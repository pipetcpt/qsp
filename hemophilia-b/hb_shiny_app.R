##############################################################################
# Hemophilia B QSP — Shiny Interactive Dashboard
# Tabs: (1) Patient Profile · (2) FIX Replacement PK · (3) Non-Factor Rebalancing
#       (4) AAV Gene Therapy · (5) PD Core Metrics · (6) Bleed Risk & Arthropathy
#       (7) Scenario Comparison · (8) Inhibitor & Biomarkers
##############################################################################

library(shiny)
library(bslib)
library(plotly)
library(dplyr)
library(tidyr)
library(ggplot2)

# ── Inline mini-ODE simulator (self-contained Euler integrator; mirrors
#    hb_mrgsolve_model.R without requiring an mrgsolve/C++ toolchain) ──────

simulate_hb <- function(
    BW           = 70,
    base_fix     = 0.5,     # IU/dL endogenous FIX (severe <1, moderate 1-5, mild >5)
    tx_type      = "shl",   # none | shl | ehl_fc | ehl_alb | ehl_gp | conc | mars | fitu | aav
    dose_fix     = 40,      # IU/kg per replacement dose
    freq_fix     = 2,       # times per week (SHL) — ignored for interval EHL regimens
    interval_ehl = 7,       # days between EHL doses
    dose_conc_maint = 15,   # mg concizumab maintenance (daily)
    dose_mars_maint = 150,  # mg marstacimab maintenance (weekly)
    dose_fitu    = 50,      # mg/month fitusiran
    aav_dose     = 100,     # relative AAV vector genome dose units
    use_steroid  = TRUE,    # corticosteroid taper for AAV transaminitis
    nab_block    = 0,       # 0-1, pre-existing anti-AAV neutralizing antibody block
    inhibitor_bu = 0,       # pre-existing inhibitor titer (BU/mL)
    null_mut     = FALSE,   # null-mutation genotype -> higher inhibitor formation risk
    duration_d   = 365
) {
  dt      <- 1.0  # hours
  n_steps <- duration_d * 24

  # --- FIX SHL PK ---
  CL_FIX <- 0.30; Vc_FIX <- 6.0; Q_FIX <- 0.15; Vp_FIX <- 4.0

  # --- FIX EHL PK variants ---
  ehl_par <- switch(tx_type,
    ehl_fc  = list(CL = 0.075, Vc = 6.0, Q = 0.05, Vp = 4.0),
    ehl_alb = list(CL = 0.048, Vc = 5.8, Q = 0.03, Vp = 3.8),
    ehl_gp  = list(CL = 0.065, Vc = 6.0, Q = 0.04, Vp = 4.0),
    list(CL = 0.075, Vc = 6.0, Q = 0.05, Vp = 4.0)
  )

  # --- AAV gene therapy ---
  ka_AAV <- 0.35; k_transduce <- 1.00; k_expr_ramp <- 0.0015
  Expr_plateau_max <- 40.0
  k_capsid_immune <- 0.006; k_immune_decay <- 0.004
  k_ALT_rise <- 0.08; k_ALT_fall <- 0.03
  steroid_suppress <- 0.75; k_vector_dilution <- 0.0000099
  k_antigen_clear <- 0.000825; k_immune_erosion <- 0.0015

  # --- Concizumab / Marstacimab PK ---
  ka_CONC <- 0.020; CL_CONC <- 0.0021; Vc_CONC <- 3.4; F_CONC <- 0.65; EC50_CONC <- 45.0
  ka_MARS <- 0.018; CL_MARS <- 0.0032; Vc_MARS <- 4.6; F_MARS <- 0.70; EC50_MARS <- 60.0
  Emax_TFPI <- 0.90; FIXeq_TFPI <- 18.0; FIXeq_AT <- 12.0

  # --- Fitusiran + AT ---
  ka_FITU <- 0.008; CL_FITU <- 0.022; Vc_FITU <- 8.5; F_FITU <- 0.75
  Emax_FITU <- 0.92; EC50_FITU <- 0.008
  ksyn_ATm <- 0.0065; kdeg_ATm <- 0.0065
  ksyn_ATp <- 0.0030; kdeg_ATp <- 0.0030

  # --- Inhibitor ---
  IC50_inh <- 1.0; Ki_max <- 150.0
  k_inhibit <- 0.00010 * (if (null_mut) 3.0 else 1.0)
  k_inh_off <- 0.0006

  # --- Thrombin generation ---
  ETP_base <- 100.0; ETP_EC50 <- 4.0; ETP_hill <- 0.85; k_ETP_up <- 0.85

  # --- Bleed model ---
  ABR_base <- 28.0; FIX_ABR_EC50 <- 3.5; FIX_ABR_hill <- 1.15; ABR_floor <- 0.4

  # --- Arthropathy ---
  k_joint_in <- 0.0007; k_joint_rep <- 0.0001
  k_syno_in  <- 0.002;  k_syno_out  <- 0.0005

  # --- QoL ---
  k_QoL_joint <- 0.004; k_QoL_ABR <- 0.010

  # --- State vector ---
  FIX_C <- base_fix; FIX_P <- base_fix
  FIXe_C <- 0.0; FIXe_P <- 0.0
  AAV_Vector <- 0.0; Transduced_Hep <- 0.0; Transgene_Expr <- 0.0; Capsid_Antigen <- 0.0
  Capsid_Immune <- 0.0; ALT_level <- 1.0
  CONC_SC <- 0.0; CONC_C <- 0.0
  MARS_SC <- 0.0; MARS_C <- 0.0
  FITU_SC <- 0.0; FITU_C <- 0.0
  AT_mRNA <- 1.0; AT_prot <- 1.0
  Inhibitor <- inhibitor_bu
  ETP <- 18.0; CumBleeds <- 0.0
  JointScore <- 6.0; Synovitis <- 0.1; QoL <- 0.78

  use_shl  <- tx_type == "shl"
  use_ehl  <- tx_type %in% c("ehl_fc", "ehl_alb", "ehl_gp")
  use_conc <- tx_type == "conc"
  use_mars <- tx_type == "mars"
  use_fitu <- tx_type == "fitu"
  use_aav  <- tx_type == "aav"

  fix_dose_iudl  <- dose_fix * 1.0        # ~1 IU/dL per IU/kg recovery for FIX
  shl_interval_h <- 168 / freq_fix
  ehl_times_h    <- seq(0, duration_d * 24, by = interval_ehl * 24)

  conc_loading_mg <- 210; conc_times_h <- seq(24, duration_d * 24, by = 24)
  mars_loading_mg <- 300; mars_times_h <- seq(24, duration_d * 24, by = 7 * 24)
  fitu_times_h    <- seq(0, duration_d * 24, by = 28 * 24)

  rec_every <- 6
  out_n <- floor(n_steps / rec_every) + 1
  rec <- matrix(NA_real_, nrow = out_n, ncol = 15)
  colnames(rec) <- c("time_days", "FIX_activity", "FIX_total", "Transgene",
                      "Concizumab", "Marstacimab", "Fitusiran", "AT_prot",
                      "ETP", "ABR", "JointScore", "QoL", "Inhibitor",
                      "ALT_level", "ImmuneActivation")
  rec_idx <- 1

  for (step in 1:n_steps) {
    t_h <- (step - 1) * dt

    if (use_shl && isTRUE(abs(t_h %% shl_interval_h) < dt / 2)) {
      FIX_C <- FIX_C + fix_dose_iudl
    }
    if (use_ehl && any(abs(ehl_times_h - t_h) < dt / 2)) {
      FIXe_C <- FIXe_C + fix_dose_iudl
    }
    if (use_aav && t_h < dt) {
      AAV_Vector <- AAV_Vector + aav_dose
    }
    if (use_conc) {
      if (t_h < dt) CONC_SC <- CONC_SC + conc_loading_mg
      if (any(abs(conc_times_h - t_h) < dt / 2)) CONC_SC <- CONC_SC + dose_conc_maint
    }
    if (use_mars) {
      if (t_h < dt) MARS_SC <- MARS_SC + mars_loading_mg
      if (any(abs(mars_times_h - t_h) < dt / 2)) MARS_SC <- MARS_SC + dose_mars_maint
    }
    if (use_fitu && any(abs(fitu_times_h - t_h) < dt / 2)) {
      FITU_SC <- FITU_SC + dose_fitu
    }

    FIX_replace_raw <- FIX_C + FIXe_C
    FIX_act <- max(FIX_replace_raw / (1.0 + Inhibitor / IC50_inh), 0)

    TFPI_effect <- 0.0
    if (use_conc) TFPI_effect <- Emax_TFPI * CONC_C / (EC50_CONC + CONC_C + 1e-9) * FIXeq_TFPI
    if (use_mars) TFPI_effect <- Emax_TFPI * MARS_C / (EC50_MARS + MARS_C + 1e-9) * FIXeq_TFPI

    # FIX SHL PK
    dFIX_C <- (-CL_FIX/Vc_FIX*FIX_C - Q_FIX/Vc_FIX*FIX_C + Q_FIX/Vp_FIX*FIX_P) * dt
    dFIX_P <- (Q_FIX/Vc_FIX*FIX_C - Q_FIX/Vp_FIX*FIX_P) * dt
    FIX_C <- max(FIX_C + dFIX_C, 0); FIX_P <- max(FIX_P + dFIX_P, 0)

    # FIX EHL PK
    dFIXe_C <- (-ehl_par$CL/ehl_par$Vc*FIXe_C - ehl_par$Q/ehl_par$Vc*FIXe_C + ehl_par$Q/ehl_par$Vp*FIXe_P) * dt
    dFIXe_P <- (ehl_par$Q/ehl_par$Vc*FIXe_C - ehl_par$Q/ehl_par$Vp*FIXe_P) * dt
    FIXe_C <- max(FIXe_C + dFIXe_C, 0); FIXe_P <- max(FIXe_P + dFIXe_P, 0)

    # AAV gene therapy: transduction flux drives both hepatocyte transduction and a
    # transient intracellular capsid-antigen pool (~5-week clearance) that resolves
    # even though the durable episomal FIX-Padua transgene persists — so immune-
    # mediated erosion is self-limiting rather than an indefinite drag on expression.
    transduce_rate <- if (use_aav) k_transduce * (1 - nab_block) else 0
    dAAV <- -ka_AAV * AAV_Vector * dt
    transduction_flux <- transduce_rate * AAV_Vector/100 * (1 - Transduced_Hep)

    dAntigen <- (transduction_flux - k_antigen_clear * Capsid_Antigen) * dt
    Capsid_Antigen <- max(Capsid_Antigen + dAntigen, 0)

    immune_drive <- k_capsid_immune * Capsid_Antigen
    if (use_steroid && use_aav) immune_drive <- immune_drive * (1 - steroid_suppress)
    dImmune <- (immune_drive * (1 - Capsid_Immune) - k_immune_decay * Capsid_Immune) * dt
    Capsid_Immune <- max(min(Capsid_Immune + dImmune, 1), 0)
    dALT <- (k_ALT_rise * Capsid_Immune * (1 - ALT_level/5) - k_ALT_fall * (ALT_level - 1)) * dt
    ALT_level <- max(ALT_level + dALT, 1)

    dTransd <- (transduction_flux - k_vector_dilution * Transduced_Hep -
                  k_immune_erosion * Capsid_Immune * Transduced_Hep) * dt
    AAV_Vector <- max(AAV_Vector + dAAV, 0)
    Transduced_Hep <- max(min(Transduced_Hep + dTransd, 1), 0)

    Expr_target <- Expr_plateau_max * Transduced_Hep
    dExpr <- k_expr_ramp * (Expr_target - Transgene_Expr) * dt
    Transgene_Expr <- max(Transgene_Expr + dExpr, 0)

    # Concizumab / Marstacimab PK
    dCONC_SC <- -ka_CONC * CONC_SC * dt
    dCONC_C  <- (F_CONC * ka_CONC * CONC_SC * 1000 / Vc_CONC - CL_CONC/Vc_CONC*CONC_C) * dt
    CONC_SC <- max(CONC_SC + dCONC_SC, 0); CONC_C <- max(CONC_C + dCONC_C, 0)

    dMARS_SC <- -ka_MARS * MARS_SC * dt
    dMARS_C  <- (F_MARS * ka_MARS * MARS_SC * 1000 / Vc_MARS - CL_MARS/Vc_MARS*MARS_C) * dt
    MARS_SC <- max(MARS_SC + dMARS_SC, 0); MARS_C <- max(MARS_C + dMARS_C, 0)

    # Fitusiran PK + AT knockdown
    dFITU_SC <- -ka_FITU * FITU_SC * dt
    dFITU_C  <- (F_FITU * ka_FITU * FITU_SC / Vc_FITU - CL_FITU/Vc_FITU*FITU_C) * dt
    FITU_SC <- max(FITU_SC + dFITU_SC, 0); FITU_C <- max(FITU_C + dFITU_C, 0)

    kd_ATm_total <- kdeg_ATm
    if (use_fitu) kd_ATm_total <- kdeg_ATm * (1 + Emax_FITU * FITU_C / (EC50_FITU + FITU_C + 1e-9))
    dAT_mRNA <- (ksyn_ATm - kd_ATm_total*AT_mRNA) * dt
    dAT_prot <- (ksyn_ATp*AT_mRNA - kdeg_ATp*AT_prot) * dt
    AT_mRNA <- max(AT_mRNA + dAT_mRNA, 0); AT_prot <- max(AT_prot + dAT_prot, 0)

    AT_rebalance_effect <- FIXeq_AT * (1 - AT_prot)
    FIX_total <- FIX_act + Transgene_Expr + TFPI_effect + AT_rebalance_effect

    # Inhibitor
    inhibit_form <- if (FIX_act > 0) k_inhibit * FIX_act * (1 - Inhibitor/Ki_max) else 0
    dInhib <- (inhibit_form - k_inh_off*Inhibitor) * dt
    Inhibitor <- max(Inhibitor + dInhib, 0)

    # ETP
    ETP_ss <- ETP_base * FIX_total^ETP_hill / (ETP_EC50^ETP_hill + FIX_total^ETP_hill + 1e-9) *
      (1 + (1 - AT_prot) * 0.5)
    dETP <- k_ETP_up * (ETP_ss - ETP) * dt
    ETP <- max(ETP + dETP, 0)

    # Bleed rate
    FIX_prot <- FIX_total^FIX_ABR_hill / (FIX_ABR_EC50^FIX_ABR_hill + FIX_total^FIX_ABR_hill + 1e-9)
    ABR_inst <- max(ABR_base * (1 - FIX_prot) + ABR_floor, ABR_floor)
    bleed_per_h <- ABR_inst / 8760

    dSyno <- (k_syno_in*bleed_per_h*(1 - Synovitis) - k_syno_out*Synovitis) * dt
    Synovitis <- max(min(Synovitis + dSyno, 1), 0)

    dJoint <- (k_joint_in*bleed_per_h*(100 - JointScore) - k_joint_rep*JointScore*(1 - JointScore/100)) * dt
    JointScore <- max(min(JointScore + dJoint, 100), 0)

    QoL_tgt <- max(1.0 - k_QoL_joint*JointScore/100 - k_QoL_ABR*ABR_inst/28, 0.1)
    dQoL <- 0.01 * (QoL_tgt - QoL) * dt
    QoL <- max(min(QoL + dQoL, 1.0), 0.1)

    if ((step - 1) %% rec_every == 0) {
      rec[rec_idx, ] <- c(t_h/24, FIX_act, FIX_total, Transgene_Expr,
                           CONC_C, MARS_C, FITU_C, AT_prot,
                           ETP, ABR_inst, JointScore, QoL, Inhibitor,
                           ALT_level, Capsid_Immune)
      rec_idx <- rec_idx + 1
    }
  }

  as.data.frame(rec[1:(rec_idx - 1), ])
}

# ── UI ───────────────────────────────────────────────────────────────────────

tx_choices <- c(
  "No prophylaxis (on-demand)"     = "none",
  "SHL-rFIX 2x/week"                = "shl",
  "EHL-rFIX-Fc (Alprolix-like)"     = "ehl_fc",
  "EHL-rFIX-Albumin (Idelvion-like)" = "ehl_alb",
  "GlycoPEG-rFIX (Refixia-like)"    = "ehl_gp",
  "Concizumab SC daily"             = "conc",
  "Marstacimab SC weekly"           = "mars",
  "Fitusiran SC monthly"            = "fitu",
  "AAV Gene Therapy (single dose)"  = "aav"
)

ui <- page_navbar(
  title = "Hemophilia B QSP Dashboard",
  theme = bs_theme(
    bootswatch = "darkly",
    primary    = "#2b9348",
    secondary  = "#52b788",
    base_font  = font_google("Roboto"),
    heading_font = font_google("Roboto Condensed")
  ),
  bg = "#1a1a2e",
  fillable = TRUE,

  # ── Tab 1: Patient Profile ────────────────────────────────────────────────
  nav_panel(
    title = "1. Patient Profile",
    icon  = icon("user-circle"),
    layout_sidebar(
      sidebar = sidebar(
        width = 320, bg = "#16213e",
        h5("Patient Parameters", class = "text-info"),
        sliderInput("BW", "Body Weight (kg)", 10, 120, 70, step = 5),
        sliderInput("base_fix", "Endogenous FIX (%)", 0, 5, 0.5, step = 0.5),
        selectInput("severity", "HB Severity",
                    choices = c("Severe (<1%)" = "severe",
                                "Moderate (1-5%)" = "moderate",
                                "Mild (>5%)" = "mild")),
        checkboxInput("null_mut", "Null-mutation genotype (higher inhibitor risk)", FALSE),
        sliderInput("inhibitor_bu", "Pre-existing Inhibitor (BU/mL)", 0, 40, 0, step = 1),
        hr(),
        h5("Treatment", class = "text-info"),
        selectInput("tx_type", "Treatment Modality", choices = tx_choices, selected = "shl"),
        conditionalPanel(
          "input.tx_type == 'shl'",
          sliderInput("dose_fix", "FIX Dose (IU/kg)", 10, 100, 40, step = 5),
          sliderInput("freq_fix", "Frequency (x/week)", 1, 7, 2, step = 1)
        ),
        conditionalPanel(
          "input.tx_type == 'ehl_fc' || input.tx_type == 'ehl_alb' || input.tx_type == 'ehl_gp'",
          sliderInput("dose_fix_ehl", "EHL-FIX Dose (IU/kg)", 25, 100, 50, step = 5),
          sliderInput("interval_ehl", "Dosing Interval (days)", 7, 14, 7, step = 7)
        ),
        conditionalPanel(
          "input.tx_type == 'aav'",
          sliderInput("nab_block", "Pre-existing Anti-AAV NAb Block", 0, 1, 0, step = 0.1),
          checkboxInput("use_steroid", "Corticosteroid taper for transaminitis", TRUE)
        ),
        sliderInput("duration_d", "Simulation Duration (days)", 90, 3650, 365, step = 30),
        actionButton("run_sim", "Run Simulation", class = "btn-primary btn-lg w-100", icon = icon("play"))
      ),
      card(
        card_header("Patient Summary"),
        layout_columns(
          col_widths = c(4, 4, 4),
          value_box("Severity", textOutput("sev_box"), theme = "primary", showcase = icon("tint")),
          value_box("Inhibitor", textOutput("inhib_box"), theme = "warning", showcase = icon("shield-virus")),
          value_box("Treatment", textOutput("tx_box"), theme = "success", showcase = icon("pills"))
        ),
        p("Hemophilia B is an X-linked recessive bleeding disorder caused by
          deficiency of coagulation factor IX (FIX). Severity is classified by
          endogenous FIX activity: severe (<1 IU/dL), moderate (1-5 IU/dL), or
          mild (>5 IU/dL). Roughly 1-3% of severe HB patients develop inhibitors
          (anti-FIX antibodies), which uniquely carry risks of anaphylaxis and
          nephrotic syndrome during immune tolerance induction.",
          class = "text-muted small mt-2"),
        tableOutput("patient_summary_tbl")
      )
    )
  ),

  # ── Tab 2: FIX Replacement PK ──────────────────────────────────────────────
  nav_panel(
    title = "2. FIX PK",
    icon  = icon("chart-line"),
    layout_sidebar(
      sidebar = sidebar(
        width = 240, bg = "#16213e",
        h5("PK Display Options", class = "text-info"),
        checkboxInput("log_scale", "Log Y-axis", FALSE),
        sliderInput("pk_days", "Display (days)", 1, 90, 28),
        hr(),
        helpText("FIX half-life (recovery ~1 IU/dL per IU/kg):"),
        tags$ul(class = "text-muted small",
          tags$li("SHL-rFIX (BeneFIX/Rixubis): ~18-24 h"),
          tags$li("EHL-rFIX-Fc (Alprolix): ~82 h"),
          tags$li("EHL-rFIX-Albumin (Idelvion): ~102-104 h"),
          tags$li("GlycoPEG-rFIX (Refixia): ~93 h")
        )
      ),
      card(
        card_header("FIX Activity vs. Time"),
        plotlyOutput("pk_plot", height = "400px")
      ),
      card(
        card_header("PK Summary Statistics"),
        tableOutput("pk_summary_tbl")
      )
    )
  ),

  # ── Tab 3: Non-Factor Rebalancing ──────────────────────────────────────────
  nav_panel(
    title = "3. Non-Factor Rebalancing",
    icon  = icon("balance-scale-right"),
    layout_columns(
      col_widths = c(6, 6),
      card(
        card_header("Concizumab / Marstacimab (Anti-TFPI) Concentration"),
        plotlyOutput("tfpi_plot", height = "320px"),
        p("Anti-TFPI monoclonal antibodies neutralize tissue factor pathway
          inhibitor, amplifying the extrinsic pathway to rebalance hemostasis
          independent of FIX replacement.", class = "text-muted small")
      ),
      card(
        card_header("Fitusiran — Antithrombin Knockdown"),
        plotlyOutput("at_plot", height = "320px"),
        p("Fitusiran (siRNA) knocks down antithrombin mRNA (t1/2 ~4.4 days) and
          protein (t1/2 ~9.6 days), reducing AT to ~15-30% of baseline and
          augmenting thrombin generation.", class = "text-muted small")
      )
    )
  ),

  # ── Tab 4: AAV Gene Therapy ────────────────────────────────────────────────
  nav_panel(
    title = "4. Gene Therapy",
    icon  = icon("dna"),
    layout_columns(
      col_widths = c(7, 5),
      card(
        card_header("Endogenous FIX-Padua Transgene Expression"),
        plotlyOutput("aav_expr_plot", height = "350px")
      ),
      card(
        card_header("Capsid Immune Response & Transaminitis"),
        plotlyOutput("aav_alt_plot", height = "350px")
      )
    ),
    card(
      p("Single-dose IV administration of an AAV5 (etranacogene dezaparvovec-like)
        or AAVRh74var (fidanacogene elaparvovec-like) vector transduces
        hepatocytes and drives episomal expression of the hyperfunctional
        FIX-Padua (R338L) variant (~8x specific activity). Capsid-specific
        CD8+ T-cell responses can cause transient transaminitis, managed with
        a reactive corticosteroid taper; durability may decline over years due
        to hepatocyte turnover diluting non-integrating episomal vector.",
        class = "text-muted small")
    )
  ),

  # ── Tab 5: PD Core Metrics ─────────────────────────────────────────────────
  nav_panel(
    title = "5. PD Core Metrics",
    icon  = icon("heartbeat"),
    layout_columns(
      col_widths = c(6, 6),
      card(
        card_header("Thrombin Generation Potential (ETP)"),
        plotlyOutput("etp_plot", height = "320px"),
        p("ETP reflects thrombin generation capacity. Normal: ~100 (normalized).
          Severe HB: ~10-20% of normal. Prophylaxis restores ETP proportionally
          to FIX trough levels.", class = "text-muted small")
      ),
      card(
        card_header("Total Effective Hemostatic FIX-Equivalent"),
        plotlyOutput("fixtotal_plot", height = "320px"),
        p("Sum of replacement-derived FIX activity, AAV transgene expression,
          and TFPI-mAb FIX-equivalent contribution.", class = "text-muted small")
      )
    )
  ),

  # ── Tab 6: Bleed Risk & Arthropathy ────────────────────────────────────────
  nav_panel(
    title = "6. Bleed Risk & Arthropathy",
    icon  = icon("droplet"),
    layout_columns(
      col_widths = c(7, 5),
      card(
        card_header("Annual Bleed Rate (ABR) Over Time"),
        plotlyOutput("abr_plot", height = "350px")
      ),
      card(
        card_header("Hemophilic Arthropathy & QoL"),
        plotlyOutput("joint_qol_plot", height = "350px")
      )
    ),
    layout_columns(
      col_widths = c(4, 4, 4),
      value_box("Mean ABR", textOutput("mean_abr_vb"), theme = "danger", showcase = icon("tint")),
      value_box("Final Joint Score", textOutput("joint_vb"), theme = "warning", showcase = icon("bone")),
      value_box("Final QoL", textOutput("qol_vb"), theme = "success", showcase = icon("smile"))
    )
  ),

  # ── Tab 7: Scenario Comparison ─────────────────────────────────────────────
  nav_panel(
    title = "7. Scenario Comparison",
    icon  = icon("chart-bar"),
    layout_sidebar(
      sidebar = sidebar(
        width = 300, bg = "#16213e",
        h5("Compare Scenarios", class = "text-info"),
        checkboxGroupInput("compare_scenarios", "Include Scenarios",
          choices = tx_choices,
          selected = c("none", "shl", "ehl_fc", "conc", "aav")
        ),
        sliderInput("cmp_duration", "Duration (years)", 1, 10, 2),
        actionButton("run_comparison", "Compare", class = "btn-info w-100", icon = icon("balance-scale"))
      ),
      card(
        card_header("ABR Comparison Across Scenarios"),
        plotlyOutput("cmp_abr", height = "320px")
      ),
      card(
        card_header("Joint Score Comparison (Long-term)"),
        plotlyOutput("cmp_joint", height = "320px")
      )
    )
  ),

  # ── Tab 8: Inhibitor & Biomarkers ──────────────────────────────────────────
  nav_panel(
    title = "8. Inhibitor & Biomarkers",
    icon  = icon("microscope"),
    layout_columns(
      col_widths = c(6, 6),
      card(
        card_header("FIX Activity vs. ETP Relationship"),
        plotlyOutput("fix_etp_scatter", height = "300px"),
        p("Non-linear relationship between FIX activity and thrombin
          generation (Hill equation). EC50 for ETP ~4 IU/dL.",
          class = "text-muted small")
      ),
      card(
        card_header("Inhibitor Titer Dynamics"),
        plotlyOutput("inhibitor_plot", height = "300px"),
        p("Anti-FIX antibody (inhibitor) titer. Null-mutation genotypes carry
          higher inhibitor risk. High titers necessitate bypassing agents and
          preclude standard immune tolerance induction due to nephrotic
          syndrome/anaphylaxis risk.", class = "text-muted small")
      )
    ),
    card(
      card_header("Clinical Outcomes Summary Table"),
      tableOutput("biomarker_summary")
    )
  )
)

# ── Server ───────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Simulating...", value = 0.5, {
      dose_val <- if (input$tx_type %in% c("ehl_fc", "ehl_alb", "ehl_gp")) input$dose_fix_ehl else input$dose_fix
      simulate_hb(
        BW = input$BW, base_fix = input$base_fix, tx_type = input$tx_type,
        dose_fix = dose_val, freq_fix = input$freq_fix,
        interval_ehl = if (!is.null(input$interval_ehl)) input$interval_ehl else 7,
        nab_block = if (!is.null(input$nab_block)) input$nab_block else 0,
        use_steroid = if (!is.null(input$use_steroid)) input$use_steroid else TRUE,
        inhibitor_bu = input$inhibitor_bu, null_mut = input$null_mut,
        duration_d = input$duration_d
      )
    })
  }, ignoreNULL = FALSE)

  observeEvent(TRUE, {
    shinyjs::runjs("document.getElementById('run_sim').click();")
  }, once = TRUE, ignoreInit = FALSE)

  # ── Tab 1 ──
  output$sev_box <- renderText({
    switch(input$severity, severe = "Severe (<1%)", moderate = "Moderate (1-5%)", mild = "Mild (>5%)")
  })
  output$inhib_box <- renderText({
    if (input$inhibitor_bu == 0) "Negative" else paste0(input$inhibitor_bu, " BU/mL")
  })
  output$tx_box <- renderText({ names(tx_choices)[tx_choices == input$tx_type] })
  output$patient_summary_tbl <- renderTable({
    data.frame(
      Parameter = c("Body Weight", "FIX Severity", "Genotype", "Inhibitor Status", "Simulation Duration"),
      Value = c(paste0(input$BW, " kg"), input$severity,
                if (input$null_mut) "Null mutation" else "Missense/other",
                if (input$inhibitor_bu == 0) "None" else paste0(input$inhibitor_bu, " BU/mL"),
                paste0(input$duration_d, " days"))
    )
  })

  # ── Tab 2 PK ──
  output$pk_plot <- renderPlotly({
    df <- sim_data() %>% filter(time_days <= input$pk_days)
    plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~FIX_total, name = "Total Effective FIX (IU/dL)", line = list(color = "#2b9348", width = 2)) %>%
      add_lines(y = ~FIX_activity, name = "Replacement FIX (IU/dL)", line = list(color = "#52b788", dash = "dash")) %>%
      layout(
        xaxis = list(title = "Time (days)", color = "white"),
        yaxis = list(title = "FIX Activity (IU/dL)", type = if (input$log_scale) "log" else "linear", color = "white"),
        shapes = list(
          list(type = "line", x0 = 0, x1 = input$pk_days, y0 = 1, y1 = 1, line = list(color = "red", dash = "dot")),
          list(type = "line", x0 = 0, x1 = input$pk_days, y0 = 15, y1 = 15, line = list(color = "orange", dash = "dot"))
        ),
        annotations = list(
          list(x = 1, y = 2, text = "1% trough", showarrow = FALSE, font = list(color = "red", size = 10)),
          list(x = 1, y = 17, text = "15% optimal", showarrow = FALSE, font = list(color = "orange", size = 10))
        ),
        paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e",
        font = list(color = "white"), legend = list(font = list(color = "white"))
      )
  })
  output$pk_summary_tbl <- renderTable({
    df <- sim_data() %>% filter(time_days > 0)
    data.frame(
      Metric = c("Mean FIX Activity (IU/dL)", "Min FIX (trough)", "Max FIX (peak)", "Mean Total FIX-Equivalent"),
      Value = round(c(mean(df$FIX_activity), min(df$FIX_activity), max(df$FIX_activity), mean(df$FIX_total)), 2)
    )
  })

  # ── Tab 3 Non-factor ──
  output$tfpi_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~Concizumab, name = "Concizumab (ng/mL)", line = list(color = "#f8961e", width = 2)) %>%
      add_lines(y = ~Marstacimab, name = "Marstacimab (ng/mL)", line = list(color = "#f94144", dash = "dash")) %>%
      layout(xaxis = list(title = "Time (days)", color = "white"),
             yaxis = list(title = "Concentration (ng/mL)", color = "white"),
             paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e",
             font = list(color = "white"), legend = list(font = list(color = "white")))
  })
  output$at_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~AT_prot * 100, name = "AT Protein (%)", line = list(color = "#f8961e", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)", color = "white"),
             yaxis = list(title = "Antithrombin Level (% of baseline)", range = c(0, 110), color = "white"),
             shapes = list(list(type = "line", x0 = 0, x1 = max(df$time_days), y0 = 25, y1 = 25,
                                 line = list(color = "red", dash = "dot"))),
             annotations = list(list(x = 5, y = 30, text = "Fitusiran target ~15-25%",
                                      showarrow = FALSE, font = list(color = "red", size = 10))),
             paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e", font = list(color = "white"))
  })

  # ── Tab 4 AAV ──
  output$aav_expr_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~Transgene, name = "Transgene FIX Expression (IU/dL)", line = list(color = "#2b9348", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)", color = "white"),
             yaxis = list(title = "Endogenous FIX Expression (IU/dL)", color = "white"),
             paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e", font = list(color = "white"))
  })
  output$aav_alt_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~ALT_level, name = "ALT (fold-elevation)", line = list(color = "#d90429", width = 2)) %>%
      add_lines(y = ~ImmuneActivation * 5, name = "Capsid Immune Activation (x5)",
                line = list(color = "#f8961e", dash = "dash")) %>%
      layout(xaxis = list(title = "Time (days)", color = "white"),
             yaxis = list(title = "ALT fold / Immune Activation", color = "white"),
             paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e",
             font = list(color = "white"), legend = list(font = list(color = "white")))
  })

  # ── Tab 5 PD ──
  output$etp_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~ETP, name = "Thrombin ETP", line = list(color = "#06d6a0", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)", color = "white"),
             yaxis = list(title = "ETP (normalized)", range = c(0, 110), color = "white"),
             paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e", font = list(color = "white"))
  })
  output$fixtotal_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~FIX_total, name = "Total FIX-Equivalent", line = list(color = "#43aa8b", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)", color = "white"),
             yaxis = list(title = "FIX-Equivalent (IU/dL)", color = "white"),
             paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e", font = list(color = "white"))
  })

  # ── Tab 6 Bleed ──
  output$abr_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~ABR, name = "Instantaneous ABR", line = list(color = "#e63946", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)", color = "white"),
             yaxis = list(title = "Annual Bleed Rate", color = "white"),
             shapes = list(list(type = "line", x0 = 0, x1 = max(df$time_days), y0 = 3, y1 = 3,
                                 line = list(color = "green", dash = "dot"))),
             paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e", font = list(color = "white"))
  })
  output$joint_qol_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~JointScore, name = "Joint Damage Score", line = list(color = "#f3722c", width = 2)) %>%
      add_lines(y = ~QoL * 100, name = "QoL (x100)", line = list(color = "#43aa8b", dash = "dash")) %>%
      layout(xaxis = list(title = "Time (days)", color = "white"),
             yaxis = list(title = "Score", color = "white"),
             paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e",
             font = list(color = "white"), legend = list(font = list(color = "white")))
  })
  output$mean_abr_vb <- renderText({ round(mean(sim_data()$ABR), 1) })
  output$joint_vb <- renderText({ round(tail(sim_data()$JointScore, 1), 1) })
  output$qol_vb <- renderText({ round(tail(sim_data()$QoL, 1), 2) })

  # ── Tab 7 Comparison ──
  cmp_data <- eventReactive(input$run_comparison, {
    withProgress(message = "Comparing scenarios...", value = 0.3, {
      bind_rows(lapply(input$compare_scenarios, function(s) {
        simulate_hb(BW = input$BW, base_fix = input$base_fix, tx_type = s,
                    duration_d = input$cmp_duration * 365) %>%
          mutate(scenario = names(tx_choices)[tx_choices == s])
      }))
    })
  }, ignoreNULL = FALSE)

  output$cmp_abr <- renderPlotly({
    df <- cmp_data()
    plot_ly(df, x = ~time_days, y = ~ABR, color = ~scenario) %>%
      add_lines() %>%
      layout(xaxis = list(title = "Time (days)", color = "white"),
             yaxis = list(title = "ABR (bleeds/year)", color = "white"),
             paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e",
             font = list(color = "white"), legend = list(font = list(color = "white")))
  })
  output$cmp_joint <- renderPlotly({
    df <- cmp_data()
    plot_ly(df, x = ~time_days, y = ~JointScore, color = ~scenario) %>%
      add_lines() %>%
      layout(xaxis = list(title = "Time (days)", color = "white"),
             yaxis = list(title = "Joint Damage Score", color = "white"),
             paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e",
             font = list(color = "white"), legend = list(font = list(color = "white")))
  })

  # ── Tab 8 Biomarkers ──
  output$fix_etp_scatter <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~FIX_total, y = ~ETP, type = "scatter", mode = "markers",
            marker = list(color = "#06d6a0", size = 5, opacity = 0.6)) %>%
      layout(xaxis = list(title = "FIX-Equivalent Activity (IU/dL)", color = "white"),
             yaxis = list(title = "Thrombin ETP", color = "white"),
             paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e", font = list(color = "white"))
  })
  output$inhibitor_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~Inhibitor, name = "Inhibitor Titer (BU/mL)", line = list(color = "#d90429", width = 2)) %>%
      layout(xaxis = list(title = "Time (days)", color = "white"),
             yaxis = list(title = "Inhibitor Titer (BU/mL)", color = "white"),
             paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e", font = list(color = "white"))
  })
  output$biomarker_summary <- renderTable({
    df <- sim_data()
    data.frame(
      Metric = c("Mean ABR (bleeds/yr)", "Final Joint Score", "Final QoL",
                 "Mean Inhibitor Titer (BU/mL)", "Mean ETP", "Mean FIX-Equivalent (IU/dL)"),
      Value = round(c(mean(df$ABR), tail(df$JointScore, 1), tail(df$QoL, 1),
                      mean(df$Inhibitor), mean(df$ETP), mean(df$FIX_total)), 2)
    )
  })
}

shinyApp(ui, server)
