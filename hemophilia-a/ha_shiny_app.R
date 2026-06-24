##############################################################################
# Hemophilia A QSP — Shiny Interactive Dashboard
# Tabs: (1) Patient Profile · (2) FVIII PK · (3) PD Core Metrics
#       (4) Bleed Risk & ABR · (5) Scenario Comparison · (6) Biomarkers
##############################################################################

library(shiny)
library(bslib)
library(plotly)
library(dplyr)
library(tidyr)
library(ggplot2)

# ── Inline mini-ODE simulator (no mrgsolve dependency for Shiny) ──────────

simulate_ha <- function(
    BW          = 70,
    base_fviii  = 0.5,   # IU/dL — endogenous FVIII (1 = 1%; 0 = severe)
    use_shl     = TRUE,
    use_emic    = FALSE,
    use_fitu    = FALSE,
    dose_fviii  = 25,    # IU/kg (SHL or EHL dose)
    freq_fviii  = 3,     # times per week (SHL)
    dose_emic   = 1.5,   # mg/kg maintenance (after loading)
    dose_fitu   = 80,    # mg/month
    inhibitor_bu = 0,    # pre-existing inhibitor titer (BU/mL)
    duration_d  = 365
) {
  # --- Parameters ---
  dt      <- 1.0  # hours
  n_steps <- duration_d * 24

  # FVIII PK
  CL_FVIII <- if (use_shl) 3.0  else 1.8
  Vc_FVIII <- 3.2
  Q_FVIII  <- if (use_shl) 1.8  else 1.0
  Vp_FVIII <- 2.5

  # Emicizumab PK
  ka_EMIC   <- 0.014
  CL_EMIC   <- 0.0024
  Vc_EMIC   <- 2.9
  Q_EMIC    <- 0.0015
  Vp_EMIC   <- 2.7
  F_EMIC    <- 0.80
  EC50_EMIC <- 0.045
  Emax_EMIC <- 0.85

  # Fitusiran PK
  ka_FITU   <- 0.008
  CL_FITU   <- 0.022
  Vc_FITU   <- 8.5
  F_FITU    <- 0.75
  Emax_FITU <- 0.92
  EC50_FITU <- 0.008
  ksyn_ATm  <- 0.0065; kdeg_ATm <- 0.0065
  ksyn_ATp  <- 0.0030; kdeg_ATp <- 0.0030

  # Inhibitor
  IC50_inh <- 1.0

  # Bleed model
  ABR_base       <- 30.0
  FVIII_ABR_EC50 <- 3.0
  FVIII_ABR_hill <- 1.2
  ABR_floor      <- 0.5

  # Arthropathy
  k_joint_in  <- 0.0008
  k_joint_rep <- 0.0001
  k_syno_in   <- 0.002
  k_syno_out  <- 0.0005

  # QoL
  k_QoL_joint <- 0.004
  k_QoL_ABR   <- 0.010

  # ETP
  ETP_base        <- 100.0
  ETP_EC50        <- 5.0
  ETP_hill        <- 0.8
  AT_inhibit_ETP  <- 0.3
  k_ETP_up        <- 0.90

  # --- State vector ---
  FVIII_C  <- base_fviii
  FVIII_P  <- base_fviii * Q_FVIII / Q_FVIII  # ~equal initial
  EMIC_SC  <- 0.0; EMIC_C <- 0.0; EMIC_P <- 0.0
  FITU_SC  <- 0.0; FITU_C <- 0.0
  AT_mRNA  <- 1.0; AT_prot <- 1.0
  Inhibitor <- inhibitor_bu
  ETP      <- 20.0
  CumBleeds <- 0.0
  JointScore <- 5.0
  QoL      <- 0.80
  Synovitis <- 0.1

  # Dosing schedules
  fviii_dose_iudl <- dose_fviii * 2.0  # 2% per IU/kg
  fviii_interval_h <- 168 / freq_fviii  # hours between doses

  # Loading doses: 3 mg/kg ×4 weeks Q1W
  emic_loading_mg  <- 3.0 * BW
  emic_maint_mg    <- dose_emic * BW
  emic_loading_times <- (0:3) * 7 * 24  # h
  emic_maint_times   <- (4:60) * 7 * 24  # up to 60 weeks maintenance

  fitu_dose_mg   <- dose_fitu
  fitu_times     <- seq(0, duration_d * 24, by = 28 * 24)

  # --- Output storage (record every 6 hours) ---
  rec_every <- 6
  out_n     <- floor(n_steps / rec_every) + 1
  times_rec  <- numeric(out_n)
  FVIII_rec  <- numeric(out_n)
  EMIC_rec   <- numeric(out_n)
  FITU_rec   <- numeric(out_n)
  AT_rec     <- numeric(out_n)
  ETP_rec    <- numeric(out_n)
  ABR_rec    <- numeric(out_n)
  Joint_rec  <- numeric(out_n)
  QoL_rec    <- numeric(out_n)
  Inhib_rec  <- numeric(out_n)
  rec_idx    <- 1

  for (step in 1:n_steps) {
    t_h <- (step - 1) * dt

    # ── Dosing events ──
    if (use_shl || !use_emic && !use_fitu) {
      if (isTRUE(abs(t_h %% fviii_interval_h) < dt / 2)) {
        FVIII_C <- FVIII_C + fviii_dose_iudl
      }
    }
    if (use_emic) {
      if (any(abs(emic_loading_times - t_h) < dt / 2)) {
        EMIC_SC <- EMIC_SC + emic_loading_mg
      }
      if (any(abs(emic_maint_times - t_h) < dt / 2)) {
        EMIC_SC <- EMIC_SC + emic_maint_mg
      }
    }
    if (use_fitu) {
      if (any(abs(fitu_times - t_h) < dt / 2)) {
        FITU_SC <- FITU_SC + fitu_dose_mg
      }
    }

    # ── Computed quantities ──
    FVIII_act <- FVIII_C / (1.0 + Inhibitor / IC50_inh)
    FVIII_act <- max(FVIII_act, 0)

    EMIC_eff <- if (use_emic) {
      Emax_EMIC * EMIC_C / (EC50_EMIC + EMIC_C + 1e-9) * 15.0
    } else 0.0

    FVIII_total <- FVIII_act + EMIC_eff

    # ── FVIII PK ──
    dFVIII_C <- (-CL_FVIII / Vc_FVIII * FVIII_C
                 - Q_FVIII / Vc_FVIII * FVIII_C
                 + Q_FVIII / Vp_FVIII * FVIII_P) * dt
    dFVIII_P <- (Q_FVIII / Vc_FVIII * FVIII_C
                 - Q_FVIII / Vp_FVIII * FVIII_P) * dt
    FVIII_C <- max(FVIII_C + dFVIII_C, 0)
    FVIII_P <- max(FVIII_P + dFVIII_P, 0)

    # ── Emicizumab PK ──
    dEMIC_SC <- -ka_EMIC * EMIC_SC * dt
    dEMIC_C  <- (F_EMIC * ka_EMIC * EMIC_SC / Vc_EMIC
                 - (CL_EMIC + Q_EMIC) / Vc_EMIC * EMIC_C
                 + Q_EMIC / Vp_EMIC * EMIC_P) * dt
    dEMIC_P  <- (Q_EMIC / Vc_EMIC * EMIC_C
                 - Q_EMIC / Vp_EMIC * EMIC_P) * dt
    EMIC_SC <- max(EMIC_SC + dEMIC_SC, 0)
    EMIC_C  <- max(EMIC_C  + dEMIC_C,  0)
    EMIC_P  <- max(EMIC_P  + dEMIC_P,  0)

    # ── Fitusiran PK ──
    dFITU_SC <- -ka_FITU * FITU_SC * dt
    dFITU_C  <- (F_FITU * ka_FITU * FITU_SC / Vc_FITU
                 - CL_FITU / Vc_FITU * FITU_C) * dt
    FITU_SC <- max(FITU_SC + dFITU_SC, 0)
    FITU_C  <- max(FITU_C  + dFITU_C,  0)

    # ── AT mRNA/Protein ──
    kd_ATm_total <- kdeg_ATm
    if (use_fitu) {
      kd_ATm_total <- kdeg_ATm * (1 + Emax_FITU * FITU_C / (EC50_FITU + FITU_C + 1e-9))
    }
    dAT_mRNA <- (ksyn_ATm - kd_ATm_total * AT_mRNA) * dt
    dAT_prot <- (ksyn_ATp * AT_mRNA - kdeg_ATp * AT_prot) * dt
    AT_mRNA  <- max(AT_mRNA + dAT_mRNA, 0)
    AT_prot  <- max(AT_prot + dAT_prot, 0)

    # ── ETP ──
    ETP_ss <- ETP_base * FVIII_total^ETP_hill /
              (ETP_EC50^ETP_hill + FVIII_total^ETP_hill + 1e-9) *
              (1.0 + (1.0 - AT_prot) * 0.5)
    dETP <- k_ETP_up * (ETP_ss - ETP) * dt
    ETP  <- max(ETP + dETP, 0)

    # ── Bleed rate ──
    FVIII_prot <- FVIII_total^FVIII_ABR_hill /
                  (FVIII_ABR_EC50^FVIII_ABR_hill +
                   FVIII_total^FVIII_ABR_hill + 1e-9)
    ABR_inst <- ABR_base * (1 - FVIII_prot) + ABR_floor
    ABR_inst <- max(ABR_inst, ABR_floor)

    dCumBleeds <- ABR_inst / 8760.0 * dt
    CumBleeds  <- CumBleeds + dCumBleeds

    # ── Synovitis ──
    bleed_per_h <- ABR_inst / 8760.0
    dSyno <- (k_syno_in * bleed_per_h * (1 - Synovitis)
              - k_syno_out * Synovitis) * dt
    Synovitis <- max(0, min(1, Synovitis + dSyno))

    # ── Joint score ──
    dJoint <- (k_joint_in * bleed_per_h * (100 - JointScore)
               - k_joint_rep * JointScore * (1 - JointScore / 100)) * dt
    JointScore <- max(0, min(100, JointScore + dJoint))

    # ── QoL ──
    QoL_tgt <- 1.0 - k_QoL_joint * JointScore / 100 - k_QoL_ABR * ABR_inst / 30
    QoL_tgt <- max(0.1, QoL_tgt)
    dQoL <- 0.01 * (QoL_tgt - QoL) * dt
    QoL  <- max(0.1, min(1.0, QoL + dQoL))

    # ── Record ──
    if ((step - 1) %% rec_every == 0) {
      times_rec[rec_idx] <- t_h / 24
      FVIII_rec[rec_idx] <- FVIII_total
      EMIC_rec[rec_idx]  <- EMIC_C
      FITU_rec[rec_idx]  <- FITU_C
      AT_rec[rec_idx]    <- AT_prot
      ETP_rec[rec_idx]   <- ETP
      ABR_rec[rec_idx]   <- ABR_inst
      Joint_rec[rec_idx] <- JointScore
      QoL_rec[rec_idx]   <- QoL
      Inhib_rec[rec_idx] <- Inhibitor
      rec_idx <- rec_idx + 1
    }
  }

  idx <- 1:(rec_idx - 1)
  data.frame(
    time_days   = times_rec[idx],
    FVIII       = FVIII_rec[idx],
    Emicizumab  = EMIC_rec[idx],
    Fitusiran   = FITU_rec[idx],
    AT_prot     = AT_rec[idx],
    ETP         = ETP_rec[idx],
    ABR         = ABR_rec[idx],
    JointScore  = Joint_rec[idx],
    QoL         = QoL_rec[idx],
    Inhibitor   = Inhib_rec[idx]
  )
}

# ── UI ───────────────────────────────────────────────────────────────────────

ui <- page_navbar(
  title = "Hemophilia A QSP Dashboard",
  theme = bs_theme(
    bootswatch = "darkly",
    primary    = "#0077b6",
    secondary  = "#48cae4",
    base_font  = font_google("Roboto"),
    heading_font = font_google("Roboto Condensed")
  ),
  bg = "#1a1a2e",
  fillable = TRUE,

  # ── Tab 1: Patient Profile ─────────────────────────────────────────────────
  nav_panel(
    title = "1. Patient Profile",
    icon  = icon("user-circle"),
    layout_sidebar(
      sidebar = sidebar(
        width = 320, bg = "#16213e",
        h5("Patient Parameters", class = "text-info"),
        sliderInput("BW",   "Body Weight (kg)", 20, 120, 70, step = 5),
        sliderInput("base_fviii", "Endogenous FVIII (%)",
                    0, 5, 0.5, step = 0.5),
        selectInput("severity", "HA Severity",
                    choices = c("Severe (<1%)"    = "severe",
                                "Moderate (1–5%)" = "moderate",
                                "Mild (>5%)"      = "mild")),
        sliderInput("inhibitor_bu", "Inhibitor Titer (BU/mL)", 0, 50, 0, step = 1),
        hr(),
        h5("Treatment", class = "text-info"),
        selectInput("tx_type", "Treatment Modality",
                    choices = c("No prophylaxis"       = "none",
                                "SHL-FVIII 3×/week"   = "shl",
                                "EHL-FVIII Q3-4 days" = "ehl",
                                "Emicizumab Q1W"      = "emic_q1w",
                                "Emicizumab Q4W"      = "emic_q4w",
                                "Fitusiran Q1M"       = "fitu"),
                    selected = "shl"),
        conditionalPanel(
          "input.tx_type == 'shl' || input.tx_type == 'ehl'",
          sliderInput("dose_fviii", "FVIII Dose (IU/kg)",
                      10, 80, 25, step = 5),
          sliderInput("freq_fviii", "Frequency (×/week)",
                      1, 7, 3, step = 1)
        ),
        conditionalPanel(
          "input.tx_type == 'emic_q1w'",
          sliderInput("dose_emic", "Emicizumab Maint (mg/kg/wk)",
                      0.75, 3.0, 1.5, step = 0.25)
        ),
        sliderInput("duration_d", "Simulation Duration (days)",
                    90, 3650, 365, step = 30),
        actionButton("run_sim", "Run Simulation",
                     class = "btn-primary btn-lg w-100", icon = icon("play"))
      ),
      card(
        card_header("Patient Summary"),
        layout_columns(
          col_widths = c(4, 4, 4),
          value_box("Severity", textOutput("sev_box"),
                    theme = "primary", showcase = icon("tint")),
          value_box("Inhibitor", textOutput("inhib_box"),
                    theme = "warning", showcase = icon("shield-virus")),
          value_box("Treatment", textOutput("tx_box"),
                    theme = "success", showcase = icon("pills"))
        ),
        p("Hemophilia A is an X-linked recessive bleeding disorder caused by
          deficiency of coagulation factor VIII (FVIII). Severity is classified
          by endogenous FVIII activity: severe (<1 IU/dL), moderate (1–5 IU/dL),
          or mild (>5 IU/dL). Approximately 30% of severe HA patients develop
          inhibitors (anti-FVIII antibodies), complicating replacement therapy.",
          class = "text-muted small mt-2"),
        tableOutput("patient_summary_tbl")
      )
    )
  ),

  # ── Tab 2: FVIII PK ────────────────────────────────────────────────────────
  nav_panel(
    title = "2. FVIII PK",
    icon  = icon("chart-line"),
    layout_sidebar(
      sidebar = sidebar(
        width = 240, bg = "#16213e",
        h5("PK Display Options", class = "text-info"),
        checkboxInput("log_scale", "Log Y-axis", FALSE),
        sliderInput("pk_days", "Display (days)", 1, 60, 28),
        hr(),
        helpText("FVIII half-life:"),
        tags$ul(class = "text-muted small",
          tags$li("SHL-FVIII: ~8–12 h"),
          tags$li("EHL-FVIII (Fc/PEG): ~18–19 h"),
          tags$li("Emicizumab: ~4–5 weeks")
        )
      ),
      card(
        card_header("FVIII Activity & Emicizumab Concentration vs. Time"),
        plotlyOutput("pk_plot", height = "400px")
      ),
      card(
        card_header("PK Summary Statistics"),
        tableOutput("pk_summary_tbl")
      )
    )
  ),

  # ── Tab 3: PD Core Metrics ─────────────────────────────────────────────────
  nav_panel(
    title = "3. PD Core Metrics",
    icon  = icon("heartbeat"),
    layout_columns(
      col_widths = c(6, 6),
      card(
        card_header("Thrombin Generation Potential (ETP)"),
        plotlyOutput("etp_plot", height = "320px"),
        p("ETP reflects thrombin generation capacity. Normal: ~100 (normalized).
          Severe HA: ~10–20% of normal. Prophylaxis restores ETP proportionally
          to FVIII trough levels.", class = "text-muted small")
      ),
      card(
        card_header("Antithrombin Level (Fitusiran Effect)"),
        plotlyOutput("at_plot", height = "320px"),
        p("Fitusiran (siRNA) knocks down antithrombin mRNA (t½ ~4.4 days) and
          protein (t½ ~9.6 days), reducing AT to ~20–30% of baseline. This
          augments thrombin generation independent of FVIII.", class = "text-muted small")
      )
    )
  ),

  # ── Tab 4: Bleed Risk & ABR ────────────────────────────────────────────────
  nav_panel(
    title = "4. Bleed Risk & ABR",
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
      value_box("Mean ABR", textOutput("mean_abr_vb"),
                theme = "danger",  showcase = icon("tint")),
      value_box("Final Joint Score", textOutput("joint_vb"),
                theme = "warning", showcase = icon("bone")),
      value_box("Final QoL (EQ-5D)", textOutput("qol_vb"),
                theme = "success", showcase = icon("smile"))
    )
  ),

  # ── Tab 5: Scenario Comparison ─────────────────────────────────────────────
  nav_panel(
    title = "5. Scenario Comparison",
    icon  = icon("chart-bar"),
    layout_sidebar(
      sidebar = sidebar(
        width = 280, bg = "#16213e",
        h5("Compare Scenarios", class = "text-info"),
        checkboxGroupInput("compare_scenarios", "Include Scenarios",
          choices = c(
            "No Prophylaxis"       = "none",
            "SHL-FVIII 3×/wk"    = "shl",
            "EHL-FVIII Q3-4d"    = "ehl",
            "Emicizumab Q1W"     = "emic_q1w",
            "Emicizumab Q4W"     = "emic_q4w",
            "Fitusiran Q1M"      = "fitu"
          ),
          selected = c("none", "shl", "emic_q1w", "fitu")
        ),
        sliderInput("cmp_duration", "Duration (years)", 1, 10, 2),
        actionButton("run_comparison", "Compare",
                     class = "btn-info w-100", icon = icon("balance-scale"))
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

  # ── Tab 6: Biomarkers ──────────────────────────────────────────────────────
  nav_panel(
    title = "6. Biomarkers",
    icon  = icon("microscope"),
    layout_columns(
      col_widths = c(6, 6),
      card(
        card_header("FVIII Activity vs. ETP Relationship"),
        plotlyOutput("fviii_etp_scatter", height = "300px"),
        p("Non-linear relationship between FVIII activity and thrombin
          generation (Hill equation). EC50 for ETP ~5 IU/dL.",
          class = "text-muted small")
      ),
      card(
        card_header("Inhibitor Titer Dynamics"),
        plotlyOutput("inhibitor_plot", height = "300px"),
        p("Anti-FVIII antibody (inhibitor) titer dynamics.
          High responders (>5 BU/mL) require alternative therapies
          (emicizumab, fitusiran) bypassing FVIII.", class = "text-muted small")
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

  # Reactive simulation
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Simulating...", value = 0.5, {
      use_shl  <- input$tx_type %in% c("shl", "none")
      use_emic <- input$tx_type %in% c("emic_q1w", "emic_q4w")
      use_fitu <- input$tx_type == "fitu"
      freq_fviii_val <- if (input$tx_type == "ehl") 2 else input$freq_fviii
      dose_emic_val  <- if (input$tx_type == "emic_q4w") 6.0 else input$dose_emic

      simulate_ha(
        BW          = input$BW,
        base_fviii  = input$base_fviii,
        use_shl     = use_shl,
        use_emic    = use_emic,
        use_fitu    = use_fitu,
        dose_fviii  = input$dose_fviii,
        freq_fviii  = freq_fviii_val,
        dose_emic   = dose_emic_val,
        dose_fitu   = 80,
        inhibitor_bu = input$inhibitor_bu,
        duration_d  = input$duration_d
      )
    })
  }, ignoreNULL = FALSE)

  # Default simulation on load
  observeEvent(TRUE, {
    shinyjs::runjs("document.getElementById('run_sim').click();")
  }, once = TRUE, ignoreInit = FALSE)

  # ── Tab 1 outputs ──────────────────────────────────────────────────────────
  output$sev_box <- renderText({
    sev <- switch(input$severity,
                  severe   = "Severe (<1%)",
                  moderate = "Moderate (1–5%)",
                  mild     = "Mild (>5%)")
    sev
  })
  output$inhib_box <- renderText({
    if (input$inhibitor_bu == 0) "Negative" else
      paste0(input$inhibitor_bu, " BU/mL")
  })
  output$tx_box <- renderText({
    switch(input$tx_type,
           none     = "None (on-demand)",
           shl      = "SHL-FVIII",
           ehl      = "EHL-FVIII",
           emic_q1w = "Emicizumab Q1W",
           emic_q4w = "Emicizumab Q4W",
           fitu     = "Fitusiran Q1M")
  })
  output$patient_summary_tbl <- renderTable({
    data.frame(
      Parameter = c("Body Weight", "FVIII Severity", "Inhibitor Status",
                    "Simulation Duration"),
      Value = c(paste0(input$BW, " kg"),
                input$severity,
                if (input$inhibitor_bu == 0) "None"
                else paste0(input$inhibitor_bu, " BU/mL"),
                paste0(input$duration_d, " days"))
    )
  })

  # ── Tab 2 PK ───────────────────────────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    df <- sim_data() %>% filter(time_days <= input$pk_days)
    p  <- plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~FVIII, name = "FVIII Activity (IU/dL)",
                line = list(color = "#0077b6", width = 2)) %>%
      add_lines(y = ~Emicizumab * 100, name = "Emicizumab (mg/L ×100)",
                line = list(color = "#90e0ef", dash = "dash")) %>%
      add_lines(y = ~Fitusiran * 1000, name = "Fitusiran (mg/L ×1000)",
                line = list(color = "#ffd166", dash = "dot")) %>%
      layout(
        xaxis = list(title = "Time (days)", color = "white"),
        yaxis = list(title = "Concentration / Activity",
                     type  = if (input$log_scale) "log" else "linear",
                     color = "white"),
        shapes = list(
          list(type = "line", x0 = 0, x1 = input$pk_days,
               y0 = 1, y1 = 1,
               line = list(color = "red", dash = "dot", width = 1.5)),
          list(type = "line", x0 = 0, x1 = input$pk_days,
               y0 = 15, y1 = 15,
               line = list(color = "orange", dash = "dot", width = 1.5))
        ),
        annotations = list(
          list(x = 1, y = 2,  text = "1% trough", showarrow = FALSE,
               font = list(color = "red",    size = 10)),
          list(x = 1, y = 17, text = "15% optimal", showarrow = FALSE,
               font = list(color = "orange", size = 10))
        ),
        paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e",
        font = list(color = "white"), legend = list(font = list(color = "white"))
      )
    p
  })

  output$pk_summary_tbl <- renderTable({
    df <- sim_data() %>% filter(time_days > 0)
    data.frame(
      Metric = c("Mean FVIII Activity (IU/dL)", "Min FVIII (trough)",
                 "Max FVIII (peak)", "Mean Emicizumab (mg/L)",
                 "Mean Fitusiran (mg/L)"),
      Value  = round(c(mean(df$FVIII), min(df$FVIII),
                       max(df$FVIII), mean(df$Emicizumab),
                       mean(df$Fitusiran)), 2)
    )
  })

  # ── Tab 3 PD ───────────────────────────────────────────────────────────────
  output$etp_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~ETP, name = "Thrombin ETP",
                line = list(color = "#06d6a0", width = 2)) %>%
      layout(
        xaxis = list(title = "Time (days)", color = "white"),
        yaxis = list(title = "ETP (normalized)", range = c(0, 110),
                     color = "white"),
        shapes = list(
          list(type = "line", x0 = 0, x1 = max(df$time_days),
               y0 = 80, y1 = 80,
               line = list(color = "green", dash = "dot"))),
        annotations = list(
          list(x = 5, y = 85, text = "Normal ETP", showarrow = FALSE,
               font = list(color = "green", size = 10))),
        paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e",
        font = list(color = "white")
      )
  })

  output$at_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~AT_prot * 100, name = "AT Protein (%)",
                line = list(color = "#f8961e", width = 2)) %>%
      layout(
        xaxis = list(title = "Time (days)", color = "white"),
        yaxis = list(title = "Antithrombin Level (% of baseline)",
                     range = c(0, 110), color = "white"),
        shapes = list(
          list(type = "line", x0 = 0, x1 = max(df$time_days),
               y0 = 25, y1 = 25,
               line = list(color = "red", dash = "dot"))),
        annotations = list(
          list(x = 5, y = 30, text = "Fitusiran target ~25%",
               showarrow = FALSE, font = list(color = "red", size = 10))),
        paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e",
        font = list(color = "white")
      )
  })

  # ── Tab 4 Bleed ────────────────────────────────────────────────────────────
  output$abr_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~ABR, name = "Instantaneous ABR",
                line = list(color = "#e63946", width = 2)) %>%
      layout(
        xaxis = list(title = "Time (days)", color = "white"),
        yaxis = list(title = "Annual Bleed Rate", color = "white"),
        shapes = list(
          list(type = "line", x0 = 0, x1 = max(df$time_days),
               y0 = 3, y1 = 3,
               line = list(color = "green", dash = "dot"))),
        paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e",
        font = list(color = "white")
      )
  })

  output$joint_qol_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~JointScore, name = "Joint Score (0-100)",
                line = list(color = "#8338ec", width = 2)) %>%
      add_lines(y = ~QoL * 100, name = "QoL × 100 (EQ-5D)",
                line = list(color = "#06d6a0", width = 2, dash = "dash")) %>%
      layout(
        xaxis = list(title = "Time (days)", color = "white"),
        yaxis = list(title = "Score", range = c(0, 105), color = "white"),
        paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e",
        font = list(color = "white"),
        legend = list(font = list(color = "white"))
      )
  })

  output$mean_abr_vb <- renderText({
    df <- sim_data() %>% filter(time_days > 7)
    round(mean(df$ABR), 1)
  })
  output$joint_vb <- renderText({
    df <- sim_data()
    round(tail(df$JointScore, 1), 1)
  })
  output$qol_vb <- renderText({
    df <- sim_data()
    round(tail(df$QoL, 3), 2)
  })

  # ── Tab 5 Scenario Comparison ──────────────────────────────────────────────
  cmp_data <- eventReactive(input$run_comparison, {
    scenarios <- input$compare_scenarios
    if (length(scenarios) == 0) scenarios <- c("none", "shl")
    dur <- input$cmp_duration * 365

    scenario_map <- list(
      none     = list(use_shl = TRUE,  use_emic = FALSE, use_fitu = FALSE,
                       dose_fviii = 25, freq_fviii = 3, dose_emic = 1.5,
                       label = "No Prophylaxis"),
      shl      = list(use_shl = TRUE,  use_emic = FALSE, use_fitu = FALSE,
                       dose_fviii = 25, freq_fviii = 3, dose_emic = 1.5,
                       label = "SHL-FVIII 3×/wk"),
      ehl      = list(use_shl = FALSE, use_emic = FALSE, use_fitu = FALSE,
                       dose_fviii = 50, freq_fviii = 2, dose_emic = 1.5,
                       label = "EHL-FVIII Q3-4d"),
      emic_q1w = list(use_shl = FALSE, use_emic = TRUE,  use_fitu = FALSE,
                       dose_fviii = 25, freq_fviii = 3, dose_emic = 1.5,
                       label = "Emicizumab Q1W"),
      emic_q4w = list(use_shl = FALSE, use_emic = TRUE,  use_fitu = FALSE,
                       dose_fviii = 25, freq_fviii = 3, dose_emic = 6.0,
                       label = "Emicizumab Q4W"),
      fitu     = list(use_shl = FALSE, use_emic = FALSE, use_fitu = TRUE,
                       dose_fviii = 25, freq_fviii = 3, dose_emic = 1.5,
                       label = "Fitusiran Q1M")
    )

    withProgress(message = "Comparing scenarios...", {
      results <- lapply(scenarios, function(s) {
        sm <- scenario_map[[s]]
        df <- simulate_ha(BW = input$BW,
                          base_fviii = input$base_fviii,
                          use_shl  = sm$use_shl,
                          use_emic = sm$use_emic,
                          use_fitu = sm$use_fitu,
                          dose_fviii = sm$dose_fviii,
                          freq_fviii = sm$freq_fviii,
                          dose_emic  = sm$dose_emic,
                          duration_d = dur)
        df$scenario <- sm$label
        df
      })
      bind_rows(results)
    })
  })

  output$cmp_abr <- renderPlotly({
    df <- cmp_data()
    colors <- c("#e63946","#457b9d","#2a9d8f","#e9c46a","#f4a261","#264653")
    scenarios <- unique(df$scenario)
    p <- plot_ly()
    for (i in seq_along(scenarios)) {
      d <- df %>% filter(scenario == scenarios[i])
      p <- add_lines(p, x = d$time_days / 365, y = d$ABR,
                     name = scenarios[i],
                     line = list(color = colors[i], width = 2))
    }
    layout(p,
           xaxis = list(title = "Time (years)", color = "white"),
           yaxis = list(title = "Annual Bleed Rate", color = "white"),
           shapes = list(list(type = "line", x0 = 0,
                               x1 = max(df$time_days) / 365,
                               y0 = 3, y1 = 3,
                               line = list(color = "green", dash = "dot"))),
           paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e",
           font = list(color = "white"),
           legend = list(font = list(color = "white")))
  })

  output$cmp_joint <- renderPlotly({
    df <- cmp_data()
    colors <- c("#e63946","#457b9d","#2a9d8f","#e9c46a","#f4a261","#264653")
    scenarios <- unique(df$scenario)
    p <- plot_ly()
    for (i in seq_along(scenarios)) {
      d <- df %>% filter(scenario == scenarios[i])
      p <- add_lines(p, x = d$time_days / 365, y = d$JointScore,
                     name = scenarios[i],
                     line = list(color = colors[i], width = 2))
    }
    layout(p,
           xaxis = list(title = "Time (years)", color = "white"),
           yaxis = list(title = "Pettersson Joint Score", color = "white"),
           paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e",
           font = list(color = "white"),
           legend = list(font = list(color = "white")))
  })

  # ── Tab 6 Biomarkers ───────────────────────────────────────────────────────
  output$fviii_etp_scatter <- renderPlotly({
    fviii_seq <- seq(0, 100, by = 0.5)
    etp_calc <- 100 * fviii_seq^0.8 / (5^0.8 + fviii_seq^0.8)
    df_curve <- data.frame(FVIII = fviii_seq, ETP = etp_calc)
    df <- sim_data()

    plot_ly() %>%
      add_lines(data = df_curve, x = ~FVIII, y = ~ETP,
                name = "Theoretical Curve",
                line = list(color = "#90e0ef", width = 2)) %>%
      add_markers(data = df %>% filter(row_number() %% 10 == 1),
                  x = ~FVIII, y = ~ETP,
                  name = "Simulation",
                  marker = list(color = "#06d6a0", size = 5, opacity = 0.5)) %>%
      layout(
        xaxis = list(title = "FVIII Activity (IU/dL)", color = "white"),
        yaxis = list(title = "Thrombin ETP (normalized)", color = "white"),
        paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e",
        font = list(color = "white"),
        legend = list(font = list(color = "white"))
      )
  })

  output$inhibitor_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(df, x = ~time_days) %>%
      add_lines(y = ~Inhibitor, name = "Inhibitor Titer (BU/mL)",
                line = list(color = "#ff6b6b", width = 2)) %>%
      layout(
        xaxis = list(title = "Time (days)", color = "white"),
        yaxis = list(title = "Inhibitor Titer (BU/mL)", color = "white"),
        shapes = list(
          list(type = "line", x0 = 0, x1 = max(df$time_days),
               y0 = 5, y1 = 5,
               line = list(color = "orange", dash = "dot"))),
        annotations = list(
          list(x = 5, y = 6, text = "High responder threshold (5 BU/mL)",
               showarrow = FALSE, font = list(color = "orange", size = 10))),
        paper_bgcolor = "#1a1a2e", plot_bgcolor = "#16213e",
        font = list(color = "white")
      )
  })

  output$biomarker_summary <- renderTable({
    df <- sim_data() %>% filter(time_days > 7)
    data.frame(
      Biomarker = c("Mean FVIII Activity (IU/dL)",
                    "FVIII Trough < 1% Time (%)",
                    "Mean ETP (normalized)",
                    "Mean AT Level (%)",
                    "Mean ABR",
                    "Final Joint Score",
                    "Final QoL (EQ-5D)",
                    "Final Inhibitor (BU/mL)"),
      Value = round(c(
        mean(df$FVIII),
        100 * mean(df$FVIII < 1),
        mean(df$ETP),
        mean(df$AT_prot) * 100,
        mean(df$ABR),
        tail(df$JointScore, 1),
        tail(df$QoL, 1),
        tail(df$Inhibitor, 1)
      ), 2),
      Target = c(">1 (trough), >15 (optimal)",
                 "< 20%",
                 "> 80",
                 "20–30% (Fitusiran) / 100 (others)",
                 "< 3",
                 "< 20",
                 "> 0.80",
                 "< 0.6 (negative)")
    )
  }, striped = TRUE, hover = TRUE)
}

shinyApp(ui = ui, server = server)
