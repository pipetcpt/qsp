## ============================================================
## IBS QSP Model — Interactive Shiny Dashboard
## Irritable Bowel Syndrome: Brain-Gut Axis QSP Model
## 6 Tabs: Patient Profile · PK · Biomarkers ·
##         Symptom Scores · Scenario Comparison · Biomarker Tracker
## ============================================================
## Dependencies: shiny, shinydashboard, plotly, deSolve, dplyr,
##               tidyr, DT, shinyWidgets
## ============================================================

library(shiny)
library(shinydashboard)
library(plotly)
library(deSolve)
library(dplyr)
library(tidyr)
library(DT)
library(shinyWidgets)

## ============================================================
## 1. ODE System (deSolve version — standalone for Shiny)
## ============================================================

ibs_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    # ---- Drug PK ----
    Alo_RO  <- ifelse(Emax_alo * Cp1 / (IC50_alo + Cp1 + 1e-9) > 0,
                      Emax_alo * Cp1 / (IC50_alo + Cp1 + 1e-9), 0)
    Pru_RO  <- ifelse(Emax_pru * Cp3 / (EC50_pru + Cp3 + 1e-9) > 0,
                      Emax_pru * Cp3 / (EC50_pru + Cp3 + 1e-9), 0)

    # ---- Stress / CRF / Cortisol ----
    dSTRESS <- -0.02 * (STRESS - stress_base)
    dCRF    <- k_crf_syn * STRESS - k_crf_deg * CRF
    dCORT   <- k_cort_syn * CRF - k_cort_deg * CORT

    # ---- 5-HT dynamics ----
    syn_5HT  <- k_5HT_syn * base_5HT_adj *
                (1 + 0.3 * CORT / 15) * (1 + 0.2 * MAST_ACT / 0.25)
    SERT_eff <- SERT_OCC + SERT_inh_ami
    if (SERT_eff > 0.95) SERT_eff <- 0.95
    reup_5HT <- k_5HT_reup * GUT_5HT * (1 - SERT_eff)
    dGUT_5HT <- syn_5HT - reup_5HT - k_5HT_deg * GUT_5HT

    dSERT_OCC <- 0.5 * (Alo_RO - SERT_OCC)

    # ---- Mast cells ----
    mast_stim  <- k_mast_activ * (0.4 * CRF / 0.8 +
                                   0.3 * MICROB / 0.35 +
                                   0.3 * max(1 - BARRIER, 0))
    mast_inhib <- k_mast_rest * (1 + 0.4 * SCFA / 10)
    dMAST_ACT  <- mast_stim * (1 - MAST_ACT) - mast_inhib * MAST_ACT

    # ---- Inflammation ----
    inflam_in  <- k_inflam_on * (0.5 * MAST_ACT + 0.3 * max(1 - BARRIER, 0) +
                                  0.2 * MICROB)
    inflam_out <- k_inflam_off * INFLAM * (1 + 0.5 * SCFA / 10)
    dINFLAM    <- max(inflam_in - inflam_out, -INFLAM)

    # ---- Barrier ----
    b_repair <- k_barrier_rep * SCFA / (5 + SCFA) * (1 - BARRIER)
    b_break  <- k_barrier_brk * INFLAM * (1 + 0.3 * CORT / 15) * BARRIER
    dBARRIER  <- b_repair - b_break

    # ---- Microbiome ----
    dysb_inc <- k_dysb_form * (0.5 * INFLAM + 0.5 * max(1 - BARRIER, 0))
    dysb_dec <- k_dysb_res * (1 + probiotic_eff * 2 + rifax_effect * 3)
    dMICROB   <- dysb_inc * (1 - MICROB) - dysb_dec * MICROB

    # ---- SCFA ----
    dSCFA <- SCFA_prod_base * (1 - MICROB) - k_SCFA_util * SCFA

    # ---- Visceral Sensitivity ----
    sens_drive <- k_vis_sens_on * (0.4 * MAST_ACT + 0.4 * INFLAM + 0.2 * DRG_ACT)
    alo_inh    <- Alo_RO * 0.5
    ami_inh    <- pain_inh_ami * 0.6
    sens_res   <- k_vis_sens_off * VIS_SENS * (1 + alo_inh + ami_inh)
    dVIS_SENS  <- sens_drive * (1 - VIS_SENS) - sens_res

    # ---- DRG Activation ----
    drg_drive <- k_DRG_on * (0.4 * GUT_5HT / 1.5 + 0.3 * MAST_ACT + 0.3 * INFLAM)
    drg_block <- k_DRG_off * DRG_ACT * (1 + Alo_RO * 0.7 + pain_inh_ami * 0.4)
    dDRG_ACT   <- drg_drive * (1 - DRG_ACT) - drg_block

    # ---- Motility ----
    sHT_motil <- ifelse(IBS_TYPE == 1, 0.3 * (GUT_5HT / 1.5 - 1), 0)
    motil_tgt  <- motil_base *
                  (1 + Pru_RO * 0.4) *
                  (1 + sHT_motil) *
                  (1 - 0.15 * MICROB)
    dMOTIL     <- k_motil_adapt * (motil_tgt - MOTIL)

    # ---- Pain ----
    pain_drv   <- 8 * (0.45 * DRG_ACT + 0.35 * VIS_SENS + 0.2 * STRESS)
    pain_inh   <- Alo_RO * 0.5 + pain_inh_ami
    pain_tgt   <- pain_drv * (1 - pain_inh * 0.6)
    dPAIN      <- 0.5 * (pain_tgt - PAIN)

    # ---- Bloating ----
    bloat_tgt <- 7 * (0.6 * MICROB + 0.4 * INFLAM) * (1 - rifax_effect * 0.6)
    dBLOAT    <- 0.4 * (bloat_tgt - BLOAT)

    # ---- Stool BSFS ----
    stool_tgt <- if (IBS_TYPE == 1) {
      4.5 + 1.5 * (MOTIL - 1) + 0.5 * (1 - BARRIER) + Pru_RO * 1.5 - Alo_RO * 1.0
    } else {
      2.5 - 1.5 * (1 - MOTIL) + 0.3 * SCFA / 10 + Pru_RO * 1.5 - Alo_RO * 1.0
    }
    stool_tgt <- max(1, min(7, stool_tgt))
    dSTOOL    <- 0.5 * (stool_tgt - STOOL)

    # ---- IBS-SSS ----
    stool_dev <- max(0, ifelse(IBS_TYPE == 1, STOOL - 4, 4 - STOOL))
    sss_tgt   <- 15 * PAIN + 10 * BLOAT + 8 * stool_dev * 10 + 12 * DRG_ACT * 10
    sss_tgt   <- max(0, min(500, sss_tgt))
    dIBS_SSS  <- 0.3 * (sss_tgt - IBS_SSS)

    # ---- Drug PK ----
    dCp1 <- -(CL1 / V1 + Q1 / V1) * Cp1 + (Q1 / V2) * Cp2
    dCp2 <- (Q1 / V1) * Cp1 - (Q1 / V2) * Cp2
    dCp3 <- -(CL3 / V3) * Cp3

    list(c(dSTRESS, dCRF, dCORT, dGUT_5HT, dSERT_OCC,
           dMAST_ACT, dINFLAM, dBARRIER, dMICROB, dSCFA,
           dVIS_SENS, dMOTIL, dPAIN, dBLOAT, dSTOOL, dIBS_SSS,
           dDRG_ACT, dCp1, dCp2, dCp3),
         Alo_RO = Alo_RO, Pru_RO = Pru_RO,
         DEFEC_FREQ = 1 + 2.5 * MOTIL,
         QoL = max(0, min(100, 100 * (1 - IBS_SSS / 500))),
         SSS_RED = (220 - IBS_SSS) / 220 * 100)
  })
}

## ============================================================
## 2. Simulation Helper
## ============================================================

run_ibs_sim <- function(
    ibs_type     = 1,      # 1=IBS-D, 2=IBS-C, 3=IBS-M
    stress_level = 0.45,
    scenario     = 1,
    weeks        = 26,
    dose_alo     = 0,      # μg alosetron
    dose_pru     = 0,      # μg prucalopride
    sert_inh_ami = 0,
    pain_inh_ami = 0,
    rifax_eff    = 0,
    prob_eff     = 0
) {
  base_5HT_adj <- ifelse(ibs_type == 1, 1.8, 0.65)

  y0 <- c(
    STRESS   = stress_level,
    CRF      = 0.8,
    CORT     = 15.0,
    GUT_5HT  = base_5HT_adj,
    SERT_OCC = 0.0,
    MAST_ACT = 0.28,
    INFLAM   = 0.22,
    BARRIER  = 0.82,
    MICROB   = 0.38,
    SCFA     = 9.5,
    VIS_SENS = 0.37,
    MOTIL    = ifelse(ibs_type == 1, 1.25, 0.75),
    PAIN     = ifelse(ibs_type == 1, 4.8, 3.5),
    BLOAT    = 3.5,
    STOOL    = ifelse(ibs_type == 1, 5.5, 2.5),
    IBS_SSS  = ifelse(ibs_type == 1, 245, 195),
    DRG_ACT  = 0.35,
    Cp1      = 0.0,
    Cp2      = 0.0,
    Cp3      = 0.0
  )

  parms <- list(
    IBS_TYPE       = ibs_type,
    stress_base    = stress_level,
    base_5HT_adj   = base_5HT_adj,
    k_crf_syn      = 0.8,   k_crf_deg     = 0.5,
    k_cort_syn     = 1.2,   k_cort_deg    = 0.35,
    k_5HT_syn      = 0.6,   k_5HT_reup    = 0.8,
    k_5HT_deg      = 0.3,   SERT_inh_ami  = sert_inh_ami,
    k_mast_activ   = 0.4,   k_mast_rest   = 0.25,
    k_inflam_on    = 0.3,   k_inflam_off  = 0.2,
    k_barrier_rep  = 0.15,  k_barrier_brk = 0.25,
    k_dysb_form    = 0.15,  k_dysb_res    = 0.10,
    SCFA_prod_base = 8.0,   k_SCFA_util   = 0.4,
    k_vis_sens_on  = 0.3,   k_vis_sens_off = 0.15,
    motil_base     = 1.0,   k_motil_adapt  = 0.4,
    k_DRG_on       = 0.5,   k_DRG_off      = 0.35,
    pain_inh_ami   = pain_inh_ami,
    rifax_effect   = rifax_eff,
    probiotic_eff  = prob_eff,
    CL1 = 30, V1 = 65, Q1 = 15, V2 = 55,
    IC50_alo = 1.2, Emax_alo = 0.90,
    CL3 = 7.0, V3 = 200,
    EC50_pru = 2.5, Emax_pru = 0.85
  )

  sim_hours <- weeks * 7 * 24
  times     <- seq(0, sim_hours, by = 4)

  # Dosing events
  dose_times_alo <- seq(0, sim_hours, by = 12)
  dose_times_pru <- seq(0, sim_hours, by = 24)
  dose_amt_alo   <- dose_alo / 65   # convert to ng/mL (Vd=65L)
  dose_amt_pru   <- dose_pru / 200  # convert to ng/mL (Vd=200L)

  # Run deSolve with event dosing
  events_alo <- data.frame(var  = "Cp1",
                            time = dose_times_alo,
                            value = dose_amt_alo,
                            method = "add")
  events_pru <- data.frame(var  = "Cp3",
                            time = dose_times_pru,
                            value = dose_amt_pru,
                            method = "add")
  all_events <- if (dose_alo > 0 && dose_pru > 0) {
    bind_rows(events_alo, events_pru) %>% arrange(time)
  } else if (dose_alo > 0) {
    events_alo
  } else if (dose_pru > 0) {
    events_pru
  } else {
    NULL
  }

  if (!is.null(all_events)) {
    out <- deSolve::ode(y = y0, times = times, func = ibs_ode,
                        parms = parms,
                        events = list(data = all_events),
                        method = "lsoda")
  } else {
    out <- deSolve::ode(y = y0, times = times, func = ibs_ode,
                        parms = parms, method = "lsoda")
  }

  as.data.frame(out) %>%
    mutate(time_days = time / 24)
}

## ============================================================
## 3. UI Definition
## ============================================================

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = "IBS QSP Dashboard",
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      menuItem("① Patient Profile",      tabName = "profile",    icon = icon("user-md")),
      menuItem("② Drug PK",              tabName = "pk",         icon = icon("pills")),
      menuItem("③ Biomarkers",           tabName = "biomarkers", icon = icon("flask")),
      menuItem("④ Symptom Scores",       tabName = "symptoms",   icon = icon("chart-line")),
      menuItem("⑤ Scenario Comparison",  tabName = "scenarios",  icon = icon("layer-group")),
      menuItem("⑥ Biomarker Tracker",    tabName = "tracker",    icon = icon("heartbeat"))
    ),

    hr(),
    h5("  Patient & Disease Settings", style = "color:#BDC3C7; margin-left:10px;"),
    selectInput("ibs_type", "IBS Subtype:",
                choices = c("IBS-D (Diarrhea)" = 1,
                            "IBS-C (Constipation)" = 2,
                            "IBS-M (Mixed)" = 3),
                selected = 1),
    sliderInput("stress", "Stress Level (0–1):",
                min = 0.1, max = 0.9, value = 0.45, step = 0.05),
    sliderInput("weeks", "Simulation Duration (weeks):",
                min = 4, max = 52, value = 24, step = 4),

    hr(),
    h5("  Drug Treatment", style = "color:#BDC3C7; margin-left:10px;"),
    selectInput("drug_choice", "Select Drug:",
                choices = c("No Drug (Placebo)"              = "none",
                            "Alosetron 0.5 mg BID (IBS-D)"   = "alosetron",
                            "Ondansetron 4 mg TID"            = "ondansetron",
                            "Prucalopride 2 mg QD (IBS-C)"   = "prucalopride",
                            "Amitriptyline 25 mg QN"          = "amitriptyline",
                            "Rifaximin + Probiotics"          = "rifaximin"),
                selected = "none"),
    actionButton("run_sim", "Run Simulation",
                 icon = icon("play"), class = "btn-success",
                 style = "margin:10px; width:220px;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #F0F4F8; }
      .small-box { border-radius: 8px; }
      .box { border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); }
    "))),

    tabItems(

      ## ----------------------------------------------------------
      ## TAB 1: Patient Profile
      ## ----------------------------------------------------------
      tabItem(tabName = "profile",
        fluidRow(
          valueBoxOutput("vb_pain",     width = 3),
          valueBoxOutput("vb_sss",      width = 3),
          valueBoxOutput("vb_stool",    width = 3),
          valueBoxOutput("vb_qol",      width = 3)
        ),
        fluidRow(
          box(title = "Disease Overview — Irritable Bowel Syndrome (IBS)",
              width = 6, solidHeader = TRUE, status = "primary",
              p("IBS is a functional gastrointestinal disorder affecting ~11% of the global population,
                characterized by chronic abdominal pain, altered bowel habits, and bloating without
                structural abnormality."),
              h4("Diagnostic Criteria (Rome IV)"),
              tags$ul(
                tags$li("Recurrent abdominal pain ≥1 day/week (average) in the last 3 months"),
                tags$li("Associated with ≥2 of: (1) related to defecation, (2) change in stool frequency,
                         (3) change in stool form/appearance"),
                tags$li("Onset ≥6 months before diagnosis")
              ),
              h4("Subtypes"),
              tags$ul(
                tags$li(strong("IBS-D:"), " Diarrhea-predominant (BSFS 6-7 >25% time)"),
                tags$li(strong("IBS-C:"), " Constipation-predominant (BSFS 1-2 >25% time)"),
                tags$li(strong("IBS-M:"), " Mixed bowel habits"),
                tags$li(strong("IBS-U:"), " Unclassified")
              )
          ),
          box(title = "IBS-SSS Calculator", width = 6,
              solidHeader = TRUE, status = "warning",
              p("IBS Symptom Severity Score (Francis 1997) — 0 to 500 points"),
              tags$table(class = "table table-striped",
                tags$thead(tags$tr(
                  tags$th("Category"), tags$th("Score Range"), tags$th("Interpretation")
                )),
                tags$tbody(
                  tags$tr(tags$td("Mild"),     tags$td("75–175"),  tags$td("Mild symptoms")),
                  tags$tr(tags$td("Moderate"), tags$td("175–300"), tags$td("Significant impact")),
                  tags$tr(tags$td("Severe"),   tags$td(">300"),    tags$td("Substantial disability")),
                  tags$tr(tags$td("Remission"), tags$td("<75"),    tags$td("Minimal symptoms"))
                )
              ),
              hr(),
              h5("Clinical Responder Definition"),
              p("≥50 point reduction from baseline IBS-SSS OR IBS-SSS < 175 at Week 12")
          )
        ),
        fluidRow(
          box(title = "Natural History vs. Treated: IBS-SSS Over Time",
              width = 12, plotlyOutput("plot_profile_sss", height = "350px"))
        )
      ),

      ## ----------------------------------------------------------
      ## TAB 2: Drug PK
      ## ----------------------------------------------------------
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Drug PK Parameters", width = 4, solidHeader = TRUE, status = "info",
              DTOutput("tbl_pk_params")
          ),
          box(title = "Plasma Concentration — Alosetron (Cp1, ng/mL)",
              width = 8, plotlyOutput("plot_pk_alo", height = "300px"))
        ),
        fluidRow(
          box(title = "5-HT3 Receptor Occupancy (Alosetron)",
              width = 6, plotlyOutput("plot_pk_ro1", height = "280px")),
          box(title = "5-HT4 Receptor Occupancy (Prucalopride)",
              width = 6, plotlyOutput("plot_pk_ro2", height = "280px"))
        )
      ),

      ## ----------------------------------------------------------
      ## TAB 3: Biomarkers
      ## ----------------------------------------------------------
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "Gut Serotonin (5-HT, nmol/g tissue)",
              width = 6, plotlyOutput("plot_bio_5ht", height = "280px")),
          box(title = "Mast Cell Activation Index (0–1)",
              width = 6, plotlyOutput("plot_bio_mast", height = "280px"))
        ),
        fluidRow(
          box(title = "Mucosal Inflammatory Index (0–1)",
              width = 6, plotlyOutput("plot_bio_inflam", height = "280px")),
          box(title = "Epithelial Barrier Integrity (0–1)",
              width = 6, plotlyOutput("plot_bio_barrier", height = "280px"))
        ),
        fluidRow(
          box(title = "Gut Microbiome Dysbiosis Index (0–1)",
              width = 6, plotlyOutput("plot_bio_microb", height = "280px")),
          box(title = "Short-Chain Fatty Acids (mmol/L)",
              width = 6, plotlyOutput("plot_bio_scfa", height = "280px"))
        )
      ),

      ## ----------------------------------------------------------
      ## TAB 4: Symptom Scores
      ## ----------------------------------------------------------
      tabItem(tabName = "symptoms",
        fluidRow(
          box(title = "Abdominal Pain (NRS 0–10)",
              width = 6, plotlyOutput("plot_sym_pain", height = "280px")),
          box(title = "Bloating/Distension Score (0–10)",
              width = 6, plotlyOutput("plot_sym_bloat", height = "280px"))
        ),
        fluidRow(
          box(title = "Bristol Stool Form Scale (1–7)",
              width = 6, plotlyOutput("plot_sym_stool", height = "280px")),
          box(title = "IBS-SSS Composite (0–500)",
              width = 6, plotlyOutput("plot_sym_sss", height = "280px"))
        ),
        fluidRow(
          box(title = "Visceral Hypersensitivity Index (0–1)",
              width = 6, plotlyOutput("plot_sym_vis", height = "280px")),
          box(title = "Gut Motility Index (1=normal)",
              width = 6, plotlyOutput("plot_sym_motil", height = "280px"))
        )
      ),

      ## ----------------------------------------------------------
      ## TAB 5: Scenario Comparison
      ## ----------------------------------------------------------
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "All Scenarios — IBS-SSS Over 26 Weeks",
              width = 12, plotlyOutput("plot_scen_sss", height = "400px"))
        ),
        fluidRow(
          box(title = "Pain Scores — All Scenarios",
              width = 6, plotlyOutput("plot_scen_pain", height = "300px")),
          box(title = "Microbiome Dysbiosis — All Scenarios",
              width = 6, plotlyOutput("plot_scen_microb", height = "300px"))
        ),
        fluidRow(
          box(title = "Scenario Summary Table (Week 12)", width = 12,
              DTOutput("tbl_scen_summary"))
        )
      ),

      ## ----------------------------------------------------------
      ## TAB 6: Biomarker Tracker
      ## ----------------------------------------------------------
      tabItem(tabName = "tracker",
        fluidRow(
          box(title = "Quality of Life Trajectory (0–100)",
              width = 6, plotlyOutput("plot_tr_qol", height = "280px")),
          box(title = "DRG Activation Index (0–1)",
              width = 6, plotlyOutput("plot_tr_drg", height = "280px"))
        ),
        fluidRow(
          box(title = "Serum Chromogranin A Proxy (EC cell marker)",
              width = 6, plotlyOutput("plot_tr_cga", height = "280px")),
          box(title = "Rectal Sensitivity Threshold (barostat mmHg)",
              width = 6, plotlyOutput("plot_tr_rectal", height = "280px"))
        ),
        fluidRow(
          box(title = "Download Simulation Data", width = 4,
              downloadButton("dl_data", "Download CSV"),
              br(), br(),
              p("Downloads the full time-course simulation data for the selected scenario.")),
          box(title = "Biomarker Summary at Week 4/8/12/24", width = 8,
              DTOutput("tbl_biomarker_summary"))
        )
      )
    )  # end tabItems
  )    # end dashboardBody
)      # end dashboardPage

## ============================================================
## 4. Server Logic
## ============================================================

server <- function(input, output, session) {

  ## Reactive: get drug parameters from selection
  drug_params <- reactive({
    switch(input$drug_choice,
      "alosetron"    = list(dose_alo=500, dose_pru=0, sert_inh_ami=0,
                            pain_inh_ami=0, rifax_eff=0, prob_eff=0),
      "ondansetron"  = list(dose_alo=250, dose_pru=0, sert_inh_ami=0,
                            pain_inh_ami=0, rifax_eff=0, prob_eff=0),
      "prucalopride" = list(dose_alo=0, dose_pru=2000, sert_inh_ami=0,
                            pain_inh_ami=0, rifax_eff=0, prob_eff=0),
      "amitriptyline"= list(dose_alo=0, dose_pru=0, sert_inh_ami=0.30,
                            pain_inh_ami=0.40, rifax_eff=0, prob_eff=0),
      "rifaximin"    = list(dose_alo=0, dose_pru=0, sert_inh_ami=0,
                            pain_inh_ami=0, rifax_eff=0.7, prob_eff=0.5),
      # none
      list(dose_alo=0, dose_pru=0, sert_inh_ami=0,
           pain_inh_ami=0, rifax_eff=0, prob_eff=0)
    )
  })

  ## Primary simulation (reactive to button)
  sim_data <- eventReactive(input$run_sim, {
    dp <- drug_params()
    run_ibs_sim(
      ibs_type     = as.integer(input$ibs_type),
      stress_level = input$stress,
      weeks        = input$weeks,
      dose_alo     = dp$dose_alo,
      dose_pru     = dp$dose_pru,
      sert_inh_ami = dp$sert_inh_ami,
      pain_inh_ami = dp$pain_inh_ami,
      rifax_eff    = dp$rifax_eff,
      prob_eff     = dp$prob_eff
    )
  }, ignoreNULL = FALSE)

  ## Untreated baseline (for comparison)
  baseline_data <- reactive({
    run_ibs_sim(
      ibs_type     = as.integer(input$ibs_type),
      stress_level = input$stress,
      weeks        = input$weeks
    )
  })

  ## ---- Value Boxes ----
  output$vb_pain <- renderValueBox({
    d <- sim_data()
    last <- tail(d, 1)
    valueBox(round(last$PAIN, 1), "Abdominal Pain (NRS)",
             icon = icon("exclamation-triangle"),
             color = ifelse(last$PAIN > 6, "red", ifelse(last$PAIN > 3, "yellow", "green")))
  })
  output$vb_sss <- renderValueBox({
    d <- sim_data()
    last <- tail(d, 1)
    valueBox(round(last$IBS_SSS, 0), "IBS-SSS Score",
             icon = icon("chart-bar"),
             color = ifelse(last$IBS_SSS > 300, "red", ifelse(last$IBS_SSS > 175, "yellow", "green")))
  })
  output$vb_stool <- renderValueBox({
    d <- sim_data()
    last <- tail(d, 1)
    valueBox(round(last$STOOL, 1), "Bristol Stool Form (1–7)",
             icon = icon("water"), color = "blue")
  })
  output$vb_qol <- renderValueBox({
    d <- sim_data()
    last <- tail(d, 1)
    qol_val <- round(max(0, min(100, 100 * (1 - last$IBS_SSS / 500))), 1)
    valueBox(paste0(qol_val, "%"), "Quality of Life (QoL)",
             icon = icon("smile"), color = ifelse(qol_val > 60, "green", "yellow"))
  })

  ## ---- Profile Tab ----
  output$plot_profile_sss <- renderPlotly({
    d_treated  <- sim_data() %>% mutate(Group = "Treated")
    d_baseline <- baseline_data() %>% mutate(Group = "Untreated")
    d <- bind_rows(d_treated, d_baseline)
    plot_ly(d, x = ~time_days, y = ~IBS_SSS, color = ~Group,
            type = "scatter", mode = "lines",
            colors = c("Treated" = "#3498DB", "Untreated" = "#E74C3C")) %>%
      add_segments(x = 0, xend = max(d$time_days), y = 175, yend = 175,
                   line = list(dash = "dash", color = "grey"), name = "Responder") %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "IBS-SSS (0–500)", range = c(0, 500)),
             title = "IBS-SSS: Natural History vs. Treatment")
  })

  ## ---- PK Tab ----
  output$tbl_pk_params <- renderDT({
    tibble(
      Drug      = c("Alosetron", "Prucalopride", "Amitriptyline", "Rifaximin"),
      Dose      = c("0.5 mg BID", "2 mg QD", "10–25 mg QN", "550 mg TID×14d"),
      t_half    = c("1.5 h", "24 h", "20 h", "6 h"),
      Vd        = c("65 L", "200 L", "1000 L", "30 L"),
      F_oral    = c("60%", "90%", "45%", "<0.4%"),
      IC50_EC50 = c("1.2 ng/mL", "2.5 ng/mL", "SERT Ki~25nM", "MIC90<0.1 mg/L")
    ) %>% datatable(options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  output$plot_pk_alo <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~Cp1, type = "scatter", mode = "lines",
            line = list(color = "#3498DB")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Alosetron Cp1 (ng/mL)"))
  })

  output$plot_pk_ro1 <- renderPlotly({
    d <- sim_data()
    ro <- 0.90 * d$Cp1 / (1.2 + d$Cp1 + 1e-9)
    plot_ly(d, x = ~time_days, y = ro, type = "scatter", mode = "lines",
            line = list(color = "#2980B9")) %>%
      add_segments(x=0, xend=max(d$time_days), y=0.8, yend=0.8,
                   line=list(dash="dash", color="red"), name="80% RO") %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "5-HT3 Receptor Occupancy", range = c(0, 1)))
  })

  output$plot_pk_ro2 <- renderPlotly({
    d <- sim_data()
    ro <- 0.85 * d$Cp3 / (2.5 + d$Cp3 + 1e-9)
    plot_ly(d, x = ~time_days, y = ro, type = "scatter", mode = "lines",
            line = list(color = "#27AE60")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "5-HT4 Receptor Occupancy (Prucalopride)", range = c(0, 1)))
  })

  ## ---- Biomarker Tab ----
  mk_bio_plot <- function(y_col, y_title, color) {
    d_t <- sim_data() %>% mutate(Group = "Treated")
    d_b <- baseline_data() %>% mutate(Group = "Untreated")
    d   <- bind_rows(d_t, d_b)
    plot_ly(d, x = ~time_days, y = as.formula(paste0("~", y_col)),
            color = ~Group, type = "scatter", mode = "lines",
            colors = c("Treated" = color, "Untreated" = "#E74C3C")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = y_title))
  }

  output$plot_bio_5ht    <- renderPlotly(mk_bio_plot("GUT_5HT", "5-HT (nmol/g)", "#F39C12"))
  output$plot_bio_mast   <- renderPlotly(mk_bio_plot("MAST_ACT", "Mast Activation (0–1)", "#9B59B6"))
  output$plot_bio_inflam <- renderPlotly(mk_bio_plot("INFLAM", "Inflammation (0–1)", "#E74C3C"))
  output$plot_bio_barrier<- renderPlotly(mk_bio_plot("BARRIER", "Barrier Integrity (0–1)", "#27AE60"))
  output$plot_bio_microb <- renderPlotly(mk_bio_plot("MICROB", "Dysbiosis (0–1)", "#E67E22"))
  output$plot_bio_scfa   <- renderPlotly(mk_bio_plot("SCFA", "SCFA (mmol/L)", "#16A085"))

  ## ---- Symptom Tab ----
  output$plot_sym_pain  <- renderPlotly(mk_bio_plot("PAIN",     "Pain NRS (0–10)", "#E74C3C"))
  output$plot_sym_bloat <- renderPlotly(mk_bio_plot("BLOAT",    "Bloating (0–10)", "#8E44AD"))
  output$plot_sym_stool <- renderPlotly(mk_bio_plot("STOOL",    "BSFS (1–7)",      "#3498DB"))
  output$plot_sym_sss   <- renderPlotly(mk_bio_plot("IBS_SSS",  "IBS-SSS (0–500)", "#E67E22"))
  output$plot_sym_vis   <- renderPlotly(mk_bio_plot("VIS_SENS", "Visceral Sens. (0–1)", "#C0392B"))
  output$plot_sym_motil <- renderPlotly(mk_bio_plot("MOTIL",    "Motility Index", "#27AE60"))

  ## ---- Scenario Comparison ----
  all_scenarios <- reactive({
    ibs_t <- as.integer(input$ibs_type)
    sl    <- input$stress
    wk    <- input$weeks

    scen_list <- list(
      list(label="① Untreated",           dose_alo=0,   dose_pru=0,    sert_inh_ami=0,    pain_inh_ami=0,    rifax_eff=0,   prob_eff=0, ibs_type=ibs_t),
      list(label="② Alosetron",            dose_alo=500, dose_pru=0,    sert_inh_ami=0,    pain_inh_ami=0,    rifax_eff=0,   prob_eff=0, ibs_type=1),
      list(label="③ Prucalopride",         dose_alo=0,   dose_pru=2000, sert_inh_ami=0,    pain_inh_ami=0,    rifax_eff=0,   prob_eff=0, ibs_type=2),
      list(label="④ Amitriptyline",        dose_alo=0,   dose_pru=0,    sert_inh_ami=0.30, pain_inh_ami=0.40, rifax_eff=0,   prob_eff=0, ibs_type=ibs_t),
      list(label="⑤ Rifaximin+Probiotics", dose_alo=0,   dose_pru=0,    sert_inh_ami=0,    pain_inh_ami=0,    rifax_eff=0.7, prob_eff=0.5, ibs_type=ibs_t)
    )

    purrr::map_dfr(scen_list, function(s) {
      run_ibs_sim(ibs_type=s$ibs_type, stress_level=sl, weeks=wk,
                  dose_alo=s$dose_alo, dose_pru=s$dose_pru,
                  sert_inh_ami=s$sert_inh_ami, pain_inh_ami=s$pain_inh_ami,
                  rifax_eff=s$rifax_eff, prob_eff=s$prob_eff) %>%
        mutate(label = s$label)
    })
  })

  output$plot_scen_sss <- renderPlotly({
    d <- all_scenarios()
    plot_ly(d, x = ~time_days, y = ~IBS_SSS, color = ~label,
            type = "scatter", mode = "lines") %>%
      add_segments(x=0, xend=max(d$time_days), y=175, yend=175,
                   line=list(dash="dash", color="grey"), name="Responder") %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "IBS-SSS (0–500)", range=c(0,500)),
             title  = "IBS Symptom Severity Score — All Treatment Scenarios")
  })

  output$plot_scen_pain <- renderPlotly({
    d <- all_scenarios()
    plot_ly(d, x = ~time_days, y = ~PAIN, color = ~label,
            type = "scatter", mode = "lines") %>%
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="Pain NRS", range=c(0,10)))
  })

  output$plot_scen_microb <- renderPlotly({
    d <- all_scenarios()
    plot_ly(d, x = ~time_days, y = ~MICROB, color = ~label,
            type = "scatter", mode = "lines") %>%
      layout(xaxis=list(title="Time (days)"), yaxis=list(title="Dysbiosis (0–1)", range=c(0,1)))
  })

  output$tbl_scen_summary <- renderDT({
    d <- all_scenarios()
    d %>%
      filter(abs(time_days - 84) == min(abs(time_days - 84))) %>%
      group_by(label) %>%
      slice(1) %>%
      summarise(
        `IBS-SSS`    = round(mean(IBS_SSS), 1),
        `Pain NRS`   = round(mean(PAIN), 2),
        `Bloating`   = round(mean(BLOAT), 2),
        `BSFS`       = round(mean(STOOL), 2),
        `Dysbiosis`  = round(mean(MICROB), 3),
        `Barrier`    = round(mean(BARRIER), 3),
        `5-HT`       = round(mean(GUT_5HT), 3),
        `QoL (%)`    = round(pmax(0, 100*(1-mean(IBS_SSS)/500)), 1),
        .groups = "drop"
      ) %>%
      datatable(rownames = FALSE, options = list(dom = "t", pageLength = 10))
  })

  ## ---- Biomarker Tracker Tab ----
  output$plot_tr_qol <- renderPlotly({
    d <- sim_data() %>%
      mutate(QoL = pmax(0, pmin(100, 100 * (1 - IBS_SSS / 500))))
    plot_ly(d, x = ~time_days, y = ~QoL, type = "scatter", mode = "lines",
            line = list(color = "#2ECC71")) %>%
      layout(xaxis=list(title="Time (days)"),
             yaxis=list(title="QoL Score (0–100)", range=c(0,100)))
  })

  output$plot_tr_drg <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~DRG_ACT, type = "scatter", mode = "lines",
            line = list(color = "#E74C3C")) %>%
      layout(xaxis=list(title="Time (days)"),
             yaxis=list(title="DRG Activation (0–1)", range=c(0,1)))
  })

  output$plot_tr_cga <- renderPlotly({
    d <- sim_data() %>% mutate(CgA = 50 + 80 * GUT_5HT)
    plot_ly(d, x = ~time_days, y = ~CgA, type = "scatter", mode = "lines",
            line = list(color = "#F39C12")) %>%
      layout(xaxis=list(title="Time (days)"),
             yaxis=list(title="Serum Chromogranin A (proxy, ng/mL)"))
  })

  output$plot_tr_rectal <- renderPlotly({
    d <- sim_data() %>% mutate(RecThr = 20 * (1 - VIS_SENS))
    plot_ly(d, x = ~time_days, y = ~RecThr, type = "scatter", mode = "lines",
            line = list(color = "#9B59B6")) %>%
      add_segments(x=0, xend=max(d$time_days), y=12, yend=12,
                   line=list(dash="dash", color="red"), name="Hypersensitivity") %>%
      layout(xaxis=list(title="Time (days)"),
             yaxis=list(title="Rectal Threshold (mmHg)"),
             title = "Higher threshold = less hypersensitivity")
  })

  output$tbl_biomarker_summary <- renderDT({
    d <- sim_data()
    key_days <- c(28, 56, 84, 168)
    purrr::map_dfr(key_days, function(kd) {
      d %>%
        filter(abs(time_days - kd) == min(abs(time_days - kd))) %>%
        slice(1) %>%
        mutate(
          Week       = kd / 7,
          `5-HT`     = round(GUT_5HT, 3),
          `Mast`     = round(MAST_ACT, 3),
          `Inflam`   = round(INFLAM, 3),
          `Barrier`  = round(BARRIER, 3),
          `Dysbiosis`= round(MICROB, 3),
          `SCFA`     = round(SCFA, 2),
          `Pain`     = round(PAIN, 2),
          `IBS-SSS`  = round(IBS_SSS, 1),
          `QoL %`    = round(pmax(0, 100*(1-IBS_SSS/500)), 1)
        ) %>%
        select(Week, `5-HT`, Mast, Inflam, Barrier, Dysbiosis, SCFA, Pain, `IBS-SSS`, `QoL %`)
    }) %>%
      datatable(rownames=FALSE, options=list(dom="t", pageLength=10))
  })

  ## ---- Download ----
  output$dl_data <- downloadHandler(
    filename = function() paste0("ibs_simulation_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(sim_data(), file, row.names = FALSE)
  )
}

## ============================================================
## 5. Launch
## ============================================================

shinyApp(ui, server)
