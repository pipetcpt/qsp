## =============================================================================
## Polymyalgia Rheumatica (PMR) — Interactive QSP Shiny Dashboard
## =============================================================================
## 6 Tabs:
##   Tab 1: Patient Profile & Disease Setup
##   Tab 2: Drug PK (Prednisolone + Tocilizumab)
##   Tab 3: Inflammatory Markers (IL-6, CRP, ESR)
##   Tab 4: Disease Activity (PMR-AS, Pain, Stiffness)
##   Tab 5: Scenario Comparison (7 treatment arms)
##   Tab 6: Biomarker Explorer (Bone, HPA, GC-sparing)
## =============================================================================

library(shiny)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(mrgsolve)

## ─────────────────────────────────────────────────────────────────────────────
## Inline mrgsolve model (same as in the .R model file)
## ─────────────────────────────────────────────────────────────────────────────

pmr_model_code <- '
$PARAM
KA_PRED=2.0 F_PRED=0.82 CL_PRED=14.0 V1_PRED=30.0 V2_PRED=50.0
Q_PRED=8.0 FU_PRED=0.28
KA_TCZ=0.012 F_TCZ=0.80 CL_TCZ=0.29 V1_TCZ=3.6 V2_TCZ=2.5
Q_TCZ=0.35 KINT=0.003
KSYN_CORT=0.14 KEL_CORT=0.46 CORT0=12.0 IC50_HPA=2.5 IMAX_HPA=0.95
KSYN_IL6=6.0 KEL_IL6=0.35 IL6_BASE=15.0 EC50_IL6=200.0 EMAX_IL6=0.90
KSYN_sR=0.08 KEL_sR=0.005 SIL6R_BASE=40.0
KSYN_CRP=0.018 KEL_CRP=0.025 CRP_BASE=35.0 STIM_CRP=0.8
KSYN_ESR=0.05 KEL_ESR=0.015 ESR_BASE=55.0
BMD_BASE=1.0 KBMD_LOSS=0.0003 KBMD_REC=0.00005
PMRAS_MAX=55.0 PMRAS_BASE=2.0 EC50_PMRAS=180.0 EMAX_PMRAS=0.92
EC50_TCZ_PMRAS=50.0 EMAX_TCZ_PMRAS=0.70
K_RELAPSE=0.0001 PRED_PROTECT=10.0

$CMT DEPOT_PRED CENT_PRED PERI_PRED DEPOT_TCZ CENT_TCZ PERI_TCZ
     CORT IL6 SIL6R CRP ESR BMD PMRAS FLARE

$MAIN
if(NEWIND<=1){
  _init_CORT=CORT0; _init_IL6=IL6_BASE; _init_SIL6R=SIL6R_BASE;
  _init_CRP=CRP_BASE; _init_ESR=ESR_BASE; _init_BMD=BMD_BASE;
  _init_PMRAS=PMRAS_MAX; _init_FLARE=0.01;
}

$ODE
double k12_P=Q_PRED/V1_PRED, k21_P=Q_PRED/V2_PRED, kel_P=CL_PRED/V1_PRED;
double Cp_FREE_h=CENT_PRED/V1_PRED*FU_PRED*1000.0;
dxdt_DEPOT_PRED=-KA_PRED*DEPOT_PRED;
dxdt_CENT_PRED=KA_PRED*DEPOT_PRED*F_PRED-(kel_P+k12_P)*CENT_PRED+k21_P*PERI_PRED;
dxdt_PERI_PRED=k12_P*CENT_PRED-k21_P*PERI_PRED;

double k12_T=Q_TCZ/V1_TCZ, k21_T=Q_TCZ/V2_TCZ, kel_T=CL_TCZ/V1_TCZ;
double kel_TMDD=KINT*SIL6R/(SIL6R_BASE+SIL6R);
double Cp_TCZnM=CENT_TCZ/V1_TCZ/148.0*1e6;
dxdt_DEPOT_TCZ=-KA_TCZ*DEPOT_TCZ;
dxdt_CENT_TCZ=KA_TCZ*DEPOT_TCZ*F_TCZ-(kel_T+kel_TMDD+k12_T)*CENT_TCZ+k21_T*PERI_TCZ;
dxdt_PERI_TCZ=k12_T*CENT_TCZ-k21_T*PERI_TCZ;

double inh_HPA=(IMAX_HPA*Cp_FREE_h)/(IC50_HPA+Cp_FREE_h);
dxdt_CORT=KSYN_CORT*(1-inh_HPA)-KEL_CORT*CORT;

double disease_factor=PMRAS/PMRAS_MAX;
double ksyn_IL6_eff=KSYN_IL6*(1+3.0*disease_factor);
double inh_IL6_GC=(EMAX_IL6*Cp_FREE_h)/(EC50_IL6+Cp_FREE_h);
double TCZ_occ_sR=(Cp_TCZnM>0)?Cp_TCZnM/(0.5+Cp_TCZnM):0.0;
double inh_IL6_TCZ=0.85*TCZ_occ_sR;
double eff_inh_IL6=1-(1-inh_IL6_GC)*(1-inh_IL6_TCZ);
dxdt_IL6=ksyn_IL6_eff*(1-eff_inh_IL6)-KEL_IL6*IL6;

double sIL6R_stim=1.0+0.4*TCZ_occ_sR;
dxdt_SIL6R=KSYN_sR*sIL6R_stim-KEL_sR*SIL6R;

double IL6_ratio=IL6/IL6_BASE;
double ksyn_CRP_eff=KSYN_CRP*(1+STIM_CRP*(IL6_ratio-1));
dxdt_CRP=ksyn_CRP_eff-KEL_CRP*CRP;
double ESR_target=ESR_BASE*(IL6/IL6_BASE)*(CRP/CRP_BASE+1)/2.0;
dxdt_ESR=KEL_ESR*(ESR_target-ESR);

double daily_pred_mg=CENT_PRED/V1_PRED*CL_PRED*24.0;
dxdt_BMD=-KBMD_LOSS*daily_pred_mg+KBMD_REC*(BMD_BASE-BMD);

double eff_GC_PMRAS=(EMAX_PMRAS*Cp_FREE_h)/(EC50_PMRAS+Cp_FREE_h);
double eff_TCZ_PMRAS=(EMAX_TCZ_PMRAS*Cp_TCZnM)/(EC50_TCZ_PMRAS+Cp_TCZnM);
double combined_eff=1-(1-eff_GC_PMRAS)*(1-eff_TCZ_PMRAS);
double PMRAS_eq=PMRAS_MAX*(1-combined_eff);
if(PMRAS_eq<PMRAS_BASE) PMRAS_eq=PMRAS_BASE;
dxdt_PMRAS=0.02*(PMRAS_eq-PMRAS);

double pred_dose_c=CENT_PRED/V1_PRED*CL_PRED*24.0;
double relapse_driver=(pred_dose_c<PRED_PROTECT)?K_RELAPSE*(PRED_PROTECT-pred_dose_c)/PRED_PROTECT:0;
dxdt_FLARE=relapse_driver-0.002*FLARE;

$TABLE
double Cp_PRED_out=CENT_PRED/V1_PRED;
double Cp_FREE_out=CENT_PRED/V1_PRED*FU_PRED*1000.0;
double Cp_TCZ_out=CENT_TCZ/V1_TCZ;
double Cp_TCZnM_out=CENT_TCZ/V1_TCZ/148.0*1e6;
double CRP_norm=CRP/CRP_BASE;

$CAPTURE Cp_PRED_out Cp_FREE_out Cp_TCZ_out Cp_TCZnM_out
         CORT IL6 SIL6R CRP ESR BMD PMRAS FLARE CRP_norm
'

## Compile model (cached with options)
options(mrgsolve.soloc = tempdir())
mod <- suppressMessages(mcode("PMR_Shiny", pmr_model_code))

## ─────────────────────────────────────────────────────────────────────────────
## Simulation helper function
## ─────────────────────────────────────────────────────────────────────────────

run_sim <- function(
    pred_dose_mg = 15,
    taper_start_wk = 4,
    taper_rate = 2.5,
    min_pred_mg = 0,
    tcz_dose_mg = 0,
    tcz_interval_h = 336,
    sim_days = 365,
    pmras_max = 50,
    il6_base = 15,
    crp_base = 35,
    esr_base = 55,
    ec50_il6 = 200,
    ec50_pmras = 180,
    patient_id = 1
) {
  # Build dosing events
  ev_list <- list()

  # Prednisolone BID tapering
  admin_times <- seq(0, (sim_days - 1) * 24, by = 12)
  taper_start_h <- taper_start_wk * 7 * 24
  doses_per_admin <- sapply(admin_times, function(h) {
    mo_after <- max(0, (h - taper_start_h) / (30.4375 * 24))
    current_total <- max(min_pred_mg, pred_dose_mg - taper_rate * mo_after)
    current_total / 2
  })
  ev_pred <- ev(cmt = 1, amt = doses_per_admin, time = admin_times)
  ev_list[["pred"]] <- ev_pred

  # Tocilizumab SC
  if (tcz_dose_mg > 0) {
    tcz_times <- seq(0, (sim_days - 1) * 24, by = tcz_interval_h)
    ev_list[["tcz"]] <- ev(cmt = 4, amt = tcz_dose_mg, time = tcz_times)
  }

  events <- Reduce(c, ev_list)

  # Update patient-specific parameters
  iparams <- list(
    PMRAS_MAX = pmras_max,
    IL6_BASE  = il6_base,
    CRP_BASE  = crp_base,
    ESR_BASE  = esr_base,
    EC50_IL6  = ec50_il6,
    EC50_PMRAS = ec50_pmras
  )
  mod_run <- param(mod, iparams)

  out <- mrgsim(mod_run, ev = events,
                end = sim_days * 24, delta = 12) %>%
    as.data.frame() %>%
    mutate(Time_days = time / 24, Time_weeks = time / 168, Patient = patient_id)

  out
}

## ─────────────────────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = "PMR QSP Dashboard",
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,

    sidebarMenu(
      id = "tabs",
      menuItem("Patient Profile", tabName = "tab_patient",
               icon = icon("user-md")),
      menuItem("Drug PK", tabName = "tab_pk",
               icon = icon("pills")),
      menuItem("Inflammatory Markers", tabName = "tab_inflam",
               icon = icon("fire")),
      menuItem("Disease Activity", tabName = "tab_disease",
               icon = icon("chart-line")),
      menuItem("Scenario Comparison", tabName = "tab_scenarios",
               icon = icon("layer-group")),
      menuItem("Biomarker Explorer", tabName = "tab_biomarker",
               icon = icon("flask"))
    ),

    hr(),
    h5("Patient Parameters", style = "color:#ECF0F1; padding-left:15px;"),

    sliderInput("pmras_max", "Initial PMR-AS (severity):",
                min = 15, max = 70, value = 45, step = 1),
    sliderInput("il6_base", "Baseline IL-6 (pg/mL):",
                min = 5, max = 120, value = 20, step = 1),
    sliderInput("crp_base", "Baseline CRP (mg/L):",
                min = 5, max = 150, value = 40, step = 1),
    sliderInput("esr_base", "Baseline ESR (mm/hr):",
                min = 20, max = 120, value = 60, step = 1),

    hr(),
    h5("Treatment: Prednisolone", style = "color:#ECF0F1; padding-left:15px;"),
    sliderInput("pred_dose", "Initial Dose (mg/day):",
                min = 0, max = 30, value = 15, step = 2.5),
    sliderInput("taper_start", "Start Taper (weeks):",
                min = 1, max = 12, value = 4, step = 1),
    sliderInput("taper_rate", "Taper Rate (mg/month):",
                min = 0, max = 5, value = 2.5, step = 0.5),
    sliderInput("min_pred", "Minimum Pred Dose (mg/day):",
                min = 0, max = 10, value = 0, step = 1),

    hr(),
    h5("Treatment: Tocilizumab", style = "color:#ECF0F1; padding-left:15px;"),
    checkboxInput("use_tcz", "Add Tocilizumab", value = FALSE),
    conditionalPanel(
      "input.use_tcz == true",
      sliderInput("tcz_dose", "TCZ Dose (mg/injection):",
                  min = 80, max = 162, value = 162, step = 40),
      radioButtons("tcz_freq", "TCZ Frequency:",
                   choices = c("Weekly (Q1W)" = 168, "Biweekly (Q2W)" = 336),
                   selected = 336)
    ),

    hr(),
    sliderInput("sim_days", "Simulation Duration (days):",
                min = 90, max = 730, value = 365, step = 30),
    actionButton("run_btn", "Run Simulation", icon = icon("play"),
                 class = "btn-success", width = "90%")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #F5F7FA; }
      .box-title { font-weight: bold; }
      .value-box .inner h3 { font-size: 26px; }
      .tab-content { padding-top: 10px; }
      .nav-tabs-custom .nav-tabs li.active a { border-top-color: #1565C0; }
    "))),

    tabItems(

      ## ── TAB 1: PATIENT PROFILE ─────────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          valueBoxOutput("vbox_pmras", width = 3),
          valueBoxOutput("vbox_crp", width = 3),
          valueBoxOutput("vbox_esr", width = 3),
          valueBoxOutput("vbox_il6", width = 3)
        ),
        fluidRow(
          box(title = "PMR Disease Overview", status = "primary", solidHeader = TRUE,
              width = 6,
              tags$div(
                tags$h4("Polymyalgia Rheumatica (PMR)"),
                tags$p("PMR is a common inflammatory condition affecting adults aged ≥50, characterized by:"),
                tags$ul(
                  tags$li("Bilateral aching and stiffness of shoulder and hip girdles"),
                  tags$li("Morning stiffness lasting >45 minutes"),
                  tags$li("Elevated acute-phase reactants (CRP, ESR, IL-6)"),
                  tags$li("Rapid response to glucocorticoids (hallmark)"),
                  tags$li("Overlap with Giant Cell Arteritis in 15–20%")
                ),
                tags$p(strong("Epidemiology:"), " Incidence 50–100/100,000/year in ≥50y; F:M ≈ 2–3:1; Nordic ancestry higher risk"),
                tags$p(strong("Treatment:"), " Prednisolone 12.5–25mg/d with slow taper; TCZ approved for GCA, RCT evidence in PMR (SEMAPHORE, NCT03600818)")
              )
          ),
          box(title = "Diagnosis Criteria (ACR/EULAR 2012)", status = "info", solidHeader = TRUE,
              width = 6,
              tableOutput("diagnosis_table")
          )
        ),
        fluidRow(
          box(title = "Baseline Parameter Summary (Current Patient)", status = "warning",
              solidHeader = TRUE, width = 12,
              DT::dataTableOutput("baseline_table"))
        )
      ),

      ## ── TAB 2: DRUG PK ─────────────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Prednisolone Plasma Concentration (Total)", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_pred_pk", height = 350)),
          box(title = "Free Prednisolone (Active Fraction)", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_pred_free", height = 350))
        ),
        fluidRow(
          box(title = "Tocilizumab Plasma Concentration (nM)", status = "info",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_tcz_pk", height = 350)),
          box(title = "PK Summary Statistics", status = "info",
              solidHeader = TRUE, width = 6,
              DT::dataTableOutput("pk_summary_table"))
        )
      ),

      ## ── TAB 3: INFLAMMATORY MARKERS ────────────────────────────────────────
      tabItem(tabName = "tab_inflam",
        fluidRow(
          box(title = "Plasma IL-6 Dynamics", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_il6", height = 350)),
          box(title = "Soluble IL-6R (sIL-6Rα)", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_sil6r", height = 350))
        ),
        fluidRow(
          box(title = "CRP Response", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_crp", height = 350)),
          box(title = "ESR Response", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_esr", height = 350))
        ),
        fluidRow(
          box(title = "Inflammatory Marker Summary by Timepoint",
              status = "primary", solidHeader = TRUE, width = 12,
              DT::dataTableOutput("inflam_summary"))
        )
      ),

      ## ── TAB 4: DISEASE ACTIVITY ─────────────────────────────────────────────
      tabItem(tabName = "tab_disease",
        fluidRow(
          box(title = "PMR Activity Score (PMR-AS) Over Time", status = "success",
              solidHeader = TRUE, width = 8,
              plotlyOutput("plot_pmras", height = 400)),
          box(title = "PMR-AS Interpretation", status = "success",
              solidHeader = TRUE, width = 4,
              tags$div(
                tags$h5("PMR-AS Score Thresholds:"),
                tags$table(class = "table table-condensed",
                  tags$thead(tags$tr(
                    tags$th("Score"), tags$th("Interpretation")
                  )),
                  tags$tbody(
                    tags$tr(tags$td("≤ 7"),    tags$td(style="color:green", "Remission")),
                    tags$tr(tags$td("7–17.5"), tags$td(style="color:orange", "Low Activity")),
                    tags$tr(tags$td(">17.5"),  tags$td(style="color:red",   "Active PMR"))
                  )
                ),
                tags$hr(),
                tags$h5("PMR-AS Formula:"),
                tags$p("PMR-AS = 2.45 × VAS_pain + 0.02 × ESR + 0.70 × PGA + 0.35 × EL_stiffness + 0.58 × HAQ-DI"),
                tags$h5("Target:"),
                tags$p("PMR-AS ≤ 7 = clinical remission; aim by week 4–8 with standard GC therapy")
              )
          )
        ),
        fluidRow(
          box(title = "HPA Axis: Endogenous Cortisol Suppression", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_cortisol", height = 300)),
          box(title = "Flare/Relapse Risk Score", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("plot_flare", height = 300))
        )
      ),

      ## ── TAB 5: SCENARIO COMPARISON ─────────────────────────────────────────
      tabItem(tabName = "tab_scenarios",
        fluidRow(
          box(title = "Select Scenarios to Compare", status = "primary",
              solidHeader = TRUE, width = 12,
              checkboxGroupInput("selected_scenarios",
                label = NULL,
                choices = list(
                  "S1: No Treatment (Natural History)" = "S1",
                  "S2: Pred 15mg → Taper 2.5mg/mo (ACR Standard)" = "S2",
                  "S3: Pred 22.5mg → Rapid Taper 4mg/mo" = "S3",
                  "S4: Pred 15mg → Slow Taper 1mg/mo" = "S4",
                  "S5: TCZ 162mg QW + Pred 12.5mg → Taper" = "S5",
                  "S6: TCZ 162mg Q2W + Pred 12.5mg → Taper" = "S6",
                  "S7: TCZ QW Only (Steroid-Free)" = "S7"
                ),
                selected = c("S2", "S5", "S6", "S7"),
                inline = FALSE
              )
          )
        ),
        fluidRow(
          box(title = "CRP Comparison by Scenario", status = "primary",
              solidHeader = TRUE, width = 6,
              plotlyOutput("comp_crp", height = 350)),
          box(title = "PMR-AS Comparison by Scenario", status = "success",
              solidHeader = TRUE, width = 6,
              plotlyOutput("comp_pmras", height = 350))
        ),
        fluidRow(
          box(title = "IL-6 Comparison by Scenario", status = "danger",
              solidHeader = TRUE, width = 6,
              plotlyOutput("comp_il6", height = 350)),
          box(title = "BMD Comparison (GC Bone Safety)", status = "warning",
              solidHeader = TRUE, width = 6,
              plotlyOutput("comp_bmd", height = 350))
        ),
        fluidRow(
          box(title = "Scenario Summary at Key Timepoints",
              status = "info", solidHeader = TRUE, width = 12,
              DT::dataTableOutput("scenario_summary_table"))
        )
      ),

      ## ── TAB 6: BIOMARKER EXPLORER ───────────────────────────────────────────
      tabItem(tabName = "tab_biomarker",
        fluidRow(
          box(title = "Bone Mineral Density (GC-Induced Osteoporosis Risk)",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("bio_bmd", height = 320)),
          box(title = "Endogenous Cortisol vs Prednisolone Free Cp",
              status = "info", solidHeader = TRUE, width = 6,
              plotlyOutput("bio_cortisol", height = 320))
        ),
        fluidRow(
          box(title = "GC Sparing Effect: Cumulative Pred Dose",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("bio_gc_sparing", height = 320)),
          box(title = "Biomarker Correlation: IL-6 vs CRP vs PMR-AS",
              status = "success", solidHeader = TRUE, width = 6,
              plotlyOutput("bio_correlation", height = 320))
        ),
        fluidRow(
          box(title = "Clinical Reference Values & PMR Risk Markers",
              status = "danger", solidHeader = TRUE, width = 12,
              DT::dataTableOutput("biomarker_ref_table"))
        )
      )

    ) # end tabItems
  )   # end dashboardBody
)     # end dashboardPage

## ─────────────────────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  ## Reactive: run main simulation
  sim_data <- eventReactive(input$run_btn, {
    withProgress(message = "Running PMR Simulation...", value = 0.5, {
      tcz_mg <- if (input$use_tcz) input$tcz_dose else 0
      tcz_iv  <- if (input$use_tcz) as.numeric(input$tcz_freq) else 336

      run_sim(
        pred_dose_mg   = input$pred_dose,
        taper_start_wk = input$taper_start,
        taper_rate     = input$taper_rate,
        min_pred_mg    = input$min_pred,
        tcz_dose_mg    = tcz_mg,
        tcz_interval_h = tcz_iv,
        sim_days       = input$sim_days,
        pmras_max      = input$pmras_max,
        il6_base       = input$il6_base,
        crp_base       = input$crp_base,
        esr_base       = input$esr_base
      )
    })
  }, ignoreNULL = FALSE)

  ## Reactive: run all 7 scenarios (for comparison tab)
  all_scenarios <- eventReactive(input$run_btn, {
    withProgress(message = "Running 7 Scenarios...", {
      sc_defs <- list(
        S1 = list(pred=0,    taper=0,   tcz=0,   tcz_h=336, label="S1: No Treatment",
                  col="#999999"),
        S2 = list(pred=15,   taper=2.5, tcz=0,   tcz_h=336, label="S2: Pred 15mg Std",
                  col="#2196F3"),
        S3 = list(pred=22.5, taper=4.0, tcz=0,   tcz_h=336, label="S3: Pred 22.5mg Rapid",
                  col="#FF9800"),
        S4 = list(pred=15,   taper=1.0, tcz=0,   tcz_h=336, label="S4: Pred 15mg Slow",
                  col="#9C27B0"),
        S5 = list(pred=12.5, taper=2.5, tcz=162, tcz_h=168, label="S5: TCZ QW + Pred",
                  col="#4CAF50"),
        S6 = list(pred=12.5, taper=2.5, tcz=162, tcz_h=336, label="S6: TCZ Q2W + Pred",
                  col="#009688"),
        S7 = list(pred=0,    taper=0,   tcz=162, tcz_h=168, label="S7: TCZ Only (GC-Free)",
                  col="#F44336")
      )

      results <- lapply(names(sc_defs), function(sname) {
        sc <- sc_defs[[sname]]
        setProgress(message = paste("Running", sc$label))
        out <- run_sim(
          pred_dose_mg   = sc$pred,
          taper_start_wk = input$taper_start,
          taper_rate     = sc$taper,
          tcz_dose_mg    = sc$tcz,
          tcz_interval_h = sc$tcz_h,
          sim_days       = input$sim_days,
          pmras_max      = input$pmras_max,
          il6_base       = input$il6_base,
          crp_base       = input$crp_base,
          esr_base       = input$esr_base
        )
        out$Scenario   <- sc$label
        out$ScenarioID <- sname
        out$Color      <- sc$col
        out
      })
      bind_rows(results)
    })
  }, ignoreNULL = FALSE)

  ## ── Value boxes ─────────────────────────────────────────────────────────────
  output$vbox_pmras <- renderValueBox({
    valueBox(
      value = input$pmras_max,
      subtitle = "Initial PMR-AS",
      icon = icon("chart-bar"),
      color = if (input$pmras_max > 35) "red" else if (input$pmras_max > 17) "orange" else "green"
    )
  })
  output$vbox_crp <- renderValueBox({
    valueBox(
      value = paste(input$crp_base, "mg/L"),
      subtitle = "Baseline CRP",
      icon = icon("fire"),
      color = if (input$crp_base > 50) "red" else if (input$crp_base > 20) "orange" else "yellow"
    )
  })
  output$vbox_esr <- renderValueBox({
    valueBox(
      value = paste(input$esr_base, "mm/hr"),
      subtitle = "Baseline ESR",
      icon = icon("tint"),
      color = if (input$esr_base > 70) "red" else if (input$esr_base > 40) "orange" else "yellow"
    )
  })
  output$vbox_il6 <- renderValueBox({
    valueBox(
      value = paste(input$il6_base, "pg/mL"),
      subtitle = "Baseline IL-6",
      icon = icon("bacterium"),
      color = if (input$il6_base > 30) "red" else "orange"
    )
  })

  ## ── Diagnosis table ──────────────────────────────────────────────────────────
  output$diagnosis_table <- renderTable({
    data.frame(
      Criterion = c("Age ≥ 50 years", "Bilateral shoulder aching",
                    "Elevated CRP and/or ESR",
                    "Morning stiffness ≥ 45 min",
                    "Hip girdle pain or limited ROM",
                    "Negative RF and anti-CCP",
                    "Absence of other diagnosis"),
      Points = c("Required", "+2", "+2", "+1", "+1", "+2", "Exclusion"),
      Note = c("Mandatory", "Key feature", "CRP>5mg/L or ESR>50",
               "Key feature", "Present in ~70%",
               "Mandatory", "")
    )
  }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "s")

  ## ── Baseline table ───────────────────────────────────────────────────────────
  output$baseline_table <- DT::renderDataTable({
    data.frame(
      Parameter = c("PMR Activity Score (PMR-AS)", "IL-6", "CRP", "ESR",
                    "Prednisolone Dose", "Tocilizumab", "Simulation Duration"),
      Value = c(input$pmras_max, paste(input$il6_base, "pg/mL"),
                paste(input$crp_base, "mg/L"), paste(input$esr_base, "mm/hr"),
                paste(input$pred_dose, "mg/day"),
                if (input$use_tcz) paste(input$tcz_dose, "mg",
                                         if (input$tcz_freq == 168) "QW" else "Q2W") else "None",
                paste(input$sim_days, "days")),
      Category = c("Disease Activity", "Biomarker", "Biomarker", "Biomarker",
                   "Treatment", "Treatment", "Simulation")
    )
  }, options = list(pageLength = 10, dom = 't'), rownames = FALSE)

  ## ── PK plots ─────────────────────────────────────────────────────────────────
  output$plot_pred_pk <- renderPlotly({
    dat <- sim_data()
    p <- ggplot(dat, aes(x = Time_days, y = Cp_PRED_out)) +
      geom_line(color = "#2196F3", linewidth = 0.9) +
      labs(x = "Time (days)", y = "Prednisolone (μg/mL)", title = "Total Prednisolone") +
      theme_bw()
    ggplotly(p, tooltip = c("x", "y"))
  })

  output$plot_pred_free <- renderPlotly({
    dat <- sim_data()
    p <- ggplot(dat, aes(x = Time_days, y = Cp_FREE_out)) +
      geom_line(color = "#1565C0", linewidth = 0.9) +
      geom_hline(yintercept = 180, linetype = "dashed", color = "red", alpha = 0.6) +
      annotate("text", x = max(dat$Time_days) * 0.7, y = 195,
               label = "EC50 for PMR-AS (180 ng/mL)", size = 3, color = "red") +
      labs(x = "Time (days)", y = "Free Prednisolone (ng/mL)", title = "Free (Active) Prednisolone") +
      theme_bw()
    ggplotly(p, tooltip = c("x", "y"))
  })

  output$plot_tcz_pk <- renderPlotly({
    dat <- sim_data()
    p <- ggplot(dat, aes(x = Time_days, y = Cp_TCZnM_out)) +
      geom_line(color = "#009688", linewidth = 0.9) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "orange", alpha = 0.7) +
      annotate("text", x = max(dat$Time_days) * 0.7, y = 55,
               label = "EC50 PMR-AS (50 nM)", size = 3, color = "orange") +
      labs(x = "Time (days)", y = "Tocilizumab (nM)", title = "Tocilizumab Cp") +
      theme_bw()
    ggplotly(p, tooltip = c("x", "y"))
  })

  output$pk_summary_table <- DT::renderDataTable({
    dat <- sim_data()
    timepoints <- c(1, 7, 14, 28, 56, 84)
    pk_df <- dat %>%
      filter(abs(Time_days - timepoints[1]) < 0.6 |
             abs(Time_days - timepoints[2]) < 0.6 |
             abs(Time_days - timepoints[3]) < 0.6 |
             abs(Time_days - timepoints[4]) < 0.6 |
             abs(Time_days - timepoints[5]) < 0.6 |
             abs(Time_days - timepoints[6]) < 0.6) %>%
      mutate(Day = round(Time_days)) %>%
      filter(Day %in% timepoints) %>%
      group_by(Day) %>%
      summarise(
        Pred_ug_mL  = round(max(Cp_PRED_out, na.rm=TRUE), 3),
        PredFree_ng = round(max(Cp_FREE_out, na.rm=TRUE), 1),
        TCZ_nM      = round(max(Cp_TCZnM_out, na.rm=TRUE), 2),
        .groups = "drop"
      )
    pk_df
  }, options = list(dom = 't'), rownames = FALSE)

  ## ── Inflammatory marker plots ────────────────────────────────────────────────
  output$plot_il6 <- renderPlotly({
    dat <- sim_data()
    p <- ggplot(dat, aes(x = Time_days, y = IL6)) +
      geom_line(color = "#F44336", linewidth = 0.9) +
      geom_hline(yintercept = 3.4, linetype = "dashed", color = "gray40") +
      annotate("text", x = max(dat$Time_days)*0.7, y = 5, label = "Normal (<3.4 pg/mL)", size=3) +
      labs(x="Time (days)", y="IL-6 (pg/mL)", title="Plasma IL-6") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_sil6r <- renderPlotly({
    dat <- sim_data()
    p <- ggplot(dat, aes(x = Time_days, y = SIL6R)) +
      geom_line(color = "#E91E63", linewidth = 0.9) +
      geom_hline(yintercept = 40, linetype = "dashed", color = "gray40") +
      annotate("text", x = max(dat$Time_days)*0.7, y = 42, label = "Baseline (40 ng/mL)", size=3) +
      labs(x="Time (days)", y="sIL-6Rα (ng/mL)", title="Soluble IL-6 Receptor") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_crp <- renderPlotly({
    dat <- sim_data()
    p <- ggplot(dat, aes(x = Time_days, y = CRP)) +
      geom_line(color = "#FF9800", linewidth = 0.9) +
      geom_hline(yintercept = 5, linetype = "dashed", color = "darkgreen") +
      annotate("text", x = max(dat$Time_days)*0.7, y = 7, label = "Normal (<5 mg/L)", size=3) +
      labs(x="Time (days)", y="CRP (mg/L)", title="C-Reactive Protein") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_esr <- renderPlotly({
    dat <- sim_data()
    p <- ggplot(dat, aes(x = Time_days, y = ESR)) +
      geom_line(color = "#FF6F00", linewidth = 0.9) +
      geom_hline(yintercept = 20, linetype = "dashed", color = "darkgreen") +
      annotate("text", x = max(dat$Time_days)*0.7, y = 22, label = "Normal (<20 mm/hr)", size=3) +
      labs(x="Time (days)", y="ESR (mm/hr)", title="Erythrocyte Sedimentation Rate") +
      theme_bw()
    ggplotly(p)
  })

  output$inflam_summary <- DT::renderDataTable({
    dat <- sim_data()
    days_check <- c(0, 14, 28, 84, 182, 365)
    dat %>%
      filter(sapply(Time_days, function(d) any(abs(d - days_check) < 1))) %>%
      mutate(Day = round(Time_days)) %>%
      filter(Day %in% days_check) %>%
      group_by(Day) %>%
      slice(1) %>%
      ungroup() %>%
      select(Day, IL6, SIL6R, CRP, ESR) %>%
      mutate(across(where(is.numeric), ~round(., 1)))
  }, options = list(dom='t'), rownames=FALSE)

  ## ── Disease activity plots ───────────────────────────────────────────────────
  output$plot_pmras <- renderPlotly({
    dat <- sim_data()
    p <- ggplot(dat, aes(x = Time_days, y = PMRAS)) +
      geom_line(color = "#4CAF50", linewidth = 1.1) +
      geom_hline(yintercept = 7,    linetype = "dashed", color = "darkgreen", linewidth = 0.7) +
      geom_hline(yintercept = 17.5, linetype = "dashed", color = "orange",    linewidth = 0.7) +
      geom_ribbon(aes(ymin = 0, ymax = 7), fill = "green",  alpha = 0.05) +
      geom_ribbon(aes(ymin = 7, ymax = 17.5), fill = "yellow", alpha = 0.05) +
      geom_ribbon(aes(ymin = 17.5, ymax = max(PMRAS, na.rm=T)*1.1), fill = "red", alpha = 0.05) +
      labs(x="Time (days)", y="PMR-AS Score", title="PMR Activity Score (PMR-AS)") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_cortisol <- renderPlotly({
    dat <- sim_data()
    p <- ggplot(dat, aes(x = Time_days, y = CORT)) +
      geom_line(color = "#FFC107", linewidth = 0.9) +
      geom_hline(yintercept = 12, linetype = "dashed", color = "gray40") +
      annotate("text", x = max(dat$Time_days)*0.7, y = 13.5,
               label = "Normal AM cortisol (~12 μg/dL)", size=3) +
      labs(x="Time (days)", y="Cortisol (μg/dL)", title="Endogenous Cortisol (HPA Suppression)") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_flare <- renderPlotly({
    dat <- sim_data()
    p <- ggplot(dat, aes(x = Time_days, y = FLARE * 100)) +
      geom_area(fill = "#F44336", alpha = 0.3) +
      geom_line(color = "#F44336", linewidth = 0.9) +
      labs(x="Time (days)", y="Relapse Risk Score (%)", title="Flare/Relapse Risk") +
      theme_bw()
    ggplotly(p)
  })

  ## ── Scenario comparison plots ────────────────────────────────────────────────
  sc_data_filtered <- reactive({
    dat <- all_scenarios()
    sel <- input$selected_scenarios
    id_map <- c(S1="S1: No Treatment", S2="S2: Pred 15mg Std",
                S3="S3: Pred 22.5mg Rapid", S4="S4: Pred 15mg Slow",
                S5="S5: TCZ QW + Pred", S6="S6: TCZ Q2W + Pred",
                S7="S7: TCZ Only (GC-Free)")
    dat %>% filter(ScenarioID %in% sel)
  })

  col_map_react <- reactive({
    c("S1: No Treatment"     = "#999999",
      "S2: Pred 15mg Std"    = "#2196F3",
      "S3: Pred 22.5mg Rapid"= "#FF9800",
      "S4: Pred 15mg Slow"   = "#9C27B0",
      "S5: TCZ QW + Pred"    = "#4CAF50",
      "S6: TCZ Q2W + Pred"   = "#009688",
      "S7: TCZ Only (GC-Free)"="#F44336")
  })

  output$comp_crp <- renderPlotly({
    dat <- sc_data_filtered()
    p <- ggplot(dat, aes(x = Time_days, y = CRP, color = Scenario)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 5, linetype="dashed", color="gray40") +
      scale_color_manual(values = col_map_react()) +
      labs(x="Time (days)", y="CRP (mg/L)") +
      theme_bw() + theme(legend.position="bottom", legend.text=element_text(size=7))
    ggplotly(p)
  })

  output$comp_pmras <- renderPlotly({
    dat <- sc_data_filtered()
    p <- ggplot(dat, aes(x = Time_days, y = PMRAS, color = Scenario)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 7,    linetype="dashed", color="darkgreen") +
      geom_hline(yintercept = 17.5, linetype="dashed", color="orange") +
      scale_color_manual(values = col_map_react()) +
      labs(x="Time (days)", y="PMR-AS Score") +
      theme_bw() + theme(legend.position="bottom", legend.text=element_text(size=7))
    ggplotly(p)
  })

  output$comp_il6 <- renderPlotly({
    dat <- sc_data_filtered()
    p <- ggplot(dat, aes(x = Time_days, y = IL6, color = Scenario)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 3.4, linetype="dashed", color="gray40") +
      scale_color_manual(values = col_map_react()) +
      coord_cartesian(ylim = c(0, pmax(50, max(dat$IL6, na.rm=T) * 1.1))) +
      labs(x="Time (days)", y="IL-6 (pg/mL)") +
      theme_bw() + theme(legend.position="bottom", legend.text=element_text(size=7))
    ggplotly(p)
  })

  output$comp_bmd <- renderPlotly({
    dat <- sc_data_filtered()
    p <- ggplot(dat, aes(x = Time_days, y = BMD * 100, color = Scenario)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 90, linetype="dashed", color="red", alpha=0.6) +
      annotate("text", x = 50, y = 91, label="Osteopenia risk", size=3, color="red") +
      scale_color_manual(values = col_map_react()) +
      labs(x="Time (days)", y="BMD (% of baseline)") +
      theme_bw() + theme(legend.position="bottom", legend.text=element_text(size=7))
    ggplotly(p)
  })

  output$scenario_summary_table <- DT::renderDataTable({
    dat <- all_scenarios()
    days_to_check <- c(14, 28, 84, 182, 364)
    dat %>%
      mutate(Day = round(Time_days)) %>%
      filter(Day %in% days_to_check) %>%
      group_by(ScenarioID, Scenario, Day) %>%
      slice(1) %>%
      ungroup() %>%
      select(Scenario, Day, CRP, ESR, IL6, PMRAS, CORT, BMD) %>%
      mutate(
        CRP    = round(CRP, 1),
        ESR    = round(ESR, 1),
        IL6    = round(IL6, 1),
        PMRAS  = round(PMRAS, 1),
        CORT   = round(CORT, 1),
        BMD    = round(BMD * 100, 1)
      ) %>%
      rename(`CRP (mg/L)`=CRP, `ESR (mm/hr)`=ESR, `IL-6 (pg/mL)`=IL6,
             `PMR-AS`=PMRAS, `Cortisol (μg/dL)`=CORT, `BMD (%)`=BMD)
  }, options = list(pageLength=15, scrollX=TRUE), rownames=FALSE)

  ## ── Biomarker explorer ────────────────────────────────────────────────────────
  output$bio_bmd <- renderPlotly({
    dat <- all_scenarios()
    p <- ggplot(dat, aes(x=Time_days, y=BMD*100, color=Scenario)) +
      geom_line(linewidth=0.8) +
      geom_hline(yintercept=90, linetype="dashed", color="red") +
      geom_hline(yintercept=80, linetype="dashed", color="darkred") +
      scale_color_manual(values=col_map_react()) +
      labs(x="Time (days)", y="BMD (%)", title="Bone Mineral Density") +
      theme_bw() + theme(legend.position="bottom", legend.text=element_text(size=7))
    ggplotly(p)
  })

  output$bio_cortisol <- renderPlotly({
    dat <- sim_data()
    p <- ggplot(dat) +
      geom_line(aes(x=Time_days, y=CORT, color="Cortisol"), linewidth=0.9) +
      geom_line(aes(x=Time_days, y=Cp_FREE_out/30, color="Free Pred (×1/30)"), linewidth=0.8) +
      scale_color_manual(values=c("Cortisol"="#FFC107", "Free Pred (×1/30)"="#2196F3"),
                         name="Marker") +
      labs(x="Time (days)", y="μg/dL / scaled ng/mL", title="HPA Axis: Cortisol & Free Prednisolone") +
      theme_bw()
    ggplotly(p)
  })

  output$bio_gc_sparing <- renderPlotly({
    dat <- all_scenarios()
    cum_gc <- dat %>%
      group_by(Scenario, ScenarioID) %>%
      arrange(Time_days) %>%
      mutate(Cum_GC = cumsum(Cp_PRED_out * 0.5)) %>%  # approx
      select(Scenario, Time_days, Cum_GC)

    p <- ggplot(cum_gc, aes(x=Time_days, y=Cum_GC, color=Scenario)) +
      geom_line(linewidth=0.8) +
      scale_color_manual(values=col_map_react()) +
      labs(x="Time (days)", y="Cumulative GC Exposure (arb. units)",
           title="GC Sparing: Cumulative Prednisolone") +
      theme_bw() + theme(legend.position="bottom", legend.text=element_text(size=7))
    ggplotly(p)
  })

  output$bio_correlation <- renderPlotly({
    dat <- sim_data()
    # Sample every 5 days for scatter
    dat_sub <- dat %>% filter(Time_days %% 5 < 0.6) %>%
      mutate(PMRAS_bin = cut(PMRAS, breaks=c(0,7,17.5,35,70),
                             labels=c("Remission","Low","Active","High")))
    p <- ggplot(dat_sub, aes(x=IL6, y=CRP, size=PMRAS, color=PMRAS_bin)) +
      geom_point(alpha=0.6) +
      scale_color_manual(values=c("Remission"="darkgreen","Low"="gold","Active"="orange","High"="red")) +
      labs(x="IL-6 (pg/mL)", y="CRP (mg/L)",
           title="Biomarker Correlation: IL-6 vs CRP (size=PMR-AS)",
           color="PMR-AS Level") +
      theme_bw()
    ggplotly(p)
  })

  output$biomarker_ref_table <- DT::renderDataTable({
    data.frame(
      Biomarker = c("CRP", "ESR", "IL-6", "PMR-AS", "sIL-6R",
                    "Cortisol (AM)", "BMD (lumbar)", "Fibrinogen",
                    "Alkaline Phosphatase", "Morning Stiffness"),
      Normal_Range = c("<5 mg/L", "<20 mm/hr (W)/<15 (M)", "<3.4 pg/mL",
                       "≤7 (remission)", "20–50 ng/mL", "6–23 μg/dL",
                       "T-score > -1.0", "200–400 mg/dL", "44–147 U/L", "<15 min"),
      Typical_PMR = c("35–100 mg/L", "55–100 mm/hr", "15–80 pg/mL",
                      "25–55 at diagnosis", "40–90 ng/mL (elevated)",
                      "Suppressed with GC", "Declines >6mo GC",
                      "400–800 mg/dL", "Mildly elevated", ">45 min"),
      Clinical_Significance = c("Primary treatment response marker",
                                  "Supports diagnosis; less specific",
                                  "Key driver of CRP and disease activity",
                                  "Composite disease activity score",
                                  "TCZ → paradoxical rise in free IL-6",
                                  "Monitor for HPA suppression",
                                  "Osteoporosis prophylaxis needed >3mo",
                                  "Acute phase; drives ESR",
                                  "Mild elevation in active PMR",
                                  "Cardinal symptom; resolves in 24–72h with GC")
    )
  }, options = list(pageLength=10, dom='t'), rownames=FALSE)

}

## ─────────────────────────────────────────────────────────────────────────────
## Launch
## ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
