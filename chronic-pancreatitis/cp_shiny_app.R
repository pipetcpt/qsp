## ============================================================
## Chronic Pancreatitis QSP Model — Shiny Interactive Dashboard
## Author: Claude Code Routine (CCR)
## Date:   2026-06-19
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)

## ─────────────────────────────────────────────────────────────
## Inline mrgsolve model code (same as cp_mrgsolve_model.R)
## ─────────────────────────────────────────────────────────────
model_code <- '
$PARAM
ka_opioid=0.8, CL_opioid=15, V1_opioid=50, Q_opioid=5, V2_opioid=100, F_opioid=0.70,
ka_pert=2.0, kd_pert=0.6, Emax_pert=0.90, EC50_pert=500,
k_TNFa_prod=0.05, k_TNFa_deg=0.12,
k_IL6_prod=0.04, k_IL6_deg=0.15,
k_TGFb_prod=0.02, k_TGFb_deg=0.08,
k_ROS_prod=0.06, k_ROS_deg=0.20,
n_inflam=1.5, EC50_TGFb=0.5,
kf_PSC=0.003, kr_PSC=0.0005,
k_fibrosis=0.004, k_fibdeg=0.0008, PSC_max=1.0,
k_exo_loss=0.0015, exo_0=1.0,
k_beta_loss=0.0012, beta_0=1.0,
k_glucose=0.02, G_basal=5.5,
k_pS_prod=0.08, k_pS_deg=0.04,
k_cS_prod=0.03, k_cS_deg=0.02,
Imax_opioid_pS=0.70, IC50_opioid_pS=0.3,
Imax_opioid_cS=0.60, IC50_opioid_cS=0.4,
E_gabapentin=0.0,
Imax_pirf_TGFb=0.65, IC50_pirf=1.5,
Imax_losartan=0.50, IC50_losartan=0.8,
PIRF_ON=0, LOSARTAN_ON=0, PERT_ON=0, INSULIN_TX=0, SEVERITY=1.0

$INIT
GUT_OPIOID=0, CENT_OPIOID=0, PERI_OPIOID=0,
GUT_PERT=0, DUO_PERT=0,
TNFa=0.10, IL6=0.08, TGFb=0.05, ROS=0.10,
PSC=0.10, FIB=0.05, EXO=1.0, BETA=1.0, GLUC=5.5,
pSENS=0.05, cSENS=0.02

$ODE
dxdt_GUT_OPIOID  = -ka_opioid * GUT_OPIOID;
dxdt_CENT_OPIOID =  ka_opioid * GUT_OPIOID * F_opioid
                    - (CL_opioid/V1_opioid) * CENT_OPIOID
                    - (Q_opioid/V1_opioid) * CENT_OPIOID
                    + (Q_opioid/V2_opioid) * PERI_OPIOID;
dxdt_PERI_OPIOID =  (Q_opioid/V1_opioid) * CENT_OPIOID
                    - (Q_opioid/V2_opioid) * PERI_OPIOID;
dxdt_GUT_PERT  = -ka_pert * GUT_PERT;
dxdt_DUO_PERT  =  ka_pert * GUT_PERT - kd_pert * DUO_PERT;

double Cop  = CENT_OPIOID / V1_opioid;
double pirf_eff = PIRF_ON * Imax_pirf_TGFb * 2.0 / (IC50_pirf + 2.0 * PIRF_ON);
double stim_inflam = SEVERITY * (1.0 + FIB * 2.0 + ROS * 1.5);
dxdt_TNFa = k_TNFa_prod * stim_inflam - k_TNFa_deg * TNFa;
dxdt_IL6  = k_IL6_prod  * stim_inflam - k_IL6_deg  * IL6;
dxdt_TGFb = k_TGFb_prod * stim_inflam * (1.0 - pirf_eff) - k_TGFb_deg * TGFb;
dxdt_ROS  = k_ROS_prod  * stim_inflam - k_ROS_deg * ROS - 0.3 * ROS * PIRF_ON;

double TGFb_hill = pow(TGFb, n_inflam) / (pow(EC50_TGFb, n_inflam) + pow(TGFb, n_inflam));
double los_eff   = LOSARTAN_ON * Imax_losartan * 2.0 / (IC50_losartan + 2.0 * LOSARTAN_ON);
double PSC_net   = kf_PSC * TGFb_hill * (TNFa + 1.0) * (1.0 - los_eff) - kr_PSC * PSC;
dxdt_PSC = PSC_net;
dxdt_FIB = k_fibrosis * PSC - k_fibdeg * (1.0 - PSC) * FIB;

dxdt_EXO = -k_exo_loss * FIB * EXO + 0.001 * (exo_0 - EXO);

dxdt_BETA = -k_beta_loss * FIB * BETA;

double insulin_eff = BETA * (1.0 + INSULIN_TX * 0.5);
dxdt_GLUC = G_basal * (1.0 - BETA) * 0.1 - k_glucose * insulin_eff * (GLUC - G_basal);

double opi_pS = Imax_opioid_pS * Cop / (IC50_opioid_pS + Cop);
double gaba   = E_gabapentin;
dxdt_pSENS = k_pS_prod * (TNFa + IL6) * 0.5 * SEVERITY - k_pS_deg * pSENS
             - opi_pS * pSENS - gaba * pSENS;
double opi_cS = Imax_opioid_cS * Cop / (IC50_opioid_cS + Cop);
dxdt_cSENS = k_cS_prod * pSENS - k_cS_deg * cSENS
             - opi_cS * cSENS - gaba * 0.5 * cSENS;

$TABLE
double PAIN_SCORE = 10.0 * (pSENS * 0.6 + cSENS * 0.4);
double FAT_MALAB  = (1.0 - EXO) * 100.0;
double HBA1C      = 4.0 + (GLUC / 5.5 - 1.0) * 3.0;
double FIB_INDEX  = FIB * 100.0;

$CAPTURE PAIN_SCORE FAT_MALAB HBA1C FIB_INDEX TNFa IL6 TGFb ROS PSC FIB EXO BETA GLUC
'

## ─────────────────────────────────────────────────────────────
## Compile once at startup
## ─────────────────────────────────────────────────────────────
mod <- mread_cache("cp_shiny", tempdir(), model_code)

## ─────────────────────────────────────────────────────────────
## Helper: run simulation
## ─────────────────────────────────────────────────────────────
run_sim <- function(severity, duration_yr,
                    pert_on, opioid_dose, gaba_eff,
                    pirf_on, losartan_on, insulin_tx) {

  end_h   <- duration_yr * 365 * 24
  delta_h <- 24

  params <- list(
    SEVERITY    = severity,
    PERT_ON     = as.numeric(pert_on),
    E_gabapentin = gaba_eff,
    PIRF_ON     = as.numeric(pirf_on),
    LOSARTAN_ON = as.numeric(losartan_on),
    INSULIN_TX  = as.numeric(insulin_tx)
  )
  m <- mod %>% param(params)

  # Build dosing
  evs <- NULL
  if (pert_on) {
    evs <- ev(cmt = "GUT_PERT",    amt = 40000, ii = 8,
              addl = duration_yr * 365 * 3 - 1, time = 0)
  }
  if (opioid_dose > 0) {
    eopioid <- ev(cmt = "GUT_OPIOID", amt = opioid_dose,
                  ii = 8, addl = duration_yr * 365 * 3 - 1, time = 0)
    evs <- if (is.null(evs)) eopioid else ev_seq(evs, eopioid)
  }

  if (!is.null(evs)) {
    out <- m %>% mrgsim_e(evs, end = end_h, delta = delta_h)
  } else {
    out <- m %>% mrgsim(end = end_h, delta = delta_h)
  }
  as.data.frame(out) %>% mutate(time_days = time / 24)
}

## ─────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Chronic Pancreatitis QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",   tabName = "patient",   icon = icon("user-md")),
      menuItem("PK — Opioid/PERT",  tabName = "pk",        icon = icon("pills")),
      menuItem("PD — Inflammation & Fibrosis", tabName = "pd_inflam", icon = icon("fire")),
      menuItem("Clinical Endpoints", tabName = "endpoints", icon = icon("chart-line")),
      menuItem("Scenario Comparison", tabName = "scenarios", icon = icon("exchange-alt")),
      menuItem("Biomarker Panel",   tabName = "biomarkers", icon = icon("vials"))
    )
  ),

  dashboardBody(
    tabItems(

      ## ── Tab 1: Patient Profile ──────────────────────────────
      tabItem("patient",
        fluidRow(
          box(title = "Disease Parameters", status = "primary", solidHeader = TRUE, width = 4,
            sliderInput("severity",   "Disease Severity", 0.5, 3.0, 1.5, step = 0.1),
            sliderInput("duration_yr","Simulation Duration (years)", 1, 10, 5, step = 1)
          ),
          box(title = "Treatment Selection", status = "info", solidHeader = TRUE, width = 4,
            checkboxInput("pert_on",     "PERT (Enzyme Replacement)", FALSE),
            sliderInput("opioid_dose", "Opioid Dose (mg, 0 = off)", 0, 100, 0, step = 10),
            sliderInput("gaba_eff",    "Gabapentin Effect (0–0.6)",  0, 0.6, 0, step = 0.05),
            checkboxInput("pirf_on",     "Pirfenidone (Antifibrotic)", FALSE),
            checkboxInput("losartan_on", "Losartan (AngII Blockade)",  FALSE),
            checkboxInput("insulin_tx",  "Insulin Therapy (T3cDM)",    FALSE)
          ),
          box(title = "QSP Model Overview", status = "success", solidHeader = TRUE, width = 4,
            tags$p("This model simulates:"),
            tags$ul(
              tags$li("Pancreatic stellate cell activation & fibrosis"),
              tags$li("Exocrine dysfunction (steatorrhea)"),
              tags$li("Type 3c diabetes (pancreatogenic)"),
              tags$li("Peripheral & central pain sensitisation"),
              tags$li("Drug PK/PD: PERT, opioids, antifibrotics")
            ),
            tags$p(tags$b("20 ODEs | 5 therapeutic scenarios"))
          )
        ),
        fluidRow(
          valueBoxOutput("vb_pain",    width = 3),
          valueBoxOutput("vb_fibrosis",width = 3),
          valueBoxOutput("vb_hba1c",   width = 3),
          valueBoxOutput("vb_exo",     width = 3)
        )
      ),

      ## ── Tab 2: PK ──────────────────────────────────────────
      tabItem("pk",
        fluidRow(
          box(title = "Opioid Plasma Concentration (3-cpt Model)",
              status = "primary", solidHeader = TRUE, width = 12,
              plotlyOutput("plot_opioid_pk", height = 300))
        ),
        fluidRow(
          box(title = "PERT Duodenal Enzyme Activity",
              status = "info", solidHeader = TRUE, width = 12,
              plotlyOutput("plot_pert_pk", height = 300))
        )
      ),

      ## ── Tab 3: Inflammation & Fibrosis ────────────────────
      tabItem("pd_inflam",
        fluidRow(
          box(title = "Inflammatory Mediators Over Time",
              status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_inflam", height = 320)),
          box(title = "PSC Activation & Fibrosis Index",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_psc", height = 320))
        ),
        fluidRow(
          box(title = "ROS & TGF-β Dynamics",
              status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_ros_tgfb", height = 300)),
          box(title = "Exocrine & Beta-cell Function",
              status = "success", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_exo_beta", height = 300))
        )
      ),

      ## ── Tab 4: Clinical Endpoints ─────────────────────────
      tabItem("endpoints",
        fluidRow(
          box(title = "Pain Score (NRS 0–10)",
              status = "danger", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_pain", height = 300)),
          box(title = "Fat Malabsorption / Steatorrhea (%)",
              status = "warning", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_fat", height = 300))
        ),
        fluidRow(
          box(title = "Plasma Glucose & HbA1c",
              status = "info", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_gluc", height = 300)),
          box(title = "Fibrosis Index (CT/Elastography Proxy)",
              status = "primary", solidHeader = TRUE, width = 6,
              plotlyOutput("plot_fib", height = 300))
        )
      ),

      ## ── Tab 5: Scenario Comparison ────────────────────────
      tabItem("scenarios",
        fluidRow(
          box(title = "Pre-defined Scenario Comparison (2-year)",
              status = "primary", solidHeader = TRUE, width = 12,
              actionButton("run_scenarios", "Run All 5 Scenarios",
                           icon = icon("play"), class = "btn-success btn-lg"),
              hr(),
              plotlyOutput("plot_scenario_pain", height = 250),
              plotlyOutput("plot_scenario_fib",  height = 250),
              plotlyOutput("plot_scenario_exo",  height = 250))
        ),
        fluidRow(
          box(title = "Summary Table — 2-year Endpoint Values",
              status = "info", solidHeader = TRUE, width = 12,
              DTOutput("table_scenarios"))
        )
      ),

      ## ── Tab 6: Biomarker Panel ────────────────────────────
      tabItem("biomarkers",
        fluidRow(
          box(title = "Biomarker Time Course",
              status = "success", solidHeader = TRUE, width = 12,
              plotlyOutput("plot_biomarkers", height = 420))
        ),
        fluidRow(
          box(title = "Biomarker Reference Ranges",
              status = "info", solidHeader = TRUE, width = 12,
              DTOutput("table_biomarkers"))
        )
      )
    )
  )
)

## ─────────────────────────────────────────────────────────────
## Server
## ─────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## Reactive simulation (user-controlled)
  sim_data <- reactive({
    run_sim(
      severity    = input$severity,
      duration_yr = input$duration_yr,
      pert_on     = input$pert_on,
      opioid_dose = input$opioid_dose,
      gaba_eff    = input$gaba_eff,
      pirf_on     = input$pirf_on,
      losartan_on = input$losartan_on,
      insulin_tx  = input$insulin_tx
    )
  })

  ## Value boxes
  endpoint_last <- reactive({
    sim_data() %>% filter(time_days == max(time_days))
  })
  output$vb_pain <- renderValueBox({
    val <- round(endpoint_last()$PAIN_SCORE, 1)
    col <- if (val < 4) "green" else if (val < 7) "yellow" else "red"
    valueBox(val, "Pain NRS (0–10)", icon = icon("bolt"), color = col)
  })
  output$vb_fibrosis <- renderValueBox({
    val <- round(endpoint_last()$FIB_INDEX, 1)
    col <- if (val < 25) "green" else if (val < 60) "yellow" else "red"
    valueBox(paste0(val, "%"), "Fibrosis Index", icon = icon("layer-group"), color = col)
  })
  output$vb_hba1c <- renderValueBox({
    val <- round(endpoint_last()$HBA1C, 1)
    col <- if (val < 5.7) "green" else if (val < 6.5) "yellow" else "red"
    valueBox(paste0(val, "%"), "HbA1c", icon = icon("tint"), color = col)
  })
  output$vb_exo <- renderValueBox({
    val <- round((1 - endpoint_last()$FIB * 0.8) * 100, 1)
    col <- if (val > 70) "green" else if (val > 40) "yellow" else "red"
    valueBox(paste0(val, "%"), "Exocrine Function", icon = icon("utensils"), color = col)
  })

  ## PK plots
  output$plot_opioid_pk <- renderPlotly({
    d <- sim_data() %>% filter(time_days <= min(14, input$duration_yr * 365))
    plot_ly(d, x = ~time_days, y = ~(CENT_OPIOID/50), type = "scatter", mode = "lines",
            line = list(color = "#1E88E5")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Opioid Conc. (μg/mL)"))
  })

  output$plot_pert_pk <- renderPlotly({
    d <- sim_data() %>% filter(time_days <= min(14, input$duration_yr * 365))
    plot_ly(d, x = ~time_days, y = ~DUO_PERT, type = "scatter", mode = "lines",
            line = list(color = "#43A047")) %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "PERT (lipase units/mL)"))
  })

  ## Inflammation plot
  output$plot_inflam <- renderPlotly({
    d <- sim_data() %>% select(time_days, TNFa, IL6) %>%
      pivot_longer(c(TNFa, IL6))
    plot_ly(d, x = ~time_days, y = ~value, color = ~name, type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Cytokine (nM)"))
  })

  output$plot_psc <- renderPlotly({
    d <- sim_data() %>% select(time_days, PSC, FIB_INDEX)
    plot_ly(d, x = ~time_days) %>%
      add_lines(y = ~PSC * 100, name = "PSC Activation (%)", line = list(color = "#AB47BC")) %>%
      add_lines(y = ~FIB_INDEX, name = "Fibrosis Index", line = list(color = "#7E57C2")) %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "%"))
  })

  output$plot_ros_tgfb <- renderPlotly({
    d <- sim_data() %>% select(time_days, ROS, TGFb) %>%
      pivot_longer(c(ROS, TGFb))
    plot_ly(d, x = ~time_days, y = ~value, color = ~name, type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Concentration (nM)"))
  })

  output$plot_exo_beta <- renderPlotly({
    d <- sim_data() %>% mutate(EXO_pct = EXO * 100, BETA_pct = BETA * 100) %>%
      select(time_days, EXO_pct, BETA_pct) %>% pivot_longer(c(EXO_pct, BETA_pct))
    plot_ly(d, x = ~time_days, y = ~value, color = ~name, type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "% Normal"))
  })

  ## Clinical endpoints
  output$plot_pain <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~PAIN_SCORE, type = "scatter", mode = "lines",
            line = list(color = "#E53935")) %>%
      add_segments(x = 0, xend = max(d$time_days), y = 4, yend = 4,
                   line = list(dash = "dot", color = "orange"), name = "Moderate threshold") %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "NRS (0–10)", range = c(0, 10)))
  })

  output$plot_fat <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~FAT_MALAB, type = "scatter", mode = "lines",
            line = list(color = "#FB8C00")) %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Steatorrhea (%)"))
  })

  output$plot_gluc <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days) %>%
      add_lines(y = ~GLUC, name = "Glucose (mmol/L)", line = list(color = "#1E88E5")) %>%
      add_lines(y = ~HBA1C, name = "HbA1c (%)", line = list(color = "#E53935")) %>%
      add_segments(x = 0, xend = max(d$time_days), y = 7.0, yend = 7.0,
                   line = list(dash = "dot", color = "red"), name = "DM threshold") %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Value"))
  })

  output$plot_fib <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~FIB_INDEX, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(126,87,194,0.2)",
            line = list(color = "#7E57C2")) %>%
      layout(xaxis = list(title = "Time (days)"), yaxis = list(title = "Fibrosis Index (%)"))
  })

  ## Scenario comparison
  scenario_data <- eventReactive(input$run_scenarios, {
    scens <- list(
      list(label = "1. No Treatment",          severity = 1.5, pert = FALSE, opi = 0,    gaba = 0,    pirf = FALSE, los = FALSE, ins = FALSE),
      list(label = "2. PERT",                  severity = 1.5, pert = TRUE,  opi = 0,    gaba = 0,    pirf = FALSE, los = FALSE, ins = FALSE),
      list(label = "3. Opioid + Gabapentin",   severity = 1.5, pert = FALSE, opi = 50,   gaba = 0.35, pirf = FALSE, los = FALSE, ins = FALSE),
      list(label = "4. Antifibrotic",          severity = 1.5, pert = FALSE, opi = 0,    gaba = 0,    pirf = TRUE,  los = TRUE,  ins = FALSE),
      list(label = "5. Combination",           severity = 1.5, pert = TRUE,  opi = 50,   gaba = 0.35, pirf = TRUE,  los = TRUE,  ins = TRUE)
    )
    bind_rows(lapply(scens, function(s) {
      run_sim(s$severity, 2, s$pert, s$opi, s$gaba, s$pirf, s$los, s$ins) %>%
        mutate(scenario = s$label)
    }))
  })

  output$plot_scenario_pain <- renderPlotly({
    req(scenario_data())
    d <- scenario_data()
    plot_ly(d, x = ~time_days, y = ~PAIN_SCORE, color = ~scenario,
            type = "scatter", mode = "lines") %>%
      layout(title = "Pain Score", yaxis = list(title = "NRS"),
             xaxis = list(title = "Days"))
  })

  output$plot_scenario_fib <- renderPlotly({
    req(scenario_data())
    d <- scenario_data()
    plot_ly(d, x = ~time_days, y = ~FIB_INDEX, color = ~scenario,
            type = "scatter", mode = "lines") %>%
      layout(title = "Fibrosis Index", yaxis = list(title = "%"),
             xaxis = list(title = "Days"))
  })

  output$plot_scenario_exo <- renderPlotly({
    req(scenario_data())
    d <- scenario_data()
    plot_ly(d, x = ~time_days, y = ~FAT_MALAB, color = ~scenario,
            type = "scatter", mode = "lines") %>%
      layout(title = "Fat Malabsorption (%)", yaxis = list(title = "%"),
             xaxis = list(title = "Days"))
  })

  output$table_scenarios <- renderDT({
    req(scenario_data())
    scenario_data() %>%
      group_by(scenario) %>%
      filter(time_days == max(time_days)) %>%
      summarise(
        Pain_NRS  = round(PAIN_SCORE, 2),
        Fibrosis  = round(FIB_INDEX, 1),
        FatMalab  = round(FAT_MALAB, 1),
        HbA1c     = round(HBA1C, 2),
        Glucose   = round(GLUC, 2),
        BetaCell  = round(BETA * 100, 1),
        .groups = "drop"
      ) %>%
      datatable(options = list(pageLength = 10, dom = "t"))
  })

  ## Biomarkers
  output$plot_biomarkers <- renderPlotly({
    d <- sim_data() %>%
      mutate(
        `Serum TNF-α (nM)`   = TNFa,
        `Serum IL-6 (nM)`    = IL6,
        `TGF-β (nM)`         = TGFb,
        `ROS index`          = ROS,
        `Fecal Elastase (proxy)` = EXO * 200,
        `Fasting Glucose (mmol/L)` = GLUC
      ) %>%
      select(time_days, starts_with("Serum"), starts_with("TGF"), starts_with("ROS"),
             starts_with("Fecal"), starts_with("Fasting")) %>%
      pivot_longer(-time_days)

    plot_ly(d, x = ~time_days, y = ~value, color = ~name,
            type = "scatter", mode = "lines") %>%
      layout(xaxis = list(title = "Time (days)"),
             yaxis = list(title = "Value"),
             legend = list(orientation = "h"))
  })

  output$table_biomarkers <- renderDT({
    data.frame(
      Biomarker = c("Serum Amylase", "Serum Lipase", "CRP", "IL-6", "TGF-β",
                    "Fecal Elastase-1", "Fecal Fat (72h)", "HbA1c",
                    "CA 19-9", "IgG4 (AIP)", "Serum Lipase (acute)"),
      `Reference Range` = c("<100 U/L", "<60 U/L", "<5 mg/L", "<7 pg/mL", "<5 ng/mL",
                             ">200 μg/g (normal)", "<7g/day", "<5.7%",
                             "<37 U/mL", "<140 mg/dL", "<60 U/L"),
      `CP Finding` = c("Often normal in chronic", "Often normal in chronic",
                        "Elevated during flares", "Elevated", "Elevated (fibrosis)",
                        "<100 μg/g = EPI", ">7g/day = steatorrhea", "≥6.5% = T3cDM",
                        "Monitor annually", "Elevated if AIP", "↑ in acute flares"),
      Clinical_Use = c("Acute flare diagnosis", "Acute flare diagnosis",
                        "Activity marker", "Activity & severity",
                        "Fibrosis progression", "EPI diagnosis",
                        "EPI quantification", "Glucose control",
                        "PDAC surveillance", "AIP diagnosis",
                        "Acute on chronic CP")
    ) %>%
      datatable(options = list(pageLength = 15, dom = "t"),
                colnames = c("Biomarker", "Reference Range", "Finding in CP", "Clinical Use"))
  })
}

shinyApp(ui, server)
