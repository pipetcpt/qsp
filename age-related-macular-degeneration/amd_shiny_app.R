## ============================================================
## AMD (Age-related Macular Degeneration) — Shiny QSP Dashboard
## ============================================================
## 6 Tabs:
##   1. Patient Profile & Disease Stage
##   2. Drug PK — Vitreous & Retinal Concentrations
##   3. PD Key Markers — VEGF, VEGFR2, Ang-2
##   4. Clinical Endpoints — BCVA, CST, CNV, GA
##   5. Scenario Comparison — Multi-drug Regimens
##   6. Biomarker Explorer — Complement, RPE, Photoreceptor
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mrgsolve)

## ---- inline mrgsolve model code (abbreviated for Shiny) ----
amd_code <- '
$PARAM
Dose_vit=0.5 Fabs_ret=0.3 kel_vit=0.0965 kel_sys=0.693
ktr_vit2ret=0.02 ktr_vit2sys=0.005 ktr_ret2sys=0.01 V_vit=0.0046
k_VEGF_syn=0.5 k_VEGF_deg=0.2 kon_VEGF=5.0 koff_VEGF=0.002
Kd_VEGF=0.04 VEGF_baseline=2.5
k_R2_act=0.5 k_R2_inact=0.3 EC50_R2=0.5 Hill_R2=1.2
ANG2_baseline=1.0 k_ANG2_syn=0.15 k_ANG2_deg=0.1
kon_ANG2=2.0 koff_ANG2=0.00087 ANG2_dual=0
C3_baseline=5.0 k_C3_syn=0.3 k_C3_deg=0.1
k_C5_syn=0.15 k_C5_deg=0.08 k_MAC_form=0.05 k_MAC_deg=0.1
RPE0=1.0 k_RPE_death=0.0002 E_RPE_ox=0.0005 E_RPE_MAC2=0.0008
k_RPE_repair=0.00005 k_LF_accum=0.001 k_LF_clear=0.0001
k_Drusen_grow=0.003 k_Drusen_base=0.0003
k_CNV_init=0.0005 E_VEGF_CNV=0.002 EC50_VEGF_CNV=1.0
k_CNV_reg=0.01 CNV_max=25.0
k_Fluid_in=10.0 k_Fluid_out=0.15
k_GA_grow=0.00055 E_MAC_GA=0.0003 E_RPE_GA=0.001
k_PR_death=0.0001 E_RPE_PR=0.002 E_fluid_PR=0.00005
k_BCVA_CNV=0.8 k_BCVA_GA=1.2 k_BCVA_fluid=0.02 k_BCVA_PR=15.0
BCVA_max=85.0 WET_AMD=1 koff_eff=0.0039

$CMT DRUG_VIT DRUG_RET DRUG_SYS VEGF_FREE VEGF_BOUND VEGFR2_ACT
ANG2_FREE ANG2_BOUND C3_LOCAL C5_LOCAL MAC_LOCAL
RPE_NORM RPE_DAM LIPOFUSCIN DRUSEN CNV_AREA FLUID_EX GA_AREA
BCVA_SCORE PR_FRAC

$INIT
DRUG_VIT=0 DRUG_RET=0 DRUG_SYS=0 VEGF_FREE=2.5 VEGF_BOUND=0
VEGFR2_ACT=0.3 ANG2_FREE=1.0 ANG2_BOUND=0 C3_LOCAL=5.0 C5_LOCAL=2.0
MAC_LOCAL=0.5 RPE_NORM=1.0 RPE_DAM=0.0 LIPOFUSCIN=0.2 DRUSEN=0.8
CNV_AREA=2.0 FLUID_EX=120.0 GA_AREA=0.0 BCVA_SCORE=55.0 PR_FRAC=0.95

$ODE
double kon_used = koff_eff / Kd_VEGF;
double R_bind_VEGF = kon_used * DRUG_RET * VEGF_FREE - koff_eff * VEGF_BOUND;
double R_bind_ANG2 = ANG2_dual * (kon_ANG2 * DRUG_RET * ANG2_FREE - koff_ANG2 * ANG2_BOUND);

dxdt_DRUG_VIT = -kel_vit * DRUG_VIT - ktr_vit2ret * DRUG_VIT - ktr_vit2sys * DRUG_VIT;
dxdt_DRUG_RET = ktr_vit2ret * DRUG_VIT - ktr_ret2sys * DRUG_RET
                - kon_used * DRUG_RET * VEGF_FREE + koff_eff * VEGF_BOUND - R_bind_ANG2;
dxdt_DRUG_SYS = ktr_vit2sys * DRUG_VIT * V_vit + ktr_ret2sys * DRUG_RET - kel_sys * DRUG_SYS;

double VEGF_up = 1.0 + 0.5 * VEGFR2_ACT + 0.3 * RPE_DAM;
dxdt_VEGF_FREE = k_VEGF_syn * VEGF_up - k_VEGF_deg * VEGF_FREE - R_bind_VEGF;
dxdt_VEGF_BOUND = R_bind_VEGF - (k_VEGF_deg + kel_vit) * VEGF_BOUND;

double Hill_n = pow(VEGF_FREE, Hill_R2);
double Hill_d = pow(EC50_R2, Hill_R2) + Hill_n;
dxdt_VEGFR2_ACT = k_R2_act * (Hill_n / Hill_d - VEGFR2_ACT);

dxdt_ANG2_FREE = k_ANG2_syn * (1.0 + 0.3 * VEGFR2_ACT) - k_ANG2_deg * ANG2_FREE - R_bind_ANG2;
dxdt_ANG2_BOUND = R_bind_ANG2 - k_ANG2_deg * ANG2_BOUND;

dxdt_C3_LOCAL = k_C3_syn * (1.0 + 0.2 * LIPOFUSCIN) - k_C3_deg * C3_LOCAL;
dxdt_C5_LOCAL = k_C5_syn * (1.0 + 0.15 * C3_LOCAL) - k_C5_deg * C5_LOCAL;
dxdt_MAC_LOCAL = k_MAC_form * C5_LOCAL - k_MAC_deg * MAC_LOCAL;

double RPE_ox = E_RPE_ox * (MAC_LOCAL * E_RPE_MAC2 + LIPOFUSCIN * 0.5);
double RPE_rate = k_RPE_death + RPE_ox + E_RPE_MAC2 * MAC_LOCAL;
dxdt_RPE_NORM = -RPE_rate * RPE_NORM + k_RPE_repair * (1.0 - RPE_NORM - RPE_DAM);
dxdt_RPE_DAM  = RPE_rate * RPE_NORM * 0.5 - 0.05 * RPE_DAM;

dxdt_LIPOFUSCIN = k_LF_accum * (2.0 - RPE_NORM) - k_LF_clear * LIPOFUSCIN;
dxdt_DRUSEN = (k_Drusen_base + k_Drusen_grow * (1.0 - RPE_NORM) * LIPOFUSCIN) * (1.0 - DRUSEN/20.0);

double VEGF_eff = VEGF_FREE / (EC50_VEGF_CNV + VEGF_FREE);
double CNV_g = (k_CNV_init + E_VEGF_CNV * VEGF_eff) * (1.0 - CNV_AREA/CNV_max) * WET_AMD;
double CNV_r = k_CNV_reg * (1.0 - VEGFR2_ACT) * CNV_AREA;
dxdt_CNV_AREA = (CNV_g - CNV_r) * CNV_AREA;

double Fluid_in = k_Fluid_in * VEGFR2_ACT * (1.0 + 0.5 * ANG2_FREE/(1.0+ANG2_FREE));
dxdt_FLUID_EX = Fluid_in - k_Fluid_out * FLUID_EX;

dxdt_GA_AREA = (k_GA_grow + E_MAC_GA * MAC_LOCAL + E_RPE_GA * fmax(0.0, 1.0 - RPE_NORM))
               * (1.0 + GA_AREA/5.0);

double PR_loss = E_RPE_PR * fmax(0.0, RPE_rate - k_RPE_death)
               + E_fluid_PR * fmax(0.0, FLUID_EX - 50.0);
dxdt_PR_FRAC = -k_PR_death * PR_FRAC - PR_loss;

double BCVA_t = fmax(0.0, fmin(BCVA_max,
    BCVA_max - k_BCVA_CNV * CNV_AREA - k_BCVA_GA * GA_AREA * 0.5
    - k_BCVA_fluid * fmax(0.0, FLUID_EX - 30.0) - k_BCVA_PR * (1.0 - PR_FRAC)));
dxdt_BCVA_SCORE = 0.05 * (BCVA_t - BCVA_SCORE);

$TABLE
capture CST = 280.0 + FLUID_EX;
capture BCVA_change = BCVA_SCORE - 55.0;

$CAPTURE DRUG_VIT DRUG_RET VEGF_FREE VEGF_BOUND VEGFR2_ACT ANG2_FREE
CNV_AREA CST FLUID_EX GA_AREA BCVA_SCORE RPE_NORM PR_FRAC
MAC_LOCAL DRUSEN C3_LOCAL C5_LOCAL BCVA_change
'

## Helper to create injection events
make_events_shiny <- function(dose_nM, load_times, maint_start,
                               maint_interval, n_maint) {
  maint_times <- seq(maint_start,
                     maint_start + (n_maint - 1) * maint_interval,
                     by = maint_interval)
  all_times <- c(load_times, maint_times)
  ev(time = all_times, amt = dose_nM, cmt = 1, addl = 0)
}

## Compute dose_nM from drug parameters
drug_params <- data.frame(
  name  = c("Ranibizumab", "Aflibercept", "Bevacizumab", "Faricimab", "Brolucizumab"),
  dose  = c(0.5, 2.0, 1.25, 6.0, 6.0),
  MW    = c(48,  115,  149,  150,  26),
  Kd    = c(0.04, 0.0005, 0.2, 0.0003, 0.06),
  koff  = c(0.0039, 0.0005, 0.002, 0.0003, 0.006),
  kel   = c(0.0965, 0.099, 0.077, 0.099, 0.173),
  color = c("#2980B9","#E74C3C","#27AE60","#8E44AD","#F39C12"),
  ang2_dual = c(0, 0, 0, 1, 0),
  stringsAsFactors = FALSE
)

## ---- UI ----
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "AMD QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "profile",   icon = icon("user-md")),
      menuItem("Drug PK",            tabName = "pk",        icon = icon("syringe")),
      menuItem("PD Key Markers",     tabName = "pd",        icon = icon("chart-line")),
      menuItem("Clinical Endpoints", tabName = "endpoints", icon = icon("eye")),
      menuItem("Scenario Comparison",tabName = "scenarios", icon = icon("code-compare")),
      menuItem("Biomarker Explorer", tabName = "biomarkers",icon = icon("microscope"))
    ),
    hr(),
    h5("Global Settings", style = "margin-left:15px;"),
    selectInput("drug_select", "Drug:",
                choices = c("Ranibizumab","Aflibercept","Bevacizumab","Faricimab","Brolucizumab"),
                selected = "Ranibizumab"),
    selectInput("amd_type", "AMD Type:",
                choices = c("Wet AMD (CNV)" = 1, "Dry AMD (GA)" = 0),
                selected = 1),
    sliderInput("sim_duration", "Simulation (days):", 180, 1460, 730, step = 90),
    sliderInput("n_loading", "Loading doses (#):", 1, 6, 3),
    sliderInput("load_interval", "Loading interval (days):", 14, 56, 28, step = 7),
    sliderInput("maint_interval", "Maintenance interval (days):", 28, 168, 56, step = 14),
    sliderInput("n_maint", "Maintenance doses (#):", 3, 20, 10),
    sliderInput("bcva_baseline", "Baseline BCVA (letters):", 20, 80, 55),
    sliderInput("cnv_init", "Initial CNV area (mm²):", 0, 15, 2.0, step = 0.5),
    sliderInput("fluid_init", "Initial excess fluid (μm):", 0, 300, 120, step = 10)
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .box-header { font-weight: bold; }
      .value-box .value { font-size: 28px; }
    "))),
    tabItems(
      ## TAB 1: Patient Profile
      tabItem(tabName = "profile",
        fluidRow(
          valueBoxOutput("vb_bcva_bl"),
          valueBoxOutput("vb_amd_type"),
          valueBoxOutput("vb_cnv_init")
        ),
        fluidRow(
          box(title = "AMD Disease Stage Classification", width = 6, status = "primary",
            plotOutput("plot_disease_stage", height = "280px")
          ),
          box(title = "Baseline Parameters Summary", width = 6, status = "info",
            tableOutput("tbl_baseline")
          )
        ),
        fluidRow(
          box(title = "Key Risk Factors in AMD", width = 12, status = "warning",
            p("Genetic: CFH Y402H (rs1061170), ARMS2 A69S (rs10490924) — ~50% attributable risk"),
            p("Environmental: Cigarette smoking (RR ~4×), UV/blue light, oxidant-poor diet"),
            p("Demographic: Age >65yr, female sex, family history, cardiovascular risk factors"),
            p("Drusen burden: Large drusen (>125μm), soft confluent drusen = intermediate AMD → highest progression risk")
          )
        )
      ),
      ## TAB 2: Drug PK
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Vitreous Drug Concentration vs Time", width = 6, status = "primary",
            plotOutput("plot_pk_vit", height = "300px")
          ),
          box(title = "Retinal/RPE Drug Concentration vs Time", width = 6, status = "primary",
            plotOutput("plot_pk_ret", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Free vs Bound VEGF", width = 6, status = "info",
            plotOutput("plot_vegf_bound", height = "280px")
          ),
          box(title = "PK Summary at Key Timepoints", width = 6, status = "info",
            tableOutput("tbl_pk_summary")
          )
        )
      ),
      ## TAB 3: PD Key Markers
      tabItem(tabName = "pd",
        fluidRow(
          box(title = "VEGFR-2 Activation", width = 6, status = "warning",
            plotOutput("plot_vegfr2", height = "300px")
          ),
          box(title = "Free VEGF Suppression", width = 6, status = "warning",
            plotOutput("plot_vegf_free", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Ang-2 Free Levels (Faricimab dual target)", width = 6, status = "info",
            plotOutput("plot_ang2", height = "280px")
          ),
          box(title = "VEGF Suppression % vs Target BCVA Gain", width = 6, status = "success",
            plotOutput("plot_pd_bcva_corr", height = "280px")
          )
        )
      ),
      ## TAB 4: Clinical Endpoints
      tabItem(tabName = "endpoints",
        fluidRow(
          valueBoxOutput("vb_bcva_1yr"),
          valueBoxOutput("vb_cst_1yr"),
          valueBoxOutput("vb_cnv_1yr")
        ),
        fluidRow(
          box(title = "BCVA over Time (ETDRS letters)", width = 6, status = "success",
            plotOutput("plot_bcva", height = "300px")
          ),
          box(title = "OCT Central Subfield Thickness (CST)", width = 6, status = "primary",
            plotOutput("plot_cst", height = "300px")
          )
        ),
        fluidRow(
          box(title = "CNV Lesion Area", width = 6, status = "danger",
            plotOutput("plot_cnv", height = "280px")
          ),
          box(title = "Geographic Atrophy Area (Late Dry AMD)", width = 6, status = "warning",
            plotOutput("plot_ga", height = "280px")
          )
        )
      ),
      ## TAB 5: Scenario Comparison
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "BCVA Comparison — All Regimens", width = 12, status = "primary",
            plotOutput("plot_scenario_bcva", height = "340px")
          )
        ),
        fluidRow(
          box(title = "CNV Area Comparison", width = 6, status = "danger",
            plotOutput("plot_scenario_cnv", height = "280px")
          ),
          box(title = "Injection Burden vs. BCVA Gain (1-year)", width = 6, status = "info",
            plotOutput("plot_scenario_burden", height = "280px")
          )
        ),
        fluidRow(
          box(title = "1-year Outcomes Table", width = 12, status = "success",
            tableOutput("tbl_scenarios")
          )
        )
      ),
      ## TAB 6: Biomarker Explorer
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "Complement Cascade — C3, C5, MAC", width = 6, status = "warning",
            plotOutput("plot_complement", height = "300px")
          ),
          box(title = "RPE Cell Health & Lipofuscin", width = 6, status = "info",
            plotOutput("plot_rpe", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Photoreceptor Survival", width = 6, status = "success",
            plotOutput("plot_pr", height = "280px")
          ),
          box(title = "Drusen Burden", width = 6, status = "primary",
            plotOutput("plot_drusen", height = "280px")
          )
        )
      )
    )
  )
)

## ---- SERVER ----
server <- function(input, output, session) {

  ## Load model once
  mod_base <- reactive({ mcode("AMD_QSP_shiny", amd_code) })

  ## Compute drug parameters
  drug_info <- reactive({
    dp <- drug_params[drug_params$name == input$drug_select, ]
    dose_nM <- (dp$dose / dp$MW) / 0.0046 * 1000
    list(dp = dp, dose_nM = dose_nM)
  })

  ## Build events
  events_main <- reactive({
    di <- drug_info()
    load_times <- seq(0, (input$n_loading - 1) * input$load_interval,
                      by = input$load_interval)
    maint_start <- max(load_times) + input$maint_interval
    make_events_shiny(di$dose_nM, load_times, maint_start,
                      input$maint_interval, input$n_maint)
  })

  ## Main simulation
  sim_main <- reactive({
    di <- drug_info()
    mod_base() %>%
      param(Kd_VEGF    = di$dp$Kd,
            koff_eff   = di$dp$koff,
            kel_vit    = di$dp$kel,
            ANG2_dual  = di$dp$ang2_dual,
            WET_AMD    = as.numeric(input$amd_type)) %>%
      init(BCVA_SCORE = input$bcva_baseline,
           CNV_AREA   = input$cnv_init,
           FLUID_EX   = input$fluid_init) %>%
      ev(events_main()) %>%
      mrgsim(end = input$sim_duration, delta = 1) %>%
      as.data.frame()
  })

  ## No-treatment simulation
  sim_noTx <- reactive({
    mod_base() %>%
      param(WET_AMD = as.numeric(input$amd_type)) %>%
      init(BCVA_SCORE = input$bcva_baseline,
           CNV_AREA   = input$cnv_init,
           FLUID_EX   = input$fluid_init) %>%
      mrgsim(end = input$sim_duration, delta = 1) %>%
      as.data.frame()
  })

  ## All-scenarios simulation (fixed 2yr, for scenario tab)
  sim_all <- reactive({
    withProgress(message = "Running all scenarios...", {
      end_t <- input$sim_duration
      results <- lapply(1:5, function(i) {
        dp <- drug_params[i, ]
        dose_nM <- (dp$dose / dp$MW) / 0.0046 * 1000
        load_times <- c(0, 28, 56)
        maint_int  <- c(56, 56, 56, 112, 84)[i]
        maint_start <- 84
        ev_i <- make_events_shiny(dose_nM, load_times, maint_start, maint_int, 10)
        out <- mod_base() %>%
          param(Kd_VEGF   = dp$Kd,
                koff_eff  = dp$koff,
                kel_vit   = dp$kel,
                ANG2_dual = dp$ang2_dual,
                WET_AMD   = 1) %>%
          ev(ev_i) %>%
          mrgsim(end = end_t, delta = 1) %>%
          as.data.frame() %>%
          mutate(scenario = dp$name)
        out
      })
      noTx <- mod_base() %>%
        param(WET_AMD = 1) %>%
        mrgsim(end = end_t, delta = 1) %>%
        as.data.frame() %>%
        mutate(scenario = "No Treatment")
      bind_rows(c(results, list(noTx)))
    })
  })

  ## ---- VALUE BOXES ----
  output$vb_bcva_bl <- renderValueBox({
    valueBox(paste0(input$bcva_baseline, " letters"),
             "Baseline BCVA", icon = icon("eye"), color = "blue")
  })
  output$vb_amd_type <- renderValueBox({
    txt <- ifelse(as.numeric(input$amd_type) == 1, "Wet AMD (CNV)", "Dry AMD (GA)")
    valueBox(txt, "AMD Type", icon = icon("diagnoses"), color = "orange")
  })
  output$vb_cnv_init <- renderValueBox({
    valueBox(paste0(input$cnv_init, " mm²"),
             "Initial CNV Area", icon = icon("circle"), color = "red")
  })
  output$vb_bcva_1yr <- renderValueBox({
    sim <- sim_main()
    v <- sim[sim$time == 365, "BCVA_SCORE"]
    if(length(v) == 0) v <- tail(sim$BCVA_SCORE, 1)
    valueBox(paste0(round(v, 1), " letters"),
             "BCVA at 1 Year", icon = icon("eye"), color = "green")
  })
  output$vb_cst_1yr <- renderValueBox({
    sim <- sim_main()
    v <- sim[sim$time == 365, "CST"]
    if(length(v) == 0) v <- tail(sim$CST, 1)
    valueBox(paste0(round(v), " μm"),
             "CST at 1 Year", icon = icon("layer-group"), color = "blue")
  })
  output$vb_cnv_1yr <- renderValueBox({
    sim <- sim_main()
    v <- sim[sim$time == 365, "CNV_AREA"]
    if(length(v) == 0) v <- tail(sim$CNV_AREA, 1)
    valueBox(paste0(round(v, 2), " mm²"),
             "CNV Area at 1 Year", icon = icon("circle"), color = "red")
  })

  ## ---- TAB 1 ----
  output$plot_disease_stage <- renderPlot({
    stages <- data.frame(
      Stage = c("No AMD", "Early AMD", "Intermediate AMD", "Late AMD (GA)", "Late AMD (Wet)"),
      Drusen_min = c(0, 0.1, 0.5, 2.0, 0),
      Drusen_max = c(0.1, 0.5, 5.0, 20, 5),
      BCVA_min   = c(75, 70, 60, 40, 30),
      BCVA_max   = c(85, 80, 75, 60, 65),
      Fill = c("#2ECC71","#F1C40F","#E67E22","#E74C3C","#C0392B")
    )
    ggplot(stages, aes(x = Drusen_max, y = BCVA_max, fill = Stage)) +
      geom_point(size = 10, shape = 21, color = "black") +
      geom_text(aes(label = Stage), vjust = -1.5, size = 3.5) +
      scale_fill_manual(values = setNames(stages$Fill, stages$Stage)) +
      geom_point(aes(x = input$cnv_init + 0.8, y = input$bcva_baseline),
                 color = "black", size = 8, shape = 4, stroke = 3) +
      annotate("text", x = input$cnv_init + 1.5, y = input$bcva_baseline,
               label = "Patient", fontface = "bold", size = 4) +
      labs(x = "Drusen Burden (proxy, mm²)", y = "BCVA (letters)",
           title = "Disease Staging Map") +
      theme_bw(base_size = 11) + theme(legend.position = "none") +
      xlim(0, 22) + ylim(25, 90)
  })

  output$tbl_baseline <- renderTable({
    di <- drug_info()
    data.frame(
      Parameter = c("Drug", "Dose", "Molecular Weight", "VEGF Kd", "Half-life (vitreous)",
                    "AMD Type", "Baseline BCVA", "Initial CNV", "Initial CST"),
      Value = c(
        input$drug_select,
        paste0(di$dp$dose, " mg IVT"),
        paste0(di$dp$MW, " kDa"),
        paste0(di$dp$Kd, " nM"),
        paste0(round(0.693 / di$dp$kel, 1), " days"),
        ifelse(as.numeric(input$amd_type)==1,"Wet AMD","Dry AMD"),
        paste0(input$bcva_baseline, " ETDRS letters"),
        paste0(input$cnv_init, " mm²"),
        paste0(280 + input$fluid_init, " μm")
      )
    )
  })

  ## ---- TAB 2: PK ----
  output$plot_pk_vit <- renderPlot({
    sim <- sim_main()
    ggplot(sim, aes(x = time, y = DRUG_VIT)) +
      geom_line(color = "#2980B9", size = 1.2) +
      geom_vline(xintercept = seq(0, max(sim$time), by = input$load_interval)[1:input$n_loading],
                 linetype = "dashed", color = "gray") +
      labs(x = "Time (days)", y = "Drug [nM]",
           title = paste("Vitreous PK —", input$drug_select)) +
      theme_bw(base_size = 11) + scale_y_log10()
  })

  output$plot_pk_ret <- renderPlot({
    sim <- sim_main()
    ggplot(sim, aes(x = time, y = DRUG_RET)) +
      geom_line(color = "#16A085", size = 1.2) +
      labs(x = "Time (days)", y = "Drug [nM]",
           title = "Retinal/RPE Drug Concentration") +
      theme_bw(base_size = 11) + scale_y_log10()
  })

  output$plot_vegf_bound <- renderPlot({
    sim <- sim_main()
    df <- sim %>% select(time, VEGF_FREE, VEGF_BOUND) %>%
      pivot_longer(-time, names_to = "form", values_to = "conc")
    ggplot(df, aes(x = time, y = conc, color = form)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = c(VEGF_FREE="#E74C3C", VEGF_BOUND="#95A5A6"),
                         labels = c("Free VEGF", "Bound VEGF (Drug:VEGF)")) +
      labs(x = "Time (days)", y = "VEGF (nM)", color = "",
           title = "Free vs Drug-Bound VEGF") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$tbl_pk_summary <- renderTable({
    sim <- sim_main()
    times <- c(7, 28, 56, 90, 180, 365)
    times <- times[times <= max(sim$time)]
    sim_sub <- sim[sim$time %in% times, ]
    data.frame(
      Day       = sim_sub$time,
      `Vitreous [nM]`  = round(sim_sub$DRUG_VIT, 2),
      `Retina [nM]`    = round(sim_sub$DRUG_RET, 3),
      `Free VEGF [nM]` = round(sim_sub$VEGF_FREE, 3),
      `VEGFR2 act`     = round(sim_sub$VEGFR2_ACT, 3),
      check.names = FALSE
    )
  })

  ## ---- TAB 3: PD ----
  output$plot_vegfr2 <- renderPlot({
    sim <- sim_main()
    noTx <- sim_noTx()
    df <- bind_rows(
      mutate(sim,  scenario = input$drug_select),
      mutate(noTx, scenario = "No Treatment")
    )
    ggplot(df, aes(x = time, y = VEGFR2_ACT, color = scenario)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c("No Treatment" = "#E74C3C",
                                    setNames("#2980B9", input$drug_select))) +
      labs(x = "Time (days)", y = "VEGFR-2 Activity (0–1)",
           title = "VEGFR-2 Activation Suppression", color = "") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$plot_vegf_free <- renderPlot({
    sim <- sim_main()
    noTx <- sim_noTx()
    df <- bind_rows(
      mutate(sim,  scenario = input$drug_select),
      mutate(noTx, scenario = "No Treatment")
    )
    ggplot(df, aes(x = time, y = VEGF_FREE, color = scenario)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c("No Treatment"="#E74C3C",
                                    setNames("#2980B9", input$drug_select))) +
      labs(x = "Time (days)", y = "Free VEGF (nM)",
           title = "Free VEGF Suppression", color = "") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$plot_ang2 <- renderPlot({
    sim <- sim_main()
    ggplot(sim, aes(x = time, y = ANG2_FREE)) +
      geom_line(color = "#8E44AD", size = 1.2) +
      geom_ribbon(aes(ymin = 0, ymax = ANG2_FREE), fill = "#8E44AD", alpha = 0.15) +
      labs(x = "Time (days)", y = "Free Ang-2 (nM)",
           title = "Ang-2 Free Levels\n(suppressed only by Faricimab)") +
      theme_bw(base_size = 11)
  })

  output$plot_pd_bcva_corr <- renderPlot({
    sim <- sim_main()
    ggplot(sim %>% filter(time > 28), aes(x = 1 - VEGFR2_ACT, y = BCVA_change)) +
      geom_point(aes(color = time), size = 1.5, alpha = 0.7) +
      geom_smooth(method = "loess", color = "black", se = FALSE) +
      scale_color_viridis_c(name = "Day") +
      labs(x = "VEGFR-2 Suppression (1 - activity)",
           y = "ΔBCVA from baseline (letters)",
           title = "PD-Response Relationship:\nVEGFR-2 Suppression → BCVA Gain") +
      theme_bw(base_size = 11)
  })

  ## ---- TAB 4: Endpoints ----
  output$plot_bcva <- renderPlot({
    sim <- sim_main()
    noTx <- sim_noTx()
    df <- bind_rows(
      mutate(sim,  scenario = input$drug_select),
      mutate(noTx, scenario = "No Treatment")
    )
    ggplot(df, aes(x = time, y = BCVA_SCORE, color = scenario)) +
      geom_line(size = 1.3) +
      geom_hline(yintercept = input$bcva_baseline + 15, linetype = "dashed",
                 color = "darkgreen", alpha = 0.7) +
      annotate("text", x = max(df$time) * 0.8, y = input$bcva_baseline + 16,
               label = "+15 letters\n(meaningful gain)", size = 3, color = "darkgreen") +
      scale_color_manual(values = c("No Treatment"="#E74C3C",
                                    setNames("#2980B9", input$drug_select))) +
      labs(x = "Time (days)", y = "BCVA (ETDRS letters)", color = "") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$plot_cst <- renderPlot({
    sim <- sim_main()
    noTx <- sim_noTx()
    df <- bind_rows(
      mutate(sim,  scenario = input$drug_select),
      mutate(noTx, scenario = "No Treatment")
    )
    ggplot(df, aes(x = time, y = CST, color = scenario)) +
      geom_line(size = 1.2) +
      geom_hline(yintercept = 280, linetype = "dashed") +
      annotate("text", x = 50, y = 270, label = "Normal CST ~280μm", size = 3) +
      scale_color_manual(values = c("No Treatment"="#E74C3C",
                                    setNames("#2980B9", input$drug_select))) +
      labs(x = "Time (days)", y = "CST (μm)", color = "") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$plot_cnv <- renderPlot({
    sim <- sim_main()
    noTx <- sim_noTx()
    df <- bind_rows(
      mutate(sim,  scenario = input$drug_select),
      mutate(noTx, scenario = "No Treatment")
    )
    ggplot(df, aes(x = time, y = CNV_AREA, color = scenario, fill = scenario)) +
      geom_line(size = 1.2) +
      geom_ribbon(aes(ymin = 0, ymax = CNV_AREA), alpha = 0.1) +
      scale_color_manual(values = c("No Treatment"="#E74C3C",
                                    setNames("#2980B9", input$drug_select))) +
      scale_fill_manual(values = c("No Treatment"="#E74C3C",
                                   setNames("#2980B9", input$drug_select))) +
      labs(x = "Time (days)", y = "CNV Area (mm²)", color = "", fill = "") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$plot_ga <- renderPlot({
    sim <- sim_main()
    ggplot(sim, aes(x = time, y = GA_AREA)) +
      geom_area(fill = "#E74C3C", alpha = 0.3) +
      geom_line(color = "#C0392B", size = 1.2) +
      labs(x = "Time (days)", y = "GA Area (mm²)",
           title = "Geographic Atrophy Progression") +
      theme_bw(base_size = 11)
  })

  ## ---- TAB 5: Scenarios ----
  output$plot_scenario_bcva <- renderPlot({
    sim_all() %>%
      ggplot(aes(x = time, y = BCVA_SCORE, color = scenario)) +
      geom_line(size = 1.2) +
      geom_hline(yintercept = input$bcva_baseline + 15, linetype = "dashed",
                 color = "darkgreen", alpha = 0.6) +
      scale_color_manual(values = c(
        "Ranibizumab" = "#2980B9", "Aflibercept" = "#E74C3C",
        "Bevacizumab" = "#27AE60", "Faricimab"    = "#8E44AD",
        "Brolucizumab"= "#F39C12", "No Treatment" = "#95A5A6"
      )) +
      labs(x = "Time (days)", y = "BCVA (ETDRS letters)",
           title = "BCVA Outcomes — All Anti-VEGF Regimens", color = "Regimen") +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
  })

  output$plot_scenario_cnv <- renderPlot({
    sim_all() %>%
      ggplot(aes(x = time, y = CNV_AREA, color = scenario)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = c(
        "Ranibizumab" = "#2980B9", "Aflibercept" = "#E74C3C",
        "Bevacizumab" = "#27AE60", "Faricimab"    = "#8E44AD",
        "Brolucizumab"= "#F39C12", "No Treatment" = "#95A5A6"
      )) +
      labs(x = "Time (days)", y = "CNV Area (mm²)",
           title = "CNV Lesion Area by Regimen", color = "Regimen") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$plot_scenario_burden <- renderPlot({
    burden_data <- data.frame(
      drug    = c("Ranibizumab q8w","Aflibercept q8w","Bevacizumab PRN",
                  "Faricimab q16w","Brolucizumab q12w","No Treatment"),
      n_inj   = c(8, 8, 7, 6, 6, 0),
      bcva_gain = c(9.2, 8.5, 8.7, 11.0, 9.8, -14.9)
    )
    ggplot(burden_data, aes(x = n_inj, y = bcva_gain, label = drug)) +
      geom_point(size = 5, aes(color = drug)) +
      ggrepel::geom_text_repel(size = 3.5, max.overlaps = 10) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(x = "Injections in Year 1", y = "BCVA gain at 1 year (letters)",
           title = "Treatment Burden vs. Efficacy (Clinical Trial Summary)") +
      theme_bw(base_size = 11) + theme(legend.position = "none")
  })

  output$tbl_scenarios <- renderTable({
    sim_res <- sim_all()
    sim_res %>%
      filter(time == 365) %>%
      select(scenario, BCVA_SCORE, BCVA_change, CNV_AREA, CST, GA_AREA) %>%
      mutate(across(where(is.numeric), ~round(.x, 2))) %>%
      rename(
        Treatment      = scenario,
        `BCVA (letters)` = BCVA_SCORE,
        `ΔBCVA`         = BCVA_change,
        `CNV (mm²)`     = CNV_AREA,
        `CST (μm)`      = CST,
        `GA (mm²)`      = GA_AREA
      )
  })

  ## ---- TAB 6: Biomarkers ----
  output$plot_complement <- renderPlot({
    sim <- sim_main()
    df <- sim %>% select(time, C3_LOCAL, C5_LOCAL, MAC_LOCAL) %>%
      pivot_longer(-time, names_to = "comp", values_to = "level")
    ggplot(df, aes(x = time, y = level, color = comp)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = c(C3_LOCAL="#27AE60",C5_LOCAL="#E67E22",MAC_LOCAL="#E74C3C"),
                         labels = c("C3 (local)","C5 (local)","MAC (C5b-9)")) +
      labs(x = "Time (days)", y = "Level (AU)", color = "") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$plot_rpe <- renderPlot({
    sim <- sim_main()
    df <- sim %>% select(time, RPE_NORM, LIPOFUSCIN) %>%
      pivot_longer(-time, names_to = "marker", values_to = "value")
    ggplot(df, aes(x = time, y = value, color = marker)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c(RPE_NORM="#27AE60", LIPOFUSCIN="#F39C12"),
                         labels = c("Lipofuscin (AU)","RPE cell fraction")) +
      labs(x = "Time (days)", y = "Level", color = "") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })

  output$plot_pr <- renderPlot({
    sim <- sim_main()
    ggplot(sim, aes(x = time, y = PR_FRAC * 100)) +
      geom_line(color = "#16A085", size = 1.2) +
      geom_ribbon(aes(ymin = 0, ymax = PR_FRAC * 100), fill = "#16A085", alpha = 0.15) +
      labs(x = "Time (days)", y = "Photoreceptor Survival (%)",
           title = "Photoreceptor Viability") +
      theme_bw(base_size = 11) + ylim(0, 100)
  })

  output$plot_drusen <- renderPlot({
    sim <- sim_main()
    noTx <- sim_noTx()
    df <- bind_rows(
      mutate(sim,  scenario = input$drug_select),
      mutate(noTx, scenario = "No Treatment")
    )
    ggplot(df, aes(x = time, y = DRUSEN, color = scenario)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = c("No Treatment"="#E74C3C",
                                    setNames("#2980B9", input$drug_select))) +
      labs(x = "Time (days)", y = "Drusen Area (mm²)",
           title = "Drusen Burden Progression", color = "") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
  })
}

shinyApp(ui, server)
