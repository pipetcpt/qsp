## =============================================================================
## Acute Intermittent Porphyria (AIP) – Interactive QSP Shiny Dashboard
## =============================================================================
## Author  : Claude Code Routine (CCR)
## Date    : 2026-06-25
## Usage   : shiny::runApp("aip_shiny_app.R")
##
## Tabs:
##   1. Patient Profile & Disease Parameters
##   2. Drug PK  (Givosiran plasma/liver + Hemin plasma)
##   3. PD Markers (ALAS1 mRNA/protein, ALA/PBG fold-change)
##   4. Clinical Endpoints (AAR, urine biomarkers, neurotoxicity)
##   5. Scenario Comparison (6 treatment arms, waterfall/spider)
##   6. Biomarkers (ALA/PBG normalization, renal function, VPop)
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ─── mrgsolve model code (embedded) ─────────────────────────────────────────
aip_code <- '
$PARAM
KA_GIV=0.50, F1_GIV=0.90, CL_GIV=0.038, V1_GIV=0.22,
Q_GIV=0.12, V2_GIV=1.60, BW=65.0,
KUP_LIV=4.0, GALN_KM=12.0, KCLEAR_LIV=0.0039, MW_GIV=14540,
EC50_siRNA=450.0, EMAX_siRNA=0.92, HILL_siRNA=1.50,
KDEG_mRNA=0.693, KDEG_PROT=0.231,
KSYN_ALA=5.0, KD_ALA=1.8, K12_ALA=0.50, K21_ALA=0.10, KREN_ALA=0.55,
ALA0_H=2.50, ALA0_P=0.50,
KSYN_PBG=1.60, PBGD_ACT=0.50, KD_PBG=1.20, K12_PBG=0.40, KREN_PBG=0.65,
PBG0_H=1.50, PBG0_P=0.25,
KSYN_HEME=3.50, KD_HEME=0.33, HEME0=8.50, KFB_HEME=2.20, EC50_HEME=8.50,
CL_HEM=1.60, V_HEM=0.09, KUP_HEM=6.0, KHO1_HEM=0.80,
KALAS1_HEM=3.50, EC50_HEM=4.50,
KNTOX_ACC=0.025, KREC_NTOX=0.06, ALA_THRESH=4.50,
CYCTRIG=0, TRIG_AMP=0.40, CYCP=28.0, CYCP_PH=21.0,
GLU_DOSE=0, GLU_IEFF=0.50, GLU_EC50=8.0

$CMT GIV_SC GIV_C GIV_P GIV_LIV ALAS1_mRNA ALAS1_PROT
     ALA_LIV ALA_PLAS PBG_LIV PBG_PLAS HEME_LIV
     HEM_C HEM_LIV NEUROTOX ATK_DAY AUC_ALA AUC_PBG

$MAIN
GIV_SC_0=0; GIV_C_0=0; GIV_P_0=0; GIV_LIV_0=0;
ALAS1_mRNA_0=1.0; ALAS1_PROT_0=1.0;
ALA_LIV_0=ALA0_H*(1.0/PBGD_ACT); ALA_PLAS_0=ALA0_P*(1.0/PBGD_ACT);
PBG_LIV_0=PBG0_H*(1.0/(PBGD_ACT*PBGD_ACT)); PBG_PLAS_0=PBG0_P*(1.0/PBGD_ACT);
HEME_LIV_0=HEME0*PBGD_ACT; HEM_C_0=0; HEM_LIV_0=0;
NEUROTOX_0=0; ATK_DAY_0=0; AUC_ALA_0=0; AUC_PBG_0=0;

$ODE
double abs_GIV=KA_GIV*GIV_SC;
dxdt_GIV_SC=-abs_GIV;
double GIV_C_ng=GIV_C*MW_GIV/1e6;
double UPT_LIV=KUP_LIV*GIV_C_ng/(GALN_KM+GIV_C_ng)*GIV_C;
dxdt_GIV_C=F1_GIV*abs_GIV/BW*1000-CL_GIV*GIV_C-Q_GIV*GIV_C+Q_GIV*GIV_P-KUP_LIV*GIV_C;
dxdt_GIV_P=Q_GIV*GIV_C-Q_GIV*GIV_P;
dxdt_GIV_LIV=UPT_LIV-KCLEAR_LIV*GIV_LIV;
double siEFF=EMAX_siRNA*pow(GIV_LIV,HILL_siRNA)/(pow(EC50_siRNA,HILL_siRNA)+pow(GIV_LIV,HILL_siRNA));
double t_cyc=fmod(SOLVERTIME,CYCP);
double HORMTRIG=1.0+CYCTRIG*TRIG_AMP*fmax(sin(2.0*3.14159265*(t_cyc-CYCP_PH)/CYCP),0.0);
double GLU_EFF=GLU_IEFF*GLU_DOSE/(GLU_EC50+GLU_DOSE);
double HEME_FB=1.0/(1.0+pow(HEME_LIV/EC50_HEME,KFB_HEME));
double HEM_EFF=1.0/(1.0+(HEM_LIV/EC50_HEM)*KALAS1_HEM);
double ALAS1_SYN=KDEG_mRNA*HORMTRIG*HEME_FB*HEM_EFF*(1.0-siEFF)*(1.0-GLU_EFF);
dxdt_ALAS1_mRNA=ALAS1_SYN-KDEG_mRNA*ALAS1_mRNA;
dxdt_ALAS1_PROT=KDEG_PROT*ALAS1_mRNA-KDEG_PROT*ALAS1_PROT;
double ALA_SYNTH=KSYN_ALA*ALAS1_PROT;
double ALA_UTIL=KD_ALA*ALA_LIV;
double ALA_EXP=K12_ALA*ALA_LIV;
double ALA_RET=K21_ALA*ALA_PLAS;
dxdt_ALA_LIV=ALA_SYNTH-ALA_UTIL-ALA_EXP+ALA_RET;
dxdt_ALA_PLAS=ALA_EXP-ALA_RET-KREN_ALA*ALA_PLAS;
double PBG_SYNTH=KSYN_PBG*ALA_LIV;
double PBG_UTIL=KD_PBG*PBGD_ACT*PBG_LIV;
double PBG_EXP=K12_PBG*PBG_LIV;
dxdt_PBG_LIV=PBG_SYNTH-PBG_UTIL-PBG_EXP;
dxdt_PBG_PLAS=PBG_EXP-KREN_PBG*PBG_PLAS;
double HEME_PROD=KSYN_HEME*PBGD_ACT*ALAS1_PROT;
dxdt_HEME_LIV=HEME_PROD-KD_HEME*HEME_LIV+KHO1_HEM*HEM_LIV*0.3;
dxdt_HEM_C=-(CL_HEM+KUP_HEM*V_HEM)*HEM_C;
dxdt_HEM_LIV=KUP_HEM*HEM_C*V_HEM/(BW*0.025)-KHO1_HEM*HEM_LIV;
double ALA_FC=ALA_PLAS/ALA0_P;
double NTOX_INPUT=fmax(0.0,ALA_FC-ALA_THRESH);
dxdt_NEUROTOX=KNTOX_ACC*NTOX_INPUT-KREC_NTOX*NEUROTOX;
double ATK_PROB=(ALA_PLAS>ALA_THRESH*ALA0_P)?1.0:0.0;
dxdt_ATK_DAY=ATK_PROB;
dxdt_AUC_ALA=ALA_PLAS;
dxdt_AUC_PBG=PBG_PLAS;

$TABLE
capture ALA_FC=ALA_PLAS/ALA0_P;
capture PBG_FC=PBG_PLAS/PBG0_P;
capture siEFF_out=EMAX_siRNA*pow(GIV_LIV,HILL_siRNA)/(pow(EC50_siRNA,HILL_siRNA)+pow(GIV_LIV,HILL_siRNA));
capture ALA_UR=KREN_ALA*ALA_PLAS*1440;
capture PBG_UR=KREN_PBG*PBG_PLAS*1440;
capture HORM_TRIG=1.0+CYCTRIG*TRIG_AMP*fmax(sin(2.0*3.14159265*(fmod(SOLVERTIME,CYCP)-CYCP_PH)/CYCP),0.0);
'

mod_global <- mcode("aip_shiny", aip_code, quiet = TRUE)

## ─── Helpers ─────────────────────────────────────────────────────────────────
run_sim <- function(mod, params_list, ev_obj = NULL, end = 365, delta = 0.5) {
  m <- do.call(param, c(list(mod), params_list))
  if (!is.null(ev_obj)) m <- ev(m, ev_obj)
  mrgsim(m, end = end, delta = delta,
         outvars = c("ALA_FC","PBG_FC","ALAS1_mRNA","ALAS1_PROT",
                     "ALA_PLAS","PBG_PLAS","HEME_LIV",
                     "GIV_C","GIV_LIV","HEM_C","HEM_LIV",
                     "siEFF_out","NEUROTOX","ATK_DAY",
                     "ALA_UR","PBG_UR","HORM_TRIG")) %>%
    as_tibble()
}

make_giv_ev <- function(bw, dose_mgkg, n_months) {
  dose_mg <- bw * dose_mgkg
  days    <- seq(0, (n_months - 1) * 28, by = 28)
  ev(amt = dose_mg, cmt = "GIV_SC", time = days)
}

make_hem_ev <- function(bw, start_day, dose_mgkg, n_days) {
  dose_ug <- bw * dose_mgkg * 1000
  ev(amt = dose_ug, cmt = "HEM_C", rate = -2,
     time = seq(start_day, start_day + n_days - 1, by = 1))
}

## ─── UI ──────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "AIP QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "tab_patient",   icon = icon("user")),
      menuItem("Drug PK",              tabName = "tab_pk",        icon = icon("pills")),
      menuItem("PD Markers",           tabName = "tab_pd",        icon = icon("chart-line")),
      menuItem("Clinical Endpoints",   tabName = "tab_clinical",  icon = icon("hospital")),
      menuItem("Scenario Comparison",  tabName = "tab_scenario",  icon = icon("balance-scale")),
      menuItem("Biomarkers & VPop",    tabName = "tab_biomarker", icon = icon("dna"))
    )
  ),

  dashboardBody(
    ## ── Shared CSS ──────────────────────────────────────────────────────────
    tags$head(tags$style(HTML(
      ".box-header{background-color:#5B2C8D!important;color:#fff!important;}
       .small-box{border-radius:8px;}
       .info-box-icon{background-color:#5B2C8D!important;}
       h4{color:#4A235A;font-weight:600;}
       .skin-purple .main-header .logo{background-color:#5B2C8D;}
       .skin-purple .main-header .navbar{background-color:#5B2C8D;}"
    ))),

    tabItems(

      ## ── Tab 1: Patient Profile ─────────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Demographics", width = 4, status = "primary",
            numericInput("bw",      "Body Weight (kg):",      value = 65, min = 40, max = 120),
            selectInput("sex",      "Sex:",                   choices = c("Female" = 1, "Male" = 0)),
            numericInput("egfr",    "eGFR (mL/min/1.73m²):", value = 90, min = 10, max = 130),
            sliderInput("pbgd_act", "PBGD Activity (fraction of normal):",
                        min = 0.1, max = 1.0, value = 0.50, step = 0.05),
            checkboxInput("cyctrig", "Enable Menstrual Cycle Trigger", value = TRUE),
            numericInput("sim_days", "Simulation Duration (days):", value = 365, min = 30, max = 730)
          ),
          box(title = "AIP Genetic Variant", width = 4, status = "warning",
            selectInput("hmbs_class", "HMBS Mutation Class:",
                        choices = c("Missense (common, ~50% activity)" = 0.50,
                                    "Nonsense/frameshift (~25% activity)" = 0.25,
                                    "Splice site (~40% activity)" = 0.40,
                                    "Exon deletion (~15% activity)" = 0.15,
                                    "Normal (wildtype control)" = 1.00)),
            actionButton("apply_genotype", "Apply Genotype", class = "btn-warning"),
            hr(),
            h4("Genotype → Phenotype"),
            p("PBGD activity directly sets the rate-limiting HMBS enzyme level,
               determining baseline ALA/PBG accumulation and attack susceptibility.")
          ),
          box(title = "Trigger Factors", width = 4, status = "danger",
            sliderInput("trig_amp", "Hormonal ALAS1 Upregulation (amplitude):",
                        min = 0, max = 1.0, value = 0.40, step = 0.05),
            sliderInput("glu_dose", "IV Glucose Rate (g/h; 0 = no glucose Rx):",
                        min = 0, max = 25, value = 0, step = 0.5),
            checkboxGroupInput("drug_triggers", "Active Drug Triggers (CYP inducers):",
                               choices = c("Barbiturates (+30% ALAS1)" = "barb",
                                           "Rifampicin (+25% ALAS1)"  = "rif",
                                           "Antiepileptics (+20% ALAS1)" = "aed"),
                               selected = NULL)
          )
        ),
        fluidRow(
          valueBoxOutput("vbox_ala_ss"),
          valueBoxOutput("vbox_pbg_ss"),
          valueBoxOutput("vbox_heme_ss")
        ),
        fluidRow(
          box(title = "Steady-State Disease Profile (Pre-Treatment)", width = 12,
            plotOutput("plot_ss_profile", height = "300px"))
        )
      ),

      ## ── Tab 2: Drug PK ─────────────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Givosiran PK Settings", width = 4, status = "primary",
            selectInput("giv_regimen", "Givosiran Regimen:",
                        choices = c("2.5 mg/kg Q28d (standard)" = 2.5,
                                    "1.25 mg/kg Q28d (reduced)" = 1.25,
                                    "5.0 mg/kg Q28d (high dose)"= 5.0,
                                    "No givosiran (control)"    = 0)),
            numericInput("giv_start_day", "First Dose Day:", value = 0, min = 0, max = 60),
            numericInput("giv_n_months",  "Number of Doses (months):", value = 12, min = 1, max = 24),
            hr(),
            h4("Givosiran PK Parameters"),
            sliderInput("cl_giv", "CL (L/kg/d):", min = 0.01, max = 0.10, value = 0.038, step = 0.002),
            sliderInput("v1_giv", "V1 (L/kg):",   min = 0.05, max = 0.50, value = 0.22,  step = 0.01),
            sliderInput("kup_liv","KUP_LIV (d⁻¹):", min = 1.0, max = 10.0, value = 4.0, step = 0.5)
          ),
          box(title = "Hemin IV Settings", width = 4, status = "danger",
            checkboxInput("use_hemin", "Add Hemin IV (acute attack treatment)", value = FALSE),
            numericInput("hem_start", "Hemin IV Start Day:", value = 30, min = 0, max = 360),
            sliderInput("hem_dose",   "Hemin Dose (mg/kg/d):", min = 1, max = 6, value = 3, step = 0.5),
            sliderInput("hem_days",   "Hemin Course Duration (days):", min = 2, max = 7, value = 4)
          ),
          box(title = "PK Summary", width = 4,
            tableOutput("table_pk_summary")
          )
        ),
        fluidRow(
          box(title = "Givosiran Plasma PK", width = 6,
            plotlyOutput("plot_giv_plasma", height = "300px")),
          box(title = "Givosiran Liver Concentration", width = 6,
            plotlyOutput("plot_giv_liver", height = "300px"))
        ),
        fluidRow(
          box(title = "Hemin IV Plasma Concentration", width = 6,
            plotlyOutput("plot_hem_plasma", height = "300px")),
          box(title = "siRNA Knockdown Efficiency", width = 6,
            plotlyOutput("plot_si_eff", height = "300px"))
        )
      ),

      ## ── Tab 3: PD Markers ──────────────────────────────────────────────────
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(title = "ALAS1 mRNA Knockdown", width = 6,
            plotlyOutput("plot_alas1_mrna", height = "300px")),
          box(title = "ALAS1 Protein Activity", width = 6,
            plotlyOutput("plot_alas1_prot", height = "300px"))
        ),
        fluidRow(
          box(title = "Plasma ALA (fold-change over normal)", width = 6,
            plotlyOutput("plot_ala_fc", height = "300px")),
          box(title = "Plasma PBG (fold-change over normal)", width = 6,
            plotlyOutput("plot_pbg_fc", height = "300px"))
        ),
        fluidRow(
          box(title = "Hepatic Free Heme Pool", width = 6,
            plotlyOutput("plot_heme_pool", height = "300px")),
          box(title = "Hormonal Trigger Pattern", width = 6,
            plotlyOutput("plot_horm_trig", height = "300px"))
        )
      ),

      ## ── Tab 4: Clinical Endpoints ──────────────────────────────────────────
      tabItem(tabName = "tab_clinical",
        fluidRow(
          valueBoxOutput("vbox_aar"),
          valueBoxOutput("vbox_ala_norm"),
          valueBoxOutput("vbox_pbg_norm")
        ),
        fluidRow(
          box(title = "Urinary ALA Rate (µmol/day)", width = 6,
            plotlyOutput("plot_urine_ala", height = "300px")),
          box(title = "Urinary PBG Rate (µmol/day)", width = 6,
            plotlyOutput("plot_urine_pbg", height = "300px"))
        ),
        fluidRow(
          box(title = "Neurotoxicity Index", width = 6,
            plotlyOutput("plot_neurotox", height = "300px")),
          box(title = "Cumulative Attack-Risk Days", width = 6,
            plotlyOutput("plot_atk_days", height = "300px"))
        ),
        fluidRow(
          box(title = "Clinical Endpoint Summary", width = 12,
            DTOutput("table_clinical_endpoints"))
        )
      ),

      ## ── Tab 5: Scenario Comparison ─────────────────────────────────────────
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "Comparison Scenarios", width = 3, status = "primary",
            checkboxGroupInput("scen_select", "Select Scenarios:",
              choices = c(
                "Placebo (No Rx)"              = "placebo",
                "Givosiran 1.25 mg/kg Q1M"     = "giv_low",
                "Givosiran 2.5 mg/kg Q1M (Std)"= "giv_std",
                "Givosiran 5.0 mg/kg Q1M (HD)" = "giv_hd",
                "Hemin IV (Acute)"             = "hemin",
                "Gene Therapy (PBGD 95%)"      = "gene"
              ),
              selected = c("placebo","giv_std","gene")
            ),
            selectInput("scen_endpoint", "Compare on:",
                        choices = c("Plasma ALA (FC)" = "ALA_FC",
                                    "Plasma PBG (FC)" = "PBG_FC",
                                    "ALAS1 mRNA"      = "ALAS1_mRNA",
                                    "Neurotoxicity"   = "NEUROTOX")),
            actionButton("run_scenarios", "Run Scenarios", class = "btn-primary btn-block")
          ),
          box(title = "Scenario Time-Course", width = 9,
            plotlyOutput("plot_scenario_tc", height = "400px"))
        ),
        fluidRow(
          box(title = "AAR Reduction Comparison", width = 6,
            plotlyOutput("plot_aar_bar", height = "350px")),
          box(title = "ALA & PBG Reduction Spider Plot", width = 6,
            plotlyOutput("plot_spider", height = "350px"))
        ),
        fluidRow(
          box(title = "Scenario Summary Table", width = 12,
            DTOutput("table_scenario_summary"))
        )
      ),

      ## ── Tab 6: Biomarkers & VPop ───────────────────────────────────────────
      tabItem(tabName = "tab_biomarker",
        fluidRow(
          box(title = "VPop Simulation Settings", width = 3, status = "info",
            numericInput("vpop_n",     "Number of Patients (N):", value = 50, min = 10, max = 200),
            sliderInput("vpop_cv_cl",  "PK CL IIV (CV%, log-normal):", min = 10, max = 60, value = 30),
            sliderInput("vpop_cv_kup", "ASGPR Uptake IIV (CV%):",      min = 10, max = 60, value = 35),
            sliderInput("vpop_female", "Female % (hormonal triggers):", min = 50, max = 100, value = 85),
            actionButton("run_vpop", "Run VPop (N patients)", class = "btn-info btn-block")
          ),
          box(title = "VPop: Plasma ALA Response (Median + 90% PI)", width = 9,
            plotlyOutput("plot_vpop_ala", height = "400px"))
        ),
        fluidRow(
          box(title = "ALA Normalization Rate Over Time", width = 6,
            plotlyOutput("plot_ala_norm_rate", height = "300px")),
          box(title = "PBG Normalization Rate Over Time", width = 6,
            plotlyOutput("plot_pbg_norm_rate", height = "300px"))
        ),
        fluidRow(
          box(title = "Exposure-Response: AUC vs AAR Reduction", width = 6,
            plotlyOutput("plot_er_auc", height = "300px")),
          box(title = "Renal Biomarker: Estimated GFR Trend", width = 6,
            plotlyOutput("plot_egfr_trend", height = "300px"))
        )
      )

    ) # end tabItems
  )   # end dashboardBody
)

## ─── Server ──────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## -- Apply genotype button --------------------------------------------------
  observeEvent(input$apply_genotype, {
    updateSliderInput(session, "pbgd_act",
                      value = as.numeric(input$hmbs_class))
  })

  ## -- Reactive: base params --------------------------------------------------
  base_params <- reactive({
    trig_mult <- 1.0
    if ("barb" %in% input$drug_triggers) trig_mult <- trig_mult + 0.30
    if ("rif"  %in% input$drug_triggers) trig_mult <- trig_mult + 0.25
    if ("aed"  %in% input$drug_triggers) trig_mult <- trig_mult + 0.20

    list(
      BW       = input$bw,
      PBGD_ACT = input$pbgd_act,
      CYCTRIG  = as.numeric(input$sex) * as.numeric(input$cyctrig),
      TRIG_AMP = input$trig_amp * trig_mult,
      CL_GIV   = input$cl_giv,
      V1_GIV   = input$v1_giv,
      KUP_LIV  = input$kup_liv,
      GLU_DOSE = input$glu_dose
    )
  })

  ## -- Reactive: event object -------------------------------------------------
  ev_object <- reactive({
    giv_dose <- as.numeric(input$giv_regimen)
    ev_out   <- NULL

    if (giv_dose > 0) {
      ev_out <- make_giv_ev(bw       = input$bw,
                            dose_mgkg= giv_dose,
                            n_months = input$giv_n_months)
    }
    if (input$use_hemin) {
      hem_ev <- make_hem_ev(bw        = input$bw,
                            start_day = input$hem_start,
                            dose_mgkg = input$hem_dose,
                            n_days    = input$hem_days)
      ev_out <- if (is.null(ev_out)) hem_ev else ev_out + hem_ev
    }
    ev_out
  })

  ## -- Reactive: simulation output --------------------------------------------
  sim_out <- reactive({
    run_sim(mod_global, base_params(), ev_object(), end = input$sim_days)
  })

  ## ── Tab 1: Value Boxes ─────────────────────────────────────────────────────
  output$vbox_ala_ss <- renderValueBox({
    d <- sim_out()
    ss_ala <- round(tail(d$ALA_FC, 20) %>% mean(), 2)
    color  <- if (ss_ala > 4.5) "red" else if (ss_ala > 2) "yellow" else "green"
    valueBox(paste0(ss_ala, "×"), "Plasma ALA (× normal at SS)", icon = icon("vial"), color = color)
  })
  output$vbox_pbg_ss <- renderValueBox({
    d <- sim_out()
    ss_pbg <- round(tail(d$PBG_FC, 20) %>% mean(), 2)
    color  <- if (ss_pbg > 4.5) "red" else if (ss_pbg > 2) "yellow" else "green"
    valueBox(paste0(ss_pbg, "×"), "Plasma PBG (× normal at SS)", icon = icon("flask"), color = color)
  })
  output$vbox_heme_ss <- renderValueBox({
    d   <- sim_out()
    ss_h <- round(tail(d$HEME_LIV, 20) %>% mean(), 2)
    valueBox(paste0(ss_h, " nmol/g"), "Hepatic Free Heme (SS)", icon = icon("heartbeat"), color = "purple")
  })

  output$plot_ss_profile <- renderPlot({
    d <- sim_out()
    d_long <- d %>%
      select(time, ALA_FC, PBG_FC, ALAS1_mRNA, siEFF_out) %>%
      pivot_longer(-time, names_to = "var", values_to = "val") %>%
      mutate(var = recode(var,
        ALA_FC    = "Plasma ALA (FC)",
        PBG_FC    = "Plasma PBG (FC)",
        ALAS1_mRNA= "ALAS1 mRNA (rel)",
        siEFF_out = "siRNA KD (fraction)"
      ))
    ggplot(d_long, aes(x = time, y = val, color = var)) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~var, scales = "free_y") +
      scale_color_brewer(palette = "Set1") +
      labs(x = "Time (days)", y = "Value", color = NULL) +
      theme_bw() + theme(legend.position = "none")
  })

  ## ── Tab 2: PK Plots ────────────────────────────────────────────────────────
  output$plot_giv_plasma <- renderPlotly({
    d <- sim_out()
    p <- ggplot(d, aes(x = time, y = GIV_C)) +
      geom_line(color = "#6A0DAD", linewidth = 0.8) +
      labs(x = "Time (days)", y = "Givosiran Plasma (µg/L)") + theme_bw()
    ggplotly(p)
  })
  output$plot_giv_liver <- renderPlotly({
    d <- sim_out()
    p <- ggplot(d, aes(x = time, y = GIV_LIV)) +
      geom_line(color = "#4C1D95", linewidth = 0.8) +
      labs(x = "Time (days)", y = "Liver Givosiran (ng/g)") + theme_bw()
    ggplotly(p)
  })
  output$plot_hem_plasma <- renderPlotly({
    d <- sim_out()
    p <- ggplot(d, aes(x = time, y = HEM_C)) +
      geom_line(color = "#DC143C", linewidth = 0.8) +
      labs(x = "Time (days)", y = "Hemin Plasma (µg/L)") + theme_bw()
    ggplotly(p)
  })
  output$plot_si_eff <- renderPlotly({
    d <- sim_out()
    p <- ggplot(d, aes(x = time, y = siEFF_out * 100)) +
      geom_line(color = "#6A0DAD", linewidth = 0.8) +
      geom_hline(yintercept = 80, linetype = "dashed", color = "grey50") +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = "Time (days)", y = "ALAS1 Knockdown (%)") + theme_bw()
    ggplotly(p)
  })
  output$table_pk_summary <- renderTable({
    d   <- sim_out()
    last <- tail(d, 1)
    tibble(
      Parameter = c("Cmax Givosiran Plasma (µg/L)",
                    "Liver Givosiran Conc. (ng/g)",
                    "siRNA KD Efficiency (%)",
                    "Hemin Plasma Peak (µg/L)",
                    "Hemin Liver Peak (nmol/g)"),
      Value = round(c(max(d$GIV_C, na.rm=T),
                      max(d$GIV_LIV, na.rm=T),
                      max(d$siEFF_out, na.rm=T)*100,
                      max(d$HEM_C, na.rm=T),
                      max(d$HEM_LIV, na.rm=T)), 3)
    )
  })

  ## ── Tab 3: PD Plots ────────────────────────────────────────────────────────
  output$plot_alas1_mrna <- renderPlotly({
    d <- sim_out()
    p <- ggplot(d, aes(x = time, y = ALAS1_mRNA)) +
      geom_line(color = "#E63946", linewidth = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      labs(x = "Time (days)", y = "ALAS1 mRNA (relative)") + theme_bw()
    ggplotly(p)
  })
  output$plot_alas1_prot <- renderPlotly({
    d <- sim_out()
    p <- ggplot(d, aes(x = time, y = ALAS1_PROT)) +
      geom_line(color = "#FF6600", linewidth = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      labs(x = "Time (days)", y = "ALAS1 Protein (relative)") + theme_bw()
    ggplotly(p)
  })
  output$plot_ala_fc <- renderPlotly({
    d <- sim_out()
    p <- ggplot(d, aes(x = time, y = ALA_FC)) +
      geom_line(color = "#DC143C", linewidth = 0.8) +
      geom_hline(yintercept = 4.5, linetype = "dashed", color = "red", alpha = 0.6) +
      geom_hline(yintercept = 1, linetype = "dotted", color = "grey50") +
      labs(x = "Time (days)", y = "ALA fold-change") + theme_bw()
    ggplotly(p)
  })
  output$plot_pbg_fc <- renderPlotly({
    d <- sim_out()
    p <- ggplot(d, aes(x = time, y = PBG_FC)) +
      geom_line(color = "#E9C46A", linewidth = 0.8) +
      geom_hline(yintercept = 1, linetype = "dotted", color = "grey50") +
      labs(x = "Time (days)", y = "PBG fold-change") + theme_bw()
    ggplotly(p)
  })
  output$plot_heme_pool <- renderPlotly({
    d <- sim_out()
    p <- ggplot(d, aes(x = time, y = HEME_LIV)) +
      geom_line(color = "#8B0000", linewidth = 0.8) +
      geom_hline(yintercept = 8.5, linetype = "dashed", color = "grey50") +
      labs(x = "Time (days)", y = "Hepatic Heme (nmol/g)") + theme_bw()
    ggplotly(p)
  })
  output$plot_horm_trig <- renderPlotly({
    d <- sim_out()
    p <- ggplot(d, aes(x = time, y = HORM_TRIG)) +
      geom_line(color = "#F4A261", linewidth = 0.8) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      labs(x = "Time (days)", y = "Hormonal Trigger (amplitude)") + theme_bw()
    ggplotly(p)
  })

  ## ── Tab 4: Clinical Endpoints ──────────────────────────────────────────────
  output$vbox_aar <- renderValueBox({
    d   <- sim_out()
    aar <- round(max(d$ATK_DAY, na.rm=T) / max(d$time, na.rm=T) * 365, 1)
    col <- if (aar > 5) "red" else if (aar > 2) "yellow" else "green"
    valueBox(aar, "Annualized Attack Rate (AAR proxy)", icon = icon("exclamation-triangle"), color = col)
  })
  output$vbox_ala_norm <- renderValueBox({
    d    <- sim_out()
    pct  <- round(mean(tail(d$ALA_FC, 50) < 2.0) * 100, 1)
    col  <- if (pct > 70) "green" else if (pct > 40) "yellow" else "red"
    valueBox(paste0(pct, "%"), "ALA Normalization Rate (SS)", icon = icon("check"), color = col)
  })
  output$vbox_pbg_norm <- renderValueBox({
    d   <- sim_out()
    pct <- round(mean(tail(d$PBG_FC, 50) < 2.0) * 100, 1)
    col <- if (pct > 60) "green" else if (pct > 30) "yellow" else "red"
    valueBox(paste0(pct, "%"), "PBG Normalization Rate (SS)", icon = icon("check-circle"), color = col)
  })

  output$plot_urine_ala <- renderPlotly({
    d <- sim_out()
    p <- ggplot(d, aes(x = time, y = ALA_UR)) +
      geom_line(color = "#E63946", linewidth = 0.8) +
      labs(x = "Time (days)", y = "Urine ALA (µmol/day)") + theme_bw()
    ggplotly(p)
  })
  output$plot_urine_pbg <- renderPlotly({
    d <- sim_out()
    p <- ggplot(d, aes(x = time, y = PBG_UR)) +
      geom_line(color = "#F4A261", linewidth = 0.8) +
      labs(x = "Time (days)", y = "Urine PBG (µmol/day)") + theme_bw()
    ggplotly(p)
  })
  output$plot_neurotox <- renderPlotly({
    d <- sim_out()
    p <- ggplot(d, aes(x = time, y = NEUROTOX)) +
      geom_line(color = "#8B0000", linewidth = 0.8) +
      labs(x = "Time (days)", y = "Neurotoxicity Index") + theme_bw()
    ggplotly(p)
  })
  output$plot_atk_days <- renderPlotly({
    d <- sim_out()
    p <- ggplot(d, aes(x = time, y = ATK_DAY)) +
      geom_line(color = "#CC0000", linewidth = 0.8) +
      labs(x = "Time (days)", y = "Cumulative Attack-Risk Days") + theme_bw()
    ggplotly(p)
  })
  output$table_clinical_endpoints <- renderDT({
    d <- sim_out()
    tab <- tibble(
      Endpoint = c("Annualized Attack Rate (proxy)",
                   "Mean Plasma ALA (fold-change)",
                   "Mean Plasma PBG (fold-change)",
                   "ALA <2× normal (normalization rate %)",
                   "PBG <2× normal (normalization rate %)",
                   "Mean ALAS1 mRNA Knockdown (%)",
                   "Peak Neurotoxicity Index",
                   "Cumulative AUC-ALA (µmol·d/L)",
                   "Cumulative AUC-PBG (µmol·d/L)"),
      Value = c(
        round(max(d$ATK_DAY,na.rm=T)/max(d$time,na.rm=T)*365, 2),
        round(mean(d$ALA_FC, na.rm=T), 3),
        round(mean(d$PBG_FC, na.rm=T), 3),
        round(mean(d$ALA_FC < 2.0, na.rm=T)*100, 1),
        round(mean(d$PBG_FC < 2.0, na.rm=T)*100, 1),
        round((1 - mean(d$ALAS1_mRNA, na.rm=T))*100, 1),
        round(max(d$NEUROTOX, na.rm=T), 4),
        round(max(d$AUC_ALA, na.rm=T), 1),
        round(max(d$AUC_PBG, na.rm=T), 1)
      )
    )
    datatable(tab, options = list(dom = 't', pageLength = 20),
              rownames = FALSE, class = "stripe hover")
  })

  ## ── Tab 5: Scenario Comparison ─────────────────────────────────────────────
  scen_data <- eventReactive(input$run_scenarios, {
    bp  <- base_params()
    bw  <- bp$BW
    end <- input$sim_days

    scen_list <- list()
    scens_sel <- input$scen_select

    if ("placebo" %in% scens_sel)
      scen_list[["Placebo"]] <-
        run_sim(mod_global, bp, NULL, end = end) %>% mutate(scenario="Placebo")

    if ("giv_low" %in% scens_sel)
      scen_list[["Giv 1.25"]] <-
        run_sim(mod_global, bp, make_giv_ev(bw, 1.25, 13), end = end) %>%
        mutate(scenario="Givosiran 1.25 mg/kg")

    if ("giv_std" %in% scens_sel)
      scen_list[["Giv 2.5"]] <-
        run_sim(mod_global, bp, make_giv_ev(bw, 2.5, 13), end = end) %>%
        mutate(scenario="Givosiran 2.5 mg/kg (Std)")

    if ("giv_hd" %in% scens_sel)
      scen_list[["Giv 5.0"]] <-
        run_sim(mod_global, bp, make_giv_ev(bw, 5.0, 13), end = end) %>%
        mutate(scenario="Givosiran 5.0 mg/kg (HD)")

    if ("hemin" %in% scens_sel)
      scen_list[["Hemin"]] <-
        run_sim(mod_global, bp, make_hem_ev(bw, 30, 3, 4), end = min(120, end)) %>%
        mutate(scenario="Hemin IV (Acute)")

    if ("gene" %in% scens_sel) {
      bp_gene <- bp; bp_gene$PBGD_ACT <- 0.95
      scen_list[["Gene"]] <-
        run_sim(mod_global, bp_gene, NULL, end = end) %>%
        mutate(scenario="Gene Therapy")
    }
    bind_rows(scen_list)
  })

  output$plot_scenario_tc <- renderPlotly({
    req(scen_data())
    d   <- scen_data()
    col_name <- input$scen_endpoint
    p   <- ggplot(d, aes_string(x = "time", y = col_name, color = "scenario")) +
      geom_line(linewidth = 0.8) +
      scale_color_brewer(palette = "Set1") +
      labs(x = "Time (days)", y = col_name, color = NULL) +
      theme_bw()
    ggplotly(p)
  })

  output$plot_aar_bar <- renderPlotly({
    req(scen_data())
    d <- scen_data()
    aar_df <- d %>%
      group_by(scenario) %>%
      summarise(AAR = max(ATK_DAY, na.rm=T) / max(time, na.rm=T) * 365, .groups="drop")
    p <- ggplot(aar_df, aes(x = reorder(scenario, -AAR), y = AAR, fill = scenario)) +
      geom_col(show.legend = FALSE) +
      scale_fill_brewer(palette = "Set1") +
      coord_flip() +
      labs(x = NULL, y = "Attack-Risk Days/Year") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_spider <- renderPlotly({
    req(scen_data())
    d <- scen_data()
    spdr <- d %>%
      group_by(scenario) %>%
      summarise(
        ALA_Mean   = mean(ALA_FC, na.rm=T),
        PBG_Mean   = mean(PBG_FC, na.rm=T),
        ALAS1_KD   = (1-mean(ALAS1_mRNA, na.rm=T))*100,
        Ntox_Peak  = max(NEUROTOX, na.rm=T)*10,
        .groups = "drop"
      )
    spdr_long <- spdr %>%
      pivot_longer(-scenario, names_to = "axis", values_to = "val")
    p <- ggplot(spdr_long, aes(x = axis, y = val, fill = scenario, group = scenario)) +
      geom_col(position = "dodge") +
      scale_fill_brewer(palette = "Set1") +
      labs(x = NULL, y = "Value", fill = NULL) +
      theme_bw() + theme(axis.text.x = element_text(angle = 25, hjust=1))
    ggplotly(p)
  })

  output$table_scenario_summary <- renderDT({
    req(scen_data())
    d <- scen_data()
    tab <- d %>%
      group_by(scenario) %>%
      summarise(
        AAR_proxy   = round(max(ATK_DAY,na.rm=T)/max(time,na.rm=T)*365, 2),
        ALA_FC_mean = round(mean(ALA_FC, na.rm=T), 3),
        PBG_FC_mean = round(mean(PBG_FC, na.rm=T), 3),
        ALAS1_KD_pct= round((1-mean(ALAS1_mRNA,na.rm=T))*100, 1),
        Ntox_peak   = round(max(NEUROTOX, na.rm=T), 4),
        .groups     = "drop"
      )
    datatable(tab, options = list(dom = 't', pageLength = 10),
              rownames = FALSE, class = "stripe hover") %>%
      formatStyle("AAR_proxy", background = styleColorBar(tab$AAR_proxy, "lightcoral"))
  })

  ## ── Tab 6: VPop ─────────────────────────────────────────────────────────────
  vpop_data <- eventReactive(input$run_vpop, {
    N   <- input$vpop_n
    bp  <- base_params()
    cv_cl  <- input$vpop_cv_cl / 100
    cv_kup <- input$vpop_cv_kup / 100
    pct_f  <- input$vpop_female / 100
    bw     <- bp$BW
    set.seed(123)

    sims <- lapply(1:N, function(i) {
      params_i <- bp
      params_i$CL_GIV  <- bp$CL_GIV * exp(rnorm(1, 0, cv_cl))
      params_i$KUP_LIV <- bp$KUP_LIV * exp(rnorm(1, 0, cv_kup))
      params_i$CYCTRIG <- as.numeric(runif(1) < pct_f) * as.numeric(input$cyctrig)
      bw_i <- rnorm(1, bw, bw * 0.15) |> pmax(40) |> pmin(110)
      params_i$BW <- bw_i
      ev_i <- make_giv_ev(bw_i, 2.5, 13)
      run_sim(mod_global, params_i, ev_i, end = input$sim_days) %>%
        mutate(ID = i)
    })
    bind_rows(sims)
  })

  output$plot_vpop_ala <- renderPlotly({
    req(vpop_data())
    d <- vpop_data() %>%
      group_by(time) %>%
      summarise(med = median(ALA_FC, na.rm=T),
                lo  = quantile(ALA_FC, 0.05, na.rm=T),
                hi  = quantile(ALA_FC, 0.95, na.rm=T), .groups="drop")
    p <- ggplot(d, aes(x = time)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#457B9D", alpha = 0.3) +
      geom_line(aes(y = med), color = "#457B9D", linewidth = 1) +
      geom_hline(yintercept = 4.5, linetype = "dashed", color = "red", alpha = 0.6) +
      geom_hline(yintercept = 1, linetype = "dotted", color = "grey50") +
      labs(x = "Time (days)", y = "Plasma ALA (fold-change)", color = NULL) +
      theme_bw()
    ggplotly(p)
  })

  output$plot_ala_norm_rate <- renderPlotly({
    req(vpop_data())
    d <- vpop_data() %>%
      group_by(time) %>%
      summarise(norm_rate = mean(ALA_FC < 2.0, na.rm=T)*100, .groups="drop")
    p <- ggplot(d, aes(x = time, y = norm_rate)) +
      geom_line(color = "#2A9D8F", linewidth = 0.8) +
      geom_hline(yintercept = 73, linetype = "dashed", color = "grey50") +
      annotate("text", x = 5, y = 75, label = "ENVISION: 73% at 6M", size = 3.5, color="grey40") +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = "Time (days)", y = "ALA Normalization Rate (%)") + theme_bw()
    ggplotly(p)
  })

  output$plot_pbg_norm_rate <- renderPlotly({
    req(vpop_data())
    d <- vpop_data() %>%
      group_by(time) %>%
      summarise(norm_rate = mean(PBG_FC < 2.0, na.rm=T)*100, .groups="drop")
    p <- ggplot(d, aes(x = time, y = norm_rate)) +
      geom_line(color = "#E9C46A", linewidth = 0.8) +
      geom_hline(yintercept = 63, linetype = "dashed", color = "grey50") +
      scale_y_continuous(limits = c(0, 100)) +
      labs(x = "Time (days)", y = "PBG Normalization Rate (%)") + theme_bw()
    ggplotly(p)
  })

  output$plot_er_auc <- renderPlotly({
    req(vpop_data())
    d <- vpop_data() %>%
      group_by(ID) %>%
      summarise(
        AUC_LIV  = max(GIV_LIV, na.rm=T),
        AAR_proxy = max(ATK_DAY,na.rm=T)/max(time,na.rm=T)*365,
        .groups = "drop"
      )
    p <- ggplot(d, aes(x = AUC_LIV, y = AAR_proxy)) +
      geom_point(color = "#6A0DAD", alpha = 0.5) +
      geom_smooth(method = "loess", color = "#4C1D95", se = TRUE, alpha = 0.2) +
      labs(x = "Liver Givosiran Cmax (ng/g)", y = "AAR Proxy") + theme_bw()
    ggplotly(p)
  })

  output$plot_egfr_trend <- renderPlotly({
    d <- sim_out()
    # Simulate eGFR decline due to chronic ALA nephrotoxicity
    egfr_data <- d %>%
      mutate(
        AUC_cum   = cumsum(ALA_PLAS) * 0.5,  # cumulative ALA × time step
        eGFR_decl = pmax(input$egfr - AUC_cum * 0.0005, 15)  # 0.05 mL/min per µmol·d
      )
    p <- ggplot(egfr_data, aes(x = time, y = eGFR_decl)) +
      geom_line(color = "#457B9D", linewidth = 0.8) +
      geom_hline(yintercept = 60, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 30, linetype = "dotted", color = "red") +
      labs(x = "Time (days)",
           y = "Estimated eGFR (mL/min/1.73m²)",
           caption = "Estimated ALA-mediated nephrotoxicity trajectory") +
      theme_bw()
    ggplotly(p)
  })

}

## ─── Launch App ──────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
