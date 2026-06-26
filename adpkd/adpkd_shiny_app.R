## ============================================================
##  ADPKD QSP Shiny Dashboard
##  Autosomal Dominant Polycystic Kidney Disease
##
##  Tabs:
##    1. Patient Profile       — Genotype, baseline TKV/eGFR, Mayo class
##    2. Drug PK               — Tolvaptan/Everolimus/Octreotide Cp curves
##    3. PD Biomarkers         — V2R occupancy, urine osmolality, cAMP, mTOR
##    4. Disease Progression   — TKV growth, eGFR decline over 3 years
##    5. Scenario Comparison   — 5 treatment arms, forest plot, table
##    6. Biomarker & Risk      — PROPKD score, Mayo class, ESRD prediction
## ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)
library(DT)
library(bslib)

## ============================================================
##  Embedded mrgsolve model
## ============================================================

adpkd_code <- '
$PARAM
KA_TOLV=1.0 CL_TOLV=4.0 V1_TOLV=10.0 Q_TOLV=2.0 V2_TOLV=20.0 F_TOLV=0.56
EC50_TOLV=50.0 EMAX_TOLV=1.0 HILL_TOLV=1.5
KA_EVER=0.4 CL_EVER=14.0 V1_EVER=20.0 Q_EVER=5.0 V2_EVER=50.0
EC50_EVER=5.0 EMAX_EVER=0.80 HILL_EVER=1.0
KREL_OCT=0.002 CL_OCT=0.5 V1_OCT=8.0
KA_ACEI=0.8 CL_ACEI=5.0 V1_ACEI=15.0 EC50_ACEI=0.5 EMAX_ACEI=0.7
KIN_CAMP=0.5 KOUT_CAMP=0.5 CAMP_SS=1.0 CA2_SCALE=1.5
KIN_MTOR=0.3 KOUT_MTOR=0.3 MTOR_SS=1.0 CAMP_MTOR=0.3
KIN_ANGII=0.2 KOUT_ANGII=0.2 ANGII_SS=1.0 RENIN_TKV=0.3
KIN_BP=26.0 KOUT_BP=0.20 ANGII_BP=5.0 BP_SS=130.0
TKV0=1500.0 KGROW_TKV=6.28e-6 EGFR0=70.0 KDECL_EGFR=4.00e-4
TKV_EGFR=0.5 BP_EGFR=0.3
FCAMP=0.50 FMTOR=0.30 FBASE=0.20
UOSM_BASE=600.0 UOSM_MIN=50.0 KOUT_UOSM=0.5

$CMT AGUT ACENT APERI EGUT ECENT EPERI OCTDEP OCTCENT ACEI_GUT ACEI_CENT
AVP_ST CAMP_ST MTOR_ST ANGII_ST BP_ST TKV_ST EGFR_ST UOSM_ST NEPH_ST

$MAIN
F_AGUT = F_TOLV;
AVP_ST_0=1.0; CAMP_ST_0=CA2_SCALE; MTOR_ST_0=1.2;
ANGII_ST_0=1.0; BP_ST_0=BP_SS; TKV_ST_0=TKV0;
EGFR_ST_0=EGFR0; UOSM_ST_0=UOSM_BASE; NEPH_ST_0=1.0;

$ODE
double Cp_tolv_ngml = (ACENT/V1_TOLV)*1000.0;
double Ce_ever_ngml = (ECENT/V1_EVER)*1000.0;
double Cp_oct_ngml  = (OCTCENT/V1_OCT)*1000.0;
double Cp_acei_mgl  = ACEI_CENT/V1_ACEI;
double INH_TOLV = EMAX_TOLV*pow(Cp_tolv_ngml,HILL_TOLV)/(pow(EC50_TOLV,HILL_TOLV)+pow(Cp_tolv_ngml,HILL_TOLV));
double INH_EVER = EMAX_EVER*Ce_ever_ngml/(EC50_EVER+Ce_ever_ngml);
double INH_OCT  = 0.40*Cp_oct_ngml/(1.0+Cp_oct_ngml);
double INH_ACEI = EMAX_ACEI*Cp_acei_mgl/(EC50_ACEI+Cp_acei_mgl);
dxdt_AGUT   = -KA_TOLV*AGUT;
dxdt_ACENT  =  KA_TOLV*AGUT-(CL_TOLV/V1_TOLV)*ACENT-(Q_TOLV/V1_TOLV)*ACENT+(Q_TOLV/V2_TOLV)*APERI;
dxdt_APERI  = (Q_TOLV/V1_TOLV)*ACENT-(Q_TOLV/V2_TOLV)*APERI;
dxdt_EGUT   = -KA_EVER*EGUT;
dxdt_ECENT  =  KA_EVER*EGUT-(CL_EVER/V1_EVER)*ECENT-(Q_EVER/V1_EVER)*ECENT+(Q_EVER/V2_EVER)*EPERI;
dxdt_EPERI  = (Q_EVER/V1_EVER)*ECENT-(Q_EVER/V2_EVER)*EPERI;
dxdt_OCTDEP  = -KREL_OCT*OCTDEP;
dxdt_OCTCENT =  KREL_OCT*OCTDEP*1000.0-(CL_OCT/V1_OCT)*OCTCENT;
dxdt_ACEI_GUT  = -KA_ACEI*ACEI_GUT;
dxdt_ACEI_CENT =  KA_ACEI*ACEI_GUT-(CL_ACEI/V1_ACEI)*ACEI_CENT;
dxdt_AVP_ST  = KIN_CAMP - KOUT_CAMP*AVP_ST;
double CAMP_IN  = KIN_CAMP*AVP_ST*CA2_SCALE*(1.0-INH_TOLV)*(1.0-INH_OCT);
dxdt_CAMP_ST = CAMP_IN - KOUT_CAMP*CAMP_ST;
double MTOR_IN  = KIN_MTOR*(1.0+CAMP_MTOR*CAMP_ST);
dxdt_MTOR_ST = MTOR_IN - KOUT_MTOR*MTOR_ST*(1.0+INH_EVER);
double RENIN_F  = 1.0+RENIN_TKV*(TKV_ST/TKV0-1.0);
dxdt_ANGII_ST = KIN_ANGII*RENIN_F*(1.0-INH_ACEI) - KOUT_ANGII*ANGII_ST;
dxdt_BP_ST = KIN_BP+ANGII_BP*ANGII_ST - KOUT_BP*BP_ST;
double TKV_MOD = FBASE + FCAMP*(CAMP_ST/CAMP_SS) + FMTOR*(MTOR_ST/MTOR_SS);
double TKV_MOD_NORM = FBASE + FCAMP*CA2_SCALE + FMTOR*1.2;
dxdt_TKV_ST = (KGROW_TKV*TKV_MOD/TKV_MOD_NORM)*TKV_ST;
double TC = 1.0+TKV_EGFR*(TKV_ST/TKV0-1.0);
double BD = 1.0+BP_EGFR*(BP_ST/BP_SS-1.0);
if(TC<0.5) TC=0.5; if(BD<0.5) BD=0.5;
dxdt_EGFR_ST = (EGFR_ST>5.0) ? -KDECL_EGFR*TC*BD : 0.0;
dxdt_NEPH_ST = -KDECL_EGFR*TC*BD/EGFR0;
double UOSM_T = UOSM_BASE*(1.0-0.90*INH_TOLV)+UOSM_MIN*0.90*INH_TOLV;
if(UOSM_T<UOSM_MIN) UOSM_T=UOSM_MIN;
dxdt_UOSM_ST = KOUT_UOSM*(UOSM_T-UOSM_ST);

$TABLE
double Cp_tolv=Cp_tolv_ngml; double Ce_ever=Ce_ever_ngml;
double V2R_OCC=INH_TOLV*100.0; double mTOR_INH=INH_EVER*100.0;
double TKV_L=TKV_ST/1000.0;
double TKV_pct=(TKV_ST/TKV0-1.0)*100.0;
double eGFR=EGFR_ST; double BP=BP_ST; double Uosm=UOSM_ST;
double CAMP_rel=CAMP_ST/CAMP_SS; double MTOR_rel=MTOR_ST/MTOR_SS;

$CAPTURE Cp_tolv Ce_ever V2R_OCC mTOR_INH TKV_L TKV_pct eGFR BP Uosm CAMP_rel MTOR_rel
'

## Compile once
mod_app <- tryCatch(mcode("adpkd_shiny", adpkd_code, quiet = TRUE),
                    error = function(e) NULL)

## ============================================================
##  Helper functions
## ============================================================

run_sim <- function(params, scenario_list, sim_hours = 3 * 365.25 * 24) {
  if (is.null(mod_app)) return(NULL)

  mod_local <- mod_app %>%
    param(TKV0        = params$TKV0,
          EGFR0       = params$EGFR0,
          KGROW_TKV   = params$kgrow,
          KDECL_EGFR  = params$kdecl,
          CA2_SCALE   = params$ca2_scale,
          BP_SS       = params$BP0)

  results <- lapply(scenario_list, function(scen) {
    ev_dose <- scen$ev
    out <- mod_local %>%
      mrgsim_e(ev_dose, end = sim_hours, delta = 24, recover = "time") %>%
      as.data.frame() %>%
      mutate(scenario = scen$label, time_yr = time / (365.25 * 24))
    out
  })
  bind_rows(results)
}

build_scenarios <- function(input, sim_hours) {
  scen <- list()
  scen[[1]] <- list(label = "Placebo",
                    ev = ev(cmt = "AGUT", amt = 0, time = 0))
  if (input$use_tolv) {
    dose_am <- input$tolv_dose * 0.75
    dose_pm <- input$tolv_dose * 0.25
    scen[[2]] <- list(
      label = paste0("Tolvaptan ", input$tolv_dose, " mg/day"),
      ev = ev(cmt = "AGUT", amt = dose_am,
              time = seq(0, sim_hours, by = 24)) +
           ev(cmt = "AGUT", amt = dose_pm,
              time = seq(8, sim_hours + 8, by = 24))
    )
  }
  if (input$use_ever) {
    scen[[length(scen) + 1]] <- list(
      label = paste0("Everolimus ", input$ever_dose, " mg/day"),
      ev = ev(cmt = "EGUT", amt = input$ever_dose,
              time = seq(0, sim_hours, by = 24))
    )
  }
  if (input$use_acei) {
    scen[[length(scen) + 1]] <- list(
      label = "ACEi/ARB",
      ev = ev(cmt = "ACEI_GUT", amt = 10,
              time = seq(0, sim_hours, by = 24))
    )
  }
  if (input$use_combo && input$use_tolv && input$use_acei) {
    dose_am <- input$tolv_dose * 0.75
    dose_pm <- input$tolv_dose * 0.25
    scen[[length(scen) + 1]] <- list(
      label = "Tolvaptan + ACEi Combo",
      ev = ev(cmt = "AGUT", amt = dose_am,
              time = seq(0, sim_hours, by = 24)) +
           ev(cmt = "AGUT", amt = dose_pm,
              time = seq(8, sim_hours + 8, by = 24)) +
           ev(cmt = "ACEI_GUT", amt = 10,
              time = seq(0, sim_hours, by = 24))
    )
  }
  scen
}

mayo_class <- function(httkv) {
  if      (httkv < 210)  "Class 1A (very low)"
  else if (httkv < 320)  "Class 1B (low)"
  else if (httkv < 480)  "Class 1C (intermediate)"
  else if (httkv < 720)  "Class 1D (high)"
  else                   "Class 1E (very high)"
}

propkd_score <- function(gene, trunc_mut, hematuria, proteinuria,
                          hypertension, tkv_mayo) {
  score <- 0
  score <- score + if (gene == "PKD1" && trunc_mut)  4 else if (gene == "PKD1") 2 else 0
  score <- score + if (hematuria)   1 else 0
  score <- score + if (proteinuria) 1 else 0
  score <- score + if (hypertension) 2 else 0
  score <- score + if (tkv_mayo == "Class 1D (high)") 1 else if (tkv_mayo == "Class 1E (very high)") 2 else 0
  score
}

esrd_estimate <- function(propkd) {
  # Based on Cornec-Le Gall E et al. JASN 2016
  if      (propkd <= 3) "~70-80 years"
  else if (propkd <= 6) "~55-65 years"
  else                  "<55 years (high risk)"
}


## ============================================================
##  UI
## ============================================================

ui <- navbarPage(
  title = "ADPKD QSP Dashboard",
  theme = bs_theme(bootswatch = "flatly"),

  ## ---- Tab 1: Patient Profile ----
  tabPanel("1. Patient Profile",
    sidebarLayout(
      sidebarPanel(
        h4("Patient Characteristics"),
        numericInput("age",    "Age (years)",           value = 38, min = 18, max = 70),
        selectInput("sex",     "Sex",                   choices = c("Female", "Male")),
        numericInput("bwt",    "Body Weight (kg)",      value = 70, min = 30, max = 150),
        numericInput("height", "Height (cm)",           value = 170, min = 140, max = 200),
        hr(),
        h4("Genotype"),
        selectInput("gene",     "Gene",
                    choices = c("PKD1 truncating", "PKD1 non-truncating", "PKD2", "Unknown")),
        checkboxInput("trunc_mut", "Truncating mutation confirmed", TRUE),
        hr(),
        h4("Baseline Values"),
        numericInput("TKV0_inp",  "Baseline TKV (mL)",  value = 1500, min = 200, max = 6000),
        numericInput("EGFR0_inp", "Baseline eGFR (mL/min/1.73m²)", value = 70, min = 15, max = 120),
        numericInput("BP0_inp",   "Baseline SBP (mmHg)", value = 135, min = 100, max = 180),
        hr(),
        h4("Clinical Features"),
        checkboxInput("hematuria",   "History of hematuria", FALSE),
        checkboxInput("proteinuria", "Proteinuria > 0.3 g/g", FALSE),
        checkboxInput("htn_diag",    "Hypertension diagnosed", TRUE)
      ),
      mainPanel(
        fluidRow(
          column(6,
            wellPanel(
              h4("Calculated Indices"),
              verbatimTextOutput("patient_summary")
            )
          ),
          column(6,
            wellPanel(
              h4("PROPKD Score"),
              verbatimTextOutput("propkd_out")
            )
          )
        ),
        fluidRow(
          column(12,
            wellPanel(
              h4("Mayo Imaging Classification"),
              plotOutput("mayo_plot", height = "250px")
            )
          )
        ),
        fluidRow(
          column(12,
            wellPanel(
              h4("Monitoring Recommendations"),
              tableOutput("monitoring_table")
            )
          )
        )
      )
    )
  ),

  ## ---- Tab 2: Drug PK ----
  tabPanel("2. Drug PK",
    sidebarLayout(
      sidebarPanel(
        h4("Pharmacokinetic Simulation"),
        selectInput("pk_drug", "Drug",
                    choices = c("Tolvaptan", "Everolimus", "Octreotide LAR")),
        conditionalPanel(
          condition = "input.pk_drug == 'Tolvaptan'",
          sliderInput("tolv_pk_dose",  "Total Daily Dose (mg)", 60, 15, 120, 15),
          helpText("Split: 75% AM + 25% PM (8h interval)")
        ),
        conditionalPanel(
          condition = "input.pk_drug == 'Everolimus'",
          sliderInput("ever_pk_dose", "Daily Dose (mg)", 2.5, 0.5, 10, 0.5)
        ),
        conditionalPanel(
          condition = "input.pk_drug == 'Octreotide LAR'",
          sliderInput("oct_dose_mg",  "Depot Dose (mg)", 30, 10, 40, 10),
          helpText("Monthly IM injection (28-day interval)")
        ),
        sliderInput("pk_days", "Simulation Duration (days)", 14, 3, 56, 1),
        actionButton("run_pk", "Run PK Simulation", class = "btn-primary btn-lg")
      ),
      mainPanel(
        plotlyOutput("pk_plot", height = "350px"),
        hr(),
        h4("PK Parameters & Target Ranges"),
        tableOutput("pk_table"),
        h4("PK Metrics (Steady-state)"),
        verbatimTextOutput("pk_metrics")
      )
    )
  ),

  ## ---- Tab 3: PD Biomarkers ----
  tabPanel("3. PD Biomarkers",
    sidebarLayout(
      sidebarPanel(
        h4("PD Simulation Settings"),
        sliderInput("pd_tolv_dose", "Tolvaptan Dose (mg/day)", 60, 15, 120, 15),
        sliderInput("pd_ever_dose", "Everolimus Dose (mg/day)", 2.5, 0, 10, 0.5),
        sliderInput("pd_sim_days",  "Duration (days)", 30, 7, 90, 7),
        actionButton("run_pd", "Run PD Simulation", class = "btn-primary btn-lg"),
        hr(),
        h5("Key PD Targets"),
        p("• V2R occupancy target: ≥90%"),
        p("• Urine Osm reduction: <300 mOsm/kg"),
        p("• mTOR inhibition (S6K1): ≥80%")
      ),
      mainPanel(
        fluidRow(
          column(6, plotlyOutput("v2r_plot",    height = "280px")),
          column(6, plotlyOutput("uosm_plot",   height = "280px"))
        ),
        fluidRow(
          column(6, plotlyOutput("camp_plot",   height = "280px")),
          column(6, plotlyOutput("mtor_pd_plot",height = "280px"))
        ),
        hr(),
        h4("PD Summary at Selected Timepoint"),
        sliderInput("pd_time_select", "Timepoint (days)", 14, 1, 90, 1),
        tableOutput("pd_summary_table")
      )
    )
  ),

  ## ---- Tab 4: Disease Progression ----
  tabPanel("4. Disease Progression",
    sidebarLayout(
      sidebarPanel(
        h4("Simulation Parameters"),
        sliderInput("dp_years", "Simulation Duration (years)", 3, 1, 10, 1),
        hr(),
        h4("Treatment"),
        checkboxInput("use_tolv",  "Tolvaptan", TRUE),
        sliderInput("tolv_dose", "Tolvaptan dose (mg/day)", 60, 15, 120, 15),
        checkboxInput("use_ever",  "Everolimus (mTOR)", FALSE),
        sliderInput("ever_dose", "Everolimus dose (mg/day)", 2.5, 0.5, 5, 0.5),
        checkboxInput("use_acei",  "ACEi/ARB", FALSE),
        checkboxInput("use_combo", "Show combination arm", FALSE),
        hr(),
        h4("Disease Parameters"),
        sliderInput("ca2_scale_sl", "ADPKD cAMP amplification (1=mild, 2=severe)",
                    1.5, 1.0, 2.5, 0.1),
        sliderInput("kgrow_sl", "TKV growth multiplier (1=average)",
                    1.0, 0.5, 2.0, 0.1),
        actionButton("run_dp", "Simulate Progression", class = "btn-primary btn-lg")
      ),
      mainPanel(
        fluidRow(
          column(6, plotlyOutput("tkv_plot",    height = "300px")),
          column(6, plotlyOutput("egfr_dp_plot",height = "300px"))
        ),
        fluidRow(
          column(6, plotlyOutput("bp_dp_plot",  height = "300px")),
          column(6, plotlyOutput("ckd_stage_plot", height = "300px"))
        )
      )
    )
  ),

  ## ---- Tab 5: Scenario Comparison ----
  tabPanel("5. Scenario Comparison",
    sidebarLayout(
      sidebarPanel(
        h4("Select Treatments to Compare"),
        checkboxInput("sc_placebo", "Placebo (control)", TRUE),
        checkboxInput("sc_tolv60",  "Tolvaptan 60 mg/day", TRUE),
        checkboxInput("sc_tolv120", "Tolvaptan 120 mg/day", TRUE),
        checkboxInput("sc_ever",    "Everolimus 2.5 mg/day", FALSE),
        checkboxInput("sc_acei",    "ACEi/ARB (ramipril)", FALSE),
        sliderInput("sc_years", "Time horizon (years)", 3, 1, 5, 1),
        actionButton("run_sc", "Run Comparison", class = "btn-primary btn-lg"),
        hr(),
        h5("Calibration Reference"),
        p("TEMPO 3:4 (Torres 2012):"),
        p("TKV: +4.4% vs +8.0%/yr"),
        p("REPRISE (Torres 2017):"),
        p("eGFR change: -2.3 vs -3.6 mL/min"),
        p("SIRENA (Serra/Walz 2010):"),
        p("TKV: -4.9% vs +9.7% (1 yr)")
      ),
      mainPanel(
        fluidRow(
          column(6, plotlyOutput("sc_tkv_plot",   height = "300px")),
          column(6, plotlyOutput("sc_egfr_plot",  height = "300px"))
        ),
        hr(),
        h4("Outcome Comparison Table (at selected time horizon)"),
        DTOutput("sc_table"),
        hr(),
        h4("Forest Plot: TKV Growth Rate vs. Placebo"),
        plotOutput("forest_plot", height = "300px")
      )
    )
  ),

  ## ---- Tab 6: Biomarker & Risk Assessment ----
  tabPanel("6. Biomarker & Risk",
    sidebarLayout(
      sidebarPanel(
        h4("Risk Assessment Inputs"),
        selectInput("risk_gene", "Gene (mutation type)",
                    choices = c("PKD1 truncating" = "PKD1T",
                                "PKD1 non-truncating" = "PKD1N",
                                "PKD2" = "PKD2",
                                "Unknown" = "UNK")),
        numericInput("risk_age",   "Age at diagnosis (years)", 35, 18, 65),
        numericInput("risk_TKV",   "TKV at diagnosis (mL)",  1200, 200, 5000),
        numericInput("risk_ht",    "Height (m)",               1.70, 1.40, 2.00, 0.01),
        numericInput("risk_eGFR",  "eGFR (mL/min/1.73m²)",     72, 20, 120),
        checkboxInput("risk_hematuria",   "History of hematuria", FALSE),
        checkboxInput("risk_proteinuria", "Proteinuria > 0.3 g/g", FALSE),
        checkboxInput("risk_htn",         "Hypertension", TRUE),
        actionButton("calc_risk", "Calculate Risk", class = "btn-primary btn-lg")
      ),
      mainPanel(
        fluidRow(
          column(6,
            wellPanel(
              h4("Mayo Imaging Classification"),
              verbatimTextOutput("mayo_risk_out")
            )
          ),
          column(6,
            wellPanel(
              h4("PROPKD Score"),
              verbatimTextOutput("propkd_risk_out")
            )
          )
        ),
        fluidRow(
          column(12,
            wellPanel(
              h4("Predicted Disease Trajectory"),
              plotOutput("risk_trajectory_plot", height = "300px")
            )
          )
        ),
        fluidRow(
          column(6,
            wellPanel(
              h4("ESRD Prediction"),
              verbatimTextOutput("esrd_out")
            )
          ),
          column(6,
            wellPanel(
              h4("Monitoring & Treatment"),
              verbatimTextOutput("treatment_rec")
            )
          )
        ),
        fluidRow(
          column(12,
            h4("Biomarker Reference Ranges"),
            tableOutput("biomarker_ref_table")
          )
        )
      )
    )
  )
)


## ============================================================
##  SERVER
## ============================================================

server <- function(input, output, session) {

  ## ---- Tab 1: Patient Profile ----

  httkv_val <- reactive({
    round(input$TKV0_inp / (input$height / 100), 1)
  })

  mayo_val <- reactive({ mayo_class(httkv_val()) })

  output$patient_summary <- renderPrint({
    httkv <- httkv_val()
    bmi   <- input$bwt / (input$height / 100)^2
    cat(sprintf("BMI:        %.1f kg/m²\n", bmi))
    cat(sprintf("htTKV:      %.0f mL/m\n",  httkv))
    cat(sprintf("Mayo Class: %s\n",  mayo_val()))
    gene_lab <- switch(input$gene,
      "PKD1 truncating"     = "PKD1 (truncating) — severe",
      "PKD1 non-truncating" = "PKD1 (non-truncating) — moderate",
      "PKD2"                = "PKD2 — milder course",
      "Unknown"             = "Gene unknown")
    cat(sprintf("Genotype:   %s\n", gene_lab))
    cat(sprintf("eGFR:       %.0f mL/min/1.73m²\n", input$EGFR0_inp))
    ckd <- if (input$EGFR0_inp >= 90) "CKD Stage 1"
           else if (input$EGFR0_inp >= 60) "CKD Stage 2"
           else if (input$EGFR0_inp >= 45) "CKD Stage 3a"
           else if (input$EGFR0_inp >= 30) "CKD Stage 3b"
           else if (input$EGFR0_inp >= 15) "CKD Stage 4"
           else "CKD Stage 5 (ESRD)"
    cat(sprintf("CKD Stage:  %s\n", ckd))
  })

  output$propkd_out <- renderPrint({
    gene_type <- if (grepl("PKD1", input$gene)) {
      if (input$trunc_mut) "PKD1" else "PKD1_NT"
    } else if (input$gene == "PKD2") "PKD2" else "UNK"
    sc <- 0
    if (gene_type == "PKD1"   ) sc <- sc + 4
    if (gene_type == "PKD1_NT") sc <- sc + 2
    if (input$hematuria)    sc <- sc + 1
    if (input$proteinuria)  sc <- sc + 1
    if (input$htn_diag)     sc <- sc + 2
    mayo_sc <- mayo_val()
    if (grepl("1D", mayo_sc)) sc <- sc + 1
    if (grepl("1E", mayo_sc)) sc <- sc + 2
    risk_cat <- if (sc <= 3) "Low risk (score 0-3)"
                else if (sc <= 6) "Intermediate risk (4-6)"
                else "High risk (7-9)"
    cat(sprintf("PROPKD Score: %d / 9\n", sc))
    cat(sprintf("Risk Category: %s\n", risk_cat))
    cat("\nEstimated ESRD:\n")
    cat(sprintf("  %s\n", esrd_estimate(sc)))
  })

  output$mayo_plot <- renderPlot({
    httkv <- httkv_val()
    breaks <- c(210, 320, 480, 720)
    classes <- c("1A", "1B", "1C", "1D", "1E")
    df_mayo <- data.frame(
      class = classes,
      lower = c(0, breaks),
      upper = c(breaks, 1500)
    )
    cur_class <- mayo_val()
    ggplot(df_mayo, aes(x = class)) +
      geom_bar(aes(y = upper - lower, fill = class),
               stat = "identity", width = 0.7) +
      geom_hline(yintercept = httkv - df_mayo$lower[1], linetype = "dashed",
                 color = "red", linewidth = 1.5) +
      scale_fill_manual(
        values = c("1A" = "#4CAF50", "1B" = "#8BC34A",
                   "1C" = "#FFC107", "1D" = "#FF9800", "1E" = "#F44336")
      ) +
      annotate("text", x = 3, y = 400,
               label = paste("Patient:", cur_class, "\nhtTKV:", httkv, "mL/m"),
               color = "red", size = 4) +
      labs(title = "Mayo Imaging Classification (htTKV bands)",
           x = "Class", y = "htTKV range (mL/m)", fill = "Class") +
      theme_bw(base_size = 12) + theme(legend.position = "none")
  })

  output$monitoring_table <- renderTable({
    data.frame(
      Parameter     = c("TKV (MRI)", "eGFR", "Blood pressure", "Urine osmolality",
                        "ALT/AST (tolvaptan)", "Urine protein", "Urinary AQP2"),
      Interval      = c("Every 12 months", "Every 6 months", "Every visit",
                        "At initiation", "Monthly × 18 months, then Q6M",
                        "Annually", "Optional (research)"),
      `Target/Reference` = c(httkv_val(), round(input$EGFR0_inp, 0),
                               "<130/80 mmHg", "<300 mOsm/kg",
                               "<3× ULN", "<0.3 g/g", "Not standardised")
    )
  })


  ## ---- Tab 2: Drug PK ----

  pk_data <- eventReactive(input$run_pk, {
    req(mod_app)
    sim_h <- input$pk_days * 24
    if (input$pk_drug == "Tolvaptan") {
      d_am <- input$tolv_pk_dose * 0.75
      d_pm <- input$tolv_pk_dose * 0.25
      ev_d <- ev(cmt = "AGUT", amt = d_am, time = seq(0, sim_h, by = 24)) +
              ev(cmt = "AGUT", amt = d_pm, time = seq(8, sim_h + 8, by = 24))
      yvar  <- "Cp_tolv"; yunits <- "ng/mL"; target_lo <- 50; target_hi <- 300
    } else if (input$pk_drug == "Everolimus") {
      ev_d <- ev(cmt = "EGUT", amt = input$ever_pk_dose, time = seq(0, sim_h, by = 24))
      yvar  <- "Ce_ever"; yunits <- "ng/mL"; target_lo <- 3; target_hi <- 8
    } else {
      ev_d <- ev(cmt = "OCTDEP", amt = input$oct_dose_mg * 1000,
                 time = seq(0, sim_h, by = 24 * 28))
      yvar  <- "Cp_oct"; yunits <- "ng/mL"; target_lo <- 0.5; target_hi <- 5
    }
    out <- mod_app %>%
      mrgsim_e(ev_d, end = sim_h, delta = 1) %>%
      as.data.frame()
    list(data = out, yvar = yvar, yunits = yunits,
         target_lo = target_lo, target_hi = target_hi)
  })

  output$pk_plot <- renderPlotly({
    req(pk_data())
    d <- pk_data()
    y <- d$data[[d$yvar]]
    p <- ggplot(d$data, aes(x = time / 24, y = .data[[d$yvar]])) +
      geom_line(color = "#1565C0", linewidth = 1) +
      geom_hline(yintercept = d$target_lo, linetype = "dashed", color = "green") +
      geom_hline(yintercept = d$target_hi, linetype = "dashed", color = "red") +
      labs(x = "Time (days)", y = paste0("Concentration (", d$yunits, ")"),
           title = paste(input$pk_drug, "PK Profile")) +
      theme_bw()
    ggplotly(p)
  })

  output$pk_table <- renderTable({
    data.frame(
      Parameter   = c("Oral bioavailability", "t½", "CL", "Vd", "Dosing"),
      Tolvaptan   = c("56%", "~8 h", "4 L/h", "30 L", "AM+PM split"),
      Everolimus  = c("16%", "~30 h", "14 L/h", "70 L", "Once daily"),
      `Octreotide LAR` = c("~80% (depot)", "~28 days", "0.5 L/h", "8 L", "Monthly IM")
    )
  })

  output$pk_metrics <- renderPrint({
    req(pk_data())
    d    <- pk_data()
    dat  <- d$data
    conc <- dat[[d$yvar]]
    t    <- dat$time
    cat(sprintf("Cmax (overall):  %.1f %s\n", max(conc),   d$yunits))
    cat(sprintf("Cmin (overall):  %.1f %s\n", min(conc),   d$yunits))
    n24 <- dat %>% filter(time >= 24 & time <= 48)
    cat(sprintf("Cmax (Day 2):    %.1f %s\n", max(n24[[d$yvar]]), d$yunits))
    cat(sprintf("Tmax (Day 2):    %.1f h\n",  n24$time[which.max(n24[[d$yvar]])] - 24))
    cat(sprintf("AUC (Day 2):     %.0f %s·h\n",
        trapz(n24$time, n24[[d$yvar]]), d$yunits))
  })

  # Simple trapezoid helper
  trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)


  ## ---- Tab 3: PD Biomarkers ----

  pd_data <- eventReactive(input$run_pd, {
    req(mod_app)
    sim_h <- input$pd_sim_days * 24
    d_am  <- input$pd_tolv_dose * 0.75
    d_pm  <- input$pd_tolv_dose * 0.25
    ev_d  <- ev(cmt = "AGUT", amt = d_am, time = seq(0, sim_h, by = 24)) +
             ev(cmt = "AGUT", amt = d_pm, time = seq(8, sim_h + 8, by = 24)) +
             ev(cmt = "EGUT", amt = input$pd_ever_dose, time = seq(0, sim_h, by = 24))
    out <- mod_app %>%
      mrgsim_e(ev_d, end = sim_h, delta = 1) %>%
      as.data.frame() %>%
      mutate(time_d = time / 24)
    out
  })

  output$v2r_plot <- renderPlotly({
    req(pd_data())
    p <- ggplot(pd_data(), aes(x = time_d, y = V2R_OCC)) +
      geom_line(color = "#1565C0", linewidth = 1.1) +
      geom_hline(yintercept = 90, linetype = "dashed", color = "green") +
      annotate("text", x = input$pd_sim_days * 0.7, y = 88, label = "Target ≥90%",
               color = "darkgreen", size = 3.5) +
      labs(title = "V2R Occupancy by Tolvaptan",
           x = "Time (days)", y = "V2R Occupancy (%)") +
      ylim(0, 105) + theme_bw()
    ggplotly(p)
  })

  output$uosm_plot <- renderPlotly({
    req(pd_data())
    p <- ggplot(pd_data(), aes(x = time_d, y = Uosm)) +
      geom_line(color = "#E65100", linewidth = 1.1) +
      geom_hline(yintercept = 300, linetype = "dashed", color = "blue") +
      annotate("text", x = input$pd_sim_days * 0.7, y = 310, label = "Target <300",
               color = "blue", size = 3.5) +
      labs(title = "Urine Osmolality (Tolvaptan PD)",
           x = "Time (days)", y = "Uosm (mOsm/kg)") +
      theme_bw()
    ggplotly(p)
  })

  output$camp_plot <- renderPlotly({
    req(pd_data())
    p <- ggplot(pd_data(), aes(x = time_d, y = CAMP_rel)) +
      geom_line(color = "#FF9800", linewidth = 1.1) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey") +
      labs(title = "Relative cAMP (Collecting Duct)",
           x = "Time (days)", y = "cAMP (ADPKD baseline = 1)") +
      theme_bw()
    ggplotly(p)
  })

  output$mtor_pd_plot <- renderPlotly({
    req(pd_data())
    p <- ggplot(pd_data(), aes(x = time_d, y = mTOR_INH)) +
      geom_line(color = "#6A1B9A", linewidth = 1.1) +
      geom_hline(yintercept = 80, linetype = "dashed", color = "green") +
      annotate("text", x = input$pd_sim_days * 0.7, y = 78, label = "Tgt ≥80%",
               color = "darkgreen", size = 3.5) +
      labs(title = "mTOR Inhibition (Everolimus PD)",
           x = "Time (days)", y = "mTOR Inhibition (%)") +
      ylim(0, 105) + theme_bw()
    ggplotly(p)
  })

  output$pd_summary_table <- renderTable({
    req(pd_data())
    t_h <- input$pd_time_select * 24
    closest <- pd_data() %>% filter(abs(time - t_h) == min(abs(time - t_h))) %>% head(1)
    data.frame(
      Biomarker = c("Tolvaptan Cp", "V2R Occupancy", "Urine Osmolality",
                    "Relative cAMP", "Everolimus Ce", "mTOR Inhibition"),
      Value     = c(sprintf("%.1f ng/mL",  closest$Cp_tolv),
                    sprintf("%.1f%%",        closest$V2R_OCC),
                    sprintf("%.0f mOsm/kg",  closest$Uosm),
                    sprintf("%.2f (rel)",    closest$CAMP_rel),
                    sprintf("%.2f ng/mL",   closest$Ce_ever),
                    sprintf("%.1f%%",        closest$mTOR_INH)),
      Target    = c("50-300 ng/mL", "≥90%", "<300 mOsm/kg",
                    "<1.0 (↓cAMP)", "3-8 ng/mL", "≥80%"),
      Status    = c(
        if (closest$Cp_tolv  >= 50 && closest$Cp_tolv <= 300) "✓ In range" else "— Out",
        if (closest$V2R_OCC  >= 90) "✓ Target met" else "— Below",
        if (closest$Uosm     <  300) "✓ Target met" else "— Above",
        if (closest$CAMP_rel <  1.0) "✓ Reduced" else "— Elevated",
        if (closest$Ce_ever  >= 3 && closest$Ce_ever <= 8) "✓ In window" else "— Adj.",
        if (closest$mTOR_INH >= 80) "✓ Target met" else "— Below"
      )
    )
  })


  ## ---- Tab 4: Disease Progression ----

  dp_data <- eventReactive(input$run_dp, {
    req(mod_app)
    sim_h <- input$dp_years * 365.25 * 24
    params <- list(
      TKV0      = input$TKV0_inp,
      EGFR0     = input$EGFR0_inp,
      kgrow     = 6.28e-6 * input$kgrow_sl,
      kdecl     = 4e-4,
      ca2_scale = input$ca2_scale_sl,
      BP0       = input$BP0_inp
    )
    scenarios <- build_scenarios(input, sim_h)
    run_sim(params, scenarios, sim_h)
  })

  output$tkv_plot <- renderPlotly({
    req(dp_data())
    p <- ggplot(dp_data(), aes(x = time_yr, y = TKV_L, color = scenario)) +
      geom_line(linewidth = 1.1) +
      labs(title = "TKV Over Time", x = "Years", y = "TKV (L)", color = NULL) +
      theme_bw() + theme(legend.position = "bottom",
                         legend.text = element_text(size = 8))
    ggplotly(p) %>% layout(legend = list(orientation = "h"))
  })

  output$egfr_dp_plot <- renderPlotly({
    req(dp_data())
    p <- ggplot(dp_data(), aes(x = time_yr, y = eGFR, color = scenario)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = c(90, 60, 45, 30, 15),
                 linetype = "dotted", color = "grey60") +
      labs(title = "eGFR Over Time", x = "Years",
           y = "eGFR (mL/min/1.73m²)", color = NULL) +
      theme_bw() + theme(legend.position = "bottom",
                         legend.text = element_text(size = 8))
    ggplotly(p) %>% layout(legend = list(orientation = "h"))
  })

  output$bp_dp_plot <- renderPlotly({
    req(dp_data())
    p <- ggplot(dp_data(), aes(x = time_yr, y = BP, color = scenario)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 130, linetype = "dashed", color = "red") +
      labs(title = "Blood Pressure", x = "Years", y = "SBP (mmHg)", color = NULL) +
      theme_bw() + theme(legend.position = "bottom",
                         legend.text = element_text(size = 8))
    ggplotly(p) %>% layout(legend = list(orientation = "h"))
  })

  output$ckd_stage_plot <- renderPlotly({
    req(dp_data())
    ckd_df <- dp_data() %>%
      mutate(CKD = case_when(
        eGFR >= 90 ~ 1, eGFR >= 60 ~ 2, eGFR >= 45 ~ 3,
        eGFR >= 30 ~ 4, eGFR >= 15 ~ 5, TRUE ~ 5
      ))
    p <- ggplot(ckd_df, aes(x = time_yr, y = CKD, color = scenario)) +
      geom_line(linewidth = 1.1) +
      scale_y_continuous(breaks = 1:5,
                         labels = c("1 (≥90)", "2 (60-89)", "3 (30-59)",
                                    "4 (15-29)", "5 (<15)")) +
      labs(title = "CKD Stage Progression",
           x = "Years", y = "CKD Stage", color = NULL) +
      theme_bw() + theme(legend.position = "bottom",
                         legend.text = element_text(size = 8))
    ggplotly(p) %>% layout(legend = list(orientation = "h"))
  })


  ## ---- Tab 5: Scenario Comparison ----

  sc_data <- eventReactive(input$run_sc, {
    req(mod_app)
    sim_h <- input$sc_years * 365.25 * 24
    params <- list(TKV0 = input$TKV0_inp, EGFR0 = input$EGFR0_inp,
                   kgrow = 6.28e-6, kdecl = 4e-4, ca2_scale = 1.5,
                   BP0 = input$BP0_inp)
    scens <- list()
    if (input$sc_placebo)
      scens[[1]] <- list(label = "Placebo",
                         ev = ev(cmt = "AGUT", amt = 0, time = 0))
    if (input$sc_tolv60) {
      scens[[length(scens)+1]] <- list(label = "Tolvaptan 60 mg/day",
        ev = ev(cmt = "AGUT", amt = 45, time = seq(0, sim_h, by = 24)) +
             ev(cmt = "AGUT", amt = 15, time = seq(8, sim_h+8, by = 24)))
    }
    if (input$sc_tolv120) {
      scens[[length(scens)+1]] <- list(label = "Tolvaptan 120 mg/day",
        ev = ev(cmt = "AGUT", amt = 90, time = seq(0, sim_h, by = 24)) +
             ev(cmt = "AGUT", amt = 30, time = seq(8, sim_h+8, by = 24)))
    }
    if (input$sc_ever) {
      scens[[length(scens)+1]] <- list(label = "Everolimus 2.5 mg/day",
        ev = ev(cmt = "EGUT", amt = 2.5, time = seq(0, sim_h, by = 24)))
    }
    if (input$sc_acei) {
      scens[[length(scens)+1]] <- list(label = "ACEi/ARB",
        ev = ev(cmt = "ACEI_GUT", amt = 10, time = seq(0, sim_h, by = 24)))
    }
    if (length(scens) == 0) {
      scens[[1]] <- list(label = "Placebo",
                         ev = ev(cmt = "AGUT", amt = 0, time = 0))
    }
    run_sim(params, scens, sim_h)
  })

  output$sc_tkv_plot <- renderPlotly({
    req(sc_data())
    p <- ggplot(sc_data(), aes(x = time_yr, y = TKV_L, color = scenario)) +
      geom_line(linewidth = 1.1) +
      labs(title = "TKV Comparison", x = "Years", y = "TKV (L)", color = NULL) +
      theme_bw() + theme(legend.position = "bottom",
                         legend.text = element_text(size = 8))
    ggplotly(p) %>% layout(legend = list(orientation = "h"))
  })

  output$sc_egfr_plot <- renderPlotly({
    req(sc_data())
    p <- ggplot(sc_data(), aes(x = time_yr, y = eGFR, color = scenario)) +
      geom_line(linewidth = 1.1) +
      labs(title = "eGFR Comparison", x = "Years",
           y = "eGFR (mL/min/1.73m²)", color = NULL) +
      theme_bw() + theme(legend.position = "bottom",
                         legend.text = element_text(size = 8))
    ggplotly(p) %>% layout(legend = list(orientation = "h"))
  })

  output$sc_table <- renderDT({
    req(sc_data())
    tbl <- sc_data() %>%
      group_by(scenario) %>%
      slice_max(order_by = time_yr, n = 1) %>%
      summarise(
        `TKV (L)` = round(last(TKV_L), 2),
        `TKV Change (%)` = round(last(TKV_pct), 1),
        `eGFR (mL/min)` = round(last(eGFR), 1),
        `eGFR Δ` = round(input$EGFR0_inp - last(eGFR), 1),
        `BP (mmHg)` = round(last(BP), 1),
        `cAMP (rel)` = round(last(CAMP_rel), 2),
        .groups = "drop"
      ) %>%
      rename(Scenario = scenario)
    datatable(tbl, options = list(dom = "t"), rownames = FALSE)
  })

  output$forest_plot <- renderPlot({
    req(sc_data())
    summary_df <- sc_data() %>%
      group_by(scenario) %>%
      slice_max(time_yr, n = 1) %>%
      summarise(TKV_pct = last(TKV_pct), .groups = "drop")
    plac_val <- summary_df %>% filter(scenario == "Placebo") %>% pull(TKV_pct)
    if (length(plac_val) == 0) plac_val <- summary_df$TKV_pct[1]
    forest_df <- summary_df %>%
      mutate(diff_vs_placebo = TKV_pct - plac_val,
             lo = diff_vs_placebo - 5,
             hi = diff_vs_placebo + 5) %>%
      filter(scenario != "Placebo")
    ggplot(forest_df, aes(x = diff_vs_placebo, y = reorder(scenario, diff_vs_placebo))) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.3, color = "#1565C0") +
      geom_point(size = 4, color = "#1565C0") +
      labs(title = "TKV Growth vs. Placebo at 3 Years",
           x = "Difference in TKV Growth (%)", y = NULL) +
      theme_bw(base_size = 12)
  })


  ## ---- Tab 6: Risk Assessment ----

  risk_result <- eventReactive(input$calc_risk, {
    httkv_r <- input$risk_TKV / input$risk_ht
    mc <- mayo_class(httkv_r)
    gene_ok <- input$risk_gene
    sc <- 0
    if (gene_ok == "PKD1T") sc <- sc + 4
    if (gene_ok == "PKD1N") sc <- sc + 2
    if (input$risk_hematuria)   sc <- sc + 1
    if (input$risk_proteinuria) sc <- sc + 1
    if (input$risk_htn)         sc <- sc + 2
    if (grepl("1D", mc)) sc <- sc + 1
    if (grepl("1E", mc)) sc <- sc + 2
    list(httkv = httkv_r, mayo = mc, propkd = sc, esrd = esrd_estimate(sc))
  })

  output$mayo_risk_out <- renderPrint({
    req(risk_result())
    r <- risk_result()
    cat(sprintf("htTKV:      %.0f mL/m\n", r$httkv))
    cat(sprintf("Mayo Class: %s\n", r$mayo))
  })

  output$propkd_risk_out <- renderPrint({
    req(risk_result())
    r <- risk_result()
    rc <- if (r$propkd <= 3) "Low (0-3)"
          else if (r$propkd <= 6) "Intermediate (4-6)"
          else "High (7-9)"
    cat(sprintf("PROPKD Score: %d / 9\n", r$propkd))
    cat(sprintf("Risk Category: %s\n", rc))
  })

  output$risk_trajectory_plot <- renderPlot({
    req(risk_result(), mod_app)
    r  <- risk_result()
    sc <- r$propkd
    kgrow_adj <- 6.28e-6 * (1 + (sc - 4) * 0.15)   # scale by risk
    kdecl_adj <- 4.00e-4 * (1 + (sc - 4) * 0.10)
    sim_h <- 10 * 365.25 * 24
    params <- list(TKV0 = input$risk_TKV, EGFR0 = input$risk_eGFR,
                   kgrow = kgrow_adj, kdecl = kdecl_adj,
                   ca2_scale = 1.5, BP0 = 135)
    ev_plac <- ev(cmt = "AGUT", amt = 0, time = 0)
    ev_tolv <- ev(cmt = "AGUT", amt = 45, time = seq(0, sim_h, by = 24)) +
               ev(cmt = "AGUT", amt = 15, time = seq(8, sim_h + 8, by = 24))
    df_p <- run_sim(params, list(list(label = "Placebo",   ev = ev_plac)), sim_h)
    df_t <- run_sim(params, list(list(label = "Tolvaptan", ev = ev_tolv)), sim_h)
    df_all <- bind_rows(df_p, df_t)
    ggplot(df_all, aes(x = time_yr, y = eGFR, color = scenario)) +
      geom_line(linewidth = 1.3) +
      geom_hline(yintercept = c(60, 30, 15), linetype = "dotted", color = "grey") +
      scale_color_manual(values = c(Placebo = "#888888", Tolvaptan = "#1565C0")) +
      labs(title = sprintf("Predicted 10-year eGFR (PROPKD=%d)", r$propkd),
           x = "Years", y = "eGFR (mL/min/1.73m²)", color = NULL) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
  })

  output$esrd_out <- renderPrint({
    req(risk_result())
    r <- risk_result()
    cat("Predicted age at ESRD:\n")
    cat(sprintf("  %s\n\n", r$esrd))
    cat(sprintf("Current age:    %d years\n", input$risk_age))
    cat(sprintf("PROPKD score:   %d / 9\n",   r$propkd))
    cat(sprintf("Gene:           %s\n", input$risk_gene))
  })

  output$treatment_rec <- renderPrint({
    req(risk_result())
    r <- risk_result()
    cat("Treatment Recommendations:\n\n")
    if (r$propkd >= 4) {
      cat("★ TOLVAPTAN: Consider if:\n")
      cat("  - Age 18-55, eGFR ≥25 mL/min\n")
      cat("  - Rapidly progressive (Mayo 1C-1E)\n")
      cat("  - Monitor LFTs monthly × 18 mo\n\n")
    } else {
      cat("  - Tolvaptan: consider watchful waiting\n\n")
    }
    cat("★ ACEi/ARB: Recommended if BP >130/80\n")
    cat("  HALT-PKD target: <110/75 mmHg (young)\n\n")
    cat("★ EVEROLIMUS: Not standard (SIRENA, RAPYD)\n")
    cat("  Use only in clinical trials\n\n")
    cat("★ Lifestyle: high water intake (≥3L/day)\n")
    cat("  (reduces endogenous AVP → ↓cAMP)")
  })

  output$biomarker_ref_table <- renderTable({
    data.frame(
      Biomarker        = c("TKV (MRI)", "htTKV", "Annual TKV growth",
                           "eGFR", "Urine Osm (on tolvaptan)", "V2R Occupancy",
                           "S6K1 phospho (everolimus)", "Urinary AQP2"),
      `Normal / Target`= c("Varies by age", ">600 mL/m: high risk", "<5%/yr on treatment",
                           ">60 mL/min", "<300 mOsm/kg", "≥90%",
                           "≥80% reduction", "Decreases with tolvaptan"),
      `Clinical Use`   = c("Primary endpoint, eligibility",
                           "Mayo classification",
                           "Treatment efficacy",
                           "CKD staging", "V2R target engagement",
                           "PD biomarker", "PD biomarker (mTOR trial)",
                           "Exploratory PD")
    )
  })

}

## ============================================================
##  Run App
## ============================================================
shinyApp(ui = ui, server = server)
