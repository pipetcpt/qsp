################################################################################
## Atrial Fibrillation QSP — Interactive Shiny Dashboard
## ============================================================
## Tabs: Patient Profile | Drug PK | Electrophysiology PD |
##       Thromboembolism Risk | Scenario Comparison | Biomarker Dashboard
##
## Dependencies: shiny, shinydashboard, mrgsolve, ggplot2, dplyr, tidyr,
##               DT, plotly
################################################################################

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(mrgsolve)

# ==============================================================================
# INLINE mrgsolve MODEL (simplified, fast version for Shiny)
# ==============================================================================

af_model_inline <- '
$PARAM
ka_AMIO=0.06 CL_AMIO=3.5 V1_AMIO=40 V2_AMIO=4200 k12_AMIO=0.02 k21_AMIO=0.002 F_AMIO=0.46
ka_APIX=1.2  CL_APIX=3.3  V1_APIX=21  F_APIX=0.50
ka_METRO=1.5 CL_METRO=65  V1_METRO=290 F_METRO=0.40
IC50_AMIO_ERP=0.5 IC50_AMIO_HR=1.0 Emax_AMIO_ERP=0.35 Emax_AMIO_HR=0.45
IC50_METRO_HR=25  Emax_METRO_HR=0.40
IC50_APIX_FXa=0.08 Emax_APIX_FXa=0.95
AF0=0.60 ERP0=180 kfib=0.001 kfib_ERP=20 kAF_remod=0.005
HR0_AF=140 QTc0=400 kQTc_ERP=0.5 AngII0=1.0 ROS0=1.0
FXa0=1.0 Thrombin0=1.0 kStroke_base=0.035 kStroke_Thr=0.04 kNE_decay=0.1 NE0=1.0
kAngII_fib=0.003 kROS_fib=0.002 kThr_FXa=0.8 SMAD_base=0.5

$CMT GI_AMIO C1_AMIO C2_AMIO GI_APIX C1_APIX GI_METRO C1_METRO
     AF_BURDEN ERP Fibrosis QTc HR_AF AngII ROS FXa Thrombin STROKE_RISK NE IL6 SMAD23

$ODE
double ke_AMIO  = CL_AMIO  / V1_AMIO;
double ke_APIX  = CL_APIX  / V1_APIX;
double ke_METRO = CL_METRO / V1_METRO;

dxdt_GI_AMIO  = -ka_AMIO  * GI_AMIO;
dxdt_C1_AMIO  = F_AMIO  * ka_AMIO  * GI_AMIO  - ke_AMIO  * C1_AMIO  - k12_AMIO * C1_AMIO  + k21_AMIO * C2_AMIO;
dxdt_C2_AMIO  = k12_AMIO * C1_AMIO - k21_AMIO * C2_AMIO;
dxdt_GI_APIX  = -ka_APIX  * GI_APIX;
dxdt_C1_APIX  = F_APIX  * ka_APIX  * GI_APIX  - ke_APIX  * C1_APIX;
dxdt_GI_METRO = -ka_METRO * GI_METRO;
dxdt_C1_METRO = F_METRO * ka_METRO * GI_METRO - ke_METRO * C1_METRO;

double Cp_AMIO  = C1_AMIO  / V1_AMIO;
double Cp_APIX_ug = C1_APIX / V1_APIX;
double Cp_METRO_ng = C1_METRO / V1_METRO * 1000.0;

double Eff_AMIO_ERP  = Emax_AMIO_ERP  * Cp_AMIO  / (IC50_AMIO_ERP  + Cp_AMIO);
double Eff_AMIO_HR   = Emax_AMIO_HR   * Cp_AMIO  / (IC50_AMIO_HR   + Cp_AMIO);
double Eff_METRO_HR  = Emax_METRO_HR  * Cp_METRO_ng / (IC50_METRO_HR  + Cp_METRO_ng);
double Eff_APIX_FXa  = Emax_APIX_FXa  * Cp_APIX_ug / (IC50_APIX_FXa  + Cp_APIX_ug);

double HR_red = Eff_AMIO_HR + Eff_METRO_HR - Eff_AMIO_HR * Eff_METRO_HR;
if(HR_red > 0.8) HR_red = 0.8;

double dERP_remodel = -kAF_remod * AF_BURDEN * ERP;
double dERP_drug    = Eff_AMIO_ERP * ERP0;
double dERP_fib     = -kfib_ERP * Fibrosis / 24.0;
dxdt_ERP = (dERP_remodel + dERP_drug / 24.0 + dERP_fib) - 0.0001 * (ERP - (ERP0 - kfib_ERP * Fibrosis));

double ERP_eff_AF = 1.0 / (1.0 + exp((ERP - 200.0) / 20.0));
double Fib_eff_AF = 1.5 * Fibrosis;
double kAF_in  = 0.003 * ERP_eff_AF * (1.0 + Fib_eff_AF);
double kAF_out = 0.002 * (1.0 - ERP_eff_AF);
dxdt_AF_BURDEN = kAF_in * (1.0 - AF_BURDEN) - kAF_out * AF_BURDEN + 0.0001 * NE * AF_BURDEN * (1.0 - AF_BURDEN);
if(AF_BURDEN < 0.001 && dxdt_AF_BURDEN < 0) dxdt_AF_BURDEN = 0;
if(AF_BURDEN > 0.999 && dxdt_AF_BURDEN > 0) dxdt_AF_BURDEN = 0;

double kFib_in = kfib * (AngII * kAngII_fib + ROS * kROS_fib + SMAD23 * 0.002);
dxdt_Fibrosis = kFib_in * (1.0 - Fibrosis) - 0.00005 * Fibrosis;
if(Fibrosis < 0.001 && dxdt_Fibrosis < 0) dxdt_Fibrosis = 0;
if(Fibrosis > 0.999 && dxdt_Fibrosis > 0) dxdt_Fibrosis = 0;

double QTc_target = QTc0 + kQTc_ERP * (ERP - ERP0) + 15.0 * Eff_AMIO_ERP;
dxdt_QTc = 0.02 * (QTc_target - QTc);

double HR_target = HR0_AF * (1.0 - HR_red) * (1.0 + 0.3 * NE);
dxdt_HR_AF = 0.05 * (HR_target - HR_AF);

dxdt_AngII = 0.02 * AF_BURDEN * (2.0 - AngII) - 0.05 * (AngII - AngII0);
dxdt_ROS   = 0.015 * AF_BURDEN * AngII - 0.04 * (ROS - ROS0);
dxdt_SMAD23= 0.03 * AngII * (1.5 - SMAD23) - 0.02 * (SMAD23 - SMAD_base);
dxdt_IL6   = 0.01 * AF_BURDEN * (3.0 - IL6) - 0.03 * (IL6 - 1.0);

double FXa_prod = 0.1 * AF_BURDEN * (1.0 + 0.5 * Thrombin);
double FXa_inh  = Eff_APIX_FXa * FXa;
dxdt_FXa = FXa_prod - 0.15 * FXa - FXa_inh + 0.05 * (FXa0 - FXa);

dxdt_Thrombin = kThr_FXa * FXa * AF_BURDEN - 0.2 * Thrombin + 0.02 * (Thrombin0 - Thrombin);

double stroke_rate = kStroke_base * Thrombin * AF_BURDEN * 100.0 + kStroke_Thr * (Thrombin - 1.0);
if(stroke_rate < 0) stroke_rate = 0;
dxdt_STROKE_RISK = 0.01 * (stroke_rate - STROKE_RISK);

dxdt_NE = 0.01 * AF_BURDEN * (2.0 - NE) - kNE_decay * (NE - NE0);

$TABLE
double Cp_AMIO_out  = C1_AMIO  / V1_AMIO;
double Cp_APIX_out  = C1_APIX  / V1_APIX * 1000.0;
double Cp_METRO_out = C1_METRO / V1_METRO * 1000.0;
double AntiXa_pct   = Emax_APIX_FXa * (C1_APIX / V1_APIX) / (IC50_APIX_FXa + (C1_APIX / V1_APIX)) * 100.0;
double NT_proBNP    = AF_BURDEN * 500.0 + 200.0;
double CRP_out      = IL6 * 3.5;
double LA_diam      = 38.0 + 10.0 * Fibrosis;

$CAPTURE AF_BURDEN ERP QTc HR_AF STROKE_RISK Fibrosis
         Cp_AMIO_out Cp_APIX_out Cp_METRO_out
         FXa Thrombin AntiXa_pct NT_proBNP CRP_out LA_diam AngII ROS NE IL6
'

# Compile model once at startup
cat("Compiling AF mrgsolve model...\n")
af_mod <- tryCatch(
  mcode("AF_Shiny_QSP", af_model_inline, quiet = TRUE),
  error = function(e) {
    message("mrgsolve compile error: ", e$message)
    NULL
  }
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

run_scenario_sim <- function(
    af_type    = "persistent",   # paroxysmal / persistent / longstanding
    dose_amio  = 0,   # mg/day (0 = off)
    dose_apix  = 0,   # mg BID (0 = off)
    dose_metro = 0,   # mg BID (0 = off)
    t_days     = 365,
    dt_h       = 12,
    age        = 65,
    chads       = 3,
    comorbid    = c()
) {
  if (is.null(af_mod)) return(NULL)

  # Set baseline based on AF type
  af_init <- switch(af_type,
    "paroxysmal"   = 0.25,
    "persistent"   = 0.60,
    "longstanding" = 0.85,
    0.60
  )

  fib_init <- switch(af_type,
    "paroxysmal"   = 0.10,
    "persistent"   = 0.22,
    "longstanding" = 0.40,
    0.22
  )

  # Comorbidity adjustments
  angII_init <- 1.0
  if ("HTN" %in% comorbid)  angII_init <- angII_init + 0.3
  if ("DM"  %in% comorbid)  angII_init <- angII_init + 0.2
  if ("HF"  %in% comorbid)  angII_init <- angII_init + 0.4

  # Age effect on ERP
  erp_init <- 180 - (max(age - 50, 0)) * 0.3

  stroke_init <- (chads * 1.2) * af_init  # simplified CHA2DS2-VASc
  if (stroke_init < 0) stroke_init <- 0.1

  init_v <- c(
    GI_AMIO=0, C1_AMIO=0, C2_AMIO=0,
    GI_APIX=0, C1_APIX=0,
    GI_METRO=0, C1_METRO=0,
    AF_BURDEN=af_init, ERP=erp_init, Fibrosis=fib_init,
    QTc=410, HR_AF=138, AngII=angII_init, ROS=1.2,
    FXa=1.0, Thrombin=1.0, STROKE_RISK=stroke_init,
    NE=1.1, IL6=1.3, SMAD23=0.55
  )

  # Build event table
  ev_list <- list()

  if (dose_amio > 0) {
    # Loading 400 mg BID x 28d then maintenance
    load_times <- seq(0, 28*24-1, by=12)
    ev_list[["amio_load"]] <- data.frame(
      time=load_times, cmt=1, amt=400, evid=1)
    maint_times <- seq(28*24, t_days*24-1, by=24)
    ev_list[["amio_maint"]] <- data.frame(
      time=maint_times, cmt=1, amt=dose_amio, evid=1)
  }

  if (dose_apix > 0) {
    apix_times <- seq(0, t_days*24-1, by=12)
    ev_list[["apix"]] <- data.frame(
      time=apix_times, cmt=4, amt=dose_apix, evid=1)
  }

  if (dose_metro > 0) {
    metro_times <- seq(0, t_days*24-1, by=12)
    ev_list[["metro"]] <- data.frame(
      time=metro_times, cmt=6, amt=dose_metro, evid=1)
  }

  ev_df <- if (length(ev_list) > 0) {
    do.call(rbind, ev_list)
  } else {
    data.frame(time=0, cmt=1, amt=0, evid=0)
  }

  ev_df <- ev_df[order(ev_df$time), ]

  tryCatch({
    out <- mrgsim(
      af_mod,
      idata  = data.frame(ID=1),
      events = ev_df,
      init   = init_v,
      tgrid  = tgrid(0, t_days*24, dt_h),
      output = "df"
    )
    df <- as.data.frame(out)
    df$time_days <- df$time / 24
    df
  }, error = function(e) {
    message("Sim error: ", e$message)
    NULL
  })
}

calc_chads <- function(age, htn, dm, hf, prior_stroke, pad, female) {
  score <- 0
  if (prior_stroke) score <- score + 2
  if (age >= 75)    score <- score + 2
  if (age >= 65 & age < 75) score <- score + 1
  if (htn)    score <- score + 1
  if (dm)     score <- score + 1
  if (hf)     score <- score + 1
  if (pad)    score <- score + 1
  if (female) score <- score + 1
  score
}

chads_stroke_risk <- function(score) {
  # Approximate annual stroke risk by score (%)
  risks <- c(0, 0.2, 0.6, 1.5, 2.8, 4.0, 5.3, 6.6, 7.9, 9.6, 10)
  idx <- min(score + 1, length(risks))
  risks[idx]
}

# ==============================================================================
# UI
# ==============================================================================

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = span(icon("heartbeat"), "AF QSP Dashboard"),
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",        tabName = "tab_patient",   icon = icon("user-md")),
      menuItem("Drug PK",                tabName = "tab_pk",        icon = icon("flask")),
      menuItem("Electrophysiology PD",   tabName = "tab_ep",        icon = icon("wave-square")),
      menuItem("Thromboembolism Risk",   tabName = "tab_thrombo",   icon = icon("tint")),
      menuItem("Scenario Comparison",    tabName = "tab_scenarios", icon = icon("chart-bar")),
      menuItem("Biomarker Dashboard",    tabName = "tab_biomarker", icon = icon("microscope"))
    ),
    hr(),
    div(style="padding:10px; color:#ccc; font-size:11px;",
      strong("AF QSP Model v1.0"),
      br(), "mrgsolve + Shiny",
      br(), "Calibrated: AFFIRM, RACE,",
      br(), "ARISTOTLE, RE-LY"
    )
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background-color: #F5F7FA; }
        .box { border-radius: 8px; }
        .traffic-green  { background: #27AE60; color: white; border-radius: 6px; padding: 8px; text-align: center; font-weight: bold; }
        .traffic-yellow { background: #F39C12; color: white; border-radius: 6px; padding: 8px; text-align: center; font-weight: bold; }
        .traffic-red    { background: #E74C3C; color: white; border-radius: 6px; padding: 8px; text-align: center; font-weight: bold; }
        .metric-box { background: white; border-radius: 8px; padding: 12px; margin: 5px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); text-align: center; }
        .metric-value { font-size: 22px; font-weight: bold; color: #2C3E50; }
        .metric-label { font-size: 11px; color: #7F8C8D; }
      "))
    ),

    tabItems(

      # ======================================================================
      # TAB 1: PATIENT PROFILE
      # ======================================================================
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title=tagList(icon("user"), " Patient Demographics & Risk Factors"),
              width=4, status="primary", solidHeader=TRUE,
            sliderInput("age", "Age (years)", min=18, max=90, value=68, step=1),
            radioButtons("sex", "Sex", choices=c("Male"="male","Female"="female"), inline=TRUE),
            radioButtons("af_type", "AF Type",
              choices=c("Paroxysmal"="paroxysmal",
                        "Persistent"="persistent",
                        "Long-standing Persistent"="longstanding"),
              selected="persistent"),
            sliderInput("lvef", "LVEF (%)", min=15, max=75, value=55, step=1),
            sliderInput("la_diam_inp", "LA Diameter (mm)", min=30, max=65, value=44, step=1),
            sliderInput("egfr", "eGFR (mL/min/1.73m²)", min=10, max=130, value=72, step=1)
          ),

          box(title=tagList(icon("stethoscope"), " Comorbidities (CHA₂DS₂-VASc)"),
              width=4, status="warning", solidHeader=TRUE,
            checkboxInput("htn",         "Hypertension",              value=TRUE),
            checkboxInput("dm",          "Diabetes Mellitus",         value=FALSE),
            checkboxInput("hf",          "Heart Failure (HFrEF/HFpEF)", value=FALSE),
            checkboxInput("prior_stroke","Prior Stroke/TIA",          value=FALSE),
            checkboxInput("pad",         "PAD / Vascular Disease",    value=FALSE),
            checkboxInput("osa",         "Obstructive Sleep Apnea",   value=FALSE),
            checkboxInput("obesity_cb",  "Obesity (BMI ≥30)",         value=FALSE),
            hr(),
            strong("HAS-BLED Factors:"),
            checkboxInput("hasbled_htn",    "Uncontrolled HTN (SBP>160)", value=FALSE),
            checkboxInput("hasbled_renal",  "Abnormal Renal Function",    value=FALSE),
            checkboxInput("hasbled_bleed",  "Prior Bleeding History",     value=FALSE),
            checkboxInput("hasbled_alcohol","Alcohol use (≥8 units/wk)",  value=FALSE)
          ),

          box(title=tagList(icon("calculator"), " Risk Score Calculation"),
              width=4, status="success", solidHeader=TRUE,
            h4("CHA₂DS₂-VASc Score"),
            div(class="metric-box",
              div(class="metric-value", textOutput("chads_score_out")),
              div(class="metric-label", "CHA₂DS₂-VASc Score")
            ),
            div(class="metric-box",
              div(class="metric-value", textOutput("chads_risk_out")),
              div(class="metric-label", "Estimated Annual Stroke Risk")
            ),
            hr(),
            h4("HAS-BLED Score"),
            div(class="metric-box",
              div(class="metric-value", textOutput("hasbled_out")),
              div(class="metric-label", "HAS-BLED Score (bleeding risk)")
            ),
            hr(),
            h4("Anticoagulation Recommendation"),
            uiOutput("anticoag_rec"),
            hr(),
            h4("Renal Dose Guidance"),
            uiOutput("renal_guidance")
          )
        )
      ),

      # ======================================================================
      # TAB 2: DRUG PK
      # ======================================================================
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title=tagList(icon("pills"), " Dosing Regimen"), width=3, status="info", solidHeader=TRUE,
            h5(strong("Amiodarone (Antiarrhythmic)")),
            checkboxInput("use_amio", "Use Amiodarone", value=TRUE),
            sliderInput("dose_amio_inp", "Maintenance Dose (mg/day)", 50, 400, 200, 50),
            helpText("Loading: 400 mg BID × 28 days automatically applied"),
            hr(),
            h5(strong("Apixaban (Anticoagulant)")),
            checkboxInput("use_apix", "Use Apixaban", value=TRUE),
            sliderInput("dose_apix_inp", "Dose (mg BID)", 2.5, 10, 5, 2.5),
            helpText("2.5 mg BID if ≥2 dose-reduction criteria"),
            hr(),
            h5(strong("Metoprolol (Rate Control)")),
            checkboxInput("use_metro", "Use Metoprolol", value=FALSE),
            sliderInput("dose_metro_inp", "Dose (mg BID)", 12.5, 200, 50, 12.5),
            hr(),
            sliderInput("pk_tmax", "Simulation Duration (days)", 30, 365, 180, 30),
            actionButton("run_pk", "Run PK Simulation", icon=icon("play"),
                         class="btn-primary btn-block")
          ),

          box(title=tagList(icon("chart-line"), " Plasma Concentration-Time Profiles"),
              width=9, status="primary", solidHeader=TRUE,
            tabsetPanel(
              tabPanel("Amiodarone",
                plotOutput("pk_amio_plot", height="300px"),
                div(style="padding:10px; background:#eef; border-radius:6px; font-size:12px;",
                  strong("PK Parameters:"), " ka=0.06/h, CL=3.5 L/h, V1=40 L, V2=4200 L (fat depot), F=46%",
                  br(), "t½α ~3.2d, t½β ~40-55 days (very long due to fat depot)",
                  br(), strong("Therapeutic range:"), " 1.0–2.5 µg/mL"
                )
              ),
              tabPanel("Apixaban",
                plotOutput("pk_apix_plot", height="300px"),
                div(style="padding:10px; background:#eef; border-radius:6px; font-size:12px;",
                  strong("PK Parameters:"), " ka=1.2/h, CL=3.3 L/h, V1=21 L, F=50%",
                  br(), "t½ ~12 hours; Cmax ~120 ng/mL (5mg BID)",
                  br(), strong("Therapeutic range:"), " 50–200 ng/mL (trough/peak)"
                )
              ),
              tabPanel("Metoprolol",
                plotOutput("pk_metro_plot", height="300px"),
                div(style="padding:10px; background:#eef; border-radius:6px; font-size:12px;",
                  strong("PK Parameters:"), " ka=1.5/h, CL=65 L/h, V1=290 L, F=40%",
                  br(), "t½ ~3.5 hours; high first-pass metabolism",
                  br(), strong("Therapeutic range:"), " 20–100 ng/mL"
                )
              ),
              tabPanel("PK Summary Table",
                tableOutput("pk_summary_table")
              )
            )
          )
        )
      ),

      # ======================================================================
      # TAB 3: ELECTROPHYSIOLOGY PD
      # ======================================================================
      tabItem(tabName = "tab_ep",
        fluidRow(
          box(title="EP Model Controls", width=3, status="info", solidHeader=TRUE,
            h5(strong("Patient & Drug Settings")),
            helpText("Uses settings from Tabs 1 & 2"),
            sliderInput("ep_tmax", "Duration (days)", 30, 365, 365, 30),
            hr(),
            h5("AF Type Zones:"),
            div(style="background:#E8F8E8;padding:8px;border-radius:4px;margin:3px;",
              strong("Paroxysmal:"), " AF burden <25%"),
            div(style="background:#FFF8E8;padding:8px;border-radius:4px;margin:3px;",
              strong("Persistent:"), " 25–75%"),
            div(style="background:#F8E8E8;padding:8px;border-radius:4px;margin:3px;",
              strong("Long-standing/Permanent:"), " >75%"),
            hr(),
            actionButton("run_ep", "Run EP Simulation", icon=icon("play"),
                         class="btn-primary btn-block")
          ),

          box(title=tagList(icon("wave-square"), " Electrophysiology Outputs"),
              width=9, status="primary", solidHeader=TRUE,
            tabsetPanel(
              tabPanel("AF Burden (%)",
                plotOutput("ep_af_plot", height="280px")
              ),
              tabPanel("ERP (ms)",
                plotOutput("ep_erp_plot", height="280px")
              ),
              tabPanel("Heart Rate Control",
                plotOutput("ep_hr_plot", height="280px")
              ),
              tabPanel("QTc Safety",
                plotOutput("ep_qtc_plot", height="280px")
              ),
              tabPanel("ERP vs AF Burden",
                plotOutput("ep_erp_af_scatter", height="280px")
              )
            )
          )
        )
      ),

      # ======================================================================
      # TAB 4: THROMBOEMBOLISM RISK
      # ======================================================================
      tabItem(tabName = "tab_thrombo",
        fluidRow(
          box(title="Coagulation Controls", width=3, status="info", solidHeader=TRUE,
            helpText("Anticoagulation settings from Tab 2"),
            sliderInput("thrombo_tmax", "Duration (days)", 30, 365, 365, 30),
            hr(),
            h5(strong("ARISTOTLE Trial Data:")),
            div(style="background:#EEF;padding:8px;border-radius:4px;font-size:12px;",
              "Apixaban 5mg BID vs Warfarin:", br(),
              "Stroke/SE: 1.27 vs 1.60%/yr", br(),
              strong("RRR: 21%, ARR: 0.33%/yr"), br(),
              "NNT: ~303 for 1 year"
            ),
            hr(),
            actionButton("run_thrombo", "Run Simulation", icon=icon("play"),
                         class="btn-primary btn-block")
          ),

          box(title=tagList(icon("tint"), " Thromboembolism Risk Outputs"),
              width=9, status="danger", solidHeader=TRUE,
            fluidRow(
              column(6, plotOutput("thrombo_fxa_plot",     height="240px")),
              column(6, plotOutput("thrombo_thr_plot",     height="240px"))
            ),
            fluidRow(
              column(6, plotOutput("thrombo_stroke_plot",  height="240px")),
              column(6, plotOutput("thrombo_pie",          height="240px"))
            ),
            fluidRow(
              box(width=12, status="warning",
                strong("NNT Calculator"),
                div(style="font-size:13px;",
                  uiOutput("nnt_output")
                )
              )
            )
          )
        )
      ),

      # ======================================================================
      # TAB 5: SCENARIO COMPARISON
      # ======================================================================
      tabItem(tabName = "tab_scenarios",
        fluidRow(
          box(title="Simulation Controls", width=3, status="info", solidHeader=TRUE,
            helpText("Run all 6 standard treatment scenarios"),
            radioButtons("highlight_scenario", "Highlight Scenario:",
              choices = c(
                "No Treatment"        = "S1",
                "Metoprolol"          = "S2",
                "Amiodarone"          = "S3",
                "Apixaban"            = "S4",
                "Metro + Apix"        = "S5",
                "Amio + Apix (SoC)"   = "S6"
              ), selected="S6"),
            sliderInput("scen_tmax", "Duration (days)", 90, 365, 365, 90),
            actionButton("run_all", "Run All Scenarios",
                         icon=icon("play-circle"),
                         class="btn-success btn-block"),
            hr(),
            downloadButton("dl_results", "Download CSV", class="btn-default btn-block")
          ),

          box(title=tagList(icon("chart-bar"), " Treatment Scenario Comparison"),
              width=9, status="primary", solidHeader=TRUE,
            tabsetPanel(
              tabPanel("AF Burden Over Time",
                plotOutput("scen_af_plot", height="300px")
              ),
              tabPanel("Stroke Risk Over Time",
                plotOutput("scen_stroke_plot", height="300px")
              ),
              tabPanel("Summary Table",
                div(style="overflow-x:auto;",
                  tableOutput("scen_summary_table")
                )
              ),
              tabPanel("Bar Comparison",
                fluidRow(
                  column(6, plotOutput("scen_bar_af",     height="280px")),
                  column(6, plotOutput("scen_bar_stroke", height="280px"))
                )
              )
            )
          )
        )
      ),

      # ======================================================================
      # TAB 6: BIOMARKER DASHBOARD
      # ======================================================================
      tabItem(tabName = "tab_biomarker",
        fluidRow(
          box(title="Biomarker Controls", width=3, status="info", solidHeader=TRUE,
            helpText("Uses current patient profile & drug settings"),
            sliderInput("bio_tmax", "Duration (days)", 30, 365, 365, 30),
            hr(),
            h5(strong("Normal Ranges:")),
            div(style="font-size:11px;",
              "NT-proBNP: <125 pg/mL", br(),
              "CRP: <3 mg/L", br(),
              "LA diameter: <40 mm", br(),
              "LVEF: ≥50%", br(),
              "Fibrosis: <15% (mild)"
            ),
            hr(),
            actionButton("run_bio", "Run Biomarker Simulation",
                         icon=icon("play"), class="btn-primary btn-block"),
            hr(),
            downloadButton("dl_report", "Download Report", class="btn-default btn-block")
          ),

          box(title=tagList(icon("microscope"), " Biomarker Dashboard"),
              width=9, status="success", solidHeader=TRUE,
            fluidRow(
              column(3, div(class="metric-box",
                div(class="metric-value", textOutput("bio_af_val")),
                div(class="metric-label", "AF Burden (%)"),
                uiOutput("traffic_af")
              )),
              column(3, div(class="metric-box",
                div(class="metric-value", textOutput("bio_erp_val")),
                div(class="metric-label", "ERP (ms)"),
                uiOutput("traffic_erp")
              )),
              column(3, div(class="metric-box",
                div(class="metric-value", textOutput("bio_hr_val")),
                div(class="metric-label", "Heart Rate (bpm)"),
                uiOutput("traffic_hr")
              )),
              column(3, div(class="metric-box",
                div(class="metric-value", textOutput("bio_stroke_val")),
                div(class="metric-label", "Stroke Risk (%/yr)"),
                uiOutput("traffic_stroke")
              ))
            ),
            hr(),
            tabsetPanel(
              tabPanel("NT-proBNP & CRP",
                plotOutput("bio_bnp_crp_plot", height="260px")
              ),
              tabPanel("Fibrosis & LA Diameter",
                plotOutput("bio_fib_la_plot", height="260px")
              ),
              tabPanel("Neurohormonal (AngII, NE, ROS)",
                plotOutput("bio_neuro_plot", height="260px")
              ),
              tabPanel("Traffic Light Summary",
                plotOutput("bio_traffic_plot", height="260px")
              )
            )
          )
        )
      )
    )  # end tabItems
  )    # end dashboardBody
)      # end dashboardPage

# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {

  # ---------- Reactive: CHA2DS2-VASc score ----------
  chads_score <- reactive({
    calc_chads(
      age          = input$age,
      htn          = input$htn,
      dm           = input$dm,
      hf           = input$hf,
      prior_stroke = input$prior_stroke,
      pad          = input$pad,
      female       = (input$sex == "female")
    )
  })

  has_bled_score <- reactive({
    s <- 0
    if (input$hasbled_htn)    s <- s + 1
    if (input$hasbled_renal)  s <- s + 1
    if (input$hasbled_bleed)  s <- s + 1
    if (input$hasbled_alcohol)s <- s + 1
    if (input$age > 65)       s <- s + 1
    s
  })

  output$chads_score_out <- renderText({
    chads_score()
  })

  output$chads_risk_out <- renderText({
    risk <- chads_stroke_risk(chads_score())
    paste0(risk, "%/year")
  })

  output$hasbled_out <- renderText({
    has_bled_score()
  })

  output$anticoag_rec <- renderUI({
    sc <- chads_score()
    if (sc >= 2) {
      div(style="background:#E8F8E8;padding:10px;border-radius:6px;",
        icon("check-circle", style="color:green"),
        strong(" RECOMMEND anticoagulation"),
        p("CHA₂DS₂-VASc ≥2: Oral anticoagulation indicated."),
        p("Preferred: ", strong("Apixaban 5mg BID"), " (ARISTOTLE: RRR 21% vs warfarin)")
      )
    } else if (sc == 1) {
      div(style="background:#FFF8E8;padding:10px;border-radius:6px;",
        icon("exclamation-circle", style="color:orange"),
        strong(" CONSIDER anticoagulation"),
        p("CHA₂DS₂-VASc = 1: Anticoagulation may be considered.")
      )
    } else {
      div(style="background:#F8E8E8;padding:10px;border-radius:6px;",
        icon("info-circle", style="color:blue"),
        strong(" No anticoagulation needed"),
        p("CHA₂DS₂-VASc = 0 (male) or 1 (female): Low risk.")
      )
    }
  })

  output$renal_guidance <- renderUI({
    egfr <- input$egfr
    if (egfr < 15) {
      div(style="color:red;", strong("eGFR <15: DOACs contraindicated. Consider warfarin."))
    } else if (egfr < 30) {
      div(style="color:orange;", strong("eGFR 15-29: Reduce dabigatran. Avoid rivaroxaban. Apixaban 2.5mg BID if criteria."))
    } else if (egfr < 50) {
      div(style="color:#aa8800;", strong("eGFR 30-49: Consider apixaban 2.5mg BID if age ≥80 or weight ≤60kg."))
    } else {
      div(style="color:green;", strong("eGFR ≥50: Standard dosing appropriate."))
    }
  })

  # ---- PK Simulations ----
  pk_data <- eventReactive(input$run_pk, {
    dose_amio  <- if (input$use_amio)  input$dose_amio_inp  else 0
    dose_apix  <- if (input$use_apix)  input$dose_apix_inp  else 0
    dose_metro <- if (input$use_metro) input$dose_metro_inp else 0

    run_scenario_sim(
      af_type    = input$af_type,
      dose_amio  = dose_amio,
      dose_apix  = dose_apix,
      dose_metro = dose_metro,
      t_days     = input$pk_tmax,
      age        = input$age,
      chads       = chads_score()
    )
  }, ignoreNULL=FALSE)

  output$pk_amio_plot <- renderPlot({
    df <- pk_data()
    if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=Cp_AMIO_out)) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=1.0, ymax=2.5,
               alpha=0.2, fill="green") +
      geom_line(color="#1A5276", size=1.2) +
      annotate("text", x=input$pk_tmax*0.8, y=1.75,
               label="Therapeutic\n1–2.5 µg/mL", size=3, color="darkgreen") +
      labs(title="Amiodarone Plasma Concentration",
           x="Time (days)", y="Concentration (µg/mL)") +
      theme_bw(base_size=12) +
      theme(panel.grid.minor=element_blank())
  })

  output$pk_apix_plot <- renderPlot({
    df <- pk_data()
    if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=Cp_APIX_out)) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=50, ymax=200,
               alpha=0.2, fill="blue") +
      geom_line(color="#7D3C98", size=1.2) +
      annotate("text", x=input$pk_tmax*0.8, y=125,
               label="Therapeutic\n50–200 ng/mL", size=3, color="navy") +
      labs(title="Apixaban Plasma Concentration",
           x="Time (days)", y="Concentration (ng/mL)") +
      theme_bw(base_size=12) +
      theme(panel.grid.minor=element_blank())
  })

  output$pk_metro_plot <- renderPlot({
    df <- pk_data()
    if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=Cp_METRO_out)) +
      annotate("rect", xmin=-Inf, xmax=Inf, ymin=20, ymax=100,
               alpha=0.2, fill="orange") +
      geom_line(color="#E67E22", size=1.2) +
      annotate("text", x=input$pk_tmax*0.8, y=60,
               label="Therapeutic\n20–100 ng/mL", size=3, color="darkorange") +
      labs(title="Metoprolol Plasma Concentration",
           x="Time (days)", y="Concentration (ng/mL)") +
      theme_bw(base_size=12) +
      theme(panel.grid.minor=element_blank())
  })

  output$pk_summary_table <- renderTable({
    df <- pk_data()
    if (is.null(df)) return(NULL)
    last <- tail(df, 5)
    data.frame(
      "Drug"          = c("Amiodarone","Apixaban","Metoprolol"),
      "Final Cp"      = c(
        round(last$Cp_AMIO_out[nrow(last)], 3),
        round(last$Cp_APIX_out[nrow(last)], 1),
        round(last$Cp_METRO_out[nrow(last)], 1)
      ),
      "Units"         = c("µg/mL","ng/mL","ng/mL"),
      "Ther.Range"    = c("1.0–2.5","50–200","20–100"),
      "In Range"      = c(
        ifelse(last$Cp_AMIO_out[nrow(last)] >= 1.0 & last$Cp_AMIO_out[nrow(last)] <= 2.5, "Yes","No"),
        ifelse(last$Cp_APIX_out[nrow(last)] >= 50 & last$Cp_APIX_out[nrow(last)] <= 200, "Yes","No"),
        ifelse(last$Cp_METRO_out[nrow(last)] >= 20 & last$Cp_METRO_out[nrow(last)] <= 100, "Yes","No")
      )
    )
  })

  # ---- EP Simulation ----
  ep_data <- eventReactive(input$run_ep, {
    dose_amio  <- if (input$use_amio)  input$dose_amio_inp  else 0
    dose_apix  <- if (input$use_apix)  input$dose_apix_inp  else 0
    dose_metro <- if (input$use_metro) input$dose_metro_inp else 0
    run_scenario_sim(
      af_type   = input$af_type, dose_amio=dose_amio,
      dose_apix=dose_apix, dose_metro=dose_metro,
      t_days=input$ep_tmax, age=input$age, chads=chads_score()
    )
  }, ignoreNULL=FALSE)

  output$ep_af_plot <- renderPlot({
    df <- ep_data()
    if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=AF_BURDEN*100)) +
      annotate("rect", xmin=-Inf,xmax=Inf, ymin=75,ymax=100, alpha=0.12, fill="red") +
      annotate("rect", xmin=-Inf,xmax=Inf, ymin=25,ymax=75,  alpha=0.12, fill="orange") +
      annotate("rect", xmin=-Inf,xmax=Inf, ymin=0, ymax=25,  alpha=0.12, fill="green") +
      annotate("text", x=5, y=88, label="Long-standing Persistent", size=3, hjust=0, color="#922B21") +
      annotate("text", x=5, y=50, label="Persistent", size=3, hjust=0, color="#D35400") +
      annotate("text", x=5, y=12, label="Paroxysmal", size=3, hjust=0, color="#1E8449") +
      geom_line(color="#1A5276", size=1.3) +
      scale_y_continuous(limits=c(0,100)) +
      labs(title="AF Burden Over Time", x="Time (days)", y="AF Burden (%)") +
      theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
  })

  output$ep_erp_plot <- renderPlot({
    df <- ep_data()
    if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=ERP)) +
      geom_hline(yintercept=200, linetype="dashed", color="blue", size=0.8) +
      geom_line(color="#1E8449", size=1.3) +
      annotate("text", x=input$ep_tmax*0.7, y=202,
               label="Reentry threshold (200ms)", size=3, color="blue") +
      labs(title="Atrial ERP Over Time", x="Time (days)", y="ERP (ms)") +
      theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
  })

  output$ep_hr_plot <- renderPlot({
    df <- ep_data()
    if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=HR_AF)) +
      annotate("rect", xmin=-Inf,xmax=Inf, ymin=60,ymax=110, alpha=0.15, fill="green") +
      geom_line(color="#E74C3C", size=1.3) +
      annotate("text", x=5, y=85, label="Rate control target (<110 bpm)", size=3, hjust=0, color="darkgreen") +
      labs(title="Ventricular Rate During AF", x="Time (days)", y="Heart Rate (bpm)") +
      theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
  })

  output$ep_qtc_plot <- renderPlot({
    df <- ep_data()
    if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=QTc)) +
      annotate("rect", xmin=-Inf,xmax=Inf, ymin=470,ymax=Inf, alpha=0.15, fill="red") +
      annotate("rect", xmin=-Inf,xmax=Inf, ymin=440,ymax=470, alpha=0.1, fill="orange") +
      geom_line(color="#8E44AD", size=1.3) +
      annotate("text", x=5, y=475, label="QTc > 470ms: High TdP Risk", size=3, hjust=0, color="red") +
      annotate("text", x=5, y=452, label="QTc 440–470ms: Monitor", size=3, hjust=0, color="orange") +
      labs(title="QTc Interval (Safety Monitoring)", x="Time (days)", y="QTc (ms)") +
      theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
  })

  output$ep_erp_af_scatter <- renderPlot({
    df <- ep_data()
    if (is.null(df)) return(NULL)
    df_sub <- df[seq(1, nrow(df), by=4), ]
    ggplot(df_sub, aes(x=ERP, y=AF_BURDEN*100, color=time_days)) +
      geom_point(size=1.5, alpha=0.7) +
      scale_color_gradient(low="red", high="blue", name="Day") +
      geom_vline(xintercept=200, linetype="dashed", color="gray50") +
      labs(title="ERP vs AF Burden (Phase Plot)",
           x="Atrial ERP (ms)", y="AF Burden (%)") +
      theme_bw(base_size=12) + theme(panel.grid.minor=element_blank())
  })

  # ---- Thromboembolism ----
  thrombo_data <- eventReactive(input$run_thrombo, {
    dose_apix  <- if (input$use_apix)  input$dose_apix_inp  else 0
    dose_amio  <- if (input$use_amio)  input$dose_amio_inp  else 0
    dose_metro <- if (input$use_metro) input$dose_metro_inp else 0
    run_scenario_sim(
      af_type=input$af_type, dose_amio=dose_amio, dose_apix=dose_apix,
      dose_metro=dose_metro, t_days=input$thrombo_tmax,
      age=input$age, chads=chads_score()
    )
  }, ignoreNULL=FALSE)

  output$thrombo_fxa_plot <- renderPlot({
    df <- thrombo_data()
    if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=FXa)) +
      geom_line(color="#C0392B", size=1.2) +
      geom_hline(yintercept=1, linetype="dashed", color="gray") +
      labs(title="FXa Activity", x="Days", y="FXa (rel. units)") +
      theme_bw(base_size=11) + theme(panel.grid.minor=element_blank())
  })

  output$thrombo_thr_plot <- renderPlot({
    df <- thrombo_data()
    if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=Thrombin)) +
      geom_line(color="#922B21", size=1.2) +
      labs(title="Thrombin Activity", x="Days", y="Thrombin (rel. units)") +
      theme_bw(base_size=11) + theme(panel.grid.minor=element_blank())
  })

  output$thrombo_stroke_plot <- renderPlot({
    df <- thrombo_data()
    if (is.null(df)) return(NULL)
    ggplot(df, aes(x=time_days, y=STROKE_RISK)) +
      geom_line(color="#17202A", size=1.2) +
      labs(title="Annual Stroke Risk", x="Days", y="Stroke Risk (%/yr)") +
      theme_bw(base_size=11) + theme(panel.grid.minor=element_blank())
  })

  output$thrombo_pie <- renderPlot({
    df <- thrombo_data()
    if (is.null(df)) return(NULL)
    last_row <- tail(df, 1)
    base_risk <- last_row$STROKE_RISK[1]
    rrr <- 0.21  # ARISTOTLE
    treated_risk <- base_risk * (1 - rrr)
    pie_data <- data.frame(
      cat   = c("Residual stroke risk\n(with anticoag)",
                "Absolute risk reduced\n(by anticoag)"),
      value = c(treated_risk, base_risk - treated_risk),
      fill  = c("#E74C3C","#27AE60")
    )
    ggplot(pie_data, aes(x="", y=value, fill=cat)) +
      geom_col(width=1) +
      coord_polar("y") +
      scale_fill_manual(values=c("#E74C3C","#27AE60")) +
      labs(title="Absolute Risk Reduction\n(Apixaban, ARISTOTLE RRR=21%)",
           fill=NULL, x=NULL, y=NULL) +
      theme_void(base_size=10) +
      theme(legend.position="bottom")
  })

  output$nnt_output <- renderUI({
    df <- thrombo_data()
    if (is.null(df)) return(NULL)
    last_row <- tail(df, 1)
    base_risk <- last_row$STROKE_RISK[1] / 100
    arr <- base_risk * 0.21
    nnt <- if (arr > 0) round(1 / arr) else "N/A"
    div(
      strong("Based on ARISTOTLE trial (apixaban 5mg BID):"),
      tags$ul(
        tags$li(paste("Baseline stroke risk: ", round(base_risk*100,2), "%/year")),
        tags$li(paste("Relative risk reduction (RRR): 21%")),
        tags$li(paste("Absolute risk reduction (ARR):", round(arr*100,3), "%/year")),
        tags$li(paste0("NNT for 1 year: ", nnt, " patients"))
      ),
      em("Reference: Granger et al. NEJM 2011;365:981. PMID: 21870978")
    )
  })

  # ---- All 6 Scenarios ----
  scen_data <- eventReactive(input$run_all, {
    withProgress(message="Running 6 scenarios...", value=0, {
      sc_list <- list()
      sc_defs <- list(
        S1=list(label="No Treatment",      amio=0,   apix=0, metro=0),
        S2=list(label="Metoprolol",        amio=0,   apix=0, metro=50),
        S3=list(label="Amiodarone",        amio=200, apix=0, metro=0),
        S4=list(label="Apixaban",          amio=0,   apix=5, metro=0),
        S5=list(label="Metro+Apix",        amio=0,   apix=5, metro=50),
        S6=list(label="Amio+Apix (SoC)",   amio=200, apix=5, metro=0)
      )
      n <- length(sc_defs)
      for (i in seq_along(sc_defs)) {
        sc  <- sc_defs[[i]]
        sid <- names(sc_defs)[i]
        incProgress(1/n, detail=paste("Running", sc$label))
        df <- run_scenario_sim(
          af_type   = input$af_type,
          dose_amio = sc$amio, dose_apix = sc$apix, dose_metro = sc$metro,
          t_days    = input$scen_tmax,
          age       = input$age, chads = chads_score()
        )
        if (!is.null(df)) {
          df$scenario_id <- sid
          df$label       <- sc$label
          sc_list[[sid]] <- df
        }
      }
      bind_rows(sc_list)
    })
  })

  scen_colors <- c(
    "No Treatment"   = "#E41A1C",
    "Metoprolol"     = "#377EB8",
    "Amiodarone"     = "#4DAF4A",
    "Apixaban"       = "#984EA3",
    "Metro+Apix"     = "#FF7F00",
    "Amio+Apix (SoC)"= "#A65628"
  )

  output$scen_af_plot <- renderPlot({
    df <- scen_data()
    if (is.null(df)) return(NULL)
    df$label <- factor(df$label, levels=names(scen_colors))
    ggplot(df, aes(x=time_days, y=AF_BURDEN*100, color=label,
                   size=ifelse(scenario_id==input$highlight_scenario, 1.8, 0.8))) +
      geom_line(alpha=0.85) +
      scale_size_identity() +
      scale_color_manual(values=scen_colors, name="Scenario") +
      labs(title="AF Burden — All Scenarios", x="Time (days)", y="AF Burden (%)") +
      theme_bw(base_size=12) + theme(legend.position="bottom", panel.grid.minor=element_blank())
  })

  output$scen_stroke_plot <- renderPlot({
    df <- scen_data()
    if (is.null(df)) return(NULL)
    df$label <- factor(df$label, levels=names(scen_colors))
    ggplot(df, aes(x=time_days, y=STROKE_RISK, color=label,
                   size=ifelse(scenario_id==input$highlight_scenario, 1.8, 0.8))) +
      geom_line(alpha=0.85) +
      scale_size_identity() +
      scale_color_manual(values=scen_colors, name="Scenario") +
      labs(title="Annual Stroke Risk — All Scenarios", x="Time (days)", y="Stroke Risk (%/yr)") +
      theme_bw(base_size=12) + theme(legend.position="bottom", panel.grid.minor=element_blank())
  })

  output$scen_summary_table <- renderTable({
    df <- scen_data()
    if (is.null(df)) return(NULL)
    df %>%
      filter(time_days > max(time_days, na.rm=TRUE) - 1) %>%
      group_by(Scenario=label) %>%
      summarise(
        "AF Burden (%)"     = round(mean(AF_BURDEN*100, na.rm=TRUE), 1),
        "ERP (ms)"          = round(mean(ERP, na.rm=TRUE), 1),
        "HR (bpm)"          = round(mean(HR_AF, na.rm=TRUE), 1),
        "QTc (ms)"          = round(mean(QTc, na.rm=TRUE), 1),
        "Stroke Risk (%/yr)"= round(mean(STROKE_RISK, na.rm=TRUE), 2),
        "Fibrosis (%)"      = round(mean(Fibrosis*100, na.rm=TRUE), 1),
        "NT-proBNP (pg/mL)" = round(mean(NT_proBNP, na.rm=TRUE), 0),
        .groups="drop"
      )
  }, striped=TRUE, bordered=TRUE, hover=TRUE, spacing="s")

  output$scen_bar_af <- renderPlot({
    df <- scen_data()
    if (is.null(df)) return(NULL)
    sumdf <- df %>%
      filter(time_days > max(time_days,na.rm=TRUE)-1) %>%
      group_by(label) %>%
      summarise(af=round(mean(AF_BURDEN*100,na.rm=TRUE),1), .groups="drop")
    sumdf$label <- factor(sumdf$label, levels=names(scen_colors))
    ggplot(sumdf, aes(x=label, y=af, fill=label)) +
      geom_col(alpha=0.85) +
      geom_text(aes(label=paste0(af,"%")), vjust=-0.3, size=3) +
      scale_fill_manual(values=scen_colors) +
      labs(title="AF Burden at End of Simulation",
           x=NULL, y="AF Burden (%)") +
      theme_bw(base_size=11) +
      theme(axis.text.x=element_text(angle=30,hjust=1), legend.position="none",
            panel.grid.minor=element_blank())
  })

  output$scen_bar_stroke <- renderPlot({
    df <- scen_data()
    if (is.null(df)) return(NULL)
    sumdf <- df %>%
      filter(time_days > max(time_days,na.rm=TRUE)-1) %>%
      group_by(label) %>%
      summarise(sr=round(mean(STROKE_RISK,na.rm=TRUE),2), .groups="drop")
    sumdf$label <- factor(sumdf$label, levels=names(scen_colors))
    ggplot(sumdf, aes(x=label, y=sr, fill=label)) +
      geom_col(alpha=0.85) +
      geom_text(aes(label=paste0(sr,"%")), vjust=-0.3, size=3) +
      scale_fill_manual(values=scen_colors) +
      labs(title="Stroke Risk at End of Simulation",
           x=NULL, y="Stroke Risk (%/year)") +
      theme_bw(base_size=11) +
      theme(axis.text.x=element_text(angle=30,hjust=1), legend.position="none",
            panel.grid.minor=element_blank())
  })

  output$dl_results <- downloadHandler(
    filename = function() paste0("AF_QSP_scenarios_", Sys.Date(), ".csv"),
    content  = function(file) {
      df <- scen_data()
      if (!is.null(df)) write.csv(df, file, row.names=FALSE)
    }
  )

  # ---- Biomarker Dashboard ----
  bio_data <- eventReactive(input$run_bio, {
    dose_amio  <- if (input$use_amio)  input$dose_amio_inp  else 0
    dose_apix  <- if (input$use_apix)  input$dose_apix_inp  else 0
    dose_metro <- if (input$use_metro) input$dose_metro_inp else 0
    run_scenario_sim(
      af_type=input$af_type, dose_amio=dose_amio, dose_apix=dose_apix,
      dose_metro=dose_metro, t_days=input$bio_tmax,
      age=input$age, chads=chads_score()
    )
  }, ignoreNULL=FALSE)

  bio_last <- reactive({
    df <- bio_data()
    if (is.null(df)) return(NULL)
    tail(df, 1)
  })

  output$bio_af_val    <- renderText({ if(!is.null(bio_last())) round(bio_last()$AF_BURDEN*100,1) else "—" })
  output$bio_erp_val   <- renderText({ if(!is.null(bio_last())) round(bio_last()$ERP,0) else "—" })
  output$bio_hr_val    <- renderText({ if(!is.null(bio_last())) round(bio_last()$HR_AF,0) else "—" })
  output$bio_stroke_val<- renderText({ if(!is.null(bio_last())) round(bio_last()$STROKE_RISK,2) else "—" })

  traffic_light <- function(val, green_max, yellow_max, label_desc) {
    if (is.null(val)) return(NULL)
    cls <- if (val <= green_max) "traffic-green" else if (val <= yellow_max) "traffic-yellow" else "traffic-red"
    stat <- if (val <= green_max) "✓ Normal" else if (val <= yellow_max) "⚠ Borderline" else "✗ High"
    div(class=cls, stat)
  }

  output$traffic_af     <- renderUI({ if(!is.null(bio_last())) traffic_light(bio_last()$AF_BURDEN*100, 25, 75) })
  output$traffic_erp    <- renderUI({ if(!is.null(bio_last())) traffic_light(bio_last()$HR_AF, 80, 110) })
  output$traffic_hr     <- renderUI({ if(!is.null(bio_last())) traffic_light(bio_last()$HR_AF, 80, 110) })
  output$traffic_stroke <- renderUI({ if(!is.null(bio_last())) traffic_light(bio_last()$STROKE_RISK, 2, 4) })

  output$bio_bnp_crp_plot <- renderPlot({
    df <- bio_data()
    if (is.null(df)) return(NULL)
    df_long <- df %>%
      select(time_days, NT_proBNP, CRP_out) %>%
      pivot_longer(-time_days, names_to="marker", values_to="value")
    ggplot(df_long, aes(x=time_days, y=value, color=marker)) +
      geom_line(size=1.2) +
      facet_wrap(~marker, scales="free_y",
                 labeller=labeller(marker=c(NT_proBNP="NT-proBNP (pg/mL)", CRP_out="CRP (mg/L)"))) +
      scale_color_manual(values=c(NT_proBNP="#E74C3C", CRP_out="#E67E22")) +
      labs(title="Cardiac & Inflammatory Biomarkers", x="Time (days)", y="Concentration") +
      theme_bw(base_size=11) + theme(legend.position="none", panel.grid.minor=element_blank())
  })

  output$bio_fib_la_plot <- renderPlot({
    df <- bio_data()
    if (is.null(df)) return(NULL)
    df_long <- df %>%
      select(time_days, Fibrosis, LA_diam) %>%
      mutate(Fibrosis_pct = Fibrosis*100) %>%
      select(time_days, Fibrosis_pct, LA_diam) %>%
      pivot_longer(-time_days, names_to="marker", values_to="value")
    ggplot(df_long, aes(x=time_days, y=value, color=marker)) +
      geom_line(size=1.2) +
      facet_wrap(~marker, scales="free_y",
                 labeller=labeller(marker=c(Fibrosis_pct="Fibrosis Score (%)", LA_diam="LA Diameter (mm)"))) +
      scale_color_manual(values=c(Fibrosis_pct="#1E8449", LA_diam="#2874A6")) +
      labs(title="Structural Remodeling Biomarkers", x="Time (days)", y="Value") +
      theme_bw(base_size=11) + theme(legend.position="none", panel.grid.minor=element_blank())
  })

  output$bio_neuro_plot <- renderPlot({
    df <- bio_data()
    if (is.null(df)) return(NULL)
    df_long <- df %>%
      select(time_days, AngII, ROS, NE) %>%
      pivot_longer(-time_days, names_to="marker", values_to="value")
    ggplot(df_long, aes(x=time_days, y=value, color=marker)) +
      geom_line(size=1.2) +
      scale_color_manual(values=c(AngII="#E74C3C", ROS="#8E44AD", NE="#F39C12"),
                         labels=c(AngII="Angiotensin II", ROS="Reactive Oxygen Species", NE="Norepinephrine")) +
      labs(title="Neurohormonal Mediators (Relative Units)",
           x="Time (days)", y="Relative Units", color="Mediator") +
      geom_hline(yintercept=1, linetype="dashed", color="gray50") +
      theme_bw(base_size=11) + theme(panel.grid.minor=element_blank(), legend.position="bottom")
  })

  output$bio_traffic_plot <- renderPlot({
    df <- bio_data()
    if (is.null(df)) return(NULL)
    last <- tail(df, 1)

    metrics <- data.frame(
      metric = c("AF Burden", "ERP (ms)", "Heart Rate", "QTc (ms)",
                 "Stroke Risk", "Fibrosis", "NT-proBNP"),
      value  = c(
        round(last$AF_BURDEN*100, 1),
        round(last$ERP, 0),
        round(last$HR_AF, 0),
        round(last$QTc, 0),
        round(last$STROKE_RISK, 2),
        round(last$Fibrosis*100, 1),
        round(last$NT_proBNP, 0)
      ),
      green_max  = c(25, Inf, 80, 440, 1.5, 15, 125),
      yellow_max = c(75, 190, 110, 470, 4.0, 35, 500)
    )
    metrics$status <- with(metrics,
      ifelse(value <= green_max, "Normal",
             ifelse(value <= yellow_max, "Borderline", "Abnormal"))
    )
    metrics$status <- factor(metrics$status, levels=c("Normal","Borderline","Abnormal"))

    ggplot(metrics, aes(x=reorder(metric, as.numeric(status)),
                        y=0.5, fill=status)) +
      geom_tile(height=0.9, color="white", size=2) +
      geom_text(aes(label=paste0(metric, "\n", value)), size=3.5, fontface="bold") +
      scale_fill_manual(values=c(Normal="#27AE60",Borderline="#F39C12",Abnormal="#E74C3C")) +
      labs(title="Biomarker Traffic Light Status (End of Simulation)",
           x=NULL, y=NULL, fill="Status") +
      coord_flip() +
      theme_minimal(base_size=12) +
      theme(axis.text=element_blank(), axis.ticks=element_blank(),
            panel.grid=element_blank(), legend.position="right")
  })

  output$dl_report <- downloadHandler(
    filename = function() paste0("AF_QSP_biomarker_report_", Sys.Date(), ".csv"),
    content  = function(file) {
      df <- bio_data()
      if (!is.null(df)) {
        out <- df %>%
          select(time_days, AF_BURDEN, ERP, QTc, HR_AF, STROKE_RISK,
                 Fibrosis, NT_proBNP, CRP_out, LA_diam, AngII, ROS, NE,
                 Cp_AMIO_out, Cp_APIX_out, Cp_METRO_out)
        write.csv(out, file, row.names=FALSE)
      }
    }
  )

}  # end server

# Run the app
shinyApp(ui = ui, server = server)
