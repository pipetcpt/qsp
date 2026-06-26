## =============================================================================
## NAFLD/MASLD QSP Interactive Shiny Dashboard
## 6 Tabs: Patient Profile · PK · PD Biomarkers · Clinical Endpoints ·
##         Scenario Comparison · Biomarker Panel
## =============================================================================

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(ggplot2)
library(dplyr)
library(tidyr)
library(plotly)
library(DT)
library(mrgsolve)

## ─────────────────────────────────────────────────────────────────────────────
## MODEL CODE (embedded for self-contained app)
## ─────────────────────────────────────────────────────────────────────────────

MODEL_CODE <- '
$PARAM
ka=0.8, CL=4.5, V1=15.0, V2=30.0, Q=2.0, F1=0.65
ka_glp1=0.005, CL_glp1=0.055, V1_glp1=8.5
kLFFA_in=0.18, kDNL_base=0.12, kBeta_ox=0.22, kVLDL_sec=0.10, kLTAG_deg=0.05
IR_base=1.0, k_IR_FFA=0.15, k_Adipo=0.30, Adipo_base=12.0, Leptin_base=18.0
k_Leptin_IR=0.04
kROS_gen=0.08, kROS_deg=0.25, kROS_FFA=0.10
kNRF2_act=0.05, kNRF2_deg=0.15
kER_FFA=0.06, kER_ROS=0.04, kER_deg=0.12
kKup_act=0.10, kKup_res=0.08
kTNF_prod=0.20, kTNF_deg=0.35
kIL6_prod=0.15, kIL6_deg=0.30
kIL1_prod=0.12, kIL1_deg=0.28
kNeutro=0.05, kMCP1_deg=0.25
kApop_TNF=0.03, kApop_ROS=0.02, kApop_ER=0.025, kApop_res=0.10
kALT_lysis=15.0, ALT_base=35.0
kHSC_act=0.006, kHSC_res=0.004
kTGFb_prod=0.08, kTGFb_deg=0.20
kCol_prod=0.012, kCol_deg=0.003
EC50_FXR=0.5, Emax_FXR_DNL=0.40, Emax_FXR_TGF=0.35, Emax_FXR_LPS=0.25
EC50_GLP1=0.3, Emax_GLP1_IR=0.45, Emax_GLP1_DNL=0.30, Emax_GLP1_Kup=0.35
kFGF19_prod=0.08, kFGF19_deg=0.20, kBA_synth=0.10

$CMT
GUT CENT PERI GUT_GLP1 CENT_GLP1
LFFA LTAG LDAG
ROS_LVR NRF2_ACT ER_STRESS
KUP_ACT TNF IL6 IL1B MCP1 NEUTRO
HEPATO_APOP
HSC_ACTIV TGF_B1 COLLAGEN

$MAIN
double Cp_FXR  = CENT / V1;
double Cp_GLP1 = CENT_GLP1 / V1_glp1;
double IR_index = IR_base
                  * (1 + k_IR_FFA * LFFA / 100.0)
                  * (1 + k_Leptin_IR * Leptin_base / 18.0)
                  * (1 / (1 + k_Adipo * Adipo_base / 12.0));
double E_FXR_DNL  = Emax_FXR_DNL * Cp_FXR / (EC50_FXR + Cp_FXR);
double E_FXR_TGF  = Emax_FXR_TGF * Cp_FXR / (EC50_FXR + Cp_FXR);
double E_FXR_LPS  = Emax_FXR_LPS * Cp_FXR / (EC50_FXR + Cp_FXR);
double E_GLP1_IR  = Emax_GLP1_IR  * Cp_GLP1 / (EC50_GLP1 + Cp_GLP1);
double E_GLP1_DNL = Emax_GLP1_DNL * Cp_GLP1 / (EC50_GLP1 + Cp_GLP1);
double E_GLP1_Kup = Emax_GLP1_Kup * Cp_GLP1 / (EC50_GLP1 + Cp_GLP1);
double kLFFA_eff = kLFFA_in * IR_index * (1 - E_GLP1_IR * 0.4);
double kDNL_eff  = kDNL_base * IR_index * (1 - E_FXR_DNL) * (1 - E_GLP1_DNL);
double kBeta_eff = kBeta_ox * (Adipo_base / 12.0) * (1 + 0.4 * NRF2_ACT);
double KUP_drive = kKup_act * (1 + 0.5 * ROS_LVR) * (1 - E_GLP1_Kup)
                 * (1 - E_FXR_LPS * 0.3);
double TGFb_drive = kTGFb_prod * KUP_ACT * (1 + 0.3 * IL1B)
                  * (1 - E_FXR_TGF);
double Col_syn = kCol_prod * HSC_ACTIV * TGF_B1;
double Col_deg_rate = kCol_deg * (1 - 0.6 * HSC_ACTIV);

$ODE
dxdt_GUT      = -ka * GUT;
dxdt_CENT     = ka * F1 * GUT - (CL/V1)*CENT - (Q/V1)*CENT + (Q/V2)*PERI;
dxdt_PERI     = (Q/V1)*CENT - (Q/V2)*PERI;
dxdt_GUT_GLP1 = -ka_glp1 * GUT_GLP1;
dxdt_CENT_GLP1 = ka_glp1 * GUT_GLP1 - (CL_glp1/V1_glp1)*CENT_GLP1;

dxdt_LFFA = kLFFA_eff + kDNL_eff * 80.0
           - kBeta_eff * LFFA - 0.15 * LFFA - 0.02 * LFFA;
dxdt_LTAG = 0.15 * LFFA - kVLDL_sec * LTAG - kLTAG_deg * LTAG;
dxdt_LDAG = 0.02 * LFFA - 0.08 * LDAG;

dxdt_ROS_LVR = kROS_gen + kROS_FFA * LFFA / 100.0 + 0.02 * KUP_ACT
              - kROS_deg * (1 + NRF2_ACT) * ROS_LVR;
dxdt_NRF2_ACT = kNRF2_act * ROS_LVR * (1 - NRF2_ACT) - kNRF2_deg * NRF2_ACT;
dxdt_ER_STRESS = kER_FFA * LFFA / 100.0 + kER_ROS * ROS_LVR - kER_deg * ER_STRESS;

dxdt_KUP_ACT = KUP_drive * (1 - KUP_ACT)
              - kKup_res * (1 + k_Adipo * Adipo_base / 12.0) * KUP_ACT;
dxdt_TNF  = kTNF_prod * KUP_ACT * (1 + 0.2 * IL1B) - kTNF_deg * TNF;
dxdt_IL6  = kIL6_prod * KUP_ACT + 0.05 * TNF - kIL6_deg * IL6;
dxdt_IL1B = kIL1_prod * KUP_ACT * (1 + 0.3 * ROS_LVR) - kIL1_deg * IL1B;
dxdt_MCP1 = 0.15 * KUP_ACT * (1 + TNF * 0.2) - kMCP1_deg * MCP1;
dxdt_NEUTRO = kNeutro * (IL1B + 0.5 * TNF) * (1 - NEUTRO * 0.3) - 0.08 * NEUTRO;

dxdt_HEPATO_APOP = kApop_TNF * TNF + kApop_ROS * ROS_LVR + kApop_ER * ER_STRESS
                  - kApop_res * HEPATO_APOP;

dxdt_HSC_ACTIV = kHSC_act * TGF_B1 * (1 - HSC_ACTIV)
               + 0.002 * HEPATO_APOP * (1 - HSC_ACTIV)
               + 0.001 * TNF * (1 - HSC_ACTIV)
               - kHSC_res * (1 + k_Adipo * Adipo_base / 12.0) * HSC_ACTIV;
dxdt_TGF_B1  = TGFb_drive + 0.05 * HSC_ACTIV - kTGFb_deg * TGF_B1;
dxdt_COLLAGEN = Col_syn - Col_deg_rate * COLLAGEN;

$TABLE
capture Cp_FXR_ug = CENT / V1;
capture Cp_GLP1_ug = CENT_GLP1 / V1_glp1;
capture Hepatic_TG = LTAG;
capture ROS_level = ROS_LVR;
capture ER_stress = ER_STRESS;
capture Kupffer_act = KUP_ACT;
capture TNFalpha = TNF;
capture IL6_level = IL6;
capture IL1beta = IL1B;
capture Apoptosis = HEPATO_APOP;
capture HSC_act = HSC_ACTIV;
capture TGFbeta1 = TGF_B1;
capture Collagen = COLLAGEN;
capture ALT = ALT_base + kALT_lysis * HEPATO_APOP;

$INIT
GUT=0, CENT=0, PERI=0, GUT_GLP1=0, CENT_GLP1=0
LFFA=120, LTAG=80, LDAG=15
ROS_LVR=0.6, NRF2_ACT=0.25, ER_STRESS=0.4
KUP_ACT=0.3, TNF=0.25, IL6=0.2, IL1B=0.15, MCP1=0.1, NEUTRO=0.1
HEPATO_APOP=0.3
HSC_ACTIV=0.2, TGF_B1=0.3, COLLAGEN=2.0
'

## ─────────────────────────────────────────────────────────────────────────────
## COMPILE MODEL
## ─────────────────────────────────────────────────────────────────────────────
nafld_mod <- mcode("nafld_shiny", MODEL_CODE, quiet = TRUE)

## ─────────────────────────────────────────────────────────────────────────────
## HELPER FUNCTIONS
## ─────────────────────────────────────────────────────────────────────────────

compute_nas <- function(df) {
  df %>% mutate(
    NAS_steatosis  = case_when(Hepatic_TG > 100 ~ 3, Hepatic_TG > 60 ~ 2,
                               Hepatic_TG > 30 ~ 1, TRUE ~ 0),
    NAS_inflam     = case_when(Kupffer_act > 0.7 ~ 3, Kupffer_act > 0.4 ~ 2,
                               Kupffer_act > 0.2 ~ 1, TRUE ~ 0),
    NAS_ballooning = case_when(Apoptosis > 1.5 ~ 2, Apoptosis > 0.5 ~ 1,
                               TRUE ~ 0),
    NAS_total      = NAS_steatosis + NAS_inflam + NAS_ballooning,
    Fibrosis_stage = case_when(Collagen > 7 ~ "F4 (Cirrhosis)",
                               Collagen > 5 ~ "F3 (Bridging)",
                               Collagen > 3.5 ~ "F2 (Perisinusoidal+Periportal)",
                               Collagen > 2 ~ "F1 (Perisinusoidal)",
                               TRUE ~ "F0 (None)"),
    Fibrosis_num   = case_when(Collagen > 7 ~ 4, Collagen > 5 ~ 3,
                               Collagen > 3.5 ~ 2, Collagen > 2 ~ 1, TRUE ~ 0)
  )
}

run_model_sim <- function(mod, fxr_dose, fxr_freq,
                          glp1_dose, glp1_freq,
                          adipo, leptin,
                          fxr_emax, glp1_emax,
                          sim_weeks) {

  mod_upd <- param(mod,
    Adipo_base = adipo,
    Leptin_base = leptin,
    Emax_FXR_DNL = fxr_emax,
    Emax_GLP1_IR = glp1_emax
  )

  events <- NULL
  if (fxr_dose > 0) {
    events <- c(events, ev(time = 0, cmt = 1, amt = fxr_dose,
                           ii = fxr_freq, addl = sim_weeks * 7 / fxr_freq - 1))
  }
  if (glp1_dose > 0) {
    events <- c(events, ev(time = 0, cmt = 4, amt = glp1_dose,
                           ii = glp1_freq, addl = sim_weeks / (glp1_freq / 168) - 1))
  }
  if (is.null(events)) events <- ev(time = 0, cmt = 1, amt = 0)

  out <- mrgsim(mod_upd,
                events = events,
                end = sim_weeks * 168,
                delta = 24)
  df <- as.data.frame(out)
  df$time_days  <- df$time / 24
  df$time_weeks <- df$time / 168
  compute_nas(df)
}

THEME_NAFLD <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        legend.position = "bottom")

SCENARIO_COLORS <- c(
  "No Treatment"         = "#E53935",
  "OCA 25mg/day"         = "#1E88E5",
  "Semaglutide 2.4mg/wk" = "#43A047",
  "OCA + Sema"           = "#8E24AA",
  "Resmetirom 80mg/day"  = "#FB8C00",
  "Custom"               = "#00ACC1"
)

## ─────────────────────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = "NAFLD/MASLD QSP",
    titleWidth = 260
  ),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "patient",   icon = icon("user-md")),
      menuItem("PK Profiles",          tabName = "pk",        icon = icon("pills")),
      menuItem("PD Biomarkers",        tabName = "pd",        icon = icon("flask")),
      menuItem("Clinical Endpoints",   tabName = "clinical",  icon = icon("chart-line")),
      menuItem("Scenario Comparison",  tabName = "scenarios", icon = icon("exchange-alt")),
      menuItem("Biomarker Panel",      tabName = "biomarker", icon = icon("vials"))
    ),
    hr(),
    h4("  Global Settings", style = "color:#ccc; padding-left:10px"),
    sliderInput("sim_weeks",  "Simulation (weeks):", 4, 208, 104, step = 4,
                width = "90%"),
    selectInput("drug_regimen", "Treatment Regimen:",
                choices = c("No Treatment",
                            "OCA 25mg/day",
                            "Semaglutide 2.4mg/wk",
                            "OCA + Sema",
                            "Resmetirom 80mg/day",
                            "Custom"),
                selected = "OCA + Sema",
                width = "90%"),
    conditionalPanel(
      condition = "input.drug_regimen == 'Custom'",
      numericInput("fxr_dose",  "FXR Agonist Dose (mg):", 25, 0, 200, step = 5,
                   width = "90%"),
      numericInput("glp1_dose", "GLP-1 RA Dose (mg):",   2.4, 0, 10, step = 0.4,
                   width = "90%")
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box-title { font-weight: bold; }
      .content-wrapper { background-color: #f4f6f9; }
      .nav-tabs-custom { background: #fff; }
    "))),

    tabItems(

      ## ── TAB 1: PATIENT PROFILE ─────────────────────────────────────────────
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "Patient Characteristics", width = 4, solidHeader = TRUE,
              status = "primary",
              sliderInput("bmi",     "BMI (kg/m²):",  20, 50, 32, step = 0.5),
              sliderInput("hba1c",   "HbA1c (%):",    4,  12,  7.2, step = 0.1),
              sliderInput("adipo",   "Adiponectin (µg/mL):", 2, 30, 8, step = 0.5),
              sliderInput("leptin",  "Leptin (ng/mL):", 2, 80, 28, step = 1),
              sliderInput("alt_init","Baseline ALT (U/L):", 20, 200, 75, step = 5),
              selectInput("fibrosis_base", "Baseline Fibrosis Stage:",
                          choices = c("F0","F1","F2","F3"), selected = "F1")
          ),
          box(title = "Disease Activity Profile", width = 4, solidHeader = TRUE,
              status = "warning",
              h4("Estimated Baseline State"),
              tableOutput("patient_summary"),
              hr(),
              h4("NAS Component Breakdown"),
              plotOutput("nas_donut", height = "220px")
          ),
          box(title = "Disease Stage & Risk", width = 4, solidHeader = TRUE,
              status = "danger",
              valueBoxOutput("nas_box",    width = 12),
              valueBoxOutput("fib_box",    width = 12),
              valueBoxOutput("cv_risk_box", width = 12),
              br(),
              h5("10-year Liver-Related Event Risk"),
              plotOutput("event_risk_bar", height = "160px")
          )
        ),
        fluidRow(
          box(title = "Mechanistic Map", width = 12, solidHeader = TRUE,
              status = "info",
              p("The NAFLD/MASLD QSP mechanistic map covers 10 biological subsystems:",
                tags$br(),
                tags$b("Adipose tissue & IR · Hepatic lipid metabolism · Oxidative/ER stress ·"),
                tags$br(),
                tags$b("Kupffer cell inflammation · Hepatocyte death · HSC fibrosis ·"),
                tags$br(),
                tags$b("Gut-liver axis · Drug PK/PD · Clinical endpoints")),
              imageOutput("mech_map", height = "500px")
          )
        )
      ),

      ## ── TAB 2: PK PROFILES ────────────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "PK Parameters", width = 3, solidHeader = TRUE,
              status = "primary",
              h4("FXR Agonist (OCA-like)"),
              numericInput("pk_CL",   "CL (L/h):",   4.5, 0.5, 20, step = 0.5),
              numericInput("pk_V1",   "V1 (L):",     15,  2,   80, step = 1),
              numericInput("pk_V2",   "V2 (L):",     30,  5,   150, step = 5),
              numericInput("pk_ka",   "ka (h⁻¹):",   0.8, 0.1, 3,  step = 0.1),
              numericInput("pk_F1",   "F (%):",       65,  10,  100, step = 5),
              hr(),
              h4("GLP-1 RA (Semaglutide-like)"),
              numericInput("pk_CLg",  "CL (L/h):",   0.055, 0.01, 0.5, step = 0.005),
              numericInput("pk_V1g",  "V1 (L):",     8.5,   2,    30,  step = 0.5),
              numericInput("pk_kag",  "ka (h⁻¹):",   0.005, 0.001, 0.05, step = 0.001)
          ),
          box(title = "FXR Agonist PK Profile", width = 9, solidHeader = TRUE,
              status = "info",
              fluidRow(
                column(4,
                  selectInput("pk_dose_fxr",  "Dose (mg):",
                              choices = c("5","10","25","50","100"), selected = "25"),
                  selectInput("pk_freq_fxr",  "Dosing Interval:",
                              choices = c("Once daily (24h)" = 24,
                                         "Twice daily (12h)" = 12,
                                         "Weekly (168h)" = 168),
                              selected = 24),
                  numericInput("pk_doses_n",  "# Doses:", 14, 1, 100)
                ),
                column(8,
                  plotlyOutput("pk_fxr_plot", height = "280px")
                )
              ),
              hr(),
              fluidRow(
                column(4,
                  selectInput("pk_dose_glp1", "GLP-1 RA Dose (mg):",
                              choices = c("0.25","0.5","1.0","2.4"), selected = "2.4"),
                  selectInput("pk_freq_glp1", "Dosing Interval:",
                              choices = c("Weekly (168h)" = 168,
                                         "Daily (24h)" = 24),
                              selected = 168)
                ),
                column(8,
                  plotlyOutput("pk_glp1_plot", height = "280px")
                )
              )
          )
        ),
        fluidRow(
          box(title = "PK Summary Table (Steady State)", width = 12,
              solidHeader = TRUE, status = "info",
              DTOutput("pk_table"))
        )
      ),

      ## ── TAB 3: PD BIOMARKERS ──────────────────────────────────────────────
      tabItem(tabName = "pd",
        fluidRow(
          box(title = "PD Display Options", width = 3, solidHeader = TRUE,
              status = "primary",
              checkboxGroupInput("pd_vars", "Show Markers:",
                choices = list(
                  "Hepatic TG"          = "Hepatic_TG",
                  "ROS Level"           = "ROS_level",
                  "ER Stress"           = "ER_stress",
                  "Kupffer Activation"  = "Kupffer_act",
                  "TNF-α"               = "TNFalpha",
                  "IL-6"                = "IL6_level",
                  "IL-1β"               = "IL1beta",
                  "HSC Activation"      = "HSC_act",
                  "TGF-β1"              = "TGFbeta1",
                  "Collagen"            = "Collagen",
                  "Apoptosis Index"     = "Apoptosis"
                ),
                selected = c("Hepatic_TG", "Kupffer_act", "Collagen",
                             "TGFbeta1", "ROS_level")
              ),
              hr(),
              selectInput("pd_scale", "Y-axis scale:",
                          choices = c("Free" = "free_y", "Fixed" = "fixed")),
              actionButton("run_pd", "Run Simulation",
                           class = "btn-primary btn-block")
          ),
          box(title = "PD Biomarker Time Courses", width = 9, solidHeader = TRUE,
              status = "warning",
              plotlyOutput("pd_plot", height = "500px")
          )
        )
      ),

      ## ── TAB 4: CLINICAL ENDPOINTS ─────────────────────────────────────────
      tabItem(tabName = "clinical",
        fluidRow(
          box(title = "NAS Score Components", width = 6, solidHeader = TRUE,
              status = "danger",
              plotlyOutput("nas_plot", height = "320px"),
              hr(),
              plotlyOutput("nas_components_plot", height = "280px")
          ),
          box(title = "Fibrosis Progression", width = 6, solidHeader = TRUE,
              status = "warning",
              plotlyOutput("fibrosis_plot", height = "320px"),
              hr(),
              fluidRow(
                valueBoxOutput("fibrosis_at_end", width = 6),
                valueBoxOutput("nas_at_end",       width = 6)
              )
          )
        ),
        fluidRow(
          box(title = "Liver Enzymes (ALT & AST)", width = 6,
              solidHeader = TRUE, status = "info",
              plotlyOutput("alt_ast_plot", height = "320px")
          ),
          box(title = "Non-Invasive Scores (FIB-4, ELF)", width = 6,
              solidHeader = TRUE, status = "primary",
              plotlyOutput("noninvasive_plot", height = "320px"),
              p("FIB-4 < 1.30: low fibrosis risk; FIB-4 > 2.67: high risk",
                style = "color: gray; font-size: 11px;")
          )
        )
      ),

      ## ── TAB 5: SCENARIO COMPARISON ────────────────────────────────────────
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "Scenarios to Compare", width = 3, solidHeader = TRUE,
              status = "primary",
              checkboxGroupInput("scen_select", "Select Scenarios:",
                choices = c("No Treatment",
                            "OCA 25mg/day",
                            "Semaglutide 2.4mg/wk",
                            "OCA + Sema",
                            "Resmetirom 80mg/day"),
                selected = c("No Treatment", "OCA 25mg/day",
                             "Semaglutide 2.4mg/wk", "OCA + Sema",
                             "Resmetirom 80mg/day")
              ),
              hr(),
              selectInput("scen_endpoint", "Primary Endpoint:",
                          choices = c("NAS Total"       = "NAS_total",
                                      "Collagen (ECM)"  = "Collagen",
                                      "Hepatic TG"      = "Hepatic_TG",
                                      "ALT"             = "ALT",
                                      "Kupffer Activ."  = "Kupffer_act",
                                      "TGF-β1"          = "TGFbeta1")),
              actionButton("run_scenarios", "Compare All",
                           class = "btn-success btn-block")
          ),
          box(title = "Comparative Time Course", width = 9, solidHeader = TRUE,
              status = "success",
              plotlyOutput("scenario_plot", height = "420px")
          )
        ),
        fluidRow(
          box(title = "Response Rate Comparison (Week 52)", width = 6,
              solidHeader = TRUE, status = "info",
              plotlyOutput("response_bar", height = "320px"),
              p("Response = NAS improvement ≥2 without fibrosis worsening,",
                "OR fibrosis improvement ≥1 stage")
          ),
          box(title = "Summary Table (Week 52 & Week 104)", width = 6,
              solidHeader = TRUE, status = "primary",
              DTOutput("scenario_table"))
        )
      ),

      ## ── TAB 6: BIOMARKER PANEL ────────────────────────────────────────────
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "Biomarker Dashboard", width = 12, solidHeader = TRUE,
              status = "primary",
              fluidRow(
                column(3,
                  h4("Lipid Markers"),
                  valueBoxOutput("tag_box",  width = 12),
                  valueBoxOutput("vldl_box", width = 12)
                ),
                column(3,
                  h4("Inflammatory"),
                  valueBoxOutput("tnf_box", width = 12),
                  valueBoxOutput("il6_box", width = 12)
                ),
                column(3,
                  h4("Fibrosis"),
                  valueBoxOutput("tgfb_box",   width = 12),
                  valueBoxOutput("col_box",    width = 12)
                ),
                column(3,
                  h4("Liver Injury"),
                  valueBoxOutput("alt2_box",   width = 12),
                  valueBoxOutput("hsc_box",    width = 12)
                )
              )
          )
        ),
        fluidRow(
          box(title = "Biomarker Heatmap Over Time", width = 8,
              solidHeader = TRUE, status = "warning",
              plotlyOutput("biomarker_heatmap", height = "400px")
          ),
          box(title = "Waterfall Plot (Week 52 vs Baseline)", width = 4,
              solidHeader = TRUE, status = "info",
              plotlyOutput("waterfall_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "Key Clinical Trial Comparison", width = 12,
              solidHeader = TRUE, status = "success",
              DTOutput("trial_table"))
        )
      )
    )
  )
)

## ─────────────────────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  ## ── REACTIVE: run model simulation ────────────────────────────────────────
  sim_data <- reactive({
    req(input$drug_regimen, input$sim_weeks)
    isolate({

      regimen <- input$drug_regimen
      fxr_dose  <- switch(regimen,
        "No Treatment"         = 0,
        "OCA 25mg/day"         = 25,
        "Semaglutide 2.4mg/wk" = 0,
        "OCA + Sema"           = 25,
        "Resmetirom 80mg/day"  = 80,
        "Custom"               = input$fxr_dose
      )
      glp1_dose <- switch(regimen,
        "No Treatment"         = 0,
        "OCA 25mg/day"         = 0,
        "Semaglutide 2.4mg/wk" = 2.4,
        "OCA + Sema"           = 2.4,
        "Resmetirom 80mg/day"  = 0,
        "Custom"               = input$glp1_dose
      )
      fxr_emax <- ifelse(regimen == "Resmetirom 80mg/day", 0.5, 0.40)
      glp1_emax <- ifelse(regimen == "Resmetirom 80mg/day", 0.0, 0.45)

      run_model_sim(
        mod       = nafld_mod,
        fxr_dose  = fxr_dose, fxr_freq = 24,
        glp1_dose = glp1_dose, glp1_freq = 168,
        adipo     = input$adipo,
        leptin    = input$leptin,
        fxr_emax  = fxr_emax,
        glp1_emax = glp1_emax,
        sim_weeks = input$sim_weeks
      )
    })
  }) %>% bindEvent(input$run_pd, input$drug_regimen, input$sim_weeks,
                   input$adipo, input$leptin, ignoreNULL = FALSE)

  ## ── TAB 1: Patient summary ──────────────────────────────────────────────
  output$patient_summary <- renderTable({
    data.frame(
      Metric    = c("BMI", "HbA1c", "Adiponectin", "Leptin", "ALT (baseline)"),
      Value     = c(input$bmi, input$hba1c, input$adipo,
                    input$leptin, input$alt_init),
      Unit      = c("kg/m²", "%", "µg/mL", "ng/mL", "U/L"),
      Reference = c("18.5–25", "4–5.6", ">10", "<15 (M) <30 (F)", "<40")
    )
  }, striped = TRUE, hover = TRUE)

  output$nas_box <- renderValueBox({
    valueBox(
      value = paste0(min(round(input$alt_init / 15), 8), "/8"),
      subtitle = "Estimated NAS Score",
      icon = icon("circle-exclamation"),
      color = ifelse(input$alt_init > 80, "red",
                     ifelse(input$alt_init > 50, "orange", "yellow"))
    )
  })

  output$fib_box <- renderValueBox({
    fib <- input$fibrosis_base
    valueBox(
      value = fib,
      subtitle = "Baseline Fibrosis Stage",
      icon = icon("layer-group"),
      color = ifelse(fib %in% c("F3","F4"), "red",
                     ifelse(fib == "F2", "orange", "yellow"))
    )
  })

  output$cv_risk_box <- renderValueBox({
    risk <- ifelse(input$bmi > 35 | input$hba1c > 8, "High", "Moderate")
    valueBox(
      value = risk,
      subtitle = "Cardiovascular Risk",
      icon = icon("heart"),
      color = ifelse(risk == "High", "red", "orange")
    )
  })

  output$nas_donut <- renderPlot({
    vals <- c(
      Steatosis  = min(floor(input$alt_init / 30), 3),
      Inflammation = min(floor(input$leptin / 20), 3),
      Ballooning = min(floor(input$alt_init / 70), 2)
    )
    par(mar = c(0, 0, 1, 0))
    pie(pmax(vals, 0.1), labels = paste0(names(vals), "\n(", vals, ")"),
        col = c("#2196F3", "#FF5722", "#9C27B0"),
        main = paste0("NAS = ", sum(vals)),
        cex = 0.9)
  })

  output$mech_map <- renderImage({
    list(src = "nafld_qsp_model.png",
         width = "100%", alt = "NAFLD/MASLD QSP Mechanistic Map")
  }, deleteFile = FALSE)

  ## ── TAB 2: PK ─────────────────────────────────────────────────────────────
  output$pk_fxr_plot <- renderPlotly({
    dose  <- as.numeric(input$pk_dose_fxr)
    freq  <- as.numeric(input$pk_freq_fxr)
    n_dos <- as.numeric(input$pk_doses_n)

    mod_pk <- param(nafld_mod,
      CL = input$pk_CL, V1 = input$pk_V1, V2 = input$pk_V2,
      ka = input$pk_ka, F1 = input$pk_F1 / 100
    )
    ev_pk <- ev(time = 0, cmt = 1, amt = dose, ii = freq, addl = n_dos - 1)
    out   <- mrgsim(mod_pk, events = ev_pk,
                    end = freq * n_dos, delta = 0.5)
    df    <- as.data.frame(out)
    df$Cp <- df$CENT / input$pk_V1

    plot_ly(df, x = ~time, y = ~Cp, type = "scatter", mode = "lines",
            line = list(color = "#1E88E5", width = 2.5)) %>%
      layout(title = "FXR Agonist Plasma Concentration",
             xaxis = list(title = "Time (h)"),
             yaxis = list(title = "Cp (µg/mL)"))
  })

  output$pk_glp1_plot <- renderPlotly({
    dose <- as.numeric(input$pk_dose_glp1)
    freq <- as.numeric(input$pk_freq_glp1)
    mod_pk <- param(nafld_mod,
      CL_glp1 = input$pk_CLg, V1_glp1 = input$pk_V1g, ka_glp1 = input$pk_kag
    )
    ev_pk <- ev(time = 0, cmt = 4, amt = dose, ii = freq, addl = 11)
    out   <- mrgsim(mod_pk, events = ev_pk,
                    end = freq * 12, delta = 1)
    df    <- as.data.frame(out)
    df$Cp <- df$CENT_GLP1 / input$pk_V1g

    plot_ly(df, x = ~time, y = ~Cp, type = "scatter", mode = "lines",
            line = list(color = "#43A047", width = 2.5)) %>%
      layout(title = "GLP-1 RA (sc) Plasma Concentration",
             xaxis = list(title = "Time (h)"),
             yaxis = list(title = "Cp (µg/mL)"))
  })

  output$pk_table <- renderDT({
    data.frame(
      Drug        = c("OCA (FXR agonist)", "Semaglutide (GLP-1 RA)",
                      "Resmetirom (THRβ agonist)"),
      `Half-life` = c("~24h (modified bile acid)", "~168h (sc weekly)",
                      "~76h (once daily)"),
      `Tmax`      = c("1–2h", "24–36h sc", "1–4h"),
      `Cmax_ss`   = c("~0.5 µg/mL", "~0.1 µg/mL", "~1.5 µg/mL"),
      `AUC_ss`    = c("~8 h·µg/mL", "~25 h·µg/mL", "~70 h·µg/mL"),
      EC50        = c("0.5 µg/mL (FXR)", "0.3 µg/mL (GLP-1R)",
                      "0.4 µg/mL (THRβ)")
    )
  }, options = list(pageLength = 5, scrollX = TRUE))

  ## ── TAB 3: PD ─────────────────────────────────────────────────────────────
  output$pd_plot <- renderPlotly({
    df <- sim_data()
    req(length(input$pd_vars) >= 1)

    df_long <- df %>%
      select(time_weeks, all_of(input$pd_vars)) %>%
      pivot_longer(cols = -time_weeks, names_to = "marker", values_to = "value")

    p <- ggplot(df_long, aes(x = time_weeks, y = value, color = marker)) +
      geom_line(linewidth = 1.2) +
      facet_wrap(~marker, scales = input$pd_scale, ncol = 2) +
      labs(x = "Time (weeks)", y = "Value", color = "Marker") +
      THEME_NAFLD +
      theme(legend.position = "none",
            strip.text = element_text(face = "bold"))
    ggplotly(p)
  })

  ## ── TAB 4: Clinical Endpoints ──────────────────────────────────────────────
  output$nas_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time_weeks, y = NAS_total)) +
      geom_line(color = "#E53935", linewidth = 1.5) +
      geom_ribbon(aes(ymin = NAS_total - 0.5, ymax = NAS_total + 0.5),
                  alpha = 0.15, fill = "#E53935") +
      geom_hline(yintercept = 5, linetype = "dashed", color = "black") +
      annotate("text", x = 1, y = 5.2, label = "MASH threshold",
               color = "black", size = 3.5) +
      scale_y_continuous(limits = c(0, 8), breaks = 0:8) +
      labs(x = "Time (weeks)", y = "NAS Score") +
      THEME_NAFLD
    ggplotly(p)
  })

  output$nas_components_plot <- renderPlotly({
    df <- sim_data()
    df_comp <- df %>%
      select(time_weeks, NAS_steatosis, NAS_inflam, NAS_ballooning) %>%
      pivot_longer(-time_weeks, names_to = "component", values_to = "score")

    p <- ggplot(df_comp, aes(x = time_weeks, y = score, fill = component)) +
      geom_area(alpha = 0.7) +
      scale_fill_manual(values = c("#64B5F6","#FFB74D","#EF9A9A"),
                        labels = c("Ballooning","Inflammation","Steatosis")) +
      labs(x = "Time (weeks)", y = "NAS Component Score", fill = "") +
      THEME_NAFLD
    ggplotly(p)
  })

  output$fibrosis_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time_weeks, y = Collagen)) +
      geom_line(color = "#5C6BC0", linewidth = 1.5) +
      geom_hline(yintercept = c(2, 3.5, 5, 7), linetype = "dotted",
                 color = "gray60") +
      annotate("text", x = max(df$time_weeks) * 0.9,
               y = c(2.2, 3.7, 5.2, 7.2),
               label = c("F1","F2","F3","F4"),
               color = "gray40", size = 3) +
      labs(x = "Time (weeks)", y = "Collagen Content (normalized)") +
      THEME_NAFLD
    ggplotly(p)
  })

  output$alt_ast_plot <- renderPlotly({
    df <- sim_data() %>%
      mutate(AST = ALT * 0.7) %>%
      select(time_weeks, ALT, AST) %>%
      pivot_longer(-time_weeks, names_to = "enzyme", values_to = "value")

    p <- ggplot(df, aes(x = time_weeks, y = value, color = enzyme)) +
      geom_line(linewidth = 1.3) +
      scale_color_manual(values = c("ALT" = "#E53935", "AST" = "#FB8C00")) +
      geom_hline(yintercept = 40, linetype = "dashed", color = "gray50") +
      labs(x = "Time (weeks)", y = "Liver Enzymes (U/L)", color = "") +
      THEME_NAFLD
    ggplotly(p)
  })

  output$noninvasive_plot <- renderPlotly({
    df <- sim_data() %>%
      mutate(
        FIB4  = pmax(1.0 + (Collagen - 2) * 0.5, 0.5),
        ELF   = pmax(7.7 + 0.5 * log(Collagen + 0.1) + 0.3 * TGFbeta1, 5)
      ) %>%
      select(time_weeks, FIB4, ELF) %>%
      pivot_longer(-time_weeks, names_to = "score", values_to = "value")

    p <- ggplot(df, aes(x = time_weeks, y = value, color = score)) +
      geom_line(linewidth = 1.3) +
      facet_wrap(~score, scales = "free_y") +
      scale_color_manual(values = c("FIB4" = "#1565C0", "ELF" = "#6A1B9A")) +
      labs(x = "Time (weeks)", y = "Score", color = "") +
      THEME_NAFLD
    ggplotly(p)
  })

  output$fibrosis_at_end <- renderValueBox({
    df <- sim_data()
    last <- tail(df, 1)
    stage <- last$Fibrosis_stage
    valueBox(stage, "Final Fibrosis Stage",
             icon = icon("layer-group"),
             color = ifelse(grepl("F4|F3", stage), "red",
                            ifelse(grepl("F2", stage), "orange", "green")))
  })

  output$nas_at_end <- renderValueBox({
    df <- sim_data()
    last_nas <- tail(df$NAS_total, 1)
    valueBox(round(last_nas, 1), "Final NAS Score",
             icon = icon("star"),
             color = ifelse(last_nas >= 5, "red",
                            ifelse(last_nas >= 3, "orange", "green")))
  })

  ## ── TAB 5: Scenario Comparison ─────────────────────────────────────────────
  scenario_data <- eventReactive(input$run_scenarios, {
    scenarios <- input$scen_select
    configs <- list(
      "No Treatment"         = list(fxr = 0,    glp1 = 0,   fxr_e = 0.40, glp1_e = 0.45),
      "OCA 25mg/day"         = list(fxr = 25,   glp1 = 0,   fxr_e = 0.40, glp1_e = 0),
      "Semaglutide 2.4mg/wk" = list(fxr = 0,    glp1 = 2.4, fxr_e = 0,    glp1_e = 0.45),
      "OCA + Sema"           = list(fxr = 25,   glp1 = 2.4, fxr_e = 0.40, glp1_e = 0.45),
      "Resmetirom 80mg/day"  = list(fxr = 80,   glp1 = 0,   fxr_e = 0.50, glp1_e = 0)
    )

    bind_rows(lapply(scenarios, function(s) {
      cfg <- configs[[s]]
      df  <- run_model_sim(
        nafld_mod,
        fxr_dose = cfg$fxr, fxr_freq = 24,
        glp1_dose = cfg$glp1, glp1_freq = 168,
        adipo = input$adipo, leptin = input$leptin,
        fxr_emax = cfg$fxr_e, glp1_emax = cfg$glp1_e,
        sim_weeks = input$sim_weeks
      )
      df$scenario <- s
      df
    }))
  }, ignoreNULL = FALSE)

  output$scenario_plot <- renderPlotly({
    df  <- scenario_data()
    var <- input$scen_endpoint

    p <- ggplot(df, aes_string(x = "time_weeks", y = var,
                               color = "scenario")) +
      geom_line(linewidth = 1.3) +
      scale_color_manual(values = SCENARIO_COLORS,
                         breaks = names(SCENARIO_COLORS)) +
      labs(x = "Time (weeks)", y = var, color = "Treatment") +
      THEME_NAFLD
    ggplotly(p) %>% layout(legend = list(orientation = "h", y = -0.2))
  })

  output$response_bar <- renderPlotly({
    df <- scenario_data()
    bas <- df %>% filter(time_weeks == 0) %>%
      select(scenario, NAS_base = NAS_total, Fib_base = Fibrosis_num)

    resp <- df %>%
      filter(abs(time_weeks - 52) == min(abs(df$time_weeks - 52))) %>%
      left_join(bas, by = "scenario") %>%
      mutate(
        NAS_resp  = (NAS_base - NAS_total) >= 2,
        Fib_resp  = (Fib_base - Fibrosis_num) >= 1,
        Any_resp  = NAS_resp | Fib_resp
      ) %>%
      group_by(scenario) %>%
      summarise(Response_pct = mean(Any_resp) * 100, .groups = "drop")

    plot_ly(resp, x = ~scenario, y = ~Response_pct, type = "bar",
            marker = list(color = unname(SCENARIO_COLORS[resp$scenario])),
            text = ~paste0(round(Response_pct, 1), "%"),
            textposition = "outside") %>%
      layout(yaxis = list(title = "Responder Rate (%)", range = c(0, 100)),
             xaxis = list(title = ""),
             title = "Response Rate at Week 52")
  })

  output$scenario_table <- renderDT({
    df <- scenario_data()
    time_pts <- c(52, 104)
    bind_rows(lapply(time_pts, function(wk) {
      df %>%
        filter(abs(time_weeks - wk) == min(abs(df$time_weeks - wk))) %>%
        group_by(scenario) %>%
        summarise(
          Week           = wk,
          NAS            = round(mean(NAS_total), 2),
          Fibrosis_norm  = round(mean(Collagen), 2),
          ALT_UL         = round(mean(ALT), 1),
          Kupffer        = round(mean(Kupffer_act), 3),
          TGFbeta1       = round(mean(TGFbeta1), 3),
          .groups = "drop"
        )
    })) %>% arrange(Week, scenario)
  }, options = list(pageLength = 10, scrollX = TRUE))

  ## ── TAB 6: Biomarker Panel ─────────────────────────────────────────────────
  output$tag_box <- renderValueBox({
    df <- sim_data()
    v  <- round(tail(df$Hepatic_TG, 1), 1)
    valueBox(v, "Hepatic TG (µmol/g)",
             icon = icon("droplet"),
             color = ifelse(v > 80, "red", ifelse(v > 40, "orange", "green")))
  })

  output$vldl_box <- renderValueBox({
    valueBox("↑", "VLDL (hypertrig.)",
             icon = icon("heart-pulse"),
             color = "orange")
  })

  output$tnf_box <- renderValueBox({
    df <- sim_data()
    v  <- round(tail(df$TNFalpha, 1), 3)
    valueBox(v, "TNF-α (AU)",
             icon = icon("fire"),
             color = ifelse(v > 0.4, "red", ifelse(v > 0.2, "orange", "green")))
  })

  output$il6_box <- renderValueBox({
    df <- sim_data()
    v  <- round(tail(df$IL6_level, 1), 3)
    valueBox(v, "IL-6 (AU)",
             icon = icon("wind"),
             color = ifelse(v > 0.3, "red", "orange"))
  })

  output$tgfb_box <- renderValueBox({
    df <- sim_data()
    v  <- round(tail(df$TGFbeta1, 1), 3)
    valueBox(v, "TGF-β1 (AU)",
             icon = icon("layer-group"),
             color = ifelse(v > 0.4, "red", ifelse(v > 0.2, "orange", "green")))
  })

  output$col_box <- renderValueBox({
    df <- sim_data()
    v  <- round(tail(df$Collagen, 1), 2)
    valueBox(v, "Collagen (normalized)",
             icon = icon("bars"),
             color = ifelse(v > 5, "red", ifelse(v > 3.5, "orange", "green")))
  })

  output$alt2_box <- renderValueBox({
    df <- sim_data()
    v  <- round(tail(df$ALT, 1), 1)
    valueBox(v, "ALT (U/L)",
             icon = icon("syringe"),
             color = ifelse(v > 80, "red", ifelse(v > 40, "orange", "green")))
  })

  output$hsc_box <- renderValueBox({
    df <- sim_data()
    v  <- round(tail(df$HSC_act, 1), 3)
    valueBox(v, "HSC Activation (0–1)",
             icon = icon("star-of-life"),
             color = ifelse(v > 0.5, "red", ifelse(v > 0.3, "orange", "green")))
  })

  output$biomarker_heatmap <- renderPlotly({
    df <- sim_data()
    vars <- c("Hepatic_TG","ROS_level","ER_stress","Kupffer_act",
              "TNFalpha","IL6_level","TGFbeta1","Collagen","Apoptosis")
    df_heat <- df %>%
      filter(time_weeks %% 4 == 0) %>%
      select(time_weeks, all_of(vars)) %>%
      mutate(across(all_of(vars), ~scale(.) %>% as.numeric())) %>%
      pivot_longer(-time_weeks, names_to = "marker", values_to = "z_score")

    plot_ly(df_heat, x = ~time_weeks, y = ~marker, z = ~z_score,
            type = "heatmap",
            colorscale = list(c(0,"#2196F3"), c(0.5,"white"), c(1,"#E53935"))) %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = ""),
             title = "Z-score Normalized Biomarker Heatmap")
  })

  output$waterfall_plot <- renderPlotly({
    df <- sim_data()
    baseline <- df %>% filter(time_weeks == 0) %>%
      select(Hepatic_TG, ROS_level, Kupffer_act, TNFalpha,
             TGFbeta1, Collagen, Apoptosis)
    wk52 <- df %>%
      filter(abs(time_weeks - 52) == min(abs(df$time_weeks - 52))) %>%
      select(Hepatic_TG, ROS_level, Kupffer_act, TNFalpha,
             TGFbeta1, Collagen, Apoptosis)

    delta <- ((as.numeric(wk52[1,]) - as.numeric(baseline[1,])) /
               (as.numeric(baseline[1,]) + 1e-9)) * 100
    wf_df <- data.frame(
      marker = names(baseline),
      pct_change = delta
    ) %>% arrange(pct_change)

    plot_ly(wf_df, x = ~pct_change, y = ~reorder(marker, pct_change),
            type = "bar", orientation = "h",
            marker = list(color = ifelse(wf_df$pct_change < 0,
                                         "#43A047", "#E53935"))) %>%
      layout(xaxis = list(title = "% Change from Baseline"),
             yaxis = list(title = ""),
             title = "Biomarker Change at Week 52")
  })

  output$trial_table <- renderDT({
    data.frame(
      Trial       = c("MAESTRO-NASH (resmetirom 80mg)",
                      "MAESTRO-NASH (resmetirom 100mg)",
                      "REGENERATE (OCA 25mg)",
                      "CENTAUR (cenicriviroc)",
                      "LEAN (liraglutide 1.8mg)",
                      "NATIVE (semaglutide 2.4mg)"),
      Duration    = c("52 wk","52 wk","18 mo","52 wk","48 wk","72 wk"),
      NAS_response = c("25.9%","29.9%","—","—","39%","59%"),
      Fibrosis_response = c("24.2%","25.9%","23%","20%","—","43%"),
      Primary_endpoint = c("NAS≥2 + F stable","NAS≥2 + F stable",
                           "F≥1 no NAS worsen","F≥1 no NAS worsen",
                           "NAS≥2 + F stable","NASH resolv"),
      Reference   = c("Harrison NEJM 2023","Harrison NEJM 2023",
                      "Sanyal Hepatol 2021","Friedman Hepatol 2018",
                      "Armstrong Lancet 2016","Newsome Lancet 2021")
    )
  }, options = list(pageLength = 6, scrollX = TRUE))
}

## ─────────────────────────────────────────────────────────────────────────────
## LAUNCH
## ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
