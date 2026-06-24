## =============================================================================
## Cholelithiasis (Gallstone Disease) QSP — Interactive Shiny Dashboard
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)
library(scales)

## ---- Embedded mrgsolve model ------------------------------------------------
model_code <- '
$PARAM
DOSE_UDCA=750, BWT=70, ka_UDCA=0.80, F_UDCA=0.50,
Vd_UDCA=12.0, CL_UDCA=40.0, kEHC_UDCA=0.35, f_hep=0.70, f_bile=0.80,
DOSE_STAT=0, ka_STAT=0.60, F_STAT=0.05, Vd_STAT=3.2, CL_STAT=70.0,
Km_STAT=0.008, Emax_STAT=0.85,
DOSE_EZET=0, ka_EZET=0.50, F_EZET=0.35, Vd_EZET=4.0, CL_EZET=8.0,
Km_EZET=0.006, Emax_EZET=0.80,
BA_synth0=0.52, BA_pool0=3.5, kBA_EHC=0.42, kBA_fecal=0.014,
E_UDCA_BA=0.25, KD_FXR=8.0,
CHOL_h0=15.0, k_CHOL_syn=1.50, k_CHOL_deg=0.10, E_FXR_CHOL=0.20,
k_CHOL_bil=0.080, k_PL_bil=0.18, PL_bil0=12.0, BA_bil0=35.0, CHOL_bil0=4.2,
GB_vol0=30.0, GB_vol_min=5.0, k_GB_fill=0.025, CCK_peak=1.0, k_GB_empty=0.30,
CSI_thresh=1.05, k_nucleat=0.0005, k_grow=0.012, E_UDCA_dis=0.025,
Stone_vol0=0.5, Stone_max=5.0,
IL6_base=2.0, CRP_base=0.5, k_IL6_prod=0.05, k_IL6_elim=0.15,
k_CRP_prod=0.30, k_CRP_elim=0.035, WLOSS=0, k_WL=0.002

$CMT
A_gut_UDCA A_plas_UDCA A_hep_UDCA A_bile_UDCA A_gb_UDCA
A_gut_STAT A_plas_STAT
A_gut_EZET A_plas_EZET
BA_pool CHOL_h CHOL_bil PL_bil
GB_vol Crystal_mass Stone_V
IL6 CRP_plas

$MAIN
double Vd_UDCA_L = Vd_UDCA * BWT;
double Vd_STAT_L = Vd_STAT * BWT;
double Vd_EZET_L = Vd_EZET * BWT;
double C_UDCA_bile = (A_bile_UDCA / 392.6) * 1000.0;
double C_STAT_plas = A_plas_STAT / Vd_STAT_L;
double C_EZET_plas = A_plas_EZET / Vd_EZET_L;
double E_STAT = (Emax_STAT * C_STAT_plas) / (Km_STAT + C_STAT_plas);
double E_EZET = (Emax_EZET * C_EZET_plas) / (Km_EZET + C_EZET_plas);
double BA_conc_portal = (BA_pool / 3.5) * 40.0;
double FXR_act = BA_conc_portal / (KD_FXR + BA_conc_portal);
double k_CHOL_syn_eff = k_CHOL_syn * (1.0 - E_STAT) * (1.0 + 0.3*(1.0 - CHOL_h/CHOL_h0));
double C_UDCA_bile_norm = C_UDCA_bile / 500.0;
double E_UDCA_CSI = E_FXR_CHOL + E_UDCA_dis * C_UDCA_bile_norm;
if(E_UDCA_CSI > 0.70) E_UDCA_CSI = 0.70;
double k_CHOL_bil_eff = k_CHOL_bil * (1.0 - E_UDCA_CSI) * (1.0 - E_STAT * 0.3);
double BA_bil = BA_pool / 0.10;
double CSI_val = CHOL_bil / (0.1875 * BA_bil + 0.1429 * PL_bil + 1e-6);
double k_dissol_eff = E_UDCA_dis * C_UDCA_bile_norm * 0.5;
double delta_CSI = (CSI_val > CSI_thresh) ? (CSI_val - CSI_thresh) : 0.0;
double nucleat_rate = k_nucleat * delta_CSI * GB_vol;
double growth_rate = k_grow * delta_CSI * Stone_V * (1.0 - Stone_V / Stone_max);
double k_IL6_stim = k_IL6_prod * Stone_V * (Stone_V > 0.1 ? 1.0 : 0.0);
if(NEWIND <= 1) {
  _INIT(Stone_V) = Stone_vol0; _INIT(Crystal_mass) = Stone_vol0 * 100.0;
  _INIT(BA_pool) = BA_pool0; _INIT(CHOL_h) = CHOL_h0;
  _INIT(CHOL_bil) = CHOL_bil0; _INIT(PL_bil) = PL_bil0;
  _INIT(GB_vol) = GB_vol0; _INIT(IL6) = IL6_base; _INIT(CRP_plas) = CRP_base;
}

$ODE
dxdt_A_gut_UDCA  = -ka_UDCA * A_gut_UDCA;
double UDCA_abs = ka_UDCA * A_gut_UDCA;
dxdt_A_plas_UDCA = UDCA_abs * (1.0 - f_hep) - (CL_UDCA / Vd_UDCA_L) * A_plas_UDCA;
dxdt_A_hep_UDCA  = UDCA_abs * f_hep - kEHC_UDCA * A_hep_UDCA;
dxdt_A_bile_UDCA = f_bile * kEHC_UDCA * A_hep_UDCA - 0.50 * kEHC_UDCA * A_bile_UDCA;
dxdt_A_gb_UDCA   = 0.40 * f_bile * kEHC_UDCA * A_hep_UDCA - kEHC_UDCA * A_gb_UDCA;
dxdt_A_gut_STAT  = -ka_STAT * A_gut_STAT;
dxdt_A_plas_STAT = ka_STAT * A_gut_STAT * F_STAT - (CL_STAT / Vd_STAT_L) * A_plas_STAT;
dxdt_A_gut_EZET  = -ka_EZET * A_gut_EZET;
dxdt_A_plas_EZET = ka_EZET * A_gut_EZET * F_EZET - (CL_EZET / Vd_EZET_L) * A_plas_EZET;
double BA_syn_rate = (BA_synth0 / 24.0) * (1.0 - 0.50 * FXR_act) * (1.0 + E_UDCA_BA * C_UDCA_bile_norm);
dxdt_BA_pool = BA_syn_rate - kBA_fecal * BA_pool;
dxdt_CHOL_h  = k_CHOL_syn_eff + 0.20 * (1.0 + E_STAT * 0.8) - k_CHOL_deg * CHOL_h - k_CHOL_bil_eff;
dxdt_CHOL_bil = k_CHOL_bil_eff - 0.08 * CHOL_bil;
dxdt_PL_bil   = k_PL_bil - 0.06 * PL_bil;
dxdt_GB_vol   = k_GB_fill * (GB_vol0 - GB_vol) - k_GB_empty * CCK_peak * GB_vol;
dxdt_Crystal_mass = nucleat_rate * 50.0 + growth_rate * 80.0 - k_dissol_eff * Crystal_mass;
dxdt_Stone_V  = growth_rate - k_dissol_eff * Stone_V;
dxdt_IL6      = IL6_base * k_IL6_elim + k_IL6_stim - k_IL6_elim * IL6;
dxdt_CRP_plas = k_CRP_prod * IL6 / (IL6_base + 1.0) * CRP_base - k_CRP_elim * CRP_plas;

$TABLE
double CSI_out = CHOL_bil / (0.1875 * (BA_pool / 0.10) + 0.1429 * PL_bil + 1e-6);
double UDCA_bile_uM = (A_bile_UDCA / 392.6) * 1000.0;
double UDCA_plas_mgL = A_plas_UDCA / (Vd_UDCA * BWT);
double STAT_plas_mgL = A_plas_STAT / (Vd_STAT * BWT);
double Stone_mm = pow(Stone_V * 6.0 / 3.14159, 0.333) * 10.0;
capture CSI_out UDCA_bile_uM UDCA_plas_mgL STAT_plas_mgL Stone_V Stone_mm
capture BA_pool CHOL_bil PL_bil CHOL_h IL6 CRP_plas GB_vol Crystal_mass
'

## ---- Compile model on app load ----------------------------------------------
mod <- mcode("chol_shiny", model_code, quiet = TRUE)

## ---- Helper function to run one scenario ------------------------------------
run_sim <- function(mod, dose_UDCA, dose_STAT, dose_EZET, bwt,
                    stone_init, wloss, dur_days = 365) {
  ev_list <- list()
  if(dose_UDCA > 0)
    ev_list$u <- ev(amt = dose_UDCA/3, cmt = "A_gut_UDCA", ii = 8,
                    addl = dur_days * 3 - 1)
  if(dose_STAT > 0)
    ev_list$s <- ev(amt = dose_STAT,   cmt = "A_gut_STAT", ii = 24,
                    addl = dur_days - 1)
  if(dose_EZET > 0)
    ev_list$e <- ev(amt = dose_EZET,   cmt = "A_gut_EZET", ii = 24,
                    addl = dur_days - 1)

  ev_final <- if(length(ev_list) == 0) ev(time=0, amt=0, cmt=1)
              else Reduce(ev_seq, ev_list)

  mod %>%
    param(DOSE_UDCA = dose_UDCA, DOSE_STAT = dose_STAT, DOSE_EZET = dose_EZET,
          BWT = bwt, Stone_vol0 = stone_init, WLOSS = wloss) %>%
    mrgsim(events = ev_final, end = dur_days * 24, delta = 24, obsonly = TRUE) %>%
    as_tibble() %>%
    mutate(Day = time / 24)
}

## ============================================================================
## UI
## ============================================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Cholelithiasis QSP Dashboard", titleWidth = 320),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "tab_patient",  icon = icon("user")),
      menuItem("Drug PK",            tabName = "tab_pk",       icon = icon("pills")),
      menuItem("Bile & CSI",         tabName = "tab_csi",      icon = icon("flask")),
      menuItem("Gallstone Dynamics", tabName = "tab_stone",    icon = icon("circle")),
      menuItem("Scenario Comparison",tabName = "tab_compare",  icon = icon("chart-bar")),
      menuItem("Biomarkers",         tabName = "tab_bio",      icon = icon("heartbeat"))
    ),
    hr(),
    h5("  Global Settings", style="color:#ccc; margin-left:10px;"),
    sliderInput("dur", "Simulation Duration (days)", 90, 730, 365, step = 30),
    br()
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f9f9f9; }
      .box { border-radius: 6px; }
    "))),

    tabItems(

      ## -------- TAB 1: Patient Profile ----------------------------------------
      tabItem("tab_patient",
        fluidRow(
          box(title = "Patient Characteristics", width = 4, status = "primary",
              solidHeader = TRUE,
              numericInput("bwt",    "Body Weight (kg)",     70, 40, 150, 5),
              selectInput("sex",     "Sex",
                          c("Male" = "M", "Female (premenopausal)" = "F_pre",
                            "Female (postmenopausal)" = "F_post")),
              numericInput("age",    "Age (years)",          50, 18, 90, 1),
              numericInput("bmi",    "BMI (kg/m²)",          28, 16, 55, 0.5),
              selectInput("diet",   "Dietary Pattern",
                          c("Western (High fat)" = "western",
                            "Mediterranean"       = "med",
                            "Low-calorie"         = "low_cal")),
              checkboxInput("rapid_wl", "Rapid Weight Loss (>1.5 kg/week)", FALSE),
              checkboxInput("crohn",   "Crohn's Disease", FALSE),
              checkboxInput("cirrhosis","Cirrhosis", FALSE)
          ),
          box(title = "Gallstone Characteristics", width = 4, status = "warning",
              solidHeader = TRUE,
              selectInput("stone_type",  "Stone Type",
                          c("Cholesterol (radiolucent)" = "chol",
                            "Mixed"                     = "mixed",
                            "Pigment (black)"           = "pigment_black",
                            "Pigment (brown)"           = "pigment_brown")),
              numericInput("stone_init_mm", "Initial Stone Diameter (mm)", 8, 0, 30, 1),
              selectInput("stone_loc",   "Stone Location",
                          c("Gallbladder"    = "gb",
                            "Cystic Duct"    = "cystic",
                            "Common Bile Duct"= "cbd")),
              checkboxInput("sx_colic",  "Symptomatic (biliary colic)", TRUE),
              checkboxInput("sx_cholecyst","History of cholecystitis", FALSE)
          ),
          box(title = "Risk Score", width = 4, status = "danger",
              solidHeader = TRUE,
              verbatimTextOutput("risk_text"),
              br(),
              h5("UDCA Candidacy for Dissolution Therapy"),
              verbatimTextOutput("udca_candidate")
          )
        ),
        fluidRow(
          box(title = "Clinical Decision Summary", width = 12, status = "info",
              solidHeader = TRUE,
              tableOutput("patient_summary")
          )
        )
      ),

      ## -------- TAB 2: Drug PK -----------------------------------------------
      tabItem("tab_pk",
        fluidRow(
          box(title = "Drug Selection & Doses", width = 4, status = "primary",
              solidHeader = TRUE,
              checkboxInput("use_udca", "Use UDCA", TRUE),
              conditionalPanel("input.use_udca == true",
                sliderInput("dose_udca", "UDCA Daily Dose (mg/day)",
                            250, 1500, 750, step = 250),
                selectInput("udca_freq", "Dosing Frequency",
                            c("TID (every 8h)"  = "TID",
                              "BID (every 12h)" = "BID",
                              "QD (once daily)" = "QD"))
              ),
              hr(),
              checkboxInput("use_stat", "Add Statin (Simvastatin)", FALSE),
              conditionalPanel("input.use_stat == true",
                selectInput("dose_stat", "Simvastatin Dose",
                            c("10 mg" = 10, "20 mg" = 20, "40 mg" = 40, "80 mg" = 80))
              ),
              hr(),
              checkboxInput("use_ezet", "Add Ezetimibe", FALSE),
              conditionalPanel("input.use_ezet == true",
                numericInput("dose_ezet", "Ezetimibe Dose (mg/day)", 10, 5, 20, 5)
              ),
              actionButton("run_sim", "Run Simulation", icon = icon("play"),
                           class = "btn-primary btn-block", style = "margin-top:15px;")
          ),
          box(title = "UDCA Plasma Concentration", width = 8, status = "info",
              solidHeader = TRUE,
              plotlyOutput("pk_udca_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "UDCA Biliary Concentration (µmol/L)", width = 6,
              status = "success", solidHeader = TRUE,
              plotlyOutput("pk_bile_udca", height = "260px")),
          box(title = "Statin & Ezetimibe Plasma", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("pk_stat_plot", height = "260px"))
        )
      ),

      ## -------- TAB 3: Bile Acid & CSI Dynamics --------------------------------
      tabItem("tab_csi",
        fluidRow(
          box(title = "Cholesterol Saturation Index (CSI)", width = 8,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("csi_plot", height = "340px")),
          box(title = "CSI Interpretation", width = 4,
              status = "info", solidHeader = TRUE,
              valueBoxOutput("csi_val_box", width = 12),
              br(),
              p("CSI = Biliary cholesterol / (0.1875 × [BA] + 0.1429 × [PL])"),
              p(strong("CSI < 1.0:"), " Unsaturated — stones cannot form"),
              p(strong("CSI = 1.0:"), " Saturation point (lithogenic threshold)"),
              p(strong("CSI > 1.0:"), " Supersaturated — nucleation can begin"),
              hr(),
              p("UDCA reduces CSI by enriching the hydrophilic BA fraction,
                 increasing micellar solubilization of cholesterol.")
          )
        ),
        fluidRow(
          box(title = "Total Bile Acid Pool (g)", width = 6,
              status = "success", solidHeader = TRUE,
              plotlyOutput("ba_pool_plot", height = "260px")),
          box(title = "Biliary Composition (mmol/L)", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("bil_comp_plot", height = "260px"))
        )
      ),

      ## -------- TAB 4: Gallstone Dynamics --------------------------------------
      tabItem("tab_stone",
        fluidRow(
          box(title = "Gallstone Volume Over Time", width = 8,
              status = "primary", solidHeader = TRUE,
              plotlyOutput("stone_vol_plot", height = "340px")),
          box(title = "Dissolution Progress", width = 4,
              status = "success", solidHeader = TRUE,
              valueBoxOutput("stone_dissolv_box", width = 12),
              br(),
              p(strong("UDCA Dissolution Criteria:")),
              tags$ul(
                tags$li("Cholesterol stones (radiolucent on X-ray)"),
                tags$li("Diameter < 5 mm: best response"),
                tags$li("Diameter 5–10 mm: good response (months)"),
                tags$li("Diameter > 10 mm: poor dissolution"),
                tags$li("Treatment duration: 6–18 months"),
                tags$li("Recurrence rate: 30–50% within 5 years")
              )
          )
        ),
        fluidRow(
          box(title = "Crystal Mass Dynamics (mg)", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("crystal_plot", height = "260px")),
          box(title = "Gallbladder Volume (mL)", width = 6,
              status = "info", solidHeader = TRUE,
              plotlyOutput("gb_vol_plot", height = "260px"))
        )
      ),

      ## -------- TAB 5: Scenario Comparison -------------------------------------
      tabItem("tab_compare",
        fluidRow(
          box(title = "Multi-Scenario Configuration", width = 12,
              status = "primary", solidHeader = TRUE, collapsible = TRUE,
              fluidRow(
                column(4,
                  h4("Scenario A"),
                  selectInput("sc_a", NULL,
                    c("No treatment"           = "none",
                      "UDCA 750 mg/day"        = "udca750",
                      "UDCA 1050 mg/day"       = "udca1050",
                      "UDCA + Simvastatin 40mg"= "udca_stat",
                      "Ezetimibe 10 mg"        = "ezet",
                      "UDCA + Weight Loss"     = "udca_wl"), selected = "none")
                ),
                column(4,
                  h4("Scenario B"),
                  selectInput("sc_b", NULL,
                    c("No treatment"           = "none",
                      "UDCA 750 mg/day"        = "udca750",
                      "UDCA 1050 mg/day"       = "udca1050",
                      "UDCA + Simvastatin 40mg"= "udca_stat",
                      "Ezetimibe 10 mg"        = "ezet",
                      "UDCA + Weight Loss"     = "udca_wl"), selected = "udca750")
                ),
                column(4,
                  h4("Scenario C"),
                  selectInput("sc_c", NULL,
                    c("No treatment"           = "none",
                      "UDCA 750 mg/day"        = "udca750",
                      "UDCA 1050 mg/day"       = "udca1050",
                      "UDCA + Simvastatin 40mg"= "udca_stat",
                      "Ezetimibe 10 mg"        = "ezet",
                      "UDCA + Weight Loss"     = "udca_wl"), selected = "udca_stat")
                )
              ),
              actionButton("run_compare", "Compare Scenarios",
                           icon = icon("chart-bar"), class = "btn-success")
          )
        ),
        fluidRow(
          box(title = "Stone Volume Comparison", width = 6,
              status = "primary", solidHeader = TRUE,
              plotlyOutput("compare_stone", height = "300px")),
          box(title = "CSI Comparison", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("compare_csi", height = "300px"))
        ),
        fluidRow(
          box(title = "Endpoint Summary Table (Day 0, 180, 365)", width = 12,
              status = "info", solidHeader = TRUE,
              DTOutput("compare_table"))
        )
      ),

      ## -------- TAB 6: Biomarkers & Clinical Endpoints -------------------------
      tabItem("tab_bio",
        fluidRow(
          box(title = "Inflammatory Biomarkers", width = 8,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("inflam_plot", height = "300px")),
          box(title = "Biomarker Interpretation", width = 4,
              status = "info", solidHeader = TRUE,
              valueBoxOutput("crp_box", width = 12),
              valueBoxOutput("il6_box", width = 12),
              br(),
              p(strong("CRP Reference Ranges:")),
              tags$ul(
                tags$li("Normal: < 1 mg/L"),
                tags$li("Low-grade inflammation: 1–10 mg/L"),
                tags$li("Acute phase (cholecystitis): > 10 mg/L"),
                tags$li("Severe infection/cholangitis: > 50 mg/L")
              )
          )
        ),
        fluidRow(
          box(title = "Hepatic Cholesterol Dynamics", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("chol_hep_plot", height = "260px")),
          box(title = "Liver Function Trends", width = 6,
              status = "success", solidHeader = TRUE,
              plotlyOutput("lft_plot", height = "260px"))
        ),
        fluidRow(
          box(title = "UDCA Dose-Response: Stone Dissolution at 12 Months",
              width = 12, status = "primary", solidHeader = TRUE,
              plotlyOutput("dose_resp_plot", height = "260px"))
        )
      )
    )  # end tabItems
  )
)

## ============================================================================
## Server
## ============================================================================
server <- function(input, output, session) {

  ## ---- Scenario parameter lookup -------------------------------------------
  get_scenario_params <- function(sc_id) {
    switch(sc_id,
      "none"      = list(DOSE_UDCA = 0,    DOSE_STAT = 0,  DOSE_EZET = 0,  WLOSS = 0),
      "udca750"   = list(DOSE_UDCA = 750,  DOSE_STAT = 0,  DOSE_EZET = 0,  WLOSS = 0),
      "udca1050"  = list(DOSE_UDCA = 1050, DOSE_STAT = 0,  DOSE_EZET = 0,  WLOSS = 0),
      "udca_stat" = list(DOSE_UDCA = 750,  DOSE_STAT = 40, DOSE_EZET = 0,  WLOSS = 0),
      "ezet"      = list(DOSE_UDCA = 0,    DOSE_STAT = 0,  DOSE_EZET = 10, WLOSS = 0),
      "udca_wl"   = list(DOSE_UDCA = 750,  DOSE_STAT = 0,  DOSE_EZET = 0,  WLOSS = 1)
    )
  }

  ## ---- Reactive simulation (main tab) -------------------------------------
  sim_data <- eventReactive(input$run_sim, {
    dose_udca <- if(input$use_udca) input$dose_udca else 0
    dose_stat <- if(input$use_stat) as.numeric(input$dose_stat) else 0
    dose_ezet <- if(input$use_ezet) input$dose_ezet else 0
    wloss     <- if(input$rapid_wl) 1 else 0

    run_sim(mod, dose_UDCA = dose_udca, dose_STAT = dose_stat,
            dose_EZET = dose_ezet, bwt = input$bwt,
            stone_init = (input$stone_init_mm / 10)^3 * pi / 6,
            wloss = wloss, dur_days = input$dur)
  }, ignoreNULL = FALSE)

  ## ---- Risk score ----------------------------------------------------------
  output$risk_text <- renderText({
    score <- 0
    if(input$sex != "M") score <- score + 2
    if(input$bmi > 30)   score <- score + 2
    if(input$age > 40)   score <- score + 1
    if(input$rapid_wl)   score <- score + 2
    if(input$crohn)      score <- score + 1
    if(input$cirrhosis)  score <- score + 2
    risk_cat <- if(score >= 6) "HIGH" else if(score >= 3) "MODERATE" else "LOW"
    paste("Risk Score:", score, "/10\nCategory:", risk_cat,
          "\n\nKey risk factors present:\n",
          if(input$sex != "M") "• Female sex (+2)\n" else "",
          if(input$bmi > 30)   "• Obesity BMI>30 (+2)\n" else "",
          if(input$rapid_wl)   "• Rapid weight loss (+2)\n" else "",
          if(input$crohn)      "• Crohn's disease (+1)\n" else "",
          if(input$cirrhosis)  "• Cirrhosis (+2)\n" else "")
  })

  output$udca_candidate <- renderText({
    if(input$stone_type == "chol" && input$stone_init_mm <= 15 &&
       input$stone_loc == "gb") {
      "✓ CANDIDATE for UDCA dissolution therapy\n(Cholesterol stone, GB location, ≤15mm)"
    } else if(input$stone_type %in% c("pigment_black","pigment_brown")) {
      "✗ NOT candidate: Pigment stone\n(UDCA ineffective for pigment stones)"
    } else {
      "△ Partial candidate: assess stone size\nand type by CT/HIDA scan"
    }
  })

  output$patient_summary <- renderTable({
    stone_vol <- (input$stone_init_mm / 10)^3 * pi / 6
    data.frame(
      Parameter = c("Body Weight", "BMI", "Stone Diameter", "Stone Volume",
                    "Stone Type", "Dosing Duration"),
      Value     = c(paste(input$bwt, "kg"),
                    paste(input$bmi, "kg/m²"),
                    paste(input$stone_init_mm, "mm"),
                    paste(round(stone_vol, 3), "mL"),
                    input$stone_type,
                    paste(input$dur, "days"))
    )
  })

  ## ---- PK plots -----------------------------------------------------------
  output$pk_udca_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = Day, y = UDCA_plas_mgL)) +
      geom_line(color = "#1565C0", linewidth = 1.2) +
      labs(x = "Day", y = "UDCA Plasma (mg/L)", title = "UDCA Plasma Concentration") +
      theme_classic()
    ggplotly(p)
  })

  output$pk_bile_udca <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = Day, y = UDCA_bile_uM)) +
      geom_line(color = "#2E7D32", linewidth = 1.2) +
      geom_hline(yintercept = 300, linetype="dashed", color="orange") +
      annotate("text", x=max(d$Day)*0.7, y=310, label="Target >300 µmol/L",
               color="orange", size=3) +
      labs(x = "Day", y = "UDCA in Bile (µmol/L)") +
      theme_classic()
    ggplotly(p)
  })

  output$pk_stat_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = Day, y = STAT_plas_mgL)) +
      geom_line(color = "#6A1B9A", linewidth = 1.2) +
      labs(x = "Day", y = "Statin Plasma (mg/L)", title = "Statin PK") +
      theme_classic()
    ggplotly(p)
  })

  ## ---- Bile/CSI plots -----------------------------------------------------
  output$csi_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = Day, y = CSI_out)) +
      geom_hline(yintercept = 1.0, linetype="dashed", color="red", linewidth=0.8) +
      annotate("text", x=5, y=1.02, label="Lithogenic threshold (CSI=1)", color="red", size=3.5) +
      geom_line(color="#C62828", linewidth=1.3) +
      labs(x = "Day", y = "CSI", title = "Cholesterol Saturation Index") +
      theme_classic()
    ggplotly(p)
  })

  output$csi_val_box <- renderValueBox({
    d <- sim_data()
    last_csi <- tail(d$CSI_out, 1)
    color <- if(last_csi > 1.2) "red" else if(last_csi > 1.0) "yellow" else "green"
    valueBox(round(last_csi, 2), "Final CSI", icon = icon("flask"), color = color)
  })

  output$ba_pool_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = Day, y = BA_pool)) +
      geom_line(color = "#1B5E20", linewidth = 1.2) +
      labs(x = "Day", y = "BA Pool (g)") +
      theme_classic()
    ggplotly(p)
  })

  output$bil_comp_plot <- renderPlotly({
    d <- sim_data() %>%
      select(Day, CHOL_bil, PL_bil) %>%
      pivot_longer(c(CHOL_bil, PL_bil), names_to = "Component", values_to = "Conc")
    p <- ggplot(d, aes(x = Day, y = Conc, color = Component)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("CHOL_bil" = "#FF8F00", "PL_bil" = "#0277BD")) +
      labs(x = "Day", y = "mmol/L") +
      theme_classic()
    ggplotly(p)
  })

  ## ---- Stone plots --------------------------------------------------------
  output$stone_vol_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = Day, y = Stone_V)) +
      geom_area(fill = "#E3F2FD", alpha = 0.5) +
      geom_line(color = "#0D47A1", linewidth = 1.3) +
      labs(x = "Day", y = "Stone Volume (mL)", title = "Gallstone Volume") +
      theme_classic()
    ggplotly(p)
  })

  output$stone_dissolv_box <- renderValueBox({
    d <- sim_data()
    init_v <- d$Stone_V[1]
    last_v <- tail(d$Stone_V, 1)
    pct    <- if(init_v > 0) (1 - last_v/init_v) * 100 else 0
    color  <- if(pct > 50) "green" else if(pct > 20) "yellow" else "red"
    valueBox(paste0(round(pct, 1), "%"), "Stone Dissolved",
             icon = icon("check-circle"), color = color)
  })

  output$crystal_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = Day, y = Crystal_mass)) +
      geom_line(color = "#7B1FA2", linewidth = 1.1) +
      labs(x = "Day", y = "Crystal Mass (mg)") +
      theme_classic()
    ggplotly(p)
  })

  output$gb_vol_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = Day, y = GB_vol)) +
      geom_line(color = "#880E4F", linewidth = 1.1) +
      labs(x = "Day", y = "GB Volume (mL)") +
      theme_classic()
    ggplotly(p)
  })

  ## ---- Scenario Comparison ------------------------------------------------
  compare_data <- eventReactive(input$run_compare, {
    sc_ids   <- c(A = input$sc_a, B = input$sc_b, C = input$sc_c)
    sc_names <- c(
      "none"      = "No Treatment",
      "udca750"   = "UDCA 750 mg",
      "udca1050"  = "UDCA 1050 mg",
      "udca_stat" = "UDCA+Statin",
      "ezet"      = "Ezetimibe",
      "udca_wl"   = "UDCA+WtLoss"
    )
    bind_rows(lapply(names(sc_ids), function(nm) {
      p <- get_scenario_params(sc_ids[nm])
      run_sim(mod, dose_UDCA = p$DOSE_UDCA, dose_STAT = p$DOSE_STAT,
              dose_EZET = p$DOSE_EZET, bwt = input$bwt,
              stone_init = (input$stone_init_mm / 10)^3 * pi / 6,
              wloss = p$WLOSS, dur_days = input$dur) %>%
        mutate(Scenario = sc_names[sc_ids[nm]])
    }))
  })

  output$compare_stone <- renderPlotly({
    d <- compare_data()
    p <- ggplot(d, aes(x = Day, y = Stone_V, color = Scenario)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Day", y = "Stone Volume (mL)") +
      theme_classic(base_size = 11)
    ggplotly(p)
  })

  output$compare_csi <- renderPlotly({
    d <- compare_data()
    p <- ggplot(d, aes(x = Day, y = CSI_out, color = Scenario)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
      geom_line(linewidth = 1.1) +
      labs(x = "Day", y = "CSI") +
      theme_classic(base_size = 11)
    ggplotly(p)
  })

  output$compare_table <- renderDT({
    d <- compare_data() %>%
      filter(Day %in% c(0, 90, 180, 365)) %>%
      group_by(Scenario, Day) %>%
      summarise(
        `Stone Vol (mL)` = round(mean(Stone_V), 3),
        `Stone Diam (mm)`= round(mean(Stone_mm), 1),
        `CSI`            = round(mean(CSI_out), 2),
        `BA Pool (g)`    = round(mean(BA_pool), 2),
        `UDCA Bile (µM)` = round(mean(UDCA_bile_uM), 0),
        `CRP (mg/L)`     = round(mean(CRP_plas), 2),
        .groups = "drop"
      )
    datatable(d, options = list(pageLength = 12, scrollX = TRUE),
              rownames = FALSE)
  })

  ## ---- Biomarker plots ----------------------------------------------------
  output$inflam_plot <- renderPlotly({
    d <- sim_data() %>%
      select(Day, IL6, CRP_plas) %>%
      pivot_longer(c(IL6, CRP_plas), names_to = "Marker", values_to = "Value")
    p <- ggplot(d, aes(x = Day, y = Value, color = Marker)) +
      geom_line(linewidth = 1.1) +
      facet_wrap(~Marker, scales = "free_y") +
      scale_color_manual(values = c("IL6" = "#C62828", "CRP_plas" = "#E65100")) +
      labs(x = "Day", y = "Concentration") +
      theme_classic()
    ggplotly(p)
  })

  output$crp_box <- renderValueBox({
    d <- sim_data()
    crp <- round(tail(d$CRP_plas, 1), 1)
    color <- if(crp > 10) "red" else if(crp > 1) "yellow" else "green"
    valueBox(paste(crp, "mg/L"), "CRP (final)", icon = icon("vial"), color = color)
  })

  output$il6_box <- renderValueBox({
    d <- sim_data()
    il6 <- round(tail(d$IL6, 1), 1)
    color <- if(il6 > 10) "red" else if(il6 > 3) "yellow" else "green"
    valueBox(paste(il6, "pg/mL"), "IL-6 (final)", icon = icon("microscope"), color = color)
  })

  output$chol_hep_plot <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x = Day, y = CHOL_h)) +
      geom_line(color = "#FF8F00", linewidth = 1.1) +
      labs(x = "Day", y = "Hepatic Cholesterol (mmol)") +
      theme_classic()
    ggplotly(p)
  })

  output$lft_plot <- renderPlotly({
    d <- sim_data() %>%
      mutate(ALP_proxy = 80 + 15 * pmax(0, CSI_out - 1) * Stone_V * 2,
             GGT_proxy = 25 + 10 * pmax(0, CSI_out - 1)) %>%
      select(Day, ALP_proxy, GGT_proxy) %>%
      pivot_longer(c(ALP_proxy, GGT_proxy), names_to = "Enzyme", values_to = "U_L")
    p <- ggplot(d, aes(x = Day, y = U_L, color = Enzyme)) +
      geom_line(linewidth = 1.0) +
      labs(x = "Day", y = "U/L (modeled proxy)") +
      theme_classic()
    ggplotly(p)
  })

  output$dose_resp_plot <- renderPlotly({
    doses <- seq(250, 1500, by = 250)
    stone_init <- (input$stone_init_mm / 10)^3 * pi / 6
    dr <- lapply(doses, function(d) {
      out <- run_sim(mod, dose_UDCA = d, dose_STAT = 0, dose_EZET = 0,
                     bwt = input$bwt, stone_init = stone_init,
                     wloss = 0, dur_days = 365)
      last <- tail(out, 1)
      pct  <- if(stone_init > 0) (1 - last$Stone_V / stone_init) * 100 else 0
      data.frame(Dose = d, Pct_Dissolved = max(0, pct))
    }) %>% bind_rows()

    p <- ggplot(dr, aes(x = Dose, y = Pct_Dissolved)) +
      geom_point(size = 3, color = "#1565C0") +
      geom_line(color = "#1565C0", linewidth = 1.1) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
      labs(x = "UDCA Daily Dose (mg/day)",
           y = "Stone Dissolution at 12 months (%)") +
      theme_classic()
    ggplotly(p)
  })
}

## ---- Launch -----------------------------------------------------------------
shinyApp(ui, server)
