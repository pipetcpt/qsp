################################################################################
# Pneumoconiosis QSP — Interactive Shiny Dashboard
# 진폐증 정량적 시스템 약리학 대시보드
#
# Tabs:
#   1. 환자 프로파일 (Patient Profile & Exposure History)
#   2. PK 시뮬레이션 (Drug PK — NAC, Pirfenidone, Nintedanib, Tetrandrine)
#   3. PD 주요 지표 (Cytokines, Oxidative Stress, Fibrotic Signaling)
#   4. 임상 엔드포인트 (FVC, DLCO, mPAP, KL-6, Dyspnea)
#   5. 시나리오 비교 (Multi-scenario Comparison)
#   6. 바이오마커 & 가상환자 (Biomarkers & Virtual Population)
################################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(DT)

# ============================================================
# QSP MODEL DEFINITION (inline, same as R model file)
# ============================================================
pnm_model <- mrgsolve::mcode("pnm_shiny", '
$PARAM
D_in=0.5, k_clr=0.05, k_phag=0.10, D_overload=20.0,
AM0_base=100.0, k_AM_act=0.04, k_AM_rec=0.02, k_AM_pyr=0.01, k_mono=0.03,
k_NLRP3=0.08, k_NLRP3_d=0.15,
k_IL1b_syn=5.0, k_IL1b_deg=0.3,
k_TNFa_syn=3.0, k_TNFa_deg=0.4,
k_TGFb_syn=2.0, k_TGFb_deg=0.20,
k_IL10_syn=1.0, k_IL10_deg=0.25,
k_ROS_syn=0.12, k_ROS_deg=0.20,
GSH_base=10.0, k_GSH_syn=0.08, k_GSH_dep=0.06,
k_Neut_rec=0.15, k_Neut_d=0.20,
k_Fibro_syn=0.05, k_Fibro_deg=0.03,
k_MyoFib=0.08, k_MyoFib_d=0.02,
k_Coll_syn=0.10, k_Coll_deg=0.01,
FVC_0=100.0, FVC_min=30.0, k_FVC_loss=0.008, k_PVR_inc=0.005, PVR_0=1.0,
NAC_ka=1.5, NAC_CL=15.0, NAC_Vd=25.0, NAC_Emax=0.80, NAC_EC50=5.0,
Pirf_ka=2.0, Pirf_CL=3.5, Pirf_Vd=70.0, Pirf_Emax=0.75, Pirf_EC50=2.0,
Nint_ka=0.8, Nint_CL=1.3, Nint_Vd=640.0, Nint_Emax=0.70, Nint_EC50=0.05,
Tetra_ka=1.2, Tetra_CL=8.0, Tetra_Vd=120.0, Tetra_Emax=0.65, Tetra_EC50=0.3,
ALAG_NAC=0, ALAG_Pirf=0, ALAG_Nint=0, ALAG_Tetra=0

$INIT
D_alv=0, D_clr=0, AM_rest=100, AM_act=0, NLRP3=0,
IL1b=0.5, TNFa=1.0, TGFb=2.0, IL10=2.0,
ROS=1.0, GSH=10.0, Neutro=5.0,
Fibro=1.0, MyoFib=0.0, Coll=1.0,
FVC=100.0, PVR=1.0,
C_NAC=0, C_Pirf=0, C_Nint=0, C_Tetra=0

$ODE
double overload_factor = 1.0 / (1.0 + D_alv / D_overload);
dxdt_D_alv = D_in - k_clr*D_alv*overload_factor - k_phag*AM_rest*D_alv/(D_alv+5.0);
dxdt_D_clr = k_clr*D_alv*overload_factor + k_phag*AM_rest*D_alv/(D_alv+5.0);
double Pirf_inh_TGFb = Pirf_Emax*C_Pirf/(Pirf_EC50+C_Pirf);
double Nint_inh = Nint_Emax*C_Nint/(Nint_EC50+C_Nint);
double NAC_eff = NAC_Emax*C_NAC/(NAC_EC50+C_NAC);
double Tetra_inh = Tetra_Emax*C_Tetra/(Tetra_EC50+C_Tetra);
double D_effect = k_AM_act*D_alv;
dxdt_AM_rest = k_mono*10.0 - D_effect*AM_rest - k_AM_pyr*AM_act;
dxdt_AM_act  = D_effect*AM_rest - k_AM_rec*AM_act - k_AM_pyr*AM_act;
double NLRP3_drive = k_NLRP3*AM_act*ROS/(1.0+IL10/5.0)*(1.0-Tetra_inh);
dxdt_NLRP3 = NLRP3_drive - k_NLRP3_d*NLRP3;
dxdt_IL1b = k_IL1b_syn*NLRP3 - k_IL1b_deg*IL1b;
dxdt_TNFa = k_TNFa_syn*AM_act*(1.0+0.5*IL1b/10.0)*(1.0-Tetra_inh) - k_TNFa_deg*TNFa;
dxdt_TGFb = k_TGFb_syn*(AM_act+0.3*IL1b/5.0)*(1.0-Pirf_inh_TGFb)*(1.0-0.3*Nint_inh) - k_TGFb_deg*TGFb;
dxdt_IL10 = k_IL10_syn*(AM_act*0.5+2.0) - k_IL10_deg*IL10;
dxdt_ROS = k_ROS_syn*(AM_act+0.3*Neutro) - k_ROS_deg*ROS*(1.0+GSH/GSH_base) - 0.5*NAC_eff*ROS;
dxdt_GSH = k_GSH_syn*GSH_base*(1.0+NAC_eff) - k_GSH_dep*ROS*GSH - 0.02*GSH;
dxdt_Neutro = k_Neut_rec*(IL1b+TNFa)/10.0 - k_Neut_d*Neutro;
double Pirf_inh_F = Pirf_Emax*0.5*C_Pirf/(Pirf_EC50+C_Pirf);
double Nint_inh_F = Nint_Emax*0.4*C_Nint/(Nint_EC50+C_Nint);
dxdt_Fibro = k_Fibro_syn*TGFb*(1.0-Pirf_inh_F-Nint_inh_F) - k_Fibro_deg*Fibro;
dxdt_MyoFib = k_MyoFib*TGFb*Fibro*(1.0-Pirf_inh_TGFb)*(1.0-0.5*Nint_inh) - k_MyoFib_d*MyoFib;
dxdt_Coll = k_Coll_syn*MyoFib*(1.0-0.6*Pirf_inh_TGFb) - k_Coll_deg*Coll;
double FVC_target = FVC_0 - k_FVC_loss*(Coll-1.0)*100.0;
if(FVC_target < FVC_min) FVC_target = FVC_min;
dxdt_FVC = 0.1*(FVC_target - FVC);
dxdt_PVR = k_PVR_inc*(Coll-1.0) + 0.02*(100.0-FVC)/100.0 - 0.01*(PVR-1.0);
dxdt_C_NAC   = NAC_ka*ALAG_NAC     - (NAC_CL/NAC_Vd)*C_NAC;
dxdt_C_Pirf  = Pirf_ka*ALAG_Pirf   - (Pirf_CL/Pirf_Vd)*C_Pirf;
dxdt_C_Nint  = Nint_ka*ALAG_Nint   - (Nint_CL/Nint_Vd)*C_Nint;
dxdt_C_Tetra = Tetra_ka*ALAG_Tetra - (Tetra_CL/Tetra_Vd)*C_Tetra;

$CAPTURE
double FEV1=FVC*0.80;
double DLCO=100.0-0.7*(100.0-FVC);
double mPAP=15.0+(PVR-1.0)*12.0;
double KL6=200.0+Coll*150.0;
double dyspnea=(100.0-FVC)/14.0;
double Pirf_inh_TGFb=Pirf_Emax*C_Pirf/(Pirf_EC50+C_Pirf);
double Nint_inh=Nint_Emax*C_Nint/(Nint_EC50+C_Nint);
double NAC_eff=NAC_Emax*C_NAC/(NAC_EC50+C_NAC);
', quiet = TRUE)

# Helper to run simulation
run_sim <- function(mod, params_list, drug_dose_nac, drug_dose_pirf, drug_dose_nint,
                    drug_dose_tetra, sim_end = 120) {
  mod <- do.call(param, c(list(mod), params_list))

  events <- ev(time = 0, amt = 0, cmt = 1)  # dummy
  if (drug_dose_nac > 0)
    events <- events + ev(amt = drug_dose_nac, cmt = "C_NAC",   time = 0, ii = 0.333, addl = sim_end * 3)
  if (drug_dose_pirf > 0)
    events <- events + ev(amt = drug_dose_pirf, cmt = "C_Pirf", time = 0, ii = 0.333, addl = sim_end * 3)
  if (drug_dose_nint > 0)
    events <- events + ev(amt = drug_dose_nint, cmt = "C_Nint", time = 0, ii = 0.5,   addl = sim_end * 2)
  if (drug_dose_tetra > 0)
    events <- events + ev(amt = drug_dose_tetra, cmt = "C_Tetra", time = 0, ii = 0.5, addl = sim_end * 2)

  mrgsim(mod, end = sim_end, delta = 0.5, events = events) %>% as_tibble()
}

# ============================================================
# UI DEFINITION
# ============================================================
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = "Pneumoconiosis QSP",
    titleWidth = 250
  ),

  dashboardSidebar(
    width = 250,
    sidebarMenu(
      menuItem("Patient Profile",       tabName = "tab_patient",   icon = icon("user-md")),
      menuItem("Drug PK",               tabName = "tab_pk",        icon = icon("pills")),
      menuItem("PD Biomarkers",         tabName = "tab_pd",        icon = icon("flask")),
      menuItem("Clinical Endpoints",    tabName = "tab_clinical",  icon = icon("lungs")),
      menuItem("Scenario Comparison",   tabName = "tab_scenario",  icon = icon("chart-line")),
      menuItem("Virtual Population",    tabName = "tab_vp",        icon = icon("users"))
    ),

    hr(),
    h5("Simulation Settings", style = "color:white; padding-left:10px;"),
    sliderInput("sim_end",  "Duration (months):", 12, 120, 60, step = 12),
    hr(),
    h5("Dust Exposure", style = "color:white; padding-left:10px;"),
    sliderInput("D_in",     "Dust Input Rate (mg/mo):", 0, 3, 0.5, step = 0.1),
    selectInput("dust_type", "Dust Type:",
                choices = c("Silica (Silicosis)"    = "silica",
                            "Coal Dust (CWP)"       = "coal",
                            "Asbestos"              = "asbestos",
                            "Mixed Dust"            = "mixed"),
                selected = "silica")
  ),

  dashboardBody(

    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f5f7fa; }
      .box { border-top-color: #2c3e50 !important; }
      .value-box { min-height: 90px; }
    "))),

    tabItems(

      # ----------------------------------------------------------
      # TAB 1: Patient Profile
      # ----------------------------------------------------------
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Characteristics", width = 4, status = "primary",
            sliderInput("age",    "Age (years):", 30, 70, 50),
            selectInput("sex",    "Sex:", choices = c("Male", "Female"), selected = "Male"),
            sliderInput("smoking","Smoking Pack-Years:", 0, 60, 20),
            sliderInput("work_yrs","Years of Occupational Exposure:", 1, 40, 15),
            selectInput("pneumo_type", "Pneumoconiosis Type:",
                        choices = c("Simple Silicosis (Category 1-2)" = "simple",
                                    "Complicated Silicosis / PMF" = "pmf",
                                    "Coal Workers' Pneumoconiosis" = "cwp",
                                    "Asbestosis" = "asbestosis"),
                        selected = "simple"),
            selectInput("stage", "Disease Stage:",
                        choices = c("Early (0-5 yr exposure)" = "early",
                                    "Intermediate (5-15 yr)" = "intermediate",
                                    "Advanced (>15 yr)" = "advanced"),
                        selected = "intermediate")
          ),
          box(title = "Baseline Clinical Status", width = 4, status = "warning",
            sliderInput("FVC_0",  "Baseline FVC (% predicted):", 50, 120, 90),
            sliderInput("DLCO_0", "Baseline DLCO (% predicted):", 40, 120, 85),
            sliderInput("KL6_0",  "Baseline KL-6 (U/mL):", 100, 1000, 250),
            sliderInput("mPAP_0", "Baseline mPAP (mmHg):", 10, 40, 17),
            checkboxInput("sil_tb", "Comorbid Silicotuberculosis", FALSE),
            checkboxInput("lung_ca","Comorbid Lung Cancer Risk", FALSE)
          ),
          box(title = "Patient Summary", width = 4, status = "info",
            valueBoxOutput("vbox_stage", width = 12),
            valueBoxOutput("vbox_fvc",   width = 12),
            valueBoxOutput("vbox_risk",  width = 12),
            br(),
            htmlOutput("patient_summary")
          )
        ),
        fluidRow(
          box(title = "Exposure Timeline & Disease Progression Conceptual Model",
              width = 12, status = "primary",
              p(strong("Pneumoconiosis Pathophysiology (Conceptual):"), style = "font-size:14px"),
              tags$ol(
                tags$li("Mineral dust (silica/coal/asbestos) inhaled → deposited in alveoli"),
                tags$li("Alveolar macrophages attempt phagocytosis → NLRP3 inflammasome activation"),
                tags$li("IL-1β, TNF-α, TGF-β1 released → sustained inflammatory cascade"),
                tags$li("Oxidative stress (ROS) amplifies macrophage injury and cytokine signaling"),
                tags$li("TGF-β1 → fibroblast → myofibroblast differentiation → collagen deposition"),
                tags$li("Progressive nodular fibrosis → PMF (coalescence of nodules)"),
                tags$li("Restrictive lung disease: FVC↓, DLCO↓, TLC↓"),
                tags$li("Pulmonary hypertension → cor pulmonale → death")
              ),
              br(),
              imageOutput("mech_img", height = "300px")
          )
        )
      ),

      # ----------------------------------------------------------
      # TAB 2: Drug PK
      # ----------------------------------------------------------
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Drug Dosing", width = 3, status = "primary",
            h4("N-Acetylcysteine (NAC)"),
            checkboxInput("use_nac", "Use NAC", FALSE),
            sliderInput("nac_dose", "NAC Dose (mg, per administration):", 200, 1800, 600, step = 200),
            selectInput("nac_freq", "Dosing Frequency:", c("TID" = "tid", "BID" = "bid", "QD" = "qd")),
            hr(),
            h4("Pirfenidone"),
            checkboxInput("use_pirf", "Use Pirfenidone", FALSE),
            sliderInput("pirf_dose", "Pirfenidone Dose (mg, per dose):", 267, 1002, 801, step = 267),
            selectInput("pirf_freq", "Dosing Frequency:", c("TID" = "tid", "BID" = "bid")),
            hr(),
            h4("Nintedanib"),
            checkboxInput("use_nint", "Use Nintedanib", FALSE),
            sliderInput("nint_dose", "Nintedanib Dose (mg, per dose):", 100, 200, 150, step = 50),
            hr(),
            h4("Tetrandrine"),
            checkboxInput("use_tetra", "Use Tetrandrine", FALSE),
            sliderInput("tetra_dose", "Tetrandrine Dose (mg, per dose):", 60, 120, 60, step = 20)
          ),
          box(title = "PK Profiles (Plasma Concentrations)", width = 9, status = "success",
            plotlyOutput("pk_plot", height = "500px"),
            hr(),
            DTOutput("pk_table")
          )
        )
      ),

      # ----------------------------------------------------------
      # TAB 3: PD Biomarkers
      # ----------------------------------------------------------
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(title = "PD Parameter Adjustments", width = 3, status = "warning",
            h4("Inflammasome"),
            sliderInput("k_NLRP3",    "NLRP3 Activation Rate:", 0.01, 0.3,  0.08, step = 0.01),
            sliderInput("k_IL1b_syn", "IL-1β Synthesis Rate:", 0.5, 15.0, 5.0, step = 0.5),
            hr(),
            h4("Fibrotic Signaling"),
            sliderInput("k_TGFb_syn", "TGF-β1 Synthesis Rate:", 0.2, 8.0, 2.0, step = 0.2),
            sliderInput("k_Coll_syn", "Collagen Synthesis Rate:", 0.02, 0.4, 0.10, step = 0.02),
            hr(),
            h4("Oxidative Stress"),
            sliderInput("k_ROS_syn",  "ROS Production Rate:", 0.02, 0.4, 0.12, step = 0.02),
            sliderInput("k_GSH_syn",  "GSH Synthesis Rate:", 0.01, 0.3, 0.08, step = 0.01)
          ),
          box(title = "Cytokine & Inflammasome Dynamics", width = 9, status = "danger",
            plotlyOutput("pd_cytokines", height = "350px"),
            hr(),
            plotlyOutput("pd_ros_gsh", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Macrophage & Cellular Dynamics", width = 6, status = "primary",
            plotlyOutput("pd_mac", height = "300px")
          ),
          box(title = "Fibroblast → Myofibroblast → Collagen Cascade", width = 6, status = "warning",
            plotlyOutput("pd_fibro", height = "300px")
          )
        )
      ),

      # ----------------------------------------------------------
      # TAB 4: Clinical Endpoints
      # ----------------------------------------------------------
      tabItem(tabName = "tab_clinical",
        fluidRow(
          valueBoxOutput("ep_fvc",     width = 3),
          valueBoxOutput("ep_dlco",    width = 3),
          valueBoxOutput("ep_mpap",    width = 3),
          valueBoxOutput("ep_kl6",     width = 3)
        ),
        fluidRow(
          box(title = "Pulmonary Function Decline (FVC, FEV1, DLCO)", width = 6, status = "danger",
            plotlyOutput("ep_lung_fn", height = "350px")
          ),
          box(title = "Pulmonary Hypertension (mPAP & PVR)", width = 6, status = "warning",
            plotlyOutput("ep_pah", height = "350px")
          )
        ),
        fluidRow(
          box(title = "KL-6 Biomarker & Dyspnea Score", width = 6, status = "info",
            plotlyOutput("ep_kl6_plot", height = "300px")
          ),
          box(title = "Clinical Outcomes Table", width = 6, status = "primary",
            DTOutput("ep_table")
          )
        )
      ),

      # ----------------------------------------------------------
      # TAB 5: Scenario Comparison
      # ----------------------------------------------------------
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "Scenario Selection", width = 3, status = "primary",
            checkboxGroupInput("scenarios", "Compare Scenarios:",
              choices = c(
                "1. No Treatment"               = "sc1",
                "2. NAC Monotherapy"            = "sc2",
                "3. Pirfenidone Monotherapy"    = "sc3",
                "4. Nintedanib Monotherapy"     = "sc4",
                "5. Dust Cessation + NAC"       = "sc5",
                "6. Pirfenidone + NAC Combo"    = "sc6"
              ),
              selected = c("sc1", "sc3", "sc4", "sc6")
            ),
            hr(),
            selectInput("sc_endpoint", "Primary Endpoint:",
                        choices = c("FVC" = "FVC", "DLCO" = "DLCO",
                                    "Collagen" = "Coll", "TGF-β" = "TGFb",
                                    "KL-6" = "KL6", "mPAP" = "mPAP",
                                    "ROS" = "ROS"),
                        selected = "FVC"),
            actionButton("run_sc", "Run Comparison", class = "btn-primary btn-block", icon = icon("play"))
          ),
          box(title = "Multi-Scenario Endpoint Comparison", width = 9, status = "success",
            plotlyOutput("sc_plot", height = "450px"),
            hr(),
            DTOutput("sc_summary_table")
          )
        )
      ),

      # ----------------------------------------------------------
      # TAB 6: Virtual Population
      # ----------------------------------------------------------
      tabItem(tabName = "tab_vp",
        fluidRow(
          box(title = "Virtual Population Settings", width = 3, status = "primary",
            sliderInput("n_vp", "Number of Patients:", 20, 200, 50, step = 10),
            sliderInput("cv_pk", "PK Variability (CV%):", 10, 60, 30, step = 5),
            sliderInput("cv_pd", "PD Variability (CV%):", 10, 60, 30, step = 5),
            selectInput("vp_treatment", "Treatment:",
                        choices = c("No Treatment" = "none",
                                    "Pirfenidone" = "pirf",
                                    "Nintedanib" = "nint")),
            actionButton("run_vp", "Run VP Simulation", class = "btn-success btn-block", icon = icon("play")),
            hr(),
            downloadButton("dl_vp", "Download VP Results", class = "btn-info btn-block")
          ),
          box(title = "Virtual Population: FVC Distribution", width = 9, status = "info",
            plotlyOutput("vp_fvc_plot", height = "400px"),
            hr(),
            plotlyOutput("vp_hist_plot", height = "250px")
          )
        ),
        fluidRow(
          box(title = "Population Responder Analysis", width = 6, status = "warning",
            plotlyOutput("vp_responder", height = "300px")
          ),
          box(title = "Biomarker Distribution at 5 Years", width = 6, status = "primary",
            plotlyOutput("vp_biomarker", height = "300px")
          )
        )
      )

    ) # end tabItems
  ) # end dashboardBody
) # end dashboardPage

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # --- Reactive: Base simulation ---
  base_sim <- reactive({
    params <- list(
      D_in       = input$D_in,
      k_NLRP3    = if (!is.null(input$k_NLRP3)) input$k_NLRP3 else 0.08,
      k_IL1b_syn = if (!is.null(input$k_IL1b_syn)) input$k_IL1b_syn else 5.0,
      k_TGFb_syn = if (!is.null(input$k_TGFb_syn)) input$k_TGFb_syn else 2.0,
      k_Coll_syn = if (!is.null(input$k_Coll_syn)) input$k_Coll_syn else 0.10,
      k_ROS_syn  = if (!is.null(input$k_ROS_syn)) input$k_ROS_syn else 0.12,
      k_GSH_syn  = if (!is.null(input$k_GSH_syn)) input$k_GSH_syn else 0.08,
      FVC_0      = if (!is.null(input$FVC_0)) input$FVC_0 else 90.0
    )
    dose_nac  <- if (isTRUE(input$use_nac))  input$nac_dose  * 0.044 else 0
    dose_pirf <- if (isTRUE(input$use_pirf)) input$pirf_dose * 0.025 else 0
    dose_nint <- if (isTRUE(input$use_nint)) input$nint_dose * 0.0017 else 0
    dose_tetra<- if (isTRUE(input$use_tetra)) input$tetra_dose * 0.012 else 0
    run_sim(pnm_model, params, dose_nac, dose_pirf, dose_nint, dose_tetra, input$sim_end)
  })

  # --- Patient Profile Value Boxes ---
  output$vbox_stage <- renderValueBox({
    valueBox(value = input$stage, subtitle = "Disease Stage", icon = icon("stethoscope"), color = "blue")
  })
  output$vbox_fvc <- renderValueBox({
    col <- if (input$FVC_0 >= 80) "green" else if (input$FVC_0 >= 60) "yellow" else "red"
    valueBox(value = paste0(input$FVC_0, "%"), subtitle = "Baseline FVC", icon = icon("lungs"), color = col)
  })
  output$vbox_risk <- renderValueBox({
    risk <- if (input$work_yrs > 20) "High" else if (input$work_yrs > 10) "Moderate" else "Low"
    col  <- if (risk == "High") "red" else if (risk == "Moderate") "yellow" else "green"
    valueBox(value = risk, subtitle = "Progression Risk", icon = icon("exclamation-triangle"), color = col)
  })
  output$patient_summary <- renderUI({
    HTML(paste0(
      "<b>Age:</b> ", input$age, " yr | <b>Sex:</b> ", input$sex, "<br>",
      "<b>Smoking:</b> ", input$smoking, " pack-yr<br>",
      "<b>Occupational exposure:</b> ", input$work_yrs, " yr<br>",
      "<b>Type:</b> ", input$pneumo_type, "<br>",
      "<b>Silicotuberculosis:</b> ", if (isTRUE(input$sil_tb)) "YES ⚠️" else "No", "<br>",
      "<b>Lung Cancer Risk:</b> ", if (isTRUE(input$lung_ca)) "Elevated ⚠️" else "Standard"
    ))
  })

  # --- PK Plot ---
  output$pk_plot <- renderPlotly({
    sim <- base_sim()
    pk_long <- sim %>%
      select(time, C_NAC, C_Pirf, C_Nint, C_Tetra) %>%
      pivot_longer(-time, names_to = "Drug", values_to = "Conc") %>%
      filter(Conc > 0.001)
    if (nrow(pk_long) == 0) {
      return(plotly_empty() %>% layout(title = "No drug selected — enable a drug to see PK"))
    }
    p <- ggplot(pk_long, aes(x = time, y = Conc, color = Drug)) +
      geom_line(linewidth = 1.0) +
      labs(x = "Time (months)", y = "Plasma Concentration (µg/mL)",
           title = "Drug Plasma Concentration Profiles") +
      theme_bw()
    ggplotly(p)
  })
  output$pk_table <- renderDT({
    sim <- base_sim()
    sim %>%
      filter(time %in% c(0, 1, 3, 6, 12, 24, 48, input$sim_end)) %>%
      select(time, C_NAC, C_Pirf, C_Nint, C_Tetra) %>%
      mutate(across(-time, ~round(., 4))) %>%
      datatable(rownames = FALSE, options = list(pageLength = 8))
  })

  # --- PD Cytokine Plot ---
  output$pd_cytokines <- renderPlotly({
    sim <- base_sim()
    p <- sim %>%
      select(time, IL1b, TNFa, TGFb, IL10) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x = time, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Time (months)", y = "Concentration (pg/mL)",
           title = "Cytokine Dynamics", color = "Cytokine") +
      theme_bw()
    ggplotly(p)
  })
  output$pd_ros_gsh <- renderPlotly({
    sim <- base_sim()
    p <- sim %>%
      select(time, ROS, GSH) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x = time, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Time (months)", y = "Normalized value",
           title = "Oxidative Stress: ROS vs GSH", color = "Variable") +
      theme_bw()
    ggplotly(p)
  })
  output$pd_mac <- renderPlotly({
    sim <- base_sim()
    p <- sim %>%
      select(time, AM_rest, AM_act, Neutro) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x = time, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Time (months)", y = "Cell index",
           title = "Macrophage & Neutrophil Dynamics", color = "Cell") +
      theme_bw()
    ggplotly(p)
  })
  output$pd_fibro <- renderPlotly({
    sim <- base_sim()
    p <- sim %>%
      select(time, Fibro, MyoFib, Coll) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x = time, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      labs(x = "Time (months)", y = "Index (normalized)",
           title = "Fibrosis Cascade: Fibroblast → Collagen", color = "Component") +
      theme_bw()
    ggplotly(p)
  })

  # --- Clinical Endpoints ---
  output$ep_fvc  <- renderValueBox({
    sim <- base_sim()
    val <- round(sim$FVC[nrow(sim)], 1)
    col <- if (val >= 80) "green" else if (val >= 60) "yellow" else "red"
    valueBox(paste0(val, "%"), "FVC at End", icon = icon("lungs"), color = col)
  })
  output$ep_dlco <- renderValueBox({
    sim <- base_sim()
    val <- round(sim$DLCO[nrow(sim)], 1)
    col <- if (val >= 75) "green" else if (val >= 55) "yellow" else "red"
    valueBox(paste0(val, "%"), "DLCO at End", icon = icon("wind"), color = col)
  })
  output$ep_mpap <- renderValueBox({
    sim <- base_sim()
    val <- round(sim$mPAP[nrow(sim)], 1)
    col <- if (val < 25) "green" else if (val < 35) "yellow" else "red"
    valueBox(paste0(val, " mmHg"), "mPAP at End", icon = icon("heart"), color = col)
  })
  output$ep_kl6  <- renderValueBox({
    sim <- base_sim()
    val <- round(sim$KL6[nrow(sim)], 0)
    col <- if (val < 500) "green" else if (val < 1000) "yellow" else "red"
    valueBox(paste0(val, " U/mL"), "KL-6 at End", icon = icon("vial"), color = col)
  })
  output$ep_lung_fn <- renderPlotly({
    sim <- base_sim()
    p <- sim %>%
      select(time, FVC, FEV1, DLCO) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x = time, y = value, color = name)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 70, linetype = "dashed", color = "red", alpha = 0.5) +
      labs(x = "Time (months)", y = "% predicted",
           title = "Pulmonary Function Decline", color = "Measure") +
      theme_bw()
    ggplotly(p)
  })
  output$ep_pah <- renderPlotly({
    sim <- base_sim()
    p <- sim %>%
      select(time, mPAP, PVR) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x = time, y = value, color = name)) +
      geom_line(linewidth = 1.2) +
      geom_hline(yintercept = 25, linetype = "dashed", color = "red", alpha = 0.6) +
      labs(x = "Time (months)", y = "Value",
           title = "Pulmonary Hypertension Risk", color = "Variable") +
      theme_bw()
    ggplotly(p)
  })
  output$ep_kl6_plot <- renderPlotly({
    sim <- base_sim()
    p <- sim %>%
      select(time, KL6, dyspnea) %>%
      pivot_longer(-time) %>%
      ggplot(aes(x = time, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      facet_wrap(~name, scales = "free_y", nrow = 1) +
      labs(x = "Time (months)", y = "Value", title = "KL-6 & Dyspnea Score") +
      theme_bw() +
      theme(legend.position = "none")
    ggplotly(p)
  })
  output$ep_table <- renderDT({
    sim <- base_sim()
    sim %>%
      filter(time %in% c(0, 6, 12, 24, 36, 48, 60, 84, 120) & time <= input$sim_end) %>%
      select(time, FVC, FEV1, DLCO, mPAP, KL6, dyspnea, Coll) %>%
      mutate(across(-time, ~round(., 2))) %>%
      datatable(rownames = FALSE, caption = "Clinical Outcomes Over Time",
                options = list(pageLength = 10))
  })

  # --- Scenario Comparison ---
  scenario_data <- eventReactive(input$run_sc, {
    base_params <- list(D_in = input$D_in, FVC_0 = input$FVC_0 %||% 90)

    sc_list <- list()
    if ("sc1" %in% input$scenarios)
      sc_list[["1. No Treatment"]] <- run_sim(pnm_model, base_params, 0, 0, 0, 0, input$sim_end)
    if ("sc2" %in% input$scenarios)
      sc_list[["2. NAC Mono"]] <- run_sim(pnm_model, base_params, 600*0.044, 0, 0, 0, input$sim_end)
    if ("sc3" %in% input$scenarios)
      sc_list[["3. Pirfenidone"]] <- run_sim(pnm_model, base_params, 0, 801*0.025, 0, 0, input$sim_end)
    if ("sc4" %in% input$scenarios)
      sc_list[["4. Nintedanib"]] <- run_sim(pnm_model, base_params, 0, 0, 150*0.0017, 0, input$sim_end)
    if ("sc5" %in% input$scenarios)
      sc_list[["5. Dust Cessation+NAC"]] <- run_sim(
        param(pnm_model, D_in = 0), list(FVC_0 = input$FVC_0 %||% 90),
        600*0.044, 0, 0, 0, input$sim_end)
    if ("sc6" %in% input$scenarios)
      sc_list[["6. Pirf+NAC Combo"]] <- run_sim(pnm_model, base_params, 600*0.044, 801*0.025, 0, 0, input$sim_end)

    bind_rows(sc_list, .id = "Scenario")
  }, ignoreNULL = FALSE)

  output$sc_plot <- renderPlotly({
    df <- scenario_data()
    ep <- input$sc_endpoint
    if (is.null(df) || nrow(df) == 0) return(plotly_empty())
    p <- ggplot(df, aes_string(x = "time", y = ep, color = "Scenario")) +
      geom_line(linewidth = 1.1) +
      labs(x = "Time (months)", y = ep, title = paste("Scenario Comparison:", ep)) +
      theme_bw()
    ggplotly(p)
  })
  output$sc_summary_table <- renderDT({
    df <- scenario_data()
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>%
      group_by(Scenario) %>%
      filter(time == max(time)) %>%
      select(Scenario, FVC, DLCO, Coll, TGFb, KL6, mPAP) %>%
      mutate(across(where(is.numeric), ~round(., 2))) %>%
      datatable(rownames = FALSE, options = list(dom = "t"))
  })

  # --- Virtual Population ---
  vp_sim <- eventReactive(input$run_vp, {
    n <- input$n_vp
    cv_pd <- input$cv_pd / 100
    set.seed(123)
    pop <- tibble(
      ID = 1:n,
      k_TGFb_syn = rlnorm(n, log(2.0), cv_pd),
      k_Coll_syn = rlnorm(n, log(0.10), cv_pd),
      k_FVC_loss = rlnorm(n, log(0.008), cv_pd),
      D_in       = rlnorm(n, log(input$D_in), cv_pd),
      FVC_0      = rnorm(n, input$FVC_0 %||% 90, 8)
    )
    dose_pirf <- if (input$vp_treatment == "pirf") 801*0.025 else 0
    dose_nint <- if (input$vp_treatment == "nint") 150*0.0017 else 0

    purrr::map_dfr(split(pop, 1:nrow(pop)), function(row) {
      mod <- pnm_model %>%
        param(k_TGFb_syn = row$k_TGFb_syn, k_Coll_syn = row$k_Coll_syn,
              k_FVC_loss = row$k_FVC_loss, D_in = row$D_in) %>%
        init(FVC = max(row$FVC_0, 50))
      run_sim(mod, list(), 0, dose_pirf, dose_nint, 0, input$sim_end) %>%
        mutate(ID = row$ID)
    })
  })

  output$vp_fvc_plot <- renderPlotly({
    df <- vp_sim()
    if (is.null(df)) return(plotly_empty())
    summary_df <- df %>%
      group_by(time) %>%
      summarise(med = median(FVC), q5 = quantile(FVC, 0.05), q95 = quantile(FVC, 0.95))
    p <- ggplot(summary_df, aes(x = time)) +
      geom_ribbon(aes(ymin = q5, ymax = q95), fill = "steelblue", alpha = 0.3) +
      geom_line(aes(y = med), color = "steelblue", linewidth = 1.5) +
      geom_hline(yintercept = 70, linetype = "dashed", color = "red") +
      labs(x = "Time (months)", y = "FVC (% predicted)",
           title = paste0("Virtual Population (n=", input$n_vp, "): FVC Median ± 90% PI")) +
      theme_bw()
    ggplotly(p)
  })
  output$vp_hist_plot <- renderPlotly({
    df <- vp_sim()
    if (is.null(df)) return(plotly_empty())
    end_df <- df %>% filter(time == max(time))
    p <- ggplot(end_df, aes(x = FVC)) +
      geom_histogram(fill = "steelblue", color = "white", bins = 20) +
      geom_vline(xintercept = 70, color = "red", linetype = "dashed") +
      labs(x = "FVC at End (%)", y = "Count",
           title = "Distribution of FVC at End of Simulation") +
      theme_bw()
    ggplotly(p)
  })
  output$vp_responder <- renderPlotly({
    df <- vp_sim()
    if (is.null(df)) return(plotly_empty())
    responders <- df %>%
      filter(time == max(time)) %>%
      summarise(
        n_normal    = sum(FVC >= 80),
        n_mild      = sum(FVC >= 60 & FVC < 80),
        n_moderate  = sum(FVC >= 40 & FVC < 60),
        n_severe    = sum(FVC < 40)
      ) %>%
      pivot_longer(everything(), names_to = "Category", values_to = "Count") %>%
      mutate(Category = gsub("n_", "", Category))
    p <- ggplot(responders, aes(x = Category, y = Count, fill = Category)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("normal" = "#27AE60", "mild" = "#F39C12",
                                    "moderate" = "#E67E22", "severe" = "#E74C3C")) +
      labs(title = "FVC Category at End of Simulation", x = "FVC Category", y = "Patients") +
      theme_bw() + theme(legend.position = "none")
    ggplotly(p)
  })
  output$vp_biomarker <- renderPlotly({
    df <- vp_sim()
    if (is.null(df)) return(plotly_empty())
    yr5 <- df %>% filter(abs(time - min(60, input$sim_end)) < 1)
    p <- ggplot(yr5, aes(x = KL6)) +
      geom_density(fill = "tomato", alpha = 0.5) +
      geom_vline(xintercept = 500, color = "red", linetype = "dashed") +
      labs(title = "KL-6 Distribution at 5 Years", x = "KL-6 (U/mL)", y = "Density") +
      theme_bw()
    ggplotly(p)
  })
  output$dl_vp <- downloadHandler(
    filename = function() paste0("pneumoconiosis_vp_sim_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(vp_sim(), file, row.names = FALSE)
  )
}

# Run App
shinyApp(ui = ui, server = server)
