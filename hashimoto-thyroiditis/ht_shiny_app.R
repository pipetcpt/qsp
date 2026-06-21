##############################################################################
# Hashimoto's Thyroiditis — Interactive QSP Shiny Dashboard
# 7 Tabs: Patient Profile · HPT Dynamics · Autoimmune Markers ·
#         Drug PK · Clinical Endpoints · Scenario Comparison · Mechanistic Map
##############################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)

# ============================================================
# mrgsolve MODEL (inline — same core as ht_mrgsolve_model.R)
# ============================================================
ht_model_code <- '
$PARAM
k_TRH_prod=0.5 k_TRH_deg=1.0 k_TSH_stim=4.0 k_TSH_deg=6.0
EC50_T3_TSH=0.005 hill_TSH=2.5
k_T4_syn=0.12 k_T3_syn=0.008 k_T4sec=1.5 k_T3sec=1.5
EC50_TSH_T4=2.0 hill_T4=1.5
k_T4_deg=0.099 k_T3_deg=0.693 k_D1=0.3 k_D2=0.15
k_T4_tissue=0.2 k_T4_tissue_ret=0.1 k_T3_tissue=0.35 k_T3_tissue_ret=0.15
k_Th1_base=0.01 k_Th1_stim=2.5 k_Th1_treg=3.0 k_Th1_deg=0.5
k_Treg_prod=0.3 k_Treg_deg=0.4 k_Treg_inh=0.8 Treg_base=0.75
k_Bcell_stim=1.5 k_Bcell_deg=0.6
k_Ab_prod=0.3 k_Ab_prod2=0.2 k_Ab_deg=0.05
Th1_0=0.15 Treg_0=0.65 Bcell_0=0.20 Ab_TPO_0=150.0 Ab_Tg_0=120.0
k_dmg=0.8 k_repair=0.2 DmgThy_0=0.25 k_dmg_scale=0.001
MW_T4=776.87 F_LT4=0.75 ka_LT4=0.48 CL_LT4=1.3 Vc_LT4=10.0
Vp_LT4=22.0 k12_LT4=0.25 k21_LT4=0.11 kel_LT4=0.13 f_T4conv=0.35
ka_LiT3=5.0 CL_LiT3=25.0 Vc_LiT3=40.0 kel_LiT3=0.625
F_Se=0.85 k_Se_abs=2.0 k_Se_elim=0.15 Vd_Se=50.0 Se_0=80.0
Emax_Se_GPx=0.6 EC50_Se_GPx=90.0 Emax_Se_Ab=0.55 EC50_Se_Ab=95.0
Se_dose=0.0 LT4_dose=0.0 LiT3_dose=0.0 MMI_block=0.0

$INIT
TRH=0.5 TSH=2.0 T4thy=0.12 T3thy=0.008
T4p=0.10 T3p=0.0045 T4t=0.15 T3t=0.008
Th1=0.15 Treg=0.65 Bcell=0.20
AntiTPOAb=150.0 AntiTgAb=120.0
DmgThy=0.25
LT4g=0.0 LT4c=0.0 LT4p=0.0 LiT3c=0.0
Se=80.0

$ODE
double T3_fb = T3p + 0.5*T4p*k_D2;
double TRH_inhib = 1.0/(1.0+pow(T3_fb/EC50_T3_TSH, hill_TSH));
dxdt_TRH = k_TRH_prod*TRH_inhib - k_TRH_deg*TRH;

double T3_pit = T3p + k_D2*T4p*2.0;
double TSH_inhib = pow(EC50_T3_TSH,hill_TSH)/(pow(EC50_T3_TSH,hill_TSH)+pow(T3_pit*0.5,hill_TSH));
dxdt_TSH = k_TSH_stim*TRH*TSH_inhib - k_TSH_deg*TSH;

double TSH_stim = pow(TSH,hill_T4)/(pow(EC50_TSH_T4,hill_T4)+pow(TSH,hill_T4));
double thyroid_func = 1.0 - DmgThy;
double TPO_activity = 1.0 - MMI_block;
dxdt_T4thy = k_T4_syn*TSH_stim*thyroid_func*TPO_activity - k_T4sec*T4thy;
dxdt_T3thy = k_T3_syn*TSH_stim*thyroid_func*TPO_activity - k_T3sec*T3thy;

double Se_eff_GPx = Emax_Se_GPx*Se/(EC50_Se_GPx+Se);
double DIO_Se_eff = 1.0 + Se_eff_GPx*0.3;
double LT4c_nmol  = LT4c/(MW_T4*1e-3*Vc_LT4);
double LiT3c_nmol = LiT3c/(651.0*1e-3*Vc_LiT3);

dxdt_T4p = k_T4sec*T4thy + LT4c_nmol*kel_LT4
           - k_D1*DIO_Se_eff*T4p - k_D2*DIO_Se_eff*T4p
           - k_T4_deg*T4p - k_T4_tissue*T4p + k_T4_tissue_ret*T4t;
dxdt_T3p = k_T3sec*T3thy + (k_D1+k_D2)*DIO_Se_eff*T4p
           + LiT3c_nmol*kel_LiT3 - k_T3_deg*T3p
           - k_T3_tissue*T3p + k_T3_tissue_ret*T3t;
dxdt_T4t = k_T4_tissue*T4p - k_T4_tissue_ret*T4t - k_D2*T4t*0.5;
dxdt_T3t = k_T3_tissue*T3p + k_D2*T4t*0.5 - k_T3_tissue_ret*T3t - k_T3_deg*T3t*0.5;

double Ag_drive  = 0.1 + DmgThy*0.9;
dxdt_Th1  = k_Th1_base + k_Th1_stim*Ag_drive*(1.0-Th1) - k_Th1_treg*Treg*Th1 - k_Th1_deg*Th1;
dxdt_Treg = k_Treg_prod*Treg_base - k_Treg_inh*Th1*Treg - k_Treg_deg*Treg;
dxdt_Bcell= k_Bcell_stim*Th1*(1.0-Bcell) - k_Bcell_deg*Bcell;

double Se_eff_Ab = 1.0 - Emax_Se_Ab*Se/(EC50_Se_Ab+Se);
dxdt_AntiTPOAb = k_Ab_prod*Bcell*1000.0*Se_eff_Ab - k_Ab_deg*AntiTPOAb;
dxdt_AntiTgAb  = k_Ab_prod2*Bcell*800.0*Se_eff_Ab  - k_Ab_deg*AntiTgAb;

double dmg_drive = k_dmg*Th1*(AntiTPOAb/500.0)*(1.0+Se_eff_GPx*(-0.4));
dxdt_DmgThy = dmg_drive*k_dmg_scale*(1.0-DmgThy) - k_repair*Treg*DmgThy;

dxdt_LT4g = LT4_dose*F_LT4 - ka_LT4*LT4g;
dxdt_LT4c = ka_LT4*LT4g - (kel_LT4+k12_LT4)*LT4c + k21_LT4*LT4p;
dxdt_LT4p = k12_LT4*LT4c - k21_LT4*LT4p;
dxdt_LiT3c= LiT3_dose*0.95 - (kel_LiT3+k12_LT4*0.5)*LiT3c;
dxdt_Se   = Se_dose*F_Se/Vd_Se - k_Se_elim*Se;

$TABLE
double fT4_pmolL  = T4p*1000.0;
double fT3_pmolL  = T3p*1000.0;
double TSH_clin   = TSH;
double T3_norm    = T3t/0.008;
double HypoScore  = 10.0*(1.0-fmin(T3_norm,1.0));
double HR_effect  = 60.0 + 20.0*fmin(T3_norm,1.5);
double LDL_rel    = 4.0/fmax(T3_norm,0.3);
double Fatigue    = 5.0*(1.0-T3_norm)+2.0;
double BMD_risk   = -0.5*(LT4c/100.0);
double ThyVol_rel = 1.0+0.3*(0.5-DmgThy)*2.0+0.1*(AntiTPOAb/300.0);

$CAPTURE fT4_pmolL fT3_pmolL TSH_clin AntiTPOAb AntiTgAb
DmgThy Th1 Treg Bcell HypoScore HR_effect LDL_rel Fatigue BMD_risk
ThyVol_rel T4p T3p LT4c LiT3c Se
'

ht_mod <- suppressMessages(mrgsolve::mcode("ht_shiny", ht_model_code))

# ============================================================
# HELPER FUNCTIONS
# ============================================================

run_sim <- function(LT4=0, Se=0, LiT3=0, MMI=0,
                    years=3, Dmg0=0.25, Ab0=150, Treg0=0.65) {
  ht_mod %>%
    param(LT4_dose=LT4, Se_dose=Se, LiT3_dose=LiT3, MMI_block=MMI) %>%
    init(DmgThy=Dmg0, AntiTPOAb=Ab0, Treg=Treg0) %>%
    mrgsim(end=years*365, delta=7, obsonly=TRUE) %>%
    as_tibble() %>%
    mutate(Time_years = time/365, Time_months = time/30)
}

theme_qsp <- function() {
  theme_bw(base_size=12) +
    theme(plot.title=element_text(face="bold", size=12),
          legend.position="bottom",
          panel.grid.minor=element_blank())
}

# ============================================================
# UI
# ============================================================

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = "Hashimoto's Thyroiditis QSP",
    titleWidth = 320
  ),

  dashboardSidebar(
    width = 320,
    sidebarMenu(
      id = "tabs",
      menuItem("Patient Profile",        tabName="tab_patient",   icon=icon("user-md")),
      menuItem("HPT Axis Dynamics",      tabName="tab_hpt",       icon=icon("chart-line")),
      menuItem("Autoimmune Markers",     tabName="tab_immune",    icon=icon("shield-virus")),
      menuItem("Drug PK",                tabName="tab_pk",        icon=icon("capsules")),
      menuItem("Clinical Endpoints",     tabName="tab_clinical",  icon=icon("heartbeat")),
      menuItem("Scenario Comparison",    tabName="tab_scenario",  icon=icon("layer-group")),
      menuItem("Mechanistic Map",        tabName="tab_map",       icon=icon("project-diagram"))
    ),

    hr(),
    h4("Treatment Settings", style="padding-left:15px;color:white;"),

    sliderInput("lt4_dose", "Levothyroxine (μg/day)",
                min=0, max=250, value=0, step=12.5),
    sliderInput("se_dose", "Selenium (μg/day)",
                min=0, max=400, value=0, step=25),
    sliderInput("lit3_dose", "Liothyronine T3 (μg/day)",
                min=0, max=20, value=0, step=2.5),
    sliderInput("mmi_block", "Methimazole Effect (0-1)",
                min=0, max=1, value=0, step=0.05),

    hr(),
    h4("Patient Parameters", style="padding-left:15px;color:white;"),
    sliderInput("dmg0", "Initial Thyroid Damage (%)",
                min=0, max=80, value=25, step=5),
    sliderInput("ab0", "Initial Anti-TPO Ab (IU/mL)",
                min=34, max=1000, value=150, step=25),
    sliderInput("treg0", "Treg Level (relative, 0-1)",
                min=0.1, max=1.0, value=0.65, step=0.05),
    sliderInput("years", "Simulation Duration (years)",
                min=1, max=10, value=5, step=1),

    hr(),
    actionButton("run_btn", "Run Simulation",
                 class="btn-success btn-lg", width="90%",
                 style="margin:5px 15px;")
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side { background-color: #F5F5F5; }
        .box { border-top: 3px solid #2471A3; }
        .value-box .icon-large { font-size: 55px; }
      "))
    ),

    tabItems(
      # ------------------------------------------------------------------
      # TAB 1: Patient Profile
      # ------------------------------------------------------------------
      tabItem(tabName="tab_patient",
        fluidRow(
          box(title="Disease Overview: Hashimoto's Thyroiditis", width=12,
              status="primary", solidHeader=TRUE,
              p(strong("Autoimmune chronic thyroiditis"), " — the most common autoimmune
              disease (prevalence ~2–5% globally, female:male ≈ 7:1). Characterized by
              lymphocytic infiltration of the thyroid gland, anti-TPO and anti-Tg antibodies,
              and progressive thyroid failure leading to hypothyroidism."),
              p("This QSP model integrates the", strong("HPT axis"), ",
              thyroid hormone synthesis (T4/T3), peripheral deiodination (DIO1/DIO2/DIO3),
              autoimmune T-cell and B-cell dynamics, thyroid damage, and pharmacology of",
              strong("levothyroxine (LT4), liothyronine (T3), and selenium."))
          )
        ),
        fluidRow(
          valueBoxOutput("vb_tsh", width=3),
          valueBoxOutput("vb_ft4", width=3),
          valueBoxOutput("vb_antitpo", width=3),
          valueBoxOutput("vb_dmg", width=3)
        ),
        fluidRow(
          box(title="Current Patient Parameters", width=6, status="info",
              tableOutput("patient_params_table")),
          box(title="Disease Status at 1 Year", width=6, status="warning",
              tableOutput("oneyear_table"))
        ),
        fluidRow(
          box(title="Model Structure Summary", width=12, status="primary",
              tags$ul(
                tags$li(strong("18 ODEs:"), "TRH, TSH, thyroidal T4/T3, plasma/tissue T4/T3, Th1, Treg, Bcell, Anti-TPO Ab, Anti-Tg Ab, thyroid damage, LT4 2-cpt PK, LiT3 1-cpt PK, Selenium"),
                tags$li(strong("HPT Axis:"), "TRH → TSH → thyroid T4/T3 synthesis with sigmoidal feedback (T3-mediated TSH suppression via pituitary DIO2)"),
                tags$li(strong("Autoimmune:"), "Th1/Treg balance drives antigen-dependent thyrocyte destruction; B-cell/plasma cell antibody production; selenium modulates GPx activity"),
                tags$li(strong("Pharmacology:"), "LT4 2-compartment oral model (F=75%, t½=7 days); selenium dose-response on Anti-TPO reduction (~50% at 200 μg/day per Gärtner 2002)"),
                tags$li(strong("Key outputs:"), "TSH, fT4, fT3, Anti-TPO Ab, thyroid damage, hypothyroid score, heart rate, LDL-C, BMD risk")
              )
          )
        )
      ),

      # ------------------------------------------------------------------
      # TAB 2: HPT Axis Dynamics
      # ------------------------------------------------------------------
      tabItem(tabName="tab_hpt",
        fluidRow(
          box(title="TSH Dynamics (mIU/L)", width=6, status="primary",
              solidHeader=TRUE,
              plotlyOutput("plt_tsh", height="300px")),
          box(title="Free T4 (pmol/L)", width=6, status="primary",
              solidHeader=TRUE,
              plotlyOutput("plt_ft4", height="300px"))
        ),
        fluidRow(
          box(title="Free T3 (pmol/L)", width=6, status="success",
              solidHeader=TRUE,
              plotlyOutput("plt_ft3", height="300px")),
          box(title="TRH Dynamics (normalized)", width=6, status="success",
              solidHeader=TRUE,
              plotlyOutput("plt_trh", height="300px"))
        ),
        fluidRow(
          box(title="Thyroid Gland — T4/T3 Synthesis", width=12, status="info",
              solidHeader=TRUE,
              plotlyOutput("plt_thyroid_synth", height="280px"))
        )
      ),

      # ------------------------------------------------------------------
      # TAB 3: Autoimmune Markers
      # ------------------------------------------------------------------
      tabItem(tabName="tab_immune",
        fluidRow(
          box(title="Anti-TPO Antibody (IU/mL)", width=6, status="danger",
              solidHeader=TRUE,
              plotlyOutput("plt_antitpo", height="300px")),
          box(title="Anti-Tg Antibody (IU/mL)", width=6, status="danger",
              solidHeader=TRUE,
              plotlyOutput("plt_antitg", height="300px"))
        ),
        fluidRow(
          box(title="Th1 vs Treg Dynamics", width=6, status="warning",
              solidHeader=TRUE,
              plotlyOutput("plt_th1treg", height="300px")),
          box(title="B Cell Activation & Th1/Treg Ratio", width=6, status="warning",
              solidHeader=TRUE,
              plotlyOutput("plt_ratio", height="300px"))
        ),
        fluidRow(
          box(title="Thyroid Damage Progression (0–100%)", width=12, status="danger",
              solidHeader=TRUE,
              plotlyOutput("plt_damage", height="280px"))
        )
      ),

      # ------------------------------------------------------------------
      # TAB 4: Drug PK
      # ------------------------------------------------------------------
      tabItem(tabName="tab_pk",
        fluidRow(
          box(title="Levothyroxine Plasma Concentration (LT4 Central, μg)", width=6,
              status="primary", solidHeader=TRUE,
              plotlyOutput("plt_lt4pk", height="300px")),
          box(title="LT4 Gut + Peripheral Compartments", width=6,
              status="primary", solidHeader=TRUE,
              plotlyOutput("plt_lt4_2cpt", height="300px"))
        ),
        fluidRow(
          box(title="Selenium Plasma (μg/L)", width=6, status="success",
              solidHeader=TRUE,
              plotlyOutput("plt_se", height="300px")),
          box(title="Liothyronine (T3) Plasma", width=6, status="success",
              solidHeader=TRUE,
              plotlyOutput("plt_lit3", height="300px"))
        ),
        fluidRow(
          box(title="LT4 Dose-Response Summary", width=12, status="info",
              solidHeader=TRUE,
              DTOutput("lt4_dr_table"))
        )
      ),

      # ------------------------------------------------------------------
      # TAB 5: Clinical Endpoints
      # ------------------------------------------------------------------
      tabItem(tabName="tab_clinical",
        fluidRow(
          box(title="Hypothyroid Symptom Score (0–10)", width=6, status="danger",
              solidHeader=TRUE,
              plotlyOutput("plt_hyposcore", height="300px")),
          box(title="Heart Rate Effect (bpm)", width=6, status="info",
              solidHeader=TRUE,
              plotlyOutput("plt_hr", height="300px"))
        ),
        fluidRow(
          box(title="LDL-Cholesterol (relative, mmol/L)", width=6, status="warning",
              solidHeader=TRUE,
              plotlyOutput("plt_ldl", height="300px")),
          box(title="Thyroid Volume (relative)", width=6, status="warning",
              solidHeader=TRUE,
              plotlyOutput("plt_thyrvol", height="300px"))
        ),
        fluidRow(
          box(title="Fatigue Score & BMD Risk", width=12, status="danger",
              solidHeader=TRUE,
              plotlyOutput("plt_fatigue_bmd", height="280px"))
        )
      ),

      # ------------------------------------------------------------------
      # TAB 6: Scenario Comparison
      # ------------------------------------------------------------------
      tabItem(tabName="tab_scenario",
        fluidRow(
          box(title="Scenario Comparison Parameters", width=12,
              status="primary", solidHeader=TRUE,
              fluidRow(
                column(3, checkboxGroupInput("scen_sel",
                  "Select Scenarios:",
                  choices=c("Untreated"="sc1",
                            "LT4 100 μg"="sc2",
                            "Se 200 μg"="sc3",
                            "LT4+Se"="sc4",
                            "LT4+T3"="sc5",
                            "High LT4"="sc6",
                            "Early LT4+Se"="sc7"),
                  selected=c("sc1","sc2","sc3","sc4"))),
                column(9,
                  plotlyOutput("plt_scen_tsh", height="250px"))
              )
          )
        ),
        fluidRow(
          box(title="Anti-TPO Ab — Scenario Comparison", width=6,
              status="danger", solidHeader=TRUE,
              plotlyOutput("plt_scen_ab", height="300px")),
          box(title="Thyroid Damage — Scenario Comparison", width=6,
              status="danger", solidHeader=TRUE,
              plotlyOutput("plt_scen_dmg", height="300px"))
        ),
        fluidRow(
          box(title="Scenario Summary Table (1-year & 5-year outcomes)", width=12,
              status="info", solidHeader=TRUE,
              DTOutput("scen_table"))
        )
      ),

      # ------------------------------------------------------------------
      # TAB 7: Mechanistic Map
      # ------------------------------------------------------------------
      tabItem(tabName="tab_map",
        fluidRow(
          box(title="QSP Mechanistic Map — Hashimoto's Thyroiditis", width=12,
              status="primary", solidHeader=TRUE,
              p("Full mechanistic map: HPT Axis · Thyroid Follicular Biology ·
                Peripheral T4/T3 Metabolism · Autoimmune Pathogenesis · Drug PK/PD ·
                Clinical Endpoints & Comorbidities"),
              tags$a("View full-size SVG →",
                     href="ht_qsp_model.svg", target="_blank"),
              br(), br(),
              imageOutput("mech_map_img", height="700px")
          )
        )
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  # Reactive simulation
  sim_data <- eventReactive(input$run_btn, {
    withProgress(message="Running ODE simulation...", value=0.5, {
      run_sim(LT4    = input$lt4_dose,
              Se     = input$se_dose,
              LiT3   = input$lit3_dose,
              MMI    = input$mmi_block,
              years  = input$years,
              Dmg0   = input$dmg0/100,
              Ab0    = input$ab0,
              Treg0  = input$treg0)
    })
  }, ignoreNULL=FALSE)

  # Default run at startup
  observe({
    if (is.null(sim_data())) {
      isolate({ run_sim() })
    }
  })

  # ---- VALUE BOXES ---------------------------------------------------
  output$vb_tsh <- renderValueBox({
    d <- sim_data()
    val <- round(tail(d$TSH_clin, 1), 2)
    status <- if (val < 0.4) "warning" else if (val > 4.0) "danger" else "success"
    icon_nm <- if (val > 4.0) "arrow-up" else if (val < 0.4) "arrow-down" else "check"
    valueBox(paste(val, "mIU/L"), "TSH (final)", icon=icon(icon_nm), color=status)
  })

  output$vb_ft4 <- renderValueBox({
    d <- sim_data()
    val <- round(tail(d$fT4_pmolL, 1), 1)
    status <- if (val >= 9 && val <= 23) "success" else "danger"
    valueBox(paste(val, "pmol/L"), "Free T4 (final)", icon=icon("flask"), color=status)
  })

  output$vb_antitpo <- renderValueBox({
    d <- sim_data()
    val <- round(tail(d$AntiTPOAb, 1), 0)
    status <- if (val > 200) "danger" else if (val > 34) "warning" else "success"
    valueBox(paste(val, "IU/mL"), "Anti-TPO Ab (final)", icon=icon("shield-virus"), color=status)
  })

  output$vb_dmg <- renderValueBox({
    d <- sim_data()
    val <- round(tail(d$DmgThy, 1)*100, 1)
    status <- if (val > 60) "danger" else if (val > 30) "warning" else "success"
    valueBox(paste(val, "%"), "Thyroid Damage (final)", icon=icon("heartbeat"), color=status)
  })

  # ---- PATIENT TABLE ---------------------------------------------------
  output$patient_params_table <- renderTable({
    data.frame(
      Parameter = c("LT4 Dose", "Selenium Dose", "Liothyronine", "MMI Effect",
                    "Initial Damage", "Baseline Anti-TPO Ab", "Treg Level"),
      Value     = c(paste(input$lt4_dose, "μg/day"),
                    paste(input$se_dose, "μg/day"),
                    paste(input$lit3_dose, "μg/day"),
                    input$mmi_block,
                    paste(input$dmg0, "%"),
                    paste(input$ab0, "IU/mL"),
                    input$treg0),
      stringsAsFactors=FALSE
    )
  }, striped=TRUE, bordered=TRUE)

  output$oneyear_table <- renderTable({
    d <- sim_data()
    d1 <- d %>% filter(abs(Time_years-1) == min(abs(Time_years-1))) %>% slice(1)
    data.frame(
      Endpoint = c("TSH (mIU/L)", "fT4 (pmol/L)", "fT3 (pmol/L)",
                   "Anti-TPO Ab (IU/mL)", "Anti-Tg Ab (IU/mL)",
                   "Thyroid Damage (%)", "Hypothyroid Score",
                   "Heart Rate (bpm)", "LDL Relative (mmol/L)"),
      Value_1yr = c(round(d1$TSH_clin, 2), round(d1$fT4_pmolL, 1), round(d1$fT3_pmolL, 2),
                    round(d1$AntiTPOAb, 0), round(d1$AntiTgAb, 0),
                    round(d1$DmgThy*100, 1), round(d1$HypoScore, 1),
                    round(d1$HR_effect, 0), round(d1$LDL_rel, 2)),
      Normal_Range = c("0.4–4.0", "9–23", "3.5–6.5", "<34", "<40",
                       "0", "0", "60–80", "<3.0", stringsAsFactors=FALSE)
    )
  }, striped=TRUE, bordered=TRUE)

  # ---- HPT PLOTS ---------------------------------------------------
  make_plotly <- function(d, x, y, title, ylab, ref_lines=NULL) {
    p <- ggplot(d, aes_string(x=x, y=y)) +
      geom_line(color="#2471A3", linewidth=1.1)
    if (!is.null(ref_lines)) {
      for (rl in ref_lines) {
        p <- p + geom_hline(yintercept=rl$y, linetype="dashed",
                            color=rl$color, linewidth=0.7)
      }
    }
    p <- p + labs(title=title, x="Time (years)", y=ylab) + theme_qsp()
    ggplotly(p) %>% layout(margin=list(t=30))
  }

  output$plt_tsh <- renderPlotly({
    d <- sim_data()
    make_plotly(d, "Time_years", "TSH_clin", "TSH (mIU/L)", "mIU/L",
                list(list(y=0.4, color="green"), list(y=4.0, color="green")))
  })

  output$plt_ft4 <- renderPlotly({
    d <- sim_data()
    make_plotly(d, "Time_years", "fT4_pmolL", "Free T4 (pmol/L)", "pmol/L",
                list(list(y=9, color="green"), list(y=23, color="green")))
  })

  output$plt_ft3 <- renderPlotly({
    d <- sim_data()
    make_plotly(d, "Time_years", "fT3_pmolL", "Free T3 (pmol/L)", "pmol/L",
                list(list(y=3.5, color="green"), list(y=6.5, color="green")))
  })

  output$plt_trh <- renderPlotly({
    d <- sim_data()
    make_plotly(d, "Time_years", "TRH", "TRH Dynamics (normalized)", "TRH (rel)")
  })

  output$plt_thyroid_synth <- renderPlotly({
    d <- sim_data()
    p <- plot_ly(d) %>%
      add_lines(x=~Time_years, y=~T4p*1000, name="fT4 (pmol/L)", line=list(color="#2980B9")) %>%
      add_lines(x=~Time_years, y=~T3p*1000*10, name="fT3 × 10 (pmol/L)", line=list(color="#E74C3C")) %>%
      layout(title="Plasma T4/T3 (fT3 × 10 for scale)",
             xaxis=list(title="Time (years)"), yaxis=list(title="pmol/L"))
    p
  })

  # ---- AUTOIMMUNE PLOTS ---------------------------------------------------
  output$plt_antitpo <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=Time_years, y=AntiTPOAb)) +
      geom_line(color="#C0392B", linewidth=1.1) +
      geom_hline(yintercept=34, linetype="dashed", color="green", linewidth=0.8) +
      annotate("text", x=max(d$Time_years)*0.8, y=50, label="Upper normal: 34 IU/mL",
               color="grey40", size=3.5) +
      labs(x="Time (years)", y="IU/mL") + theme_qsp()
    ggplotly(p)
  })

  output$plt_antitg <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=Time_years, y=AntiTgAb)) +
      geom_line(color="#E74C3C", linewidth=1.1) +
      geom_hline(yintercept=40, linetype="dashed", color="green", linewidth=0.8) +
      labs(x="Time (years)", y="IU/mL") + theme_qsp()
    ggplotly(p)
  })

  output$plt_th1treg <- renderPlotly({
    d <- sim_data()
    p <- plot_ly(d) %>%
      add_lines(x=~Time_years, y=~Th1, name="Th1", line=list(color="#E74C3C")) %>%
      add_lines(x=~Time_years, y=~Treg, name="Treg", line=list(color="#27AE60")) %>%
      layout(xaxis=list(title="Time (years)"), yaxis=list(title="Relative units"))
    p
  })

  output$plt_ratio <- renderPlotly({
    d <- sim_data()
    p <- plot_ly(d) %>%
      add_lines(x=~Time_years, y=~Th1/Treg, name="Th1/Treg Ratio",
                line=list(color="#9B59B6")) %>%
      add_lines(x=~Time_years, y=~Bcell*5, name="B Cells × 5",
                line=list(color="#3498DB")) %>%
      layout(xaxis=list(title="Time (years)"), yaxis=list(title="Ratio / Relative"))
    p
  })

  output$plt_damage <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=Time_years, y=DmgThy*100)) +
      geom_area(fill="#E74C3C", alpha=0.3) +
      geom_line(color="#C0392B", linewidth=1.2) +
      geom_hline(yintercept=50, linetype="dashed", color="red") +
      annotate("text", x=max(d$Time_years)*0.8, y=55, label="50% damage",
               color="red", size=3.5) +
      labs(x="Time (years)", y="Damage (%)") + theme_qsp()
    ggplotly(p)
  })

  # ---- DRUG PK PLOTS ---------------------------------------------------
  output$plt_lt4pk <- renderPlotly({
    d <- sim_data()
    make_plotly(d, "Time_years", "LT4c", "LT4 Central Concentration", "LT4c (μg)")
  })

  output$plt_lt4_2cpt <- renderPlotly({
    d <- sim_data()
    p <- plot_ly(d) %>%
      add_lines(x=~Time_years, y=~LT4g, name="Gut", line=list(color="#E67E22")) %>%
      add_lines(x=~Time_years, y=~LT4c, name="Central", line=list(color="#2980B9")) %>%
      layout(xaxis=list(title="Time (years)"), yaxis=list(title="Amount (μg)"))
    p
  })

  output$plt_se <- renderPlotly({
    d <- sim_data()
    p <- ggplot(d, aes(x=Time_years, y=Se)) +
      geom_line(color="#27AE60", linewidth=1.1) +
      geom_hline(yintercept=90, linetype="dashed", color="blue") +
      annotate("text", x=max(d$Time_years)*0.7, y=95, label="Normal Se ~90 μg/L",
               color="blue", size=3.5) +
      labs(x="Time (years)", y="Selenium (μg/L)") + theme_qsp()
    ggplotly(p)
  })

  output$plt_lit3 <- renderPlotly({
    d <- sim_data()
    make_plotly(d, "Time_years", "LiT3c", "Liothyronine Plasma (μg)", "LiT3c (μg)")
  })

  output$lt4_dr_table <- renderDT({
    doses <- c(0, 25, 50, 75, 100, 125, 150, 175)
    res <- lapply(doses, function(ld) {
      withProgress(message=paste("LT4 =", ld), value=ld/175, {
        ht_mod %>%
          param(LT4_dose=ld) %>%
          init(DmgThy=input$dmg0/100, AntiTPOAb=input$ab0) %>%
          mrgsim(end=365, delta=365, obsonly=TRUE) %>%
          as_tibble() %>% slice_tail(n=1) %>%
          transmute(`LT4 (μg/day)`=ld,
                    `TSH (mIU/L)`=round(TSH_clin,2),
                    `fT4 (pmol/L)`=round(fT4_pmolL,1),
                    `fT3 (pmol/L)`=round(fT3_pmolL,2),
                    `Hypo Score`=round(HypoScore,1),
                    `HR (bpm)`=round(HR_effect,0),
                    `LDL (mmol/L)`=round(LDL_rel,2))
      })
    }) %>% bind_rows()
    datatable(res, options=list(pageLength=10, dom="t"),
              rownames=FALSE) %>%
      formatStyle("TSH (mIU/L)",
        backgroundColor = styleInterval(c(0.4, 4.0),
          c("#F5B7B1","#ABEBC6","#F9E79F")))
  })

  # ---- CLINICAL ENDPOINTS ---------------------------------------------------
  output$plt_hyposcore <- renderPlotly({
    d <- sim_data()
    make_plotly(d, "Time_years", "HypoScore", "Hypothyroid Symptom Score", "Score (0-10)",
                list(list(y=2, color="green"), list(y=6, color="red")))
  })

  output$plt_hr <- renderPlotly({
    d <- sim_data()
    make_plotly(d, "Time_years", "HR_effect", "Heart Rate (bpm)", "bpm",
                list(list(y=60, color="green"), list(y=80, color="green")))
  })

  output$plt_ldl <- renderPlotly({
    d <- sim_data()
    make_plotly(d, "Time_years", "LDL_rel", "LDL-C (relative, mmol/L)", "mmol/L",
                list(list(y=3.0, color="orange")))
  })

  output$plt_thyrvol <- renderPlotly({
    d <- sim_data()
    make_plotly(d, "Time_years", "ThyVol_rel", "Thyroid Volume (relative)", "Rel. volume")
  })

  output$plt_fatigue_bmd <- renderPlotly({
    d <- sim_data()
    p <- plot_ly(d) %>%
      add_lines(x=~Time_years, y=~Fatigue, name="Fatigue Score",
                line=list(color="#E74C3C")) %>%
      add_lines(x=~Time_years, y=~BMD_risk*(-20), name="BMD Risk × (-20)",
                line=list(color="#8E44AD", dash="dot")) %>%
      layout(xaxis=list(title="Time (years)"),
             yaxis=list(title="Score / BMD risk (scaled)"))
    p
  })

  # ---- SCENARIO COMPARISON ---------------------------------------------------
  scen_colors <- c(sc1="#E74C3C", sc2="#3498DB", sc3="#2ECC71",
                   sc4="#9B59B6", sc5="#E67E22", sc6="#1ABC9C", sc7="#F39C12")
  scen_labels <- c(sc1="Untreated", sc2="LT4 100 μg", sc3="Se 200 μg",
                   sc4="LT4+Se", sc5="LT4+T3", sc6="High LT4 175 μg",
                   sc7="Early LT4+Se")
  scen_params <- list(
    sc1=c(0,0,0,0), sc2=c(100,0,0,0), sc3=c(0,200,0,0),
    sc4=c(100,200,0,0), sc5=c(100,0,7.5,0),
    sc6=c(175,0,0,0), sc7=c(75,200,0,0)
  )

  all_scen_data <- reactive({
    selected <- input$scen_sel
    if (length(selected) == 0) return(NULL)
    withProgress(message="Running scenario comparison...", value=0, {
      lapply(selected, function(sc) {
        incProgress(1/length(selected))
        p <- scen_params[[sc]]
        run_sim(LT4=p[1], Se=p[2], LiT3=p[3], MMI=p[4],
                years=input$years,
                Dmg0=input$dmg0/100, Ab0=input$ab0, Treg0=input$treg0) %>%
          mutate(Scenario=scen_labels[sc], sc_id=sc)
      }) %>% bind_rows()
    })
  })

  output$plt_scen_tsh <- renderPlotly({
    d <- all_scen_data()
    if (is.null(d)) return(plotly_empty())
    p <- plot_ly()
    for (sc in unique(d$sc_id)) {
      dd <- d %>% filter(sc_id == sc)
      p <- p %>% add_lines(data=dd, x=~Time_years, y=~TSH_clin,
                            name=unique(dd$Scenario),
                            line=list(color=scen_colors[sc]))
    }
    p %>% layout(xaxis=list(title="Time (years)"), yaxis=list(title="TSH (mIU/L)"),
                 shapes=list(
                   list(type="rect", x0=0, x1=input$years, y0=0.4, y1=4.0,
                        fillcolor="rgba(0,200,0,0.05)", line=list(width=0))
                 ))
  })

  output$plt_scen_ab <- renderPlotly({
    d <- all_scen_data()
    if (is.null(d)) return(plotly_empty())
    p <- plot_ly()
    for (sc in unique(d$sc_id)) {
      dd <- d %>% filter(sc_id == sc)
      p <- p %>% add_lines(data=dd, x=~Time_years, y=~AntiTPOAb,
                            name=unique(dd$Scenario),
                            line=list(color=scen_colors[sc]))
    }
    p %>% layout(xaxis=list(title="Time (years)"), yaxis=list(title="IU/mL"))
  })

  output$plt_scen_dmg <- renderPlotly({
    d <- all_scen_data()
    if (is.null(d)) return(plotly_empty())
    p <- plot_ly()
    for (sc in unique(d$sc_id)) {
      dd <- d %>% filter(sc_id == sc)
      p <- p %>% add_lines(data=dd, x=~Time_years, y=~DmgThy*100,
                            name=unique(dd$Scenario),
                            line=list(color=scen_colors[sc]))
    }
    p %>% layout(xaxis=list(title="Time (years)"), yaxis=list(title="Damage (%)"))
  })

  output$scen_table <- renderDT({
    d <- all_scen_data()
    if (is.null(d)) return(datatable(data.frame()))
    tbl <- d %>%
      group_by(Scenario) %>%
      summarise(
        `TSH at 1yr`      = round(mean(TSH_clin[abs(Time_years-1)<0.1]),2),
        `fT4 at 1yr`      = round(mean(fT4_pmolL[abs(Time_years-1)<0.1]),1),
        `Anti-TPO 1yr`    = round(mean(AntiTPOAb[abs(Time_years-1)<0.1]),0),
        `Damage 5yr (%)`  = round(mean(DmgThy[Time_years>=(max(Time_years)-0.1)])*100,1),
        `HypoScore 1yr`   = round(mean(HypoScore[abs(Time_years-1)<0.1]),1),
        .groups="drop"
      )
    datatable(tbl, rownames=FALSE, options=list(dom="t")) %>%
      formatStyle("TSH at 1yr",
                  backgroundColor=styleInterval(c(0.4,4.0),c("#F5B7B1","#ABEBC6","#F9E79F")))
  })

  # ---- MECHANISTIC MAP ---------------------------------------------------
  output$mech_map_img <- renderImage({
    list(src="ht_qsp_model.png", contentType="image/png",
         width="100%", alt="Hashimoto's Thyroiditis QSP Map")
  }, deleteFile=FALSE)
}

# ============================================================
# RUN APP
# ============================================================
shinyApp(ui=ui, server=server)
