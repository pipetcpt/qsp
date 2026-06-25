# ============================================================
# Ischemic Stroke QSP — Interactive Shiny Dashboard
# 허혈성 뇌졸중 QSP 대시보드
#
# Tabs:
#   1. Patient Profile & Risk Factors
#   2. Acute Treatment PK (tPA)
#   3. Ischemic Cascade (CBF / ATP / Excitotoxicity)
#   4. Neuroinflammation & BBB
#   5. Clinical Endpoints (NIHSS / Infarct / mRS)
#   6. Secondary Prevention PK (Aspirin / Apixaban)
#   7. Scenario Comparison
#   8. Biomarker Dynamics
# ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(DT)
library(plotly)

# ---- mrgsolve model (embedded) --------------------------------
is_code <- '
$PARAM
k_tpa_fibrinol=8.0, k_spont_lysis=0.003,
CBF_normal=55.0, k_cbf_restore=0.20,
CL_tpa=550.0, V1_tpa=3500, Q_tpa=650.0, V2_tpa=4200,
ka_asp=2.0, CL_asp=10.0, V_asp=12.0, F_asp=0.70,
IC50_asp_cox=0.025,
ka_noac=3.3, CL_noac=3.3, V1_noac=21.0, Q_noac=3.7,
V2_noac=25.0, F_noac=0.50, IC50_noac_xa=0.10,
k_atp_deplete=0.50, k_atp_recover=0.35,
k_glut_release=2.5, k_glut_clear=0.80,
k_ca2_influx=1.20, k_ca2_efflux=0.60, CA2_norm=0.0001,
k_ros_ca2=0.40, k_ros_clear=0.50, k_ros_reperfu=2.0,
k_il6_prod=0.30, k_il6_clear=0.154,
k_bbb_damage=0.08, k_bbb_repair=0.015,
k_pen_convert=0.15, k_pen_salvage=0.30,
k_nihss_worsen=0.04, k_nihss_improve=0.025, k_neuropl=0.002,
NIHSS_init=14.0

$INIT
THROMBUS=1,CBF_CORE=8,CBF_PEN=15,TPA_CENT=0,TPA_PERI=0,
ASP_GUT=0,ASP_CENT=0,NOAC_GUT=0,NOAC_CENT=0,NOAC_PERI=0,
ATP_PEN=1,GLUT=0.1,CA2=0.0001,ROS=0.1,IL6=2,BBB=1,INFARCT=5,NIHSS=14

$ODE
double tPA_Cp=TPA_CENT/V1_tpa;
double recanal=1.0-THROMBUS;
dxdt_THROMBUS=-k_tpa_fibrinol*tPA_Cp*THROMBUS-k_spont_lysis*THROMBUS;
if(THROMBUS<0) THROMBUS=0;
dxdt_CBF_CORE=k_cbf_restore*recanal*(8.0*2-CBF_CORE)-0.01;
dxdt_CBF_PEN=k_cbf_restore*recanal*(CBF_normal-CBF_PEN);
double k10t=(CL_tpa/V1_tpa)*60;double k12t=(Q_tpa/V1_tpa)*60;double k21t=(Q_tpa/V2_tpa)*60;
dxdt_TPA_CENT=-k10t*TPA_CENT-k12t*TPA_CENT+k21t*TPA_PERI;
dxdt_TPA_PERI=k12t*TPA_CENT-k21t*TPA_PERI;
double kael=CL_asp/V_asp;
dxdt_ASP_GUT=-ka_asp*ASP_GUT;
dxdt_ASP_CENT=ka_asp*F_asp*ASP_GUT-kael*ASP_CENT;
double k10n=CL_noac/V1_noac;double k12n=Q_noac/V1_noac;double k21n=Q_noac/V2_noac;
dxdt_NOAC_GUT=-ka_noac*NOAC_GUT;
dxdt_NOAC_CENT=ka_noac*F_noac*NOAC_GUT-k10n*NOAC_CENT-k12n*NOAC_CENT+k21n*NOAC_PERI;
dxdt_NOAC_PERI=k12n*NOAC_CENT-k21n*NOAC_PERI;
double cbf_f=CBF_PEN/CBF_normal;
dxdt_ATP_PEN=k_atp_recover*cbf_f*(1-ATP_PEN)-k_atp_deplete*(1-cbf_f)*ATP_PEN;
dxdt_GLUT=k_glut_release*(1-ATP_PEN)*10-k_glut_clear*ATP_PEN*GLUT;
dxdt_CA2=k_ca2_influx*GLUT*CA2_norm-k_ca2_efflux*ATP_PEN*CA2;
double rrf=(recanal>0.3)?k_ros_reperfu:1.0;
dxdt_ROS=k_ros_ca2*(CA2/CA2_norm)*0.01*rrf-k_ros_clear*ROS;
dxdt_IL6=k_il6_prod*ROS*10-k_il6_clear*IL6;
double mmp9=IL6/(IL6+20.0);
dxdt_BBB=-k_bbb_damage*mmp9*BBB+k_bbb_repair*(1-BBB);
if(BBB<0) BBB=0; if(BBB>1) BBB=1;
dxdt_INFARCT=k_pen_convert*(1-ATP_PEN)*5.0*(1-0.4*(1-BBB));
dxdt_NIHSS=k_nihss_worsen*(INFARCT/60)-k_nihss_improve*recanal*NIHSS-k_neuropl*NIHSS;

$TABLE
double Cp_tpa=TPA_CENT/V1_tpa*1000;
double Cp_asp=ASP_CENT/V_asp;
double Cp_noac=NOAC_CENT/V1_noac*1000;
double recanalization=1.0-THROMBUS;
double mRS_est=(NIHSS<=1)?0:(NIHSS<=4)?1:(NIHSS<=9)?2:(NIHSS<=15)?3:(NIHSS<=20)?4:(NIHSS<=25)?5:6;
double COX1_inh=ASP_CENT/(ASP_CENT+IC50_asp_cox*V_asp);
double Xa_inh=NOAC_CENT/(NOAC_CENT+IC50_noac_xa*V1_noac);
'

mod <- mcode("is_shiny", is_code, quiet = TRUE)

# ---- Helpers --------------------------------------------------
tpa_ev <- function(onset_h, dose_mg) {
  ev(time = onset_h, amt = dose_mg * 0.1,  cmt = "TPA_CENT") +
  ev(time = onset_h, amt = dose_mg * 0.9,  cmt = "TPA_CENT", rate = dose_mg * 0.9)
}
asp_ev_fn <- function(start_h, dose_mg, days) {
  ev(time = seq(start_h, start_h + (days-1)*24, 24), amt = dose_mg, cmt = "ASP_GUT")
}
noac_ev_fn <- function(start_h, dose_mg, days) {
  ts <- sort(as.vector(outer(seq(0, days-1)*24, c(0,12), "+")) + start_h)
  ev(time = ts, amt = dose_mg, cmt = "NOAC_GUT")
}

run_sim <- function(p, ev_list, end_h = 2160, delta = 0.5) {
  tryCatch({
    all_ev <- Reduce(function(a,b) as.ev(a)+as.ev(b), ev_list)
    mrgsim(param(mod, .list = p), events = all_ev, end = end_h, delta = delta) %>%
      as_tibble()
  }, error = function(e) tibble(time=0, error=TRUE))
}

# ============================================================
# UI
# ============================================================
ui <- dashboardPage(
  skin = "red",
  dashboardHeader(title = "Ischemic Stroke QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",       tabName="tab_patient",  icon=icon("user-injured")),
      menuItem("Acute Treatment PK",    tabName="tab_acute",    icon=icon("syringe")),
      menuItem("Ischemic Cascade",      tabName="tab_cascade",  icon=icon("brain")),
      menuItem("Neuroinflamm & BBB",    tabName="tab_inflam",   icon=icon("fire")),
      menuItem("Clinical Endpoints",    tabName="tab_clinical", icon=icon("chart-line")),
      menuItem("Secondary Prevention",  tabName="tab_secondary",icon=icon("pills")),
      menuItem("Scenario Comparison",   tabName="tab_compare",  icon=icon("balance-scale")),
      menuItem("Biomarker Dynamics",    tabName="tab_biom",     icon=icon("vial"))
    ),
    hr(),
    # Global patient parameters
    h5("Patient Parameters", style="margin-left:10px;color:#aaa"),
    sliderInput("bw",     "Body weight (kg)",    min=40, max=120, value=70, step=5),
    sliderInput("age",    "Age (years)",          min=30, max=90,  value=65, step=1),
    selectInput("stroke_type", "Stroke Etiology",
                choices=c("Large artery atherosclerosis","Cardioembolic (AF)","Small vessel","Cryptogenic"),
                selected="Large artery atherosclerosis"),
    sliderInput("nihss0", "Initial NIHSS",        min=1,  max=35,  value=14, step=1),
    sliderInput("infarct0","Initial core (mL)",   min=0,  max=50,  value=5,  step=1),
    sliderInput("pen0",   "Initial penumbra (mL)",min=10, max=100, value=45, step=5)
  ),

  dashboardBody(
    tabItems(

      # ======================== Tab 1: Patient Profile ========================
      tabItem("tab_patient",
        fluidRow(
          box(title="Stroke Parameters", width=6, status="danger",
            sliderInput("onset2treat", "Onset-to-treatment time (hours)", 0.5, 6, 2, 0.5),
            checkboxGroupInput("comorbid", "Comorbidities",
              choices=c("Hypertension","Diabetes","Dyslipidemia","Atrial Fibrillation","Prior Stroke","Smoking"),
              selected=c("Hypertension","Dyslipidemia")),
            numericInput("systolic_bp",  "Systolic BP (mmHg)", 160, 80, 250),
            numericInput("glucose",      "Blood glucose (mg/dL)", 130, 60, 400),
            numericInput("ldl",          "LDL-C (mg/dL)", 120, 30, 300)
          ),
          box(title="Risk Factor Summary", width=6, status="warning",
            DTOutput("risk_table"),
            br(),
            valueBoxOutput("vb_nihss"),
            valueBoxOutput("vb_infarct")
          )
        ),
        fluidRow(
          box(title="Stroke Severity vs. Treatment Window", width=12,
            plotlyOutput("p_severity_window"))
        )
      ),

      # ======================== Tab 2: Acute Treatment PK ========================
      tabItem("tab_acute",
        fluidRow(
          box(title="tPA Dosing", width=4, status="success",
            numericInput("tpa_dose_mg", "tPA dose (mg = 0.9 mg/kg × BW)", value=63, min=10, max=100),
            sliderInput("tpa_onset",   "Onset-to-needle time (hours)", 0.5, 6, 2, 0.5),
            checkboxInput("evt_on",    "Add EVT (endovascular thrombectomy)", FALSE),
            sliderInput("evt_onset",   "EVT onset (hours)",  1, 24, 3, 0.5),
            actionButton("run_acute",  "Run Simulation", class="btn-success btn-block")
          ),
          box(title="PK Parameters", width=4, status="info",
            numericInput("cl_tpa",  "tPA CL (mL/min)", 550, 100, 1500),
            numericInput("v1_tpa",  "tPA V1 (mL)",    3500, 500, 8000),
            numericInput("pai1",    "PAI-1 (relative inhibition, 0–1)", 0.0, 0, 1)
          ),
          box(title="Fibrinolysis Efficacy", width=4, status="warning",
            sliderInput("k_tpa_eff", "tPA fibrinolysis rate", 1, 20, 8, 0.5),
            plotlyOutput("p_recanal_gauge", height="180px")
          )
        ),
        fluidRow(
          box(title="tPA Plasma Concentration — PK Curve", width=6,
            plotlyOutput("p_tpa_pk")),
          box(title="Thrombus Dissolution & Recanalization", width=6,
            plotlyOutput("p_recanal"))
        )
      ),

      # ======================== Tab 3: Ischemic Cascade ========================
      tabItem("tab_cascade",
        fluidRow(
          box(title="Ischemia Parameters", width=3, status="danger",
            sliderInput("cbf_pen",  "Penumbra CBF (mL/100g/min)", 5, 30, 15, 1),
            sliderInput("cbf_core", "Core CBF (mL/100g/min)", 1, 15, 8, 1),
            checkboxInput("collateral_good", "Good collateral circulation", FALSE),
            actionButton("run_cascade", "Run", class="btn-danger btn-block")
          ),
          box(title="CBF Over Time", width=9,
            plotlyOutput("p_cbf"))
        ),
        fluidRow(
          box(title="ATP Depletion (Penumbra)", width=6, plotlyOutput("p_atp")),
          box(title="Glutamate & Ca²⁺",         width=6, plotlyOutput("p_glut_ca2"))
        )
      ),

      # ======================== Tab 4: Neuroinflammation & BBB ========================
      tabItem("tab_inflam",
        fluidRow(
          box(title="Inflammation Parameters", width=3, status="warning",
            sliderInput("k_il6_prod",  "IL-6 production rate", 0.05, 1.0, 0.30, 0.05),
            sliderInput("k_bbb_dam",   "BBB damage rate",       0.01, 0.20, 0.08, 0.01),
            sliderInput("k_bbb_rep",   "BBB repair rate",       0.005,0.05, 0.015,0.005),
            actionButton("run_inflam", "Run", class="btn-warning btn-block")
          ),
          box(title="Neuroinflammation Timeline", width=9,
            plotlyOutput("p_inflam"))
        ),
        fluidRow(
          box(title="BBB Integrity", width=6, plotlyOutput("p_bbb")),
          box(title="Edema Risk (BBB × Infarct Volume)", width=6, plotlyOutput("p_edema"))
        )
      ),

      # ======================== Tab 5: Clinical Endpoints ========================
      tabItem("tab_clinical",
        fluidRow(
          box(title="NIHSS Trajectory (90 days)", width=8, plotlyOutput("p_nihss")),
          box(title="Outcome KPIs", width=4, status="danger",
            valueBoxOutput("vb_90d_nihss",  width=12),
            valueBoxOutput("vb_90d_mrs",    width=12),
            valueBoxOutput("vb_90d_infvol", width=12)
          )
        ),
        fluidRow(
          box(title="Infarct Volume Over Time", width=6, plotlyOutput("p_infarct")),
          box(title="mRS Distribution at 90 Days", width=6, plotlyOutput("p_mrs"))
        )
      ),

      # ======================== Tab 6: Secondary Prevention ========================
      tabItem("tab_secondary",
        fluidRow(
          box(title="Antiplatelet Settings", width=4, status="info",
            selectInput("antiplatelet", "Antiplatelet Agent",
                        c("Aspirin 100 mg QD","Aspirin 300 mg QD","Clopidogrel 75 mg QD",
                          "Dual (Aspirin + Clopidogrel)","None")),
            sliderInput("asp_start", "Start time after stroke (hours)", 0, 72, 24, 6)
          ),
          box(title="Anticoagulation (AF)", width=4, status="info",
            checkboxInput("use_noac", "Use NOAC (apixaban 5 mg BID)", FALSE),
            sliderInput("noac_start","NOAC start (hours post-stroke)", 24, 168, 48, 12),
            selectInput("noac_dose", "Apixaban dose", c("2.5 mg BID","5 mg BID"), "5 mg BID")
          ),
          box(title="Statin / Antihypertensive", width=4, status="info",
            checkboxInput("use_statin","High-intensity statin", TRUE),
            selectInput("statin_type","Statin",
                        c("Atorvastatin 80 mg","Rosuvastatin 40 mg","None")),
            actionButton("run_secondary","Run", class="btn-info btn-block")
          )
        ),
        fluidRow(
          box(title="Aspirin PK (72 hours)", width=6, plotlyOutput("p_asp_pk")),
          box(title="Apixaban PK (7 days)", width=6, plotlyOutput("p_noac_pk"))
        ),
        fluidRow(
          box(title="COX-1 Inhibition (Aspirin)", width=6, plotlyOutput("p_cox1")),
          box(title="Factor Xa Inhibition (Apixaban)", width=6, plotlyOutput("p_xa"))
        )
      ),

      # ======================== Tab 7: Scenario Comparison ========================
      tabItem("tab_compare",
        fluidRow(
          box(title="Select Scenarios", width=3, status="primary",
            checkboxGroupInput("scen_sel", "Scenarios to plot",
              choices = c("Standard tPA 2h","Late tPA 4.5h","Antiplatelet Only",
                          "tPA + NOAC","EVT @ 3h"),
              selected = c("Standard tPA 2h","Late tPA 4.5h","Antiplatelet Only")),
            selectInput("compare_var","Outcome Variable",
              c("NIHSS","INFARCT","mRS_est","BBB","IL6","recanalization","ATP_PEN")),
            actionButton("run_compare","Run All Scenarios", class="btn-primary btn-block")
          ),
          box(title="Scenario Comparison Plot", width=9,
            plotlyOutput("p_compare", height="400px"))
        ),
        fluidRow(
          box(title="90-Day Summary Table", width=12,
            DTOutput("compare_table"))
        )
      ),

      # ======================== Tab 8: Biomarker Dynamics ========================
      tabItem("tab_biom",
        fluidRow(
          box(title="Biomarker Simulation (0–7 days)", width=8,
            plotlyOutput("p_biomarkers", height="400px")),
          box(title="Biomarker Interpretation", width=4,
            h4("Key Biomarkers"),
            tags$ul(
              tags$li(tags$b("GFAP:"), " Astrocyte injury — rises within 1h, peaks ~24h"),
              tags$li(tags$b("UCH-L1:"), " Neuronal damage — early (minutes–hours)"),
              tags$li(tags$b("NSE:"), " Neuronal death — peaks 24–72h"),
              tags$li(tags$b("S100β:"), " Glial/BBB injury — peaks 48–72h"),
              tags$li(tags$b("IL-6:"), " Acute phase — modeled directly in ODE"),
              tags$li(tags$b("D-Dimer:"), " Coagulation — rises with thrombus lysis"),
              tags$li(tags$b("hsCRP:"), " Systemic inflammation — day 1–7")
            ),
            hr(),
            selectInput("biom_var", "Biomarker to highlight",
                        c("IL6","ROS","BBB","GLUT","CA2")),
            plotlyOutput("p_biom_detail", height="200px")
          )
        )
      )
    ) # end tabItems
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # ---- Reactive: update tPA dose based on body weight ----
  observe({
    updateNumericInput(session, "tpa_dose_mg", value = round(input$bw * 0.9))
    updateNumericInput(session, "nihss0", value = input$nihss0)
  })

  # ---- Base simulation (reactive, triggered by main params) ----
  base_sim <- reactive({
    bw <- input$bw
    dose_tpa <- bw * 0.9
    onset <- input$onset2treat
    p_override <- list(
      k_tpa_fibrinol   = isolate(input$k_tpa_eff) %||% 8,
      k_nihss_worsen   = 0.04,
      NIHSS_init       = input$nihss0
    )
    ev_list <- list(tpa_ev(onset, dose_tpa), asp_ev_fn(24, 100, 90))
    run_sim(p_override, ev_list, end_h = 2160)
  }) %>% bindEvent(input$run_acute, input$run_cascade,
                   input$run_inflam, input$run_secondary,
                   ignoreNULL = FALSE)

  # ---- Tab 1: Risk table ----
  output$risk_table <- renderDT({
    tibble(
      Factor    = c("Age","Blood Pressure","LDL-C","Glucose","Comorbidities"),
      Value     = c(input$age, input$systolic_bp, input$ldl, input$glucose,
                    paste(input$comorbid, collapse=", ")),
      `Target`  = c("<75y","<130 mmHg","<70 mg/dL","<140 mg/dL","Controlled"),
      `Status`  = c(
        ifelse(input$age < 75, "OK", "High Risk"),
        ifelse(input$systolic_bp < 130, "OK", "Elevated"),
        ifelse(input$ldl < 70, "At Target", "Above Target"),
        ifelse(input$glucose < 140, "OK", "Elevated"),
        ifelse(length(input$comorbid) <= 1, "Low", "Multiple")
      )
    ) %>% datatable(options = list(dom="t", pageLength=10), rownames=FALSE)
  })

  output$vb_nihss   <- renderValueBox(valueBox(input$nihss0, "Initial NIHSS", icon=icon("chart-bar"), color="red"))
  output$vb_infarct <- renderValueBox(valueBox(paste0(input$infarct0," mL"), "Initial Core", icon=icon("brain"), color="orange"))

  output$p_severity_window <- renderPlotly({
    df <- tibble(
      onset_h = seq(0.5, 6, 0.5),
      nihss_improve = pmax(0, input$nihss0 * (1 - onset_h/8)),
      pen_saved = pmax(0, input$pen0 * (1 - onset_h/7))
    )
    plot_ly(df, x=~onset_h) %>%
      add_lines(y=~nihss_improve, name="Projected NIHSS improvement", line=list(color="steelblue")) %>%
      add_lines(y=~pen_saved,     name="Penumbra salvaged (mL)",      yaxis="y2",
                line=list(color="orange", dash="dot")) %>%
      layout(title="Benefit vs. Treatment Window",
             xaxis=list(title="Onset-to-treatment (h)"),
             yaxis=list(title="NIHSS improvement"),
             yaxis2=list(title="Penumbra salvaged (mL)", overlaying="y", side="right"))
  })

  # ---- Tab 2: tPA PK ----
  output$p_tpa_pk <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim) || !"Cp_tpa" %in% names(sim)) return(plotly_empty())
    d <- sim %>% filter(time <= 24)
    plot_ly(d, x=~time, y=~Cp_tpa, type="scatter", mode="lines",
            line=list(color="#e41a1c", width=2)) %>%
      layout(title="tPA Plasma Concentration",
             xaxis=list(title="Time (hours)"),
             yaxis=list(title="tPA (ng/mL)"))
  })

  output$p_recanal <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    d <- sim %>% filter(time <= 24)
    plot_ly(d, x=~time) %>%
      add_lines(y=~recanalization, name="Recanalization", line=list(color="#2ca02c")) %>%
      add_lines(y=~THROMBUS,       name="Thrombus burden", line=list(color="#d62728", dash="dot")) %>%
      layout(title="Thrombus Dissolution",
             xaxis=list(title="Time (hours)"),
             yaxis=list(title="Fraction (0–1)", range=c(0,1)))
  })

  output$p_recanal_gauge <- renderPlotly({
    sim <- base_sim()
    val <- if (!is.null(sim)) round(max(sim$recanalization, na.rm=TRUE)*100) else 50
    plot_ly(type="indicator", mode="gauge+number",
            value=val,
            gauge=list(axis=list(range=c(0,100)),
                       bar=list(color="#2ca02c"))) %>%
      layout(margin=list(t=10,b=10))
  })

  # ---- Tab 3: Cascade ----
  output$p_cbf <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    d <- sim %>% filter(time <= 48)
    plot_ly(d, x=~time) %>%
      add_lines(y=~CBF_CORE, name="Core CBF",     line=list(color="#d62728")) %>%
      add_lines(y=~CBF_PEN,  name="Penumbra CBF", line=list(color="#ff7f0e")) %>%
      add_lines(y=rep(55, nrow(d)), x=~time, name="Normal CBF",
                line=list(color="grey50", dash="dot")) %>%
      layout(title="Cerebral Blood Flow (CBF)",
             xaxis=list(title="Time (hours)"),
             yaxis=list(title="CBF (mL/100g/min)"))
  })

  output$p_atp <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    d <- sim %>% filter(time <= 72)
    plot_ly(d, x=~time, y=~ATP_PEN, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(44,160,44,0.2)",
            line=list(color="#2ca02c")) %>%
      layout(title="Penumbral ATP (normalized)",
             xaxis=list(title="Time (hours)"),
             yaxis=list(title="ATP (0–1)", range=c(0,1)))
  })

  output$p_glut_ca2 <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    d <- sim %>% filter(time <= 48)
    plot_ly(d, x=~time) %>%
      add_lines(y=~GLUT,         name="Glutamate (mmol/L)", line=list(color="#d62728")) %>%
      add_lines(y=~CA2*10000,    name="Ca2+ (×10⁻⁴ mmol/L)", yaxis="y2",
                line=list(color="#9467bd", dash="dot")) %>%
      layout(title="Excitotoxicity Cascade",
             xaxis=list(title="Time (hours)"),
             yaxis=list(title="Glutamate (mmol/L)"),
             yaxis2=list(title="Ca2+ (×10⁻⁴)", overlaying="y", side="right"))
  })

  # ---- Tab 4: Inflammation ----
  output$p_inflam <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    d <- sim %>% filter(time <= 336)
    plot_ly(d, x=~time/24) %>%
      add_lines(y=~IL6, name="IL-6 (pg/mL)", line=list(color="#d62728", width=2)) %>%
      add_lines(y=~ROS*50, name="ROS (×50 a.u.)", yaxis="y2",
                line=list(color="#ff7f0e", dash="dot")) %>%
      layout(title="Neuroinflammation Dynamics (0–14 days)",
             xaxis=list(title="Day"),
             yaxis=list(title="IL-6 (pg/mL)"),
             yaxis2=list(title="ROS (a.u.)", overlaying="y", side="right"))
  })

  output$p_bbb <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    d <- sim %>% filter(time <= 168)
    plot_ly(d, x=~time/24, y=~BBB, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(31,119,180,0.2)",
            line=list(color="#1f77b4", width=2)) %>%
      layout(title="BBB Integrity (0–7 days)",
             xaxis=list(title="Day"),
             yaxis=list(title="BBB integrity (0–1)", range=c(0,1)))
  })

  output$p_edema <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    d <- sim %>% filter(time <= 168) %>%
         mutate(edema_risk = (1-BBB) * INFARCT)
    plot_ly(d, x=~time/24, y=~edema_risk, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(214,39,40,0.2)",
            line=list(color="#d62728")) %>%
      layout(title="Edema Risk Index [(1-BBB) × InfarctVol]",
             xaxis=list(title="Day"),
             yaxis=list(title="Edema Risk (a.u.)"))
  })

  # ---- Tab 5: Clinical ----
  output$p_nihss <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    d <- sim %>% filter(time <= 2160)
    plot_ly(d, x=~time/24, y=~NIHSS, type="scatter", mode="lines",
            line=list(color="#e41a1c", width=2)) %>%
      layout(title="NIHSS Trajectory (90 days)",
             xaxis=list(title="Day"), yaxis=list(title="NIHSS Score"))
  })

  output$p_infarct <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    d <- sim %>% filter(time <= 2160)
    plot_ly(d, x=~time/24, y=~INFARCT, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(255,127,14,0.2)",
            line=list(color="#ff7f0e", width=2)) %>%
      layout(title="Infarct Volume (mL) — 90 days",
             xaxis=list(title="Day"), yaxis=list(title="Volume (mL)"))
  })

  output$p_mrs <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    d90 <- sim %>% filter(time == max(sim$time[sim$time <= 2160]))
    val <- round(mean(d90$mRS_est, na.rm=TRUE))
    mrs_df <- tibble(mRS=0:6, prob=c(0,0,0,0,0,0,0))
    mrs_df$prob[val+1] <- 1
    plot_ly(mrs_df, x=~mRS, y=~prob, type="bar",
            marker=list(color=ifelse(mrs_df$mRS==val,"#e41a1c","#aec7e8"))) %>%
      layout(title=paste0("90-Day mRS Estimate: ", val),
             xaxis=list(title="mRS", dtick=1),
             yaxis=list(title="Probability"))
  })

  output$vb_90d_nihss  <- renderValueBox({
    sim <- base_sim()
    val <- if (!is.null(sim)) round(tail(sim$NIHSS[sim$time<=2160],1),1) else "N/A"
    valueBox(val, "90d NIHSS", icon=icon("chart-bar"), color="red")
  })
  output$vb_90d_mrs    <- renderValueBox({
    sim <- base_sim()
    val <- if (!is.null(sim)) round(tail(sim$mRS_est[sim$time<=2160],1),0) else "N/A"
    valueBox(val, "90d mRS", icon=icon("wheelchair"), color="orange")
  })
  output$vb_90d_infvol <- renderValueBox({
    sim <- base_sim()
    val <- if (!is.null(sim)) paste0(round(tail(sim$INFARCT[sim$time<=2160],1),1)," mL") else "N/A"
    valueBox(val, "90d Infarct Vol", icon=icon("brain"), color="maroon")
  })

  # ---- Tab 6: Secondary Prevention ----
  output$p_asp_pk <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    d <- sim %>% filter(time <= 72)
    plot_ly(d, x=~time, y=~Cp_asp, type="scatter", mode="lines",
            line=list(color="#2ca02c", width=2)) %>%
      layout(title="Aspirin Plasma Conc (72h)",
             xaxis=list(title="Hours"), yaxis=list(title="Aspirin (mg/L)"))
  })

  output$p_noac_pk <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    d <- sim %>% filter(time <= 168)
    plot_ly(d, x=~time/24, y=~Cp_noac, type="scatter", mode="lines",
            line=list(color="#1f77b4", width=2)) %>%
      layout(title="Apixaban Plasma Conc (7 days)",
             xaxis=list(title="Day"), yaxis=list(title="Apixaban (ng/mL)"))
  })

  output$p_cox1 <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    d <- sim %>% filter(time <= 72)
    plot_ly(d, x=~time, y=~COX1_inh*100, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(44,160,44,0.2)",
            line=list(color="#2ca02c", width=2)) %>%
      layout(title="COX-1 Inhibition (Aspirin)",
             xaxis=list(title="Hours"), yaxis=list(title="% Inhibition", range=c(0,100)))
  })

  output$p_xa <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    d <- sim %>% filter(time <= 168)
    plot_ly(d, x=~time/24, y=~Xa_inh*100, type="scatter", mode="lines",
            fill="tozeroy", fillcolor="rgba(31,119,180,0.2)",
            line=list(color="#1f77b4", width=2)) %>%
      layout(title="Factor Xa Inhibition (Apixaban)",
             xaxis=list(title="Day"), yaxis=list(title="% Inhibition", range=c(0,100)))
  })

  # ---- Tab 7: Scenario Comparison ----
  scen_data <- eventReactive(input$run_compare, {
    bw <- input$bw
    dose_tpa <- bw * 0.9
    base_p <- list()

    sims <- list()
    if ("Standard tPA 2h" %in% input$scen_sel)
      sims[["Standard tPA 2h"]] <- run_sim(base_p,
        list(tpa_ev(2, dose_tpa), asp_ev_fn(24,100,90)), 2160)

    if ("Late tPA 4.5h" %in% input$scen_sel)
      sims[["Late tPA 4.5h"]] <- run_sim(base_p,
        list(tpa_ev(4.5, dose_tpa), asp_ev_fn(24,100,90)), 2160)

    if ("Antiplatelet Only" %in% input$scen_sel)
      sims[["Antiplatelet Only"]] <- run_sim(base_p,
        list(asp_ev_fn(0,100,90)), 2160)

    if ("tPA + NOAC" %in% input$scen_sel)
      sims[["tPA + NOAC"]] <- run_sim(base_p,
        list(tpa_ev(2, dose_tpa), noac_ev_fn(48, 5, 88)), 2160)

    if ("EVT @ 3h" %in% input$scen_sel) {
      mod_evt <- param(mod, k_tpa_fibrinol=60, k_spont_lysis=1.0)
      sims[["EVT @ 3h"]] <- tryCatch(
        mrgsim(mod_evt, events=asp_ev_fn(24,100,90), end=2160, delta=0.5) %>%
          as_tibble(), error=function(e) NULL)
    }

    bind_rows(lapply(names(sims), function(n) sims[[n]] %>% mutate(scenario=n)))
  }, ignoreNULL=FALSE)

  output$p_compare <- renderPlotly({
    df <- scen_data()
    if (is.null(df) || nrow(df)==0) return(plotly_empty())
    var <- input$compare_var
    df$yval <- df[[var]]
    p <- ggplot(df, aes(time/24, yval, color=scenario)) +
      geom_line(linewidth=1) +
      labs(x="Day", y=var, color=NULL) + theme_bw()
    ggplotly(p)
  })

  output$compare_table <- renderDT({
    df <- scen_data()
    if (is.null(df) || nrow(df)==0) return(datatable(tibble()))
    df %>%
      filter(time %in% c(0,24,168,720,2160)) %>%
      mutate(Day=time/24) %>%
      group_by(scenario, Day) %>%
      summarise(NIHSS=round(mean(NIHSS),1), mRS=round(mean(mRS_est),2),
                Infarct_mL=round(mean(INFARCT),1), BBB=round(mean(BBB),2),
                .groups="drop") %>%
      datatable(options=list(pageLength=20), rownames=FALSE)
  })

  # ---- Tab 8: Biomarkers ----
  output$p_biomarkers <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    d <- sim %>% filter(time <= 168) %>%
      mutate(
        GFAP_proxy   = IL6 * 0.8 + ROS * 2,
        NSE_proxy    = (1 - ATP_PEN) * 30 + INFARCT * 0.3,
        S100B_proxy  = (1 - BBB) * 25,
        DimDimer_proxy = THROMBUS * 500 + (1 - THROMBUS) * 3000
      )
    plot_ly(d, x=~time/24) %>%
      add_lines(y=~IL6,          name="IL-6 (pg/mL)",        line=list(color="#e41a1c")) %>%
      add_lines(y=~GFAP_proxy,   name="GFAP proxy (ng/mL)",  line=list(color="#ff7f0e")) %>%
      add_lines(y=~NSE_proxy,    name="NSE proxy (ng/mL)",   line=list(color="#2ca02c")) %>%
      add_lines(y=~S100B_proxy,  name="S100β proxy (μg/L)",  line=list(color="#9467bd")) %>%
      layout(title="Biomarker Dynamics (0–7 days)",
             xaxis=list(title="Day"),
             yaxis=list(title="Concentration (proxy units)"))
  })

  output$p_biom_detail <- renderPlotly({
    sim <- base_sim()
    if (is.null(sim)) return(plotly_empty())
    var <- input$biom_var
    d <- sim %>% filter(time <= 336)
    plot_ly(d, x=~time/24, y=d[[var]], type="scatter", mode="lines",
            line=list(color="#e41a1c", width=2)) %>%
      layout(xaxis=list(title="Day"), yaxis=list(title=var))
  })
}

# ============================================================
# RUN
# ============================================================
`%||%` <- function(a, b) if (!is.null(a)) a else b

shinyApp(ui = ui, server = server)
