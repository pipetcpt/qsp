## ============================================================================
## AAA QSP Shiny Application
## Abdominal Aortic Aneurysm — Interactive PK/PD Dashboard
## Tabs: Patient Profile · Drug PK · MMP Biomarkers · Aortic Dynamics
##        Scenario Comparison · Risk Assessment
## ============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)
library(shinycssloaders)

## ── Inline Model Code (same as aaa_mrgsolve_model.R) ─────────────────────────
aaa_code <- '
$PARAM
ka_d=0.50, F_d=0.93, CL_d=3.90, Vd_d=110.0, Q_d=15.0, Vp_d=300.0,
kta_d=0.05, kte_d=0.04,
ka_s=0.50, F_s=0.05, CL_s=36.0, Vd_s=245.0, Q_s=8.0, Vp_s=120.0,
ka_bb=0.90, F_bb=0.30, CL_bb=50.0, Vd_bb=300.0, Q_bb=10.0, Vp_bb=200.0,
MAC0=100, k_MAC_in=0.10, k_MAC_out=0.05,
TNF0=10, k_TNF_syn=1.00, k_TNF_deg=0.50,
ROS0=1.00, k_ROS_syn=0.20, k_ROS_deg=0.15,
MMP9_0=1.00, k_MMP9_syn=0.50, k_MMP9_deg=0.30,
Imax_d9=0.85, IC50_d9=0.20, Imax_s9=0.40, IC50_s9=0.05,
MMP2_0=1.00, k_MMP2_syn=0.30, k_MMP2_deg=0.25, Imax_d2=0.70, IC50_d2=0.25,
ELAST0=100, k_ELAST_deg=0.020, k_ELAST_syn=0.001,
COL0=100, k_COL_deg=0.015, k_COL_syn=0.025,
VSMC0=100, k_VSMC_apop=0.008, k_VSMC_prol=0.005,
ILT0=0, k_ILT_grow=0.002, k_ILT_lyse=0.001,
DIAM0=30, k_diam_grow=0.005, ECM_weight=0.60, VSMC_weight=0.20, ILT_weight=0.10,
SBP0=145, Imax_bb=0.20, IC50_bb=0.02, HR0=75, Imax_bb_hr=0.30, IC50_bb_hr=0.015,
ACEi_dose=0.0, Imax_acei=0.18, IC50_acei=0.005

$CMT DGUT DCENT DPERIPH DAORTA SGUT SCENT SPERIPH BGUT BCENT BPERIPH
     MAC TNF ROSO MMP9 MMP2 ELAST COLLAG VSMC ILT DIAM

$MAIN
DGUT_0=0; DCENT_0=0; DPERIPH_0=0; DAORTA_0=0;
SGUT_0=0; SCENT_0=0; SPERIPH_0=0;
BGUT_0=0; BCENT_0=0; BPERIPH_0=0;
MAC_0=MAC0; TNF_0=TNF0; ROSO_0=ROS0;
MMP9_0=MMP9_0; MMP2_0=MMP2_0;
ELAST_0=ELAST0; COLLAG_0=COL0; VSMC_0=VSMC0; ILT_0=ILT0; DIAM_0=DIAM0;

$ODE
double ke_d=CL_d/Vd_d, k12_d=Q_d/Vd_d, k21_d=Q_d/Vp_d;
dxdt_DGUT=-ka_d*DGUT;
dxdt_DCENT=ka_d*F_d*DGUT-(ke_d+k12_d)*DCENT+k21_d*DPERIPH;
dxdt_DPERIPH=k12_d*DCENT-k21_d*DPERIPH;
dxdt_DAORTA=kta_d*(DCENT/Vd_d)-kte_d*DAORTA;
double Cp_doxy=DCENT/Vd_d, Ct_doxy=DAORTA;

double ke_s=CL_s/Vd_s, k12_s=Q_s/Vd_s, k21_s=Q_s/Vp_s;
dxdt_SGUT=-ka_s*SGUT;
dxdt_SCENT=ka_s*F_s*SGUT-(ke_s+k12_s)*SCENT+k21_s*SPERIPH;
dxdt_SPERIPH=k12_s*SCENT-k21_s*SPERIPH;
double Cp_stat=SCENT/Vd_s;

double ke_bb=CL_bb/Vd_bb, k12_bb=Q_bb/Vd_bb, k21_bb=Q_bb/Vp_bb;
dxdt_BGUT=-ka_bb*BGUT;
dxdt_BCENT=ka_bb*F_bb*BGUT-(ke_bb+k12_bb)*BCENT+k21_bb*BPERIPH;
dxdt_BPERIPH=k12_bb*BCENT-k21_bb*BPERIPH;
double Cp_prop=BCENT/Vd_bb;

double SBP_red_bb=Imax_bb*Cp_prop/(IC50_bb+Cp_prop);
double SBP_red_acei=Imax_acei*ACEi_dose/(IC50_acei+ACEi_dose);
double SBP_curr=SBP0*(1.0-SBP_red_bb-SBP_red_acei);
if(SBP_curr<90) SBP_curr=90;
double WALL_STR=SBP_curr*DIAM/(DIAM0*145.0);

double MAC_stim=WALL_STR*ROSO/ROS0;
dxdt_MAC=k_MAC_in*MAC0*MAC_stim-k_MAC_out*MAC;

double statin_NF_inh=0.30*Cp_stat/(0.05+Cp_stat);
dxdt_TNF=k_TNF_syn*(MAC/MAC0)*(1.0-statin_NF_inh)-k_TNF_deg*TNF;

double statin_ROS_inh=0.35*Cp_stat/(0.08+Cp_stat);
dxdt_ROSO=k_ROS_syn*(MAC/MAC0)-k_ROS_deg*ROSO*(1.0+statin_ROS_inh);

double InhMMP9_doxy=Imax_d9*Ct_doxy/(IC50_d9+Ct_doxy);
double InhMMP9_stat=Imax_s9*Cp_stat/(IC50_s9+Cp_stat);
double InhMMP9_total=1.0-(1.0-InhMMP9_doxy)*(1.0-InhMMP9_stat);
dxdt_MMP9=k_MMP9_syn*(MAC/MAC0)*(1.0-InhMMP9_total)-k_MMP9_deg*MMP9;

double InhMMP2_doxy=Imax_d2*Ct_doxy/(IC50_d2+Ct_doxy);
dxdt_MMP2=k_MMP2_syn*(1.0-InhMMP2_doxy)-k_MMP2_deg*MMP2;

double VSMC_apop_rate=k_VSMC_apop*(MMP9/MMP9_0)*(TNF/TNF0)*(ROSO/ROS0);
dxdt_VSMC=k_VSMC_prol*(VSMC/VSMC0)*VSMC0-VSMC_apop_rate*VSMC;
if(VSMC<10) dxdt_VSMC=0;

double elast_deg_rate=k_ELAST_deg*(MMP9/MMP9_0)*(MMP2/MMP2_0);
dxdt_ELAST=k_ELAST_syn*ELAST0-elast_deg_rate*ELAST;
if(ELAST<0) dxdt_ELAST=0;

double col_deg_rate=k_COL_deg*(MMP9/MMP9_0);
dxdt_COLLAG=k_COL_syn*(VSMC/VSMC0)*COL0-col_deg_rate*COLLAG;
if(COLLAG<5) dxdt_COLLAG=0;

double turbulence_factor=pow(DIAM/DIAM0,2.0);
dxdt_ILT=k_ILT_grow*turbulence_factor-k_ILT_lyse*ILT;
if(ILT<0) dxdt_ILT=0;

double ECM_loss_idx=(1.0-ELAST/ELAST0)*ECM_weight;
double VSMC_loss_idx=(1.0-VSMC/VSMC0)*VSMC_weight;
double ILT_contrib=(ILT/50.0)*ILT_weight;
double expansion_driver=1.0+ECM_loss_idx+VSMC_loss_idx+ILT_contrib;
double stress_accel=(WALL_STR>1.2)?(WALL_STR-1.2)*0.5:0.0;
dxdt_DIAM=k_diam_grow*expansion_driver*(1.0+stress_accel);
if(DIAM>120) dxdt_DIAM=0;

$TABLE
double Rupture_P=1.0/(1.0+exp(-(0.12*(DIAM-55.0))));
double SBP_mmHg=SBP_curr;
double Wall_Stress_idx=WALL_STR;

$CAPTURE Cp_doxy Ct_doxy Cp_stat Cp_prop MMP9 MMP2 ELAST COLLAG VSMC ILT DIAM
         TNF ROSO MAC SBP_mmHg Wall_Stress_idx Rupture_P
'
mod_global <- mcode("AAA_shiny", aaa_code, quiet = TRUE)

## ── Color palette ─────────────────────────────────────────────────────────────
scenario_colors <- c(
  "No Treatment"         = "#E41A1C",
  "Doxycycline"          = "#377EB8",
  "Statin"               = "#4DAF4A",
  "Propranolol"          = "#FF7F00",
  "Doxycy + Statin"      = "#984EA3",
  "Triple Therapy"       = "#A65628"
)

## ============================================================================
## UI
## ============================================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "AAA QSP Model"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "profile",    icon = icon("user-md")),
      menuItem("Drug PK",             tabName = "pk",         icon = icon("pills")),
      menuItem("MMP Biomarkers",      tabName = "biomarkers", icon = icon("vials")),
      menuItem("Aortic Wall Dynamics",tabName = "aortic",     icon = icon("heartbeat")),
      menuItem("Scenario Comparison", tabName = "scenarios",  icon = icon("chart-line")),
      menuItem("Rupture Risk",        tabName = "risk",       icon = icon("exclamation-triangle"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f8f9fa; }
      .box { border-radius: 8px; }
    "))),
    tabItems(

      ## ── TAB 1: Patient Profile ──────────────────────────────────────────────
      tabItem(tabName = "profile",
        fluidRow(
          box(title = "Patient Characteristics", width = 4, status = "primary",
            numericInput("DIAM0",    "Initial Aortic Diameter (mm)", 35, 20, 100, 1),
            numericInput("SBP0",     "Baseline Systolic BP (mmHg)", 145, 90, 200, 5),
            numericInput("MMP9_0_i", "Baseline MMP-9 (ng/mL)",      2.5, 0.5, 20, 0.5),
            numericInput("ELAST0_i", "Baseline Elastin (%)",        100, 10, 100, 5),
            numericInput("VSMC0_i",  "Baseline VSMC Density (%)",   100, 10, 100, 5)
          ),
          box(title = "Risk Profile", width = 4, status = "warning",
            checkboxInput("smoke",   "Current/Ex-Smoker",   TRUE),
            checkboxInput("htn",     "Hypertension",        TRUE),
            checkboxInput("dyslip",  "Dyslipidemia",        FALSE),
            checkboxInput("family",  "Family History",      FALSE),
            numericInput("age_pt",   "Age (years)",         68, 40, 95, 1),
            selectInput("sex_pt",    "Sex",  choices = c("Male", "Female"), selected = "Male")
          ),
          box(title = "Patient Risk Summary", width = 4, status = "danger",
            verbatimTextOutput("riskSummary")
          )
        ),
        fluidRow(
          box(title = "AAA Staging Reference", width = 12, status = "info",
            tableOutput("staging_table")
          )
        )
      ),

      ## ── TAB 2: Drug PK ──────────────────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Doxycycline Dosing", width = 3, status = "primary",
            numericInput("doxy_dose", "Dose (mg)", 100, 20, 300, 10),
            numericInput("doxy_int",  "Interval (h)", 12, 8, 24, 4),
            numericInput("doxy_days", "Duration (days)", 2, 1, 7, 1),
            actionButton("run_pk_doxy", "Simulate Doxy PK", class = "btn-primary")
          ),
          box(title = "Statin Dosing", width = 3, status = "success",
            numericInput("stat_dose", "Dose (mg)", 40, 10, 80, 5),
            numericInput("stat_int",  "Interval (h)", 24, 12, 24, 12),
            numericInput("stat_days", "Duration (days)", 3, 1, 7, 1),
            actionButton("run_pk_stat", "Simulate Statin PK", class = "btn-success")
          ),
          box(title = "Propranolol Dosing", width = 3, status = "warning",
            numericInput("bb_dose", "Dose (mg)", 40, 10, 160, 10),
            numericInput("bb_int",  "Interval (h)", 8, 6, 12, 2),
            numericInput("bb_days", "Duration (days)", 2, 1, 7, 1),
            actionButton("run_pk_bb", "Simulate BB PK", class = "btn-warning")
          ),
          box(title = "PK Metrics", width = 3, status = "info",
            verbatimTextOutput("pk_summary")
          )
        ),
        fluidRow(
          box(title = "Doxycycline PK — Plasma vs Aortic Tissue", width = 12,
            withSpinner(plotlyOutput("pk_doxy_plot", height = "350px"))
          )
        ),
        fluidRow(
          box(title = "Statin Plasma PK", width = 6,
            withSpinner(plotlyOutput("pk_stat_plot", height = "300px"))
          ),
          box(title = "Propranolol Plasma PK", width = 6,
            withSpinner(plotlyOutput("pk_bb_plot", height = "300px"))
          )
        )
      ),

      ## ── TAB 3: MMP Biomarkers ───────────────────────────────────────────────
      tabItem(tabName = "biomarkers",
        fluidRow(
          box(title = "Doxycycline Settings", width = 3, status = "primary",
            numericInput("bm_doxy_dose", "Doxy Dose (mg)", 100, 0, 200, 25),
            numericInput("bm_stat_dose", "Statin Dose (mg)", 40, 0, 80, 10),
            numericInput("bm_sim_days",  "Simulation Days", 180, 30, 730, 30),
            actionButton("run_bm", "Run Biomarker Sim", class = "btn-primary btn-block")
          ),
          box(title = "Biomarker Targets", width = 9,
            fluidRow(
              valueBoxOutput("mmp9_box", width = 4),
              valueBoxOutput("mmp2_box", width = 4),
              valueBoxOutput("tnf_box",  width = 4)
            )
          )
        ),
        fluidRow(
          box(title = "MMP-9 Activity Over Time", width = 6,
            withSpinner(plotlyOutput("mmp9_plot", height = "300px"))
          ),
          box(title = "MMP-2 Activity Over Time", width = 6,
            withSpinner(plotlyOutput("mmp2_plot", height = "300px"))
          )
        ),
        fluidRow(
          box(title = "TNF-α and ROS Dynamics", width = 6,
            withSpinner(plotlyOutput("inflam_plot", height = "300px"))
          ),
          box(title = "Macrophage Activity Index", width = 6,
            withSpinner(plotlyOutput("mac_plot", height = "300px"))
          )
        )
      ),

      ## ── TAB 4: Aortic Wall Dynamics ─────────────────────────────────────────
      tabItem(tabName = "aortic",
        fluidRow(
          box(title = "Treatment Settings", width = 3, status = "primary",
            sliderInput("aw_doxy", "Doxy Dose (mg BID)", 0, 200, 100, 25),
            sliderInput("aw_stat", "Statin Dose (mg QD)", 0, 80, 40, 10),
            sliderInput("aw_bb",   "Propranolol (mg TID)", 0, 160, 40, 20),
            sliderInput("aw_days", "Duration (days)", 90, 730, 365, 30),
            actionButton("run_aortic", "Simulate", class = "btn-primary btn-block")
          ),
          box(title = "Aortic Diameter Trajectory", width = 9,
            withSpinner(plotlyOutput("diam_plot", height = "350px"))
          )
        ),
        fluidRow(
          box(title = "Elastin Content", width = 4,
            withSpinner(plotlyOutput("elast_plot", height = "280px"))
          ),
          box(title = "Collagen Content", width = 4,
            withSpinner(plotlyOutput("collag_plot", height = "280px"))
          ),
          box(title = "VSMC Density & ILT Volume", width = 4,
            withSpinner(plotlyOutput("vsmc_ilt_plot", height = "280px"))
          )
        )
      ),

      ## ── TAB 5: Scenario Comparison ──────────────────────────────────────────
      tabItem(tabName = "scenarios",
        fluidRow(
          box(title = "Simulation Settings", width = 3, status = "primary",
            numericInput("sc_diam0", "Initial Diameter (mm)", 35, 20, 60, 1),
            numericInput("sc_sbp0",  "Baseline SBP (mmHg)", 145, 110, 180, 5),
            numericInput("sc_days",  "Simulation (days)", 365, 180, 730, 30),
            checkboxGroupInput("sc_scenarios", "Scenarios:",
              choices = c("No Treatment", "Doxycycline", "Statin",
                          "Propranolol", "Doxycy + Statin", "Triple Therapy"),
              selected = c("No Treatment", "Doxycycline", "Statin", "Triple Therapy")
            ),
            actionButton("run_scenarios", "Run Comparison", class = "btn-primary btn-block")
          ),
          box(title = "Diameter Comparison", width = 9,
            withSpinner(plotlyOutput("sc_diam_plot", height = "350px"))
          )
        ),
        fluidRow(
          box(title = "MMP-9 Comparison", width = 6,
            withSpinner(plotlyOutput("sc_mmp9_plot", height = "280px"))
          ),
          box(title = "Rupture Risk Comparison", width = 6,
            withSpinner(plotlyOutput("sc_risk_plot", height = "280px"))
          )
        ),
        fluidRow(
          box(title = "1-Year Summary Table", width = 12,
            withSpinner(DTOutput("sc_table"))
          )
        )
      ),

      ## ── TAB 6: Rupture Risk Assessment ──────────────────────────────────────
      tabItem(tabName = "risk",
        fluidRow(
          box(title = "Current Patient State", width = 4, status = "danger",
            numericInput("rk_diam",  "Current Diameter (mm)", 45, 20, 120, 1),
            numericInput("rk_sbp",   "Systolic BP (mmHg)", 145, 90, 200, 5),
            numericInput("rk_mmp9",  "Plasma MMP-9 (ng/mL)", 5, 0.5, 50, 0.5),
            numericInput("rk_elast", "Elastin Integrity (%)", 60, 0, 100, 5),
            numericInput("rk_ilt",   "ILT Volume (mL)", 10, 0, 100, 5),
            actionButton("calc_risk", "Calculate Risk", class = "btn-danger btn-block")
          ),
          box(title = "Rupture Risk Estimate", width = 4, status = "danger",
            valueBoxOutput("rupture_prob_box", width = 12),
            br(),
            htmlOutput("risk_interpretation")
          ),
          box(title = "Annual Rupture Risk by Diameter", width = 4, status = "info",
            tableOutput("rupture_ref_table")
          )
        ),
        fluidRow(
          box(title = "Diameter-Based Rupture Risk Curve", width = 6,
            withSpinner(plotlyOutput("risk_curve_plot", height = "320px"))
          ),
          box(title = "5-Year Projection under Treatment", width = 6,
            withSpinner(plotlyOutput("longterm_plot", height = "320px"))
          )
        ),
        fluidRow(
          box(title = "Clinical Recommendations", width = 12, status = "warning",
            htmlOutput("clinical_recommendations")
          )
        )
      )

    )  # end tabItems
  )
)

## ============================================================================
## SERVER
## ============================================================================
server <- function(input, output, session) {

  ## ── Reactive: simulation function ──────────────────────────────────────────
  run_sim <- function(doxy_dose = 0, stat_dose = 0, bb_dose = 0,
                      diam0 = 35, sbp0 = 145, mmp9_0 = 1, elast0 = 100, vsmc0 = 100,
                      days = 365) {
    end_h <- days * 24
    events <- NULL
    if (doxy_dose > 0) {
      ev_d <- ev(cmt = "DGUT", amt = doxy_dose, ii = 12,
                 addl = days * 2 - 1, time = 0)
      events <- if (is.null(events)) ev_d else events + ev_d
    }
    if (stat_dose > 0) {
      ev_s <- ev(cmt = "SGUT", amt = stat_dose, ii = 24,
                 addl = days - 1, time = 0)
      events <- if (is.null(events)) ev_s else events + ev_s
    }
    if (bb_dose > 0) {
      ev_b <- ev(cmt = "BGUT", amt = bb_dose, ii = 8,
                 addl = days * 3 - 1, time = 0)
      events <- if (is.null(events)) ev_b else events + ev_b
    }
    m <- mod_global %>%
      param(DIAM0 = diam0, SBP0 = sbp0, MMP9_0 = mmp9_0,
            ELAST0 = elast0, VSMC0 = vsmc0)
    if (is.null(events)) {
      out <- m %>% mrgsim(end = end_h, delta = 24)
    } else {
      out <- m %>% mrgsim(events = events, end = end_h, delta = 24)
    }
    as_tibble(out) %>% mutate(Day = time / 24)
  }

  ## ── Tab 1: Patient Risk Summary ────────────────────────────────────────────
  output$riskSummary <- renderText({
    score <- 0
    if (input$smoke)   score <- score + 2
    if (input$htn)     score <- score + 2
    if (input$dyslip)  score <- score + 1
    if (input$family)  score <- score + 2
    if (input$age_pt >= 65) score <- score + 1
    if (input$sex_pt == "Male") score <- score + 1
    if (input$DIAM0 >= 55) score <- score + 3
    risk_level <- ifelse(score >= 7, "HIGH", ifelse(score >= 4, "MODERATE", "LOW"))
    paste0("Risk Score: ", score, "/12\n",
           "Risk Level: ", risk_level, "\n\n",
           "Recommended Action:\n",
           if (input$DIAM0 >= 55) "  ► Urgent surgical evaluation\n" else
           if (input$DIAM0 >= 45) "  ► 6-monthly surveillance\n" else
           "  ► Annual ultrasound surveillance\n",
           if (score >= 7) "  ► Consider preventive pharmacotherapy\n" else "")
  })

  output$staging_table <- renderTable({
    data.frame(
      Stage       = c("I (Small)", "II (Medium)", "III (Large)", "IV (Very Large)", "V (Giant)"),
      `Diameter (mm)` = c("<40", "40–45", "45–55", "55–70", ">70"),
      `Rupture Risk (%/yr)` = c("<1", "1–3", "3–15", "10–25", ">25"),
      Surveillance = c("Annual US", "6-month US", "3–6 month + CT", "Surgical eval.", "Emergency"),
      Treatment    = c("Risk factor control", "Pharmacotherapy", "Pharmacotherapy + plan surgery",
                       "EVAR or open repair", "Emergency repair"),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  ## ── Tab 2: PK Simulations ──────────────────────────────────────────────────
  pk_doxy_data <- eventReactive(input$run_pk_doxy, {
    ev_d <- ev(cmt = "DGUT", amt = input$doxy_dose, ii = 12,
               addl = input$doxy_days * 2 - 1, time = 0)
    mod_global %>%
      mrgsim(events = ev_d, end = input$doxy_days * 24, delta = 0.5) %>%
      as_tibble()
  }, ignoreNULL = FALSE)

  output$pk_doxy_plot <- renderPlotly({
    df <- pk_doxy_data()
    plot_ly(df) %>%
      add_lines(x = ~time, y = ~Cp_doxy, name = "Plasma (Cp)", line = list(color = "#2C7BB6")) %>%
      add_lines(x = ~time, y = ~Ct_doxy, name = "Aortic Tissue (Ct)", line = list(color = "#D7191C")) %>%
      layout(xaxis = list(title = "Time (h)"), yaxis = list(title = "Doxycycline (mg/L)"),
             legend = list(orientation = "h"))
  })

  output$pk_stat_plot <- renderPlotly({
    ev_s <- ev(cmt = "SGUT", amt = input$stat_dose, ii = 24,
               addl = input$stat_days - 1, time = 0)
    df <- mod_global %>%
      mrgsim(events = ev_s, end = input$stat_days * 24, delta = 0.5) %>%
      as_tibble()
    plot_ly(df, x = ~time, y = ~Cp_stat, type = "scatter", mode = "lines",
            line = list(color = "#4DAF4A")) %>%
      layout(xaxis = list(title = "Time (h)"), yaxis = list(title = "Statin (mg/L)"))
  })

  output$pk_bb_plot <- renderPlotly({
    ev_b <- ev(cmt = "BGUT", amt = input$bb_dose, ii = 8,
               addl = input$bb_days * 3 - 1, time = 0)
    df <- mod_global %>%
      mrgsim(events = ev_b, end = input$bb_days * 24, delta = 0.5) %>%
      as_tibble()
    plot_ly(df, x = ~time, y = ~Cp_prop, type = "scatter", mode = "lines",
            line = list(color = "#FF7F00")) %>%
      layout(xaxis = list(title = "Time (h)"), yaxis = list(title = "Propranolol (mg/L)"))
  })

  output$pk_summary <- renderText({
    df <- pk_doxy_data()
    cmax <- max(df$Cp_doxy, na.rm = TRUE)
    tmax_i <- which.max(df$Cp_doxy)
    tmax <- df$time[tmax_i]
    auc  <- sum(diff(df$time) * (head(df$Cp_doxy, -1) + tail(df$Cp_doxy, -1)) / 2, na.rm = TRUE)
    sprintf("Doxycycline PK:\n Cmax = %.3f mg/L\n Tmax = %.1f h\n AUC0-t = %.2f mg·h/L\n Aortic Ct max = %.4f mg/L",
            cmax, tmax, auc, max(df$Ct_doxy, na.rm = TRUE))
  })

  ## ── Tab 3: Biomarkers ─────────────────────────────────────────────────────
  bm_data <- eventReactive(input$run_bm, {
    run_sim(doxy_dose = input$bm_doxy_dose, stat_dose = input$bm_stat_dose,
            days = input$bm_sim_days)
  }, ignoreNULL = FALSE)

  output$mmp9_box <- renderValueBox({
    df <- bm_data()
    final <- tail(df, 1)
    valueBox(sprintf("%.2f ng/mL", final$MMP9), "MMP-9 at End",
             icon = icon("flask"), color = if (final$MMP9 > 5) "red" else "green")
  })
  output$mmp2_box <- renderValueBox({
    df <- bm_data()
    final <- tail(df, 1)
    valueBox(sprintf("%.2f ng/mL", final$MMP2), "MMP-2 at End",
             icon = icon("flask"), color = if (final$MMP2 > 3) "orange" else "green")
  })
  output$tnf_box <- renderValueBox({
    df <- bm_data()
    final <- tail(df, 1)
    valueBox(sprintf("%.1f pg/mL", final$TNF), "TNF-α at End",
             icon = icon("biohazard"), color = if (final$TNF > 15) "red" else "yellow")
  })

  output$mmp9_plot <- renderPlotly({
    df <- bm_data()
    plot_ly(df, x = ~Day, y = ~MMP9, type = "scatter", mode = "lines",
            line = list(color = "#FF8C00", width = 2)) %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Active MMP-9 (ng/mL)"))
  })
  output$mmp2_plot <- renderPlotly({
    df <- bm_data()
    plot_ly(df, x = ~Day, y = ~MMP2, type = "scatter", mode = "lines",
            line = list(color = "#E41A1C", width = 2)) %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Active MMP-2 (ng/mL)"))
  })
  output$inflam_plot <- renderPlotly({
    df <- bm_data()
    plot_ly(df) %>%
      add_lines(x = ~Day, y = ~TNF, name = "TNF-α (pg/mL)", line = list(color = "#E41A1C")) %>%
      add_lines(x = ~Day, y = ~ROSO * 10, name = "ROS ×10 (nmol/mg)", line = list(color = "#FF8C00")) %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Level"),
             legend = list(orientation = "h"))
  })
  output$mac_plot <- renderPlotly({
    df <- bm_data()
    plot_ly(df, x = ~Day, y = ~MAC, type = "scatter", mode = "lines",
            line = list(color = "#9370DB", width = 2)) %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Macrophage Activity (AU)"))
  })

  ## ── Tab 4: Aortic Wall Dynamics ───────────────────────────────────────────
  aw_data <- eventReactive(input$run_aortic, {
    run_sim(doxy_dose = input$aw_doxy, stat_dose = input$aw_stat, bb_dose = input$aw_bb,
            diam0 = input$sc_diam0, sbp0 = input$sc_sbp0, days = input$aw_days)
  }, ignoreNULL = FALSE)

  output$diam_plot <- renderPlotly({
    df <- aw_data()
    plot_ly(df, x = ~Day, y = ~DIAM, type = "scatter", mode = "lines",
            line = list(color = "#E41A1C", width = 2)) %>%
      add_segments(x = 0, xend = max(df$Day), y = 55, yend = 55,
                   line = list(color = "red", dash = "dash"), name = "Surgical Threshold") %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Aortic Diameter (mm)"))
  })
  output$elast_plot <- renderPlotly({
    df <- aw_data()
    plot_ly(df, x = ~Day, y = ~ELAST, type = "scatter", mode = "lines",
            line = list(color = "#DAA520", width = 2)) %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Elastin Content (%)"))
  })
  output$collag_plot <- renderPlotly({
    df <- aw_data()
    plot_ly(df, x = ~Day, y = ~COLLAG, type = "scatter", mode = "lines",
            line = list(color = "#8B4513", width = 2)) %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Collagen Content (%)"))
  })
  output$vsmc_ilt_plot <- renderPlotly({
    df <- aw_data()
    plot_ly(df) %>%
      add_lines(x = ~Day, y = ~VSMC, name = "VSMC (%)", line = list(color = "#228B22")) %>%
      add_lines(x = ~Day, y = ~ILT, name = "ILT (mL)", line = list(color = "#9370DB")) %>%
      layout(xaxis = list(title = "Day"), yaxis = list(title = "Level"),
             legend = list(orientation = "h"))
  })

  ## ── Tab 5: Scenario Comparison ────────────────────────────────────────────
  sc_data <- eventReactive(input$run_scenarios, {
    sc_map <- list(
      "No Treatment"    = c(0,   0,  0),
      "Doxycycline"     = c(100, 0,  0),
      "Statin"          = c(0,   40, 0),
      "Propranolol"     = c(0,   0,  40),
      "Doxycy + Statin" = c(100, 40, 0),
      "Triple Therapy"  = c(100, 40, 40)
    )
    selected <- input$sc_scenarios
    bind_rows(lapply(selected, function(sc) {
      doses <- sc_map[[sc]]
      run_sim(doxy_dose = doses[1], stat_dose = doses[2], bb_dose = doses[3],
              diam0 = input$sc_diam0, sbp0 = input$sc_sbp0, days = input$sc_days) %>%
        mutate(Scenario = sc)
    }))
  }, ignoreNULL = FALSE)

  output$sc_diam_plot <- renderPlotly({
    df <- sc_data()
    p <- ggplot(df, aes(Day, DIAM, color = Scenario)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 55, linetype = "dashed", color = "red") +
      scale_color_manual(values = scenario_colors) +
      labs(x = "Day", y = "Aortic Diameter (mm)") + theme_bw()
    ggplotly(p)
  })
  output$sc_mmp9_plot <- renderPlotly({
    df <- sc_data()
    p <- ggplot(df, aes(Day, MMP9, color = Scenario)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = scenario_colors) +
      labs(x = "Day", y = "MMP-9 (ng/mL)") + theme_bw()
    ggplotly(p)
  })
  output$sc_risk_plot <- renderPlotly({
    df <- sc_data()
    p <- ggplot(df, aes(Day, Rupture_P * 100, color = Scenario)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = scenario_colors) +
      labs(x = "Day", y = "Rupture Probability (%)") + theme_bw()
    ggplotly(p)
  })
  output$sc_table <- renderDT({
    df <- sc_data() %>%
      group_by(Scenario) %>%
      filter(abs(Day - max(Day)) < 0.1) %>%
      summarise(
        `Diameter (mm)`   = round(mean(DIAM), 2),
        `MMP-9 (ng/mL)`   = round(mean(MMP9), 3),
        `Elastin (%)`     = round(mean(ELAST), 1),
        `VSMC (%)`        = round(mean(VSMC), 1),
        `ILT (mL)`        = round(mean(ILT), 2),
        `Rupture Risk (%)` = round(mean(Rupture_P) * 100, 2),
        `SBP (mmHg)`      = round(mean(SBP_mmHg), 1),
        .groups = "drop"
      )
    datatable(df, options = list(pageLength = 10, dom = "t"), rownames = FALSE) %>%
      formatStyle("Rupture Risk (%)", backgroundColor = styleInterval(c(5, 15), c("lightgreen", "yellow", "salmon")))
  })

  ## ── Tab 6: Risk Assessment ────────────────────────────────────────────────
  risk_val <- eventReactive(input$calc_risk, {
    d <- input$rk_diam
    rp <- 1 / (1 + exp(-(0.12 * (d - 55))))
    rp * 100
  })

  output$rupture_prob_box <- renderValueBox({
    rp <- risk_val()
    color <- if (rp > 20) "red" else if (rp > 5) "orange" else "green"
    valueBox(sprintf("%.1f%%", rp), "Rupture Probability",
             icon = icon("exclamation-triangle"), color = color)
  })

  output$risk_interpretation <- renderUI({
    rp <- risk_val()
    msg <- if (rp > 20) {
      HTML("<div class='alert alert-danger'><b>VERY HIGH RISK:</b> Immediate surgical consultation recommended. Consider emergency EVAR or open repair.</div>")
    } else if (rp > 10) {
      HTML("<div class='alert alert-warning'><b>HIGH RISK:</b> Urgent surgical evaluation. Elective repair indicated if diameter >55mm.</div>")
    } else if (rp > 3) {
      HTML("<div class='alert alert-info'><b>MODERATE RISK:</b> Close surveillance (3–6 months). Optimise pharmacotherapy (antihypertensives, statin).</div>")
    } else {
      HTML("<div class='alert alert-success'><b>LOW RISK:</b> Annual surveillance. Continue risk factor modification.</div>")
    }
    msg
  })

  output$rupture_ref_table <- renderTable({
    data.frame(
      `Diameter (mm)` = c("<40", "40–44", "45–49", "50–54", "55–59", "60–69", "≥70"),
      `Annual Rupture Risk (%)` = c("<1", "1", "1–3", "3–5", "5–10", "10–20", ">20"),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, bordered = TRUE)

  output$risk_curve_plot <- renderPlotly({
    diam_seq <- seq(25, 100, by = 1)
    rp_seq <- 1 / (1 + exp(-(0.12 * (diam_seq - 55)))) * 100
    df_curve <- data.frame(Diameter = diam_seq, Rupture_Prob = rp_seq)
    plot_ly(df_curve, x = ~Diameter, y = ~Rupture_Prob, type = "scatter", mode = "lines",
            line = list(color = "#E41A1C", width = 2.5)) %>%
      add_segments(x = input$rk_diam, xend = input$rk_diam, y = 0, yend = risk_val(),
                   line = list(color = "blue", dash = "dash"), name = "Current Patient") %>%
      add_segments(x = 55, xend = 55, y = 0, yend = 100,
                   line = list(color = "orange", dash = "dot"), name = "Surgical Threshold (55mm)") %>%
      layout(xaxis = list(title = "Aortic Diameter (mm)"),
             yaxis = list(title = "Rupture Probability (%)"),
             legend = list(orientation = "h"))
  })

  output$longterm_plot <- renderPlotly({
    df_no <- run_sim(diam0 = input$rk_diam, sbp0 = input$rk_sbp, days = 5 * 365) %>%
      mutate(Scenario = "No Treatment")
    df_tx <- run_sim(doxy_dose = 100, stat_dose = 40, bb_dose = 40,
                     diam0 = input$rk_diam, sbp0 = input$rk_sbp, days = 5 * 365) %>%
      mutate(Scenario = "Triple Therapy")
    df_all <- bind_rows(df_no, df_tx)
    p <- ggplot(df_all, aes(Day / 365, DIAM, color = Scenario)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 55, linetype = "dashed", color = "red") +
      scale_color_manual(values = c("No Treatment" = "#E41A1C", "Triple Therapy" = "#377EB8")) +
      labs(x = "Year", y = "Aortic Diameter (mm)") + theme_bw()
    ggplotly(p)
  })

  output$clinical_recommendations <- renderUI({
    d <- input$rk_diam
    HTML(paste0(
      "<ul>",
      "<li><b>Surveillance:</b> ", if (d >= 55) "Surgical evaluation" else if (d >= 45) "CT/MRI every 6 months" else "Annual ultrasound", "</li>",
      "<li><b>Blood Pressure:</b> Target SBP &lt;130 mmHg. Use β-blockers (e.g., propranolol 40mg TID) or ACE inhibitors.</li>",
      "<li><b>MMP Inhibition:</b> ", if (d >= 40) "Consider doxycycline 100mg BID (off-label, shown to reduce MMP-9 and elastin degradation in PHAST trial)" else "Doxycycline may be considered in high-MMP patients", "</li>",
      "<li><b>Lipid Control:</b> High-intensity statin therapy (atorvastatin 40–80mg) for pleiotropic anti-inflammatory and antioxidant effects.</li>",
      "<li><b>Smoking Cessation:</b> CRITICAL — smoking accelerates AAA growth by ~0.4 mm/year.</li>",
      "<li><b>Physical Activity:</b> Avoid vigorous exercise (Valsalva maneuver) if diameter &gt;5 cm.</li>",
      "</ul>"
    ))
  })

}

## ============================================================================
## LAUNCH
## ============================================================================
shinyApp(ui, server)
