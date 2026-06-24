##############################################################################
#  PAH QSP — Shiny Interactive Dashboard
#  Requires: shiny, mrgsolve, dplyr, ggplot2, patchwork, shinydashboard,
#            DT, plotly, shinyjs
#
#  Run: shiny::runApp("pah_shiny_app.R")
##############################################################################

library(shiny)
library(shinydashboard)
library(shinyjs)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(patchwork)

# ── Inline mrgsolve model (same ODE as in pah_mrgsolve_model.R) ─────────────
PAH_MODEL_CODE <- '
$PROB PAH QSP Shiny Model

$PARAM
ERA_F=0.50 ERA_ka=0.693 ERA_CL=15 ERA_Vc=35 ERA_Q=9 ERA_Vp=50 ERA_ke0=0.3
ERA_IC50=1.0 ERA_Emax=0.85 ERA_n_hill=1.5
PDE5_F=0.80 PDE5_ka=0.35 PDE5_CL=3.5 PDE5_Vc=63 PDE5_ke0=0.2
PDE5_IC50=0.5 PDE5_Emax=0.90 PDE5_n=1.2
PGI2_CL=800 PGI2_Vc=30 PGI2_ke0=2.0 PGI2_EC50=0.1 PGI2_Emax=0.95 PGI2_n=1.0
ET1_PAH=9.0 ET1_kprod=0.15 ET1_kdeg=0.15 ET1_stimHyp=1.8
cGMP_baseline=5.0 cGMP_kprod=1.0 cGMP_kdeg_0=0.8 NO_synth_eff=0.6 sGC_baseline=1.0
cAMP_baseline=4.0 cAMP_kprod=0.9 cAMP_kdeg=0.7 IP_r_eff=1.0
VRI_0=1.0 VRI_kgrowth=0.015 VRI_kdrug_ERA=0.30 VRI_kdrug_PDE=0.20 VRI_kdrug_PGI=0.45 VRI_max=3.0
PVR_normal=100 PVR_PAH0=900 PVR_tone_frac=0.45 PVR_remod_frac=0.55
PAWP_baseline=10 CO_baseline=3.5
Ees_0=0.8 Ees_max=2.5 Ees_kHyp=0.003 Ees_kFail=0.006 RV_fail_thresh=1.5
BNP_normal=20 BNP_kprod=0.1 BNP_kdeg=0.08
sixMWD_normal=580 sixMWD_PAH0=330

$CMT ERA_gut ERA_central ERA_periph ERA_effect PDE5_gut PDE5_central PDE5_effect
     PGI2_central PGI2_effect ET1 cGMP cAMP VRI Ees_RV BNP_conc

$ODE
double ERA_abs   = ERA_ka*ERA_gut;
double ERA_distr = ERA_Q*(ERA_central/ERA_Vc - ERA_periph/ERA_Vp);
double ERA_elim  = (ERA_CL/ERA_Vc)*ERA_central;
dxdt_ERA_gut     = -ERA_abs;
dxdt_ERA_central =  ERA_abs*ERA_F - ERA_distr - ERA_elim;
dxdt_ERA_periph  =  ERA_distr;
dxdt_ERA_effect  =  ERA_ke0*(ERA_central/ERA_Vc - ERA_effect);

double PDE5_abs  = PDE5_ka*PDE5_gut;
double PDE5_elim = (PDE5_CL/PDE5_Vc)*PDE5_central;
dxdt_PDE5_gut     = -PDE5_abs;
dxdt_PDE5_central =  PDE5_abs*PDE5_F - PDE5_elim;
dxdt_PDE5_effect  =  PDE5_ke0*(PDE5_central/PDE5_Vc - PDE5_effect);

double PGI2_elim = (PGI2_CL/PGI2_Vc)*PGI2_central;
dxdt_PGI2_central = -PGI2_elim;
dxdt_PGI2_effect  = PGI2_ke0*(PGI2_central/PGI2_Vc - PGI2_effect);

double ERA_Ce  = ERA_effect;
double ERA_Inh = ERA_Emax*pow(ERA_Ce,ERA_n_hill)/(pow(ERA_IC50,ERA_n_hill)+pow(ERA_Ce,ERA_n_hill));
double PDE5_Inh= PDE5_Emax*pow(PDE5_effect,PDE5_n)/(pow(PDE5_IC50,PDE5_n)+pow(PDE5_effect,PDE5_n));
double PGI2_Act= PGI2_Emax*PGI2_effect/(PGI2_EC50+PGI2_effect);

dxdt_ET1  = ET1_kprod*ET1_stimHyp*(1.0-0.15*ERA_Inh)*(ET1_PAH-ET1) - ET1_kdeg*ET1;
dxdt_cGMP = cGMP_kprod*NO_synth_eff*sGC_baseline - cGMP_kdeg_0*(1.0-PDE5_Inh)*cGMP;
dxdt_cAMP = cAMP_kprod*IP_r_eff*(1.0+PGI2_Act) - cAMP_kdeg*cAMP;

double VRI_drug = VRI_kdrug_ERA*ERA_Inh + VRI_kdrug_PDE*(cGMP/cGMP_baseline-1.0) + VRI_kdrug_PGI*PGI2_Act;
dxdt_VRI  = VRI_kgrowth*VRI*(1.0-VRI/VRI_max) - VRI_drug*VRI;

double PVR_c = PVR_normal+(PVR_PAH0-PVR_normal)*(PVR_tone_frac*(ET1/ET1_PAH)*(1.0-ERA_Inh)*(cGMP_baseline/cGMP)*(cAMP_baseline/cAMP)+PVR_remod_frac*(VRI/VRI_0));
double CO_c  = CO_baseline*(PVR_PAH0/PVR_c);
double mPAP_c= CO_c*PVR_c/80.0+PAWP_baseline;
double Ea_c  = mPAP_c/(CO_c*1000.0/60.0);
double dEes  = (Ea_c < RV_fail_thresh) ? Ees_kHyp*(Ees_max-Ees_RV) : -Ees_kFail*(Ees_RV-Ees_0*0.5);
dxdt_Ees_RV = dEes;
dxdt_BNP_conc = BNP_kprod*(mPAP_c/25.0)*BNP_normal - BNP_kdeg*BNP_conc;

$TABLE
double ERA_ef = ERA_Emax*pow(ERA_effect,ERA_n_hill)/(pow(ERA_IC50,ERA_n_hill)+pow(ERA_effect,ERA_n_hill));
double P5_ef  = PDE5_Emax*pow(PDE5_effect,PDE5_n)/(pow(PDE5_IC50,PDE5_n)+pow(PDE5_effect,PDE5_n));
double PG_ef  = PGI2_Emax*PGI2_effect/(PGI2_EC50+PGI2_effect);
double PVR_s  = PVR_normal+(PVR_PAH0-PVR_normal)*(PVR_tone_frac*(ET1/ET1_PAH)*(1.0-ERA_ef)*(cGMP_baseline/cGMP)*(cAMP_baseline/cAMP)+PVR_remod_frac*(VRI/VRI_0));
double CO_s   = CO_baseline*(PVR_PAH0/PVR_s);
double mPAP_s = CO_s*PVR_s/80.0+PAWP_baseline;
double sixMWD_s = sixMWD_normal-(PVR_s-PVR_normal)*(sixMWD_normal-sixMWD_PAH0)/(PVR_PAH0-PVR_normal);
if(sixMWD_s<50) sixMWD_s=50;
if(sixMWD_s>sixMWD_normal) sixMWD_s=sixMWD_normal;
double WHO;
if(sixMWD_s>440) WHO=1.0;
else if(sixMWD_s>300) WHO=1.0+(440-sixMWD_s)/140.0;
else if(sixMWD_s>150) WHO=2.0+(300-sixMWD_s)/150.0;
else WHO=3.0+(150-sixMWD_s)/150.0;
if(WHO>4.0) WHO=4.0;
double Ea_s  = mPAP_s/(CO_s*1000.0/60.0);
double coup  = Ees_RV/Ea_s;
capture PVR_dyn=PVR_s; capture mPAP_mmHg=mPAP_s; capture CO_Lmin=CO_s;
capture sixMWD_m=sixMWD_s; capture WHO_FC=WHO; capture BNP_pg=BNP_conc;
capture ET1_pg=ET1; capture cGMP_nM=cGMP; capture cAMP_nM=cAMP;
capture VRI_idx=VRI; capture Ees_mmHg=Ees_RV; capture Ea_mmHg=Ea_s;
capture RV_PA_coup=coup;
capture ERA_Inh_pct=ERA_ef*100; capture PDE5_Inh_pct=P5_ef*100; capture PGI2_Act_pct=PG_ef*100;
'

# Compile model on app startup
mod_global <- mread_cache("pah_shiny", tempdir(), PAH_MODEL_CODE)

# ─────────────────────────────────────────────────────────────────────────────
#  HELPER: Run simulation for given parameters
# ─────────────────────────────────────────────────────────────────────────────
run_sim <- function(
    mod,
    era_dose   = 125,   # mg BID
    pde5_dose  = 40,    # mg QD
    pgi2_rate  = 0,     # ng/kg/min
    bw         = 70,
    duration_wk= 12,
    vri_init   = 1.0,   # baseline remodelling severity (1=mild, 2=moderate, 3=severe)
    bnp_init   = 200
) {
  duration_h <- duration_wk * 7 * 24
  events <- list()

  if (era_dose > 0) {
    events$era <- ev(cmt="ERA_gut", amt=era_dose*1000, ii=12,
                     addl=floor(duration_h/12)-1)
  }
  if (pde5_dose > 0) {
    events$pde5 <- ev(cmt="PDE5_gut", amt=pde5_dose*1000, ii=24,
                      addl=floor(duration_h/24)-1)
  }
  if (pgi2_rate > 0) {
    total_ng_h <- pgi2_rate * bw * 60
    events$pgi2 <- ev(cmt="PGI2_central", amt=total_ng_h*duration_h,
                      rate=total_ng_h, time=0)
  }

  ev_combined <- if (length(events) > 0) Reduce(c, events) else ev(time=0, amt=0, cmt=1)

  mod %>%
    init(ET1=9, cGMP=5, cAMP=4, VRI=vri_init, Ees_RV=0.8, BNP_conc=bnp_init) %>%
    ev(ev_combined) %>%
    mrgsim(end=duration_h, delta=24, obsonly=TRUE) %>%
    as.data.frame() %>%
    mutate(time_days = time/24, time_weeks = time/(24*7))
}

# ─────────────────────────────────────────────────────────────────────────────
#  UI
# ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = "PAH QSP Model",
    titleWidth = 260
  ),

  dashboardSidebar(
    width = 260,
    useShinyjs(),
    sidebarMenu(
      id = "sidebar",
      menuItem("Overview",         tabName = "tab_overview",  icon = icon("info-circle")),
      menuItem("Patient Profile",  tabName = "tab_patient",   icon = icon("user")),
      menuItem("Drug Dosing",      tabName = "tab_dosing",    icon = icon("pills")),
      menuItem("Simulation",       tabName = "tab_sim",       icon = icon("chart-line")),
      menuItem("PK Profiles",      tabName = "tab_pk",        icon = icon("wave-square")),
      menuItem("Dose-Response",    tabName = "tab_dr",        icon = icon("sliders-h")),
      menuItem("Risk Assessment",  tabName = "tab_risk",      icon = icon("heartbeat")),
      menuItem("Mechanistic Map",  tabName = "tab_map",       icon = icon("project-diagram"))
    ),
    hr(),
    tags$div(
      style = "padding: 10px; font-size: 11px; color: #aaa;",
      "PAH QSP v1.0", br(),
      "mrgsolve + Shiny"
    )
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background-color: #f4f7fb; }
        .box { border-radius: 8px; }
        .value-box .inner h3 { font-size: 22px; }
        .kpi-label { font-size: 12px; color: #888; }
      "))
    ),

    tabItems(

      # ── TAB 1: OVERVIEW ─────────────────────────────────────────────────
      tabItem(tabName = "tab_overview",
        h2("Pulmonary Arterial Hypertension — QSP Model"),
        fluidRow(
          box(
            width = 12, title = "Model Description", status = "primary",
            solidHeader = TRUE,
            p("This interactive QSP dashboard simulates the pharmacology of",
              "pulmonary arterial hypertension (PAH) using a multi-scale",
              "mechanistic model. Key pathways modelled include:"),
            tags$ul(
              tags$li(strong("Endothelin-1 (ET-1) axis:"),
                      " ETA/ETB receptor signalling → ROCK → vasoconstriction"),
              tags$li(strong("NO-cGMP-PDE5 axis:"),
                      " eNOS → NO → sGC → cGMP → PKG → vasodilation"),
              tags$li(strong("PGI2-cAMP axis:"),
                      " Prostacyclin synthase → IP receptor → AC → cAMP → vasodilation"),
              tags$li(strong("Vascular remodelling:"),
                      " BMPR2/TGF-β/PDGF signalling → PASMC proliferation → VRI index"),
              tags$li(strong("Right ventricular function:"),
                      " RV-PA coupling (Ees/Ea), BNP, TAPSE"),
              tags$li(strong("Clinical endpoints:"),
                      " 6MWD, WHO functional class, mPAP, PVR")
            ),
            hr(),
            h4("Drug classes simulated:"),
            tags$table(
              class = "table table-bordered table-sm",
              tags$thead(tags$tr(
                tags$th("Class"), tags$th("Drugs"), tags$th("Target"), tags$th("Key Trial")
              )),
              tags$tbody(
                tags$tr(tags$td("ERA"),  tags$td("Bosentan, Ambrisentan, Macitentan"),
                        tags$td("ETA/ETB receptor block"), tags$td("SERAPHIN, ARIES")),
                tags$tr(tags$td("PDE5i"), tags$td("Sildenafil, Tadalafil"),
                        tags$td("PDE5 inhibition → ↑cGMP"), tags$td("SUPER-1, PHIRST")),
                tags$tr(tags$td("sGC"),  tags$td("Riociguat"),
                        tags$td("sGC stimulation → ↑cGMP"), tags$td("PATENT-1")),
                tags$tr(tags$td("PGI2"), tags$td("Epoprostenol, Treprostinil, Selexipag"),
                        tags$td("IP receptor agonism → ↑cAMP"), tags$td("GRIPHON"))
              )
            )
          )
        )
      ),

      # ── TAB 2: PATIENT PROFILE ──────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        h3("Patient Profile & Baseline"),
        fluidRow(
          box(
            width = 4, title = "Demographics", status = "info", solidHeader = TRUE,
            numericInput("bw",   "Body weight (kg)",  value = 70,  min = 30, max = 150),
            numericInput("age",  "Age (years)",        value = 45,  min = 18, max = 90),
            selectInput("sex",   "Sex",                choices = c("Female", "Male")),
            selectInput("who_fc_base", "Baseline WHO FC",
                        choices = c("FC I"="1","FC II"="2","FC III"="3","FC IV"="4"),
                        selected = "3")
          ),
          box(
            width = 4, title = "Haemodynamics at Diagnosis", status = "warning",
            solidHeader = TRUE,
            numericInput("base_mPAP", "Baseline mPAP (mmHg)", value = 45, min = 25, max = 100),
            numericInput("base_PVR",  "Baseline PVR (dyn·s·cm⁻⁵)", value = 900, min = 200, max = 3000),
            numericInput("base_CO",   "Baseline CO (L/min)",  value = 3.5, min = 1, max = 8, step = 0.1),
            numericInput("base_6MWD", "Baseline 6MWD (m)",   value = 330, min = 50, max = 580)
          ),
          box(
            width = 4, title = "Biomarkers", status = "danger", solidHeader = TRUE,
            numericInput("base_BNP",  "Baseline BNP (pg/mL)", value = 200, min = 5, max = 5000),
            numericInput("base_ET1",  "Baseline ET-1 (pg/mL)",value = 9.0, min = 1, max = 25),
            numericInput("vri_init",  "Remodelling severity\n(VRI: 0.5 mild → 3 severe)",
                         value = 1.0, min = 0.1, max = 3.0, step = 0.1),
            selectInput("disease_sub", "PAH Subtype",
                        choices = c("Idiopathic (IPAH)",
                                    "Heritable (BMPR2 mutation)",
                                    "CTD-associated",
                                    "CHD-associated",
                                    "Drug/toxin-induced"))
          )
        )
      ),

      # ── TAB 3: DRUG DOSING ──────────────────────────────────────────────
      tabItem(tabName = "tab_dosing",
        h3("Drug Dosing & Regimen"),
        fluidRow(
          box(
            width = 4, title = "ERA (Endothelin Receptor Antagonist)",
            status = "primary", solidHeader = TRUE,
            checkboxInput("use_era", "Include ERA", value = TRUE),
            conditionalPanel("input.use_era",
              selectInput("era_drug", "Drug",
                          choices = c("Bosentan (dual ERA)" = "bosentan",
                                      "Ambrisentan (ETA-selective)" = "ambrisentan",
                                      "Macitentan (tissue-penetrant)" = "macitentan")),
              sliderInput("era_dose", "Dose (mg BID)",
                          min = 0, max = 250, value = 125, step = 31.25),
              tags$small("Target trough: 0.1–10 ng/mL | IC50 ≈ 1 ng/mL (ETA)")
            )
          ),
          box(
            width = 4, title = "PDE5 Inhibitor / sGC Stimulator",
            status = "success", solidHeader = TRUE,
            checkboxInput("use_pde5", "Include PDE5i/sGC", value = TRUE),
            conditionalPanel("input.use_pde5",
              selectInput("pde5_drug", "Drug",
                          choices = c("Tadalafil 40 mg QD" = "tadalafil",
                                      "Sildenafil 20 mg TID" = "sildenafil",
                                      "Riociguat 2.5 mg TID (sGC)" = "riociguat")),
              sliderInput("pde5_dose", "Daily dose (mg)",
                          min = 0, max = 120, value = 40, step = 10),
              tags$small("Tadalafil IC50 ≈ 0.5 ng/mL | PDE5 inhibition ↑ cGMP")
            )
          ),
          box(
            width = 4, title = "Prostacyclin Analogue",
            status = "danger", solidHeader = TRUE,
            checkboxInput("use_pgi2", "Include Prostacyclin", value = FALSE),
            conditionalPanel("input.use_pgi2",
              selectInput("pgi2_drug", "Drug",
                          choices = c("Epoprostenol IV (t½ ~3 min)" = "epoprostenol",
                                      "Treprostinil SC/IV (t½ ~4 h)" = "treprostinil",
                                      "Selexipag oral (IP agonist)" = "selexipag")),
              sliderInput("pgi2_rate", "Infusion rate (ng/kg/min)",
                          min = 0, max = 20, value = 2, step = 0.5),
              tags$small("Starting dose 2 ng/kg/min; titrate by 1–2 ng/kg/min q2wk")
            )
          )
        ),
        fluidRow(
          box(
            width = 6, title = "Simulation Duration", status = "warning",
            solidHeader = TRUE,
            sliderInput("duration_wk", "Duration (weeks)", min = 4, max = 52, value = 12),
            actionButton("run_sim", "Run Simulation", icon = icon("play"),
                         class = "btn-primary btn-lg", width = "100%")
          ),
          box(
            width = 6, title = "Current Regimen", status = "info",
            solidHeader = TRUE,
            verbatimTextOutput("regimen_summary")
          )
        )
      ),

      # ── TAB 4: SIMULATION RESULTS ────────────────────────────────────────
      tabItem(tabName = "tab_sim",
        h3("Simulation Results"),
        fluidRow(
          valueBoxOutput("kpi_pvr",   width = 3),
          valueBoxOutput("kpi_mPAP",  width = 3),
          valueBoxOutput("kpi_6MWD",  width = 3),
          valueBoxOutput("kpi_BNP",   width = 3)
        ),
        fluidRow(
          box(
            width = 6, title = "Pulmonary Vascular Resistance (PVR)",
            status = "danger", solidHeader = TRUE,
            plotlyOutput("plot_PVR", height = "280px")
          ),
          box(
            width = 6, title = "Mean Pulmonary Artery Pressure (mPAP)",
            status = "warning", solidHeader = TRUE,
            plotlyOutput("plot_mPAP", height = "280px")
          )
        ),
        fluidRow(
          box(
            width = 6, title = "6-Minute Walk Distance (6MWD)",
            status = "success", solidHeader = TRUE,
            plotlyOutput("plot_6MWD", height = "280px")
          ),
          box(
            width = 6, title = "BNP / NT-proBNP",
            status = "primary", solidHeader = TRUE,
            plotlyOutput("plot_BNP", height = "280px")
          )
        ),
        fluidRow(
          box(
            width = 6, title = "RV-PA Coupling (Ees/Ea ratio)",
            status = "danger", solidHeader = TRUE,
            plotlyOutput("plot_coupling", height = "280px")
          ),
          box(
            width = 6, title = "Vascular Remodelling Index (VRI)",
            status = "warning", solidHeader = TRUE,
            plotlyOutput("plot_VRI", height = "280px")
          )
        ),
        fluidRow(
          box(
            width = 12, title = "Full Results Table", status = "info",
            solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
            DTOutput("sim_table")
          )
        )
      ),

      # ── TAB 5: PK PROFILES ──────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        h3("Drug PK Profiles"),
        fluidRow(
          box(
            width = 4, title = "ERA Effect-Site Concentration",
            status = "primary", solidHeader = TRUE,
            plotlyOutput("plot_pk_era", height = "300px"),
            tags$small("Blue = plasma | Red = effect-site | Dashed = IC50")
          ),
          box(
            width = 4, title = "PDE5i Effect-Site Concentration",
            status = "success", solidHeader = TRUE,
            plotlyOutput("plot_pk_pde5", height = "300px")
          ),
          box(
            width = 4, title = "Prostacyclin Effect-Site",
            status = "danger", solidHeader = TRUE,
            plotlyOutput("plot_pk_pgi2", height = "300px")
          )
        ),
        fluidRow(
          box(
            width = 6, title = "Second Messenger Dynamics",
            status = "warning", solidHeader = TRUE,
            plotlyOutput("plot_2nd_mess", height = "300px")
          ),
          box(
            width = 6, title = "ET-1 Dynamics",
            status = "danger", solidHeader = TRUE,
            plotlyOutput("plot_ET1", height = "300px")
          )
        )
      ),

      # ── TAB 6: DOSE-RESPONSE ────────────────────────────────────────────
      tabItem(tabName = "tab_dr",
        h3("Dose-Response Analysis"),
        fluidRow(
          box(
            width = 4, title = "DR Parameters", status = "primary",
            solidHeader = TRUE,
            selectInput("dr_drug", "Drug to vary",
                        choices = c("ERA dose (mg BID)" = "era",
                                    "PDE5i dose (mg/day)" = "pde5",
                                    "PGI2 rate (ng/kg/min)" = "pgi2")),
            sliderInput("dr_wk", "Evaluation time (weeks)", min=4, max=24, value=12),
            actionButton("run_dr", "Run DR Analysis", icon=icon("chart-bar"),
                         class="btn-warning btn-lg", width="100%")
          ),
          box(
            width = 8, title = "Dose-Response: 6MWD & PVR",
            status = "success", solidHeader = TRUE,
            plotlyOutput("plot_DR", height = "380px")
          )
        )
      ),

      # ── TAB 7: RISK ASSESSMENT ──────────────────────────────────────────
      tabItem(tabName = "tab_risk",
        h3("ESC/ERS PAH Risk Assessment (Adapted)"),
        fluidRow(
          box(
            width = 6, title = "Risk Score at Baseline vs End of Simulation",
            status = "danger", solidHeader = TRUE,
            tableOutput("risk_table")
          ),
          box(
            width = 6, title = "Risk Score Interpretation",
            status = "info", solidHeader = TRUE,
            tags$table(
              class = "table table-bordered",
              tags$thead(tags$tr(
                tags$th("Parameter"), tags$th("Low risk"), tags$th("Intermediate"), tags$th("High risk")
              )),
              tags$tbody(
                tags$tr(tags$td("WHO-FC"),        tags$td("I-II"),   tags$td("III"),   tags$td("IV")),
                tags$tr(tags$td("6MWD (m)"),      tags$td(">440"),   tags$td("165–440"), tags$td("<165")),
                tags$tr(tags$td("BNP (pg/mL)"),   tags$td("<50"),    tags$td("50–300"), tags$td(">300")),
                tags$tr(tags$td("mPAP (mmHg)"),   tags$td("<25"),    tags$td("25–38"), tags$td(">38")),
                tags$tr(tags$td("PVR (WU)"),       tags$td("<4"),     tags$td("4–8"),   tags$td(">8")),
                tags$tr(tags$td("Ees/Ea"),         tags$td(">1.5"),   tags$td("1–1.5"), tags$td("<1"))
              )
            )
          )
        )
      ),

      # ── TAB 8: MECHANISTIC MAP ──────────────────────────────────────────
      tabItem(tabName = "tab_map",
        h3("PAH Mechanistic Map"),
        fluidRow(
          box(
            width = 12, title = "QSP Pathway Map (>100 nodes)",
            status = "primary", solidHeader = TRUE,
            tags$img(
              src = "pah_qsp_model.svg",
              style = "width:100%; border:1px solid #ddd; border-radius:4px;"
            ),
            hr(),
            downloadButton("dl_svg",  "Download SVG"),
            downloadButton("dl_png",  "Download PNG"),
            downloadButton("dl_dot",  "Download DOT source")
          )
        )
      )

    )  # end tabItems
  )  # end dashboardBody
)  # end dashboardPage


# ─────────────────────────────────────────────────────────────────────────────
#  SERVER
# ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Regimen summary ────────────────────────────────────────────────────────
  output$regimen_summary <- renderText({
    lines <- c("=== Current Regimen ===")
    if (isTRUE(input$use_era))
      lines <- c(lines, sprintf("ERA: %s %.0f mg BID", input$era_drug, input$era_dose))
    if (isTRUE(input$use_pde5))
      lines <- c(lines, sprintf("PDE5i: %s %.0f mg/day", input$pde5_drug, input$pde5_dose))
    if (isTRUE(input$use_pgi2))
      lines <- c(lines, sprintf("PGI2: %s %.1f ng/kg/min", input$pgi2_drug, input$pgi2_rate))
    if (length(lines) == 1) lines <- c(lines, "[No drugs selected]")
    paste(lines, collapse = "\n")
  })

  # ── Reactive: run simulation ───────────────────────────────────────────────
  sim_results <- eventReactive(input$run_sim, {
    withProgress(message = "Running simulation...", value = 0.5, {
      run_sim(
        mod         = mod_global,
        era_dose    = if (input$use_era)  input$era_dose  else 0,
        pde5_dose   = if (input$use_pde5) input$pde5_dose else 0,
        pgi2_rate   = if (input$use_pgi2) input$pgi2_rate else 0,
        bw          = input$bw,
        duration_wk = input$duration_wk,
        vri_init    = input$vri_init,
        bnp_init    = input$base_BNP
      )
    })
  }, ignoreNULL = FALSE)

  # ── KPI value boxes ────────────────────────────────────────────────────────
  output$kpi_pvr <- renderValueBox({
    df <- sim_results()
    val <- round(tail(df$PVR_dyn, 1))
    delta <- round(val - df$PVR_dyn[1])
    valueBox(
      value    = paste0(val, " dyn·s·cm⁻⁵"),
      subtitle = paste0("PVR (Δ", delta, " from baseline)"),
      icon     = icon("tachometer-alt"),
      color    = if (val < 400) "green" else if (val < 700) "yellow" else "red"
    )
  })

  output$kpi_mPAP <- renderValueBox({
    df <- sim_results()
    val <- round(tail(df$mPAP_mmHg, 1), 1)
    valueBox(
      value    = paste0(val, " mmHg"),
      subtitle = "Mean PAP (end of simulation)",
      icon     = icon("heartbeat"),
      color    = if (val < 30) "green" else if (val < 40) "yellow" else "red"
    )
  })

  output$kpi_6MWD <- renderValueBox({
    df <- sim_results()
    val <- round(tail(df$sixMWD_m, 1))
    delta <- round(val - df$sixMWD_m[1])
    valueBox(
      value    = paste0(val, " m"),
      subtitle = paste0("6MWD (Δ+", delta, " m)"),
      icon     = icon("walking"),
      color    = if (val > 440) "green" else if (val > 300) "yellow" else "red"
    )
  })

  output$kpi_BNP <- renderValueBox({
    df <- sim_results()
    val <- round(tail(df$BNP_pg, 1))
    valueBox(
      value    = paste0(val, " pg/mL"),
      subtitle = "BNP (end of simulation)",
      icon     = icon("vial"),
      color    = if (val < 50) "green" else if (val < 300) "yellow" else "red"
    )
  })

  # ── Simulation plots ───────────────────────────────────────────────────────
  make_plotly <- function(df, yvar, ylab, color, ref_line = NULL, ref_label = NULL) {
    p <- plot_ly(df, x = ~time_weeks, y = as.formula(paste0("~", yvar)),
                 type = "scatter", mode = "lines",
                 line = list(color = color, width = 2.5)) %>%
      layout(xaxis = list(title = "Time (weeks)"),
             yaxis = list(title = ylab),
             hovermode = "x unified",
             margin = list(t = 10, b = 40))
    if (!is.null(ref_line))
      p <- p %>%
        add_lines(y = rep(ref_line, nrow(df)), x = df$time_weeks,
                  line = list(dash = "dot", color = "grey", width = 1),
                  name = ref_label)
    p
  }

  output$plot_PVR     <- renderPlotly({ df <- sim_results(); make_plotly(df, "PVR_dyn",   "PVR (dyn·s·cm⁻⁵)", "#e74c3c") })
  output$plot_mPAP    <- renderPlotly({ df <- sim_results(); make_plotly(df, "mPAP_mmHg", "mPAP (mmHg)",       "#f39c12", 25, "PAH threshold") })
  output$plot_6MWD    <- renderPlotly({ df <- sim_results(); make_plotly(df, "sixMWD_m",  "6MWD (m)",          "#27ae60", 440, "FC I threshold") })
  output$plot_BNP     <- renderPlotly({ df <- sim_results(); make_plotly(df, "BNP_pg",    "BNP (pg/mL)",       "#2980b9", 300, "High risk >300") })
  output$plot_coupling<- renderPlotly({ df <- sim_results(); make_plotly(df, "RV_PA_coup","Ees/Ea",            "#8e44ad", 1.0, "Decompensation") })
  output$plot_VRI     <- renderPlotly({ df <- sim_results(); make_plotly(df, "VRI_idx",   "VRI (dimensionless)","#e67e22") })

  # ── PK plots ───────────────────────────────────────────────────────────────
  output$plot_pk_era <- renderPlotly({
    df <- sim_results()
    plot_ly(df, x=~time_weeks) %>%
      add_lines(y=~ERA_Ce_ng, name="ERA effect-site", line=list(color="red", width=2)) %>%
      layout(xaxis=list(title="Time (weeks)"), yaxis=list(title="Conc (ng/mL)"))
  })

  output$plot_pk_pde5 <- renderPlotly({
    df <- sim_results()
    plot_ly(df, x=~time_weeks) %>%
      add_lines(y=~PDE5_Ce_ng, name="PDE5i effect-site", line=list(color="green", width=2)) %>%
      add_lines(y=~PDE5_Inh_pct, name="PDE5 inhibition (%)", yaxis="y2",
                line=list(color="darkgreen", dash="dash")) %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="Conc (ng/mL)"),
             yaxis2=list(overlaying="y", side="right", title="Inhibition (%)"))
  })

  output$plot_pk_pgi2 <- renderPlotly({
    df <- sim_results()
    plot_ly(df, x=~time_weeks) %>%
      add_lines(y=~PGI2_Ce_ng, name="PGI2 effect-site", line=list(color="blue", width=2)) %>%
      add_lines(y=~PGI2_Act_pct, name="IP activation (%)", yaxis="y2",
                line=list(color="navy", dash="dash")) %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="Conc (ng/mL)"),
             yaxis2=list(overlaying="y", side="right", title="Activation (%)"))
  })

  output$plot_2nd_mess <- renderPlotly({
    df <- sim_results()
    plot_ly(df, x=~time_weeks) %>%
      add_lines(y=~cGMP_nM, name="cGMP (nM)", line=list(color="green", width=2)) %>%
      add_lines(y=~cAMP_nM, name="cAMP (nM)", line=list(color="blue", width=2)) %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="Conc (nM)"))
  })

  output$plot_ET1 <- renderPlotly({
    df <- sim_results()
    plot_ly(df, x=~time_weeks) %>%
      add_lines(y=~ET1_pg, name="ET-1 (pg/mL)", line=list(color="red", width=2)) %>%
      add_lines(y=~ERA_Inh_pct, name="ERA inhibition (%)", yaxis="y2",
                line=list(color="purple", dash="dash")) %>%
      layout(xaxis=list(title="Weeks"), yaxis=list(title="ET-1 (pg/mL)"),
             yaxis2=list(overlaying="y", side="right", title="ERA Inhibition (%)"))
  })

  # ── Dose-response ──────────────────────────────────────────────────────────
  dr_results <- eventReactive(input$run_dr, {
    withProgress(message = "Running dose-response...", value = 0, {
      drug <- input$dr_drug
      wk   <- input$dr_wk

      doses <- switch(drug,
        era  = c(0, 31.25, 62.5, 125, 250),
        pde5 = c(0, 10, 20, 40, 80, 120),
        pgi2 = c(0, 0.5, 1, 2, 4, 8, 12)
      )

      lapply(seq_along(doses), function(i) {
        incProgress(1/length(doses))
        d <- doses[i]
        df <- run_sim(
          mod       = mod_global,
          era_dose  = if (drug == "era")  d else 0,
          pde5_dose = if (drug == "pde5") d else 0,
          pgi2_rate = if (drug == "pgi2") d else 0,
          duration_wk = wk,
          vri_init  = input$vri_init,
          bnp_init  = input$base_BNP
        )
        df %>% filter(time == max(time)) %>%
          mutate(dose = d)
      }) %>% bind_rows()
    })
  })

  output$plot_DR <- renderPlotly({
    df <- dr_results()
    subplot(
      plot_ly(df, x=~dose, y=~sixMWD_m, type="scatter", mode="lines+markers",
              name="6MWD (m)", line=list(color="green")) %>%
        layout(xaxis=list(title="Dose"), yaxis=list(title="6MWD (m)")),
      plot_ly(df, x=~dose, y=~PVR_dyn, type="scatter", mode="lines+markers",
              name="PVR", line=list(color="red")) %>%
        layout(xaxis=list(title="Dose"), yaxis=list(title="PVR")),
      plot_ly(df, x=~dose, y=~mPAP_mmHg, type="scatter", mode="lines+markers",
              name="mPAP", line=list(color="orange")) %>%
        layout(xaxis=list(title="Dose"), yaxis=list(title="mPAP (mmHg)")),
      nrows = 1, shareX = FALSE, titleX = TRUE
    ) %>% layout(title = paste("Dose-Response at", input$dr_wk, "weeks"))
  })

  # ── Risk assessment table ──────────────────────────────────────────────────
  output$risk_table <- renderTable({
    df <- sim_results()
    base <- df[1, ]
    last <- df[nrow(df), ]

    risk_classify <- function(val, breaks, labels) {
      cut(val, breaks = c(-Inf, breaks, Inf), labels = labels,
          right = TRUE, include.lowest = TRUE)[1]
    }

    data.frame(
      Parameter = c("mPAP (mmHg)", "PVR (WU)", "6MWD (m)", "BNP (pg/mL)",
                    "WHO-FC", "Ees/Ea"),
      Baseline  = c(
        round(base$mPAP_mmHg, 1),
        round(base$PVR_dyn/80, 1),
        round(base$sixMWD_m),
        round(base$BNP_pg),
        round(base$WHO_FC, 1),
        round(base$RV_PA_coup, 2)
      ),
      `End of Sim` = c(
        round(last$mPAP_mmHg, 1),
        round(last$PVR_dyn/80, 1),
        round(last$sixMWD_m),
        round(last$BNP_pg),
        round(last$WHO_FC, 1),
        round(last$RV_PA_coup, 2)
      ),
      `Risk (end)` = c(
        ifelse(last$mPAP_mmHg < 30, "Low", ifelse(last$mPAP_mmHg < 40, "Intermediate", "High")),
        ifelse(last$PVR_dyn/80 < 4,  "Low", ifelse(last$PVR_dyn/80 < 8,  "Intermediate", "High")),
        ifelse(last$sixMWD_m > 440,  "Low", ifelse(last$sixMWD_m > 165,  "Intermediate", "High")),
        ifelse(last$BNP_pg < 50,     "Low", ifelse(last$BNP_pg < 300,    "Intermediate", "High")),
        ifelse(last$WHO_FC < 2,      "Low", ifelse(last$WHO_FC < 3.5,    "Intermediate", "High")),
        ifelse(last$RV_PA_coup > 1.5,"Low", ifelse(last$RV_PA_coup > 1.0,"Intermediate", "High"))
      ),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # ── Full results table ─────────────────────────────────────────────────────
  output$sim_table <- renderDT({
    df <- sim_results() %>%
      filter(time %% 24 < 1) %>%
      select(time_days, time_weeks, PVR_dyn, mPAP_mmHg, CO_Lmin,
             sixMWD_m, WHO_FC, BNP_pg, ET1_pg, cGMP_nM, cAMP_nM,
             VRI_idx, Ees_mmHg, Ea_mmHg, RV_PA_coup) %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
    datatable(df, options = list(pageLength = 15, scrollX = TRUE),
              class = "compact stripe hover")
  })

  # ── Downloads ──────────────────────────────────────────────────────────────
  output$dl_svg <- downloadHandler(
    filename = "pah_qsp_model.svg",
    content  = function(file) {
      file.copy("pah_qsp_model.svg", file)
    }
  )
  output$dl_png <- downloadHandler(
    filename = "pah_qsp_model.png",
    content  = function(file) {
      file.copy("pah_qsp_model.png", file)
    }
  )
  output$dl_dot <- downloadHandler(
    filename = "pah_qsp_model.dot",
    content  = function(file) {
      file.copy("pah_qsp_model.dot", file)
    }
  )
}

# ── Launch ───────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
