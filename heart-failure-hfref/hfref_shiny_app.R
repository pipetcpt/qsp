################################################################################
# Heart Failure with Reduced Ejection Fraction (HFrEF)
# Interactive Shiny QSP Dashboard
# 8 Tabs: Patient Profile · PK · RAAS+SNS · NPS+Hemodynamics ·
#         Remodeling · Clinical Endpoints · Scenario Comparison · Biomarker Risk
################################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(scales)

# ─────────────────────────────────────────────────────────────────────────────
# MODEL CODE (embedded for standalone app)
# ─────────────────────────────────────────────────────────────────────────────
hfref_model_code <- '
$PARAM
kAngI_prod=0.12, kAngI_deg=0.05, kACE_Vmax=0.20, kACE_Km=2.0,
kACE2_Vmax=0.04, kACE2_Km=2.0, kAngII_deg=0.30, kAng17_deg=0.25,
EC50_AngII=40.0, kAldo_max=0.60, kAldo_deg=0.18, Aldo_base=8.0,
kNE_prod=180.0, kNE_deg=0.55, NE_base=325.0, NE_HF_factor=2.2,
kBNP_prod=48.0, kBNP_NEPdeg=0.30, kBNP_other=0.12, BNP_base=35.0,
NTproBNP_ratio=8.5, kcGMP_prod=2.5, kcGMP_PDE=0.80, kcGMP_other=0.20,
LVEDV_base=280.0, HR_base=82.0, HR_setpoint=82.0,
SVR_base=1400.0, EF_base=0.27, EF_max=0.65, kSVR_decay=0.015, MAP_target=85.0,
kFib_prod=0.002, kFib_deg=0.0005, Fib_base=0.35, Fib_max=0.85,
kHyp_prod=0.004, kHyp_deg=0.001, Hyp_base=1.45,
kTGFb1_prod=1.5, kTGFb1_deg=0.25, TGFb1_base=12.0,
kTNFa_prod=0.8, kTNFa_deg=0.35, TNFa_base=6.5,
kIL6_prod=1.2, kIL6_deg=0.40, IL6_base=5.8,
CL_LBQ=4.6, Vd_LBQ=26.0, EC50_NEP_inhib=80.0,
CL_Valsar=1.9, Vd_Valsar=75.0, EC50_AT1R_val=200.0,
CL_BB=55.0, Vd_BB=220.0, EC50_BB=50.0, Emax_BB_HR=0.30, Emax_BB_NE=0.45,
CL_MRA=7.0, Vd_MRA=43.0, EC50_MRA=50.0, Emax_MRA_Ald=0.80,
CL_SGLT2=16.1, Vd_SGLT2=73.8, EC50_SGLT2=40.0, Emax_SGLT2_V=0.18, Emax_SGLT2_C=0.12,
CL_IVA=10.0, Vd_IVA=100.0, EC50_IVA=20.0, Emax_IVA_HR=0.25,
kEF_ARNI=0.003, kEF_BB=0.002, kEF_SGLT2=0.0015

$INIT
AngI=2.0, AngII=55.0, Ang17=3.5, Aldo=22.0,
NE=715.0, BNP=380.0, NTpBNP=3230.0, cGMP=12.0,
LVEDV=280.0, HR=82.0, SVR=1400.0, LVEF=0.27,
Fib=0.35, TGFb1=12.0, Hyp=1.45, TNFa=6.5, IL6=5.8,
LBQ_C=0.0, Valsar_C=0.0, BB_C=0.0, MRA_C=0.0, SGLT2_C=0.0, IVA_C=0.0

$ODE
double E_NEP   = LBQ_C    / (LBQ_C    + EC50_NEP_inhib);
double E_AT1R  = Valsar_C  / (Valsar_C  + EC50_AT1R_val);
double E_BB_HR = Emax_BB_HR  * BB_C    / (BB_C    + EC50_BB);
double E_BB_NE = Emax_BB_NE  * BB_C    / (BB_C    + EC50_BB);
double E_MRA   = Emax_MRA_Ald * MRA_C  / (MRA_C   + EC50_MRA);
double E_S_V   = Emax_SGLT2_V * SGLT2_C / (SGLT2_C + EC50_SGLT2);
double E_S_C   = Emax_SGLT2_C * SGLT2_C / (SGLT2_C + EC50_SGLT2);
double E_IVA   = Emax_IVA_HR  * IVA_C  / (IVA_C   + EC50_IVA);

double ACE_r   = kACE_Vmax  * AngI / (kACE_Km  + AngI);
double ACE2_r  = kACE2_Vmax * AngI / (kACE2_Km + AngI) * (1 + E_NEP*0.5);
double AT1R_FB = 1.0 + E_AT1R * 1.5;

dxdt_AngI   = kAngI_prod - ACE_r - ACE2_r - kAngI_deg*AngI;
dxdt_AngII  = ACE_r * AT1R_FB - kAngII_deg*AngII;
dxdt_Ang17  = ACE2_r - kAng17_deg*Ang17*(1.0-E_NEP*0.4);
dxdt_Aldo   = kAldo_max*pow(AngII,1.2)/(pow(EC50_AngII,1.2)+pow(AngII,1.2))*(1.0-E_AT1R*0.7) - kAldo_deg*Aldo;
dxdt_NE     = kNE_prod*NE_HF_factor - kNE_deg*(1.0+E_BB_NE*0.6)*NE - 0.015*(cGMP-8.0)*NE;
double WS = LVEDV*SVR/(LVEF*1e6+0.001);
double BNP_pr = kBNP_prod*(1.0+0.8*(WS-0.02)/0.02);
if(BNP_pr<0) BNP_pr=0;
dxdt_BNP    = BNP_pr - kBNP_NEPdeg*(1.0-E_NEP*0.9)*BNP - kBNP_other*BNP;
double NTP_EF = (1.0-(LVEF-0.27)*3.0); if(NTP_EF<0.2) NTP_EF=0.2;
dxdt_NTpBNP = 0.55*BNP_pr*NTproBNP_ratio*NTP_EF - 0.008*NTpBNP;
dxdt_cGMP   = kcGMP_prod*BNP/(BNP+200.0) - (kcGMP_PDE+kcGMP_other)*cGMP;
double LVEDV_load = 0.03*Aldo/Aldo_base;
double LVEDV_ul   = 0.025*cGMP/10.0 + E_S_V*0.008*LVEDV + 0.0010*(E_BB_NE+E_AT1R)*(LVEDV-120.0);
dxdt_LVEDV  = LVEDV_load - LVEDV_ul;
double HR_t = HR_setpoint*(1.0+0.25*NE/NE_base) - HR_setpoint*(1.0+0.25*NE/NE_base)*E_BB_HR - E_IVA*25.0;
dxdt_HR     = 0.15*(HR_t - HR);
double SVR_t = SVR_base*(1.0+0.40*(AngII-15.0)/15.0)*(1.0+0.20*(NE/NE_base-1.0))*(1.0-E_AT1R*0.35-E_S_C*0.08);
dxdt_SVR    = 0.05*(SVR_t - SVR);
double EF_w  = 0.0012*(Fib-0.1) + 0.0008*(Hyp-1.0);
double EF_b  = kEF_BB*E_BB_NE*(EF_max-LVEF) + kEF_ARNI*E_AT1R*(EF_max-LVEF) + kEF_SGLT2*E_S_C*(EF_max-LVEF);
dxdt_LVEF   = EF_b - EF_w; if(LVEF>EF_max) dxdt_LVEF=0; if(LVEF<0.05) dxdt_LVEF=0;
double TGF_p = kTGFb1_prod*(AngII/55.0)*(Aldo/22.0)*(TNFa/6.5)*(1.0-E_AT1R*0.55-E_MRA*0.40);
dxdt_TGFb1  = TGF_p - kTGFb1_deg*TGFb1;
dxdt_Fib    = kFib_prod*TGFb1*(1.0-Fib/Fib_max) - kFib_deg*Fib - E_MRA*0.008*Fib - E_AT1R*0.006*Fib;
double Hyp_p = kHyp_prod*(AngII/55.0)*(NE/715.0)*(LVEDV/280.0);
dxdt_Hyp    = Hyp_p - kHyp_deg*Hyp - kHyp_deg*E_BB_NE*Hyp*2.5; if(Hyp<1.0) dxdt_Hyp=0;
dxdt_TNFa   = kTNFa_prod*(1.0+0.5*(NE/NE_base-1.0)) - kTNFa_deg*TNFa - E_S_C*0.012*TNFa;
dxdt_IL6    = kIL6_prod*(TNFa/TNFa_base)*(AngII/55.0) - kIL6_deg*IL6;
dxdt_LBQ_C  = -(CL_LBQ/Vd_LBQ)*LBQ_C;
dxdt_Valsar_C = -(CL_Valsar/Vd_Valsar)*Valsar_C;
dxdt_BB_C   = -(CL_BB/Vd_BB)*BB_C;
dxdt_MRA_C  = -(CL_MRA/Vd_MRA)*MRA_C;
dxdt_SGLT2_C = -(CL_SGLT2/Vd_SGLT2)*SGLT2_C;
dxdt_IVA_C  = -(CL_IVA/Vd_IVA)*IVA_C;

$TABLE
double SV        = LVEDV*LVEF;
double CO        = HR*SV/1000.0;
double MAP       = SVR*CO/80.0;
double PCWP      = 8.0+(LVEDV-120.0)*0.12; if(PCWP<5)PCWP=5; if(PCWP>45)PCWP=45;
double LVEF_pct  = LVEF*100.0;
double NYHA_score = 4.0;
if(CO>=5.0 && PCWP<=12.0 && LVEF_pct>=40.0) NYHA_score=1.0;
else if(CO>=4.0 && PCWP<=18.0 && LVEF_pct>=30.0) NYHA_score=2.0;
else if(CO>=3.0 && PCWP<=25.0 && LVEF_pct>=20.0) NYHA_score=3.0;
double eGFR = 65.0*(1.0-0.25*(AngII-15.0)/55.0+0.10*(Valsar_C/(Valsar_C+EC50_AT1R_val)));
if(eGFR<10) eGFR=10; if(eGFR>90) eGFR=90;
capture SV, CO, MAP, PCWP, LVEF_pct, NYHA_score, eGFR;
'

# ─────────────────────────────────────────────────────────────────────────────
# COMPILE MODEL
# ─────────────────────────────────────────────────────────────────────────────
mod <- mread_cache("hfref_shiny", tempdir(), hfref_model_code)

# ─────────────────────────────────────────────────────────────────────────────
# DARK GGPLOT THEME
# ─────────────────────────────────────────────────────────────────────────────
theme_hf <- function() {
  theme_dark(base_size=13) +
    theme(
      plot.background  = element_rect(fill="#1a1a2e", color=NA),
      panel.background = element_rect(fill="#16213e", color=NA),
      panel.grid.major = element_line(color="#2a2a4a"),
      panel.grid.minor = element_line(color="#1e1e3a"),
      axis.text  = element_text(color="#cccccc"),
      axis.title = element_text(color="#eeeeee", face="bold"),
      plot.title = element_text(color="white", face="bold", size=14),
      legend.background = element_rect(fill="#1a1a2e"),
      legend.text = element_text(color="#cccccc"),
      strip.background = element_rect(fill="#0d3460"),
      strip.text = element_text(color="white", face="bold")
    )
}

pal_scenario <- c(
  "No Therapy"                            = "#e74c3c",
  "ACEi + BB"                             = "#e67e22",
  "ARNI + BB + MRA"                       = "#3498db",
  "ARNI + BB + MRA + SGLT2i"             = "#27ae60",
  "ARNI + BB + MRA + SGLT2i + IVA"       = "#9b59b6"
)

# ─────────────────────────────────────────────────────────────────────────────
# SIMULATION HELPER
# ─────────────────────────────────────────────────────────────────────────────
run_sim <- function(mod_local,
                    baseline_LVEF, baseline_NE_factor, baseline_Fib,
                    duration_months,
                    use_ARNI, arni_dose, use_BB, bb_dose,
                    use_MRA, mra_dose, use_SGLT2, sglt2_dose,
                    use_IVA, iva_dose,
                    use_ACEi, acei_dose) {

  duration_h <- duration_months * 30.44 * 24

  mod2 <- mod_local %>%
    init(LVEF   = baseline_LVEF / 100,
         NE     = 325 * baseline_NE_factor,
         Fib    = baseline_Fib) %>%
    param(EF_base      = baseline_LVEF / 100,
          NE_HF_factor = baseline_NE_factor,
          Fib_base     = baseline_Fib)

  ev_list <- list()

  if (use_ARNI && arni_dose > 0) {
    lbq_amt    <- arni_dose * 0.75 * 0.60
    valsar_amt <- arni_dose * 1.06 * 0.94
    ev_list$lbq    <- ev(amt=lbq_amt,    cmt="LBQ_C",    ii=12, addl=floor(duration_h/12), time=0, rate=-2)
    ev_list$valsar <- ev(amt=valsar_amt, cmt="Valsar_C",  ii=12, addl=floor(duration_h/12), time=0, rate=-2)
  }
  if (use_ACEi && acei_dose > 0 && !use_ARNI) {
    # Represent ACEi via AT1R block approximation using valsartan compartment
    acei_equiv <- acei_dose * 0.50 * 0.50   # rough enalaprilat equivalent for AT1R block
    ev_list$acei <- ev(amt=acei_equiv, cmt="Valsar_C", ii=12, addl=floor(duration_h/12), time=0, rate=-2)
  }
  if (use_BB && bb_dose > 0) {
    ev_list$bb <- ev(amt=bb_dose*0.45, cmt="BB_C", ii=24, addl=floor(duration_h/24), time=0, rate=-2)
  }
  if (use_MRA && mra_dose > 0) {
    ev_list$mra <- ev(amt=mra_dose*0.69, cmt="MRA_C", ii=24, addl=floor(duration_h/24), time=0, rate=-2)
  }
  if (use_SGLT2 && sglt2_dose > 0) {
    ev_list$sglt2 <- ev(amt=sglt2_dose*0.86, cmt="SGLT2_C", ii=24, addl=floor(duration_h/24), time=0, rate=-2)
  }
  if (use_IVA && iva_dose > 0) {
    ev_list$iva <- ev(amt=iva_dose, cmt="IVA_C", ii=12, addl=floor(duration_h/12), time=0, rate=-2)
  }

  dosing <- if (length(ev_list) > 0) Reduce(c, ev_list) else ev()

  out <- mod2 %>%
    ev(dosing) %>%
    mrgsim(end=duration_h, delta=24) %>%
    as.data.frame() %>%
    mutate(time_months = time / (30.44 * 24))

  return(out)
}

# ─────────────────────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "black",
  dashboardHeader(title = span(icon("heartbeat"), "HFrEF QSP Dashboard"),
                  titleWidth = 300),
  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",     tabName="tab_patient",  icon=icon("user-md")),
      menuItem("Drug Pharmacokinetics", tabName="tab_pk",    icon=icon("pills")),
      menuItem("RAAS & SNS",          tabName="tab_raas",    icon=icon("dna")),
      menuItem("NPS & Hemodynamics",  tabName="tab_nps",     icon=icon("heartbeat")),
      menuItem("Cardiac Remodeling",  tabName="tab_remodel", icon=icon("heart")),
      menuItem("Clinical Endpoints",  tabName="tab_endpoints",icon=icon("chart-line")),
      menuItem("Scenario Comparison", tabName="tab_compare", icon=icon("exchange-alt")),
      menuItem("Biomarker Risk",      tabName="tab_risk",    icon=icon("exclamation-triangle"))
    ),
    br(), hr(),
    h5("Patient Characteristics", style="color:#aaa; padding-left:15px"),
    sliderInput("init_LVEF",   "Baseline LVEF (%)",   min=10, max=45, value=27, step=1),
    sliderInput("init_NE_fac", "SNS Activation (×NE baseline)", min=1.0, max=4.0, value=2.2, step=0.1),
    sliderInput("init_Fib",    "Baseline Fibrosis Score (0-1)",  min=0.05, max=0.80, value=0.35, step=0.05),
    sliderInput("duration",    "Simulation Duration (months)", min=1, max=24, value=12, step=1),
    hr(),
    h5("Drug Regimen", style="color:#aaa; padding-left:15px"),
    checkboxInput("use_ARNI",  "ARNI (Sacubitril/Valsartan)", value=TRUE),
    conditionalPanel("input.use_ARNI",
      sliderInput("arni_dose", "Sacubitril dose (mg BID)", min=24, max=200, value=97, step=1)),
    checkboxInput("use_ACEi",  "ACEi (Enalapril, if no ARNI)", value=FALSE),
    conditionalPanel("input.use_ACEi && !input.use_ARNI",
      sliderInput("acei_dose", "Enalapril dose (mg BID)", min=2.5, max=20, value=10, step=2.5)),
    checkboxInput("use_BB",    "β-Blocker (Metoprolol Succ.)",  value=TRUE),
    conditionalPanel("input.use_BB",
      sliderInput("bb_dose",   "β-Blocker dose (mg/day)", min=12.5, max=200, value=200, step=12.5)),
    checkboxInput("use_MRA",   "MRA (Eplerenone)",   value=TRUE),
    conditionalPanel("input.use_MRA",
      sliderInput("mra_dose",  "MRA dose (mg/day)",  min=12.5, max=50, value=50, step=12.5)),
    checkboxInput("use_SGLT2", "SGLT2i (Empagliflozin)", value=TRUE),
    conditionalPanel("input.use_SGLT2",
      sliderInput("sglt2_dose","SGLT2i dose (mg/day)",  min=5, max=25, value=10, step=5)),
    checkboxInput("use_IVA",   "Ivabradine (if HR ≥70 bpm)", value=FALSE),
    conditionalPanel("input.use_IVA",
      sliderInput("iva_dose",  "Ivabradine dose (mg BID)", min=2.5, max=7.5, value=7.5, step=2.5)),
    br(),
    actionButton("run_sim", "Run Simulation", icon=icon("play"),
                 style="background-color:#27ae60; color:white; width:90%; margin-left:5%")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper, .right-side { background-color:#0d1117; }
      .box { background-color:#161b22; border-top:3px solid #3498db; }
      .box-title { color:#58a6ff; }
      h3, h4 { color:#c9d1d9; }
      .nav-tabs-custom .nav-tabs li.active a { background-color:#161b22; color:#58a6ff; }
    "))),
    tabItems(

      # ── TAB 1: Patient Profile ──────────────────────────────────────────
      tabItem("tab_patient",
        fluidRow(
          box(width=12, title="Patient Baseline Profile & Disease Severity",
              solidHeader=TRUE, status="primary",
              fluidRow(
                valueBoxOutput("vbox_LVEF",    width=3),
                valueBoxOutput("vbox_NYHA",    width=3),
                valueBoxOutput("vbox_BNP",     width=3),
                valueBoxOutput("vbox_NTpBNP",  width=3)
              ),
              fluidRow(
                valueBoxOutput("vbox_CO",   width=3),
                valueBoxOutput("vbox_HR",   width=3),
                valueBoxOutput("vbox_PCWP", width=3),
                valueBoxOutput("vbox_NE",   width=3)
              )
          )
        ),
        fluidRow(
          box(width=6, title="Hemodynamic State (Baseline)", solidHeader=TRUE, status="info",
              plotlyOutput("plot_baseline_hemo", height=350)),
          box(width=6, title="Neurohormonal Activation Profile", solidHeader=TRUE, status="warning",
              plotlyOutput("plot_baseline_neuro", height=350))
        ),
        fluidRow(
          box(width=12, title="Guideline-Directed Medical Therapy (GDMT) Pillar Status",
              solidHeader=TRUE, status="success",
              tableOutput("tbl_gdmt_status"))
        )
      ),

      # ── TAB 2: Pharmacokinetics ─────────────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(width=12, title="Plasma Concentration–Time Profiles (First 48h)",
              solidHeader=TRUE, status="primary",
              plotlyOutput("plot_pk_profiles", height=400))
        ),
        fluidRow(
          box(width=6, title="Drug Effect (Emax) Summary",
              solidHeader=TRUE, status="info",
              plotlyOutput("plot_effect_summary", height=350)),
          box(width=6, title="Steady-State PK Parameters",
              solidHeader=TRUE, status="warning",
              tableOutput("tbl_pk_ss"))
        )
      ),

      # ── TAB 3: RAAS & SNS ──────────────────────────────────────────────
      tabItem("tab_raas",
        fluidRow(
          box(width=6, title="Angiotensin II", solidHeader=TRUE, status="danger",
              plotlyOutput("plot_AngII", height=300)),
          box(width=6, title="Aldosterone", solidHeader=TRUE, status="danger",
              plotlyOutput("plot_Aldo", height=300))
        ),
        fluidRow(
          box(width=6, title="Plasma Norepinephrine", solidHeader=TRUE, status="warning",
              plotlyOutput("plot_NE", height=300)),
          box(width=6, title="Angiotensin 1-7 (Cardioprotective)", solidHeader=TRUE, status="success",
              plotlyOutput("plot_Ang17", height=300))
        )
      ),

      # ── TAB 4: NPS & Hemodynamics ──────────────────────────────────────
      tabItem("tab_nps",
        fluidRow(
          box(width=6, title="BNP & NT-proBNP", solidHeader=TRUE, status="primary",
              plotlyOutput("plot_BNP", height=300)),
          box(width=6, title="cGMP (Downstream of NPS)", solidHeader=TRUE, status="success",
              plotlyOutput("plot_cGMP", height=300))
        ),
        fluidRow(
          box(width=6, title="Cardiac Output & Stroke Volume", solidHeader=TRUE, status="info",
              plotlyOutput("plot_CO", height=300)),
          box(width=6, title="PCWP & SVR (Congestion + Afterload)", solidHeader=TRUE, status="warning",
              plotlyOutput("plot_PCWP_SVR", height=300))
        )
      ),

      # ── TAB 5: Cardiac Remodeling ──────────────────────────────────────
      tabItem("tab_remodel",
        fluidRow(
          box(width=6, title="LVEF Trajectory", solidHeader=TRUE, status="primary",
              plotlyOutput("plot_LVEF", height=300)),
          box(width=6, title="LVEDV (Reverse Remodeling)", solidHeader=TRUE, status="info",
              plotlyOutput("plot_LVEDV", height=300))
        ),
        fluidRow(
          box(width=6, title="Myocardial Fibrosis Score", solidHeader=TRUE, status="danger",
              plotlyOutput("plot_Fib", height=300)),
          box(width=6, title="TGF-β1 & Hypertrophy Index", solidHeader=TRUE, status="warning",
              plotlyOutput("plot_TGF_Hyp", height=300))
        )
      ),

      # ── TAB 6: Clinical Endpoints ──────────────────────────────────────
      tabItem("tab_endpoints",
        fluidRow(
          box(width=6, title="NYHA Functional Class", solidHeader=TRUE, status="primary",
              plotlyOutput("plot_NYHA", height=300)),
          box(width=6, title="Heart Rate", solidHeader=TRUE, status="info",
              plotlyOutput("plot_HR", height=300))
        ),
        fluidRow(
          box(width=6, title="Estimated GFR", solidHeader=TRUE, status="warning",
              plotlyOutput("plot_eGFR", height=300)),
          box(width=6, title="Inflammation (TNF-α, IL-6)", solidHeader=TRUE, status="danger",
              plotlyOutput("plot_inflam", height=300))
        )
      ),

      # ── TAB 7: Scenario Comparison ─────────────────────────────────────
      tabItem("tab_compare",
        fluidRow(
          box(width=12, title="Multi-Scenario Comparison — 5 Treatment Strategies",
              solidHeader=TRUE, status="primary",
              p("Comparing: No Therapy | ACEi+BB | ARNI+BB+MRA | +SGLT2i | +Ivabradine (max GDMT)"),
              plotlyOutput("plot_compare_LVEF", height=350))
        ),
        fluidRow(
          box(width=6, solidHeader=TRUE, status="info",
              title="NT-proBNP Comparison",
              plotlyOutput("plot_compare_NTpBNP", height=300)),
          box(width=6, solidHeader=TRUE, status="warning",
              title="Cardiac Output Comparison",
              plotlyOutput("plot_compare_CO", height=300))
        ),
        fluidRow(
          box(width=12, title="Summary Table at 12 Months", solidHeader=TRUE, status="success",
              DTOutput("tbl_compare_12mo"))
        )
      ),

      # ── TAB 8: Biomarker Risk ──────────────────────────────────────────
      tabItem("tab_risk",
        fluidRow(
          box(width=12, title="Biomarker-Based Risk Classification & Response Assessment",
              solidHeader=TRUE, status="danger",
              fluidRow(
                infoBoxOutput("ibox_BNP_risk",  width=4),
                infoBoxOutput("ibox_EF_class",  width=4),
                infoBoxOutput("ibox_NYHA_risk",  width=4)
              )
          )
        ),
        fluidRow(
          box(width=6, title="Risk Zone: BNP vs. LVEF (12-month trajectory)",
              solidHeader=TRUE, status="primary",
              plotlyOutput("plot_risk_scatter", height=400)),
          box(width=6, title="ESC 2021 HF Guideline Response Criteria",
              solidHeader=TRUE, status="info",
              tableOutput("tbl_esc_criteria"))
        ),
        fluidRow(
          box(width=12, title="Relative Risk Reduction Estimates (vs. No Therapy)",
              solidHeader=TRUE, status="success",
              plotlyOutput("plot_RRR", height=300))
        )
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# SERVER
# ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Reactive simulation
  sim_result <- eventReactive(input$run_sim, {
    withProgress(message="Running QSP simulation...", value=0, {
      incProgress(0.2)
      out <- run_sim(
        mod_local          = mod,
        baseline_LVEF      = input$init_LVEF,
        baseline_NE_factor = input$init_NE_fac,
        baseline_Fib       = input$init_Fib,
        duration_months    = input$duration,
        use_ARNI           = input$use_ARNI,
        arni_dose          = if(input$use_ARNI) input$arni_dose else 0,
        use_BB             = input$use_BB,
        bb_dose            = if(input$use_BB)   input$bb_dose   else 0,
        use_MRA            = input$use_MRA,
        mra_dose           = if(input$use_MRA)  input$mra_dose  else 0,
        use_SGLT2          = input$use_SGLT2,
        sglt2_dose         = if(input$use_SGLT2) input$sglt2_dose else 0,
        use_IVA            = input$use_IVA,
        iva_dose           = if(input$use_IVA)  input$iva_dose  else 0,
        use_ACEi           = input$use_ACEi,
        acei_dose          = if(input$use_ACEi) input$acei_dose else 0
      )
      incProgress(0.8)
      out
    })
  }, ignoreNULL=FALSE)

  # Multi-scenario simulation (fixed 5 scenarios)
  scenarios_result <- reactive({
    withProgress(message="Building comparison scenarios...", value=0, {
      scenarios <- list(
        list(label="No Therapy",   arni=FALSE, acei=FALSE, bb=FALSE, mra=FALSE, sglt2=FALSE, iva=FALSE),
        list(label="ACEi + BB",    arni=FALSE, acei=TRUE,  bb=TRUE,  mra=FALSE, sglt2=FALSE, iva=FALSE),
        list(label="ARNI + BB + MRA", arni=TRUE, acei=FALSE, bb=TRUE, mra=TRUE,  sglt2=FALSE, iva=FALSE),
        list(label="ARNI + BB + MRA + SGLT2i", arni=TRUE, acei=FALSE, bb=TRUE, mra=TRUE, sglt2=TRUE, iva=FALSE),
        list(label="ARNI + BB + MRA + SGLT2i + IVA", arni=TRUE, acei=FALSE, bb=TRUE, mra=TRUE, sglt2=TRUE, iva=TRUE)
      )
      do.call(rbind, lapply(seq_along(scenarios), function(i) {
        s <- scenarios[[i]]
        incProgress(0.15)
        run_sim(mod, input$init_LVEF, input$init_NE_fac, input$init_Fib,
                duration_months=12,
                use_ARNI=s$arni,  arni_dose=97,
                use_ACEi=s$acei,  acei_dose=10,
                use_BB=s$bb,      bb_dose=200,
                use_MRA=s$mra,    mra_dose=50,
                use_SGLT2=s$sglt2, sglt2_dose=10,
                use_IVA=s$iva,    iva_dose=7.5) %>%
          mutate(scenario=s$label)
      }))
    })
  })

  # ── TAB 1 VALUE BOXES ────────────────────────────────────────────────────
  init_data <- reactive({
    d <- sim_result()
    d[1, ]
  })

  output$vbox_LVEF <- renderValueBox({
    d <- init_data()
    valueBox(sprintf("%.0f%%", d$LVEF_pct), "Baseline LVEF",
             icon=icon("heart"), color=if(d$LVEF_pct<30)"red" else "yellow")
  })
  output$vbox_NYHA <- renderValueBox({
    d <- init_data()
    valueBox(sprintf("Class %.0f", d$NYHA_score), "NYHA Class",
             icon=icon("walking"), color="orange")
  })
  output$vbox_BNP <- renderValueBox({
    d <- init_data()
    valueBox(sprintf("%.0f pg/mL", d$BNP), "BNP",
             icon=icon("tachometer-alt"), color=if(d$BNP>400)"red" else "yellow")
  })
  output$vbox_NTpBNP <- renderValueBox({
    d <- init_data()
    valueBox(sprintf("%.0f pg/mL", d$NTpBNP), "NT-proBNP",
             icon=icon("chart-bar"), color=if(d$NTpBNP>3000)"red" else "yellow")
  })
  output$vbox_CO <- renderValueBox({
    d <- init_data()
    valueBox(sprintf("%.1f L/min", d$CO), "Cardiac Output",
             icon=icon("heartbeat"), color=if(d$CO<3.5)"red" else "blue")
  })
  output$vbox_HR <- renderValueBox({
    d <- init_data()
    valueBox(sprintf("%.0f bpm", d$HR), "Heart Rate",
             icon=icon("clock"), color=if(d$HR>80)"orange" else "green")
  })
  output$vbox_PCWP <- renderValueBox({
    d <- init_data()
    valueBox(sprintf("%.0f mmHg", d$PCWP), "PCWP",
             icon=icon("lungs"), color=if(d$PCWP>18)"red" else "yellow")
  })
  output$vbox_NE <- renderValueBox({
    d <- init_data()
    ne_val <- d$NE
    valueBox(sprintf("%.0f pg/mL", ne_val), "Plasma NE",
             icon=icon("bolt"), color=if(ne_val>600)"red" else "orange")
  })

  output$tbl_gdmt_status <- renderTable({
    data.frame(
      Pillar = c("RAAS Inhibition (ARNI/ACEi)", "β-Blocker", "MRA",
                 "SGLT2 Inhibitor", "Ivabradine (if HR ≥70)"),
      Status = c(
        ifelse(input$use_ARNI || input$use_ACEi, "✓ Active", "✗ Not prescribed"),
        ifelse(input$use_BB,    "✓ Active", "✗ Not prescribed"),
        ifelse(input$use_MRA,   "✓ Active", "✗ Not prescribed"),
        ifelse(input$use_SGLT2, "✓ Active", "✗ Not prescribed"),
        ifelse(input$use_IVA,   "✓ Active", "✗ Not prescribed")
      ),
      `Class/Level of Evidence` = c("Class I / Level A","Class I / Level A",
                                      "Class I / Level A","Class I / Level A",
                                      "Class IIa / Level B"),
      `Clinical Trial` = c("PARADIGM-HF / CONSENSUS","MERIT-HF / CIBIS-II",
                             "RALES / EMPHASIS-HF","EMPEROR-Reduced / DAPA-HF",
                             "SHIFT")
    )
  }, striped=TRUE, bordered=TRUE, spacing="s")

  # ── TAB 2: PK ───────────────────────────────────────────────────────────
  output$plot_pk_profiles <- renderPlotly({
    d <- sim_result() %>% filter(time_months <= 2/30.44*30)  # first 48h
    d_long <- d %>%
      select(time, LBQ_C, Valsar_C, BB_C, MRA_C, SGLT2_C, IVA_C) %>%
      pivot_longer(-time, names_to="Drug", values_to="Cp")
    p <- ggplot(d_long, aes(time/1, Cp, color=Drug)) +
      geom_line(linewidth=1.2) +
      labs(x="Time (hours)", y="Plasma Concentration (ng/mL)",
           title="Drug PK Profiles — First 48 Hours") +
      theme_hf() + scale_color_brewer(palette="Set1")
    ggplotly(p)
  })

  output$plot_effect_summary <- renderPlotly({
    d_end <- tail(sim_result(), 1)
    effects <- data.frame(
      Drug = c("NEP Inhibition\n(LBQ657)", "AT1R Block\n(Valsartan)",
               "β1-Block\n(BB)", "MR Block\n(MRA)", "SGLT2i\n(Diuresis)", "If Block\n(IVA)"),
      Effect = c(d_end$NEP_E * 100, d_end$AT1R_E * 100, d_end$BB_E * 100,
                 d_end$MRA_E * 100, d_end$SGLT2_E * 100,
                 (d_end$IVA_C / (d_end$IVA_C + 20.0)) * 100)
    )
    p <- ggplot(effects, aes(Drug, Effect, fill=Drug)) +
      geom_col(alpha=0.85) +
      geom_hline(yintercept=80, linetype="dashed", color="yellow") +
      scale_fill_brewer(palette="Set2") +
      labs(x=NULL, y="Effect (% of Emax)", title="Drug Effect at End of Simulation") +
      theme_hf() + theme(legend.position="none")
    ggplotly(p)
  })

  output$tbl_pk_ss <- renderTable({
    data.frame(
      Drug = c("LBQ657 (NEPi)", "Valsartan (ARB)", "Metoprolol (BB)",
               "Eplerenone (MRA)", "Empagliflozin (SGLT2i)", "Ivabradine"),
      `t½ (h)` = round(c(26/4.6*0.693, 75/1.9*0.693, 220/55*0.693,
                           43/7*0.693, 73.8/16.1*0.693, 100/10*0.693), 1),
      `Dose Interval` = c("BID","BID","QD","QD","QD","BID"),
      `Target Therapeutic Range` = c("LBQ657 >80 ng/mL", "Valsartan 1-6 µg/mL",
                                       "Metoprolol 50-200 ng/mL", "Eplerenone 200-800 ng/mL",
                                       "Empagliflozin 100-500 ng/mL", "Ivabradine 20-100 ng/mL")
    )
  }, striped=TRUE, bordered=TRUE, spacing="s")

  # ── PLOTTING HELPERS ─────────────────────────────────────────────────────
  make_ggplotly <- function(d, y_var, y_label, title_str, ref_line=NULL,
                             ref_label=NULL, log_y=FALSE) {
    p <- ggplot(d, aes(time_months, .data[[y_var]])) +
      geom_line(color="#3498db", linewidth=1.3) +
      labs(x="Time (months)", y=y_label, title=title_str) +
      theme_hf()
    if (!is.null(ref_line)) {
      p <- p + geom_hline(yintercept=ref_line, linetype="dashed", color="#e74c3c", alpha=0.8)
      if (!is.null(ref_label))
        p <- p + annotate("text", x=max(d$time_months)*0.75, y=ref_line*1.05,
                          label=ref_label, color="#e74c3c", size=3)
    }
    if (log_y) p <- p + scale_y_log10()
    ggplotly(p)
  }

  # ── TAB 3: RAAS & SNS ───────────────────────────────────────────────────
  output$plot_AngII  <- renderPlotly({ make_ggplotly(sim_result(), "AngII", "AngII (pg/mL)", "Angiotensin II", 55, "HFrEF baseline") })
  output$plot_Aldo   <- renderPlotly({ make_ggplotly(sim_result(), "Aldo",  "Aldosterone (ng/dL)", "Aldosterone", 22, "HFrEF baseline") })
  output$plot_NE     <- renderPlotly({ make_ggplotly(sim_result(), "NE",    "Plasma NE (pg/mL)", "Norepinephrine", 325, "Normal (<325 pg/mL)") })
  output$plot_Ang17  <- renderPlotly({ make_ggplotly(sim_result(), "Ang17", "Ang1-7 (pg/mL)", "Angiotensin 1-7 (Cardioprotective)") })

  # ── TAB 4: NPS & HEMO ───────────────────────────────────────────────────
  output$plot_BNP    <- renderPlotly({
    d <- sim_result() %>%
      select(time_months, BNP, NTpBNP) %>%
      pivot_longer(-time_months)
    p <- ggplot(d, aes(time_months, value, color=name)) +
      geom_line(linewidth=1.2) +
      geom_hline(yintercept=100,  linetype="dashed", color="#e74c3c", alpha=0.7) +
      geom_hline(yintercept=900,  linetype="dashed", color="#e67e22", alpha=0.7) +
      scale_color_manual(values=c(BNP="#3498db", NTpBNP="#e74c3c"),
                         labels=c("BNP (pg/mL)", "NT-proBNP (pg/mL)")) +
      scale_y_log10() +
      labs(x="Time (months)", y="Concentration (pg/mL, log)", title="BNP & NT-proBNP", color=NULL) +
      theme_hf()
    ggplotly(p)
  })
  output$plot_cGMP   <- renderPlotly({ make_ggplotly(sim_result(), "cGMP",  "cGMP (pmol/mL)", "Cyclic GMP") })
  output$plot_CO     <- renderPlotly({
    d <- sim_result() %>% select(time_months, CO, SV) %>% pivot_longer(-time_months)
    p <- ggplot(d, aes(time_months, value, color=name)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c(CO="#27ae60", SV="#3498db"),
                         labels=c("CO (L/min)", "SV (mL)")) +
      labs(x="Time (months)", y=NULL, title="Cardiac Output & Stroke Volume", color=NULL) +
      theme_hf()
    ggplotly(p)
  })
  output$plot_PCWP_SVR <- renderPlotly({
    d <- sim_result() %>%
      mutate(SVR_norm = SVR / 100) %>%
      select(time_months, PCWP, SVR_norm) %>% pivot_longer(-time_months)
    p <- ggplot(d, aes(time_months, value, color=name)) +
      geom_line(linewidth=1.2) +
      geom_hline(yintercept=18, linetype="dashed", color="#e74c3c", alpha=0.7) +
      scale_color_manual(values=c(PCWP="#e74c3c", SVR_norm="#f39c12"),
                         labels=c("PCWP (mmHg)", "SVR (×100 dynes·s/cm5)")) +
      labs(x="Time (months)", y=NULL, title="PCWP & SVR", color=NULL) +
      theme_hf()
    ggplotly(p)
  })

  # ── TAB 5: REMODELING ───────────────────────────────────────────────────
  output$plot_LVEF   <- renderPlotly({ make_ggplotly(sim_result(), "LVEF_pct", "LVEF (%)", "LVEF Trajectory", 40, "HFmrEF threshold (40%)") })
  output$plot_LVEDV  <- renderPlotly({ make_ggplotly(sim_result(), "LVEDV", "LVEDV (mL)", "LV End-Diastolic Volume (Reverse Remodeling)", 120, "Normal LVEDV (~120 mL)") })
  output$plot_Fib    <- renderPlotly({ make_ggplotly(sim_result(), "Fib",   "Fibrosis Score (0–1)", "Myocardial Fibrosis") })
  output$plot_TGF_Hyp <- renderPlotly({
    d <- sim_result() %>% select(time_months, TGFb1, Hyp) %>% pivot_longer(-time_months)
    p <- ggplot(d, aes(time_months, value, color=name)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c(TGFb1="#e74c3c", Hyp="#f39c12"),
                         labels=c("TGF-β1 (pg/mL)", "Hypertrophy Index")) +
      labs(x="Time (months)", y=NULL, title="TGF-β1 & Hypertrophy Index", color=NULL) +
      theme_hf()
    ggplotly(p)
  })

  # ── TAB 6: CLINICAL ENDPOINTS ───────────────────────────────────────────
  output$plot_NYHA   <- renderPlotly({ make_ggplotly(sim_result(), "NYHA_score", "NYHA Class (1–4)", "NYHA Functional Class", 2.5, "Class III boundary") })
  output$plot_HR     <- renderPlotly({ make_ggplotly(sim_result(), "HR", "Heart Rate (bpm)", "Heart Rate", 70, "Target HR <70 bpm (SHIFT)") })
  output$plot_eGFR   <- renderPlotly({ make_ggplotly(sim_result(), "eGFR", "eGFR (mL/min/1.73m²)", "Estimated GFR", 60, "CKD stage 3 threshold") })
  output$plot_inflam <- renderPlotly({
    d <- sim_result() %>% select(time_months, TNFa, IL6) %>% pivot_longer(-time_months)
    p <- ggplot(d, aes(time_months, value, color=name)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=c(TNFa="#e74c3c", IL6="#e67e22"),
                         labels=c("TNF-α (pg/mL)", "IL-6 (pg/mL)")) +
      labs(x="Time (months)", y="Cytokine (pg/mL)", title="Inflammatory Cytokines", color=NULL) +
      theme_hf()
    ggplotly(p)
  })

  # ── TAB 7: SCENARIO COMPARISON ──────────────────────────────────────────
  output$plot_compare_LVEF <- renderPlotly({
    d <- scenarios_result()
    p <- ggplot(d, aes(time_months, LVEF_pct, color=scenario)) +
      geom_line(linewidth=1.2) +
      geom_hline(yintercept=40, linetype="dashed", color="grey60") +
      scale_color_manual(values=pal_scenario) +
      labs(x="Time (months)", y="LVEF (%)", title="LVEF — All Scenarios", color="Scenario") +
      theme_hf() + theme(legend.position="right")
    ggplotly(p)
  })
  output$plot_compare_NTpBNP <- renderPlotly({
    d <- scenarios_result()
    p <- ggplot(d, aes(time_months, NTpBNP, color=scenario)) +
      geom_line(linewidth=1.2) + scale_y_log10() +
      scale_color_manual(values=pal_scenario) +
      labs(x="Time (months)", y="NT-proBNP (pg/mL, log)", title="NT-proBNP", color="Scenario") +
      theme_hf() + theme(legend.position="right")
    ggplotly(p)
  })
  output$plot_compare_CO <- renderPlotly({
    d <- scenarios_result()
    p <- ggplot(d, aes(time_months, CO, color=scenario)) +
      geom_line(linewidth=1.2) +
      scale_color_manual(values=pal_scenario) +
      labs(x="Time (months)", y="CO (L/min)", title="Cardiac Output", color="Scenario") +
      theme_hf() + theme(legend.position="right")
    ggplotly(p)
  })
  output$tbl_compare_12mo <- renderDT({
    d <- scenarios_result() %>%
      filter(abs(time_months - 12) < 0.1) %>%
      select(scenario, LVEF_pct, NTpBNP, BNP, CO, HR, PCWP, Fib, NYHA_score, eGFR) %>%
      mutate(across(where(is.numeric), ~round(.x, 1))) %>%
      rename(Scenario=scenario, `LVEF (%)`=LVEF_pct, `NT-proBNP`=NTpBNP,
             `BNP`=BNP, `CO (L/min)`=CO, `HR (bpm)`=HR,
             `PCWP (mmHg)`=PCWP, `Fibrosis`=Fib,
             `NYHA`=NYHA_score, `eGFR`=eGFR)
    datatable(d, options=list(dom='t', pageLength=5, scrollX=TRUE),
              rownames=FALSE, class="compact stripe hover")
  })

  # ── TAB 8: RISK ─────────────────────────────────────────────────────────
  end_data <- reactive({ tail(sim_result(), 1) })

  output$ibox_BNP_risk <- renderInfoBox({
    d <- end_data()
    bnp_val <- d$NTpBNP
    status <- if(bnp_val > 2000) "High Risk" else if(bnp_val > 1000) "Intermediate" else "Low Risk"
    color  <- if(bnp_val > 2000) "red" else if(bnp_val > 1000) "yellow" else "green"
    infoBox("NT-proBNP Risk", sprintf("%.0f pg/mL\n%s", bnp_val, status),
            icon=icon("chart-bar"), color=color)
  })
  output$ibox_EF_class <- renderInfoBox({
    d <- end_data()
    ef <- d$LVEF_pct
    cls <- if(ef >= 50) "HFpEF (≥50%)" else if(ef >= 40) "HFmrEF (40-49%)" else "HFrEF (<40%)"
    col <- if(ef >= 50) "green" else if(ef >= 40) "yellow" else "red"
    infoBox("EF Classification", sprintf("%.0f%%\n%s", ef, cls),
            icon=icon("heart"), color=col)
  })
  output$ibox_NYHA_risk <- renderInfoBox({
    d <- end_data()
    nyha <- round(d$NYHA_score)
    col <- c("green","yellow","orange","red")[pmax(1,pmin(4,nyha))]
    infoBox("NYHA Class", sprintf("Class %d", nyha),
            icon=icon("walking"), color=col)
  })

  output$plot_risk_scatter <- renderPlotly({
    d <- sim_result() %>%
      filter(row_number() %% 5 == 0) %>%
      mutate(risk_zone = case_when(
        LVEF_pct < 30 & NTpBNP > 2000 ~ "Very High Risk",
        LVEF_pct < 35 | NTpBNP > 1000 ~ "High Risk",
        LVEF_pct < 40 | NTpBNP > 500  ~ "Intermediate",
        TRUE ~ "Lower Risk"))
    p <- ggplot(d, aes(LVEF_pct, NTpBNP, color=risk_zone,
                       text=sprintf("Month: %.1f\nLVEF: %.1f%%\nNT-proBNP: %.0f", time_months, LVEF_pct, NTpBNP))) +
      geom_path(color="grey50", linewidth=0.5, alpha=0.5) +
      geom_point(size=2.5) +
      scale_y_log10() +
      geom_vline(xintercept=c(40,50), linetype="dashed", color="white", alpha=0.4) +
      geom_hline(yintercept=c(900,3000), linetype="dashed", color="white", alpha=0.4) +
      scale_color_manual(values=c("Very High Risk"="#e74c3c","High Risk"="#e67e22",
                                   "Intermediate"="#f1c40f","Lower Risk"="#27ae60")) +
      labs(x="LVEF (%)", y="NT-proBNP (log)", title="Risk Zone Trajectory", color="Risk Zone") +
      theme_hf()
    ggplotly(p, tooltip="text")
  })

  output$tbl_esc_criteria <- renderTable({
    d <- end_data()
    data.frame(
      Criterion = c("LVEF ≥50%","LVEF 40-49%","NT-proBNP <900 pg/mL",
                    "BNP <100 pg/mL","HR <70 bpm (if sinus)","NYHA Class I-II",
                    "LVEDV normalization","CO ≥5 L/min"),
      Target = c("≥50%","40-49% (HFmrEF)","<900","<100","<70","I or II",
                 "<150 mL","≥5.0"),
      Current = c(sprintf("%.1f%%", d$LVEF_pct),
                  ifelse(d$LVEF_pct>=40 & d$LVEF_pct<50,"Yes","No"),
                  sprintf("%.0f", d$NTpBNP),
                  sprintf("%.0f", d$BNP),
                  sprintf("%.0f", d$HR),
                  ifelse(d$NYHA_score<=2,"Class I-II","Class III-IV"),
                  sprintf("%.0f mL", d$LVEDV),
                  sprintf("%.1f", d$CO)),
      Met = c(ifelse(d$LVEF_pct>=50,"✓","✗"),
              ifelse(d$LVEF_pct>=40 & d$LVEF_pct<50,"✓","–"),
              ifelse(d$NTpBNP<900,"✓","✗"),
              ifelse(d$BNP<100,"✓","✗"),
              ifelse(d$HR<70,"✓","✗"),
              ifelse(d$NYHA_score<=2,"✓","✗"),
              ifelse(d$LVEDV<150,"✓","✗"),
              ifelse(d$CO>=5,"✓","✗"))
    )
  }, striped=TRUE, bordered=TRUE, spacing="s")

  output$plot_RRR <- renderPlotly({
    d_cmp <- scenarios_result() %>%
      filter(abs(time_months - 12) < 0.1)
    baseline_LVEF   <- d_cmp$LVEF_pct[d_cmp$scenario == "No Therapy"]
    baseline_NTpBNP <- d_cmp$NTpBNP[d_cmp$scenario == "No Therapy"]
    baseline_CO     <- d_cmp$CO[d_cmp$scenario == "No Therapy"]

    d_rrr <- d_cmp %>%
      filter(scenario != "No Therapy") %>%
      mutate(
        LVEF_gain    = LVEF_pct - baseline_LVEF,
        NTpBNP_red   = (1 - NTpBNP / baseline_NTpBNP) * 100,
        CO_gain      = CO - baseline_CO
      ) %>%
      select(scenario, LVEF_gain, NTpBNP_red, CO_gain) %>%
      pivot_longer(-scenario, names_to="Metric", values_to="Value")

    p <- ggplot(d_rrr, aes(scenario, Value, fill=scenario)) +
      geom_col(alpha=0.85, show.legend=FALSE) +
      facet_wrap(~Metric, scales="free_y",
                 labeller=labeller(Metric=c(LVEF_gain="LVEF Gain (%pt)",
                                            NTpBNP_red="NT-proBNP Reduction (%)",
                                            CO_gain="CO Increase (L/min)"))) +
      scale_fill_manual(values=pal_scenario) +
      labs(x=NULL, y="Improvement vs. No Therapy",
           title="Treatment Benefit at 12 Months (vs. No Therapy)") +
      theme_hf() + theme(axis.text.x=element_text(angle=30, hjust=1, size=8))
    ggplotly(p)
  })
}

# ─────────────────────────────────────────────────────────────────────────────
# LAUNCH
# ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui=ui, server=server)
