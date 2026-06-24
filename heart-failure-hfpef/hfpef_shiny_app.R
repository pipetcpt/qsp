##############################################################
#  HFpEF QSP Interactive Shiny Dashboard
#  Heart Failure with Preserved Ejection Fraction
#  6+ Tabs: Patient Profile | PK | PD Biomarkers |
#           Clinical Endpoints | Scenario Comparison | Biomarkers
#  Author: Claude Code Routine (CCR)
#  Date: 2026-06-17
##############################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)

# ─────────────────────────────────────────────────────────────
# MODEL CODE (identical to standalone model, embedded for Shiny)
# ─────────────────────────────────────────────────────────────
hfpef_code <- '
$PROB HFpEF QSP Model — Shiny Interface

$PARAM
Age=72 Weight=95 BSA=2.1 BMI=33
base_HTN=1 base_T2DM=1 base_AF=0 base_CKD_sev=0.3
empa_F=0.78 empa_ka=1.5 empa_CL=9.4 empa_V1=73.8 empa_Q=18.0 empa_V2=36.0
empa_EC50=5.0 empa_Emax=1.0
sac_F_sac=0.60 sac_ka=1.0 sac_CL_sac=4.5 sac_V_sac=35.0
sac_kconv=0.25 val_F=0.23 val_CL=1.3 val_V=17.0
val_EC50=0.10 sac_EC50_NEP=5.0
fin_F=0.43 fin_ka=1.2 fin_CL=3.2 fin_V=52.0 fin_EC50=0.9 fin_Emax=0.95
furo_F=0.60 furo_ka=2.0 furo_CL=8.0 furo_V=12.0 furo_Vt=0.35
furo_EC50=2.0 furo_Emax=0.80
base_AngII=1.0 base_Aldo=1.0 kout_AngII=0.5 kin_AngII=0.5 AngII_Aldo_EC50=1.5
base_ANP=1.0 base_BNP=1.0 kout_ANP=0.35 kout_BNP=0.20 kin_ANP=0.35 kin_BNP=0.20
LVEDP_ANP_EC50=20 BNP_hill=2.0
base_LVM_idx=110 base_fibrosis=0.30 base_Ecol=1.3
kgrowth_LVM=0.001 kdecay_LVM=0.002 kfib_up=0.005 kfib_down=0.003
base_titin_stiff=0.45 ktitin_PKG=0.10 ktitin_age=0.001
base_HR=78 base_SV=70 base_SVR=1400 base_LVEDP=16 base_PCWP=18
LVEDP_sv_sens=-0.3
base_GFR=62 base_Na_exc=1.0 kGFR_aldo=0.05 base_UricAcid=7.2
base_CRP=3.5 base_IL6=1.0 base_TNFa=1.0 base_sST2=1.0
kout_CRP=0.02 kin_CRP=0.07
base_cGMP=1.0 base_PKG=1.0 kout_cGMP=0.3 kin_cGMP=0.3 PKG_ox_EC50=2.0
use_empa=0 use_arni=0 use_fin=0 use_furo=0

$CMT
EMPA_GUT EMPA_C EMPA_P
SAC_GUT SAC_C VAL_GUT VAL_C
FIN_GUT FIN_C
FURO_GUT FURO_C
ANGII ALDO ANP BNP
LVM_IDX FIBROSIS TITIN_STIFF
LVEDP_dyn SVR_dyn
CRP_dyn IL6_dyn
cGMP_dyn PKG_dyn
GFR_dyn Na_exc_dyn
NT_proBNP_dyn

$MAIN
F_EMPA_GUT = empa_F;
F_SAC_GUT  = sac_F_sac;
F_VAL_GUT  = val_F;
F_FIN_GUT  = fin_F;
F_FURO_GUT = furo_F;

$ODE
double empa_Cp = EMPA_C / empa_V1;
dxdt_EMPA_GUT = -empa_ka * EMPA_GUT;
dxdt_EMPA_C   =  empa_ka * EMPA_GUT - (empa_CL/empa_V1)*EMPA_C - (empa_Q/empa_V1)*EMPA_C + (empa_Q/empa_V2)*EMPA_P;
dxdt_EMPA_P   =  (empa_Q/empa_V1)*EMPA_C - (empa_Q/empa_V2)*EMPA_P;
double SGLT2_inh = use_empa * empa_Emax * empa_Cp / (empa_EC50 + empa_Cp);

double sac_Cp = SAC_C / sac_V_sac;
double val_Cp = VAL_C / val_V;
dxdt_SAC_GUT = -sac_ka * SAC_GUT;
dxdt_SAC_C   =  sac_ka * SAC_GUT * sac_kconv - (sac_CL_sac/sac_V_sac)*SAC_C;
dxdt_VAL_GUT = -sac_ka * VAL_GUT;
dxdt_VAL_C   =  sac_ka * VAL_GUT - (val_CL/val_V)*VAL_C;
double NEP_inh = use_arni * sac_Cp / (sac_EC50_NEP + sac_Cp);
double AT1_blk = use_arni * val_Cp  / (val_EC50     + val_Cp);

double fin_Cp = FIN_C / fin_V;
dxdt_FIN_GUT = -fin_ka * FIN_GUT;
dxdt_FIN_C   =  fin_ka * FIN_GUT - (fin_CL/fin_V)*FIN_C;
double MR_blk = use_fin * fin_Emax * fin_Cp / (fin_EC50 + fin_Cp);

double furo_Cp_tubular = (FURO_C/furo_V)*furo_Vt;
dxdt_FURO_GUT = -furo_ka * FURO_GUT;
dxdt_FURO_C   =  furo_ka * FURO_GUT - (furo_CL/furo_V)*FURO_C;
double furo_nat = use_furo * furo_Emax * furo_Cp_tubular / (furo_EC50 + furo_Cp_tubular);

double AngII_current = base_AngII + ANGII;
double AT1_eff = 1.0 - AT1_blk;
double kin_ANGII_adj = kin_AngII * AT1_eff * (1.0 + 0.3*base_HTN);
dxdt_ANGII = kin_ANGII_adj - kout_AngII * AngII_current;

double Aldo_stim = AngII_current*AngII_current/(AngII_Aldo_EC50*AngII_Aldo_EC50 + AngII_current*AngII_current);
double Aldo_current = base_Aldo + ALDO;
double MR_suppress = 1.0 - 0.5*MR_blk;
dxdt_ALDO = kin_AngII * Aldo_stim * MR_suppress - kout_AngII * Aldo_current;

double LVEDP_current = base_LVEDP + LVEDP_dyn;
double ANP_current = base_ANP + ANP;
double BNP_current = base_BNP + BNP;
double LVEDP_stim = pow(LVEDP_current, BNP_hill)/(pow(LVEDP_ANP_EC50, BNP_hill) + pow(LVEDP_current, BNP_hill));
double kout_ANP_eff = kout_ANP * (1.0 - 0.7*NEP_inh);
double kout_BNP_eff = kout_BNP * (1.0 - 0.7*NEP_inh);
dxdt_ANP = kin_ANP*(1.0 + 2.0*LVEDP_stim) - kout_ANP_eff*ANP_current;
dxdt_BNP = kin_BNP*(1.0 + 3.0*LVEDP_stim) - kout_BNP_eff*BNP_current;

double cGMP_current = base_cGMP + cGMP_dyn;
double PKG_current  = base_PKG  + PKG_dyn;
double cGMP_synth = kin_cGMP*(1.0 + 0.8*(ANP_current-1.0) + 0.5*(BNP_current-1.0));
double PDE_activity = 1.0 + 0.3*base_CKD_sev;
double cGMP_degrad  = kout_cGMP*PDE_activity*cGMP_current;
dxdt_cGMP_dyn = cGMP_synth - cGMP_degrad - base_cGMP*kout_cGMP;

double ROS_proxy = 1.0 + 0.4*base_T2DM + 0.3*base_HTN + 0.2*base_CKD_sev;
double PKG_inhibit = ROS_proxy/(PKG_ox_EC50 + ROS_proxy);
double PKG_activation = cGMP_current/(1.0 + cGMP_current);
dxdt_PKG_dyn = 0.5*PKG_activation*(1.0-PKG_inhibit) - 0.3*PKG_current + 0.3*base_PKG;

double LVM_current  = base_LVM_idx  + LVM_IDX;
double fibr_current = base_fibrosis + FIBROSIS;
double titin_current= base_titin_stiff + TITIN_STIFF;
double LVM_AngII_drive = kgrowth_LVM * AngII_current * (1.0 - 0.4*AT1_blk);
double LVM_regress  = kdecay_LVM * (LVM_current - 95.0);
dxdt_LVM_IDX = LVM_AngII_drive - LVM_regress;

double TGFb_proxy = 1.0 + 0.3*Aldo_current + 0.2*AngII_current;
double fib_drive  = kfib_up  * TGFb_proxy * (1.0-0.5*MR_blk) * (1.0-0.3*SGLT2_inh);
double fib_regress= kfib_down * fibr_current;
dxdt_FIBROSIS = fib_drive - fib_regress;

double titin_age_effect = ktitin_age*(Age-60.0)*0.01;
double titin_PKG_soften = ktitin_PKG*PKG_current;
dxdt_TITIN_STIFF = titin_age_effect - titin_PKG_soften*titin_current + ktitin_PKG*base_titin_stiff*0.1;

double titin_LVEDP = 8.0*titin_current;
double fibr_LVEDP  = 6.0*(fibr_current/0.30);
double vol_LVEDP   = 4.0*(Aldo_current-1.0)*(1.0-furo_nat*0.8)*(1.0-SGLT2_inh*0.4);
double target_LVEDP= titin_LVEDP + fibr_LVEDP + vol_LVEDP + 4.0;
dxdt_LVEDP_dyn = 0.05*(target_LVEDP - LVEDP_current);

double SVR_current = base_SVR + SVR_dyn;
double target_SVR  = base_SVR*(1.0+0.15*(AngII_current-1.0))*(1.0-0.20*AT1_blk)*(1.0-0.10*(cGMP_current-1.0));
dxdt_SVR_dyn = 0.03*(target_SVR - SVR_current);

double CRP_current = base_CRP + CRP_dyn;
double IL6_current = base_IL6 + IL6_dyn;
double IL6_drive = 0.1*(1.0+0.3*BMI/30.0+0.2*base_T2DM+0.2*AngII_current)*(1.0-0.25*SGLT2_inh)*(1.0-0.15*MR_blk);
double IL6_decay = 0.05*IL6_current;
dxdt_IL6_dyn = IL6_drive - IL6_decay;
double CRP_synth = kin_CRP*IL6_current;
double CRP_decay = kout_CRP*CRP_current;
dxdt_CRP_dyn = CRP_synth - CRP_decay;

double GFR_current = base_GFR + GFR_dyn;
double Na_current  = base_Na_exc + Na_exc_dyn;
double GFR_aldo_effect = -kGFR_aldo*(Aldo_current-1.0);
double GFR_sglt2_init  = 0.05*SGLT2_inh*(-1.0);
double GFR_target  = base_GFR*(1.0-0.2*base_CKD_sev) + GFR_aldo_effect + GFR_sglt2_init*0.5;
dxdt_GFR_dyn = 0.02*(GFR_target - GFR_current);

double Na_exc_target = base_Na_exc + furo_nat*0.5 + SGLT2_inh*0.15 + NEP_inh*0.10 - (Aldo_current-1.0)*0.20;
dxdt_Na_exc_dyn = 0.1*(Na_exc_target - Na_current);

double NTpBNP_current = 400.0 + NT_proBNP_dyn;
double NTpBNP_target  = 400.0 + 15.0*(LVEDP_current-base_LVEDP) + 25.0*(fibr_current-base_fibrosis)*1000.0
                        - 80.0*furo_nat - 60.0*SGLT2_inh;
dxdt_NT_proBNP_dyn = 0.008*(NTpBNP_target - NTpBNP_current);

$TABLE
capture Empa_Cp   = EMPA_C/empa_V1;
capture Sac_Cp    = SAC_C/sac_V_sac;
capture Val_Cp    = VAL_C/val_V;
capture Fin_Cp    = FIN_C/fin_V;
capture Furo_Cp   = FURO_C/furo_V;
capture SGLT2_inh_out = use_empa * empa_Emax * Empa_Cp/(empa_EC50+Empa_Cp);
capture NEP_inh_out   = use_arni * Sac_Cp/(sac_EC50_NEP+Sac_Cp);
capture AT1_blk_out   = use_arni * Val_Cp/(val_EC50+Val_Cp);
capture MR_blk_out    = use_fin  * fin_Emax * Fin_Cp/(fin_EC50+Fin_Cp);
capture AngII_n   = base_AngII + ANGII;
capture Aldo_n    = base_Aldo  + ALDO;
capture ANP_lv    = base_ANP   + ANP;
capture BNP_lv    = base_BNP   + BNP;
capture LVM       = base_LVM_idx   + LVM_IDX;
capture Fibr      = base_fibrosis  + FIBROSIS;
capture Titin     = base_titin_stiff + TITIN_STIFF;
capture cGMP_lv   = base_cGMP + cGMP_dyn;
capture PKG_lv    = base_PKG  + PKG_dyn;
capture LVEDP_out = base_LVEDP + LVEDP_dyn;
capture SVR_out   = base_SVR   + SVR_dyn;
capture CO_out    = (base_SV/1000.0)*base_HR;
capture MAP_out   = CO_out * SVR_out / 80.0;
capture eGFR_out  = base_GFR + GFR_dyn;
capture NaExc     = base_Na_exc + Na_exc_dyn;
capture hsCRP_out = base_CRP + CRP_dyn;
capture IL6_out   = base_IL6 + IL6_dyn;
capture NTpBNP    = 400.0 + NT_proBNP_dyn;
capture RiskScore = 20.0*(LVEDP_out/20.0) + 20.0*(NTpBNP/400.0) + 15.0*Fibr + 15.0*(1.0-eGFR_out/60.0) + 10.0*(hsCRP_out/3.0) + 10.0*Titin + 10.0*(LVM/110.0);

$CAPTURE
Empa_Cp Sac_Cp Val_Cp Fin_Cp Furo_Cp
SGLT2_inh_out NEP_inh_out AT1_blk_out MR_blk_out
AngII_n Aldo_n ANP_lv BNP_lv
LVM Fibr Titin cGMP_lv PKG_lv
LVEDP_out SVR_out CO_out MAP_out
eGFR_out NaExc hsCRP_out IL6_out NTpBNP RiskScore
'

# Compile model once globally
MOD <- mread_cache("hfpef_shiny", temp = TRUE, code = hfpef_code, quiet = TRUE)

# ─────────────────────────────────────────────────────────────
# HELPER: run single simulation
# ─────────────────────────────────────────────────────────────
run_sim <- function(mod, params, evs, n_weeks = 52) {
  mod %>%
    param(params) %>%
    ev(evs) %>%
    mrgsim(end = n_weeks * 7 * 24, delta = 24, obsonly = TRUE) %>%
    as_tibble() %>%
    mutate(day = time / 24, week = day / 7)
}

build_events <- function(use_empa, empa_dose,
                         use_arni, arni_dose,
                         use_fin,  fin_dose,
                         use_furo, furo_dose,
                         n_weeks = 52) {
  n_days <- n_weeks * 7
  ev_list <- ev()
  if (use_empa > 0 && empa_dose > 0)
    ev_list <- ev_list + ev(amt = empa_dose * 1000, cmt = "EMPA_GUT",
                            ii = 24, addl = n_days - 1)
  if (use_arni > 0 && arni_dose > 0) {
    ev_list <- ev_list + ev(amt = arni_dose * 1000, cmt = "SAC_GUT",
                            ii = 12, addl = n_days * 2 - 1)
    ev_list <- ev_list + ev(amt = 103 * 1000, cmt = "VAL_GUT",
                            ii = 12, addl = n_days * 2 - 1)
  }
  if (use_fin > 0 && fin_dose > 0)
    ev_list <- ev_list + ev(amt = fin_dose * 1000, cmt = "FIN_GUT",
                            ii = 24, addl = n_days - 1)
  if (use_furo > 0 && furo_dose > 0)
    ev_list <- ev_list + ev(amt = furo_dose * 1000, cmt = "FURO_GUT",
                            ii = 12, addl = n_days * 2 - 1)
  ev_list
}

theme_hfpef <- theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title    = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#1a3a5c"),
    strip.text = element_text(color = "white", face = "bold")
  )

# ─────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "HFpEF QSP Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "tab_patient",   icon = icon("user")),
      menuItem("Drug PK",              tabName = "tab_pk",        icon = icon("pills")),
      menuItem("PD Biomarkers",        tabName = "tab_pd",        icon = icon("heartbeat")),
      menuItem("Clinical Endpoints",   tabName = "tab_clinical",  icon = icon("hospital")),
      menuItem("Scenario Comparison",  tabName = "tab_compare",   icon = icon("chart-bar")),
      menuItem("Biomarker Panel",      tabName = "tab_biomarker", icon = icon("vials"))
    ),
    hr(),
    h5("Simulation Settings", style = "padding-left:15px; color:#aaa"),
    sliderInput("n_weeks", "Duration (weeks):", min = 4, max = 104, value = 52, step = 4)
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box-header { background: #1a3a5c !important; color: white !important; }
      .box-header .box-title { color: white !important; }
    "))),

    tabItems(

      # ── TAB 1: Patient Profile ──────────────────────────────
      tabItem("tab_patient",
        fluidRow(
          box(title = "Demographics", width = 4, solidHeader = TRUE, status = "primary",
            sliderInput("pt_age",    "Age (years):",   min=50, max=90, value=72),
            sliderInput("pt_bmi",    "BMI (kg/m²):",   min=20, max=55, value=33),
            sliderInput("pt_weight", "Weight (kg):",   min=50, max=160, value=95),
            selectInput("pt_sex",    "Sex:", c("Female","Male"), selected="Female")
          ),
          box(title = "Comorbidities", width = 4, solidHeader = TRUE, status = "warning",
            checkboxInput("pt_htn",  "Hypertension (HTN)",              value = TRUE),
            checkboxInput("pt_t2dm", "Type 2 Diabetes (T2DM)",          value = TRUE),
            checkboxInput("pt_af",   "Atrial Fibrillation (AF)",         value = FALSE),
            sliderInput("pt_ckd",    "CKD Severity (0=none, 1=severe):",
                        min=0, max=1, value=0.3, step=0.1)
          ),
          box(title = "Baseline Cardiac Status", width = 4, solidHeader = TRUE, status = "danger",
            sliderInput("pt_lvedp",  "Baseline LVEDP (mmHg):", min=10, max=30, value=16),
            sliderInput("pt_gfr",    "Baseline eGFR (mL/min/1.73m²):", min=20, max=100, value=62),
            sliderInput("pt_fibr",   "LV Fibrosis (fraction):", min=0.1, max=0.6, value=0.3, step=0.05),
            sliderInput("pt_lvm",    "LV Mass Index (g/m²):", min=80, max=160, value=110)
          )
        ),
        fluidRow(
          box(title = "Treatment Selection", width = 6, solidHeader = TRUE, status = "success",
            checkboxInput("tx_empa", "Empagliflozin (SGLT2i)", value = TRUE),
            conditionalPanel("input.tx_empa",
              sliderInput("dose_empa", "Empagliflozin dose (mg QD):",
                          min=2.5, max=25, value=10, step=2.5)
            ),
            checkboxInput("tx_arni", "Sacubitril/Valsartan (ARNI)", value = FALSE),
            conditionalPanel("input.tx_arni",
              selectInput("dose_arni", "ARNI dose (sacubitril mg BID):",
                          c("24/26" = 24, "49/51" = 49, "97/103" = 97), selected = "97")
            ),
            checkboxInput("tx_fin",  "Finerenone (MRA)", value = FALSE),
            conditionalPanel("input.tx_fin",
              sliderInput("dose_fin", "Finerenone dose (mg QD):", min=10, max=40, value=20, step=10)
            ),
            checkboxInput("tx_furo", "Furosemide (Loop Diuretic)", value = TRUE),
            conditionalPanel("input.tx_furo",
              sliderInput("dose_furo", "Furosemide dose (mg BID):", min=20, max=160, value=40, step=20)
            ),
            actionButton("run_sim", "Run Simulation", class="btn btn-primary btn-lg",
                         icon = icon("play"))
          ),
          box(title = "Patient Risk Summary", width = 6, solidHeader = TRUE, status = "info",
            valueBoxOutput("vb_nyha",   width = 6),
            valueBoxOutput("vb_ntpbnp", width = 6),
            valueBoxOutput("vb_lvedp",  width = 6),
            valueBoxOutput("vb_egfr",   width = 6),
            br(),
            helpText("Risk classification is updated after simulation run.")
          )
        )
      ),

      # ── TAB 2: Drug PK ─────────────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(title = "Plasma Concentration — Time Profile", width = 12,
              solidHeader = TRUE, status = "primary",
            plotOutput("pk_plot", height = "420px"),
            helpText("Steady-state profiles shown after initial dosing period.")
          )
        ),
        fluidRow(
          box(title = "PK Parameters Summary", width = 6, solidHeader = TRUE, status = "info",
            tableOutput("pk_table")
          ),
          box(title = "Drug-Target Engagement (%)", width = 6, solidHeader = TRUE, status = "warning",
            plotOutput("te_plot", height = "280px")
          )
        )
      ),

      # ── TAB 3: PD Biomarkers ───────────────────────────────
      tabItem("tab_pd",
        fluidRow(
          box(title = "Neurohumoral Markers over Time", width = 6,
              solidHeader = TRUE, status = "primary",
            plotOutput("neuro_plot", height = "320px")
          ),
          box(title = "cGMP-PKG Pathway Activation", width = 6,
              solidHeader = TRUE, status = "success",
            plotOutput("cgmp_plot", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Cardiac Structure — LVM & Fibrosis", width = 6,
              solidHeader = TRUE, status = "warning",
            plotOutput("struct_plot", height = "320px")
          ),
          box(title = "Inflammation Markers", width = 6,
              solidHeader = TRUE, status = "danger",
            plotOutput("inflam_plot", height = "320px")
          )
        )
      ),

      # ── TAB 4: Clinical Endpoints ──────────────────────────
      tabItem("tab_clinical",
        fluidRow(
          box(title = "LV End-Diastolic Pressure (LVEDP)", width = 6,
              solidHeader = TRUE, status = "primary",
            plotOutput("lvedp_plot", height = "320px"),
            helpText("Normal LVEDP < 12 mmHg. HFpEF typically 15-25 mmHg at rest.")
          ),
          box(title = "NT-proBNP Trajectory", width = 6,
              solidHeader = TRUE, status = "danger",
            plotOutput("ntpbnp_plot", height = "320px"),
            helpText("Diagnostic threshold ≥125 pg/mL (HFpEF 2021 ESC Guidelines).")
          )
        ),
        fluidRow(
          box(title = "Renal Function — eGFR", width = 6,
              solidHeader = TRUE, status = "warning",
            plotOutput("egfr_plot", height = "320px"),
            helpText("SGLT2 inhibitors cause an initial hemodynamic eGFR dip (~3-5 mL/min) that is reversible and protective long-term.")
          ),
          box(title = "HF Hospitalization Risk Score", width = 6,
              solidHeader = TRUE, status = "info",
            plotOutput("risk_plot", height = "320px"),
            helpText("Composite risk integrating LVEDP, NT-proBNP, fibrosis, eGFR, CRP, Titin, LVM.")
          )
        )
      ),

      # ── TAB 5: Scenario Comparison ─────────────────────────
      tabItem("tab_compare",
        fluidRow(
          box(title = "Multi-Scenario Comparison at 52 Weeks",
              width = 12, solidHeader = TRUE, status = "primary",
            p("Compare 5 treatment scenarios: Placebo, Empagliflozin, ARNI, Finerenone, Combination."),
            actionButton("run_compare", "Run All Scenarios", class="btn btn-success btn-lg",
                         icon = icon("layer-group")),
            hr(),
            plotOutput("compare_plot", height = "500px")
          )
        ),
        fluidRow(
          box(title = "52-Week Outcomes Table", width = 12,
              solidHeader = TRUE, status = "info",
            DTOutput("compare_table")
          )
        )
      ),

      # ── TAB 6: Biomarker Panel ─────────────────────────────
      tabItem("tab_biomarker",
        fluidRow(
          box(title = "Biomarker Dashboard — Radar / Spider Chart",
              width = 8, solidHeader = TRUE, status = "primary",
            plotOutput("radar_plot", height = "450px")
          ),
          box(title = "Biomarker Reference Ranges", width = 4,
              solidHeader = TRUE, status = "info",
            tableOutput("biomarker_ref"),
            hr(),
            h5("Week-52 Values"),
            tableOutput("biomarker_w52")
          )
        ),
        fluidRow(
          box(title = "Uric Acid & Renal Biomarkers", width = 6,
              solidHeader = TRUE, status = "warning",
            plotOutput("uricacid_plot", height = "280px")
          ),
          box(title = "Dose-Response: Empagliflozin", width = 6,
              solidHeader = TRUE, status = "danger",
            sliderInput("dr_weeks", "Evaluate at week:", min=4, max=104, value=52, step=4),
            actionButton("run_dr", "Run Dose-Response", class="btn btn-warning",
                         icon = icon("chart-line")),
            plotOutput("dr_plot", height = "240px")
          )
        )
      )
    )
  )
)

# ─────────────────────────────────────────────────────────────
# SERVER
# ─────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Reactive: simulation results for single patient ────────
  sim_data <- reactiveVal(NULL)

  observeEvent(input$run_sim, {
    withProgress(message = "Running HFpEF simulation...", value = 0.5, {
      pars <- list(
        Age           = input$pt_age,
        BMI           = input$pt_bmi,
        Weight        = input$pt_weight,
        base_HTN      = as.integer(input$pt_htn),
        base_T2DM     = as.integer(input$pt_t2dm),
        base_AF       = as.integer(input$pt_af),
        base_CKD_sev  = input$pt_ckd,
        base_LVEDP    = input$pt_lvedp,
        base_GFR      = input$pt_gfr,
        base_fibrosis = input$pt_fibr,
        base_LVM_idx  = input$pt_lvm,
        use_empa      = as.integer(input$tx_empa),
        use_arni      = as.integer(input$tx_arni),
        use_fin       = as.integer(input$tx_fin),
        use_furo      = as.integer(input$tx_furo)
      )

      evs <- build_events(
        use_empa  = pars$use_empa,  empa_dose = if(input$tx_empa) input$dose_empa else 0,
        use_arni  = pars$use_arni,  arni_dose = if(input$tx_arni) as.numeric(input$dose_arni) else 0,
        use_fin   = pars$use_fin,   fin_dose  = if(input$tx_fin)  input$dose_fin  else 0,
        use_furo  = pars$use_furo,  furo_dose = if(input$tx_furo) input$dose_furo else 0,
        n_weeks   = input$n_weeks
      )

      out <- run_sim(MOD, pars, evs, n_weeks = input$n_weeks)
      sim_data(out)
      setProgress(1)
    })
  })

  # Auto-run on startup with defaults
  observe({
    if (is.null(sim_data())) {
      isolate({
        pars <- list(Age=72, BMI=33, Weight=95, base_HTN=1, base_T2DM=1,
                     base_AF=0, base_CKD_sev=0.3, base_LVEDP=16, base_GFR=62,
                     base_fibrosis=0.30, base_LVM_idx=110,
                     use_empa=1, use_arni=0, use_fin=0, use_furo=1)
        evs <- build_events(1, 10, 0, 0, 0, 0, 1, 40, n_weeks=52)
        out <- run_sim(MOD, pars, evs, n_weeks=52)
        sim_data(out)
      })
    }
  })

  # ── Value boxes ────────────────────────────────────────────
  get_w52 <- function() {
    req(sim_data())
    sim_data() %>% filter(week >= max(week) - 0.5) %>%
      summarise_if(is.numeric, mean, na.rm = TRUE)
  }

  output$vb_nyha <- renderValueBox({
    d <- get_w52()
    lvedp <- d$LVEDP_out[1]
    nyha <- ifelse(lvedp < 14, "I", ifelse(lvedp < 18, "II",
            ifelse(lvedp < 22, "III", "IV")))
    clr  <- ifelse(lvedp < 14, "green", ifelse(lvedp < 18, "yellow", "red"))
    valueBox(paste("NYHA", nyha), "Functional Class", icon=icon("heartbeat"), color=clr)
  })

  output$vb_ntpbnp <- renderValueBox({
    d <- get_w52()
    val <- round(d$NTpBNP[1])
    clr <- ifelse(val < 125, "green", ifelse(val < 400, "yellow", "red"))
    valueBox(paste(val, "pg/mL"), "NT-proBNP", icon=icon("vials"), color=clr)
  })

  output$vb_lvedp <- renderValueBox({
    d <- get_w52()
    val <- round(d$LVEDP_out[1], 1)
    clr <- ifelse(val < 12, "green", ifelse(val < 18, "yellow", "red"))
    valueBox(paste(val, "mmHg"), "LVEDP", icon=icon("compress"), color=clr)
  })

  output$vb_egfr <- renderValueBox({
    d <- get_w52()
    val <- round(d$eGFR_out[1], 1)
    clr <- ifelse(val >= 60, "green", ifelse(val >= 30, "yellow", "red"))
    valueBox(paste(val, "mL/min"), "eGFR", icon=icon("tint"), color=clr)
  })

  # ── PK plot ────────────────────────────────────────────────
  output$pk_plot <- renderPlot({
    req(sim_data())
    d <- sim_data() %>%
      select(day, Empa_Cp, Sac_Cp, Val_Cp, Fin_Cp, Furo_Cp) %>%
      pivot_longer(-day, names_to = "Drug", values_to = "Concentration") %>%
      filter(Concentration > 0.01)

    ggplot(d, aes(day, Concentration, color = Drug)) +
      geom_line(linewidth = 1) +
      facet_wrap(~Drug, scales = "free_y", ncol = 2) +
      scale_color_brewer(palette = "Set1") +
      labs(x = "Day", y = "Concentration (ng/mL or μg/mL)",
           title = "Drug Plasma Concentration — Time Profiles") +
      theme_hfpef
  })

  output$te_plot <- renderPlot({
    req(sim_data())
    d <- sim_data() %>%
      filter(abs(week - max(week)) < 0.5) %>%
      summarise(SGLT2 = mean(SGLT2_inh_out)*100,
                NEP   = mean(NEP_inh_out)*100,
                AT1R  = mean(AT1_blk_out)*100,
                MR    = mean(MR_blk_out)*100) %>%
      pivot_longer(everything(), names_to="Target", values_to="Engagement_pct")

    ggplot(d, aes(Target, Engagement_pct, fill=Target)) +
      geom_col(show.legend=FALSE) +
      geom_text(aes(label=sprintf("%.1f%%", Engagement_pct)), vjust=-0.3, size=4) +
      ylim(0, 110) +
      scale_fill_brewer(palette="Dark2") +
      labs(x = "Drug Target", y = "Engagement (%)",
           title = "Target Engagement at Steady State") +
      theme_hfpef
  })

  output$pk_table <- renderTable({
    data.frame(
      Drug        = c("Empagliflozin","Sacubitrilat","Valsartan","Finerenone","Furosemide"),
      `t½ (h)`    = c(12, 12, 10, 20, 2),
      `Vd (L)`    = c("73.8", "35.0", "17.0", "52.0", "12.0"),
      `CL (L/h)`  = c(9.4, 4.5, 1.3, 3.2, 8.0),
      `F (%)`     = c(78, 60, 23, 43, 60),
      `EC50 (ng/mL)` = c(5.0, 5.0, 0.1, 0.9, "2.0 μg/mL"),
      check.names = FALSE
    )
  })

  # ── Neurohumoral plot ──────────────────────────────────────
  output$neuro_plot <- renderPlot({
    req(sim_data())
    d <- sim_data() %>%
      select(week, AngII_n, Aldo_n, ANP_lv, BNP_lv) %>%
      pivot_longer(-week, names_to="Marker", values_to="Value")

    ggplot(d, aes(week, Value, color=Marker)) +
      geom_line(linewidth=1) +
      geom_hline(yintercept=1, linetype=2, color="gray50") +
      scale_color_brewer(palette="Set2") +
      labs(x="Week", y="Normalized Level (1 = normal)",
           title="Neurohumoral Biomarkers") +
      theme_hfpef
  })

  output$cgmp_plot <- renderPlot({
    req(sim_data())
    d <- sim_data() %>%
      select(week, cGMP_lv, PKG_lv) %>%
      pivot_longer(-week, names_to="Marker", values_to="Value")

    ggplot(d, aes(week, Value, color=Marker)) +
      geom_line(linewidth=1.2) +
      geom_hline(yintercept=1, linetype=2, color="gray50") +
      scale_color_manual(values=c(cGMP_lv="#7B1FA2", PKG_lv="#4CAF50")) +
      labs(x="Week", y="Normalized Activity (1=normal)",
           title="cGMP-PKG Pathway Activation") +
      theme_hfpef
  })

  output$struct_plot <- renderPlot({
    req(sim_data())
    d <- sim_data() %>%
      select(week, LVM, Fibr) %>%
      mutate(LVM_norm = LVM/110, Fibr_norm = Fibr/0.30) %>%
      select(week, LVM_norm, Fibr_norm) %>%
      pivot_longer(-week)

    ggplot(d, aes(week, value, color=name)) +
      geom_line(linewidth=1.2) +
      geom_hline(yintercept=1, linetype=2, color="gray50") +
      scale_color_manual(
        values=c(LVM_norm="#FF5722", Fibr_norm="#795548"),
        labels=c("LV Mass (normalized)","LV Fibrosis (normalized)")
      ) +
      labs(x="Week", y="Normalized (1 = baseline)",
           title="Cardiac Structural Remodeling") +
      theme_hfpef
  })

  output$inflam_plot <- renderPlot({
    req(sim_data())
    d <- sim_data() %>%
      select(week, hsCRP_out, IL6_out) %>%
      mutate(CRP_norm = hsCRP_out/3.5, IL6_norm = IL6_out) %>%
      select(week, CRP_norm, IL6_norm) %>%
      pivot_longer(-week)

    ggplot(d, aes(week, value, color=name)) +
      geom_line(linewidth=1.2) +
      geom_hline(yintercept=1, linetype=2, color="gray50") +
      scale_color_manual(
        values=c(CRP_norm="#F44336", IL6_norm="#E91E63"),
        labels=c("hsCRP (normalized)","IL-6 (normalized)")
      ) +
      labs(x="Week", y="Normalized (1 = baseline)",
           title="Inflammation Markers") +
      theme_hfpef
  })

  # ── Clinical endpoint plots ────────────────────────────────
  output$lvedp_plot <- renderPlot({
    req(sim_data())
    ggplot(sim_data(), aes(week, LVEDP_out)) +
      geom_line(color="#1a3a5c", linewidth=1.3) +
      geom_hline(yintercept=12, linetype=2, color="#F44336") +
      annotate("text", x=max(sim_data()$week)*0.8, y=11,
               label="Normal (<12 mmHg)", color="#F44336", size=3.5) +
      labs(x="Week", y="LVEDP (mmHg)",
           title="LV End-Diastolic Pressure") +
      theme_hfpef
  })

  output$ntpbnp_plot <- renderPlot({
    req(sim_data())
    ggplot(sim_data(), aes(week, NTpBNP)) +
      geom_line(color="#E91E63", linewidth=1.3) +
      geom_hline(yintercept=125, linetype=2, color="#FF9800") +
      annotate("text", x=max(sim_data()$week)*0.75, y=115,
               label="HFpEF threshold 125 pg/mL", color="#FF9800", size=3.5) +
      labs(x="Week", y="NT-proBNP (pg/mL)",
           title="NT-proBNP Biomarker") +
      theme_hfpef
  })

  output$egfr_plot <- renderPlot({
    req(sim_data())
    ggplot(sim_data(), aes(week, eGFR_out)) +
      geom_line(color="#4CAF50", linewidth=1.3) +
      geom_hline(yintercept=60, linetype=2, color="#FF9800") +
      geom_hline(yintercept=30, linetype=2, color="#F44336") +
      annotate("text", x=5, y=62, label="CKD G2 boundary", color="#FF9800", size=3) +
      annotate("text", x=5, y=32, label="CKD G3b boundary", color="#F44336", size=3) +
      labs(x="Week", y="eGFR (mL/min/1.73m²)",
           title="Renal Function — eGFR") +
      theme_hfpef
  })

  output$risk_plot <- renderPlot({
    req(sim_data())
    ggplot(sim_data(), aes(week, RiskScore)) +
      geom_area(fill="#90CAF9", alpha=0.4) +
      geom_line(color="#1565C0", linewidth=1.3) +
      labs(x="Week", y="Composite Risk Score (0-100)",
           title="HF Hospitalization Risk Score",
           subtitle="Lower = Better Prognosis") +
      theme_hfpef
  })

  # ── Scenario Comparison ────────────────────────────────────
  compare_data <- reactiveVal(NULL)

  observeEvent(input$run_compare, {
    withProgress(message = "Running 5 treatment scenarios...", value = 0, {
      scn_defs <- list(
        list(nm="Placebo + Furosemide",    col="#888888", empa=0, arni=0, fin=0, furo=1),
        list(nm="Empagliflozin 10mg",      col="#E53935", empa=1, arni=0, fin=0, furo=1),
        list(nm="Sacubitril/Valsartan",    col="#1E88E5", empa=0, arni=1, fin=0, furo=1),
        list(nm="Finerenone 20mg",         col="#43A047", empa=0, arni=0, fin=1, furo=1),
        list(nm="Empa + ARNI (Combo)",     col="#FB8C00", empa=1, arni=1, fin=0, furo=1)
      )
      out_list <- lapply(seq_along(scn_defs), function(i) {
        s <- scn_defs[[i]]
        evs <- build_events(s$empa,10, s$arni,97, s$fin,20, s$furo,40, n_weeks=52)
        pars <- list(use_empa=s$empa, use_arni=s$arni, use_fin=s$fin, use_furo=s$furo)
        run_sim(MOD, pars, evs, n_weeks=52) %>%
          mutate(scenario=s$nm, color=s$col)
      })
      compare_data(bind_rows(out_list))
      setProgress(1)
    })
  })

  output$compare_plot <- renderPlot({
    req(compare_data())
    d <- compare_data()
    colors <- setNames(unique(d$color), unique(d$scenario))

    vars <- c("LVEDP_out"="LVEDP (mmHg)", "NTpBNP"="NT-proBNP (pg/mL)",
              "Fibr"="LV Fibrosis", "RiskScore"="Risk Score")

    dlong <- d %>%
      select(week, scenario, color, all_of(names(vars))) %>%
      pivot_longer(all_of(names(vars)), names_to="Variable", values_to="Value") %>%
      mutate(Variable = vars[Variable])

    ggplot(dlong, aes(week, Value, color=scenario)) +
      geom_line(linewidth=1) +
      facet_wrap(~Variable, scales="free_y", ncol=2) +
      scale_color_manual(values=colors) +
      labs(x="Week", y="", title="52-Week Treatment Scenario Comparison") +
      theme_hfpef
  })

  output$compare_table <- renderDT({
    req(compare_data())
    compare_data() %>%
      filter(abs(week-52) < 0.5) %>%
      group_by(scenario) %>%
      summarise(
        `LVEDP (mmHg)`     = round(mean(LVEDP_out),   1),
        `NT-proBNP (pg/mL)` = round(mean(NTpBNP),     0),
        `LV Fibrosis`      = round(mean(Fibr),         3),
        `eGFR (mL/min)`    = round(mean(eGFR_out),     1),
        `hsCRP (mg/L)`     = round(mean(hsCRP_out),    2),
        `Risk Score`       = round(mean(RiskScore),    1),
        .groups = "drop"
      ) %>%
      datatable(options=list(dom="t", pageLength=10), rownames=FALSE) %>%
      formatStyle("Risk Score",
        background = styleColorBar(c(0,80), "#90CAF9"),
        backgroundSize = "100% 90%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      )
  })

  # ── Biomarker Panel ────────────────────────────────────────
  output$biomarker_ref <- renderTable({
    data.frame(
      Biomarker  = c("NT-proBNP","BNP","hsCRP","eGFR","LV Fibrosis","LVEDP"),
      Normal     = c("<125 pg/mL","<35 pg/mL","<1 mg/L",">60","<15%","<12 mmHg"),
      `HFpEF`    = c("125-500","35-200",">2","30-60","20-40%","15-25"),
      check.names= FALSE
    )
  })

  output$biomarker_w52 <- renderTable({
    req(sim_data())
    d <- sim_data() %>% filter(abs(week - max(week)) < 0.5) %>%
      summarise_if(is.numeric, mean) %>%
      select(NTpBNP, hsCRP_out, eGFR_out, Fibr, LVEDP_out, cGMP_lv, PKG_lv)

    data.frame(
      Biomarker = c("NT-proBNP (pg/mL)","hsCRP (mg/L)","eGFR","Fibrosis","LVEDP (mmHg)","cGMP (norm)","PKG (norm)"),
      Value     = round(c(d$NTpBNP, d$hsCRP_out, d$eGFR_out, d$Fibr*100,
                           d$LVEDP_out, d$cGMP_lv, d$PKG_lv), 2)
    )
  })

  output$radar_plot <- renderPlot({
    req(sim_data())
    d <- sim_data() %>% filter(abs(week - max(week)) < 0.5) %>%
      summarise_if(is.numeric, mean)

    vars_radar <- data.frame(
      name  = c("NT-proBNP","hsCRP","LV Fibrosis","LVEDP","Titin Stiffness",
                "Risk Score","cGMP Act.","eGFR"),
      value = c(
        pmin(d$NTpBNP/400, 2),        # normalized 0-2
        pmin(d$hsCRP_out/6, 2),
        pmin(d$Fibr/0.30, 2),
        pmin(d$LVEDP_out/16, 2),
        pmin(d$Titin/0.45, 2),
        pmin(d$RiskScore/50, 2),
        pmin(d$cGMP_lv, 2),
        pmax(1 - d$eGFR_out/62, 0)
      )
    )
    vars_radar <- rbind(vars_radar, vars_radar[1,])

    ggplot(vars_radar, aes(name, value)) +
      geom_col(aes(fill=value), show.legend=FALSE) +
      coord_polar() +
      scale_fill_gradient(low="#4CAF50", high="#F44336") +
      ylim(0, 2) +
      labs(title="Biomarker Profile (Normalized to Baseline; 1=baseline, >1=worse)") +
      theme_hfpef +
      theme(axis.title=element_blank(), axis.text.y=element_blank())
  })

  output$uricacid_plot <- renderPlot({
    req(sim_data())
    d <- sim_data() %>%
      mutate(UricAcid_proxy = 7.2 * (1 + 0.05*(1 - eGFR_out/62)),
             Creatinine     = 1.2 * (60/eGFR_out))

    dlong <- d %>% select(week, UricAcid_proxy, Creatinine) %>%
      pivot_longer(-week, names_to="Marker", values_to="Value")

    ggplot(dlong, aes(week, Value, color=Marker)) +
      geom_line(linewidth=1.2) +
      scale_color_brewer(palette="Dark2") +
      labs(x="Week", y="Level (mg/dL)",
           title="Uric Acid & Creatinine") +
      theme_hfpef
  })

  # ── Dose-Response ──────────────────────────────────────────
  dr_data <- reactiveVal(NULL)

  observeEvent(input$run_dr, {
    withProgress(message="Running empagliflozin dose-response...", value=0.5, {
      doses <- c(0, 2.5, 5, 10, 25)
      out <- map_dfr(doses, function(d) {
        evs <- build_events(as.integer(d>0), d, 0,0,0,0, 1, 40,
                            n_weeks = input$dr_weeks)
        pars <- list(use_empa=as.integer(d>0), use_arni=0, use_fin=0, use_furo=1)
        run_sim(MOD, pars, evs, n_weeks=input$dr_weeks) %>%
          filter(abs(week - input$dr_weeks) < 0.5) %>%
          summarise_if(is.numeric, mean) %>%
          mutate(dose=d)
      })
      dr_data(out)
      setProgress(1)
    })
  })

  output$dr_plot <- renderPlot({
    req(dr_data())
    ggplot(dr_data(), aes(dose, RiskScore)) +
      geom_line(color="#E53935", linewidth=1.2) +
      geom_point(color="#E53935", size=4) +
      geom_text(aes(label=sprintf("%.1f", RiskScore)), vjust=-1, size=3.5) +
      scale_x_continuous(breaks=c(0,2.5,5,10,25)) +
      labs(x="Empagliflozin Dose (mg)", y="Composite Risk Score",
           title=paste("Dose-Response at Week", input$dr_weeks)) +
      theme_hfpef
  })
}

# ─────────────────────────────────────────────────────────────
# LAUNCH
# ─────────────────────────────────────────────────────────────
if (interactive()) {
  shinyApp(ui, server)
}
