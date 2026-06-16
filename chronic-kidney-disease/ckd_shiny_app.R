################################################################################
# CKD QSP Shiny Dashboard
# Tabs:
#   1. Patient Profile & CKD Stage
#   2. Drug PK Profiles
#   3. Renal Function (eGFR, Proteinuria, Nephron Mass)
#   4. Key PD Biomarkers (RAAS, Inflammation, Fibrosis, MBD, Anemia)
#   5. Clinical Endpoints & Scenario Comparison
#   6. Cardiovascular Risk & CKD-MBD Panel
################################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(patchwork)
library(tidyr)

# ─────────────────────────────────────────────────────────────
# INLINE MODEL
# ─────────────────────────────────────────────────────────────
ckd_model_code <- '
$PARAM @annotated
eGFR0:45 : Baseline eGFR UACR0:300 : Baseline UACR Hb0:11.0 : Baseline Hb
SBP0:150 : Baseline SBP T2DM_flag:1 : T2DM flag
kNephron:0.0004 : Nephron loss rate kFib:0.0008 : Fibrosis rate
kGS:0.0003 : Glomerulosclerosis kEGFR_loss:0.0006 : eGFR decline
kAngII_prod:0.5 : AngII prod kAngII_deg:0.5 : AngII deg
kAldo_prod:0.3 : Aldo prod kAldo_deg:0.3 : Aldo deg
AngII_ss:1.0 : AngII SS Aldo_ss:1.0 : Aldo SS
IL6_ss:1.0 : IL6 SS TNFa_ss:1.0 : TNFa SS
kIL6_prod:0.1 : IL6 prod kIL6_deg:0.15 : IL6 deg
kTNFa_prod:0.08 : TNFa prod kTNFa_deg:0.12 : TNFa deg
kMacro:0.05 : Macrophage in kMacro_out:0.08 : Macrophage out NF_kB_base:1.0 : NF-kB
TGFb1_ss:1.0 : TGFb SS kTGF_prod:0.06 : TGFb prod kTGF_deg:0.08 : TGFb deg
kCol_prod:0.04 : Collagen prod kCol_deg:0.01 : Collagen deg Smad7_inh:0.3 : Smad7
kEMT:0.02 : EMT rate
FGF23_ss:100 : FGF23 SS PTH_ss:65 : PTH SS Klotho_ss:1.0 : Klotho SS
Phos_ss:3.8 : Phos SS kFGF23_prod:0.08 : FGF23 prod kFGF23_deg:0.08 : FGF23 deg
kPTH_prod:0.12 : PTH prod kPTH_deg:0.15 : PTH deg
kKlotho_syn:0.05 : Klotho syn kKlotho_deg:0.05 : Klotho deg
kPhos_retent:0.02 : Phos retention kVitD_prod:0.04 : VitD prod kVitD_deg:0.06 : VitD deg
EPO_ss:12 : EPO SS Hep_ss:25 : Hep SS
kEPO_prod:0.10 : EPO prod kEPO_deg:0.30 : EPO deg
kHb_prod:0.15 : Hb prod kHb_deg:0.008 : Hb deg
kHep_prod:0.08 : Hep prod kHep_deg:0.05 : Hep deg
kLVH_prod:0.001 : LVH prod kLVH_reg:0.0005 : LVH reg kVC:0.002 : VC rate
F_ACEi:0.28 : F ACEi ka_ACEi:0.8 : ka ACEi CL_ACEi:2.5 : CL ACEi V_ACEi:12 : V ACEi
IC50_ACEi:0.5 : IC50 ACEi Imax_ACEi:0.90 : Imax ACEi
F_ARB:0.33 : F ARB ka_ARB:0.9 : ka ARB CL_ARB:7.0 : CL ARB V_ARB:34 : V ARB
IC50_ARB:0.3 : IC50 ARB Imax_ARB:0.85 : Imax ARB
F_Fin:0.44 : F Fin ka_Fin:1.2 : ka Fin CL_Fin:5.0 : CL Fin V1_Fin:22 : V1 Fin
V2_Fin:18 : V2 Fin Q_Fin:3.0 : Q Fin IC50_Fin:0.8 : IC50 Fin Imax_Fin:0.90 : Imax Fin
F_SGLT2:0.78 : F SGLT2 ka_SGLT2:1.5 : ka SGLT2 CL_SGLT2:7.2 : CL SGLT2 V_SGLT2:22 : V SGLT2
IC50_SGLT2:0.6 : IC50 SGLT2 Imax_SGLT2:0.88 : Imax SGLT2
F_ESA:0.62 : F ESA ka_ESA:0.016 : ka ESA CL_ESA:0.12 : CL ESA V_ESA:8.0 : V ESA
Emax_ESA:3.0 : Emax ESA EC50_ESA:10 : EC50 ESA
F_PHI:0.88 : F PHI ka_PHI:1.1 : ka PHI CL_PHI:3.5 : CL PHI V_PHI:16 : V PHI
Imax_PHI:0.80 : Imax PHI IC50_PHI:1.5 : IC50 PHI

$CMT ACEi_gut ACEi_c ARB_gut ARB_c Fin_gut Fin_c Fin_p
SGLT2_gut SGLT2_c ESA_sc ESA_c PHI_gut PHI_c
Nephron eGFR UACR_st AngII Aldo IL6 TNFa Macro TGFb Collagen
FGF23 Klotho Phos PTH VitD EPO Hepcidin Hemoglobin LVH_idx VC_idx BP

$MAIN
double C_ACEi  = ACEi_c  / V_ACEi;
double C_ARB   = ARB_c   / V_ARB;
double C_Fin   = Fin_c   / V1_Fin;
double C_SGLT2 = SGLT2_c / V_SGLT2;
double C_ESA   = ESA_c;
double C_PHI   = PHI_c   / V_PHI;
double E_ACEi  = Imax_ACEi  * C_ACEi  / (IC50_ACEi  + C_ACEi);
double E_ARB   = Imax_ARB   * C_ARB   / (IC50_ARB   + C_ARB);
double E_Fin   = Imax_Fin   * C_Fin   / (IC50_Fin   + C_Fin);
double E_SGLT2 = Imax_SGLT2 * C_SGLT2 / (IC50_SGLT2 + C_SGLT2);
double E_ESA   = Emax_ESA   * C_ESA   / (EC50_ESA   + C_ESA);
double E_PHI   = Imax_PHI   * C_PHI   / (IC50_PHI   + C_PHI);
double E_RAS   = 1.0 - (1.0 - E_ACEi) * (1.0 - E_ARB);
double AngII_norm = AngII / AngII_ss;
double Aldo_norm  = Aldo  / Aldo_ss;
double IL6_norm   = IL6   / IL6_ss;
double TGFb_norm  = TGFb  / TGFb1_ss;
double FGF23_norm = FGF23 / FGF23_ss;
double PTH_norm   = PTH   / PTH_ss;
double Uremia_tox = (1.0 - eGFR / eGFR0);
double GFP_effect = 1.0 + 0.5 * AngII_norm * (1.0 - E_RAS);
double SGLT2_GH_red = 0.3 * E_SGLT2;
double NF_kB_act = NF_kB_base * (1.0 + 0.5 * AngII_norm + 0.3 * Aldo_norm + 0.4 * Uremia_tox)
                 * (1.0 - 0.5 * E_RAS) * (1.0 - 0.3 * E_Fin);
double HIF_stim  = 1.0 + 2.5 * E_PHI;
double Smad23_act = TGFb_norm * (1.0 - Smad7_inh) * (1.0 - 0.4 * E_Fin);

$ODE
dxdt_ACEi_gut  = -ka_ACEi * ACEi_gut;
dxdt_ACEi_c    = F_ACEi * ka_ACEi * ACEi_gut - CL_ACEi * C_ACEi;
dxdt_ARB_gut   = -ka_ARB * ARB_gut;
dxdt_ARB_c     = F_ARB * ka_ARB * ARB_gut - CL_ARB * C_ARB;
dxdt_Fin_gut   = -ka_Fin * Fin_gut;
dxdt_Fin_c     = F_Fin * ka_Fin * Fin_gut - (CL_Fin + Q_Fin) * C_Fin + Q_Fin * (Fin_p / V2_Fin);
dxdt_Fin_p     = Q_Fin * C_Fin * V1_Fin - Q_Fin * (Fin_p / V2_Fin);
dxdt_SGLT2_gut = -ka_SGLT2 * SGLT2_gut;
dxdt_SGLT2_c   = F_SGLT2 * ka_SGLT2 * SGLT2_gut - CL_SGLT2 * C_SGLT2;
dxdt_ESA_sc    = -ka_ESA * ESA_sc;
dxdt_ESA_c     = F_ESA * ka_ESA * ESA_sc - CL_ESA * C_ESA;
dxdt_PHI_gut   = -ka_PHI * PHI_gut;
dxdt_PHI_c     = F_PHI * ka_PHI * PHI_gut - CL_PHI * C_PHI;
double nephron_loss = kNephron * Nephron * (GFP_effect - SGLT2_GH_red) * (1.0 + 0.5 * Collagen);
dxdt_Nephron   = -nephron_loss;
double eGFR_target = eGFR0 * Nephron * (1.0 + 0.12 * E_SGLT2 + 0.08 * E_RAS);
dxdt_eGFR      = 0.05 * (eGFR_target - eGFR);
double UACR_ss_cur = UACR0 * (GFP_effect - SGLT2_GH_red)
                   * (1.0 - 0.45 * E_RAS) * (1.0 - 0.30 * E_SGLT2) * (1.0 - 0.25 * E_Fin);
dxdt_UACR_st   = 0.02 * (UACR_ss_cur - UACR_st);
double AngII_ss_new = AngII_ss * (1.0 - E_RAS);
dxdt_AngII     = kAngII_prod * AngII_ss_new - kAngII_deg * AngII;
double Aldo_drive = Aldo_ss * AngII_norm * (1.0 - 0.5 * E_ACEi);
dxdt_Aldo      = kAldo_prod * Aldo_drive - kAldo_deg * Aldo;
double BP_ss_target = SBP0 - 18.0 * E_RAS - 5.0 * E_SGLT2 - 3.0 * E_Fin
                    + 20.0 * (AngII_norm - 1.0) + 5.0 * Uremia_tox;
dxdt_BP        = 0.005 * (BP_ss_target - BP);
double Mac_drive  = kMacro * (1.0 + 0.5 * AngII_norm + 0.5 * Uremia_tox);
dxdt_Macro     = Mac_drive - kMacro_out * Macro;
double IL6_prod   = kIL6_prod * NF_kB_act * (1.0 + 0.3 * Macro);
dxdt_IL6       = IL6_prod - kIL6_deg * IL6;
double TNFa_prod  = kTNFa_prod * NF_kB_act * (1.0 + 0.4 * Macro);
dxdt_TNFa      = TNFa_prod - kTNFa_deg * TNFa;
double TGFb_prod  = kTGF_prod * (AngII_norm + 0.3 * IL6_norm + 0.2 * Macro)
                  * (1.0 - 0.35 * E_Fin) * (1.0 - 0.15 * E_RAS);
dxdt_TGFb      = TGFb_prod - kTGF_deg * TGFb;
dxdt_Collagen  = kCol_prod * Smad23_act - kCol_deg * Collagen;
double Phos_target = Phos_ss * (1.0 + kPhos_retent * (eGFR0 / (eGFR + 1.0) - 1.0));
dxdt_Phos      = 0.01 * (Phos_target - Phos);
double Klotho_target = Klotho_ss * (eGFR / eGFR0);
dxdt_Klotho    = kKlotho_syn * Klotho_target - kKlotho_deg * Klotho;
double FGF23_drive = kFGF23_prod * (Phos / Phos_ss) * (eGFR0 / (eGFR + 0.1));
dxdt_FGF23     = FGF23_drive - kFGF23_deg * FGF23;
double VitD_ss_cur = 30.0 * (eGFR / eGFR0) / (1.0 + 0.02 * (FGF23 - FGF23_ss));
dxdt_VitD      = kVitD_prod * VitD_ss_cur - kVitD_deg * VitD;
double VitD_norm   = VitD / 30.0;
double PTH_stim    = kPTH_prod * (1.0 + 0.5 * (Phos / Phos_ss - 1.0))
                               / (1.0 + 0.5 * VitD_norm) / (1.0 + 0.3 * Klotho);
dxdt_PTH       = PTH_stim - kPTH_deg * PTH;
double EPO_prod   = kEPO_prod * EPO_ss * (eGFR / eGFR0) * HIF_stim + E_ESA * kEPO_prod * EPO_ss;
dxdt_EPO       = EPO_prod - kEPO_deg * EPO;
double Hep_prod   = kHep_prod * (1.0 + 1.5 * IL6_norm);
dxdt_Hepcidin  = Hep_prod - kHep_deg * Hepcidin;
double EPO_stim   = (EPO / EPO_ss) + E_ESA + 0.8 * E_PHI;
double Hb_prod    = kHb_prod * EPO_stim / (1.0 + 0.2 * (Hepcidin / Hep_ss));
double Hb_loss    = kHb_deg  * Hemoglobin * (1.0 + 0.3 * Uremia_tox);
dxdt_Hemoglobin = Hb_prod - Hb_loss;
double LVH_drive  = kLVH_prod * (BP / SBP0) * Aldo_norm * (1.0 - 0.4 * E_Fin);
dxdt_LVH_idx   = LVH_drive - kLVH_reg * LVH_idx * (E_RAS + 0.5 * E_Fin + 0.3 * E_SGLT2);
double VC_drive   = kVC * (Phos / Phos_ss) * PTH_norm * (1.0 + 0.5 * Uremia_tox);
dxdt_VC_idx    = VC_drive * (1.0 - VC_idx);

$TABLE
double eGFR_out   = eGFR;
double UACR_out   = UACR_st;
double Hb_out     = Hemoglobin;
double PTH_out    = PTH;
double SBP_out    = BP;
double FGF23_out  = FGF23;
double Phos_out   = Phos;
double Col_out    = Collagen;
double LVH_out    = LVH_idx;
double CV_risk    = 0.3 * (BP / SBP0 - 0.5) + 0.3 * VC_idx + 0.2 * (LVH_idx - 0.5) + 0.2 * (UACR_st / UACR0);
double CKD_stage;
if      (eGFR >= 90) CKD_stage = 1;
else if (eGFR >= 60) CKD_stage = 2;
else if (eGFR >= 30) CKD_stage = 3;
else if (eGFR >= 15) CKD_stage = 4;
else                  CKD_stage = 5;

$CAPTURE eGFR_out UACR_out Hb_out PTH_out SBP_out FGF23_out Phos_out
Col_out LVH_out VC_idx Nephron CV_risk CKD_stage AngII Aldo IL6 TNFa TGFb
Klotho VitD EPO Hepcidin
'

ckd_mod <- suppressMessages(mrgsolve::mcode("CKD_Shiny", ckd_model_code))

# ─────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────
ckd_init <- function(egfr0, uacr0, hb0, sbp0) {
  list(
    Nephron = 1.0, eGFR = egfr0, UACR_st = uacr0,
    AngII = 1.0, Aldo = 1.0, IL6 = 1.0, TNFa = 1.0, Macro = 1.0,
    TGFb = 1.0, Collagen = 0.2, FGF23 = 100, Klotho = 1.0,
    Phos = 3.8, PTH = 65, VitD = 30.0, EPO = 12, Hepcidin = 25,
    Hemoglobin = hb0, LVH_idx = 1.0, VC_idx = 0.05, BP = sbp0
  )
}

run_sim <- function(egfr0, uacr0, hb0, sbp0, duration_yr,
                    use_acei, dose_acei,
                    use_arb,  dose_arb,
                    use_fin,  dose_fin,
                    use_sglt2, dose_sglt2,
                    use_esa,  dose_esa,
                    use_phi,  dose_phi) {

  evs <- list()
  if (use_acei  && dose_acei  > 0) evs[["acei"]]  <- mrgsolve::ev(amt = dose_acei,  cmt = "ACEi_gut",  ii = 24, addl = duration_yr * 365 - 1)
  if (use_arb   && dose_arb   > 0) evs[["arb"]]   <- mrgsolve::ev(amt = dose_arb,   cmt = "ARB_gut",   ii = 24, addl = duration_yr * 365 - 1)
  if (use_fin   && dose_fin   > 0) evs[["fin"]]   <- mrgsolve::ev(amt = dose_fin,   cmt = "Fin_gut",   ii = 24, addl = duration_yr * 365 - 1)
  if (use_sglt2 && dose_sglt2 > 0) evs[["sglt2"]] <- mrgsolve::ev(amt = dose_sglt2, cmt = "SGLT2_gut", ii = 24, addl = duration_yr * 365 - 1)
  if (use_esa   && dose_esa   > 0) evs[["esa"]]   <- mrgsolve::ev(amt = dose_esa,   cmt = "ESA_sc",    ii = 24 * 7 / 3, addl = duration_yr * 52 * 3 - 1)
  if (use_phi   && dose_phi   > 0) evs[["phi"]]   <- mrgsolve::ev(amt = dose_phi,   cmt = "PHI_gut",   ii = 24 * 7 / 3, addl = duration_yr * 52 * 3 - 1)

  if (length(evs) == 0) evs[["null"]] <- mrgsolve::ev(amt = 0, cmt = "ACEi_gut", time = 0)
  combined_ev <- Reduce(mrgsolve::ev_c, evs)

  ckd_mod %>%
    mrgsolve::init(ckd_init(egfr0, uacr0, hb0, sbp0)) %>%
    mrgsolve::ev(combined_ev) %>%
    mrgsolve::mrgsim(end = duration_yr * 365 * 24, delta = 24, obsonly = TRUE) %>%
    as.data.frame() %>%
    mutate(Day = time / 24, Year = Day / 365)
}

theme_ckd <- function() {
  theme_bw(base_size = 12) +
    theme(
      legend.position   = "bottom",
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "#1A5276"),
      strip.text        = element_text(color = "white", face = "bold")
    )
}

ckd_stage_color <- function(egfr) {
  case_when(
    egfr >= 90 ~ "#27AE60",
    egfr >= 60 ~ "#2ECC71",
    egfr >= 45 ~ "#F39C12",
    egfr >= 30 ~ "#E67E22",
    egfr >= 15 ~ "#E74C3C",
    TRUE       ~ "#922B21"
  )
}

# ─────────────────────────────────────────────────────────────
# UI
# ─────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = "CKD QSP Dashboard",
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",         tabName = "tab_patient", icon = icon("user-md")),
      menuItem("Drug PK Profiles",        tabName = "tab_pk",      icon = icon("flask")),
      menuItem("Renal Function",          tabName = "tab_renal",   icon = icon("filter")),
      menuItem("PD Biomarkers",           tabName = "tab_pd",      icon = icon("chart-line")),
      menuItem("Clinical Endpoints",      tabName = "tab_clinical",icon = icon("stethoscope")),
      menuItem("CV Risk & CKD-MBD",      tabName = "tab_cv",      icon = icon("heartbeat"))
    ),
    hr(),
    h5("Patient Parameters", style = "padding-left:15px; color:#BDC3C7;"),
    sliderInput("egfr0",  "Baseline eGFR (mL/min)", 10, 90, 45, step = 5),
    sliderInput("uacr0",  "Baseline UACR (mg/g)",  30, 3000, 300, step = 30),
    sliderInput("hb0",    "Baseline Hb (g/dL)",    7, 14, 11.0, step = 0.5),
    sliderInput("sbp0",   "Baseline SBP (mmHg)",   110, 200, 150, step = 5),
    sliderInput("dur_yr", "Simulation (years)",    1, 5, 3),
    hr(),
    h5("Drug Selection", style = "padding-left:15px; color:#BDC3C7;"),

    fluidRow(
      column(6, checkboxInput("use_acei",  "ACEi",   value = TRUE)),
      column(6, sliderInput( "dose_acei", NULL,       1, 20, 10, step = 1))
    ),
    fluidRow(
      column(6, checkboxInput("use_arb",   "ARB",    value = FALSE)),
      column(6, sliderInput( "dose_arb",  NULL,      25, 200, 100, step = 25))
    ),
    fluidRow(
      column(6, checkboxInput("use_fin",   "Finerenone",  value = FALSE)),
      column(6, sliderInput( "dose_fin",  NULL,  5, 40, 20, step = 5))
    ),
    fluidRow(
      column(6, checkboxInput("use_sglt2", "SGLT2i", value = FALSE)),
      column(6, sliderInput( "dose_sglt2", NULL, 5, 25, 10, step = 5))
    ),
    fluidRow(
      column(6, checkboxInput("use_esa",   "ESA",    value = FALSE)),
      column(6, sliderInput( "dose_esa",  NULL, 1000, 8000, 4000, step = 1000))
    ),
    fluidRow(
      column(6, checkboxInput("use_phi",   "HIF-PHI", value = FALSE)),
      column(6, sliderInput( "dose_phi",  NULL,  50, 200, 100, step = 50))
    ),
    hr(),
    actionButton("run_sim", "Run Simulation", class = "btn-primary btn-block",
                 icon = icon("play"))
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box-header { background-color: #2C3E50 !important; color: white !important; }
      .content-wrapper { background-color: #F0F3F4; }
    "))),

    tabItems(

      # ── TAB 1: Patient Profile ──────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          valueBoxOutput("vb_egfr",  width = 3),
          valueBoxOutput("vb_stage", width = 3),
          valueBoxOutput("vb_uacr",  width = 3),
          valueBoxOutput("vb_hb",    width = 3)
        ),
        fluidRow(
          valueBoxOutput("vb_sbp",   width = 3),
          valueBoxOutput("vb_pth",   width = 3),
          valueBoxOutput("vb_fgf23", width = 3),
          valueBoxOutput("vb_phos",  width = 3)
        ),
        fluidRow(
          box(title = "Patient Summary", width = 6, solidHeader = TRUE, status = "primary",
            tableOutput("patient_summary_tbl")
          ),
          box(title = "CKD Classification (KDIGO 2024)", width = 6, solidHeader = TRUE, status = "info",
            plotOutput("ckd_kdigo_plot", height = "280px")
          )
        ),
        fluidRow(
          box(title = "Trajectory Overview", width = 12, solidHeader = TRUE, status = "primary",
            plotlyOutput("overview_plot", height = "350px")
          )
        )
      ),

      # ── TAB 2: Drug PK Profiles ────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "RAAS Blocker PK (ACEi / ARB)", width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("pk_raas", height = "300px")
          ),
          box(title = "Finerenone PK (2-compartment)", width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("pk_fin", height = "300px")
          )
        ),
        fluidRow(
          box(title = "SGLT2 Inhibitor PK (Dapagliflozin)", width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("pk_sglt2", height = "300px")
          ),
          box(title = "Anemia Treatment PK (ESA / HIF-PHI)", width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("pk_anemia", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Drug Effect Summary at Steady-State", width = 12, solidHeader = TRUE, status = "info",
            DT::dataTableOutput("pk_effect_tbl")
          )
        )
      ),

      # ── TAB 3: Renal Function ──────────────────────────────
      tabItem(tabName = "tab_renal",
        fluidRow(
          box(title = "eGFR Trajectory", width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("renal_egfr", height = "320px")
          ),
          box(title = "UACR / Proteinuria", width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("renal_uacr", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Nephron Mass (Fraction Remaining)", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("renal_nephron", height = "300px")
          ),
          box(title = "Renal Fibrosis Index (Collagen)", width = 6, solidHeader = TRUE, status = "danger",
            plotlyOutput("renal_fibrosis", height = "300px")
          )
        )
      ),

      # ── TAB 4: PD Biomarkers ──────────────────────────────
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(title = "RAAS Markers (Ang II, Aldosterone)", width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("pd_raas", height = "280px")
          ),
          box(title = "Inflammation (IL-6, TNF-α, TGF-β1)", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("pd_inflam", height = "280px")
          )
        ),
        fluidRow(
          box(title = "CKD-MBD (FGF-23, Klotho, PTH, Phosphate, VitD)", width = 6, solidHeader = TRUE, status = "danger",
            plotlyOutput("pd_mbd", height = "280px")
          ),
          box(title = "Anemia (EPO, Hepcidin, Hemoglobin)", width = 6, solidHeader = TRUE, status = "success",
            plotlyOutput("pd_anemia", height = "280px")
          )
        )
      ),

      # ── TAB 5: Clinical Endpoints ─────────────────────────
      tabItem(tabName = "tab_clinical",
        fluidRow(
          box(title = "Scenario Comparison Table", width = 12, solidHeader = TRUE, status = "primary",
            DT::dataTableOutput("scenario_tbl")
          )
        ),
        fluidRow(
          box(title = "eGFR Change from Baseline (%)", width = 6, solidHeader = TRUE, status = "info",
            plotlyOutput("endpoint_egfr_pct", height = "300px")
          ),
          box(title = "UACR Change from Baseline (%)", width = 6, solidHeader = TRUE, status = "info",
            plotlyOutput("endpoint_uacr_pct", height = "300px")
          )
        ),
        fluidRow(
          box(title = "CKD Stage Progression Timeline", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("ckd_stage_plot", height = "300px")
          ),
          box(title = "Composite Endpoint Risk Score", width = 6, solidHeader = TRUE, status = "danger",
            plotlyOutput("cv_risk_plot", height = "300px")
          )
        )
      ),

      # ── TAB 6: CV Risk & CKD-MBD ─────────────────────────
      tabItem(tabName = "tab_cv",
        fluidRow(
          box(title = "Blood Pressure Trajectory", width = 6, solidHeader = TRUE, status = "danger",
            plotlyOutput("cv_bp", height = "280px")
          ),
          box(title = "LVH Index (Left Ventricular Hypertrophy)", width = 6, solidHeader = TRUE, status = "danger",
            plotlyOutput("cv_lvh", height = "280px")
          )
        ),
        fluidRow(
          box(title = "Vascular Calcification Index", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("cv_vc", height = "280px")
          ),
          box(title = "PTH (Secondary Hyperparathyroidism)", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("cv_pth", height = "280px")
          )
        ),
        fluidRow(
          box(title = "FGF-23 & Klotho Axis", width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("mbd_fgf23_klotho", height = "280px")
          ),
          box(title = "Active Vitamin D (Calcitriol)", width = 6, solidHeader = TRUE, status = "success",
            plotlyOutput("mbd_vitd", height = "280px")
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

  # Reactive simulation result
  sim_result <- eventReactive(input$run_sim, {
    withProgress(message = "Running CKD QSP simulation...", {
      run_sim(
        egfr0     = input$egfr0,  uacr0   = input$uacr0,
        hb0       = input$hb0,   sbp0    = input$sbp0,
        duration_yr = input$dur_yr,
        use_acei  = input$use_acei,  dose_acei  = input$dose_acei,
        use_arb   = input$use_arb,   dose_arb   = input$dose_arb,
        use_fin   = input$use_fin,   dose_fin   = input$dose_fin,
        use_sglt2 = input$use_sglt2, dose_sglt2 = input$dose_sglt2,
        use_esa   = input$use_esa,   dose_esa   = input$dose_esa,
        use_phi   = input$use_phi,   dose_phi   = input$dose_phi
      )
    })
  }, ignoreNULL = FALSE)

  # Multi-scenario comparison (fixed scenarios)
  multi_results <- reactive({
    req(input$run_sim >= 0)
    e0 <- input$egfr0; u0 <- input$uacr0; h0 <- input$hb0; s0 <- input$sbp0
    dy <- input$dur_yr

    scens <- list(
      "No Treatment"          = list(use_acei=F,dose_acei=0, use_arb=F,dose_arb=0,
                                     use_fin=F,dose_fin=0,   use_sglt2=F,dose_sglt2=0,
                                     use_esa=F,dose_esa=0,   use_phi=F,dose_phi=0),
      "ACEi Monotherapy"      = list(use_acei=T,dose_acei=10,use_arb=F,dose_arb=0,
                                     use_fin=F,dose_fin=0,   use_sglt2=F,dose_sglt2=0,
                                     use_esa=F,dose_esa=0,   use_phi=F,dose_phi=0),
      "ACEi + Finerenone"     = list(use_acei=T,dose_acei=10,use_arb=F,dose_arb=0,
                                     use_fin=T,dose_fin=20,  use_sglt2=F,dose_sglt2=0,
                                     use_esa=F,dose_esa=0,   use_phi=F,dose_phi=0),
      "ACEi + SGLT2i"         = list(use_acei=T,dose_acei=10,use_arb=F,dose_arb=0,
                                     use_fin=F,dose_fin=0,   use_sglt2=T,dose_sglt2=10,
                                     use_esa=F,dose_esa=0,   use_phi=F,dose_phi=0),
      "Triple Therapy"        = list(use_acei=T,dose_acei=10,use_arb=F,dose_arb=0,
                                     use_fin=T,dose_fin=20,  use_sglt2=T,dose_sglt2=10,
                                     use_esa=F,dose_esa=0,   use_phi=F,dose_phi=0)
    )
    colors <- c("#E74C3C","#3498DB","#9B59B6","#27AE60","#F39C12")

    bind_rows(lapply(seq_along(scens), function(i) {
      s <- scens[[i]]
      run_sim(e0,u0,h0,s0,dy,
              s$use_acei,s$dose_acei, s$use_arb,s$dose_arb,
              s$use_fin,s$dose_fin,   s$use_sglt2,s$dose_sglt2,
              s$use_esa,s$dose_esa,   s$use_phi,s$dose_phi) %>%
        mutate(Scenario = names(scens)[i], Color = colors[i])
    }))
  })

  # Current sim data shortcuts
  d <- reactive({ sim_result() })
  last_row <- reactive({ d() %>% slice_tail(n = 1) })

  # ── VALUE BOXES ─────────────────────────────────────────────
  output$vb_egfr  <- renderValueBox({
    v <- round(last_row()$eGFR_out, 1)
    valueBox(paste0(v, " mL/min"), "Current eGFR",
             color = if (v >= 60) "green" else if (v >= 30) "yellow" else "red",
             icon = icon("filter"))
  })
  output$vb_stage <- renderValueBox({
    s <- last_row()$CKD_stage
    valueBox(paste0("G", round(s)), "CKD Stage",
             color = c("green","lime","yellow","orange","red","maroon")[min(round(s),5)],
             icon = icon("tachometer-alt"))
  })
  output$vb_uacr  <- renderValueBox({
    v <- round(last_row()$UACR_out, 0)
    valueBox(paste0(v, " mg/g"), "UACR",
             color = if (v < 30) "green" else if (v < 300) "yellow" else "red",
             icon = icon("droplet"))
  })
  output$vb_hb    <- renderValueBox({
    v <- round(last_row()$Hb_out, 1)
    valueBox(paste0(v, " g/dL"), "Hemoglobin",
             color = if (v >= 12) "green" else if (v >= 10) "yellow" else "red",
             icon = icon("tint"))
  })
  output$vb_sbp   <- renderValueBox({
    v <- round(last_row()$SBP_out, 0)
    valueBox(paste0(v, " mmHg"), "SBP",
             color = if (v < 130) "green" else if (v < 160) "yellow" else "red",
             icon = icon("heartbeat"))
  })
  output$vb_pth   <- renderValueBox({
    v <- round(last_row()$PTH_out, 0)
    valueBox(paste0(v, " pg/mL"), "Intact PTH",
             color = if (v < 150) "green" else if (v < 300) "yellow" else "red",
             icon = icon("bone"))
  })
  output$vb_fgf23 <- renderValueBox({
    v <- round(last_row()$FGF23_out, 0)
    valueBox(paste0(v, " pg/mL"), "FGF-23",
             color = if (v < 200) "green" else if (v < 500) "yellow" else "red",
             icon = icon("dna"))
  })
  output$vb_phos  <- renderValueBox({
    v <- round(last_row()$Phos_out, 1)
    valueBox(paste0(v, " mg/dL"), "Serum Phosphate",
             color = if (v < 4.5) "green" else if (v < 5.5) "yellow" else "red",
             icon = icon("atom"))
  })

  # ── PATIENT SUMMARY TABLE ───────────────────────────────────
  output$patient_summary_tbl <- renderTable({
    r <- last_row()
    data.frame(
      Parameter = c("eGFR", "CKD Stage", "UACR", "Hemoglobin",
                    "SBP", "PTH", "FGF-23", "Phosphate",
                    "Collagen/Fibrosis", "LVH Index", "VC Index", "CV Risk"),
      Value = c(
        paste0(round(r$eGFR_out, 1), " mL/min/1.73m²"),
        paste0("G", round(r$CKD_stage)),
        paste0(round(r$UACR_out, 0), " mg/g Cr"),
        paste0(round(r$Hb_out, 1), " g/dL"),
        paste0(round(r$SBP_out, 0), " mmHg"),
        paste0(round(r$PTH_out, 0), " pg/mL"),
        paste0(round(r$FGF23_out, 0), " pg/mL"),
        paste0(round(r$Phos_out, 1), " mg/dL"),
        round(r$Col_out, 3),
        round(r$LVH_out, 3),
        round(r$VC_idx, 3),
        round(r$CV_risk, 3)
      )
    )
  }, striped = TRUE, hover = TRUE)

  # ── KDIGO PLOT ──────────────────────────────────────────────
  output$ckd_kdigo_plot <- renderPlot({
    egfr_now <- last_row()$eGFR_out
    uacr_now <- last_row()$UACR_out

    kdigo_df <- expand.grid(
      eGFR_cat = c("G1 (≥90)", "G2 (60-89)", "G3a (45-59)", "G3b (30-44)", "G4 (15-29)", "G5 (<15)"),
      UACR_cat = c("A1 (<30)", "A2 (30-300)", "A3 (>300)")
    )
    risk_mat <- matrix(
      c("Low","Mod","High",
        "Low","Mod","High",
        "Mod","High","VH",
        "High","VH","VH",
        "VH","VH","VH",
        "VH","VH","VH"),
      nrow = 6, ncol = 3, byrow = TRUE
    )
    kdigo_df$Risk <- c(risk_mat)
    kdigo_df$Fill <- case_when(
      kdigo_df$Risk == "Low"  ~ "#2ECC71",
      kdigo_df$Risk == "Mod"  ~ "#F1C40F",
      kdigo_df$Risk == "High" ~ "#E67E22",
      TRUE                    ~ "#E74C3C"
    )

    egfr_cat <- case_when(
      egfr_now >= 90 ~ "G1 (≥90)",
      egfr_now >= 60 ~ "G2 (60-89)",
      egfr_now >= 45 ~ "G3a (45-59)",
      egfr_now >= 30 ~ "G3b (30-44)",
      egfr_now >= 15 ~ "G4 (15-29)",
      TRUE           ~ "G5 (<15)"
    )
    uacr_cat <- case_when(
      uacr_now < 30  ~ "A1 (<30)",
      uacr_now < 300 ~ "A2 (30-300)",
      TRUE           ~ "A3 (>300)"
    )

    ggplot(kdigo_df, aes(x = UACR_cat, y = eGFR_cat, fill = Fill)) +
      geom_tile(color = "white", linewidth = 1.5) +
      geom_tile(data = subset(kdigo_df, UACR_cat == uacr_cat & eGFR_cat == egfr_cat),
                fill = NA, color = "#2C3E50", linewidth = 3) +
      geom_text(aes(label = Risk), size = 4, fontface = "bold") +
      scale_fill_identity() +
      scale_y_discrete(limits = rev) +
      labs(x = "Albuminuria Category", y = "GFR Category",
           title = "KDIGO Risk Heat Map (■ = current patient)") +
      theme_ckd() +
      theme(legend.position = "none")
  })

  # ── OVERVIEW PLOT ───────────────────────────────────────────
  output$overview_plot <- renderPlotly({
    df <- d()
    p <- ggplot(df, aes(x = Year)) +
      geom_line(aes(y = eGFR_out,  color = "eGFR (mL/min)"),       linewidth = 1.0) +
      geom_line(aes(y = Hb_out * 5, color = "Hb × 5 (g/dL)"),      linewidth = 1.0, linetype = "dashed") +
      geom_line(aes(y = SBP_out / 10, color = "SBP ÷ 10 (mmHg)"),  linewidth = 1.0, linetype = "dotdash") +
      geom_hline(yintercept = 15, linetype = "dashed", color = "gray60", alpha = 0.7) +
      scale_color_manual(values = c("eGFR (mL/min)"="#2980B9",
                                    "Hb × 5 (g/dL)"="#1ABC9C",
                                    "SBP ÷ 10 (mmHg)"="#E74C3C")) +
      labs(x = "Year", y = "Value (scaled)", title = "Patient Trajectory Overview") +
      theme_ckd()
    ggplotly(p) %>% layout(legend = list(orientation = "h", x = 0, y = -0.2))
  })

  # ── PK PLOTS ────────────────────────────────────────────────
  pk_24h <- reactive({
    # Single-dose PK over 72h for visualization
    run_sim(input$egfr0, input$uacr0, input$hb0, input$sbp0,
            duration_yr = 0.02,
            use_acei = input$use_acei,   dose_acei  = input$dose_acei,
            use_arb  = input$use_arb,    dose_arb   = input$dose_arb,
            use_fin  = input$use_fin,    dose_fin   = input$dose_fin,
            use_sglt2= input$use_sglt2,  dose_sglt2 = input$dose_sglt2,
            use_esa  = input$use_esa,    dose_esa   = input$dose_esa,
            use_phi  = input$use_phi,    dose_phi   = input$dose_phi) %>%
      filter(Day <= 7)
  })

  make_pk_plot <- function(df, vars, labels, colors, title, ylab) {
    long_df <- df %>%
      select(Day, all_of(vars)) %>%
      pivot_longer(-Day, names_to = "Compound", values_to = "Conc") %>%
      mutate(Compound = factor(Compound, levels = vars, labels = labels))
    p <- ggplot(long_df, aes(x = Day * 24, y = Conc, color = Compound)) +
      geom_line(linewidth = 1.0) +
      scale_color_manual(values = setNames(colors, labels)) +
      labs(x = "Hour", y = ylab, title = title) +
      theme_ckd()
    ggplotly(p)
  }

  output$pk_raas <- renderPlotly({
    make_pk_plot(pk_24h(), c("AngII","Aldo"), c("Ang II (norm)","Aldosterone (norm)"),
                 c("#27AE60","#F39C12"), "RAAS Dynamics (First 7 Days)", "Normalized Conc.")
  })

  output$pk_fin <- renderPlotly({
    df <- pk_24h()
    p <- ggplot(df, aes(x = Day * 24)) +
      geom_line(aes(y = Fin_c, color = "Central"), linewidth = 1.0) +
      geom_line(aes(y = Fin_p, color = "Peripheral"), linewidth = 1.0, linetype = "dashed") +
      scale_color_manual(values = c("Central"="#9B59B6","Peripheral"="#D7BDE2")) +
      labs(x = "Hour", y = "Amount (mg)", title = "Finerenone 2-Compartment PK") +
      theme_ckd()
    ggplotly(p)
  })

  output$pk_sglt2 <- renderPlotly({
    make_pk_plot(pk_24h(), c("SGLT2_c"), c("Dapagliflozin (central)"),
                 c("#3498DB"), "SGLT2i PK (First 7 Days)", "Amount (mg)")
  })

  output$pk_anemia <- renderPlotly({
    make_pk_plot(pk_24h(), c("EPO","Hepcidin","Hemoglobin"),
                 c("EPO (mIU/mL)","Hepcidin (ng/mL)","Hemoglobin (g/dL)"),
                 c("#1ABC9C","#E74C3C","#F39C12"),
                 "Anemia Dynamics (First 7 Days)", "Concentration")
  })

  output$pk_effect_tbl <- DT::renderDataTable({
    lr <- last_row()
    data.frame(
      Drug = c("Ramipril (ACEi)","Losartan (ARB)","Finerenone (nsMRA)",
               "Dapagliflozin (SGLT2i)","Epoetin alfa (ESA)","Roxadustat (HIF-PHI)"),
      `Mechanism` = c("ACE inhibition","AT1R blockade","MR blockade (nsSMRA)",
                      "SGLT2 inhibition","EPOR stimulation","PHD inhibition → ↑HIF"),
      `PK Model` = c("1-cmt + prodrug","1-cmt + prodrug","2-cmt","1-cmt","1-cmt SC","1-cmt"),
      `Target Effect` = c("↓Ang II, ↓Aldosterone","↓AT1R activation","↓MR activation, ↓fibrosis",
                          "↓SGLT2, TGF, ↓UACR","↑Erythropoiesis","↑HIF-1α → ↑EPO"),
      `Clinical Benefit` = c("↓UACR, ↓BP","↓UACR, ↓BP","↓UACR, ↓CV events","↓eGFR loss, ↓UACR",
                              "↑Hb","↑Hb, ↓hepcidin dependence")
    )
  }, options = list(pageLength = 10, dom = 't'), rownames = FALSE)

  # ── RENAL FUNCTION PLOTS ────────────────────────────────────
  output$renal_egfr <- renderPlotly({
    df <- d()
    p <- ggplot(df, aes(x = Year, y = eGFR_out)) +
      geom_line(color = "#2980B9", linewidth = 1.2) +
      geom_ribbon(aes(ymin = eGFR_out * 0.9, ymax = eGFR_out * 1.1),
                  fill = "#2980B9", alpha = 0.15) +
      geom_hline(yintercept = c(15, 30, 45, 60, 90),
                 linetype = "dashed", color = "gray60", alpha = 0.5) +
      labs(x = "Year", y = "eGFR (mL/min/1.73m²)", title = "eGFR Over Time") +
      theme_ckd()
    ggplotly(p)
  })

  output$renal_uacr <- renderPlotly({
    df <- d()
    p <- ggplot(df, aes(x = Year, y = UACR_out)) +
      geom_line(color = "#E74C3C", linewidth = 1.2) +
      geom_hline(yintercept = 30,  linetype = "dashed", color = "#27AE60") +
      geom_hline(yintercept = 300, linetype = "dashed", color = "#E67E22") +
      scale_y_log10() +
      labs(x = "Year", y = "UACR (mg/g Cr, log scale)", title = "Proteinuria (UACR)") +
      theme_ckd()
    ggplotly(p)
  })

  output$renal_nephron <- renderPlotly({
    df <- d()
    p <- ggplot(df, aes(x = Year, y = Nephron * 100)) +
      geom_area(fill = "#F39C12", alpha = 0.4) +
      geom_line(color = "#F39C12", linewidth = 1.2) +
      labs(x = "Year", y = "Remaining Nephrons (%)", title = "Functional Nephron Mass") +
      ylim(0, 100) +
      theme_ckd()
    ggplotly(p)
  })

  output$renal_fibrosis <- renderPlotly({
    df <- d()
    p <- ggplot(df, aes(x = Year, y = Col_out)) +
      geom_area(fill = "#8E44AD", alpha = 0.4) +
      geom_line(color = "#8E44AD", linewidth = 1.2) +
      labs(x = "Year", y = "Fibrosis Index (0–1)", title = "Renal Fibrosis (TGF-β/Collagen)") +
      ylim(0, 1) +
      theme_ckd()
    ggplotly(p)
  })

  # ── PD BIOMARKER PLOTS ──────────────────────────────────────
  pd_long <- reactive({
    d() %>%
      select(Year, AngII, Aldo, IL6, TNFa, TGFb,
             FGF23_out, Klotho, PTH_out, Phos_out, VitD,
             EPO, Hepcidin, Hemoglobin = Hb_out)
  })

  output$pd_raas <- renderPlotly({
    df <- pd_long() %>%
      select(Year, `Ang II (norm)` = AngII, `Aldosterone (norm)` = Aldo) %>%
      pivot_longer(-Year)
    p <- ggplot(df, aes(x = Year, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray60") +
      scale_color_manual(values = c("#27AE60","#F39C12")) +
      labs(x = "Year", y = "Normalized Level", title = "RAAS Cascade") +
      theme_ckd()
    ggplotly(p)
  })

  output$pd_inflam <- renderPlotly({
    df <- pd_long() %>%
      select(Year, `IL-6 (norm)` = IL6, `TNF-α (norm)` = TNFa, `TGF-β1 (norm)` = TGFb) %>%
      pivot_longer(-Year)
    p <- ggplot(df, aes(x = Year, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("#E74C3C","#F39C12","#8E44AD")) +
      labs(x = "Year", y = "Normalized Level", title = "Inflammatory Mediators") +
      theme_ckd()
    ggplotly(p)
  })

  output$pd_mbd <- renderPlotly({
    df <- pd_long() %>%
      select(Year, `FGF-23 (pg/mL)` = FGF23_out, `PTH (pg/mL)` = PTH_out,
             `Phos (mg/dL×10)` = Phos_out, `VitD (pg/mL)` = VitD,
             `Klotho (norm×50)` = Klotho) %>%
      mutate(`Phos (mg/dL×10)` = `Phos (mg/dL×10)` * 10,
             `Klotho (norm×50)` = `Klotho (norm×50)` * 50) %>%
      pivot_longer(-Year)
    p <- ggplot(df, aes(x = Year, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("#E74C3C","#F1948A","#922B21","#82E0AA","#5DADE2")) +
      labs(x = "Year", y = "Value (mixed units)", title = "CKD-MBD Markers") +
      theme_ckd()
    ggplotly(p)
  })

  output$pd_anemia <- renderPlotly({
    df <- pd_long() %>%
      select(Year, `EPO (mIU/mL)` = EPO, `Hepcidin (ng/mL)` = Hepcidin,
             `Hemoglobin (g/dL×5)` = Hemoglobin) %>%
      mutate(`Hemoglobin (g/dL×5)` = `Hemoglobin (g/dL×5)` * 5) %>%
      pivot_longer(-Year)
    p <- ggplot(df, aes(x = Year, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("#1ABC9C","#E74C3C","#3498DB")) +
      labs(x = "Year", y = "Concentration", title = "Anemia of CKD") +
      theme_ckd()
    ggplotly(p)
  })

  # ── SCENARIO COMPARISON ─────────────────────────────────────
  output$scenario_tbl <- DT::renderDataTable({
    mr <- multi_results()
    end_time <- max(mr$Day)
    mr %>%
      filter(abs(Day - end_time) < 2) %>%
      group_by(Scenario) %>%
      summarise(
        `eGFR (mL/min)` = round(mean(eGFR_out), 1),
        `UACR (mg/g)`   = round(mean(UACR_out), 0),
        `SBP (mmHg)`    = round(mean(SBP_out),  0),
        `Hb (g/dL)`     = round(mean(Hb_out),   1),
        `PTH (pg/mL)`   = round(mean(PTH_out),  0),
        `CKD Stage`     = round(mean(CKD_stage), 1),
        `CV Risk`       = round(mean(CV_risk),   3),
        .groups = "drop"
      )
  }, options = list(pageLength = 10, dom = 't'), rownames = FALSE)

  output$endpoint_egfr_pct <- renderPlotly({
    df <- multi_results() %>%
      group_by(Scenario, Color) %>%
      mutate(eGFR_pct = 100 * (eGFR_out - first(eGFR_out)) / first(eGFR_out)) %>%
      ungroup()
    cols <- df %>% distinct(Scenario, Color) %>% deframe()
    p <- ggplot(df, aes(x = Year, y = eGFR_pct, color = Scenario)) +
      geom_line(linewidth = 1.0) +
      geom_hline(yintercept = -40, linetype = "dashed", color = "gray60") +
      scale_color_manual(values = cols) +
      labs(x = "Year", y = "eGFR Change from Baseline (%)",
           title = "eGFR Change (%, threshold -40%)") +
      theme_ckd()
    ggplotly(p)
  })

  output$endpoint_uacr_pct <- renderPlotly({
    df <- multi_results() %>%
      group_by(Scenario, Color) %>%
      mutate(UACR_pct = 100 * (UACR_out - first(UACR_out)) / first(UACR_out)) %>%
      ungroup()
    cols <- df %>% distinct(Scenario, Color) %>% deframe()
    p <- ggplot(df, aes(x = Year, y = UACR_pct, color = Scenario)) +
      geom_line(linewidth = 1.0) +
      scale_color_manual(values = cols) +
      labs(x = "Year", y = "UACR Change from Baseline (%)", title = "UACR % Change") +
      theme_ckd()
    ggplotly(p)
  })

  output$ckd_stage_plot <- renderPlotly({
    df <- multi_results()
    cols <- df %>% distinct(Scenario, Color) %>% deframe()
    p <- ggplot(df, aes(x = Year, y = CKD_stage, color = Scenario)) +
      geom_step(linewidth = 1.0) +
      scale_color_manual(values = cols) +
      scale_y_continuous(breaks = 1:5, labels = paste0("G", 1:5)) +
      labs(x = "Year", y = "CKD Stage", title = "CKD Stage Progression") +
      theme_ckd()
    ggplotly(p)
  })

  output$cv_risk_plot <- renderPlotly({
    df <- multi_results()
    cols <- df %>% distinct(Scenario, Color) %>% deframe()
    p <- ggplot(df, aes(x = Year, y = CV_risk, color = Scenario)) +
      geom_line(linewidth = 1.0) +
      scale_color_manual(values = cols) +
      labs(x = "Year", y = "Composite CV Risk Score", title = "CV Risk Composite") +
      theme_ckd()
    ggplotly(p)
  })

  # ── CV & MBD PLOTS ──────────────────────────────────────────
  output$cv_bp <- renderPlotly({
    df <- d()
    p <- ggplot(df, aes(x = Year, y = SBP_out)) +
      geom_line(color = "#E74C3C", linewidth = 1.2) +
      geom_hline(yintercept = 130, linetype = "dashed", color = "#27AE60") +
      labs(x = "Year", y = "SBP (mmHg)", title = "Blood Pressure") +
      theme_ckd()
    ggplotly(p)
  })

  output$cv_lvh <- renderPlotly({
    df <- d()
    p <- ggplot(df, aes(x = Year, y = LVH_out)) +
      geom_line(color = "#C0392B", linewidth = 1.2) +
      labs(x = "Year", y = "LVH Index (normalized)", title = "Left Ventricular Hypertrophy") +
      theme_ckd()
    ggplotly(p)
  })

  output$cv_vc <- renderPlotly({
    df <- d()
    p <- ggplot(df, aes(x = Year, y = VC_idx)) +
      geom_line(color = "#922B21", linewidth = 1.2) +
      geom_area(fill = "#922B21", alpha = 0.2) +
      labs(x = "Year", y = "VC Index (0–1)", title = "Vascular Calcification") +
      ylim(0, 1) +
      theme_ckd()
    ggplotly(p)
  })

  output$cv_pth <- renderPlotly({
    df <- d()
    p <- ggplot(df, aes(x = Year, y = PTH_out)) +
      geom_line(color = "#E74C3C", linewidth = 1.2) +
      geom_hline(yintercept = 150, linetype = "dashed", color = "#F39C12") +
      geom_hline(yintercept = 300, linetype = "dashed", color = "#E74C3C") +
      labs(x = "Year", y = "Intact PTH (pg/mL)", title = "Secondary Hyperparathyroidism") +
      theme_ckd()
    ggplotly(p)
  })

  output$mbd_fgf23_klotho <- renderPlotly({
    df <- d() %>%
      select(Year, `FGF-23 (pg/mL ÷ 2)` = FGF23_out, `Klotho × 100` = Klotho) %>%
      mutate(`FGF-23 (pg/mL ÷ 2)` = `FGF-23 (pg/mL ÷ 2)` / 2,
             `Klotho × 100` = `Klotho × 100` * 100) %>%
      pivot_longer(-Year)
    p <- ggplot(df, aes(x = Year, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("#E74C3C","#27AE60")) +
      labs(x = "Year", y = "Scaled Value", title = "FGF-23 vs. Klotho Axis") +
      theme_ckd()
    ggplotly(p)
  })

  output$mbd_vitd <- renderPlotly({
    df <- d()
    p <- ggplot(df, aes(x = Year, y = VitD)) +
      geom_line(color = "#F39C12", linewidth = 1.2) +
      geom_hline(yintercept = 20, linetype = "dashed", color = "#E74C3C") +
      labs(x = "Year", y = "Calcitriol (pg/mL)", title = "Active Vitamin D (1,25(OH)₂D₃)") +
      theme_ckd()
    ggplotly(p)
  })

  # Auto-run on startup
  observeEvent(input$run_sim, {}, ignoreNULL = TRUE, ignoreInit = FALSE)
}

# ─────────────────────────────────────────────────────────────
# LAUNCH
# ─────────────────────────────────────────────────────────────
shinyApp(ui, server)
