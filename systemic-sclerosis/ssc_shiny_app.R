## ============================================================
## Systemic Sclerosis (SSc) QSP — Interactive Shiny Application
## ============================================================
## Tabs:
##  1. Patient Profile & Disease Subtype
##  2. Drug PK (Nintedanib, Tocilizumab, MMF, Bosentan, Iloprost)
##  3. Fibrosis Biomarkers (TGF-β, IL-6, mRSS, ECM, collagen)
##  4. Pulmonary Endpoints (FVC, DLCO, HRCT, ILD severity)
##  5. Vascular Endpoints (ET-1, PVR, mPAP, 6MWD, WHO-FC)
##  6. Scenario Comparison & Virtual Population
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(tidyr)
library(DT)
library(shinyWidgets)
library(scales)

## ── Inline mrgsolve model (same as in ssc_mrgsolve_model.R) ──
ssc_code <- '
$PARAM
  BW=65, F_NINT=0.047, KA_NINT=1.2, CL_NINT=1050, V1_NINT=985,
  Q_NINT=420, V2_NINT=2100,
  CL_TCZ=0.22, V1_TCZ=3.7, Q_TCZ=0.09, V2_TCZ=3.0,
  F_MMF=0.94, KA_MMF=1.5, CL_MMF=15.6, V_MMF=30.0,
  F_BOSEN=0.50, KA_BOSEN=0.8, CL_BOSEN=15.5, V_BOSEN=78.0,
  F_ILOP=0.17, KA_ILOP=5.0, CL_ILOP=75.0, V_ILOP=20.0,
  kprod_TGFB=0.05, kdeg_TGFB=0.08, TGFB_base=0.625,
  kprod_IL6=0.20, kdeg_IL6=0.25, IL6_base=0.80,
  kprod_Th17=0.004, kdeg_Th17=0.010, Th17_base=0.40,
  IL6_EC50_Th17=1.5,
  kprod_Bnv=0.008, kdeg_Bnv=0.005, Bnv_base=1.6,
  BAFF_conc=1.0,
  kact_FAct=0.015, kdeg_FAct=0.008, FAct_base=0.30,
  kform_Myo=0.012, kdeg_Myo=0.006, Myo_base=0.15,
  ksynth_Col1=0.020, kdeg_Col1=0.003, Col1_base=6.67,
  ksynth_Col3=0.012, kdeg_Col3=0.003, Col3_base=4.0,
  ECM_0=1.0, kECM=0.005, kdegECM=0.004,
  mRSS_0=20.0, kmRSS_up=0.0020, kmRSS_dn=0.0008, mRSS_max=51.0,
  FVC_0=78.0, kFVC_loss=0.0012, kFVC_rec=0.0003, FVC_min=30.0,
  DLCO_0=62.0, kDLCO_loss=0.0015,
  ET1_0=3.5, kprod_ET1=0.015, kdeg_ET1=0.04,
  NO_0=0.50, kprod_NO=0.02, kdeg_NO=0.04,
  PGI2_0=0.30, kprod_PGI2=0.015, kdeg_PGI2=0.05,
  Endo_0=0.55, kdmg_Endo=0.003, krep_Endo=0.002,
  PVR_0=350.0, kPVR_up=0.0015, kPVR_dn=0.0005,
  mPAP_0=30.0,
  SixMWD_0=380.0, kSixMWD_loss=0.0010, SixMWD_min=100.0,
  EMAX_NINT_FAct=0.60, EC50_NINT_FAct=85.0,
  EMAX_NINT_FVC=0.55, EC50_NINT_FVC=80.0,
  EMAX_TCZ_IL6=0.85, EC50_TCZ_IL6=2.5,
  EMAX_TCZ_mRSS=0.40, EC50_TCZ_mRSS=3.0,
  EMAX_MMF_Th17=0.50, EC50_MMF_Th17=1.8,
  EMAX_BOSEN_ET1=0.75, EC50_BOSEN_ET1=1.2,
  EMAX_BOSEN_PVR=0.60, EC50_BOSEN_PVR=1.5,
  EMAX_ILOP_PGI2=0.80, EC50_ILOP_PGI2=0.5,
  EMAX_ILOP_PVR=0.50, EC50_ILOP_PVR=0.6

$INIT
  A_NINT1=0, A_NINT2=0, A_NINT_GUT=0,
  A_TCZ1=0, A_TCZ2=0,
  A_MMF_GUT=0, A_MPA=0,
  A_BOSEN_GUT=0, A_BOSEN=0,
  A_ILOP_LUNG=0, A_ILOP=0,
  TGFB=0.625, IL6=0.80, Th17=0.40, Bnv=1.60,
  FAct=0.30, Myo=0.15, Col1=6.67, Col3=4.00, ECM=1.00,
  mRSS=20.0, FVC=78.0, DLCO=62.0,
  ET1=3.50, NO=0.50, PGI2_c=0.30, Endo=0.55,
  PVR=350.0, SixMWD=380.0

$ODE
  double C_NINT  = A_NINT1  / V1_NINT  * 1000;
  double C_TCZ   = A_TCZ1   / V1_TCZ;
  double C_MPA   = A_MPA    / V_MMF;
  double C_BOSEN = A_BOSEN  / V_BOSEN;
  double C_ILOP  = A_ILOP   / V_ILOP;

  dxdt_A_NINT_GUT = -KA_NINT * A_NINT_GUT;
  dxdt_A_NINT1    =  KA_NINT * F_NINT * A_NINT_GUT
                   - (CL_NINT + Q_NINT) / V1_NINT * A_NINT1
                   + Q_NINT / V2_NINT * A_NINT2;
  dxdt_A_NINT2    =  Q_NINT / V1_NINT * A_NINT1
                   - Q_NINT / V2_NINT * A_NINT2;

  dxdt_A_TCZ1 = -(CL_TCZ + Q_TCZ) * A_TCZ1 / V1_TCZ
                + Q_TCZ * A_TCZ2 / V2_TCZ;
  dxdt_A_TCZ2 =  Q_TCZ * A_TCZ1 / V1_TCZ
               - Q_TCZ * A_TCZ2 / V2_TCZ;

  dxdt_A_MMF_GUT = -KA_MMF * A_MMF_GUT;
  dxdt_A_MPA     =  KA_MMF * F_MMF * A_MMF_GUT - CL_MMF * C_MPA;

  dxdt_A_BOSEN_GUT = -KA_BOSEN * A_BOSEN_GUT;
  dxdt_A_BOSEN     =  KA_BOSEN * F_BOSEN * A_BOSEN_GUT
                    - CL_BOSEN * C_BOSEN;

  dxdt_A_ILOP_LUNG = -KA_ILOP * A_ILOP_LUNG;
  dxdt_A_ILOP      =  KA_ILOP * F_ILOP * A_ILOP_LUNG
                    - CL_ILOP * C_ILOP;

  double E_NINT_FAct  = EMAX_NINT_FAct * C_NINT  / (EC50_NINT_FAct + C_NINT);
  double E_NINT_FVC   = EMAX_NINT_FVC  * C_NINT  / (EC50_NINT_FVC  + C_NINT);
  double E_TCZ_IL6    = EMAX_TCZ_IL6   * C_TCZ   / (EC50_TCZ_IL6   + C_TCZ);
  double E_TCZ_mRSS   = EMAX_TCZ_mRSS  * C_TCZ   / (EC50_TCZ_mRSS  + C_TCZ);
  double E_MMF_Th17   = EMAX_MMF_Th17  * C_MPA   / (EC50_MMF_Th17  + C_MPA);
  double E_BOSEN_ET1  = EMAX_BOSEN_ET1 * C_BOSEN / (EC50_BOSEN_ET1 + C_BOSEN);
  double E_BOSEN_PVR  = EMAX_BOSEN_PVR * C_BOSEN / (EC50_BOSEN_PVR + C_BOSEN);
  double E_ILOP_PGI2  = EMAX_ILOP_PGI2 * C_ILOP  / (EC50_ILOP_PGI2 + C_ILOP);
  double E_ILOP_PVR   = EMAX_ILOP_PVR  * C_ILOP  / (EC50_ILOP_PVR  + C_ILOP);

  double TGFB_stim = 1.0 + 0.3 * (ECM - 1.0);
  dxdt_TGFB = kprod_TGFB * TGFB_stim * (1 - 0.2 * E_NINT_FAct) - kdeg_TGFB * TGFB;

  double IL6_eff  = IL6 * (1 - E_TCZ_IL6);
  double IL6_stim = 1.0 + 0.5 * (TGFB / TGFB_base - 1.0);
  dxdt_IL6 = kprod_IL6 * IL6_stim - kdeg_IL6 * IL6;

  double Th17_stim = IL6_eff / (IL6_EC50_Th17 + IL6_eff) + 0.3 * TGFB / TGFB_base;
  dxdt_Th17 = kprod_Th17 * Th17_stim * (1 - E_MMF_Th17) - kdeg_Th17 * Th17;

  dxdt_Bnv = kprod_Bnv * BAFF_conc * (1 - 0.3 * E_MMF_Th17) - kdeg_Bnv * Bnv;

  double FAct_stim = (TGFB / TGFB_base) * (1 + 0.2 * (Th17 / 0.4 - 1.0));
  dxdt_FAct = kact_FAct * FAct_stim * (1 - E_NINT_FAct) - kdeg_FAct * FAct;

  dxdt_Myo  = kform_Myo * FAct - kdeg_Myo * Myo;
  double col_stim = Myo / Myo_base;
  dxdt_Col1 = ksynth_Col1 * col_stim - kdeg_Col1 * Col1;
  dxdt_Col3 = ksynth_Col3 * col_stim - kdeg_Col3 * Col3;
  dxdt_ECM  = kECM * (Col1 / Col1_base + Col3 / Col3_base) / 2.0 - kdegECM * ECM;

  double mRSS_prog_rate = kmRSS_up * (ECM - 1.0) * (1 - E_TCZ_mRSS);
  dxdt_mRSS = (mRSS_prog_rate - kmRSS_dn) * mRSS * (1 - mRSS / mRSS_max);
  if (mRSS < 0) dxdt_mRSS = 0;

  double FVC_loss_rate = kFVC_loss * (ECM / 1.0) * (1 - E_NINT_FVC);
  dxdt_FVC  = -(FVC_loss_rate - kFVC_rec) * (FVC - FVC_min);
  dxdt_DLCO = -kDLCO_loss * (ECM / 1.0 + 0.5 * (PVR / PVR_0 - 1.0)) * (DLCO - 30.0);

  double Endo_dmg = kdmg_Endo * (TGFB / TGFB_base + Th17 / 0.4) / 2.0;
  double Endo_rep = krep_Endo * (NO / NO_0 + PGI2_c / PGI2_0) / 2.0
                  * (1 + E_ILOP_PGI2 * 0.5);
  dxdt_Endo = Endo_rep * (1 - Endo) - Endo_dmg * Endo;

  double ET1_prod = kprod_ET1 * (1 / (Endo + 0.01)) * (1 - E_BOSEN_ET1);
  dxdt_ET1 = ET1_prod - kdeg_ET1 * ET1;

  double NO_prod  = kprod_NO  * Endo * (1 + 0.3 * E_ILOP_PGI2);
  dxdt_NO = NO_prod - kdeg_NO * NO;

  double PGI2_prod = kprod_PGI2 * Endo * (1 + E_ILOP_PGI2);
  dxdt_PGI2_c = PGI2_prod - kdeg_PGI2 * PGI2_c;

  double PVR_ET1_drive = ET1 / ET1_0;
  double PVR_NO_brake  = NO   / NO_0;
  double PVR_PGI_brake = PGI2_c / PGI2_0;
  double PVR_prog = kPVR_up * PVR_ET1_drive * (1 - E_BOSEN_PVR) * (ECM / 1.0);
  double PVR_reg  = kPVR_dn * (PVR_NO_brake + PVR_PGI_brake) / 2.0 * (1 + E_ILOP_PVR);
  dxdt_PVR = (PVR_prog - PVR_reg) * PVR * 0.10;

  double SixMWD_loss = kSixMWD_loss * (PVR / PVR_0)
                     * (1 - E_ILOP_PVR * 0.5 - E_BOSEN_PVR * 0.3);
  dxdt_SixMWD = -SixMWD_loss * (SixMWD - SixMWD_min);

$TABLE
  double C_NINT_obs  = A_NINT1 / V1_NINT * 1000;
  double C_TCZ_obs   = A_TCZ1  / V1_TCZ;
  double C_MPA_obs   = A_MPA   / V_MMF;
  double C_BOSEN_obs = A_BOSEN / V_BOSEN;
  double C_ILOP_obs  = A_ILOP  / V_ILOP;
  double mPAP_est    = 0.61 * PVR / 80.0 + 2.0;
  double WHO_FC_est  = (SixMWD > 440) ? 1 :
                       (SixMWD > 315) ? 2 :
                       (SixMWD > 165) ? 3 : 4;

$CAPTURE
  C_NINT_obs C_TCZ_obs C_MPA_obs C_BOSEN_obs C_ILOP_obs
  TGFB IL6 Th17 Bnv FAct Myo Col1 Col3 ECM
  mRSS FVC DLCO ET1 NO PGI2_c Endo PVR SixMWD mPAP_est WHO_FC_est
'

ssc_model <- mcode("ssc_shiny", ssc_code, quiet = TRUE)

## ── Build dosing events ────────────────────────────────────
make_events <- function(use_nint, dose_nint,
                        use_tcz,  dose_tcz,
                        use_mmf,  dose_mmf,
                        use_bosen, dose_bosen,
                        use_ilop,  dose_ilop,
                        dur_yr) {
  dur_h <- dur_yr * 8760
  n_doses <- function(ii) floor(dur_h / ii)

  evs <- list()
  if (use_nint)
    evs[["nint"]] <- ev(cmt = "A_NINT_GUT", amt = dose_nint,
                        ii = 12, addl = n_doses(12) - 1)
  if (use_tcz)
    evs[["tcz"]]  <- ev(cmt = "A_TCZ1", amt = dose_tcz * 1000,
                        ii = 672, addl = n_doses(672) - 1)
  if (use_mmf)
    evs[["mmf"]]  <- ev(cmt = "A_MMF_GUT", amt = dose_mmf,
                        ii = 12, addl = n_doses(12) - 1)
  if (use_bosen)
    evs[["bosen"]] <- ev(cmt = "A_BOSEN_GUT", amt = dose_bosen,
                         ii = 12, addl = n_doses(12) - 1)
  if (use_ilop)
    evs[["ilop"]] <- ev(cmt = "A_ILOP_LUNG", amt = dose_ilop,
                        ii = 8, addl = n_doses(8) - 1)

  if (length(evs) == 0) return(ev(time = 0, cmt = 1, amt = 0))
  Reduce(`+`, evs)
}

run_sim_shiny <- function(model, events, fvc_bl, mrss_bl, pvr_bl, dur_yr,
                          kfvc_loss = 0.0012, emax_nint = 0.55) {
  model %>%
    param(kFVC_loss = kfvc_loss, EMAX_NINT_FVC = emax_nint,
          PVR_0 = pvr_bl) %>%
    init(FVC = fvc_bl, mRSS = mrss_bl, PVR = pvr_bl,
         SixMWD = 380) %>%
    ev(events) %>%
    mrgsim(end = dur_yr * 8760, delta = 24) %>%
    as_tibble() %>%
    mutate(time_yr = time / 8760)
}

## ── UI ────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",

  dashboardHeader(
    title = tags$span(
      style = "font-size:14px; font-weight:bold;",
      "Systemic Sclerosis QSP Dashboard"
    ),
    titleWidth = 320
  ),

  dashboardSidebar(
    width = 320,
    sidebarMenu(
      id = "tabs",
      menuItem("Patient Profile",      tabName = "tab_patient",   icon = icon("user")),
      menuItem("Drug PK",              tabName = "tab_pk",        icon = icon("capsules")),
      menuItem("Fibrosis Biomarkers",  tabName = "tab_fibro",     icon = icon("dna")),
      menuItem("Pulmonary Endpoints",  tabName = "tab_lung",      icon = icon("lungs")),
      menuItem("Vascular Endpoints",   tabName = "tab_vasc",      icon = icon("heartbeat")),
      menuItem("Scenario Comparison",  tabName = "tab_compare",   icon = icon("chart-bar"))
    ),

    hr(),
    h5("  Simulation Settings", style = "padding-left:12px; color:#aaa;"),

    sliderInput("dur_yr", "Simulation Duration (years)",
                min = 0.5, max = 5, value = 2, step = 0.5),

    hr(),
    h5("  Patient Baseline", style = "padding-left:12px; color:#aaa;"),
    sliderInput("fvc_bl", "Baseline FVC % predicted",
                min = 40, max = 95, value = 78, step = 1),
    sliderInput("mrss_bl", "Baseline mRSS (0-51)",
                min = 5, max = 45, value = 20, step = 1),
    sliderInput("pvr_bl", "Baseline PVR (dyn·s·cm⁻⁵)",
                min = 150, max = 800, value = 350, step = 25),

    hr(),
    h5("  Drug Treatment", style = "padding-left:12px; color:#aaa;"),

    checkboxInput("use_nint", "Nintedanib", value = FALSE),
    conditionalPanel("input.use_nint",
      sliderInput("dose_nint", "Dose (mg, BID)",
                  min = 50, max = 200, value = 150, step = 50)),

    checkboxInput("use_tcz", "Tocilizumab (IV q4w)", value = FALSE),
    conditionalPanel("input.use_tcz",
      sliderInput("dose_tcz_mgkg", "Dose (mg/kg)",
                  min = 4, max = 12, value = 8, step = 4)),

    checkboxInput("use_mmf", "MMF (Mycophenolate)", value = FALSE),
    conditionalPanel("input.use_mmf",
      sliderInput("dose_mmf", "Dose (mg, BID)",
                  min = 500, max = 2000, value = 1500, step = 250)),

    checkboxInput("use_bosen", "Bosentan (ERA)", value = FALSE),
    conditionalPanel("input.use_bosen",
      sliderInput("dose_bosen", "Dose (mg, BID)",
                  min = 62.5, max = 250, value = 125, step = 62.5)),

    checkboxInput("use_ilop", "Iloprost (inhaled)", value = FALSE),
    conditionalPanel("input.use_ilop",
      sliderInput("dose_ilop", "Dose (μg, 6x/day)",
                  min = 2.5, max = 5, value = 2.5, step = 2.5)),

    hr(),
    actionButton("run_sim", "Run Simulation",
                 class = "btn-primary", style = "margin:10px; width:90%;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .skin-purple .main-header .logo { background-color:#4a235a; }
      .skin-purple .main-header .navbar { background-color:#5b2c6f; }
      .skin-purple .main-sidebar { background-color:#2c1254; }
      .box.box-solid.box-purple { border-top-color:#7d3c98; }
      .value-box .icon { font-size:45px!important; }
      .summary-box { border-left:4px solid #7d3c98; padding:8px 12px;
                     background:#f9f0ff; border-radius:4px; margin-bottom:8px; }
    "))),

    tabItems(

      ## ── TAB 1: Patient Profile ──────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Disease Overview: Systemic Sclerosis (SSc)", width = 12,
              solidHeader = TRUE, status = "purple",
              fluidRow(
                column(6,
                  h4("Disease Classification"),
                  tags$ul(
                    tags$li(strong("Limited cutaneous SSc (lcSSc):"),
                      " skin thickening limited to hands/face; anti-centromere Ab; PAH risk"),
                    tags$li(strong("Diffuse cutaneous SSc (dcSSc):"),
                      " widespread skin/organ involvement; anti-Scl-70; ILD/SRC risk"),
                    tags$li(strong("SSc sine scleroderma:"),
                      " internal organ fibrosis without skin thickening")
                  ),
                  h4("Key Pathophysiology Triad"),
                  tags$ol(
                    tags$li(strong("Autoimmunity:"), " T/B cell dysregulation, autoantibodies"),
                    tags$li(strong("Vasculopathy:"), " Raynaud's, endothelial injury, ET-1↑"),
                    tags$li(strong("Fibrosis:"), " TGF-β, myofibroblasts, ECM accumulation")
                  )
                ),
                column(6,
                  h4("QSP Model Components"),
                  tags$table(class = "table table-condensed table-bordered",
                    tags$thead(tags$tr(
                      tags$th("Layer"), tags$th("Components"), tags$th("Drugs Modeled")
                    )),
                    tags$tbody(
                      tags$tr(tags$td("Immune"), tags$td("TGF-β, IL-6, Th17, B cells"),
                               tags$td("TCZ, MMF, Rituximab")),
                      tags$tr(tags$td("Fibrosis"), tags$td("FAct, Myo, Col1/3, ECM, mRSS"),
                               tags$td("Nintedanib, Pirfenidone")),
                      tags$tr(tags$td("Pulmonary"), tags$td("FVC, DLCO, ILD severity"),
                               tags$td("Nintedanib, TCZ")),
                      tags$tr(tags$td("Vascular"), tags$td("ET-1, PVR, mPAP, 6MWD"),
                               tags$td("Bosentan, Iloprost, Sildenafil")),
                      tags$tr(tags$td("Renal"), tags$td("RAAS, SRC risk"),
                               tags$td("ACEi (captopril)"))
                    )
                  )
                )
              )
          )
        ),
        fluidRow(
          valueBoxOutput("vbox_fvc",    width = 3),
          valueBoxOutput("vbox_mrss",   width = 3),
          valueBoxOutput("vbox_pvr",    width = 3),
          valueBoxOutput("vbox_6mwd",   width = 3)
        ),
        fluidRow(
          box(title = "Simulated Disease Trajectory", width = 12,
              solidHeader = TRUE, status = "purple",
              plotlyOutput("profile_plot", height = "350px"))
        )
      ),

      ## ── TAB 2: Drug PK ──────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Plasma Concentration — Nintedanib (ng/mL)",
              width = 6, solidHeader = TRUE, status = "info",
              plotlyOutput("pk_nint_plot", height = "280px")),
          box(title = "Plasma Concentration — Tocilizumab (μg/mL)",
              width = 6, solidHeader = TRUE, status = "info",
              plotlyOutput("pk_tcz_plot",  height = "280px"))
        ),
        fluidRow(
          box(title = "Plasma Concentration — MPA (μg/mL; MMF active moiety)",
              width = 4, solidHeader = TRUE, status = "warning",
              plotlyOutput("pk_mpa_plot",   height = "250px")),
          box(title = "Plasma Concentration — Bosentan (μg/mL)",
              width = 4, solidHeader = TRUE, status = "warning",
              plotlyOutput("pk_bosen_plot", height = "250px")),
          box(title = "Plasma Concentration — Iloprost (ng/mL)",
              width = 4, solidHeader = TRUE, status = "warning",
              plotlyOutput("pk_ilop_plot",  height = "250px"))
        ),
        fluidRow(
          box(title = "PK Parameters Reference", width = 12,
              solidHeader = TRUE,
              DT::dataTableOutput("pk_params_table"))
        )
      ),

      ## ── TAB 3: Fibrosis Biomarkers ───────────────────────
      tabItem(tabName = "tab_fibro",
        fluidRow(
          box(title = "TGF-β1 (nmol/L) — Master Fibrotic Cytokine",
              width = 6, solidHeader = TRUE, status = "danger",
              plotlyOutput("fibro_tgfb_plot", height = "260px")),
          box(title = "IL-6 (pmol/L) — Pro-inflammatory / Pro-fibrotic",
              width = 6, solidHeader = TRUE, status = "danger",
              plotlyOutput("fibro_il6_plot",  height = "260px"))
        ),
        fluidRow(
          box(title = "ECM Accumulation (normalized)",
              width = 4, solidHeader = TRUE, status = "warning",
              plotlyOutput("fibro_ecm_plot",  height = "240px")),
          box(title = "Myofibroblasts (relative)",
              width = 4, solidHeader = TRUE, status = "warning",
              plotlyOutput("fibro_myo_plot",  height = "240px")),
          box(title = "Collagen I & III (normalized)",
              width = 4, solidHeader = TRUE, status = "warning",
              plotlyOutput("fibro_col_plot",  height = "240px"))
        ),
        fluidRow(
          box(title = "Modified Rodnan Skin Score (mRSS, 0-51)",
              width = 6, solidHeader = TRUE, status = "danger",
              plotlyOutput("fibro_mrss_plot", height = "260px")),
          box(title = "Th17 Cells & B Cells (relative)",
              width = 6, solidHeader = TRUE, status = "info",
              plotlyOutput("fibro_immune_plot", height = "260px"))
        )
      ),

      ## ── TAB 4: Pulmonary Endpoints ───────────────────────
      tabItem(tabName = "tab_lung",
        fluidRow(
          box(title = "FVC % Predicted (key ILD endpoint)",
              width = 6, solidHeader = TRUE, status = "success",
              plotlyOutput("lung_fvc_plot",  height = "280px")),
          box(title = "DLCO % Predicted",
              width = 6, solidHeader = TRUE, status = "success",
              plotlyOutput("lung_dlco_plot", height = "280px"))
        ),
        fluidRow(
          box(title = "Annual FVC Change (% per year)",
              width = 6, solidHeader = TRUE, status = "warning",
              plotlyOutput("lung_fvc_rate_plot", height = "260px")),
          box(title = "FVC Thresholds: Clinical Significance",
              width = 6, solidHeader = TRUE,
              tags$div(class = "summary-box",
                tags$h5("FVC % Thresholds (ATS/ERS)"),
                tags$table(class = "table table-condensed",
                  tags$tr(tags$td(strong("≥80%")), tags$td("Normal/mild restriction")),
                  tags$tr(tags$td(strong("70-79%")), tags$td("Mild ILD — monitor closely")),
                  tags$tr(tags$td(strong("55-69%")), tags$td("Moderate ILD — treatment indicated")),
                  tags$tr(tags$td(strong("<55%")), tags$td("Severe ILD — transplant evaluation"))
                ),
                tags$h5("SENSCIS Trial (Nintedanib, NEJM 2019)"),
                tags$p("• Nintedanib reduced FVC decline by 41% vs placebo"),
                tags$p("• Annual rate: −52.4 mL/yr (nintedanib) vs −93.3 mL/yr (placebo)"),
                tags$p("• Primary endpoint: rate of FVC decline over 52 weeks")
              )
          )
        )
      ),

      ## ── TAB 5: Vascular Endpoints ────────────────────────
      tabItem(tabName = "tab_vasc",
        fluidRow(
          box(title = "Pulmonary Vascular Resistance (dyn·s·cm⁻⁵)",
              width = 6, solidHeader = TRUE, status = "info",
              plotlyOutput("vasc_pvr_plot",  height = "270px")),
          box(title = "6-Minute Walk Distance (meters)",
              width = 6, solidHeader = TRUE, status = "info",
              plotlyOutput("vasc_6mwd_plot", height = "270px"))
        ),
        fluidRow(
          box(title = "Endothelin-1 (pg/mL)",
              width = 4, solidHeader = TRUE, status = "warning",
              plotlyOutput("vasc_et1_plot",  height = "240px")),
          box(title = "Nitric Oxide & PGI2 (relative)",
              width = 4, solidHeader = TRUE, status = "warning",
              plotlyOutput("vasc_no_plot",   height = "240px")),
          box(title = "Endothelial Integrity (0-1)",
              width = 4, solidHeader = TRUE, status = "warning",
              plotlyOutput("vasc_endo_plot", height = "240px"))
        ),
        fluidRow(
          box(title = "Estimated mean PAP & WHO Functional Class",
              width = 6, solidHeader = TRUE, status = "danger",
              plotlyOutput("vasc_mpap_plot", height = "250px")),
          box(title = "PAH Vascular Pathology Notes", width = 6, solidHeader = TRUE,
              tags$div(class = "summary-box",
                tags$h5("SSc-PAH Clinical Context"),
                tags$ul(
                  tags$li("SSc-PAH accounts for ~30% of PAH cases"),
                  tags$li("5-year survival: ~50% (worse than IPAH)"),
                  tags$li("mPAP >20 mmHg + PVR >2 WU = diagnostic criteria (ESC 2022)"),
                  tags$li("Anti-U3-RNP antibodies: high PAH risk marker"),
                  tags$li("Annual echo screening recommended in SSc")
                ),
                tags$h5("ERA Therapy (Bosentan RAPIDS-1)"),
                tags$p("Bosentan reduced new digital ulcers by 48% vs placebo"),
                tags$p("Macitentan showed mortality reduction in SSc-PAH")
              )
          )
        )
      ),

      ## ── TAB 6: Scenario Comparison ───────────────────────
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title = "Multi-Scenario Outcome Comparison", width = 12,
              solidHeader = TRUE, status = "purple",
              fluidRow(
                column(4,
                  checkboxGroupInput("compare_scenarios",
                    "Select Scenarios to Compare:",
                    choices = c(
                      "Untreated",
                      "Nintedanib 150 mg BID",
                      "Tocilizumab 8 mg/kg q4w",
                      "MMF 1500 mg BID",
                      "Bosentan 125 mg BID",
                      "Nintedanib + Tocilizumab",
                      "Full Combination (Nint+TCZ+MMF+Bosen)"
                    ),
                    selected = c("Untreated", "Nintedanib 150 mg BID",
                                 "Nintedanib + Tocilizumab",
                                 "Full Combination (Nint+TCZ+MMF+Bosen)")
                  )
                ),
                column(8,
                  selectInput("compare_endpoint", "Primary Endpoint:",
                    choices = c("FVC % predicted" = "FVC",
                                "mRSS"            = "mRSS",
                                "PVR"             = "PVR",
                                "SixMWD"          = "SixMWD",
                                "ET-1 (pg/mL)"    = "ET1",
                                "TGF-β1"          = "TGFB",
                                "IL-6"            = "IL6",
                                "ECM"             = "ECM"),
                    selected = "FVC"
                  ),
                  plotlyOutput("compare_main_plot", height = "340px")
                )
              )
          )
        ),
        fluidRow(
          box(title = "Outcomes at Simulation End", width = 7,
              solidHeader = TRUE, status = "primary",
              DT::dataTableOutput("compare_summary_table")),
          box(title = "Relative Treatment Effect vs Untreated", width = 5,
              solidHeader = TRUE, status = "primary",
              plotlyOutput("compare_waterfall", height = "300px"))
        )
      )
    )
  )
)

## ── SERVER ────────────────────────────────────────────────
server <- function(input, output, session) {

  ## ── Run simulation (reactive) ───────────────────────────
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running SSc QSP simulation...", {

      events <- make_events(
        input$use_nint, input$dose_nint,
        input$use_tcz,  input$dose_tcz_mgkg * 65,
        input$use_mmf,  input$dose_mmf,
        input$use_bosen, input$dose_bosen,
        input$use_ilop, if (input$use_ilop) input$dose_ilop else 0,
        input$dur_yr
      )

      run_sim_shiny(
        ssc_model, events,
        fvc_bl  = input$fvc_bl,
        mrss_bl = input$mrss_bl,
        pvr_bl  = input$pvr_bl,
        dur_yr  = input$dur_yr
      )
    })
  }, ignoreNULL = FALSE)

  ## Run default simulation on load
  default_data <- reactive({
    run_sim_shiny(ssc_model, ev(time = 0, cmt = 1, amt = 0),
                  78, 20, 350, 2)
  })

  get_data <- reactive({
    tryCatch(sim_data(), error = function(e) default_data())
  })

  ## ── Value boxes ─────────────────────────────────────────
  output$vbox_fvc <- renderValueBox({
    d <- get_data()
    last_fvc <- tail(d$FVC, 1)
    color <- if (last_fvc >= 70) "green" else if (last_fvc >= 55) "yellow" else "red"
    valueBox(sprintf("%.1f%%", last_fvc), "FVC % Predicted",
             icon = icon("lungs"), color = color)
  })
  output$vbox_mrss <- renderValueBox({
    d <- get_data()
    last_mrss <- tail(d$mRSS, 1)
    color <- if (last_mrss < 14) "green" else if (last_mrss < 22) "yellow" else "red"
    valueBox(sprintf("%.1f", last_mrss), "mRSS (skin score)",
             icon = icon("hand"), color = color)
  })
  output$vbox_pvr <- renderValueBox({
    d <- get_data()
    last_pvr <- tail(d$PVR, 1)
    color <- if (last_pvr < 300) "green" else if (last_pvr < 600) "yellow" else "red"
    valueBox(sprintf("%.0f", last_pvr), "PVR (dyn·s·cm⁻⁵)",
             icon = icon("heartbeat"), color = color)
  })
  output$vbox_6mwd <- renderValueBox({
    d <- get_data()
    last_6mwd <- tail(d$SixMWD, 1)
    color <- if (last_6mwd > 440) "green" else if (last_6mwd > 315) "yellow" else "red"
    valueBox(sprintf("%.0f m", last_6mwd), "6-Minute Walk Distance",
             icon = icon("walking"), color = color)
  })

  ## ── Profile overview ────────────────────────────────────
  output$profile_plot <- renderPlotly({
    d <- get_data()
    fig <- plot_ly(d) %>%
      add_trace(x = ~time_yr, y = ~FVC,   name = "FVC %",  type = "scatter",
                mode = "lines", line = list(color = "#27ae60")) %>%
      add_trace(x = ~time_yr, y = ~mRSS,  name = "mRSS",   type = "scatter",
                mode = "lines", line = list(color = "#e74c3c")) %>%
      layout(title = "Key Disease Endpoints Over Time",
             xaxis = list(title = "Time (years)"),
             yaxis = list(title = "Value"),
             legend = list(orientation = "h"))
    fig
  })

  ## ── PK plots ─────────────────────────────────────────────
  make_pk_plot <- function(d, col, title, unit, n_days = 3) {
    pd <- d %>% filter(time <= n_days * 24)
    if (max(pd[[col]], na.rm = TRUE) < 1e-6) {
      return(plot_ly() %>% add_text(x = 0.5, y = 0.5,
        text = "Drug not selected", textposition = "middle center"))
    }
    plot_ly(pd, x = ~time, y = as.formula(paste0("~", col)),
            type = "scatter", mode = "lines",
            line = list(color = "#2980b9")) %>%
      layout(title = title,
             xaxis = list(title = "Time (hours)"),
             yaxis = list(title = unit))
  }

  output$pk_nint_plot  <- renderPlotly({
    make_pk_plot(get_data(), "C_NINT_obs", "Nintedanib PK", "ng/mL")
  })
  output$pk_tcz_plot   <- renderPlotly({
    make_pk_plot(get_data(), "C_TCZ_obs", "Tocilizumab PK", "μg/mL", n_days = 28)
  })
  output$pk_mpa_plot   <- renderPlotly({
    make_pk_plot(get_data(), "C_MPA_obs",   "MPA PK", "μg/mL")
  })
  output$pk_bosen_plot <- renderPlotly({
    make_pk_plot(get_data(), "C_BOSEN_obs", "Bosentan PK", "μg/mL")
  })
  output$pk_ilop_plot  <- renderPlotly({
    make_pk_plot(get_data(), "C_ILOP_obs",  "Iloprost PK", "ng/mL", n_days = 1)
  })

  output$pk_params_table <- DT::renderDataTable({
    data.frame(
      Drug = c("Nintedanib", "Tocilizumab", "MPA (MMF)", "Bosentan", "Iloprost"),
      Route = c("Oral", "IV", "Oral", "Oral", "Inhaled"),
      `F (%)` = c("4.7", "100", "94", "50", "17"),
      `CL (L/h)` = c("1050", "0.22", "15.6", "15.5", "75"),
      `V1 (L)` = c("985", "3.7", "30", "78", "20"),
      `t½ (h)` = c("0.65", "11.6", "12.0", "5.4", "0.2"),
      Reference = c(
        "Dallmann 2020", "Frey 2010",
        "Kiberd 2006", "van Giersbergen 2003",
        "Hoeper PAH")
    )
  }, options = list(dom = "t", pageLength = 10), rownames = FALSE)

  ## ── Fibrosis plots ──────────────────────────────────────
  make_ts_plot <- function(d, col, title, ylab, color = "#e74c3c",
                           hline = NULL, hline_label = NULL) {
    p <- plot_ly(d, x = ~time_yr, y = as.formula(paste0("~", col)),
                 type = "scatter", mode = "lines",
                 line = list(color = color)) %>%
      layout(title = title,
             xaxis = list(title = "Time (years)"),
             yaxis = list(title = ylab))
    if (!is.null(hline)) {
      p <- p %>% add_segments(x = 0, xend = max(d$time_yr),
                              y = hline, yend = hline,
                              line = list(dash = "dash", color = "gray"))
    }
    p
  }

  output$fibro_tgfb_plot <- renderPlotly({
    make_ts_plot(get_data(), "TGFB", "TGF-β1", "nmol/L", "#c0392b")
  })
  output$fibro_il6_plot  <- renderPlotly({
    make_ts_plot(get_data(), "IL6", "IL-6", "pmol/L", "#e74c3c")
  })
  output$fibro_ecm_plot  <- renderPlotly({
    make_ts_plot(get_data(), "ECM", "ECM Accumulation", "normalized", "#8e44ad")
  })
  output$fibro_myo_plot  <- renderPlotly({
    make_ts_plot(get_data(), "Myo", "Myofibroblasts", "relative units", "#d35400")
  })
  output$fibro_col_plot  <- renderPlotly({
    d <- get_data()
    plot_ly(d) %>%
      add_trace(x = ~time_yr, y = ~Col1, name = "Collagen I",
                type = "scatter", mode = "lines",
                line = list(color = "#f39c12")) %>%
      add_trace(x = ~time_yr, y = ~Col3, name = "Collagen III",
                type = "scatter", mode = "lines",
                line = list(color = "#e67e22")) %>%
      layout(title = "Collagen I & III",
             xaxis = list(title = "Time (years)"),
             yaxis = list(title = "normalized"),
             legend = list(orientation = "h"))
  })
  output$fibro_mrss_plot <- renderPlotly({
    d <- get_data()
    plot_ly(d, x = ~time_yr, y = ~mRSS,
            type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(231,76,60,0.2)",
            line = list(color = "#e74c3c")) %>%
      layout(title = "mRSS Over Time",
             xaxis = list(title = "Time (years)"),
             yaxis = list(title = "mRSS (0-51)", range = c(0, 51)))
  })
  output$fibro_immune_plot <- renderPlotly({
    d <- get_data()
    plot_ly(d) %>%
      add_trace(x = ~time_yr, y = ~Th17, name = "Th17 cells",
                type = "scatter", mode = "lines",
                line = list(color = "#3498db")) %>%
      add_trace(x = ~time_yr, y = ~Bnv,  name = "B cells (naive)",
                type = "scatter", mode = "lines",
                line = list(color = "#2ecc71")) %>%
      layout(title = "Immune Cells",
             xaxis = list(title = "Time (years)"),
             yaxis = list(title = "relative units"),
             legend = list(orientation = "h"))
  })

  ## ── Lung plots ──────────────────────────────────────────
  output$lung_fvc_plot <- renderPlotly({
    d <- get_data()
    plot_ly(d, x = ~time_yr, y = ~FVC,
            type = "scatter", mode = "lines",
            line = list(color = "#27ae60")) %>%
      add_segments(x = 0, xend = max(d$time_yr), y = 70, yend = 70,
                   line = list(dash = "dash", color = "#e74c3c"),
                   name = "FVC 70% threshold") %>%
      layout(title = "FVC % Predicted",
             xaxis = list(title = "Time (years)"),
             yaxis = list(title = "FVC %", range = c(30, 100)))
  })
  output$lung_dlco_plot <- renderPlotly({
    make_ts_plot(get_data(), "DLCO", "DLCO % Predicted", "DLCO %", "#1abc9c")
  })
  output$lung_fvc_rate_plot <- renderPlotly({
    d <- get_data() %>%
      mutate(fvc_rate = c(NA, diff(FVC) / diff(time_yr)))
    plot_ly(d %>% filter(!is.na(fvc_rate)),
            x = ~time_yr, y = ~fvc_rate,
            type = "scatter", mode = "lines",
            line = list(color = "#f39c12")) %>%
      add_segments(x = 0, xend = max(d$time_yr), y = 0, yend = 0,
                   line = list(dash = "dash", color = "gray")) %>%
      layout(title = "FVC Annual Change Rate",
             xaxis = list(title = "Time (years)"),
             yaxis = list(title = "ΔFVC / year (%)"))
  })

  ## ── Vascular plots ──────────────────────────────────────
  output$vasc_pvr_plot  <- renderPlotly({
    make_ts_plot(get_data(), "PVR", "PVR", "dyn·s·cm⁻⁵", "#8e44ad")
  })
  output$vasc_6mwd_plot <- renderPlotly({
    d <- get_data()
    plot_ly(d, x = ~time_yr, y = ~SixMWD,
            type = "scatter", mode = "lines",
            line = list(color = "#2980b9")) %>%
      add_segments(x = 0, xend = max(d$time_yr), y = 165, yend = 165,
                   line = list(dash = "dash", color = "#e74c3c"),
                   name = "WHO-FC IV threshold") %>%
      layout(title = "6-Minute Walk Distance",
             xaxis = list(title = "Time (years)"),
             yaxis = list(title = "6MWD (meters)"))
  })
  output$vasc_et1_plot  <- renderPlotly({
    make_ts_plot(get_data(), "ET1", "Endothelin-1", "pg/mL", "#c0392b")
  })
  output$vasc_no_plot   <- renderPlotly({
    d <- get_data()
    plot_ly(d) %>%
      add_trace(x = ~time_yr, y = ~NO,    name = "NO",   type = "scatter",
                mode = "lines", line = list(color = "#3498db")) %>%
      add_trace(x = ~time_yr, y = ~PGI2_c, name = "PGI2", type = "scatter",
                mode = "lines", line = list(color = "#9b59b6")) %>%
      layout(title = "NO & Prostacyclin",
             xaxis = list(title = "Time (years)"),
             yaxis = list(title = "relative units"),
             legend = list(orientation = "h"))
  })
  output$vasc_endo_plot <- renderPlotly({
    make_ts_plot(get_data(), "Endo", "Endothelial Integrity", "0-1", "#1abc9c")
  })
  output$vasc_mpap_plot <- renderPlotly({
    d <- get_data()
    plot_ly(d, x = ~time_yr, y = ~mPAP_est,
            type = "scatter", mode = "lines",
            line = list(color = "#e74c3c")) %>%
      add_segments(x = 0, xend = max(d$time_yr), y = 20, yend = 20,
                   line = list(dash = "dash", color = "gray")) %>%
      layout(title = "Estimated mPAP",
             xaxis = list(title = "Time (years)"),
             yaxis = list(title = "mPAP (mmHg)"))
  })

  ## ── Scenario comparison ─────────────────────────────────
  scenario_sims <- reactive({
    dur_yr <- input$dur_yr
    fvc_bl  <- input$fvc_bl
    mrss_bl <- input$mrss_bl
    pvr_bl  <- input$pvr_bl

    scenario_list <- list(
      "Untreated" = ev(time = 0, cmt = 1, amt = 0),
      "Nintedanib 150 mg BID" =
        ev(cmt = "A_NINT_GUT", amt = 150, ii = 12, addl = 2919),
      "Tocilizumab 8 mg/kg q4w" =
        ev(cmt = "A_TCZ1", amt = 520000, ii = 672, addl = 51),
      "MMF 1500 mg BID" =
        ev(cmt = "A_MMF_GUT", amt = 1500, ii = 12, addl = 2919),
      "Bosentan 125 mg BID" =
        ev(cmt = "A_BOSEN_GUT", amt = 125, ii = 12, addl = 2919),
      "Nintedanib + Tocilizumab" =
        ev(cmt = "A_NINT_GUT", amt = 150, ii = 12, addl = 2919) +
        ev(cmt = "A_TCZ1", amt = 520000, ii = 672, addl = 51),
      "Full Combination (Nint+TCZ+MMF+Bosen)" =
        ev(cmt = "A_NINT_GUT", amt = 150, ii = 12, addl = 2919) +
        ev(cmt = "A_TCZ1", amt = 520000, ii = 672, addl = 51) +
        ev(cmt = "A_MMF_GUT", amt = 1500, ii = 12, addl = 2919) +
        ev(cmt = "A_BOSEN_GUT", amt = 125, ii = 12, addl = 2919)
    )

    selected <- input$compare_scenarios
    if (is.null(selected)) selected <- c("Untreated", "Nintedanib 150 mg BID")

    bind_rows(lapply(selected, function(sc) {
      run_sim_shiny(ssc_model, scenario_list[[sc]],
                   fvc_bl, mrss_bl, pvr_bl, dur_yr) %>%
        mutate(scenario = sc)
    }))
  })

  colors_palette <- c(
    "Untreated"                              = "#e74c3c",
    "Nintedanib 150 mg BID"                  = "#2980b9",
    "Tocilizumab 8 mg/kg q4w"               = "#27ae60",
    "MMF 1500 mg BID"                        = "#f39c12",
    "Bosentan 125 mg BID"                    = "#8e44ad",
    "Nintedanib + Tocilizumab"               = "#1abc9c",
    "Full Combination (Nint+TCZ+MMF+Bosen)" = "#d35400"
  )

  output$compare_main_plot <- renderPlotly({
    d  <- scenario_sims()
    ep <- input$compare_endpoint
    ylab_map <- c(FVC = "FVC %", mRSS = "mRSS (0-51)",
                  PVR = "dyn·s·cm⁻⁵", SixMWD = "meters",
                  ET1 = "pg/mL", TGFB = "nmol/L",
                  IL6 = "pmol/L", ECM = "normalized")

    p <- plot_ly()
    for (sc in unique(d$scenario)) {
      d_sc <- d %>% filter(scenario == sc)
      p <- p %>% add_trace(
        data = d_sc, x = ~time_yr,
        y = as.formula(paste0("~", ep)),
        name = sc, type = "scatter", mode = "lines",
        line = list(color = colors_palette[sc])
      )
    }
    p %>% layout(
      title = paste("Comparison:", ep),
      xaxis = list(title = "Time (years)"),
      yaxis = list(title = ylab_map[ep]),
      legend = list(orientation = "h", y = -0.3)
    )
  })

  output$compare_summary_table <- DT::renderDataTable({
    d <- scenario_sims()
    d %>%
      group_by(scenario) %>%
      filter(time_yr == max(time_yr)) %>%
      summarise(
        `FVC %`   = round(mean(FVC), 1),
        mRSS      = round(mean(mRSS), 1),
        PVR       = round(mean(PVR), 0),
        `6MWD (m)` = round(mean(SixMWD), 0),
        `ET-1`    = round(mean(ET1), 2),
        `mPAP`    = round(mean(mPAP_est), 1),
        .groups = "drop"
      )
  }, options = list(dom = "t"), rownames = FALSE)

  output$compare_waterfall <- renderPlotly({
    d <- scenario_sims()
    baseline_fvc <- d %>%
      filter(scenario == "Untreated") %>%
      filter(time_yr == max(time_yr)) %>%
      pull(FVC) %>% mean()

    wf <- d %>%
      filter(time_yr == max(time_yr)) %>%
      group_by(scenario) %>%
      summarise(FVC_end = mean(FVC), .groups = "drop") %>%
      mutate(FVC_diff = FVC_end - baseline_fvc) %>%
      filter(scenario != "Untreated") %>%
      arrange(desc(FVC_diff))

    plot_ly(wf, x = ~reorder(scenario, FVC_diff), y = ~FVC_diff,
            type = "bar",
            marker = list(color = "#2980b9")) %>%
      layout(title = "FVC Gain vs Untreated",
             xaxis = list(title = "", tickangle = -30),
             yaxis = list(title = "ΔFVC % vs Untreated"),
             showlegend = FALSE)
  })
}

shinyApp(ui = ui, server = server)
