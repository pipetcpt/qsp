## =============================================================================
## Microscopic Polyangiitis (MPA) — Shiny Interactive Dashboard
## QSP Model: 21 ODE compartments, 7 treatment scenarios
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

## ── Inline mrgsolve model (same as mpa_mrgsolve_model.R) ──────────────────
mpa_code <- '
$PROB MPA QSP Model — Shiny Version

$PARAM
  CL_CY=8.6, V1_CY=31, Q_CY=18, V2_CY=38, KA_CY=1.2, F_CY=0.75,
  CLm_CY=3.2, Vm_CY=12, FRAC_OH=0.10, IC50_CY_BC=0.5, IC50_CY_PC=1.2,
  CL_RTX=0.014, V1_RTX=3.6, Q_RTX=0.22, V2_RTX=5.1, KDE_RTX=0.0008,
  KA_PR=3.6, F_PR=0.85, CL_PR=14.5, V_PR=46,
  KA_AVA=0.8, F_AVA=0.45, CL_AVA=4.2, V_AVA=60,
  Ks_BC=0.012, Kd_BC=0.012, Kp_BC=0.08, Kd_PC=0.003,
  KS_ANCA=0.002, KD_ANCA=0.015,
  KS_C5a=0.04, KD_C5a=0.35, KANCA_C5=0.15, IC50_AVA_C5=0.03,
  KN_act=0.08, KN_res=0.06, C5a_EC50=0.5,
  KE_inj=0.05, KE_res=0.03,
  KR_inflam=0.04, KR_res=0.02, KR_fibr=0.002, KD_fibr=0.0003,
  GFR0=90, KG_loss=0.08, GFR_min=5,
  KP_inflam=0.06, KP_res=0.04, DLCO0=100,
  KS_CRP=2.5, KD_CRP=0.06,
  wt_ANCA=0.30, wt_GFR=0.35, wt_PULM=0.20, wt_CRP=0.15,
  EMAX_GC=0.75, EC50_GC=0.08, HILL_GC=1.5,
  EMAX_AZA=0.55, EC50_AZA=150, EMAX_MMF=0.60, EC50_MMF=2.5,
  CY_dose=0, RTX_mg=0, PRED_mg=60, AZA_mg=0, MMF_mg=0,
  PLEX_on=0, AVA_mg=0

$CMT DEPOT_CY CY_C CY_P OH_CY RTX_C RTX_P PRED_DEP PRED_C AVA_DEP AVA_C
     B_CELL PLASMA ANCA C5a PMN_ACT ENDO_INJ RENAL_I RENAL_F GFR_C PULM_I CRP_C

$INIT
  DEPOT_CY=0, CY_C=0, CY_P=0, OH_CY=0, RTX_C=0, RTX_P=0,
  PRED_DEP=0, PRED_C=0, AVA_DEP=0, AVA_C=0,
  B_CELL=1.5, PLASMA=2.5, ANCA=8.0, C5a=2.0,
  PMN_ACT=2.5, ENDO_INJ=0.6,
  RENAL_I=3.0, RENAL_F=0.3, GFR_C=45.0, PULM_I=1.5, CRP_C=85.0

$ODE
  double dose_CY_cont = CY_dose / 24.0;
  dxdt_DEPOT_CY = F_CY * dose_CY_cont - KA_CY * DEPOT_CY;
  dxdt_CY_C = KA_CY * DEPOT_CY / V1_CY - (CL_CY+Q_CY)/V1_CY*CY_C + Q_CY/V2_CY*CY_P;
  dxdt_CY_P = Q_CY/V1_CY*CY_C - Q_CY/V2_CY*CY_P;
  dxdt_OH_CY = FRAC_OH*CL_CY/Vm_CY*CY_C - CLm_CY/Vm_CY*OH_CY;
  dxdt_RTX_C = -(CL_RTX+Q_RTX)/V1_RTX*RTX_C + Q_RTX/V2_RTX*RTX_P - KDE_RTX*B_CELL*RTX_C;
  dxdt_RTX_P = Q_RTX/V1_RTX*RTX_C - Q_RTX/V2_RTX*RTX_P;
  double dose_PR_cont = PRED_mg / 24.0;
  dxdt_PRED_DEP = F_PR*dose_PR_cont - KA_PR*PRED_DEP;
  dxdt_PRED_C = KA_PR*PRED_DEP/V_PR - CL_PR/V_PR*PRED_C;
  double dose_AVA_cont = AVA_mg*2.0/24.0;
  dxdt_AVA_DEP = F_AVA*dose_AVA_cont - KA_AVA*AVA_DEP;
  dxdt_AVA_C = KA_AVA*AVA_DEP/V_AVA - CL_AVA/V_AVA*AVA_C;
  double INH_OH_CY_BC = OH_CY/(OH_CY+IC50_CY_BC);
  double INH_OH_CY_PC = OH_CY/(OH_CY+IC50_CY_PC);
  double INH_GC = EMAX_GC*pow(PRED_C,HILL_GC)/(pow(PRED_C,HILL_GC)+pow(EC50_GC,HILL_GC));
  double INH_AVA = AVA_C/(AVA_C+IC50_AVA_C5);
  double INH_RTX = (RTX_C>0.001) ? RTX_C/(RTX_C+0.1) : 0.0;
  double EFF_AZA = (AZA_mg>0) ? EMAX_AZA*AZA_mg/(AZA_mg+EC50_AZA*24.0/2.0) : 0.0;
  double EFF_MMF = (MMF_mg>0) ? EMAX_MMF*MMF_mg/(MMF_mg*24.0+EC50_MMF*V_PR*24.0) : 0.0;
  double PLEX_kd = (PLEX_on>0.5) ? 0.08 : 0.0;
  double B_death = Kd_BC + INH_OH_CY_BC*0.15 + INH_RTX*0.12 + EFF_AZA*0.06 + EFF_MMF*0.05;
  dxdt_B_CELL = Ks_BC - B_death*B_CELL - Kp_BC*(1-INH_GC*0.3)*B_CELL;
  dxdt_PLASMA = Kp_BC*(1-INH_GC*0.3)*B_CELL - (Kd_PC+INH_OH_CY_PC*0.10+EFF_MMF*0.08)*PLASMA;
  dxdt_ANCA = KS_ANCA*PLASMA - (KD_ANCA+PLEX_kd)*ANCA;
  dxdt_C5a = KS_C5a*(1+KANCA_C5*PMN_ACT) - KD_C5a*(1+INH_AVA)*C5a;
  double ANCA_drive = (ANCA>1.0) ? (ANCA-1.0)/(ANCA+2.0) : 0.0;
  double C5a_drive = C5a/(C5a+C5a_EC50);
  dxdt_PMN_ACT = KN_act*(ANCA_drive+C5a_drive)*(1-INH_GC*0.4)*(1-INH_AVA*0.5) - KN_res*PMN_ACT;
  dxdt_ENDO_INJ = KE_inj*PMN_ACT*(1-INH_GC*0.5) - KE_res*ENDO_INJ;
  double RI_drive = ENDO_INJ*(ANCA>1.0?1.0:0.3)*(1-INH_GC*0.6)*(1-INH_OH_CY_BC*0.3);
  dxdt_RENAL_I = KR_inflam*RI_drive - KR_res*RENAL_I;
  dxdt_RENAL_F = KR_fibr*RENAL_I - KD_fibr*RENAL_F;
  double GFR_loss = KG_loss*RENAL_I*(1-INH_GC*0.3);
  dxdt_GFR_C = (GFR_C>GFR_min) ? -GFR_loss : 0.0;
  dxdt_PULM_I = KP_inflam*PMN_ACT*(1-INH_GC*0.65)*(1-INH_OH_CY_BC*0.25) - KP_res*PULM_I;
  dxdt_CRP_C = KS_CRP*(RENAL_I+PULM_I)*(1-INH_GC*0.8) - KD_CRP*CRP_C;

$TABLE
  double ANCA_titer = ANCA;
  double GFR_obs = GFR_C;
  double Creat_est = (GFR_C>1) ? 8100.0/(GFR_C*1.1) : 80.0;
  double DLCO_obs = DLCO0 - 15.0*PULM_I;
  if(DLCO_obs<20) DLCO_obs=20;
  double BVAS_analog = wt_ANCA*(ANCA/8.0)*10.0 + wt_GFR*((GFR0-GFR_C)/GFR0)*20.0
                     + wt_PULM*PULM_I*5.0 + wt_CRP*log(CRP_C+1.0)*2.0;
  if(BVAS_analog<0) BVAS_analog=0;
  double B_depletion = (RTX_C>0.001) ? RTX_C/(RTX_C+0.2)*0.95 : 0.0;
  double in_remission = (BVAS_analog<2.0&&ANCA<1.5) ? 1.0 : 0.0;

$CAPTURE CY_C OH_CY RTX_C PRED_C AVA_C B_CELL PLASMA ANCA_titer C5a PMN_ACT
         ENDO_INJ RENAL_I RENAL_F GFR_obs Creat_est PULM_I DLCO_obs
         CRP_C BVAS_analog B_depletion in_remission
'

## ── Compile once ────────────────────────────────────────────────────────────
mod <- mcode("MPA_Shiny", mpa_code)

## ── Scenario definitions ────────────────────────────────────────────────────
scenario_defs <- data.frame(
  key       = c("1_Untreated","2_CY_GC","3_RTX_GC","4_CY_PLEX","5_RTX_AVA","6_AZA_Maint","7_RTX_Maint"),
  label     = c(
    "Untreated (natural history)",
    "CY oral + Prednisolone (CYCLOPS)",
    "Rituximab + Prednisolone (RAVE)",
    "CY + Prednisolone + Plasma Exchange (PEXIVAS)",
    "Rituximab + Avacopan — GC-free (ADVOCATE)",
    "AZA Maintenance post-induction (IMPROVE)",
    "Rituximab Maintenance 500mg Q6M (MAINRITSAN)"
  ),
  CY_dose   = c(0,  140,  0, 140,  0,   0,  0),
  RTX_mg    = c(0,    0, 375,   0,375,   0,500),
  PRED_mg   = c(5,   60,  60,  60,  0,  10, 10),
  AZA_mg    = c(0,    0,   0,   0,  0, 150,  0),
  MMF_mg    = c(0,    0,   0,   0,  0,   0,  0),
  PLEX_on   = c(0,    0,   0,   1,  0,   0,  0),
  AVA_mg    = c(0,    0,   0,   0, 30,   0,  0),
  stringsAsFactors = FALSE
)

sc_colors <- c(
  "1_Untreated"    = "#CC0000",
  "2_CY_GC"        = "#FF8800",
  "3_RTX_GC"       = "#0055CC",
  "4_CY_PLEX"      = "#9900CC",
  "5_RTX_AVA"      = "#00AA44",
  "6_AZA_Maint"    = "#888888",
  "7_RTX_Maint"    = "#006699"
)

run_sim <- function(sc_rows, dur_days=90) {
  results <- lapply(seq_len(nrow(sc_rows)), function(i) {
    r <- sc_rows[i,]
    mod %>%
      param(CY_dose=r$CY_dose, RTX_mg=r$RTX_mg, PRED_mg=r$PRED_mg,
            AZA_mg=r$AZA_mg, MMF_mg=r$MMF_mg, PLEX_on=r$PLEX_on,
            AVA_mg=r$AVA_mg) %>%
      mrgsim(end=dur_days*24, delta=6) %>%
      as.data.frame() %>%
      mutate(scenario=r$key, label=r$label, time_days=time/24)
  })
  bind_rows(results)
}

## ============================================================================
## UI
## ============================================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "MPA QSP Dashboard", titleWidth=300),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",      tabName="tab_profile",  icon=icon("user-md")),
      menuItem("Drug PK",              tabName="tab_pk",       icon=icon("pills")),
      menuItem("Immunology & ANCA",    tabName="tab_immuno",   icon=icon("microscope")),
      menuItem("Renal Endpoints",      tabName="tab_renal",    icon=icon("kidney")),
      menuItem("Pulmonary Endpoints",  tabName="tab_pulm",     icon=icon("lungs")),
      menuItem("Scenario Comparison",  tabName="tab_scenario", icon=icon("chart-line")),
      menuItem("References",           tabName="tab_refs",     icon=icon("book"))
    ),
    hr(),
    h4("Simulation Controls", style="color:white; padding-left:15px;"),
    sliderInput("sim_duration", "Duration (days):", min=30, max=365, value=90, step=15),
    sliderInput("init_gfr",  "Baseline GFR:", min=5, max=90, value=45, step=5),
    sliderInput("init_anca", "Baseline ANCA titer:", min=1, max=20, value=8, step=0.5),
    sliderInput("init_bvas", "Initial BVAS severity (1-3=mild, 4=severe):",
                min=1, max=4, value=3, step=0.5),
    hr(),
    checkboxGroupInput("sc_select", "Select Scenarios:",
      choiceNames  = scenario_defs$label,
      choiceValues = scenario_defs$key,
      selected     = c("1_Untreated","2_CY_GC","3_RTX_GC","5_RTX_AVA")
    ),
    actionButton("run_sim", "Run Simulation", icon=icon("play"),
                 class="btn btn-success btn-block", style="margin:10px 15px;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f7f9; }
      .box-title { font-weight: bold; }
    "))),
    tabItems(

      ## ----------------------------------------------------------------
      ## TAB 1: PATIENT PROFILE
      ## ----------------------------------------------------------------
      tabItem(tabName="tab_profile",
        fluidRow(
          box(title="Microscopic Polyangiitis (MPA) — Disease Overview",
              width=12, status="primary", solidHeader=TRUE,
            fluidRow(
              column(6,
                h4("Key Features"),
                tags$ul(
                  tags$li("Small-vessel ANCA vasculitis (arterioles, venules, capillaries)"),
                  tags$li("Dominant ANCA: Anti-MPO IgG (pANCA) in ~75% of cases"),
                  tags$li("Pauci-immune glomerulonephritis (no immune deposits by IF)"),
                  tags$li("Rapidly progressive GN (RPGN) → crescentic GN"),
                  tags$li("Diffuse alveolar hemorrhage (DAH) in 30-40%"),
                  tags$li("No granuloma formation (distinguishes from GPA)"),
                  tags$li("Annual incidence: 2-10 per million (higher in Japan)")
                ),
                h4("Diagnostic Criteria"),
                tags$ul(
                  tags$li("ANCA positivity (MPO-ANCA by ELISA)"),
                  tags$li("Biopsy: pauci-immune necrotizing vasculitis"),
                  tags$li("RPGN with crescent formation on renal biopsy"),
                  tags$li("Pulmonary capillaritis ± DAH"),
                  tags$li("Exclusion of GPA, EGPA, PAN (no large-vessel disease)")
                )
              ),
              column(6,
                h4("Disease Activity Score (BVAS)"),
                tableOutput("bvas_table"),
                h4("Key Biomarkers"),
                tableOutput("biomarker_table")
              )
            )
          )
        ),
        fluidRow(
          box(title="MPA Pathophysiology Pathway",
              width=6, status="warning",
            p("Core disease mechanism:"),
            tags$ol(
              tags$li("Genetic susceptibility (HLA-DQ, PTPN22) + environmental trigger"),
              tags$li("Loss of B cell tolerance → anti-MPO IgG production (pANCA)"),
              tags$li("TNF-α/GM-CSF priming causes surface expression of MPO on neutrophils"),
              tags$li("ANCA binds MPO on primed neutrophils → FcγR activation"),
              tags$li("Neutrophil respiratory burst → ROS, degranulation, NETosis"),
              tags$li("Small vessel wall destruction → pauci-immune vasculitis"),
              tags$li("Glomerular capillaritis → RPGN, crescent GN, GFR↓"),
              tags$li("Pulmonary capillaritis → DAH, hemoptysis, DLCO↓")
            )
          ),
          box(title="Treatment Algorithm",
              width=6, status="success",
            h5("Induction Therapy (0-3 months)"),
            tags$ul(
              tags$li("Rituximab 375mg/m² × 4 doses + Prednisolone OR"),
              tags$li("Cyclophosphamide IV 15mg/kg q3w × 6 + Prednisolone"),
              tags$li("PLEX for severe renal (Cr>500μmol/L) or DAH (controversial)"),
              tags$li("Avacopan 30mg BID may replace prednisolone (ADVOCATE 2021)")
            ),
            h5("Maintenance Therapy (3-24 months)"),
            tags$ul(
              tags$li("Rituximab 500mg q6m (preferred, MAINRITSAN) OR"),
              tags$li("Azathioprine 2mg/kg/day OR"),
              tags$li("Mycophenolate mofetil 3g/day (inferior to RTX)")
            ),
            h5("Prognosis"),
            tags$ul(
              tags$li("5-year survival ~80% with treatment"),
              tags$li("Relapse rate: 25-35% at 5 years (higher than GPA)"),
              tags$li("ESRD risk: 20-30% at 5 years without adequate treatment")
            )
          )
        )
      ),

      ## ----------------------------------------------------------------
      ## TAB 2: DRUG PK
      ## ----------------------------------------------------------------
      tabItem(tabName="tab_pk",
        fluidRow(
          box(title="Cyclophosphamide & 4-OH-CY PK", width=6, status="warning", solidHeader=TRUE,
            plotlyOutput("plot_cy_pk", height="320px")
          ),
          box(title="Rituximab PK (Central Compartment)", width=6, status="primary", solidHeader=TRUE,
            plotlyOutput("plot_rtx_pk", height="320px")
          )
        ),
        fluidRow(
          box(title="Prednisolone PK", width=6, status="danger", solidHeader=TRUE,
            plotlyOutput("plot_pred_pk", height="320px")
          ),
          box(title="Avacopan PK (C5aR1 Inhibitor)", width=6, status="success", solidHeader=TRUE,
            plotlyOutput("plot_ava_pk", height="320px")
          )
        ),
        fluidRow(
          box(title="PK Parameter Summary", width=12, status="info",
            DTOutput("pk_param_table")
          )
        )
      ),

      ## ----------------------------------------------------------------
      ## TAB 3: IMMUNOLOGY & ANCA
      ## ----------------------------------------------------------------
      tabItem(tabName="tab_immuno",
        fluidRow(
          box(title="ANCA Titer (Anti-MPO IgG)", width=6, status="danger", solidHeader=TRUE,
            plotlyOutput("plot_anca", height="320px")
          ),
          box(title="B Cell & Plasma Cell Dynamics", width=6, status="primary", solidHeader=TRUE,
            plotlyOutput("plot_bcell", height="320px")
          )
        ),
        fluidRow(
          box(title="Complement C5a & Neutrophil Activation", width=6, status="warning", solidHeader=TRUE,
            plotlyOutput("plot_c5a_pmn", height="320px")
          ),
          box(title="Endothelial Injury Index", width=6, status="info", solidHeader=TRUE,
            plotlyOutput("plot_endo", height="320px")
          )
        )
      ),

      ## ----------------------------------------------------------------
      ## TAB 4: RENAL ENDPOINTS
      ## ----------------------------------------------------------------
      tabItem(tabName="tab_renal",
        fluidRow(
          box(title="GFR Trajectory", width=6, status="primary", solidHeader=TRUE,
            plotlyOutput("plot_gfr", height="320px")
          ),
          box(title="Estimated Serum Creatinine", width=6, status="warning", solidHeader=TRUE,
            plotlyOutput("plot_creat", height="320px")
          )
        ),
        fluidRow(
          box(title="Renal Inflammation Index", width=6, status="danger", solidHeader=TRUE,
            plotlyOutput("plot_renal_i", height="320px")
          ),
          box(title="Renal Fibrosis Index", width=6, status="info", solidHeader=TRUE,
            plotlyOutput("plot_renal_f", height="320px")
          )
        ),
        fluidRow(
          box(title="Renal Outcome Summary", width=12, status="success",
            DTOutput("renal_summary_tbl")
          )
        )
      ),

      ## ----------------------------------------------------------------
      ## TAB 5: PULMONARY ENDPOINTS
      ## ----------------------------------------------------------------
      tabItem(tabName="tab_pulm",
        fluidRow(
          box(title="Pulmonary Inflammation (DAH Index)", width=6, status="danger", solidHeader=TRUE,
            plotlyOutput("plot_pulm_i", height="320px")
          ),
          box(title="DLCO (% Predicted)", width=6, status="primary", solidHeader=TRUE,
            plotlyOutput("plot_dlco", height="320px")
          )
        ),
        fluidRow(
          box(title="CRP (Systemic Inflammation)", width=6, status="warning", solidHeader=TRUE,
            plotlyOutput("plot_crp", height="320px")
          ),
          box(title="DAH Severity Interpretation", width=6, status="info",
            h4("DAH Clinical Staging (Pulmonary Index)"),
            tableOutput("dah_stage_table"),
            hr(),
            h4("Pulmonary Management"),
            tags$ul(
              tags$li("Mild DAH (index 1-2): high-dose steroids"),
              tags$li("Moderate DAH (index 2-4): add CY or RTX"),
              tags$li("Severe DAH (index >4): consider PLEX + mechanical ventilation"),
              tags$li("Recurrent DAH: maintenance RTX preferred over AZA")
            )
          )
        )
      ),

      ## ----------------------------------------------------------------
      ## TAB 6: SCENARIO COMPARISON
      ## ----------------------------------------------------------------
      tabItem(tabName="tab_scenario",
        fluidRow(
          box(title="BVAS Disease Activity Score", width=6, status="primary", solidHeader=TRUE,
            plotlyOutput("comp_bvas", height="320px")
          ),
          box(title="GFR Comparison", width=6, status="success", solidHeader=TRUE,
            plotlyOutput("comp_gfr", height="320px")
          )
        ),
        fluidRow(
          box(title="ANCA Titer Comparison", width=6, status="warning", solidHeader=TRUE,
            plotlyOutput("comp_anca", height="320px")
          ),
          box(title="DLCO Comparison", width=6, status="info", solidHeader=TRUE,
            plotlyOutput("comp_dlco", height="320px")
          )
        ),
        fluidRow(
          box(title="Outcome Summary at End of Simulation", width=12, status="danger",
            DTOutput("outcome_summary_tbl")
          )
        )
      ),

      ## ----------------------------------------------------------------
      ## TAB 7: REFERENCES
      ## ----------------------------------------------------------------
      tabItem(tabName="tab_refs",
        fluidRow(
          box(title="Key Clinical Trials", width=6, status="primary", solidHeader=TRUE,
            tags$ul(
              tags$li(tags$b("RAVE (2010):"), " RTX vs CY for ANCA vasculitis induction — PMID 20647198"),
              tags$li(tags$b("RITUXVAS (2010):"), " RTX + CY IV vs CY IV for renal disease — PMID 20647199"),
              tags$li(tags$b("MYCYC (2019):"), " MMF vs CY for induction — PMID 31373547"),
              tags$li(tags$b("IMPROVE (2010):"), " MMF vs AZA for maintenance — PMID 20647200"),
              tags$li(tags$b("MAINRITSAN (2014):"), " RTX vs AZA maintenance — PMID 25372085"),
              tags$li(tags$b("MAINRITSAN2 (2018):"), " Tailored vs fixed RTX maintenance — PMID 30396304"),
              tags$li(tags$b("PEXIVAS (2020):"), " Plasma exchange in ANCA vasculitis — PMID 32053298"),
              tags$li(tags$b("ADVOCATE (2021):"), " Avacopan vs prednisone — PMID 33596356")
            )
          ),
          box(title="Mechanistic References", width=6, status="warning", solidHeader=TRUE,
            tags$ul(
              tags$li("Jennette & Falk. N Engl J Med 2013; PMID: 23984731 (ANCA disease mechanisms)"),
              tags$li("Specks et al. Arthritis Rheum 2013; PMID: 23666854 (RTX remission biomarkers)"),
              tags$li("Mukhtyar et al. Ann Rheum Dis 2009; PMID: 18413398 (BVAS scoring)"),
              tags$li("Schreiber et al. J Clin Invest 2016 (NET in ANCA vasculitis)"),
              tags$li("Hu et al. Kidney Int 2020 (complement in ANCA vasculitis)")
            )
          )
        ),
        fluidRow(
          box(title="About This Model", width=12, status="info",
            p("This QSP model integrates 21 ODE compartments covering drug PK (CY, RTX, PRED, avacopan),",
              "B cell/ANCA dynamics, complement (C5a), neutrophil-mediated vasculitis, renal pauci-immune GN,",
              "and pulmonary DAH to simulate MPA disease progression and treatment response."),
            p("Calibration targets: RAVE (PMID 20647198), PEXIVAS (PMID 32053298), ADVOCATE (PMID 33596356)."),
            p("Model version 1.0 | Generated by QSP Library CCR | 2026-06-20")
          )
        )
      )

    ) # end tabItems
  )
)

## ============================================================================
## SERVER
## ============================================================================
server <- function(input, output, session) {

  ## ── Simulation reactive ──────────────────────────────────────────────────
  sim_data <- eventReactive(input$run_sim, {
    sc_sel <- scenario_defs %>% filter(key %in% input$sc_select)
    if (nrow(sc_sel) == 0) return(NULL)

    # Update initial conditions based on sliders
    mod2 <- mod %>%
      init(GFR_C = input$init_gfr,
           ANCA  = input$init_anca,
           RENAL_I = input$init_bvas * 0.8,
           PULM_I  = input$init_bvas * 0.4)

    results <- lapply(seq_len(nrow(sc_sel)), function(i) {
      r <- sc_sel[i,]
      mod2 %>%
        param(CY_dose=r$CY_dose, RTX_mg=r$RTX_mg, PRED_mg=r$PRED_mg,
              AZA_mg=r$AZA_mg, MMF_mg=r$MMF_mg, PLEX_on=r$PLEX_on,
              AVA_mg=r$AVA_mg) %>%
        mrgsim(end=input$sim_duration*24, delta=6) %>%
        as.data.frame() %>%
        mutate(scenario=r$key, label=r$label, time_days=time/24)
    })
    bind_rows(results)
  }, ignoreNULL=FALSE)

  ## Initialize with default run
  observe({
    isolate({ input$run_sim })
  })

  ## ── Helper: plotly line chart ────────────────────────────────────────────
  make_plotly <- function(df, y_col, y_label, title, hline=NULL, hline_lab=NULL) {
    if (is.null(df) || nrow(df)==0) return(plotly_empty())
    df$color_key <- df$scenario
    p <- plot_ly(df, x=~time_days, y=~.data[[y_col]], color=~label,
                 colors=unname(sc_colors[unique(df$scenario)]),
                 type="scatter", mode="lines",
                 line=list(width=2)) %>%
      layout(title=list(text=title, font=list(size=13)),
             xaxis=list(title="Time (days)"),
             yaxis=list(title=y_label),
             legend=list(orientation="h", y=-0.3, font=list(size=9)))
    if (!is.null(hline)) {
      p <- p %>% add_lines(x=range(df$time_days),
                           y=c(hline,hline), showlegend=FALSE,
                           line=list(dash="dash", color="gray50", width=1))
    }
    p
  }

  ## ── Tab 1 tables ─────────────────────────────────────────────────────────
  output$bvas_table <- renderTable({
    data.frame(
      Component     = c("Renal","General","CNS","Cardiovascular","Pulmonary","ENT","Skin","Eyes"),
      `Max Score`   = c(12, 12, 9, 6, 6, 6, 6, 6),
      `MPA Typical` = c("High","Moderate","Low","Low","Moderate","Low","Low","Low")
    )
  }, striped=TRUE, bordered=TRUE, hover=TRUE)

  output$biomarker_table <- renderTable({
    data.frame(
      Biomarker   = c("MPO-ANCA","CRP","ESR","Creatinine","Hematuria","Proteinuria","DLCO","C5a"),
      `Normal`    = c("<3.5 EU","<10 mg/L","<20 mm/h","<110 μmol/L","0 RBC/HPF","<0.3 g/d","75-140%","<0.1 relative"),
      `Active MPA`= c(">100 EU",">50 mg/L",">50 mm/h",">200 μmol/L",">10 RBC/HPF",">3 g/d","<60%",">0.5 relative")
    )
  }, striped=TRUE, bordered=TRUE)

  ## ── PK plots (Tab 2) ─────────────────────────────────────────────────────
  output$plot_cy_pk <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "OH_CY", "4-OH-CY (mg/L)", "Cyclophosphamide Active Metabolite")
  })
  output$plot_rtx_pk <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "RTX_C", "RTX Concentration (mg/L)", "Rituximab PK")
  })
  output$plot_pred_pk <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "PRED_C", "Prednisolone (mg/L)", "Prednisolone PK", hline=0.08)
  })
  output$plot_ava_pk <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "AVA_C", "Avacopan (mg/L)", "Avacopan PK", hline=0.03)
  })
  output$pk_param_table <- renderDT({
    datatable(data.frame(
      Drug          = c("Cyclophosphamide","4-OH-CY","Rituximab","Prednisolone","Avacopan"),
      Route         = c("Oral/IV","Metabolite","IV","Oral","Oral"),
      Dose          = c("2mg/kg/d or 15mg/kg IV","–","375mg/m²×4","1mg/kg→taper","30mg BID"),
      `t1/2 (h)`    = c("3-10","~4","~350h","2-4","~16"),
      `Vd (L)`      = c("31","12","3.6","46","60"),
      `CL (L/h)`    = c("8.6","3.2","0.014","14.5","4.2"),
      Target        = c("DNA alkylation","DNA alkylation","CD20 (ADCC/CDC)","GR (transrepression)","C5aR1 blockade")
    ), rownames=FALSE, options=list(dom="t", paging=FALSE))
  })

  ## ── Immunology plots (Tab 3) ──────────────────────────────────────────────
  output$plot_anca <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "ANCA_titer", "ANCA Titer (relative units)",
                "Anti-MPO ANCA Titer", hline=1.0)
  })
  output$plot_bcell <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "PLASMA", "Plasma Cell Index", "Plasma Cell Dynamics")
  })
  output$plot_c5a_pmn <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "PMN_ACT", "PMN Activation Index", "Neutrophil Activation")
  })
  output$plot_endo <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "ENDO_INJ", "Endothelial Injury Index", "Endothelial Injury")
  })

  ## ── Renal plots (Tab 4) ───────────────────────────────────────────────────
  output$plot_gfr <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "GFR_obs", "GFR (mL/min/1.73m²)",
                "Glomerular Filtration Rate", hline=15)
  })
  output$plot_creat <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "Creat_est", "Creatinine (μmol/L)", "Serum Creatinine", hline=500)
  })
  output$plot_renal_i <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "RENAL_I", "Renal Inflammation Index", "Renal Inflammation (RPGN Activity)")
  })
  output$plot_renal_f <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "RENAL_F", "Renal Fibrosis Index", "Renal Fibrosis Progression")
  })
  output$renal_summary_tbl <- renderDT({
    df <- sim_data()
    if (is.null(df)) return(datatable(data.frame()))
    df %>%
      filter(time_days >= max(time_days) - 1) %>%
      group_by(label) %>%
      summarise(
        `GFR (final)` = round(mean(GFR_obs),1),
        `Creatinine (final)` = round(mean(Creat_est),0),
        `Renal Inflammation` = round(mean(RENAL_I),2),
        `Fibrosis Index` = round(mean(RENAL_F),3),
        .groups="drop"
      ) %>%
      datatable(rownames=FALSE, options=list(dom="t",paging=FALSE))
  })

  ## ── Pulmonary plots (Tab 5) ───────────────────────────────────────────────
  output$plot_pulm_i <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "PULM_I", "DAH Index", "Pulmonary/DAH Index")
  })
  output$plot_dlco <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "DLCO_obs", "DLCO (% predicted)", "DLCO Trajectory", hline=70)
  })
  output$plot_crp <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "CRP_C", "CRP (mg/L)", "C-Reactive Protein", hline=10)
  })
  output$dah_stage_table <- renderTable({
    data.frame(
      Stage  = c("Subclinical","Mild","Moderate","Severe"),
      Index  = c("< 0.5","0.5–2.0","2.0–4.0","≥ 4.0"),
      DLCO   = c("> 80%","70–80%","50–70%","< 50%"),
      Action = c("Monitor","Steroid","CY/RTX","PLEX+ICU")
    )
  }, striped=TRUE, bordered=TRUE)

  ## ── Scenario comparison plots (Tab 6) ─────────────────────────────────────
  output$comp_bvas <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "BVAS_analog", "BVAS (continuous)", "Disease Activity (BVAS)")
  })
  output$comp_gfr <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "GFR_obs", "GFR (mL/min/1.73m²)", "GFR Comparison", hline=15)
  })
  output$comp_anca <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "ANCA_titer", "ANCA Titer", "ANCA Comparison", hline=1.0)
  })
  output$comp_dlco <- renderPlotly({
    df <- sim_data()
    make_plotly(df, "DLCO_obs", "DLCO (% predicted)", "DLCO Comparison", hline=70)
  })
  output$outcome_summary_tbl <- renderDT({
    df <- sim_data()
    if (is.null(df)) return(datatable(data.frame()))
    df %>%
      filter(time_days >= max(time_days) - 1) %>%
      group_by(label) %>%
      summarise(
        `ANCA Titer`   = round(mean(ANCA_titer), 2),
        `GFR`          = round(mean(GFR_obs), 1),
        `DLCO (%)`     = round(mean(DLCO_obs), 1),
        `BVAS`         = round(mean(BVAS_analog), 2),
        `CRP (mg/L)`   = round(mean(CRP_C), 1),
        `Remission`    = ifelse(mean(in_remission) > 0.5, "Yes", "No"),
        .groups="drop"
      ) %>%
      datatable(rownames=FALSE,
                options=list(dom="t",paging=FALSE),
                rowCallback=JS("function(row,data){if(data[6]=='Yes'){$('td',row).css('background-color','#d4edda');}else{$('td',row).css('background-color','#f8d7da');}}"))
  })

}

## ── Launch ───────────────────────────────────────────────────────────────────
shinyApp(ui=ui, server=server)
