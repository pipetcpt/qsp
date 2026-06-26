## =============================================================================
## Preeclampsia QSP Shiny Dashboard
## =============================================================================
## Interactive simulation & visualization for the Preeclampsia QSP model
## Tabs: Patient Profile · Drug PK · Angiogenic Balance · Cardio-Renal ·
##       HELLP & Neuro · Scenario Comparison · Biomarkers · About
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)

## ─────────────────────────────────────────────────────────────────────────────
## COMPILE mrgsolve MODEL (inline)
## ─────────────────────────────────────────────────────────────────────────────
pe_code <- '
$PARAM
  GA0=84, PE_severity=0.7,
  kprod_sFlt1=0.015, kel_sFlt1=0.04,
  kprod_PlGF=0.08, kel_PlGF=0.1,
  kprod_sEng=0.01, kel_sEng=0.05,
  k_sFlt1_rise=0.025, wk20=140,
  kNO_base=1.0, kNO_inhib=0.4,
  kET1_base=1.0, kET1_stim=0.3,
  kROS_base=1.0, kROS_stim=0.25,
  SBP_base=118, DBP_base=74,
  kSVR_NO=15, kSVR_ET1=8, kSVR_AngII=5,
  kBP_NO_damp=0.1, kBP_Asp=3,
  GFR_base=140, kGFR_sFlt1=0.15,
  kProt_base=150, kProt_rise=8, kProt_GFR=0.3,
  PLT_base=250, kPLT_TXA2=0.15, kPLT_regen=0.03,
  LDH_base=180, kLDH_mth=0.02,
  SeizThresh_base=1.0, kSeiz_MAP=0.15, kMg_protect=0.3,
  ka_Asp=2.0, F_Asp=0.4, Vd_Asp=0.14, CL_Asp=0.6,
  ke_Sal=0.04, Emax_COX1=1.0, EC50_COX1=0.1, koff_COX1=0.007,
  ka_Lab=0.5, F_Lab=0.25, Vd_Lab=9.4, CL_Lab=1.8,
  Emax_Lab=25, EC50_Lab=200,
  ka_Nif=0.2, F_Nif=0.75, Vd_Nif=0.8, CL_Nif=0.5,
  Emax_Nif=20, EC50_Nif=20,
  Vd_Mg=0.55, CL_Mg=0.12, Mg_baseline=0.8,
  Emax_NMDA=0.8, EC50_NMDA=2.0, Mg_tox1=3.5, Mg_tox2=5.0,
  BWT=70

$CMT DEPOT_ASP ASPIRIN SALICYLATE COX1_INH
     DEPOT_LAB LABETALOL
     DEPOT_NIF NIFEDIPINE
     MG_PLASMA
     SFLT1 PLGF SENG
     NO_EA ET1 ROS
     SBP DBP GFR_C PROTEINURIA PLATELET LDH_MK SEIZURE_RISK

$INIT
  MG_PLASMA=38.5, SFLT1=1000, PLGF=80, SENG=2.0,
  NO_EA=1.0, ET1=1.5, ROS=1.0,
  SBP=118, DBP=74, GFR_C=140, PROTEINURIA=150,
  PLATELET=250, LDH_MK=180, SEIZURE_RISK=0.05

$ODE
  double GA = GA0 + TIME;
  double fGest_sFlt1 = (GA > wk20) ?
    1.0 + PE_severity * k_sFlt1_rise * (GA - wk20) : 1.0;

  double C_Asp  = ASPIRIN   / (Vd_Asp * BWT);
  double C_Lab  = LABETALOL / (Vd_Lab * BWT);
  double C_Nif  = NIFEDIPINE/ (Vd_Nif * BWT);
  double Mg_conc = MG_PLASMA / (Vd_Mg * BWT);
  if (C_Asp  < 0) C_Asp  = 0;
  if (C_Lab  < 0) C_Lab  = 0;
  if (C_Nif  < 0) C_Nif  = 0;
  if (Mg_conc < 0) Mg_conc = 0;

  dxdt_DEPOT_ASP  = -ka_Asp * DEPOT_ASP;
  dxdt_ASPIRIN    =  F_Asp * ka_Asp * DEPOT_ASP - (CL_Asp/Vd_Asp) * C_Asp * Vd_Asp * BWT;
  dxdt_SALICYLATE = (CL_Asp/Vd_Asp) * C_Asp * Vd_Asp * BWT - ke_Sal * SALICYLATE;
  double Imax_COX1 = Emax_COX1 * C_Asp / (EC50_COX1 + C_Asp);
  dxdt_COX1_INH   = Imax_COX1 * (1.0 - COX1_INH) - koff_COX1 * COX1_INH;

  dxdt_DEPOT_LAB  = -ka_Lab * DEPOT_LAB;
  dxdt_LABETALOL  = F_Lab * ka_Lab * DEPOT_LAB - CL_Lab * BWT * C_Lab;

  dxdt_DEPOT_NIF  = -ka_Nif * DEPOT_NIF;
  dxdt_NIFEDIPINE = F_Nif * ka_Nif * DEPOT_NIF - CL_Nif * BWT * C_Nif;

  double GFR_frac  = GFR_C / 120.0;
  double CL_Mg_adj = CL_Mg * GFR_frac;
  dxdt_MG_PLASMA  = -(CL_Mg_adj * BWT) * (Mg_conc - Mg_baseline);

  double prod_sFlt1 = kprod_sFlt1 * fGest_sFlt1;
  dxdt_SFLT1 = prod_sFlt1 * 1000 - kel_sFlt1 * SFLT1;
  double PlGF_sup = 1.0 / (1.0 + (SFLT1/5000.0) * PE_severity * 0.8);
  dxdt_PLGF  = kprod_PlGF * 100 * PlGF_sup - kel_PlGF * PLGF;
  dxdt_SENG  = kprod_sEng * fGest_sFlt1 * 3 - kel_sEng * SENG;

  double NO_target = kNO_base - kNO_inhib * (SFLT1/5000.0) * PE_severity;
  if (NO_target < 0.1) NO_target = 0.1;
  dxdt_NO_EA = 0.2 * (NO_target - NO_EA);

  double ET1_target = kET1_base + kET1_stim * (1.0 - NO_EA);
  dxdt_ET1 = 0.3 * (ET1_target - ET1);

  double ROS_target = kROS_base + kROS_stim * (SFLT1/5000.0) * PE_severity;
  dxdt_ROS = 0.2 * (ROS_target - ROS);

  double E_Lab = Emax_Lab * C_Lab / (EC50_Lab + C_Lab);
  double E_Nif = Emax_Nif * C_Nif / (EC50_Nif + C_Nif);
  double E_Asp_BP = kBP_Asp * COX1_INH;

  double SBP_target = SBP_base + kSVR_NO*(1.0-NO_EA) + kSVR_ET1*(ET1-kET1_base) +
    kSVR_AngII*PE_severity - E_Lab - E_Nif - E_Asp_BP;
  double DBP_target = DBP_base + 0.7*kSVR_NO*(1.0-NO_EA) + 0.6*kSVR_ET1*(ET1-kET1_base) +
    0.5*kSVR_AngII*PE_severity - 0.7*E_Lab - 0.7*E_Nif - 0.5*E_Asp_BP;

  dxdt_SBP = kBP_NO_damp * 24 * (SBP_target - SBP);
  dxdt_DBP = kBP_NO_damp * 24 * (DBP_target - DBP);

  double MAP_val = DBP + (SBP-DBP)/3.0;
  double MAP_excess = MAP_val - 90.0;
  double GFR_target = GFR_base - kGFR_sFlt1*(SFLT1-2000.0)/1000.0*PE_severity -
    (MAP_excess > 0 ? 0.5*MAP_excess : 0);
  if (GFR_target < 15) GFR_target = 15;
  dxdt_GFR_C = 0.1 * (GFR_target - GFR_C);

  double Prot_target = kProt_base + kProt_rise*(SFLT1-2000.0)/1000.0*PE_severity +
    kProt_GFR*(GFR_base-GFR_C);
  if (Prot_target < 0) Prot_target = 0;
  dxdt_PROTEINURIA = 0.15 * (Prot_target - PROTEINURIA);

  double TXA2_idx = (1.0 - COX1_INH) * PE_severity * 0.5;
  double PLT_target = PLT_base * (1.0 - TXA2_idx*kPLT_TXA2*10);
  if (PLT_target < 10) PLT_target = 10;
  dxdt_PLATELET = kPLT_regen * (PLT_target - PLATELET);
  double LDH_target = LDH_base + 200*TXA2_idx*kLDH_mth*100;
  dxdt_LDH_MK = 0.05 * (LDH_target - LDH_MK);

  double MAP_neurisk = MAP_val - 100.0;
  double Mg_protect_fx = Emax_NMDA*(Mg_conc-Mg_baseline)/(EC50_NMDA+(Mg_conc-Mg_baseline));
  if (Mg_protect_fx < 0) Mg_protect_fx = 0;
  double Seizure_target = 0.05 + 0.003*(MAP_neurisk>0?MAP_neurisk:0)*kSeiz_MAP -
    kMg_protect*Mg_protect_fx;
  if (Seizure_target < 0) Seizure_target = 0;
  if (Seizure_target > 1.0) Seizure_target = 1.0;
  dxdt_SEIZURE_RISK = 0.2 * (Seizure_target - SEIZURE_RISK);

$TABLE
  double GA_weeks = (GA0 + TIME) / 7.0;
  double MAP_out  = DBP + (SBP - DBP) / 3.0;
  double sFlt1_PlGF_ratio = (PLGF > 0) ? SFLT1/PLGF : 9999;
  double Mg_plasma_mmolL  = MG_PLASMA / (Vd_Mg * BWT);
  double Lab_ngmL = LABETALOL / (Vd_Lab * BWT);
  double Nif_ngmL = NIFEDIPINE / (Vd_Nif * BWT);
  double Asp_mgL  = ASPIRIN   / (Vd_Asp * BWT);
  capture GA_weeks MAP_out sFlt1_PlGF_ratio Mg_plasma_mmolL Lab_ngmL Nif_ngmL Asp_mgL;
'
base_mod <- mcode("pe_shiny", pe_code, quiet = TRUE)

## ─────────────────────────────────────────────────────────────────────────────
## SIMULATION FUNCTION
## ─────────────────────────────────────────────────────────────────────────────
run_simulation <- function(severity, dose_asp, dose_lab, dose_nif,
                            mg_load, mg_maint_gh, start_asp_wk,
                            start_trt_wk, bwt = 70) {

  sim_end_days <- 196  # 12 → 40 wk
  delta_t <- 0.5

  mod_u <- param(base_mod,
                  PE_severity = severity,
                  SBP_base = 118, DBP_base = 74,
                  BWT = bwt)

  ev_list <- list()

  # Aspirin
  if (dose_asp > 0) {
    t0 <- max(0, (start_asp_wk - 12) * 7 * 24)
    addl_asp <- max(0, floor((sim_end_days - t0/24) * 1))
    ev_list[["asp"]] <- ev(amt = dose_asp, cmt = 1, time = t0, ii = 24, addl = addl_asp)
  }

  # Labetalol
  if (dose_lab > 0) {
    t_lab <- max(0, (start_trt_wk - 12) * 7 * 24)
    addl_lab <- max(0, floor((sim_end_days - t_lab/24) * 2))
    ev_list[["lab"]] <- ev(amt = dose_lab, cmt = 5, time = t_lab, ii = 12, addl = addl_lab)
  }

  # Nifedipine
  if (dose_nif > 0) {
    t_nif <- max(0, (start_trt_wk - 12) * 7 * 24)
    addl_nif <- max(0, floor(sim_end_days - t_nif/24))
    ev_list[["nif"]] <- ev(amt = dose_nif, cmt = 7, time = t_nif, ii = 24, addl = addl_nif)
  }

  # Magnesium sulfate
  if (mg_load > 0 || mg_maint_gh > 0) {
    t_mg <- max(0, (start_trt_wk - 12) * 7 * 24)
    ev_mg <- list()
    if (mg_load > 0)
      ev_mg[["load"]] <- ev(amt = mg_load * 1000, cmt = 9, time = t_mg)
    if (mg_maint_gh > 0)
      ev_mg[["maint"]] <- ev(amt = mg_maint_gh * 1000, cmt = 9,
                              time = t_mg + 1, ii = 1, addl = 47)
    ev_list <- c(ev_list, ev_mg)
  }

  if (length(ev_list) > 0) {
    ev_all <- Reduce(c, ev_list)
    out <- mrgsim(mod_u, ev = ev_all, end = sim_end_days, delta = delta_t)
  } else {
    out <- mrgsim(mod_u, end = sim_end_days, delta = delta_t)
  }

  as.data.frame(out)
}

## ─────────────────────────────────────────────────────────────────────────────
## UI
## ─────────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "purple",

  dashboardHeader(title = "Preeclampsia QSP Dashboard", titleWidth = 320),

  dashboardSidebar(
    width = 290,
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("Drug PK",              tabName = "tab_pk",       icon = icon("pills")),
      menuItem("Angiogenic Balance",   tabName = "tab_angio",    icon = icon("dna")),
      menuItem("Cardio-Renal",         tabName = "tab_cardio",   icon = icon("heartbeat")),
      menuItem("HELLP & Neuro",        tabName = "tab_hellp",    icon = icon("brain")),
      menuItem("Scenario Comparison",  tabName = "tab_scenario", icon = icon("chart-bar")),
      menuItem("Biomarker Panel",      tabName = "tab_bio",      icon = icon("vials")),
      menuItem("About",                tabName = "tab_about",    icon = icon("info-circle"))
    ),
    br(),
    h5("  Disease Parameters", style = "color:#CCC; padding-left:10px"),
    sliderInput("severity", "PE Severity (0=none, 1=severe)",
                min = 0, max = 1, value = 0.7, step = 0.05),
    sliderInput("bwt", "Body Weight (kg)", min = 45, max = 120, value = 70),
    br(),
    h5("  Drug Dosing", style = "color:#CCC; padding-left:10px"),
    sliderInput("dose_asp",  "Aspirin (mg/day; 0=off)", min = 0, max = 150, value = 75, step = 25),
    sliderInput("start_asp", "Aspirin start (wk GA)",   min = 12, max = 20, value = 12),
    sliderInput("dose_lab",  "Labetalol (mg BID; 0=off)", min = 0, max = 400, value = 0, step = 50),
    sliderInput("dose_nif",  "Nifedipine MR (mg/day; 0=off)", min = 0, max = 90, value = 0, step = 10),
    sliderInput("start_trt", "Antihypertensive start (wk)", min = 20, max = 36, value = 24),
    sliderInput("mg_load",   "MgSO4 Loading dose (g IV; 0=off)", min = 0, max = 6, value = 0, step = 1),
    sliderInput("mg_maint",  "MgSO4 Maintenance (g/h; 0=off)",   min = 0, max = 3, value = 0, step = 0.5),
    br(),
    actionButton("run_sim", "Run Simulation", icon = icon("play"),
                 style = "background-color:#7B1FA2; color:white; width:90%; margin:5%")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box { border-top: 3px solid #7B1FA2; }
      .nav-tabs-custom > .nav-tabs > li.active > a { border-top-color: #7B1FA2; }
      .badge-danger { background-color: #E63946; }
    "))),

    tabItems(

      ## ── TAB 1: Patient Profile ──────────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(width = 12, title = "Patient Profile & Disease Overview",
              status = "purple", solidHeader = TRUE,
              fluidRow(
                column(4,
                  h4("Current Parameters"),
                  tableOutput("tbl_params")
                ),
                column(4,
                  h4("Clinical Status at 36 Weeks"),
                  uiOutput("clinical_status_box")
                ),
                column(4,
                  h4("Key Thresholds"),
                  tableOutput("tbl_thresholds")
                )
              )
          )
        ),
        fluidRow(
          box(width = 6, title = "sFlt-1/PlGF Ratio Trajectory",
              plotOutput("p_ratio_overview", height = 300)),
          box(width = 6, title = "Blood Pressure Trajectory",
              plotOutput("p_bp_overview", height = 300))
        )
      ),

      ## ── TAB 2: Drug PK ────────────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(width = 12, title = "Drug Pharmacokinetics",
              status = "blue", solidHeader = TRUE,
              tabsetPanel(
                tabPanel("Aspirin",   plotOutput("p_asp_pk",  height = 350)),
                tabPanel("Labetalol", plotOutput("p_lab_pk",  height = 350)),
                tabPanel("Nifedipine",plotOutput("p_nif_pk",  height = 350)),
                tabPanel("Magnesium", plotOutput("p_mg_pk",   height = 350)),
                tabPanel("PK Summary Table", DTOutput("tbl_pk_summary"))
              )
          )
        )
      ),

      ## ── TAB 3: Angiogenic Balance ──────────────────────────────────────────
      tabItem(tabName = "tab_angio",
        fluidRow(
          box(width = 6, title = "sFlt-1 Plasma Levels",
              plotOutput("p_sflt1", height = 300)),
          box(width = 6, title = "PlGF Plasma Levels",
              plotOutput("p_plgf", height = 300))
        ),
        fluidRow(
          box(width = 6, title = "sFlt-1 / PlGF Ratio",
              plotOutput("p_ratio", height = 300)),
          box(width = 6, title = "Soluble Endoglin (sEng)",
              plotOutput("p_seng", height = 300))
        )
      ),

      ## ── TAB 4: Cardio-Renal ──────────────────────────────────────────────
      tabItem(tabName = "tab_cardio",
        fluidRow(
          box(width = 6, title = "Systolic & Diastolic BP",
              plotOutput("p_bp", height = 300)),
          box(width = 6, title = "Endothelial Markers (NO, ET-1, ROS)",
              plotOutput("p_endothel", height = 300))
        ),
        fluidRow(
          box(width = 6, title = "GFR Trajectory",
              plotOutput("p_gfr", height = 300)),
          box(width = 6, title = "Proteinuria",
              plotOutput("p_prot", height = 300))
        )
      ),

      ## ── TAB 5: HELLP & Neuro ─────────────────────────────────────────────
      tabItem(tabName = "tab_hellp",
        fluidRow(
          box(width = 6, title = "Platelet Count",
              plotOutput("p_plt", height = 300)),
          box(width = 6, title = "LDH (Hemolysis Marker)",
              plotOutput("p_ldh", height = 300))
        ),
        fluidRow(
          box(width = 6, title = "Seizure Risk & Eclampsia",
              plotOutput("p_seiz", height = 300)),
          box(width = 6, title = "Magnesium Plasma Level & Safety",
              plotOutput("p_mg_safety", height = 300))
        )
      ),

      ## ── TAB 6: Scenario Comparison ───────────────────────────────────────
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(width = 12, title = "Six Treatment Scenario Comparison",
              status = "purple", solidHeader = TRUE,
              plotOutput("p_scenarios", height = 600))
        ),
        fluidRow(
          box(width = 12, title = "Comparative Summary at GA 36 weeks",
              DTOutput("tbl_scenario_comparison"))
        )
      ),

      ## ── TAB 7: Biomarker Panel ───────────────────────────────────────────
      tabItem(tabName = "tab_bio",
        fluidRow(
          box(width = 12, title = "Biomarker Heatmap Over Gestation",
              status = "purple", solidHeader = TRUE,
              plotOutput("p_biomarker_heatmap", height = 500))
        ),
        fluidRow(
          box(width = 12, title = "Biomarker Data Table",
              DTOutput("tbl_biomarkers"))
        )
      ),

      ## ── TAB 8: About ──────────────────────────────────────────────────────
      tabItem(tabName = "tab_about",
        fluidRow(
          box(width = 12, title = "About This Model", status = "purple", solidHeader = TRUE,
            HTML("
              <h4>Preeclampsia QSP Model</h4>
              <p><strong>Disease:</strong> Preeclampsia (PE) — a pregnancy complication
              affecting 2–8% of pregnancies worldwide, characterized by new-onset
              hypertension (≥140/90 mmHg) after 20 weeks gestation with proteinuria
              or end-organ damage.</p>

              <h5>Model Structure (20 ODEs)</h5>
              <ul>
                <li><strong>Angiogenic imbalance:</strong> sFlt-1, PlGF, sEng dynamics</li>
                <li><strong>Endothelial dysfunction:</strong> NO bioavailability, ET-1, ROS</li>
                <li><strong>Cardiovascular:</strong> SBP, DBP (SVR-driven, drug-modulated)</li>
                <li><strong>Renal:</strong> GFR decline, proteinuria (glomerular endotheliosis)</li>
                <li><strong>Coagulation/HELLP:</strong> Platelet count, LDH</li>
                <li><strong>Neurological:</strong> Seizure risk (Mg-NMDA axis)</li>
                <li><strong>Drug PK:</strong> Aspirin, Labetalol, Nifedipine, MgSO4</li>
              </ul>

              <h5>Clinical Calibration</h5>
              <ul>
                <li>ASPRE trial (Rolnik et al., <em>Lancet</em> 2017): aspirin 75–150 mg from
                11–14 wk reduces early PE by 62%</li>
                <li>Maynard et al. (<em>J Clin Invest</em> 2003): sFlt-1 overexpression causes
                PE phenotype in rats</li>
                <li>Verlohren et al. (<em>Am J Obstet Gynecol</em> 2010): sFlt-1/PlGF ratio
                cutoff ≥38 (sensitivity 82%, specificity 95%) predicts PE within 4 weeks</li>
                <li>CHIPS trial (Magee et al., <em>NEJM</em> 2015): tight BP control
                (target DBP 85) vs. less tight (target 100)</li>
                <li>Magpie trial (Altman et al., <em>Lancet</em> 2002): MgSO4 reduces
                eclampsia risk by 58%</li>
              </ul>

              <h5>Key Biomarker Thresholds</h5>
              <table class='table table-bordered' style='width:50%'>
                <tr><th>Biomarker</th><th>Normal</th><th>PE Risk</th><th>Severe</th></tr>
                <tr><td>SBP (mmHg)</td><td>&lt;140</td><td>≥140</td><td>≥160</td></tr>
                <tr><td>DBP (mmHg)</td><td>&lt;90</td><td>≥90</td><td>≥110</td></tr>
                <tr><td>sFlt-1/PlGF ratio</td><td>&lt;38</td><td>38–85</td><td>&gt;85</td></tr>
                <tr><td>Proteinuria (mg/24h)</td><td>&lt;300</td><td>≥300</td><td>≥5000</td></tr>
                <tr><td>Platelets (×10³/µL)</td><td>&gt;150</td><td>100–150</td><td>&lt;100</td></tr>
                <tr><td>LDH (IU/L)</td><td>&lt;600</td><td>600–800</td><td>&gt;800</td></tr>
                <tr><td>Mg²⁺ plasma (mmol/L)</td><td>0.7–1.0</td><td>1.7–3.5 (therapeutic)</td><td>&gt;3.5 (toxic)</td></tr>
              </table>

              <p><em>Model by Claude Code Routine (CCR) — QSP Disease Model Library</em></p>
            ")
          )
        )
      )
    )
  )
)

## ─────────────────────────────────────────────────────────────────────────────
## SERVER
## ─────────────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  ## Reactive simulation data
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running QSP simulation...", value = 0.5, {
      run_simulation(
        severity     = input$severity,
        dose_asp     = input$dose_asp,
        dose_lab     = input$dose_lab,
        dose_nif     = input$dose_nif,
        mg_load      = input$mg_load,
        mg_maint_gh  = input$mg_maint,
        start_asp_wk = input$start_asp,
        start_trt_wk = input$start_trt,
        bwt          = input$bwt
      )
    })
  }, ignoreNULL = FALSE)

  ## Scenario comparison (pre-defined 6 scenarios)
  scenario_data <- reactive({
    sev <- input$severity
    bwt <- input$bwt

    s1 <- run_simulation(sev, 0,  0,  0,  0,   0, 12, 24, bwt); s1$scen <- "1. No treatment"
    s2 <- run_simulation(sev, 75, 0,  0,  0,   0, 12, 24, bwt); s2$scen <- "2. Aspirin 75mg/d"
    s3 <- run_simulation(sev, 0,  200,0,  0,   0, 12, 24, bwt); s3$scen <- "3. Labetalol 200mg BID"
    s4 <- run_simulation(sev, 0,  0, 30,  0,   0, 12, 24, bwt); s4$scen <- "4. Nifedipine 30mg/d"
    s5 <- run_simulation(sev, 0,  0,  0,  4,   1, 12, 30, bwt); s5$scen <- "5. MgSO4 (load+maint)"
    s6 <- run_simulation(sev, 75,200,  0,  4,   1, 12, 24, bwt); s6$scen <- "6. Aspirin+Lab+Mg"
    bind_rows(s1, s2, s3, s4, s5, s6)
  })

  ## Shared ggplot theme
  thm <- theme_bw(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  pal <- c("#E63946","#457B9D","#2A9D8F","#E9C46A","#F4A261","#264653")

  ## ── TAB 1 Outputs ────────────────────────────────────────────────────────
  output$tbl_params <- renderTable({
    data.frame(
      Parameter = c("PE Severity", "Body Weight (kg)",
                    "Aspirin dose (mg/d)", "Labetalol (mg BID)",
                    "Nifedipine (mg/d)", "MgSO4 load (g)"),
      Value = c(input$severity, input$bwt,
                input$dose_asp, input$dose_lab,
                input$dose_nif, input$mg_load)
    )
  })

  output$tbl_thresholds <- renderTable({
    data.frame(
      Marker = c("SBP PE", "SBP Severe", "DBP PE", "DBP Severe",
                 "sFlt1/PlGF", "Proteinuria", "Platelets"),
      Threshold = c("≥140 mmHg", "≥160 mmHg", "≥90 mmHg", "≥110 mmHg",
                    ">38 (>85 severe)", "≥300 mg/24h", "<100 k/µL (HELLP)")
    )
  })

  output$clinical_status_box <- renderUI({
    d <- sim_data()
    row <- d %>% filter(abs(GA_weeks - 36) == min(abs(GA_weeks - 36))) %>% slice(1)
    pe_flag    <- row$SBP >= 140 || row$DBP >= 90
    sev_flag   <- row$SBP >= 160 || row$DBP >= 110
    hellp_flag <- row$PLATELET < 100 && row$LDH_MK > 600
    ratio_val  <- round(row$sFlt1_PlGF_ratio, 1)

    status_color <- if (sev_flag || hellp_flag) "red"
                    else if (pe_flag) "orange" else "green"
    status_text  <- if (sev_flag) "SEVERE PREECLAMPSIA"
                    else if (hellp_flag) "HELLP SYNDROME"
                    else if (pe_flag) "PREECLAMPSIA"
                    else "NORMOTENSIVE"

    tagList(
      tags$div(style = paste0("background:", status_color, "; color:white; padding:10px;
                               border-radius:5px; font-weight:bold; font-size:1.1em; margin-bottom:10px"),
               status_text),
      tags$p(paste0("SBP: ", round(row$SBP, 0), " | DBP: ", round(row$DBP, 0), " mmHg")),
      tags$p(paste0("sFlt-1/PlGF ratio: ", ratio_val)),
      tags$p(paste0("Proteinuria: ", round(row$PROTEINURIA, 0), " mg/24h")),
      tags$p(paste0("Platelets: ", round(row$PLATELET, 0), " ×10³/µL")),
      tags$p(paste0("Seizure Risk: ", round(row$SEIZURE_RISK * 100, 1), "%"))
    )
  })

  output$p_ratio_overview <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = GA_weeks, y = sFlt1_PlGF_ratio)) +
      geom_line(color = "#7B1FA2", linewidth = 1) +
      geom_hline(yintercept = 38, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 85, linetype = "dotted", color = "red") +
      scale_y_log10() +
      labs(x = "Gestational Age (weeks)", y = "sFlt-1/PlGF ratio (log)") + thm
  })

  output$p_bp_overview <- renderPlot({
    d <- sim_data()
    d %>% select(GA_weeks, SBP, DBP) %>%
      pivot_longer(c(SBP, DBP)) %>%
      ggplot(aes(x = GA_weeks, y = value, color = name)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 140, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 160, linetype = "dotted", color = "red") +
      scale_color_manual(values = c("SBP" = "#E63946", "DBP" = "#457B9D")) +
      labs(x = "Gestational Age (weeks)", y = "Blood Pressure (mmHg)") + thm
  })

  ## ── TAB 2: Drug PK ───────────────────────────────────────────────────────
  output$p_asp_pk <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = GA_weeks)) +
      geom_line(aes(y = Asp_mgL, color = "Aspirin plasma (mg/L)"), linewidth = 1) +
      geom_line(aes(y = COX1_INH, color = "COX-1 Inhibition (0-1)"), linewidth = 1) +
      scale_color_manual(values = c("Aspirin plasma (mg/L)" = "#E63946",
                                    "COX-1 Inhibition (0-1)" = "#2A9D8F")) +
      labs(x = "GA (weeks)", y = "Concentration / Effect", title = "Aspirin PK & COX-1 Inhibition") + thm
  })

  output$p_lab_pk <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = GA_weeks, y = Lab_ngmL)) +
      geom_line(color = "#E9C46A", linewidth = 1) +
      geom_hline(yintercept = 200, linetype = "dashed", color = "gray50") +
      annotate("text", x = 30, y = 205, label = "EC50 (200 ng/mL)", size = 4) +
      labs(x = "GA (weeks)", y = "Labetalol (ng/mL)", title = "Labetalol Plasma Concentration") + thm
  })

  output$p_nif_pk <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = GA_weeks, y = Nif_ngmL)) +
      geom_line(color = "#2A9D8F", linewidth = 1) +
      geom_hline(yintercept = 20, linetype = "dashed", color = "gray50") +
      annotate("text", x = 30, y = 21, label = "EC50 (20 ng/mL)", size = 4) +
      labs(x = "GA (weeks)", y = "Nifedipine (ng/mL)", title = "Nifedipine Plasma Concentration") + thm
  })

  output$p_mg_pk <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = GA_weeks, y = Mg_plasma_mmolL)) +
      geom_line(color = "#457B9D", linewidth = 1) +
      geom_hline(yintercept = 1.7, linetype = "dashed", color = "green4") +
      geom_hline(yintercept = 3.5, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 5.0, linetype = "dotted", color = "red") +
      annotate("text", x = 30, y = 1.75, label = "Therapeutic lower (1.7)", size = 3.5) +
      annotate("text", x = 30, y = 3.55, label = "NMJ block risk (3.5)",    size = 3.5) +
      annotate("text", x = 30, y = 5.05, label = "Resp. depress. (5.0)",    size = 3.5) +
      labs(x = "GA (weeks)", y = "Mg²⁺ (mmol/L)", title = "Magnesium Plasma Level") + thm
  })

  output$tbl_pk_summary <- renderDT({
    d <- sim_data()
    d %>%
      filter(GA_weeks %in% c(20, 24, 28, 32, 36, 38)) %>%
      select(GA_weeks, Asp_mgL, Lab_ngmL, Nif_ngmL, Mg_plasma_mmolL, COX1_INH) %>%
      mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
      datatable(options = list(pageLength = 10))
  })

  ## ── TAB 3: Angiogenic Balance ─────────────────────────────────────────────
  output$p_sflt1 <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = GA_weeks, y = SFLT1)) +
      geom_line(color = "#E63946", linewidth = 1) +
      labs(x = "GA (wk)", y = "sFlt-1 (pg/mL)", title = "sFlt-1 Plasma Levels") + thm
  })

  output$p_plgf <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = GA_weeks, y = PLGF)) +
      geom_line(color = "#2A9D8F", linewidth = 1) +
      labs(x = "GA (wk)", y = "PlGF (pg/mL)", title = "Placental Growth Factor (PlGF)") + thm
  })

  output$p_ratio <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = GA_weeks, y = sFlt1_PlGF_ratio)) +
      geom_line(color = "#7B1FA2", linewidth = 1.2) +
      geom_hline(yintercept = 38, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 85, linetype = "dotted", color = "red") +
      scale_y_log10() +
      labs(x = "GA (wk)", y = "Ratio (log scale)", title = "sFlt-1/PlGF Ratio (Diagnostic)") + thm
  })

  output$p_seng <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = GA_weeks, y = SENG)) +
      geom_line(color = "#E9C46A", linewidth = 1) +
      labs(x = "GA (wk)", y = "sEng (ng/mL)", title = "Soluble Endoglin (Anti-angiogenic)") + thm
  })

  ## ── TAB 4: Cardio-Renal ─────────────────────────────────────────────────
  output$p_bp <- renderPlot({
    d <- sim_data()
    d %>% select(GA_weeks, SBP, DBP) %>%
      pivot_longer(c(SBP, DBP)) %>%
      ggplot(aes(x = GA_weeks, y = value, color = name)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 140, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 160, linetype = "dotted", color = "red") +
      scale_color_manual(values = c(SBP = "#E63946", DBP = "#457B9D")) +
      labs(x = "GA (wk)", y = "mmHg", title = "Blood Pressure") + thm
  })

  output$p_endothel <- renderPlot({
    d <- sim_data()
    d %>% select(GA_weeks, NO_EA, ET1, ROS) %>%
      pivot_longer(c(NO_EA, ET1, ROS)) %>%
      ggplot(aes(x = GA_weeks, y = value, color = name)) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = c(NO_EA = "#2A9D8F", ET1 = "#E63946", ROS = "#E9C46A")) +
      labs(x = "GA (wk)", y = "Relative / Actual", title = "Endothelial Markers") + thm
  })

  output$p_gfr <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = GA_weeks, y = GFR_C)) +
      geom_line(color = "#457B9D", linewidth = 1) +
      geom_hline(yintercept = 60, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 30, linetype = "dotted", color = "red") +
      labs(x = "GA (wk)", y = "GFR (mL/min/1.73m²)", title = "Glomerular Filtration Rate") + thm
  })

  output$p_prot <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = GA_weeks, y = PROTEINURIA)) +
      geom_line(color = "#F4A261", linewidth = 1) +
      geom_hline(yintercept = 300, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 5000, linetype = "dotted", color = "red") +
      labs(x = "GA (wk)", y = "mg/24h", title = "Proteinuria") + thm
  })

  ## ── TAB 5: HELLP & Neuro ─────────────────────────────────────────────────
  output$p_plt <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = GA_weeks, y = PLATELET)) +
      geom_line(color = "#7B1FA2", linewidth = 1) +
      geom_hline(yintercept = 150, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 100, linetype = "dotted", color = "red") +
      labs(x = "GA (wk)", y = "×10³/µL", title = "Platelet Count") + thm
  })

  output$p_ldh <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = GA_weeks, y = LDH_MK)) +
      geom_line(color = "#E63946", linewidth = 1) +
      geom_hline(yintercept = 600, linetype = "dashed", color = "orange") +
      labs(x = "GA (wk)", y = "IU/L", title = "LDH (Hemolysis Marker)") + thm
  })

  output$p_seiz <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = GA_weeks, y = SEIZURE_RISK * 100)) +
      geom_line(color = "#9B2226", linewidth = 1) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
      labs(x = "GA (wk)", y = "Seizure Risk (%)", title = "Eclampsia Risk") + thm
  })

  output$p_mg_safety <- renderPlot({
    d <- sim_data()
    ggplot(d, aes(x = GA_weeks, y = Mg_plasma_mmolL)) +
      geom_line(color = "#264653", linewidth = 1) +
      geom_ribbon(aes(ymin = 1.7, ymax = 3.5), fill = "green", alpha = 0.1) +
      geom_hline(yintercept = 1.7, color = "green4", linetype = "dashed") +
      geom_hline(yintercept = 3.5, color = "orange", linetype = "dashed") +
      geom_hline(yintercept = 5.0, color = "red", linetype = "dotted") +
      labs(x = "GA (wk)", y = "Mg²⁺ (mmol/L)", title = "Plasma Magnesium (Therapeutic Window)") + thm
  })

  ## ── TAB 6: Scenario Comparison ───────────────────────────────────────────
  output$p_scenarios <- renderPlot({
    d <- scenario_data()
    p1 <- ggplot(d, aes(x = GA_weeks, y = SBP, color = scen)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 140, linetype = "dashed") +
      scale_color_manual(values = pal) +
      labs(x = "GA (wk)", y = "SBP (mmHg)", title = "SBP") + thm +
      theme(legend.position = "none")

    p2 <- ggplot(d, aes(x = GA_weeks, y = sFlt1_PlGF_ratio, color = scen)) +
      geom_line(linewidth = 0.9) +
      scale_y_log10() +
      geom_hline(yintercept = 38, linetype = "dashed") +
      scale_color_manual(values = pal) +
      labs(x = "GA (wk)", y = "sFlt-1/PlGF (log)", title = "Angiogenic Ratio") + thm +
      theme(legend.position = "none")

    p3 <- ggplot(d, aes(x = GA_weeks, y = PROTEINURIA, color = scen)) +
      geom_line(linewidth = 0.9) +
      geom_hline(yintercept = 300, linetype = "dashed") +
      scale_color_manual(values = pal) +
      labs(x = "GA (wk)", y = "Proteinuria (mg/24h)", title = "Proteinuria") + thm +
      theme(legend.position = "none")

    p4 <- ggplot(d, aes(x = GA_weeks, y = SEIZURE_RISK * 100, color = scen)) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = pal) +
      labs(x = "GA (wk)", y = "Seizure Risk (%)", title = "Eclampsia Risk") + thm +
      theme(legend.position = "bottom")

    gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)
  })

  output$tbl_scenario_comparison <- renderDT({
    d <- scenario_data()
    d %>%
      group_by(scen) %>%
      filter(abs(GA_weeks - 36) == min(abs(GA_weeks - 36))) %>%
      slice(1) %>%
      select(Scenario = scen, SBP, DBP,
             `sFlt1/PlGF` = sFlt1_PlGF_ratio,
             `GFR` = GFR_C,
             `Proteinuria` = PROTEINURIA,
             `Platelets` = PLATELET,
             `LDH` = LDH_MK,
             `Seizure Risk` = SEIZURE_RISK) %>%
      mutate(across(where(is.numeric), ~ round(.x, 1))) %>%
      datatable(options = list(pageLength = 10))
  })

  ## ── TAB 7: Biomarker Heatmap ─────────────────────────────────────────────
  output$p_biomarker_heatmap <- renderPlot({
    d <- sim_data()
    markers <- c("SFLT1","PLGF","SENG","NO_EA","ET1","ROS",
                 "SBP","DBP","GFR_C","PROTEINURIA",
                 "PLATELET","LDH_MK","SEIZURE_RISK")

    d_heat <- d %>%
      filter(GA_weeks %% 2 < 0.2) %>%
      select(GA_weeks, all_of(markers)) %>%
      mutate(across(all_of(markers), ~ (.x - min(.x)) / (max(.x) - min(.x) + 1e-9))) %>%
      pivot_longer(all_of(markers), names_to = "Marker", values_to = "Normalized")

    ggplot(d_heat, aes(x = GA_weeks, y = Marker, fill = Normalized)) +
      geom_tile() +
      scale_fill_gradientn(colors = c("#2A9D8F","#E9C46A","#E63946"),
                           name = "Normalized\nIntensity") +
      labs(x = "Gestational Age (weeks)", y = NULL,
           title = "Biomarker Heatmap (Normalized) Over Gestation") +
      theme_bw(base_size = 12) + theme(panel.grid = element_blank())
  })

  output$tbl_biomarkers <- renderDT({
    d <- sim_data()
    d %>%
      filter(GA_weeks %in% c(12, 16, 20, 24, 28, 32, 36, 38)) %>%
      select(GA_weeks, SBP, DBP, MAP_out, sFlt1_PlGF_ratio,
             GFR_C, PROTEINURIA, PLATELET, LDH_MK, SEIZURE_RISK,
             Mg_plasma_mmolL) %>%
      mutate(across(where(is.numeric), ~ round(.x, 1))) %>%
      datatable(options = list(pageLength = 10))
  })
}

## ─────────────────────────────────────────────────────────────────────────────
## LAUNCH
## ─────────────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
