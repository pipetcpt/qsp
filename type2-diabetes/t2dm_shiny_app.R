## =============================================================================
##  Type 2 Diabetes QSP — Interactive Shiny Dashboard
##  7 tabs: Patient Profile · Drug PK · Glucose/Insulin · β-cell ·
##          Clinical Endpoints · Scenario Comparison · Biomarkers & Safety
## =============================================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(plotly)
library(mrgsolve)

# ============================================================================
# Inline model (same as t2dm_mrgsolve_model.R but compact for Shiny)
# ============================================================================
t2dm_code <- '
$PARAM
BW=90, HbA1c0=9.0, Gp0=180, Ip0=20, Gc0=150, eGFR0=80, UACR0=50,
USE_MET=0, USE_EMPA=0, USE_SEMA=0, USE_DPP4=0, USE_SU=0, USE_INS=0, USE_PIOG=0,
ka_met=1.2, F_met=0.55, Vc_met=300, k12_met=0.15, k21_met=0.08, CL_met=35,
ka_empa=0.60, Vc_empa=73, CL_empa=10.6,
ka_sema=0.0085, Vc_sema=12.5, CL_sema=0.053,
ka_dpp4=0.80, Vc_dpp4=198, CL_dpp4=12.4,
ka_su=0.50, Vc_su=12.6, CL_su=3.1,
ka_ins=0.048, Vc_ins=8.0, CL_ins=0.96,
ka_piog=0.70, Vc_piog=89, CL_piog=5.6,
SI=8e-4, Sg=0.01, EGP0=2.4, Rd0=2.4, p2=0.05, Vg=1.5, Vi=0.05,
beta_sens=0.6, beta_M0=1.0, k_prolif=1e-4, k_apop=1e-3,
Gth=90, phi_max=900,
kout_Gc=0.3, Gc_Gp50=100, EGP_Gc=0.015,
kin_GLP1=15, kout_GLP1=4.0, kDPP4=3.5, GLP1_Gp=0.05,
Imax_dpp4=0.80, IC50_dpp4=100,
TmG_base=340, RGT0=180, Imax_empa=0.55, IC50_empa=30000, GFR_val=100,
Emax_su=2.5, EC50_su=50, n_su=1.5,
Emax_sema_ins=0.6, EC50_sema=5, Emax_sema_wt=4.5, ksema_wt=0.004,
IR_H0=2.5, IR_P0=2.0, k_IR_FFA=0.005, k_IR_rec=0.001,
FFA0=0.8, kFFA_rel=0.15, kFFA_up=0.20, FFA_Ip50=15,
kHbA1c=0.0084, HbA1c_ss=0.165,
k_eGFR_decline=2e-5, k_UACR_rise=0.004, k_eGFR_empa=1.5e-5, k_UACR_empa=0.003,
Emax_piog_IR=0.35, EC50_piog=200000, Emax_piog_wt=3.0,
BW_rate=0

$CMT
MET_GUT MET_C MET_P EMPA_C SEMA_SC SEMA_C DPP4I_C SU_C INS_SC INS_C PIOG_C
Gp Gt Ip X_action Gc GLP1 beta_mass IR_H IR_P FFA BW_t
HbA1c_cmpt eGFR_cmpt UACR_cmpt

$INIT
MET_GUT=0, MET_C=0, MET_P=0, EMPA_C=0, SEMA_SC=0, SEMA_C=0,
DPP4I_C=0, SU_C=0, INS_SC=0, INS_C=0, PIOG_C=0,
Gp=180, Gt=180, Ip=20, X_action=0.016, Gc=150, GLP1=5,
beta_mass=0.7, IR_H=2.5, IR_P=2.0, FFA=0.8, BW_t=90,
HbA1c_cmpt=9.0, eGFR_cmpt=80, UACR_cmpt=50

$MAIN
double Cmet   = MET_C  / Vc_met  * 1000;
double Cempa  = EMPA_C / Vc_empa * 1000;
double Csema  = SEMA_C / Vc_sema;
double Cdpp4  = DPP4I_C/ Vc_dpp4 * 1000;
double Csu    = SU_C   / Vc_su   * 1000;
double Cins   = INS_C  / Vc_ins  * 1000;
double Cpiog  = PIOG_C / Vc_piog * 1000;

double Emet_EGP  = USE_MET  * 0.30 * Cmet / (Cmet + 800);
double Emet_Rd   = USE_MET  * 0.12 * Cmet / (Cmet + 800);
double E_empa    = USE_EMPA * Imax_empa * Cempa / (Cempa + IC50_empa);
double TmG       = TmG_base * (1 - E_empa);
double UGE_loss  = fmax(0.0, GFR_val * Gp / 100 - TmG) * 60 / BW_t;
double E_sema_ins= USE_SEMA * Emax_sema_ins * pow(Csema,1.5) / (pow(Csema,1.5) + pow(EC50_sema,1.5));
double E_sema_Gc = USE_SEMA * 0.40 * Csema / (Csema + EC50_sema);
double DPP4_inh  = USE_DPP4 * Imax_dpp4 * Cdpp4 / (Cdpp4 + IC50_dpp4);
double E_su      = USE_SU   * Emax_su * pow(Csu, n_su) / (pow(Csu,n_su) + pow(EC50_su,n_su));
double E_piog_IR = USE_PIOG * Emax_piog_IR * Cpiog / (Cpiog + EC50_piog);

double GLP1_total = GLP1 + USE_SEMA * Csema * 50;
double GSIS_base  = phi_max * beta_mass * beta_sens
                    * pow(fmax(0.0, Gp - Gth), 1.5)
                    / (pow(fmax(0.0, Gp - Gth), 1.5) + pow(90.0, 1.5));
double GSIS_inc   = GSIS_base * (1 + 0.35 * GLP1_total/(GLP1_total + 10.0) + E_sema_ins);
double GSIS_su    = GSIS_inc * (1 + E_su);
double dIp_sec    = GSIS_su;

double EGP_ins_s  = 1.0 / (1.0 + 3.0 * X_action);
double EGP_Gc_d   = 1.0 + EGP_Gc * fmax(0.0, Gc - Gc0);
double EGP_IR_d   = IR_H / 2.5;
double EGP_val    = fmax(0.5, fmin(EGP0 * EGP_ins_s * EGP_Gc_d * EGP_IR_d * (1 - Emet_EGP), 8.0));
double SI_eff     = SI / IR_P;
double Rd_ins     = SI_eff * X_action * Gt;
double Rd_val     = Rd0 + Rd_ins + Emet_Rd * Rd0;
double Gc_ss      = Gc0 * (Gc_Gp50 / fmax(Gp,50.0)) * (1-E_sema_Gc)
                    * (1.0 / (1.0 + 0.02 * fmax(0.0, Ip-Ip0)));

double beta_prol  = k_prolif * beta_mass * fmax(0.0, Gp-90)/90.0;
double beta_apo   = k_apop * beta_mass * (1+0.5*fmax(0.0,FFA-0.6))
                    * (1+0.3*fmax(0.0,Gp-180)/180.0)
                    * (1 - 0.3*GLP1_total/(GLP1_total+10.0));
double Ip_al      = 1.0 / (1.0 + Ip / FFA_Ip50);
double FFA_rel    = kFFA_rel * Ip_al * BW_t / 90.0;
double FFA_up     = kFFA_up * FFA;
double HbA1c_tgt  = 5.0 + HbA1c_ss * fmax(0.0, fmax(Gp,70.0) - 90.0);
double HbA1c_exc  = fmax(0.0, HbA1c_cmpt - 7.0);
double eGFR_dec   = k_eGFR_decline * HbA1c_exc * 24;
double eGFR_prot  = USE_EMPA * k_eGFR_empa * 24 * E_empa;
double UACR_drv   = k_UACR_rise * HbA1c_exc * 24;
double UACR_prot  = USE_EMPA * k_UACR_empa * 24 * E_empa;
double BW_UGE     = UGE_loss / 40.0 * 0.001;
double BW_sema_r  = USE_SEMA * ksema_wt * fmax(0.0, Emax_sema_wt - fmax(0.0,90-BW_t)) / 24.0;
double dBW        = BW_rate - BW_UGE - BW_sema_r;

$ODE
dxdt_MET_GUT = -ka_met * MET_GUT;
dxdt_MET_C   = ka_met * MET_GUT - (CL_met/Vc_met + k12_met) * MET_C + k21_met * MET_P;
dxdt_MET_P   = k12_met * MET_C - k21_met * MET_P;
dxdt_EMPA_C  = -CL_empa/Vc_empa * EMPA_C;
dxdt_SEMA_SC = -ka_sema * SEMA_SC;
dxdt_SEMA_C  = ka_sema * SEMA_SC - CL_sema/Vc_sema * SEMA_C;
dxdt_DPP4I_C = -CL_dpp4/Vc_dpp4 * DPP4I_C;
dxdt_SU_C    = -CL_su/Vc_su * SU_C;
dxdt_INS_SC  = -ka_ins * INS_SC;
dxdt_INS_C   = ka_ins * INS_SC - CL_ins/Vc_ins * INS_C;
dxdt_PIOG_C  = -CL_piog/Vc_piog * PIOG_C;
dxdt_Gp      = (EGP_val - Rd_val) / (Vg * BW_t) * 10 - Sg*(Gp-90) - UGE_loss/(Vg*BW_t);
dxdt_Gt      = Sg*(Gp - Gt) - Rd_ins/(Vg*BW_t)*5;
dxdt_Ip      = (dIp_sec - CL_ins*Ip/Vc_ins) / Vi / BW_t * 10;
dxdt_X_action= -p2*X_action + p2*SI*Ip;
dxdt_Gc      = kout_Gc * (Gc_ss - Gc);
dxdt_GLP1    = kin_GLP1 + GLP1_Gp*fmax(0.0,Gp-90) - (kDPP4*(1-DPP4_inh) + kout_GLP1)*GLP1;
dxdt_beta_mass = (beta_prol - beta_apo) / (24.0*365.25);
dxdt_IR_H    = k_IR_FFA*fmax(0.0,FFA-0.5) - k_IR_rec*(IR_H-1.0) - 0.001*Emet_EGP - 0.003*E_piog_IR*IR_H;
dxdt_IR_P    = k_IR_FFA*fmax(0.0,FFA-0.5) - k_IR_rec*(IR_P-1.0) - 0.003*E_piog_IR*IR_P;
dxdt_FFA     = FFA_rel - FFA_up;
dxdt_BW_t    = dBW;
dxdt_HbA1c_cmpt = kHbA1c * (HbA1c_tgt - HbA1c_cmpt);
dxdt_eGFR_cmpt  = -(eGFR_dec - eGFR_prot);
dxdt_UACR_cmpt  = UACR_drv - UACR_prot;

$CAPTURE
Cmet Cempa Csema Cdpp4 Csu Cins Cpiog
EGP_val Rd_val UGE_loss E_empa E_sema_ins DPP4_inh E_su E_piog_IR
beta_prol beta_apo GLP1 Gc Ip
IR_H IR_P FFA BW_t
HbA1c_cmpt eGFR_cmpt UACR_cmpt
'

mod_shiny <- mcode("T2DM_Shiny", t2dm_code, quiet = TRUE)

# ============================================================================
# Scenario runner
# ============================================================================
run_sim <- function(params_list, drugs, days = 365) {
  events <- list()
  if ("met"  %in% drugs) events <- c(events, list(ev(amt=550,   cmt="MET_GUT", ii=12, addl=2*days-1, time=0)))
  if ("empa" %in% drugs) events <- c(events, list(ev(amt=0.19e6/450.9*8.6, cmt="EMPA_C", ii=24, addl=days-1, time=0)))
  if ("sema" %in% drugs) events <- c(events, list(ev(amt=1000/4114*1000*0.89, cmt="SEMA_SC", ii=168, addl=ceiling(days/7)-1, time=0)))
  if ("dpp4" %in% drugs) events <- c(events, list(ev(amt=87e3/407.5, cmt="DPP4I_C", ii=24, addl=days-1, time=0)))
  if ("su"   %in% drugs) events <- c(events, list(ev(amt=4*1000/490.6, cmt="SU_C", ii=24, addl=days-1, time=0)))
  if ("ins"  %in% drugs) events <- c(events, list(ev(amt=20*0.0347*0.91*1000/6103, cmt="INS_SC", ii=24, addl=days-1, time=0)))
  if ("piog" %in% drugs) events <- c(events, list(ev(amt=30*830/356.5, cmt="PIOG_C", ii=24, addl=days-1, time=0)))

  drug_flags <- list(
    USE_MET  = as.numeric("met"  %in% drugs),
    USE_EMPA = as.numeric("empa" %in% drugs),
    USE_SEMA = as.numeric("sema" %in% drugs),
    USE_DPP4 = as.numeric("dpp4" %in% drugs),
    USE_SU   = as.numeric("su"   %in% drugs),
    USE_INS  = as.numeric("ins"  %in% drugs),
    USE_PIOG = as.numeric("piog" %in% drugs)
  )
  all_params <- c(params_list, drug_flags)

  if (length(events) == 0) events <- list(ev(amt=0, cmt="MET_GUT", time=0))
  dose_ev <- Reduce(c, events)

  mod_shiny %>%
    init(Gp       = params_list$Gp0,
         Gt       = params_list$Gp0,
         HbA1c_cmpt = params_list$HbA1c0,
         eGFR_cmpt  = params_list$eGFR0,
         BW_t     = params_list$BW) %>%
    param(all_params) %>%
    ev(dose_ev) %>%
    mrgsim(end = days * 24, delta = 6, rtol = 1e-4, atol = 1e-6) %>%
    as.data.frame() %>%
    mutate(time_days = time / 24)
}

# ============================================================================
# UI
# ============================================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "T2DM QSP Dashboard", titleWidth = 280),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "tab_patient",  icon = icon("user-md")),
      menuItem("Drug PK",            tabName = "tab_pk",       icon = icon("capsules")),
      menuItem("Glucose & Insulin",  tabName = "tab_gi",       icon = icon("tint")),
      menuItem("β-cell Dynamics",    tabName = "tab_beta",     icon = icon("dna")),
      menuItem("Clinical Endpoints", tabName = "tab_endo",     icon = icon("chart-line")),
      menuItem("Scenario Comparison",tabName = "tab_compare",  icon = icon("balance-scale")),
      menuItem("Biomarkers & Safety",tabName = "tab_bio",      icon = icon("flask"))
    ),
    hr(),
    h5("Patient Parameters", style = "color:white; padding-left:15px"),
    sliderInput("BW",     "Body Weight (kg)", 50, 150, 90, step = 5),
    sliderInput("HbA1c0", "Baseline HbA1c (%)", 7, 13, 9.0, step = 0.1),
    sliderInput("eGFR0",  "Baseline eGFR",      20, 120, 80, step = 5),
    sliderInput("Gp0",    "Baseline Glucose (mg/dL)", 100, 350, 180, step = 10),
    sliderInput("sim_days", "Simulation Duration (days)", 90, 730, 365, step = 30),
    hr(),
    h5("Drug Selection", style = "color:white; padding-left:15px"),
    checkboxGroupInput("drugs", NULL,
                       choices  = c("Metformin 1g BID" = "met",
                                    "Empagliflozin 10mg QD" = "empa",
                                    "Semaglutide 1mg QW" = "sema",
                                    "Sitagliptin 100mg QD" = "dpp4",
                                    "Glimepiride 4mg QD" = "su",
                                    "Insulin Degludec 20U QD" = "ins",
                                    "Pioglitazone 30mg QD" = "piog"),
                       selected = c("met")),
    actionButton("run_btn", "Run Simulation", class = "btn-success btn-block",
                 icon = icon("play"))
  ),

  dashboardBody(
    tags$head(tags$style(HTML(
      ".box { border-top-color: #3c8dbc; }
       .nav-tabs-custom>.tab-content { padding: 8px; }"
    ))),
    tabItems(

      # ----------------------------------------------------------------
      # TAB 1: Patient Profile
      # ----------------------------------------------------------------
      tabItem("tab_patient",
        fluidRow(
          box(title = "Patient Characteristics", width = 4, status = "primary",
              h4("Current Settings"),
              tableOutput("tbl_patient")),
          box(title = "Disease Classification", width = 4, status = "warning",
              h4("Glycemic Staging"),
              htmlOutput("html_staging")),
          box(title = "Cardiovascular Risk Profile", width = 4, status = "danger",
              h4("Risk Factors"),
              tableOutput("tbl_risk"))
        ),
        fluidRow(
          box(title = "Pharmacological Overview", width = 12,
              plotlyOutput("plot_drug_radar", height = "350px"))
        )
      ),

      # ----------------------------------------------------------------
      # TAB 2: Drug PK
      # ----------------------------------------------------------------
      tabItem("tab_pk",
        fluidRow(
          box(title = "Plasma Drug Concentrations", width = 8, status = "primary",
              plotlyOutput("plot_pk", height = "420px")),
          box(title = "PK Summary", width = 4,
              tableOutput("tbl_pk_summary"))
        ),
        fluidRow(
          box(title = "Insulin (Exogenous + Endogenous)", width = 6,
              plotlyOutput("plot_insulin_pk", height = "300px")),
          box(title = "GLP-1 with DPP-4i Effect", width = 6,
              plotlyOutput("plot_glp1_pk", height = "300px"))
        )
      ),

      # ----------------------------------------------------------------
      # TAB 3: Glucose & Insulin Dynamics
      # ----------------------------------------------------------------
      tabItem("tab_gi",
        fluidRow(
          box(title = "Plasma Glucose (Gp)", width = 6, status = "primary",
              plotlyOutput("plot_glucose", height = "300px")),
          box(title = "Plasma Insulin (Ip)", width = 6,
              plotlyOutput("plot_Ip", height = "300px"))
        ),
        fluidRow(
          box(title = "Endogenous Glucose Production (EGP)", width = 6,
              plotlyOutput("plot_EGP", height = "300px")),
          box(title = "Glucose Disposal Rate (Rd)", width = 6,
              plotlyOutput("plot_Rd", height = "300px"))
        ),
        fluidRow(
          box(title = "Glucagon (Gc)", width = 6,
              plotlyOutput("plot_glucagon", height = "300px")),
          box(title = "Insulin Resistance Indices", width = 6,
              plotlyOutput("plot_IR", height = "300px"))
        )
      ),

      # ----------------------------------------------------------------
      # TAB 4: β-cell Dynamics
      # ----------------------------------------------------------------
      tabItem("tab_beta",
        fluidRow(
          box(title = "β-cell Mass (Normalized)", width = 6, status = "warning",
              plotlyOutput("plot_beta_mass", height = "300px")),
          box(title = "Proliferation vs Apoptosis Rate", width = 6,
              plotlyOutput("plot_beta_rates", height = "300px"))
        ),
        fluidRow(
          box(title = "β-cell Mass vs HbA1c Phase Portrait", width = 6,
              plotlyOutput("plot_phase", height = "300px")),
          box(title = "FFA (Lipotoxicity on β-cells)", width = 6,
              plotlyOutput("plot_FFA", height = "300px"))
        )
      ),

      # ----------------------------------------------------------------
      # TAB 5: Clinical Endpoints
      # ----------------------------------------------------------------
      tabItem("tab_endo",
        fluidRow(
          box(title = "HbA1c Trajectory", width = 6, status = "success",
              plotlyOutput("plot_hba1c", height = "300px")),
          box(title = "Body Weight Change", width = 6,
              plotlyOutput("plot_bw", height = "300px"))
        ),
        fluidRow(
          box(title = "eGFR (Renal Protection)", width = 6,
              plotlyOutput("plot_egfr", height = "300px")),
          box(title = "UACR (Albuminuria)", width = 6,
              plotlyOutput("plot_uacr", height = "300px"))
        ),
        fluidRow(
          box(title = "SGLT2i: Urinary Glucose Excretion", width = 6,
              plotlyOutput("plot_uge", height = "300px")),
          box(title = "52-Week Outcomes Table", width = 6,
              DT::dataTableOutput("tbl_52w"))
        )
      ),

      # ----------------------------------------------------------------
      # TAB 6: Scenario Comparison
      # ----------------------------------------------------------------
      tabItem("tab_compare",
        fluidRow(
          box(title = "Multi-Scenario: Select Comparators", width = 12, status = "info",
              fluidRow(
                column(3, checkboxGroupInput("cmp1", "Scenario 1", choices = c("met","empa","sema","dpp4","su","ins","piog"), selected = "met")),
                column(3, checkboxGroupInput("cmp2", "Scenario 2", choices = c("met","empa","sema","dpp4","su","ins","piog"), selected = c("met","empa"))),
                column(3, checkboxGroupInput("cmp3", "Scenario 3", choices = c("met","empa","sema","dpp4","su","ins","piog"), selected = c("met","sema"))),
                column(3, checkboxGroupInput("cmp4", "Scenario 4", choices = c("met","empa","sema","dpp4","su","ins","piog"), selected = c("met","empa","sema")))
              ),
              actionButton("run_compare", "Compare Scenarios", class = "btn-info"))
        ),
        fluidRow(
          box(title = "ΔHbA1c at 52 Weeks", width = 6,
              plotlyOutput("plot_cmp_hba1c", height = "350px")),
          box(title = "ΔBody Weight at 52 Weeks", width = 6,
              plotlyOutput("plot_cmp_bw", height = "350px"))
        ),
        fluidRow(
          box(title = "Scenario Summary Table", width = 12,
              DT::dataTableOutput("tbl_compare"))
        )
      ),

      # ----------------------------------------------------------------
      # TAB 7: Biomarkers & Safety
      # ----------------------------------------------------------------
      tabItem("tab_bio",
        fluidRow(
          box(title = "Active GLP-1 (pmol/L)", width = 6, status = "info",
              plotlyOutput("plot_glp1", height = "300px")),
          box(title = "Plasma FFA (Lipotoxicity)", width = 6,
              plotlyOutput("plot_ffa", height = "300px"))
        ),
        fluidRow(
          box(title = "Empagliflozin: SGLT2 Inhibition Profile", width = 6,
              plotlyOutput("plot_sglt2_inhib", height = "300px")),
          box(title = "Insulin Resistance Trajectory", width = 6,
              plotlyOutput("plot_IR2", height = "300px"))
        ),
        fluidRow(
          box(title = "Safety Monitoring Summary", width = 12,
              status = "warning",
              fluidRow(
                column(6, h4("SGLT2i Considerations"),
                       tags$ul(
                         tags$li("Genital mycotic infections: ~5-10% risk"),
                         tags$li("DKA risk (low, monitor if unwell)"),
                         tags$li("Mild osmotic diuresis — hydration important"),
                         tags$li("eGFR dip (acute, reversible, ~3 mL/min/1.73m²)"),
                         tags$li("Bone fracture risk with canagliflozin (less so empa)")
                       )),
                column(6, h4("GLP-1RA Considerations"),
                       tags$ul(
                         tags$li("GI side effects: nausea/vomiting (transient)"),
                         tags$li("Pancreatitis risk (rare, monitor amylase)"),
                         tags$li("Thyroid C-cell hyperplasia (rodents only)"),
                         tags$li("Injection site reactions"),
                         tags$li("Contraindicated in MEN2 / medullary thyroid Ca")
                       ))
              ))
        )
      )

    )  # end tabItems
  )   # end dashboardBody
)     # end dashboardPage

# ============================================================================
# SERVER
# ============================================================================
server <- function(input, output, session) {

  # Reactive simulation
  sim_data <- eventReactive(input$run_btn, {
    withProgress(message = "Simulating T2DM model...", value = 0, {
      params <- list(
        BW = input$BW, HbA1c0 = input$HbA1c0,
        eGFR0 = input$eGFR0, UACR0 = 50,
        Gp0 = input$Gp0, Ip0 = 20, Gc0 = 150,
        SI = 8e-4, Sg = 0.01, EGP0 = 2.4
      )
      incProgress(0.3)
      out <- run_sim(params, input$drugs %||% character(0), days = input$sim_days)
      incProgress(0.7)
      out
    })
  }, ignoreNULL = FALSE)

  # Comparison simulation
  cmp_data <- eventReactive(input$run_compare, {
    params <- list(BW=input$BW, HbA1c0=input$HbA1c0, eGFR0=input$eGFR0,
                   UACR0=50, Gp0=input$Gp0, Ip0=20, Gc0=150,
                   SI=8e-4, Sg=0.01, EGP0=2.4)
    scen_list <- list(
      list(name=paste(input$cmp1, collapse="+"), drugs=input$cmp1),
      list(name=paste(input$cmp2, collapse="+"), drugs=input$cmp2),
      list(name=paste(input$cmp3, collapse="+"), drugs=input$cmp3),
      list(name=paste(input$cmp4, collapse="+"), drugs=input$cmp4)
    )
    bind_rows(lapply(scen_list, function(s) {
      tryCatch({
        d <- run_sim(params, s$drugs, days = input$sim_days)
        d$scenario <- s$name
        d
      }, error = function(e) NULL)
    }))
  })

  # ---- Tab 1: Patient ----
  output$tbl_patient <- renderTable({
    data.frame(
      Parameter = c("Body Weight","BMI","HbA1c","FPG","eGFR","Active Drugs"),
      Value     = c(paste0(input$BW, " kg"),
                    paste0(round(input$BW / (1.70^2), 1), " kg/m²"),
                    paste0(input$HbA1c0, " %"),
                    paste0(input$Gp0, " mg/dL"),
                    paste0(input$eGFR0, " mL/min/1.73m²"),
                    paste0(length(input$drugs %||% 0), " drug(s)"))
    )
  }, striped = TRUE, bordered = TRUE)

  output$html_staging <- renderUI({
    h <- input$HbA1c0
    stage <- if (h < 7.0) "<span style='color:green'>Well-controlled</span>"
             else if (h < 8.5) "<span style='color:orange'>Suboptimal</span>"
             else "<span style='color:red'>Poorly controlled</span>"
    eGFR <- input$eGFR0
    ckd_stage <- if (eGFR >= 90) "G1" else if (eGFR >= 60) "G2"
                  else if (eGFR >= 45) "G3a" else if (eGFR >= 30) "G3b"
                  else if (eGFR >= 15) "G4" else "G5"
    HTML(paste0(
      "<b>Glycemic Control:</b> ", stage, "<br>",
      "<b>HbA1c:</b> ", h, "% → Estimated avg glucose: ", round(28.7*h - 46.7, 0), " mg/dL<br>",
      "<b>CKD Stage:</b> ", ckd_stage, " (eGFR=", eGFR, ")<br>",
      "<b>ADA Intensification:</b> ",
      if (h > 9) "Triple/Insulin therapy recommended" else if (h > 7.5) "Dual therapy" else "Lifestyle ± Metformin"
    ))
  })

  output$tbl_risk <- renderTable({
    bmi <- round(input$BW / (1.70^2), 1)
    data.frame(
      Factor = c("Obesity","Hyperglycemia","CKD Risk","BMI Category"),
      Status = c(if (bmi > 30) "Obese" else if (bmi > 25) "Overweight" else "Normal",
                 if (input$HbA1c0 > 9) "High" else if (input$HbA1c0 > 7) "Moderate" else "Low",
                 if (input$eGFR0 < 60) "Elevated" else "Normal",
                 paste0(bmi, " kg/m²"))
    )
  }, striped = TRUE)

  output$plot_drug_radar <- renderPlotly({
    drug_data <- data.frame(
      drug  = c("Metformin","Empagliflozin","Semaglutide","Sitagliptin","Glimepiride","Insulin","Pioglitazone"),
      HbA1c = c(1.4, 0.7, 1.4, 0.7, 1.0, 2.0, 0.9),
      Weight= c(-2, -3, -4.5, 0, 1.0, 2.0, 3.0),
      CV_prot=c(1, 3, 3, 1, 0, 1, 1),
      Renal = c(1, 3, 2, 1, 0, 1, 1),
      Hypo  = c(0, 0, 0, 0, 3, 4, 0)
    )
    plot_ly(drug_data, type = "bar",
            x = ~drug, y = ~HbA1c, name = "HbA1c reduction (%)",
            marker = list(color = "steelblue")) %>%
      add_trace(y = ~Weight, name = "Weight effect (kg)", marker = list(color = "coral")) %>%
      layout(title = "Drug Comparison Overview",
             barmode = "group",
             xaxis = list(title = ""), yaxis = list(title = "Effect size"))
  })

  # ---- Tab 2: PK ----
  output$plot_pk <- renderPlotly({
    d <- sim_data()
    p <- plot_ly()
    if ("met"  %in% input$drugs) p <- add_trace(p, data=d, x=~time_days, y=~Cmet,   name="Metformin (ng/mL)",      type="scatter", mode="lines")
    if ("empa" %in% input$drugs) p <- add_trace(p, data=d, x=~time_days, y=~Cempa,  name="Empagliflozin (ng/mL)",  type="scatter", mode="lines")
    if ("sema" %in% input$drugs) p <- add_trace(p, data=d, x=~time_days, y=~Csema,  name="Semaglutide (nmol/L)",   type="scatter", mode="lines")
    if ("dpp4" %in% input$drugs) p <- add_trace(p, data=d, x=~time_days, y=~Cdpp4,  name="Sitagliptin (ng/mL)",    type="scatter", mode="lines")
    if ("su"   %in% input$drugs) p <- add_trace(p, data=d, x=~time_days, y=~Csu,    name="Glimepiride (ng/mL)",    type="scatter", mode="lines")
    if ("ins"  %in% input$drugs) p <- add_trace(p, data=d, x=~time_days, y=~Cins,   name="Ins Degludec (ng/mL)",   type="scatter", mode="lines")
    if ("piog" %in% input$drugs) p <- add_trace(p, data=d, x=~time_days, y=~Cpiog,  name="Pioglitazone (ng/mL)",   type="scatter", mode="lines")
    p %>% layout(title = "Drug Plasma Concentrations",
                 xaxis = list(title = "Time (days)"),
                 yaxis = list(title = "Concentration"),
                 legend = list(orientation = "h"))
  })

  output$plot_insulin_pk <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~Ip, name = "Insulin (mU/L)", type = "scatter", mode = "lines") %>%
      layout(title = "Plasma Insulin", xaxis = list(title = "Days"), yaxis = list(title = "mU/L"))
  })

  output$plot_glp1_pk <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~GLP1, name = "Active GLP-1 (pmol/L)", type = "scatter", mode = "lines",
            line = list(color = "purple")) %>%
      layout(title = "Active GLP-1", xaxis = list(title = "Days"), yaxis = list(title = "pmol/L"))
  })

  output$tbl_pk_summary <- renderTable({
    d <- sim_data()
    rows <- list()
    if ("met"  %in% input$drugs) rows[[length(rows)+1]] <- data.frame(Drug="Metformin",   Cmax=round(max(d$Cmet, na.rm=T),0),   Tmax=round(d$time_days[which.max(d$Cmet)],1), Css=round(mean(tail(d$Cmet,20)),0), Unit="ng/mL")
    if ("empa" %in% input$drugs) rows[[length(rows)+1]] <- data.frame(Drug="Empagliflozin", Cmax=round(max(d$Cempa,na.rm=T),0), Tmax=round(d$time_days[which.max(d$Cempa)],1), Css=round(mean(tail(d$Cempa,20)),0), Unit="ng/mL")
    if ("sema" %in% input$drugs) rows[[length(rows)+1]] <- data.frame(Drug="Semaglutide",  Cmax=round(max(d$Csema,na.rm=T),3),  Tmax=round(d$time_days[which.max(d$Csema)],1), Css=round(mean(tail(d$Csema,20)),3), Unit="nmol/L")
    if (length(rows) == 0) return(data.frame(Message="No drug selected"))
    bind_rows(rows)
  }, striped = TRUE, bordered = TRUE)

  # ---- Tab 3: Glucose & Insulin ----
  make_plotly_ts <- function(d, y_var, title, ylab, ref_line = NULL, color = "steelblue") {
    p <- plot_ly(d, x = ~time_days, y = as.formula(paste0("~",y_var)),
                 type = "scatter", mode = "lines",
                 line = list(color = color)) %>%
      layout(title = title, xaxis = list(title = "Days"), yaxis = list(title = ylab))
    if (!is.null(ref_line)) p <- add_hline(p, y = ref_line, line = list(dash = "dash", color = "red", width = 1))
    p
  }

  output$plot_glucose  <- renderPlotly({ make_plotly_ts(sim_data(), "Gp",      "Plasma Glucose",      "mg/dL", 126, "coral") })
  output$plot_Ip       <- renderPlotly({ make_plotly_ts(sim_data(), "Ip",      "Plasma Insulin",      "mU/L",  NULL, "steelblue") })
  output$plot_EGP      <- renderPlotly({ make_plotly_ts(sim_data(), "EGP_val", "Endogenous Glucose Production", "mg/kg/min", NULL, "#e67e22") })
  output$plot_Rd       <- renderPlotly({ make_plotly_ts(sim_data(), "Rd_val",  "Glucose Disposal (Rd)", "mg/kg/min", NULL, "#27ae60") })
  output$plot_glucagon <- renderPlotly({ make_plotly_ts(sim_data(), "Gc",      "Glucagon",             "pg/mL", NULL, "#8e44ad") })

  output$plot_IR <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~IR_H, name = "Hepatic IR", type = "scatter", mode = "lines") %>%
      add_trace(y = ~IR_P, name = "Peripheral IR") %>%
      layout(title = "Insulin Resistance Index", xaxis = list(title = "Days"), yaxis = list(title = "IR index (1=normal)"))
  })

  # ---- Tab 4: β-cell ----
  output$plot_beta_mass  <- renderPlotly({ make_plotly_ts(sim_data(), "beta_mass", "β-cell Mass", "Relative mass (1=normal)", NULL, "#f39c12") })

  output$plot_beta_rates <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~beta_prol, name = "Proliferation", type = "scatter", mode = "lines") %>%
      add_trace(y = ~beta_apo, name = "Apoptosis", line = list(dash = "dash")) %>%
      layout(title = "β-cell Proliferation vs Apoptosis", xaxis = list(title = "Days"), yaxis = list(title = "Rate (1/year)"))
  })

  output$plot_phase <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~HbA1c_cmpt, y = ~beta_mass, type = "scatter", mode = "lines+markers",
            marker = list(size = 3), color = ~time_days, colors = "Blues") %>%
      layout(title = "β-cell Mass vs HbA1c Phase Portrait",
             xaxis = list(title = "HbA1c (%)"), yaxis = list(title = "β-cell Mass"))
  })

  output$plot_FFA <- renderPlotly({ make_plotly_ts(sim_data(), "FFA", "Plasma FFA (Lipotoxicity)", "mmol/L", NULL, "#e74c3c") })

  # ---- Tab 5: Clinical Endpoints ----
  output$plot_hba1c <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~HbA1c_cmpt, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(41,182,246,0.2)",
            line = list(color = "#0288d1", width = 2)) %>%
      add_hline(y = 7.0, line = list(dash = "dash", color = "green", width = 1)) %>%
      layout(title = "HbA1c Trajectory", xaxis = list(title = "Days"), yaxis = list(title = "HbA1c (%)"))
  })

  output$plot_bw   <- renderPlotly({ make_plotly_ts(sim_data(), "BW_t",     "Body Weight",    "kg", NULL, "#16a085") })
  output$plot_egfr <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~eGFR_cmpt, type = "scatter", mode = "lines",
            line = list(color = "#2ecc71", width = 2)) %>%
      add_hline(y = 60, line = list(dash = "dash", color = "orange")) %>%
      layout(title = "eGFR Trajectory", xaxis = list(title = "Days"), yaxis = list(title = "eGFR (mL/min/1.73m²)"))
  })
  output$plot_uacr <- renderPlotly({ make_plotly_ts(sim_data(), "UACR_cmpt","UACR",          "μg/mg Cr", 30, "#e74c3c") })
  output$plot_uge  <- renderPlotly({
    d <- sim_data() %>% mutate(UGE_daily = UGE_loss * 24)
    make_plotly_ts(d, "UGE_daily", "Urinary Glucose Excretion", "mg/kg/day", NULL, "#3498db")
  })

  output$tbl_52w <- DT::renderDataTable({
    d <- sim_data()
    last <- tail(d, 4) %>% slice(1)
    tibble::tibble(
      Endpoint = c("HbA1c (%)","ΔHbA1c (%)","Fasting Glucose (mg/dL)","Body Weight (kg)","ΔBW (kg)",
                   "eGFR (mL/min/1.73m²)","UACR (μg/mg)","β-cell Mass"),
      Value    = c(round(last$HbA1c_cmpt,1), round(last$HbA1c_cmpt - input$HbA1c0, 1),
                   round(last$Gp,0), round(last$BW_t,1), round(last$BW_t - input$BW, 1),
                   round(last$eGFR_cmpt,1), round(last$UACR_cmpt,0),
                   round(last$beta_mass, 3))
    )
  }, options = list(pageLength = 10, dom = "t"))

  # ---- Tab 6: Comparison ----
  output$plot_cmp_hba1c <- renderPlotly({
    d <- cmp_data()
    req(d, nrow(d) > 0)
    plot_ly(d, x = ~time_days, y = ~HbA1c_cmpt, color = ~scenario,
            type = "scatter", mode = "lines") %>%
      add_hline(y = 7.0, line = list(dash = "dash", color = "gray")) %>%
      layout(title = "HbA1c by Scenario", xaxis = list(title = "Days"), yaxis = list(title = "HbA1c (%)"))
  })

  output$plot_cmp_bw <- renderPlotly({
    d <- cmp_data()
    req(d, nrow(d) > 0)
    plot_ly(d, x = ~time_days, y = ~BW_t, color = ~scenario,
            type = "scatter", mode = "lines") %>%
      layout(title = "Body Weight by Scenario", xaxis = list(title = "Days"), yaxis = list(title = "kg"))
  })

  output$tbl_compare <- DT::renderDataTable({
    d <- cmp_data()
    req(d, nrow(d) > 0)
    d %>%
      group_by(scenario) %>%
      filter(time_days == max(time_days)) %>%
      slice(1) %>%
      summarise(
        `HbA1c (%)` = round(mean(HbA1c_cmpt),1),
        `ΔHbA1c (%)` = round(mean(HbA1c_cmpt) - input$HbA1c0, 1),
        `FPG (mg/dL)` = round(mean(Gp),0),
        `BW (kg)` = round(mean(BW_t),1),
        `ΔBW (kg)` = round(mean(BW_t) - input$BW, 1),
        `eGFR` = round(mean(eGFR_cmpt),1),
        .groups = "drop"
      )
  }, options = list(pageLength = 10))

  # ---- Tab 7: Biomarkers ----
  output$plot_glp1    <- renderPlotly({ make_plotly_ts(sim_data(), "GLP1",       "Active GLP-1",      "pmol/L", NULL, "#9b59b6") })
  output$plot_ffa     <- renderPlotly({ make_plotly_ts(sim_data(), "FFA",        "Plasma FFA",        "mmol/L", NULL, "#e74c3c") })
  output$plot_sglt2_inhib <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~E_empa, type = "scatter", mode = "lines",
            line = list(color = "#1abc9c", width = 2)) %>%
      layout(title = "Empagliflozin SGLT2 Inhibition", xaxis = list(title = "Days"), yaxis = list(title = "Fraction inhibited (0-1)"))
  })
  output$plot_IR2 <- renderPlotly({
    d <- sim_data()
    plot_ly(d, x = ~time_days, y = ~IR_H, name = "Hepatic IR", type = "scatter", mode = "lines") %>%
      add_trace(y = ~IR_P, name = "Peripheral IR") %>%
      layout(title = "Insulin Resistance", xaxis = list(title = "Days"), yaxis = list(title = "IR index"))
  })
}

# Null-coalescing helper
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# ============================================================================
shinyApp(ui = ui, server = server)
