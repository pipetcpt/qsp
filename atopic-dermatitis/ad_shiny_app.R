## =============================================================================
## Atopic Dermatitis QSP — Interactive Shiny Application
## =============================================================================
## Tabs:
##  1. Patient Profile & Disease Severity Configuration
##  2. Drug PK Profiles (Dupilumab / Upadacitinib / Nemolizumab)
##  3. PD Biomarkers (TARC, Eosinophils, IL-31, IgE, STAT6)
##  4. Clinical Endpoints (EASI, IGA, NRS itch)
##  5. Treatment Scenario Comparison
##  6. Epidermal Barrier & Pruritus Deep Dive
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(DT)
library(plotly)

## ---------------------------------------------------------------------------
## Inline mrgsolve model
## ---------------------------------------------------------------------------

ad_model_code <- '
$PARAM @annotated
  KA_DUP    : 0.12  : /day  : SC absorption dupilumab
  CL_DUP    : 0.21  : L/day : Clearance dupilumab
  V1_DUP    : 3.5   : L     : Central volume
  V2_DUP    : 4.2   : L     : Peripheral volume
  Q_DUP     : 0.8   : L/day : Inter-compartmental Q
  F_SC_DUP  : 0.64  :       : SC bioavailability
  MW_DUP    : 148000: g/mol : Dupilumab MW
  Kd_IL4Ra  : 0.027 : nM   : KD for IL-4Ra
  Rtot_IL4Ra: 2.5   : nM   : Total IL-4Ra
  kon_IL4Ra : 0.26  :/nM/d : on-rate
  koff_IL4Ra: 0.007 : /day : off-rate
  kSTAT6_in : 0.5   : /day : pSTAT6 induction
  kSTAT6_out: 0.8   : /day : pSTAT6 elimination
  EC50_STAT6: 0.45  :      : EC50 RO for STAT6 inhibition
  Emax_STAT6: 0.95  :      : Emax STAT6 inhibition
  KA_UPA    : 3.2   : /day : UPA oral absorption
  CL_UPA    : 38.0  : L/day: UPA clearance
  Vd_UPA    : 166   : L    : UPA volume
  F_UPA     : 0.79  :      : UPA bioavailability
  IC50_JAK1 : 0.045 : ng/mL: JAK1 IC50
  Imax_JAK1 : 0.97  :      : JAK1 Imax
  KA_NEMO   : 0.09  : /day : Nemolizumab absorption
  CL_NEMO   : 0.15  : L/day: Nemolizumab clearance
  V1_NEMO   : 3.8   : L    : Nemolizumab volume
  F_SC_NEMO : 0.72  :      : Nemolizumab bioavailability
  IL31_base : 45.0  :pg/mL : Baseline IL-31
  kIL31_prod: 0.18  : /day : IL-31 production
  kIL31_deg : 0.22  : /day : IL-31 degradation
  kTh2_in   : 8.5   : /day : Th2 influx
  kTh2_out  : 0.09  : /day : Th2 efflux
  kTARC_in  : 2.0   : /day : TARC production
  kTARC_out : 0.35  : /day : TARC clearance
  kILC2_in  : 3.0   : /day : ILC2 activation
  kILC2_out : 0.12  : /day : ILC2 deactivation
  kEos_in   : 25.0  :c/uL/d: Eos production
  kEos_out  : 0.07  : /day : Eos clearance
  kIgE_in   : 0.05  : /day : IgE production
  kIgE_out  : 0.005 : /day : IgE catabolism
  kFLG_rest : 0.02  : /day : FLG restoration
  kFLG_sup  : 0.04  : /day : STAT6 FLG suppression
  kInfl_prod: 0.8   : /day : Inflammation production
  kInfl_elim: 0.15  : /day : Inflammation resolution
  kEASI_resp: 0.10  : /day : EASI response rate
  EASI_base : 29.0  :      : Baseline EASI
  EASI_min  : 0.5   :      : Minimum EASI
  TEWL_base : 25.0  :g/m2/h: Baseline TEWL
  TEWL_norm : 8.0   :g/m2/h: Normal TEWL
  TCS_effect: 0.0   :      : TCS additional STAT6 inh.

$CMT DUP_SC DUP_C DUP_P Rfree RC pSTAT6 UPA_GI UPA_C NEMO_SC NEMO_C
     IL31 Th2 TARC ILC2 Eos_blood IgE FLG TEWL SkinInfl EASI

$INIT DUP_SC=0 DUP_C=0 DUP_P=0 Rfree=2.5 RC=0 pSTAT6=100
      UPA_GI=0 UPA_C=0 NEMO_SC=0 NEMO_C=0 IL31=45 Th2=100
      TARC=3500 ILC2=50 Eos_blood=350 IgE=2000 FLG=0.40
      TEWL=25 SkinInfl=100 EASI=29

$ODE
  dxdt_DUP_SC = -KA_DUP * DUP_SC;
  double Cp_dup = DUP_C / V1_DUP;
  dxdt_DUP_C = KA_DUP * DUP_SC * F_SC_DUP
               - (CL_DUP/V1_DUP)*DUP_C - (Q_DUP/V1_DUP)*DUP_C
               + (Q_DUP/V2_DUP)*DUP_P;
  dxdt_DUP_P = (Q_DUP/V1_DUP)*DUP_C - (Q_DUP/V2_DUP)*DUP_P;
  double Cp_dup_nM = (Cp_dup * 1000.0)/MW_DUP * 1e6;
  dxdt_Rfree = koff_IL4Ra*RC - kon_IL4Ra*Rfree*Cp_dup_nM;
  dxdt_RC    = kon_IL4Ra*Rfree*Cp_dup_nM - koff_IL4Ra*RC;
  double RO  = RC/(Rtot_IL4Ra + 1e-10);
  double Cp_upa_ng = (UPA_C/Vd_UPA)*1000.0;
  double JAK1_inh = Imax_JAK1*Cp_upa_ng/(IC50_JAK1+Cp_upa_ng+1e-10);
  double STAT6_inh = 1.0-(1.0-Emax_STAT6*RO/(EC50_STAT6+RO+1e-10))
                     *(1.0-JAK1_inh)*(1.0-TCS_effect);
  if(STAT6_inh>0.99) STAT6_inh=0.99;
  dxdt_pSTAT6 = kSTAT6_in*(1.0-STAT6_inh)*100.0 - kSTAT6_out*pSTAT6;
  dxdt_UPA_GI = -KA_UPA * UPA_GI;
  dxdt_UPA_C  = KA_UPA*F_UPA*UPA_GI - (CL_UPA/Vd_UPA)*UPA_C;
  dxdt_NEMO_SC = -KA_NEMO * NEMO_SC;
  double Cp_nemo = NEMO_C/V1_NEMO;
  dxdt_NEMO_C  = KA_NEMO*F_SC_NEMO*NEMO_SC - (CL_NEMO/V1_NEMO)*NEMO_C;
  double RO_NEMO = Cp_nemo/(Cp_nemo + 0.015);
  double IL31_block = RO_NEMO*0.85 + JAK1_inh*0.60;
  if(IL31_block>0.95) IL31_block=0.95;
  double ILC2_dr = ILC2/50.0; double Th2_dr = Th2/100.0;
  dxdt_IL31 = kIL31_prod*(Th2_dr*0.7+ILC2_dr*0.3)*45.0*(1.0-IL31_block) - kIL31_deg*IL31;
  dxdt_ILC2 = kILC2_in*(1.0-JAK1_inh*0.5)*50.0/(1.0+ILC2/50.0) - kILC2_out*ILC2;
  dxdt_Th2  = kTh2_in*(TARC/3500.0)*(1.0-STAT6_inh*0.7) - kTh2_out*Th2;
  dxdt_TARC = kTARC_in*(pSTAT6/100.0)*3500.0 - kTARC_out*TARC;
  double IL5eq = Th2_dr*ILC2_dr;
  dxdt_Eos_blood = kEos_in*IL5eq*(1.0-STAT6_inh*0.8) - kEos_out*Eos_blood;
  dxdt_IgE  = kIgE_in*(Th2_dr)*(1.0-STAT6_inh*0.6)*2000.0 - kIgE_out*IgE;
  dxdt_FLG  = kFLG_rest*(1.0-FLG) - kFLG_sup*(pSTAT6/100.0)*FLG;
  double TEWL_tgt = 8.0+(25.0-8.0)*(1.0-FLG)/0.60;
  dxdt_TEWL = 0.15*(TEWL_tgt-TEWL);
  double Infl_dr = Th2_dr*(pSTAT6/100.0)*(Eos_blood/350.0);
  dxdt_SkinInfl = kInfl_prod*Infl_dr*100.0 - kInfl_elim*SkinInfl;
  dxdt_EASI = kEASI_resp*(EASI_min+(EASI_base-EASI_min)*(SkinInfl/100.0)-EASI);

$TABLE
  double RO_pct = RC/(Rtot_IL4Ra+1e-10)*100.0;
  double JAK1pct = Imax_JAK1*(UPA_C/Vd_UPA*1000.0)/(IC50_JAK1+UPA_C/Vd_UPA*1000.0+1e-10)*100.0;
  double itch_raw = 10.0*(IL31/45.0)*(1.0-NEMO_C/(NEMO_C+0.015*3.8)*0.85);
  double NRS = itch_raw*(1.0-JAK1pct/100.0*0.7);
  if(NRS<0.2) NRS=0.2; if(NRS>10.0) NRS=10.0;
  double IGA = (EASI>21)?4:(EASI>14)?3:(EASI>7)?2:(EASI>2)?1:0;
  double EASI75 = (EASI<=EASI_base*0.25)?1.0:0.0;
  capture RO_pct, JAK1pct, NRS, IGA, EASI75, TARC, Eos_blood, IgE, IL31, FLG, TEWL
$CAPTURE EASI NRS IGA RO_pct JAK1pct
'

ad_mod <<- mcode("ad_shiny", ad_model_code, quiet = TRUE)

## ---------------------------------------------------------------------------
## Helper: run single simulation
## ---------------------------------------------------------------------------

run_sim <- function(params_list, dosing_ev, sim_days = 364) {
  base_params <- list(
    KA_DUP=0.12, CL_DUP=0.21, V1_DUP=3.5, V2_DUP=4.2, Q_DUP=0.8, F_SC_DUP=0.64,
    MW_DUP=148000, Kd_IL4Ra=0.027, Rtot_IL4Ra=2.5, kon_IL4Ra=0.26, koff_IL4Ra=0.007,
    kSTAT6_in=0.5, kSTAT6_out=0.8, EC50_STAT6=0.45, Emax_STAT6=0.95,
    KA_UPA=3.2, CL_UPA=38.0, Vd_UPA=166, F_UPA=0.79, IC50_JAK1=0.045, Imax_JAK1=0.97,
    KA_NEMO=0.09, CL_NEMO=0.15, V1_NEMO=3.8, F_SC_NEMO=0.72,
    IL31_base=45, kIL31_prod=0.18, kIL31_deg=0.22,
    kTh2_in=8.5, kTh2_out=0.09, kTARC_in=2.0, kTARC_out=0.35,
    kILC2_in=3.0, kILC2_out=0.12, kEos_in=25.0, kEos_out=0.07,
    kIgE_in=0.05, kIgE_out=0.005, kFLG_rest=0.02, kFLG_sup=0.04,
    kInfl_prod=0.8, kInfl_elim=0.15, kEASI_resp=0.10,
    EASI_base=29, EASI_min=0.5, TEWL_base=25, TEWL_norm=8, TCS_effect=0
  )
  final_params <- modifyList(base_params, params_list)
  m <- do.call(param, c(list(ad_mod), final_params))
  if (!is.null(dosing_ev)) {
    out <- m %>% ev(dosing_ev) %>% mrgsim(delta=1, end=sim_days)
  } else {
    out <- m %>% mrgsim(delta=1, end=sim_days)
  }
  as_tibble(out)
}

## ---------------------------------------------------------------------------
## UI
## ---------------------------------------------------------------------------

ui <- dashboardPage(
  skin = "purple",

  dashboardHeader(
    title = "AD QSP Model",
    titleWidth = 260
  ),

  dashboardSidebar(
    width = 260,
    sidebarMenu(
      menuItem("Patient Profile",      tabName = "tab_patient",  icon = icon("user")),
      menuItem("Drug PK Profiles",     tabName = "tab_pk",       icon = icon("pills")),
      menuItem("PD Biomarkers",        tabName = "tab_pd",       icon = icon("vial")),
      menuItem("Clinical Endpoints",   tabName = "tab_endpoints",icon = icon("chart-line")),
      menuItem("Scenario Comparison",  tabName = "tab_compare",  icon = icon("layer-group")),
      menuItem("Barrier & Pruritus",   tabName = "tab_barrier",  icon = icon("allergies"))
    ),

    hr(),

    ## Global Treatment Selection
    h5("Treatment Selection", style="margin-left:15px; color:#ddd;"),

    checkboxInput("use_dupilumab", "Dupilumab (Q2W SC)", value = TRUE),
    conditionalPanel("input.use_dupilumab",
      sliderInput("dup_load", "Loading Dose (mg):", 300, 900, 600, 50),
      sliderInput("dup_maint","Maintenance (mg Q2W):", 100, 600, 300, 50)
    ),

    checkboxInput("use_upadacitinib", "Upadacitinib (QD oral)", value = FALSE),
    conditionalPanel("input.use_upadacitinib",
      radioButtons("upa_dose", "Dose:", c("15mg" = 15, "30mg" = 30), selected = 30, inline = TRUE)
    ),

    checkboxInput("use_nemolizumab", "Nemolizumab (Q4W SC)", value = FALSE),

    checkboxInput("use_tcs", "Topical CS (TCS)", value = FALSE),
    conditionalPanel("input.use_tcs",
      sliderInput("tcs_potency", "TCS Potency (STAT6 supp.):", 0.1, 0.6, 0.35, 0.05)
    ),

    sliderInput("sim_weeks", "Simulation Duration (wk):", 16, 52, 52, 4),
    actionButton("run_sim", "Run Simulation", icon = icon("play"),
                 class = "btn-success", style = "width:90%; margin:5px 5% 0 5%;")
  ),

  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background: #f4f6f9; }
        .box { border-top: 3px solid #6c3483; }
        .box-header .box-title { font-size: 15px; }
        .ggplot-output { background: white; border-radius: 6px; padding: 8px; }
      "))
    ),

    tabItems(

      ## =====================================================================
      ## TAB 1: Patient Profile
      ## =====================================================================
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title="Disease Severity Configuration", width=6, status="primary", solidHeader=TRUE,
            sliderInput("EASI_baseline", "Baseline EASI Score:", 7, 72, 29, 1),
            selectInput("disease_subtype", "AD Subtype:",
                        c("Western AD (Th2-dominant)" = "western",
                          "Asian AD (Th2+Th17 mixed)" = "asian",
                          "Intrinsic AD (IgE normal)" = "intrinsic")),
            sliderInput("IgE_baseline", "Baseline IgE (IU/mL):", 100, 10000, 2000, 100),
            sliderInput("Eos_baseline", "Baseline Eosinophil (/μL):", 100, 2000, 350, 50),
            sliderInput("TARC_baseline", "Baseline TARC (pg/mL):", 500, 15000, 3500, 100)
          ),

          box(title="Patient Demographics", width=6, status="info", solidHeader=TRUE,
            sliderInput("patient_age", "Age (years):", 18, 80, 35, 1),
            radioButtons("patient_sex", "Sex:", c("Male", "Female"), inline = TRUE),
            sliderInput("bsa_affected", "BSA Affected (%):", 5, 90, 40, 5),
            checkboxGroupInput("comorbidities", "Comorbidities:",
                               c("Allergic Asthma", "Allergic Rhinitis",
                                 "Food Allergy", "Anxiety/Depression")),
            hr(),
            h5("Predicted IGA at Baseline:"),
            verbatimTextOutput("baseline_IGA")
          )
        ),

        fluidRow(
          box(title="AD Severity Classification Reference", width=12, status="warning",
            DTOutput("severity_table")
          )
        )
      ),

      ## =====================================================================
      ## TAB 2: Drug PK Profiles
      ## =====================================================================
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title="Dupilumab PK — Plasma Concentration", width=6, status="primary", solidHeader=TRUE,
            plotlyOutput("plot_dup_pk", height = "320px")
          ),
          box(title="IL-4Rα Receptor Occupancy (%)", width=6, status="primary", solidHeader=TRUE,
            plotlyOutput("plot_ro", height = "320px")
          )
        ),
        fluidRow(
          box(title="Upadacitinib PK — Plasma Concentration", width=6, status="warning", solidHeader=TRUE,
            plotlyOutput("plot_upa_pk", height = "320px")
          ),
          box(title="JAK1 Inhibition (%)", width=6, status="warning", solidHeader=TRUE,
            plotlyOutput("plot_jak1", height = "320px")
          )
        ),
        fluidRow(
          box(title="PK Summary Table", width=12, status="info",
            DTOutput("pk_summary_table")
          )
        )
      ),

      ## =====================================================================
      ## TAB 3: PD Biomarkers
      ## =====================================================================
      tabItem(tabName = "tab_pd",
        fluidRow(
          box(title="TARC/CCL17 (Th2 Biomarker)", width=6, status="purple", solidHeader=TRUE,
            plotlyOutput("plot_tarc", height = "300px")
          ),
          box(title="Blood Eosinophil Count", width=6, status="success", solidHeader=TRUE,
            plotlyOutput("plot_eos", height = "300px")
          )
        ),
        fluidRow(
          box(title="IL-31 (Pruritus Cytokine)", width=6, status="warning", solidHeader=TRUE,
            plotlyOutput("plot_il31", height = "300px")
          ),
          box(title="Total IgE Level", width=6, status="danger", solidHeader=TRUE,
            plotlyOutput("plot_ige", height = "300px")
          )
        ),
        fluidRow(
          box(title="pSTAT6 Phosphorylation (%)", width=12, status="info",
            plotlyOutput("plot_stat6", height = "280px")
          )
        )
      ),

      ## =====================================================================
      ## TAB 4: Clinical Endpoints
      ## =====================================================================
      tabItem(tabName = "tab_endpoints",
        fluidRow(
          box(title="EASI Score Over Time", width=8, status="primary", solidHeader=TRUE,
            plotlyOutput("plot_easi", height = "350px")
          ),
          box(title="EASI Response Rates", width=4, status="success", solidHeader=TRUE,
            plotlyOutput("plot_easi_resp", height = "350px")
          )
        ),
        fluidRow(
          box(title="NRS Itch Score (0-10)", width=6, status="warning", solidHeader=TRUE,
            plotlyOutput("plot_nrs", height = "300px")
          ),
          box(title="IGA Score Over Time", width=6, status="info", solidHeader=TRUE,
            plotlyOutput("plot_iga", height = "300px")
          )
        ),
        fluidRow(
          box(title="Simulated Response Endpoints", width=12,
            DTOutput("endpoint_table")
          )
        )
      ),

      ## =====================================================================
      ## TAB 5: Scenario Comparison
      ## =====================================================================
      tabItem(tabName = "tab_compare",
        fluidRow(
          box(title="Scenario Comparison Settings", width=4, status="primary",
            checkboxGroupInput("compare_scenarios",
                               "Scenarios to Compare:",
                               c("No Treatment"       = "none",
                                 "TCS Only"           = "tcs",
                                 "Dupilumab Q2W"      = "dup",
                                 "Upadacitinib 30mg"  = "upa",
                                 "Nemolizumab Q4W"    = "nemo",
                                 "Dupilumab + TCS"    = "dup_tcs"),
                               selected = c("none", "dup", "upa")),
            actionButton("run_compare", "Run Comparison", icon = icon("play"),
                         class = "btn-primary", style = "width:100%;")
          ),
          box(title="EASI Score Comparison", width=8, status="primary", solidHeader=TRUE,
            plotlyOutput("plot_compare_easi", height = "380px")
          )
        ),
        fluidRow(
          box(title="NRS Itch Comparison", width=6, status="warning", solidHeader=TRUE,
            plotlyOutput("plot_compare_nrs", height = "320px")
          ),
          box(title="TARC Reduction Comparison", width=6, status="info", solidHeader=TRUE,
            plotlyOutput("plot_compare_tarc", height = "320px")
          )
        ),
        fluidRow(
          box(title="Response Endpoint Summary (Week 16 & 52)", width=12,
            DTOutput("compare_summary_table")
          )
        )
      ),

      ## =====================================================================
      ## TAB 6: Barrier & Pruritus
      ## =====================================================================
      tabItem(tabName = "tab_barrier",
        fluidRow(
          box(title="Filaggrin (FLG) Expression Recovery", width=6, status="success", solidHeader=TRUE,
            plotlyOutput("plot_flg", height = "300px")
          ),
          box(title="Trans-Epidermal Water Loss (TEWL)", width=6, status="warning", solidHeader=TRUE,
            plotlyOutput("plot_tewl", height = "300px")
          )
        ),
        fluidRow(
          box(title="Itch-Inflammation Feedback Loop", width=8, status="danger", solidHeader=TRUE,
            plotlyOutput("plot_itch_infl", height = "350px")
          ),
          box(title="Barrier Metrics Summary", width=4, status="info",
            h5("Barrier Recovery at Weeks:"),
            DTOutput("barrier_table")
          )
        ),
        fluidRow(
          box(title="Mechanistic Summary", width=12, status="primary",
            HTML("
              <div style='font-size:14px; line-height:1.8;'>
              <b>Epidermal Barrier in Atopic Dermatitis:</b><br>
              • FLG mutations present in ~30% European, ~15% Asian AD patients<br>
              • STAT6 activation (IL-4/IL-13) suppresses FLG transcription even without genetic mutation<br>
              • Dupilumab restores FLG expression via pSTAT6 inhibition within 4-8 weeks<br>
              • TEWL normalization lags FLG recovery by ~2-4 weeks due to barrier reconstruction<br>
              • IL-31 pathway (Nemolizumab target): direct neural itch signaling via IL-31Ra on DRG neurons<br>
              • JAK inhibitors (Upadacitinib) block TSLP, IL-4, IL-13, IL-31 simultaneously for rapid itch relief<br>
              • Emollients provide adjunct barrier support independent of immunological pathways
              </div>
            ")
          )
        )
      )
    )
  )
)

## ---------------------------------------------------------------------------
## Server
## ---------------------------------------------------------------------------

server <- function(input, output, session) {

  ## -----------------------------------------------------------------------
  ## Reactive: Build Dosing Events
  ## -----------------------------------------------------------------------
  build_dosing <- function(dup = FALSE, upa = FALSE, nemo = FALSE,
                            dup_load = 600, dup_maint = 300,
                            upa_dose = 30, tcs_eff = 0,
                            sim_wk = 52) {
    events_list <- list()
    sim_days <- sim_wk * 7

    if (dup) {
      ev_load <- ev(amt = dup_load, cmt = "DUP_SC", time = 0)
      n_maint <- floor(sim_days / 14)
      ev_m <- ev(amt = dup_maint, cmt = "DUP_SC",
                 time = seq(14, by = 14, length.out = n_maint))
      events_list <- c(events_list, list(ev_load, ev_m))
    }
    if (upa) {
      ev_upa <- ev(amt = as.numeric(upa_dose), cmt = "UPA_GI",
                   time = seq(0, by = 1, length.out = sim_days))
      events_list <- c(events_list, list(ev_upa))
    }
    if (nemo) {
      n_nemo <- ceiling(sim_days / 28)
      ev_n <- ev(amt = 60, cmt = "NEMO_SC",
                 time = seq(0, by = 28, length.out = n_nemo))
      events_list <- c(events_list, list(ev_n))
    }

    if (length(events_list) == 0) return(NULL)
    do.call(c, events_list)
  }

  ## -----------------------------------------------------------------------
  ## Reactive: Run Primary Simulation
  ## -----------------------------------------------------------------------
  sim_result <- eventReactive(input$run_sim, {
    withProgress(message = "Running simulation...", value = 0.3, {
      ev_obj <- build_dosing(
        dup      = input$use_dupilumab,
        upa      = input$use_upadacitinib,
        nemo     = input$use_nemolizumab,
        dup_load = input$dup_load,
        dup_maint = input$dup_maint,
        upa_dose  = input$upa_dose,
        tcs_eff   = if (input$use_tcs) input$tcs_potency else 0,
        sim_wk   = input$sim_weeks
      )
      params_extra <- list(TCS_effect = if (input$use_tcs) input$tcs_potency else 0)
      setProgress(0.6)
      result <- run_sim(params_extra, ev_obj, sim_days = input$sim_weeks * 7)
      setProgress(1.0)
      result
    })
  }, ignoreNULL = FALSE)

  # Run with default on start
  observe({ sim_result() })

  ## -----------------------------------------------------------------------
  ## TAB 1 OUTPUTS
  ## -----------------------------------------------------------------------
  output$baseline_IGA <- renderText({
    easi <- input$EASI_baseline
    iga <- if (easi > 21) "4 — Severe"
           else if (easi > 14) "3 — Moderate"
           else if (easi > 7) "2 — Mild"
           else if (easi > 2) "1 — Almost Clear"
           else "0 — Clear"
    paste("IGA:", iga)
  })

  output$severity_table <- renderDT({
    df <- data.frame(
      IGA   = c("0 — Clear", "1 — Almost Clear", "2 — Mild",
                "3 — Moderate", "4 — Severe"),
      EASI  = c("0", "1-6", "7-14", "15-21", ">21"),
      SCORAD = c("0", "<15", "15-40", "40-60", ">60"),
      BSA   = c("<1%", "1-5%", "6-15%", "16-30%", ">30%"),
      NRS   = c("0", "1-3", "3-5", "5-7", ">7")
    )
    datatable(df, options = list(dom = "t", pageLength = 5),
              rownames = FALSE,
              caption = "AD Severity Classification (IGA, EASI, SCORAD, BSA, NRS)")
  })

  ## -----------------------------------------------------------------------
  ## TAB 2 — PK PLOTS
  ## -----------------------------------------------------------------------
  output$plot_dup_pk <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>%
      mutate(Cp_dup = DUP_C / 3.5, week = time / 7)
    plot_ly(df, x = ~week, y = ~Cp_dup, type = "scatter", mode = "lines",
            line = list(color = "#2980b9", width = 2)) %>%
      layout(title = "Dupilumab Cp (μg/mL)",
             xaxis = list(title = "Week"),
             yaxis = list(title = "Cp (μg/mL)"),
             shapes = list(list(type="line", x0=0, x1=max(df$week),
                               y0=2, y1=2, line=list(dash="dash", color="green"))))
  })

  output$plot_ro <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>% mutate(week = time / 7)
    plot_ly(df, x = ~week, y = ~RO_pct, type = "scatter", mode = "lines",
            line = list(color = "#27ae60", width = 2)) %>%
      layout(title = "IL-4Rα Receptor Occupancy",
             xaxis = list(title = "Week"),
             yaxis = list(title = "RO (%)", range = c(0, 100)),
             shapes = list(
               list(type="line", x0=0, x1=max(df$week), y0=70, y1=70,
                    line=list(dash="dash", color="orange")),
               list(type="line", x0=0, x1=max(df$week), y0=90, y1=90,
                    line=list(dash="dot", color="green"))))
  })

  output$plot_upa_pk <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>%
      mutate(Cp_upa = UPA_C / 166 * 1000, week = time / 7)
    plot_ly(df, x = ~week, y = ~Cp_upa, type = "scatter", mode = "lines",
            line = list(color = "#8e44ad", width = 2)) %>%
      layout(title = "Upadacitinib Cp (ng/mL)",
             xaxis = list(title = "Week"),
             yaxis = list(title = "Cp (ng/mL)"))
  })

  output$plot_jak1 <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>% mutate(week = time / 7)
    plot_ly(df, x = ~week, y = ~JAK1pct, type = "scatter", mode = "lines",
            line = list(color = "#7d3c98", width = 2)) %>%
      layout(title = "JAK1 Inhibition (%)",
             xaxis = list(title = "Week"),
             yaxis = list(title = "JAK1 Inhibition (%)", range = c(0, 100)))
  })

  output$pk_summary_table <- renderDT({
    req(sim_result())
    df <- sim_result() %>%
      filter(time %in% c(14, 28, 56, 112, 224, 364)) %>%
      mutate(Week = time / 7,
             `Dup Cp (μg/mL)` = round(DUP_C / 3.5, 3),
             `RO (%)` = round(RO_pct, 1),
             `UPA Cp (ng/mL)` = round(UPA_C / 166 * 1000, 2),
             `JAK1 Inh (%)` = round(JAK1pct, 1)) %>%
      select(Week, `Dup Cp (μg/mL)`, `RO (%)`, `UPA Cp (ng/mL)`, `JAK1 Inh (%)`)
    datatable(df, options = list(dom = "t"), rownames = FALSE)
  })

  ## -----------------------------------------------------------------------
  ## TAB 3 — PD BIOMARKERS
  ## -----------------------------------------------------------------------
  make_plotly_line <- function(df, x, y, title, ylab, color,
                               hline = NULL, hline_color = "gray") {
    p <- plot_ly(df, x = as.formula(paste0("~", x)),
                 y = as.formula(paste0("~", y)),
                 type = "scatter", mode = "lines",
                 line = list(color = color, width = 2)) %>%
      layout(title = title,
             xaxis = list(title = "Week"),
             yaxis = list(title = ylab))
    if (!is.null(hline)) {
      p <- p %>% layout(shapes = list(
        list(type="line", x0=0, x1=max(df[[x]]),
             y0=hline, y1=hline, line=list(dash="dash", color=hline_color))))
    }
    p
  }

  output$plot_tarc <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>% mutate(week = time / 7)
    make_plotly_line(df, "week", "TARC", "TARC/CCL17 (pg/mL)",
                     "TARC (pg/mL)", "#9b59b6", hline = 450, hline_color = "green")
  })

  output$plot_eos <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>% mutate(week = time / 7)
    make_plotly_line(df, "week", "Eos_blood", "Blood Eosinophil Count",
                     "Eosinophils (/μL)", "#27ae60", hline = 500, hline_color = "orange")
  })

  output$plot_il31 <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>% mutate(week = time / 7)
    make_plotly_line(df, "week", "IL31", "IL-31 Level (pg/mL)",
                     "IL-31 (pg/mL)", "#f39c12", hline = 20, hline_color = "green")
  })

  output$plot_ige <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>% mutate(week = time / 7)
    make_plotly_line(df, "week", "IgE", "Total IgE (IU/mL)",
                     "IgE (IU/mL)", "#e74c3c", hline = 200, hline_color = "darkgreen")
  })

  output$plot_stat6 <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>% mutate(week = time / 7)
    make_plotly_line(df, "week", "pSTAT6", "pSTAT6 (AU)",
                     "pSTAT6 (AU)", "#1a5276", hline = 50, hline_color = "orange")
  })

  ## -----------------------------------------------------------------------
  ## TAB 4 — CLINICAL ENDPOINTS
  ## -----------------------------------------------------------------------
  output$plot_easi <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>% mutate(week = time / 7)
    plot_ly(df, x = ~week, y = ~EASI, type = "scatter", mode = "lines",
            line = list(color = "#6c3483", width = 2.5)) %>%
      layout(title = "EASI Score", xaxis = list(title = "Week"),
             yaxis = list(title = "EASI (0-72)", range = c(0, max(df$EASI) + 3)),
             shapes = list(
               list(type="line", x0=0, x1=max(df$week), y0=7, y1=7,
                    line=list(dash="dash", color="green")),
               list(type="line", x0=0, x1=max(df$week), y0=14, y1=14,
                    line=list(dash="dash", color="orange")),
               list(type="line", x0=0, x1=max(df$week), y0=21, y1=21,
                    line=list(dash="dash", color="red"))))
  })

  output$plot_easi_resp <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>%
      filter(time %in% c(0, 28, 56, 84, 112, 168, 224, 280, 336, 364)) %>%
      mutate(week = time / 7,
             EASI75 = EASI <= input$EASI_baseline * 0.25,
             EASI90 = EASI <= input$EASI_baseline * 0.10,
             IGA01  = IGA <= 1)
    plot_ly(df, x = ~week) %>%
      add_bars(y = ~as.numeric(EASI75) * 100, name = "EASI-75",
               marker = list(color = "#27ae60")) %>%
      add_bars(y = ~as.numeric(EASI90) * 100, name = "EASI-90",
               marker = list(color = "#2980b9")) %>%
      add_bars(y = ~as.numeric(IGA01) * 100, name = "IGA 0/1",
               marker = list(color = "#8e44ad")) %>%
      layout(title = "Response (%)", xaxis = list(title = "Week"),
             yaxis = list(title = "Response Rate (%)", range = c(0, 105)),
             barmode = "group")
  })

  output$plot_nrs <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>% mutate(week = time / 7)
    make_plotly_line(df, "week", "NRS", "Pruritus NRS (0-10)",
                     "NRS Score", "#e67e22", hline = 4, hline_color = "blue")
  })

  output$plot_iga <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>% mutate(week = time / 7)
    make_plotly_line(df, "week", "IGA", "IGA Score (0-4)",
                     "IGA", "#2c3e50")
  })

  output$endpoint_table <- renderDT({
    req(sim_result())
    df <- sim_result() %>%
      filter(time %in% c(0, 56, 112, 168, 224, 280, 336, 364)) %>%
      mutate(Week = time / 7,
             EASI_r  = round(EASI, 1),
             NRS_r   = round(NRS, 1),
             IGA_r   = round(IGA, 0),
             TARC_r  = round(TARC, 0),
             Eos_r   = round(Eos_blood, 0),
             IgE_r   = round(IgE, 0),
             EASI75  = ifelse(EASI <= input$EASI_baseline * 0.25, "Yes", "No"),
             IGA01   = ifelse(IGA <= 1, "Yes", "No")) %>%
      select(Week, EASI_r, NRS_r, IGA_r, EASI75, IGA01, TARC_r, Eos_r, IgE_r)
    names(df) <- c("Week", "EASI", "NRS", "IGA", "EASI-75?", "IGA0/1?",
                   "TARC(pg/mL)", "Eos(/μL)", "IgE(IU/mL)")
    datatable(df, options = list(dom = "t", pageLength = 10),
              rownames = FALSE)
  })

  ## -----------------------------------------------------------------------
  ## TAB 5 — SCENARIO COMPARISON
  ## -----------------------------------------------------------------------
  compare_data <- eventReactive(input$run_compare, {
    withProgress(message = "Comparing scenarios...", value = 0.1, {
      scenarios <- input$compare_scenarios
      map_dfr(scenarios, function(sc) {
        dup_on  <- sc %in% c("dup", "dup_tcs")
        upa_on  <- sc == "upa"
        nemo_on <- sc == "nemo"
        tcs_eff <- if (sc %in% c("tcs", "dup_tcs")) 0.35 else 0.0

        ev_obj <- build_dosing(dup=dup_on, upa=upa_on, nemo=nemo_on,
                               dup_load=600, dup_maint=300,
                               upa_dose=30, tcs_eff=tcs_eff,
                               sim_wk=52)
        res <- run_sim(list(TCS_effect=tcs_eff), ev_obj, sim_days=364)
        label <- switch(sc,
          "none" = "No Treatment",
          "tcs"  = "TCS Only",
          "dup"  = "Dupilumab Q2W",
          "upa"  = "Upadacitinib 30mg",
          "nemo" = "Nemolizumab Q4W",
          "dup_tcs" = "Dupilumab + TCS"
        )
        res %>% mutate(scenario = label)
      })
    })
  })

  sc_colors <- c(
    "No Treatment"      = "#e74c3c",
    "TCS Only"          = "#f39c12",
    "Dupilumab Q2W"     = "#2980b9",
    "Upadacitinib 30mg" = "#8e44ad",
    "Nemolizumab Q4W"   = "#27ae60",
    "Dupilumab + TCS"   = "#1a5276"
  )

  output$plot_compare_easi <- renderPlotly({
    req(compare_data())
    df <- compare_data() %>% mutate(week = time / 7)
    p <- plot_ly()
    for (sc in unique(df$scenario)) {
      sub <- df %>% filter(scenario == sc)
      p <- add_trace(p, x=sub$week, y=sub$EASI, type="scatter", mode="lines",
                     name=sc, line=list(width=2, color=sc_colors[[sc]]))
    }
    p %>% layout(title="EASI Scenario Comparison",
                 xaxis=list(title="Week"),
                 yaxis=list(title="EASI"),
                 legend=list(orientation="h", y=-0.25))
  })

  output$plot_compare_nrs <- renderPlotly({
    req(compare_data())
    df <- compare_data() %>% mutate(week = time / 7)
    p <- plot_ly()
    for (sc in unique(df$scenario)) {
      sub <- df %>% filter(scenario == sc)
      p <- add_trace(p, x=sub$week, y=sub$NRS, type="scatter", mode="lines",
                     name=sc, line=list(width=2, color=sc_colors[[sc]]))
    }
    p %>% layout(title="NRS Itch Comparison",
                 xaxis=list(title="Week"),
                 yaxis=list(title="NRS (0-10)"),
                 legend=list(orientation="h", y=-0.25))
  })

  output$plot_compare_tarc <- renderPlotly({
    req(compare_data())
    df <- compare_data() %>% mutate(week = time / 7)
    p <- plot_ly()
    for (sc in unique(df$scenario)) {
      sub <- df %>% filter(scenario == sc)
      p <- add_trace(p, x=sub$week, y=sub$TARC, type="scatter", mode="lines",
                     name=sc, line=list(width=2, color=sc_colors[[sc]]))
    }
    p %>% layout(title="TARC/CCL17 Comparison",
                 xaxis=list(title="Week"),
                 yaxis=list(title="TARC (pg/mL)"),
                 legend=list(orientation="h", y=-0.25))
  })

  output$compare_summary_table <- renderDT({
    req(compare_data())
    df <- compare_data() %>%
      filter(time %in% c(112, 364)) %>%
      mutate(week = time / 7) %>%
      group_by(scenario, week) %>%
      summarise(
        EASI     = round(mean(EASI), 1),
        NRS      = round(mean(NRS), 1),
        `IGA0/1` = ifelse(mean(IGA) <= 1, "Yes", "No"),
        `EASI-75` = ifelse(mean(EASI75) > 0.5, "Yes", "No"),
        TARC     = round(mean(TARC), 0),
        Eos      = round(mean(Eos_blood), 0),
        .groups  = "drop"
      )
    datatable(df, options=list(dom="t", pageLength=20), rownames=FALSE,
              caption="Scenario Response Summary at Weeks 16 & 52")
  })

  ## -----------------------------------------------------------------------
  ## TAB 6 — BARRIER & PRURITUS
  ## -----------------------------------------------------------------------
  output$plot_flg <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>% mutate(week = time / 7)
    make_plotly_line(df, "week", "FLG", "Filaggrin Expression (relative)",
                     "FLG (rel, normal=1.0)", "#27ae60", hline=0.85, hline_color="green")
  })

  output$plot_tewl <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>% mutate(week = time / 7)
    make_plotly_line(df, "week", "TEWL", "Trans-Epidermal Water Loss (g/m²/h)",
                     "TEWL", "#e67e22", hline=10, hline_color="green")
  })

  output$plot_itch_infl <- renderPlotly({
    req(sim_result())
    df <- sim_result() %>% mutate(week = time / 7)
    plot_ly(df) %>%
      add_trace(x=~week, y=~NRS, name="NRS Itch", type="scatter", mode="lines",
                line=list(color="#e74c3c", width=2), yaxis="y") %>%
      add_trace(x=~week, y=~SkinInfl, name="Skin Infl. Index", type="scatter", mode="lines",
                line=list(color="#1a5276", width=2, dash="dash"), yaxis="y2") %>%
      layout(title="Itch-Inflammation Feedback",
             xaxis=list(title="Week"),
             yaxis=list(title="NRS Itch (0-10)", range=c(0,11)),
             yaxis2=list(title="Skin Inflammation (AU)", overlaying="y",
                         side="right", range=c(0,150)),
             legend=list(orientation="h", y=-0.25))
  })

  output$barrier_table <- renderDT({
    req(sim_result())
    df <- sim_result() %>%
      filter(time %in% c(0, 28, 56, 112, 224, 364)) %>%
      mutate(Week = time / 7,
             FLG_r  = round(FLG, 3),
             TEWL_r = round(TEWL, 1)) %>%
      select(Week, FLG_r, TEWL_r)
    names(df) <- c("Week", "FLG (rel.)", "TEWL (g/m²/h)")
    datatable(df, options=list(dom="t", pageLength=10), rownames=FALSE)
  })
}

## ---------------------------------------------------------------------------
## Launch
## ---------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
