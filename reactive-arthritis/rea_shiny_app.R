# ============================================================
# Reactive Arthritis (ReA) QSP — Interactive Shiny Dashboard
# ============================================================
# Tabs: 1) Patient Profile  2) Drug PK  3) Disease Dynamics
#       4) Clinical Endpoints  5) Scenario Comparison  6) Biomarker Panel
# ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(patchwork)
library(scales)

# ============================================================
# EMBED MODEL CODE
# ============================================================
rea_code <- '
$PARAM @annotated
k_path_decay  : 0.008  : Pathogen decay rate (hr-1)
k_path_innate : 0.05   : Innate-mediated pathogen killing (hr-1)
k_abx_gi      : 0      : GI antibiotic effect
k_abx_chlamyd : 0      : Anti-chlamydial antibiotic effect
k_chron_decay : 0.001  : Chronic pathogen decay (hr-1)
k_innate_act  : 0.3    : Innate activation rate
k_innate_decay: 0.05   : Innate resolution (hr-1)
k_neut_prod   : 0.2    : Neutrophil production rate
k_neut_decay  : 0.03   : Neutrophil decay (hr-1)
k_macro_act   : 0.15   : Macrophage activation rate
k_macro_res   : 0.03   : Macrophage resolution (hr-1)
k_Th1_diff    : 0.08   : Th1 differentiation (hr-1)
k_Th17_diff   : 0.07   : Th17 differentiation (hr-1)
k_Treg_prod   : 0.03   : Treg production (hr-1)
k_Th1_decay   : 0.02   : Th1 decay (hr-1)
k_Th17_decay  : 0.025  : Th17 decay (hr-1)
k_Treg_decay  : 0.01   : Treg decay (hr-1)
HLAB27_eff    : 1.6    : HLA-B27 Th17 fold-amplification
USE_HLAB27    : 1      : 1=HLA-B27+, 0=HLA-B27-
k_TNF_prod    : 4.0    : TNF-alpha production (pg/mL/hr)
k_TNF_decay   : 0.4    : TNF-alpha clearance (hr-1)
k_IL17_prod   : 2.5    : IL-17A production (pg/mL/hr)
k_IL17_decay  : 0.25   : IL-17A clearance (hr-1)
k_IL6_prod    : 3.5    : IL-6 production (pg/mL/hr)
k_IL6_decay   : 0.35   : IL-6 clearance (hr-1)
k_IL10_prod   : 1.2    : IL-10 production (pg/mL/hr)
k_IL10_decay  : 0.2    : IL-10 clearance (hr-1)
k_IFNg_prod   : 1.5    : IFN-gamma production (pg/mL/hr)
k_IFNg_decay  : 0.3    : IFN-gamma clearance (hr-1)
k_synov_act   : 0.12   : Synovitis activation
k_synov_res   : 0.015  : Synovitis resolution (hr-1)
k_cart_dmg    : 0.00003: Cartilage damage rate
k_pain_act    : 0.25   : Pain activation
k_pain_res    : 0.08   : Pain resolution (hr-1)
k_CRP_prod    : 1.8    : CRP production by IL-6
k_CRP_decay   : 0.036  : CRP clearance (hr-1)
k_jcount_rate : 0.004  : Joint count equilibration (hr-1)
ka_NSAID  : 1.2   : NSAID absorption (hr-1)
ke_NSAID  : 0.33  : NSAID elimination (hr-1)
F_NSAID   : 0.95  : NSAID bioavailability
IC50_NSAID: 1.8   : NSAID COX-2 IC50 (µg/mL)
Emax_NSAID: 0.85  : NSAID max COX-2 inhibition
ka_SSZ    : 0.25  : SSZ absorption (hr-1)
ke_SSZ    : 0.033 : SSZ elimination (hr-1)
F_SSZ     : 0.10  : SSZ bioavailability
IC50_SSZ  : 10.0  : SSZ IC50 (µg/mL)
Emax_SSZ  : 0.55  : SSZ max effect
ka_MTX    : 0.7   : MTX absorption (hr-1)
ke_MTX    : 0.065 : MTX elimination (hr-1)
F_MTX     : 0.78  : MTX bioavailability
IC50_MTX  : 0.45  : MTX IC50 (µmol/L)
Emax_MTX  : 0.72  : MTX max effect
ka_TNFi   : 0.009 : TNFi SC absorption (hr-1)
ke_TNFi   : 0.0055: TNFi elimination (hr-1)
F_TNFi    : 0.76  : TNFi bioavailability
IC50_TNFi : 1.2   : TNFi IC50 (µg/mL)
Emax_TNFi : 0.95  : TNFi max neutralization
ka_IL17i  : 0.006 : IL17i SC absorption (hr-1)
ke_IL17i  : 0.0025: IL17i elimination (hr-1)
F_IL17i   : 0.73  : IL17i bioavailability
IC50_IL17i: 0.7   : IL17i IC50 (µg/mL)
Emax_IL17i: 0.93  : IL17i max neutralization

$CMT PATH CHRON_PATH INNATE NEUT MACRO TH1 TH17 TREG
     TNF_c IL17_c IL6_c IL10_c IFNg_c SYNOV CARTDMG PAIN CRP_c JCOUNT
     NSAID_DEPOT NSAID_C SSZ_DEPOT SSZ_C MTX_DEPOT MTX_C
     TNFi_DEPOT TNFi_C IL17i_DEPOT IL17i_C

$INIT
PATH=1, CHRON_PATH=0, INNATE=0, NEUT=0.5, MACRO=0.5,
TH1=0.3, TH17=0.2, TREG=0.2,
TNF_c=3, IL17_c=1.5, IL6_c=2, IL10_c=2, IFNg_c=1,
SYNOV=0.2, CARTDMG=0, PAIN=5, CRP_c=1, JCOUNT=0.5,
NSAID_DEPOT=0, NSAID_C=0, SSZ_DEPOT=0, SSZ_C=0,
MTX_DEPOT=0, MTX_C=0, TNFi_DEPOT=0, TNFi_C=0,
IL17i_DEPOT=0, IL17i_C=0

$ODE
double E_NSAID  = Emax_NSAID  * NSAID_C  / (IC50_NSAID  + NSAID_C  + 1e-12);
double E_SSZ    = Emax_SSZ    * SSZ_C    / (IC50_SSZ    + SSZ_C    + 1e-12);
double E_MTX    = Emax_MTX    * MTX_C    / (IC50_MTX    + MTX_C    + 1e-12);
double E_TNFi   = Emax_TNFi   * TNFi_C   / (IC50_TNFi   + TNFi_C   + 1e-12);
double E_IL17i  = Emax_IL17i  * IL17i_C  / (IC50_IL17i  + IL17i_C  + 1e-12);
double HB27 = 1.0 + (HLAB27_eff - 1.0) * USE_HLAB27;
double IL10_supp = 1.0 / (1.0 + IL10_c / 6.0);
dxdt_PATH = -(k_path_decay + k_path_innate * INNATE + k_abx_gi) * PATH;
dxdt_CHRON_PATH = 0.002 * PATH - k_chron_decay * CHRON_PATH - k_abx_chlamyd * CHRON_PATH;
dxdt_INNATE = k_innate_act * (PATH + CHRON_PATH * 0.5) - k_innate_decay * INNATE;
dxdt_NEUT   = k_neut_prod * (INNATE + IL17_c / 20.0) - k_neut_decay * NEUT;
dxdt_MACRO  = k_macro_act * (INNATE + TH1 * 0.4 + IFNg_c / 20.0) - k_macro_res * MACRO;
dxdt_TH1    = k_Th1_diff  * MACRO * IL10_supp * (1.0 - E_MTX) - k_Th1_decay  * TH1;
dxdt_TH17   = k_Th17_diff * MACRO * HB27 * IL10_supp * (1.0 - E_MTX * 0.8) - k_Th17_decay * TH17;
dxdt_TREG   = k_Treg_prod * (1.0 + TH1 * 0.2) - k_Treg_decay * TREG;
dxdt_TNF_c  = k_TNF_prod  * (MACRO + TH1 * 0.7) * (1.0 - E_TNFi) * (1.0 - E_SSZ * 0.5) * (1.0 - E_MTX * 0.6) - k_TNF_decay * TNF_c;
dxdt_IL17_c = k_IL17_prod * TH17 * (1.0 - E_IL17i) - k_IL17_decay * IL17_c;
dxdt_IL6_c  = k_IL6_prod  * MACRO * (1.0 - E_SSZ * 0.45) * (1.0 - E_MTX * 0.4) - k_IL6_decay * IL6_c;
dxdt_IL10_c = k_IL10_prod * TREG - k_IL10_decay * IL10_c;
dxdt_IFNg_c = k_IFNg_prod * TH1  - k_IFNg_decay * IFNg_c;
double synov_drive = k_synov_act * (TNF_c / 5.0 + IL17_c / 3.0 + IL6_c / 5.0 + CHRON_PATH * 2.0);
dxdt_SYNOV = synov_drive - k_synov_res * SYNOV;
dxdt_CARTDMG = k_cart_dmg * SYNOV * (1.0 - CARTDMG);
double pain_drive = k_pain_act * (SYNOV * 8.0 + TNF_c / 4.0 + IL6_c / 6.0);
dxdt_PAIN = pain_drive - k_pain_res * PAIN - E_NSAID * k_pain_res * PAIN;
dxdt_CRP_c = k_CRP_prod * IL6_c - k_CRP_decay * CRP_c;
dxdt_JCOUNT = k_jcount_rate * (12.0 * SYNOV / (SYNOV + 1.5) - JCOUNT);
dxdt_NSAID_DEPOT = -ka_NSAID * NSAID_DEPOT;
dxdt_NSAID_C     =  ka_NSAID * F_NSAID * NSAID_DEPOT - ke_NSAID * NSAID_C;
dxdt_SSZ_DEPOT   = -ka_SSZ * SSZ_DEPOT;
dxdt_SSZ_C       =  ka_SSZ * F_SSZ * SSZ_DEPOT - ke_SSZ * SSZ_C;
dxdt_MTX_DEPOT   = -ka_MTX * MTX_DEPOT;
dxdt_MTX_C       =  ka_MTX * F_MTX * MTX_DEPOT - ke_MTX * MTX_C;
dxdt_TNFi_DEPOT  = -ka_TNFi * TNFi_DEPOT;
dxdt_TNFi_C      =  ka_TNFi * F_TNFi * TNFi_DEPOT - ke_TNFi * TNFi_C;
dxdt_IL17i_DEPOT = -ka_IL17i * IL17i_DEPOT;
dxdt_IL17i_C     =  ka_IL17i * F_IL17i * IL17i_DEPOT - ke_IL17i * IL17i_C;

$TABLE
capture VAS_pain=PAIN; capture swollen_jt=JCOUNT; capture CRP_obs=CRP_c;
capture TNF_obs=TNF_c; capture IL17_obs=IL17_c; capture IL6_obs=IL6_c;
capture IL10_obs=IL10_c; capture IFNg_obs=IFNg_c;
capture SYNOV_obs=SYNOV; capture CARTDMG_pct=CARTDMG*100;
capture PATH_obs=PATH; capture CHRON_obs=CHRON_PATH;
capture TH1_obs=TH1; capture TH17_obs=TH17; capture TREG_obs=TREG;
capture NSAID_obs=NSAID_C; capture SSZ_obs=SSZ_C; capture MTX_obs=MTX_C;
capture TNFi_obs=TNFi_C; capture IL17i_obs=IL17i_C;
'

mod <- mcode("rea_qsp", rea_code, quiet = TRUE)

# ============================================================
# HELPER: Build dosing events
# ============================================================
build_events <- function(tx_type, nsaid_dur, nsaid_dose,
                         ssz_on, ssz_dose, mtx_on, mtx_dose,
                         tnfi_on, tnfi_dose, il17i_on, il17i_dose,
                         abx_on, sim_days) {
  sim_hr <- sim_days * 24
  evs <- list()

  if (nsaid_dur > 0) {
    times_nsaid <- seq(0, min(nsaid_dur * 24 - 12, sim_hr), by = 12)
    evs[["nsaid"]] <- ev(cmt = "NSAID_DEPOT", time = times_nsaid, amt = nsaid_dose)
  }
  if (ssz_on) {
    times_ssz <- seq(0, sim_hr - 12, by = 12)
    evs[["ssz"]] <- ev(cmt = "SSZ_DEPOT", time = times_ssz, amt = ssz_dose)
  }
  if (mtx_on) {
    times_mtx <- seq(168, sim_hr - 168, by = 168)
    evs[["mtx"]] <- ev(cmt = "MTX_DEPOT", time = times_mtx, amt = mtx_dose)
  }
  if (tnfi_on) {
    times_tnfi <- seq(2016, sim_hr - 168, by = 168)
    evs[["tnfi"]] <- ev(cmt = "TNFi_DEPOT", time = times_tnfi, amt = tnfi_dose)
  }
  if (il17i_on) {
    load_t <- seq(2016, 2016 + 4 * 168, by = 168)
    maint_t <- seq(2016 + 4 * 168 + 720, sim_hr - 720, by = 720)
    evs[["il17i"]] <- ev(cmt = "IL17i_DEPOT", time = c(load_t, maint_t), amt = il17i_dose)
  }

  if (length(evs) == 0) return(ev(time = 0, amt = 0, cmt = 1))
  do.call(c, evs)
}

run_simulation <- function(hlab27, abx_gi, abx_chlamyd, sim_days,
                           nsaid_dur, nsaid_dose, ssz_on, ssz_dose,
                           mtx_on, mtx_dose, tnfi_on, tnfi_dose,
                           il17i_on, il17i_dose) {
  params <- list(
    USE_HLAB27 = as.numeric(hlab27),
    k_abx_gi = ifelse(abx_gi, 0.5, 0),
    k_abx_chlamyd = ifelse(abx_chlamyd, 0.3, 0)
  )
  ev_obj <- build_events("custom", nsaid_dur, nsaid_dose, ssz_on, ssz_dose,
                          mtx_on, mtx_dose, tnfi_on, tnfi_dose,
                          il17i_on, il17i_dose, abx_gi, sim_days)
  mod %>%
    param(params) %>%
    mrgsim(ev = ev_obj, end = sim_days * 24, delta = 12, carry_out = "time") %>%
    as.data.frame() %>%
    mutate(day = time / 24)
}

# ============================================================
# UI
# ============================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Reactive Arthritis (ReA) QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "patient",   icon = icon("user")),
      menuItem("Drug PK",             tabName = "pk",        icon = icon("capsules")),
      menuItem("Disease Dynamics",    tabName = "disease",   icon = icon("dna")),
      menuItem("Clinical Endpoints",  tabName = "clinical",  icon = icon("chart-line")),
      menuItem("Scenario Comparison", tabName = "scenario",  icon = icon("balance-scale")),
      menuItem("Biomarker Panel",     tabName = "biomarker", icon = icon("flask"))
    ),
    hr(),
    h5("Global Parameters", style = "padding-left:15px; color:white;"),
    selectInput("hlab27", "HLA-B27 Status:",
                choices = c("Positive (+)" = "1", "Negative (-)" = "0"), selected = "1"),
    selectInput("infect_type", "Triggering Infection:",
                choices = c("Chlamydia trachomatis (GU)" = "chlamydia",
                            "Salmonella (GI)"             = "salmonella",
                            "Campylobacter (GI)"          = "campylobacter",
                            "Yersinia (GI)"               = "yersinia",
                            "Shigella (GI)"               = "shigella"),
                selected = "chlamydia"),
    checkboxInput("abx_gi",      "GI Antibiotic (Fluoroquinolone)", value = FALSE),
    checkboxInput("abx_chlamyd", "Anti-Chlamydial Abx (Doxycycline)", value = FALSE),
    sliderInput("sim_days", "Simulation Duration (days):",
                min = 90, max = 730, value = 365, step = 30)
  ),

  dashboardBody(
    tabItems(

      # ----------------------------------------------------------
      # TAB 1: Patient Profile
      # ----------------------------------------------------------
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "Demographic & Risk Inputs", width = 4, solidHeader = TRUE, status = "primary",
            sliderInput("age",    "Age (years):", 20, 70, 35, 1),
            selectInput("sex",   "Sex:", choices = c("Male" = "M", "Female" = "F"), selected = "M"),
            sliderInput("weight","Weight (kg):", 40, 130, 70, 1),
            selectInput("disease_duration", "Disease Duration at Presentation:",
                        choices = c("Acute (<4 weeks)"     = "acute",
                                    "Subacute (4–12 weeks)" = "subacute",
                                    "Chronic (>12 weeks)"   = "chronic"),
                        selected = "subacute"),
            sliderInput("baseline_crp",  "Baseline CRP (mg/L):", 1, 100, 25, 1),
            sliderInput("baseline_pain", "Baseline VAS Pain (0–100):", 0, 100, 60, 1),
            sliderInput("baseline_jcount","Baseline Swollen Joint Count:", 0, 28, 4, 1)
          ),
          box(title = "Risk Stratification", width = 4, solidHeader = TRUE, status = "warning",
            verbatimTextOutput("risk_score_text"),
            plotOutput("radar_plot", height = "280px")
          ),
          box(title = "ReA Pathophysiology Overview", width = 4, solidHeader = TRUE, status = "info",
            HTML("
              <h5>Reactive Arthritis Key Features</h5>
              <ul>
                <li><b>Triggering infections</b>: Chlamydia (GU), Salmonella/Campylobacter/
                    Yersinia/Shigella (GI)</li>
                <li><b>Triad</b>: Urethritis + Arthritis + Conjunctivitis</li>
                <li><b>HLA-B27</b>: +ve in ~75% ReA (vs. 8% general pop.)</li>
                <li><b>Chronicity</b>: ~15–30% persist >6 months</li>
                <li><b>Axial involvement</b>: Sacroiliitis in 10–40% HLA-B27+</li>
                <li><b>Key cytokines</b>: TNF-α, IL-17A, IL-6</li>
                <li><b>Treatment</b>: NSAIDs → DMARDs → Biologics</li>
              </ul>
            "),
            hr(),
            tableOutput("patient_summary_table")
          )
        ),
        fluidRow(
          box(title = "Disease Course Prediction (Based on Input Profile)", width = 12,
              solidHeader = TRUE, status = "success",
            plotlyOutput("course_prediction_plot", height = "300px")
          )
        )
      ),

      # ----------------------------------------------------------
      # TAB 2: Drug PK
      # ----------------------------------------------------------
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Drug Selection & Dosing", width = 3, solidHeader = TRUE, status = "primary",
            checkboxInput("use_nsaid",  "NSAID (Naproxen)", value = TRUE),
            conditionalPanel("input.use_nsaid",
              sliderInput("nsaid_dose", "Naproxen Dose (mg):", 250, 750, 500, 50),
              sliderInput("nsaid_dur",  "NSAID Duration (days):", 14, 365, 84, 7)
            ),
            hr(),
            checkboxInput("use_ssz", "Sulfasalazine", value = FALSE),
            conditionalPanel("input.use_ssz",
              sliderInput("ssz_dose", "SSZ Dose per dose (mg):", 250, 1500, 1000, 250)
            ),
            hr(),
            checkboxInput("use_mtx", "Methotrexate (weekly SC)", value = FALSE),
            conditionalPanel("input.use_mtx",
              sliderInput("mtx_dose", "MTX Weekly Dose (mg):", 7.5, 25, 15, 2.5)
            ),
            hr(),
            checkboxInput("use_tnfi", "TNF Inhibitor (Etanercept)", value = FALSE),
            conditionalPanel("input.use_tnfi",
              sliderInput("tnfi_dose", "Etanercept Dose (mg/wk):", 25, 50, 50, 25),
              helpText("Started at week 12 (after DMARD trial)")
            ),
            hr(),
            checkboxInput("use_il17i", "IL-17i (Secukinumab)", value = FALSE),
            conditionalPanel("input.use_il17i",
              sliderInput("il17i_dose", "Secukinumab Dose (mg):", 150, 300, 300, 150)
            ),
            actionButton("run_pk", "Run Simulation", class = "btn-primary btn-block")
          ),
          box(title = "Plasma Concentration–Time Profiles", width = 9, solidHeader = TRUE, status = "info",
            plotlyOutput("pk_plot", height = "450px"),
            hr(),
            dataTableOutput("pk_table")
          )
        )
      ),

      # ----------------------------------------------------------
      # TAB 3: Disease Dynamics
      # ----------------------------------------------------------
      tabItem(tabName = "disease",
        fluidRow(
          box(title = "Immune Cell Dynamics", width = 6, solidHeader = TRUE, status = "warning",
            plotlyOutput("immune_cell_plot", height = "380px")
          ),
          box(title = "Cytokine Time Course", width = 6, solidHeader = TRUE, status = "danger",
            plotlyOutput("cytokine_plot", height = "380px")
          )
        ),
        fluidRow(
          box(title = "Pathogen Clearance & Synovial Inflammation", width = 6,
              solidHeader = TRUE, status = "primary",
            plotlyOutput("pathogen_synov_plot", height = "300px")
          ),
          box(title = "T-cell Subset Ratios & Regulatory Balance", width = 6,
              solidHeader = TRUE, status = "success",
            plotlyOutput("tcell_ratio_plot", height = "300px"),
            helpText("Th17/Treg ratio: elevated ratio indicates pro-inflammatory state")
          )
        )
      ),

      # ----------------------------------------------------------
      # TAB 4: Clinical Endpoints
      # ----------------------------------------------------------
      tabItem(tabName = "clinical",
        fluidRow(
          box(title = "VAS Pain Score (0–100 mm)", width = 4, solidHeader = TRUE, status = "danger",
            plotlyOutput("pain_plot", height = "300px"),
            verbatimTextOutput("pain_milestones")
          ),
          box(title = "CRP (mg/L)", width = 4, solidHeader = TRUE, status = "warning",
            plotlyOutput("crp_plot", height = "300px"),
            verbatimTextOutput("crp_milestones")
          ),
          box(title = "Swollen Joint Count", width = 4, solidHeader = TRUE, status = "info",
            plotlyOutput("jcount_plot", height = "300px"),
            verbatimTextOutput("jcount_milestones")
          )
        ),
        fluidRow(
          box(title = "Cartilage Damage (cumulative %)", width = 6,
              solidHeader = TRUE, status = "primary",
            plotlyOutput("cart_damage_plot", height = "250px"),
            helpText("Irreversible damage accumulates with sustained synovitis")
          ),
          box(title = "Clinical Outcome Summary", width = 6, solidHeader = TRUE, status = "success",
            dataTableOutput("clinical_outcome_table")
          )
        )
      ),

      # ----------------------------------------------------------
      # TAB 5: Scenario Comparison
      # ----------------------------------------------------------
      tabItem(tabName = "scenario",
        fluidRow(
          box(title = "Select Scenarios to Compare", width = 3, solidHeader = TRUE, status = "primary",
            checkboxGroupInput("compare_scenarios", "Treatment Arms:",
              choices = c(
                "No Treatment"             = "s1",
                "NSAIDs Only"              = "s2",
                "Antibiotics + NSAIDs"     = "s3",
                "NSAIDs + Sulfasalazine"   = "s4",
                "NSAIDs + SSZ → Etanercept"= "s5"
              ),
              selected = c("s1", "s2", "s4", "s5")
            ),
            selectInput("compare_endpoint", "Primary Endpoint:",
                        choices = c("VAS Pain"     = "VAS_pain",
                                    "CRP (mg/L)"   = "CRP_obs",
                                    "Joint Count"  = "swollen_jt",
                                    "IL-17A"       = "IL17_obs",
                                    "TNF-α"        = "TNF_obs",
                                    "Synovitis"    = "SYNOV_obs"),
                        selected = "CRP_obs"),
            actionButton("run_compare", "Compare Scenarios", class = "btn-success btn-block"),
            hr(),
            h5("Outcome at Selected Timepoint"),
            sliderInput("compare_week", "Week:", 4, 52, 24, 4)
          ),
          box(title = "Treatment Comparison — Time Course", width = 9,
              solidHeader = TRUE, status = "info",
            plotlyOutput("compare_plot", height = "350px"),
            hr(),
            plotlyOutput("compare_forest", height = "220px")
          )
        ),
        fluidRow(
          box(title = "Outcome Table at Key Timepoints", width = 12,
              solidHeader = TRUE, status = "success",
            dataTableOutput("compare_table")
          )
        )
      ),

      # ----------------------------------------------------------
      # TAB 6: Biomarker Panel
      # ----------------------------------------------------------
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "HLA-B27 Effect Modifier Analysis", width = 6,
              solidHeader = TRUE, status = "warning",
            sliderInput("hlab27_fold", "HLA-B27 Th17 Amplification Factor:",
                        1.0, 3.0, 1.6, 0.1),
            plotlyOutput("hlab27_effect_plot", height = "320px"),
            helpText("HLA-B27 amplifies Th17 differentiation via UPR/ER-stress → IL-23")
          ),
          box(title = "Cytokine Correlation Map", width = 6,
              solidHeader = TRUE, status = "danger",
            plotlyOutput("cytokine_heatmap", height = "350px"),
            helpText("Pearson correlation of cytokine dynamics across simulation period")
          )
        ),
        fluidRow(
          box(title = "Drug Effect Sensitivity Analysis (IC50 Perturbation)", width = 6,
              solidHeader = TRUE, status = "primary",
            selectInput("sa_drug", "Drug for Sensitivity Analysis:",
                        choices = c("NSAID"         = "NSAID",
                                    "Sulfasalazine" = "SSZ",
                                    "Etanercept"    = "TNFi"),
                        selected = "TNFi"),
            sliderInput("sa_ic50_mult", "IC50 Multiplier (fold):", 0.1, 5.0, 1.0, 0.1),
            plotlyOutput("sa_plot", height = "250px")
          ),
          box(title = "Chronicity Risk & Biomarker Profile", width = 6,
              solidHeader = TRUE, status = "success",
            plotlyOutput("chronicity_plot", height = "280px"),
            verbatimTextOutput("chronicity_text"),
            hr(),
            helpText("Risk factors: HLA-B27+, Chlamydia persistence, elevated IL-17, delayed treatment")
          )
        )
      )

    )  # end tabItems
  )  # end dashboardBody
)  # end dashboardPage

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {

  # ---- Reactive: run single-arm simulation ----
  sim_data <- eventReactive(input$run_pk, {
    req(input$sim_days)
    run_simulation(
      hlab27       = input$hlab27,
      abx_gi       = input$abx_gi,
      abx_chlamyd  = input$abx_chlamyd,
      sim_days     = input$sim_days,
      nsaid_dur    = if (input$use_nsaid) input$nsaid_dur else 0,
      nsaid_dose   = input$nsaid_dose,
      ssz_on       = input$use_ssz,
      ssz_dose     = input$ssz_dose,
      mtx_on       = input$use_mtx,
      mtx_dose     = input$mtx_dose,
      tnfi_on      = input$use_tnfi,
      tnfi_dose    = input$tnfi_dose,
      il17i_on     = input$use_il17i,
      il17i_dose   = input$il17i_dose
    )
  }, ignoreNULL = FALSE)

  # Also run on app start with defaults
  default_sim <- reactive({
    run_simulation(
      hlab27 = input$hlab27, abx_gi = input$abx_gi,
      abx_chlamyd = input$abx_chlamyd, sim_days = input$sim_days,
      nsaid_dur = 84, nsaid_dose = 500,
      ssz_on = FALSE, ssz_dose = 1000,
      mtx_on = FALSE, mtx_dose = 15,
      tnfi_on = FALSE, tnfi_dose = 50,
      il17i_on = FALSE, il17i_dose = 300
    )
  })

  active_sim <- reactive({
    if (is.null(sim_data())) default_sim() else sim_data()
  })

  # ---- TAB 1: Patient Profile ----
  output$risk_score_text <- renderText({
    hb27  <- input$hlab27 == "1"
    chron <- input$disease_duration == "chronic"
    crp_h <- input$baseline_crp > 30
    pain_h <- input$baseline_pain > 60
    score <- sum(c(hb27, chron, crp_h, pain_h, input$infect_type == "chlamydia"))
    risk_cat <- if (score <= 1) "LOW (Self-limiting likely)"
                else if (score <= 3) "MODERATE (DMARD needed)"
                else "HIGH (Biologic candidacy)"
    paste0("Chronicity Risk Score: ", score, "/5\nClassification: ", risk_cat,
           "\n\nRisk Factors Present:\n",
           if (hb27) "  [+] HLA-B27 positive\n" else "  [-] HLA-B27 negative\n",
           if (input$infect_type == "chlamydia") "  [+] Chlamydia (synovial persistence risk)\n"
           else "  [-] Non-Chlamydia (lower persistence)\n",
           if (crp_h) "  [+] High baseline CRP (>30 mg/L)\n" else "  [-] Moderate CRP\n",
           if (chron) "  [+] Chronic presentation (>12 wk)\n" else "  [-] Acute/subacute onset\n",
           if (pain_h) "  [+] Severe baseline pain (>60 VAS)\n" else "  [-] Moderate pain\n")
  })

  output$patient_summary_table <- renderTable({
    data.frame(
      Parameter = c("Age", "Sex", "Weight", "HLA-B27", "Infection", "CRP", "VAS Pain", "Swollen Joints"),
      Value = c(
        paste(input$age, "years"),
        ifelse(input$sex == "M", "Male", "Female"),
        paste(input$weight, "kg"),
        ifelse(input$hlab27 == "1", "Positive", "Negative"),
        input$infect_type,
        paste(input$baseline_crp, "mg/L"),
        paste(input$baseline_pain, "/100 mm"),
        input$baseline_jcount
      )
    )
  })

  output$course_prediction_plot <- renderPlotly({
    d <- default_sim()
    plot_ly(d, x = ~day) %>%
      add_lines(y = ~VAS_pain, name = "VAS Pain", line = list(color = "red")) %>%
      add_lines(y = ~CRP_obs, name = "CRP (mg/L)", line = list(color = "orange")) %>%
      add_lines(y = ~swollen_jt * 5, name = "Joint Count ×5", line = list(color = "blue", dash = "dot")) %>%
      layout(title = "Predicted Disease Course (NSAID baseline)",
             xaxis = list(title = "Day"),
             yaxis = list(title = "Score / Concentration"),
             shapes = list(
               list(type = "line", x0 = 84, x1 = 84, y0 = 0, y1 = 100,
                    line = list(color = "gray", dash = "dash")),
               list(type = "line", x0 = 168, x1 = 168, y0 = 0, y1 = 100,
                    line = list(color = "gray", dash = "dash"))
             ),
             legend = list(orientation = "h"))
  })

  # ---- TAB 2: Drug PK ----
  output$pk_plot <- renderPlotly({
    d <- active_sim()
    p <- plot_ly(d, x = ~day)
    if (input$use_nsaid)
      p <- p %>% add_lines(y = ~NSAID_obs, name = "Naproxen (µg/mL)", line = list(color = "#FF8C00"))
    if (input$use_ssz)
      p <- p %>% add_lines(y = ~SSZ_obs, name = "Sulfasalazine (µg/mL)", line = list(color = "#4169E1"))
    if (input$use_mtx)
      p <- p %>% add_lines(y = ~MTX_obs * 10, name = "MTX ×10 (µmol/L)", line = list(color = "#228B22"))
    if (input$use_tnfi)
      p <- p %>% add_lines(y = ~TNFi_obs, name = "Etanercept (µg/mL)", line = list(color = "#9400D3"))
    if (input$use_il17i)
      p <- p %>% add_lines(y = ~IL17i_obs, name = "Secukinumab (µg/mL)", line = list(color = "#00CED1"))

    p %>% layout(title = "Drug Plasma Concentration–Time Profiles",
                 xaxis = list(title = "Day"),
                 yaxis = list(title = "Concentration"),
                 legend = list(orientation = "h"))
  })

  output$pk_table <- renderDataTable({
    d <- active_sim()
    wks <- c(1, 2, 4, 8, 12, 24, 36, 52)
    tps <- sapply(wks, function(w) which.min(abs(d$day - w * 7)))
    tps <- tps[tps <= nrow(d)]
    res <- d[tps, c("day", "NSAID_obs", "SSZ_obs", "MTX_obs", "TNFi_obs", "IL17i_obs")]
    res$day <- round(res$day, 0)
    names(res) <- c("Day", "Naproxen (µg/mL)", "SSZ (µg/mL)", "MTX (µmol/L)",
                    "Etanercept (µg/mL)", "Secukinumab (µg/mL)")
    datatable(round(res, 3), options = list(pageLength = 8, dom = "t"))
  })

  # ---- TAB 3: Disease Dynamics ----
  output$immune_cell_plot <- renderPlotly({
    d <- active_sim()
    plot_ly(d, x = ~day) %>%
      add_lines(y = ~TH1_obs, name = "Th1 Cells", line = list(color = "red")) %>%
      add_lines(y = ~TH17_obs, name = "Th17 Cells", line = list(color = "darkorange")) %>%
      add_lines(y = ~TREG_obs, name = "Treg Cells", line = list(color = "green")) %>%
      add_lines(y = ~NEUT / 3, name = "Neutrophils /3", line = list(color = "purple", dash = "dot")) %>%
      layout(title = "Immune Cell Dynamics",
             xaxis = list(title = "Day"), yaxis = list(title = "Normalized Units"),
             legend = list(orientation = "h"))
  })

  output$cytokine_plot <- renderPlotly({
    d <- active_sim()
    plot_ly(d, x = ~day) %>%
      add_lines(y = ~TNF_obs,  name = "TNF-α (pg/mL)",  line = list(color = "#CC0000")) %>%
      add_lines(y = ~IL17_obs, name = "IL-17A (pg/mL)", line = list(color = "#FF6600")) %>%
      add_lines(y = ~IL6_obs,  name = "IL-6 (pg/mL)",   line = list(color = "#0066CC")) %>%
      add_lines(y = ~IL10_obs, name = "IL-10 (pg/mL)",  line = list(color = "#009933")) %>%
      add_lines(y = ~IFNg_obs, name = "IFN-γ (pg/mL)",  line = list(color = "#9900CC", dash = "dash")) %>%
      layout(title = "Cytokine Time Course",
             xaxis = list(title = "Day"), yaxis = list(title = "pg/mL"),
             legend = list(orientation = "h"))
  })

  output$pathogen_synov_plot <- renderPlotly({
    d <- active_sim()
    plot_ly(d, x = ~day) %>%
      add_lines(y = ~PATH_obs, name = "Acute Pathogen Load", line = list(color = "red")) %>%
      add_lines(y = ~CHRON_obs * 10, name = "Chronic Synovial Pathogen ×10",
                line = list(color = "darkred", dash = "dot")) %>%
      add_lines(y = ~SYNOV_obs, name = "Synovitis Score", line = list(color = "blue")) %>%
      layout(title = "Pathogen Clearance & Synovitis",
             xaxis = list(title = "Day"), yaxis = list(title = "Normalized Score"),
             legend = list(orientation = "h"))
  })

  output$tcell_ratio_plot <- renderPlotly({
    d <- active_sim()
    d$th17_treg <- d$TH17_obs / (d$TREG_obs + 0.01)
    d$th1_treg  <- d$TH1_obs  / (d$TREG_obs + 0.01)
    plot_ly(d, x = ~day) %>%
      add_lines(y = ~th17_treg, name = "Th17/Treg ratio", line = list(color = "orange")) %>%
      add_lines(y = ~th1_treg,  name = "Th1/Treg ratio",  line = list(color = "red")) %>%
      add_lines(y = ~TREG_obs, name = "Treg (abs)", line = list(color = "green", dash = "dot")) %>%
      layout(title = "T-cell Ratios (pro-inflammatory balance)",
             xaxis = list(title = "Day"), yaxis = list(title = "Ratio / AU"),
             legend = list(orientation = "h"))
  })

  # ---- TAB 4: Clinical Endpoints ----
  output$pain_plot <- renderPlotly({
    d <- active_sim()
    plot_ly(d, x = ~day, y = ~VAS_pain, type = "scatter", mode = "lines",
            line = list(color = "red", width = 2)) %>%
      add_lines(y = rep(30, nrow(d)), name = "Mild threshold",
                line = list(color = "orange", dash = "dash")) %>%
      layout(title = "VAS Pain (0–100 mm)",
             xaxis = list(title = "Day"),
             yaxis = list(title = "mm", range = c(0, 100)),
             showlegend = FALSE)
  })

  output$pain_milestones <- renderText({
    d <- active_sim()
    t50 <- d$day[which(d$VAS_pain <= 50)[1]]
    t30 <- d$day[which(d$VAS_pain <= 30)[1]]
    paste0("Pain <50 (moderate): day ", ifelse(is.na(t50), ">365", round(t50)),
           "\nPain <30 (mild): day ", ifelse(is.na(t30), ">365", round(t30)))
  })

  output$crp_plot <- renderPlotly({
    d <- active_sim()
    plot_ly(d, x = ~day, y = ~CRP_obs, type = "scatter", mode = "lines",
            line = list(color = "orange", width = 2)) %>%
      add_lines(y = rep(5, nrow(d)), name = "ULN 5 mg/L",
                line = list(color = "red", dash = "dash")) %>%
      layout(title = "CRP (mg/L)", xaxis = list(title = "Day"),
             yaxis = list(title = "mg/L"), showlegend = FALSE)
  })

  output$crp_milestones <- renderText({
    d <- active_sim()
    t_norm <- d$day[which(d$CRP_obs <= 5)[1]]
    final_crp <- round(tail(d$CRP_obs, 1), 1)
    paste0("CRP normalizes (<5 mg/L): day ", ifelse(is.na(t_norm), ">365", round(t_norm)),
           "\nFinal CRP: ", final_crp, " mg/L")
  })

  output$jcount_plot <- renderPlotly({
    d <- active_sim()
    plot_ly(d, x = ~day, y = ~swollen_jt, type = "scatter", mode = "lines",
            line = list(color = "steelblue", width = 2)) %>%
      layout(title = "Swollen Joint Count", xaxis = list(title = "Day"),
             yaxis = list(title = "Joints", range = c(0, 12)))
  })

  output$jcount_milestones <- renderText({
    d <- active_sim()
    t0j <- d$day[which(d$swollen_jt <= 1)[1]]
    paste0("Joint count ≤1: day ", ifelse(is.na(t0j), ">365", round(t0j)),
           "\nFinal joint count: ", round(tail(d$swollen_jt, 1), 1))
  })

  output$cart_damage_plot <- renderPlotly({
    d <- active_sim()
    plot_ly(d, x = ~day, y = ~CARTDMG_pct, type = "scatter", mode = "lines",
            fill = "tozeroy", fillcolor = "rgba(255,50,50,0.2)",
            line = list(color = "darkred")) %>%
      layout(title = "Cumulative Cartilage Damage (%)",
             xaxis = list(title = "Day"), yaxis = list(title = "%", range = c(0, 15)))
  })

  output$clinical_outcome_table <- renderDataTable({
    d <- active_sim()
    wks <- c(4, 8, 12, 24, 36, 52)
    rows <- lapply(wks, function(w) {
      tp <- which.min(abs(d$day - w * 7))
      data.frame(
        Week = w,
        `VAS Pain` = round(d$VAS_pain[tp], 1),
        `CRP (mg/L)` = round(d$CRP_obs[tp], 1),
        `Joint Count` = round(d$swollen_jt[tp], 1),
        `IL-17A (pg/mL)` = round(d$IL17_obs[tp], 1),
        `TNF-α (pg/mL)` = round(d$TNF_obs[tp], 1),
        `Cart Dmg (%)` = round(d$CARTDMG_pct[tp], 2),
        check.names = FALSE
      )
    })
    datatable(do.call(rbind, rows), options = list(pageLength = 7, dom = "t"))
  })

  # ---- TAB 5: Scenario Comparison ----
  compare_data <- eventReactive(input$run_compare, {
    scenarios <- list()
    params_base <- list(USE_HLAB27 = as.numeric(input$hlab27),
                        k_abx_gi = 0, k_abx_chlamyd = 0)
    sim_hr <- input$sim_days * 24

    if ("s1" %in% input$compare_scenarios) {
      scenarios[["S1: No Treatment"]] <- mod %>%
        param(params_base) %>%
        mrgsim(end = sim_hr, delta = 12) %>%
        as.data.frame() %>%
        mutate(day = time / 24, scenario = "S1: No Treatment")
    }
    if ("s2" %in% input$compare_scenarios) {
      nsaid_ev <- ev(cmt = "NSAID_DEPOT", time = seq(0, 2016, by = 12), amt = 500)
      scenarios[["S2: NSAIDs Only"]] <- mod %>%
        param(params_base) %>%
        mrgsim(ev = nsaid_ev, end = sim_hr, delta = 12) %>%
        as.data.frame() %>%
        mutate(day = time / 24, scenario = "S2: NSAIDs Only")
    }
    if ("s3" %in% input$compare_scenarios) {
      nsaid_ev <- ev(cmt = "NSAID_DEPOT", time = seq(0, 2016, by = 12), amt = 500)
      scenarios[["S3: ABX + NSAIDs"]] <- mod %>%
        param(c(params_base, list(k_abx_gi = 0.5, k_abx_chlamyd = 0.3))) %>%
        mrgsim(ev = nsaid_ev, end = sim_hr, delta = 12) %>%
        as.data.frame() %>%
        mutate(day = time / 24, scenario = "S3: ABX + NSAIDs")
    }
    if ("s4" %in% input$compare_scenarios) {
      nsaid_ev <- ev(cmt = "NSAID_DEPOT", time = seq(0, 2016, by = 12), amt = 500)
      ssz_ev   <- ev(cmt = "SSZ_DEPOT",   time = seq(0, sim_hr - 12, by = 12), amt = 1000)
      scenarios[["S4: NSAIDs + SSZ"]] <- mod %>%
        param(params_base) %>%
        mrgsim(ev = c(nsaid_ev, ssz_ev), end = sim_hr, delta = 12) %>%
        as.data.frame() %>%
        mutate(day = time / 24, scenario = "S4: NSAIDs + SSZ")
    }
    if ("s5" %in% input$compare_scenarios) {
      nsaid_ev <- ev(cmt = "NSAID_DEPOT", time = seq(0, 2016, by = 12), amt = 500)
      ssz_ev   <- ev(cmt = "SSZ_DEPOT",   time = seq(0, sim_hr - 12, by = 12), amt = 1000)
      etn_ev   <- ev(cmt = "TNFi_DEPOT",  time = seq(2016, sim_hr - 168, by = 168), amt = 50)
      scenarios[["S5: NSAIDs + SSZ → ETN"]] <- mod %>%
        param(params_base) %>%
        mrgsim(ev = c(nsaid_ev, ssz_ev, etn_ev), end = sim_hr, delta = 12) %>%
        as.data.frame() %>%
        mutate(day = time / 24, scenario = "S5: NSAIDs + SSZ → ETN")
    }
    bind_rows(scenarios)
  }, ignoreNULL = FALSE)

  output$compare_plot <- renderPlotly({
    d <- compare_data()
    req(nrow(d) > 0)
    endpoint <- input$compare_endpoint
    pal <- c("#8B0000","#FF8C00","#228B22","#4169E1","#9400D3")
    scens <- unique(d$scenario)
    p <- plot_ly()
    for (i in seq_along(scens)) {
      sub <- filter(d, scenario == scens[i])
      p <- p %>% add_lines(data = sub, x = ~day, y = as.formula(paste0("~", endpoint)),
                           name = scens[i], line = list(color = pal[i], width = 1.8))
    }
    p %>% layout(title = paste("Comparison:", endpoint),
                 xaxis = list(title = "Day"),
                 yaxis = list(title = endpoint),
                 legend = list(orientation = "h"))
  })

  output$compare_forest <- renderPlotly({
    d <- compare_data()
    req(nrow(d) > 0)
    endpoint <- input$compare_endpoint
    target_day <- input$compare_week * 7
    res <- d %>%
      group_by(scenario) %>%
      summarise(val = get(endpoint)[which.min(abs(day - target_day))], .groups = "drop")
    plot_ly(res, x = ~val, y = ~scenario, type = "bar", orientation = "h",
            marker = list(color = c("#8B0000","#FF8C00","#228B22","#4169E1","#9400D3"))) %>%
      layout(title = paste0(endpoint, " at Week ", input$compare_week),
             xaxis = list(title = endpoint), yaxis = list(title = ""))
  })

  output$compare_table <- renderDataTable({
    d <- compare_data()
    req(nrow(d) > 0)
    wks <- c(4, 12, 24, 52)
    rows <- lapply(unique(d$scenario), function(sc) {
      sub <- filter(d, scenario == sc)
      lapply(wks, function(w) {
        tp <- which.min(abs(sub$day - w * 7))
        data.frame(Scenario = sc, Week = w,
                   `VAS Pain` = round(sub$VAS_pain[tp], 1),
                   `CRP (mg/L)` = round(sub$CRP_obs[tp], 1),
                   `Swollen Jts` = round(sub$swollen_jt[tp], 1),
                   `IL-17A` = round(sub$IL17_obs[tp], 1),
                   `Cart Dmg (%)` = round(sub$CARTDMG_pct[tp], 2),
                   check.names = FALSE)
      }) %>% bind_rows()
    }) %>% bind_rows()
    datatable(rows, options = list(pageLength = 20, scrollX = TRUE))
  })

  # ---- TAB 6: Biomarker Panel ----
  output$hlab27_effect_plot <- renderPlotly({
    fold <- input$hlab27_fold
    sim_pos <- mod %>%
      param(list(USE_HLAB27 = 1, HLAB27_eff = fold)) %>%
      mrgsim(end = 365 * 24, delta = 12) %>% as.data.frame() %>%
      mutate(day = time / 24, grp = paste0("HLA-B27+ (×", fold, ")"))
    sim_neg <- mod %>%
      param(list(USE_HLAB27 = 0)) %>%
      mrgsim(end = 365 * 24, delta = 12) %>% as.data.frame() %>%
      mutate(day = time / 24, grp = "HLA-B27–")

    d <- bind_rows(sim_pos, sim_neg)
    plot_ly(d, x = ~day, y = ~IL17_obs, color = ~grp,
            type = "scatter", mode = "lines") %>%
      layout(title = "HLA-B27 Effect on IL-17A Dynamics",
             xaxis = list(title = "Day"), yaxis = list(title = "IL-17A (pg/mL)"),
             legend = list(orientation = "h"))
  })

  output$cytokine_heatmap <- renderPlotly({
    d <- active_sim()
    cytokines <- c("TNF_obs", "IL17_obs", "IL6_obs", "IL10_obs", "IFNg_obs")
    cmat <- cor(d[, cytokines], use = "complete.obs")
    rownames(cmat) <- colnames(cmat) <- c("TNF-α", "IL-17A", "IL-6", "IL-10", "IFN-γ")
    plot_ly(x = colnames(cmat), y = rownames(cmat), z = cmat,
            type = "heatmap", colorscale = "RdBu", reversescale = TRUE,
            zmin = -1, zmax = 1) %>%
      layout(title = "Cytokine Correlation Matrix")
  })

  output$sa_plot <- renderPlotly({
    mult <- input$sa_ic50_mult
    drug <- input$sa_drug
    param_name <- paste0("IC50_", drug)

    base_val <- switch(drug, NSAID = 1.8, SSZ = 10.0, TNFi = 1.2)
    new_val <- base_val * mult

    sim_new <- mod %>%
      param(setNames(list(new_val), param_name)) %>%
      mrgsim(ev = ev(cmt = if (drug == "NSAID") "NSAID_DEPOT"
                     else if (drug == "SSZ") "SSZ_DEPOT" else "TNFi_DEPOT",
                     time = seq(0, 8760 - 12, by = if (drug == "TNFi") 168 else 12),
                     amt = switch(drug, NSAID = 500, SSZ = 1000, TNFi = 50)),
             end = 8760, delta = 12) %>%
      as.data.frame() %>%
      mutate(day = time / 24, grp = paste0("IC50 ×", mult))

    sim_base <- mod %>%
      mrgsim(ev = ev(cmt = if (drug == "NSAID") "NSAID_DEPOT"
                     else if (drug == "SSZ") "SSZ_DEPOT" else "TNFi_DEPOT",
                     time = seq(0, 8760 - 12, by = if (drug == "TNFi") 168 else 12),
                     amt = switch(drug, NSAID = 500, SSZ = 1000, TNFi = 50)),
             end = 8760, delta = 12) %>%
      as.data.frame() %>%
      mutate(day = time / 24, grp = "IC50 Base")

    d <- bind_rows(sim_new, sim_base)
    plot_ly(d, x = ~day, y = ~CRP_obs, color = ~grp,
            type = "scatter", mode = "lines") %>%
      layout(title = paste0("CRP Sensitivity: ", drug, " IC50"),
             xaxis = list(title = "Day"), yaxis = list(title = "CRP (mg/L)"),
             legend = list(orientation = "h"))
  })

  output$chronicity_plot <- renderPlotly({
    d <- active_sim()
    d6m <- filter(d, day <= 180)
    plot_ly(d6m, x = ~day) %>%
      add_lines(y = ~SYNOV_obs,       name = "Synovitis",      line = list(color = "blue")) %>%
      add_lines(y = ~CHRON_obs * 20,  name = "Chronic Pathogen ×20", line = list(color = "red", dash = "dot")) %>%
      add_lines(y = ~CARTDMG_pct * 3, name = "Cart Damage ×3", line = list(color = "darkred")) %>%
      layout(title = "Chronicity Indicators (First 6 Months)",
             xaxis = list(title = "Day"),
             yaxis = list(title = "Score (normalized)"),
             legend = list(orientation = "h"))
  })

  output$chronicity_text <- renderText({
    d <- active_sim()
    d6m <- d[which.min(abs(d$day - 180)), ]
    chronic_risk <- (as.numeric(input$hlab27) * 0.3 +
                     (input$infect_type == "chlamydia") * 0.25 +
                     (d6m$SYNOV_obs > 2) * 0.25 +
                     (d6m$CHRON_obs > 0.01) * 0.2)
    paste0("6-Month Chronicity Risk: ", round(chronic_risk * 100, 0), "%\n",
           "Synovitis at 6 months: ", round(d6m$SYNOV_obs, 2), "\n",
           "IL-17A at 6 months: ", round(d6m$IL17_obs, 1), " pg/mL\n",
           "Cart Damage at 6 months: ", round(d6m$CARTDMG_pct, 2), "%")
  })

}  # end server

# ============================================================
# LAUNCH
# ============================================================
shinyApp(ui = ui, server = server)
