## =============================================================================
## Congenital Adrenal Hyperplasia (CAH) – Interactive QSP Shiny Dashboard
## 선천성 부신증식증 (21-수산화효소 결핍증) 정량적 시스템 약리학 대시보드
## =============================================================================
## Tabs:
##   1. Patient Profile     – Genotype, phenotype, body weight
##   2. Drug PK             – Plasma concentration-time profiles
##   3. Steroid Biomarkers  – 17-OHP, ACTH, Androstenedione, Cortisol
##   4. Clinical Endpoints  – Height SDS, Bone Age, BMD, Renin
##   5. Scenario Comparison – Side-by-side comparison of 6 treatment arms
##   6. Biomarker Dashboard – Heatmap, target attainment, radar chart
## =============================================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

# ---- Inline mrgsolve model (compact version for Shiny) ----------------------
cah_model_code <- '
$PARAM @annotated
CYP21A2_res : 0.01  : Residual 21-OH activity (0=null, 0.01=SW, 0.02=SV, 0.20=NC)
k_CRH_prod  : 0.30  : CRH production (1/h)
k_CRH_deg   : 7.00  : CRH degradation (1/h)
CRH_ss      : 1.0   : CRH baseline
k_ACTH_prod : 2.0   : ACTH production (pmol/L/h)
k_ACTH_deg  : 0.90  : ACTH degradation (1/h)
ACTH_ss     : 15.0  : ACTH baseline (pmol/L)
IC50_GC     : 50.0  : GC IC50 ACTH suppression (nmol/L)
n_GC        : 2.0   : GC Hill n
EC50_CRH    : 0.50  : CRH EC50 for ACTH (normalized)
n_CRH       : 2.0   : Hill n CRH->ACTH
CHOL_ss     : 100.0 : Cholesterol pool
k_PREG      : 0.20  : Cholesterol->Pregnenolone
k_PROG      : 0.40  : Preg->Prog
k_17OHP     : 0.45  : Prog->17-OHP
k_DHEA      : 0.25  : 17-OH Preg->DHEA
k_A4        : 0.12  : DHEA->A4 + 17-OHP shunt
k_21OH      : 0.55  : CYP21A2 flux
k_CORT      : 0.40  : Compound S->Cortisol
k_ALD       : 0.30  : DOC->Aldosterone
k_A4_T      : 0.20  : A4->Testosterone
k_17OHP_deg : 0.30  : 17-OHP clearance
k_A4_deg    : 0.50  : A4 clearance
k_T_deg     : 0.25  : T clearance
k_CORT_deg  : 0.45  : Cortisol clearance
k_ALD_deg   : 0.40  : Aldosterone clearance
k_RENIN_prod: 0.50  : Renin production
k_RENIN_deg : 0.30  : Renin clearance
IC50_ALD_REN: 0.50  : Ald IC50 renin suppression
k_BA_adv    : 0.15  : Bone age advance rate
k_bone_GC   : 0.05  : GC->BMD loss
BMD_ss      : 1.0   : Normal BMD
k_BMD_rec   : 0.002 : BMD recovery
HC_F : 0.95 : HC bioavailability
HC_ka: 2.50 : HC absorption (1/h)
HC_CL: 17.0 : HC clearance (L/h)
HC_V1: 15.0 : HC central volume (L)
HC_Q : 5.00 : HC Q (L/h)
HC_V2: 20.0 : HC peripheral volume (L)
PRED_F : 0.82 : PRED bioavailability
PRED_ka: 2.00 : PRED absorption (1/h)
PRED_CL: 10.5 : PRED clearance (L/h)
PRED_V : 35.0 : PRED volume (L)
FC_F : 0.90 : FC bioavailability
FC_ka: 3.00 : FC absorption (1/h)
FC_CL: 8.00 : FC clearance (L/h)
FC_V : 25.0 : FC volume (L)
TILD_F  : 0.65  : Tildacerfont bioavailability
TILD_ka : 0.80  : Tildacerfont absorption (1/h)
TILD_CL : 12.0  : Tildacerfont clearance (L/h)
TILD_V  : 80.0  : Tildacerfont volume (L)
TILD_IC50: 0.004: Tildacerfont IC50 (mg/L)
TILD_n  : 1.20  : Tildacerfont Hill n
CRINE_F  : 0.50  : Crinecerfont bioavailability
CRINE_ka : 1.20  : Crinecerfont absorption (1/h)
CRINE_CL : 18.0  : Crinecerfont clearance (L/h)
CRINE_V  : 95.0  : Crinecerfont volume (L)
CRINE_IC50: 0.0003: Crinecerfont IC50 (mg/L)
CRINE_n  : 1.50  : Crinecerfont Hill n

$CMT @annotated
CRH : CRH
ACTH: ACTH (pmol/L)
OHP17: 17-OHP (nmol/L)
A4: Androstenedione (nmol/L)
TESTO: Testosterone (nmol/L)
CORTISOL: Cortisol (nmol/L)
ALDOST: Aldosterone (pmol/L)
RENIN: Renin (normalized)
HEIGHT_SDS: Height SDS
BONE_AGE: Bone age advance (yr)
BMD: BMD (normalized)
HC_GUT: HC gut (mg)
HC_CENT: HC central (mg/L)
HC_PERI: HC peripheral (mg)
PRED_GUT: PRED gut (mg)
PRED_CENT: PRED central (mg/L)
FC_GUT: FC gut (mg)
FC_CENT: FC central (mg/L)
TILD_GUT: Tild gut (mg)
TILD_CENT: Tild central (mg/L)
CRINE_GUT: Crine gut (mg)
CRINE_CENT: Crine central (mg/L)

$MAIN
double GC_eff = HC_CENT * 1.0 + PRED_CENT * 4.0;
double GC_nmol = GC_eff * 2760.0;
double GC_inh  = pow(GC_nmol,n_GC) / (pow(IC50_GC,n_GC) + pow(GC_nmol,n_GC));
double GC_inh_t = 1.0 - GC_inh;
double TILD_occ = (TILD_CENT > 0) ?
  pow(TILD_CENT,TILD_n) / (pow(TILD_IC50,TILD_n)+pow(TILD_CENT,TILD_n)) : 0.0;
double CRINE_occ = (CRINE_CENT > 0) ?
  pow(CRINE_CENT,CRINE_n)/(pow(CRINE_IC50,CRINE_n)+pow(CRINE_CENT,CRINE_n)) : 0.0;
double CRF1_blk = std::min(1.0, TILD_occ + CRINE_occ);
double CRH_eff = CRH * (1.0 - CRF1_blk);
double ACTH_ratio = ACTH / ACTH_ss;
double F21 = CYP21A2_res;
double A4_ex = std::max(0.0, A4 - 7.0);
double T_ex  = std::max(0.0, TESTO - 1.5);

$ODE
double crh_hill = pow(CRH_eff,n_CRH)/(pow(EC50_CRH,n_CRH)+pow(CRH_eff,n_CRH));
dxdt_CRH  = k_CRH_prod * CRH_ss * (1.0 - GC_inh * 0.5) - k_CRH_deg * CRH;
dxdt_ACTH = k_ACTH_prod * crh_hill * GC_inh_t - k_ACTH_deg * ACTH;

double preg  = k_PREG * ACTH_ratio * CHOL_ss;
double prog  = k_PROG * ACTH_ratio;
double A4_17OHP_shunt = 0.10 * OHP17 * (1.0 - F21);

dxdt_OHP17   = k_17OHP * prog * ACTH_ratio -
               k_21OH * F21 * OHP17 - k_17OHP_deg * OHP17;
dxdt_A4      = k_A4 * ACTH_ratio + A4_17OHP_shunt -
               (k_A4_T + k_A4_deg) * A4;
dxdt_TESTO   = k_A4_T * A4 - k_T_deg * TESTO;
dxdt_CORTISOL= k_CORT * k_21OH * F21 * OHP17 * 0.3 -
               k_CORT_deg * CORTISOL +
               HC_CENT * 1.0 * 2760.0 * 0.005;
dxdt_ALDOST  = k_ALD * k_21OH * F21 * OHP17 * 0.1 +
               40.0 * RENIN - k_ALD_deg * ALDOST;
double ALD_n  = ALDOST / 200.0;
double ren_s  = 1.0 / (1.0 + pow(ALD_n/IC50_ALD_REN,2.0));
dxdt_RENIN   = k_RENIN_prod * ren_s - k_RENIN_deg * RENIN;
dxdt_HEIGHT_SDS = -k_BA_adv * A4_ex - 0.005 * std::max(0.0, GC_nmol - 300.0);
dxdt_BONE_AGE   = k_BA_adv * (A4_ex + T_ex * 2.0);
dxdt_BMD        = -k_bone_GC * std::max(0.0, GC_nmol - 300.0) * 0.0005 +
                  k_BMD_rec * (BMD_ss - BMD);
dxdt_HC_GUT  = -HC_ka * HC_GUT;
dxdt_HC_CENT = HC_F * HC_ka * HC_GUT / HC_V1 -
               (HC_CL/HC_V1 + HC_Q/HC_V1)*HC_CENT + HC_Q/HC_V2*HC_PERI;
dxdt_HC_PERI = HC_Q/HC_V1 * HC_CENT - HC_Q/HC_V2 * HC_PERI;
dxdt_PRED_GUT  = -PRED_ka * PRED_GUT;
dxdt_PRED_CENT = PRED_F * PRED_ka * PRED_GUT / PRED_V - PRED_CL/PRED_V*PRED_CENT;
dxdt_FC_GUT  = -FC_ka * FC_GUT;
dxdt_FC_CENT = FC_F * FC_ka * FC_GUT / FC_V - FC_CL/FC_V * FC_CENT;
dxdt_TILD_GUT  = -TILD_ka * TILD_GUT;
dxdt_TILD_CENT = TILD_F*TILD_ka*TILD_GUT/TILD_V - TILD_CL/TILD_V*TILD_CENT;
dxdt_CRINE_GUT  = -CRINE_ka * CRINE_GUT;
dxdt_CRINE_CENT = CRINE_F*CRINE_ka*CRINE_GUT/CRINE_V - CRINE_CL/CRINE_V*CRINE_CENT;

$TABLE
capture ACTH_pgmL = ACTH * 22.0;
capture A4_nmol   = A4;
capture T_nmol    = TESTO;
capture CRF1_block_pct = CRF1_blk * 100.0;
capture TILD_pct  = TILD_occ * 100.0;
capture CRINE_pct = CRINE_occ * 100.0;
capture GC_nmol_val = GC_nmol;
capture HC_nmol   = HC_CENT * 2760.0;
capture PRED_nmol = PRED_CENT * 2950.0;
capture TILD_conc_mgL = TILD_CENT;
capture CRINE_conc_mgL = CRINE_CENT;

$INIT @annotated
CRH = 1.0 : CRH
ACTH = 50.0: ACTH (pmol/L; elevated in CAH)
OHP17 = 120.0: 17-OHP nmol/L
A4 = 25.0  : A4 nmol/L
TESTO = 4.0: T nmol/L
CORTISOL = 50.0: Cortisol nmol/L
ALDOST = 50.0  : Aldosterone pmol/L
RENIN = 3.0    : Renin normalized
HEIGHT_SDS = 1.5: Height SDS
BONE_AGE = 2.0  : Bone age advance (yr)
BMD = 0.95      : BMD
HC_GUT=0 HC_CENT=0 HC_PERI=0
PRED_GUT=0 PRED_CENT=0
FC_GUT=0 FC_CENT=0
TILD_GUT=0 TILD_CENT=0
CRINE_GUT=0 CRINE_CENT=0

$SET delta=1 end=8760
'

mod <- mcode("cah_shiny", cah_model_code)

# ---- Helper: build dosing events ----
build_events <- function(HC_dose, HC_freq, PRED_dose, PRED_freq,
                         FC_dose, TILD_dose, CRINE_dose, CRINE_freq,
                         days = 365) {
  ev_list <- list()
  # HC
  if (HC_dose > 0) {
    HC_times <- switch(as.character(HC_freq),
                       "1" = 8, "2" = c(8, 20), "3" = c(7, 13, 19))
    for (d in 0:(days-1)) for (t in HC_times)
      ev_list[[length(ev_list)+1]] <- ev(amt = HC_dose/HC_freq,
                                          time = d*24+t, cmt="HC_GUT")
  }
  # PRED
  if (PRED_dose > 0) {
    PRED_times <- if (PRED_freq == 2) c(8,20) else 8
    for (d in 0:(days-1)) for (t in PRED_times)
      ev_list[[length(ev_list)+1]] <- ev(amt = PRED_dose/PRED_freq,
                                          time = d*24+t, cmt="PRED_GUT")
  }
  # FC
  if (FC_dose > 0) {
    for (d in 0:(days-1))
      ev_list[[length(ev_list)+1]] <- ev(amt = FC_dose/1000,
                                          time = d*24+8, cmt="FC_GUT")
  }
  # Tildacerfont
  if (TILD_dose > 0) {
    for (d in 0:(days-1))
      ev_list[[length(ev_list)+1]] <- ev(amt = TILD_dose,
                                          time = d*24+8, cmt="TILD_GUT")
  }
  # Crinecerfont
  if (CRINE_dose > 0) {
    CRINE_times <- if (CRINE_freq == 2) c(8,20) else 8
    for (d in 0:(days-1)) for (t in CRINE_times)
      ev_list[[length(ev_list)+1]] <- ev(amt = CRINE_dose/CRINE_freq,
                                          time = d*24+t, cmt="CRINE_GUT")
  }
  if (length(ev_list) == 0) return(NULL)
  do.call(c, ev_list)
}

run_sim <- function(mut_type, HC_dose, HC_freq, PRED_dose, PRED_freq,
                    FC_dose, TILD_dose, CRINE_dose, CRINE_freq, days = 365) {
  cyp21 <- switch(mut_type,
                  "SW (Null)"           = 0.005,
                  "SV (Simple Virilizing)" = 0.020,
                  "NC (Non-Classical)"  = 0.200)
  evts <- build_events(HC_dose, HC_freq, PRED_dose, PRED_freq,
                       FC_dose, TILD_dose, CRINE_dose, CRINE_freq, days)
  sim <- mod %>%
    param(CYP21A2_res = cyp21) %>%
    mrgsim(events = evts, end = days*24, delta = 4) %>%
    as.data.frame() %>%
    mutate(time_days = time / 24)
  sim
}

# ---- UI -----------------------------------------------------------------------
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "CAH QSP Dashboard | 선천성 부신증식증"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "patient",   icon = icon("user")),
      menuItem("Drug PK",             tabName = "pk",        icon = icon("pills")),
      menuItem("Steroid Biomarkers",  tabName = "steroids",  icon = icon("vial")),
      menuItem("Clinical Endpoints",  tabName = "endpoints", icon = icon("chart-line")),
      menuItem("Scenario Comparison", tabName = "scenarios", icon = icon("layer-group")),
      menuItem("Biomarker Dashboard", tabName = "heatmap",   icon = icon("th"))
    ),

    hr(),
    h5("Patient & Disease", style = "padding-left:15px; color:#ccc"),
    selectInput("mut_type", "Mutation Type (Phenotype):",
                choices = c("SW (Null)", "SV (Simple Virilizing)", "NC (Non-Classical)"),
                selected = "SW (Null)"),
    numericInput("pt_weight", "Body Weight (kg):", value = 70, min = 5, max = 150),
    numericInput("pt_age", "Age (years):", value = 25, min = 1, max = 60),
    selectInput("pt_sex", "Sex:", choices = c("Female (46,XX)", "Male (46,XY)")),
    numericInput("sim_days", "Simulation duration (days):", value = 365, min = 30, max = 730),

    hr(),
    h5("Glucocorticoid Therapy", style = "padding-left:15px; color:#ccc"),
    numericInput("HC_dose", "Hydrocortisone (mg/day):", value = 20, min = 0, max = 60, step = 2.5),
    selectInput("HC_freq", "HC Frequency:", choices = c("TID (3×)"=3, "BID (2×)"=2, "QD (1×)"=1), selected = 3),
    numericInput("PRED_dose", "Prednisolone (mg/day):", value = 0, min = 0, max = 20, step = 0.5),
    selectInput("PRED_freq", "PRED Frequency:", choices = c("BID (2×)"=2, "QD (1×)"=1), selected = 2),

    hr(),
    h5("Mineralocorticoid Therapy", style = "padding-left:15px; color:#ccc"),
    numericInput("FC_dose", "Fludrocortisone (mcg/day):", value = 100, min = 0, max = 300, step = 50),

    hr(),
    h5("Novel Therapies (CRF1 Antagonists)", style = "padding-left:15px; color:#ccc"),
    numericInput("TILD_dose", "Tildacerfont (mg QD):", value = 0, min = 0, max = 300, step = 25),
    numericInput("CRINE_dose", "Crinecerfont (mg/day):", value = 0, min = 0, max = 400, step = 50),
    selectInput("CRINE_freq", "Crinecerfont Freq:", choices = c("BID (2×)"=2, "QD (1×)"=1), selected = 2),

    hr(),
    actionButton("run_sim", "Run Simulation", class = "btn-primary btn-block",
                 icon = icon("play"))
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #F8F9FA; }
      .box { border-radius: 8px; }
      .info-box { border-radius: 8px; }
    "))),

    tabItems(
      # =========================================================
      # TAB 1: Patient Profile
      # =========================================================
      tabItem("patient",
        fluidRow(
          box(width=12, status="primary", solidHeader=TRUE,
              title="CAH Patient Profile & Disease Background",
              fluidRow(
                column(4,
                  h4("Mutation Classification"),
                  tableOutput("mut_table")
                ),
                column(4,
                  h4("Expected Steroid Profile (Untreated)"),
                  tableOutput("baseline_steroids")
                ),
                column(4,
                  h4("Patient Summary"),
                  uiOutput("patient_summary")
                )
              ))
        ),
        fluidRow(
          box(width=6, status="info", solidHeader=TRUE,
              title="Pathophysiology – CYP21A2 Deficiency",
              tags$p("The CYP21A2 (21-hydroxylase) enzyme converts progesterone
              → 11-deoxycorticosterone (DOC) and 17-OHP → 11-deoxycortisol.
              In its absence, substrate accumulates and is shunted to
              androgen biosynthesis. ACTH rises due to absent cortisol feedback."),
              tags$ul(
                tags$li("Salt-Wasting (SW): >95% loss of function, aldosterone deficiency"),
                tags$li("Simple Virilizing (SV): ~1-2% residual, virilization only"),
                tags$li("Non-Classical (NC): ~20-50% residual, milder symptoms")
              ),
              tags$p(strong("Key drugs:"), " Hydrocortisone (glucocorticoid replacement),
              Fludrocortisone (mineralocorticoid replacement), Tildacerfont/Crinecerfont
              (CRF1 antagonists – reduce ACTH drive on adrenal).")
          ),
          box(width=6, status="warning", solidHeader=TRUE,
              title="Treatment Goals",
              tableOutput("treatment_goals")
          )
        )
      ),

      # =========================================================
      # TAB 2: Drug PK
      # =========================================================
      tabItem("pk",
        fluidRow(
          box(width=12, status="success", solidHeader=TRUE,
              title="Plasma Concentration-Time Profiles",
              plotlyOutput("pk_plot", height = "500px"))
        ),
        fluidRow(
          box(width=6, status="info",
              title="PK Parameters",
              tableOutput("pk_params_table")),
          box(width=6, status="info",
              title="Drug Interaction Notes",
              tags$ul(
                tags$li("Hydrocortisone t½ ≈ 1.5 h – requires TID dosing for physiologic coverage"),
                tags$li("Prednisolone t½ ≈ 2.5 h – BID dosing common; 4× GC potency"),
                tags$li("Dexamethasone t½ ≈ 3.8 h – 25× potency; bedtime dosing to suppress morning ACTH"),
                tags$li("Tildacerfont t½ ≈ 12-14 h – once-daily; IC50 ≈ 4 nM"),
                tags$li("Crinecerfont t½ ≈ 8-10 h – BID; IC50 ≈ 0.5 nM (10× more potent than TILD)")
              ))
        )
      ),

      # =========================================================
      # TAB 3: Steroid Biomarkers
      # =========================================================
      tabItem("steroids",
        fluidRow(
          valueBoxOutput("vb_17OHP"),
          valueBoxOutput("vb_ACTH"),
          valueBoxOutput("vb_A4")
        ),
        fluidRow(
          box(width=6, status="primary", solidHeader=TRUE,
              title="17-OHP (nmol/L) – Key Biomarker",
              plotlyOutput("plot_17OHP", height = "300px")),
          box(width=6, status="danger", solidHeader=TRUE,
              title="ACTH (pg/mL)",
              plotlyOutput("plot_ACTH", height = "300px"))
        ),
        fluidRow(
          box(width=6, status="warning", solidHeader=TRUE,
              title="Androstenedione (nmol/L)",
              plotlyOutput("plot_A4", height = "300px")),
          box(width=6, status="success", solidHeader=TRUE,
              title="Cortisol (nmol/L) – Replacement Adequacy",
              plotlyOutput("plot_CORT", height = "300px"))
        )
      ),

      # =========================================================
      # TAB 4: Clinical Endpoints
      # =========================================================
      tabItem("endpoints",
        fluidRow(
          valueBoxOutput("vb_HT"),
          valueBoxOutput("vb_BA"),
          valueBoxOutput("vb_BMD")
        ),
        fluidRow(
          box(width=6, status="info", solidHeader=TRUE,
              title="Height SDS (z-score)",
              plotlyOutput("plot_HT", height = "300px")),
          box(width=6, status="warning", solidHeader=TRUE,
              title="Bone Age Advancement (years)",
              plotlyOutput("plot_BA", height = "300px"))
        ),
        fluidRow(
          box(width=6, status="success", solidHeader=TRUE,
              title="Bone Mineral Density (normalized)",
              plotlyOutput("plot_BMD", height = "300px")),
          box(width=6, status="primary", solidHeader=TRUE,
              title="Renin / Mineralocorticoid Axis",
              plotlyOutput("plot_RENIN", height = "300px"))
        )
      ),

      # =========================================================
      # TAB 5: Scenario Comparison
      # =========================================================
      tabItem("scenarios",
        fluidRow(
          box(width=12, status="primary", solidHeader=TRUE,
              title="Treatment Scenario Comparison (Fixed 6 Arms)",
              plotlyOutput("scenario_plot", height = "650px"))
        ),
        fluidRow(
          box(width=12, status="info",
              title="12-Month Outcome Table",
              DTOutput("scenario_table"))
        )
      ),

      # =========================================================
      # TAB 6: Biomarker Dashboard
      # =========================================================
      tabItem("heatmap",
        fluidRow(
          box(width=6, status="primary", solidHeader=TRUE,
              title="Target Attainment (% time in target range)",
              plotlyOutput("target_plot", height = "400px")),
          box(width=6, status="info", solidHeader=TRUE,
              title="CRF1 Receptor Occupancy Over Time",
              plotlyOutput("crf1_plot", height = "400px"))
        ),
        fluidRow(
          box(width=12, status="success", solidHeader=TRUE,
              title="Biomarker Summary Table (Mean ± Last 30 days)",
              DTOutput("biomarker_table"))
        )
      )
    )
  )
)

# ---- SERVER -------------------------------------------------------------------
server <- function(input, output, session) {

  # Reactive simulation
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Running CAH simulation...", {
      run_sim(
        mut_type   = input$mut_type,
        HC_dose    = input$HC_dose,
        HC_freq    = as.numeric(input$HC_freq),
        PRED_dose  = input$PRED_dose,
        PRED_freq  = as.numeric(input$PRED_freq),
        FC_dose    = input$FC_dose,
        TILD_dose  = input$TILD_dose,
        CRINE_dose = input$CRINE_dose,
        CRINE_freq = as.numeric(input$CRINE_freq),
        days       = input$sim_days
      )
    })
  }, ignoreNULL = FALSE)

  # Scenario comparison (pre-computed fixed arms)
  scenarios_data <- reactive({
    withProgress(message = "Running scenario comparison (6 arms)...", {
      mut <- input$mut_type
      base_days <- 365

      s1 <- run_sim(mut, 0,  3, 0, 2, 0,   0,   0, 2, base_days) %>%
        mutate(Scenario = "Untreated")
      s2 <- run_sim(mut, 20, 3, 0, 2, 100, 0,   0, 2, base_days) %>%
        mutate(Scenario = "HC TID + FC")
      s3 <- run_sim(mut, 0,  3, 5, 2, 100, 0,   0, 2, base_days) %>%
        mutate(Scenario = "Prednisolone BID + FC")
      s4 <- run_sim("NC (Non-Classical)", 0, 3, 0, 2, 0, 0, 0, 2, base_days) %>%
        mutate(Scenario = "NC-CAH Untreated")
      s5 <- run_sim(mut, 15, 3, 0, 2, 100, 100, 0, 2, base_days) %>%
        mutate(Scenario = "Tildacerfont + HC + FC")
      s6 <- run_sim(mut, 15, 3, 0, 2, 100, 0, 200, 2, base_days) %>%
        mutate(Scenario = "Crinecerfont + HC + FC")

      bind_rows(s1, s2, s3, s4, s5, s6)
    })
  })

  # ---- TAB 1: Patient Profile ----
  output$mut_table <- renderTable({
    data.frame(
      Type     = c("SW (Salt-Wasting)", "SV (Simple Virilizing)", "NC (Non-Classical)"),
      Mutation = c("Null (I2G, del30kb, Q318X)", "I172N, P30L", "V281L, R339H"),
      CYP21A2  = c("<1%", "1-2%", "20-50%"),
      Aldosterone = c("Deficient", "Normal", "Normal"),
      Frequency= c("~75% of CAH", "~25%", "~0.1-0.2% pop")
    )
  })

  output$baseline_steroids <- renderTable({
    res <- switch(input$mut_type,
                  "SW (Null)"              = c("17-OHP: >120 nmol/L", "ACTH: >1000 pg/mL",
                                               "A4: >25 nmol/L", "Cortisol: <50 nmol/L",
                                               "Aldosterone: <50 pmol/L", "Renin: elevated 3-5×"),
                  "SV (Simple Virilizing)" = c("17-OHP: 30-120 nmol/L", "ACTH: 100-500 pg/mL",
                                               "A4: 15-30 nmol/L", "Cortisol: 50-150 nmol/L",
                                               "Aldosterone: normal", "Renin: normal"),
                  "NC (Non-Classical)"     = c("17-OHP: 6-36 nmol/L (stimulated)", "ACTH: 60-200 pg/mL",
                                               "A4: 7-15 nmol/L", "Cortisol: normal",
                                               "Aldosterone: normal", "Renin: normal")
    )
    data.frame(Biomarker = res)
  }, colnames = FALSE)

  output$patient_summary <- renderUI({
    tags$ul(
      tags$li(paste("Phenotype:", input$mut_type)),
      tags$li(paste("Age:", input$pt_age, "years")),
      tags$li(paste("Weight:", input$pt_weight, "kg")),
      tags$li(paste("Sex:", input$pt_sex)),
      tags$li(paste("HC dose:", input$HC_dose, "mg/day TID")),
      tags$li(paste("FC dose:", input$FC_dose, "mcg/day")),
      if (input$TILD_dose > 0) tags$li(paste("Tildacerfont:", input$TILD_dose, "mg QD")),
      if (input$CRINE_dose > 0) tags$li(paste("Crinecerfont:", input$CRINE_dose, "mg BID"))
    )
  })

  output$treatment_goals <- renderTable({
    data.frame(
      Biomarker = c("17-OHP", "ACTH", "Androstenedione", "Testosterone",
                    "Renin", "Height SDS"),
      Target    = c("<36 nmol/L (adults)", "<100 pg/mL",
                    "<7 nmol/L (adults)", "Age-appropriate",
                    "Normal range", "0 ± 2"),
      Rationale = c("Androgen suppression", "HPA suppression",
                    "Androgen excess marker", "Virilization marker",
                    "MC adequacy", "Normal linear growth")
    )
  })

  # ---- TAB 2: PK ----
  output$pk_plot <- renderPlotly({
    df <- sim_data()
    # Focus on day 1-7 for PK profile
    df7 <- filter(df, time_days <= 7)
    p <- ggplot(df7) +
      geom_line(aes(x = time_days, y = HC_nmol, color = "Hydrocortisone (nmol/L)"), linewidth=0.9) +
      geom_line(aes(x = time_days, y = TILD_conc_mgL * 1e4, color = "Tildacerfont ×10⁴ mg/L"), linewidth=0.9, linetype="dashed") +
      geom_line(aes(x = time_days, y = CRINE_conc_mgL * 1e6, color = "Crinecerfont ×10⁶ mg/L"), linewidth=0.9, linetype="dotted") +
      scale_color_manual(values = c("Hydrocortisone (nmol/L)"="#2171B5",
                                     "Tildacerfont ×10⁴ mg/L"="#FE9929",
                                     "Crinecerfont ×10⁶ mg/L"="#41AB5D")) +
      labs(x = "Time (days)", y = "Concentration", color = "Drug",
           title = "Drug PK – First Week") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_params_table <- renderTable({
    data.frame(
      Drug = c("Hydrocortisone", "Prednisolone", "Fludrocortisone",
               "Tildacerfont", "Crinecerfont"),
      `F (%)` = c(95, 82, 90, 65, 50),
      `t½ (h)` = c(1.5, 2.5, 3.5, 13, 9),
      `GC Potency` = c("1×", "4×", "—", "—", "—"),
      `IC50 CRF1R` = c("—", "—", "—", "4 nM", "0.5 nM")
    )
  })

  # ---- TAB 3: Steroid Biomarkers ----
  last_17OHP  <- reactive(round(mean(tail(sim_data()$OHP17, 100)), 1))
  last_ACTH   <- reactive(round(mean(tail(sim_data()$ACTH_pgmL, 100)), 1))
  last_A4     <- reactive(round(mean(tail(sim_data()$A4_nmol, 100)), 1))

  output$vb_17OHP <- renderValueBox({
    val <- last_17OHP()
    col <- if (val < 36) "green" else if (val < 100) "yellow" else "red"
    valueBox(paste(val, "nmol/L"), "17-OHP (target <36)", icon = icon("vial"), color = col)
  })
  output$vb_ACTH <- renderValueBox({
    val <- last_ACTH()
    col <- if (val < 100) "green" else if (val < 300) "yellow" else "red"
    valueBox(paste(val, "pg/mL"), "ACTH (target <100)", icon = icon("arrow-up"), color = col)
  })
  output$vb_A4 <- renderValueBox({
    val <- last_A4()
    col <- if (val < 7) "green" else if (val < 15) "yellow" else "red"
    valueBox(paste(val, "nmol/L"), "Androstenedione (target <7)", icon = icon("exclamation"), color = col)
  })

  mk_plotly <- function(df, y, ylab, target = NULL, color = "#2171B5") {
    p <- ggplot(df, aes_string(x = "time_days", y = y)) +
      geom_line(color = color, linewidth = 0.9) +
      labs(x = "Time (days)", y = ylab) +
      theme_bw()
    if (!is.null(target))
      p <- p + geom_hline(yintercept = target, linetype = "dashed", color = "red")
    ggplotly(p)
  }

  output$plot_17OHP <- renderPlotly(mk_plotly(sim_data(), "OHP17",
                                               "17-OHP (nmol/L)", 36, "#E6550D"))
  output$plot_ACTH  <- renderPlotly(mk_plotly(sim_data(), "ACTH_pgmL",
                                               "ACTH (pg/mL)", 100, "#CB181D"))
  output$plot_A4    <- renderPlotly(mk_plotly(sim_data(), "A4_nmol",
                                               "A4 (nmol/L)", 7, "#FE9929"))
  output$plot_CORT  <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time_days, y = CORTISOL)) +
      geom_line(color = "#2171B5", linewidth = 0.9) +
      geom_hline(yintercept = c(138, 690), linetype = "dashed", color = "gray50") +
      labs(x = "Time (days)", y = "Cortisol (nmol/L)") +
      theme_bw()
    ggplotly(p)
  })

  # ---- TAB 4: Clinical Endpoints ----
  output$vb_HT  <- renderValueBox({
    val <- round(tail(sim_data()$HEIGHT_SDS, 1), 2)
    col <- if (abs(val) < 2) "green" else "yellow"
    valueBox(paste(val, "SDS"), "Height SDS", icon = icon("child"), color = col)
  })
  output$vb_BA  <- renderValueBox({
    val <- round(tail(sim_data()$BONE_AGE, 1), 1)
    col <- if (val < 1) "green" else if (val < 3) "yellow" else "red"
    valueBox(paste(val, "years ahead"), "Bone Age Advance", icon = icon("bone"), color = col)
  })
  output$vb_BMD <- renderValueBox({
    val <- round(tail(sim_data()$BMD, 1), 3)
    col <- if (val > 0.9) "green" else if (val > 0.8) "yellow" else "red"
    valueBox(val, "BMD (normalized)", icon = icon("chart-bar"), color = col)
  })

  output$plot_HT  <- renderPlotly(mk_plotly(sim_data(), "HEIGHT_SDS",
                                              "Height SDS", 0, "#74C476"))
  output$plot_BA  <- renderPlotly(mk_plotly(sim_data(), "BONE_AGE",
                                              "Bone Age Adv (yr)", 0, "#FD8D3C"))
  output$plot_BMD <- renderPlotly(mk_plotly(sim_data(), "BMD",
                                              "BMD (normalized)", 1.0, "#9E9AC8"))
  output$plot_RENIN <- renderPlotly(mk_plotly(sim_data(), "RENIN",
                                               "Renin (normalized)", 1.0, "#6BAED6"))

  # ---- TAB 5: Scenario Comparison ----
  output$scenario_plot <- renderPlotly({
    df <- scenarios_data()
    colors6 <- c("Untreated"="#CB181D", "HC TID + FC"="#2171B5",
                  "Prednisolone BID + FC"="#238B45",
                  "NC-CAH Untreated"="#6A51A3",
                  "Tildacerfont + HC + FC"="#FE9929",
                  "Crinecerfont + HC + FC"="#41AB5D")

    p <- ggplot(df %>% filter(time_days <= 365),
                aes(x = time_days, y = OHP17, color = Scenario)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 36, linetype = "dashed", color = "gray40") +
      scale_color_manual(values = colors6) +
      scale_y_log10() +
      labs(x = "Time (days)", y = "17-OHP (nmol/L, log scale)",
           title = "Six Treatment Scenarios – 17-OHP Trajectory",
           color = "Treatment") +
      theme_bw(base_size = 12)
    ggplotly(p)
  })

  output$scenario_table <- renderDT({
    df <- scenarios_data() %>%
      filter(time_days >= 300) %>%
      group_by(Scenario) %>%
      summarise(
        `17-OHP (nmol/L)` = round(mean(OHP17), 1),
        `ACTH (pg/mL)`    = round(mean(ACTH_pgmL), 1),
        `A4 (nmol/L)`     = round(mean(A4_nmol), 2),
        `Cortisol (nmol/L)` = round(mean(CORTISOL), 1),
        `Height SDS`      = round(mean(HEIGHT_SDS), 2),
        `Bone Age Adv`    = round(mean(BONE_AGE), 2),
        `BMD`             = round(mean(BMD), 3),
        `Renin`           = round(mean(RENIN), 2),
        .groups = "drop"
      )
    datatable(df, options = list(pageLength = 10, scrollX = TRUE)) %>%
      formatStyle("17-OHP (nmol/L)", backgroundColor =
                    styleInterval(c(36, 100), c("#c8f7c5","#ffe599","#f7c6c6")))
  })

  # ---- TAB 6: Biomarker Dashboard ----
  output$target_plot <- renderPlotly({
    df <- sim_data()
    last30 <- filter(df, time_days >= max(time_days) - 30)
    targets <- list(
      `17-OHP <36` = mean(last30$OHP17 < 36) * 100,
      `ACTH <100 pg/mL` = mean(last30$ACTH_pgmL < 100) * 100,
      `A4 <7 nmol/L` = mean(last30$A4_nmol < 7) * 100,
      `Cortisol 138-690` = mean(last30$CORTISOL >= 138 & last30$CORTISOL <= 690) * 100,
      `Renin 0.5-2` = mean(last30$RENIN >= 0.5 & last30$RENIN <= 2) * 100,
      `BMD >0.85` = mean(last30$BMD > 0.85) * 100
    )
    tdf <- data.frame(Biomarker = names(targets),
                      `Pct_In_Target` = unlist(targets))
    p <- ggplot(tdf, aes(x = reorder(Biomarker, Pct_In_Target), y = Pct_In_Target,
                         fill = Pct_In_Target)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 80, linetype = "dashed", color = "green") +
      coord_flip() +
      scale_fill_gradient(low = "#CB181D", high = "#41AB5D") +
      labs(x = "", y = "% Time in Target Range",
           title = "Target Attainment (Last 30 days)") +
      theme_bw() + theme(legend.position = "none")
    ggplotly(p)
  })

  output$crf1_plot <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = time_days)) +
      geom_line(aes(y = TILD_pct, color = "Tildacerfont"), linewidth = 0.9) +
      geom_line(aes(y = CRINE_pct, color = "Crinecerfont"), linewidth = 0.9) +
      geom_line(aes(y = CRF1_block_pct, color = "Combined"), linetype = "dashed", linewidth = 0.9) +
      scale_color_manual(values = c("Tildacerfont"="#FE9929",
                                     "Crinecerfont"="#41AB5D",
                                     "Combined"="#CB181D")) +
      labs(x = "Time (days)", y = "CRF1R Occupancy (%)",
           color = "Agent", title = "CRF1 Receptor Occupancy Over Time") +
      theme_bw()
    ggplotly(p)
  })

  output$biomarker_table <- renderDT({
    df <- sim_data()
    last30 <- filter(df, time_days >= max(time_days) - 30)
    bm <- data.frame(
      Biomarker = c("17-OHP", "ACTH", "Androstenedione", "Testosterone",
                    "Cortisol", "Aldosterone", "Renin",
                    "Height SDS", "Bone Age Adv", "BMD", "CRF1 Block"),
      Unit = c("nmol/L","pg/mL","nmol/L","nmol/L","nmol/L","pmol/L",
               "norm","SDS","years","norm","%"),
      Mean_Last30 = round(c(mean(last30$OHP17), mean(last30$ACTH_pgmL),
                             mean(last30$A4_nmol), mean(last30$T_nmol),
                             mean(last30$CORTISOL), mean(last30$ALDOST),
                             mean(last30$RENIN), mean(last30$HEIGHT_SDS),
                             mean(last30$BONE_AGE), mean(last30$BMD),
                             mean(last30$CRF1_block_pct)), 2),
      Target = c("<36","<100","<7","Age-approp","138-690","Normal",
                 "0.5-2","±2","<1","≥0.9","—"),
      Status = c(
        if (mean(last30$OHP17) < 36) "✓ On Target" else "✗ Elevated",
        if (mean(last30$ACTH_pgmL) < 100) "✓ On Target" else "✗ Elevated",
        if (mean(last30$A4_nmol) < 7) "✓ On Target" else "✗ Elevated",
        if (mean(last30$T_nmol) < 2) "✓ On Target" else "✗ Elevated",
        if (mean(last30$CORTISOL) >= 138) "✓ On Target" else "✗ Low",
        if (mean(last30$ALDOST) > 100) "✓ On Target" else "✗ Low",
        if (mean(last30$RENIN) <= 2) "✓ On Target" else "✗ Elevated",
        if (abs(mean(last30$HEIGHT_SDS)) < 2) "✓ Normal" else "✗ Abnormal",
        if (mean(last30$BONE_AGE) < 1) "✓ Normal" else "✗ Elevated",
        if (mean(last30$BMD) > 0.9) "✓ Normal" else "✗ Low",
        "—"
      )
    )
    datatable(bm, options = list(pageLength = 15)) %>%
      formatStyle("Status",
                  backgroundColor = styleEqual(
                    c("✓ On Target","✓ Normal","✗ Elevated","✗ Low","✗ Abnormal"),
                    c("#c8f7c5","#c8f7c5","#f7c6c6","#f7c6c6","#ffe599")
                  ))
  })
}

# ---- Launch ----
shinyApp(ui = ui, server = server)
