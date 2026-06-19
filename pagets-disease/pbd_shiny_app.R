## ============================================================
## Paget's Disease of Bone — Interactive QSP Shiny Dashboard
## 파젯병 정량적 시스템 약리학 인터랙티브 대시보드
##
## Author: CCR (Claude Code Routine)
## Date  : 2026-06-19
##
## Tabs:
##   1. Patient Profile     — demographics & disease severity
##   2. PK Profile          — drug concentration-time curves
##   3. Biomarkers          — ALP, CTx, P1NP, TRAP5b
##   4. Clinical Endpoints  — bone quality, pain, fracture risk
##   5. Scenario Comparison — side-by-side treatment arms
##   6. Mechanism Explorer  — interactive pathway diagram
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)
library(plotly)

## ---- Embedded mrgsolve model code -------------------------
pbd_model_code <- '
$PARAM @annotated
CL_bp    = 5.0    : L/h     BP systemic clearance
V1_bp    = 8.0    : L       BP central volume
V2_bp    = 40.0   : L       BP peripheral volume
Q_bp     = 2.0    : L/h     BP intercompartment Q
Ka_bp    = 0.5    : h-1     Oral absorption (oral route)
F_bp     = 0.006  : -       Oral bioavailability
Kbone_on = 0.05   : h-1     Bone uptake rate
Kbone_off= 0.0001 : h-1     Bone release rate
CL_bone  = 0.002  : h-1     Bone-local clearance
Ka_den   = 0.003  : h-1     Denosumab SC absorption
F_den    = 0.62   : -       Denosumab bioavailability
CL_den   = 0.013  : L/h     Denosumab CL
V1_den   = 3.0    : L       Denosumab Vd
V2_den   = 2.5    : L       Denosumab peripheral Vd
Q_den    = 0.08   : L/h     Denosumab inter-Q
KD_den   = 0.008  : nM      Denosumab KD
RANKL_syn = 0.08  : nmol/L/h RANKL synthesis
RANKL_deg = 0.04  : h-1     RANKL degradation
OPG_syn   = 0.06  : nmol/L/h OPG synthesis
OPG_deg   = 0.03  : h-1     OPG degradation
PBD_fold  = 3.0   : -       RANKL/OPG elevation fold
OCP_syn   = 500   : cells/mL/h OCP production
OCP_deg   = 0.04  : h-1     OCP degradation
OC_diff   = 0.001 : h-1     OC differentiation
OC_max    = 2000  : cells/mL Max OC density
OCN_life  = 0.008 : h-1     OC apoptosis rate
OCN_base  = 300   : cells/mL Baseline OC
EC50_oc   = 1.0   : nmol/L  EC50 RANKL for OC diff
OBP_syn   = 600   : cells/mL/h OBP production
OBP_deg   = 0.03  : h-1     OBP degradation
OB_diff   = 0.0008: h-1     OB differentiation
OBN_life  = 0.006 : h-1     OB apoptosis rate
OBN_base  = 250   : cells/mL Baseline OB
OB_coupling = 1.5 : -       OC-OB coupling coefficient
CTX_prod  = 0.0002: ng/mL/cell/h CTx per OC
CTX_elim  = 0.15  : h-1     CTx elimination
P1NP_prod = 0.0003: ng/mL/cell/h P1NP per OB
P1NP_elim = 0.04  : h-1     P1NP elimination
ALP_prod  = 0.0006: IU/L/cell/h ALP per OB
ALP_elim  = 0.01  : h-1     ALP elimination
TRAP_prod = 0.0001: U/L/cell/h TRAP5b per OC
TRAP_elim = 0.08  : h-1     TRAP5b elimination
BMD_form  = 0.00001: g/cm2/h BMD formation rate
BMD_resorb= 0.00002: g/cm2/h BMD resorption rate
BMD_base  = 1.5   : g/cm2   BMD baseline pagetic
PAIN_base = 6.0   : VAS     Baseline pain
PAIN_slope= 0.5   : -       Pain ALP sensitivity
PAIN_recov= 0.001 : h-1     Pain recovery
Emax_bp   = 0.95  : -       Max BP inhibition
EC50_bp   = 0.01  : nmol/L  EC50 bone BP
Emax_den  = 0.90  : -       Max denosumab RANKL-neut
EC50_rankl_den = 0.05 : nmol/L Den EC50

$INIT
DRUG1 = 0.0  DRUG2 = 0.0  DRUGB = 0.0
DEN_SC = 0.0 DEN_C = 0.0  DEN_P = 0.0
OCP = 12500.0 OCN = 300.0
OBP = 20000.0 OBN = 250.0
RANKL_f = 2.0 OPG_f = 1.5
CTX = 3.0  P1NP = 80.0  ALP_b = 350.0
TRAP_s = 5.5 BMD = 1.5  PAIN = 6.0

$ODE
double C1_bp  = DRUG1 / V1_bp;
double C2_bp  = DRUG2 / V2_bp;
double C_den  = DEN_C / V1_den;
double CP_den = DEN_P / V2_den;
double Den_RANKL_inh = Emax_den * DEN_C / (EC50_rankl_den + DEN_C);
double RANKL_free_eff = RANKL_f * (1 - Den_RANKL_inh);
double OC_diff_stim = OC_diff * RANKL_free_eff / (EC50_oc + RANKL_free_eff);
double Cbone_eff = DRUGB / (DRUGB + 0.001);
double BP_oc_inh = Emax_bp * Cbone_eff / (EC50_bp + Cbone_eff);
dxdt_DRUG1  = -CL_bp*C1_bp - Q_bp*(C1_bp-C2_bp) - Kbone_on*DRUG1 + Kbone_off*DRUGB;
dxdt_DRUG2  =  Q_bp*(C1_bp - C2_bp);
dxdt_DRUGB  =  Kbone_on*DRUG1 - Kbone_off*DRUGB - CL_bone*DRUGB;
dxdt_DEN_SC = -Ka_den * F_den * DEN_SC;
dxdt_DEN_C  =  Ka_den*F_den*DEN_SC - CL_den*C_den - Q_den*(C_den-CP_den) - KD_den*DEN_C*RANKL_f;
dxdt_DEN_P  =  Q_den*(C_den - CP_den);
dxdt_RANKL_f= RANKL_syn*PBD_fold - RANKL_deg*RANKL_f - KD_den*DEN_C*RANKL_f/100.0;
dxdt_OPG_f  = OPG_syn - OPG_deg*OPG_f;
dxdt_OCP    = OCP_syn - OCP_deg*OCP - OC_diff_stim*OCP;
dxdt_OCN    = OC_diff_stim*OCP - OCN_life*OCN*(1 + BP_oc_inh*5.0);
double OB_drive = OB_coupling * OCN / OBN_base;
dxdt_OBP    = OBP_syn - OBP_deg*OBP - OB_diff*OBP;
dxdt_OBN    = OB_diff*OBP*OB_drive - OBN_life*OBN;
dxdt_CTX    = CTX_prod*OCN - CTX_elim*CTX;
dxdt_P1NP   = P1NP_prod*OBN - P1NP_elim*P1NP;
dxdt_ALP_b  = ALP_prod*OBN - ALP_elim*ALP_b;
dxdt_TRAP_s = TRAP_prod*OCN - TRAP_elim*TRAP_s;
dxdt_BMD    = BMD_form*(OBN/OBN_base) - BMD_resorb*(OCN/OCN_base);
dxdt_PAIN   = PAIN_slope*((ALP_b-100)/300.0)*0.001 - PAIN_recov*(PAIN-1.0);

$TABLE
double C1_out = DRUG1/V1_bp;
double DEN_out= DEN_C/V1_den;
double RANKL_OPG_r = RANKL_f/(OPG_f+0.001);
double ALP_pct = ALP_b/350.0*100.0;
double FractureRisk = (BMD < 1.2) ? 3.0 : (BMD < 1.4 ? 2.0 : 1.0);

$CAPTURE C1_out DEN_out RANKL_OPG_r ALP_pct FractureRisk
'

## ---- Build model at startup ----
pbd_mod <- mcode("pbd_shiny", pbd_model_code, quiet = TRUE)

## ---- Helper: run simulation ----
run_sim <- function(mod, drug, dose_mg, route, duration_days,
                    pbd_fold, ob_coupling, pain_base,
                    end_days = 365) {
  t_grid <- seq(0, end_days * 24, by = 24)

  params_list <- list(
    PBD_fold    = pbd_fold,
    OB_coupling = ob_coupling,
    PAIN_base   = pain_base
  )

  dose_units <- dose_mg * 1000   # mg -> mcg (arbitrary unit scale)

  if (drug == "None") {
    ev_obj <- ev(time = 0, amt = 0, cmt = 1)
  } else if (drug == "Zoledronate") {
    ev_obj <- ev(time = 0, amt = dose_units, cmt = 1, rate = dose_units / 0.25)
  } else if (drug == "Pamidronate") {
    n_inf <- max(1, ceiling(duration_days / 1))
    ev_obj <- ev(time = 0, amt = dose_units, cmt = 1, rate = dose_units / 4,
                 addl = n_inf - 1, ii = 24)
  } else if (drug == "Alendronate") {
    n_doses <- floor(duration_days)
    ev_obj <- ev(time = 0, amt = dose_units, cmt = 1,
                 addl = n_doses - 1, ii = 24)
  } else if (drug == "Denosumab") {
    n_doses <- max(1, floor(end_days / 180))
    ev_obj <- ev(time = 0, amt = dose_mg, cmt = 4,
                 addl = n_doses - 1, ii = 180 * 24)
  }

  mod %>%
    param(params_list) %>%
    mrgsim(ev = ev_obj, tgrid = t_grid) %>%
    as_tibble() %>%
    mutate(day = time / 24)
}

## ============================================================
## UI
## ============================================================

ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(
    title = "PBD QSP Dashboard",
    titleWidth = 280
  ),
  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "patient",    icon = icon("user-md")),
      menuItem("PK Profile",          tabName = "pk",         icon = icon("flask")),
      menuItem("Biomarkers",          tabName = "biomarkers", icon = icon("chart-line")),
      menuItem("Clinical Endpoints",  tabName = "clinical",   icon = icon("bone")),
      menuItem("Scenario Comparison", tabName = "compare",    icon = icon("balance-scale")),
      menuItem("Mechanism Explorer",  tabName = "mechanism",  icon = icon("project-diagram"))
    ),
    hr(),
    h5("  Treatment Settings", style = "color:#ccc; padding-left:15px;"),
    selectInput("drug", "Drug",
                choices = c("None", "Zoledronate", "Pamidronate",
                            "Alendronate", "Denosumab"),
                selected = "Zoledronate"),
    sliderInput("dose_mg", "Dose (mg)", min = 0.5, max = 100, value = 5, step = 0.5),
    sliderInput("duration_d", "Duration (days)", min = 1, max = 365, value = 1),
    sliderInput("end_days", "Simulation (days)", min = 90, max = 730, value = 365),
    hr(),
    h5("  Disease Parameters", style = "color:#ccc; padding-left:15px;"),
    sliderInput("pbd_fold", "RANKL/OPG Fold (disease severity)",
                min = 1, max = 8, value = 3, step = 0.5),
    sliderInput("ob_coupling", "OC→OB Coupling",
                min = 0.5, max = 3.0, value = 1.5, step = 0.1),
    sliderInput("pain_base", "Baseline Pain (VAS)",
                min = 1, max = 10, value = 6, step = 0.5),
    actionButton("run_sim", "Run Simulation",
                 class = "btn-primary btn-block", style = "margin:10px;")
  ),
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper { background-color: #f8f9fa; }
        .box { border-radius: 8px; }
        .value-box .inner { padding: 10px 15px; }
      "))
    ),
    tabItems(
      ## =============================================
      ## TAB 1: PATIENT PROFILE
      ## =============================================
      tabItem(
        tabName = "patient",
        h2("Patient Profile & Disease Overview"),
        fluidRow(
          valueBoxOutput("vbox_alp",    width = 3),
          valueBoxOutput("vbox_ctx",    width = 3),
          valueBoxOutput("vbox_pain",   width = 3),
          valueBoxOutput("vbox_fxrisk", width = 3)
        ),
        fluidRow(
          box(title = "Paget's Disease of Bone — Overview", width = 6,
              status = "purple", solidHeader = TRUE,
              HTML("
                <b>Paget's disease of bone (PBD)</b> is the second most common
                metabolic bone disease after osteoporosis. It is characterised by
                focal areas of greatly accelerated bone remodelling, with markedly
                increased and dysregulated osteoclastic bone resorption coupled
                with compensatory excessive osteoblastic bone formation.
                <br><br>
                <b>Key features:</b>
                <ul>
                  <li>Abnormal giant osteoclasts (up to 100+ nuclei; normal 3–5)</li>
                  <li>RANKL/OPG ratio markedly elevated → OC hyperactivation</li>
                  <li>Woven, disorganised bone → mosaic cement-line pattern</li>
                  <li>Bone ALP: primary biochemical marker (↑2–100× normal)</li>
                  <li>Paramyxovirus inclusion bodies in pagetic OC (MV, CDV)</li>
                  <li>SQSTM1/p62 mutations in ~25–40% of familial cases</li>
                </ul>
                <b>Primary treatment goal:</b> biochemical remission (ALP normalisation)
              ")
          ),
          box(title = "Disease Severity — Current Parameters", width = 6,
              status = "info", solidHeader = TRUE,
              tableOutput("disease_summary_table")
          )
        ),
        fluidRow(
          box(title = "Baseline Disease Activity (ALP distribution in PBD)",
              width = 12, status = "warning", solidHeader = TRUE,
              plotlyOutput("alp_distribution_plot", height = "300px"))
        )
      ),

      ## =============================================
      ## TAB 2: PK PROFILE
      ## =============================================
      tabItem(
        tabName = "pk",
        h2("Pharmacokinetic Profile"),
        fluidRow(
          box(title = "Drug Plasma Concentration–Time", width = 8,
              status = "blue", solidHeader = TRUE,
              plotlyOutput("pk_plot", height = "400px")),
          box(title = "PK Parameters", width = 4,
              status = "info", solidHeader = TRUE,
              tableOutput("pk_params_table")
          )
        ),
        fluidRow(
          box(title = "Bone-Bound Drug Accumulation", width = 6,
              status = "purple", solidHeader = TRUE,
              plotlyOutput("bone_pk_plot", height = "300px")),
          box(title = "RANKL/OPG Ratio over Time", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("rankl_opg_plot", height = "300px"))
        )
      ),

      ## =============================================
      ## TAB 3: BIOMARKERS
      ## =============================================
      tabItem(
        tabName = "biomarkers",
        h2("Bone Turnover Biomarkers"),
        fluidRow(
          box(title = "Bone ALP — Primary Endpoint", width = 6,
              status = "purple", solidHeader = TRUE,
              plotlyOutput("alp_plot", height = "320px")),
          box(title = "Serum CTx — Resorption Marker", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("ctx_plot", height = "320px"))
        ),
        fluidRow(
          box(title = "Serum P1NP — Formation Marker", width = 6,
              status = "success", solidHeader = TRUE,
              plotlyOutput("p1np_plot", height = "320px")),
          box(title = "TRAP5b — Osteoclast Activity", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("trap_plot", height = "320px"))
        ),
        fluidRow(
          box(title = "Biomarker Summary Table (Day 0, 30, 90, 180, 365)",
              width = 12, status = "info", solidHeader = TRUE,
              DTOutput("biomarker_table"))
        )
      ),

      ## =============================================
      ## TAB 4: CLINICAL ENDPOINTS
      ## =============================================
      tabItem(
        tabName = "clinical",
        h2("Clinical Endpoints & Bone Quality"),
        fluidRow(
          box(title = "BMD at Pagetic Site", width = 6,
              status = "success", solidHeader = TRUE,
              plotlyOutput("bmd_plot", height = "320px")),
          box(title = "Pain VAS Score", width = 6,
              status = "danger", solidHeader = TRUE,
              plotlyOutput("pain_plot", height = "320px"))
        ),
        fluidRow(
          box(title = "Active Osteoclast Count", width = 6,
              status = "warning", solidHeader = TRUE,
              plotlyOutput("ocn_plot", height = "320px")),
          box(title = "Active Osteoblast Count", width = 6,
              status = "info", solidHeader = TRUE,
              plotlyOutput("obn_plot", height = "320px"))
        ),
        fluidRow(
          box(title = "Fracture Risk Classification",
              width = 12, status = "danger", solidHeader = TRUE,
              fluidRow(
                column(6, plotlyOutput("fracture_risk_plot", height = "300px")),
                column(6, HTML("
                  <h4>Fracture Risk Criteria in PBD</h4>
                  <table class='table table-bordered'>
                    <tr><th>BMD (g/cm²)</th><th>Risk Level</th><th>Action</th></tr>
                    <tr class='danger'><td>&lt;1.2</td><td>High (3×)</td><td>Immediate treatment + orthopedic referral</td></tr>
                    <tr class='warning'><td>1.2–1.4</td><td>Moderate (2×)</td><td>Antiresorptive + weight-bearing protection</td></tr>
                    <tr class='success'><td>&gt;1.4</td><td>Normal</td><td>Monitor, treat if symptomatic</td></tr>
                  </table>
                  <p><i>Note: BMD in PBD reflects quantity not quality; pagetic bone has
                  increased volume but disorganised microarchitecture.</i></p>
                "))
              )
          )
        )
      ),

      ## =============================================
      ## TAB 5: SCENARIO COMPARISON
      ## =============================================
      tabItem(
        tabName = "compare",
        h2("Multi-Scenario Comparison"),
        fluidRow(
          box(title = "Select Scenarios to Compare", width = 12,
              status = "primary", solidHeader = TRUE,
              fluidRow(
                column(3,
                  checkboxGroupInput("comp_scenarios", "Treatment Arms:",
                    choices = c("No Treatment", "Zoledronate 5mg", "Pamidronate 60mg",
                                "Alendronate 40mg", "Denosumab 60mg"),
                    selected = c("No Treatment", "Zoledronate 5mg", "Denosumab 60mg"))
                ),
                column(3,
                  sliderInput("comp_pbd_fold", "RANKL/OPG Fold", 1, 8, 3, 0.5),
                  sliderInput("comp_days", "Days to simulate", 90, 730, 365)
                ),
                column(3,
                  radioButtons("comp_endpoint", "Primary Endpoint:",
                    choices = c("Bone ALP" = "ALP_b",
                                "Serum CTx" = "CTX",
                                "P1NP"      = "P1NP",
                                "Pain VAS"  = "PAIN"),
                    selected = "ALP_b")
                ),
                column(3,
                  actionButton("run_compare", "Compare Scenarios",
                               class = "btn-success btn-lg",
                               style = "margin-top:25px;")
                )
              )
          )
        ),
        fluidRow(
          box(title = "Endpoint Trajectory Comparison", width = 8,
              status = "primary", solidHeader = TRUE,
              plotlyOutput("compare_plot", height = "420px")),
          box(title = "At Day 365 Summary", width = 4,
              status = "info", solidHeader = TRUE,
              DTOutput("compare_table"))
        ),
        fluidRow(
          box(title = "Responder Analysis (ALP normalisation rate)",
              width = 12, status = "success", solidHeader = TRUE,
              plotlyOutput("responder_plot", height = "300px"))
        )
      ),

      ## =============================================
      ## TAB 6: MECHANISM EXPLORER
      ## =============================================
      tabItem(
        tabName = "mechanism",
        h2("Mechanism Explorer — Pagetic Pathways"),
        fluidRow(
          box(title = "Mechanistic Map (SVG)", width = 8,
              status = "purple", solidHeader = TRUE,
              tags$iframe(
                src = "pbd_qsp_model.svg",
                width = "100%", height = "600px",
                style = "border:none;"
              )
          ),
          box(title = "Pathway Summary", width = 4,
              status = "info", solidHeader = TRUE,
              HTML("
                <h4>Key Mechanistic Nodes</h4>
                <dl>
                  <dt>RANKL/OPG Axis</dt>
                  <dd>Central regulator of osteoclastogenesis. RANKL/OPG ratio ↑3–6× in PBD drives giant OC formation.</dd>
                  <dt>NF-κB / NFATc1</dt>
                  <dd>Downstream of RANK; NFATc1 is the master transcription factor for osteoclast differentiation. p62/SQSTM1 mutations enhance NF-κB activation.</dd>
                  <dt>FPP Synthase (Bisphosphonate target)</dt>
                  <dd>Inhibition of FPP synthase blocks protein prenylation (Rab, Rho GTPases) → cytoskeletal disruption → OC apoptosis.</dd>
                  <dt>Wnt / β-catenin</dt>
                  <dd>Drives osteoblast differentiation; Runx2/Osterix are downstream effectors. Sclerostin ↓ in PBD allows exaggerated bone formation.</dd>
                  <dt>OC–OB Coupling</dt>
                  <dd>TGF-β1, IGF-1, BMP-2 released from resorbed matrix drive compensatory (but disorganised) bone formation.</dd>
                  <dt>Inclusion Bodies</dt>
                  <dd>Paramyxovirus (measles) nucleocapsid proteins found in pagetic OC; may enhance RANKL signalling and fusion.</dd>
                </dl>
              ")
          )
        ),
        fluidRow(
          box(title = "Drug Mechanism Summary", width = 12,
              status = "warning", solidHeader = TRUE,
              fluidRow(
                column(3, wellPanel(
                  h4("Nitrogen-containing BPs"),
                  p("Zoledronate, pamidronate, alendronate, risedronate"),
                  tags$ul(
                    tags$li("Bind hydroxyapatite with high affinity"),
                    tags$li("Taken up by osteoclasts during bone resorption"),
                    tags$li("Inhibit FPP synthase → block protein prenylation"),
                    tags$li("Cytoskeletal disruption → enhanced OC apoptosis"),
                    tags$li("Zoledronate: single IV 5 mg achieves >90% ALP reduction in 6 months")
                  )
                )),
                column(3, wellPanel(
                  h4("Denosumab"),
                  p("Anti-RANKL fully human monoclonal antibody (IgG2)"),
                  tags$ul(
                    tags$li("Binds sRANKL/mRANKL with picomolar affinity (KD ~0.001 nM)"),
                    tags$li("Prevents RANK–RANKL interaction"),
                    tags$li("Rapid suppression of CTx within 1–3 days"),
                    tags$li("SC 60 mg Q6M → sustained OC suppression"),
                    tags$li("Reversible: rebound effect on RANKL cessation")
                  )
                )),
                column(3, wellPanel(
                  h4("Calcitonin (historical)"),
                  p("Salmon calcitonin — now largely replaced by BPs"),
                  tags$ul(
                    tags$li("Binds calcitonin receptor on OC surface"),
                    tags$li("Inhibits OC motility and ruffled border"),
                    tags$li("Short-lived effect (escape phenomenon)"),
                    tags$li("Mainly for acute hypercalcaemia/pain")
                  )
                )),
                column(3, wellPanel(
                  h4("Treatment Goals"),
                  p("Biochemical remission = primary objective"),
                  tags$ul(
                    tags$li("ALP normalisation (< ULN)"),
                    tags$li("CTx < 0.6 ng/mL"),
                    tags$li("P1NP < 40 ng/mL"),
                    tags$li("Pain ≤ 2 VAS"),
                    tags$li("Prevention of complications: fracture, deformity, sarcoma")
                  )
                ))
              )
          )
        )
      )  # end tabItem mechanism
    )  # end tabItems
  )  # end dashboardBody
)  # end dashboardPage

## ============================================================
## SERVER
## ============================================================

server <- function(input, output, session) {

  ## ---- Reactive: run simulation ----
  sim_result <- eventReactive(input$run_sim, {
    withProgress(message = "Running QSP simulation...", {
      run_sim(
        mod          = pbd_mod,
        drug         = input$drug,
        dose_mg      = input$dose_mg,
        route        = ifelse(input$drug %in% c("Alendronate"), "oral", "IV"),
        duration_days= input$duration_d,
        pbd_fold     = input$pbd_fold,
        ob_coupling  = input$ob_coupling,
        pain_base    = input$pain_base,
        end_days     = input$end_days
      )
    })
  }, ignoreNULL = FALSE)

  ## Run with defaults on startup
  observe({
    if (is.null(isolate(sim_result()))) input$run_sim
  })

  ## ---- Value Boxes ----
  output$vbox_alp <- renderValueBox({
    df <- sim_result()
    alp_last <- round(tail(df$ALP_b, 1), 0)
    color <- if (alp_last < 120) "green" else if (alp_last < 300) "yellow" else "red"
    valueBox(paste0(alp_last, " IU/L"), "Bone ALP", icon = icon("flask"), color = color)
  })

  output$vbox_ctx <- renderValueBox({
    df <- sim_result()
    ctx_last <- round(tail(df$CTX, 1), 2)
    color <- if (ctx_last < 0.6) "green" else if (ctx_last < 2) "yellow" else "red"
    valueBox(paste0(ctx_last, " ng/mL"), "Serum CTx", icon = icon("vial"), color = color)
  })

  output$vbox_pain <- renderValueBox({
    df <- sim_result()
    pain_last <- round(tail(df$PAIN, 1), 1)
    color <- if (pain_last <= 3) "green" else if (pain_last <= 6) "yellow" else "red"
    valueBox(paste0(pain_last, " VAS"), "Pain Score", icon = icon("heartbeat"), color = color)
  })

  output$vbox_fxrisk <- renderValueBox({
    df <- sim_result()
    fx_last <- tail(df$FractureRisk, 1)
    risk_label <- if (fx_last >= 3) "High Risk" else if (fx_last >= 2) "Moderate" else "Low Risk"
    color <- if (fx_last >= 3) "red" else if (fx_last >= 2) "yellow" else "green"
    valueBox(risk_label, "Fracture Risk", icon = icon("bone"), color = color)
  })

  ## ---- Disease summary table ----
  output$disease_summary_table <- renderTable({
    tibble(
      Parameter   = c("RANKL/OPG Fold", "OC–OB Coupling",
                      "Baseline ALP", "Baseline CTx", "Baseline Pain"),
      Value       = c(input$pbd_fold, input$ob_coupling,
                      "350 IU/L", "3.0 ng/mL", paste(input$pain_base, "VAS")),
      Status      = c(
        ifelse(input$pbd_fold >= 3, "Elevated (active PBD)", "Mild elevation"),
        ifelse(input$ob_coupling >= 1.5, "Normal coupling", "Impaired coupling"),
        "Elevated (PBD active)", "Elevated (>6×)", "Moderate-severe"
      )
    )
  })

  ## ---- ALP distribution plot ----
  output$alp_distribution_plot <- renderPlotly({
    set.seed(123)
    alp_vals <- c(
      rnorm(100, mean = 80,   sd = 20),    # normal controls
      rnorm(200, mean = 250,  sd = 100),   # mild PBD
      rnorm(150, mean = 600,  sd = 200),   # moderate PBD
      rnorm(80,  mean = 1500, sd = 500)    # severe PBD
    )
    df_alp <- tibble(ALP = pmax(30, alp_vals),
                     Group = c(rep("Normal", 100), rep("Mild PBD", 200),
                               rep("Moderate PBD", 150), rep("Severe PBD", 80)))
    p <- ggplot(df_alp, aes(ALP, fill = Group)) +
      geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
      geom_vline(xintercept = 120, color = "red", linetype = "dashed") +
      scale_fill_manual(values = c("Normal"="#4CAF50","Mild PBD"="#FFC107",
                                   "Moderate PBD"="#FF5722","Severe PBD"="#B71C1C")) +
      scale_x_log10() +
      labs(x = "Bone ALP (IU/L, log scale)", y = "Count") +
      theme_minimal()
    ggplotly(p)
  })

  ## ---- PK plot ----
  output$pk_plot <- renderPlotly({
    df <- sim_result()
    col_pk <- if (input$drug == "Denosumab") "DEN_out" else "C1_out"
    y_lab  <- if (input$drug == "Denosumab") "Denosumab (mg/L)" else "BP Plasma (nmol/L)"
    p <- ggplot(df, aes(day, .data[[col_pk]])) +
      geom_line(color = "#1E88E5", linewidth = 1.2) +
      labs(x = "Day", y = y_lab,
           title = paste(input$drug, "— Plasma Concentration")) +
      theme_classic()
    ggplotly(p)
  })

  output$bone_pk_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(day, DRUGB)) +
      geom_line(color = "#7B1FA2", linewidth = 1.2) +
      labs(x = "Day", y = "Bone-Bound Drug (nmol/g)",
           title = "Bone-Bound Drug Accumulation") +
      theme_classic()
    ggplotly(p)
  })

  output$rankl_opg_plot <- renderPlotly({
    df <- sim_result()
    p <- ggplot(df, aes(day, RANKL_OPG_r)) +
      geom_line(color = "#E65100", linewidth = 1.2) +
      geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray") +
      labs(x = "Day", y = "RANKL / OPG Ratio",
           title = "RANKL/OPG Ratio") +
      theme_classic()
    ggplotly(p)
  })

  output$pk_params_table <- renderTable({
    tibble(
      Parameter   = c("Drug", "Route", "Dose (mg)", "CL (L/h)", "V1 (L)", "Bone t½"),
      Value       = c(input$drug,
                      ifelse(input$drug %in% c("Alendronate","Risedronate"), "Oral", "IV/SC"),
                      as.character(input$dose_mg), "5.0", "8.0", "~10 years")
    )
  })

  ## ---- Biomarker plots ----
  make_bio_plot <- function(df, yvar, ylabel, color, ref_line = NULL) {
    p <- ggplot(df, aes(day, .data[[yvar]])) +
      geom_line(color = color, linewidth = 1.2)
    if (!is.null(ref_line)) {
      p <- p + geom_hline(yintercept = ref_line, linetype = "dashed", color = "gray50")
    }
    p + labs(x = "Day", y = ylabel) + theme_classic()
  }

  output$alp_plot   <- renderPlotly({ ggplotly(make_bio_plot(sim_result(), "ALP_b", "Bone ALP (IU/L)", "#7B1FA2", 120)) })
  output$ctx_plot   <- renderPlotly({ ggplotly(make_bio_plot(sim_result(), "CTX",   "CTx (ng/mL)",     "#C62828", 0.6)) })
  output$p1np_plot  <- renderPlotly({ ggplotly(make_bio_plot(sim_result(), "P1NP",  "P1NP (ng/mL)",    "#2E7D32", 40)) })
  output$trap_plot  <- renderPlotly({ ggplotly(make_bio_plot(sim_result(), "TRAP_s","TRAP5b (U/L)",     "#F57F17", 3.0)) })

  output$biomarker_table <- renderDT({
    df <- sim_result()
    df %>%
      filter(day %in% c(0, 30, 90, 180, 365)) %>%
      transmute(
        Day      = day,
        `ALP (IU/L)`   = round(ALP_b, 1),
        `CTx (ng/mL)`  = round(CTX,   2),
        `P1NP (ng/mL)` = round(P1NP,  1),
        `TRAP5b (U/L)` = round(TRAP_s, 2),
        `ALP % of baseline` = round(ALP_pct, 1)
      ) %>%
      datatable(options = list(pageLength = 6, searching = FALSE),
                rownames = FALSE) %>%
      formatStyle("ALP (IU/L)",
        backgroundColor = styleInterval(c(120, 300), c("#C8E6C9","#FFF9C4","#FFCDD2")))
  })

  ## ---- Clinical Endpoint plots ----
  output$bmd_plot <- renderPlotly({
    p <- ggplot(sim_result(), aes(day, BMD)) +
      geom_line(color = "#388E3C", linewidth = 1.2) +
      geom_hline(yintercept = 1.4, linetype = "dashed", color = "orange", alpha = 0.7) +
      geom_hline(yintercept = 1.2, linetype = "dashed", color = "red", alpha = 0.7) +
      labs(x = "Day", y = "BMD (g/cm²)", title = "BMD — Pagetic Site") +
      theme_classic()
    ggplotly(p)
  })

  output$pain_plot <- renderPlotly({
    p <- ggplot(sim_result(), aes(day, PAIN)) +
      geom_line(color = "#C62828", linewidth = 1.2) +
      scale_y_continuous(limits = c(0, 10)) +
      geom_hline(yintercept = 3, linetype = "dashed", color = "orange") +
      labs(x = "Day", y = "Pain VAS (0–10)", title = "Pain Score") +
      theme_classic()
    ggplotly(p)
  })

  output$ocn_plot <- renderPlotly({
    p <- ggplot(sim_result(), aes(day, OCN)) +
      geom_line(color = "#E65100", linewidth = 1.2) +
      geom_hline(yintercept = 100, linetype = "dashed", color = "gray50") +
      labs(x = "Day", y = "OC Count (cells/mL)", title = "Active Osteoclast Number") +
      theme_classic()
    ggplotly(p)
  })

  output$obn_plot <- renderPlotly({
    p <- ggplot(sim_result(), aes(day, OBN)) +
      geom_line(color = "#1565C0", linewidth = 1.2) +
      geom_hline(yintercept = 80, linetype = "dashed", color = "gray50") +
      labs(x = "Day", y = "OB Count (cells/mL)", title = "Active Osteoblast Number") +
      theme_classic()
    ggplotly(p)
  })

  output$fracture_risk_plot <- renderPlotly({
    p <- ggplot(sim_result(), aes(day, FractureRisk)) +
      geom_line(color = "#B71C1C", linewidth = 1.4) +
      scale_y_continuous(breaks = c(1, 2, 3), labels = c("Low", "Moderate", "High"),
                         limits = c(0.5, 3.5)) +
      labs(x = "Day", y = "Fracture Risk Level", title = "Fracture Risk Trajectory") +
      theme_classic()
    ggplotly(p)
  })

  ## ---- Scenario Comparison ----
  scenario_map <- list(
    "No Treatment"     = list(drug = "None",        dose = 0,  dur = 1),
    "Zoledronate 5mg"  = list(drug = "Zoledronate", dose = 5,  dur = 1),
    "Pamidronate 60mg" = list(drug = "Pamidronate", dose = 60, dur = 3),
    "Alendronate 40mg" = list(drug = "Alendronate", dose = 40, dur = 180),
    "Denosumab 60mg"   = list(drug = "Denosumab",   dose = 60, dur = 1)
  )

  compare_result <- eventReactive(input$run_compare, {
    scen_list <- input$comp_scenarios
    pbd_f     <- input$comp_pbd_fold
    days      <- input$comp_days

    results <- lapply(scen_list, function(sc) {
      cfg <- scenario_map[[sc]]
      run_sim(pbd_mod, cfg$drug, cfg$dose, "IV", cfg$dur,
              pbd_f, 1.5, 6.0, days) %>%
        mutate(Scenario = sc)
    })
    bind_rows(results)
  })

  compare_colors <- c(
    "No Treatment"     = "#E53935",
    "Zoledronate 5mg"  = "#1E88E5",
    "Pamidronate 60mg" = "#43A047",
    "Alendronate 40mg" = "#FB8C00",
    "Denosumab 60mg"   = "#8E24AA"
  )

  output$compare_plot <- renderPlotly({
    req(compare_result())
    df   <- compare_result()
    yvar <- input$comp_endpoint
    ylabs <- c(ALP_b = "Bone ALP (IU/L)", CTX = "CTx (ng/mL)",
                P1NP = "P1NP (ng/mL)",    PAIN = "Pain VAS")
    p <- ggplot(df, aes(day, .data[[yvar]], color = Scenario)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = compare_colors) +
      labs(x = "Day", y = ylabs[yvar]) +
      theme_classic() +
      theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$compare_table <- renderDT({
    req(compare_result())
    compare_result() %>%
      filter(day == max(day)) %>%
      group_by(Scenario) %>%
      summarize(
        `ALP (IU/L)` = round(mean(ALP_b), 0),
        `CTx`        = round(mean(CTX), 2),
        `P1NP`       = round(mean(P1NP), 0),
        `Pain VAS`   = round(mean(PAIN), 1),
        .groups = "drop"
      ) %>%
      datatable(options = list(pageLength = 8, searching = FALSE), rownames = FALSE)
  })

  output$responder_plot <- renderPlotly({
    req(compare_result())
    df_resp <- compare_result() %>%
      group_by(Scenario, day) %>%
      summarize(ALP_norm_frac = mean(ALP_b < 120) * 100, .groups = "drop")
    p <- ggplot(df_resp, aes(day, ALP_norm_frac, color = Scenario)) +
      geom_line(linewidth = 1.1) +
      scale_y_continuous(limits = c(0, 100)) +
      scale_color_manual(values = compare_colors) +
      labs(x = "Day", y = "% Patients with ALP Normalisation",
           title = "Responder Rate — ALP Normalisation (< 120 IU/L)") +
      theme_classic()
    ggplotly(p)
  })
}

## ============================================================
## LAUNCH
## ============================================================
shinyApp(ui = ui, server = server)
