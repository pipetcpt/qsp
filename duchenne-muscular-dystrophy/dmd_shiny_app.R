## ============================================================
## DMD QSP Shiny Application
## Duchenne Muscular Dystrophy — Interactive QSP Dashboard
## 6 Tabs: Patient Profile · PK · Dystrophin & Membrane ·
##         Inflammation & Fibrosis · Clinical Endpoints ·
##         Scenario Comparison & Biomarkers
## ============================================================

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(DT)

## ===== Inline mrgsolve model code =====
CODE_DMD <- '
$PARAM @annotated
KA_CS    : 1.2    : Steroid absorption rate (h-1)
CL_CS    : 25.0   : Steroid clearance (L/h)
V1_CS    : 35.0   : Steroid central Vd (L)
V2_CS    : 80.0   : Steroid peripheral Vd (L)
Q_CS     : 8.0    : Steroid Q (L/h)
CL_ASO   : 210.0  : ASO clearance (mL/h/kg)
V1_ASO   : 290.0  : ASO central V (mL/kg)
V2_ASO   : 6500.0 : ASO tissue V (mL/kg)
Q_ASO    : 18.0   : ASO Q (mL/h/kg)
KUP_ASO  : 0.0015 : ASO muscle uptake (h-1)
KOUT_ASO : 0.0008 : ASO intracellular loss (h-1)
KAAV     : 0.003  : AAV transduction rate (h-1)
DEGRAD   : 5e-5   : AAV vector degradation (h-1)
KEXP_DYS : 0.004  : Micro-dys expression (h-1)
KDEC_DYS : 0.00008: Micro-dys decay (h-1)
DYS_BASE : 0.001  : Baseline dystrophin DMD (%)
EC50_DYS : 5.0    : EC50 dystrophin membrane (%)
HILL_DYS : 2.0    : Hill coefficient
KD_MEM   : 0.15   : Membrane damage rate (h-1)
KR_MEM   : 0.08   : Membrane repair rate (h-1)
KD_CS_MEM: 0.012  : CS membrane stabilization
CA_IN    : 2.5    : Ca influx rate
CA_OUT   : 1.0    : Ca efflux rate
KMEM_CA  : 0.8    : Membrane-Ca coupling
KROS_IN  : 0.5    : ROS generation rate
KROS_CA  : 0.3    : Ca-ROS coefficient
KROS_EL  : 0.4    : ROS elimination
KNF_IN   : 1.2    : NFkB activation rate
KNF_EL   : 0.5    : NFkB deactivation rate
EC50_CS_NF: 15.0  : CS IC50 NFkB (ng/mL)
KM1_IN   : 0.08   : M1 recruitment (h-1)
KM1_EL   : 0.04   : M1 elimination (h-1)
KM2_IN   : 0.03   : M1->M2 switch (h-1)
KM2_EL   : 0.03   : M2 elimination (h-1)
M1_0     : 15.0   : Initial M1 (cells/µL)
KTGF_IN  : 0.2    : TGFb secretion (pg/mL/cell/h)
KTGF_EL  : 0.15   : TGFb clearance (h-1)
KFIB_IN  : 0.002  : Fibrosis progression
KFIB_EL  : 0.00015: Fibrosis resolution
FIB_MAX  : 100.0  : Max fibrosis score
KSC_REGEN: 0.01   : SC self-renewal
KSC_EXHST: 0.002  : SC exhaustion rate
KMF_DEC  : 0.0004 : MF decline rate
KMF_REGEN: 0.0002 : MF recovery rate
KSWMD_DC : 0.00006: 6MWD decline rate
KFVC_DC  : 0.00004: FVC decline rate
KLVEF_DC : 0.000025: LVEF decline rate
EFF_ASO  : 0.04   : ASO dystrophin efficiency
HDAC_FIB : 0.25   : Givinostat fibrosis reduction
WT       : 22.0   : Body weight (kg)

$CMT DEPOT_CS CENT_CS PERIPH_CS CENT_ASO MUS_ASO IC_ASO
     AAV_CIRC AAV_MUS DYS MEMI CAI ROS NFkB M1 M2 TGFb
     FIB SC MF SWD FVC_pct LVEF_pct

$INIT
DEPOT_CS=0 CENT_CS=0 PERIPH_CS=0 CENT_ASO=0 MUS_ASO=0 IC_ASO=0
AAV_CIRC=0 AAV_MUS=0 DYS=0.1 MEMI=0.35 CAI=2.5 ROS=3.0
NFkB=3.5 M1=15.0 M2=8.0 TGFb=12.0 FIB=10.0 SC=100.0
MF=80.0 SWD=380.0 FVC_pct=95.0 LVEF_pct=62.0

$MAIN
double C_CS = CENT_CS / V1_CS * 1000.0;
double DYS_TOTAL = DYS_BASE + DYS + EFF_ASO * IC_ASO;
double DYS_EFF = pow(DYS_TOTAL, HILL_DYS) / (pow(EC50_DYS, HILL_DYS) + pow(DYS_TOTAL, HILL_DYS));
double CS_NF_INH = C_CS / (EC50_CS_NF + C_CS);
double CA_INFLUX = CA_IN * (1.0 - MEMI) * KMEM_CA;
double NECRO_RATE = 0.1 * CAI * ROS;
double SC_EFF = SC / 100.0;
double FIB_EFF_NEG = 1.0 - (FIB / FIB_MAX) * 0.7;

$ODE
dxdt_DEPOT_CS  = -KA_CS * DEPOT_CS;
dxdt_CENT_CS   = KA_CS * DEPOT_CS - (CL_CS/V1_CS)*CENT_CS - (Q_CS/V1_CS)*CENT_CS + (Q_CS/V2_CS)*PERIPH_CS;
dxdt_PERIPH_CS = (Q_CS/V1_CS)*CENT_CS - (Q_CS/V2_CS)*PERIPH_CS;
dxdt_CENT_ASO  = -(CL_ASO/1000.0)*CENT_ASO - (Q_ASO/1000.0)*CENT_ASO + (Q_ASO/1000.0)*MUS_ASO;
dxdt_MUS_ASO   = (Q_ASO/1000.0)*CENT_ASO - (Q_ASO/1000.0)*MUS_ASO - KUP_ASO*MUS_ASO;
dxdt_IC_ASO    = KUP_ASO*MUS_ASO*1000.0 - KOUT_ASO*IC_ASO;
dxdt_AAV_CIRC  = -KAAV*AAV_CIRC - DEGRAD*AAV_CIRC;
dxdt_AAV_MUS   = KAAV*AAV_CIRC - KDEC_DYS*AAV_MUS;
dxdt_DYS       = 0.0;
dxdt_MEMI      = KR_MEM*DYS_EFF*(1.0-MEMI)*SC_EFF - KD_MEM*(1.0-DYS_EFF)*CAI + KD_CS_MEM*C_CS*(1.0-MEMI);
dxdt_CAI       = CA_INFLUX*(1.0-DYS_EFF) + 0.1*(1.0-MEMI) - CA_OUT*CAI;
dxdt_ROS       = KROS_IN*(1.0+KROS_CA*(CAI-1.0)) - KROS_EL*ROS;
double NF_INPUT= 1.0 + 2.0*(ROS-1.0)/3.0 + 0.5*(M1/M1_0-1.0);
dxdt_NFkB      = KNF_IN*(NF_INPUT>1.0?NF_INPUT:1.0) - KNF_EL*NFkB*(1.0+2.0*CS_NF_INH);
double M1_IN   = KM1_IN*NFkB*(1.0-CS_NF_INH*0.7);
dxdt_M1        = M1_IN - KM1_EL*M1 - KM2_IN*M1;
dxdt_M2        = KM2_IN*M1 - KM2_EL*M2;
dxdt_TGFb      = KTGF_IN*M2 - KTGF_EL*TGFb;
dxdt_FIB       = KFIB_IN*TGFb*(1.0-FIB/FIB_MAX)*24.0*(1.0-HDAC_FIB*0.0) - KFIB_EL*FIB;
dxdt_SC        = KSC_REGEN*SC*(1.0-SC/100.0)*DYS_EFF - KSC_EXHST*NECRO_RATE*SC - 0.0001*FIB*SC/100.0;
dxdt_MF        = KMF_REGEN*SC_EFF*DYS_EFF*(100.0-MF) - KMF_DEC*(1.0-DYS_EFF)*(1.0-MEMI)*MF - KMF_DEC*0.5*(FIB/FIB_MAX)*MF;
dxdt_SWD       = -KSWMD_DC*(1.0+(1.0-DYS_EFF)*2.0)*SWD + KSWMD_DC*0.5*(MF/80.0)*DYS_EFF*SWD;
dxdt_FVC_pct   = -KFVC_DC*(1.0-DYS_EFF*0.5)*FVC_pct;
dxdt_LVEF_pct  = -KLVEF_DC*(1.0-DYS_EFF*0.3)*LVEF_pct;

$TABLE
double CK_serum = 20000.0*(1.0-MEMI)*(MF/80.0);
double Dystrophin_pct = DYS_BASE + DYS + EFF_ASO*IC_ASO;
double C_CS_ngmL = CENT_CS / V1_CS * 1000.0;
double NSAA = 34.0*(SWD/400.0)*(MF/100.0);
if(NSAA>34) NSAA=34; if(NSAA<0) NSAA=0;

$CAPTURE C_CS_ngmL CK_serum Dystrophin_pct NSAA
         DYS TGFb FIB SC M1 M2 MEMI CAI ROS NFkB
         SWD FVC_pct LVEF_pct MF IC_ASO
'

## ===== Compile model =====
mod_base <- mcode("dmd_shiny", CODE_DMD)

## ===== Simulation helper =====
run_sim <- function(params, dose_cs = 0, freq_cs = 24, dose_aso = 0,
                    dose_aav = 0, sim_years = 6, age_start = 8,
                    hdac_eff = 0) {
  t_end <- sim_years * 8760
  mod <- param(mod_base, c(params, HDAC_FIB = hdac_eff))

  events <- c()
  if (dose_cs > 0)
    events <- c(events, ev(time=0, amt=dose_cs, cmt="DEPOT_CS", ii=freq_cs, addl=t_end/freq_cs-1))
  if (dose_aso > 0)
    events <- c(events, ev(time=0, amt=dose_aso, cmt="CENT_ASO", ii=168, addl=t_end/168-1))
  if (dose_aav > 0)
    events <- c(events, ev(time=0, amt=dose_aav, cmt="AAV_CIRC"))

  if (length(events) == 0)
    events <- ev(time=0, amt=0, cmt="DEPOT_CS")

  ev_obj <- do.call(rbind, events)
  out <- mrgsim(mod, ev_obj, end = t_end, delta = 24, obsonly = TRUE)
  df <- as.data.frame(out)
  df$age_yr <- age_start + df$time / 8760
  df
}

## ===== UI =====
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "DMD QSP Model"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "tab_patient",    icon = icon("user")),
      menuItem("PK Panel",            tabName = "tab_pk",         icon = icon("chart-line")),
      menuItem("Dystrophin & Membrane",tabName = "tab_dys",       icon = icon("dna")),
      menuItem("Inflammation & Fibrosis", tabName = "tab_inflam", icon = icon("fire")),
      menuItem("Clinical Endpoints",  tabName = "tab_clinical",   icon = icon("heartbeat")),
      menuItem("Scenario Comparison", tabName = "tab_scenarios",  icon = icon("bar-chart"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML(
      ".box-title { font-weight: bold; }
       .scenario-badge { padding: 3px 8px; border-radius: 4px; color: white; margin-right: 4px; }
       .scenario-1 { background-color: #E74C3C; }
       .scenario-2 { background-color: #3498DB; }
       .scenario-3 { background-color: #2ECC71; }
       .scenario-4 { background-color: #9B59B6; }
       .scenario-5 { background-color: #E67E22; }
       .scenario-6 { background-color: #1ABC9C; }"
    ))),

    tabItems(

      ## ============================================================
      ## TAB 1: PATIENT PROFILE
      ## ============================================================
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Characteristics", width = 4, status = "primary", solidHeader = TRUE,
            sliderInput("age_start", "Starting Age (years):", min = 4, max = 12, value = 8, step = 1),
            numericInput("weight",   "Body Weight (kg):", value = 22, min = 10, max = 60),
            selectInput("genotype", "DMD Mutation Type:",
                        choices = c("Exon deletion (amenable to exon 51 skip)" = "del_51",
                                    "Exon deletion (amenable to exon 53 skip)" = "del_53",
                                    "Exon deletion (amenable to exon 45 skip)" = "del_45",
                                    "Nonsense mutation (amenable to ataluren)" = "nonsense",
                                    "Large deletion (not amenable to exon skip)" = "large_del",
                                    "Duplication" = "dup")),
            numericInput("baseline_6mwd",   "Baseline 6MWD (m):", value = 380, min = 100, max = 600),
            numericInput("baseline_fvc",    "Baseline FVC (% pred):", value = 95, min = 40, max = 120),
            numericInput("baseline_lvef",   "Baseline LVEF (%):", value = 62, min = 20, max = 75),
            numericInput("sim_years",       "Simulation Duration (years):", value = 6, min = 1, max = 10)
          ),
          box(title = "Drug Eligibility & Treatment Selection", width = 4, status = "warning", solidHeader = TRUE,
            h4("Corticosteroid"),
            selectInput("cs_drug", "Corticosteroid:",
                        choices = c("None" = "none",
                                    "Deflazacort (0.9 mg/kg/d)" = "dfz",
                                    "Prednisone (0.75 mg/kg/d)" = "pred",
                                    "Vamorolone (6 mg/kg/d)" = "vamo")),
            h4("Exon-Skipping ASO"),
            selectInput("aso_drug", "Exon-Skipping ASO:",
                        choices = c("None" = "none",
                                    "Eteplirsen 30 mg/kg/wk (Exon 51)" = "etepl",
                                    "Golodirsen 30 mg/kg/wk (Exon 53)" = "golod",
                                    "Casimersen 30 mg/kg/wk (Exon 45)" = "casim")),
            h4("Gene / Other Therapy"),
            checkboxInput("use_aav",       "Gene Therapy (Elevidys, single dose)", FALSE),
            checkboxInput("use_givinostat","Givinostat (HDAC inhibitor, 2×5 mg/kg/d)", FALSE),
            checkboxInput("use_acei",      "ACEi + β-Blocker (cardiac protection)", FALSE)
          ),
          box(title = "Disease Overview — Duchenne Muscular Dystrophy", width = 4, status = "info", solidHeader = TRUE,
            HTML('<div style="font-size:13px;">
            <b>Gene:</b> DMD (Xp21.2), 2.4 Mb, 79 exons<br/>
            <b>Protein:</b> Dystrophin (427 kDa, rod-shaped)<br/>
            <b>Prevalence:</b> 1 in 3,500–5,000 live male births<br/>
            <b>Pathophysiology:</b><br/>
            &nbsp;• No dystrophin → DAPC destabilization<br/>
            &nbsp;• Sarcolemmal fragility → Ca²⁺ overload<br/>
            &nbsp;• Necrosis → chronic inflammation<br/>
            &nbsp;• TGF-β → fibrosis → SC exhaustion<br/>
            <b>Clinical Milestones (untreated):</b><br/>
            &nbsp;• Gower sign: ~5yr<br/>
            &nbsp;• LoA: median ~12yr (untreated)<br/>
            &nbsp;• NIV dependence: ~19yr<br/>
            &nbsp;• Cardiomyopathy: ~20yr<br/>
            <b>Approved drugs (2024):</b><br/>
            &nbsp;• Deflazacort (2017), Vamorolone (2023)<br/>
            &nbsp;• Eteplirsen/Golodirsen/Casimersen (2016-21)<br/>
            &nbsp;• Elevidys (2023), Givinostat (2024)<br/>
            </div>')
          )
        ),
        fluidRow(
          valueBoxOutput("vbox_6mwd"),
          valueBoxOutput("vbox_fvc"),
          valueBoxOutput("vbox_lvef"),
          valueBoxOutput("vbox_ck"),
          valueBoxOutput("vbox_fibrosis"),
          valueBoxOutput("vbox_dystrophin")
        )
      ),

      ## ============================================================
      ## TAB 2: PK PANEL
      ## ============================================================
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Corticosteroid Plasma Concentration", width = 6, status = "primary", solidHeader = TRUE,
            plotlyOutput("plot_cs_pk", height = "320px")
          ),
          box(title = "ASO Plasma & Intracellular Concentrations", width = 6, status = "success", solidHeader = TRUE,
            plotlyOutput("plot_aso_pk", height = "320px")
          )
        ),
        fluidRow(
          box(title = "AAV Vector Distribution", width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("plot_aav_pk", height = "320px")
          ),
          box(title = "PK Summary Table", width = 6, status = "info", solidHeader = TRUE,
            DTOutput("dt_pk_summary")
          )
        )
      ),

      ## ============================================================
      ## TAB 3: DYSTROPHIN & MEMBRANE INTEGRITY
      ## ============================================================
      tabItem(tabName = "tab_dys",
        fluidRow(
          box(title = "Dystrophin Level (% of Normal)", width = 6, status = "success", solidHeader = TRUE,
            plotlyOutput("plot_dystrophin", height = "320px")
          ),
          box(title = "Membrane Integrity (fraction, 0=broken, 1=intact)", width = 6, status = "primary", solidHeader = TRUE,
            plotlyOutput("plot_memi", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Intracellular Ca²⁺ (relative to normal=1)", width = 6, status = "danger", solidHeader = TRUE,
            plotlyOutput("plot_cai", height = "300px")
          ),
          box(title = "Reactive Oxygen Species (ROS, relative units)", width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("plot_ros", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Calcium & Membrane Mechanistic Notes", width = 12, status = "info",
            HTML('<b>Pathophysiology of Ca²⁺ Dysregulation in DMD:</b>
            <ul>
            <li>Loss of DAPC → sarcolemmal fragility → microtears during eccentric contraction</li>
            <li>TRPC1/4/6 stretch-activated channels chronically open → Ca²⁺ influx ↑</li>
            <li>SERCA pump inhibited → SR Ca²⁺ reuptake ↓ → cytosolic overload</li>
            <li>Calpain activation → cytoskeletal cleavage → irreversible necrosis</li>
            <li>Mitochondrial Ca²⁺ overload → PTP opening → ROS burst</li>
            <li><b>Dystrophin restoration</b> → DAPC re-assembly → membrane stabilization → Ca²⁺ normalization</li>
            </ul>')
          )
        )
      ),

      ## ============================================================
      ## TAB 4: INFLAMMATION & FIBROSIS
      ## ============================================================
      tabItem(tabName = "tab_inflam",
        fluidRow(
          box(title = "NF-κB Activity (fold over basal)", width = 6, status = "danger", solidHeader = TRUE,
            plotlyOutput("plot_nfkb", height = "300px")
          ),
          box(title = "M1 / M2 Macrophage Dynamics", width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("plot_macro", height = "300px")
          )
        ),
        fluidRow(
          box(title = "TGF-β1 Concentration (pg/mL)", width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("plot_tgfb", height = "300px")
          ),
          box(title = "Muscle Fibrosis Score (0-100)", width = 6, status = "danger", solidHeader = TRUE,
            plotlyOutput("plot_fibrosis", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Satellite Cell Pool (% of Normal)", width = 6, status = "success", solidHeader = TRUE,
            plotlyOutput("plot_sc", height = "300px")
          ),
          box(title = "Inflammation → Fibrosis Pathway Notes", width = 6, status = "info",
            HTML('<b>Inflammation → Fibrosis Cascade:</b>
            <ul>
            <li>DAMP release (HMGB1, ATP, mtDNA) → TLR4/MyD88 → NF-κB</li>
            <li>M1 macrophages: pro-inflammatory (TNF-α, IL-1β, IL-6)</li>
            <li>M2 macrophages: anti-inflammatory but <b>pro-fibrotic</b> (TGF-β1)</li>
            <li>TGF-β1 → Smad2/3 → myofibroblast activation → Collagen I/III ↑</li>
            <li>FAPs (fibro-adipogenic progenitors) → fat infiltration</li>
            <li><b>Givinostat</b> (HDAC inhibitor): redirects FAPs from adipogenic to myogenic fate</li>
            <li>SC exhaustion: Repeated necrosis/regeneration cycles → telomere shortening → senescence</li>
            </ul>')
          )
        )
      ),

      ## ============================================================
      ## TAB 5: CLINICAL ENDPOINTS
      ## ============================================================
      tabItem(tabName = "tab_clinical",
        fluidRow(
          box(title = "6-Minute Walk Distance (6MWD, meters)", width = 6, status = "primary", solidHeader = TRUE,
            plotlyOutput("plot_6mwd", height = "320px")
          ),
          box(title = "North Star Ambulatory Assessment (NSAA, 0-34)", width = 6, status = "success", solidHeader = TRUE,
            plotlyOutput("plot_nsaa", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Respiratory Function — FVC (% predicted)", width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("plot_fvc", height = "300px")
          ),
          box(title = "Cardiac Function — LVEF (%)", width = 6, status = "danger", solidHeader = TRUE,
            plotlyOutput("plot_lvef", height = "300px")
          )
        ),
        fluidRow(
          box(title = "Serum CK Biomarker (U/L)", width = 6, status = "info", solidHeader = TRUE,
            plotlyOutput("plot_ck", height = "300px")
          ),
          box(title = "Clinical Endpoint Summary Table", width = 6, status = "primary",
            DTOutput("dt_clinical")
          )
        )
      ),

      ## ============================================================
      ## TAB 6: SCENARIO COMPARISON & BIOMARKERS
      ## ============================================================
      tabItem(tabName = "tab_scenarios",
        fluidRow(
          box(title = "Run All 6 Scenarios", width = 3, status = "primary", solidHeader = TRUE,
            actionButton("run_all", "Run All Scenarios", icon = icon("play"),
                         class = "btn-primary btn-block"),
            br(),
            h4("Scenario Legend:"),
            HTML('<span class="scenario-badge scenario-1">1</span> Natural History<br/>
                  <span class="scenario-badge scenario-2">2</span> Deflazacort<br/>
                  <span class="scenario-badge scenario-3">3</span> Prednisone<br/>
                  <span class="scenario-badge scenario-4">4</span> Eteplirsen<br/>
                  <span class="scenario-badge scenario-5">5</span> Gene Therapy<br/>
                  <span class="scenario-badge scenario-6">6</span> Dfz + Eteplirsen'),
            br(), br(),
            sliderInput("compare_year", "Year for Comparison Bar Chart:", min = 1, max = 10, value = 6, step = 1)
          ),
          box(title = "6MWD — All Scenarios", width = 9, status = "primary", solidHeader = TRUE,
            plotlyOutput("plot_compare_6mwd", height = "350px")
          )
        ),
        fluidRow(
          box(title = "FVC & Fibrosis Comparison", width = 6, status = "warning", solidHeader = TRUE,
            plotlyOutput("plot_compare_fvc_fib", height = "320px")
          ),
          box(title = "Scenario Comparison at Selected Year (Bar)", width = 6, status = "success", solidHeader = TRUE,
            plotlyOutput("plot_bar_compare", height = "320px")
          )
        ),
        fluidRow(
          box(title = "Comparative Outcomes Table (at selected year)", width = 12, status = "info",
            DTOutput("dt_compare")
          )
        )
      )
    )
  )
)

## ===== SERVER =====
server <- function(input, output, session) {

  ## ---- Reactive: build parameter overrides ----
  patient_params <- reactive({
    list(WT = input$weight)
  })

  ## ---- Reactive: dose inputs ----
  dose_inputs <- reactive({
    cs_dose <- 0
    if (input$cs_drug == "dfz")  cs_dose <- input$weight * 0.9
    if (input$cs_drug == "pred") cs_dose <- input$weight * 0.75
    if (input$cs_drug == "vamo") cs_dose <- input$weight * 6.0

    aso_dose <- 0
    if (input$aso_drug %in% c("etepl","golod","casim")) aso_dose <- 30  # mg/kg equiv.

    aav_dose <- if (input$use_aav) 1.33e14 else 0

    hdac <- if (input$use_givinostat) 0.25 else 0

    list(cs = cs_dose, aso = aso_dose, aav = aav_dose, hdac = hdac)
  })

  ## ---- Main simulation (single run) ----
  sim_data <- reactive({
    doses <- dose_inputs()
    params <- patient_params()
    df <- run_sim(
      params    = params,
      dose_cs   = doses$cs,
      dose_aso  = doses$aso,
      dose_aav  = doses$aav,
      sim_years = input$sim_years,
      age_start = input$age_start,
      hdac_eff  = doses$hdac
    )
    df
  })

  ## ---- All scenarios (run_all button) ----
  all_scenarios <- reactiveVal(NULL)

  observeEvent(input$run_all, {
    withProgress(message = "Running 6 scenarios...", {
      t_end <- input$sim_years * 8760
      w     <- input$weight
      age_s <- input$age_start

      doses_list <- list(
        list(name="1. Natural History",         cs=0,       aso=0,  aav=0,      hdac=0),
        list(name="2. Deflazacort 0.9 mg/kg/d", cs=w*0.9,  aso=0,  aav=0,      hdac=0),
        list(name="3. Prednisone 0.75 mg/kg/d", cs=w*0.75, aso=0,  aav=0,      hdac=0),
        list(name="4. Eteplirsen 30 mg/kg/wk",  cs=0,       aso=30, aav=0,      hdac=0),
        list(name="5. Gene Therapy (Elevidys)",  cs=0,       aso=0,  aav=1.33e14,hdac=0),
        list(name="6. Deflazacort + Eteplirsen", cs=w*0.9,  aso=30, aav=0,      hdac=0)
      )

      all_df <- bind_rows(lapply(doses_list, function(d) {
        df <- run_sim(params=list(WT=w), dose_cs=d$cs, dose_aso=d$aso,
                      dose_aav=d$aav, sim_years=input$sim_years,
                      age_start=age_s, hdac_eff=d$hdac)
        df$scenario <- d$name
        df
      }))

      all_df$scenario <- factor(all_df$scenario, levels = sapply(doses_list, `[[`, "name"))
      all_scenarios(all_df)
    })
  })

  ## ---- Value boxes ----
  last_row <- reactive({ tail(sim_data(), 1) })

  output$vbox_6mwd <- renderValueBox({
    v <- round(last_row()$SWD, 0)
    valueBox(paste0(v, " m"), "6MWD at End", icon = icon("walking"),
             color = ifelse(v > 300, "green", "red"))
  })
  output$vbox_fvc <- renderValueBox({
    v <- round(last_row()$FVC_pct, 1)
    valueBox(paste0(v, " %"), "FVC % Predicted", icon = icon("lungs"),
             color = ifelse(v > 60, "green", ifelse(v > 40, "yellow", "red")))
  })
  output$vbox_lvef <- renderValueBox({
    v <- round(last_row()$LVEF_pct, 1)
    valueBox(paste0(v, " %"), "LVEF", icon = icon("heartbeat"),
             color = ifelse(v > 50, "green", "red"))
  })
  output$vbox_ck <- renderValueBox({
    v <- round(last_row()$CK_serum, 0)
    valueBox(format(v, big.mark=","), "Serum CK (U/L)", icon = icon("vial"),
             color = ifelse(v < 5000, "yellow", "red"))
  })
  output$vbox_fibrosis <- renderValueBox({
    v <- round(last_row()$FIB, 1)
    valueBox(v, "Fibrosis Score (0-100)", icon = icon("layer-group"),
             color = ifelse(v < 40, "green", ifelse(v < 70, "yellow", "red")))
  })
  output$vbox_dystrophin <- renderValueBox({
    v <- round(last_row()$Dystrophin_pct, 2)
    valueBox(paste0(v, " %"), "Dystrophin Level", icon = icon("dna"),
             color = ifelse(v > 4, "green", ifelse(v > 1, "yellow", "red")))
  })

  ## ---- PK plots ----
  output$plot_cs_pk <- renderPlotly({
    df <- sim_data() %>% filter(age_yr <= input$age_start + min(input$sim_years, 0.05))
    if (nrow(df) == 0) df <- sim_data()
    p <- ggplot(sim_data() %>% filter(age_yr <= input$age_start + 0.2),
                aes(x = time/24, y = C_CS_ngmL)) +
      geom_line(color = "#3498DB", linewidth = 1.1) +
      labs(x = "Time (days)", y = "Concentration (ng/mL)",
           title = "Corticosteroid Plasma Conc.") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_aso_pk <- renderPlotly({
    df <- sim_data() %>% filter(age_yr <= input$age_start + 0.1)
    p <- ggplot(df, aes(x = time/24)) +
      geom_line(aes(y = IC_ASO, color = "Intracellular (nmol/L)"), linewidth = 1.1) +
      scale_color_manual(values = c("Intracellular (nmol/L)" = "#9B59B6")) +
      labs(x = "Time (days)", y = "Conc", color = "",
           title = "ASO Concentrations") +
      theme_bw() + theme(legend.position = "top")
    ggplotly(p)
  })

  output$plot_aav_pk <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = age_yr, y = AAV_MUS)) +
      geom_line(color = "#E67E22", linewidth = 1.1) +
      labs(x = "Age (years)", y = "AAV in Muscle (vg/cell scale)",
           title = "AAV Vector in Muscle") +
      theme_bw()
    ggplotly(p)
  })

  output$dt_pk_summary <- renderDT({
    df <- sim_data()
    sum_tbl <- data.frame(
      Parameter     = c("Steroid Cmax (ng/mL)", "Steroid AUC (first day)",
                        "ASO Intracellular (peak, nmol/L)", "AAV Muscle (steady-state)"),
      Value         = c(
        round(max(df$C_CS_ngmL), 1),
        round(sum(df$C_CS_ngmL[df$time <= 24]) * 24, 0),
        round(max(df$IC_ASO), 2),
        round(tail(df$AAV_MUS, 1), 2)
      )
    )
    datatable(sum_tbl, options = list(dom = "t"), rownames = FALSE)
  })

  ## ---- Dystrophin & membrane plots ----
  output$plot_dystrophin <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = age_yr, y = Dystrophin_pct)) +
      geom_line(color = "#27AE60", linewidth = 1.2) +
      geom_hline(yintercept = 4, linetype = "dashed", color = "gray50") +
      annotate("text", x = min(df$age_yr)+0.1, y = 4.5, label = "4% functional threshold",
               hjust = 0, size = 3, color = "gray40") +
      labs(x = "Age (years)", y = "Dystrophin (%)", title = "Dystrophin Level") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_memi <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = age_yr, y = MEMI)) +
      geom_line(color = "#2980B9", linewidth = 1.2) +
      geom_hline(yintercept = 0.35, linetype = "dashed", color = "gray50") +
      ylim(0, 1) +
      labs(x = "Age (years)", y = "Membrane Integrity (0-1)",
           title = "Sarcolemmal Membrane Integrity") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_cai <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = age_yr, y = CAI)) +
      geom_line(color = "#E74C3C", linewidth = 1.2) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      labs(x = "Age (years)", y = "Cytosolic Ca²⁺ (relative)", title = "Intracellular Calcium") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_ros <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = age_yr, y = ROS)) +
      geom_line(color = "#F39C12", linewidth = 1.2) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      labs(x = "Age (years)", y = "ROS (relative)", title = "Reactive Oxygen Species") +
      theme_bw()
    ggplotly(p)
  })

  ## ---- Inflammation & fibrosis plots ----
  output$plot_nfkb <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = age_yr, y = NFkB)) +
      geom_line(color = "#8E44AD", linewidth = 1.2) +
      labs(x = "Age (years)", y = "NF-κB (fold)", title = "NF-κB Inflammatory Activity") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_macro <- renderPlotly({
    df <- sim_data() %>% select(age_yr, M1, M2) %>%
      pivot_longer(cols = c(M1, M2), names_to = "Type", values_to = "Density")
    p <- ggplot(df, aes(x = age_yr, y = Density, color = Type)) +
      geom_line(linewidth = 1.2) +
      scale_color_manual(values = c("M1" = "#E74C3C", "M2" = "#27AE60")) +
      labs(x = "Age (years)", y = "Cells/µL", color = "",
           title = "Macrophage Dynamics (M1/M2)") +
      theme_bw() + theme(legend.position = "top")
    ggplotly(p)
  })

  output$plot_tgfb <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = age_yr, y = TGFb)) +
      geom_line(color = "#E67E22", linewidth = 1.2) +
      labs(x = "Age (years)", y = "TGF-β1 (pg/mL)", title = "TGF-β1 Fibrogenic Signaling") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_fibrosis <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = age_yr, y = FIB)) +
      geom_line(color = "#A04000", linewidth = 1.2) +
      ylim(0, 100) +
      labs(x = "Age (years)", y = "Fibrosis Score (0-100)", title = "Muscle Fibrosis Progression") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_sc <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = age_yr, y = SC)) +
      geom_line(color = "#27AE60", linewidth = 1.2) +
      geom_hline(yintercept = 20, linetype = "dashed", color = "red", alpha = 0.5) +
      ylim(0, 110) +
      labs(x = "Age (years)", y = "SC Pool (% of normal)",
           title = "Satellite Cell Pool Dynamics") +
      theme_bw()
    ggplotly(p)
  })

  ## ---- Clinical endpoint plots ----
  output$plot_6mwd <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = age_yr, y = SWD)) +
      geom_line(color = "#2980B9", linewidth = 1.3) +
      geom_hline(yintercept = 300, linetype = "dashed", color = "red", alpha = 0.6) +
      annotate("text", x = min(df$age_yr)+0.1, y = 310, label = "~LoA threshold",
               hjust = 0, size = 3, color = "red") +
      labs(x = "Age (years)", y = "6MWD (m)", title = "6-Minute Walk Distance") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_nsaa <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = age_yr, y = NSAA)) +
      geom_line(color = "#16A085", linewidth = 1.3) +
      ylim(0, 34) +
      labs(x = "Age (years)", y = "NSAA Score (0-34)",
           title = "North Star Ambulatory Assessment") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_fvc <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = age_yr, y = FVC_pct)) +
      geom_line(color = "#E67E22", linewidth = 1.3) +
      geom_hline(yintercept = 50, linetype = "dashed", color = "red", alpha = 0.6) +
      annotate("text", x = min(df$age_yr)+0.1, y = 51, label = "NIV threshold <50%",
               hjust = 0, size = 3, color = "red") +
      labs(x = "Age (years)", y = "FVC (% predicted)", title = "Respiratory Function (FVC)") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_lvef <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = age_yr, y = LVEF_pct)) +
      geom_line(color = "#E74C3C", linewidth = 1.3) +
      geom_hline(yintercept = 55, linetype = "dashed", color = "gray50") +
      annotate("text", x = min(df$age_yr)+0.1, y = 56, label = "Normal LVEF ≥55%",
               hjust = 0, size = 3, color = "gray40") +
      labs(x = "Age (years)", y = "LVEF (%)", title = "Cardiac Function (LVEF)") +
      theme_bw()
    ggplotly(p)
  })

  output$plot_ck <- renderPlotly({
    df <- sim_data()
    p <- ggplot(df, aes(x = age_yr, y = CK_serum)) +
      geom_line(color = "#9B59B6", linewidth = 1.3) +
      labs(x = "Age (years)", y = "CK (U/L)", title = "Serum CK Biomarker") +
      theme_bw()
    ggplotly(p)
  })

  output$dt_clinical <- renderDT({
    df <- sim_data()
    tpts <- seq(0, max(df$time), by = 8760)
    sum_tbl <- df %>% filter(time %in% tpts) %>%
      mutate(Age = paste0("Age ", input$age_start + round(time/8760), "yr")) %>%
      select(Age, SWD, NSAA, FVC_pct, LVEF_pct, FIB, CK_serum, Dystrophin_pct) %>%
      mutate(across(where(is.numeric), ~round(., 1)))
    datatable(sum_tbl, rownames = FALSE,
              colnames = c("Age", "6MWD(m)", "NSAA", "FVC%", "LVEF%", "Fibrosis", "CK(U/L)", "Dys%"),
              options = list(dom = "t", pageLength = 15))
  })

  ## ---- Scenario comparison plots ----
  scenario_colors <- c(
    "1. Natural History"          = "#E74C3C",
    "2. Deflazacort 0.9 mg/kg/d"  = "#3498DB",
    "3. Prednisone 0.75 mg/kg/d"  = "#2ECC71",
    "4. Eteplirsen 30 mg/kg/wk"   = "#9B59B6",
    "5. Gene Therapy (Elevidys)"  = "#E67E22",
    "6. Deflazacort + Eteplirsen" = "#1ABC9C"
  )

  output$plot_compare_6mwd <- renderPlotly({
    df <- all_scenarios()
    req(df)
    p <- ggplot(df, aes(x = age_yr, y = SWD, color = scenario)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = scenario_colors) +
      geom_hline(yintercept = 300, linetype = "dashed", color = "gray60") +
      labs(x = "Age (years)", y = "6MWD (m)", color = "",
           title = "6MWD — All Treatment Scenarios") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p) %>% layout(legend = list(orientation = "h"))
  })

  output$plot_compare_fvc_fib <- renderPlotly({
    df <- all_scenarios()
    req(df)
    p1 <- ggplot(df, aes(x = age_yr, y = FVC_pct, color = scenario)) +
      geom_line(linewidth = 0.9) +
      scale_color_manual(values = scenario_colors) +
      labs(x = "Age (yr)", y = "FVC%", title = "FVC & Fibrosis") +
      theme_bw() + theme(legend.position = "none")
    ggplotly(p1)
  })

  output$plot_bar_compare <- renderPlotly({
    df <- all_scenarios()
    req(df)
    yr <- input$compare_year
    t_select <- yr * 8760
    df_bar <- df %>%
      filter(abs(time - t_select) == min(abs(time - t_select))) %>%
      group_by(scenario) %>% slice(1) %>% ungroup() %>%
      select(scenario, SWD, FVC_pct, LVEF_pct, Dystrophin_pct) %>%
      pivot_longer(-scenario, names_to = "Metric", values_to = "Value")

    p <- ggplot(df_bar %>% filter(Metric == "SWD"),
                aes(x = scenario, y = Value, fill = scenario)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = scenario_colors) +
      coord_flip() +
      labs(x = "", y = "6MWD (m)",
           title = paste0("6MWD at Year +", yr)) +
      theme_bw() + theme(legend.position = "none")
    ggplotly(p)
  })

  output$dt_compare <- renderDT({
    df <- all_scenarios()
    req(df)
    yr <- input$compare_year
    t_select <- yr * 8760
    df_tbl <- df %>%
      filter(abs(time - t_select) == min(abs(time - t_select))) %>%
      group_by(scenario) %>% slice(1) %>% ungroup() %>%
      select(scenario, SWD, NSAA, FVC_pct, LVEF_pct, FIB, Dystrophin_pct, SC, CK_serum) %>%
      mutate(across(where(is.numeric), ~round(., 1)))
    datatable(df_tbl, rownames = FALSE,
              colnames = c("Scenario","6MWD(m)","NSAA","FVC%","LVEF%","Fibrosis","Dys%","SC%","CK(U/L)"),
              options = list(dom = "t", pageLength = 10))
  })
}

## ===== Launch =====
shinyApp(ui = ui, server = server)
