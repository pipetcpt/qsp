## ============================================================
## Myotonic Dystrophy Type 1 (DM1) — Shiny QSP Dashboard
## 근긴장성 이영양증 제1형 인터랙티브 시뮬레이션 대시보드
##
## 7 Tabs:
##  1. Patient Profile & Genetics
##  2. Drug PK (Mexiletine / ASO)
##  3. Muscle & Myotonia PD
##  4. Cardiac Safety & Endpoints
##  5. Multi-Scenario Comparison
##  6. Biomarker Panel (Splicing Indices)
##  7. CNS & Metabolic Outcomes
## ============================================================

# install.packages(c("shiny","shinydashboard","dplyr","ggplot2",
#                    "tidyr","mrgsolve","DT","scales","plotly"))
library(shiny)
library(shinydashboard)
library(dplyr)
library(ggplot2)
library(tidyr)
library(mrgsolve)
library(DT)
library(scales)
library(plotly)

## ============================================================
## Inline mrgsolve model code (condensed)
## ============================================================
dm1_code <- '
$PARAM
CTG_repeat=400 BW=70
Ka_mex=1.2 F_mex=0.85 Vc_mex=65 Vp_mex=195 CL_mex=32 Q_mex=12
MW_mex=179.26 fu_mex=0.30
k_aso_abs=0.30 V_aso_p=5 k_aso_dist=0.05 k_aso_nuc=0.15
k_aso_elim=0.003 k_aso_pelim=0.10
KD_MBNL1=0.25 k_MBNL1_eq=0.05 MBNL1_total=1.0
CUGBP1_max=2.5 EC50_CUG_CUGBP=0.5 hill_CUG_CUGBP=1.5
CLCN1_fetal_norm=0.05 CLCN1_fetal_DM1=0.80
SERCA1_fetal_norm=0.10 SERCA1_fetal_DM1=0.60
INSR_fetal_norm=0.30 INSR_fetal_DM1=0.75
k_splice_eq=0.008 k_foci_form=0.002
Myotonia_max=7.0 k_myot_eq=0.50
EC50_ClC1_myo=0.40 hill_myo=2.0
IC50_mex_nav=3.0 hill_mex_nav=1.8
Grip_DM1_base=18.0 k_grip_eq=0.001
Muscle_mass_0=24.0 k_muscle_loss=0.0002
PR_baseline=210 QTc_baseline=440 k_PR_prog=0.0001 k_QTc_eq=0.05
HOMA_IR_0=3.5 k_HOMA_eq=0.01 FVC_0=75.0 k_FVC_prog=0.0001
MEX_ON=0 ASO_ON=0 GT_ON=0

$CMT MEX_GUT MEX_CENT MEX_PERI ASO_PLASMA ASO_MUSCLE ASO_NUCL
     CUG_FOCI MBNL1_FREE CUGBP1_ACT CLCN1_FETAL SERCA_FETAL INSR_FETAL
     MYOTONIA GRIP_STR MUSCLE_MASS PR_INT QTc_INT HOMA_IR FVC_PCT

$MAIN
double CTG_ratio = CTG_repeat/400.0;
double CTG_factor = pow(CTG_ratio, 0.6);
double ASO_eff = ASO_ON * ASO_NUCL/(ASO_NUCL+0.5);
double foci_target = 1.0 * CTG_factor * (1.0-ASO_eff);
if(GT_ON>0) foci_target *= 0.15;
double MBNL1_free_target = MBNL1_total/(1.0+CUG_FOCI/KD_MBNL1);
if(GT_ON>0) MBNL1_free_target = 0.85;
double CUGBP1_target = 1.0 + (CUGBP1_max-1.0)*pow(CUG_FOCI,hill_CUG_CUGBP)/
                       (pow(EC50_CUG_CUGBP,hill_CUG_CUGBP)+pow(CUG_FOCI,hill_CUG_CUGBP));
double mbnl1_eff = MBNL1_FREE;
double cugbp1_eff = CUGBP1_ACT/CUGBP1_max;
double CLCN1_tgt = CLCN1_fetal_norm + (CLCN1_fetal_DM1-CLCN1_fetal_norm)*(1-mbnl1_eff)*cugbp1_eff*CTG_factor;
if(CLCN1_tgt>0.95) CLCN1_tgt=0.95;
if(CLCN1_tgt<CLCN1_fetal_norm) CLCN1_tgt=CLCN1_fetal_norm;
double SERCA_tgt = SERCA1_fetal_norm + (SERCA1_fetal_DM1-SERCA1_fetal_norm)*(1-mbnl1_eff)*cugbp1_eff*CTG_factor;
if(SERCA_tgt>0.90) SERCA_tgt=0.90;
double INSR_tgt = INSR_fetal_norm + (INSR_fetal_DM1-INSR_fetal_norm)*(1-mbnl1_eff)*cugbp1_eff*CTG_factor;
if(INSR_tgt>0.90) INSR_tgt=0.90;
double ClC1_act = 1.0 - CLCN1_FETAL;
double Cp_mex = (MEX_CENT/Vc_mex)*1000.0/MW_mex;
double Cp_free = Cp_mex*fu_mex;
double nav_block = MEX_ON*pow(Cp_free,hill_mex_nav)/(pow(IC50_mex_nav,hill_mex_nav)+pow(Cp_free,hill_mex_nav));
double ClC1_myo_red = pow(ClC1_act,hill_myo)/(pow(EC50_ClC1_myo,hill_myo)+pow(ClC1_act,hill_myo));
double myot_tgt = Myotonia_max*(1-ClC1_myo_red)*(1-0.80*nav_block);
if(myot_tgt<0) myot_tgt=0;
double grip_tgt = Grip_DM1_base*(1-0.5*SERCA_FETAL)*(1+0.2*ClC1_act);
double PR_tgt = PR_baseline + 20*CTG_factor*(1-0.4*ASO_eff);
double QTc_mex_eff = 15*nav_block;
double QTc_tgt = QTc_baseline + 15*CTG_factor - QTc_mex_eff;
if(QTc_tgt<380) QTc_tgt=380;
double HOMA_tgt = HOMA_IR_0*(1+1.5*(INSR_FETAL-INSR_fetal_norm));
double FVC_tgt = FVC_0*(MUSCLE_MASS/Muscle_mass_0);
if(FVC_tgt<20) FVC_tgt=20;

$ODE
dxdt_MEX_GUT  = -Ka_mex*MEX_GUT;
dxdt_MEX_CENT = Ka_mex*MEX_GUT-(CL_mex/Vc_mex)*MEX_CENT-(Q_mex/Vc_mex)*MEX_CENT+(Q_mex/Vp_mex)*MEX_PERI;
dxdt_MEX_PERI = (Q_mex/Vc_mex)*MEX_CENT-(Q_mex/Vp_mex)*MEX_PERI;
dxdt_ASO_PLASMA = -k_aso_pelim*ASO_PLASMA - k_aso_dist*ASO_PLASMA;
dxdt_ASO_MUSCLE = k_aso_dist*ASO_PLASMA - k_aso_nuc*ASO_MUSCLE - k_aso_elim*ASO_MUSCLE;
dxdt_ASO_NUCL   = k_aso_nuc*ASO_MUSCLE - k_aso_elim*ASO_NUCL;
dxdt_CUG_FOCI   = k_foci_form*(foci_target-CUG_FOCI);
dxdt_MBNL1_FREE = k_MBNL1_eq*(MBNL1_free_target-MBNL1_FREE);
dxdt_CUGBP1_ACT = 0.10*(CUGBP1_target-CUGBP1_ACT);
dxdt_CLCN1_FETAL = k_splice_eq*(CLCN1_tgt-CLCN1_FETAL);
dxdt_SERCA_FETAL = k_splice_eq*(SERCA_tgt-SERCA_FETAL);
dxdt_INSR_FETAL  = k_splice_eq*(INSR_tgt-INSR_FETAL);
dxdt_MYOTONIA   = k_myot_eq*(myot_tgt-MYOTONIA);
dxdt_GRIP_STR   = k_grip_eq*(grip_tgt-GRIP_STR);
dxdt_MUSCLE_MASS = -k_muscle_loss*MUSCLE_MASS + 0.00005*ClC1_act*MUSCLE_MASS;
dxdt_PR_INT = k_PR_prog*(PR_tgt-PR_INT);
dxdt_QTc_INT = k_QTc_eq*(QTc_tgt-QTc_INT);
dxdt_HOMA_IR = k_HOMA_eq*(HOMA_tgt-HOMA_IR);
dxdt_FVC_PCT = 0.001*(FVC_tgt-FVC_PCT);

$INIT
MEX_GUT=0 MEX_CENT=0 MEX_PERI=0 ASO_PLASMA=0 ASO_MUSCLE=0 ASO_NUCL=0
CUG_FOCI=1.0 MBNL1_FREE=0.20 CUGBP1_ACT=2.0
CLCN1_FETAL=0.78 SERCA_FETAL=0.55 INSR_FETAL=0.72
MYOTONIA=5.5 GRIP_STR=18.0 MUSCLE_MASS=22.0
PR_INT=215 QTc_INT=445 HOMA_IR=3.8 FVC_PCT=73.0

$CAPTURE Cp_mex nav_block ClC1_act ASO_eff
'

mod_dm1 <- mcode("DM1_Shiny", dm1_code)

## ============================================================
## Helper: run simulation
## ============================================================
run_sim <- function(ctg, mex_dose, mex_on, aso_dose, aso_on, gt_on,
                    sim_days = 365, delta = 6) {
  idata <- data.frame(CTG_repeat = ctg, MEX_ON = mex_on,
                      ASO_ON = aso_on, GT_ON = gt_on)

  ev_list <- list()
  if (mex_on == 1 && mex_dose > 0) {
    times <- seq(0, sim_days * 24, by = 8)
    ev_list <- c(ev_list, list(ev(cmt = "MEX_GUT", amt = mex_dose, time = times)))
  }
  if (aso_on == 1 && aso_dose > 0) {
    times <- seq(0, sim_days * 24, by = 28 * 24)
    ev_list <- c(ev_list, list(ev(cmt = "ASO_PLASMA", amt = aso_dose, time = times)))
  }
  events <- if (length(ev_list) > 0) do.call(c, ev_list) else ev(amt = 0, time = 0, cmt = 1)

  out <- mrgsim(mod_dm1, idata = idata, events = events,
                end = sim_days * 24, delta = delta,
                carry_out = c("CTG_repeat"))
  df <- as.data.frame(out)
  df$Day <- df$time / 24
  df
}

theme_dm1 <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "#E3F2FD"),
        plot.title = element_text(face = "bold", size = 13))

## ============================================================
## UI
## ============================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "DM1 QSP Dashboard", titleWidth = 300),

  dashboardSidebar(
    width = 270,
    sidebarMenu(
      menuItem("Patient Profile", tabName = "patient", icon = icon("user-circle")),
      menuItem("Drug PK", tabName = "pk", icon = icon("pills")),
      menuItem("Muscle & Myotonia", tabName = "muscle", icon = icon("hand-paper")),
      menuItem("Cardiac Safety", tabName = "cardiac", icon = icon("heartbeat")),
      menuItem("Scenario Comparison", tabName = "scenario", icon = icon("chart-bar")),
      menuItem("Biomarker Panel", tabName = "biomarker", icon = icon("flask")),
      menuItem("CNS & Metabolic", tabName = "cns", icon = icon("brain"))
    ),
    hr(),
    h5("Patient Parameters", style="padding-left:15px; color:#90CAF9;"),
    sliderInput("ctg_repeat", "CTG Repeat Size",
                min = 100, max = 2000, value = 400, step = 50),
    hr(),
    h5("Mexiletine", style="padding-left:15px; color:#90CAF9;"),
    checkboxInput("mex_on", "Mexiletine ON", value = FALSE),
    sliderInput("mex_dose", "Dose (mg TID)",
                min = 50, max = 400, value = 300, step = 50),
    hr(),
    h5("ASO Therapy", style="padding-left:15px; color:#90CAF9;"),
    checkboxInput("aso_on", "ASO ON (DYNE-101 / AOC 1001)", value = FALSE),
    sliderInput("aso_dose", "ASO Dose (mg / 4 wks)",
                min = 50, max = 400, value = 200, step = 50),
    hr(),
    h5("Gene Therapy", style="padding-left:15px; color:#90CAF9;"),
    checkboxInput("gt_on", "AAV-MBNL1 (Experimental)", value = FALSE),
    hr(),
    sliderInput("sim_days", "Simulation (days)",
                min = 7, max = 730, value = 365, step = 7),
    actionButton("run", "▶ Run Simulation", class = "btn-primary btn-block")
  ),

  dashboardBody(
    tags$head(tags$style(HTML(
      ".content-wrapper { background-color: #F8F9FA; }
       .box { border-top: 3px solid #1565C0; }
       .info-box { min-height: 80px; }
       .info-box-icon { height: 80px; line-height: 80px; }
       .info-box-content { padding-top: 10px; }"
    ))),

    tabItems(

      ## ---- TAB 1: Patient Profile ----
      tabItem("patient",
        fluidRow(
          infoBoxOutput("box_severity", width = 4),
          infoBoxOutput("box_mbnl1",   width = 4),
          infoBoxOutput("box_foci",    width = 4)
        ),
        fluidRow(
          box(title = "DM1 Disease Mechanism", width = 6, solidHeader = TRUE,
              status = "primary",
              HTML('<img src="" alt="" style="width:100%"/>
              <table class="table table-bordered" style="font-size:11px;">
              <thead><tr><th>Component</th><th>Normal</th><th>DM1</th><th>Consequence</th></tr></thead>
              <tbody>
              <tr><td>CTG Repeat</td><td>5–37</td><td>>50 → &gt;10,000</td><td>RNA gain-of-function</td></tr>
              <tr><td>MBNL1 (free)</td><td>~95%</td><td>10–30%</td><td>Fetal splicing reversion</td></tr>
              <tr><td>CUGBP1/CELF1</td><td>Normal</td><td>2–3× hyperactive</td><td>Competing splicing</td></tr>
              <tr><td>CLCN1 fetal %</td><td>~5%</td><td>60–90%</td><td>Myotonia</td></tr>
              <tr><td>INSR fetal %</td><td>~30%</td><td>70–80%</td><td>Insulin resistance</td></tr>
              <tr><td>SERCA1 fetal %</td><td>~10%</td><td>50–70%</td><td>Ca²⁺ dysregulation</td></tr>
              </tbody></table>'
              )
          ),
          box(title = "Severity Staging (CTG-based)", width = 6, solidHeader = TRUE,
              status = "warning",
              tableOutput("severity_table"),
              hr(),
              HTML('<div style="font-size:11px; color:#555;">
              <b>MIRS Grading:</b> 1=None · 2=Mild (CTG <400) · 3=Moderate (400–800) ·
              4=Severe (>800) · 5=Requiring assistance<br>
              <b>Risk categories:</b> Cardiac (annual ECG/Holter) · Respiratory (spirometry) ·
              Cognitive (INI-Brief) · Metabolic (HbA1c/HOMA-IR)
              </div>')
          )
        ),
        fluidRow(
          box(title = "Initial Conditions at Simulation Start",
              width = 12, solidHeader = TRUE, status = "info",
              tableOutput("init_conditions"))
        )
      ),

      ## ---- TAB 2: Drug PK ----
      tabItem("pk",
        fluidRow(
          box(title = "Mexiletine Plasma Concentration (μM)",
              width = 7, solidHeader = TRUE, status = "primary",
              plotlyOutput("plot_mex_pk", height = 350)),
          box(title = "PK Summary", width = 5, solidHeader = TRUE, status = "info",
              tableOutput("pk_summary"),
              hr(),
              HTML('<div style="font-size:11px;">
              <b>Key PK parameters (Warner 2015):</b><br>
              Ka = 1.2 h⁻¹ · Vc = 65 L · CL = 32 L/h<br>
              t½ = 10–12 h · F = 85% · fu = 30%<br>
              <b>IC₅₀ Nav1.4</b> = 3 μM (muscle)<br>
              <b>IC₅₀ Nav1.5</b> = 12 μM (cardiac — monitoring needed)<br>
              <b>CYP2D6 PM</b>: dose 50% reduction recommended
              </div>'))
        ),
        fluidRow(
          box(title = "ASO Tissue Concentration (Muscle)",
              width = 7, solidHeader = TRUE, status = "success",
              plotlyOutput("plot_aso_pk", height = 350)),
          box(title = "ASO PK Properties", width = 5, solidHeader = TRUE, status = "success",
              HTML('<div style="font-size:12px;">
              <b>DYNE-101 / AOC 1001 Pharmacokinetics:</b><br><br>
              ● Subcutaneous (SC) or IV delivery<br>
              ● Tissue half-life: ~2–4 weeks (muscle)<br>
              ● Nuclear uptake: 15–20% of muscle conc.<br>
              ● Dosing: every 4 weeks (monthly)<br>
              ● DMPK mRNA knockdown: 50–80% at wk 12<br>
              ● CLCN1 splicing rescue: ~40 pp improvement<br><br>
              <b>Mechanism:</b> RNase H-mediated degradation of<br>
              mutant DMPK mRNA → dissolution of nuclear<br>
              CUGn RNA foci → release of sequestered MBNL1
              </div>'))
        )
      ),

      ## ---- TAB 3: Muscle & Myotonia ----
      tabItem("muscle",
        fluidRow(
          valueBoxOutput("vbox_myotonia", width = 4),
          valueBoxOutput("vbox_grip",     width = 4),
          valueBoxOutput("vbox_nav",      width = 4)
        ),
        fluidRow(
          box(title = "Myotonia VAS Score Over Time",
              width = 6, solidHeader = TRUE, status = "primary",
              plotlyOutput("plot_myotonia", height = 320)),
          box(title = "Hand Grip Strength (kg)",
              width = 6, solidHeader = TRUE, status = "success",
              plotlyOutput("plot_grip", height = 320))
        ),
        fluidRow(
          box(title = "Nav1.4 Blockade & ClC-1 Activity",
              width = 6, solidHeader = TRUE, status = "warning",
              plotlyOutput("plot_nav_clc1", height = 300)),
          box(title = "Muscle Mass Over Time",
              width = 6, solidHeader = TRUE, status = "danger",
              plotlyOutput("plot_muscle", height = 300))
        )
      ),

      ## ---- TAB 4: Cardiac Safety ----
      tabItem("cardiac",
        fluidRow(
          valueBoxOutput("vbox_PR", width = 3),
          valueBoxOutput("vbox_QTc", width = 3),
          valueBoxOutput("vbox_cardiac_risk", width = 3),
          valueBoxOutput("vbox_mex_cardiac", width = 3)
        ),
        fluidRow(
          box(title = "QTc Interval Over Time",
              width = 6, solidHeader = TRUE, status = "danger",
              plotlyOutput("plot_qtc", height = 320)),
          box(title = "PR Interval Over Time",
              width = 6, solidHeader = TRUE, status = "warning",
              plotlyOutput("plot_pr", height = 320))
        ),
        fluidRow(
          box(title = "Cardiac Risk Stratification (DM1)", width = 6,
              solidHeader = TRUE, status = "primary",
              HTML('<table class="table table-sm" style="font-size:11px;">
              <thead><tr><th>Parameter</th><th>Low Risk</th><th>Intermediate</th><th>High Risk</th></tr></thead>
              <tbody>
              <tr><td>PR interval</td><td>&lt;200 ms</td><td>200–240 ms</td><td>&gt;240 ms</td></tr>
              <tr><td>QRS duration</td><td>&lt;120 ms</td><td>120–140 ms</td><td>&gt;140 ms</td></tr>
              <tr><td>HV interval</td><td>&lt;55 ms</td><td>55–70 ms</td><td>&gt;70 ms</td></tr>
              <tr><td>QTc</td><td>&lt;440 ms</td><td>440–460 ms</td><td>&gt;460 ms</td></tr>
              </tbody></table>
              <small>Groh et al. 2008 NEJM: HV > 70 ms → 5× SCD risk.
              Annual ECG + 24h Holter recommended for all DM1 patients.</small>')),
          box(title = "Mexiletine Cardiac Safety Notes", width = 6,
              solidHeader = TRUE, status = "warning",
              HTML('<div style="font-size:12px;">
              <b>Class Ib Antiarrhythmic Properties:</b><br>
              ● Nav1.5 blockade: IC₅₀ ~12 μM (less than Nav1.4 ~3 μM)<br>
              ● QTc shortening: 10–20 ms at therapeutic doses<br>
              ● Use-dependent block: faster at high HR<br><br>
              <b>Contraindications / Cautions:</b><br>
              ● Existing significant conduction disease (PR >300 ms)<br>
              ● 2nd/3rd degree AV block without pacemaker<br>
              ● Sick sinus syndrome<br><br>
              <b>Monitoring (MELT Trial Protocol):</b><br>
              ● ECG at baseline, 2 wk, 4 wk, each dose change<br>
              ● HV interval via EPS if significant QRS widening
              </div>'))
        )
      ),

      ## ---- TAB 5: Scenario Comparison ----
      tabItem("scenario",
        fluidRow(
          box(title = "Define Comparison Scenarios", width = 12,
              solidHeader = TRUE, status = "primary",
              column(3, checkboxGroupInput("sc_list",
                label = "Scenarios to compare:",
                choices = list(
                  "Natural History" = 1,
                  "Mexiletine 200 mg TID" = 2,
                  "Mexiletine 300 mg TID" = 3,
                  "ASO monthly (DYNE-101)" = 4,
                  "Mexiletine + ASO (combo)" = 5,
                  "Gene Therapy (AAV-MBNL1)" = 6,
                  "Severe DM1 (CTG=1200)" = 7
                ),
                selected = c(1, 3, 4, 5, 6)
              )),
              column(9,
                plotlyOutput("plot_scenario_myotonia", height = 250),
                plotlyOutput("plot_scenario_clcn1", height = 250)
              )
          )
        ),
        fluidRow(
          box(title = "1-Year Outcomes Summary Table",
              width = 12, solidHeader = TRUE, status = "success",
              DTOutput("scenario_table"))
        )
      ),

      ## ---- TAB 6: Biomarker Panel ----
      tabItem("biomarker",
        fluidRow(
          box(title = "CLCN1 Fetal Splicing Index (%)",
              width = 6, solidHeader = TRUE, status = "primary",
              plotlyOutput("plot_clcn1", height = 300)),
          box(title = "INSR Fetal Splicing Index (%)",
              width = 6, solidHeader = TRUE, status = "warning",
              plotlyOutput("plot_insr", height = 300))
        ),
        fluidRow(
          box(title = "SERCA1 Fetal Splicing Index (%)",
              width = 6, solidHeader = TRUE, status = "danger",
              plotlyOutput("plot_serca", height = 300)),
          box(title = "MBNL1 Free Fraction & CUG Foci",
              width = 6, solidHeader = TRUE, status = "success",
              plotlyOutput("plot_mbnl1", height = 300))
        ),
        fluidRow(
          box(title = "Splicing Biomarker Reference Values", width = 12,
              solidHeader = TRUE, status = "info",
              HTML('<table class="table table-bordered" style="font-size:11px;">
              <thead><tr><th>Biomarker</th><th>Normal Adult</th><th>DM1 (moderate)</th>
              <th>After ASO</th><th>Method</th></tr></thead>
              <tbody>
              <tr><td>CLCN1 exon 7a %</td><td>~5%</td><td>60–90%</td><td>30–50%</td>
                  <td>RT-PCR / muscle biopsy</td></tr>
              <tr><td>INSR exon 11 excl. %</td><td>~30%</td><td>65–80%</td><td>45–60%</td>
                  <td>RT-PCR / blood</td></tr>
              <tr><td>SERCA1 exon 22 excl. %</td><td>~10%</td><td>40–70%</td><td>25–45%</td>
                  <td>RT-PCR / muscle biopsy</td></tr>
              <tr><td>RNA foci count</td><td>0 foci/nucleus</td><td>10–50 foci</td><td>2–10 foci</td>
                  <td>FISH / muscle biopsy</td></tr>
              <tr><td>DMPK mRNA</td><td>1.0 (relative)</td><td>1.0 (retained)</td><td>0.2–0.5</td>
                  <td>RT-qPCR / blood</td></tr>
              </tbody></table>'))
        )
      ),

      ## ---- TAB 7: CNS & Metabolic ----
      tabItem("cns",
        fluidRow(
          box(title = "HOMA-IR (Insulin Resistance Index)",
              width = 6, solidHeader = TRUE, status = "warning",
              plotlyOutput("plot_homa", height = 300)),
          box(title = "FVC % Predicted (Pulmonary Function)",
              width = 6, solidHeader = TRUE, status = "danger",
              plotlyOutput("plot_fvc", height = 300))
        ),
        fluidRow(
          box(title = "CNS Manifestations in DM1", width = 6,
              solidHeader = TRUE, status = "primary",
              HTML('<table class="table table-sm" style="font-size:11px;">
              <thead><tr><th>Manifestation</th><th>Prevalence</th><th>Mechanism</th></tr></thead>
              <tbody>
              <tr><td>Excessive daytime sleepiness</td><td>70–80%</td><td>Orexin ↓, MBNL2 CNS sequestration</td></tr>
              <tr><td>Cognitive impairment</td><td>60–70%</td><td>Tau 4R/3R imbalance, white matter lesions</td></tr>
              <tr><td>Depression / anxiety</td><td>60–70%</td><td>GABAa splicing, frontal atrophy</td></tr>
              <tr><td>Executive dysfunction</td><td>50–60%</td><td>Prefrontal white matter</td></tr>
              <tr><td>Personality changes</td><td>40–50%</td><td>Frontal-temporal involvement</td></tr>
              <tr><td>Central sleep apnea</td><td>20–30%</td><td>CNS respiratory control</td></tr>
              </tbody></table>')),
          box(title = "Metabolic & Systemic Features", width = 6,
              solidHeader = TRUE, status = "warning",
              HTML('<table class="table table-sm" style="font-size:11px;">
              <thead><tr><th>Feature</th><th>Prevalence</th><th>Mechanism</th></tr></thead>
              <tbody>
              <tr><td>Insulin resistance</td><td>75%</td><td>INSR-A fetal isoform, GLUT4 ↓</td></tr>
              <tr><td>Type 2 Diabetes</td><td>20–30%</td><td>Chronic IR progression</td></tr>
              <tr><td>Testicular atrophy</td><td>50–80% males</td><td>DMPK kinase haploinsufficiency</td></tr>
              <tr><td>Dysphagia</td><td>60–70%</td><td>Pharyngeal smooth muscle involvement</td></tr>
              <tr><td>GI dysmotility</td><td>50–60%</td><td>Smooth muscle splicing defects</td></tr>
              <tr><td>Posterior cataracts</td><td>50–80%</td><td>Lens MBNL1/2 loss</td></tr>
              <tr><td>Hypothyroidism</td><td>10–20%</td><td>Thyroid smooth muscle / direct</td></tr>
              </tbody></table>'))
        )
      )
    )
  )
)

## ============================================================
## SERVER
## ============================================================
server <- function(input, output, session) {

  ## ---- Reactive simulation ----
  sim_data <- eventReactive(input$run, {
    run_sim(
      ctg      = input$ctg_repeat,
      mex_dose = input$mex_dose,
      mex_on   = as.integer(input$mex_on),
      aso_dose = input$aso_dose,
      aso_on   = as.integer(input$aso_on),
      gt_on    = as.integer(input$gt_on),
      sim_days = input$sim_days
    )
  }, ignoreNULL = FALSE)

  ## Multi-scenario data
  scenario_defs <- list(
    list(label="Natural History", ctg=400, mex_on=0, mex_dose=0, aso_on=0, aso_dose=0, gt_on=0, color="#555555"),
    list(label="Mexiletine 200mg", ctg=400, mex_on=1, mex_dose=200, aso_on=0, aso_dose=0, gt_on=0, color="#42A5F5"),
    list(label="Mexiletine 300mg", ctg=400, mex_on=1, mex_dose=300, aso_on=0, aso_dose=0, gt_on=0, color="#1565C0"),
    list(label="ASO monthly", ctg=400, mex_on=0, mex_dose=0, aso_on=1, aso_dose=200, gt_on=0, color="#FF6F00"),
    list(label="Mex+ASO combo", ctg=400, mex_on=1, mex_dose=300, aso_on=1, aso_dose=200, gt_on=0, color="#6A1B9A"),
    list(label="Gene Therapy", ctg=400, mex_on=0, mex_dose=0, aso_on=0, aso_dose=0, gt_on=1, color="#1B5E20"),
    list(label="Severe (CTG=1200)", ctg=1200, mex_on=0, mex_dose=0, aso_on=0, aso_dose=0, gt_on=0, color="#C62828")
  )

  sc_data_all <- reactive({
    req(input$sc_list)
    sel <- as.integer(input$sc_list)
    bind_rows(lapply(sel, function(i) {
      sc <- scenario_defs[[i]]
      df <- run_sim(sc$ctg, sc$mex_dose, sc$mex_on, sc$aso_dose, sc$aso_on, sc$gt_on, 365)
      df$Label <- sc$label
      df$Color <- sc$color
      df
    }))
  })

  ## ---- InfoBoxes: Patient Profile ----
  output$box_severity <- renderInfoBox({
    ctg <- input$ctg_repeat
    sev <- if(ctg < 200) "Mild" else if(ctg < 400) "Moderate" else if(ctg < 800) "Severe" else "Very Severe"
    infoBox("CTG Repeat / Severity", paste0(ctg, " repeats — ", sev),
            icon = icon("dna"), color = if(ctg>800) "red" else if(ctg>400) "orange" else "yellow")
  })
  output$box_mbnl1 <- renderInfoBox({
    ctg <- input$ctg_repeat
    mbnl1_pct <- round(100 / (1 + (ctg/400)^0.6 / 0.25))
    infoBox("Est. Free MBNL1 %", paste0(mbnl1_pct, "% (normal ~95%)"),
            icon = icon("microscope"), color = if(mbnl1_pct < 20) "red" else "orange")
  })
  output$box_foci <- renderInfoBox({
    ctg <- input$ctg_repeat
    foci_rel <- round((ctg/400)^0.6, 2)
    infoBox("RNA Foci Burden", paste0(foci_rel, "× (relative to CTG=400)"),
            icon = icon("atom"), color = if(foci_rel > 2) "red" else "yellow")
  })

  ## Severity table
  output$severity_table <- renderTable({
    data.frame(
      CTG_Range = c("50–150","150–400","400–800","800–2000",">2000"),
      Severity = c("Mild DM2-like","Moderate (classic)","Severe","Very severe","Congenital-DM1 overlap"),
      Myotonia_VAS = c("1–2","3–5","5–7","7–9","8–10"),
      Cardiac_Risk = c("Low","Moderate","High","Very High","Extreme")
    )
  }, striped = TRUE, bordered = TRUE, small = TRUE)

  ## Initial conditions
  output$init_conditions <- renderTable({
    ctg <- input$ctg_repeat
    data.frame(
      Parameter = c("CLCN1 fetal %","INSR fetal %","SERCA1 fetal %",
                    "MBNL1 free %","Myotonia VAS","Grip (kg)",
                    "PR interval (ms)","QTc (ms)","HOMA-IR","FVC %"),
      Value = c("~78%","~72%","~55%","~20%","5.5","18.0","215","445","3.8","73%"),
      Normal = c("~5%","~30%","~10%","~95%","0","36","<200","<440","<2.5",">80%")
    )
  }, striped = TRUE, bordered = TRUE, small = TRUE)

  ## ---- PK Plots ----
  output$plot_mex_pk <- renderPlotly({
    df <- sim_data()
    p <- df %>% filter(Day <= min(10, max(df$Day))) %>%
      ggplot(aes(Day, Cp_mex)) +
      geom_line(color = "#1565C0", linewidth = 1.2) +
      geom_hline(yintercept = 3, lty = 2, color = "red") +
      geom_hline(yintercept = 10, lty = 2, color = "darkred") +
      annotate("label", x = 8, y = 3.5, label = "IC50 Nav1.4 (3μM)", size = 2.5) +
      annotate("label", x = 8, y = 10.5, label = "Upper target (10μM)", size = 2.5) +
      labs(x = "Day", y = "Cp Mexiletine (μM)") +
      theme_dm1
    ggplotly(p)
  })

  output$pk_summary <- renderTable({
    df <- sim_data()
    ss <- df %>% filter(Day >= 3, Day <= 4)
    data.frame(
      Metric = c("Cmax (μM)", "Cmin (μM)", "Css,avg (μM)", "Nav1.4 Block (avg%)"),
      Value = c(round(max(ss$Cp_mex), 2),
                round(min(ss$Cp_mex), 2),
                round(mean(ss$Cp_mex), 2),
                round(mean(ss$nav_block) * 100, 1))
    )
  })

  output$plot_aso_pk <- renderPlotly({
    df <- sim_data()
    p <- df %>%
      ggplot(aes(Day, ASO_MUSCLE)) +
      geom_line(color = "#00695C", linewidth = 1.2) +
      geom_line(aes(y = ASO_NUCL), color = "#004D40", lty = 2, linewidth = 1) +
      annotate("label", x = max(df$Day)*0.85, y = max(df$ASO_MUSCLE)*0.8,
               label = "Muscle tissue", size = 2.5, color = "#00695C") +
      annotate("label", x = max(df$Day)*0.85, y = max(df$ASO_NUCL)*0.8,
               label = "Nuclear (active)", size = 2.5, color = "#004D40") +
      labs(x = "Day", y = "ASO Concentration (μg/g equivalent)") +
      theme_dm1
    ggplotly(p)
  })

  ## ---- Muscle & Myotonia Plots ----
  output$vbox_myotonia <- renderValueBox({
    df <- sim_data()
    val <- round(tail(df$MYOTONIA, 1), 1)
    valueBox(val, "Final Myotonia VAS", icon = icon("hand-paper"),
             color = if(val > 5) "red" else if(val > 2) "orange" else "green")
  })
  output$vbox_grip <- renderValueBox({
    df <- sim_data()
    val <- round(tail(df$GRIP_STR, 1), 1)
    valueBox(paste0(val, " kg"), "Final Grip Strength", icon = icon("dumbbell"),
             color = if(val < 12) "red" else if(val < 20) "orange" else "green")
  })
  output$vbox_nav <- renderValueBox({
    df <- sim_data()
    val <- round(mean(df$nav_block[df$Day >= (max(df$Day)-2)]) * 100, 0)
    valueBox(paste0(val, "%"), "Nav1.4 Blockade (SS)", icon = icon("bolt"),
             color = if(val > 70) "green" else if(val > 30) "yellow" else "red")
  })

  output$plot_myotonia <- renderPlotly({
    df <- sim_data()
    p <- df %>% ggplot(aes(Day, MYOTONIA)) +
      geom_line(color = "#C62828", linewidth = 1.2) +
      geom_hline(yintercept = 0, lty = 2, color = "green4") +
      scale_y_continuous(limits = c(0, 8)) +
      labs(x = "Day", y = "Myotonia VAS (0–10)") + theme_dm1
    ggplotly(p)
  })

  output$plot_grip <- renderPlotly({
    df <- sim_data()
    p <- df %>% ggplot(aes(Day, GRIP_STR)) +
      geom_line(color = "#2E7D32", linewidth = 1.2) +
      geom_hline(yintercept = 36, lty = 2, color = "blue", alpha = 0.7) +
      annotate("label", x = max(df$Day)*0.8, y = 37, label = "Normal ~36 kg", size = 2.5) +
      labs(x = "Day", y = "Grip Strength (kg)") + theme_dm1
    ggplotly(p)
  })

  output$plot_nav_clc1 <- renderPlotly({
    df <- sim_data()
    p <- df %>% select(Day, nav_block, ClC1_act) %>%
      pivot_longer(-Day, names_to = "Variable", values_to = "Value") %>%
      mutate(Variable = recode(Variable,
                               nav_block = "Nav1.4 Blockade",
                               ClC1_act  = "ClC-1 Activity")) %>%
      ggplot(aes(Day, Value, color = Variable)) +
      geom_line(linewidth = 1.1) +
      scale_color_manual(values = c("Nav1.4 Blockade" = "#1565C0", "ClC-1 Activity" = "#2E7D32")) +
      scale_y_continuous(labels = percent_format()) +
      labs(x = "Day", y = "Fractional Activity") + theme_dm1
    ggplotly(p)
  })

  output$plot_muscle <- renderPlotly({
    df <- sim_data()
    p <- df %>% ggplot(aes(Day, MUSCLE_MASS)) +
      geom_line(color = "#6A1B9A", linewidth = 1.2) +
      labs(x = "Day", y = "Lean Muscle Mass (kg)") + theme_dm1
    ggplotly(p)
  })

  ## ---- Cardiac Plots ----
  output$vbox_PR <- renderValueBox({
    df <- sim_data()
    val <- round(tail(df$PR_INT, 1), 0)
    valueBox(paste0(val, " ms"), "PR Interval", icon = icon("heartbeat"),
             color = if(val > 240) "red" else if(val > 200) "orange" else "green")
  })
  output$vbox_QTc <- renderValueBox({
    df <- sim_data()
    val <- round(tail(df$QTc_INT, 1), 0)
    valueBox(paste0(val, " ms"), "QTc Interval", icon = icon("wave-square"),
             color = if(val > 460) "red" else if(val > 440) "orange" else "green")
  })
  output$vbox_cardiac_risk <- renderValueBox({
    ctg <- input$ctg_repeat
    risk <- if(ctg > 1000) "Very High" else if(ctg > 400) "High" else "Moderate"
    valueBox(risk, "SCD Risk Level", icon = icon("exclamation-triangle"),
             color = if(risk == "Very High") "red" else if(risk == "High") "orange" else "yellow")
  })
  output$vbox_mex_cardiac <- renderValueBox({
    df <- sim_data()
    nav <- mean(df$nav_block[df$Day >= (max(df$Day)-2)])
    qtc_short <- round(15 * nav, 0)
    valueBox(paste0("-", qtc_short, " ms"), "QTc Shortening (Mexiletine)",
             icon = icon("shield-alt"), color = "blue")
  })

  output$plot_qtc <- renderPlotly({
    df <- sim_data()
    p <- df %>% ggplot(aes(Day, QTc_INT)) +
      geom_line(color = "#C62828", linewidth = 1.2) +
      geom_hline(yintercept = 450, lty = 2, color = "red") +
      geom_hline(yintercept = 440, lty = 3, color = "orange") +
      annotate("label", x = max(df$Day)*0.7, y = 452, label = "Alert >450 ms", size = 2.5) +
      labs(x = "Day", y = "QTc (ms)") + theme_dm1
    ggplotly(p)
  })

  output$plot_pr <- renderPlotly({
    df <- sim_data()
    p <- df %>% ggplot(aes(Day, PR_INT)) +
      geom_line(color = "#FF6F00", linewidth = 1.2) +
      geom_hline(yintercept = 200, lty = 2, color = "orange") +
      geom_hline(yintercept = 240, lty = 2, color = "red") +
      annotate("label", x = max(df$Day)*0.7, y = 242, label = "High risk >240 ms", size = 2.5) +
      labs(x = "Day", y = "PR Interval (ms)") + theme_dm1
    ggplotly(p)
  })

  ## ---- Scenario Comparison ----
  output$plot_scenario_myotonia <- renderPlotly({
    df <- sc_data_all()
    colors <- setNames(sapply(as.integer(input$sc_list), function(i) scenario_defs[[i]]$color),
                       sapply(as.integer(input$sc_list), function(i) scenario_defs[[i]]$label))
    p <- df %>% ggplot(aes(Day, MYOTONIA, color = Label)) +
      geom_line(linewidth = 1.0) +
      scale_color_manual(values = colors) +
      labs(x = "Day", y = "Myotonia VAS") + theme_dm1 +
      theme(legend.position = "right", legend.text = element_text(size=8))
    ggplotly(p)
  })

  output$plot_scenario_clcn1 <- renderPlotly({
    df <- sc_data_all()
    colors <- setNames(sapply(as.integer(input$sc_list), function(i) scenario_defs[[i]]$color),
                       sapply(as.integer(input$sc_list), function(i) scenario_defs[[i]]$label))
    p <- df %>% ggplot(aes(Day, CLCN1_FETAL * 100, color = Label)) +
      geom_line(linewidth = 1.0) +
      geom_hline(yintercept = 5, lty = 2, color = "green4", alpha = 0.7) +
      scale_color_manual(values = colors) +
      labs(x = "Day", y = "CLCN1 Fetal % (splicing index)") + theme_dm1 +
      theme(legend.position = "right", legend.text = element_text(size=8))
    ggplotly(p)
  })

  output$scenario_table <- renderDT({
    df <- sc_data_all()
    tbl <- df %>% filter(Day >= 364) %>%
      group_by(Label) %>%
      summarise(
        `Myotonia VAS` = round(mean(MYOTONIA), 2),
        `Grip (kg)` = round(mean(GRIP_STR), 1),
        `CLCN1 fetal %` = round(mean(CLCN1_FETAL*100), 1),
        `INSR fetal %` = round(mean(INSR_FETAL*100), 1),
        `PR (ms)` = round(mean(PR_INT), 0),
        `QTc (ms)` = round(mean(QTc_INT), 0),
        `HOMA-IR` = round(mean(HOMA_IR), 2),
        `FVC %` = round(mean(FVC_PCT), 1),
        .groups = "drop"
      )
    datatable(tbl, options = list(pageLength = 10, dom = "t"), rownames = FALSE) %>%
      formatStyle("Myotonia VAS", backgroundColor = styleInterval(c(2, 5), c("#C8E6C9","#FFF9C4","#FFCDD2")))
  })

  ## ---- Biomarker Plots ----
  output$plot_clcn1 <- renderPlotly({
    df <- sim_data()
    p <- df %>% ggplot(aes(Day, CLCN1_FETAL*100)) +
      geom_line(color = "#1565C0", linewidth=1.2) +
      geom_hline(yintercept=5, lty=2, color="green4") +
      labs(x="Day", y="CLCN1 Fetal Isoform (%)") + theme_dm1
    ggplotly(p)
  })
  output$plot_insr <- renderPlotly({
    df <- sim_data()
    p <- df %>% ggplot(aes(Day, INSR_FETAL*100)) +
      geom_line(color="#F9A825", linewidth=1.2) +
      geom_hline(yintercept=30, lty=2, color="green4") +
      labs(x="Day", y="INSR Fetal Isoform (%)") + theme_dm1
    ggplotly(p)
  })
  output$plot_serca <- renderPlotly({
    df <- sim_data()
    p <- df %>% ggplot(aes(Day, SERCA_FETAL*100)) +
      geom_line(color="#C62828", linewidth=1.2) +
      geom_hline(yintercept=10, lty=2, color="green4") +
      labs(x="Day", y="SERCA1 Fetal Isoform (%)") + theme_dm1
    ggplotly(p)
  })
  output$plot_mbnl1 <- renderPlotly({
    df <- sim_data()
    p <- df %>%
      select(Day, MBNL1_FREE, CUG_FOCI) %>%
      pivot_longer(-Day) %>%
      mutate(name = recode(name, MBNL1_FREE="MBNL1 Free Fraction", CUG_FOCI="CUG Foci Burden")) %>%
      ggplot(aes(Day, value, color=name)) +
      geom_line(linewidth=1.1) +
      scale_color_manual(values=c("MBNL1 Free Fraction"="#2E7D32","CUG Foci Burden"="#C62828")) +
      labs(x="Day", y="Relative Level (0–1)") + theme_dm1 +
      theme(legend.position="bottom")
    ggplotly(p)
  })

  ## ---- CNS & Metabolic Plots ----
  output$plot_homa <- renderPlotly({
    df <- sim_data()
    p <- df %>% ggplot(aes(Day, HOMA_IR)) +
      geom_line(color="#E65100", linewidth=1.2) +
      geom_hline(yintercept=2.5, lty=2, color="green4") +
      annotate("label", x=max(df$Day)*0.7, y=2.7, label="Normal <2.5", size=2.5) +
      labs(x="Day", y="HOMA-IR") + theme_dm1
    ggplotly(p)
  })
  output$plot_fvc <- renderPlotly({
    df <- sim_data()
    p <- df %>% ggplot(aes(Day, FVC_PCT)) +
      geom_line(color="#880E4F", linewidth=1.2) +
      geom_hline(yintercept=50, lty=2, color="red") +
      geom_hline(yintercept=70, lty=2, color="orange") +
      annotate("label", x=max(df$Day)*0.7, y=52, label="Ventilator threshold <50%", size=2.5) +
      labs(x="Day", y="FVC % Predicted") + theme_dm1
    ggplotly(p)
  })
}

## ============================================================
## RUN
## ============================================================
shinyApp(ui, server)
