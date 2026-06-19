## ============================================================
## Addison's Disease (Primary Adrenal Insufficiency)
## Interactive QSP Dashboard — Shiny App
##
## Tabs:
##   1. Patient Profile & Disease Overview
##   2. Pharmacokinetics (HC, FC, DHEA)
##   3. HPA Axis & Cortisol Dynamics
##   4. Electrolytes & Hemodynamics
##   5. Clinical Endpoints & Biomarkers
##   6. Treatment Scenario Comparison
##   7. Adrenal Crisis Risk
##   8. Dose Optimization
## ============================================================

library(shiny)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(shinydashboard)
library(DT)

# ──────────────────────────────────────────────────────────────
# MODEL CODE (embedded)
# ──────────────────────────────────────────────────────────────
code <- '
$PARAM
Ka_HC=1.20, F_HC=0.96, CL_HC=90.0, Vc_HC=15.0, Vp_HC=32.0, Q_HC=25.0, fu_HC=0.05,
Ka_FC=1.50, CL_FC=4.0, Vc_FC=80.0, Vp_FC=30.0, Q_FC=5.0,
Ka_DHEA=0.80, CL_DHEA=12.0, Vc_DHEA=50.0,
k_CRH_syn=0.50, k_CRH_deg=0.60, CRH0=0.833,
k_ACTH_syn=2.40, k_ACTH_deg=0.30, ACTH0=8.0,
ACTH_EC50=0.4, ACTH_Emax=3.0,
GC_IC50=15.0, GC_hill=1.5,
k_Cort_prod=0.80, k_Cort_deg=0.40, ACTH_Cort_EC50=5.0,
Circ_amp=0.60, Circ_phase=6.0, Circ_period=24.0,
GR_tot=1.0, kon_GR=0.50, koff_GR=0.10, k_GR_nuc=0.20, k_GR_ret=0.05,
k_mRNA_syn=0.30, k_mRNA_deg=0.20,
k_Na_in=0.40, k_Na_out=0.40, Na_SS=140.0,
FC_Na_Emax=0.60, FC_Na_EC50=0.008,
k_K_in=0.30, k_K_out=0.30, K_SS=4.0,
FC_K_Emax=0.50, FC_K_EC50=0.008,
MAP_SS=90.0, k_MAP_vol=0.08, MAP_Na_coeff=0.20, k_MAP_ret=0.10,
Gluc0=90.0, k_Gluc_GC=0.05, k_Gluc_base=0.03, Gluc_ins_sens=0.04,
BMD0=1.0, k_BMD_deg_GC=0.002, k_BMD_regen=0.001,
AR0=0.05, k_AR_prog=0.0002,
k_ACTH_sig=0.10, k_MSH_deg=0.15,
k_crisis=0.0010, cortisol_safe=5.0, k_crisis_decay=0.005,
DHEA_FLAG=0

$CMT HC_GUT HC_CENT HC_PERI FC_GUT FC_CENT DHEA_CENT
     CRH ACTH_CMT Cort_endo
     GR_FREE GR_BOUND GR_MRNA
     Na_CMT K_CMT MAP_CMT
     Gluc_CMT BMD_CMT ACTH_SIG AR_CMT Crisis_risk

$INIT
HC_GUT=0, HC_CENT=0, HC_PERI=0,
FC_GUT=0, FC_CENT=0, DHEA_CENT=0,
CRH=0.833, ACTH_CMT=8.0, Cort_endo=0.0,
GR_FREE=0.90, GR_BOUND=0.05, GR_MRNA=1.50,
Na_CMT=135.0, K_CMT=5.2, MAP_CMT=75.0,
Gluc_CMT=70.0, BMD_CMT=0.95, ACTH_SIG=1.5,
AR_CMT=0.05, Crisis_risk=0.0

$ODE
double HC_Cp = HC_CENT / Vc_HC;
double HC_free = HC_Cp * fu_HC;
double HC_nmol = HC_free * 2.76;
double FC_Cp = FC_CENT / Vc_FC;

dxdt_HC_GUT  = -Ka_HC * HC_GUT;
dxdt_HC_CENT =  Ka_HC * F_HC * HC_GUT - (CL_HC + Q_HC)/Vc_HC*HC_CENT + Q_HC/Vp_HC*HC_PERI;
dxdt_HC_PERI =  Q_HC/Vc_HC*HC_CENT - Q_HC/Vp_HC*HC_PERI;
dxdt_FC_GUT  = -Ka_FC * FC_GUT;
dxdt_FC_CENT =  Ka_FC * FC_GUT - (CL_FC+Q_FC)/Vc_FC*FC_CENT + Q_FC/Vp_FC*FC_PERI;
dxdt_DHEA_CENT = Ka_DHEA * DHEA_FLAG * 1000 - CL_DHEA/Vc_DHEA*DHEA_CENT;

double t_mod = fmod(SOLVERTIME, Circ_period);
double circ  = 1.0 + Circ_amp * cos(2.0*3.14159*(t_mod - Circ_phase)/Circ_period);
double Cort_total = Cort_endo + HC_Cp/2.76;
double GC_free_nmol = HC_nmol + Cort_endo * 27.6;
double GC_fb = 1.0 / (1.0 + pow(GC_free_nmol/GC_IC50, GC_hill));
double CRH_norm = CRH/CRH0;
double ACTH_stim = ACTH_Emax * CRH_norm / (ACTH_EC50 + CRH_norm);

dxdt_CRH      = k_CRH_syn * circ * GC_fb - k_CRH_deg * CRH;
dxdt_ACTH_CMT = k_ACTH_syn * ACTH_stim * GC_fb - k_ACTH_deg * ACTH_CMT;
double AR = AR_CMT < 0.001 ? 0.001 : AR_CMT;
dxdt_Cort_endo = k_Cort_prod * AR * ACTH_CMT/(ACTH_Cort_EC50+ACTH_CMT) - k_Cort_deg*Cort_endo;

double GR_bind = kon_GR * HC_nmol * GR_FREE;
double GR_dis  = koff_GR * GR_BOUND;
dxdt_GR_FREE  = GR_tot * k_GR_ret - GR_bind + GR_dis;
dxdt_GR_BOUND = GR_bind - GR_dis - k_GR_nuc * GR_BOUND;
double GR_nuc_eq = GR_BOUND * k_GR_nuc / (k_GR_ret + k_mRNA_deg);
dxdt_GR_MRNA  = k_mRNA_syn*(1.0+GR_nuc_eq) - k_mRNA_deg*GR_MRNA;

double FC_eff_Na = FC_Na_Emax * FC_Cp/(FC_Na_EC50+FC_Cp);
dxdt_Na_CMT = k_Na_in*(1.0+FC_eff_Na) - k_Na_out*Na_CMT/Na_SS;
double FC_eff_K  = FC_K_Emax  * FC_Cp/(FC_K_EC50+FC_Cp);
dxdt_K_CMT  = k_K_in - k_K_out*(1.0+FC_eff_K)*K_CMT/K_SS;

double Na_dev = (Na_CMT-Na_SS)/Na_SS;
dxdt_MAP_CMT = k_MAP_ret*(MAP_SS + (MAP_Na_coeff*Na_dev + k_MAP_vol*FC_Cp/0.01)*MAP_SS - MAP_CMT);

dxdt_Gluc_CMT = k_Gluc_GC*Cort_total + Gluc0*k_Gluc_base - k_Gluc_base*(1.0+Gluc_ins_sens)*Gluc_CMT/Gluc0;
double GC_excess = Cort_total > 20.0 ? (Cort_total-20.0)/20.0 : 0.0;
dxdt_BMD_CMT = k_BMD_regen*(1.0-BMD_CMT) - k_BMD_deg_GC*GC_excess*BMD_CMT;
dxdt_ACTH_SIG = k_ACTH_sig*ACTH_CMT - k_MSH_deg*ACTH_SIG;

double immune_act = 1.0 - GC_free_nmol/(GC_free_nmol+50.0);
dxdt_AR_CMT = -k_AR_prog * immune_act * AR_CMT;

double crisis_d = (Cort_total < cortisol_safe) ? k_crisis*(cortisol_safe-Cort_total)/cortisol_safe : 0.0;
dxdt_Crisis_risk = crisis_d - k_crisis_decay*Crisis_risk;

$TABLE
double Total_Cort = Cort_endo + HC_CENT/Vc_HC/2.76;
double Cp_HC_ugdL = HC_CENT/Vc_HC/27.6;
double Cp_FC_ngmL = FC_CENT/Vc_FC*1000.0;
double Cp_DHEA    = DHEA_CENT/Vc_DHEA;

$CAPTURE Total_Cort Cp_HC_ugdL Cp_FC_ngmL Cp_DHEA
         ACTH_CMT Cort_endo
         Na_CMT K_CMT MAP_CMT
         Gluc_CMT BMD_CMT ACTH_SIG AR_CMT Crisis_risk
         GR_FREE GR_BOUND GR_MRNA
'

mod <- mcode("pai_shiny", code, quiet = TRUE)

# ──────────────────────────────────────────────────────────────
# SIMULATION HELPER
# ──────────────────────────────────────────────────────────────
run_sim <- function(hc_am, hc_nn, hc_pm, fc_mcg, use_dhea,
                    hc_mode = "IR", days = 90, stress = FALSE,
                    stress_start = 30, stress_days = 7, stress_mult = 2) {

  end_t <- days * 24
  Ka_use <- ifelse(hc_mode == "MR", 0.30, 1.20)

  ev_HC_am <- ev(amt = hc_am * 1000, cmt = "HC_GUT", time = 8,  ii = 24, addl = days - 1)
  ev_HC_nn <- ev(amt = hc_nn * 1000, cmt = "HC_GUT", time = 12, ii = 24, addl = days - 1)
  ev_HC_pm <- ev(amt = hc_pm * 1000, cmt = "HC_GUT", time = 18, ii = 24, addl = days - 1)
  ev_FC    <- ev(amt = fc_mcg,        cmt = "FC_GUT", time = 8,  ii = 24, addl = days - 1)

  evts <- c(ev_HC_am, ev_HC_nn, ev_HC_pm, ev_FC)
  if (use_dhea) {
    ev_DHEA <- ev(amt = 25000, cmt = "DHEA_CENT", time = 8, ii = 24, addl = days - 1)
    evts <- c(evts, ev_DHEA)
  }

  if (stress) {
    t0 <- stress_start * 24
    ev_st_am <- ev(amt = hc_am * stress_mult * 1000, cmt = "HC_GUT",
                   time = 8 + t0,  ii = 24, addl = stress_days - 1)
    ev_st_nn <- ev(amt = hc_nn * stress_mult * 1000, cmt = "HC_GUT",
                   time = 12 + t0, ii = 24, addl = stress_days - 1)
    ev_st_pm <- ev(amt = hc_pm * stress_mult * 1000, cmt = "HC_GUT",
                   time = 18 + t0, ii = 24, addl = stress_days - 1)
    evts <- c(evts, ev_st_am, ev_st_nn, ev_st_pm)
  }

  mod %>%
    param(Ka_HC = Ka_use,
          DHEA_FLAG = ifelse(use_dhea, 1, 0)) %>%
    ev(evts) %>%
    mrgsim(end = end_t, delta = 0.5) %>%
    as.data.frame() %>%
    mutate(Day = time / 24, TOD = time %% 24)
}

# ──────────────────────────────────────────────────────────────
# UI
# ──────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Addison's Disease QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("1. Patient Profile",         tabName = "tab_profile",   icon = icon("user-md")),
      menuItem("2. Pharmacokinetics",         tabName = "tab_pk",        icon = icon("pills")),
      menuItem("3. HPA Axis & Cortisol",      tabName = "tab_hpa",       icon = icon("brain")),
      menuItem("4. Electrolytes & BP",        tabName = "tab_elec",      icon = icon("heartbeat")),
      menuItem("5. Clinical Endpoints",       tabName = "tab_clinical",  icon = icon("chart-line")),
      menuItem("6. Scenario Comparison",      tabName = "tab_scenario",  icon = icon("balance-scale")),
      menuItem("7. Adrenal Crisis Risk",      tabName = "tab_crisis",    icon = icon("exclamation-triangle")),
      menuItem("8. Dose Optimization",        tabName = "tab_dose",      icon = icon("sliders-h"))
    ),

    hr(),
    h4("  Dosing Parameters", style = "color:white;"),
    sliderInput("hc_am",  "HC Morning (mg)", 5, 20, 10, 1),
    sliderInput("hc_nn",  "HC Noon (mg)",    0, 10,  5, 1),
    sliderInput("hc_pm",  "HC Evening (mg)", 0, 10,  5, 1),
    sliderInput("fc_mcg", "FC Dose (μg)",    25, 300, 100, 25),
    radioButtons("hc_mode", "HC Formulation",
                 choices = c("Immediate-Release" = "IR",
                             "Modified-Release"  = "MR"),
                 selected = "IR"),
    checkboxInput("use_dhea", "Add DHEA 25 mg/day", FALSE),
    sliderInput("sim_days", "Simulation Duration (days)", 30, 365, 90, 15),
    hr(),
    checkboxInput("stress_on", "Simulate Illness Episode", FALSE),
    sliderInput("stress_start", "Illness Start (day)", 5, 60, 20, 1),
    sliderInput("stress_days",  "Illness Duration (days)", 1, 14, 7, 1),
    sliderInput("stress_mult",  "Stress Dose Multiplier", 1, 4, 2, 0.5),
    hr(),
    actionButton("run_sim", "Run Simulation", icon = icon("play"),
                 style = "color:white;background-color:#1976D2;width:90%")
  ),

  dashboardBody(
    tabItems(

      # ── Tab 1: Patient Profile ─────────────────────────────
      tabItem("tab_profile",
        fluidRow(
          box(title = "Disease Overview — Primary Adrenal Insufficiency", width = 12, status = "primary",
            HTML('<p><b>Addison\'s Disease</b> (Primary Adrenal Insufficiency, PAI) is characterized by
                  destruction of the adrenal cortex, leading to deficiency of <b>cortisol</b> (glucocorticoid),
                  <b>aldosterone</b> (mineralocorticoid), and <b>DHEA</b> (androgen).</p>
                  <ul>
                    <li><b>Prevalence:</b> ~100–140 per million (UK/Scandinavia)</li>
                    <li><b>Autoimmune:</b> 80–90% of cases (anti-21-hydroxylase antibodies)</li>
                    <li><b>Peak incidence:</b> 3rd–4th decade; F:M ≈ 2:1</li>
                    <li><b>APS-2:</b> 50% have associated autoimmune thyroid disease</li>
                  </ul>
                  <p><b>Key pathophysiology:</b> CD4+ Th1/Th17 and CD8+ CTL attack adrenocortical cells →
                  cortisol/aldosterone/DHEA deficiency → HPA axis disinhibition → ACTH ↑↑ (→ hyperpigmentation).</p>')
          )
        ),
        fluidRow(
          valueBoxOutput("vbox_cortisol", width = 3),
          valueBoxOutput("vbox_acth",     width = 3),
          valueBoxOutput("vbox_na",       width = 3),
          valueBoxOutput("vbox_k",        width = 3)
        ),
        fluidRow(
          box(title = "Clinical Presentation", width = 6, status = "warning",
            tableOutput("dx_table")
          ),
          box(title = "Diagnostic Criteria (Endocrine Society 2016)", width = 6, status = "info",
            tableOutput("dx_criteria")
          )
        ),
        fluidRow(
          box(title = "QSP Model Structure", width = 12,
            HTML('<img src="https://via.placeholder.com/1000x300?text=QSP+Model+Mechanistic+Map" width="100%">
                  <p>Model: 20 ODE compartments covering HC/FC/DHEA PK, HPA axis (CRH-ACTH-cortisol),
                  GR/MR receptor signaling, Na+/K+/MAP, BMD, blood glucose, adrenal reserve, crisis risk.</p>')
          )
        )
      ),

      # ── Tab 2: Pharmacokinetics ────────────────────────────
      tabItem("tab_pk",
        fluidRow(
          box(title = "Hydrocortisone Plasma Concentration", width = 6, status = "primary",
            plotlyOutput("pk_hc_plot", height = 350)
          ),
          box(title = "Fludrocortisone Plasma Concentration", width = 6, status = "info",
            plotlyOutput("pk_fc_plot", height = 350)
          )
        ),
        fluidRow(
          box(title = "PK Parameters Summary", width = 12,
            DTOutput("pk_table")
          )
        ),
        fluidRow(
          box(title = "Cortisol AUC₀₋₂₄ Analysis", width = 6, status = "success",
            plotlyOutput("pk_auc_plot", height = 300)
          ),
          box(title = "HC Formulation Comparison (IR vs MR)", width = 6, status = "warning",
            plotlyOutput("pk_form_plot", height = 300)
          )
        )
      ),

      # ── Tab 3: HPA Axis ────────────────────────────────────
      tabItem("tab_hpa",
        fluidRow(
          box(title = "Plasma ACTH Trajectory", width = 6, status = "danger",
            plotlyOutput("hpa_acth_plot", height = 350)
          ),
          box(title = "Circadian Cortisol Profile (Steady State)", width = 6, status = "primary",
            plotlyOutput("hpa_circ_plot", height = 350)
          )
        ),
        fluidRow(
          box(title = "CRH Dynamics", width = 4, status = "info",
            plotlyOutput("hpa_crh_plot", height = 280)
          ),
          box(title = "GR Occupancy", width = 4, status = "success",
            plotlyOutput("hpa_gr_plot", height = 280)
          ),
          box(title = "Hyperpigmentation Index (ACTH Signal)", width = 4, status = "warning",
            plotlyOutput("hpa_pig_plot", height = 280)
          )
        )
      ),

      # ── Tab 4: Electrolytes & Blood Pressure ──────────────
      tabItem("tab_elec",
        fluidRow(
          box(title = "Serum Sodium (Na⁺)", width = 4, status = "primary",
            plotlyOutput("elec_na_plot", height = 300)
          ),
          box(title = "Serum Potassium (K⁺)", width = 4, status = "warning",
            plotlyOutput("elec_k_plot", height = 300)
          ),
          box(title = "Mean Arterial Pressure", width = 4, status = "danger",
            plotlyOutput("elec_map_plot", height = 300)
          )
        ),
        fluidRow(
          box(title = "Electrolyte Reference Ranges", width = 12,
            HTML('<table style="width:100%; border-collapse:collapse;">
              <tr style="background:#1976D2;color:white;">
                <th>Parameter</th><th>Normal Range</th><th>Addison\'s (untreated)</th><th>Target (treated)</th>
              </tr>
              <tr><td>Na⁺ (mEq/L)</td><td>135–145</td><td>125–134 (hyponatremia)</td><td>138–142</td></tr>
              <tr><td>K⁺ (mEq/L)</td><td>3.5–5.0</td><td>5.5–7.0 (hyperkalemia)</td><td>3.8–4.5</td></tr>
              <tr><td>MAP (mmHg)</td><td>70–100</td><td>55–70 (hypotension)</td><td>80–95</td></tr>
              <tr><td>Glucose (mg/dL)</td><td>70–100</td><td>50–70 (hypoglycemia)</td><td>80–110</td></tr>
            </table>')
          )
        )
      ),

      # ── Tab 5: Clinical Endpoints ─────────────────────────
      tabItem("tab_clinical",
        fluidRow(
          box(title = "Blood Glucose", width = 6, status = "warning",
            plotlyOutput("clin_gluc_plot", height = 300)
          ),
          box(title = "Bone Mineral Density (normalized)", width = 6, status = "success",
            plotlyOutput("clin_bmd_plot", height = 300)
          )
        ),
        fluidRow(
          box(title = "Functional Adrenal Reserve", width = 6, status = "info",
            plotlyOutput("clin_ar_plot", height = 280)
          ),
          box(title = "DHEA-S Levels", width = 6, status = "primary",
            plotlyOutput("clin_dhea_plot", height = 280)
          )
        ),
        fluidRow(
          box(title = "Clinical Outcomes Summary at End of Simulation", width = 12,
            DTOutput("outcomes_table")
          )
        )
      ),

      # ── Tab 6: Scenario Comparison ────────────────────────
      tabItem("tab_scenario",
        fluidRow(
          box(title = "Treatment Scenario Definitions", width = 12, status = "primary",
            HTML('<table style="width:100%;border-collapse:collapse;">
              <tr style="background:#1976D2;color:white;">
                <th>#</th><th>Scenario</th><th>HC Dose</th><th>FC Dose</th><th>DHEA</th><th>Notes</th>
              </tr>
              <tr><td>1</td><td>No Treatment</td><td>—</td><td>—</td><td>—</td><td>Natural disease course</td></tr>
              <tr><td>2</td><td>Standard HC IR + FC</td><td>10+5+5 mg/day</td><td>100 μg</td><td>—</td><td>Standard of care</td></tr>
              <tr><td>3</td><td>Modified-Release HC + FC</td><td>20 mg/day (MR)</td><td>100 μg</td><td>—</td><td>Plenadren/Efmody</td></tr>
              <tr><td>4</td><td>Triple Replacement</td><td>10+5+5 mg/day</td><td>100 μg</td><td>25 mg</td><td>DHEA added (Arlt 2006)</td></tr>
              <tr><td>5</td><td>Stress Dosing</td><td>2× during illness</td><td>100 μg</td><td>—</td><td>Sick-day rules</td></tr>
            </table>')
          )
        ),
        fluidRow(
          box(title = "Cortisol Comparison Across Scenarios", width = 6, status = "primary",
            plotlyOutput("scen_cort_plot", height = 350)
          ),
          box(title = "ACTH Normalization", width = 6, status = "danger",
            plotlyOutput("scen_acth_plot", height = 350)
          )
        ),
        fluidRow(
          box(title = "Scenario Outcomes at Day 90", width = 12,
            DTOutput("scen_summary_table")
          )
        )
      ),

      # ── Tab 7: Adrenal Crisis Risk ────────────────────────
      tabItem("tab_crisis",
        fluidRow(
          box(title = "Cumulative Adrenal Crisis Risk", width = 8, status = "danger",
            plotlyOutput("crisis_risk_plot", height = 400)
          ),
          box(title = "Crisis Key Facts", width = 4, status = "warning",
            HTML('<h4>Adrenal Crisis</h4>
                  <ul>
                    <li>Incidence: <b>5.2/100 patient-years</b></li>
                    <li>Mortality: <b>6–17%</b> per episode</li>
                    <li>Main triggers: infection (30%), GI illness, surgery</li>
                    <li>Prevention: sick-day rules, emergency kit</li>
                    <li>Treatment: HC 100 mg IV bolus + saline</li>
                  </ul>
                  <h4>Cortisol Threshold</h4>
                  <p>Crisis risk increases when total cortisol < 5 μg/dL</p>
                  <p>Emergency requirement: 200–400 mg/day HC IV</p>')
          )
        ),
        fluidRow(
          box(title = "Crisis Risk vs Cortisol Level", width = 6, status = "danger",
            plotlyOutput("crisis_dose_plot", height = 300)
          ),
          box(title = "Stress Dosing Simulation", width = 6, status = "info",
            plotlyOutput("crisis_stress_plot", height = 300)
          )
        )
      ),

      # ── Tab 8: Dose Optimization ──────────────────────────
      tabItem("tab_dose",
        fluidRow(
          box(title = "Cortisol AUC₀₋₂₄ Heatmap (HC Dose Optimization)", width = 12, status = "success",
            plotlyOutput("dose_auc_heatmap", height = 400)
          )
        ),
        fluidRow(
          box(title = "Optimal HC Split Ratio", width = 6, status = "primary",
            plotlyOutput("dose_ratio_plot", height = 300)
          ),
          box(title = "FC Dose — Electrolyte Response", width = 6, status = "info",
            plotlyOutput("dose_fc_plot", height = 300)
          )
        ),
        fluidRow(
          box(title = "Treatment Monitoring Targets", width = 12,
            HTML('<table style="width:100%;border-collapse:collapse;">
              <tr style="background:#1976D2;color:white;">
                <th>Parameter</th><th>Monitoring Interval</th><th>Target</th><th>Action if Out of Range</th>
              </tr>
              <tr><td>Cortisol day curve (CDC)</td><td>Annually / dose change</td><td>AUC 35–70 μg·h/dL</td><td>Adjust HC total dose/timing</td></tr>
              <tr><td>Plasma renin activity</td><td>Every 3–6 months</td><td>1–4 ng/mL/h</td><td>Adjust FC dose</td></tr>
              <tr><td>Serum electrolytes</td><td>Every 3–6 months</td><td>Na 138–142, K 3.5–5.0</td><td>Adjust FC dose</td></tr>
              <tr><td>DHEA-S</td><td>Annually</td><td>50–200 μg/dL</td><td>Add/adjust DHEA</td></tr>
              <tr><td>DXA (BMD)</td><td>Every 2–3 years</td><td>T-score > –1</td><td>Minimize HC dose, add bisphosphonate?</td></tr>
              <tr><td>HbA1c / fasting glucose</td><td>Annually</td><td>HbA1c <5.7%</td><td>Minimize HC dose</td></tr>
            </table>')
          )
        )
      )
    )
  )
)

# ──────────────────────────────────────────────────────────────
# SERVER
# ──────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Reactive simulation ──────────────────────────────────
  sim_data <- eventReactive(input$run_sim, {
    run_sim(
      hc_am      = input$hc_am,
      hc_nn      = input$hc_nn,
      hc_pm      = input$hc_pm,
      fc_mcg     = input$fc_mcg,
      use_dhea   = input$use_dhea,
      hc_mode    = input$hc_mode,
      days       = input$sim_days,
      stress     = input$stress_on,
      stress_start = input$stress_start,
      stress_days  = input$stress_days,
      stress_mult  = input$stress_mult
    )
  }, ignoreNULL = FALSE)

  # Run default on load
  observe({
    isolate({ input$run_sim })
  })

  # ── Value Boxes ──────────────────────────────────────────
  output$vbox_cortisol <- renderValueBox({
    df <- sim_data()
    val <- tail(df$Total_Cort, 1) %>% round(1)
    status <- ifelse(val < 5, "danger", ifelse(val > 20, "warning", "success"))
    valueBox(paste0(val, " μg/dL"), "Cortisol (final)", icon = icon("flask"), color = status)
  })

  output$vbox_acth <- renderValueBox({
    df <- sim_data()
    val <- tail(df$ACTH_CMT, 1) %>% round(1)
    status <- ifelse(val > 46, "danger", "success")
    valueBox(paste0(val, " pg/mL"), "ACTH (final)", icon = icon("dna"), color = status)
  })

  output$vbox_na <- renderValueBox({
    df <- sim_data()
    val <- tail(df$Na_CMT, 1) %>% round(1)
    status <- ifelse(val < 135, "danger", ifelse(val > 145, "warning", "success"))
    valueBox(paste0(val, " mEq/L"), "Serum Na⁺", icon = icon("tint"), color = status)
  })

  output$vbox_k <- renderValueBox({
    df <- sim_data()
    val <- tail(df$K_CMT, 1) %>% round(2)
    status <- ifelse(val > 5.5, "danger", ifelse(val < 3.5, "warning", "success"))
    valueBox(paste0(val, " mEq/L"), "Serum K⁺", icon = icon("bolt"), color = status)
  })

  # ── Static tables ────────────────────────────────────────
  output$dx_table <- renderTable({
    data.frame(
      Symptom = c("Fatigue/Weakness","Weight Loss","Hyperpigmentation","Salt Craving",
                  "Hypotension","Nausea/Vomiting","Abdominal Pain","Decreased Libido"),
      Frequency = c(">99%","87%","~80%","~80%","~90%","~75%","~31%","Female ~50%"),
      Mechanism = c("↓Cortisol","↓Cortisol/Aldo","↑ACTH→MSH","↓Aldosterone",
                    "↓Aldo/GC → ↓volume","↓GC","↓GC/electrolytes","↓DHEA")
    )
  })

  output$dx_criteria <- renderTable({
    data.frame(
      Test = c("8am Cortisol","Plasma ACTH","Cosyntropin Stim","Low-dose Stim","Anti-21OH Ab","Plasma Renin"),
      Diagnostic_Threshold = c("<3 μg/dL (= PAI)",">2× ULN (>80 pg/mL)","Peak <18 μg/dL","Peak <18 μg/dL",
                                "Positive (80% sensitivity)","↑ (>4 ng/mL/h)"),
      Interpretation = c("Confirms GC def","Primary vs secondary","Standard confirm","Sensitive confirm",
                          "Autoimmune cause","Mineralocorticoid def")
    )
  })

  # ── PK Plots ─────────────────────────────────────────────
  output$pk_hc_plot <- renderPlotly({
    df <- sim_data() %>% filter(Day <= 3)
    p <- ggplot(df, aes(x = time, y = Cp_HC_ugdL)) +
      geom_line(color = "#1976D2", linewidth = 1) +
      geom_hline(yintercept = c(3, 18), linetype = "dashed", color = "grey60") +
      labs(x = "Time (h)", y = "HC Plasma (μg/dL)",
           title = "HC Plasma Concentration (First 3 Days)") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_fc_plot <- renderPlotly({
    df <- sim_data() %>% filter(Day <= 10)
    p <- ggplot(df, aes(x = Day, y = Cp_FC_ngmL)) +
      geom_line(color = "#388E3C", linewidth = 1) +
      labs(x = "Day", y = "FC Plasma (ng/mL)",
           title = "Fludrocortisone Plasma Concentration") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_table <- renderDT({
    data.frame(
      Parameter    = c("HC Ka (h⁻¹)","HC Bioavailability","HC CL (L/h)","HC Vc (L)","HC t½ (h)",
                       "FC Ka (h⁻¹)","FC CL (L/h)","FC Vc (L)","FC t½ (h)"),
      Value        = c("1.20","96%","90","15","1.7","1.50","4.0","80","18–36"),
      Source       = c("Johansson 2009","—","Johansson 2009","—","Johansson 2009",
                       "Estimated","Estimated","—","Briggs 2006"),
      stringsAsFactors = FALSE
    ), options = list(pageLength = 12, dom = "t")
  })

  output$pk_auc_plot <- renderPlotly({
    df <- sim_data() %>% filter(Day >= max(Day) - 1)
    auc_approx <- sum(df$Total_Cort * 0.5)
    p <- ggplot(df, aes(x = TOD, y = Total_Cort)) +
      geom_area(fill = "#64B5F6", alpha = 0.4) +
      geom_line(color = "#1565C0", linewidth = 1) +
      geom_hline(yintercept = c(5, 25), linetype = "dashed", color = "grey60") +
      labs(x = "Time of Day (h)", y = "Cortisol (μg/dL)",
           title = paste0("24h Cortisol Profile (AUC≈", round(auc_approx,0), " μg·h/dL)")) +
      theme_bw()
    ggplotly(p)
  })

  output$pk_form_plot <- renderPlotly({
    df_ir <- run_sim(10,5,5,100,FALSE,"IR",3) %>% mutate(Form="IR")
    df_mr <- run_sim(20,0,0,100,FALSE,"MR",3) %>% mutate(Form="MR")
    df_both <- bind_rows(df_ir, df_mr) %>% filter(Day <= 3)
    p <- ggplot(df_both, aes(x=time, y=Total_Cort, color=Form)) +
      geom_line(linewidth=1) +
      scale_color_manual(values=c("IR"="#E53935","MR"="#1E88E5")) +
      labs(x="Time (h)", y="Cortisol (μg/dL)", title="IR vs MR-HC: Cortisol Profile") +
      theme_bw()
    ggplotly(p)
  })

  # ── HPA Axis Plots ────────────────────────────────────────
  output$hpa_acth_plot <- renderPlotly({
    df <- sim_data() %>% filter(TOD < 0.6)
    p <- ggplot(df, aes(x = Day, y = ACTH_CMT)) +
      geom_line(color = "#E53935", linewidth = 1) +
      geom_hline(yintercept = c(10, 46), linetype = "dashed", color = "grey60") +
      annotate("text", x = max(df$Day)*0.1, y = 47, label = "ULN = 46 pg/mL", size = 3) +
      labs(x = "Day", y = "ACTH (pg/mL)", title = "Plasma ACTH (Morning)") +
      theme_bw()
    ggplotly(p)
  })

  output$hpa_circ_plot <- renderPlotly({
    df <- sim_data() %>%
      filter(Day >= max(Day) - 2)
    p <- ggplot(df, aes(x = TOD, y = Total_Cort)) +
      geom_line(color = "#1976D2", linewidth = 1.2) +
      geom_ribbon(aes(ymin = 3, ymax = 20), fill = "#64B5F6", alpha = 0.1) +
      scale_x_continuous(breaks = seq(0, 24, 4)) +
      labs(x = "Time of Day (h)", y = "Cortisol (μg/dL)",
           title = "Circadian Cortisol Profile (Steady State)") +
      theme_bw()
    ggplotly(p)
  })

  output$hpa_crh_plot <- renderPlotly({
    df <- sim_data() %>% filter(TOD < 0.6)
    p <- ggplot(df, aes(x = Day, y = CRH)) +
      geom_line(color = "#FB8C00", linewidth = 1) +
      labs(x = "Day", y = "CRH (pmol/mL)", title = "Hypothalamic CRH") +
      theme_bw()
    ggplotly(p)
  })

  output$hpa_gr_plot <- renderPlotly({
    df <- sim_data() %>% filter(Day <= 7)
    p <- ggplot(df, aes(x = time)) +
      geom_line(aes(y = GR_FREE,  color = "Free GR"),  linewidth = 0.8) +
      geom_line(aes(y = GR_BOUND, color = "Bound GR"), linewidth = 0.8) +
      scale_color_manual(values = c("Free GR" = "#42A5F5", "Bound GR" = "#EF5350")) +
      labs(x = "Time (h)", y = "GR (AU)", title = "Glucocorticoid Receptor Occupancy") +
      theme_bw()
    ggplotly(p)
  })

  output$hpa_pig_plot <- renderPlotly({
    df <- sim_data() %>% filter(TOD < 0.6)
    p <- ggplot(df, aes(x = Day, y = ACTH_SIG)) +
      geom_line(color = "#8D6E63", linewidth = 1) +
      geom_hline(yintercept = 1.5, linetype = "dashed", color = "grey60") +
      labs(x = "Day", y = "ACTH/MSH Signal (AU)",
           title = "Hyperpigmentation Index\n(ACTH/MSH activity)") +
      theme_bw()
    ggplotly(p)
  })

  # ── Electrolyte Plots ─────────────────────────────────────
  output$elec_na_plot <- renderPlotly({
    df <- sim_data() %>% filter(TOD < 0.6)
    p <- ggplot(df, aes(x = Day, y = Na_CMT)) +
      geom_line(color = "#1976D2", linewidth = 1) +
      geom_hline(yintercept = c(135, 145), linetype = "dashed", color = "grey60") +
      annotate("rect", xmin=-Inf,xmax=Inf,ymin=135,ymax=145, fill="#81D4FA", alpha=0.15) +
      labs(x = "Day", y = "Na⁺ (mEq/L)", title = "Serum Sodium") +
      theme_bw()
    ggplotly(p)
  })

  output$elec_k_plot <- renderPlotly({
    df <- sim_data() %>% filter(TOD < 0.6)
    p <- ggplot(df, aes(x = Day, y = K_CMT)) +
      geom_line(color = "#F57C00", linewidth = 1) +
      geom_hline(yintercept = c(3.5, 5.0), linetype = "dashed", color = "grey60") +
      annotate("rect", xmin=-Inf,xmax=Inf,ymin=3.5,ymax=5.0, fill="#FFE082", alpha=0.15) +
      labs(x = "Day", y = "K⁺ (mEq/L)", title = "Serum Potassium") +
      theme_bw()
    ggplotly(p)
  })

  output$elec_map_plot <- renderPlotly({
    df <- sim_data() %>% filter(TOD < 0.6)
    p <- ggplot(df, aes(x = Day, y = MAP_CMT)) +
      geom_line(color = "#E53935", linewidth = 1) +
      geom_hline(yintercept = c(70, 100), linetype = "dashed", color = "grey60") +
      annotate("rect", xmin=-Inf,xmax=Inf,ymin=70,ymax=100, fill="#FFCDD2", alpha=0.15) +
      labs(x = "Day", y = "MAP (mmHg)", title = "Mean Arterial Pressure") +
      theme_bw()
    ggplotly(p)
  })

  # ── Clinical Endpoints ────────────────────────────────────
  output$clin_gluc_plot <- renderPlotly({
    df <- sim_data() %>% filter(TOD < 0.6)
    p <- ggplot(df, aes(x = Day, y = Gluc_CMT)) +
      geom_line(color = "#FBC02D", linewidth = 1) +
      geom_hline(yintercept = c(70, 110), linetype = "dashed", color = "grey60") +
      labs(x = "Day", y = "Glucose (mg/dL)", title = "Blood Glucose") +
      theme_bw()
    ggplotly(p)
  })

  output$clin_bmd_plot <- renderPlotly({
    df <- sim_data() %>% filter(TOD < 0.6)
    p <- ggplot(df, aes(x = Day, y = BMD_CMT)) +
      geom_line(color = "#7B1FA2", linewidth = 1) +
      geom_hline(yintercept = c(0.85, 1.0), linetype = "dashed", color = "grey60") +
      labs(x = "Day", y = "BMD (normalized)", title = "Bone Mineral Density") +
      theme_bw()
    ggplotly(p)
  })

  output$clin_ar_plot <- renderPlotly({
    df <- sim_data() %>% filter(TOD < 0.6)
    p <- ggplot(df, aes(x = Day, y = AR_CMT * 100)) +
      geom_line(color = "#00897B", linewidth = 1) +
      labs(x = "Day", y = "Adrenal Reserve (%)", title = "Functional Adrenal Reserve") +
      theme_bw()
    ggplotly(p)
  })

  output$clin_dhea_plot <- renderPlotly({
    df <- sim_data() %>% filter(TOD < 0.6)
    p <- ggplot(df, aes(x = Day, y = Cp_DHEA)) +
      geom_line(color = "#6A1B9A", linewidth = 1) +
      geom_hline(yintercept = c(50, 200), linetype = "dashed", color = "grey60") +
      labs(x = "Day", y = "DHEA-S (arb. units)", title = "DHEA Plasma Level") +
      theme_bw()
    ggplotly(p)
  })

  output$outcomes_table <- renderDT({
    df <- sim_data()
    last <- tail(df, 1)
    data.frame(
      Endpoint = c("Total Cortisol (μg/dL)","ACTH (pg/mL)","Na⁺ (mEq/L)",
                   "K⁺ (mEq/L)","MAP (mmHg)","Glucose (mg/dL)",
                   "BMD (normalized)","Crisis Risk (%)","Adrenal Reserve (%)"),
      Value    = round(c(last$Total_Cort, last$ACTH_CMT, last$Na_CMT,
                          last$K_CMT, last$MAP_CMT, last$Gluc_CMT,
                          last$BMD_CMT, last$Crisis_risk * 100,
                          last$AR_CMT * 100), 2),
      Target   = c("5–20","10–46","135–145","3.5–5.0","70–100",
                   "70–110","≥0.90","<5%","—"),
      Status   = c(
        ifelse(last$Total_Cort >= 5 & last$Total_Cort <= 20, "✅ OK", "⚠️ Check"),
        ifelse(last$ACTH_CMT < 46, "✅ OK", "⚠️ Elevated"),
        ifelse(last$Na_CMT >= 135 & last$Na_CMT <= 145, "✅ OK", "⚠️ Check"),
        ifelse(last$K_CMT >= 3.5 & last$K_CMT <= 5.0, "✅ OK", "⚠️ Check"),
        ifelse(last$MAP_CMT >= 70, "✅ OK", "⚠️ Low"),
        ifelse(last$Gluc_CMT >= 70, "✅ OK", "⚠️ Low"),
        ifelse(last$BMD_CMT >= 0.85, "✅ OK", "⚠️ Low"),
        ifelse(last$Crisis_risk < 0.05, "✅ Low", "⚠️ High"),
        "—"
      ), stringsAsFactors = FALSE
    )
  }, options = list(pageLength = 12, dom = "t"))

  # ── Scenario Comparison ───────────────────────────────────
  all_scenarios <- reactive({
    list(
      run_sim(0, 0, 0, 0,   FALSE, "IR", 90) %>% mutate(Scenario = "1. No Treatment"),
      run_sim(10,5, 5, 100, FALSE, "IR", 90) %>% mutate(Scenario = "2. Std HC+FC"),
      run_sim(20,0, 0, 100, FALSE, "MR", 90) %>% mutate(Scenario = "3. MR-HC+FC"),
      run_sim(10,5, 5, 100, TRUE,  "IR", 90) %>% mutate(Scenario = "4. HC+FC+DHEA"),
      run_sim(10,5, 5, 100, FALSE, "IR", 90, stress=TRUE, stress_start=20) %>%
                                               mutate(Scenario = "5. Stress Dosing")
    ) %>% bind_rows()
  })

  output$scen_cort_plot <- renderPlotly({
    df <- all_scenarios() %>% filter(TOD < 0.6)
    p <- ggplot(df, aes(x = Day, y = Total_Cort, color = Scenario)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = c(5, 20), linetype="dashed", color="grey60") +
      labs(x="Day", y="Cortisol (μg/dL)", title="Cortisol — 5 Scenarios") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$scen_acth_plot <- renderPlotly({
    df <- all_scenarios() %>% filter(TOD < 0.6)
    p <- ggplot(df, aes(x = Day, y = ACTH_CMT, color = Scenario)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 46, linetype="dashed", color="grey60") +
      labs(x="Day", y="ACTH (pg/mL)", title="ACTH Normalization") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$scen_summary_table <- renderDT({
    df <- all_scenarios()
    df %>%
      filter(abs(Day - 90) < 0.1) %>%
      group_by(Scenario) %>%
      summarise(
        `Cortisol (μg/dL)` = round(mean(Total_Cort), 1),
        `ACTH (pg/mL)`     = round(mean(ACTH_CMT), 1),
        `Na (mEq/L)`       = round(mean(Na_CMT), 1),
        `K (mEq/L)`        = round(mean(K_CMT), 2),
        `MAP (mmHg)`       = round(mean(MAP_CMT), 1),
        `Glucose (mg/dL)`  = round(mean(Gluc_CMT), 1),
        `Crisis Risk (%)`  = round(mean(Crisis_risk * 100), 1),
        .groups = "drop"
      )
  }, options = list(dom = "t"))

  # ── Adrenal Crisis Risk ────────────────────────────────────
  output$crisis_risk_plot <- renderPlotly({
    df <- all_scenarios() %>% filter(TOD < 0.6)
    p <- ggplot(df, aes(x = Day, y = Crisis_risk * 100, color = Scenario)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 5, linetype = "dashed", color = "red") +
      annotate("text", x = 5, y = 5.5, label = "5% threshold", size = 3, color = "red") +
      labs(x = "Day", y = "Crisis Risk (%)", title = "Cumulative Adrenal Crisis Risk") +
      theme_bw() + theme(legend.position = "bottom")
    ggplotly(p)
  })

  output$crisis_dose_plot <- renderPlotly({
    cort_seq <- seq(0, 25, 0.5)
    k_c <- 0.001; cs <- 5.0; kd <- 0.005
    risk <- sapply(cort_seq, function(c) {
      if (c < cs) k_c * (cs - c) / cs / kd else 0
    })
    p <- ggplot(data.frame(Cortisol = cort_seq, Risk = risk * 100),
                aes(x = Cortisol, y = Risk)) +
      geom_line(color = "#E53935", linewidth = 1.5) +
      geom_vline(xintercept = 5, linetype = "dashed", color = "grey60") +
      labs(x = "Total Cortisol (μg/dL)", y = "Steady-State Crisis Risk (%)",
           title = "Crisis Risk vs Cortisol Level") +
      theme_bw()
    ggplotly(p)
  })

  output$crisis_stress_plot <- renderPlotly({
    df_no_stress <- run_sim(10,5,5,100,FALSE,"IR",60,FALSE) %>%
      mutate(Type = "Normal Dosing") %>% filter(TOD < 0.6)
    df_stress <- run_sim(10,5,5,100,FALSE,"IR",60,TRUE,20,7,2) %>%
      mutate(Type = "Stress Dosing (2×)") %>% filter(TOD < 0.6)
    df_both <- bind_rows(df_no_stress, df_stress)
    p <- ggplot(df_both, aes(x = Day, y = Total_Cort, color = Type)) +
      geom_line(linewidth = 1) +
      geom_vline(xintercept = c(20,27), linetype="dashed", color="grey40") +
      annotate("rect", xmin=20,xmax=27,ymin=-Inf,ymax=Inf, fill="yellow", alpha=0.2) +
      annotate("text", x=23.5, y=max(df_both$Total_Cort)*0.9,
               label="Illness\nEpisode", size=3) +
      scale_color_manual(values = c("Normal Dosing"="#1E88E5","Stress Dosing (2×)"="#E53935")) +
      labs(x="Day", y="Cortisol (μg/dL)", title="Stress Dosing Protocol: Cortisol Coverage") +
      theme_bw()
    ggplotly(p)
  })

  # ── Dose Optimization ──────────────────────────────────────
  output$dose_auc_heatmap <- renderPlotly({
    doses <- c(10, 15, 20, 25, 30)
    splits <- list(c(10,0,0), c(6,2,2), c(5,2.5,2.5), c(4,2,4),
                   c(7,4,4), c(8,4,4), c(10,4,6), c(10,5,5),
                   c(12,5,5), c(12,6,7))
    auc_res <- lapply(splits, function(sp) {
      total_hc <- sum(sp)
      d <- run_sim(sp[1], sp[2], sp[3], 100, FALSE, "IR", 10)
      d_ss <- d %>% filter(Day >= 8)
      auc <- sum(d_ss$Total_Cort) * 0.5 / 2  # AU/day
      data.frame(AM = sp[1], Noon = sp[2], PM = sp[3], Total = total_hc, AUC24 = auc)
    }) %>% bind_rows()

    p <- ggplot(auc_res, aes(x = factor(paste0(AM,"/",Noon,"/",PM)), y = AUC24, fill = AUC24)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = c(35, 70), linetype = "dashed", color = "darkred") +
      scale_fill_gradient2(low="#1976D2", mid="#FFEE58", high="#E53935", midpoint=52.5) +
      labs(x = "HC Dose Split (AM/Noon/PM, mg)", y = "Cortisol AUC₀₋₂₄ (μg·h/dL)",
           title = "Cortisol AUC by HC Dose Split — Target 35–70 μg·h/dL") +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggplotly(p)
  })

  output$dose_ratio_plot <- renderPlotly({
    df <- sim_data() %>% filter(Day >= max(Day) - 2)
    p <- ggplot(df, aes(x = TOD, y = Total_Cort)) +
      geom_line(color = "#1976D2", linewidth = 1.2) +
      geom_hline(yintercept = c(3, 20), linetype = "dashed", color = "grey60") +
      scale_x_continuous(breaks = seq(0, 24, 4)) +
      labs(x = "Time of Day (h)", y = "Cortisol (μg/dL)",
           title = "Cortisol Profile — Current Dosing") +
      theme_bw()
    ggplotly(p)
  })

  output$dose_fc_plot <- renderPlotly({
    fc_doses <- c(0, 50, 100, 150, 200, 300)
    fc_res <- lapply(fc_doses, function(fc) {
      d <- run_sim(10,5,5,fc,FALSE,"IR",30) %>% tail(48)
      data.frame(FC_dose = fc,
                 Na_mean = mean(d$Na_CMT),
                 K_mean  = mean(d$K_CMT),
                 MAP_mean= mean(d$MAP_CMT))
    }) %>% bind_rows()

    p <- ggplot(fc_res) +
      geom_line(aes(x = FC_dose, y = Na_mean, color = "Na⁺"), linewidth=1) +
      geom_line(aes(x = FC_dose, y = K_mean*10, color = "K⁺ ×10"), linewidth=1) +
      geom_line(aes(x = FC_dose, y = MAP_mean/2, color = "MAP/2"), linewidth=1) +
      scale_color_manual(values = c("Na⁺"="#1976D2","K⁺ ×10"="#F57C00","MAP/2"="#E53935")) +
      geom_hline(yintercept = c(135, 145), linetype="dashed", color="lightblue") +
      labs(x = "FC Dose (μg/day)", y = "Scaled Values",
           title = "FC Dose-Response: Na, K, MAP") +
      theme_bw()
    ggplotly(p)
  })
}

shinyApp(ui, server)
