##############################################################################
# EoE (Eosinophilic Esophagitis) Interactive QSP Shiny Dashboard
#
# Tabs:
#   1. Patient Profile & Disease Inputs
#   2. Drug PK — Concentration-Time Profiles
#   3. Cytokine & Immune Response (IL-13, IL-5, Eotaxin-3, IgE)
#   4. Eosinophil Dynamics (Blood & Tissue)
#   5. Clinical Endpoints (Dysphagia, EREFS, Histology, Barrier)
#   6. Treatment Scenario Comparison (Side-by-Side)
#   7. Biomarker Panel & Heatmap
#
# Author: Claude Code Routine — QSP Disease Model Library
# Date: 2026-06-24
##############################################################################

library(shiny)
library(shinydashboard)
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(tidyr)
library(scales)

# ============================================================
# EMBEDDED MODEL CODE (same as standalone R model)
# ============================================================

eoe_model_code <- '
$PARAM @annotated
kin_IL13  : 320   : IL-13 production rate (pg/mL/day)
kout_IL13 : 4.0   : IL-13 elimination (1/day)
kin_IL5   : 60    : IL-5 production (pg/mL/day)
kout_IL5  : 4.0   : IL-5 elimination (1/day)
kin_eotax : 1600  : Eotaxin-3 production (pg/mL/day)
kout_eotax: 4.0   : Eotaxin-3 elimination (1/day)
hill_eotax: 1.5   : IL-13→eotaxin-3 Hill coefficient
kin_eosbl : 2400  : Blood EOS production (cells/µL/day)
kout_eosbl: 4.0   : Blood EOS turnover (1/day)
k_recruit : 0.06  : Tissue recruitment (1/day)
kout_eos_eso: 0.05: Tissue EOS apoptosis (1/day)
kin_mast  : 0.10  : Mast cell influx (cells·mm⁻²/day)
kout_mast : 0.002 : Mast cell turnover (1/day)
kin_fibro : 0.001 : Fibrosis accumulation (1/day)
kout_fibro: 0.0008: Fibrosis resolution (1/day)
kin_IgE   : 0.15  : IgE production (IU/mL/day)
kout_IgE  : 0.0005: IgE elimination (1/day)
tau_EPBAR : 7     : Epithelial barrier time constant (days)
Emax_IL13_eotax : 0.95 : Max IL-13→eotaxin-3
IC50_IL13_eotax : 40   : IL-13 EC50 for eotaxin-3 (pg/mL)
Emax_eotax_recruit: 0.90: Max eotaxin-3→recruitment
IC50_eotax_recruit: 200 : Eotaxin-3 EC50 (pg/mL)
Emax_IL5_eosbl  : 0.80  : Max IL-5→BM eosinopoiesis
IC50_IL5_eosbl  : 8     : IL-5 EC50 (pg/mL)
Emax_eos_fibro  : 0.70  : Max eos→fibrosis
IC50_eos_fibro  : 40    : Eos EC50 for fibrosis (eos/hpf)
Emax_IL13_bar   : 0.70  : Max IL-13 barrier disruption
IC50_IL13_bar   : 40    : IL-13 EC50 barrier (pg/mL)
Emax_eos_bar    : 0.50  : Max eos barrier disruption
IC50_eos_bar    : 50    : Eos EC50 barrier (eos/hpf)
ka_bud_sys  : 0.6   : BUD esoph→systemic (1/day)
ke_bud_eso  : 2.4   : BUD esoph elim (1/day)
ke_bud_sys  : 8.0   : BUD systemic elim (1/day)
Emax_bud_IL13: 0.85 : BUD max IL-13 suppression
IC50_bud_IL13: 0.1  : BUD IC50 IL-13 (mg/L)
Emax_bud_eotax:0.90 : BUD max eotaxin suppression
IC50_bud_eotax:0.08 : BUD IC50 eotaxin (mg/L)
Emax_bud_eos_eso: 0.85: BUD max eos apoptosis
IC50_bud_eos_eso:0.1 : BUD IC50 eos esoph (mg/L)
ka_dup    : 0.18  : DUP ka SC (1/day)
F_dup     : 0.64  : DUP SC bioavailability
CL_dup    : 0.21  : DUP clearance (L/day)
Vd_dup_C  : 3.5   : DUP central volume (L)
Vd_dup_P  : 2.8   : DUP peripheral volume (L)
Q_dup     : 1.5   : DUP Q (L/day)
Emax_dup_STAT6: 0.95: DUP max STAT6 blockade
IC50_dup_STAT6: 2.0 : DUP IC50 STAT6 (mg/L)
Emax_dup_IgE: 0.70  : DUP max IgE reduction
IC50_dup_IgE: 10    : DUP IC50 IgE (mg/L)
ka_mepo   : 0.34  : MEPO ka SC (1/day)
F_mepo    : 0.81  : MEPO bioavailability
CL_mepo   : 0.28  : MEPO clearance (L/day)
Vd_mepo   : 3.6   : MEPO volume (L)
Emax_mepo_IL5: 0.95: MEPO max IL-5 block
IC50_mepo_IL5: 1.0  : MEPO IC50 IL-5 (mg/L)
ka_cenda  : 14.4  : CENDA ka (1/day)
F_cenda   : 0.35  : CENDA bioavailability
CL_cenda  : 360   : CENDA apparent CL (L/day)
Vd_cenda  : 100   : CENDA apparent Vd (L)
Emax_cenda_IL13: 0.90: CENDA max IL-13 block
IC50_cenda_IL13: 0.05 : CENDA IC50 IL-13 (mg/L)
DIETARY : 0 : Dietary effect (0=none, 1=active)
Emax_diet_antigen: 0.80: Max dietary antigen reduction

$CMT @annotated
BUD_ESO  : Budesonide esoph (mg/L)
BUD_SYS  : Budesonide systemic (mg/L)
DUP_SC   : Dupilumab SC depot (mg)
DUP_C    : Dupilumab central (mg/L)
DUP_P    : Dupilumab peripheral (mg/L)
MEPO_SC  : Mepolizumab SC depot (mg)
MEPO_C   : Mepolizumab central (mg/L)
CENDA_GUT: Cendakimab gut (mg)
CENDA_C  : Cendakimab central (mg/L)
IL13     : Esophageal IL-13 (pg/mL)
IL5      : Circulating IL-5 (pg/mL)
EOTAX3   : Eotaxin-3/CCL26 (pg/mL)
EOS_BL   : Blood eosinophils (cells/µL)
EOS_ESO  : Tissue eosinophils (eos/hpf)
MAST_ESO : Esoph mast cells (per mm²)
FIBRO    : LP fibrosis (0-1)
IGE_TOT  : Total IgE (IU/mL)
EPBAR    : Barrier integrity (0-1)

$MAIN
double Inh_bud_IL13  = Emax_bud_IL13  * BUD_ESO / (BUD_ESO + IC50_bud_IL13);
double Inh_bud_eotax = Emax_bud_eotax * BUD_ESO / (BUD_ESO + IC50_bud_eotax);
double Enh_bud_apop  = Emax_bud_eos_eso * BUD_ESO / (BUD_ESO + IC50_bud_eos_eso);
double Inh_dup_STAT6 = Emax_dup_STAT6 * DUP_C / (DUP_C + IC50_dup_STAT6);
double Inh_dup_IgE   = Emax_dup_IgE   * DUP_C / (DUP_C + IC50_dup_IgE);
double Inh_mepo_IL5 = Emax_mepo_IL5 * MEPO_C / (MEPO_C + IC50_mepo_IL5);
double Inh_cenda_IL13 = Emax_cenda_IL13 * CENDA_C / (CENDA_C + IC50_cenda_IL13);
double Emax_diet = DIETARY * Emax_diet_antigen;
double Prod_IL13_factor = (1.0 - Inh_bud_IL13) * (1.0 - Emax_diet);
double IL13_signal = IL13 * (1.0 - Inh_dup_STAT6) * (1.0 - Inh_cenda_IL13);
double Stim_eotax = Emax_IL13_eotax * pow(IL13_signal, hill_eotax) /
                    (pow(IC50_IL13_eotax, hill_eotax) + pow(IL13_signal, hill_eotax));
double Prod_eotax = kin_eotax * (1.0 + Stim_eotax) * (1.0 - Inh_bud_eotax);
double IL5_effective = IL5 * (1.0 - Inh_mepo_IL5);
double Stim_eos_IL5 = Emax_IL5_eosbl * IL5_effective / (IC50_IL5_eosbl + IL5_effective);
double Stim_recruit = Emax_eotax_recruit * EOTAX3 / (IC50_eotax_recruit + EOTAX3);
double k_recruit_eff = k_recruit * (1.0 + Stim_recruit);
double Stim_fibro = Emax_eos_fibro * EOS_ESO / (IC50_eos_fibro + EOS_ESO);
double Disrupt_IL13 = Emax_IL13_bar * IL13 / (IC50_IL13_bar + IL13);
double Disrupt_eos  = Emax_eos_bar  * EOS_ESO / (IC50_eos_bar + EOS_ESO);
double EPBAR_ss = (1.0 - Disrupt_IL13 * (1.0 - Inh_dup_STAT6) * (1.0 - Inh_cenda_IL13) * (1.0 - Inh_bud_IL13)) *
                  (1.0 - 0.5 * Disrupt_eos);
EPBAR_ss = EPBAR_ss < 0.05 ? 0.05 : EPBAR_ss;
double FIBRO_FRAC = FIBRO / (FIBRO + 0.5);
double EREFS_val = 3.0 * EOS_ESO/(EOS_ESO+50.0) + 1.6 * EOS_ESO/(EOS_ESO+50.0) +
                   1.2 * (1.0-EPBAR) + 1.8 * FIBRO_FRAC + 1.4 * FIBRO_FRAC;
if(EREFS_val > 18.0) EREFS_val = 18.0;

$ODE
dxdt_BUD_ESO = -ka_bud_sys * BUD_ESO - ke_bud_eso * BUD_ESO;
dxdt_BUD_SYS = ka_bud_sys * BUD_ESO - ke_bud_sys * BUD_SYS;
dxdt_DUP_SC = -ka_dup * DUP_SC;
dxdt_DUP_C  = ka_dup * F_dup * DUP_SC / Vd_dup_C - (CL_dup/Vd_dup_C)*DUP_C - (Q_dup/Vd_dup_C)*DUP_C + (Q_dup/Vd_dup_P)*DUP_P;
dxdt_DUP_P  = (Q_dup/Vd_dup_C)*DUP_C - (Q_dup/Vd_dup_P)*DUP_P;
dxdt_MEPO_SC = -ka_mepo * MEPO_SC;
dxdt_MEPO_C  = ka_mepo * F_mepo * MEPO_SC / Vd_mepo - (CL_mepo/Vd_mepo)*MEPO_C;
dxdt_CENDA_GUT = -ka_cenda * CENDA_GUT;
dxdt_CENDA_C   = ka_cenda * F_cenda * CENDA_GUT / Vd_cenda - (CL_cenda/Vd_cenda)*CENDA_C;
dxdt_IL13 = kin_IL13 * Prod_IL13_factor - kout_IL13 * IL13 * (1.0 + Inh_cenda_IL13 * 2.0);
dxdt_IL5  = kin_IL5 * (1.0 - Inh_bud_IL13*0.5) * (1.0 - Emax_diet*0.7) - kout_IL5 * IL5 * (1.0 + Inh_mepo_IL5*3.0);
dxdt_EOTAX3 = Prod_eotax * (1.0 - Inh_dup_STAT6) - kout_eotax * EOTAX3;
dxdt_EOS_BL = kin_eosbl * (1.0 + Stim_eos_IL5) * (1.0 - Inh_bud_IL13*0.3) - kout_eosbl * EOS_BL - k_recruit_eff * EOS_BL;
dxdt_EOS_ESO = k_recruit_eff * EOS_BL - kout_eos_eso * EOS_ESO * (1.0 + Enh_bud_apop*3.0);
dxdt_MAST_ESO = kin_mast * (1.0 + IL13_signal/80.0*0.5) * (1.0 - Emax_diet*0.5) - kout_mast * MAST_ESO * (1.0 + Inh_dup_STAT6*0.5);
dxdt_FIBRO = kin_fibro * (1.0 + Stim_fibro) * (1.0 - Inh_dup_STAT6*0.4) * (1.0 - Inh_cenda_IL13*0.3) - kout_fibro * FIBRO;
dxdt_IGE_TOT = kin_IgE * (1.0 - Emax_diet*0.3) - kout_IgE * IGE_TOT * (1.0 + Inh_dup_IgE);
dxdt_EPBAR = (EPBAR_ss - EPBAR) / tau_EPBAR;

$TABLE
double EREFS_SCORE = EREFS_val;
double HISTO_REMIS = (EOS_ESO < 15.0) ? 1.0 : 0.0;
double DYSPHAG_SCORE = 10.0 * (0.5*EOS_ESO/(EOS_ESO+60.0) + 0.3*FIBRO/(FIBRO+0.5) + 0.2*(1.0-EPBAR));
double IL13_SIGNAL = IL13_signal;

$CAPTURE BUD_ESO BUD_SYS DUP_C DUP_P MEPO_C CENDA_C
IL13 IL5 EOTAX3 IL13_SIGNAL
EOS_BL EOS_ESO MAST_ESO EPBAR FIBRO IGE_TOT
EREFS_SCORE HISTO_REMIS DYSPHAG_SCORE
'

# ============================================================
# SIMULATION HELPERS
# ============================================================

compile_model <- function() {
  mcode("EoE_QSP_Shiny", eoe_model_code)
}

build_dose_events <- function(drug, dose, interval_days, sim_days) {
  if (drug == "none") return(data.frame())
  dose_times <- seq(0, sim_days, by = interval_days)
  cmt_map <- c(budesonide = 1, dupilumab = 3, mepolizumab = 6, cendakimab = 8)
  data.frame(time = dose_times, amt = dose,
             cmt = cmt_map[drug], evid = 1, rate = 0, ii = 0)
}

sim_scenario <- function(mod, drug1, dose1, int1,
                          drug2, dose2, int2,
                          dietary, sim_days,
                          bl_IL13 = 80, bl_EOS_ESO = 80,
                          bl_FIBRO = 0.4, bl_IGE = 300) {
  ev1 <- build_dose_events(drug1, dose1, int1, sim_days)
  ev2 <- if (drug2 != "none") build_dose_events(drug2, dose2, int2, sim_days)
          else data.frame()
  events <- bind_rows(ev1, ev2)
  events <- events[order(events$time), ]

  init_vals <- list(IL13 = bl_IL13, EOS_ESO = bl_EOS_ESO,
                    FIBRO = bl_FIBRO, IGE_TOT = bl_IGE,
                    EOTAX3 = 400 * bl_IL13 / 80,
                    EOS_BL  = 600 * bl_EOS_ESO / 80,
                    EPBAR = 1.0 - 0.6 * bl_EOS_ESO / 80)

  mod2 <- mod %>%
    param(DIETARY = as.numeric(dietary)) %>%
    init(IL13 = init_vals$IL13,
         EOS_ESO = init_vals$EOS_ESO,
         FIBRO = init_vals$FIBRO,
         IGE_TOT = init_vals$IGE_TOT,
         EOTAX3 = init_vals$EOTAX3,
         EOS_BL = init_vals$EOS_BL,
         EPBAR = init_vals$EPBAR,
         MAST_ESO = 50 * bl_EOS_ESO / 80)

  if (nrow(events) == 0) {
    out <- mrgsim(mod2, end = sim_days, delta = 7)
  } else {
    out <- mrgsim_df(mod2, events = as.data.frame(events),
                     end = sim_days, add = seq(0, sim_days, by = 7))
  }
  as.data.frame(out) %>% mutate(week = time / 7)
}

# ============================================================
# UI
# ============================================================

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = "EoE QSP Dashboard",
    titleWidth = 280
  ),

  dashboardSidebar(
    width = 280,
    sidebarMenu(
      menuItem("Patient Profile", tabName = "patient", icon = icon("user")),
      menuItem("Drug PK Profiles", tabName = "pk", icon = icon("flask")),
      menuItem("Cytokine & Immune Response", tabName = "cytokines", icon = icon("dna")),
      menuItem("Eosinophil Dynamics", tabName = "eosinophil", icon = icon("microscope")),
      menuItem("Clinical Endpoints", tabName = "endpoints", icon = icon("chart-line")),
      menuItem("Scenario Comparison", tabName = "comparison", icon = icon("columns")),
      menuItem("Biomarker Panel", tabName = "biomarker", icon = icon("vial"))
    ),

    hr(),

    # PRIMARY TREATMENT
    h5("Primary Treatment", style = "padding-left:15px; color:white; font-weight:bold"),
    selectInput("drug1", "Drug 1",
                choices = c("None" = "none",
                            "Budesonide ODT (1 mg BID)" = "budesonide",
                            "Dupilumab 300 mg SC q2w" = "dupilumab",
                            "Mepolizumab 300 mg SC q4w" = "mepolizumab",
                            "Cendakimab 160 mg PO QD" = "cendakimab"),
                selected = "dupilumab"),
    conditionalPanel(
      condition = "input.drug1 != 'none'",
      numericInput("dose1", "Dose (mg)", value = 300, min = 1, max = 1000),
      numericInput("interval1", "Dosing Interval (days)", value = 14, min = 1, max = 28)
    ),

    # COMBINATION DRUG
    h5("Combination Drug", style = "padding-left:15px; color:white; font-weight:bold"),
    selectInput("drug2", "Drug 2 (optional)",
                choices = c("None" = "none",
                            "Budesonide ODT (1 mg BID)" = "budesonide",
                            "Dupilumab 300 mg SC q2w" = "dupilumab",
                            "Mepolizumab 300 mg SC q4w" = "mepolizumab",
                            "Cendakimab 160 mg PO QD" = "cendakimab"),
                selected = "none"),

    checkboxInput("dietary", "Dietary Elimination (SFED)", value = FALSE),

    hr(),

    # PATIENT PARAMETERS
    h5("Patient Parameters", style = "padding-left:15px; color:white; font-weight:bold"),
    sliderInput("bl_eos", "Baseline Tissue Eos (eos/hpf)",
                min = 15, max = 300, value = 80, step = 5),
    sliderInput("bl_IL13", "Baseline IL-13 (pg/mL)",
                min = 10, max = 200, value = 80, step = 5),
    sliderInput("bl_fibro", "Baseline Fibrosis (0-1)",
                min = 0, max = 1, value = 0.4, step = 0.05),
    sliderInput("bl_IgE", "Baseline Total IgE (IU/mL)",
                min = 50, max = 1000, value = 300, step = 25),

    sliderInput("sim_weeks", "Simulation Duration (weeks)",
                min = 12, max = 104, value = 52, step = 4),

    actionButton("run", "Run Simulation",
                 class = "btn-primary btn-block",
                 style = "margin:10px; width:90%")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f4f6f9; }
      .box-header { font-weight: bold; }
    "))),

    tabItems(

      # ──────────────────────────────────────────────────────────────
      # TAB 1: PATIENT PROFILE
      # ──────────────────────────────────────────────────────────────
      tabItem(tabName = "patient",
        fluidRow(
          box(title = "EoE Disease Overview", width = 12, solidHeader = TRUE,
              status = "primary",
              HTML("
                <div style='font-size:14px'>
                <h4>Eosinophilic Esophagitis (EoE)</h4>
                <p><b>Definition:</b> A chronic, immune-mediated disease characterized by
                eosinophilic infiltration of the esophagus (≥15 eos/hpf) causing symptoms of
                esophageal dysfunction (dysphagia, food impaction, heartburn).</p>
                <div style='display:grid; grid-template-columns:1fr 1fr 1fr; gap:20px'>
                  <div style='background:#e3f2fd; padding:15px; border-radius:8px'>
                    <b>Pathophysiology</b><br>
                    Th2-dominant inflammation driven by food allergens → IL-4/IL-5/IL-13 axis →
                    eotaxin-3 (CCL26) production → CCR3+ eosinophil recruitment → epithelial
                    disruption (↓DSG1, ↓filaggrin) → fibrosis and stricture
                  </div>
                  <div style='background:#fce4ec; padding:15px; border-radius:8px'>
                    <b>Epidemiology</b><br>
                    Prevalence ~1/2000 in Western countries<br>
                    Male predominance (M:F = 3:1)<br>
                    Strong atopic comorbidity (75% have asthma/rhinitis)<br>
                    Peak incidence ages 20-40
                  </div>
                  <div style='background:#e8f5e9; padding:15px; border-radius:8px'>
                    <b>Treatment Goals</b><br>
                    • Histological remission: peak eos/hpf &lt;15<br>
                    • Dysphagia improvement (DSQ)<br>
                    • Endoscopic improvement (EREFS)<br>
                    • Quality of life improvement<br>
                    • Prevention of stricture progression
                  </div>
                </div>
                <br>
                <table border='1' style='width:100%; border-collapse:collapse; font-size:13px'>
                  <tr style='background:#1565C0; color:white'>
                    <th style='padding:8px'>Treatment</th>
                    <th style='padding:8px'>Mechanism</th>
                    <th style='padding:8px'>Histological Remission</th>
                    <th style='padding:8px'>Key Trial</th>
                  </tr>
                  <tr><td style='padding:6px'>Swallowed Budesonide</td>
                    <td>GR-mediated pan-cytokine suppression</td>
                    <td>~58% (ODT, 6 weeks)</td><td>ApplE (Lucendo 2022)</td></tr>
                  <tr style='background:#f5f5f5'><td style='padding:6px'>Dupilumab 300 mg SC q2w</td>
                    <td>Anti-IL-4Rα (blocks IL-4 + IL-13)</td>
                    <td>~60% (24 weeks)</td><td>MATS (Hirano 2022)</td></tr>
                  <tr><td style='padding:6px'>Cendakimab 160 mg QD</td>
                    <td>Anti-IL-13 (oral, direct neutralization)</td>
                    <td>~64% (24 weeks)</td><td>CACTUS (Hirano 2023)</td></tr>
                  <tr style='background:#f5f5f5'><td style='padding:6px'>Mepolizumab 300 mg SC q4w</td>
                    <td>Anti-IL-5 (blood EOS depletion)</td>
                    <td>~26% (phase 2)</td><td>Stein 2006, Rothenberg 2014</td></tr>
                  <tr><td style='padding:6px'>6-Food Elimination (SFED)</td>
                    <td>Allergen avoidance</td>
                    <td>~72%</td><td>Molina-Infante 2018</td></tr>
                </table>
                </div>
              ")
          )
        ),
        fluidRow(
          valueBoxOutput("vb_eos_esophagus"),
          valueBoxOutput("vb_dysphag"),
          valueBoxOutput("vb_erefs")
        )
      ),

      # ──────────────────────────────────────────────────────────────
      # TAB 2: DRUG PK
      # ──────────────────────────────────────────────────────────────
      tabItem(tabName = "pk",
        fluidRow(
          box(title = "Drug Concentration–Time Profiles", width = 8,
              solidHeader = TRUE, status = "primary",
              plotlyOutput("pk_plot", height = "450px")),
          box(title = "PK Parameters", width = 4,
              solidHeader = TRUE, status = "info",
              tableOutput("pk_table"))
        ),
        fluidRow(
          box(title = "Budesonide — Esophageal vs. Systemic Exposure", width = 6,
              solidHeader = TRUE, status = "warning",
              plotlyOutput("bud_pk_plot", height = "350px")),
          box(title = "Cendakimab — Oral PK (Daily Dosing)", width = 6,
              solidHeader = TRUE, status = "success",
              plotlyOutput("cenda_pk_plot", height = "350px"))
        )
      ),

      # ──────────────────────────────────────────────────────────────
      # TAB 3: CYTOKINE & IMMUNE RESPONSE
      # ──────────────────────────────────────────────────────────────
      tabItem(tabName = "cytokines",
        fluidRow(
          box(title = "IL-13 Dynamics (Key EoE Driver)", width = 6,
              solidHeader = TRUE, status = "warning",
              plotlyOutput("il13_plot", height = "350px")),
          box(title = "IL-5 Dynamics (Eosinopoietin)", width = 6,
              solidHeader = TRUE, status = "warning",
              plotlyOutput("il5_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "Eotaxin-3 / CCL26 (Tissue Eosinophil Recruiter)", width = 6,
              solidHeader = TRUE, status = "primary",
              plotlyOutput("eotax3_plot", height = "350px")),
          box(title = "Total Serum IgE", width = 6,
              solidHeader = TRUE, status = "info",
              plotlyOutput("ige_plot", height = "350px"))
        )
      ),

      # ──────────────────────────────────────────────────────────────
      # TAB 4: EOSINOPHIL DYNAMICS
      # ──────────────────────────────────────────────────────────────
      tabItem(tabName = "eosinophil",
        fluidRow(
          box(title = "Tissue Eosinophils (eos/hpf) — Primary Histological Endpoint",
              width = 8, solidHeader = TRUE, status = "danger",
              plotlyOutput("eos_tissue_plot", height = "400px")),
          box(title = "Histological Remission",
              width = 4, solidHeader = TRUE, status = "success",
              plotlyOutput("histo_pie", height = "250px"),
              hr(),
              h5("Time to Histological Remission (<15 eos/hpf)"),
              verbatimTextOutput("remission_time"))
        ),
        fluidRow(
          box(title = "Blood Eosinophils (Absolute Count)", width = 6,
              solidHeader = TRUE, status = "warning",
              plotlyOutput("eos_blood_plot", height = "350px")),
          box(title = "Esophageal Mast Cells", width = 6,
              solidHeader = TRUE, status = "info",
              plotlyOutput("mast_plot", height = "350px"))
        )
      ),

      # ──────────────────────────────────────────────────────────────
      # TAB 5: CLINICAL ENDPOINTS
      # ──────────────────────────────────────────────────────────────
      tabItem(tabName = "endpoints",
        fluidRow(
          box(title = "Dysphagia Score (0–10)", width = 6,
              solidHeader = TRUE, status = "primary",
              plotlyOutput("dysphag_plot", height = "350px")),
          box(title = "EREFS Score (0–18)", width = 6,
              solidHeader = TRUE, status = "primary",
              plotlyOutput("erefs_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "Epithelial Barrier Integrity (0=disrupted, 1=intact)", width = 6,
              solidHeader = TRUE, status = "warning",
              plotlyOutput("epbar_plot", height = "350px")),
          box(title = "Lamina Propria Fibrosis (0–1)", width = 6,
              solidHeader = TRUE, status = "danger",
              plotlyOutput("fibro_plot", height = "350px"))
        ),
        fluidRow(
          box(title = "Summary Table — Weeks 12 & 24", width = 12,
              solidHeader = TRUE, status = "info",
              DTOutput("endpoints_table"))
        )
      ),

      # ──────────────────────────────────────────────────────────────
      # TAB 6: SCENARIO COMPARISON
      # ──────────────────────────────────────────────────────────────
      tabItem(tabName = "comparison",
        fluidRow(
          box(title = "Tissue Eosinophils — All Treatments", width = 12,
              solidHeader = TRUE, status = "primary",
              plotlyOutput("compare_eos_plot", height = "450px"))
        ),
        fluidRow(
          box(title = "Dysphagia Comparison", width = 6,
              solidHeader = TRUE, status = "info",
              plotlyOutput("compare_dysphag", height = "350px")),
          box(title = "Fibrosis Comparison", width = 6,
              solidHeader = TRUE, status = "warning",
              plotlyOutput("compare_fibro", height = "350px"))
        ),
        fluidRow(
          box(title = "Head-to-Head at Week 24", width = 12,
              solidHeader = TRUE, status = "success",
              plotlyOutput("compare_bar", height = "350px"))
        )
      ),

      # ──────────────────────────────────────────────────────────────
      # TAB 7: BIOMARKER PANEL
      # ──────────────────────────────────────────────────────────────
      tabItem(tabName = "biomarker",
        fluidRow(
          box(title = "Biomarker Response Heatmap (% Change from Baseline at Week 24)",
              width = 12, solidHeader = TRUE, status = "primary",
              plotlyOutput("biomarker_heatmap", height = "450px"))
        ),
        fluidRow(
          box(title = "Biomarker Correlation Matrix", width = 6,
              solidHeader = TRUE, status = "info",
              plotlyOutput("biomarker_corr", height = "400px")),
          box(title = "Biomarker Trajectories", width = 6,
              solidHeader = TRUE, status = "success",
              plotlyOutput("biomarker_all", height = "400px"))
        )
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  # Compile model once
  mod <- reactive({
    withProgress(message = "Compiling EoE QSP model...", {
      compile_model()
    })
  })

  # Primary simulation
  sim_data <- eventReactive(input$run, {
    req(mod())
    withProgress(message = "Running simulation...", {
      sim_scenario(
        mod = mod(),
        drug1 = input$drug1,
        dose1 = input$dose1,
        int1  = input$interval1,
        drug2 = input$drug2,
        dose2 = if (input$drug2 == "budesonide") 1 else
                if (input$drug2 == "cendakimab") 160 else 300,
        int2  = if (input$drug2 == "budesonide") 0.5 else
                if (input$drug2 == "cendakimab") 1 else
                if (input$drug2 == "mepolizumab") 28 else 14,
        dietary = input$dietary,
        sim_days = input$sim_weeks * 7,
        bl_IL13 = input$bl_IL13,
        bl_EOS_ESO = input$bl_eos,
        bl_FIBRO = input$bl_fibro,
        bl_IGE = input$bl_IgE
      )
    })
  })

  # Multi-scenario comparison simulation
  compare_data <- eventReactive(input$run, {
    req(mod())
    withProgress(message = "Running comparison scenarios...", {
      scenarios_list <- list(
        list(d1 = "none",        d2 = "none",       lbl = "No Treatment"),
        list(d1 = "budesonide",  d2 = "none",       lbl = "Budesonide ODT"),
        list(d1 = "dupilumab",   d2 = "none",       lbl = "Dupilumab"),
        list(d1 = "mepolizumab", d2 = "none",       lbl = "Mepolizumab"),
        list(d1 = "cendakimab",  d2 = "none",       lbl = "Cendakimab"),
        list(d1 = "dupilumab",   d2 = "budesonide", lbl = "Dupilumab + Bud")
      )
      results <- lapply(scenarios_list, function(sc) {
        out <- tryCatch(
          sim_scenario(mod = mod(),
                       drug1 = sc$d1, dose1 = if (sc$d1 == "budesonide") 1 else
                                               if (sc$d1 == "cendakimab") 160 else 300,
                       int1 = if (sc$d1 == "budesonide") 0.5 else
                              if (sc$d1 == "cendakimab") 1 else
                              if (sc$d1 == "mepolizumab") 28 else 14,
                       drug2 = sc$d2, dose2 = 1, int2 = 0.5,
                       dietary = input$dietary,
                       sim_days = input$sim_weeks * 7,
                       bl_IL13 = input$bl_IL13,
                       bl_EOS_ESO = input$bl_eos,
                       bl_FIBRO = input$bl_fibro,
                       bl_IGE = input$bl_IgE),
          error = function(e) NULL)
        if (!is.null(out)) out$scenario <- sc$lbl
        out
      })
      bind_rows(Filter(Negate(is.null), results))
    })
  })

  # ── VALUE BOXES ──
  output$vb_eos_esophagus <- renderValueBox({
    df <- sim_data()
    val <- round(tail(df$EOS_ESO, 1), 0)
    col <- if (val < 15) "green" else if (val < 50) "yellow" else "red"
    valueBox(paste0(val, " eos/hpf"), "Peak Tissue Eos (Final Week)",
             icon = icon("microscope"), color = col)
  })
  output$vb_dysphag <- renderValueBox({
    df <- sim_data()
    val <- round(tail(df$DYSPHAG_SCORE, 1), 1)
    col <- if (val < 3) "green" else if (val < 6) "yellow" else "red"
    valueBox(paste0(val, " / 10"), "Dysphagia Score (Final Week)",
             icon = icon("utensils"), color = col)
  })
  output$vb_erefs <- renderValueBox({
    df <- sim_data()
    val <- round(tail(df$EREFS_SCORE, 1), 1)
    col <- if (val < 5) "green" else if (val < 10) "yellow" else "red"
    valueBox(paste0(val, " / 18"), "EREFS Score (Final Week)",
             icon = icon("eye"), color = col)
  })

  # ── PK PLOTS ──
  output$pk_plot <- renderPlotly({
    df <- sim_data()
    pk_cols <- c("DUP_C", "MEPO_C", "CENDA_C")
    pk_labels <- c("Dupilumab (mg/L)", "Mepolizumab (mg/L)", "Cendakimab (mg/L)")
    p <- plot_ly()
    for (i in seq_along(pk_cols)) {
      if (max(df[[pk_cols[i]]], na.rm = TRUE) > 0.001) {
        p <- add_trace(p, x = df$week, y = df[[pk_cols[i]]],
                       name = pk_labels[i], type = "scatter", mode = "lines",
                       line = list(width = 2))
      }
    }
    p %>% layout(title = "Biologic Drug PK Profiles",
                 xaxis = list(title = "Week"),
                 yaxis = list(title = "Concentration (mg/L)"),
                 hovermode = "x unified")
  })

  output$bud_pk_plot <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      add_trace(x = df$week, y = df$BUD_ESO, name = "Esophageal",
                type = "scatter", mode = "lines",
                line = list(color = "#FF6D00", width = 2)) %>%
      add_trace(x = df$week, y = df$BUD_SYS, name = "Systemic",
                type = "scatter", mode = "lines",
                line = list(color = "#0288D1", width = 2)) %>%
      layout(xaxis = list(title = "Week"),
             yaxis = list(title = "Budesonide (mg/L)"),
             hovermode = "x unified")
  })

  output$cenda_pk_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(x = df$week, y = df$CENDA_C, type = "scatter", mode = "lines",
            line = list(color = "#00695C", width = 2)) %>%
      layout(xaxis = list(title = "Week"),
             yaxis = list(title = "Cendakimab (mg/L)"),
             hovermode = "x unified")
  })

  output$pk_table <- renderTable({
    data.frame(
      Parameter = c("Dupilumab t½", "Mepolizumab t½", "Cendakimab t½",
                    "Budesonide esoph t½", "DUP Css,trough", "MEPO Css,trough"),
      Value = c("~23 days", "~20 days", "~4.6 h",
                "~6 h (esoph)", "~65-90 µg/mL", "~30 µg/mL")
    )
  })

  # ── CYTOKINE PLOTS ──
  gg_to_plotly <- function(df, y_col, title, y_lab, color = "#E91E63") {
    plot_ly(x = df$week, y = df[[y_col]], type = "scatter", mode = "lines",
            line = list(color = color, width = 2.5)) %>%
      layout(title = title,
             xaxis = list(title = "Week"),
             yaxis = list(title = y_lab),
             hovermode = "x unified")
  }

  output$il13_plot <- renderPlotly({
    gg_to_plotly(sim_data(), "IL13", "IL-13 (pg/mL)", "IL-13 (pg/mL)", "#FF6D00")
  })
  output$il5_plot <- renderPlotly({
    gg_to_plotly(sim_data(), "IL5", "IL-5 (pg/mL)", "IL-5 (pg/mL)", "#6A1B9A")
  })
  output$eotax3_plot <- renderPlotly({
    gg_to_plotly(sim_data(), "EOTAX3", "Eotaxin-3/CCL26 (pg/mL)", "CCL26 (pg/mL)", "#1565C0")
  })
  output$ige_plot <- renderPlotly({
    gg_to_plotly(sim_data(), "IGE_TOT", "Total IgE (IU/mL)", "IgE (IU/mL)", "#37474F")
  })

  # ── EOSINOPHIL PLOTS ──
  output$eos_tissue_plot <- renderPlotly({
    df <- sim_data()
    plot_ly() %>%
      add_trace(x = df$week, y = df$EOS_ESO, name = "Tissue Eos",
                type = "scatter", mode = "lines",
                line = list(color = "#D32F2F", width = 2.5)) %>%
      add_segments(x = min(df$week), xend = max(df$week),
                   y = 15, yend = 15,
                   line = list(color = "gray", dash = "dash", width = 1),
                   name = "Remission threshold (<15 eos/hpf)") %>%
      layout(xaxis = list(title = "Week"),
             yaxis = list(title = "Tissue Eosinophils (eos/hpf)"),
             hovermode = "x unified")
  })

  output$eos_blood_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(x = df$week, y = df$EOS_BL, type = "scatter", mode = "lines",
            line = list(color = "#E57373", width = 2)) %>%
      layout(xaxis = list(title = "Week"),
             yaxis = list(title = "Blood Eosinophils (cells/µL)"))
  })

  output$mast_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(x = df$week, y = df$MAST_ESO, type = "scatter", mode = "lines",
            line = list(color = "#827717", width = 2)) %>%
      layout(xaxis = list(title = "Week"),
             yaxis = list(title = "Mast Cells (per mm²)"))
  })

  output$histo_pie <- renderPlotly({
    df <- sim_data()
    last_val <- tail(df$EOS_ESO, 1)
    remiss <- if (last_val < 15) "Remission (<15 eos/hpf)" else "Active (≥15 eos/hpf)"
    plot_ly(labels = c("Remission", "Active"),
            values = c(ifelse(last_val < 15, 1, 0), ifelse(last_val >= 15, 1, 0)),
            type = "pie",
            marker = list(colors = c("#4CAF50", "#F44336"))) %>%
      layout(showlegend = TRUE)
  })

  output$remission_time <- renderText({
    df <- sim_data()
    remiss_rows <- which(df$EOS_ESO < 15)
    if (length(remiss_rows) == 0) {
      "No histological remission achieved within simulation period"
    } else {
      first_wk <- df$week[remiss_rows[1]]
      sprintf("First remission: Week %.0f\n(EOS_ESO = %.1f eos/hpf)",
              first_wk, df$EOS_ESO[remiss_rows[1]])
    }
  })

  # ── CLINICAL ENDPOINT PLOTS ──
  output$dysphag_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(x = df$week, y = df$DYSPHAG_SCORE, type = "scatter", mode = "lines",
            line = list(color = "#1565C0", width = 2)) %>%
      layout(xaxis = list(title = "Week"),
             yaxis = list(title = "Dysphagia Score (0-10)", range = c(0, 10)))
  })

  output$erefs_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(x = df$week, y = df$EREFS_SCORE, type = "scatter", mode = "lines",
            line = list(color = "#7B1FA2", width = 2)) %>%
      layout(xaxis = list(title = "Week"),
             yaxis = list(title = "EREFS Score (0-18)", range = c(0, 18)))
  })

  output$epbar_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(x = df$week, y = df$EPBAR, type = "scatter", mode = "lines",
            line = list(color = "#2E7D32", width = 2)) %>%
      layout(xaxis = list(title = "Week"),
             yaxis = list(title = "Barrier Integrity (0-1)", range = c(0, 1)))
  })

  output$fibro_plot <- renderPlotly({
    df <- sim_data()
    plot_ly(x = df$week, y = df$FIBRO, type = "scatter", mode = "lines",
            line = list(color = "#37474F", width = 2)) %>%
      layout(xaxis = list(title = "Week"),
             yaxis = list(title = "LP Fibrosis Score (0-1)", range = c(0, 1)))
  })

  output$endpoints_table <- renderDT({
    df <- sim_data()
    wk12 <- df[which.min(abs(df$week - 12)), ]
    wk24 <- df[which.min(abs(df$week - 24)), ]
    tbl <- data.frame(
      Timepoint = c("Baseline", "Week 12", "Week 24"),
      `Eos/hpf` = c(input$bl_eos, round(wk12$EOS_ESO, 0), round(wk24$EOS_ESO, 0)),
      `Blood EOS` = c(round(600*input$bl_eos/80, 0), round(wk12$EOS_BL, 0), round(wk24$EOS_BL, 0)),
      `IL-13 (pg/mL)` = c(input$bl_IL13, round(wk12$IL13, 1), round(wk24$IL13, 1)),
      `Dysphagia` = c(round(10*(0.5*input$bl_eos/(input$bl_eos+60)+0.2*0.6), 1),
                      round(wk12$DYSPHAG_SCORE, 1), round(wk24$DYSPHAG_SCORE, 1)),
      EREFS = c(round(10*input$bl_eos/(input$bl_eos+50)*3, 1),
                round(wk12$EREFS_SCORE, 1), round(wk24$EREFS_SCORE, 1)),
      `Histologic Remission` = c("No", ifelse(wk12$EOS_ESO < 15, "YES", "No"),
                                  ifelse(wk24$EOS_ESO < 15, "YES", "No"))
    )
    datatable(tbl, options = list(dom = "t", pageLength = 5), rownames = FALSE)
  })

  # ── SCENARIO COMPARISON ──
  SCEN_COLORS <- c("No Treatment" = "#E53935",
                   "Budesonide ODT" = "#FF8F00",
                   "Dupilumab" = "#1565C0",
                   "Mepolizumab" = "#6A1B9A",
                   "Cendakimab" = "#00695C",
                   "Dupilumab + Bud" = "#37474F")

  output$compare_eos_plot <- renderPlotly({
    df <- compare_data()
    p <- plot_ly()
    for (sc in unique(df$scenario)) {
      d <- filter(df, scenario == sc)
      col <- SCEN_COLORS[sc]
      p <- add_trace(p, x = d$week, y = d$EOS_ESO,
                     name = sc, type = "scatter", mode = "lines",
                     line = list(color = col, width = 2))
    }
    p %>%
      add_segments(x = 0, xend = max(df$week), y = 15, yend = 15,
                   line = list(color = "gray", dash = "dot"),
                   name = "Remission (<15 eos/hpf)",
                   inherit = FALSE) %>%
      layout(xaxis = list(title = "Week"),
             yaxis = list(title = "Tissue Eosinophils (eos/hpf)"),
             hovermode = "x unified")
  })

  output$compare_dysphag <- renderPlotly({
    df <- compare_data()
    p <- plot_ly()
    for (sc in unique(df$scenario)) {
      d <- filter(df, scenario == sc)
      p <- add_trace(p, x = d$week, y = d$DYSPHAG_SCORE,
                     name = sc, type = "scatter", mode = "lines",
                     line = list(color = SCEN_COLORS[sc], width = 2))
    }
    p %>% layout(xaxis = list(title = "Week"),
                 yaxis = list(title = "Dysphagia Score (0-10)"))
  })

  output$compare_fibro <- renderPlotly({
    df <- compare_data()
    p <- plot_ly()
    for (sc in unique(df$scenario)) {
      d <- filter(df, scenario == sc)
      p <- add_trace(p, x = d$week, y = d$FIBRO,
                     name = sc, type = "scatter", mode = "lines",
                     line = list(color = SCEN_COLORS[sc], width = 2))
    }
    p %>% layout(xaxis = list(title = "Week"),
                 yaxis = list(title = "LP Fibrosis (0-1)"))
  })

  output$compare_bar <- renderPlotly({
    df <- compare_data()
    wk24 <- df %>%
      group_by(scenario) %>%
      filter(abs(week - 24) == min(abs(week - 24))) %>%
      slice(1) %>%
      ungroup()

    metrics <- c("EOS_ESO", "DYSPHAG_SCORE", "EREFS_SCORE", "FIBRO")
    labels <- c("Tissue Eos (eos/hpf)", "Dysphagia (0-10)", "EREFS (0-18)", "Fibrosis (0-1)")

    wk24_long <- wk24 %>%
      select(scenario, all_of(metrics)) %>%
      pivot_longer(cols = all_of(metrics), names_to = "metric", values_to = "value") %>%
      mutate(metric_label = factor(metric, levels = metrics, labels = labels))

    plot_ly(wk24_long, x = ~metric_label, y = ~value, color = ~scenario,
            colors = SCEN_COLORS,
            type = "bar") %>%
      layout(barmode = "group",
             xaxis = list(title = "Endpoint"),
             yaxis = list(title = "Value at Week 24"),
             title = "Head-to-Head Comparison at Week 24")
  })

  # ── BIOMARKER PANEL ──
  output$biomarker_heatmap <- renderPlotly({
    df <- compare_data()
    bl <- df %>% filter(week == 0) %>% group_by(scenario) %>% slice(1)
    wk24 <- df %>% group_by(scenario) %>%
      filter(abs(week - 24) == min(abs(week - 24))) %>% slice(1) %>% ungroup()

    markers <- c("EOS_ESO", "EOS_BL", "IL13", "IL5", "EOTAX3", "IGE_TOT", "MAST_ESO",
                 "EPBAR", "FIBRO")
    marker_labels <- c("Tissue Eos", "Blood Eos", "IL-13", "IL-5",
                        "Eotaxin-3", "Total IgE", "Mast Cells", "Barrier", "Fibrosis")

    pct_chg <- lapply(unique(df$scenario), function(sc) {
      bl_row <- filter(bl, scenario == sc)[1, ]
      wk_row <- filter(wk24, scenario == sc)[1, ]
      sapply(markers, function(m) {
        bl_val <- bl_row[[m]]
        if (is.null(bl_val) || is.na(bl_val) || bl_val == 0) return(0)
        (wk_row[[m]] - bl_val) / abs(bl_val) * 100
      })
    })
    names(pct_chg) <- unique(df$scenario)
    mat <- do.call(rbind, pct_chg)
    colnames(mat) <- marker_labels

    plot_ly(z = mat,
            x = marker_labels,
            y = rownames(mat),
            type = "heatmap",
            colorscale = list(c(0, "#1565C0"), c(0.5, "white"), c(1, "#D32F2F")),
            zmid = 0,
            text = round(mat, 1),
            texttemplate = "%{text}%",
            showscale = TRUE) %>%
      layout(title = "Biomarker % Change from Baseline at Week 24",
             xaxis = list(title = "Biomarker"),
             yaxis = list(title = "Treatment"))
  })

  output$biomarker_all <- renderPlotly({
    df <- sim_data()
    markers <- list(
      list(col = "EOS_ESO", label = "Tissue Eos", color = "#D32F2F"),
      list(col = "IL13", label = "IL-13", color = "#FF6D00"),
      list(col = "EOTAX3", label = "Eotaxin-3 (÷10)", color = "#1565C0"),
      list(col = "EPBAR", label = "Barrier (×50)", color = "#2E7D32"),
      list(col = "FIBRO", label = "Fibrosis (×100)", color = "#37474F")
    )
    p <- plot_ly()
    for (m in markers) {
      scale <- if (m$col == "EOTAX3") 0.1 else if (m$col %in% c("EPBAR", "FIBRO")) 50 else 1
      p <- add_trace(p, x = df$week,
                     y = df[[m$col]] * (if (m$col == "EPBAR") 50 else
                                         if (m$col == "FIBRO") 100 else
                                         if (m$col == "EOTAX3") 0.1 else 1),
                     name = m$label, type = "scatter", mode = "lines",
                     line = list(color = m$color, width = 2))
    }
    p %>% layout(xaxis = list(title = "Week"),
                 yaxis = list(title = "Scaled Value"),
                 hovermode = "x unified",
                 title = "Scaled Biomarker Trajectories")
  })

  output$biomarker_corr <- renderPlotly({
    df <- sim_data()
    vars <- c("EOS_ESO", "EOS_BL", "IL13", "IL5", "EOTAX3",
              "IGE_TOT", "MAST_ESO", "EPBAR", "FIBRO")
    lbls <- c("Tissue Eos", "Blood Eos", "IL-13", "IL-5",
               "Eotaxin-3", "IgE", "Mast Cells", "Barrier", "Fibrosis")
    mat <- cor(df[, vars], use = "complete.obs")
    rownames(mat) <- colnames(mat) <- lbls

    plot_ly(z = mat, x = lbls, y = lbls,
            type = "heatmap",
            colorscale = "RdBu",
            zmid = 0,
            text = round(mat, 2),
            texttemplate = "%{text}") %>%
      layout(title = "Biomarker Correlation Matrix")
  })
}

# ============================================================
# LAUNCH
# ============================================================

shinyApp(ui = ui, server = server)
