## ============================================================
## Cystic Fibrosis (CF) QSP — Interactive Shiny Dashboard
## ============================================================
## Tabs:
##   1. Patient Profile & Mutation Class
##   2. CFTR Modulator PK Profiles
##   3. CFTR Function & Biomarkers
##   4. Lung Function (ppFEV1 / Exacerbations)
##   5. Treatment Scenario Comparison
##   6. ASL, Infection & Inflammation
## ============================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(DT)
library(plotly)
library(mrgsolve)

# ── Embed the mrgsolve model code ───────────────────────────
cf_code <- '
$PARAM
IVA_dose_flag = 0; LUM_dose_flag = 0; TEZ_dose_flag = 0; ELX_dose_flag = 0
TOBRA_flag = 0; mutation_type = 2
IVA_F = 0.67; IVA_Ka = 0.75; IVA_Vc = 97.1; IVA_Vp = 255.9
IVA_CL = 17.3; IVA_Q = 5.2; IVA_EC50 = 0.1; IVA_Emax = 0.85; IVA_hill = 1.5
ELX_F = 0.80; ELX_Ka = 0.56; ELX_Vc = 193.0; ELX_CL = 4.1; ELX_EC50 = 0.05; ELX_Emax = 0.55
TEZ_F = 0.70; TEZ_Ka = 0.65; TEZ_Vc = 271.0; TEZ_CL = 13.5; TEZ_EC50 = 0.08; TEZ_Emax = 0.25
CFTR_synth = 1.0; CFTR_deg_WT = 0.05; CFTR_fold_EF = 0.01; CFTR_fold_WT = 0.75
CFTR_deg_B = 0.5; CFTR_traffic = 0.3; CFTR_endo = 0.15; CFTR_Po_WT = 0.45
CFTR_Po_EF = 0.04; CFTR_Po_G551 = 0.02
ASL_normal = 7.0; ASL_min = 1.0; kASL_secrete = 0.05; kASL_absorb = 0.08
kASL_restore = 0.02; ASL_MCC_thresh = 4.0
Pa_growth = 0.35; Pa_Kmax = 1e8; Pa_kill_host = 0.15; Pa_biofilm_k = 0.08
Pa_biofilm_d = 0.02; Pa_biofilm_prot = 0.05; TOBRA_Emax = 0.90; TOBRA_EC50 = 2.0
IL8_baseline = 100.0; IL8_bacteria = 0.5; IL8_deg = 0.2
Neu_recruit = 0.8; Neu_halflife = 12.0; Damage_rate = 0.003; Damage_repair = 0.001
FEV1_0 = 90.0; FEV1_decline_noRx = 1.5; FEV1_decline_ETI = 0.2
kFEV1_damage = 0.8; kFEV1_recover = 0.1; FEV1_exac_loss = 1.2
Exac_rate_base = 1.8; Exac_bacteria_k = 3.0; Exac_ASL_k = 2.5
Pancreas_0 = 15.0; BMI_0 = 18.5; BMI_target = 22.0; kBMI_ETI = 0.005

$CMT IVA_gut IVA_central IVA_periph ELX_gut ELX_central TEZ_gut TEZ_central
     CFTR_bandB CFTR_mem ASL Pa_free Pa_film IL8 Neutrophil Damage FEV1 Exac_cumul BMI_state Panc_fn

$MAIN
double fold_eff; double gating_eff;
if (mutation_type == 1) { fold_eff = 0.0; gating_eff = 0.0; }
else if (mutation_type == 2) { fold_eff = CFTR_fold_EF; gating_eff = CFTR_Po_EF; }
else if (mutation_type == 3) { fold_eff = 0.60; gating_eff = CFTR_Po_G551; }
else { fold_eff = 0.25; gating_eff = 0.20; }
IVA_gut_0 = 0; IVA_central_0 = 0; IVA_periph_0 = 0;
ELX_gut_0 = 0; ELX_central_0 = 0; TEZ_gut_0 = 0; TEZ_central_0 = 0;
CFTR_bandB_0 = fold_eff * 100.0; CFTR_mem_0 = fold_eff * 40.0;
ASL_0 = ASL_normal * (0.3 + 0.7 * fold_eff / CFTR_fold_WT);
Pa_free_0 = 0.1; Pa_film_0 = 0.01;
IL8_0 = IL8_baseline; Neutrophil_0 = 0.5; Damage_0 = 0.0;
FEV1_0 = FEV1_0; Exac_cumul_0 = 0.0;
BMI_state_0 = BMI_0; Panc_fn_0 = (mutation_type >= 3) ? 70.0 : 15.0;

$ODE
dxdt_IVA_gut = -IVA_Ka * IVA_gut;
dxdt_IVA_central = IVA_Ka * IVA_gut - (IVA_CL/IVA_Vc)*IVA_central - (IVA_Q/IVA_Vc)*IVA_central + (IVA_Q/IVA_Vp)*IVA_periph;
dxdt_IVA_periph = (IVA_Q/IVA_Vc)*IVA_central - (IVA_Q/IVA_Vp)*IVA_periph;
dxdt_ELX_gut = -ELX_Ka * ELX_gut;
dxdt_ELX_central = ELX_Ka * ELX_gut - (ELX_CL/ELX_Vc)*ELX_central;
dxdt_TEZ_gut = -TEZ_Ka * TEZ_gut;
dxdt_TEZ_central = TEZ_Ka * TEZ_gut - (TEZ_CL/TEZ_Vc)*TEZ_central;
double IVA_Cp = IVA_central / IVA_Vc;
double ELX_Cp = ELX_central / ELX_Vc;
double TEZ_Cp = TEZ_central / TEZ_Vc;
double fold_eff2; double gating_eff2;
if (mutation_type == 1) { fold_eff2 = 0.0; gating_eff2 = 0.0; }
else if (mutation_type == 2) { fold_eff2 = CFTR_fold_EF; gating_eff2 = CFTR_Po_EF; }
else if (mutation_type == 3) { fold_eff2 = 0.60; gating_eff2 = CFTR_Po_G551; }
else { fold_eff2 = 0.25; gating_eff2 = 0.20; }
double ELX_effect = (ELX_Cp > 0) ? ELX_Emax * pow(ELX_Cp, IVA_hill) / (pow(ELX_EC50, IVA_hill) + pow(ELX_Cp, IVA_hill)) : 0.0;
double TEZ_effect = (TEZ_Cp > 0) ? TEZ_Emax * TEZ_Cp / (TEZ_EC50 + TEZ_Cp) : 0.0;
double LUM_effect = LUM_dose_flag * 0.10;
double correction = 1.0 - (1.0 - ELX_effect)*(1.0 - TEZ_effect)*(1.0 - LUM_effect);
double potentiation = (IVA_Cp > 0) ? IVA_Emax * pow(IVA_Cp, IVA_hill) / (pow(IVA_EC50, IVA_hill) + pow(IVA_Cp, IVA_hill)) : 0.0;
double fold_eff_c = fold_eff2 + correction * (CFTR_fold_WT - fold_eff2);
double gating_c   = gating_eff2 + potentiation * (CFTR_Po_WT - gating_eff2);
dxdt_CFTR_bandB = CFTR_synth * 10.0 - (CFTR_traffic + correction*0.2)*CFTR_bandB - CFTR_deg_B*(1.0-correction*0.7)*CFTR_bandB;
dxdt_CFTR_mem   = (CFTR_traffic + correction*0.2)*CFTR_bandB - CFTR_deg_WT*CFTR_mem - CFTR_endo*CFTR_mem;
double CFTR_fn_v = (fold_eff_c * gating_c) / (CFTR_fold_WT * CFTR_Po_WT) * 100.0;
if (CFTR_fn_v > 100.0) CFTR_fn_v = 100.0;
double ASL_sec = kASL_secrete * (CFTR_fn_v / 100.0) * (ASL_normal - ASL);
double ASL_abs = kASL_absorb * (1.0 - (CFTR_fn_v/100.0)*0.5);
dxdt_ASL = ASL_sec - ASL_abs * ASL + kASL_restore*(ASL_normal - ASL);
if (ASL < ASL_min) dxdt_ASL = 0.0;
double Pa_tot = Pa_free + Pa_film;
double ASL_kf = (ASL > ASL_MCC_thresh) ? 1.0 : ASL / ASL_MCC_thresh;
double host_k = Pa_kill_host * Neutrophil * ASL_kf * Pa_free;
double TOBRA_Cp2 = TOBRA_flag * 100.0;
double TOBRA_k = TOBRA_Emax * TOBRA_Cp2 / (TOBRA_EC50 + TOBRA_Cp2) * Pa_free;
dxdt_Pa_free = Pa_growth * Pa_free * (1.0 - Pa_free/Pa_Kmax) - host_k - TOBRA_k - Pa_biofilm_k * Pa_free + Pa_biofilm_d * Pa_film;
dxdt_Pa_film = Pa_biofilm_k * Pa_free - Pa_biofilm_d * Pa_film - TOBRA_k * Pa_biofilm_prot;
double Pa_sig = log1p(Pa_tot);
dxdt_IL8       = IL8_bacteria * Pa_sig * IL8_baseline - IL8_deg * IL8;
dxdt_Neutrophil= Neu_recruit * IL8 / (IL8_baseline*5.0) - Neutrophil / Neu_halflife;
dxdt_Damage    = Damage_rate * Neutrophil * (1.0 + Pa_sig/3.0) - Damage_repair * Damage;
double FEV1_nat = FEV1_decline_noRx / 8760.0;
double FEV1_ben = (ELX_dose_flag > 0.5) ? (FEV1_decline_noRx - FEV1_decline_ETI)/8760.0 : 0.0;
dxdt_FEV1 = -FEV1_nat*FEV1 + FEV1_ben*FEV1 - kFEV1_damage*dxdt_Damage + kFEV1_recover*(CFTR_fn_v/100.0)*(95.0 - FEV1);
if (FEV1 < 20.0) dxdt_FEV1 = 0.0;
double exac_h = (Exac_rate_base/8760.0) * (1.0 + Exac_bacteria_k*log1p(Pa_tot)) * (1.0 + Exac_ASL_k*fmax(0.0,ASL_MCC_thresh - ASL)/ASL_MCC_thresh) * fmax(0.0, 1.0 - CFTR_fn_v/150.0);
dxdt_Exac_cumul = exac_h;
dxdt_BMI_state  = kBMI_ETI * ELX_dose_flag * (BMI_target - BMI_state) * 24.0;
dxdt_Panc_fn    = 0.0;

$TABLE
capture IVA_Cp = IVA_central / IVA_Vc;
capture ELX_Cp = ELX_central / ELX_Vc;
capture TEZ_Cp = TEZ_central / TEZ_Vc;
double fold2; double gate2;
if (mutation_type == 1) { fold2 = 0.0; gate2 = 0.0; }
else if (mutation_type == 2) { fold2 = CFTR_fold_EF; gate2 = CFTR_Po_EF; }
else if (mutation_type == 3) { fold2 = 0.60; gate2 = CFTR_Po_G551; }
else { fold2 = 0.25; gate2 = 0.20; }
double ELXe = (ELX_central/ELX_Vc > 0) ? ELX_Emax * pow(ELX_central/ELX_Vc, IVA_hill)/(pow(ELX_EC50,IVA_hill)+pow(ELX_central/ELX_Vc,IVA_hill)) : 0.0;
double TEZe = (TEZ_central/TEZ_Vc > 0) ? TEZ_Emax*(TEZ_central/TEZ_Vc)/(TEZ_EC50+TEZ_central/TEZ_Vc) : 0.0;
double corr2 = 1.0 - (1.0-ELXe)*(1.0-TEZe)*(1.0-LUM_dose_flag*0.10);
double pote2 = (IVA_central/IVA_Vc > 0) ? IVA_Emax*pow(IVA_central/IVA_Vc,IVA_hill)/(pow(IVA_EC50,IVA_hill)+pow(IVA_central/IVA_Vc,IVA_hill)) : 0.0;
double fce2 = fold2 + corr2*(CFTR_fold_WT - fold2);
double gce2 = gate2 + pote2*(CFTR_Po_WT - gate2);
capture CFTR_fn = fmin(100.0, (fce2*gce2)/(CFTR_fold_WT*CFTR_Po_WT)*100.0);
capture sweat_Cl = 110.0 - 70.0*(CFTR_fn/100.0);
capture ppFEV1 = FEV1;
capture exac_total = Exac_cumul;
capture ASL_ht = ASL;
capture logPa = log10(fmax(Pa_free, 1e-6));
capture BMI_val = BMI_state;
capture IL8_val = IL8;
capture Neu_val = Neutrophil;
capture Dam_val = Damage;

$CAPTURE IVA_Cp ELX_Cp TEZ_Cp CFTR_fn sweat_Cl ppFEV1 exac_total ASL_ht logPa BMI_val IL8_val Neu_val Dam_val
'

# Load model
cf_mod <- mrgsolve::mcode("CF_shiny", cf_code)

# Helper to run simulation
run_cf_sim <- function(mod, params, events = NULL, end_h = 8760, delta_h = 24) {
  m2 <- param(mod, params)
  if (!is.null(events)) {
    out <- mrgsim(m2, events = events, end = end_h, delta = delta_h)
  } else {
    out <- mrgsim(m2, end = end_h, delta = delta_h)
  }
  df <- as.data.frame(out)
  df$time_days <- df$time / 24
  return(df)
}

make_ev_shiny <- function(give_IVA, give_TEZ, give_ELX, give_LUM, end_h) {
  evs <- list()
  if (give_IVA) evs[["IVA"]] <- ev(amt = 150*0.67, cmt = "IVA_gut", time = 0, ii = 12, addl = ceiling(end_h/12))
  if (give_TEZ) evs[["TEZ"]] <- ev(amt = 100*0.70, cmt = "TEZ_gut", time = 0, ii = 24, addl = ceiling(end_h/24))
  if (give_ELX) evs[["ELX"]] <- ev(amt = 200*0.80, cmt = "ELX_gut", time = 0, ii = 24, addl = ceiling(end_h/24))
  if (give_LUM) evs[["LUM"]] <- ev(amt = 200*0.65, cmt = "IVA_gut", time = 0, ii = 12, addl = ceiling(end_h/12))  # simplified
  if (length(evs) > 0) return(do.call(c, evs)) else return(NULL)
}

# ── Mutation info table ──────────────────────────────────────
mutation_info <- data.frame(
  Class    = c("I","II","III","IV","V","VI"),
  Defect   = c("No protein (NMD/PTC)","Processing (ERAD)",
               "Gating defect","Conductance defect",
               "Reduced quantity","Membrane stability"),
  Examples = c("W1282X, G542X, R553X", "ΔF508 (~70% CF alleles)",
               "G551D (~4%)", "R117H",
               "3849+10kbC→T","120del23"),
  Drug     = c("Ataluren (stop readthrough; investigational)",
               "ELX/TEZ correctors + IVA potentiator (Trikafta)",
               "IVA potentiator alone (Kalydeco)",
               "IVA potentiator ± corrector",
               "TEZ/IVA (Symdeko) or ETI",
               "TEZ/IVA or ETI"),
  stringsAsFactors = FALSE
)

# ── Clinical trial benchmark table ──────────────────────────
trial_bench <- data.frame(
  Trial    = c("STRIVE (IVA)","EVOLENT (TEZ/IVA)","VX-445-102","VX-445-103 (AURORA)"),
  Drug     = c("Ivacaftor 150mg q12h","TEZ/IVA","ELX/TEZ/IVA","ELX/TEZ/IVA"),
  Mutation = c("G551D/WT","ΔF508/ΔF508","ΔF508/MF","ΔF508/ΔF508"),
  N        = c(144, 504, 403, 107),
  ppFEV1   = c("+10.6 pp", "+3.4 pp", "+14.3 pp", "+13.8 pp"),
  SweatCl  = c("-48 mmol/L", "-9.7 mmol/L", "-41.8 mmol/L", "-45.1 mmol/L"),
  ExacRR   = c("-55%", "-35%", "-63%", "-60%"),
  CFQR     = c("+8.1", "+4.0", "+17.4", "+20.2"),
  stringsAsFactors = FALSE
)

# ════════════════════════════════════════════════════════════
# SHINY UI
# ════════════════════════════════════════════════════════════
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "CF QSP Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",    tabName = "tab_patient",   icon = icon("user")),
      menuItem("CFTR Modulator PK",  tabName = "tab_pk",        icon = icon("pills")),
      menuItem("CFTR Function",      tabName = "tab_cftr",      icon = icon("dna")),
      menuItem("Lung Function",      tabName = "tab_lung",      icon = icon("lungs")),
      menuItem("Scenario Comparison",tabName = "tab_scenario",  icon = icon("chart-bar")),
      menuItem("ASL & Infection",    tabName = "tab_asl",       icon = icon("bacteria"))
    )
  ),
  dashboardBody(
    tabItems(
      # ──────────────────────────────────────────────────────
      # TAB 1: Patient Profile
      # ──────────────────────────────────────────────────────
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient & Mutation Configuration", status = "primary",
              solidHeader = TRUE, width = 4,
              radioButtons("mutation_type", "CFTR Mutation Class:",
                choices = list(
                  "Class I — No protein (W1282X)" = 1,
                  "Class II — ΔF508/ΔF508 (processing)" = 2,
                  "Class III — G551D (gating defect)" = 3,
                  "Class IV/Mixed — R117H / compound het." = 4
                ), selected = 2),
              sliderInput("age_start", "Age at Treatment Start (years):", 2, 40, 18),
              sliderInput("FEV1_base", "Baseline ppFEV1 (%):", 20, 110, 90),
              sliderInput("sim_years", "Simulation Duration (years):", 1, 10, 2),
              sliderInput("exac_rate_base", "Baseline Exacerbation Rate/yr:", 0.5, 5, 1.8, step = 0.1)
          ),
          box(title = "CFTR Mutation Classes — Summary", status = "info",
              solidHeader = TRUE, width = 8,
              DTOutput("mutation_table"),
              br(),
              p("Epidemiology: ~90% of CF patients carry at least one ΔF508 allele (Class II). ",
                "Approximately 70% are ΔF508 homozygous."),
              p("Trikafta (ELX/TEZ/IVA) is approved for patients aged ≥2 years with at least one ",
                "ΔF508 allele OR one of 177 responsive mutations (~90% of CF patients).")
          )
        ),
        fluidRow(
          box(title = "Clinical Trial Benchmarks", status = "success",
              solidHeader = TRUE, width = 12,
              DTOutput("trial_table")
          )
        )
      ),
      # ──────────────────────────────────────────────────────
      # TAB 2: CFTR Modulator PK
      # ──────────────────────────────────────────────────────
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "PK Simulation Settings", status = "primary",
              solidHeader = TRUE, width = 3,
              checkboxInput("pk_IVA", "Ivacaftor 150mg q12h", TRUE),
              checkboxInput("pk_ELX", "Elexacaftor 200mg q24h", TRUE),
              checkboxInput("pk_TEZ", "Tezacaftor 100mg q24h", TRUE),
              sliderInput("pk_days", "Simulation Days:", 1, 14, 7),
              hr(),
              h5("Key PK Parameters:"),
              tableOutput("pk_params_table")
          ),
          box(title = "Plasma Concentration Profiles", status = "primary",
              solidHeader = TRUE, width = 9,
              plotlyOutput("pk_plot", height = "400px"),
              br(),
              verbatimTextOutput("pk_summary")
          )
        )
      ),
      # ──────────────────────────────────────────────────────
      # TAB 3: CFTR Function & Biomarkers
      # ──────────────────────────────────────────────────────
      tabItem(tabName = "tab_cftr",
        fluidRow(
          box(title = "Treatment Selection", status = "primary",
              solidHeader = TRUE, width = 3,
              checkboxInput("cftr_IVA", "Ivacaftor", TRUE),
              checkboxInput("cftr_LUM", "Lumacaftor (Orkambi)", FALSE),
              checkboxInput("cftr_TEZ", "Tezacaftor", TRUE),
              checkboxInput("cftr_ELX", "Elexacaftor (ETI)", TRUE),
              hr(),
              h5("Expected CFTR Function:"),
              tableOutput("cftr_expected_table")
          ),
          box(title = "CFTR Function (% WT) over Time", status = "success",
              solidHeader = TRUE, width = 9,
              plotlyOutput("cftr_fn_plot", height = "280px"),
              br(),
              plotlyOutput("sweat_cl_plot", height = "250px")
          )
        ),
        fluidRow(
          valueBoxOutput("cftr_fn_box", width = 3),
          valueBoxOutput("sweat_cl_box", width = 3),
          valueBoxOutput("correction_box", width = 3),
          valueBoxOutput("potentiation_box", width = 3)
        )
      ),
      # ──────────────────────────────────────────────────────
      # TAB 4: Lung Function
      # ──────────────────────────────────────────────────────
      tabItem(tabName = "tab_lung",
        fluidRow(
          box(title = "Lung Function Simulation", status = "primary",
              solidHeader = TRUE, width = 3,
              checkboxGroupInput("lung_scenarios", "Scenarios to show:",
                choices = list(
                  "Untreated" = "none",
                  "Ivacaftor (G551D)" = "iva",
                  "LUM/IVA (Orkambi)" = "lum",
                  "TEZ/IVA (Symdeko)" = "tez",
                  "ETI (Trikafta)" = "eti"
                ),
                selected = c("none", "eti")
              ),
              sliderInput("lung_years", "Duration (years):", 1, 10, 5)
          ),
          box(title = "ppFEV1 Trajectory", status = "success",
              solidHeader = TRUE, width = 9,
              plotlyOutput("fev1_plot", height = "350px"),
              br(),
              plotlyOutput("exac_plot", height = "220px")
          )
        ),
        fluidRow(
          box(title = "FEV1 Severity Classification", status = "info",
              solidHeader = TRUE, width = 12,
              fluidRow(
                valueBoxOutput("fev1_mild_box", width = 3),
                valueBoxOutput("fev1_mod_box", width = 3),
                valueBoxOutput("fev1_sev_box", width = 3),
                valueBoxOutput("fev1_vsev_box", width = 3)
              )
          )
        )
      ),
      # ──────────────────────────────────────────────────────
      # TAB 5: Scenario Comparison
      # ──────────────────────────────────────────────────────
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "Simulation Settings", status = "primary",
              solidHeader = TRUE, width = 3,
              sliderInput("sc_years", "Simulation Duration (years):", 1, 5, 2),
              radioButtons("sc_endpoint", "Primary Endpoint:",
                choices = list(
                  "ppFEV1 (% predicted)" = "ppFEV1",
                  "Sweat Chloride (mmol/L)" = "sweat_Cl",
                  "CFTR Function (%)" = "CFTR_fn",
                  "Annual Exacerbations" = "exac_total",
                  "ASL Height (μm)" = "ASL_ht",
                  "BMI (kg/m²)" = "BMI_val"
                ), selected = "ppFEV1"),
              checkboxInput("sc_tobra", "Add Tobramycin (Scenario 6)", FALSE)
          ),
          box(title = "Multi-scenario Time Course", status = "success",
              solidHeader = TRUE, width = 9,
              plotlyOutput("scenario_plot", height = "400px")
          )
        ),
        fluidRow(
          box(title = "52-Week Endpoint Summary Table", status = "warning",
              solidHeader = TRUE, width = 12,
              DTOutput("scenario_summary_table")
          )
        )
      ),
      # ──────────────────────────────────────────────────────
      # TAB 6: ASL, Infection & Inflammation
      # ──────────────────────────────────────────────────────
      tabItem(tabName = "tab_asl",
        fluidRow(
          box(title = "ASL & Infection Settings", status = "primary",
              solidHeader = TRUE, width = 3,
              checkboxInput("asl_eti",  "Trikafta (ETI)", TRUE),
              checkboxInput("asl_tobra","Inhaled Tobramycin", FALSE),
              checkboxInput("asl_dnase","Dornase alfa (DNase)", FALSE),
              hr(),
              h5("Mucociliary Clearance:"),
              p("ASL must exceed 4 μm for effective MCC."),
              p("ΔF508 CF: ASL typically collapses to 1–3 μm."),
              p("ETI restores ASL toward 5–7 μm range.")
          ),
          box(title = "ASL Height & MCC Status", status = "info",
              solidHeader = TRUE, width = 9,
              plotlyOutput("asl_plot", height = "280px"),
              br(),
              plotlyOutput("bacteria_plot", height = "220px")
          )
        ),
        fluidRow(
          box(title = "Inflammatory Markers", status = "danger",
              solidHeader = TRUE, width = 12,
              plotlyOutput("inflam_plot", height = "300px")
          )
        )
      )
    )
  )
)

# ════════════════════════════════════════════════════════════
# SHINY SERVER
# ════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  # ── Static tables ────────────────────────────────────────
  output$mutation_table <- renderDT({
    datatable(mutation_info, options = list(pageLength = 6, dom = "t"),
              rownames = FALSE) |>
      formatStyle("Class", fontWeight = "bold")
  })

  output$trial_table <- renderDT({
    datatable(trial_bench, options = list(pageLength = 5, dom = "t"),
              rownames = FALSE) |>
      formatStyle("ppFEV1", color = "green", fontWeight = "bold")
  })

  output$pk_params_table <- renderTable({
    data.frame(
      Drug = c("Ivacaftor","Elexacaftor","Tezacaftor"),
      F    = c("67%","80%","70%"),
      Vc   = c("97 L","193 L","271 L"),
      CL   = c("17.3 L/h","4.1 L/h","13.5 L/h"),
      t_half = c("12 h","27 h","14 h")
    )
  })

  # ── PK Simulation ────────────────────────────────────────
  pk_sim <- reactive({
    end_h <- input$pk_days * 24
    evs <- make_ev_shiny(input$pk_IVA, input$pk_TEZ, input$pk_ELX, FALSE, end_h)
    params <- list(mutation_type = 2, IVA_dose_flag = as.numeric(input$pk_IVA),
                   TEZ_dose_flag = as.numeric(input$pk_TEZ),
                   ELX_dose_flag = as.numeric(input$pk_ELX),
                   LUM_dose_flag = 0, TOBRA_flag = 0)
    run_cf_sim(cf_mod, params, evs, end_h = end_h, delta_h = 0.5)
  })

  output$pk_plot <- renderPlotly({
    df <- pk_sim()
    p <- ggplot(df, aes(x = time / 24)) +
      geom_line(aes(y = IVA_Cp, color = "Ivacaftor"), size = 1) +
      geom_line(aes(y = ELX_Cp, color = "Elexacaftor"), size = 1) +
      geom_line(aes(y = TEZ_Cp, color = "Tezacaftor"), size = 1) +
      scale_color_manual(values = c("Ivacaftor" = "#FF7F00",
                                    "Elexacaftor" = "#00B050",
                                    "Tezacaftor" = "#7030A0")) +
      labs(x = "Time (days)", y = "Plasma Concentration (μg/mL)", color = "Drug") +
      theme_bw()
    ggplotly(p)
  })

  output$pk_summary <- renderText({
    df <- pk_sim()
    paste0(
      "IVA: Cmax = ", round(max(df$IVA_Cp), 3), " μg/mL | ",
      "ELX: Cmax = ", round(max(df$ELX_Cp), 3), " μg/mL | ",
      "TEZ: Cmax = ", round(max(df$TEZ_Cp), 3), " μg/mL\n",
      "IVA Tmax = ", round(df$time[which.max(df$IVA_Cp)], 1), " h"
    )
  })

  # ── CFTR Function Simulation ──────────────────────────────
  cftr_sim <- reactive({
    params <- list(mutation_type = as.integer(input$mutation_type),
                   IVA_dose_flag = as.numeric(input$cftr_IVA),
                   LUM_dose_flag = as.numeric(input$cftr_LUM),
                   TEZ_dose_flag = as.numeric(input$cftr_TEZ),
                   ELX_dose_flag = as.numeric(input$cftr_ELX),
                   TOBRA_flag = 0, FEV1_0 = input$FEV1_base,
                   Exac_rate_base = input$exac_rate_base)
    evs <- make_ev_shiny(input$cftr_IVA, input$cftr_TEZ, input$cftr_ELX,
                         input$cftr_LUM, 8760)
    run_cf_sim(cf_mod, params, evs)
  })

  output$cftr_fn_plot <- renderPlotly({
    df <- cftr_sim()
    p <- ggplot(df, aes(x = time_days, y = CFTR_fn)) +
      geom_line(color = "#00B050", size = 1.3) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
      annotate("text", x = 5, y = 11.5, label = "~10% threshold (minimal disease)",
               hjust = 0, size = 3, color = "red") +
      labs(x = "Days", y = "CFTR Function (% of WT)") +
      theme_bw()
    ggplotly(p)
  })

  output$sweat_cl_plot <- renderPlotly({
    df <- cftr_sim()
    p <- ggplot(df, aes(x = time_days, y = sweat_Cl)) +
      geom_line(color = "#0070C0", size = 1.3) +
      geom_hline(yintercept = 60, linetype = "dashed", color = "red") +
      geom_hline(yintercept = 30, linetype = "dashed", color = "green3") +
      labs(x = "Days", y = "Sweat Chloride (mmol/L)") +
      theme_bw()
    ggplotly(p)
  })

  output$cftr_fn_box <- renderValueBox({
    df <- cftr_sim()
    val <- round(mean(tail(df$CFTR_fn, 5)), 1)
    valueBox(paste0(val, "%"), "CFTR Function", icon = icon("dna"),
             color = if (val > 30) "green" else if (val > 10) "yellow" else "red")
  })

  output$sweat_cl_box <- renderValueBox({
    df <- cftr_sim()
    val <- round(mean(tail(df$sweat_Cl, 5)), 1)
    valueBox(paste0(val, " mmol/L"), "Sweat Chloride", icon = icon("droplet"),
             color = if (val < 30) "green" else if (val < 60) "yellow" else "red")
  })

  output$correction_box <- renderValueBox({
    df <- cftr_sim()
    mtype <- as.integer(input$mutation_type)
    fold_base <- if (mtype == 1) 0 else if (mtype == 2) 0.01 else 0.60
    correction_eff <- round((mean(tail(df$CFTR_fn,5))/100 * 0.75*0.45 - fold_base*0.04)/(0.75*0.45 - fold_base*0.04)*100, 1)
    valueBox(paste0(max(0, correction_eff), "%"), "Corrector Effect", icon = icon("wrench"),
             color = "blue")
  })

  output$potentiation_box <- renderValueBox({
    df <- cftr_sim()
    pot_val <- if (input$cftr_IVA) round(mean(tail(df$CFTR_fn,5)) * 0.6, 1) else 0
    valueBox(paste0(max(0, pot_val), "%"), "Potentiator Effect", icon = icon("bolt"),
             color = "purple")
  })

  output$cftr_expected_table <- renderTable({
    data.frame(
      Regimen = c("Untreated ΔF508","IVA alone","LUM/IVA","TEZ/IVA","ETI (Trikafta)"),
      CFTR_fn = c("<1%","~1%","~10–15%","~15–25%","~40–60%"),
      Sweat_Cl = c("~110","~100","~95","~85","~65–70")
    )
  })

  # ── Lung Function Simulation ──────────────────────────────
  lung_data <- reactive({
    end_h <- input$lung_years * 8760
    sc_map <- list(
      none = list(IVA = FALSE, LUM = FALSE, TEZ = FALSE, ELX = FALSE, mut = 2, col = "#CC0000", lbl = "Untreated"),
      iva  = list(IVA = TRUE,  LUM = FALSE, TEZ = FALSE, ELX = FALSE, mut = 3, col = "#0070C0", lbl = "IVA (G551D)"),
      lum  = list(IVA = TRUE,  LUM = TRUE,  TEZ = FALSE, ELX = FALSE, mut = 2, col = "#FF7F00", lbl = "LUM/IVA"),
      tez  = list(IVA = TRUE,  LUM = FALSE, TEZ = TRUE,  ELX = FALSE, mut = 2, col = "#7030A0", lbl = "TEZ/IVA"),
      eti  = list(IVA = TRUE,  LUM = FALSE, TEZ = TRUE,  ELX = TRUE,  mut = 2, col = "#00B050", lbl = "ETI/Trikafta")
    )
    results <- lapply(input$lung_scenarios, function(sc_name) {
      sc <- sc_map[[sc_name]]
      evs <- make_ev_shiny(sc$IVA, sc$TEZ, sc$ELX, sc$LUM, end_h)
      params <- list(mutation_type = sc$mut, IVA_dose_flag = as.numeric(sc$IVA),
                     LUM_dose_flag = as.numeric(sc$LUM), TEZ_dose_flag = as.numeric(sc$TEZ),
                     ELX_dose_flag = as.numeric(sc$ELX), TOBRA_flag = 0,
                     FEV1_0 = input$FEV1_base, Exac_rate_base = input$exac_rate_base)
      df <- run_cf_sim(cf_mod, params, evs, end_h, 24)
      df$label <- sc$lbl
      df$color <- sc$col
      df
    })
    do.call(rbind, results)
  })

  output$fev1_plot <- renderPlotly({
    df <- lung_data()
    p <- ggplot(df, aes(x = time_days, y = ppFEV1, color = label)) +
      geom_line(size = 1.2) +
      geom_hline(yintercept = c(70, 40), linetype = "dashed", color = c("orange","red")) +
      labs(x = "Days", y = "ppFEV1 (% predicted)", color = "Scenario") +
      theme_bw()
    ggplotly(p)
  })

  output$exac_plot <- renderPlotly({
    df <- lung_data()
    p <- ggplot(df, aes(x = time_days, y = exac_total, color = label)) +
      geom_line(size = 1.2) +
      labs(x = "Days", y = "Cumulative Exacerbations", color = "Scenario") +
      theme_bw()
    ggplotly(p)
  })

  output$fev1_mild_box <- renderValueBox({
    valueBox("ppFEV1 ≥ 70%", "Mild CF", icon = icon("smile"), color = "green")})
  output$fev1_mod_box  <- renderValueBox({
    valueBox("40–69%", "Moderate CF", icon = icon("meh"), color = "yellow")})
  output$fev1_sev_box  <- renderValueBox({
    valueBox("30–39%", "Severe CF", icon = icon("frown"), color = "orange")})
  output$fev1_vsev_box <- renderValueBox({
    valueBox("< 30%", "Very Severe\n(Transplant eval)", icon = icon("hospital"), color = "red")})

  # ── Scenario Comparison ───────────────────────────────────
  scenario_sim_all <- reactive({
    end_h <- input$sc_years * 8760
    sc_defs <- list(
      list(name="Untreated ΔF508",color="#CC0000",IVA=F,LUM=F,TEZ=F,ELX=F,mut=2,tobra=F),
      list(name="Ivacaftor G551D",color="#0070C0",IVA=T,LUM=F,TEZ=F,ELX=F,mut=3,tobra=F),
      list(name="LUM/IVA Orkambi",color="#FF7F00",IVA=T,LUM=T,TEZ=F,ELX=F,mut=2,tobra=F),
      list(name="TEZ/IVA Symdeko",color="#7030A0",IVA=T,LUM=F,TEZ=T,ELX=F,mut=2,tobra=F),
      list(name="ETI Trikafta",   color="#00B050",IVA=T,LUM=F,TEZ=T,ELX=T,mut=2,tobra=F),
      list(name="ETI + Tobramycin",color="#00B0F0",IVA=T,LUM=F,TEZ=T,ELX=T,mut=2,tobra=T)
    )
    results <- lapply(sc_defs, function(sc) {
      evs <- make_ev_shiny(sc$IVA, sc$TEZ, sc$ELX, sc$LUM, end_h)
      params <- list(mutation_type = sc$mut, IVA_dose_flag = as.numeric(sc$IVA),
                     LUM_dose_flag = as.numeric(sc$LUM), TEZ_dose_flag = as.numeric(sc$TEZ),
                     ELX_dose_flag = as.numeric(sc$ELX),
                     TOBRA_flag = as.numeric(sc$tobra || input$sc_tobra),
                     FEV1_0 = input$FEV1_base, Exac_rate_base = input$exac_rate_base)
      df <- run_cf_sim(cf_mod, params, evs, end_h, 24)
      df$Scenario <- sc$name
      df$Color    <- sc$color
      df
    })
    do.call(rbind, results)
  })

  output$scenario_plot <- renderPlotly({
    df <- scenario_sim_all()
    ep <- input$sc_endpoint
    ylabel <- switch(ep,
      ppFEV1    = "ppFEV1 (% predicted)",
      sweat_Cl  = "Sweat Chloride (mmol/L)",
      CFTR_fn   = "CFTR Function (% WT)",
      exac_total= "Cumulative Exacerbations",
      ASL_ht    = "ASL Height (μm)",
      BMI_val   = "BMI (kg/m²)")
    cols <- setNames(unique(df$Color), unique(df$Scenario))
    p <- ggplot(df, aes_string(x = "time_days", y = ep, color = "Scenario")) +
      geom_line(size = 1.2) +
      scale_color_manual(values = cols) +
      labs(x = "Days", y = ylabel, color = "Scenario") +
      theme_bw()
    ggplotly(p)
  })

  output$scenario_summary_table <- renderDT({
    df <- scenario_sim_all()
    ep_52 <- df |>
      filter(time_days >= input$sc_years * 365 * 0.95) |>
      group_by(Scenario) |>
      summarise(
        ppFEV1    = round(mean(ppFEV1), 1),
        Sweat_Cl  = round(mean(sweat_Cl), 1),
        CFTR_fn   = round(mean(CFTR_fn), 1),
        Exac_yr   = round(mean(exac_total) / input$sc_years, 2),
        ASL_um    = round(mean(ASL_ht), 2),
        BMI       = round(mean(BMI_val), 1),
        .groups = "drop"
      )
    datatable(ep_52, options = list(dom = "t", pageLength = 10), rownames = FALSE) |>
      formatStyle("CFTR_fn",
        backgroundColor = styleInterval(c(10, 30), c("#FFCCCC", "#FFFFCC", "#CCFFCC"))) |>
      formatStyle("ppFEV1",
        backgroundColor = styleInterval(c(40, 70), c("#FFCCCC", "#FFE5CC", "#CCFFCC")))
  })

  # ── ASL & Infection tab ───────────────────────────────────
  asl_sim <- reactive({
    params_no  <- list(mutation_type = 2, IVA_dose_flag = 0, LUM_dose_flag = 0,
                       TEZ_dose_flag = 0, ELX_dose_flag = 0, TOBRA_flag = 0,
                       FEV1_0 = input$FEV1_base)
    params_eti <- list(mutation_type = 2, IVA_dose_flag = 1, LUM_dose_flag = 0,
                       TEZ_dose_flag = 1, ELX_dose_flag = 1,
                       TOBRA_flag = as.numeric(input$asl_tobra),
                       FEV1_0 = input$FEV1_base)
    ev_eti <- make_ev_shiny(TRUE, TRUE, TRUE, FALSE, 8760)

    df_no <- run_cf_sim(cf_mod, params_no, NULL)
    df_no$trt <- "Untreated"

    df_eti <- run_cf_sim(cf_mod, params_eti, ev_eti)
    df_eti$trt <- if (input$asl_eti) "ETI (Trikafta)" else "Untreated"

    rbind(df_no, if (input$asl_eti) df_eti else df_no[0,])
  })

  output$asl_plot <- renderPlotly({
    df <- asl_sim()
    p <- ggplot(df, aes(x = time_days, y = ASL_ht, color = trt)) +
      geom_line(size = 1.2) +
      geom_hline(yintercept = 7, linetype = "dashed", color = "green3") +
      geom_hline(yintercept = 4, linetype = "dashed", color = "orange") +
      scale_color_manual(values = c("Untreated" = "#CC0000", "ETI (Trikafta)" = "#00B050")) +
      labs(x = "Days", y = "ASL Height (μm)", color = "Treatment") +
      theme_bw()
    ggplotly(p)
  })

  output$bacteria_plot <- renderPlotly({
    df <- asl_sim()
    p <- ggplot(df, aes(x = time_days, y = logPa, color = trt)) +
      geom_line(size = 1.2) +
      scale_color_manual(values = c("Untreated" = "#CC0000", "ETI (Trikafta)" = "#00B050")) +
      labs(x = "Days", y = "log₁₀(P. aeruginosa CFU/mL)", color = "Treatment") +
      theme_bw()
    ggplotly(p)
  })

  output$inflam_plot <- renderPlotly({
    df <- asl_sim()
    p1 <- ggplot(df, aes(x = time_days, color = trt)) +
      geom_line(aes(y = IL8_val / 100, linetype = "IL-8 (×100 pg/mL)"), size = 1) +
      geom_line(aes(y = Neu_val, linetype = "Neutrophils (×10⁶/mL)"), size = 1) +
      geom_line(aes(y = Dam_val, linetype = "Damage Score"), size = 1) +
      scale_color_manual(values = c("Untreated" = "#CC0000", "ETI (Trikafta)" = "#00B050")) +
      labs(x = "Days", y = "Biomarker (scaled)", color = "Treatment", linetype = "Biomarker") +
      theme_bw()
    ggplotly(p1)
  })
}

shinyApp(ui, server)
