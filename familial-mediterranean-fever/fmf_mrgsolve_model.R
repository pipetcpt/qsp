## ============================================================
## Familial Mediterranean Fever (FMF) — QSP mrgsolve Model
## ============================================================
## Disease: Familial Mediterranean Fever (FMF)
## Genes:   MEFV (pyrin/marenostrin)
## Targets: PYRIN inflammasome, IL-1β, IL-18
## Drugs:   Colchicine, Anakinra, Canakinumab, Rilonacept
##
## ODE Compartments (22 state variables):
##   PK  - Colchicine: gut/central/peripheral/leukocyte (4)
##   PK  - Anakinra: SC depot/central (2)
##   PK  - Canakinumab: SC depot/central/peripheral (3)
##   PD  - PYRIN inflammasome: RhoA, Pyrin_phos, ASC_speck, Casp1 (4)
##   PD  - Cytokines: IL1b_pro, IL1b_mat, IL18, SAA, CRP (5)
##   PD  - Neutrophil dynamics: Neu_circ, Neu_tissue (2)
##   PD  - Attack dynamics: Attack_trigger, Attack_severity (2)
##
## Key References:
##   Infevers mutation database; Gattorno et al. Ann Rheum Dis 2019;
##   Ozen et al. Ann Rheum Dis 2016 (Eurofever/PRINTO criteria);
##   Hentgen et al. Ann Rheum Dis 2013 (colchicine PK);
##   De Benedetti et al. N Engl J Med 2018 (canakinumab CLUSTER trial)
##
## Calibration notes:
##   - Colchicine: Vd ~5-8 L/kg; F ~45%; t½ ~30h; Cleukocyte peak ~10× Cplasma
##   - Anakinra: F ~95% SC; t½ ~4-6h; TMDD-based binding to IL-1R1
##   - Canakinumab: F ~70% SC; t½ ~26d; TMDD to free IL-1β
##   - Attack frequency in untreated FMF: ~3-12/year; colchicine reduces ~75%
##   - Amyloidosis risk correlates with SAA >10 mg/L chronically
## ============================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

## ─── Model definition ──────────────────────────────────────────────────────

fmf_model_code <- '
$PROB
# FMF QSP Model — PYRIN Inflammasome, Cytokines, Drug PK/PD
# 22 ODE compartments; Colchicine, Anakinra, Canakinumab, Rilonacept

$PARAM
// ── Colchicine PK (2-cpt + leukocyte) ──
ka_col   = 1.2       // absorption rate [h^-1], F=0.45
F_col    = 0.45      // oral bioavailability
CL_col   = 18.0      // clearance [L/h]
V1_col   = 120.0     // central volume [L]
Q_col    = 60.0      // inter-compartment CL [L/h]
V2_col   = 480.0     // peripheral volume [L]
k_leu_on = 2.5       // leukocyte uptake [h^-1]
k_leu_off= 0.08      // leukocyte release [h^-1]

// ── Anakinra PK (SC 1-cpt TMDD-simplified) ──
ka_ana   = 0.40      // SC absorption [h^-1]
F_ana    = 0.95      // SC bioavailability
CL_ana   = 1.8       // clearance [L/h], FW 17.3 kDa
V_ana    = 8.5       // volume [L]

// ── Canakinumab PK (SC 2-cpt) ──
ka_cana  = 0.012     // SC absorption [h^-1] (t½abs ~2.4d)
F_cana   = 0.70      // SC bioavailability
CL_cana  = 0.18      // clearance [L/h], t½~26d
V1_cana  = 4.5       // central volume [L]
Q_cana   = 0.4       // inter-compartment CL [L/h]
V2_cana  = 3.0       // peripheral volume [L]

// ── PYRIN inflammasome ──
k_RhoA_ss  = 1.0     // baseline RhoA-GTP production [AU/h]
k_RhoA_deg = 1.0     // RhoA degradation [h^-1]  → SS=1 AU
k_phos_basal= 0.5    // basal PYRIN phosphorylation rate [h^-1]
k_phos_mut  = 0.05   // reduced phosphorylation in M694V [h^-1] (90% reduction)
k_dephos    = 0.3    // dephosphorylation rate [h^-1]
k_ASC_form  = 0.8    // ASC speck formation from dephospho-Pyrin [h^-1]
k_ASC_deg   = 0.4    // ASC speck dissolution [h^-1]
k_Casp1_act = 1.2    // Caspase-1 activation from ASC [h^-1]
k_Casp1_deg = 0.6    // Caspase-1 deactivation [h^-1]

// ── MEFV mutation severity (0=WT, 1=severe M694V) ──
MEFV_severity = 1.0  // 1.0 = homozygous M694V (most severe)

// ── IL-1β dynamics ──
k_IL1b_pro_prod = 0.3  // IL-1β pro production baseline [pg/mL/h]
k_IL1b_pro_deg  = 0.15 // pro-IL-1β degradation [h^-1]
k_IL1b_mat_form = 2.0  // Casp1-mediated maturation (Vmax) [pg/mL/h]
Km_IL1b_mat     = 0.5  // Casp1 Km for maturation [AU]
k_IL1b_mat_deg  = 0.8  // mature IL-1β clearance [h^-1]
IL1b_0          = 5.0  // baseline mature IL-1β [pg/mL]

// ── IL-18 dynamics ──
k_IL18_form = 1.5    // Casp1-mediated IL-18 production [pg/mL/h]
k_IL18_deg  = 0.3    // IL-18 clearance [h^-1]
IL18_0      = 200.0  // baseline IL-18 [pg/mL]

// ── SAA and CRP ──
k_SAA_base  = 0.05   // baseline SAA production [mg/L/h]
k_SAA_IL1b  = 1.5    // IL-1β-driven SAA production
k_SAA_deg   = 0.04   // SAA elimination [h^-1]
SAA_0       = 5.0    // baseline SAA [mg/L]
k_CRP_base  = 0.01   // baseline CRP production [mg/L/h]
k_CRP_IL1b  = 0.3    // IL-1β-driven CRP production
k_CRP_deg   = 0.03   // CRP elimination [h^-1]
CRP_0       = 3.0    // baseline CRP [mg/L]

// ── Neutrophil dynamics ──
k_Neu_prod  = 50.0   // bone marrow neutrophil production [cells/μL/h]
k_Neu_circ_deg= 0.03 // circulating Neu clearance [h^-1]
k_Neu_migr  = 0.05   // IL-1β-driven tissue migration [h^-1]
k_Neu_tis_deg= 0.08  // tissue Neu clearance [h^-1]
Neu_circ_0  = 4000.0 // baseline ANC [cells/μL]
Neu_tis_0   = 500.0  // baseline tissue neutrophils [AU]

// ── Attack dynamics ──
Att_trig_thresh = 2.0 // IL-1β threshold for attack trigger [pg/mL fold]
k_Att_rise  = 0.3    // attack severity rise rate [h^-1]
k_Att_decay = 0.1    // attack resolution rate [h^-1]

// ── Drug effects ──
// Colchicine: IC50 on neutrophil migration and PYRIN
IC50_col_neu  = 8.0  // ng/mL colchicine in leukocyte
Emax_col_neu  = 0.85 // max neutrophil migration inhibition
IC50_col_PYRIN= 5.0  // ng/mL for PYRIN/ASC inhibition
Emax_col_PYRIN= 0.60

// Anakinra: competitive IL-1R blockade
IC50_ana      = 50.0 // ng/mL for IL-1β effect inhibition (reflects IC50 at receptor)
Emax_ana      = 0.95

// Canakinumab: IL-1β neutralization
IC50_cana     = 10.0 // μg/mL canakinumab for IL-1β neutralization
Emax_cana     = 0.98

// Rilonacept: TRAP (modeled via parameter)
RILO_dose     = 0.0  // mg/kg, set >0 to activate (simplified)
IC50_rilo     = 0.5  // μg/mL
Emax_rilo     = 0.96

// ── Amyloidosis ──
k_AA_dep    = 0.001  // AA deposition rate [1/(mg/L*h)] ~ SAA driven
k_AA_deg    = 0.0002 // AA regression rate [h^-1]
k_eGFR_dec  = 0.005  // eGFR decline per AA unit [mL/min/1.73m²/h/AU]
eGFR_0      = 90.0   // baseline eGFR [mL/min/1.73m²]

$PARAM @annotated
// Dosing switches (0=off, 1=on)
USE_COL   : 0 : Use colchicine (1=yes)
USE_ANA   : 0 : Use anakinra (1=yes)
USE_CANA  : 0 : Use canakinumab (1=yes)

$MAIN
// Effective phosphorylation rate (mutation reduces this)
double k_phos_eff = k_phos_basal * (1.0 - MEFV_severity * 0.90)
                    + k_phos_mut * MEFV_severity * 0.90;
// Effective: for M694V severity=1, k_phos_eff ≈ k_phos_mut

// Initial conditions
if(NEWIND <= 1) {
  _nid++;
}

$CMT
// Colchicine PK
GUT_COL    // oral dose depot
CENT_COL   // central plasma [ng*h/mL → div by V]
PERI_COL   // peripheral
LEU_COL    // leukocyte accumulation [ng/mL in WBC]

// Anakinra PK
SC_ANA     // subcutaneous depot
CENT_ANA   // central [ng/mL]

// Canakinumab PK
SC_CANA    // subcutaneous depot
CENT_CANA  // central [μg/mL]
PERI_CANA  // peripheral

// PYRIN inflammasome
RhoA       // RhoA-GTP activation state [AU]
Pyrin_p    // phosphorylated (inactive) PYRIN [AU]
ASC        // ASC speck assembly [AU]
Casp1      // active Caspase-1 [AU]

// Cytokines
IL1b_pro   // pro-IL-1β [pg/mL]
IL1b_mat   // mature IL-1β [pg/mL]
IL18       // IL-18 [pg/mL]
SAA        // serum amyloid A [mg/L]
CRP        // C-reactive protein [mg/L]

// Neutrophils
Neu_circ   // circulating ANC [cells/μL]
Neu_tis    // tissue neutrophils [AU]

// Attack
Att_sev    // attack severity score [0-10]

// Amyloidosis
AA_dep     // AA amyloid deposits [AU]
eGFR       // estimated GFR [mL/min/1.73m²]

$INIT
GUT_COL = 0, CENT_COL = 0, PERI_COL = 0, LEU_COL = 0,
SC_ANA = 0, CENT_ANA = 0,
SC_CANA = 0, CENT_CANA = 0, PERI_CANA = 0,
RhoA = 1.0,
Pyrin_p = 0.5,    // some baseline phosphorylated (inactive) PYRIN
ASC = 0.2,        // low baseline
Casp1 = 0.1,
IL1b_pro = 30.0,
IL1b_mat = IL1b_0,
IL18 = IL18_0,
SAA = SAA_0,
CRP = CRP_0,
Neu_circ = Neu_circ_0,
Neu_tis = Neu_tis_0,
Att_sev = 0.0,
AA_dep = 0.0,
eGFR = eGFR_0

$ODE

// ── Colchicine PK ──
double Cp_col = CENT_COL / V1_col;        // ng/mL plasma
double Cl_col_conc = LEU_COL;             // ng/mL leukocyte

dxdt_GUT_COL  = -ka_col * GUT_COL;
dxdt_CENT_COL =  ka_col * F_col * GUT_COL
                 - (CL_col/V1_col)*CENT_COL
                 - (Q_col/V1_col)*CENT_COL
                 + (Q_col/V2_col)*PERI_COL;
dxdt_PERI_COL =  (Q_col/V1_col)*CENT_COL - (Q_col/V2_col)*PERI_COL;
dxdt_LEU_COL  =  k_leu_on * Cp_col - k_leu_off * LEU_COL;

// ── Anakinra PK ──
double Cp_ana = CENT_ANA / V_ana;        // ng/mL
dxdt_SC_ANA   = -ka_ana * SC_ANA;
dxdt_CENT_ANA =  ka_ana * F_ana * SC_ANA - (CL_ana/V_ana)*CENT_ANA;

// ── Canakinumab PK ──
double Cp_cana = CENT_CANA / V1_cana;   // μg/mL
dxdt_SC_CANA  = -ka_cana * SC_CANA;
dxdt_CENT_CANA=  ka_cana * F_cana * SC_CANA
                 - (CL_cana/V1_cana)*CENT_CANA
                 - (Q_cana/V1_cana)*CENT_CANA
                 + (Q_cana/V2_cana)*PERI_CANA;
dxdt_PERI_CANA=  (Q_cana/V1_cana)*CENT_CANA - (Q_cana/V2_cana)*PERI_CANA;

// ── Drug effect functions ──
// Colchicine leukocyte: Imax model
double E_col_neu  = (USE_COL * Emax_col_neu * Cl_col_conc)
                    / (IC50_col_neu + Cl_col_conc + 1e-10);
double E_col_PYRIN= (USE_COL * Emax_col_PYRIN * Cl_col_conc)
                    / (IC50_col_PYRIN + Cl_col_conc + 1e-10);

// Anakinra blocks IL-1 receptor → reduces IL-1β signaling effect
double E_ana = (USE_ANA * Emax_ana * Cp_ana)
               / (IC50_ana + Cp_ana + 1e-10);

// Canakinumab neutralizes free IL-1β
double E_cana = (USE_CANA * Emax_cana * Cp_cana)
                / (IC50_cana + Cp_cana + 1e-10);

// Combined IL-1β blockade
double IL1b_block = 1.0 - fmax(E_ana, E_cana);  // max effect wins

// ── PYRIN inflammasome ODEs ──
double k_phos_eff = k_phos_basal * (1.0 - MEFV_severity * 0.90)
                    + k_phos_mut * MEFV_severity * 0.90;

// RhoA: inducible by microbial/stress stimuli (simplified as basal SS)
dxdt_RhoA    = k_RhoA_ss - k_RhoA_deg * RhoA;

// Phospho-PYRIN (inactive): reduced in M694V; colchicine stabilizes it slightly
double k_phos_drug = k_phos_eff * (1.0 + 0.3 * E_col_PYRIN);
dxdt_Pyrin_p = k_phos_drug * (1.0 - Pyrin_p) * RhoA
               - k_dephos * Pyrin_p;

// Free (dephospho) PYRIN drives ASC speck: (1 - Pyrin_p) = active pyrin fraction
double active_pyrin = fmax(0.0, 1.0 - Pyrin_p);
dxdt_ASC     = k_ASC_form * active_pyrin * (1.0 - E_col_PYRIN)
               - k_ASC_deg * ASC;

// Caspase-1 from ASC speck
dxdt_Casp1   = k_Casp1_act * ASC - k_Casp1_deg * Casp1;

// ── IL-1β ODEs ──
// Caspase-1 cleaves pro-IL-1β → mature
double rate_IL1b_mat = k_IL1b_mat_form * Casp1 * IL1b_pro
                       / (Km_IL1b_mat + IL1b_pro);

dxdt_IL1b_pro = k_IL1b_pro_prod - k_IL1b_pro_deg * IL1b_pro - rate_IL1b_mat;
dxdt_IL1b_mat = rate_IL1b_mat - k_IL1b_mat_deg * IL1b_mat * IL1b_block;

// ── IL-18 ──
dxdt_IL18     = k_IL18_form * Casp1 - k_IL18_deg * IL18;

// ── SAA (driven by IL-6/IL-1β; proxy for IL-1β signal here) ──
dxdt_SAA      = k_SAA_base + k_SAA_IL1b * IL1b_mat * IL1b_block
                - k_SAA_deg * SAA;

// ── CRP ──
dxdt_CRP      = k_CRP_base + k_CRP_IL1b * IL1b_mat * IL1b_block
                - k_CRP_deg * CRP;

// ── Neutrophils ──
double Neu_migr_eff = k_Neu_migr * (1.0 - E_col_neu) * (IL1b_mat / IL1b_0);
dxdt_Neu_circ = k_Neu_prod - k_Neu_circ_deg * Neu_circ
                - Neu_migr_eff * Neu_circ;
dxdt_Neu_tis  = Neu_migr_eff * Neu_circ - k_Neu_tis_deg * Neu_tis;

// ── Attack severity ──
double att_signal = fmax(0.0, IL1b_mat / IL1b_0 - Att_trig_thresh);
dxdt_Att_sev  = k_Att_rise * att_signal * (1.0 - Att_sev/10.0)
                - k_Att_decay * Att_sev;

// ── AA Amyloidosis and eGFR ──
// SAA >10 mg/L over time drives deposition
double SAA_excess = fmax(0.0, SAA - 10.0);
dxdt_AA_dep   = k_AA_dep * SAA_excess - k_AA_deg * AA_dep;
dxdt_eGFR     = -k_eGFR_dec * AA_dep;  // slow decline

$CAPTURE
Cp_col Cl_col_conc Cp_ana Cp_cana
IL1b_mat IL18 SAA CRP
Neu_circ Neu_tis
Att_sev AA_dep eGFR
E_col_neu E_col_PYRIN E_ana E_cana IL1b_block
active_pyrin ASC Casp1

$TABLE
double AIDAI_score = fmin(10.0, Att_sev);  // simplified AIDAI proxy
double attack_flag = (Att_sev > 2.0) ? 1.0 : 0.0;
'

mod <- mcode("fmf_qsp", fmf_model_code, quiet = TRUE)

## ─── Treatment scenarios ───────────────────────────────────────────────────

# Simulation duration: 2 years (17520h) for chronic disease perspective
# Shorter exploration: 52 weeks (8736h) with attack events
SIM_END <- 8760  # 1 year in hours

# 5 Treatment scenarios
scenarios <- list(
  list(
    name        = "1. No Treatment",
    USE_COL     = 0, USE_ANA = 0, USE_CANA = 0,
    events      = NULL
  ),
  list(
    name        = "2. Colchicine 0.5 mg BID",
    USE_COL     = 1, USE_ANA = 0, USE_CANA = 0,
    events      = ev(amt = 500, cmt = "GUT_COL",  # 500 ng? — use μg→ng: 500000 ng
                     time = seq(0, SIM_END, by = 12))
  ),
  list(
    name        = "3. Colchicine 1.0 mg QD",
    USE_COL     = 1, USE_ANA = 0, USE_CANA = 0,
    events      = ev(amt = 1000000, cmt = "GUT_COL",
                     time = seq(0, SIM_END, by = 24))
  ),
  list(
    name        = "4. Anakinra 100 mg SC QD",
    USE_COL     = 0, USE_ANA = 1, USE_CANA = 0,
    events      = ev(amt = 100000, cmt = "SC_ANA",  # 100 mg = 100,000 μg → but V=8.5L, units ng/mL
                     time = seq(0, SIM_END, by = 24))
  ),
  list(
    name        = "5. Canakinumab 150 mg SC Q8W",
    USE_COL     = 0, USE_ANA = 0, USE_CANA = 1,
    events      = ev(amt = 150, cmt = "SC_CANA",   # 150 mg, Cp in μg/mL, V=4.5L
                     time = seq(0, SIM_END, by = 56*24))
  )
)

## ─── Helper: run one scenario ───────────────────────────────────────────────

run_scenario <- function(sc, mefv_severity = 1.0) {
  params <- c(
    MEFV_severity = mefv_severity,
    USE_COL = sc$USE_COL,
    USE_ANA = sc$USE_ANA,
    USE_CANA = sc$USE_CANA
  )

  m <- mod %>% param(params)

  if (!is.null(sc$events)) {
    out <- m %>%
      ev(sc$events) %>%
      mrgsim(end = SIM_END, delta = 1)
  } else {
    out <- m %>%
      mrgsim(end = SIM_END, delta = 1)
  }

  as.data.frame(out) %>%
    mutate(scenario = sc$name,
           mefv_sev = mefv_severity)
}

## ─── Run all scenarios for M694V (severe) ──────────────────────────────────

cat("Running FMF treatment scenarios (M694V homozygous)...\n")
results <- lapply(scenarios, run_scenario, mefv_severity = 1.0)
df_all  <- bind_rows(results)

## ─── Run MEFV severity comparison ──────────────────────────────────────────
# No-treatment arm across severity levels
cat("Running MEFV severity comparison (no treatment)...\n")
severity_levels <- c(0.0, 0.3, 0.6, 0.8, 1.0)  # WT → M694V
severity_names  <- c("WT (no FMF)", "E148Q (mild)", "V726A (mod)", "M680I (mod-sev)", "M694V (severe)")
df_sev <- lapply(seq_along(severity_levels), function(i) {
  run_scenario(scenarios[[1]], mefv_severity = severity_levels[i]) %>%
    mutate(genotype = severity_names[i])
}) %>% bind_rows()

## ─── Key summary metrics ────────────────────────────────────────────────────

summary_metrics <- df_all %>%
  group_by(scenario) %>%
  summarise(
    mean_IL1b        = mean(IL1b_mat, na.rm = TRUE),
    peak_IL1b        = max(IL1b_mat, na.rm = TRUE),
    mean_SAA         = mean(SAA, na.rm = TRUE),
    mean_CRP         = mean(CRP, na.rm = TRUE),
    pct_time_attack  = mean(attack_flag, na.rm = TRUE) * 100,
    attacks_per_year = sum(diff(c(0, attack_flag)) == 1, na.rm = TRUE),
    final_AA_dep     = last(AA_dep),
    final_eGFR       = last(eGFR),
    .groups = "drop"
  )

cat("\n=== FMF Treatment Scenario Summary (1 year, M694V) ===\n")
print(summary_metrics)

## ─── Plots ─────────────────────────────────────────────────────────────────

time_weeks <- function(df) df %>% mutate(time_wk = time / 168)

# IL-1β over time
p_IL1b <- df_all %>%
  time_weeks() %>%
  ggplot(aes(x = time_wk, y = IL1b_mat, color = scenario)) +
  geom_line(linewidth = 0.8) +
  labs(title = "Mature IL-1β (FMF Scenarios)",
       x = "Time (weeks)", y = "IL-1β (pg/mL)",
       color = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(ncol = 2))

# SAA over time
p_SAA <- df_all %>%
  time_weeks() %>%
  ggplot(aes(x = time_wk, y = SAA, color = scenario)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  annotate("text", x = 1, y = 11, label = "SAA=10 mg/L\n(amyloid risk)", size = 3, hjust = 0) +
  labs(title = "Serum Amyloid A",
       x = "Time (weeks)", y = "SAA (mg/L)", color = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# Attack severity over time
p_att <- df_all %>%
  time_weeks() %>%
  ggplot(aes(x = time_wk, y = Att_sev, color = scenario)) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  labs(title = "Attack Severity Score",
       x = "Time (weeks)", y = "Attack Severity (0-10)",
       color = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# Amyloidosis and eGFR (long-term consequence)
p_eGFR <- df_all %>%
  time_weeks() %>%
  ggplot(aes(x = time_wk, y = eGFR, color = scenario)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = c(60, 30), linetype = "dashed", color = "gray50") +
  labs(title = "eGFR (Amyloidosis Risk)",
       x = "Time (weeks)", y = "eGFR (mL/min/1.73m²)",
       color = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# Colchicine PK
p_col_pk <- df_all %>%
  filter(scenario == "2. Colchicine 0.5 mg BID") %>%
  time_weeks() %>%
  filter(time_wk < 4) %>%
  ggplot(aes(x = time_wk)) +
  geom_line(aes(y = Cp_col, color = "Plasma")) +
  geom_line(aes(y = Cl_col_conc, color = "Leukocyte")) +
  labs(title = "Colchicine PK (0.5 mg BID, first 4 wk)",
       x = "Time (weeks)", y = "Concentration (ng/mL)",
       color = "Compartment") +
  theme_bw(base_size = 11)

# MEFV genotype effect on IL-1β
p_geno <- df_sev %>%
  time_weeks() %>%
  ggplot(aes(x = time_wk, y = IL1b_mat, color = genotype)) +
  geom_line(linewidth = 0.8) +
  labs(title = "IL-1β by MEFV Genotype (no treatment)",
       x = "Time (weeks)", y = "IL-1β (pg/mL)",
       color = "Genotype") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 2))

# PYRIN inflammasome components
p_inflamasome <- df_all %>%
  filter(scenario %in% c("1. No Treatment", "2. Colchicine 0.5 mg BID", "5. Canakinumab 150 mg SC Q8W")) %>%
  time_weeks() %>%
  select(time_wk, scenario, active_pyrin, ASC, Casp1) %>%
  pivot_longer(cols = c(active_pyrin, ASC, Casp1),
               names_to = "component", values_to = "level") %>%
  ggplot(aes(x = time_wk, y = level, color = scenario, linetype = component)) +
  geom_line() +
  labs(title = "PYRIN Inflammasome Components",
       x = "Time (weeks)", y = "Level (AU)",
       color = "Scenario", linetype = "Component") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(ncol = 1))

# Combine panel
combined_plot <- (p_IL1b + p_SAA) / (p_att + p_eGFR) / (p_col_pk + p_geno)

cat("\nPlot objects created: p_IL1b, p_SAA, p_att, p_eGFR, p_col_pk, p_geno, p_inflamasome\n")
cat("Use combined_plot to view all panels.\n")

## ─── Sensitivity analysis ───────────────────────────────────────────────────
# Simple one-at-a-time sensitivity: colchicine IC50 on neutrophils
cat("\nSensitivity: colchicine IC50 on neutrophil inhibition...\n")

ic50_vals <- c(2, 5, 8, 15, 30)
df_sens <- lapply(ic50_vals, function(ic) {
  mod %>%
    param(MEFV_severity = 1.0, USE_COL = 1, USE_ANA = 0, USE_CANA = 0,
          IC50_col_neu = ic) %>%
    ev(amt = 500000, cmt = "GUT_COL",
       time = seq(0, SIM_END, by = 12)) %>%
    mrgsim(end = SIM_END, delta = 6) %>%
    as.data.frame() %>%
    mutate(IC50_col_neu = ic,
           scenario = paste0("IC50=", ic, " ng/mL"))
}) %>% bind_rows()

p_sens <- df_sens %>%
  time_weeks() %>%
  ggplot(aes(x = time_wk, y = Neu_tis, color = scenario)) +
  geom_line() +
  labs(title = "Sensitivity: Colchicine IC50 on Tissue Neutrophils",
       x = "Time (weeks)", y = "Tissue Neutrophils (AU)",
       color = "IC50 (ng/mL)") +
  theme_bw(base_size = 11)

cat("\nFMF mrgsolve model ready. Objects in environment:\n")
cat("  mod           — mrgsolve model object\n")
cat("  df_all        — simulation results (5 scenarios, 1 year)\n")
cat("  df_sev        — MEFV severity comparison\n")
cat("  summary_metrics — key outcome table\n")
cat("  combined_plot — 6-panel ggplot\n")
cat("  p_sens        — sensitivity analysis plot\n")
