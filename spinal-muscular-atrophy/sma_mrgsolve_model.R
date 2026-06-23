## =============================================================================
## Spinal Muscular Atrophy (SMA) — mrgsolve QSP Model
## =============================================================================
## Disease: Spinal Muscular Atrophy (5q-SMA)
## Model covers:
##   1. SMN2 pre-mRNA alternative splicing (exon 7 inclusion dynamics)
##   2. SMN protein homeostasis (FL-SMN pool)
##   3. Motor neuron survival/degeneration (alpha-MN pool)
##   4. Neuromuscular junction maturation
##   5. Skeletal muscle mass dynamics
##   6. Clinical endpoint proxies (CMAP, HFMSE/CHOP-INTEND, FVC)
##   7. Nusinersen intrathecal PK (2-compartment CSF + CNS tissue)
##   8. Risdiplam oral PK (2-compartment with CNS penetration)
##   9. Onasemnogene (Zolgensma) IV gene therapy PK/transduction
##
## Calibration references:
##   - Darras et al. (2019) NEJM — ENDEAR trial (nusinersen type I)
##   - Mercuri et al. (2018) NEJM — CHERISH trial (nusinersen type II/III)
##   - Baranello et al. (2021) NEJM — firefish/sunfish (risdiplam)
##   - Day et al. (2021) NEJM — STR1VE (Zolgensma)
##   - Kletzl et al. (2019) J Pharmacokinet Pharmacodyn — nusinersen PK
##   - Poirier et al. (2021) CPT:PSP — risdiplam PK
##   - Al-Zaidy et al. (2019) Mol Ther — Zolgensma biodistribution
##
## Usage:
##   library(mrgsolve); library(tidyverse)
##   source("sma_mrgsolve_model.R")
##   mod <- sma_model()
##   out <- mrgsim(mod, ev_nusinersen(), delta=1, end=730)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

## ─────────────────────────────────────────────────────────────
## 1. Model Definition
## ─────────────────────────────────────────────────────────────
sma_model <- function() {
  mrgsolve::mcode("sma_qsp", '
$PROB
SMA QSP model: SMN biology, motor neuron degeneration,
nusinersen / risdiplam / Zolgensma PK-PD

$PARAM @annotated
// --- SMN2 Splicing Parameters ---
Emax_NUS    : 0.60 : Nusinersen max fractional increase in exon-7 inclusion (0-1)
EC50_NUS    : 5.0  : Nusinersen CNS tissue EC50 for exon-7 inclusion (ng/g)
hill_NUS    : 1.5  : Hill coefficient, nusinersen splicing effect
Emax_RIS    : 0.50 : Risdiplam max fractional increase in exon-7 inclusion
EC50_RIS    : 80.0 : Risdiplam plasma EC50 for splicing effect (ng/mL)
hill_RIS    : 1.2  : Hill coefficient, risdiplam splicing effect
E7I_base    : 0.10 : Baseline exon-7 inclusion from SMN2 (fraction, ~10%)
E7I_max     : 0.90 : Maximum achievable exon-7 inclusion fraction

// --- SMN mRNA & Protein Dynamics ---
k_SMN2_txn  : 1.0  : SMN2 transcription rate constant (relative units/day)
k_FL_deg    : 0.693: FL-SMN mRNA degradation rate (t½ ~1 day; /day)
k_d7_deg    : 2.079: SMN-Δ7 mRNA degradation rate (t½ ~8 h; /day)
k_prot_syn  : 2.0  : SMN protein synthesis rate per FL-mRNA unit
k_prot_deg  : 0.231: SMN protein degradation rate (t½ ~3 days; /day)
SMN_thresh  : 0.30 : SMN protein threshold (fraction of normal) below which MN death begins

// --- Motor Neuron Pool ---
MN0         : 1.0  : Initial motor neuron pool (normalized to 1.0 = 100%)
k_MN_death  : 0.002: Daily MN death rate in SMA (fraction/day)
k_MN_rescue : 0.04 : Rate at which SMN restoration rescues MNs (/day)
k_MN_spont  : 0.0001: Spontaneous MN death rate (aging; /day)
MN_min      : 0.05 : Minimum MN pool fraction (5% irreducible)
d_MN_hill   : 2.0  : Hill coefficient for SMN-dependent MN death rate

// --- Neuromuscular Junction ---
k_NMJ_mat   : 0.01 : Rate of NMJ maturation (/day)
k_NMJ_det   : 0.005: Rate of NMJ deterioration with MN loss (/day)
NMJ_max     : 1.0  : Maximum NMJ maturity score

// --- Skeletal Muscle ---
Muscle0     : 1.0  : Initial muscle mass (normalized)
k_muscle_at : 0.004: Daily muscle atrophy rate (denervation; /day)
k_muscle_gr : 0.002: Daily muscle recovery/growth rate (/day)
muscle_min  : 0.10 : Minimum muscle mass fraction

// --- Clinical Endpoint Scaling ---
CMAP_max    : 10.0 : Maximum CMAP amplitude (mV, normal adult)
CHOP_max    : 64.0 : Maximum CHOP-INTEND score
HFMSE_max   : 66.0 : Maximum HFMSE score
FVC_max     : 100.0: Maximum FVC % predicted
RULM_max    : 37.0 : Maximum RULM score

// ------- Nusinersen PK Parameters -------
// Intrathecal 2-compartment (CSF lumbar → cervical → CNS tissue)
V_CSF_L     : 30.0 : CSF lumbar volume (mL)
V_CSF_C     : 70.0 : CSF cervical+cranial volume (mL)
V_CNS       : 1500.0: CNS tissue volume (g, for ng/g)
Q_CSF       : 0.43 : CSF bulk flow rate (mL/min = 619 mL/day)
CL_ASO_CSF  : 0.05 : ASO clearance from CSF lumbar (mL/min = 72 mL/day)
CL_ASO_CNS  : 0.001: ASO elimination from CNS tissue (fraction/day)
ka_CNS_NUS  : 0.05 : Uptake rate from CSF to CNS tissue (/day)
Vmax_NUS    : 50.0 : Max uptake capacity from CSF to CNS (ng/day)
Km_NUS      : 5.0  : Michaelis constant for CNS uptake (ng/mL)

// ------- Risdiplam PK Parameters -------
// 1-compartment with CNS penetration
ka_RIS      : 1.5  : Oral absorption rate constant (/h, converted in ODEs)
F_RIS       : 0.99 : Oral bioavailability
Vd_RIS      : 70.0 : Volume of distribution (L)
CL_RIS      : 1.0  : Clearance (L/h)
Kp_brain    : 0.5  : Brain:plasma Kp ratio for risdiplam
fu_RIS      : 0.07 : Unbound fraction in plasma

// ------- Onasemnogene (Zolgensma) PK Parameters -------
k_vg_clear  : 2.31 : Vector genome clearance rate from plasma (/day, t½ ~7h)
k_MN_trans  : 0.001: Rate of motor neuron transduction (vg/nucleus/day)
k_tg_txn    : 0.5  : Transgene transcription rate (relative units/day)
k_tg_mRNA_deg : 0.693: Transgene mRNA degradation rate (t½ ~1 day; /day)
k_tg_prot   : 2.0  : Transgene protein synthesis rate
eff_tg      : 1.0  : Transgene expression efficiency (0–1)
AAV9_Ab_block: 0.0 : Anti-AAV9 antibody neutralization (0=no block, 1=full block)

// ------- Disease Subtype Parameters -------
SMN2_copies : 2.0  : Number of SMN2 gene copies (1–4)
disease_onset_age: 0.0: Age at symptom onset (months, for natural history)

$CMT @annotated
// Nusinersen PK compartments
A_CSF_L     : CSF lumbar nusinersen (ng)
A_CSF_C     : CSF cervical nusinersen (ng)
A_CNS_NUS   : CNS tissue nusinersen (ng/g equivalent)

// Risdiplam PK compartments
A_gut_RIS   : GI tract risdiplam (mg)
A_plasma_RIS: Plasma risdiplam (mg)
A_CNS_RIS   : CNS risdiplam (mg)

// Zolgensma PK
A_plasma_ZOL: Plasma vector genome concentration (vg/mL * Vdist)
A_MN_ZOL    : Transduced MN vg load (vg/nucleus)
A_tg_mRNA   : Transgene-derived FL-SMN mRNA (relative units)

// Disease biology
FL_SMN_mRNA : Full-length SMN mRNA (relative units)
dSMN_mRNA   : SMN-delta7 mRNA (relative units)
SMN_prot    : Full-length SMN protein pool (normalized, 0-1)
MN_pool     : Alpha motor neuron pool (normalized, 0-1)
NMJ_score   : NMJ maturation score (0-1)
Muscle_mass : Skeletal muscle mass (normalized, 0-1)

// Cumulative endpoints
AUC_SMN     : AUC of SMN protein (for efficacy assessment)
MN_lost     : Cumulative MN loss (fraction)

$MAIN
// SMN2 copy-number scaling (more SMN2 copies = more baseline FL-SMN)
double SMN2_scale = SMN2_copies / 2.0;  // normalized to 2 copies

// Exon-7 inclusion at t=0 (drug-naive)
double E7I_ss = E7I_base * SMN2_scale;

// Initial SMN mRNA at steady state (no drug)
double FL_mRNA_ss = k_SMN2_txn * E7I_ss / k_FL_deg;
double dSMN_mRNA_ss = k_SMN2_txn * (1.0 - E7I_ss) / k_d7_deg;

// Initial SMN protein (normalized: normal = 1.0 for SMN2_copies=4)
double SMN_prot_ss = k_prot_syn * FL_mRNA_ss / k_prot_deg;

// Initial conditions
FL_SMN_mRNA_0  = FL_mRNA_ss;
dSMN_mRNA_0    = dSMN_mRNA_ss;
SMN_prot_0     = SMN_prot_ss;
MN_pool_0      = MN0;
NMJ_score_0    = 0.8 * MN0;   // NMJ roughly proportional to MN at baseline
Muscle_mass_0  = Muscle0;
AUC_SMN_0      = 0.0;
MN_lost_0      = 0.0;

$ODE
// ─────────────────────────────────────────────────────────
// Nusinersen PK (intrathecal)
// C_CSF_L (ng/mL) = A_CSF_L / V_CSF_L
// ─────────────────────────────────────────────────────────
double C_CSF_L = A_CSF_L / V_CSF_L;
double C_CSF_C = A_CSF_C / V_CSF_C;
double C_CNS_N = A_CNS_NUS / V_CNS;   // ng/g in CNS tissue

double flow_L_to_C = (Q_CSF * 1440.0) * C_CSF_L;  // ng/day (Q in mL/min → mL/day)
double CL_L = CL_ASO_CSF * 1440.0;                 // mL/day
double uptake_NUS = (Vmax_NUS * C_CSF_L) / (Km_NUS + C_CSF_L);  // nonlinear CNS uptake

dxdt_A_CSF_L  = - flow_L_to_C - CL_L * C_CSF_L - uptake_NUS;
dxdt_A_CSF_C  = flow_L_to_C - CL_L * C_CSF_C;
dxdt_A_CNS_NUS = uptake_NUS - CL_ASO_CNS * A_CNS_NUS;

// ─────────────────────────────────────────────────────────
// Risdiplam PK (oral, per-day dosing)
// ─────────────────────────────────────────────────────────
double C_plasma_RIS = A_plasma_RIS / Vd_RIS;  // mg/L = ug/mL (mg/L * 1000 = ng/mL)
double C_plasma_ng  = C_plasma_RIS * 1000.0;  // convert to ng/mL for PD
double C_CNS_ng     = (A_CNS_RIS / Vd_RIS) * 1000.0 * Kp_brain;  // ng/mL CNS

double ka_day = ka_RIS * 24.0;  // convert /h to /day
double CL_day = CL_RIS * 24.0;

dxdt_A_gut_RIS    = -ka_day * A_gut_RIS;
dxdt_A_plasma_RIS =  ka_day * F_RIS * A_gut_RIS - CL_day * C_plasma_RIS;
dxdt_A_CNS_RIS    =  ka_day * Kp_brain * A_gut_RIS - CL_day * Kp_brain * (A_CNS_RIS / Vd_RIS);

// ─────────────────────────────────────────────────────────
// Zolgensma PK (IV AAV9 gene therapy)
// ─────────────────────────────────────────────────────────
double C_vg_plasma = A_plasma_ZOL;  // vector genomes in plasma compartment

dxdt_A_plasma_ZOL = -k_vg_clear * A_plasma_ZOL;
dxdt_A_MN_ZOL     = k_MN_trans * A_plasma_ZOL * (1.0 - AAV9_Ab_block);   // transduction of MNs
dxdt_A_tg_mRNA    = k_tg_txn * eff_tg * A_MN_ZOL - k_tg_mRNA_deg * A_tg_mRNA;

// ─────────────────────────────────────────────────────────
// Drug Effects on Exon-7 Inclusion
// ─────────────────────────────────────────────────────────
// Nusinersen effect (ISS-N1 blockade increases E7I)
double Emax_eff_NUS = Emax_NUS * pow(C_CNS_N, hill_NUS) /
                      (pow(EC50_NUS, hill_NUS) + pow(C_CNS_N, hill_NUS));

// Risdiplam effect (SRSF1/Tra2b enhancement)
double Emax_eff_RIS = Emax_RIS * pow(C_plasma_ng, hill_RIS) /
                      (pow(EC50_RIS, hill_RIS) + pow(C_plasma_ng, hill_RIS));

// Combined exon-7 inclusion rate (max at E7I_max)
double E7I_current = E7I_base + (E7I_max - E7I_base) * (Emax_eff_NUS + Emax_eff_RIS -
                     Emax_eff_NUS * Emax_eff_RIS);  // independent combination
E7I_current = fmin(E7I_current, E7I_max);
E7I_current = fmax(E7I_current, E7I_base);

// Scale by SMN2 copy number
E7I_current *= SMN2_scale;

// ─────────────────────────────────────────────────────────
// SMN mRNA Dynamics
// ─────────────────────────────────────────────────────────
double SMN_from_ZOL = k_tg_prot > 0 ? A_tg_mRNA : 0.0;  // transgene contribution

dxdt_FL_SMN_mRNA = k_SMN2_txn * E7I_current - k_FL_deg * FL_SMN_mRNA;
dxdt_dSMN_mRNA   = k_SMN2_txn * (1.0 - fmin(E7I_current, 1.0)) - k_d7_deg * dSMN_mRNA;

// ─────────────────────────────────────────────────────────
// SMN Protein Dynamics
// ─────────────────────────────────────────────────────────
double SMN_synthesis = k_prot_syn * (FL_SMN_mRNA + A_tg_mRNA);  // includes transgene
dxdt_SMN_prot = SMN_synthesis - k_prot_deg * SMN_prot;
double SMN_norm = fmax(SMN_prot, 0.0);  // bounded at 0

// ─────────────────────────────────────────────────────────
// Motor Neuron Pool Dynamics
// ─────────────────────────────────────────────────────────
// MN death rate depends on SMN protein (sigmoidal, below threshold)
double SMN_ratio = SMN_norm / SMN_thresh;
double death_rate_factor = 1.0 / (1.0 + pow(SMN_ratio, d_MN_hill));
double k_death_effective = k_MN_death * death_rate_factor + k_MN_spont;

// MN rescue proportional to SMN recovery above threshold
double rescue_factor = fmax(0.0, (SMN_norm - SMN_thresh) / (1.0 - SMN_thresh));
double k_rescue_effective = k_MN_rescue * rescue_factor;

dxdt_MN_pool = k_rescue_effective * (MN0 - MN_pool) - k_death_effective * MN_pool;
dxdt_MN_pool = fmax(dxdt_MN_pool, -(MN_pool - MN_min));  // don't go below minimum

// ─────────────────────────────────────────────────────────
// NMJ Maturation
// ─────────────────────────────────────────────────────────
dxdt_NMJ_score = k_NMJ_mat * MN_pool * (NMJ_max - NMJ_score)
                 - k_NMJ_det * (1.0 - MN_pool) * NMJ_score;

// ─────────────────────────────────────────────────────────
// Skeletal Muscle Mass
// ─────────────────────────────────────────────────────────
double innervation = NMJ_score * MN_pool;
dxdt_Muscle_mass = k_muscle_gr * innervation * (Muscle0 - Muscle_mass)
                   - k_muscle_at * (1.0 - innervation) * Muscle_mass;
dxdt_Muscle_mass = fmax(dxdt_Muscle_mass, -(Muscle_mass - muscle_min));

// ─────────────────────────────────────────────────────────
// Cumulative Endpoints
// ─────────────────────────────────────────────────────────
dxdt_AUC_SMN = SMN_prot;
dxdt_MN_lost = fmax(0.0, -dxdt_MN_pool);

$TABLE
// Derived clinical variables
double CMAP      = CMAP_max * MN_pool * NMJ_score;
double FVC_pct   = FVC_max * Muscle_mass * 0.85 + 15.0;  // respiratory muscle proxy
double HFMSE     = HFMSE_max * pow(Muscle_mass, 1.5) * MN_pool;
double CHOP      = CHOP_max * Muscle_mass * NMJ_score;
double RULM      = RULM_max * Muscle_mass * MN_pool;

// PK readouts
double C_CSF_L_out  = A_CSF_L / V_CSF_L;
double C_CNS_NUS_out = A_CNS_NUS / V_CNS;
double C_plas_RIS   = (A_plasma_RIS / Vd_RIS) * 1000.0;  // ng/mL

// Exon-7 inclusion (approximation from current drug levels)
double E7I_out = E7I_base + (E7I_max - E7I_base) *
  (Emax_NUS * pow(C_CNS_NUS_out, hill_NUS) / (pow(EC50_NUS, hill_NUS) + pow(C_CNS_NUS_out, hill_NUS)) +
   Emax_RIS * pow(C_plas_RIS, hill_RIS) / (pow(EC50_RIS, hill_RIS) + pow(C_plas_RIS, hill_RIS)));

capture CMAP      = CMAP;
capture FVC       = FVC_pct;
capture HFMSE     = HFMSE;
capture CHOP_INTEND = CHOP;
capture RULM      = RULM;
capture E7_inclusion = E7I_out;
capture SMN_protein  = SMN_prot;
capture MN_fraction  = MN_pool;
capture NMJ_maturity = NMJ_score;
capture C_CSF_lumbar = C_CSF_L_out;
capture C_CNS_nusinersen = C_CNS_NUS_out;
capture C_plasma_risdiplam = C_plas_RIS;
capture Transgene_mRNA = A_tg_mRNA;
')
}

## ─────────────────────────────────────────────────────────────
## 2. Dosing Event Builders
## ─────────────────────────────────────────────────────────────

# Nusinersen (intrathecal): ENDEAR/CHERISH loading + maintenance
# cmt=1 (A_CSF_L), intrathecal 12 mg = 12000 ng injected
ev_nusinersen <- function(start_day = 0) {
  loading <- c(0, 14, 28, 63)
  maintenance <- seq(63 + 120, 730, by = 120)  # every 4 months
  days <- unique(c(loading + start_day, maintenance + start_day))
  days <- days[days <= 730]
  mrgsolve::ev(
    ID    = 1,
    time  = days,
    amt   = 12000,    # 12 mg = 12000 ng in CSF lumbar (cmt=1)
    cmt   = 1,
    evid  = 1,
    addl  = 0
  )
}

# Risdiplam oral 5 mg daily (fixed dose adult)
# amt in mg, cmt=4 (A_gut_RIS)
ev_risdiplam <- function(start_day = 0, end_day = 730, dose_mg = 5) {
  mrgsolve::ev(
    ID    = 1,
    time  = seq(start_day, end_day, by = 1),
    amt   = dose_mg,
    cmt   = 4,
    evid  = 1
  )
}

# Risdiplam weight-based (pediatric) 0.2 mg/kg
ev_risdiplam_peds <- function(weight_kg = 15, start_day = 0, end_day = 730) {
  ev_risdiplam(start_day, end_day, dose_mg = 0.2 * weight_kg)
}

# Zolgensma single IV dose (1.1 × 10^14 vg/kg, 15 kg child = 1.65 × 10^15 vg)
# Expressed as normalized units; cmt=7 (A_plasma_ZOL)
ev_zolgensma <- function(dose_vg = 1.65e15) {
  mrgsolve::ev(
    ID   = 1,
    time = 0,
    amt  = dose_vg / 1e14,  # normalized to 10^14 units
    cmt  = 7,
    evid = 1
  )
}

## ─────────────────────────────────────────────────────────────
## 3. Simulation Scenarios
## ─────────────────────────────────────────────────────────────
run_scenarios <- function() {
  mod <- sma_model()

  # Parameter sets for disease subtypes
  params_type1 <- list(SMN2_copies = 2, MN0 = 1.0, k_MN_death = 0.004)
  params_type2 <- list(SMN2_copies = 3, MN0 = 1.0, k_MN_death = 0.002)
  params_type3 <- list(SMN2_copies = 4, MN0 = 1.0, k_MN_death = 0.001)

  scenarios <- list(
    # ── Scenario 1: Untreated SMA Type I natural history ──────
    list(
      name   = "SMA Type I — No Treatment",
      params = params_type1,
      events = mrgsolve::ev(time = 9999, amt = 0, cmt = 1)  # no drug
    ),

    # ── Scenario 2: Nusinersen in SMA Type I (ENDEAR) ─────────
    list(
      name   = "SMA Type I — Nusinersen",
      params = params_type1,
      events = ev_nusinersen(start_day = 0)
    ),

    # ── Scenario 3: Risdiplam in SMA Type II (SUNFISH) ────────
    list(
      name   = "SMA Type II — Risdiplam",
      params = params_type2,
      events = ev_risdiplam(dose_mg = 5)
    ),

    # ── Scenario 4: Zolgensma in presymptomatic SMA type I ───
    list(
      name   = "Presymptomatic SMA — Zolgensma",
      params = modifyList(params_type1, list(MN0 = 0.95)),
      events = ev_zolgensma(dose_vg = 1.65e15)
    ),

    # ── Scenario 5: Nusinersen late start (type II, 2 yrs) ───
    list(
      name   = "SMA Type II — Nusinersen Late Start",
      params = params_type2,
      events = ev_nusinersen(start_day = 365)  # start after 1 year
    ),

    # ── Scenario 6: Risdiplam pediatric weight-based ─────────
    list(
      name   = "SMA Type II — Risdiplam Pediatric (15 kg)",
      params = params_type2,
      events = ev_risdiplam_peds(weight_kg = 15)
    )
  )

  results <- lapply(seq_along(scenarios), function(i) {
    sc  <- scenarios[[i]]
    mod2 <- mrgsolve::param(mod, sc$params)
    out  <- mrgsolve::mrgsim(mod2, sc$events, delta = 1, end = 730, obsonly = TRUE)
    df   <- as.data.frame(out)
    df$Scenario <- sc$name
    df
  })

  bind_rows(results)
}

## ─────────────────────────────────────────────────────────────
## 4. Visualization
## ─────────────────────────────────────────────────────────────
plot_results <- function(df) {
  pal <- c(
    "SMA Type I — No Treatment"                     = "#B71C1C",
    "SMA Type I — Nusinersen"                       = "#4CAF50",
    "SMA Type II — Risdiplam"                       = "#2196F3",
    "Presymptomatic SMA — Zolgensma"                = "#9C27B0",
    "SMA Type II — Nusinersen Late Start"           = "#FF9800",
    "SMA Type II — Risdiplam Pediatric (15 kg)"    = "#00BCD4"
  )
  lt <- c(
    "SMA Type I — No Treatment"                     = "dashed",
    "SMA Type I — Nusinersen"                       = "solid",
    "SMA Type II — Risdiplam"                       = "solid",
    "Presymptomatic SMA — Zolgensma"                = "solid",
    "SMA Type II — Nusinersen Late Start"           = "dotdash",
    "SMA Type II — Risdiplam Pediatric (15 kg)"    = "dotted"
  )

  p1 <- ggplot(df, aes(time, SMN_protein, color = Scenario, linetype = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = pal) + scale_linetype_manual(values = lt) +
    labs(title = "SMN Protein Pool Over Time",
         x = "Days", y = "SMN Protein (normalized)") +
    theme_bw(14) + theme(legend.position = "bottom", legend.text = element_text(size = 8))

  p2 <- ggplot(df, aes(time, MN_fraction, color = Scenario, linetype = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = pal) + scale_linetype_manual(values = lt) +
    labs(title = "Motor Neuron Pool Over Time",
         x = "Days", y = "MN Pool (fraction of baseline)") +
    theme_bw(14) + theme(legend.position = "bottom", legend.text = element_text(size = 8))

  p3 <- ggplot(df, aes(time, CMAP, color = Scenario, linetype = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = pal) + scale_linetype_manual(values = lt) +
    labs(title = "CMAP Amplitude Over Time",
         x = "Days", y = "CMAP (mV)") +
    theme_bw(14) + theme(legend.position = "bottom", legend.text = element_text(size = 8))

  p4 <- ggplot(df, aes(time, CHOP_INTEND, color = Scenario, linetype = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = pal) + scale_linetype_manual(values = lt) +
    labs(title = "CHOP-INTEND Score Over Time",
         x = "Days", y = "CHOP-INTEND (0–64)") +
    theme_bw(14) + theme(legend.position = "bottom", legend.text = element_text(size = 8))

  p5 <- ggplot(df, aes(time, FVC, color = Scenario, linetype = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = pal) + scale_linetype_manual(values = lt) +
    labs(title = "FVC % Predicted Over Time",
         x = "Days", y = "FVC (% predicted)") +
    theme_bw(14) + theme(legend.position = "bottom", legend.text = element_text(size = 8))

  p6 <- ggplot(df, aes(time, E7_inclusion, color = Scenario, linetype = Scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = pal) + scale_linetype_manual(values = lt) +
    labs(title = "Exon 7 Inclusion Rate (SMN2)",
         x = "Days", y = "Exon 7 Inclusion (fraction)") +
    theme_bw(14) + theme(legend.position = "bottom", legend.text = element_text(size = 8))

  cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3)
}

## ─────────────────────────────────────────────────────────────
## 5. Sensitivity Analysis
## ─────────────────────────────────────────────────────────────
sensitivity_analysis <- function() {
  mod <- sma_model()

  # Vary SMN2 copy number (1–4)
  sa_smn2 <- lapply(1:4, function(copies) {
    m <- mrgsolve::param(mod, list(SMN2_copies = copies, k_MN_death = 0.003))
    ev_null <- mrgsolve::ev(time = 9999, amt = 0, cmt = 1)
    df <- as.data.frame(mrgsolve::mrgsim(m, ev_null, delta = 1, end = 365))
    df$SMN2_copies <- copies
    df
  }) |> bind_rows()

  # Vary nusinersen EC50
  ec50_vals <- c(2.5, 5.0, 10.0, 20.0)
  sa_ec50 <- lapply(ec50_vals, function(ec50) {
    m <- mrgsolve::param(mod, list(EC50_NUS = ec50, SMN2_copies = 2))
    df <- as.data.frame(mrgsolve::mrgsim(m, ev_nusinersen(), delta = 1, end = 730))
    df$EC50_NUS <- ec50
    df
  }) |> bind_rows()

  list(smn2_copies = sa_smn2, ec50_nusinersen = sa_ec50)
}

## ─────────────────────────────────────────────────────────────
## 6. Virtual Population (PopPK)
## ─────────────────────────────────────────────────────────────
virtual_population <- function(n = 100, seed = 42) {
  set.seed(seed)
  mod <- sma_model()

  # Sample interindividual variability
  idata <- data.frame(
    ID           = 1:n,
    SMN2_copies  = sample(2:3, n, replace = TRUE, prob = c(0.6, 0.4)),
    k_MN_death   = rlnorm(n, log(0.002), 0.3),
    EC50_NUS     = rlnorm(n, log(5.0), 0.4),
    Emax_NUS     = rnorm(n, 0.60, 0.05) |> pmin(0.85) |> pmax(0.35),
    k_prot_deg   = rlnorm(n, log(0.231), 0.2)
  )

  ev_nus <- ev_nusinersen()
  out <- mrgsolve::mrgsim(mod, idata = idata, events = ev_nus,
                          delta = 7, end = 730, obsonly = TRUE)
  as.data.frame(out)
}

## ─────────────────────────────────────────────────────────────
## 7. Main: Run Everything
## ─────────────────────────────────────────────────────────────
if (FALSE) {  # set TRUE to run interactively
  library(mrgsolve); library(dplyr); library(ggplot2); library(cowplot)

  cat("Building SMA QSP model...\n")
  mod <- sma_model()
  print(mod)

  cat("Running 6 treatment scenarios...\n")
  df_all <- run_scenarios()

  cat("Plotting results...\n")
  p <- plot_results(df_all)
  print(p)

  cat("Running sensitivity analysis...\n")
  sa <- sensitivity_analysis()
  print(sa$smn2_copies |> filter(time == 365) |>
    group_by(SMN2_copies) |> summarise(SMN_prot = last(SMN_prot), MN = last(MN_pool)))

  cat("Generating virtual population (n=100)...\n")
  vpc <- virtual_population(100)
  cat("VPC summary:\n")
  print(vpc |> filter(time == 365) |>
    group_by(ID) |> slice(n()) |>
    summarise(across(c(MN_fraction, CMAP, CHOP_INTEND, FVC),
                     list(median = median, q5 = ~quantile(., 0.05), q95 = ~quantile(., 0.95)))))
}
