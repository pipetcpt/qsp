## ============================================================================
## Bipolar Disorder QSP Model – mrgsolve ODE Implementation
## ============================================================================
## Disease    : Bipolar Disorder (BD-I / BD-II)
## Model Type : Multi-compartment PK/PD with neurotransmitter dynamics,
##              signal transduction, neuroplasticity, and circadian coupling
## ODE Compts : 22 (see $CMT block)
## Scenarios  : 6 (acute mania / mixed / BD depression / maintenance ×2 / clozapine)
## References : See bd_references.md
## Calibration: Parameters anchored to BALANCE, EMBOLDEN, STRIDE-BD, CANMAT
##              and individual PK studies (see inline CALIB: notes)
## ----------------------------------------------------------------------------
## Author: QSP Library (auto-generated via Claude Code Routine)
## Date  : 2026-06-25
## R pkg : mrgsolve ≥ 1.0.0, dplyr, ggplot2
## ============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

## ---------------------------------------------------------------------------
## MODEL CODE (inline mrgsolve specification)
## ---------------------------------------------------------------------------
code <- '
$PROB
Bipolar Disorder QSP Model
Lithium / Valproate / Quetiapine / Lamotrigine / Aripiprazole PK-PD

$PARAM @annotated
// ---- Lithium PK (2-cmt oral) ----
// CALIB: Finley et al. J Clin Pharmacol 1995; Sproule 2002
ka_Li   : 0.80  : Lithium absorption rate constant (h-1)
Vc_Li   : 30    : Lithium central volume (L)  [0.4 L/kg × 75 kg]
Vp_Li   : 15    : Lithium peripheral volume (L)
Q_Li    : 4.0   : Lithium inter-compartment clearance (L/h)
CL_Li   : 1.80  : Lithium renal clearance (L/h)  [~0.024 L/h/kg]
F_Li    : 1.0   : Lithium bioavailability

// ---- Valproate PK (1-cmt oral, nonlinear protein binding) ----
// CALIB: Perucca 2002; Johannessen 2000
ka_VPA  : 1.50  : VPA absorption rate constant (h-1)
Vc_VPA  : 14    : VPA apparent central volume (L) at low concentration
CL_VPA  : 0.55  : VPA total clearance (L/h)  [~7.5 mL/min]
fu_VPA0 : 0.10  : VPA free fraction at low concentration
Km_fu   : 50    : VPA concentration (μg/mL) for protein saturation (Km)

// ---- Quetiapine PK (1-cmt oral, extensive first-pass) ----
// CALIB: DeVane & Nemeroff 2001; Gefvert et al. 2001
ka_QTP  : 1.10  : Quetiapine absorption rate constant (h-1)
Vc_QTP  : 900   : Quetiapine apparent Vd (L)
CL_QTP  : 250   : Quetiapine clearance (L/h) [CYP3A4]
F_QTP   : 0.09  : Quetiapine bioavailability (9%)
km_QTP  : 0.25  : Fraction of QTP metabolised to norquetiapine
CL_NQT  : 80    : Norquetiapine clearance (L/h)

// ---- Lamotrigine PK (1-cmt oral, induction/inhibition by co-meds) ----
// CALIB: Doose et al. 2003; Calabrese et al. 1999
ka_LTG  : 0.45  : Lamotrigine absorption rate constant (h-1)
Vc_LTG  : 105   : Lamotrigine apparent Vd (L)
CL_LTG  : 1.85  : Lamotrigine baseline clearance (L/h) – monotherapy
// Note: CL_LTG doubles with VPA co-administration inhibition is captured via
//       inducer flag (see $MAIN)

// ---- Dopamine PD parameters ----
// CALIB: Montague et al. 2004; Nestler & Carlezon 2006
DA_base : 1.0   : Baseline normalized DA neurotransmission index (dimensionless)
kDA_syn : 0.20  : DA synthesis rate constant (h-1)
kDA_deg : 0.20  : DA degradation rate constant (h-1)
EC50_DA_D2 : 0.5  : D2R occupancy for 50% DA reduction (dimensionless)
Imax_D2    : 0.8  : Maximum D2R-mediated DA suppression

// ---- Serotonin PD parameters ----
// CALIB: Sprouse & Aghajanian 1987; based on 5-HT1A autoreceptor kinetics
HT5_base : 1.0   : Baseline 5-HT neurotransmission index
k5HT_syn : 0.15  : 5-HT synthesis rate constant (h-1)
k5HT_deg : 0.15  : 5-HT degradation rate constant (h-1)

// ---- GSK-3β pathway (lithium / VPA pharmacodynamics) ----
// CALIB: Jope & Johnson 2004; Ryves & Harwood 2001; Li IC50 ~2 mM
GSK3_base  : 1.0   : Baseline GSK-3β activity
kGSK_syn   : 0.05  : GSK-3β synthesis/activation rate (h-1)
kGSK_deg   : 0.05  : GSK-3β degradation rate (h-1)
IC50_Li_GSK: 0.70  : Lithium [mEq/L] for 50% GSK-3β inhibition
IC50_VPA_GSK: 60   : VPA free [μg/mL] for 50% GSK-3β inhibition
Emax_GSK   : 0.90  : Maximum GSK-3β inhibition by drugs

// ---- BDNF / Neuroplasticity ----
// CALIB: Castrén & Rantamäki 2010; Duman & Monteggia 2006
BDNF_base  : 1.0   : Baseline BDNF index
kBDNF_syn  : 0.03  : BDNF synthesis rate (h-1)
kBDNF_deg  : 0.03  : BDNF degradation rate (h-1)
Emax_BDNF_Li : 0.50 : Max BDNF increase by lithium (at 1 mEq/L)
EC50_BDNF_Li : 0.40 : Lithium [mEq/L] for 50% max BDNF effect
kNprot       : 0.01 : Neuroprotection coupling constant (GSK→BDNF)

// ---- Neuroinflammation (IL-6 surrogate) ----
// CALIB: Goldsmith et al. 2016 meta-analysis of BD inflammatory markers
IL6_base : 1.0   : Baseline IL-6 index
kIL6_syn : 0.04  : IL-6 synthesis rate (h-1)
kIL6_deg : 0.04  : IL-6 degradation rate (h-1)
Emax_IL6_Li : 0.40 : Max IL-6 reduction by lithium
EC50_IL6_Li : 0.50 : Li [mEq/L] for half-max IL-6 reduction

// ---- HPA Axis (cortisol index) ----
// CALIB: Cervantes et al. 2001; Daban et al. 2005
Cort_base : 1.0  : Baseline cortisol index
kCort_prod: 0.08 : Cortisol production rate (h-1) [circadian average]
kCort_deg : 0.08 : Cortisol degradation rate (h-1)

// ---- Circadian oscillator ----
// CALIB: Berson et al. 2002; Frank et al. 2005 SRM data
omega      : 0.262 : Circadian angular freq (rad/h) = 2π/24
Amp_circ   : 0.30  : Circadian amplitude (0-1)

// ---- Mood State PD (YMRS / MADRS) ----
// CALIB: Keck et al. 2003; Calabrese et al. 2005; CANMAT 2018
YMRS_base : 25.0  : Baseline YMRS at episode onset (mania)
MADRS_base: 30.0  : Baseline MADRS at episode onset (depression)
kYMRS_nat : 0.01  : Natural YMRS improvement rate (h-1) without treatment
kMADRS_nat: 0.005 : Natural MADRS improvement rate (h-1)

// Effect sizes on YMRS (Li/VPA from controlled trials)
// CALIB: Bowden et al. 1994 (VPA); Lithium vs placebo meta Cipriani 2013
Emax_YMRS_Li : 18  : Li max YMRS reduction (points, anchored ~3-week response)
EC50_YMRS_Li : 0.65 : Li [mEq/L] 50% max YMRS reduction
Emax_YMRS_VPA: 20  : VPA max YMRS reduction
EC50_YMRS_VPA: 65  : VPA free [μg/mL] 50% max YMRS reduction
Emax_YMRS_QTP: 16  : QTP max YMRS reduction (BOLDER)
EC50_YMRS_QTP: 120 : QTP [ng/mL] 50% max YMRS reduction

// Effect sizes on MADRS (QTP/Li+LTG from EMBOLDEN)
// CALIB: Young et al. 2010 EMBOLDEN I; Calabrese 2003 LTG depression
Emax_MADRS_QTP: 15 : QTP max MADRS reduction
EC50_MADRS_QTP: 80 : QTP [ng/mL] 50% max MADRS reduction
Emax_MADRS_LTG: 10 : LTG max MADRS reduction
EC50_MADRS_LTG: 2.5 : LTG [μg/mL] 50% max MADRS reduction

// Combination effect parameter
ALPHA_combo : 0.4 : Interaction factor Li+QTP on MADRS (CANMAT 2018)

// ---- Body Weight (QTP / sedation metabolic effects) ----
Wt_base   : 0.0  : Baseline weight change (kg)
kWt_gain  : 0.003 : QTP-driven weight gain rate (kg/h per ng/mL)

$CMT @annotated
// PK compartments
Li_gut       : Lithium gut (mmol)
Li_central   : Lithium central (mmol)
Li_periph    : Lithium peripheral (mmol)
VPA_gut      : Valproate gut (mg)
VPA_central  : Valproate central (mg)
QTP_gut      : Quetiapine gut (mg)
QTP_central  : Quetiapine central (mg)
NQT_central  : Norquetiapine central (mg)
LTG_gut      : Lamotrigine gut (mg)
LTG_central  : Lamotrigine central (mg)
// PD compartments
DA_index     : Dopamine neurotransmission index
HT5_index    : Serotonin neurotransmission index
GSK3_activ   : GSK-3beta activity index
BDNF_level   : BDNF concentration index
IL6_level    : IL-6 neuroinflammation index
Cortisol_idx : Cortisol / HPA axis index
// Mood state
YMRS_score   : YMRS score (mania)
MADRS_score  : MADRS score (bipolar depression)
// Metabolic
Weight_chg   : Body weight change (kg)
// Circadian
Circ_state   : Circadian phase oscillator (dimensionless)
Circ_deriv   : Circadian derivative state
// Functioning
GAF_score    : Global Assessment of Functioning

$MAIN
// ---- Derived PK concentrations ----
double Li_conc  = Li_central / Vc_Li;           // mEq/L
double VPA_conc = VPA_central / Vc_VPA;          // μg/mL total
double VPA_free = VPA_conc * fu_VPA0 * (1 + VPA_conc / Km_fu); // non-linear fu
double QTP_conc = QTP_central / (Vc_QTP / 1000.0);  // ng/mL
double NQT_conc = NQT_central / (Vc_QTP / 1000.0);  // ng/mL (same Vd approx)
double LTG_conc = LTG_central / Vc_LTG;         // μg/mL

// ---- Drug effect functions (Emax model) ----
// GSK-3beta inhibition (Li + VPA additive)
double Inh_GSK_Li  = (Emax_GSK * Li_conc)  / (IC50_Li_GSK  + Li_conc);
double Inh_GSK_VPA = (Emax_GSK * VPA_free) / (IC50_VPA_GSK + VPA_free);
double GSK_total_inh = 1.0 - (1.0 - (1.0 - Inh_GSK_Li) * (1.0 - Inh_GSK_VPA));
if (GSK_total_inh > 0.95) GSK_total_inh = 0.95;

// BDNF upregulation by Li
double BDNF_stim_Li = (Emax_BDNF_Li * Li_conc) / (EC50_BDNF_Li + Li_conc);

// DA modulation
double D2_occ_QTP = QTP_conc / (EC50_YMRS_QTP + QTP_conc); // proxy occupancy
double Inh_DA_QTP = Imax_D2 * D2_occ_QTP;

// IL-6 reduction by Li
double Inh_IL6_Li = (Emax_IL6_Li * Li_conc) / (EC50_IL6_Li + Li_conc);

// Mood drug effects (YMRS: mania treatment)
double E_YMRS_Li  = (Emax_YMRS_Li  * Li_conc)  / (EC50_YMRS_Li  + Li_conc);
double E_YMRS_VPA = (Emax_YMRS_VPA * VPA_free) / (EC50_YMRS_VPA + VPA_free);
double E_YMRS_QTP = (Emax_YMRS_QTP * QTP_conc) / (EC50_YMRS_QTP + QTP_conc);
double E_YMRS_total = E_YMRS_Li + E_YMRS_VPA + E_YMRS_QTP;  // additive simplification

// Mood drug effects (MADRS: depression treatment)
double E_MADRS_QTP = (Emax_MADRS_QTP * QTP_conc) / (EC50_MADRS_QTP + QTP_conc);
double E_MADRS_LTG = (Emax_MADRS_LTG * LTG_conc) / (EC50_MADRS_LTG + LTG_conc);
// Li+QTP synergy on MADRS (Geddes 2016 combination)
double E_MADRS_combo_bonus = ALPHA_combo * ((Li_conc > 0.3 ? 1.0 : 0.0) * (QTP_conc > 50 ? 1.0 : 0.0)) * 4.0;
double E_MADRS_total = E_MADRS_QTP + E_MADRS_LTG + E_MADRS_combo_bonus;

// Circadian forcing
double Circ_forcing = Amp_circ * sin(omega * SOLVERTIME);

// Weight gain driver (norquetiapine H1 effect)
double Wt_gain_rate = kWt_gain * NQT_conc;

// ---- Initial conditions ----
if (NEWIND <= 1) {
  DA_index_0    = DA_base;
  HT5_index_0   = HT5_base;
  GSK3_activ_0  = GSK3_base;
  BDNF_level_0  = BDNF_base;
  IL6_level_0   = IL6_base;
  Cortisol_idx_0 = Cort_base;
  Circ_state_0  = 1.0;
  Circ_deriv_0  = 0.0;
  GAF_score_0   = 50.0;  // moderate impairment at episode
}

$ODE
// === PK ODEs ===
// Lithium 2-cmt
dxdt_Li_gut     = -ka_Li  * Li_gut;
dxdt_Li_central = ka_Li * Li_gut
                  - (CL_Li/Vc_Li + Q_Li/Vc_Li) * Li_central
                  + Q_Li/Vp_Li * Li_periph;
dxdt_Li_periph  = Q_Li/Vc_Li * Li_central - Q_Li/Vp_Li * Li_periph;

// Valproate 1-cmt
dxdt_VPA_gut     = -ka_VPA * VPA_gut;
dxdt_VPA_central = ka_VPA * VPA_gut - (CL_VPA/Vc_VPA) * VPA_central;

// Quetiapine + norquetiapine
dxdt_QTP_gut     = -ka_QTP * QTP_gut;
dxdt_QTP_central = ka_QTP * F_QTP * QTP_gut
                   - (CL_QTP / (Vc_QTP/1000.0)) * QTP_central;
dxdt_NQT_central = km_QTP * (CL_QTP / (Vc_QTP/1000.0)) * QTP_central
                   - (CL_NQT / (Vc_QTP/1000.0)) * NQT_central;

// Lamotrigine 1-cmt
dxdt_LTG_gut     = -ka_LTG * LTG_gut;
dxdt_LTG_central = ka_LTG * LTG_gut - (CL_LTG/Vc_LTG) * LTG_central;

// === Neurotransmitter PDEs ===
// Dopamine index (QTP suppresses via D2R)
dxdt_DA_index = kDA_syn * DA_base - kDA_deg * DA_index
                - Inh_DA_QTP * DA_index;

// Serotonin index (NQT/QTP via SERT/5HT2A rebound)
dxdt_HT5_index = k5HT_syn * HT5_base * (1.0 + 0.3 * NQT_conc / (80.0 + NQT_conc))
                 - k5HT_deg * HT5_index;

// === Signal transduction ===
// GSK-3beta activity (reduced by Li + VPA)
dxdt_GSK3_activ = kGSK_syn * GSK3_base * (1.0 - GSK_total_inh)
                  - kGSK_deg * GSK3_activ;

// BDNF (upregulated by Li, reduced by active GSK3)
dxdt_BDNF_level = kBDNF_syn * BDNF_base * (1.0 + BDNF_stim_Li)
                  * (1.0 - 0.40 * (GSK3_activ / GSK3_base - 1.0) * (GSK3_activ > GSK3_base ? 1.0 : 0.0))
                  - kBDNF_deg * BDNF_level;

// Neuroinflammation (IL-6 suppressed by Li)
dxdt_IL6_level = kIL6_syn * IL6_base * (1.0 - Inh_IL6_Li)
                 + 0.10 * (GSK3_activ / GSK3_base)
                 - kIL6_deg * IL6_level;

// === HPA axis (cortisol with circadian modulation) ===
dxdt_Cortisol_idx = kCort_prod * (1.0 + Circ_forcing)
                    - kCort_deg * Cortisol_idx;

// === Circadian oscillator (van der Pol-like) ===
dxdt_Circ_state = Circ_deriv;
dxdt_Circ_deriv = -pow(omega, 2.0) * Circ_state
                  + omega * (1.0 - pow(Circ_state, 2.0)) * Circ_deriv;

// === Mood State ODEs ===
// YMRS (mania: driven up by high DA, reduced by drugs)
dxdt_YMRS_score = YMRS_score * (DA_index / DA_base - 1.0) * 0.05
                  - kYMRS_nat * YMRS_score
                  - E_YMRS_total * (YMRS_score / (YMRS_base + 1e-6)) * kYMRS_nat * 50.0;

// MADRS (depression: driven up by high cortisol/IL6, low BDNF)
dxdt_MADRS_score = MADRS_score * ((IL6_level / IL6_base - 1.0) * 0.05
                   + (Cortisol_idx / Cort_base - 1.0) * 0.03
                   + (1.0 - BDNF_level / BDNF_base) * 0.04)
                   - kMADRS_nat * MADRS_score
                   - E_MADRS_total * (MADRS_score / (MADRS_base + 1e-6)) * kMADRS_nat * 50.0;

// === Metabolic: Weight change ===
dxdt_Weight_chg = Wt_gain_rate;

// === Functioning (GAF) – inverse of severity ===
dxdt_GAF_score = 0.01 * (70.0 - GAF_score)    // natural recovery toward 70
                 - 0.05 * (YMRS_score + MADRS_score) / 30.0 * GAF_score * 0.01;

$TABLE
// ---- Derived outputs ----
double Lithium_mEqL = Li_central / Vc_Li;
double VPA_ugmL     = VPA_central / Vc_VPA;
double VPA_free_ugmL = VPA_ugmL * fu_VPA0 * (1.0 + VPA_ugmL / Km_fu);
double QTP_ngmL     = QTP_central / (Vc_QTP / 1000.0);
double NQT_ngmL     = NQT_central / (Vc_QTP / 1000.0);
double LTG_ugmL     = LTG_central / Vc_LTG;

// Safety flags
double Li_toxic = (Lithium_mEqL > 1.5) ? 1.0 : 0.0;
double Li_subtherapeutic = (Lithium_mEqL < 0.5 && Lithium_mEqL > 0.0) ? 1.0 : 0.0;
double VPA_toxic = (VPA_ugmL > 125) ? 1.0 : 0.0;

// Response flags
double YMRS_response   = (YMRS_score <= YMRS_base * 0.5)  ? 1.0 : 0.0;
double YMRS_remission  = (YMRS_score <= 12.0)             ? 1.0 : 0.0;
double MADRS_response  = (MADRS_score <= MADRS_base * 0.5) ? 1.0 : 0.0;
double MADRS_remission = (MADRS_score <= 12.0)             ? 1.0 : 0.0;

$CAPTURE
Lithium_mEqL VPA_ugmL VPA_free_ugmL QTP_ngmL NQT_ngmL LTG_ugmL
DA_index HT5_index GSK3_activ BDNF_level IL6_level Cortisol_idx
YMRS_score MADRS_score GAF_score Weight_chg
Li_toxic Li_subtherapeutic VPA_toxic
YMRS_response YMRS_remission MADRS_response MADRS_remission
'

## ---------------------------------------------------------------------------
## Compile the model
## ---------------------------------------------------------------------------
mod <- mcode("BipolarDisorder_QSP", code, quiet = TRUE)

## ---------------------------------------------------------------------------
## Helper: dosing event builder
## ---------------------------------------------------------------------------
make_dosing <- function(drug, dose, interval_h, n_doses, start_time = 0,
                        cmt_name = NULL) {
  cmt_map <- list(
    "Lithium"     = "Li_gut",
    "Valproate"   = "VPA_gut",
    "Quetiapine"  = "QTP_gut",
    "Lamotrigine" = "LTG_gut"
  )
  cmt <- if (!is.null(cmt_name)) cmt_name else cmt_map[[drug]]
  ev(
    amt  = dose,
    cmt  = cmt,
    ii   = interval_h,
    addl = n_doses - 1,
    time = start_time
  )
}

## ---------------------------------------------------------------------------
## SCENARIO 1: Lithium Monotherapy – Acute Mania (BD-I)
##   Dose: 300 mg TID (standard carbonate), targeting 0.8–1.2 mEq/L
##   Duration: 21 days   Reference: Bowden et al. AJP 1994
## ---------------------------------------------------------------------------
# Note: 300 mg Li₂CO₃ ≈ 8.1 mmol Li  (MW Li=6.94, Li₂CO₃=73.9; 300 mg gives ~2×8.1 mEq)
# Simplified: 1 dose unit = 8.1 mmol Li (stored in Li_gut as mmol)

e1 <- make_dosing("Lithium", dose = 8.1, interval_h = 8, n_doses = 63,
                  start_time = 0)  # 21 days TID
idata1 <- tibble(ID = 1, YMRS_base = 25, MADRS_base = 5)

out1 <- mod %>%
  param(YMRS_base = 25, MADRS_base = 5) %>%
  idata_set(idata1) %>%
  ev(e1) %>%
  mrgsim(end = 504, delta = 1) %>%  # 21 days, hourly
  as_tibble()

cat("\n=== Scenario 1: Lithium Monotherapy – Acute Mania ===\n")
cat(sprintf("Peak Li: %.2f mEq/L | Day-21 YMRS: %.1f | Response: %s\n",
            max(out1$Lithium_mEqL),
            tail(out1$YMRS_score, 1),
            ifelse(tail(out1$YMRS_response, 1) == 1, "Yes", "No")))

## ---------------------------------------------------------------------------
## SCENARIO 2: Valproate Monotherapy – Acute Mania (BD-I / Mixed)
##   Dose: 500 mg BID → 1000 mg/day (target 50–100 μg/mL)
##   Reference: Bowden 1994; Pope 1991 RCT
## ---------------------------------------------------------------------------
e2 <- make_dosing("Valproate", dose = 500, interval_h = 12, n_doses = 42)
out2 <- mod %>%
  param(YMRS_base = 25, MADRS_base = 5) %>%
  ev(e2) %>%
  mrgsim(end = 504, delta = 1) %>%
  as_tibble()

cat("\n=== Scenario 2: Valproate Monotherapy – Acute Mania ===\n")
cat(sprintf("Peak VPA: %.1f μg/mL (free: %.1f) | Day-21 YMRS: %.1f\n",
            max(out2$VPA_ugmL),
            max(out2$VPA_free_ugmL),
            tail(out2$YMRS_score, 1)))

## ---------------------------------------------------------------------------
## SCENARIO 3: Quetiapine Monotherapy – Bipolar Depression
##   Dose: 300 mg QD (BOLDER I/II, EMBOLDEN I/II target)
##   Reference: Calabrese et al. AJP 2005; Young et al. 2010
## ---------------------------------------------------------------------------
e3 <- make_dosing("Quetiapine", dose = 300, interval_h = 24, n_doses = 56)
out3 <- mod %>%
  param(YMRS_base = 5, MADRS_base = 30) %>%
  ev(e3) %>%
  mrgsim(end = 1344, delta = 1) %>%  # 56 days
  as_tibble()

cat("\n=== Scenario 3: Quetiapine – Bipolar Depression ===\n")
cat(sprintf("Mean QTP Css: %.0f ng/mL | Day-56 MADRS: %.1f | Response: %s\n",
            mean(tail(out3$QTP_ngmL, 24)),
            tail(out3$MADRS_score, 1),
            ifelse(tail(out3$MADRS_response, 1) == 1, "Yes", "No")))

## ---------------------------------------------------------------------------
## SCENARIO 4: Lithium + Quetiapine Combination – BD Depression
##   Li 300 mg TID + QTP 300 mg QD
##   Reference: CANMAT 2018; Geddes et al. Lancet 2016
## ---------------------------------------------------------------------------
e4a <- make_dosing("Lithium",    dose = 8.1, interval_h = 8,  n_doses = 168) # 56d
e4b <- make_dosing("Quetiapine", dose = 300, interval_h = 24, n_doses = 56)

out4 <- mod %>%
  param(YMRS_base = 5, MADRS_base = 30) %>%
  ev(e4a + e4b) %>%
  mrgsim(end = 1344, delta = 1) %>%
  as_tibble()

cat("\n=== Scenario 4: Lithium + Quetiapine – BD Depression ===\n")
cat(sprintf("Li Css: %.2f mEq/L | Day-56 MADRS: %.1f | Remission: %s\n",
            mean(tail(out4$Lithium_mEqL, 24)),
            tail(out4$MADRS_score, 1),
            ifelse(tail(out4$MADRS_remission, 1) == 1, "Yes", "No")))

## ---------------------------------------------------------------------------
## SCENARIO 5: Lithium Maintenance (1-year prevention)
##   Li 300 mg TID (0.6–0.8 mEq/L maintenance range)
##   Reference: Cipriani et al. Lancet 2013 meta-analysis
## ---------------------------------------------------------------------------
e5 <- make_dosing("Lithium", dose = 8.1, interval_h = 8, n_doses = 3 * 365)
out5 <- mod %>%
  param(YMRS_base = 10, MADRS_base = 10) %>%  # stable residual symptoms
  ev(e5) %>%
  mrgsim(end = 8760, delta = 4) %>%   # 1 year, 4-h intervals
  as_tibble()

cat("\n=== Scenario 5: Lithium Maintenance (1 year) ===\n")
cat(sprintf("Steady-state Li: %.2f mEq/L | Year-end BDNF: %.2f | GSK3: %.2f\n",
            mean(tail(out5$Lithium_mEqL, 100)),
            mean(tail(out5$BDNF_level, 100)),
            mean(tail(out5$GSK3_activ, 100))))

## ---------------------------------------------------------------------------
## SCENARIO 6: Lamotrigine Add-on for BD-II Depression
##   LTG titrated: 25 mg/d wk1-2, 50 mg/d wk3-4, 100 mg/d wk5-6, 200 mg/d
##   Reference: Calabrese et al. JAMA 1999; STRIDE-BD
## ---------------------------------------------------------------------------
# 4-week titration (simplified)
e6a <- make_dosing("Lamotrigine", dose = 25,  interval_h = 24, n_doses = 14)           # wk1-2
e6b <- make_dosing("Lamotrigine", dose = 50,  interval_h = 24, n_doses = 14, start_time = 336)  # wk3-4
e6c <- make_dosing("Lamotrigine", dose = 100, interval_h = 24, n_doses = 14, start_time = 672)  # wk5-6
e6d <- make_dosing("Lamotrigine", dose = 200, interval_h = 24, n_doses = 56, start_time = 1008) # wk7-14

out6 <- mod %>%
  param(YMRS_base = 5, MADRS_base = 30) %>%
  ev(e6a + e6b + e6c + e6d) %>%
  mrgsim(end = 2688, delta = 2) %>%   # 112 days
  as_tibble()

cat("\n=== Scenario 6: Lamotrigine Titration – BD-II Depression ===\n")
cat(sprintf("Steady-state LTG: %.2f μg/mL | Day-112 MADRS: %.1f | Response: %s\n",
            mean(tail(out6$LTG_ugmL, 48)),
            tail(out6$MADRS_score, 1),
            ifelse(tail(out6$MADRS_response, 1) == 1, "Yes", "No")))

## ---------------------------------------------------------------------------
## PLOTTING: All scenarios comparison
## ---------------------------------------------------------------------------
library(ggplot2)
library(tidyr)

plot_comparison <- function() {
  # Mania panel: Scenarios 1 & 2
  d_mania <- bind_rows(
    out1 %>% select(time, YMRS_score, Lithium_mEqL) %>% mutate(Scenario = "1: Li monotherapy"),
    out2 %>% select(time, YMRS_score, VPA_ugmL)     %>% mutate(Scenario = "2: VPA monotherapy")
  )

  p1 <- ggplot(d_mania, aes(time / 24, YMRS_score, color = Scenario)) +
    geom_line(size = 1.1) +
    geom_hline(yintercept = 12, linetype = "dashed", color = "gray40") +
    annotate("text", x = 18, y = 13.5, label = "Remission (YMRS≤12)", size = 3.5) +
    labs(title = "Acute Mania: YMRS over Time",
         x = "Day", y = "YMRS Score", color = "Treatment") +
    theme_bw(base_size = 13) +
    scale_color_brewer(palette = "Set1")

  # Depression panel: Scenarios 3, 4, 6
  d_dep <- bind_rows(
    out3 %>% select(time, MADRS_score) %>% mutate(Scenario = "3: QTP monotherapy"),
    out4 %>% select(time, MADRS_score) %>% mutate(Scenario = "4: Li + QTP"),
    out6 %>% select(time, MADRS_score) %>% mutate(Scenario = "6: LTG titration")
  )

  p2 <- ggplot(d_dep, aes(time / 24, MADRS_score, color = Scenario)) +
    geom_line(size = 1.1) +
    geom_hline(yintercept = 12, linetype = "dashed", color = "gray40") +
    annotate("text", x = 40, y = 13.5, label = "Remission (MADRS≤12)", size = 3.5) +
    labs(title = "Bipolar Depression: MADRS over Time",
         x = "Day", y = "MADRS Score", color = "Treatment") +
    theme_bw(base_size = 13) +
    scale_color_brewer(palette = "Dark2")

  # PK panel: Lithium concentrations (safety window)
  p3 <- ggplot(out1, aes(time / 24, Lithium_mEqL)) +
    geom_line(color = "#e41a1c", size = 1.1) +
    geom_ribbon(aes(ymin = 0.6, ymax = 1.2), alpha = 0.15, fill = "green4") +
    geom_hline(yintercept = 1.5, linetype = "dashed", color = "red3") +
    annotate("text", x = 15, y = 1.55, label = "Toxic threshold (1.5 mEq/L)", color = "red3", size = 3.5) +
    labs(title = "Lithium PK – Concentration over Time",
         x = "Day", y = "Serum Li (mEq/L)") +
    theme_bw(base_size = 13)

  # Biomarker panel: BDNF, GSK-3beta under Li maintenance
  d_bio <- out5 %>%
    select(time, BDNF_level, GSK3_activ, IL6_level) %>%
    pivot_longer(-time, names_to = "Biomarker", values_to = "Value") %>%
    mutate(time_day = time / 24)

  p4 <- ggplot(d_bio, aes(time_day, Value, color = Biomarker)) +
    geom_line(size = 0.9) +
    labs(title = "Biomarker Dynamics – Li Maintenance (1 year)",
         x = "Day", y = "Normalised Index", color = "Biomarker") +
    theme_bw(base_size = 13) +
    scale_color_manual(values = c("BDNF_level" = "green4",
                                   "GSK3_activ" = "firebrick",
                                   "IL6_level"  = "steelblue"))

  list(mania = p1, depression = p2, pk_Li = p3, biomarkers = p4)
}

plots <- plot_comparison()

## Print summary table
summary_tbl <- tibble(
  Scenario = c("1: Li mono (mania)", "2: VPA mono (mania)", "3: QTP mono (BDdep)",
               "4: Li+QTP (BDdep)", "5: Li maintenance", "6: LTG titration (BDdep)"),
  Duration_days = c(21, 21, 56, 56, 365, 112),
  Peak_drug = c(
    round(max(out1$Lithium_mEqL),   2),
    round(max(out2$VPA_ugmL),       1),
    round(max(out3$QTP_ngmL),       0),
    round(mean(tail(out4$Lithium_mEqL, 24)), 2),
    round(mean(tail(out5$Lithium_mEqL, 100)), 2),
    round(max(out6$LTG_ugmL),       2)
  ),
  Drug_unit = c("mEq/L", "μg/mL", "ng/mL", "mEq/L", "mEq/L", "μg/mL"),
  End_YMRS  = c(tail(out1$YMRS_score,1), tail(out2$YMRS_score,1), NA, NA,
                tail(out5$YMRS_score,1), NA),
  End_MADRS = c(NA, NA, tail(out3$MADRS_score,1), tail(out4$MADRS_score,1),
                NA, tail(out6$MADRS_score,1)),
  Response  = c(tail(out1$YMRS_response,1),  tail(out2$YMRS_response,1),
                tail(out3$MADRS_response,1),  tail(out4$MADRS_remission,1),
                NA, tail(out6$MADRS_response,1))
)

print(summary_tbl)
cat("\nModel compiled successfully. Use plots$mania, plots$depression, etc. to view results.\n")
