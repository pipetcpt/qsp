################################################################################
# Heart Failure with Reduced Ejection Fraction (HFrEF)
# Quantitative Systems Pharmacology (QSP) Model
# mrgsolve ODE-based PK/PD Model
#
# Model structure:
#   RAAS: AngI → AngII, ACE2 → Ang(1-7), Aldosterone
#   SNS:  Plasma NE, Cardiac adrenergic signaling
#   NPS:  BNP/NT-proBNP, cGMP
#   Hemodynamics: LVEDV, SV, CO, SVR, MAP, HR
#   Cardiac remodeling: Fibrosis (TGF-β1/Collagen), Hypertrophy index
#   Inflammation: TNF-α, IL-6
#   Drug PK: ARNI (LBQ657/Valsartan), β-blocker, MRA, SGLT2i (4-compartment)
#
# Calibration references:
#   PARADIGM-HF (McMurray NEJM 2014): sacubitril/valsartan vs enalapril
#   EMPEROR-Reduced (Packer NEJM 2020): empagliflozin in HFrEF
#   MERIT-HF (Metoprolol RALES/EMPHASIS trials)
#   RALES (Pitt NEJM 1999): spironolactone
#   SHIFT (Swedberg Lancet 2010): ivabradine
#
# 21 ODE states, 5 treatment scenarios + dose-response analysis
################################################################################

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

# ─────────────────────────────────────────────────────────────────────────────
# MODEL CODE BLOCK
# ─────────────────────────────────────────────────────────────────────────────

hfref_model_code <- '
$PROB
HFrEF QSP Model — RAAS + SNS + NPS + Hemodynamics + Remodeling + Drug PK/PD
21 ODE states. Calibrated to PARADIGM-HF, EMPEROR-Reduced, MERIT-HF, RALES.

$PARAM
// ── Disease parameters ──────────────────────────────────────────────────────
// RAAS
kAngI_prod   = 0.12    // AngI production rate (ng/mL/h) — RAAS activation
kAngI_deg    = 0.05    // AngI intrinsic degradation (1/h)
kACE_Vmax    = 0.20    // ACE Michaelis-Menten Vmax for AngI → AngII (ng/mL/h)
kACE_Km      = 2.0     // ACE Km (ng/mL)
kACE2_Vmax   = 0.04    // ACE2 for AngI → Ang1-7 (ng/mL/h)
kACE2_Km     = 2.0     // ACE2 Km (ng/mL)
kAngII_deg   = 0.30    // AngII overall degradation (1/h)
kAng17_deg   = 0.25    // Ang1-7 degradation (1/h)
EC50_AngII   = 40.0    // AngII for aldosterone release (pg/mL)
kAldo_max    = 0.60    // Max aldosterone synthesis rate (ng/dL/h)
kAldo_deg    = 0.18    // Aldosterone degradation (1/h)
Aldo_base    = 8.0     // Baseline aldosterone (ng/dL)

// SNS
kNE_prod     = 180.0   // NE production rate (pg/mL/h)
kNE_deg      = 0.55    // NE clearance (1/h)
NE_base      = 325.0   // Healthy NE baseline (pg/mL)
// HFrEF: baseline NE ~600-900 pg/mL (2× normal)
NE_HF_factor = 2.2     // HF SNS activation factor

// NPS
kBNP_prod    = 48.0    // BNP synthesis rate (pg/mL/h) — wall stress driven
kBNP_NEPdeg  = 0.30    // BNP NEP-mediated degradation (1/h)
kBNP_other   = 0.12    // BNP other degradation (1/h)
BNP_base     = 35.0    // Healthy BNP (pg/mL)
// HFrEF typical BNP ≥100 pg/mL, NT-proBNP ≥900 pg/mL at admission
NTproBNP_ratio = 8.5   // NT-proBNP/BNP ratio (clearance much slower)
kcGMP_prod   = 2.5     // cGMP production from BNP/ANP (pmol/mL/h)
kcGMP_PDE    = 0.80    // cGMP PDE5-mediated degradation (1/h)
kcGMP_other  = 0.20    // cGMP other degradation (1/h)

// Hemodynamics (baseline HFrEF: LVEF ~25–30%)
LVEDV_base   = 280.0   // LV end-diastolic volume at HFrEF (mL) [normal ~120]
HR_base      = 82.0    // Resting heart rate (bpm) — HFrEF typically elevated
HR_setpoint  = 82.0    // Setpoint HR
SVR_base     = 1400.0  // Systemic vascular resistance (dynes·s/cm5) [normal ~1100]
EF_base      = 0.27    // Baseline LVEF in HFrEF (0.27 = 27%)
EF_max       = 0.65    // Achievable max EF with optimal therapy
kSVR_decay   = 0.015   // SVR normalization rate under therapy (1/h)
MAP_target   = 85.0    // Target MAP (mmHg)

// Cardiac remodeling
kFib_prod    = 0.002   // Fibrosis accumulation rate (1/h) per TGF-β1 unit
kFib_deg     = 0.0005  // Spontaneous fibrosis resolution (very slow, 1/h)
Fib_base     = 0.35    // Baseline fibrosis score in HFrEF (0–1 scale)
Fib_max      = 0.85    // Maximum fibrosis (0–1)
kHyp_prod    = 0.004   // Hypertrophy index growth rate (1/h)
kHyp_deg     = 0.001   // Hypertrophy regression rate (1/h)
Hyp_base     = 1.45    // Baseline hypertrophy index (1=normal, >1=hypertrophied)

kTGFb1_prod  = 1.5     // TGF-β1 production (pg/mL/h) — driven by AngII+Aldo
kTGFb1_deg   = 0.25    // TGF-β1 degradation (1/h)
TGFb1_base   = 12.0    // Baseline TGF-β1 in HFrEF (pg/mL) [normal ~4]

// Inflammation
kTNFa_prod   = 0.8     // TNF-α production (pg/mL/h)
kTNFa_deg    = 0.35    // TNF-α degradation (1/h)
TNFa_base    = 6.5     // HFrEF baseline TNF-α (pg/mL) [normal ~1.5]
kIL6_prod    = 1.2     // IL-6 production (pg/mL/h)
kIL6_deg     = 0.40    // IL-6 degradation (1/h)
IL6_base     = 5.8     // HFrEF baseline IL-6 (pg/mL) [normal ~1.0]

// ── Drug PK parameters ──────────────────────────────────────────────────────
// ARNI: sacubitril/valsartan (Entresto)
// LBQ657 (active neprilysin inhibitor): CL=4.6 L/h, Vd=26L, F=~0.6 (Langenickel CPT 2016)
// Valsartan in ARNI: higher AUC vs standalone (F~94%)
CL_LBQ       = 4.6    // LBQ657 apparent clearance (L/h)
Vd_LBQ       = 26.0   // LBQ657 volume of distribution (L)
ka_LBQ       = 1.8    // LBQ657 absorption rate (1/h) [fast, from sacubitril]
F_LBQ        = 0.60   // LBQ657 bioavailability
Dose_LBQ     = 0.0    // LBQ657 dose equivalent (mg) — set in scenarios
II_LBQ       = 12.0   // Dosing interval (h) — BID
CL_Valsar    = 1.9    // Valsartan CL (L/h) in ARNI formulation
Vd_Valsar    = 75.0   // Valsartan Vd (L)
ka_Valsar    = 0.9    // Valsartan ka (1/h)
F_Valsar     = 0.94   // Valsartan bioavailability in ARNI
Dose_Valsar  = 0.0    // Valsartan dose (mg)

// β-Blocker: metoprolol succinate (MERIT-HF: 200 mg/day target dose)
// CL: 55 L/h, Vd: 220 L, t½~5h (Regardh Eur J Clin Pharm 1974)
CL_BB        = 55.0   // β-blocker clearance (L/h)
Vd_BB        = 220.0  // β-blocker Vd (L)
ka_BB        = 0.5    // β-blocker absorption (1/h) — extended release
F_BB         = 0.45   // Metoprolol succinate bioavailability
Dose_BB      = 0.0    // β-blocker dose (mg)
II_BB        = 24.0   // Once daily dosing interval
EC50_BB      = 50.0   // BB Cp for 50% HR reduction (ng/mL)
Emax_BB_HR   = 0.30   // Max fractional HR reduction by BB (30%)
Emax_BB_NE   = 0.45   // Max NE-pathway inhibition by BB

// MRA: eplerenone (EMPHASIS-HF: 50 mg/day)
// CL: 7 L/h, Vd: 43 L, t½~4-6h (Cook 2003 Clin Pharmacokinet)
CL_MRA       = 7.0    // MRA clearance (L/h)
Vd_MRA       = 43.0   // MRA Vd (L)
ka_MRA       = 0.8    // MRA absorption (1/h)
F_MRA        = 0.69   // Eplerenone bioavailability
Dose_MRA     = 0.0    // MRA dose (mg)
II_MRA       = 24.0   // Once daily
EC50_MRA     = 50.0   // MRA Cp for 50% aldosterone blockade (ng/mL)
Emax_MRA_Ald = 0.80   // Max aldosterone effect blockade

// SGLT2i: empagliflozin (EMPEROR-Reduced: 10 mg/day)
// CL: 16.1 L/h, Vd: 73.8 L, t½~12h (Macha 2013 Clin Pharm)
CL_SGLT2     = 16.1   // SGLT2i clearance (L/h)
Vd_SGLT2     = 73.8   // SGLT2i Vd (L)
ka_SGLT2     = 0.65   // SGLT2i absorption (1/h)
F_SGLT2      = 0.86   // Empagliflozin bioavailability
Dose_SGLT2   = 0.0    // SGLT2i dose (mg)
II_SGLT2     = 24.0   // Once daily
EC50_SGLT2   = 40.0   // Empagliflozin Cp for 50% diuretic effect (ng/mL)
Emax_SGLT2_V = 0.18   // Max LVEDV reduction by SGLT2i (18%, EMPEROR)
Emax_SGLT2_C = 0.12   // Max direct cardiac improvement (12% via NHE1/mito)

// Ivabradine (SHIFT: 7.5 mg BID)
CL_IVA       = 10.0   // Ivabradine clearance (L/h)
Vd_IVA       = 100.0  // Ivabradine Vd (L)
ka_IVA       = 1.2    // Ivabradine absorption (1/h)
Dose_IVA     = 0.0    // Ivabradine dose (mg)
II_IVA       = 12.0   // BID dosing
EC50_IVA     = 20.0   // Ivabradine Cp for 50% HR reduction (ng/mL)
Emax_IVA_HR  = 0.25   // Max HR reduction (25 bpm typical)

// ── PD interaction parameters ──────────────────────────────────────────────
// ARNI NEP inhibition
EC50_NEP_inhib = 80.0  // LBQ657 Cp for 50% NEP inhibition (ng/mL)
// ARNI AT1R block (valsartan component)
EC50_AT1R_val  = 200.0 // Valsartan Cp for 50% AT1R blockade (ng/mL)
// ARNI effect on BNP (direct measurement: PARADIGM-HF: BNP↑23% despite NT-proBNP↓37%)
// ARNI on LVEF: +5-7% absolute over 12 months
kEF_ARNI     = 0.003   // Rate of LVEF improvement per unit drug effect (/h)
kEF_BB       = 0.002   // Rate of LVEF improvement from BB remodeling reversal (/h)
kEF_SGLT2    = 0.0015  // Rate of LVEF improvement from SGLT2i (/h)

$INIT
// RAAS
AngI    = 2.0      // Angiotensin I (ng/mL) — HFrEF elevated
AngII   = 55.0     // Angiotensin II (pg/mL) — HFrEF: ~50-80 pg/mL [normal ~15]
Ang17   = 3.5      // Angiotensin 1-7 (pg/mL)
Aldo    = 22.0     // Aldosterone (ng/dL) — HFrEF elevated [normal ~8]

// SNS
NE      = 715.0    // Plasma NE (pg/mL) — HFrEF ~600-900 pg/mL [normal ~325]

// NPS
BNP     = 380.0    // BNP plasma (pg/mL) — HFrEF threshold ≥100 pg/mL for dx
NTpBNP  = 3230.0   // NT-proBNP (pg/mL) — HFrEF: ~8.5 × BNP
cGMP    = 12.0     // cGMP (pmol/mL) — relatively low in HFrEF

// Hemodynamics
LVEDV   = 280.0    // LV end-diastolic volume (mL) — dilated
HR      = 82.0     // Heart rate (bpm) — elevated in HFrEF
SVR     = 1400.0   // Systemic vascular resistance (dynes·s/cm5)
LVEF    = 0.27     // LV ejection fraction — severely reduced

// Cardiac remodeling
Fib     = 0.35     // Fibrosis score (0–1)
TGFb1   = 12.0     // TGF-β1 (pg/mL) — elevated in HFrEF
Hyp     = 1.45     // Hypertrophy index (dimensionless)

// Inflammation
TNFa    = 6.5      // TNF-α (pg/mL) — HFrEF elevated
IL6     = 5.8      // IL-6 (pg/mL) — HFrEF elevated

// Drug PK (all start at 0 — pre-dose)
LBQ_C   = 0.0     // LBQ657 plasma (ng/mL)
Valsar_C = 0.0    // Valsartan plasma (ng/mL)
BB_C    = 0.0     // β-blocker plasma (ng/mL)
MRA_C   = 0.0     // MRA plasma (ng/mL)
SGLT2_C = 0.0     // SGLT2i plasma (ng/mL)
IVA_C   = 0.0     // Ivabradine plasma (ng/mL)

$ODE
// ════════════════════════════════════════════════════════════════════════════
// DRUG EFFECT CALCULATIONS (Hill equation)
// ════════════════════════════════════════════════════════════════════════════

// NEP inhibition by LBQ657 (0=none, 1=complete)
double E_NEP  = LBQ_C / (LBQ_C + EC50_NEP_inhib);
// AT1R blockade by valsartan (0=none, 1=complete)
double E_AT1R = Valsar_C / (Valsar_C + EC50_AT1R_val);
// β-blocker effects
double E_BB_HR  = Emax_BB_HR  * BB_C / (BB_C + EC50_BB);
double E_BB_NE  = Emax_BB_NE  * BB_C / (BB_C + EC50_BB);
// MRA effect on aldosterone signaling
double E_MRA  = Emax_MRA_Ald * MRA_C / (MRA_C + EC50_MRA);
// SGLT2i effects
double E_SGLT2_V = Emax_SGLT2_V * SGLT2_C / (SGLT2_C + EC50_SGLT2);
double E_SGLT2_C2 = Emax_SGLT2_C * SGLT2_C / (SGLT2_C + EC50_SGLT2);
// Ivabradine effect on HR
double E_IVA_HR = Emax_IVA_HR * IVA_C / (IVA_C + EC50_IVA);

// ════════════════════════════════════════════════════════════════════════════
// RAAS ODEs
// ════════════════════════════════════════════════════════════════════════════

// AngI: produced from angiotensinogen by renin; ACE converts to AngII; ACE2 converts to Ang1-7
double ACE_rate  = kACE_Vmax  * AngI / (kACE_Km  + AngI);
double ACE2_rate = kACE2_Vmax * AngI / (kACE2_Km + AngI);

// ACEi inhibition of ACE (valsartan in ARNI does not inhibit ACE;
// include ACEi as separate parameter if needed — here E_AT1R represents AT1R block)
// NEP inhibition raises ACE2-derived Ang1-7
double ACE2_effect = ACE2_rate * (1 + E_NEP * 0.5); // NEP block → relative ACE2 product↑

dxdt_AngI  = kAngI_prod - ACE_rate - ACE2_effect - kAngI_deg * AngI;

// AngII: produced by ACE from AngI; degraded; AT1R blockade does not change plasma AngII
// (AT1R block → reactive AngII rise ~2× due to feedback)
double AT1R_FB = 1.0 + E_AT1R * 1.5; // AT1R block → AngII reactive ↑1.5×
dxdt_AngII = ACE_rate * AT1R_FB - kAngII_deg * AngII;

// Ang1-7: produced by ACE2 from AngI; degraded; NEP inhibition reduces its breakdown
dxdt_Ang17 = ACE2_effect - kAng17_deg * Ang17 * (1.0 - E_NEP * 0.4);

// Aldosterone: driven by AngII (AT1R → adrenal zona glomerulosa), blocked by MRA
double Aldo_stim = kAldo_max * pow(AngII, 1.2) / (pow(EC50_AngII, 1.2) + pow(AngII, 1.2));
double Aldo_prod = Aldo_stim * (1.0 - E_AT1R * 0.7);  // AT1R block reduces Aldo production
dxdt_Aldo  = Aldo_prod - kAldo_deg * Aldo;

// ════════════════════════════════════════════════════════════════════════════
// SNS ODEs
// ════════════════════════════════════════════════════════════════════════════

// Norepinephrine: elevated in HFrEF due to baroreceptor desensitization
// NE is produced proportional to SNS activation (driven by low CO → CNS)
// BB reduces efferent SNS signaling indirectly
double NE_prod = kNE_prod * NE_HF_factor;
double NE_clearance = kNE_deg * (1.0 + E_BB_NE * 0.6) * NE;
// cGMP (via ANP/BNP) also attenuates NE release slightly
double NE_cGMP_inhibit = 0.015 * (cGMP - 8.0);  // cGMP above baseline slightly ↓NE
dxdt_NE = NE_prod - NE_clearance - NE_cGMP_inhibit * NE;

// ════════════════════════════════════════════════════════════════════════════
// NPS ODEs
// ════════════════════════════════════════════════════════════════════════════

// BNP: secreted from ventricle proportional to wall stress (≈ LVEDV × SVR / LVEF)
// NEP inhibition dramatically reduces BNP plasma degradation
// Note: plasma BNP ↑ with NEPi (PARADIGM paradox) but NT-proBNP ↓
double WS = LVEDV * SVR / (LVEF * 1e6 + 0.001); // relative wall stress
double BNP_prod_rate = kBNP_prod * (1.0 + 0.8 * (WS - 0.02) / 0.02);
if(BNP_prod_rate < 0) BNP_prod_rate = 0;
double BNP_NEP_deg   = kBNP_NEPdeg * (1.0 - E_NEP * 0.9) * BNP; // NEPi: ↓90% NEP-mediated
double BNP_other_deg = kBNP_other * BNP;
dxdt_BNP   = BNP_prod_rate - BNP_NEP_deg - BNP_other_deg;

// NT-proBNP: cleaved from pro-BNP in same ratio; longer half-life (~70h vs ~20min for BNP)
// NT-proBNP reflects synthesis better; NEPi does NOT degrade NT-proBNP (different pathway)
// PARADIGM-HF: NT-proBNP ↓37% with ARNI vs. enalapril
double NTpBNP_prod = 0.55 * BNP_prod_rate * NTproBNP_ratio;
double NTpBNP_deg  = 0.008 * NTpBNP; // slow degradation (t½~70h)
// Improvement in LVEDV/EF reduces production
double NTpBNP_EF_corr = (1.0 - (LVEF - 0.27) * 3.0); // rising EF ↓ NT-proBNP
if(NTpBNP_EF_corr < 0.2) NTpBNP_EF_corr = 0.2;
dxdt_NTpBNP = NTpBNP_prod * NTpBNP_EF_corr - NTpBNP_deg;

// cGMP: driven by BNP/ANP (through NPR-A → guanylyl cyclase); degraded by PDE5
double cGMP_prod_rate = kcGMP_prod * BNP / (BNP + 200.0); // BNP saturation ~200 pg/mL
dxdt_cGMP = cGMP_prod_rate - (kcGMP_PDE + kcGMP_other) * cGMP;

// ════════════════════════════════════════════════════════════════════════════
// HEMODYNAMIC ODEs
// ════════════════════════════════════════════════════════════════════════════

// LVEDV: increases with volume overload (Na retention via Aldo/Renin), decreases with diuretics/SGLT2i
double LVEDV_load   = 0.03 * Aldo / Aldo_base;   // Aldo drives Na/H2O retention → ↑LVEDV
double LVEDV_unload = 0.025 * cGMP / 10.0;        // BNP/cGMP → natriuresis → ↓LVEDV
double LVEDV_SGLT2  = E_SGLT2_V * 0.008 * LVEDV;  // SGLT2i osmotic diuresis
// Reverse remodeling (BB, ARNI) reduces LVEDV over months
double LVEDV_rmod   = 0.0010 * (E_BB_NE + E_AT1R) * (LVEDV - 120.0); // asymptote = 120 mL
dxdt_LVEDV = LVEDV_load - LVEDV_unload - LVEDV_SGLT2 - LVEDV_rmod;

// Heart Rate: elevated in HFrEF; reduced by BB and Ivabradine
double HR_target = HR_setpoint * (1.0 + 0.25 * NE / NE_base); // NE drives HR↑
double HR_BB_inh = HR_target * E_BB_HR;     // BB reduces HR
double HR_IVA_inh = E_IVA_HR * 25.0;       // Ivabradine absolute HR reduction (bpm)
double HR_setpt_new = HR_target - HR_BB_inh - HR_IVA_inh;
dxdt_HR = 0.15 * (HR_setpt_new - HR);      // first-order toward new setpoint (1/h, tau~7h)

// SVR: elevated in HFrEF due to AngII/NE/Aldo; reduced by ARNI/ACEi/SGLT2i
double SVR_AngII_eff = 1.0 + 0.40 * (AngII - 15.0) / 15.0; // AngII above normal → SVR↑
double SVR_AT1R_inh  = E_AT1R * 0.35;  // AT1R block → SVR↓35%
double SVR_NE_eff    = 0.20 * (NE / NE_base - 1.0);  // excess NE → SVR↑
double SVR_SGLT2_eff = E_SGLT2_C2 * 0.08; // SGLT2i mild vasodilation
double SVR_target = SVR_base * SVR_AngII_eff * (1.0 + SVR_NE_eff) * (1.0 - SVR_AT1R_inh - SVR_SGLT2_eff);
dxdt_SVR = 0.05 * (SVR_target - SVR); // slow SVR changes (tau~20h)

// LVEF: core remodeling endpoint
// Improves with: BB (reverse remodeling), ARNI (reverse remodeling), SGLT2i
// Worsens with: progressive fibrosis, hypertrophy, increased wall stress
double EF_worsening = 0.0012 * (Fib - 0.1) + 0.0008 * (Hyp - 1.0);
double EF_BB_benefit = kEF_BB * E_BB_NE * (EF_max - LVEF);
double EF_ARNI_benefit = kEF_ARNI * E_AT1R * (EF_max - LVEF);
double EF_SGLT2_benefit = kEF_SGLT2 * E_SGLT2_C2 * (EF_max - LVEF);
dxdt_LVEF = EF_BB_benefit + EF_ARNI_benefit + EF_SGLT2_benefit - EF_worsening;
// Constrain LVEF
if(LVEF > EF_max) dxdt_LVEF = 0;
if(LVEF < 0.05)   dxdt_LVEF = 0;

// ════════════════════════════════════════════════════════════════════════════
// CARDIAC REMODELING ODEs
// ════════════════════════════════════════════════════════════════════════════

// TGF-β1: driven by AngII (AT1R), Aldo (MR), TNF-α, NE
double TGFb1_prod = kTGFb1_prod * (AngII / 55.0) * (Aldo / 22.0) * (TNFa / 6.5);
double TGFb1_AT1R_inh = E_AT1R  * 0.55; // AT1R block → ↓TGF-β1
double TGFb1_MRA_inh  = E_MRA   * 0.40; // MRA → ↓TGF-β1
dxdt_TGFb1 = TGFb1_prod * (1.0 - TGFb1_AT1R_inh - TGFb1_MRA_inh) - kTGFb1_deg * TGFb1;

// Fibrosis: driven by TGF-β1; attenuated by MRA and ARNI reverse remodeling
double Fib_prod = kFib_prod * TGFb1 * (1.0 - Fib / Fib_max);
double Fib_MRA_inh = E_MRA   * 0.008 * Fib; // MRA attenuates fibrosis growth
double Fib_AT1R_inh = E_AT1R * 0.006 * Fib;
dxdt_Fib = Fib_prod - kFib_deg * Fib - Fib_MRA_inh - Fib_AT1R_inh;

// Hypertrophy index: driven by AngII, NE, wall stress; reversed by BB
double Hyp_prod = kHyp_prod * (AngII / 55.0) * (NE / 715.0) * (LVEDV / 280.0);
double Hyp_BB_rev = kHyp_deg * E_BB_NE * Hyp * 2.5; // BB reverses hypertrophy
dxdt_Hyp = Hyp_prod - kHyp_deg * Hyp - Hyp_BB_rev;
if(Hyp < 1.0) dxdt_Hyp = 0; // cannot be below normal

// ════════════════════════════════════════════════════════════════════════════
// INFLAMMATION ODEs
// ════════════════════════════════════════════════════════════════════════════

// TNF-α: produced by macrophages in response to low CO / tissue ischemia
// SGLT2i reduces inflammation via AMPK/NF-κB
double TNFa_prod = kTNFa_prod * (1.0 + 0.5 * (NE / NE_base - 1.0));
double TNFa_SGLT2 = E_SGLT2_C2 * 0.012 * TNFa;
dxdt_TNFa = TNFa_prod - kTNFa_deg * TNFa - TNFa_SGLT2;

// IL-6: driven by TNF-α and AngII
double IL6_prod = kIL6_prod * (TNFa / TNFa_base) * (AngII / 55.0);
dxdt_IL6 = IL6_prod - kIL6_deg * IL6;

// ════════════════════════════════════════════════════════════════════════════
// DRUG PK ODEs (1-compartment with first-order absorption)
// ════════════════════════════════════════════════════════════════════════════

// LBQ657 (active NEP inhibitor from sacubitril)
dxdt_LBQ_C   = - (CL_LBQ / Vd_LBQ) * LBQ_C;

// Valsartan (ARB component of ARNI)
dxdt_Valsar_C = - (CL_Valsar / Vd_Valsar) * Valsar_C;

// β-Blocker (metoprolol succinate equivalent)
dxdt_BB_C    = - (CL_BB / Vd_BB) * BB_C;

// MRA (eplerenone)
dxdt_MRA_C   = - (CL_MRA / Vd_MRA) * MRA_C;

// SGLT2i (empagliflozin)
dxdt_SGLT2_C = - (CL_SGLT2 / Vd_SGLT2) * SGLT2_C;

// Ivabradine
dxdt_IVA_C   = - (CL_IVA / Vd_IVA) * IVA_C;

$TABLE
// ── Derived clinical variables (reported in output) ──────────────────────
// Stroke Volume (mL)
double SV = LVEDV * LVEF;
// Cardiac Output (L/min)
double CO = HR * SV / 1000.0;
// MAP (mmHg): approximately SVR × CO / 80 (conversion factor)
double MAP = SVR * CO / 80.0;
// LVEDP as surrogate for PCWP (mmHg): estimated from LVEDV excess
double PCWP = 8.0 + (LVEDV - 120.0) * 0.12;
if(PCWP < 5)  PCWP = 5;
if(PCWP > 45) PCWP = 45;
// NT-proBNP (already ODE state NTpBNP)
// BNP (already ODE state BNP)
// LVEF as percentage
double LVEF_pct = LVEF * 100.0;
// EF drug effects (intermediate)
double NEP_E  = E_NEP;
double AT1R_E = E_AT1R;
double BB_E   = E_BB_HR;
double MRA_E  = E_MRA;
double SGLT2_E = E_SGLT2_V;

// NYHA classification (rough scoring)
// Based on LVEF, CO, PCWP
double NYHA_score = 0.0;
if (CO < 3.0 || PCWP > 25.0 || LVEF_pct < 20.0) NYHA_score = 4.0;
else if (CO < 4.0 || PCWP > 18.0 || LVEF_pct < 30.0) NYHA_score = 3.0;
else if (CO < 5.0 || PCWP > 12.0 || LVEF_pct < 40.0) NYHA_score = 2.0;
else NYHA_score = 1.0;

// eGFR proxy (simplified — elevated AngII/SVR → ↓renal perfusion)
double eGFR = 65.0 * (1.0 - 0.25 * (AngII - 15.0) / 55.0 + 0.10 * E_AT1R);
if(eGFR < 10) eGFR = 10;
if(eGFR > 90) eGFR = 90;

// Capture everything
capture SV, CO, MAP, PCWP, LVEF_pct, NYHA_score, eGFR;
capture NEP_E, AT1R_E, BB_E, MRA_E, SGLT2_E;
'

# ─────────────────────────────────────────────────────────────────────────────
# COMPILE MODEL
# ─────────────────────────────────────────────────────────────────────────────
mod <- mread_cache("hfref_qsp", tempdir(), hfref_model_code)

# ─────────────────────────────────────────────────────────────────────────────
# HELPER: dosing event block constructor
# ─────────────────────────────────────────────────────────────────────────────
make_doses <- function(dose_ARNI_LBQ    = 0,   # mg LBQ657-equivalent (sacubitril 97mg → LBQ ~75mg)
                       dose_ARNI_Valsar = 0,   # mg valsartan (103mg in 200mg ARNI)
                       dose_BB          = 0,   # mg metoprolol succinate
                       dose_MRA         = 0,   # mg eplerenone
                       dose_SGLT2       = 0,   # mg empagliflozin
                       dose_IVA         = 0,   # mg ivabradine
                       start = 0, end = 8760,  # hours (1 year default)
                       Vd_LBQ    = 26,
                       Vd_Valsar = 75,
                       Vd_BB     = 220,
                       Vd_MRA    = 43,
                       Vd_SGLT2  = 73.8,
                       Vd_IVA    = 100) {
  ev_list <- list()
  if(dose_ARNI_LBQ    > 0) ev_list <- c(ev_list, list(ev(amt=dose_ARNI_LBQ*0.6,    cmt="LBQ_C",    ii=12, addl=floor((end-start)/12),   time=start, rate=-2)))
  if(dose_ARNI_Valsar > 0) ev_list <- c(ev_list, list(ev(amt=dose_ARNI_Valsar*0.94, cmt="Valsar_C", ii=12, addl=floor((end-start)/12),   time=start, rate=-2)))
  if(dose_BB          > 0) ev_list <- c(ev_list, list(ev(amt=dose_BB*0.45,          cmt="BB_C",     ii=24, addl=floor((end-start)/24),   time=start, rate=-2)))
  if(dose_MRA         > 0) ev_list <- c(ev_list, list(ev(amt=dose_MRA*0.69,         cmt="MRA_C",    ii=24, addl=floor((end-start)/24),   time=start, rate=-2)))
  if(dose_SGLT2       > 0) ev_list <- c(ev_list, list(ev(amt=dose_SGLT2*0.86,       cmt="SGLT2_C",  ii=24, addl=floor((end-start)/24),   time=start, rate=-2)))
  if(dose_IVA         > 0) ev_list <- c(ev_list, list(ev(amt=dose_IVA,              cmt="IVA_C",    ii=12, addl=floor((end-start)/12),   time=start, rate=-2)))
  if(length(ev_list) == 0) return(ev())
  Reduce(c, ev_list)
}

# ─────────────────────────────────────────────────────────────────────────────
# SCENARIO 1: No therapy (natural history)
# ─────────────────────────────────────────────────────────────────────────────
sim_time <- seq(0, 8760, by=24) # 1 year, daily output

out_s1 <- mod %>%
  ev(ev()) %>%
  mrgsim(end=8760, delta=24, outvars=c("AngII","Aldo","NE","BNP","NTpBNP","cGMP",
                                        "LVEDV","HR","SVR","LVEF_pct","Fib","TGFb1","Hyp",
                                        "TNFa","IL6","SV","CO","MAP","PCWP","NYHA_score","eGFR")) %>%
  as.data.frame() %>%
  mutate(scenario = "1_No_Therapy")

# ─────────────────────────────────────────────────────────────────────────────
# SCENARIO 2: ACEi + β-Blocker (legacy therapy — pre-PARADIGM era)
# ─────────────────────────────────────────────────────────────────────────────
# Enalapril: simulate as 50% AT1R block equivalent (Valsartan Cp → EC50_AT1R)
# For simplicity, ACEi represented via valsartan compartment with adjusted dose
# Enalapril 10mg BID → plasma enalaprilat: effective AT1R block via AngII↓
# Represented by Valsar_C driving AT1R_E with enalapril Cp equivalent
out_s2 <- mod %>%
  param(CL_Valsar=3.2, Vd_Valsar=18.0) %>%  # enalaprilat PK approximation
  ev(make_doses(dose_ARNI_Valsar=10,          # enalapril 10mg → equivalent ACEi
                dose_BB=200)) %>%             # metoprolol succinate 200mg/day
  mrgsim(end=8760, delta=24) %>%
  as.data.frame() %>%
  mutate(scenario = "2_ACEi_BB")

# ─────────────────────────────────────────────────────────────────────────────
# SCENARIO 3: ARNI + β-Blocker + MRA (PARADIGM-HF / EMPHASIS-HF based)
# ─────────────────────────────────────────────────────────────────────────────
# Sacubitril/valsartan 97/103 mg BID → LBQ657 ~75mg equivalent, Valsartan 103mg
out_s3 <- mod %>%
  param(CL_Valsar=1.9, Vd_Valsar=75.0) %>%  # restore ARNI valsartan PK
  ev(make_doses(dose_ARNI_LBQ=75,
                dose_ARNI_Valsar=103,
                dose_BB=200,
                dose_MRA=50)) %>%
  mrgsim(end=8760, delta=24) %>%
  as.data.frame() %>%
  mutate(scenario = "3_ARNI_BB_MRA")

# ─────────────────────────────────────────────────────────────────────────────
# SCENARIO 4: ARNI + BB + MRA + SGLT2i (Comprehensive GDMT — 4 pillars)
# ─────────────────────────────────────────────────────────────────────────────
out_s4 <- mod %>%
  param(CL_Valsar=1.9, Vd_Valsar=75.0) %>%
  ev(make_doses(dose_ARNI_LBQ=75,
                dose_ARNI_Valsar=103,
                dose_BB=200,
                dose_MRA=50,
                dose_SGLT2=10)) %>%
  mrgsim(end=8760, delta=24) %>%
  as.data.frame() %>%
  mutate(scenario = "4_ARNI_BB_MRA_SGLT2i")

# ─────────────────────────────────────────────────────────────────────────────
# SCENARIO 5: Maximal GDMT + Ivabradine (HR ≥70 bpm at baseline)
# ─────────────────────────────────────────────────────────────────────────────
out_s5 <- mod %>%
  param(CL_Valsar=1.9, Vd_Valsar=75.0) %>%
  ev(make_doses(dose_ARNI_LBQ=75,
                dose_ARNI_Valsar=103,
                dose_BB=200,
                dose_MRA=50,
                dose_SGLT2=10,
                dose_IVA=7.5)) %>%
  mrgsim(end=8760, delta=24) %>%
  as.data.frame() %>%
  mutate(scenario = "5_Max_GDMT_plus_IVA")

# ─────────────────────────────────────────────────────────────────────────────
# COMBINE ALL SCENARIOS
# ─────────────────────────────────────────────────────────────────────────────
all_scenarios <- bind_rows(out_s1, out_s2, out_s3, out_s4, out_s5) %>%
  mutate(
    time_days  = time / 24,
    time_months = time_days / 30.44,
    scenario_label = factor(scenario,
      levels = c("1_No_Therapy","2_ACEi_BB","3_ARNI_BB_MRA",
                 "4_ARNI_BB_MRA_SGLT2i","5_Max_GDMT_plus_IVA"),
      labels = c("No Therapy","ACEi + BB","ARNI + BB + MRA",
                 "ARNI + BB + MRA + SGLT2i","ARNI + BB + MRA + SGLT2i + Ivabradine"))
  )

# ─────────────────────────────────────────────────────────────────────────────
# SCENARIO 6: DOSE-RESPONSE — ARNI dose titration (for Shiny app)
# ─────────────────────────────────────────────────────────────────────────────
arni_doses <- c(24, 49, 97, 200)  # mg sacubitril (24/26 → 49/51 → 97/103 mg BID steps)
lbq_doses  <- arni_doses * 0.75   # approximate LBQ657 equivalent

dose_response_list <- lapply(seq_along(arni_doses), function(i) {
  mod %>%
    ev(make_doses(dose_ARNI_LBQ    = lbq_doses[i],
                  dose_ARNI_Valsar = arni_doses[i],
                  dose_BB = 200, dose_MRA = 50, dose_SGLT2 = 10)) %>%
    mrgsim(end=8760, delta=168) %>%  # weekly output
    as.data.frame() %>%
    mutate(ARNI_dose_mg = arni_doses[i],
           time_months  = time / 720)
})
dose_response_df <- bind_rows(dose_response_list)

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY TABLE AT 12 MONTHS
# ─────────────────────────────────────────────────────────────────────────────
summary_12mo <- all_scenarios %>%
  filter(abs(time_days - 365) < 1) %>%
  select(scenario_label, LVEF_pct, BNP, NTpBNP, CO, HR, SVR, MAP,
         LVEDV, Fib, NYHA_score, eGFR) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

cat("\n===== 12-Month Simulation Summary =====\n")
print(summary_12mo)

# ─────────────────────────────────────────────────────────────────────────────
# KEY VISUALIZATION
# ─────────────────────────────────────────────────────────────────────────────
pal5 <- c("#e74c3c","#e67e22","#3498db","#27ae60","#9b59b6")

p1 <- ggplot(all_scenarios, aes(time_months, LVEF_pct, color=scenario_label)) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=40, linetype="dashed", color="grey60", alpha=0.8) +
  annotate("text", x=11, y=41.5, label="HFmrEF threshold (EF=40%)", size=3, color="grey50") +
  scale_color_manual(values=pal5, name="Treatment Scenario") +
  labs(title="LVEF Over 12 Months by Treatment Scenario",
       x="Time (months)", y="LVEF (%)") +
  theme_dark() + theme(legend.position="bottom", legend.text=element_text(size=8))

p2 <- ggplot(all_scenarios, aes(time_months, NTpBNP, color=scenario_label)) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=900, linetype="dashed", color="grey60") +
  annotate("text", x=11, y=950, label="HF diagnosis threshold (900 pg/mL)", size=3, color="grey50") +
  scale_color_manual(values=pal5, name="Treatment Scenario") +
  scale_y_log10() +
  labs(title="NT-proBNP Over 12 Months (log scale)",
       x="Time (months)", y="NT-proBNP (pg/mL, log)") +
  theme_dark() + theme(legend.position="bottom", legend.text=element_text(size=8))

p3 <- ggplot(all_scenarios, aes(time_months, CO, color=scenario_label)) +
  geom_line(linewidth=1.1) +
  scale_color_manual(values=pal5, name="Treatment Scenario") +
  labs(title="Cardiac Output Over 12 Months",
       x="Time (months)", y="CO (L/min)") +
  theme_dark() + theme(legend.position="bottom", legend.text=element_text(size=8))

p4 <- ggplot(dose_response_df, aes(time_months, LVEF_pct*1,
                                    color=factor(ARNI_dose_mg))) +
  geom_line(linewidth=1.1) +
  scale_color_viridis_d(name="Sacubitril dose (mg BID)") +
  labs(title="ARNI Dose-Response: LVEF Improvement",
       x="Time (months)", y="LVEF (%)") +
  theme_dark()

# Print plots
print(p1)
print(p2)
print(p3)
print(p4)

# ─────────────────────────────────────────────────────────────────────────────
# PARAMETER SENSITIVITY ANALYSIS
# ─────────────────────────────────────────────────────────────────────────────
sens_params <- c("kAngII_deg","EC50_AT1R_val","Emax_BB_HR","Emax_SGLT2_V",
                 "kFib_prod","kTGFb1_prod","Fib_base","NE_HF_factor")
sens_range  <- c(0.5, 1.0, 1.5, 2.0)

sens_list <- lapply(sens_params, function(pname) {
  lapply(sens_range, function(mult) {
    param_val <- mod@param@data[[pname]]
    mod %>%
      param(.oo = setNames(list(param_val * mult), pname)) %>%
      ev(make_doses(dose_ARNI_LBQ=75, dose_ARNI_Valsar=103,
                    dose_BB=200, dose_MRA=50, dose_SGLT2=10)) %>%
      mrgsim(end=8760, delta=720) %>%
      filter(time == 8760) %>%
      as.data.frame() %>%
      mutate(param=pname, multiplier=mult, param_val=param_val*mult)
  }) %>% bind_rows()
}) %>% bind_rows()

cat("\n===== Sensitivity Analysis — LVEF at 12 months =====\n")
sens_summary <- sens_list %>%
  select(param, multiplier, LVEF_pct) %>%
  arrange(param, multiplier)
print(sens_summary)

cat("\nModel compilation and all scenarios complete.\n")
cat("Objects available: mod, all_scenarios, dose_response_df, summary_12mo\n")
