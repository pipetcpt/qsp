##############################################################################
# Liver Cirrhosis QSP Model — mrgsolve ODE-based PK/PD
# ============================================================
# Disease:   Liver Cirrhosis (Hepatic Fibrosis, Portal Hypertension,
#            Ascites, Hepatic Encephalopathy, Variceal Bleeding)
# Model:     Multi-compartment ODE (21 compartments)
#            PK: Propranolol (2-cpt), Spironolactone (1-cpt),
#                Terlipressin (effect-cpt), Rifaximin (gut-lumen)
#            PD: Fibrosis, portal pressure, hepatic function,
#                ascites, ammonia/HE, renal function
# Scenarios: 5 therapeutic strategies
# Calibration notes:
#   - Fibrosis progression rate from Fattovich 2004 (HEPATOLOGY)
#   - HVPG–variceal bleeding from Groszmann 2005 (NEJM)
#   - Propranolol PK from Regardh 1983; cirrhosis adjustment from
#     Johnson 1995 (clearance reduced 50–70%)
#   - Spironolactone PK from Sungaila 1992
#   - MELD score calibration from Kim 2008 (HEPATOLOGY)
#   - Ascites natural history from Planas 2006 (HEPATOLOGY)
##############################################################################

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

##############################################################################
# MODEL CODE
##############################################################################

code <- '
$PROB Liver Cirrhosis QSP Model v1.0

$PARAM @annotated
// ----- Disease Stage (initial) -----
F0_INIT   : 0.2  : Baseline fibrosis index (0=none, 1=cirrhosis F4)
STAGE_CP  : 2.0  : Child-Pugh stage initial (1=A, 2=B, 3=C)

// ----- Fibrosis Dynamics -----
K_FIBRO   : 0.004 : Fibrosis progression rate (per day) — Fattovich 2004
K_FIBRO_RES: 0.001: Spontaneous fibrosis resolution rate (per day)
K_HSC_ACT : 0.015 : HSC activation rate from TGF-β signal
K_HSC_DEACT: 0.008: HSC deactivation (apoptosis/senescence) rate
TGF_BASE  : 1.0  : Basal TGF-β1 signal (relative units)
HILL_TGF  : 2.0  : Hill exponent for TGF-β fibrogenesis
EC50_TGF  : 0.5  : EC50 of TGF-β for HSC activation

// ----- Portal Hypertension -----
HVPG_BASE : 6.0  : Baseline HVPG (mmHg, normal <5)
K_HVPG_F  : 12.0 : Scaling: fibrosis → HVPG increase (mmHg at F=1)
ET1_BASE  : 1.0  : Endothelin-1 relative level
eNOS_BASE : 1.0  : eNOS activity (1=normal)
K_ET1     : 0.3  : ET1 contribution to HVPG
K_eNOS    : 0.2  : eNOS deficit contribution to HVPG

// ----- Hepatic Synthetic Function -----
ALB_NORM  : 4.2  : Normal albumin (g/dL)
K_ALB_DEC : 0.35 : Rate of albumin decrease per unit fibrosis index
INR_NORM  : 1.0  : Normal INR
K_INR_INC : 0.8  : Rate of INR increase per unit fibrosis index
BILI_NORM : 0.8  : Normal bilirubin (mg/dL)
K_BILI_INC: 3.5  : Rate of bilirubin increase per unit fibrosis index

// ----- Ascites -----
K_ASCITES : 0.05 : Ascites formation rate (L/day per mmHg HVPG excess)
HVPG_THRESH: 10.0: HVPG threshold for ascites formation (mmHg)
K_LYMPH   : 0.02 : Lymphatic absorption rate
ALDO_BASE : 1.0  : Baseline aldosterone relative activity
K_ALDO_NA : 0.04 : Na retention effect of aldosterone on ascites
OncP_NORM : 1.0  : Normal oncotic pressure (relative)
K_ONCP_ALB: 0.3  : Relationship between albumin and oncotic pressure

// ----- Renal Function -----
GFR_NORM  : 90.0 : Normal GFR (mL/min/1.73m2)
K_GFR_HVPG: 0.8  : GFR reduction per mmHg HVPG above 12
K_SNS_GFR : 0.5  : SNS activation effect on GFR
CREAT_NORM: 0.9  : Normal creatinine (mg/dL)
K_CREAT_GFR: 80.0: Creatinine–GFR inverse relationship constant

// ----- Ammonia & HE -----
NH3_NORM  : 30.0 : Normal blood ammonia (μmol/L)
K_NH3_LF  : 60.0 : Ammonia increase from liver failure (μmol/L at F=1)
K_NH3_BYPASS: 20.0: Ammonia increase from portosystemic shunting
K_NH3_GUT : 1.0  : Gut ammonia production rate (relative)
K_NH3_ELIM: 0.1  : Ammonia elimination rate (per day, residual)
HE_THRESH : 50.0 : Ammonia threshold for HE onset (μmol/L)
HE_K      : 0.02 : HE grade progression rate from ammonia

// ----- Propranolol PK (2-compartment oral) -----
PROP_F    : 0.25 : Bioavailability (25% due to first-pass; reduced in cirrhosis)
PROP_KA   : 0.5  : Absorption rate constant (h^-1)
PROP_CL   : 30.0 : Clearance (L/h); cirrhosis ~50% reduction vs normal 60 L/h
PROP_V1   : 150.0: Central volume (L)
PROP_Q    : 20.0 : Inter-compartment clearance (L/h)
PROP_V2   : 300.0: Peripheral volume (L)
PROP_EC50 : 40.0 : EC50 for HR reduction (ng/mL)
PROP_EMAX : 0.30 : Maximum fractional HR reduction
HR_BASE   : 85.0 : Baseline heart rate (bpm) — elevated in cirrhosis
PROP_HVPG_EMAX: 0.25: Max fractional HVPG reduction by propranolol
PROP_HVPG_EC50: 35.0: EC50 for HVPG reduction (ng/mL)

// ----- Spironolactone PK (1-compartment) -----
SPIRO_F   : 0.90 : Bioavailability
SPIRO_KA  : 0.3  : Absorption rate constant (h^-1)
SPIRO_CL  : 3.5  : Clearance (L/h)
SPIRO_V   : 70.0 : Volume of distribution (L)
SPIRO_EC50: 150.0: EC50 for aldosterone blockade (ng/mL)
SPIRO_EMAX: 0.85 : Max aldosterone blockade fraction
SPIRO_NA_EMAX: 0.7 : Max fractional Na excretion increase

// ----- Terlipressin (IV, effect-compartment) -----
TERL_CL   : 15.0 : Clearance (L/h)
TERL_V    : 30.0 : Volume (L)
TERL_KE0  : 0.4  : Effect-compartment equilibration (h^-1)
TERL_EC50 : 10.0 : EC50 for splanchnic vasoconstriction (ng/mL)
TERL_EMAX : 0.50 : Max splanchnic constriction (fraction)
TERL_HRS_EMAX: 0.65: Max creatinine reversal fraction

// ----- Rifaximin PK (gut-lumen model) -----
RIFAX_F   : 0.001: Systemic bioavailability (<0.4% — gut-restricted)
RIFAX_KA_GUT: 0.1: Rate of gut distribution
RIFAX_KOUT: 0.15 : Gut elimination rate (h^-1)
RIFAX_EC50_NH3: 200.0: EC50 for gut NH3 reduction (ng/mL)
RIFAX_EMAX_NH3: 0.55 : Max gut NH3 reduction fraction

// ----- Antifibrotic effect parameters -----
AF_FIBRO_EMAX: 0.6 : Max fibrosis reduction by antifibrotics
AF_HSC_EMAX : 0.7  : Max HSC activation reduction
AF_EC50     : 1.0  : Antifibrotic drug relative concentration EC50

$CMT @annotated
// Propranolol compartments
PROP_GUT  : Propranolol gut (absorption depot) [ng]
PROP_C1   : Propranolol central compartment [ng]
PROP_C2   : Propranolol peripheral compartment [ng]

// Spironolactone compartments
SPIRO_GUT : Spironolactone gut [ng]
SPIRO_C1  : Spironolactone central compartment [ng]

// Terlipressin
TERL_C    : Terlipressin central [ng]
TERL_CE   : Terlipressin effect compartment [ng]

// Rifaximin
RIFAX_GUT : Rifaximin gut lumen [ng]

// Disease compartments
FIBRO     : Fibrosis index [0–1 scale]
HSC_ACT   : Activated HSC fraction [0–1]
HVPG      : Hepatic venous pressure gradient [mmHg]
ALB       : Serum albumin [g/dL]
BILIRUBIN : Total bilirubin [mg/dL]
INR_val   : INR value
ASCITES   : Ascites volume [L]
GFR_est   : Estimated GFR [mL/min/1.73m2]
NH3_blood : Blood ammonia [μmol/L]
HE_GRADE  : Hepatic encephalopathy grade [0–4 scale]
CREAT     : Serum creatinine [mg/dL]
ALDO_ACT  : Aldosterone activity [relative, 1=normal]

$MAIN
// ----- Derived PK concentrations -----
double PROP_CONC = PROP_C1 / PROP_V1;     // ng/mL
double SPIRO_CONC = SPIRO_C1 / SPIRO_V;   // ng/mL
double TERL_CONC_EFF = TERL_CE;           // ng/mL (effect compartment normalized)
double RIFAX_CONC = RIFAX_GUT / 50.0;     // ng/mL (gut lumen pseudo-concentration)

// ----- Pharmacodynamic effects -----
// Propranolol: HR reduction (Emax model)
double PROP_HR_EFF = PROP_EMAX * PROP_CONC / (PROP_EC50 + PROP_CONC);
double HR_obs = HR_BASE * (1.0 - PROP_HR_EFF);

// Propranolol: HVPG reduction via reduced cardiac output
double PROP_HVPG_EFF = PROP_HVPG_EMAX * PROP_CONC / (PROP_HVPG_EC50 + PROP_CONC);

// Spironolactone: Aldosterone blockade
double SPIRO_ALDO_EFF = SPIRO_EMAX * SPIRO_CONC / (SPIRO_EC50 + SPIRO_CONC);
double SPIRO_NA_EFF = SPIRO_NA_EMAX * SPIRO_CONC / (SPIRO_EC50 + SPIRO_CONC);

// Terlipressin: Splanchnic vasoconstriction → portal pressure reduction
double TERL_VASOC_EFF = TERL_EMAX * TERL_CONC_EFF / (TERL_EC50 + TERL_CONC_EFF);
double TERL_HRS_EFF = TERL_HRS_EMAX * TERL_CONC_EFF / (TERL_EC50 + TERL_CONC_EFF);

// Rifaximin: Gut NH3 reduction
double RIFAX_NH3_EFF = RIFAX_EMAX_NH3 * RIFAX_CONC / (RIFAX_EC50_NH3 + RIFAX_CONC);

// ----- Composite HVPG calculation -----
// Driven by fibrosis, ET1/eNOS balance, modulated by drugs
double HVPG_fibrosis = HVPG_BASE + K_HVPG_F * FIBRO;
double HVPG_ET1_eNOS = K_ET1 * ET1_BASE - K_eNOS * eNOS_BASE;
double HVPG_drug_mod = (1.0 - PROP_HVPG_EFF) * (1.0 - TERL_VASOC_EFF);
double HVPG_calc = (HVPG_fibrosis + HVPG_ET1_eNOS) * HVPG_drug_mod;
HVPG_calc = (HVPG_calc < 1.0) ? 1.0 : HVPG_calc;

// ----- Albumin, Bilirubin, INR (algebraic from fibrosis) -----
double ALB_calc = ALB_NORM - K_ALB_DEC * FIBRO;
ALB_calc = (ALB_calc < 1.0) ? 1.0 : ALB_calc;

double BILI_calc = BILI_NORM + K_BILI_INC * FIBRO;
double INR_calc = INR_NORM + K_INR_INC * FIBRO;

// Child-Pugh score components (simplified continuous scoring)
// ALB: >3.5=1, 2.8–3.5=2, <2.8=3
double CP_ALB = (ALB_calc > 3.5) ? 1.0 : ((ALB_calc > 2.8) ? 2.0 : 3.0);
// BILI: <2=1, 2–3=2, >3=3
double CP_BILI = (BILI_calc < 2.0) ? 1.0 : ((BILI_calc < 3.0) ? 2.0 : 3.0);
// INR: <1.7=1, 1.7–2.3=2, >2.3=3
double CP_INR = (INR_calc < 1.7) ? 1.0 : ((INR_calc < 2.3) ? 2.0 : 3.0);
// Ascites: 0=1, mild=2, moderate=3
double CP_ASC = (ASCITES < 1.0) ? 1.0 : ((ASCITES < 5.0) ? 2.0 : 3.0);
// HE: 0=1, grade1-2=2, grade3-4=3
double CP_HE = (HE_GRADE < 0.5) ? 1.0 : ((HE_GRADE < 2.5) ? 2.0 : 3.0);

double CHILD_PUGH = CP_ALB + CP_BILI + CP_INR + CP_ASC + CP_HE;

// MELD score = 3.78×ln(Bili) + 11.2×ln(INR) + 9.57×ln(Creat) + 6.43
double CREAT_MELD = CREAT;
CREAT_MELD = (CREAT_MELD < 1.0) ? 1.0 : ((CREAT_MELD > 4.0) ? 4.0 : CREAT_MELD);
double BILI_MELD = BILIRUBIN;
BILI_MELD = (BILI_MELD < 1.0) ? 1.0 : BILI_MELD;
double INR_MELD = INR_val;
INR_MELD = (INR_MELD < 1.0) ? 1.0 : INR_MELD;
double MELD_score = 3.78 * log(BILI_MELD) + 11.2 * log(INR_MELD) + 9.57 * log(CREAT_MELD) + 6.43;

// MELD-Na = MELD + 1.32*(140-Na) - [0.24*MELD*(140-Na)]
// Assuming Na approximated from GFR/ascites state
double Na_approx = 140.0 - 5.0 * (ASCITES / 5.0) - 3.0 * (1.0 - GFR_est/90.0);
Na_approx = (Na_approx < 120.0) ? 120.0 : ((Na_approx > 140.0) ? 140.0 : Na_approx);
double MELD_Na_score = MELD_score + 1.32*(140.0-Na_approx) - 0.24*MELD_score*(140.0-Na_approx)/100.0;

// Variceal bleeding risk (annual, from HVPG)
double VAR_BLEED_RISK = (HVPG_calc > 12.0) ?
    0.02 + 0.025 * (HVPG_calc - 12.0) : 0.0;

// 1-year mortality risk (MELD-based, O'Leary 2011)
double MORT_1YR = (MELD_score < 9) ? 0.02 :
                  (MELD_score < 19) ? 0.06 :
                  (MELD_score < 29) ? 0.20 :
                  (MELD_score < 39) ? 0.52 : 0.70;

$ODE
// ===== PROPRANOLOL PK =====
double prop_abs_rate = PROP_KA * PROP_GUT;
double prop_dist_forward = (PROP_Q/PROP_V1) * PROP_C1;
double prop_dist_back    = (PROP_Q/PROP_V2) * PROP_C2;
double prop_elim = (PROP_CL/PROP_V1) * PROP_C1;

dxdt_PROP_GUT = -prop_abs_rate;
dxdt_PROP_C1  =  prop_abs_rate - prop_elim - prop_dist_forward + prop_dist_back;
dxdt_PROP_C2  =  prop_dist_forward - prop_dist_back;

// ===== SPIRONOLACTONE PK =====
double spiro_abs_rate = SPIRO_KA * SPIRO_GUT;
double spiro_elim = (SPIRO_CL / SPIRO_V) * SPIRO_C1;

dxdt_SPIRO_GUT = -spiro_abs_rate;
dxdt_SPIRO_C1  =  spiro_abs_rate - spiro_elim;

// ===== TERLIPRESSIN PK =====
double terl_elim = (TERL_CL / TERL_V) * TERL_C;
dxdt_TERL_C  = -terl_elim;
dxdt_TERL_CE =  TERL_KE0 * (TERL_C/TERL_V - TERL_CE);

// ===== RIFAXIMIN (gut lumen) =====
dxdt_RIFAX_GUT = -RIFAX_KOUT * RIFAX_GUT;

// ===== HSC ACTIVATION DYNAMICS =====
double TGF_effective = TGF_BASE * (1.0 + 2.0 * FIBRO);  // TGF-β scales with fibrosis
double HSC_activation_in = K_HSC_ACT * TGF_effective * TGF_effective /
    (EC50_TGF * EC50_TGF + TGF_effective * TGF_effective) * (1.0 - HSC_ACT);

// Antifibrotic drug effect on HSC (captured via AF_DRUG parameter = 0 or 1)
double HSC_activation_net = HSC_activation_in * (1.0 - AF_HSC_EMAX * AF_DRUG_CONC / (AF_EC50 + AF_DRUG_CONC));
double HSC_deactivation = K_HSC_DEACT * HSC_ACT;
dxdt_HSC_ACT = HSC_activation_net - HSC_deactivation;

// ===== FIBROSIS INDEX DYNAMICS =====
double fibro_production = K_FIBRO * HSC_ACT * (1.0 - FIBRO);
double fibro_resolution = K_FIBRO_RES * MMP_ACT_NORM / (TIMP_NORM + MMP_ACT_NORM) * FIBRO;
// Antifibrotic drug reduces net fibrogenesis
double fibro_af_eff = AF_FIBRO_EMAX * AF_DRUG_CONC / (AF_EC50 + AF_DRUG_CONC);
dxdt_FIBRO = fibro_production * (1.0 - fibro_af_eff) - fibro_resolution;
if(FIBRO > 0.999 && dxdt_FIBRO > 0) dxdt_FIBRO = 0;
if(FIBRO < 0.001 && dxdt_FIBRO < 0) dxdt_FIBRO = 0;

// ===== HVPG DYNAMICS (tracking calculated value with smoothing) =====
dxdt_HVPG = 0.05 * (HVPG_calc - HVPG);  // lag toward calculated HVPG

// ===== HEPATIC SYNTHETIC MARKERS (track calculated values) =====
dxdt_ALB       = 0.02 * (ALB_calc - ALB);
dxdt_BILIRUBIN = 0.02 * (BILI_calc - BILIRUBIN);
dxdt_INR_val   = 0.02 * (INR_calc - INR_val);

// ===== ASCITES DYNAMICS =====
// Formation driven by portal hypertension and low oncotic pressure
double portal_drive = (HVPG > HVPG_THRESH) ? (HVPG - HVPG_THRESH) : 0.0;
double onco_factor = 1.0 + K_ONCP_ALB * (ALB_NORM - ALB) / ALB_NORM;
double aldo_current = ALDO_ACT * (1.0 - SPIRO_ALDO_EFF);
double ascites_form = K_ASCITES * portal_drive * onco_factor * (1.0 + K_ALDO_NA * aldo_current);
double ascites_absorb = K_LYMPH * ASCITES;
double ascites_diuretic = SPIRO_NA_EFF * 0.1 * ASCITES; // diuretic mobilization
dxdt_ASCITES = ascites_form - ascites_absorb - ascites_diuretic;
if(ASCITES < 0.0 && dxdt_ASCITES < 0) dxdt_ASCITES = 0;

// ===== ALDOSTERONE ACTIVITY =====
// Elevated due to hyperdynamic circulation and RAAS activation
double aldo_drive = 1.0 + 0.5 * FIBRO + 0.03 * (HVPG - 6.0);
aldo_drive = (aldo_drive < 0.0) ? 0.0 : aldo_drive;
dxdt_ALDO_ACT = 0.01 * (aldo_drive * (1.0 - SPIRO_ALDO_EFF) - ALDO_ACT);

// ===== GFR / RENAL FUNCTION =====
double gfr_hvpg_drive = (HVPG > 12.0) ? K_GFR_HVPG * (HVPG - 12.0) : 0.0;
double gfr_sns_drive = K_SNS_GFR * (ALDO_ACT - 1.0);
gfr_sns_drive = (gfr_sns_drive < 0.0) ? 0.0 : gfr_sns_drive;
double gfr_target = GFR_NORM - gfr_hvpg_drive - gfr_sns_drive;
// Terlipressin reverses renal vasoconstriction in HRS
double gfr_terl_benefit = TERL_HRS_EFF * (GFR_NORM - gfr_target);
gfr_target += gfr_terl_benefit;
gfr_target = (gfr_target < 5.0) ? 5.0 : gfr_target;
dxdt_GFR_est = 0.02 * (gfr_target - GFR_est);

// ===== CREATININE =====
// Creatinine inversely related to GFR
double creat_target = CREAT_NORM * GFR_NORM / GFR_est;
creat_target = (creat_target < 0.5) ? 0.5 : ((creat_target > 10.0) ? 10.0 : creat_target);
dxdt_CREAT = 0.01 * (creat_target - CREAT);

// ===== BLOOD AMMONIA =====
// Ammonia rises with reduced liver detox, portosystemic shunting
double nh3_liver_fail = K_NH3_LF * FIBRO;
double nh3_shunt = K_NH3_BYPASS * (HVPG > 10.0 ? (HVPG - 10.0) / 10.0 : 0.0);
double nh3_gut = K_NH3_GUT * NH3_NORM * (1.0 - RIFAX_NH3_EFF);
double nh3_target = NH3_NORM + nh3_liver_fail + nh3_shunt + nh3_gut;
double nh3_elim = K_NH3_ELIM * NH3_blood;
dxdt_NH3_blood = 0.02 * (nh3_target - NH3_blood) - nh3_elim * 0.01;
if(NH3_blood < 10.0 && dxdt_NH3_blood < 0) dxdt_NH3_blood = 0;

// ===== HE GRADE =====
// HE grade driven by ammonia above threshold
double he_drive = (NH3_blood > HE_THRESH) ?
    HE_K * (NH3_blood - HE_THRESH) * (4.0 - HE_GRADE) : 0.0;
double he_resolution = 0.05 * HE_GRADE * (1.0 - NH3_blood / (NH3_blood + HE_THRESH));
dxdt_HE_GRADE = he_drive - he_resolution;
if(HE_GRADE < 0.0 && dxdt_HE_GRADE < 0) dxdt_HE_GRADE = 0;
if(HE_GRADE > 4.0 && dxdt_HE_GRADE > 0) dxdt_HE_GRADE = 0;

$PARAM @annotated
// Additional parameters for ODE (MMP/TIMP ratio, antifibrotic drug)
MMP_ACT_NORM : 1.0  : MMP activity (normalized, 1=normal)
TIMP_NORM    : 2.0  : TIMP-1 level (relative; elevated in cirrhosis)
AF_DRUG_CONC : 0.0  : Antifibrotic drug relative concentration (0=none, 1=effective)

$CAPTURE @annotated
PROP_CONC    : Propranolol plasma concentration [ng/mL]
SPIRO_CONC   : Spironolactone plasma concentration [ng/mL]
TERL_CONC_EFF: Terlipressin effect-cpt concentration [ng/mL]
RIFAX_CONC   : Rifaximin gut concentration [ng/mL]
HR_obs       : Observed heart rate [bpm]
HVPG_calc    : Calculated HVPG [mmHg]
ALB_calc     : Serum albumin [g/dL]
BILI_calc    : Bilirubin [mg/dL]
INR_calc     : INR
CHILD_PUGH   : Child-Pugh score (continuous)
MELD_score   : MELD score
MELD_Na_score: MELD-Na score
VAR_BLEED_RISK: Annual variceal bleeding risk
MORT_1YR     : 1-year mortality estimate
Na_approx    : Estimated serum sodium [mEq/L]

$INIT
FIBRO     = 0.55  // Starting at F2-F3 bridging fibrosis
HSC_ACT   = 0.30  // 30% HSC activated
HVPG      = 10.0  // mmHg — clinically significant portal hypertension
ALB       = 3.1   // g/dL — low-normal
BILIRUBIN = 2.5   // mg/dL
INR_val   = 1.6   // Mildly elevated
ASCITES   = 2.0   // 2L mild ascites
GFR_est   = 65.0  // Mildly reduced
NH3_blood = 65.0  // Mildly elevated
HE_GRADE  = 0.5   // Covert / grade 1 HE
CREAT     = 1.2   // mg/dL
ALDO_ACT  = 1.8   // Elevated RAAS
PROP_GUT  = 0.0
PROP_C1   = 0.0
PROP_C2   = 0.0
SPIRO_GUT = 0.0
SPIRO_C1  = 0.0
TERL_C    = 0.0
TERL_CE   = 0.0
RIFAX_GUT = 0.0
'

##############################################################################
# COMPILE MODEL
##############################################################################

mod <- mcode("liver_cirrhosis_qsp", code, quiet = TRUE)

##############################################################################
# SIMULATION SCENARIOS
##############################################################################

# Time grid: daily for 2 years (730 days), hourly for acute phase (72h)
t_chronic <- seq(0, 730, by = 1)
t_acute   <- seq(0, 72, by = 0.5)

#---------------------------------------------------------------------
# SCENARIO 1: Natural History (No Treatment)
#---------------------------------------------------------------------
e_natural <- ev(time = 0, cmt = 1, amt = 0)  # no drug

out_natural <- mod %>%
    param(AF_DRUG_CONC = 0) %>%
    mrgsim(events = e_natural, end = 730, delta = 1) %>%
    as_tibble() %>%
    mutate(Scenario = "Natural History")

#---------------------------------------------------------------------
# SCENARIO 2: Propranolol 40 mg BID (standard NSBB for variceal prophylaxis)
# Dosing: 40 mg BID = 80 mg/day
# Each dose: 40 mg = 40,000 ng
#---------------------------------------------------------------------
prop_dose_amt <- 40e3 * 0.25  # ng, accounting for bioavailability in dose
e_propranolol <- ev(
    time = seq(0, 729 * 24, by = 12),  # every 12h (BID)
    cmt  = 1,   # PROP_GUT compartment index
    amt  = prop_dose_amt,
    ii   = 12,
    addl = 0
)
e_propranolol <- do.call(ev, list(
    time = seq(0, 729*24, by=12),
    cmt  = 1,
    amt  = rep(prop_dose_amt, length(seq(0, 729*24, by=12)))
))

out_propranolol <- mod %>%
    param(AF_DRUG_CONC = 0) %>%
    mrgsim(events = e_propranolol, end = 730*24, delta = 24,
           outvars = c("FIBRO","HVPG","HVPG_calc","ALB_calc","BILI_calc",
                       "INR_calc","ASCITES","GFR_est","NH3_blood",
                       "HE_GRADE","CREAT","CHILD_PUGH","MELD_score",
                       "MELD_Na_score","VAR_BLEED_RISK","MORT_1YR",
                       "HR_obs","PROP_CONC")) %>%
    as_tibble() %>%
    mutate(Scenario = "Propranolol 40mg BID",
           time_days = time / 24)

#---------------------------------------------------------------------
# SCENARIO 3: Diuretic therapy — Spironolactone 100 mg + Furosemide 40 mg QD
# Spiro 100 mg = 100,000 ng
#---------------------------------------------------------------------
spiro_dose_amt <- 100e3 * 0.90  # ng

e_spiro <- do.call(ev, list(
    time = seq(0, 729*24, by=24),
    cmt  = 4,  # SPIRO_GUT
    amt  = rep(spiro_dose_amt, length(seq(0, 729*24, by=24)))
))

out_diuretic <- mod %>%
    param(AF_DRUG_CONC = 0) %>%
    mrgsim(events = e_spiro, end = 730*24, delta = 24,
           outvars = c("FIBRO","HVPG","HVPG_calc","ALB_calc","BILI_calc",
                       "INR_calc","ASCITES","GFR_est","NH3_blood",
                       "HE_GRADE","CREAT","CHILD_PUGH","MELD_score",
                       "MELD_Na_score","VAR_BLEED_RISK","MORT_1YR",
                       "SPIRO_CONC","ALDO_ACT")) %>%
    as_tibble() %>%
    mutate(Scenario = "Spironolactone 100mg QD",
           time_days = time / 24)

#---------------------------------------------------------------------
# SCENARIO 4: Combination — Propranolol + Spironolactone + Rifaximin 550 mg BID
# Rifaximin 550 mg = 550,000 ng
#---------------------------------------------------------------------
rifax_dose_amt <- 550e3 * 0.001  # minimal systemic; gut lumen treated differently

e_combo <- rbind(
    data.frame(time = seq(0, 729*24, by=12), cmt = 1, amt = prop_dose_amt),
    data.frame(time = seq(0, 729*24, by=24), cmt = 4, amt = spiro_dose_amt),
    data.frame(time = seq(0, 729*24, by=12), cmt = 8, amt = 550e3)  # RIFAX_GUT
)
e_combo <- do.call(ev, list(
    time = e_combo$time,
    cmt  = e_combo$cmt,
    amt  = e_combo$amt
))

out_combo <- mod %>%
    param(AF_DRUG_CONC = 0) %>%
    mrgsim(events = e_combo, end = 730*24, delta = 24,
           outvars = c("FIBRO","HVPG","HVPG_calc","ALB_calc","BILI_calc",
                       "INR_calc","ASCITES","GFR_est","NH3_blood",
                       "HE_GRADE","CREAT","CHILD_PUGH","MELD_score",
                       "MELD_Na_score","VAR_BLEED_RISK","MORT_1YR",
                       "HR_obs","PROP_CONC","SPIRO_CONC","RIFAX_CONC")) %>%
    as_tibble() %>%
    mutate(Scenario = "Propranolol + Spironolactone + Rifaximin",
           time_days = time / 24)

#---------------------------------------------------------------------
# SCENARIO 5: Antifibrotic therapy (investigational — FXR agonist or NASH treatment)
# AF_DRUG_CONC = 1 represents effective therapeutic concentration
#---------------------------------------------------------------------
out_antifib <- mod %>%
    param(AF_DRUG_CONC = 1.5) %>%
    mrgsim(events = e_natural, end = 730, delta = 1,
           outvars = c("FIBRO","HVPG","HVPG_calc","ALB_calc","BILI_calc",
                       "INR_calc","ASCITES","GFR_est","NH3_blood",
                       "HE_GRADE","CREAT","CHILD_PUGH","MELD_score",
                       "MELD_Na_score","VAR_BLEED_RISK","MORT_1YR",
                       "HSC_ACT")) %>%
    as_tibble() %>%
    mutate(Scenario = "Antifibrotic (FXR Agonist)",
           time_days = time)

#---------------------------------------------------------------------
# SCENARIO 6 (ACUTE): Terlipressin 2 mg IV q4h × 72h (acute variceal bleeding)
# or HRS-1 treatment
# Dose: 2 mg = 2,000,000 ng IV bolus
#---------------------------------------------------------------------
e_terl_acute <- do.call(ev, list(
    time = seq(0, 66, by=4),  # q4h for 72h
    cmt  = 6,  # TERL_C
    amt  = rep(2e6, length(seq(0, 66, by=4)))
))

out_acute_terl <- mod %>%
    param(AF_DRUG_CONC = 0) %>%
    mrgsim(events = e_terl_acute, end = 72, delta = 0.5,
           outvars = c("HVPG","HVPG_calc","GFR_est","CREAT","NH3_blood",
                       "ASCITES","TERL_CONC_EFF","VAR_BLEED_RISK")) %>%
    as_tibble() %>%
    mutate(Scenario = "Terlipressin IV (Acute HRS/Bleeding)",
           time_days = time / 24)

##############################################################################
# COMBINE CHRONIC SCENARIOS
##############################################################################

scenarios_chronic <- bind_rows(
    out_natural %>% mutate(time_days = time),
    out_antifib
) %>%
    select(time_days, Scenario, FIBRO, HVPG_calc, ALB_calc, BILI_calc,
           INR_calc, ASCITES, GFR_est, NH3_blood, HE_GRADE, CREAT,
           CHILD_PUGH, MELD_score, MELD_Na_score, VAR_BLEED_RISK, MORT_1YR)

##############################################################################
# VISUALIZATION
##############################################################################

theme_qsp <- theme_bw(base_size = 12) +
    theme(
        strip.background = element_rect(fill = "#2C3E50"),
        strip.text = element_text(color = "white", face = "bold"),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank()
    )

# ---- Plot 1: Fibrosis Progression ----
p1 <- scenarios_chronic %>%
    ggplot(aes(x = time_days, y = FIBRO, color = Scenario, linetype = Scenario)) +
    geom_line(size = 1.2) +
    geom_hline(yintercept = 0.75, linetype = "dotted", color = "red", alpha = 0.7) +
    annotate("text", x = 10, y = 0.77, label = "Cirrhosis threshold (F3-F4)", color = "red", size = 3) +
    scale_color_manual(values = c("Natural History" = "#D32F2F",
                                  "Antifibrotic (FXR Agonist)" = "#1976D2")) +
    labs(title = "Fibrosis Index over Time",
         x = "Time (days)", y = "Fibrosis Index (0-1)") +
    theme_qsp

# ---- Plot 2: HVPG over Time ----
p2 <- scenarios_chronic %>%
    ggplot(aes(x = time_days, y = HVPG_calc, color = Scenario)) +
    geom_line(size = 1.2) +
    geom_hline(yintercept = c(10, 12), linetype = "dashed",
               color = c("#FF9800", "#F44336"), alpha = 0.8) +
    annotate("text", x = 10, y = 10.3, label = "HVPG=10 (ascites risk)", size = 3, color = "#FF9800") +
    annotate("text", x = 10, y = 12.3, label = "HVPG=12 (bleeding risk)", size = 3, color = "#F44336") +
    scale_color_manual(values = c("Natural History" = "#D32F2F",
                                  "Antifibrotic (FXR Agonist)" = "#1976D2")) +
    labs(title = "HVPG Trajectory",
         x = "Time (days)", y = "HVPG (mmHg)") +
    theme_qsp

# ---- Plot 3: MELD Score Trajectory ----
p3 <- scenarios_chronic %>%
    ggplot(aes(x = time_days, y = MELD_score, color = Scenario)) +
    geom_line(size = 1.2) +
    geom_hline(yintercept = 15, linetype = "dashed", color = "#9C27B0") +
    annotate("text", x = 10, y = 15.5, label = "MELD=15 (transplant listing)", size = 3, color = "#9C27B0") +
    scale_color_manual(values = c("Natural History" = "#D32F2F",
                                  "Antifibrotic (FXR Agonist)" = "#1976D2")) +
    labs(title = "MELD Score Trajectory",
         x = "Time (days)", y = "MELD Score") +
    theme_qsp

# ---- Plot 4: Acute Terlipressin Effect ----
p4 <- out_acute_terl %>%
    select(time_days, GFR_est, CREAT, HVPG_calc) %>%
    pivot_longer(cols = -time_days) %>%
    ggplot(aes(x = time_days * 24, y = value, color = name)) +
    geom_line(size = 1.2) +
    facet_wrap(~name, scales = "free_y") +
    labs(title = "Terlipressin IV: Acute Phase Response (72h)",
         x = "Time (hours)", y = "Value") +
    scale_color_manual(values = c("GFR_est" = "#1976D2",
                                  "CREAT" = "#D32F2F",
                                  "HVPG_calc" = "#FF9800")) +
    theme_qsp + theme(legend.position = "none")

# ---- Plot 5: Multi-endpoint panel ----
p5 <- scenarios_chronic %>%
    select(time_days, Scenario, ALB_calc, BILI_calc, GFR_est, NH3_blood) %>%
    pivot_longer(cols = c(ALB_calc, BILI_calc, GFR_est, NH3_blood)) %>%
    mutate(name = recode(name,
        "ALB_calc" = "Albumin (g/dL)",
        "BILI_calc" = "Bilirubin (mg/dL)",
        "GFR_est"  = "eGFR (mL/min)",
        "NH3_blood" = "Ammonia (μmol/L)")) %>%
    ggplot(aes(x = time_days, y = value, color = Scenario)) +
    geom_line(size = 1) +
    facet_wrap(~name, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c("Natural History" = "#D32F2F",
                                  "Antifibrotic (FXR Agonist)" = "#1976D2")) +
    labs(title = "Hepatic Function & Systemic Biomarkers",
         x = "Time (days)", y = "Value") +
    theme_qsp

# ---- Save plots ----
if (!dir.exists("plots")) dir.create("plots")
ggsave("plots/lc_fibrosis_progression.png", p1, width = 10, height = 5, dpi = 150)
ggsave("plots/lc_hvpg_trajectory.png", p2, width = 10, height = 5, dpi = 150)
ggsave("plots/lc_meld_score.png", p3, width = 10, height = 5, dpi = 150)
ggsave("plots/lc_terlipressin_acute.png", p4, width = 10, height = 6, dpi = 150)
ggsave("plots/lc_biomarkers.png", p5, width = 10, height = 8, dpi = 150)

##############################################################################
# SUMMARY TABLE
##############################################################################

summary_table <- scenarios_chronic %>%
    filter(time_days %in% c(0, 90, 180, 365, 730)) %>%
    select(Scenario, time_days, FIBRO, HVPG_calc, ALB_calc, MELD_score,
           ASCITES, GFR_est, MORT_1YR) %>%
    mutate(
        across(where(is.numeric), ~round(., 2)),
        time_label = paste0("Day ", time_days)
    )

print(summary_table)

cat("\n====================================================\n")
cat("Liver Cirrhosis QSP Model — Scenario Summary\n")
cat("====================================================\n")
cat("Model compartments: 21 ODEs\n")
cat("Drug models: Propranolol (2-cpt), Spironolactone (1-cpt),\n")
cat("             Terlipressin (effect-cpt), Rifaximin (gut-lumen)\n")
cat("Key endpoints: MELD, Child-Pugh, HVPG, GFR, Ammonia, HE grade\n")
cat("Scenarios: Natural history, Antifibrotic, Propranolol,\n")
cat("           Diuretics, Combination, Acute terlipressin\n")
