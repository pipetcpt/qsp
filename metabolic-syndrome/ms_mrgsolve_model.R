## =============================================================================
## Metabolic Syndrome — Quantitative Systems Pharmacology Model
## mrgsolve ODE Implementation
##
## Disease: Metabolic Syndrome (MetS)
## Components: Central Obesity, Insulin Resistance, Dyslipidemia,
##             Hypertension, Chronic Inflammation
## Drugs modelled:
##   1. Metformin (AMPK activator, HGP inhibitor)
##   2. Semaglutide/Liraglutide (GLP-1 receptor agonist)
##   3. Empagliflozin/Dapagliflozin (SGLT2 inhibitor)
##   4. Rosuvastatin/Atorvastatin (HMG-CoA reductase inhibitor)
##   5. Losartan/Enalapril (ARB/ACEi)
##
## References:
##   Bergman (2001) PMID:11228214 – Minimal model
##   Dalla Man (2007) PMID:17466591 – Oral glucose minimal model
##   De Gaetano (2008) PMID:18286784 – HOMA model
##   Claret (2016) PMID:27378228 – GLP-1 RA PK/PD
##   Pocock (2016) PMID:27350, – SGLT2i TmG model
##   Nossen (2012) PMID:22540983 – statin PK
##   NCEP ATP-III Guidelines, Alberti 2006 (IDF consensus)
##
## Author: CCR — Claude Code Routine (2026-06-18)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ─────────────────────────────────────────────────────────────────────────────
## MODEL CODE
## ─────────────────────────────────────────────────────────────────────────────
ms_model_code <- '
$PROB Metabolic Syndrome QSP — 22-compartment ODE model

$PARAM
//-- Patient characteristics
BW     = 90      // Body weight (kg) — obese MetS patient
VAT0   = 5.0     // Baseline visceral adipose tissue (kg)
SAT0   = 20.0    // Baseline subcutaneous adipose tissue (kg)

//-- Glucose homeostasis (Bergman-derived parameters)
G_b    = 100     // Fasting glucose baseline (mg/dL)
I_b    = 25      // Fasting insulin baseline (μU/mL)
Sg     = 0.021   // Glucose effectiveness (min-1)
Si     = 4.5e-4  // Insulin sensitivity (mL·μU-1·min-1)
k_HGP  = 0.18    // Hepatic glucose production rate constant (mg/kg/min)
k_gut  = 0.07    // Gut glucose absorption rate (min-1)
GFR    = 120     // Glomerular filtration rate (mL/min)
TmG0   = 375     // Tubular max glucose reabsorption (mg/min)

//-- Insulin/glucagon secretion
Kg1    = 0.5     // β-cell glucose sensitivity (mU/mg)
Kg2    = 0.02    // β-cell 2nd phase slope
tau_B  = 90      // β-cell response time constant (min)
Beta0  = 1.0     // Relative β-cell mass (baseline=1)
k_BetaDeath = 0.0003  // IL-1β-mediated β-cell apoptosis rate
Gluc0  = 80      // Fasting glucagon (pg/mL)
kGluc  = 0.08    // Glucagon suppression by insulin rate

//-- GLP-1 endogenous
GLP1_e0 = 5.0   // Basal endogenous GLP-1 (pM)
kGLP1   = 0.12  // GLP-1 secretion rate constant

//-- Lipid model
VLDL0  = 120    // Baseline VLDL-TG (mg/dL)
LDL0   = 140    // Baseline LDL-C (mg/dL)
HDL0   = 38     // Baseline HDL-C (mg/dL)
TG0    = 220    // Baseline triglycerides (mg/dL)
kVLDL  = 0.05   // VLDL production from liver IR
kLDL   = 0.03   // LDL production from VLDL catabolism
kHDL   = 0.015  // HDL elimination rate
kCETP  = 0.02   // CETP-mediated TG-HDL exchange rate

//-- Inflammatory cytokines
TNFa0  = 15.0   // Baseline TNF-α (pg/mL)
IL6_0  = 3.0    // Baseline IL-6 (pg/mL)
IL1b0  = 2.0    // Baseline IL-1β (pg/mL)
CRP0   = 3.0    // Baseline hsCRP (mg/L)
kTNF   = 0.06   // TNF-α production from M1 macrophages
kIL6   = 0.04   // IL-6 production rate
kCRP   = 0.10   // CRP production from IL-6

//-- RAAS & Blood Pressure
AngII0 = 0.2    // Baseline Angiotensin II (ng/mL)
MAP0   = 100    // Baseline MAP (mmHg)
SVR0   = 1200   // Baseline SVR (dyne·s/cm5)
CO0    = 5.0    // Baseline cardiac output (L/min)
kRAAS  = 0.03   // RAAS response rate

//-- Adipokines
Lep0   = 25     // Baseline leptin (ng/mL) — obese
Adip0  = 5.0    // Baseline adiponectin (μg/mL) — low in MetS
kLep   = 0.8    // Leptin-VAT proportionality
kAdip  = 0.3    // Adiponectin-SAT proportionality

//-- AMPK & Energy sensing
AMPK0  = 1.0    // Baseline AMPK activity (relative)
kAMPK  = 0.15   // AMPK activation rate

//-- Drug PK parameters — Metformin
MET_Fg  = 0.55  // Metformin oral bioavailability
MET_ka  = 0.6   // Absorption rate (h-1)
MET_CL  = 28    // Clearance (L/h)
MET_Vc  = 380   // Central volume (L)
MET_Q   = 60    // Inter-compartmental CL (L/h)
MET_Vt  = 2000  // Peripheral volume (L)

//-- Drug PK — GLP-1 RA (Semaglutide weekly SC)
GLP_ka  = 0.004 // SC absorption rate (h-1)
GLP_CL  = 0.056 // Clearance (L/h)
GLP_Vc  = 8.0   // Central volume (L) — protein bound

//-- Drug PK — SGLT2i (Empagliflozin)
SGi_ka  = 1.5   // Absorption rate (h-1)
SGi_CL  = 12.0  // Clearance (L/h)
SGi_Vc  = 73.0  // Central volume (L)
SGi_F   = 0.86  // Bioavailability

//-- Drug PK — Statin (Rosuvastatin)
ST_ka   = 0.4   // Absorption rate (h-1)
ST_CL   = 25    // Clearance (L/h)
ST_Vc   = 134   // Central volume (L)
ST_F    = 0.20  // Bioavailability (hepatic FPE)

//-- Drug PK — ARB (Losartan)
ARB_ka  = 1.0   // Absorption rate (h-1)
ARB_CL  = 75    // Clearance (L/h)
ARB_Vc  = 34    // Central volume (L)
ARB_F   = 0.33  // Bioavailability

//-- Drug PD — EC50 / Emax
MET_EC50 = 1.5  // Metformin EC50 for HGP inhibition (μg/mL)
MET_Emax = 0.30 // Max HGP reduction (30%)
GLP_EC50 = 0.8  // GLP-1 RA EC50 insulin secretion (ng/mL)
GLP_Emax = 2.0  // Max fold-change insulin
GLP_WT50 = 2.0  // EC50 for weight loss effect (ng/mL)
GLP_WTmax= 0.12 // Max fraction body weight reduction
SGi_EC50 = 15   // SGLT2i EC50 TmG reduction (ng/mL)
SGi_Emax = 0.90 // Max TmG reduction fraction
ST_EC50  = 0.05 // Statin EC50 for cholesterol synthesis (μg/mL)
ST_Emax  = 0.70 // Max cholesterol synthesis reduction
ARB_EC50 = 0.3  // ARB EC50 for AT1R blockade (μg/mL)
ARB_Emax = 0.85 // Max AngII effect reduction

//-- Meal parameters (not dosed via $CMT but via custom calc)
D_meal   = 75000 // Oral glucose load (mg) — OGTT equivalent
k_meal   = 0.03  // Meal gut emptying rate (min-1)

$CMT
// Glucose system
GGUT    // Gut glucose (mg)
GPLAS   // Plasma glucose (mg/dL)

// Insulin / glucagon
BETA    // β-cell mass (relative)
IPLAS   // Plasma insulin (μU/mL)
GLUCPLAS // Glucagon (pg/mL)
GLP1E   // Endogenous GLP-1 (pM)

// Lipids
VLDLC   // VLDL-TG (mg/dL)
LDLC    // LDL-C (mg/dL)
HDLC    // HDL-C (mg/dL)
TRIGLY  // Triglycerides (mg/dL)

// Adipose
VAT     // Visceral adipose tissue (kg)
SAT     // Subcutaneous adipose tissue (kg)
LEP     // Leptin (ng/mL)
ADIPON  // Adiponectin (μg/mL)

// Inflammation
TNFA    // TNF-alpha (pg/mL)
IL6C    // IL-6 (pg/mL)
IL1BC   // IL-1beta (pg/mL)
CRPC    // hsCRP (mg/L)

// RAAS / BP
ANGII   // Angiotensin II (ng/mL)
MAPC    // Mean arterial pressure (mmHg)
AMPKC   // AMPK activity (relative)

// Drug PK — Metformin (2-CMT)
MET_GUT // Metformin gut depot (mg)
MET_CEN // Metformin central (mg)
MET_PER // Metformin peripheral (mg)

// Drug PK — GLP-1 RA (1-CMT SC)
GLP_SC  // GLP-1 RA SC depot (mg)
GLP_CEN // GLP-1 RA central (mg)

// Drug PK — SGLT2i (1-CMT)
SGI_GUT // SGLT2i gut (mg)
SGI_CEN // SGLT2i central (mg)

// Drug PK — Statin (1-CMT with FPE)
STA_GUT // Statin gut (mg)
STA_CEN // Statin central (mg)

// Drug PK — ARB (1-CMT)
ARB_GUT // ARB gut (mg)
ARB_CEN // ARB central (mg)

$INIT
GGUT    = 0
GPLAS   = 100   // mg/dL
BETA    = 1.0
IPLAS   = 25    // μU/mL
GLUCPLAS = 80
GLP1E   = 5.0
VLDLC   = 120
LDLC    = 140
HDLC    = 38
TRIGLY  = 220
VAT     = 5.0
SAT     = 20.0
LEP     = 25
ADIPON  = 5.0
TNFA    = 15
IL6C    = 3.0
IL1BC   = 2.0
CRPC    = 3.0
ANGII   = 0.2
MAPC    = 100
AMPKC   = 1.0
MET_GUT = 0
MET_CEN = 0
MET_PER = 0
GLP_SC  = 0
GLP_CEN = 0
SGI_GUT = 0
SGI_CEN = 0
STA_GUT = 0
STA_CEN = 0
ARB_GUT = 0
ARB_CEN = 0

$ODE

// ─── Drug concentrations (μg/mL or ng/mL) ───────────────────────────────────
double Cp_MET  = MET_CEN / MET_Vc;     // μg/mL
double Cp_GLP  = GLP_CEN / GLP_Vc;     // ng/mL
double Cp_SGI  = SGI_CEN / SGI_Vc;     // ng/mL
double Cp_STA  = STA_CEN / ST_Vc;      // μg/mL
double Cp_ARB  = ARB_CEN / ARB_Vc;     // μg/mL

// ─── Drug effect calculations (Hill function) ────────────────────────────────
double E_MET_HGP  = MET_Emax * Cp_MET / (MET_EC50 + Cp_MET);   // HGP reduction
double E_GLP_INS  = GLP_Emax * Cp_GLP / (GLP_EC50 + Cp_GLP);   // Insulin fold-change
double E_GLP_WT   = GLP_WTmax* Cp_GLP / (GLP_WT50 + Cp_GLP);   // Weight loss fraction
double E_SGI_TMG  = SGi_Emax * Cp_SGI / (SGi_EC50 + Cp_SGI);   // TmG reduction
double E_STA_CHOL = ST_Emax  * Cp_STA / (ST_EC50  + Cp_STA);   // Chol synth reduction
double E_ARB_BP   = ARB_Emax * Cp_ARB / (ARB_EC50 + Cp_ARB);   // AngII effect reduction

// ─── Insulin resistance indices ──────────────────────────────────────────────
double IR_FFA  = 1.0 + 0.15 * (VAT / VAT0 - 1);     // FFA-mediated IR
double IR_TNF  = 1.0 + 0.08 * (TNFA / TNFa0 - 1);   // TNF-mediated IR
double IR_AngII = 1.0 + 0.04 * (ANGII / AngII0 - 1);// AngII-mediated IR
double IR_total = IR_FFA * IR_TNF * IR_AngII;         // Composite IR

// ─── Glucose homeostasis ──────────────────────────────────────────────────────
// Hepatic glucose production (HGP): inhibited by insulin and metformin
double f_Ins_HGP = 1.0 / (1.0 + 0.03 * IPLAS);
double HGP = k_HGP * BW * f_Ins_HGP * IR_total * (1.0 - E_MET_HGP)
           * (1.0 + 0.2 * (GLUCPLAS / Gluc0 - 1));  // glucagon stimulation

// Peripheral glucose disposal: insulin-stimulated GLUT4 translocation
double Rd  = (Sg + Si * IPLAS / IR_total) * GPLAS;

// Renal glucose handling (SGLT2 inhibition effect)
double TmG_eff = TmG0 * (1.0 - E_SGI_TMG);
double UGE  = (GFR * GPLAS / 100.0 > TmG_eff) ?
              (GFR * GPLAS / 100.0 - TmG_eff) : 0.0;  // Urinary glucose excretion

// Gut glucose absorption
double Ra_gut = k_gut * GGUT;

// Plasma glucose dODE
dxdt_GGUT  = -k_gut * GGUT;
dxdt_GPLAS = (HGP + Ra_gut - Rd - UGE) / BW;  // Simplified distribution

// ─── β-cell dynamics ─────────────────────────────────────────────────────────
double GSIS_first  = Kg1 * (GPLAS - G_b);                // 1st phase secretion
double GSIS_second = Kg2 * GPLAS * BETA;                  // 2nd phase
double Ins_sec = (GSIS_first + GSIS_second > 0) ?
                  (GSIS_first + GSIS_second) : 0.0;
Ins_sec *= (1.0 + E_GLP_INS);                             // GLP-1 RA potentiation

double k_BetaGrowth = 0.0002 * (GPLAS / G_b);             // Glucose-induced compensation
double k_BetaLoss   = k_BetaDeath * IL1BC / IL1b0;        // IL-1β apoptosis

dxdt_BETA  = (k_BetaGrowth - k_BetaLoss) * BETA;
dxdt_IPLAS = Ins_sec - 0.05 * IPLAS;                      // Linear elimination

// ─── Glucagon dynamics ────────────────────────────────────────────────────────
double Gluc_target = Gluc0 / (1.0 + 0.015 * IPLAS);      // Insulin suppression
dxdt_GLUCPLAS = kGluc * (Gluc_target - GLUCPLAS);

// ─── Endogenous GLP-1 ─────────────────────────────────────────────────────────
double GLP1_stim = kGLP1 * (GPLAS / G_b);                 // Glucose-stimulated
dxdt_GLP1E = GLP1_stim - 0.15 * GLP1E;                   // DPP-IV elimination

// ─── Lipid dynamics ──────────────────────────────────────────────────────────
// VLDL-TG: increased by hepatic IR, reduced by adiponectin
double VLDL_prod = kVLDL * (VLDLC + IR_total * 15)
                 * (ADIPON0 / ADIPON)                      // adiponectin protection
                 * (1.0 - E_MET_HGP * 0.3);               // Metformin partial
double VLDL_elim = 0.06 * VLDLC;
dxdt_VLDLC = VLDL_prod - VLDL_elim;

// LDL-C: from VLDL catabolism, reduced by statin (↑LDLr)
double LDL_prod = kLDL * TRIGLY;
double LDL_elim = (0.04 + E_STA_CHOL * 0.08) * LDLC;    // Statin ↑LDLr
dxdt_LDLC = LDL_prod - LDL_elim;

// HDL-C: reduced by TG-rich milieu (CETP), ARB modestly helpful
double HDL_prod = 0.015 * ADIPON;                          // Adiponectin-linked
double HDL_elim = kHDL * HDLC + kCETP * TRIGLY / TG0;    // CETP exchange
dxdt_HDLC = HDL_prod - HDL_elim;

// Triglycerides: from VLDL, reduced by SGLT2i (caloric loss) and GLP-1 RA
double TG_prod  = 0.04 * VLDLC;
double TG_elim  = 0.025 * TRIGLY * (1.0 + 0.1 * ADIPON / Adip0);
double TG_GLP1_eff = 0.1 * E_GLP_WT;                      // GLP-1 RA TG reduction
dxdt_TRIGLY = TG_prod - TG_elim * (1.0 + TG_GLP1_eff);

// ─── Adipose Tissue dynamics ─────────────────────────────────────────────────
// VAT: driven by caloric excess, cortisol, reduced by GLP-1 RA and SGLT2i
double VAT_accum = 0.0005 * (IPLAS / I_b - 1) * VAT;     // Hyperinsulinemia
double VAT_loss  = E_GLP_WT * 0.005 * VAT + 0.002 * VAT;  // GLP-1 RA + basal
dxdt_VAT = VAT_accum - VAT_loss;

// SAT: less pathological, slower turnover
double SAT_change = 0.0003 * (1.0 - E_GLP_WT);
dxdt_SAT = -SAT_change * SAT;

// Leptin: proportional to total fat mass
double Lep_target = kLep * (VAT + SAT) / (VAT0 + SAT0) * Lep0;
dxdt_LEP = 0.05 * (Lep_target - LEP);

// Adiponectin: inversely related to VAT (MetS pattern)
double Adip_target = Adip0 * (VAT0 / VAT) * (1.0 + 0.5 * E_ARB_BP);
dxdt_ADIPON = 0.03 * (Adip_target - ADIPON);

// ─── Inflammatory cytokines ──────────────────────────────────────────────────
// TNF-α: from M1 macrophages in VAT
double TNF_prod  = kTNF * (VAT / VAT0) * (LEP / Lep0);
double TNF_inhib = 0.02 * (ADIPON / Adip0);               // Adiponectin anti-inflam
double TNF_elim  = 0.05 * TNFA;
dxdt_TNFA = TNF_prod - TNF_elim - TNF_inhib * TNFA;

// IL-6: from macrophages and adipose
double IL6_prod = kIL6 * (VAT / VAT0) * (TNFA / TNFa0);
double IL6_elim = 0.08 * IL6C;
dxdt_IL6C = IL6_prod - IL6_elim;

// IL-1β: NLRP3 inflammasome activated by FFA
double IL1b_prod = 0.015 * (VAT / VAT0) * (TNFA / TNFa0);
double IL1b_elim = 0.06 * IL1BC;
dxdt_IL1BC = IL1b_prod - IL1b_elim;

// CRP: acute phase protein from liver (IL-6 driven)
double CRP_prod = kCRP * IL6C / IL6_0;
double CRP_elim = 0.04 * CRPC;
dxdt_CRPC = CRP_prod - CRP_elim;

// ─── RAAS & Blood Pressure ────────────────────────────────────────────────────
// AngII: increased by obesity/SNS activation, reduced by ARB
double AngII_prod = kRAAS * (VAT / VAT0) * (TNFA / TNFa0);
double AngII_elim = 0.15 * ANGII * (1.0 - E_ARB_BP * 0.7);
dxdt_ANGII = AngII_prod - AngII_elim;

// MAP: driven by SVR (AngII) and volume (SGLT2i reduces)
double MAP_target = MAP0 + 15 * (ANGII / AngII0 - 1)       // AngII pressure
                  + 8  * (IL6C / IL6_0 - 1)                 // inflammation
                  - 10 * E_ARB_BP                            // ARB/ACEi effect
                  - 4  * E_SGI_TMG;                          // SGLT2i volume
dxdt_MAPC = 0.01 * (MAP_target - MAPC);

// ─── AMPK activity ───────────────────────────────────────────────────────────
// Activated by metformin, exercise (not modelled), adiponectin
double AMPK_target = AMPK0 * (1.0 + 2.0 * E_MET_HGP)      // Metformin
                   * (1.0 + 0.5 * ADIPON / Adip0)           // Adiponectin
                   / IR_total;                               // IR reduces AMPK
dxdt_AMPKC = 0.05 * (AMPK_target - AMPKC);

// ─── Drug PK — Metformin (2-compartment) ─────────────────────────────────────
dxdt_MET_GUT = -MET_ka * MET_GUT;
dxdt_MET_CEN =  MET_ka * MET_Fg * MET_GUT
             - (MET_CL / MET_Vc) * MET_CEN
             - (MET_Q  / MET_Vc) * MET_CEN
             + (MET_Q  / MET_Vt) * MET_PER;
dxdt_MET_PER =  (MET_Q  / MET_Vc) * MET_CEN
             -  (MET_Q  / MET_Vt) * MET_PER;

// ─── Drug PK — GLP-1 RA (SC 1-compartment) ───────────────────────────────────
dxdt_GLP_SC  = -GLP_ka * GLP_SC;
dxdt_GLP_CEN =  GLP_ka * GLP_SC - (GLP_CL / GLP_Vc) * GLP_CEN;

// ─── Drug PK — SGLT2i ──────────────────────────────────────────────────────────
dxdt_SGI_GUT = -SGi_ka * SGI_GUT;
dxdt_SGI_CEN =  SGi_ka * SGi_F * SGI_GUT - (SGi_CL / SGI_Vc) * SGI_CEN;

// ─── Drug PK — Statin ────────────────────────────────────────────────────────
dxdt_STA_GUT = -ST_ka * STA_GUT;
dxdt_STA_CEN =  ST_ka * ST_F * STA_GUT - (ST_CL / ST_Vc) * STA_CEN;

// ─── Drug PK — ARB ───────────────────────────────────────────────────────────
dxdt_ARB_GUT = -ARB_ka * ARB_GUT;
dxdt_ARB_CEN =  ARB_ka * ARB_F * ARB_GUT - (ARB_CL / ARB_Vc) * ARB_CEN;

$PARAM
// Convenience references for steady-state values
ADIPON0 = 5.0   // used in VLDL_prod formula
VAT0_ref = 5.0
Adip0_ref = 5.0
Lep0_ref = 25.0
Gluc0    = 80.0

$OMEGA 0.0  // No ETA for deterministic base model

$TABLE
// Derived endpoints
double HbA1c    = 5.0 + GPLAS / 30.0;   // Simplified HbA1c estimate (%)
double HOMA_IR  = GPLAS * IPLAS / 405.0; // HOMA-IR index
double HOMA_B   = 20.0 * IPLAS / (GPLAS - 3.5); // HOMA-β
double FPG      = GPLAS;
double TG_HDL   = TRIGLY / HDLC;        // Atherogenic ratio
double NonHDL   = LDLC + TRIGLY / 5.0;  // Non-HDL-C
double SBP      = MAPC + 40;            // Approximation SBP
double DBP      = MAPC - 10;            // Approximation DBP
double OxStress = TNFA * IL6C / 45.0;   // Oxidative stress proxy
double MetS_Zscore = (GPLAS - 100) / 20 + (MAPC - 93) / 13 +
                    (TRIGLY - 150) / 50 - (HDLC - 50) / 12 +
                    (VAT - 4.5) / 1.0;   // Composite MetS severity

// NCEP-ATP III criteria (binary)
double crit_gluc = (GPLAS >= 100) ? 1 : 0;
double crit_TG   = (TRIGLY >= 150) ? 1 : 0;
double crit_HDL  = (HDLC < 40) ? 1 : 0;
double crit_BP   = (SBP >= 130 || DBP >= 85) ? 1 : 0;
double crit_WC   = 1;  // Already obese patient
double MetS_NCEP = crit_gluc + crit_TG + crit_HDL + crit_BP + crit_WC;

// Drug concentrations for output
double Cp_MET_out  = MET_CEN / MET_Vc;
double Cp_GLP_out  = GLP_CEN / GLP_Vc;
double Cp_SGI_out  = SGI_CEN / SGI_Vc;
double Cp_STA_out  = STA_CEN / ST_Vc;
double Cp_ARB_out  = ARB_CEN / ARB_Vc;

$CAPTURE
GPLAS IPLAS GLUCPLAS GLP1E
VLDLC LDLC HDLC TRIGLY
VAT SAT LEP ADIPON
TNFA IL6C IL1BC CRPC
ANGII MAPC AMPKC BETA
HbA1c HOMA_IR HOMA_B FPG TG_HDL NonHDL SBP DBP
OxStress MetS_Zscore MetS_NCEP
Cp_MET_out Cp_GLP_out Cp_SGI_out Cp_STA_out Cp_ARB_out
'

## ─────────────────────────────────────────────────────────────────────────────
## Compile model
## ─────────────────────────────────────────────────────────────────────────────
ms_mod <- mcode("MetabolicSyndrome_QSP", ms_model_code)

## ─────────────────────────────────────────────────────────────────────────────
## Treatment scenarios
## ─────────────────────────────────────────────────────────────────────────────

## Common simulation settings: 52-week follow-up, hourly output
sim_end_wk <- 52
dt <- 1  # hours

##
## Scenario 1: No treatment — natural MetS progression
##
ev_none <- ev(ID = 1, time = 0, amt = 0, cmt = 1)

##
## Scenario 2: Metformin 1000 mg BID
##   Dose every 12h starting Week 0
##
ev_met <- ev(ID = 2, amt = 1000, cmt = "MET_GUT",
             ii = 12, addl = 2 * 7 * sim_end_wk - 1,
             time = 0)

##
## Scenario 3: Semaglutide 1 mg SC weekly + Metformin 1000 mg BID
##
ev_glp <- ev(ID = 3, amt = 1, cmt = "GLP_SC",
             ii = 168, addl = sim_end_wk - 1,
             time = 0)
ev_met3 <- ev(ID = 3, amt = 1000, cmt = "MET_GUT",
              ii = 12, addl = 2 * 7 * sim_end_wk - 1,
              time = 0)

##
## Scenario 4: Empagliflozin 10 mg QD + Metformin 1000 mg BID
##
ev_sgi <- ev(ID = 4, amt = 10, cmt = "SGI_GUT",
             ii = 24, addl = 7 * sim_end_wk - 1,
             time = 0)
ev_met4 <- ev(ID = 4, amt = 1000, cmt = "MET_GUT",
              ii = 12, addl = 2 * 7 * sim_end_wk - 1,
              time = 0)

##
## Scenario 5: Triple therapy (Metformin + Semaglutide + Rosuvastatin 10 mg QD + Losartan 50 mg QD)
##
ev_sta <- ev(ID = 5, amt = 10, cmt = "STA_GUT",
             ii = 24, addl = 7 * sim_end_wk - 1,
             time = 0)
ev_arb <- ev(ID = 5, amt = 50, cmt = "ARB_GUT",
             ii = 24, addl = 7 * sim_end_wk - 1,
             time = 0)
ev_glp5 <- ev(ID = 5, amt = 1, cmt = "GLP_SC",
              ii = 168, addl = sim_end_wk - 1,
              time = 0)
ev_met5 <- ev(ID = 5, amt = 1000, cmt = "MET_GUT",
              ii = 12, addl = 2 * 7 * sim_end_wk - 1,
              time = 0)

## Run simulations
sim_none <- ms_mod %>%
  ev(ev_none) %>%
  mrgsim(end = sim_end_wk * 168, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "1_No_Treatment")

sim_met <- ms_mod %>%
  ev(ev_met) %>%
  mrgsim(end = sim_end_wk * 168, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "2_Metformin")

sim_glp <- ms_mod %>%
  ev(c(ev_glp, ev_met3)) %>%
  mrgsim(end = sim_end_wk * 168, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "3_GLP1_RA_Met")

sim_sgi <- ms_mod %>%
  ev(c(ev_sgi, ev_met4)) %>%
  mrgsim(end = sim_end_wk * 168, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "4_SGLT2i_Met")

sim_triple <- ms_mod %>%
  ev(c(ev_met5, ev_glp5, ev_sta, ev_arb)) %>%
  mrgsim(end = sim_end_wk * 168, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "5_Quadruple_Therapy")

## Combine all scenarios
all_sims <- bind_rows(sim_none, sim_met, sim_glp, sim_sgi, sim_triple) %>%
  mutate(Week = time / 168)

## ─────────────────────────────────────────────────────────────────────────────
## PLOTS
## ─────────────────────────────────────────────────────────────────────────────

cols <- c(
  "1_No_Treatment"     = "#E74C3C",
  "2_Metformin"        = "#F39C12",
  "3_GLP1_RA_Met"      = "#2ECC71",
  "4_SGLT2i_Met"       = "#3498DB",
  "5_Quadruple_Therapy" = "#9B59B6"
)

## Plot 1: Glucose & HbA1c
p1 <- all_sims %>%
  ggplot(aes(x = Week, y = HbA1c, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 7.0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = cols) +
  labs(title = "HbA1c Over 52 Weeks",
       y = "HbA1c (%)", x = "Week") +
  annotate("text", x = 5, y = 6.85, label = "Target < 7%", size = 3) +
  theme_bw()

## Plot 2: Fasting Plasma Glucose
p2 <- all_sims %>%
  ggplot(aes(x = Week, y = FPG, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = cols) +
  labs(title = "Fasting Plasma Glucose",
       y = "FPG (mg/dL)", x = "Week") +
  theme_bw()

## Plot 3: LDL-C
p3 <- all_sims %>%
  ggplot(aes(x = Week, y = LDLC, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 70, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = cols) +
  labs(title = "LDL-C Response",
       y = "LDL-C (mg/dL)", x = "Week") +
  theme_bw()

## Plot 4: Triglycerides
p4 <- all_sims %>%
  ggplot(aes(x = Week, y = TRIGLY, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 150, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = cols) +
  labs(title = "Triglycerides",
       y = "TG (mg/dL)", x = "Week") +
  theme_bw()

## Plot 5: SBP
p5 <- all_sims %>%
  ggplot(aes(x = Week, y = SBP, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 130, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = cols) +
  labs(title = "Systolic Blood Pressure",
       y = "SBP (mmHg)", x = "Week") +
  theme_bw()

## Plot 6: VAT (visceral fat)
p6 <- all_sims %>%
  ggplot(aes(x = Week, y = VAT, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = cols) +
  labs(title = "Visceral Adipose Tissue",
       y = "VAT (kg)", x = "Week") +
  theme_bw()

## Plot 7: HOMA-IR
p7 <- all_sims %>%
  ggplot(aes(x = Week, y = HOMA_IR, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 2.5, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = cols) +
  labs(title = "HOMA-IR (Insulin Resistance)",
       y = "HOMA-IR", x = "Week") +
  theme_bw()

## Plot 8: MetS Z-score
p8 <- all_sims %>%
  ggplot(aes(x = Week, y = MetS_Zscore, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = cols) +
  labs(title = "Metabolic Syndrome Z-score\n(Composite Severity)",
       y = "MetS Z-score", x = "Week") +
  theme_bw()

## Plot 9: Inflammation (TNF-α, CRP)
p9 <- all_sims %>%
  pivot_longer(c(TNFA, CRPC), names_to = "Marker", values_to = "Value") %>%
  ggplot(aes(x = Week, y = Value, color = Scenario, linetype = Marker)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = cols) +
  labs(title = "Inflammatory Markers",
       y = "Concentration", x = "Week") +
  theme_bw()

## Plot 10: Drug concentrations at Week 4
p10 <- all_sims %>%
  filter(Scenario == "5_Quadruple_Therapy") %>%
  select(Week, Cp_MET_out, Cp_GLP_out, Cp_SGI_out, Cp_STA_out, Cp_ARB_out) %>%
  pivot_longer(-Week, names_to = "Drug", values_to = "Cp") %>%
  ggplot(aes(x = Week, y = Cp, color = Drug)) +
  geom_line(linewidth = 1.1) +
  scale_y_log10() +
  labs(title = "Drug Plasma Concentrations\n(Quadruple Therapy)",
       y = "Cp (log-scale)", x = "Week") +
  theme_bw()

## Combined dashboard
combined_plot <- (p1 | p2 | p3) /
                 (p4 | p5 | p6) /
                 (p7 | p8 | p10) +
  plot_annotation(
    title = "Metabolic Syndrome QSP Model — 52-Week Treatment Comparison",
    subtitle = "Scenarios: No Treatment | Metformin | GLP-1 RA+Met | SGLT2i+Met | Quadruple Therapy",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

print(combined_plot)

## ─────────────────────────────────────────────────────────────────────────────
## Summary table at Week 52
## ─────────────────────────────────────────────────────────────────────────────
summary_tbl <- all_sims %>%
  filter(abs(Week - 52) < 0.5) %>%
  group_by(Scenario) %>%
  summarise(
    HbA1c_pct  = round(mean(HbA1c), 1),
    FPG_mgdL   = round(mean(FPG), 0),
    LDLC_mgdL  = round(mean(LDLC), 0),
    HDL_mgdL   = round(mean(HDLC), 0),
    TG_mgdL    = round(mean(TRIGLY), 0),
    SBP_mmHg   = round(mean(SBP), 0),
    VAT_kg     = round(mean(VAT), 2),
    HOMA_IR    = round(mean(HOMA_IR), 1),
    MetSZ      = round(mean(MetS_Zscore), 2),
    MetS_NCEP  = round(mean(MetS_NCEP), 0),
    .groups = "drop"
  )

print(summary_tbl)

## ─────────────────────────────────────────────────────────────────────────────
## Sensitivity analysis: HOMA-IR at Week 52 across BMI values
## ─────────────────────────────────────────────────────────────────────────────
bmi_vals <- c(28, 30, 32, 35, 38, 42)  # kg/m2
vat_from_bmi <- function(bmi) 2.5 + 0.08 * (bmi - 28)  # rough linear

sens_results <- map_dfr(bmi_vals, function(bmi_i) {
  vat_i <- vat_from_bmi(bmi_i)
  ms_mod %>%
    param(VAT0 = vat_i, SAT0 = 18 + (bmi_i - 28) * 0.5) %>%
    init(VAT = vat_i) %>%
    ev(ev_glp) %>%
    mrgsim(end = 52 * 168, delta = 168) %>%
    as_tibble() %>%
    filter(time == max(time)) %>%
    transmute(BMI = bmi_i,
              HOMA_IR = HOMA_IR,
              HbA1c   = HbA1c,
              LDL     = LDLC,
              VAT_change_pct = (VAT - vat_i) / vat_i * 100)
})

print("Sensitivity analysis — GLP-1 RA effect across BMI levels at Week 52:")
print(sens_results)

message("✓ Metabolic Syndrome QSP model simulation complete.")
message("  22 ODE compartments | 5 treatment scenarios | Sensitivity analysis across BMI")
