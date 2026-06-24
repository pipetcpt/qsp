##############################################################
#  Heart Failure with Preserved Ejection Fraction (HFpEF)
#  QSP / PK-PD mrgsolve Model
#  Language: R + mrgsolve
#  Author: Claude Code Routine (CCR)
#  Date: 2026-06-17
#
#  Key clinical references:
#  - EMPEROR-Preserved (Empagliflozin, NEJM 2021)
#  - DELIVER (Dapagliflozin, NEJM 2022)
#  - PARAGON-HF (Sacubitril/Valsartan, NEJM 2019)
#  - TOPCAT (Spironolactone, NEJM 2014)
##############################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

# ============================================================
# MODEL DEFINITION
# ============================================================

hfpef_model_code <- '
$PROB
HFpEF Quantitative Systems Pharmacology Model
Compartments: Drug PK (4 drugs) + Neurohumoral (RAAS/SNS) +
              Cardiac Structure + Hemodynamics + Inflammation +
              Renal + Metabolic + Clinical Endpoints

$PARAM
// ── Patient demographics ──────────────────────────────────
Age      = 72      // years
Weight   = 95      // kg (typical obese HFpEF patient)
BSA      = 2.1     // m²
BMI      = 33      // kg/m²

// ── Baseline comorbidity scores ───────────────────────────
base_HTN     = 1   // 0/1 flag
base_T2DM    = 1
base_AF      = 0
base_CKD_sev = 0.3 // 0=none to 1=severe (eGFR <30)

// ── Empagliflozin PK (2-compartment, oral) ───────────────
empa_F    = 0.78   // Bioavailability
empa_ka   = 1.5    // Absorption rate (1/h)
empa_CL   = 9.4    // Clearance (L/h)
empa_V1   = 73.8   // Central volume (L)
empa_Q    = 18.0   // Inter-compartment CL (L/h)
empa_V2   = 36.0   // Peripheral volume (L)
empa_EC50 = 5.0    // SGLT2 inhibition EC50 (ng/mL)
empa_Emax = 1.0    // Max inhibition

// ── Sacubitril/Valsartan PK ──────────────────────────────
sac_F_sac = 0.60   // Sacubitril bioavailability
sac_ka    = 1.0    // Absorption rate (1/h)
sac_CL_sac = 4.5  // CL sacubitrilat (L/h)
sac_V_sac  = 35.0  // V sacubitrilat (L)
sac_kconv  = 0.25  // Sacubitril → Sacubitrilat (1/h)
val_F     = 0.23   // Valsartan bioavailability
val_CL    = 1.3    // Valsartan CL (L/h)
val_V     = 17.0   // Valsartan V (L)
val_EC50  = 0.10   // AT1R blockade EC50 (ng/mL Valsartan)
sac_EC50_NEP = 5.0 // NEP inhibition EC50 (ng/mL Sacubitrilat)

// ── Finerenone (MRA) PK ──────────────────────────────────
fin_F     = 0.43   // Bioavailability
fin_ka    = 1.2    // Absorption (1/h)
fin_CL    = 3.2    // CL (L/h)
fin_V     = 52.0   // V (L)
fin_EC50  = 0.9    // MR blockade EC50 (ng/mL)
fin_Emax  = 0.95   // Max MR blockade

// ── Furosemide PK ────────────────────────────────────────
furo_F    = 0.60   // Oral bioavailability
furo_ka   = 2.0    // Absorption rate (1/h)
furo_CL   = 8.0    // CL (L/h)
furo_V    = 12.0   // V (L)
furo_Vt   = 0.35   // Tubular secretion CL fraction
furo_EC50 = 2.0    // Natriuresis EC50 (μg/mL tubular)
furo_Emax = 0.80   // Max natriuresis effect

// ── RAAS baseline & dynamics ─────────────────────────────
base_AngII    = 1.0    // normalized (1 = normal)
base_Aldo     = 1.0    // normalized
kout_AngII    = 0.5    // degradation (1/h)
kin_AngII     = 0.5    // synthesis (normalized units/h)
AngII_Aldo_EC50 = 1.5  // AngII for 50% Aldo stimulation

// ── ANP/BNP / Natriuretic peptide system ─────────────────
base_ANP   = 1.0   // normalized (1 = normal; HFpEF ~3-5x elevated)
base_BNP   = 1.0
kout_ANP   = 0.35  // degradation (1/h, by NEP & NPR-C)
kout_BNP   = 0.20  // degradation (1/h)
kin_ANP    = 0.35
kin_BNP    = 0.20
LVEDP_ANP_EC50 = 20  // LVEDP (mmHg) for 50% ANP stimulation
BNP_hill   = 2.0   // Hill coefficient for pressure-ANP/BNP

// ── Cardiac structure ─────────────────────────────────────
base_LVM_idx  = 110  // LV mass index (g/m²); normal <95 (F), <115 (M)
base_fibrosis = 0.30 // Fractional fibrosis (0-1); HFpEF typical 0.25-0.40
base_Ecol     = 1.3  // Collagen I/III ratio (elevated in HFpEF)
kgrowth_LVM   = 0.001  // LVM hypertrophy rate per unit AngII×ET1
kdecay_LVM    = 0.002  // LVM regression rate
kfib_up       = 0.005  // Fibrosis accrual per unit TGFb
kfib_down     = 0.003  // Fibrosis resolution rate

// ── Titin stiffness & LV filling ─────────────────────────
base_titin_stiff = 0.45  // 0=soft, 1=maximal stiffness
ktitin_PKG       = 0.10  // PKG phosphorylation softening rate
ktitin_age       = 0.001 // Age-related stiffening per year above 60

// ── Hemodynamics ─────────────────────────────────────────
base_HR    = 78    // bpm
base_SV    = 70    // mL
base_SVR   = 1400  // dyne·s·cm⁻⁵
base_LVEDP = 16    // mmHg (elevated in HFpEF; normal <12)
base_PCWP  = 18    // mmHg
LVEDP_sv_sens = -0.3  // SV sensitivity to LVEDP change (mL/mmHg excess)

// ── Renal dynamics ───────────────────────────────────────
base_GFR     = 62  // mL/min/1.73m²  (CKD stage 2-3 typical in HFpEF)
base_Na_exc  = 1.0 // Normalized Na excretion
kGFR_aldo    = 0.05 // Aldosterone effect on GFR (subtle)
base_UricAcid = 7.2  // mg/dL (elevated in HFpEF)

// ── Inflammation biomarkers ───────────────────────────────
base_CRP    = 3.5   // mg/L (hsCRP; HFpEF > 2 mg/L typical)
base_IL6    = 1.0   // normalized (1=normal ~2 pg/mL)
base_TNFa   = 1.0   // normalized
base_sST2   = 1.0   // normalized (HFpEF elevated sST2 ~40 ng/mL)
kout_CRP    = 0.02  // CRP half-life ~19h → kout=0.036/h
kin_CRP     = 0.07  // hepatic synthesis stimulated by IL-6

// ── cGMP-PKG pathway ─────────────────────────────────────
base_cGMP = 1.0    // normalized
base_PKG   = 1.0   // normalized (1 = normal PKG activity)
kout_cGMP  = 0.3   // PDE5/9 degradation
kin_cGMP   = 0.3   // ANP/BNP/NO stimulated synthesis
PKG_ox_EC50 = 2.0  // Oxidative stress EC50 for PKG inactivation

// ── Treatment flags ──────────────────────────────────────
use_empa  = 0  // 1 = empagliflozin 10mg QD
use_arni  = 0  // 1 = sacubitril/valsartan 97/103mg BID
use_fin   = 0  // 1 = finerenone 20mg QD
use_furo  = 0  // 1 = furosemide 40mg BID

$CMT
// Empagliflozin PK
EMPA_GUT EMPA_C EMPA_P

// ARNI PK
SAC_GUT SAC_C VAL_GUT VAL_C

// Finerenone PK
FIN_GUT FIN_C

// Furosemide PK
FURO_GUT FURO_C

// Neurohumoral
ANGII ALDO

// Natriuretic peptides
ANP BNP

// Cardiac structure
LVM_IDX FIBROSIS TITIN_STIFF

// Hemodynamics (state variables tracking deviation from baseline)
LVEDP_dyn SVR_dyn

// Inflammation
CRP_dyn IL6_dyn

// cGMP-PKG pathway
cGMP_dyn PKG_dyn

// Renal
GFR_dyn Na_exc_dyn

// Clinical composite
NT_proBNP_dyn   // NT-proBNP (pg/mL)

$MAIN
// Empagliflozin steady-state dosing (QD → ADDL+II handle repeat)
F_EMPA_GUT = empa_F;

// ARNI
F_SAC_GUT = sac_F_sac;
F_VAL_GUT = val_F;

// MRA
F_FIN_GUT = fin_F;

// Furosemide
F_FURO_GUT = furo_F;

$ODE
//============================================================
// 1. EMPAGLIFLOZIN PK (2-compartment)
//============================================================
double empa_Cp = EMPA_C / empa_V1;  // ng/mL

dxdt_EMPA_GUT = -empa_ka * EMPA_GUT;
dxdt_EMPA_C   =  empa_ka * EMPA_GUT
                - (empa_CL / empa_V1) * EMPA_C
                - (empa_Q  / empa_V1) * EMPA_C
                + (empa_Q  / empa_V2) * EMPA_P;
dxdt_EMPA_P   =  (empa_Q  / empa_V1) * EMPA_C
                - (empa_Q  / empa_V2) * EMPA_P;

// SGLT2 inhibition effect (Emax model)
double SGLT2_inh = empa_Emax * empa_Cp / (empa_EC50 + empa_Cp);
SGLT2_inh = use_empa * SGLT2_inh;

//============================================================
// 2. SACUBITRIL/VALSARTAN PK
//============================================================
double sac_Cp  = SAC_C  / sac_V_sac;   // Sacubitrilat (ng/mL)
double val_Cp  = VAL_C  / val_V;        // Valsartan (ng/mL)

dxdt_SAC_GUT = -sac_ka * SAC_GUT;
dxdt_SAC_C   =  sac_ka * SAC_GUT * sac_kconv
               - (sac_CL_sac / sac_V_sac) * SAC_C;
dxdt_VAL_GUT = -sac_ka * VAL_GUT;
dxdt_VAL_C   =  sac_ka * VAL_GUT
               - (val_CL / val_V) * VAL_C;

// NEP inhibition and AT1R blockade effects
double NEP_inh = use_arni * sac_Cp / (sac_EC50_NEP + sac_Cp);
double AT1_blk = use_arni * val_Cp  / (val_EC50     + val_Cp);

//============================================================
// 3. FINERENONE PK (1-compartment)
//============================================================
double fin_Cp = FIN_C / fin_V;   // ng/mL

dxdt_FIN_GUT = -fin_ka * FIN_GUT;
dxdt_FIN_C   =  fin_ka * FIN_GUT - (fin_CL / fin_V) * FIN_C;

// MR blockade effect
double MR_blk = use_fin * fin_Emax * fin_Cp / (fin_EC50 + fin_Cp);

//============================================================
// 4. FUROSEMIDE PK (tubular secretion model)
//============================================================
double furo_Cp_plasma  = FURO_C / furo_V;           // μg/mL
double furo_Cp_tubular = furo_Cp_plasma * furo_Vt;  // tubular conc (proxy)

dxdt_FURO_GUT = -furo_ka  * FURO_GUT;
dxdt_FURO_C   =  furo_ka  * FURO_GUT
                - (furo_CL / furo_V) * FURO_C;

// Furosemide natriuresis effect
double furo_nat = use_furo * furo_Emax * furo_Cp_tubular /
                  (furo_EC50 + furo_Cp_tubular);

//============================================================
// 5. RAAS — AngII & Aldosterone
//============================================================
// AngII increases with HTN, decreases with AT1R blockade (neg feedback)
double AngII_current = base_AngII + ANGII;

// AT1R blockade reduces renin-independent AngII production
double AT1_eff   = 1.0 - AT1_blk;  // fraction of AT1R active
double kin_ANGII_adj = kin_AngII * AT1_eff * (1.0 + 0.3 * base_HTN);

dxdt_ANGII = kin_ANGII_adj - kout_AngII * AngII_current;

// Aldosterone driven by AngII; inhibited by MR blockade (feedback on synthesis)
double Aldo_stim = AngII_current * AngII_current /
                   (AngII_Aldo_EC50 * AngII_Aldo_EC50 +
                    AngII_current * AngII_current);
double Aldo_current = base_Aldo + ALDO;
double MR_suppress  = 1.0 - 0.5 * MR_blk;  // partial suppression via feedback

dxdt_ALDO = kin_AngII * Aldo_stim * MR_suppress - kout_AngII * Aldo_current;

//============================================================
// 6. NATRIURETIC PEPTIDES — ANP & BNP
//============================================================
double LVEDP_current = base_LVEDP + LVEDP_dyn;
double ANP_current   = base_ANP   + ANP;
double BNP_current   = base_BNP   + BNP;

// Pressure-dependent NP synthesis (Hill equation)
double LVEDP_stim = pow(LVEDP_current, BNP_hill) /
                    (pow(LVEDP_ANP_EC50, BNP_hill) +
                     pow(LVEDP_current,  BNP_hill));

// NEP inhibition increases ANP/BNP bioavailability (reduces kout)
double kout_ANP_eff = kout_ANP * (1.0 - 0.7 * NEP_inh);
double kout_BNP_eff = kout_BNP * (1.0 - 0.7 * NEP_inh);

dxdt_ANP = kin_ANP * (1.0 + 2.0 * LVEDP_stim) - kout_ANP_eff * ANP_current;
dxdt_BNP = kin_BNP * (1.0 + 3.0 * LVEDP_stim) - kout_BNP_eff * BNP_current;

//============================================================
// 7. cGMP-PKG PATHWAY
//============================================================
double cGMP_current = base_cGMP + cGMP_dyn;
double PKG_current  = base_PKG  + PKG_dyn;

// cGMP synthesis: driven by ANP/BNP (NPR-A) + NO (sGC)
// cGMP degradation: PDE5/9 (reduced by oxidative stress on PKG)
double cGMP_synth = kin_cGMP * (1.0 + 0.8 * (ANP_current - 1.0) +
                                 0.5 * (BNP_current - 1.0));
double PDE_activity = 1.0 + 0.3 * base_CKD_sev;  // PDE9 up in HFpEF comorbid
double cGMP_degrad  = kout_cGMP * PDE_activity * cGMP_current;

dxdt_cGMP_dyn = cGMP_synth - cGMP_degrad - base_cGMP * kout_cGMP;

// PKG activity: activated by cGMP, inhibited by oxidative stress (ROS/ONOO)
double ROS_proxy   = 1.0 + 0.4 * base_T2DM + 0.3 * base_HTN +
                     0.2 * base_CKD_sev;
double PKG_inhibit = ROS_proxy / (PKG_ox_EC50 + ROS_proxy);
double PKG_activation = cGMP_current / (1.0 + cGMP_current);

dxdt_PKG_dyn = 0.5 * PKG_activation * (1.0 - PKG_inhibit) -
               0.3 * PKG_current + 0.3 * base_PKG;  // equilibrate

//============================================================
// 8. CARDIAC STRUCTURE
//    LVM Index (g/m²), Fibrosis (fraction), Titin stiffness
//============================================================
double LVM_current    = base_LVM_idx   + LVM_IDX;
double fibr_current   = base_fibrosis  + FIBROSIS;
double titin_current  = base_titin_stiff + TITIN_STIFF;

// LV hypertrophy: driven by AngII × ET-1 signal (proxy = AngII_current)
// Regressed by ARNI/MRA reducing AngII/Aldo
double LVM_AngII_drive = kgrowth_LVM * AngII_current * (1.0 - 0.4 * AT1_blk);
double LVM_regress     = kdecay_LVM  * (LVM_current - 95.0);  // toward normal

dxdt_LVM_IDX = LVM_AngII_drive - LVM_regress;

// Fibrosis: driven by Aldo→MR and TGF-β; reduced by MR blockade + SGLT2i
double TGFb_proxy  = 1.0 + 0.3 * Aldo_current + 0.2 * AngII_current;
double fib_drive   = kfib_up   * TGFb_proxy * (1.0 - 0.5 * MR_blk)
                               * (1.0 - 0.3 * SGLT2_inh);
double fib_regress = kfib_down * fibr_current;

dxdt_FIBROSIS = fib_drive - fib_regress;

// Titin stiffness: increases with aging; reduced by PKG phosphorylation
double titin_age_effect = ktitin_age * (Age - 60.0) * 0.01;
double titin_PKG_soften = ktitin_PKG * PKG_current;

dxdt_TITIN_STIFF = titin_age_effect - titin_PKG_soften * titin_current
                   + ktitin_PKG * base_titin_stiff * 0.1;

//============================================================
// 9. HEMODYNAMICS
//    LVEDP & SVR as dynamic deviation from baseline
//============================================================
// LVEDP depends on titin stiffness, fibrosis, volume load
double titin_LVEDP = 8.0 * titin_current;         // stiffness contribution (mmHg)
double fibr_LVEDP  = 6.0 * (fibr_current / 0.30); // fibrosis contribution (mmHg)

// Volume effect: Aldosterone/RAAS increases preload; furosemide/SGLT2i reduces it
double vol_LVEDP   = 4.0 * (Aldo_current - 1.0) * (1.0 - furo_nat * 0.8)
                   * (1.0 - SGLT2_inh * 0.4);

double target_LVEDP = titin_LVEDP + fibr_LVEDP + vol_LVEDP + 4.0;
// Equilibrium rate toward target LVEDP
dxdt_LVEDP_dyn = 0.05 * (target_LVEDP - LVEDP_current);

// SVR: increases with AngII, ET-1, decreases with NP/cGMP and ARNI
double target_SVR = base_SVR * (1.0 + 0.15 * (AngII_current - 1.0))
                  * (1.0 - 0.20 * AT1_blk)
                  * (1.0 - 0.10 * (cGMP_current - 1.0));
double SVR_current = base_SVR + SVR_dyn;

dxdt_SVR_dyn = 0.03 * (target_SVR - SVR_current);

//============================================================
// 10. INFLAMMATION
//============================================================
double CRP_current = base_CRP + CRP_dyn;
double IL6_current = base_IL6 + IL6_dyn;

// IL-6 driven by obesity + T2DM + AngII; reduced by SGLT2i and MRA
double IL6_drive  = 0.1 * (1.0 + 0.3 * BMI / 30.0 + 0.2 * base_T2DM
                           + 0.2 * AngII_current)
                  * (1.0 - 0.25 * SGLT2_inh) * (1.0 - 0.15 * MR_blk);
double IL6_decay  = 0.05 * IL6_current;

dxdt_IL6_dyn = IL6_drive - IL6_decay;

// CRP synthesized by liver in response to IL-6
double CRP_synth  = kin_CRP * IL6_current;
double CRP_decay  = kout_CRP * CRP_current;

dxdt_CRP_dyn = CRP_synth - CRP_decay;

//============================================================
// 11. RENAL FUNCTION
//============================================================
double GFR_current   = base_GFR  + GFR_dyn;
double Na_current    = base_Na_exc + Na_exc_dyn;

// GFR: reduced by efferent dilation (SGLT2i initial hemodynamic effect)
//      then partially restored as sodium balance normalizes
double GFR_aldo_effect = -kGFR_aldo * (Aldo_current - 1.0);
double GFR_sglt2_init  = 0.05 * SGLT2_inh * (-1.0); // initial dip
double GFR_target      = base_GFR * (1.0 - 0.2 * base_CKD_sev)
                       + GFR_aldo_effect + GFR_sglt2_init * 0.5;

dxdt_GFR_dyn = 0.02 * (GFR_target - GFR_current);

// Na excretion: increased by furosemide, SGLT2i, ARNI (NP effect)
// decreased by Aldosterone
double Na_exc_target = base_Na_exc
                     + furo_nat   * 0.5
                     + SGLT2_inh  * 0.15
                     + NEP_inh    * 0.10
                     - (Aldo_current - 1.0) * 0.20;

dxdt_Na_exc_dyn = 0.1 * (Na_exc_target - Na_current);

//============================================================
// 12. NT-proBNP (clinical biomarker, pg/mL)
//    NT-proBNP ≈ 4-6x BNP; NEP does NOT degrade NT-proBNP
//    but rises with filling pressure; ARNI increases BNP not NT-proBNP
//============================================================
double NTpBNP_current = 400.0 + NT_proBNP_dyn;  // baseline 400 pg/mL (HFpEF)

double NTpBNP_target  = 400.0
  + 15.0 * (LVEDP_current - base_LVEDP)     // pressure-driven rise
  + 25.0 * (fibr_current  - base_fibrosis) * 1000.0
  - 80.0 * furo_nat                          // relief of congestion
  - 60.0 * SGLT2_inh;                        // SGLT2i effect (EMPEROR-Pres)

dxdt_NT_proBNP_dyn = 0.008 * (NTpBNP_target - NTpBNP_current);

$TABLE
// ── PK outputs ───────────────────────────────────────────
capture Empa_Cp     = EMPA_C / empa_V1;      // ng/mL
capture Sac_Cp      = SAC_C  / sac_V_sac;    // ng/mL (sacubitrilat)
capture Val_Cp      = VAL_C  / val_V;         // ng/mL
capture Fin_Cp      = FIN_C  / fin_V;         // ng/mL
capture Furo_Cp     = FURO_C / furo_V;        // μg/mL

// ── PD / effect ───────────────────────────────────────────
capture SGLT2_inhibition = use_empa * empa_Emax * Empa_Cp /
                           (empa_EC50 + Empa_Cp);
capture NEP_inhibition   = use_arni * Sac_Cp / (sac_EC50_NEP + Sac_Cp);
capture AT1R_blockade    = use_arni * Val_Cp  / (val_EC50     + Val_Cp);
capture MR_blockade      = use_fin  * fin_Emax * Fin_Cp / (fin_EC50 + Fin_Cp);
capture Furosemide_Natriuresis = use_furo * furo_Emax * (FURO_C / furo_V * furo_Vt) /
                                 (furo_EC50 + FURO_C / furo_V * furo_Vt);

// ── Neurohumoral ─────────────────────────────────────────
capture AngII_norm    = base_AngII + ANGII;
capture Aldosterone_n = base_Aldo  + ALDO;
capture ANP_level     = base_ANP   + ANP;
capture BNP_level     = base_BNP   + BNP;

// ── Cardiac / structural ─────────────────────────────────
capture LV_Mass_idx   = base_LVM_idx   + LVM_IDX;   // g/m²
capture LV_Fibrosis   = base_fibrosis  + FIBROSIS;  // fraction (0-1)
capture Titin_S       = base_titin_stiff + TITIN_STIFF;
capture cGMP_level    = base_cGMP + cGMP_dyn;
capture PKG_activity  = base_PKG  + PKG_dyn;

// ── Hemodynamics ─────────────────────────────────────────
capture LVEDP_mmHg    = base_LVEDP + LVEDP_dyn;
capture SVR_val       = base_SVR   + SVR_dyn;
capture CO_val        = (base_SV / 1000.0) * base_HR;  // rough (L/min)
// simplified MAP: CO × SVR / 80
capture MAP_mmHg      = CO_val * SVR_val / 80.0;

// ── Renal ────────────────────────────────────────────────
capture eGFR          = base_GFR + GFR_dyn;
capture Na_excretion  = base_Na_exc + Na_exc_dyn;

// ── Inflammation ─────────────────────────────────────────
capture hsCRP         = base_CRP + CRP_dyn;
capture IL6_level     = base_IL6 + IL6_dyn;

// ── Biomarker ─────────────────────────────────────────────
capture NT_proBNP_pgmL = 400.0 + NT_proBNP_dyn;

// ── Composite clinical risk score (0–100) ────────────────
// Higher = worse prognosis
capture HFH_risk_score = 20.0 * (LVEDP_mmHg / 20.0)
                       + 20.0 * (NT_proBNP_pgmL / 400.0)
                       + 15.0 * LV_Fibrosis
                       + 15.0 * (1.0 - eGFR / 60.0)
                       + 10.0 * (hsCRP / 3.0)
                       + 10.0 * Titin_S
                       + 10.0 * (LV_Mass_idx / 110.0);

$CAPTURE
Empa_Cp Sac_Cp Val_Cp Fin_Cp Furo_Cp
SGLT2_inhibition NEP_inhibition AT1R_blockade MR_blockade
Furosemide_Natriuresis
AngII_norm Aldosterone_n ANP_level BNP_level
LV_Mass_idx LV_Fibrosis Titin_S cGMP_level PKG_activity
LVEDP_mmHg SVR_val CO_val MAP_mmHg
eGFR Na_excretion hsCRP IL6_level
NT_proBNP_pgmL HFH_risk_score
'

# Compile the model
mod <- mread_cache("hfpef_qsp", temp = TRUE, code = hfpef_model_code, quiet = TRUE)

# ============================================================
# DOSING EVENTS HELPER FUNCTION
# ============================================================

make_dosing <- function(
    empa_dose = 10,    # mg, QD
    arni_dose = 97,    # mg sacubitril, BID
    fin_dose  = 20,    # mg, QD
    furo_dose = 40,    # mg, BID
    use_empa = TRUE, use_arni = FALSE, use_fin = FALSE, use_furo = FALSE,
    start_day = 0, n_weeks = 52) {

  n_hours <- n_weeks * 7 * 24
  ev <- ev()

  if (use_empa && empa_dose > 0) {
    ev <- ev + ev(amt = empa_dose * 1000,  # convert mg → μg → ng (ng/mL×L=μg×1000)
                  cmt = "EMPA_GUT", ii = 24, addl = ceiling(n_hours/24) - 1,
                  time = start_day * 24)
  }
  if (use_arni && arni_dose > 0) {
    ev <- ev + ev(amt = arni_dose * 1000,
                  cmt = "SAC_GUT", ii = 12, addl = ceiling(n_hours/12) - 1,
                  time = start_day * 24)
    ev <- ev + ev(amt = 103 * 1000,  # valsartan 103mg
                  cmt = "VAL_GUT", ii = 12, addl = ceiling(n_hours/12) - 1,
                  time = start_day * 24)
  }
  if (use_fin && fin_dose > 0) {
    ev <- ev + ev(amt = fin_dose * 1000,
                  cmt = "FIN_GUT", ii = 24, addl = ceiling(n_hours/24) - 1,
                  time = start_day * 24)
  }
  if (use_furo && furo_dose > 0) {
    ev <- ev + ev(amt = furo_dose * 1000,
                  cmt = "FURO_GUT", ii = 12, addl = ceiling(n_hours/12) - 1,
                  time = start_day * 24)
  }
  ev
}

# ============================================================
# SCENARIO DEFINITIONS (5 treatment scenarios)
# ============================================================

scenarios <- list(
  list(
    name = "1. Placebo (Standard of Care)",
    use_empa = 0, use_arni = 0, use_fin = 0, use_furo = 1,
    label  = "Placebo + Furosemide",
    color  = "#888888"
  ),
  list(
    name = "2. Empagliflozin 10 mg QD (EMPEROR-Preserved)",
    use_empa = 1, use_arni = 0, use_fin = 0, use_furo = 1,
    label  = "Empagliflozin 10 mg",
    color  = "#e63946"
  ),
  list(
    name = "3. Sacubitril/Valsartan 97/103 mg BID (PARAGON-HF)",
    use_empa = 0, use_arni = 1, use_fin = 0, use_furo = 1,
    label  = "Sacubitril/Valsartan",
    color  = "#2196F3"
  ),
  list(
    name = "4. Finerenone 20 mg QD (non-diabetic HFpEF)",
    use_empa = 0, use_arni = 0, use_fin = 1, use_furo = 1,
    label  = "Finerenone 20 mg",
    color  = "#4CAF50"
  ),
  list(
    name = "5. Combination (Empa + Sacubitril/Val)",
    use_empa = 1, use_arni = 1, use_fin = 0, use_furo = 1,
    label  = "Empa + ARNI (Combo)",
    color  = "#FF9800"
  )
)

# ============================================================
# SIMULATION FUNCTION
# ============================================================

run_scenario <- function(sc, n_weeks = 52) {
  dose_ev <- make_dosing(
    use_empa = as.logical(sc$use_empa),
    use_arni = as.logical(sc$use_arni),
    use_fin  = as.logical(sc$use_fin),
    use_furo = as.logical(sc$use_furo),
    n_weeks  = n_weeks
  )

  sim_times <- seq(0, n_weeks * 7 * 24, by = 12)  # every 12 h

  out <- mod %>%
    param(use_empa = sc$use_empa,
          use_arni = sc$use_arni,
          use_fin  = sc$use_fin,
          use_furo = sc$use_furo) %>%
    ev(dose_ev) %>%
    mrgsim(end = max(sim_times), delta = 12, obsonly = TRUE) %>%
    as_tibble() %>%
    mutate(
      day      = time / 24,
      week     = day / 7,
      scenario = sc$label,
      color    = sc$color
    )
  out
}

# Run all scenarios
message("Running HFpEF QSP simulations across 5 treatment scenarios...")
results <- map_dfr(scenarios, run_scenario)
message("Simulations complete.")

# ============================================================
# DOSE-RESPONSE ANALYSIS — Empagliflozin
# ============================================================

empa_doses <- c(0, 2.5, 5, 10, 25)

dr_results <- map_dfr(empa_doses, function(d) {
  ev_d <- ev(amt = d * 1000, cmt = "EMPA_GUT", ii = 24, addl = 364 - 1)

  mod %>%
    param(use_empa = as.integer(d > 0),
          use_arni = 0, use_fin = 0, use_furo = 1) %>%
    ev(ev_d) %>%
    mrgsim(end = 52 * 7 * 24, delta = 24, obsonly = TRUE) %>%
    as_tibble() %>%
    filter(abs(time - 52 * 7 * 24) < 1) %>%
    mutate(empa_dose_mg = d)
}) %>%
  select(empa_dose_mg, LVEDP_mmHg, NT_proBNP_pgmL, LV_Fibrosis,
         SGLT2_inhibition, eGFR, HFH_risk_score)

# ============================================================
# REPORTING: KEY SUMMARY AT WEEK 52
# ============================================================

summary_tbl <- results %>%
  filter(abs(week - 52) < 0.5) %>%
  group_by(scenario) %>%
  summarise(
    Week           = round(mean(week), 1),
    LVEDP_mmHg     = round(mean(LVEDP_mmHg),     1),
    NT_proBNP_pgmL = round(mean(NT_proBNP_pgmL), 0),
    LV_Fibrosis    = round(mean(LV_Fibrosis),     3),
    LV_Mass_idx    = round(mean(LV_Mass_idx),     1),
    eGFR           = round(mean(eGFR),            1),
    hsCRP          = round(mean(hsCRP),           2),
    HFH_risk_score = round(mean(HFH_risk_score),  1),
    .groups = "drop"
  )

cat("\n====================================================\n")
cat("   HFpEF QSP MODEL — 52-WEEK OUTCOMES SUMMARY\n")
cat("====================================================\n")
print(summary_tbl, width = 120)

cat("\n--- Dose-Response: Empagliflozin at 52 Weeks ---\n")
print(dr_results, digits = 3, width = 100)

# ============================================================
# PLOTS
# ============================================================

theme_qsp <- theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "#1a3a5c", color = NA),
    strip.text  = element_text(color = "white", face = "bold"),
    legend.position = "bottom",
    legend.title    = element_blank(),
    panel.grid.minor = element_blank()
  )

results_weekly <- results %>% filter(week %% 1 < 0.1)

# Plot 1: LVEDP over time
p1 <- ggplot(results_weekly, aes(week, LVEDP_mmHg, color = scenario)) +
  geom_hline(yintercept = 12, linetype = 2, color = "gray50") +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = setNames(
    sapply(scenarios, `[[`, "color"), sapply(scenarios, `[[`, "label"))) +
  annotate("text", x = 50, y = 11.5, label = "Normal LVEDP (<12 mmHg)",
           size = 3, color = "gray40") +
  labs(title = "LV End-Diastolic Pressure over 52 Weeks",
       x = "Week", y = "LVEDP (mmHg)") +
  theme_qsp

# Plot 2: NT-proBNP biomarker
p2 <- ggplot(results_weekly, aes(week, NT_proBNP_pgmL, color = scenario)) +
  geom_hline(yintercept = 125, linetype = 2, color = "gray50") +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = setNames(
    sapply(scenarios, `[[`, "color"), sapply(scenarios, `[[`, "label"))) +
  annotate("text", x = 48, y = 115, label = "HFpEF threshold 125 pg/mL",
           size = 3, color = "gray40") +
  labs(title = "NT-proBNP Biomarker Trajectory",
       x = "Week", y = "NT-proBNP (pg/mL)") +
  theme_qsp

# Plot 3: LV Fibrosis
p3 <- ggplot(results_weekly, aes(week, LV_Fibrosis * 100, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = setNames(
    sapply(scenarios, `[[`, "color"), sapply(scenarios, `[[`, "label"))) +
  labs(title = "LV Interstitial Fibrosis",
       x = "Week", y = "Fibrosis (%)") +
  theme_qsp

# Plot 4: HF Hospitalization Risk Score
p4 <- ggplot(results_weekly, aes(week, HFH_risk_score, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = setNames(
    sapply(scenarios, `[[`, "color"), sapply(scenarios, `[[`, "label"))) +
  labs(title = "Composite HF Hospitalization Risk Score",
       subtitle = "Lower = Better Prognosis",
       x = "Week", y = "Risk Score (0-100)") +
  theme_qsp

# Plot 5: PK profiles at steady state (Day 7)
pk_ss <- results %>% filter(week >= 1, week <= 2)

p5 <- ggplot(pk_ss %>% filter(Empa_Cp > 0 | Sac_Cp > 0),
             aes(week * 7 * 24 - 7 * 24, color = scenario)) +
  geom_line(aes(y = Empa_Cp), linetype = 1, linewidth = 1) +
  labs(title = "Empagliflozin Plasma Concentration (Steady-State, Day 7-14)",
       x = "Hours post first dose", y = "Concentration (ng/mL)") +
  theme_qsp

# Plot 6: Dose-Response at 52 weeks
p6 <- ggplot(dr_results, aes(empa_dose_mg, HFH_risk_score)) +
  geom_line(color = "#e63946", linewidth = 1.2) +
  geom_point(color = "#e63946", size = 3) +
  geom_text(aes(label = paste0(empa_dose_mg, " mg")),
            vjust = -1, size = 3.5) +
  labs(title = "Empagliflozin Dose-Response at 52 Weeks",
       x = "Empagliflozin Dose (mg)", y = "HFH Risk Score") +
  theme_qsp

# Arrange and save
if (requireNamespace("gridExtra", quietly = TRUE)) {
  library(gridExtra)
  pdf("hfpef_qsp_results.pdf", width = 14, height = 18)
  grid.arrange(p1, p2, p3, p4, p6, ncol = 2, nrow = 3)
  dev.off()
  message("Plots saved to hfpef_qsp_results.pdf")
} else {
  ggsave("hfpef_lvedp.png", p1, width = 10, height = 5, dpi = 150)
  ggsave("hfpef_ntprobnp.png", p2, width = 10, height = 5, dpi = 150)
  ggsave("hfpef_fibrosis.png", p3, width = 10, height = 5, dpi = 150)
  ggsave("hfpef_risk.png", p4, width = 10, height = 5, dpi = 150)
  ggsave("hfpef_doseresponse.png", p6, width = 8, height = 5, dpi = 150)
  message("Individual plots saved as PNG files.")
}

# ============================================================
# VIRTUAL PATIENT POPULATION (n=500)
# ============================================================

message("Running virtual patient population simulation (n=500)...")

set.seed(42)
n_pop <- 500

vp_params <- tibble(
  ID       = 1:n_pop,
  Age      = rnorm(n_pop, 72, 8) %>% pmax(50) %>% pmin(90),
  Weight   = rnorm(n_pop, 95, 18) %>% pmax(55) %>% pmin(160),
  BMI      = rnorm(n_pop, 33, 6)  %>% pmax(20) %>% pmin(55),
  base_CKD_sev = rbeta(n_pop, 1.5, 4),
  base_T2DM    = rbinom(n_pop, 1, 0.60),
  base_HTN     = rbinom(n_pop, 1, 0.85),
  base_AF      = rbinom(n_pop, 1, 0.35),
  base_LVEDP   = rnorm(n_pop, 16, 3) %>% pmax(10),
  base_GFR     = rnorm(n_pop, 62, 18) %>% pmax(20) %>% pmin(100),
  base_fibrosis = rbeta(n_pop, 3, 6)
)

# Simulate empagliflozin vs placebo arms
vp_empa <- vp_params %>%
  mutate(use_empa = 1L, use_arni = 0L, use_fin = 0L, use_furo = 1L,
         arm = "Empagliflozin")
vp_pbo  <- vp_params %>%
  mutate(use_empa = 0L, use_arni = 0L, use_fin = 0L, use_furo = 1L,
         arm = "Placebo")
vp_all  <- bind_rows(vp_empa, vp_pbo)

# Run population simulation
pop_results <- mod %>%
  idata_set(vp_all) %>%
  ev(ev(amt = 10000, cmt = "EMPA_GUT", ii = 24, addl = 365 - 1)) %>%
  mrgsim(end = 52 * 168, delta = 168, obsonly = TRUE) %>%  # every week
  as_tibble()

# Population summary at week 52
pop_w52 <- pop_results %>%
  filter(abs(time - 52 * 168) < 1) %>%
  left_join(vp_all %>% select(ID, arm), by = "ID") %>%
  group_by(arm) %>%
  summarise(
    n = n(),
    LVEDP_mean   = round(mean(LVEDP_mmHg),     1),
    LVEDP_sd     = round(sd(LVEDP_mmHg),       1),
    NTpBNP_mean  = round(mean(NT_proBNP_pgmL), 0),
    NTpBNP_sd    = round(sd(NT_proBNP_pgmL),   0),
    eGFR_mean    = round(mean(eGFR),            1),
    eGFR_sd      = round(sd(eGFR),             1),
    Risk_mean    = round(mean(HFH_risk_score),  1),
    Risk_sd      = round(sd(HFH_risk_score),    1),
    .groups = "drop"
  )

cat("\n====================================================\n")
cat("   VIRTUAL POPULATION (n=500) — WEEK 52 OUTCOMES\n")
cat("====================================================\n")
print(pop_w52, width = 120)

message("\nHFpEF QSP simulation complete.")
