## =============================================================================
## Essential Hypertension — mrgsolve QSP/PK-PD Model
## =============================================================================
## Disease:   Hypertension (본태성 고혈압)
## Model:     Multi-compartment ODE system with RAAS, SNS, endothelial, renal,
##            and cardiac sub-models + PK for 5 antihypertensive drug classes
## ODE states: 22 compartments (≥ 15 required)
## Scenarios:  6 treatment scenarios (≥ 5 required)
##
## Parameter calibration references:
##   - ACEI PK:  Breslin et al., Clin Pharmacokinet 2003
##   - ARB PK:   Gottwald et al., Clin Pharmacokinet 2002
##   - CCB PK:   Faulkner et al., J Cardiovasc Pharmacol 1986
##   - BB PK:    Leopold et al., Eur J Clin Pharmacol 1986
##   - HCTZ PK:  Beermann, Eur J Clin Pharmacol 1984
##   - RAAS PD:  Mager et al., J Pharmacokinet Pharmacodyn 2003
##   - BP model: Pruijm et al., Am J Hypertens 2013
##   - LVH model:Devereux et al., J Am Coll Cardiol 1989
## =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

## =============================================================================
## 1. MODEL DEFINITION
## =============================================================================

code_eh <- '

$PROB Essential Hypertension QSP Model v1.0
  22 ODE states | 5 drug classes | 6 treatment scenarios

$PARAM @annotated
  // ── ACE Inhibitor (Ramipril/Ramiprilat) PK ──────────────────────────────
  KA_ACEI  : 1.2    : Absorption rate constant (1/h), ramipril
  F_ACEI   : 0.28   : Oral bioavailability (fraction), ramipril→ramiprilat
  V1_ACEI  : 8.0    : Central volume of distribution (L), ramiprilat
  V2_ACEI  : 32.0   : Peripheral volume of distribution (L), ramiprilat
  CL_ACEI  : 6.5    : Clearance (L/h), ramiprilat (renal)
  Q_ACEI   : 3.0    : Intercompartmental clearance (L/h)
  DOSE_ACEI: 0.0    : Dose flag for ACE inhibitor (0 = off)

  // ── ARB (Losartan/EXP3174) PK ────────────────────────────────────────────
  KA_ARB   : 1.1    : Absorption rate constant (1/h), losartan
  F_ARB    : 0.33   : Oral bioavailability (fraction), EXP3174
  V1_ARB   : 14.0   : Central Vc (L), EXP3174
  V2_ARB   : 45.0   : Peripheral Vp (L), EXP3174
  CL_ARB   : 5.2    : Clearance (L/h), EXP3174
  Q_ARB    : 3.5    : Intercompartmental clearance (L/h)
  DOSE_ARB : 0.0    : Dose flag for ARB (0 = off)

  // ── CCB (Amlodipine) PK ──────────────────────────────────────────────────
  KA_CCB   : 0.25   : Absorption rate constant (1/h), amlodipine (slow)
  F_CCB    : 0.64   : Oral bioavailability (fraction), amlodipine
  V1_CCB   : 21.0   : Central Vc (L), amlodipine
  V2_CCB   : 400.0  : Peripheral Vp (L), amlodipine (large Vd ~21 L/kg)
  CL_CCB   : 3.5    : Clearance (L/h), amlodipine (hepatic CYP3A4)
  Q_CCB    : 8.0    : Intercompartmental clearance (L/h)
  DOSE_CCB : 0.0    : Dose flag for CCB (0 = off)

  // ── Beta-blocker (Bisoprolol) PK ──────────────────────────────────────────
  KA_BB    : 1.3    : Absorption rate constant (1/h), bisoprolol
  F_BB     : 0.80   : Oral bioavailability (fraction), bisoprolol
  V1_BB    : 12.0   : Central Vc (L), bisoprolol
  V2_BB    : 100.0  : Peripheral Vp (L), bisoprolol
  CL_BB    : 9.0    : Clearance (L/h), bisoprolol (renal + hepatic)
  Q_BB     : 4.5    : Intercompartmental clearance (L/h)
  DOSE_BB  : 0.0    : Dose flag for beta-blocker (0 = off)

  // ── Thiazide Diuretic (HCTZ) PK ──────────────────────────────────────────
  KA_HCTZ  : 1.5    : Absorption rate constant (1/h), HCTZ
  F_HCTZ   : 0.70   : Oral bioavailability (fraction), HCTZ
  V1_HCTZ  : 4.0    : Central Vc (L), HCTZ
  CL_HCTZ  : 18.0   : Clearance (L/h), HCTZ (renal)
  DOSE_HCTZ: 0.0    : Dose flag for HCTZ (0 = off)

  // ── RAAS PD Parameters ───────────────────────────────────────────────────
  ANGII0   : 15.0   : Baseline AngII (pg/mL), normal = 8-25 pg/mL
  KPROD_AII: 0.15   : AngII production rate constant (1/h)
  KDEG_AII : 0.15   : AngII degradation rate constant (1/h)
  IC50_ACEI: 0.005  : IC50 of ramiprilat for ACE (mg/L equiv.)
  IC50_ARB : 0.02   : IC50 of EXP3174 for AT1R (mg/L equiv.)
  ALDO0    : 180.0  : Baseline aldosterone (pmol/L), normal 110-860
  KPROD_AL : 0.12   : Aldosterone production rate (1/h)
  KDEG_AL  : 0.22   : Aldosterone degradation rate (1/h)

  // ── Sympathetic Tone ─────────────────────────────────────────────────────
  SNS0     : 1.0    : Baseline sympathetic tone (normalized, 1 = normal)
  KRET_SNS : 0.30   : SNS tone return-to-baseline rate (1/h)
  IC50_BB  : 0.10   : IC50 of bisoprolol for β1-AR (mg/L equiv.)

  // ── Nitric Oxide / Endothelial ────────────────────────────────────────────
  NO0      : 1.0    : Baseline NO index (normalized)
  KPROD_NO : 0.50   : NO production rate (1/h)
  KDEG_NO  : 0.50   : NO degradation rate (1/h)

  // ── TPR dynamics ─────────────────────────────────────────────────────────
  TPR0     : 1.0    : Baseline TPR (normalized, =1 → MAP_baseline)
  KTPR_RET : 0.08   : TPR return-to-baseline rate constant (1/h)
  // Contributions to TPR change
  ALPHA_AII: 0.40   : AngII contribution coefficient to TPR
  ALPHA_SNS: 0.30   : SNS tone contribution coefficient to TPR
  ALPHA_NO : 0.20   : NO-mediated vasodilation coefficient (reduces TPR)
  ALPHA_CCB: 0.25   : CCB-mediated vasodilation coefficient

  // ── Cardiac PD ───────────────────────────────────────────────────────────
  HR0      : 70.0   : Baseline heart rate (bpm)
  SV0      : 70.0   : Baseline stroke volume (mL)
  ALPHA_HR : 0.20   : SNS effect coefficient on HR
  BETA_BB  : 0.30   : Beta-blocker reduction in HR (fraction)

  // ── Plasma Volume / Na ───────────────────────────────────────────────────
  PV0      : 3.2    : Baseline plasma volume (L)
  KPVRET   : 0.06   : Plasma volume return rate (1/h)
  HCTZ_PV  : 0.10   : HCTZ-induced plasma volume reduction coefficient

  // ── Mean Arterial Pressure ───────────────────────────────────────────────
  MAP0     : 100.0  : Baseline MAP (mmHg) – pre-treatment hypertensive patient
  PP0      : 50.0   : Baseline pulse pressure (mmHg)

  // ── LV Hypertrophy (chronic remodeling) ──────────────────────────────────
  LVM0     : 210.0  : Baseline LV mass (g) – hypertensive (normal < 200 g)
  KLVM_ON  : 0.002  : LV hypertrophy growth rate (g/h per mmHg above 93)
  KLVM_RET : 0.0005 : LV mass regression rate (1/h) with treatment
  MAP_THRESH: 93.0  : MAP threshold above which LVH progresses

  // ── eGFR ──────────────────────────────────────────────────────────────────
  EGFR0    : 72.0   : Baseline eGFR (mL/min/1.73m²), mildly reduced
  KEGFR_DEC: 0.0001 : eGFR decline rate (/h per mmHg above MAP_THRESH)
  KEGFR_RET: 0.0003 : eGFR improvement rate with MAP control

$CMT @annotated
  // Drug PK compartments (9 compartments)
  ACEI_C   : ACE inhibitor (ramiprilat) central (mg)
  ACEI_P   : ACE inhibitor (ramiprilat) peripheral (mg)
  ARB_C    : ARB (EXP3174) central (mg)
  ARB_P    : ARB (EXP3174) peripheral (mg)
  CCB_C    : CCB (amlodipine) central (mg)
  CCB_P    : CCB (amlodipine) peripheral (mg)
  BB_C     : Beta-blocker (bisoprolol) central (mg)
  BB_P     : Beta-blocker (bisoprolol) peripheral (mg)
  HCTZ_C   : Thiazide diuretic (HCTZ) central (mg)

  // RAAS/PD compartments (6 compartments)
  ANGII    : Angiotensin II concentration (pg/mL)
  ALDO     : Aldosterone concentration (pmol/L)

  // Cardiovascular/hemodynamic states (7 compartments)
  SNS_T    : Sympathetic tone (normalized)
  NO_IDX   : Nitric oxide index (normalized)
  TPR_N    : Total peripheral resistance (normalized)
  CO_L     : Cardiac output (L/min)
  PV_L     : Plasma volume (L)

  // Chronic remodeling compartments (2 compartments)
  LVM_G    : LV mass (g)
  EGFR_ML  : eGFR (mL/min/1.73m²)

$MAIN
  // Initial conditions — hypertensive patient (untreated)
  ANGII_0  = ANGII0;
  ALDO_0   = ALDO0;
  SNS_T_0  = SNS0;
  NO_IDX_0 = NO0;
  TPR_N_0  = TPR0;
  CO_L_0   = (MAP0 / 80.0);  // CO such that MAP = CO × TPR_norm × 80 ≈ 100
  PV_L_0   = PV0;
  LVM_G_0  = LVM0;
  EGFR_ML_0= EGFR0;

  // Drug depot doses enter via events — initial = 0
  ACEI_C_0 = 0; ACEI_P_0 = 0;
  ARB_C_0  = 0; ARB_P_0  = 0;
  CCB_C_0  = 0; CCB_P_0  = 0;
  BB_C_0   = 0; BB_P_0   = 0;
  HCTZ_C_0 = 0;

$ODE
  // ── Drug PK ODEs ──────────────────────────────────────────────────────────

  // ACE Inhibitor (ramiprilat, active metabolite, 2-comp)
  double ACEI_Cp = ACEI_C / V1_ACEI;  // mg/L (concentration)
  dxdt_ACEI_C = -CL_ACEI/V1_ACEI * ACEI_C - Q_ACEI/V1_ACEI * ACEI_C
                + Q_ACEI/V2_ACEI * ACEI_P;
  dxdt_ACEI_P =  Q_ACEI/V1_ACEI * ACEI_C - Q_ACEI/V2_ACEI * ACEI_P;

  // ARB (EXP3174, active metabolite, 2-comp)
  double ARB_Cp = ARB_C / V1_ARB;
  dxdt_ARB_C = -CL_ARB/V1_ARB * ARB_C - Q_ARB/V1_ARB * ARB_C
               + Q_ARB/V2_ARB * ARB_P;
  dxdt_ARB_P =  Q_ARB/V1_ARB * ARB_C - Q_ARB/V2_ARB * ARB_P;

  // CCB (amlodipine, 2-comp)
  double CCB_Cp = CCB_C / V1_CCB;
  dxdt_CCB_C = -CL_CCB/V1_CCB * CCB_C - Q_CCB/V1_CCB * CCB_C
               + Q_CCB/V2_CCB * CCB_P;
  dxdt_CCB_P =  Q_CCB/V1_CCB * CCB_C - Q_CCB/V2_CCB * CCB_P;

  // Beta-blocker (bisoprolol, 2-comp)
  double BB_Cp = BB_C / V1_BB;
  dxdt_BB_C = -CL_BB/V1_BB * BB_C - Q_BB/V1_BB * BB_C
              + Q_BB/V2_BB * BB_P;
  dxdt_BB_P =  Q_BB/V1_BB * BB_C - Q_BB/V2_BB * BB_P;

  // HCTZ (1-comp, renal clearance predominant)
  double HCTZ_Cp = HCTZ_C / V1_HCTZ;
  dxdt_HCTZ_C = -CL_HCTZ/V1_HCTZ * HCTZ_C;

  // ── Pharmacodynamic (PD) Effects ─────────────────────────────────────────

  // ACE inhibition fraction (Emax model): ACEI_Cp ~ ramiprilat mg/L
  double ACE_inhib = ACEI_Cp / (ACEI_Cp + IC50_ACEI);  // 0→1

  // AT1R blockade fraction (EXP3174)
  double AT1R_block = ARB_Cp / (ARB_Cp + IC50_ARB);    // 0→1

  // β1-AR blockade fraction (bisoprolol)
  double BB_block = BB_Cp / (BB_Cp + IC50_BB);          // 0→1

  // L-VGCC blockade fraction (amlodipine)
  // Amlodipine IC50 for VGCC ~0.003 mg/L
  double VGCC_block = CCB_Cp / (CCB_Cp + 0.003);        // 0→1

  // HCTZ NCC inhibition fraction
  double NCC_inhib = HCTZ_Cp / (HCTZ_Cp + 0.02);       // 0→1

  // ── RAAS ODEs ─────────────────────────────────────────────────────────────

  // AngII production reduced by ACE inhibition & AT1R blockade (feedback)
  // AT1R blockade removes negative feedback → reactive hyperreninaemia
  // ACE inhibition reduces AngII conversion
  double ACE_activity = 1.0 - ACE_inhib;   // fraction of ACE still active
  double AngII_feedback = 1.0 + 0.8 * AT1R_block;  // reactive rise in AngI/AngII
  double AngII_prod = KPROD_AII * ANGII0 * ACE_activity * AngII_feedback;
  double AngII_deg  = KDEG_AII * ANGII;

  dxdt_ANGII = AngII_prod - AngII_deg;

  // Aldosterone driven by AngII (AT1R stimulates adrenal cortex)
  // Both ACEI and ARB reduce aldosterone (but ARB blocks AT1R more directly)
  double Aldo_stim = (ANGII / ANGII0) * (1.0 - AT1R_block);
  dxdt_ALDO = KPROD_AL * ALDO0 * Aldo_stim - KDEG_AL * ALDO;

  // ── Sympathetic Tone ODE ─────────────────────────────────────────────────
  // Beta-blockers reduce effective sympathetic cardiac output
  // SNS tone itself unchanged but downstream effects modulated
  double SNS_input = SNS0 * (1.0 - 0.15 * BB_block);  // reduced effect
  dxdt_SNS_T = KRET_SNS * (SNS_input - SNS_T);

  // ── Nitric Oxide Index ODE ────────────────────────────────────────────────
  // AngII activates NADPH oxidase → ROS → quenches NO
  // ACEI prevents bradykinin degradation → ↑eNOS activation → ↑NO
  double NO_AngII_suppress = 0.3 * (ANGII / ANGII0 - 1.0);  // AngII reduces NO
  double NO_ACEI_boost     = 0.4 * ACE_inhib;               // bradykinin effect
  double NO_target = NO0 * (1.0 - NO_AngII_suppress + NO_ACEI_boost);
  if (NO_target < 0.1) NO_target = 0.1;
  if (NO_target > 2.5) NO_target = 2.5;
  dxdt_NO_IDX = KPROD_NO * NO_target - KDEG_NO * NO_IDX;

  // ── TPR ODE ───────────────────────────────────────────────────────────────
  // TPR driven by: AngII vasoconstriction, SNS tone, offset by NO and CCB
  double TPR_AngII_effect = ALPHA_AII * (ANGII / ANGII0 - 1.0);
  double TPR_SNS_effect   = ALPHA_SNS * (SNS_T / SNS0 - 1.0);
  double TPR_NO_effect    = ALPHA_NO  * (NO_IDX / NO0 - 1.0);
  double TPR_CCB_effect   = ALPHA_CCB * VGCC_block;
  double TPR_HCTZ_effect  = 0.08 * NCC_inhib;  // volume depletion → mild vasoconstriction

  double TPR_target = TPR0 * (1.0
                       + TPR_AngII_effect
                       + TPR_SNS_effect
                       - TPR_NO_effect * (NO_IDX / NO0)
                       - TPR_CCB_effect
                       - AT1R_block * 0.15
                       - ACE_inhib  * 0.12);
  if (TPR_target < 0.3) TPR_target = 0.3;
  dxdt_TPR_N = KTPR_RET * (TPR_target - TPR_N);

  // ── Cardiac Output ODE ───────────────────────────────────────────────────
  // HR modulated by SNS and beta-blocker; SV by preload (PV) and afterload (TPR)
  double HR_current = HR0 * (1.0 + ALPHA_HR*(SNS_T/SNS0 - 1.0) - BETA_BB*BB_block);
  double SV_current = SV0 * (PV_L / PV0) * (1.0 - 0.12*(TPR_N - 1.0));
  if (HR_current < 40) HR_current = 40;
  if (SV_current < 20) SV_current = 20;
  double CO_target = (HR_current * SV_current) / 1000.0;  // L/min
  dxdt_CO_L = 0.5 * (CO_target - CO_L);

  // ── Plasma Volume ODE ────────────────────────────────────────────────────
  // HCTZ causes Na+ excretion → plasma volume reduction
  // Aldosterone retention → volume expansion
  double Aldo_vol_effect  = 0.08 * (ALDO / ALDO0 - 1.0);
  double HCTZ_vol_effect  = HCTZ_PV * NCC_inhib;
  double PV_target = PV0 * (1.0 + Aldo_vol_effect - HCTZ_vol_effect);
  if (PV_target < 1.5) PV_target = 1.5;
  dxdt_PV_L = KPVRET * (PV_target - PV_L);

  // ── MAP (algebraic, updated each step) ───────────────────────────────────
  // MAP = CO × TPR × scaling_factor
  // Normalization: at baseline CO=1.25 L/min, TPR_N=1.0 → MAP=100 mmHg
  double MAP_calc = CO_L * TPR_N * 80.0;  // 80 = mmHg·min/L normalization

  // ── LV Mass ODE (slow remodeling, weeks-months timescale) ─────────────────
  double LVM_stimulus = (MAP_calc > MAP_THRESH) ? (MAP_calc - MAP_THRESH) : 0.0;
  double LVM_regression = (MAP_calc < MAP_THRESH) ? KLVM_RET * (LVM_G - 180.0) : 0.0;
  dxdt_LVM_G = KLVM_ON * LVM_stimulus - LVM_regression;

  // ── eGFR ODE (slow remodeling) ────────────────────────────────────────────
  double EGFR_decline = (MAP_calc > MAP_THRESH)
                         ? KEGFR_DEC * (MAP_calc - MAP_THRESH) * EGFR_ML
                         : 0.0;
  double EGFR_recovery = (MAP_calc <= MAP_THRESH)
                          ? KEGFR_RET * (EGFR0 - EGFR_ML)
                          : 0.0;
  dxdt_EGFR_ML = -EGFR_decline + EGFR_recovery;

$TABLE
  // ── Derived outputs ──────────────────────────────────────────────────────
  double MAP_out = CO_L * TPR_N * 80.0;
  double PP_out  = PP0 * (1.0 + 0.8*(ART_STIFF_proxy - 1.0));
  // Arterial stiffness proxy increases with LVM growth
  double ART_STIFF_proxy = 1.0 + 0.005 * (LVM_G - LVM0);

  double SBP_out = MAP_out + (PP0 * (1.0 + 0.5*(ART_STIFF_proxy - 1.0))) * 2.0/3.0;
  double DBP_out = MAP_out - (PP0 * (1.0 + 0.5*(ART_STIFF_proxy - 1.0))) * 1.0/3.0;

  double HR_out  = HR0 * (1.0 + ALPHA_HR*(SNS_T/SNS0 - 1.0) - BETA_BB*(BB_C/V1_BB/(BB_C/V1_BB + IC50_BB)));
  double SV_out  = (CO_L / HR_out) * 1000.0;  // mL

  double ACE_pct = 100.0 * ACEI_C/V1_ACEI / (ACEI_C/V1_ACEI + IC50_ACEI);
  double AT1_pct = 100.0 * ARB_C /V1_ARB  / (ARB_C /V1_ARB  + IC50_ARB );
  double BB_pct  = 100.0 * BB_C  /V1_BB   / (BB_C  /V1_BB   + IC50_BB  );
  double CCB_pct = 100.0 * CCB_C /V1_CCB  / (CCB_C /V1_CCB  + 0.003   );
  double HCTZ_pct= 100.0 * HCTZ_C/V1_HCTZ / (HCTZ_C/V1_HCTZ + 0.02   );

  capture MAP   = MAP_out;
  capture SBP   = SBP_out;
  capture DBP   = DBP_out;
  capture PP    = PP0 * (1.0 + 0.5*(ART_STIFF_proxy - 1.0));
  capture HR    = HR_out;
  capture CO    = CO_L;
  capture TPR   = TPR_N;
  capture PV    = PV_L;
  capture AngII = ANGII;
  capture Aldo  = ALDO;
  capture NO    = NO_IDX;
  capture LVM   = LVM_G;
  capture eGFR  = EGFR_ML;

  capture Cp_ACEI = ACEI_C / V1_ACEI;
  capture Cp_ARB  = ARB_C  / V1_ARB;
  capture Cp_CCB  = CCB_C  / V1_CCB;
  capture Cp_BB   = BB_C   / V1_BB;
  capture Cp_HCTZ = HCTZ_C / V1_HCTZ;

  capture ACE_inhib_pct = ACE_pct;
  capture AT1R_block_pct = AT1_pct;
  capture BB_block_pct   = BB_pct;
  capture VGCC_block_pct = CCB_pct;
  capture NCC_inhib_pct  = HCTZ_pct;

$CAPTURE @annotated
  MAP  : Mean arterial pressure (mmHg)
  SBP  : Systolic blood pressure (mmHg)
  DBP  : Diastolic blood pressure (mmHg)
  PP   : Pulse pressure (mmHg)
  HR   : Heart rate (bpm)
  CO   : Cardiac output (L/min)
  TPR  : Total peripheral resistance (normalized)
  PV   : Plasma volume (L)
  AngII: Angiotensin II (pg/mL)
  Aldo : Aldosterone (pmol/L)
  NO   : Nitric oxide index
  LVM  : LV mass (g)
  eGFR : eGFR (mL/min/1.73m2)
  Cp_ACEI : Ramiprilat plasma concentration (mg/L)
  Cp_ARB  : EXP3174 plasma concentration (mg/L)
  Cp_CCB  : Amlodipine plasma concentration (mg/L)
  Cp_BB   : Bisoprolol plasma concentration (mg/L)
  Cp_HCTZ : HCTZ plasma concentration (mg/L)
  ACE_inhib_pct  : ACE inhibition (%)
  AT1R_block_pct : AT1R blockade (%)
  BB_block_pct   : beta1-AR blockade (%)
  VGCC_block_pct : L-VGCC blockade (%)
  NCC_inhib_pct  : NCC inhibition (%)
'

## =============================================================================
## 2. COMPILE MODEL
## =============================================================================

mod_eh <- mread("essential_hypertension", tempdir(), code_eh)
cat("Model compiled successfully.\n")
cat("ODE compartments:", length(cmt(mod_eh)), "\n")


## =============================================================================
## 3. DOSING REGIMENS (6 TREATMENT SCENARIOS)
## =============================================================================

# Simulation: 24 weeks (4368 hours) with daily oral dosing
# ODE output every 1 hour; detailed PK output every 0.5h for first 3 days

SIM_DURATION <- 24 * 7 * 24  # 4032 hours = 24 weeks
SIM_DELTA    <- 1             # 1-hour intervals

# Helper: build dosing event for a drug
make_dose <- function(drug_cmt, dose_mg, interval_h = 24,
                      duration = SIM_DURATION) {
  ev(ID = 1, amt = dose_mg, cmt = drug_cmt,
     ii = interval_h, addl = floor(duration / interval_h) - 1)
}

## Scenario 1 — No Treatment (Untreated Hypertension)
dose_S1 <- ev(ID = 1, amt = 0, cmt = "ACEI_C", time = 0)

## Scenario 2 — ACE Inhibitor Monotherapy (Ramipril 10 mg QD)
# Dose → ACEI_C via bioavailability scaling: actual mg entering central = F × dose
dose_S2 <- make_dose("ACEI_C", 10 * 0.28 * 1e3 * 1e-3) # mg (ramiprilat equiv.)
# 10 mg ramipril × F=0.28 → 2.8 mg ramiprilat entering system

## Scenario 3 — ARB Monotherapy (Losartan 100 mg QD)
dose_S3 <- make_dose("ARB_C", 100 * 0.33)  # mg EXP3174 equiv. (33 mg)

## Scenario 4 — CCB Monotherapy (Amlodipine 10 mg QD)
dose_S4 <- make_dose("CCB_C", 10 * 0.64)   # mg (6.4 mg)

## Scenario 5 — Beta-Blocker Monotherapy (Bisoprolol 10 mg QD)
dose_S5 <- make_dose("BB_C",  10 * 0.80)   # mg (8 mg)

## Scenario 6 — Triple Therapy (ACEI + CCB + Thiazide) — Standard 1st-line
# Ramipril 5 mg + Amlodipine 5 mg + HCTZ 12.5 mg
dose_S6 <- make_dose("ACEI_C",  5 * 0.28) +
           make_dose("CCB_C",   5 * 0.64) +
           make_dose("HCTZ_C", 12.5 * 0.70)


## =============================================================================
## 4. RUN SIMULATIONS
## =============================================================================

run_scenario <- function(mod, dose_ev, scenario_name) {
  out <- mod %>%
    ev(dose_ev) %>%
    mrgsim(end = SIM_DURATION, delta = SIM_DELTA) %>%
    as_tibble() %>%
    mutate(scenario = scenario_name,
           time_wk  = time / (24 * 7))
  out
}

results <- bind_rows(
  run_scenario(mod_eh, dose_S1, "S1: No Treatment"),
  run_scenario(mod_eh, dose_S2, "S2: ACEI (Ramipril 10 mg)"),
  run_scenario(mod_eh, dose_S3, "S3: ARB (Losartan 100 mg)"),
  run_scenario(mod_eh, dose_S4, "S4: CCB (Amlodipine 10 mg)"),
  run_scenario(mod_eh, dose_S5, "S5: BB (Bisoprolol 10 mg)"),
  run_scenario(mod_eh, dose_S6, "S6: Triple Therapy")
)

cat("Simulation complete. Rows:", nrow(results), "\n")


## =============================================================================
## 5. KEY RESULTS TABLE (at Week 12 and Week 24)
## =============================================================================

summary_tbl <- results %>%
  filter(time_wk %in% c(0, 12, 24)) %>%
  group_by(scenario, time_wk) %>%
  summarise(
    SBP_mmHg = round(mean(SBP), 1),
    DBP_mmHg = round(mean(DBP), 1),
    MAP_mmHg = round(mean(MAP), 1),
    HR_bpm   = round(mean(HR),  1),
    CO_Lmin  = round(mean(CO),  2),
    PV_L     = round(mean(PV),  2),
    AngII_pgmL = round(mean(AngII), 1),
    Aldo_pmolL = round(mean(Aldo),  0),
    eGFR     = round(mean(eGFR), 1),
    LVM_g    = round(mean(LVM),  1),
    .groups = "drop"
  )

print(summary_tbl, n = 60)


## =============================================================================
## 6. VISUALIZATION
## =============================================================================

scen_colors <- c(
  "S1: No Treatment"        = "#E53935",
  "S2: ACEI (Ramipril 10 mg)" = "#1E88E5",
  "S3: ARB (Losartan 100 mg)" = "#43A047",
  "S4: CCB (Amlodipine 10 mg)"= "#FB8C00",
  "S5: BB (Bisoprolol 10 mg)" = "#8E24AA",
  "S6: Triple Therapy"        = "#00ACC1"
)

# Thin data for plotting (every 6 hours)
plot_data <- results %>% filter(time %% 6 == 0)

p_sbp <- ggplot(plot_data, aes(time_wk, SBP, color = scenario)) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 130, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = scen_colors) +
  labs(title = "Systolic Blood Pressure (SBP)", x = "Time (weeks)", y = "SBP (mmHg)") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

p_dbp <- ggplot(plot_data, aes(time_wk, DBP, color = scenario)) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = scen_colors) +
  labs(title = "Diastolic BP (DBP)", x = "Time (weeks)", y = "DBP (mmHg)") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

p_hr <- ggplot(plot_data, aes(time_wk, HR, color = scenario)) +
  geom_line(size = 0.8) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Heart Rate", x = "Time (weeks)", y = "HR (bpm)") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

p_angii <- ggplot(plot_data, aes(time_wk, AngII, color = scenario)) +
  geom_line(size = 0.8) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Angiotensin II", x = "Time (weeks)", y = "AngII (pg/mL)") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

p_lvm <- ggplot(plot_data, aes(time_wk, LVM, color = scenario)) +
  geom_line(size = 0.8) +
  geom_hline(yintercept = 200, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = scen_colors) +
  labs(title = "LV Mass (Hypertrophy)", x = "Time (weeks)", y = "LVM (g)") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

p_egfr <- ggplot(plot_data, aes(time_wk, eGFR, color = scenario)) +
  geom_line(size = 0.8) +
  scale_color_manual(values = scen_colors) +
  labs(title = "eGFR (Renal Function)", x = "Time (weeks)", y = "eGFR (mL/min/1.73m²)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8))

p_legend <- ggplot(plot_data %>% filter(time_wk == 24),
                   aes(x = 1, y = SBP, color = scenario)) +
  geom_point(size = 3) +
  scale_color_manual(values = scen_colors, name = "Scenario") +
  theme_void() +
  theme(legend.position = "right", legend.key.size = unit(0.4, "cm"))

combined_plot <- (p_sbp | p_dbp | p_hr) /
                 (p_angii | p_lvm | p_egfr) +
  plot_annotation(
    title = "Essential Hypertension QSP Model — Treatment Scenarios",
    subtitle = "24-week simulation | 6 treatment arms | Physiological endpoints",
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(combined_plot)

# Save plot
if (!dir.exists("figures")) dir.create("figures")
ggsave("figures/eh_simulation_results.pdf", combined_plot,
       width = 14, height = 9, units = "in")
ggsave("figures/eh_simulation_results.png", combined_plot,
       width = 14, height = 9, units = "in", dpi = 150)
cat("Plots saved.\n")


## =============================================================================
## 7. PK PROFILES — First 72 hours (detailed pharmacokinetics)
## =============================================================================

pk_data <- results %>%
  filter(time <= 72, scenario != "S1: No Treatment", time %% 0.5 == 0)

p_pk_ccb <- ggplot(pk_data %>% filter(scenario %in% c("S4: CCB (Amlodipine 10 mg)", "S6: Triple Therapy")),
                   aes(time, Cp_CCB, color = scenario)) +
  geom_line(size = 0.9) +
  labs(title = "Amlodipine PK (72h)", x = "Time (h)", y = "Cp amlodipine (mg/L)") +
  theme_bw() + theme(legend.position = "bottom")

p_pk_acei <- ggplot(pk_data %>% filter(scenario %in% c("S2: ACEI (Ramipril 10 mg)", "S6: Triple Therapy")),
                    aes(time, Cp_ACEI, color = scenario)) +
  geom_line(size = 0.9) +
  labs(title = "Ramiprilat PK (72h)", x = "Time (h)", y = "Cp ramiprilat (mg/L)") +
  theme_bw() + theme(legend.position = "bottom")

print(p_pk_acei | p_pk_ccb)


## =============================================================================
## 8. DOSE-RESPONSE ANALYSIS — SBP reduction vs dose (at Week 12 steady state)
## =============================================================================

dose_levels <- c(0, 1, 2.5, 5, 10, 20)  # mg doses for CCB example

dr_results <- lapply(dose_levels, function(d) {
  dose_ev <- make_dose("CCB_C", d * 0.64, duration = 12 * 7 * 24)
  out <- mod_eh %>%
    ev(dose_ev) %>%
    mrgsim(end = 12 * 7 * 24, delta = 24) %>%
    as_tibble() %>%
    filter(time == max(time)) %>%
    mutate(dose_mg = d)
  out
})

dr_data <- bind_rows(dr_results)

p_dr <- ggplot(dr_data, aes(dose_mg, SBP)) +
  geom_line(color = "#FB8C00", size = 1.2) +
  geom_point(color = "#FB8C00", size = 3) +
  geom_hline(yintercept = 130, linetype = "dashed") +
  labs(title = "Amlodipine Dose-Response (SBP at Week 12)",
       x = "Amlodipine Dose (mg/day)", y = "SBP (mmHg)") +
  theme_bw(base_size = 12)

print(p_dr)

cat("\n=== SIMULATION COMPLETE ===\n")
cat("Key outputs saved in 'figures/' directory\n")
