## =============================================================================
## Pulmonary Arterial Hypertension (PAH) — QSP Model
## mrgsolve implementation
##
## Model scope:
##   PK  : sildenafil (PDE5i, 2-cmt oral)
##          bosentan   (ERA dual, 2-cmt oral)
##          treprostinil (SC, 1-cmt)
##   PD  : ET-1 dynamics
##          NO / cGMP pathway
##          Prostacyclin / cAMP pathway
##          Vascular tone (composite)
##          Vascular remodeling (slow, weeks)
##          Pulmonary hemodynamics (PVR, mPAP, CO)
##          Right-ventricular mechanics (Ees, Ea, coupling)
##          Biomarkers  : NT-proBNP, 6MWD
##          WHO FC score
##
## Reference physiology (baseline PAH patient):
##   mPAP  = 40 mmHg (normal ≤20 mmHg)
##   PVR   = 800 dyne·s·cm⁻⁵ (normal ~70-100)
##   CO    = 4.0 L/min (reduced from ~5.5 normal)
##   PAWP  = 8 mmHg
##   6MWD  = 350 m (normal ~550 m)
##   NT-proBNP = 1000 pg/mL (elevated)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

## ─────────────────────────────────────────────────────────────────────────────
## mrgsolve model code block
## ─────────────────────────────────────────────────────────────────────────────
pah_code <- '
$PROB
  PAH QSP Model — ET-1 / NO-cGMP / PGI2-cAMP / Remodeling / Hemodynamics
  Version 1.0

$PARAM @annotated
  // ── Dosing flags (1=on, 0=off) ──────────────────────────────
  DOSE_SIL : 0   : Sildenafil dosing flag (1=yes)
  DOSE_BOS : 0   : Bosentan dosing flag (1=yes)
  DOSE_TRE : 0   : Treprostinil dosing flag (1=yes)

  // ── Sildenafil PK (2-cmt oral) ───────────────────────────────
  // Ref: Muirhead GJ, Br J Clin Pharmacol 2002
  KA_SIL   : 0.90  : Absorption rate constant (1/h)
  CL_SIL   : 40.0  : Clearance (L/h)
  V2_SIL   : 100.0 : Central volume (L)
  Q_SIL    : 18.0  : Inter-compartmental clearance (L/h)
  V3_SIL   : 200.0 : Peripheral volume (L)
  F1_SIL   : 0.40  : Oral bioavailability
  // PD parameters
  EC50_SIL : 30.0  : EC50 for PDE5 inhibition (ng/mL)
  EMAX_SIL : 0.95  : Emax for PDE5 inhibition (fraction)

  // ── Bosentan PK (2-cmt oral) ─────────────────────────────────
  // Ref: van Giersbergen PL, Clin Pharmacokinet 2002
  KA_BOS   : 0.60  : Absorption rate constant (1/h)
  CL_BOS   : 18.0  : Clearance (L/h)
  V2_BOS   : 42.0  : Central volume (L)
  Q_BOS    : 8.0   : Inter-compartmental clearance (L/h)
  V3_BOS   : 50.0  : Peripheral volume (L)
  F1_BOS   : 0.50  : Oral bioavailability
  // PD parameters (ETA blockade)
  EC50_BOS : 150.0 : EC50 for ETA/ETB block (ng/mL)
  EMAX_BOS : 0.90  : Emax for receptor blockade (fraction)

  // ── Treprostinil PK (SC, 1-cmt) ──────────────────────────────
  // Ref: Wade M, J Clin Pharmacol 2004
  KA_TRE   : 0.80  : SC absorption rate (1/h)
  CL_TRE   : 38.0  : Clearance (L/h)
  V2_TRE   : 14.0  : Volume of distribution (L)
  F1_TRE   : 1.00  : SC bioavailability
  // PD parameters (IP receptor agonism → cAMP)
  EC50_TRE : 0.50  : EC50 for cAMP stimulation (ng/mL)
  EMAX_TRE : 0.80  : Emax for cAMP effect (fraction)

  // ── ET-1 dynamics ────────────────────────────────────────────
  KSY_ET1  : 0.60  : ET-1 synthesis rate (pmol/L/h) at baseline
  KDG_ET1  : 0.12  : ET-1 degradation/clearance rate (1/h)
  ET1_0    : 5.0   : Baseline plasma ET-1 (pmol/L)
  STIM_ET1 : 1.50  : Disease-driven ET-1 synthesis multiplier (PAH vs normal)
  ETB_FRAC : 0.30  : Fraction of ET-1 cleared by ETB-EC receptor

  // ── cGMP dynamics ────────────────────────────────────────────
  KPROD_CGMP : 0.50  : Basal sGC activity → cGMP production rate (nmol/L/h)
  KDEG_CGMP  : 0.40  : PDE5-mediated cGMP degradation rate (1/h)
  CGMP_0     : 1.25  : Baseline cGMP (nmol/L, PAH reduced)
  NO_EFF     : 0.60  : Relative NO-driven sGC stimulation in PAH (0–1, reduced)
  SGC_STIM   : 1.0   : Riociguat sGC stimulation multiplier (1 = none, >1 = drug)

  // ── cAMP dynamics ────────────────────────────────────────────
  KPROD_CAMP : 0.80  : Basal AC activity → cAMP production (nmol/L/h)
  KDEG_CAMP  : 0.60  : PDE3/4-mediated cAMP degradation rate (1/h)
  CAMP_0     : 1.33  : Baseline cAMP (nmol/L, PAH reduced)
  PGI2_EFF   : 0.50  : Relative PGI2-driven cAMP stimulation in PAH (0–1)

  // ── Vascular tone (dimensionless, 0=fully dilated, 1=max constrict) ─────
  VT_BASE    : 0.70  : Baseline vascular tone in PAH (0.35 = normal)
  KET1_VT    : 0.15  : ET-1 contribution coefficient to vascular tone
  KCGMP_VT   : 0.20  : cGMP vasodilatory coefficient
  KCAMP_VT   : 0.15  : cAMP vasodilatory coefficient
  ET1_REF    : 5.0   : Reference ET-1 level for tone calculation (pmol/L)
  CGMP_REF   : 1.25  : Reference cGMP for tone calculation (nmol/L)
  CAMP_REF   : 1.33  : Reference cAMP for tone calculation (nmol/L)

  // ── Vascular Remodeling Index (VRI, slow dynamics) ───────────
  VRI_0      : 0.40  : Baseline VRI in PAH (0=no remod, 1=max)
  KGROW_VRI  : 0.002 : VRI growth rate (1/h) driven by pathological stimuli
  KREG_VRI   : 0.001 : VRI regression rate (1/h) driven by treatment
  VRI_MAX    : 0.95  : Maximum VRI (near complete occlusion)

  // ── Hemodynamics ─────────────────────────────────────────────
  PVR_BASE   : 800.0  : Baseline PVR in PAH (dyne·s·cm⁻⁵)
  PVR_NORM   : 80.0   : Normal PVR (used for tone-scaling anchor)
  PAWP_0     : 8.0    : Pulmonary artery wedge pressure (mmHg, fixed)
  CO_0       : 4.0    : Baseline CO in PAH (L/min)
  ALPHA_CO   : 0.30   : Sensitivity of CO to RV coupling change
  HR_0       : 82.0   : Baseline heart rate (bpm)

  // ── RV mechanics ─────────────────────────────────────────────
  Ees_0      : 0.70   : Baseline RV end-systolic elastance (mmHg/mL)
  Ea_BASE    : 0.75   : Baseline effective arterial elastance (mmHg/mL, PAH)
  Ees_HYP    : 0.40   : Max Ees increase from hypertrophy (fractional)
  TAU_HYP    : 720.0  : Time constant of RV hypertrophy (h, ~30 days)
  Ees_DIL    : 0.30   : Ees reduction when decompensated (fractional)
  COUP_THRESH: 0.80   : Ees/Ea threshold below which dilation begins

  // ── NT-proBNP ─────────────────────────────────────────────────
  BNP_BASE   : 1000.0 : Baseline NT-proBNP in PAH (pg/mL)
  BNP_SENS   : 2.5    : Sensitivity of NT-proBNP to RVSP change
  RVSP_REF   : 60.0   : Reference RVSP for NT-proBNP (mmHg)

  // ── 6MWD ─────────────────────────────────────────────────────
  WALK_BASE  : 350.0  : Baseline 6MWD (m)
  WALK_MAX   : 550.0  : Maximum achievable 6MWD (m, healthy normal)
  WALK_SENS  : 1.8    : Sensitivity of 6MWD to PVR change (fractional)

$CMT @annotated
  // Sildenafil compartments
  DEPOT_SIL : Sildenafil GI depot (mg)
  CENT_SIL  : Sildenafil central plasma (mg)
  PERI_SIL  : Sildenafil peripheral (mg)

  // Bosentan compartments
  DEPOT_BOS : Bosentan GI depot (mg)
  CENT_BOS  : Bosentan central plasma (mg)
  PERI_BOS  : Bosentan peripheral (mg)

  // Treprostinil compartments
  SC_TRE    : Treprostinil SC depot (mcg)
  CENT_TRE  : Treprostinil central plasma (mcg)

  // PD state variables
  ET1_PD    : Plasma ET-1 (pmol/L)
  CGMP_PD   : Intracellular cGMP (nmol/L)
  CAMP_PD   : Intracellular cAMP (nmol/L)
  VRI_PD    : Vascular Remodeling Index (0–1)
  RV_HYP    : RV hypertrophy index (0=none, 1=max)

$MAIN
  // ── Derived PK concentrations (ng/mL or mcg/L) ──────────────
  double C_SIL = CENT_SIL / V2_SIL * 1000.0;    // ng/mL
  double C_BOS = CENT_BOS / V2_BOS * 1000.0;    // ng/mL
  double C_TRE = CENT_TRE / V2_TRE * 1000.0;    // ng/mL

  // ── Drug effects (Emax models) ───────────────────────────────
  double EFF_PDE5  = EMAX_SIL * pow(C_SIL, 1.5) /
                     (pow(EC50_SIL, 1.5) + pow(C_SIL, 1.5));   // PDE5 inhib
  double EFF_ERA   = EMAX_BOS * C_BOS / (EC50_BOS + C_BOS);    // ERA block
  double EFF_TRE   = EMAX_TRE * C_TRE / (EC50_TRE + C_TRE);   // IP agonism

  // ── ET-1 equation helpers ───────────────────────────────────
  // ERA reduces ETB clearance and ETA-mediated vascular tone
  // ETB blockade (BOS) reduces plasma ET-1 clearance → transient ET-1 ↑
  double ETB_eff  = ETB_FRAC * (1.0 - EFF_ERA * 0.5);  // partial ETB block
  double SYN_ET1  = KSY_ET1 * STIM_ET1;

  // ── Vascular tone (dimensionless) ───────────────────────────
  // ET-1 increases tone; cGMP and cAMP reduce tone
  double ET1_ratio  = ET1_PD  / ET1_REF;
  double CGMP_ratio = CGMP_PD / CGMP_REF;
  double CAMP_ratio = CAMP_PD / CAMP_REF;
  // ERA blocks ETA-mediated vasoconstriction (reduces KET1_VT contribution)
  double ETA_block  = 1.0 - EFF_ERA;
  double VT = VT_BASE
              + KET1_VT  * (ET1_ratio  - 1.0) * ETA_block
              - KCGMP_VT * (CGMP_ratio - 1.0)
              - KCAMP_VT * (CAMP_ratio - 1.0);
  VT = (VT < 0.05) ? 0.05 : (VT > 1.0) ? 1.0 : VT;  // bound [0.05, 1.0]

  // ── PVR ─────────────────────────────────────────────────────
  // PVR = basal × vascular_tone_ratio × (1 + structural_remodeling)
  double PVR = PVR_BASE * (VT / VT_BASE) * (1.0 + 0.50 * VRI_PD);

  // ── mPAP (Ohm's law: mPAP = PAWP + PVR × CO / 80) ──────────
  // CO_curr updated from RV coupling below; initialise from CO_0
  double PVR_mmHg = PVR / 80.0;  // convert dyne·s·cm⁻⁵ → mmHg·min/L
  // CO estimated dynamically (see below)

  // ── RV end-systolic elastance (Ees) ─────────────────────────
  double Ees = Ees_0 * (1.0 + Ees_HYP * RV_HYP);

  // ── Effective arterial elastance (Ea) ────────────────────────
  // Ea ≈ mPAP / SV ≈ (PVR_mmHg × HR_0 / 60)
  double Ea = Ea_BASE * (PVR / PVR_BASE);

  // ── RV-PA coupling ratio ─────────────────────────────────────
  double COUP = (Ea > 0) ? Ees / Ea : 1.0;
  COUP = (COUP > 3.0) ? 3.0 : COUP;

  // ── CO (cardiac output) ─────────────────────────────────────
  // Compensated: CO maintained when COUP ≥ 1.0
  // Decompensated: CO falls as COUP < COUP_THRESH
  double CO_adj = CO_0;
  if (COUP < 1.0) {
    CO_adj = CO_0 * (0.4 + 0.6 * COUP);  // linear decline below coupling=1
  } else {
    CO_adj = CO_0 * (1.0 + ALPHA_CO * (COUP - 1.0));  // slight ↑ when compensated
  }
  CO_adj = (CO_adj < 1.0) ? 1.0 : CO_adj;

  // ── mPAP ─────────────────────────────────────────────────────
  double mPAP = PAWP_0 + PVR_mmHg * CO_adj;

  // ── RVSP ─────────────────────────────────────────────────────
  double RVSP = mPAP + 5.0;  // RVSP ≈ mPAP + ~5 mmHg (tricuspid grad)

  // ── NT-proBNP ────────────────────────────────────────────────
  double NT_proBNP = BNP_BASE * exp(BNP_SENS * (RVSP - RVSP_REF) / RVSP_REF);
  NT_proBNP = (NT_proBNP < 50) ? 50 : NT_proBNP;

  // ── 6MWD ─────────────────────────────────────────────────────
  double WALK6 = WALK_MAX - (WALK_MAX - WALK_BASE) *
                 pow(PVR / PVR_BASE, WALK_SENS) /
                 (1.0 + ALPHA_CO * fmax(COUP - 1.0, 0));
  WALK6 = (WALK6 < 50.0) ? 50.0 : (WALK6 > 650.0) ? 650.0 : WALK6;

  // ── WHO Functional Class (numeric 1–4) ────────────────────────
  double WHO_FC;
  if      (mPAP < 25 && WALK6 > 500)    WHO_FC = 1;
  else if (mPAP < 35 && WALK6 > 400)    WHO_FC = 2;
  else if (mPAP < 50 && WALK6 > 250)    WHO_FC = 3;
  else                                    WHO_FC = 4;

  // ── ET-1 steady-state (for initial conditions) ───────────────
  ET1_PD_0  = SYN_ET1 / (KDG_ET1 + ETB_FRAC * KDG_ET1);
  CGMP_PD_0 = CGMP_0;
  CAMP_PD_0 = CAMP_0;
  VRI_PD_0  = VRI_0;
  RV_HYP_0  = 0.3;  // PAH patients typically have some RV hypertrophy

$ODE
  // ════════════════════════════════════════════════════════════
  // SILDENAFIL PK (2-cmt first-order oral)
  // ════════════════════════════════════════════════════════════
  double C_SIL_od = CENT_SIL / V2_SIL * 1000.0;
  double EFF_SIL  = EMAX_SIL * pow(C_SIL_od, 1.5) /
                    (pow(EC50_SIL, 1.5) + pow(C_SIL_od, 1.5));

  dxdt_DEPOT_SIL = -KA_SIL * DEPOT_SIL;
  dxdt_CENT_SIL  =  KA_SIL * DEPOT_SIL * F1_SIL
                    - (CL_SIL + Q_SIL) / V2_SIL * CENT_SIL
                    + Q_SIL / V3_SIL * PERI_SIL;
  dxdt_PERI_SIL  =  Q_SIL / V2_SIL * CENT_SIL
                    - Q_SIL / V3_SIL * PERI_SIL;

  // ════════════════════════════════════════════════════════════
  // BOSENTAN PK (2-cmt first-order oral)
  // ════════════════════════════════════════════════════════════
  double C_BOS_od = CENT_BOS / V2_BOS * 1000.0;
  double EFF_BOS  = EMAX_BOS * C_BOS_od / (EC50_BOS + C_BOS_od);

  dxdt_DEPOT_BOS = -KA_BOS * DEPOT_BOS;
  dxdt_CENT_BOS  =  KA_BOS * DEPOT_BOS * F1_BOS
                    - (CL_BOS + Q_BOS) / V2_BOS * CENT_BOS
                    + Q_BOS / V3_BOS * PERI_BOS;
  dxdt_PERI_BOS  =  Q_BOS / V2_BOS * CENT_BOS
                    - Q_BOS / V3_BOS * PERI_BOS;

  // ════════════════════════════════════════════════════════════
  // TREPROSTINIL PK (1-cmt SC)
  // ════════════════════════════════════════════════════════════
  double C_TRE_od = CENT_TRE / V2_TRE * 1000.0;
  double EFF_TRE2 = EMAX_TRE * C_TRE_od / (EC50_TRE + C_TRE_od);

  dxdt_SC_TRE   = -KA_TRE * SC_TRE;
  dxdt_CENT_TRE =  KA_TRE * SC_TRE * F1_TRE
                   - CL_TRE / V2_TRE * CENT_TRE;

  // ════════════════════════════════════════════════════════════
  // ET-1 DYNAMICS
  // Synthesis stimulated by disease state; degradation via NEP
  // and ETB-EC receptor clearance (ETB blocked by ERA)
  // ════════════════════════════════════════════════════════════
  double SYN_ode   = KSY_ET1 * STIM_ET1;
  double ETB_clr   = ETB_FRAC * KDG_ET1 * (1.0 - EFF_BOS * 0.5);
  double DG_ET1    = (KDG_ET1 * (1.0 - ETB_FRAC) + ETB_clr) * ET1_PD;

  dxdt_ET1_PD  = SYN_ode - DG_ET1;

  // ════════════════════════════════════════════════════════════
  // cGMP DYNAMICS
  // Production: sGC activated by NO (reduced in PAH due to ↓NO/↑ROS)
  //             and by riociguat (SGC_STIM parameter)
  // Degradation: PDE5 (inhibited by sildenafil/tadalafil)
  // ════════════════════════════════════════════════════════════
  double PROD_cGMP = KPROD_CGMP * NO_EFF * SGC_STIM;
  double DEG_cGMP  = KDEG_CGMP * (1.0 - EFF_SIL) * CGMP_PD;

  dxdt_CGMP_PD = PROD_cGMP - DEG_cGMP;

  // ════════════════════════════════════════════════════════════
  // cAMP DYNAMICS
  // Production: IP receptor (PGI2 endogenous + prostacyclin drugs)
  //             AC activated by IP agonism
  // Degradation: PDE3 / PDE4
  // ════════════════════════════════════════════════════════════
  double IP_stim   = PGI2_EFF + EFF_TRE2;          // endogenous + drug
  double PROD_cAMP = KPROD_CAMP * IP_stim / PGI2_EFF;  // scaled so baseline OK
  double DEG_cAMP  = KDEG_CAMP * CAMP_PD;

  dxdt_CAMP_PD = PROD_cAMP - DEG_cAMP;

  // ════════════════════════════════════════════════════════════
  // VASCULAR REMODELING INDEX (VRI) — slow dynamics (weeks)
  // Growth driven by: high ET-1, low cGMP/cAMP, inflammation
  // Regression driven by: treatment (ERA + PDE5i + IP agonist)
  // ════════════════════════════════════════════════════════════
  // Proliferative stimulus (0→1: 0=treated, 1=maximal disease)
  double ET1_ratio2  = ET1_PD  / ET1_REF;
  double CGMP_ratio2 = CGMP_PD / CGMP_REF;
  double CAMP_ratio2 = CAMP_PD / CAMP_REF;

  double PROL_stim = fmax(0.0, 0.5 * (ET1_ratio2 - 1.0)
                         - 0.3 * (CGMP_ratio2 - 1.0)
                         - 0.2 * (CAMP_ratio2 - 1.0));
  double REGR_drug = EFF_BOS * 0.4 + EFF_SIL * 0.3 + EFF_TRE2 * 0.3;

  dxdt_VRI_PD = KGROW_VRI * PROL_stim * (VRI_MAX - VRI_PD)
               - KREG_VRI  * REGR_drug  * VRI_PD;

  // ════════════════════════════════════════════════════════════
  // RV HYPERTROPHY INDEX (slow dynamics, weeks–months)
  // Driven by: Ees/Ea ratio < 1 (increased afterload)
  // ════════════════════════════════════════════════════════════
  // Recalculate coupling in ODE context
  double PVR_ode    = PVR_BASE * (VT / VT_BASE) * (1.0 + 0.5 * VRI_PD);
  double Ea_ode     = Ea_BASE * (PVR_ode / PVR_BASE);
  double Ees_ode    = Ees_0 * (1.0 + Ees_HYP * RV_HYP);
  double COUP_ode   = (Ea_ode > 0) ? Ees_ode / Ea_ode : 1.0;
  COUP_ode = (COUP_ode > 3.0) ? 3.0 : COUP_ode;

  // Hypertrophy stimulus: mismatch between afterload and contractility
  double HYP_stim = fmax(0.0, 1.0 - COUP_ode);
  double HYP_regr = fmax(0.0, COUP_ode - 1.0) * 0.5;

  dxdt_RV_HYP = (1.0 / TAU_HYP) * (HYP_stim * (1.0 - RV_HYP)
                                   - HYP_regr * RV_HYP);

$CAPTURE
  // Capture key derived variables for output
  PVR CGMP_PD CAMP_PD ET1_PD VRI_PD RV_HYP
  mPAP RVSP CO_adj COUP NT_proBNP WALK6 WHO_FC
  EFF_SIL EFF_BOS EFF_TRE2
  C_SIL_od C_BOS_od C_TRE_od
'

## ─────────────────────────────────────────────────────────────────────────────
## Compile model
## ─────────────────────────────────────────────────────────────────────────────
mod <- mcode("PAH_QSP", pah_code)

## ─────────────────────────────────────────────────────────────────────────────
## Helper: PVR (needs to be computed outside ODE for initial capture)
## ─────────────────────────────────────────────────────────────────────────────
## Because MAIN and CAPTURE don't retain VT without an explicit $MAIN capture
## we compute VT in $MAIN and tag it, but for simplicity here we rely on
## the $CAPTURE block above and verify outputs in simulation.

## ─────────────────────────────────────────────────────────────────────────────
## Simulation 1: No treatment (natural disease course over 52 weeks)
## ─────────────────────────────────────────────────────────────────────────────
sim_untreated <- mod %>%
  param(DOSE_SIL = 0, DOSE_BOS = 0, DOSE_TRE = 0, STIM_ET1 = 1.5) %>%
  mrgsim(end = 52 * 7 * 24, delta = 24) %>%  # 52 weeks, 24h steps
  as_tibble()

## ─────────────────────────────────────────────────────────────────────────────
## Simulation 2: Sildenafil monotherapy (20 mg TID)
## TID dosing via evid=1, amt=20mg, ii=8h
## ─────────────────────────────────────────────────────────────────────────────
e_sil <- ev(amt = 20, ii = 8, addl = 52*7*3 - 1, cmt = "DEPOT_SIL", evid = 1)

sim_sil <- mod %>%
  param(DOSE_SIL = 1, DOSE_BOS = 0, DOSE_TRE = 0) %>%
  mrgsim(ev = e_sil, end = 52 * 7 * 24, delta = 24) %>%
  as_tibble()

## ─────────────────────────────────────────────────────────────────────────────
## Simulation 3: Bosentan monotherapy (125 mg BID, after 4-week titration)
## ─────────────────────────────────────────────────────────────────────────────
# 4-week titration at 62.5 mg BID, then 125 mg BID
e_bos_titrate <- ev(amt = 62.5, ii = 12, addl = 4*7*2 - 1, cmt = "DEPOT_BOS")
e_bos_maint   <- ev(amt = 125.0, ii = 12, addl = 48*7*2 - 1,
                    cmt = "DEPOT_BOS", time = 4*7*24)
e_bos <- ev_seq(e_bos_titrate, e_bos_maint)

sim_bos <- mod %>%
  param(DOSE_SIL = 0, DOSE_BOS = 1, DOSE_TRE = 0) %>%
  mrgsim(ev = e_bos, end = 52 * 7 * 24, delta = 24) %>%
  as_tibble()

## ─────────────────────────────────────────────────────────────────────────────
## Simulation 4: Combination — ambrisentan + tadalafil (AMBITION-like)
## ERA approximated by bosentan params; tadalafil by sildenafil params (t½ diff)
## ─────────────────────────────────────────────────────────────────────────────
e_era_combo <- ev(amt = 125.0, ii = 24, addl = 52*7 - 1, cmt = "DEPOT_BOS")
e_pde5_combo <- ev(amt = 40.0, ii = 24, addl = 52*7 - 1, cmt = "DEPOT_SIL")
e_combo <- ev_seq(e_era_combo, e_pde5_combo)

sim_combo <- mod %>%
  param(DOSE_SIL = 1, DOSE_BOS = 1, DOSE_TRE = 0,
        # Tadalafil-like PK (longer t½)
        KA_SIL = 0.70, CL_SIL = 7.0, V2_SIL = 63.0,
        Q_SIL = 3.5, V3_SIL = 100.0) %>%
  mrgsim(ev = e_combo, end = 52 * 7 * 24, delta = 24) %>%
  as_tibble()

## ─────────────────────────────────────────────────────────────────────────────
## Simulation 5: Triple therapy — ERA + PDE5i + treprostinil SC
## ─────────────────────────────────────────────────────────────────────────────
e_triple_era  <- ev(amt = 125.0, ii = 24, addl = 52*7 - 1, cmt = "DEPOT_BOS")
e_triple_sil  <- ev(amt = 40.0,  ii = 24, addl = 52*7 - 1, cmt = "DEPOT_SIL")
# Treprostinil SC: ~40 ng/mL target. 1.25 mcg/kg/min → ~90 mcg/h for 70kg patient
e_triple_tre  <- ev(amt = 90*24, ii = 24, addl = 52*7 - 1, cmt = "SC_TRE")
e_triple      <- ev_seq(e_triple_era, e_triple_sil, e_triple_tre)

sim_triple <- mod %>%
  param(DOSE_SIL = 1, DOSE_BOS = 1, DOSE_TRE = 1) %>%
  mrgsim(ev = e_triple, end = 52 * 7 * 24, delta = 24) %>%
  as_tibble()

## ─────────────────────────────────────────────────────────────────────────────
## Combine results for plotting
## ─────────────────────────────────────────────────────────────────────────────
results <- bind_rows(
  sim_untreated %>% mutate(Group = "No treatment"),
  sim_sil       %>% mutate(Group = "Sildenafil 20mg TID"),
  sim_bos       %>% mutate(Group = "Bosentan 125mg BID"),
  sim_combo     %>% mutate(Group = "ERA + PDE5i (combo)"),
  sim_triple    %>% mutate(Group = "Triple therapy")
) %>%
  mutate(
    Week = time / (7 * 24),
    Group = factor(Group, levels = c(
      "No treatment", "Sildenafil 20mg TID",
      "Bosentan 125mg BID", "ERA + PDE5i (combo)", "Triple therapy"
    ))
  )

## ─────────────────────────────────────────────────────────────────────────────
## ggplot2 visualisation
## ─────────────────────────────────────────────────────────────────────────────
color_pal <- c(
  "No treatment"         = "#C0392B",
  "Sildenafil 20mg TID"  = "#1ABC9C",
  "Bosentan 125mg BID"   = "#2980B9",
  "ERA + PDE5i (combo)"  = "#8E44AD",
  "Triple therapy"       = "#27AE60"
)

theme_qsp <- theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title    = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#2C3E50"),
    strip.text       = element_text(colour = "white", face = "bold")
  )

## --- Plot 1: PVR over time ---
p_pvr <- ggplot(results, aes(x = Week, y = PVR, colour = Group)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 51, y = 85, label = "Normal PVR ≈ 80", size = 3.5, colour = "grey50") +
  scale_colour_manual(values = color_pal) +
  labs(title = "Pulmonary Vascular Resistance over Time",
       x = "Week", y = "PVR (dyne·s·cm⁻⁵)") +
  theme_qsp

## --- Plot 2: mPAP over time ---
p_mpap <- ggplot(results, aes(x = Week, y = mPAP, colour = Group)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 20, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 51, y = 21.5, label = "Normal mPAP ≤ 20 mmHg", size = 3.5, colour = "grey50") +
  scale_colour_manual(values = color_pal) +
  labs(title = "Mean Pulmonary Arterial Pressure over Time",
       x = "Week", y = "mPAP (mmHg)") +
  theme_qsp

## --- Plot 3: 6MWD over time ---
p_6mwd <- ggplot(results, aes(x = Week, y = WALK6, colour = Group)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 550, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 51, y = 558, label = "Normal ≈ 550 m", size = 3.5, colour = "grey50") +
  scale_colour_manual(values = color_pal) +
  labs(title = "6-Minute Walk Distance over Time",
       x = "Week", y = "6MWD (m)") +
  theme_qsp

## --- Plot 4: NT-proBNP over time ---
p_bnp <- ggplot(results, aes(x = Week, y = NT_proBNP, colour = Group)) +
  geom_line(size = 1) +
  scale_y_log10() +
  geom_hline(yintercept = 300, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 51, y = 280, label = "Target <300 pg/mL", size = 3.5, colour = "grey50") +
  scale_colour_manual(values = color_pal) +
  labs(title = "NT-proBNP Biomarker over Time",
       x = "Week", y = "NT-proBNP (pg/mL, log scale)") +
  theme_qsp

## --- Plot 5: cGMP and cAMP (pharmacodynamic markers) ---
p_second_mess <- results %>%
  filter(Group %in% c("No treatment", "Sildenafil 20mg TID",
                      "ERA + PDE5i (combo)", "Triple therapy")) %>%
  tidyr::pivot_longer(cols = c(CGMP_PD, CAMP_PD),
                      names_to  = "Messenger",
                      values_to = "Conc") %>%
  mutate(Messenger = recode(Messenger,
                             CGMP_PD = "cGMP (nmol/L)",
                             CAMP_PD = "cAMP (nmol/L)")) %>%
  ggplot(aes(x = Week, y = Conc, colour = Group, linetype = Messenger)) +
  geom_line(size = 1) +
  scale_colour_manual(values = color_pal) +
  labs(title = "Second Messengers (cGMP & cAMP) over Time",
       x = "Week", y = "Concentration (nmol/L)") +
  theme_qsp

## --- Plot 6: RV coupling ---
p_coup <- ggplot(results, aes(x = Week, y = COUP, colour = Group)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1.0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = COUP_THRESH, linetype = "dotted", colour = "#E74C3C") +
  annotate("text", x = 51, y = 1.05, label = "Coupling = 1.0", size = 3.5) +
  annotate("text", x = 51, y = 0.75, label = "Decompensation threshold",
           size = 3.5, colour = "#E74C3C") +
  scale_colour_manual(values = color_pal) +
  labs(title = "RV-PA Coupling (Ees/Ea) over Time",
       x = "Week", y = "Ees/Ea Ratio") +
  theme_qsp

## --- Plot 7: ET-1 plasma ---
p_et1 <- ggplot(results, aes(x = Week, y = ET1_PD, colour = Group)) +
  geom_line(size = 1) +
  scale_colour_manual(values = color_pal) +
  labs(title = "Plasma ET-1 over Time",
       x = "Week", y = "ET-1 (pmol/L)") +
  theme_qsp

## --- Arrange all plots ---
if (requireNamespace("gridExtra", quietly = TRUE)) {
  library(gridExtra)
  main_fig <- grid.arrange(p_pvr, p_mpap, p_6mwd, p_bnp,
                            p_second_mess, p_coup,
                            ncol = 2,
                            top = "PAH QSP Model — Treatment Comparison (52 weeks)")
}

## ─────────────────────────────────────────────────────────────────────────────
## Summary table at Week 12 and Week 52
## ─────────────────────────────────────────────────────────────────────────────
make_summary <- function(df, week_target) {
  df %>%
    filter(abs(Week - week_target) < 0.6) %>%
    group_by(Group) %>%
    slice(1) %>%
    select(Group, PVR, mPAP, WALK6, NT_proBNP, COUP, WHO_FC) %>%
    rename(
      "PVR (dyn·s/cm⁵)" = PVR,
      "mPAP (mmHg)"      = mPAP,
      "6MWD (m)"         = WALK6,
      "NT-proBNP (pg/mL)" = NT_proBNP,
      "Ees/Ea"            = COUP,
      "WHO FC"            = WHO_FC
    ) %>%
    mutate(across(where(is.numeric), ~round(.x, 1)))
}

cat("\n====  Week 12 Summary  ====\n")
print(as.data.frame(make_summary(results, 12)))

cat("\n====  Week 52 Summary  ====\n")
print(as.data.frame(make_summary(results, 52)))

## ─────────────────────────────────────────────────────────────────────────────
## Sensitivity Analysis: EC50 sildenafil effect on cGMP at steady state
## ─────────────────────────────────────────────────────────────────────────────
ec50_range <- c(10, 20, 30, 50, 100)
sa_list <- lapply(ec50_range, function(ec50) {
  e_sa <- ev(amt = 20, ii = 8, addl = 52*7*3, cmt = "DEPOT_SIL")
  mod %>%
    param(EC50_SIL = ec50) %>%
    mrgsim(ev = e_sa, end = 52 * 7 * 24, delta = 24) %>%
    as_tibble() %>%
    mutate(EC50 = paste0("EC50=", ec50, " ng/mL"),
           Week = time / (7 * 24))
})
sa_df <- bind_rows(sa_list)

p_sa <- ggplot(sa_df, aes(x = Week, y = CGMP_PD, colour = EC50)) +
  geom_line(size = 1) +
  scale_colour_viridis_d(option = "plasma") +
  labs(title = "Sensitivity Analysis: EC50 Sildenafil → cGMP Response",
       x = "Week", y = "cGMP (nmol/L)") +
  theme_qsp

print(p_sa)

message("PAH QSP simulation complete.")
