## ============================================================
## Systemic Sclerosis (SSc) QSP Model — mrgsolve Implementation
## ============================================================
## Disease: Systemic Sclerosis (Scleroderma)
## Model scope:
##   PK: Nintedanib (FGFR/PDGFR/VEGFRi), Tocilizumab (anti-IL-6R),
##       Mycophenolate mofetil (MMF), Bosentan (ERA), Iloprost (PGI2-analog)
##   PD: Fibrosis (mRSS, FVC), Vascular (ET-1, PVR, mPAP, 6MWD),
##       Immune (TGF-β, IL-6, B cells, Th17), PAH progression
## Parameter calibration references:
##   - SENSCIS trial (Distler JHW, NEJM 2019) – nintedanib in SSc-ILD
##   - FOCUSSS / TRANSFORM trial (tocilizumab) – FVC preservation
##   - faSScinate trial (tocilizumab, mRSS)
##   - Bosentan RAPIDS-1/2 (digital ulcers)
##   - STEP trial (iloprost)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ============================================================
## MODEL CODE BLOCK
## ============================================================

ssc_code <- '
$PROB Systemic Sclerosis QSP Model
  PK compartments: Nintedanib (oral 2-cpt), Tocilizumab (IV/SC mAb 2-cpt),
                   MMF/MPA (oral 1-cpt), Bosentan (oral 1-cpt), Iloprost (inh 1-cpt)
  PD compartments: TGF-beta, IL-6, Th17 cells, B cells (naive/memory),
                   Fibroblast activation (FAct), Myofibroblasts (Myo),
                   Collagen (Col1, Col3), ECM accumulation,
                   Skin fibrosis (mRSS), FVC % predicted, DLCO,
                   ET-1, PVR, mPAP, RV function, 6MWD, WHO-FC proxy,
                   Endothelial integrity (Endo), NO, Prostacyclin (PGI2)

$PARAM
  // ── Body & dosing ──────────────────────────────────────
  BW        = 65       // body weight (kg)

  // ── Nintedanib PK (2-cpt oral, based on SENSCIS pop-PK) ──
  // Ref: Dallmann A et al, Clin Pharmacokinet 2020
  F_NINT    = 0.047    // bioavailability (oral)
  KA_NINT   = 1.2      // absorption rate constant (1/h)
  CL_NINT   = 1050     // clearance (L/h)
  V1_NINT   = 985      // central volume (L)
  Q_NINT    = 420      // inter-compartment CL (L/h)
  V2_NINT   = 2100     // peripheral volume (L)

  // ── Tocilizumab PK (2-cpt IV mAb) ─────────────────────
  // Ref: Frey N et al, J Clin Pharmacol 2010
  CL_TCZ    = 0.22     // linear CL (L/h) — dose-dependent target-mediated
  V1_TCZ    = 3.7      // central volume (L)
  Q_TCZ     = 0.09     // inter-cpt CL (L/h)
  V2_TCZ    = 3.0      // peripheral volume (L)

  // ── MMF/MPA PK (1-cpt, active moiety MPA) ─────────────
  // Ref: Kiberd BA et al, Transplantation 2006
  F_MMF     = 0.94     // bioavailability
  KA_MMF    = 1.5      // absorption (1/h)
  CL_MMF    = 15.6     // CL (L/h)
  V_MMF     = 30.0     // volume (L)

  // ── Bosentan PK (1-cpt oral) ──────────────────────────
  // Ref: van Giersbergen PLM, Clin Pharmacokinet 2003
  F_BOSEN   = 0.50
  KA_BOSEN  = 0.8
  CL_BOSEN  = 15.5
  V_BOSEN   = 78.0

  // ── Iloprost PK (inhaled, 1-cpt) ──────────────────────
  // Ref: Hoeper MM et al. PAH inh iloprost data
  F_ILOP    = 0.17     // inhaled bioavailability
  KA_ILOP   = 5.0
  CL_ILOP   = 75.0
  V_ILOP    = 20.0

  // ── Immune/Cytokine PD ─────────────────────────────────
  kprod_TGFB  = 0.05   // TGF-beta basal production (nmol/h)
  kdeg_TGFB   = 0.08   // TGF-beta degradation (1/h)
  TGFB_base   = 0.625  // basal TGF-beta (nmol/L) ~5 ng/mL

  kprod_IL6   = 0.20   // IL-6 basal prod (pmol/h)
  kdeg_IL6    = 0.25   // IL-6 degradation (1/h)
  IL6_base    = 0.80   // basal IL-6 (pmol/L) ~8 pg/mL (elevated in SSc)

  // Th17 cells
  kprod_Th17  = 0.004
  kdeg_Th17   = 0.010
  Th17_base   = 0.40   // relative units
  IL6_EC50_Th17 = 1.5  // IL-6 EC50 for Th17 induction

  // B cells (naive)
  kprod_Bnv   = 0.008
  kdeg_Bnv    = 0.005
  Bnv_base    = 1.6    // relative units
  BAFF_conc   = 1.0    // BAFF (assumed constant unless belimumab added)

  // ── Fibrosis PD ───────────────────────────────────────
  kact_FAct   = 0.015  // fibroblast activation rate by TGF-beta
  kdeg_FAct   = 0.008  // fibroblast deactivation
  FAct_base   = 0.30   // basal activated fibroblast fraction

  kform_Myo   = 0.012  // myofibroblast formation
  kdeg_Myo    = 0.006
  Myo_base    = 0.15

  ksynth_Col1 = 0.020  // collagen I synthesis
  kdeg_Col1   = 0.003  // collagen I degradation
  Col1_base   = 6.67   // basal normalized col1

  ksynth_Col3 = 0.012
  kdeg_Col3   = 0.003
  Col3_base   = 4.0

  ECM_0       = 1.0    // baseline ECM (normalized)
  kECM        = 0.005  // ECM accumulation from collagen
  kdegECM     = 0.004

  // ── Skin fibrosis (mRSS) ──────────────────────────────
  // mRSS range 0-51; typical dcSSc at diagnosis ~20
  mRSS_0      = 20.0
  kmRSS_up    = 0.0020  // mRSS progression rate driven by ECM
  kmRSS_dn    = 0.0008  // spontaneous improvement
  mRSS_max    = 51.0

  // ── Lung (FVC % predicted) ───────────────────────────
  FVC_0       = 78.0    // baseline FVC % predicted (SSc-ILD typical)
  kFVC_loss   = 0.0012  // FVC annual decline rate (/h) untreated ~2.7%/yr
  kFVC_rec    = 0.0003  // partial recovery rate
  FVC_min     = 30.0    // minimum (death threshold)

  DLCO_0      = 62.0    // baseline DLCO % predicted
  kDLCO_loss  = 0.0015

  // ── Vascular PD ───────────────────────────────────────
  ET1_0       = 3.5     // basal ET-1 (pg/mL) — elevated in SSc ~3-4 pg/mL
  kprod_ET1   = 0.015
  kdeg_ET1    = 0.04

  NO_0        = 0.50    // relative NO (reduced in SSc)
  kprod_NO    = 0.02
  kdeg_NO     = 0.04

  PGI2_0      = 0.30    // relative prostacyclin (reduced in SSc)
  kprod_PGI2  = 0.015
  kdeg_PGI2   = 0.05

  Endo_0      = 0.55    // endothelial integrity (0-1); SSc ~55% baseline
  kdmg_Endo   = 0.003   // endothelial damage rate
  krep_Endo   = 0.002   // endothelial repair rate

  PVR_0       = 350.0   // pulmonary vascular resistance (dyn*s*cm-5)
                        // SSc-PAH typical baseline ~400
  kPVR_up     = 0.0015  // PVR progression
  kPVR_dn     = 0.0005  // PVR improvement

  mPAP_0      = 30.0    // mean PAP mmHg (SSc-PAH borderline elevated)

  SixMWD_0    = 380.0   // 6MWD meters — SSc-PAH at diagnosis
  kSixMWD_loss = 0.0010
  SixMWD_min  = 100.0

  // ── Drug PD effect parameters ─────────────────────────
  // Nintedanib: anti-fibrotic via PDGFR/FGFR inh
  EMAX_NINT_FAct = 0.60  // max inhibition of fibroblast activation
  EC50_NINT_FAct = 85.0  // EC50 nintedanib for fibroblast inh (ng/mL)
                         // Free Cmax ~120 ng/mL after 150 mg BID
  EMAX_NINT_FVC  = 0.55  // max effect on FVC decline reduction
  EC50_NINT_FVC  = 80.0

  // Tocilizumab: anti-IL-6R → ↓Th17, ↓IL-6 effects
  EMAX_TCZ_IL6   = 0.85  // max IL-6 signaling blockade
  EC50_TCZ_IL6   = 2.5   // EC50 (μg/mL) — typical Ctrough 5-10 μg/mL
  EMAX_TCZ_mRSS  = 0.40
  EC50_TCZ_mRSS  = 3.0

  // MMF: ↓Th17/B cell proliferation
  EMAX_MMF_Th17  = 0.50
  EC50_MMF_Th17  = 1.8   // MPA μg/mL — typical trough 1-3.5 μg/mL

  // Bosentan: ↓ET-1 effects via ERA
  EMAX_BOSEN_ET1 = 0.75
  EC50_BOSEN_ET1 = 1.2   // μg/mL — typical Ctrough ~1.5 μg/mL
  EMAX_BOSEN_PVR = 0.60
  EC50_BOSEN_PVR = 1.5

  // Iloprost: ↑PGI2-like effects, ↓PVR
  EMAX_ILOP_PGI2 = 0.80
  EC50_ILOP_PGI2 = 0.5   // ng/mL — inhaled iloprost Cmax ~0.1-2 ng/mL
  EMAX_ILOP_PVR  = 0.50
  EC50_ILOP_PVR  = 0.6

$INIT
  // PK compartments
  A_NINT1 = 0, A_NINT2 = 0, A_NINT_GUT = 0
  A_TCZ1  = 0, A_TCZ2  = 0
  A_MMF_GUT = 0, A_MPA = 0
  A_BOSEN_GUT = 0, A_BOSEN = 0
  A_ILOP_LUNG = 0, A_ILOP = 0

  // Immune/cytokine
  TGFB   = 0.625
  IL6    = 0.80
  Th17   = 0.40
  Bnv    = 1.60

  // Fibrosis
  FAct   = 0.30
  Myo    = 0.15
  Col1   = 6.67
  Col3   = 4.00
  ECM    = 1.00

  // Skin
  mRSS   = 20.0

  // Lung
  FVC    = 78.0
  DLCO   = 62.0

  // Vascular
  ET1    = 3.50
  NO     = 0.50
  PGI2_c = 0.30
  Endo   = 0.55
  PVR    = 350.0
  SixMWD = 380.0

$ODE
  // ────────────────────────────────────────────────────────────
  // PK: Nintedanib (2-cpt oral, MW ~539.6)
  // Dose in mg → conc in ng/mL (V in L/kg approximate total)
  double C_NINT = A_NINT1 / V1_NINT * 1000; // ng/mL

  dxdt_A_NINT_GUT = -KA_NINT * A_NINT_GUT;
  dxdt_A_NINT1    =  KA_NINT * F_NINT * A_NINT_GUT
                   - (CL_NINT + Q_NINT) / V1_NINT * A_NINT1
                   + Q_NINT / V2_NINT * A_NINT2;
  dxdt_A_NINT2    =  Q_NINT / V1_NINT * A_NINT1
                   - Q_NINT / V2_NINT * A_NINT2;

  // ────────────────────────────────────────────────────────────
  // PK: Tocilizumab (2-cpt IV, conc in μg/mL)
  double C_TCZ = A_TCZ1 / V1_TCZ;  // μg/mL

  dxdt_A_TCZ1 = -(CL_TCZ + Q_TCZ) * A_TCZ1 / V1_TCZ
                + Q_TCZ * A_TCZ2 / V2_TCZ;
  dxdt_A_TCZ2 =  Q_TCZ * A_TCZ1 / V1_TCZ
               - Q_TCZ * A_TCZ2 / V2_TCZ;

  // ────────────────────────────────────────────────────────────
  // PK: MPA (active moiety of MMF, 1-cpt)
  double C_MPA = A_MPA / V_MMF;   // μg/mL

  dxdt_A_MMF_GUT = -KA_MMF * A_MMF_GUT;
  dxdt_A_MPA     =  KA_MMF * F_MMF * A_MMF_GUT - CL_MMF * C_MPA;

  // ────────────────────────────────────────────────────────────
  // PK: Bosentan (1-cpt oral)
  double C_BOSEN = A_BOSEN / V_BOSEN;  // μg/mL

  dxdt_A_BOSEN_GUT = -KA_BOSEN * A_BOSEN_GUT;
  dxdt_A_BOSEN     =  KA_BOSEN * F_BOSEN * A_BOSEN_GUT
                    - CL_BOSEN * C_BOSEN;

  // ────────────────────────────────────────────────────────────
  // PK: Iloprost (inhaled, 1-cpt, ng/mL)
  double C_ILOP = A_ILOP / V_ILOP;  // ng/mL

  dxdt_A_ILOP_LUNG = -KA_ILOP * A_ILOP_LUNG;
  dxdt_A_ILOP      =  KA_ILOP * F_ILOP * A_ILOP_LUNG
                    - CL_ILOP * C_ILOP;

  // ────────────────────────────────────────────────────────────
  // Drug effect functions (inhibition Hill equations)
  double E_NINT_FAct  = EMAX_NINT_FAct * C_NINT / (EC50_NINT_FAct + C_NINT);
  double E_NINT_FVC   = EMAX_NINT_FVC  * C_NINT / (EC50_NINT_FVC  + C_NINT);
  double E_TCZ_IL6    = EMAX_TCZ_IL6   * C_TCZ  / (EC50_TCZ_IL6   + C_TCZ);
  double E_TCZ_mRSS   = EMAX_TCZ_mRSS  * C_TCZ  / (EC50_TCZ_mRSS  + C_TCZ);
  double E_MMF_Th17   = EMAX_MMF_Th17  * C_MPA  / (EC50_MMF_Th17  + C_MPA);
  double E_BOSEN_ET1  = EMAX_BOSEN_ET1 * C_BOSEN / (EC50_BOSEN_ET1 + C_BOSEN);
  double E_BOSEN_PVR  = EMAX_BOSEN_PVR * C_BOSEN / (EC50_BOSEN_PVR + C_BOSEN);
  double E_ILOP_PGI2  = EMAX_ILOP_PGI2 * C_ILOP  / (EC50_ILOP_PGI2 + C_ILOP);
  double E_ILOP_PVR   = EMAX_ILOP_PVR  * C_ILOP  / (EC50_ILOP_PVR  + C_ILOP);

  // ────────────────────────────────────────────────────────────
  // TGF-β: driven by ECM feedback + FAct; reduced by nintedanib indirectly
  double TGFB_stim  = 1.0 + 0.3 * (ECM - 1.0);  // ECM feedback
  double TGFB_prod  = kprod_TGFB * TGFB_stim * (1 - 0.2 * E_NINT_FAct);
  dxdt_TGFB = TGFB_prod - kdeg_TGFB * TGFB;

  // ────────────────────────────────────────────────────────────
  // IL-6: elevated in SSc; inhibited by TCZ (receptor blockade)
  // TCZ blocks IL-6 receptor — represented as ↓ effective IL-6
  double IL6_eff    = IL6 * (1 - E_TCZ_IL6);
  double IL6_stim   = 1.0 + 0.5 * (TGFB / TGFB_base - 1.0);
  dxdt_IL6 = kprod_IL6 * IL6_stim - kdeg_IL6 * IL6;

  // ────────────────────────────────────────────────────────────
  // Th17: driven by IL-6 + TGF-β; inhibited by MMF, TCZ
  double Th17_stim = IL6_eff / (IL6_EC50_Th17 + IL6_eff) + 0.3 * TGFB / TGFB_base;
  dxdt_Th17 = kprod_Th17 * Th17_stim * (1 - E_MMF_Th17)
            - kdeg_Th17 * Th17;

  // ────────────────────────────────────────────────────────────
  // B cells (naive): BAFF-dependent survival; inhibited by MMF
  double Bnv_stim = BAFF_conc * (1 - 0.3 * E_MMF_Th17);
  dxdt_Bnv = kprod_Bnv * Bnv_stim - kdeg_Bnv * Bnv;

  // ────────────────────────────────────────────────────────────
  // Fibroblast activation: TGF-β + IL-6 + Th17 driven
  double FAct_stim = (TGFB / TGFB_base) * (1 + 0.2 * (Th17 / 0.4 - 1.0));
  dxdt_FAct = kact_FAct * FAct_stim * (1 - E_NINT_FAct)
            - kdeg_FAct * FAct;

  // ────────────────────────────────────────────────────────────
  // Myofibroblasts
  dxdt_Myo = kform_Myo * FAct - kdeg_Myo * Myo;

  // ────────────────────────────────────────────────────────────
  // Collagen I & III synthesis
  double col_stim = Myo / Myo_base;
  dxdt_Col1 = ksynth_Col1 * col_stim - kdeg_Col1 * Col1;
  dxdt_Col3 = ksynth_Col3 * col_stim - kdeg_Col3 * Col3;

  // ────────────────────────────────────────────────────────────
  // ECM accumulation (normalized)
  double ECM_input  = kECM * (Col1 / Col1_base + Col3 / Col3_base) / 2.0;
  dxdt_ECM = ECM_input - kdegECM * ECM;

  // ────────────────────────────────────────────────────────────
  // mRSS (skin score 0-51)
  // Driven by ECM; limited by mRSS_max; TCZ reduces progression
  double mRSS_prog_rate = kmRSS_up * (ECM - 1.0) * (1 - E_TCZ_mRSS);
  double mRSS_reg_rate  = kmRSS_dn;
  dxdt_mRSS = (mRSS_prog_rate - mRSS_reg_rate) * mRSS
              * (1 - mRSS / mRSS_max);
  // Protect from going negative
  if (mRSS < 0) dxdt_mRSS = 0;

  // ────────────────────────────────────────────────────────────
  // FVC % predicted — decline driven by lung fibrosis (ECM×IL6)
  // Nintedanib slows decline (E_NINT_FVC)
  double FVC_loss_rate = kFVC_loss * (ECM / 1.0) * (1 - E_NINT_FVC);
  double FVC_rec_rate  = kFVC_rec;
  dxdt_FVC  = -(FVC_loss_rate - FVC_rec_rate) * (FVC - FVC_min);

  // DLCO: more sensitive to vascular disease
  double DLCO_loss_rate = kDLCO_loss * (ECM / 1.0 + 0.5 * (PVR / PVR_0 - 1.0));
  dxdt_DLCO = -DLCO_loss_rate * (DLCO - 30.0);

  // ────────────────────────────────────────────────────────────
  // Endothelial integrity
  double Endo_dmg = kdmg_Endo * (TGFB / TGFB_base + Th17 / 0.4) / 2.0;
  double Endo_rep = krep_Endo * (NO / NO_0 + PGI2_c / PGI2_0) / 2.0
                  * (1 + E_ILOP_PGI2 * 0.5);
  dxdt_Endo = Endo_rep * (1 - Endo) - Endo_dmg * Endo;

  // ────────────────────────────────────────────────────────────
  // ET-1: elevated in SSc; bosentan reduces ET-1 signaling
  double ET1_prod = kprod_ET1 * (1 / (Endo + 0.01)) * (1 - E_BOSEN_ET1);
  dxdt_ET1 = ET1_prod - kdeg_ET1 * ET1;

  // ────────────────────────────────────────────────────────────
  // NO: reduced by endothelial injury; iloprost partially restores
  double NO_prod = kprod_NO * Endo * (1 + 0.3 * E_ILOP_PGI2);
  dxdt_NO = NO_prod - kdeg_NO * NO;

  // ────────────────────────────────────────────────────────────
  // Prostacyclin (PGI2): reduced in SSc; iloprost analog effect
  double PGI2_prod = kprod_PGI2 * Endo * (1 + E_ILOP_PGI2);
  dxdt_PGI2_c = PGI2_prod - kdeg_PGI2 * PGI2_c;

  // ────────────────────────────────────────────────────────────
  // PVR: driven by ET-1, reduced NO/PGI2, ECM-driven remodeling
  // Bosentan: ERA blockade; Iloprost: PGI2 analog vasodilation
  double PVR_ET1_drive = ET1 / ET1_0;
  double PVR_NO_brake  = NO   / NO_0;
  double PVR_PGI_brake = PGI2_c / PGI2_0;
  double PVR_prog  = kPVR_up * PVR_ET1_drive * (1 - E_BOSEN_PVR)
                   * (ECM / 1.0);
  double PVR_reg   = kPVR_dn * (PVR_NO_brake + PVR_PGI_brake) / 2.0
                   * (1 + E_ILOP_PVR);
  dxdt_PVR = (PVR_prog - PVR_reg) * PVR * 0.10;

  // ────────────────────────────────────────────────────────────
  // 6-Minute Walk Distance
  double SixMWD_loss = kSixMWD_loss * (PVR / PVR_0)
                     * (1 - E_ILOP_PVR * 0.5 - E_BOSEN_PVR * 0.3);
  dxdt_SixMWD = -SixMWD_loss * (SixMWD - SixMWD_min);

  // mPAP: estimated from PVR (Chemla formula approximation)
  // mPAP = 0.61 * PVR_indexed + 2 (simplified)

$TABLE
  double C_NINT_obs  = C_NINT;
  double C_TCZ_obs   = C_TCZ;
  double C_MPA_obs   = C_MPA;
  double C_BOSEN_obs = C_BOSEN;
  double C_ILOP_obs  = C_ILOP;

  // Derived endpoints
  double mPAP_est    = 0.61 * PVR / 80.0 + 2.0;  // simplified approximation
  double WHO_FC_est  = (SixMWD > 440) ? 1 :
                       (SixMWD > 315) ? 2 :
                       (SixMWD > 165) ? 3 : 4;

  // CRISS score approximation (ACR/EULAR SSc response index)
  // CRISS based on: mRSS, FVC, DLCO, physician VAS, patient VAS
  double CRISS_approx = 1.0 / (1.0 + exp(-(
    -8.15
    - 0.089 * (mRSS - 20.0)
    + 0.049 * (FVC  - 78.0)
    + 0.032 * (DLCO - 62.0)
  )));

$CAPTURE
  C_NINT_obs C_TCZ_obs C_MPA_obs C_BOSEN_obs C_ILOP_obs
  TGFB IL6 Th17 Bnv
  FAct Myo Col1 Col3 ECM
  mRSS FVC DLCO
  ET1 NO PGI2_c Endo
  PVR SixMWD mPAP_est WHO_FC_est
  CRISS_approx
'

## ============================================================
## Compile model
## ============================================================

ssc_model <- mcode("ssc_qsp", ssc_code)
cat("Model compiled successfully.\n")

## ============================================================
## Dosing Events — 5 Treatment Scenarios
## ============================================================

## Time scale: hours; 1 year = 8760 h; simulate 2 years = 17520 h
SIM_DUR <- 17520   # 2 years in hours
STEP    <- 24      # output every 24 h

# Base disease events (no drug)
ev_none <- ev(time = 0, cmt = 1, amt = 0)

# Scenario 1: Nintedanib 150 mg BID (every 12h)
# amt in mg → compartment A_NINT_GUT
ev_nint <- ev(cmt = "A_NINT_GUT", amt = 150, ii = 12, addl = 1459)

# Scenario 2: Tocilizumab 8 mg/kg IV q4w
# 65 kg × 8 mg/kg = 520 mg = 520,000 μg → directly into A_TCZ1
# iv bolus: bioavailability=1, so rate input into central compartment
ev_tcz <- ev(cmt = "A_TCZ1", amt = 520 * 1000, ii = 672, addl = 25)

# Scenario 3: MMF 3000 mg/day (1500 mg BID)
ev_mmf <- ev(cmt = "A_MMF_GUT", amt = 1500, ii = 12, addl = 1459)

# Scenario 4: Bosentan 125 mg BID (ERA for PAH)
ev_bosen <- ev(cmt = "A_BOSEN_GUT", amt = 125, ii = 12, addl = 1459)

# Scenario 5: Combination Nintedanib + Tocilizumab + MMF
ev_combo <- ev_nint + ev_tcz + ev_mmf

## ============================================================
## Run simulations
## ============================================================

run_sim <- function(model, events, label, dur = SIM_DUR, step = STEP) {
  out <- model %>%
    ev(events) %>%
    mrgsim(end = dur, delta = step) %>%
    as_tibble() %>%
    mutate(scenario = label,
           time_yr  = time / 8760)
  return(out)
}

cat("Running simulations...\n")

sim_none  <- run_sim(ssc_model, ev_none,  "Untreated")
sim_nint  <- run_sim(ssc_model, ev_nint,  "Nintedanib 150 mg BID")
sim_tcz   <- run_sim(ssc_model, ev_tcz,   "Tocilizumab 8 mg/kg q4w")
sim_mmf   <- run_sim(ssc_model, ev_mmf,   "MMF 3000 mg/day")
sim_bosen <- run_sim(ssc_model, ev_bosen, "Bosentan 125 mg BID")
sim_combo <- run_sim(ssc_model, ev_combo, "Nintedanib + TCZ + MMF")

all_sims <- bind_rows(sim_none, sim_nint, sim_tcz, sim_mmf, sim_bosen, sim_combo)

## ============================================================
## Plotting — Key Outcomes
## ============================================================

plot_colors <- c(
  "Untreated"                = "#e74c3c",
  "Nintedanib 150 mg BID"    = "#2980b9",
  "Tocilizumab 8 mg/kg q4w"  = "#27ae60",
  "MMF 3000 mg/day"          = "#f39c12",
  "Bosentan 125 mg BID"      = "#8e44ad",
  "Nintedanib + TCZ + MMF"   = "#1abc9c"
)

# FVC % predicted
p_fvc <- ggplot(all_sims, aes(x = time_yr, y = FVC, color = scenario)) +
  geom_line(size = 1.0) +
  scale_color_manual(values = plot_colors) +
  labs(title = "FVC % Predicted (SSc-ILD)",
       x = "Time (years)", y = "FVC % predicted",
       color = "Treatment") +
  geom_hline(yintercept = 70, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.05, y = 68.5, label = "FVC 70% threshold",
           hjust = 0, size = 3, color = "gray50") +
  ylim(40, 85) + theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

# mRSS
p_mrss <- ggplot(all_sims, aes(x = time_yr, y = mRSS, color = scenario)) +
  geom_line(size = 1.0) +
  scale_color_manual(values = plot_colors) +
  labs(title = "Modified Rodnan Skin Score (mRSS)",
       x = "Time (years)", y = "mRSS (0-51)",
       color = "Treatment") +
  ylim(0, 40) + theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

# PVR
p_pvr <- ggplot(all_sims, aes(x = time_yr, y = PVR, color = scenario)) +
  geom_line(size = 1.0) +
  scale_color_manual(values = plot_colors) +
  labs(title = "Pulmonary Vascular Resistance",
       x = "Time (years)", y = "PVR (dyn·s·cm⁻⁵)",
       color = "Treatment") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

# 6MWD
p_6mwd <- ggplot(all_sims, aes(x = time_yr, y = SixMWD, color = scenario)) +
  geom_line(size = 1.0) +
  scale_color_manual(values = plot_colors) +
  labs(title = "6-Minute Walk Distance",
       x = "Time (years)", y = "6MWD (meters)",
       color = "Treatment") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

# ET-1
p_et1 <- ggplot(all_sims, aes(x = time_yr, y = ET1, color = scenario)) +
  geom_line(size = 1.0) +
  scale_color_manual(values = plot_colors) +
  labs(title = "Endothelin-1 (ET-1)",
       x = "Time (years)", y = "ET-1 (pg/mL)",
       color = "Treatment") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

# TGF-β
p_tgfb <- ggplot(all_sims, aes(x = time_yr, y = TGFB, color = scenario)) +
  geom_line(size = 1.0) +
  scale_color_manual(values = plot_colors) +
  labs(title = "TGF-β1 (fibrosis driver)",
       x = "Time (years)", y = "TGF-β1 (nmol/L)",
       color = "Treatment") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

# PK: Nintedanib
pk_nint <- all_sims %>% filter(scenario == "Nintedanib 150 mg BID")
p_pk_nint <- ggplot(pk_nint %>% filter(time_yr < 0.05),
                    aes(x = time * 1, y = C_NINT_obs)) +
  geom_line(color = "#2980b9", size = 1.2) +
  labs(title = "Nintedanib PK (first 18h)",
       x = "Time (hours)", y = "Plasma Conc (ng/mL)") +
  theme_bw()

# Combine plots
combined_plot <- (p_fvc | p_mrss) / (p_pvr | p_6mwd) / (p_et1 | p_tgfb)
print(combined_plot)

## ============================================================
## Summary table at 2 years
## ============================================================

summary_2yr <- all_sims %>%
  filter(abs(time_yr - max(time_yr)) < 0.01) %>%
  group_by(scenario) %>%
  summarise(
    FVC_2yr    = round(mean(FVC), 1),
    mRSS_2yr   = round(mean(mRSS), 1),
    PVR_2yr    = round(mean(PVR), 0),
    SixMWD_2yr = round(mean(SixMWD), 0),
    ET1_2yr    = round(mean(ET1), 2),
    TGFB_2yr   = round(mean(TGFB), 3),
    mPAP_2yr   = round(mean(mPAP_est), 1),
    WHO_FC_2yr = round(mean(WHO_FC_est), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(FVC_2yr))

cat("\n=== 2-Year Outcome Summary ===\n")
print(summary_2yr)

## ============================================================
## Virtual Patient Population (uncertainty analysis)
## ============================================================

set.seed(42)
n_vp <- 200  # 200 virtual patients

vp_params <- data.frame(
  ID        = 1:n_vp,
  FVC_0     = rnorm(n_vp, 78, 8),      # FVC at baseline
  mRSS_0    = rnorm(n_vp, 20, 6),      # mRSS at baseline
  PVR_0     = rnorm(n_vp, 350, 80),    # PVR at baseline
  kFVC_loss = rnorm(n_vp, 0.0012, 0.0003),
  EMAX_NINT_FVC = rnorm(n_vp, 0.55, 0.10)
)

# Clamp to physiological range
vp_params <- vp_params %>%
  mutate(
    FVC_0     = pmax(50, pmin(95, FVC_0)),
    mRSS_0    = pmax(5,  pmin(40, mRSS_0)),
    PVR_0     = pmax(150, pmin(800, PVR_0)),
    kFVC_loss = pmax(0.0005, pmin(0.003, kFVC_loss))
  )

run_vp_sim <- function(vp_row, model, events, label) {
  model %>%
    param(
      FVC_0        = vp_row$FVC_0,
      mRSS_0       = vp_row$mRSS_0,
      PVR_0        = vp_row$PVR_0,
      kFVC_loss    = vp_row$kFVC_loss,
      EMAX_NINT_FVC = vp_row$EMAX_NINT_FVC
    ) %>%
    init(FVC = vp_row$FVC_0, mRSS = vp_row$mRSS_0, PVR = vp_row$PVR_0) %>%
    ev(events) %>%
    mrgsim(end = SIM_DUR, delta = 24 * 30) %>%  # monthly
    as_tibble() %>%
    mutate(ID = vp_row$ID, scenario = label, time_yr = time / 8760)
}

cat("\nRunning virtual patient simulations...\n")
vp_none_list <- lapply(1:n_vp, function(i)
  run_vp_sim(vp_params[i, ], ssc_model, ev_none, "Untreated"))
vp_nint_list <- lapply(1:n_vp, function(i)
  run_vp_sim(vp_params[i, ], ssc_model, ev_nint, "Nintedanib"))

vp_none_df <- bind_rows(vp_none_list)
vp_nint_df <- bind_rows(vp_nint_list)
vp_all    <- bind_rows(vp_none_df, vp_nint_df)

# FVC distribution at 2 years
vp_2yr <- vp_all %>%
  filter(abs(time_yr - max(time_yr)) < 0.05)

p_vp <- ggplot(vp_2yr, aes(x = FVC, fill = scenario)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("Untreated" = "#e74c3c",
                               "Nintedanib" = "#2980b9")) +
  labs(title = "Virtual Patient FVC Distribution at 2 Years",
       x = "FVC % predicted", y = "Density",
       fill = "Treatment") +
  theme_bw()

print(p_vp)

cat("\n=== Virtual Patient FVC Summary at 2 Years ===\n")
vp_summary <- vp_2yr %>%
  group_by(scenario) %>%
  summarise(
    FVC_median = round(median(FVC), 1),
    FVC_q25    = round(quantile(FVC, 0.25), 1),
    FVC_q75    = round(quantile(FVC, 0.75), 1),
    prop_below70 = round(mean(FVC < 70) * 100, 1),
    .groups = "drop"
  )
print(vp_summary)

cat("\nSSc QSP model simulation complete.\n")
