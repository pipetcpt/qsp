## =============================================================================
## IBS QSP Model — mrgsolve ODE-based PK/PD Simulation
## Irritable Bowel Syndrome: Brain-Gut Axis, 5-HT Signaling,
##   Visceral Hypersensitivity, Mucosal Immunity, Drug PK/PD
##
## Compartments (23 ODEs):
##   Drug PK  (8): Alosetron (2-CMT oral), Linaclotide (luminal),
##                 Rifaximin (luminal), Loperamide (2-CMT oral)
##   5-HT     (3): EC-cell 5-HT release, SERT-mediated reuptake, 5-HT3R occupancy
##   Motility (3): Colonic transit index, intestinal secretion, MMC score
##   Pain/Sens(4): Visceral nociception, spinal sensitization,
##                 central sensitization, pain NRS
##   Immune   (3): Mast cell activity, low-grade inflammation, gut permeability
##   Disease  (2): Dysbiosis score, HPA stress index
##
## Treatment scenarios:
##   1. No treatment (natural history)
##   2. Alosetron 1 mg BID (IBS-D)
##   3. Linaclotide 290 µg QD (IBS-C)
##   4. Rifaximin 550 mg TID × 14 days (IBS-D/SIBO)
##   5. Loperamide 2 mg PRN (IBS-D, symptom rescue)
##   6. Alosetron + Low-dose TCA (combination, IBS-D refractory)
##   7. Linaclotide + Tegaserod 6 mg BID (combination, IBS-C refractory)
##
## Key calibration references:
##   Camilleri 2001 (Alosetron IBS-D, NEJM): pain ↓40%, urgency ↓50%
##   Chey 2012 (Linaclotide IBS-C, NEJM): SBM+1/wk responder 33%
##   Pimentel 2011 (Rifaximin IBS-D, NEJM): global relief 40.8% vs 31.2%
##   Chang 2010 (Loperamide): stool consistency BSS ↓1.2
##   Drossman 2003 (TCA/SSRI meta-analysis): IBS-SSS ↓50–70 pts
## =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ─────────────────────────────────────────────────────────────────────────────
# 1. MODEL DEFINITION
# ─────────────────────────────────────────────────────────────────────────────

ibs_model_code <- '
$PROB IBS QSP mrgsolve model v1.0

$PARAM
// ── Alosetron PK (2-CMT oral) ─────────────────────────────────────────────
  KA_ALO   = 0.70,   // h⁻¹  absorption rate (Tmax ~1 h)
  CL_ALO   = 50.0,   // L/h  total clearance (t½ ~1.5 h)
  VC_ALO   = 77.0,   // L    central volume
  VP_ALO   = 52.0,   // L    peripheral volume
  Q_ALO    = 18.0,   // L/h  inter-compartmental clearance
  F_ALO    = 0.60,   // bioavailability 60%
  IC50_ALO = 0.002,  // µg/mL  5-HT3R IC50 (1.6 nM × 294 g/mol)

// ── Linaclotide PK (luminal, minimal systemic) ────────────────────────────
  KA_LIN   = 0.10,   // h⁻¹  luminal absorption (very slow)
  CL_LIN   = 200.0,  // L/h  rapid systemic clearance (negligible plasma)
  VC_LIN   = 5.0,    // L    systemic volume (not relevant clinically)
  F_LIN    = 0.02,   // bioavailability <2%
  EC50_LIN = 0.05,   // nM   GC-C EC50

// ── Rifaximin PK (luminal, <0.4% systemic) ───────────────────────────────
  KA_RIF   = 0.05,   // h⁻¹  minimal GI absorption
  CL_RIF   = 300.0,  // L/h
  VC_RIF   = 10.0,   // L
  F_RIF    = 0.004,  // bioavailability <0.4%
  KLUM_RIF = 0.08,   // h⁻¹  luminal elimination (excretion)
  EC50_RIF = 0.10,   // µg/mL luminal MIC90 for gut flora

// ── Loperamide PK (2-CMT oral, BBB-excluded) ─────────────────────────────
  KA_LOP   = 0.30,   // h⁻¹  absorption (Tmax ~2.5 h)
  CL_LOP   = 14.0,   // L/h  (t½ ~10 h)
  VC_LOP   = 100.0,  // L
  VP_LOP   = 80.0,   // L
  Q_LOP    = 12.0,   // L/h
  F_LOP    = 0.40,   // bioavailability 40%
  IC50_LOP = 0.001,  // µg/mL µ-OR IC50

// ── TCA (Amitriptyline) PK (simplified 1-CMT) ────────────────────────────
  KA_TCA   = 0.35,   // h⁻¹
  CL_TCA   = 22.0,   // L/h  (t½ ~20 h)
  VC_TCA   = 1400.0, // L    large Vd (lipophilic)
  F_TCA    = 0.50,   // bioavailability 50%
  EC50_TCA = 0.08,   // µg/mL SERT inhibition

// ── Tegaserod PK (simplified 1-CMT) ──────────────────────────────────────
  KA_TEGA  = 0.55,   // h⁻¹
  CL_TEGA  = 77.0,   // L/h  (t½ ~11 h)
  VC_TEGA  = 368.0,  // L
  F_TEGA   = 0.10,   // bioavailability ~10%
  EC50_TEGA= 0.015,  // µg/mL 5-HT4R EC50

// ── 5-HT signaling ───────────────────────────────────────────────────────
  KSY_5HT  = 0.50,   // h⁻¹  5-HT synthesis/release rate constant
  KRU_5HT  = 1.20,   // h⁻¹  SERT-mediated reuptake rate
  KDG_5HT  = 0.30,   // h⁻¹  MAO degradation
  HT3R_BASE= 0.20,   // baseline 5-HT3R tone (0–1 normalized)
  HT4R_BASE= 0.30,   // baseline 5-HT4R tone

// ── Motility / secretion ─────────────────────────────────────────────────
  CTT_BASE = 36.0,   // h    baseline colonic transit time (normal)
  K_CTT    = 0.05,   // h⁻¹  motility return-to-baseline rate
  SEC_BASE = 1.0,    // arbitrary units, baseline secretion
  K_SEC    = 0.10,   // h⁻¹  secretion rate constant

// ── Visceral hypersensitivity & pain ─────────────────────────────────────
  P_BASE   = 5.5,    // baseline pain NRS (IBS typical, 0–10)
  K_PAIN   = 0.15,   // h⁻¹  pain return-to-baseline
  VH_BASE  = 0.60,   // baseline visceral hypersensitivity index (0–1)
  K_VH     = 0.008,  // h⁻¹  visceral hypersensitivity change rate
  K_CSS    = 0.004,  // h⁻¹  central sensitization rate

// ── Mucosal immunity ─────────────────────────────────────────────────────
  MC_BASE  = 0.40,   // baseline mast cell activity (0–1)
  K_MC     = 0.12,   // h⁻¹  MC activation rate
  LGI_BASE = 0.45,   // baseline low-grade inflammation (0–1)
  K_LGI    = 0.05,   // h⁻¹
  PERM_BASE= 0.50,   // baseline permeability index (0–1)
  K_PERM   = 0.03,   // h⁻¹

// ── Dysbiosis & HPA ──────────────────────────────────────────────────────
  DYS_BASE = 0.55,   // baseline dysbiosis score (0–1)
  K_DYS    = 0.02,   // h⁻¹  dysbiosis change rate
  STRESS_BASE = 0.60,// baseline HPA stress (0–1)
  K_STRESS = 0.04,   // h⁻¹

// ── IBS subtype flags ─────────────────────────────────────────────────────
  IBSD_FLAG = 1,     // 1=IBS-D, 0=IBS-C  (changes baseline CTT/secretion)
  IBS_SEVERITY = 2,  // 1=mild, 2=moderate, 3=severe (modifies baseline pain)

// ── Drug dosing flags ─────────────────────────────────────────────────────
  USE_ALO  = 0,      // alosetron on/off
  USE_LIN  = 0,      // linaclotide on/off
  USE_RIF  = 0,      // rifaximin on/off
  USE_LOP  = 0,      // loperamide on/off
  USE_TCA  = 0,      // TCA on/off
  USE_TEGA = 0       // tegaserod on/off

$CMT
// Drug PK compartments
  GUT_ALO, CENT_ALO, PERIPH_ALO,       // alosetron (oral depot + 2-CMT)
  LUM_LIN,  CENT_LIN,                   // linaclotide (luminal + systemic)
  LUM_RIF,  CENT_RIF,                   // rifaximin (luminal + systemic)
  GUT_LOP,  CENT_LOP, PERIPH_LOP,       // loperamide (oral + 2-CMT)
  CENT_TCA,                              // TCA (amitriptyline)
  CENT_TEGA,                             // tegaserod

// Serotonin signaling
  HT5_EC,    // EC-cell 5-HT pool (nmol/g tissue)
  HT3R_OCC,  // 5-HT3R fractional occupancy (0–1)
  HT4R_TON,  // 5-HT4R tone (0–1)

// Motility & secretion
  CTT,       // colonic transit time (h)
  SEC_IDX,   // intestinal secretion index (a.u.)
  MMC_IDX,   // migrating motor complex activity (0–1)

// Visceral hypersensitivity & pain
  VH_IDX,    // visceral hypersensitivity index (0–1)
  SP_SENS,   // spinal sensitization state (0–1)
  CENT_SENS, // central sensitization (0–1)
  PAIN_NRS,  // pain score (0–10)

// Mucosal immunity
  MC_ACT,    // mast cell activity (0–1)
  LGI_IDX,   // low-grade inflammation (0–1)
  PERM_IDX,  // gut permeability index (0–1)

// Disease state
  DYS_IDX,   // dysbiosis score (0–1)
  STRESS_IDX // HPA stress index (0–1)

$MAIN
// IBS subtype–specific adjustments at initialization
  double CTT_init    = IBSD_FLAG == 1 ? CTT_BASE * 0.72 : CTT_BASE * 1.35;
  double PAIN_init   = P_BASE + (IBS_SEVERITY - 1) * 1.5;
  double STRESS_init = STRESS_BASE + (IBS_SEVERITY - 1) * 0.10;

  if (NEWIND <= 1) {
    CTT_0      = CTT_init;
    PAIN_NRS_0 = PAIN_init;
    STRESS_IDX_0 = STRESS_init;
    SEC_IDX_0  = IBSD_FLAG == 1 ? SEC_BASE * 1.4 : SEC_BASE * 0.7;
    MC_ACT_0   = MC_BASE + (IBS_SEVERITY - 1) * 0.08;
    LGI_IDX_0  = LGI_BASE + (IBS_SEVERITY - 1) * 0.08;
    VH_IDX_0   = VH_BASE + (IBS_SEVERITY - 1) * 0.10;
    DYS_IDX_0  = DYS_BASE;
    PERM_IDX_0 = PERM_BASE;
    SP_SENS_0  = 0.45;
    CENT_SENS_0= 0.35;
    HT5_EC_0   = IBSD_FLAG == 1 ? 1.2 : 0.8;
    HT3R_OCC_0 = HT3R_BASE;
    HT4R_TON_0 = HT4R_BASE;
    MMC_IDX_0  = 0.60;
  }

$ODE
// ── Plasma concentrations (µg/mL) ────────────────────────────────────────
  double Cp_ALO  = CENT_ALO  / VC_ALO;
  double Cp_LIN  = CENT_LIN  / VC_LIN;
  double Cp_LUM_RIF = LUM_RIF;           // luminal concentration (a.u.)
  double Cp_LOP  = CENT_LOP  / VC_LOP;
  double Cp_TCA  = CENT_TCA  / VC_TCA;
  double Cp_TEGA = CENT_TEGA / VC_TEGA;

// ── PD effect functions (Hill, Emax model) ────────────────────────────────
  // Alosetron: 5-HT3R occupancy (0–1), Emax=1
  double E_ALO   = USE_ALO  * (Cp_ALO  / (Cp_ALO  + IC50_ALO));
  // Linaclotide: luminal GC-C activation
  double E_LIN   = USE_LIN  * (LUM_LIN / (LUM_LIN + EC50_LIN));
  // Rifaximin: luminal antibacterial effect (on dysbiosis)
  double E_RIF   = USE_RIF  * (Cp_LUM_RIF / (Cp_LUM_RIF + EC50_RIF));
  // Loperamide: µ-OR occupancy
  double E_LOP   = USE_LOP  * (Cp_LOP  / (Cp_LOP  + IC50_LOP));
  // TCA: SERT/NET inhibition
  double E_TCA   = USE_TCA  * (Cp_TCA  / (Cp_TCA  + EC50_TCA));
  // Tegaserod: 5-HT4R partial agonism
  double E_TEGA  = USE_TEGA * (Cp_TEGA / (Cp_TEGA + EC50_TEGA)) * 0.6; // partial: Emax=0.6

// ── Alosetron PK ─────────────────────────────────────────────────────────
  dxdt_GUT_ALO   = -KA_ALO * GUT_ALO;
  dxdt_CENT_ALO  = KA_ALO * F_ALO * GUT_ALO - (CL_ALO + Q_ALO) / VC_ALO * CENT_ALO
                   + Q_ALO / VP_ALO * PERIPH_ALO;
  dxdt_PERIPH_ALO= Q_ALO / VC_ALO * CENT_ALO - Q_ALO / VP_ALO * PERIPH_ALO;

// ── Linaclotide PK ────────────────────────────────────────────────────────
  dxdt_LUM_LIN   = -KA_LIN * LUM_LIN - 0.30 * LUM_LIN; // luminal → absorbed + degraded
  dxdt_CENT_LIN  = KA_LIN * F_LIN * LUM_LIN - CL_LIN / VC_LIN * CENT_LIN;

// ── Rifaximin PK ──────────────────────────────────────────────────────────
  dxdt_LUM_RIF   = -KLUM_RIF * LUM_RIF - KA_RIF * LUM_RIF;
  dxdt_CENT_RIF  = KA_RIF * F_RIF * LUM_RIF - CL_RIF / VC_RIF * CENT_RIF;

// ── Loperamide PK ─────────────────────────────────────────────────────────
  dxdt_GUT_LOP   = -KA_LOP * GUT_LOP;
  dxdt_CENT_LOP  = KA_LOP * F_LOP * GUT_LOP - (CL_LOP + Q_LOP) / VC_LOP * CENT_LOP
                   + Q_LOP / VP_LOP * PERIPH_LOP;
  dxdt_PERIPH_LOP= Q_LOP / VC_LOP * CENT_LOP - Q_LOP / VP_LOP * PERIPH_LOP;

// ── TCA PK ────────────────────────────────────────────────────────────────
  dxdt_CENT_TCA  = KA_TCA * F_TCA * 0 - CL_TCA / VC_TCA * CENT_TCA;
  // Note: TCA dose injected via event table

// ── Tegaserod PK ──────────────────────────────────────────────────────────
  dxdt_CENT_TEGA = KA_TEGA * F_TEGA * 0 - CL_TEGA / VC_TEGA * CENT_TEGA;

// ── EC-cell 5-HT pool ────────────────────────────────────────────────────
  // LGI drives TPH1 upregulation → more 5-HT; SERT clearance
  double TPH1_up  = 1.0 + 0.6 * LGI_IDX;
  double SERT_fun = 1.0 - 0.85 * E_TCA;   // TCA blocks SERT
  dxdt_HT5_EC = KSY_5HT * TPH1_up - KRU_5HT * SERT_fun * HT5_EC - KDG_5HT * HT5_EC;

// ── 5-HT3R occupancy ─────────────────────────────────────────────────────
  // HT5_EC drives 5-HT3R activation; alosetron blocks it
  double HT3R_drive = HT5_EC / (HT5_EC + 1.0);   // normalized
  dxdt_HT3R_OCC = 0.30 * (HT3R_drive - HT3R_OCC) - E_ALO * HT3R_OCC * 1.5;

// ── 5-HT4R tone ──────────────────────────────────────────────────────────
  // Baseline 5-HT4R; tegaserod partial agonism increases it
  double HT4R_drive = HT5_EC / (HT5_EC + 1.5);
  dxdt_HT4R_TON = 0.20 * (HT4R_drive - HT4R_TON) + E_TEGA * (1.0 - HT4R_TON) * 0.5;

// ── Colonic transit time ─────────────────────────────────────────────────
  // 5-HT3R ↑ → faster transit (IBS-D); alosetron blocks → normalises
  // 5-HT4R ↑ → faster (prokinetic); tegaserod drives IBS-C
  // Loperamide → slows; linaclotide → faster (IBS-C)
  double CTT_target = CTT_init
                     * (1.0 - 0.40 * E_ALO)          // alosetron: slow
                     * (1.0 - 0.25 * E_LOP)           // loperamide: slow
                     * (1.0 - 0.20 * E_TEGA * (1.0 - IBSD_FLAG)) // tega: IBS-C only
                     * (1.0 - 0.10 * E_TCA);           // TCA: slow (anticholinergic)
  // Linaclotide: speeds transit in IBS-C
  double CTT_lina_eff = USE_LIN * IBSD_FLAG == 0 ? -0.15 * E_LIN * CTT : 0;
  dxdt_CTT = K_CTT * (CTT_target - CTT) + CTT_lina_eff;

// ── Secretion index ───────────────────────────────────────────────────────
  // Linaclotide/lubiprostone → ↑ secretion (IBS-C benefit)
  // Loperamide → ↓ secretion
  double SEC_target = (IBSD_FLAG == 1 ? SEC_BASE * 1.4 : SEC_BASE * 0.7)
                     * (1.0 + 0.50 * E_LIN * (1 - IBSD_FLAG)) // lina → ↑ sec
                     * (1.0 - 0.35 * E_LOP);                   // lope → ↓ sec
  dxdt_SEC_IDX = K_SEC * (SEC_target - SEC_IDX);

// ── MMC activity ──────────────────────────────────────────────────────────
  // MMC impaired in IBS-D (rapid); rifaximin restores by clearing SIBO
  double MMC_target = 0.60 + 0.25 * E_RIF; // rifaximin restores MMC
  dxdt_MMC_IDX = 0.04 * (MMC_target - MMC_IDX);

// ── Visceral hypersensitivity ─────────────────────────────────────────────
  // 5-HT3R drives sensitization; alosetron reduces; TCA reduces pain signal
  // Linaclotide: cGMP → ↓ TRPV1 activation (reduces pain afferent firing)
  double VH_target = VH_BASE + 0.20 * HT3R_OCC
                   - 0.30 * E_ALO
                   - 0.15 * E_TCA
                   - 0.12 * E_LIN;
  VH_target = fmax(0.10, fmin(1.0, VH_target));
  dxdt_VH_IDX = K_VH * (VH_target - VH_IDX);

// ── Spinal sensitization ──────────────────────────────────────────────────
  // VH drives spinal sensitization (wind-up); LGI (cytokines) amplifies
  double SP_target = 0.40 * VH_IDX + 0.20 * LGI_IDX + 0.10 * STRESS_IDX
                   - 0.20 * E_TCA; // TCA centrally modulates pain
  SP_target = fmax(0.05, fmin(0.95, SP_target));
  dxdt_SP_SENS = K_CSS * (SP_target - SP_SENS);

// ── Central sensitization ─────────────────────────────────────────────────
  double CS_target = 0.50 * SP_SENS + 0.25 * STRESS_IDX
                   - 0.25 * E_TCA;
  CS_target = fmax(0.05, fmin(0.95, CS_target));
  dxdt_CENT_SENS = K_CSS * 0.8 * (CS_target - CENT_SENS);

// ── Pain NRS ──────────────────────────────────────────────────────────────
  // Driven by VH, spinal & central sensitization; modulated by drugs
  double P_target = P_BASE
                   + 3.0 * VH_IDX
                   + 2.5 * SP_SENS
                   + 2.0 * CENT_SENS
                   + 1.5 * MC_ACT
                   + (IBS_SEVERITY - 1) * 1.5
                   - 3.5 * E_ALO          // alosetron: pain ↓40%
                   - 2.0 * E_TCA          // TCA: centrally
                   - 1.5 * E_LIN;         // linaclotide: visceral afferent
  P_target = fmax(0.0, fmin(10.0, P_target));
  dxdt_PAIN_NRS = K_PAIN * (P_target - PAIN_NRS);

// ── Mast cell activity ────────────────────────────────────────────────────
  // Stress (HPA/CRH) drives MC; LGI amplifies; drugs modulate weakly
  double MC_target = MC_BASE + 0.35 * STRESS_IDX + 0.20 * DYS_IDX
                   - 0.10 * E_RIF;  // rifaximin reduces dysbiosis → ↓ MC
  MC_target = fmax(0.10, fmin(1.0, MC_target));
  dxdt_MC_ACT = K_MC * (MC_target - MC_ACT);

// ── Low-grade inflammation ────────────────────────────────────────────────
  // MC + dysbiosis + permeability → LGI; butyrate (inverse dysbiosis) → ↓
  double LGI_target = LGI_BASE + 0.30 * MC_ACT + 0.25 * DYS_IDX
                    + 0.15 * PERM_IDX
                    - 0.15 * E_RIF;
  LGI_target = fmax(0.10, fmin(1.0, LGI_target));
  dxdt_LGI_IDX = K_LGI * (LGI_target - LGI_IDX);

// ── Gut permeability ──────────────────────────────────────────────────────
  // LGI + stress → leaky gut; rifaximin (↓ dysbiosis) → repairs
  double PERM_target = PERM_BASE + 0.25 * LGI_IDX + 0.20 * STRESS_IDX
                     - 0.15 * E_RIF;
  PERM_target = fmax(0.10, fmin(1.0, PERM_target));
  dxdt_PERM_IDX = K_PERM * (PERM_target - PERM_IDX);

// ── Dysbiosis ─────────────────────────────────────────────────────────────
  // Rifaximin reduces; stress chronically shifts; microbiome slow dynamics
  double DYS_target = DYS_BASE + 0.20 * STRESS_IDX - 0.45 * E_RIF;
  DYS_target = fmax(0.10, fmin(1.0, DYS_target));
  dxdt_DYS_IDX = K_DYS * (DYS_target - DYS_IDX);

// ── HPA stress index ──────────────────────────────────────────────────────
  // Exogenous stress input (chronic); TCA modestly dampens HPA
  double STRESS_target = STRESS_BASE - 0.10 * E_TCA;
  dxdt_STRESS_IDX = K_STRESS * (STRESS_target - STRESS_IDX);

$TABLE
  capture Cp_ALO  = CENT_ALO  / VC_ALO;
  capture Cp_LOP  = CENT_LOP  / VC_LOP;
  capture Cp_TCA  = CENT_TCA  / VC_TCA;
  capture Cp_TEGA = CENT_TEGA / VC_TEGA;
  capture LumRIF  = LUM_RIF;
  capture LumLIN  = LUM_LIN;
  capture serotonin = HT5_EC;
  capture pain    = PAIN_NRS;
  capture ctransit= CTT;
  capture secretion = SEC_IDX;
  capture mmc     = MMC_IDX;
  capture viscHyp = VH_IDX;
  capture spinal  = SP_SENS;
  capture central = CENT_SENS;
  capture mast    = MC_ACT;
  capture inflam  = LGI_IDX;
  capture permeab = PERM_IDX;
  capture dysbiosis = DYS_IDX;
  capture stress  = STRESS_IDX;
  capture ht3r    = HT3R_OCC;
  capture ht4r    = HT4R_TON;
  capture ibs_sss = fmin(500, (pain / 10.0) * 100 + (viscHyp * 150)
                         + (fabs(ctransit - 36.0) / 36.0) * 100
                         + (inflam * 75) + (dysbiosis * 75));
'

# Compile model
mod <- mcode("ibs_qsp", ibs_model_code, quiet = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# 2. DOSING HELPERS
# ─────────────────────────────────────────────────────────────────────────────

# Alosetron 1 mg BID × 12 weeks
dose_alo <- function(dose_mg = 1, freq_h = 12, duration_d = 84) {
  ev(amt   = dose_mg * 1000 / 294.35 * VC_ALO_val,  # rough dose in model units
     cmt   = "GUT_ALO",
     ii    = freq_h, addl = as.integer(duration_d * 24 / freq_h) - 1)
}

# Simpler: direct amount in mg/L equivalents
make_ev_mg <- function(amt_mg, cmt_name, freq_h = 12, duration_d = 84) {
  ev(amt  = amt_mg, cmt = cmt_name,
     ii   = freq_h, addl = as.integer(duration_d * 24 / freq_h) - 1)
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. SCENARIO DEFINITIONS
# ─────────────────────────────────────────────────────────────────────────────

scenario_list <- list(
  # 1 — No treatment (IBS-D, moderate)
  list(
    name = "1. No Treatment (IBS-D)", col = "#aaaaaa",
    p    = list(IBSD_FLAG=1, IBS_SEVERITY=2,
                USE_ALO=0, USE_LIN=0, USE_RIF=0,
                USE_LOP=0, USE_TCA=0, USE_TEGA=0),
    ev   = ev(amt=0, cmt="GUT_ALO", time=0)
  ),
  # 2 — Alosetron 1 mg BID (IBS-D)
  list(
    name = "2. Alosetron 1 mg BID (IBS-D)", col = "#2a9aff",
    p    = list(IBSD_FLAG=1, IBS_SEVERITY=2,
                USE_ALO=1, USE_LIN=0, USE_RIF=0,
                USE_LOP=0, USE_TCA=0, USE_TEGA=0),
    ev   = ev(amt=1, cmt="GUT_ALO", ii=12,
              addl=13*7*2-1) # 13 weeks BID
  ),
  # 3 — Linaclotide 290 µg QD (IBS-C)
  list(
    name = "3. Linaclotide 290 µg QD (IBS-C)", col = "#ff8844",
    p    = list(IBSD_FLAG=0, IBS_SEVERITY=2,
                USE_ALO=0, USE_LIN=1, USE_RIF=0,
                USE_LOP=0, USE_TCA=0, USE_TEGA=0),
    ev   = ev(amt=0.29, cmt="LUM_LIN", ii=24,
              addl=13*7-1)
  ),
  # 4 — Rifaximin 550 mg TID × 14 d (IBS-D/SIBO)
  list(
    name = "4. Rifaximin 550 mg TID ×14 d (IBS-D)", col = "#44cc44",
    p    = list(IBSD_FLAG=1, IBS_SEVERITY=2,
                USE_ALO=0, USE_LIN=0, USE_RIF=1,
                USE_LOP=0, USE_TCA=0, USE_TEGA=0),
    ev   = ev(amt=550, cmt="LUM_RIF", ii=8, addl=3*14-1)
  ),
  # 5 — Loperamide 2 mg PRN BID (IBS-D)
  list(
    name = "5. Loperamide 2 mg BID (IBS-D)", col = "#cc44cc",
    p    = list(IBSD_FLAG=1, IBS_SEVERITY=2,
                USE_ALO=0, USE_LIN=0, USE_RIF=0,
                USE_LOP=1, USE_TCA=0, USE_TEGA=0),
    ev   = ev(amt=2, cmt="GUT_LOP", ii=12, addl=13*7*2-1)
  ),
  # 6 — Alosetron + TCA combination (IBS-D refractory)
  list(
    name = "6. Alosetron + TCA 25 mg (IBS-D refractory)", col = "#ff4444",
    p    = list(IBSD_FLAG=1, IBS_SEVERITY=3,
                USE_ALO=1, USE_LIN=0, USE_RIF=0,
                USE_LOP=0, USE_TCA=1, USE_TEGA=0),
    ev   = c(ev(amt=1, cmt="GUT_ALO", ii=12, addl=13*7*2-1),
             ev(amt=25, cmt="CENT_TCA", ii=24, addl=13*7-1))
  ),
  # 7 — Linaclotide + Tegaserod (IBS-C refractory)
  list(
    name = "7. Linaclotide + Tegaserod 6 mg BID (IBS-C refractory)", col = "#ffcc00",
    p    = list(IBSD_FLAG=0, IBS_SEVERITY=3,
                USE_ALO=0, USE_LIN=1, USE_RIF=0,
                USE_LOP=0, USE_TCA=0, USE_TEGA=1),
    ev   = c(ev(amt=0.29, cmt="LUM_LIN", ii=24, addl=13*7-1),
             ev(amt=6, cmt="CENT_TEGA", ii=12, addl=13*7*2-1))
  )
)

# ─────────────────────────────────────────────────────────────────────────────
# 4. RUN SIMULATIONS
# ─────────────────────────────────────────────────────────────────────────────

sim_time <- seq(0, 91 * 24, by = 6)   # 91 days (13 weeks) at 6-h resolution

run_scenario <- function(sc) {
  p_update <- do.call(param, c(list(mod), sc$p))
  out <- mrgsim(p_update, events = sc$ev, tgrid = sim_time, obsonly = TRUE)
  as.data.frame(out) %>%
    mutate(scenario = sc$name, color = sc$col,
           time_d = time / 24,
           time_wk = time / 168)
}

cat("Running IBS QSP simulations...\n")
results_all <- bind_rows(lapply(scenario_list, run_scenario))
cat("Done. Rows:", nrow(results_all), "\n")

# Weekly summary
results_wk <- results_all %>%
  mutate(week = floor(time_wk)) %>%
  group_by(scenario, color, week) %>%
  summarise(
    pain_mean     = mean(pain, na.rm = TRUE),
    ctransit_mean = mean(ctransit, na.rm = TRUE),
    inflam_mean   = mean(inflam, na.rm = TRUE),
    dysbiosis_mean= mean(dysbiosis, na.rm = TRUE),
    permeab_mean  = mean(permeab, na.rm = TRUE),
    mast_mean     = mean(mast, na.rm = TRUE),
    serotonin_mean= mean(serotonin, na.rm = TRUE),
    ht3r_mean     = mean(ht3r, na.rm = TRUE),
    ht4r_mean     = mean(ht4r, na.rm = TRUE),
    ibs_sss_mean  = mean(ibs_sss, na.rm = TRUE),
    .groups = "drop"
  )

# ─────────────────────────────────────────────────────────────────────────────
# 5. KEY PLOTS
# ─────────────────────────────────────────────────────────────────────────────

theme_qsp <- theme_bw(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "#0a0a14", color = NA),
    panel.background = element_rect(fill = "#121230"),
    panel.grid       = element_line(color = "#2a2a50"),
    axis.text        = element_text(color = "white"),
    axis.title       = element_text(color = "white"),
    plot.title       = element_text(color = "white", face = "bold"),
    plot.subtitle    = element_text(color = "#aaaacc"),
    legend.background= element_rect(fill = "#0a0a20"),
    legend.text      = element_text(color = "white"),
    legend.title     = element_text(color = "white"),
    strip.background = element_rect(fill = "#1a1a40"),
    strip.text       = element_text(color = "white")
  )

sc_colors <- setNames(
  sapply(scenario_list, `[[`, "col"),
  sapply(scenario_list, `[[`, "name")
)

# Plot 1: Abdominal Pain NRS over time
p1 <- ggplot(results_wk, aes(x = week, y = pain_mean,
                              color = scenario, group = scenario)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 3.5, linetype = "dashed", color = "#ffcc00", size = 0.8) +
  annotate("text", x = 2, y = 3.7, label = "Mild–Moderate threshold",
           color = "#ffcc00", size = 3.5) +
  scale_color_manual(values = sc_colors) +
  labs(title = "IBS QSP — Abdominal Pain NRS (0–10)",
       subtitle = "Weekly mean, all 7 treatment scenarios",
       x = "Week", y = "Pain NRS", color = "Scenario") +
  scale_x_continuous(breaks = seq(0, 13, 2)) +
  theme_qsp

# Plot 2: IBS-SSS composite score
p2 <- ggplot(results_wk, aes(x = week, y = ibs_sss_mean,
                               color = scenario, group = scenario)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 175, linetype = "dashed", color = "#44ff44") +
  annotate("text", x = 1.5, y = 185, label = "Mild (<175)",
           color = "#44ff44", size = 3.5) +
  scale_color_manual(values = sc_colors) +
  labs(title = "IBS-SSS Composite Score (0–500)",
       subtitle = "Severe >300, Moderate 175–300, Mild <175",
       x = "Week", y = "IBS-SSS", color = "Scenario") +
  theme_qsp

# Plot 3: Colonic transit time
p3 <- ggplot(results_wk, aes(x = week, y = ctransit_mean,
                               color = scenario, group = scenario)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = 30, ymax = 42), fill = "#44aa44", alpha = 0.1) +
  scale_color_manual(values = sc_colors) +
  labs(title = "Colonic Transit Time (h)",
       subtitle = "Normal range: 30–42 h (green band)",
       x = "Week", y = "CTT (h)", color = "Scenario") +
  theme_qsp

# Plot 4: Gut microbiome dysbiosis & low-grade inflammation
p4 <- results_wk %>%
  select(scenario, color, week, dysbiosis_mean, inflam_mean) %>%
  pivot_longer(cols = c(dysbiosis_mean, inflam_mean),
               names_to = "biomarker", values_to = "value") %>%
  mutate(biomarker = ifelse(biomarker == "dysbiosis_mean",
                            "Dysbiosis Score", "Low-Grade Inflammation")) %>%
  ggplot(aes(x = week, y = value, color = scenario, linetype = biomarker)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = sc_colors) +
  scale_linetype_manual(values = c("Dysbiosis Score" = "solid",
                                   "Low-Grade Inflammation" = "dashed")) +
  labs(title = "Gut Dysbiosis & Low-Grade Inflammation (0–1 index)",
       x = "Week", y = "Index (0–1)", color = "Scenario",
       linetype = "Biomarker") +
  theme_qsp

# Plot 5: Serotonin & 5-HT3R occupancy
p5 <- results_wk %>%
  filter(scenario %in% c("1. No Treatment (IBS-D)",
                          "2. Alosetron 1 mg BID (IBS-D)",
                          "6. Alosetron + TCA 25 mg (IBS-D refractory)")) %>%
  select(scenario, color, week, serotonin_mean, ht3r_mean) %>%
  pivot_longer(cols = c(serotonin_mean, ht3r_mean)) %>%
  mutate(name = ifelse(name == "serotonin_mean", "EC-cell 5-HT (a.u.)", "5-HT3R Occupancy")) %>%
  ggplot(aes(x = week, y = value, color = scenario, linetype = name)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = sc_colors) +
  scale_linetype_manual(values = c("EC-cell 5-HT (a.u.)" = "solid",
                                   "5-HT3R Occupancy" = "dashed")) +
  labs(title = "5-HT Signaling: EC-cell Serotonin & 5-HT3R Occupancy",
       subtitle = "IBS-D scenarios (alosetron effect)",
       x = "Week", y = "Value (a.u.)", color = "Scenario",
       linetype = "Biomarker") +
  theme_qsp

# Plot 6: Week 13 summary bar chart
wk13_summary <- results_wk %>%
  filter(week == 13) %>%
  select(scenario, pain_mean, ibs_sss_mean, ctransit_mean,
         inflam_mean, dysbiosis_mean) %>%
  mutate(
    pain_pct_change = 100 * (pain_mean - results_wk %>%
                               filter(week==0, scenario==scenario[1]) %>%
                               pull(pain_mean) %>% mean()),
    pain_label = round(pain_mean, 1)
  )

# Print clinical endpoint table
cat("\n=== 13-Week Clinical Endpoint Summary ===\n")
endpoint_table <- results_wk %>%
  filter(week %in% c(0, 4, 8, 13)) %>%
  select(scenario, week, pain_mean, ctransit_mean,
         inflam_mean, ibs_sss_mean) %>%
  mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  rename(
    `Week`    = week,
    `Scenario`= scenario,
    `Pain NRS`= pain_mean,
    `CTT (h)` = ctransit_mean,
    `Inflammation` = inflam_mean,
    `IBS-SSS` = ibs_sss_mean
  )
print(endpoint_table, n = Inf)

# ─────────────────────────────────────────────────────────────────────────────
# 6. SENSITIVITY ANALYSIS: Stress level vs Pain/SSS
# ─────────────────────────────────────────────────────────────────────────────

stress_levels <- seq(0.2, 0.9, by = 0.1)
sens_results <- lapply(stress_levels, function(sl) {
  p_s <- param(mod, IBSD_FLAG=1, IBS_SEVERITY=2,
               USE_ALO=0, USE_LIN=0, USE_RIF=0,
               USE_LOP=0, USE_TCA=0, USE_TEGA=0,
               STRESS_BASE=sl)
  out <- mrgsim(p_s,
                events = ev(amt=0, cmt="GUT_ALO"),
                tgrid = seq(0, 13*168, 168),
                obsonly = TRUE)
  df <- as.data.frame(out)
  wk13 <- tail(df, 1)
  data.frame(stress = sl,
             pain   = wk13$pain,
             ibs_sss= wk13$ibs_sss,
             inflam = wk13$inflam,
             dysbiosis = wk13$dysbiosis)
})
sens_df <- bind_rows(sens_results)

p_sens <- ggplot(sens_df) +
  geom_line(aes(x = stress, y = pain, color = "Pain NRS"), size = 1.3) +
  geom_line(aes(x = stress, y = ibs_sss / 50, color = "IBS-SSS / 50"), size = 1.3) +
  geom_line(aes(x = stress, y = inflam * 10, color = "Inflammation × 10"), size = 1.3, linetype = "dashed") +
  scale_color_manual(values = c("Pain NRS" = "#ff4444",
                                 "IBS-SSS / 50" = "#ffcc00",
                                 "Inflammation × 10" = "#44ccff")) +
  labs(title = "Sensitivity: HPA Stress Level vs IBS Endpoints (Week 13)",
       x = "Baseline Stress Index (0–1)",
       y = "Scaled outcome", color = "Endpoint") +
  theme_qsp

# ─────────────────────────────────────────────────────────────────────────────
# 7. DOSE–RESPONSE ANALYSIS (Alosetron)
# ─────────────────────────────────────────────────────────────────────────────

alo_doses <- c(0.25, 0.5, 1.0, 2.0)
dose_resp <- lapply(alo_doses, function(d) {
  ev_d <- ev(amt = d, cmt = "GUT_ALO", ii = 12, addl = 13*7*2-1)
  p_d  <- param(mod, IBSD_FLAG=1, IBS_SEVERITY=2,
                USE_ALO=1, USE_LIN=0, USE_RIF=0,
                USE_LOP=0, USE_TCA=0, USE_TEGA=0)
  out  <- mrgsim(p_d, events = ev_d,
                 tgrid = seq(0, 13*168, 168), obsonly = TRUE)
  df   <- as.data.frame(out)
  wk13 <- tail(df, 1)
  data.frame(dose = d,
             pain = wk13$pain,
             ctransit = wk13$ctransit,
             ht3r_occ = wk13$ht3r,
             ibs_sss  = wk13$ibs_sss)
})
dose_resp_df <- bind_rows(dose_resp)

cat("\n=== Alosetron Dose–Response at Week 13 ===\n")
print(dose_resp_df)

p_dr <- ggplot(dose_resp_df) +
  geom_line(aes(x = dose, y = pain, color = "Pain NRS"), size = 1.4) +
  geom_point(aes(x = dose, y = pain, color = "Pain NRS"), size = 3) +
  geom_line(aes(x = dose, y = ht3r_occ * 10, color = "5-HT3R Occ. × 10"), size = 1.4) +
  geom_point(aes(x = dose, y = ht3r_occ * 10, color = "5-HT3R Occ. × 10"), size = 3) +
  scale_x_log10() +
  scale_color_manual(values = c("Pain NRS" = "#ff4444",
                                 "5-HT3R Occ. × 10" = "#2a9aff")) +
  labs(title = "Alosetron Dose–Response (IBS-D, Week 13)",
       x = "Alosetron Dose per BID dose (mg, log scale)",
       y = "Endpoint value", color = "Endpoint") +
  theme_qsp

# ─────────────────────────────────────────────────────────────────────────────
# 8. SAVE RESULTS
# ─────────────────────────────────────────────────────────────────────────────

cat("\n=== IBS QSP Model: Key Model Parameters ===\n")
cat("  ODE compartments: 23\n")
cat("  Drug classes modeled: 6 (alosetron, linaclotide, rifaximin,\n")
cat("                           loperamide, TCA, tegaserod)\n")
cat("  Treatment scenarios: 7\n")
cat("  Simulation duration: 13 weeks\n")
cat("  Sensitivity analysis: HPA stress 0.2–0.9\n")
cat("  Dose-response: alosetron 0.25–2 mg BID\n")
cat("\n  Key disease drivers:\n")
cat("    Brain-Gut Axis: HPA stress → CRH → MC activation → LGI → TJ disruption\n")
cat("    5-HT axis: EC cell release → 5-HT3R (pain/transit) + 5-HT4R (prokinetic)\n")
cat("    Dysbiosis: ↑LPS → TLR4 → NFκB → MC degranulation → pain sensitization\n")
cat("    VH: TRPV1 sensitization → spinal → central sensitization → pain\n")

# Return all key objects for interactive use
invisible(list(
  mod           = mod,
  results_all   = results_all,
  results_wk    = results_wk,
  sens_df       = sens_df,
  dose_resp_df  = dose_resp_df,
  plots = list(pain = p1, sss = p2, ctt = p3,
               microbiome = p4, serotonin = p5,
               sensitivity = p_sens, doseresponse = p_dr)
))
