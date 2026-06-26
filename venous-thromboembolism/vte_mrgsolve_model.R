##############################################################################
## Venous Thromboembolism (VTE) QSP Model — mrgsolve
## Coagulation Cascade · Fibrinolysis · Multi-Drug PK/PD
##
## Compartments (19 ODEs):
##   PK  : RIV_GUT, RIV_CENT, RIV_PERIPH (Rivaroxaban)
##          APIX_CENT (Apixaban)
##          DABI_CENT (Dabigatran)
##          WARF_CENT (Warfarin)
##          ENOX_CENT (Enoxaparin anti-Xa)
##   PD  : FXa_FREE, FIIa_FREE (coagulation activity)
##          FIBRIN (fibrin accumulation)
##          CLOT (thrombus burden)
##          PLASMIN (fibrinolysis)
##          DDIMER (biomarker)
##          VK_OX, VK_RED (Vitamin K cycle)
##          FVII_pool, FX_pool, FII_pool (factor synthesis/depletion)
##
## Key References:
##   - Mueck W et al. J Thromb Haemost 2011 (Rivaroxaban PK/PD)
##   - Frost C et al. Clin Pharmacokinet 2015 (Apixaban PK)
##   - Stangier J et al. Clin Pharmacokinet 2008 (Dabigatran PK)
##   - Holford NHG. Clin Pharmacokinet 1986 (Warfarin model)
##   - Walenga JM et al. Semin Thromb Hemost 1999 (LMWH)
##   - Leidel BA et al. Eur J Emerg Med 2010 (VTE treatment)
##############################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ══════════════════════════════════════════════════════════════════════════
## MODEL CODE
## ══════════════════════════════════════════════════════════════════════════

vte_model_code <- '
$PROB
VTE QSP Model — Coagulation/Fibrinolysis/Multi-Drug PK-PD
Rivaroxaban | Apixaban | Dabigatran | Warfarin | Enoxaparin

$PARAM
// ── Rivaroxaban PK (2-compartment, oral) ────────────────────────────────
// Mueck W et al. Clin Pharmacokinet 2011; Kubitza D et al.
KA_RIV   = 1.2    // h-1 absorption rate constant (fed state)
F_RIV    = 0.93   // bioavailability (with food)
CL_RIV   = 4.8    // L/h total clearance (renal+metabolic)
V1_RIV   = 33.0   // L central volume of distribution
Q_RIV    = 3.2    // L/h intercompartmental clearance
V2_RIV   = 20.5   // L peripheral volume
WT_ref   = 70     // kg reference weight

// ── Apixaban PK (1-compartment) ─────────────────────────────────────────
// Frost C et al. Clin Pharmacokinet 2015
KA_APIX  = 0.78   // h-1
F_APIX   = 0.50   // ~50% bioavailability
CL_APIX  = 3.3    // L/h
V1_APIX  = 23.0   // L

// ── Dabigatran PK (1-compartment, prodrug) ──────────────────────────────
// Stangier J et al. Clin Pharmacokinet 2008
KA_DABI  = 0.35   // h-1 (slow absorption, acid-labile)
F_DABI   = 0.065  // ~6.5% bioavailability (pH dependent)
CL_DABI  = 8.5    // L/h (primarily renal)
V1_DABI  = 80.0   // L (high Vd)

// ── Warfarin PK + VK cycle (indirect) ───────────────────────────────────
// Holford 1986; Goulooze SC Clin Pharmacokinet 2020
KA_WARF  = 0.9    // h-1 (rapid oral absorption)
F_WARF   = 1.0    // ~100% bioavailability
CL_WARF  = 0.20   // L/h (CYP2C9-mediated)
V1_WARF  = 9.5    // L
IC50_WARF= 0.65   // mg/L (EC50 for VKORC1 inhibition)
HILL_W   = 1.0    // Hill coefficient (warfarin-VKORC1)
KIN_VK   = 0.14   // 1/h VitK oxidation to reduced form rate
KOUT_VK  = 0.18   // 1/h VitK reduced → oxidized (baseline)
VK0_ox   = 1.0    // relative units baseline VitKox
VK0_red  = VK0_ox * KIN_VK / KOUT_VK  // steady-state ratio

// ── Enoxaparin PK (SC, anti-Xa) ─────────────────────────────────────────
// Hulot JS et al. Clin Pharmacol Ther 2004; Becker RC et al.
KA_ENOX  = 0.23   // h-1 SC absorption rate (Tmax ~3-5h)
F_ENOX   = 0.92   // SC bioavailability ~92%
CL_ENOX  = 0.82   // L/h (renal + saturable non-renal)
V1_ENOX  = 4.5    // L (low Vd — stays in plasma)

// ── Coagulation PD Parameters ────────────────────────────────────────────
// Thrombin generation model adapted from:
// Wajima T et al. Clin Pharmacol Ther 2009 (mechanism-based coag model)
FXa_base  = 1.0   // nM baseline FXa activity (relative units)
KF_FXa   = 0.12   // h-1 FXa clearance (inactivation by AT etc.)
KG_FXa   = 0.08   // h-1 FXa generation rate at baseline
EMAX_RIV = 0.97   // maximum inhibition by Rivaroxaban
EC50_RIV = 12.0   // ng/mL (EC50 for FXa inhibition)
HILL_RIV = 1.3    // Hill coefficient
EMAX_APIX= 0.97   // maximum inhibition by Apixaban
EC50_APIX= 5.0    // ng/mL
HILL_APX = 1.2
EMAX_ENOX= 0.85   // Anti-Xa via LMWH-AT
EC50_ENOX= 0.35   // IU/mL

FIIa_base = 1.0   // nM baseline thrombin (relative)
KF_FIIa  = 0.20   // h-1 FIIa clearance
KG_FIIa  = 0.15   // h-1 FIIa generation (proportional to FXa)
EMAX_DABI= 0.95   // direct thrombin inhibitor Emax
EC50_DABI= 35.0   // ng/mL dabigatran EC50
HILL_DABI= 1.0

// ── Fibrin Dynamics ─────────────────────────────────────────────────────
FIBRIN_base=0.0   // initial clot = 0
KF_FBR   = 0.05   // h-1 fibrin formation rate (proportional to FIIa)
KL_FBR   = 0.03   // h-1 fibrin spontaneous lysis (plasmin-dependent)
CLOT_max = 100.0  // maximum clot size (arbitrary units)

// ── Clot (Thrombus) Dynamics ────────────────────────────────────────────
KG_CLOT  = 0.04   // h-1 clot growth rate
KD_CLOT  = 0.008  // h-1 spontaneous clot decay
CLOT_init= 80.0   // initial clot burden (DVT/PE scenario)

// ── Fibrinolysis (Plasmin) ───────────────────────────────────────────────
// tPA-plasminogen-plasmin system
PLASMIN_base = 0.5// relative plasmin activity
KP_FORM  = 0.06   // h-1 plasmin formation rate
KP_DECAY = 0.15   // h-1 plasmin decay (alpha2-antiplasmin)
PAI1_eff = 0.7    // PAI-1 inhibitory effect (fraction)
K_DDIMER = 0.10   // h-1 D-dimer generation from fibrinolysis
K_DDIMER_CL=0.025 // h-1 D-dimer clearance

// ── Factor Synthesis (Warfarin PD — indirect response) ──────────────────
KSYN_FII  = 0.012 // h-1 FII (t1/2 ~58h)
KDEG_FII  = 0.012 // h-1 FII degradation
FII_init  = 100.0 // % normal FII at baseline
KSYN_FVII = 0.12  // h-1 FVII (t1/2 ~5.7h, earliest drop)
KDEG_FVII = 0.12  // h-1 FVII degradation
FVII_init = 100.0 // % normal FVII
KSYN_FX   = 0.017 // h-1 FX (t1/2 ~41h)
KDEG_FX   = 0.017 // h-1 FX degradation
FX_init   = 100.0 // % normal FX

// ── Renal Function Adjustment (eGFR) ────────────────────────────────────
eGFR_pat = 90.0   // mL/min/1.73m2 (patient eGFR)
eGFR_ref = 90.0   // reference eGFR

// ── Patient Parameters ──────────────────────────────────────────────────
BWT      = 70.0   // kg body weight
AGE_pat  = 55.0   // years
CRCL_adj = 1.0    // creatinine clearance adjustment (1=normal)

$CMT
// PK compartments
RIV_GUT RIV_CENT RIV_PERIPH
APIX_CENT
DABI_CENT
WARF_CENT
ENOX_CENT
// PD compartments (coagulation/fibrinolysis)
FXa_ACT
FIIa_ACT
FIBRIN_FORM
CLOT_SIZE
PLASMIN_ACT
DDIMER_CONC
// Vitamin K cycle (Warfarin PD)
VK_OX
VK_RED
// Factor pools (Warfarin indirect effect)
FVII_POOL
FX_POOL
FII_POOL

$INIT
// PK: zero dose at start
RIV_GUT   = 0
RIV_CENT  = 0
RIV_PERIPH= 0
APIX_CENT = 0
DABI_CENT = 0
WARF_CENT = 0
ENOX_CENT = 0
// PD: at steady-state (no drug, with clot present)
FXa_ACT   = 1.0
FIIa_ACT  = 1.0
FIBRIN_FORM = 0.0
CLOT_SIZE = CLOT_init   // Start with established clot
PLASMIN_ACT = PLASMIN_base
DDIMER_CONC = 2.5       // ng/mL elevated at VTE diagnosis
// Vitamin K cycle (baseline)
VK_OX     = VK0_ox
VK_RED    = VK0_red
// Factor pools at 100%
FVII_POOL = FVII_init
FX_POOL   = FX_init
FII_POOL  = FII_init

$MAIN
// ── Renal function adjustment ────────────────────────────────────────────
double RF_adj = eGFR_pat / eGFR_ref;   // renal function ratio
double RF_adj_DABI = pow(RF_adj, 0.85); // stronger renal effect on dabi
double RF_adj_ENOX = pow(RF_adj, 0.65); // renal effect on LMWH

// ── Rivaroxaban Cp (ng/mL) ───────────────────────────────────────────────
double Cp_RIV = RIV_CENT / (V1_RIV * (BWT / WT_ref));

// ── Apixaban Cp (ng/mL) ─────────────────────────────────────────────────
double Cp_APIX = APIX_CENT / V1_APIX;

// ── Dabigatran Cp (ng/mL) ───────────────────────────────────────────────
double Cp_DABI = DABI_CENT / (V1_DABI * (70.0 / BWT));

// ── Warfarin Cp (mg/L) ──────────────────────────────────────────────────
double Cp_WARF = WARF_CENT / V1_WARF;

// ── Enoxaparin Anti-Xa (IU/mL) ──────────────────────────────────────────
double Cp_ENOX = ENOX_CENT / V1_ENOX;

// ── FXa inhibition fraction ─────────────────────────────────────────────
double INH_FXa_RIV  = (EMAX_RIV * pow(Cp_RIV, HILL_RIV)) /
                      (pow(EC50_RIV, HILL_RIV) + pow(Cp_RIV, HILL_RIV));
double INH_FXa_APIX = (EMAX_APIX * pow(Cp_APIX, HILL_APX)) /
                      (pow(EC50_APIX, HILL_APX) + pow(Cp_APIX, HILL_APX));
double INH_FXa_ENOX = (EMAX_ENOX * Cp_ENOX) / (EC50_ENOX + Cp_ENOX);
double INH_FXa_TOT  = 1.0 - (1.0 - INH_FXa_RIV) *
                             (1.0 - INH_FXa_APIX) *
                             (1.0 - INH_FXa_ENOX);  // Bliss independence

// ── FIIa inhibition fraction (Dabigatran + indirect from FXa inhib) ─────
double INH_FIIa_DABI= (EMAX_DABI * pow(Cp_DABI, HILL_DABI)) /
                      (pow(EC50_DABI, HILL_DABI) + pow(Cp_DABI, HILL_DABI));
double INH_FIIa_TOT = 1.0 - (1.0 - INH_FIIa_DABI) *
                             (1.0 - INH_FXa_TOT * 0.7); // FXa inh → ↓ FIIa gen

// ── Warfarin indirect effect on factor synthesis ─────────────────────────
double WARF_INH = (EMAX_RIV * pow(Cp_WARF, HILL_W)) /
                  (pow(IC50_WARF, HILL_W) + pow(Cp_WARF, HILL_W));
// Factor synthesis inhibited proportional to VK_RED depletion
double VK_ratio = VK_RED / VK0_red;
double INH_FVII = 1.0 - VK_ratio;    // % inhibition of FVII synthesis
double INH_FX   = 1.0 - VK_ratio;
double INH_FII  = 1.0 - VK_ratio;

// ── INR calculation (empirical, based on factor pools) ──────────────────
// INR ≈ 1 when all factors 100%, increases as FVII/FX/FII drop
double PT_norm  = 12.0; // s
double PT_pct   = (FVII_POOL / 100.0) * 0.5 +
                  (FX_POOL   / 100.0) * 0.3 +
                  (FII_POOL  / 100.0) * 0.2;
double PT_pat   = PT_norm / (PT_pct + 0.001);
double INR      = PT_pat / PT_norm;

// ── Anti-Xa level (for enoxaparin monitoring) ───────────────────────────
double ANTI_XA = Cp_ENOX;

// ── aPTT (for heparin monitoring) ───────────────────────────────────────
// (simplified: aPTT increases with thrombin inhibition)
double aPTT     = 32.0 * (1.0 + 2.5 * INH_FIIa_DABI);

// ── Thrombus resolution score (0-100%) ──────────────────────────────────
double CLOT_PCT = 100.0 * CLOT_SIZE / CLOT_init;

$ODE
// ════════════════════════════════════════════════════════════════════════
// PHARMACOKINETICS
// ════════════════════════════════════════════════════════════════════════

// ── Rivaroxaban 2-compartment ────────────────────────────────────────────
dxdt_RIV_GUT   = -KA_RIV * RIV_GUT;
dxdt_RIV_CENT  =  KA_RIV * RIV_GUT
                 - (CL_RIV / V1_RIV) * RIV_CENT
                 - (Q_RIV  / V1_RIV) * RIV_CENT
                 + (Q_RIV  / V2_RIV) * RIV_PERIPH;
dxdt_RIV_PERIPH=  (Q_RIV  / V1_RIV) * RIV_CENT
                 - (Q_RIV  / V2_RIV) * RIV_PERIPH;

// ── Apixaban 1-compartment ───────────────────────────────────────────────
dxdt_APIX_CENT = -KA_APIX * APIX_CENT
                 - (CL_APIX / V1_APIX) * APIX_CENT;

// Note: dosing enters APIX_CENT directly as bolus with F_APIX applied
// (handled via F in $CMT or explicit dose × F at time of dosing)

// ── Dabigatran 1-compartment (prodrug → active) ──────────────────────────
dxdt_DABI_CENT = -KA_DABI * DABI_CENT
                 - (CL_DABI / V1_DABI) / RF_adj_DABI * DABI_CENT;

// ── Warfarin 1-compartment ───────────────────────────────────────────────
dxdt_WARF_CENT = -KA_WARF * WARF_CENT
                 - (CL_WARF / V1_WARF) * WARF_CENT;

// ── Enoxaparin SC 1-compartment ──────────────────────────────────────────
dxdt_ENOX_CENT = -KA_ENOX * ENOX_CENT
                 - (CL_ENOX / V1_ENOX) / RF_adj_ENOX * ENOX_CENT;

// ════════════════════════════════════════════════════════════════════════
// COAGULATION PHARMACODYNAMICS
// ════════════════════════════════════════════════════════════════════════

// ── Factor Xa Activity ───────────────────────────────────────────────────
// FXa generation: baseline + stimulus (VTE) − inhibition
double FXa_gen = KG_FXa * FXa_base * (FX_POOL / 100.0);
double FXa_clr = KF_FXa * FXa_ACT;
dxdt_FXa_ACT   = FXa_gen * (1.0 - INH_FXa_TOT) - FXa_clr;

// ── Thrombin (FIIa) Activity ─────────────────────────────────────────────
// FIIa generation proportional to FXa and FII pool
double FIIa_gen = KG_FIIa * FXa_ACT * (FII_POOL / 100.0);
double FIIa_clr = KF_FIIa * FIIa_ACT;
dxdt_FIIa_ACT  = FIIa_gen * (1.0 - INH_FIIa_TOT) - FIIa_clr;

// ── Fibrin Formation ─────────────────────────────────────────────────────
// Fibrin formation = proportional to FIIa; lysis = proportional to plasmin
double FIBRIN_form_rate = KF_FBR * FIIa_ACT * FIIa_base;
double FIBRIN_lysis_rate= KL_FBR * PLASMIN_ACT * FIBRIN_FORM;
dxdt_FIBRIN_FORM= FIBRIN_form_rate - FIBRIN_lysis_rate;

// ── Thrombus/Clot Burden ─────────────────────────────────────────────────
// Clot grows with ongoing fibrin formation; decays with plasmin lysis
double CLOT_growth = KG_CLOT * FIBRIN_FORM;
double CLOT_decay  = KD_CLOT * PLASMIN_ACT * CLOT_SIZE;
dxdt_CLOT_SIZE = CLOT_growth - CLOT_decay;

// ── Fibrinolysis: Plasmin Activity ───────────────────────────────────────
// Plasmin forms from tPA-mediated activation of plasminogen
// Thrombolytic therapy increases Plasmin_base by ~20×
double PLASMIN_form = KP_FORM * (1.0 - PAI1_eff * 0.3);  // PAI-1 partial inhib
double PLASMIN_inhib= KP_DECAY * PLASMIN_ACT;
dxdt_PLASMIN_ACT  = PLASMIN_form - PLASMIN_inhib;

// ── D-dimer Dynamics ─────────────────────────────────────────────────────
// D-dimer generated by plasmin degradation of cross-linked fibrin
double DDIMER_gen = K_DDIMER * PLASMIN_ACT * FIBRIN_FORM;
double DDIMER_clr = K_DDIMER_CL * DDIMER_CONC;
dxdt_DDIMER_CONC  = DDIMER_gen - DDIMER_clr;

// ════════════════════════════════════════════════════════════════════════
// VITAMIN K CYCLE (WARFARIN INDIRECT EFFECT)
// ════════════════════════════════════════════════════════════════════════
// VKORC1 inhibition: ↓ VitK_ox → VitK_red conversion
double WARF_effect = WARF_INH;  // 0 to ~0.97 inhibition
dxdt_VK_OX  =  KOUT_VK * VK_RED - KIN_VK  * VK_OX;
dxdt_VK_RED =  KIN_VK  * VK_OX  * (1.0 - WARF_effect) - KOUT_VK * VK_RED;

// ════════════════════════════════════════════════════════════════════════
// FACTOR SYNTHESIS (INDIRECT RESPONSE via VK_RED depletion)
// ════════════════════════════════════════════════════════════════════════
// Factor synthesis rate reduced when VK_RED depleted (warfarin effect)
double STIM_VK = VK_RED / VK0_red;   // 0→1: 0 = full warfarin effect
dxdt_FVII_POOL =  KSYN_FVII * STIM_VK * FVII_init
                 - KDEG_FVII * FVII_POOL;
dxdt_FX_POOL   =  KSYN_FX   * STIM_VK * FX_init
                 - KDEG_FX   * FX_POOL;
dxdt_FII_POOL  =  KSYN_FII  * STIM_VK * FII_init
                 - KDEG_FII  * FII_POOL;

$TABLE
// ── PK outputs ──────────────────────────────────────────────────────────
capture Cp_RIV   = Cp_RIV;     // ng/mL Rivaroxaban
capture Cp_APIX  = Cp_APIX;    // ng/mL Apixaban
capture Cp_DABI  = Cp_DABI;    // ng/mL Dabigatran
capture Cp_WARF  = Cp_WARF;    // mg/L Warfarin
capture ANTI_XA  = ANTI_XA;    // IU/mL Enoxaparin anti-Xa

// ── PD biomarkers ───────────────────────────────────────────────────────
capture INR       = INR;            // PT/INR ratio
capture aPTT_out  = aPTT;          // aPTT (s)
capture INH_FXa   = INH_FXa_TOT*100; // % FXa inhibition
capture INH_FIIa  = INH_FIIa_TOT*100;// % FIIa inhibition
capture TG_ETP    = FIIa_ACT * 100;  // % Thrombin generation ETP proxy
capture DDIMER    = DDIMER_CONC;    // ng/mL D-dimer
capture CLOT_PCT2 = CLOT_SIZE / CLOT_init * 100; // % residual clot
capture FIBRIN_mg = FIBRIN_FORM;    // mg/dL fibrin pool
capture FVII_pct  = FVII_POOL;     // % normal FVII
capture FX_pct    = FX_POOL;       // % normal FX
capture FII_pct   = FII_POOL;      // % normal FII

$CAPTURE
Cp_RIV Cp_APIX Cp_DABI Cp_WARF ANTI_XA
INR aPTT_out INH_FXa INH_FIIa TG_ETP
DDIMER CLOT_PCT2 FIBRIN_mg
FVII_pct FX_pct FII_pct
'

## ══════════════════════════════════════════════════════════════════════════
## Compile model
## ══════════════════════════════════════════════════════════════════════════

mod <- mcode("vte_qsp", vte_model_code)

## ══════════════════════════════════════════════════════════════════════════
## Helper: build dosing regimen
## ══════════════════════════════════════════════════════════════════════════

make_dosing <- function(drug, scenario, duration_days = 90) {
  dur_h <- duration_days * 24

  if (drug == "rivaroxaban" && scenario == "DVT_treatment") {
    # Phase 1: 15 mg BID × 21 days (with food), then 20 mg QD
    d1 <- ev(amt = 15 * 0.93, cmt = "RIV_GUT", ii = 12, addl = 41, time = 0)
    d2 <- ev(amt = 20 * 0.93, cmt = "RIV_GUT", ii = 24, addl = (duration_days - 21) - 1,
             time = 21 * 24)
    c(d1, d2)

  } else if (drug == "rivaroxaban" && scenario == "prophylaxis") {
    # 10 mg QD
    ev(amt = 10 * 0.93, cmt = "RIV_GUT", ii = 24,
       addl = duration_days - 1, time = 0)

  } else if (drug == "apixaban" && scenario == "PE_treatment") {
    # 10 mg BID × 7 days, then 5 mg BID × 83 days
    d1 <- ev(amt = 10 * 0.50, cmt = "APIX_CENT", ii = 12, addl = 13, time = 0)
    d2 <- ev(amt = 5  * 0.50, cmt = "APIX_CENT", ii = 12, addl = (duration_days - 7)*2 - 1,
             time = 7 * 24)
    c(d1, d2)

  } else if (drug == "dabigatran" && scenario == "DVT_treatment") {
    # After 5-10d parenteral: 150 mg BID
    ev(amt = 150 * 0.065, cmt = "DABI_CENT", ii = 12,
       addl = (duration_days - 0) * 2 - 1, time = 0)

  } else if (drug == "warfarin" && scenario == "long_term") {
    # 5 mg QD (loading), adjust to INR 2-3
    ev(amt = 5 * 1.0, cmt = "WARF_CENT", ii = 24,
       addl = duration_days - 1, time = 0)

  } else if (drug == "enoxaparin" && scenario == "bridge") {
    # 1 mg/kg BID SC (70 kg patient)
    ev(amt = 70 * 1.0 * 0.92, cmt = "ENOX_CENT", ii = 12,
       addl = 9 * 2 - 1, time = 0)   # 10 days bridging

  } else if (drug == "enoxaparin" && scenario == "prophylaxis") {
    # 40 mg QD SC (surgical prophylaxis)
    ev(amt = 40 * 0.92, cmt = "ENOX_CENT", ii = 24,
       addl = 13, time = 0)  # 14 days post-op
  }
}

## ══════════════════════════════════════════════════════════════════════════
## SCENARIO 1: DVT Treatment — Rivaroxaban (Phase I: 15mg BID / Phase II: 20mg QD)
## ══════════════════════════════════════════════════════════════════════════

cat("\n═══ Scenario 1: DVT Treatment — Rivaroxaban ═══\n")
dose_riva <- make_dosing("rivaroxaban", "DVT_treatment", duration_days = 90)

sim1 <- mod %>%
  param(CLOT_init = 80, CLOT_SIZE = 80) %>%
  mrgsim(events = dose_riva, end = 90 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "DVT: Rivaroxaban (15 BID→20 QD)", time_d = time / 24)

cat("Peak Rivaroxaban Cp (Day 1, BID):", round(max(sim1$Cp_RIV[sim1$time_d < 1]), 1), "ng/mL\n")
cat("Day 7 Clot Reduction:", round(100 - sim1$CLOT_PCT2[sim1$time == 168], 1), "%\n")
cat("Day 30 Clot Reduction:", round(100 - sim1$CLOT_PCT2[sim1$time == 720], 1), "%\n")

## ══════════════════════════════════════════════════════════════════════════
## SCENARIO 2: PE Treatment — Apixaban (10 mg BID × 7d, then 5 mg BID)
## ══════════════════════════════════════════════════════════════════════════

cat("\n═══ Scenario 2: PE Treatment — Apixaban ═══\n")
dose_apix <- make_dosing("apixaban", "PE_treatment", duration_days = 90)

sim2 <- mod %>%
  param(CLOT_init = 80, CLOT_SIZE = 80, PLASMIN_base = 0.4) %>%
  mrgsim(events = dose_apix, end = 90 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "PE: Apixaban (10 BID×7d→5 BID)", time_d = time / 24)

cat("Peak Apixaban Cp (Day 1):", round(max(sim2$Cp_APIX[sim2$time_d < 1]), 1), "ng/mL\n")

## ══════════════════════════════════════════════════════════════════════════
## SCENARIO 3: Warfarin + LMWH Bridge (Target INR 2-3)
## ══════════════════════════════════════════════════════════════════════════

cat("\n═══ Scenario 3: Warfarin + Enoxaparin Bridge ═══\n")
dose_warf <- make_dosing("warfarin", "long_term", duration_days = 90)
dose_enox <- make_dosing("enoxaparin", "bridge")

dose_combo <- c(dose_warf, dose_enox)

sim3 <- mod %>%
  param(CLOT_init = 80, CLOT_SIZE = 80) %>%
  mrgsim(events = dose_combo, end = 90 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Warfarin + LMWH Bridge", time_d = time / 24)

cat("Day 7 INR:", round(sim3$INR[sim3$time == 168], 2), "\n")
cat("Day 30 INR:", round(sim3$INR[sim3$time == 720], 2), "\n")
cat("Day 7 FVII%:", round(sim3$FVII_pct[sim3$time == 168], 1), "%\n")

## ══════════════════════════════════════════════════════════════════════════
## SCENARIO 4: Surgical VTE Prophylaxis — Enoxaparin 40 mg QD
## ══════════════════════════════════════════════════════════════════════════

cat("\n═══ Scenario 4: Surgical Prophylaxis — Enoxaparin 40 mg QD ═══\n")
dose_enox_proph <- make_dosing("enoxaparin", "prophylaxis")

sim4 <- mod %>%
  param(CLOT_init = 0, CLOT_SIZE = 0, FIBRIN_FORM = 0) %>%
  mrgsim(events = dose_enox_proph, end = 14 * 24, delta = 0.5) %>%
  as.data.frame() %>%
  mutate(scenario = "Prophylaxis: Enoxaparin 40 mg QD", time_d = time / 24)

cat("Day 1 Peak Anti-Xa:", round(max(sim4$ANTI_XA[sim4$time_d <= 1]), 2), "IU/mL\n")
cat("Day 7 Peak Anti-Xa:", round(max(sim4$ANTI_XA[sim4$time_d >= 6 & sim4$time_d <= 7]), 2), "IU/mL\n")

## ══════════════════════════════════════════════════════════════════════════
## SCENARIO 5: Extended VTE Prevention — Rivaroxaban 10 mg QD
## ══════════════════════════════════════════════════════════════════════════

cat("\n═══ Scenario 5: Extended Prevention — Rivaroxaban 10 mg QD ═══\n")
dose_riva_ext <- make_dosing("rivaroxaban", "prophylaxis", duration_days = 180)

sim5 <- mod %>%
  param(CLOT_init = 0, CLOT_SIZE = 0) %>%
  mrgsim(events = dose_riva_ext, end = 180 * 24, delta = 2) %>%
  as.data.frame() %>%
  mutate(scenario = "Extended: Rivaroxaban 10 mg QD", time_d = time / 24)

cat("Steady-state Cp_RIV trough:", round(mean(tail(sim5$Cp_RIV[sim5$time_d %% 1 < 0.1], 30)), 1), "ng/mL\n")
cat("Steady-state FXa inhibition:", round(mean(tail(sim5$INH_FXa, 500)), 1), "%\n")

## ══════════════════════════════════════════════════════════════════════════
## SCENARIO 6: Renally Impaired Patient — Dabigatran 110 mg BID (CKD Stage 3)
## ══════════════════════════════════════════════════════════════════════════

cat("\n═══ Scenario 6: CKD3 — Dabigatran 110 mg BID ═══\n")
dose_dabi_ckd <- ev(amt = 110 * 0.065, cmt = "DABI_CENT", ii = 12,
                    addl = 89 * 2, time = 0)
sim6_norm <- mod %>%
  param(eGFR_pat = 90, CLOT_init = 80, CLOT_SIZE = 80) %>%
  mrgsim(events = dose_dabi_ckd, end = 90 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Dabigatran 110mg BID (Normal GFR)", time_d = time / 24)

sim6_ckd <- mod %>%
  param(eGFR_pat = 30, CLOT_init = 80, CLOT_SIZE = 80) %>%
  mrgsim(events = dose_dabi_ckd, end = 90 * 24, delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "Dabigatran 110mg BID (CKD3 GFR=30)", time_d = time / 24)

cat("Normal GFR SS Cp_DABI:", round(mean(tail(sim6_norm$Cp_DABI, 200)), 1), "ng/mL\n")
cat("CKD3 GFR=30 SS Cp_DABI:", round(mean(tail(sim6_ckd$Cp_DABI, 200)), 1), "ng/mL\n")

## ══════════════════════════════════════════════════════════════════════════
## COMBINE & PLOT
## ══════════════════════════════════════════════════════════════════════════

all_sims <- bind_rows(sim1, sim2, sim3)

# ── Plot 1: PK profiles ──────────────────────────────────────────────────
p_pk <- all_sims %>%
  filter(time_d <= 10) %>%
  select(time_d, scenario, Cp_RIV, Cp_APIX, INR) %>%
  pivot_longer(c(Cp_RIV, Cp_APIX, INR)) %>%
  ggplot(aes(x = time_d, y = value, color = scenario, linetype = name)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~name, scales = "free_y",
             labeller = labeller(name = c(Cp_RIV = "Rivaroxaban (ng/mL)",
                                           Cp_APIX= "Apixaban (ng/mL)",
                                           INR    = "INR"))) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "VTE Drug PK: First 10 Days",
       x = "Time (days)", y = "Concentration / Ratio",
       color = "Treatment", linetype = "Analyte") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ── Plot 2: Clot resolution ──────────────────────────────────────────────
p_clot <- all_sims %>%
  ggplot(aes(x = time_d, y = CLOT_PCT2, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  annotate("text", x = 45, y = 52, label = "50% Clot Reduction", size = 3) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "VTE Thrombus Resolution Over 90 Days",
       x = "Time (days)", y = "Residual Clot (%)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ── Plot 3: D-dimer dynamics ─────────────────────────────────────────────
p_ddimer <- all_sims %>%
  ggplot(aes(x = time_d, y = DDIMER, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "darkgreen",
             linewidth = 0.8) +
  annotate("text", x = 60, y = 0.6, label = "Normal upper limit (500 ng/mL)",
           size = 3, color = "darkgreen") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "D-dimer Dynamics During VTE Treatment",
       x = "Time (days)", y = "D-dimer (μg/mL FEU)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ── Plot 4: Warfarin factor pool depletion ───────────────────────────────
p_warf <- sim3 %>%
  select(time_d, FVII_pct, FX_pct, FII_pct, INR) %>%
  pivot_longer(c(FVII_pct, FX_pct, FII_pct, INR)) %>%
  ggplot(aes(x = time_d, y = value, color = name)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = c(2, 3), linetype = "dashed",
             color = c("red", "red"), alpha = 0.7) +
  annotate("text", x = 50, y = 2.1, label = "INR 2.0", size = 3, color = "red") +
  annotate("text", x = 50, y = 3.1, label = "INR 3.0", size = 3, color = "red") +
  scale_color_manual(values = c(FVII_pct = "#e53935", FX_pct = "#1e88e5",
                                 FII_pct = "#43a047", INR = "#8e24aa"),
                     labels = c("FVII (%)", "FX (%)", "FII (%)", "INR")) +
  labs(title = "Warfarin — Vitamin K-Dependent Factor Depletion & INR",
       x = "Time (days)", y = "Factor Level (%) or INR",
       color = "Parameter") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# ── Plot 5: Renally impaired comparison ──────────────────────────────────
p_renal <- bind_rows(sim6_norm, sim6_ckd) %>%
  filter(time_d <= 7) %>%
  ggplot(aes(x = time_d, y = Cp_DABI, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = c(30, 200), linetype = c("dotted", "dashed"),
             color = c("blue", "red")) +
  annotate("text", x = 4, y = 35,  label = "Trough target >30 ng/mL", size = 3, color = "blue") +
  annotate("text", x = 4, y = 210, label = "Peak concern >200 ng/mL", size = 3, color = "red") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Dabigatran: Normal vs. CKD Stage 3 (GFR 30)",
       x = "Time (days)", y = "Dabigatran Cp (ng/mL)",
       color = "Renal Function") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

# Print plots
print(p_pk)
print(p_clot)
print(p_ddimer)
print(p_warf)
print(p_renal)

## ══════════════════════════════════════════════════════════════════════════
## Summary Table
## ══════════════════════════════════════════════════════════════════════════

summary_tbl <- bind_rows(sim1, sim2, sim3) %>%
  group_by(scenario) %>%
  summarise(
    Day7_Clot_Pct   = round(CLOT_PCT2[which.min(abs(time_d - 7))], 1),
    Day30_Clot_Pct  = round(CLOT_PCT2[which.min(abs(time_d - 30))], 1),
    Day90_Clot_Pct  = round(CLOT_PCT2[which.min(abs(time_d - 90))], 1),
    Peak_Ddimer     = round(max(DDIMER), 2),
    Day30_Ddimer    = round(DDIMER[which.min(abs(time_d - 30))], 2),
    .groups = "drop"
  )

cat("\n═══ Clot Resolution Summary ═══\n")
print(summary_tbl)

cat("\nModel run complete.\n")
cat("5 treatment scenarios simulated:\n")
cat("  1. DVT: Rivaroxaban 15mg BID×21d → 20mg QD\n")
cat("  2. PE : Apixaban 10mg BID×7d → 5mg BID\n")
cat("  3. Bridge: Warfarin 5mg QD + Enoxaparin 1mg/kg BID\n")
cat("  4. Prophylaxis: Enoxaparin 40mg QD\n")
cat("  5. Extended: Rivaroxaban 10mg QD\n")
cat("  6. CKD: Dabigatran 110mg BID (renal impairment sensitivity)\n")
