# =====================================================================
# Von Willebrand Disease (VWD) — mrgsolve QSP Model
#   Author : Claude Code Routine (2026-07-01)
#   Scope  : VWF gene (12p13.31) quantitative (Type 1/3) or qualitative
#            (Type 2A/2B/2M/2N) deficiency -> impaired platelet-collagen/
#            GPIba adhesion (primary hemostasis) and/or FVIII stabilization
#            (secondary hemostasis) -> mucocutaneous bleeding, menorrhagia,
#            GI angiodysplasia bleeding (Type 2A/3/acquired), postpartum
#            hemorrhage. ADAMTS13 cleaves the shear-unfolded A2 domain,
#            regulating high-molecular-weight multimer (HMWM) content.
#   PK/PD  : Desmopressin (DDAVP, IV/SC/intranasal; V2-like endothelial
#            receptor -> acute Weibel-Palade body release of VWF+FVIII+
#            VWFpp; tachyphylaxis; contraindicated in Type 2B/platelet-type)
#            . Plasma-derived VWF/FVIII concentrate (Humate-P/Wilate-like).
#            Recombinant VWF (vonicog alfa/Vonvendi-like; VWF-only, delayed
#            endogenous FVIII stabilization rise). Tranexamic acid
#            (antifibrinolytic, menorrhagia/surgical adjunct).
#   Outputs: VWF:Ag, VWF:RCo, HMWM fraction, FVIII:C, platelet count,
#            composite bleeding score (ISTH-BAT-like), menstrual blood loss,
#            GI blood loss, hemoglobin, serum sodium (DDAVP safety),
#            thrombotic-risk index (factor overcorrection safety signal).
#   References (calibration): Mannucci PM. Blood 1997 (PMID 9326215; DDAVP
#            pharmacology/kinetics, first 20 years); Federici AB, et al.
#            Blood 2004 (PMID 14630825; DDAVP biologic-response multicenter
#            study by VWD type); Gill JC, et al. Blood 2015 (PMID 26239086;
#            rVWF/vonicog alfa pivotal phase 3 PK/hemostatic efficacy);
#            Mannucci PM, et al. Blood 2013 (PMID 23777763; rVWF phase 1 PK,
#            VWF:RCo t1/2 ~21h); Dobrkovska A, et al. Haemophilia 1998
#            (PMID 10028316; Humate-P PK, VWF:RCo t1/2 ~12-20h, FVIII:C
#            t1/2 ~8-12h); James PD, Connell NT, et al. Blood Adv 2021
#            (PMID 33570651; ASH/ISTH/NHF/WFH 2021 diagnosis guideline);
#            Connell NT, Flood VH, et al. Blood Adv 2021 (PMID 33570647;
#            ASH/ISTH/NHF/WFH 2021 management guideline, DDAVP trial-of-
#            response, PK-guided perioperative dosing); Leebeek FW,
#            Eikenboom JC. N Engl J Med 2016 (PMID 27959741; comprehensive
#            VWD review); Lukes AS, et al. Obstet Gynecol 2010 (PMID
#            20859150; tranexamic acid for heavy menstrual bleeding RCT).
# =====================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

vwd_code <- '
$PROB
# Von Willebrand Disease (VWD) QSP model
# 21 ODE compartments: 6 drug PK + 15 disease/PD/clinical

$PLUGIN autodec

$PARAM @annotated
// ============================================
// Desmopressin (DDAVP) PK (1-cpt, IV/SC/IN) — Mannucci 1997
// ============================================
KA_DDAVP   : 4.5    : DDAVP absorption rate, IV~instant/IN·SC slower (1/h)
KE_DDAVP   : 0.23   : DDAVP elimination rate (1/h)          // t1/2 ~3 h
V_DDAVP    : 30.0   : DDAVP apparent Vd (L)
F_DDAVP    : 0.15   : Bioavailability (IN/SC route; IV=1 handled via dosing)

// ============================================
// Recombinant VWF (vonicog alfa-like) PK — Mannucci 2013 / Gill 2015
// ============================================
KE_RVWF    : 0.033  : rVWF:RCo elimination rate (1/h)        // t1/2 ~21 h
V_RVWF     : 35.0   : rVWF apparent Vd (dL; ~plasma volume, ~2 IU/dL rise per IU/kg dosed)

// ============================================
// Plasma-derived VWF/FVIII concentrate PK — Dobrkovska 1998 (Humate-P-like)
// ============================================
KE_PDVWF     : 0.045  : PD-concentrate VWF:RCo elimination rate (1/h)  // t1/2 ~15 h
V_PDVWF      : 35.0   : PD-concentrate VWF apparent Vd (dL; ~plasma volume)
VWF_FVIII_RATIO : 2.4 : VWF:RCo to FVIII:C ratio in PD concentrate (IU:IU)

// ============================================
// Tranexamic acid PK (oral/IV, antifibrinolytic)
// ============================================
KA_TXA     : 1.1    : TXA oral absorption rate (1/h)
KE_TXA     : 0.29   : TXA elimination rate (1/h)             // t1/2 ~2.4 h
V_TXA      : 25.0   : TXA apparent Vd (L)
F_TXA      : 0.35   : TXA oral bioavailability

// ============================================
// Disease baseline (genotype-dependent; overridden per scenario)
// ============================================
VWFAG_BASE   : 35.0  : Baseline VWF:Ag (IU/dL); Type1 mild~35, Type3~<3
VWFRCO_BASE  : 30.0  : Baseline VWF:RCo (IU/dL)
FVIII_BASE   : 40.0  : Baseline FVIII:C (IU/dL)
HMWM_BASE    : 0.55  : Baseline HMWM fraction (0-1); normal~1.0, Type2A/3 low
PLT_BASE     : 250   : Baseline platelet count (x10^9/L)
NA_BASE      : 140   : Baseline serum sodium (mmol/L)
HB_BASE      : 13.0  : Baseline hemoglobin (g/dL)
ADAMTS13_BASE: 1.0   : Baseline ADAMTS13 activity (relative, normal=1.0)
TYPE2B_FLAG  : 0     : Flag (0/1): Type 2B / platelet-type GPIba gain-of-function
ACQUIRED_FLAG: 0     : Flag (0/1): acquired VWS (shear/AS or lymphoproliferative)
TYPE2N_FLAG  : 0     : Flag (0/1): Type 2N (defective FVIII binding)

// ============================================
// WPB releasable-store & DDAVP pharmacodynamics — Federici 2008
// ============================================
WPB0          : 1.0   : Baseline normalized WPB releasable store (0-1)
EC50_DDAVP    : 0.8   : DDAVP conc. for half-max WPB release (ng/mL)
HILL_DDAVP    : 1.3   : Hill coefficient
EMAX_WPBREL   : 3.5   : Max fold-rise in VWF:Ag/RCo from full WPB release
K_WPBDEPLETE  : 0.15  : Rate constant, WPB store depletion per release event (1/h)
K_WPBREGEN    : 0.010 : Rate constant, WPB store regeneration (1/h, slow, days)
K_VWFPP_REL   : 2.0   : Rate constant, acute VWFpp release with WPB (1/h, fast)
K_VWFPP_CLEAR : 0.35  : VWFpp plasma clearance rate (1/h) // faster than mature VWF

// ============================================
// VWF:Ag / VWF:RCo / HMWM turnover
// ============================================
KOUT_VWFAG    : 0.058 : VWF:Ag turnover/clearance rate (1/h) // t1/2 ~12h endogenous
KOUT_VWFRCO   : 0.058 : VWF:RCo turnover/clearance rate (1/h)
KOUT_HMWM     : 0.10  : HMWM fraction relaxation rate (1/h)
K_ADAMTS13_CLEAVE : 0.35 : ADAMTS13-mediated HMWM cleavage rate constant (1/h, scaled by activity)
K_ADAMTS13_CONSUME: 0.20 : ADAMTS13 transient consumption rate post-WPB surge (1/h)
K_ADAMTS13_REGEN  : 0.05 : ADAMTS13 activity regeneration rate (1/h)

// ============================================
// FVIII stabilization
// ============================================
KOUT_FVIII    : 0.075 : Endogenous FVIII:C turnover rate (1/h) // t1/2 ~9h free, longer if VWF-bound
FVIII_STAB_GAIN: 1.15 : Max FVIII:C gain per unit VWF:Ag (stabilization coupling)
TYPE2N_FVIII_PENALTY: 0.35 : Fractional FVIII:C ceiling reduction in Type 2N

// ============================================
// Platelet count dynamics (Type 2B / platelet-type clearance)
// ============================================
KOUT_PLT      : 0.010 : Baseline platelet turnover rate (1/h)
K_PLT_CLEAR_2B: 0.30  : Type-2B/platelet-type aggregation-clearance rate constant (1/h)
EC50_2B_CLEAR : 1.0   : VWF surge (fold-baseline) for half-max Type-2B clearance
PLT_FLOOR_2B  : 40    : Minimum platelet count under severe 2B clearance (x10^9/L)

// ============================================
// Bleeding score / clinical endpoints
// ============================================
BLEED_BASE     : 4.0   : Baseline composite bleeding score (ISTH-BAT-like, a.u.)
BLEED_GAIN_VWF : 8.0    : Bleeding-score increment per unit VWF:RCo deficit (normalized)
BLEED_GAIN_FVIII: 3.0   : Bleeding-score increment per unit FVIII:C deficit
BLEED_GAIN_PLT : 4.0    : Bleeding-score increment per unit platelet deficit
KOUT_BLEED     : 0.25   : Bleeding-score equilibration rate (1/h)
EMAX_TXA_BLEED : 0.30   : Max fractional bleeding-score reduction from TXA
EC50_TXA       : 8.0    : TXA conc. for half-max antifibrinolytic effect (mcg/mL)

// ============================================
// Menorrhagia / GI bleed / hemoglobin
// ============================================
MENS_BASE      : 40    : Baseline menstrual blood loss (mL/cycle, normal ~30-40)
MENS_GAIN      : 220   : Max additional menstrual loss from severe deficiency (mL/cycle)
KOUT_MENS      : 0.02  : Menstrual-loss compartment relaxation rate (1/h, cycle-scale)
HORMONAL_ON    : 0     : Flag (0/1): hormonal therapy (COC/LNG-IUS) active
HORMONAL_REDUX : 0.55  : Fractional menstrual-loss reduction with hormonal therapy
GI_BASE        : 2.0   : Baseline chronic GI blood loss (mL/day)
GI_GAIN_HMWM   : 25.0  : Max additional GI loss from HMWM deficit (angiodysplasia)
KOUT_GI        : 0.04  : GI blood-loss compartment relaxation rate (1/h)
KOUT_HB        : 0.006 : Hemoglobin relaxation rate to chronic-loss-adjusted target (1/h, slow)
HB_LOSS_SENS   : 0.006 : Hemoglobin sensitivity to cumulative blood loss (g/dL per mL/day equiv)

// ============================================
// Safety: hyponatremia (DDAVP), thrombotic risk (over-correction)
// ============================================
EMAX_NA_DROP   : 8.0   : Max serum sodium drop with repeated DDAVP + free water (mmol/L)
EC50_NA        : 1.5   : DDAVP cumulative-exposure index for half-max Na effect
KOUT_NA        : 0.03  : Sodium re-equilibration rate (1/h)
FLUID_RESTRICT : 0     : Flag (0/1): post-dose fluid restriction (mitigates Na drop)
THROMB_GAIN    : 6.0   : Max thrombotic-risk index from factor-concentrate overcorrection
EC50_THROMB    : 200   : VWF:RCo (IU/dL) for half-max thrombotic-risk signal
KOUT_THROMB    : 0.08  : Thrombotic-risk index equilibration rate (1/h)

// ============================================
// Pregnancy / physiologic modulator
// ============================================
PREG_BOOST     : 1.0   : Multiplier on VWF:Ag/RCo/FVIII synthesis (1.0=non-pregnant, ~2.5=3rd trimester)
POSTPARTUM_ON  : 0     : Flag (0/1): postpartum decline phase active

$CMT @annotated
DDAVP_DEPOT : DDAVP SC/IN depot (mcg)
DDAVP_CP    : DDAVP central concentration (ng/mL equiv)
RVWF_CP     : Recombinant VWF circulating amount (IU; divide by V_RVWF for conc.)
PDVWF_CP    : Plasma-derived VWF/FVIII concentrate circulating amount (IU; divide by V_PDVWF for conc.)
TXA_DEPOT   : Tranexamic acid oral depot (mg)
TXA_CP      : Tranexamic acid central concentration (mcg/mL)
WPB_STORE   : Normalized Weibel-Palade body releasable store (0-1)
VWFPP       : VWF propeptide (IU/dL equiv, acute release marker)
VWF_AG      : Plasma VWF antigen (IU/dL)
VWF_RCO     : Plasma VWF ristocetin-cofactor / GPIbM activity (IU/dL)
HMWM        : High-molecular-weight multimer fraction (0-1)
ADAMTS13_ACT: ADAMTS13 activity (relative, normal=1.0)
FVIII_C     : Plasma FVIII coagulant activity (IU/dL)
PLT_COUNT   : Platelet count (x10^9/L)
BLEED_SCORE : Composite bleeding-severity index (ISTH-BAT-like, a.u.)
MENS_LOSS   : Menstrual blood loss (mL/cycle, quasi-steady index)
GI_LOSS     : Chronic GI blood loss (mL/day, quasi-steady index)
HB          : Hemoglobin (g/dL)
NA_SERUM    : Serum sodium (mmol/L)
THROMB_RISK : Thrombotic-risk index (a.u., factor-overcorrection safety signal)

$MAIN
F_DDAVP_DEPOT = F_DDAVP;
F_TXA_DEPOT   = F_TXA;

if (NEWIND <= 1) {
  WPB_STORE   = WPB0;
  VWFPP       = VWFAG_BASE * 0.10;
  VWF_AG      = VWFAG_BASE;
  VWF_RCO     = VWFRCO_BASE;
  HMWM        = HMWM_BASE;
  ADAMTS13_ACT= ADAMTS13_BASE;
  FVIII_C     = FVIII_BASE;
  PLT_COUNT   = PLT_BASE;
  BLEED_SCORE = BLEED_BASE;
  MENS_LOSS   = MENS_BASE;
  GI_LOSS     = GI_BASE;
  HB          = HB_BASE;
  NA_SERUM    = NA_BASE;
  THROMB_RISK = 0;
}

$ODE
// ---- PK ----
dxdt_DDAVP_DEPOT = -KA_DDAVP * DDAVP_DEPOT;
dxdt_DDAVP_CP    =  KA_DDAVP * DDAVP_DEPOT / V_DDAVP - KE_DDAVP * DDAVP_CP;

dxdt_RVWF_CP     = -KE_RVWF * RVWF_CP;     // dosed directly (IU) via bolus into RVWF_CP
dxdt_PDVWF_CP    = -KE_PDVWF * PDVWF_CP;   // dosed directly (IU) via bolus into PDVWF_CP
double RVWF_CONC  = RVWF_CP / V_RVWF;      // IU/dL
double PDVWF_CONC = PDVWF_CP / V_PDVWF;    // IU/dL

dxdt_TXA_DEPOT   = -KA_TXA * TXA_DEPOT;
dxdt_TXA_CP      =  KA_TXA * TXA_DEPOT / V_TXA - KE_TXA * TXA_CP;

// ---- DDAVP -> WPB release (V2-like endothelial receptor) ----
double DDAVP_DRIVE = pow(DDAVP_CP, HILL_DDAVP) / (pow(EC50_DDAVP, HILL_DDAVP) + pow(DDAVP_CP, HILL_DDAVP));
double WPB_RELEASE_RATE = K_WPBDEPLETE * WPB_STORE * DDAVP_DRIVE;
dxdt_WPB_STORE = K_WPBREGEN * (WPB0 - WPB_STORE) - WPB_RELEASE_RATE;

// ---- VWFpp: acute co-release marker, fast clearance ----
dxdt_VWFPP = K_VWFPP_REL * WPB_RELEASE_RATE * EMAX_WPBREL * VWFAG_BASE - K_VWFPP_CLEAR * VWFPP;

// ---- ADAMTS13: transient consumption after large WPB surge (ULVWF substrate), slow regen ----
dxdt_ADAMTS13_ACT = K_ADAMTS13_REGEN * (ADAMTS13_BASE - ADAMTS13_ACT) - K_ADAMTS13_CONSUME * WPB_RELEASE_RATE * ADAMTS13_ACT;

// ---- Endogenous VWF:Ag: synthesis (pregnancy-boosted) + acute WPB release - clearance ----
double VWFAG_SYNTH_TARGET = VWFAG_BASE * PREG_BOOST;
dxdt_VWF_AG = KOUT_VWFAG * (VWFAG_SYNTH_TARGET - VWF_AG) + WPB_RELEASE_RATE * EMAX_WPBREL * VWFAG_BASE;

double EXO_VWF_TOTAL = RVWF_CONC + PDVWF_CONC;    // exogenous VWF:RCo-equivalent contribution (own PK decay)
double VWF_AG_TOTAL  = VWF_AG + EXO_VWF_TOTAL;    // reported total antigen = endogenous + replacement

// ---- Endogenous VWF:RCo activity: tracks Ag but discounted by HMWM deficit (qualitative defects reduce RCo/Ag ratio) ----
double VWFRCO_TARGET = VWF_AG * (0.4 + 0.6 * HMWM) * (VWFRCO_BASE / (VWFAG_BASE + 1e-6));
dxdt_VWF_RCO = KOUT_VWFRCO * (VWFRCO_TARGET - VWF_RCO);
double VWF_RCO_TOTAL = VWF_RCO + EXO_VWF_TOTAL;   // replacement product multimers intact -> full RCo credit

// ---- HMWM: shear/ADAMTS13 cleavage of endogenous pool; effective multimer competency boosted by exogenous multimer-replete product ----
double HMWM_TARGET_ENDO = HMWM_BASE * (1 - 0.5 * ACQUIRED_FLAG);
double HMWM_CLEAVE  = K_ADAMTS13_CLEAVE * ADAMTS13_ACT * HMWM * (WPB_RELEASE_RATE * 5.0 + 0.05);
double HMWM_RESTORE = KOUT_HMWM * (EXO_VWF_TOTAL / (VWF_AG + EXO_VWF_TOTAL + 1e-6));
dxdt_HMWM = KOUT_HMWM * (HMWM_TARGET_ENDO - HMWM) - HMWM_CLEAVE + HMWM_RESTORE;
double HMWM_EFFECTIVE = fmin(1.0, HMWM + EXO_VWF_TOTAL / (VWFAG_BASE * 3 + 1e-6));

// ---- FVIII:C: endogenous stabilization target raised by BOTH endogenous VWF:Ag and exogenous rVWF (RVWF_CONC);
//      because FVIII_C only relaxes toward its target at rate KOUT_FVIII (t1/2 ~9h), a pure rVWF dose (no
//      co-formulated FVIII) produces a genuinely delayed/gradual endogenous FVIII:C climb -- reproducing the
//      delayed/gradual rise reported for VWF-only rVWF product (Gill 2015) without an artificial bypass term.
//      Type-2N caps the achievable ceiling (defective D-prime-D3/FVIII binding). ----
double FVIII_CEILING = FVIII_BASE * (1 - TYPE2N_FVIII_PENALTY * TYPE2N_FLAG);
double FVIII_STAB_DRIVE = fmin(1.0, (VWF_AG + FVIII_STAB_GAIN * RVWF_CONC) / (VWFAG_BASE + 1e-6));
double FVIII_STAB_TARGET = FVIII_CEILING * (0.3 + 0.7 * FVIII_STAB_DRIVE) * PREG_BOOST;
dxdt_FVIII_C = KOUT_FVIII * (FVIII_STAB_TARGET - FVIII_C);

double FVIII_EXO_DIRECT = PDVWF_CONC / VWF_FVIII_RATIO;   // co-formulated FVIII in PD concentrate: immediate, own PK decay
double FVIII_C_TOTAL = FVIII_C + FVIII_EXO_DIRECT;

// ---- Platelet count: baseline homeostasis; Type-2B/platelet-type surge-triggered clearance ----
double VWF_FOLD = VWF_AG_TOTAL / (VWFAG_BASE + 1e-6);
double CLEAR_2B = TYPE2B_FLAG * K_PLT_CLEAR_2B * pow(VWF_FOLD, 2) / (pow(EC50_2B_CLEAR, 2) + pow(VWF_FOLD, 2));
dxdt_PLT_COUNT = KOUT_PLT * (PLT_BASE - PLT_COUNT) - CLEAR_2B * (PLT_COUNT - PLT_FLOOR_2B);

// ---- Composite bleeding score: rises with VWF:RCo/FVIII:C/platelet deficits, reduced by TXA ----
double VWF_DEFICIT    = fmax(0.0, 1 - VWF_RCO_TOTAL / 100.0);
double FVIII_DEFICIT  = fmax(0.0, 1 - FVIII_C_TOTAL / 100.0);
double PLT_DEFICIT    = fmax(0.0, 1 - PLT_COUNT / PLT_BASE);
double BLEED_TARGET_RAW = BLEED_BASE + BLEED_GAIN_VWF * VWF_DEFICIT + BLEED_GAIN_FVIII * FVIII_DEFICIT + BLEED_GAIN_PLT * PLT_DEFICIT;
double TXA_EFFECT = EMAX_TXA_BLEED * TXA_CP / (EC50_TXA + TXA_CP);
double BLEED_TARGET = BLEED_TARGET_RAW * (1 - TXA_EFFECT);
dxdt_BLEED_SCORE = KOUT_BLEED * (BLEED_TARGET - BLEED_SCORE);

// ---- Menstrual blood loss: driven by bleeding score, reduced by hormonal therapy ----
double MENS_TARGET_RAW = MENS_BASE + MENS_GAIN * VWF_DEFICIT * (0.5 + 0.5*FVIII_DEFICIT);
double MENS_TARGET = MENS_TARGET_RAW * (1 - HORMONAL_REDUX * HORMONAL_ON) * (1 - TXA_EFFECT);
dxdt_MENS_LOSS = KOUT_MENS * (MENS_TARGET - MENS_LOSS);

// ---- GI blood loss: driven by effective HMWM deficit (angiodysplasia, Type 2A/3/acquired) ----
double GI_TARGET = GI_BASE + GI_GAIN_HMWM * fmax(0.0, HMWM_BASE - HMWM_EFFECTIVE) / (HMWM_BASE + 1e-6) * (0.5 + 0.5*ACQUIRED_FLAG);
dxdt_GI_LOSS = KOUT_GI * (GI_TARGET - GI_LOSS);

// ---- Hemoglobin: slow relaxation to chronic-loss-adjusted target ----
double HB_TARGET = HB_BASE - HB_LOSS_SENS * (MENS_LOSS/30.0 + GI_LOSS);
dxdt_HB = KOUT_HB * (HB_TARGET - HB);

// ---- Serum sodium: DDAVP cumulative free-water retention risk ----
double NA_DROP = EMAX_NA_DROP * (1 - 0.6*FLUID_RESTRICT) * DDAVP_CP / (EC50_NA + DDAVP_CP);
dxdt_NA_SERUM = KOUT_NA * ( (NA_BASE - NA_DROP) - NA_SERUM );

// ---- Thrombotic-risk index: factor-concentrate overcorrection safety signal ----
double THROMB_TARGET = THROMB_GAIN * pow(VWF_RCO_TOTAL, 2) / (pow(EC50_THROMB, 2) + pow(VWF_RCO_TOTAL, 2));
dxdt_THROMB_RISK = KOUT_THROMB * (THROMB_TARGET - THROMB_RISK);

$CAPTURE VWF_AG_TOTAL VWF_RCO_TOTAL HMWM_EFFECTIVE ADAMTS13_ACT FVIII_C_TOTAL PLT_COUNT BLEED_SCORE MENS_LOSS GI_LOSS HB NA_SERUM THROMB_RISK WPB_STORE VWFPP DDAVP_CP RVWF_CONC PDVWF_CONC TXA_CP
'

vwd_mod <- mcode("vwd_qsp", vwd_code)

# =====================================================================
# Treatment scenarios (10)
#   All simulations: adult VWD patient, 30-day horizon unless noted.
#   Genotype set via param() overrides (VWFAG_BASE, VWFRCO_BASE, HMWM_BASE,
#   FVIII_BASE, TYPE2B_FLAG, TYPE2N_FLAG, ACQUIRED_FLAG).
# =====================================================================
WT <- 70  # kg, representative adult body weight

make_ev <- function(amt, ii, addl, cmt, time = 0) {
  ev(time = time, amt = amt, ii = ii, addl = addl, cmt = cmt)
}

scenarios <- list(
  "1_Type1_Mild_Untreated"          = NULL,
  "2_Type1_DDAVP_IV_single"         = make_ev(0.3 * WT, 0, 0, "DDAVP_DEPOT"),
  "3_Type1_DDAVP_Intranasal_Repeat" = make_ev(300, 12, 5, "DDAVP_DEPOT"),
  "4_Type2B_DDAVP_Contraindicated"  = make_ev(0.3 * WT, 0, 0, "DDAVP_DEPOT"),
  "5_Type3_Severe_rVWF_Vonvendi"    = make_ev(50 * WT, 24, 6, "RVWF_CP"),
  "6_Type3_Severe_PDConcentrate"    = make_ev(50 * WT, 24, 6, "PDVWF_CP"),
  "7_Menorrhagia_TXA_Hormonal"      = make_ev(1300, 8, 89, "TXA_DEPOT"),
  "8_Acquired_VWS_DDAVP_transient"  = make_ev(0.3 * WT, 0, 0, "DDAVP_DEPOT"),
  "9_Pregnancy_Type1_Peripartum"    = NULL,
  "10_Surgical_Major_PKguided_PDConc" = make_ev(60 * WT, 12, 20, "PDVWF_CP")
)

run_scenario <- function(name, ev_obj, params = list(), end = 720) {
  m <- vwd_mod
  if (length(params) > 0) m <- m %>% param(params)
  if (!is.null(ev_obj)) {
    out <- m %>% ev(ev_obj) %>% mrgsim(end = end, delta = 1) %>% as_tibble()
  } else {
    out <- m %>% mrgsim(end = end, delta = 1) %>% as_tibble()
  }
  out$scenario <- name
  out
}

# Genotype/parameter presets used by the example runs below:
p_type1   <- list(VWFAG_BASE = 35, VWFRCO_BASE = 30, FVIII_BASE = 42, HMWM_BASE = 0.85)
p_type2b  <- list(VWFAG_BASE = 45, VWFRCO_BASE = 22, FVIII_BASE = 55, HMWM_BASE = 0.55,
                   TYPE2B_FLAG = 1, PLT_BASE = 130)
p_type3   <- list(VWFAG_BASE = 2, VWFRCO_BASE = 1, FVIII_BASE = 6, HMWM_BASE = 0.02)
p_acquired<- list(VWFAG_BASE = 40, VWFRCO_BASE = 18, FVIII_BASE = 45, HMWM_BASE = 0.30,
                   ACQUIRED_FLAG = 1, ADAMTS13_BASE = 0.6)
p_preg    <- list(VWFAG_BASE = 35, VWFRCO_BASE = 30, FVIII_BASE = 42, HMWM_BASE = 0.85,
                   PREG_BOOST = 2.5)

# Example runs (uncomment to execute):
# results <- bind_rows(
#   run_scenario("1_Type1_Mild_Untreated", NULL, p_type1),
#   run_scenario("2_Type1_DDAVP_IV_single", scenarios[["2_Type1_DDAVP_IV_single"]], p_type1, end = 48),
#   run_scenario("3_Type1_DDAVP_Intranasal_Repeat", scenarios[["3_Type1_DDAVP_Intranasal_Repeat"]], p_type1, end = 168),
#   run_scenario("4_Type2B_DDAVP_Contraindicated", scenarios[["4_Type2B_DDAVP_Contraindicated"]], p_type2b, end = 48),
#   run_scenario("5_Type3_Severe_rVWF_Vonvendi", scenarios[["5_Type3_Severe_rVWF_Vonvendi"]], p_type3, end = 168),
#   run_scenario("6_Type3_Severe_PDConcentrate", scenarios[["6_Type3_Severe_PDConcentrate"]], p_type3, end = 168),
#   run_scenario("7_Menorrhagia_TXA_Hormonal", scenarios[["7_Menorrhagia_TXA_Hormonal"]],
#                c(p_type1, list(HORMONAL_ON = 1)), end = 720),
#   run_scenario("8_Acquired_VWS_DDAVP_transient", scenarios[["8_Acquired_VWS_DDAVP_transient"]], p_acquired, end = 48),
#   run_scenario("9_Pregnancy_Type1_Peripartum", NULL, p_preg, end = 720),
#   run_scenario("10_Surgical_Major_PKguided_PDConc", scenarios[["10_Surgical_Major_PKguided_PDConc"]], p_type3, end = 240)
# )
#
# ggplot(results, aes(time, VWF_RCO, color = scenario)) + geom_line(linewidth=1) +
#   labs(x = "Hour", y = "VWF:RCo (IU/dL)", title = "VWD: VWF:RCo trajectory by scenario")

# =====================================================================
# Calibration notes:
#  - VWFAG_BASE/VWFRCO_BASE Type-1 defaults (35/30 IU/dL) reflect the
#    ISTH/ASH/NHF/WFH 2021 diagnostic threshold zone (<30 IU/dL definite,
#    30-50 IU/dL "low VWF"); Type 3 defaults (<3 IU/dL) per James/Connell
#    Blood Adv 2021 diagnosis guideline (PMID 33570651).
#  - DDAVP EC50/EMAX_WPBREL calibrated so a single 0.3 mcg/kg IV dose in a
#    Type-1 patient produces a ~3-4x rise in VWF:Ag/RCo/FVIII:C peaking at
#    ~60-90 min, matching the pooled multicenter response-rate data in
#    Federici 2004 Blood (rise magnitude and duration ~6-8h before decline
#    toward KOUT_VWFAG/RCO half-life ~12h).
#  - Tachyphylaxis (K_WPBDEPLETE vs slow K_WPBREGEN) reproduces the
#    well-described blunted 2nd/3rd-dose response to repeated DDAVP over
#    24-48h dosing intervals (Mannucci 1997 Blood).
#  - rVWF (RVWF_CP) KE_RVWF corresponds to VWF:RCo t1/2 ~21h (Mannucci 2013
#    Blood phase-1 PK); because rVWF only raises the FVIII_STAB_TARGET (and
#    FVIII_C relaxes toward it at KOUT_FVIII, t1/2 ~9h) rather than adding
#    directly to plasma FVIII, the model reproduces the genuinely delayed/
#    gradual post-infusion FVIII:C climb reported for a VWF-only product in
#    the Gill 2015 Blood pivotal trial (vonicog alfa has no co-formulated
#    FVIII, unlike plasma-derived concentrate).
#  - Plasma-derived concentrate (PDVWF_CP) KE_PDVWF and VWF_FVIII_RATIO=2.4
#    approximate Humate-P-class VWF:RCo t1/2 ~12-20h / FVIII:C t1/2 ~8-12h
#    (Dobrkovska 1998 Haemophilia).
#  - Type-2B platelet clearance (K_PLT_CLEAR_2B, EC50_2B_CLEAR) illustrates
#    why DDAVP and high-dose replacement are used cautiously in Type 2B:
#    a VWF surge (endogenous or exogenous) transiently worsens
#    thrombocytopenia via spontaneous platelet aggregation/clearance.
#  - Menstrual/GI blood-loss and hemoglobin submodels are illustrative
#    physiologic proxies (not trial-calibrated point estimates) intended
#    to demonstrate qualitative direction and relative magnitude of
#    treatment effect (TXA, hormonal therapy, factor replacement).
# =====================================================================
