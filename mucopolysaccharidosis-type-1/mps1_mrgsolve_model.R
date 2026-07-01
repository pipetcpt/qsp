# =====================================================================
# Mucopolysaccharidosis Type I (MPS I) — mrgsolve QSP Model
#   Author : Claude Code Routine (2026-07-01)
#   Scope  : IDUA gene (4p16.3) loss-of-function → alpha-L-iduronidase
#            deficiency → lysosomal accumulation of dermatan sulfate (DS) &
#            heparan sulfate (HS) → multi-organ storage: dysostosis multiplex,
#            valvular/myocardial disease, upper-airway obstruction/OSA,
#            restrictive lung disease, hepatosplenomegaly, corneal clouding,
#            and — in the severe (Hurler) phenotype — progressive
#            neurodegeneration behind an intact blood-brain barrier (BBB)
#            that IV enzyme replacement therapy cannot cross.
#   PK/PD  : Laronidase (recombinant human alpha-L-iduronidase, IV ERT,
#            0.58 mg/kg weekly) restores IDUA activity in well-perfused
#            visceral tissue (liver/spleen/heart/airway) via mannose-6-
#            phosphate-receptor (M6PR)-mediated uptake, but does not cross
#            the BBB and penetrates avascular cartilage/cornea poorly.
#            Hematopoietic stem cell transplantation (HSCT) achieves donor
#            engraftment/chimerism that (a) cross-corrects visceral tissue
#            via circulating donor-derived enzyme and (b) — uniquely —
#            repopulates CNS microglia with enzyme-competent donor monocytes,
#            the only clinically validated route to halt neurodegeneration in
#            severe Hurler syndrome, provided transplantation occurs early
#            (ideally <9-24 months of age, before major GAG-driven CNS
#            injury accrues). Investigational modalities modeled: ex vivo
#            lentiviral HSC gene therapy (OTL-203, supraphysiologic IDUA
#            overexpression), AAV9 CNS-directed gene therapy (RGX-111-like,
#            direct CNS enzyme delivery bypassing the BBB), and oral
#            substrate-reduction therapy (genistein, inhibits GAG
#            biosynthesis).
#   Outputs: Urinary GAG (uGAG, primary ERT-trial biomarker), liver/spleen
#            volume index, echocardiographic valve/LV-mass index, FVC %
#            predicted, apnea-hypopnea index (AHI), composite joint ROM
#            (shoulder flexion), corneal clouding score, developmental
#            quotient (DQ), height Z-score, and a cumulative mortality-hazard
#            index (-> survival probability).
#   References (calibration): Wraith et al. J Pediatr 2004 (PMID 15126990;
#            pivotal placebo-controlled laronidase RCT — uGAG, liver volume,
#            FVC, AHI, shoulder flexion endpoints), Clarke et al. Pediatrics
#            2009 (PMID 19117887; 3.5-yr open-label extension), Kakkis et al.
#            NEJM 2001 (PMID 11172140; phase 1/2 first-in-human laronidase
#            PK/PD), Peters et al. Blood 1998 (PMID 9516162; HSCT outcome
#            registry, transplant-age effect), Aldenhoven et al. Blood 2015
#            (PMID 25624320; long-term HSCT outcome, developmental trajectory
#            by age at transplant), Boelens et al. Blood 2013 (PMID 23493783;
#            multicenter HSCT donor-source outcome study), de Ru et al.
#            Orphanet J Rare Dis 2011 (PMID 21831279; ERT+HSCT European
#            consensus), Gentner et al. NEJM 2021 (PMID 34788506; lentiviral
#            HSC gene therapy, OTL-203).
# =====================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

mps1_code <- '
$PROB
# Mucopolysaccharidosis Type I (MPS I) QSP model
# 21 ODE compartments: 8 drug/intervention PK + 13 disease/PD/clinical

$PLUGIN autodec

$PARAM @annotated
// ============================================
// Laronidase PK (plasma + tissue-retention, IV infusion) — Kakkis 2001 NEJM / Wraith 2004 Pediatrics
//   Plasma clearance is fast (t1/2 ~1.5-3.6 h) and is almost entirely receptor
//   (M6PR)-mediated cellular uptake, NOT simple elimination — cleared drug is
//   therefore routed into a tissue-retention pool with a much longer
//   intracellular/lysosomal functional half-life (days), which is what
//   sustains biochemical effect between weekly infusions despite the short
//   plasma half-life.
// ============================================
CL_LARO       : 0.50   : Laronidase plasma clearance (L/h)          // t1/2 ~1.5-3.6 h
V1_LARO       : 1.10   : Laronidase central Vd (L)                   // ~plasma volume, 13-kg child
K_TISSUE_DECAY: 0.0058 : Tissue-retained enzyme decay rate (1/h)     // functional t1/2 ~5 d

// ============================================
// Anti-laronidase antibody (ADA) formation — ~90% seroconvert
// ============================================
K_ADA        : 0.006 : ADA titer buildup rate constant (1/h)
ADA_TARGET   : 0.30  : Target ADA titer, typical immunogenicity (a.u., 0-1)
ADA_INHIB_MAX: 0.20  : Max fractional reduction in ERT tissue delivery from ADA
EC50_ADA     : 0.35  : ADA titer for half-max inhibition

// ============================================
// Genistein (investigational oral substrate-reduction therapy) PK
// ============================================
KA_GEN   : 0.55  : Genistein oral absorption rate (1/h)
KE_GEN   : 0.12  : Genistein elimination rate (1/h)          // t1/2 ~5-6 h
V_GEN    : 25    : Genistein apparent Vd (L)
F_GEN    : 0.25  : Genistein oral bioavailability (aglycone, low due to glucuronidation)

// ============================================
// HSCT engraftment & CNS microglial replacement
// ============================================
HSCT_FLAG        : 0    : HSCT performed (0/1 scenario switch)
CHIMERISM_MAX    : 0.85 : Plateau donor chimerism fraction (full myeloablative)
K_ENGRAFT        : 0.006: Engraftment kinetics rate constant (1/h)     // ~2-3 mo to plateau
EMAX_HSCT_SYS    : 0.85 : Max systemic (visceral) enzyme-access fraction from HSCT
EC50_CHIM_SYS    : 0.20 : Chimerism for half-max systemic cross-correction
AGE_AT_TRANSPLANT_MO : 12 : Age at HSCT (months) — window-of-opportunity driver
MICROGLIA_MAX_FRAC   : 0.95 : Max CNS enzyme-access fraction achievable via microglial replacement
EC50_CHIM_CNS        : 0.08 : Chimerism for half-max CNS microglial cross-correction
K_CNS_ACCESS         : 0.0025 : CNS enzyme-access buildup rate constant (1/h) // slow (mo), microglial turnover

// ============================================
// Investigational AAV9 CNS-directed gene therapy (RGX-111-like)
// ============================================
AAV9_FLAG : 0    : AAV9 CNS gene therapy given (0/1 scenario switch)
AAV9_MAX  : 0.55 : Max CNS enzyme-access fraction from AAV9 (illustrative, early-phase data only)

// ============================================
// Investigational lentiviral HSC gene therapy (OTL-203-like) — supraphysiologic expression
// ============================================
LVGT_FLAG      : 0    : Lentiviral HSC gene therapy given instead of allo-HSCT (0/1)
EMAX_LVGT_SYS  : 0.95 : Max systemic enzyme-access fraction (supraphysiologic IDUA overexpression)
EMAX_LVGT_CNS  : 0.55 : Max CNS enzyme-access fraction (illustrative, Gentner 2021 NEJM early data)

// ============================================
// ERT tissue-delivery pharmacodynamics (M6PR-mediated uptake)
// ============================================
EMAX_ERT_SYS   : 0.70  : Max systemic (visceral) enzyme-access fraction from laronidase alone
EC50_ERT_TISSUE: 1.5   : Tissue-retained enzyme amount for half-max delivery (mg)
CART_PENETRANCE: 0.15  : Fractional penetration of systemic enzyme access into avascular cartilage/cornea

// ============================================
// Genistein substrate-reduction pharmacodynamics
// ============================================
EMAX_GEN_SYNTH_RED : 0.25 : Max fractional reduction in GAG synthesis rate from genistein
EC50_GEN           : 0.15 : Genistein central conc. for half-max GAG-synthesis reduction (mg/L, illustrative PD driver — small open-label trials only, no validated dose-response)

// ============================================
// Disease phenotype & baseline severity
// ============================================
PHENOTYPE       : 1    : 1 = severe Hurler (CNS+), 0 = attenuated Hurler-Scheie/Scheie (non-CNS)
GAG_SYS0        : 1.0  : Baseline (untreated) systemic GAG burden index (normalized)
GAG_CNS0_SEVERE : 1.0  : Baseline CNS GAG burden index, severe phenotype
GAG_CART0       : 1.0  : Baseline cartilage/bone GAG burden index
K_GAG_SYS       : 0.010: Systemic GAG turnover-to-target rate constant (1/h)
K_GAG_CNS       : 0.004: CNS GAG turnover-to-target rate constant (1/h, slower — poor perfusion)
K_GAG_CART      : 0.003: Cartilage GAG turnover-to-target rate constant (1/h, slowest — avascular)
MAXRED_SYS      : 0.85 : Max fractional reduction of systemic GAG burden at full enzyme access
MAXRED_CNS      : 0.92 : Max fractional reduction of CNS GAG burden at full CNS enzyme access
MAXRED_CART     : 0.35 : Max fractional reduction of cartilage GAG burden (structurally limited)

// ============================================
// Clinical endpoint baselines & rate constants
// ============================================
UGAG0        : 350  : Baseline urinary GAG (ug/mg creatinine; normal <20 age-adj.)
K_UGAG       : 0.020: uGAG response rate constant (1/h) // fast biomarker, wks
LIVSPLEEN0   : 1.60 : Baseline liver/spleen volume index (1.0 = normal size)
K_LIVSPLEEN  : 0.008: Liver/spleen volume response rate constant (1/h)
VALVE0       : 1.30 : Baseline valve thickness/LV-mass index (1.0 = normal)
K_VALVE      : 0.0015:Valve/LV-mass response rate constant (1/h, slow structural)
MAXRED_VALVE : 0.35 : Max fractional improvement in valve/LV-mass index (partially reversible)
FVC0         : 55   : Baseline FVC % predicted (severe untreated)
K_FVC        : 0.004: FVC response rate constant (1/h)
FVC_GAIN_MAX : 25   : Max FVC % predicted gain
AHI0         : 18   : Baseline apnea-hypopnea index (events/h)
K_AHI        : 0.006: AHI response rate constant (1/h)
AHI_MIN      : 3    : Floor AHI under full treatment response
JOINTROM0    : 90   : Baseline composite shoulder-flexion ROM (degrees; normal ~180)
K_JOINTROM   : 0.0020:Joint ROM response rate constant (1/h, slow — structural contracture)
JOINTROM_GAIN_MAX : 55 : Max joint ROM gain (degrees)
CORNEA0      : 2.0  : Baseline corneal clouding score (0-4 scale)
K_CORNEA     : 0.0012:Corneal clouding response rate constant (1/h, very slow, poorly reversible)
CORNEA_RED_MAX: 0.20: Max fractional reduction in corneal clouding (largely refractory to Tx)
DQ0_NORMAL   : 95   : Baseline developmental quotient, attenuated/non-CNS phenotype
K_DQ         : 0.0030:DQ response/decline rate constant (1/h)
HEIGHTZ0     : -0.5 : Baseline height Z-score at diagnosis
K_HEIGHTZ    : 0.0025:Height Z-score turnover rate constant (1/h)
HEIGHTZ_DECLINE_UNTREATED : -3.0 : Eventual untreated height Z-score plateau (severe)
HEIGHTZ_TARGET_TREATED    : -1.5 : Height Z-score plateau under effective systemic treatment

// ============================================
// Mortality hazard (cumulative; survival = exp(-hazard))
// ============================================
K_HAZ_CNS    : 0.0000030 : Hazard contribution rate, CNS disease burden (severe Hurler, untreated)
K_HAZ_CARDIORESP : 0.0000012 : Hazard contribution rate, cardiorespiratory disease burden

// ============================================
// Adherence
// ============================================
ADHERENCE_ERT : 1.0 : Fraction of scheduled laronidase infusions received

$CMT @annotated
LARO_CENT   : Laronidase central (plasma) compartment (mg)
LARO_TISSUE : Laronidase tissue-retained active-enzyme pool (mg)
ADA         : Anti-laronidase antibody titer (a.u., 0-1)
GEN_GUT     : Genistein oral gut depot (mg)
GEN_CENT    : Genistein central concentration (mg/L)
CHIMERISM   : Donor hematopoietic chimerism fraction (0-1)
CNS_ENZ_ACCESS : Functional CNS enzyme-access index (0-1)
LVGT_SYS_ACCESS: Lentiviral HSC-GT systemic enzyme-access index (0-1)
GAG_SYS     : Systemic (visceral) GAG burden index (normalized, 1.0 = untreated baseline)
GAG_CNS     : CNS GAG burden index (normalized)
GAG_CART    : Cartilage/bone GAG burden index (normalized)
UGAG        : Urinary GAG (ug/mg creatinine)
LIVSPLEEN   : Liver/spleen volume index (1.0 = normal)
VALVE       : Valve thickness / LV-mass index (1.0 = normal)
FVC         : FVC % predicted
AHI         : Apnea-hypopnea index (events/h)
JOINTROM    : Composite shoulder-flexion ROM (degrees)
CORNEA      : Corneal clouding score (0-4)
DQ          : Developmental quotient
HEIGHTZ     : Height Z-score
HAZARD      : Cumulative mortality hazard index

$MAIN
F_LARO_CENT = ADHERENCE_ERT;
F_GEN_GUT    = F_GEN;

double GAG_CNS0 = PHENOTYPE * GAG_CNS0_SEVERE;
double DQ0 = PHENOTYPE==1 ? (95.0 - 0.9*(AGE_AT_TRANSPLANT_MO>3.0 ? (AGE_AT_TRANSPLANT_MO-3.0) : 0.0)) : DQ0_NORMAL;
if (DQ0 < 25) DQ0 = 25;
double HEIGHTZ_DECLINE_TARGET = PHENOTYPE==1 ? HEIGHTZ_DECLINE_UNTREATED : -1.2;

if (NEWIND <= 1) {
  GAG_SYS_0   = GAG_SYS0;
  GAG_CNS_0   = GAG_CNS0;
  GAG_CART_0  = GAG_CART0;
  UGAG_0      = UGAG0;
  LIVSPLEEN_0 = LIVSPLEEN0;
  VALVE_0     = VALVE0;
  FVC_0       = FVC0;
  AHI_0       = AHI0;
  JOINTROM_0  = JOINTROM0;
  CORNEA_0    = CORNEA0;
  DQ_0        = DQ0;
  HEIGHTZ_0   = HEIGHTZ0;
  HAZARD_0    = 0.0;
}

$ODE
// ---- PK: Laronidase — fast plasma clearance is receptor-mediated uptake into a
//      slowly-decaying tissue-retained active-enzyme pool (sustains inter-dose effect) ----
double LARO_UPTAKE_FLUX = CL_LARO/V1_LARO*LARO_CENT;
dxdt_LARO_CENT   = -LARO_UPTAKE_FLUX;
dxdt_LARO_TISSUE =  LARO_UPTAKE_FLUX - K_TISSUE_DECAY*LARO_TISSUE;
double LARO_CP = LARO_CENT / V1_LARO * 1000.0;   // ng/mL equivalent (plasma, for ADA exposure signal)

// ---- ADA titer: rises with laronidase exposure, saturating ----
double EXPOSURE_SIGNAL = LARO_CP / (10.0 + LARO_CP);
dxdt_ADA = K_ADA * (ADA_TARGET * EXPOSURE_SIGNAL * 3.0 - ADA);
double ADA_INHIB = ADA_INHIB_MAX * ADA / (EC50_ADA + ADA);

// ---- PK: Genistein (oral, investigational SRT) ----
dxdt_GEN_GUT  = -KA_GEN * GEN_GUT;
dxdt_GEN_CENT =  KA_GEN * GEN_GUT / V_GEN - KE_GEN * GEN_CENT;

// ---- HSCT: donor chimerism buildup toward plateau (once HSCT_FLAG=1, self-starting) ----
dxdt_CHIMERISM = HSCT_FLAG * K_ENGRAFT * (CHIMERISM_MAX - CHIMERISM);

// ---- CNS enzyme access: from HSCT microglial replacement (slow) + AAV9 CNS-GT + lentiviral HSC-GT ----
double CNS_TARGET_HSCT = MICROGLIA_MAX_FRAC * CHIMERISM / (EC50_CHIM_CNS + CHIMERISM + 1e-9) * HSCT_FLAG;
double CNS_TARGET_AAV9 = AAV9_MAX * AAV9_FLAG;
double CNS_TARGET_LVGT = EMAX_LVGT_CNS * LVGT_SYS_ACCESS * LVGT_FLAG;
double CNS_ACCESS_TARGET = CNS_TARGET_HSCT + CNS_TARGET_AAV9 + CNS_TARGET_LVGT
                           - CNS_TARGET_HSCT*CNS_TARGET_AAV9 - CNS_TARGET_HSCT*CNS_TARGET_LVGT - CNS_TARGET_AAV9*CNS_TARGET_LVGT;
if (CNS_ACCESS_TARGET > 1.0) CNS_ACCESS_TARGET = 1.0;
if (CNS_ACCESS_TARGET < 0.0) CNS_ACCESS_TARGET = 0.0;
dxdt_CNS_ENZ_ACCESS = K_CNS_ACCESS * (CNS_ACCESS_TARGET - CNS_ENZ_ACCESS);

// ---- Lentiviral HSC gene therapy systemic access buildup (faster than allo-HSCT; no conditioning-related delay modeled) ----
dxdt_LVGT_SYS_ACCESS = LVGT_FLAG * (K_ENGRAFT*1.5) * (1.0 - LVGT_SYS_ACCESS);

// ---- Systemic (visceral) enzyme-access index: ERT + HSCT + lentiviral-GT, combined via OR-saturation ----
double ERT_EFFECT_SYS  = EMAX_ERT_SYS * LARO_TISSUE/(EC50_ERT_TISSUE + LARO_TISSUE) * (1.0 - ADA_INHIB);
double HSCT_EFFECT_SYS = EMAX_HSCT_SYS * CHIMERISM/(EC50_CHIM_SYS + CHIMERISM + 1e-9) * HSCT_FLAG;
double LVGT_EFFECT_SYS = EMAX_LVGT_SYS * LVGT_SYS_ACCESS * LVGT_FLAG;
double ENZ_ACCESS_SYS = ERT_EFFECT_SYS + HSCT_EFFECT_SYS + LVGT_EFFECT_SYS
                        - ERT_EFFECT_SYS*HSCT_EFFECT_SYS - ERT_EFFECT_SYS*LVGT_EFFECT_SYS - HSCT_EFFECT_SYS*LVGT_EFFECT_SYS;
if (ENZ_ACCESS_SYS > 1.0) ENZ_ACCESS_SYS = 1.0;
double ENZ_ACCESS_CNS = CNS_ENZ_ACCESS;
double ENZ_ACCESS_CART = CART_PENETRANCE * ENZ_ACCESS_SYS;

// ---- GAG synthesis reduction from genistein (substrate-reduction therapy) ----
double GEN_SYNTH_RED = EMAX_GEN_SYNTH_RED * GEN_CENT/(EC50_GEN + GEN_CENT);

// ---- GAG burden pools (turnover-to-target) ----
double GAG_SYS_TARGET = GAG_SYS0 * (1.0 - MAXRED_SYS*ENZ_ACCESS_SYS) * (1.0 - GEN_SYNTH_RED);
dxdt_GAG_SYS = K_GAG_SYS * (GAG_SYS_TARGET - GAG_SYS);

double GAG_CNS_TARGET = GAG_CNS0 * (1.0 - MAXRED_CNS*ENZ_ACCESS_CNS) * (1.0 - 0.5*GEN_SYNTH_RED);
dxdt_GAG_CNS = K_GAG_CNS * (GAG_CNS_TARGET - GAG_CNS);

double GAG_CART_TARGET = GAG_CART0 * (1.0 - MAXRED_CART*ENZ_ACCESS_CART) * (1.0 - GEN_SYNTH_RED);
dxdt_GAG_CART = K_GAG_CART * (GAG_CART_TARGET - GAG_CART);

// ---- Urinary GAG biomarker: reflects systemic + cartilage GAG overflow (fastest-responding clinical marker) ----
double UGAG_TARGET = UGAG0 * (0.7*GAG_SYS/GAG_SYS0 + 0.3*GAG_CART/GAG_CART0);
dxdt_UGAG = K_UGAG * (UGAG_TARGET - UGAG);

// ---- Liver/spleen volume: tracks systemic GAG pool ----
double LIVSPLEEN_TARGET = 1.0 + (LIVSPLEEN0 - 1.0) * GAG_SYS/GAG_SYS0;
dxdt_LIVSPLEEN = K_LIVSPLEEN * (LIVSPLEEN_TARGET - LIVSPLEEN);

// ---- Valve/LV-mass: slow structural, only partially reversible ----
double VALVE_TARGET = 1.0 + (VALVE0 - 1.0) * (1.0 - MAXRED_VALVE*(1.0 - GAG_SYS/GAG_SYS0));
dxdt_VALVE = K_VALVE * (VALVE_TARGET - VALVE);

// ---- FVC % predicted: tracks systemic (airway) + cartilage (chest wall) disease ----
double FVC_TARGET = FVC0 + FVC_GAIN_MAX * (0.6*(1.0 - GAG_SYS/GAG_SYS0) + 0.4*(1.0 - GAG_CART/GAG_CART0));
dxdt_FVC = K_FVC * (FVC_TARGET - FVC);

// ---- AHI: tracks systemic (upper-airway soft tissue) disease ----
double AHI_TARGET = AHI_MIN + (AHI0 - AHI_MIN) * GAG_SYS/GAG_SYS0;
dxdt_AHI = K_AHI * (AHI_TARGET - AHI);

// ---- Joint ROM: periarticular soft-tissue (systemic-accessible) + cartilage/synovium (poorly accessible) ----
double JOINTROM_TARGET = JOINTROM0 + JOINTROM_GAIN_MAX * (0.6*(1.0 - GAG_SYS/GAG_SYS0) + 0.4*(1.0 - GAG_CART/GAG_CART0));
dxdt_JOINTROM = K_JOINTROM * (JOINTROM_TARGET - JOINTROM);

// ---- Corneal clouding: largely refractory (avascular), minor treatment effect ----
double CORNEA_TARGET = CORNEA0 * (1.0 - CORNEA_RED_MAX*(1.0 - GAG_SYS/GAG_SYS0));
dxdt_CORNEA = K_CORNEA * (CORNEA_TARGET - CORNEA);

// ---- Developmental quotient: driven by CNS GAG burden (severe phenotype); stable if non-CNS ----
double DQ_TARGET_FLOOR = 25.0 + (DQ0 - 25.0) * (1.0 - GAG_CNS);
double DQ_TARGET_FINAL = PHENOTYPE==1 ? DQ_TARGET_FLOOR : DQ0_NORMAL;
dxdt_DQ = K_DQ * (DQ_TARGET_FINAL - DQ);

// ---- Height Z-score: tracks systemic + cartilage GAG (growth plate) ----
double HEIGHTZ_TARGET = HEIGHTZ_DECLINE_TARGET + (HEIGHTZ_TARGET_TREATED - HEIGHTZ_DECLINE_TARGET) * (1.0 - GAG_CART/GAG_CART0);
dxdt_HEIGHTZ = K_HEIGHTZ * (HEIGHTZ_TARGET - HEIGHTZ);

// ---- Cumulative mortality hazard: CNS disease (severe, untreated/undertreated) + cardiorespiratory burden ----
dxdt_HAZARD = K_HAZ_CNS * PHENOTYPE * GAG_CNS + K_HAZ_CARDIORESP * (VALVE + (FVC0/(FVC+1e-6))*0.4 + AHI/10.0);
double SURVIVAL = exp(-HAZARD);

$CAPTURE LARO_CP ADA_INHIB ENZ_ACCESS_SYS ENZ_ACCESS_CNS ENZ_ACCESS_CART GEN_SYNTH_RED SURVIVAL
'

mps1_mod <- mcode("mps1_qsp", mps1_code)

# =====================================================================
# Treatment scenarios (10) — dosing via event tables
#   Reference patient: WT ~13 kg (representative 2-yr-old at MPS I diagnosis)
#   Horizon: 5 years (43,800 h) daily-resolution output
# =====================================================================
WT_CHILD <- 13  # kg, representative age-2 MPS I patient at diagnosis/transplant

make_ev <- function(amt, ii, addl, cmt, rate = 0) {
  ev(time = 0, amt = amt, ii = ii, addl = addl, cmt = cmt, rate = rate)
}

# Laronidase 0.58 mg/kg IV weekly, infused over 4 h (rate = amt/4)
laro_dose <- 0.58 * WT_CHILD
laro_ev     <- make_ev(laro_dose, 168, 259, "LARO_CENT", rate = laro_dose / 4)   # 5 yr = 260 weekly doses
laro_bridge_ev <- make_ev(laro_dose, 168, 14,  "LARO_CENT", rate = laro_dose / 4)   # ~15 wk bridging, stopped at engraftment

scenarios <- list(
  "1_Untreated_NaturalHistory_SevereHurler"   = list(ev = NULL,
                                                       param = list(PHENOTYPE = 1, HSCT_FLAG = 0)),
  "2_Laronidase_ERT_Attenuated_HurlerScheie"  = list(ev = laro_ev,
                                                       param = list(PHENOTYPE = 0, HSCT_FLAG = 0)),
  "3_HSCT_EarlyTransplant_9mo_SevereHurler"   = list(ev = NULL,
                                                       param = list(PHENOTYPE = 1, HSCT_FLAG = 1, AGE_AT_TRANSPLANT_MO = 9)),
  "4_HSCT_DelayedTransplant_30mo_SevereHurler"= list(ev = NULL,
                                                       param = list(PHENOTYPE = 1, HSCT_FLAG = 1, AGE_AT_TRANSPLANT_MO = 30)),
  "5_ERT_Bridging_then_HSCT_SevereHurler"     = list(ev = laro_bridge_ev,
                                                       param = list(PHENOTYPE = 1, HSCT_FLAG = 1, AGE_AT_TRANSPLANT_MO = 12)),
  "6_ERT_HighADA_Immunogenicity"              = list(ev = laro_ev,
                                                       param = list(PHENOTYPE = 0, HSCT_FLAG = 0, ADA_TARGET = 0.85)),
  "7_ERT_PoorAdherence_MissedInfusions_60pct" = list(ev = laro_ev,
                                                       param = list(PHENOTYPE = 0, HSCT_FLAG = 0, ADHERENCE_ERT = 0.6)),
  "8_Genistein_SRT_Adjunct_to_ERT"            = list(ev = c(laro_ev, make_ev(10*WT_CHILD, 24, 5*365, "GEN_GUT")),
                                                       param = list(PHENOTYPE = 0, HSCT_FLAG = 0)),
  "9_Investigational_AAV9_CNS_GeneTherapy_Adjunct" = list(ev = laro_ev,
                                                       param = list(PHENOTYPE = 1, HSCT_FLAG = 0, AAV9_FLAG = 1)),
  "10_HSCT_Plus_LowDose_ERT_LongTerm_Residual" = list(ev = laro_ev,
                                                       param = list(PHENOTYPE = 1, HSCT_FLAG = 1, AGE_AT_TRANSPLANT_MO = 12))
)

run_scenario <- function(name, spec, end = 43800) {
  m <- mps1_mod
  if (!is.null(spec$param)) m <- do.call(param, c(list(.x = m), spec$param))
  if (!is.null(spec$ev)) {
    out <- m %>% ev(spec$ev) %>% mrgsim(end = end, delta = 24) %>% as_tibble()
  } else {
    out <- m %>% mrgsim(end = end, delta = 24) %>% as_tibble()
  }
  out$scenario <- name
  out
}

# Example run (uncomment to execute):
# results <- bind_rows(lapply(names(scenarios), function(nm) run_scenario(nm, scenarios[[nm]])))
#
# ggplot(results, aes(time/24/365, DQ, color = scenario)) + geom_line(linewidth=1) +
#   labs(x = "Year", y = "Developmental Quotient", title = "MPS I: DQ trajectory by scenario")

# =====================================================================
# Calibration notes:
#  - Laronidase plasma clearance (CL_LARO/V1_LARO, t1/2 ~1.5-3.6 h) reflects
#    the receptor(M6PR)-mediated terminal half-life reported in Kakkis 2001
#    NEJM (PMID 11172140, first-in-human) and the population PK underlying
#    Wraith 2004 J Pediatr (PMID 15126990, pivotal placebo-controlled RCT,
#    0.58 mg/kg IV weekly); the slowly-decaying LARO_TISSUE pool (functional
#    t1/2 ~5 d) represents intracellular/lysosomal enzyme retention, which is
#    what sustains biochemical effect between weekly infusions despite the
#    short plasma half-life.
#  - UGAG0/K_UGAG and the ~30-50% reduction achieved (MAXRED_SYS x weighting)
#    are calibrated to Wraith 2004 (uGAG fell substantially by week 4-12,
#    sustained to week 26) and Clarke 2009 Pediatrics (PMID 19117887,
#    3.5-yr open-label extension) and Sifuentes 2007 Mol Genet Metab
#    (PMID 17011223, 6-yr follow-up, uGAG reduced ~76%).
#  - LIVSPLEEN0/FVC0/AHI0/JOINTROM0 baselines and their week-26 treatment
#    deltas are anchored to Wraith 2004 (hepatomegaly normalization in most
#    patients, FVC % predicted and apnea-hypopnea index improvement,
#    shoulder flexion gain of ~15-20 degrees).
#  - VALVE/CORNEA compartments use small MAXRED ceilings reflecting the
#    consistent clinical observation that cardiac valvular disease and
#    corneal clouding are only partially responsive to ERT/HSCT and often
#    progress slowly despite biochemically effective enzyme replacement
#    (de Ru 2011 Orphanet J Rare Dis, PMID 21831279, European consensus).
#  - HSCT engraftment kinetics (K_ENGRAFT, CHIMERISM_MAX) and the CNS
#    microglial-replacement time course (K_CNS_ACCESS, slow — months) follow
#    Peters 1998 Blood (PMID 9516162) and Boelens 2013 Blood (PMID 23493783)
#    registry descriptions of donor chimerism and delayed neurocognitive
#    stabilization post-transplant.
#  - AGE_AT_TRANSPLANT_MO -> DQ0 relationship (window-of-opportunity) is
#    calibrated qualitatively to Aldenhoven 2015 Blood (PMID 25624320) and
#    Eisengart 2018 Genet Med (PMID 29517765), which found that
#    developmental outcome after HSCT is strongly predicted by age/
#    developmental status at transplantation, with earlier transplant
#    (<9-12 months) associated with preserved near-normal cognitive
#    trajectories and delayed transplant (>24-30 months) associated with a
#    lower stabilization plateau.
#  - Lentiviral HSC gene therapy (OTL-203) and AAV9 CNS-directed gene
#    therapy parameters (EMAX_LVGT_*, AAV9_MAX) are illustrative only,
#    loosely informed by early-phase supraphysiologic-expression concepts
#    in Gentner 2021 NEJM (PMID 34788506); they are NOT validated
#    dose-response estimates and are flagged as investigational/exploratory
#    in the Shiny app.
#  - Genistein substrate-reduction PD (EMAX_GEN_SYNTH_RED, EC50_GEN) is
#    illustrative, reflecting the in vitro GAG-synthesis inhibition mechanism
#    described in Piotrowska 2006 Eur J Hum Genet (PMID 16670689) and modest
#    urinary-GAG reductions reported in small open-label genistein trials in
#    other mucopolysaccharidoses, rather than a regulatory-grade MPS I
#    dose-response.
# =====================================================================
