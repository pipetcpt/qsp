## ===========================================================================
##  mmd_mrgsolve_model.R  —  Moyamoya Disease (MMD) QSP model for mrgsolve
## ===========================================================================
##
##  THESIS
##  ------
##  In moyamoya disease the MEASURED quantity is cerebral blood flow.  The
##  STATE variable is how much dilatory range the cortical arteriole has left.
##  One hemisphere is written as a two-node resistive network:
##
##      Pa --[ R_ica(stenosis) ]--+
##      Pa --[ R_moya           ]--+--> P_A --[ R_artA ]--> Pv
##      Pa --[ R_pva            ]--+        (autoregulator, has a FLOOR)
##      Pa --[ R_byp (surgery)  ]--+
##                                 |
##                          [ R_coll ]        <-- donated by the PCA
##                                 |
##      Pa --[ R_pca ]---------> P_B --[ R_artB ]--> Pv
##
##  Nodal balance (Pv = ICP):
##      (gi+gm+gv+gb)(Pa-P_A) + gc(P_B-P_A) = gA(P_A-Pv)
##                 gp(Pa-P_B)               = gB(P_B-Pv) + gc(P_B-P_A)
##
##  Autoregulation picks gA and gB to meet demand, CLAMPED to
##  [1/R_art_max_eff, 1/R_art_min].  Everything else in the model is a
##  consequence of that clamp.  In particular:
##
##   *  A CRITICAL INLET CONDUCTANCE gS* = 2.201 mL/min/mmHg exists.  Above it
##      CBF is defended and dCBF/dMAP = 0; below it the arteriole is on its
##      floor, CBF becomes pressure-passive, and the measured CO2/acetazolamide
##      response turns NEGATIVE.  All three switch at the same point.
##   *  THE PENUMBRAL THRESHOLD IS DERIVED, not fitted:
##          CBF_crit = CMRO2 / (CaO2 * OEF_max)
##      = 19.71 mL/100g/min at Hb 15 g/dL (the textbook 20), and
##      = 36.96 mL/100g/min at Hb 8 g/dL.  Sickle-cell moyamoya therefore
##      infarcts at flows a normal brain tolerates, and TRANSFUSION — which
##      touches no vessel — is the intervention that moves the threshold.
##   *  ACETAZOLAMIDE AND PaCO2 ARE NOT THE SAME PROBE.  Acetazolamide acts on
##      the arteriole only; PaCO2 acts on the arteriole AND on the pial
##      collateral CONDUITS.  In a normal brain the distinction is invisible;
##      in moyamoya the conduits are the circulation.
##   *  TWO COLLATERAL ROUTES WITH OPPOSITE SAFETY PROFILES.  The
##      leptomeningeal route (g_coll) is donated by the PCA and is safe.  The
##      periventricular / choroidal route (g_pva) is thin-walled and is the
##      haemorrhage source.  PCA involvement removes the safe donor AND forces
##      the dangerous route — one flag, two opposite-signed consequences.
##   *  RNF213 ENTERS TWICE WITH OPPOSITE SIGN: it drives the intimal lesion
##      and it CAPS the collateral ceiling.
##
##  STRUCTURE:  35 ODE compartments (lesion, two collateral routes, surgical
##  graft, arteriolar remodelling, angiogenic signalling, injury, haemorrhage,
##  and a 7-agent PK layer) + a fully algebraic quasi-steady haemodynamic
##  layer solved by fixed-point iteration inside $ODE.
##
##  CALIBRATION / PROVENANCE
##  ------------------------
##  Anchors used to SET parameters (see mmd_references.md for PMIDs):
##    - Normal cortical CBF 50 mL/100g/min, CMRO2 3.3, OEF 0.33; large-artery
##      share of cerebrovascular resistance ~25%.
##    - Acetazolamide 1 g IV raises normal CBF 30-40%  -> AZ_EMAX 0.45.
##    - CO2 reactivity ~3.5%/mmHg PaCO2                -> K_CO2 0.035.
##    - Symptomatic adult MMD affected-hemisphere CBF ~15-25% below normal
##      with OEF 0.42-0.50 (15O-PET misery perfusion) -> collateral ceilings.
##    - JAM Trial CONSERVATIVE arm 5-y rebleeding 31.6% -> HEM_HAZ0 0.414.
##      The bypass HR is then a PREDICTION, not a fit (see $PROB note).
##    - Post-bypass symptomatic hyperperfusion peaks day 2-7, resolves over
##      2-3 weeks -> K_REM_OFF 0.05 /d (tau ~ 20 d).
##    - Direct STA-MCA patency > 95%; peri-operative stroke 4-7%.
##
##  Every equation below is a line-by-line translation of
##  mmd_reference_model.py, whose printed output (mmd_reference_output.txt)
##  is the numerical truth this file is checked against.
##
##  Author: QSP Disease Model Library (Claude Code Routine), 2026-08-06
##  EDUCATIONAL / RESEARCH USE ONLY — not for clinical decision-making.
## ===========================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)

code <- '
$PROB
# Moyamoya Disease QSP model
Reserve, not flow, is the state variable.  35 ODEs + an algebraic
quasi-steady two-node haemodynamic network with a dilatory floor.

$PARAM @annotated
// ---- systemic / geometry -------------------------------------------------
MAP0      :  90.0  : baseline mean arterial pressure (mmHg)
MAP_REF   :  90.0  : fixed normalising pressure for wall stress (mmHg)
ICP       :  10.0  : venous outflow pressure (mmHg)
MASS_A    : 250.0  : MCA (anterior) territory mass (g)
MASS_B    : 120.0  : PCA (posterior) territory mass (g)
HB        :  15.0  : haemoglobin (g/dL)
SAO2      :   0.98 : arterial oxygen saturation (fraction)
CMRO2_0   :   3.30 : cerebral metabolic rate for O2 (mL O2/100g/min)
OEF_MAX   :   0.85 : hard ceiling on oxygen extraction fraction
OEF_BASE  :   0.335: resting oxygen extraction fraction
CBF_A0    :  50.0  : normal CBF territory A (mL/100g/min)
CBF_B0    :  50.0  : normal CBF territory B (mL/100g/min)
FRAC_PROX :   0.25 : large-artery share of territory resistance (fraction)
DIL_MAX   :   2.50 : maximal arteriolar dilatation (R_art0/R_art_min)
CON_MAX   :   2.20 : maximal arteriolar constriction (R_art_max/R_art0)
GCOLL0    :   0.050: native leptomeningeal conductance (mL/min/mmHg)
GMOYA0    :   0.004: native basal perforator conductance (mL/min/mmHg)
GPVA0     :   0.002: native periventricular conductance (mL/min/mmHg)

// ---- RNF213 / intimal hyperplasia (the lesion) ---------------------------
RNF       :   1.0  : RNF213 R4810K dose (0 wt, 1 het, 1.6 hom)
SEC       :   0.0  : RNF213-independent intimal drive (quasi-moyamoya)
K_PDGF    :   0.90 : PDGF-BB production gain (/d)
K_PDGF_D  :   0.60 : PDGF-BB decay (/d)
K_SMC     :   2.6e-4 : intimal SMC accumulation rate constant (/d)
SMC_MAX   :   1.0  : intimal SMC saturation (dimensionless)
STEN_MAX  :   0.945: asymptotic fractional DIAMETER stenosis
K_STATIN  :   0.28 : max fractional suppression of K_SMC by statin

// ---- collateral / angiogenesis -------------------------------------------
ANGIO     :   0.58 : angiogenic capacity (child ~1.0, adult ~0.35-0.6)
K_HIF_ON  :   6.0  : HIF-1a induction gain (/d)
K_HIF_OFF :   3.0  : HIF-1a decay (/d)
K_VEGF    :   2.0  : VEGF-A production gain (/d)
K_VEGF_D  :   1.4  : VEGF-A decay (/d)
K_MMP     :   1.5  : MMP-9 production gain (/d)
K_MMP_D   :   1.0  : MMP-9 decay (/d)
K_MOYA    :   0.018: basal perforator recruitment (/d)
GMOYA_CAP :   3.40 : basal perforator ceiling before RNF213 penalty
RNF_CAP   :   0.42 : collateral ceiling penalty per unit RNF (fraction)
K_MOYA_REG:   0.030: perforator regression on VEGF withdrawal (/d)
K_COLL    :   0.009: leptomeningeal recruitment (/d)
GCOLL_CAP :   1.00 : leptomeningeal ceiling before penalties
K_COLL_REG:   0.012: leptomeningeal regression (/d)
VEGF_REG  :   0.60 : VEGF below which the collateral bed is pruned

// ---- the periventricular / choroidal anastomosis: THE BLEEDER ------------
PCA_INV   :   0.0  : posterior cerebral artery involvement (0-1)
PCA_STEN  :   0.25 : PCA diameter stenosis at PCA_INV = 1
PCA_COLL  :   0.60 : fractional loss of leptomeningeal ceiling at PCA_INV = 1
PCA_PVA   :   2.20 : gain on periventricular ceiling at PCA_INV = 1
GPVA_CAP  :   0.55 : periventricular ceiling before penalties
K_PVA     :   0.014: periventricular recruitment (/d)
K_PVA_REG :   0.030: periventricular regression (/d)
GPVA_RREF :   0.050: conductance above which the perforator counts as dilated
GMOYA_RREF:   0.20 : same, for the basal perforator bed
Q_LOAD_REF: 125.0  : native design flow of territory A (mL/min)

// ---- surgery -------------------------------------------------------------
GBYP_DIR  :   1.05 : mature direct STA-MCA conductance (mL/min/mmHg)
TAU_BYP   :   6.0  : direct graft maturation time constant (d)
GBYP_IND  :   1.15 : indirect (EDAS/EMS) construct ceiling
K_BYP_IND :   0.055: VEGF- and ANGIO-gated indirect ingrowth (/d)
SURG_KIND :   0    : 0 none, 1 direct, 2 indirect, 3 combined
SURG_T    : 1e9    : day of operation
SURG_SCALE:   1.0  : construct size multiplier

// ---- arteriolar remodelling (the hyperperfusion mechanism) ---------------
K_REM_ON  :   0.075: gain of chronic-dilatation remodelling (/d)
K_REM_OFF :   0.050: recovery of constrictive range (/d)
// VASOPARALYSIS: REM_CON = DIL_MAX*CON_MAX = 5.50 collapses the WHOLE
// [g_con, g_dil] interval at REMOD = 1, which is what the hyperperfusion
// literature means by the word.
REM_CON   :   5.50 : g_art constrictive-ceiling multiplier at REMOD = 1
F_FRAC    :   0.10 : peri-anastomotic share of territory A (fraction)
G_LEAK    :   0.80 : pial/M3 coupling, focal cortex <-> territory (mL/min/mmHg)
TAU_ADAPT :  22.0  : barrier / autoregulatory re-adaptation time constant (d)

// ---- watershed geometry --------------------------------------------------
WS_FRAC   :   0.30 : borderzone fraction of territory A
WS_K      :   0.55 : max watershed penalty at fully collateral supply
W_BYP     :   0.80 : a direct graft is 80% as good as antegrade flow

// ---- injury / event hazards ----------------------------------------------
K_INF     :   0.0165 : infarct accrual gain (/d)
K_INF_EMB :   4.0e-5 : infarct accrual per unit embolic hazard (/d)
ISCH_HAZ0 :   0.0480 : ischaemic-event hazard at zero margin (/yr)
CBF_SD0   :   4.50   : autoregulating perfusion variability (mL/100g/min)
K_PASS    :   0.70   : extra variability when the arteriole is on its floor
ISCH_EMB  :   0.075  : ischaemic hazard per unit embolic hazard (/yr)
K_EMB     :   0.35   : embolic hazard gain from low-shear perforators
K_EMB_ASA :   0.55   : max fractional embolic reduction by antiplatelet
K_TIA     :   1.0    : accumulator for time below threshold

// ---- perforator wall / haemorrhage --------------------------------------
K_ANEU    :   0.0110 : microaneurysm formation gain (/d)
K_ANEU_REP:   0.0060 : microaneurysm repair (/d)
K_MINO    :   0.45   : minocycline boost to wall repair (fraction)
HEM_HAZ0  :   0.4140 : haemorrhage hazard at ANEU=1, sigma=1 (/yr)
HEM_POW   :   2.0    : wall-stress exponent of the haemorrhage hazard
HEM_VOL   :  26.0    : mean haematoma volume per event (mL)

// ---- cognition / BBB ----------------------------------------------------
K_COG     :   0.00090: z-score loss per unit relative hypoperfusion (/d)
K_COG_INF :   3.10   : z-score loss per unit infarct fraction
K_BBB     :   0.45   : BBB breakdown gain (/d)
K_BBB_REP :   0.25   : BBB repair (/d)
HYPER_THR :   1.35   : CBF/CBF0 above which hyperperfusion injures
K_EDEMA   :   0.85   : oedema formation (/d)
K_EDEMA_R :   0.30   : oedema resolution (/d)

// ---- drug PK ------------------------------------------------------------
ASA_KA    :  12.0  : aspirin absorption (/d)
ASA_KE    :   8.3  : aspirin elimination (/d)
ASA_V     :  12.0  : aspirin volume (L)
K_PLT_ON  :  45.0  : platelet COX-1 inactivation gain (L/mg/d)
K_PLT_OFF :   0.10 : platelet pool turnover (/d)
CILO_KA   :   4.2  : cilostazol absorption (/d)
CILO_KE   :   1.51 : cilostazol elimination (/d)
CILO_V    :  95.0  : cilostazol volume (L)
CILO_EC50 :   0.55 : cilostazol EC50 (mg/L)
CILO_EMAX :   0.22 : cilostazol max arteriolar dilatory stimulus
CILO_APL  :   0.45 : cilostazol antiplatelet fraction
NIF_KA    :   0.55 : nifedipine GITS absorption (/d)
NIF_KE    :   2.31 : nifedipine elimination (/d)
NIF_V     :  55.0  : nifedipine volume (L)
NIF_EC50  :   0.030: nifedipine EC50 (mg/L)
NIF_EMAX_CBR: 0.30 : nifedipine max cerebral dilatory stimulus
NIF_EMAX_MAP: 0.16 : nifedipine max fractional MAP reduction
STAT_KA   :   6.0  : atorvastatin absorption (/d)
STAT_KE   :   0.50 : atorvastatin effect-site turnover (/d)
STAT_V    : 380.0  : atorvastatin volume (L)
STAT_EC50 :   0.012: atorvastatin EC50 (mg/L)
MIN_KA    :   5.0  : minocycline absorption (/d)
MIN_KE    :   0.77 : minocycline elimination (/d)
MIN_V     : 100.0  : minocycline volume (L)
MIN_EC50  :   1.20 : minocycline EC50 (mg/L)
AZ_KE     :   3.30 : acetazolamide elimination (/d)
AZ_V      :  18.0  : acetazolamide volume (L)
AZ_EC50   :   8.0  : acetazolamide EC50 (mg/L)
AZ_EMAX   :   0.45 : acetazolamide max arteriolar dilatory stimulus
AHT_KA    :   8.0  : antihypertensive absorption (/d)
AHT_KE    :   1.0  : antihypertensive elimination (/d)
AHT_V     :  60.0  : antihypertensive volume (L)
AHT_EC50  :   0.10 : antihypertensive EC50 (mg/L)
AHT_EMAX  :   0.22 : antihypertensive max fractional MAP reduction

// ---- CO2 and pressure forcings -----------------------------------------
PACO2_0   :  40.0  : reference arterial PaCO2 (mmHg)
PACO2     :  40.0  : current arterial PaCO2 (mmHg) — forcing input
K_CO2     :   0.0350: arteriolar DEMAND change per mmHg PaCO2
K_CO2_COLL:   0.0200: leptomeningeal CONDUCTANCE change per mmHg PaCO2
K_CO2_MOYA:   0.0110: perforator CONDUCTANCE change per mmHg PaCO2
MAP_MULT  :   1.0  : multiplicative MAP policy (e.g. post-op BP control)
MAP_MULT_T:  -1e9  : day the MAP policy starts (so it cannot rewrite history)
PROG      :   1.0  : lesion-progression switch (0 freezes the lesion)

$CMT @annotated
PDGF   : PDGF-BB signal (a.u.)
SMC    : intimal smooth-muscle mass, terminal ICA (0-1)
GMOYA  : basal moyamoya perforator conductance (mL/min/mmHg)
GPVA   : periventricular/choroidal conductance (mL/min/mmHg)
GCOLL  : leptomeningeal collateral conductance (mL/min/mmHg)
GBYP   : surgical bypass conductance (mL/min/mmHg)
REMOD  : arteriolar outward remodelling, loss of constrictive range (0-1)
HIF    : HIF-1a activity (a.u.)
VEGF   : VEGF-A activity (a.u.)
MMP9   : MMP-9 activity (a.u.)
INFA   : infarcted fraction of territory A (0-1)
INFB   : infarcted fraction of territory B (0-1)
ANEU   : perforator microaneurysm burden (0-1)
HEMV   : cumulative haemorrhage volume (mL)
COG    : cognitive z-score change (z)
BBB    : blood-brain barrier disruption (a.u.)
EDEMA  : vasogenic oedema (a.u.)
EMBH   : cumulative embolic exposure (a.u.-d)
TIAB   : cumulative days with CBF_ws below threshold (d)
HYPOX  : cumulative relative hypoperfusion exposure (a.u.-d)
HEMH   : cumulative haemorrhage hazard (expected events)
ISCH   : cumulative ischaemic-event hazard (expected events)
CBFAD  : flow the peri-anastomotic barrier has adapted to (mL/100g/min)
ASA_G  : aspirin gut (mg)
ASA_C  : aspirin central (mg)
PLTI   : platelet COX-1 inhibition (0-1)
CILO_G : cilostazol gut (mg)
CILO_C : cilostazol central (mg)
NIF_G  : nifedipine gut (mg)
NIF_C  : nifedipine central (mg)
STAT_G : atorvastatin gut (mg)
STAT_C : atorvastatin effect site (mg)
MIN_G  : minocycline gut (mg)
MIN_C  : minocycline central (mg)
AZ_C   : acetazolamide central (mg)
AHT_G  : antihypertensive gut (mg)
AHT_C  : antihypertensive central (mg)

$GLOBAL
#define _EPS 1e-12

namespace mmd {

  inline double clamp(double x, double lo, double hi) {
    return (x < lo) ? lo : ((x > hi) ? hi : x);
  }

  struct Net { double PA; double PF; double PB; };

  // Three-node linear solve.
  //
  //     Pa --[gSA]--+                 Pa --[gSF + gb]--+
  //                 v                                  v
  //     P_B --[gc]-> P_A <---[gl]---> P_F      (peri-anastomotic cortex)
  //                  |                 |
  //                [gA]              [gF]
  //                  v                 v
  //                  Pv                Pv
  //
  //   A: gSA(Pa-P_A) + gc(P_B-P_A) + gl(P_F-P_A) = gA(P_A-Pv)
  //   F: (gSF+gb)(Pa-P_F) + gl(P_A-P_F)          = gF(P_F-Pv)
  //   B: gp(Pa-P_B) + gc(P_A-P_B)                = gB(P_B-Pv)
  //
  // With gb = 0, and gSA:gSF and gA:gF both split by mass, P_A == P_F
  // exactly -- so the focal compartment is invisible until an operation.
  // Solved by Cramer determinants (3x3, always well conditioned here).
  inline Net solve3(double Pa, double Pv, double gSA, double gSF, double gb,
                    double gp, double gc, double gl,
                    double gA, double gF, double gB) {
    double gSFb = gSF + gb;
    double a11 = -(gSA + gc + gl + gA), a12 = gl,                  a13 = gc;
    double a21 = gl,                    a22 = -(gSFb + gl + gF),   a23 = 0.0;
    double a31 = gc,                    a32 = 0.0,                 a33 = -(gp + gc + gB);
    double b1 = -(gSA * Pa + gA * Pv);
    double b2 = -(gSFb * Pa + gF * Pv);
    double b3 = -(gp * Pa + gB * Pv);
    double det = a11 * (a22 * a33 - a23 * a32)
               - a12 * (a21 * a33 - a23 * a31)
               + a13 * (a21 * a32 - a22 * a31);
    double d1  = b1  * (a22 * a33 - a23 * a32)
               - a12 * (b2  * a33 - a23 * b3)
               + a13 * (b2  * a32 - a22 * b3);
    double d2  = a11 * (b2  * a33 - a23 * b3)
               - b1  * (a21 * a33 - a23 * a31)
               + a13 * (a21 * b3  - b2  * a31);
    double d3  = a11 * (a22 * b3  - b2  * a32)
               - a12 * (a21 * b3  - b2  * a31)
               + b1  * (a21 * a32 - a22 * a31);
    Net n;
    n.PA = d1 / det; n.PF = d2 / det; n.PB = d3 / det;
    return n;
  }

  struct Auto { double gA, gF, gB, PA, PF, PB, QA, QF, QB; };

  // Quasi-steady autoregulation on the three-node network.  Each compartment
  // picks the arteriolar conductance that would deliver its demanded flow,
  // CLAMPED to its own dilatory floor and constrictive ceiling.  Fixed-point
  // iterated because the compartments share conduits.  This clamp is the ONE
  // place the central assumption of the model lives.
  inline Auto autoreg(double Pa, double Pv, double gSA, double gSF, double gb,
                      double gp, double gc, double gl,
                      double QtA, double QtF, double QtB,
                      double gAc, double gAd, double gFc, double gFd,
                      double gBc, double gBd,
                      double gA0, double gF0, double gB0, int iters) {
    double gA = gA0, gF = gF0, gB = gB0;
    Net n = solve3(Pa, Pv, gSA, gSF, gb, gp, gc, gl, gA, gF, gB);
    for (int i = 0; i < iters; ++i) {
      double PAr = (gSA * Pa + gc * n.PB + gl * n.PF - QtA) / (gSA + gc + gl);
      gA = (PAr - Pv > 1e-9) ? QtA / (PAr - Pv) : gAd;
      gA = clamp(gA, gAc, gAd);
      double PFr = ((gSF + gb) * Pa + gl * n.PA - QtF) / (gSF + gb + gl);
      gF = (PFr - Pv > 1e-9) ? QtF / (PFr - Pv) : gFd;
      gF = clamp(gF, gFc, gFd);
      double PBr = (gp * Pa + gc * n.PA - QtB) / (gp + gc);
      gB = (PBr - Pv > 1e-9) ? QtB / (PBr - Pv) : gBd;
      gB = clamp(gB, gBc, gBd);
      n = solve3(Pa, Pv, gSA, gSF, gb, gp, gc, gl, gA, gF, gB);
    }
    Auto a;
    a.gA = gA; a.gF = gF; a.gB = gB;
    a.PA = n.PA; a.PF = n.PF; a.PB = n.PB;
    a.QA = gA * (n.PA - Pv);
    a.QF = gF * (n.PF - Pv);
    a.QB = gB * (n.PB - Pv);
    return a;
  }

  // ---- parameter block, filled once per record in $MAIN ----------------
  struct Par {
    double MAP0, MAP_REF, ICP, MASS_A, MASS_B, HB, SAO2, CMRO2_0, OEF_MAX,
           OEF_BASE, CBF_A0, CBF_B0, FRAC_PROX, DIL_MAX, CON_MAX,
           GCOLL0, GMOYA0, GPVA0, STEN_MAX, REM_CON, WS_K, W_BYP,
           GMOYA_RREF, GPVA_RREF, Q_LOAD_REF, PCA_STEN, PCA_INV,
           F_FRAC, G_LEAK,
           PACO2, PACO2_0, K_CO2, K_CO2_COLL, K_CO2_MOYA, MAP_MULT,
           MAP_MULT_ON,
           CILO_V, CILO_EC50, CILO_EMAX, NIF_V, NIF_EC50, NIF_EMAX_CBR,
           NIF_EMAX_MAP, AZ_V, AZ_EC50, AZ_EMAX, AHT_V, AHT_EC50, AHT_EMAX;
    // derived
    double CAO2, CBF_CRIT, QA0, QB0, GA0, GB0, GA_DIL, GA_CON,
           GB_DIL, GB_CON, GI0, GP0, SIG_REF;
    // focal (peri-anastomotic) split, all derived from F_FRAC
    double MASS_F, MASS_At, QF0, QAt0, GF0, GAt0,
           GF_DIL, GF_CON, GAt_DIL, GAt_CON;
  };

  // ---- the complete algebraic haemodynamic layer -----------------------
  // Called IDENTICALLY from $ODE and $TABLE, so the reported values are the
  // ones the derivatives were built from -- no stale-global hazard.
  struct Hemo {
    double STEN, gi, gm, gv, gb, gc, gS, gp, gl, gSA, gSF, Pa, dem;
    double PA, PF, PB, gA, gF, gB, gAc, gAd, gFc, gFd, QA, QF, QB, ar_pos;
    double CBFA, CBFF, CBFT, CBFB, CBF_WS;
    double OEFA, OEFB, OEF_WS, CMROA, CMROB;
    double ante, k_ws, sigma, sig_pva, load, load_pva, P_perf;
    double q_moya, q_pva, q_coll, q_ica, q_byp, q_leak;
    double e_cilo, e_nif, e_az, e_aht, on_floor;
  };

  inline Hemo hemo(const Par& p,
                   double SMC, double GMOYA, double GPVA, double GCOLL,
                   double GBYP, double REMOD, double INFA, double INFB,
                   double CILO_C, double NIF_C, double AZ_C, double AHT_C) {
    Hemo h;
    h.STEN = p.STEN_MAX * clamp(SMC, 0.0, 1.0);

    // CO2 acts on the COLLATERAL CONDUITS as well as on the arteriole.
    double dco2   = p.PACO2 - p.PACO2_0;
    double f_coll = clamp(1.0 + p.K_CO2_COLL * dco2, 0.25, 1.80);
    double f_moya = clamp(1.0 + p.K_CO2_MOYA * dco2, 0.35, 1.60);

    // the four inlets, in parallel
    h.gi = p.GI0 * pow(1.0 - h.STEN, 4.0);            // antegrade, Poiseuille
    h.gm = (GMOYA > 0.0 ? GMOYA : 0.0) * f_moya;      // basal perforators
    h.gv = (GPVA  > 0.0 ? GPVA  : 0.0) * f_moya;      // periventricular route
    h.gb = (GBYP  > 0.0 ? GBYP  : 0.0);               // graft: not CO2-gated
    h.gc = (GCOLL > 0.0 ? GCOLL : 0.0) * f_coll;      // leptomeningeal (PCA)
    h.gS = h.gi + h.gm + h.gv + h.gb;
    h.gp = p.GP0;
    h.gl = p.G_LEAK;
    // native inlets are shared between the territorial and peri-anastomotic
    // cortex in proportion to mass, so the focal compartment is invisible
    // until a graft is added to it
    double gnat = h.gi + h.gm + h.gv;
    h.gSA = gnat * (1.0 - p.F_FRAC);
    h.gSF = gnat * p.F_FRAC;

    // drug effects
    double c_cilo = CILO_C / p.CILO_V, c_nif = NIF_C / p.NIF_V;
    double c_az   = AZ_C   / p.AZ_V,   c_aht = AHT_C / p.AHT_V;
    h.e_cilo = p.CILO_EMAX * c_cilo / (p.CILO_EC50 + c_cilo);
    h.e_nif  = p.NIF_EMAX_CBR * c_nif / (p.NIF_EC50 + c_nif);
    double e_nif_map = p.NIF_EMAX_MAP * c_nif / (p.NIF_EC50 + c_nif);
    h.e_az  = p.AZ_EMAX  * c_az  / (p.AZ_EC50  + c_az);
    h.e_aht = p.AHT_EMAX * c_aht / (p.AHT_EC50 + c_aht);

    // a MAP policy applies only from the day it is started -- otherwise a
    // post-operative BP target would silently rewrite the whole pre-operative
    // history and make the hemisphere a different one
    double mm = (p.MAP_MULT_ON > 0.5) ? p.MAP_MULT : 1.0;
    h.Pa = p.MAP0 * (1.0 - e_nif_map) * (1.0 - h.e_aht) * mm;

    // arteriolar demand
    double fco2 = 1.0 + p.K_CO2 * dco2;
    if (fco2 < 0.35) fco2 = 0.35;
    h.dem = fco2 * (1.0 + h.e_az) * (1.0 + h.e_cilo) * (1.0 + h.e_nif);
    double QtA = p.QAt0 * h.dem * (1.0 - INFA);  // infarct removes DEMAND
    double QtF = p.QF0  * h.dem;
    double QtB = p.QB0  * h.dem * (1.0 - INFB);

    // the clamps; the CONSTRICTIVE one is widened by chronic remodelling
    // (vasoparalysis), which is what makes hyperperfusion possible at all
    double rem = clamp(REMOD, 0.0, 1.0);
    h.gAd = p.GAt_DIL;
    h.gFd = p.GF_DIL;
    h.gAc = p.GAt_CON * (1.0 + (p.REM_CON - 1.0) * rem);
    h.gFc = p.GF_CON  * (1.0 + (p.REM_CON - 1.0) * rem);
    if (h.gAc > h.gAd * 0.999) h.gAc = h.gAd * 0.999;
    if (h.gFc > h.gFd * 0.999) h.gFc = h.gFd * 0.999;

    Auto A = autoreg(h.Pa, p.ICP, h.gSA, h.gSF, h.gb, h.gp, h.gc, h.gl,
                     QtA, QtF, QtB, h.gAc, h.gAd, h.gFc, h.gFd,
                     p.GB_CON, p.GB_DIL, p.GAt0, p.GF0, p.GB0, 30);
    h.PA = A.PA; h.PF = A.PF; h.PB = A.PB;
    h.gA = A.gA; h.gF = A.gF; h.gB = A.gB;
    h.QA = A.QA; h.QF = A.QF; h.QB = A.QB;
    h.ar_pos = h.gA / h.gAd;              // 1.0 == on the dilatory floor
    h.on_floor = 1.0 / (1.0 + exp(-(h.ar_pos - 0.93) / 0.02));

    double massA = p.MASS_At * (1.0 - INFA) + 1e-6;
    double massF = p.MASS_F;
    double massB = p.MASS_B  * (1.0 - INFB) + 1e-6;
    h.CBFA = h.QA / massA * 100.0;     // per 100 g of SURVIVING tissue
    h.CBFF = h.QF / massF * 100.0;     // peri-anastomotic cortex
    h.CBFB = h.QB / massB * 100.0;
    // what a SPECT of the WHOLE MCA territory would report
    h.CBFT = (h.QA + h.QF) / (massA + massF) * 100.0;

    // oxygen: OEF rises to a ceiling, then CMRO2 falls
    double dA = h.CBFA * p.CAO2;
    h.OEFA = clamp(p.CMRO2_0 / (dA > _EPS ? dA : _EPS), p.OEF_BASE, p.OEF_MAX);
    h.CMROA = dA * h.OEFA;
    double dB = h.CBFB * p.CAO2;
    h.OEFB = clamp(p.CMRO2_0 / (dB > _EPS ? dB : _EPS), p.OEF_BASE, p.OEF_MAX);
    h.CMROB = dB * h.OEFB;

    // watershed: collateral flow reaches the CORTEX first, so the deep
    // borderzone is supplied last.  The penalty scales with how
    // COLLATERAL-dependent the supply is, not with how low the mean flow is.
    // Antegrade FLOW fraction, not conductance fraction: the graft reaches the
    // territorial cortex through the pial coupling, so its contribution is
    // whatever actually crosses that connection.
    h.q_leak = h.gl * (h.PF - h.PA);
    double q_in_A = h.gSA * (h.Pa - h.PA) + h.gc * (h.PB - h.PA) + h.q_leak;
    if (q_in_A < 1e-9) q_in_A = 1e-9;
    double q_ante = h.gi * (1.0 - p.F_FRAC) * (h.Pa - h.PA)
                    + p.W_BYP * (h.q_leak > 0.0 ? h.q_leak : 0.0);
    h.ante   = clamp(q_ante / q_in_A, 0.0, 1.0);
    h.k_ws   = 1.0 - p.WS_K * (1.0 - h.ante);
    h.CBF_WS = h.CBFA * h.k_ws;
    double dW = h.CBF_WS * p.CAO2;
    h.OEF_WS = clamp(p.CMRO2_0 / (dW > _EPS ? dW : _EPS),
                     p.OEF_BASE, p.OEF_MAX);

    // perforator wall stress.  G ~ N*r^4 with both N and r growing, so
    // r ~ (G/Gref)^(1/8); Laplace with a wall that thins as it dilates
    // (h ~ r^-1/2) gives sigma ~ P_transmural * dilatation^1.5.
    h.P_perf = 0.5 * (h.Pa + h.PA);
    double d1 = h.gm / p.GMOYA_RREF; if (d1 < 1.0) d1 = 1.0;
    double d2 = h.gv / p.GPVA_RREF;  if (d2 < 1.0) d2 = 1.0;
    h.sigma   = h.P_perf * pow(pow(d1, 0.125), 1.5) / p.SIG_REF;
    h.sig_pva = h.P_perf * pow(pow(d2, 0.125), 1.5) / p.SIG_REF;
    h.q_moya = h.gm * (h.Pa - h.PA);
    h.q_pva  = h.gv * (h.Pa - h.PA);
    h.q_coll = h.gc * (h.PB - h.PA);
    h.q_ica  = h.gi * (h.Pa - h.PA);
    h.q_byp  = h.gb * (h.Pa - h.PF);      // the graft feeds the FOCAL node
    h.load     = h.q_moya / p.Q_LOAD_REF;
    h.load_pva = h.q_pva  / p.Q_LOAD_REF;
    return h;
  }
}

mmd::Par PB_;   // parameter block, refreshed in $MAIN

$PREAMBLE
// nothing: every derived constant depends on a parameter and so is set in
// $MAIN, which runs once per record and therefore once per ID.

$MAIN
PB_.MAP0 = MAP0;   PB_.MAP_REF = MAP_REF; PB_.ICP = ICP;
PB_.MASS_A = MASS_A; PB_.MASS_B = MASS_B; PB_.HB = HB; PB_.SAO2 = SAO2;
PB_.CMRO2_0 = CMRO2_0; PB_.OEF_MAX = OEF_MAX; PB_.OEF_BASE = OEF_BASE;
PB_.CBF_A0 = CBF_A0; PB_.CBF_B0 = CBF_B0; PB_.FRAC_PROX = FRAC_PROX;
PB_.DIL_MAX = DIL_MAX; PB_.CON_MAX = CON_MAX;
PB_.GCOLL0 = GCOLL0; PB_.GMOYA0 = GMOYA0; PB_.GPVA0 = GPVA0;
PB_.STEN_MAX = STEN_MAX; PB_.REM_CON = REM_CON;
PB_.WS_K = WS_K; PB_.W_BYP = W_BYP;
PB_.GMOYA_RREF = GMOYA_RREF; PB_.GPVA_RREF = GPVA_RREF;
PB_.Q_LOAD_REF = Q_LOAD_REF; PB_.PCA_STEN = PCA_STEN; PB_.PCA_INV = PCA_INV;
PB_.F_FRAC = F_FRAC; PB_.G_LEAK = G_LEAK;
PB_.PACO2 = PACO2; PB_.PACO2_0 = PACO2_0; PB_.K_CO2 = K_CO2;
PB_.K_CO2_COLL = K_CO2_COLL; PB_.K_CO2_MOYA = K_CO2_MOYA;
PB_.MAP_MULT = MAP_MULT;
PB_.MAP_MULT_ON = (TIME >= MAP_MULT_T) ? 1.0 : 0.0;
PB_.CILO_V = CILO_V; PB_.CILO_EC50 = CILO_EC50; PB_.CILO_EMAX = CILO_EMAX;
PB_.NIF_V = NIF_V; PB_.NIF_EC50 = NIF_EC50;
PB_.NIF_EMAX_CBR = NIF_EMAX_CBR; PB_.NIF_EMAX_MAP = NIF_EMAX_MAP;
PB_.AZ_V = AZ_V; PB_.AZ_EC50 = AZ_EC50; PB_.AZ_EMAX = AZ_EMAX;
PB_.AHT_V = AHT_V; PB_.AHT_EC50 = AHT_EC50; PB_.AHT_EMAX = AHT_EMAX;

// ---- derived oxygen constants: THE THRESHOLD IS DERIVED -----------------
PB_.CAO2     = 1.34 * HB * SAO2 / 100.0;             // mL O2 / mL blood
PB_.CBF_CRIT = CMRO2_0 / (PB_.CAO2 * OEF_MAX);       // mL/100g/min

// ---- derived resistive geometry ----------------------------------------
double cpp   = MAP0 - ICP;
PB_.QA0      = CBF_A0 * MASS_A / 100.0;              // mL/min
PB_.QB0      = CBF_B0 * MASS_B / 100.0;
double R_ICA0  = FRAC_PROX * (cpp / PB_.QA0);
double R_ARTA0 = (1.0 - FRAC_PROX) * (cpp / PB_.QA0);
double R_PCA0  = FRAC_PROX * (cpp / PB_.QB0);
double R_ARTB0 = (1.0 - FRAC_PROX) * (cpp / PB_.QB0);
PB_.GA0    = 1.0 / R_ARTA0;
PB_.GB0    = 1.0 / R_ARTB0;
PB_.GA_DIL = DIL_MAX / R_ARTA0;
PB_.GA_CON = 1.0 / (R_ARTA0 * CON_MAX);
PB_.GB_DIL = DIL_MAX / R_ARTB0;
PB_.GB_CON = 1.0 / (R_ARTB0 * CON_MAX);
PB_.GI0    = 1.0 / R_ICA0;
// focal (peri-anastomotic) split of masses, demands and arteriolar limits
PB_.MASS_F  = MASS_A * F_FRAC;
PB_.MASS_At = MASS_A * (1.0 - F_FRAC);
PB_.QF0     = PB_.QA0 * F_FRAC;
PB_.QAt0    = PB_.QA0 * (1.0 - F_FRAC);
PB_.GF0     = PB_.GA0 * F_FRAC;
PB_.GAt0    = PB_.GA0 * (1.0 - F_FRAC);
PB_.GF_DIL  = PB_.GA_DIL * F_FRAC;
PB_.GF_CON  = PB_.GA_CON * F_FRAC;
PB_.GAt_DIL = PB_.GA_DIL * (1.0 - F_FRAC);
PB_.GAt_CON = PB_.GA_CON * (1.0 - F_FRAC);
PB_.GP0    = (1.0 / R_PCA0) * pow(1.0 - PCA_STEN * PCA_INV, 4.0);

// MAP-INDEPENDENT wall-stress normaliser, at the healthy operating point
double cpp_ref = MAP_REF - ICP;
PB_.SIG_REF = 0.5 * (MAP_REF + (ICP + cpp_ref * (1.0 - FRAC_PROX)));

$ODE
mmd::Hemo H = mmd::hemo(PB_, SMC, GMOYA, GPVA, GCOLL, GBYP, REMOD,
                        INFA, INFB, CILO_C, NIF_C, AZ_C, AHT_C);

// ---- drug PK -----------------------------------------------------------
dxdt_ASA_G  = -ASA_KA * ASA_G;
dxdt_ASA_C  =  ASA_KA * ASA_G - ASA_KE * ASA_C;
double c_asa = ASA_C / ASA_V;
dxdt_PLTI   = K_PLT_ON * c_asa * (1.0 - PLTI) - K_PLT_OFF * PLTI;
dxdt_CILO_G = -CILO_KA * CILO_G;
dxdt_CILO_C =  CILO_KA * CILO_G - CILO_KE * CILO_C;
dxdt_NIF_G  = -NIF_KA  * NIF_G;
dxdt_NIF_C  =  NIF_KA  * NIF_G  - NIF_KE  * NIF_C;
dxdt_STAT_G = -STAT_KA * STAT_G;
dxdt_STAT_C =  STAT_KA * STAT_G - STAT_KE * STAT_C;
dxdt_MIN_G  = -MIN_KA  * MIN_G;
dxdt_MIN_C  =  MIN_KA  * MIN_G  - MIN_KE  * MIN_C;
dxdt_AHT_G  = -AHT_KA  * AHT_G;
dxdt_AHT_C  =  AHT_KA  * AHT_G  - AHT_KE  * AHT_C;
dxdt_AZ_C   = -AZ_KE   * AZ_C;

double c_stat = STAT_C / STAT_V;
double c_min  = MIN_C / MIN_V;
double e_stat = c_stat / (STAT_EC50 + c_stat);
double e_min  = c_min  / (MIN_EC50  + c_min);

// ---- lesion: (RNF213 + secondary drive) -> PDGF -> intimal SMC ---------
dxdt_PDGF = K_PDGF * (RNF + SEC) * (1.0 + 0.6 * MMP9) - K_PDGF_D * PDGF;
dxdt_SMC  = K_SMC * (1.0 - K_STATIN * e_stat) * PDGF
            * (1.0 - SMC / SMC_MAX) * PROG;

// ---- hypoxic drive -----------------------------------------------------
double hypo = 1.0 - H.CBFA / CBF_A0; if (hypo < 0.0) hypo = 0.0;
dxdt_HIF  = K_HIF_ON * hypo - K_HIF_OFF * HIF;
dxdt_VEGF = K_VEGF   * HIF  - K_VEGF_D  * VEGF;
dxdt_MMP9 = K_MMP    * VEGF - K_MMP_D   * MMP9;

// ---- collaterals: RNF213 caps the ceiling it also creates --------------
double cap_m = GMOYA_CAP * (1.0 - RNF_CAP * RNF) * ANGIO;
double cap_c = GCOLL_CAP * (1.0 - RNF_CAP * RNF) * ANGIO
               * (1.0 - PCA_COLL * PCA_INV);
double cap_v = GPVA_CAP  * (1.0 - RNF_CAP * RNF) * ANGIO
               * (1.0 + PCA_PVA  * PCA_INV);
if (cap_m < GMOYA0) cap_m = GMOYA0;
if (cap_c < GCOLL0) cap_c = GCOLL0;
if (cap_v < GPVA0 ) cap_v = GPVA0;
double reg = 1.0 - VEGF / VEGF_REG; if (reg < 0.0) reg = 0.0;
double hm = 1.0 - GMOYA / cap_m; if (hm < 0.0) hm = 0.0;
double hc = 1.0 - GCOLL / cap_c; if (hc < 0.0) hc = 0.0;
double hv = 1.0 - GPVA  / cap_v; if (hv < 0.0) hv = 0.0;
double xm = GMOYA - GMOYA0; if (xm < 0.0) xm = 0.0;
double xc = GCOLL - GCOLL0; if (xc < 0.0) xc = 0.0;
double xv = GPVA  - GPVA0 ; if (xv < 0.0) xv = 0.0;
dxdt_GMOYA = K_MOYA * VEGF * ANGIO * hm - K_MOYA_REG * reg * xm;
dxdt_GCOLL = K_COLL * VEGF * ANGIO * hc - K_COLL_REG * reg * xc;
// the periventricular route is recruited by the DEEP borderzone deficit
// specifically -- which is what makes it a marker of watershed failure
double dvw = 1.0 - H.CBF_WS / CBF_A0; if (dvw < 0.0) dvw = 0.0;
dxdt_GPVA  = K_PVA * VEGF * dvw * ANGIO * hv - K_PVA_REG * reg * xv;

// ---- surgery: the only intervention that changes gS -------------------
dxdt_GBYP = 0.0;
if (SURG_KIND > 0.5 && SOLVERTIME >= SURG_T) {
  if (SURG_KIND < 1.5) {                                   // direct
    dxdt_GBYP = (GBYP_DIR * SURG_SCALE - GBYP) / TAU_BYP;
  } else if (SURG_KIND < 2.5) {                            // indirect
    double gap = GBYP_IND * SURG_SCALE - GBYP; if (gap < 0.0) gap = 0.0;
    dxdt_GBYP = K_BYP_IND * VEGF * ANGIO * gap;
  } else {                                                 // combined
    double gd = GBYP_DIR * SURG_SCALE - GBYP; if (gd < 0.0) gd = 0.0;
    double gt = (GBYP_DIR + GBYP_IND) * SURG_SCALE - GBYP;
    if (gt < 0.0) gt = 0.0;
    dxdt_GBYP = gd / TAU_BYP + K_BYP_IND * VEGF * ANGIO * gt;
  }
}

// ---- arteriolar remodelling: the price of chronic maximal dilatation --
dxdt_REMOD = K_REM_ON * H.on_floor * (1.0 - REMOD) - K_REM_OFF * REMOD;

// ---- ischaemic injury -------------------------------------------------
double CBFCRIT_ = PB_.CBF_CRIT;
double sevA = (CBFCRIT_ - H.CBF_WS) / CBFCRIT_; if (sevA < 0.0) sevA = 0.0;
double sevB = (CBFCRIT_ - H.CBFB)   / CBFCRIT_; if (sevB < 0.0) sevB = 0.0;
double apl  = PLTI + CILO_APL * (H.e_cilo / (CILO_EMAX > _EPS ? CILO_EMAX : _EPS));
if (apl > 1.0) apl = 1.0;
double emb  = K_EMB * GMOYA * (1.0 - K_EMB_ASA * apl); if (emb < 0.0) emb = 0.0;
double roomA = WS_FRAC - INFA; if (roomA < 0.0) roomA = 0.0;

dxdt_INFA  = K_INF * roomA * sevA * sevA + K_INF_EMB * emb * (1.0 - INFA);
dxdt_INFB  = K_INF * (1.0 - INFB) * sevB * sevB;
dxdt_EMBH  = emb;
dxdt_HYPOX = hypo;
dxdt_TIAB  = (H.CBF_WS < CBFCRIT_) ? K_TIA : 0.0;

// Clinical ischaemic events are threshold CROSSINGS by a fluctuating
// perfusion, not a state of permanent sub-threshold flow.  The natural form
// is therefore exponential in the MARGIN to CBF_crit, with a scale equal to
// the perfusion variability -- which is itself WIDER when the arteriole is on
// its floor, because there MAP noise passes straight through to the tissue.
double sd     = CBF_SD0 * (1.0 + K_PASS * H.on_floor);
double margin = mmd::clamp((H.CBF_WS - CBFCRIT_) / sd, -6.0, 25.0);
dxdt_ISCH = (ISCH_HAZ0 * exp(-margin) + ISCH_EMB * emb) / 365.0;

// ---- perforator microaneurysm and haemorrhage -------------------------
double lpv = H.load_pva > 0.0 ? H.load_pva : 0.0;
double stress = H.sig_pva * lpv;
dxdt_ANEU = K_ANEU * stress * stress * (1.0 - ANEU)
            - K_ANEU_REP * (1.0 + K_MINO * e_min) * ANEU;
double haz = HEM_HAZ0 * ANEU * pow(H.sig_pva, HEM_POW);      // per year
dxdt_HEMH = haz / 365.0;
dxdt_HEMV = haz / 365.0 * HEM_VOL;

// ---- hyperperfusion -> BBB -> oedema ---------------------------------
// FOCAL and RELATIVE: read on the peri-anastomotic cortex, against the flow
// that cortex had adapted to.  50 mL/100g/min is normal, and it injures a
// barrier that has lived on 18 for five years.  REMOD sets the HEIGHT of the
// surge; TAU_ADAPT sets its DURATION.
dxdt_CBFAD = (H.CBFF - CBFAD) / TAU_ADAPT;
double hyper_rel = H.CBFF / (CBFAD > 5.0 ? CBFAD : 5.0);
double over = hyper_rel - HYPER_THR; if (over < 0.0) over = 0.0;
dxdt_BBB   = K_BBB   * over - K_BBB_REP * BBB;
dxdt_EDEMA = K_EDEMA * BBB  - K_EDEMA_R * EDEMA;

// ---- cognition -------------------------------------------------------
dxdt_COG = -K_COG * hypo - K_COG_INF * dxdt_INFA;

$TABLE
// Recomputed with the SAME function $ODE used, so every reported number is
// consistent with the derivatives at this time point.
mmd::Hemo H = mmd::hemo(PB_, SMC, GMOYA, GPVA, GCOLL, GBYP, REMOD,
                        INFA, INFB, CILO_C, NIF_C, AZ_C, AHT_C);

double CBFA    = H.CBFA;
double CBFF    = H.CBFF;      // peri-anastomotic cortex
double CBFT    = H.CBFT;      // what a territorial SPECT would report
double CBFB    = H.CBFB;
double CBFWS   = H.CBF_WS;
double OEFA    = H.OEFA;
double OEFWS   = H.OEF_WS;
double CMROA   = H.CMROA;
double CBFCRIT = PB_.CBF_CRIT;
double STEN    = H.STEN;
double PA      = H.PA;
double PF      = H.PF;
double PB      = H.PB;
double MAPx    = H.Pa;
double gS      = H.gS;
double gICA    = H.gi;
double gMOYA   = H.gm;
double gPVA    = H.gv;
double gCOLL   = H.gc;
double gBYP    = H.gb;
double gART    = H.gA;
double ARPOS   = H.ar_pos;        // 1.0 == arteriole on its dilatory floor
double ANTE    = H.ante;
double SIGPVA  = H.sig_pva;
double QMOYA   = H.q_moya;
double QPVA    = H.q_pva;
double QCOLL   = H.q_coll;
double QICA    = H.q_ica;
double QBYP    = H.q_byp;
double MARGIN  = H.CBF_WS - PB_.CBF_CRIT;
double PISCH   = 1.0 - exp(-ISCH);          // cumulative event probability
double PHEM    = 1.0 - exp(-HEMH);
double HYPER   = H.CBFF / CBF_A0;           // absolute, vs textbook normal
double HYPERREL = H.CBFF / (CBFAD > 5.0 ? CBFAD : 5.0);  // vs adapted flow
double QLEAK   = H.q_leak;
double INFPCT  = INFA * 100.0;

// ---- INTRINSIC reserve: dilate A maximally, hold B and the conduits ---
mmd::Net rn = mmd::solve3(H.Pa, ICP, H.gSA, H.gSF, H.gb, H.gp, H.gc, H.gl,
                          H.gAd, H.gF, H.gB);
double massA2 = PB_.MASS_At * (1.0 - INFA) + 1e-6;
double CBFAMAX = H.gAd * (rn.PA - ICP) / massA2 * 100.0;
double CVRINT  = (CBFAMAX - H.CBFA) / (H.CBFA > _EPS ? H.CBFA : _EPS) * 100.0;

// ---- MEASURED reactivity: an acetazolamide-like stimulus dilates BOTH -
mmd::Auto A2 = mmd::autoreg(H.Pa, ICP, H.gSA, H.gSF, H.gb, H.gp, H.gc, H.gl,
                            PB_.QAt0 * H.dem * (1.0 - INFA) * 1.45,
                            PB_.QF0  * H.dem * 1.45,
                            PB_.QB0  * H.dem * (1.0 - INFB) * 1.45,
                            H.gAc, H.gAd, H.gFc, H.gFd,
                            PB_.GB_CON, PB_.GB_DIL,
                            PB_.GAt0, PB_.GF0, PB_.GB0, 30);
double CVRMEAS = (A2.QA / massA2 * 100.0 - H.CBFA)
                 / (H.CBFA > _EPS ? H.CBFA : _EPS) * 100.0;

// ---- pressure-passivity dCBF/dMAP ------------------------------------
mmd::Auto A3 = mmd::autoreg(H.Pa + 5.0, ICP, H.gSA, H.gSF, H.gb, H.gp,
                            H.gc, H.gl,
                            PB_.QAt0 * H.dem * (1.0 - INFA),
                            PB_.QF0  * H.dem,
                            PB_.QB0  * H.dem * (1.0 - INFB),
                            H.gAc, H.gAd, H.gFc, H.gFd,
                            PB_.GB_CON, PB_.GB_DIL,
                            PB_.GAt0, PB_.GF0, PB_.GB0, 30);
double DCBFDMAP = (A3.QA / massA2 * 100.0 - H.CBFA) / 5.0;

$CAPTURE
CBFA CBFF CBFT CBFB CBFWS OEFA OEFWS CMROA CBFCRIT STEN PA PF PB MAPx
gS gICA gMOYA gPVA gCOLL gBYP gART ARPOS ANTE SIGPVA
QMOYA QPVA QCOLL QICA QBYP QLEAK MARGIN PISCH PHEM HYPER HYPERREL INFPCT
CVRINT CVRMEAS DCBFDMAP
'

mod <- mcode_cache("mmd_qsp", code)

## ===========================================================================
##  ARCHETYPES
##  Only ANGIO (angiogenic capacity), K_SMC (lesion growth rate), PCA_INV
##  (posterior-circulation involvement), Hb and the RNF213/secondary drive
##  differ between these patients.  Everything else in the results below is a
##  CONSEQUENCE of the network, not a separate assumption.
## ===========================================================================
ARCHETYPE <- list(
  paediatric   = list(RNF = 1.0, SEC = 0.0, ANGIO = 1.00, K_SMC = 5.5e-4,
                      PCA_INV = 0.0, HB = 15.0),
  adult_isch   = list(RNF = 1.0, SEC = 0.0, ANGIO = 0.58, K_SMC = 2.6e-4,
                      PCA_INV = 0.0, HB = 15.0),
  adult_haem   = list(RNF = 1.0, SEC = 0.0, ANGIO = 0.86, K_SMC = 2.3e-4,
                      PCA_INV = 1.0, HB = 15.0),
  sickle_cell  = list(RNF = 0.0, SEC = 1.3, ANGIO = 0.85, K_SMC = 3.4e-4,
                      PCA_INV = 0.0, HB =  8.0),
  asymptomatic = list(RNF = 1.0, SEC = 0.0, ANGIO = 0.66, K_SMC = 1.05e-4,
                      PCA_INV = 0.0, HB = 15.0)
)

`%||%` <- function(a, b) if (is.null(a)) b else a

init_state <- function(p) {
  c(PDGF = unname(p$K_PDGF %||% 0.90) * (p$RNF + p$SEC) / 0.60,
    SMC = 0, GMOYA = 0.004, GPVA = 0.002, GCOLL = 0.050, GBYP = 0,
    REMOD = 0, HIF = 0, VEGF = 0, MMP9 = 0, INFA = 0, INFB = 0,
    ANEU = 0.02, HEMV = 0, COG = 0, BBB = 0, EDEMA = 0, EMBH = 0,
    TIAB = 0, HYPOX = 0, HEMH = 0, ISCH = 0, CBFAD = 50,
    ASA_G = 0, ASA_C = 0, PLTI = 0, CILO_G = 0, CILO_C = 0,
    NIF_G = 0, NIF_C = 0, STAT_G = 0, STAT_C = 0, MIN_G = 0, MIN_C = 0,
    AZ_C = 0, AHT_G = 0, AHT_C = 0)
}

## Convenience: build a mrgsolve run for one archetype + one intervention
run_mmd <- function(arch = "adult_isch", days = 3650, delta = 1,
                    drugs = NULL, surgery = NULL, pars = NULL,
                    map_mult = 1.0) {
  p <- ARCHETYPE[[arch]]
  m <- mod %>% param(p) %>% param(MAP_MULT = map_mult) %>% init(init_state(p))
  if (map_mult != 1.0 && !is.null(surgery))
    m <- m %>% param(MAP_MULT_T = surgery$day)   # post-op policy, not history
  if (!is.null(pars))   m <- m %>% param(pars)
  if (!is.null(surgery)) {
    kind <- match(surgery$kind, c("direct", "indirect", "combined"))
    m <- m %>% param(SURG_KIND = kind, SURG_T = surgery$day,
                     SURG_SCALE = surgery$scale %||% 1.0)
  }
  ## dosing
  ev <- NULL
  if (!is.null(drugs)) for (nm in names(drugs)) {
    d <- drugs[[nm]]
    cmt <- if (nm == "AZ") "AZ_C" else paste0(nm, "_G")
    e <- ev(amt = d$amt, cmt = cmt, time = d$start %||% 0,
            ii = d$ii %||% 0, addl = d$addl %||% 0)
    ev <- if (is.null(ev)) e else c(ev, e)
  }
  ## NOTE: PaCO2 is a PARAMETER, not an event.  Challenge scenarios set it
  ## through `pars = list(PACO2 = ...)` together with PROG = 0 to freeze the
  ## lesion, which is the right way to read an instantaneous probe.
  out <- if (is.null(ev)) m %>% mrgsim(end = days, delta = delta)
         else             m %>% mrgsim(events = ev, end = days, delta = delta)
  as_tibble(out)
}

## ===========================================================================
##  SCENARIO LIBRARY  (18 scenarios)
##  Each entry states WHAT IT IS FOR, i.e. which of the model's claims it
##  tests, so the file reads as an experiment list rather than a dose table.
## ===========================================================================
SCENARIOS <- list(

  ## --- natural history: the silent phase and the reserve cliff ----------
  S01_natural_paediatric = list(
    what = "Untreated child.  Fast lesion, high angiogenic capacity.  Shows
            reserve lost years before flow moves.",
    arch = "paediatric", days = 3650),

  S02_natural_adult_isch = list(
    what = "Untreated adult, ischaemic phenotype.  The reference natural
            history: CVRINT reaches 0 around y2, CBFA does not fall until y3-5,
            and the steady state is CBFWS pinned at CBFCRIT by infarction
            removing demand.",
    arch = "adult_isch", days = 3650),

  S03_natural_adult_haem = list(
    what = "Untreated adult, haemorrhagic phenotype (PCA involved).  Same
            equations; the periventricular route carries the load and gPVA,
            SIGPVA, ANEU and PHEM all rise with it.",
    arch = "adult_haem", days = 3650),

  S04_natural_sickle = list(
    what = "Sickle-cell moyamoya, Hb 8 g/dL.  The vascular lesion is ordinary;
            CBFCRIT is 36.96 instead of 19.71, so the SAME flow is ischaemic.",
    arch = "sickle_cell", days = 1825),

  S05_natural_asymptomatic = list(
    what = "Slowly progressive / incidentally found disease.  Reserve is still
            positive at y5 -- the AMORE-type patient.",
    arch = "asymptomatic", days = 3650),

  ## --- diagnostics: the two probes are not the same --------------------
  S06_acetazolamide_challenge = list(
    what = "Acetazolamide 1000 mg IV at day 1460 (reserve lost).  Territory B
            answers +39%; territory A answers -1.7% -- STEAL, because B is A's
            supply.  This is the acetazolamide-SPECT stage-II signature.",
    arch = "adult_isch", days = 1460.2, delta = 0.002,
    drugs = list(AZ = list(amt = 1000, start = 1460))),

  S07_hypercapnia_challenge = list(
    what = "PaCO2 40 -> 50 mmHg.  Opens the pial CONDUITS as well as the
            arteriole, so territory A gains +4.4% where acetazolamide lost
            1.7%.  An acetazolamide study and a breath-hold study are NOT
            interchangeable in moyamoya.",
    arch = "adult_isch", days = 1460, pars = list(PACO2 = 50, PROG = 0)),

  S08_crying_child = list(
    what = "Hyperventilation to PaCO2 25 mmHg -- crying, blowing, wind
            instrument.  Demand falls (protective) while the conduits constrict
            (harmful).  In a hemisphere with reserve this is safe; in a
            decompensated one it crosses CBFCRIT.",
    arch = "paediatric", days = 1460, pars = list(PACO2 = 25, PROG = 0)),

  ## --- the blood-pressure paradox --------------------------------------
  S09_hypotension_ischaemic = list(
    what = "MAP -20% in a pressure-passive ischaemic hemisphere.  dCBF/dMAP is
            ~0.5 mL/100g/min per mmHg here, so this is a stroke.",
    arch = "adult_isch", days = 2555, map_mult = 0.80),

  S10_bp_lowering_haemorrhagic = list(
    what = "The SAME MAP reduction in the haemorrhagic phenotype, where it
            lowers periventricular wall stress.  Opposite sign of benefit from
            identical equations.",
    arch = "adult_haem", days = 2555, map_mult = 0.80,
    drugs = list(AHT = list(amt = 10, ii = 1, addl = 2554))),

  ## --- surgery ---------------------------------------------------------
  S11_direct_bypass_isch = list(
    what = "Direct STA-MCA at y5 in the ischaemic phenotype.  gBYP is a
            PARALLEL CONDUCTANCE -- the only intervention in the model that
            changes gS.",
    arch = "adult_isch", days = 3650,
    surgery = list(kind = "direct", day = 1825)),

  S12_direct_bypass_haem_JAM = list(
    what = "Direct bypass at y7 in the haemorrhagic phenotype: the JAM-trial
            emulation.  Watch QPVA collapse as PA rises -- the operation
            RETIRES the fragile vessel rather than reinforcing it.  Predicted
            rebleed HR ~0.36 vs the observed 0.355.",
    arch = "adult_haem", days = 2555 + 1825,
    surgery = list(kind = "direct", day = 2555)),

  S13_indirect_child = list(
    what = "Indirect EDAS in a child.  VEGF- and ANGIO-gated ingrowth reaches
            most of a direct graft within a year.",
    arch = "paediatric", days = 2555,
    surgery = list(kind = "indirect", day = 1460)),

  S14_indirect_adult = list(
    what = "The SAME indirect operation in an adult.  One parameter (ANGIO)
            separates a near-cure from a near-miss.",
    arch = "adult_isch", days = 2555,
    surgery = list(kind = "indirect", day = 1460)),

  S15_combined_adult = list(
    what = "Combined direct + indirect in an adult -- the immediate graft plus
            whatever angiogenesis the patient can still supply.",
    arch = "adult_isch", days = 2555,
    surgery = list(kind = "combined", day = 1460)),

  S16_hyperperfusion = list(
    what = "Direct bypass onto a hemisphere whose arteriole has been on its
            floor for 5 y (REMOD ~ 0.59).  Normal pressure arrives at a vessel
            that cannot constrict: HYPER peaks ~day 2-7 and resolves over 2-3
            weeks, with BBB and EDEMA following.  A reserve-intact hemisphere
            does not do this at all (compare S17).",
    arch = "adult_isch", days = 3285 + 90, delta = 0.05,
    surgery = list(kind = "direct", day = 3285)),

  S17_bypass_preserved_reserve = list(
    what = "Identical operation on a hemisphere that still has reserve.  No
            hyperperfusion, and little to gain -- the model's account of why
            asymptomatic reserve-intact disease is a different decision.",
    arch = "asymptomatic", days = 1095 + 90, delta = 0.05,
    surgery = list(kind = "direct", day = 1095)),

  S18_postop_bp_control = list(
    what = "The sign flip: the SAME hemisphere that wanted a HIGHER MAP before
            surgery wants a lower one for two weeks afterwards.  Nothing about
            the patient changed except which side of gS* they are on.",
    arch = "adult_isch", days = 3285 + 60, delta = 0.05, map_mult = 0.85,
    surgery = list(kind = "direct", day = 3285)),

  ## --- medical therapy -------------------------------------------------
  S19_aspirin = list(
    what = "Aspirin 100 mg/d.  Acts on the EMBOLIC hazard from low-shear
            perforators; gS is untouched.",
    arch = "adult_isch", days = 1825,
    drugs = list(ASA = list(amt = 100, ii = 1, addl = 1824))),

  S20_cilostazol = list(
    what = "Cilostazol 100 mg b.d.  Antiplatelet plus a mild dilatory stimulus
            that a floored arteriole cannot use.",
    arch = "adult_isch", days = 1825,
    drugs = list(CILO = list(amt = 100, ii = 0.5, addl = 3649))),

  S21_nifedipine_steal = list(
    what = "Nifedipine GITS 60 mg/d -- the instructive failure.  A real
            cerebral vasodilator that cannot dilate a vessel already on its
            floor, that dilates the DONOR territory instead, and that lowers
            the MAP a pressure-passive territory is living on.",
    arch = "adult_isch", days = 1825,
    drugs = list(NIF = list(amt = 60, ii = 1, addl = 1824))),

  S22_statin = list(
    what = "Atorvastatin 20 mg/d -- the only agent in the model that touches
            the LESION (K_SMC), and even so it slows rather than stops it.",
    arch = "adult_isch", days = 1825,
    drugs = list(STAT = list(amt = 20, ii = 1, addl = 1824))),

  S23_minocycline = list(
    what = "Minocycline 100 mg b.d.  MMP-9 inhibition raises perforator wall
            repair, so it moves PHEM without moving CBF.",
    arch = "adult_haem", days = 2555,
    drugs = list(MIN = list(amt = 100, ii = 0.5, addl = 5109))),

  S24_transfusion_sickle = list(
    what = "Sickle-cell moyamoya transfused to Hb 11.  No vessel is touched;
            CBFCRIT falls from 36.96 to 26.88 and the infarct stops.  This is
            the model's account of STOP/SIT.",
    arch = "sickle_cell", days = 1825, pars = list(HB = 11.0))
)

## ===========================================================================
##  RUNNERS
## ===========================================================================
run_scenario <- function(nm) {
  s <- SCENARIOS[[nm]]
  run_mmd(arch = s$arch, days = s$days, delta = s$delta %||% 1,
          drugs = s$drugs, surgery = s$surgery, pars = s$pars,
          map_mult = s$map_mult %||% 1.0) %>% mutate(scenario = nm)
}

run_all <- function() bind_rows(lapply(names(SCENARIOS), run_scenario))

## ---------------------------------------------------------------------------
##  VALIDATION TABLE — reproduce the numbers in mmd_reference_output.txt.
##  Run this after any parameter change.  Expected values are quoted from the
##  Python reference implementation, which is the numerical truth for this
##  model; agreement to ~1% is expected (both solve the same equations, but
##  LSODA and mrgsolve's ODEPACK build differ in step control).
## ---------------------------------------------------------------------------
validate <- function() {
  chk <- function(label, got, want, tol = 0.02) {
    ## relative tolerance, except that a target of exactly zero is checked
    ## absolutely (otherwise the test is unsatisfiable)
    lim <- if (isTRUE(all.equal(want, 0))) max(tol, 1e-6) else tol * abs(want)
    ok <- is.finite(got) && abs(got - want) <= lim
    cat(sprintf("  %-46s got %10.3f  want %10.3f  %s\n",
                label, got, want, if (ok) "OK" else "** CHECK **"))
    invisible(ok)
  }
  cat("\n== derived constants ==\n")
  b <- run_mmd("adult_isch", days = 1, delta = 1)
  chk("CBF_crit, Hb 15 (mL/100g/min)", b$CBFCRIT[1], 19.71)
  chk("baseline CBF_A (mL/100g/min)",  b$CBFA[1],    50.00)
  chk("baseline OEF_A",                b$OEFA[1],     0.335)
  chk("baseline CMRO2_A",              b$CMROA[1],    3.300)
  chk("baseline P_A (mmHg)",           b$PA[1],      70.01)
  chk("baseline dCBF/dMAP",            b$DCBFDMAP[1], 0.000, tol = 1e-3)
  s <- run_mmd("sickle_cell", days = 1, delta = 1)
  chk("CBF_crit, Hb 8 (mL/100g/min)",  s$CBFCRIT[1], 36.96)

  cat("\n== adult-ischaemic natural history ==\n")
  o <- run_mmd("adult_isch", days = 3650, delta = 5)
  g <- function(d, v) approx(o$time, o[[v]], d)$y
  chk("y1  CVR_intrinsic (%)",  g(365,  "CVRINT"), 38.3, tol = 0.06)
  chk("y2  CVR_intrinsic (%)",  g(730,  "CVRINT"),  0.0, tol = 0.05)
  chk("y1  CBF_A",              g(365,  "CBFA"),   50.0)
  chk("y5  CBF_A",              g(1825, "CBFA"),   39.6, tol = 0.04)
  chk("y5  CBF_watershed",      g(1825, "CBFWS"),  21.5, tol = 0.06)
  chk("y10 CBF_watershed",      g(3650, "CBFWS"),  17.9, tol = 0.06)
  chk("y10 infarct fraction %", g(3650, "INFPCT"),  6.71, tol = 0.10)
  chk("y10 OEF_A",              g(3650, "OEFA"),    0.423, tol = 0.03)

  cat("\n== haemorrhagic phenotype ==\n")
  h <- run_mmd("adult_haem", days = 3650, delta = 5)
  gh <- function(d, v) approx(h$time, h[[v]], d)$y
  chk("y10 g_PVA (mL/min/mmHg)", gh(3650, "gPVA"), 0.524, tol = 0.06)
  chk("y10 P(haemorrhage)",      gh(3650, "PHEM"), 0.290, tol = 0.10)
  ph <- run_mmd("adult_isch", days = 3650, delta = 5)
  chk("y10 P(haem), PCA spared", approx(ph$time, ph$PHEM, 3650)$y, 0.031,
      tol = 0.20)

  cat("\n== the critical inlet conductance ==\n")
  cat("  gS* = 2.201 mL/min/mmHg: above it ARPOS < 1 and DCBFDMAP = 0,\n")
  cat("  below it ARPOS = 1 and DCBFDMAP > 0.  Check both sides:\n")
  for (st in c(0.20, 0.24)) {
    q <- mod %>% param(ARCHETYPE$adult_isch) %>% init(init_state(ARCHETYPE$adult_isch)) %>%
      init(SMC = st / 0.945) %>% param(PROG = 0) %>% mrgsim(end = 1, delta = 1) %>% as_tibble()
    cat(sprintf("    STEN %.2f -> gS %.3f  ARPOS %.3f  dCBF/dMAP %.3f  CVRmeas %+.1f%%\n",
                st, q$gS[1], q$ARPOS[1], q$DCBFDMAP[1], q$CVRMEAS[1]))
  }
  invisible(TRUE)
}

## ---------------------------------------------------------------------------
##  JAM-trial emulation: the haemorrhage HR is a PREDICTION.  Only the
##  CONSERVATIVE 5-y rebleeding rate (31.6%) was used to set HEM_HAZ0.
## ---------------------------------------------------------------------------
jam_emulation <- function(op_day = 2555, follow = 1825) {
  out <- lapply(c(conservative = 0, bypass = 1), function(k) {
    p <- ARCHETYPE$adult_haem
    m <- mod %>% param(p) %>% init(init_state(p)) %>%
      param(SURG_KIND = k, SURG_T = op_day)
    as_tibble(mrgsim(m, end = op_day + follow, delta = 5))
  })
  haz <- sapply(out, function(d) {
    i0 <- which.min(abs(d$time - op_day))
    tail(d$HEMH, 1) - d$HEMH[i0]
  })
  isc <- sapply(out, function(d) {
    i0 <- which.min(abs(d$time - op_day))
    tail(d$ISCH, 1) - d$ISCH[i0]
  })
  cat("\n== JAM-trial emulation (haemorrhagic phenotype, 5 y) ==\n")
  cat(sprintf("  5-y rebleeding  conservative %.3f   bypass %.3f   HR %.3f\n",
              1 - exp(-haz[["conservative"]]), 1 - exp(-haz[["bypass"]]),
              haz[["bypass"]] / haz[["conservative"]]))
  cat(sprintf("  5-y ischaemic   conservative %.3f   bypass %.3f   HR %.3f\n",
              1 - exp(-isc[["conservative"]]), 1 - exp(-isc[["bypass"]]),
              isc[["bypass"]] / isc[["conservative"]]))
  cat("  observed (Miyamoto 2014, PMID 24668203): 31.6% vs 11.9%, HR 0.355\n")
  invisible(out)
}

## ---------------------------------------------------------------------------
##  MAP optimisation: the two phenotypes have different optima from the SAME
##  equations.  This is the blood-pressure paradox as calculus.
## ---------------------------------------------------------------------------
map_optimum <- function(arch = "adult_isch", day = 2555,
                        maps = seq(65, 125, by = 5)) {
  p <- ARCHETYPE[[arch]]
  base <- mod %>% param(p) %>% init(init_state(p)) %>%
    mrgsim(end = day, delta = 5) %>% as_tibble()
  st <- as.numeric(tail(base, 1)[names(init_state(p))])
  names(st) <- names(init_state(p))
  res <- lapply(maps, function(mp) {
    q <- mod %>% param(p) %>% init(st) %>%
      param(MAP_MULT = mp / 90) %>% param(PROG = 0) %>%
      mrgsim(end = 1, delta = 1) %>% as_tibble()
    ## read the hazard constants from the model rather than hard-coding them,
    ## so this reporting helper cannot drift away from the ODE it summarises
    pp <- as.list(param(mod))
    sd <- pp$CBF_SD0 * (1 + pp$K_PASS * (q$ARPOS[1] > 0.99))
    ih <- pp$ISCH_HAZ0 *
      exp(-pmax(pmin((q$CBFWS[1] - q$CBFCRIT[1]) / sd, 25), -6)) +
      pp$ISCH_EMB * pp$K_EMB * st[["GMOYA"]]
    hh <- pp$HEM_HAZ0 * st[["ANEU"]] * q$SIGPVA[1]^pp$HEM_POW
    tibble(MAP = mp, CBFA = q$CBFA[1], CBFWS = q$CBFWS[1],
           SIGPVA = q$SIGPVA[1], isch = ih, haem = hh, total = ih + hh)
  })
  r <- bind_rows(res)
  cat(sprintf("\n== MAP optimum, %s at day %d ==\n", arch, day))
  print(as.data.frame(round(r, 4)))
  cat(sprintf("  hazard-minimising MAP = %d mmHg\n", r$MAP[which.min(r$total)]))
  invisible(r)
}

## ---------------------------------------------------------------------------
##  Virtual population: what actually predicts the outcome?
## ---------------------------------------------------------------------------
virtual_population <- function(n = 200, days = 1825, seed = 20260806) {
  set.seed(seed)
  idata <- tibble(
    ID      = seq_len(n),
    RNF     = sample(c(0, 1, 1.6), n, TRUE, prob = c(0.20, 0.72, 0.08)),
    ANGIO   = pmin(pmax(rnorm(n, 0.62, 0.22), 0.12), 1.05),
    K_SMC   = pmin(pmax(rlnorm(n, log(2.6e-4), 0.42), 5e-5), 9e-4),
    HB      = pmin(pmax(rnorm(n, 13.8, 1.5), 7.0), 16.5),
    MAP0    = pmin(pmax(rnorm(n, 92, 11), 68), 125),
    PCA_INV = sample(c(0, 1), n, TRUE, prob = c(0.70, 0.30)),
    ## pial coupling: the least well constrained parameter in the model, and
    ## the one that decides whether a bypass hyperperfuses (see README S8)
    G_LEAK  = pmin(pmax(rlnorm(n, log(0.80), 0.55), 0.12), 3.0)
  ) %>% mutate(SEC = ifelse(RNF == 0, 1.0, 0.0))
  out <- mod %>% idata_set(idata) %>%
    init(init_state(list(RNF = 1, SEC = 0, K_PDGF = 0.9))) %>%
    mrgsim(end = days, delta = 25) %>% as_tibble()
  last <- out %>% group_by(ID) %>% slice_tail(n = 1) %>% ungroup() %>%
    left_join(idata, by = "ID")
  cat("\n== virtual population, 5 y ==\n")
  for (v in c("INFPCT", "PISCH", "PHEM")) {
    cat(sprintf("  correlates of %s:\n", v))
    for (x in c("RNF", "ANGIO", "K_SMC", "HB", "MAP0", "PCA_INV",
                "STEN", "gS", "CBFWS", "CVRINT")) {
      r <- suppressWarnings(cor(last[[x]], last[[v]]))
      if (is.finite(r) && abs(r) > 0.12)
        cat(sprintf("      %-8s r = %+.3f\n", x, r))
    }
  }
  invisible(last)
}

## ---------------------------------------------------------------------------
if (identical(environment(), globalenv()) && !interactive()) {
  validate()
  jam_emulation()
  map_optimum("adult_isch", 2555)
  map_optimum("adult_haem", 3285)
  invisible(virtual_population(120))
}
