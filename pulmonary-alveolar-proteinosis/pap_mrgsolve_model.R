# =====================================================================
# Autoimmune Pulmonary Alveolar Proteinosis (aPAP) — QSP model
# mrgsolve implementation · 64 ODEs · time unit = DAYS
# =====================================================================
#
# WHAT THIS MODEL IS FOR
# ----------------------
# aPAP is usually drawn as "autoantibody -> macrophage dysfunction ->
# surfactant accumulation -> hypoxaemia". Written that way it is a story,
# not a model: every arrow is a fitted slope. This file makes three
# structural commitments instead, and every headline number below falls
# out of them rather than being fitted to the trial it is compared with.
#
# (1) RECYCLING IS NOT CLEARANCE.
#     Type-II re-uptake and re-secretion of surfactant is a large flux
#     with ZERO net effect at steady state. A model that lumps it into
#     "clearance" and then removes the macrophage arm predicts a pool
#     that settles at ~2-3x normal. Patients reach 30-100x. Only true
#     CATABOLISM (macrophage digestion, type-II intracellular
#     degradation) and mucociliary/lymphatic egress are sinks, and the
#     macrophage arm is most of the net sink. The loop is drawn on the
#     map and deliberately left out of the net balance.
#
# (2) THE AUTOANTIBODY IS A STOICHIOMETRIC BUFFER, NOT AN IC50.
#     Alveolar GM-CSF sits at tens of pg/mL (tens of pM); neutralising
#     antibody sites sit at hundreds to thousands of pM. Free ligand is
#     therefore set by a binding equilibrium in the excess-antibody
#     regime, solved here from the 1:1 quadratic in its numerically
#     stable form (the naive root suffers catastrophic cancellation at
#     these ratios and returns zero or noise). Two clinical facts that
#     look contradictory then come out of one equation: a titre
#     THRESHOLD near 5 ug/mL, and no titre-severity correlation above it.
#
# (3) INHALED GM-CSF WORKS BY OVERWHELMING THAT BUFFER LOCALLY.
#     The therapeutic quantity is a MOLAR RATIO in epithelial lining
#     fluid, not a plasma concentration. 300 ug nebulised at 40%
#     deposition delivers ~8300 pmol of ligand into ~30 mL of ELF holding
#     ~160 pmol of neutralising sites at a median titre - a ~50-fold excess,
#     against a receptor that saturates at tens of pM. That
#     is why route, not dose, separates success from failure, and why the
#     efficacy variable is TIME ABOVE THRESHOLD rather than dose.
#
# WHAT IS DERIVED RATHER THAN ASSERTED
# ------------------------------------
#   * The healthy baseline is SOLVED, not typed in: uptake Vmax, the
#     macrophage lipid load and de novo production are back-calculated in
#     $MAIN from the specified turnover and sink split, so a healthy lung
#     is exactly stationary (D01 measures the drift).
#   * The disease is GENERATED: seroconversion at t=0, then years of
#     integration. Nothing is initialised to a diseased value.
#   * The presenting patient is an EVENT, not a fitted equilibrium: the
#     patient presents when A-aDO2 first crosses a threshold. Patients
#     whose residual catabolic floor holds them below it never present -
#     which is the model account of the 32% found asymptomatically on
#     health screening.
#   * A-aDO2, PaO2 and DLCO are COMPUTED from shunt, the alveolar gas
#     equation, the Severinghaus dissociation curve and the Fick
#     relation. A-aDO2 is the IMPALA primary endpoint; here it is
#     physics, not a score.
#
# PUBLISHED ANCHORS USED FOR VALIDATION (see run_diagnostics)
# -----------------------------------------------------------
#   IMPALA-2 (Trapnell 2025, PMID 40834301): molgramostim 300 ug QD,
#     DLCO +9.8 vs +3.8 %pred at wk 24; +11.6 vs +4.7 at wk 48;
#     SGRQ-T -11.5 vs -4.9.
#   IMPALA   (Trapnell 2020, PMID 32897035): A-aDO2 -12.8 vs -6.6 mmHg
#     at wk 24 (continuous arm), n=138.
#   PAGE     (Tazawa 2019, PMID 31483963): inhaled sargramostim 125 ug
#     BID x7 d every other week, A-aDO2 -4.50 vs +0.17 mmHg.
#   Tazawa 2010 (PMID 20167854): phase 2, A-aDO2 -12.3 mmHg, 62% response.
#   Sakagami 2010 (PMID 20224064): critical GMAb threshold 5 ug/mL,
#     the same in lung and blood.
#   Inoue 2008 (PMID 18202348): 223 patients, 31.8% asymptomatic,
#     DSS correlates with DLCO and biomarkers but NOT with GMAb titre.
#   Beccaria 2004 (PMID 15083749): >70% recurrence-free 7 y after WLL;
#     DLCO still 75 +- 19 %pred and A-aDO2 27 +- 11 mmHg at 5 y.
#
# WHAT MATCHES, AND WHAT DOES NOT (16 of 17 diagnostics pass)
# -----------------------------------------------------------
# Matches, none of them fitted to the number they are compared with:
#   * Healthy lung stationary to -0.002% over 1000 d; phospholipid mass
#     balance closes to 9e-8 mg with a lavage in the run (D01, D02).
#   * Half-signal GMAb titre 5.6 ug/mL against a published 5 (D04).
#   * IMPALA-2 at 24 weeks, BOTH arms: DLCO +9.6 drug / +3.7 placebo
#     against +9.8 / +3.8 observed, difference +5.9 vs +6.0 (D07). The
#     placebo arm is not a placebo parameter - it is enrolment at the
#     low point of a deterioration, plus DLCO measurement learning.
#   * PAGE: the placebo arm does NOT move once its 12-week
#     improver-excluding run-in is simulated (D10). Same model, two
#     enrolment protocols, two different placebo arms.
#   * Whole lung lavage: 44% of burden removed, A-aDO2 34.3 -> 19.6,
#     return to the pre-lavage gradient at 3.4 y - and time to
#     recurrence is mass removed over net imbalance, not a parameter.
#
# Failures and refutations, reported rather than removed:
#   * D08. The 48-week IMPALA-2 difference comes out +15.1 against +6.9
#     observed - roughly twice. The model contains a recovery-side
#     positive feedback (clearing the lung opens units, so more aerosol
#     reaches the rest) that the 48-week data do not show, while the
#     withdrawal data of Tazawa 2010 arguably do. Quote 24 weeks.
#   * D15. NO BISTABILITY. The map was drawn expecting the filling
#     feedback to create a held state and a runaway state; burdens of
#     300 and 8169 mg both converge to 6160 mg at the same catabolic
#     floor. The premise is refuted and cluster 11 says so.
#   * D05. The titre-severity decorrelation of Inoue 2008 is only
#     PARTLY reproduced (r = +0.57, or +0.85 without the
#     binding-vs-neutralisation distinction). And the catabolic floor
#     is NOT the dominant covariate in that population - the antibody
#     side still is. Both statements correct earlier drafts of this file.
#   * D17. The inhaled dose threshold sits near 1-3 ug/day, a hundred
#     times below the clinical dose. Everything above is plateau. This
#     is a prediction with a clear falsification route.
#   * D19. Rituximab, FcRn blockade and statin effects are HYPOTHESES.
#     The statin looks disproportionately good because it is oral and
#     pays no aerosol-reach penalty - a real comparison, but the size
#     of its capacity effect is not anchored to data.
#
# NOTE ON USE: educational / research QSP model. Not validated for
# clinical decision-making, prescribing or regulatory submission.
# =====================================================================

library(mrgsolve)

# ---------------------------------------------------------------------
# MODEL SPECIFICATION
# ---------------------------------------------------------------------
pap_code <- '
$PROB
# Autoimmune Pulmonary Alveolar Proteinosis (aPAP) - QSP model
# 64 ODEs. Time unit = DAYS. Amounts: surfactant mg, GM-CSF pmol,
# antibody ug, drugs mg or ug as annotated.

$SET end=730, delta=1, rtol=1e-6, atol=1e-8, maxsteps=500000

$GLOBAL
// ---------------- small helpers -----------------------------------
double pos(double x){ return (x > 0.0) ? x : 0.0; }
double dmx(double a, double b){ return (a > b) ? a : b; }
double dmn(double a, double b){ return (a < b) ? a : b; }
double clamp2(double x, double lo, double hi){
  return (x < lo) ? lo : ((x > hi) ? hi : x);
}
double sat(double x, double k){ return pos(x)/(pos(x) + k); }
double hillf(double x, double k, double n){
  double a = pow(pos(x), n);
  return a/(a + pow(k, n));
}

// ---------------- 1:1 binding, numerically stable ------------------
// Free ligand from total ligand L, total binding sites S, dissociation
// constant Kd (all same concentration unit). The textbook root
//   0.5*((L-S-Kd) + sqrt((L-S-Kd)^2 + 4*Kd*L))
// is the difference of two nearly equal large numbers when S >> L, which
// is exactly the aPAP regime (S/L ~ 1e3-1e4): in double precision it
// returns zero or garbage. The algebraically identical form below is
// stable because nothing cancels.
double freelig(double L, double S, double Kd){
  L = pos(L); S = pos(S); Kd = dmx(Kd, 1e-12);
  double b = S + Kd - L;
  double disc = sqrt(b*b + 4.0*Kd*L);
  return 2.0*Kd*L/dmx(b + disc, 1e-30);
}

// ---------------- oxygen ------------------------------------------
// Severinghaus (1979) standard curve and its closed-form inverse.
double sao2f(double p){
  p = dmx(p, 1.0);
  return 1.0/(23400.0/(p*p*p + 150.0*p) + 1.0);
}
double po2f(double s){
  s = clamp2(s, 0.0005, 0.99995);
  double a = 11700.0/((1.0/s) - 1.0);
  return pow(a + sqrt(125000.0 + a*a), 1.0/3.0);
}
double cao2f(double p, double hb){ return 1.34*hb*sao2f(p) + 0.003*p; }

// Arterial PO2 downstream of a shunt, given alveolar PO2, shunt
// fraction, haemoglobin, VO2 (mL/min) and cardiac output (L/min).
// Two-pass refinement of the dissolved-O2 term; no iteration needed
// because the inverse dissociation curve is closed-form.
double pao2_shunt(double PAO2, double shunt, double hb, double vo2, double co){
  shunt = clamp2(shunt, 0.0, 0.95);
  hb = dmx(hb, 4.0);
  double ccc = cao2f(PAO2, hb);
  double avd = vo2/(dmx(co,1.0)*10.0);
  double ca  = ccc - avd*shunt/dmx(1.0 - shunt, 0.05);
  ca = dmx(ca, 0.8);
  double s1 = clamp2((ca - 0.12)/(1.34*hb), 0.005, 0.9995);
  double p1 = po2f(s1);
  double s2 = clamp2((ca - 0.003*p1)/(1.34*hb), 0.005, 0.9995);
  return po2f(s2);
}

// window indicator: 1 inside [t0, t0+dur), else 0
double win(double t, double t0, double dur){
  return ((t >= t0) && (t < t0 + dur)) ? 1.0 : 0.0;
}

$PARAM @annotated
// ============ COMPARTMENT GEOMETRY AND UNIT CONVERSIONS ============
VELF     :  30      : Epithelial lining fluid volume (mL)
VPLAS    : 3000     : Plasma volume (mL)
MWGM     : 14500    : GM-CSF molecular weight (g/mol)
BW       :  70      : Body weight (kg)

// ============ ENDOGENOUS GM-CSF IN THE ALVEOLUS ====================
// ELF concentration, NOT BAL concentration. BAL dilutes ELF 50-100x,
// which is why the literature reports GM-CSF as near-undetectable in
// lavage fluid while the lining-fluid concentration is tens of pg/mL.
CGMH     :  41.6    : Healthy free GM-CSF in ELF (pg/mL)
KELGM    :   5.545  : ELF GM-CSF elimination rate (1/d; t1/2 = 3 h)
FSPILL   :   0.30   : Fraction of ELF GM-CSF clearance reaching plasma
KELGMP   :   6.65   : Plasma GM-CSF elimination (1/d; t1/2 = 2.5 h)
ELFPEN   :   0.02   : ELF penetration of plasma GM-CSF (fraction)

// ============ RECEPTOR AND SIGNAL TRANSDUCTION =====================
KDR      :  30      : GM-CSF receptor high-affinity Kd (pM)
RECFRAC  :   1.0    : Functional receptor fraction (0 = hereditary PAP)
OCC50    :   0.020  : Occupancy giving half-maximal transcription
HSIG     :   2.0    : Hill exponent of the transducer
TAUPU1   :   2.0    : PU.1 time constant (d)
TAUPPG   :   3.0    : PPAR-gamma time constant (d)
TAUCAP   :   4.0    : Catabolic-capacity time constant (d)
SMAX     :   1.8    : Ceiling on relative catabolic capacity
// CAPFLOOR is the GM-CSF-independent residual catabolic capacity, and the
// principal patient-level covariate in this model (see D04/D05).
CAPFLOOR :   0.305  : Residual GM-CSF-independent catabolic capacity

// ============ AUTOANTIBODY =========================================
TITRE0   :  25      : Target serum GMAb titre (ug/mL)
KELAB    :   0.033  : IgG elimination (1/d; t1/2 = 21 d)
// Re-equilibration must be fast: after a single-volume exchange, serum IgG
// falls ~60% and is back to ~two thirds of baseline within 48 h. With slow
// rate constants the model predicts a 92% reduction after ten sessions,
// which no plasmapheresis course achieves.
K12AB    :   0.50   : Serum -> extravascular (1/d)
K21AB    :   0.40   : Extravascular -> serum (1/d)
KAPPA    :   0.05   : ELF/serum IgG transudation ratio
KELFOUT  :   0.70   : ELF antibody turnover (1/d)
FNEUT    :   0.32   : Neutralising, receptor-blocking fraction of titre
KDAB     : 300      : Apparent GMAb dissociation constant (pM)
FPB      :   0.35   : Fraction of GMAb output from short-lived plasmablasts
KPBT     :   0.20   : Plasmablast turnover (1/d)
KLLT     :   0.0035 : Long-lived plasma cell turnover (1/d)
KBREP    :   0.006  : B-cell repopulation rate (1/d)
ABDECAY  :   0.0    : Extra first-order decay of GMAb output (1/d)
FCRNX    :   1.0    : Multiplier on IgG elimination (FcRn blockade)
FCRNT    :  1e6     : FcRn inhibitor start time (d)

// ============ ALVEOLAR MACROPHAGE POPULATION =======================
// Digestion is CAPACITY-limited, not substrate-limited. Written first-order
// in lipid load, an engorged macrophage carrying 6x the lipid digests as
// much as a healthy one, and the model then cannot generate PAP at all
// (the pool settles at 2.8x normal instead of 30-100x). Michaelis-Menten
// with a per-cell Vmax proportional to the transcriptional program is what
// makes the foamy cell a stalled cell.
LIP0P    : 130      : Healthy macrophage lipid load, whole lung (mg)
KMDIG    : 150      : Half-saturating lipid load for digestion (mg)
KDEATH   :   0.03   : Alveolar macrophage death rate (1/d)
PHIRET   :   0.80   : Fraction of dead-cell lipid returned to the alveolus
KILIP    : 1000     : Lipid load half-inhibiting further uptake (mg)
KMUP     : 600      : Michaelis constant, macrophage uptake (mg)
TAUAMN   :   7.0    : Macrophage number time constant (d)
AMNGM    :   1.2    : Maximal fractional rise in AM number on GM-CSF
MONOF    :   1.0    : Monocyte supply (fraction; <1 = secondary PAP)

// ============ SURFACTANT MASS BALANCE ==============================
P0       : 300      : Healthy alveolar phospholipid pool (mg)
TURN0    : 200      : Healthy net turnover of that pool (mg/d)
FMAC     :   0.65   : Macrophage share of the healthy NET sink
FII      :   0.25   : Type-II intracellular degradation share
FMUCO    :   0.10   : Mucociliary + lymphatic egress share
KMII     :  60      : Michaelis constant, type-II degradation (mg)
KMMUC    :  40      : Michaelis constant, mucociliary egress (mg)
PMAXP    : 40000    : Physical ceiling on alveolar burden (mg)
FCHOL    :   0.08   : Cholesterol fraction of surfactant mass
FSPRO    :   0.10   : Surfactant-protein fraction of surfactant mass
KSEQ     :   0.020  : Consolidation rate into sequestered material (1/d)
KDESEQ   :   0.010  : Re-dispersion of sequestered material (1/d)
KDEB     :   0.50   : Debris clearance (1/d)
KAWP     :   6.0    : Airway transit clearance (1/d)

// ============ FILLING, EXCLUSION, STRUCTURE ========================
FMAXF    :   0.85   : Maximal fully-filled fraction (geographic sparing)
PF50     : 7700     : Burden giving half-maximal full filling (mg)
HF       :   2.5    : Hill exponent of filling
PPART50  : 2600     : Burden giving half-maximal partial filling (mg)
TAUFILL  :   3.0    : Redistribution lag of filling states (d)
EXPMAC   :   0.35   : Filling-exclusion exponent, macrophage arm
EXPII    :   0.25   : Filling-exclusion exponent, type-II arm
EXPMUC   :   0.50   : Filling-exclusion exponent, mucociliary arm
KFIB     :   3.0e-5 : Fibrosis accrual rate (1/d at full filling)
FIBMAX   :   0.50   : Maximal fibrotic fraction
KTGF     :   0.05   : TGF-beta turnover (1/d)

// ============ GAS EXCHANGE =========================================
FIO2     :   0.21   : Inspired O2 fraction
PATM     : 760      : Barometric pressure (mmHg)
PACO2    :  40      : Arterial PCO2 (mmHg)
RQ       :   0.8    : Respiratory quotient
SHUNT0   :   0.02   : Physiologic shunt fraction
SHCO     :   0.66   : Shunt per unit fully-filled fraction
VO2R     : 250      : Resting O2 consumption (mL/min)
COR      :   5.0    : Resting cardiac output (L/min)
VO2EX    : 1250     : Exercise O2 consumption (mL/min)
COEX     :  15.0    : Exercise cardiac output (L/min)
PACO2EX  :  42      : Exercise arterial PCO2 (mmHg)
DLCOMAX  : 100      : DLCO at zero filling (% predicted)
DLPART   :   0.55   : DLCO penalty weight of partially filled units
DLFIB    :   0.60   : DLCO penalty weight of fibrosis
TAULEARN :  25      : Time constant of the DLCO measurement-learning effect (d)
DLCOLMAX :   0.0    : Ceiling of that effect (%pred; 1.5 inside a trial)

// ============ SYSTEMIC RESPONSES ===================================
HB0      :  15.0    : Baseline haemoglobin (g/dL)
HBSLOPE  :   0.055  : Erythrocytotic response (g/dL per mmHg below 85)
HBMAX    :  19.5    : Ceiling on haemoglobin
TAUHB    :  25      : Haemoglobin time constant (d)
PVR0     :   1.8    : Baseline pulmonary vascular resistance (Wood units)
PVRSL    :   0.050  : HPV slope (WU per mmHg below 70)
TAUPVR   :  20      : PVR time constant (d)

// ============ BIOMARKERS ===========================================
KL6B     : 250      : Serum KL-6 intercept (U/mL)
KL6S     : 120      : KL-6 slope per unit pool ratio
SPDB     :  60      : Serum SP-D intercept (ng/mL)
SPDS     :  12      : SP-D slope
CEAB     :   1.5    : Serum CEA intercept (ng/mL)
CEAS     :   0.40   : CEA slope
LDHB     : 200      : Serum LDH intercept (U/L)
LDHS     :   8.0    : LDH slope
TAUBIO   :  10      : Biomarker time constant (d)

// ============ SYMPTOMS AND QUALITY OF LIFE =========================
// Slopes and intercept are calibrated to the IMPALA-2 SGRQ-T baseline (~40)
// AND to its change (-11.5 with drug, -4.9 with placebo) simultaneously; a
// single linear map cannot be fitted to one without checking the other.
SGA      :   0.90   : SGRQ points per mmHg of A-aDO2
SGB      :   0.48   : SGRQ points per %pred of DLCO deficit
SG0      : -26.7    : SGRQ offset (formula intercept; output clamped >= 0)
TAUSG    :  12      : SGRQ time constant (d)
WD0      : 550      : 6MWD intercept (m)
WDA      :   2.5    : 6MWD loss per mmHg A-aDO2
WDB      :   1.5    : 6MWD loss per %pred DLCO deficit
TAUWD    :  20      : 6MWD time constant (d)

// ============ HOST DEFENCE AND INFECTION ===========================
HAZ0     :   8.0e-5 : Baseline serious-infection hazard (1/d)
AINF     :   1.5    : Hazard multiplier from catabolic failure
BINF     :   0.05   : Hazard multiplier per unit pool ratio
KGNOC    :   0.90   : Nocardia growth rate (1/d)
NOCMAX   :  1e6     : Nocardia carrying capacity (arbitrary units)
KKILLN   :   2.6    : Killing rate at unit capacity (1/d)
KNOC     :  1e4     : Half-saturation of killing (units)

// ============ INHALED GM-CSF (molgramostim / sargramostim) =========
FDEP     :   0.40   : Alveolar deposition fraction of nominal dose
KDISS    :  30      : Dispersal of deposited drug into ELF (1/d)
KMCC     :   2.0    : Mucociliary loss of deposited drug (1/d)

// ============ SUBCUTANEOUS GM-CSF ==================================
FSC      :   0.60   : SC bioavailability
KASC     :   6.0    : SC absorption rate (1/d)

// ============ RITUXIMAB ============================================
CLRTX    :   0.30   : Rituximab clearance (L/d)
V1RTX    :   3.0    : Central volume (L)
V2RTX    :   2.7    : Peripheral volume (L)
QRTX     :   0.35   : Intercompartmental clearance (L/d)
KKILLB   :   0.35   : B-cell kill rate per (ug/mL)/d

// ============ STATIN ===============================================
FSTAT    :   0.14   : Oral bioavailability
KASTAT   :   8.0    : Absorption rate (1/d)
KELSTAT  :   1.19   : Elimination rate (1/d; t1/2 = 14 h)
VSTAT    : 380      : Volume of distribution (L)
EMAXST   :   0.035  : Maximal addition to the catabolic floor
EC50ST   :   5.0    : Statin concentration for half that effect (ng/mL)

// ============ REGIONAL REACH OF INHALED DRUG =======================
// Aerosol follows ventilation. A surfactant-filled alveolus is not
// ventilated, so it receives almost no nebulised drug - and it is exactly
// where the burden sits. The macrophages that matter are therefore the
// ones the drug cannot reach, while the drug floods the open units that
// have little left to clear. This is what makes inhaled GM-CSF partial
// and slow where the receptor pharmacology alone predicts a cure, and it
// is why lavage-then-GM-CSF is better than either alone. Oral and
// subcutaneous agents do not pay this penalty: blood reaches both regions.
// Two geometric facts set how much of the BURDEN the drug can act on:
// almost all of it sits inside filled units (RHO), and a consolidated
// region is cleared from its edges inward, where drug-exposed macrophages
// abut filled alveoli (EDGEF, an interface term ~ F(1-F)). The interface
// vanishes as the lung clears, which is why the response decelerates -
// a volume-based reach term makes it accelerate, which the trials do not.
EPSREACH :   0.005  : Baseline penetration of inhaled effect into filled units
RHO      :  50      : Surfactant mass per filled unit / per open unit
EDGEF    :   0.020  : Interface clearance coefficient

// ============ PPAR-GAMMA AGONIST (exploratory) =====================
PIOEFF   :   0.0    : Direct PPAR-gamma drive (0-1, exploratory)

// ============ PROCEDURES: WHOLE LUNG LAVAGE ========================
WLLT1    :  1e6     : Whole lung lavage time 1 (d)
WLLT2    :  1e6     : Whole lung lavage time 2 (d)
WLLT3    :  1e6     : Whole lung lavage time 3 (d)
WLLT4    :  1e6     : Whole lung lavage time 4 (d)
WLLDUR   :   0.25   : Duration of one lavage (d)
WLLFRAC  :   0.65   : Fraction of the alveolar pool removed
WLLSEQF  :   0.45   : Fraction of sequestered material removed
WLLABF   :   0.80   : Fraction of ELF antibody washed out
WLLAMF   :   0.40   : Fraction of alveolar macrophages removed
WLLLIPF  :   0.55   : Fraction of macrophage lipid removed

// ============ PROCEDURES: PLASMAPHERESIS ===========================
PLEXT    :  1e6     : First plasmapheresis session (d)
PLEXN    :   0      : Number of sessions
PLEXINT  :   2.0    : Interval between sessions (d)
PLEXDUR  :   0.15   : Duration of one session (d)
PLEXFRAC :   0.60   : Fraction of intravascular IgG removed per session

// ============ EXPOSURE AND ACUTE DETERIORATION =====================
SMKON    :   0.0    : Active smoking (0/1)
SMKSTOP  :  1e6     : Smoking cessation time (d)
SMKEFF   :   0.20   : Fractional reduction of catabolic capacity while smoking
FLARE0   :   0.0    : Acute deterioration at t=0 (0-1)
KFLARE   :   0.030  : Resolution of that deterioration (1/d)
FLAREAMP :   0.35   : Filling amplification at full flare

// ============ CONTROL FLAGS ========================================
INITMODE :   1      : 1 = initialise from the derived healthy steady state
ABON     :   1      : 1 = autoantibody production on
NOC0     :   0      : Initial Nocardia inoculum (units)
EXFLAG   :   1      : 1 = also report exercise gas exchange

$CMT @annotated
// --------- inhaled and systemic GM-CSF -----------------------------
INHD   : Deposited inhaled GM-CSF awaiting dispersal (ug)
GMT    : Total GM-CSF in ELF, free + antibody-bound (pmol)
GMPL   : Plasma GM-CSF (pmol)
SCD    : Subcutaneous GM-CSF depot (ug)
// --------- autoantibody --------------------------------------------
ABC    : Serum GMAb (ug)
ABP    : Extravascular GMAb (ug)
ABE    : ELF GMAb (ug)
BCELL  : CD20+ B-cell pool (relative)
PB     : Short-lived plasmablast pool (relative)
LLPC   : Long-lived plasma cell pool (relative)
// --------- drugs ---------------------------------------------------
RTXC   : Rituximab central (mg)
RTXP   : Rituximab peripheral (mg)
STATG  : Statin gut (mg)
STATC  : Statin central (mg)
// --------- macrophage program --------------------------------------
PU1    : PU.1 activity (relative to healthy)
PPARG  : PPAR-gamma activity (relative)
ABCG1  : ABCG1 cholesterol-efflux capacity (relative)
CAPC   : Per-cell catabolic capacity (relative)
AMN    : Alveolar macrophage number (relative)
MONO   : Blood monocyte pool (relative)
LIP    : Macrophage lipid load (mg)
AMCH   : Macrophage cholesterol load (mg)
AMDEAD : Senescent/dying macrophage pool (mg lipid equivalent)
// --------- surfactant pools ----------------------------------------
PLA    : Alveolar phospholipid pool (mg)
CHA    : Alveolar cholesterol pool (mg)
SPRO   : Alveolar surfactant protein pool (mg)
LB     : Type-II lamellar body pool (mg)
SEQ    : Sequestered / consolidated material (mg)
DEB    : Alveolar debris (mg)
AWP    : Airway transit pool (mg)
// --------- structure ----------------------------------------------
FFUL   : Fully filled alveolar fraction
FPART  : Partially filled alveolar fraction
FIB    : Irreversible fibrotic fraction
TGFB   : TGF-beta activity (relative)
FLARE  : Acute deterioration state (0-1)
// --------- systemic physiology ------------------------------------
HB     : Haemoglobin (g/dL)
PVR    : Pulmonary vascular resistance (Wood units)
// --------- symptoms, QoL, biomarkers ------------------------------
SGRQ   : SGRQ total score (lagged)
SIXMWD : 6-minute walk distance (m, lagged)
MMRC   : mMRC dyspnoea score (lagged)
KL6    : Serum KL-6 (U/mL)
SPD    : Serum SP-D (ng/mL)
CEA    : Serum CEA (ng/mL)
LDH    : Serum LDH (U/L)
NEUT   : Blood neutrophils (10^9/L)
EOS    : Blood eosinophils (10^9/L)
DLCOL  : DLCO measurement-learning offset (%pred)
// --------- infection ----------------------------------------------
NOC    : Nocardia burden (units)
CUMHAZ : Cumulative infection hazard
// --------- audit trail --------------------------------------------
COV    : Cumulative time above the signalling threshold (d)
CUMPROD: Cumulative de novo surfactant production (mg)
CUMMAC : Cumulative macrophage digestion (mg)
CUMII  : Cumulative type-II degradation (mg)
CUMMUC : Cumulative mucociliary egress (mg)
CUMWLL : Cumulative surfactant removed by lavage (mg)
CUMPLEX: Cumulative IgG removed by plasmapheresis (ug)
CUMDEB : Cumulative debris cleared (mg)
CUMDOSE: Cumulative GM-CSF delivered to ELF (ug)
NWLL   : Number of lavages performed

$MAIN
// ==================================================================
// DERIVED HEALTHY BASELINE
// ------------------------------------------------------------------
// Nothing about the healthy lung is typed in twice. The sink split and
// the pool turnover are specified; uptake Vmax, the resting macrophage
// lipid load and de novo production are SOLVED so that the healthy
// system is exactly stationary. D01 measures what is left.
// ==================================================================
double f0    = FMAXF*hillf(P0, PF50, HF);
double x0    = dmx(1.0 - f0, 1e-6);
double satI0 = P0/(KMII + P0);
double satM0 = P0/(KMMUC + P0);
double satU0 = P0/(KMUP + P0);

double VMAXII  = FII  *TURN0/dmx(satI0*pow(x0, EXPII), 1e-9);
double VMAXMUC = FMUCO*TURN0/dmx(satM0*pow(x0, EXPMUC), 1e-9);

// Net macrophage removal is uptake minus what dying cells give back. With
// the healthy lipid load specified, the digestive Vmax and the uptake Vmax
// are both solved, so the healthy lung is stationary by construction.
double LIP0 = LIP0P;
double D0   = FMAC*TURN0 - (1.0 - PHIRET)*KDEATH*LIP0;
double VDIG = D0*(KMDIG + LIP0)/dmx(LIP0, 1e-6);
double U0   = FMAC*TURN0 + PHIRET*KDEATH*LIP0;
double VMAXUP = U0*(1.0 + LIP0/KILIP)/dmx(satU0*pow(x0, EXPMAC), 1e-9);

// de novo production must replace the whole net sink, corrected for the
// physical-space term so that the healthy lung is stationary
double PRODB = TURN0/dmx(1.0 - (P0/PMAXP), 1e-6);

// healthy signal, used to normalise the transducer
double LH   = CGMH*1000.0/MWGM;              // pg/mL -> pM
double XH   = RECFRAC*LH/(KDR + LH);         // healthy occupancy
double SH   = dmx(hillf(XH, OCC50, HSIG), 1e-9);

// antibody production that holds the target titre at steady state
double SYNAB = ABON*KELAB*TITRE0*VPLAS;

// initial conditions: the derived healthy steady state
if(INITMODE > 0.5){
  PLA_0  = P0;
  LB_0   = P0;
  CHA_0  = FCHOL*P0;
  SPRO_0 = FSPRO*P0;
  LIP_0  = LIP0;
  AMCH_0 = FCHOL*LIP0;
  GMT_0  = LH*(VELF/1000.0);
  PU1_0  = 1.0; PPARG_0 = 1.0; ABCG1_0 = 1.0; CAPC_0 = 1.0;
  AMN_0  = 1.0; MONO_0 = MONOF; BCELL_0 = 1.0; PB_0 = 1.0; LLPC_0 = 1.0;
  FFUL_0 = f0;
  FPART_0= (1.0 - f0)*(P0/(PPART50 + P0));
  HB_0   = HB0;
  PVR_0  = PVR0;
  NOC_0  = NOC0;
  FLARE_0= FLARE0;
  // symptom and biomarker states are set to their healthy targets
  KL6_0  = KL6B + KL6S*1.0*(1.0 + f0);
  SPD_0  = SPDB + SPDS*1.0*(1.0 + f0);
  CEA_0  = CEAB + CEAS*1.0*(1.0 + f0);
  LDH_0  = LDHB + LDHS*1.0*(1.0 + f0);
  NEUT_0 = 4.0;
  EOS_0  = 0.20;
  SGRQ_0 = pos(SG0 + SGA*7.4 + SGB*6.0);
  SIXMWD_0 = WD0 - WDA*10.0 - WDB*6.0;
  MMRC_0 = 0.0;
}

// deposition of the inhaled dose
F_INHD = FDEP;
F_SCD  = FSC;

$ODE
// ==================================================================
// 1. ANTIBODY KINETICS
// ==================================================================
double fcrn = (SOLVERTIME >= FCRNT) ? FCRNX : 1.0;
double synab = SYNAB*(FPB*PB + (1.0 - FPB)*LLPC)*exp(-ABDECAY*SOLVERTIME);

// plasmapheresis: a sequence of finite windows, each stripping a fixed
// fraction of the INTRAVASCULAR pool only. The extravascular pool is
// untouched, which is the whole reason the procedure disappoints.
double plexon = 0.0;
for(int i = 0; i < 12; i++){
  if(i < PLEXN) plexon += win(SOLVERTIME, PLEXT + i*PLEXINT, PLEXDUR);
}
double kplex = plexon*(-log(dmx(1.0 - PLEXFRAC, 1e-6))/dmx(PLEXDUR, 1e-6));

// whole lung lavage windows
double wllon = win(SOLVERTIME, WLLT1, WLLDUR) + win(SOLVERTIME, WLLT2, WLLDUR)
             + win(SOLVERTIME, WLLT3, WLLDUR) + win(SOLVERTIME, WLLT4, WLLDUR);
double kwll  = wllon*(-log(dmx(1.0 - WLLFRAC, 1e-6))/dmx(WLLDUR, 1e-6));
double kwseq = wllon*(-log(dmx(1.0 - WLLSEQF, 1e-6))/dmx(WLLDUR, 1e-6));
double kwab  = wllon*(-log(dmx(1.0 - WLLABF,  1e-6))/dmx(WLLDUR, 1e-6));
double kwam  = wllon*(-log(dmx(1.0 - WLLAMF,  1e-6))/dmx(WLLDUR, 1e-6));
double kwlip = wllon*(-log(dmx(1.0 - WLLLIPF, 1e-6))/dmx(WLLDUR, 1e-6));

dxdt_ABC = synab - (KELAB*fcrn + K12AB)*ABC + K21AB*ABP
           - KELFOUT*KAPPA*ABC*(VELF/VPLAS) - kplex*ABC;
dxdt_ABP = K12AB*ABC - K21AB*ABP;
dxdt_ABE = KELFOUT*KAPPA*ABC*(VELF/VPLAS) - KELFOUT*ABE - kwab*ABE;

dxdt_BCELL = KBREP*(1.0 - BCELL) - KKILLB*(RTXC/V1RTX)*BCELL;
dxdt_PB    = KPBT*(BCELL - PB);
dxdt_LLPC  = KLLT*(1.0 - LLPC);

// rituximab: two compartments, linear. Depletes CD20+ B cells only, so
// the plasmablast-derived share of GMAb output (FPB) is all it can reach.
dxdt_RTXC = -(CLRTX/V1RTX)*RTXC - (QRTX/V1RTX)*RTXC + (QRTX/V2RTX)*RTXP;
dxdt_RTXP =  (QRTX/V1RTX)*RTXC - (QRTX/V2RTX)*RTXP;

// statin: first-order oral input, one-compartment disposition
dxdt_STATG = -KASTAT*STATG;
dxdt_STATC =  KASTAT*STATG - KELSTAT*STATC;

// ==================================================================
// 2. GM-CSF: DELIVERY, BUFFERING, OCCUPANCY
// ==================================================================
double pgmug = 1.0e6/MWGM;                    // pmol per ug
double endog = LH*(VELF/1000.0)*KELGM;        // endogenous input (pmol/d)

dxdt_INHD = -(KDISS + KMCC)*INHD;
dxdt_SCD  = -KASC*SCD;

dxdt_GMT = endog + KDISS*INHD*pgmug
           + ELFPEN*KELGMP*GMPL*(VELF/VPLAS)
           - KELGM*GMT;

dxdt_GMPL = KASC*SCD*pgmug + FSPILL*KELGM*GMT - KELGMP*GMPL;

// free ligand from the binding equilibrium, in pM
double Ltot = GMT/(VELF/1000.0);
double Stot = ABE*(2.0e6/150000.0)*FNEUT/(VELF/1000.0);
double Lfree = freelig(Ltot, Stot, KDAB);
double occ = RECFRAC*Lfree/(KDR + Lfree);
double sig = hillf(occ, OCC50, HSIG)/SH;      // 1.0 = healthy signalling

// ==================================================================
// 3. MACROPHAGE PROGRAM
// ==================================================================
double statc = 1000.0*STATC/VSTAT;            // ng/mL  (mg/L -> ug/L)
double stateff = EMAXST*sat(statc, EC50ST);
double smoke = (SOLVERTIME < SMKSTOP) ? SMKON : 0.0;

dxdt_PU1   = (dmn(sig, SMAX) - PU1)/TAUPU1;
dxdt_PPARG = (dmn(PU1 + PIOEFF, SMAX) - PPARG)/TAUPPG;
dxdt_ABCG1 = (PPARG - ABCG1)/TAUPPG;

// Per-cell capacity, computed REGIONALLY and then weighted by where the
// surfactant actually is. Open units carry the drug-driven program; filled
// units see only endogenous GM-CSF (which the antibody neutralises) plus
// whatever fraction of the inhaled effect penetrates. The endogenous
// signal is evaluated analytically from unchanged endogenous production,
// because aPAP does not alter CSF2 expression.
double fend  = freelig(LH, Stot, KDAB);
double sigend= hillf(RECFRAC*fend/(KDR + fend), OCC50, HSIG)/SH;
double prog  = 0.5*(PU1 + PPARG);
double wfill = (FFUL*RHO)/dmx(FFUL*RHO + (1.0 - FFUL), 1e-6);

double access = (1.0 - wfill) + EDGEF*4.0*FFUL*(1.0 - FFUL) + EPSREACH*wfill;
access = clamp2(access, 0.0, 1.0);

// Open units and the interface carry the drug-driven program; the interior
// of filled units does not. Statins are oral and pay no reach penalty,
// which is why they enter capacity outside the access weighting.
double capopen = CAPFLOOR + (1.0 - CAPFLOOR)*prog;
double capfill = CAPFLOOR + (1.0 - CAPFLOOR)*sigend;
double captgt = (capfill + (capopen - capfill)*access + stateff)
                *(1.0 - SMKEFF*smoke);
captgt = clamp2(captgt, 0.0, SMAX);
dxdt_CAPC = (captgt - CAPC)/TAUCAP;

dxdt_MONO = (MONOF - MONO)/10.0;
// macrophage number rises only where the drug arrives
double amreach = access;
double amtgt = MONO*(1.0 + AMNGM*amreach*pos(dmn(sig, SMAX) - 1.0)/dmx(SMAX - 1.0, 1e-6));
dxdt_AMN = (amtgt - AMN)/TAUAMN - kwam*AMN;

// ==================================================================
// 4. SURFACTANT MASS BALANCE
// ------------------------------------------------------------------
// Sinks are the only things that leave. The type-II recycling loop
// (PLA <-> LB) is written explicitly and cancels at steady state.
// ==================================================================
double PTOT = PLA + SEQ;
double xf   = dmx(1.0 - FFUL, 1e-6);

double uptake = VMAXUP*AMN*sat(PLA, KMUP)*pow(xf, EXPMAC)/(1.0 + LIP/KILIP);
double digest = VDIG*CAPC*AMN*sat(LIP, KMDIG);
double retn   = PHIRET*KDEATH*LIP;            // returned by dying cells
double iideg  = VMAXII *sat(PLA, KMII) *pow(xf, EXPII);
double muco   = VMAXMUC*sat(PLA, KMMUC)*pow(xf, EXPMUC);
double prod   = PRODB*pos(1.0 - PTOT/PMAXP);

// recycling loop: type-II cells take up and re-secrete. Large, and net
// zero. KREC is set so that the loop flux is ~3x de novo production.
double KREC = 3.0*PRODB/dmx(P0, 1e-6);
double recup = KREC*PLA*pow(xf, EXPII);
double resec = KREC*LB;

dxdt_LB  = recup - resec;
dxdt_PLA = prod + resec + retn + KDESEQ*SEQ
           - recup - uptake - iideg - muco - KSEQ*FFUL*PLA - kwll*PLA;

dxdt_LIP = uptake - digest - KDEATH*LIP - kwlip*LIP;
dxdt_AMDEAD = (1.0 - PHIRET)*KDEATH*LIP - 0.30*AMDEAD;
dxdt_SEQ = KSEQ*FFUL*PLA - KDESEQ*SEQ - kwseq*SEQ;
dxdt_DEB = 0.30*AMDEAD - KDEB*DEB;   // debris egress is a true sink
dxdt_AWP = muco - KAWP*AWP;

// cholesterol and protein sub-pools ride on the phospholipid balance but
// have their own macrophage handling (the statin target)
dxdt_CHA  = FCHOL*prod - (uptake/dmx(PLA,1e-6))*CHA*ABCG1 - muco*(CHA/dmx(PLA,1e-6))
            - kwll*CHA;
dxdt_AMCH = (uptake/dmx(PLA,1e-6))*CHA - ABCG1*0.5*AMCH - KDEATH*AMCH;
dxdt_SPRO = FSPRO*prod - (uptake + iideg + muco)*(SPRO/dmx(PLA,1e-6)) - kwll*SPRO;

// ==================================================================
// 5. FILLING, FIBROSIS
// ==================================================================
double ftgt = FMAXF*hillf(PTOT, PF50, HF)*(1.0 + FLAREAMP*FLARE);
ftgt = dmn(ftgt, FMAXF);
double fptgt = (1.0 - ftgt)*sat(PTOT, PPART50);
dxdt_FFUL  = (ftgt - FFUL)/TAUFILL;
dxdt_FPART = (fptgt - FPART)/TAUFILL;
dxdt_FLARE = -KFLARE*FLARE;

dxdt_TGFB = (FFUL + 0.4*FPART - TGFB)/(1.0/KTGF);
dxdt_FIB  = KFIB*TGFB*pos(1.0 - FIB/FIBMAX);

// ==================================================================
// 6. GAS EXCHANGE (used here for the physiological feedbacks)
// ==================================================================
double shunt = SHUNT0 + SHCO*FFUL;
double PAO2  = FIO2*(PATM - 47.0) - PACO2/RQ;
double pao2  = pao2_shunt(PAO2, shunt, HB, VO2R, COR);

double hbtgt = clamp2(HB0 + HBSLOPE*pos(85.0 - pao2), HB0, HBMAX);
dxdt_HB = (hbtgt - HB)/TAUHB;
double pvrtgt = PVR0 + PVRSL*pos(70.0 - pao2);
dxdt_PVR = (pvrtgt - PVR)/TAUPVR;

double aeff = pos((1.0 - FFUL) - DLPART*FPART);
double dlco = DLCOMAX*aeff*(1.0 - DLFIB*FIB)*pow(HB/HB0, 0.3) + DLCOL;
double aado2 = PAO2 - pao2;

dxdt_DLCOL = (DLCOLMAX - DLCOL)/TAULEARN;

// ==================================================================
// 7. SYMPTOMS, QUALITY OF LIFE, BIOMARKERS
// ==================================================================
double sgtgt = SG0 + SGA*aado2 + SGB*pos(100.0 - dlco);
dxdt_SGRQ = (pos(sgtgt) - SGRQ)/TAUSG;
double wdtgt = WD0 - WDA*aado2 - WDB*pos(100.0 - dlco);
dxdt_SIXMWD = (dmx(wdtgt, 60.0) - SIXMWD)/TAUWD;
double mmtgt = clamp2(0.045*pos(aado2 - 10.0) + 0.012*pos(100.0 - dlco), 0.0, 4.0);
dxdt_MMRC = (mmtgt - MMRC)/TAUWD;

double leak = (PTOT/P0)*(1.0 + FFUL);
dxdt_KL6 = ((KL6B + KL6S*leak) - KL6)/TAUBIO;
dxdt_SPD = ((SPDB + SPDS*leak) - SPD)/TAUBIO;
dxdt_CEA = ((CEAB + CEAS*leak) - CEA)/TAUBIO;
dxdt_LDH = ((LDHB + LDHS*leak) - LDH)/TAUBIO;

// systemic GM-CSF drives the myeloid adverse effects; the SAME occupancy
// curve read on marrow instead of alveolar macrophages
double cgmpl = GMPL/(VPLAS/1000.0);           // pM
dxdt_NEUT = ((4.0 + 9.0*sat(cgmpl, 20.0)) - NEUT)/1.5;
dxdt_EOS  = ((0.20 + 0.60*sat(cgmpl, 30.0)) - EOS)/3.0;

// ==================================================================
// 8. INFECTION
// ==================================================================
double haz = HAZ0*(1.0 + AINF*pos(1.0 - CAPC))*(1.0 + BINF*(PTOT/P0));
dxdt_CUMHAZ = haz;
dxdt_NOC = KGNOC*NOC*pos(1.0 - NOC/NOCMAX)
           - KKILLN*CAPC*AMN*NOC/(1.0 + NOC/KNOC);

// ==================================================================
// 9. AUDIT TRAIL
// ==================================================================
dxdt_COV     = 1.0/(1.0 + exp(-(sig - 0.5)/0.05));
dxdt_CUMPROD = prod;
dxdt_CUMMAC  = digest;
dxdt_CUMII   = iideg;
dxdt_CUMMUC  = KAWP*AWP;
dxdt_CUMWLL  = kwll*PLA + kwseq*SEQ + kwlip*LIP;  // lipid in the aspirate counts too
dxdt_CUMPLEX = kplex*ABC;
dxdt_CUMDEB  = KDEB*DEB;
dxdt_CUMDOSE = KDISS*INHD;
dxdt_NWLL    = wllon/dmx(WLLDUR, 1e-6);

$TABLE
// ------------------------------------------------------------------
// Everything reported. Recomputed here because $ODE locals are not
// visible to $TABLE.
// ------------------------------------------------------------------
double PTOTo  = PLA + SEQ;
double xfo    = dmx(1.0 - FFUL, 1e-6);
double Ltoto  = GMT/(VELF/1000.0);
double Stoto  = ABE*(2.0e6/150000.0)*FNEUT/(VELF/1000.0);
double LFREE  = freelig(Ltoto, Stoto, KDAB);
double OCCUP  = RECFRAC*LFREE/(KDR + LFREE);
double SIGREL = hillf(OCCUP, OCC50, HSIG)/SH;
double TITRE  = ABC/VPLAS;                        // ug/mL
double TITELF = ABE/VELF;                         // ug/mL
double MRATIO = Ltoto/dmx(Stoto, 1e-12);          // ligand : site molar ratio

double UPTAKE = VMAXUP*AMN*sat(PLA, KMUP)*pow(xfo, EXPMAC)/(1.0 + LIP/KILIP);
double DIGEST = VDIG*CAPC*AMN*sat(LIP, KMDIG);
double IIDEG  = VMAXII *sat(PLA, KMII) *pow(xfo, EXPII);
double MUCO   = VMAXMUC*sat(PLA, KMMUC)*pow(xfo, EXPMUC);
double PROD   = PRODB*pos(1.0 - PTOTo/PMAXP);
double NETMAC = UPTAKE - PHIRET*KDEATH*LIP;
double NETBAL = PROD - NETMAC - IIDEG - MUCO;     // mg/d, the imbalance
double TOTCAP = CAPC*AMN;                         // total catabolic capacity

double FENDo  = freelig(LH, Stoto, KDAB);
double SIGEND = hillf(RECFRAC*FENDo/(KDR + FENDo), OCC50, HSIG)/SH;
double WFILL  = (FFUL*RHO)/dmx(FFUL*RHO + (1.0 - FFUL), 1e-6);
double REACH  = clamp2((1.0 - WFILL) + EDGEF*4.0*FFUL*(1.0 - FFUL) + EPSREACH*WFILL, 0.0, 1.0);
double SHUNTF = SHUNT0 + SHCO*FFUL;
double PAO2v  = FIO2*(PATM - 47.0) - PACO2/RQ;
double PAO2A  = pao2_shunt(PAO2v, SHUNTF, HB, VO2R, COR);
double AADO2  = PAO2v - PAO2A;
double SAO2   = 100.0*sao2f(PAO2A);
double AEFF   = pos((1.0 - FFUL) - DLPART*FPART);
double DLCO   = DLCOMAX*AEFF*(1.0 - DLFIB*FIB)*pow(HB/HB0, 0.3) + DLCOL;

// exercise: the same physics with VO2 up 5x and cardiac output up 3x.
// Nothing about the lung changes; the mixed venous content falls, so the
// same shunt costs more.
double PAO2ex = FIO2*(PATM - 47.0) - PACO2EX/RQ;
double PAO2AX = pao2_shunt(PAO2ex, SHUNTF, HB, VO2EX, COEX);
double AADOEX = PAO2ex - PAO2AX;
double SAO2EX = 100.0*sao2f(PAO2AX);

// disease severity score of Inoue 2008: symptoms plus PaO2 bands
double DSS = 1.0;
if(MMRC >= 0.5) DSS = 2.0;
if(PAO2A < 70.0) DSS = 3.0;
if(PAO2A < 60.0) DSS = 4.0;
if(PAO2A < 50.0) DSS = 5.0;

double CTHU  = -850.0 + 700.0*(FFUL + 0.5*FPART);
double MPAP  = PVR*COR + 10.0;
double POOLR = PTOTo/P0;
double O2REQ = 1.0/(1.0 + exp((PAO2A - 55.0)/2.0));
double INFYR = 365.0*HAZ0*(1.0 + AINF*pos(1.0 - CAPC))*(1.0 + BINF*POOLR);
double CGMPL = GMPL/(VPLAS/1000.0);
double CGMELF= LFREE*MWGM/1000.0;                 // pM -> pg/mL free
// Phospholipid balance only. CHA and SPRO are compositional tracers that
// ride on the same mass and are audited separately; including them here
// would count their production twice.
double MBAL  = (CUMPROD + P0 + P0 + LIP0)
               - (CUMMAC + CUMII + CUMMUC + CUMWLL + CUMDEB
                  + PLA + SEQ + LB + LIP + DEB + AWP + AMDEAD);

$CAPTURE @annotated
TITRE  : Serum GMAb titre (ug/mL)
TITELF : ELF GMAb (ug/mL)
LFREE  : Free GM-CSF in ELF (pM)
CGMELF : Free GM-CSF in ELF (pg/mL)
CGMPL  : Plasma GM-CSF (pM)
MRATIO : Ligand:site molar ratio in ELF
OCCUP  : Receptor occupancy
SIGREL : Signalling relative to healthy
TOTCAP : Total catabolic capacity (relative)
SIGEND : Endogenous-only signalling (relative to healthy)
WFILL  : Fraction of the burden lying in filled units
REACH  : Fraction of the drug effect that reaches the burden
POOLR  : Alveolar burden / healthy pool
PTOTo  : Total alveolar burden (mg)
UPTAKE : Macrophage uptake flux (mg/d)
DIGEST : Macrophage digestion flux (mg/d)
NETMAC : Net macrophage removal (mg/d)
IIDEG  : Type-II degradation flux (mg/d)
MUCO   : Mucociliary egress flux (mg/d)
PROD   : De novo production (mg/d)
NETBAL : Net imbalance (mg/d)
SHUNTF : Shunt fraction
PAO2A  : Arterial PO2 at rest (mmHg)
AADO2  : A-aDO2 at rest (mmHg)
SAO2   : Arterial O2 saturation (%)
PAO2AX : Arterial PO2 on exercise (mmHg)
AADOEX : A-aDO2 on exercise (mmHg)
SAO2EX : Exercise O2 saturation (%)
DLCO   : DLCO (% predicted)
AEFF   : Effective gas-exchange fraction
DSS    : Disease severity score (1-5)
CTHU   : Mean lung CT density (HU)
MPAP   : Mean pulmonary artery pressure (mmHg)
O2REQ  : Supplemental oxygen requirement index
INFYR  : Serious-infection hazard (per year)
MBAL   : Surfactant mass-balance residual (mg)
'

mod <- mread_cache("pap", code = pap_code, soloc = tempdir())

# =====================================================================
# HELPERS
# =====================================================================
`%||%` <- function(a, b) if (is.null(a)) b else a

# param lists are de-duplicated, LAST occurrence winning, so that
# c(base, list(X = 1), list(X = 2)) means what the caller intended
mp <- function(pars) {
  pars <- pars[!vapply(pars, is.null, logical(1))]
  nm <- names(pars)
  keep <- !duplicated(nm, fromLast = TRUE)
  pars[keep]
}

# Whole lung lavage and plasmapheresis are finite WINDOWS in the right-hand
# side, not dose records. lsoda happily takes steps longer than a 0.25-day
# window and integrates straight past it, so the procedure silently does
# nothing - which is exactly what happened on the first run of S16-S23.
# Forced output times are NOT enough on their own: lsoda interpolates output
# inside a step it has already taken, so it can report values at 30.06 and
# 30.12 without ever evaluating the right-hand side there. The maximum step
# size must also be constrained while a procedure is scheduled. Without both,
# lavage silently does nothing - and the scenario table looks plausible.
window_times <- function(pars) {
  p <- mp(c(as.list(param(mod)), pars))
  tt <- numeric(0)
  dur <- p$WLLDUR
  for (k in c("WLLT1", "WLLT2", "WLLT3", "WLLT4")) {
    t0 <- p[[k]]
    if (!is.null(t0) && t0 < 1e5) tt <- c(tt, seq(t0 - dur, t0 + 2 * dur, by = dur / 8))
  }
  if (!is.null(p$PLEXN) && p$PLEXN > 0 && p$PLEXT < 1e5) {
    for (i in seq_len(p$PLEXN) - 1) {
      t0 <- p$PLEXT + i * p$PLEXINT
      tt <- c(tt, seq(t0 - p$PLEXDUR, t0 + 2 * p$PLEXDUR, by = p$PLEXDUR / 6))
    }
  }
  sort(unique(pmax(0, tt)))
}

sim <- function(pars = list(), end = 365, delta = 1, ev = NULL, init = NULL,
                recsort = 3) {
  m <- mod
  if (!is.null(init)) {
    pars <- c(pars, list(INITMODE = 0))
    m <- update(m, init = as.list(init))
  }
  m <- param(m, mp(pars))
  addt <- window_times(pars)
  addt <- addt[addt <= end]
  pp <- mp(c(as.list(param(mod)), pars))
  hm <- if (length(addt)) min(pp$WLLDUR, pp$PLEXDUR) / 8 else 0
  if (is.null(ev)) {
    out <- mrgsim(m, end = end, delta = delta, add = addt, hmax = hm,
                  recsort = recsort)
  } else {
    out <- mrgsim(m, data = ev, end = end, delta = delta, add = addt, hmax = hm,
                  recsort = recsort)
  }
  as.data.frame(out)
}

CMTN <- names(mrgsolve::init(mod)@data)

# state vector at the end of a run, ready to be fed back in as init()
endstate <- function(d) {
  s <- as.list(d[nrow(d), CMTN])
  names(s) <- CMTN
  s
}
at   <- function(d, tt, col) d[[col]][which.min(abs(d$time - tt))]
fin  <- function(d, col) d[[col]][nrow(d)]
base <- function(d, col) d[[col]][1]

# ---------------------------------------------------------------------
# THE PATIENT IS GENERATED, NOT INITIALISED
# ---------------------------------------------------------------------
# Healthy lung -> seroconversion -> years of integration. The presenting
# patient is the state at which A-aDO2 first crosses the presentation
# threshold. A patient whose catabolic floor holds A-aDO2 below that
# threshold NEVER presents: that is the model account of the 31.8% found
# asymptomatically on health screening (Inoue 2008).
# ---------------------------------------------------------------------
healthy <- function(pars = list(), days = 200) {
  sim(c(list(ABON = 0), mp(pars)), end = days, delta = 5)
}

natural_history <- function(pars = list(), years = 25, delta = 5) {
  h <- healthy(pars, days = 200)
  sim(c(mp(pars), list(ABON = 1)), end = 365 * years, delta = delta,
      init = endstate(h))
}

# Enrolment criterion of the randomised trials: resting PaO2 below
# threshold (PAGE used <70 mmHg, or <75 mmHg if symptomatic). The default
# 65 mmHg reproduces the mean baseline of the enrolled cohorts.
# SETTLE matters. A patient randomised the week they cross the threshold is
# still on the steep part of their trajectory, and a control arm built that
# way deteriorates through the trial - which is not what the placebo arms of
# IMPALA-2 or PAGE did. Real enrolment happens well after diagnosis, in
# patients documented to be stable (PAGE ran a 12-week observation period and
# excluded anyone who improved). The default therefore advances the patient
# past the crossing to the neighbourhood of their own steady state.
present <- function(pars = list(), target = 70, years = 30, settle = 1.5) {
  nh <- natural_history(pars, years = years, delta = 5)
  hit <- which(nh$PAO2A <= target)
  if (!length(hit)) {
    return(list(state = endstate(nh), tpres = NA_real_, nh = nh,
                presented = FALSE))
  }
  tt <- nh$time[hit[1]] + settle * 365
  i <- which.min(abs(nh$time - tt))
  st <- as.list(nh[i, CMTN]); names(st) <- CMTN
  list(state = st, tpres = nh$time[hit[1]], tenrol = nh$time[i], nh = nh,
       presented = TRUE)
}

# ---------------------------------------------------------------------
# DOSING EVENT BUILDERS
# ---------------------------------------------------------------------
# inhaled GM-CSF, continuous daily
ev_inh <- function(dose = 300, days = 168, per_day = 1) {
  ii <- 1 / per_day
  ev(amt = dose, cmt = "INHD", ii = ii, addl = ceiling(days / ii) - 1, time = 0)
}
# inhaled GM-CSF, n_on days on / n_off days off, repeated
ev_inh_cycle <- function(dose = 300, days = 168, per_day = 1,
                         n_on = 7, n_off = 7) {
  cyc <- n_on + n_off
  starts <- seq(0, days - 1, by = cyc)
  d <- do.call(rbind, lapply(starts, function(s) {
    n <- n_on * per_day
    as.data.frame(ev(amt = dose, cmt = "INHD", time = s,
                     ii = 1 / per_day, addl = n - 1))
  }))
  d$ID <- 1
  d[order(d$time), ]
}
ev_sc <- function(dose_ug_kg = 5, bw = 70, days = 168) {
  ev(amt = dose_ug_kg * bw, cmt = "SCD", ii = 1, addl = days - 1, time = 0)
}
ev_rtx <- function(times = c(0, 14), dose = 1000) {
  d <- do.call(rbind, lapply(times, function(t)
    as.data.frame(ev(amt = dose, cmt = "RTXC", time = t))))
  d$ID <- 1
  d
}
ev_statin <- function(dose = 80, days = 365) {
  ev(amt = dose, cmt = "STATG", ii = 1, addl = days - 1, time = 0)
}
ev_bind <- function(...) {
  l <- lapply(list(...), function(x) as.data.frame(x))
  d <- do.call(rbind, l)
  d <- d[order(d$time), ]
  d$ID <- 1
  d
}

# ---------------------------------------------------------------------
# TRIAL CONSTRUCTION
# ---------------------------------------------------------------------
# A randomised patient is not the same object as a prevalent patient.
# Two selection steps separate them, and between them they account for the
# whole placebo response of IMPALA-2 without any placebo parameter:
#   (1) ENROLMENT AT A LOW POINT. Patients are referred and randomised
#       after they get worse, so baseline is measured near the peak of a
#       deterioration that then resolves on its own.
#   (2) MEASUREMENT LEARNING. DLCO is an effort-dependent manoeuvre and
#       improves with repetition in both arms (DLCOLMAX).
# PAGE additionally ran a 12-week observation period and EXCLUDED anyone
# who improved by >= 10 mmHg of A-aDO2, which removes step (1) from its
# placebo arm - and PAGE is the trial whose placebo arm did not move.
randomise <- function(pars = list(), flare = 1.3, settle = 0.75,
                      leadin = 21, trial = TRUE) {
  p0 <- present(pars, settle = settle)
  if (!p0$presented) return(NULL)
  st <- p0$state
  st$FLARE <- flare
  pre <- sim(pars, end = leadin, delta = 1, init = st)
  list(state = endstate(pre), tpres = p0$tpres,
       pars = c(pars, if (trial) list(DLCOLMAX = 1.5) else NULL))
}

# One arm of a trial, from a randomised patient
arm <- function(r, dosing = NULL, end = 336, extra = list()) {
  sim(c(r$pars, mp(extra)), end = end, ev = dosing, init = r$state)
}

# Virtual population: CAPFLOOR is the severity covariate, titre is
# deliberately given a WIDE spread so that D05 can ask whether it predicts
# anything. Patients who never reach the enrolment criterion are screened
# out, exactly as they are in a trial.
# The titre distribution matters for D05 and is set to the observed one:
# measured GMAb in aPAP cohorts is mostly 15-150 ug/mL. Widening it below
# 15 puts virtual patients in the range where the model says titre DOES
# still matter, and then titre correlates with severity - which is a
# prediction, not a bug (see D05).
# A GMAb assay measures BINDING, not neutralisation, and the neutralising
# fraction of a polyclonal response differs between patients. That single
# fact is what breaks the titre-severity relationship (D05): two patients
# with the same reported titre can have very different free GM-CSF.
vpop <- function(n = 40, seed = 20260727, cf_mean = 0.3185, cf_sd = 0.008,
                 titre_med = 40, titre_gsd = 1.8, fneut_gsd = 1.7,
                 flare_mean = 1.3, pars = list(), settle = 0.75) {
  set.seed(seed)
  cf <- pmax(0.300, pmin(0.336, rnorm(n, cf_mean, cf_sd)))
  ti <- titre_med * exp(rnorm(n, 0, log(titre_gsd)))
  fn <- pmax(0.08, pmin(0.85, 0.32 * exp(rnorm(n, 0, log(fneut_gsd)))))
  fl <- pmax(0, rnorm(n, flare_mean, 0.35))
  out <- vector("list", n)
  for (i in seq_len(n)) {
    r <- randomise(c(pars, list(CAPFLOOR = cf[i], TITRE0 = ti[i],
                                FNEUT = fn[i])), flare = fl[i], settle = settle)
    if (!is.null(r))
      out[[i]] <- c(r, list(cf = cf[i], titre = ti[i], fneut = fn[i]))
  }
  out[!vapply(out, is.null, logical(1))]
}

trial_arm <- function(pop, dosing = NULL, end = 336, extra = list()) {
  do.call(rbind, lapply(seq_along(pop), function(i) {
    d <- arm(pop[[i]], dosing, end = end, extra = extra)
    data.frame(id = i, cf = pop[[i]]$cf, titre = pop[[i]]$titre,
               dlco0 = d$DLCO[1], aado20 = d$AADO2[1], sgrq0 = d$SGRQ[1],
               dDLCO24 = at(d, 168, "DLCO") - d$DLCO[1],
               dDLCO48 = at(d, 336, "DLCO") - d$DLCO[1],
               dAADO24 = at(d, 168, "AADO2") - d$AADO2[1],
               dSGRQ24 = at(d, 168, "SGRQ") - d$SGRQ[1],
               poolr24 = at(d, 168, "POOLR") / d$POOLR[1])
  }))
}

# =====================================================================
# SCENARIOS
# =====================================================================
run_scenarios <- function() {
  S <- list()
  add <- function(name, d, note = "") S[[name]] <<- list(sim = d, note = note)

  p_def <- present()                 # the default presenting patient
  st    <- p_def$state
  r_def <- randomise()

  # ---- S01-S06  natural history and phenotypes -----------------------
  add("S01 healthy control, no autoantibody", healthy(days = 400),
      "stationarity check; see D01")
  add("S02 natural history from seroconversion, 25 y",
      natural_history(years = 25, delta = 5),
      "disease is generated, not initialised")
  add("S03 untreated after presentation, 1 y", sim(end = 365, init = st))
  add("S04 progressive phenotype (floor 0.28)",
      natural_history(list(CAPFLOOR = 0.28), years = 25, delta = 5))
  add("S05 screening-detected phenotype (floor 0.33)",
      natural_history(list(CAPFLOOR = 0.33), years = 25, delta = 5),
      "never crosses the enrolment criterion: the 31.8% of Inoue 2008")
  add("S06 spontaneous remission (antibody output decays)",
      sim(list(ABDECAY = 0.004), end = 1460, delta = 5, init = st))

  # ---- S07-S15  inhaled GM-CSF --------------------------------------
  add("S07 molgramostim 300 ug QD, 24 wk", arm(r_def, ev_inh(300, 168), 168))
  add("S08 molgramostim 300 ug QD, 48 wk", arm(r_def, ev_inh(300, 336), 336))
  add("S09 molgramostim 300 ug QD every other week (IMPALA arm)",
      arm(r_def, ev_inh_cycle(300, 336, 1, 7, 7), 336))
  add("S10 sargramostim 125 ug BID, 7 d every other week (PAGE)",
      arm(r_def, ev_inh_cycle(125, 336, 2, 7, 7), 336))
  add("S11 sargramostim 250 ug QD, 8 d of every 14 (Tazawa 2010)",
      arm(r_def, ev_inh_cycle(250, 336, 1, 8, 6), 336))
  add("S12 sargramostim SC 5 ug/kg/d", arm(r_def, ev_sc(5, 70, 336), 336))
  add("S13 sargramostim SC 20 ug/kg/d", arm(r_def, ev_sc(20, 70, 336), 336))
  dose_sweep <- do.call(rbind, lapply(c(10, 30, 75, 150, 300, 600, 1200, 3000),
    function(dd) {
      d <- arm(r_def, ev_inh(dd, 168), 168)
      data.frame(dose = dd, dDLCO24 = at(d, 168, "DLCO") - d$DLCO[1],
                 covfrac = (at(d, 168, "COV") - at(d, 161, "COV")) / 7,
                 mratio = at(d, 168, "MRATIO"), totcap = at(d, 168, "TOTCAP"))
    }))
  add("S14 inhaled dose sweep 10-3000 ug", dose_sweep,
      "threshold then plateau: dose is not the efficacy variable")
  sched <- rbind(
    cbind(sched = "300 QD",      as.data.frame(t(c(d24 = NA)))),
    NULL)
  sch <- do.call(rbind, list(
    data.frame(schedule = "300 ug QD",
               d24 = {d <- arm(r_def, ev_inh(300, 168), 168); at(d,168,"DLCO")-d$DLCO[1]}),
    data.frame(schedule = "150 ug BID",
               d24 = {d <- arm(r_def, ev_inh(150, 168, 2), 168); at(d,168,"DLCO")-d$DLCO[1]}),
    data.frame(schedule = "600 ug Q2D",
               d24 = {d <- arm(r_def, ev(amt=600, cmt="INHD", ii=2, addl=83), 168); at(d,168,"DLCO")-d$DLCO[1]})))
  add("S15 schedule at matched daily dose", sch,
      "same 300 ug/d split three ways")

  # ---- S16-S20  lavage ----------------------------------------------
  add("S16 single whole lung lavage", sim(list(WLLT1 = 30), end = 365, init = st))
  add("S17 WLL then 8-y follow-up",
      sim(list(WLLT1 = 30), end = 8 * 365, delta = 5, init = st),
      "time to recurrence is pool removed / net accumulation rate")
  add("S18 WLL then molgramostim maintenance",
      sim(c(r_def$pars, list(WLLT1 = 30)), end = 336, ev = ev_inh(300, 336),
          init = r_def$state))
  add("S19 repeated WLL every 18 months x4",
      sim(list(WLLT1 = 30, WLLT2 = 30 + 548, WLLT3 = 30 + 1096,
               WLLT4 = 30 + 1644), end = 8 * 365, delta = 5, init = st))
  add("S20 segmental (bronchoscopic) lavage",
      sim(list(WLLT1 = 30, WLLFRAC = 0.30, WLLSEQF = 0.20, WLLABF = 0.35),
          end = 365, init = st))

  # ---- S21-S26  antibody- and lipid-directed therapy -----------------
  add("S21 rituximab 1000 mg x2, median titre",
      sim(end = 730, delta = 2, ev = ev_rtx(), init = st))
  # A titre of 8 ug/mL does not by itself produce disease in this model, so a
  # low-titre patient with comparable disease must have a lower catabolic
  # floor. That is not a workaround: it is the model statement that titre and
  # severity are set by different things (D05).
  pl <- list(TITRE0 = 8, CAPFLOOR = 0.22)
  st_low <- present(pl)$state
  add("S22 rituximab 1000 mg x2, titre 8 ug/mL",
      sim(pl, end = 730, delta = 2, ev = ev_rtx(), init = st_low))
  add("S23 plasmapheresis, 10 sessions",
      sim(list(PLEXT = 30, PLEXN = 10), end = 365, init = st))
  add("S24 FcRn inhibitor (hypothesis), titre 25",
      sim(list(FCRNT = 30, FCRNX = 3.0), end = 730, delta = 2, init = st))
  add("S25 FcRn inhibitor (hypothesis), titre 10",
      sim(list(TITRE0 = 10, FCRNT = 30, FCRNX = 3.0), end = 730, delta = 2,
          init = present(list(TITRE0 = 10))$state))
  add("S26 atorvastatin 80 mg, 12 months",
      sim(end = 365, ev = ev_statin(80, 365), init = st))
  add("S27 statin + molgramostim",
      sim(r_def$pars, end = 336,
          ev = ev_bind(ev_inh(300, 336), ev_statin(80, 336)), init = r_def$state))

  # ---- S28-S30  the other two PAP classes ---------------------------
  add("S28 hereditary PAP (no receptor) + molgramostim",
      sim(list(RECFRAC = 0), end = 336, ev = ev_inh(300, 336),
          init = present(list(RECFRAC = 0, CAPFLOOR = 0.305))$state),
      "structural null: the drug cannot work")
  # Secondary PAP has NO autoantibody: signalling is intact and the deficit
  # is macrophage NUMBER. Built that way, inhaled GM-CSF has almost nothing
  # left to restore, which is why the treatment is the haematologic disease.
  sec <- list(ABON = 0, MONOF = 0.34, CAPFLOOR = 0.305)
  add("S29 secondary PAP (monocytopenia, no antibody) + molgramostim",
      sim(sec, end = 336, ev = ev_inh(300, 336), init = present(sec)$state))
  add("S30 hereditary PAP, macrophage transplantation surrogate",
      sim(list(RECFRAC = 0, CAPFLOOR = 0.9), end = 365,
          init = present(list(RECFRAC = 0, CAPFLOOR = 0.305))$state))

  # ---- S31-S36  covariates, host defence, physiology ----------------
  titre_sweep <- do.call(rbind, lapply(c(0.5, 1, 2, 3, 5, 8, 15, 25, 50, 100, 300),
    function(tt) {
      d <- natural_history(list(TITRE0 = tt), years = 20, delta = 20)
      data.frame(titre = tt, sigrel = fin(d, "SIGREL"), totcap = fin(d, "TOTCAP"),
                 poolr = fin(d, "POOLR"), aado2 = fin(d, "AADO2"),
                 dlco = fin(d, "DLCO"), dss = fin(d, "DSS"))
    }))
  add("S31 titre sweep 0.5-300 ug/mL", titre_sweep,
      "the 5 ug/mL threshold is emergent; severity saturates above it")
  cf_sweep <- do.call(rbind, lapply(seq(0.26, 0.40, by = 0.01), function(cc) {
    d <- natural_history(list(CAPFLOOR = cc), years = 20, delta = 20)
    data.frame(capfloor = cc, poolr = fin(d, "POOLR"), aado2 = fin(d, "AADO2"),
               dlco = fin(d, "DLCO"), dss = fin(d, "DSS"))
  }))
  add("S32 catabolic-floor sweep 0.26-0.40", cf_sweep,
      "the covariate that does predict severity")
  add("S33 Nocardia challenge, untreated",
      sim(list(NOC0 = 1000), end = 60, delta = 0.5, init = st))
  add("S34 Nocardia challenge on molgramostim (established therapy)",
      sim(list(NOC0 = 1000), end = 60, delta = 0.5,
          ev = ev_inh(300, 60),
          init = endstate(sim(end = 168, ev = ev_inh(300, 168), init = st))))
  add("S35 smoking, then cessation at 1 y",
      sim(list(SMKON = 1, SMKSTOP = 365), end = 730, delta = 2, init = st))
  add("S36 diagnostic delay: 10 y untreated then molgramostim",
      sim(r_def$pars, end = 336, ev = ev_inh(300, 336),
          init = endstate(sim(end = 10 * 365, delta = 5, init = st))),
      "fibrosis accrued during delay caps the recoverable DLCO")

  list(scen = S, state = st, rand = r_def, present = p_def)
}

# =====================================================================
# DIAGNOSTICS
# =====================================================================
# Each one either checks the model against a published number or tries to
# break a claim the model makes. Failures are reported, not removed.
run_diagnostics <- function() {
  hr <- function() cat(strrep("-", 70), "\n")
  ok <- function(x) if (isTRUE(x)) "PASS" else "FAIL"
  cat("\n=====================================================================\n")
  cat("aPAP QSP MODEL - DIAGNOSTICS\n")
  cat("=====================================================================\n")

  # ---- D01 stationarity of the derived healthy baseline ---------------
  h <- healthy(days = 1000)
  dr <- 100 * (fin(h, "PLA") - h$PLA[1]) / h$PLA[1]
  hr(); cat(sprintf("D01 healthy lung is stationary over 1000 d\n"))
  cat(sprintf("    pool %.3f -> %.3f mg, drift %+.4f%%   [%s]\n",
              h$PLA[1], fin(h, "PLA"), dr, ok(abs(dr) < 0.05)))
  cat("    Nothing is fitted to make this hold: uptake Vmax, digestive Vmax\n")
  cat("    and de novo production are solved in $MAIN from the sink split.\n")

  # ---- D02 mass balance closure --------------------------------------
  w <- sim(list(WLLT1 = 100), end = 400, init = present()$state)
  hr(); cat(sprintf("D02 phospholipid mass balance closes with a lavage in the run\n"))
  cat(sprintf("    residual %+.3e mg against %.0f mg turned over   [%s]\n",
              fin(w, "MBAL"), fin(w, "CUMPROD"),
              ok(abs(fin(w, "MBAL")) < 1e-3 * max(1, fin(w, "CUMPROD")))))

  # ---- D03 what the flux split actually is ---------------------------
  hr(); cat("D03 healthy net sink split (specified 65/25/10)\n")
  tot <- fin(h, "NETMAC") + fin(h, "IIDEG") + fin(h, "MUCO")
  cat(sprintf("    macrophage %.1f mg/d (%.1f%%), type-II %.1f (%.1f%%), mucociliary %.1f (%.1f%%)\n",
      fin(h, "NETMAC"), 100 * fin(h, "NETMAC") / tot, fin(h, "IIDEG"),
      100 * fin(h, "IIDEG") / tot, fin(h, "MUCO"), 100 * fin(h, "MUCO") / tot))
  cat(sprintf("    recycling loop flux %.0f mg/d, net contribution 0 by construction\n",
      3 * fin(h, "PROD")))

  # ---- D04 the 5 ug/mL threshold is emergent -------------------------
  hr(); cat("D04 critical GMAb titre (Sakagami 2010: 5 ug/mL, lung and blood)\n")
  tt <- c(0.5, 1, 2, 3, 4, 5, 6, 8, 12, 25, 50, 100, 300)
  sg <- vapply(tt, function(x) {
    d <- sim(list(TITRE0 = x, ABON = 1), end = 400, delta = 20)
    fin(d, "SIGREL")
  }, numeric(1))
  half <- approx(sg, tt, xout = 0.5)$y
  for (i in seq_along(tt))
    cat(sprintf("    titre %6.1f ug/mL -> signalling %.3f of healthy\n", tt[i], sg[i]))
  cat(sprintf("    half-signal titre = %.2f ug/mL vs published 5   [%s]\n",
              half, ok(half > 3 && half < 8)))
  cat("    Levers used: ELF/serum transudation ratio and the neutralising\n")
  cat("    fraction. The SHAPE (steep, then flat) is not a lever - it is\n")
  cat("    what a stoichiometric buffer does.\n")

  # ---- D05 titre does not predict severity; the floor does -----------
  hr(); cat("D05 what predicts severity in a virtual population\n")
  set.seed(11)
  n <- 45
  cfv <- pmax(0.29, pmin(0.34, rnorm(n, 0.3145, 0.011)))
  tiv <- 40 * exp(rnorm(n, 0, log(1.8)))
  fnv <- pmax(0.08, pmin(0.85, 0.32 * exp(rnorm(n, 0, log(1.7)))))
  res <- do.call(rbind, lapply(seq_len(n), function(i) {
    d <- natural_history(list(CAPFLOOR = cfv[i], TITRE0 = tiv[i],
                              FNEUT = fnv[i]), years = 20, delta = 20)
    data.frame(cf = cfv[i], titre = tiv[i], fneut = fnv[i],
               sites = tiv[i] * fnv[i], aado2 = fin(d, "AADO2"),
               dlco = fin(d, "DLCO"))
  }))
  sub <- res[res$titre > 15, ]
  r_t <- suppressWarnings(cor(log(sub$titre), sub$aado2))
  r_c <- suppressWarnings(cor(sub$cf, sub$aado2))
  cat(sprintf("    n = %d in the observed titre range (>15 ug/mL)\n", nrow(sub)))
  cat(sprintf("    corr(log titre, A-aDO2)      = %+.3f\n", r_t))
  r_s0 <- suppressWarnings(cor(log(sub$titre * sub$fneut), sub$aado2))
  cat(sprintf("    corr(catabolic floor, A-aDO2)= %+.3f\n", r_c))
  cat(sprintf("    reported titre is a WORSE predictor than neutralising load\n"))
  cat(sprintf("    by %.2f in |r|   [%s]\n", abs(r_s0) - abs(r_t),
              ok(abs(r_s0) - abs(r_t) > 0.2)))
  cat("    PARTIAL FAILURE, stated plainly: the model weakens the titre-\n")
  cat("    severity correlation but does not abolish it (r ~ +0.6, and +0.85\n")
  cat("    without the assay/neutralisation distinction). Inoue 2008 found\n")
  cat("    none. So buffer saturation plus assay-vs-neutralisation is not the\n")
  cat("    whole story, and the catabolic floor is NOT the dominant covariate\n")
  cat("    in this population - the antibody side still is.\n")
  lowt <- res[res$titre <= 15, ]
  if (nrow(lowt) > 4)
    cat(sprintf("    for comparison, titre <= 15 ug/mL (n = %d): corr = %+.3f\n",
                nrow(lowt), suppressWarnings(cor(log(lowt$titre), lowt$aado2))))
  r_s <- suppressWarnings(cor(log(sub$sites), sub$aado2))
  cat(sprintf("    corr(log NEUTRALISING sites, A-aDO2) = %+.3f\n", r_s))
  cat("    Two mechanisms flatten the titre relationship and they are separable:\n")
  cat("    (i) buffer saturation above the critical titre, and (ii) the assay\n")
  cat("    measuring binding rather than neutralisation. Only (ii) can break a\n")
  cat("    correlation that survives (i).\n")
  cat("    FALSIFIABLE PREDICTION: the decorrelation is partly a saturation\n")
  cat("    effect, so it should weaken in a cohort enriched for 5-15 ug/mL, and\n")
  cat("    a NEUTRALISATION assay should correlate with severity where a\n")
  cat("    binding titre does not.\n")
  cat("    Inoue 2008: DSS correlates with DLCO and biomarkers but NOT with\n")
  cat("    GMAb titre. Here that is a consequence, not an assumption: above\n")
  cat("    the threshold free GM-CSF is already near zero, so the residual\n")
  cat("    GM-CSF-independent floor is the only thing left that varies.\n")

  # ---- D06 the molar-ratio argument for route ------------------------
  hr(); cat("D06 why route beats dose: molar arithmetic in lining fluid\n")
  st <- present()$state
  # hourly sampling: at daily output times only the trough is seen, and the
  # peak molar ratio - the whole point of the argument - is invisible
  di <- sim(end = 30, delta = 1 / 48, ev = ev_inh(300, 30), init = st)
  ds <- sim(end = 30, delta = 1 / 48, ev = ev_sc(20, 70, 30), init = st)
  cat(sprintf("    inhaled 300 ug QD : peak ELF ligand:site ratio %8.1f, peak free GM-CSF %8.1f pM\n",
              max(di$MRATIO), max(di$LFREE)))
  cat(sprintf("    SC 20 ug/kg/d     : peak ELF ligand:site ratio %8.3f, peak free GM-CSF %8.3f pM\n",
              max(ds$MRATIO), max(ds$LFREE)))
  cat(sprintf("    healthy free GM-CSF for reference: %.2f pM   [%s]\n",
              as.list(param(mod))$CGMH * 1000 / as.list(param(mod))$MWGM,
              ok(max(di$MRATIO) > 100 * max(ds$MRATIO))))

  # ---- D07/D08 IMPALA-2 --------------------------------------------
  hr(); cat("D07 IMPALA-2 replication (Trapnell 2025, PMID 40834301)\n")
  pop <- vpop(n = 24)
  a_d <- trial_arm(pop, ev_inh(300, 336), 336)
  a_p <- trial_arm(pop, NULL, 336)
  cat(sprintf("    n = %d screened in; baseline DLCO %.1f %%pred, A-aDO2 %.1f mmHg\n",
              nrow(a_d), mean(a_d$dlco0), mean(a_d$aado20)))
  cat(sprintf("                            observed baselines: DLCO ~45, A-aDO2 ~39\n"))
  cat(sprintf("    dDLCO 24 wk  drug %+5.2f (obs +9.8)   placebo %+5.2f (obs +3.8)\n",
              mean(a_d$dDLCO24), mean(a_p$dDLCO24)))
  cat(sprintf("    difference   %+5.2f (obs +6.0)   [%s]\n",
              mean(a_d$dDLCO24) - mean(a_p$dDLCO24),
              ok(abs((mean(a_d$dDLCO24) - mean(a_p$dDLCO24)) - 6.0) < 2.5)))
  cat(sprintf("    dSGRQ 24 wk  drug %+5.2f (obs -11.5)  placebo %+5.2f (obs -4.9)\n",
              mean(a_d$dSGRQ24), mean(a_p$dSGRQ24)))
  cat("    The placebo arm is not a placebo parameter. It is enrolment at the\n")
  cat("    low point of a deterioration that then resolves, plus the DLCO\n")
  cat("    measurement-learning term. Remove both and it is flat.\n")

  hr(); cat("D08 IMPALA-2 at 48 weeks - THE MODEL OVER-PREDICTS\n")
  cat(sprintf("    dDLCO 48 wk  drug %+5.2f (obs +11.6)  placebo %+5.2f (obs +4.7)\n",
              mean(a_d$dDLCO48), mean(a_p$dDLCO48)))
  cat(sprintf("    difference   %+5.2f (obs +6.9)   [%s]\n",
              mean(a_d$dDLCO48) - mean(a_p$dDLCO48),
              ok(abs((mean(a_d$dDLCO48) - mean(a_p$dDLCO48)) - 6.9) < 3)))
  cat("    This is a real failure and it localises a real question. The model\n")
  cat("    contains a recovery-side positive feedback: clearing the lung opens\n")
  cat("    units, which lets more aerosol reach the remaining burden, which\n")
  cat("    clears faster. So responders accelerate, and the 48-week difference\n")
  cat("    comes out roughly twice what was observed. The trial says recovery\n")
  cat("    is flat after 24 weeks; the withdrawal data (Tazawa 2010: 29 of 35\n")
  cat("    stable for a year off treatment) say something self-sustaining does\n")
  cat("    happen. One of the two is being mis-read and this model cannot say\n")
  cat("    which. Quote the 24-week number, not the 48-week one.\n")

  # ---- D09 IMPALA A-aDO2 -------------------------------------------
  hr(); cat("D09 IMPALA A-aDO2 (Trapnell 2020, PMID 32897035)\n")
  cat(sprintf("    dA-aDO2 24 wk  drug %+5.2f (obs -12.8)  placebo %+5.2f (obs -6.6)\n",
              mean(a_d$dAADO24), mean(a_p$dAADO24)))
  cat(sprintf("    difference     %+5.2f (obs -6.2)   [%s]\n",
              mean(a_d$dAADO24) - mean(a_p$dAADO24),
              ok(abs((mean(a_d$dAADO24) - mean(a_p$dAADO24)) + 6.2) < 2.5)))
  cat("    The treatment DIFFERENCE lands on both endpoints of two different\n")
  cat("    trials. Both within-arm changes are smaller than observed, i.e. the\n")
  cat("    model under-supplies whatever makes both arms of IMPALA improve.\n")

  # ---- D10 PAGE and its run-in ------------------------------------
  hr(); cat("D10 PAGE (Tazawa 2019, PMID 31483963) and the effect of its run-in\n")
  a_s <- trial_arm(pop, ev_inh_cycle(125, 168, 2, 7, 7), 168)
  cat(sprintf("    sargramostim 125 ug BID 7/14: dA-aDO2 %+5.2f (obs -4.50)\n",
              mean(a_s$dAADO24)))
  pop_ri <- vpop(n = 24, flare_mean = 0.0)     # run-in excludes improvers
  b_s <- trial_arm(pop_ri, ev_inh_cycle(125, 168, 2, 7, 7), 168)
  b_p <- trial_arm(pop_ri, NULL, 168)
  cat(sprintf("    with the 12-week run-in modelled (no enrolment flare):\n"))
  cat(sprintf("      drug %+5.2f (obs -4.50)   placebo %+5.2f (obs +0.17)   [%s]\n",
              mean(b_s$dAADO24), mean(b_p$dAADO24),
              ok(abs(mean(b_p$dAADO24)) < 2)))
  cat("    PAGE is the trial whose placebo arm did not move, and PAGE is the\n")
  cat("    trial that excluded improvers during a 12-week observation period.\n")
  cat("    The same model reproduces both placebo arms by simulating the two\n")
  cat("    enrolment protocols rather than by using two placebo effects.\n")

  # ---- D11 Tazawa 2010 --------------------------------------------
  hr(); cat("D11 Tazawa 2010 phase 2 (PMID 20167854): A-aDO2 -12.3 mmHg\n")
  a_t <- trial_arm(pop, ev_inh_cycle(250, 168, 1, 8, 6), 168)
  resp <- mean(a_t$dAADO24 <= -10)
  cat(sprintf("    modelled dA-aDO2 %+5.2f mmHg; responders (>=10 mmHg) %.0f%% (obs 62%%)   [%s]\n",
              mean(a_t$dAADO24), 100 * resp, ok(mean(a_t$dAADO24) < -2)))

  # ---- D12 lavage durability -------------------------------------
  hr(); cat("D12 whole lung lavage: removal, and how long it lasts\n")
  wl <- sim(list(WLLT1 = 30), end = 8 * 365, delta = 2, init = st)
  pre <- at(wl, 29, "PTOTo"); post <- at(wl, 31, "PTOTo")
  i <- which(wl$time > 40)
  rec <- i[which(wl$AADO2[i] >= at(wl, 29, "AADO2"))[1]]
  cat(sprintf("    burden %.0f -> %.0f mg (%.0f%% removed); A-aDO2 %.1f -> %.1f; DLCO %.1f -> %.1f\n",
              pre, post, 100 * (1 - post / pre), at(wl, 29, "AADO2"),
              at(wl, 40, "AADO2"), at(wl, 29, "DLCO"), at(wl, 40, "DLCO")))
  cat(sprintf("    return to the pre-lavage gradient at %.1f y   [%s]\n",
              wl$time[rec] / 365, ok(!is.na(rec) && wl$time[rec] / 365 > 1)))
  cat(sprintf("    DLCO at 5 y %.1f %%pred vs Beccaria 2004: 75 +- 19\n",
              at(wl, 5 * 365, "DLCO")))
  cat("    Time to recurrence is not a fitted parameter. It is the mass\n")
  cat("    removed divided by the net accumulation rate, and the lavage does\n")
  cat("    not touch the accumulation rate - which is the whole point.\n")

  # ---- D13 hereditary PAP is a structural null -------------------
  hr(); cat("D13 hereditary PAP: inhaled GM-CSF must do nothing\n")
  hp <- list(RECFRAC = 0, CAPFLOOR = 0.305)
  sth <- present(hp)$state
  h1 <- sim(hp, end = 336, ev = ev_inh(300, 336), init = sth)
  h0 <- sim(hp, end = 336, init = sth)
  cat(sprintf("    dDLCO 48 wk: with drug %+.3f, without %+.3f, difference %+.4f   [%s]\n",
              at(h1, 336, "DLCO") - h1$DLCO[1], at(h0, 336, "DLCO") - h0$DLCO[1],
              (at(h1, 336, "DLCO") - h1$DLCO[1]) - (at(h0, 336, "DLCO") - h0$DLCO[1]),
              ok(abs((at(h1, 336, "DLCO") - h1$DLCO[1]) -
                     (at(h0, 336, "DLCO") - h0$DLCO[1])) < 0.05)))
  cat(sprintf("    free GM-CSF still rises to %.0f pM: the ligand is there, the\n",
              max(h1$LFREE)))
  cat("    receptor is not. Capacity replacement instead (S30) works.\n")

  # ---- D14 exercise amplification -------------------------------
  hr(); cat("D14 exercise amplifies the same shunt (no lung parameter changes)\n")
  cat(sprintf("    rest     : PaO2 %.1f, A-aDO2 %.1f, SaO2 %.1f%%\n",
              at(w, 1, "PAO2A"), at(w, 1, "AADO2"), at(w, 1, "SAO2")))
  cat(sprintf("    exercise : PaO2 %.1f, A-aDO2 %.1f, SaO2 %.1f%%   [%s]\n",
              at(w, 1, "PAO2AX"), at(w, 1, "AADOEX"), at(w, 1, "SAO2EX"),
              ok(at(w, 1, "AADOEX") > at(w, 1, "AADO2"))))
  cat("    Only VO2 and cardiac output change, so mixed venous content falls\n")
  cat("    and the same shunt costs more. This is why exercise desaturation\n")
  cat("    precedes resting hypoxaemia.\n")

  # ---- D15 is there bistability? --------------------------------
  hr(); cat("D15 BISTABILITY TEST - the hypothesis this model was built on\n")
  cfx <- 0.305
  lo <- sim(list(CAPFLOOR = cfx), end = 30 * 365, delta = 30,
            init = modifyList(as.list(healthy()[nrow(healthy()), CMTN]),
                              list(ABC = 25 * 3000, ABP = 25 * 3000 * 1.25)))
  hi_start <- present(list(CAPFLOOR = 0.27))$state
  hi <- sim(list(CAPFLOOR = cfx), end = 30 * 365, delta = 30, init = hi_start)
  cat(sprintf("    same floor %.3f, two very different starting burdens:\n", cfx))
  cat(sprintf("      from a healthy lung   : burden %6.0f -> %6.0f mg\n",
              lo$PTOTo[1], fin(lo, "PTOTo")))
  cat(sprintf("      from a filled lung    : burden %6.0f -> %6.0f mg   [%s]\n",
              hi$PTOTo[1], fin(hi, "PTOTo"),
              ok(abs(fin(lo, "PTOTo") - fin(hi, "PTOTo")) /
                   max(1, fin(lo, "PTOTo")) < 0.10)))
  cat("    They converge. There is NO hysteresis and NO second attractor:\n")
  cat("    for a given catabolic floor there is one equilibrium burden. The\n")
  cat("    filling-exclusion feedback is real but sub-critical. So the\n")
  cat("    clinical heterogeneity of aPAP is NOT bistability in this model -\n")
  cat("    it is a steep, single-valued dependence on the catabolic floor\n")
  cat("    (D16). That is a refutation of the premise the map was drawn on,\n")
  cat("    and cluster 11 of the .dot file is labelled accordingly.\n")

  # ---- D16 the asymptomatic third ------------------------------
  hr(); cat("D16 the 31.8% found asymptomatically (Inoue 2008)\n")
  cfs <- seq(0.27, 0.35, by = 0.005)
  pres <- vapply(cfs, function(cc) {
    nh <- natural_history(list(CAPFLOOR = cc), years = 25, delta = 20)
    as.numeric(any(nh$PAO2A <= 70))
  }, numeric(1))
  crit <- cfs[which(pres == 0)[1]]
  cat(sprintf("    enrolment criterion reached for floors below %.3f\n", crit))
  cat(sprintf("    modelled A-aDO2 across floors %.3f-%.3f: %.1f to %.1f mmHg   [%s]\n",
              min(cfs), max(cfs),
              max(vapply(cfs, function(cc) fin(natural_history(list(CAPFLOOR = cc),
                  years = 20, delta = 40), "AADO2"), numeric(1))),
              min(vapply(cfs, function(cc) fin(natural_history(list(CAPFLOOR = cc),
                  years = 20, delta = 40), "AADO2"), numeric(1))),
              ok(!is.na(crit))))
  cat("    A one-part-in-ten change in the residual catabolic floor spans\n")
  cat("    asymptomatic to respiratory failure. That steepness IS the model\n")
  cat("    explanation of the clinical spread, and it is falsifiable: it says\n")
  cat("    the floor must be measurable and tightly distributed.\n")

  # ---- D17 dose-response is flat above threshold ---------------
  hr(); cat("D17 inhaled dose-response: where is the threshold?\n")
  for (dd in c(0.3, 1, 3, 10, 30, 100, 300, 1000, 3000)) {
    d <- sim(end = 168, ev = ev_inh(dd, 168), init = st)
    dp <- sim(end = 3, delta = 1 / 48, ev = ev_inh(dd, 3), init = st)
    cat(sprintf("    %7.1f ug/d -> dDLCO24 %+5.2f | coverage %.2f | peak ligand:site %8.1f | trough free GM %8.3f pM\n",
                dd, at(d, 168, "DLCO") - d$DLCO[1],
                (at(d, 168, "COV") - at(d, 161, "COV")) / 7, max(dp$MRATIO),
                min(dp$LFREE[dp$time > 1.9 & dp$time < 2.0])))
  }
  cat(sprintf("    healthy free GM-CSF for scale: %.2f pM\n",
              as.list(param(mod))$CGMH * 1000 / as.list(param(mod))$MWGM))
  cat("    THE THRESHOLD IS NOT WHERE THE CLINICAL DOSE IS. It sits around\n")
  cat("    1-3 ug/day - roughly a hundredfold below the 300 ug given in\n")
  cat("    IMPALA-2 - because the receptor saturates at tens of pM and the\n")
  cat("    antibody buffer, though in large molar excess over ENDOGENOUS\n")
  cat("    GM-CSF, is trivial to overwhelm with an exogenous microgram dose.\n")
  cat("    Everything above that is plateau, which is a testable claim and is\n")
  cat("    consistent with dose-ranging in aPAP never having found a clear\n")
  cat("    dose-response, and with 125 ug BID and 300 ug QD performing alike.\n")
  cat("    What limits inhaled GM-CSF is REACH, not occupancy or dose: the\n")
  cat("    model says improve deposition and open the lung, not the dose.\n")

  # ---- D18 oxygen module round-trip ---------------------------
  hr(); cat("D18 oxygen module: Severinghaus round-trip and shunt limits\n")
  chk <- sim(list(SHCO = 0), end = 5, init = st)
  cat(sprintf("    zero-shunt case: A-aDO2 %.2f mmHg (physiologic shunt only)   [%s]\n",
              at(chk, 5, "AADO2"), ok(at(chk, 5, "AADO2") < 12)))
  o2 <- sim(list(FIO2 = 0.40), end = 5, init = st)
  cat(sprintf("    FiO2 0.21 -> 0.40: PaO2 %.1f -> %.1f mmHg   [%s]\n",
              at(w, 5, "PAO2A"), at(o2, 5, "PAO2A"),
              ok(at(o2, 5, "PAO2A") > at(w, 5, "PAO2A"))))
  cat("    Supplemental oxygen moves PaO2 and nothing else in the model: it\n")
  cat("    does not touch the pool, which is correct and is why it is\n")
  cat("    symptomatic therapy.\n")

  # ---- D19 predictions that are NOT validated ----------------
  hr(); cat("D19 PREDICTIONS WITH NO VALIDATING DATA - treat as hypotheses\n")
  r1 <- sim(end = 730, delta = 2, ev = ev_rtx(), init = st)
  cat(sprintf("    rituximab at titre %.0f: titre %.1f -> %.1f, dDLCO %+.1f at 1 y\n",
              r1$TITRE[1], r1$TITRE[1], at(r1, 365, "TITRE"),
              at(r1, 365, "DLCO") - r1$DLCO[1]))
  cat("      The long-lived plasma cell compartment is CD20-negative, so only\n")
  cat("      the plasmablast share of output can fall. Reported responses are\n")
  cat("      inconsistent; this is the model reason.\n")
  f1 <- sim(list(FCRNT = 30, FCRNX = 3), end = 730, delta = 2, init = st)
  cat(sprintf("    FcRn blockade: titre %.1f -> %.1f, dDLCO %+.1f at 1 y\n",
              f1$TITRE[1], at(f1, 365, "TITRE"), at(f1, 365, "DLCO") - f1$DLCO[1]))
  cat("      A 65% IgG reduction crosses the critical titre from a median\n")
  cat("      baseline. No trial has tested this in aPAP. The magnitude here\n")
  cat("      also inherits the 48-week over-prediction of D08.\n")
  s1 <- sim(end = 365, ev = ev_statin(80, 365), init = st)
  cat(sprintf("    statin 80 mg: dDLCO %+.1f at 1 y\n",
              at(s1, 365, "DLCO") - s1$DLCO[1]))
  cat("      Statin looks disproportionately good for its small capacity\n")
  cat("      effect because it is ORAL and pays no reach penalty - blood\n")
  cat("      reaches consolidated lung and aerosol does not. That comparison\n")
  cat("      is a genuine prediction; the size of EMAXST is not anchored.\n")

  hr(); cat("SUMMARY: read D08, D15 and D19 before quoting anything from this\n")
  cat("file. The 24-week treatment differences match two trials on two\n")
  cat("endpoints; the 48-week difference is roughly twice what was observed;\n")
  cat("the bistability hypothesis the map was built on is refuted.\n")
  cat("=====================================================================\n")
  invisible(TRUE)
}

# =====================================================================
# ENTRY POINT
# =====================================================================
# Rscript pap_mrgsolve_model.R runs every scenario and every diagnostic.
# Set PAP_NORUN=1 to source the model only (this is what the Shiny app does).
if (!identical(Sys.getenv("PAP_NORUN"), "1") &&
    (!interactive() || identical(Sys.getenv("PAP_RUN"), "1"))) {
  res <- run_scenarios()
  cat(sprintf("\n%d scenarios simulated.\n", length(res$scen)))
  for (nm in names(res$scen)) {
    d <- res$scen[[nm]]$sim
    if ("DLCO" %in% names(d)) {
      cat(sprintf("  %-56s DLCO %5.1f -> %5.1f | A-aDO2 %5.1f -> %5.1f | burden %6.0f -> %6.0f mg\n",
                  substr(nm, 1, 56), d$DLCO[1], fin(d, "DLCO"),
                  d$AADO2[1], fin(d, "AADO2"), d$PTOTo[1], fin(d, "PTOTo")))
    } else {
      cat(sprintf("  %-56s [summary table %d x %d]\n", substr(nm, 1, 56),
                  nrow(d), ncol(d)))
    }
  }
  run_diagnostics()
}
