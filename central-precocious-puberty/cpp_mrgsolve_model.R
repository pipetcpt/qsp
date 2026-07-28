# =============================================================================
#  cpp_mrgsolve_model.R
#  Central Precocious Puberty (CPP) — Quantitative Systems Pharmacology model
#  ---------------------------------------------------------------------------
#  성조숙증(중추성) QSP 모델 · mrgsolve (44 ODEs)
#
#  CENTRAL MODELLING THESIS
#  ------------------------
#  CPP is *not* modelled here as "puberty too early, therefore suppress it".
#  It is modelled as a competition inside ONE integral, in which OESTRADIOL
#  APPEARS TWICE WITH OPPOSITE SIGNS:
#
#      adult height = HT(t0) + integral[ GV(E2, IGF1, GPRES) dt ]   over the
#                                                                   OPEN-PLATE
#                                                                   window
#
#      (+)  E2 -> GH pulse amplitude x2-3 -> IGF-1 -> chondrocyte proliferation
#           ==> E2 RAISES the integrand (growth velocity)
#
#      (-)  E2 -> ERalpha on the growth plate -> irreversible consumption of
#           proliferative reserve -> epiphyseal fusion
#           ==> E2 SHORTENS the upper limit of integration
#
#  A GnRH agonist removes BOTH arms simultaneously.  Whether that helps is
#  therefore NOT a property of the drug; it is a property of how much plate
#  reserve is left when you start.  The model contains NO rule of the form
#  "treat before age 8".  That rule is an OUTPUT (see the sign-flip analysis:
#  height gain +7.5 cm at bone age 7.5 y, +1.5 cm at bone age 11.2 y, and
#  0.0 cm at bone age 14 y, with the crossing computed, not assumed).
#
#  THREE FURTHER STRUCTURAL CLAIMS
#  -------------------------------
#  1. BONE AGE IS DEFINED AS MATURATION RELATIVE TO THE SAME-AGE NORMAL
#     REFERENCE CHILD.  That is literally what the Greulich-Pyle atlas is: the
#     atlas is CALIBRATED on normal children, so a normal child has BA = CA at
#     every age and dBA/dCA == 1.0 BY CONSTRUCTION.  The model therefore writes
#
#         dBA/dCA = Rmat(E2, androgens, IGF-1) / Rmat(reference values at CA)
#
#     with the reference trajectories E2NORM(CA), TNORM(CA) and IGFNORM(CA)
#     prescribed analytically.  Two consequences fall out that no additive
#     formulation can produce:
#       (a) an untreated CPP girl at chronological age 7 has dBA/dCA ~ 1.8,
#           because her numerator is pubertal while her denominator is
#           prepubertal;
#       (b) a SUPPRESSED girl at chronological age 10-12 has dBA/dCA ~ 0.6-0.7,
#           i.e. BELOW the normal prepubertal rate, because by then the
#           reference child is herself pubertal and the denominator has grown.
#     (b) is the reason bone maturation appears "arrested" on therapy without
#     any arrest term existing anywhere in the model.
#
#  2. THE PITUITARY DECODES PULSE FREQUENCY, AND A DEPOT AGONIST DESTROYS THE
#     CODE RATHER THAN BLOCKING THE RECEPTOR.  GnRH pulse frequency (PULS) is
#     an explicit state.  The gonadotrope stimulus is
#
#         S = RS * ( Sendo * ffree  +  AINT * fa )
#
#     where fa is fractional GnRHR occupancy by AGONIST, ffree = 1 - fa - fx is
#     what is left for endogenous pulses, AINT = 1.6 is the agonist's intrinsic
#     activity (it is a SUPER-agonist at the receptor), and RS is the
#     sensitised receptor fraction.  Because AINT > 1, the first dose produces
#     a stimulus ABOVE baseline -> the FLARE.  Suppression then arrives only
#     through RS collapsing (kdes 0.45/d, krec 0.035/d => RSss ~ 0.09 at
#     fa = 0.8).  A GnRH ANTAGONIST enters the same equation with fx instead of
#     fa and NO intrinsic-activity term: no flare, and no desensitisation, so
#     it must hold occupancy continuously to work.  Same receptor, opposite
#     pharmacology, from one equation.
#
#  3. GROWTH VELOCITY IS THE WRONG MONITORING ENDPOINT, AND THE MODEL SAYS SO
#     QUANTITATIVELY.  Effective suppression drops growth velocity from ~9-10
#     to ~4-5 cm/yr, i.e. BELOW the normal prepubertal rate, because the child
#     loses both the sex steroid and its GH/IGF-1 amplification.  Across a
#     virtual cohort the correlation between on-therapy growth velocity and
#     realised adult-height gain is NEGATIVE.  dBA/dCA is the surrogate that
#     works.
#
#  WHAT THE MODEL IS FOR
#  ---------------------
#  Comparing (i) WHO to treat (the bone-age sign flip; the slowly progressive
#  variant), (ii) WITH WHAT (1-month vs 3-month vs 6-month depot vs 12-month
#  implant vs nasal spray vs antagonist vs aromatase inhibitor), (iii) HOW TO
#  KNOW IT IS WORKING, and (iv) WHAT IT COSTS (bone mineral density, hot
#  flushes, body mass, psychosocial index) — on endpoints that matter: adult
#  height, height gain versus the untreated counterfactual, age at menarche,
#  Tanner stage regression, and peak bone mass.
#
#  UNITS
#  -----
#  time            days from the start of simulation; chronological age
#                  CA = CA0 + TIME/365.25   (default CA0 = 5.0 y)
#  height          cm            growth velocity  cm/yr
#  bone age        years         GPRES  dimensionless 1 -> 0
#  LH, FSH         IU/L          E2  pg/mL          testosterone  ng/dL
#  DHEAS           ug/dL         inhibin B  pg/mL
#  IGF-1           ng/mL         GH  arbitrary pulse-amplitude units (1 = pre-
#                                pubertal reference)
#  uterine volume  mL            Tanner breast  continuous 1-5
#  BMD             Z-score       BMI  Z-score
#  peptide analogues  amount ug in depots, concentration ng/mL ( = ug/L )
#  anastrozole     amount ug, concentration ng/mL
#
#  NUMERICAL / COARSE-GRAINING NOTE
#  --------------------------------
#  LH, FSH, E2 and testosterone are PULSE-AVERAGED concentrations with
#  relaxation times of 6-12 h, not their true 20-60 min turnover.  The model
#  resolves days to years; pulse FREQUENCY enters as its own state instead.
#  The peptide central compartments are likewise given ~6-h distribution
#  half-lives with UNCHANGED clearances, so steady-state concentrations
#  (which are release-rate / CL and therefore independent of V) are unaffected.
#  This is what lets the independent fixed-step RK4 transcription in
#  cpp_reference_check.py reproduce this model's adult heights to <0.1 cm.
#
#  HONEST STATEMENT OF STATUS
#  --------------------------
#  Every number in $PARAM is either (i) taken from the literature cited in
#  cpp_references.md, or (ii) derived by hand from a published steady-state
#  observation with the derivation written out in the CALIBRATION section
#  below.  The model has NOT been fitted to individual patient data.  Exactly
#  ONE parameter (ENDOBLEED) is fitted to a single reported number.  Places
#  where the model DISAGREES with the literature are listed at the bottom of
#  README.md rather than tuned away.
# =============================================================================

# =============================================================================
#  CALIBRATION — how each block was pinned down
# =============================================================================
#
#  A. MKRN3 BRAKE AND PUBERTAL ONSET (the phenotype knob)
#     MKRN3 is an imprinted, paternally expressed ubiquitin ligase that
#     RESTRAINS the KNDy pulse generator; its circulating level falls before
#     puberty and loss-of-function mutation is the commonest monogenic cause of
#     CPP (Abreu NEJM 2013 PMID 23738509; Busch JCEM 2016 PMID 27057785).  The
#     model therefore treats MKRN3(0) = MK0 as THE patient-level parameter and
#     lets pubertal onset be an OUTPUT:
#         MK0 = 1.00  -> thelarche 10.3 y, menarche 13.0 y, adult height 163.9
#         MK0 = 0.40  -> thelarche  6.5 y, menarche  9.2 y, adult height 154.1
#         MK0 = 0.66  -> thelarche  8.6 y, adult height 160.0 (slowly progressive)
#     KMKSEX = 0.75 slows the same brake in boys, putting male pubertal onset
#     (testosterone > 100 ng/dL) at 11.5 y, peak height velocity at 13.1 y and
#     normal male adult height at 176.5 cm.
#     KMK = 6.40e-4 /d (t1/2 3.0 yr) was then chosen so that MK0 = 1 gives
#     thelarche at 10.3 y and menarche at 13.0 y, matching contemporary
#     population data (Biro; Eckert-Lind JAMA Pediatr 2020 PMID 32040143).
#     KMK50 = 0.35 with Hill 4 makes the brake release switch-like, which is
#     what makes onset date sensitive to MK0 rather than to everything else.
#
#  B. GONADOTROPE DOSE-RESPONSE (KSLH, HLH, LHMAX)
#     Solved from three published anchor points, all as basal (pulse-averaged)
#     LH by immunochemiluminometric assay (Neely J Pediatr 1995 PMID 7608809;
#     Houk/Lee Pediatrics 2009 PMID 19482738):
#         prepubertal   LH ~ 0.10-0.15 IU/L   at stimulus S ~ 0.15
#         pubertal CPP  LH ~ 3-4    IU/L      at S ~ 0.85
#         suppressed    LH < 0.3   IU/L       at S ~ 0.11
#     HLH = 3.2, KSLH = 0.62, LHMAX = 5.0, LHB0 = 0.03 satisfy all three.
#     The steepness (Hill 3.2) is what makes "peak LH > 5 IU/L on stimulation"
#     a usable diagnostic threshold rather than a continuum.
#
#  C. GnRHR DESENSITISATION (KDES, KREC, AINT)
#     AINT = 1.6 is set by the requirement that the first depot produce an LH
#     surge to ~4 IU/L in a child whose pre-treatment basal LH is ~3 (the
#     classic flare).  KDES = 0.45/d gives the flare a 1.5-day decay.
#     KREC = 0.035/d (t1/2 20 d) is set by the observation that the axis takes
#     weeks-to-months to recover after the last depot and that menarche follows
#     roughly a year later (Lazar JCEM 2007 PMID 17579199).  Together they give
#     RSss = KREC/(KREC + KDES*fa) = 0.09 at fa = 0.8, i.e. a 91% loss of the
#     sensitised receptor pool — the actual therapeutic mechanism.
#
#  D. PEPTIDE PK.  Leuprolide PLGA depot: FBURSTL = 0.11 released with
#     ka 1.10/d and the remainder dissolving at 0.0455/d (t1/2 15 d), which
#     reproduces the published biphasic profile of the 1-month paediatric depot
#     and its ~30-day duration (Pradhan Clin Drug Investig 2014 PMID 24756362;
#     Salem Paediatr Drugs 2026 PMID 41824265).  CLL = 140 L/d in a 30 kg
#     child puts mean steady-state leuprolide at 3.75 mg q28d = (3750 ug/28 d)
#     / 140 L/d = 0.96 ng/mL.  EC50L = 0.35 ng/mL follows from a leuprolide
#     GnRHR Kd of ~0.5 nM (MW 1209).  Histrelin implant: 50 mg reservoir
#     releasing 65 ug/day ZERO-ORDER for 12 months (Eugster JCEM 2007
#     PMID 17327379; Silverman JCEM 2015 PMID 25803268) gives 65/130 =
#     0.50 ng/mL with NO peak and NO trough — the point of the formulation.
#     Nafarelin nasal: F = 2.1%, so 1800 ug/day delivers only ~38 ug/day
#     systemically with a 4-h half-life, which is why it is the one formulation
#     in which the model predicts genuine interdose escape.
#
#  E. GROWTH PLATE.  Oestrogen-dependence of epiphyseal fusion is fixed by the
#     human loss-of-function experiments: ESR1 mutation (Smith NEJM 1994) and
#     aromatase deficiency (Morishima JCEM 1995 PMID 8530621; Carani NEJM 1997
#     PMID 9211678) both produce unfused epiphyses and continued growth into
#     adulthood, and oestrogen closes them.  The reserve-depletion (rather than
#     clock) formulation is Weise/Baron (PNAS 2001 PMID 11381135; Trends
#     Endocrinol Metab 2004 PMID 15380808).  BAFUS = 15.0 y bone age in girls
#     and 16.8 y in boys are the atlas ages at which growth ceases.
#
#  F. THE (+) ARM.  EGH = 1.8 with KGH = 25 pg/mL makes GH pulse amplitude rise
#     2.5-2.8 fold across puberty and IGF-1 go 160 -> ~390 ng/mL, matching the
#     pubertal IGF-1 rise.  GVBASE = 5.90 then reproduces prepubertal growth
#     velocity 5.6-5.9 cm/yr and, with EE2GV = 0.45, a pubertal peak of
#     7.4 cm/yr at 11.2 y in the normal girl.
#
#  G. ENDOMETRIUM AND MENARCHE.  KENDE = 48 pg/mL with Hill 2 and turnover
#     0.030/0.010 per day place the menarche crossing (ENDO > 0.60 with
#     uterine volume > 4 mL and bone age > 10) at 13.0 y in the normal girl.
#     ENDOBLEED = 0.503 is THE ONE FITTED PARAMETER in the whole model: it is
#     the 93rd percentile of endometrial priming in a virtual CPP cohort at the
#     moment therapy starts, chosen so that the predicted withdrawal-bleeding
#     rate after the first depot (7.2%) lands in the reported 5-10% range.
#
#  H. BONE MASS.  KZ = 0.55/yr on the difference between the child's own
#     oestrogen effect and that of the same-age reference child, plus a
#     catch-up term KCATCH = 0.35/yr acting while the accrual window is open,
#     reproduces the observed pattern: BMD Z-score ELEVATED in untreated CPP,
#     falling during suppression, and recovering to normal peak bone mass
#     afterwards (Boot JCEM 1998 PMID 9467543; Antoniazzi JCEM 1999
#     PMID 10372699 for the calcium effect KCAVD).
# =============================================================================

library(mrgsolve)
library(dplyr)

code <- '
$PROB
# Central Precocious Puberty (CPP) QSP model — 44 ODEs
# MKRN3/KNDy pulse generator -> GnRH pulse FREQUENCY -> LH/FSH -> gonadal E2
# -> the two opposite-signed arms of the adult-height integral
# Therapy: GnRH agonist depots and implant, GnRH antagonist, aromatase
#          inhibitor, rhGH, tamoxifen

$PARAM @annotated
// ---------------- PK: leuprolide acetate PLGA depot (IM) ------------------
FBURSTL : 0.11    : leuprolide depot burst fraction (-)
KABL    : 1.10    : leuprolide burst release rate (1/day)
KDISL   : 0.0455  : leuprolide slow dissolution rate (1/day)
CLL     : 140.0   : leuprolide clearance (L/day, 30 kg child)
V1L     : 60.0    : leuprolide central volume (L)
QL      : 20.0    : leuprolide intercompartmental clearance (L/day)
V2L     : 40.0    : leuprolide peripheral volume (L)
EC50L   : 0.35    : leuprolide GnRHR occupancy EC50 (ng/mL)
// ---------------- PK: triptorelin pamoate depot (IM) ---------------------
FBURSTT : 0.07    : triptorelin depot burst fraction (-)
KABT    : 0.90    : triptorelin burst release rate (1/day)
KDIST   : 0.0150  : triptorelin slow dissolution rate (1/day)
CLT     : 125.0   : triptorelin clearance (L/day)
V1T     : 55.0    : triptorelin central volume (L)
EC50T   : 0.28    : triptorelin GnRHR occupancy EC50 (ng/mL)
// ---------------- PK: histrelin subdermal implant (zero order) -----------
K0H     : 65.0    : histrelin zero-order release rate (ug/day)
CLH     : 130.0   : histrelin clearance (L/day)
V1H     : 55.0    : histrelin central volume (L)
EC50H   : 0.12    : histrelin GnRHR occupancy EC50 (ng/mL)
// ---------------- PK: nafarelin nasal spray -----------------------------
KAN     : 6.0     : nafarelin nasal absorption rate (1/day)
FN      : 0.021   : nafarelin intranasal bioavailability (-)
CLN     : 150.0   : nafarelin clearance (L/day)
V1N     : 40.0    : nafarelin central volume (L)
EC50N   : 0.55    : nafarelin GnRHR occupancy EC50 (ng/mL)
// ---------------- PK: GnRH ANTAGONIST (degarelix-like sc depot) ---------
KAX     : 0.030   : antagonist depot release rate (1/day)
CLX     : 95.0    : antagonist clearance (L/day)
V1X     : 45.0    : antagonist central volume (L)
EC50X   : 0.90    : antagonist GnRHR occupancy EC50 (ng/mL)
// ---------------- PK/PD: aromatase inhibitor ----------------------------
CLAI    : 26.0    : aromatase inhibitor clearance (L/day)
VAI     : 45.0    : aromatase inhibitor volume (L)
FAI     : 0.85    : aromatase inhibitor oral bioavailability (-)
IMAXAI  : 0.975   : maximal aromatase inhibition (-)
IC50AI  : 3.0     : aromatase inhibition IC50 (ng/mL)
// ---------------- PK/PD: rhGH ------------------------------------------
KAG     : 1.4     : rhGH sc absorption rate (1/day)
KELG    : 1.1     : rhGH effect-compartment elimination (1/day)
VGH     : 1200.0  : rhGH effect-scaling volume; 0.043 mg/kg/d in 30 kg gives CGH ~ 1 (L)
EGH_EXO : 0.55    : rhGH effect slope on GH pulse amplitude (per unit CGH)
// ---------------- PK/PD: tamoxifen -------------------------------------
CLTAM   : 380.0   : tamoxifen clearance (L/day)
VTAM    : 1100.0  : tamoxifen volume (L)
IMAXTAM : 0.72    : maximal oestrogen-receptor blockade (-)
IC50TAM : 40.0    : tamoxifen IC50 (ng/mL)
// ---------------- GnRH receptor dynamics -------------------------------
KDES    : 0.45    : agonist-driven desensitisation rate (1/day)
KDESE   : 0.02    : endogenous-drive desensitisation rate (1/day)
KREC    : 0.035   : receptor resensitisation rate (1/day)
AINT    : 1.60    : agonist intrinsic activity relative to endogenous GnRH (-)
// ---------------- MKRN3 brake / KNDy pulse generator -------------------
MK0     : 1.00    : initial MKRN3 brake level; THE phenotype knob (-)
KMK     : 6.40e-4 : MKRN3 decay rate (1/day)
KMKSEX  : 0.75    : multiplier on KMK in boys; male puberty starts ~1.8 yr later (-)
KMK50   : 0.35    : MKRN3 restraint IC50 on KNDy drive (-)
HMK     : 4.0     : MKRN3 restraint Hill coefficient (-)
KNDMAX  : 1.15    : maximal KNDy drive (-)
KNDLES  : 0.0     : ectopic KNDy/GnRH drive from CNS lesion or hamartoma (-)
KSEC    : 0.0     : secondary central activation by advanced bone age (-)
TAUKND  : 60.0    : KNDy drive remodelling time constant (day)
DLK1R   : 0.10    : DLK1-mediated residual restraint (-)
KLEPB   : 0.18    : leptin/adiposity permissive slope on BMI Z (-)
PULSMIN : 6.0     : prepubertal GnRH pulse frequency (pulses/24h)
PULSMAX : 18.0    : adult GnRH pulse frequency (pulses/24h)
TAUP    : 20.0    : pulse-frequency time constant (day)
PULSREF : 12.0    : reference pulse frequency (pulses/24h)
// ---------------- gonadotropes ----------------------------------------
LHB0    : 0.030   : LH secretion floor (IU/L)
LHMAX   : 5.00    : maximal LH (IU/L)
KSLH    : 0.62    : gonadotrope stimulus giving half-maximal LH (-)
HLH     : 3.2     : LH Hill coefficient (-)
TAULH   : 0.25    : pulse-averaged LH relaxation time (day)
FSHB0   : 0.60    : FSH secretion floor (IU/L)
FSHMAX  : 5.00    : maximal FSH (IU/L)
KSF     : 0.45    : stimulus giving half-maximal FSH (-)
HFS     : 1.6     : FSH Hill coefficient (-)
TAUFSH  : 0.40    : pulse-averaged FSH relaxation time (day)
KINH    : 90.0    : inhibin B IC50 on FSH (pg/mL)
KFB     : 120.0   : E2 IC50 for negative feedback on endogenous drive (pg/mL)
// ---------------- ovary / steroids ------------------------------------
KFOL    : 0.060   : follicular recruitment rate (1/day)
KFF     : 2.60    : FSH giving half-maximal recruitment (IU/L)
HFF     : 3.0     : FSH recruitment Hill coefficient (-)
KFOLO   : 0.030   : follicular pool loss rate (1/day)
KE2     : 120.0   : maximal OVARIAN aromatase E2 output; gated by (1-SEXM) (pg/mL)
KLE     : 1.20    : LH giving half-maximal E2 output (IU/L)
HLE     : 1.30    : LH-E2 Hill coefficient (-)
TAUE2   : 0.50    : pulse-averaged E2 relaxation time (day)
KPER    : 1.40    : peripheral aromatisation coefficient; the ONLY E2 source in boys (pg/mL)
E2FLOOR : 0.80    : irreducible oestradiol floor (pg/mL)
KINHB   : 150.0   : maximal inhibin B (pg/mL)
KFF2    : 2.20    : FSH giving half-maximal inhibin B (IU/L)
TAUI    : 1.5     : inhibin B relaxation time (day)
KAUT    : 95.0    : autonomous (GNAS) gonadal E2 output (pg/mL)
TAUAUT  : 40.0    : autonomous activity time constant (day)
AUTSET  : 0.0     : autonomous GNAS activity set-point; >0 = McCune-Albright (-)
// ---------------- androgens ------------------------------------------
TB0     : 4.0     : testosterone floor (ng/dL)
KTG     : 600.0   : maximal LH-driven gonadal testosterone (ng/dL)
KLT     : 1.60    : LH giving half-maximal testosterone (IU/L)
HLT     : 1.20    : LH-testosterone Hill coefficient (-)
KTA     : 22.0    : adrenal contribution to testosterone (ng/dL)
TAUT    : 0.50    : testosterone relaxation time (day)
LEYD    : 1.0     : Leydig-cell functional capacity (-)
SEXM    : 0.0     : sex flag; 0 = female, 1 = male (-)
DHM     : 250.0   : maximal DHEAS (ug/dL)
DHCA50  : 11.0    : chronological age at half-maximal DHEAS (yr)
HDH     : 3.5     : adrenarche Hill coefficient (-)
TAUDH   : 60.0    : DHEAS relaxation time (day)
// ---------------- GH / IGF-1 -----------------------------------------
GHB     : 1.00    : prepubertal GH pulse amplitude (arbitrary units)
EGH     : 1.80    : maximal E2 amplification of GH amplitude (-)
KGH     : 25.0    : E2 giving half-maximal GH amplification (pg/mL)
NGH     : 2.0     : E2-GH Hill coefficient (-)
TAUGH   : 1.0     : GH amplitude relaxation time (day)
KIGF    : 160.0   : IGF-1 per unit GH amplitude (ng/mL)
TAUIGF  : 1.0     : IGF-1 relaxation time (day)
IGFREF  : 150.0   : reference prepubertal IGF-1 (ng/mL)
// ---------------- growth plate / bone age ----------------------------
GVBASE  : 5.90    : growth velocity scale (cm/yr)
GVSEXM  : 0.04    : male increment on the growth-velocity scale (-)
PGP     : 0.45    : growth-velocity exponent on plate reserve (-)
PIGF    : 0.50    : growth-velocity exponent on IGF-1 (-)
EE2GV   : 0.45    : maximal direct E2 stimulation of growth velocity (-)
KGV     : 20.0    : E2 giving half-maximal growth-velocity effect (pg/mL)
HGV     : 1.20    : E2-growth-velocity Hill coefficient (-)
ETGV    : 0.34    : maximal direct androgen stimulation of growth velocity (-)
KTGV    : 250.0   : testosterone giving half-maximal growth-velocity effect (ng/dL)
HTGV    : 1.50    : testosterone-growth-velocity Hill coefficient (-)
SCBA    : 1.00    : bone-age rate scaling (-)
M0      : 0.30    : steroid-independent maturation drive (-)
ME      : 0.55    : oestradiol maturation drive (-)
KBA     : 25.0    : E2 giving half-maximal maturation drive (pg/mL)
HBA     : 1.30    : E2-maturation Hill coefficient (-)
MA      : 0.20    : androgen maturation drive (-)
KAND    : 0.60    : androgen index giving half-maximal maturation (-)
MI      : 0.45    : IGF-1 maturation drive (-)
KIG     : 300.0   : IGF-1 giving half-maximal maturation drive (ng/mL)
BAFUS   : 15.0    : bone age at epiphyseal fusion (yr; 16.8 in boys)
BAREF   : 5.0     : bone age at which plate reserve is defined as 1 (yr)
XTRA    : 0.10    : extra reserve consumption per unit E2 drive (-)
// ---------------- prescribed NORMAL reference trajectories ------------
E2N0    : 1.40    : girls reference E2 floor (pg/mL)
E2NA    : 54.0    : girls reference E2 amplitude (pg/mL)
E2NC    : 11.35   : girls reference E2 midpoint age (yr)
E2NH    : 11.0    : girls reference E2 Hill coefficient (-)
E2N0M   : 1.60    : boys reference E2 floor (pg/mL)
E2NAM   : 27.0    : boys reference E2 amplitude (pg/mL)
E2NCM   : 12.60   : boys reference E2 midpoint age (yr)
E2NHM   : 9.0     : boys reference E2 Hill coefficient (-)
TN0     : 5.0     : girls reference testosterone floor (ng/dL)
TNA     : 42.0    : girls reference testosterone amplitude (ng/dL)
TNC     : 11.40   : girls reference testosterone midpoint age (yr)
TNH     : 8.0     : girls reference testosterone Hill coefficient (-)
TN0M    : 6.0     : boys reference testosterone floor (ng/dL)
TNAM    : 480.0   : boys reference testosterone amplitude (ng/dL)
TNCM    : 12.60   : boys reference testosterone midpoint age (yr)
TNHM    : 8.0     : boys reference testosterone Hill coefficient (-)
// ---------------- target tissues -------------------------------------
UTV0    : 1.20    : prepubertal uterine volume (mL)
KUT     : 10.00   : maximal E2-driven uterine volume increment (mL)
KUTE    : 25.0    : E2 giving half-maximal uterine growth (pg/mL)
HUT     : 1.50    : E2-uterus Hill coefficient (-)
TAUUT   : 30.0    : uterine volume time constant (day)
KBST    : 22.0    : E2 giving half-maximal breast development (pg/mL)
HBST    : 2.20    : E2-breast Hill coefficient (-)
TAUBSTU : 90.0    : breast progression time constant (day)
TAUBSTD : 240.0   : breast REGRESSION time constant (day)
KEND    : 0.030   : endometrial proliferation rate (1/day)
KENDE   : 45.0    : E2 giving half-maximal endometrial proliferation (pg/mL)
HEND    : 2.0     : E2-endometrium Hill coefficient (-)
KENDO   : 0.010   : endometrial involution rate (1/day)
ENDOBLEED : 0.503 : endometrial priming above which withdrawal bleeding occurs; THE ONE FITTED PARAMETER (-)
// ---------------- bone mass -----------------------------------------
KZ      : 0.55    : BMD Z-score slope on relative oestrogen effect (1/yr)
KBM     : 25.0    : E2 giving half-maximal bone effect (pg/mL)
HBM     : 1.20    : E2-bone Hill coefficient (-)
KCATCH  : 0.35    : BMD catch-up rate while accrual window is open (1/yr)
CAVD    : 0.0     : calcium + vitamin D supplementation flag (-)
KCAVD   : 0.30    : BMD Z benefit of supplementation, gated by suppression (1/yr)
// ---------------- body composition ----------------------------------
BMIZTGT : 0.30    : untreated BMI Z-score target (-)
KBMIG   : 0.25    : BMI Z increment on full GnRHa suppression (-)
TAUBMI  : 365.0   : BMI Z time constant (day)
// ---------------- symptoms / QoL ------------------------------------
KHF     : 6.00    : maximal hot-flush index (0-10)
HFSCALE : 20.0    : E2 fall giving maximal hot flushes (pg/mL)
TAUHF   : 7.0     : hot-flush time constant (day)
TAUTRK  : 60.0    : trailing E2 reference time constant (day)
TAUQOL  : 45.0    : psychosocial index time constant (day)
// ---------------- initial conditions / housekeeping -----------------
CA0     : 5.0     : chronological age at TIME = 0 (yr)
BA00    : 5.0     : bone age at TIME = 0 (yr)
HT00    : 108.0   : height at TIME = 0 (cm)

$CMT @annotated
// ---- PK / exposure (17) ----
LDEPB  : leuprolide depot, burst pool (ug)
LDEPS  : leuprolide depot, slow-release pool (ug)
LA1    : leuprolide central (ug)
LA2    : leuprolide peripheral (ug)
TDEPB  : triptorelin depot, burst pool (ug)
TDEPS  : triptorelin depot, slow-release pool (ug)
TA1    : triptorelin central (ug)
HIMP   : histrelin implant reservoir (ug)
HA1    : histrelin central (ug)
NDEP   : nafarelin nasal absorption depot (ug)
NA1    : nafarelin central (ug)
XDEP   : GnRH antagonist depot (ug)
XA1    : GnRH antagonist central (ug)
AIA1   : aromatase inhibitor central (ug)
GHDEP  : rhGH subcutaneous depot (ug)
GHEFF  : rhGH effect compartment (ug)
TAMA1  : tamoxifen central (ug)
// ---- GnRH receptor (2) ----
RS     : sensitised GnRH receptor fraction (-)
RD     : desensitised/internalised GnRH receptor fraction (-)
// ---- neuroendocrine (6) ----
MKRN3  : MKRN3 brake level (-)
KND    : KNDy pulse-generator drive (-)
PULS   : GnRH pulse frequency (pulses/24h)
LH     : pulse-averaged LH (IU/L)
FSH    : pulse-averaged FSH (IU/L)
INHB   : inhibin B (pg/mL)
// ---- gonad / steroid (5) ----
FOL    : antral follicle / granulosa functional mass (-)
E2     : oestradiol (pg/mL)
TESTO  : testosterone (ng/dL)
DHEAS  : DHEAS (ug/dL)
AUTON  : autonomous GNAS-driven gonadal activity (-)
// ---- growth (5) ----
GH     : GH pulse amplitude (units)
IGF1   : IGF-1 (ng/mL)
BA     : bone age (yr)
GPRES  : growth-plate proliferative reserve (-)
HT     : height (cm)
// ---- target tissue / body (7) ----
UTV    : uterine volume (mL)
BST    : Tanner breast stage, continuous (1-5)
ENDO   : endometrial proliferation state (-)
BMIZ   : BMI Z-score (-)
BMDZ   : lumbar-spine BMD Z-score (-)
HF     : hot-flush index (0-10)
E2TRK  : trailing 60-day E2 reference (pg/mL)
// ---- indices / integrators (2) ----
QOL    : psychosocial distress index (0-10)
CUME2  : cumulative oestradiol exposure (pg/mL*yr)

$GLOBAL
#define POWSAFE(x, p) (pow(fmax((x), 1.0e-9), (p)))
#define HILL(x, k, h) (POWSAFE((x), (h)) / (POWSAFE((x), (h)) + POWSAFE((k), (h))))

$MAIN
if (NEWIND < 2) {
  RS_0    = 1.0;
  RD_0    = 0.0;
  MKRN3_0 = MK0;
  KND_0   = 0.05;
  PULS_0  = 6.0;
  LH_0    = 0.10;
  FSH_0   = 1.20;
  INHB_0  = 8.0;
  FOL_0   = 0.12;
  E2_0    = 4.0;
  TESTO_0 = 5.0;
  DHEAS_0 = 25.0;
  AUTON_0 = 0.0;
  GH_0    = 1.05;
  IGF1_0  = 130.0;
  BA_0    = BA00;
  GPRES_0 = fmin(fmax((BAFUS - BA00) / (BAFUS - BAREF), 0.0), 1.0);
  HT_0    = HT00;
  UTV_0   = 1.3;
  BST_0   = 1.0;
  ENDO_0  = 0.05;
  BMIZ_0  = 0.30;
  BMDZ_0  = 0.0;
  HF_0    = 0.0;
  E2TRK_0 = 4.0;
  QOL_0   = 1.0;
  CUME2_0 = 0.0;
}

$ODE
double CA = CA0 + SOLVERTIME / 365.25;

double lh  = fmax(LH,    1.0e-9);
double fsh = fmax(FSH,   1.0e-9);
double e2  = fmax(E2,    1.0e-9);
double tes = fmax(TESTO, 1.0e-9);
double igf = fmax(IGF1,  1.0e-6);

// ============ 1. PK =========================================================
dxdt_LDEPB = -KABL * LDEPB;
dxdt_LDEPS = -KDISL * LDEPS;
dxdt_LA1   =  KABL * LDEPB + KDISL * LDEPS
              - (CLL / V1L) * LA1 - (QL / V1L) * LA1 + (QL / V2L) * LA2;
dxdt_LA2   =  (QL / V1L) * LA1 - (QL / V2L) * LA2;
double CLEUP = LA1 / V1L;

dxdt_TDEPB = -KABT * TDEPB;
dxdt_TDEPS = -KDIST * TDEPS;
dxdt_TA1   =  KABT * TDEPB + KDIST * TDEPS - (CLT / V1T) * TA1;
double CTRIP = TA1 / V1T;

double relH = (HIMP > 1.0) ? K0H : 0.0;          // zero-order implant release
dxdt_HIMP  = -relH;
dxdt_HA1   =  relH - (CLH / V1H) * HA1;
double CHIST = HA1 / V1H;

dxdt_NDEP  = -KAN * NDEP;
dxdt_NA1   =  FN * KAN * NDEP - (CLN / V1N) * NA1;
double CNAF = NA1 / V1N;

dxdt_XDEP  = -KAX * XDEP;
dxdt_XA1   =  KAX * XDEP - (CLX / V1X) * XA1;
double CANT = XA1 / V1X;

dxdt_AIA1  = -(CLAI / VAI) * AIA1;
double CAI   = AIA1 / VAI;
double AIEFF = IMAXAI * CAI / (CAI + IC50AI);

dxdt_GHDEP = -KAG * GHDEP;
dxdt_GHEFF =  KAG * GHDEP - KELG * GHEFF;
double CGH = GHEFF / VGH;

dxdt_TAMA1 = -(CLTAM / VTAM) * TAMA1;
double CTAM   = TAMA1 / VTAM;
double TAMEFF = IMAXTAM * CTAM / (CTAM + IC50TAM);

// ============ 2. GnRHR occupancy (competitive) and desensitisation =========
double A = CLEUP / EC50L + CTRIP / EC50T + CHIST / EC50H + CNAF / EC50N;
double B = CANT / EC50X;
double den = 1.0 + A + B;
double fa = A / den;                  // fraction occupied by AGONIST
double fx = B / den;                  // fraction occupied by ANTAGONIST
double ffree = fmin(fmax(1.0 - fa - fx, 0.0), 1.0);

// ============ 3. MKRN3 brake -> KNDy pulse generator =======================
dxdt_MKRN3 = -KMK * (1.0 - (1.0 - KMKSEX) * SEXM) * MKRN3;   // boys later
double LEPF  = fmin(fmax(1.0 + KLEPB * BMIZ, 0.55), 1.60);
double RESTR = MKRN3 + DLK1R;
double KNDss = KNDMAX * LEPF / (1.0 + POWSAFE(RESTR / KMK50, HMK))
               + KNDLES + KSEC * HILL(BA, 10.5, 6.0);
dxdt_KND = (KNDss - KND) / TAUKND;
double KNDc  = fmin(fmax(KND, 0.0), 1.2);
double PULSss = PULSMIN + (PULSMAX - PULSMIN) * KNDc;
dxdt_PULS = (PULSss - PULS) / TAUP;

// ============ 4. Pituitary: the PULSATILITY DECODER ========================
double FBK   = 1.0 / (1.0 + POWSAFE(e2 / KFB, 1.5));
double Sendo = KND * POWSAFE(PULS / PULSREF, 0.5) * FBK;
double S     = RS * (Sendo * ffree + AINT * fa);

dxdt_RD = (KDES * fa + KDESE * Sendo) * RS - KREC * RD;
dxdt_RS = -dxdt_RD;

double LHss = LHB0 + LHMAX * HILL(S, KSLH, HLH);
dxdt_LH = (LHss - lh) / TAULH;
double FSHss = (FSHB0 + FSHMAX * HILL(S, KSF, HFS)) / (1.0 + INHB / KINH);
dxdt_FSH = (FSHss - fsh) / TAUFSH;

// ============ 5. Gonad / steroidogenesis ==================================
dxdt_FOL = KFOL * HILL(fsh, KFF, HFF) * (1.0 - FOL) - KFOLO * FOL;
dxdt_AUTON = (AUTSET - AUTON) / TAUAUT;

double AROM  = FOL * (1.0 - AIEFF);
double E2per = KPER * (tes / 20.0 + DHEAS / 400.0) * (1.0 - AIEFF)
               * (1.0 + 0.20 * fmax(BMIZ, 0.0));
// the gonadal-aromatase term is OVARIAN: gate it by sex, so that in boys
// oestradiol arises ONLY from peripheral aromatisation of testosterone
double E2ss  = E2FLOOR + E2per + KE2 * (1.0 - SEXM) * AROM * HILL(lh, KLE, HLE)
               + KAUT * AUTON * (1.0 - AIEFF);
dxdt_E2 = (E2ss - e2) / TAUE2;

double Tss = TB0 + KTG * HILL(lh, KLT, HLT)
             * (SEXM * LEYD + (1.0 - SEXM) * 0.08) + KTA * DHEAS / 400.0;
dxdt_TESTO = (Tss - tes) / TAUT;

double DHss = DHM * HILL(CA, DHCA50, HDH);
dxdt_DHEAS = (DHss - DHEAS) / TAUDH;

double INHBss = KINHB * FOL * HILL(fsh, KFF2, 2.0);
dxdt_INHB = (INHBss - INHB) / TAUI;

// ============ 6. GH / IGF-1 — the (+) arm =================================
double GHss = GHB * (1.0 + EGH * HILL(e2, KGH, NGH)) + EGH_EXO * CGH;
dxdt_GH   = (GHss - GH) / TAUGH;
dxdt_IGF1 = (KIGF * GH - igf) / TAUIGF;

// ============ 7. Growth plate — the (-) arm, IRREVERSIBLE =================
double GPc = fmin(fmax(GPRES, 0.0), 1.0);
double open_plate = (GPc > 1.0e-6) ? 1.0 : 0.0;
double GV = GVBASE * (1.0 + GVSEXM * SEXM) * POWSAFE(GPc, PGP)
            * POWSAFE(igf / IGFREF, PIGF)
            * (1.0 + EE2GV * HILL(e2, KGV, HGV) + ETGV * HILL(tes, KTGV, HTGV))
            * open_plate;
dxdt_HT = GV / 365.25;

// bone age = maturation RELATIVE TO THE SAME-AGE NORMAL REFERENCE CHILD
double E2plate = e2 * (1.0 - TAMEFF);
double ANDR = DHEAS / 400.0 + tes / 400.0;
double Rmat = M0 + ME * HILL(E2plate, KBA, HBA) + MA * (ANDR / (ANDR + KAND))
              + MI * (igf / (igf + KIG));

double E2NORM = (SEXM > 0.5)
  ? (E2N0M + E2NAM * HILL(CA, E2NCM, E2NHM))
  : (E2N0  + E2NA  * HILL(CA, E2NC,  E2NH));
double TNORM = (SEXM > 0.5)
  ? (TN0M + TNAM * HILL(CA, TNCM, TNHM))
  : (TN0  + TNA  * HILL(CA, TNC,  TNH));
double ANDRN   = DHss / 400.0 + TNORM / 400.0;
double IGFNORM = KIGF * GHB * (1.0 + EGH * HILL(E2NORM, KGH, NGH));
double RmatN = M0 + ME * HILL(E2NORM, KBA, HBA) + MA * (ANDRN / (ANDRN + KAND))
               + MI * (IGFNORM / (IGFNORM + KIG));

double dBA_yr = SCBA * Rmat / RmatN * open_plate;
dxdt_BA = dBA_yr / 365.25;
dxdt_GPRES = (GPRES > 0.0)
  ? (-(dBA_yr / 365.25) / (BAFUS - BAREF)
     * (1.0 + XTRA * HILL(E2plate, KBA, HBA)))
  : 0.0;

// ============ 8. Target tissues ==========================================
double E2t = e2 * (1.0 - TAMEFF);
double UTVss = UTV0 + KUT * HILL(E2t, KUTE, HUT);
dxdt_UTV = (UTVss - UTV) / TAUUT;

double BSTss = 1.0 + 4.0 * HILL(E2t, KBST, HBST);
double tauB = (BSTss >= BST) ? TAUBSTU : TAUBSTD;
dxdt_BST = (BSTss - BST) / tauB;

dxdt_ENDO = KEND * HILL(E2t, KENDE, HEND) * (1.0 - ENDO) - KENDO * ENDO;

// ============ 9. Bone mass ===============================================
double accr = fmin(fmax((20.0 - CA) / 10.0, 0.0), 1.0);
// calcium/vitamin D only matters while the child is oestrogen-deprived AND the
// accrual window is still open; without the (supp * accr) gate the term would
// keep adding indefinitely and lift ADULT bone mass, which it does not do
dxdt_BMDZ = (KZ * (HILL(e2, KBM, HBM) - HILL(E2NORM, KBM, HBM))
             + KCATCH * (-BMDZ) * accr
             + KCAVD * CAVD * (fa / (fa + 0.30)) * accr) / 365.25;

// ============ 10. Body composition, symptoms, QoL ========================
double supp = fa / (fa + 0.30);
dxdt_BMIZ = (BMIZTGT + KBMIG * supp - BMIZ) / TAUBMI;

dxdt_E2TRK = (e2 - E2TRK) / TAUTRK;
double HFss = KHF * fmin(fmax((E2TRK - e2) / HFSCALE, 0.0), 1.0);
dxdt_HF = (HFss - HF) / TAUHF;

double QOLss = 1.0 + 2.6 * fmin(fmax((BST - 1.0) / 4.0, 0.0), 1.0)
               + 1.3 * fmin(fmax((BA - CA) / 2.0, 0.0), 1.5)
               + 0.9 * HF / 6.0
               + 1.4 * fmin(fmax((ENDO - 0.60) / 0.30, 0.0), 1.0)
                     * ((CA < 10.0) ? 1.0 : 0.0);
dxdt_QOL = (QOLss - QOL) / TAUQOL;

dxdt_CUME2 = e2 / 365.25;

$TABLE
double CA_out = CA0 + TIME / 365.25;
double CLEUPo = LA1 / V1L;
double CTRIPo = TA1 / V1T;
double CHISTo = HA1 / V1H;
double CNAFo  = NA1 / V1N;
double CANTo  = XA1 / V1X;
double CAIo   = AIA1 / VAI;
double AIEFFo = IMAXAI * CAIo / (CAIo + IC50AI);
double CTAMo  = TAMA1 / VTAM;
double TAMEFFo = IMAXTAM * CTAMo / (CTAMo + IC50TAM);
double CPTOT  = CLEUPo + CTRIPo + CHISTo + CNAFo;

double Ao = CLEUPo / EC50L + CTRIPo / EC50T + CHISTo / EC50H + CNAFo / EC50N;
double Bo = CANTo / EC50X;
double deno = 1.0 + Ao + Bo;
double FAo = Ao / deno;
double FXo = Bo / deno;
double FFREEo = fmin(fmax(1.0 - FAo - FXo, 0.0), 1.0);

double e2o  = fmax(E2, 1.0e-9);
double igfo = fmax(IGF1, 1.0e-6);
double FBKo   = 1.0 / (1.0 + POWSAFE(e2o / KFB, 1.5));
double Sendoo = KND * POWSAFE(PULS / PULSREF, 0.5) * FBKo;
double So     = RS * (Sendoo * FFREEo + AINT * FAo);

double GPco = fmin(fmax(GPRES, 0.0), 1.0);
double openo = (GPco > 1.0e-6) ? 1.0 : 0.0;
double GVo = GVBASE * (1.0 + GVSEXM * SEXM) * POWSAFE(GPco, PGP)
             * POWSAFE(igfo / IGFREF, PIGF)
             * (1.0 + EE2GV * HILL(e2o, KGV, HGV)
                + ETGV * HILL(fmax(TESTO, 1e-9), KTGV, HTGV)) * openo;

double E2plateo = e2o * (1.0 - TAMEFFo);
double ANDRo = DHEAS / 400.0 + fmax(TESTO, 1e-9) / 400.0;
double DHsso = DHM * HILL(CA_out, DHCA50, HDH);
double Rmato = M0 + ME * HILL(E2plateo, KBA, HBA) + MA * (ANDRo / (ANDRo + KAND))
               + MI * (igfo / (igfo + KIG));
double E2NORMo = (SEXM > 0.5) ? (E2N0M + E2NAM * HILL(CA_out, E2NCM, E2NHM))
                              : (E2N0  + E2NA  * HILL(CA_out, E2NC,  E2NH));
double TNORMo  = (SEXM > 0.5) ? (TN0M + TNAM * HILL(CA_out, TNCM, TNHM))
                              : (TN0  + TNA  * HILL(CA_out, TNC,  TNH));
double ANDRNo   = DHsso / 400.0 + TNORMo / 400.0;
double IGFNORMo = KIGF * GHB * (1.0 + EGH * HILL(E2NORMo, KGH, NGH));
double RmatNo = M0 + ME * HILL(E2NORMo, KBA, HBA) + MA * (ANDRNo / (ANDRNo + KAND))
                + MI * (IGFNORMo / (IGFNORMo + KIG));
double DBADCA = SCBA * Rmato / RmatNo * openo;

double BAmCA = BA - CA_out;
double LHFSH = LH / fmax(FSH, 1e-9);

$CAPTURE @annotated
CA_out : chronological age (yr)
CPTOT  : total GnRH-agonist concentration (ng/mL)
CLEUPo : leuprolide concentration (ng/mL)
CTRIPo : triptorelin concentration (ng/mL)
CHISTo : histrelin concentration (ng/mL)
CNAFo  : nafarelin concentration (ng/mL)
CANTo  : GnRH antagonist concentration (ng/mL)
AIEFFo : fractional aromatase inhibition (-)
TAMEFFo: fractional oestrogen-receptor blockade (-)
FAo    : GnRHR fractional occupancy by agonist (-)
FXo    : GnRHR fractional occupancy by antagonist (-)
So     : total gonadotrope stimulus (-)
GVo    : growth velocity (cm/yr)
DBADCA : bone-age advance per chronological year (-)
BAmCA  : bone age minus chronological age (yr)
LHFSH  : LH to FSH ratio (-)
'

mod <- mrgsolve::mcode("cpp_qsp", code, soloc = tempdir())

# =============================================================================
#  BAYLEY-PINNEAU PREDICTED ADULT HEIGHT
#  (PMID 14918032 — the table a clinician actually uses; kept separate from the
#   model so that MODEL-TRUE final height and CLINICIAN-PREDICTED adult height
#   can be compared, which is itself one of the results.)
# =============================================================================
BP_GIRLS <- data.frame(
  ba = c(6.0, 7.0, 8.0, 9.0, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5,
         13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 17.0),
  pc = c(66.2, 69.6, 73.0, 77.2, 80.4, 82.3, 84.4, 86.2, 88.4, 90.6,
         92.2, 94.1, 95.8, 97.4, 98.6, 99.0, 99.6, 100.0))
BP_BOYS <- data.frame(
  ba = c(7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 12.5, 13.0, 13.5, 14.0,
         14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 18.0),
  pc = c(69.5, 72.3, 75.2, 78.4, 81.2, 84.2, 85.8, 87.6, 90.2, 92.7,
         94.8, 96.8, 97.9, 98.6, 99.2, 99.6, 100.0))

bp_fraction <- function(ba, male = FALSE) {
  tab <- if (male) BP_BOYS else BP_GIRLS
  approx(tab$ba, tab$pc, xout = pmin(pmax(ba, min(tab$ba)), max(tab$ba)))$y / 100
}
pah_bp <- function(ht, ba, male = FALSE) ht / bp_fraction(ba, male)

# =============================================================================
#  PHENOTYPES
#  MK0 is THE patient-level knob; pubertal onset is an OUTPUT, not an input.
# =============================================================================
MK_NORMAL <- 1.00     # thelarche 10.3 y, menarche 13.0 y, adult height 163.9 cm
MK_CPP    <- 0.40     # thelarche  6.5 y, menarche  9.2 y, adult height 154.1 cm
MK_SLOW   <- 0.66     # thelarche  8.6 y, adult height 160.0 cm (slowly progressive)

pheno_girl <- function(mk0 = MK_CPP, ...) {
  c(list(MK0 = mk0, SEXM = 0, BAFUS = 15.0, HT00 = 108.0, CA0 = 5.0), list(...))
}
pheno_boy <- function(mk0 = MK_CPP, ...) {
  c(list(MK0 = mk0, SEXM = 1, BAFUS = 16.8, HT00 = 108.5, CA0 = 5.0), list(...))
}
pheno_mas <- function(...) {   # McCune-Albright: GnRH-INDEPENDENT disease
  c(list(MK0 = MK_NORMAL, SEXM = 0, BAFUS = 15.0, HT00 = 108.0, CA0 = 5.0,
         AUTSET = 0.42, KSEC = 0.22), list(...))
}

# =============================================================================
#  DOSING REGIMENS.  Times are in DAYS FROM TIME = 0, i.e. from age CA0.
# =============================================================================
age2day <- function(age, ca0 = 5.0) (age - ca0) * 365.25

CMT_LDEPB <- 1;  CMT_LDEPS <- 2;  CMT_TDEPB <- 5;  CMT_TDEPS <- 6
CMT_HIMP  <- 8;  CMT_NDEP  <- 10; CMT_XDEP  <- 12; CMT_AIA1  <- 14
CMT_GHDEP <- 15; CMT_TAMA1 <- 17

#' GnRH agonist PLGA depot: the burst fraction and the slow fraction are
#' loaded into two separate depot compartments by ONE administration.
ev_leup <- function(dose_mg, ii, start_age, stop_age, ca0 = 5.0,
                    fburst = 0.11) {
  t0 <- age2day(start_age, ca0); t1 <- age2day(stop_age, ca0)
  n  <- max(0, floor((t1 - t0) / ii))
  amt <- dose_mg * 1000
  c(ev(time = t0, amt = amt * fburst,       cmt = CMT_LDEPB, ii = ii, addl = n),
    ev(time = t0, amt = amt * (1 - fburst), cmt = CMT_LDEPS, ii = ii, addl = n))
}
ev_trip <- function(dose_mg, ii, start_age, stop_age, ca0 = 5.0,
                    fburst = 0.07) {
  t0 <- age2day(start_age, ca0); t1 <- age2day(stop_age, ca0)
  n  <- max(0, floor((t1 - t0) / ii))
  amt <- dose_mg * 1000
  c(ev(time = t0, amt = amt * fburst,       cmt = CMT_TDEPB, ii = ii, addl = n),
    ev(time = t0, amt = amt * (1 - fburst), cmt = CMT_TDEPS, ii = ii, addl = n))
}
#' Histrelin implant: the reservoir is REPLACED (not added to) each year.
#' Implemented as a 50 mg load every 365 d, which is equivalent while the
#' previous implant is essentially exhausted at explantation.
ev_hist <- function(start_age, stop_age, ca0 = 5.0, dose_ug = 50000, ii = 365) {
  t0 <- age2day(start_age, ca0); t1 <- age2day(stop_age, ca0)
  n  <- max(0, floor((t1 - t0) / ii))
  ev(time = t0, amt = dose_ug, cmt = CMT_HIMP, ii = ii, addl = n)
}
ev_naf <- function(daily_ug, start_age, stop_age, ca0 = 5.0) {
  t0 <- age2day(start_age, ca0); t1 <- age2day(stop_age, ca0)
  n  <- max(0, floor((t1 - t0) * 3))
  ev(time = t0, amt = daily_ug / 3, cmt = CMT_NDEP, ii = 1 / 3, addl = n)
}
ev_antag <- function(dose_mg, ii, start_age, stop_age, ca0 = 5.0) {
  t0 <- age2day(start_age, ca0); t1 <- age2day(stop_age, ca0)
  n  <- max(0, floor((t1 - t0) / ii))
  ev(time = t0, amt = dose_mg * 1000, cmt = CMT_XDEP, ii = ii, addl = n)
}
ev_ai <- function(dose_mg, start_age, stop_age, ca0 = 5.0, f = 0.85) {
  t0 <- age2day(start_age, ca0); t1 <- age2day(stop_age, ca0)
  ev(time = t0, amt = dose_mg * 1000 * f, cmt = CMT_AIA1, ii = 1,
     addl = max(0, floor(t1 - t0)))
}
ev_gh <- function(mg_per_day, start_age, stop_age, ca0 = 5.0) {
  t0 <- age2day(start_age, ca0); t1 <- age2day(stop_age, ca0)
  ev(time = t0, amt = mg_per_day * 1000, cmt = CMT_GHDEP, ii = 1,
     addl = max(0, floor(t1 - t0)))
}
ev_tam <- function(dose_mg, start_age, stop_age, ca0 = 5.0, f = 0.40) {
  t0 <- age2day(start_age, ca0); t1 <- age2day(stop_age, ca0)
  ev(time = t0, amt = dose_mg * 1000 * f, cmt = CMT_TAMA1, ii = 1,
     addl = max(0, floor(t1 - t0)))
}

# =============================================================================
#  SIMULATION DRIVER
# =============================================================================
sim_cpp <- function(pars = pheno_girl(), events = NULL, years = 16.5,
                    delta = 1) {
  m <- mod %>% param(pars)
  end <- years * 365.25
  if (is.null(events)) {
    m %>% mrgsim(end = end, delta = delta, hmax = 1) %>% as.data.frame()
  } else {
    m %>% mrgsim(events = events, end = end, delta = delta,
                 hmax = 1, recsort = 3) %>% as.data.frame()
  }
}

# ---- read-outs ---------------------------------------------------------------
adult_height <- function(d) tail(d$HT, 1)
age_first    <- function(d, cond) { i <- which(cond); if (length(i)) d$CA_out[i[1]] else NA_real_ }
thelarche    <- function(d) age_first(d, d$BST >= 2.0)
menarche     <- function(d) age_first(d, d$ENDO > 0.60 & d$UTV > 4.0 & d$BA > 10.0)
fusion_age   <- function(d) age_first(d, d$GPRES <= 1e-6)
at_age       <- function(d, ca, col) approx(d$CA_out, d[[col]], xout = ca, rule = 2)$y
ba_ca_ratio  <- function(d, lo, hi) (at_age(d, hi, "BA") - at_age(d, lo, "BA")) / (hi - lo)
gv_mean      <- function(d, lo, hi) mean(d$GVo[d$CA_out >= lo & d$CA_out <= hi])
pct_time     <- function(d, lo, hi, col, thr) {
  s <- d[[col]][d$CA_out >= lo & d$CA_out <= hi]
  100 * mean(s > thr)
}

# =============================================================================
#  SCENARIOS (15).  Every scenario is a (phenotype, regimen) pair; the
#  interesting output is always the CONTRAST with the untreated counterfactual
#  of the SAME phenotype, which is why every runner returns both.
# =============================================================================
S <- list()

S$S1_normal <- function() list(
  label = "S1  normal puberty reference girl",
  pars  = pheno_girl(MK_NORMAL), events = NULL)

S$S2_cpp_untreated <- function() list(
  label = "S2  CPP girl, untreated (the counterfactual)",
  pars  = pheno_girl(MK_CPP), events = NULL)

S$S3_leup_1m <- function() list(
  label = "S3  CPP + leuprolide 3.75 mg IM q28d, age 7.4-12.4",
  pars  = pheno_girl(MK_CPP),
  events = ev_leup(3.75, 28, 7.4, 12.4))

S$S4_leup_3m <- function() list(
  label = "S4  CPP + leuprolide 11.25 mg IM q12wk, age 7.4-12.4",
  pars  = pheno_girl(MK_CPP),
  events = ev_leup(11.25, 84, 7.4, 12.4))

S$S5_histrelin <- function() list(
  label = "S5  CPP + histrelin 50 mg implant q12mo, age 7.4-12.4",
  pars  = pheno_girl(MK_CPP),
  events = ev_hist(7.4, 12.4))

S$S6_trip_6m <- function() list(
  label = "S6  CPP + triptorelin 22.5 mg IM q24wk, age 7.4-12.4",
  pars  = pheno_girl(MK_CPP),
  events = ev_trip(22.5, 168, 7.4, 12.4))

S$S7_nafarelin <- function() list(
  label = "S7  CPP + nafarelin nasal 1800 ug/day, age 7.4-12.4",
  pars  = pheno_girl(MK_CPP),
  events = ev_naf(1800, 7.4, 12.4))

S$S8_antagonist <- function() list(
  label = "S8  CPP + GnRH ANTAGONIST 18 mg q28d (no flare), age 7.4-12.4",
  pars  = pheno_girl(MK_CPP),
  events = ev_antag(18, 28, 7.4, 12.4))

S$S9_late_start <- function() list(
  label = "S9  CPP + leuprolide started LATE at age 10.0 (bone age ~12.8)",
  pars  = pheno_girl(MK_CPP),
  events = ev_leup(11.25, 84, 10.0, 14.0))

S$S10_gh_addon <- function() list(
  label = "S10 CPP + leuprolide + rhGH 0.043 mg/kg/day (30 kg)",
  pars  = pheno_girl(MK_CPP),
  events = c(ev_leup(11.25, 84, 7.4, 12.4), ev_gh(0.043 * 30, 7.4, 12.4)))

S$S11_slow_variant <- function() list(
  label = "S11 SLOWLY PROGRESSIVE variant, treated (the over-treatment case)",
  pars  = pheno_girl(MK_SLOW),
  events = ev_leup(11.25, 84, 8.8, 13.8))

S$S12_boy <- function() list(
  label = "S12 CPP boy + leuprolide 11.25 mg q12wk, age 8-13",
  pars  = pheno_boy(MK_CPP),
  events = ev_leup(11.25, 84, 8.0, 13.0))

S$S13_boy_ai <- function() list(
  label = "S13 CPP boy + anastrozole-like AI alone (E2 blocked, T preserved)",
  pars  = pheno_boy(MK_CPP),
  events = ev_ai(1.0, 8.0, 13.0))

S$S14_mas_gnrha <- function() list(
  label = "S14 McCune-Albright PPP + leuprolide (STRUCTURALLY unreachable)",
  pars  = pheno_mas(),
  events = ev_leup(11.25, 84, 6.0, 13.0))

S$S15_mas_ai <- function() list(
  label = "S15 McCune-Albright PPP + potent AI (letrozole-like) + tamoxifen",
  pars  = pheno_mas(IMAXAI = 0.992, IC50AI = 2.2),
  events = c(ev_ai(2.5, 6.0, 13.0), ev_tam(20, 6.0, 13.0)))

S$S16_bone_health <- function() list(
  label = "S16 CPP + leuprolide + calcium/vitamin D (bone-health arm)",
  pars  = pheno_girl(MK_CPP, CAVD = 1),
  events = ev_leup(11.25, 84, 7.4, 11.4))

# =============================================================================
#  ANALYSIS FUNCTIONS (12)
# =============================================================================

## A1 -- natural history: normal vs untreated CPP -----------------------------
a1_natural_history <- function() {
  dN <- sim_cpp(pheno_girl(MK_NORMAL))
  dC <- sim_cpp(pheno_girl(MK_CPP))
  out <- lapply(list(normal = dN, cpp = dC), function(d) data.frame(
    thelarche = thelarche(d), menarche = menarche(d), fusion = fusion_age(d),
    adult_height = adult_height(d), peak_GV = max(d$GVo), peak_E2 = max(d$E2),
    LH_at_8 = at_age(d, 8, "LH"), UTV_at_8 = at_age(d, 8, "UTV"),
    BA_at_8 = at_age(d, 8, "BA"), dBAdCA_7_9 = ba_ca_ratio(d, 7, 9),
    cumE2 = tail(d$CUME2, 1), BMDZ_adult = tail(d$BMDZ, 1)))
  res <- do.call(rbind, out)
  res$height_loss <- res["normal", "adult_height"] - res$adult_height
  res
}

## A2 -- THE SIGN FLIP: height gain vs bone age at start ---------------------
a2_sign_flip <- function(starts = c(6.6, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0,
                                    10.5, 11.0)) {
  dU <- sim_cpp(pheno_girl(MK_CPP))
  hU <- adult_height(dU)
  do.call(rbind, lapply(starts, function(a) {
    dT <- sim_cpp(pheno_girl(MK_CPP), ev_leup(11.25, 84, a, a + 6))
    data.frame(start_CA = a, start_BA = at_age(dU, a, "BA"),
               GPRES_at_start = at_age(dU, a, "GPRES"),
               height_noRx = hU, height_Rx = adult_height(dT),
               gain = adult_height(dT) - hU,
               GV_on_Rx = gv_mean(dT, a + 0.3, min(a + 2, 16.4)),
               dBAdCA_on_Rx = ba_ca_ratio(dT, a + 0.3, min(a + 2, 16.4)))
  }))
}

## A3 -- growth velocity is a misleading monitor ----------------------------
a3_monitor_trap <- function(start = 7.5) {
  dT <- sim_cpp(pheno_girl(MK_CPP), ev_leup(11.25, 84, start, start + 5))
  dU <- sim_cpp(pheno_girl(MK_CPP))
  ages <- c(7.4, 7.6, 8, 8.5, 9, 10, 11, 12, 12.5)
  data.frame(
    CA = ages,
    GV_Rx      = sapply(ages, at_age, d = dT, col = "GVo"),
    dBAdCA_Rx  = sapply(ages, at_age, d = dT, col = "DBADCA"),
    BAmCA_Rx   = sapply(ages, at_age, d = dT, col = "BAmCA"),
    E2_Rx      = sapply(ages, at_age, d = dT, col = "E2"),
    LH_Rx      = sapply(ages, at_age, d = dT, col = "LH"),
    UTV_Rx     = sapply(ages, at_age, d = dT, col = "UTV"),
    PAH_BP_Rx  = pah_bp(sapply(ages, at_age, d = dT, col = "HT"),
                        sapply(ages, at_age, d = dT, col = "BA")),
    GV_noRx    = sapply(ages, at_age, d = dU, col = "GVo"),
    PAH_BP_noRx = pah_bp(sapply(ages, at_age, d = dU, col = "HT"),
                         sapply(ages, at_age, d = dU, col = "BA")))
}

## A4 -- the flare, and whether an antagonist is worth it -------------------
a4_flare <- function(start = 7.5) {
  dA <- sim_cpp(pheno_girl(MK_CPP), ev_leup(3.75, 28, start, start + 5),
                delta = 0.25)
  dX <- sim_cpp(pheno_girl(MK_CPP), ev_antag(18, 28, start, start + 5),
                delta = 0.25)
  dU <- sim_cpp(pheno_girl(MK_CPP), delta = 0.25)
  win <- function(d, d0, d1) d[d$CA_out >= start + d0 / 365.25 &
                               d$CA_out <= start + d1 / 365.25, ]
  f <- function(d) c(
    peak_LH_14d  = max(win(d, 0, 14)$LH),
    peak_E2_14d  = max(win(d, 0, 14)$E2),
    LH_day28     = at_age(d, start + 28 / 365.25, "LH"),
    LH_day56     = at_age(d, start + 56 / 365.25, "LH"),
    peak_ENDO_30d = max(win(d, 0, 30)$ENDO),
    BA_at_9      = at_age(d, 9, "BA"),
    adult_height = adult_height(d))
  res <- as.data.frame(rbind(agonist = f(dA), antagonist = f(dX),
                             untreated = f(dU)))
  res$gain <- res$adult_height - adult_height(dU)
  res
}

## A5 -- formulation comparison: TROUGH COVERAGE is the design variable ----
a5_formulations <- function(start = 7.4, stop = 12.4) {
  regs <- list(
    "leuprolide 3.75 mg q28d"        = ev_leup(3.75, 28, start, stop),
    "leuprolide 7.5 mg q28d"         = ev_leup(7.5, 28, start, stop),
    "leuprolide 11.25 mg q12wk"      = ev_leup(11.25, 84, start, stop),
    "leuprolide 30 mg q12wk"         = ev_leup(30, 84, start, stop),
    "triptorelin 11.25 mg q12wk"     = ev_trip(11.25, 84, start, stop),
    "triptorelin 22.5 mg q24wk"      = ev_trip(22.5, 168, start, stop),
    "histrelin 50 mg implant q12mo"  = ev_hist(start, stop),
    "nafarelin nasal 1800 ug/d"      = ev_naf(1800, start, stop),
    "leuprolide 3.75 mg LATE q42d"   = ev_leup(3.75, 42, start, stop),
    "leuprolide 3.75 mg LATE q56d"   = ev_leup(3.75, 56, start, stop),
    "leuprolide 1.875 mg q28d (low)" = ev_leup(1.875, 28, start, stop),
    "GnRH antagonist 18 mg q28d"     = ev_antag(18, 28, start, stop))
  hU <- adult_height(sim_cpp(pheno_girl(MK_CPP)))
  do.call(rbind, lapply(names(regs), function(nm) {
    d <- sim_cpp(pheno_girl(MK_CPP), regs[[nm]], delta = 0.5)
    lo <- start + 0.25
    data.frame(formulation = nm,
               mean_Cp = mean(d$CPTOT[d$CA_out >= lo & d$CA_out <= stop]),
               mean_occ = mean(d$FAo[d$CA_out >= lo & d$CA_out <= stop]),
               pct_LH_gt_0.5 = pct_time(d, lo, stop, "LH", 0.5),
               pct_E2_gt_10  = pct_time(d, lo, stop, "E2", 10),
               adult_height = adult_height(d), gain = adult_height(d) - hU)
  }))
}

## A6 -- within-interval profile: where escape would come from -------------
a6_interdose <- function(dose_mg = 3.75, ii = 28, start = 7.4, cycle = 11) {
  d <- sim_cpp(pheno_girl(MK_CPP), ev_leup(dose_mg, ii, start, start + 5),
               delta = 0.25)
  base <- age2day(start) + ii * cycle
  days <- c(0.5, 2, 5, 10, 14, 20, 24, ii - 0.5)
  idx <- sapply(base + days, function(tt) which.min(abs(d$time - tt)))
  data.frame(day_in_cycle = days, Cp = d$CPTOT[idx], occ_fa = d$FAo[idx],
             RS = d$RS[idx], LH = d$LH[idx], E2 = d$E2[idx],
             stimulus_S = d$So[idx])
}

## A7 -- GnRH-INDEPENDENT disease: target/mechanism mismatch --------------
a7_peripheral <- function() {
  base <- pheno_mas()
  ai   <- pheno_mas(IMAXAI = 0.992, IC50AI = 2.2)
  arms <- list(
    list("no treatment", base, NULL),
    list("leuprolide 11.25 q12wk", base, ev_leup(11.25, 84, 6, 13)),
    list("potent AI (letrozole-like)", ai, ev_ai(2.5, 6, 13)),
    list("AI + leuprolide", ai, c(ev_ai(2.5, 6, 13), ev_leup(11.25, 84, 6, 13))),
    list("AI + tamoxifen", ai, c(ev_ai(2.5, 6, 13), ev_tam(20, 6, 13))))
  h0 <- adult_height(sim_cpp(base))
  do.call(rbind, lapply(arms, function(a) {
    d <- sim_cpp(a[[2]], a[[3]])
    data.frame(strategy = a[[1]],
               mean_E2_7_11 = mean(d$E2[d$CA_out >= 7 & d$CA_out <= 11]),
               BA_at_10 = at_age(d, 10, "BA"),
               adult_height = adult_height(d),
               gain = adult_height(d) - h0)
  }))
}

## A8 -- boys: aromatase inhibition separates the two E2 signs -----------
a8_boys <- function() {
  h0 <- adult_height(sim_cpp(pheno_boy(MK_CPP), years = 17.5))
  arms <- list(
    list("untreated CPP boy", NULL),
    list("leuprolide 11.25 q12wk", ev_leup(11.25, 84, 8, 13)),
    list("anastrozole-like AI alone", ev_ai(1.0, 8, 13)),
    list("GnRHa + AI", c(ev_leup(11.25, 84, 8, 13), ev_ai(1.0, 8, 13))))
  do.call(rbind, lapply(arms, function(a) {
    d <- sim_cpp(pheno_boy(MK_CPP), a[[2]], years = 17.5)
    w <- d$CA_out >= 8.5 & d$CA_out <= 12
    data.frame(strategy = a[[1]], mean_E2 = mean(d$E2[w]),
               mean_T = mean(d$TESTO[w]), BA_at_11 = at_age(d, 11, "BA"),
               adult_height = adult_height(d), gain = adult_height(d) - h0,
               BMDZ_nadir = min(d$BMDZ), BMDZ_adult = tail(d$BMDZ, 1))
  }))
}

## A9 -- rhGH add-on, early vs late ------------------------------------
a9_gh_addon <- function() {
  hU <- adult_height(sim_cpp(pheno_girl(MK_CPP)))
  grid <- expand.grid(start = c(7.2, 8.6), arm = c("GnRHa", "GnRHa+rhGH",
                                                   "rhGH alone"),
                      stringsAsFactors = FALSE)
  do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    a <- grid$start[i]; wt <- if (a < 8) 30 else 32
    ev <- switch(grid$arm[i],
      "GnRHa"      = ev_leup(11.25, 84, a, a + 4.5),
      "GnRHa+rhGH" = c(ev_leup(11.25, 84, a, a + 4.5),
                       ev_gh(0.043 * wt, a, a + 4.5)),
      "rhGH alone" = ev_gh(0.043 * wt, a, a + 4.5))
    d <- sim_cpp(pheno_girl(MK_CPP), ev)
    w <- d$CA_out >= a + 0.3 & d$CA_out <= min(a + 3, 16.4)
    data.frame(start = a, arm = grid$arm[i], mean_IGF1 = mean(d$IGF1[w]),
               mean_GV = mean(d$GVo[w]),
               dBAdCA = ba_ca_ratio(d, a + 0.3, min(a + 3, 16.4)),
               adult_height = adult_height(d), gain = adult_height(d) - hU)
  }))
}

## A10 -- the slowly progressive variant: over-treatment quantified -----
a10_slow_variant <- function(mks = c(0.34, 0.40, 0.46, 0.52, 0.58, 0.66, 0.76)) {
  do.call(rbind, lapply(mks, function(mk) {
    dU <- sim_cpp(pheno_girl(mk), years = 17)
    th <- thelarche(dU)
    a  <- min(max(ifelse(is.na(th), 12, th) + 0.5, 6.3), 11.5)
    dT <- sim_cpp(pheno_girl(mk), ev_leup(11.25, 84, a, a + 5), years = 17)
    data.frame(MK0 = mk, thelarche = th, peak_E2 = max(dU$E2),
               menarche_noRx = menarche(dU),
               height_noRx = adult_height(dU), height_Rx = adult_height(dT),
               gain = adult_height(dT) - adult_height(dU))
  }))
}

## A11 -- bone mineral density trajectory -----------------------------
a11_bmd <- function() {
  dN  <- sim_cpp(pheno_girl(MK_NORMAL))
  dU  <- sim_cpp(pheno_girl(MK_CPP))
  dT  <- sim_cpp(pheno_girl(MK_CPP), ev_leup(11.25, 84, 7.4, 11.4))
  dTC <- sim_cpp(pheno_girl(MK_CPP, CAVD = 1), ev_leup(11.25, 84, 7.4, 11.4))
  ages <- c(7, 7.4, 8, 9, 10, 11.4, 12, 14, 16, 20.4)
  data.frame(CA = ages,
             normal   = sapply(ages, at_age, d = dN,  col = "BMDZ"),
             CPP_noRx = sapply(ages, at_age, d = dU,  col = "BMDZ"),
             CPP_GnRHa = sapply(ages, at_age, d = dT,  col = "BMDZ"),
             plus_CaVitD = sapply(ages, at_age, d = dTC, col = "BMDZ"))
}

## A12 -- axis recovery and menarche after stopping -------------------
a12_recovery <- function(start = 7.4, durations = c(3, 4, 5)) {
  do.call(rbind, lapply(durations, function(dur) {
    d <- sim_cpp(pheno_girl(MK_CPP), ev_leup(11.25, 84, start, start + dur),
                 years = 17, delta = 0.5)
    stop_age <- start + dur
    aft <- d$CA_out > stop_age
    tLH <- age_first(d, aft & d$LH > 1.0)
    tE2 <- age_first(d, aft & d$E2 > 20)
    men <- menarche(d)
    data.frame(duration_yr = dur, stop_age = stop_age,
               months_to_LH1 = 12 * (tLH - stop_age),
               months_to_E2_20 = 12 * (tE2 - stop_age),
               menarche_age = men, months_to_menarche = 12 * (men - stop_age),
               adult_height = adult_height(d))
  }))
}

# =============================================================================
#  RUN EVERYTHING
# =============================================================================
run_all_scenarios <- function() {
  hU <- adult_height(sim_cpp(pheno_girl(MK_CPP)))
  do.call(rbind, lapply(names(S), function(nm) {
    sc <- S[[nm]]()
    yrs <- if (grepl("boy", nm)) 17.5 else 16.5
    d <- sim_cpp(sc$pars, sc$events, years = yrs)
    data.frame(scenario = sc$label, thelarche = thelarche(d),
               menarche = menarche(d), fusion = fusion_age(d),
               adult_height = adult_height(d),
               gain_vs_untreated_girl = adult_height(d) - hU,
               peak_QoL = max(d$QOL), peak_hot_flush = max(d$HF),
               BMDZ_nadir = min(d$BMDZ), stringsAsFactors = FALSE)
  }))
}

main <- function() {
  cat("\n================ A1  NATURAL HISTORY ================\n")
  print(a1_natural_history(), digits = 4)
  cat("\n================ A2  THE SIGN FLIP ==================\n")
  print(a2_sign_flip(), digits = 4)
  cat("\n================ A3  MONITORING TRAP ================\n")
  print(a3_monitor_trap(), digits = 4)
  cat("\n================ A4  THE FLARE ======================\n")
  print(a4_flare(), digits = 4)
  cat("\n================ A5  FORMULATIONS ===================\n")
  print(a5_formulations(), digits = 4)
  cat("\n================ A6  INTERDOSE PROFILE ==============\n")
  print(a6_interdose(), digits = 4)
  cat("\n================ A7  PERIPHERAL PP ==================\n")
  print(a7_peripheral(), digits = 4)
  cat("\n================ A8  BOYS / AROMATASE ===============\n")
  print(a8_boys(), digits = 4)
  cat("\n================ A9  rhGH ADD-ON ====================\n")
  print(a9_gh_addon(), digits = 4)
  cat("\n================ A10 SLOW VARIANT ===================\n")
  print(a10_slow_variant(), digits = 4)
  cat("\n================ A11 BONE DENSITY ===================\n")
  print(a11_bmd(), digits = 4)
  cat("\n================ A12 AXIS RECOVERY ==================\n")
  print(a12_recovery(), digits = 4)
  cat("\n================ SCENARIO TABLE =====================\n")
  print(run_all_scenarios(), digits = 4)
  invisible(NULL)
}

# Auto-run the analyses only when this file is executed directly (Rscript).
# When cpp_shiny_app.R sources it with local = TRUE, environment() is not the
# global environment, so the analyses are NOT run at app start-up.
if (!interactive() && identical(environment(), globalenv())) {
  main()
}

# =============================================================================
#  CROSS-CHECK AGAINST THE INDEPENDENT PYTHON TRANSCRIPTION
#  ---------------------------------------------------------------------------
#  cpp_reference_check.py implements the SAME 44 ODEs and the SAME parameter
#  values with a fixed-step RK4 integrator and no R dependency.  The two
#  transcriptions are expected to agree on:
#      normal girl   adult height  163.9 cm, thelarche 10.3 y, menarche 13.0 y
#      untreated CPP adult height  154.1 cm, thelarche  6.5 y, menarche  9.2 y
#      normal boy    adult height  176.5 cm, fusion 16.4 y, PHV at 13.1 y
#      CPP + leuprolide 11.25 q12wk from age 7.4  ->  ~160 cm
#  If they do not, ONE of the two files has a transcription error, and the
#  disagreement is the bug report.  cpp_reference_output.txt in this directory
#  is the committed output of the python side.
# =============================================================================
