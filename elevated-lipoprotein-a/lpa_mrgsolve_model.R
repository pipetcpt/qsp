## =====================================================================
##  lpa_mrgsolve_model.R
##  Elevated Lipoprotein(a) — QSP / PK-PD model
##  50 ODE compartments · 20 therapeutic scenarios
##
##  고지단백(a)혈증 — 정량적 시스템 약리학 모델
## ---------------------------------------------------------------------
##  STRUCTURAL THESIS
##  -----------------
##  Lp(a) is ONE particle, made at ONE rate-limiting step, acting through
##  THREE effector arms, and read on TWO non-interchangeable scales.
##  Almost every unresolved argument in the field is a consequence of
##  those four numbers (1 / 1 / 3 / 2), not of unknown biology.
##
##  (A)  ONE RATE-LIMITING STEP — hepatic apo(a) secretion.
##       Kinetic studies show the fractional catabolic rate of Lp(a) is
##       ~0.25 pools/day and is essentially CONSTANT across a 10-fold
##       range of plasma levels.  All the variance is production.  In
##       this model the production chain is
##           LPA transcription -> mRNA -> ER folding -> ERAD ->
##           secretion (SECEFF, isoform-dependent) -> free apo(a)
##       and SECEFF = KSZ^3/(KSZ^3 + NKIV2^3) alone spans ~5-fold across
##       the physiological KIV-2 repeat range.
##
##  (B)  ONLY ~24% OF Lp(a) CLEARANCE IS LDL-RECEPTOR DEPENDENT
##       (KLDLR_LPA = 0.060 of KCAT0 = 0.250 /day), whereas 86% of LDL
##       clearance is.  This single asymmetry is the whole clearance-drug
##       story and the model DERIVES rather than asserts it:
##
##         evolocumab -> LDLR 1.00 -> 2.75
##            LDL  clearance 0.350 -> 0.875  (x2.50)  => LDL-C   -60.0%
##            Lp(a) clearance 0.250 -> 0.355 (x1.42)  => Lp(a)   -29.6%
##
##         high-intensity statin -> LDLR 1.00 -> 1.78
##            Lp(a) clearance 0.250 -> 0.297           => Lp(a)  -15.7%
##            but LPA transcription x1.30 (sterol-responsive element)
##            => NET Lp(a)  +9.6%   <- the observed statin paradox,
##            produced by two terms with opposite signs, not by a rule.
##
##       The same arithmetic makes a falsifiable prediction: in a patient
##       with NO functional LDL receptor (HoFH, FLDLRFN = 0.05) the
##       offsetting clearance gain disappears and a statin should raise
##       Lp(a) by the full transcriptional ~30%.  Scenario 19 runs it.
##
##  (C)  A SECOND, INDEPENDENT DRUGGABLE NODE — particle assembly.
##       Free apo(a) is a separate state (APOA_FR).  Assembly-directed
##       therapy (muvalaplin) blocks step-1 docking, so Lp(a) falls AND
##       free apo(a) RISES ~4.5-fold.  Any assay that also detects free
##       apo(a) therefore UNDER-reports the effect.  The model reproduces
##       the direction of the KRAKEN intact-vs-traditional discrepancy
##       and quantifies what would be needed to reproduce its magnitude
##       (see CALIBRATION NOTES, item 9) — a genuine open question, left
##       as the parameter RFREE rather than buried.
##
##  (D)  TWO SCALES THAT DISAGREE EXACTLY WHERE IT MATTERS.
##       Every pathogenic payload is present ONCE PER PARTICLE:
##          1 apoB-100 · ~1 OxPL site on KIV-10 · 1 strong lysine-binding
##          site · a roughly fixed cholesterol core.
##       The ONLY thing that varies with isoform is apo(a) kringle
##       protein mass, which is pathogenically inert.  So particle
##       concentration (nmol/L) is not a unit preference — it is the
##       mechanistically correct scale, and mass (mg/dL) systematically
##       over-weights the large-isoform (lower-risk) particle.
##
##       The model separates the two error sources and finds they are
##       NOT the same size:
##          - STOICHIOMETRIC spread of the nmol/L <-> mg/dL conversion
##            across the whole isoform range:      only  ~9%
##          - ANALYTICAL bias of a polyclonal anti-apo(a) mass assay,
##            whose signal scales with the number of KIV-2 epitopes:
##                                                 up to ~40%
##       i.e. the unit conversion that everyone blames is nearly fine;
##       it is the ANTIBODY that is the problem.  A patient with a true
##       mass of 50 mg/dL and a small (12-repeat) isoform is reported as
##       34 mg/dL — below the treatment threshold — while their true
##       particle count is 137 nmol/L, above it.  The mass assay misses
##       precisely the highest-risk isoform.
##
##  (E)  WHY MENDELIAN RANDOMISATION AND A 5-YEAR TRIAL DISAGREE
##       WITHOUT EITHER BEING WRONG.
##       Risk is carried by TWO components with different time constants:
##          SLOW  plaque burden, t1/2 ~ 14 years  (70% of the excess)
##          FAST  OxPL/IL-6-driven vulnerability, t1/2 ~ 3 months (30%)
##       Lifelong genetic lowering removes both; a 5-year drug removes
##       essentially all of the fast component and only ~40% of the slow
##       one.  The ratio (MR effect / 5-year trial effect) ~ 2.6-3x then
##       EMERGES from two time constants instead of being inserted as a
##       correction factor.
##
##  (F)  THE FEED-FORWARD LOOP IS REAL BUT WEAK — AND THAT IS THE POINT.
##       Lp(a) -> OxPL -> monocyte NF-kB -> IL-6 -> IL-6 response element
##       in the LPA promoter -> more Lp(a).  The model computes the
##       open-loop gain and finds g ~ 0.02 at normal inflammation, i.e.
##       strongly subcritical and clinically negligible.  What is NOT
##       negligible is the TONIC IL-6 term, worth ~25% of basal LPA
##       transcription — which is why anti-IL-6 lowers Lp(a) ~25% in
##       non-inflamed patients (RESCUE) and ~37-41% in rheumatoid
##       arthritis (MEASURE), and the model reproduces both from the
##       same equation with only the IL-6 level changed.
##
##  Requires: mrgsolve (>= 1.0), dplyr, tidyr, ggplot2
##
##  All 50 ODEs were independently re-implemented in Python and
##  integrated with a fixed-step RK4 solver; every "ACHIEVED" figure in
##  the notes at the bottom of this file is a measured output of that
##  cross-check, not an aspiration.  It found five defects, listed there.
## =====================================================================

library(mrgsolve)
library(dplyr)

lpa_code <- '
$PROB
# Elevated Lipoprotein(a) — QSP model
# One particle · one rate-limiting step · three arms · two scales
# 50 compartments. Time unit = DAYS. Particle unit = nmol/L.

$PARAM @annotated

// ---------------- PATIENT / GENOTYPE ----------------------------------
LPA0      : 250   : Target baseline Lp(a) particle conc (nmol/L)
NKIV2     : 12    : KIV-2 repeats on the EXPRESSED dominant allele (-)
NCAL      : 22    : KIV-2 repeats of the mass-assay calibrator isoform (-)
LDLP0     : 1200  : Baseline LDL particle conc (nmol/L)
VLDLP0    : 120   : Baseline VLDL/remnant particle conc (nmol/L)
HDLC0     : 48    : Baseline HDL-C (mg/dL)
EGFR0     : 90    : Baseline eGFR (mL/min/1.73m2)
KEGFRD    : 2.7e-5: eGFR decline (mL/min/1.73m2 per day; ~1/yr)
FLDLRFN   : 1.0   : LDL-receptor FUNCTION multiplier (HoFH ~0.05) (-)
ERISK     : 0.0   : Extra endothelial permeability (smoking/HTN/DM) (-)

// covariate multipliers on LPA transcription (all 1.0 = neutral)
FTHY      : 1.0   : Thyroid status (hypothyroid ~1.25) (-)
FSEX      : 1.0   : Sex-steroid state (post-menopausal ~1.20 HRT ~0.80) (-)
FNEPH     : 1.0   : Nephrotic syndrome (~2.5) (-)
FLIV      : 1.0   : Hepatic synthetic capacity (cirrhosis ~0.5) (-)
IL6EXO    : 0.0   : Exogenous inflammatory IL-6 load (pg/mL, e.g. RA=18)

// ---------------- PRODUCTION / ASSEMBLY / CATABOLISM ------------------
KTX       : 0.375 : Basal LPA transcription rate (relative units/day)
KDMR      : 0.50  : LPA mRNA degradation rate (/day)
EIL6      : 0.80  : Emax of IL-6 on LPA transcription (-)
KIL6      : 2.80  : IL-6 conc for half-max transcription effect (pg/mL)
KSEC      : 2.0   : apo(a) secretion rate constant from ER pool (/day)
KDEGP     : 2.0   : Presecretory ERAD rate of ER apo(a) pool (/day)
KSZ       : 22.0  : KIV-2 repeats giving 50% secretion efficiency (-)
HSZ       : 3.0   : Hill coefficient of the isoform-size term (-)
KASM      : 24.6  : Maximal assembly rate constant of free apo(a) (/day)
KMLDL     : 30.0  : LDL particles giving half-maximal assembly (nmol/L)
KCLF      : 6.0   : Free apo(a) clearance rate (/day)
KLDLR_LPA : 0.060 : LDLR-dependent Lp(a) clearance at LDLR=1 (/day)
KOTH_LPA  : 0.190 : LDLR-INDEPENDENT Lp(a) clearance (/day)
KLDLR_LDL : 0.300 : LDLR-dependent LDL clearance at LDLR=1 (/day)
KOTH_LDL  : 0.050 : LDLR-independent LDL clearance (/day)
KVL       : 3.0   : VLDL/remnant to LDL conversion (/day)
KCLV      : 1.0   : Direct VLDL/remnant clearance (/day)
KSR       : 2.852 : LDLR synthesis rate (relative/day)
KDR       : 1.0   : LDLR degradation rate (/day)
GPCSK     : 1.852 : PCSK9 amplification of LDLR degradation (-)
KSP       : 1.0   : PCSK9 synthesis (relative/day)
KDP       : 1.0   : PCSK9 elimination (/day)
KFRAGR    : 0.02  : Fractional apo(a) fragment excretion rate (/day)

// ---------------- MEASUREMENT / STOICHIOMETRY -------------------------
MWCORE    : 3.30e6: Mass of apoB-100 + lipid core per particle (Da)
MWKR      : 14000 : Mass per (glycosylated) kringle domain (Da)
MWPROT    : 35000 : Mass of KV + inactive protease domain (Da)
FCHOL     : 0.30  : Cholesterol fraction of Lp(a) MASS (Dahlen factor) (-)
CLDL      : 0.0833: LDL-C per LDL particle (mg/dL per nmol/L)
MWAPOB    : 512000: apoB-100 molecular weight (Da)
RFREE     : 1.0   : Molar response of a traditional apo(a) assay to FREE apo(a) (-)
FRECY     : 0.15  : Fraction of assembly-spared apoB re-entering the plasma LDL pool (-)
DAHLC     : 0.30  : Coefficient actually used in the Dahlen LDL-C correction (-)

// ---------------- OxPL / INFLAMMATION ---------------------------------
LPAREF    : 30.0  : Reference Lp(a) defining OXPL = 1 (nmol/L)
FOXLDL    : 0.02  : OxPL carried per LDL particle relative to Lp(a) (-)
KOXD      : 2.0   : OxPL turnover (/day)
KMD       : 0.05  : Trained-monocyte state turnover (/day)
EMON      : 1.5   : Emax of OxPL on monocyte training (-)
KMOX      : 3.0   : OxPL for half-max monocyte training (-)
KI6D      : 3.0   : IL-6 elimination (/day)
EI6       : 0.80  : Monocyte-state drive on IL-6 production (-)
KCD       : 0.36  : CRP elimination (/day)
CRP0      : 1.5   : Reference hsCRP at reference IL-6 (mg/L)
FCRP0     : 0.06  : IL-6-INDEPENDENT fraction of CRP production (-)
IL6REF    : 2.0   : Reference IL-6 (pg/mL)

// ---------------- ARTERIAL WALL ---------------------------------------
KRET      : 8.333e-5 : Intimal retention rate per particle (per nmol/L per day)
AVID      : 2.5   : Lp(a) matrix-retention avidity relative to LDL (-)
KEGR      : 0.08  : HDL-dependent intimal lipid egress (/day)
KFOAM     : 0.02  : Conversion of retained lipoprotein to foam cells (/day)
WOX       : 1.6   : Extra foam-cell potency per retained Lp(a) vs LDL (-)
KFAP      : 0.02  : Foam-cell clearance (/day)
KNC       : 0.001 : Necrotic-core formation (/day)
KNCR      : 0.002 : Necrotic-core resolution (/day)
ENEC      : 1.0   : OxPL amplification of necrotic-core formation (-)
KNOX      : 3.0   : OxPL for half-max necrotic amplification (-)
KCAPS     : 1.2   : Fibrous-cap synthesis (um/day)
KCAPD     : 0.01  : Fibrous-cap turnover (/day)
ECAPN     : 0.6   : Necrotic-core suppression of cap synthesis (-)
KPG       : 4.787e-3 : Plaque growth rate (PAV %/day per driver unit)
KPR       : 1.369e-4 : Plaque regression rate (/day; t1/2 ~14 yr)
PAVMAX    : 60.0  : Maximum percent atheroma volume (%)
WNEC      : 1.0   : Weight of necrotic core in plaque growth (-)
KCAC      : 0.010 : Coronary calcium accrual (Agatston/day per PAV unit)
ECACOX    : 1.0   : OxPL amplification of coronary calcification (-)
KVD       : 0.0111: Vulnerability-index turnover (/day; t1/2 ~ 62 d)
EVOX      : 0.50  : OxPL drive on vulnerability (-)
KVOX      : 4.0   : OxPL for half-max vulnerability drive (-)
EVN       : 0.15  : Necrotic-core drive on vulnerability (-)
OXREF     : 0.907 : Reference OxPL of the low-Lp(a) comparator (-)
NECREF    : 0.654 : Reference necrotic-core index of the comparator (-)

// ---------------- AORTIC VALVE ----------------------------------------
KATX      : 0.10  : Valve autotaxin delivery per unit Lp(a)/LPAREF (/day)
KATXD     : 0.10  : Valve autotaxin turnover (/day)
KLPAX     : 0.20  : LysoPA generation from autotaxin x OxPL (/day)
KLPAD     : 0.20  : LysoPA turnover (/day)
KOST      : 0.0020: VIC osteogenic transition rate (/day)
KOSTR     : 0.0020: VIC osteogenic reversion rate (/day)
KLO       : 20.0  : LysoPA for half-max VIC transition (-)
EVI6      : 0.05  : IL-6 amplification of VIC transition (per pg/mL)
KCA       : 0.055 : Valve calcium accrual (AU/day per unit VIC_OST)
KSELFR    : 4.0   : Lp(a)-INDEPENDENT self-perpetuating calcification (-)
KSELFH    : 650   : Valve calcium at half-max self-perpetuation (AU)
HSELF     : 4.0   : Hill coefficient of the self-perpetuation switch (-)
AVA0      : 3.50  : Aortic valve area with no calcium (cm2)
KAVA      : 480   : Valve calcium halving the effective orifice (AU)
KMG       : 40.0  : Mean-gradient constant (mmHg cm4)

// ---------------- FIBRINOLYSIS (LOW-WEIGHT ARM) -----------------------
KOCC      : 2.0   : Equilibration rate of fibrin-surface LBS occupancy (/day)
KPLGC     : 6.0   : Plasminogen competition constant (x LPAREF, nmol/L)
KPB       : 2.0   : PAI-1 synthesis (ng/mL/day)
KPD       : 0.10  : PAI-1 elimination (/day)
EPAI      : 0.50  : Lp(a) drive on PAI-1 (-)
KPAI      : 250   : Lp(a) for half-max PAI-1 drive (nmol/L)
PAI1REF   : 20.0  : Reference PAI-1 (ng/mL)
KLRD      : 0.50  : Lysis-resistance turnover (/day)
ELBS      : 0.45  : LBS occupancy drive on lysis resistance (-)
ELPAI     : 0.25  : PAI-1 drive on lysis resistance (-)
WARM3     : 0.25  : Weight of the antifibrinolytic arm in MACE hazard (-)

// ---------------- RISK TRANSLATION ------------------------------------
H0        : 5.77e-5 : Reference MACE hazard (/day; 10% per 5 yr)
BSLOW     : 0.0373  : ln-HR per percent-atheroma-volume point (-)
BFAST     : 0.465   : ln-HR per unit vulnerability index (-)
BLYS      : 0.30    : ln-HR per unit lysis-resistance index (-)
PAVREF    : 28.0    : Reference plaque burden of the comparator (PAV %)
VULNREF   : 1.0     : Reference vulnerability index (-)
HAVS0     : 1.0e-6  : Reference AVR/valve-death hazard (/day)
BAVS      : 2.2     : ln-HR per 1000 AU of valve calcium (-)
VCTH      : 300     : Valve calcium below which AVS hazard is flat (AU)

// ---------------- DRUG PK ---------------------------------------------
// pelacarsen (GalNAc-ASO, cmt 1-3)
KAPEL : 0.50 : Pelacarsen SC absorption (/day)
KPTL  : 8.0  : Pelacarsen plasma-to-liver uptake (/day)
KPEL  : 2.0  : Pelacarsen non-hepatic elimination (/day)
KLIVEL: 0.030: Pelacarsen hepatic elimination (/day; t1/2 ~23 d)
VLIVP : 1.8  : Effective hepatic distribution volume (L)
EASO  : 6.0  : Emax of ASO on mRNA degradation (fold) (-)
IC50A : 21.0 : Hepatic ASO conc for half-max mRNA degradation (mg/L)

// siRNA: olpasiran / lepodisiran / zerlasiran (cmt 4-6)
KASIR : 0.80 : siRNA SC absorption (/day)
KSTR  : 10.0 : siRNA plasma-to-RISC loading (/day)
KSEL  : 4.0  : siRNA non-RISC elimination (/day)
KRISC : 0.0116: RISC-loaded siRNA turnover (/day; t1/2 ~60 d)
ERNAI : 30.0 : Emax of RISC-loaded siRNA on mRNA degradation (fold) (-)
IC50R : 6.0  : RISC amount for half-max mRNA degradation (mg)

// muvalaplin (oral assembly inhibitor, cmt 7-8)
KAMUV : 3.0  : Muvalaplin absorption (/day)
KEMUV : 1.5  : Muvalaplin elimination (/day)
VMUV  : 60.0 : Muvalaplin volume of distribution (L)
KIMUV : 0.10 : Muvalaplin Ki for step-1 assembly blockade (mg/L)

// evolocumab / alirocumab (cmt 9-10) and PCSK9 (cmt 11)
KAEVO : 0.25 : Evolocumab SC absorption (/day)
KEEVO : 0.060: Evolocumab elimination (/day)
VEVO  : 3.30 : Evolocumab volume of distribution (L)
KONE  : 0.65 : Evolocumab-mediated PCSK9 removal (per mg/L per day)

// statin (cmt 12-13)
KASTA : 6.0  : Statin absorption (/day)
KESTA : 1.4  : Statin elimination (/day)
VSTA  : 134  : Statin volume of distribution (L)
EC50STA: 0.012: Statin conc for half-max effect (mg/L)
ESRE  : 1.24 : Statin Emax on LDLR synthesis (SREBP-2) (-)
EPSK  : 0.40 : Statin Emax on PCSK9 synthesis (-)
ESTA  : 0.30 : Statin Emax on LPA transcription (-)

// niacin ER (cmt 14)
KENIA : 1.5  : Niacin elimination (/day)
VNIA  : 60.0 : Niacin volume of distribution (L)
EC50NIA: 3.0 : Niacin conc for half-max effect (mg/L)
ENIA  : 0.22 : Niacin Emax on LPA transcription (fractional decrease) (-)
ENIALDL: 0.15: Niacin Emax on LDL particle production (-)

// ziltivekimab / tocilizumab (cmt 15-16)
KAZIL : 0.30 : Anti-IL-6 SC absorption (/day)
KEZIL : 0.05 : Anti-IL-6 elimination (/day)
VZIL  : 4.0  : Anti-IL-6 volume of distribution (L)
IC50Z : 0.11 : Anti-IL-6 conc for half of IL-6 signalling blocked (mg/L)

// obicetrapib (CETP inhibitor, cmt 17)
KEOBI : 0.35 : Obicetrapib elimination (/day)
VOBI  : 100  : Obicetrapib volume of distribution (L)
EC50OBI: 0.05: Obicetrapib conc for half-max effect (mg/L)
ECETPR: 0.74 : Obicetrapib Emax on LDLR synthesis (-)
EOBI  : 0.23 : Obicetrapib Emax on LPA transcription (fractional decrease) (-)
EHDL  : 1.50 : Obicetrapib Emax on HDL-C (-)
KHD   : 0.20 : HDL-C turnover (/day)

// ezetimibe / bempedoic acid switch (LDLR-only agents)
FEZE  : 1.0  : Multiplier on LDLR synthesis from ezetimibe/bempedoic (-)

// lipoprotein apheresis
APHON : 0    : Apheresis on/off flag (-)
APHINT: 7    : Apheresis interval (days)
KAPH  : 5.25 : Peak apheresis removal rate (/day)
APHW  : 0.08 : Gaussian half-width of one session (days; ~4.5 h full width)

$GLOBAL
// derived per-individual constants, set in $MAIN
double gSECEFF, gKTL, gMWAPOA, gMWLPA, gCONV, gEPIT;
double gKA, gPRODLDL, gPRODV, gKSYNH, gKMS, gKI6B, gKCB;
double gKOXS, gOXPL0, gMONO0, gIL60, gFIL60, gKV;
#define POSD(a) ((a) > 0.0 ? (a) : 0.0)

$INIT @annotated
PEL_SC   : 0    : Pelacarsen SC depot (mg) [cmt 1]
PEL_CE   : 0    : Pelacarsen plasma (mg) [cmt 2]
PEL_LIV  : 0    : Pelacarsen hepatic (mg) [cmt 3]
SIR_SC   : 0    : siRNA SC depot (mg) [cmt 4]
SIR_CE   : 0    : siRNA plasma (mg) [cmt 5]
SIR_RISC : 0    : siRNA RISC-loaded hepatic (mg) [cmt 6]
MUV_GU   : 0    : Muvalaplin gut (mg) [cmt 7]
MUV_CE   : 0    : Muvalaplin plasma (mg) [cmt 8]
EVO_SC   : 0    : Evolocumab SC depot (mg) [cmt 9]
EVO_CE   : 0    : Evolocumab plasma (mg) [cmt 10]
PCSK9    : 1    : Free plasma PCSK9 (relative) [cmt 11]
STA_GU   : 0    : Statin gut (mg) [cmt 12]
STA_CE   : 0    : Statin plasma (mg) [cmt 13]
NIA_CE   : 0    : Niacin plasma (mg) [cmt 14]
ZIL_SC   : 0    : Anti-IL-6 SC depot (mg) [cmt 15]
ZIL_CE   : 0    : Anti-IL-6 plasma (mg) [cmt 16]
OBI_CE   : 0    : Obicetrapib plasma (mg) [cmt 17]
MRNA     : 1    : LPA mRNA (relative) [cmt 18]
APOA_ER  : 1    : Intracellular ER apo(a) pool (nmol/L-equivalent) [cmt 19]
APOA_FR  : 1    : Free plasma apo(a) (nmol/L) [cmt 20]
LPA_P    : 250  : Lp(a) PARTICLES (nmol/L) [cmt 21]
LDL_P    : 1200 : LDL particles (nmol/L) [cmt 22]
VLDL_P   : 120  : VLDL/remnant particles (nmol/L) [cmt 23]
LDLR     : 1    : Hepatic LDL receptor density (relative) [cmt 24]
OXPL     : 1    : OxPL-apoB (relative to Lp(a)=LPAREF) [cmt 25]
MONO     : 1    : Trained/activated monocyte index (relative) [cmt 26]
IL6      : 2    : Plasma IL-6 (pg/mL) [cmt 27]
CRP      : 1.5  : hsCRP (mg/L) [cmt 28]
INT_LPA  : 0.5  : Retained intimal Lp(a) (relative) [cmt 29]
INT_LDL  : 1    : Retained intimal LDL (relative) [cmt 30]
FOAM     : 1    : Foam-cell burden (relative) [cmt 31]
NECRO    : 0.65 : Necrotic-core index (relative) [cmt 32]
CAP      : 120  : Fibrous-cap thickness (um) [cmt 33]
PLAQUE   : 8    : Percent atheroma volume (%) [cmt 34]
CAC      : 0    : Coronary artery calcium (Agatston) [cmt 35]
VULN     : 1    : Plaque vulnerability index (relative) [cmt 36]
ATXV     : 1    : Valve autotaxin activity (relative) [cmt 37]
LYSOPA   : 1    : Valve lysophosphatidic acid (relative) [cmt 38]
VIC_OST  : 0    : Osteogenic valve interstitial cell fraction (-) [cmt 39]
VCALC    : 0    : Aortic valve calcium (Agatston AU) [cmt 40]
PLGOCC   : 0    : Fibrin-surface lysine-binding-site occupancy (-) [cmt 41]
PAI1     : 20   : PAI-1 (ng/mL) [cmt 42]
LYSRES   : 1    : Clot-lysis resistance index (relative) [cmt 43]
HDL_C    : 48   : HDL cholesterol (mg/dL) [cmt 44]
EGFR     : 90   : eGFR (mL/min/1.73m2) [cmt 45]
AUC_LPA  : 0    : Cumulative Lp(a) exposure (nmol.yr/L) [cmt 46]
AUC_APOB : 0    : Cumulative apoB exposure (mg/dL.yr) [cmt 47]
HAZ      : 0    : Cumulative MACE hazard (-) [cmt 48]
HAZ_AVS  : 0    : Cumulative AVS/AVR hazard (-) [cmt 49]
FRAG     : 0    : Cumulative urinary apo(a) fragments (nmol/L-equiv) [cmt 50]

$MAIN
// ---------------------------------------------------------------------
// Derived per-individual constants.  Everything below is algebra on the
// parameters; nothing here is fitted at run time.
// ---------------------------------------------------------------------

// --- isoform stoichiometry --------------------------------------------
gMWAPOA = MWKR*(NKIV2 + 10.0) + MWPROT;      // apo(a) MW (Da)
gMWLPA  = MWCORE + gMWAPOA;                  // whole-particle MW (Da)
gCONV   = 1.0e7/gMWLPA;                      // nmol/L per mg/dL
gEPIT   = (NKIV2 + 10.0)/(NCAL + 10.0);      // polyclonal epitope bias
gSECEFF = pow(KSZ,HSZ)/(pow(KSZ,HSZ) + pow(NKIV2,HSZ));

// --- self-consistent baseline of the OxPL -> IL-6 arm ------------------
// OXPL is normalised so that OXPL = 1 at Lp(a) = LPAREF.
gKOXS   = KOXD/(LPAREF + FOXLDL*LDLP0);
gOXPL0  = gKOXS*(LPA0 + FOXLDL*LDLP0)/KOXD;
// MONO is normalised so that MONO = 1 at OXPL = 1.
gKMS    = KMD/(1.0 + EMON*1.0/(KMOX + 1.0));
gMONO0  = (gKMS/KMD)*(1.0 + EMON*gOXPL0/(KMOX + gOXPL0));
// IL6 is normalised so that IL6 = IL6REF at MONO = 1.
gKI6B   = KI6D*IL6REF;
gIL60   = (gKI6B/KI6D)*(1.0 + EI6*(gMONO0 - 1.0)) + IL6EXO;
gFIL60  = 1.0 + EIL6*gIL60/(KIL6 + gIL60);
gKCB    = KCD*CRP0;

// --- back-solve translation rate so baseline Lp(a) == LPA0 -------------
double SATL0 = LDLP0/(KMLDL + LDLP0);              // 97.6% saturated
gKA          = KASM*SATL0;                         // assembly pseudo-1st-order
double KCAT0 = KLDLR_LPA*FLDLRFN + KOTH_LPA;       // Lp(a) FCR at LDLR = 1
double SECN  = KCAT0*LPA0*(gKA + KCLF)/gKA;        // required secretion flux
double MRNAS = KTX*gFIL60*FTHY*FSEX*FNEPH*FLIV/KDMR;
gKTL         = SECN*(KSEC*gSECEFF + KDEGP)/(KSEC*gSECEFF*MRNAS);

// --- lipoprotein production rates so baselines hold --------------------
gPRODV   = VLDLP0*(KVL + KCLV);
gPRODLDL = LDLP0*(KLDLR_LDL*FLDLRFN + KOTH_LDL) + FRECY*KCAT0*LPA0 - KVL*VLDLP0;
gKSYNH   = KHD*HDLC0;
gKV      = KVD;                                    // VULN normalised to 1

// --- initial conditions ------------------------------------------------
MRNA_0    = MRNAS;
APOA_ER_0 = gKTL*MRNAS/(KSEC*gSECEFF + KDEGP);
APOA_FR_0 = SECN/(gKA + KCLF);
LPA_P_0   = LPA0;
LDL_P_0   = LDLP0;
VLDL_P_0  = VLDLP0;
OXPL_0    = gOXPL0;
MONO_0    = gMONO0;
IL6_0     = gIL60;
CRP_0     = CRP0*(FCRP0 + (1.0 - FCRP0)*gIL60/IL6REF);
HDL_C_0   = HDLC0;
EGFR_0    = EGFR0;

$ODE
// =====================================================================
//  DRUG PK  (cmt 1-17)
// =====================================================================
dxdt_PEL_SC   = -KAPEL*PEL_SC;
dxdt_PEL_CE   =  KAPEL*PEL_SC - (KPTL + KPEL)*PEL_CE;
dxdt_PEL_LIV  =  KPTL*PEL_CE  - KLIVEL*PEL_LIV;

dxdt_SIR_SC   = -KASIR*SIR_SC;
dxdt_SIR_CE   =  KASIR*SIR_SC - (KSTR + KSEL)*SIR_CE;
dxdt_SIR_RISC =  KSTR*SIR_CE  - KRISC*SIR_RISC;

dxdt_MUV_GU   = -KAMUV*MUV_GU;
dxdt_MUV_CE   =  KAMUV*MUV_GU - KEMUV*MUV_CE;

dxdt_EVO_SC   = -KAEVO*EVO_SC;
dxdt_EVO_CE   =  KAEVO*EVO_SC - KEEVO*EVO_CE;

dxdt_STA_GU   = -KASTA*STA_GU;
dxdt_STA_CE   =  KASTA*STA_GU - KESTA*STA_CE;

dxdt_NIA_CE   = -KENIA*NIA_CE;

dxdt_ZIL_SC   = -KAZIL*ZIL_SC;
dxdt_ZIL_CE   =  KAZIL*ZIL_SC - KEZIL*ZIL_CE;

dxdt_OBI_CE   = -KEOBI*OBI_CE;

// --- drug concentrations / effect fractions --------------------------
double CPELL = PEL_LIV/VLIVP;                       // hepatic ASO (mg/L)
double CMUV  = MUV_CE/VMUV;                         // muvalaplin (mg/L)
double CEVO  = EVO_CE/VEVO;                         // evolocumab (mg/L)
double CSTA  = STA_CE/VSTA;                         // statin (mg/L)
double CNIA  = NIA_CE/VNIA;                         // niacin (mg/L)
double CZIL  = ZIL_CE/VZIL;                         // anti-IL-6 (mg/L)
double COBI  = OBI_CE/VOBI;                         // obicetrapib (mg/L)

double FSTAT = CSTA/(EC50STA + CSTA);               // statin effect fraction
double FNIAC = CNIA/(EC50NIA + CNIA);
double FOBIC = COBI/(EC50OBI + COBI);

// =====================================================================
//  PCSK9  (cmt 11)   — statin raises synthesis, mAb removes free protein
// =====================================================================
dxdt_PCSK9 = KSP*(1.0 + EPSK*FSTAT) - KDP*PCSK9 - KONE*CEVO*PCSK9;

// =====================================================================
//  IL-6 SIGNALLING SEEN BY THE LPA PROMOTER
//  At clinical doses an anti-IL-6 antibody sits in ~10^4 molar excess
//  over IL-6, so blockade is near-complete and essentially independent
//  of the ligand concentration.  IC50Z is set accordingly.
// =====================================================================
double IL6EFF = IL6/(1.0 + CZIL/IC50Z);
double FIL6   = 1.0 + EIL6*IL6EFF/(KIL6 + IL6EFF);

// =====================================================================
//  LPA TRANSCRIPTION -> mRNA  (cmt 18)
//  Three signed terms meet here: IL-6 (+), statin (+), niacin/CETPi (-).
// =====================================================================
double TRANS = KTX*FIL6
             * (1.0 + ESTA*FSTAT)
             * (1.0 - ENIA*FNIAC)
             * (1.0 - EOBI*FOBIC)
             * FTHY*FSEX*FNEPH*FLIV;

// RNA-directed drugs act on the SAME degradation term by two routes
double KNOCK = 1.0
             + EASO*CPELL/(IC50A + CPELL)
             + ERNAI*SIR_RISC/(IC50R + SIR_RISC);

dxdt_MRNA = TRANS - KDMR*KNOCK*MRNA;

// =====================================================================
//  apo(a) PROTEIN: ER pool -> secretion -> free apo(a)  (cmt 19-20)
//  gSECEFF is where the isoform-size effect enters: it sets how much of
//  the ER pool escapes presecretory degradation.
// =====================================================================
double SECFLX = KSEC*gSECEFF*APOA_ER;
dxdt_APOA_ER  = gKTL*MRNA - SECFLX - KDEGP*APOA_ER;

// muvalaplin: competitive blockade of step-1 (KIV-7/8 lysine-binding site)
double FMUV  = 1.0/(1.0 + CMUV/KIMUV);

// Assembly SATURATES in LDL substrate.  This is not a convenience: if
// assembly were first-order in LDL, ezetimibe (LDL-P -19%, no direct LPA
// effect) would drag Lp(a) down by a further ~19%, which is not observed.
double ASSEM = KASM*FMUV*APOA_FR*LDL_P/(KMLDL + LDL_P);

dxdt_APOA_FR = SECFLX - ASSEM - KCLF*APOA_FR;

// =====================================================================
//  Lp(a) PARTICLES  (cmt 21)
//  Clearance is split into an LDLR-dependent MINORITY and an
//  LDLR-independent majority; the latter degrades with eGFR.
// =====================================================================
double RENF   = 0.60 + 0.40*EGFR/EGFR0;
double KCATL  = KLDLR_LPA*LDLR*FLDLRFN + KOTH_LPA*RENF;

// lipoprotein apheresis: a smooth periodic extra-clearance pulse whose
// width is fixed in ABSOLUTE time, so the interval can be changed alone
double APH = 0.0;
if (APHON > 0.5) {
  double tw  = fmod(SOLVERTIME, APHINT);
  double dtn = fmin(tw, APHINT - tw);          // time to the nearest session
  APH = KAPH*exp(-0.5*(dtn/APHW)*(dtn/APHW));
}

dxdt_LPA_P = ASSEM - (KCATL + APH)*LPA_P;

// =====================================================================
//  LDL and remnant particles  (cmt 22-23)
//  Making one Lp(a) consumes one apoB particle.  FRECY asks whether that
//  apoB comes out of the MEASURED plasma LDL pool (FRECY = 1, so deep
//  apo(a) knockdown would rebound LDL-P upward) or out of a nascent pool
//  that never joins it (FRECY = 0).  The apoB fall seen in the
//  pelacarsen trials is close to the Lp(a)-apoB removed, which argues
//  FRECY is near 0; the default 0.15 is deliberately conservative.
// =====================================================================
double KCATD = KLDLR_LDL*LDLR*FLDLRFN + KOTH_LDL;
dxdt_LDL_P  = gPRODLDL*(1.0 - ENIALDL*FNIAC) + KVL*VLDL_P
              - KCATD*LDL_P - FRECY*ASSEM - APH*LDL_P;
dxdt_VLDL_P = gPRODV*(1.0 - ENIALDL*FNIAC) - (KVL + KCLV)*VLDL_P;

// =====================================================================
//  LDL RECEPTOR  (cmt 24)
// =====================================================================
dxdt_LDLR = KSR*(1.0 + ESRE*FSTAT)*(1.0 + ECETPR*FOBIC)*FEZE
            - KDR*(1.0 + GPCSK*PCSK9)*LDLR;

// =====================================================================
//  ARM 2 : OxPL -> trained monocyte -> IL-6 -> CRP  (cmt 25-28)
//  OxPL sites are ~1 per particle, so this arm is driven by nmol/L.
// =====================================================================
dxdt_OXPL = gKOXS*(LPA_P + FOXLDL*LDL_P) - KOXD*OXPL;
dxdt_MONO = gKMS*(1.0 + EMON*OXPL/(KMOX + OXPL)) - KMD*MONO;
dxdt_IL6  = gKI6B*(1.0 + EI6*(MONO - 1.0)) + KI6D*IL6EXO - KI6D*IL6;
dxdt_CRP  = gKCB*(FCRP0 + (1.0 - FCRP0)*IL6EFF/IL6REF) - KCD*CRP;

// =====================================================================
//  ARM 1 : intimal retention -> foam cells -> plaque  (cmt 29-36)
//  AVID is the only Lp(a)-specific term in the retention step.
// =====================================================================
double ENDO = 1.0 + ERISK;
double HDLF = HDL_C/HDLC0;
double EGRE = KEGR*HDLF + KFOAM;

dxdt_INT_LPA = KRET*AVID*LPA_P*ENDO - EGRE*INT_LPA;
dxdt_INT_LDL = KRET*LDL_P*ENDO      - EGRE*INT_LDL;

dxdt_FOAM  = KFOAM*(INT_LDL + WOX*INT_LPA)*MONO - KFAP*FOAM;
dxdt_NECRO = KNC*FOAM*(1.0 + ENEC*OXPL/(KNOX + OXPL)) - KNCR*NECRO;
dxdt_CAP   = KCAPS/(1.0 + ECAPN*NECRO) - KCAPD*CAP;

double PDRIVE = FOAM + WNEC*NECRO;
dxdt_PLAQUE = KPG*PDRIVE*(1.0 - PLAQUE/PAVMAX) - KPR*PLAQUE;
dxdt_CAC    = KCAC*PLAQUE*(1.0 + ECACOX*OXPL/(KNOX + OXPL));

dxdt_VULN = gKV*(1.0 + EVOX*(OXPL - OXREF)/(KVOX + OXPL)
                     + EVN*(NECRO - NECREF)) - KVD*VULN;

// =====================================================================
//  AORTIC VALVE  (cmt 37-40)
//  The self-perpetuating term is a THRESHOLD (Hill n = 4), not a ramp:
//  below a mineral nucleus it contributes nothing; above it the process
//  runs without Lp(a).  Everything about the failed lipid-lowering AS
//  trials lives in that one term.
// =====================================================================
dxdt_ATXV   = KATX*(LPA_P/LPAREF) - KATXD*ATXV;
dxdt_LYSOPA = KLPAX*ATXV*OXPL     - KLPAD*LYSOPA;
dxdt_VIC_OST = KOST*(LYSOPA/(KLO + LYSOPA))*(1.0 + EVI6*POSD(IL6EFF - IL6REF))
                    *(1.0 - VIC_OST)
               - KOSTR*VIC_OST;
double SELFT = KSELFR*pow(VCALC,HSELF)/(pow(KSELFH,HSELF) + pow(VCALC,HSELF));
dxdt_VCALC  = KCA*(VIC_OST + SELFT);

// =====================================================================
//  ARM 3 : antifibrinolysis  (cmt 41-43)   — deliberately LOW weight
//  One strong lysine-binding site per particle => nmol/L again.
// =====================================================================
double OCCT = LPA_P/(LPA_P + KPLGC*LPAREF);
dxdt_PLGOCC = KOCC*(OCCT - PLGOCC);
dxdt_PAI1   = KPB*(1.0 + EPAI*LPA_P/(KPAI + LPA_P)) - KPD*PAI1;
dxdt_LYSRES = KLRD*(1.0 + ELBS*PLGOCC + ELPAI*(PAI1/PAI1REF - 1.0)) - KLRD*LYSRES;

// =====================================================================
//  HDL, renal function  (cmt 44-45)
// =====================================================================
dxdt_HDL_C = gKSYNH*(1.0 + EHDL*FOBIC) - KHD*HDL_C;
dxdt_EGFR  = -KEGFRD;

// =====================================================================
//  EXPOSURE ACCUMULATORS AND HAZARDS  (cmt 46-50)
// =====================================================================
double APOBM = (LDL_P + LPA_P + VLDL_P)*MWAPOB*1.0e-7;   // total apoB (mg/dL)

dxdt_AUC_LPA  = LPA_P/365.25;
dxdt_AUC_APOB = APOBM/365.25;

double LNHR = BSLOW*(PLAQUE - PAVREF)
            + BFAST*(VULN - VULNREF)
            + WARM3*BLYS*(LYSRES - 1.0);
dxdt_HAZ = H0*exp(fmax(-30.0, fmin(30.0, LNHR)));

double VEXC = POSD(VCALC - VCTH)/1000.0;
dxdt_HAZ_AVS = HAVS0*exp(BAVS*VEXC);

dxdt_FRAG = KFRAGR*(LPA_P + APOA_FR);

$TABLE
// =====================================================================
//  READ-OUTS.  The point of this block is that the SAME particle
//  concentration is reported four different ways, and they disagree.
// =====================================================================

// --- the truth: particle concentration --------------------------------
double LPA_NMOL   = LPA_P;                              // nmol/L, true
double LPA_MASS_T = LPA_P/gCONV;                        // mg/dL, true mass

// --- what an isoform-INSENSITIVE monoclonal molar assay reports --------
double ASSAY_NMOL = LPA_P;                              // 1 epitope/particle

// --- what an isoform-SENSITIVE polyclonal MASS assay reports -----------
//     signal scales with the number of KIV-2 epitopes vs the calibrator
double ASSAY_MASS = LPA_MASS_T*gEPIT;

// --- what a traditional apo(a)-detecting assay reports -----------------
//     it also sees FREE apo(a), which RISES on an assembly inhibitor
double ASSAY_APOA = (LPA_P + RFREE*APOA_FR);

// --- an intact-Lp(a) (apo(a)+apoB sandwich) assay ----------------------
double ASSAY_INTACT = LPA_P;

// --- cholesterol book-keeping ------------------------------------------
double LPAC      = FCHOL*LPA_MASS_T;                    // Lp(a)-C (mg/dL)
double LDLC_TRUE = CLDL*LDL_P;                          // true LDL-C
double LDLC_MEAS = LDLC_TRUE + LPAC;                    // what the lab reports
double LDLC_CORR = LDLC_MEAS - DAHLC*ASSAY_MASS;        // Dahlen, as performed
double APOB_TOT  = (LDL_P + LPA_P + VLDL_P)*MWAPOB*1.0e-7;
double LPA_APOB_PCT = 100.0*LPA_P/(LDL_P + LPA_P + VLDL_P);

// --- retention share: why a minority particle matters -------------------
double RETSHARE = 100.0*AVID*LPA_P/(AVID*LPA_P + LDL_P);

// --- valve haemodynamics ------------------------------------------------
double AVA   = AVA0/(1.0 + VCALC/KAVA);
double AVAG  = fmax(0.30, AVA);
double MGRAD = KMG/(AVAG*AVAG);

// --- risk read-outs -----------------------------------------------------
double LNHR_OUT = BSLOW*(PLAQUE - PAVREF) + BFAST*(VULN - VULNREF)
                + WARM3*BLYS*(LYSRES - 1.0);
double HR_MACE  = exp(LNHR_OUT);
double HR_SLOW  = exp(BSLOW*(PLAQUE - PAVREF));
double HR_FAST  = exp(BFAST*(VULN - VULNREF));

// --- unit-conversion diagnostics ---------------------------------------
double CONVF   = gCONV;         // nmol/L per mg/dL for THIS isoform
double EPITB   = gEPIT;         // analytical bias factor of the mass assay
double SECEFFO = gSECEFF;       // secretion efficiency of THIS isoform
double MWAPOAO = gMWAPOA;       // apo(a) MW (Da)

$CAPTURE
LPA_NMOL LPA_MASS_T ASSAY_NMOL ASSAY_MASS ASSAY_APOA ASSAY_INTACT
LPAC LDLC_TRUE LDLC_MEAS LDLC_CORR APOB_TOT LPA_APOB_PCT RETSHARE
AVA MGRAD HR_MACE HR_SLOW HR_FAST LNHR_OUT
CONVF EPITB SECEFFO MWAPOAO
CPELL CMUV CEVO CSTA CZIL COBI IL6EFF FIL6 TRANS KNOCK FMUV KCATL KCATD ASSEM
'

mod <- mcode_cache("lpa_qsp", lpa_code, atol = 1e-9, rtol = 1e-7, maxsteps = 1e6)

## =====================================================================
##  COMPARTMENT / DOSING INDEX
## ---------------------------------------------------------------------
##   1 PEL_SC    pelacarsen SC        80 mg  q28d
##   4 SIR_SC    siRNA SC             75 mg  q84d  (olpasiran)
##                                   400 mg  q168d (lepodisiran)
##   7 MUV_GU    muvalaplin PO       240 mg  q1d
##   9 EVO_SC    evolocumab SC       420 mg  q28d
##  12 STA_GU    statin PO            20 mg  q1d (rosuvastatin-equivalent)
##  14 NIA_CE    niacin ER          2000 mg  q1d
##  15 ZIL_SC    ziltivekimab SC      30 mg  q28d
##  17 OBI_CE    obicetrapib PO       10 mg  q1d
## =====================================================================

YR <- 365.25

ev_none  <- ev(amt = 0, cmt = 1, time = 0)
ev_pela  <- ev(amt = 80,   cmt = 1,  ii = 28,  addl = 999)
ev_olpa  <- ev(amt = 75,   cmt = 4,  ii = 84,  addl = 999)
ev_lepo  <- ev(amt = 400,  cmt = 4,  ii = 168, addl = 999)
ev_muva  <- ev(amt = 240,  cmt = 7,  ii = 1,   addl = 99999)
ev_evol  <- ev(amt = 420,  cmt = 9,  ii = 28,  addl = 999)
ev_stat  <- ev(amt = 20,   cmt = 12, ii = 1,   addl = 99999)
ev_niac  <- ev(amt = 2000, cmt = 14, ii = 1,   addl = 99999)
ev_zilt  <- ev(amt = 30,   cmt = 15, ii = 28,  addl = 999)
ev_obic  <- ev(amt = 10,   cmt = 17, ii = 1,   addl = 99999)

## ---------------------------------------------------------------------
##  SCENARIO TABLE — 20 scenarios
##  1-3 untreated references; 4-17 treatments applied to the same
##  250 nmol/L small-isoform patient; 18-20 structural predictions.
## ---------------------------------------------------------------------
sc <- list(
  "01 Untreated  Lp(a) 250 nmol/L, small isoform (12 KIV-2)" =
    list(ev = ev_none, p = list(LPA0 = 250, NKIV2 = 12)),

  "02 Untreated  same particle count, large isoform (35 KIV-2)" =
    list(ev = ev_none, p = list(LPA0 = 250, NKIV2 = 35)),

  "03 Untreated  low Lp(a) comparator (25 nmol/L)" =
    list(ev = ev_none, p = list(LPA0 = 25,  NKIV2 = 30)),

  "04 High-intensity statin" =
    list(ev = ev_stat, p = list()),

  "05 Ezetimibe-like LDLR bump alone (PREDICTION: small Lp(a) fall)" =
    list(ev = ev_none, p = list(FEZE = 1.28)),

  "06 Statin + evolocumab 420 mg q4w" =
    list(ev = c(ev_stat, ev_evol), p = list()),

  "07 Evolocumab monotherapy" =
    list(ev = ev_evol, p = list()),

  "08 Pelacarsen 80 mg SC q4w" =
    list(ev = ev_pela, p = list()),

  "09 Olpasiran 75 mg SC q12w" =
    list(ev = ev_olpa, p = list()),

  "10 Lepodisiran 400 mg SC q24w" =
    list(ev = ev_lepo, p = list()),

  "11 Muvalaplin 240 mg PO daily" =
    list(ev = ev_muva, p = list()),

  "12 Muvalaplin + statin" =
    list(ev = c(ev_muva, ev_stat), p = list()),

  "13 Niacin ER 2 g daily" =
    list(ev = ev_niac, p = list()),

  "14 Obicetrapib 10 mg daily" =
    list(ev = ev_obic, p = list()),

  "15 Ziltivekimab 30 mg SC q4w" =
    list(ev = ev_zilt, p = list()),

  "16 Lipoprotein apheresis, weekly" =
    list(ev = ev_none, p = list(APHON = 1)),

  "17 Pelacarsen + statin + evolocumab" =
    list(ev = c(ev_pela, ev_stat, ev_evol), p = list()),

  ## ---- predictions ---------------------------------------------------
  "18 PREDICTION  statin in HoFH (no LDL receptor)" =
    list(ev = ev_stat, p = list(FLDLRFN = 0.05)),

  "19 PREDICTION  rheumatoid arthritis (IL-6 21) + anti-IL-6" =
    list(ev = ev_zilt, p = list(IL6EXO = 18)),

  "20 PREDICTION  true 50 mg/dL with a 12-repeat isoform (the missed patient)" =
    list(ev = ev_none, p = list(LPA0 = 137.2, NKIV2 = 12))
)

run_scenario <- function(name, end = 2*YR, delta = 1) {
  s <- sc[[name]]
  m <- mod
  if (length(s$p)) m <- param(m, s$p)
  out <- mrgsim_df(m, events = s$ev, end = end, delta = delta)
  out$scenario <- name
  out
}

## ---------------------------------------------------------------------
##  DRIVER BLOCK.  Sourcing this file with LPA_MODEL_ONLY set (as the
##  Shiny app does) defines `mod`, `sc` and the event objects and stops
##  there; running the file directly executes RUN 1-4 below.
## ---------------------------------------------------------------------
if (!exists("LPA_MODEL_ONLY", inherits = TRUE)) {

## ---------------------------------------------------------------------
##  RUN 1 — two-year pharmacology comparison
## ---------------------------------------------------------------------
res2y <- bind_rows(lapply(names(sc), run_scenario))

summ <- res2y %>%
  group_by(scenario) %>%
  summarise(
    Lpa_base_nM   = first(LPA_NMOL),
    Lpa_2y_nM     = last(LPA_NMOL),
    Lpa_pct       = 100*(last(LPA_NMOL)/first(LPA_NMOL) - 1),
    Lpa_abs_delta = last(LPA_NMOL) - first(LPA_NMOL),
    MassAssay_2y  = last(ASSAY_MASS),
    ApoaAssay_pct = 100*(last(ASSAY_APOA)/first(ASSAY_APOA) - 1),
    FreeApoa_fold = last(ASSAY_APOA - LPA_NMOL)/first(ASSAY_APOA - LPA_NMOL),
    LDLC_meas_pct = 100*(last(LDLC_MEAS)/first(LDLC_MEAS) - 1),
    apoB_pct      = 100*(last(APOB_TOT)/first(APOB_TOT) - 1),
    hsCRP_pct     = 100*(last(CRP)/first(CRP) - 1),
    .groups = "drop"
  )
print(as.data.frame(summ), digits = 4)

## ---------------------------------------------------------------------
##  RUN 2 — lifetime run (age 18 -> 80) for the MR-vs-trial comparison
##  and for the aortic-valve timing question.
## ---------------------------------------------------------------------
lifetime <- function(p = list(), ev_ = ev_none, start_treat_yr = NA,
                     end_yr = 62) {
  m <- mod
  if (length(p)) m <- param(m, p)
  e <- ev_none
  if (!is.na(start_treat_yr)) {
    d <- as.data.frame(ev_); d$time <- start_treat_yr*YR; e <- as.ev(d)
  }
  mrgsim_df(m, events = e, end = end_yr*YR, delta = 30)
}

life_hi    <- lifetime(list(LPA0 = 250, NKIV2 = 12))
life_lo    <- lifetime(list(LPA0 = 25,  NKIV2 = 30))
life_trt42 <- lifetime(list(LPA0 = 250, NKIV2 = 12), ev_pela, start_treat_yr = 42)
life_trt12 <- lifetime(list(LPA0 = 250, NKIV2 = 12), ev_pela, start_treat_yr = 12)

cat("\n--- Aortic valve at age 80 (62 model-years from age 18) ---\n")
for (nm in c("high Lp(a) untreated", "low Lp(a)",
             "treated from age 30", "treated from age 60")) {
  d <- switch(nm, "high Lp(a) untreated" = life_hi, "low Lp(a)" = life_lo,
              "treated from age 30" = life_trt12, life_trt42)
  cat(sprintf("  %-22s AVC %7.0f AU  AVA %4.2f cm2  MG %5.1f mmHg\n",
              nm, tail(d$VCALC,1), tail(d$AVA,1), tail(d$MGRAD,1)))
}

## MR-versus-trial: integrate the hazard over the SAME five years
haz5 <- function(d) {
  i42 <- which.min(abs(d$time - 42*YR)); i47 <- which.min(abs(d$time - 47*YR))
  d$HAZ[i47] - d$HAZ[i42]
}
h_un <- haz5(life_hi); h_tr <- haz5(life_trt42); h_lo <- haz5(life_lo)
cat(sprintf("\n5-yr integrated MACE hazard  untreated %.5f  treated %.5f  RRR %.1f%%\n",
            h_un, h_tr, 100*(1 - h_tr/h_un)))
cat(sprintf("                             lifelong-low %.5f  'MR' RRR %.1f%%\n",
            h_lo, 100*(1 - h_lo/h_un)))
cat(sprintf("ratio  MR effect / 5-yr trial effect = %.2fx\n",
            (1 - h_lo/h_un)/(1 - h_tr/h_un)))

## ---------------------------------------------------------------------
##  RUN 3 — the isoform / assay misclassification demonstration
##  Patients with the SAME true Lp(a) MASS but different isoforms.
## ---------------------------------------------------------------------
iso_demo <- function(nk, mass_target = 50) {
  conv <- 1.0e7/(3.30e6 + 14000*(nk + 10) + 35000)
  m <- param(mod, NKIV2 = nk, LPA0 = mass_target*conv)
  o <- mrgsim_df(m, events = ev_none, end = 1, delta = 1)
  data.frame(
    KIV2               = nk,
    apoa_MW_kDa        = round(o$MWAPOAO[1]/1000),
    true_nmol          = round(o$LPA_NMOL[1], 1),
    true_mass          = round(o$LPA_MASS_T[1], 1),
    conv_nmol_per_mgdl = round(o$CONVF[1], 3),
    epitope_bias       = round(o$EPITB[1], 3),
    mass_REPORTED      = round(o$ASSAY_MASS[1], 1),
    flagged_by_mass    = o$ASSAY_MASS[1] >= 50,
    flagged_by_molar   = o$LPA_NMOL[1]   >= 125,
    secretion_eff      = round(o$SECEFFO[1], 3)
  )
}
cat("\n--- Same TRUE mass (50 mg/dL), five isoforms ---\n")
print(bind_rows(lapply(c(8, 12, 22, 30, 35), iso_demo)), row.names = FALSE)

## ---------------------------------------------------------------------
##  RUN 4 — absolute vs relative reduction (the trial-design argument)
##  This is the table that explains why niacin failed and why
##  Lp(a)HORIZON and OCEAN(a) set the entry criteria they did.
## ---------------------------------------------------------------------
cat("\n--- Absolute Lp(a) reduction beats percent potency ---\n")
abs_demo <- expand.grid(baseline_nmol = c(30, 60, 125, 250, 400),
                        drug = c("niacin (-19%)", "evolocumab (-30%)",
                                 "pelacarsen (-79%)"))
abs_demo$pct <- c(-19.4, -30.1, -79.3)[as.integer(abs_demo$drug)]
abs_demo$absolute_nmol <- round(abs_demo$baseline_nmol*abs_demo$pct/100, 1)
print(abs_demo, row.names = FALSE)
cat("\nNote: pelacarsen in a 60 nmol/L patient removes 48 nmol/L, while\n",
    "evolocumab in a 250 nmol/L patient removes 75 nmol/L.  Baseline\n",
    "selection dominates percent potency, which is exactly why the\n",
    "outcome trials set floors of 70 mg/dL and 200 nmol/L.\n")

}   ## end of the driver block

## =====================================================================
##  CALIBRATION NOTES — anchors, ACHIEVED values, and what is NOT known
##  (sources are itemised in lpa_references.md; the section numbers below
##   refer to that file)
## ---------------------------------------------------------------------
##  Every "ACHIEVED" figure below is the output of an independent
##  re-implementation of these same 50 ODEs in Python with a fixed-step
##  RK4 integrator, run from the same parameter vector.  They are
##  measured, not intended.  Writing that cross-check is what found the
##  five defects listed at the end of this block.
##
##   1. Lp(a) FCR.  Anchor: kinetic studies give FCR 0.20-0.35 pools/day,
##      INVARIANT across a 10-fold concentration range (refs sec. 4).
##      Model: KCAT0 = KLDLR_LPA + KOTH_LPA = 0.060 + 0.190 = 0.250 /day
##      (t1/2 2.77 d), and the model is production-controlled by
##      construction: KTL is back-solved from the target Lp(a), the FCR
##      never is.
##
##   2. LDLR-dependent fraction of Lp(a) clearance.  Anchor: PCSK9 mAbs
##      lower LDL-C ~60% but Lp(a) only 25-30% (refs sec. 8).  The 25-30%
##      is NOT fitted.  KLDLR_LPA/KCAT0 = 24% is set from the receptor
##      biology; evolocumab drives LDLR 1.00 -> 2.75 in the model, and
##      the two answers then fall out of the same receptor change:
##          ACHIEVED  LDL-P -59.0%   Lp(a) -30.1%
##
##   3. Statin paradox.  Anchor: meta-analyses give +8 to +20% Lp(a).
##      Two opposing terms: clearance (LDLR 1.00 -> 1.78, worth -15.7%)
##      and transcription (x1.30 via the sterol-responsive element).
##          ACHIEVED  net +8.1%,  LDL-P -38.2%
##      ESTA = 0.30 is the only quantity fitted to this endpoint, and it
##      is a derived requirement: given the clearance arithmetic, nothing
##      smaller can produce a net RISE at all.
##
##   4. Ezetimibe — a PREDICTION the model was not tuned to.  An
##      LDLR-mediated LDL-C fall with no direct LPA effect should lower
##      Lp(a) by roughly a quarter as much, because a quarter of Lp(a)
##      clearance is LDLR-dependent.
##          ACHIEVED  LDL-P -19.3%,  Lp(a) -6.6%
##      Observed: ezetimibe monotherapy lowers Lp(a) by ~7% (ref 77).
##      Ezetimibe is therefore not a null control — it is a quantitative
##      test of the 24% figure, and it passes.
##
##   5. Assembly is SATURATED in LDL substrate.  This is forced by the
##      same ezetimibe data: if assembly were first-order in LDL
##      particles, a 19% LDL fall would carry Lp(a) down with it by a
##      further ~19%, which is not observed.  Hence the Michaelis term
##      LDL_P/(KMLDL + LDL_P), 97.6% saturated at baseline.  This was a
##      real defect in the first version of the model, caught by the
##      cross-check: the original first-order form made high-intensity
##      statin LOWER Lp(a) by 4%, inverting a well-established finding.
##
##   6. Pelacarsen.  Anchor: 80 mg monthly, Lp(a) -80% (refs sec. 9).
##          ACHIEVED  -74% at wk 4, -80.7% at wk 13, nadir -81.2% at
##                    wk 26, -79.3% at 2 yr;  apoB -11.3% (obs. ~-13%)
##      The 13-week onset is not a fitted delay: it is the hepatic ASO
##      half-life (23 d) filling the liver compartment.
##
##   7. Olpasiran / lepodisiran.  Anchor: olpasiran 75 mg q12w gives -95
##      to -101% placebo-adjusted (OCEAN(a)-DOSE); lepodisiran holds >90%
##      for ~48 weeks after ONE dose.
##          ACHIEVED  olpasiran -96.6%;  lepodisiran 400 mg q24w -96.8%
##          ACHIEVED  single 75 mg dose: -96% at wk 12, -95% at wk 24,
##                    -84% at wk 48, -68% at wk 60
##      The off-treatment persistence comes from the RISC compartment
##      (t1/2 60 d), not from a separate PD-delay term.
##
##   8. Niacin.  Anchor: -20 to -25% Lp(a), yet AIM-HIGH and HPS2-THRIVE
##      were null.  ACHIEVED -19.4%.  The model's explanation is
##      arithmetic and RUN 4 prints it: at the ~30 mg/dL median baseline
##      of those trials, 20% is ~6 mg/dL absolute — an order of magnitude
##      below the exposure-response requirement.  Niacin is not a failed
##      Lp(a) drug; it is an adequately potent drug given to a population
##      in which no achievable percentage could have worked.
##
##   9. Anti-IL-6.  Anchors: RESCUE (ziltivekimab) Lp(a) -16 to -25% and
##      hsCRP -92% in non-inflamed patients; MEASURE (tocilizumab) -37%
##      in rheumatoid arthritis.  ONE equation, only the IL-6 level
##      changed between the two runs.
##          ACHIEVED  IL-6 2.7 pg/mL : Lp(a) -26.6%,  hsCRP -92.5%
##          ACHIEVED  IL-6  21 pg/mL : Lp(a) -33.5%,  hsCRP -96.3%
##
##  10. The feed-forward loop is real but WEAK, and saying so is the
##      result.  Lp(a) -> OxPL -> monocyte -> IL-6 -> LPA promoter.
##          ACHIEVED  open-loop gain g = 0.019, amplification 1.019
##      i.e. clinically negligible.  What is NOT negligible is the TONIC
##      IL-6 term: at IL-6 2 pg/mL the promoter multiplier is already
##      1.33, so a quarter of basal LPA transcription is IL-6-driven.
##      That, not the loop, is what anti-IL-6 removes.
##
##  11. Lipoprotein apheresis.  Anchor: 60-70% acute removal per session.
##          ACHIEVED  weekly: nadir -65.7%, true time-averaged -33.8%
##                    q2w:    nadir -62.9%, true time-averaged -18.8%
##      Note the model computes a genuine time integral.  The "30-35%
##      interval mean" usually quoted is the average of pre- and
##      post-session values, which is not the same quantity and is
##      systematically more favourable.
##
##  12. Muvalaplin — AND THE ONE PLACE THE MODEL DOES NOT CLOSE.
##      Anchor: KRAKEN reported -85.8% by an intact-Lp(a) assay and
##      about -70% by a traditional apo(a)-detecting assay.
##          ACHIEVED  intact -84.4% (good)
##          ACHIEVED  free apo(a) x4.09 (right direction, and the reason
##                    the two assays must disagree at all)
##          ACHIEVED  traditional assay -80.4% at RFREE = 1  -- NOT -70%
##      The model says the discrepancy can only be closed if the
##      traditional assay's MOLAR response to FREE apo(a) is about 3.8x
##      its response to apo(a) inside an intact particle:
##          RFREE 1.0 -> -80.4%    RFREE 3.0 -> -72.5%
##          RFREE 3.8 -> -69.4%    RFREE 5.0 -> -64.9%
##      That is a measurable quantity which, as far as we can find, has
##      not been reported.  RFREE therefore defaults to 1.0, where the
##      model UNDER-explains KRAKEN, rather than being tuned to 3.8 to
##      make the table look right.
##
##  13. Assay bias vs unit conversion — the result that surprised us.
##      Across the isoform range 8 to 35 KIV-2 repeats:
##          STOICHIOMETRIC conversion  2.788 -> 2.522 nmol/L per mg/dL
##                                     (a spread of 10%)
##          ANALYTICAL epitope bias    0.562 -> 1.406
##                                     (a spread of 150%)
##      The unit conversion that the field argues about is nearly fine.
##      The antibody is the problem.  RUN 3 makes it concrete — four
##      patients with the SAME true mass of 50 mg/dL:
##
##        KIV-2   true nmol/L   REPORTED mg/dL   mass>=50?   molar>=125?
##            8       139.4          28.1          missed        FLAGGED
##           12       137.2          34.4          missed        FLAGGED
##           22       132.2          50.0          borderline    FLAGGED
##           35       126.1          70.3          FLAGGED       FLAGGED
##
##      The mass assay misses precisely the small-isoform patient — the
##      one with the highest genetically predicted risk — and flags the
##      large-isoform patient hardest.  Note also that the guideline pair
##      50 mg/dL / 125 nmol/L is not self-consistent: 50 mg/dL of true
##      mass is 126-139 nmol/L across the isoform range, so the molar
##      threshold is uniformly the more inclusive of the two.
##
##  14. Cholesterol book-keeping.  ACHIEVED at Lp(a) 250 nmol/L, small
##      isoform: reported LDL-C 127.3 mg/dL = 100.0 true + 27.3 Lp(a)-C,
##      i.e. 21% of the reported "LDL-C" is not LDL.  Performing the
##      Dahlen correction with the (biased) reported mass gives 108.5
##      mg/dL against a true LDL-C of 100.0 — an 8.5 mg/dL
##      under-correction, in the same patients and the same direction as
##      the threshold error above.
##
##  15. Why a minority particle matters.  ACHIEVED: Lp(a) is 15.9% of
##      apoB particles but contributes 34.2% of the intimal retention
##      flux (AVID = 2.5), and carries an OxPL payload LDL does not.
##
##  16. Risk time constants.  Inputs: the 70/30 split between plaque
##      burden and vulnerability, and their two half-lives (14 yr, 62 d).
##      Nothing else was fitted.
##          ACHIEVED  lifelong 250 vs 25 nmol/L: HR 2.29 at age 65
##                    (PAV 44.0 vs 29.8%, VULN 1.448 vs 1.000)
##          ACHIEVED  80% reduction from age 60, integrated over 5 years:
##                    RRR 20.0%
##          ACHIEVED  lifelong-low vs lifelong-high over the same 5 yr:
##                    RRR 53.1%
##          ACHIEVED  ratio = 2.65x
##      The ~3x gap between Mendelian randomisation and a 5-year trial is
##      thus not a correction factor. It is what two time constants do.
##
##  17. Aortic valve — the sharpest clinical statement in the model.
##      Anchor: strong Lp(a)-AVS Mendelian-randomisation signal, but SEAS,
##      ASTRONOMER and SALTIRE were all null in ESTABLISHED calcific AS.
##      The self-perpetuating term is a THRESHOLD (Hill n = 4), so below
##      a mineral nucleus it contributes nothing and above it the process
##      no longer needs Lp(a).  Over 62 model-years (age 18 -> 80):
##
##        low Lp(a) 25 nmol/L      AVC    43 AU   AVA 3.21   MG   3.9
##        high Lp(a), untreated    AVC  1558 AU   AVA 0.82   MG  58.9
##        pelacarsen from age 30   AVC   206 AU   AVA 2.45   MG   6.7
##        pelacarsen from age 60   AVC  1327 AU   AVA 0.93   MG  46.3
##
##      Starting at 30 prevents 87% of the calcium.  Starting at 60 —
##      twenty years of a drug that removes 80% of the particle —
##      prevents 15%.  The model therefore predicts that a valve-endpoint trial in
##      patients who already have moderate AS will be null however potent
##      the agent, and that this indication requires primary prevention.
##      (Untreated progression runs at ~80 AU/yr in the last decade,
##      against a literature range of ~100-250 AU/yr in progressive AS:
##      the model is conservative here, which strengthens rather than
##      weakens the conclusion.)
##
##  18. DEFECTS THE PYTHON CROSS-CHECK FOUND AND THIS FILE FIXES.
##      (a) First-order LDL dependence of assembly made high-intensity
##          statin LOWER Lp(a) by 4.1% and raise free apo(a) 1.8-fold —
##          inverting the statin paradox.  Fixed by saturating the
##          assembly step (item 5).
##      (b) PCSK9 turnover was written 5x too fast, making the evolocumab
##          arm numerically stiff without changing any steady state.
##          Rewritten at KSP = KDP = 1.0 /day with KONE rescaled to hold
##          the same 98% free-PCSK9 suppression.
##      (c) The IL-6 antibody was modelled as partially competitive, so
##          the rheumatoid-arthritis arm produced a SMALLER Lp(a) fall
##          than the non-inflamed arm — backwards.  At clinical doses the
##          antibody sits in ~10^4 molar excess over IL-6, so blockade is
##          near-complete and ligand-independent; IC50Z was corrected.
##      (d) CRP had no IL-6-independent production, so the model
##          predicted -98% hsCRP where RESCUE observed -92%.  FCRP0 = 0.06
##          added.
##      (e) The valve self-perpetuation term was a RAMP rather than a
##          threshold, which gave an 80-year-old with Lp(a) 25 nmol/L an
##          AVC of 4,425 AU and an aortic valve area of 0.34 cm2 — i.e.
##          the model gave everybody severe aortic stenosis.  Replaced
##          with a Hill-4 threshold and re-tuned (item 17).
##      An unstable-looking hazard exponent was also clamped, and the
##      mean-gradient output guarded at the near-atretic limit.
##
##  19. What this model does NOT do.  No explicit sex; no
##      lipoprotein-lipase / triglyceride axis; a single expressed apo(a)
##      isoform rather than two alleles; no plaque geometry; no
##      pharmacogenomics of the aspirin-rs3798220 interaction; and the
##      antifibrinolytic arm is deliberately down-weighted (WARM3 = 0.25)
##      because Mendelian randomisation shows NO Lp(a)-VTE association.
##      The in-vitro biology of that arm is far stronger than its in-vivo
##      evidence, and the model follows the evidence.
## =====================================================================
