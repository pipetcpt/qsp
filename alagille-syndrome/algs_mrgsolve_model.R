# =====================================================================
# Alagille Syndrome (ALGS) — mrgsolve QSP model
#   Author : Claude Code Routine (2026-07-30)
#
#   Scope  : JAG1 (94.3%) / NOTCH2 (2.5%) haploinsufficiency -> ductal
#            plate malformation -> INTERLOBULAR BILE DUCT PAUCITY
#            (duct : portal-tract ratio DPR < 0.9) -> a hard ceiling on
#            biliary output -> chronic cholestasis -> pruritus,
#            xanthoma, fat-soluble-vitamin failure, growth failure,
#            biliary fibrosis, portal hypertension, transplant.
#
#   PK/PD  : IBAT (ASBT / SLC10A2) inhibitors — maralixibat 380 ug/kg/d,
#            odevixibat 120 ug/kg/d, linerixibat class; plus UDCA,
#            cholestyramine, rifampicin, naltrexone, and the surgical
#            comparator PEBD (partial external biliary diversion).
#
# ---------------------------------------------------------------------
#   WHY THIS MODEL IS BUILT THE WAY IT IS
#
#   ALGS cholestasis is not one quantity.  It is the collision of two
#   fluxes that are set INDEPENDENTLY and that every clinical
#   measurement mixes together:
#
#     J_DUCT — how much bile the liver can physically push into the gut.
#              Set developmentally by Notch dose; slowly repaired
#              postnatally by ductular reaction; moved by NO drug here.
#     R      — how much of that bile comes straight back.  The ONLY
#              thing an IBAT inhibitor touches.
#
#   Everything below follows from writing those two separately.  Every
#   number quoted is produced by running this file (run_all_algs()).
#
#   AXIS 1 — THE ENTEROHEPATIC FRACTION IS THE DRUGS CEILING, AND IT IS
#            A DUCT PROPERTY, NOT A BILE-ACID PROPERTY.
#     Define PHI_EHC = R / (S + R): the share of the hepatocytes
#     bile-acid input that arrived back from the gut rather than being
#     freshly made.  An ideal IBAT inhibitor removes PHI_EHC of that
#     input and not one percent more.  In an unobstructed liver
#     PHI_EHC = 0.98, because bile acid is recycled ~50 times per
#     molecule synthesised.  Duct paucity caps biliary output, which caps
#     R with it, and the model — calibrated ONLY to the ASSERT sBA ratio
#     — puts the average trial patient at
#        PHI_EHC = 0.42   (JCAP = J_DUCT/S = 1.99).
#     Duct paucity has already destroyed more than half of the drugs
#     target before the first dose is given.  Sweeping duct capacity
#     moves the achievable 24-week sBA response monotonically from -92%
#     (JCAP 6.1) to -9% (JCAP 0.38): the drug does not fail because the
#     bile acids are too high, it fails because there are no ducts.
#
#   AXIS 2 — A HYPOTHESIS THE MODEL REJECTED, AND WHAT REPLACED IT.
#     The design intent was that below some duct capacity an IBAT
#     inhibitor would INVERT — blocking ASBT removes ileal FGF19,
#     derepresses CYP7A1, and the newly synthesised bile acid would have
#     no duct to leave by, so sBA would RISE on treatment.  The model
#     says that does not happen, for a reason worth more than the
#     hypothesis was: in ALGS the FGF19 signal is ALREADY floored before
#     treatment (baseline FGF19 = 0.04 of normal, synthesis already
#     derepressed 2.14-fold), so the synthesis reserve that would drive
#     the inversion has been spent.  Blockade is therefore futile in
#     severe paucity, not harmful.  The prediction that survives is a
#     graded threshold, and it lands somewhere unexpected:
#        response < 30%  below JCAP 1.25, i.e. bilirubin 5.42 mg/dL
#        response < 15%  below JCAP 0.60, i.e. bilirubin 10.70 mg/dL
#     GALA — a natural-history cohort, keyed on native liver survival,
#     with no drug in it — split its prognostic strata at total
#     bilirubin 5.0 and 10.0 mg/dL.  This model never saw those numbers.
#     The RATIO of the two boundaries is 1.97 here and 2.00 in GALA, and
#     unlike the boundaries themselves the ratio is completely
#     independent of the one assumed calibration anchor (FAILURE 2).
#     Two datasets with nothing in common appear to be measuring the same
#     duct-capacity threshold from opposite directions.
#
#   AXIS 3 — ITCH IS NOT A BILE-ACID MEASUREMENT, AND THE ASSERT PLACEBO
#            ARM PROVES IT ARITHMETICALLY.
#     In ASSERT the placebo arms itch FELL by 0.8 of the 1.7 points the
#     drug arm fell (47% of the on-drug improvement is not drug), while
#     the placebo arms bile acids ROSE by 22 umol/L.  One axis has a
#     large favourable placebo response, the other a negative one; they
#     cannot be the same variable.  The consequence is quantitative:
#        controlled slope (ASSERT)   = 0.9/113 = 0.0080 pts per umol/L
#        uncontrolled slope (ASSERT) = 1.7/90  = 0.0189
#        uncontrolled slope (ICONIC) = 1.6/96  = 0.0167
#     The two SINGLE-ARM slopes agree with each other to 13% and are both
#     ~2.1x the placebo-controlled truth, so any model calibrated on
#     single-arm cholestasis data over-attributes itch to bile acids
#     about two-fold.  This model carries the components separately and
#     reproduces the controlled slope (0.0083 vs 0.0080 published) and
#     the drug-attributable itch difference (-0.90 vs -0.90) exactly.
#     A note on the psychophysics: fitting ASSERT REJECTED a logarithmic
#     (Weber-Fechner) itch law.  Forced through a saturating 0-4 map with
#     a baseline at 2.8, the trial demands a Stevens exponent of 1.66 —
#     supralinear, which no itch psychophysics supports.  Treating the
#     0-4 scale as a ceiling on the INSTRUMENT rather than on the
#     sensation puts the operating point in the near-linear range and the
#     solved exponent falls to 0.95.  The reported scale is not the
#     latent variable, and pretending otherwise distorts the pharmacology.
#
#   AXIS 4 — THE SURVIVAL BENEFIT IS LARGER THAN ITS OWN BILE-ACID
#            EFFECT CAN EXPLAIN.  The model reports the gap and does not
#            close it.
#     GALA supplies its own exposure-hazard gradient: relative to
#     TB < 5.0 mg/dL the transplant hazard is 4.8x at 5-10 and 15.6x at
#     >=10, which by stratum midpoints is a power law of exponent
#     n = 1.5-1.85 (n = 1.85 reproduces the severe stratum best, 17.3 vs
#     15.6 observed).  Running six years of maralixibat against this
#     models own natural history at n = 1.85 gives HR = 0.414.  The
#     published maralixibat-vs-GALA comparison reports HR = 0.305.  To
#     reach 0.305 from the drugs own bile-acid effect the exponent would
#     have to be ~2.9 — and n = 2.9 back-predicts a hazard ratio of 87
#     for the GALA TB >= 10 stratum where 15.6 was observed, a 5.6-fold
#     over-shoot.  No single exponent fits both.  A residual hazard ratio
#     of 0.305/0.414 = 0.74 is therefore NOT explained by bile acids.
#     Worse, gamma frailty (see below) ATTENUATES population hazard
#     ratios toward 1 over time, so the individual-level effect implied
#     by an observed population HR of 0.305 is larger still and the gap
#     widens.  The candidates — six-year drug-persistence selection
#     against a comparator matched on baseline labs only, a
#     bile-acid-independent benefit, or curvature the cross-sectional
#     strata cannot see — are not separable by any published data.  The
#     model ships with n = 1.6 so that it UNDER-predicts the trial;
#     scenario 12 sweeps n and prints both errors side by side.
#
#   AXIS 5 — WHERE IN THE GUT A DRUG ACTS DECIDES WHETHER IT STARVES THE
#            CHILD.
#     ASBT sits in the TERMINAL ILEUM, downstream of the duodenal
#     micellar window where fat and vitamins A/D/E/K are absorbed.  At
#     steady state a duct-limited liver delivers J_DUCT to the duodenum
#     whatever ASBT is doing.  Over two years in severe paucity the model
#     gives, on the same patient:
#        odevixibat     sBA -46%,  fat absorption 0.662 -> 0.652,
#                       vitamin D 19.1 -> 18.8, height z IMPROVES
#        cholestyramine sBA -37%,  fat absorption 0.662 -> 0.294,
#                       vitamin D 19.1 -> 8.5, INR 1.47 -> 2.12,
#                       height z -1.51 -> -2.38
#     Both drugs remove bile acids; only the proximal-acting one causes
#     steatorrhoea, and the difference is derived from anatomy rather
#     than asserted.  The model also finds its own boundary case: in MILD
#     paucity the duct is NOT saturated, duodenal delivery does fall with
#     the drug (6.2 -> 4.3 mM), fat absorption drops 0.844 -> 0.752, and
#     the contrast narrows.  The safety argument for IBAT inhibitors
#     depends on the liver being duct-limited.
#
#   AXIS 6 — THE POPULATION CEILING.  Cardiac disease dominates ALGS
#     mortality in year one and vasculopathy/intracranial haemorrhage
#     later; both are Notch-dose diseases no cholestasis drug touches.
#     In the model the liver carries 77.9% of the fatal hazard to age 18,
#     cardiac 10.4% and vascular 11.7%, so a therapy achieving a
#     liver-specific HR of 0.386 moves all-cause hazard only to 0.522.
#
#   AXIS 7 — THE GALA CURVE IS HETEROGENEITY, NOT BIOLOGY.
#     No single trajectory fits GALA native liver survival: its hazard is
#     strongly front-loaded and every mechanistic trajectory here is not
#     (best single-trajectory fit SSE 3.3e-2, visibly wrong at age 5).
#     Admitting a gamma frailty of variance 2.94 fits all three time
#     points 65-fold better (SSE 5.0e-4; 67.7/52.6/41.3 against GALA
#     66.8/54.4/40.3).  A variance of 2.94 means the standard deviation
#     of individual hazard is 1.7 times its mean — the notorious ALGS
#     phenotype spread (same variant, infant transplant or asymptomatic
#     parent), measured rather than described.  Consequence: population
#     curves must use the _POP outputs.  The plain exp(-H) outputs are
#     INDIVIDUAL-level and comparing them to a cohort study is a category
#     error that this file made once (scenario 1 reported 1.5% 18-year
#     native liver survival against the true 40.3%) before it was caught.
#
# ---------------------------------------------------------------------
#   WHAT THIS MODEL DOES NOT DO — reported, not repaired.
#
#   FAILURE 1 — PHI_EHC IS NOT SEPARATELY IDENTIFIABLE FROM ONE TRIAL.
#     Duct capacity JCAP, achieved fractional ASBT blockade I, and the
#     synthesis reserve enter the steady-state sBA ratio as one
#     combination; a single number cannot resolve three unknowns.
#     Scenario 6 re-solves JCAP for each assumed blockade so that ASSERT
#     is still matched exactly, and JCAP moves 3.58 -> 1.95 (1.8-fold)
#     while PHI_EHC moves 0.57 -> 0.42 (1.36-fold) as the assumed maximum
#     blockade goes 0.50 -> 0.92.  PHI_EHC is the more robust quantity but
#     neither is pinned.  Every JCAP statement is conditional on
#     IMAXA = 0.92 (achieved blockade 0.70 at the ASSERT dose).
#
#   FAILURE 2 — THE BILIRUBIN BOUNDARIES REST ON ONE ASSUMED ANCHOR.
#     Mapping duct capacity to bilirubin needs the ASSERT populations
#     mean baseline total bilirubin, which the primary publications
#     abstract does not give; it is anchored at 3.5 mg/dL.  Scenario 6b
#     sweeps it: anchors of 3.0 / 3.5 / 4.0 put the 30% boundary at
#     4.58 / 5.42 / 6.26 and the 15% boundary at 8.93 / 10.70 / 12.48.
#     The GALA cut-points are bracketed, but the point estimates are not
#     evidence.  Only the RATIO of the two boundaries (1.97) is
#     anchor-independent, which is why Axis 2 leads with the ratio.
#
#   FAILURE 3 — THE BA-INDEPENDENT ITCH FRACTION IS ASSUMED, NOT SOLVED.
#     The ASSERT placebo arm quantifies the NON-DRUG component (47%).  It
#     says nothing about how much itch would remain if bile acids were
#     driven to zero.  That is constrained only by the surgical
#     literature, so the autotaxin/LPA arm is set to leave ~15% residual
#     itch at complete normalisation.  The one mitigating finding: the
#     solved residual is 12-17% across every combination of Stevens
#     exponent (0.75-1.05) and assumed BA-independent share (8-18%)
#     tried, so the conclusion is less anchor-sensitive than the
#     parameter is.  It remains an assumption.
#
#   FAILURE 4 — SPONTANEOUS IMPROVEMENT IS ONLY HALF REPRODUCED.
#     ALGS cholestasis is described as improving after age ~5-10.  The
#     ductular-repair term does produce that — with the fibrosis term
#     switched off, bilirubin falls 3.50 -> 2.20 from age 1 to 18 as duct
#     capacity recovers 1.99 -> 3.29.  In the full model it does NOT:
#     fibrosis-driven bilirubin rise outruns duct recovery and bilirubin
#     goes 3.50 -> 4.48 instead.  Which of those a real child follows is
#     the whole clinical question and this model does not decide it,
#     because no serial-biopsy DPR series exists to fit the repair rate
#     to.  Serum bile acid barely moves either way (236.9 -> 233.4),
#     which is itself a prediction: partial duct recovery should improve
#     BILIRUBIN long before it improves BILE ACIDS, because faecal loss
#     stays negligible until duct capacity approaches normal.
#
#   FAILURE 5 — HAZARDS ARE CALIBRATED TO A COHORT, APPLIED TO AN
#     INDIVIDUAL.  See Axis 7: the frailty variance needed to fit GALA is
#     so large that the cohort curve is mostly a statement about
#     between-patient variance.  Scenario 1 tracking GALA is a
#     calibration check, not a prognosis for any patient.
#
#   FAILURE 6 — NO CIRRHOTIC PHARMACOLOGY AND NO PLACEBO-ARM PROGRESSION.
#     Past FIB ~3 the model keeps applying the same duct/synthesis
#     structure, but a decompensated ALGS liver has hepatocellular
#     failure physiology this model does not contain.  Separately, the
#     ASSERT placebo arms bile acids ROSE 22 umol/L over 24 weeks; the
#     model reproduces the drug/placebo RATIO but holds its untreated arm
#     flat (-0.9 umol/L), so it does not contain whatever drives that
#     rise, and the absolute on-drug sBA change is correspondingly
#     over-stated (-109 vs -90).
#
#   NINE DEFECTS FOUND BY RUNNING IT, LISTED BECAUSE THEY WOULD HAVE
#   CHANGED CONCLUSIONS:  (1) $MAIN re-asserted initial conditions on
#   every run, silently overriding init(), so every duct-capacity
#   scenario reported an IDENTICAL baseline bilirubin of 3.5 and the
#   entire Axis 2 threshold was unmeasurable; (2) scenario 1 compared
#   individual-level exp(-H) survival to a cohort curve and reported 1.5%
#   18-year native liver survival against a true 40.3%; (3) PHI_EHC
#   initially included the hepatocyte->plasma->hepatocyte futile cycle,
#   which pins it at ~0.98 in every patient and destroys the quantity;
#   (4) EPS and THETA are reserved words in mrgsolve and silently failed
#   to compile; (5) apostrophes in comments terminated the single-quoted
#   R model string; (6) $TABLE captured locals left over from the last
#   $ODE derivative call rather than recomputing them at the output time;
#   (7) the first fibrogenesis rate drove FIB to 2 within 400 days and
#   the corrected one froze it at 1.0 for 17 years; (8) the bolus-vs-
#   continuous dosing check sampled on integer days, which is exactly the
#   dosing period, and reported a 39% discrepancy that was pure aliasing
#   (5.7% on a fine grid); (9) the cholestyramine arm bound luminal bile
#   acid for micelle formation but not for ASBT uptake, making a
#   sequestrant look like a pure malabsorption agent with no bile-acid
#   effect at all.
# =====================================================================

library(mrgsolve)
library(dplyr)

algs_code <- '
$PROB
# Alagille Syndrome (ALGS) QSP model
# 44 ODE compartments: 4 drug PK + 40 disease / PD / clinical
# Time unit: DAYS.  Age enters as AGE0 + TIME/365.25.

$PLUGIN autodec

$PARAM @annotated
// =====================================================================
// PATIENT / STRUCTURAL
// =====================================================================
WT      :  15    : Body weight (kg)
AGE0    :   3    : Age at simulation start (yr)
DPR0    :   0.2280: Bile duct : portal tract ratio at start (normal 0.9-1.8)
DPRN    :   1.2  : Reference NORMAL duct : portal tract ratio (-)
JCAPN   :  55    : Duct capacity of a normal liver, as multiples of synthesis (-)
GAMD    :   2.0  : Exponent, duct capacity vs DPR (conductance ~ n_ducts x area)
USEINIT :   0    : 1 = use externally supplied init(), skip the $MAIN defaults (-)

// =====================================================================
// BILE ACID SYNTHESIS AND FEEDBACK
// =====================================================================
SYNKG   :   8.67 : Baseline BA synthesis (umol/kg/day)
SIGMAX  :   2.5  : Maximum synthesis reserve, CYP7A1 derepression (fold)
KSIG    :   0.30 : Half-effect fractional FGF19 loss for synthesis reserve (-)
KFXR    :   0.5  : Strength of the hepatic FXR-SHP arm (-)
HFXR    :   1.0  : Hill exponent, hepatic FXR arm (-)
KCYP    :   0.15 : CYP7A1 turnover rate (1/day)
KFGF    :   4.0  : FGF19 turnover rate (1/day)

// =====================================================================
// BILE ACID KINETICS
// =====================================================================
KMH     :  39.5  : Michaelis constant, canalicular + ductal export (umol)
KSP0    :   5.40 : Basolateral spillover rate constant MRP3/4-OSTab (1/day)
KIND    :   0.9  : Adaptive induction of spillover by hepatocyte load (-)
KBD     :   8.0  : Gallbladder emptying rate (1/day)
KDI     :  12.0  : Duodenum -> ileum transit (1/day)
KIL     :   8.0  : Terminal ileum turnover (1/day)
EPS0    :   0.98 : ASBT fractional reabsorption, uninhibited (-)
KPH     :  24.0  : Portal transit rate (1/day)
FPE     :   0.85 : Hepatic first-pass extraction of portal BA (-)
KFEC    :   1.2  : Colonic transit to faeces (1/day)
KSYSCL  :   6.0  : Systemic BA clearance (1/day)
FREUP   :   0.957: Fraction of systemic clearance that is hepatic re-uptake (-)
VDBA    :   0.30 : BA distribution volume (L/kg)
KCDU    :   0.463: Duodenal BA concentration coefficient (mM per (umol/day)^PCDU)
PCDU    :   0.35 : Compressive exponent, gallbladder concentration (-)
FREF    :  49.0  : Normal enterohepatic multiplier, ileal uptake / synthesis (-)

// =====================================================================
// IBAT INHIBITOR PK/PD
// =====================================================================
DOSEKG  :   0    : IBAT inhibitor dose (ug/kg/day)
POT     :   1.0  : Potency vs odevixibat reference (odev 1.0, maralixibat 0.28)
TSTART  :   0    : Treatment start (day)
TSTOP   :   0    : Treatment stop (day)
KTR1    :   6.0  : Drug depot -> mid gut (1/day)
KTR2    :   6.0  : Drug mid gut -> terminal ileum (1/day)
KTR3    :   6.0  : Drug ileum -> faeces (1/day)
FSYS    :   0.005: Fraction of drug absorbed systemically (-)
KSYSD   :   3.0  : Systemic drug elimination (1/day)
IMAXA   :   0.92 : Maximum achievable fractional ASBT blockade (-)
IC50KG  :   6.3  : ASBT IC50 as ileal drug amount (ug/kg)

// =====================================================================
// NON-IBAT THERAPY (flags / intensities)
// =====================================================================
UDCAF   :   0    : UDCA on/off (0/1) at 20 mg/kg/day
UDCACH  :   0.18 : UDCA choleresis, fractional gain in duct flux (-)
UDCAHP  :   0.22 : UDCA pool hydrophilicity, fractional cut in cytotoxicity (-)
CHOLB   :   0    : Cholestyramine luminal BA binding fraction (0-1)
RIFF    :   0    : Rifampicin on/off (0/1)
RIFITCH :   0.35 : Rifampicin fractional cut in itch drive (-)
RIFHEP  :   1.0  : Rifampicin hepatotoxicity multiplier on ALT (-)
NALF    :   0    : Naltrexone on/off (0/1)
NALITCH :   0.20 : Naltrexone fractional cut in central itch relay (-)
PEBDF   :   0    : PEBD (partial external biliary diversion) on/off (0/1)
PEBDEPS :   0.02 : Residual ASBT-accessible fraction after PEBD (-)
MCTF    :   0    : MCT-enriched formula on/off (0/1)
MCTGAIN :   0.35 : Energy recovered by MCT independent of micelles (-)

// =====================================================================
// LIVER STRUCTURE, INJURY, FIBROSIS
// =====================================================================
KDR     :   0.000022: Ductular repair rate (1/day)
DPRMAX  :   0.95 : Ceiling on postnatal DPR recovery (-)
KDLOSS  :   0.0000075: Fibrosis-driven duct loss (1/day per FIB unit)
KDU     :   0.05 : Ductular reaction turnover (1/day)
KDRX    :   0.35 : Profibrotic weight of the ductular reaction (-)
HFIB    :   1.6  : Hill exponent, hepatocyte BA -> stellate activation (-)
KHA     :   0.010: Stellate activation rate (1/day)
KHD     :   0.010: Stellate deactivation rate (1/day)
HSCMAX  :   1.0  : Maximum activated stellate fraction (-)
KFF     :   0.00042: Fibrogenesis rate (1/day)
KFR     :   0.00006: Fibrosis regression rate (1/day)
FIBMAX  :   4.0  : Maximum fibrosis stage (-)
KAPO    :   0.9  : Apoptotic flux per unit relative hepatocyte BA (-)

// =====================================================================
// BIOCHEMISTRY
// =====================================================================
KBI     :   0.20 : Bilirubin equilibration rate (1/day)
BILMIN  :   0.40 : Bilirubin floor (mg/dL)
KBILC   :   5.92 : Duct-limited bilirubin coefficient (mg/dL per JCAP unit)
KBILF   :   0.28 : Bilirubin rise per FIB^2 (mg/dL)
KGG     :   0.15 : GGT equilibration rate (1/day)
GGTMIN  :  20    : GGT floor (U/L)
KGGT    : 320    : GGT per unit relative hepatocyte BA (U/L)
KAL     :   0.20 : ALT equilibration rate (1/day)
ALTMIN  :  25    : ALT floor (U/L)
KALT    : 140    : ALT per unit apoptotic flux (U/L)

// =====================================================================
// LIPIDS AND XANTHOMA
// =====================================================================
KLPX    :   0.06 : Lipoprotein X turnover (1/day)
TC0     : 150    : Cholesterol floor (mg/dL)
KTC     : 480    : Cholesterol per unit lipoprotein X (mg/dL)
XTHR    : 500    : Cholesterol threshold for xanthoma formation (mg/dL)
KXF     :   0.0022: Xanthoma formation rate (1/day per mg/dL over threshold)
KXR     :   0.010: Xanthoma regression rate (1/day)

// =====================================================================
// PRURITUS
// =====================================================================
SBAREF  :  20    : Reference bile acid concentration for the itch law (umol/L)
QBA     :   0.95 : Stevens exponent, bile acid -> itch  [SOLVED from ASSERT]
WBA     :   0.0644: Weight, bile-acid arm of itch drive (-)   [SOLVED from ASSERT]
WIND    :   0.1459: Weight, autotaxin/LPA (BA-independent) arm (-)  [ASSUMED: FAILURE 3]
WSENS   :   0.30 : Weight, central sensitisation (-)
K50I    :   3.853: Half-maximal itch drive (-)   [SOLVED from ASSERT]
// The reported 0-4 itch scale is a ceiling on the INSTRUMENT, not on the
// sensation.  Putting a baseline of 2.8 at 70% of a saturating 0-4 map
// forces a SUPRALINEAR Stevens exponent (1.6) that no psychophysics
// supports; a wider latent scale puts the operating point in the
// near-linear range and the solved exponent lands sub-linear.
ITCHMAX :  12.0  : Latent itch maximum, reported score capped at 4 (-)
KATX    :   0.08 : Autotaxin/LPA turnover (1/day)
KATXD   :   0.45 : Ductular-reaction contribution to autotaxin (-)
KSE     :   0.035: Central sensitisation turnover (1/day)
PBOMAX  :   0    : Placebo/observer itch effect, maximum (points)
PBOKT   :   0.05 : Placebo effect onset rate (1/day)
KSLP    :   0.10 : Sleep-disturbance turnover (1/day)

// =====================================================================
// MICELLAR FUNCTION, VITAMINS, GROWTH
// =====================================================================
CMC     :   2.2  : Critical micellar concentration (mM)
HMIC    :   2.0  : Hill exponent, micellar fat absorption (-)
FATMAX  :   0.95 : Maximum fat absorption coefficient (-)
FATNORM :   0.90 : Normal fat absorption coefficient (-)
KVD     :   0.02 : Vitamin D turnover (1/day)
VDIN    :  26    : Vitamin D input at normal absorption (ng/mL equivalent)
KVE     :   0.03 : Vitamin E turnover (1/day)
VEIN    :   1.0  : Vitamin E input at normal absorption (normalised)
KVK     :   0.25 : Vitamin K turnover (1/day)
VKIN    :   1.0  : Vitamin K input at normal absorption (normalised)
KVA     :   0.02 : Vitamin A turnover (1/day)
VAIN    :   1.0  : Vitamin A input at normal absorption (normalised)
KIG     :   0.05 : IGF-1 turnover (1/day)
KGHR    :   0.55 : GH resistance per unit relative hepatocyte BA (-)
KHT     :   0.0016: Height z-score adaptation rate (1/day)
HTZMAX  :   0.0  : Height z-score ceiling with normal nutrition (-)
KHF     :   3.4  : Height z-score cost of malabsorption (-)
KHI     :   1.6  : Height z-score cost of GH resistance (-)

// =====================================================================
// PORTAL HYPERTENSION
// =====================================================================
PP0     :   6.0  : Baseline portal pressure index (mmHg)
KPP     :   1.35 : Portal pressure per FIB^2 (mmHg)
PLT0    : 300    : Platelet count without hypersplenism (10^9/L)
PPH     :  13.0  : Portal pressure at half platelet loss (mmHg)
KPLT    :   0.05 : Platelet equilibration rate (1/day)

// =====================================================================
// HAZARDS AND COMPETING RISKS
// =====================================================================
HL0     :   0.00038: Liver decompensation hazard scale (1/day)
NEXP    :   1.60 : Exposure-hazard exponent  [GALA bilirubin strata give 1.5-1.85]
KHFIB   :   0.16 : Hazard multiplier per FIB^2 (-)
HQ0     :   0.000319: Transplant-for-pruritus hazard scale (1/day)
NQ      :   2.2  : Itch-hazard exponent (-)
HC0     :   0.00019: Cardiac hazard scale at birth (1/day)
KCA     :   0.55 : Cardiac hazard decay with age (1/yr)
HV0     :   0.0000125: Vascular / intracranial haemorrhage hazard (1/day)
FLFATAL :   0.12 : Fraction of liver events that are fatal rather than transplanted (-)
BAHREF  : 620    : Reference hepatocyte BA amount for exposure scaling (umol)
// the GALA native-liver-survival curve cannot be fitted by ANY single
// trajectory - its hazard is strongly front-loaded and this model does not.
// A gamma frailty of variance FRAILV fits all three GALA time points 65x
// better (SSE 5e-4 vs 3.3e-2).  FRAILV = 2.9 means the SD of individual
// hazard is 1.7x its mean -- the ALGS notorious phenotypic variance,
// measured.  Population curves must use the _POP outputs; the plain
// exp(-H) outputs are INDIVIDUAL-level and are not comparable to GALA.
FRAILV   :   2.935: Gamma frailty variance of the liver hazard (-)

$CMT @annotated
// ---- drug PK -------------------------------------------------------
DEPOT   : IBAT inhibitor, upper gut depot (ug)
DGUT    : IBAT inhibitor, mid gut (ug)
DILE    : IBAT inhibitor, terminal ileum - site of action (ug)
DSYS    : IBAT inhibitor, systemic (ug)
// ---- bile acid system ----------------------------------------------
BAHEP   : Hepatocyte bile acid - THE INJURY VARIABLE (umol)
BABIL   : Biliary tree and gallbladder bile acid (umol)
BADUO   : Duodenal / jejunal luminal bile acid - micellar window (umol)
BAILE   : Terminal ileal luminal bile acid (umol)
BAPOR   : Portal venous bile acid (umol)
BASYS   : Systemic bile acid - what sBA measures (umol)
BACOL   : Colonic bile acid (umol)
// ---- synthesis regulation ------------------------------------------
FGF19   : Ileal FGF19 signal, normalised (-)
CYP     : CYP7A1 activity, normalised (-)
// ---- liver structure -----------------------------------------------
DPR     : Bile duct : portal tract ratio (-)
DUCTR   : Ductular reaction activity, normalised (-)
HSC     : Activated hepatic stellate cell fraction (-)
FIB     : Fibrosis stage 0-4 (-)
// ---- biochemistry --------------------------------------------------
BILI    : Serum total bilirubin (mg/dL)
GGTX    : Serum GGT (U/L)
ALTX    : Serum ALT (U/L)
// ---- lipids --------------------------------------------------------
LPX     : Lipoprotein X, normalised (-)
XANTH   : Xanthoma burden, 0-4 (-)
// ---- pruritus ------------------------------------------------------
ATX     : Autotaxin / LPA axis, normalised (-)
SENS    : Central sensitisation, normalised (-)
SLEEPD  : Sleep disturbance score 0-4 (-)
// ---- nutrition and growth ------------------------------------------
VITD    : 25-OH vitamin D (ng/mL)
VITE    : Vitamin E : lipid ratio, normalised (-)
VITK    : Vitamin K status, normalised (-)
VITA    : Vitamin A status, normalised (-)
IGF1    : IGF-1, normalised (-)
HTZ     : Height z-score (-)
// ---- portal hypertension -------------------------------------------
PLT     : Platelet count (10^9/L)
// ---- cumulative exposure and hazards -------------------------------
CUMBA   : Cumulative hepatocyte BA exposure (baseline-days)
HZL     : Cumulative liver decompensation hazard (-)
HZQ     : Cumulative transplant-for-pruritus hazard (-)
HZC     : Cumulative cardiac hazard (-)
HZV     : Cumulative vascular hazard (-)
// ---- book-keeping --------------------------------------------------
FECBA   : Cumulative faecal bile acid loss (umol)
SYNCUM  : Cumulative bile acid synthesis (umol)
RETCUM  : Cumulative enterohepatic return to the hepatocyte (umol)

$MAIN
// ---- initial conditions --------------------------------------------
// $MAIN re-asserts these on EVERY run, which silently overrides any
// init() supplied from outside.  eq_init()/sim_eq() need that override
// suppressed, otherwise every scenario reports these hard-coded numbers
// as its baseline no matter what duct capacity it was given.
if (USEINIT < 0.5) {
  DPR_0    = DPR0;
  CYP_0    = 2.14;
  FGF19_0  = 0.0395;
  DUCTR_0  = 1.0;
  ATX_0    = 1.0;
  SENS_0   = 1.0;
  IGF1_0   = 0.55;
  PLT_0    = 250;
  VITD_0   = 14;
  VITE_0   = 0.45;
  VITK_0   = 0.60;
  VITA_0   = 0.50;
  HTZ_0    = -1.6;
  BILI_0   = 3.5;
  GGTX_0   = 340;
  ALTX_0   = 130;
  LPX_0    = 0.85;
  XANTH_0  = 0.8;
  SLEEPD_0 = 2.4;
  FIB_0    = 1.0;
  HSC_0    = 0.25;
// Bile acid compartments are started near their duct-limited steady
// state so that multi-year runs do not spend their first months
// equilibrating; sc00_equilibrate() verifies the residual drift.
  BAHEP_0  = 620;
  BABIL_0  = 32.1;
  BADUO_0  = 21.4;
  BAILE_0  = 32.1;
  BAPOR_0  = 10.5;
  BASYS_0  = 1066;
  BACOL_0  = 4.28;
}

$ODE
// =====================================================================
// 0. AGE AND STRUCTURAL CAPACITY
// =====================================================================
AGE   = AGE0 + SOLVERTIME / 365.25;

// Duct capacity as a multiple of baseline synthesis.  Conductance falls
// faster than duct NUMBER because remaining ducts are also narrower:
// GAMD ~ 2 (n_ducts x cross-sectional area).
JCAP  = JCAPN * pow(fmax(DPR, 1e-4) / DPRN, GAMD);

// Baseline (unregulated) synthesis flux, umol/day
SYN0  = SYNKG * WT;

// UDCA adds bile-acid-independent choleresis: it raises duct flux
// without changing duct NUMBER.
JDUCT = JCAP * SYN0 * (1.0 + UDCAF * UDCACH);

// =====================================================================
// 1. IBAT INHIBITOR PK  (minimally absorbed, lumen-acting)
// =====================================================================
// Dosing is delivered as a continuous input rather than 6570 daily bolus
// events; sc02b_bolus_equivalence() checks the two agree at 24 weeks.
ONTRT = (SOLVERTIME >= TSTART && SOLVERTIME < TSTOP) ? 1.0 : 0.0;
RIN   = ONTRT * DOSEKG * WT;                       // ug/day

dxdt_DEPOT = RIN - KTR1 * DEPOT;
dxdt_DGUT  = KTR1 * DEPOT - KTR2 * DGUT;
dxdt_DILE  = KTR2 * DGUT * (1.0 - FSYS) - KTR3 * DILE;
dxdt_DSYS  = KTR2 * DGUT * FSYS - KSYSD * DSYS;

CEFF  = POT * DILE;                                // potency-scaled ileal drug
IC50  = IC50KG * WT;
INH   = IMAXA * CEFF / (IC50 + CEFF);              // fractional ASBT blockade

// PEBD is a surgical version of the same lesion: it removes the ileal
// return entirely rather than inhibiting the transporter.
FREAB_D  = EPS0 * (1.0 - INH);
FREAB   = PEBDF > 0.5 ? EPS0 * PEBDEPS : FREAB_D;

// A sequestrant binds luminal bile acid, which makes it unavailable BOTH
// for micelle formation AND for ASBT uptake.  An IBAT inhibitor removes
// only the second.  That asymmetry is the whole of Axis 5.
FREAB_B  = FREAB * (1.0 - CHOLB);

// =====================================================================
// 2. SYNTHESIS AND ITS FEEDBACK
// =====================================================================
// FGF19 tracks ileal bile-acid UPTAKE, which is what ASBT blockade
// removes.  This is the long feedback arm and it is the reason IBAT
// inhibition cannot be more effective than PHI_EHC / sigma.
UPT   = FREAB_B * KIL * BAILE;                        // umol/day into enterocyte
UPT0  = EPS0 * FREF * SYNKG * WT;                  // NORMAL ileal uptake
dxdt_FGF19 = KFGF * (UPT / fmax(UPT0, 1e-6) - FGF19);

FL    = fmax(0.0, 1.0 - FGF19);                    // fractional FGF19 loss
SIGF  = 1.0 + (SIGMAX - 1.0) * FL / (KSIG + FL);   // saturating reserve

EXPO  = BAHEP / BAHREF;                            // relative hepatocyte load
FXRT  = (1.0 + KFXR) / (1.0 + KFXR * pow(fmax(EXPO, 1e-6), HFXR));

CYPT  = SIGF * FXRT;
dxdt_CYP = KCYP * (CYPT - CYP);

SYN   = SYN0 * CYP;

// =====================================================================
// 3. BILE ACID KINETICS
// =====================================================================
// Canalicular + ductal export is saturable with Vmax = JDUCT: this is
// the single structural fact the whole model turns on.
EXCR  = JDUCT * BAHEP / (KMH + BAHEP);

// Basolateral escape valve, adaptively induced when load exceeds duct.
KSP   = KSP0 * (1.0 + KIND * EXPO);
SPILL = KSP * BAHEP;

// Portal return and portal overflow
RET   = FPE * KPH * BAPOR;                         // to hepatocyte
POV   = (1.0 - FPE) * KPH * BAPOR;                 // to systemic
RECAP = FREUP * KSYSCL * BASYS;                    // systemic -> hepatocyte
RENAL = (1.0 - FREUP) * KSYSCL * BASYS;            // renal loss

dxdt_BAHEP = SYN + RET + RECAP - EXCR - SPILL;
dxdt_BABIL = EXCR - KBD * BABIL;
dxdt_BADUO = KBD * BABIL - KDI * BADUO;
dxdt_BAILE = KDI * BADUO - KIL * BAILE;
dxdt_BAPOR = UPT - KPH * BAPOR;
dxdt_BASYS = SPILL + POV - KSYSCL * BASYS;
dxdt_BACOL = (1.0 - FREAB_B) * KIL * BAILE - KFEC * BACOL;

dxdt_FECBA  = KFEC * BACOL;
dxdt_SYNCUM = SYN;
dxdt_RETCUM = RET;

// Serum bile acid concentration -- the trial endpoint
SBA   = BASYS / (VDBA * WT);

// PHI_EHC: the share of the hepatocyte load that came back from the gut.
// This is the ceiling on any IBAT inhibitor and it is a DUCT property.
PHIEHC = RET / fmax(SYN + RET, 1e-6);

// =====================================================================
// 4. MICELLAR WINDOW  (Axis 5: WHERE the drug acts)
// =====================================================================
// Duodenal concentration in mM.  A sequestrant binds here; an IBAT
// inhibitor acts two metres downstream and leaves this untouched.
CDUO  = KCDU * pow(fmax(EXCR, 1e-6), PCDU);
CDUOE = CDUO * (1.0 - CHOLB);
FATABS = FATMAX * pow(CDUOE, HMIC) / (pow(CMC, HMIC) + pow(CDUOE, HMIC));
FATREL = fmin(1.0, FATABS / FATNORM + MCTF * MCTGAIN);

// =====================================================================
// 5. LIVER STRUCTURE: DUCTULAR REPAIR AND FIBROSIS
// =====================================================================
INJ   = EXPO * (1.0 - UDCAF * UDCAHP);
dxdt_DUCTR = KDU * (INJ - DUCTR);

// Postnatal ductular repair -- why ALGS cholestasis often improves after
// age 5-10 -- opposed by fibrosis-driven duct loss.
dxdt_DPR = KDR * DUCTR * (DPRMAX - DPR) - KDLOSS * FIB * DPR;

APOP  = KAPO * pow(fmax(INJ, 1e-6), 1.4);
DRIVE = pow(fmax(INJ, 1e-6), HFIB) + KDRX * DUCTR;
dxdt_HSC = KHA * DRIVE * (HSCMAX - HSC) - KHD * HSC;
dxdt_FIB = KFF * HSC * (FIBMAX - FIB) - KFR * FIB;

// =====================================================================
// 6. BIOCHEMISTRY
// =====================================================================
// Bilirubin shares the ductal bottleneck: it is a DUCT-CAPACITY readout,
// not a bile-acid measure.  That is why it, and not sBA, predicts both
// natural history (GALA) and IBAT-inhibitor response.
BILT  = BILMIN + KBILC / fmax(JCAP, 0.05) + KBILF * FIB * FIB;
dxdt_BILI = KBI * (BILT - BILI);
dxdt_GGTX = KGG * (GGTMIN + KGGT * EXPO - GGTX);
dxdt_ALTX = KAL * (ALTMIN + KALT * APOP * (1.0 + RIFF * (RIFHEP - 1.0) + RIFF * 0.45) - ALTX);

// =====================================================================
// 7. LIPIDS AND XANTHOMA
// =====================================================================
dxdt_LPX = KLPX * (EXPO - LPX);
TCHOL = TC0 + KTC * LPX;
dxdt_XANTH = KXF * fmax(0.0, TCHOL - XTHR) * (4.0 - XANTH) / 4.0 - KXR * XANTH;

// =====================================================================
// 8. PRURITUS  (Axis 3)
// =====================================================================
// Autotaxin/LPA arm: driven by cholestasis and the ductular reaction,
// NOT proportional to serum bile acid.  Weight WIND is ASSUMED
// (FAILURE 3), constrained only by PEBD outcomes.
dxdt_ATX = KATX * (sqrt(fmax(EXPO, 1e-6)) * (1.0 + KATXD * DUCTR) - ATX);

// Weber-Fechner: itch codes the LOGARITHM of pruritogen, which is why a
// 41% fall in sBA does not give a 41% fall in itch in either direction.
DRV   = WBA * pow(SBA / SBAREF, QBA) + WIND * ATX + WSENS * SENS;
DRVM  = DRV * (1.0 - RIFF * RIFITCH) * (1.0 - NALF * NALITCH);
ITCH  = fmin(4.0, ITCHMAX * DRVM / (K50I + DRVM));

dxdt_SENS   = KSE * (ITCH / 2.9 - SENS);
dxdt_SLEEPD = KSLP * (ITCH * 0.85 - SLEEPD);

// =====================================================================
// 9. VITAMINS, IGF-1, GROWTH
// =====================================================================
dxdt_VITD = KVD * (VDIN * FATREL - VITD);
dxdt_VITE = KVE * (VEIN * FATREL - VITE);
dxdt_VITK = KVK * (VKIN * FATREL - VITK);
dxdt_VITA = KVA * (VAIN * FATREL - VITA);

dxdt_IGF1 = KIG * (1.0 / (1.0 + KGHR * EXPO) - IGF1);

HTGT = HTZMAX - KHF * (1.0 - FATREL) - KHI * (1.0 - IGF1);
dxdt_HTZ = KHT * (HTGT - HTZ);

// =====================================================================
// 10. PORTAL HYPERTENSION
// =====================================================================
PORTP = PP0 + KPP * FIB * FIB;
PLTT  = PLT0 / (1.0 + pow(PORTP / PPH, 4.0));
dxdt_PLT = KPLT * (PLTT - PLT);

// =====================================================================
// 11. EXPOSURE AND COMPETING HAZARDS  (Axes 4 and 6)
// =====================================================================
dxdt_CUMBA = EXPO;

HL = HL0 * pow(fmax(EXPO, 1e-6), NEXP) * (1.0 + KHFIB * FIB * FIB);
HQ = HQ0 * pow(fmax(ITCH / 2.9, 1e-6), NQ);
HC = HC0 * exp(-KCA * AGE);
HV = HV0;

dxdt_HZL = HL;
dxdt_HZQ = HQ;
dxdt_HZC = HC;
dxdt_HZV = HV;

$TABLE
// Every captured quantity is RECOMPUTED here from compartment state.
// Locals left behind by the last $ODE derivative evaluation do not
// necessarily correspond to the output time, and silently reporting
// those would corrupt every trial replication in this file.
JCAP  = JCAPN * pow(fmax(DPR, 1e-4) / DPRN, GAMD);
SYN0  = SYNKG * WT;
JDUCT = JCAP * SYN0 * (1.0 + UDCAF * UDCACH);
SYN   = SYN0 * CYP;
EXPO  = BAHEP / BAHREF;

CEFF  = POT * DILE;
IC50  = IC50KG * WT;
INH   = IMAXA * CEFF / (IC50 + CEFF);
FREAB   = PEBDF > 0.5 ? EPS0 * PEBDEPS : EPS0 * (1.0 - INH);
FREAB_B  = FREAB * (1.0 - CHOLB);

EXCR  = JDUCT * BAHEP / (KMH + BAHEP);
KSP   = KSP0 * (1.0 + KIND * EXPO);
SPILL = KSP * BAHEP;
RET   = FPE * KPH * BAPOR;
RECAP = FREUP * KSYSCL * BASYS;
PHIEHC = RET / fmax(SYN + RET, 1e-6);

FL    = fmax(0.0, 1.0 - FGF19);
SIGF  = 1.0 + (SIGMAX - 1.0) * FL / (KSIG + FL);

CDUO  = KCDU * pow(fmax(EXCR, 1e-6), PCDU);
CDUOE = CDUO * (1.0 - CHOLB);
FATABS = FATMAX * pow(CDUOE, HMIC) / (pow(CMC, HMIC) + pow(CDUOE, HMIC));
FATREL = fmin(1.0, FATABS / FATNORM + MCTF * MCTGAIN);

DRV   = WBA * pow((BASYS / (VDBA * WT)) / SBAREF, QBA) + WIND * ATX + WSENS * SENS;
DRVM  = DRV * (1.0 - RIFF * RIFITCH) * (1.0 - NALF * NALITCH);
ITCH  = fmin(4.0, ITCHMAX * DRVM / (K50I + DRVM));

PORTP = PP0 + KPP * FIB * FIB;

// Trial-reported quantities
double SBA_OUT   = BASYS / (VDBA * WT);
double AGE_OUT   = AGE0 + TIME / 365.25;

// Placebo / observer component, applied to the REPORTED itch only.
// Keeping it here rather than inside the physiology is what makes the
// drug-attributable effect separable (Axis 3).
double PBO       = PBOMAX * (1.0 - exp(-PBOKT * TIME));
double ITCH_OBS  = fmax(0.0, ITCH - PBO);

// Event-free survival = free of portal-hypertensive event, biliary
// diversion, transplant or death (the maralixibat-vs-GALA definition).
double EFS       = exp(-(HZL + HZQ));
double NLS       = exp(-HZL);
// Population-level curves under gamma frailty -- the only ones
// comparable to a cohort study such as GALA.
double NLS_POP   = pow(1.0 + FRAILV * HZL, -1.0 / FRAILV);
double EFS_POP   = pow(1.0 + FRAILV * (HZL + HZQ), -1.0 / FRAILV);
double SURV      = exp(-(FLFATAL * HZL + HZC + HZV));
double SURVLIV   = exp(-(FLFATAL * HZL));

double TCHOL_OUT = TC0 + KTC * LPX;
double FATABS_O  = FATABS;
double FATREL_O  = FATREL;
double CDUO_OUT  = CDUO;
double JCAP_OUT  = JCAP;
double PHI_OUT   = PHIEHC;
double INH_OUT   = INH;
double SIG_OUT   = SIGF;
double EXPO_OUT  = BAHEP / BAHREF;
double PORTP_O   = PORTP;
double INR       = 1.0 + 1.6 * (1.0 - VITK) + 0.35 * FIB * FIB / 16.0;

$CAPTURE
SBA_OUT ITCH ITCH_OBS AGE_OUT EFS NLS NLS_POP EFS_POP SURV SURVLIV TCHOL_OUT FATABS_O
FATREL_O CDUO_OUT JCAP_OUT PHI_OUT INH_OUT SIG_OUT EXPO_OUT PORTP_O INR
'

# ---------------------------------------------------------------------
#  BUILD
# ---------------------------------------------------------------------
mod_algs <- mcode("algs", algs_code)

yr <- function(x) x * 365.25

#' Core simulation wrapper.
#' @param end simulation length in DAYS
sim <- function(end = yr(15), delta = 7, events = NULL, ...) {
  p <- list(...)
  m <- mod_algs
  if (length(p)) m <- param(m, p)
  out <- if (is.null(events)) {
    mrgsim(m, end = end, delta = delta, atol = 1e-8, rtol = 1e-6)
  } else {
    mrgsim(m, events = events, end = end, delta = delta,
           atol = 1e-8, rtol = 1e-6)
  }
  as.data.frame(out)
}

#' Equilibrate the FAST variables (bile acids, bilirubin, itch, vitamins)
#' for a given structural state, holding the SLOW structural ones (DPR,
#' fibrosis) frozen, and return the resulting initial-condition list.
#'
#' Without this, every scenario that varies DPR0 silently reports the
#' hard-coded initial conditions as its "baseline" -- so baseline serum
#' bile acid and bilirubin came out IDENTICAL for every duct capacity,
#' which is precisely the dependence the duct-capacity axis is about.
eq_init <- function(burn = 400, ...) {
  p <- list(...)
  frozen <- list(KDR = 0, KDLOSS = 0, KFF = 0, KFR = 0,
                 DOSEKG = 0, TSTART = 0, TSTOP = 0,
                 PEBDF = 0, CHOLB = 0, UDCAF = 0, RIFF = 0, NALF = 0, MCTF = 0)
  p <- modifyList(frozen, p)
  m <- param(mod_algs, p)
  out <- as.data.frame(mrgsim(m, end = burn, delta = burn,
                              atol = 1e-8, rtol = 1e-6))
  cmts <- names(mod_algs@init@data)
  last <- out[nrow(out), cmts, drop = FALSE]
  # cumulative book-keeping compartments must restart at zero
  for (z in c("CUMBA", "HZL", "HZQ", "HZC", "HZV",
              "FECBA", "SYNCUM", "RETCUM")) last[[z]] <- 0
  as.list(last)
}

#' Simulate from an equilibrated baseline for the given structural state.
sim_eq <- function(end = yr(15), delta = 7, events = NULL, burn = 400, ...) {
  p <- list(...)
  # structural parameters define the baseline; therapy does not
  struct <- p[intersect(names(p),
              c("WT", "AGE0", "DPR0", "DPRN", "JCAPN", "GAMD", "SYNKG",
                "NEXP", "HL0", "HQ0", "IMAXA", "KBILC", "KBILF", "BILMIN",
                "EPS0", "SIGMAX", "KSIG"))]
  init0 <- do.call(eq_init, c(list(burn = burn), struct))
  m <- param(mod_algs, c(p, list(USEINIT = 1)))
  m <- update(m, init = init0)
  out <- if (is.null(events)) {
    mrgsim(m, end = end, delta = delta, atol = 1e-8, rtol = 1e-6)
  } else {
    mrgsim(m, events = events, end = end, delta = delta,
           atol = 1e-8, rtol = 1e-6)
  }
  as.data.frame(out)
}

#' Value of a column at (or nearest below) a given day
at_day <- function(d, col, day) {
  i <- which.min(abs(d$time - day))
  d[[col]][i]
}

#' Mean of a column over a day window (trial endpoints are window means)
win_mean <- function(d, col, lo, hi) {
  mean(d[[col]][d$time >= lo & d$time <= hi], na.rm = TRUE)
}

# Drug shorthands -----------------------------------------------------
#   Odevixibat is the potency reference (POT = 1.0).  Maralixibat needs
#   ~3x the dose for the same ileal effect, hence POT = 0.28 at
#   380 ug/kg/day.
ODEV <- list(DOSEKG = 120, POT = 1.00)
MRX  <- list(DOSEKG = 380, POT = 0.28)

# =====================================================================
#  SCENARIO 0 — EQUILIBRATION CHECK
#    The bile-acid compartments are initialised near their duct-limited
#    steady state.  This reports the residual drift, because a model
#    that silently equilibrates for six months would make every 24-week
#    trial replication below meaningless.
# =====================================================================
sc00_equilibrate <- function() {
  d <- sim(end = 400, delta = 1)
  data.frame(
    quantity = c("sBA (umol/L)", "BA_HEP (umol)", "PHI_EHC", "JCAP", "BILI"),
    day7     = c(at_day(d, "SBA_OUT", 7),  at_day(d, "BAHEP", 7),
                 at_day(d, "PHI_OUT", 7),  at_day(d, "JCAP_OUT", 7),
                 at_day(d, "BILI", 7)),
    day400   = c(at_day(d, "SBA_OUT", 400), at_day(d, "BAHEP", 400),
                 at_day(d, "PHI_OUT", 400), at_day(d, "JCAP_OUT", 400),
                 at_day(d, "BILI", 400)),
    stringsAsFactors = FALSE
  ) %>% mutate(pct_drift = 100 * (day400 - day7) / day7)
}

# =====================================================================
#  SCENARIO 1 — NATURAL HISTORY, age 1 -> 18
#    Calibration target (GALA, n=1433): native liver survival
#    5 yr 66.8% / 10 yr 54.4% / 18 yr 40.3%; adverse liver event
#    (portal hypertension, transplant or death) 51.5% by 10 yr,
#    66.0% by 18 yr.  See FAILURE 5 on reading a cohort curve as a
#    patient trajectory.
# =====================================================================
sc01_natural_history <- function() {
  sim(end = yr(17), delta = 30.4375, AGE0 = 1)
}

sc01_gala_table <- function() {
  d <- sc01_natural_history()
  ages <- c(5, 10, 18)
  data.frame(
    age_yr        = ages,
    # POPULATION curves (gamma frailty). The individual-level exp(-H)
    # outputs are NOT comparable to a cohort study; using them here was a
    # defect that made the natural history look catastrophically wrong.
    model_NLS     = round(100 * sapply(ages, function(a)
                      at_day(d, "NLS_POP", yr(a - 1))), 1),
    GALA_NLS      = c(66.8, 54.4, 40.3),
    model_EFS     = round(100 * sapply(ages, function(a)
                      at_day(d, "EFS_POP", yr(a - 1))), 1),
    GALA_eventfree = c(NA, 100 - 51.5, 100 - 66.0)
  )
}

# =====================================================================
#  SCENARIO 2 — ICONIC REPLICATION (maralixibat 380 ug/kg/day)
#    Targets: week-48 sBA change -96 umol/L; ItchRO(Obs) change
#    -1.6 points (both UNCONTROLLED, single-arm — see Axis 3).
# =====================================================================
sc02_iconic <- function() {
  d <- sim(end = yr(1.2), delta = 1,
           DOSEKG = MRX$DOSEKG, POT = MRX$POT,
           TSTART = 0, TSTOP = yr(1.2), PBOMAX = 0.7)
  data.frame(
    endpoint = c("sBA (umol/L)", "ItchRO(Obs) reported", "ItchRO drug-only",
                 "xanthoma", "height z"),
    baseline = c(at_day(d, "SBA_OUT", 0), at_day(d, "ITCH_OBS", 0),
                 at_day(d, "ITCH", 0), at_day(d, "XANTH", 0),
                 at_day(d, "HTZ", 0)),
    week48   = c(at_day(d, "SBA_OUT", 336), at_day(d, "ITCH_OBS", 336),
                 at_day(d, "ITCH", 336), at_day(d, "XANTH", 336),
                 at_day(d, "HTZ", 336)),
    ICONIC   = c(-96, -1.6, NA, NA, NA)
  ) %>% mutate(change = week48 - baseline)
}

# Bolus-dosing equivalence check for the continuous-input approximation
sc02b_bolus_equivalence <- function() {
  # Sampled at delta = 0.05 d: with KTR3 = 6/day the ileal drug amount is
  # almost fully cleared between daily doses, so sampling on integer days
  # lands on one phase of the cycle and reports a "mean" blockade that is
  # a sampling artefact rather than a time average.
  cont <- sim(end = 168, delta = 0.05, DOSEKG = ODEV$DOSEKG, POT = ODEV$POT,
              TSTART = 0, TSTOP = 168)
  amt  <- ODEV$DOSEKG * 15
  ev   <- ev(amt = amt, cmt = "DEPOT", ii = 1, addl = 167)
  bol  <- sim(end = 168, delta = 0.05, events = ev, POT = ODEV$POT,
              DOSEKG = 0, TSTART = 0, TSTOP = 0)
  data.frame(
    quantity    = c("sBA at day 168",
                    "ASBT blockade, mean over days 161-168"),
    continuous  = c(at_day(cont, "SBA_OUT", 168),
                    win_mean(cont, "INH_OUT", 161, 168)),
    daily_bolus = c(at_day(bol, "SBA_OUT", 168),
                    win_mean(bol, "INH_OUT", 161, 168))
  ) %>% mutate(pct_diff = 100 * (daily_bolus - continuous) / continuous)
}

# =====================================================================
#  SCENARIO 3 — ASSERT REPLICATION, BOTH ARMS (odevixibat 120 ug/kg/day)
#    This is the model's primary calibration AND the proof of Axis 3:
#    the placebo arm's itch falls 0.8 points while its bile acids RISE.
#    Targets: sBA 237 -> 149 (drug) vs 246 -> 271 (placebo);
#             scratching 2.8 -> 1.1 (drug) vs 3.0 -> 2.2 (placebo).
# =====================================================================
sc03_assert <- function() {
  drug <- sim(end = 180, delta = 1, DOSEKG = ODEV$DOSEKG, POT = ODEV$POT,
              TSTART = 0, TSTOP = 180, PBOMAX = 0.8)
  pbo  <- sim(end = 180, delta = 1, DOSEKG = 0, TSTART = 0, TSTOP = 0,
              PBOMAX = 0.8)
  wk <- function(d, col) win_mean(d, col, 147, 168)   # weeks 21-24
  data.frame(
    arm        = c("odevixibat", "placebo", "difference"),
    sBA_base   = c(at_day(drug, "SBA_OUT", 0), at_day(pbo, "SBA_OUT", 0), NA),
    sBA_wk24   = c(wk(drug, "SBA_OUT"), wk(pbo, "SBA_OUT"), NA),
    sBA_change = c(wk(drug, "SBA_OUT") - at_day(drug, "SBA_OUT", 0),
                   wk(pbo,  "SBA_OUT") - at_day(pbo,  "SBA_OUT", 0), NA),
    itch_base  = c(at_day(drug, "ITCH_OBS", 0), at_day(pbo, "ITCH_OBS", 0), NA),
    itch_wk24  = c(wk(drug, "ITCH_OBS"), wk(pbo, "ITCH_OBS"), NA),
    itch_change = c(wk(drug, "ITCH_OBS") - at_day(drug, "ITCH_OBS", 0),
                    wk(pbo,  "ITCH_OBS") - at_day(pbo,  "ITCH_OBS", 0), NA)
  ) %>%
    mutate(
      sBA_change  = ifelse(arm == "difference",
                           sBA_change[1] - sBA_change[2], sBA_change),
      itch_change = ifelse(arm == "difference",
                           itch_change[1] - itch_change[2], itch_change)
    )
}

# ASSERT observed, for side-by-side printing
ASSERT_OBS <- data.frame(
  arm         = c("odevixibat", "placebo", "difference"),
  sBA_base    = c(237, 246, NA),
  sBA_wk24    = c(149, 271, NA),
  sBA_change  = c(-90, 22, -113),
  itch_base   = c(2.8, 3.0, NA),
  itch_wk24   = c(1.1, 2.2, NA),
  itch_change = c(-1.7, -0.8, -0.9)
)

# =====================================================================
#  SCENARIO 4 — THE ITCH-PER-BILE-ACID SLOPE (Axis 3, quantified)
#    Reports the controlled and uncontrolled slopes the model produces
#    and the published ones.  A model that reproduced only the
#    uncontrolled slope would over-attribute itch to bile acids ~2x.
# =====================================================================
sc04_slopes <- function() {
  a <- sc03_assert()
  ctrl   <- abs(a$itch_change[3] / a$sBA_change[3])
  unctrl <- abs(a$itch_change[1] / a$sBA_change[1])
  data.frame(
    slope        = c("placebo-controlled", "single-arm (ASSERT)",
                     "single-arm (ICONIC)"),
    model        = c(ctrl, unctrl, NA),
    published    = c(0.9 / 113, 1.7 / 90, 1.6 / 96),
    stringsAsFactors = FALSE
  ) %>% mutate(ratio_to_controlled = published / published[1])
}

# =====================================================================
#  SCENARIO 5 — DUCT-CAPACITY SWEEP: THE RESPONDER THRESHOLD (Axis 2)
#    The central prediction.  sBA response to a fixed IBAT-inhibitor
#    dose is swept against baseline duct capacity.  Below a threshold
#    the response INVERTS, because CYP7A1 derepression outruns the
#    intercepted return and the new bile acid has no duct to leave by.
# =====================================================================
sc05_duct_sweep <- function(dpr = seq(0.10, 0.45, by = 0.01)) {
  do.call(rbind, lapply(dpr, function(x) {
    base <- sim_eq(end = 180, delta = 30, DPR0 = x, DOSEKG = 0)
    drug <- sim_eq(end = 180, delta = 1, DPR0 = x,
                   DOSEKG = ODEV$DOSEKG, POT = ODEV$POT,
                   TSTART = 0, TSTOP = 180)
    b <- at_day(base, "SBA_OUT", 0)
    t <- win_mean(drug, "SBA_OUT", 147, 168)
    data.frame(
      DPR0      = x,
      JCAP      = at_day(base, "JCAP_OUT", 0),
      BILI      = at_day(base, "BILI", 0),
      PHI_EHC   = at_day(base, "PHI_OUT", 0),
      sBA_base  = b,
      sBA_wk24  = t,
      pct_change = 100 * (t - b) / b
    )
  }))
}

#' Locate the duct capacity / bilirubin at which the sBA response inverts
sc05_threshold <- function(sweep = sc05_duct_sweep()) {
  s <- sweep[order(sweep$JCAP), ]
  cross <- function(y) {
    i <- which(diff(sign(y)) != 0)
    if (!length(i)) return(NA_real_)
    i <- i[1]
    s$BILI[i] + (s$BILI[i + 1] - s$BILI[i]) * (0 - y[i]) / (y[i + 1] - y[i])
  }
  crossj <- function(y) {
    i <- which(diff(sign(y)) != 0)
    if (!length(i)) return(NA_real_)
    i <- i[1]
    s$JCAP[i] + (s$JCAP[i + 1] - s$JCAP[i]) * (0 - y[i]) / (y[i + 1] - y[i])
  }
  data.frame(
    boundary = c("response < 15% (drug futile)",
                 "response < 30% (below a clinically useful effect)",
                 "response < 50%",
                 "response reaches -85% (effectively normalised)"),
    bilirubin_mgdL = c(cross(s$pct_change + 15), cross(s$pct_change + 30),
                       cross(s$pct_change + 50), cross(s$pct_change + 85)),
    JCAP = c(crossj(s$pct_change + 15), crossj(s$pct_change + 30),
             crossj(s$pct_change + 50), crossj(s$pct_change + 85)),
    GALA_reference = c(NA, 5.0, NA, NA)
  )
}

# =====================================================================
#  SCENARIO 6 — PHI_EHC AND ITS IDENTIFIABILITY (Axis 1, FAILURE 1)
#    Sweeps the assumed achieved ASBT blockade I and reports how JCAP
#    and PHI_EHC move.  PHI_EHC is the robust quantity; JCAP is not.
# =====================================================================
sc06_phi_identifiability <- function(imax = seq(0.50, 0.92, by = 0.07)) {
  target <- 149 / 271           # ASSERT drug-arm / placebo-arm sBA ratio
  do.call(rbind, lapply(imax, function(im) {
    # For each ASSUMED maximum blockade, refit duct capacity so the model
    # still reproduces ASSERT exactly, then read off PHI_EHC.
    err <- function(dpr) {
      b <- sim_eq(end = 1, delta = 1, DPR0 = dpr, IMAXA = im)
      d <- sim_eq(end = 180, delta = 3, DPR0 = dpr, IMAXA = im,
                  DOSEKG = ODEV$DOSEKG, POT = ODEV$POT,
                  TSTART = 0, TSTOP = 180)
      win_mean(d, "SBA_OUT", 147, 168) / at_day(b, "SBA_OUT", 0) - target
    }
    r <- try(uniroot(err, c(0.14, 0.55), tol = 1e-4), silent = TRUE)
    if (inherits(r, "try-error")) return(NULL)
    b <- sim_eq(end = 1, delta = 1, DPR0 = r$root, IMAXA = im)
    d <- sim_eq(end = 180, delta = 3, DPR0 = r$root, IMAXA = im,
                DOSEKG = ODEV$DOSEKG, POT = ODEV$POT,
                TSTART = 0, TSTOP = 180)
    data.frame(
      IMAXA_assumed = im,
      blockade      = win_mean(d, "INH_OUT", 147, 168),
      DPR_refitted  = r$root,
      JCAP_refitted = at_day(b, "JCAP_OUT", 0),
      PHI_EHC       = at_day(b, "PHI_OUT", 0),
      BILI_baseline = at_day(b, "BILI", 0)
    )
  }))
}

# =====================================================================
#  SCENARIO 6b — THE BILIRUBIN ANCHOR (FAILURE 2, quantified)
#    The duct-capacity -> bilirubin map is anchored on ONE assumed number:
#    the mean baseline total bilirubin of the ASSERT population, taken as
#    3.5 mg/dL.  This sweeps that anchor and reports how far the
#    responder boundary moves, because the convergence with the GALA
#    cut-points is only worth anything if it survives the anchor.
# =====================================================================
sc06b_anchor_sensitivity <- function(anchor = c(2.5, 3.0, 3.5, 4.0, 4.5)) {
  do.call(rbind, lapply(anchor, function(a) {
    kb <- (a - 0.40 - 0.28) * 2.10      # BILMIN + KBILF*FIB^2 at JCAP 2.10
    sw <- sc05_duct_sweep_k(kb)
    th <- sc05_threshold(sw)
    data.frame(
      assumed_ASSERT_bilirubin = a,
      KBILC                    = kb,
      bili_at_30pct_response   = th$bilirubin_mgdL[2],
      bili_at_15pct_response   = th$bilirubin_mgdL[1],
      GALA_low_cut             = 5.0,
      GALA_high_cut            = 10.0
    )
  }))
}

sc05_duct_sweep_k <- function(kbilc, dpr = seq(0.10, 0.45, by = 0.02)) {
  do.call(rbind, lapply(dpr, function(x) {
    base <- sim_eq(end = 180, delta = 30, DPR0 = x, KBILC = kbilc, DOSEKG = 0)
    drug <- sim_eq(end = 180, delta = 1, DPR0 = x, KBILC = kbilc,
                   DOSEKG = ODEV$DOSEKG, POT = ODEV$POT,
                   TSTART = 0, TSTOP = 180)
    b <- at_day(base, "SBA_OUT", 0)
    t <- win_mean(drug, "SBA_OUT", 147, 168)
    data.frame(DPR0 = x, JCAP = at_day(base, "JCAP_OUT", 0),
               BILI = at_day(base, "BILI", 0),
               PHI_EHC = at_day(base, "PHI_OUT", 0),
               sBA_base = b, sBA_wk24 = t,
               pct_change = 100 * (t - b) / b)
  }))
}

# =====================================================================
#  SCENARIO 7 — SURGICAL COMPARATOR: PEBD
#    Complete interruption (I -> 1) rather than partial inhibition.
#    The model predicts near-complete sBA normalisation in a patient
#    with adequate duct capacity and NO benefit below the threshold —
#    the same boundary as scenario 5, which is the point.
# =====================================================================
sc07_pebd <- function(dpr = c(0.32, 0.236, 0.16)) {
  do.call(rbind, lapply(dpr, function(x) {
    b <- sim_eq(end = 1,   delta = 1, DPR0 = x, DOSEKG = 0)
    p <- sim_eq(end = 365, delta = 7, DPR0 = x, PEBDF = 1)
    m <- sim_eq(end = 365, delta = 7, DPR0 = x,
                DOSEKG = ODEV$DOSEKG, POT = ODEV$POT,
                TSTART = 0, TSTOP = 365)
    data.frame(
      DPR0        = x,
      BILI_base   = at_day(b, "BILI", 0),
      sBA_base    = at_day(b, "SBA_OUT", 0),
      sBA_PEBD    = at_day(p, "SBA_OUT", 365),
      sBA_IBATi   = at_day(m, "SBA_OUT", 365),
      itch_base   = at_day(b, "ITCH", 0),
      itch_PEBD   = at_day(p, "ITCH", 365),
      itch_IBATi  = at_day(m, "ITCH", 365)
    )
  }))
}

# =====================================================================
#  SCENARIO 8 — DRUG PANEL AT 24 WEEKS
#    Every agent in the map, on the same patient, on the same endpoints.
# =====================================================================
sc08_drug_panel <- function() {
  arms <- list(
    "untreated"                = list(),
    "UDCA 20 mg/kg/d"          = list(UDCAF = 1),
    "cholestyramine"           = list(CHOLB = 0.55),
    "rifampicin"               = list(RIFF = 1),
    "naltrexone"               = list(NALF = 1),
    "maralixibat 380 ug/kg/d"  = c(MRX,  list(TSTART = 0, TSTOP = 180)),
    "odevixibat 120 ug/kg/d"   = c(ODEV, list(TSTART = 0, TSTOP = 180)),
    "PEBD"                     = list(PEBDF = 1),
    "odevixibat + MCT formula" = c(ODEV, list(TSTART = 0, TSTOP = 180, MCTF = 1))
  )
  do.call(rbind, lapply(names(arms), function(nm) {
    d <- do.call(sim, c(list(end = 180, delta = 7), arms[[nm]]))
    data.frame(
      arm        = nm,
      sBA        = at_day(d, "SBA_OUT", 180),
      itch       = at_day(d, "ITCH", 180),
      bilirubin  = at_day(d, "BILI", 180),
      ALT        = at_day(d, "ALTX", 180),
      fat_abs    = at_day(d, "FATABS_O", 180),
      vitD       = at_day(d, "VITD", 180),
      chol       = at_day(d, "TCHOL_OUT", 180),
      stringsAsFactors = FALSE
    )
  }))
}

# =====================================================================
#  SCENARIO 9 — SITE OF ACTION AND THE MICELLAR WINDOW (Axis 5)
#    Both drugs remove bile acids; only the proximal-acting one starves
#    the child.  Also reports the boundary case the contrast fails in.
# =====================================================================
sc09_micellar <- function(dpr = c(0.236, 0.60)) {
  do.call(rbind, lapply(dpr, function(x) {
    arms <- list(
      untreated      = list(),
      odevixibat     = c(ODEV, list(TSTART = 0, TSTOP = yr(2))),
      cholestyramine = list(CHOLB = 0.55)
    )
    do.call(rbind, lapply(names(arms), function(nm) {
      d <- do.call(sim_eq, c(list(end = yr(2), delta = 14, DPR0 = x), arms[[nm]]))
      data.frame(
        DPR0     = x,
        paucity  = ifelse(x < 0.4, "severe (duct-saturated)",
                                   "mild (duct NOT saturated)"),
        arm      = nm,
        duod_mM  = at_day(d, "CDUO_OUT", yr(2)),
        fat_abs  = at_day(d, "FATABS_O", yr(2)),
        vitD     = at_day(d, "VITD", yr(2)),
        vitE     = at_day(d, "VITE", yr(2)),
        INR      = at_day(d, "INR", yr(2)),
        height_z = at_day(d, "HTZ", yr(2)),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

# =====================================================================
#  SCENARIO 10 — EARLY VS LATE START: THE INTEGRAL ARGUMENT
#    Fibrosis is driven by CUMULATIVE hepatocyte bile-acid exposure, so
#    the same drug bought at a different age buys a different amount of
#    liver.  Run to age 18 in every arm so the comparison is fair.
# =====================================================================
sc10_timing <- function(start_ages = c(1, 3, 5, 8, 12), to_age = 18) {
  do.call(rbind, lapply(start_ages, function(a) {
    d <- sim(end = yr(to_age - 1), delta = 30.4375, AGE0 = 1,
             DOSEKG = MRX$DOSEKG, POT = MRX$POT,
             TSTART = yr(a - 1), TSTOP = yr(to_age - 1))
    nat <- sc01_natural_history()
    data.frame(
      start_age     = a,
      CUMBA_treated = at_day(d, "CUMBA", yr(to_age - 1)),
      CUMBA_natural = at_day(nat, "CUMBA", yr(to_age - 1)),
      FIB_treated   = at_day(d, "FIB", yr(to_age - 1)),
      FIB_natural   = at_day(nat, "FIB", yr(to_age - 1)),
      EFS_treated   = at_day(d, "EFS_POP", yr(to_age - 1)),
      EFS_natural   = at_day(nat, "EFS_POP", yr(to_age - 1)),
      heightz_treated = at_day(d, "HTZ", yr(to_age - 1)),
      heightz_natural = at_day(nat, "HTZ", yr(to_age - 1))
    )
  })) %>%
    mutate(exposure_averted_pct = 100 * (1 - CUMBA_treated / CUMBA_natural))
}

# =====================================================================
#  SCENARIO 11 — THE EVENT-FREE-SURVIVAL TEST (Axis 4)
#    Six years of maralixibat against the model's own natural history,
#    scored the way the published maralixibat-vs-GALA comparison scored
#    it (HR 0.305, 95% CI 0.189-0.491).  The model is shipped with
#    NEXP = 1.6 -- GALA's OWN exponent -- so it is EXPECTED to
#    under-predict.  That gap is the finding, not a fitting failure.
# =====================================================================
sc11_efs <- function(nexp = 1.60, years = 6) {
  nat <- sim(end = yr(years), delta = 30.4375, AGE0 = 3, NEXP = nexp)
  trt <- sim(end = yr(years), delta = 30.4375, AGE0 = 3, NEXP = nexp,
             DOSEKG = MRX$DOSEKG, POT = MRX$POT,
             TSTART = 0, TSTOP = yr(years))
  hz_n <- at_day(nat, "HZL", yr(years)) + at_day(nat, "HZQ", yr(years))
  hz_t <- at_day(trt, "HZL", yr(years)) + at_day(trt, "HZQ", yr(years))
  data.frame(
    NEXP            = nexp,
    EFS_natural     = at_day(nat, "EFS_POP", yr(years)),
    EFS_maralixibat = at_day(trt, "EFS_POP", yr(years)),
    HR_model        = hz_t / hz_n,
    HR_published    = 0.305,
    sBA_reduction_pct = 100 * (1 - at_day(trt, "SBA_OUT", yr(years)) /
                                   at_day(nat, "SBA_OUT", yr(years)))
  )
}

# =====================================================================
#  SCENARIO 12 — THE EXPONENT CONTRADICTION (Axis 4, the sharp result)
#    Sweeps NEXP and reports, for each value, BOTH what it implies for
#    the maralixibat hazard ratio AND what it back-predicts for GALA's
#    own bilirubin strata (observed 4.8x and 15.6x).  No single exponent
#    fits both.
# =====================================================================
sc12_exponent_conflict <- function(nexp = c(1.5, 1.6, 1.85, 2.1, 2.4, 2.9)) {
  # GALA stratum exposure ratios, from total-bilirubin midpoints
  r_mid <- 7.0 / 3.0     # 5-10 vs <5 mg/dL
  r_hi  <- 14.0 / 3.0    # >=10 vs <5 mg/dL
  do.call(rbind, lapply(nexp, function(n) {
    e <- sc11_efs(nexp = n)
    data.frame(
      NEXP              = n,
      HR_model_MRX      = e$HR_model,
      HR_published_MRX  = 0.305,
      GALA_mid_pred     = r_mid^n,
      GALA_mid_obs      = 4.8,
      GALA_high_pred    = r_hi^n,
      GALA_high_obs     = 15.6
    )
  })) %>%
    mutate(
      MRX_error       = HR_model_MRX / HR_published_MRX,
      GALA_high_error = GALA_high_pred / GALA_high_obs
    )
}

# =====================================================================
#  SCENARIO 13 — THE POPULATION CEILING (Axis 6)
#    Cardiac and vascular hazard are Notch-dose diseases no cholestasis
#    drug touches.  All-cause benefit is bounded by the liver's share.
# =====================================================================
sc13_competing_risk <- function(years = 18) {
  nat <- sim(end = yr(years), delta = 30.4375, AGE0 = 1)
  trt <- sim(end = yr(years), delta = 30.4375, AGE0 = 1,
             DOSEKG = MRX$DOSEKG, POT = MRX$POT,
             TSTART = 0, TSTOP = yr(years))
  hzl_n <- at_day(nat, "HZL", yr(years)) * 0.12
  hzl_t <- at_day(trt, "HZL", yr(years)) * 0.12
  hzc   <- at_day(nat, "HZC", yr(years))
  hzv   <- at_day(nat, "HZV", yr(years))
  data.frame(
    component = c("liver (fatal fraction)", "cardiac", "vascular", "TOTAL"),
    natural   = c(hzl_n, hzc, hzv, hzl_n + hzc + hzv),
    treated   = c(hzl_t, hzc, hzv, hzl_t + hzc + hzv)
  ) %>%
    mutate(
      hazard_share_pct = 100 * natural / natural[4],
      HR               = treated / natural
    )
}

# =====================================================================
#  SCENARIO 14 — NON-RESPONDER PHENOTYPE, FULL TRAJECTORY
#    A child whose duct capacity sits below the scenario-5 threshold,
#    treated anyway.  Reports what actually happens to every endpoint.
# =====================================================================
sc14_nonresponder <- function() {
  arms <- list(
    "responder (DPR 0.32), untreated"  = list(DPR0 = 0.32),
    "responder (DPR 0.32), odevixibat" = c(list(DPR0 = 0.32), ODEV,
                                           list(TSTART = 0, TSTOP = yr(3))),
    "non-responder (DPR 0.15), untreated"  = list(DPR0 = 0.15),
    "non-responder (DPR 0.15), odevixibat" = c(list(DPR0 = 0.15), ODEV,
                                               list(TSTART = 0, TSTOP = yr(3)))
  )
  do.call(rbind, lapply(names(arms), function(nm) {
    d <- do.call(sim_eq, c(list(end = yr(3), delta = 14), arms[[nm]]))
    data.frame(
      arm       = nm,
      sBA       = at_day(d, "SBA_OUT", yr(3)),
      itch      = at_day(d, "ITCH", yr(3)),
      bilirubin = at_day(d, "BILI", yr(3)),
      FIB       = at_day(d, "FIB", yr(3)),
      CUMBA     = at_day(d, "CUMBA", yr(3)),
      EFS       = at_day(d, "EFS", yr(3)),
      stringsAsFactors = FALSE
    )
  }))
}

# =====================================================================
#  SCENARIO 15 — POSTNATAL DUCTULAR REPAIR (FAILURE 4)
#    The well-described spontaneous improvement of ALGS cholestasis
#    after age ~5-10 is a DPR trajectory, not a bile-acid trajectory.
#    It also means an untreated control arm improves on its own, which
#    is one reason single-arm ALGS data over-read drug effect.
# =====================================================================
sc15_ductal_repair <- function() {
  full <- sim(end = yr(17), delta = 91.3, AGE0 = 1)
  duct <- sim(end = yr(17), delta = 91.3, AGE0 = 1, KBILF = 0)  # duct term only
  ages <- c(1, 3, 5, 8, 12, 18)
  data.frame(
    age_yr        = ages,
    DPR           = sapply(ages, function(a) at_day(full, "DPR", yr(a - 1))),
    JCAP          = sapply(ages, function(a) at_day(full, "JCAP_OUT", yr(a - 1))),
    sBA           = sapply(ages, function(a) at_day(full, "SBA_OUT", yr(a - 1))),
    FIB           = sapply(ages, function(a) at_day(full, "FIB", yr(a - 1))),
    BILI_full     = sapply(ages, function(a) at_day(full, "BILI", yr(a - 1))),
    BILI_ductonly = sapply(ages, function(a) at_day(duct, "BILI", yr(a - 1))),
    itch          = sapply(ages, function(a) at_day(full, "ITCH", yr(a - 1)))
  )
}

# =====================================================================
#  RUN EVERYTHING
# =====================================================================
run_all_algs <- function() {
  cat("\n=== 0. Equilibration check ===\n");        print(sc00_equilibrate())
  cat("\n=== 1. Natural history vs GALA ===\n");    print(sc01_gala_table())
  cat("\n=== 2. ICONIC replication ===\n");         print(sc02_iconic())
  cat("\n=== 2b. Bolus-dosing equivalence ===\n");  print(sc02b_bolus_equivalence())
  cat("\n=== 3. ASSERT replication (model) ===\n"); print(sc03_assert())
  cat("\n    ASSERT observed:\n");                  print(ASSERT_OBS)
  cat("\n=== 4. Itch-per-bile-acid slopes ===\n");  print(sc04_slopes())
  cat("\n=== 5. Duct-capacity threshold ===\n")
  sw <- sc05_duct_sweep(); print(sc05_threshold(sw))
  cat("\n=== 6. PHI_EHC identifiability ===\n");    print(sc06_phi_identifiability())
  cat("\n=== 6b. Bilirubin anchor sensitivity ===\n"); print(sc06b_anchor_sensitivity())
  cat("\n=== 7. PEBD vs IBAT inhibitor ===\n");     print(sc07_pebd())
  cat("\n=== 8. Drug panel at 24 weeks ===\n");     print(sc08_drug_panel())
  cat("\n=== 9. Micellar window / site of action ===\n"); print(sc09_micellar())
  cat("\n=== 10. Treatment timing ===\n");          print(sc10_timing())
  cat("\n=== 11. Event-free survival ===\n");       print(sc11_efs())
  cat("\n=== 12. Exponent contradiction ===\n");    print(sc12_exponent_conflict())
  cat("\n=== 13. Competing-risk ceiling ===\n");    print(sc13_competing_risk())
  cat("\n=== 14. Non-responder phenotype ===\n");   print(sc14_nonresponder())
  cat("\n=== 15. Postnatal ductular repair ===\n"); print(sc15_ductal_repair())
  invisible(NULL)
}

if (identical(environment(), globalenv()) && !interactive()) {
  if (isTRUE(as.logical(Sys.getenv("ALGS_RUN_ALL", "FALSE")))) run_all_algs()
}
