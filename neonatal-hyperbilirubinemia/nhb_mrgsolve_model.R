# =============================================================================
#  nhb_mrgsolve_model.R
#  Neonatal hyperbilirubinaemia / kernicterus — Quantitative Systems Pharmacology
#  ---------------------------------------------------------------------------
#  신생아 고빌리루빈혈증(신생아 황달) QSP 모델 · mrgsolve · 34 ODEs
#
#  ===========================================================================
#  CENTRAL MODELLING THESIS
#  ===========================================================================
#  This model is NOT organised around "bilirubin goes up, phototherapy brings
#  it down".  It is organised around a single structural fact:
#
#      THE QUANTITY THAT IS MEASURED IS NOT THE QUANTITY THAT INJURES,
#      AND THE MAP BETWEEN THEM IS A SATURABLE BINDING ISOTHERM.
#
#  Total serum bilirubin (TSB) is what the laboratory reports and what every
#  guideline threshold is written against.  Unbound ("free") bilirubin Bf is
#  the only species that crosses the blood-brain barrier.  They are related by
#
#      TSB = Bf + A*K1*Bf/(1+K1*Bf) + A*n2*K2*Bf/(1+K2*Bf)
#
#  with A = albumin (uM) and K1 >> K2.  Four things follow, and all four are
#  OUTPUTS of this model rather than rules written into it.
#
#  ---------------------------------------------------------------------------
#  THESIS 1.  THE AAP "NEUROTOXICITY RISK FACTOR" THRESHOLD REDUCTION IS
#             COMPUTED, AND THE BILIRUBIN/ALBUMIN RATIO IS DERIVED, NOT
#             POSTULATED.
#
#  Inverting the isotherm at a fixed Bf = 30 nM gives the TSB that carries
#  equal risk in different infants (reference output A1b):
#
#      albumin 3.5 g/dL, pH 7.40 ............ TSB 18.58   (reference)
#      albumin 3.0 g/dL (an AAP risk factor)  TSB 15.93   -2.65 mg/dL
#      albumin 2.5 g/dL ..................... TSB 13.27   -5.31 mg/dL
#      acidaemia pH 7.15 (K1 x 0.70) ........ TSB 15.59   -2.99 mg/dL
#      ceftriaxone-type displacement ........ TSB 15.59   -2.99 mg/dL
#      sepsis + acidosis + albumin 2.6 ...... TSB 10.59   -7.99 mg/dL
#
#  The AAP 2022 curves lower the phototherapy threshold by about 2 mg/dL when
#  a neurotoxicity risk factor is present.  The model reproduces -2.65 mg/dL
#  for albumin 3.0 g/dL from binding stoichiometry alone, having never been
#  shown a guideline.
#
#  The same inversion explains WHY the bilirubin/albumin molar ratio is used
#  as a risk index and WHERE it fails.  At fixed Bf,
#
#      B/A = K1*Bf/(1+K1*Bf) + n2*K2*Bf/(1+K2*Bf)
#
#  which contains no albumin term at all: the first three rows above have
#  DIFFERENT albumin and DIFFERENT TSB but an IDENTICAL B/A of 0.604.  That is
#  the mathematical content of "B/A is better than TSB".  But B/A does contain
#  K1 — so in the acidaemia and displacement rows, at the same Bf = 30 nM, B/A
#  falls to 0.506.  B/A therefore fails precisely in the infants in whom
#  binding AFFINITY, rather than binding CAPACITY, is the problem: the septic,
#  acidotic, ceftriaxone-exposed infant.  No amount of albumin measurement
#  recovers that.
#
#  ---------------------------------------------------------------------------
#  THESIS 2.  TSB IS THE INTEGRAL OF A DIFFERENCE OF FLUXES, SO IT IS NOT
#             IDENTIFIABLE ON ITS OWN — AND END-TIDAL CO MAKES IT SO.
#
#  d(burden)/dt = production - clearance.  A doubled production and a halved
#  clearance produce indistinguishable TSB curves but respond to completely
#  different drugs.  Because haem oxygenase releases exactly one CO per
#  bilirubin, carboxyhaemoglobin is a direct read-out of the PRODUCTION arm.
#  In the reference output A3 the physiologic infant runs ETCOc 1.38 ppm and
#  the isoimmune infant 2.48-3.37 ppm at identical peak TSB values.  The model
#  therefore says which lever can work: stannsoporfin and IVIG act only on the
#  production arm, phenobarbital and gene transfer only on the clearance arm,
#  and phototherapy on neither (see Thesis 4).
#
#  ---------------------------------------------------------------------------
#  THESIS 3.  THE ASSAY COUNTS PHOTOISOMERS AS BILIRUBIN, SO PART OF THE
#             PHOTOTHERAPY RESPONSE AND PART OF THE REBOUND ARE ARTEFACTS OF
#             THE MEASUREMENT.
#
#  Phototherapy is modelled as the two reactions it actually is: a fast,
#  REVERSIBLE (4Z,15Z) <-> (4Z,15E) configurational isomerisation that reaches
#  a photostationary state, and an IRREVERSIBLE cyclisation of the E-isomer to
#  lumirubin, which is the true excretory route.  Both photoisomers are
#  measured as "bilirubin" by routine assays.  Consequence (reference A5):
#  during intensive phototherapy the photoisomers are 26 % of the reported
#  TSB, the reported number rises by ~11 % in the first hour while the native
#  pigment is already falling, and when the lamps go off the reported TSB
#  first DIPS (lumirubin, t1/2 1 h, washes out) and only then rebounds as the
#  stored E-isomer reverts thermally to Z.  Part of "rebound" is photochemical
#  bookkeeping, not new pigment.
#
#  ---------------------------------------------------------------------------
#  THESIS 4.  PHOTOTHERAPY IS THE ONLY BILIRUBIN-LOWERING ROUTE THAT DOES NOT
#             PASS THROUGH UGT1A1 — AND THE USUAL EXPLANATION FOR WHY IT STOPS
#             WORKING AS THE CHILD GROWS IS QUANTITATIVELY WRONG.
#
#  Lumirubin needs no conjugation, which is why phototherapy works at all in
#  Crigler-Najjar type I, where the model's phenobarbital arm is exactly
#  superimposable on no treatment (induction enters as GENO*ont*(1+E*C/(EC50+C))
#  and GENO = 0 annihilates the induction along with the enzyme — you cannot
#  induce a gene product that does not exist).
#
#  Every photoreaction rate in this model is (photon delivery) x (CONCENTRATION
#  in the irradiated layer), never x (amount), because a photon reaction happens
#  in the volume the light reaches.  Photon delivery therefore scales as body
#  surface area, and photo-clearance per kilogram scales as BSA/W ~ W^-0.30.
#  The textbook explanation for phototherapy failure in Crigler-Najjar is
#  exactly that term.  Running the ceiling at frozen body size from birth to
#  18 years in three nested models (reference A11b) shows it is not enough:
#
#      [G]     geometry only ...................... ceiling rises 1.89-fold
#      [G+P]   + production falls 7.7 -> 3.4 mg/kg/d  ceiling FLAT (0.84-fold)
#      [G+P+O] + age-dependent skin optics ....... ceiling rises 1.72-fold
#
#  BSA/W falls by a factor 0.42 from birth to 18 years; bilirubin production
#  per kilogram falls by 3.80/8.50 = 0.45 over the same interval.  The two
#  cancel.  Geometry plus physiology predicts that the lamps should keep
#  working, which contradicts the clinical course.  What is left to carry the
#  failure is skin optics (thicker, more pigmented skin puts less of the
#  pigment pool inside the blue-light penetration depth) and achievable
#  exposure hours.  That reframing is actionable in a way the geometric story
#  is not: defend irradiance and hours per day, and replace the missing enzyme
#  before the optical term wins.
#
#  ===========================================================================
#  PROVENANCE OF EVERY NUMBER, AND WHAT IS NOT TRUSTWORTHY
#  ===========================================================================
#  There is no R runtime in the environment in which this file was written, so
#  the equations and parameters below were developed and calibrated in an
#  INDEPENDENT pure-python RK4 implementation of the identical system,
#  nhb_reference_check.py, whose verbatim output is nhb_reference_output.txt.
#  Every quantitative claim in this header and in README.md is from that file.
#  If this R model and that output ever disagree, the python file is the one
#  that was actually executed.
#
#  Calibration targets actually hit (see nhb_reference_output.txt):
#    * bilirubin production 7.68 mg/kg/day, ~2x the adult 3.8 ................ A14a
#    * bilirubin distribution space 2.5 dL/kg, 80 % extravascular ............ A14b
#    * physiologic peak TSB 8.54 mg/dL at 108 h, 7.29 at day 7, 2.89 at
#      day 14 in an exclusively breast-fed term infant ....................... A2
#    * free bilirubin 8-21 nM at TSB 8-15 with albumin 3.5 g/dL .............. A1a
#    * ETCOc 1.38 ppm physiologic, 2.5-3.4 ppm isoimmune, ~10 ppm in a
#      fulminant G6PD crisis ................................................ A3, A8
#    * intensive phototherapy: 39 % fall in reported TSB in 24 h ............. A4
#    * photostationary E-isomer 20.6 % and lumirubin 5.5 % of reported TSB ... A5
#    * double-volume exchange: 40 % immediate fall in TSB, 82 % fall in Bf,
#      rebound to 56 % of pre-exchange ...................................... A7
#    * Crigler-Najjar II untreated plateau 23 mg/dL; ~10 % of adult UGT1A1
#      activity from gene transfer holds TSB near 9 mg/dL off the lamps ..... A11a,c
#
#  KNOWN BIASES — please read before using any absolute number:
#    1. A double-volume exchange removes ~40-45 % of the body burden here,
#       against the ~25 % of classical teaching.  The single well-mixed
#       extravascular compartment is the reason: it lets the whole
#       extravascular pool feed plasma at one fast rate constant.  A deep,
#       slowly-exchanging pool would reduce it.  The same compartment makes
#       the phototherapy ceiling in Crigler-Najjar optimistic (4.9-8.4 mg/dL
#       against a published 15-25).  These two biases have the SAME cause and
#       are linked: PSXKG cannot be tuned to fix one without breaking the
#       other, and it was fixed to the phototherapy calibration.
#    2. TAUUGT lumps UGT1A1 protein ontogeny with OATP1B1 and ligandin
#       maturation into one time constant calibrated on the TSB trajectory.
#       It is not a measurement of UGT1A1 protein.
#    3. The AAP 2022 thresholds are an ANALYTIC APPROXIMATION to the published
#       curves (plateau by gestational age, 36 h rise constant), not the
#       tabulated values.  Do not use them for care.
#    4. The skin-optics term (FOPTMIN, TAUOPT) is directionally established
#       but its parameter values are illustrative.  Thesis 4 rests on the
#       CANCELLATION of the two well-measured terms, not on these two.
#
#  TRANSCRIPTION VERIFICATION.  Because mrgsolve could not be executed here,
#  the C++ of the $GLOBAL/$ODE/$TABLE blocks was extracted, compiled standalone
#  with g++ -Wall (clean), and evaluated at a deliberately non-trivial state
#  (postnatal 31.5 h, phototherapy running, a double-volume exchange in
#  progress, a G6PD oxidant challenge active, and IVIG + stannsoporfin +
#  phenobarbital + ursodeoxycholic acid + an intraluminal binder all on board)
#  against the python reference at the identical state.  All 34 derivatives and
#  the derived outputs TSB, Bf, B/A, TcB and ETCOc agreed to better than 1e-8
#  relative, the single exception being brain bilirubin at 1.3e-9, which comes
#  from the free-bilirubin solver using 40 bisection steps here and 60 there.
#  That establishes that the two implementations are the same model; it does
#  NOT establish that mrgsolve will parse the R plumbing around them, which
#  remains unexecuted.
#
#  NOT FOR CLINICAL USE.  Educational and research model only.
# =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# =============================================================================
#  MODEL CODE
# =============================================================================
nhb_code <- '
$PROB
# Neonatal hyperbilirubinaemia QSP model (34 ODEs)
# See the R file header for the four theses and the calibration provenance.

$PARAM @annotated
// ---- anthropometry / geometry ---------------------------------------------
GA      :  40.0  : completed weeks of gestation (wk)
W0      :   3.40 : birth weight (kg)
KW      :   0.010: weight tracking rate toward the growth reference (1/h)
VP_KG   :   0.50 : plasma volume (dL/kg)
VEX_KG  :   2.00 : extravascular albumin/bilirubin space (dL/kg)
VLIV_KG :   0.40 : hepatocyte water (dL/kg)
VBL_KG  :   0.85 : whole-blood volume (dL/kg)
PSXKG   :   0.50 : plasma<->interstitium permeability (dL/h/kg)

// ---- erythron / haemolysis ------------------------------------------------
HB0     :  17.0  : birth haemoglobin (g/dL)
LRBC    :  80.0  : neonatal RBC lifespan, 80 term / 50 preterm (d)
RET0    :   3.0  : birth reticulocyte fraction (%)
TAURET  :  36.0  : marrow response time constant (h)
RETMAX  :  12.0  : maximal reticulocyte response (%)
KIMM    :   0.0045: immune haemolysis rate at ABMAT = 1 (1/h)
KOX     :   0.0040: oxidative (G6PD) haemolysis rate coefficient (1/h)
KABCL   :   0.00138: alloantibody clearance, t1/2 21 d (1/h)
KRESORB :   0.0060: extravasated-blood resorption, t1/2 115 h (1/h)
G6PD    :   0.0  : 0 normal, 1 G6PD-deficient (-)
ABMAT0  :   0.0  : birth alloantibody load, 0 = no isoimmunisation (-)
BLEX0   :   0.0  : extravasated haemoglobin at birth, cephalohaematoma (g)
FIXHB   :   0.0  : >0 clamps haemoglobin (decomposition analyses only) (g/dL)

// ---- bilirubin production -------------------------------------------------
BILPERHB:  34.0  : mg bilirubin per g haemoglobin catabolised (mg/g)
FEARLY  :   0.25 : early-labelled fraction of the RBC flux (-)
KHO1    :   0.030: HO-1 turnover (1/h)
KISNMP  :   8.0  : stannsoporfin Ki against HO-1 (ug/mL)
KCO     :   0.4694: COHb formation per bilirubin flux (%/h per mg/h/kg)
KCOOUT  :   0.1386: COHb elimination, t1/2 5 h (1/h)
PRODSCALE:  1.0  : explicit multiplier on bilirubin production (-)

// ---- albumin binding ------------------------------------------------------
ALB0    :   3.40 : birth albumin (g/dL)
ALBSET  :   3.60 : albumin set point (g/dL)
KALB    :   0.004: albumin turnover toward the set point (1/h)
ALBDON  :   4.00 : donor-blood albumin during exchange (g/dL)
KA1     :  45.0  : primary-site association constant, term neonate (1/uM)
KA2     :   1.00 : secondary-site association constant (1/uM)
NSITE2  :   1.00 : secondary sites per albumin molecule (-)
FMATK   :   1.00 : binding-affinity maturity factor, preterm < 1 (-)
FACID   :   1.00 : acidaemia factor on KA, ~0.70 at pH 7.15 (-)
FDISP   :   1.00 : displacer factor, ceftriaxone/ibuprofen/FFA ~0.6-0.8 (-)

// ---- hepatic handling -----------------------------------------------------
CLUP    :   0.800: sinusoidal uptake clearance at 3.4 kg (dL/h)
FOATP   :   1.00 : SLCO1B1 genotype factor on uptake (-)
KBACK   :   0.500: hepatocyte -> plasma efflux (1/h)
VMAXU   :  14.0  : UGT1A1 Vmax at adult activity, at 3.4 kg (mg/h)
KMU     :   5.00 : UGT1A1 Km, hepatocyte water (mg/dL)
ALLO    :   0.75 : allometric exponent for clearances (-)
KALT    :   0.015: UGT1A1-INDEPENDENT elimination of unconjugated bilirubin (1/h)
F0UGT   :   0.008: UGT1A1 activity at birth as a fraction of adult (-)
TAUUGT  : 330.0  : postnatal maturation time constant (h)
KUGT    :   0.030: UGT1A1 protein turnover, t1/2 23 h (1/h)
GENO    :   1.00 : UGT1A1 genotype activity factor (-)
EPB     :   2.50 : maximal phenobarbital induction, fold-1 (-)
EC50PB  :  20.0  : phenobarbital EC50 (mg/L)
KMRP2   :   0.400: MRP2 canalicular export (1/h)
FCHOL   :   1.00 : cholestasis factor on MRP2, bronze baby < 1 (-)

// ---- enterohepatic shunt --------------------------------------------------
KBGLUC  :   0.250: beta-glucuronidase deconjugation at BGA = 1 (1/h)
BGA0    :   1.00 : beta-glucuronidase activity at birth (-)
BGABM   :   1.00 : extra beta-glucuronidase from human milk (-)
BGAMIN  :   0.15 : floor once the flora establish (-)
TAUBGA  : 504.0  : flora-maturation time constant (h)
KREAB   :   0.200: jejunal reabsorption of unconjugated bilirubin (1/h)
KTRANS  :   0.100: gut transit / stool output at full enteral intake (1/h)
FTRANS0 :   0.35 : transit floor at zero enteral intake (-)
KB50    : 300.0  : intraluminal binder amount giving 50 % block (mg)
BGU0    :   5.0  : meconium unconjugated bilirubin load (mg)
BGC0    :  10.0  : meconium conjugated bilirubin load (mg)
BREAST  :   1.0  : 1 human milk, 0 formula (-)

// ---- phototherapy: photon delivery x CONCENTRATION, never x amount -------
KZE     :   5.44e-5: Z to E photoisomerisation (mg/h per uW/nm per mg/dL)
KEZP    :   1.02e-4: E to Z photoreversion (mg/h per uW/nm per mg/dL)
KLUMF   :   4.675e-5: E to lumirubin photocyclisation (mg/h per uW/nm per mg/dL)
I50     :  45.0  : irradiance at which the skin optically saturates (uW/cm2/nm)
KLUM    :   0.700: lumirubin elimination, t1/2 1 h, no UGT1A1 needed (1/h)
KEZBIL  :   0.020: E-isomer biliary excretion, no UGT1A1 needed (1/h)
KREV    :   0.010: E -> Z thermal (dark) reversion (1/h)
FOPTMIN :   0.30 : asymptotic optical accessibility of adult skin (-)
TAUOPT  :   3.00 : optical-accessibility decay constant (years)
AGEOFF  :   0.0  : hours added to t for the skin-optics term only (h)

// ---- exchange transfusion -------------------------------------------------
HBDON   :  13.0  : reconstituted donor-blood haemoglobin (g/dL)
ETEFF   :   0.45 : exchange efficiency, mixing/recirculation losses (-)
FABREM  :   1.00 : efficiency of antibody / coated-RBC removal (-)

// ---- drug PK --------------------------------------------------------------
KASNMP  :   0.173: stannsoporfin IM absorption, t1/2 4 h (1/h)
VSNMP   :   3.00 : stannsoporfin volume of distribution (dL/kg)
KESNMP  :   0.023: stannsoporfin elimination, t1/2 30 h (1/h)
KAPB    :   0.400: phenobarbital absorption (1/h)
VPB     :   9.00 : phenobarbital volume of distribution (dL/kg)
KEPB    :   0.00693: phenobarbital elimination, t1/2 100 h (1/h)
VIGG    :   0.50 : IVIG central volume (dL/kg)
KIGGCP  :   0.010: IVIG central -> peripheral (1/h)
KIGGPC  :   0.008: IVIG peripheral -> central (1/h)
KIGGEL  :   0.00120: IVIG elimination, t1/2 24 d (1/h)
IGG0    :  10.0  : endogenous transplacental IgG (g/L)
IMAXIVIG:   0.65 : maximal fractional block of immune haemolysis (-)
IC50IVIG:   6.00 : IVIG EC50 above baseline (g/L)
KAUDCA  :   0.500: ursodeoxycholic acid absorption (1/h)
VUDCA   :   3.00 : ursodeoxycholic acid volume of distribution (dL/kg)
KEUDCA  :   0.140: ursodeoxycholic acid elimination (1/h)
EMAXU   :   0.80 : maximal ursodeoxycholic acid effect (-)
EC50UDCA:   4.00 : ursodeoxycholic acid EC50 (umol/L)
KTGON   :   0.00289: transgene expression onset, t1/2 10 d (1/h)
TGXMAX  :   0.100: plateau transgene activity as a fraction of adult (-)
GTSTART : 1.0e9  : time of AAV8-hUGT1A1 administration (h)

// ---- blood-brain barrier and injury --------------------------------------
KINBBB  :   0.100: free bilirubin -> brain influx (1/h)
KOUTBBB :   0.100: brain efflux including P-glycoprotein (1/h)
FBBB    :   1.00 : BBB permeability multiplier, prematurity/sepsis (-)
BBRTHR  :  35.0  : injury threshold for brain bilirubin (nM-eq)
KINJ    :   5.0e-4: injury accrual (1/(nM*h))
KREP    :   0.0020: repair of the reversible injury component (1/h)
BBRTHRA :  25.0  : ABR threshold, lower than the injury threshold (nM-eq)
KABR    :   4.0e-4: ABR deficit accrual (1/(nM*h))
KABRREC :   0.020: ABR recovery (1/h)
INJ50   :   0.50 : kernicterus logistic midpoint (-)
INJSL   :   0.08 : kernicterus logistic slope (-)

// ---- AAP 2022 thresholds (analytic approximation) ------------------------
RF      :   0.0  : 1 if any neurotoxicity risk factor is present (-)

// ---- scenario controls ---------------------------------------------------
IRRSET  :   0.0  : phototherapy spectral irradiance when on (uW/cm2/nm)
FBSASET :   0.80 : irradiated fraction of body surface area (-)
PTSTART : 1.0e9  : phototherapy start time (h)
PTSTOP  : 1.0e9  : phototherapy stop time (h)
PTDUTY  :  24.0  : hours of phototherapy per 24 h in open-loop mode (h)
AUTOPT  :   0.0  : 1 = closed-loop phototherapy on the AAP threshold (-)
PTOFFM  :   2.0  : closed-loop switch-off margin below the threshold (mg/dL)
ETRATE  :   0.0  : exchange rate while running (mL/kg/h)
ETSTART : 1.0e9  : exchange transfusion start time (h)
ETDUR   :   3.0  : exchange transfusion duration (h)
OXSET   :   0.0  : oxidant challenge multiplier while running (-)
OXSTART : 1.0e9  : oxidant challenge start time (h)
OXDUR   :  24.0  : oxidant challenge duration (h)
FIEARLY :   1.00 : enteral intake adequacy before FISW (-)
FILATE  :   1.00 : enteral intake adequacy after FISW (-)
FISW    :   0.0  : time at which enteral intake changes (h)

$CMT @annotated
W     : body weight (kg)
HB    : haemoglobin concentration (g/dL)
RET   : reticulocyte fraction, marrow output signal (%)
ABMAT : maternal alloantibody load on neonatal RBC (-)
HO1   : haem oxygenase-1 activity (-)
COHB  : carboxyhaemoglobin (%)
BLEX  : extravasated blood haemoglobin (g)
BP    : plasma unconjugated bilirubin (mg)
BEX   : extravascular unconjugated bilirubin, the photon target (mg)
BLIV  : hepatocyte unconjugated bilirubin, ligandin pool (mg)
BCON  : hepatocyte conjugated bilirubin (mg)
BGC   : gut-lumen conjugated bilirubin (mg)
BGU   : gut-lumen unconjugated bilirubin, reabsorbable (mg)
BSTL  : cumulative faecal + renal bilirubin output (mg)
LUMI  : plasma lumirubin (mg)
EZ    : plasma (4Z,15E) configurational photoisomer (mg)
ALB   : serum albumin (g/dL)
UGT   : UGT1A1 activity as a fraction of adult (-)
TGX   : AAV8-hUGT1A1 transgene-derived activity (-)
BBR   : basal-ganglia bilirubin (nM-eq)
INJ   : cumulative neuronal injury index (-)
ABRD  : auditory brainstem response deficit (-)
ASNMP : stannsoporfin IM depot (mg)
CSNMP : stannsoporfin plasma (ug/mL)
APB   : phenobarbital gut depot (mg)
CPB   : phenobarbital plasma (mg/L)
IGGC  : IVIG central (g)
IGGP  : IVIG peripheral (g)
AUDCA : ursodeoxycholic acid gut depot (mg)
CUDCA : ursodeoxycholic acid plasma (umol/L)
GBIND : intraluminal bilirubin binder, agar/charcoal/zinc (mg)
AUCX  : AUC of TSB above the AAP phototherapy threshold (mg*h/dL)
PTH   : cumulative phototherapy exposure (h)
ETV   : cumulative exchange-transfused volume (mL/kg)

$GLOBAL
namespace {

// --- WHO median growth reference, scaled to birth weight -------------------
// The first fortnight is written out explicitly so that the physiologic
// postnatal weight nadir (~6 % at day 3) is part of the reference and not a
// separate term.
const int NW = 15;
const double WT_AGE[NW] = {0.0, 3.0/365.0, 7.0/365.0, 14.0/365.0, 28.0/365.0,
                           0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 14.0, 18.0};
const double WT_REL[NW] = {1.000, 0.939, 0.985, 1.059, 1.191, 1.765, 2.294,
                           2.824, 3.588, 4.794, 6.029, 7.441, 9.382, 14.706,
                           19.412};
const int NH = 11;
const double HT_AGE[NH] = {0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0,
                           14.0, 18.0};
const double HT_CM[NH]  = {50.0, 61.4, 67.6, 75.7, 87.1, 103.3, 116.0, 127.3,
                           137.8, 163.2, 176.1};

double interp(const double *x, const double *y, int n, double xx) {
  if (xx <= x[0]) return y[0];
  if (xx >= x[n-1]) return y[n-1];
  for (int i = 0; i < n-1; ++i)
    if (xx >= x[i] && xx <= x[i+1])
      return y[i] + (y[i+1]-y[i]) * (xx-x[i]) / (x[i+1]-x[i]);
  return y[n-1];
}

double wref(double t_h, double w0) {
  return w0 * interp(WT_AGE, WT_REL, NW, t_h/8760.0);
}

double href(double t_h, double w0) {
  double sc = pow(w0/3.40, 1.0/3.0);
  double h  = interp(HT_AGE, HT_CM, NH, t_h/8760.0);
  return (t_h < 8760.0) ? h*sc : h;
}

// Haycock body surface area in cm2: BSA(m2) = 0.024265 * W^0.5378 * H^0.3964
double bsa_cm2(double w, double h) {
  return 0.024265 * pow(w, 0.5378) * pow(h, 0.3964) * 1.0e4;
}

// --- free (unbound) bilirubin from a TWO-CLASS saturable isotherm ----------
// BT = Bf + A*K1*Bf/(1+K1*Bf) + A*n2*K2*Bf/(1+K2*Bf)
// Monotone in Bf, so 40 bisection steps; there is no closed form and a
// one-site isotherm predicts physically impossible Bf above B/A ~ 1.
double bfree_uM(double cp_mgdl, double alb_gdl, double ka1, double ka2,
                double n2, double fk) {
  double bt = cp_mgdl * 17.1;        // mg/dL -> umol/L
  double a  = alb_gdl * 150.4;       // g/dL  -> umol/L (albumin MW 66.5 kDa)
  if (bt <= 0.0) return 0.0;
  double k1 = ka1 * fk, k2 = ka2 * fk;
  double lo = 0.0, hi = bt, bf, tot;
  for (int i = 0; i < 40; ++i) {
    bf  = 0.5*(lo+hi);
    tot = bf + a*k1*bf/(1.0+k1*bf) + a*n2*k2*bf/(1.0+k2*bf);
    if (tot < bt) lo = bf; else hi = bf;
  }
  return 0.5*(lo+hi);
}

// --- AAP 2022 thresholds, ANALYTIC APPROXIMATION (not the tabulated curves)
double thr_pt(double t_h, double ga, double rf) {
  double gaeff   = (ga > 40.0) ? 40.0 : ga;
  double plateau = 16.5 + 1.0*(gaeff - 35.0);
  if (plateau > 21.0) plateau = 21.0;
  plateau -= 2.0*rf;
  return plateau * (0.42 + 0.58*(1.0 - exp(-t_h/36.0)));
}
double thr_et(double t_h, double ga, double rf) {
  return thr_pt(t_h, ga, rf) + (3.5 - 1.0*rf);
}

// --- haemoglobin set point sensed by the renal/hepatic oxygen sensor -------
// Falls after birth (arterial pO2 jumps, HbF -> HbA) and recovers.  This is
// the mechanism of the physiologic nadir; there is no separate "nadir" term.
double hbset(double t_h) {
  double d = t_h/24.0;
  if (d <  3.0)   return 17.0;
  if (d < 60.0)   return 17.0 - 6.0*(d-3.0)/57.0;
  if (d < 120.0)  return 11.0 + 1.5*(d-60.0)/60.0;
  return 12.5;
}

double sig(double x) {
  if (x >  40.0) return 1.0;
  if (x < -40.0) return 0.0;
  return 1.0/(1.0+exp(-x));
}

double win(double t, double t0, double dur, double val) {
  return (t >= t0 && t < t0 + dur) ? val : 0.0;
}

} // namespace

double CPNOW, CISONOW, BFNOW, BANOW, TSBNOW, THRPTNOW, THRETNOW, PRODNOW;
double ETCONOW, TCBNOW, IRRNOW, KERNNOW, PHOTONOW, UPNOW, CONJNOW, REABNOW;
double ALBFNOW, WNOW, BSANOW, FIVIGNOW, FSNMPNOW, EZPCTNOW;

$MAIN
W_0     = W0;
HB_0    = (FIXHB > 0.0) ? FIXHB : HB0;
RET_0   = RET0;
ABMAT_0 = ABMAT0;
HO1_0   = 1.0;
COHB_0  = 1.20;
BLEX_0  = BLEX0;
// cord TSB 1.5 mg/dL, already equilibrated across plasma and interstitium
BP_0    = 1.5 * VP_KG  * W0;
BEX_0   = 1.5 * VEX_KG * W0;
ALB_0   = ALB0;
UGT_0   = GENO * F0UGT;
BGU_0   = BGU0;
BGC_0   = BGC0;
IGGC_0  = IGG0 * VIGG * W0 / 10.0;

$ODE
double Wc   = (W    > 0.4) ? W   : 0.4;
double ALBc = (ALB  > 0.5) ? ALB : 0.5;
double t    = SOLVERTIME;

// ---- geometry -------------------------------------------------------------
double HT   = href(t, W0);
double BSA  = bsa_cm2(Wc, HT);
double VP   = VP_KG   * Wc;
double VEX  = VEX_KG  * Wc;
double VLIV = VLIV_KG * Wc;
double VBL  = VBL_KG  * Wc;
double VEZ  = VP + VEX;
double allo = pow(Wc/3.40, ALLO);

// ---- concentrations and the binding isotherm -----------------------------
double CP   = BP   / VP;                      // native (4Z,15Z) UCB, mg/dL
double CEX  = BEX  / VEX;
double CLIV = BLIV / VLIV;
double CEZ  = EZ   / VEZ;
double CISO = (LUMI + EZ) / VEZ;              // photoisomers, mg/dL
double TSB  = CP + CISO;                      // what the laboratory reports
// The photoisomers occupy albumin as well, so the site pool left for the
// native pigment is smaller than total albumin.
double ALBF = ALBc - CISO*17.1/150.4;
if (ALBF < 0.30) ALBF = 0.30;
double fk   = FMATK * FACID * FDISP;
double BFuM = bfree_uM(CP, ALBF, KA1, KA2, NSITE2, fk);
double BFnM = 1000.0 * BFuM;
double BA   = (CP*17.1) / (ALBc*150.4);

double THRPT = thr_pt(t, GA, RF);
double THRET = thr_et(t, GA, RF);

// ---- drug effects --------------------------------------------------------
double CIGG   = 10.0 * IGGC / (VIGG*Wc);
double iggex  = CIGG - IGG0; if (iggex < 0.0) iggex = 0.0;
double F_IVIG = IMAXIVIG*iggex/(IC50IVIG+iggex);
double F_SNMP = 1.0/(1.0 + CSNMP/KISNMP);     // competitive HO-1 inhibition
double F_UDCA = EMAXU*CUDCA/(EC50UDCA+CUDCA);
double F_PB   = 1.0 + EPB*CPB/(EC50PB+CPB);

// ---- haemolysis ----------------------------------------------------------
double k_basal = 1.0/(LRBC*24.0);
double k_imm   = KIMM*ABMAT*(1.0-F_IVIG);
double ox      = win(t, OXSTART, OXDUR, OXSET);
double k_ox    = KOX*ox*G6PD;
double k_hem   = k_basal + k_imm + k_ox;

// ---- exchange transfusion: a real 2-4 h procedure, not an instant jump ---
double etr  = win(t, ETSTART, ETDUR, ETRATE);
double QET  = etr*Wc/100.0*ETEFF;             // dL/h of plasma turned over
double eton = (etr > 0.0) ? 1.0 : 0.0;

// ---- enteral intake -----------------------------------------------------
double fi = (t < FISW) ? FIEARLY : FILATE;

// ---- erythron -----------------------------------------------------------
double hset = hbset(t);
double hgap = (hset - HB)/((hset > 1.0) ? hset : 1.0); if (hgap < 0.0) hgap = 0.0;
double rtgt = RET0*exp(-t/1000.0) + RETMAX*hgap/(0.15+hgap);
double kprod = k_basal*HB0/RET0;
dxdt_RET = (FIXHB > 0.0) ? 0.0 : (rtgt - RET)/TAURET;
dxdt_HB  = (FIXHB > 0.0) ? 0.0 :
           (-k_hem*HB + kprod*RET + QET/VBL*(HBDON-HB)*eton
            - HB*(KW*(wref(t,W0)-Wc))/Wc);
dxdt_ABMAT = -KABCL*ABMAT - QET/VP*ABMAT*FABREM*eton;
dxdt_BLEX  = -KRESORB*BLEX;

// ---- bilirubin production ----------------------------------------------
double hbflux    = k_hem*HB*VBL;               // g Hb/h from circulating RBC
double prod_rbc  = BILPERHB*hbflux;
double prod_earl = FEARLY*prod_rbc;
double prod_ex   = BILPERHB*KRESORB*BLEX;
double PROD      = (prod_rbc + prod_earl + prod_ex)*HO1*F_SNMP*PRODSCALE;

double himm = k_imm/0.01; if (himm > 1.0) himm = 1.0;
dxdt_HO1  = KHO1*(1.0 + 0.30*himm - HO1);
dxdt_COHB = KCO*PROD/Wc - KCOOUT*COHB;

// ---- phototherapy: photon delivery and the two photoreactions -----------
double irr;
if (AUTOPT > 0.5) {
  // closed loop on the infant own AAP threshold, smoothed for the solver
  irr = IRRSET * sig((TSB - (THRPT - 0.5*PTOFFM))/0.35);
} else {
  double ph = fmod(t - PTSTART, 24.0);
  irr = (t >= PTSTART && t < PTSTOP && ph < PTDUTY) ? IRRSET : 0.0;
}
double ieff = irr/(1.0 + irr/I50);             // optical saturation of skin
double fbsa = (irr > 0.0) ? FBSASET : 0.0;
double fopt = FOPTMIN + (1.0-FOPTMIN)*exp(-((t+AGEOFF)/8760.0)/TAUOPT);
double u    = ieff*fbsa*BSA*fopt;              // photon delivery
// Rates are photon delivery x CONCENTRATION in the irradiated layer.  Writing
// them against the AMOUNT makes delivery scale as BSA*W and silently inverts
// every size and growth conclusion.
double PZE  = KZE  *u*CEX;                     // Z -> E in irradiated skin
double PEZ  = KEZP *u*CEZ;                     // E -> Z photoreversion
double PLUM = KLUMF*u*CEZ;                     // E -> lumirubin, one-way
double TREV = KREV *EZ;                        // E -> Z thermal, in the dark
double EBIL = KEZBIL*EZ;                       // E-isomer into bile
dxdt_PTH = (irr > 0.0) ? 1.0 : 0.0;

// ---- hepatic handling --------------------------------------------------
double UP   = CLUP*allo*FOATP*CP;
double BACK = KBACK*BLIV;
double CONJ = VMAXU*allo*UGT*CLIV/(KMU+CLIV);
double EXPB = KMRP2*FCHOL*BCON;
double ALTE = KALT*BLIV;                       // non-UGT1A1 escape route

// ---- enterohepatic shunt ----------------------------------------------
double bga = BGAMIN + (BGA0*(1.0+BGABM*BREAST) - BGAMIN)*exp(-t/TAUBGA);
double ftr = (FTRANS0 + (1.0-FTRANS0)*fi)*(1.0+0.5*F_UDCA);
if (t < 12.0) ftr *= 0.35;                     // meconium has not passed
double occ  = GBIND/(KB50+GBIND);
double DECON = KBGLUC*bga*BGC;
double REAB  = KREAB*BGU*(1.0-occ)*(1.0-0.25*F_UDCA);
double TRC   = KTRANS*ftr*BGC;
double TRU   = KTRANS*ftr*BGU;

// ---- bilirubin mass balance ------------------------------------------
double PSX = PSXKG*Wc;
double fx  = PSX*(CP - CEX);
// The photoisomers are distributed over the WHOLE albumin space, so the
// reverse reactions must return pigment to plasma and interstitium in
// proportion to their volumes.  Returning all of it to plasma would turn the
// futile Z->E->Z cycle into a photon-driven pump out of the tissues.
double RBP = (PEZ+TREV)*VP /VEZ;
double RBX = (PEZ+TREV)*VEX/VEZ;

dxdt_BP   = PROD - fx - UP + BACK + RBP - QET*CP*eton;
dxdt_BEX  = fx - PZE + RBX;
dxdt_BLIV = UP - BACK - CONJ - ALTE + REAB;
dxdt_BCON = CONJ - EXPB;
dxdt_BGC  = EXPB - DECON - TRC;
dxdt_BGU  = DECON - REAB - TRU;
dxdt_LUMI = PLUM - KLUM*LUMI - QET*(LUMI/VEZ)*eton;
dxdt_EZ   = PZE - PEZ - PLUM - TREV - EBIL - QET*(EZ/VEZ)*eton;
dxdt_BSTL = TRC + TRU + KLUM*LUMI + EBIL + ALTE;
dxdt_ETV  = etr;

// ---- albumin ---------------------------------------------------------
dxdt_ALB = KALB*(ALBSET-ALB) + QET/VP*(ALBDON-ALB)*eton;

// ---- UGT1A1: ontogeny x genotype x induction, with protein turnover ---
double ont  = F0UGT + (1.0-F0UGT)*(1.0-exp(-t/TAUUGT));
double utgt = GENO*ont*F_PB + TGX;
dxdt_UGT = KUGT*(utgt - UGT);
dxdt_TGX = (t >= GTSTART) ? KTGON*(TGXMAX - TGX) : 0.0;

// ---- neurotoxicity ---------------------------------------------------
dxdt_BBR = KINBBB*FBBB*BFnM - KOUTBBB*BBR;
double exc = BBR - BBRTHR;  if (exc < 0.0) exc = 0.0;
double exa = BBR - BBRTHRA; if (exa < 0.0) exa = 0.0;
double injc = (INJ < 1.0) ? INJ : 1.0;
dxdt_INJ  = KINJ*exc*(1.0-injc) - KREP*injc;
dxdt_ABRD = KABR*exa - KABRREC*ABRD;

// ---- drug PK ---------------------------------------------------------
dxdt_ASNMP = -KASNMP*ASNMP;
dxdt_CSNMP = KASNMP*ASNMP/(VSNMP*Wc)*10.0 - KESNMP*CSNMP;
dxdt_APB   = -KAPB*APB;
dxdt_CPB   = KAPB*APB/(VPB*Wc)*10.0 - KEPB*CPB;
dxdt_IGGC  = -KIGGCP*IGGC + KIGGPC*IGGP - KIGGEL*IGGC - QET/VP*IGGC*eton;
dxdt_IGGP  = KIGGCP*IGGC - KIGGPC*IGGP;
dxdt_AUDCA = -KAUDCA*AUDCA;
dxdt_CUDCA = KAUDCA*AUDCA/(VUDCA*Wc)*10.0 - KEUDCA*CUDCA;
dxdt_GBIND = -KTRANS*ftr*GBIND;

// ---- growth and threshold bookkeeping -------------------------------
dxdt_W = KW*(wref(t,W0)*(1.0-0.12*(1.0-fi)) - Wc);
double over = TSB - THRPT; if (over < 0.0) over = 0.0;
dxdt_AUCX = over;

// ---- hand the derived quantities to the output block ---------------
CPNOW = CP; CISONOW = CISO; BFNOW = BFnM; BANOW = BA; TSBNOW = TSB;
THRPTNOW = THRPT; THRETNOW = THRET; PRODNOW = PROD; ETCONOW = 1.15*COHB;
TCBNOW = 0.92*CEX; IRRNOW = irr; WNOW = Wc; BSANOW = BSA; ALBFNOW = ALBF;
PHOTONOW = PLUM + EBIL; UPNOW = UP; CONJNOW = CONJ; REABNOW = REAB;
FIVIGNOW = F_IVIG; FSNMPNOW = F_SNMP;
EZPCTNOW = (TSB > 1e-9) ? 100.0*CISO/TSB : 0.0;

$TABLE
double TSBOUT   = TSBNOW;          // total serum bilirubin as measured (mg/dL)
double UCBNAT   = CPNOW;           // native (4Z,15Z) pigment only (mg/dL)
double ISOMER   = CISONOW;         // photoisomers measured as bilirubin
double ISOPCT   = EZPCTNOW;
double BF       = BFNOW;           // FREE bilirubin (nM) - the toxic species
double BARATIO  = BANOW;
double TCB      = TCBNOW;
double THRESHPT = THRPTNOW;
double THRESHET = THRETNOW;
double ESCAL    = THRETNOW - 2.0;
double PROD24   = PRODNOW*24.0/WNOW;   // bilirubin production (mg/kg/day)
double ETCOC    = ETCONOW;             // end-tidal CO (ppm)
double KERN     = sig((INJ-INJ50)/INJSL);
double PHOTOFLX = PHOTONOW;            // irreversible photochemical removal
double IRRAD    = IRRNOW;
double BSAOUT   = BSANOW;
double WTOUT    = WNOW;
double ALBFREE  = ALBFNOW;
double UGTPCT   = 100.0*UGT;
double OVERPT   = (TSBNOW > THRPTNOW) ? 1.0 : 0.0;
double OVERET   = (TSBNOW > THRETNOW) ? 1.0 : 0.0;

$CAPTURE
TSBOUT UCBNAT ISOMER ISOPCT BF BARATIO TCB THRESHPT THRESHET ESCAL
PROD24 ETCOC KERN PHOTOFLX IRRAD BSAOUT WTOUT ALBFREE UGTPCT OVERPT OVERET
'

mod <- mcode("nhb", nhb_code)

# =============================================================================
#  GENOTYPES — UGT1A1 activity as a fraction of wild-type adult activity
# =============================================================================
GENOTYPES <- c(
  wild          = 1.00,   # *1/*1
  UGT1A1_28_het = 0.65,   # (TA)7 heterozygote
  Gilbert       = 0.35,   # *28/*28
  UGT1A1_6      = 0.40,   # G71R homozygote, common in East Asia
  CN2           = 0.05,   # Crigler-Najjar type II (Arias)
  CN1           = 0.00    # Crigler-Najjar type I
)

# =============================================================================
#  DOSING HELPERS
#  Every dose is an amount into the relevant depot compartment.  Weights are
#  the scenario birth weight, so amounts are computed per scenario.
# =============================================================================
dose_ivig  <- function(wt, time = 20, g_per_kg = 1.0)
  ev(time = time, amt = g_per_kg*wt, cmt = "IGGC")

dose_snmp  <- function(wt, time = 24, mg_per_kg = 4.5)
  ev(time = time, amt = mg_per_kg*wt, cmt = "ASNMP")

dose_pheno <- function(wt, start = 0, days = 14, mg_per_kg = 5.0)
  ev(time = start, amt = mg_per_kg*wt, cmt = "APB", ii = 24, addl = days-1)

dose_udca  <- function(wt, start = 0, days = 14, mg_per_kg = 10.0)
  ev(time = start, amt = mg_per_kg*wt, cmt = "AUDCA", ii = 12, addl = 2*days-1)

dose_agar  <- function(wt, start = 0, days = 7, mg_per_kg_day = 250)
  ev(time = start, amt = mg_per_kg_day*wt/4, cmt = "GBIND", ii = 6,
     addl = 4*days-1)

# =============================================================================
#  THE TEN SHIPPED SCENARIOS
#  Phototherapy is CLOSED-LOOP (AUTOPT = 1) in the scenario library: it comes
#  on when TSB crosses that infant's own AAP-2022 threshold and goes off
#  PTOFFM below it, which is how it is actually given.  Phototherapy hours are
#  therefore an OUTPUT of the model rather than an assumption.
# =============================================================================
AUTO <- list(AUTOPT = 1, IRRSET = 30, FBSASET = 0.80, PTOFFM = 2.0)

scenarios <- list(

  S1 = list(
    label = "S1 physiologic term, exclusively breast-fed",
    par   = c(AUTO),
    ev    = NULL, wt = 3.40),

  S2 = list(
    label = "S2 suboptimal intake, delayed stooling (AAP 2022 term)",
    par   = c(AUTO, list(FIEARLY = 0.30, FISW = 96, FILATE = 0.95)),
    ev    = NULL, wt = 3.40),

  S3 = list(
    # isoimmune haemolytic disease is itself an AAP neurotoxicity risk factor,
    # so RF = 1 lowers this infant's own thresholds
    label = "S3 ABO isoimmune haemolytic disease, DAT+",
    par   = c(AUTO, list(ABMAT0 = 0.12, RF = 1)),
    ev    = NULL, wt = 3.40),

  S4 = list(
    label = "S4 Rh(D) disease + IVIG 1 g/kg",
    par   = c(AUTO, list(ABMAT0 = 0.30, HB0 = 13.5, ALB0 = 3.0, RF = 1)),
    ev    = dose_ivig(3.40, 20), wt = 3.40),

  S5 = list(
    label = "S5 Rh(D) disease + IVIG + double-volume exchange at 30 h",
    par   = c(AUTO, list(ABMAT0 = 0.30, HB0 = 13.5, ALB0 = 3.0, RF = 1,
                         ETSTART = 30, ETDUR = 3, ETRATE = 56.7)),
    ev    = dose_ivig(3.40, 20), wt = 3.40),

  S6 = list(
    label = "S6 G6PD deficiency, 24 h oxidant challenge from 60 h",
    par   = c(AUTO, list(G6PD = 1, RF = 1, OXSTART = 60, OXDUR = 24,
                         OXSET = 1.0)),
    ev    = NULL, wt = 3.40),

  S7 = list(
    label = "S7 late preterm 35 wk with neurotoxicity risk factors",
    par   = c(AUTO, list(GA = 35, W0 = 2.40, ALB0 = 2.8, ALBSET = 3.0,
                         LRBC = 50, FMATK = 0.85, RF = 1, FBBB = 1.6,
                         ABMAT0 = 0.12)),
    ev    = NULL, wt = 2.40),

  S8 = list(
    label = "S8 UGT1A1*6/*6 prolonged jaundice (East Asian)",
    par   = c(AUTO, list(GENO = GENOTYPES[["UGT1A1_6"]])),
    ev    = NULL, wt = 3.40),

  S9 = list(
    label = "S9 Crigler-Najjar type II + phenobarbital 5 mg/kg/day",
    par   = c(AUTO, list(GENO = GENOTYPES[["CN2"]])),
    ev    = dose_pheno(3.40, 0, 14), wt = 3.40),

  S10 = list(
    label = "S10 ABO isoimmune + stannsoporfin 4.5 mg/kg IM",
    par   = c(AUTO, list(ABMAT0 = 0.12, RF = 1)),
    ev    = dose_snmp(3.40, 24), wt = 3.40)
)

run_scenario <- function(s, end = 336, delta = 0.5) {
  m <- param(mod, s$par)
  if (is.null(s$ev)) {
    out <- mrgsim(m, end = end, delta = delta, hmax = 0.25)
  } else {
    out <- mrgsim(m, events = s$ev, end = end, delta = delta, hmax = 0.25)
  }
  as_tibble(out) %>% mutate(scenario = s$label)
}

run_all <- function(end = 336) bind_rows(lapply(scenarios, run_scenario,
                                                end = end))

# =============================================================================
#  SUMMARY TABLE — the shape of nhb_reference_output.txt section A13
# =============================================================================
summarise_scenarios <- function(sims) {
  sims %>%
    group_by(scenario) %>%
    summarise(
      peak_TSB   = max(TSBOUT),
      t_peak_h   = time[which.max(TSBOUT)],
      peak_Bf_nM = max(BF),
      PT_hours   = max(PTH),
      h_over_PT  = sum(OVERPT)*(time[2]-time[1]),
      h_over_ET  = sum(OVERET)*(time[2]-time[1]),
      AUC_over   = max(AUCX),
      peak_ETCOc = max(ETCOC),
      Hb_day7    = HB[which.min(abs(time-168))],
      injury     = last(INJ),
      P_kern     = last(KERN),
      .groups    = "drop") %>%
    arrange(scenario)
}

# =============================================================================
#  ANALYSIS 1 — the iso-Bf contour: the AAP risk-factor threshold reduction
#  falls out of binding stoichiometry.  Reproduces reference output A1b.
# =============================================================================
free_bilirubin <- function(tsb, alb, ka1 = 45, ka2 = 1, n2 = 1, fk = 1) {
  bt <- tsb*17.1; a <- rep_len(alb*150.4, length(bt))
  k1 <- ka1*fk;   k2 <- ka2*fk
  vapply(seq_along(bt), function(i) {
    lo <- 0; hi <- bt[i]
    for (j in 1:60) {
      bf  <- (lo+hi)/2
      tot <- bf + a[i]*k1*bf/(1+k1*bf) + a[i]*n2*k2*bf/(1+k2*bf)
      if (tot < bt[i]) lo <- bf else hi <- bf
    }
    1000*(lo+hi)/2
  }, numeric(1))
}

iso_bf_table <- function(target_nM = 30) {
  cond <- tibble::tribble(
    ~condition,                              ~alb, ~fk,
    "term, pH 7.40, no displacer",            3.5, 1.00,
    "albumin 3.0 g/dL (an AAP risk factor)",  3.0, 1.00,
    "albumin 2.5 g/dL",                       2.5, 1.00,
    "acidaemia pH 7.15 (KA x 0.70)",          3.5, 0.70,
    "ceftriaxone / ibuprofen displacement",   3.5, 0.70,
    "preterm 30 wk, albumin 2.8, KA x 0.75",  2.8, 0.75,
    "sepsis + acidosis + albumin 2.6",        2.6, 0.595)
  cond %>%
    rowwise() %>%
    mutate(TSB_at_target = {
      lo <- 0; hi <- 60
      for (j in 1:80) {
        mid <- (lo+hi)/2
        if (free_bilirubin(mid, alb, fk = fk) < target_nM) lo <- mid else hi <- mid
      }
      (lo+hi)/2
    }) %>%
    ungroup() %>%
    mutate(dTSB   = TSB_at_target - first(TSB_at_target),
           B_A    = TSB_at_target*17.1/(alb*150.4))
}
# The B_A column is IDENTICAL for rows 1-3 and lower for rows 4-7: the
# bilirubin/albumin ratio is an albumin-free surrogate for Bf but NOT an
# affinity-free one.

# =============================================================================
#  ANALYSIS 2 — phototherapy dose-response: irradiance versus exposed area
# =============================================================================
pt_dose_response <- function() {
  grid <- expand.grid(IRRSET = c(0, 8, 15, 30, 50),
                      FBSASET = c(0.35, 0.80, 1.00))
  bind_rows(lapply(seq_len(nrow(grid)), function(i) {
    m <- param(mod, ABMAT0 = 0.12, RF = 1, PTSTART = 24, PTSTOP = 1e9,
               IRRSET = grid$IRRSET[i], FBSASET = grid$FBSASET[i])
    as_tibble(mrgsim(m, end = 96, delta = 1, hmax = 0.25)) %>%
      mutate(irr = grid$IRRSET[i], fbsa = grid$FBSASET[i])
  })) %>%
    group_by(irr, fbsa) %>%
    summarise(TSB24 = TSBOUT[time == 24],
              TSB48 = TSBOUT[time == 48],
              TSB72 = TSBOUT[time == 72],
              pct_change_24h = 100*(TSBOUT[time == 48]/TSBOUT[time == 24] - 1),
              .groups = "drop")
}

# =============================================================================
#  ANALYSIS 3 — photoisomers are measured as bilirubin
# =============================================================================
photoisomer_run <- function() {
  m <- param(mod, ABMAT0 = 0.12, RF = 1, PTSTART = 24, PTSTOP = 72,
             IRRSET = 30, FBSASET = 0.80)
  as_tibble(mrgsim(m, end = 144, delta = 0.5, hmax = 0.1)) %>%
    select(time, TSBOUT, UCBNAT, ISOMER, ISOPCT, BF)
}

# =============================================================================
#  ANALYSIS 4 — exchange transfusion mechanics.  Nothing tells the model to
#  remove "50 %": the removal is the integral of QET*Cp over the procedure,
#  limited by refilling from the extravascular pool.
# =============================================================================
exchange_run <- function() {
  m <- param(mod, ABMAT0 = 0.30, HB0 = 13.5, ALB0 = 2.9, RF = 1,
             PTSTART = 18, PTSTOP = 1e9, IRRSET = 30, FBSASET = 0.80,
             ETSTART = 30, ETDUR = 3, ETRATE = 56.7)
  as_tibble(mrgsim(m, end = 96, delta = 0.25, hmax = 0.05)) %>%
    select(time, TSBOUT, BF, ALB, HB, ABMAT, BARATIO)
}

# =============================================================================
#  ANALYSIS 5 — Crigler-Najjar: the phototherapy ceiling versus body size,
#  decomposed into geometry [G], + production ontogeny [G+P], and + skin
#  optics [G+P+O].  Reproduces reference output A11b, whose conclusion is that
#  the surface-to-mass term is CANCELLED by production ontogeny.
# =============================================================================
PROD_ONTOGENY <- tibble::tibble(
  age_y = c(0, 0.25, 0.5, 1, 2, 4, 6, 8, 10, 14, 18),
  rel   = c(1.000, 0.850, 0.780, 0.700, 0.620, 0.550, 0.520, 0.500,
            0.480, 0.460, 0.447))
W_REL <- tibble::tibble(
  age_y = c(0, 0.25, 0.5, 1, 2, 4, 6, 8, 10, 14, 18),
  rel   = c(1.000, 1.765, 2.294, 2.824, 3.588, 4.794, 6.029, 7.441,
            9.382, 14.706, 19.412))

cn_ceiling <- function(days = 35) {
  bind_rows(lapply(seq_len(nrow(W_REL)), function(i) {
    age <- W_REL$age_y[i]
    w   <- 3.40*W_REL$rel[i]
    ps  <- approx(PROD_ONTOGENY$age_y, PROD_ONTOGENY$rel, age)$y
    variants <- list(
      G      = list(PRODSCALE = 1,  FOPTMIN = 1,   AGEOFF = 0),
      GP     = list(PRODSCALE = ps, FOPTMIN = 1,   AGEOFF = 0),
      GPO    = list(PRODSCALE = ps, FOPTMIN = 0.30, AGEOFF = age*8760))
    vals <- vapply(variants, function(v) {
      m <- param(mod, c(list(GENO = 0, W0 = w, KW = 0, FIXHB = 15,
                             BGU0 = 0, BGC0 = 0, IRRSET = 30, FBSASET = 0.80,
                             PTSTART = 0, PTSTOP = 1e9, PTDUTY = 12), v))
      o <- as_tibble(mrgsim(m, end = 24*days, delta = 1, hmax = 0.25))
      mean(tail(o$TSBOUT, 24))
    }, numeric(1))
    tibble(age_y = age, W = w, TSB_G = vals[["G"]], TSB_GP = vals[["GP"]],
           TSB_GPO = vals[["GPO"]])
  }))
}
# cn_ceiling() runs the same 12 h/day duty cycle as the reference output
# (PTDUTY = 12).  FIXHB = 15 clamps haemoglobin so that the three nested
# variants differ ONLY in the term being tested: without the clamp the erythron
# drifts during each plateau run and the decomposition is confounded.

# =============================================================================
#  ANALYSIS 6 — AAV8-hUGT1A1 gene transfer in Crigler-Najjar type I
# =============================================================================
gene_therapy_run <- function() {
  m <- param(mod, GENO = 0, GTSTART = 1440, TGXMAX = 0.10,
             IRRSET = 30, FBSASET = 0.80, PTSTART = 24, PTSTOP = 1560)
  as_tibble(mrgsim(m, end = 3360, delta = 12, hmax = 1)) %>%
    mutate(day = time/24) %>%
    select(day, TSBOUT, TGX, BF, UGTPCT)
}

# =============================================================================
#  ANALYSIS 7 — genotype x feeding interaction.  Poor enteral intake raises
#  the peak without touching UGT1A1: the shunt term competes with clearance,
#  so lactation support has the same units as a drug effect.
# =============================================================================
genotype_feeding <- function() {
  gts  <- c(wild = 1.00, het = 0.65, Gilbert = 0.35, star6 = 0.40)
  feed <- list(`human milk`        = list(BREAST = 1, FIEARLY = 1.0,  FILATE = 1.0),
               `formula`           = list(BREAST = 0, FIEARLY = 1.0,  FILATE = 1.0),
               `milk, poor intake` = list(BREAST = 1, FIEARLY = 0.35, FILATE = 0.35))
  bind_rows(lapply(names(gts), function(g) {
    bind_rows(lapply(names(feed), function(f) {
      m <- param(mod, c(list(GENO = gts[[g]]), feed[[f]]))
      o <- as_tibble(mrgsim(m, end = 672, delta = 6, hmax = 0.5))
      tibble(genotype = g, feeding = f,
             peak = max(o$TSBOUT), t_peak = o$time[which.max(o$TSBOUT)],
             d7 = o$TSBOUT[o$time == 168], d14 = o$TSBOUT[o$time == 336],
             d28 = o$TSBOUT[o$time == 672])
    }))
  }))
}

# =============================================================================
#  ANALYSIS 8 — interrupting the enterohepatic shunt as pharmacology
# =============================================================================
ehc_interventions <- function() {
  arms <- list(
    `none`                        = list(par = list(), ev = NULL),
    `oral agar 250 mg/kg/day`     = list(par = list(), ev = dose_agar(3.40)),
    `UDCA 10 mg/kg q12h`          = list(par = list(), ev = dose_udca(3.40)),
    `formula supplementation`      = list(par = list(BREAST = 0), ev = NULL),
    `phenobarbital 5 mg/kg/day`   = list(par = list(), ev = dose_pheno(3.40)))
  bind_rows(lapply(names(arms), function(a) {
    m <- param(mod, arms[[a]]$par)
    o <- if (is.null(arms[[a]]$ev)) mrgsim(m, end = 240, delta = 6, hmax = 0.5)
         else mrgsim(m, events = arms[[a]]$ev, end = 240, delta = 6, hmax = 0.5)
    o <- as_tibble(o)
    tibble(intervention = a, peak = max(o$TSBOUT),
           d7 = o$TSBOUT[o$time == 168], faecal_mg = o$BSTL[o$time == 168])
  }))
}

# =============================================================================
#  PLOTS
# =============================================================================
plot_scenarios <- function(sims) {
  sims %>%
    ggplot(aes(time/24, TSBOUT, colour = scenario)) +
    geom_line(linewidth = 0.7) +
    geom_line(aes(y = THRESHPT), linetype = 2, colour = "grey40",
              show.legend = FALSE) +
    geom_line(aes(y = THRESHET), linetype = 3, colour = "grey20",
              show.legend = FALSE) +
    facet_wrap(~scenario, ncol = 2) +
    labs(x = "postnatal age (days)", y = "total serum bilirubin (mg/dL)",
         title = "Neonatal hyperbilirubinaemia: ten scenarios",
         subtitle = paste("dashed = AAP-2022 phototherapy threshold,",
                          "dotted = exchange threshold (approximation)")) +
    theme_bw(base_size = 10) + theme(legend.position = "none")
}

plot_tsb_vs_free <- function(sims) {
  sims %>%
    select(time, scenario, TSBOUT, BF) %>%
    ggplot(aes(TSBOUT, BF, colour = scenario)) +
    geom_path(linewidth = 0.7) +
    geom_hline(yintercept = 30, linetype = 2) +
    annotate("text", x = 5, y = 33, label = "Bf ~ 30 nM injury threshold",
             size = 3, hjust = 0) +
    labs(x = "total serum bilirubin, what is measured (mg/dL)",
         y = "free bilirubin, what injures (nM)",
         title = "The same TSB is not the same risk",
         subtitle = "each trajectory has its own albumin, affinity and barrier") +
    theme_bw(base_size = 10)
}

plot_photoisomers <- function() {
  photoisomer_run() %>%
    pivot_longer(c(TSBOUT, UCBNAT, ISOMER)) %>%
    mutate(name = recode(name,
                         TSBOUT = "reported TSB (assay)",
                         UCBNAT = "native (4Z,15Z) pigment",
                         ISOMER = "photoisomers")) %>%
    ggplot(aes(time, value, colour = name)) +
    geom_line(linewidth = 0.8) +
    annotate("rect", xmin = 24, xmax = 72, ymin = -Inf, ymax = Inf,
             alpha = 0.08, fill = "blue") +
    labs(x = "postnatal age (h)", y = "mg/dL", colour = NULL,
         title = "The assay counts photoisomers as bilirubin",
         subtitle = "shaded = lamps on; note the DIP at lamps-off as lumirubin washes out") +
    theme_bw(base_size = 10)
}

# =============================================================================
#  DRIVER
# =============================================================================
if (interactive()) {
  sims <- run_all(end = 336)
  print(summarise_scenarios(sims), n = 20, width = 200)
  print(iso_bf_table(), n = 20, width = 200)
  print(pt_dose_response(), n = 20, width = 200)
  print(genotype_feeding(), n = 20, width = 200)
  print(ehc_interventions(), n = 20, width = 200)
  print(gene_therapy_run(), n = 30, width = 200)
  print(plot_scenarios(sims))
  print(plot_tsb_vs_free(sims))
  print(plot_photoisomers())
}

# =============================================================================
#  END.  Cross-check every number against nhb_reference_output.txt, which was
#  produced by the independent implementation in nhb_reference_check.py.
# =============================================================================
