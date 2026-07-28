## =============================================================================
##  ERYTHROPOIETIC PROTOPORPHYRIA (EPP) / X-LINKED PROTOPORPHYRIA (XLP)
##  Quantitative Systems Pharmacology model — mrgsolve implementation
## =============================================================================
##
##  THE IDEA THE MODEL IS BUILT AROUND
##  ----------------------------------
##  Protoporphyrin IX is, uniquely, BOTH the product of a biosynthetic pathway
##  AND the substrate of the enzyme that is broken.  Everything clinical is
##  therefore a RATIO — supply divided by disposal — and the SIGN of any
##  intervention is a property of WHICH ARM IS RATE-LIMITING, not of the drug.
##
##  Three results that are computed rather than assumed fall out of that single
##  structural fact.  (Every number below is reproduced by the dependency-free
##  Python mirror `epp_reference_check.py`; see `epp_reference_output.txt`.)
##
##  (1) THE ~35% PENETRANCE THRESHOLD IS DERIVED, NOT TYPED IN.
##      Normal ferrochelatase CAPACITY is set to 2.9x the normal ALAS2 flux
##      (VFECH x Fe-saturation / vALAS = 2.900).  The terminal step therefore
##      stops keeping up below FRES* = 1/2.9 = 34.5% residual activity — which
##      is the ~35% threshold that makes EPP behave recessively despite two
##      dominant-looking alleles.  Erythrocyte PPIX over the FRES sweep:
##          FECH 100%  0.13 umol/L      FECH 30%   6.12
##          FECH  50%  0.53             FECH 20%  16.11
##          FECH  40%  1.24             FECH 15%  21.42   <- typical EPP
##          FECH  35%  2.66  <- knee    FECH 10%  26.77
##
##  (2) THE ZINC-PROTOPORPHYRIN SIGNATURE SEPARATES EPP FROM XLP BY ITSELF.
##      FECH is also the enzyme that inserts ZINC.  Break it (EPP) and the
##      substrate piles up while the metal cannot go in either; leave it intact
##      and merely flood it (XLP) and zinc insertion rises with the substrate:
##          normal   free PPIX  0.13   Zn-PP 0.85   metal-free 13%
##          EPP      free PPIX 21.42   Zn-PP 1.42   metal-free 94%   (x164 / x1.7)
##          XLP      free PPIX 17.49   Zn-PP 8.13   metal-free 68%   (x134 / x9.6)
##      Nothing in the code special-cases the two diseases: only FRES and GOF
##      differ between those two runs.
##
##  (3) THE TWO DRUG AXES MULTIPLY, THEY DO NOT ADD.
##      Photodynamic dose D = E x T(melanin) x [PPIX]skin x t is a PRODUCT, so
##      tolerance time t_tol = D_pro/(E x T x PPIX) responds to the product of
##      the two available interventions:
##          afamelanotide 16 mg q60d   t_tol x1.93   (melanin up, PPIX UNCHANGED)
##          bitopertin 60 mg od        t_tol x1.59   (PPIX -39%, melanin UNCHANGED)
##          additive expectation       t_tol x2.52
##          Bliss/multiplicative       t_tol x3.07
##          OBSERVED in the model      t_tol x3.08   (+22% over additive)
##      A 180-day season: 89 pain-free sun hours untreated, 166 on afamelanotide,
##      256 on the combination.
##
##  (4) PROTOPORPHYRIC HEPATOPATHY IS A SADDLE-NODE, NOT AN ACCUMULATION.
##      PPIX has exactly one exit (bile), PPIX crystals reduce bile flow, and
##      reduced bile flow raises hepatic PPIX.  Positive feedback on a
##      single-exit system gives MORE THAN ONE steady state.  Continuation over
##      residual FECH activity locates two folds:
##          FECH > 16.1%  monostable healthy
##          16.1% - 7.4%  BISTABLE (e.g. at FECH 10%: healthy hepatic PPIX 0.88,
##                        saddle at 1.46, disease branch at 5.32)
##          FECH < 7.4%   the healthy branch CEASES TO EXIST
##      A typical EPP patient sits on the healthy branch with a finite margin to
##      a threshold — which is why protoporphyric liver failure is abrupt and
##      usually triggered, and why the lifetime incidence is a few percent
##      rather than universal.
##
##  (5) THE PATIENT IS A FEEDBACK CONTROLLER AND THE LOOP DELAY IS THE DISEASE.
##      Prodrome = sensor, retreat = actuator, human reaction time = loop delay.
##      One identical sunny day, one identical PPIX, varying only the delay:
##          2 min  -> peak pain NRS 1.2      20 min -> NRS 3.5
##          10 min -> NRS 1.6                40 min -> NRS 8.3 (peaks 18 h later)
##      A shielding drug does not only buy minutes of sunlight: by slowing the
##      RATE of dose accrual it makes the SAME human delay non-limiting.  That
##      benefit is invisible to an "hours in sunlight" endpoint.
##
##  MODEL SIZE
##  ----------
##  37 ODE compartments: heme biosynthesis (9) - erythrocyte pools (2) -
##  distribution/hepatobiliary (6) - photobiology & injury cascade (7) -
##  MC1R/melanin (2) - drug PK (9) - accumulators (2).
##  9 built-in therapeutic scenarios (see the driver at the bottom).
##
##  UNITS: time = HOURS.  Porphyrins = umol/L of the stated compartment.
##  Glycine = mmol/L.  Iron = umol/L.  Drug amounts = ug, concentrations ug/L.
##
##  CALIBRATION NOTES (targets, and where they come from)
##  -----------------------------------------------------
##  * normal erythrocyte protoporphyrin < 1.5 umol/L, predominantly ZINC-
##    chelated; model gives total 0.98 umol/L, 87% Zn-PP.
##  * clinically manifest EPP 10-40 umol/L, >85% metal-free; model gives
##    21.4 umol/L, 94% metal-free at FRES = 0.15.
##  * sunlight tolerance to first prodrome of minutes, not hours; model gives
##    10.2 min of direct summer sun untreated.
##  * afamelanotide roughly doubles time in sunlight (Langendonk et al., NEJM
##    2015: ~116 vs ~61 h over the European 270-day study); model gives a
##    season ratio of 166/89 = 1.87.
##  * bitopertin lowers whole-blood metal-free PPIX by roughly 40% (BEACON /
##    APOLLO phase 2); model gives -39% at 60 mg once daily.
##  * beta-carotene is included and is deliberately WEAK (t_tol x1.09): its
##    quenching must compete with a ~3 us singlet-oxygen lifetime at plasma
##    concentrations of 1-3 mg/L, and the clinical evidence is correspondingly
##    unimpressive.  A model that made it work would be the wrong model.
##
##  IMPORTANT: this is an educational / research QSP model. It is not validated
##  for clinical decision-making.  See ../README.md disclaimer.
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

code <- '
$PROB
# Erythropoietic Protoporphyria / X-Linked Protoporphyria QSP model
# 37 ODEs. FECH loss-of-function or ALAS2 gain-of-function -> PPIX overflow
# -> 400-410 nm photoexcitation -> singlet oxygen -> burning pain.

$PARAM @annotated
// ---------------- genotype -------------------------------------------------
FRES   : 0.15    : Residual ferrochelatase activity (fraction of normal)
GOF    : 1.0     : ALAS2 gain-of-function multiplier (XLP ~3.5)
FMC1R  : 1.0     : MC1R functionality (~0.3 for red-hair loss-of-function)

// ---------------- erythroid glycine / ALAS2 --------------------------------
KINGLY : 4.0     : Erythroblast glycine influx via GlyT1/SLC6A9 (mmol/L/h)
KOUTGLY: 5.0     : Glycine consumption other than ALAS2 (1/h)
KMGLY  : 0.8     : ALAS2 Km for glycine - half-saturated at baseline (mmol/L)
VALAS  : 201.4   : ALAS2 Vmax (umol/L erythron/h)
KHFB   : 1200    : Heme feedback constant on ALAS2 (weak in the erythron)
KIRE   : 12.0    : IRE/IRP half-point on the labile iron pool (umol/L)
NIRE   : 2.0     : IRE/IRP Hill coefficient

// ---------------- pathway transfer (all enzymes in large excess) -----------
KALAD  : 6.0     : ALA dehydratase transfer (1/h)
KHMBS  : 6.0     : PBG deaminase transfer (1/h)
KUROD  : 6.0     : Uroporphyrinogen decarboxylase transfer (1/h)
KCPOX  : 6.0     : Coproporphyrinogen oxidase transfer (1/h)
KPPOX  : 6.0     : Protoporphyrinogen oxidase transfer (1/h)

// ---------------- ferrochelatase -------------------------------------------
VFECH  : 435.0   : FECH Vmax; normal CAPACITY = 2.9 x normal flux (umol/L/h)
KMPPIX : 1.0     : FECH Km for protoporphyrin IX (umol/L)
KMFE   : 10.0    : FECH Km for Fe2+ (umol/L)
KISC   : 20.0    : 2Fe-2S cluster assembly half-point on labile iron (umol/L)
VZNREL : 0.02    : Zn-insertion efficiency relative to Fe-insertion
KMZN   : 5.0     : FECH Km for Zn2+ (umol/L)
ZN     : 8.0     : Erythroblast labile zinc (umol/L)
KIFZ   : 3.0     : Fe competition constant at the FECH metal site (umol/L)
KESC   : 0.672   : PPIX escape erythron -> plasma (1/h) [THE OVERFLOW ROUTE]
KHEMEO : 1.0     : Heme utilisation (1/h)

// ---------------- iron ------------------------------------------------------
FE0    : 20.0    : Reference labile iron pool (umol/L)
KFEIN  : 20.0    : Baseline iron supply (umol/L/h)
KFEOUT : 1.0     : Iron efflux / storage (1/h)
FEDOSE : 0.0     : Iron supplementation as a fraction of KFEIN

// ---------------- erythrocyte pools -----------------------------------------
KPROD  : 2.4056E-4 : Erythropoietic release, ln2/(120 d) (1/h)
KSEN   : 2.4056E-4 : Erythrocyte senescence (1/h)
KEFFL  : 7.2E-4    : PPIX efflux out of circulating erythrocytes (1/h)
ALPHA  : 1.0       : Erythroblast -> erythrocyte PPIX partition
BETAZN : 8.09E-4   : Zn-PP packaging coefficient
KZN2   : 8.0E-5    : Late-reticulocyte Zn insertion (requires INTACT FECH)
ERYSUP : 1.0       : Erythropoietic suppression (transfusion < 1)

// ---------------- distribution ----------------------------------------------
KTR    : 0.0521  : Erythron escape -> plasma scaling
KPLLIV : 3.0     : Plasma -> hepatocyte uptake (1/h)
KLIVPL : 0.1     : Hepatocyte -> plasma back-flux (1/h)
KPLSK  : 2.0     : Plasma -> skin (1/h)
KSKPL  : 1.667   : Skin -> plasma (1/h); skin:plasma partition 1.2
KEHC   : 0.05    : Enterohepatic reabsorption (1/h)
KFEC   : 0.30    : Faecal loss from the gut lumen (1/h)
CHOLEFF: 0.0     : Cholestyramine effect on reabsorption (0-0.9)
KPHER  : 0.0     : Plasmapheresis clearance (1/h)

// ---------------- hepatobiliary ---------------------------------------------
VBILE  : 20.0    : Canalicular Vmax - the SATURABLE choke point (umol/L/h)
KMBILE : 3.0     : Canalicular Km (umol/L)
KCHIN  : 0.010   : Cholestasis induction (1/h)
KCHOUT : 0.004   : Cholestasis resolution (1/h)
KCH    : 1.8     : Hepatic PPIX at which crystals drive cholestasis (umol/L)
NCH    : 6.0     : Hill coefficient of crystal -> cholestasis
FCHOL  : 0.92    : Maximum fractional loss of bile flow
KINJ   : 0.012   : Hepatocyte injury accrual (1/h)
KINJOUT: 0.010   : Hepatocyte injury resolution (1/h)
KILIV  : 5.0     : Hepatic PPIX scale for direct hepatotoxicity (umol/L)

// ---------------- photobiology -----------------------------------------------
ODMEL  : 0.693   : Epidermal optical density per melanin unit at 405 nm
KPHOT  : 8.20    : Photodynamic dose rate constant
EREF   : 1.0     : Reference irradiance for the tolerance-time read-out
EMAX   : 1.0     : Peak ambient violet irradiance (normalised)
FCLOUD : 1.0     : Cloud transmission at 400-410 nm (overcast ~0.5)
FGLASS : 1.0     : Window-glass transmission at 405 nm (~0.9 - it does NOT block)
FCLOTH : 1.0     : Clothing / opaque screen transmission (0.05 = full cover)
DPRO   : 1.0     : Prodrome threshold (dose units)
DCRIT  : 3.0     : Full phototoxic-reaction threshold (dose units)
KDB    : 0.15    : Decay of the prodromal dose memory (1/h)
IC50BC : 12000   : Beta-carotene conc. for 50% singlet-oxygen quenching (ug/L)
KOXOUT : 0.70    : Oxidative-signal clearance (1/h)
OUTS   : 10.0    : Start of the intended outdoor window (hour of day)
OUTE   : 18.0    : End of the intended outdoor window (hour of day)
EXPOSE : 1.0     : Master switch for sun exposure (0 = strict avoidance)
TAUDEL : 0.033   : Behavioural retreat delay (h); ~2 min adult, ~0.7 h child
NBOUT  : 3.0     : Outdoor bouts per day (for the censored endpoint)
WINDOW : 8.0     : 10:00-18:00 registrational window (h)

// ---------------- injury cascade ---------------------------------------------
KMAST  : 1.0     : Mast-cell degranulation gain
KMASTO : 0.25    : Mast-cell signal decay (1/h)
KC5    : 1.0     : Complement activation gain
KC5O   : 0.35    : Complement signal decay (1/h)
KED    : 0.05    : Oedema formation gain
KEDO   : 0.06    : Oedema resolution (1/h)
KNOC   : 0.22    : Nociceptor sensitisation gain
KNOCO  : 0.020   : Nociceptor desensitisation (1/h) - the days-long tail
KPAIN  : 6.0     : Nociceptor drive giving NRS 5 (Hill, n = 2)
AMAX   : 8.0     : Maximum supra-threshold inflammatory amplification
KAMP   : 3.0     : Dose at which amplification is half-maximal
NAMP   : 4.0     : Steepness of the all-or-nothing reaction gate

// ---------------- MC1R / melanin ----------------------------------------------
EMAXMC : 12.0    : Maximal fold-increase in melanogenic drive over basal
EC50AF : 0.30    : Afamelanotide plasma EC50 (ug/L)
EC50DE : 30000   : Dersimelagon plasma EC50 (ug/L)
KTYRIN : 4.13E-3 : Tyrosinase/MITF program formation (1/h), t1/2 7 d
KTYRO  : 4.13E-3 : Tyrosinase/MITF program decay (1/h)
KMELIN : 8.25E-4 : Melanin formation (1/h), t1/2 35 d
KMELO  : 8.25E-4 : Melanin turnover (1/h)

// ---------------- drug PK -------------------------------------------------------
KREL   : 0.030   : Afamelanotide implant release (1/h) - ~5 days
FAFA   : 1.0     : Afamelanotide implant bioavailability
VAFA   : 20.0    : Afamelanotide volume of distribution (L)
KEAFA  : 1.386   : Afamelanotide elimination (1/h), t1/2 ~30 min
KADER  : 0.50    : Dersimelagon absorption (1/h)
VDER   : 200.0   : Dersimelagon V/F (L)
KEDER  : 0.0289  : Dersimelagon elimination (1/h), t1/2 ~24 h
KABIT  : 0.40    : Bitopertin absorption (1/h)
VBIT   : 700.0   : Bitopertin V/F (L)
KEBIT  : 0.01386 : Bitopertin elimination (1/h), t1/2 ~50 h
IMAXBIT: 0.85    : Maximum GlyT1 inhibition of glycine influx
IC50BIT: 328.0   : Bitopertin IC50 at GlyT1 (ug/L)
KABC   : 0.05    : Beta-carotene absorption (1/h)
VBC    : 3000.0  : Beta-carotene V/F (L) - low bioavailability, large V
KEBC   : 0.00206 : Beta-carotene elimination (1/h), t1/2 ~14 d
KEHEM  : 0.10    : Hemin elimination (1/h)
IMAXHEM: 0.35    : Maximum ALAS suppression by hemin (weak: ALAS1 >> ALAS2)
IC50HEM: 5.0     : Hemin IC50 (ug/L)

$CMT @annotated
GLY   : Erythroblast glycine (mmol/L)
ALA   : 5-aminolevulinic acid (umol/L)
PBG   : Porphobilinogen (umol/L)
UPG   : Uroporphyrinogen III (umol/L)
CPG   : Coproporphyrinogen III (umol/L)
PPG   : Protoporphyrinogen IX (umol/L)
PPIXE : Erythroblast protoporphyrin IX (umol/L)
HEME  : Erythroid heme (umol/L)
FE    : Labile iron pool (umol/L)
PRBC  : Erythrocyte metal-free PPIX (umol/L)   [THE clinical biomarker]
ZRBC  : Erythrocyte zinc protoporphyrin (umol/L)
PPL   : Plasma PPIX (umol/L)
PSK   : Dermal / perivascular PPIX (umol/L)    [THE photoactive pool]
PLIV  : Hepatocyte PPIX (umol/L)
PGUT  : Intestinal luminal PPIX (umol/L)
CHOL  : Cholestasis index (0-1)
LINJ  : Hepatocyte injury index
OX    : Oxidative / 4-HNE damage signal
MAST  : Mast-cell degranulation signal
C5A   : Complement C3a/C5a signal
EDEMA : Erythema / oedema
NOCI  : Nociceptor sensitisation
AVOID : Fraction of time in shade (behavioural actuator, 0-1)
DBOUT : Prodromal dose memory (the patient SENSOR)
TYR   : Tyrosinase / MITF melanogenic program
MEL   : Epidermal eumelanin optical density (1 = untreated)
ADEP  : Afamelanotide implant depot (ug)
AC    : Afamelanotide plasma (ug/L)
DG    : Dersimelagon gut (ug)
DC    : Dersimelagon plasma (ug/L)
BG    : Bitopertin gut (ug)
BC    : Bitopertin plasma (ug/L)
CG    : Beta-carotene gut (ug)
CC    : Beta-carotene plasma (ug/L)
HEMC  : Hemin plasma (ug/L)
SUNCUM: Cumulative pain-free sun exposure (h)
RXNCUM: Cumulative phototoxic reaction-days

$MAIN
// baseline initial conditions (the driver overwrites these with the solved
// steady state; these are only a sensible starting point)
GLY_0   = 0.8;
FE_0    = FE0;
HEME_0  = 100.0;
TYR_0   = 1.0;
MEL_0   = 1.0;

$ODE
// ---- guards ---------------------------------------------------------------
double gly  = (GLY   > 1e-9) ? GLY   : 1e-9;
double fe   = (FE    > 1e-9) ? FE    : 1e-9;
double ppix = (PPIXE > 0.0)  ? PPIXE : 0.0;
double pliv = (PLIV  > 0.0)  ? PLIV  : 0.0;
double psk  = (PSK   > 0.0)  ? PSK   : 0.0;

// ===========================================================================
// 1. GLYCINE  — the substrate arm, and the bitopertin target
// ===========================================================================
double IBIT = IMAXBIT*BC/(IC50BIT + BC);
dxdt_GLY = KINGLY*(1.0 - IBIT) - KOUTGLY*gly;

// ===========================================================================
// 2. ALAS2 — substrate x iron(IRE/IRP) x heme feedback x genotype
//    The IRE/IRP term is the (+) arm of the iron paradox: iron REPLETION
//    de-represses ALAS2 translation and so raises the supply of substrate.
// ===========================================================================
double GLYSAT = gly/(KMGLY + gly);
double IREF   = (pow(fe,NIRE)/(pow(KIRE,NIRE) + pow(fe,NIRE))) /
                (pow(FE0,NIRE)/(pow(KIRE,NIRE) + pow(FE0,NIRE)));
double FHEME  = 1.0/(1.0 + pow(HEME/KHFB, 2.0));
double IHEM   = IMAXHEM*HEMC/(IC50HEM + HEMC);
double vALAS  = VALAS*GOF*GLYSAT*IREF*FHEME*(1.0 - IHEM)*ERYSUP;

dxdt_ALA = vALAS       - KALAD*ALA;
dxdt_PBG = KALAD*ALA   - KHMBS*PBG;
dxdt_UPG = KHMBS*PBG   - KUROD*UPG;
dxdt_CPG = KUROD*UPG   - KCPOX*CPG;
dxdt_PPG = KCPOX*CPG   - KPPOX*PPG;

// ===========================================================================
// 3. FERROCHELATASE — ONE enzyme, TWO metals.
//    Iron enters the (-) disposal arm twice: as the co-substrate (KMFE) and
//    through 2Fe-2S cluster availability (KISC), which sets Vmax.  The net
//    sign of an iron intervention is therefore emergent, not coded.
//    Zinc insertion is the SAME enzyme, which is why zinc protoporphyrin
//    reports on residual FECH activity rather than on the size of the overload.
// ===========================================================================
double FISC  = (fe/(KISC + fe))/(FE0/(KISC + FE0));
double fFECH = FRES*FISC;
double sat   = ppix/(KMPPIX + ppix);
double vFe   = VFECH*fFECH*sat*(fe/(KMFE + fe));
double vZn   = VFECH*fFECH*VZNREL*sat*(ZN/(KMZN + ZN))*(KIFZ/(KIFZ + fe));
double Jesc  = KESC*ppix;              // the OVERFLOW route out of the erythron

dxdt_PPIXE = KPPOX*PPG - vFe - vZn - Jesc;
dxdt_HEME  = vFe - KHEMEO*HEME;
dxdt_FE    = KFEIN*(1.0 + FEDOSE) - KFEOUT*fe;

// ===========================================================================
// 4. ERYTHROCYTE POOLS
//    Zn-PP has two sources: enzymatic insertion during erythropoiesis, and
//    late-reticulocyte insertion, which needs an INTACT enzyme (KZN2*fFECH).
//    That second term is what makes XLP zinc-rich and EPP metal-free.
// ===========================================================================
dxdt_PRBC = KPROD*ALPHA*ppix*ERYSUP - (KSEN + KEFFL)*PRBC;
dxdt_ZRBC = BETAZN*vZn + KZN2*fFECH*PRBC - KSEN*ZRBC;

// ===========================================================================
// 5. DISTRIBUTION — a hydrophobic molecule with ZERO renal clearance and
//    exactly ONE exit: bile.  That single-exit topology is what puts the
//    catastrophic complication in the liver.
// ===========================================================================
double reabs = KEHC*(1.0 - CHOLEFF)*PGUT;
dxdt_PPL = KTR*Jesc + KEFFL*PRBC + KSKPL*psk + reabs + KLIVPL*pliv
           - (KPLLIV + KPLSK + KPHER)*PPL;
dxdt_PSK = KPLSK*PPL - KSKPL*psk;

double BF   = 1.0 - FCHOL*CHOL;                       // bile-flow factor
double vcan = VBILE*BF*pliv/(KMBILE + pliv);          // SATURABLE
dxdt_PLIV = KPLLIV*PPL - vcan - KLIVPL*pliv;
dxdt_PGUT = vcan - KFEC*PGUT - reabs;

// ---- the vicious cycle: crystals -> cholestasis -> less exit -> crystals ---
double hillc = pow(pliv,NCH)/(pow(KCH,NCH) + pow(pliv,NCH));
dxdt_CHOL = KCHIN*hillc*(1.0 - CHOL) - KCHOUT*CHOL;
dxdt_LINJ = KINJ*(CHOL*CHOL + pow(pliv/KILIV, 2.0)) - KINJOUT*LINJ;

// ===========================================================================
// 6. MC1R -> EUMELANIN — the only axis that changes the OPTICAL DENOMINATOR.
//    Note that the melanogenic program (TYR) has a 7-day half-life and melanin
//    a 35-day one, while afamelanotide plasma has a 30-MINUTE half-life. The
//    PD therefore outlasts the PK by more than an order of magnitude, which is
//    exactly why a 5-day implant is dosed every 60 days.
// ===========================================================================
double occ = AC/EC50AF + DC/EC50DE;
double SIG = FMC1R*EMAXMC*occ/(1.0 + occ);
dxdt_TYR = KTYRIN*(1.0 + SIG) - KTYRO*TYR;
dxdt_MEL = KMELIN*TYR - KMELO*MEL;

// ===========================================================================
// 7. PHOTOBIOLOGY — the disease is, at bottom, an absorption spectrum.
//    Ambient irradiance x epidermal transmittance x photosensitiser
//    concentration.  A PRODUCT — which is why interventions multiply.
// ===========================================================================
double TEPI = exp(-ODMEL*MEL);
double QUEN = 1.0/(1.0 + CC/IC50BC);
double tod  = fmod(SOLVERTIME, 24.0);
double sarg = sin(3.14159265358979*(tod - 6.0)/12.0);
double sun  = (tod >= 6.0 && tod <= 18.0 && sarg > 0.0) ? pow(sarg, 1.3) : 0.0;
double fout = (EXPOSE > 0.5 && tod >= OUTS && tod <= OUTE) ? 1.0 : 0.0;
double Eamb = EMAX*sun*FCLOUD*FGLASS*FCLOTH;
double R    = KPHOT*Eamb*TEPI*psk*QUEN*fout*(1.0 - AVOID);

// ---- the behavioural control loop -----------------------------------------
// DBOUT is the SENSOR the patient actually has: it integrates dose with a
// multi-hour memory, which is why a single prodrome ends the whole outing
// rather than a minute of it.  TAUDEL is the LOOP DELAY, and the dose accrued during it decides
// prodrome-only versus a three-day reaction.
dxdt_DBOUT = R - KDB*DBOUT;
double SW  = 1.0/(1.0 + exp(-(DBOUT - DPRO)/0.08));
dxdt_AVOID = (SW - AVOID)/TAUDEL;

// ===========================================================================
// 8. CUTANEOUS INJURY CASCADE
//    Singlet oxygen lives ~3 us; the pain lasts ~3 days.  The persistence is
//    in the slow state variables below, not in the photon.
//    Above DCRIT the mast-cell/complement arm becomes self-amplifying, which
//    is what makes the full phototoxic reaction all-or-nothing.
// ===========================================================================
dxdt_OX = R - KOXOUT*OX;
double dbn = pow(DBOUT/KAMP, NAMP);
double AMP = 1.0 + AMAX*dbn/(1.0 + dbn);
dxdt_MAST  = KMAST*OX*AMP - KMASTO*MAST;
dxdt_C5A   = KC5*OX*AMP   - KC5O*C5A;
dxdt_EDEMA = KED*(MAST + C5A) - KEDO*EDEMA;
dxdt_NOCI  = KNOC*(1.5*OX + 0.5*EDEMA) - KNOCO*NOCI;

// ===========================================================================
// 9. ACCUMULATORS — the registrational endpoints
// ===========================================================================
dxdt_SUNCUM = fout*(1.0 - AVOID);
dxdt_RXNCUM = (1.0/24.0)/(1.0 + exp(-(DBOUT - DCRIT)/0.20));

// ===========================================================================
// 10. DRUG PK
// ===========================================================================
dxdt_ADEP = -KREL*ADEP;
dxdt_AC   = KREL*ADEP*FAFA/VAFA - KEAFA*AC;
dxdt_DG   = -KADER*DG;
dxdt_DC   = KADER*DG/VDER - KEDER*DC;
dxdt_BG   = -KABIT*BG;
dxdt_BC   = KABIT*BG/VBIT - KEBIT*BC;
dxdt_CG   = -KABC*CG;
dxdt_CC   = KABC*CG/VBC - KEBC*CC;
dxdt_HEMC = -KEHEM*HEMC;

$TABLE
double TOTEP   = PRBC + ZRBC;                       // total erythrocyte porphyrin
double ZNFRAC  = (TOTEP > 0) ? 100.0*ZRBC/TOTEP : 0.0;
double MFFRAC  = 100.0 - ZNFRAC;                    // metal-free % (diagnostic)
double TEPIo   = exp(-ODMEL*MEL);
double QUENo   = 1.0/(1.0 + CC/IC50BC);
double PSKo    = (PSK > 1e-12) ? PSK : 1e-12;
double TTOL_H  = DPRO/(KPHOT*EREF*TEPIo*PSKo*QUENo);
double TTOLMIN = 60.0*TTOL_H;                       // minutes to first prodrome
double SUNHRD  = WINDOW*(1.0 - exp(-NBOUT*TTOL_H/WINDOW));   // censored endpoint
double PAIN    = 10.0*NOCI*NOCI/(KPAIN*KPAIN + NOCI*NOCI);
double ALT     = 25.0*(1.0 + 3.0*LINJ);
double TBIL    = 0.6*(1.0 + 6.0*CHOL*LINJ);
double MELIDX  = MEL;
double BILEFLOW= 1.0 - FCHOL*CHOL;

$CAPTURE @annotated
TOTEP   : Total erythrocyte protoporphyrin (umol/L)
ZNFRAC  : Zinc protoporphyrin fraction (%)
MFFRAC  : Metal-free protoporphyrin fraction (%) - EPP >85, XLP 50-85
TEPIo   : Epidermal transmittance at 405 nm
TTOLMIN : Sunlight tolerance time to first prodrome (min)
SUNHRD  : Pain-free sun exposure per day, censored at the window (h)
PAIN    : Pain numerical rating scale 0-10
ALT     : Alanine aminotransferase (U/L)
TBIL    : Total bilirubin (mg/dL)
MELIDX  : Epidermal melanin optical density (1 = untreated)
BILEFLOW: Bile flow as a fraction of normal
'

mod <- mcode("epp_qsp", code)

## =============================================================================
##  STEADY STATE
##  The slow pools (erythrocyte PPIX t1/2 ~40 d, melanin ~35 d, cholestasis
##  ~7 d) take a long time to settle, so every scenario is started from a
##  pre-equilibrated state rather than from the nominal initials.
##  NOTE the second argument: equilibration must be run with EXPOSE = 0, so
##  that the biochemical steady state is not contaminated by sunlight.
## =============================================================================
epp_steady <- function(mod, pars = list(), days = 1200) {
  m <- mod %>% param(pars) %>% param(EXPOSE = 0)
  out <- m %>% mrgsim(end = days * 24, delta = 24, hmax = 0.5)
  tail(as.data.frame(out), 1)
}

init_from <- function(mod, ss) {
  cmts <- names(mod@init@data)
  vals <- as.list(ss[, intersect(cmts, names(ss)), drop = FALSE])
  mod %>% init(vals)
}

## =============================================================================
##  DOSING EVENT BUILDERS
## =============================================================================
ev_afamelanotide <- function(n = 3, interval_d = 60, dose_mg = 16)
  ev(amt = dose_mg * 1000, cmt = "ADEP", ii = interval_d * 24, addl = n - 1)

ev_dersimelagon  <- function(days = 180, dose_mg = 300)
  ev(amt = dose_mg * 1000, cmt = "DG", ii = 24, addl = days - 1)

ev_bitopertin    <- function(days = 180, dose_mg = 60)
  ev(amt = dose_mg * 1000, cmt = "BG", ii = 24, addl = days - 1)

ev_betacarotene  <- function(days = 180, dose_mg = 180)
  ev(amt = dose_mg * 1000, cmt = "CG", ii = 24, addl = days - 1)

ev_hemin         <- function(n = 4, dose_mg = 250)
  ev(amt = dose_mg * 1000, cmt = "HEMC", ii = 24, addl = n - 1)

## =============================================================================
##  SCENARIO 1 — UNTREATED EPP THROUGH A SUMMER
##  Establishes the reference: ~21 umol/L erythrocyte PPIX, ~10 min of direct
##  summer sun before the prodrome, ~0.5 pain-free sun hours a day.
## =============================================================================
scenario_untreated <- function(mod) {
  ss <- epp_steady(mod, list(FRES = 0.15))
  init_from(mod, ss) %>%
    param(FRES = 0.15, EXPOSE = 1) %>%
    mrgsim(end = 180 * 24, delta = 6, hmax = 0.05) %>%
    as.data.frame() %>% mutate(scenario = "1. Untreated EPP")
}

## =============================================================================
##  SCENARIO 2 — AFAMELANOTIDE 16 mg IMPLANT q60d x 3
##  The PK/PD hysteresis is the point: plasma is gone by day 5, protection is
##  still rising at day 30, and the dose interval is set by melanin turnover.
##  Erythrocyte PPIX must NOT move — that is the model's falsifiable claim.
## =============================================================================
scenario_afamelanotide <- function(mod) {
  ss <- epp_steady(mod, list(FRES = 0.15))
  init_from(mod, ss) %>%
    param(FRES = 0.15, EXPOSE = 1) %>%
    mrgsim(events = ev_afamelanotide(3, 60, 16),
           end = 180 * 24, delta = 6, hmax = 0.05) %>%
    as.data.frame() %>% mutate(scenario = "2. Afamelanotide 16 mg q60d")
}

## =============================================================================
##  SCENARIO 3 — DERSIMELAGON 300 mg ORALLY ONCE DAILY
##  Same receptor, continuous rather than pulsed stimulation: no sawtooth, but
##  a smaller time-averaged melanin optical density than the implant.
## =============================================================================
scenario_dersimelagon <- function(mod, dose_mg = 300) {
  ss <- epp_steady(mod, list(FRES = 0.15))
  init_from(mod, ss) %>%
    param(FRES = 0.15, EXPOSE = 1) %>%
    mrgsim(events = ev_dersimelagon(180, dose_mg),
           end = 180 * 24, delta = 6, hmax = 0.05) %>%
    as.data.frame() %>%
    mutate(scenario = sprintf("3. Dersimelagon %d mg od", dose_mg))
}

## =============================================================================
##  SCENARIO 4 — BITOPERTIN 60 mg ORALLY ONCE DAILY
##  The other axis entirely: GlyT1 blockade starves ALAS2 of glycine, lowering
##  the photosensitiser itself.  Melanin must NOT move.
## =============================================================================
scenario_bitopertin <- function(mod, dose_mg = 60) {
  ss <- epp_steady(mod, list(FRES = 0.15))
  init_from(mod, ss) %>%
    param(FRES = 0.15, EXPOSE = 1) %>%
    mrgsim(events = ev_bitopertin(180, dose_mg),
           end = 180 * 24, delta = 6, hmax = 0.05) %>%
    as.data.frame() %>%
    mutate(scenario = sprintf("4. Bitopertin %d mg od", dose_mg))
}

## =============================================================================
##  SCENARIO 5 — COMBINATION, AND THE BLISS TEST
##  Because photodynamic dose is a PRODUCT, shielding and source-reduction
##  compose multiplicatively.  This scenario exists to be compared against the
##  additive expectation computed in `combination_index()` below.
## =============================================================================
scenario_combination <- function(mod) {
  ss <- epp_steady(mod, list(FRES = 0.15))
  init_from(mod, ss) %>%
    param(FRES = 0.15, EXPOSE = 1) %>%
    mrgsim(events = c(ev_afamelanotide(3, 60, 16), ev_bitopertin(180, 60)),
           end = 180 * 24, delta = 6, hmax = 0.05) %>%
    as.data.frame() %>% mutate(scenario = "5. Afamelanotide + bitopertin")
}

combination_index <- function(mod) {
  season_mean <- function(df) mean(df$TTOLMIN)
  base <- season_mean(scenario_untreated(mod))
  fa   <- season_mean(scenario_afamelanotide(mod)) / base
  fb   <- season_mean(scenario_bitopertin(mod))    / base
  fc   <- season_mean(scenario_combination(mod))   / base
  data.frame(
    metric = c("afamelanotide alone", "bitopertin alone",
               "ADDITIVE expectation (fa + fb - 1)",
               "MULTIPLICATIVE / Bliss (fa * fb)",
               "OBSERVED combination",
               "excess over additive (%)"),
    fold   = c(fa, fb, fa + fb - 1, fa * fb, fc,
               100 * (fc / (fa + fb - 1) - 1))
  )
}

## =============================================================================
##  SCENARIO 6 — X-LINKED PROTOPORPHYRIA
##  Identical equations; only ALAS2 gain-of-function (GOF) and intact FECH
##  (FRES = 1) differ.  The zinc-protoporphyrin fraction separates it from EPP
##  without anything in the model being told that these are different diseases.
## =============================================================================
scenario_xlp <- function(mod) {
  ss <- epp_steady(mod, list(FRES = 1.0, GOF = 3.5))
  init_from(mod, ss) %>%
    param(FRES = 1.0, GOF = 3.5, EXPOSE = 1) %>%
    mrgsim(end = 180 * 24, delta = 6, hmax = 0.05) %>%
    as.data.frame() %>% mutate(scenario = "6. X-linked protoporphyria")
}

## =============================================================================
##  SCENARIO 7 — THE IRON SIGN FLIP
##  Iron supplementation given to EPP and to XLP.  Same nutrient, same
##  equations, opposite usefulness: in XLP the disposal arm is intact and wins;
##  in EPP it is broken, so all that is left is the IRE/IRP supply arm.
## =============================================================================
scenario_iron <- function(mod, fedose = 0.5) {
  grid <- expand.grid(FRES = c(0.15, 1.0), FEDOSE = c(0, fedose))
  out <- lapply(seq_len(nrow(grid)), function(i) {
    gof <- if (grid$FRES[i] > 0.5) 3.5 else 1.0
    ss  <- epp_steady(mod, list(FRES = grid$FRES[i], GOF = gof,
                                FEDOSE = grid$FEDOSE[i]))
    data.frame(disease = if (gof > 1) "XLP" else "EPP",
               FEDOSE  = grid$FEDOSE[i],
               PPIX    = ss$PRBC, ZNPP = ss$ZRBC,
               ZNFRAC  = ss$ZNFRAC, TTOLMIN = ss$TTOLMIN)
  })
  do.call(rbind, out)
}

## =============================================================================
##  SCENARIO 8 — PROTOPORPHYRIC HEPATOPATHY AND ITS RESCUE
##  Continuation over residual FECH activity locates the two folds of the
##  cholestasis feedback loop; then a patient who has fallen onto the
##  cholestatic branch is treated.  Liver transplantation is deliberately
##  absent: it changes no parameter here, which is why the graft re-accumulates.
## =============================================================================
hepatic_continuation <- function(mod, fres_grid = seq(0.30, 0.05, by = -0.01)) {
  do.call(rbind, lapply(fres_grid, function(f) {
    ss <- epp_steady(mod, list(FRES = f), days = 2000)
    data.frame(FRES = f, PPIXRBC = ss$PRBC, PLIV = ss$PLIV,
               CHOL = ss$CHOL, ALT = ss$ALT, TBIL = ss$TBIL)
  }))
}

scenario_hepatopathy_rescue <- function(mod) {
  # start ON the cholestatic branch: equilibrate a very severe genotype
  ss <- epp_steady(mod, list(FRES = 0.06), days = 3000)
  base <- init_from(mod, ss) %>% param(FRES = 0.08, EXPOSE = 0)
  arms <- list(
    "no treatment"                = list(),
    "cholestyramine 16 g/d"       = list(CHOLEFF = 0.85),
    "RBC transfusion"             = list(ERYSUP  = 0.40),
    "plasmapheresis"              = list(KPHER   = 0.50),
    "cholestyramine + transfusion"= list(CHOLEFF = 0.85, ERYSUP = 0.40),
    "all three"                   = list(CHOLEFF = 0.85, ERYSUP = 0.40,
                                         KPHER = 0.50),
    "bone-marrow transplant"      = list(FRES    = 1.00))
  do.call(rbind, lapply(names(arms), function(a) {
    o <- base %>% param(arms[[a]]) %>%
      mrgsim(end = 90 * 24, delta = 24, hmax = 0.5) %>% as.data.frame()
    tail(o, 1) %>% mutate(arm = a)
  }))
}

## =============================================================================
##  SCENARIO 9 — THE BEHAVIOURAL CONTROL LOOP
##  Identical photobiology, identical PPIX; only the retreat delay changes.
##  This is where a shielding drug shows a benefit that the registrational
##  "hours in sunlight" endpoint cannot see.
## =============================================================================
scenario_control_loop <- function(mod, delays_min = c(0.5, 1, 2, 5, 10, 20, 40)) {
  ss   <- epp_steady(mod, list(FRES = 0.15))
  base <- init_from(mod, ss) %>% param(FRES = 0.15, EXPOSE = 1)
  # afamelanotide arm sampled 30 days after an implant (near the melanin peak)
  ss30 <- init_from(mod, ss) %>% param(EXPOSE = 0) %>%
    mrgsim(events = ev_afamelanotide(1), end = 30 * 24, delta = 24) %>%
    as.data.frame() %>% tail(1)
  afa  <- init_from(mod, ss30) %>% param(FRES = 0.15, EXPOSE = 1)
  do.call(rbind, lapply(delays_min, function(d) {
    u <- base %>% param(TAUDEL = d/60) %>%
      mrgsim(end = 14 * 24, delta = 0.25, hmax = 0.02) %>% as.data.frame()
    a <- afa  %>% param(TAUDEL = d/60) %>%
      mrgsim(end = 14 * 24, delta = 0.25, hmax = 0.02) %>% as.data.frame()
    data.frame(delay_min = d,
               peak_pain_untreated = max(u$PAIN),
               rxn_days_untreated  = tail(u$RXNCUM, 1),
               peak_pain_afa       = max(a$PAIN),
               rxn_days_afa        = tail(a$RXNCUM, 1))
  }))
}

## =============================================================================
##  ENVIRONMENTAL MODIFIERS
##  Why window glass, cloud cover and SPF sunscreen all fail: the action
##  spectrum is 400-410 nm, which is visible violet, not ultraviolet.
## =============================================================================
environmental_table <- function(mod) {
  ss <- epp_steady(mod, list(FRES = 0.15))
  conds <- list(
    "open summer sun"              = list(),
    "behind ordinary window glass" = list(FGLASS = 0.90),
    "overcast sky"                 = list(FCLOUD = 0.50),
    "SPF 50 chemical sunscreen"    = list(),   # deliberately no violet effect
    "iron-oxide tinted sunscreen"  = list(FCLOTH = 0.35),
    "long sleeves + wide-brim hat" = list(FCLOTH = 0.15),
    "full opaque cover"            = list(FCLOTH = 0.05))
  do.call(rbind, lapply(names(conds), function(cn) {
    p <- conds[[cn]]
    Eeff <- prod(unlist(c(list(FGLASS = 1, FCLOUD = 1, FCLOTH = 1)[
      setdiff(c("FGLASS","FCLOUD","FCLOTH"), names(p))], p)))
    ttol <- ss$TTOLMIN / Eeff
    data.frame(condition = cn, transmitted = Eeff, tolerance_min = ttol,
               sun_h_day = 8 * (1 - exp(-3 * (ttol/60) / 8)))
  }))
}

## =============================================================================
##  RUN EVERYTHING
## =============================================================================
run_all <- function(mod) {
  sims <- bind_rows(
    scenario_untreated(mod),
    scenario_afamelanotide(mod),
    scenario_dersimelagon(mod, 300),
    scenario_bitopertin(mod, 60),
    scenario_combination(mod),
    scenario_xlp(mod)
  )

  cat("\n== SEASON SUMMARY (means over 180 days) ==\n")
  print(sims %>% group_by(scenario) %>%
          summarise(RBC_PPIX = last(PRBC), ZnPP = last(ZRBC),
                    metal_free_pct = last(MFFRAC), melanin = mean(MELIDX),
                    tolerance_min  = mean(TTOLMIN),
                    sun_h_per_day  = mean(SUNHRD),
                    season_sun_h   = mean(SUNHRD) * 180,
                    .groups = "drop"), n = 20)

  cat("\n== COMBINATION INDEX (is the interaction additive or multiplicative?) ==\n")
  print(combination_index(mod))

  cat("\n== IRON: THE SIGN FLIP ==\n")
  print(scenario_iron(mod))

  cat("\n== HEPATIC CONTINUATION (locate the folds) ==\n")
  print(hepatic_continuation(mod))

  cat("\n== HEPATOPATHY RESCUE ==\n")
  print(scenario_hepatopathy_rescue(mod) %>%
          select(arm, PLIV, CHOL, ALT, TBIL, PRBC))

  cat("\n== BEHAVIOURAL CONTROL LOOP ==\n")
  print(scenario_control_loop(mod))

  cat("\n== ENVIRONMENTAL MODIFIERS ==\n")
  print(environmental_table(mod))

  invisible(sims)
}

## =============================================================================
##  PLOTS
## =============================================================================
plot_biomarkers <- function(sims) {
  sims %>%
    select(time, scenario, `RBC PPIX (umol/L)` = PRBC,
           `Zn-PP (umol/L)` = ZRBC, `melanin OD` = MELIDX,
           `tolerance (min)` = TTOLMIN) %>%
    pivot_longer(-c(time, scenario)) %>%
    ggplot(aes(time / 24, value, colour = scenario)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~name, scales = "free_y") +
    labs(x = "day", y = NULL,
         title = "EPP QSP model — the two therapeutic axes are orthogonal",
         subtitle = paste("MC1R agonists move melanin and never PPIX;",
                          "bitopertin moves PPIX and never melanin")) +
    theme_bw() + theme(legend.position = "bottom")
}

plot_hysteresis <- function(mod) {
  ss <- epp_steady(mod, list(FRES = 0.15))
  init_from(mod, ss) %>% param(FRES = 0.15, EXPOSE = 0) %>%
    mrgsim(events = ev_afamelanotide(3, 60, 16),
           end = 240 * 24, delta = 2) %>%
    as.data.frame() %>%
    select(time, `afamelanotide plasma (ug/L)` = AC,
           `melanin OD` = MELIDX, `tolerance (min)` = TTOLMIN) %>%
    pivot_longer(-time) %>%
    ggplot(aes(time / 24, value)) +
    geom_line(linewidth = 0.7, colour = "#8040a0") +
    facet_wrap(~name, scales = "free_y", ncol = 1) +
    labs(x = "day", y = NULL,
         title = "A five-day drug with a sixty-day effect",
         subtitle = paste("afamelanotide plasma half-life ~30 min; melanin",
                          "half-life ~35 d — the PD outlasts the PK by",
                          "orders of magnitude")) +
    theme_bw()
}

## =============================================================================
##  ENTRY POINT
## =============================================================================
if (sys.nframe() == 0L) {
  sims <- run_all(mod)
  print(plot_biomarkers(sims))
}
