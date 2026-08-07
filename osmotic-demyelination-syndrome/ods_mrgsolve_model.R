## =====================================================================
##  ods_mrgsolve_model.R
##  Osmotic Demyelination Syndrome (ODS) — QSP / PK-PD model
##  40 ODE compartments · 20 therapeutic scenarios
##
##  삼투성 탈수초 증후군 — 정량적 시스템 약리학 모델
## ---------------------------------------------------------------------
##  STRUCTURAL THESIS
##  -----------------
##  Osmotic demyelination is not caused by hyponatraemia, and it is not
##  caused by sodium.  It is caused by an ASYMMETRY BETWEEN TWO TRANSPORT
##  TIME CONSTANTS inside the astrocyte:
##
##     EFFLUX of organic osmolytes   VRAC / LRRC8A, an ION CHANNEL   t½ ~ 15 h
##     INFLUX of organic osmolytes   SMIT1 / TauT / BGT1, CARRIERS
##                                   that TonEBP must first transcribe  t½ ~ 2-3 d
##
##  A channel opens in milliseconds.  A transporter has to be transcribed.
##  The brain can give solute away about six times faster than it can take
##  it back, and every clinical rule about correcting hyponatraemia is a
##  statement about that ratio.
##
##  The disease is carried by ONE derived quantity:
##
##      OMEGA(t) = ORG_set( Osm_eff(t) )  −  ORG(t)        [mOsm/kg brain water]
##
##  the ORGANIC OSMOLYTE DEFICIT — what the brain OUGHT to be holding at
##  the tonicity it now finds itself in, minus what it DOES hold.  OMEGA is
##  zero in the normal brain; it is zero in the chronically ADAPTED
##  hyponatraemic brain (which is why a sodium of 110 is not itself an
##  injury); and it becomes positive only when plasma tonicity moves faster
##  than transcription.
##
##  SIX RESULTS THE MODEL DERIVES RATHER THAN ASSUMES
##  --------------------------------------------------
##  (1) THE TWO GUIDELINE LIMITS ARE ONE NUMBER PLUS A TRANSPORTER.
##      A single universal astrocyte threshold, OMEGA* = 8 mOsm/kg, is used
##      for every patient.  With normal organic-transport capacity the
##      threshold is first crossed between +10 and +12 mmol/L/24 h.  With
##      the capacity of an alcoholic/malnourished patient (FOSM 0.55) it is
##      first crossed between +6 and +8.  Those are exactly the two numbers
##      in every guideline ("≤10-12 normal risk, ≤8 high risk"), and here
##      they are not two rules — they are one threshold read through two
##      transporters.  Risk factors in this model act ONLY on FOSM.
##
##  (2) THE KIDNEY, NOT THE PRESCRIPTION, SETS THE CORRECTION RATE.
##      The central experiment is one switch.  A hypovolaemic hyponatraemic
##      patient is given 0.9% saline; in arm A the AVP axis responds
##      physiologically, in arm B AVP is frozen at its presentation value.
##      Same patient, same fluid, no hypertonic saline in either arm.  In
##      arm A, volume repletion removes the non-osmotic AVP stimulus, urine
##      osmolality collapses from >1000 to ~50 mOsm/kg, electrolyte-free
##      water clearance opens, and the sodium climbs on its own; arm B is
##      quiet.  See ods_verification_output.txt §4 for the numbers.  The
##      overcorrection is made by the kidney being RELEASED, not by
##      anything that was infused.
##
##  (3) THEREFORE THE TREATMENT IS TO REMOVE THE KIDNEY FROM THE LOOP.
##      A proactive desmopressin clamp (2 µg IV q8h from the outset) plus
##      titrated 3% saline turns the correction rate back into something a
##      prescription can set.  DDAVP does not lower the sodium; it makes the
##      sodium obey the order.  The model then adds something that was not
##      designed in.  A clamp STORES WATER, and when the drug finally clears
##      the kidney excretes all of it at once: at a fluid intake of 1.5 L/d
##      the sodium rebounds 17 mmol/L in 24 h after the clamp is stopped and
##      the patient demyelinates — the clamp caused the injury it prevented.
##      At 0.5 L/d the rebound is 4.8 and the osmotic stress never approaches
##      threshold.  Tapering the DOSE from 2 µg to 0.25 µg changes nothing,
##      because 0.25 µg still gives ~10 pg/mL against a V2 EC50 of 1.6 —
##      the receptor is saturated until the drug is gone, whatever the
##      schedule.  The derived rule is therefore not "wean the desmopressin"
##      but "hold the free-water balance neutral while it is running": every
##      litre retained is a litre that leaves the moment the clamp comes off.
##
##  (4) POTASSIUM IS SODIUM.  Edelman's relation puts K_e in the numerator,
##      so 40 mmol of KCl into 42 L of body water is 1.11 × 40 / 42 =
##      1.06 mmol/L of sodium correction that appears on no fluid chart.
##      Repleting a hypokalaemic hyponatraemic patient at 120 mmol/day
##      spends a large part of the day's allowance before any saline is
##      hung.  Hypokalaemia is a risk factor by TWO independent routes in
##      this model — that one, and the loss of the Na+ gradient that the
##      Na+-coupled osmolyte carriers run on.
##
##  (5) ACUTE AND CHRONIC HYPONATRAEMIA ARE OPPOSITE DISEASES AT THE SAME
##      SODIUM.  The model was never told this.  An 8-hour water
##      intoxication to [Na] 110 leaves the organic osmolyte pool 93% full
##      and the brain swollen by 9% (herniation risk); a 21-day adaptation to
##      the SAME [Na] 110 leaves brain water normal and the pool 40% empty.
##      Correct both at +14 mmol/L/24 h and the chronic patient ends with a
##      peak deficit of 32 while the acute patient ends with 0.6.  The danger
##      in one is oedema and in the other is the treatment.
##
##  (6) THE MRI IS LATE ON PURPOSE.  Astrocytes die within a day or two,
##      myelin is lost over the next few days, the patient deteriorates
##      after that, and the radiological lesion appears later still.  The
##      model reproduces the biphasic course — sodium normal, patient
##      better, deficit still zero — and predicts a negative scan at the
##      moment of clinical deterioration.
##
##  WHAT WAS FITTED AND WHAT WAS NOT
##  ---------------------------------
##  Fitted to NORMAL physiology (not to ODS):  the water/solute steady
##  state (WIN solved so that the healthy model is an exact steady state,
##  max|dy/dt| = 2.5e-11), Edelman's published regression, urine
##  concentrating range, plasma urea, AVP osmotic threshold and gain.
##  Fitted to ADAPTATION data:  the osmoresponsiveness coefficients β_i,
##  chosen so that chronic adaptation to [Na] 110 costs ~40% of the total
##  organic pool and ~66% of myo-inositol (Verbalis & Gullans 1991), and
##  the influx/efflux time constants, chosen so that myo-inositol takes
##  ~5 days to come back (Verbalis & Gullans 1993).
##  Fitted to the INJURY dose-response:  four constants (KINJ, KAST, KOLI,
##  KDEM) chosen so that the peak clinical deficit is graded across
##  correction rates.  The THRESHOLD itself (OMEGA* = 8) is calibrated so
##  that the normal-risk limit lands between 10 and 12 — one number, and
##  the high-risk limit of 6-8 is then a PREDICTION of FOSM = 0.55, not a
##  second fit.
##  NOT fitted, i.e. predictions:  the acute/chronic asymmetry, the
##  autonomous overcorrection of §4, the size of the DDAVP effect, the
##  relowering deadline, the Furst urine/serum ratio rule, the failure of
##  0.9% saline when urine osmolality exceeds 308, and the MRI lag.
##
##  LIMITATIONS — stated, not buried
##  ---------------------------------
##  a.  R was not available in the build container, so this file is
##      EQUATION-VERIFIED but NOT COMPILE-VERIFIED.  Every ODE here was
##      independently re-implemented in Python/scipy
##      (ods_verify_python.py) and integrated; that file is the source of
##      every number quoted in the documentation.
##  b.  OMEGA has never been measured in a living human brain.  What has
##      been measured is ¹H-MRS myo-inositol (Videen 1995, Restuccia 2004),
##      which is one term of it.  The falsifiable prediction is that the
##      MRS myo-inositol deficit at 24-48 h after correction, not the
##      24-hour sodium rise, is what separates the patients who
##      demyelinate from those who do not.
##  c.  Urea is given no osmotic protective role.  Under a pure
##      effective-tonicity accounting it cannot have one — it crosses the
##      barrier and cancels from both sides.  The animal protection data
##      are real, so urea here acts on the blood-brain barrier and
##      microglia, where those studies actually measured an effect.  If a
##      study ever shows urea blunting the osmotic insult itself, this
##      model is wrong in a specific and checkable way.
##  d.  The clinical deficit scale (0-100) is an ordinal construct mapped
##      onto lesion burden.  It is not an mRS and should not be read as one.
##  e.  Serum potassium is an algebraic function of the total-body deficit,
##      not a transcellular-shift model; acid-base, insulin and beta-agonist
##      effects on [K] are absent.
##  f.  Cerebral blood flow, ICP and the glymphatic system are not modelled;
##      the herniation index is a monotone proxy for brain swelling only.
##
##  Units: TIME = DAYS.  Brain solutes = mOsm (or mmol) per kg of BASELINE
##  brain water.  Body solutes = mmol.  Volumes = L.  Requires mrgsolve.
## =====================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

code <- '
$PROB
# Osmotic Demyelination Syndrome — QSP model (40 ODEs)

$PARAM @annotated
// ---- body fluid / Edelman -------------------------------------------
EDA      : 1.11    : Edelman slope
EDB      : 25.6    : Edelman intercept (mmol/L)
FECF     : 0.65    : osmotically active fraction of exchangeable Na
OSMX     : 5.0     : non-Na effective osmoles (mOsm/kg)

// ---- renal solute handling ------------------------------------------
NAINT    : 120     : dietary Na intake (mmol/d)
KINT     : 60      : dietary K intake (mmol/d)
NAESET   : 3000    : exchangeable Na set point (mmol)
KESET    : 3266    : exchangeable K set point (mmol)
KNAEX    : 1.50    : Na excretion gain (/d)
KKEX     : 1.00    : K excretion gain (/d)
PUREA    : 400     : urea generation (mmol/d)
PUREA2   : 400     : urea generation after TSOLUP (mmol/d)
TSOLUP   : 1e9     : time of solute-intake step (d)
CLUREA   : 55      : renal urea clearance (L/d)
OTHOSM   : 60      : other urinary osmoles (mOsm/d)
UOSMMIN  : 50      : maximally dilute urine (mOsm/kg)
UOSMMAX  : 1150    : maximally concentrated urine (mOsm/kg)
SOL0     : 900     : solute load at which the medulla starts to wash out
SOLREF   : 900     : washout scale (mOsm/d)

// ---- water ------------------------------------------------------------
WIN      : 1.638295 : oral water intake (L/d) - solved for exact steady state
WMET     : 0.30    : metabolic water (L/d)
WINSENS  : 0.90    : insensible loss (L/d)
OSMTHIRST: 292     : thirst threshold (mOsm/kg)
KTHIRST  : 0.25    : L/d per mOsm/kg above thirst threshold

// ---- AVP axis ----------------------------------------------------------
OSMTHR   : 280     : osmotic threshold for AVP (mOsm/kg)
GAINOSM  : 0.45    : pg/mL per mOsm/kg
AVPVMAX  : 12      : maximum non-osmotic AVP (pg/mL)
VOLEC50  : 0.08    : ECF deficit for half-maximal AVP (fraction)
VOLHILL  : 3.0     : Hill coefficient, baroreceptor arm
AVPSIADH : 0       : autonomous AVP secretion (pg/mL)
TAVPOFF  : 1e9     : time the AVP stimulus is removed (d)
AVPFREEZE: 0       : 1 = clamp AVP (counterfactual)
TFREEZE  : 1e9     : time the clamp starts (d)
AVPFRZV  : 0       : clamped AVP value (pg/mL)
TAUAVP   : 0.015   : AVP turnover (d)
EC50AVP  : 1.20    : V2 EC50 for AVP (pg/mL)
EC50DDA  : 1.60    : V2 EC50 for desmopressin (pg/mL)
KITLV    : 30      : apparent V2 Ki for tolvaptan, total drug (ng/mL)
TAUAQP   : 0.060   : AQP2 trafficking (d)
THIAZ    : 0       : thiazide flag

// ---- brain baseline composition (mOsm per kg baseline brain water) -----
IMP      : 10.0    : impermeant solute
INS0     : 7.0     : myo-inositol
TAU0     : 2.0     : taurine
GLX0     : 15.0    : glutamate + glutamine
CRE0     : 10.0    : creatine + phosphocreatine
GPC0     : 2.5     : GPC + PC
OTH0     : 11.5    : NAA + glycine + betaine + others
BINS     : 3.2     : osmoresponsiveness, myo-inositol
BTAU     : 3.2     : osmoresponsiveness, taurine
BGLX     : 2.6     : osmoresponsiveness, Glx
BCRE     : 1.0     : osmoresponsiveness, creatine
BGPC     : 2.4     : osmoresponsiveness, GPC
BOTH     : 0.8     : osmoresponsiveness, others
ORGFLOOR : 0.12    : floor as fraction of the baseline pool

// ---- brain transport kinetics -------------------------------------------
TAUEFF   : 0.90    : VRAC efflux time constant (d), t1/2 15 h
TAUINF   : 1.25    : carrier influx time constant (d)
TAUTON   : 0.25    : TonEBP activation (d)
TAUSMIT  : 1.00    : carrier protein turnover (d)
ETON     : 0.9     : TonEBP gain
TONMAX   : 2.6     : TonEBP ceiling
FOSM     : 1.0     : organic transport capacity (THE risk factor)
TAUELEC  : 0.15    : RVI/RVD time constant (d)
JRVIMAX  : 60      : maximum inorganic accumulation (mOsm/kg/d)
JRVDMAX  : 150     : maximum inorganic loss (mOsm/kg/d)
BELCEIL  : 1.10    : inorganic ceiling (x baseline)
TAUBUR   : 0.10    : brain-plasma urea equilibration (d)

// ---- injury --------------------------------------------------------------
WSHR     : 1.0     : weight of shrinkage in the stress term
OMSTAR   : 8.0     : injury threshold (mOsm/kg) - ONE number for everybody
RISK     : 1.0     : threshold multiplier (kept at 1 in every scenario)
KINJ     : 0.300   : injury rate constant
HINJ     : 1.00    : injury exponent above threshold
FNUT     : 1.0     : astrocyte energy reserve
TAUATP   : 0.10    : ATP recovery (d)
KATPU    : 0.030   : ATP consumption per unit stress
KAST     : 0.30    : astrocyte death gain
ASTTHR   : 0.030   : astrocyte loss tolerated by the gap-junction syncytium
KASTREP  : 0.250   : astrocyte repopulation (/d)
KBBB     : 0.60    : BBB opening gain
TAUBBB   : 2.50    : BBB resealing (d)
KMG      : 0.50    : microglial activation gain
TAUMG    : 4.00    : microglial decay (d)
KCYT     : 0.60    : cytokine production gain
TAUCYT   : 1.50    : cytokine decay (d)
KIGG     : 0.50    : IgG/complement entry gain
TAUIGG   : 3.00    : IgG clearance (d)
KOLI     : 0.20    : oligodendrocyte death gain
KHIT     : 1.20    : saturation of the oligodendrocyte hit
KSEV     : 0.42    : lesion burden for half-maximal deficit
WOA      : 1.00    : weight, astrocyte loss
WOC      : 0.45    : weight, cytokines
WOI      : 0.35    : weight, IgG/complement
KOPCP    : 0.070   : OPC proliferation (/d)
KOPCD    : 0.055   : OPC differentiation (/d)
KDEM     : 0.40    : acute demyelination gain
KMYE     : 0.060   : myelin-oligodendrocyte coupling (/d)
WPONS    : 1.00    : topographic weight, central pons
WEXP     : 0.55    : topographic weight, extrapontine
TAULES   : 45      : lesion resolution (d)
TAUMRI   : 3.5     : MRI signal lag (d)
TAUDEF   : 1.50    : clinical deficit lag (d)
DEFMAX   : 100     : deficit scale

// ---- drug PK ---------------------------------------------------------------
VDDA     : 25      : desmopressin volume (L)
KDDA     : 5.55    : desmopressin elimination (/d), t1/2 3 h
KADDA    : 36      : desmopressin SC absorption (/d)
VTLV     : 210     : tolvaptan V/F (L)
KTLV     : 2.08    : tolvaptan elimination (/d), t1/2 8 h
KATLV    : 12      : tolvaptan absorption (/d)
FTLV     : 0.50    : tolvaptan bioavailability
KAUREA   : 48      : oral urea absorption (/d)
UREADOSE : 0       : oral urea (g/d)
VDEX     : 60      : dexamethasone volume (L)
KDEX     : 4.16    : dexamethasone elimination (/d)
DEXON    : 0       : dexamethasone 16 mg/d flag
VMIN     : 80      : minocycline volume (L)
KMIN     : 1.04    : minocycline elimination (/d)
MINOON   : 0       : minocycline 200 mg/d flag
EMAXUREA : 0.55    : urea Emax on BBB/microglia
EC50UREA : 12      : urea EC50 (mmol/L above baseline)
EMAXDEX  : 0.45    : dexamethasone Emax on BBB
EC50DEX  : 15      : dexamethasone EC50 (ng/mL)
EMAXMIN  : 0.50    : minocycline Emax on microglia
EC50MIN  : 2.0     : minocycline EC50 (ug/mL)

// ---- fluids, losses and the correction controller ---------------------------
R3PCT    : 0       : fixed 3% NaCl rate (L/d)
R09      : 0       : fixed 0.9% NaCl rate (L/d)
RD5W     : 0       : fixed 5% dextrose rate (L/d)
RKCL     : 0       : KCl infusion (mmol/d)
NA3      : 513     : [Na] of 3% saline (mmol/L)
NA09     : 154     : [Na] of 0.9% saline (mmol/L)
NALOSS   : 0       : extrarenal Na loss (mmol/d)
KLOSS    : 0       : extrarenal K loss (mmol/d)
WLOSS    : 0       : extrarenal water loss (L/d)
TLOSSEND : 0       : losses stop here (d)
ACUTE    : 0       : acute water-loading flag
WLOAD    : 0       : free-water load rate (L/d)
TWLEND   : 0       : end of the water load (d)
CTRLON   : 0       : closed-loop 3% NaCl titration flag
TCORR    : 0       : time correction starts (d)
RATETGT  : 6.0     : prescribed correction rate (mmol/L/d)
NASTART  : 110     : [Na] at the start of correction
NACAP    : 135     : [Na] at which the prescription ends
KP3      : 2.0     : controller gain (L/d per mmol/L)
R3MAX    : 4.0     : 3% NaCl clamp (L/d)
RESCUE   : 0       : relowering flag
TRESCUE  : 1e9     : time relowering starts (d)
NARES    : 118     : relowering target [Na]
KPD5W    : 0.8     : relowering controller gain (L/d per mmol/L)
DURRES   : 1.0     : relowering treatment window (d)
D5WMAX   : 5.0     : 5% dextrose clamp (L/d), ~3 mL/kg/h

$CMT @annotated
TBW    : total body water (L)
NAE    : exchangeable sodium (mmol)
KE     : exchangeable potassium (mmol)
UREAB  : total body urea (mmol)
AVP    : plasma vasopressin (pg/mL)
AQP2   : apical AQP2 water permeability (0-1)
BELEC  : brain inorganic solute (mOsm/kg brain water)
INS    : brain myo-inositol
TAUR   : brain taurine
GLX    : brain glutamate + glutamine
CRE    : brain creatine + phosphocreatine
GPC    : brain GPC + PC
OTH    : brain NAA + others
BURE   : brain urea
TONEBP : TonEBP/NFAT5 nuclear activity
SMIT   : SMIT1/TauT/BGT1 carrier protein
ATP    : astrocyte energy charge
AST    : astrocyte viability
BBBP   : blood-brain barrier permeability (1 = normal)
MG     : microglial activation
CYT    : TNF/IL-1b composite
IGG    : brain IgG / complement
OLI    : oligodendrocyte viability
OPC    : oligodendrocyte precursor pool
MYE    : myelin content
LESP   : pontine lesion burden
LESE   : extrapontine lesion burden
MRI    : radiological signal
DEF    : clinical deficit (0-100)
DDAD   : desmopressin SC depot (ug)
DDAC   : desmopressin central (ug)
TLVD   : tolvaptan gut (mg)
TLVC   : tolvaptan central (mg)
UREAG  : oral urea in gut (mmol)
DEXC   : dexamethasone (mg)
MINC   : minocycline (mg)
CUMI   : cumulative injury integral
CUMNA  : cumulative sodium infused (mmol)
CUMV   : cumulative volume infused (L)
CUMEFW : cumulative electrolyte-free water clearance (L)

$GLOBAL
// --------------------------------------------------------------------
//  SMOOTH SWITCHES.  Every hinge in this model is C1 and EXACTLY zero
//  below the knee.  Two numerical lessons are baked in here and both were
//  found the hard way in the Python twin:
//
//   (i)  The popular smooth hinge 0.5*(x + sqrt(x*x + e*e)) is NOT zero
//        below the knee; it leaves a floor of e^2/(4|x|).  Harmless in a
//        rate equation, fatal in an integrated one — fed through the
//        astrocyte -> microglia -> oligodendrocyte chain that floor
//        demyelinated a COMPLETELY HEALTHY brain to 75% of normal over 60
//        simulated days.  SPOS below is exactly zero for x <= 0.
//
//   (ii) A switch on a MEASURED state ("stop the infusion when [Na]
//        reaches the cap") makes the controller hunt at the cap; LSODA
//        spent 7.7e5 function evaluations there and never finished.  The
//        controller below stops at a TIME computed from the prescription,
//        which is deterministic and is what a prescription actually is.
// --------------------------------------------------------------------
#define SPOS(x, e)  ((x) <= 0.0 ? 0.0 : ((x) >= (e) ? (x) - 0.5*(e) : (x)*(x)/(2.0*(e))))
#define SCLAMP(x, lo, hi, e)  ((lo) + SPOS((x)-(lo), e) - SPOS((x)-(hi), e))
#define GATE(x, w)  (0.5*(1.0 + tanh((x)/(w))))

// Variables that $TABLE captures must live at FILE scope: a `double`
// declared inside $ODE is destroyed before $TABLE runs, so capturing it
// would silently record whatever the compiler left behind.
double UOSM, VU, UNAK, EFWC, ORG, ORGSET, OMEGA, STRESS, BWREL, SHRINK, SWELLB, INJ, DDACP, TLVCP;

$MAIN
// initial conditions: the exact healthy steady state (max|dy/dt| = 2.5e-11)
TBW_0    = 42.0;
NAE_0    = NAESET;
KE_0     = KESET;
UREAB_0  = PUREA/CLUREA*42.0;
AVP_0    = 2.25;
AQP2_0   = 0.652;
BELEC_0  = 285.0 - IMP - (INS0+TAU0+GLX0+CRE0+GPC0+OTH0);
INS_0    = INS0;
TAUR_0   = TAU0;
GLX_0    = GLX0;
CRE_0    = CRE0;
GPC_0    = GPC0;
OTH_0    = OTH0;
BURE_0   = PUREA/CLUREA;
TONEBP_0 = 1.0;
SMIT_0   = 1.0;
ATP_0    = 1.0;
AST_0    = 1.0;
BBBP_0   = 1.0;
OLI_0    = 1.0;
MYE_0    = 1.0;

$ODE
// =====================================================================
//  ALGEBRAIC LAYER  (identical to algebra() in ods_verify_python.py)
// =====================================================================
double TBWx = TBW > 5.0 ? TBW : 5.0;

// ---- Edelman ---------------------------------------------------------
double SNA    = EDA*(NAE + KE)/TBWx - EDB;
double OSMEFF = 2.0*SNA + OSMX;
double UREAP  = UREAB/TBWx;
double SK     = 4.0 + (KE - KESET)/300.0;   if(SK < 1.5) SK = 1.5;

// ---- ECF volume and the non-osmotic AVP drive -------------------------
double SNAx  = SNA > 60.0 ? SNA : 60.0;
double ECFV  = FECF*NAE/SNAx;
double ECFV0 = FECF*NAESET/140.0;
double VD    = SPOS((ECFV0 - ECFV)/ECFV0, 2e-3);
double AVPOSM = GAINOSM*SPOS(OSMEFF - OSMTHR, 0.30);
double VDH    = pow(VD, VOLHILL);
double AVPVOL = AVPVMAX*VDH/(pow(VOLEC50, VOLHILL) + VDH);
double AVPTGT = AVPOSM + AVPVOL + AVPSIADH*GATE(TAVPOFF - SOLVERTIME, 0.01);
if(AVPFREEZE > 0.5){
  double wfz = GATE(SOLVERTIME - TFREEZE, 0.01);
  AVPTGT = (1.0 - wfz)*AVPTGT + wfz*AVPFRZV;
}

// ---- V2 receptor: AVP + desmopressin against tolvaptan -----------------
DDACP = (DDAC > 0 ? DDAC : 0.0)/VDDA*1e6;     // ug   -> pg/mL
TLVCP = (TLVC > 0 ? TLVC : 0.0)/VTLV*1e6;     // mg   -> ng/mL
double AGO   = (AVP > 0 ? AVP : 0.0)/EC50AVP + DDACP/EC50DDA;
double OCC   = AGO/(1.0 + AGO + TLVCP/KITLV);

// ---- urine -------------------------------------------------------------
double UOMIN = UOSMMIN + 100.0*THIAZ;
double ENA   = 0.05*NAINT + SPOS(KNAEX*(NAE - NAESET) + 0.95*NAINT, 0.5);
double EK    = 0.05*KINT  + SPOS(KKEX *(KE  - KESET ) + 0.95*KINT , 0.5);
double EUREA = CLUREA*UREAP;
double WINT  = WIN + KTHIRST*SPOS(OSMEFF - OSMTHIRST, 0.3);
double PU    = PUREA + (PUREA2 - PUREA)*GATE(SOLVERTIME - TSOLUP, 0.02);
double LS    = GATE(TLOSSEND - SOLVERTIME, 0.004);
double SOLEXC = 2.0*(ENA + EK) + EUREA + OTHOSM;
// flow-dependent medullary washout: the countercurrent multiplier cannot
// hold its gradient against a large solute load.  Without this term the
// model produces urine [Na+K] above 400 mmol/L, which does not exist.
double UOMAX = UOSMMAX/(1.0 + SPOS(SOLEXC - SOL0, 1.0)/SOLREF);
UOSM  = UOMIN + (UOMAX - UOMIN)*AQP2;
VU    = SOLEXC/(UOSM + 30.0);
UNAK  = (ENA + EK)/(VU > 1e-6 ? VU : 1e-6);
EFWC  = VU*(1.0 - UNAK/SNAx);

// ---- infusions and the controller --------------------------------------
double R3 = R3PCT, R09v = R09, RD5 = RD5W;
if(ACUTE > 0.5)  RD5 += WLOAD*GATE(TWLEND - SOLVERTIME, 0.004);
if(CTRLON > 0.5){
  double rt    = RATETGT > 0.1 ? RATETGT : 0.1;
  double tstop = TCORR + (NACAP - NASTART)/rt;
  double ramp  = NASTART + RATETGT*(SOLVERTIME - TCORR);
  R3 += SCLAMP(KP3*(ramp - SNA), 0.0, R3MAX, 0.02)
        * GATE(SOLVERTIME - TCORR, 0.01) * GATE(tstop - SOLVERTIME, 0.01);
}
if(RESCUE > 0.5)
  RD5 += SCLAMP(KPD5W*(SNA - NARES), 0.0, D5WMAX, 0.02)
         * GATE(SOLVERTIME - TRESCUE, 0.01);
double VINF  = R3 + R09v + RD5;
double NAINF = R3*NA3 + R09v*NA09;

// ---- brain set points ---------------------------------------------------
double FR = (OSMEFF - 285.0)/285.0;
#define SETP(base, beta) ((base)*ORGFLOOR + SPOS((base)*(1.0 + (beta)*FR) - (base)*ORGFLOOR, 0.02))
double INSs = SETP(INS0, BINS);
double TAUs = SETP(TAU0, BTAU);
double GLXs = SETP(GLX0, BGLX);
double CREs = SETP(CRE0, BCRE);
double GPCs = SETP(GPC0, BGPC);
double OTHs = SETP(OTH0, BOTH);
ORGSET = INSs + TAUs + GLXs + CREs + GPCs + OTHs;
double ORG0   = INS0 + TAU0 + GLX0 + CRE0 + GPC0 + OTH0;
ORG    = (INS>0?INS:0) + (TAUR>0?TAUR:0) + (GLX>0?GLX:0)
              + (CRE>0?CRE:0) + (GPC>0?GPC:0)   + (OTH>0?OTH:0);

double BELVOL = OSMEFF - IMP - ORG;                    // what RVI is asked for
double BELCAP = (285.0 - IMP - ORG0)*BELCEIL;          // what the pump can hold

// ---- brain water: osmotic equilibrium is instantaneous ------------------
double SOLUTE = (BELEC > 0 ? BELEC : 0.0) + ORG + IMP;
BWREL  = SOLUTE/OSMEFF;
SHRINK = SPOS(1.0 - BWREL, 5e-4);
SWELLB = SPOS(BWREL - 1.0, 5e-4);

// ---- THE DRIVER ---------------------------------------------------------
OMEGA  = SPOS(ORGSET - ORG, 0.03);
STRESS = OMEGA + WSHR*SHRINK*OSMEFF;
double EXC    = SPOS(STRESS - OMSTAR*RISK, 0.03);
double ATPx   = ATP > 0.05 ? ATP : 0.05;
INJ    = KINJ*pow(EXC, HINJ)/ATPx;

// =====================================================================
//  DERIVATIVES
// =====================================================================
dxdt_TBW   = WINT + WMET - WINSENS - VU + VINF - WLOSS*LS;
dxdt_NAE   = NAINT + NAINF - ENA - NALOSS*LS;
dxdt_KE    = KINT  + RKCL  - EK   - KLOSS*LS;
dxdt_UREAB = PU + KAUREA*UREAG - EUREA;

dxdt_AVP   = (AVPTGT - AVP)/TAUAVP;
dxdt_AQP2  = (OCC - AQP2)/TAUAQP;

// ---- brain inorganic: RVI / RVD (channels, fast both ways) --------------
double BTGT = BELCAP - SPOS(BELCAP - BELVOL, 0.05);      // = min(BELVOL, BELCAP)
dxdt_BELEC = SCLAMP((BTGT - BELEC)/TAUELEC, -JRVDMAX, JRVIMAX, 0.5);

// ---- brain organic osmolytes: THE ASYMMETRY -----------------------------
//   efflux = VRAC/LRRC8A, a channel      -> fast, energy-independent
//   influx = SMIT1/TauT/BGT1, carriers   -> slow, transcription-gated and
//                                           dependent on the Na+ gradient
double KFAC = SK/4.0;  if(KFAC > 1.10) KFAC = 1.10;  if(KFAC < 0.40) KFAC = 0.40;
double GINF = FOSM*SMIT*KFAC*(ATP > 0.2 ? ATP : 0.2);
#define ORGDX(state, setp) { double df = (setp) - (state); \
    double up = SPOS(df, 0.03); double dn = up - df;                       \
    dxdt_##state = up/TAUINF*GINF - dn/TAUEFF; }
ORGDX(INS,  INSs)
ORGDX(TAUR, TAUs)
ORGDX(GLX,  GLXs)
ORGDX(CRE,  CREs)
ORGDX(GPC,  GPCs)
ORGDX(OTH,  OTHs)
dxdt_BURE = (UREAP - BURE)/TAUBUR;

// ---- TonEBP -> carrier transcription -------------------------------------
double TONT = 1.0 + (TONMAX - 1.0)
            - SPOS((TONMAX - 1.0) - ETON*OMEGA/ORG0*10.0, 0.02);
dxdt_TONEBP = (TONT - TONEBP)/TAUTON;
dxdt_SMIT   = (TONEBP - SMIT)/TAUSMIT;

// ---- astrocyte energetics and viability -----------------------------------
double SKx = SK > 1.5 ? SK : 1.5;
dxdt_ATP = (FNUT - ATP)/TAUATP - KATPU*STRESS*ATP*(4.0/SKx);
dxdt_AST = -KAST*INJ*AST + KASTREP*(1.0 - AST);
// The astrocyte syncytium is gap-junction coupled (Cx43/Cx30): scattered
// loss is covered by neighbours, and oligodendrocyte trophic support fails
// only once loss exceeds ASTTHR.  Without this threshold a very small but
// very PROLONGED astrocyte deficit integrates into total demyelination.
double ASTLOSS = SPOS((1.0 - AST) - ASTTHR, 0.01);

// ---- BBB, microglia, cytokines, humoral entry ------------------------------
double UREAX = SPOS(BURE - 7.3, 0.05);
double PU_   = EMAXUREA*UREAX/(EC50UREA + UREAX);
double DEXCP = (DEXC > 0 ? DEXC : 0.0)/VDEX*1e6;
double PD_   = EMAXDEX*DEXCP/(EC50DEX + DEXCP);
double MINCP = (MINC > 0 ? MINC : 0.0)/VMIN*1e3;
double PM_   = EMAXMIN*MINCP/(EC50MIN + MINCP);

dxdt_BBBP = KBBB*ASTLOSS*(1.0 - PU_)*(1.0 - PD_) - (BBBP - 1.0)/TAUBBB;
dxdt_MG   = KMG*(ASTLOSS + 0.5*SPOS(BBBP - 1.0, 1e-3))*(1.0 - PM_)*(1.0 - PU_)
            - MG/TAUMG;
dxdt_CYT  = KCYT*MG - CYT/TAUCYT;
dxdt_IGG  = KIGG*SPOS(BBBP - 1.0, 1e-3) - IGG/TAUIGG;

// ---- oligodendrocyte and myelin --------------------------------------------
double HIT0 = WOA*ASTLOSS + WOC*CYT + WOI*IGG;
double HIT  = HIT0/(KHIT + HIT0);            // the pons cannot demyelinate twice
dxdt_OLI = -KOLI*HIT*OLI + KOPCD*OPC*(1.0 - OLI);
dxdt_OPC = KOPCP*(1.0 - OLI)*(1.0 - OPC) - KOPCD*OPC*(1.0 - OLI);
// myelin TRACKS the oligodendrocyte population (a 10% cell deficit is a ~10%
// myelin deficit, not a collapse) plus an acute stripping term.
dxdt_MYE = KMYE*(OLI - MYE) - KDEM*HIT*MYE;

// ---- topography, imaging, clinical deficit ----------------------------------
double DM = SPOS(1.0 - MYE, 1e-4);
dxdt_LESP = (WPONS*DM - LESP)/TAULES*12.0;
dxdt_LESE = (WEXP *DM - LESE)/TAULES*12.0;
double LMAX = LESP > LESE ? LESP : LESE;
dxdt_MRI  = (LMAX - MRI)/TAUMRI;
double LB  = LESP + 0.5*LESE;
double SEV = pow(LB, 1.4)/(pow(KSEV, 1.4) + pow(LB, 1.4));
dxdt_DEF  = (DEFMAX*SEV - DEF)/TAUDEF;

// ---- drug PK -----------------------------------------------------------------
dxdt_DDAD  = -KADDA*DDAD;
dxdt_DDAC  = KADDA*DDAD - KDDA*DDAC;
dxdt_TLVD  = -KATLV*TLVD;
dxdt_TLVC  = FTLV*KATLV*TLVD - KTLV*TLVC;
dxdt_UREAG = UREADOSE/0.060 - KAUREA*UREAG;
dxdt_DEXC  = DEXON*16.0  - KDEX*DEXC;
dxdt_MINC  = MINOON*200.0 - KMIN*MINC;

// ---- bookkeeping ---------------------------------------------------------------
dxdt_CUMI   = INJ;
dxdt_CUMNA  = NAINF;
dxdt_CUMV   = VINF;
dxdt_CUMEFW = EFWC;

$TABLE
capture SODIUM   = EDA*(NAE + KE)/(TBW > 5 ? TBW : 5) - EDB;
capture TONICITY = 2.0*SODIUM + OSMX;
capture POTASSIUM= (4.0 + (KE - KESET)/300.0 < 1.5) ? 1.5 : 4.0 + (KE - KESET)/300.0;
capture BUNmgdL  = UREAB/(TBW > 5 ? TBW : 5)*2.80;
capture AVPpg    = AVP;
capture UOSMOL   = UOSM;
capture UVOL     = VU;
capture UNAKc    = UNAK;
capture EFWCL    = EFWC;
capture FURST    = UNAK/(SODIUM > 60 ? SODIUM : 60);
capture ORGTOT   = ORG;
capture ORGSETC  = ORGSET;
capture OMEGAc   = OMEGA;
capture STRESSc  = STRESS;
capture THRESH   = OMSTAR*RISK;
capture BRAINH2O = 80.0*BWREL;
capture SHRINKc  = SHRINK*100.0;
capture SWELLc   = SWELLB*100.0;
capture HERNIDX  = 100.0/(1.0 + exp(-(SWELLB - 0.070)/0.012));
capture INJRATE  = INJ;
capture MYOINOS  = INS;
capture DDAVPpg  = DDACP;
capture TLVngmL  = TLVCP;
capture DEFICIT  = DEF;
capture MRIPOS   = MRI;
'

mod <- mcode("ods", code)

## =====================================================================
##  HELPERS
## =====================================================================

## Build a chronically ADAPTED hyponatraemic patient BY SIMULATION.
## The adapted brain composition is an OUTPUT of the model, never an input:
## the osmolyte pool empties because the model empties it.
make_chronic <- function(mod, pheno = c("siadh","hypovol","thiazide","potomania","adrenal"),
                         win, days = 21) {
  pheno <- match.arg(pheno)
  pset <- switch(pheno,
    siadh     = list(AVPSIADH = 6.0, WIN = win),
    hypovol   = list(NALOSS = 250, KLOSS = 60, WLOSS = 1.6, TLOSSEND = days, WIN = win),
    thiazide  = list(THIAZ = 1, NAINT = 40, KINT = 20, AVPSIADH = 2.0, WIN = win),
    potomania = list(NAINT = 15, KINT = 15, PUREA = 110, UOSMMAX = 600, WIN = win),
    adrenal   = list(AVPSIADH = 4.0, NAINT = 30, KINT = 90, WIN = win))
  out <- mod %>% param(pset) %>% mrgsim(end = days, delta = 0.25)
  list(param = pset, state = as.numeric(tail(out, 1)[, names(init(mod))]), out = out)
}

## Solve for the water intake that gives a target sodium at `days`.
solve_win <- function(mod, pheno, target = 110, days = 21, lo = 1.0, hi = 14.0) {
  f <- function(w) {
    r <- make_chronic(mod, pheno, win = w, days = days)
    tail(r$out$SODIUM, 1) - target
  }
  if (f(hi) > 0) return(hi)
  uniroot(f, c(lo, hi), tol = 1e-4)$root
}

## Desmopressin 2 ug IV q8h
ddavp_q8 <- function(from, to, amt = 2, ID = 1)
  ev(ID = ID, time = from, cmt = "DDAC", amt = amt, ii = 8/24,
     addl = max(0, floor((to - from)/(8/24)) - 1))

## Tolvaptan 15 mg once daily
tolvaptan <- function(days = 7, amt = 15, ID = 1)
  ev(ID = ID, time = 0, cmt = "TLVD", amt = amt, ii = 1, addl = days - 1)

## Largest rise in [Na] over any 24 h window
na_rise_24 <- function(time, na) {
  best <- 0
  for (i in seq_along(time)) {
    j <- which(time >= time[i] + 1)[1]
    if (is.na(j)) break
    best <- max(best, na[j] - na[i])
  }
  best
}

## =====================================================================
##  20 SCENARIOS
##  Every scenario starts from a patient the model built for itself.
##  Reference values for each are in ods_verification_output.txt.
## =====================================================================
run_scenarios <- function(mod) {

  ## --- build the source patients -------------------------------------
  w_siadh <- solve_win(mod, "siadh",     110, 21)
  w_hypo  <- solve_win(mod, "hypovol",   110,  7)
  w_pot   <- solve_win(mod, "potomania", 110, 21)
  w_adr   <- solve_win(mod, "adrenal",   110, 21)

  SI <- make_chronic(mod, "siadh",     w_siadh, 21)
  HY <- make_chronic(mod, "hypovol",   w_hypo,   7)
  PO <- make_chronic(mod, "potomania", w_pot,   21)
  AD <- make_chronic(mod, "adrenal",   w_adr,   21)

  ## An ACUTE patient: the same sodium reached in 8 hours, so the
  ## osmolytes have not left yet.  This is a different disease.
  w_ac <- uniroot(function(x) {
      o <- mod %>% param(AVPSIADH = 6, ACUTE = 1, WLOAD = x, TWLEND = 8/24) %>%
             mrgsim(end = 8/24, delta = 1/96)
      tail(o$SODIUM, 1) - 110
    }, c(1, 60), tol = 1e-4)$root
  AC <- mod %>% param(AVPSIADH = 6, ACUTE = 1, WLOAD = w_ac, TWLEND = 8/24) %>%
          mrgsim(end = 8/24, delta = 1/96)
  ac_state <- as.numeric(tail(AC, 1)[, names(init(mod))])

  start_at <- function(mod, state) do.call(init, c(list(mod), as.list(setNames(state, names(init(mod))))))
  corr <- function(rate, na0 = 110, cap = 140)
    list(CTRLON = 1, RATETGT = rate, NASTART = na0, NACAP = cap, WIN = 1.0)
  stop_losses <- list(NALOSS = 0, KLOSS = 0, WLOSS = 0, TLOSSEND = 0)

  S <- list()

  ## S01  normal control ------------------------------------------------
  S$S01 <- mod %>% mrgsim(end = 30, delta = 0.25)

  ## S02  chronic SIADH, untreated — ADAPTED, not injured ---------------
  S$S02 <- start_at(mod, SI$state) %>% param(SI$param) %>% mrgsim(end = 30, delta = 0.25)

  ## S03  ACUTE hyponatraemia corrected fast (+20/24 h) — safe ----------
  S$S03 <- start_at(mod, ac_state) %>%
    param(c(list(AVPSIADH = 6, ACUTE = 0), corr(20))) %>% mrgsim(end = 30, delta = 0.25)

  ## S04  chronic, guideline correction +6/24 h -------------------------
  S$S04 <- start_at(mod, SI$state) %>% param(c(SI$param, corr(6))) %>%
    mrgsim(end = 30, delta = 0.25)

  ## S05  chronic, +12/24 h (the normal-risk ceiling) -------------------
  S$S05 <- start_at(mod, SI$state) %>% param(c(SI$param, corr(12))) %>%
    mrgsim(end = 30, delta = 0.25)

  ## S06  ** THE CENTRAL SCENARIO **  hypovolaemic patient given 0.9%
  ##      saline.  No hypertonic saline is prescribed at any point.
  S$S06 <- start_at(mod, HY$state) %>%
    param(c(HY$param, stop_losses, list(R09 = 2.0))) %>% mrgsim(end = 30, delta = 0.25)

  ## S07  ** THE COUNTERFACTUAL **  identical, but AVP is frozen at its
  ##      presentation value.  The difference between S06 and S07 is the
  ##      whole disease.
  S$S07 <- start_at(mod, HY$state) %>%
    param(c(HY$param, stop_losses,
            list(R09 = 2.0, AVPFREEZE = 1, TFREEZE = 0,
                 AVPFRZV = HY$state[which(names(init(mod)) == "AVP")]))) %>%
    mrgsim(end = 30, delta = 0.25)

  ## S08  proactive desmopressin clamp + titrated 3% saline -------------
  S$S08 <- start_at(mod, HY$state) %>%
    param(c(HY$param, stop_losses, list(R09 = 0.5), corr(6, cap = 130))) %>%
    mrgsim(events = c(ddavp_q8(0, 5),
                      ev(time = 5,   cmt = "DDAC", amt = 2, ii = 0.5, addl = 3),
                      ev(time = 7,   cmt = "DDAC", amt = 2, ii = 1.0, addl = 2)),
           end = 30, delta = 0.25)

  ## S09-S11  relowering rescue at 12, 24, 48 h -------------------------
  ## A rescue is a change of PRESCRIPTION, not just a change of infusion
  ## rate: the saline stops, the fluid restriction and the desmopressin
  ## start, and the dextrose is titrated.  Running those parameters from
  ## t = 0 would abolish the overcorrection they are meant to rescue, so
  ## each rescue is simulated as two phases.
  for (T in c(12, 24, 48)) {
    nm <- sprintf("S%02d_relower_%dh", 8 + which(c(12,24,48) == T), T)
    th <- T/24
    a  <- start_at(mod, HY$state) %>%
      param(c(HY$param, stop_losses, list(R09 = 2.0))) %>%
      mrgsim(end = th, delta = 0.05)
    st2 <- as.numeric(tail(a, 1)[, names(init(mod))])
    S[[nm]] <- start_at(mod, st2) %>%
      param(c(HY$param, stop_losses,
              list(R09 = 0, WIN = 0.5, RESCUE = 1, TRESCUE = th,
                   NARES = 118, DURRES = 1.0))) %>%
      mrgsim(events = ddavp_q8(th, th + 2.5), end = 60, delta = 0.25)
  }

  ## S12  potassium repletion alone — sodium rises with no sodium given -
  k_state <- SI$state
  k_state[which(names(init(mod)) == "KE")] <-
    k_state[which(names(init(mod)) == "KE")] - 450
  S$S12 <- start_at(mod, k_state) %>% param(c(SI$param, list(RKCL = 120, WIN = 1.0))) %>%
    mrgsim(end = 14, delta = 0.25)

  ## S13  tolvaptan 15 mg daily — the aquaresis you cannot switch off ---
  S$S13 <- start_at(mod, SI$state) %>% param(SI$param) %>%
    mrgsim(events = tolvaptan(7), end = 14, delta = 0.25)

  ## S14  alcoholic / malnourished at a "safe" +8/24 h ------------------
  S$S14 <- start_at(mod, SI$state) %>%
    param(c(SI$param, corr(8), list(FOSM = 0.55, FNUT = 0.85))) %>%
    mrgsim(end = 60, delta = 0.25)

  ## S15  normal-risk patient at the same +8/24 h (comparator) ----------
  S$S15 <- start_at(mod, SI$state) %>% param(c(SI$param, corr(8))) %>%
    mrgsim(end = 60, delta = 0.25)

  ## S16  oral urea 30 g/day --------------------------------------------
  S$S16 <- start_at(mod, SI$state) %>%
    param(c(SI$param, list(UREADOSE = 30, WIN = 1.0))) %>% mrgsim(end = 14, delta = 0.25)

  ## S17  beer potomania, solute intake restored (refeeding) ------------
  S$S17 <- start_at(mod, PO$state) %>%
    param(c(PO$param, list(TSOLUP = 0, PUREA2 = 400))) %>% mrgsim(end = 30, delta = 0.25)

  ## S18  slow but large: +8/day for four days (32 mmol/L total) --------
  S$S18 <- start_at(mod, SI$state) %>% param(c(SI$param, corr(8, cap = 142))) %>%
    mrgsim(end = 30, delta = 0.25)

  ## S19  adrenal insufficiency: hydrocortisone opens the aquaresis -----
  S$S19 <- start_at(mod, AD$state) %>%
    param(c(AD$param, list(TAVPOFF = 0, WIN = 1.0))) %>% mrgsim(end = 30, delta = 0.25)

  ## S20  cirrhosis / pre-transplant: the osmolyte buffer is already spent
  CI <- mod %>% param(c(SI$param, list(FOSM = 0.70, INS0 = 4.2, GLX0 = 19.5))) %>%
    mrgsim(end = 21, delta = 0.25)
  S$S20 <- start_at(mod, as.numeric(tail(CI, 1)[, names(init(mod))])) %>%
    param(c(SI$param, corr(8), list(FOSM = 0.70, INS0 = 4.2, GLX0 = 19.5))) %>%
    mrgsim(end = 60, delta = 0.25)

  S
}

## =====================================================================
##  THE ONE PLOT THAT MATTERS
##  S06 against S07: same patient, same fluid, one switch.
## =====================================================================
plot_central_experiment <- function(S) {
  d <- bind_rows(
    as_tibble(S$S06) %>% mutate(arm = "AVP responds (physiological)"),
    as_tibble(S$S07) %>% mutate(arm = "AVP frozen (counterfactual)")) %>%
    filter(time <= 7) %>%
    select(time, arm, SODIUM, UOSMOL, EFWCL, STRESSc, DEFICIT) %>%
    tidyr::pivot_longer(-c(time, arm))
  ggplot(d, aes(time, value, colour = arm)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~name, scales = "free_y") +
    labs(x = "days", y = NULL,
         title = "Who sets the correction rate?",
         subtitle = paste("Hypovolaemic hyponatraemia, 0.9% saline only.",
                          "No hypertonic saline in either arm.")) +
    theme_bw() + theme(legend.position = "bottom")
}

## =====================================================================
##  USAGE
##    S <- run_scenarios(mod)
##    plot_central_experiment(S)
##    plot(S$S06, SODIUM + UOSMOL + EFWCL + OMEGAc + STRESSc + DEFICIT ~ time)
##
##  Reference values for every scenario: ods_verification_output.txt
##  Interactive exploration: ods_shiny_app.R
## =====================================================================
