## ============================================================================
##  rfs_mrgsolve_model.R
##  REFEEDING SYNDROME  -  Quantitative Systems Pharmacology model (mrgsolve)
## ============================================================================
##
##  THE THESIS
##  ----------
##  Refeeding syndrome is not a disease of low serum phosphate.  It is a
##  FLUX MISMATCH, and the two fluxes are set by two different people.
##
##      Lambda_P  =  J_demand / J_supply
##
##      J_demand  the phosphate that must move into cells because the
##                glycolytic flux has raised the cellular organic-phosphate
##                set-point (G6P, F1,6-bisP, ATP, phosphocreatine, 2,3-DPG),
##                plus the phosphorus laid down in new lean tissue.
##                Set by the CLINICIAN, through the glucose infusion rate.
##
##      J_supply  absorbed dietary phosphate + prescribed intravenous
##                phosphate + net bone efflux + the renal phosphate saved by
##                switching urinary excretion off.
##                Set by the PATIENT'S HISTORY.
##
##  Three pieces of arithmetic make that ratio dangerous, and all three are
##  stoichiometry rather than fitted parameters.
##
##  (1) THE MEASURED COMPARTMENT IS NOT THE COMPARTMENT AT RISK
##          ECF inorganic phosphate  1.15 mmol/L x 12.4 L  =      14 mmol
##          total body phosphorus                          =  22,780 mmol
##      The laboratory samples 0.06 % of the pool.  A refeeding demand of
##      60 mmol/d is FOUR TIMES the entire measured compartment, per day.
##      This is why a normal admission phosphate is not reassurance - it is
##      the expected finding - and why NICE triages on history rather than
##      on the blood test.
##
##  (2) THE RENAL RESERVE IS SINGLE-USE.  Urinary phosphate excretion is a
##      THRESHOLD, GFR x max(0, Pser - TmP/GFR).  The kidney's entire
##      contribution is to stop excreting, which is worth about 14 mmol/d
##      and is fully spent within 24 hours.  After that it has nothing more
##      to offer however low the phosphate goes.
##
##  (3) BONE IS SLOW AND MAGNESIUM CAN SWITCH IT OFF.  Net skeletal efflux
##      is ~10 mmol/d and PTH-driven; hypomagnesaemia suppresses PTH
##      secretion, so the last endogenous supply line is disabled by a
##      second electrolyte that the same syndrome depletes.
##
##  Thiamine runs on a DIFFERENT CLOCK.  Both pools take about three weeks
##  to empty, but they REFILL two orders of magnitude apart, because at
##  supra-physiological plasma concentrations thiamine bypasses the
##  saturable ThTR-1/ThTR-2 carriers by passive diffusion, whereas
##  phosphate has to be pumped into 40 kg of cells against a saturable
##  transporter.
##
##  CALIBRATION POLICY
##  ------------------
##  Spent on NORMAL physiology (nothing on refeeding syndrome): total body
##  P/K/Mg and their bone/ICF/ECF partition; 67 mmol P and 72 mmol K per kg
##  fat-free mass; serum reference intervals; TmP/GFR; GFR 6 L/h; urinary
##  P 29, K 70, Mg 4.9 mmol/d in balance with a normal diet; a 26.5 mg
##  thiamine store with a 14 d biological half-life - which REPRODUCES the
##  1.1-1.4 mg/d RDA as an OUTPUT rather than taking it as an input;
##  Cunningham REE; the Forbes fat/lean partition; insulin t1/2 5 min.
##  Spent on NON-refeeding pharmacology: enteral formula content (22 mmol P,
##  38 mmol K, 6 mmol Mg per 1000 kcal); saturable oral thiamine uptake;
##  the 0.6-1.0 mmol/L fall in serum potassium produced by insulin-dextrose
##  given for hyperkalaemia; the Ca x P solubility product; QTc sensitivity.
##  Spent on REFEEDING SYNDROME ITSELF: THREE numbers - sIns (the size of
##  the flux-driven rise in the cellular phosphate set-point), kFill (how
##  fast cells move toward it), and HAZ (one global mortality-hazard scale).
##
##  VERIFICATION
##  ------------
##  Every equation below was independently re-implemented in Python/scipy
##  (`rfs_reference_model.py`, output in `rfs_reference_output.txt`).  That
##  exercise found and fixed five real defects, all of which are corrected
##  here:
##    1. the cellular phosphate set-point was referenced to a WELL-FED
##       control, so a starved patient refed at 30 kcal/kg/d never crossed
##       it and phosphate never fell - the disease mechanism was inert;
##    2. the set-point was then driven by INSULIN, which has the wrong sign:
##       when hypophosphataemia throttles glucose disposal, glucose and
##       insulin RISE, so the demand grew at exactly the moment the cell
##       could no longer meet it.  It is driven by glycolytic FLUX here;
##    3. renal potassium excretion had no depletion adaptation, which drove
##       intracellular potassium to 23 % of normal during the starvation
##       phase - incompatible with life;
##    4. PaCO2 was an unbounded ratio inside an exponential hazard, which
##       produced 100 % mortality in two scenarios as a pure artefact;
##    5. 75 % of insulin secretion sat in a glucose-independent basal term,
##       leaving insulin - and therefore the whole phosphate demand -
##       almost unresponsive to feeding.
##
##  48 ODEs, 16 treatment scenarios.
##
##  EDUCATIONAL / RESEARCH MODEL.  NOT FOR CLINICAL USE.
## ============================================================================

library(mrgsolve)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

## ============================================================================
##  THE MODEL
## ============================================================================

rfs_code <- '
$PARAM
// ---- phase and regimen switches -------------------------------------------
PHASE     = 2         // 0 healthy, 1 starvation, 2 refeeding
STARVE_FRAC = 0.32      // intake as a fraction of requirement during starvation ()
QUALITY   = 0.25      // micronutrient quality of the starvation diet (0-1)
ALCOHOL   = 0.00      // alcohol use (0-1)
DIURETIC  = 0.00      // loop/thiazide diuretic exposure (0-1)
GI_LOSS   = 0.20      // vomiting/purging/diarrhoeal losses (0-1)
GUT_FAIL  = 0.00      // fractional loss of intestinal uptake during refeeding ()

KCAL_START = 10.0      // starting energy delivery (kcal/kg/d)
KCAL_GOAL = 30.0      // target energy delivery (kcal/kg/d)
ADV_FRAC  = 0.20      // fraction of goal added per day ()
ADV_START = 2.00      // day on which advancement begins (d)
CHO_FRAC  = 0.50      // energy fraction from carbohydrate ()
FAT_FRAC  = 0.32      // energy fraction from fat ()
PRO_FRAC  = 0.18      // energy fraction from protein ()
ROUTE     = 0         // 0 enteral formula, 1 intravenous dextrose
FORM_P    = 22.0      // formula phosphate (mmol/1000 kcal)
FORM_K    = 38.0      // formula potassium (mmol/1000 kcal)
FORM_MG   = 6.00      // formula magnesium (mmol/1000 kcal)
FORM_TH   = 5.70      // formula thiamine (umol/1000 kcal)

P_DOSE    = 0.00      // intravenous phosphate (mmol/kg/d)
K_DOSE    = 0.00      // intravenous potassium (mmol/kg/d)
MG_DOSE   = 0.00      // intravenous magnesium (mmol/kg/d)
MG_BOLUS  = 0         // 1 gives the daily magnesium over 2 h instead of 24 h
TH_IV     = 0.00      // intravenous thiamine (mg/d)
TH_PO     = 0.00      // oral thiamine (mg/d)
TH_LEAD   = 0.00      // hours of thiamine given BEFORE the first glucose (h)
TH_DELAY  = 0.00      // hours by which thiamine is started LATE (h)
TH_DAYS   = 10.0      // duration of thiamine course (d)
P_DELAY   = 0.00      // days by which phosphate is started late (d)
RESCUE    = 0         // 1 enables reactive repletion at the usual thresholds
NA_INTAKE = 100       // sodium intake (mmol/d)

// ---- anthropometry ---------------------------------------------------------
BW0       = 62.0      // pre-morbid body weight (kg)
HT        = 1.68      // height (m)
FM0       = 14.0      // pre-morbid fat mass (kg)
FFM0      = 48.0      // pre-morbid fat-free mass (kg)
FECF      = 0.200     // ECF volume per kg body weight (L/kg)
FVG       = 0.160     // glucose distribution volume (L/kg)
FVP       = 0.045     // plasma volume (L/kg)
GFR0      = 6.00      // glomerular filtration rate (L/h)

// ---- phosphate -------------------------------------------------------------
RHOP      = 67.0      // intracellular phosphate per kg fat-free mass (mmol/kg)
PB0       = 19550     // bone phosphorus (mmol)
PSER0     = 1.15      // reference serum phosphate (mmol/L)
PTH50TM   = 4.20      // PTH for half-maximal TmP suppression (pmol/L)
FGF50TM   = 45.0      // FGF23 for half-maximal TmP suppression (pg/mL)
HTMP      = 1.30      // Hill coefficient of TmP suppression ()
OBLIGP    = 0.0833    // obligate non-threshold phosphate loss (mmol/h)
FABSP     = 0.65      // fractional intestinal phosphate absorption ()
KAP       = 0.25      // gut phosphate transit (1/h)
KMPUP     = 0.30      // Km of cellular Na-Pi cotransport on serum P (mmol/L)
SINS      = 0.750     // scale of the flux-driven cellular P set-point rise ()
USAT      = 0.50      // normalised glycolytic flux at half the set-point rise ()
KFILL     = 0.0145    // rate at which cells approach the set-point (1/h)
HBRES     = 1.60      // Hill coefficient of PTH on bone resorption ()
KPRECIP   = 0.55      // calcium-phosphate precipitation rate (1/h)
KSP       = 4.40      // calcium-phosphate solubility product (mmol2/L2)
ETOHP     = 0.32      // fractional fall of TmP/GFR at full alcohol use ()
KMGLY     = 0.25      // phosphate Km of glycolytic (GAPDH) flux (mmol/L)

// ---- potassium -------------------------------------------------------------
RHOK      = 72.0      // intracellular potassium per kg fat-free mass (mmol/kg)
KSER0     = 4.20      // reference serum potassium (mmol/L)
KKP       = 8.60      // Na/K-ATPase pump flux scale (L/h)
EMAXK     = 0.34      // maximal fractional rise of pump activity ()
TAUX      = 1.20      // insulin effect scale of the pump response ()
FABSK     = 0.90      // fractional intestinal potassium absorption ()
KAK       = 0.40      // gut potassium transit (1/h)
BETAMGK   = 1.90      // renal K wasting per unit fractional hypomagnesaemia ()
KOBL      = 0.25      // irreducible potassium loss (mmol/h)
HROMK     = 6.00      // steepness of ROMK down-regulation in K depletion ()
ROMKMIN   = 0.05      // floor on distal potassium secretion ()
RALD      = 0.30      // aldosterone contribution to sodium retention ()

// ---- magnesium -------------------------------------------------------------
RHOMG     = 7.30      // exchangeable magnesium per kg fat-free mass (mmol/kg)
MGSER0    = 0.85      // reference serum magnesium (mmol/L)
FUFMG     = 0.70      // ultrafilterable magnesium fraction ()
FABSMG    = 0.40      // fractional intestinal magnesium absorption ()
KAMG      = 0.35      // gut magnesium transit (1/h)
KMGP      = 1.30      // insulin-stimulated magnesium uptake scale (L/h)
EMAXMG    = 0.20      // maximal fractional rise of magnesium uptake ()
ETOHMG    = 0.28      // fractional fall of TmMg at full alcohol use ()

// ---- calcium / PTH / vitamin D / FGF23 -------------------------------------
CA0       = 1.20      // reference ionised calcium (mmol/L)
PTH0      = 4.20      // reference PTH (pmol/L)
PTHMAX    = 26.0      // maximal PTH secretion (pmol/L)
PTHMIN    = 0.60      // minimal PTH secretion (pmol/L)
HCA       = 6.00      // steepness of the PTH-calcium sigmoid ()
KPTHEL    = 6.93      // PTH elimination (1/h)
MGPTH50   = 0.22      // magnesium for half-maximal PTH secretion (mmol/L)
CTD0      = 120       // reference 1,25(OH)2 vitamin D (pmol/L)
KCTD      = 0.10      // 1,25D turnover (1/h)
FGF0      = 45.0      // reference FGF23 (pg/mL)
KFGF      = 0.14      // FGF23 turnover (1/h)

// ---- thiamine --------------------------------------------------------------
THT0      = 100.0     // whole-body thiamine store (umol)
KTHDEG    = 0.00206   // thiamine degradation, ln2/14 d (1/h)
KTHCHO    = 2.05e-4   // thiamine consumed per mmol glucose disposed (umol/mmol)
THP0C     = 10.0      // reference plasma thiamine (nmol/L)
CLTH      = 0.42      // plasma thiamine clearance above threshold (L/h)
THRENALC  = 30.0      // plasma thiamine above which renal loss is free (nmol/L)
VMAXTHA   = 7.50      // saturable intestinal thiamine uptake (umol/h)
KMTHA     = 6.00      // gut thiamine Km (umol)
KMTHT     = 12.0      // tissue carrier Km (nmol/L)
KTHPASS   = 0.0042    // PASSIVE tissue uptake clearance (L/h)
KTHOFF    = 0.020     // tissue thiamine release (1/h)
KMTPP     = 0.35      // fractional store for half enzyme activity ()
HTPP      = 3.00      // Hill coefficient of TPP-dependent activity ()
ETOHABS   = 0.45      // fractional thiamine absorption at full alcohol use ()
ETOHTPK   = 0.55      // thiamine pyrophosphokinase activity at full alcohol use ()

// ---- glucose / insulin ------------------------------------------------------
G0        = 5.00      // reference plasma glucose (mmol/L)
INS0      = 45.0      // reference plasma insulin (pmol/L)
INSLOW    = 15.0      // insulin at the starvation glycaemia (pmol/L)
GLOW      = 3.50      // glycaemia of prolonged starvation (mmol/L)
EGP0      = 44.0      // maximal hepatic glucose output (mmol/h)
X50       = 1.00      // insulin effect for half EGP suppression ()
RD0       = 5.00      // insulin-independent glucose clearance (L/h)
GTHR      = 3.00      // glucose secretion threshold (mmol/L)
HSEC      = 2.00      // Hill coefficient of glucose-stimulated secretion ()
CLINS     = 8.32      // insulin elimination (1/h)
P2        = 1.50      // insulin effect-compartment rate (1/h)
KACHO     = 0.55      // gut carbohydrate emptying (1/h)
KAFAT     = 0.28      // gut fat emptying (1/h)
KAPRO     = 0.38      // gut protein emptying (1/h)
GLY0      = 2200      // reference glycogen (mmol glucosyl)
GLYMAX    = 3300      // maximal glycogen (mmol glucosyl)
JGLY      = 12.0      // basal glycogen cycling flux (mmol/h)
GCG0      = 60.0      // reference glucagon (ng/L)
KGCG      = 0.60      // glucagon turnover (1/h)
SIMIN     = 0.45      // insulin sensitivity floor of the adapted starved state ()
FINC      = 0.35      // incretin share of fed insulin secretion ()
INCIV     = 0.45      // residual incretin drive of intravenous dextrose ()
RALOWR    = 0.32      // intake fraction used for the two-point insulin calibration ()

// ---- lipid / lactate --------------------------------------------------------
FFA0      = 0.45      // reference plasma NEFA (mmol/L)
KFFAOX    = 1.65      // NEFA oxidation (1/h)
BHB0      = 0.10      // reference ketones (mmol/L)
KBHBOX    = 0.62      // ketone oxidation (1/h)
LAC0      = 1.00      // reference lactate (mmol/L)
CLLAC     = 22.0      // lactate clearance (L/h)
KLAC2     = 1.50      // lactate produced per unit blocked pyruvate flux ()

// ---- organ energetics -------------------------------------------------------
KATPM     = 0.14      // myocardial ATP index relaxation (1/h)
KPATP     = 0.35      // serum P for half-maximal oxidative phosphorylation (mmol/L)
HPATP     = 2.00      // Hill coefficient of the phosphate energetic terms ()
KPCR      = 0.30      // phosphocreatine turnover (1/h)
KPPCR     = 0.50      // serum P for half-maximal phosphocreatine (mmol/L)
KDPG      = 0.10      // 2,3-DPG turnover (1/h)
KPDPG     = 0.55      // serum P for half-maximal 2,3-DPG (mmol/L)
KATPD     = 0.16      // diaphragm energetic turnover (1/h)

// ---- composition / heart / volume -------------------------------------------
RHOFAT    = 9440      // energy density of adipose tissue (kcal/kg)
RHOLEAN   = 1800      // energy density of lean tissue (kcal/kg)
PAL       = 1.30      // physical activity level ()
KAT       = 0.0060    // adaptive thermogenesis rate (1/h)
ATMIN     = 0.78      // floor of adaptive thermogenesis ()
LVM0      = 150       // pre-morbid left ventricular mass (g)
KLVM      = 0.00198   // LV remodelling rate, 3-week constant (1/h)
NAREF     = 140       // plasma sodium (mmol/L)
RINS      = 0.62      // maximal fractional sodium retention from insulin ()
KNAB      = 0.0139    // natriuresis of the retained sodium load (1/h)
KALD      = 0.06      // aldosterone turnover (1/h)

// ---- hazards / gas exchange -------------------------------------------------
QTC0      = 400       // reference QTc (ms)
AK        = 26.0      // QTc per mmol/L of hypokalaemia (ms)
AMG       = 22.0      // QTc per mmol/L of hypomagnesaemia (ms)
ACA       = 95.0      // QTc per mmol/L of hypocalcaemia (ms)
H0ARR     = 5.0e-6    // baseline arrhythmia hazard (1/h)
CQT       = 0.045     // hazard per ms of QTc above 440 ()
CATPARR   = 3.20      // hazard amplification by myocardial ATP loss ()
H0CF      = 3.0e-6    // baseline cardiac-failure hazard (1/h)
PCONG     = 2.00      // exponent of the congestion index ()
H0RF      = 2.0e-6    // baseline respiratory-failure hazard (1/h)
PACO2REF  = 55.0      // PaCO2 above which the hazard rises (mmHg)
QCO2      = 0.075     // hazard per mmHg above the reference ()
H0WE      = 3.0e-4    // Wernicke incidence hazard scale (1/h)
HAZMAX    = 0.020     // ceiling on the instantaneous total hazard (1/h)
VAMIN     = 0.35      // floor on relative alveolar ventilation capacity ()
CONTMIN   = 0.20      // floor on the myocardial contractility index ()
HAZ       = 1.00      // global mortality hazard scale ()
KCALO2    = 4.83      // energy per litre of oxygen (kcal/L)
PACO2_0   = 40.0      // reference PaCO2 (mmHg)
CHOOXMAX  = 4.60      // maximal carbohydrate oxidation (mg/kg/min)



$CMT // 48 compartments
// A_CHO    luminal carbohydrate (g)
// A_FAT    luminal fat (g)
// A_PRO    luminal protein (g)
// GLU      plasma glucose (mmol/L)
// INS      plasma insulin (pmol/L)
// X        insulin effect compartment ()
// GCG      glucagon (ng/L)
// GLY      glycogen (mmol glucosyl)
// FFA      plasma NEFA (mmol/L)
// BHB      ketones (mmol/L)
// LAC      lactate (mmol/L)
// PGUT     luminal phosphate (mmol)
// PE       ECF phosphate (mmol)
// PI       intracellular phosphate (mmol)
// PB       bone phosphorus (mmol)
// KGUT     luminal potassium (mmol)
// KE       ECF potassium (mmol)
// KI       intracellular potassium (mmol)
// MGGUT    luminal magnesium (mmol)
// MGE      ECF magnesium (mmol)
// MGI      exchangeable magnesium (mmol)
// CAE      ECF ionised calcium (mmol)
// CAB      precipitated / deposited calcium (mmol)
// PTH      parathyroid hormone (pmol/L)
// CTD      1,25(OH)2 vitamin D (pmol/L)
// FGF      FGF23 (pg/mL)
// THGUT    luminal thiamine (umol)
// THP      plasma thiamine (umol)
// THT      tissue thiamine store (umol)
// ATPM     myocardial ATP index ()
// PCRM     myocardial phosphocreatine index ()
// DPG      erythrocyte 2,3-DPG index ()
// ATPD     diaphragm energetic index ()
// FM       fat mass (kg)
// FFM      fat-free mass (kg)
// AT       adaptive thermogenesis ()
// LVM      left ventricular mass (g)
// NAB      cumulative sodium balance (mmol)
// ALD      aldosterone ()
// H_ARR    cumulative arrhythmia hazard ()
// H_CF     cumulative cardiac-failure hazard ()
// H_RF     cumulative respiratory-failure hazard ()
// H_WE     cumulative Wernicke incidence ()
// SURV     survival fraction ()
// AUCP     exposure below 0.65 mmol/L (mmol/L*h)
// CUMP     cumulative phosphate delivered (mmol)
// CUMK     cumulative potassium delivered (mmol)
// CUMKCAL  cumulative energy delivered (kcal)
  A_CHO A_FAT A_PRO GLU INS X
  GCG GLY FFA BHB LAC PGUT
  PE PI PB KGUT KE KI
  MGGUT MGE MGI CAE CAB PTH
  CTD FGF THGUT THP THT ATPM
  PCRM DPG ATPD FM FFM AT
  LVM NAB ALD H_ARR H_CF H_RF
  H_WE SURV AUCP CUMP CUMK CUMKCAL

$GLOBAL
#define SQ(a) ((a)*(a))
namespace {
  inline double posv(double x){ return x > 0.0 ? x : 0.0; }
  inline double hillf(double x, double k, double h){
    double xx = posv(x);
    return pow(xx,h) / (pow(k,h) + pow(xx,h));
  }
  // Smooth "is x below thr".  A hard threshold makes the integrator chatter
  // on the crossing, and real prescribing is not instantaneous either.
  inline double swbelow(double x, double thr, double w){
    double z = (x - thr)/w;
    if(z >  40.0) return 0.0;
    if(z < -40.0) return 1.0;
    return 1.0/(1.0 + exp(z));
  }
  inline double pumpm(double X, double Emax, double tau){
    if(X >= 1.0) return 1.0 + Emax*(1.0 - exp(-(X-1.0)/tau));
    return 1.0 - 0.6*Emax*(1.0 - exp(-(1.0-X)/tau));
  }
  inline double reekcal(double FFM, double AT){ return (370.0 + 21.6*FFM)*AT; }
}
// derived constants, recomputed in $MAIN so they always track $PARAM
double VECF0, VG0, VP0, VLAC, LVMFRAC, PI0, KI0, MGI0, PE0, KE0, MGE0, CAE0;
double REE0, TEE0, KCALH0, RA0, SI, UINS0, KSEC, SECB, KINC, KGLYS, KGLYB;
double KLIPO, KKETO, CA50, TMPMAX, TMP0, KBRES, KBFORM, JBONE0, KKLEAK, KKSEC;
double KMGLEAK, TMMG, KCAU, KCABONE, VMAXTHT, THP0A, THABS0, FPDH0, FORMTH0;
double LACPROD0, VO2_0, VCO2REF, FU1;

$MAIN
// ---------------------------------------------------------------------------
//  DERIVED CONSTANTS.  Every one of these is fixed by a BALANCE CONDITION on
//  normal physiology, not fitted to refeeding data: the healthy state is a
//  numerically exact steady state of the 44 physiological equations below
//  (verified at max|dy/dt| = 6e-15 in the Python reference).
// ---------------------------------------------------------------------------
VECF0 = FECF*BW0;  VG0 = FVG*BW0;  VP0 = FVP*BW0;  VLAC = 0.55*BW0;
LVMFRAC = LVM0/FFM0;
PI0 = RHOP*FFM0;  KI0 = RHOK*FFM0;  MGI0 = RHOMG*FFM0;
PE0 = PSER0*VECF0; KE0 = KSER0*VECF0; MGE0 = MGSER0*VECF0; CAE0 = CA0*VECF0;

REE0 = reekcal(FFM0, 1.0);
TEE0 = REE0*PAL/0.90;              // includes 10 % diet-induced thermogenesis
KCALH0 = TEE0/24.0;

RA0 = (KCALH0*0.50/4.0)*1000.0/180.0;
double EGPSS = EGP0/(1.0 + 1.0/X50);
SI    = ((RA0 + EGPSS)/G0 - RD0);
UINS0 = SI*G0;
FU1   = 1.0/(USAT + 1.0);

double CVP = CLINS*VP0;
double A1  = pow(posv(G0   - GTHR), HSEC);
double A2  = pow(posv(GLOW - GTHR), HSEC);
KINC = FINC*CVP*INS0/RA0;
double B1 = CVP*INS0   - KINC*RA0;
double B2 = CVP*INSLOW - KINC*RALOWR*RA0;
KSEC = (B1 - B2)/(A1 - A2);
SECB = B1 - KSEC*A1;

KGLYS = JGLY/(G0*(1.0 - GLY0/GLYMAX));
KGLYB = JGLY*(1.0 + 3.0)/GLY0;
KLIPO = KFFAOX*FFA0*(1.0 + 3.2);
KKETO = KBHBOX*BHB0*(1.0 + 6.0)/FFA0;

double RR = (PTHMAX - PTHMIN)/(PTH0 - PTHMIN) - 1.0;
CA50 = CA0/pow(RR, 1.0/HCA);

double FP0 = 22.0*KCALH0/1000.0;
double JPABS0 = FABSP*FP0;
TMP0 = PSER0 - (JPABS0 - OBLIGP)/GFR0;
TMPMAX = TMP0*((1.0 + hillf(PTH0, PTH50TM, HTMP))
              *(1.0 + 0.55*hillf(FGF0, FGF50TM, HTMP)));
JBONE0 = 10.0/24.0;
KBRES  = JBONE0/(PB0*hillf(PTH0, PTH0, HBRES));
KBFORM = JBONE0/PB0;

KKLEAK = (KKP*KSER0)/KI0;
double FK0 = 38.0*KCALH0/1000.0;
double EXK0 = FABSK*FK0;
KKSEC = (EXK0 - KOBL)/KSER0;

KMGLEAK = (KMGP*MGSER0)/MGI0;
double FMG0 = 6.0*KCALH0/1000.0;
double EXMG0 = FABSMG*FMG0;
TMMG = FUFMG*MGSER0 - (EXMG0 - 0.5/24.0)/GFR0;

KCAU = (5.0/24.0)/CA0;
KCABONE = (5.0/24.0)/CAE0;

THP0A = THP0C*VP0/1000.0;
double LOSSTH = KTHDEG*THT0 + KTHCHO*(RA0 + EGPSS);
double NEEDTH = LOSSTH + KTHOFF*THT0;
VMAXTHT = (NEEDTH - KTHPASS*THP0C)*(KMTHT + THP0C)/THP0C;
double JUP0 = VMAXTHT*THP0C/(KMTHT + THP0C) + KTHPASS*THP0C;
THABS0 = JUP0 + CLTH*posv(THP0C - THRENALC)/1000.0 - KTHOFF*THT0;
FORMTH0 = THABS0/(KCALH0/1000.0);
FPDH0 = hillf(1.0, KMTPP, HTPP);
LACPROD0 = CLLAC*LAC0;

VO2_0 = TEE0/(1440.0*KCALO2);
VCO2REF = 0.85*VO2_0;

// ---------------------------------------------------------------------------
//  HEALTHY INITIAL CONDITIONS
// ---------------------------------------------------------------------------
A_CHO_0 = (KCALH0*0.50/4.0)/KACHO;
A_FAT_0 = (KCALH0*0.32/9.0)/KAFAT;
A_PRO_0 = (KCALH0*0.18/4.0)/KAPRO;
GLU_0 = G0; INS_0 = INS0; X_0 = 1.0; GCG_0 = GCG0; GLY_0 = GLY0;
FFA_0 = FFA0; BHB_0 = BHB0; LAC_0 = LAC0;
PGUT_0 = (22.0*KCALH0/1000.0)/KAP;  PE_0 = PE0; PI_0 = PI0; PB_0 = PB0;
KGUT_0 = (38.0*KCALH0/1000.0)/KAK;  KE_0 = KE0; KI_0 = KI0;
MGGUT_0 = (6.0*KCALH0/1000.0)/KAMG; MGE_0 = MGE0; MGI_0 = MGI0;
CAE_0 = CAE0; PTH_0 = PTH0; CTD_0 = CTD0; FGF_0 = FGF0;
double THGF = FORMTH0*KCALH0/1000.0;
THGUT_0 = KMTHA*THGF/((VMAXTHA - THGF) > 1e-9 ? (VMAXTHA - THGF) : 1e-9);
THP_0 = THP0A; THT_0 = THT0;
ATPM_0 = 1.0; PCRM_0 = 1.0; DPG_0 = 1.0; ATPD_0 = 1.0;
FM_0 = FM0; FFM_0 = FFM0; AT_0 = 1.0; LVM_0 = LVM0;
NAB_0 = 0.0; ALD_0 = 1.0; SURV_0 = 1.0;

$ODE
// ---------------------------------------------------------------------------
double fm = posv(FM), ffm = posv(FFM);
double BW = fm + ffm;  if(BW < 20.0) BW = 20.0;
double Vecf = FECF*BW0 + NAB/NAREF;  if(Vecf < 1.0) Vecf = 1.0;
double Vg = FVG*BW, Vp = FVP*BW;

double Pser  = posv(PE)/Vecf;
double Kser  = posv(KE)/Vecf;
double MGser = posv(MGE)/Vecf;
double Caser = posv(CAE)/Vecf;
double Cth   = posv(THP)/Vp*1000.0;
double day   = SOLVERTIME/24.0;

// ---- 1. nutrition ----------------------------------------------------------
double kcal_d, cho_f, fat_f, pro_f, quality;
int route = 0;
if(PHASE < 0.5){
  kcal_d = TEE0; cho_f = 0.50; fat_f = 0.32; pro_f = 0.18; quality = 1.0;
} else if(PHASE < 1.5){
  kcal_d = STARVE_FRAC*TEE0; cho_f = 0.50; fat_f = 0.32; pro_f = 0.18;
  quality = QUALITY;
} else {
  double kk;
  if(day < ADV_START) kk = KCAL_START;
  else {
    double n = floor(day - ADV_START) + 1.0;
    kk = KCAL_START + n*ADV_FRAC*KCAL_GOAL;
    if(kk > KCAL_GOAL) kk = KCAL_GOAL;
  }
  kcal_d = kk*BW; cho_f = CHO_FRAC; fat_f = FAT_FRAC; pro_f = PRO_FRAC;
  quality = 1.0; route = (ROUTE > 0.5) ? 1 : 0;
  if(SOLVERTIME < TH_LEAD) kcal_d = 0.0;   // thiamine lead-in, no glucose yet
}
double kcal_h = kcal_d/24.0;
double g_cho = kcal_h*cho_f/4.0, g_fat = kcal_h*fat_f/9.0, g_pro = kcal_h*pro_f/4.0;
double sc = kcal_h/1000.0;

// Diet QUALITY hits the nutrients unequally: refined carbohydrate and alcohol
// are almost devoid of thiamine and poor in phosphate, but any real food
// still carries potassium roughly in proportion to its mass.
double qP = 0.40 + 0.60*quality, qK = 0.70 + 0.30*quality;
double qMg = 0.45 + 0.55*quality, qTh = quality;
double fP, fK, fMg, fTh;
if(route == 1){ fP = 0.0; fK = 0.0; fMg = 0.0; fTh = 0.0; }
else { fP = FORM_P*sc*qP; fK = FORM_K*sc*qK; fMg = FORM_MG*sc*qMg; fTh = FORM_TH*sc*qTh; }
if(PHASE < 0.5){ fP = 22.0*sc; fK = 38.0*sc; fMg = 6.0*sc; fTh = FORMTH0*sc; }

double alc = (PHASE > 0.5 && PHASE < 1.5) ? ALCOHOL : 0.0;
double absTh = 1.0 - alc*(1.0 - ETOHABS);
double tpk   = 1.0 - alc*(1.0 - ETOHTPK);
if(PHASE > 1.5) absTh *= (1.0 - GUT_FAIL);

// ---- 2. prescribed repletion -----------------------------------------------
double ivP = 0.0, ivK = 0.0, ivMg = 0.0, ivTh = 0.0, poTh = 0.0;
if(PHASE > 1.5){
  ivP = (day >= P_DELAY) ? P_DOSE*BW/24.0 : 0.0;
  ivK = K_DOSE*BW/24.0;
  if(MG_BOLUS > 0.5) ivMg = (fmod(SOLVERTIME,24.0) < 2.0) ? MG_DOSE*BW/2.0 : 0.0;
  else               ivMg = MG_DOSE*BW/24.0;
  if(SOLVERTIME >= TH_DELAY && SOLVERTIME < TH_DELAY + TH_DAYS*24.0){
    ivTh = TH_IV/265.4*1000.0/24.0;
    poTh = TH_PO/265.4*1000.0/24.0;
  }
  if(RESCUE > 0.5){
    ivP  += 0.50*BW/24.0*swbelow(Pser , 0.65, 0.030);
    ivK  += 1.00*BW/24.0*swbelow(Kser , 3.30, 0.100);
    ivMg += 0.20*BW/24.0*swbelow(MGser, 0.65, 0.030);
  }
}

// ---- 3. glucose / insulin ---------------------------------------------------
double Ra_gut;
if(route == 1 && PHASE > 1.5){ Ra_gut = g_cho*1000.0/180.0; dxdt_A_CHO = -KACHO*A_CHO; }
else { dxdt_A_CHO = g_cho - KACHO*A_CHO; Ra_gut = KACHO*posv(A_CHO)*1000.0/180.0; }
dxdt_A_FAT = g_fat - KAFAT*A_FAT;
dxdt_A_PRO = g_pro - KAPRO*A_PRO;

double gly_f = posv(GLY)/(0.35*GLY0);  if(gly_f > 1.0) gly_f = 1.0;
double EGP = EGP0/(1.0 + X/X50)*(0.7 + 0.3*GCG/GCG0)*(0.45 + 0.55*gly_f);
if(EGP < 0.0) EGP = 0.0;

// Prolonged starvation down-regulates insulin-stimulated glucose disposal.
// AT is the slow marker of that adapted state, so sensitivity recovers over
// about a week - which is why refeeding hyperglycaemia is real and transient.
double ATn = (AT - ATMIN)/(1.0 - ATMIN);  if(ATn < 0.0) ATn = 0.0; if(ATn > 1.0) ATn = 1.0;
double SIe = SI*(SIMIN + (1.0 - SIMIN)*ATn);
// Glycolysis needs inorganic phosphate: GAPDH consumes Pi stoichiometrically,
// so profound hypophosphataemia throttles glucose disposal itself.  This is
// the brake that stops serum phosphate reaching absurd values - the demand is
// extinguished by its own consequence.
double fPgly = (Pser/(KMGLY + Pser))/(PSER0/(KMGLY + PSER0));
if(fPgly > 1.0) fPgly = 1.0;
SIe *= fPgly;

double Uins = SIe*X*GLU;
double Rd = RD0*GLU + Uins;
double glyS = KGLYS*X*GLU*posv(1.0 - GLY/GLYMAX);
double glyB = KGLYB*GLY/(1.0 + 3.0*X);
double renG = GFR0*posv(GLU - 10.0);
dxdt_GLU = (Ra_gut + EGP + glyB - Rd - glyS - renG)/Vg;

double beta = 0.55 + 0.45*((ffm/FFM0 < 1.0) ? ffm/FFM0 : 1.0);
double inc = KINC*Ra_gut*((route == 1 && PHASE > 1.5) ? INCIV : 1.0);
double sec = SECB + KSEC*pow(posv(GLU - GTHR), HSEC) + inc;
dxdt_INS = (sec*beta - CLINS*INS*Vp)/Vp;
dxdt_X   = P2*(INS/INS0 - X);
double gcg_t = GCG0*(1.0 + 1.6*posv(G0 - GLU)/G0)*(1.0 + 0.9)/(1.0 + 0.9*X);
dxdt_GCG = KGCG*(gcg_t - GCG);
dxdt_GLY = glyS - glyB;
dxdt_FFA = KLIPO/(1.0 + 3.2*X)*(fm/FM0 + 0.15) - KFFAOX*FFA;
dxdt_BHB = KKETO*FFA/(1.0 + 6.0*X) - KBHBOX*BHB;

// ---- 4. thiamine -------------------------------------------------------------
double Jth_gut = VMAXTHA*posv(THGUT)/(KMTHA + posv(THGUT));
dxdt_THGUT = fTh + poTh - Jth_gut;
double Jth_abs = absTh*Jth_gut;
double Jth_up = (VMAXTHT*Cth/(KMTHT + Cth) + KTHPASS*Cth)*tpk;
double ren_th = CLTH*posv(Cth - THRENALC)/1000.0;
dxdt_THP = Jth_abs + ivTh - Jth_up + KTHOFF*THT - ren_th;
dxdt_THT = Jth_up - KTHOFF*THT - KTHDEG*THT - KTHCHO*Rd;

double thn = posv(THT)/THT0;
double fPDH = hillf(thn, KMTPP, HTPP)/FPDH0;  if(fPDH > 1.0) fPDH = 1.0;
double fTK = fPDH;

double lac_prod = (LACPROD0 + KLAC2*Rd*(1.0 - fPDH))*(1.0 + 0.85*posv(1.0 - DPG));
dxdt_LAC = (lac_prod - CLLAC*LAC)/VLAC;

// ---- 5. energy balance and body composition ----------------------------------
double REE = reekcal(ffm, AT);
double TEE = REE*PAL + 0.10*kcal_d;
double Ebal = kcal_d - TEE;
double cat = 1.0 - kcal_d/((TEE > 1.0) ? TEE : 1.0);
if(cat < 0.0) cat = 0.0;  if(cat > 1.0) cat = 1.0;
double fLean = 0.10 + 0.50*exp(-fm/8.0);
if(fLean < 0.08) fLean = 0.08;  if(fLean > 0.55) fLean = 0.55;
double dFFM_d = fLean*Ebal/RHOLEAN;
double dFM_d  = (1.0 - fLean)*Ebal/RHOFAT;
if(Ebal > 0.0){
  double a1 = Pser/0.70;  if(a1 > 1.0) a1 = 1.0;
  double a2 = Kser/3.50;  if(a2 > 1.0) a2 = 1.0;
  dFFM_d *= a1*a2;
}
dxdt_FFM = dFFM_d/24.0;
dxdt_FM  = dFM_d/24.0;
double rat = kcal_d/((TEE > 1.0) ? TEE : 1.0);  if(rat > 1.0) rat = 1.0;
dxdt_AT  = KAT*(ATMIN + (1.0 - ATMIN)*rat - AT);
dxdt_LVM = KLVM*(LVMFRAC*ffm - LVM);

// ---- 6. PHOSPHATE  -- the engine of the syndrome ------------------------------
// The set-point is driven by the GLYCOLYTIC FLUX actually achieved, not by the
// insulin signal: the pool being refilled IS the pool of phosphorylated
// intermediates, so its target must track the flux through them.
double un = Uins/UINS0;
double fU = un/(USAT + un);
double PI_t = RHOP*ffm*(1.0 + SINS*(fU - FU1));
double gap = PI_t - PI;
double fUp = Pser/(KMPUP + Pser);
double JPnet = KFILL*gap*((gap > 0.0) ? fUp : 1.0);

double Jres = KBRES*posv(PB)*hillf(PTH, PTH0, HBRES);
double Jform = KBFORM*posv(PB)
             * ((Pser/PSER0 < 2.0) ? Pser/PSER0 : 2.0)
             * ((Caser/CA0 < 2.0) ? Caser/CA0 : 2.0);

// Renal excretion is a THRESHOLD.  The kidney's whole contribution is to stop
// excreting, and a starved patient has already done that.
double TmP = TMPMAX/((1.0 + hillf(PTH, PTH50TM, HTMP))
                    *(1.0 + 0.55*hillf(FGF, FGF50TM, HTMP)));
TmP *= (1.0 - ETOHP*alc);
double ExP = GFR0*posv(Pser - TmP) + OBLIGP*Pser/PSER0;
double gi_P = (PHASE > 0.5 && PHASE < 1.5) ? GI_LOSS*4.0/24.0 : 0.0;

double JPabs = FABSP*KAP*posv(PGUT)*(0.7 + 0.3*CTD/CTD0);
dxdt_PGUT = fP - KAP*PGUT;

double sup = posv(Caser*Pser - KSP);
double Jprec = KPRECIP*sup*Vecf;

dxdt_PE = JPabs + ivP + Jres - JPnet - Jform - ExP - Jprec - gi_P;
dxdt_PI = JPnet;
dxdt_PB = Jform - Jres;

// ---- 7. POTASSIUM --------------------------------------------------------------
double depK = posv(1.0 - KI/(RHOK*((ffm > 1.0) ? ffm : 1.0)));
double JKin = KKP*pumpm(X, EMAXK, TAUX)*Kser*(1.0 + 0.8*depK);
double JKout = KKLEAK*posv(KI)*(1.0 + 0.9*cat);
double JKnew = dxdt_FFM*RHOK;
double fMgK = 1.0 + BETAMGK*posv(1.0 - MGser/MGSER0);
// ROMK/BK down-regulation: the distal nephron adapts to potassium depletion
// over days, but never all the way to zero - which is why a depleted patient
// still loses potassium in the urine.
double romk = pow(KI/(RHOK*((ffm > 1.0) ? ffm : 1.0)), HROMK);
if(romk > 1.0) romk = 1.0;  if(romk < ROMKMIN) romk = ROMKMIN;
double ExK = KOBL + KKSEC*Kser*ALD*fMgK*(1.0 + 1.4*DIURETIC)*romk;
double JKabs = FABSK*KAK*posv(KGUT);
dxdt_KGUT = fK - KAK*KGUT;
double gi_K = (PHASE > 0.5 && PHASE < 1.5) ? GI_LOSS*25.0/24.0 : 0.0;
dxdt_KE = JKabs + ivK + JKout - JKin - ExK - gi_K - JKnew;
dxdt_KI = JKin - JKout + JKnew;

// ---- 8. MAGNESIUM ---------------------------------------------------------------
double depMg = posv(1.0 - MGI/(RHOMG*((ffm > 1.0) ? ffm : 1.0)));
double JMgin = KMGP*pumpm(X, EMAXMG, TAUX)*MGser*(1.0 + 0.7*depMg);
double JMgout = KMGLEAK*posv(MGI)*(1.0 + 0.6*cat);
double JMgnew = dxdt_FFM*RHOMG;
double TmMg = TMMG*(1.0 - ETOHMG*alc);
double ExMg = (GFR0*posv(FUFMG*MGser - TmMg) + 0.5/24.0*MGser/MGSER0)
            * (1.0 + 1.8*DIURETIC);
double JMgabs = FABSMG*KAMG*posv(MGGUT);
dxdt_MGGUT = fMg - KAMG*MGGUT;
double gi_Mg = (PHASE > 0.5 && PHASE < 1.5) ? GI_LOSS*6.0/24.0 : 0.0;
dxdt_MGE = JMgabs + ivMg + JMgout - JMgin - ExMg - gi_Mg - JMgnew;
dxdt_MGI = JMgin - JMgout + JMgnew;

// ---- 9. CALCIUM / PTH / 1,25D / FGF23 --------------------------------------------
double fMgPTH = (MGser/(MGPTH50 + MGser))/(MGSER0/(MGPTH50 + MGSER0));
if(fMgPTH > 1.05) fMgPTH = 1.05;
double secPTH = (PTHMIN + (PTHMAX - PTHMIN)/(1.0 + pow(Caser/CA50, HCA)))*fMgPTH;
dxdt_PTH = KPTHEL*(secPTH - PTH);

double CaRes = KCABONE*CAE0*hillf(PTH, PTH0, HBRES)*2.0;
double CaForm = KCABONE*posv(CAE)*((Pser/PSER0 < 2.0) ? Pser/PSER0 : 2.0);
double kr = kcal_h/KCALH0;  if(kr > 2.0) kr = 2.0;
double CaAbs = (5.0/24.0)*(CTD/CTD0)*(0.4 + 0.6*kr);
double CaU = KCAU*Caser*(1.0 + 0.4*posv(1.0 - PTH/PTH0));
dxdt_CAE = CaAbs + CaRes - CaForm - CaU - 1.5*Jprec;
dxdt_CAB = Jprec + CaForm - CaRes;
dxdt_CTD = KCTD*(CTD0*(0.5 + 0.5*PTH/PTH0)*(1.0 + 0.5*posv(1.0 - Pser/PSER0)) - CTD);
dxdt_FGF = KFGF*(FGF0*(0.25 + 0.75*Pser/PSER0)*(1.0 + 0.4*(CTD/CTD0 - 1.0)) - FGF);

// ---- 10. ORGAN ENERGETICS ---------------------------------------------------------
double nP = hillf(PSER0, KPATP, HPATP);
double fPox = hillf(Pser, KPATP, HPATP)/nP;
double oxid = 0.35 + 0.65*fPDH;  if(oxid > 1.0) oxid = 1.0;
double tA = fPox*oxid;  if(tA > 1.0) tA = 1.0;
dxdt_ATPM = KATPM*(tA - ATPM);
double nPc = hillf(PSER0, KPPCR, HPATP);
double tC = (hillf(Pser, KPPCR, HPATP)/nPc)*ATPM;  if(tC > 1.0) tC = 1.0;
dxdt_PCRM = KPCR*(tC - PCRM);
double nPd = hillf(PSER0, KPDPG, HPATP);
double tD = hillf(Pser, KPDPG, HPATP)/nPd;  if(tD > 1.0) tD = 1.0;
dxdt_DPG = KDPG*(tD - DPG);
dxdt_ATPD = KATPD*(tA - ATPD);

// ---- 11. SODIUM / VOLUME / ALDOSTERONE -----------------------------------------------
// Potassium is the dominant aldosterone secretagogue, so hypokalaemia brakes
// its own renal loss; volume depletion pushes the other way.
double ald_t = (0.15 + 0.85*SQ(Kser/KSER0))*(1.0 + 1.8*posv(1.0 - BW/BW0));
dxdt_ALD = KALD*(ald_t - ALD);
double Na_in = (PHASE > 1.5) ? NA_INTAKE/24.0 : NA_INTAKE*0.5/24.0;
double ret = RINS*posv(X - 1.0)/(1.0 + posv(X - 1.0))
           + RALD*posv(ALD - 1.0)/(1.0 + posv(ALD - 1.0));
if(ret > 0.92) ret = 0.92;
dxdt_NAB = Na_in*ret - KNAB*NAB;

// ---- 12. GAS EXCHANGE --------------------------------------------------------------
double choMax = CHOOXMAX*BW*60.0/180.0;
double choAct = (Rd < choMax) ? Rd : choMax;
double dnl = posv(Rd - choMax)*((GLY > 0.92*GLYMAX) ? 1.0 : 0.25);
double EE_h = (REE*PAL + 0.10*kcal_d)/24.0;
double e_cho = choAct*0.180*4.0;  if(e_cho > EE_h) e_cho = EE_h;
double e_fat = posv(EE_h - e_cho);
double RQ = (1.00*e_cho + 0.71*e_fat)/((e_cho + e_fat > 1e-9) ? (e_cho + e_fat) : 1e-9);
double VO2 = EE_h*24.0/(1440.0*KCALO2);
double VCO2 = RQ*VO2 + 2.44*dnl*22.4/1000.0/60.0;
// Alveolar ventilation cannot fall to zero in a living patient, and PaCO2 is
// bounded on both sides.  Leaving this ratio unbounded puts an unbounded
// argument inside the respiratory-failure exponential.
double VAcap = ATPD*sqrt(ffm/FFM0);  if(VAcap < VAMIN) VAcap = VAMIN;
double PaCO2 = PACO2_0*(VCO2/VCO2REF)/VAcap;
if(PaCO2 < 25.0) PaCO2 = 25.0;  if(PaCO2 > 120.0) PaCO2 = 120.0;

// ---- 13. CARDIAC FUNCTION AND HAZARDS -----------------------------------------------
double CONT = ATPM*sqrt(posv(PCRM));  if(CONT < CONTMIN) CONT = CONTMIN;
double lvr = LVM/LVM0;  if(lvr < 0.25) lvr = 0.25;
double cong = posv(Vecf/VECF0 - 1.0)/lvr/CONT;
double QTc = QTC0 + AK*posv(KSER0 - Kser) + AMG*posv(MGSER0 - MGser)
           + ACA*posv(CA0 - Caser);
double h_arr = HAZ*H0ARR*exp(CQT*posv(QTc - 440.0))*(1.0 + CATPARR*posv(1.0 - ATPM));
double h_cf  = HAZ*H0CF*pow(1.0 + cong, PCONG)/CONT;
double advd = (ATPD > VAMIN) ? ATPD : VAMIN;
double h_rf  = HAZ*H0RF*exp(QCO2*posv(PaCO2 - PACO2REF))/advd;
double rr2 = Rd/40.0;  if(rr2 > 2.0) rr2 = 2.0;
double h_we = H0WE*SQ(posv(1.0 - fTK))*(1.0 + 2.0*rr2);
double h_tot = h_arr + h_cf + h_rf;  if(h_tot > HAZMAX) h_tot = HAZMAX;

dxdt_H_ARR = h_arr;
dxdt_H_CF  = h_cf;
dxdt_H_RF  = h_rf;
dxdt_H_WE  = h_we*(1.0 - ((H_WE < 0.999) ? H_WE : 0.999));
dxdt_SURV  = -SURV*h_tot;
dxdt_AUCP  = posv(0.65 - Pser);
dxdt_CUMP  = ivP + JPabs;
dxdt_CUMK  = ivK + JKabs;
dxdt_CUMKCAL = kcal_h;

$TABLE
double BWo   = posv(FM) + posv(FFM);
double Vecfo = FECF*BW0 + NAB/NAREF;  if(Vecfo < 1.0) Vecfo = 1.0;
double Psero = posv(PE)/Vecfo;
double Ksero = posv(KE)/Vecfo;
double MGsero = posv(MGE)/Vecfo;
double Casero = posv(CAE)/Vecfo;
double thno  = posv(THT)/THT0;
double fTKo  = hillf(thno, KMTPP, HTPP)/FPDH0;  if(fTKo > 1.0) fTKo = 1.0;
double QTco  = QTC0 + AK*posv(KSER0 - Ksero) + AMG*posv(MGSER0 - MGsero)
             + ACA*posv(CA0 - Casero);

capture PSER  = Psero;
capture KSER  = Ksero;
capture MGSER = MGsero;
capture CASER = Casero;
capture CAXP  = Casero*Psero;
capture BWkg  = BWo;
capture BMI   = BWo/(HT*HT);
capture PICF  = posv(PI)/(RHOP*((FFM > 1.0) ? FFM : 1.0));
capture KICF  = posv(KI)/(RHOK*((FFM > 1.0) ? FFM : 1.0));
capture MGICF = posv(MGI)/(RHOMG*((FFM > 1.0) ? FFM : 1.0));
capture TKACT = fTKo;
capture QTC   = QTco;
capture MORT  = 100.0*(1.0 - SURV);
capture WERN  = 100.0*H_WE;
capture OEDEMA = posv(Vecfo - FECF*BW0);
capture TBP   = posv(PE) + posv(PI) + posv(PB);
'

rfs <- mcode_cache("rfs", rfs_code)

## ============================================================================
##  DRIVER:  starvation phase, then the refeeding protocol
## ============================================================================
##  The deficit is SIMULATED rather than asserted.  Phase 1 runs the patient
##  through their starvation history; phase 2 refeeds whatever state that
##  produced.  The model's first real test is that serum phosphate, potassium
##  and magnesium stay inside (or at the edge of) their reference intervals
##  throughout phase 1 while the stores empty - because if the blood tests
##  moved, admission bloods would be a useful triage test and this syndrome
##  would not be dangerous.

STATES <- c("A_CHO","A_FAT","A_PRO","GLU","INS","X","GCG","GLY","FFA","BHB",
            "LAC","PGUT","PE","PI","PB","KGUT","KE","KI","MGGUT","MGE","MGI",
            "CAE","CAB","PTH","CTD","FGF","THGUT","THP","THT","ATPM","PCRM",
            "DPG","ATPD","FM","FFM","AT","LVM","NAB","ALD","H_ARR","H_CF",
            "H_RF","H_WE","SURV","AUCP","CUMP","CUMK","CUMKCAL")

#' Run one virtual patient: starvation history, then a refeeding protocol.
rfs_run <- function(starve_days = 60, sim_days = 14, ...) {
  args <- list(...)
  starve_par <- args[intersect(names(args),
                               c("STARVE_FRAC","QUALITY","ALCOHOL",
                                 "DIURETIC","GI_LOSS"))]

  ## ---- phase 1: starvation ------------------------------------------------
  m1 <- rfs
  if (length(starve_par)) m1 <- do.call(param, c(list(m1), starve_par))
  m1 <- param(m1, PHASE = 1)
  o1 <- mrgsim(m1, end = starve_days * 24, delta = 6, hmax = 6)
  last <- as.data.frame(o1)[nrow(as.data.frame(o1)), ]

  ## clocks that only make sense during refeeding are reset at admission
  st <- as.list(last[STATES])
  for (k in c("H_ARR","H_CF","H_RF","H_WE","AUCP","CUMP","CUMK","CUMKCAL"))
    st[[k]] <- 0
  st[["SURV"]] <- 1

  ## ---- phase 2: refeeding -------------------------------------------------
  m2 <- do.call(param, c(list(rfs), args))
  m2 <- param(m2, PHASE = 2)
  m2 <- do.call(init, c(list(m2), st))
  o2 <- mrgsim(m2, end = sim_days * 24, delta = 0.25, hmax = 0.5)

  list(starve = as.data.frame(o1),
       refeed = as.data.frame(o2),
       admission = last)
}

#' Condense one run into the numbers a clinician would actually look at.
rfs_summary <- function(r) {
  d <- r$refeed
  at <- function(day, col) d[[col]][which.min(abs(d$time - day * 24))]
  data.frame(
    P_nadir   = min(d$PSER),
    P_nadir_d = d$time[which.min(d$PSER)] / 24,
    P_d1 = at(1, "PSER"), P_d3 = at(3, "PSER"), P_d7 = at(7, "PSER"),
    K_d3 = at(3, "KSER"), Mg_d3 = at(3, "MGSER"), Ca_d3 = at(3, "CASER"),
    QTc_max   = max(d$QTC),
    Lac_max   = max(d$LAC),
    TK_min    = min(d$TKACT),
    CaxP_max  = max(d$CAXP),
    oedema_L  = max(d$OEDEMA),
    mortality = tail(d$MORT, 1),
    wernicke  = tail(d$WERN, 1),
    kcal      = tail(d$CUMKCAL, 1),
    dFFM      = tail(d$FFM, 1) - d$FFM[1]
  )
}

## ============================================================================
##  SIXTEEN TREATMENT SCENARIOS
## ============================================================================
##  Every scenario shares the same simulated starvation history unless its
##  entry says otherwise, so each difference is caused by the REGIMEN.

HX <- list(STARVE_FRAC = 0.32, QUALITY = 0.25, GI_LOSS = 0.20)
hx <- function(...) modifyList(HX, list(...))

rfs_scenarios <- list(

  ## -- guideline comparators ------------------------------------------------
  "01 NICE 2006 cautious ramp" = hx(
    KCAL_START = 10, KCAL_GOAL = 30, ADV_FRAC = 0.20, ADV_START = 2,
    P_DOSE = 0.5, K_DOSE = 2.5, MG_DOSE = 0.3,
    TH_IV = 200, TH_LEAD = 0.5, TH_DAYS = 10, RESCUE = 1),

  "02 ASPEN 2020 faster ramp" = hx(
    KCAL_START = 15, KCAL_GOAL = 30, ADV_FRAC = 0.33, ADV_START = 1,
    P_DOSE = 0.6, K_DOSE = 3.0, MG_DOSE = 0.4,
    TH_IV = 100, TH_LEAD = 0.5, TH_DAYS = 7, RESCUE = 1),

  ## -- the 2x2: which prophylaxis is doing the work? ------------------------
  "03 full feed, no prophylaxis" = hx(
    KCAL_START = 30, KCAL_GOAL = 30, ADV_START = 99),

  "04 full feed + electrolytes only" = hx(
    KCAL_START = 30, KCAL_GOAL = 30, ADV_START = 99,
    P_DOSE = 0.8, K_DOSE = 3.0, MG_DOSE = 0.4, RESCUE = 1),

  "05 full feed + thiamine only" = hx(
    KCAL_START = 30, KCAL_GOAL = 30, ADV_START = 99,
    TH_IV = 300, TH_LEAD = 0.5),

  "06 full feed + both" = hx(
    KCAL_START = 30, KCAL_GOAL = 30, ADV_START = 99,
    P_DOSE = 0.8, K_DOSE = 3.0, MG_DOSE = 0.4,
    TH_IV = 300, TH_LEAD = 0.5, RESCUE = 1),

  ## -- does the calorie cap help on its own? --------------------------------
  "07 slow ramp, no prophylaxis" = hx(
    KCAL_START = 10, KCAL_GOAL = 30, ADV_FRAC = 0.20, ADV_START = 2),

  ## -- the feed IS a phosphate prescription ---------------------------------
  "08 IV dextrose only, 100 % CHO" = hx(
    KCAL_START = 20, KCAL_GOAL = 20, ADV_START = 99, ROUTE = 1,
    CHO_FRAC = 1.0, FAT_FRAC = 0.0, PRO_FRAC = 0.0),

  "09 enteral formula, 100 % CHO" = hx(
    KCAL_START = 20, KCAL_GOAL = 20, ADV_START = 99, ROUTE = 0,
    CHO_FRAC = 1.0, FAT_FRAC = 0.0, PRO_FRAC = 0.0),

  "10 enteral formula, 40 % CHO" = hx(
    KCAL_START = 20, KCAL_GOAL = 20, ADV_START = 99, ROUTE = 0,
    CHO_FRAC = 0.40, FAT_FRAC = 0.42, PRO_FRAC = 0.18),

  ## -- the alcohol phenotype: a thiamine disease, not a phosphate one -------
  "11 alcohol, no thiamine" = hx(
    ALCOHOL = 0.75, QUALITY = 0.10, STARVE_FRAC = 0.40, GI_LOSS = 0.05,
    KCAL_START = 25, KCAL_GOAL = 25, ADV_START = 99),

  "12 alcohol + IV thiamine" = hx(
    ALCOHOL = 0.75, QUALITY = 0.10, STARVE_FRAC = 0.40, GI_LOSS = 0.05,
    KCAL_START = 25, KCAL_GOAL = 25, ADV_START = 99,
    TH_IV = 500, TH_LEAD = 1.0),

  "13 alcohol + oral thiamine" = hx(
    ALCOHOL = 0.75, QUALITY = 0.10, STARVE_FRAC = 0.40, GI_LOSS = 0.05,
    KCAL_START = 25, KCAL_GOAL = 25, ADV_START = 99, TH_PO = 300),

  ## -- magnesium: how you give it, and what it unlocks ----------------------
  "14 Mg as a 2 h bolus" = hx(
    DIURETIC = 0.6, KCAL_START = 20, KCAL_GOAL = 20, ADV_START = 99,
    MG_DOSE = 0.4, MG_BOLUS = 1, K_DOSE = 2.0),

  "15 Mg as a 24 h infusion" = hx(
    DIURETIC = 0.6, KCAL_START = 20, KCAL_GOAL = 20, ADV_START = 99,
    MG_DOSE = 0.4, MG_BOLUS = 0, K_DOSE = 2.0),

  "16 potassium without magnesium" = hx(
    DIURETIC = 0.6, KCAL_START = 20, KCAL_GOAL = 20, ADV_START = 99,
    MG_DOSE = 0.0, K_DOSE = 2.0)
)

#' Run all sixteen and return one tidy table.
rfs_run_all <- function(starve_days = 60, sim_days = 14) {
  do.call(rbind, lapply(names(rfs_scenarios), function(nm) {
    r <- do.call(rfs_run, c(list(starve_days = starve_days,
                                 sim_days = sim_days),
                            rfs_scenarios[[nm]]))
    cbind(scenario = nm, rfs_summary(r), row.names = NULL)
  }))
}

## ============================================================================
##  PARAMETER SWEEPS  -  the results the model exists to produce
## ============================================================================

#' RESULT: the demand tracks GLUCOSE, not calories.  Energy is held fixed at
#' 25 kcal/kg/d and only the carbohydrate fraction moves.
rfs_sweep_gir <- function(cho = c(.2,.3,.4,.5,.6,.75,.9,1.0)) {
  do.call(rbind, lapply(cho, function(cf) {
    r <- do.call(rfs_run, c(list(starve_days = 60, sim_days = 14),
                            hx(KCAL_START = 25, KCAL_GOAL = 25, ADV_START = 99,
                               CHO_FRAC = cf,
                               FAT_FRAC = (1 - cf) * 0.64,
                               PRO_FRAC = (1 - cf) * 0.36)))
    cbind(CHO_pct = 100 * cf, rfs_summary(r), row.names = NULL)
  }))
}

#' RESULT: admission phosphate is nearly uninformative.  Six histories, one
#' regimen; the admission panel barely moves while the outcome does.
rfs_sweep_history <- function() {
  H <- list(c(30,.50), c(45,.40), c(60,.32), c(75,.28), c(90,.25), c(120,.22))
  do.call(rbind, lapply(H, function(h) {
    r <- do.call(rfs_run, c(list(starve_days = h[1], sim_days = 14),
                            hx(STARVE_FRAC = h[2], KCAL_START = 30,
                               KCAL_GOAL = 30, ADV_START = 99)))
    a <- r$admission
    cbind(days = h[1], frac = h[2],
          adm_P = a$PSER, adm_K = a$KSER, adm_ICF_K = a$KICF,
          rfs_summary(r), row.names = NULL)
  }))
}

#' RESULT: phosphate repletion - where the benefit saturates, and where
#' (far outside clinical practice) the calcium-phosphate product finally bites.
rfs_sweep_pdose <- function(doses = c(0,.15,.3,.5,.8,1.2,2,3,5,8,12)) {
  do.call(rbind, lapply(doses, function(pd) {
    r <- do.call(rfs_run, c(list(starve_days = 60, sim_days = 14),
                            hx(KCAL_START = 30, KCAL_GOAL = 30,
                               ADV_START = 99, P_DOSE = pd)))
    cbind(P_dose = pd, rfs_summary(r), row.names = NULL)
  }))
}

#' RESULT: the two clocks.  Thiamine timing on one axis, phosphate on the other.
rfs_sweep_timing <- function() {
  th <- do.call(rbind, lapply(list(c(1,0), c(0,0), c(0,6), c(0,24), c(0,72)),
    function(z) {
      r <- do.call(rfs_run, c(list(starve_days = 60, sim_days = 14),
                              hx(KCAL_START = 30, KCAL_GOAL = 30, ADV_START = 99,
                                 TH_IV = 300, TH_LEAD = z[1], TH_DELAY = z[2])))
      cbind(agent = "thiamine", lead_h = z[1], delay_h = z[2],
            rfs_summary(r), row.names = NULL)
    }))
  ph <- do.call(rbind, lapply(c(0, 1, 2, 4), function(pd) {
      r <- do.call(rfs_run, c(list(starve_days = 60, sim_days = 14),
                              hx(KCAL_START = 30, KCAL_GOAL = 30, ADV_START = 99,
                                 P_DOSE = 0.5, P_DELAY = pd)))
      cbind(agent = "phosphate", lead_h = 0, delay_h = pd * 24,
            rfs_summary(r), row.names = NULL)
    }))
  rbind(th, ph)
}

#' RESULT (a prediction that failed, reported rather than tuned away): oral
#' and intravenous thiamine CONVERGE within a day when the gut works.  The
#' intravenous route earns its place only when absorption is impaired.
rfs_sweep_thiamine <- function() {
  R <- list(
    list(lbl = "no thiamine",        p = list()),
    list(lbl = "oral 100 mg/d",      p = list(TH_PO = 100)),
    list(lbl = "oral 300 mg/d",      p = list(TH_PO = 300)),
    list(lbl = "oral 1500 mg/d",     p = list(TH_PO = 1500)),
    list(lbl = "IV 100 mg/d",        p = list(TH_IV = 100)),
    list(lbl = "IV 500 mg/d",        p = list(TH_IV = 500)),
    list(lbl = "oral 300, vomiting", p = list(TH_PO = 300, GUT_FAIL = 0.92)),
    list(lbl = "IV 100, vomiting",   p = list(TH_IV = 100, GUT_FAIL = 0.92)))
  do.call(rbind, lapply(R, function(z) {
    r <- do.call(rfs_run, c(list(starve_days = 60, sim_days = 14),
                            do.call(hx, c(list(ALCOHOL = 0.75, QUALITY = 0.10,
                                               STARVE_FRAC = 0.40, GI_LOSS = 0.05,
                                               KCAL_START = 25, KCAL_GOAL = 25,
                                               ADV_START = 99), z$p))))
    d <- r$refeed
    tk <- function(h) d$TKACT[which.min(abs(d$time - h))]
    data.frame(regimen = z$lbl, TK_3h = tk(3), TK_12h = tk(12),
               TK_24h = tk(24), TK_d3 = tk(72),
               wernicke = tail(d$WERN, 1), row.names = NULL)
  }))
}

## ============================================================================
##  Example use
## ============================================================================
if (identical(environment(), globalenv()) && !interactive()) {
  cat("\n== sixteen scenarios ==\n");  print(rfs_run_all(), digits = 3)
  cat("\n== glucose infusion rate sweep ==\n"); print(rfs_sweep_gir(), digits = 3)
  cat("\n== starvation history sweep ==\n");    print(rfs_sweep_history(), digits = 3)
  cat("\n== phosphate dose sweep ==\n");        print(rfs_sweep_pdose(), digits = 3)
  cat("\n== timing: the two clocks ==\n");      print(rfs_sweep_timing(), digits = 3)
  cat("\n== oral vs IV thiamine ==\n");         print(rfs_sweep_thiamine(), digits = 3)
}
