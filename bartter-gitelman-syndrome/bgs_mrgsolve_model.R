## =============================================================================
##  bgs_mrgsolve_model.R
##  QSP model of the SALT-LOSING TUBULOPATHIES:
##  Bartter syndrome (types I, II, III, IVa/IVb, V) and Gitelman syndrome
##
##  38 ODEs · nephron-segment NaCl/K/Cl/HCO3/Mg/Ca handling · macula densa
##  NaCl sensing · COX-2/PGE2 amplification · RAAS · drug PK/PD · growth ·
##  nephrocalcinosis · CKD progression · 18 therapeutic and diagnostic
##  scenarios.
##
## -----------------------------------------------------------------------------
##  THE MODELLING IDEA
## -----------------------------------------------------------------------------
##  Bartter and Gitelman syndrome are given exactly ONE structural difference:
##  WHERE the NaCl transport lesion sits relative to the macula densa.
##
##      FTAL  = fractional thick-ascending-limb NaCl transport capacity  [0,1]
##      FDCT  = fractional distal-convoluted-tubule NCC capacity         [0,1]
##
##  The macula densa is anatomically at the END of the TAL, so in the model it
##  reads a signal that depends on FTAL and NOT on FDCT — and it reads it
##  through its OWN apical NKCC2, which carries the same lesion. Everything
##  clinically diagnostic then falls out with no genotype-specific switches
##  anywhere in the code:
##
##    * urinary PGE2 is 5-10x elevated when the lesion is in the TAL and
##      NORMAL when it is in the DCT  -> "hyperprostaglandin E syndrome";
##    * COX inhibition is therefore disease-modifying in Bartter and nearly
##      inert in Gitelman — the SAME drug, the SAME parameters, opposite
##      answers, because the drug acts on a loop that only one genotype runs;
##    * urinary calcium CHANGES SIGN: TAL lesions collapse the lumen-positive
##      potential and cause hypercalciuria -> nephrocalcinosis -> CKD, while
##      DCT lesions leave the TAL intact and let volume contraction drive
##      proximal Ca reabsorption -> hypocalciuria;
##    * serum magnesium is lowest in Gitelman, near-normal in ANTENATAL
##      Bartter (the huge distal Mg load recruits a hypertrophied,
##      TRPM6-upregulated DCT) and intermediate in type III — the well-known
##      but counter-intuitive ordering, reproduced without being coded in;
##    * hypomagnesaemia releases the Mg-block of ROMK, so K+ cannot be
##      repleted until Mg2+ is — the model will not let you fix potassium
##      with potassium alone.
##
## -----------------------------------------------------------------------------
##  UNITS
## -----------------------------------------------------------------------------
##  time      days                     amounts (electrolyte)  mmol
##  volume    L                        concentration          mmol/L
##  drug amt  mg                       drug conc              mg/L
##  flux      mmol/day                 urine volume           L/day
##
## -----------------------------------------------------------------------------
##  CALIBRATION TARGETS (adult, 60 kg, eGFR 100, unless stated)
## -----------------------------------------------------------------------------
##  Phenotype        K      Mg     HCO3   uCa/Cr   uPGE2   PRA    Uvol
##                 mmol/L  mmol/L  mmol/L mmol/mmol  xULN   xULN   L/d
##  healthy         4.2    0.85    25     0.30-0.40  1.0    1.0    1.5
##  Gitelman        2.8-3.1 0.48-0.55 30-32 <0.07    1.0-1.3 2-3   1.9
##  Bartter III     2.7-3.1 0.50-0.62 31-34 0.25-0.45 3-4    4-5   3.0
##  Bartter I/II    2.5-3.0 0.72-0.85 30-34 0.7-1.1   6-9    8-10  >6
##  (Bettinelli 1992 J Pediatr; Peters 2002 Am J Med; Seyberth 2011 Pediatr
##   Nephrol; Blanchard 2017 Kidney Int KDIGO consensus; Konrad 2021 Kidney
##   Int consensus; Reinalter 2002 Kidney Int for the indomethacin/PGE2 data.)
##
## -----------------------------------------------------------------------------
##  IMPORTANT — the model is initialised at the HEALTHY steady state.
##  Setting FTAL/FDCT away from 1 starts a transient. Always burn in
##  (>= 300 days) before reading a phenotype, or use bgs_steady().
## =============================================================================

library(mrgsolve)
library(dplyr)

bgs_code <- '

$PROB
# Bartter & Gitelman syndrome QSP model (bgs)
- 38 ODEs; nephron-segment transport with positional macula-densa sensing.

$GLOBAL
// ---- body-size and reference constants (assigned in $MAIN, used in $ODE) ---
// These MUST live at file scope: $MAIN locals are not visible inside $ODE.
double FBSA, FBW, VIC, ECFR, FILTNAR, LMDR, LDCTR, LASDNR, LMGR;

// ---- numerically safe softplus, used for the bicarbonate threshold --------
inline double softplus(double x, double s) {
  double z = x / s;
  return (z > 30.0) ? x : s * log1p(exp(z));
}

// ---- captured intermediates (declared here so $CAPTURE can see them) -------
double CNA, CK, CCL, CHC, CMG, CCA, CKIC;
double GFRML, GFRLD, FILTNA, RNAPT, LMD, RNATAL, LDCT, RNADCT, LASDN, RNACD;
double UNA, UK, UCL, UHCO3, UCA, UMG, UVOL, UOSMA, SOLOUT;
double MDSENSE, LUMREL, ATAL, ANCC, FPGTAL, PGREL, ANGREL, ALDREL;
double ENACF, MRACT, INH2, INH1, AMIB, MRB, ACEI, THZ, FUR;
double CIND, CCEL, CAMI, CCAN, CACE;
double JH, JSHIFT, BICTH, PDREL, PDCA, FPTCA, FTALCA, DELCA, DELMG;
double FABSMG, MGABS, MGSTR, TRPMR, TRPVR, LMGREL;
double UCACR, FEMG, FEK, QTC, CRAMP, FATIG, VOLIDX, WIN, GFRPG;
double SATIDX, KSEC, ROMKF, EGFR, HTVEL, UPGE2, NAKRATIO, TTKG;
double ULCER, AKIFLAG, DELIVREL, ADHER;

$PARAM @annotated
// ---------------- patient / genotype ---------------------------------------
BW      :  60   : body weight (kg)
HTCM    : 168   : height (cm), for BSA
PEDS    :   0   : 1 = growing child (growth endpoint active)
GFR0    : 100   : baseline GFR (mL/min/1.73m2)
FTAL    :   1   : TAL NaCl transport capacity fraction [0-1]
FDCT    :   1   : DCT NCC transport capacity fraction [0-1]
FROMKCD :   1   : collecting-duct ROMK capacity (0.30 in Bartter type II)
FCLDN   :   1   : TAL paracellular claudin-16/19 function
CASRGF  :   0   : CaSR gain-of-function activity (ADH type 1), 0-1

// ---------------- dietary intake -------------------------------------------
NAIN    : 150   : dietary Na intake (mmol/day)
KIN     :  70   : dietary K intake (mmol/day)
MGIN    :  12   : dietary Mg intake (mmol/day)
CAIN    :   6   : net intestinal Ca absorption (mmol/day)
FKABS   : 0.90  : fractional K absorption
FCLDIET : 0.15  : fraction of dietary K accompanied by chloride

// ---------------- proximal tubule ------------------------------------------
FPT0    : 0.65  : baseline PT fractional Na reabsorption
EPTANG  : 0.060 : Ang II effect on PT Na reabsorption
FPTMIN  : 0.45  : lower bound on PT fraction
FPTMAX  : 0.80  : upper bound on PT fraction

// ---------------- thick ascending limb -------------------------------------
FTALF0  : 0.80  : fraction of delivered Na reabsorbed by a healthy TAL
IPGTAL  : 2.50  : PGE2 excess giving 50% NKCC2 inhibition (xULN above 1)
ECASR   : 0.55  : CaSR gain-of-function inhibition of NKCC2/ROMK

// ---------------- macula densa / COX-2 / PGE2 ------------------------------
EMD     : 0.20  : exponent, luminal NaCl -> sensed signal
ECOX    : 7.00  : max fold COX-2 induction by loss of sensed NaCl
KOCOX   : 1.20  : COX-2 turnover (1/day)
KOPG    : 6.00  : renal PGE2 turnover (1/day)

// ---------------- distal convoluted tubule ---------------------------------
FDCTF0  : 0.55  : fraction of delivered Na reabsorbed by a healthy DCT
ENCC    : 0.35  : compensatory NCC up-phosphorylation by distal load
KNCC    : 0.20  : NCC adaptation rate (1/day)
KDCTM   : 0.020 : DCT mass remodelling rate (1/day)
EDCTHYP : 0.20  : DCT hypertrophy gain from luminal load

// ---------------- aldosterone-sensitive distal nephron ---------------------
VMCD    : 970   : max ASDN Na reabsorption (mmol/day, 60 kg)
KMCD    : 400   : ASDN Na reabsorption Km (mmol/day)
KENAC   : 0.35  : ENaC abundance adaptation rate (1/day)
EENACA  : 0.70  : MR-activity exponent on ENaC abundance
KS0     :  63   : baseline distal K secretion (mmol/day)
EKENAC  : 0.60  : ENaC exponent on K secretion
EKFLOW  : 0.80  : distal-flow exponent on K secretion
KALKK   : 0.55  : alkalosis effect on K secretion
KMGROMK : 1.20  : hypomagnesaemia release of ROMK Mg-block (KEY)

// ---------------- acid-base -------------------------------------------------
JH0     : 82.3  : baseline distal H+ secretion (mmol/day)
NEAP    :  60   : net endogenous acid production (mmol/day)
EJHENAC : 0.35  : ENaC exponent on H+ secretion
EJHFLOW : 0.20  : distal-flow exponent on H+ secretion
EJHK    : 0.35  : hypokalaemia exponent on H+ secretion
JHMAX   : 2.20  : ceiling on distal H+ secretion (fold)
KBIC    : 3.00  : bicarbonaturia gain (1/day)
BICTH0  :  26   : bicarbonate threshold at normal chloride (mmol/L)
KCLTH   : 0.55  : chloride-depletion elevation of the HCO3 threshold
UOA0    : 74.5  : unmeasured urinary anions (mmol/day)
FNH4    : 0.55  : fraction of distal H+ excreted as NH4+

// ---------------- potassium distribution -----------------------------------
KSH     : 8.00  : transcellular K exchange rate (1/day)
KALKSH  : 0.90  : alkalosis effect on the IC/EC K ratio
KICANCH : 0.20  : intracellular K anchoring rate (1/day)

// ---------------- magnesium -------------------------------------------------
FMGUF   : 0.70  : ultrafilterable fraction of plasma Mg
FMGPT   : 0.15  : PT fractional Mg reabsorption
FMGTAL  : 0.70  : TAL fractional (paracellular) Mg reabsorption
PDMIN   : 0.35  : PD-independent component of TAL paracellular transport
KMGCASR : 0.10  : CaSR-mediated TAL Mg reabsorption boost when Mg is low
FDCTMG0 : 0.667 : fraction of DELIVERED Mg reabsorbed by a healthy DCT
KMGLOAD : 0.20  : load-recruitment of DCT Mg reabsorption (KEY for Bartter I)
ETRPM6  : 0.80  : TRPM6 upregulation by hypomagnesaemia
KTRPM   : 0.05  : TRPM6 adaptation rate (1/day)
FABS0   : 0.35  : basal fractional intestinal Mg absorption
FABSMX  : 0.75  : maximal fractional intestinal Mg absorption
KMGSAT  :  30   : oral Mg dose giving half-maximal absorption saturation (mmol)
KMS     : 1.20  : bone/intracellular Mg exchange (mmol/day)
MGS0    : 800   : exchangeable Mg store (mmol)

// ---------------- calcium ---------------------------------------------------
FCAPT0  : 0.65  : baseline PT fractional Ca reabsorption
ECAVOL  : 0.22  : volume-contraction boost to PT Ca reabsorption
ECAANG  : 0.10  : Ang II boost to PT Ca reabsorption
FPTCAMX : 0.76  : ceiling on PT fractional Ca reabsorption
FCATAL  : 0.22  : TAL fractional (paracellular) Ca reabsorption
PDMINCA : 0.12  : PD-independent component of TAL paracellular Ca transport
FDCTCA  : 0.868 : fraction of DELIVERED Ca reabsorbed by DCT2/CNT
TMCA    :  29   : Tmax of DCT2/CNT transcellular Ca reabsorption (mmol/day)
ETRPV5  : 0.45  : TRPV5 upregulation by volume contraction
KTRPV   : 0.05  : TRPV5 adaptation rate (1/day)
KCABONE : 8.00  : strength of the serum-Ca homeostat (1/day)
CACA0   : 1.46  : ultrafilterable Ca set-point (mmol/L)

// ---------------- RAAS ------------------------------------------------------
KOREN   :  40   : renin turnover (1/day)
GMD     : 0.90  : macula-densa exponent on renin
GVOL    : 8.00  : volume exponent on renin
GPG     : 0.35  : PGE2 (EP4) gain on renin
SRENMAX :  12   : soft ceiling on renin stimulation
KOANG   :  60   : Ang II turnover (1/day)
ANG0    :  10   : baseline Ang II (pmol/L)
KOALD   :  12   : aldosterone turnover (1/day)
ALD0    : 300   : baseline aldosterone (pmol/L)
EAA     : 0.90  : Ang II gain on aldosterone
EKALD   : 0.80  : plasma-K exponent on aldosterone

// ---------------- water / vasopressin ---------------------------------------
UOSMB   : 600   : achievable urine osmolality at baseline (mosm/kg)
ETALCON : 0.55  : TAL contribution to the medullary gradient
EPGAQP  : 0.30  : PGE2 (EP3) antagonism of AQP2
KAVP    : 4.00  : AVP adaptation rate (1/day)
WINS    : 0.80  : insensible water loss (L/day, 60 kg)
WMET    : 0.35  : metabolic water (L/day)
WMAXK   : 0.17  : maximum achievable water intake (L/day per kg)
KTHIRST : 3.00  : tonicity gain on thirst
UREAOSM : 400   : urea + non-electrolyte urinary solute (mosm/day)

// ---------------- structural progression ------------------------------------
KNC     : 0.0016: nephrocalcinosis accrual rate (1/day)
KNCOUT  : 0.0002: nephrocalcinosis regression rate (1/day)
SSTHR   : 0.65  : supersaturation index threshold for crystal accrual
UCAREF  : 5.00  : reference urinary Ca (mmol/day)
UVREF   : 1.50  : reference urine volume (L/day)
KFIB    : 4.0e-5: nephrocalcinosis-driven GFR loss (1/day)
KHYPOK  : 2.0e-5: hypokalaemic nephropathy GFR loss (1/day)
KNSAIDK : 3.0e-5: NSAID-driven chronic GFR loss (1/day)
KGPG    : 0.35  : max acute PG-dependent GFR loss

// ---------------- growth ----------------------------------------------------
KIGF    : 0.05  : IGF-1 adaptation rate (1/day)
KHT     : 1.20  : height-SDS gain (SDS/year at IGF1 = 1)
IGFTH   : 0.60  : IGF-1 level giving zero SDS drift
EGHK    : 0.80  : hypokalaemia exponent on the GH/IGF-1 axis
EGHALK  : 0.70  : alkalosis effect on the GH/IGF-1 axis
EGHPG   : 0.30  : PGE2 effect on the GH/IGF-1 axis
EGHVOL  : 0.50  : volume-contraction effect on the GH/IGF-1 axis
RHGH    :   0   : recombinant GH effect (0 = off, 0.6 = typical dose)

// ---------------- drug PK ---------------------------------------------------
KAIND   : 17.0  : indomethacin absorption (1/day)
CLIND   : 55.0  : indomethacin CL/F (L/day, 60 kg)
VIND    : 17.4  : indomethacin V/F (L)
FIND    : 0.98  : indomethacin bioavailability
KACEL   :  8.0  : celecoxib absorption (1/day)
CLCEL   :  665  : celecoxib CL/F (L/day)
VCEL    :  400  : celecoxib V/F (L)
FCEL    : 1.00  : celecoxib relative bioavailability (folded into CL/F)
KAAMI   :  6.0  : amiloride absorption (1/day)
CLAMI   :  624  : amiloride CL/F (L/day)
VAMI    :  350  : amiloride V/F (L)
KASPI   :  8.0  : spironolactone absorption (1/day)
FCAN    : 0.60  : spironolactone -> canrenone conversion
CLCAN   :  104  : canrenone CL/F (L/day)
VCAN    :  100  : canrenone V/F (L)
KAACE   :  6.0  : enalapril absorption (1/day)
FACE    : 0.40  : enalapril -> enalaprilat conversion
CLACE   :  110  : enalaprilat CL/F (L/day)
VACE    :   60  : enalaprilat V/F (L)
KAKCL   :  8.0  : oral KCl absorption (1/day)
KAMGO   :  5.0  : oral Mg absorption (1/day)
KANAC   : 10.0  : oral NaCl absorption (1/day)

// ---------------- drug PD ---------------------------------------------------
IC50IND : 0.30  : indomethacin IC50 for renal COX-2/PGE2 (mg/L)
EMXIND  : 0.95  : indomethacin Emax on renal PGE2
IC50IN1 : 0.15  : indomethacin IC50 for COX-1 (gastric) (mg/L)
IC50CEL : 0.30  : celecoxib IC50 for renal COX-2/PGE2 (mg/L)
EMXCEL  : 0.90  : celecoxib Emax on renal PGE2
IC50CE1 : 8.00  : celecoxib IC50 for COX-1 (mg/L)
IC50AMI : 0.012 : amiloride IC50 for ENaC (mg/L)
EMXAMI  : 0.85  : amiloride Emax on ENaC
IC50CAN : 0.10  : canrenone IC50 for MR (mg/L)
EMXCAN  : 0.80  : canrenone Emax on MR
IC50ACE : 0.010 : enalaprilat IC50 for ACE (mg/L)
EMXACE  : 0.85  : enalaprilat Emax on ACE
KGIREP  : 0.08  : gastric mucosal repair rate (1/day)
KGIDAM  : 0.30  : COX-1-driven gastric mucosal damage rate (1/day)

// ---------------- diagnostic challenges & stressors -------------------------
DXHCTZ  :   0   : thiazide challenge, fractional NCC block (0-0.95)
DXFURO  :   0   : furosemide challenge, fractional NKCC2 block (0-0.95)
DXT0    : 1e9   : challenge start time (day)
DXDUR   : 0.25  : challenge duration (day)
WEXTRA  :   0   : extrarenal water loss, e.g. gastroenteritis (L/day)
WEXT0   : 1e9   : extrarenal loss start time (day)
WEXTDUR :   3   : extrarenal loss duration (day)
ADHER0  : 1.00  : adherence multiplier applied to absorbed drug/supplement

// ---------------- endpoints -------------------------------------------------
QTCB    : 400   : baseline QTc (ms)
AKQT    :  24   : QTc prolongation per mmol/L of K below 4.0 (ms)
AMQT    :   7   : QTc prolongation per 0.1 mmol/L of Mg below 0.85 (ms)
UCRE    :  13   : urinary creatinine excretion (mmol/day)

$CMT @annotated
// --- drug PK ---------------------------------------------------------------
GIND : indomethacin gut depot (mg)
CIN  : indomethacin central (mg)
GCEL : celecoxib gut depot (mg)
CCE  : celecoxib central (mg)
GAMI : amiloride gut depot (mg)
CAM  : amiloride central (mg)
GSPI : spironolactone gut depot (mg)
SCAN : canrenone central (mg)
GACE : enalapril gut depot (mg)
CAC  : enalaprilat central (mg)
GKCL : oral KCl gut depot (mmol)
GMGO : oral Mg gut depot (mmol)
GNAC : oral NaCl gut depot (mmol)
// --- body fluid and electrolytes -------------------------------------------
ECFV : extracellular fluid volume (L)
ANA  : ECF sodium content (mmol)
AKE  : ECF potassium content (mmol)
AKI  : intracellular potassium pool (mmol)
ACL  : ECF chloride content (mmol)
AHC  : ECF bicarbonate content (mmol)
AMGE : ECF magnesium content (mmol)
AMGS : exchangeable Mg store, bone/intracellular (mmol)
ACAE : ECF ultrafilterable calcium content (mmol)
// --- signalling -------------------------------------------------------------
PRA  : plasma renin activity (ng/mL/h)
ANGT : angiotensin II (pmol/L)
ALDO : aldosterone (pmol/L)
COX2 : renal COX-2 expression (relative)
PGE  : renal PGE2 (relative)
AVP  : vasopressin (relative)
// --- tubular adaptation -----------------------------------------------------
NCCP : NCC phosphorylation / abundance (relative)
ENAC : ENaC abundance (relative)
TRPV : TRPV5 abundance (relative)
TRPM : TRPM6 abundance (relative)
DCTM : DCT mass / integrity (relative)
// --- long-term outcome ------------------------------------------------------
NEPH : nephrocalcinosis burden (ultrasound grade 0-3)
GFRF : functional GFR fraction (relative)
IGF1 : IGF-1 (relative)
HTSD : height SDS
GIM  : gastric mucosal integrity (0-1)

$MAIN
// ---------------------------------------------------------------------------
// scaling: renal fluxes by BSA, intakes and pools by body weight
// ---------------------------------------------------------------------------
double BSA = sqrt(HTCM * BW / 3600.0);
FBSA  = BSA / 1.73;
FBW   = BW / 60.0;
VIC   = 0.40 * BW;                        // intracellular water (L)
ECFR  = 0.20 * BW;                        // reference ECF volume (L)

// healthy reference luminal loads (derived, so the model self-scales)
FILTNAR = GFR0 * 1.44 * FBSA * 140.0;
LMDR    = FILTNAR * (1.0 - FPT0);
LDCTR   = LMDR  * (1.0 - FTALF0);
LASDNR  = LDCTR * (1.0 - FDCTF0);
LMGR    = GFR0 * 1.44 * FBSA * 0.85 * FMGUF * (1.0 - FMGPT - FMGTAL);

// ---------------------------------------------------------------------------
// initial conditions = healthy steady state
// ---------------------------------------------------------------------------
ECFV_0 = ECFR;
ANA_0  = 140.0 * ECFR;
AKE_0  = 4.2   * ECFR;
AKI_0  = 122.9 * VIC;
ACL_0  = 103.0 * ECFR;
AHC_0  = 25.0  * ECFR;
AMGE_0 = 0.85  * ECFR;
AMGS_0 = MGS0 * FBW;
ACAE_0 = CACA0 * ECFR;
PRA_0  = 1.0;   ANGT_0 = ANG0;  ALDO_0 = ALD0;
COX2_0 = 1.0;   PGE_0  = 1.0;   AVP_0  = 1.0;
NCCP_0 = 1.0;   ENAC_0 = 1.0;   TRPV_0 = 1.0;  TRPM_0 = 1.0;  DCTM_0 = 1.0;
NEPH_0 = 0.0;   GFRF_0 = 1.0;   IGF1_0 = 1.0;  HTSD_0 = 0.0;  GIM_0  = 1.0;

$ODE
// ===========================================================================
// 0. concentrations, adherence, active challenges
// ===========================================================================
double VOL = fmax(ECFV, 0.30 * ECFR);
CNA  = ANA  / VOL;
CK   = fmax(AKE  / VOL, 0.5);
CCL  = ACL  / VOL;
CHC  = AHC  / VOL;
CMG  = fmax(AMGE / VOL, 0.10);
CCA  = ACAE / VOL;
CKIC = AKI / VIC;
ADHER = ADHER0;

THZ = ((SOLVERTIME >= DXT0) && (SOLVERTIME < DXT0 + DXDUR)) ? DXHCTZ : 0.0;
FUR = ((SOLVERTIME >= DXT0) && (SOLVERTIME < DXT0 + DXDUR)) ? DXFURO : 0.0;
double WEX = ((SOLVERTIME >= WEXT0) && (SOLVERTIME < WEXT0 + WEXTDUR)) ? WEXTRA : 0.0;

// ===========================================================================
// 1. drug concentrations and target engagement
// ===========================================================================
CIND = CIN / VIND;
CCEL = CCE / VCEL;
CAMI = CAM / VAMI;
CCAN = SCAN / VCAN;
CACE = CAC / VACE;

double IIND = EMXIND * CIND / (IC50IND + CIND);
double ICEL = EMXCEL * CCEL / (IC50CEL + CCEL);
INH2 = 1.0 - (1.0 - IIND) * (1.0 - ICEL);          // renal COX-2 inhibition
INH1 = 1.0 - (1.0 - CIND / (IC50IN1 + CIND)) *
             (1.0 - CCEL / (IC50CE1 + CCEL));      // COX-1 (gastric)
AMIB = EMXAMI * CAMI / (IC50AMI + CAMI);
MRB  = EMXCAN * CCAN / (IC50CAN + CCAN);
ACEI = EMXACE * CACE / (IC50ACE + CACE);

// ===========================================================================
// 2. glomerular filtration (chronic loss + acute PG-dependent component)
// ===========================================================================
VOLIDX = fmin(1.5, fmax(0.0, (1.0 - VOL / ECFR) / 0.10));
PGREL  = PGE;
GFRPG  = 1.0 - KGPG * fmax(0.0, 1.0 - PGREL) * fmin(1.0, VOLIDX);
GFRML  = GFR0 * GFRF * GFRPG;
GFRLD  = GFRML * 1.44 * FBSA;
EGFR   = GFRML;

// ===========================================================================
// 3. nephron-segment sodium handling
// ===========================================================================
FILTNA = GFRLD * CNA;
ANGREL = ANGT / ANG0;
ALDREL = ALDO / ALD0;

double fpt = FPT0 * (1.0 + EPTANG * (ANGREL - 1.0));
fpt   = fmin(FPTMAX, fmax(FPTMIN, fpt));
RNAPT = FILTNA * fpt;
LMD   = FILTNA - RNAPT;                     // delivered to the TAL

// --- TAL: capacity x PGE2 inhibition x CaSR brake x acute furosemide -------
FPGTAL = 1.0 / (1.0 + fmax(0.0, PGREL - 1.0) / IPGTAL);
double casr = 1.0 - ECASR * CASRGF;
ATAL   = fmax(0.0, FTAL * FPGTAL * casr * (1.0 - FUR));
double ftalf = fmin(0.95, FTALF0 * ATAL);
RNATAL = LMD * ftalf;
LDCT   = LMD - RNATAL;                      // luminal load AT the macula densa

// --- THE POSITIONAL SENSOR -------------------------------------------------
// The macula densa sits at the end of the TAL and reads luminal NaCl through
// its OWN NKCC2, which carries the same lesion. A TAL lesion therefore blinds
// it despite a HIGH luminal NaCl (exactly as furosemide does). A DCT lesion
// is downstream and leaves the signal untouched.
LUMREL  = LDCT / LDCTR;
MDSENSE = fmin(2.0, ATAL * pow(fmax(LUMREL, 0.05), EMD));

// --- DCT --------------------------------------------------------------------
ANCC   = fmax(0.0, FDCT * NCCP * (1.0 - THZ));
double fdctf = fmin(0.95, FDCTF0 * ANCC);
RNADCT = LDCT * fdctf;
LASDN  = LDCT - RNADCT;                     // distal delivery: the coupler
DELIVREL = LASDN / LASDNR;

// --- ASDN -------------------------------------------------------------------
MRACT = ALDREL * (1.0 - MRB);
ENACF = fmax(0.0, ENAC * (1.0 - AMIB));
RNACD = VMCD * FBW * ENACF * LASDN / (KMCD * FBW + LASDN);
RNACD = fmin(RNACD, 0.97 * LASDN);
UNA   = LASDN - RNACD;

// ===========================================================================
// 4. potassium: distal secretion, the Mg-ROMK brake, transcellular shift
// ===========================================================================
// Hypomagnesaemia releases the intracellular Mg2+ block of ROMK, so K+ leaks
// out of the principal cell faster. This is why potassium cannot be repleted
// until magnesium is.
ROMKF = FROMKCD * (1.0 + KMGROMK * fmax(0.0, 1.0 - CMG / 0.85));
KSEC  = KS0 * FBW
        * pow(fmax(ENACF, 0.02), EKENAC)
        * pow(fmax(DELIVREL, 0.05), EKFLOW)
        * ROMKF
        * (1.0 + KALKK * (CHC / 25.0 - 1.0))
        * pow(CK / 4.2, 0.5);
UK = fmax(0.05 * KS0 * FBW, KSEC);

double ratio = (122.9 / 4.2) * (1.0 + KALKSH * (CHC / 25.0 - 1.0));
double ckset = CKIC / fmax(ratio, 5.0);
JSHIFT = KSH * VOL * (CK - ckset);                       // ECF -> cells
double akit = 122.9 * VIC * (0.55 + 0.45 * pow(CK / 4.2, 0.6));

// ===========================================================================
// 5. acid-base
// ===========================================================================
JH = JH0 * FBW * fmin(JHMAX,
        pow(fmax(ENACF, 0.02), EJHENAC)
      * pow(fmax(DELIVREL, 0.05), EJHFLOW)
      * pow(CK / 4.2, -EJHK));
// Chloride depletion RAISES the renal bicarbonate threshold — this is why
// alkalosis only corrects when the CHLORIDE salt is given.
BICTH = BICTH0 + KCLTH * 100.0 * fmax(0.0, 1.0 - CCL / 103.0);
UHCO3 = KBIC * VOL * softplus(CHC - BICTH, 1.5);

// urinary electroneutrality closes the chloride balance
UCL = fmax(5.0, UNA + UK + FNH4 * JH - UHCO3 - UOA0 * FBW);

// ===========================================================================
// 6. magnesium — the ordering Gitelman < Bartter III < Bartter I emerges here
// ===========================================================================
double mgdef = fmax(0.0, 1.0 - CMG / 0.85);
double filtmg = GFRLD * CMG * FMGUF;
PDREL  = PDMIN + (1.0 - PDMIN) * sqrt(fmax(ATAL, 0.0));      // lumen-+ PD
double fmgpt  = FMGPT * (1.0 + 0.05 * VOLIDX);
double fmgtal = FMGTAL * PDREL * FCLDN * (1.0 + KMGCASR * mgdef);
fmgtal = fmin(fmgtal, 0.90 - fmgpt);
DELMG  = filtmg * (1.0 - fmgpt - fmgtal);
LMGREL = DELMG / fmax(LMGR, 1e-6);
TRPMR  = TRPM * (1.0 + ETRPM6 * mgdef);
double fdctmg = fmin(0.95, FDCTMG0 * TRPMR *
                     (1.0 + KMGLOAD * sqrt(fmax(0.0, LMGREL - 1.0))));
UMG    = fmax(0.02 * filtmg, DELMG * (1.0 - fdctmg));
FEMG   = 100.0 * UMG / fmax(filtmg, 1e-6);

// adaptive intestinal absorption + saturable oral supplementation
FABSMG = FABS0 + (FABSMX - FABS0) * pow(mgdef, 0.7);
MGABS  = MGIN * FBW * FABSMG;               // dietary Mg absorbed (mmol/day)
MGSTR  = KMS * FBW * (AMGS / (MGS0 * FBW)) * (1.0 - CMG / 0.85);

// ===========================================================================
// 7. calcium — where the SIGN of urinary Ca is decided
// ===========================================================================
double filtca = GFRLD * CCA;
FPTCA = fmin(FPTCAMX, FCAPT0 * (1.0 + ECAVOL * VOLIDX + ECAANG * (ANGREL - 1.0)));
PDCA  = PDMINCA + (1.0 - PDMINCA) * pow(fmax(ATAL, 0.0), 0.45);
FTALCA = FCATAL * PDCA * FCLDN;
DELCA  = fmax(0.0, filtca * (1.0 - FPTCA - FTALCA));
TRPVR  = TRPV * (1.0 + ETRPV5 * VOLIDX);
double rcadct = fmin(TMCA * FBSA * (0.80 + 0.20 * TRPVR), FDCTCA * DELCA);
UCA    = fmax(0.05, DELCA - rcadct);
UCACR  = UCA / (UCRE * FBW);

// ===========================================================================
// 8. water, urine volume and concentrating ability
// ===========================================================================
SOLOUT = 2.0 * (UNA + UK) + UOA0 * FBW + UREAOSM * FBW;
UOSMA  = UOSMB * (1.0 - ETALCON * (1.0 - FTAL))
         / (1.0 + EPGAQP * fmax(0.0, PGREL - 1.0))
         * fmin(1.8, fmax(0.5, pow(AVP, 0.35)));
UVOL   = fmin(25.0, fmax(0.30 * FBW, SOLOUT / fmax(UOSMA, 60.0)));
WIN    = fmin(WMAXK * BW,
              (UVOL + WINS * FBW - WMET * FBW)
              * (1.0 + KTHIRST * fmax(0.0, (CNA - 140.0) / 5.0))
              + 1.2 * VOLIDX);

// ===========================================================================
// 9. renin -> Ang II -> aldosterone
// ===========================================================================
double smd  = pow(1.0 / (MDSENSE + 0.05), GMD) / pow(1.0 / 1.05, GMD);
double svol = exp(GVOL * fmax(-0.05, 1.0 - VOL / ECFR));
double spg  = 1.0 + GPG * fmax(0.0, PGREL - 1.0);
double sraw = smd * svol * spg;
double sren = sraw / (1.0 + sraw / SRENMAX) * (1.0 + 1.0 / SRENMAX);

// ===========================================================================
// 10. structural progression, growth, safety
// ===========================================================================
SATIDX = (UCA / (UCAREF * FBW)) * pow(UVREF * FBW / fmax(UVOL, 0.2), 0.30)
         * (1.0 + 0.40 * fmax(0.0, CHC / 25.0 - 1.0));
double ghf = pow(fmin(1.0, CK / 3.6), EGHK)
             / (1.0 + EGHALK * fmax(0.0, CHC / 28.0 - 1.0))
             / (1.0 + EGHPG * fmax(0.0, PGREL - 1.0))
             * fmax(0.20, 1.0 - EGHVOL * fmin(1.0, VOLIDX))
             * (1.0 + RHGH);
HTVEL  = (PEDS > 0.5) ? KHT * (IGF1 - IGFTH) / (1.0 - IGFTH) : 0.0;
ULCER  = 1.0 - GIM;

// ===========================================================================
// 11. clinical endpoint read-outs
// ===========================================================================
QTC   = QTCB + AKQT * fmax(0.0, 4.0 - CK) + AMQT * 10.0 * fmax(0.0, 0.85 - CMG);
CRAMP = fmin(10.0, 6.0 * pow(fmax(0.0, 1.0 - CMG / 0.85), 0.7)
                 + 4.0 * pow(fmax(0.0, 1.0 - CK / 4.0), 0.7));
FATIG = fmin(10.0, 7.0 * pow(fmax(0.0, 1.0 - CK / 4.0), 0.7)
                 + 3.0 * fmin(1.0, VOLIDX));
FEK   = 100.0 * UK / fmax(GFRLD * CK, 1e-6);
TTKG  = (CK > 0.1) ? (UK / fmax(UVOL, 0.2)) / CK * (290.0 / fmax(UOSMA, 100.0)) : 0.0;
UPGE2 = PGREL;
NAKRATIO = UNA / fmax(UK, 1.0);
AKIFLAG = (GFRPG < 0.80) ? 1.0 : 0.0;

// ===========================================================================
// 12. DIFFERENTIAL EQUATIONS
// ===========================================================================
// --- PK ---------------------------------------------------------------------
dxdt_GIND = -KAIND * GIND;
dxdt_CIN  =  KAIND * GIND * FIND * ADHER - (CLIND * FBW / VIND) * CIN;
dxdt_GCEL = -KACEL * GCEL;
dxdt_CCE  =  KACEL * GCEL * FCEL * ADHER - (CLCEL / VCEL) * CCE;
dxdt_GAMI = -KAAMI * GAMI;
dxdt_CAM  =  KAAMI * GAMI * ADHER - (CLAMI / VAMI) * CAM;
dxdt_GSPI = -KASPI * GSPI;
dxdt_SCAN =  KASPI * GSPI * FCAN * ADHER - (CLCAN / VCAN) * SCAN;
dxdt_GACE = -KAACE * GACE;
dxdt_CAC  =  KAACE * GACE * FACE * ADHER - (CLACE / VACE) * CAC;
dxdt_GKCL = -KAKCL * GKCL;
dxdt_GMGO = -KAMGO * GMGO;
dxdt_GNAC = -KANAC * GNAC;

// oral supplement absorption (Mg absorption is saturable AND adaptive)
double kclabs = KAKCL * GKCL * FKABS * ADHER;
double mgrate = KAMGO * GMGO;
double fmgeff = FABSMG / (1.0 + mgrate / (KMGSAT * FBW));
double mgabsr = mgrate * fmgeff * ADHER;
double mgunab = mgrate - mgabsr;
double gidiar = fmax(0.0, mgunab - 15.0 * FBW);   // osmotic diarrhoea
double nacabs = KANAC * GNAC * ADHER;

// --- water and electrolytes -------------------------------------------------
dxdt_ECFV = WIN + WMET * FBW - UVOL - WINS * FBW - WEX;
dxdt_ANA  = NAIN * FBW + nacabs - UNA;
dxdt_AKE  = KIN * FBW * FKABS + kclabs - UK - JSHIFT - 0.25 * gidiar;
dxdt_AKI  = JSHIFT + KICANCH * (akit - AKI);
dxdt_ACL  = NAIN * FBW + nacabs + FCLDIET * KIN * FBW + kclabs - UCL;
dxdt_AHC  = JH - NEAP * FBW - UHCO3;
dxdt_AMGE = MGABS + mgabsr + MGSTR - UMG - 0.25 * gidiar;
dxdt_AMGS = -MGSTR;
dxdt_ACAE = CAIN * FBW + KCABONE * VOL * (CACA0 - CCA) - UCA;

// --- signalling --------------------------------------------------------------
dxdt_PRA  = KOREN * (sren - PRA);
dxdt_ANGT = KOANG * (ANG0 * PRA * (1.0 - ACEI) - ANGT);
dxdt_ALDO = KOALD * (ALD0 * (1.0 + EAA * (ANGREL - 1.0))
                     * pow(CK / 4.2, EKALD) - ALDO);
dxdt_COX2 = KOCOX * (1.0 + ECOX * fmax(0.0, 1.0 - MDSENSE) - COX2);
dxdt_PGE  = KOPG  * (COX2 * (1.0 - INH2) - PGE);
dxdt_AVP  = KAVP  * ((1.0 + 4.0 * fmax(0.0, (CNA - 140.0) / 5.0))
                     * (1.0 + 2.5 * fmin(1.0, VOLIDX)) - AVP);

// --- tubular adaptation -------------------------------------------------------
dxdt_NCCP = KNCC * (1.0 + ENCC * pow(fmax(0.0, LUMREL - 1.0), 0.5)
                    * pow(fmax(ALDREL, 0.1), 0.2) - NCCP);
dxdt_ENAC = KENAC * (pow(fmax(MRACT, 0.05), EENACA)
                     * (1.0 + 0.15 * (ANGREL - 1.0)) - ENAC);
dxdt_TRPV = KTRPV * (pow(fmax(DCTM, 0.05), 0.5) - TRPV);
dxdt_TRPM = KTRPM * (pow(fmax(DCTM, 0.05), 1.3) - TRPM);
dxdt_DCTM = KDCTM * (fmin(1.30, pow(fmax(FDCT, 0.02), 0.75)
                    * (1.0 + EDCTHYP * pow(fmax(0.0, LUMREL - 1.0), 0.4))) - DCTM);

// --- long-term outcome ---------------------------------------------------------
dxdt_NEPH = KNC * fmax(0.0, SATIDX - SSTHR) - KNCOUT * NEPH;
dxdt_GFRF = -KFIB * pow(NEPH / 3.0, 1.2)
            - KHYPOK * fmax(0.0, 1.0 - CK / 3.0)
            - KNSAIDK * INH2 * fmin(1.0, VOLIDX);
dxdt_IGF1 = KIGF * (ghf - IGF1);
dxdt_HTSD = HTVEL / 365.0;
dxdt_GIM  = KGIREP * (1.0 - GIM) - KGIDAM * INH1 * GIM;

$CAPTURE @annotated
CNA    : serum sodium (mmol/L)
CK     : serum potassium (mmol/L)
CCL    : serum chloride (mmol/L)
CHC    : serum bicarbonate (mmol/L)
CMG    : serum magnesium (mmol/L)
CCA    : ultrafilterable calcium (mmol/L)
CKIC   : intracellular potassium (mmol/L)
EGFR   : eGFR (mL/min/1.73m2)
MDSENSE: macula-densa sensed NaCl signal (relative)
LUMREL : luminal NaCl at the macula densa (relative)
DELIVREL: distal Na delivery to the ASDN (relative)
UPGE2  : urinary PGE2 (x upper limit of normal)
UNA    : urinary sodium (mmol/day)
UK     : urinary potassium (mmol/day)
UCL    : urinary chloride (mmol/day)
UCA    : urinary calcium (mmol/day)
UMG    : urinary magnesium (mmol/day)
UCACR  : urinary Ca/creatinine (mmol/mmol)
FEMG   : fractional Mg excretion (%)
FEK    : fractional K excretion (%)
TTKG   : transtubular K gradient
UVOL   : urine volume (L/day)
UOSMA  : achieved urine osmolality (mosm/kg)
QTC    : QTc interval (ms)
CRAMP  : cramp/tetany score (0-10)
FATIG  : fatigue score (0-10)
HTVEL  : height-SDS velocity (SDS/year)
SATIDX : tubular Ca supersaturation index
ULCER  : NSAID gastropathy risk (0-1)
AKIFLAG: acute PG-dependent GFR loss > 20% (0/1)
CIND   : indomethacin concentration (mg/L)
CCEL   : celecoxib concentration (mg/L)
CAMI   : amiloride concentration (mg/L)
CCAN   : canrenone concentration (mg/L)
CACE   : enalaprilat concentration (mg/L)
INH2   : renal COX-2 inhibition (fraction)
INH1   : COX-1 (gastric) inhibition (fraction)
AMIB   : ENaC blockade (fraction)
MRB    : MR blockade (fraction)
ATAL   : effective TAL activity (relative)
ANCC   : effective NCC activity (relative)
PDREL  : TAL lumen-positive PD (relative)
VOLIDX : volume-contraction index (0-1.5)
'

bgs <- mcode_cache("bgs", bgs_code)

## =============================================================================
##  GENOTYPE LIBRARY
##  Only FTAL / FDCT / FROMKCD move. Nothing else is genotype-specific.
## =============================================================================
bgs_genotypes <- function() {
  tibble::tribble(
    ~id, ~genotype,               ~gene,      ~FTAL, ~FDCT, ~FROMKCD, ~PEDS,
    1L,  "Healthy control",       "-",         1.00,  1.00,     1.00,    0,
    2L,  "Bartter I (antenatal)", "SLC12A1",   0.10,  1.00,     1.00,    1,
    3L,  "Bartter II (antenatal)","KCNJ1",     0.25,  1.00,     0.30,    1,
    4L,  "Bartter III (classic)", "CLCNKB",    0.45,  0.55,     1.00,    0,
    5L,  "Bartter IVa (+deaf)",   "BSND",      0.10,  0.35,     1.00,    1,
    6L,  "Bartter IVb (+deaf)",   "CLCNKA/KB", 0.15,  0.45,     1.00,    1,
    7L,  "Bartter V (transient)", "MAGED2",    0.20,  0.70,     1.00,    1,
    8L,  "Gitelman",              "SLC12A3",   1.00,  0.15,     1.00,    0,
    9L,  "Gitelman (severe)",     "SLC12A3",   1.00,  0.05,     1.00,    0,
    10L, "EAST/SeSAME",           "KCNJ10",    0.90,  0.30,     1.00,    1
  )
}

## Convenience: build a parameter list for a named genotype
bgs_pheno <- function(name, ...) {
  g <- bgs_genotypes()
  row <- g[g$genotype == name, ]
  if (nrow(row) == 0) stop("unknown genotype: ", name)
  c(list(FTAL = row$FTAL, FDCT = row$FDCT,
         FROMKCD = row$FROMKCD, PEDS = row$PEDS), list(...))
}

## =============================================================================
##  Burn-in to the genotype steady state (always do this before reading a
##  phenotype or starting therapy).
## =============================================================================
bgs_steady <- function(mod = bgs, pars = list(), days = 600) {
  m <- mod
  if (length(pars)) m <- mrgsolve::param(m, pars)
  out <- mrgsolve::mrgsim(m, end = days, delta = 5, hmax = 0.5)
  tail(as.data.frame(out), 1)
}

## Return a model already sitting at a genotype steady state.
bgs_init_at <- function(mod = bgs, pars = list(), days = 600) {
  ss <- bgs_steady(mod, pars, days)
  cmts <- mrgsolve::cmt(mod)
  inits <- as.list(ss[, cmts])
  m <- mod
  if (length(pars)) m <- mrgsolve::param(m, pars)
  mrgsolve::init(m, inits)
}

## =============================================================================
##  DOSING BUILDERS  (compartment numbers follow the $CMT order)
## =============================================================================
.cmtn <- function(nm) which(mrgsolve::cmt(bgs) == nm)

dose_indomethacin <- function(mg_per_kg_day = 2, bw = 60, n_daily = 3,
                              start = 0, dur_days = 365) {
  tot <- mg_per_kg_day * bw
  mrgsolve::ev(amt = tot / n_daily, cmt = .cmtn("GIND"), time = start,
               ii = 1 / n_daily, addl = n_daily * dur_days - 1)
}
dose_celecoxib <- function(mg_per_kg_day = 4, bw = 60, n_daily = 2,
                           start = 0, dur_days = 365) {
  tot <- mg_per_kg_day * bw
  mrgsolve::ev(amt = tot / n_daily, cmt = .cmtn("GCEL"), time = start,
               ii = 1 / n_daily, addl = n_daily * dur_days - 1)
}
dose_amiloride <- function(mg_day = 10, n_daily = 2, start = 0, dur_days = 365) {
  mrgsolve::ev(amt = mg_day / n_daily, cmt = .cmtn("GAMI"), time = start,
               ii = 1 / n_daily, addl = n_daily * dur_days - 1)
}
dose_spironolactone <- function(mg_day = 100, n_daily = 1, start = 0,
                                dur_days = 365) {
  mrgsolve::ev(amt = mg_day / n_daily, cmt = .cmtn("GSPI"), time = start,
               ii = 1 / n_daily, addl = n_daily * dur_days - 1)
}
dose_enalapril <- function(mg_day = 5, start = 0, dur_days = 365) {
  mrgsolve::ev(amt = mg_day, cmt = .cmtn("GACE"), time = start,
               ii = 1, addl = dur_days - 1)
}
dose_kcl <- function(mmol_day = 60, n_daily = 3, start = 0, dur_days = 365) {
  mrgsolve::ev(amt = mmol_day / n_daily, cmt = .cmtn("GKCL"), time = start,
               ii = 1 / n_daily, addl = n_daily * dur_days - 1)
}
dose_mg <- function(mmol_day = 24, n_daily = 4, start = 0, dur_days = 365) {
  ## Splitting the dose matters: intestinal Mg absorption saturates, and the
  ## unabsorbed remainder causes osmotic diarrhoea that wastes more Mg and K.
  mrgsolve::ev(amt = mmol_day / n_daily, cmt = .cmtn("GMGO"), time = start,
               ii = 1 / n_daily, addl = n_daily * dur_days - 1)
}
dose_nacl <- function(mmol_day = 100, n_daily = 3, start = 0, dur_days = 365) {
  mrgsolve::ev(amt = mmol_day / n_daily, cmt = .cmtn("GNAC"), time = start,
               ii = 1 / n_daily, addl = n_daily * dur_days - 1)
}

## =============================================================================
##  SCENARIOS
## =============================================================================
##  S1  healthy control
##  S2  Gitelman, untreated natural history
##  S3  Bartter III, untreated natural history
##  S4  Bartter I (antenatal), untreated — growth, polyuria, nephrocalcinosis
##  S5  Bartter III + indomethacin 2 mg/kg/day
##  S6  Gitelman + indomethacin 2 mg/kg/day     <- the divergence experiment
##  S7  Bartter I + celecoxib vs indomethacin (renal benefit / GI trade-off)
##  S8  Gitelman + amiloride 10 mg/day
##  S9  Gitelman + spironolactone 100 mg/day
##  S10 Gitelman: KCl alone vs KCl + Mg          <- the Mg-before-K experiment
##  S11 Bartter III combination therapy
##  S12 Bartter I + indomethacin + rhGH (growth endpoint)
##  S13 diagnostic HCTZ challenge (flat in Gitelman)
##  S14 diagnostic furosemide challenge (flat in Bartter I)
##  S15 gastroenteritis on NSAID -> acute-on-chronic AKI
##  S16 Mg dosing: 4 divided doses vs 1 daily dose
##  S17 adherence 100% vs 60%
##  S18 virtual population across the genotype library
## =============================================================================

bgs_run <- function(genotype = "Gitelman", regimen = NULL, days = 730,
                    pars = list(), burnin = 600, delta = 1) {
  p <- c(bgs_pheno(genotype), pars)
  m <- bgs_init_at(bgs, p, days = burnin)
  if (is.null(regimen)) {
    mrgsolve::mrgsim(m, end = days, delta = delta, hmax = 0.25)
  } else {
    mrgsolve::mrgsim(m, events = regimen, end = days, delta = delta, hmax = 0.25)
  }
}

bgs_scenarios <- function(days = 730) {
  bw <- 60
  list(
    S1  = list(label = "Healthy control",
               geno = "Healthy control", rx = NULL),
    S2  = list(label = "Gitelman, untreated",
               geno = "Gitelman", rx = NULL),
    S3  = list(label = "Bartter III, untreated",
               geno = "Bartter III (classic)", rx = NULL),
    S4  = list(label = "Bartter I (antenatal), untreated",
               geno = "Bartter I (antenatal)", rx = NULL,
               pars = list(BW = 12, HTCM = 85, GFR0 = 90)),
    S5  = list(label = "Bartter III + indomethacin 2 mg/kg/d",
               geno = "Bartter III (classic)",
               rx = dose_indomethacin(2, bw, 3, 0, days)),
    S6  = list(label = "Gitelman + indomethacin 2 mg/kg/d (no benefit)",
               geno = "Gitelman",
               rx = dose_indomethacin(2, bw, 3, 0, days)),
    S7  = list(label = "Bartter I + celecoxib 4 mg/kg/d",
               geno = "Bartter I (antenatal)",
               rx = dose_celecoxib(4, 12, 2, 0, days),
               pars = list(BW = 12, HTCM = 85, GFR0 = 90)),
    S8  = list(label = "Gitelman + amiloride 10 mg/d",
               geno = "Gitelman", rx = dose_amiloride(10, 2, 0, days)),
    S9  = list(label = "Gitelman + spironolactone 100 mg/d",
               geno = "Gitelman", rx = dose_spironolactone(100, 1, 0, days)),
    S10 = list(label = "Gitelman + KCl 60 + Mg 24 mmol/d",
               geno = "Gitelman",
               rx = c(dose_kcl(60, 3, 0, days), dose_mg(24, 4, 0, days))),
    S11 = list(label = "Bartter III combination (indo+amil+KCl+Mg+NaCl)",
               geno = "Bartter III (classic)",
               rx = c(dose_indomethacin(2, bw, 3, 0, days),
                      dose_amiloride(10, 2, 0, days),
                      dose_kcl(60, 3, 0, days),
                      dose_mg(20, 4, 0, days),
                      dose_nacl(120, 3, 0, days))),
    S12 = list(label = "Bartter I + indomethacin + rhGH",
               geno = "Bartter I (antenatal)",
               rx = dose_indomethacin(2, 12, 3, 0, days),
               pars = list(BW = 12, HTCM = 85, GFR0 = 90, RHGH = 0.6)),
    S13 = list(label = "Diagnostic HCTZ challenge",
               geno = "Gitelman", rx = NULL,
               pars = list(DXHCTZ = 0.90, DXT0 = 1, DXDUR = 0.25)),
    S14 = list(label = "Diagnostic furosemide challenge",
               geno = "Bartter I (antenatal)", rx = NULL,
               pars = list(DXFURO = 0.90, DXT0 = 1, DXDUR = 0.25)),
    S15 = list(label = "Gastroenteritis on NSAID -> AKI",
               geno = "Bartter III (classic)",
               rx = dose_indomethacin(2, bw, 3, 0, days),
               pars = list(WEXTRA = 2.5, WEXT0 = 30, WEXTDUR = 3)),
    S16 = list(label = "Mg 24 mmol/d as 1 dose (vs 4 in S10)",
               geno = "Gitelman",
               rx = c(dose_kcl(60, 3, 0, days), dose_mg(24, 1, 0, days))),
    S17 = list(label = "Combination at 60% adherence",
               geno = "Bartter III (classic)",
               rx = c(dose_indomethacin(2, bw, 3, 0, days),
                      dose_kcl(60, 3, 0, days), dose_mg(20, 4, 0, days)),
               pars = list(ADHER0 = 0.60))
  )
}

bgs_run_scenario <- function(key, days = 730, delta = 1) {
  sc <- bgs_scenarios(days)[[key]]
  if (is.null(sc)) stop("unknown scenario: ", key)
  out <- bgs_run(sc$geno, sc$rx, days = days,
                 pars = if (is.null(sc$pars)) list() else sc$pars,
                 delta = delta)
  df <- as.data.frame(out)
  df$scenario <- sc$label
  df$key <- key
  df
}

bgs_run_all <- function(keys = names(bgs_scenarios()), days = 730, delta = 7) {
  do.call(rbind, lapply(keys, bgs_run_scenario, days = days, delta = delta))
}

## S18 — virtual population sweep over the genotype library
bgs_vpop <- function(days = 600) {
  g <- bgs_genotypes()
  do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
    ss <- bgs_steady(bgs, list(FTAL = g$FTAL[i], FDCT = g$FDCT[i],
                               FROMKCD = g$FROMKCD[i], PEDS = g$PEDS[i]),
                     days = days)
    data.frame(genotype = g$genotype[i], gene = g$gene[i],
               FTAL = g$FTAL[i], FDCT = g$FDCT[i],
               K = round(ss$CK, 2), Mg = round(ss$CMG, 2),
               HCO3 = round(ss$CHC, 1), Cl = round(ss$CCL, 1),
               uCaCr = round(ss$UCACR, 3), FEMg = round(ss$FEMG, 1),
               uPGE2 = round(ss$UPGE2, 2), PRA = round(ss$PRA, 2),
               MDsense = round(ss$MDSENSE, 3),
               Uvol = round(ss$UVOL, 2), QTc = round(ss$QTC, 0),
               stringsAsFactors = FALSE)
  }))
}

## =============================================================================
##  QUICK-LOOK HELPERS
## =============================================================================
bgs_summary <- function(df) {
  last <- df[df$time == max(df$time), ]
  last[, c("scenario", "CK", "CMG", "CHC", "CCL", "UCACR", "FEMG",
           "UPGE2", "UVOL", "QTC", "CRAMP", "HTSD", "NEPH", "EGFR", "ULCER")]
}

## Example:
##   src <- bgs_run_all(c("S2","S3","S5","S6"), days = 365)
##   bgs_summary(src)
##   print(bgs_vpop())
##
## The key result to look for in bgs_vpop(): urinary Ca/Cr crosses from
## ~0.03 (Gitelman) through ~0.30 (healthy) to ~1.0 (Bartter I), and urinary
## PGE2 rises only in the FTAL-lesioned rows — neither is coded anywhere.
