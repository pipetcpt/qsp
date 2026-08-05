# =============================================================================
#  rvo_mrgsolve_model.R
#  Retinal vein occlusion (CRVO / BRVO) macular oedema — QSP model for mrgsolve
#  망막정맥폐쇄(중심/분지) 황반부종 — mrgsolve 기반 정량적 시스템 약리 모델
# =============================================================================
#
#  THE ONE STRUCTURAL CLAIM
#  ------------------------
#  Retinal thickening in vein occlusion is a Starling flux
#
#      Jv = Lp * S * [ (Pc - Pt) - sigma * (pi_c - pi_t) ]
#           \_____/   \________/   \___________________/
#         PERMEABILITY  PRESSURE          ONCOTIC
#             ARM         ARM               ARM
#
#  and every anti-VEGF term in this model multiplies the FIRST factor.  Not one
#  anti-VEGF term appears inside the bracket.  Pc is set by the occlusion:
#
#      Pc = (Pa * Rv + Pv * Ra) / (Ra + Rv),     Rv = Rv0 (1+OCC)/(1+COLL)
#
#  so the following are consequences of the algebra rather than assertions:
#
#   1. A CRITICAL CAPILLARY PRESSURE EXISTS.  Setting Jv = 0 with the
#      permeability arm at its floor gives
#          Pc* = Pt + sigma_max * (pi_c - pi_t0)  =  31.3 mmHg
#      Above Pc* the bracket is positive with Lp already basal, so no dose of
#      any anti-VEGF agent can dry the macula.  Non-ischaemic CRVO starts above
#      Pc* and crosses below it as collaterals mature; ischaemic CRVO never
#      crosses.  That crossing time, not the choice of agent, separates the two
#      clinical courses.
#
#   2. BECAUSE Pc IS BOUNDED BY THE ARTERIAL INLET PRESSURE, SYSTEMIC BLOOD
#      PRESSURE ENTERS THE PRESSURE ARM NON-LINEARLY.  dPc/dPa = Rv/(Ra+Rv),
#      which tends to 1 as the occlusion tightens: in an occluded eye arterial
#      pressure arrives at the capillary almost undamped.  Lowering Pa from 40
#      to 32 mmHg raises the critical venous resistance Rv* from 1.65 to 21.1.
#
#   3. DURATION OF VEGF SUPPRESSION AND DURATION OF DRYNESS ARE DIFFERENT
#      NUMBERS.  Suppression duration is pharmacology,
#          t_sup = (t_half / ln2) * ln( C0_eff / [(R-1) * alpha * KD] )
#      and for aflibercept it is ~96 days.  Dryness duration is the disease.
#      An eye that is WET while VEGF is still fully suppressed is a pressure-arm
#      eye; switching agent cannot help it, and the model says so explicitly.
#
#   4. OEDEMA IS A STATE, VISION IS AN INTEGRAL.  L_ed (letters lost to
#      thickening) is a reversible function of the current state.  L_pr (letters
#      lost to ellipsoid-zone loss) is the time-integral of thickening times
#      ischaemia and does not come back.  Delay is therefore priced in permanent
#      letters whether or not the macula is eventually dried.
#
#  NOTE ON EXECUTION
#  -----------------
#  The build environment that produced this repository has no R runtime, so
#  this file has not been executed here.  Every equation below is mirrored term
#  for term in `rvo_reference_model.py`, which IS executed; the numbers quoted
#  in README.md come from that reference implementation and the captured log is
#  `rvo_reference_output.txt`.  If you find a discrepancy between the two files,
#  the Python file is the one that was run.
#
#  CALIBRATION TARGETS (see rvo_references.md for the citations)
#  ------------------------------------------------------------
#    CRUISE      ranibizumab 0.5 mg q4w x6, CRVO ..... +14.9 letters, CST -452
#    CRUISE      sham ................................ +0.8 letters
#    CRUISE      sham -> ranibizumab at m6 ........... +7.3 letters at m12
#    COPERNICUS  aflibercept 2 mg q4w x6 ............. +17.3 letters at w24
#    COPERNICUS  sham -> aflibercept ................. +1.5 letters at w100
#    GALILEO     aflibercept ......................... +18.0 w24, +13.7 w76
#    SCORE2      aflibercept +18.9 / bevacizumab +18.6 (no affinity gap)
#    BRAVO       ranibizumab, BRVO ................... +18.3 letters at m6
#    VIBRANT     aflibercept, BRVO ................... +17.0 vs laser +12.2
#    BALATON /
#    COMINO      faricimab ........................... non-inferior to aflibercept
#    GENEVA      dexamethasone implant ............... IOP rise 16%, cataract 30%/y
#    CVOS        ischaemic CRVO -> INV/ANV ........... ~35% within 3 years
#    FRB!        registry injection counts ........... about a quarter of the
#                                                     trial letter gain
# =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

# =============================================================================
#  1.  MODEL SPECIFICATION
# =============================================================================

rvo_code <- '
$PROB
# Retinal vein occlusion (CRVO / BRVO) macular oedema
# 27 ODEs: drug PK (vitreous / systemic), VEGF-PlGF-Ang2-IL6 tone, occlusion and
# collateral haemodynamics, non-perfusion, tight junction, Starling fluid
# balance, irreversible photoreceptor damage, neovascularisation, IOP, lens,
# dexamethasone implant, and visual acuity.

$PARAM @annotated
// ---- ocular geometry -------------------------------------------------------
VVIT   :  4.0   : vitreous volume (mL)
FPEN   :  0.080 : retina-to-vitreous drug partition coefficient (-)
VP     :  4.0   : systemic volume of distribution (L)

// ---- drug identity (set per simulation; defaults = aflibercept 2 mg) --------
TVIT   :  9.10  : apparent intraocular half-life of the drug (d)
TSYS   :  0.50  : systemic half-life of free drug (d)
KDV    :  0.49  : intrinsic KD for VEGF-A165 (pM)
KDP    : 39.0   : intrinsic KD for PlGF (pM; 1e12 = does not bind)
KDA    :  1e12  : intrinsic KD for Ang-2 (pM; 1e12 = does not bind)
MWDRUG : 115000 : drug molecular weight (Da), for reporting only

// ---- the single calibrated PD scale factor ---------------------------------
// ALPHA converts an INTRINSIC KD measured by SPR in buffer into an apparent
// in-vivo IC50 at the inner blood-retina barrier.  It exceeds 1 because the
// drug must out-compete VEGFR2 (itself sub-nM), because the VEGF concentration
// that matters is a local interstitial flux rather than a vitreous average,
// and because of diffusional limitation across the retina.  ONE value is used
// for every agent, so the RELATIVE ranking of agents is fixed entirely by
// their measured KD, molar dose and half-life.
ALPHA  : 120.0  : apparent-IC50 multiplier on intrinsic KD (-)
RADEQ  :   4.9  : VEGF suppression ratio counted as adequate (-)
PERF   :   0    : 1 = infinite blockade experiment (-)
PFROM  :  42.0  : time at which infinite blockade begins (d)

// ---- retinal haemodynamics -------------------------------------------------
PA     : 40.0   : post-arteriolar inlet pressure (mmHg)
RA0    :  1.0   : arteriolar resistance, relative (-)
RV0    :  0.25  : venular plus venous resistance, relative (-)
PT     : 12.0   : retinal interstitial hydrostatic pressure (mmHg)
PIC    : 25.0   : plasma colloid osmotic pressure (mmHg)
PIT0   :  4.0   : interstitial colloid osmotic pressure when dry (mmHg)
SIGMAX :  0.92  : reflection coefficient of an intact inner BRB (-)
AUTOREG:  0.35  : fraction of arteriolar resistance lost to hypoxic dilatation (-)
IOP0   : 15.0   : baseline intraocular pressure (mmHg)

// ---- occlusion natural history ---------------------------------------------
OCC0   : 10.0   : occlusion severity multiplier at onset (-)
OCCRES :  0.55  : non-lysable residual fraction of the occlusion (-)
KLYS   :  0.0045: thrombus lysis rate constant (1/d)
KCOLL  :  0.0035: collateral growth rate constant (1/d)
COLLMAX:  0.70  : ceiling on collateral conductance (-)
HIFCOLL:  2.0   : hypoxic amplification of collateral growth (-)

// ---- hypoxia and ligand tone -----------------------------------------------
TAUHIF :  2.0   : HIF-1alpha time constant (d)
RVVEGF : 24.0   : basal VEGF-A production (pM/d)
KDEGV  :  8.0   : VEGF-A turnover (1/d)
EHIFV  : 12.0   : fold VEGF-A induction at HIF = 1 (-)
EABLAV :  0.72  : fraction of VEGF production removable by PRP (-)
EDEXV  :  0.70  : fraction of VEGF production suppressed by dexamethasone (-)
RPPLGF :  6.0   : basal PlGF production (pM/d)
KDEGP  :  3.0   : PlGF turnover (1/d)
EHIFP  :  5.0   : fold PlGF induction at HIF = 1 (-)
RAANG2 : 10.0   : basal Ang-2 production (pM/d)
KDEGA  :  2.0   : Ang-2 turnover (1/d)
EHIFA  :  7.0   : fold Ang-2 induction at HIF = 1 (-)
EDEXA  :  0.50  : fraction of Ang-2 suppressed by dexamethasone (-)
RIIL6  : 30.0   : basal IL-6 production (pg/mL/d)
KDEGI  :  6.0   : IL-6 turnover (1/d)
EHIFI  : 16.0   : fold IL-6 induction at HIF = 1 (-)
EDEXI  :  0.80  : fraction of IL-6 suppressed by dexamethasone (-)

// ---- non-perfusion and leukostasis -----------------------------------------
KNP    :  0.30  : non-perfusion growth gain (DA/d)
NPMAX  : 75.0   : retinal area available to close (disc areas)
QNP    :  0.45  : relative perfusion below which capillaries close (-)
NPLEUK :  0.80  : leukostasis amplification of capillary closure (-)
KNPREP :  0.0022: collateral-driven reperfusion rate (1/d)
KMI    :  0.020 : macular ischaemia growth gain (1/d)
QMI    :  0.35  : relative perfusion below which the fovea decompensates (-)
KMIREP :  0.0020: macular ischaemia reversal (1/d)
TAULEUK:  5.0   : leukostasis time constant (d)
WLEUKI6:  1.0   : IL-6 weight on leukostasis (-)
WLEUKV :  0.50  : VEGF weight on leukostasis (-)

// ---- tight junction and permeability ---------------------------------------
KREPTJ :  0.160 : junction re-assembly rate constant (1/d)
KDISTJ :  0.750 : junction disassembly rate constant at DRIVE = 1 (1/d)
TIE2B  :  0.50  : extra repair when Tie2 is unopposed (-)
EDEXTJ :  0.35  : direct steroid up-regulation of occludin / claudin-5 (-)
KTIE2  : 15.0   : Ang-2 for half inhibition of Tie2 (pM)
WV     :  0.80  : VEGF-A weight in the permeability drive (-)
KV     : 12.0   : VEGF-A potency in the permeability drive (pM)
WI6    :  0.08  : IL-6 weight (-)
KI6    : 45.0   : IL-6 potency (pg/mL)
WA2    :  0.06  : Ang-2 weight (-)
KA2    : 15.0   : Ang-2 potency (pM)
WPL    :  0.04  : PlGF weight (-)
KPL    :  8.0   : PlGF potency (pM)
WPRESS :  0.06  : mechanotransduction weight (-)
PPRESS0: 28.0   : capillary pressure at which mechanotransduction starts (mmHg)
PPRESSW: 12.0   : width of the mechanotransduction ramp (mmHg)
LPMAX  :  6.5   : maximum relative hydraulic conductivity (-)
LPN    :  1.30  : exponent linking junction loss to conductivity (-)
SIGN   :  0.55  : exponent linking junction integrity to sigma (-)
PITGAIN:  0.75  : protein-leak gain on interstitial oncotic pressure (-)

// ---- chronic remodelling ---------------------------------------------------
KCHRON :  0.0016: chronic remodelling rate (1/d)
WCHRON0: 70.0   : oedema above which remodelling starts (um)
WCHRONW: 300.0  : width of the remodelling ramp (um)
KCHRONR:  0.0006: chronic remodelling reversal (1/d)
TJMAXL :  0.45  : ceiling lost to chronic remodelling (-)

// ---- fluid balance ---------------------------------------------------------
KF     :  0.870 : filtration coefficient (um per d per mmHg per unit Lp)
KOUTLIN:  0.100 : linear water clearance (1/d)
VMAXPMP: 32.0   : saturable Mueller-AQP4 / RPE transport (um/d)
KMPUMP : 45.0   : Michaelis constant of the water pump (um)
FSRF   :  0.18  : fraction of filtered flux entering the subretinal space (-)
KSRFLIN:  0.060 : linear subretinal fluid clearance (1/d)
VMAXSRF: 14.0   : saturable subretinal fluid clearance (um/d)
KMSRF  : 30.0   : Michaelis constant of subretinal clearance (um)

// ---- irreversible damage ---------------------------------------------------
KEZ    :  0.0018: ellipsoid-zone damage rate (1/d)
KMEZ   : 250.0  : oedema for half maximal damage (um)
WEZMI  :  1.10  : macular ischaemia amplification of damage (-)
KEZREP :  0.0004: ellipsoid-zone repair rate (1/d)
KMEZR  : 150.0  : oedema above which repair is suppressed (um)
KDRIL  :  0.0035: DRIL formation rate (1/d)
KDRILR :  0.00025: DRIL reversal (1/d)
CSTCUM : 320.0  : thickness above which cumulative exposure accrues (um)

// ---- neovascularisation, IOP, lens -----------------------------------------
KNVI   :  0.030 : anterior neovascularisation growth (1/d)
KMNVINP: 10.0   : non-perfusion for half maximal NV drive (DA)
KMNVIV : 16.0   : VEGF for half maximal NV drive (pM)
KNVIREG:  0.012 : NV regression rate when VEGF is suppressed (1/d)
IOPNVI : 24.0   : IOP rise at full angle involvement (mmHg)
NVITHR :  0.40  : NVI above which the angle starts to fail (-)
TAUIOP :  4.0   : IOP time constant (d)
IOPDEX : 10.0   : steroid IOP rise in a responder (mmHg)
KMIOPD : 200.0  : dexamethasone for half maximal IOP rise (ng/mL)
STRESP :  1.0   : steroid responder indicator (0 or 1)
KCATDEX:  0.0024: steroid cataractogenesis rate (1/d)
KMCATD : 150.0  : dexamethasone for half maximal cataractogenesis (ng/mL)
KCATAGE:  0.00008: age-related lens opacification (1/d)

// ---- dexamethasone implant -------------------------------------------------
KRELDEX:  0.0128: implant release rate constant (1/d)
KELDEX :  2.50  : vitreous dexamethasone clearance (1/d)
KMDEXE : 120.0  : dexamethasone for half maximal steroid effect (ng/mL)

// ---- visual function -------------------------------------------------------
BCVACEI: 80.0   : letters achievable with a dry intact macula (letters)
EED    : 36.0   : maximum reversible letters lost to thickening (letters)
KMED   : 230.0  : thickening for half maximal reversible loss (um)
EPR    : 42.0   : letters lost per unit ellipsoid-zone loss (letters)
EISCH  : 22.0   : maximum letters lost to macular ischaemia (letters)
KMISCH :  1.20  : macular ischaemia index for half maximal loss (-)
ECAT   : 12.0   : letters lost at full lens opacity (letters)
EPRP   :  6.0   : letters lost at full panretinal photocoagulation (letters)
TAUBCVA: 10.0   : functional lag on visual acuity (d)

// ---- sector involvement (1.0 = CRVO, < 1 = BRVO) ---------------------------
// Two different fractions are required.  FSECFLX is the share of the central
// subfield inside the occluded sector (drives fluid and the VEGF source).
// FSECDMG is the share of the FOVEA functionally at risk, which is smaller: a
// BRVO can flood the central subfield while sparing part of the foveal
// photoreceptor mosaic.  Collapsing them into one number is what makes naive
// BRVO models predict impossibly good vision.
FSECFLX:  1.0   : sector fraction driving flux and VEGF production (-)
FSECDMG:  1.0   : sector fraction of the fovea at functional risk (-)

// ---- laser -----------------------------------------------------------------
BASECST: 250.0  : dry central subfield thickness (um)

$CMT @annotated
DVIT  : anti-VEGF drug in the vitreous (nmol)
ASYS  : anti-VEGF drug in the systemic circulation (nmol)
VTONE : VEGF-A tone, free VEGF absent drug (pM)
PTONE : PlGF tone (pM)
ATONE : Angiopoietin-2 tone (pM)
IL6   : interleukin-6 (pg/mL)
HIF   : HIF-1alpha signalling (0-1)
OCC   : occlusion severity multiplier (-)
COLL  : collateral conductance (-)
NP    : retinal capillary non-perfusion (disc areas)
MI    : macular ischaemia index (0-3)
LEUK  : leukostasis (0-1)
TJ    : tight junction integrity (0-1)
CHRON : chronic remodelling (0-1)
W     : intraretinal excess water (um)
SRF   : subretinal fluid (um)
EZ    : ellipsoid-zone integrity (0-1)
DRIL  : disorganisation of retinal inner layers (0-1)
CUMED : cumulative oedema exposure (um.d)
NVI   : anterior segment neovascularisation (0-1)
IOP   : intraocular pressure (mmHg)
CAT   : lens opacity (0-1)
IMP   : dexamethasone remaining in the implant (ng)
CDEX  : vitreous dexamethasone concentration (ng/mL)
BCVAO : observed best-corrected visual acuity (letters)
ABLA  : fraction of retina ablated by photocoagulation (0-1)
TSUP  : cumulative days of adequate VEGF suppression (d)

$GLOBAL
#define HILL(x, k) ( ((x) <= 0.0) ? 0.0 : (x) / ((k) + (x)) )
#define CLAMP(x, a, b) ( (x) < (a) ? (a) : ( (x) > (b) ? (b) : (x) ) )
#define POSPOW(x, p) ( ((x) <= 0.0) ? 0.0 : pow((x), (p)) )

// File-scope, because a local declared in $MAIN is NOT visible inside $ODE.
double KELVIT;
double KELSYS;

$MAIN
// ---- initial conditions ----------------------------------------------------
VTONE_0 = RVVEGF / KDEGV;      // 3 pM  ~ 135 pg/mL in a normal eye
PTONE_0 = RPPLGF / KDEGP;      // 2 pM
ATONE_0 = RAANG2 / KDEGA;      // 5 pM
IL6_0   = RIIL6  / KDEGI;      // 5 pg/mL
HIF_0   = 0.0;
OCC_0   = OCC0;                // the occlusion is present at t = 0
COLL_0  = 0.0;
TJ_0    = 1.0;
EZ_0    = 1.0;
IOP_0   = IOP0;
CAT_0   = 0.0;
BCVAO_0 = BCVACEI;

KELVIT = log(2.0) / TVIT;
KELSYS = log(2.0) / TSYS;

$ODE
// ===========================================================================
//  A.  DRUG CONCENTRATIONS AND COMPETITIVE INHIBITION
// ===========================================================================
double CVIT = DVIT / VVIT * 1.0e6;          // nmol/mL -> pM
double CRET = CVIT * FPEN;                  // effect-site concentration (pM)
double CPNM = ASYS / VP * 1.0e3;            // systemic concentration (nM)

int    PB   = (PERF > 0.5 && SOLVERTIME >= PFROM) ? 1 : 0;
double RINV = PB ? 1.0e12 : CRET / (KDV * ALPHA);
double RINP = PB ? 1.0e12 : CRET / (KDP * ALPHA);
double RINA = PB ? 1.0e12 : CRET / (KDA * ALPHA);

double SUPPR = 1.0 + RINV;                  // VEGF suppression ratio
double VACT  = VTONE / SUPPR;               // VEGF signalling actually delivered
double PACT  = PTONE / (1.0 + RINP);
double AACT  = ATONE / (1.0 + RINA);

double ESTER = HILL(CDEX, KMDEXE);          // steroid effect (0-1)

// ===========================================================================
//  B.  HAEMODYNAMICS — THE PRESSURE ARM
// ===========================================================================
double RARES = RA0 * (1.0 - AUTOREG * HIF);            // hypoxic dilatation
double RVRES = RV0 * (1.0 + OCC) / (1.0 + COLL);       // occlusion vs collateral
double PV    = IOP + 2.0;
double PC    = (PA * RVRES + PV * RARES) / (RARES + RVRES);
double QQ    = (PA - PV) / (RARES + RVRES);
double Q0    = (PA - (IOP0 + 2.0)) / (RA0 + RV0);
double QREL  = QQ / Q0;
double HYPOX = CLAMP(1.0 - POSPOW(QREL, 0.80), 0.0, 1.0);

// ===========================================================================
//  C.  BARRIER — THE PERMEABILITY AND ONCOTIC ARMS
// ===========================================================================
double TJs   = CLAMP(TJ, 1.0e-3, 1.0);
double LPREL = 1.0 + (LPMAX - 1.0) * POSPOW(1.0 - TJs, LPN);
double SIGMA = SIGMAX * POSPOW(TJs, SIGN);
double LEAK  = 1.0 - SIGMA / SIGMAX;
double PITT  = PIT0 + (PIC - PIT0) * LEAK * PITGAIN;

double DRIVE = WV  * HILL(VACT, KV)
             + WI6 * HILL(IL6,  KI6)
             + WA2 * HILL(AACT, KA2)
             + WPL * HILL(PACT, KPL)
             + WPRESS * CLAMP((PC - PPRESS0) / PPRESSW, 0.0, 1.0);

// ===========================================================================
//  D.  THE STARLING NODE
// ===========================================================================
double BRACKET = (PC - PT) - SIGMA * (PIC - PITT);
double JV      = KF * LPREL * BRACKET * FSECFLX;

// the same flux with the permeability arm forced to its floor: the component
// that no anti-VEGF agent can remove
double JVFLOOR = KF * 1.0 * ((PC - PT) - SIGMAX * (PIC - PIT0)) * FSECFLX;
double PCCRIT  = PT + SIGMA * (PIC - PITT);

double Wp   = (W > 0.0) ? W : 0.0;
double JOUT = KOUTLIN * Wp + VMAXPMP * HILL(Wp, KMPUMP);
double SRFp = (SRF > 0.0) ? SRF : 0.0;
double WTOT = Wp + 0.6 * SRFp;
double CST  = BASECST + WTOT;

// ===========================================================================
//  E.  DIFFERENTIAL EQUATIONS
// ===========================================================================
// --- drug disposition ------------------------------------------------------
dxdt_DVIT = -KELVIT * DVIT;
dxdt_ASYS =  KELVIT * DVIT - KELSYS * ASYS;

// --- hypoxia signalling ----------------------------------------------------
double HIFT = CLAMP(HYPOX * (1.0 + 0.60 * NP / 30.0), 0.0, 1.0);
dxdt_HIF = (HIFT - HIF) / TAUHIF;

// --- ligand tone -----------------------------------------------------------
double DEXV = 1.0 - EDEXV * ESTER;
dxdt_VTONE = RVVEGF * (1.0 + EHIFV * HIF * FSECFLX)
             * (1.0 - EABLAV * ABLA) * DEXV - KDEGV * VTONE;
dxdt_PTONE = RPPLGF * (1.0 + EHIFP * HIF * FSECFLX) * DEXV - KDEGP * PTONE;
dxdt_ATONE = RAANG2 * (1.0 + EHIFA * HIF * FSECFLX)
             * (1.0 - EDEXA * ESTER) - KDEGA * ATONE;
dxdt_IL6   = RIIL6  * (1.0 + EHIFI * HIF * FSECFLX)
             * (1.0 - EDEXI * ESTER) - KDEGI * IL6;

// --- occlusion and collaterals ---------------------------------------------
// The residual (non-lysable) stenosis is the single most consequential
// natural-history parameter: an occlusion that lyses completely cures itself
// and makes every therapeutic comparison meaningless.
double OCCFLOOR = OCCRES * OCC0;
dxdt_OCC  = -KLYS * ((OCC > OCCFLOOR) ? (OCC - OCCFLOOR) : 0.0);
dxdt_COLL =  KCOLL * (1.0 + HIFCOLL * HIF)
             * ((COLLMAX > COLL) ? (COLLMAX - COLL) : 0.0);

// --- non-perfusion, leukostasis, macular ischaemia -------------------------
double LEUKT = CLAMP(WLEUKI6 * HILL(IL6, KI6) + WLEUKV * HILL(VACT, KV), 0.0, 1.0);
dxdt_LEUK = (LEUKT - LEUK) / TAULEUK;

double STASIS = (QNP > QREL) ? (QNP - QREL) : 0.0;
dxdt_NP = KNP * POSPOW(STASIS, 1.20) * (1.0 + NPLEUK * LEUK)
          * ((NP < NPMAX) ? (1.0 - NP / NPMAX) : 0.0)
          - KNPREP * NP * HILL(COLL, 1.0);

double STASISM = (QMI > QREL) ? (QMI - QREL) : 0.0;
dxdt_MI = KMI * STASISM * (1.0 + 0.5 * NP / 20.0) - KMIREP * MI;

// --- tight junctions and chronic remodelling -------------------------------
double TIE2FREE = 1.0 / (1.0 + AACT / KTIE2);
// The steroid term is a junction effect anti-VEGF has no mechanism for, and it
// is why a steroid can work in an eye that has failed anti-VEGF.
double KREP     = KREPTJ * (1.0 + TIE2B * TIE2FREE + EDEXTJ * ESTER);
double TJMAX    = 1.0 - TJMAXL * CHRON;
dxdt_TJ = KREP * (TJMAX - TJ) - KDISTJ * DRIVE * TJ;

double CHRD = CLAMP((WTOT - WCHRON0) / WCHRONW, 0.0, 1.5);
dxdt_CHRON = KCHRON * CHRD * (1.0 - CHRON) - KCHRONR * CHRON;

// --- fluid -----------------------------------------------------------------
dxdt_W = JV - JOUT;
dxdt_SRF = FSRF * ((JV > 0.0) ? JV : 0.0) * (1.0 - TJs)
           - (KSRFLIN * SRFp + VMAXSRF * HILL(SRFp, KMSRF));

// --- irreversible damage ---------------------------------------------------
double DMG = KEZ * HILL(WTOT, KMEZ) * (1.0 + WEZMI * MI) * EZ * FSECDMG;
double REP = KEZREP * (1.0 - EZ) * (1.0 - HILL(WTOT, KMEZR));
dxdt_EZ   = REP - DMG;
dxdt_DRIL = KDRIL * HILL(WTOT, KMEZ) * (1.0 - DRIL) - KDRILR * DRIL;
dxdt_CUMED = (CST > CSTCUM) ? (CST - CSTCUM) : 0.0;

// --- neovascularisation, IOP, lens -----------------------------------------
dxdt_NVI = KNVI * HILL(NP, KMNVINP) * HILL(VACT, KMNVIV) * (1.0 - NVI)
           - KNVIREG * NVI * (1.0 - HILL(VACT, KMNVIV));
double IOPT = IOP0
            + IOPNVI * CLAMP((NVI - NVITHR) / (1.0 - NVITHR), 0.0, 1.0)
            + IOPDEX * HILL(CDEX, KMIOPD) * STRESP;
dxdt_IOP = (IOPT - IOP) / TAUIOP;
dxdt_CAT = (KCATDEX * HILL(CDEX, KMCATD) + KCATAGE) * (1.0 - CAT);

// --- dexamethasone implant -------------------------------------------------
dxdt_IMP  = -KRELDEX * IMP;
dxdt_CDEX =  KRELDEX * IMP / VVIT - KELDEX * CDEX;

// --- visual acuity, ablation, trackers -------------------------------------
double LED   = EED * HILL(WTOT, KMED);
double LPR   = EPR * (1.0 - EZ) * FSECDMG;
double LISCH = EISCH * HILL(MI, KMISCH) * FSECDMG;
double LCAT  = ECAT * CAT;
double LPRP  = EPRP * ABLA;
double BALG  = BCVACEI - LED - LPR - LISCH - LCAT - LPRP;
if (BALG < 0.0) BALG = 0.0;
dxdt_BCVAO = (BALG - BCVAO) / TAUBCVA;

dxdt_ABLA = 0.0;                       // moved only by discrete PRP events
dxdt_TSUP = (SUPPR >= RADEQ) ? 1.0 : 0.0;

$TABLE
double CVITt = DVIT / VVIT * 1.0e6;
double CRETt = CVITt * FPEN;
double CPNMt = ASYS / VP * 1.0e3;
int    PBt   = (PERF > 0.5 && TIME >= PFROM) ? 1 : 0;
double SUPPRt = PBt ? 1.0e12 : 1.0 + CRETt / (KDV * ALPHA);
double VACTt  = VTONE / SUPPRt;
double AACTt  = PBt ? 0.0 : ATONE / (1.0 + CRETt / (KDA * ALPHA));

double RAREt = RA0 * (1.0 - AUTOREG * HIF);
double RVREt = RV0 * (1.0 + OCC) / (1.0 + COLL);
double PVt   = IOP + 2.0;
double PCt   = (PA * RVREt + PVt * RAREt) / (RAREt + RVREt);
double QQt   = (PA - PVt) / (RAREt + RVREt);
double Q0t   = (PA - (IOP0 + 2.0)) / (RA0 + RV0);
double QRELt = QQt / Q0t;
double DPCDPA = RVREt / (RAREt + RVREt);

double TJt    = CLAMP(TJ, 1.0e-3, 1.0);
double LPRELt = 1.0 + (LPMAX - 1.0) * POSPOW(1.0 - TJt, LPN);
double SIGMAt = SIGMAX * POSPOW(TJt, SIGN);
double PITTt  = PIT0 + (PIC - PIT0) * (1.0 - SIGMAt / SIGMAX) * PITGAIN;
double BRACKt = (PCt - PT) - SIGMAt * (PIC - PITTt);
double JVt    = KF * LPRELt * BRACKt * FSECFLX;
double JVFLRt = KF * ((PCt - PT) - SIGMAX * (PIC - PIT0)) * FSECFLX;
double PCSTAR = PT + SIGMAX * (PIC - PIT0);       // 31.32 mmHg
double PCCRIt = PT + SIGMAt * (PIC - PITTt);

double WTOTt = ((W > 0.0) ? W : 0.0) + 0.6 * ((SRF > 0.0) ? SRF : 0.0);
double CSTt  = BASECST + WTOTt;
double LEDt  = EED * HILL(WTOTt, KMED);
double LPRt  = EPR * (1.0 - EZ) * FSECDMG;
double LISCt = EISCH * HILL(MI, KMISCH) * FSECDMG;
double BCVA  = BCVAO;
double SNELL = 20.0 * pow(10.0, (85.0 - BCVA) / 50.0);
double NVG   = (NVI > 0.50 && IOP > 25.0) ? 1.0 : 0.0;
double DRY   = (CSTt < 300.0) ? 1.0 : 0.0;
double FLOOR = (JVFLRt > 0.0) ? 1.0 : 0.0;
double CVITUG = DVIT / VVIT * MWDRUG / 1.0e6;     // ug/mL
double CPUG   = CPNMt * MWDRUG / 1.0e6;           // ug/mL

$CAPTURE @annotated
CSTt   : central subfield thickness (um)
BCVA   : best-corrected visual acuity (ETDRS letters)
SNELL  : Snellen denominator (20/x)
PCt    : retinal capillary hydrostatic pressure (mmHg)
PCSTAR : critical capillary pressure with an intact barrier (mmHg)
PCCRIt : critical capillary pressure at the current barrier state (mmHg)
FLOOR  : 1 if a pressure-arm oedema floor exists (-)
DRY    : 1 if the central subfield is below 300 um (-)
JVt    : Starling filtration flux (um/d)
JVFLRt : irreducible pressure-arm flux (um/d)
BRACKt : the Starling bracket (mmHg)
LPRELt : relative hydraulic conductivity (-)
SIGMAt : reflection coefficient (-)
QRELt  : relative retinal blood flow (-)
RVREt  : venous resistance (-)
DPCDPA : transmission of arterial pressure to the capillary (-)
VACTt  : VEGF signalling delivered (pM)
AACTt  : Ang-2 signalling delivered (pM)
SUPPRt : VEGF suppression ratio (-)
CRETt  : effect-site drug concentration (pM)
CVITUG : vitreous drug concentration (ug/mL)
CPUG   : plasma drug concentration (ug/mL)
LEDt   : reversible letters lost to thickening (letters)
LPRt   : irreversible letters lost to photoreceptor loss (letters)
LISCt  : letters lost to macular ischaemia (letters)
NVG    : 1 if neovascular glaucoma criteria are met (-)
'

mod <- mcode("rvo", rvo_code, end = 1095, delta = 1, atol = 1e-8, rtol = 1e-6)

# =============================================================================
#  2.  DRUG LIBRARY
# =============================================================================
#  MW    Da
#  dose  mg per intravitreal injection (label dose)
#  tvit  apparent intraocular half-life, days (aqueous-sampling estimates)
#  kdv   intrinsic KD for VEGF-A165, pM (Papadopoulos 2012 for AFL/RBZ/BEV)
#  kdp   KD for PlGF, pM (1e12 = does not bind)
#  kda   KD for Ang-2, pM (1e12 = does not bind)
#  tsys  systemic half-life of free drug, days
# -----------------------------------------------------------------------------
DRUGS <- tibble::tribble(
  ~drug,          ~label,                  ~mw,     ~dose, ~tvit, ~kdv,  ~kdp,  ~kda,  ~tsys,
  "ranibizumab",  "Ranibizumab 0.5 mg",    48000,   0.50,  7.19,  46.0,  1e12,  1e12,  0.083,
  "bevacizumab",  "Bevacizumab 1.25 mg",   149000,  1.25,  9.82,  58.0,  1e12,  1e12,  20.0,
  "aflibercept",  "Aflibercept 2 mg",      115000,  2.00,  9.10,  0.49,  39.0,  1e12,  0.50,
  "aflibercept8", "Aflibercept 8 mg",      115000,  8.00,  9.10,  0.49,  39.0,  1e12,  0.50,
  "brolucizumab", "Brolucizumab 6 mg",     26000,   6.00,  5.10,  28.4,  1e12,  1e12,  4.40,
  "faricimab",    "Faricimab 6 mg",        149000,  6.00,  7.50,  30.0,  1e12,  1500,  7.50
)

molar_dose <- function(drug) {
  d <- DRUGS[DRUGS$drug == drug, ]
  d$dose / d$mw * 1e6            # nmol
}

drug_params <- function(drug) {
  d <- DRUGS[DRUGS$drug == drug, ]
  list(TVIT = d$tvit, TSYS = d$tsys, KDV = d$kdv, KDP = d$kdp,
       KDA = d$kda, MWDRUG = d$mw)
}

# =============================================================================
#  3.  PHENOTYPES
# =============================================================================
#  OCC0     occlusion severity multiplier on venous resistance at t = 0
#  OCCRES   fraction of OCC0 that never lyses (the permanent stenosis)
#  COLLMAX  ceiling on collateral conductance — the BRVO / CRVO asymmetry
# -----------------------------------------------------------------------------
PHENOTYPES <- list(
  crvo_nonisch = list(OCC0 = 10.0, OCCRES = 0.55, KLYS = 0.0045,
                      COLLMAX = 0.70, FSECFLX = 1.00, FSECDMG = 1.00),
  crvo_isch    = list(OCC0 = 32.0, OCCRES = 0.80, KLYS = 0.0028,
                      COLLMAX = 0.45, FSECFLX = 1.00, FSECDMG = 1.00),
  brvo         = list(OCC0 = 22.0, OCCRES = 0.32, KLYS = 0.0060,
                      COLLMAX = 2.20, FSECFLX = 0.85, FSECDMG = 0.60),
  crvo_mild    = list(OCC0 = 5.0,  OCCRES = 0.20, KLYS = 0.0080,
                      COLLMAX = 1.50, FSECFLX = 1.00, FSECDMG = 1.00)
)

BASE_START <- 42      # day of presentation / randomisation (6 weeks post-onset)
HORIZON    <- 1095    # 36 months

# =============================================================================
#  4.  EVENT BUILDERS
# =============================================================================
ivt_events <- function(drug, times) {
  if (length(times) == 0) return(NULL)
  ev(time = times, amt = molar_dose(drug), cmt = "DVIT")
}

dex_events <- function(times) {
  if (length(times) == 0) return(NULL)
  ev(time = times, amt = 700000, cmt = "IMP")     # 0.7 mg Ozurdex, in ng
}

prp_event <- function(time, frac = 0.35) {
  if (is.null(time)) return(NULL)
  ev(time = time, amt = frac, cmt = "ABLA")
}

q4w <- function(n, start = BASE_START, q = 28) start + (seq_len(n) - 1) * q

load_then <- function(n_load, q_load, n_maint, q_maint, start = BASE_START) {
  t1 <- q4w(n_load, start, q_load)
  c(t1, tail(t1, 1) + seq_len(n_maint) * q_maint)
}

run_case <- function(phenotype = "crvo_nonisch", drug = "aflibercept",
                     ivt_times = numeric(0), dex_times = numeric(0),
                     prp_time = NULL, perfect = FALSE, pfrom = BASE_START,
                     extra = list(), end = HORIZON) {
  p <- c(PHENOTYPES[[phenotype]], drug_params(drug), extra,
         list(PERF = as.numeric(perfect), PFROM = pfrom))
  e <- Reduce(function(a, b) if (is.null(a)) b else if (is.null(b)) a else a + b,
              list(ivt_events(drug, ivt_times), dex_events(dex_times),
                   prp_event(prp_time)))
  m <- mod %>% param(p)
  out <- if (is.null(e)) m %>% mrgsim(end = end, delta = 1)
         else m %>% ev(e) %>% mrgsim(end = end, delta = 1)
  as_tibble(out)
}

# =============================================================================
#  5.  PRN AND TREAT-AND-EXTEND CONTROLLERS
# =============================================================================
#  A fixed event list cannot express "treat if wet at the monthly visit", so the
#  regimen is simulated visit by visit: integrate to the next visit, read CST,
#  decide, carry the state forward.  This is the same logic as the `Regimen`
#  class in rvo_reference_model.py.
# -----------------------------------------------------------------------------
run_reactive <- function(phenotype = "crvo_nonisch", drug = "aflibercept",
                         mode = c("tae", "prn"), load_n = 6, load_q = 28,
                         start = BASE_START, interval0 = 56, imin = 28,
                         imax = 112, istep = 28, cst_thr = 310,
                         dex_times = numeric(0), max_inj = Inf,
                         stop_at = Inf, extra = list(), end = HORIZON) {
  mode <- match.arg(mode)
  p <- c(PHENOTYPES[[phenotype]], drug_params(drug), extra,
         list(PERF = 0, PFROM = 1e9))
  m <- mod %>% param(p)
  dose <- molar_dose(drug)

  state <- NULL
  t_now <- 0
  interval <- interval0
  n_inj <- 0
  inj_times <- numeric(0)
  chunks <- list()
  next_due <- start

  repeat {
    t_next <- min(next_due, end)
    seg <- if (is.null(state)) {
      m %>% mrgsim(end = t_next, delta = 1)
    } else {
      m %>% init(state) %>% mrgsim(start = t_now, end = t_next, delta = 1)
    }
    seg <- as_tibble(seg)
    chunks[[length(chunks) + 1]] <- seg
    last <- tail(seg, 1)
    state <- as.list(last[, names(init(mod))])
    t_now <- t_next
    if (t_now >= end) break

    treat <- FALSE
    if (n_inj < max_inj && t_now <= stop_at) {
      if (n_inj < load_n) {
        treat <- TRUE
        next_due <- t_now + load_q
      } else {
        wet <- last$CSTt > cst_thr
        if (mode == "prn") {
          treat <- wet
          next_due <- t_now + load_q             # monthly monitoring
        } else {
          interval <- if (wet) max(imin, interval - istep)
                      else      min(imax, interval + istep)
          treat <- TRUE
          next_due <- t_now + interval
        }
      }
    } else {
      next_due <- t_now + load_q
    }
    if (treat) {
      state$DVIT <- state$DVIT + dose
      n_inj <- n_inj + 1
      inj_times <- c(inj_times, t_now)
    }
    for (dt in dex_times) if (abs(dt - t_now) < 0.5) state$IMP <- state$IMP + 700000
  }
  out <- bind_rows(chunks) %>% distinct(time, .keep_all = TRUE) %>% arrange(time)
  attr(out, "n_inj") <- n_inj
  attr(out, "inj_times") <- inj_times
  out
}

# =============================================================================
#  6.  SCENARIOS (31 arms)
# =============================================================================
SCENARIOS <- list(

  # --- natural history ------------------------------------------------------
  NAT_NONISCH = function() run_case("crvo_nonisch"),
  NAT_ISCH    = function() run_case("crvo_isch"),
  NAT_BRVO    = function() run_case("brvo"),

  # --- trial-like anti-VEGF regimens ---------------------------------------
  # CRUISE / HORIZON: ranibizumab 0.5 mg monthly x6 then monthly PRN
  RBZ_M6_PRN  = function() run_reactive("crvo_nonisch", "ranibizumab",
                                       mode = "prn", load_n = 6),
  # COPERNICUS / GALILEO: aflibercept monthly x6 then extend
  AFL_M6_Q8   = function() run_reactive("crvo_nonisch", "aflibercept",
                                       mode = "tae", load_n = 6),
  # SCORE2: bevacizumab, to show the affinity gap does not translate
  BEV_M6_PRN  = function() run_reactive("crvo_nonisch", "bevacizumab",
                                       mode = "prn", load_n = 6),
  # BALATON / COMINO: faricimab, Ang-2 arm stabilises the junction
  FAR_TAE     = function() run_reactive("crvo_nonisch", "faricimab",
                                       mode = "tae", load_n = 6, imax = 112),
  BRO_TAE     = function() run_reactive("crvo_nonisch", "brolucizumab",
                                       mode = "tae", load_n = 6),
  AFL8_TAE    = function() run_reactive("crvo_nonisch", "aflibercept8",
                                       mode = "tae", load_n = 3, imax = 140),

  # --- deferral: the permanent price of waiting ----------------------------
  DEFER3      = function() run_reactive("crvo_nonisch", "aflibercept",
                                       mode = "tae", load_n = 6,
                                       start = BASE_START + 3 * 30.4),
  DEFER6      = function() run_reactive("crvo_nonisch", "aflibercept",
                                       mode = "tae", load_n = 6,
                                       start = BASE_START + 6 * 30.4),
  DEFER12     = function() run_reactive("crvo_nonisch", "aflibercept",
                                       mode = "tae", load_n = 6,
                                       start = BASE_START + 12 * 30.4),

  # --- steroid --------------------------------------------------------------
  DEX_PHAKIC  = function() run_case("crvo_nonisch", dex_times =
                                      BASE_START + 182.5 * (0:5)),
  DEX_PSEUDO  = function() run_case("crvo_nonisch", dex_times =
                                      BASE_START + 182.5 * (0:5),
                                    extra = list(KCATDEX = 0, ECAT = 0)),
  AFL_PLUS_DEX = function() run_reactive("crvo_nonisch", "aflibercept",
                                        mode = "tae", load_n = 6,
                                        dex_times = c(BASE_START + 91,
                                                      BASE_START + 274),
                                        extra = list(KCATDEX = 0, ECAT = 0)),

  # --- real-world undertreatment -------------------------------------------
  # Fight Retinal Blindness! registry-like intensity: 6 in y1, 3 in y2, 2 in y3
  REALWORLD   = function() run_case("crvo_nonisch", ivt_times = BASE_START +
                                      c(0, 28, 56, 112, 182, 258, 378, 478, 598,
                                        738, 908)),
  REALWORLD_LOW = function() run_case("crvo_nonisch", ivt_times =
                                        BASE_START + c(0, 28, 56, 140, 400, 560, 800)),

  # --- THE PRESSURE-ARM EXPERIMENT -----------------------------------------
  #  Infinite blockade of everything the drug binds, starting at presentation.
  #  Whatever oedema remains is the pressure arm, by construction.
  PERFECT_NONISCH = function() run_case("crvo_nonisch", ivt_times = BASE_START,
                                        perfect = TRUE),
  PERFECT_ISCH    = function() run_case("crvo_isch", ivt_times = BASE_START,
                                        perfect = TRUE),
  PERFECT_MILD    = function() run_case("crvo_mild", ivt_times = BASE_START,
                                        perfect = TRUE),
  PERFECT_BRVO    = function() run_case("brvo", ivt_times = BASE_START,
                                        perfect = TRUE),

  # --- ischaemic CRVO and neovascular glaucoma -----------------------------
  ISCH_AFL_CONT   = function() run_case("crvo_isch", ivt_times = q4w(38)),
  ISCH_AFL_STOP6  = function() run_case("crvo_isch", ivt_times = q4w(6)),
  ISCH_AFL_PRP    = function() run_case("crvo_isch", ivt_times = q4w(6),
                                        prp_time = BASE_START + 91),
  ISCH_AFL_TAE    = function() run_reactive("crvo_isch", "aflibercept",
                                            mode = "tae", load_n = 6),

  # --- BRVO -----------------------------------------------------------------
  BRVO_AFL    = function() run_reactive("brvo", "aflibercept",
                                        mode = "prn", load_n = 6),
  BRVO_RBZ    = function() run_reactive("brvo", "ranibizumab",
                                        mode = "prn", load_n = 6),

  # --- blood pressure: the only arrow into the pressure arm ----------------
  AFL_BP_LOW      = function() run_reactive("crvo_nonisch", "aflibercept",
                                            mode = "tae", load_n = 6,
                                            extra = list(PA = 32)),
  BP_LOW_ALONE    = function() run_case("crvo_nonisch", extra = list(PA = 32)),
  ISCH_AFL_BP_LOW = function() run_reactive("crvo_isch", "aflibercept",
                                            mode = "tae", load_n = 6,
                                            extra = list(PA = 32))
)

run_all <- function(which = names(SCENARIOS)) {
  setNames(lapply(which, function(n) {
    message("running ", n)
    out <- SCENARIOS[[n]]()
    out$scenario <- n
    out
  }), which)
}

# =============================================================================
#  7.  CLOSED-FORM ANALYSES
# =============================================================================

#  7a.  The critical capillary pressure and the critical venous resistance.
#       Pc* = Pt + sigma*(pi_c - pi_t).  Setting Pc = Pc* and solving
#       Pc = (Pa Rv + Pv Ra)/(Ra + Rv) for Rv gives Rv* = (Pc* - Pv) Ra /
#       (Pa - Pc*), which diverges as Pa approaches Pc* from above: if the
#       arterial inlet pressure is below Pc*, NO occlusion can create a
#       permanent oedema floor.
critical_pressure <- function(pt = 12, sigma = 0.92, pic = 25, pit = 4) {
  pt + sigma * (pic - pit)
}

critical_rv <- function(pcstar, pa, ra = 1.0, pv = 17) {
  ifelse(pa <= pcstar, NA_real_, (pcstar - pv) * ra / (pa - pcstar))
}

analyse_critical_pressure <- function() {
  pcs <- critical_pressure()
  tbl <- tibble(Pa = c(48, 44, 40, 36, 34, 32, 31.5, 31.0)) %>%
    mutate(Rv_star   = critical_rv(pcs, Pa),
           x_base    = Rv_star / 0.25,
           dPc_dPa   = Rv_star / (1 + Rv_star))
  chronic <- critical_pressure(sigma = 0.92 * 0.70^0.55,
                              pit = 4 + 21 * (1 - 0.70^0.55) * 0.75)
  list(Pc_star = pcs, Pc_star_chronic = chronic, table = tbl,
       Rv_star_intact = critical_rv(pcs, 40),
       Rv_star_chronic = critical_rv(chronic, 40))
}

#  7b.  Suppression duration.  With C_ret(t) = C0_eff exp(-k t) and a
#       suppression ratio R = 1 + C_ret/(alpha KD), the time for which R stays
#       above a target R* is
#            t_sup = (1/k) ln[ C0_eff / ((R*-1) alpha KD) ]
#       which is LOGARITHMIC in molar dose and affinity but LINEAR in the
#       half-life.  That is why a 94-fold affinity difference between
#       ranibizumab and aflibercept becomes only a ~4-fold difference in time.
suppression_duration <- function(drug, alpha = 120, fpen = 0.08, vvit = 4,
                                 target_ratio = 4.9) {
  d <- DRUGS[DRUGS$drug == drug, ]
  c0  <- molar_dose(drug) / vvit * 1e6 * fpen            # pM at the effect site
  thr <- (target_ratio - 1) * alpha * d$kdv
  kel <- log(2) / d$tvit
  tibble(drug = drug, label = d$label, nmol = molar_dose(drug),
         C0_vit_uM = molar_dose(drug) / vvit, C0_eff_pM = c0,
         KD_pM = d$kdv, threshold_pM = thr, t_half = d$tvit,
         t_sup = ifelse(c0 <= thr, 0, log(c0 / thr) / kel),
         reservoir_term = log(c0), affinity_term = -log(thr),
         halflife_mult = d$tvit / log(2))
}

analyse_suppression <- function() {
  bind_rows(lapply(DRUGS$drug, suppression_duration))
}

#  7c.  Days suppressed versus days dry.  A ratio below 1 means the macula was
#       WET while VEGF was still adequately suppressed.  That gap is the
#       pressure arm; switching agent cannot close it.
analyse_suppressed_vs_dry <- function(res) {
  bind_rows(lapply(names(res), function(n) {
    d <- res[[n]]
    tibble(scenario = n,
           days_suppressed = tail(d$TSUP, 1),
           days_dry = sum(d$CSTt < 300, na.rm = TRUE),
           ratio = ifelse(tail(d$TSUP, 1) > 0,
                          sum(d$CSTt < 300) / tail(d$TSUP, 1), NA))
  }))
}

#  7d.  When does the eye cross below its own critical pressure?
analyse_crossing <- function() {
  bind_rows(lapply(names(PHENOTYPES), function(ph) {
    d <- run_case(ph)
    d <- d[d$time >= BASE_START, ]
    below <- which(d$PCt < d$PCCRIt)
    tibble(phenotype = ph,
           Pc_d42 = d$PCt[1], Pc_crit_d42 = d$PCCRIt[1],
           first_below = if (length(below)) d$time[below[1]] else NA_real_,
           days_above = sum(d$PCt >= d$PCCRIt))
  }))
}

#  7e.  The pressure-arm decomposition.  CST under a real regimen, CST under
#       infinite blockade, and the dry baseline: the three numbers that split
#       residual thickness into a pharmacological miss and a pressure floor.
analyse_pressure_arm <- function(res) {
  pick <- function(n, day) {
    d <- res[[n]]; d$CSTt[which.min(abs(d$time - day))]
  }
  tibble(
    phenotype = c("non-ischaemic CRVO", "ischaemic CRVO"),
    real_regimen_m36 = c(pick("AFL_M6_Q8", 1095), pick("ISCH_AFL_TAE", 1095)),
    infinite_block_m36 = c(pick("PERFECT_NONISCH", 1095), pick("PERFECT_ISCH", 1095)),
    dry_baseline = 250) %>%
    mutate(pharmacological_miss = real_regimen_m36 - infinite_block_m36,
           pressure_floor = infinite_block_m36 - dry_baseline)
}

#  7f.  Delay-response: permanent letters forfeited per month of deferral.
analyse_delay <- function(months = c(0, 0.5, 1, 2, 3, 4, 6, 9, 12)) {
  base <- NULL
  bind_rows(lapply(months, function(mo) {
    d <- run_reactive("crvo_nonisch", "aflibercept", mode = "tae",
                      load_n = 6, start = BASE_START + mo * 30.4)
    m36 <- d$BCVA[which.min(abs(d$time - 1095))]
    if (is.null(base)) base <<- m36
    tibble(delay_months = mo,
           bcva_at_start = d$BCVA[which.min(abs(d$time - (BASE_START + mo * 30.4)))],
           peak_bcva = max(d$BCVA), bcva_m36 = m36,
           ez_m36 = d$EZ[which.min(abs(d$time - 1095))],
           permanent_loss = m36 - base, n_inj = attr(d, "n_inj"))
  }))
}

# =============================================================================
#  8.  VIRTUAL POPULATION
# =============================================================================
virtual_population <- function(n = 400, seed = 20260805, end = 730) {
  set.seed(seed)
  lnorm <- function(median, cv) median * exp(rnorm(1, 0, sqrt(log(1 + cv^2))))
  bind_rows(lapply(seq_len(n), function(i) {
    occ0 <- min(40, max(1.5, lnorm(11.0, 0.55)))
    pa   <- rnorm(1, 40, 4.5)
    coll <- lnorm(1.3, 0.60)
    kez  <- lnorm(0.0018, 0.40)
    vg   <- lnorm(12.0, 0.35)
    ores <- min(0.95, max(0.05, rnorm(1, 0.55, 0.16)))
    delay <- max(7, rnorm(1, 42, 25))
    nmax  <- as.integer(min(24, max(4, rnorm(1, 11, 4))))
    d <- run_reactive("crvo_nonisch", "aflibercept", mode = "tae",
                      load_n = min(6, nmax), start = delay, max_inj = nmax,
                      extra = list(OCC0 = occ0, PA = pa, COLLMAX = coll,
                                   KEZ = kez, EHIFV = vg, OCCRES = ores),
                      end = end)
    i0 <- which.min(abs(d$time - delay)); iE <- nrow(d)
    tibble(id = i, occ0 = occ0, occ_res = ores, pa = pa, collmax = coll,
           delay = delay, n_inj = attr(d, "n_inj"),
           pc = d$PCt[iE], pc_crit = d$PCCRIt[iE],
           cst_base = d$CSTt[i0], cst = d$CSTt[iE],
           bcva_base = d$BCVA[i0], bcva = d$BCVA[iE],
           gain = d$BCVA[iE] - d$BCVA[i0], ez = d$EZ[iE],
           np = d$NP[iE], nvi = d$NVI[iE],
           dry = as.integer(d$DRY[iE] > 0.5),
           floor = as.integer(d$FLOOR[iE] > 0.5))
  }))
}

summarise_population <- function(pop) {
  pop %>% summarise(
    n = n(),
    cst_base = mean(cst_base), bcva_base = mean(bcva_base),
    cst_m24 = mean(cst), gain = mean(gain), n_inj = mean(n_inj),
    pct_dry = 100 * mean(dry), pct_floor = 100 * mean(floor),
    pct_floor_wet = 100 * mean(floor == 1 & dry == 0),
    pct_nofloor_wet = 100 * mean(floor == 0 & dry == 0),
    pct_gain15 = 100 * mean(gain >= 15), pct_lose15 = 100 * mean(gain <= -15),
    pct_nvi = 100 * mean(nvi > 0.5),
    r_pc_cst = cor(pc, cst), r_inj_cst = cor(n_inj, cst))
}

# =============================================================================
#  9.  PLOTS FOR THE FOUR CLAIMS
# =============================================================================

# Claim 1: the pressure arm is a floor, not a slope.
plot_claim1 <- function(res) {
  bind_rows(
    res$PERFECT_NONISCH %>% mutate(arm = "non-ischaemic, infinite blockade"),
    res$PERFECT_ISCH    %>% mutate(arm = "ischaemic, infinite blockade"),
    res$AFL_M6_Q8       %>% mutate(arm = "non-ischaemic, aflibercept T&E"),
    res$ISCH_AFL_TAE    %>% mutate(arm = "ischaemic, aflibercept T&E")) %>%
    ggplot(aes(time / 30.4, CSTt, colour = arm)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = c(250, 300), linetype = c("solid", "dashed"),
               colour = "grey50") +
    labs(x = "months from occlusion", y = "central subfield thickness (um)",
         title = "Infinite VEGF blockade does not dry an eye above Pc*",
         subtitle = "the gap between each solid line and 250 um is the pressure arm") +
    theme_bw() + theme(legend.position = "bottom")
}

# Claim 2: Pc and Pc* cross, or they do not.
plot_claim2 <- function() {
  bind_rows(lapply(names(PHENOTYPES), function(ph)
    run_case(ph) %>% mutate(phenotype = ph))) %>%
    select(time, phenotype, PCt, PCCRIt) %>%
    pivot_longer(c(PCt, PCCRIt)) %>%
    ggplot(aes(time / 30.4, value, colour = name)) +
    geom_line(linewidth = 0.8) + facet_wrap(~phenotype) +
    labs(x = "months from occlusion", y = "mmHg",
         title = "The crossing of Pc below Pc* is the clinical course",
         subtitle = "non-ischaemic CRVO crosses as collaterals mature; ischaemic CRVO never does") +
    theme_bw() + theme(legend.position = "bottom")
}

# Claim 3: suppression duration is logarithmic in dose, linear in half-life.
plot_claim3 <- function() {
  analyse_suppression() %>%
    ggplot(aes(reorder(label, t_sup), t_sup, fill = log10(KD_pM))) +
    geom_col() + coord_flip() +
    labs(x = NULL, y = "days of adequate VEGF suppression",
         title = "A 94-fold affinity gap is a 4-fold time gap",
         fill = "log10 KD (pM)") +
    theme_bw()
}

# Claim 4: vision is an integral, so delay is priced permanently.
plot_claim4 <- function(res) {
  bind_rows(
    res$AFL_M6_Q8 %>% mutate(arm = "immediate"),
    res$DEFER3    %>% mutate(arm = "3-month deferral"),
    res$DEFER6    %>% mutate(arm = "6-month deferral"),
    res$DEFER12   %>% mutate(arm = "12-month deferral")) %>%
    select(time, arm, BCVA, EZ) %>% pivot_longer(c(BCVA, EZ)) %>%
    ggplot(aes(time / 30.4, value, colour = arm)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~name, scales = "free_y") +
    labs(x = "months from occlusion", y = NULL,
         title = "Deferral is paid in ellipsoid zone, and the ellipsoid zone does not come back") +
    theme_bw() + theme(legend.position = "bottom")
}

# =============================================================================
# 10.  DEMONSTRATION
# =============================================================================
if (interactive()) {
  print(analyse_critical_pressure())
  print(analyse_suppression(), width = Inf)
  print(analyse_crossing())

  res <- run_all()
  summary_tbl <- bind_rows(lapply(names(res), function(n) {
    d <- res[[n]]
    pick <- function(day, col) d[[col]][which.min(abs(d$time - day))]
    tibble(scenario = n,
           cst_d42 = pick(42, "CSTt"), bcva_d42 = pick(42, "BCVA"),
           cst_m6 = pick(224, "CSTt"), bcva_m6 = pick(224, "BCVA"),
           cst_m36 = pick(1095, "CSTt"), bcva_m36 = pick(1095, "BCVA"),
           ez_m36 = pick(1095, "EZ"), np_m36 = pick(1095, "NP"),
           nvi_m36 = pick(1095, "NVI"), iop_m36 = pick(1095, "IOP"),
           n_inj = if (is.null(attr(d, "n_inj"))) NA else attr(d, "n_inj"))
  }))
  print(summary_tbl, n = 40, width = Inf)

  print(analyse_pressure_arm(res))
  print(analyse_suppressed_vs_dry(res))
  print(analyse_delay())

  pop <- virtual_population(n = 100)
  print(summarise_population(pop), width = Inf)

  print(plot_claim1(res)); print(plot_claim2())
  print(plot_claim3());    print(plot_claim4(res))
}
