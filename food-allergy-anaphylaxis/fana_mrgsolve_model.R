## =====================================================================
##  IgE-mediated food allergy and anaphylaxis (peanut as the exemplar)
##  Quantitative Systems Pharmacology model — mrgsolve implementation
## ---------------------------------------------------------------------
##  Companion files:
##     fana_qsp_model.dot / .svg / .png    mechanistic map (20 modules)
##     fana_references.md                  literature
##     fana_shiny_app.R                    interactive dashboard
##     README.md                           narrative, results, caveats
##
##  WHAT THIS MODEL IS FOR
##  ----------------------
##  Food allergy has, unusually among chronic diseases, one numerical
##  endpoint that everybody agrees on: the ELICITING DOSE (ED) — the
##  milligrams of allergen protein that provoke an objective reaction.
##  Labelling law is written from it (VITAL reference doses), trials are
##  powered on it (tolerating a single 600 mg dose of peanut protein),
##  and patients live inside it.
##
##  So this model is built to compute ONE thing well: where the ED is,
##  and what moves it. Everything else exists to support or check that
##  computation.
##
##  FIVE STRUCTURAL COMMITMENTS
##  ---------------------------
##  S1. THE SPECIFICITY ENTERS TWICE (THE SQUARE LAW).
##      FcepsilonRI binds the Fc of IgE, not its paratope, so the mast
##      cell surface samples the serum IgE pool WITHOUT REGARD TO
##      SPECIFICITY: the surface carries specific IgE at the same
##      fraction f = sIgE/total IgE as the serum does. Degranulation
##      needs ONE allergen molecule to BRIDGE TWO adjacent receptors, so
##      the density of bridgeable pairs goes as
##
##            X  ~  (rho * L * f)^2  *  [allergen]/(K + [allergen])
##
##      rho = receptor density, L = fractional occupancy. The exponent 2
##      is not fitted; it is combinatorics. It is also the most
##      consequential number in the file: anything acting on the SURFACE
##      buys twice its log-effect on the ED, anything acting on the
##      ALLERGEN buys exactly its log-effect.
##
##  S2. RECEPTOR DENSITY IS A SLOW STATE, NOT A CONSTANT.
##      Free IgE stabilises FcepsilonRI; remove the ligand and the
##      receptor is degraded. rho falls with a time constant of ~3 weeks
##      once free IgE is suppressed, giving anti-IgE a FAST arm
##      (occupancy, days) and a SLOW arm (density, months) that multiply
##      INSIDE the square. Hence the week-16 challenge (diagnostic D4).
##
##  S3. IMMUNOTHERAPY AND ANTI-IgE ACT AT DIFFERENT POINTS IN SERIES.
##      IgG4 intercepts allergen in the interstitium BEFORE the surface
##      (it divides [allergen]: LINEAR in the ED). Anti-IgE changes the
##      surface (QUADRATIC in the ED). Being in series, they MULTIPLY.
##      That is the model's main prospective claim (diagnostic D7).
##
##  S4. THE CIRCULATION HAS A RESERVE, SO SEVERITY IS A CLIFF.
##      About a third of plasma volume can move into the interstitium
##      before mean arterial pressure changes at all. MAP is therefore
##      flat-then-collapsing (Hill exponent 6 on fractional plasma
##      volume loss), not a graded function of mediator load. Grading
##      scales look ordinal because a reserve is being spent.
##
##  S5. THE MEASURED THRESHOLD IS A RANDOM VARIABLE.
##      Repeat challenges in one patient scatter by ~0.4-0.5 log10 and
##      most of that variance is cofactor-driven — exercise, NSAID,
##      alcohol, fever — acting on allergen DELIVERY. Cofactors here
##      multiply absorption and releasability and are deliberately NOT
##      allowed to touch sIgE, because a patient's antibodies do not
##      change between Tuesday and Wednesday but their threshold does.
##
##  WHAT THE MODEL DELIBERATELY DOES NOT DO
##  ---------------------------------------
##    * non-IgE food allergy (FPIES, EoE, coeliac) — different biology
##    * alpha-gal syndrome — IgE-mediated, but the kinetics are those of
##      chylomicron trafficking, not of protein absorption
##    * MRGPRX2 / idiopathic anaphylaxis / clonal mast cell disease as a
##      PRIMARY driver (hereditary alpha-tryptasaemia enters only as a
##      severity modifier, through MCBURDEN)
##    * the accuracy of self-reported food allergy, which is the largest
##      single source of error in the real clinic
##
##  UNITS
##  -----
##    time            hours (the file spans 2 minutes to 6 years)
##    IgE, IgG4       nM   (1 IU/mL IgE = 2.4 ng/mL = 0.01263 nM)
##    omalizumab      nmol in the compartment, nM in concentration terms
##    allergen        nmol / nM of effector allergen (Ara h 2-equivalent)
##    allergen DOSE   mg of PEANUT PROTEIN (the clinical unit)
##    adrenaline      ng in the compartment, ng/mL in concentration
##    volume L, pressure mmHg, FEV1 % predicted
##
##  Requires: mrgsolve (>= 1.0), dplyr, tidyr.
## =====================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
})

code <- '
$PROB IgE food allergy / anaphylaxis — surface square law, IgG4 interception, adrenaline rescue

$PARAM @annotated
// ---- patient descriptors --------------------------------------------
WT        :  70     : Body weight (kg)
SIGE0     :  40     : Baseline allergen-specific IgE (kU_A/L = IU/mL)
TIGE0     : 300     : Baseline TOTAL serum IgE (IU/mL)
MCBURDEN  :   1     : Mast cell burden multiplier (1 normal, 2-3 HaT/mastocytosis)
RELEASE   :   1     : Intrinsic releasability (0.2 nonreleaser .. 1.5)
ASTHMA    :   0     : Uncontrolled asthma flag (0/1) — SEVERITY not threshold modifier
BBLOCK    :   0     : Beta-blocker / ACEi flag (0/1)
PAFAH     :   1     : PAF-acetylhydrolase activity (1 normal, <1 = fatal phenotype)
AGEY      :  10     : Age (years) — gates the slow tolerance arm of immunotherapy

// ---- allergen delivery ------------------------------------------------
FARAH2    :   0.10  : Fraction of peanut protein that is Ara h 2/6 (effector allergen)
MWALG     :  17.0   : Effector allergen MW (kDa)
KA_GUT    :   1.5   : Absorption rate constant of intact allergen (1/h)
KDEG_GUT  :   2.5   : Luminal proteolysis rate constant (1/h)
FINTACT   :   0.0022: Fraction of absorbed allergen reaching interstitium INTACT
VDIST     :  15.0   : Allergen distribution volume, plasma+interstitium (L)
KEL_ALG   :   1.2   : Allergen elimination from interstitium (1/h)
FATMATRIX :   1     : Food matrix factor on absorption rate (1 aqueous, 0.5 high fat)

// ---- cofactors: they multiply DELIVERY, never immunology (S5) ---------
COF_EX    :   0     : Exercise flag (0/1)
COF_NSAID :   0     : NSAID / aspirin flag (0/1)
COF_ETOH  :   0     : Alcohol flag (0/1)
COF_INF   :   0     : Acute infection / fever flag (0/1)
COF_PPI   :   0     : Proton pump inhibitor flag (0/1)
FEX       :   3.0   : Exercise multiplier on intact-protein delivery
FNSAID    :   3.5   : NSAID multiplier on intact-protein delivery
FETOH     :   2.0   : Alcohol multiplier
FINF      :   1.8   : Infection multiplier on releasability
FPPI      :   2.2   : PPI multiplier on intact-protein survival
COF_CAP   :  12     : Cap on the combined delivery multiplier

// ---- the surface (S1, S2) ---------------------------------------------
KD_FCERI  :   0.10  : IgE-FcepsilonRI dissociation constant (nM)
RHO_FLOOR :   0.18  : Irreducible tissue mast cell FcepsilonRI fraction
KD_RHO    :   0.35  : Free IgE giving half-maximal receptor stabilisation (nM)
KOUT_RHO  :   0.00138: FcepsilonRI turnover (1/h; t1/2 ~21 d) — the SLOW arm
SURFREF   :   0.1209: Reference surface index rho*L*f of the calibration patient
HXL       :   2.0   : THE SQUARE LAW EXPONENT — read the README before changing

// ---- degranulation transfer function ------------------------------------
KA_ALG    :   1.0   : Allergen concentration half-saturating surface IgE (nM)
XL50      :   0.0122: Cross-link index giving half-maximal degranulation
HDG       :   1.7   : Hill exponent of the degranulation transfer function
KGRAN_DEP :   1.2   : Granule depletion rate per unit degranulation (1/h)
KGRAN_REC :   0.06  : Granule replenishment rate constant (1/h)
KHOOK     : 900     : Allergen concentration at which monovalent blockade starts (nM)

// ---- mediators ------------------------------------------------------------
KEL_HIST  :  20.8   : Histamine elimination (1/h; t1/2 2 min)
KEL_TRYP  :   0.347 : Tryptase elimination (1/h; t1/2 2 h)
KEL_PAF   :  14.0   : PAF elimination (1/h), scaled by PAFAH
KEL_CLT   :   1.4   : Cysteinyl leukotriene elimination (1/h)
FHIST     :   1.00  : Histamine pool at full degranulation (1 = full anaphylactic level)
FPAF      :   0.55  : PAF pool at full degranulation
FCLT      :   0.75  : cysLT pool at full degranulation
KTRYP_GEN :  60.0   : Tryptase appearance rate at full degranulation (ng/mL/h)
TRYP0     :   5.0   : Baseline serum tryptase (ng/mL)
KLATE     :   0.055 : Late-phase cell recruitment rate (1/h)
KELLATE   :   0.045 : Late-phase resolution rate (1/h)
FBIPH     :   0.30  : Coupling of late-phase cells back to mediator release

// ---- cardiovascular (S4) -----------------------------------------------------
PV0       :   3.0   : Baseline plasma volume (L)
KLEAK     :   1.05  : Capillary leak rate constant (1/h per unit permeability)
KRET      :   0.35  : Interstitial fluid return rate constant (1/h)
WH_LEAK   :   1.0   : Weight of histamine on permeability
WP_LEAK   :   2.6   : Weight of PAF on permeability (the dominant term)
WL_LEAK   :   1.1   : Weight of cysLT on permeability
FRAC50    :   0.36  : Fractional plasma volume loss at half-maximal MAP fall
HFRAC     :   7.0   : Hill exponent of the MAP cliff
WA1_RET   :   2.0   : alpha1 acceleration of interstitial fluid reabsorption
MAP0      :  88     : Baseline mean arterial pressure (mmHg)
KVASO     :   0.30  : Vasodilator potency of the mediator pool on SVR
CRIT_MAP  :  55     : MAP defining circulatory failure (mmHg)

// ---- respiratory ---------------------------------------------------------------
KON_BRO   :   2.6   : Bronchoconstriction onset rate (1/h per unit cysLT)
KOFF_BRO  :   1.1   : Bronchoconstriction offset rate (1/h)
WH_BRO    :   0.30  : Weight of histamine on bronchoconstriction
FEV0      : 100     : Baseline FEV1 (% predicted)
AMP_ASTHMA:   1.7   : Amplification of flow loss by pre-existing asthma

// ---- skin / GI --------------------------------------------------------------------
KON_URT   :   9.0   : Urticaria onset rate (1/h)
KOFF_URT  :   0.55  : Urticaria offset rate (1/h)
KON_GI    :   7.0   : GI symptom onset rate (1/h)
KOFF_GI   :   0.90  : GI symptom offset rate (1/h)

// ---- IgE / B cell compartment ---------------------------------------------------
KEL_IGE   :   0.0144 : Free IgE elimination (1/h; t1/2 ~2 d)
KEL_CPX   :   0.0016 : Omalizumab:IgE complex elimination (1/h; t1/2 ~18 d)
KDEATH_PC :   0.00048: IgE plasma cell death rate (1/h; t1/2 ~60 d)
KTH2_PC   :   0.00048: Th2-driven IgE plasma cell recruitment (1/h)
KBMEM     :   0.00035: Memory B recall into the IgE plasma cell pool (1/h)
KTH2_ON   :   0.0016 : Th2 relaxation rate towards its antigen set point (1/h)
ATH2      :   1.2    : Maximal antigen-driven Th2 expansion above baseline
KA_IGE    :   3.0    : Antigen sensitivity of the IgE arm (deliberately LOW)
KBM_ON    :   0.0004 : Memory B relaxation rate towards its set point (1/h)
ABM       :   1.0    : Maximal antigen-driven memory B expansion above baseline
INITMODE  :   1      : 1 = derive initial state from descriptors; 0 = use supplied init()
KPCNS     :   0.0008 : Non-specific IgE plasma cell turnover (1/h; t1/2 ~36 d)
NATRES    :   0.0    : Natural-resolution drive (0 = peanut, ~0.9 = milk/egg)

// ---- IgG4 / anergy / Treg — the immunotherapy arm (S3) --------------------------
KPC_G4    :   0.00048: IgG4 plasma cell induction per unit chronic antigen (1/h)
KDPC_G4   :   0.00048: IgG4 plasma cell death rate (1/h; t1/2 ~60 d)
KSEC_G4   :   0.55  : IgG4 secretion per plasma cell unit (nM/h)
KEL_G4    :   0.00138: IgG4 elimination (1/h; t1/2 ~21 d)
KI_G4     :  25.0   : IgG4 concentration intercepting half the allergen (nM)
KON_ANG   :  11.0   : Mast cell anergy induction (1/h per nM antigen)
KOFF_ANG  :   0.0096: Anergy loss (1/h; t1/2 ~3 d) — why missed doses matter
ANG_MAX   :   0.72  : Maximum attainable anergy
KTREG_ON  :   0.00016: Treg induction per unit chronic antigen (1/h)
KTREG_OFF :   0.00004: Treg decay (1/h) — the SLOW arm, years
TREGMAX   :   4.0   : Treg carrying capacity (relative)
AGE50_TOL :   4.0   : Age (y) at which the slow tolerance arm is half-closed
HAGE      :   3.0   : Steepness of the age gate
WTREG_TH2 :   0.85  : Treg suppression weight on Th2
WTREG_G4  :   0.50  : Treg (IL-10) boost of the IgG4 arm

// ---- omalizumab ------------------------------------------------------------------
F_OMA     :   0.62  : Omalizumab subcutaneous bioavailability
KA_OMA    :   0.0115: Subcutaneous absorption rate constant (1/h)
CL_OMA    :   0.0070: Omalizumab clearance of FREE drug (L/h at 70 kg)
CLC_OMA   :   0.0175: Clearance of the omalizumab:IgE complex (L/h at 70 kg)
VC_OMA    :   5.46  : Central volume (L at 70 kg)
Q_OMA     :   0.012 : Intercompartmental clearance (L/h)
VP_OMA    :   3.20  : Peripheral volume (L)
KD_OMA    :   1.0   : Omalizumab-IgE dissociation constant (nM)
NVAL_OMA  :   1.7   : Effective omalizumab valence for IgE (1:1 / 2:1 mixture)

// ---- dupilumab --------------------------------------------------------------------
F_DUP     :   0.64  : Dupilumab subcutaneous bioavailability
KA_DUP    :   0.010 : Absorption rate constant (1/h)
CL_DUP    :   0.0090: Clearance (L/h at 70 kg)
VC_DUP    :   4.80  : Central volume (L at 70 kg)
IC50_DUP  :  12.0   : Concentration for 50% IL-4Ra blockade (nM)
EMAX_DUP  :   0.80  : Maximal suppression of IgE class switching

// ---- adrenaline --------------------------------------------------------------------
KA_EPI_IM :   3.5   : IM absorption rate constant, anterolateral thigh (1/h)
KA_EPI_IN :   1.9   : Intranasal absorption rate constant (1/h)
KEL_EPI   :  16.6   : Adrenaline elimination (1/h; t1/2 ~2.5 min)
V_EPI     : 100     : Adrenaline distribution volume (L)
EC50_A1   :   0.42  : alpha1 EC50 (ng/mL)
EC50_B2   :   0.16  : beta2 EC50 (ng/mL)
EMAX_A1   :   0.80  : Maximal alpha1 reversal of leak and vasodilation
EMAX_B2   :   0.75  : Maximal beta2 bronchodilation / mast cell stabilisation
EPI_ENDO  :   0.05  : Endogenous adrenaline baseline (ng/mL)

// ---- H1 antihistamine (cetirizine) ---------------------------------------------------
KA_CET    :   1.3   : Absorption rate constant (1/h)
CL_CET    :   3.3   : Clearance (L/h)
VC_CET    :  35.0   : Central volume (L)
IC50_CET  :  80.0   : Concentration for 50% H1 blockade (ng/mL)
EMAX_CET  :   0.85  : Maximal H1 blockade — of the SKIN endpoint ONLY

// ---- other interventions ---------------------------------------------------------------
BTKI      :   0     : BTK inhibitor (blocks signal downstream of cross-linking)
EMAX_BTKI :   0.90  : Maximal suppression of the degranulation transfer function
ANTIKIT   :   0     : Anti-KIT mast cell depletion (0/1)
EMAX_KIT  :   0.85  : Maximal mast cell depletion
MONTELUK  :   0     : Leukotriene receptor antagonist (0/1)
EMAX_MONT :   0.55  : Maximal cysLT receptor blockade
IVFLUID   :   0     : IV crystalloid infusion rate (L/h)
SUPINE    :   0     : Supine with legs elevated (0/1) — autotransfusion

// ---- immunotherapy switch ----------------------------------------------------------------
OITDOSE   :   0     : Daily oral immunotherapy maintenance dose (mg peanut protein)
OITON     :   0     : Immunotherapy switch (0/1)

$CMT @annotated
ALGGUT  : Luminal allergen (nmol)
ALGSYS  : Interstitial effector allergen (nmol)
HIST    : Histamine pool (1 = full anaphylactic level)
TRYP    : Serum tryptase above baseline (ng/mL)
PAF     : PAF pool (relative)
CYSLT   : Cysteinyl leukotriene pool (relative)
LATEC   : Late-phase recruited cells (relative)
GRAN    : Mast cell granule content (1 = full)
ANERG   : Mast cell desensitisation / anergy (0-1)
RHO     : FcepsilonRI surface density (relative)
LEAKV   : Volume shifted from plasma to interstitium (L)
PV      : Plasma volume (L)
BRO     : Bronchoconstriction state (0-1)
URT     : Urticaria / skin state (0-1)
GISY    : Gastrointestinal symptom state (0-1)
AUCLEAK : Cumulative plasma volume deficit (L*h) — the harm integral
TH2     : Allergen-specific Th2 / Tfh13 pool (relative)
PCIGE   : Specific IgE plasma cell pool (relative)
PCNS    : Non-specific IgE plasma cell pool (relative)
BMEM    : IgG1+ memory B reservoir (relative)
SIGEt   : TOTAL specific IgE, free + omalizumab-bound (nM)
NIGEt   : TOTAL non-specific IgE, free + omalizumab-bound (nM)
PCG4    : IgG4 plasma cell pool (relative)
IGG4    : Allergen-specific IgG4 (nM)
TREG    : Allergen-specific Treg pool (relative)
OMASC   : Omalizumab subcutaneous depot (nmol)
OMAC    : Omalizumab central, free + IgE-bound (nmol)
OMAP    : Omalizumab peripheral (nmol)
DUPSC   : Dupilumab subcutaneous depot (nmol)
DUPC    : Dupilumab central (nmol)
EPIIM   : Adrenaline intramuscular depot (ng)
EPIIN   : Adrenaline intranasal depot (ng)
EPIC    : Adrenaline central (ng)
CETD    : Cetirizine gut depot (ng)
CETC    : Cetirizine central (ng)

$GLOBAL
#define IUML_TO_NM (0.01263)      // 1 IU/mL IgE = 2.4 ng/mL / 190 kDa

double posf(double x) { return (x > 1e-14) ? x : 1e-14; }

double hillf(double x, double x50, double h) {
  double xx = posf(x);
  double num = pow(xx, h);
  return num / (pow(x50, h) + num);
}

// closed-form 1:1 quasi-equilibrium binding; returns the complex
double bind11(double Dtot, double Ltot, double KD) {
  double s = Dtot + Ltot + KD;
  double disc = s*s - 4.0*Dtot*Ltot;
  if (disc < 0.0) disc = 0.0;
  return 0.5*(s - sqrt(disc));
}

$MAIN
// Everything below the patient descriptors is DERIVED, so changing
// SIGE0 / TIGE0 / AGEY re-generates a patient rather than switching
// between stored ones.
F_OMASC = F_OMA;
F_DUPSC = F_DUP;

double sige_ref = SIGE0 * IUML_TO_NM;
double tige_ref = TIGE0 * IUML_TO_NM;
double nige_ref = tige_ref - sige_ref;
if (nige_ref < 0.0) nige_ref = 0.0;

// INITMODE = 1 (default): derive the whole initial state from the
// patient descriptors. INITMODE = 0: leave the state alone, because the
// caller has supplied it with init() — this is what lets a challenge be
// run on a patient who is already 16 weeks into a drug course without
// re-simulating the course inside every bisection step.
if (INITMODE > 0.5) {
  SIGEt_0 = sige_ref;
  NIGEt_0 = nige_ref;
  RHO_0   = RHO_FLOOR + (1.0 - RHO_FLOOR) * tige_ref/(KD_RHO + tige_ref);
  GRAN_0  = 1.0;
  PV_0    = PV0;
  TH2_0   = 1.0;
  PCIGE_0 = 1.0;
  PCNS_0  = 1.0;
  BMEM_0  = 1.0;
  PCG4_0  = 0.02;
  IGG4_0  = KSEC_G4 * 0.02 / KEL_G4;  // ~8 nM = 1.2 ug/mL, a normal low level
  TREG_0  = 1.0;
}

$ODE
// =====================================================================
//  0.  COFACTORS — delivery multipliers, never immunology (S5)
// =====================================================================
double cof_deliv = 1.0;
if (COF_EX    > 0.5) cof_deliv *= FEX;
if (COF_NSAID > 0.5) cof_deliv *= FNSAID;
if (COF_ETOH  > 0.5) cof_deliv *= FETOH;
if (COF_PPI   > 0.5) cof_deliv *= FPPI;
if (cof_deliv > COF_CAP) cof_deliv = COF_CAP;

double releas = RELEASE;
if (COF_INF > 0.5) releas *= FINF;

// =====================================================================
//  1.  ALLERGEN KINETICS
//      Only the small INTACT fraction that survives digestion can
//      bridge two receptors. Everything else is nutrition.
// =====================================================================
double ka_eff = KA_GUT * FATMATRIX;
double f_int  = FINTACT * cof_deliv;
if (f_int > 0.5) f_int = 0.5;

dxdt_ALGGUT = -ka_eff*ALGGUT - KDEG_GUT*ALGGUT;
dxdt_ALGSYS =  ka_eff*ALGGUT*f_int - KEL_ALG*ALGSYS;

double Calg = posf(ALGSYS)/VDIST;                       // nM

// =====================================================================
//  2.  ANTIBODY COMPARTMENT AND OMALIZUMAB TMDD
//      Omalizumab cannot distinguish specific from non-specific IgE.
//      That is precisely why f is preserved (S1), and why TOTAL IgE
//      RISES 5-10x on treatment while FREE IgE falls.
// =====================================================================
double vco = VC_OMA * (WT/70.0);
double vpo = VP_OMA * (WT/70.0);
double clo = CL_OMA * pow(WT/70.0, 0.75);
double clc = CLC_OMA* pow(WT/70.0, 0.75);
double qo  = Q_OMA  * pow(WT/70.0, 0.75);

double omatot = posf(OMAC)/vco;                          // nM, free + bound
double igetot = posf(SIGEt) + posf(NIGEt);               // nM, free + bound

double sites  = NVAL_OMA * omatot;
double cpx    = bind11(sites, igetot, KD_OMA);           // nM IgE captured
double igefree = igetot - cpx;
if (igefree < 0.0) igefree = 0.0;
double frac_free = igefree/posf(igetot);

double omabound = cpx/NVAL_OMA;                          // nM drug engaged
double omafree  = omatot - omabound;
if (omafree < 0.0) omafree = 0.0;

double fspec = posf(SIGEt)/posf(igetot);                 // THE RATIO f

dxdt_OMASC = -KA_OMA*OMASC;
dxdt_OMAC  =  KA_OMA*OMASC
              - clo*omafree - clc*omabound
              - qo*omafree + qo*(posf(OMAP)/vpo);
dxdt_OMAP  =  qo*omafree - qo*(posf(OMAP)/vpo);

// =====================================================================
//  3.  THE SURFACE (S1, S2)
// =====================================================================
double Locc = igefree/(KD_FCERI + igefree);
double kit_eff = 1.0 - ((ANTIKIT > 0.5) ? EMAX_KIT : 0.0);
double mcpool  = MCBURDEN * kit_eff;

double SURF = RHO * Locc * fspec;
double surf_rel = SURF/SURFREF;

double rho_ss = RHO_FLOOR + (1.0 - RHO_FLOOR) * igefree/(KD_RHO + igefree);
dxdt_RHO = KOUT_RHO * (rho_ss - RHO);

// =====================================================================
//  4.  IgG4 INTERCEPTION — upstream of the surface, therefore LINEAR
// =====================================================================
double intercept = 1.0 + posf(IGG4)/KI_G4;
double Calg_eff  = Calg/intercept;

// =====================================================================
//  5.  THE ENGINE — cross-link density, QUADRATIC in surface IgE
// =====================================================================
double occ_alg = Calg_eff/(KA_ALG + Calg_eff);
double hookf   = 1.0/(1.0 + Calg_eff/KHOOK);
double XL = pow(posf(surf_rel), HXL) * occ_alg * hookf
            * (1.0 - ANERG) * releas * mcpool;

double btk_eff  = 1.0 - ((BTKI > 0.5) ? EMAX_BTKI : 0.0);
double activate = hillf(XL, XL50, HDG) * btk_eff;

// =====================================================================
//  6.  ADRENALINE AND ANTIHISTAMINE PD
// =====================================================================
double Cepi = posf(EPIC)/(V_EPI*1000.0) + EPI_ENDO;      // ng/mL
double bb   = (BBLOCK > 0.5) ? 0.35 : 1.0;
double E_a1 = EMAX_A1 * Cepi/(EC50_A1 + Cepi);
double E_b2 = EMAX_B2 * bb * Cepi/(EC50_B2 + Cepi);

double Ccet = posf(CETC)/(VC_CET*1000.0);                // ng/mL
double E_h1 = EMAX_CET * Ccet/(IC50_CET + Ccet);

dxdt_EPIIM = -KA_EPI_IM*EPIIM;
dxdt_EPIIN = -KA_EPI_IN*EPIIN;
dxdt_EPIC  =  KA_EPI_IM*EPIIM + KA_EPI_IN*EPIIN - KEL_EPI*EPIC;
dxdt_CETD  = -KA_CET*CETD;
dxdt_CETC  =  KA_CET*CETD - (CL_CET/VC_CET)*CETC;

// beta2 raises the calcium threshold: it stops release that has not
// yet happened. It does nothing about mediators already in the blood.
double degran = activate * GRAN * (1.0 - E_b2);

dxdt_GRAN = -KGRAN_DEP*degran + KGRAN_REC*(1.0 - GRAN);

// =====================================================================
//  7.  MEDIATORS
// =====================================================================
double late_boost = 1.0 + FBIPH*LATEC;
double kel_paf    = KEL_PAF*PAFAH;
double dg = degran*late_boost;

dxdt_HIST  = KEL_HIST*(FHIST*dg - HIST);
dxdt_PAF   = kel_paf *(FPAF *dg - PAF);
dxdt_CYSLT = KEL_CLT *(FCLT *dg - CYSLT);
dxdt_TRYP  = KTRYP_GEN*dg*MCBURDEN - KEL_TRYP*TRYP;
dxdt_LATEC = KLATE*degran - KELLATE*LATEC;

double mont_eff = 1.0 - ((MONTELUK > 0.5) ? EMAX_MONT : 0.0);
double cysltA = CYSLT*mont_eff;

// =====================================================================
//  8.  CARDIOVASCULAR — the reserve, then the cliff (S4)
// =====================================================================
double perm = WH_LEAK*HIST + WP_LEAK*PAF + WL_LEAK*cysltA;
double leak_in = KLEAK*perm*(1.0 - E_a1)*(PV/PV0);
if (leak_in < 0.0) leak_in = 0.0;

// alpha1 does two things to the volume balance and both matter: it
// closes the leak, and by lowering capillary hydrostatic pressure it
// speeds reabsorption of what has already gone.
double ret = KRET*(1.0 + WA1_RET*E_a1)*LEAKV;
dxdt_LEAKV = leak_in - ret;
dxdt_PV    = -leak_in + ret + IVFLUID;

double deficit = PV0 - PV;
dxdt_AUCLEAK = (deficit > 0.0) ? deficit : 0.0;

// =====================================================================
//  9.  RESPIRATORY / SKIN / GI
// =====================================================================
dxdt_BRO  = KON_BRO*(cysltA + WH_BRO*HIST)*(1.0 - BRO)
            - KOFF_BRO*BRO*(1.0 + 2.4*E_b2);
dxdt_URT  = KON_URT*HIST*(1.0 - URT) - KOFF_URT*URT*(1.0 + 1.6*E_a1);
dxdt_GISY = KON_GI*(HIST + 0.4*cysltA)*(1.0 - GISY) - KOFF_GI*GISY;

// =====================================================================
// 10.  CHRONIC IMMUNOLOGY
//      The chronic antigen signal is the DAILY-AVERAGED exposure: the
//      immunotherapy maintenance dose for a treated patient, and
//      essentially nothing for an avoiding one.
// =====================================================================
double oit_sig = (OITON > 0.5) ? OITDOSE/300.0 : 0.0;
double antigen_chronic = oit_sig + 0.02;

double Cdup = posf(DUPC)/(VC_DUP*(WT/70.0));
double dup_block = EMAX_DUP * Cdup/(IC50_DUP + Cdup);

dxdt_DUPSC = -KA_DUP*DUPSC;
dxdt_DUPC  =  KA_DUP*DUPSC - (CL_DUP*pow(WT/70.0,0.75)/(VC_DUP*(WT/70.0)))*DUPC;

double treg_ex = (TREG - 1.0 > 0.0) ? (TREG - 1.0) : 0.0;
double treg_supp = 1.0/(1.0 + WTREG_TH2*treg_ex);
double agegate = 1.0 - hillf(AGEY, AGE50_TOL, HAGE);

// Th2 and memory B are SET-POINT states, not exponential growers: the
// repertoire is antigen-driven but bounded. Writing them as unbounded
// proliferation makes any chronic-antigen scenario (i.e. all of
// immunotherapy) diverge, which is a modelling error and not biology.
double ant_ex = antigen_chronic - 0.02;
if (ant_ex < 0.0) ant_ex = 0.0;

// THE IgE ARM AND THE IgG4 ARM DO NOT SEE THE SAME ANTIGEN SIGNAL.
// That asymmetry is the immunological content of immunotherapy: oral
// antigen under tolerogenic conditions re-routes the response towards
// IgG4 and Treg rather than driving IgE in proportion. Empirically the
// sIgE rise on oral immunotherapy is transient and modest (~1.5-2x)
// while sIgG4 rises 10-100x, so the IgE arm is given a deliberately
// LOW antigen sensitivity (KA_IGE) and the IgG4 arm a linear one.
double ant_ige = ant_ex/(KA_IGE + ant_ex);

double th2_tgt = (1.0 + ATH2*ant_ige)
                 * treg_supp * (1.0 - 0.5*dup_block) / (1.0 + NATRES);
double bmem_tgt = (1.0 + ABM*ant_ige) / (1.0 + 2.0*NATRES);

dxdt_TH2   = KTH2_ON*(th2_tgt - TH2);
dxdt_BMEM  = KBM_ON*(bmem_tgt - BMEM);
dxdt_PCIGE = KTH2_PC*TH2*(1.0 - dup_block)
             + KBMEM*BMEM*ant_ige
             - KDEATH_PC*PCIGE*(1.0 + NATRES);
dxdt_PCNS  = KPCNS*(1.0 - dup_block) - KPCNS*PCNS;

// IgE turnover: free IgE is cleared fast, omalizumab-bound IgE slowly.
// This single line is why total IgE rises on anti-IgE therapy.
double kel_ige_eff = frac_free*KEL_IGE + (1.0 - frac_free)*KEL_CPX;
dxdt_SIGEt = KEL_IGE*sige_ref*PCIGE - kel_ige_eff*SIGEt;
dxdt_NIGEt = KEL_IGE*nige_ref*PCNS  - kel_ige_eff*NIGEt;

// IgG4 arm (fast, reversible) and Treg arm (slow, age-gated)
dxdt_PCG4 = KPC_G4*antigen_chronic*(1.0 + WTREG_G4*treg_ex) - KDPC_G4*PCG4;
dxdt_IGG4 = KSEC_G4*PCG4 - KEL_G4*IGG4;
dxdt_TREG = KTREG_ON*antigen_chronic*agegate*(TREGMAX - TREG)
            - KTREG_OFF*(TREG - 1.0);

// mast cell anergy: driven by the ACTUAL antigen the cell sees, lost
// with a 3-day half-life — the reason a missed week of dosing matters
dxdt_ANERG = KON_ANG*Calg*(ANG_MAX - ANERG) - KOFF_ANG*ANERG;

$TABLE
// ---- re-derive the observables (kept explicit rather than cached) -----
double vcoT = VC_OMA*(WT/70.0);
double omatotT = posf(OMAC)/vcoT;
double igetotT = posf(SIGEt) + posf(NIGEt);
double cpxT = bind11(NVAL_OMA*omatotT, igetotT, KD_OMA);
double igefreeT = igetotT - cpxT; if (igefreeT < 0.0) igefreeT = 0.0;
double fspecT = posf(SIGEt)/posf(igetotT);
double LoccT  = igefreeT/(KD_FCERI + igefreeT);
double SURFT  = RHO*LoccT*fspecT;
double CalgT  = posf(ALGSYS)/VDIST;
double interT = 1.0 + posf(IGG4)/KI_G4;
double CalgE  = CalgT/interT;
double occalgT= CalgE/(KA_ALG + CalgE);
double relT   = RELEASE*((COF_INF > 0.5) ? FINF : 1.0);
double mcT    = MCBURDEN*(1.0 - ((ANTIKIT > 0.5) ? EMAX_KIT : 0.0));
double XLT = pow(posf(SURFT/SURFREF), HXL)*occalgT
             *(1.0/(1.0 + CalgE/KHOOK))*(1.0 - ANERG)*relT*mcT;

double CepiT = posf(EPIC)/(V_EPI*1000.0) + EPI_ENDO;
double bbT   = (BBLOCK > 0.5) ? 0.35 : 1.0;
double Ea1T  = EMAX_A1*CepiT/(EC50_A1 + CepiT);
double CcetT = posf(CETC)/(VC_CET*1000.0);
double Eh1T  = EMAX_CET*CcetT/(IC50_CET + CcetT);
double montT = 1.0 - ((MONTELUK > 0.5) ? EMAX_MONT : 0.0);
double cysltT= CYSLT*montT;

double frac_lostT = (PV0 - PV)/PV0; if (frac_lostT < 0.0) frac_lostT = 0.0;
double reserveT = hillf(frac_lostT, FRAC50*bbT, HFRAC);
double vasoT = KVASO*(HIST + 0.8*PAF + 0.5*cysltT);
double toneT = (1.0 + 1.35*Ea1T)/(1.0 + vasoT);
double supT  = (SUPINE > 0.5) ? 0.13 : 0.0;
double MAPv  = MAP0*toneT*(1.0 - reserveT)*(1.0 + supT);

double ampT = 1.0 + ASTHMA*(AMP_ASTHMA - 1.0);
double FEV1v = FEV0*(1.0 - BRO*ampT);
if (FEV1v < 12.0) FEV1v = 12.0;

double URTobs = URT*(1.0 - Eh1T);

// PRACTALL-like severity grade, built from the physiology rather than
// asserted. Note that grade 1-2 are DOMINATED by skin, which is exactly
// the weakness the map annotates (cluster 10).
double SEVv = 0.0;
if (URTobs > 0.12 || GISY > 0.12) SEVv = 1.0;
if (URTobs > 0.45 || GISY > 0.45) SEVv = 2.0;
if (FEV1v < 0.80*FEV0 || MAPv < 0.86*MAP0) SEVv = 3.0;
if (FEV1v < 0.60*FEV0 || MAPv < 0.76*MAP0) SEVv = 4.0;
if (MAPv  < CRIT_MAP  || FEV1v < 0.40*FEV0) SEVv = 5.0;

capture TotIgE_IU  = igetotT/IUML_TO_NM;
capture FreeIgE_IU = igefreeT/IUML_TO_NM;
capture FreeIgE_ng = igefreeT*190.0;
capture sIgE_IU    = posf(SIGEt)/IUML_TO_NM;
capture Fspec      = fspecT;
capture Locc_occ   = LoccT;
capture SURFidx    = SURFT;
capture SURFrel    = SURFT/SURFREF;
capture XLidx      = XLT;
capture IgG4_ug    = posf(IGG4)*146.0/1000.0;
capture G4E_ratio  = (posf(IGG4)*146.0)/posf(posf(SIGEt)*190.0);
capture Interc     = interT;
capture Calg_nM    = CalgT;
capture Oma_ug     = omatotT/6.7114;
capture Cepi_ngmL  = CepiT;
capture Tryptase   = TRYP0 + TRYP;
capture MAP        = MAPv;
capture FEV1       = FEV1v;
capture URTskin    = URTobs;
capture SEV        = SEVv;
capture PVdef      = (PV0 - PV)/PV0;
capture Reserve    = reserveT;
capture MAPtone    = toneT;
capture MAPvol     = MAP0*(1.0 - reserveT);
'

mod <- mcode_cache("fana", code)

## =====================================================================
##  UNIT HELPERS
## =====================================================================
YR   <- 8766          # hours per year
WK   <- 168           # hours per week
DAY  <- 24

##' mg of peanut PROTEIN -> nmol of effector (Ara h 2-equivalent) allergen
mg_protein_to_nmol <- function(mg, f_arah2 = 0.10, mw_kda = 17) {
  mg * f_arah2 * 1e-3 / (mw_kda * 1000) * 1e9
}

##' The omalizumab label dosing table, as a function rather than a lookup.
##' Returns list(dose_mg, interval_h). NA dose = "do not dose" (outside
##' the licensed IgE/weight envelope) — a real and clinically important
##' outcome that the model reproduces as non-response.
omalizumab_dose <- function(wt, total_ige) {
  ## Approximation of the label table by the rule that actually underlies
  ## it: the drug must supply enough BINDING SITES for the IgE present,
  ## which is ~0.008 mg/kg per IU/mL per fortnight. Above the licensed
  ## envelope the dose is capped at 600 mg q2wk and the free-IgE target
  ## is simply missed — a real and clinically important outcome that the
  ## model reproduces as non-response rather than as an error.
  need2 <- 0.008 * wt * total_ige          # mg per 2 weeks, site-matched
  if (need2 <= 187.5) {
    list(dose_mg = min(600, max(75, ceiling(need2*2/75)*75)), interval_h = 4*WK)
  } else {
    list(dose_mg = min(600, max(75, ceiling(need2/75)*75)),   interval_h = 2*WK)
  }
}

## =====================================================================
##  EVENT BUILDERS
## =====================================================================
ev_allergen <- function(mg, time = 0, f_arah2 = 0.10) {
  ev(time = time, amt = mg_protein_to_nmol(mg, f_arah2), cmt = "ALGGUT")
}

ev_challenge_practall <- function(start = 0, gap = 0.5,
                                  doses = c(3, 10, 30, 100, 300, 600, 1000)) {
  ## PRACTALL escalation: semi-log steps every 30 min. Doses are SINGLE,
  ## not cumulative — the trial endpoint is a single tolerated dose.
  Reduce(`+`, lapply(seq_along(doses), function(i)
    ev_allergen(doses[i], time = start + (i - 1)*gap)))
}

ev_omalizumab <- function(dose_mg, interval_h, n, start = 0) {
  ev(time = start, amt = dose_mg*6.7114, cmt = "OMASC",
     ii = interval_h, addl = n - 1)
}

ev_dupilumab <- function(dose_mg = 300, interval_h = 2*WK, n = 12, start = 0) {
  ev(time = start, amt = dose_mg*6.8027, cmt = "DUPSC",
     ii = interval_h, addl = n - 1)
}

ev_epi_im <- function(mg = 0.3, time = 0, n = 1, ii = 0.0833) {
  ev(time = time, amt = mg*1e6, cmt = "EPIIM", ii = ii, addl = n - 1)
}

ev_epi_in <- function(mg = 2.0, time = 0) {
  ev(time = time, amt = mg*1e6, cmt = "EPIIN")
}

ev_cetirizine <- function(mg = 10, time = -1) {
  ev(time = time, amt = mg*1e6, cmt = "CETD")
}

##' Daily OIT dosing. The chronic antigen signal is carried by the
##' OITDOSE/OITON parameters; this adds the ACUTE daily spike so that
##' anergy and dosing reactions are driven by real exposure.
ev_oit_daily <- function(dose_mg, days, start = 0) {
  ev(time = start, amt = mg_protein_to_nmol(dose_mg), cmt = "ALGGUT",
     ii = 24, addl = days - 1)
}

## =====================================================================
##  PATIENT CONSTRUCTOR
## =====================================================================
##' The calibration patient: a 10-year-old, peanut allergic, sIgE
##' 40 kU/L on total IgE 300 IU/mL, no asthma, no cofactors. Every
##' scenario below is a modification of this one object.
patient <- function(...) {
  base <- list(WT = 35, AGEY = 10, SIGE0 = 40, TIGE0 = 300,
               MCBURDEN = 1, RELEASE = 1, ASTHMA = 0, BBLOCK = 0,
               PAFAH = 1, NATRES = 0)
  modifyList(base, list(...))
}

## =====================================================================
##  THE CENTRAL COMPUTATION: WHERE IS THE ELICITING DOSE?
## =====================================================================
##' find_ED(): bisect on single-dose allergen challenge to locate the
##' smallest dose producing an objective reaction (severity grade >= 2).
##'
##' This is the model's primary output. It is deliberately defined on
##' the OBSERVABLE (a graded clinical reaction) and not on any internal
##' variable, so that antihistamine masking, asthma amplification and
##' the MAP cliff all enter it the way they enter a real challenge.
find_ED <- function(par = patient(), init0 = NULL, grade = 2,
                    lo = 0.03, hi = 2e5, tol = 0.02, maxit = 40,
                    horizon = 6, delta = 0.02) {

  react <- function(dose) {
    m <- mod %>% param(par)
    if (!is.null(init0)) {
      ## INITMODE = 0 stops $MAIN from overwriting the supplied state
      m <- m %>% param(INITMODE = 0) %>% mrgsolve::init(init0)
    }
    out <- m %>%
      ev(ev_allergen(dose)) %>%
      mrgsim(end = horizon, delta = delta, atol = 1e-10, rtol = 1e-8) %>%
      as.data.frame()
    max(out$SEV, na.rm = TRUE) >= grade
  }

  ## An ED above ~10 g of protein is not reachable by eating; the number
  ## is reported anyway because the ratio is the object of interest, but
  ## it should be read as "no achievable dietary exposure".
  if (react(lo)) return(lo)
  if (!react(hi)) return(NA_real_)

  it <- 0
  while (log10(hi/lo) > tol && it < maxit) {
    mid <- sqrt(lo*hi)                    # bisect in log space
    if (react(mid)) hi <- mid else lo <- mid
    it <- it + 1
  }
  sqrt(lo*hi)
}

##' Steady-state / end-of-treatment initial conditions, so that a
##' challenge can be run on a treated patient without re-simulating the
##' whole treatment course inside the bisection.
state_at <- function(par, events = NULL, tend, delta = 24) {
  m <- mod %>% param(par)
  if (!is.null(events)) m <- m %>% ev(events)
  out <- m %>% mrgsim(end = tend, delta = delta, atol = 1e-10, rtol = 1e-8) %>%
    as.data.frame()
  last <- out[nrow(out), ]
  cn <- names(mrgsolve::init(mod))
  st <- as.list(last[, cn])
  ## the acute compartments must be reset: we are asking "what is this
  ## patient's threshold now", not "what happens if we dose mid-reaction"
  for (z in c("ALGGUT","ALGSYS","HIST","TRYP","PAF","CYSLT","LATEC",
              "LEAKV","BRO","URT","GISY","AUCLEAK","EPIIM","EPIIN",
              "EPIC","CETD","CETC")) st[[z]] <- 0
  st$GRAN <- 1
  st$PV   <- unname(param(mod)$PV0)
  st
}

## =====================================================================
##  SCENARIOS
## =====================================================================

##' SC1 — a DBPCFC in an untreated patient (PRACTALL escalation)
sc_dbpcfc <- function(par = patient()) {
  mod %>% param(par) %>%
    ev(ev_challenge_practall()) %>%
    mrgsim(end = 12, delta = 0.01) %>% as.data.frame()
}

##' SC2 — accidental ingestion of a full serving (2 g peanut protein),
##' with adrenaline given at a chosen delay. The comparison is the
##' HARM INTEGRAL (AUCLEAK), not the peak.
sc_accidental <- function(par = patient(), dose_mg = 2000,
                          epi_delay_min = 10, supine = 1,
                          n_epi = 3, epi_repeat_min = 8) {
  ## The realistic protocol is 0.3 mg IM REPEATED every 5-15 min while
  ## the reaction continues, not one shot. Modelling a single dose makes
  ## early adrenaline look useless, because with a 2.5-min half-life it
  ## has washed out before the leak peaks — which is a true statement
  ## about a single dose and a false statement about the guideline.
  p <- modifyList(par, list(SUPINE = supine))
  e <- ev_allergen(dose_mg)
  if (!is.na(epi_delay_min)) {
    e <- e + ev_epi_im(0.3, time = epi_delay_min/60, n = n_epi,
                       ii = epi_repeat_min/60)
  }
  mod %>% param(p) %>% ev(e) %>%
    mrgsim(end = 8, delta = 0.005) %>% as.data.frame()
}

##' Minutes spent below a mean arterial pressure threshold — a more
##' clinically meaningful integral than the nadir, and the quantity the
##' timing argument is actually about.
hypotension_min <- function(out, cut = 65) {
  d <- diff(out$time)
  sum(d * head(out$MAP < cut, -1)) * 60
}

##' SC3 — omalizumab, 20 weeks (OUtMATCH stage 1 geometry)
sc_omalizumab <- function(par = patient(), weeks = 20) {
  dt <- omalizumab_dose(par$WT, par$TIGE0)
  n  <- ceiling(weeks*WK/dt$interval_h)
  mod %>% param(par) %>%
    ev(ev_omalizumab(dt$dose_mg, dt$interval_h, n)) %>%
    mrgsim(end = weeks*WK, delta = 12) %>% as.data.frame()
}

##' SC4 — peanut OIT: 6 months up-dosing to 300 mg/d, 12 months
##' maintenance, then WITHDRAWAL (the PALISADE / ARC008 geometry)
sc_oit <- function(par = patient(), maint_mg = 300,
                   build_wk = 26, maint_wk = 52, off_wk = 26) {
  p <- modifyList(par, list(OITON = 1, OITDOSE = maint_mg))
  tot_d <- (build_wk + maint_wk)*7
  e <- ev_oit_daily(maint_mg, tot_d)
  ## withdrawal handled by switching OITON off through an idata split
  out1 <- mod %>% param(p) %>% ev(e) %>%
    mrgsim(end = (build_wk + maint_wk)*WK, delta = 24) %>% as.data.frame()
  st <- as.list(out1[nrow(out1), names(mrgsolve::init(mod))])
  p2 <- modifyList(par, list(OITON = 0, OITDOSE = 0, INITMODE = 0))
  out2 <- mod %>% param(p2) %>% mrgsolve::init(st) %>%
    mrgsim(end = off_wk*WK, delta = 24) %>% as.data.frame()
  out2$time <- out2$time + (build_wk + maint_wk)*WK
  rbind(out1, out2[-1, ])
}

##' SC5 — omalizumab lead-in then combined omalizumab + OIT (S3)
sc_combo <- function(par = patient(), lead_wk = 8, total_wk = 26,
                     maint_mg = 300) {
  dt <- omalizumab_dose(par$WT, par$TIGE0)
  n  <- ceiling(total_wk*WK/dt$interval_h)
  p  <- modifyList(par, list(OITON = 1, OITDOSE = maint_mg))
  e  <- ev_omalizumab(dt$dose_mg, dt$interval_h, n) +
        ev_oit_daily(maint_mg, (total_wk - lead_wk)*7, start = lead_wk*WK)
  mod %>% param(p) %>% ev(e) %>%
    mrgsim(end = total_wk*WK, delta = 24) %>% as.data.frame()
}

##' SC6 — dupilumab monotherapy, 24 weeks (the informative negative)
sc_dupilumab <- function(par = patient(), weeks = 24) {
  n <- ceiling(weeks*WK/(2*WK))
  mod %>% param(par) %>%
    ev(ev_dupilumab(300, 2*WK, n)) %>%
    mrgsim(end = weeks*WK, delta = 12) %>% as.data.frame()
}

##' SC7 — the same dose on a quiet day and after exercise + an NSAID
sc_cofactor <- function(par = patient(), dose_mg = 100) {
  quiet <- mod %>% param(par) %>% ev(ev_allergen(dose_mg)) %>%
    mrgsim(end = 8, delta = 0.01) %>% as.data.frame()
  hot <- mod %>% param(modifyList(par, list(COF_EX = 1, COF_NSAID = 1))) %>%
    ev(ev_allergen(dose_mg)) %>%
    mrgsim(end = 8, delta = 0.01) %>% as.data.frame()
  quiet$arm <- "quiet day"; hot$arm <- "exercise + NSAID"
  rbind(quiet, hot)
}

##' SC8 — adrenaline timing sweep on an identical exposure
sc_epi_timing <- function(par = patient(), dose_mg = 2000,
                          delays_min = c(2, 5, 10, 20, 30, 45, NA)) {
  do.call(rbind, lapply(delays_min, function(d) {
    o <- sc_accidental(par, dose_mg, epi_delay_min = d)
    o$delay <- ifelse(is.na(d), "none", paste0(d, " min"))
    o
  }))
}

##' SC9 — H1 antihistamine pre-treatment: the score falls, the danger
##' does not (the map's cluster-10 warning, made quantitative)
sc_antihistamine <- function(par = patient(), dose_mg = 1000) {
  plain <- mod %>% param(par) %>% ev(ev_allergen(dose_mg)) %>%
    mrgsim(end = 8, delta = 0.01) %>% as.data.frame()
  cet <- mod %>% param(par) %>%
    ev(ev_cetirizine(10, time = -2) + ev_allergen(dose_mg)) %>%
    mrgsim(start = -2, end = 8, delta = 0.01) %>% as.data.frame()
  plain$arm <- "no premedication"; cet$arm <- "cetirizine 10 mg"
  rbind(plain, cet)
}

##' SC10 — the high-total-IgE anti-IgE non-responder.
##' sIgE is raised in proportion so that f, and therefore the untreated
##' threshold, is comparable to the reference patient. The ONLY thing
##' different about this patient is how much IgE the drug has to mop up.
sc_high_ige <- function(weeks = 20, sige = 250, tige = 2000) {
  p <- patient(WT = 70, AGEY = 30, SIGE0 = sige, TIGE0 = tige)
  sc_omalizumab(p, weeks)
}

##' SC11 — natural history: milk (resolving) against peanut (persistent)
sc_natural <- function(years = 6) {
  milk   <- mod %>% param(patient(NATRES = 0.9)) %>%
    mrgsim(end = years*YR, delta = 168) %>% as.data.frame()
  peanut <- mod %>% param(patient(NATRES = 0.0)) %>%
    mrgsim(end = years*YR, delta = 168) %>% as.data.frame()
  milk$arm <- "milk (resolving)"; peanut$arm <- "peanut (persistent)"
  rbind(milk, peanut)
}

##' SC12 — the biphasic reaction
sc_biphasic <- function(par = patient(), dose_mg = 3000,
                        epi_delay_min = 25) {
  sc_accidental(par, dose_mg, epi_delay_min = epi_delay_min) %>%
    filter(time <= 8)
}

## =====================================================================
##  VIRTUAL POPULATION
##  Trials report BINARY endpoints (tolerates 600 mg: yes/no). A model
##  that only produces a median cannot be compared with them. This
##  builds a population, applies the trial entry criterion, and reports
##  the response rate the way the trial does.
## =====================================================================
vpop <- function(n = 300, seed = 20260728,
                 sige_gm = 26, sige_gsd = 2.3,
                 tige_gm = 340, tige_gsd = 2.1,
                 rel_gsd = 1.25) {
  set.seed(seed)
  sige <- rlnorm(n, log(sige_gm), log(sige_gsd))
  tige <- pmax(rlnorm(n, log(tige_gm), log(tige_gsd)), sige*1.3)
  data.frame(ID = seq_len(n), SIGE0 = sige, TIGE0 = tige,
             RELEASE = rlnorm(n, 0, log(rel_gsd)),
             WT = 35, AGEY = 10)
}

##' Compute the ED for every subject in a virtual population.
##'
##' `jitter` is the WITHIN-SUBJECT log10 standard deviation of the
##' threshold (S5) and it is applied as a multiplicative factor on
##' allergen DELIVERY, which is the mechanism cofactors actually use.
##' Every call draws independently, so entry and exit challenges are two
##' separate draws from the same patient — which is what they are in a
##' real trial, and which is the model'"'"'s explanation of the placebo
##' response rate.
vpop_ED <- function(vp, mod_par = list(), init_fun = NULL,
                    jitter = 0, seed = NULL, verbose = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  fint0 <- unname(param(mod)$FINTACT)
  out <- numeric(nrow(vp))
  for (i in seq_len(nrow(vp))) {
    jf <- if (jitter > 0) 10^rnorm(1, 0, jitter) else 1
    p <- modifyList(patient(SIGE0 = vp$SIGE0[i], TIGE0 = vp$TIGE0[i],
                            RELEASE = vp$RELEASE[i], WT = vp$WT[i],
                            AGEY = vp$AGEY[i], FINTACT = fint0*jf), mod_par)
    ini <- if (is.null(init_fun)) NULL else init_fun(p)
    out[i] <- find_ED(p, init0 = ini, tol = 0.04)
    if (verbose && i %% 25 == 0) message("  ", i, "/", nrow(vp))
  }
  out
}

## =====================================================================
##  DIAGNOSTICS
##  run_diagnostics() reproduces every number quoted in README.md.
##  Each block prints MODEL and OBSERVED side by side. Where the model
##  misses, it is printed as a MISS, not smoothed over.
## =====================================================================
fmt <- function(x, d = 2) formatC(x, format = "f", digits = d, big.mark = ",")

run_diagnostics <- function(quick = FALSE) {

  cat("\n=====================================================================\n")
  cat(" IgE-mediated food allergy / anaphylaxis — QSP model diagnostics\n")
  cat("=====================================================================\n")

  p0 <- patient()

  ## ---- D1 : the calibration patient's threshold -----------------------
  cat("\n[D1] Baseline eliciting dose of the calibration patient\n")
  ed0 <- find_ED(p0)
  cat("     model ED (grade>=2)        : ", fmt(ed0, 1), " mg peanut protein\n")
  cat("     observed (typical challenge-\n")
  cat("       proven peanut allergy)   :  ~10-100 mg, median ~30 mg\n")

  ## ---- D2 : the RATIO, not the titre (S1) -----------------------------
  cat("\n[D2] Two patients with IDENTICAL sIgE and different TOTAL IgE\n")
  edA <- find_ED(patient(SIGE0 = 40, TIGE0 = 150))
  edB <- find_ED(patient(SIGE0 = 40, TIGE0 = 300))
  edC <- find_ED(patient(SIGE0 = 40, TIGE0 = 2000))
  cat("     sIgE 40, total  150 IU/mL  (f=0.34) : ED ", fmt(edA, 1), " mg\n")
  cat("     sIgE 40, total  300 IU/mL  (f=0.17) : ED ", fmt(edB, 1), " mg\n")
  cat("     sIgE 40, total 2000 IU/mL  (f=0.025): ED ", fmt(edC, 1), " mg\n")
  cat("     PREDICTION: the same sIgE is LESS dangerous on a high total-IgE\n")
  cat("     background. This is why sIgE/total-IgE outperforms sIgE alone\n")
  cat("     (Sindher 2018; Gupta 2018). It is a consequence of f, not a fit.\n")

  ## ---- D3 : the square law, verified numerically ----------------------
  cat("\n[D3] Is the ED really quadratic in surface IgE? (S1)\n")
  base_surf <- function(p) {
    o <- mod %>% param(p) %>% mrgsim(end = 1, delta = 1) %>% as.data.frame()
    o$SURFrel[1]
  }
  ps <- list(patient(SIGE0 = 10), patient(SIGE0 = 20),
             patient(SIGE0 = 40), patient(SIGE0 = 80))
  sv <- sapply(ps, base_surf); ev_ <- sapply(ps, find_ED)
  sl <- coef(lm(log10(ev_) ~ log10(sv)))[2]
  sl_loc <- diff(log10(ev_[2:3]))/diff(log10(sv[2:3]))
  cat("     surface index (rel) : ", paste(fmt(sv, 3), collapse = "  "), "\n")
  cat("     ED (mg)             : ", paste(fmt(ev_, 1), collapse = "  "), "\n")
  cat("     slope over the whole sweep  : ", fmt(sl, 2),
      "   (theory: -2.00)\n")
  cat("     slope near the reference    : ", fmt(sl_loc, 2), "\n")
  cat("     The full-sweep slope is steeper than -2 because at the high end\n")
  cat("     the ED is large enough for the allergen term [A]/(K+[A]) to\n")
  cat("     saturate, so more than proportionally more dose is needed. The\n")
  cat("     square law is exact only while allergen is below its Kd, which\n")
  cat("     is the whole clinical range.\n")

  ## ---- D4 : the fast arm and the slow arm of anti-IgE (S2) ------------
  cat("\n[D4] Omalizumab: why the challenge is at week 16 and not week 4\n")
  o <- sc_omalizumab(p0, weeks = 20)
  pick <- function(wk) o[which.min(abs(o$time - wk*WK)), ]
  for (wk in c(0, 1, 2, 4, 8, 12, 16, 20)) {
    r <- pick(wk)
    cat(sprintf("     wk %2d  freeIgE %6.1f ng/mL   totalIgE %6.0f IU/mL   rho %.3f   SURFrel %.4f\n",
                wk, r$FreeIgE_ng, r$TotIgE_IU, r$RHO, r$SURFrel))
  }
  st4  <- state_at(p0, ev_omalizumab(omalizumab_dose(p0$WT,p0$TIGE0)$dose_mg,
                                     omalizumab_dose(p0$WT,p0$TIGE0)$interval_h,
                                     20), 4*WK)
  st16 <- state_at(p0, ev_omalizumab(omalizumab_dose(p0$WT,p0$TIGE0)$dose_mg,
                                     omalizumab_dose(p0$WT,p0$TIGE0)$interval_h,
                                     20), 16*WK)
  ed4  <- find_ED(p0, init0 = st4)
  ed16 <- find_ED(p0, init0 = st16)
  cat("     ED at week  4 : ", fmt(ed4, 0),  " mg  (", fmt(ed4/ed0, 1), "x baseline)\n")
  cat("     ED at week 16 : ", fmt(ed16, 0), " mg  (", fmt(ed16/ed0, 1), "x baseline)\n")
  cat("     observed (OUtMATCH, wk16-20): median tolerated dose 10 mg -> >1000 mg\n")
  cat("     The week-4/week-16 gap is the receptor arm. It is not a\n")
  cat("     pharmacokinetic delay: free IgE is already at target by day 7.\n")

  ## ---- D5 : total IgE rises while free IgE falls ----------------------
  cat("\n[D5] The total-IgE artefact\n")
  r20 <- pick(20)
  cat("     total IgE  baseline -> wk20 : ", fmt(o$TotIgE_IU[1], 0), " -> ",
      fmt(r20$TotIgE_IU, 0), " IU/mL  (", fmt(r20$TotIgE_IU/o$TotIgE_IU[1], 1), "x)\n")
  cat("     free  IgE  baseline -> wk20 : ", fmt(o$FreeIgE_ng[1], 0), " -> ",
      fmt(r20$FreeIgE_ng, 1), " ng/mL\n")
  cat("     observed: total IgE rises 5-10x, free IgE target <25 ng/mL\n")
  cat("     Consequence: a total-IgE assay drawn on therapy is meaningless\n")
  cat("     for re-dosing, which is why the label forbids it.\n")

  ## ---- D6 : the informative negative — dupilumab ----------------------
  cat("\n[D6] Dupilumab monotherapy: the model's explanation of a failure\n")
  d <- sc_dupilumab(p0, 24)
  dl <- d[nrow(d), ]
  cat("     sIgE      : ", fmt(d$sIgE_IU[1], 1), " -> ", fmt(dl$sIgE_IU, 1), " IU/mL\n")
  cat("     total IgE : ", fmt(d$TotIgE_IU[1], 0), " -> ", fmt(dl$TotIgE_IU, 0), " IU/mL\n")
  cat("     f         : ", fmt(d$Fspec[1], 3), " -> ", fmt(dl$Fspec, 3),
      "   <-- BARELY MOVES\n")
  stD <- state_at(p0, ev_dupilumab(300, 2*WK, 12), 24*WK)
  edD <- find_ED(p0, init0 = stD)
  cat("     ED        : ", fmt(ed0, 1), " -> ", fmt(edD, 1), " mg  (",
      fmt(edD/ed0, 2), "x)\n")
  cat("     observed  : threshold moved in ~2/24 patients (Sindher 2024-type\n")
  cat("                 dupilumab monotherapy data) — essentially no effect\n")
  cat("     WHY: dupilumab suppresses specific AND non-specific IgE together,\n")
  cat("     so f is preserved, and f is the term that is squared. It removes\n")
  cat("     the numerator and the denominator at the same time.\n")

  ## ---- D7 : the series architecture (S3) ------------------------------
  cat("\n[D7] Do immunotherapy and anti-IgE multiply?\n")
  dt <- omalizumab_dose(p0$WT, p0$TIGE0)
  st_oit <- state_at(modifyList(p0, list(OITON = 1, OITDOSE = 300)),
                     ev_oit_daily(300, 26*7), 26*WK)
  ed_oit <- find_ED(p0, init0 = st_oit)
  st_oma <- state_at(p0, ev_omalizumab(dt$dose_mg, dt$interval_h, 14), 26*WK)
  ed_oma <- find_ED(p0, init0 = st_oma)
  st_cmb <- state_at(modifyList(p0, list(OITON = 1, OITDOSE = 300)),
                     ev_omalizumab(dt$dose_mg, dt$interval_h, 14) +
                     ev_oit_daily(300, 18*7, start = 8*WK), 26*WK)
  ed_cmb <- find_ED(p0, init0 = st_cmb)
  cat("     ED baseline        : ", fmt(ed0, 1), " mg\n")
  cat("     ED OIT alone       : ", fmt(ed_oit, 0), " mg  (", fmt(ed_oit/ed0, 1), "x)\n")
  cat("     ED omalizumab alone: ", fmt(ed_oma, 0), " mg  (", fmt(ed_oma/ed0, 1), "x)\n")
  cat("     ED combination     : ", fmt(ed_cmb, 0), " mg  (", fmt(ed_cmb/ed0, 1), "x)\n")
  cat("     product of the two singles : ",
      fmt((ed_oit/ed0)*(ed_oma/ed0), 1), "x\n")
  cat("     If the combination number tracks the PRODUCT rather than the\n")
  cat("     SUM, the series architecture (S3) holds inside the model.\n")

  ## ---- D8 : the adrenaline timing integral ----------------------------
  cat("\n[D8] Adrenaline timing — an integral, not a switch\n")
  cat("     protocol: 0.3 mg IM at the stated delay, repeated at +8 and +16 min\n")
  for (d_ in c(2, 5, 10, 20, 30, 45)) {
    s <- sc_accidental(p0, 2000, epi_delay_min = d_)
    cat(sprintf("     delay %2d min : MAP nadir %5.1f   min below 65 mmHg %5.1f   peak PV deficit %4.1f%%   harm integral %5.3f L*h\n",
                d_, min(s$MAP), hypotension_min(s), 100*max(s$PVdef), max(s$AUCLEAK)))
  }
  s_none <- sc_accidental(p0, 2000, epi_delay_min = NA)
  cat(sprintf("     no adrenaline: MAP nadir %5.1f   min below 65 mmHg %5.1f   peak PV deficit %4.1f%%   harm integral %5.3f L*h\n",
              min(s_none$MAP), hypotension_min(s_none), 100*max(s_none$PVdef),
              max(s_none$AUCLEAK)))
  cat("     A single 0.3 mg dose is a ~15-minute window (t1/2 2.5 min), so\n")
  cat("     the efficacy of the DRUG does not change with delay. What grows\n")
  cat("     with delay is the deficit it has to repair, and the time spent\n")
  cat("     below the pressure at which organs are perfused. Run the same\n")
  cat("     sweep with n_epi = 1 and early adrenaline looks WORSE than late,\n")
  cat("     which is a true statement about one shot and a false statement\n")
  cat("     about the guideline.\n")

  ## ---- D9 : cofactors are threshold multipliers (S5) ------------------
  cat("\n[D9] Cofactors move the threshold, not the immunology\n")
  combos <- list("none" = list(),
                 "exercise" = list(COF_EX = 1),
                 "NSAID" = list(COF_NSAID = 1),
                 "alcohol" = list(COF_ETOH = 1),
                 "exercise+NSAID" = list(COF_EX = 1, COF_NSAID = 1),
                 "infection" = list(COF_INF = 1),
                 "PPI" = list(COF_PPI = 1))
  for (nm in names(combos)) {
    e_ <- find_ED(modifyList(p0, combos[[nm]]))
    cat(sprintf("     %-16s ED %8s mg   (%.2f log10 shift)\n",
                nm, fmt(e_, 1), log10(e_/ed0)))
  }
  cat("     observed: within-patient threshold varies ~0.4-0.5 log10 across\n")
  cat("     repeat challenges; exercise/NSAID are the strongest documented\n")
  cat("     cofactors and can lower the threshold by roughly a log.\n")

  ## ---- D10 : antihistamine masking ------------------------------------
  cat("\n[D10] H1 antihistamine: the score falls, the danger does not\n")
  for (dz in c(round(ed0*2), 1000)) {
    a <- sc_antihistamine(p0, dz)
    cat(sprintf("     dose %5d mg\n", dz))
    for (arm in unique(a$arm)) {
      s <- a[a$arm == arm, ]
      cat(sprintf("       %-18s peak skin %.2f   MAP nadir %5.1f   FEV1 nadir %5.1f   max grade %d\n",
                  arm, max(s$URTskin), min(s$MAP), min(s$FEV1), as.integer(max(s$SEV))))
    }
  }
  cat("     Near the threshold the premedicated patient is downgraded while\n")
  cat("     the physiology is unchanged: the warning sign is removed and the\n")
  cat("     lesion is left. That is the argument against premedicating a\n")
  cat("     challenge, and it is a model OUTPUT, not an assumption.\n")

  ## ---- D11 : the MAP cliff (S4) ---------------------------------------
  cat("\n[D11] Is severity graded, or is it a cliff?\n")
  for (dz in c(30, 100, 300, 1000, 3000, 10000)) {
    s <- mod %>% param(p0) %>% ev(ev_allergen(dz)) %>%
      mrgsim(end = 6, delta = 0.01) %>% as.data.frame()
    cat(sprintf("     dose %6d mg : PV deficit %4.1f%%   reserve spent %.3f   tone %.2f   MAP %5.1f   FEV1 %5.1f   tryptase %5.1f   grade %d\n",
                dz, 100*max(s$PVdef), max(s$Reserve), min(s$MAPtone),
                min(s$MAP), min(s$FEV1), max(s$Tryptase), as.integer(max(s$SEV))))
  }
  cat("     Two terms multiply into MAP: vascular TONE, which is graded, and\n")
  cat("     the RESERVE, which is a cliff. Watch the reserve column: it stays\n")
  cat("     near zero while the plasma deficit climbs, then moves suddenly.\n")
  cat("     That is S4, and it is why grading scales look ordinal.\n")

  ## ---- D12 : tryptase as a marker -------------------------------------
  cat("\n[D12] Tryptase: peak, timing, and the consensus rise criterion\n")
  s <- sc_accidental(p0, 2000, epi_delay_min = 10)
  tp <- s$time[which.max(s$Tryptase)]
  base_t <- unname(param(mod)$TRYP0)
  crit <- 1.2*base_t + 2
  cat("     peak tryptase        : ", fmt(max(s$Tryptase), 1), " ng/mL at ",
      fmt(tp*60, 0), " min\n")
  cat("     consensus rise cut-off: ", fmt(crit, 1), " ng/mL (1.2 x baseline + 2)\n")
  cat("     model exceeds cut-off : ", max(s$Tryptase) > crit, "\n")
  cat("     observed: peak 15-120 min, t1/2 ~2 h.\n")
  cat("     KNOWN MISS: in the clinic a large minority of food-triggered\n")
  cat("     anaphylaxis episodes do NOT breach this cut-off, whereas the\n")
  cat("     model breaches it for essentially any objective reaction. The\n")
  cat("     model has one well-mixed mast cell compartment; the real reason\n")
  cat("     peripheral tryptase under-reads food anaphylaxis is that the\n")
  cat("     responding cells are in gut and skin and their tryptase is\n")
  cat("     diluted before it reaches an antecubital vein. Fixing this\n")
  cat("     needs spatial structure the model does not have.\n")
  for (dz in c(10, 30, 100, 300, 1000)) {
    ss <- mod %>% param(p0) %>% ev(ev_allergen(dz)) %>%
      mrgsim(end = 6, delta = 0.02) %>% as.data.frame()
    cat(sprintf("        dose %5d mg -> peak tryptase %5.1f ng/mL  (breaches: %s)\n",
                dz, max(ss$Tryptase), max(ss$Tryptase) > crit))
  }

  ## ---- D13 : OIT time course and withdrawal ---------------------------
  cat("\n[D13] Oral immunotherapy: on over months, off over weeks\n")
  oo <- sc_oit(p0)
  for (wk in c(0, 8, 26, 52, 78, 84, 91, 104)) {
    r <- oo[which.min(abs(oo$time - wk*WK)), ]
    cat(sprintf("     wk %3d  IgG4 %6.1f ug/mL   anergy %.3f   Treg %.2f   sIgE %5.1f IU/mL\n",
                wk, r$IgG4_ug, r$ANERG, r$TREG, r$sIgE_IU))
  }
  cat("     (dosing stops at week 78)\n")
  cat("     observed: sIgG4 rises 10-100x over 6-12 months; desensitisation\n")
  cat("     is lost within weeks-to-months of stopping (ARC008, PALISADE\n")
  cat("     follow-on). The two decay rates here are the IgG4 half-life\n")
  cat("     (21 d) and the anergy half-life (3 d).\n")

  ## ---- D14 : the age gate on remission --------------------------------
  cat("\n[D14] Why age matters for OIT remission (IMPACT)\n")
  for (ag in c(1, 2, 3, 4, 6, 10, 20)) {
    pa <- patient(AGEY = ag)
    st <- state_at(modifyList(pa, list(OITON = 1, OITDOSE = 300)),
                   ev_oit_daily(300, 78*7), 78*WK)
    cat(sprintf("     age %2d y : Treg %.2f   IgG4 %5.1f ug/mL\n",
                ag, st$TREG, st$IGG4*146/1000))
  }
  cat("     observed: IMPACT remission 20%% under age 4 vs 7%% overall.\n")
  cat("     The model makes the SLOW arm age-gated and the FAST arm not,\n")
  cat("     which predicts equal desensitisation but unequal remission.\n")

  ## ---- D16 : the anti-IgE non-responder ------------------------------
  cat("\n[D16] Who does anti-IgE fail, and why\n")
  cat("     Both patients below have the SAME specific fraction f = 0.125 and\n")
  cat("     therefore comparable untreated thresholds. They differ only in\n")
  cat("     how much total IgE the drug has to neutralise.\n")
  for (cfg in list(list(w = 70, sg =  37.5, tg =  300),
                   list(w = 70, sg = 125,   tg = 1000),
                   list(w = 70, sg = 250,   tg = 2000),
                   list(w = 70, sg = 500,   tg = 4000))) {
    pp <- patient(WT = cfg$w, AGEY = 30, SIGE0 = cfg$sg, TIGE0 = cfg$tg)
    dtp <- omalizumab_dose(cfg$w, cfg$tg)
    e0  <- find_ED(pp, tol = 0.04)
    stp <- state_at(pp, ev_omalizumab(dtp$dose_mg, dtp$interval_h, 20), 16*WK)
    e1  <- find_ED(pp, init0 = stp, tol = 0.04)
    o   <- sc_omalizumab(pp, 16); l <- o[nrow(o), ]
    cat(sprintf("     total IgE %5.0f  dose %3d mg q%dwk  free IgE %5.1f ng/mL %-9s ED %8s -> %9s  (%6.1fx)\n",
                cfg$tg, dtp$dose_mg, round(dtp$interval_h/WK), l$FreeIgE_ng,
                ifelse(l$FreeIgE_ng < 25, "(at target)", "(MISSED)"),
                fmt(e0, 0), fmt(e1, 0), e1/e0))
  }
  cat("     The licensed envelope caps at 600 mg q2wk. Above it the drug can\n")
  cat("     no longer supply binding sites in stoichiometric excess, free IgE\n")
  cat("     never reaches target, the receptor arm never engages, and the\n")
  cat("     benefit collapses. This is a structural non-response, not a\n")
  cat("     pharmacodynamic one, and it is the model's single most testable\n")
  cat("     prediction about who should NOT be given anti-IgE alone.\n")

  if (quick) { cat("\n(quick mode: virtual population skipped)\n"); return(invisible(NULL)) }

  ## ---- D15 : the trial endpoint, reproduced as a rate -----------------
  cat("\n[D15] Virtual population: reproducing a BINARY trial endpoint\n")
  cat("     (this is the slow diagnostic; several minutes)\n")
  JIT <- 0.45                       # within-subject log10 SD of the threshold
  vp <- vpop(n = 160)
  ed_entry <- vpop_ED(vp, jitter = JIT, seed = 11)
  keep <- which(!is.na(ed_entry) & ed_entry <= 100)      # the entry criterion
  cat("     population n = ", nrow(vp), "; meeting entry ED<=100 mg: ",
      length(keep), "\n")
  cat("     entry ED  median ", fmt(median(ed_entry[keep]), 1), " mg   IQR ",
      fmt(quantile(ed_entry[keep], .25), 1), "-",
      fmt(quantile(ed_entry[keep], .75), 1), " mg\n")
  cat("     observed entry (OUtMATCH): reacted at <=100 mg; median maximum\n")
  cat("       TOLERATED dose 10 mg, i.e. eliciting dose ~30 mg\n")

  ## PLACEBO ARM: a second challenge in the same patients, nothing else
  ## changed. Any responder here is pure within-subject variance (S5).
  ed_pbo <- vpop_ED(vp[keep, ], jitter = JIT, seed = 22)
  cat("     placebo arm, second challenge, tolerating 600 mg : ",
      fmt(100*mean(ed_pbo >= 600, na.rm = TRUE), 1), "%\n")
  cat("     observed placebo arm                             :  7%\n")

  ## OMALIZUMAB ARM
  ed_oma_v <- vpop_ED(vp[keep, ], jitter = JIT, seed = 33,
    init_fun = function(p)
      state_at(p, ev_omalizumab(omalizumab_dose(p$WT, p$TIGE0)$dose_mg,
                                omalizumab_dose(p$WT, p$TIGE0)$interval_h, 20),
               16*WK))
  cat("     omalizumab wk16, tolerating 600 mg               : ",
      fmt(100*mean(ed_oma_v >= 600, na.rm = TRUE), 1), "%\n")
  cat("     observed (OUtMATCH stage 1)                      :  67%\n")
  cat("     median ED  entry ", fmt(median(ed_entry[keep]), 1),
      " mg  ->  omalizumab ", fmt(median(ed_oma_v, na.rm = TRUE), 0), " mg\n")
  cat("     observed    10 mg tolerated -> >1000 mg tolerated\n")
  cat("     NOTE: the placebo rate is NOT fitted. It is what falls out of\n")
  cat("     giving the same patient two challenges when the threshold has a\n")
  cat("     0.45 log10 within-subject SD (S5). A model with a deterministic\n")
  cat("     threshold cannot produce a placebo responder at all.\n")

  cat("\n=====================================================================\n")
  cat(" End of diagnostics.\n")
  cat("=====================================================================\n")
  invisible(NULL)
}

## =====================================================================
##  CALIBRATION NOTES  (what each number came from)
## ---------------------------------------------------------------------
##  ALLERGEN
##    FINTACT 0.0022  — the intact-protein fraction crossing the gut is
##       not directly measurable in humans; it is set so that the
##       calibration patient's ED lands at ~30 mg peanut protein, which
##       is the median of challenge-proven cohorts. It is the model's
##       single free scaling parameter and it is confounded with XL50.
##       Everything downstream is a RATIO and is therefore insensitive
##       to it, which is the point of building the model on thresholds.
##    KDEG_GUT / KA_GUT — Ara h 2 is a pepsin-resistant 2S albumin, so
##       the surviving fraction is high relative to labile allergens;
##       ARA h 8 (PR-10) is represented by lowering FARAH2.
##
##  THE SURFACE
##    KD_FCERI 0.1 nM — IgE-FcepsilonRI affinity, ~1e-10 M (Kinet 1999).
##       At normal total IgE (3-4 nM) the receptor is >95% occupied,
##       which is why occupancy only becomes a lever once free IgE is
##       driven into the sub-nM range. This is a structural reason why
##       anti-IgE has to be dosed to a TARGET, not to a percentage.
##    RHO_FLOOR 0.18 / KD_RHO 0.35 — calibrated so that omalizumab
##       produces a ~3.6-fold fall in TISSUE mast cell FcepsilonRI. Note
##       this is deliberately SMALLER than the ~10-30-fold fall measured
##       on BASOPHILS (MacGlashan 1997; Beck 2004): the basophil is the
##       convenient cell, not the effector cell, and using the basophil
##       number would over-predict efficacy by an order of magnitude.
##    KOUT_RHO — set from the 8-12 week time course of basophil
##       FcepsilonRI loss.
##
##  ANTI-IgE PK
##    CL/V/ka from the omalizumab population PK literature (Lowe 2009).
##    KEL_CPX/KEL_IGE ratio (9x) is set by the observed 5-10x rise in
##    TOTAL IgE on therapy — that ratio is the observable, not a fit.
##
##  IMMUNOTHERAPY
##    KSEC_G4 / KDPC_G4 — set so maintenance 300 mg/d gives peanut-
##       specific IgG4 of ~50-60 ug/mL at plateau (PALISADE, POISED).
##    KI_G4 40 nM — polyclonal IgG4 has lower average affinity than a
##       monoclonal; this is an effective, not a thermodynamic, constant.
##    KON_ANG / KOFF_ANG — the 3-day anergy half-life is set by the
##       clinical observation that missed doses restore reactivity within
##       days, which is the basis of every OIT missed-dose protocol.
##    AGE50_TOL 4 y — from IMPACT (Jones 2022): remission is strongly
##       age-dependent and the effect is concentrated under age 4.
##
##  ANAPHYLAXIS PHYSIOLOGY
##    35% plasma volume shift in ~10 min (Fisher 1986) sets KLEAK and,
##       with FRAC50 0.42 and HFRAC 6, produces the flat-then-cliff MAP.
##    Adrenaline IM 0.3 mg: Cmax ~0.4 ng/mL at ~8 min, t1/2 ~2.5 min
##       (Simons 1998/2001) — reproduced by KA_EPI_IM/KEL_EPI/V_EPI.
##    Tryptase t1/2 ~2 h and the 1.2 x baseline + 2 rise criterion
##       (Valent 2012) are used as a check, not as a fit.
##
##  KNOWN MISSES (stated here rather than hidden)
##    * The model over-predicts the median omalizumab threshold shift
##      relative to the OUtMATCH median; the trial's binary endpoint is
##      reproduced better than its median, because the endpoint is
##      dominated by the low tail of the entry ED distribution.
##    * The combination arm (D7) is almost certainly over-predicted: a
##      strictly multiplicative series with no shared saturation is the
##      model'"'"'s most aggressive assumption and the one most worth
##      testing against OUtMATCH stage 2 when it reports.
##    * Basophil and mast cell are collapsed into one surface. The
##      basophil activation test is therefore a proxy in the model in
##      the same way it is a proxy in the clinic — the model cannot
##      explain BAT/challenge discordance because it has no mechanism
##      for it.
##    * No spatial structure: gut, skin and airway mast cells see the
##      same allergen concentration, which is why the model gets the
##      ORDER of symptom onset only through differing rate constants.
## =====================================================================

if (identical(environment(), globalenv()) && !interactive()) {
  ## sourcing the file from Rscript runs the fast diagnostics
  run_diagnostics(quick = TRUE)
}
