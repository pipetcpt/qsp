## =============================================================================
##  dka_mrgsolve_model.R
##  Diabetic Ketoacidosis / Hyperglycaemic Hyperosmolar State — QSP model
##  42 ordinary differential equations, mrgsolve (C++ $ODE block)
## =============================================================================
##
##  WHAT IS UNUSUAL ABOUT THIS MODEL, AND WHY
##  -----------------------------------------
##  1. BICARBONATE IS NOT A STATE VARIABLE.  Plasma [HCO3-] is obtained at every
##     derivative evaluation by solving the physicochemical (Stewart) condition
##
##         [HCO3-](pH,PCO2) + [A-](pH) + [ketoanion] + [lactate] + [SIG] = SID
##
##     for pH by bisection inside $ODE, where SID = Na + K + (Ca,Mg) - Cl -
##     (unmetabolised infused organic anion).  This is more work per step than a
##     bicarbonate ODE, and it buys the following for free, with no book-keeping
##     term anywhere in the code:
##       * making a ketoacid from a neutral triglyceride adds a strong anion, so
##         HCO3- falls 1:1 with the ketoanion;
##       * OXIDISING a ketoanion removes it, so HCO3- is regenerated 1:1 —
##         retained ketoanions are literally "potential bicarbonate";
##       * excreting a ketoanion with Na+/K+ removes the anion AND a strong
##         cation: the anion gap closes but the bicarbonate does not recover, so
##         an organic acidosis is silently converted into a hyperchloraemic one;
##       * excreting it with NH4+ removes only the anion, so base is preserved,
##         and NH4Cl excretion (which also removes Cl-) is how the kidney repairs
##         a hyperchloraemic acidosis;
##       * 0.9% saline (SID 0) diluting a plasma SID of ~45 is acidifying,
##         automatically, with no "dilutional acidosis" term;
##       * acetate/gluconate/lactate in balanced crystalloids count as strong
##         anions until metabolised, then vanish: transient dip, then base gain.
##
##  2. TWO INSULIN EFFECT SITES.  Peripheral actions (lipolysis, glucose
##     disposal, the potassium shift) follow plasma insulin with tau 0.35 h
##     (GLUT4 translocation); hepatic glucose output follows with tau 1.5 h
##     (gluconeogenic wind-down).  That separation, not any fitted rate, is what
##     gives the observed shape of the treated glucose curve.
##
##  3. ENDOGENOUS INSULIN IS PRIVILEGED AT THE LIVER.  Secreted insulin passes
##     the liver first (~50% first-pass extraction), so portal concentration is
##     several-fold peripheral.  The CPT-1 gate therefore sees PORTF*20*SEC more
##     insulin than the adipocyte does.  Exogenous insulin has no such privilege.
##     This one asymmetry is why residual beta-cell function prevents ketosis
##     while permitting extreme hyperglycaemia — i.e. it is what makes HHS.
##
##  4. THE RENAL ESCAPE VALVE.  Renal glucose loss is filtered load minus a
##     splayed, glomerulotubular-balanced Tm.  Because filtration falls with ECF
##     volume, hyperglycaemia is partly a VOLUME disease: the conductance of
##     glucose's second exit is the extracellular volume.
##
##  5. THE PRESENTING STATE IS COMPUTED, NOT TYPED IN.  Every treatment scenario
##     below begins by integrating an untreated lead-in from a healthy steady
##     state with insulin withdrawn, and starts therapy from wherever that lands.
##
##  Companion file dka_reference_check.py is an independent, dependency-free
##  transcription of these same 42 equations; dka_reference_output.txt is its
##  committed output and the source of every number quoted in README.md.
##
##  EDUCATIONAL / RESEARCH USE ONLY.  Not validated for clinical decisions.
## =============================================================================

library(mrgsolve)
library(dplyr)

code <- '
$PARAM @annotated
// ---- anthropometry --------------------------------------------------------
BW      :  70.0  : body weight (kg)
FTBW    :   0.60 : total body water fraction of BW
FECF    :   0.20 : extracellular fluid fraction of BW
CAMG    :   6.0  : Ca2+ + Mg2+ contribution to SID (mEq/L at normal ECF volume)
ALB     :  42.0  : albumin (g/L at normal ECF volume)
SIGO    :   5.5  : baseline unmeasured strong anions (mEq/L at normal ECF volume)
OSMBASE : 288.33 : baseline total osmolality, ICF = ECF at rest (mOsm/kg)

// ---- insulin pharmacokinetics ---------------------------------------------
VINS    :   6.0  : insulin distribution volume (L)
CLINS   :  50.0  : insulin clearance (L/h) -- 1 U/h gives 20 uU/mL at steady state
KA_SC   :   1.50 : s.c. absorption, rapid-acting analogue (1/h)
KA_REG  :   0.50 : s.c. absorption, regular human insulin (1/h)
INSB    :  10.0  : basal plasma insulin in a normal subject (uU/mL)
SECMAX  :   8.0  : maximal endogenous secretion at full beta-cell mass (U/h)
KGSEC   :   8.0  : glucose EC50 for secretion (mmol/L, Hill 2)
PORTF   :   3.0  : portal:peripheral increment for ENDOGENOUS insulin
TINSE_F :   0.35 : peripheral insulin effect-site time constant (h)
TINSE_S :   1.50 : hepatic insulin effect-site time constant (h)

// ---- insulin pharmacodynamics: THE ASYMMETRY ------------------------------
IC50_LIP:  15.0  : half-maximal suppression of lipolysis (uU/mL)
EMAX_LIP:   0.92 : maximal fractional suppression of lipolysis
IC50_HGP:  30.0  : half-maximal suppression of hepatic glucose output (uU/mL)
EMAX_HGP:   0.85 : maximal fractional suppression of hepatic output
EC50_UP :  60.0  : half-maximal stimulation of glucose disposal (uU/mL)
IC50_KSH:  25.0  : half-maximal Na/K-ATPase (potassium shift) (uU/mL)

// ---- glucose --------------------------------------------------------------
HGP0    :  46.6  : basal hepatic glucose output (mmol/h)
GNGF    :   0.75 : gluconeogenic (glycogen-independent) fraction of output
GUCNS   :  52.2  : obligatory CNS/erythrocyte uptake Vmax (mmol/h)
KM_CNS  :   1.5  : Km of obligatory uptake (mmol/L) -- saturated at normoglycaemia
CLNI    :   0.524: insulin-independent mass-action clearance (L/h)
VMI_UP  :  70.0  : insulin-dependent disposal Vmax (mmol/h)
KM_UP   :   8.0  : Km of insulin-dependent disposal (mmol/L)
KGLYSYN :   0.25 : glycogen synthesis scaling

// ---- lipolysis and NEFA ---------------------------------------------------
LIPMAX  :  26.0  : maximal whole-body NEFA release (mmol/h)
KFFA    :  18.0  : NEFA disposal rate constant (1/h)
KIFFA   :   2.0  : saturation constant of NEFA disposal (mmol/L)
FHEP    :   0.50 : fraction of NEFA disposal taken up by the liver
VFFA    :   3.2  : NEFA distribution volume (L)
ADIPOSE :   1.0  : adipose lipolytic capacity factor

// ---- ketogenesis and ketone disposal --------------------------------------
KGSCALE :   4.0  : mol ketoacid per mol hepatic NEFA at a fully open CPT-1 gate
MAL0    :   1.0  : malonyl-CoA baseline
KMAL    :   1.0  : CPT-1 inhibition constant
EMAL    :   5.0  : malonyl-CoA response to portal insulin
IC50_MAL:  12.0  : portal insulin for half-maximal malonyl-CoA rise (uU/mL)
VMAX_KOX:  75.0  : saturable peripheral ketone oxidation Vmax (mmol/h)
KM_KOX  :   4.0  : Km of ketone oxidation (mmol/L)
KLIN_KOX:   1.20 : non-saturable (brain/heart) ketone extraction (L/h)
KOX_INS :   0.5  : insulin enhancement of ketone oxidation (Emax)
KCONV   :  10.0  : AcAc <-> BHB interconversion (1/h)
KEQ_BHB :   3.0  : BHB/AcAc ratio at normal hepatic redox
KACETONE:   0.020: spontaneous AcAc decarboxylation (1/h)
KACETCL :   0.030: acetone elimination (1/h) -- pulmonary, t1/2 ~23 h
TREDOX  :   1.0  : hepatic redox time constant (h)
ALCOHOL :   0.0  : ethanol-induced cytosolic redox load (0/1)

// ---- lactate --------------------------------------------------------------
LACP0   :  60.0  : basal lactate production (mmol/h)
KLACCL  :  60.0  : lactate clearance (L/h)

// ---- renal ----------------------------------------------------------------
GFRMAX  :   7.2  : GFR at full ECF volume (L/h = 120 mL/min)
TGFR    :   0.5  : GFR equilibration time constant (h)
VREL_CRIT:  0.55 : relative ECF volume at which GFR reaches zero
TMGLU   :  75.0  : maximal tubular glucose reabsorption at normal GFR (mmol/h)
KSPLAY  :   3.0  : glucose titration-curve splay (mmol/L)
FTM_GFR :   0.65 : fraction of Tm that tracks GFR (glomerulotubular balance)
SGLT2   :   0.0  : SGLT2 inhibitor on board (0/1) -- lowers Tm by 80%
FEKET0  :   0.05 : baseline fractional excretion of ketoanion
FEKET1  :   0.10 : load-dependent increment in ketoanion FE
KMFEK   :   8.0  : half-saturation of that increment (mmol/L)
FECL    :   0.014: baseline fractional excretion of chloride
FENA    :   0.011: baseline fractional excretion of sodium
FEK0    :   0.08 : baseline fractional excretion of potassium
FEK_ALDO:   0.30 : aldosterone increment in potassium FE
NH4MAX  :  12.0  : maximal renal ammoniagenesis (mmol/h)
NH4B    :   1.7  : baseline renal ammoniagenesis (mmol/h)
TNH4    :  12.0  : ammoniagenic adaptation time constant (h)
FEUREA  :   0.45 : fractional excretion of urea
PRODUREA:  12.0  : basal urea production (mmol/h)
PCR     :   6.48 : creatinine production (mg/dL * L/h)
UOSM_MIN: 300.0  : minimum urine osmolality (mOsm/L)
UOSM_MAX: 900.0  : maximum urine osmolality (mOsm/L)

// ---- potassium ------------------------------------------------------------
K0      :   4.2  : normal plasma potassium (mmol/L)
KIC0    : 125.0  : normal intracellular potassium (mmol/L)
GK      :  25.0  : transmembrane potassium conductance (L/h)
BINS_K  :   0.90 : set-point shift per unit fractional insulin-effect loss
BPH_K   :   1.20 : set-point shift per pH unit
BOSM_K  :   0.015: set-point shift per mOsm/kg above 285

// ---- water and volume -----------------------------------------------------
TSHIFT  :   0.25 : ICF<->ECF osmotic equilibration (h)
INSENS  :   0.040: insensible water loss (L/h)
POMAX   :   0.65 : maximal sustainable voluntary drinking (L/h)
OSM_THIRST: 292.0: osmotic threshold for thirst (mOsm/kg)
KVOMIT  :   0.85 : maximal fractional suppression of intake by ketosis
KMVOMIT :   8.0  : ketone concentration for half that suppression (mmol/L)
WATER   :   1.0  : access to / drive for oral water (1 = intact, 0 = none)

// ---- counter-regulatory hormones ------------------------------------------
GCG0    :  25.0  : basal glucagon (pmol/L)
TGCG    :   0.5  : glucagon time constant (h)
CORT0   :  12.0  : basal cortisol (ug/dL)
TCORT   :   1.5  : cortisol time constant (h)
EPI0    :   0.30 : basal epinephrine (nmol/L)
TEPI    :   0.25 : epinephrine time constant (h)
TIR     :   3.0  : insulin-resistance time constant (h)
TALDO   :   1.0  : aldosterone time constant (h)
TILL    :  30.0  : resolution of the precipitating illness (h)
ILL0    :   0.55 : severity of the precipitating illness at t = 0

// ---- brain ----------------------------------------------------------------
KOSMB   :   0.85 : idiogenic osmole per mOsm/kg of plasma excess
TOSMB_UP:   8.0  : osmolyte accumulation time constant (h)
TOSMB_DN:  14.0  : osmolyte washout time constant (h) -- deliberately slower
TVBR    :   0.30 : brain water equilibration (h)
COMPL   :   0.35 : intracranial compliance damping of swelling
KICP    : 900.0  : mmHg per unit relative brain volume
KINJ    :   0.075: ischaemic-injury accrual (1/h)
KRES    :   0.060: injury resolution (1/h)
KVASO   :   0.090: brain-volume contribution per unit injury

// ---- respiratory and CNS --------------------------------------------------
TPCO2   :   0.40 : PCO2 equilibration (h)
PCO2_FLOOR: 14.0 : minimum achievable PCO2 (mmHg)
PCO2_NORM :  41.5: normal PCO2 (mmHg)
TMENT   :   0.5  : mental-status time constant (h)
W_OSM   :   0.025: GCS weight on hyperosmolality
W_PH    :   2.2  : GCS weight on acidaemia
W_VBR   :  12.0  : GCS weight on brain swelling

// ---- disease knobs --------------------------------------------------------
BETA    :   0.0  : residual beta-cell function (0 = type 1 diabetes)
GLYCO0  : 300.0  : initial hepatic glycogen (mmol glucosyl units)

// ---- treatment inputs (piecewise-constant, set by the scenario) -----------
RATE_FL :   0.0  : intravenous fluid rate (L/h)
FL_NA   : 154.0  : fluid sodium (mmol/L)
FL_CL   : 154.0  : fluid chloride (mmol/L)
FL_K    :   0.0  : fluid potassium from the bag itself (mmol/L)
FL_ORG  :   0.0  : fluid metabolisable organic anion, acetate/gluconate (mmol/L)
FL_GLC  :   0.0  : fluid dextrose (mmol/L) -- D5W = 278
FL_LAC  :   0.0  : fluid lactate (mmol/L) -- Ringer lactate = 28
KCL     :   0.0  : added KCl (mmol per L of fluid)
KPO4    :   0.0  : potassium phosphate (mmol/h)
INS_IV  :   0.0  : intravenous insulin infusion (U/h)
BICARB  :   0.0  : 8.4% NaHCO3 (L/h; 1000 mmol/L)

$CMT @annotated
VECF    : extracellular fluid volume (L)
VICF    : intracellular fluid volume (L)
NAE     : extracellular sodium (mmol)
CLE     : extracellular chloride (mmol)
KE      : extracellular potassium (mmol)
KI      : intracellular potassium (mmol)
PHOSE   : extracellular phosphate (mmol)
ORGA    : infused metabolisable organic anion (mmol)
LAC     : lactate (mmol, total body)
ACAC    : acetoacetate (mmol, total body)
BHB     : beta-hydroxybutyrate (mmol, total body)
ACET    : acetone (mmol, total body)
GLU     : extracellular glucose (mmol)
GLYCO   : hepatic glycogen (mmol glucosyl units)
FFA     : plasma non-esterified fatty acid (mmol)
UREA    : urea (mmol, total body water)
CREA    : plasma creatinine (mg/dL)
PCO2    : arterial PCO2 (mmHg)
INSSC   : subcutaneous insulin depot (U)
INSP    : plasma insulin (uU/mL)
INSEF   : peripheral insulin effect site (uU/mL)
INSES   : hepatic insulin effect site (uU/mL)
GCG     : glucagon (pmol/L)
CORT    : cortisol (ug/dL)
EPI     : epinephrine (nmol/L)
IR      : insulin-resistance index
REDOX   : hepatic NADH/NAD+ index (1 = normal)
GFRR    : relative glomerular filtration rate
OSMB    : brain idiogenic osmoles above normal (mOsm/kg)
VBR     : brain water volume, relative
ALDO    : aldosterone activity index
NH4C    : renal ammoniagenic capacity (mmol/h)
BETAF   : beta-cell functional suppression factor
ILL     : precipitating illness severity
MENT    : mental status index (1 = alert)
INJ     : cerebral ischaemic injury index
UKET    : cumulative urinary ketoanion with Na+/K+ (mmol) -- base lost
UKETN   : cumulative urinary ketoanion with NH4+ (mmol) -- base kept
UGLU    : cumulative urinary glucose (mmol)
UKCUM   : cumulative urinary potassium (mmol)
UVOL    : cumulative urine volume (L)
CLIN    : cumulative chloride infused (mmol)

$MAIN
double TBW0 = FTBW * BW;
double VE0  = FECF * BW;
double VI0  = TBW0 - VE0;

VECF_0  = VE0;
VICF_0  = VI0;
NAE_0   = 140.0 * VE0;
CLE_0   = 105.0 * VE0;
KE_0    = K0 * VE0;
KI_0    = KIC0 * VI0;
PHOSE_0 = 1.15 * VE0;
ORGA_0  = 0.0;
LAC_0   = 1.0 * (VE0 + 0.5 * VI0);
ACAC_0  = 0.025 * (VE0 + 0.6 * VI0);
BHB_0   = 0.075 * (VE0 + 0.6 * VI0);
ACET_0  = 0.02 * TBW0;
GLU_0   = 5.0 * VE0;
GLYCO_0 = GLYCO0;
FFA_0   = 0.35 * VFFA;
UREA_0  = 3.33 * TBW0;
CREA_0  = 0.90;
PCO2_0  = PCO2_NORM;
INSSC_0 = 0.0;
INSP_0  = INSB;
INSEF_0 = INSB;
INSES_0 = INSB;
GCG_0   = GCG0;
CORT_0  = CORT0;
EPI_0   = EPI0;
IR_0    = 1.0;
REDOX_0 = 1.0;
GFRR_0  = 1.0;
OSMB_0  = 0.0;
VBR_0   = 1.0;
ALDO_0  = 1.0;
NH4C_0  = NH4B;
BETAF_0 = 1.0;
ILL_0   = ILL0;
MENT_0  = 1.0;
INJ_0   = 0.0;

$ODE
// ===========================================================================
//  0.  DERIVED CONCENTRATIONS
// ===========================================================================
double VE   = (VECF  > 4.0 ? VECF  : 4.0);
double VI   = (VICF  > 8.0 ? VICF  : 8.0);
double VE0d = FECF * BW;
double VI0d = (FTBW - FECF) * BW;
double TBW  = VE + VI;
double VKET = VE + 0.6 * VI;
double VLAC = VE + 0.5 * VI;
double vrel = VE / VE0d;
double conc = VE0d / VE;                 // ECF contraction factor

double cNa   = NAE / VE;
double cCl   = CLE / VE;
double cK    = KE  / VE;
double cPhos = PHOSE / VE;
double cOrg  = ORGA / VE;
double cAcAc = ACAC / VKET;
double cBHB  = BHB  / VKET;
double cKet  = cAcAc + cBHB;
double cLac  = LAC / VLAC;
double Gp    = GLU / VE;
double cUrea = UREA / TBW;

// ===========================================================================
//  1.  ACID-BASE:  physicochemical (Stewart) solve for pH, then [HCO3-]
//      Albumin and the unmeasured strong anions CONCENTRATE with the ECF;
//      holding them fixed while Na and Cl concentrate manufactures an alkalosis.
// ===========================================================================
double SIDapp = cNa + cK + CAMG * conc - cCl - cOrg;
double albc   = ALB  * conc;
double sigoc  = SIGO * conc;
double lo = 5.60, hi = 8.30, pHm, hres;
for (int it = 0; it < 45; ++it) {
  pHm  = 0.5 * (lo + hi);
  hres = 0.0301 * PCO2 * pow(10.0, pHm - 6.1)
       + albc  * (0.123 * pHm - 0.631)
       + cPhos * (0.309 * pHm - 0.469)
       + cKet + cLac + sigoc - SIDapp;
  if (hres > 0.0) hi = pHm; else lo = pHm;
}
double pH   = 0.5 * (lo + hi);
double HCO3 = 0.0301 * PCO2 * pow(10.0, pH - 6.1);
double AG   = cNa - cCl - HCO3;
double OSM_EFF = 2.0 * cNa + Gp;
double OSM_TOT = OSM_EFF + cUrea;

// ===========================================================================
//  2.  INSULIN:  two effect sites, and a portal privilege for endogenous insulin
// ===========================================================================
double IRs   = (IR > 0.3 ? IR : 0.3);
double Ieff  = INSEF / IRs;              // peripheral: minutes
double IeffS = INSES / IRs;              // hepatic: hours
double fLIP  = 1.0 - EMAX_LIP * Ieff  / (Ieff  + IC50_LIP);
double fHGP  = 1.0 - EMAX_HGP * IeffS / (IeffS + IC50_HGP);
double fUP   = Ieff / (Ieff + EC50_UP);
double fKSH  = Ieff / (Ieff + IC50_KSH);
double fKSH0 = INSB / (INSB + IC50_KSH);

double sec = SECMAX * BETA * BETAF * Gp * Gp / (Gp * Gp + KGSEC * KGSEC);
double Riv = INS_IV + sec;                                       // U/h
double abs_sc = KA_SC * INSSC;
dxdt_INSSC = -KA_SC * INSSC;
dxdt_INSP  = ((Riv + abs_sc) * 1.0e6 - CLINS * 1000.0 * INSP) / (VINS * 1000.0);
dxdt_INSEF = (INSP - INSEF) / TINSE_F;
dxdt_INSES = (INSP - INSES) / TINSE_S;

double acid_sev = (7.30 - pH) / 0.30;
if (acid_sev < 0.0) acid_sev = 0.0;
if (acid_sev > 1.5) acid_sev = 1.5;
double gtox = (Gp - 12.0) / 20.0;
if (gtox < 0.0) gtox = 0.0; if (gtox > 1.0) gtox = 1.0;
double asup = (7.30 - pH) / 0.30;
if (asup < 0.0) asup = 0.0; if (asup > 1.0) asup = 1.0;
double betaf_t = 1.0 - 0.55 * gtox - 0.60 * asup;
if (betaf_t < 0.05) betaf_t = 0.05;
dxdt_BETAF = (betaf_t - BETAF) / 6.0;

// ===========================================================================
//  3.  COUNTER-REGULATORY HORMONES (capped: an unbounded loop is not physiology)
// ===========================================================================
double hypovol = (1.0 - vrel) / 0.25;
if (hypovol < 0.0) hypovol = 0.0; if (hypovol > 1.0) hypovol = 1.0;
double hypo = (4.0 - Gp) / 2.0; if (hypo < 0.0) hypo = 0.0;

double gcg_t = GCG0 * (1.0 + 2.2 * (1.0 - Ieff / (Ieff + 40.0))
                       + 0.8 * ILL + 1.5 * hypo);
dxdt_GCG = (gcg_t - GCG) / TGCG;
double cort_t = CORT0 * (1.0 + 1.8 * ILL + 0.9 * acid_sev + 0.6 * hypovol);
dxdt_CORT = (cort_t - CORT) / TCORT;
double epi_t = EPI0 * (1.0 + 3.0 * ILL + 2.5 * hypovol + 6.0 * hypo
                       + 1.2 * acid_sev);
dxdt_EPI = (epi_t - EPI) / TEPI;

double epirel = EPI / EPI0 - 1.0;  if (epirel > 2.0) epirel = 2.0;
double ffarel = FFA / VFFA / 0.35 - 1.0; if (ffarel > 2.0) ffarel = 2.0;
if (ffarel < 0.0) ffarel = 0.0;
double ir_t = 1.0 + 0.90 * ILL + 0.55 * acid_sev + 0.30 * epirel + 0.90 * ffarel;
dxdt_IR  = (ir_t - IR) / TIR;
dxdt_ILL = -ILL / TILL;
double aldo_t = 1.0 + 2.5 * hypovol;
dxdt_ALDO = (aldo_t - ALDO) / TALDO;

double cr_gluc = 1.0 + 0.90 * (GCG / GCG0 - 1.0) + 0.50 * (EPI / EPI0 - 1.0)
               + 0.30 * (CORT / CORT0 - 1.0);
if (cr_gluc < 0.4) cr_gluc = 0.4; if (cr_gluc > 3.0) cr_gluc = 3.0;
double cr_lip = 1.0 + 0.55 * (EPI / EPI0 - 1.0) + 0.30 * (CORT / CORT0 - 1.0);
if (cr_lip < 0.4) cr_lip = 0.4; if (cr_lip > 2.2) cr_lip = 2.2;

// ===========================================================================
//  4.  LIPOLYSIS, HEPATIC NEFA UPTAKE, KETOGENESIS
// ===========================================================================
double LIP  = LIPMAX * fLIP * cr_lip * ADIPOSE;
double cffa = FFA / VFFA;
double kffa = KFFA / (1.0 + cffa / KIFFA);
double disp_ffa = kffa * FFA;
double HFFA = FHEP * disp_ffa;
dxdt_FFA = LIP - disp_ffa;

// the CPT-1 gate sees PORTAL insulin: endogenous secretion is privileged
double IeffP = (INSEF + PORTF * 20.0 * sec) / IRs;
double gcgx  = GCG / GCG0 - 1.0; if (gcgx < 0.0) gcgx = 0.0;
double mal   = MAL0 * (1.0 + EMAL * IeffP / (IeffP + IC50_MAL))
             / (1.0 + 1.2 * gcgx);
double cpt1  = 1.0 / (1.0 + mal / KMAL);
double KGEN  = KGSCALE * HFFA * cpt1;                 // mmol ketoacid / h

double hfx = HFFA / 8.2 - 1.0; if (hfx < 0.0) hfx = 0.0;
double redox_t = 1.0 + 0.75 * hfx + 3.0 * ALCOHOL;
dxdt_REDOX = (redox_t - REDOX) / TREDOX;

double KOX = VMAX_KOX * (1.0 + KOX_INS * fUP) * cKet / (cKet + KM_KOX)
           + KLIN_KOX * cKet;
double facac = 1.0 / (1.0 + KEQ_BHB * REDOX);
double Jconv = KCONV * (KEQ_BHB * REDOX * ACAC - BHB);
double Jacet = KACETONE * ACAC;

// ===========================================================================
//  5.  RENAL HANDLING
// ===========================================================================
double GFRL = GFRR * GFRMAX;                             // L/h
double tmg  = TMGLU * (1.0 - 0.80 * SGLT2)
            * ((1.0 - FTM_GFR) + FTM_GFR * GFRR);
double UGLUr = GFRL * Gp - tmg * Gp / (Gp + KSPLAY);
if (UGLUr < 0.0) UGLUr = 0.0;

double feket = FEKET0 + FEKET1 * cKet / (cKet + KMFEK);
double UKETr = GFRL * cKet * feket;
double nh4act = 0.25 + 1.4 * acid_sev; if (nh4act > 1.0) nh4act = 1.0;
double UNH4 = NH4C * nh4act * GFRR;
double nh4_t = NH4B + (NH4MAX - NH4B) * (acid_sev > 1.0 ? 1.0 : acid_sev);
dxdt_NH4C = (nh4_t - NH4C) / TNH4;

// osmotic diuresis washes out the corticomedullary gradient: fractional
// excretion of Na, Cl and K rises steeply with non-reabsorbed solute load
double osmload = UGLUr + UKETr;
double fdiur = osmload / 60.0; if (fdiur > 1.0) fdiur = 1.0;

// chloride handling is keyed to the Cl/Na RATIO, not to absolute [Cl-]: the
// kidney defends the strong-ion difference, and it must be able to work in
// both directions or a hypernatraemic patient develops a spurious alkalosis
double clx = 1.1 * (cCl / cNa - 0.750) / 0.060;
if (clx < -1.2) clx = -1.2; if (clx > 1.5) clx = 1.5;
double UCL = GFRL * cCl * FECL * (1.0 + 3.5 * fdiur)
           / (1.0 + 0.9 * (ALDO - 1.0)) * exp(clx);
// NH4+ leaves as NH4Cl whenever chloride is available: net acid excretion,
// and in the strong-ion framework removing Cl- RAISES SID, i.e. makes base
UCL = UCL + UNH4 * (cCl / 100.0 < 1.0 ? cCl / 100.0 : 1.0);

double UNA = GFRL * cNa * FENA * (1.0 + 4.0 * fdiur) / (1.0 + 0.9 * (ALDO - 1.0));
double ukx = UKETr / 10.0; if (ukx > 1.0) ukx = 1.0;
double fek = (FEK0 + FEK_ALDO * (ALDO - 1.0) / 2.5)
           * (1.0 + 1.2 * fdiur) * (1.0 + 0.8 * ukx);
double UK  = GFRL * cK * fek;
double UOA = GFRL * cOrg * 0.35;
double UUREA = GFRL * cUrea * FEUREA * (0.55 + 0.45 * GFRR);

double sol  = UGLUr + UUREA + UNA + UK + UNH4 + UCL + UKETr + UOA;
double uosm = 320.0 + 340.0 * exp(-sol / 35.0);
if (uosm < UOSM_MIN) uosm = UOSM_MIN; if (uosm > UOSM_MAX) uosm = UOSM_MAX;
double UV = sol / uosm;                                  // L/h urine flow

// ===========================================================================
//  6.  GLUCOSE
// ===========================================================================
double glyco_av = GLYCO / 120.0; if (glyco_av > 1.0) glyco_av = 1.0;
if (glyco_av < 0.0) glyco_av = 0.0;
double HGP = HGP0 * fHGP * cr_gluc * (GNGF + (1.0 - GNGF) * glyco_av);
double gsy = KGLYSYN * fUP * (Gp > 4.0 ? Gp - 4.0 : 0.0)
           * (1.0 - GLYCO / 450.0 > 0.0 ? 1.0 - GLYCO / 450.0 : 0.0);
double UPT = GUCNS * Gp / (Gp + KM_CNS)        // obligatory, already saturated
           + CLNI * Gp                        // mass-action GLUT1
           + VMI_UP * fUP * Gp / (Gp + KM_UP); // insulin-dependent GLUT4
double ginf = RATE_FL * FL_GLC;
dxdt_GLU   = HGP + ginf - UPT - UGLUr - gsy;
dxdt_GLYCO = gsy - HGP * (1.0 - GNGF);

// ===========================================================================
//  7.  KETONES AND LACTATE
// ===========================================================================
dxdt_ACAC = KGEN * facac - Jconv - Jacet - (KOX + UKETr) * facac;
dxdt_BHB  = KGEN * (1.0 - facac) + Jconv - (KOX + UKETr) * (1.0 - facac);
dxdt_ACET = Jacet - KACETCL * ACET;
double lacp = LACP0 * (1.0 + 0.9 * (1.0 - vrel > 0.0 ? (1.0 - vrel) / 0.2 : 0.0)
                       + 0.4 * (EPI / EPI0 - 1.0));
dxdt_LAC = lacp - KLACCL * cLac * (0.5 + 0.5 * GFRR) + RATE_FL * FL_LAC;

// ===========================================================================
//  8.  ELECTROLYTES
// ===========================================================================
double na_in = RATE_FL * FL_NA + BICARB * 1000.0;
double cl_in = RATE_FL * (FL_CL + KCL);           // KCL is mmol per L of fluid
double k_in  = RATE_FL * (FL_K + KCL) + KPO4;
double org_in = RATE_FL * FL_ORG;

dxdt_NAE  = na_in - UNA;
dxdt_CLE  = cl_in - UCL;
dxdt_ORGA = org_in - 0.9 * ORGA - UOA;
dxdt_PHOSE = KPO4 * 0.6 - GFRL * cPhos * 0.15 - 2.5 * fUP * cPhos;

double kset = K0 + BINS_K * (fKSH0 - fKSH) / fKSH0
            + BPH_K * (7.40 - pH) + BOSM_K * (OSM_EFF - 285.0);
double Jkout = GK * (kset - cK);
dxdt_KE = Jkout + k_in - UK;
dxdt_KI = -Jkout;

// ===========================================================================
//  9.  WATER
// ===========================================================================
double thirst = (OSM_EFF - OSM_THIRST) / 12.0;
if (thirst < 0.0) thirst = 0.0; if (thirst > 1.0) thirst = 1.0;
double vomit = KVOMIT * cKet / (cKet + KMVOMIT);
double mm = (MENT > 0.0 ? MENT : 0.0);
double PO = POMAX * WATER * thirst * (1.0 - vomit) * mm;
double insens = INSENS * (1.0 + 0.5 * ILL)
              * (1.0 + 0.35 * (41.5 - PCO2 > 0.0 ? (41.5 - PCO2) / 25.0 : 0.0));

double osmi = OSMBASE * VI0d + 2.0 * (KI - KIC0 * VI0d) + OSMB * VI0d;
double VI_target = osmi / (OSM_TOT > 150.0 ? OSM_TOT : 150.0);
double Jshift = (VI_target - VI) / TSHIFT;                 // L/h ECF -> ICF
dxdt_VECF = RATE_FL + PO + BICARB - UV - insens - Jshift;
dxdt_VICF = Jshift;

// ===========================================================================
// 10.  RENAL FUNCTION, UREA, CREATININE
// ===========================================================================
double gvr = (vrel - VREL_CRIT) / (1.0 - VREL_CRIT); if (gvr < 0.0) gvr = 0.0;
double gfr_t = pow(gvr, 1.5); if (gfr_t > 1.15) gfr_t = 1.15;
dxdt_GFRR = (gfr_t - GFRR) / TGFR;
dxdt_UREA = PRODUREA * (1.0 + 0.8 * ILL + 0.6 * fLIP) - UUREA;
dxdt_CREA = (PCR - GFRL * CREA) / TBW;

// ===========================================================================
// 11.  RESPIRATORY COMPENSATION
//      Winter's formula overestimates compensation at very low bicarbonate;
//      the flatter relation below matches observed DKA gas values.
// ===========================================================================
double pco2_t = 13.0 + 1.15 * HCO3;
if (pco2_t > 46.0) pco2_t = 46.0;
if (pco2_t < PCO2_FLOOR) pco2_t = PCO2_FLOOR;
double fatigue = 1.0 - 0.6 * (1.0 - MENT > 0.0 ? 1.0 - MENT : 0.0);
pco2_t = pco2_t + (1.0 - fatigue) * (PCO2_NORM - pco2_t);
dxdt_PCO2 = (pco2_t - PCO2) / TPCO2;

// ===========================================================================
// 12.  BRAIN
// ===========================================================================
double osmb_t = KOSMB * (OSM_TOT - 292.0 > 0.0 ? OSM_TOT - 292.0 : 0.0);
double tb = (osmb_t > OSMB ? TOSMB_UP : TOSMB_DN);
dxdt_OSMB = (osmb_t - OSMB) / tb;
double vbr_raw = (292.0 + OSMB) / (OSM_TOT > 200.0 ? OSM_TOT : 200.0);
double vbr_t = 1.0 + (vbr_raw - 1.0) * COMPL + KVASO * INJ;
dxdt_VBR = (vbr_t - VBR) / TVBR;
double isch = (25.0 - PCO2 > 0.0 ? (25.0 - PCO2) / 25.0 : 0.0)
            * (0.3 + 0.7 * (1.0 - GFRR > 0.0 ? 1.0 - GFRR : 0.0))
            * (0.3 + 0.7 * (acid_sev < 1.0 ? acid_sev : 1.0));
dxdt_INJ = KINJ * isch * 10.0 - KRES * INJ;

double ment_t = 1.0
  - W_OSM * (OSM_EFF - 310.0 > 0.0 ? OSM_EFF - 310.0 : 0.0)
  - W_PH  * (7.20 - pH > 0.0 ? 7.20 - pH : 0.0)
  - W_VBR * (VBR - 1.010 > 0.0 ? VBR - 1.010 : 0.0)
  - 0.15  * (3.3 - Gp > 0.0 ? 3.3 - Gp : 0.0);
if (ment_t < 0.0) ment_t = 0.0; if (ment_t > 1.0) ment_t = 1.0;
dxdt_MENT = (ment_t - MENT) / TMENT;

// ===========================================================================
// 13.  CUMULATIVE LEDGERS
// ===========================================================================
double frac_nh4 = (UKETr > 1e-9 ? UNH4 / UKETr : 0.0);
if (frac_nh4 > 1.0) frac_nh4 = 1.0;
dxdt_UKET  = UKETr * (1.0 - frac_nh4);       // base LOST with Na+/K+
dxdt_UKETN = UKETr * frac_nh4;               // base KEPT with NH4+
dxdt_UGLU  = UGLUr;
dxdt_UKCUM = UK;
dxdt_UVOL  = UV;
dxdt_CLIN  = cl_in;

$TABLE
capture GLUmM   = Gp;
capture GLUmgdl= Gp * 18.0182;
capture PHa     = pH;
capture BICARB_mM = HCO3;
capture ANIONGAP= AG;
capture NAmM    = cNa;
capture NA_CORR = cNa + 0.024 * (Gp * 18.0182 - 100.0);
capture CLmM    = cCl;
capture KmM     = cK;
capture BHBmM   = cBHB;
capture ACACmM  = cAcAc;
capture KETmM   = cKet;
capture BHB_ACAC= (cAcAc > 1e-9 ? cBHB / cAcAc : 0.0);
capture ACETmM  = ACET / TBW;
capture LACmM   = cLac;
capture PHOSmM  = cPhos;
capture OSMeff  = OSM_EFF;
capture OSMtot  = OSM_TOT;
capture BUN     = cUrea * 2.80;
capture CREA_TRUE = CREA;
// acetoacetate cross-reacts in the alkaline-picrate (Jaffe) creatinine assay
capture CREA_JAFFE = CREA + 0.090 * cAcAc;
capture GFRmlmin= GFRR * GFRMAX * 1000.0 / 60.0;
capture UVLh    = UV;
capture UGLUrate= UGLUr;
capture UKETrate= UKETr;
capture UNArate = UNA;
capture UKrate  = UK;
capture UNH4rate= UNH4;
capture HGPflux = HGP;
capture UPTflux = UPT;
capture KGENflux= KGEN;
capture KOXflux = KOX;
capture LIPflux = LIP;
capture CPT1gate= cpt1;
capture RENAL_CL= (Gp > 1e-9 ? UGLUr / Gp : 0.0);
capture INSPuU  = INSP;
capture INSEFuU = INSEF;
capture IPORTuU = IeffP * IRs;
capture fLIPsupp= 1.0 - fLIP;
capture fUPstim = fUP;
capture KTOT    = KE + KI;
capture POTBICARB = cKet * VKET;      // retained ketoanion = potential base
capture ICPmmHg = 10.0 + KICP * (VBR - 1.0 > 0.0 ? VBR - 1.0 : 0.0);
capture GCS     = 3.0 + 12.0 * MENT;
capture PO_intake = PO;
capture SIDmEq  = SIDapp;
capture VECFL   = VE;
capture TBWL    = TBW;
'

dka <- mcode("dka", code, atol = 1e-8, rtol = 1e-8, maxsteps = 200000)

## =============================================================================
##  FLUID LIBRARY
## =============================================================================
FLUIDS <- list(
  NS         = c(FL_NA = 154, FL_CL = 154, FL_K = 0, FL_ORG = 0,  FL_GLC = 0,   FL_LAC = 0),
  HALF_NS    = c(FL_NA = 77,  FL_CL = 77,  FL_K = 0, FL_ORG = 0,  FL_GLC = 0,   FL_LAC = 0),
  PLASMALYTE = c(FL_NA = 140, FL_CL = 98,  FL_K = 5, FL_ORG = 50, FL_GLC = 0,   FL_LAC = 0),
  LR         = c(FL_NA = 130, FL_CL = 109, FL_K = 4, FL_ORG = 0,  FL_GLC = 0,   FL_LAC = 28),
  D5NS       = c(FL_NA = 154, FL_CL = 154, FL_K = 0, FL_ORG = 0,  FL_GLC = 278, FL_LAC = 0),
  D5HALF     = c(FL_NA = 77,  FL_CL = 77,  FL_K = 0, FL_ORG = 0,  FL_GLC = 278, FL_LAC = 0),
  D10HALF    = c(FL_NA = 77,  FL_CL = 77,  FL_K = 0, FL_ORG = 0,  FL_GLC = 556, FL_LAC = 0),
  D5PL       = c(FL_NA = 140, FL_CL = 98,  FL_K = 5, FL_ORG = 50, FL_GLC = 278, FL_LAC = 0)
)

fluid_pars <- function(name) as.list(FLUIDS[[name]])

## Parameters with the dimensions of a flux, clearance or volume scale with body
## size; rate constants, concentrations, IC50s and stoichiometries do not.
SIZE_PARAMS <- c("SECMAX", "VINS", "CLINS", "HGP0", "GUCNS", "CLNI", "VMI_UP",
                 "GLYCO0", "LIPMAX", "VFFA", "VMAX_KOX", "KLIN_KOX", "LACP0",
                 "KLACCL", "GFRMAX", "TMGLU", "NH4MAX", "NH4B", "PRODUREA",
                 "PCR", "GK", "INSENS", "POMAX")

patient <- function(...) {
  ov <- list(...)
  base <- as.list(param(dka))
  q <- modifyList(base, ov)
  fac <- q$BW / base$BW
  if (abs(fac - 1) > 1e-9)
    for (k in SIZE_PARAMS) if (is.null(ov[[k]])) q[[k]] <- base[[k]] * fac
  q
}

## =============================================================================
##  LEAD-IN:  the presenting state is COMPUTED, never typed in
##  A healthy steady state with insulin withdrawn and a precipitating illness,
##  integrated for `hours`.  Whatever it lands on is the "arrival" state.
## =============================================================================
leadin <- function(p = patient(), hours = 24, delta = 0.25) {
  pp <- modifyList(p, list(RATE_FL = 0, INS_IV = 0, KCL = 0, KPO4 = 0,
                           BICARB = 0))
  pp <- modifyList(pp, fluid_pars("NS"))
  dka %>% param(pp) %>% mrgsim(end = hours, delta = delta, digits = 8) %>% as_tibble()
}

final_state <- function(sim) {
  cmts <- as.character(cmt(dka))
  as.numeric(sim[nrow(sim), cmts])
}

start_from <- function(sim) {
  cmts <- as.character(cmt(dka))
  setNames(as.list(as.numeric(sim[nrow(sim), cmts])), cmts)
}

## =============================================================================
##  TREATMENT:  a closed-loop protocol runner.
##  The dextrose switch, the insulin taper and the ADA potassium rule are
##  FEEDBACK on the simulated patient, not events on a clock.  The protocol is
##  therefore an output of the model as much as an input to it.
## =============================================================================
run_protocol <- function(p, init, hours = 24, step = 0.25,
                         fluid = "NS", dex_fluid = "D5HALF",
                         first_hour = NULL, maint = NULL,
                         ins = 0.10, ins_taper = 0.05, kcl = 40,
                         bolus = 0, dex_trigger = 13.9,
                         taper_after = 12, k_rule = TRUE, bicarb_mmol = 0) {
  bw <- p$BW
  if (is.null(first_hour)) first_hour <- 0.015 * bw
  if (is.null(maint)) maint <- 0.250 * bw / 70
  st <- init
  out <- list()
  t <- 0
  if (bolus > 0) st$INSP <- st$INSP + bolus * 1e6 / (p$VINS * 1000)
  while (t < hours - 1e-9) {
    ## --- read the patient ------------------------------------------------
    probe <- dka %>% param(modifyList(p, c(list(RATE_FL = 0, INS_IV = 0),
                                          fluid_pars(fluid)))) %>%
      init(st) %>% mrgsim(end = 1e-6, delta = 1e-6) %>% as_tibble()
    g   <- tail(probe$GLUmM, 1)
    kk  <- tail(probe$KmM, 1)

    ## --- decide the next 15 minutes of therapy ---------------------------
    fl <- fluid; rate <- if (t < 1) first_hour else if (t < 4) 2 * maint else
      if (t < taper_after) maint else 0.5 * maint
    iv <- ins * bw; kc <- if (t < 1) 0 else kcl
    if (t >= 1 && g < dex_trigger) {
      fl <- dex_fluid
      rate <- 0.250 * max(0.20, min(1.6, (dex_trigger - g) / 4))
      iv <- ins_taper * bw
    }
    if (g < 3.9) fl <- "D10HALF"
    if (k_rule && kk < 3.3) { iv <- 0; kc <- 40 }
    if (kk > 5.3) kc <- 0
    bic <- if (t < 2 && bicarb_mmol > 0) bicarb_mmol / 2000 else 0

    pp <- modifyList(p, c(list(RATE_FL = rate, INS_IV = iv, KCL = kc,
                              BICARB = bic), fluid_pars(fl)))
    seg <- dka %>% param(pp) %>% init(st) %>%
      mrgsim(start = 0, end = step, delta = step, digits = 8) %>% as_tibble()
    seg$time <- seg$time + t
    out[[length(out) + 1]] <- if (t == 0) seg else seg[-1, ]
    st <- start_from(seg)
    t <- t + step
  }
  bind_rows(out)
}

crossing_time <- function(d, col, thr, above = TRUE) {
  v <- d[[col]]; tt <- d$time
  for (i in 2:length(v)) {
    if (above && v[i - 1] < thr && v[i] >= thr)
      return(tt[i - 1] + (thr - v[i - 1]) * (tt[i] - tt[i - 1]) / (v[i] - v[i - 1]))
    if (!above && v[i - 1] > thr && v[i] <= thr)
      return(tt[i - 1] + (thr - v[i - 1]) * (tt[i] - tt[i - 1]) / (v[i] - v[i - 1]))
  }
  NA_real_
}

## =============================================================================
##  SCENARIOS
## =============================================================================

## --- S0: natural history, untreated -----------------------------------------
sc0_untreated <- function(hours = 48) leadin(patient(), hours)

## --- the reference presentation used by every treatment scenario -------------
present <- function(p = patient(), hours = 24) {
  s <- leadin(p, hours)
  list(sim = s, init = start_from(s), obs = s[nrow(s), ])
}

## --- S1: standard ADA / JBDS protocol ---------------------------------------
sc1_standard <- function() {
  pr <- present()
  run_protocol(patient(), pr$init, hours = 30)
}

## --- S2: fluid alone vs insulin alone vs both -------------------------------
sc2_levers <- function() {
  pr <- present()
  list(
    fluid_only   = run_protocol(patient(), pr$init, hours = 6, ins = 0, ins_taper = 0),
    insulin_only = run_protocol(patient(), pr$init, hours = 6, first_hour = 0,
                                maint = 0, kcl = 0),
    both         = run_protocol(patient(), pr$init, hours = 6)
  )
}

## --- S3: insulin dose-response (the low-dose insulin result) ----------------
sc3_insulin_dose <- function(doses = c(0.025, 0.05, 0.10, 0.14, 0.20, 0.40)) {
  pr <- present()
  lapply(setNames(doses, paste0(doses, " U/kg/h")), function(d)
    run_protocol(patient(), pr$init, hours = 30, ins = d, ins_taper = d / 2))
}

## --- S4: crystalloid choice, with and without acute kidney injury -----------
sc4_fluids <- function() {
  pr <- present()
  list(
    saline        = run_protocol(patient(), pr$init, 30, fluid = "NS"),
    balanced      = run_protocol(patient(), pr$init, 30, fluid = "PLASMALYTE",
                                 dex_fluid = "D5PL"),
    ringers       = run_protocol(patient(), pr$init, 30, fluid = "LR",
                                 dex_fluid = "D5PL"),
    saline_AKI    = run_protocol(patient(GFRMAX = 3.0), pr$init, 30, fluid = "NS"),
    balanced_AKI  = run_protocol(patient(GFRMAX = 3.0), pr$init, 30,
                                 fluid = "PLASMALYTE", dex_fluid = "D5PL")
  )
}

## --- S5: potassium replacement, and what the ADA hold rule buys -------------
sc5_potassium <- function(kcls = c(0, 10, 20, 40, 60)) {
  pr <- present()
  c(lapply(setNames(kcls, paste0("KCl ", kcls, " mmol/L")), function(k)
      run_protocol(patient(), pr$init, 30, kcl = k)),
    lapply(setNames(kcls, paste0("KCl ", kcls, ", no hold rule")), function(k)
      run_protocol(patient(), pr$init, 30, kcl = k, k_rule = FALSE)))
}

## --- S6: fluid rate and the cerebral-oedema question (PECARN FLUID) --------
sc6_fluid_rate <- function() {
  pr <- present()
  bw <- 70
  list(
    slow      = run_protocol(patient(), pr$init, 24, first_hour = 0.010 * bw, maint = 0.125),
    standard  = run_protocol(patient(), pr$init, 24),
    fast      = run_protocol(patient(), pr$init, 24, first_hour = 0.020 * bw, maint = 0.500),
    very_fast = run_protocol(patient(), pr$init, 24, first_hour = 0.030 * bw, maint = 0.750)
  )
}

## --- S7: bicarbonate therapy in severe acidaemia ---------------------------
sc7_bicarbonate <- function() {
  p <- patient(ILL0 = 0.95)
  pr <- present(p, 30)
  list(
    none  = run_protocol(p, pr$init, 24, bicarb_mmol = 0),
    b100  = run_protocol(p, pr$init, 24, bicarb_mmol = 100),
    b200  = run_protocol(p, pr$init, 24, bicarb_mmol = 200)
  )
}

## --- S8: the transition off the infusion -----------------------------------
sc8_transition <- function(sc_units = c(0, 10), delay_h = c(0, 2)) {
  pr <- present()
  base <- run_protocol(patient(), pr$init, 30)
  tclose <- crossing_time(base, "ANIONGAP", 12, above = FALSE)
  out <- list()
  for (u in sc_units) for (d in delay_h) {
    idx <- which.min(abs(base$time - (tclose + d)))
    st <- setNames(as.list(as.numeric(base[idx, as.character(cmt(dka))])),
                   as.character(cmt(dka)))
    st$INSSC <- st$INSSC + u
    pp <- modifyList(patient(), c(list(RATE_FL = 0.15, INS_IV = 0, KCL = 40),
                                 fluid_pars("D5HALF")))
    out[[sprintf("SC %g U, stop +%g h", u, d)]] <-
      dka %>% param(pp) %>% init(st) %>% mrgsim(end = 8, delta = 0.25) %>% as_tibble()
  }
  out
}

## --- S9: subcutaneous rapid-analogue protocol ------------------------------
sc9_subcutaneous <- function() {
  p <- patient(ILL0 = 0.30)
  pr <- present(p, 20)
  st <- pr$init; st$INSSC <- st$INSSC + 0.3 * p$BW
  out <- list(); t <- 0
  for (k in 1:12) {
    pp <- modifyList(p, c(list(RATE_FL = if (t < 1) 0.015 * p$BW else 0.25,
                               INS_IV = 0, KCL = if (t < 1) 0 else 40),
                          fluid_pars("NS")))
    seg <- dka %>% param(pp) %>% init(st) %>%
      mrgsim(start = 0, end = 2, delta = 0.25, digits = 8) %>% as_tibble()
    seg$time <- seg$time + t
    out[[k]] <- if (k == 1) seg else seg[-1, ]
    st <- start_from(seg)
    st$INSSC <- st$INSSC + if (tail(seg$GLUmM, 1) >= 13.9) 0.2 * p$BW else 0.1 * p$BW
    t <- t + 2
  }
  list(sc = bind_rows(out), iv = run_protocol(p, pr$init, 24))
}

## --- S10: the phenotype panel — one model, six diseases -------------------
sc10_phenotypes <- function() {
  cases <- list(
    "DKA, type 1, moderate illness" = list(patient(), 24),
    "DKA, severe sepsis"            = list(patient(ILL0 = 0.95), 24),
    "DKA, child 30 kg"              = list(patient(BW = 30), 24),
    "DKA in pregnancy"              = list(patient(ILL0 = 0.5, SIGO = 3.0,
                                                   ADIPOSE = 1.2), 22),
    "Euglycaemic DKA (SGLT2i)"      = list(patient(BETA = 0.10, SGLT2 = 1,
                                                   GLYCO0 = 60), 36),
    "Alcoholic ketoacidosis"        = list(patient(BETA = 0.02, ALCOHOL = 1.2,
                                                   ILL0 = 0.25, WATER = 0.6,
                                                   GLYCO0 = 15, HGP0 = 24), 30),
    "HHS, early"                    = list(patient(BETA = 0.26, WATER = 0.08,
                                                   ILL0 = 0.5), 48),
    "HHS, advanced"                 = list(patient(BETA = 0.22, WATER = 0.06,
                                                   ILL0 = 0.5), 54),
    "DKA with established AKI"      = list(patient(GFRMAX = 3.0), 24)
  )
  lapply(cases, function(cs) leadin(cs[[1]], cs[[2]]))
}

## --- S11: the renal escape valve — hyperglycaemia as a volume disease -----
sc11_escape_valve <- function(waters = c(1, 0.7, 0.5, 0.35, 0.25, 0.15, 0.05, 0)) {
  lapply(setNames(waters, paste0("WATER = ", waters)),
         function(w) leadin(patient(WATER = w), 36))
}

## --- S12: euglycaemic DKA — the valve held open pharmacologically ---------
sc12_euglycaemic <- function() list(
  off = leadin(patient(BETA = 0.10, GLYCO0 = 60), 36),
  on  = leadin(patient(BETA = 0.10, GLYCO0 = 60, SGLT2 = 1), 36)
)

## =============================================================================
##  SUMMARY HELPERS
## =============================================================================
resolution_table <- function(d) {
  tibble::tibble(
    criterion = c("glucose < 250 mg/dL", "pH > 7.30", "HCO3 >= 15",
                  "HCO3 >= 18", "anion gap <= 12", "BHB < 0.6 mmol/L"),
    hours = c(crossing_time(d, "GLUmgdl", 250, FALSE),
              crossing_time(d, "PHa", 7.30, TRUE),
              crossing_time(d, "BICARB_mM", 15, TRUE),
              crossing_time(d, "BICARB_mM", 18, TRUE),
              crossing_time(d, "ANIONGAP", 12, FALSE),
              crossing_time(d, "BHBmM", 0.6, FALSE))
  ) %>% arrange(hours)
}

presentation_row <- function(d) {
  r <- d[nrow(d), ]
  tibble::tibble(glucose_mgdl = r$GLUmgdl, pH = r$PHa, HCO3 = r$BICARB_mM,
                 anion_gap = r$ANIONGAP, BHB = r$BHBmM, AcAc = r$ACACmM,
                 BHB_AcAc = r$BHB_ACAC, Na = r$NAmM, Na_corrected = r$NA_CORR,
                 K = r$KmM, Cl = r$CLmM, BUN = r$BUN, creat_Jaffe = r$CREA_JAFFE,
                 osm_eff = r$OSMeff, GFR = r$GFRmlmin, GCS = r$GCS,
                 ECF_L = r$VECFL, potential_bicarb = r$POTBICARB)
}

## =============================================================================
##  DEMONSTRATION
## =============================================================================
if (interactive()) {

  cat("\n== The presenting state is an OUTPUT of the model ==\n")
  pr <- present()
  print(t(presentation_row(pr$sim)))

  cat("\n== Standard protocol: resolution criteria, in the order they fall ==\n")
  s1 <- sc1_standard()
  print(resolution_table(s1))
  cat("\nGlucose reaches target roughly ",
      round(crossing_time(s1, "BHBmM", 0.6, FALSE) -
            crossing_time(s1, "GLUmgdl", 250, FALSE), 1),
      " h before the ketosis clears.\n", sep = "")

  cat("\n== Insulin dose-response: the lipolysis arm saturates, disposal does not ==\n")
  s3 <- sc3_insulin_dose()
  print(dplyr::bind_rows(lapply(names(s3), function(n) tibble::tibble(
    dose = n,
    ss_uUmL = max(s3[[n]]$INSPuU),
    tBHB_0.6 = crossing_time(s3[[n]], "BHBmM", 0.6, FALSE),
    tGlu_250 = crossing_time(s3[[n]], "GLUmgdl", 250, FALSE),
    K_nadir = min(s3[[n]]$KmM),
    min_glucose = min(s3[[n]]$GLUmgdl)))))

  cat("\n== One model, several diseases ==\n")
  s10 <- sc10_phenotypes()
  print(dplyr::bind_rows(lapply(names(s10), function(n)
    dplyr::mutate(presentation_row(s10[[n]]), phenotype = n, .before = 1))))

  cat("\n== The renal escape valve ==\n")
  s11 <- sc11_escape_valve()
  print(dplyr::bind_rows(lapply(names(s11), function(n) {
    r <- s11[[n]][nrow(s11[[n]]), ]
    tibble::tibble(arm = n, glucose = r$GLUmgdl, osm_eff = r$OSMeff,
                   GFR = r$GFRmlmin, urine_glucose = r$UGLUrate,
                   BHB = r$BHBmM, GCS = r$GCS)
  })))

  cat("\n== Fluid rate barely moves the brain; the presentation does ==\n")
  s6 <- sc6_fluid_rate()
  print(dplyr::bind_rows(lapply(names(s6), function(n) tibble::tibble(
    arm = n, peak_swelling_pct = 100 * (max(s6[[n]]$ICPmmHg) - 10) / 900,
    peak_ICP = max(s6[[n]]$ICPmmHg), min_GCS = min(s6[[n]]$GCS)))))
}

## =============================================================================
##  CALIBRATION NOTES — what each group of parameters was set against
## =============================================================================
##
##  INSULIN PK.  CL 50 L/h and V 6 L give a plasma half-life of 5 min and a
##  steady state of 20 uU/mL per U/h infused, so 0.1 U/kg/h in a 70 kg adult
##  reaches ~140 uU/mL — the range measured during low-dose infusion protocols
##  (Kitabchi 1976; Alberti 1973).
##
##  THE FOUR HALF-MAXIMAL CONSTANTS are the load-bearing parameters of the whole
##  model and are taken from insulin dose-response studies rather than fitted:
##  suppression of lipolysis is half-maximal near 15 uU/mL, suppression of
##  hepatic glucose output near 30, stimulation of peripheral disposal near 60,
##  and the potassium shift near 25 (Zierler & Rabinowitz 1964; Rizza 1981;
##  Nurjhan 1986; Jensen 1989).  Everything the model says about why low-dose
##  insulin is as effective as high-dose, why glucose normalises long before the
##  ketosis, and why an intermittent subcutaneous protocol works in mild disease
##  is a consequence of the ORDERING of those four numbers.
##
##  KETONE KINETICS.  KGSCALE, VMAX_KOX, KM_KOX and KLIN_KOX are set so that the
##  fasting steady state is a total ketone body concentration near 0.8 mmol/L,
##  the untreated diabetic steady state is 16-17 mmol/L with a BHB:AcAc ratio of
##  8-9, and the treated decline is ~1 mmol/L/h — matching the JBDS target rate
##  of >0.5 mmol/L/h and the ketone turnover measurements of Owen and Hall.
##  The non-saturable KLIN_KOX arm represents concentration-driven cerebral and
##  myocardial oxidation and is what keeps the system bounded.
##
##  RENAL GLUCOSE.  TMGLU/GFRMAX gives a threshold near 10 mmol/L (185 mg/dL)
##  with splay, and FTM_GFR encodes glomerulotubular balance so that reabsorptive
##  capacity falls less than filtration does in a prerenal state.  SGLT2 = 1
##  lowers Tm by 80%, dropping the threshold to ~4 mmol/L, which is what makes
##  euglycaemic ketoacidosis possible.
##
##  ELECTROLYTE DEFICITS.  FENA, FECL and FEK0/FEK_ALDO, with the osmotic-diuresis
##  and non-reabsorbable-anion multipliers, are set so that a 24 h untreated
##  course produces the classical deficits: 6-10 mmol/kg of sodium, 3-5 mmol/kg
##  of chloride, 3-5 mmol/kg of potassium, and ~100 mL/kg of water, at a serum
##  potassium that is normal or high (Kitabchi 2009 consensus statement).
##
##  ACID-BASE.  The Figge weak-acid coefficients are used as published.  The
##  respiratory relation PCO2 = 13 + 1.15*HCO3 replaces Winter's formula, which
##  overestimates compensation at bicarbonate concentrations below ~10 mmol/L;
##  the flatter line reproduces measured DKA gas values (pH 7.24 at HCO3 10 and
##  PCO2 25; pH 7.07 at HCO3 7 and PCO2 21).
##
##  A CAVEAT WORTH READING.  Because [HCO3-] is here a residual, it is the small
##  difference of large numbers and is correspondingly ill-conditioned: a 30%
##  change in the fractional excretion of chloride alone moves the presenting
##  bicarbonate by several mmol/L and the pH by more than 0.1, while the anion
##  gap and the BHB — which involve no such cancellation — barely move.  This is
##  a property of the chemistry rather than of the implementation, and it is the
##  reason the model's bicarbonate and pH predictions deserve less trust than its
##  gap and ketone predictions.  It is also, arguably, the reason the anion gap
##  and the point-of-care ketone are the robust variables at the bedside.
## =============================================================================
