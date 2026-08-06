## =====================================================================
##  Congenital Hyperinsulinism (CHI) — mrgsolve QSP model
##  ---------------------------------------------------------------------
##  THE ORGANISING IDEA
##  -------------------
##  The disease is written as ONE number: g, the fraction of beta-cell
##  K_ATP surface conductance that survives.  g is not fitted to any CHI
##  patient.  It enters the model in exactly one place —
##
##        G_KATP = g * gKmax * P_open(ATP/ADP)
##
##  — and everything else follows as arithmetic on a conductance divider:
##
##        V_m = (G_KATP + G_GIRK)*E_K + g_leak*E_leak
##              ---------------------------------------
##                  G_KATP + G_GIRK + g_leak
##
##  Because V_m is a DIVIDER, the three classes of drug separate by
##  themselves, and this is the model's central claim:
##
##    * DIAZOXIDE is MULTIPLICATIVE in g.  It shifts the channel's
##      ATP/ADP set-point (R50) rightward, but the term it acts on is
##      g * gKmax * P_open.  If g = 0 there is nothing to open and the
##      predicted effect is EXACTLY zero — not small, zero.
##    * OCTREOTIDE is ADDITIVE in conductance.  SSTR2 opens a SECOND,
##      G-protein-gated K+ channel that does not carry the mutation, so
##      G_GIRK enters the same divider independently of g and still
##      hyperpolarises the cell at g = 0.
##    * GLUCAGON, DEXTROSE, ERSODETUG and SURGERY are g-INDEPENDENT:
##      they act on substrate supply, on the insulin receptor, or on
##      beta-cell MASS (B), which is a different symbol entirely.
##
##  The second, independent axis is that insulin does not merely lower
##  glucose — it deletes the brain's backup fuel.  Because the insulin
##  IC50 for lipolysis (12 uU/mL) lies BELOW the IC50 for glycogenolysis
##  (30) and gluconeogenesis (45), any insulin level high enough to cause
##  hypoglycaemia has already switched ketogenesis off.  Cerebral fuel is
##  therefore written as a SUM (glucose + ketone + lactate) against a
##  demand, which is what makes the CHI-specific glucose target fall out
##  of arithmetic rather than assertion.
##
## =====================================================================
##  WHAT WAS CALIBRATED (and what therefore is NOT a prediction)
##  ---------------------------------------------------------------------
##  Nine numbers are spent, and EIGHT of them are normal-neonate
##  physiology, not CHI:
##
##   N1  basal insulin 5 uU/mL at glucose 75 mg/dL           -> Smax_ins
##   N2..N5  the SHAPE of the normal glucose->insulin dose-response,
##       expressed as secretion ratios relative to glucose 75 mg/dL:
##           0.25x at 45 · 1.8x at 90 · 12x at 150 · 24x at 250
##       -> fitted (R50, nR, gKmax, kAA) by grid search; achieved
##          0.254 / 2.04 / 12.76 / 23.25
##   N6  neonatal cerebral glucose oxidation 4.2 mg/kg/min at G=75,
##       BBB GLUT1 Km 40 mg/dL, cerebral demand 4.0 mg/kg/min
##   N7  term hepatic glycogen ~2000 mg/kg (achieved 1885 on q3h feeds)
##   N8  insulin IC50 separation lipolysis 12 < ketogenesis 15 <
##       glycogenolysis 30 < gluconeogenesis 45 uU/mL (clamp literature)
##   C1  ONE CHI number: gGIRK_max = 1.2, set so that octreotide
##       30 ug/kg/day reduces the dextrose requirement of severe diffuse
##       CHI by roughly half (achieved 64%).  Without this the octreotide
##       arm would be unanchored.
##
##  Diazoxide potency (kdzx = 0.60) was set on a NORMAL beta-cell, not on
##  a CHI one: 10 mg/kg/day must produce mild iatrogenic hyperglycaemia
##  rather than abolish secretion (achieved mean glucose 134 mg/dL).
##
##  EVERYTHING BELOW IS THEREFORE A PREDICTION
##  ---------------------------------------------------------------------
##  All values verified against an independent pure-Python RK4
##  re-implementation of these same 36 ODEs (see README).
##
##  P1  Severe recessive diffuse CHI (g=0.02) needs 11.2 mg/kg/min of IV
##      dextrose on top of 4.7 mg/kg/min enteral = 15.9 mg/kg/min total
##      to hold 70 mg/dL, with insulin 83 uU/mL and BOHB 0.00 mM.
##      Observed in severe diffuse CHI: 15-20 mg/kg/min.
##  P2  Diazoxide 15 mg/kg/day reduces that requirement by 0.1%.  Across
##      the g axis the reduction is 100% (g>=0.3), 49% (g=0.2), 4.9%
##      (g=0.1), 1.1% (g=0.05), 0.6% (g=0.02) — i.e. the clinical fact
##      that biallelic recessive K_ATP CHI is diazoxide-unresponsive
##      while GDH-HI / HNF4A / dominant K_ATP are not, is ARITHMETIC.
##  P3  Octreotide 30 ug/kg/day still delivers 59% at g=0.02, where
##      diazoxide delivers 0.6%.  Two K-channel drugs, opposite verdicts,
##      one equation.
##  P4  SEVERITY SATURATES below g ~ 0.1: the dextrose requirement is
##      11.07 at g=0.10 and 11.21 at g=0.02, because V_m is already
##      clamped at the leak potential.  Falsifiable prediction: there
##      should be NO genotype-severity gradient among null and severely
##      hypomorphic K_ATP alleles, which is what is observed.
##  P5  Ketone-free brain fails at 43.7 mg/dL; at BOHB 2 mM it tolerates
##      25.3 mg/dL.  Hyperinsulinism therefore costs 18.4 mg/dL of
##      cerebral fuel margin BEFORE any glucose is lost — the reason the
##      CHI target is 70 mg/dL when ketotic hypoglycaemia is tolerated
##      into the 40s.
##  P6  Glucagon stimulation test: +62 mg/dL in CHI (glycogen 2599 mg/kg
##      preserved) vs +16 mg/dL after a 20 h normal fast (glycogen 2).
##      The test's discriminating power comes from glycogen being a STATE
##      VARIABLE that insulin protects.
##  P7  Resection is a narrow window.  Remnant 0.30-0.20 of beta-cell
##      mass = euglycaemic (mean 78-97 mg/dL); remnant <=0.10 = diabetic
##      (177-342 mg/dL); remnant 0.50 = still dextrose-dependent.  g is
##      unchanged in the remnant, so the SAME equation gives both the
##      persistent hypoglycaemia and the surgical diabetes.
##  P8  GDH-HI: leucine 150 mg/kg drops glucose 37 mg/dL vs 16 in
##      dominant K_ATP; ammonia 150 umol/L and EXACTLY unchanged by
##      diazoxide (150.0 -> 150.0), because the hepatic GDH flux is not
##      glucose-linked.  Falsifiable and clinically observed.
##  P9  Nifedipine gives 5.8%: the mechanism is real but its beta-cell
##      EC50 (~1200 ng/mL) sits ~8x above the antihypertensive Cmax, so
##      the window is closed by vasodilatation, not by biology.
##  P10 A FOCAL lesion reproduces severe disease only when
##      (lesion fraction x local beta-cell density) approaches 1.0 —
##      i.e. a small lesion must carry ~10x normal secretory density.
##      Severity depends ONLY on that product (0.10x5 and 0.05x10 both
##      give 4.2-4.3 mg/kg/min).  This is a quantitative demand on the
##      histology of focal adenomatous hyperplasia.
##
##  KNOWN TENSIONS (stated, not hidden)
##  ---------------------------------------------------------------------
##  T1  Predicted plasma insulin in untreated severe diffuse CHI is
##      ~83 uU/mL, ABOVE the 10-50 uU/mL usually reported in critical
##      samples.  This is forced arithmetic, not a free choice: a
##      dextrose requirement of 16 mg/kg/min cannot be reconciled with
##      20 uU/mL given neonatal insulin sensitivity.  Either the reported
##      critical-sample insulins understate the true 24 h exposure
##      (single timed samples, variable assay cross-reactivity), or
##      neonatal insulin sensitivity is higher than modelled.  This is
##      the model's most exposed parameter.
##  T2  Normal fasting glucose plateaus at ~72 mg/dL at 18-24 h, higher
##      than the 55-65 mg/dL commonly measured, because gluconeogenesis
##      has no substrate ceiling here.  This makes the model CONSERVATIVE
##      about normal fasting hypoglycaemia and does not affect the CHI
##      arms, which are insulin-driven.
##  T3  Ersodetug 9 mg/kg gives 75% — a larger effect than published
##      phase-2 experience.  Emax/EC50 for an allosteric insulin-receptor
##      antibody are not publicly established; treat this arm as
##      structural (it shows a g-INDEPENDENT mechanism) not quantitative.
##  T4  Sirolimus (19%) is modelled as a mass/secretion effect because
##      its mechanism in CHI is genuinely unresolved; the number carries
##      no more weight than that assumption.
##
##  DISCLAIMER: educational / research QSP model.  Not validated for
##  clinical use.  Do not use for dosing decisions.
## =====================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

chi_code <- '
$PARAM @annotated
// ---------------- beta-cell metabolic signal (ATP/ADP proxy) ----------
Rbas    : 1.0   : basal ATP/ADP signal (-)
Rmax    : 8.0   : maximal glucose-driven ATP/ADP signal (-)
KG      : 8.0   : glucokinase S0.5 (mM)
hG      : 1.7   : glucokinase Hill coefficient (-)
kAA     : 2.0   : amino-acid (GDH) anaplerotic gain (-)
KAA     : 0.5   : leucine Km for GDH activation (mM)
aGDH    : 3.0   : maximal GDH disinhibition when GTP brake is lost (-)
kMCT1   : 1.5   : beta-cell MCT1 lactate gain if SLC16A1 active (per mM)

// ---------------- K_ATP channel : THE DISEASE AXIS --------------------
R50     : 2.0   : ATP/ADP giving half-maximal channel closure (-)
nR      : 2.5   : channel closure Hill coefficient (-)
gKmax   : 12.0  : maximal K_ATP conductance relative to leak (-)
EK      : -75.0 : potassium reversal potential (mV)
Eleak   : -20.0 : background leak reversal potential (mV)
gleak   : 1.0   : background depolarising leak conductance (-)
V50     : -50.0 : Cav half-activation potential (mV)
kV      : 4.0   : Cav activation slope (mV)
tauCa   : 0.02  : cytosolic Ca time constant (h)

// ---------------- amplifying pathway / secretion ----------------------
ampmin  : 0.35  : floor of the amplifying pathway at low glucose (-)
KA      : 7.0   : glucose S0.5 of the amplifying pathway (mM)
hA      : 1.5   : amplifying pathway Hill coefficient (-)
acamp   : 0.6   : cAMP amplification gain (-)
Smax_ins: 1290  : maximal insulin secretion scale (uU/mL/h)
Sleak   : 0.0015: Ca-independent basal secretion (-)
kelI    : 8.318 : insulin elimination rate, t1/2 5 min (1/h)
aCP     : 0.0333: C-peptide secretion scale (ng/mL per unit)
kelCP   : 1.386 : C-peptide elimination, t1/2 30 min (1/h)
kx      : 6.0   : insulin action equilibration rate (1/h)

// ---------------- SSTR2 / octreotide ---------------------------------
gGIRK   : 1.2   : maximal GIRK conductance from SSTR2 (-)  [CHI-calibrated]
EC50oct : 1.0   : octreotide EC50 (ng/mL)
kdesen  : 0.010 : SSTR2 desensitisation rate (1/h)
krecov  : 0.020 : SSTR2 resensitisation rate (1/h)

// ---------------- whole-body glucose (fluxes in mg/kg/min) -----------
Vg      : 0.22  : glucose distribution volume (L/kg)
UbrMax  : 6.44  : maximal cerebral glucose uptake (mg/kg/min)
Kbr     : 40.0  : BBB GLUT1 Km (mg/dL)
acbf    : 0.35  : maximal CBF autoregulatory gain (-)
Gcbf    : 70.0  : glucose below which CBF starts to rise (mg/dL)
CMRreq  : 4.00  : cerebral energy demand, glucose-equivalents (mg/kg/min)
VketMax : 3.0   : maximal cerebral ketone oxidation (mg/kg/min glc-eq)
Kket    : 6.0   : brain MCT1 Km for BOHB (mM)
VlacMax : 0.8   : maximal cerebral lactate oxidation (mg/kg/min glc-eq)
Klac    : 4.0   : brain lactate Km (mM)
Uii0    : 2.0   : insulin-independent non-cerebral uptake scale (mg/kg/min)
KUii    : 90.0  : Km of insulin-independent uptake (mg/dL)
kren    : 0.02  : renal glucose loss above threshold (mg/kg/min per mg/dL)
Uidmax  : 20.0  : maximal insulin-dependent disposal (mg/kg/min)
KI      : 40.0  : insulin EC50 for peripheral disposal (uU/mL)
nI      : 2.0   : Hill coefficient for insulin-dependent disposal (-)
G0id    : 90.0  : reference glucose for mass-action disposal (mg/dL)

// ---------------- liver ----------------------------------------------
kglyMax : 3.60  : maximal glycogenolysis (mg/kg/min)
KGLY    : 800.0 : glycogen Km for glycogenolysis (mg/kg)
kgngMax : 3.00  : maximal gluconeogenesis (mg/kg/min)
IC50gly : 30.0  : insulin IC50 for glycogenolysis (uU/mL)
IC50gng : 45.0  : insulin IC50 for gluconeogenesis (uU/mL)
ng      : 2.0   : insulin Hill coefficient on hepatic output (-)
ksyn    : 2.5   : glycogen synthesis gain on portal glucose (-)
KsynI   : 10.0  : insulin EC50 for glycogen synthesis (uU/mL)
Ksyng   : 100.0 : glucose EC50 for glycogen synthesis (mg/dL)
GLYmax  : 4000  : maximal hepatic glycogen (mg/kg)
kabs    : 0.012 : enteral glucose absorption rate (1/min)

// ---------------- lipolysis / ketogenesis ----------------------------
kLipo   : 7.25  : maximal lipolysis (mmol/L/h)
IC50lip : 12.0  : insulin IC50 for lipolysis (uU/mL)
nl      : 1.8   : lipolysis Hill coefficient (-)
kFFAox  : 10.0  : FFA clearance (1/h)
aepiL   : 0.5   : epinephrine gain on lipolysis (-)
kket    : 0.4167: ketogenesis rate constant (1/h per mM FFA)
IC50ket : 15.0  : insulin IC50 for ketogenesis (uU/mL)
nk      : 2.0   : ketogenesis Hill coefficient (-)
kketox  : 1.5   : BOHB clearance (1/h)
agcgK   : 3.0   : glucagon gain on ketogenesis (-)
akg     : 2.5   : gain of the glycogen-depletion ketogenic switch (-)

// ---------------- counterregulation ----------------------------------
kgcgS   : 563.0 : glucagon secretion scale (pg/mL/h)
kgcgE   : 5.199 : glucagon elimination, t1/2 8 min (1/h)
agcg    : 3.0   : hypoglycaemic gain on glucagon (-)
Ghalf   : 70.0  : glucose midpoint of glucagon response (mg/dL)
kgs     : 8.0   : slope of glucagon glucose response (mg/dL)
IC50a   : 60.0  : insulin IC50 for alpha-cell suppression (uU/mL)
GCG0    : 100.0 : basal glucagon (pg/mL)
KGCG    : 400.0 : glucagon receptor Kd (pg/mL)
ashift  : 6.0   : cAMP-driven rightward shift of the insulin IC50 (-)
fgcgconv: 1000  : ug/kg -> pg/mL conversion for exogenous glucagon (-)
kepiS   : 695.0 : epinephrine secretion scale (pg/mL/h)
kepiE   : 13.90 : epinephrine elimination (1/h)
aepi    : 8.0   : hypoglycaemic gain on epinephrine (-)
Gepi    : 60.0  : glucose midpoint of adrenal response (mg/dL)
kes     : 6.0   : slope of adrenal response (mg/dL)
EPI0    : 50.0  : basal epinephrine (pg/mL)
kcortS  : 3.70  : cortisol secretion scale (ug/dL/h)
kcortE  : 0.462 : cortisol elimination (1/h)
acort   : 2.0   : hypoglycaemic gain on cortisol (-)
CORT0   : 8.0   : basal cortisol (ug/dL)

// ---------------- ammonia / amino acids / lactate --------------------
kNH3    : 40.0  : ammonia production scale (umol/L/h)
kNH3e   : 1.0   : ammonia clearance (1/h)
aNH3    : 0.917 : GDH gain on ammonia production (-)
kAAel   : 1.2   : leucine elimination (1/h)
AA0     : 0.12  : basal plasma leucine (mmol/L)
klacp   : 1.30  : lactate production scale (mmol/L/h)
klacc   : 1.0   : lactate clearance (1/h)
LAC0    : 1.3   : reference lactate (mmol/L)

// ---------------- drug PK / PD ---------------------------------------
ka_dzx  : 0.60  : diazoxide absorption (1/h)
Vc_dzx  : 0.25  : diazoxide central volume (L/kg)
CL_dzx  : 0.00866: diazoxide clearance, t1/2 ~20 h (L/kg/h)
Q_dzx   : 0.010 : diazoxide intercompartmental clearance (L/kg/h)
Vp_dzx  : 0.20  : diazoxide peripheral volume (L/kg)
Emax_dzx: 0.95  : maximal diazoxide occupancy (-)
EC50_dzx: 25.0  : diazoxide EC50 (ug/mL)
kdzx    : 0.60  : rightward shift of R50 at full occupancy (-)
ka_oct  : 1.80  : octreotide SC absorption (1/h)
V_oct   : 0.27  : octreotide volume (L/kg)
CL_oct  : 0.110 : octreotide clearance, t1/2 ~1.7 h (L/kg/h)
ka_gcg  : 6.0   : glucagon SC absorption (1/h)
V_gcg   : 0.25  : glucagon volume (L/kg)
Vc_ers  : 0.055 : ersodetug central volume (L/kg)
CL_ers  : 0.000160: ersodetug clearance, t1/2 ~2 wk (L/kg/h)
Q_ers   : 0.00020: ersodetug intercompartmental clearance (L/kg/h)
Vp_ers  : 0.045 : ersodetug peripheral volume (L/kg)
Emax_ers: 0.70  : maximal insulin-receptor blockade (-)
EC50_ers: 45.0  : ersodetug EC50 (ug/mL)
ka_sir  : 1.0   : sirolimus absorption (1/h)
V_sir   : 6.0   : sirolimus volume (L/kg)
CL_sir  : 0.070 : sirolimus clearance (L/kg/h)
Emax_sir: 0.35  : maximal mTOR effect on secretion (-)
EC50_sir: 0.012 : sirolimus EC50 (mg/L)
V_nif   : 1.2   : nifedipine volume (L/kg)
CL_nif  : 0.35  : nifedipine clearance (L/kg/h)
Emax_nif: 0.60  : maximal Cav blockade (-)
EC50_nif: 1.20  : nifedipine beta-cell EC50 (mg/L) -- far above plasma Cmax
CL_ex9  : 1.20  : exendin(9-39) clearance (L/kg/h)
V_ex9   : 0.20  : exendin(9-39) volume (L/kg)
Emax_ex9: 0.55  : maximal GLP-1R blockade (-)
EC50_ex9: 30.0  : exendin(9-39) EC50 (ng/mL)

// ---------------- genotype / disease switches ------------------------
g_ab    : 1.00  : residual K_ATP conductance of the ABNORMAL population
w_ab    : 0.00  : fraction of beta-cell mass that is abnormal (1=diffuse)
dens_a  : 1.00  : local beta-cell density of the abnormal clone (x normal)
sGTP    : 1.00  : residual GTP inhibition of GDH (0 = GDH-HI)
KGshift : 1.00  : multiplier on glucokinase S0.5 (<1 = activating GCK)
mct1    : 0.00  : beta-cell MCT1 expression (1 = SLC16A1 exercise-induced)
BMASS0  : 1.00  : beta-cell mass at time 0 (1 = intact pancreas)
kprol   : 0.0010: beta-cell proliferation (1/h)
kapo    : 0.00050: beta-cell apoptosis (1/h)
Bmax    : 1.30  : beta-cell mass carrying capacity (-)
kdev    : 1.0   : neurodevelopmental deficit accrual scale (1/h)

// ---------------- interventions --------------------------------------
GIR_fix : 0.0   : fixed IV dextrose rate (mg/kg/min)
feed    : 0.0   : continuous enteral glucose (mg/kg/h)
GCG_inf : 0.0   : continuous glucagon infusion (ug/kg/h)
EX9inf  : 0.0   : exendin(9-39) infusion (ug/kg/h)
loop    : 0.0   : 1 = closed-loop dextrose controller ON
Gtarget : 70.0  : closed-loop glucose target (mg/dL)
Kp      : 0.10  : proportional gain (mg/kg/min per mg/dL)
Ki      : 0.60  : integral gain (mg/kg/min per mg/dL per h)
GIRmax  : 40.0  : maximum deliverable dextrose (mg/kg/min)

$INIT @annotated
GLU     :  75.0 : plasma glucose (mg/dL)
GLUi    :  75.0 : interstitial (CGM) glucose (mg/dL)
INS     :   5.0 : plasma insulin (uU/mL)
X       :   5.0 : insulin at the receptor / effect site (uU/mL)
CPEP    :   1.0 : plasma C-peptide (ng/mL)
GLY     : 2000.0: hepatic glycogen (mg/kg)
FFA     :   0.60: plasma free fatty acids (mmol/L)
BOHB    :   0.15: beta-hydroxybutyrate (mmol/L)
LAC     :   1.30: plasma lactate (mmol/L)
GCG     : 100.0 : plasma glucagon (pg/mL)
EPI     :  50.0 : plasma epinephrine (pg/mL)
CORT    :   8.0 : plasma cortisol (ug/dL)
NH3     :  40.0 : plasma ammonia (umol/L)
CABn    :   0.04: Ca signal, NORMAL beta-cell population (-)
CABa    :   0.04: Ca signal, ABNORMAL beta-cell population (-)
CAMP    :   1.0 : beta-cell cAMP (relative to basal)
BMASS   :   1.0 : beta-cell mass (fraction of normal)
AA      :   0.12: plasma leucine (mmol/L)
GGUT    :   0.0 : enteral glucose in gut lumen (mg/kg)
DZXg    :   0.0 : diazoxide gut depot (mg/kg)
DZX     :   0.0 : diazoxide central (ug/mL)
DZXp    :   0.0 : diazoxide peripheral (ug/mL)
OCTs    :   0.0 : octreotide SC depot (ug/kg)
OCT     :   0.0 : octreotide plasma (ng/mL)
ROCT    :   1.0 : SSTR2 responsiveness (1 = naive)
GCGs    :   0.0 : exogenous glucagon SC depot (ug/kg)
ERS     :   0.0 : ersodetug central (ug/mL)
ERSp    :   0.0 : ersodetug peripheral (ug/mL)
SIRg    :   0.0 : sirolimus gut depot (mg/kg)
SIR     :   0.0 : sirolimus plasma (mg/L)
NIF     :   0.0 : nifedipine plasma (mg/L)
EX9     :   0.0 : exendin(9-39) plasma (ng/mL)
AUCHYPO :   0.0 : cumulative area below 70 mg/dL (mg/dL*h)
TFUEL   :   0.0 : cumulative time in cerebral fuel deficit (h)
DEV     :   0.0 : neurodevelopmental deficit index (-)
GIRi    :   0.0 : dextrose controller integral (mg/kg/min)

$GLOBAL
#define SIGLOW(x, x50, k) (1.0/(1.0+exp(((x)-(x50))/(k))))

// $ODE locals are scoped to the ODE function, so anything $TABLE needs to
// report has to live at file scope.  These are written by $ODE and read by
// $TABLE at each output time.
double gVa, gVn, gGIR, gRa, gFuel, gVgly, gVgng;

$ODE
// ================= 1. BETA-CELL: g becomes a voltage =================
double Gb   = fmax(GLU, 1e-6)/18.016;               // plasma glucose in mM
double fGDH = aGDH*(1.0 - sGTP);                    // GDH disinhibition
double gdh  = kAA*(AA/(KAA + AA))*(1.0 + fGDH);
double KGe  = KG*KGshift;                           // GCK-activating shifts S0.5
double Rt   = Rbas + (Rmax - Rbas)*pow(Gb,hG)/(pow(KGe,hG) + pow(Gb,hG))
              + gdh + kMCT1*mct1*LAC;

// diazoxide is a K_ATP OPENER: it shifts the ATP/ADP set-point rightward.
// Note it acts on P_open, which is MULTIPLIED by g -> zero effect at g=0.
double Fdzx = Emax_dzx*DZX/(EC50_dzx + DZX);
double R50e = R50*(1.0 + kdzx*Fdzx);
double Po   = 1.0/(1.0 + pow(Rt/R50e, nR));

// octreotide opens a SECOND, g-INDEPENDENT K conductance (GIRK)
double Foct  = (OCT/(EC50oct + OCT))*ROCT;
double GGIRK = gGIRK*Foct;

double GKn = gKmax*1.0  *Po;                        // normal population
double GKa = gKmax*g_ab*Po;                         // abnormal population
double Vn  = ((GKn + GGIRK)*EK + gleak*Eleak)/(GKn + GGIRK + gleak);
double Va  = ((GKa + GGIRK)*EK + gleak*Eleak)/(GKa + GGIRK + gleak);

double Fnif = Emax_nif*NIF/(EC50_nif + NIF);
double fCan = (1.0 - Fnif)/(1.0 + exp(-(Vn - V50)/kV));
double fCaa = (1.0 - Fnif)/(1.0 + exp(-(Va - V50)/kV));

double amp   = ampmin + (1.0 - ampmin)*pow(Gb,hA)/(pow(KA,hA) + pow(Gb,hA));
double Fsir  = Emax_sir*SIR/(EC50_sir + SIR);
double fcamp = 1.0 + acamp*(CAMP - 1.0);
double massT = (1.0 - w_ab)*CABn + w_ab*dens_a*CABa;
double Ssec  = Smax_ins*BMASS*(massT*amp*fcamp + Sleak)*(1.0 - Fsir);

double Fers = Emax_ers*ERS/(EC50_ers + ERS);        // insulin-receptor blockade

// ================= 2. HEPATIC OUTPUT =================================
double fgcg  = GCG/GCG0;
double Egcg  = GCG/(KGCG + GCG);                    // glucagon occupancy
double Egcg0 = GCG0/(KGCG + GCG0);
// PKA (glucagon) and PP1 (insulin) compete on phosphorylase, so raising
// cAMP shifts the insulin IC50 for glycogenolysis to the right.  This is
// what lets a pharmacological glucagon bolus overcome the insulin block
// and makes the glucagon stimulation test work.
double IC50ge = IC50gly*(1.0 + ashift*fmax(0.0, Egcg - Egcg0)/fmax(1e-9, 1.0 - Egcg0));
double finsG  = 1.0/(1.0 + pow(X/IC50ge, ng));
double finsN  = 1.0/(1.0 + pow(X/IC50gng, ng));
double Vgly   = kglyMax*(GLY/(KGLY + GLY))*(Egcg/Egcg0)
                *(1.0 + 0.5*(EPI/EPI0 - 1.0))*finsG;
double Vgng   = kgngMax*(0.4 + 0.6*fmin(fgcg,3.0))*finsN
                *(1.0 + 0.25*(CORT/CORT0 - 1.0))*pow(fmax(LAC,0.1)/LAC0, 0.3);
double EGP    = fmax(Vgly + Vgng, 0.0);
double Ra     = kabs*GGUT;                          // enteral appearance

// ================= 3. CEREBRAL FUEL: a SUM, not a glucose ============
double Ubr   = UbrMax*GLU/(Kbr + GLU)
               *(1.0 + acbf*fmax(0.0, Gcbf - GLU)/Gcbf);   // CBF autoregulation
double Uket  = VketMax*BOHB/(Kket + BOHB);
double Ulac  = VlacMax*LAC/(Klac + LAC);
double FUEL  = Ubr + Uket + Ulac;
double fuelr = FUEL/CMRreq;                         // <1 = neuroglycopenia

// ================= 4. PERIPHERAL DISPOSAL ============================
double Uii  = Uii0*GLU/(KUii + GLU);
double Uren = kren*fmax(0.0, GLU - 180.0);
double Uid  = Uidmax*pow(X,nI)/(pow(KI,nI) + pow(X,nI))*(GLU/G0id);

// ================= 5. DEXTROSE (fixed or closed loop) ================
double GIR = GIR_fix;
if (loop > 0.5) GIR = fmax(0.0, fmin(GIRmax, Kp*(Gtarget - GLU) + GIRi));

// ================= 6. ODEs ===========================================
double kconv = 6.0/Vg;      // mg/dL per h per (mg/kg/min)

dxdt_GLU  = kconv*(EGP + Ra + GIR - Uii - Uid - Uren - Ubr);
dxdt_GLUi = (GLU - GLUi)/0.17;
dxdt_INS  = Ssec - kelI*INS;
dxdt_X    = kx*(INS*(1.0 - Fers) - X);
dxdt_CPEP = aCP*Ssec - kelCP*CPEP;

// glycogen synthesis is driven by PORTAL glucose delivery, gated by
// insulin -- so it is ~0 in a true fast even though basal insulin is not.
double Vsyn = ksyn*Ra*(X/(KsynI + X))*(GLU/(Ksyng + GLU))
              *fmax(0.0, 1.0 - GLY/GLYmax);
dxdt_GLY  = 60.0*(Vsyn - Vgly);

double fl = 1.0/(1.0 + pow(X/IC50lip, nl));
dxdt_FFA  = kLipo*fl*(1.0 + aepiL*(EPI/EPI0 - 1.0)) - kFFAox*FFA;

// two independent reasons CHI is ketone-free: insulin acts directly on
// ketogenesis (IC50 15), AND the hepatic fasting switch never arms
// because insulin stops the glycogen from ever being spent.
double fk     = 1.0/(1.0 + pow(X/IC50ket, nk));
double switch_= 1.0 + akg*fmax(0.0, 1.0 - GLY/2000.0);
dxdt_BOHB = kket*FFA*fk*(1.0 + agcgK*(fgcg - 1.0))*switch_ - kketox*BOHB;

dxdt_LAC  = klacp*(GLU/90.0) - klacc*LAC;

double falpha = 1.0/(1.0 + X/IC50a);                // insulin brakes alpha-cell
dxdt_GCG  = kgcgS*(1.0 + agcg*SIGLOW(GLU,Ghalf,kgs))*falpha - kgcgE*GCG
            + fgcgconv*(ka_gcg*GCGs + GCG_inf)/V_gcg;
dxdt_EPI  = kepiS*(1.0 + aepi*SIGLOW(GLU,Gepi,kes)) - kepiE*EPI;
dxdt_CORT = kcortS*(1.0 + acort*SIGLOW(GLU,55.0,6.0)) - kcortE*CORT;
dxdt_NH3  = kNH3*(1.0 + aNH3*fGDH) - kNH3e*NH3;     // NOT glucose-linked

dxdt_CABn = (fCan - CABn)/tauCa;
dxdt_CABa = (fCaa - CABa)/tauCa;
double Fex9 = Emax_ex9*EX9/(EC50_ex9 + EX9);
dxdt_CAMP = (1.0 - Fex9 - 0.55*Foct - CAMP)/0.05;
dxdt_BMASS= kprol*BMASS*fmax(0.0, 1.0 - BMASS/Bmax) - kapo*BMASS;
dxdt_AA   = -kAAel*(AA - AA0);
dxdt_GGUT = feed - 60.0*kabs*GGUT;

// ---- drug PK ----
dxdt_DZXg = -ka_dzx*DZXg;
dxdt_DZX  = 0.9*ka_dzx*DZXg/Vc_dzx - CL_dzx*DZX/Vc_dzx - Q_dzx*(DZX - DZXp)/Vc_dzx;
dxdt_DZXp = Q_dzx*(DZX - DZXp)/Vp_dzx;
dxdt_OCTs = -ka_oct*OCTs;
dxdt_OCT  = ka_oct*OCTs/V_oct - CL_oct*OCT/V_oct;
dxdt_ROCT = krecov*(1.0 - ROCT) - kdesen*ROCT*Foct;
dxdt_GCGs = -ka_gcg*GCGs;
dxdt_ERS  = -CL_ers*ERS/Vc_ers - Q_ers*(ERS - ERSp)/Vc_ers;
dxdt_ERSp = Q_ers*(ERS - ERSp)/Vp_ers;
dxdt_SIRg = -ka_sir*SIRg;
dxdt_SIR  = ka_sir*SIRg/V_sir - CL_sir*SIR/V_sir;
dxdt_NIF  = -CL_nif*NIF/V_nif;
dxdt_EX9  = (EX9inf - CL_ex9*EX9)/V_ex9;

// ---- outcome integrals ----
dxdt_AUCHYPO = fmax(0.0, 70.0 - GLU);
dxdt_TFUEL   = SIGLOW(fuelr, 1.0, 0.02);
dxdt_DEV     = kdev*fmax(0.0, 1.0 - fuelr);
dxdt_GIRi    = (loop > 0.5) ? Ki*(Gtarget - GLU) : 0.0;

// hand the reported intermediates out to $TABLE
gVa = Va; gVn = Vn; gGIR = GIR; gRa = Ra;
gFuel = fuelr; gVgly = Vgly; gVgng = Vgng;

$TABLE
double Vm_ab   = gVa;
double Vm_norm = gVn;
double GIRout  = gGIR;
double TOTGLC  = gGIR + gRa;               // total glucose delivery
double fuelrat = gFuel;
double ICratio = (CPEP > 1e-6) ? INS/CPEP : 0.0;
double Egly    = gVgly;
double Egng    = gVgng;
double insact  = X;

$CAPTURE @annotated
Vm_ab   : resting membrane potential, abnormal beta-cells (mV)
Vm_norm : resting membrane potential, normal beta-cells (mV)
GIRout  : IV dextrose rate (mg/kg/min)
TOTGLC  : total glucose delivery, IV + enteral (mg/kg/min)
fuelrat : cerebral fuel supply / demand (-)
ICratio : insulin : C-peptide ratio (-)
Egly    : glycogenolysis flux (mg/kg/min)
Egng    : gluconeogenesis flux (mg/kg/min)
insact  : insulin at the receptor (uU/mL)
'

mod <- mcode("chi_qsp", chi_code, atol = 1e-8, rtol = 1e-8, maxsteps = 100000)

## =====================================================================
##  GENOTYPES — every one of these is a value of g (plus, for the
##  non-K_ATP forms, a shift in what feeds the same channel)
## =====================================================================
geno <- list(
  normal            = list(g_ab = 1.00, w_ab = 0),
  katp_dominant     = list(g_ab = 0.60, w_ab = 1),        # ABCC8 dominant
  katp_recess_diff  = list(g_ab = 0.02, w_ab = 1),         # biallelic, diffuse
  katp_focal        = list(g_ab = 0.02, w_ab = 0.10, dens_a = 10),
  gdh_hi            = list(g_ab = 1.00, w_ab = 1, sGTP = 0.0),   # GLUD1 / HI-HA
  gck_activating    = list(g_ab = 1.00, w_ab = 1, KGshift = 0.45),
  hnf4a             = list(g_ab = 0.70, w_ab = 1),
  schad             = list(g_ab = 1.00, w_ab = 1, sGTP = 0.35),
  slc16a1           = list(g_ab = 1.00, w_ab = 1, mct1 = 1.0),
  post_pancreatect  = list(g_ab = 0.02, w_ab = 1, BMASS0 = 0.05)
)

## ---------------------------------------------------------------------
##  Dosing helpers.  Enteral feeds q3h; 845 mg/kg systemic per feed
##  (=4.69 mg/kg/min average) represents ~150 mL/kg/day of milk after
##  splanchnic first-pass extraction.
## ---------------------------------------------------------------------
feeds     <- function(days = 4)  ev(amt = 845, cmt = "GGUT", time = 0,
                                    ii = 3, addl = days*8 - 1)
diazoxide <- function(mgkgday, days = 6) ev(amt = mgkgday/3, cmt = "DZXg",
                                    time = 0, ii = 8, addl = days*3 - 1)
octreotide<- function(ugkgday, days = 6) ev(amt = ugkgday/4, cmt = "OCTs",
                                    time = 0, ii = 6, addl = days*4 - 1)
glucagon_bolus <- function(mgkg = 0.03) ev(amt = mgkg*1000, cmt = "GCGs", time = 0)
ersodetug <- function(mgkg = 9)  ev(amt = mgkg/0.055, cmt = "ERS", time = 0)
sirolimus <- function(mgkgday, days = 10) ev(amt = mgkgday, cmt = "SIRg",
                                    time = 0, ii = 24, addl = days - 1)
nifedipine<- function(mgkg = 0.5, days = 6) ev(amt = mgkg/1.2, cmt = "NIF",
                                    time = 0, ii = 8, addl = days*3 - 1)
leucine   <- function(mmol = 1.10) ev(amt = mmol, cmt = "AA", time = 0)

run_chi <- function(gt = "normal", rx = NULL, end = 72, delta = 0.05,
                    closed_loop = FALSE, target = 70, extra = list()) {
  p <- c(geno[[gt]], extra)
  if (closed_loop) p <- c(p, list(loop = 1, Gtarget = target))
  e <- feeds(ceiling(end/24) + 1)
  if (!is.null(rx)) e <- c(e, rx)
  mod %>% param(p) %>% ev(e) %>% mrgsim(end = end, delta = delta) %>% as_tibble()
}

## =====================================================================
##  SCENARIO 1 — the normal neonate (the calibration target)
## =====================================================================
s1 <- run_chi("normal", end = 72)
## expect: glucose 85-133 mg/dL, insulin 10.6-43.9 uU/mL, glycogen ~1885 mg/kg

## =====================================================================
##  SCENARIO 2 — THE DISEASE AXIS.  Sweep g and read off the dextrose
##  requirement.  This is the model's spine: nothing about CHI was
##  fitted, so the whole severity gradient is a prediction.
##  Verified: g=1.00/0.60/0.40/0.20/0.10/0.02 -> IV dextrose
##            0.00/0.00/1.62/10.05/11.07/11.21 mg/kg/min
##  Note the SATURATION below g~0.1 (11.07 vs 11.21): once V_m is
##  clamped at the leak potential the phenotype cannot get worse.
## =====================================================================
s2 <- lapply(c(1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05, 0.02), function(g) {
  run_chi("katp_recess_diff", end = 48, closed_loop = TRUE,
          extra = list(g_ab = g)) %>%
    filter(time > 36) %>%
    summarise(g = g, GIR = mean(GIRout), TOTAL = mean(TOTGLC),
              glucose = mean(GLU), insulin = mean(INS), BOHB = mean(BOHB),
              Vm = mean(Vm_ab))
}) %>% bind_rows()

## =====================================================================
##  SCENARIO 3 — DIAZOXIDE vs OCTREOTIDE ACROSS g
##  The central pharmacological claim.  Same endpoint, same doses, two
##  drugs whose g-dependence differs because one multiplies g and the
##  other adds to the divider.
##  Verified reduction in dextrose requirement:
##     g      diazoxide 15 mg/kg/d     octreotide 30 ug/kg/d
##    0.30          100 %                    100 %
##    0.20           49 %                    100 %
##    0.10          4.9 %                    100 %
##    0.05          1.1 %                     71 %
##    0.02          0.6 %                     59 %
## =====================================================================
s3 <- lapply(c(0.30, 0.20, 0.10, 0.05, 0.02), function(g) {
  base <- run_chi("katp_recess_diff", end = 72, closed_loop = TRUE,
                  extra = list(g_ab = g))
  dzx  <- run_chi("katp_recess_diff", rx = diazoxide(15), end = 72,
                  closed_loop = TRUE, extra = list(g_ab = g))
  oct  <- run_chi("katp_recess_diff", rx = octreotide(30), end = 72,
                  closed_loop = TRUE, extra = list(g_ab = g))
  m <- function(d) mean(filter(d, time > 60)$GIRout)
  tibble(g = g, GIR0 = m(base),
         dzx_pct = 100*(m(base) - m(dzx))/m(base),
         oct_pct = 100*(m(base) - m(oct))/m(base))
}) %>% bind_rows()

## =====================================================================
##  SCENARIO 4 — DIAZOXIDE TITRATION BY GENOTYPE (open loop)
##  Verified mean glucose (mg/dL) at 0 / 5 / 10 / 15 mg/kg/day:
##    normal            99.8 / 122.2 / 134.5 / 142.1   (iatrogenic hyperglycaemia)
##    KATP dominant     80.4 /  93.5 / 100.3 / 104.5   (normalises)
##    GDH-HI            76.3 /  91.9 / 100.2 / 105.5   (normalises)
##    KATP recessive    31.8 /  31.9 /  31.9 /  31.9   (FLAT — no response)
## =====================================================================
s4 <- lapply(c("normal","katp_dominant","gdh_hi","katp_recess_diff"), function(gt) {
  lapply(c(0, 5, 10, 15), function(d) {
    r <- run_chi(gt, rx = if (d > 0) diazoxide(d) else NULL, end = 72)
    tibble(genotype = gt, dose = d, glucose = mean(filter(r, time > 60)$GLU))
  }) %>% bind_rows()
}) %>% bind_rows()

## =====================================================================
##  SCENARIO 5 — GLUCAGON STIMULATION TEST
##  A positive test (rise >30 mg/dL) is not a quirk: it is proof that
##  the glycogen was never spent, which only insulin can arrange.
##  Verified: CHI +62 mg/dL (glycogen 2599) vs normal after a 20 h fast
##  +16 mg/dL (glycogen 2).  Dextrose must be held FIXED during the test;
##  a closed loop silently cancels the rise.
## =====================================================================
s5_chi <- run_chi("katp_recess_diff", rx = glucagon_bolus(0.03), end = 4,
                  extra = list(GIR_fix = 5.0))

## =====================================================================
##  SCENARIO 6 — AGENT PANEL in severe recessive diffuse CHI (g = 0.02)
##  Verified reduction in IV dextrose requirement:
##    diazoxide 15 mg/kg/d ......  0.1 %   <- the whole point
##    nifedipine 0.5 mg/kg q8h ...  5.8 %   (EC50 unreachable in plasma)
##    sirolimus (trough 16) ...... 19.3 %
##    exendin(9-39) .............. 25.8 %
##    glucagon 15 ug/kg/h ........ 35.4 %
##    octreotide 30 ug/kg/d ...... 63.9 %
##    ersodetug 9 mg/kg .......... 74.7 %
##    octreotide + glucagon ...... 100  %
## =====================================================================
s6 <- list(
  none        = run_chi("katp_recess_diff", end = 48, closed_loop = TRUE),
  diazoxide   = run_chi("katp_recess_diff", rx = diazoxide(15),  end = 48, closed_loop = TRUE),
  octreotide  = run_chi("katp_recess_diff", rx = octreotide(30), end = 48, closed_loop = TRUE),
  glucagon    = run_chi("katp_recess_diff", end = 48, closed_loop = TRUE,
                        extra = list(GCG_inf = 15)),
  ersodetug   = run_chi("katp_recess_diff", rx = ersodetug(9),   end = 48, closed_loop = TRUE),
  sirolimus   = run_chi("katp_recess_diff", rx = sirolimus(0.06),end = 48, closed_loop = TRUE),
  nifedipine  = run_chi("katp_recess_diff", rx = nifedipine(0.5),end = 48, closed_loop = TRUE),
  oct_plus_gcg= run_chi("katp_recess_diff", rx = octreotide(30), end = 48, closed_loop = TRUE,
                        extra = list(GCG_inf = 15))
)

## =====================================================================
##  SCENARIO 7 — RESECTION EXTENT.  Surgery changes B, never g, so ONE
##  equation produces both failure modes.  Verified mean glucose:
##    remnant 0.50 -> still needs 3.4 mg/kg/min dextrose
##    remnant 0.30 -> 78 mg/dL   (euglycaemic)
##    remnant 0.20 -> 97 mg/dL   (euglycaemic)
##    remnant 0.10 -> 177 mg/dL  (DIABETIC)
##    remnant 0.02 -> 342 mg/dL  (DIABETIC)
##  Growth then re-reads the same remnant: a 2 % remnant that is
##  adequate at 3.5 kg is 0.35 % of normal per kg at 20 kg, which is why
##  ~half of near-totally resected children become diabetic by
##  adolescence without anything new happening to the pancreas.
## =====================================================================
s7 <- lapply(c(1.0, 0.5, 0.3, 0.2, 0.1, 0.05, 0.02), function(b) {
  r <- run_chi("katp_recess_diff", end = 48, closed_loop = TRUE,
               extra = list(BMASS0 = b))
  filter(r, time > 36) %>%
    summarise(remnant = b, GIR = mean(GIRout), glucose = mean(GLU),
              insulin = mean(INS))
}) %>% bind_rows()

## =====================================================================
##  SCENARIO 8 — LEUCINE LOAD: the GDH-HI signature
##  Verified glucose fall after leucine 150 mg/kg PO:
##    normal 20 mg/dL · GDH-HI 37 mg/dL · dominant K_ATP 16 mg/dL
##  and ammonia 150 umol/L in GDH-HI, EXACTLY unchanged by diazoxide.
## =====================================================================
s8 <- lapply(c("normal","gdh_hi","katp_dominant"), function(gt) {
  run_chi(gt, rx = leucine(1.10), end = 6) %>% mutate(genotype = gt)
}) %>% bind_rows()

## =====================================================================
##  SCENARIO 9 — CEREBRAL FUEL MARGIN.  Withdraw feeds from a normal
##  neonate and from a CHI neonate and watch the two fuel curves.  The
##  normal baby's fuel ratio RISES on fasting (ketones); the CHI baby's
##  falls, because the same insulin that takes the glucose has already
##  taken the ketones.
## =====================================================================
s9 <- bind_rows(
  mod %>% param(geno$normal) %>% mrgsim(end = 24, delta = 0.1) %>%
    as_tibble() %>% mutate(arm = "normal, fasting"),
  mod %>% param(geno$katp_recess_diff) %>% mrgsim(end = 24, delta = 0.1) %>%
    as_tibble() %>% mutate(arm = "CHI g=0.02, fasting")
)

## =====================================================================
##  SCENARIO 10 — FOCAL vs DIFFUSE.  Severity depends only on the
##  PRODUCT (lesion fraction x local density): 0.10x10 reproduces diffuse
##  disease (11.6 vs 11.2 mg/kg/min) while 0.10x5 and 0.05x10 are
##  indistinguishable from each other (4.23 vs 4.29).  Lesionectomy is
##  modelled by setting w_ab = 0, which is curative; near-total
##  pancreatectomy is modelled by BMASS0, which is not.
## =====================================================================
s10 <- lapply(list(c(0.02,1), c(0.05,3), c(0.10,5), c(0.05,10), c(0.10,10), c(1.0,1)),
  function(wd) {
    r <- run_chi("katp_focal", end = 48, closed_loop = TRUE,
                 extra = list(w_ab = wd[1], dens_a = wd[2]))
    filter(r, time > 36) %>%
      summarise(w_ab = wd[1], density = wd[2], product = wd[1]*wd[2],
                GIR = mean(GIRout), insulin = mean(INS))
  }) %>% bind_rows()

## ---------------------------------------------------------------------
##  A quick look at the central claim
## ---------------------------------------------------------------------
if (interactive()) {
  print(s2); print(s3); print(s4); print(s7); print(s10)
  ggplot(s3) +
    geom_line(aes(g, dzx_pct, colour = "diazoxide"), linewidth = 1.2) +
    geom_line(aes(g, oct_pct, colour = "octreotide"), linewidth = 1.2) +
    scale_x_reverse() +
    labs(x = "residual K_ATP conductance  g", y = "reduction in dextrose need (%)",
         colour = NULL,
         title = "One equation, two verdicts",
         subtitle = "diazoxide multiplies g; octreotide adds to the divider") +
    theme_minimal(base_size = 12)
}
