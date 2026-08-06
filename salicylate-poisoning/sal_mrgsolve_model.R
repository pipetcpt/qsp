## =====================================================================
##  sal_mrgsolve_model.R
##  Salicylate (Aspirin) Poisoning — QSP / PK-PD model
##  28 ODE compartments · 14 therapeutic scenarios
##
##  살리실산(아스피린) 중독 — 정량적 시스템 약리학 모델
## ---------------------------------------------------------------------
##  STRUCTURAL THESIS
##  -----------------
##  Salicylate poisoning is graded, triaged and dialysed on the strength of
##  ONE number — the TOTAL plasma salicylate concentration.  That number is
##  separated from the quantity which actually does the killing (salicylate
##  inside brain cells) by TWO multiplications, neither of which is measured:
##
##      C_brain = C_total x fu(C_total, pH) x f_n(pH_plasma)/f_n(pH_brain)
##                         \______________/  \___________________________/
##                          MULTIPLIER 1      MULTIPLIER 2
##                          albumin saturates the pH partition of a weak
##                          fu: 8.5% at 150   acid, pKa 3.0
##                          mg/L -> 43% at
##                          800 mg/L
##
##  MULTIPLIER 2 is where the clinical paradoxes live, and it turns on a
##  piece of physiology that is usually left out of pharmacokinetic models:
##  BRAIN INTRACELLULAR pH IS BUFFERED AND PLASMA pH IS NOT.  Carbon dioxide
##  crosses cell membranes instantly and is then absorbed by non-bicarbonate
##  intracellular buffers (beta ~ 35 mmol/L per pH unit), so raising PaCO2
##  from 25 to 60 mmHg moves plasma pH by 0.38 units but brain intracellular
##  pH by only 0.14.  The GRADIENT therefore swings, and with it the
##  partition coefficient — by a factor of 1.7 in this model.
##
##  THREE CONSEQUENCES THE MODEL DERIVES RATHER THAN ASSUMES
##  --------------------------------------------------------
##  (1) THE VENTILATOR IS A DRUG WITH A NARROW THERAPEUTIC INDEX.
##      Intubating the hyperventilating patient and setting a conventional
##      minute ventilation raises PaCO2 25 -> 57, drops arterial pH 7.47 ->
##      7.11, and RAISES brain salicylate — while the MEASURED plasma level
##      FALLS by 28% because the same acidaemia drives drug into muscle as
##      well.  The laboratory report improves as the patient is poisoned.
##      Scenario S07 vs S08 is that experiment with one switch changed.
##
##  (2) BICARBONATE IS THE ONLY AGENT THAT ACTS ON BOTH MULTIPLIERS AT ONCE.
##      Raising plasma pH lowers the brain partition within minutes (drug
##      leaves the brain), and raising urine pH removes the tubular
##      back-diffusion that returns 97% of the filtered load at urine pH 5.
##      Renal clearance rises 11-fold; the plasma-pH arm acts hours earlier
##      than the urine arm.  Neither arm is available to acetazolamide,
##      which alkalinises the urine while ACIDIFYING the blood.
##
##  (3) THE THERAPY IS RATE-LIMITED BY POTASSIUM.  The urine cannot be
##      alkalinised until plasma bicarbonate exceeds the renal bicarbonate
##      threshold, and hypokalaemia RAISES that threshold (5.5 mmol/L per
##      mmol/L of K deficit).  Withholding potassium delays urine pH 7.5 by
##      about 10 h and requires 3 mmol/L more plasma bicarbonate to get
##      there — and leaves brain salicylate 3-fold higher at 48 h.
##
##  WHAT IS FITTED AND WHAT IS PREDICTED
##  ------------------------------------
##  FITTED (to normal or therapeutic human data only):
##      - two-class albumin binding (fu 8.5% at 150, 29% at 500, 43% at 800)
##      - Levy capacity-limited metabolism (Vmax/Km for salicylurate and
##        phenolic glucuronide), giving t1/2 3.6 h after 650 mg and a
##        steady state of 200 mg/L on 3.9 g/day
##      - the tubular back-diffusion constant, set so that 1-6% of a
##        therapeutic dose appears unchanged in acid urine
##      - intracellular buffer powers, set to the 31P-MRS observation that
##        brain pHi is defended against acute PaCO2 change
##  PREDICTED (nothing was tuned to these):
##      - overdose t1/2 37 h, apparent Vd rising 0.11 -> 0.32 L/kg
##      - respiratory alkalosis first, mixed disorder second
##      - the intubation catastrophe and its falling plasma level
##      - the 10-hour potassium delay
##      - neuroglycopenia at a normal serum glucose
##      - the misclassification pattern of the Done nomogram
##
##  VERIFICATION
##  ------------
##  All 28 ODEs were independently re-implemented in Python/scipy
##  (sal_verify_python.py) and integrated against 54 published anchors;
##  the report is sal_verification_output.txt (54/54 pass).  That exercise
##  found and fixed four defects in earlier drafts of THIS file:
##    (a) brain intracellular pH computed from a FIXED intracellular
##        bicarbonate, which made acute hypercapnia change brain pH as much
##        as plasma pH and abolished the entire partition effect;
##    (b) single-class albumin binding, which saturated so completely that
##        the free fraction reached 64% at 800 mg/L and acidaemia could no
##        longer displace anything;
##    (c) no maintenance water intake, so every long simulation ended in a
##        spurious contraction acidosis with GFR at the floor;
##    (d) an untitrated bicarbonate infusion, which drove plasma HCO3 to 45
##        mmol/L and hid the potassium failure mode completely.
##
##  A KNOWN LIMITATION, STATED RATHER THAN HIDDEN
##  ---------------------------------------------
##  The model asserts the brain concentration that no clinical study has
##  ever measured in a living patient.  Its brain compartment is calibrated
##  to animal partition data (Hill 1971) and to the pH physics, not to human
##  outcomes.  The falsifiable prediction is the DIRECTION AND SIZE of the
##  plasma/brain divergence after a PaCO2 change: plasma -28%, brain +12%
##  within one hour, brain +21% by 24 h.
##
##  DISCLAIMER: educational and research model.  Not for clinical use.
## =====================================================================

library(mrgsolve)
library(dplyr)

sal_code <- '
$PROB
# Salicylate poisoning QSP model (28 ODEs)
# C_brain = C_total x fu x f_n(pH_plasma)/f_n(pH_brain)

$PARAM @annotated
BW      : 70    : body weight (kg)

// ---------------- absorption ----------------
KAST    : 0.55  : gastric emptying rate constant (1/h)
KAGUT   : 1.60  : intestinal absorption rate constant (1/h)
KACONC  : 0.030 : concretion dissolution rate constant (1/h)
FCONC   : 0     : fraction of dose forming a concretion (-)
FBIO    : 0.90  : oral bioavailability (-)
KHYD    : 2.10  : systemic hydrolysis aspirin to salicylate (1/h)
KAC     : 6.0   : charcoal adsorption constant (1/h per g)
KACOUT  : 0.35  : charcoal gut transit (1/h)

// ---------------- distribution ----------------
VCEN    : 5.0   : central measured volume (L)
VPER    : 28.0  : peripheral tissue water (L)
VCNS    : 1.20  : brain water (L)
QPER    : 45.0  : central-peripheral exchange of free drug (L/h)
PSCNS   : 0.95  : blood-brain permeability-surface product (L/h)
KPT     : 1.50  : tissue trapping factor at reference pH pair (-)

// ---------------- plasma protein binding ----------------
BSITE1  : 232.5 : high-affinity albumin capacity (mg/L)
KDB1    : 11.67 : high-affinity dissociation constant (mg/L)
BSITE2  : 491.4 : low-affinity albumin capacity (mg/L)
KDB2    : 380.8 : low-affinity dissociation constant (mg/L)
ALB     : 40    : serum albumin (g/L)
KDPH    : 2.00  : acidaemic weakening of binding (per pH unit)

// ---------------- ionisation ----------------
PKA     : 3.00  : salicylic acid pKa (-)

// ---------------- hepatic metabolism ----------------
VMAXSU  : 100   : Vmax salicylurate formation (mg/h)
KMSU    : 34    : Km salicylurate formation (mg/L)
VMAXPG  : 27    : Vmax phenolic glucuronide (mg/h)
KMPG    : 550   : Km phenolic glucuronide (mg/L)
CLAG    : 0.045 : acyl glucuronide clearance (L/h)
CLGA    : 0.025 : gentisic acid clearance (L/h)
GLYT    : 9.0   : glycine pool turnover time (h)
GLYUSE  : 0.0030: glycine drain per 100 mg/h salicylurate (-)

// ---------------- renal ----------------
GFR0    : 7.2   : normal GFR (L/h)
AREAB   : 260   : tubular back-diffusion constant (L/h)
HCO3TH  : 25.0  : renal bicarbonate threshold (mmol/L)
THK     : 0.22  : threshold rise per mmol/L hypokalaemia (-)
THVOL   : 0.55  : threshold rise per fractional ECF contraction (-)
THCO2   : 0.20  : threshold shift per mmHg PaCO2 (mmol/L per mmHg)
UOB     : 0.050 : basal urine output (L/h)
UOG     : 1.20  : urine output per L of ECF expansion (L/h per L)
UKB     : 20    : basal urinary potassium (mmol/L)
UKHCO3  : 0.35  : urinary K per urinary bicarbonate (-)
UPHOFF  : 1.35  : urinary non-bicarbonate buffer floor (mmol/L)

// ---------------- acid-base ----------------
VBIC    : 35.0  : bicarbonate space (L)
TAUCO2  : 0.12  : lung time constant (h)
EMAXR   : 1.65  : maximal fractional rise in alveolar ventilation (-)
EC50R   : 22.0  : brain salicylate for half-maximal drive (mg/L)
KCHEM   : 7.0   : ventilatory drive per unit fall in brain pH (-)
KFON    : 0.030 : respiratory muscle fatigue accrual (1/h)
KFOFF   : 0.10  : fatigue recovery (1/h)
NAG0    : 2.9   : dietary fixed acid load (mmol/h)
NAGC    : 20.0  : renal acid excretion per pH unit below 7.4 (mmol/h)
PHCAP   : 7.60  : arterial pH ceiling for the bicarbonate infusion (-)
PHCAPW  : 0.015 : width of the titration switch (-)

// ---------------- intracellular pH ----------------
BBRATIO : 0.536 : brain buffer base / plasma bicarbonate (-)
TAUBB   : 5.8   : brain equilibration time constant (h)
BETAB   : 35.0  : brain non-bicarbonate buffer power (mmol/L per pH)
PHREFB  : 7.05  : reference brain intracellular pH (-)
DPCO2B  : 8.0   : brain tissue PCO2 above arterial (mmHg)
TBRATIO : 0.458 : muscle buffer base / plasma bicarbonate (-)
TAUTB   : 4.3   : muscle equilibration time constant (h)
BETAT   : 30.0  : muscle non-bicarbonate buffer power (mmol/L per pH)
PHREFT  : 7.00  : reference muscle intracellular pH (-)
DPCO2T  : 6.0   : muscle tissue PCO2 above arterial (mmHg)

// ---------------- uncoupling ----------------
EC50U   : 200   : tissue salicylate for half-maximal uncoupling (mg/L)
KLAC    : 26.0  : lactate production at full uncoupling (mmol/h)
CLLAC   : 0.35  : lactate clearance (1/h)
KKET    : 12.0  : ketoacid production at full uncoupling (mmol/h)
CLKET   : 0.15  : ketoacid clearance (1/h)

// ---------------- potassium and volume ----------------
VKAP    : 250   : apparent potassium distribution volume (L)
KDIET   : 3.6   : dietary potassium intake (mmol/h)
KPHSH   : 3.0   : transcellular K shift per pH unit (mmol/L)
VECF0   : 14.0  : normal ECF volume (L)
RINTK   : 0.110 : oral maintenance water intake (L/h)
INSENS  : 0.045 : insensible water loss (L/h)
SWEAT   : 0.16  : sweat at full uncoupling (L/h)
GFRVOL  : 3.0   : GFR sensitivity to ECF contraction (-)

// ---------------- brain glucose ----------------
GLCP    : 5.5   : plasma glucose (mmol/L)
TMAXG   : 60.0  : GLUT1 transport capacity (mmol/L/h)
KTG     : 6.0   : GLUT1 Km (mmol/L)
CMR0    : 18.0  : cerebral glucose consumption (mmol/L/h)
CMRU    : 0.50  : extra consumption at full uncoupling (-)

// ---------------- thermogenesis and lung ----------------
KTEMP   : 0.90  : heat production at full uncoupling (degC/h)
CLTEMP  : 0.50  : heat dissipation (1/h)
KLUNG   : 0.055 : lung injury accrual (1/h)
CLLUNG  : 0.030 : lung injury resolution (1/h)
EC50LU  : 150   : brain salicylate for half-maximal lung injury (mg/L)

// ---------------- toxicodynamics ----------------
EC50CNS : 95.0  : brain salicylate for half-maximal CNS depression (mg/L)
HILLCNS : 2.2   : Hill coefficient for CNS depression (-)
CNSTHR  : 45.0  : brain salicylate threshold for cumulative injury (mg/L)
GLCTHR  : 0.80  : brain glucose below which injury accelerates (mmol/L)

// ---------------- interventions ----------------
RBIC    : 0     : sodium bicarbonate infusion (mmol/h)
RFLU    : 0     : crystalloid infusion (L/h)
RKCL    : 0     : potassium replacement (mmol/h)
CLHD    : 0     : extracorporeal clearance (L/h)
HDBIC   : 0     : bicarbonate gained from dialysate (mmol/h)
VENT    : 0     : 0 spontaneous / 1 controlled ventilation (-)
VASET   : 1.0   : relative alveolar ventilation when ventilated (-)
SED     : 0     : fractional respiratory depression (-)
AGEF    : 0     : elderly flag for lung injury susceptibility (-)
HCO30   : 24.0  : baseline plasma bicarbonate (mmol/L)
PACO20  : 40.0  : baseline arterial PCO2 (mmHg)

$CMT @annotated
AST   : aspirin in the stomach (mg)
AGUT  : aspirin in the small intestine (mg)
ACONC : aspirin in a concretion (mg)
AASA  : unhydrolysed aspirin, systemic (mg)
ACENT : salicylate, central compartment (mg)
APER  : salicylate, peripheral tissue (mg)
ACNS  : salicylate, brain (mg)
ASU   : salicylurate formed (mg salicylate-equivalent)
APG   : phenolic glucuronide formed (mg)
AAG   : acyl glucuronide formed (mg)
AGA   : gentisate formed (mg)
AUR   : salicylate excreted in urine (mg)
AHD   : salicylate removed by dialysis (mg)
AAC   : activated charcoal in the gut (g)

$INIT @annotated
HCO3  : 24.0  : plasma bicarbonate (mmol/L)
PACO2 : 40.0  : arterial PCO2 (mmHg)
BBB   : 12.86 : brain buffer base (mmol/L)
BBT   : 10.99 : muscle buffer base (mmol/L)
LAC   : 1.0   : lactate (mmol/L)
KET   : 0.10  : ketoacid (mmol/L)
KBAL  : 4.0   : potassium balance level (mmol/L)
VECF  : 14.0  : extracellular fluid volume (L)
GLY   : 1.0   : glycine pool (fraction of normal)
TCORE : 37.0  : core temperature (degC)
GLUB  : 1.30  : brain glucose (mmol/L)
LUNG  : 0     : lung injury (0-1)
CNSI  : 0     : cumulative CNS injury (units)
FATIG : 0     : respiratory muscle fatigue (0-1)

$GLOBAL
#define MWSAL 138.12
#define MWASA 180.16

// ---- unbound salicylate from total plasma concentration -------------
// bound(Cf) = B1*Cf/(K1+Cf) + B2*Cf/(K2+Cf).  Newton iteration; the same
// routine is used in sal_verify_python.py.
double freeSal(double ctot, double ph, double b1, double k1, double b2,
               double k2, double alb, double kdph) {
  if (ctot <= 1e-12) return 0.0;
  double sc = alb / 40.0;
  double B1 = b1 * sc, B2 = b2 * sc;
  double dph = 7.40 - ph; if (dph < 0) dph = 0;
  double m  = 1.0 + kdph * dph;
  double K1 = k1 * m, K2 = k2 * m;
  double cf = 0.3 * ctot + 1e-6;
  for (int i = 0; i < 60; i++) {
    double f = cf + B1*cf/(K1+cf) + B2*cf/(K2+cf) - ctot;
    double d = 1.0 + B1*K1/((K1+cf)*(K1+cf)) + B2*K2/((K2+cf)*(K2+cf));
    double nn = cf - f/d;
    if (nn <= 0.0) nn = 0.5*cf;
    if (fabs(nn - cf) < 1e-10*(1.0 > cf ? 1.0 : cf)) { cf = nn; break; }
    cf = nn;
  }
  if (cf < 0) cf = 0; if (cf > ctot) cf = ctot;
  return cf;
}

// ---- un-ionised (membrane-permeant) fraction of a weak acid ----------
double fNeutral(double ph, double pka) {
  return 1.0 / (1.0 + pow(10.0, ph - pka));
}

// ---- plasma pH ------------------------------------------------------
double phPlasma(double hco3, double pco2) {
  double h = hco3 > 0.30 ? hco3 : 0.30;
  double p = pco2 > 5.0  ? pco2 : 5.0;
  return 6.10 + log10(h / (0.0301 * p));
}

// ---- intracellular pH with a non-bicarbonate buffer power ------------
// [HCO3-]i(pH) = BUF - beta*(pH - phref);  pH = 6.1 + log10(HCO3i/(0.0301*PCO2))
// Newton on pH.  This is the step that makes acute hypercapnia move plasma
// pH far more than it moves intracellular pH.
double phIntra(double buf, double pco2, double beta, double phref) {
  double x = phref;
  double p = pco2 > 5.0 ? pco2 : 5.0;
  for (int i = 0; i < 40; i++) {
    double h = buf - beta * (x - phref);
    if (h < 0.20) h = 0.20;
    double g  = 6.10 + log10(h / (0.0301 * p));
    double fx = g - x;
    double dfx = -0.434294 * beta / h - 1.0;
    double xn = x - fx/dfx;
    if (xn < 5.0) xn = 5.0;
    if (xn > 8.5) xn = 8.5;
    if (fabs(xn - x) < 1e-10) { x = xn; break; }
    x = xn;
  }
  return x;
}

// ---- urine pH from urinary bicarbonate concentration ------------------
double phUrine(double uhco3, double off) {
  double v = 6.10 + log10((uhco3 + off) / 1.354);
  if (v > 8.20) v = 8.20;
  if (v < 4.60) v = 4.60;
  return v;
}

$MAIN
HCO3_0  = HCO30;
PACO2_0 = PACO20;
BBB_0   = BBRATIO * HCO30;
BBT_0   = TBRATIO * HCO30;
VECF_0  = VECF0;

$ODE
// ================= algebra =================
double ctot  = (ACENT > 0 ? ACENT : 0) / VCEN;
double hco3  = HCO3  > 1.0 ? HCO3  : 1.0;
double paco2 = PACO2 > 5.0 ? PACO2 : 5.0;
double ph    = phPlasma(hco3, paco2);
double cf    = freeSal(ctot, ph, BSITE1, KDB1, BSITE2, KDB2, ALB, KDPH);
double fu    = ctot > 1e-9 ? cf/ctot : 0.085;

double phb = phIntra(BBB > 1.0 ? BBB : 1.0, paco2 + DPCO2B, BETAB, PHREFB);
double pht = phIntra(BBT > 1.0 ? BBT : 1.0, paco2 + DPCO2T, BETAT, PHREFT);
double ccns = (ACNS > 0 ? ACNS : 0) / VCNS;
double ctis = (APER > 0 ? APER : 0) / VPER;

double kpl = KBAL - KPHSH * (ph - 7.40); if (kpl < 0.8) kpl = 0.8;

double vecf = VECF > 6.0 ? VECF : 6.0;
double contract = (VECF0 - vecf) / VECF0; if (contract < 0) contract = 0;
double gfrf = 1.0 - GFRVOL * contract;
if (gfrf > 1.10) gfrf = 1.10; if (gfrf < 0.12) gfrf = 0.12;
double gfr = GFR0 * gfrf;

double uo = (UOB + UOG * (vecf - VECF0)) * gfrf;
if (uo > 0.80) uo = 0.80; if (uo < 0.004) uo = 0.004;

double kdef = 4.0 - kpl; if (kdef < 0) kdef = 0;
double thr = HCO3TH * (1.0 + THK * kdef + THVOL * contract)
             + THCO2 * (paco2 - 40.0);
double juh = gfr * (hco3 - thr); if (juh < 0) juh = 0;
double uhco3 = juh / uo;
double uph = phUrine(uhco3, UPHOFF);

double fnu = fNeutral(uph, PKA);
double freab = 1.0 - exp(-AREAB * fnu / uo);
double clren = gfr * fu * (1.0 - freab);

double unc = ctis / (EC50U + ctis);

// ---- the two partition coefficients ----
double fnp = fNeutral(ph,  PKA);
double fnb = fNeutral(phb, PKA);
double fnt = fNeutral(pht, PKA);
double kpbrain  = fnp / fnb;                       // brain water : free plasma
double kptissue = fnp / fnt;

// ---- ventilation ----
double drives = EMAXR * ccns / (EC50R + ccns);
double drivea = KCHEM * (PHREFB - phb);
if (drivea < 0) drivea = 0; if (drivea > 1.5) drivea = 1.5;
double vaspont = (1.0 + drives + drivea) * (1.0 - 0.85*FATIG)
                 * (1.0 - 0.45*LUNG) * (1.0 - SED);
double va = VENT > 0.5 ? VASET : vaspont;
if (va < 0.15) va = 0.15;
double paco2ss = 40.0 * (1.0 + 0.60*unc) / va;
if (paco2ss > 120.0) paco2ss = 120.0; if (paco2ss < 8.0) paco2ss = 8.0;

// ---- elimination fluxes ----
double gly = GLY > 0.05 ? GLY : 0.05;
double vsu = VMAXSU * gly * ctot / (KMSU + ctot);
double vpg = VMAXPG * ctot / (KMPG + ctot);
double vag = CLAG * ctot;
double vga = CLGA * ctot;
double vre = clren * ctot;
double vhd = CLHD * ctot;

// ================= gut =================
double ac = AAC > 0 ? AAC : 0;
double kbind = KAC * ac / (1.0 + KAC * ac / 3.0);
dxdt_AST   = -KAST*AST - kbind*AST;
dxdt_ACONC = -KACONC*ACONC;
dxdt_AGUT  = KAST*AST + KACONC*ACONC - KAGUT*AGUT - kbind*AGUT;
dxdt_AAC   = -KACOUT*ac;

double absorbed = FBIO * KAGUT * AGUT;              // mg/h as aspirin
dxdt_AASA = absorbed - KHYD*AASA;
double tosal = KHYD * AASA * MWSAL / MWASA;         // mg/h as salicylate

// ================= disposition =================
double kpt = KPT * kptissue / (fNeutral(7.40,PKA) / fNeutral(7.00,PKA));
if (kpt < 1e-3) kpt = 1e-3;
double jper = QPER * (cf - APER / (VPER * kpt));
double jcns = PSCNS * (cf - ccns / (kpbrain > 1e-3 ? kpbrain : 1e-3));

dxdt_APER  = jper;
dxdt_ACNS  = jcns;
dxdt_ACENT = tosal - jper - jcns - vsu - vpg - vag - vga - vre - vhd;
dxdt_ASU   = vsu;
dxdt_APG   = vpg;
dxdt_AAG   = vag;
dxdt_AGA   = vga;
dxdt_AUR   = vre;
dxdt_AHD   = vhd;
dxdt_GLY   = (1.0 - GLY)/GLYT - GLYUSE * vsu / 100.0;

// ================= organic acids and bicarbonate =================
double dlac = KLAC*unc/VBIC - CLLAC*(LAC - 1.0);
double dket = KKET*unc/VBIC - CLKET*(KET - 0.10);
dxdt_LAC = dlac;
dxdt_KET = dket;

double salacid = tosal / MWSAL;                     // mmol/h of acid
double dph74 = 7.40 - ph; if (dph74 < 0) dph74 = 0;
double renalcomp = NAGC * dph74 * gfrf;
// clinical titration: hold the bicarbonate infusion at the arterial pH ceiling
double bicgate = 1.0 / (1.0 + exp((ph - PHCAP)/PHCAPW));
dxdt_HCO3 = -dlac - dket
            + (RBIC*bicgate + HDBIC + renalcomp - salacid - juh
               + NAG0*(gfrf - 1.0)) / VBIC;

// ================= ventilation =================
dxdt_PACO2 = (paco2ss - PACO2)/TAUCO2;
double exc = vaspont - 2.4; if (exc < 0) exc = 0;
dxdt_FATIG = KFON*exc*(1.0 - FATIG) - KFOFF*FATIG;

// ================= intracellular buffer bases =================
dxdt_BBB = (BBRATIO*hco3 - BBB)/TAUBB;
dxdt_BBT = (TBRATIO*hco3 - BBT)/TAUTB;

// ================= potassium and volume =================
double uk = UKB + UKHCO3*uhco3; if (uk > 120.0) uk = 120.0;
dxdt_KBAL = (KDIET + RKCL - uo*uk)/VKAP;
dxdt_VECF = RINTK + RFLU - uo - INSENS - SWEAT*unc;

// ================= brain glucose =================
double gb = GLUB > 0 ? GLUB : 0;
dxdt_GLUB = TMAXG*(GLCP/(KTG+GLCP) - gb/(KTG+gb)) - CMR0*(1.0 + CMRU*unc);

// ================= temperature, lung, CNS injury =================
dxdt_TCORE = KTEMP*unc - CLTEMP*(TCORE - 37.0);
dxdt_LUNG  = KLUNG*(1.0+AGEF)*ccns/(EC50LU+ccns)*(1.0-LUNG) - CLLUNG*LUNG;
double over = ccns - CNSTHR; if (over < 0) over = 0;
double hypo = (GLCTHR - gb)/GLCTHR; if (hypo < 0) hypo = 0;
dxdt_CNSI = over/100.0 * (1.0 + 2.0*hypo);

$TABLE
double ctotT  = (ACENT > 0 ? ACENT : 0)/VCEN;
double hco3T  = HCO3  > 1.0 ? HCO3  : 1.0;
double paco2T = PACO2 > 5.0 ? PACO2 : 5.0;
capture PH    = phPlasma(hco3T, paco2T);
capture CTOT  = ctotT;
capture CTOTD = ctotT / 10.0;                       // mg/dL, as reported
capture CFREE = freeSal(ctotT, PH, BSITE1, KDB1, BSITE2, KDB2, ALB, KDPH);
capture FU    = ctotT > 1e-9 ? CFREE/ctotT : 0.085;
capture PHB   = phIntra(BBB > 1.0 ? BBB : 1.0, paco2T + DPCO2B, BETAB, PHREFB);
capture PHT   = phIntra(BBT > 1.0 ? BBT : 1.0, paco2T + DPCO2T, BETAT, PHREFT);
capture CCNS  = (ACNS > 0 ? ACNS : 0)/VCNS;
capture CTIS  = (APER > 0 ? APER : 0)/VPER;

// the two multipliers, reported explicitly
capture KPBR  = fNeutral(PH,PKA) / fNeutral(PHB,PKA);      // MULTIPLIER 2
capture RATIO = ctotT > 1e-9 ? CCNS/ctotT : 0;             // brain : reported level

capture KPL   = (KBAL - KPHSH*(PH - 7.40)) > 0.8 ? (KBAL - KPHSH*(PH-7.40)) : 0.8;
double contractT = (VECF0 - (VECF > 6 ? VECF : 6))/VECF0; if (contractT < 0) contractT = 0;
double gfrfT = 1.0 - GFRVOL*contractT;
if (gfrfT > 1.10) gfrfT = 1.10; if (gfrfT < 0.12) gfrfT = 0.12;
capture GFRC  = GFR0*gfrfT*1000.0/60.0;                    // mL/min
double uoT = (UOB + UOG*((VECF > 6 ? VECF : 6) - VECF0))*gfrfT;
if (uoT > 0.80) uoT = 0.80; if (uoT < 0.004) uoT = 0.004;
capture UO    = uoT*1000.0;                                // mL/h
double kdefT = 4.0 - KPL; if (kdefT < 0) kdefT = 0;
double thrT = HCO3TH*(1.0 + THK*kdefT + THVOL*contractT) + THCO2*(paco2T - 40.0);
double juhT = GFR0*gfrfT*(hco3T - thrT); if (juhT < 0) juhT = 0;
capture UPH   = phUrine(juhT/uoT, UPHOFF);
double freabT = 1.0 - exp(-AREAB*fNeutral(UPH,PKA)/uoT);
capture CLREN = GFR0*gfrfT*FU*(1.0-freabT)*1000.0/60.0;    // mL/min
capture HCO3C = hco3T;
capture PACO2C = paco2T;
capture UNC   = CTIS/(EC50U + CTIS);
capture AGAP  = 12.0 + (LAC - 1.0) + (KET - 0.10) + ctotT/MWSAL;
capture CNS   = 100.0*pow(CCNS,HILLCNS)/(pow(EC50CNS,HILLCNS)+pow(CCNS,HILLCNS));
capture TINN  = 100.0*CCNS*CCNS/(18.0*18.0 + CCNS*CCNS);
capture VDAPP = ctotT > 1e-6 ? (ACENT+APER+ACNS)/ctotT/BW : VCEN/BW;
double sev = 0.45*CNS/100.0
           + 0.25*((7.35-PH)/0.30 > 0 ? ((7.35-PH)/0.30 < 1 ? (7.35-PH)/0.30 : 1) : 0)
           + 0.15*((GLCTHR-GLUB)/GLCTHR > 0 ? ((GLCTHR-GLUB)/GLCTHR < 1 ? (GLCTHR-GLUB)/GLCTHR : 1) : 0)
           + 0.15*((TCORE-37.5)/2.5 > 0 ? ((TCORE-37.5)/2.5 < 1 ? (TCORE-37.5)/2.5 : 1) : 0);
capture SEV   = 100.0*(sev < 1 ? sev : 1);
capture VASP  = (1.0 + EMAXR*CCNS/(EC50R+CCNS)
                 + ((KCHEM*(PHREFB-PHB) > 0) ? (KCHEM*(PHREFB-PHB) < 1.5 ?
                     KCHEM*(PHREFB-PHB) : 1.5) : 0))
                * (1.0-0.85*FATIG) * (1.0-0.45*LUNG) * (1.0-SED);
'

mod <- mcode_cache("sal_qsp", sal_code, atol = 1e-9, rtol = 1e-7, maxsteps = 1e6)

WT <- 70   # reference adult

## =====================================================================
##  Dosing and intervention helpers
## =====================================================================
##  Every event is expressed as an explicit (time, amt, cmt) row.  Parameter
##  changes over time are handled by simulating stage by stage and carrying
##  the state vector forward, so the event list must be filtered per stage —
##  which is why `ii`/`addl` are expanded here rather than left to mrgsolve.

## Aspirin ingestion: `g` grams into the stomach compartment.
asa_dose <- function(g, time = 0, cmt = "AST") {
  data.frame(time = time, amt = g * 1000, cmt = cmt)
}

## Repeated oral dosing: n doses of `mg` every `ii` hours from `start`.
dose_seq <- function(mg, ii, n, start = 0, cmt = "AST") {
  data.frame(time = start + ii * (seq_len(n) - 1), amt = mg, cmt = cmt)
}

## Activated charcoal, grams into the gut charcoal compartment.
charcoal <- function(g = 50, time = 0) {
  data.frame(time = time, amt = g, cmt = "AAC")
}

## A sodium bicarbonate BOLUS is a step change in the HCO3 state
## (mmol divided by the bicarbonate space).
bic_bolus <- function(mmol = 100, time = 0, vbic = 35) {
  data.frame(time = time, amt = mmol / vbic, cmt = "HCO3")
}

## Build an mrgsolve event object from a (time, amt, cmt) data frame.
make_ev <- function(df) {
  if (is.null(df) || !nrow(df)) return(NULL)
  Reduce(`+`, lapply(seq_len(nrow(df)), function(i)
    mrgsolve::ev(time = df$time[i], amt = df$amt[i], cmt = as.character(df$cmt[i]))))
}

## Simulate a scenario as a sequence of parameter stages.
##   stages : list of list(t = start_time, par = list(...)); parameters are
##            cumulative, so each stage only names what changes.
##   events : data frame of (time, amt, cmt) covering the whole run.
run_stages <- function(model, stages, end = 72, delta = 0.05, events = NULL) {
  bounds <- c(vapply(stages, function(s) s$t, numeric(1)), end)
  state_names <- names(as.list(mrgsolve::init(model)))
  out <- NULL
  cur_par <- list()
  cur_init <- NULL
  for (i in seq_along(stages)) {
    cur_par <- utils::modifyList(cur_par, stages[[i]]$par)
    t0 <- bounds[i]; t1 <- bounds[i + 1]
    if (t1 <= t0) next
    m <- model %>% param(cur_par)
    if (!is.null(cur_init)) m <- m %>% init(cur_init)
    ee <- NULL
    if (!is.null(events)) {
      sel <- events[events$time >= t0 - 1e-9 & events$time < t1 - 1e-9, , drop = FALSE]
      if (nrow(sel)) { sel$time <- sel$time - t0; ee <- make_ev(sel) }
    }
    seg <- if (is.null(ee)) m %>% mrgsim(end = t1 - t0, delta = delta)
           else             m %>% mrgsim(events = ee, end = t1 - t0, delta = delta)
    d <- as.data.frame(seg)
    d$time <- d$time + t0
    cur_init <- as.list(d[nrow(d), state_names])
    out <- if (is.null(out)) d else dplyr::bind_rows(out, d[-1, ])
  }
  out
}

## =====================================================================
##  FOURTEEN SCENARIOS
## =====================================================================
##  All return the same columns, so they bind together and compare directly.
##  The comparisons that carry the argument of the model are:
##      S07 vs S08   the ventilator experiment (one switch changed)
##      S05 vs S06   potassium replaced or not
##      S04 vs S05   saline diuresis vs bicarbonate at equal volume
##      S09 vs S10   dialysis early or late
##      S02 vs S12   the same plasma level, acute or chronic
## =====================================================================

SUPPORTIVE <- list(RFLU = 0.125, RKCL = 4.0)
ALKALI_1   <- list(RBIC = 75.0, RFLU = 0.25, RKCL = 10.0)   # loading, 4 h
ALKALI_2   <- list(RBIC = 37.5, RFLU = 0.25, RKCL = 10.0)   # maintenance
MASSIVE    <- list(KAST = 0.40, KAGUT = 1.20)               # large tablet load
ENTERIC    <- list(KAST = 0.16, KAGUT = 0.40)

SCEN <- list(

  S01 = list(
    label  = "Therapeutic: aspirin 650 mg q4h x 5 days",
    end    = 132,
    stages = list(list(t = 0, par = list())),
    events = dose_seq(650, ii = 4, n = 30)),

  S02 = list(
    label  = "Acute 30 g, supportive care only (crystalloid)",
    end    = 72,
    stages = list(list(t = 0, par = MASSIVE), list(t = 4, par = SUPPORTIVE)),
    events = asa_dose(30)),

  S03 = list(
    label  = "Acute 30 g, no fluid resuscitation",
    end    = 72,
    stages = list(list(t = 0, par = MASSIVE)),
    events = asa_dose(30)),

  S04 = list(
    label  = "Acute 30 g + saline diuresis 250 mL/h",
    end    = 72,
    stages = list(list(t = 0, par = MASSIVE),
                  list(t = 4, par = list(RFLU = 0.25, RKCL = 10.0))),
    events = asa_dose(30)),

  S05 = list(
    label  = "Acute 30 g + NaHCO3 with potassium (standard care)",
    end    = 72,
    stages = list(list(t = 0, par = MASSIVE), list(t = 4, par = ALKALI_1),
                  list(t = 8, par = ALKALI_2)),
    events = rbind(asa_dose(30), bic_bolus(100, time = 4))),

  S06 = list(
    label  = "Acute 30 g + NaHCO3 WITHOUT potassium",
    end    = 72,
    stages = list(list(t = 0, par = MASSIVE),
                  list(t = 4, par = utils::modifyList(ALKALI_1, list(RKCL = 0))),
                  list(t = 8, par = utils::modifyList(ALKALI_2, list(RKCL = 0)))),
    events = rbind(asa_dose(30), bic_bolus(100, time = 4))),

  S07 = list(
    label  = "Acute 30 g, intubated at 8 h to a conventional PaCO2",
    end    = 36,
    stages = list(list(t = 0, par = MASSIVE), list(t = 4, par = SUPPORTIVE),
                  list(t = 8, par = list(VENT = 1, VASET = 1.0))),
    events = asa_dose(30)),

  S08 = list(
    label  = "Acute 30 g, intubated at 8 h with MATCHED ventilation",
    end    = 36,
    stages = list(list(t = 0, par = MASSIVE), list(t = 4, par = SUPPORTIVE),
                  list(t = 8, par = list(VENT = 1, VASET = 2.30))),
    events = asa_dose(30)),

  S09 = list(
    label  = "Acute 30 g + NaHCO3, haemodialysis at 6 h",
    end    = 48,
    stages = list(list(t = 0, par = MASSIVE), list(t = 4, par = ALKALI_1),
                  list(t = 6, par = c(ALKALI_1, list(CLHD = 6.0, HDBIC = 45))),
                  list(t = 10, par = c(ALKALI_2, list(CLHD = 0, HDBIC = 0)))),
    events = rbind(asa_dose(30), bic_bolus(100, time = 4))),

  S10 = list(
    label  = "Acute 30 g + NaHCO3, haemodialysis delayed to 18 h",
    end    = 48,
    stages = list(list(t = 0, par = MASSIVE), list(t = 4, par = ALKALI_1),
                  list(t = 8, par = ALKALI_2),
                  list(t = 18, par = c(ALKALI_2, list(CLHD = 6.0, HDBIC = 45))),
                  list(t = 22, par = c(ALKALI_2, list(CLHD = 0, HDBIC = 0)))),
    events = rbind(asa_dose(30), bic_bolus(100, time = 4))),

  S11 = list(
    label  = "Enteric-coated 30 g with a 45% concretion, MDAC",
    end    = 96,
    stages = list(list(t = 0, par = ENTERIC), list(t = 4, par = SUPPORTIVE)),
    events = rbind(asa_dose(16.5), asa_dose(13.5, cmt = "ACONC"),
                   charcoal(50, 2), charcoal(50, 6), charcoal(50, 10))),

  S12 = list(
    label  = "Chronic salicylism: 1.3 g q6h, 80 y, GFR 60, albumin 28",
    end    = 340,
    stages = list(list(t = 0, par = list(GFR0 = 3.6, ALB = 28, AGEF = 1)),
                  list(t = 168, par = list(RINTK = 0.055))),
    events = dose_seq(1300, ii = 6, n = 56)),

  S13 = list(
    label  = "Acute 30 g with opioid co-ingestion (45% depression of drive)",
    end    = 36,
    stages = list(list(t = 0, par = MASSIVE), list(t = 4, par = SUPPORTIVE),
                  list(t = 8, par = list(SED = 0.45))),
    events = asa_dose(30)),

  S14 = list(
    label  = "Acute 30 g + NaHCO3 + dextrose (plasma glucose 11 mmol/L)",
    end    = 72,
    stages = list(list(t = 0, par = MASSIVE),
                  list(t = 4, par = c(ALKALI_1, list(GLCP = 11.0))),
                  list(t = 8, par = c(ALKALI_2, list(GLCP = 11.0)))),
    events = rbind(asa_dose(30), bic_bolus(100, time = 4)))
)

run_scenario <- function(key, delta = 0.05) {
  s <- SCEN[[key]]
  d <- run_stages(mod, s$stages, end = s$end, delta = delta, events = s$events)
  d$scenario <- key
  d$label <- s$label
  d
}

run_all <- function() dplyr::bind_rows(lapply(names(SCEN), run_scenario))

## =====================================================================
##  THE FOUR COMPARISONS THAT CARRY THE ARGUMENT
## =====================================================================
##  (Values below are those produced by the independent Python
##   implementation; see sal_verification_output.txt.)
##
##  1. THE VENTILATOR EXPERIMENT               S07 (PaCO2 normalised) vs
##                                             S08 (ventilation matched)
##       t = 9 h   plasma  582 vs 812 mg/L     (-28%, the report improves)
##                 brain   162 vs 145 mg/L     (+12%, the patient worsens)
##                 pH     7.11 vs 7.47
##                 brain:plasma ratio  0.279 vs 0.178   (x1.57)
##       t = 24 h  brain   129 vs 107 mg/L     (+21%, because renal
##                 clearance has also collapsed into acid urine)
##       cumulative CNS injury at 36 h:  x1.31
##
##  2. POTASSIUM                                S05 (with K) vs S06 (without)
##       hours to urine pH 7.5      11.7  vs  21.9
##       plasma HCO3 needed to get there  +2.9 mmol/L more without K
##       renal clearance at 12 h    38.5  vs  17.9 mL/min
##       brain salicylate at 48 h    x2.97 higher without K
##
##  3. WHICH FLUID                              S04 (saline) vs S05 (NaHCO3)
##       renal clearance at 12 h    16.6  vs  38.5 mL/min   (x2.3)
##       versus supportive care alone            3.5 mL/min (x11)
##       and NaHCO3 lowers the BRAIN level within 40 min, before the
##       urine has responded at all - that is the plasma-pH arm.
##
##  4. DIALYSIS TIMING                          S09 (6 h) vs S10 (18 h)
##       cumulative CNS injury at 48 h   4.41  vs  12.55   (x2.85)
##       extracorporeal clearance 6 L/h is 5.7x the best achievable
##       alkaline-urine clearance, which is the pharmacological content
##       of the EXTRIP recommendation.
##
##  5. THE SAME NUMBER, TWO DIFFERENT PATIENTS  S02 (acute) vs S12 (chronic)
##       chronic day 14: plasma 216 mg/L (21.6 mg/dL - a "moderate" level)
##                       brain 109 mg/L, CNS score 58/100
##       acute patient first reaching 216 mg/L (t = 0.9 h):
##                       brain 3 mg/L, CNS score 0
##       brain:plasma ratio differs by a factor of 36.
## =====================================================================

## Quick numerical summary of any scenario.
summarise_scenario <- function(d) {
  data.frame(
    scenario = d$scenario[1],
    peak_plasma_mgL = max(d$CTOT),
    peak_brain_mgL  = max(d$CCNS),
    min_pH          = min(d$PH),
    max_PaCO2       = max(d$PACO2C),
    min_HCO3        = min(d$HCO3C),
    max_urine_pH    = max(d$UPH),
    max_CLren_mLmin = max(d$CLREN),
    min_K           = min(d$KPL),
    max_temp        = max(d$TCORE),
    min_brain_gluc  = min(d$GLUB),
    cum_CNS_injury  = max(d$CNSI),
    peak_CNS_score  = max(d$CNS),
    peak_severity   = max(d$SEV)
  )
}

if (identical(environment(), globalenv()) && !interactive()) {
  res <- run_all()
  print(dplyr::bind_rows(lapply(split(res, res$scenario), summarise_scenario)))
}
