## ===========================================================================
##  tap_mrgsolve_model.R
##  Toxic alcohol poisoning — methanol and ethylene glycol
##  Quantitative Systems Pharmacology model (mrgsolve)
##
##  독성 알코올 중독 (메탄올 · 에틸렌글리콜) QSP 모델
##
##  ---------------------------------------------------------------------------
##  ORGANISING THESIS
##  ---------------------------------------------------------------------------
##  The parent alcohol is not the poison.  Methanol and ethylene glycol are
##  ethanol-like sedatives with a large osmotic footprint and almost no organ
##  toxicity of their own; every lesion in the disease belongs to an ACID
##  METABOLITE that the patient's own liver manufactures.  The disease is
##  therefore not a concentration but the INTEGRAL OF A FLUX, and that flux is
##  one Michaelis-Menten expression in which every therapy owns exactly one
##  term:
##
##                        Vmax_i * (S_i/Km_i)
##      v_i  =  ------------------------------------------      FORMATION
##               1 + SUM_j (S_j/Km_j) + [FOM]/Ki                fomepizole and
##                                                             ethanol own the
##                                                             DENOMINATOR
##
##           -  CLren * (1 - f_reabs(pH_urine)) * ACID          RENAL
##                                                             bicarbonate owns
##                                                             f_reabs
##
##           -  CLhd * ACID                                     DIALYSIS
##                                                             the machine owns
##                                                             this term alone
##
##           -  Vmax_THF * THF * ACID/(Km + ACID)               FOLATE (formate)
##                                                             folinic acid owns
##                                                             Vmax_THF
##
##  and the damage is the CNS integral, gated by a term nobody prescribes:
##
##      d(ACID_cns)/dt = PS * ( f_HA(pH_plasma)*C_plasma
##                              - f_HA(pH_brain)*C_cns )
##            f_HA(pH) = 1/(1 + 10^(pH - pKa))    pKa formic 3.75, glycolic 3.83
##
##  ---------------------------------------------------------------------------
##  WHAT THIS FILE CONTAINS
##  ---------------------------------------------------------------------------
##    · 50 ODEs (methanol arm, ethylene glycol arm, antidotes, acid-base,
##      calcium/oxalate, organ injury, hazards, mass-balance bookkeeping)
##    · 25 therapy scenarios
##    · every derived read-out a clinician actually looks at, including the
##      osmolal-gap / anion-gap clock, the ADH inhibition factor, the CaOx
##      supersaturation product, the "lactate gap" assay artefact, logMAR
##      visual acuity, and hazard-based P(death) / P(blindness)
##    · a self-contained acceptance test (`tap_verify()`) against the 16-row
##      reference trajectory produced by tap_python_reference.py
##
##  ---------------------------------------------------------------------------
##  UNITS  (stated once, obeyed everywhere)
##  ---------------------------------------------------------------------------
##    TIME          hours
##    AMOUNTS       mmol
##    CONCENTRATION mmol/L (mM); mg/dL appears only in $TABLE
##    VOLUMES       L            CLEARANCES  L/h
##    PaCO2         mmHg         GFR         fraction of the patient's baseline
##
##  ---------------------------------------------------------------------------
##  CALIBRATION AND VERIFICATION NOTES  (see README.md for the full account)
##  ---------------------------------------------------------------------------
##  Anchored on:
##    · methanol elimination ~8.5 mg/dL/h (range 4.4-25) at saturating
##      concentrations; terminal t1/2 ~2.5 h unblocked
##    · methanol t1/2 43-52 h with ADH blocked (no renal exit)
##    · ethylene glycol t1/2 ~17 h with ADH blocked and intact kidneys
##      (20-30% excreted unchanged) -> the model reproduces 17.3 h
##    · fomepizole 15 mg/kg load then 10 mg/kg q12h, therapeutic >8.2 ug/mL,
##      Ki 0.1-0.5 uM; q4h during haemodialysis
##    · ethanol target 100-150 mg/dL; apparent Km ~1-2 mM
##    · haemodialysis clearance 200-300 mL/min for all four small solutes
##    · Ksp(calcium oxalate monohydrate) 2.32e-9 M^2 at 37 C
##    · apparent bicarbonate space (0.40 + 2.6/[HCO3]) L/kg
##
##  ALL 50 ODEs WERE INDEPENDENTLY RE-IMPLEMENTED IN PYTHON/SCIPY
##  (tap_python_reference.py).  That exercise found and fixed eight defects in
##  this model, three of which changed a conclusion:
##
##   1. Renal bicarbonate regeneration written at 0.10/h corrected any acidosis
##      in a few hours: a lethal methanol dose produced a pH of 7.33.  It is a
##      maximal-NH4+-excretion term and belongs at ~0.010/h.
##   2. Tubular calcium-oxalate deposition written as a rate proportional to the
##      supersaturation EXCESS deposited 50 mmol/h of crystal and drove every
##      simulated ethylene-glycol patient's ionised calcium to zero within
##      hours.  Crystal formation cannot exceed oxalate DELIVERY; it is now a
##      saturating fraction of the filtered load.
##   3. Optic injury driven linearly by vitreous formate blinded every patient
##      who was exposed at a low level for a long time - including the whole
##      ethanol-antidote arm.  The retina fails for the same reason the putamen
##      does, so the driver is a retinal ATP deficit, which makes visual loss a
##      threshold function instead of an integral.
##   4. Proximal tubular injury with no flux threshold destroyed the kidney of
##      patients who received fomepizole at two hours (GFR 50%), because the
##      residual non-ADH flux never stops.
##   5. Haemodialysis bicarbonate transfer written at the small-solute clearance
##      delivered >300 mmol/h and produced HCO3 47 mM and PaCO2 79 mmHg.
##   6. Winter's formula extrapolated above a normal bicarbonate predicted a
##      PaCO2 of 79 mmHg for an alkalotic patient.
##   7. Exogenous NaHCO3 raises sodium as well as bicarbonate; ignoring the
##      sodium drove the computed anion gap to -11 on a bicarbonate infusion.
##      With sodium counted, bicarbonate therapy leaves the anion gap almost
##      unchanged - which is the clinically important behaviour.
##   8. Integrating past pH 6.6 produced numbers (pH 5.55, CNS formate 120 mM)
##      that describe a corpse; runs are now truncated at the first
##      non-survivable point.
##
##  AND VERIFICATION REFUTED TWO OF THE AUTHOR'S OWN CLAIMS, WHICH ARE LEFT
##  STATED AS REFUTED RATHER THAN QUIETLY REMOVED:
##
##   A. "Dialysis removes fomepizole, so the acid flux spikes mid-session."
##      It does not.  What matters is the margin C_therapeutic/K_inhibition:
##      812 for fomepizole, 22 for ethanol.  A 6 h session cannot spend
##      fomepizole's margin (16.4 h of continuous dialysis would be needed) and
##      omitting the q4h top-up changes the dose oxidised by +0.3%.  It CAN
##      spend most of ethanol's, and there the un-boosted arm oxidises at 3.74
##      vs 1.36 mmol/h by the end of the session.  The loop is real only for the
##      low-margin antidote.
##   B. "pH is the best predictor of outcome."  Only against the ALCOHOL level,
##      where it wins decisively (rho -0.84 vs -0.13, and the alcohol level runs
##      the WRONG WAY).  Plasma formate is better still (rho +0.94).  pH rules
##      the bedside because formate cannot be measured in the four hours during
##      which the decision must be made, not because it is the better variable.
##
##  DISCLAIMER: educational / research model.  Not validated for clinical use.
## ===========================================================================

library(mrgsolve)
library(dplyr)

code <- '
$PROB
# Toxic alcohol poisoning (methanol / ethylene glycol) QSP model
- 50 ODEs; time in hours; amounts in mmol; concentrations in mM

$GLOBAL
#define MWMEOH 32.042
#define MWFORM 46.025
#define MWEG   62.068
#define MWGLYC 76.051
#define MWOXAL 90.034
#define MWETOH 46.068
#define MWFOM  82.104

// Fraction of a weak acid present as the un-ionised, membrane-permeant species.
// This is the single most load-bearing function in the model: it is what makes
// acidaemia a DELIVERY mechanism and not merely a number to be corrected.
double fHA(double pH, double pKa) {
  return 1.0 / (1.0 + pow(10.0, pH - pKa));
}

// Apparent bicarbonate space, litres.  Garella/Fernandez-type expansion:
// space (L/kg) = 0.40 + 2.6/[HCO3].  The space GROWS as the acidosis deepens,
// which is why plasma bicarbonate decelerates on its way down - and why the
// two-gap invariant is only approximately conserved.
double bicspace(double hco3, double wt) {
  double h = hco3 > 3.0 ? hco3 : 3.0;
  return (0.40 + 2.6 / h) * wt;
}

double bloodpH(double hco3, double paco2) {
  double h = hco3 > 0.5 ? hco3 : 0.5;
  double p = paco2 > 8.0 ? paco2 : 8.0;
  return 6.10 + log10(h / (0.0301 * p));
}

// Urine pH: a saturating function of plasma bicarbonate plus a switch at the
// renal bicarbonate threshold (bicarbonaturia turns the urine alkaline).
double urinepH(double hco3) {
  double base = 5.0 + 3.0 * (hco3 / (hco3 + 20.0));
  double extra = 1.30 / (1.0 + exp(-(hco3 - 27.0) / 1.2));
  double v = base + extra;
  return v > 8.0 ? 8.0 : v;
}

// Is a time-window therapy active now?
double inwin(double t, double a, double b) {
  return (t >= a && t < b) ? 1.0 : 0.0;
}

double pos(double x) { return x > 0.0 ? x : 0.0; }

$PARAM @annotated
// ---------------------------------------------------------------- patient ----
WT       :  70.0 : body weight (kg)
GFRBASE  :   1.0 : baseline GFR as a fraction of 120 mL/min (CKD < 1)

// ------------------------------------------------------------- absorption ----
KA       :   1.8 : first-order absorption of the alcohols (1/h)
KAFOM    :   2.5 : first-order absorption of oral fomepizole (1/h)
SLOWE    :  0.65 : maximal slowing of gastric emptying by ethanol
KSLOWE   :  15.0 : ethanol concentration for half that effect (mM)

// -------------------------------- ADH: one enzyme, competing substrates -----
VMAX_M   :  95.0 : methanol oxidation capacity (mmol/h)
KM_M     :   8.0 : ADH Km for methanol (mM)
VMAX_E   :  52.0 : ethylene glycol oxidation capacity (mmol/h)
KM_E     :   6.0 : ADH Km for ethylene glycol (mM)
VMAX_ET  : 150.0 : ethanol oxidation capacity (mmol/h)
KM_ET    :   1.0 : ADH Km for ethanol (mM)
KI_FOM   : 1.5E-4: fomepizole competitive Ki (mM) = 0.15 uM
FCAT_M   :  0.06 : catalase/CYP2E1 route for methanol (fraction of Vmax)
FCAT_E   :  0.01 : non-ADH route for ethylene glycol (essentially nil)

// -------------------------------------------------- formate disposition -----
VMAX_THF :  18.0 : folate-dependent formate oxidation capacity (mmol/h)
KM_THF   :   0.5 : Km of that step (mM)
KTHF_USE : 0.0018: THF consumed per mmol of formate oxidised
KTHF_REG : 0.030 : THF pool regeneration (1/h)
CLF_R0   :  0.45 : renal formate clearance scale
FR_MAX   :  0.99 : maximal fractional tubular reabsorption of the acid
KFR      : 2.5E-4: f_HA at which reabsorption is half-maximal
PKA_F    :  3.75 : formic acid pKa
PKA_G    :  3.83 : glycolic acid pKa

// ---------------------------------------------- CNS and vitreous transfer ---
PSF      : 900.0 : permeability x area for un-ionised formic acid, brain (L/h)
PSV      :   8.0 : same, vitreous/retina (L/h)
PHB_SLOPE:  0.50 : fraction of the plasma pH shift seen inside the brain

// --------------------------------------------------- glyoxylate branching ---
K_GLX    :   4.0 : glyoxylate turnover (1/h)
KM_GLXO  :   3.0 : glycolate oxidase Km (mM)
VMAX_GLXO:   5.0 : glycolate oxidase capacity (mmol/h) - deliberately SLOW
BR_OX    :  0.25 : oxalate branch weight (the MINOR route)
BR_TH    :  0.45 : AST / thiamine-dependent branch weight
BR_PY    :  0.30 : ALT / pyridoxine-dependent branch weight
ETH_TH   :  1.10 : maximal thiamine effect on the AST branch
ETH_PY   :  1.00 : maximal pyridoxine effect on the ALT branch
FRAC_EGF :  0.04 : minor ethylene glycol -> formate route

// ------------------------------------------------------ calcium oxalate ----
KSP      : 2.32E-3: Ksp of CaOx monohydrate (mM^2) = 2.32e-9 M^2
CFTUB    :  50.0 : tubular concentration factor for the filtrate
SSCRIT_K :   8.0 : metastable supersaturation limit, tubular fluid
SSCRIT_T :   4.0 : metastable supersaturation limit, soft tissue
FP_MAX   :  0.60 : maximal precipitating fraction of the filtered oxalate load
KSS_K    : 250.0 : tubular SS excess giving half of FP_MAX
KDEP_T   :  0.10 : soft-tissue precipitation rate (1/h)
KSS_T    :  20.0 : plasma SS excess giving half of that
BUFCA    :  0.40 : fraction of deposited calcium seen as an ionised decrement
KCA_REP  : 0.055 : bone / PTH / gut restoration of ionised calcium (1/h)
CA0      :  1.18 : baseline ionised calcium (mM)
DCA_DPH  :  0.36 : ionised calcium fall per unit rise in pH (mM)

// ------------------------------------------------------------- acid-base ---
NAB      : 140.0 : baseline plasma sodium (mM)
CLB      : 104.0 : plasma chloride, held fixed (mM)
HCO30    :  24.0 : baseline bicarbonate (mM)
AG0      :  12.0 : baseline anion gap (mEq/L)
TAU_CO2  :  0.40 : respiratory response time constant (h)
PACO20   :  40.0 : baseline PaCO2 (mmHg)
KHCO3_REN: 0.010 : renal bicarbonate regeneration, NH4+ limited (1/h)
HCO3_DIA :  35.0 : dialysate bicarbonate (mM)
CL_HD_BIC:   3.0 : EFFECTIVE whole-body bicarbonate transfer on HD (L/h)
PH_TARGET:  7.35 : pH at which the bicarbonate infusion is titrated off
PH_TAU   :  0.03 : sharpness of that titration
VNA      :  0.60 : L/kg over which an exogenous sodium load spreads
KNA_OUT  : 0.014 : renal washout of the exogenous alkali load (1/h)

// ---------------------------------------------------------------- lactate ---
LACT0    :   1.0 : baseline lactate (mM)
KLACT_OUT:  0.60 : lactate elimination (1/h)
KLACT_ATP:   5.0 : anaerobic production gain per unit ATP deficit
KLACT_ETH: 0.020 : NADH-driven rise per mM ethanol

// ---------------------------------------------- mitochondria / CNS injury ---
KI_CCO   :  12.0 : formate Ki at cytochrome c oxidase (mM) - deliberately WEAK
TAU_ATP  :  0.50 : ATP pool time constant (h)
THR_PUT  :  0.25 : ATP deficit threshold for putaminal injury
K_PUT    :  0.55 : putaminal injury rate (1/h)
K_PUTREC : 0.004 : putaminal recovery (1/h)
THR_OPTA :  0.18 : retinal ATP deficit threshold for optic injury
K_OPT    :  0.40 : optic injury rate (1/h)
K_OPTREC : 0.006 : optic recovery (1/h)
LOGMARMX :   2.0 : logMAR at complete optic injury

// -------------------------------------------------------------- kidney ------
K_PTGALD : 0.0025: tubular injury per (mmol/h) of EG oxidation flux
VE_THR   :   6.0 : EG flux the tubule detoxifies without injury (mmol/h)
K_PTXTAL :  0.35 : tubular injury per (mmol/h) of CaOx deposition
DEP_THR  :  0.70 : CaOx deposition the tubule passes without injury (mmol/h)
K_PTREC  : 0.040 : tubular repair (1/h)
K_GFRINJ :  0.22 : GFR loss rate at full tubular injury (1/h)
K_GFRREP : 0.035 : GFR recovery (1/h)
GFR_MIN  :  0.05 : anuric floor

// ------------------------------------------------- sedation / respiration ---
KSED_M   :  60.0 : methanol concentration giving unit sedation score (mM)
KSED_E   :  55.0 : same, ethylene glycol
KSED_ET  :  45.0 : same, ethanol
KSED_PH  :  0.55 : pH units below 7.20 giving unit sedation score
KSED_ATP :  0.45 : ATP deficit giving unit sedation score
PACO2FAIL:  65.0 : PaCO2 target in the fully obtunded patient (mmHg)
COMA_RESP:  0.72 : coma index at which respiratory failure sets in

// -------------------------------------------------------------- hazards -----
H_PH     :  0.16 : mortality hazard per (7.20-pH)^2 unit (1/h)
H_PUT    : 0.030 : mortality hazard per unit putaminal injury (1/h)
H_RESP   : 0.055 : mortality hazard per unit respiratory failure (1/h)
H_BLIND  :  0.85 : blindness hazard per unit optic injury above threshold (1/h)
BLIND_THR:  0.18 : that threshold
PH_FATAL :  6.60 : pH below which the run is flagged non-survivable

// --------------------------------------------------------- fomepizole PK ----
VMAX_FOM :  0.90 : fomepizole elimination capacity (mmol/h)
KM_FOM   :  0.05 : Km of that step (mM)
AUTOIND  :  0.60 : fractional rise in that capacity by autoinduction
T_AUTO   :  30.0 : time at which autoinduction is half-on (h)
TAU_AUTO :   8.0 : sharpness of the autoinduction switch (h)

// --------------------------------------------------------- cofactor PD ------
KOUT_FOL :  0.14 : folinic-acid effect washout (1/h)
EMAX_FOL :  0.45 : maximal fractional rise in VMAX_THF
KOUT_TH  :  0.10 : thiamine effect washout (1/h)
KOUT_PY  :  0.10 : pyridoxine effect washout (1/h)

// -------------------------------------------------------- extracorporeal ----
CL_HD    :  14.4 : intermittent HD clearance, small solutes (L/h) = 240 mL/min
CL_HD_FOM:  12.0 : HD clearance of fomepizole (L/h)
CL_HD_ET :  13.0 : HD clearance of ethanol (L/h)
CL_CRRT  :   2.4 : CVVHD/CVVHDF clearance (L/h) = 40 mL/min

$PARAM @annotated
// ============================ THERAPY WINDOWS (a scenario is a row) =========
// Fomepizole and ethanol BOLUSES are given as dosing events into the
// corresponding compartments; everything continuous is a window below.
HD1ON    : 1E6 : start of the first haemodialysis session (h)
HD1OFF   : 1E6 : end of it (h)
HD2ON    : 1E6 : start of a second session (h)
HD2OFF   : 1E6 : end of it (h)
CRON     : 1E6 : start of CRRT (h)
CROFF    : 1E6 : end of CRRT (h)
BICON    : 1E6 : start of the bicarbonate infusion (h)
BICOFF   : 1E6 : end of it (h)
BICRATE  : 0.0 : bicarbonate infusion rate (mmol/h)
ETHON    : 1E6 : start of the ethanol maintenance infusion (h)
ETHOFF   : 1E6 : end of it (h)
ETHRATE  : 0.0 : ethanol maintenance rate (mmol/h)
ETHHDMLT : 2.5 : multiplier applied to that rate while on dialysis
FOMION   : 1E6 : start of a continuous fomepizole infusion (h)
FOMIOFF  : 1E6 : end of it (h)
FOMIRATE : 0.0 : that rate (mmol/h)
FOLON    : 1E6 : start of folinic acid (h)
FOLOFF   : 1E6 : end of it (h)
FOLRATE  : 0.20: folinic-acid input (arbitrary effect units/h)
THION    : 1E6 : start of thiamine (h)
THIOFF   : 1E6 : end of it (h)
THIRATE  : 0.16: thiamine input (effect units/h)
PYRON    : 1E6 : start of pyridoxine (h)
PYROFF   : 1E6 : end of it (h)
PYRRATE  : 0.16: pyridoxine input (effect units/h)
CAION_ON : 1E6 : start of calcium replacement (h)
CAION_OF : 1E6 : end of it (h)
CARATE   : 0.0 : elemental calcium infusion (mmol/h)
NOACID   : 0.0 : 1 = isopropanol control (ADH runs, metabolite is a KETONE)

$CMT @annotated
// ------------------------------------------------------------ methanol arm --
A_GUTM   : methanol in the gut lumen (mmol)
A_MEOH1  : methanol, central compartment (mmol)
A_MEOH2  : methanol, peripheral compartment (mmol)
A_FORM1  : formate, central (mmol) - THE methanol toxin
A_FORM2  : formate, peripheral (mmol)
A_FCNS   : formate in brain water (mmol) - the quantity that injures
A_FVIT   : formate in vitreous/retina (mmol)
THF      : hepatic tetrahydrofolate pool, relative (1 = normal)
// ----------------------------------------------------- ethylene glycol arm --
A_GUTE   : ethylene glycol in the gut lumen (mmol)
A_EG1    : ethylene glycol, central (mmol)
A_EG2    : ethylene glycol, peripheral (mmol)
A_GLYC1  : glycolate, central (mmol) - THE ethylene glycol acidosis
A_GLYC2  : glycolate, peripheral (mmol)
A_GLX    : glyoxylate (mmol) - the branch node
A_OXAL   : oxalate in the ECF (mmol)
A_CAOXK  : calcium oxalate deposited in the kidney (mmol, cumulative)
A_CAOXT  : calcium oxalate deposited in soft tissue (mmol, cumulative)
// ------------------------------------------------------------- antidotes ----
A_GUTF   : fomepizole in the gut lumen (mmol)
A_FOM1   : fomepizole, central (mmol)
A_FOM2   : fomepizole, peripheral (mmol)
A_GUTETH : ethanol in the gut lumen (mmol)
A_ETOH1  : ethanol, central (mmol)
A_ETOH2  : ethanol, peripheral (mmol)
FOLEF    : folinic-acid effect on VMAX_THF
THIEF    : thiamine effect on the AST branch
PYREF    : pyridoxine effect on the ALT branch
// --------------------------------------------------- acid-base/electrolyte --
HCO3     : plasma bicarbonate (mM)
LACT     : plasma L-lactate (mM)
PACO2    : arterial PaCO2 (mmHg)
CAION    : ionised calcium pool, before the pH correction (mM)
// ------------------------------------------------------------ organ state ---
GFRF     : GFR as a fraction of this patient baseline
PTINJ    : proximal tubular injury (0-1)
ATPC     : CNS ATP / energy charge (1 = normal)
PUT      : putaminal injury index (0-1)
OPTIC    : optic / retinal injury index (0-1)
// ------------------------------------------- cumulative and hazard states ---
AUC_FCNS : CNS formate exposure integral (mM*h)
AUC_GLYC : plasma glycolate exposure integral (mM*h)
HAZ_MORT : cumulative mortality hazard
HAZ_BLND : cumulative blindness hazard
CUM_OX_M : methanol oxidised (mmol)
CUM_OX_E : ethylene glycol oxidised (mmol)
CUM_HD_M : methanol removed by dialysis (mmol)
CUM_HD_F : formate removed by dialysis (mmol)
CUM_HD_E : ethylene glycol removed by dialysis (mmol)
CUM_HD_G : glycolate removed by dialysis (mmol)
CUM_HD_FM: fomepizole removed by dialysis (mmol)
CUM_UR_F : formate excreted in urine (mmol)
CUM_UR_G : glycolate excreted in urine (mmol)
CUM_UR_OX: oxalate excreted in urine (mmol)
CUM_BIC  : exogenous NaHCO3 delivered (mmol)

$MAIN
THF_0    = 1.0;
HCO3_0   = HCO30;
LACT_0   = LACT0;
PACO2_0  = PACO20;
CAION_0  = CA0;
GFRF_0   = 1.0;
ATPC_0   = 1.0;
A_OXAL_0 = 0.002 * 0.25 * WT;   // ~2 uM endogenous plasma oxalate

$ODE
// =========================================================================
//  0 · volumes and concentrations
// =========================================================================
double V1M   = 0.35 * WT;      double V2M   = 0.25 * WT;   // total 0.60 L/kg
double V1E   = 0.38 * WT;      double V2E   = 0.27 * WT;   // total 0.65 L/kg
double V1F   = 0.30 * WT;      double V2F   = 0.20 * WT;   // total 0.50 L/kg
double V1G   = 0.28 * WT;      double V2G   = 0.20 * WT;
double VCNS  = 0.020 * WT;     double VVIT  = 0.05;
double VECF  = 0.20 * WT;      double VOX   = 0.25 * WT;
double V1FOM = 0.42 * WT;      double V2FOM = 0.28 * WT;   // total 0.70 L/kg
double V1ET  = 0.35 * WT;      double V2ET  = 0.25 * WT;
double VLACT = 0.35 * WT;
double GFR0  = 7.2 * GFRBASE;                              // L/h

double C_M   = A_MEOH1 / V1M;   double C_M2  = A_MEOH2 / V2M;
double C_E   = A_EG1   / V1E;   double C_E2  = A_EG2   / V2E;
double C_F   = A_FORM1 / V1F;   double C_F2  = A_FORM2 / V2F;
double C_FC  = A_FCNS  / VCNS;  double C_FV  = A_FVIT  / VVIT;
double C_G   = A_GLYC1 / V1G;   double C_G2  = A_GLYC2 / V2G;
double C_GX  = A_GLX   / VECF;  double C_OX  = A_OXAL  / VOX;
double C_FOM = A_FOM1  / V1FOM; double C_FOM2= A_FOM2  / V2FOM;
double C_ET  = A_ETOH1 / V1ET;  double C_ET2 = A_ETOH2 / V2ET;

double hco3  = HCO3  > 0.5 ? HCO3  : 0.5;
double paco2 = PACO2 > 8.0 ? PACO2 : 8.0;
double pH    = bloodpH(hco3, paco2);
double pHu   = urinepH(hco3);
double gfrf  = GFRF > 0.02 ? GFRF : 0.02;
double gfr   = GFR0 * gfrf;

// ionised calcium AS MEASURED: albumin binding shifts with pH
double ca_meas = CAION - DCA_DPH * (pH - 7.40);
if (ca_meas < 0.30) ca_meas = 0.30;

// =========================================================================
//  1 · therapy switches
// =========================================================================
double onHD   = inwin(SOLVERTIME, HD1ON, HD1OFF) + inwin(SOLVERTIME, HD2ON, HD2OFF);
double onCR   = inwin(SOLVERTIME, CRON, CROFF);
if (onHD > 1.0) onHD = 1.0;
double cl_hd     = onHD * CL_HD     + (onHD > 0 ? 0.0 : onCR * CL_CRRT);
double cl_hd_fom = onHD * CL_HD_FOM + (onHD > 0 ? 0.0 : onCR * CL_CRRT);
double cl_hd_et  = onHD * CL_HD_ET  + (onHD > 0 ? 0.0 : onCR * CL_CRRT);
double anyEC     = (cl_hd > 0.0) ? 1.0 : 0.0;

double bic_inf = BICRATE  * inwin(SOLVERTIME, BICON,    BICOFF);
double fom_inf = FOMIRATE * inwin(SOLVERTIME, FOMION,   FOMIOFF);
double eth_inf = ETHRATE  * inwin(SOLVERTIME, ETHON,    ETHOFF)
                 * (anyEC > 0 ? ETHHDMLT : 1.0);
double ca_inf  = CARATE   * inwin(SOLVERTIME, CAION_ON, CAION_OF);
double fol_inf = FOLRATE  * inwin(SOLVERTIME, FOLON,    FOLOFF);
double thi_inf = THIRATE  * inwin(SOLVERTIME, THION,    THIOFF);
double pyr_inf = PYRRATE  * inwin(SOLVERTIME, PYRON,    PYROFF);

// =========================================================================
//  2 · ADH — ONE enzyme, THREE competing substrates, ONE inhibitor
//      v_i = Vmax_i*(S_i/Km_i) / (1 + SUM_j S_j/Km_j + F/Ki)
//      Fomepizole and ethanol own the DENOMINATOR and nothing else.  Because
//      the block is COMPETITIVE it is diluted by its own substrate, which is
//      why "fomepizole is 1000x ethanol" is false at the bedside: at a
//      methanol of 100 mg/dL an ethanol infusion at target still leaves ~18%
//      of the flux running and fomepizole leaves 0.6%.
// =========================================================================
double sM  = C_M  / KM_M;
double sE  = C_E  / KM_E;
double sET = C_ET / KM_ET;
double inh = C_FOM / KI_FOM;
double den = 1.0 + sM + sE + sET + inh;

double v_M  = VMAX_M  * sM  / den;
double v_E  = VMAX_E  * sE  / den;
double v_ET = VMAX_ET * sET / den;
// non-ADH routes are NOT blocked by fomepizole - this is why blockade is never
// absolute, and why 14% of a dose is still oxidised over four days of a
// perfect block
v_M += FCAT_M * VMAX_M * C_M / (KM_M * 6.0 + C_M);
v_E += FCAT_E * VMAX_E * C_E / (KM_E * 6.0 + C_E);

// =========================================================================
//  3 · absorption and parent-alcohol disposition
//      THE ONE CLEARANCE TERM THAT SPLITS THE TWO DISEASES is CLRE: 20-30% of
//      ethylene glycol leaves unchanged in urine, methanol essentially none.
//      Same antidote, same mechanism, opposite management.
// =========================================================================
double ka    = KA / (1.0 + SLOWE * C_ET / (KSLOWE + C_ET));
double abs_M = ka * A_GUTM;
double abs_E = ka * A_GUTE;

double q_M = 30.0 * (C_M - C_M2);
double q_E = 30.0 * (C_E - C_E2);
double cl_ren_M  = 0.050 * gfr;      // freely filtered, then reabsorbed
double cl_pulm_M = 0.15;             // blood:air ~1350:1
double cl_ren_E  = 0.260 * gfr;      // the term that changes management

dxdt_A_GUTM  = -abs_M;
dxdt_A_GUTE  = -abs_E;
dxdt_A_MEOH1 = abs_M - q_M - v_M - (cl_ren_M + cl_pulm_M) * C_M - cl_hd * C_M;
dxdt_A_MEOH2 = q_M;
dxdt_A_EG1   = abs_E - q_E - v_E - cl_ren_E * C_E - cl_hd * C_E;
dxdt_A_EG2   = q_E;

// =========================================================================
//  4 · formate: production, folate-dependent oxidation, renal, CNS, vitreous
// =========================================================================
double vmax_thf = VMAX_THF * (1.0 + EMAX_FOL * FOLEF);
double ox_form  = vmax_thf * THF * C_F / (KM_THF + C_F);

// renal handling is non-ionic back-diffusion of the UN-IONISED acid, so the
// reabsorbed fraction is a function of URINE pH.  This is the second half of
// what bicarbonate does, and it runs in the same direction as the first.
double fha_u  = fHA(pHu, PKA_F);
double fr_F   = FR_MAX * fha_u / (fha_u + KFR);
double cl_rF  = gfr * (1.0 - fr_F) * CLF_R0;
double exc_F  = cl_rF * C_F;

double fha_p = fHA(pH, PKA_F);
double phc   = pH > 6.60 ? pH : 6.60;         // do not extrapolate the offsets
double pH_b  = 7.00 + PHB_SLOPE * (phc - 7.40);
double pH_v  = 7.10 + PHB_SLOPE * (phc - 7.40);
double j_cns = PSF * (fha_p * C_F - fHA(pH_b, PKA_F) * C_FC);
double j_vit = PSV * (fha_p * C_F - fHA(pH_v, PKA_F) * C_FV);

double q_F         = 8.0 * (C_F - C_F2);
double form_from_e = FRAC_EGF * v_E;
double acidM       = (NOACID > 0.5) ? 0.0 : 1.0;   // isopropanol control

dxdt_A_FORM1 = acidM * (v_M + form_from_e) - q_F - ox_form - exc_F
               - j_cns - j_vit - cl_hd * C_F;
dxdt_A_FORM2 = q_F;
dxdt_A_FCNS  = j_cns;
dxdt_A_FVIT  = j_vit;
dxdt_THF     = -KTHF_USE * ox_form + KTHF_REG * (1.0 - THF);

// =========================================================================
//  5 · glycolate, the glyoxylate BRANCH RATIO, oxalate
//      Thiamine and pyridoxine widen two escape routes out of ONE node, so
//      what they own is a RATIO applied to a flux.  Close the flux and the
//      ratio has nothing to divide - the predicted non-additivity with PROMPT
//      fomepizole (and its survival alongside a LATE one).
// =========================================================================
double q_G    = 8.0 * (C_G - C_G2);
double v_glxo = VMAX_GLXO * C_G / (KM_GLXO + C_G);
double fha_ug = fHA(pHu, PKA_G);
double fr_G   = FR_MAX * fha_ug / (fha_ug + KFR);
double exc_G  = gfr * (1.0 - fr_G) * C_G;

dxdt_A_GLYC1 = v_E * (1.0 - FRAC_EGF) - q_G - v_glxo - exc_G - cl_hd * C_G;
dxdt_A_GLYC2 = q_G;

double b_ox = BR_OX;
double b_th = BR_TH * (1.0 + ETH_TH * THIEF);
double b_py = BR_PY * (1.0 + ETH_PY * PYREF);
double btot = b_ox + b_th + b_py;
double out_glx = K_GLX * A_GLX;
dxdt_A_GLX = v_glxo - out_glx;
double to_oxal     = out_glx * b_ox / btot;
double neutralised = out_glx * (b_th + b_py) / btot;

// Calcium oxalate: the model s ONLY threshold, and therefore the only place
// where an equal dose delivered SLOWLY does less harm.  Deposition is written
// as a saturating fraction of the FILTERED LOAD, because crystal formation
// cannot outrun oxalate delivery.
double ss_p    = ca_meas * C_OX / KSP;
double ss_k    = ss_p * CFTUB;
double ex_k    = pos(ss_k - SSCRIT_K);
double ex_t    = pos(ss_p - SSCRIT_T);
double filt_ox = gfr * C_OX;
double fprec   = FP_MAX * ex_k / (ex_k + KSS_K);
double dep_k   = fprec * filt_ox;
double exc_OX  = (1.0 - fprec) * filt_ox;
double dep_t   = KDEP_T * ex_t / (ex_t + KSS_T) * A_OXAL;

dxdt_A_OXAL  = to_oxal - dep_k - dep_t - exc_OX;
dxdt_A_CAOXK = dep_k;
dxdt_A_CAOXT = dep_t;

// =========================================================================
//  6 · antidote pharmacokinetics
// =========================================================================
double autoi  = 1.0 + AUTOIND / (1.0 + exp(-(SOLVERTIME - T_AUTO) / TAU_AUTO));
double el_fom = VMAX_FOM * autoi * C_FOM / (KM_FOM + C_FOM);
double abs_F  = KAFOM * A_GUTF;
double q_FOM  = 25.0 * (C_FOM - C_FOM2);
dxdt_A_GUTF = -abs_F;
dxdt_A_FOM1 = abs_F + fom_inf - q_FOM - el_fom - cl_hd_fom * C_FOM;
dxdt_A_FOM2 = q_FOM;

double abs_ET = ka * A_GUTETH;
double q_ET   = 40.0 * (C_ET - C_ET2);
dxdt_A_GUTETH = -abs_ET;
dxdt_A_ETOH1  = abs_ET + eth_inf - q_ET - v_ET - cl_hd_et * C_ET;
dxdt_A_ETOH2  = q_ET;

dxdt_FOLEF = fol_inf - KOUT_FOL * FOLEF;
dxdt_THIEF = thi_inf - KOUT_TH  * THIEF;
dxdt_PYREF = pyr_inf - KOUT_PY  * PYREF;

// =========================================================================
//  7 · acid-base — organic anion equivalents RETAINED
//      +1 per formate produced, +1 per glycolate produced, +1 more per oxalate
//      produced (glyoxylate 1- -> oxalate 2-); -1 per formate oxidised to CO2,
//      -1 per anion excreted or dialysed, -1 per glyoxylate routed to glycine
//      or alpha-hydroxy-beta-ketoadipate, -2 per oxalate excreted.
//      Oxalate that PRECIPITATES is deliberately not credited back: its proton
//      has already gone as CO2, so the bicarbonate deficit stays.
// =========================================================================
double lact_prod = KLACT_OUT * LACT0 * VLACT
                   + KLACT_ATP * pos(1.0 - ATPC) * VLACT
                   + KLACT_ETH * C_ET * VLACT;
double lact_out  = KLACT_OUT * LACT * VLACT;
dxdt_LACT = (lact_prod - lact_out) / VLACT;

double acid_in  = acidM * (v_M + form_from_e)
                  + v_E * (1.0 - FRAC_EGF) + to_oxal;
double acid_out = ox_form + exc_F + exc_G + cl_hd * (C_F + C_G)
                  + neutralised + 2.0 * exc_OX;
// In the isopropanol CONTROL the metabolite is a ketone, so acidM has already
// removed BOTH the anion and the proton from acid_in: same ADH flux, same
// osmolal gap, no acid, and every downstream lesion must therefore vanish.
double net_acid = acid_in - acid_out + (lact_prod - lact_out);

double Vb  = bicspace(hco3, WT);
double _e  = (pH - PH_TARGET) / PH_TAU;
if (_e >  50.0) _e =  50.0;
if (_e < -50.0) _e = -50.0;
double bic_eff = bic_inf / (1.0 + exp(_e));            // titrated, not poured
double hd_bic  = anyEC * CL_HD_BIC * (HCO3_DIA - hco3);
double ren_bic = KHCO3_REN * gfrf * (HCO30 - hco3) * Vb;
dxdt_HCO3 = (-net_acid + bic_eff + hd_bic + ren_bic) / Vb;
dxdt_CUM_BIC = bic_eff + pos(hd_bic) - KNA_OUT * gfrf * CUM_BIC;

// respiratory compensation, and its failure as obtundation deepens.
// Winter s formula is an ACIDOSIS rule; above a normal bicarbonate the
// alkalosis rule (~0.7 mmHg per mM) applies instead.
double sed = C_M / KSED_M + C_E / KSED_E + C_ET / KSED_ET
             + pos(7.20 - pH) / KSED_PH + pos(1.0 - ATPC) / KSED_ATP;
double coma  = sed / (1.0 + sed);
double respf = 1.0 / (1.0 + exp(-(coma - COMA_RESP) / 0.06));
double tgt_comp;
if (hco3 <= HCO30) {
  tgt_comp = 1.5 * hco3 + 8.0;
  if (tgt_comp < 12.0) tgt_comp = 12.0;
} else {
  tgt_comp = 40.0 + 0.7 * (hco3 - HCO30);
  if (tgt_comp > 55.0) tgt_comp = 55.0;
}
double tgt = tgt_comp * (1.0 - respf) + PACO2FAIL * respf;
dxdt_PACO2 = (tgt - paco2) / TAU_CO2;

// calcium: the STOICHIOMETRIC SHADOW of the crystal burden
dxdt_CAION = -(dep_k + dep_t) * BUFCA / VECF + ca_inf * BUFCA / VECF
             + KCA_REP * (CA0 - CAION);

// =========================================================================
//  8 · organ injury
// =========================================================================
double hl_ald = pow(v_E, 3.0) / (pow(v_E, 3.0) + pow(VE_THR, 3.0));
double hl_xt  = pow(dep_k, 3.0) / (pow(dep_k, 3.0) + pow(DEP_THR, 3.0));
double inj    = K_PTGALD * hl_ald * v_E + K_PTXTAL * hl_xt * dep_k;
dxdt_PTINJ = inj * (1.0 - PTINJ) - K_PTREC * PTINJ;
dxdt_GFRF  = -K_GFRINJ * PTINJ * PTINJ * (gfrf - GFR_MIN)
             + K_GFRREP * (1.0 - gfrf) * (1.0 - 0.70 * PTINJ);

double atp_tgt = 1.0 / (1.0 + C_FC / KI_CCO + C_GX / 8.0);
dxdt_ATPC = (atp_tgt - ATPC) / TAU_ATP;

double deficit = pos((1.0 - ATPC) - THR_PUT);
dxdt_PUT = K_PUT * deficit * (1.0 - PUT) - K_PUTREC * PUT;
// the retina fails for the same reason the putamen does, which makes visual
// loss a THRESHOLD function of vitreous formate rather than an integral of it
double atp_ret = 1.0 / (1.0 + C_FV / KI_CCO);
double def_opt = pos((1.0 - atp_ret) - THR_OPTA);
dxdt_OPTIC = K_OPT * def_opt * (1.0 - OPTIC) - K_OPTREC * OPTIC;

// =========================================================================
//  9 · hazards and mass-balance bookkeeping
// =========================================================================
dxdt_AUC_FCNS = C_FC;
dxdt_AUC_GLYC = C_G;
dxdt_HAZ_MORT = H_PH * pos(7.20 - pH) * pos(7.20 - pH) * 25.0
                + H_PUT * PUT + H_RESP * respf;
dxdt_HAZ_BLND = H_BLIND * pos(OPTIC - BLIND_THR);
dxdt_CUM_OX_M  = v_M;
dxdt_CUM_OX_E  = v_E;
dxdt_CUM_HD_M  = cl_hd * C_M;
dxdt_CUM_HD_F  = cl_hd * C_F;
dxdt_CUM_HD_E  = cl_hd * C_E;
dxdt_CUM_HD_G  = cl_hd * C_G;
dxdt_CUM_HD_FM = cl_hd_fom * C_FOM;
dxdt_CUM_UR_F  = exc_F;
dxdt_CUM_UR_G  = exc_G;
dxdt_CUM_UR_OX = exc_OX;

$TABLE
double V1M_   = 0.35 * WT;   double V1E_  = 0.38 * WT;
double V1F_   = 0.30 * WT;   double V1G_  = 0.28 * WT;
double VCNS_  = 0.020 * WT;  double VVIT_ = 0.05;
double VOX_   = 0.25 * WT;   double V1FM_ = 0.42 * WT;
double V1ET_  = 0.35 * WT;

double cM  = A_MEOH1 / V1M_;
double cE  = A_EG1   / V1E_;
double cF  = A_FORM1 / V1F_;
double cG  = A_GLYC1 / V1G_;
double cFC = A_FCNS  / VCNS_;
double cFV = A_FVIT  / VVIT_;
double cOX = A_OXAL  / VOX_;
double cFM = A_FOM1  / V1FM_;
double cET = A_ETOH1 / V1ET_;

double hc = HCO3  > 0.5 ? HCO3  : 0.5;
double pc = PACO2 > 8.0 ? PACO2 : 8.0;
double PHA = bloodpH(hc, pc);

// --------------------------------------------- clinical concentrations -----
capture MEOH   = cM  * MWMEOH / 10.0;          // mg/dL
capture EG     = cE  * MWEG   / 10.0;          // mg/dL
capture ETOH   = cET * MWETOH / 10.0;          // mg/dL
capture FORMmM = cF;
capture FORMmgL= cF * MWFORM;
capture GLYCmM = cG;
capture FCNS   = cFC;
capture FVIT   = cFV;
capture OXALuM = cOX * 1000.0;
capture FOMug  = cFM * MWFOM;

// -------------------------------------------------------- acid-base --------
capture pHart  = PHA;
capture HCO3o  = hc;
capture PACO2o = pc;
capture pHur   = urinepH(hc);
capture LACTo  = LACT;

// ------------------------------ THE TWO GAPS, AND THE CLOCK BETWEEN THEM ---
// Exogenous NaHCO3 raises sodium as well as bicarbonate, so with the sodium
// counted the anion gap is almost UNCHANGED by bicarbonate therapy: alkali
// fixes the pH without erasing the diagnostic gap.
double dNa = CUM_BIC / (VNA * WT);
capture NAo = NAB + dNa;
capture AG  = NAo - CLB - hc;
double unmeas = cM + cE + cET + cF + cG + pos(LACT - LACT0) + cOX;
capture OG        = unmeas + (hc - HCO30) - dNa;
capture OG_ETcorr = OG - cET;
capture DAG       = AG - AG0;
capture GAPSUM    = OG + DAG;                  // the approximate dose invariant
capture OGoverAG  = OG / (DAG > 0.5 ? DAG : 0.5);   // THE CLOCK

// --------------------------------- how complete is the blockade, really? ---
double sM_ = cM / KM_M;  double sE_ = cE / KM_E;
double sET_= cET / KM_ET; double inh_ = cFM / KI_FOM;
capture ADHDEN  = 1.0 + sM_ + sE_ + sET_ + inh_;
double base_    = 1.0 + sM_ + sE_;
capture INHFAC  = ADHDEN / base_;
capture FLUXM   = VMAX_M * sM_ / ADHDEN;
capture FLUXE   = VMAX_E * sE_ / ADHDEN;
capture FLUXLEFT= 100.0 / INHFAC;              // % of the flux still running
capture FOMKIRAT= cFM / KI_FOM;                // the margin that matters on HD

// ------------------------------------------- ion trapping: the mechanism ---
capture fHAp    = fHA(PHA, PKA_F);
capture fHAratio= fHAp / fHA(7.40, PKA_F);
double phc_     = PHA > 6.60 ? PHA : 6.60;
capture CNSEQ   = fHAp / fHA(7.00 + PHB_SLOPE * (phc_ - 7.40), PKA_F);

// ------------------------------------- calcium, crystals, kidney, heart ----
double cameas = CAION - DCA_DPH * (PHA - 7.40);
if (cameas < 0.30) cameas = 0.30;
capture CAIONo  = cameas;
capture CATOTmg = cameas / 0.5 * 4.008;
capture SSPLAS  = cameas * cOX / KSP;
capture SSTUB   = SSPLAS * CFTUB;
capture CAOXKID = A_CAOXK;
capture GFRpct  = 100.0 * GFRF;
capture PTINJo  = PTINJ;
capture QTCms   = 400.0 + 95.0 * pos(CA0 - cameas);

// -------------------------------------------- CNS, eye, and the outcomes ---
capture ATPo    = ATPC;
capture PUTo    = PUT;
capture OPTICo  = OPTIC;
capture LOGMAR  = LOGMARMX * OPTIC;
double sed_ = cM / KSED_M + cE / KSED_E + cET / KSED_ET
              + pos(7.20 - PHA) / KSED_PH + pos(1.0 - ATPC) / KSED_ATP;
capture COMA    = sed_ / (1.0 + sed_);
capture PDEATH  = 1.0 - exp(-HAZ_MORT);
capture PBLIND  = 1.0 - exp(-HAZ_BLND);
capture FATAL   = (PHA < PH_FATAL) ? 1.0 : 0.0;

// ---------------------------------------------- diagnostic curiosities -----
// Glycolate cross-reacts on some lactate-oxidase point-of-care electrodes, so
// the blood-gas lactate reads far above the laboratory lactate.  That
// "lactate gap" is a real, free clue to an ethylene glycol ingestion.
capture LACTPOC = LACT + 0.62 * cG;
capture LACTGAP = LACTPOC - LACT;
capture THFo    = THF;
capture MEOHOX  = CUM_OX_M;
capture EGOX    = CUM_OX_E;
'

## ---------------------------------------------------------------------------
tap_model <- function() mcode("tap", code)

## ===========================================================================
##  SCENARIO LIBRARY (25 scenarios)
##
##  Fomepizole and ethanol boluses are dosing events; everything continuous is
##  a parameter window.  The fomepizole schedule is generated by a helper so
##  that the q4h compression during dialysis can be switched off, which is what
##  refuted thesis statement (A) above.
## ===========================================================================
BIG <- 1e6

## methanol / EG / ethanol / fomepizole loading doses as mmol
mmol_meoh <- function(g_per_kg, wt = 70) g_per_kg * wt * 1000 / 32.042
mmol_eg   <- function(mL)                mL * 1.1132 * 1000 / 62.068
mmol_etoh <- function(g_per_kg, wt = 70) g_per_kg * wt * 1000 / 46.068
mmol_fom  <- function(mg_per_kg, wt = 70) mg_per_kg * wt / 82.104

#' Build the fomepizole event table.
#' @param start first (loading) dose time, h
#' @param hd list of c(on, off) dialysis windows
#' @param boost TRUE = compress the interval to q4h during dialysis
fom_events <- function(start, tend, hd = list(), boost = TRUE, wt = 70,
                       load_mgkg = 15, maint_mgkg = 10, q = 12, hq = 4) {
  ev <- data.frame(time = start, amt = mmol_fom(load_mgkg, wt), cmt = "A_FOM1")
  on_hd <- function(t) any(vapply(hd, function(w) t >= w[1] && t < w[2], TRUE))
  tt <- start + q; n <- 0
  while (tt < tend) {
    d <- mmol_fom(maint_mgkg, wt) * if (n >= 4) 1.5 else 1.0
    ev <- rbind(ev, data.frame(time = tt, amt = d, cmt = "A_FOM1"))
    n <- n + 1
    tt <- tt + if (on_hd(tt) && boost) hq else q
  }
  if (boost) for (w in hd) {
    tt <- w[1] + hq
    while (tt < w[2]) {
      ev <- rbind(ev, data.frame(time = tt, amt = mmol_fom(maint_mgkg, wt),
                                 cmt = "A_FOM1"))
      tt <- tt + hq
    }
  }
  ev[order(ev$time), ]
}

#' One scenario = a parameter list + an event table + a run length.
scen <- function(label, tend = 96, wt = 70, gfr = 1,
                 meoh_gkg = 0, eg_mL = 0, etoh_gkg = 0,
                 fom_start = NA, fom_boost = TRUE,
                 eth_start = NA, eth_load_gkg = 0.7, eth_rate_mgkgh = 100,
                 eth_stop = NA, eth_hd_mult = 2.5,
                 hd = list(), crrt = NULL,
                 bic = NULL, fol = NULL, thi = NULL, pyr = NULL, ca = NULL,
                 noacid = 0) {
  p <- list(WT = wt, GFRBASE = gfr, NOACID = noacid,
            HD1ON = BIG, HD1OFF = BIG, HD2ON = BIG, HD2OFF = BIG,
            CRON = BIG, CROFF = BIG, BICON = BIG, BICOFF = BIG, BICRATE = 0,
            ETHON = BIG, ETHOFF = BIG, ETHRATE = 0, ETHHDMLT = eth_hd_mult,
            FOMION = BIG, FOMIOFF = BIG, FOMIRATE = 0,
            FOLON = BIG, FOLOFF = BIG, THION = BIG, THIOFF = BIG,
            PYRON = BIG, PYROFF = BIG, CAION_ON = BIG, CAION_OF = BIG,
            CARATE = 0)
  if (length(hd) >= 1) { p$HD1ON <- hd[[1]][1]; p$HD1OFF <- hd[[1]][2] }
  if (length(hd) >= 2) { p$HD2ON <- hd[[2]][1]; p$HD2OFF <- hd[[2]][2] }
  if (!is.null(crrt))  { p$CRON  <- crrt[1];    p$CROFF  <- crrt[2] }
  if (!is.null(bic))   { p$BICON <- bic$start; p$BICOFF <- bic$stop
                         p$BICRATE <- bic$rate }
  if (!is.null(fol))   { p$FOLON <- fol[1]; p$FOLOFF <- fol[2] }
  if (!is.null(thi))   { p$THION <- thi[1]; p$THIOFF <- thi[2] }
  if (!is.null(pyr))   { p$PYRON <- pyr[1]; p$PYROFF <- pyr[2] }
  if (!is.null(ca))    { p$CAION_ON <- ca$start; p$CAION_OF <- ca$stop
                         p$CARATE <- ca$rate }

  ev <- data.frame(time = numeric(0), amt = numeric(0), cmt = character(0))
  if (meoh_gkg > 0) ev <- rbind(ev, data.frame(time = 0,
                        amt = mmol_meoh(meoh_gkg, wt), cmt = "A_GUTM"))
  if (eg_mL > 0)    ev <- rbind(ev, data.frame(time = 0,
                        amt = mmol_eg(eg_mL),          cmt = "A_GUTE"))
  if (etoh_gkg > 0) ev <- rbind(ev, data.frame(time = 0,
                        amt = mmol_etoh(etoh_gkg, wt), cmt = "A_GUTETH"))
  if (!is.na(fom_start))
    ev <- rbind(ev, fom_events(fom_start, tend, hd, fom_boost, wt))
  if (!is.na(eth_start)) {
    ev <- rbind(ev, data.frame(time = eth_start,
                               amt = mmol_etoh(eth_load_gkg, wt),
                               cmt = "A_ETOH1"))
    p$ETHON  <- eth_start
    p$ETHOFF <- if (is.na(eth_stop)) tend else eth_stop
    p$ETHRATE <- eth_rate_mgkgh * wt / 46.068
  }
  ev$evid <- 1; ev$ii <- 0; ev$addl <- 0
  list(label = label, tend = tend, param = p,
       ev = if (nrow(ev)) ev[order(ev$time), ] else ev)
}

tap_scenarios <- function() list(
  ## ------------------------------- methanol -------------------------------
  M1_untreated = scen("M1 · methanol 0.7 g/kg, untreated",
                      tend = 72, meoh_gkg = 0.7),
  M2_fom_early = scen("M2 · methanol 0.7 g/kg, fomepizole at 2 h",
                      meoh_gkg = 0.7, fom_start = 2),
  M3_fom_late  = scen("M3 · methanol 0.7 g/kg, fomepizole at 14 h",
                      meoh_gkg = 0.7, fom_start = 14),
  M4_fom_hd    = scen("M4 · methanol 0.7 g/kg, fomepizole 8 h + HD 9-15 h",
                      meoh_gkg = 0.7, fom_start = 8, hd = list(c(9, 15)),
                      bic = list(start = 8, stop = 26, rate = 30)),
  M5_hd_only   = scen("M5 · methanol 0.7 g/kg, HD alone at 15 h",
                      tend = 72, meoh_gkg = 0.7, hd = list(c(15, 21))),
  M6_coingest  = scen("M6 · methanol 0.7 g/kg + ethanol 0.8 g/kg co-ingested",
                      meoh_gkg = 0.7, etoh_gkg = 0.8),
  M7_ethanol   = scen("M7 · methanol 0.7 g/kg, ethanol antidote from 4 h",
                      meoh_gkg = 0.7, eth_start = 4, eth_stop = 60),
  M8_bicarb    = scen("M8 · methanol 0.7 g/kg, bicarbonate only from 8 h",
                      tend = 72, meoh_gkg = 0.7,
                      bic = list(start = 8, stop = 48, rate = 35)),
  M9_full      = scen("M9 · methanol 0.7 g/kg, full care at 6 h",
                      meoh_gkg = 0.7, fom_start = 6, hd = list(c(7, 13)),
                      bic = list(start = 6, stop = 24, rate = 30),
                      fol = c(6, 48)),
  M10_folinate = scen("M10 · methanol 0.7 g/kg, fomepizole 14 h + folinic acid",
                      meoh_gkg = 0.7, fom_start = 14, fol = c(14, 60)),
  M11_nodose   = scen("M11 · as M4 but fomepizole NOT re-dosed on dialysis",
                      meoh_gkg = 0.7, fom_start = 8, fom_boost = FALSE,
                      hd = list(c(9, 15)),
                      bic = list(start = 8, stop = 26, rate = 30)),
  M12_massive  = scen("M12 · methanol 1.5 g/kg, fomepizole + 2 x HD from 4 h",
                      tend = 120, meoh_gkg = 1.5, fom_start = 4,
                      hd = list(c(5, 13), c(20, 26)),
                      bic = list(start = 4, stop = 30, rate = 40)),
  M13_crrt     = scen("M13 · as M4 but 30 h of CRRT instead of 6 h of IHD",
                      tend = 120, meoh_gkg = 0.7, fom_start = 8,
                      crrt = c(9, 39),
                      bic = list(start = 8, stop = 34, rate = 30)),
  M14_eth_hdup = scen("M14 · ethanol antidote 8 h + HD, rate x2.5 on dialysis",
                      meoh_gkg = 0.7, eth_start = 8, eth_stop = 60,
                      eth_hd_mult = 2.5, hd = list(c(9, 15)),
                      bic = list(start = 8, stop = 26, rate = 30)),
  M15_eth_flat = scen("M15 · as M14 but the ethanol rate is NOT increased",
                      meoh_gkg = 0.7, eth_start = 8, eth_stop = 60,
                      eth_hd_mult = 1.0, hd = list(c(9, 15)),
                      bic = list(start = 8, stop = 26, rate = 30)),
  M16_crrt_nd  = scen("M16 · as M13 but fomepizole NOT re-dosed for the circuit",
                      tend = 120, meoh_gkg = 0.7, fom_start = 8,
                      fom_boost = FALSE, crrt = c(9, 39),
                      bic = list(start = 8, stop = 34, rate = 30)),
  M17_isoprop  = scen("M17 · CONTROL: same ADH flux, metabolite is a KETONE",
                      tend = 72, meoh_gkg = 0.7, noacid = 1),
  ## --------------------------- ethylene glycol ----------------------------
  E1_untreated = scen("E1 · ethylene glycol 90 mL, untreated", eg_mL = 90),
  E2_fom_early = scen("E2 · EG 90 mL, fomepizole at 2 h, no dialysis",
                      eg_mL = 90, fom_start = 2),
  E3_fom_hd    = scen("E3 · EG 90 mL, fomepizole 12 h + HD 13-19 h",
                      eg_mL = 90, fom_start = 12, hd = list(c(13, 19)),
                      bic = list(start = 12, stop = 30, rate = 30)),
  E4_cofactors = scen("E4 · EG 90 mL, thiamine + pyridoxine only",
                      eg_mL = 90, thi = c(2, 60), pyr = c(2, 60)),
  E5_fom_cofac = scen("E5 · EG 90 mL, fomepizole 2 h + thiamine/pyridoxine",
                      eg_mL = 90, fom_start = 2, thi = c(2, 60), pyr = c(2, 60)),
  E6_ckd       = scen("E6 · EG 90 mL in CKD (GFR 30%), fomepizole at 2 h",
                      tend = 120, eg_mL = 90, gfr = 0.30, fom_start = 2),
  E7_ethanol   = scen("E7 · EG 90 mL, ethanol antidote from 3 h + HD",
                      eg_mL = 90, eth_start = 3, eth_stop = 40,
                      hd = list(c(6, 12)),
                      bic = list(start = 3, stop = 24, rate = 30)),
  E8_calcium   = scen("E8 · EG 90 mL, fomepizole 12 h + HD + calcium loading",
                      eg_mL = 90, fom_start = 12, hd = list(c(13, 19)),
                      bic = list(start = 12, stop = 30, rate = 30),
                      ca = list(start = 12, stop = 36, rate = 3))
)

## ===========================================================================
##  RUNNER
## ===========================================================================
#' Simulate one scenario.
#' @param mod the mrgsolve model object
#' @param s   one element of tap_scenarios()
#' @param dt  output interval, h
tap_run <- function(mod, s, dt = 0.05) {
  m <- mod %>% param(s$param)
  if (nrow(s$ev) == 0) {
    out <- m %>% mrgsim(end = s$tend, delta = dt, hmax = 0.05)
  } else {
    ## data_set wants a NUMERIC cmt; resolve the names against the model
    d <- s$ev
    d$cmt <- match(as.character(d$cmt), mrgsolve::cmt(mod))
    stopifnot(!any(is.na(d$cmt)))
    d$ID <- 1
    out <- m %>% data_set(d) %>% mrgsim(end = s$tend, delta = dt, hmax = 0.05)
  }
  d <- as.data.frame(out)
  ## Truncate at the first non-survivable pH.  A model that keeps integrating
  ## past pH 6.6 reports numbers that describe a corpse, not a patient.
  bad <- which(d$pHart < s$param$PH_FATAL %||% 6.60)
  if (length(bad)) d <- d[seq_len(bad[1]), ]
  d
}
`%||%` <- function(a, b) if (is.null(a)) b else a

#' One-line summary of a run, matching the reference implementation.
tap_summary <- function(d) {
  fin <- d[nrow(d), ]
  data.frame(
    t_end   = fin$time,
    pH_min  = min(d$pHart),
    HCO3min = min(d$HCO3o),
    AG_max  = max(d$AG),
    OG_max  = max(d$OG),
    FORMmax = max(d$FORMmM),
    FCNSmax = max(d$FCNS),
    FVITmax = max(d$FVIT),
    GLYCmax = max(d$GLYCmM),
    OXALmax = max(d$OXALuM),
    CAOXkid = max(d$CAOXKID),
    CA_min  = min(d$CAIONo),
    GFR_min = min(d$GFRpct),
    PUT     = fin$PUTo,
    logMAR  = fin$LOGMAR,
    P_death = fin$PDEATH,
    P_blind = fin$PBLIND,
    MEOH_ox = fin$MEOHOX,
    EG_ox   = fin$EGOX
  )
}

#' Run every scenario and return the summary table.
tap_all <- function(mod = tap_model(), dt = 0.1) {
  S <- tap_scenarios()
  do.call(rbind, lapply(names(S), function(k) {
    d <- tap_run(mod, S[[k]], dt = dt)
    cbind(data.frame(scenario = k, label = S[[k]]$label), tap_summary(d))
  }))
}

## ===========================================================================
##  ACCEPTANCE TEST against tap_python_reference.py section 15 (scenario M4)
##
##  The two implementations share only the equations: different language,
##  different integrator, independently written.  Agreement to 3 significant
##  figures is therefore a real check and not a tautology.
## ===========================================================================
TAP_REFERENCE_M4 <- read.csv(text = '
time,MEOH,FORMmM,FCNS,FVIT,FOMug,HCO3o,pHart,PACO2o,NAo,AG,OG,LACTo,ATPo,PUTo,OPTICo,LOGMAR,PDEATH,GFRpct,CAIONo
0.0,0,0,0,0,0,24,7.3996,40,140,12,0.002,1,1,0,0,0,0,100,1.1801
2.0,101.44,4.4341,0.59506,0.1793,0,20.667,7.334,40.058,140,15.333,32.832,1.0722,0.96912,0,0,0,0.00021536,100,1.2037
4.0,87.461,7.8596,1.7933,0.64604,0,17.373,7.315,35.184,140,18.627,28.996,1.4665,0.89084,0,0,0,0.00068072,100,1.2106
8.0,61.944,13.91,4.69,2.2351,0,11.647,7.2589,26.84,140,24.353,22.577,2.688,0.73625,0.0014934,0,0,0.0030532,100,1.2308
12.0,21.212,2.6173,3.0693,2.681,6.7797,27.11,7.4343,41.716,145.97,14.861,8.4944,3.1178,0.7774,0.050571,0.0066322,0.013264,0.0097057,100,1.1676
14.0,11.05,1.059,1.8578,2.4453,13.32,30.667,7.462,44.274,146.69,12.024,6.1615,2.6781,0.84951,0.050168,0.0067897,0.013579,0.012873,100,1.1577
15.0,7.9749,0.66298,1.4168,2.3138,8.8711,31.787,7.469,45.155,146.88,11.091,5.4828,2.422,0.88019,0.049968,0.0067491,0.013498,0.01438,100,1.1551
18.0,8.3372,0.39516,0.65054,1.9469,6.3611,32.778,7.4736,46.072,146.63,9.8542,5.9138,1.7708,0.94151,0.049372,0.0066287,0.013257,0.018801,100,1.1535
21.0,7.9649,0.12812,0.29123,1.6188,18.319,33.22,7.476,46.43,146.39,9.1707,5.8272,1.384,0.97266,0.048783,0.0065105,0.013021,0.023133,100,1.1526
24.0,7.6141,0.043976,0.12414,1.3365,13.023,33.278,7.4762,46.497,146.16,8.8816,5.7171,1.1785,0.98807,0.048201,0.0063943,0.012789,0.02739,100,1.1526
30.0,6.9501,0.013267,0.02293,0.90533,5.4431,32.924,7.4737,46.269,145.68,8.7586,5.4579,1.0341,0.9978,0.047058,0.0061682,0.012336,0.035697,100,1.1535
36.0,6.344,0.0086485,0.0062127,0.61175,11.071,32.431,7.4703,45.925,145.22,8.7937,5.2023,1.0074,0.99943,0.045942,0.0059501,0.0119,0.043738,100,1.1547
48.0,5.2711,0.0067167,0.0028453,0.27869,8.7575,31.484,7.4638,45.26,144.42,8.9329,4.7211,1.0022,0.99976,0.043789,0.0055368,0.011074,0.059062,100,1.157
60.0,4.3575,0.0054965,0.00246,0.12713,7.6982,30.639,7.4577,44.666,143.73,9.0945,4.2729,1.0019,0.99979,0.041737,0.0051521,0.010304,0.07344,100,1.1592
72.0,3.5902,0.0042529,0.0020823,0.058432,13.647,29.89,7.4521,44.14,143.16,9.2663,3.8601,1.0017,0.99982,0.039781,0.0047942,0.0095885,0.086939,100,1.1612
96.0,2.4877,0.0027105,0.0011189,0.012796,16.192,28.636,7.4423,43.258,142.26,9.6201,3.1599,1.0008,0.99991,0.036139,0.0041513,0.0083025,0.11155,100,1.1648
')

#' Compare this model against the Python reference trajectory.
#' @param tol relative tolerance; 3 significant figures = 1e-3
tap_verify <- function(mod = tap_model(), tol = 2e-3) {
  s <- tap_scenarios()$M4_fom_hd
  d <- tap_run(mod, s, dt = 0.02)
  ref <- TAP_REFERENCE_M4
  cols <- setdiff(names(ref), "time")
  res <- do.call(rbind, lapply(seq_len(nrow(ref)), function(i) {
    j <- which.min(abs(d$time - ref$time[i]))
    data.frame(time = ref$time[i], variable = cols,
               ref = as.numeric(ref[i, cols]),
               got = as.numeric(d[j, cols]))
  }))
  res$rel <- abs(res$got - res$ref) / pmax(abs(res$ref), 1e-3)
  res$ok  <- res$rel < tol
  cat(sprintf("acceptance test: %d/%d comparisons within %.2g relative\n",
              sum(res$ok), nrow(res), tol))
  bad <- res[!res$ok, ]
  if (nrow(bad)) {
    cat("worst mismatches:\n")
    print(utils::head(bad[order(-bad$rel), ], 15), row.names = FALSE)
  } else {
    cat("PASS — the R and Python implementations agree everywhere.\n")
  }
  invisible(res)
}

## ===========================================================================
##  DERIVATION HELPERS — the arithmetic that the map claims, computed here
## ===========================================================================

#' How complete is each antidote, as a function of substrate concentration?
#' Reproduces section 1 of the Python reference.  A competitive block is
#' DILUTED BY ITS OWN SUBSTRATE, which is why the slogan is wrong.
tap_blockade_table <- function(mod = tap_model()) {
  p <- as.list(param(mod))
  grid <- expand.grid(meoh_mgdl = c(20, 100, 400),
                      etoh_mgdl = c(0, 100, 150),
                      fom_ugml  = c(0, 6, 10))
  grid <- subset(grid, !(etoh_mgdl > 0 & fom_ugml > 0))
  grid$M   <- grid$meoh_mgdl * 10 / 32.042
  grid$E   <- grid$etoh_mgdl * 10 / 46.068
  grid$F   <- grid$fom_ugml / 82.104
  grid$den <- 1 + grid$M / p$KM_M + grid$E / p$KM_ET + grid$F / p$KI_FOM
  grid$base<- 1 + grid$M / p$KM_M
  grid$inh_factor  <- grid$den / grid$base
  grid$pct_flux_left <- 100 * grid$base / grid$den
  grid$flux_mmol_h   <- p$VMAX_M * (grid$M / p$KM_M) / grid$den
  grid[order(grid$meoh_mgdl, grid$etoh_mgdl, grid$fom_ugml),
       c("meoh_mgdl", "etoh_mgdl", "fom_ugml", "inh_factor",
         "pct_flux_left", "flux_mmol_h")]
}

#' Ion trapping: why pH, not concentration, is the prognosticator.
tap_trapping_table <- function(pKa = 3.75, slope = 0.50) {
  ph <- c(7.45, 7.40, 7.30, 7.20, 7.10, 7.00, 6.90, 6.80, 6.70)
  f  <- function(x) 1 / (1 + 10 ^ (x - pKa))
  phb <- 7.00 + slope * (ph - 7.40)
  data.frame(pH = ph, fHA_plasma = f(ph), rel_entry = f(ph) / f(7.40),
             fHA_brain = f(phb), eq_CNS_over_plasma = f(ph) / f(phb))
}

#' The margin that decides whether dialysis can un-block an antidote.
tap_margin_table <- function(mod = tap_model()) {
  p <- as.list(param(mod))
  d <- data.frame(
    antidote     = c("fomepizole", "ethanol"),
    therapeutic  = c(10 / 82.104, 100 * 10 / 46.068),
    K_inhibition = c(p$KI_FOM, p$KM_ET),
    CL_HD        = c(p$CL_HD_FOM, p$CL_HD_ET),
    V_central    = c(0.42 * p$WT, 0.35 * p$WT))
  d$margin  <- d$therapeutic / d$K_inhibition
  d$k_HD    <- d$CL_HD / d$V_central
  d$hours_to_spend_margin <- log(d$margin) / d$k_HD
  d
}

#' The two-gap clock, read off an untreated methanol run.
tap_clock <- function(mod = tap_model(),
                      times = c(0.5, 1, 2, 4, 6, 8, 12, 16, 20, 24, 30, 36, 48)) {
  d <- tap_run(mod, tap_scenarios()$M1_untreated, dt = 0.05)
  i <- vapply(times, function(t) which.min(abs(d$time - t)), 1L)
  d[i, c("time", "MEOH", "FORMmM", "OG", "DAG", "OGoverAG", "GAPSUM", "pHart")]
}

## ===========================================================================
##  DEMONSTRATION (source this file, then call these)
## ===========================================================================
if (identical(environment(), globalenv()) && !interactive()) {
  mod <- tap_model()
  cat("\n=== ADH blockade arithmetic ===\n");   print(tap_blockade_table(mod))
  cat("\n=== ion trapping ===\n");              print(tap_trapping_table())
  cat("\n=== antidote margins on dialysis ===\n"); print(tap_margin_table(mod))
  cat("\n=== the two-gap clock ===\n");         print(tap_clock(mod))
  cat("\n=== all scenarios ===\n");             print(tap_all(mod))
  cat("\n=== acceptance test ===\n");           tap_verify(mod)
}
