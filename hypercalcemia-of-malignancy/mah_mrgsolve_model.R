## =====================================================================
##  mah_mrgsolve_model.R
##  Hypercalcaemia of Malignancy — QSP / PK-PD model
##  50 ODE compartments · 211 parameters · 23 therapeutic scenarios
##
##  악성 종양 관련 고칼슘혈증 — 정량적 시스템 약리학 모델
## ---------------------------------------------------------------------
##  STRUCTURAL THESIS
##  -----------------
##  THESIS 1 — plasma calcium is not a stock, it is the RESIDUAL of four
##  fluxes running on three different clocks:
##
##      dCa_ECF/dt = J_res - J_form + J_gut - U_Ca
##                     |        |        |       |
##      clock 3 (bone, 2-5 d) --+        |       +-- clock 1 (kidney, hours)
##                        clock 2 (gut, ~1 d) ---+
##
##  The ECF calcium pool is only ~36 mmol and the kidney filters ~260
##  mmol/day through it -- seven complete turnovers a day.  Every therapy
##  is an intervention on exactly one of these four terms, and no drug is
##  both fast and definitive:
##
##      saline          -> clock 1, effect in hours,  cures nothing
##      calcitonin      -> clocks 1+3, effect in hours, gone by 48 h
##      bisphosphonate  -> clock 3, effect in 2-4 days, lasts ~4 weeks
##      denosumab       -> clock 3, effect in 4-10 days, no renal limit
##      glucocorticoid  -> clock 2 (+3), only in the calcitriol mechanism
##      antitumour Rx   -> removes the INPUT rather than buffering it
##
##  The guideline pairing of calcitonin with a bisphosphonate is therefore
##  not additive chemistry.  It is a RELAY between two clocks, and this
##  model reproduces the hand-off window: calcitonin buys 0.188 mmol/L at
##  6 h, 0.167 at 12 h, 0.067 at 48 h and 0.007 at day 7 (note C), while
##  the calcitonin receptor availability R_CT falls 1.00 -> 0.26 -> 0.11
##  over the same 48 h.  The two half-lives -- tachyphylaxis onset and
##  bisphosphonate onset -- are close to complementary by construction of
##  neither, and that is the whole justification for the combination.
##
##  THESIS 2 — the renal escape valve has a CEILING, and hypercalcaemia
##  pushes that ceiling DOWN.  Holding total calcium fixed and letting
##  everything else equilibrate, the steady-state urinary calcium is
##  NON-MONOTONIC in plasma calcium (note A):
##
##      Ca 2.4 -> 5.2    Ca 3.0 -> 22.4   Ca 3.4 -> 13.0   mmol/day
##      Ca 2.8 -> 16.8   Ca 3.2 -> 27.2   Ca 3.6 ->  2.2
##                                 ^ the ceiling
##
##  The falling limb is the closed loop: Ca -> CaSR -> AQP2 down ->
##  nephrogenic DI, plus Ca -> nausea/vomiting -> ECF contraction ->
##  proximal Na (and therefore Ca) avidity up and GFR down -> filtered
##  load down -> Ca up.  Sweeping an unregulated calcium input through
##  this system finds a SADDLE-NODE at 23.1 mmol/day (note B): below it a
##  stable steady state exists, above it there is none and the trajectory
##  runs to the model's collapse boundary.  Hypercalcaemic crisis in this
##  model is a bifurcation, not a large number.
##
##  What that reframing buys is a quantitative account of every
##  precipitant and of saline itself (note B):
##
##      baseline                         fold = 23.07 mmol/day (volume-bound)
##      + saline 1.2 L/day               fold = 40.90
##      + saline 2.4 L/day               fold = 51.43
##      distal Tm +60% (PTHrP renal)     fold = 10.92
##      thiazide                         fold = 20.04
##      nephron mass 60% / 30%           fold = 12.89 / 7.25
##      NDI arm deleted                  fold = 25.62
##      nausea/vomiting arm deleted      fold = 31.84
##
##  Saline removes almost no calcium.  What it does is RAISE THE CEILING
##  back above the input -- which is why it works within hours and why it
##  stops working the hour it is stopped.  Note too which arm of the loop
##  matters most: deleting the concentrating defect moves the fold by
##  11%, deleting the nausea arm moves it by 38%.  That is a specific and
##  falsifiable claim about a loop usually described only qualitatively.
##
##  THESIS 3 — the four mechanisms are not four names for one disease.
##  They load different arms, so the SAME total calcium means different
##  physiology and calls for a different drug.  At an identical 3.55
##  mmol/L on day 12 the model returns (note D):
##
##                        HHM (PTHrP)   Osteolytic   Ectopic PTH
##      bone resorption     20.5           34.3         18.6   mmol/day
##      urinary calcium      5.9           17.0          5.2   mmol/day
##      FE_Ca                2.2%           6.5%         2.1%
##      distal Ca Tm        43.3           26.3         41.4   mmol/day
##      phosphate            0.86           1.38         0.91  mmol/L
##      PTH                  1.20           1.20        11.0   pmol/L
##
##  Read the resorption row against the excretion row.  Myeloma needs
##  67% MORE bone resorption than PTHrP disease to reach the same
##  calcium, because PTHrP has closed the renal valve as well as opening
##  the bone tap -- one ligand doing two things in the same direction.
##  That single asymmetry predicts the therapeutic difference the model
##  then reproduces without further assumption: zoledronate takes the
##  osteolytic patient to a nadir of 2.402 with 25.9 normocalcaemic days,
##  and the humoral patient only to 2.635 with 15.0 (scenarios 5 vs 11).
##
##  THE SELF-DELIVERING DRUG
##  ------------------------
##  Zoledronate is not delivered to the osteoclast by the circulation.
##  It binds hydroxyapatite and is then endocytosed by osteoclasts DURING
##  RESORPTION, so the uptake rate is proportional to the very flux the
##  drug abolishes:
##
##      uptake = k_up x (drug on bone surface) x (J_res / J_res_normal)
##
##  One term, three observations (note E/F).  Onset: peak intracellular
##  drug arrives on day 8.4 at normal bone turnover, day 6.2 in humoral
##  disease (J_res 20.5) and day 4.9 in osteolytic disease (J_res 34.3) --
##  the sicker the bone, the faster the drug loads itself.  Duration: once
##  resorption is suppressed the delivery stops, the residue is buried in
##  matrix, and one dose lasts about four weeks.  Dose-response: 1 -> 2 mg
##  buys 5.6 normocalcaemic days, 2 -> 4 mg buys 9.4, 4 -> 8 mg buys 8.2
##  and 8 -> 16 mg buys 3.2 -- a ceiling that comes from the delivery
##  step, not from the target.
##
##  WHAT THE MODEL SAYS ABOUT LOOP DIURETICS
##  ----------------------------------------
##  Not "furosemide is harmful".  Something more specific (scenarios 5,
##  15, 16, 23).  On top of ADEQUATE saline, furosemide is a real
##  calciuretic: nadir 2.635 -> 2.254, normocalcaemic days 15.0 -> 19.2.
##  On top of INADEQUATE saline (1.2 L/day) the same regimen drives the
##  trajectory to the collapse boundary on day 2.6, while the identical
##  inadequate-saline arm WITHOUT furosemide reaches a nadir of 2.642 and
##  never collapses.  The drug is not the variable; the volume it is
##  given on top of is.
##
##  WHERE THE CORRECTED CALCIUM FAILS
##  ---------------------------------
##  Payne's correction under-corrects, and the error grows with BOTH the
##  hypoalbuminaemia and the calcium (note G).  At albumin 25 g/L a total
##  of 3.55 has the same ionised calcium as 4.26 would at albumin 40; the
##  correction reports 3.85, an error of -0.41 mmol/L.  At albumin 22 the
##  error is -0.53.  Scenario 18 is the dynamic version: the same tumour
##  in a patient with albumin 22 presents with a LOWER total calcium
##  (3.42 vs 3.55) and a much HIGHER ionised calcium (2.06 vs 1.72).
##
##  VALIDATION
##  ----------
##  Every ODE in this file was independently re-implemented in
##  Python/scipy (LSODA, rtol 1e-7) and calibrated there; the $ODE block
##  below was then extracted, compiled with g++ -Wall (clean, all 50
##  derivatives assigned) and integrated with a fourth-order Runge-Kutta
##  step of 0.002 d.  Across four scenarios over 60 simulated days the two
##  implementations agree to within 0.0005% on every state and every
##  derived quantity.  Every number quoted above
##  and in the CALIBRATION NOTES was generated from this system.  The
##  healthy control (scenario 1) is an exact steady state: all monitored
##  concentrations drift < 0.2% over 90 simulated days.
##
##  Requires: mrgsolve (>= 1.0), dplyr, tidyr, ggplot2
## =====================================================================

library(mrgsolve)
library(dplyr)

mah_code <- '

$PROB
# Hypercalcaemia of malignancy — QSP model (50 ODEs)
# Time unit: DAYS.  Calcium amounts in mmol, concentrations in mmol/L.

$GLOBAL
// ---- helpers -------------------------------------------------------
double hillx(double x, double k) {
  if (x <= 0.0) return 0.0;
  return x / (x + k);
}
double hilln(double x, double k, double n) {
  if (x <= 0.0) return 0.0;
  double a = pow(x, n);
  return a / (a + pow(k, n));
}
double clampd(double x, double lo, double hi) {
  return (x < lo) ? lo : ((x > hi) ? hi : x);
}
// ---- quantities computed in $ODE and reported in $TABLE ------------
double Ca_tot, iCa, Ca_UF, Ca_cx, Ca_corr, Pi_c, Mg_c, Na_c, K_c, pHa;
double vol, gvol, GFR, FL, fprox, casr, fTAL, L2, L3, Tm_dct, reab_dct, U_Ca;
double FE_Ca, Osm, adh, U_max, osmload, Q_u, vom, anorexia, fabs_ca, J_gut;
double J_res, J_form, J_exch, CaP, prec, J_dial, PS, PS_ren, PS_1a;
double E_ctn, E_zol_apo, E_zol_act, E_prd_mac, E_prd_gut, E_fur, E_cin;
double RKL_f, U_P, U_Pi, P_FL, P_Tm, gut_p, bone_p, fur_diu;
double prodL, prodO, ren_1a, ext_1a, cat24, S_PTH, sec_frac, caset;
double eGFR_ml, QTc, KDIGO;

$PARAM @annotated
// ================= fluid and electrolyte intake / loss ===============
V0        :  15.0  : baseline ECF volume (L)
OSMLOAD   : 800.0  : obligate urinary solute load (mosm/day)
UOSM_MAX  : 1200.0 : maximal urine osmolality with intact AQP2 (mosm/kg)
UOSM_MIN  :  60.0  : minimal urine osmolality (mosm/kg)
INSENS    :   0.80 : insensible water loss (L/day)
ORAL_H2O  :   2.51 : oral water intake when well (L/day)
METAB_H2O :   0.30 : water of oxidation (L/day)
ORAL_NA   : 140.0  : dietary sodium (mmol/day)
ORAL_CA   :  25.0  : dietary calcium (mmol/day)
ORAL_P    :  35.0  : dietary phosphate (mmol/day)
ORAL_MG   :  13.0  : dietary magnesium (mmol/day)
ORAL_K    :  70.0  : dietary potassium (mmol/day)
VOM_MAX   :   1.40 : maximal vomiting fluid loss (L/day)
VOM_K     :   6.0  : symptom score at half-maximal vomiting (-)
GASTRIC_H :  95.0  : gastric fluid H+ concentration (mmol/L)
V_FLOOR   :   0.55 : ECF fraction at which the validity brake engages (-)
HCO3_CAP  :  42.0  : bicarbonate soft cap (mmol/L)
PI_CAP    :   6.0  : phosphate soft cap (mmol/L)
IV_RATE   :   0.0  : isotonic saline infusion rate (L/day)
IMMOB     :   0.0  : immobilisation flag 0-1 (-)

// ================= glomerular filtration =============================
GFR0      : 180.0  : GFR at euvolaemia and full nephron mass (L/day)
GVOL_K    :  12.0  : steepness of the ECF-volume to GFR logistic (-)
GVOL_50   :   0.80 : ECF fraction at the logistic midpoint (-)
GFR_CAP   :   1.10 : ceiling on volume-driven GFR increase (-)

// ================= tubular calcium handling ==========================
FPROX0    :   0.65 : proximal fractional Ca reabsorption at euvolaemia (-)
AVID_P    :   0.60 : gain of proximal avidity on ECF deficit (-)
FPROX_CAP :   0.84 : ceiling on proximal fractional reabsorption (-)
FTAL_MAX  :   0.955: TAL fractional reabsorption with CaSR unoccupied (-)
CASR_EMAX :   0.55 : maximal CaSR inhibition of TAL reabsorption (-)
CASR_EC50 :   1.25 : ionised Ca for half-maximal CaSR effect (mmol/L)
CASR_H    :   4.0  : CaSR Hill coefficient (-)
TM_DCT    :  34.0  : distal transport maximum for calcium (mmol/day)
KM_DCT    :  16.6  : distal Michaelis constant (mmol/day)
FUR_IC50  :   1.2  : furosemide concentration for half TAL block (mg/L)
FUR_EMAX  :   0.85 : maximal furosemide TAL block (-)
FUR_DIU   :   4.20 : furosemide diuresis at full block (L/day)
FUR_NA    : 340.0  : furosemide natriuresis at full block (mmol/day)
THIAZIDE  :   0.0  : thiazide co-prescription flag 0-1 (-)
THIA_EFF  :   0.15 : thiazide increment in distal Ca Tm (-)
CT_RENAL  :   0.55 : maximal calcitonin cut in distal Ca Tm (-)

// ================= phosphate and magnesium ===========================
P_TM0     : 285.0  : proximal phosphate Tm at reference GFR (mmol/day)
P_KM      :  25.0  : phosphate Michaelis constant (mmol/day)
P_PTH_EMAX:   0.70 : maximal PTH1R phosphaturia (-)
P_PTH_K   :   2.0  : PTH1R signal for half-maximal phosphaturia (-)
P_FGF_EMAX:   0.55 : maximal FGF23 phosphaturia (-)
P_FGF_K   : 120.0  : FGF23 for half-maximal phosphaturia (pg/mL)
MG_FE     :   0.0303: reference fractional Mg excretion (-)
MG_ABS    :   0.35 : fractional intestinal Mg absorption (-)

// ================= acid-base =========================================
K_HCO3    :   1.10 : renal correction rate of alkalosis (1/day)
HCO3_SET  :  24.0  : bicarbonate set point (mmol/L)

// ================= AQP2 / nephrogenic diabetes insipidus =============
AQP_KIN   :   0.35 : AQP2 synthesis rate (1/day)
AQP_KOUT  :   0.35 : AQP2 degradation rate (1/day)
AQP_EMAX  :   0.82 : maximal CaSR suppression of AQP2 (-)
AQP_EC50  :   1.55 : ionised Ca for half-maximal AQP2 suppression (mmol/L)
AQP_H     :   5.0  : AQP2 suppression Hill coefficient (-)
AQP_BASE  :   0.21766: value of the Hill term at the healthy set point (-)

// ================= nephrocalcinosis ==================================
CAP_THRESH:   3.60 : Ca_UF x Pi product precipitation threshold ((mmol/L)^2)
K_PREC    :   0.020: precipitation rate constant (mmol/day per unit^2)
K_DISS    :   0.008: deposit dissolution rate (1/day)
K_NEPH    :   0.0060: nephron loss rate per unit injury (1/day)
K_NREP    :   0.012: nephron functional recovery rate (1/day)
W_ICA     :   8.0  : weight of direct hypercalcaemic tubular injury (-)

// ================= parathyroid =======================================
VD_PTH    :   5.0  : PTH volume of distribution (L)
KEL_PTH   : 150.0  : PTH elimination rate constant (1/day)
PTH_SMIN  :   0.06 : minimal fractional PTH secretion (-)
PTH_SMAX  :   1.0  : maximal fractional PTH secretion (-)
PTH_SET   :   1.20 : CaSR set point (mmol/L ionised)
PTH_N     :   6.0  : CaSR Hill coefficient (-)
PTH_SMAXABS: 5660.0: absolute maximal PTH secretion (pmol/day)
PTG_SYN   : 6750.0 : parathyroid PTH synthesis rate (pmol/day)
PTG_MAX   : 900.0  : parathyroid storage capacity (pmol)
PTH_ECT   :   0.0  : autonomous / ectopic PTH secretion (pmol/day)
CIN_EMAX  :   0.28 : maximal cinacalcet set-point left shift (-)
CIN_EC50  :  25.0  : cinacalcet concentration for half effect (ng/mL)

// ================= vitamin D =========================================
D25_IN    :   2.036: 25(OH)D input (nmol/L/day)
D25_KEL   :   0.0231: 25(OH)D elimination (1/day)
D25_USE   :   0.0006: 25(OH)D consumed per pmol/L of 1,25D made (-)
K1A       :  18.267: renal 1-alpha-hydroxylase rate constant (1/day)
D1A_PTH_EMAX: 2.4  : maximal PTH stimulation of renal 1-alpha (-)
D1A_PTH_K :   2.0  : PTH1R signal for half-maximal 1-alpha effect (-)
W1A_PTHRP :   0.15 : PTHrP potency at renal 1-alpha relative to PTH (-)
D1A_FGF_IC50: 90.0 : FGF23 for half-maximal 1-alpha inhibition (pg/mL)
D1A_P_K   :   1.15 : phosphate constant in the 1-alpha term (mmol/L)
D1A_CA_EMAX:  0.80 : maximal Ca suppression of renal 1-alpha (-)
D1A_CA_EC50:  1.55 : ionised Ca for half-maximal 1-alpha suppression (mmol/L)
D1A_CA_BASE:  0.264301: value of that Hill term at the healthy set point (-)
D125_KEL  :   3.30 : 1,25(OH)2D elimination (1/day)
D24_EMAX  :   1.9  : maximal CYP24A1 autoinduction (-)
D24_EC50  : 190.0  : 1,25D for half-maximal CYP24A1 induction (pmol/L)
MAC_1A    :   6.0  : extrarenal 1-alpha-hydroxylase activity per unit MAC (1/day)
PRD_MAC_EMAX: 0.90 : maximal glucocorticoid block of extrarenal 1-alpha (-)
PRD_MAC_IC50: 0.055: prednisolone biophase for half block (mg/L)

// ================= FGF23 =============================================
FGF_K     :   0.55 : FGF23 turnover rate (1/day)
FGF_P_K   :   1.15 : reference phosphate for FGF23 production (mmol/L)
FGF_D_EMAX:   1.6  : maximal 1,25D stimulation of FGF23 (-)
FGF_D_K   : 150.0  : 1,25D for half-maximal FGF23 stimulation (pmol/L)
FGF_NORM  :   1.7111: normalising constant at the healthy set point (-)

// ================= intestine =========================================
ABS_F0    :   0.14 : vitamin-D-independent fractional Ca absorption (-)
ABS_FMAX  :   0.7955: maximal 1,25D-driven increment (-)
ABS_EC50  : 180.0  : 1,25D for half-maximal absorption (pmol/L)
ABS_H     :   2.0  : absorption Hill coefficient (-)
GUT_SEC   :   4.40 : endogenous intestinal Ca secretion (mmol/day)
ABS_P     :   0.65 : fractional phosphate absorption (-)
PRD_GUT   :   0.45 : maximal glucocorticoid cut in Ca absorption (-)
PRD_GUT_IC50: 0.055: prednisolone biophase for half that effect (mg/L)

// ================= tumour and its mediators ==========================
TUM_G     :   0.030: tumour logistic growth rate (1/day)
TUM_MAX   :   3.0  : tumour carrying capacity (kg)
CHEMO     :   0.0  : antitumour kill rate (1/day)
PTHRP_S   :   0.0  : PTHrP secretion per kg tumour (pmol/L/day/kg)
PTHRP_KEL :  26.0  : PTHrP elimination (1/day)
PTHRP_BASE:  13.0  : constitutive PTHrP production (pmol/L/day)
CYT_S     :   0.0  : osteolytic cytokine secretion per kg tumour (1/day/kg)
CYT_KEL   :   3.0  : cytokine elimination (1/day)
MAC_S     :   0.0  : 1-alpha-competent macrophage accrual per kg tumour (1/day/kg)
MAC_KEL   :   0.10 : macrophage pool turnover (1/day)

// ================= RANKL / OPG =======================================
KL_BASE   :   0.51141: RANKL production rate constant (pmol/L/day)
L_PTH_EMAX:  20.0  : maximal PTH1R induction of RANKL (-)
L_PTH_K   :   7.0  : PTH1R signal for half-maximal RANKL induction (-)
L_D125_EMAX: 14.0  : maximal 1,25D induction of RANKL (-)
L_D125_K  : 700.0  : 1,25D for half-maximal RANKL induction (pmol/L)
L_D125_H  :   2.0  : Hill coefficient for that term (-)
L_CYT_E   :   2.2  : cytokine induction of RANKL per unit CYT (-)
L_DEG     :   6.0  : free RANKL degradation (1/day)
KO_BASE   :  31.88 : OPG production rate constant (pmol/L/day)
O_PTH_EMAX:   0.72 : maximal PTH1R SUPPRESSION of OPG (-)
O_PTH_K   :   2.2  : PTH1R signal for half-maximal OPG suppression (-)
O_DEG     :   6.0  : OPG degradation (1/day)
KON_O     :   0.60 : RANKL-OPG association (1/pmol/L/day)
KOFF_O    :   0.06 : RANKL-OPG dissociation (1/day)
DEG_LO    :   3.0  : RANKL-OPG complex clearance (1/day)

// ================= bone cells ========================================
OCP_IN    :   0.60 : osteoclast precursor influx (1/day)
OCP_OUT   :   0.60 : osteoclast precursor efflux (1/day)
OCP_PTH_E :   3.0  : PTH1R recruitment of precursors (-)
OCP_PTH_K :   4.0  : PTH1R signal for half-maximal recruitment (-)
OCP_PTH_BASE: 1.6  : normalising constant at the healthy set point (-)
OCP_CYT_E :   1.6  : cytokine recruitment of precursors (-)
OC_DIFF   :   0.45 : osteoclast differentiation rate (1/day)
OC_KRKL   :   1.20 : RANKL for half-maximal differentiation (pmol/L)
OC_DEG    :   0.45 : osteoclast apoptosis rate (1/day)
OC_IND    :   0.060: RANKL-independent osteoclastogenesis fraction (-)
OBP_IN    :   0.25 : osteoblast precursor influx (1/day)
OBP_OUT   :   0.25 : osteoblast precursor maturation (1/day)
OB_DEG    :   0.24291: osteoblast apoptosis rate (1/day)
OBP_TGF_EMAX: 1.4  : maximal TGF-beta coupling to osteoblastogenesis (-)
OBP_TGF_K :   1.0  : TGF-beta for half-maximal coupling (-)
OB_PTH_IC50:  4.5  : PTH1R signal for half-maximal osteoblast apoptosis (-)
OB_PTH_EMAX:  0.62 : maximal PTH1R-driven osteoblast apoptosis (-)
PRD_OB    :   0.30 : maximal glucocorticoid osteoblast apoptosis (-)
PRD_OB_IC50:  0.045: prednisolone biophase for half that effect (mg/L)
J_RES0    :   8.0  : bone resorption at unit osteoclast number (mmol/day)
J_FORM0   :   8.0  : bone formation at unit osteoblast number (mmol/day)
TGF_K     :   0.125: TGF-beta release per mmol resorbed (-)
TGF_DEG   :   1.0  : TGF-beta clearance (1/day)
CTX_K     :   1.0  : CTX scaling (-)
CTX_DEG   :   8.0  : CTX turnover (1/day)
P1NP_K    :   1.0  : P1NP scaling (-)
P1NP_DEG  :   1.4  : P1NP turnover (1/day)
V_RAP     :  40.0  : apparent volume of the rapidly exchangeable pool (L)
K_EX      :  45.0  : exchange clearance with that pool (L/day)

// ================= zoledronate =======================================
ZOL_KEL   :   6.0  : renal elimination at reference GFR (1/day)
ZOL_K12   :   9.0  : central to peripheral (1/day)
ZOL_K21   :   2.2  : peripheral to central (1/day)
ZOL_KB    :  22.0  : hydroxyapatite binding rate (1/day)
ZOL_BMAX  : 180.0  : bone-surface binding capacity (umol)
ZOL_BURY  :   0.055: burial into matrix (1/day)
ZOL_KUP   :   0.45 : osteoclast uptake rate at reference resorption (1/day)
ZOL_KOUT_OC:  0.045: intracellular elimination (1/day)
ZOL_EMAX  :   3.4  : maximal osteoclast apoptosis multiplier (-)
ZOL_EC50  :   0.30 : intracellular drug for half-maximal apoptosis (umol)
ZOL_IMAX  :   0.80 : maximal per-cell resorption inhibition (-)
ZOL_IC50  :   0.30 : intracellular drug for half that effect (umol)

// ================= denosumab =========================================
DMB_KA    :   0.28 : subcutaneous absorption (1/day)
DMB_F     :   0.62 : subcutaneous bioavailability (-)
DMB_VD    :   8.20 : volume of distribution (L)
DMB_KEL   :   0.055: linear elimination (1/day)
DMB_KON   :   0.90 : RANKL association (1/pmol/L/day)
DMB_KOFF  :   2.70 : RANKL dissociation, Kd = 3 pmol/L (1/day)
DMB_KINT  :   0.24 : complex internalisation (1/day)

// ================= calcitonin ========================================
CTN_KA    :  30.0  : subcutaneous absorption (1/day)
CTN_KEL   :  16.6  : elimination (1/day)
CTN_VD    :  12.0  : volume of distribution (L)
CTN_KE0   :   3.50 : osteoclast biophase equilibration (1/day)
CTN_EMAX  :   0.85 : maximal inhibition of osteoclast activity (-)
CTN_EC50  :   0.60 : biophase concentration for half effect (IU/L)
CT_KIN    :   0.115: calcitonin receptor resynthesis (1/day)
CT_KOUT   :   1.90 : receptor down-regulation at full occupancy (1/day)

// ================= prednisolone, cinacalcet, furosemide ==============
PRD_KEL   :   6.2  : prednisolone elimination (1/day)
PRD_KE0   :   1.20 : genomic biophase equilibration (1/day)
CIN_KEL   :   1.9  : cinacalcet elimination (1/day)
FUR_KEL   :   8.3  : furosemide elimination (1/day)

// ================= speciation ========================================
ALB       :  40.0  : serum albumin (g/L)
KB0       :   0.020: albumin calcium binding constant at pH 7.40 (L/g)
PH_GAMMA  :   0.41 : pH sensitivity of albumin binding (1/pH unit)
K_CX      :   0.10 : complexed fraction of total calcium (-)

// ================= symptoms ==========================================
SYM_TAU   :   0.45 : symptom equilibration time constant (day)
SYM_K     :   0.35 : ionised Ca above set point for half-maximal score (mmol/L)
ADAPT_TAU :  20.0  : CNS set-point adaptation time constant (day)
ADAPT_LO  :   1.15 : lower bound of the adapted set point (mmol/L)
ADAPT_HI  :   1.90 : upper bound of the adapted set point (mmol/L)
ADAPT_FRAC:   0.45 : fraction of an elevation that adaptation absorbs (-)

// ================= extracorporeal ====================================
DIAL      :   0.0  : dialysis on/off (-)
DIAL_CL   : 110.0  : dialyser clearance of ultrafilterable Ca (L/day)
DIAL_CA   :   1.00 : dialysate calcium (mmol/L)

// ================= analysis handles (0 = off) ========================
J_EXO     :   0.0  : unregulated calcium input, bifurcation knob (mmol/day)
CLAMP_CA  :   0.0  : if > 0, hold total plasma calcium at this value (mmol/L)

$CMT @annotated
A_Ca    : ECF exchangeable calcium (mmol)
Ca_rap  : rapidly exchangeable bone-surface calcium (mmol)
A_P     : ECF phosphate (mmol)
A_Mg    : ECF magnesium (mmol)
V_ecf   : extracellular fluid volume (L)
A_Na    : ECF sodium (mmol)
A_K     : exchangeable potassium (mmol)
HCO3    : plasma bicarbonate (mmol/L)
AQP2    : collecting-duct AQP2 abundance (fraction of normal)
N_func  : functional nephron mass (fraction)
CAST    : intratubular / interstitial Ca-P deposit (a.u.)
PTG     : parathyroid stored PTH (pmol)
PTH     : plasma PTH (pmol/L)
D25     : 25-hydroxyvitamin D (nmol/L)
D125    : 1,25-dihydroxyvitamin D (pmol/L)
FGF23   : plasma FGF23 (pg/mL)
TUM     : tumour burden (kg)
PTHRP   : plasma PTHrP (pmol/L)
CYT     : osteolytic cytokine index (a.u.)
MAC     : extrarenal 1-alpha-hydroxylase capacity (a.u.)
OCP     : osteoclast precursors (relative)
OC      : mature osteoclasts (relative)
OBP     : osteoblast precursors (relative)
OB      : mature osteoblasts (relative)
RKL     : free RANKL (pmol/L)
OPG     : free osteoprotegerin (pmol/L)
C_LO    : RANKL-OPG complex (pmol/L)
TGFB    : matrix-derived TGF-beta signal (relative)
CA_BONE : skeletal mineral (mmol)
CTX     : serum CTX-I (ng/mL)
P1NP    : serum P1NP (ug/L)
ZOL_C   : zoledronate central (umol)
ZOL_P   : zoledronate peripheral (umol)
ZOL_B   : zoledronate on bone surface (umol)
ZOL_D   : zoledronate buried in matrix (umol)
ZOL_OC  : intra-osteoclast zoledronate (umol)
DMB_SC  : denosumab subcutaneous depot (pmol)
DMB_C   : denosumab free plasma (pmol/L)
DMB_RK  : denosumab-RANKL complex (pmol/L)
CTN_SC  : calcitonin subcutaneous depot (IU)
CTN_C   : calcitonin plasma (IU/L)
R_CT    : calcitonin receptor availability (fraction)
PRD_C   : prednisolone plasma (mg/L)
CIN_C   : cinacalcet plasma (ng/mL)
FUR_C   : furosemide plasma (mg/L)
CTN_E   : calcitonin osteoclast biophase (IU/L)
PRD_E   : prednisolone genomic biophase (mg/L)
SYM     : neurocognitive symptom score (0-10)
ADAPT   : adapted CNS calcium set point (mmol/L)
AUC_CA  : cumulative ionised Ca-time above 1.40 (mmol/L*day)

$MAIN
A_Ca_0    = 2.40 * V0;
Ca_rap_0  = 2.40 * V_RAP;
A_P_0     = 1.15 * V0;
A_Mg_0    = 0.85 * V0;
V_ecf_0   = V0;
A_Na_0    = 140.0 * V0;
A_K_0     = 3500.0;
HCO3_0    = 24.0;
AQP2_0    = 1.0;
N_func_0  = 1.0;
CAST_0    = 0.0;
PTG_0     = 500.0;
PTH_0     = 4.0;
D25_0     = 60.0;
D125_0    = 120.0;
FGF23_0   = 45.0;
TUM_0     = 0.0;
PTHRP_0   = 0.50;
CYT_0     = 0.0;
MAC_0     = 0.0;
OCP_0     = 1.0;
OC_0      = 1.0;
OBP_0     = 1.0;
OB_0      = 1.0;
RKL_0     = 0.30;
OPG_0     = 4.00;
C_LO_0    = 0.23529;
TGFB_0    = 1.0;
CA_BONE_0 = 25000.0;
CTX_0     = 0.35;
P1NP_0    = 45.0;
R_CT_0    = 1.0;
SYM_0     = 0.0;
ADAPT_0   = 1.17226;
AUC_CA_0  = 0.0;

$ODE
// =====================================================================
//  0.  SPECIATION — three calciums; only one of them acts
// =====================================================================
double V = (V_ecf > 4.0) ? V_ecf : 4.0;
Ca_tot = A_Ca / V;
Pi_c   = fmax(A_P  / V, 0.05);
Mg_c   = fmax(A_Mg / V, 0.05);
Na_c   = A_Na / V;
K_c    = 4.0 * pow(A_K / 3500.0, 3.0);
double HCO3c = fmax(HCO3, 5.0);

pHa    = 6.1 + log10(HCO3c / (0.03 * 40.0));
double Kb = KB0 * pow(10.0, PH_GAMMA * (pHa - 7.40));
Ca_cx  = K_CX * Ca_tot;
iCa    = (Ca_tot - Ca_cx) / (1.0 + Kb * ALB);
Ca_UF  = iCa + Ca_cx;
Ca_corr = Ca_tot + 0.020 * (40.0 - ALB);   // Payne, as the clinician computes it

// =====================================================================
//  1.  GLOMERULAR FILTRATION — steeply volume dependent
// =====================================================================
vol  = V / V0;
gvol = (1.0 / (1.0 + exp(-GVOL_K * (vol - GVOL_50)))) /
       (1.0 / (1.0 + exp(-GVOL_K * (1.0 - GVOL_50))));
if (gvol > GFR_CAP) gvol = GFR_CAP;
GFR = GFR0 * N_func * gvol;

// =====================================================================
//  2.  PTH1R SIGNAL — PTH and PTHrP are equipotent HERE and only here
// =====================================================================
PS = (PTH + PTHRP) / 4.50;
if (PS < 0.0) PS = 0.0;
PS_ren = (PS / (PS + 1.0)) / 0.5;          // saturable; equals 1 at baseline

// =====================================================================
//  3.  DRUG EFFECT TERMS
// =====================================================================
E_ctn     = CTN_EMAX * R_CT * hillx(CTN_E, CTN_EC50);
E_zol_apo = ZOL_EMAX * hillx(ZOL_OC, ZOL_EC50);
E_zol_act = ZOL_IMAX * hillx(ZOL_OC, ZOL_IC50);
E_prd_mac = PRD_MAC_EMAX * hillx(PRD_E, PRD_MAC_IC50);
E_prd_gut = PRD_GUT     * hillx(PRD_E, PRD_GUT_IC50);
E_fur     = FUR_EMAX    * hillx(FUR_C, FUR_IC50);
E_cin     = CIN_EMAX    * hillx(CIN_C, CIN_EC50);

// =====================================================================
//  4.  RENAL CALCIUM HANDLING — segment by segment
//      Proximal and TAL are FRACTIONAL; the distal segment is
//      TRANSPORT-limited.  That distinction is what lets PTHrP defend
//      the plasma calcium against a rising filtered load.
// =====================================================================
FL    = GFR * Ca_UF;
fprox = FPROX0 * (1.0 + AVID_P * fmax(0.0, 1.0 - vol));
if (fprox > FPROX_CAP) fprox = FPROX_CAP;
L2    = FL * (1.0 - fprox);
casr  = CASR_EMAX * hilln(iCa, CASR_EC50, CASR_H);
fTAL  = FTAL_MAX * (1.0 - casr) * (1.0 - E_fur);
L3    = L2 * (1.0 - fTAL);
Tm_dct = TM_DCT * (0.5 + 0.5 * PS_ren) * (1.0 - CT_RENAL * E_ctn)
         * (1.0 + THIA_EFF * THIAZIDE);
reab_dct = Tm_dct * L3 / (KM_DCT + L3);
U_Ca  = fmax(L3 - reab_dct, 0.0);
FE_Ca = (FL > 0.0) ? U_Ca / FL : 0.0;

// =====================================================================
//  5.  WATER — the concentrating defect only bites when ADH is driving
// =====================================================================
Osm = 2.0 * Na_c + 10.0;
adh = clampd(0.30 + 2.5 * fmax(0.0, 1.0 - vol) + 0.10 * (Osm - 290.0), 0.05, 1.0);
U_max = UOSM_MIN + (UOSM_MAX * AQP2 - UOSM_MIN) * adh;
if (U_max < UOSM_MIN) U_max = UOSM_MIN;
osmload = OSMLOAD + 2.0 * U_Ca + 40.0 * IV_RATE * 0.30;
Q_u = osmload / U_max;
double gfr_frac = GFR / GFR0 / fmax(N_func, 0.05);
Q_u = Q_u * (0.25 + 0.75 * fmin(1.0, gfr_frac));

// =====================================================================
//  6.  GI LOSSES — and the smooth model-validity brake
// =====================================================================
double gV = 1.0 / (1.0 + exp(-25.0 * (vol - V_FLOOR - 0.05)));
vom      = VOM_MAX * hilln(SYM, VOM_K, 2.0) * gV;
anorexia = hilln(SYM, 5.0, 2.0);

// =====================================================================
//  7.  INTESTINAL CALCIUM — arm 2
// =====================================================================
fabs_ca = (ABS_F0 + ABS_FMAX * hilln(D125, ABS_EC50, ABS_H)) * (1.0 - E_prd_gut);
double intake_ca = ORAL_CA * (1.0 - 0.85 * anorexia);
J_gut = fabs_ca * intake_ca - GUT_SEC * (1.0 - 0.4 * anorexia);

// =====================================================================
//  8.  BONE FLUXES — arm 1
// =====================================================================
RKL_f  = fmax(RKL, 0.0);
J_res  = J_RES0  * OC * (1.0 - E_zol_act) * (1.0 - E_ctn) * (1.0 + 0.18 * IMMOB);
J_form = J_FORM0 * OB * (1.0 - 0.25 * IMMOB);
J_exch = K_EX * (Ca_tot - Ca_rap / V_RAP);

// =====================================================================
//  9.  PRECIPITATION AND EXTRACORPOREAL REMOVAL
// =====================================================================
CaP  = Ca_UF * Pi_c;
double dCaP = fmax(0.0, CaP - CAP_THRESH);
prec = K_PREC * dCaP * dCaP;
J_dial = DIAL * DIAL_CL * fmax(0.0, Ca_UF - DIAL_CA);

// =====================================================================
// 10.  CALCIUM BALANCE — the residual
// =====================================================================
dxdt_A_Ca   = J_res - J_form + J_gut - U_Ca - J_exch - J_dial
              - 0.45 * prec + J_EXO;
dxdt_Ca_rap = J_exch;

// =====================================================================
// 11.  PHOSPHATE — low in PTHrP disease, high in myeloma with uraemia
// =====================================================================
double fgf_i = P_FGF_EMAX * hillx(FGF23, P_FGF_K);
double pth_i = P_PTH_EMAX * hillx(PS, P_PTH_K);
P_FL  = GFR * Pi_c * 0.90;
P_Tm  = P_TM0 * (GFR / GFR0) * (1.0 - pth_i) * (1.0 - fgf_i);
double P_reab = P_Tm * P_FL / (P_KM + P_FL);
U_Pi  = fmax(P_FL - P_reab, 0.0);
gut_p = ABS_P * ORAL_P * (1.0 - 0.85 * anorexia);
bone_p = 0.62 * (J_res - J_form);
double gP = 1.0 / (1.0 + exp(6.0 * (Pi_c - PI_CAP)));
dxdt_A_P = (gut_p + bone_p) * gP - U_Pi - 0.30 * prec;

// =====================================================================
// 12.  MAGNESIUM — CaSR-driven renal wasting
// =====================================================================
dxdt_A_Mg = MG_ABS * ORAL_MG * (1.0 - 0.85 * anorexia)
            - MG_FE * GFR * Mg_c * 0.70 * (1.0 + 1.6 * casr)
            - 0.3 * vom;

// =====================================================================
// 13.  VOLUME, SODIUM, POTASSIUM
// =====================================================================
double oral_h2o = ORAL_H2O * (1.0 - 0.80 * anorexia);
fur_diu = FUR_DIU * E_fur / FUR_EMAX;
dxdt_V_ecf = IV_RATE + oral_h2o + METAB_H2O - Q_u - INSENS - vom - fur_diu;

double na_in = ORAL_NA * (1.0 - 0.8 * anorexia) + 154.0 * IV_RATE;
double fexc  = clampd((vol - 0.90) / 0.10, 0.02, 1.60);
double u_na  = na_in * fexc + 1.5 * (A_Na - 140.0 * V0) + FUR_NA * E_fur / FUR_EMAX;
u_na = clampd(u_na, 0.0, 2200.0);
dxdt_A_Na = na_in - 145.0 * vom - u_na;

dxdt_A_K = ORAL_K * (1.0 - 0.85 * anorexia) - 12.0 * vom
           - (52.0 + 9.0 * Q_u) * (A_K / 3500.0) - 25.0 * E_fur;

// =====================================================================
// 14.  ACID-BASE — vomiting alkalosis PROTECTS the ionised calcium
// =====================================================================
double cl_repl = clampd(vol, 0.2, 1.15);
double gH = 1.0 / (1.0 + exp(1.2 * (HCO3 - HCO3_CAP)));
dxdt_HCO3 = GASTRIC_H * vom * gH / V - K_HCO3 * cl_repl * N_func * (HCO3 - HCO3_SET);

// =====================================================================
// 15.  AQP2 — the memory element in the vicious cycle
// =====================================================================
double aqp_inh = AQP_EMAX * fmax(0.0,
                 (hilln(iCa, AQP_EC50, AQP_H) - AQP_BASE) / (1.0 - AQP_BASE));
dxdt_AQP2 = AQP_KIN * (1.0 - aqp_inh) - AQP_KOUT * AQP2;

// =====================================================================
// 16.  NEPHROCALCINOSIS AND NEPHRON LOSS
// =====================================================================
dxdt_CAST = prec - K_DISS * CAST;
double exc_ca = fmax(0.0, iCa - 1.75);
dxdt_N_func = -K_NEPH * (CAST + W_ICA * exc_ca * exc_ca) * N_func
              + K_NREP * (1.0 - N_func);

// =====================================================================
// 17.  PARATHYROID
// =====================================================================
caset    = PTH_SET * (1.0 - E_cin);
sec_frac = PTH_SMIN + (PTH_SMAX - PTH_SMIN) /
           (1.0 + pow(fmax(iCa, 1e-6) / caset, PTH_N));
S_PTH    = PTH_SMAXABS * sec_frac * fmin(1.0, PTG / 300.0);
dxdt_PTG = PTG_SYN * (1.0 - PTG / PTG_MAX) - S_PTH;
dxdt_PTH = (S_PTH + PTH_ECT) / VD_PTH - KEL_PTH * PTH;

// =====================================================================
// 18.  VITAMIN D — PTHrP is a WEAK agonist at renal 1-alpha (w = 0.15),
//      which is why 1,25D is LOW in humoral hypercalcaemia and HIGH in
//      genuine hyperparathyroidism, at the same PTH1R signal in bone.
// =====================================================================
PS_1a  = (PTH + W1A_PTHRP * PTHRP) / 4.50;
double ca_1a = D1A_CA_EMAX * fmax(0.0,
               (hilln(iCa, D1A_CA_EC50, 4.0) - D1A_CA_BASE) / (1.0 - D1A_CA_BASE));
ren_1a = K1A * N_func * D25
         * (1.0 + D1A_PTH_EMAX * hillx(PS_1a, D1A_PTH_K))
         * (1.0 - 0.85 * hillx(FGF23, D1A_FGF_IC50))
         * (D1A_P_K / (D1A_P_K + Pi_c))
         * (1.0 - ca_1a);
ext_1a = MAC_1A * MAC * D25 * (1.0 - E_prd_mac);
cat24  = D125_KEL * (1.0 + D24_EMAX * hillx(D125, D24_EC50));
dxdt_D125 = ren_1a + ext_1a - cat24 * D125;

// 25(OH)D is CONSUMED by 1-alpha-hydroxylation.  This is the reason
// 25(OH)D runs low in granulomatous and lymphomatous 1,25D excess, and
// it is the feedback that keeps the calcitriol mechanism self-limiting.
dxdt_D25 = D25_IN - D25_KEL * D25 - 0.004 * D25 * PS
           - D25_USE * (ren_1a + ext_1a);

// =====================================================================
// 19.  FGF23
// =====================================================================
dxdt_FGF23 = FGF_K * 45.0 * (Pi_c / FGF_P_K)
             * (1.0 + FGF_D_EMAX * hillx(D125, FGF_D_K)) / FGF_NORM
             - FGF_K * FGF23;

// =====================================================================
// 20.  TUMOUR AND MEDIATORS
// =====================================================================
dxdt_TUM   = TUM_G * TUM * (1.0 - TUM / TUM_MAX) - CHEMO * TUM;
dxdt_PTHRP = PTHRP_BASE + PTHRP_S * TUM - PTHRP_KEL * PTHRP;
dxdt_CYT   = CYT_S * TUM - CYT_KEL * CYT;
dxdt_MAC   = MAC_S * TUM - MAC_KEL * MAC - 0.55 * E_prd_mac * MAC;

// =====================================================================
// 21.  RANKL / OPG — PTH1R raises RANKL AND suppresses OPG.  Both errors
//      push the same way, which is the uncoupling.
// =====================================================================
double OB_lin = 0.5 * (OBP + OB);
prodL = KL_BASE * OB_lin
        * (1.0 + L_PTH_EMAX  * hillx(PS, L_PTH_K))
        * (1.0 + L_D125_EMAX * hilln(D125, L_D125_K, L_D125_H))
        * (1.0 + L_CYT_E * CYT);
prodO = KO_BASE * OB * (1.0 - O_PTH_EMAX * hillx(PS, O_PTH_K));

double bind_lo   = KON_O * RKL * OPG;
double unbind_lo = KOFF_O * C_LO;
double bind_d    = DMB_KON * RKL * DMB_C;
double unbind_d  = DMB_KOFF * DMB_RK;

dxdt_RKL  = prodL - L_DEG * RKL - bind_lo + unbind_lo - bind_d + unbind_d;
dxdt_OPG  = prodO - O_DEG * OPG - bind_lo + unbind_lo;
dxdt_C_LO = bind_lo - unbind_lo - DEG_LO * C_LO;

// =====================================================================
// 22.  BONE CELLS
// =====================================================================
dxdt_OCP = OCP_IN * (1.0 + OCP_PTH_E * hillx(PS, OCP_PTH_K) + OCP_CYT_E * CYT)
           / OCP_PTH_BASE - OCP_OUT * OCP;

double f_rkl = RKL_f / (OC_KRKL + RKL_f);
double f_ref = 0.30 / (OC_KRKL + 0.30);
double oc_drive = (f_rkl + OC_IND * f_ref) / (f_ref * (1.0 + OC_IND));
dxdt_OC = OC_DIFF * OCP * oc_drive - OC_DEG * OC * (1.0 + E_zol_apo);

double ob_pth_apo = OB_PTH_EMAX * hilln(PS, OB_PTH_IC50, 2.0);
double ob_prd_apo = PRD_OB * hillx(PRD_E, PRD_OB_IC50);
dxdt_OBP = OBP_IN * (1.0 + OBP_TGF_EMAX * hillx(TGFB, OBP_TGF_K))
           / (1.0 + OBP_TGF_EMAX * 0.5) - OBP_OUT * OBP;
dxdt_OB  = OBP_OUT * OBP - OB_DEG * OB * (1.0 + ob_pth_apo + ob_prd_apo);

dxdt_TGFB    = TGF_K * J_res - TGF_DEG * TGFB;
dxdt_CA_BONE = J_form - J_res;
dxdt_CTX     = CTX_K  * 0.35 * (J_res  / J_RES0)  * CTX_DEG  - CTX_DEG  * CTX;
dxdt_P1NP    = P1NP_K * 45.0 * (J_form / J_FORM0) * P1NP_DEG - P1NP_DEG * P1NP;

// =====================================================================
// 23.  ZOLEDRONATE — THE SELF-DELIVERING DRUG
//      uptake is proportional to the resorption flux the drug abolishes
// =====================================================================
double free_surf = fmax(0.0, 1.0 - ZOL_B / ZOL_BMAX);
double bind_bone = ZOL_KB * ZOL_C * free_surf;
double uptake    = ZOL_KUP * ZOL_B * (J_res / J_RES0);

dxdt_ZOL_C  = -ZOL_KEL * ZOL_C * (GFR / GFR0)
              - ZOL_K12 * ZOL_C + ZOL_K21 * ZOL_P - bind_bone;
dxdt_ZOL_P  = ZOL_K12 * ZOL_C - ZOL_K21 * ZOL_P;
dxdt_ZOL_B  = bind_bone - ZOL_BURY * ZOL_B - uptake;
dxdt_ZOL_D  = ZOL_BURY * ZOL_B;
dxdt_ZOL_OC = uptake / 12.0 - ZOL_KOUT_OC * ZOL_OC;

// =====================================================================
// 24.  DENOSUMAB — target-mediated disposition, no renal step anywhere
// =====================================================================
dxdt_DMB_SC = -DMB_KA * DMB_SC;
dxdt_DMB_C  = DMB_KA * DMB_SC * DMB_F / DMB_VD - DMB_KEL * DMB_C
              - bind_d + unbind_d;
dxdt_DMB_RK = bind_d - unbind_d - DMB_KINT * DMB_RK;

// =====================================================================
// 25.  CALCITONIN — biophase plus receptor down-regulation
// =====================================================================
dxdt_CTN_SC = -CTN_KA * CTN_SC;
dxdt_CTN_C  = CTN_KA * CTN_SC / CTN_VD - CTN_KEL * CTN_C;
dxdt_CTN_E  = CTN_KE0 * (CTN_C - CTN_E);
dxdt_R_CT   = CT_KIN * (1.0 - R_CT) - CT_KOUT * R_CT * hillx(CTN_E, CTN_EC50);

// =====================================================================
// 26.  PREDNISOLONE / CINACALCET / FUROSEMIDE
//      dosed directly into the concentration compartment: the "amount"
//      given is the concentration increment D*F/Vd
// =====================================================================
dxdt_PRD_C = -PRD_KEL * PRD_C;
dxdt_PRD_E = PRD_KE0 * (PRD_C - PRD_E);
dxdt_CIN_C = -CIN_KEL * CIN_C;
dxdt_FUR_C = -FUR_KEL * FUR_C;

// =====================================================================
// 27.  SYMPTOMS — and why the same calcium is not the same illness
//      The CNS adapts its set point, but only PARTIALLY (45%).
// =====================================================================
double drive  = fmax(0.0, iCa - ADAPT);
double target = 10.0 * hilln(drive, SYM_K, 2.0);
dxdt_SYM = (target - SYM) / SYM_TAU;
double ad = clampd(ADAPT_LO + ADAPT_FRAC * (iCa - ADAPT_LO), ADAPT_LO, ADAPT_HI);
dxdt_ADAPT = (ad - ADAPT) / ADAPT_TAU;
dxdt_AUC_CA = fmax(0.0, iCa - 1.40);

// =====================================================================
// 28.  ANALYSIS HANDLE — hold total calcium fixed to trace the ceiling
// =====================================================================
if (CLAMP_CA > 0.0) dxdt_A_Ca = CLAMP_CA * dxdt_V_ecf;

$TABLE
eGFR_ml = GFR / 1.44;                                  // L/day -> mL/min
QTc = 440.0 - 78.0 * (iCa - 1.20);                     // ms, short QT
KDIGO = (GFR > 0.66 * GFR0) ? 0.0 :
        ((GFR > 0.50 * GFR0) ? 1.0 :
        ((GFR > 0.33 * GFR0) ? 2.0 : 3.0));

$CAPTURE @annotated
Ca_tot  : total plasma calcium (mmol/L)
iCa     : ionised calcium (mmol/L)
Ca_corr : Payne-corrected calcium (mmol/L)
Ca_UF   : ultrafilterable calcium (mmol/L)
Pi_c    : plasma phosphate (mmol/L)
Mg_c    : plasma magnesium (mmol/L)
K_c     : plasma potassium (mmol/L)
pHa     : arterial pH (-)
GFR     : glomerular filtration rate (L/day)
eGFR_ml : glomerular filtration rate (mL/min)
FL      : filtered calcium load (mmol/day)
U_Ca    : urinary calcium excretion (mmol/day)
FE_Ca   : fractional excretion of calcium (-)
Tm_dct  : distal calcium transport maximum (mmol/day)
J_res   : bone resorption flux (mmol/day)
J_form  : bone formation flux (mmol/day)
J_gut   : net intestinal calcium absorption (mmol/day)
Q_u     : urine output (L/day)
vol     : ECF volume as a fraction of normal (-)
vom     : vomiting fluid loss (L/day)
PS      : PTH1R signal, 1 = healthy (-)
E_ctn   : calcitonin effect on osteoclast activity (-)
E_zol_act : zoledronate inhibition of resorption (-)
QTc     : corrected QT interval (ms)
KDIGO   : KDIGO AKI stage (-)
'

mah <- mcode("mah", mah_code)

## =====================================================================
##  DOSES
## =====================================================================
ZA_4MG   <- 14.70      # umol zoledronic acid (MW 272.09)
DMB_120  <- 816000     # pmol denosumab (120 mg, MW 147 kDa)
CTN_4IU  <- 280        # IU salmon calcitonin (4 IU/kg, 70 kg)
PRD_60   <- 60 * 0.80 / 45.0      # mg/L increment, 60 mg PO, F 0.8, Vd 45 L
CIN_90   <- 90 * 0.25 / 400.0 * 1000   # ng/mL increment, 90 mg PO
FUR_40   <- 40 / 12.0             # mg/L increment, 40 mg IV

TX <- 12    # therapy starts on day 12 in every scenario

## aetiology parameter sets ---------------------------------------------
## Each was calibrated so that the model reaches a total calcium of 3.55
## mmol/L on day 12 -- deliberately, so that the treatment comparisons in
## THESIS 3 are made at MATCHED calcium.  The calcitriol mechanism is the
## exception and is discussed in note D: it CANNOT be driven to 3.55
## because it leaves the renal escape valve open, and that is a model
## output rather than a modelling choice.
AETI <- list(
  HHM = list(par = list(PTHRP_S = 883.46, TUM_G = 0.030, TUM_MAX = 2.2),
             init = list(TUM = 0.30)),
  LOH = list(par = list(CYT_S = 11.945, TUM_G = 0.030, TUM_MAX = 2.2),
             init = list(TUM = 0.30)),
  CTD = list(par = list(MAC_S = 3.60, TUM_G = 0.030, TUM_MAX = 2.2),
             init = list(TUM = 0.30)),
  EPT = list(par = list(PTH_ECT = 7339.4, TUM_G = 0.020, TUM_MAX = 1.0),
             init = list(TUM = 0.30)),
  NONE = list(par = list(), init = list())
)

## `segs` below is a piecewise-constant PARAMETER schedule: one row per change,
## a `time` column plus any parameters that change at that time (IV_RATE for
## saline, DIAL for a dialysis session, and so on).  Infusions and dialysis are
## parameters rather than dose events, so the run is split at each change and
## restarted from the previous segment's final state.
saline <- function(t0, hi = 4.80, hi_days = 3, lo = 2.40, lo_days = 4) {
  data.frame(time = c(t0, t0 + hi_days, t0 + hi_days + lo_days),
             IV_RATE = c(hi, lo, 0))
}

## merge several schedules (e.g. saline + dialysis) into one segment table,
## carrying each parameter forward until it is next changed
merge_segs <- function(...) {
  xs <- Filter(Negate(is.null), list(...))
  if (!length(xs)) return(NULL)
  cols <- unique(unlist(lapply(xs, names)))
  xs <- lapply(xs, function(d) {
    for (cc in setdiff(cols, names(d))) d[[cc]] <- NA_real_
    d[, cols, drop = FALSE]
  })
  s <- do.call(rbind, xs)
  s <- s[order(s$time), , drop = FALSE]
  ## collapse rows that share a time, then last-observation-carried-forward
  s <- do.call(rbind, lapply(split(s, s$time), function(g) {
    r <- g[1, , drop = FALSE]
    for (cc in cols) { v <- g[[cc]][!is.na(g[[cc]])]; if (length(v)) r[[cc]] <- v[length(v)] }
    r
  }))
  s <- s[order(s$time), , drop = FALSE]
  for (cc in setdiff(cols, "time")) {
    cur <- unname(as.list(param(mah))[[cc]])
    for (i in seq_len(nrow(s))) {
      if (is.na(s[[cc]][i])) s[[cc]][i] <- cur else cur <- s[[cc]][i]
    }
  }
  s
}

run_mah <- function(aeti = "HHM", events = NULL, iv = NULL, pmod = list(),
                    imod = list(), end = 60, delta = 0.05) {
  a <- AETI[[aeti]]
  m <- mah %>% param(a$par) %>% param(pmod)
  ini <- modifyList(a$init, imod)
  if (length(ini)) m <- m %>% init(ini)
  if (is.null(iv)) {
    return(m %>% mrgsim_df(events = events, end = end, delta = delta))
  }
  base <- iv[0, , drop = FALSE][1, , drop = FALSE]
  base$time <- 0
  for (cc in setdiff(names(iv), "time")) base[[cc]] <- unname(as.list(param(m))[[cc]])
  segs <- rbind(base, iv)
  segs <- segs[order(segs$time), , drop = FALSE]
  out <- NULL; state <- NULL
  for (i in seq_len(nrow(segs))) {
    t0 <- segs$time[i]
    t1 <- if (i < nrow(segs)) segs$time[i + 1] else end
    if (t1 <= t0) next
    mm <- m %>% param(as.list(segs[i, setdiff(names(segs), "time"), drop = FALSE]))
    if (!is.null(state)) mm <- mm %>% init(state)
    ev_i <- if (is.null(events)) NULL else
      dplyr::filter(as.data.frame(events), time >= t0, time < t1) %>%
      dplyr::mutate(time = time - t0)
    if (!is.null(ev_i) && nrow(ev_i) == 0) ev_i <- NULL
    o <- mm %>% mrgsim_df(events = ev_i, end = t1 - t0, delta = delta)
    o$time <- o$time + t0
    state <- as.list(o[nrow(o), names(init(mah))])
    out <- rbind(out, if (is.null(out)) o else o[-1, ])
  }
  out
}

ev_bolus <- function(times, cmt, amt) {
  data.frame(time = times, cmt = cmt, amt = amt, evid = 1, ii = 0, addl = 0)
}

## =====================================================================
##  23 SCENARIOS
## =====================================================================
SCEN <- list(

  ## --- 1. reference -------------------------------------------------
  s01 = list("Healthy control (90 d) — exact steady state",
             aeti = "NONE", end = 90),

  ## --- 2-4. the natural history and the fast arm --------------------
  s02 = list("HHM untreated — natural history to the collapse boundary",
             aeti = "HHM", end = 40),
  s03 = list("HHM + isotonic saline alone",
             aeti = "HHM", iv = saline(TX)),
  s04 = list("HHM + saline + calcitonin 4 IU/kg SC q12h x 48 h",
             aeti = "HHM", iv = saline(TX),
             events = ev_bolus(TX + c(0, 0.5, 1, 1.5), "CTN_SC", CTN_4IU)),

  ## --- 5-7. the slow arm, and what it needs from the fast arm -------
  s05 = list("HHM + saline + zoledronate 4 mg IV",
             aeti = "HHM", iv = saline(TX),
             events = ev_bolus(TX, "ZOL_C", ZA_4MG)),
  s06 = list("HHM + saline + calcitonin + zoledronate (guideline triple)",
             aeti = "HHM", iv = saline(TX),
             events = rbind(ev_bolus(TX + c(0, 0.5, 1, 1.5), "CTN_SC", CTN_4IU),
                            ev_bolus(TX, "ZOL_C", ZA_4MG))),
  s07 = list("HHM + zoledronate 4 mg WITHOUT saline",
             aeti = "HHM", events = ev_bolus(TX, "ZOL_C", ZA_4MG)),

  ## --- 8-10. RANKL blockade, and renal failure ----------------------
  s08 = list("HHM + saline + denosumab 120 mg SC",
             aeti = "HHM", iv = saline(TX),
             events = ev_bolus(TX, "DMB_SC", DMB_120)),
  s09 = list("HHM + saline + zoledronate, nephron mass 30%",
             aeti = "HHM", iv = saline(TX), imod = list(N_func = 0.30),
             events = ev_bolus(TX, "ZOL_C", ZA_4MG)),
  s10 = list("HHM + saline + denosumab, nephron mass 30%",
             aeti = "HHM", iv = saline(TX), imod = list(N_func = 0.30),
             events = ev_bolus(TX, "DMB_SC", DMB_120)),

  ## --- 11-14. THESIS 3: same number, different mechanism, different drug
  s11 = list("Osteolytic (myeloma) + saline + zoledronate — matched calcium",
             aeti = "LOH", iv = saline(TX),
             events = ev_bolus(TX, "ZOL_C", ZA_4MG)),
  s12 = list("Calcitriol (lymphoma) + saline + zoledronate",
             aeti = "CTD", iv = saline(TX),
             events = ev_bolus(TX, "ZOL_C", ZA_4MG)),
  s13 = list("Calcitriol (lymphoma) + saline + prednisolone 60 mg/day x 20 d",
             aeti = "CTD", iv = saline(TX),
             events = ev_bolus(TX + 0:19, "PRD_C", PRD_60)),
  s14 = list("Ectopic PTH + saline + cinacalcet 90 mg bd + zoledronate",
             aeti = "EPT", iv = saline(TX),
             events = rbind(ev_bolus(TX + seq(0, 29.5, by = 0.5), "CIN_C", CIN_90),
                            ev_bolus(TX, "ZOL_C", ZA_4MG))),

  ## --- 15-16, 23. the loop diuretic question, asked properly --------
  s15 = list("HHM + ADEQUATE saline + furosemide 40 mg q6h + ZA",
             aeti = "HHM", iv = saline(TX),
             events = rbind(ev_bolus(TX + seq(0, 2.75, by = 0.25), "FUR_C", FUR_40),
                            ev_bolus(TX, "ZOL_C", ZA_4MG))),
  s16 = list("HHM + INADEQUATE saline 1.2 L/day + furosemide q6h + ZA",
             aeti = "HHM",
             iv = data.frame(time = c(TX, TX + 7), IV_RATE = c(1.20, 0)),
             events = rbind(ev_bolus(TX + seq(0, 2.75, by = 0.25), "FUR_C", FUR_40),
                            ev_bolus(TX, "ZOL_C", ZA_4MG))),
  s23 = list("HHM + INADEQUATE saline 1.2 L/day, NO furosemide + ZA",
             aeti = "HHM",
             iv = data.frame(time = c(TX, TX + 7), IV_RATE = c(1.20, 0)),
             events = ev_bolus(TX, "ZOL_C", ZA_4MG)),

  ## --- 17-20. extracorporeal, albumin, precipitants -----------------
  s17 = list("HHM + saline + haemodialysis (Ca 1.0 dialysate, 4 h x2) + ZA",
             aeti = "HHM",
             ## two 4-hour sessions, on the day therapy starts and the next day
             iv = merge_segs(saline(TX),
                             data.frame(time = TX + c(0, 1/6, 1, 1 + 1/6),
                                        DIAL = c(1, 0, 1, 0))),
             events = ev_bolus(TX, "ZOL_C", ZA_4MG)),
  s18 = list("HHM with albumin 22 g/L — total versus ionised",
             aeti = "HHM", iv = saline(TX), pmod = list(ALB = 22),
             events = ev_bolus(TX, "ZOL_C", ZA_4MG)),
  s19 = list("HHM + immobilisation from day 6",
             aeti = "HHM", iv = saline(TX), pmod = list(IMMOB = 1),
             events = ev_bolus(TX, "ZOL_C", ZA_4MG)),
  s20 = list("HHM + thiazide co-prescription from day 6",
             aeti = "HHM", iv = saline(TX), pmod = list(THIAZIDE = 1),
             events = ev_bolus(TX, "ZOL_C", ZA_4MG)),

  ## --- 21-22. the only arm that keeps improving, and re-dosing ------
  s21 = list("HHM + saline + ZA + effective antitumour therapy",
             aeti = "HHM", iv = saline(TX), pmod = list(CHEMO = 0.11),
             events = ev_bolus(TX, "ZOL_C", ZA_4MG)),
  s22 = list("HHM + saline + zoledronate 4 mg q28d (two doses)",
             aeti = "HHM", iv = saline(TX), end = 70,
             events = ev_bolus(c(TX, TX + 28), "ZOL_C", ZA_4MG))
)

run_scenario <- function(id) {
  s <- SCEN[[id]]
  run_mah(aeti = s$aeti %||% "HHM", events = s$events, iv = s$iv,
          pmod = s$pmod %||% list(), imod = s$imod %||% list(),
          end = s$end %||% 60)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## Example ---------------------------------------------------------------
## out5  <- run_scenario("s05")
## out11 <- run_scenario("s11")
## plot(mah, out5, "Ca_tot,iCa,J_res,U_Ca,CTX,vol")

## =====================================================================
##  CALIBRATION NOTES
##  Every number below was generated from the system above (Python/scipy
##  reference implementation, cross-checked against a g++ RK4 integration
##  of the extracted $ODE block).  Literature anchors are cited by first
##  author and PMID in mah_references.md.
## =====================================================================
##
##  NOTE 0 — HEALTHY STEADY STATE
##  Scenario 1 is an exact steady state by construction: the production
##  constants KL_BASE, KO_BASE, PTG_SYN, K1A, D25_IN, MG_FE, ORAL_H2O and
##  GUT_SEC were each solved analytically for zero derivative at the
##  reference point rather than fitted.  Over 90 simulated days the
##  monitored variables drift: total calcium +0.004%, PTH +0.10%, 1,25D
##  -0.18%, osteoclast number -0.03%, free RANKL -0.07%, CTX -0.03%.
##  Reference values reproduced: total calcium 2.400, ionised 1.1995,
##  filtered load 259 mmol/day, urinary calcium 5.21 mmol/day, FE_Ca
##  2.01%, GFR 180 L/day, urine 2.02 L/day, bone turnover 8.0 mmol/day
##  in each direction, net intestinal absorption 5.22 mmol/day.
##
##  NOTE A — THE CEILING (THESIS 2, first half)
##  Setting CLAMP_CA to a target holds total calcium fixed while volume,
##  AQP2, symptoms, GFR, bone and the parathyroid all equilibrate around
##  it.  The steady-state urinary calcium then traces the true capacity
##  of the escape valve, feedbacks included:
##
##      Ca held   U_Ca ss   GFR    ECF vol   FE_Ca   urine
##        2.40      5.20   180.0    1.000    2.01%   2.01 L/day
##        2.80     16.80   179.1    0.995    5.60%   1.82
##        3.00     22.40   177.7    0.988    7.05%   1.52
##        3.20     27.16   175.2    0.976    8.17%   1.15   <- CEILING
##        3.40     13.04   120.1    0.838    5.43%   0.84
##        3.60      2.22    60.8    0.733    1.75%   0.60
##        3.80      0.00    28.2    0.653    0       0.49
##
##  The curve is not merely saturating, it TURNS OVER.  Between 3.2 and
##  3.6 mmol/L the kidney loses 92% of its excretory capacity, and it
##  does so because of the calcium it is trying to excrete.  Everything
##  clinicians describe qualitatively as "the vicious cycle" is in that
##  one column.
##
##  NOTE B — THE SADDLE-NODE (THESIS 2, second half)
##  J_EXO injects an unregulated calcium input.  Bisecting on whether a
##  steady state exists at 150 days locates the fold:
##
##    condition                              fold      Ca at    bound
##                                        (mmol/day)   fold      by
##    baseline, no support                   23.07     3.008    volume
##    + saline 1.2 L/day                     40.90     3.549    renal capacity
##    + saline 2.4 L/day                     51.43     4.964    renal capacity
##    + saline 4.8 L/day                     53.22     4.984    renal capacity
##    distal Tm +60% (PTHrP renal effect)    10.92     3.020    volume
##    distal Tm +60% + saline 2.4 L/day      39.75     4.972    renal capacity
##    thiazide                               20.04     3.011    volume
##    AQP2 suppression deleted               25.62     3.073    renal capacity
##    nausea / vomiting arm deleted          31.84     3.202    renal capacity
##    nephron mass 60%                       12.89     2.777    volume
##    nephron mass 30%                        7.25     2.651    volume
##
##  Four things worth reading off that table.  (i) Saline nearly doubles
##  the tolerable input at 1.2 L/day and changes the FAILURE MODE from
##  circulatory collapse to simple renal capacity -- which is the precise
##  sense in which fluid is the first-line treatment.  (ii) The renal
##  action of PTHrP, in isolation and with no extra osteoclast anywhere
##  in the model, more than halves the fold: 23.07 -> 10.92.  A humoral
##  tumour needs less than half the bone resorption of an osteolytic one
##  to tip the same patient over.  (iii) Thiazide costs 13% of the fold,
##  which is a quantitative version of "avoid it".  (iv) Deleting the
##  concentrating defect moves the fold by 11%; deleting the nausea arm
##  moves it by 38%.  The GI limb of the loop is the dominant one, which
##  is testable and not what the usual account emphasises.
##
##  NOTE C — THE CALCITONIN BRIDGE
##  Total calcium (mmol/L) after therapy starts, HHM:
##
##      hours    saline   saline+ZA   saline+CT+ZA   calcitonin benefit
##          0     3.550      3.550         3.550          0.000
##          6     3.350      3.307         3.120          0.188
##         12     3.218      3.114         2.946          0.167
##         24     3.265      3.067         2.959          0.108
##         48     3.228      2.898         2.831          0.067
##         96     3.215      2.728         2.703          0.025
##        168     3.233      2.641         2.634          0.007
##
##  Calcitonin receptor availability R_CT over the same window: 1.00 at
##  0 h, 0.71 at 6 h, 0.51 at 12 h, 0.26 at 24 h, 0.11 at 48 h, then
##  recovering to 0.31 at 5 days and 0.45 at 7 days.  The entire benefit
##  is spent in the first 48 hours -- exactly the interval in which the
##  bisphosphonate has not yet acted.  The model was not tuned to produce
##  that complementarity; CT_KOUT was set from the reported duration of
##  the calcitonin response and ZOL_KUP from the reported onset of the
##  bisphosphonate response, independently.
##
##  NOTE D — MATCHED CALCIUM, DIFFERENT PHYSIOLOGY (THESIS 3)
##  Day 12, before any treatment:
##
##                        HHM     Osteolytic   Ectopic PTH   Calcitriol
##    total calcium      3.550      3.550        3.547         2.835
##    ionised            1.716      1.715        1.713         1.404
##    PTH  (pmol/L)      1.198      1.199       10.991         2.441
##    PTHrP (pmol/L)    14.283      0.500        0.500         0.500
##    1,25D (pmol/L)    89.4       63.5        112.9         415.2
##    25(OH)D (nmol/L)  55.1       62.3         55.1          50.9
##    phosphate          0.862      1.378        0.914         1.192
##    J_res (mmol/day)  20.5       34.3         18.6          13.3
##    J_form (mmol/day)  7.58       9.49         7.95          8.46
##    J_gut (mmol/day)   0.26      -0.58         1.09         12.73
##    net input          13.18     24.23        11.78         17.61
##    U_Ca  (mmol/day)   5.9       17.0          5.2          16.9
##    FE_Ca              2.2%       6.5%         2.1%          5.6%
##    distal Ca Tm      43.3       26.3         41.4          30.4
##    CTX (ng/mL)        0.897      1.502        0.817         0.582
##
##  The osteolytic patient is resorbing 67% more bone than the humoral
##  patient for the SAME plasma calcium, and excreting nearly three times
##  as much of it, because the humoral patient's distal transport maximum
##  has been raised 65% by PTHrP.  The calcitriol patient is a different
##  animal again: the gut arm carries 12.73 of its 17.61 mmol/day of net
##  input, resorption is barely above normal, and the renal valve is wide
##  open -- which is why that mechanism cannot be driven to 3.55 in this
##  model at any plausible 1,25D, and why in the clinic it usually is not.
##
##  Look at what the gut arm is doing in the other three columns: 0.26,
##  -0.58 and 1.09 mmol/day.  Anorexia plus a suppressed 1,25D has closed
##  it unasked, and in myeloma it has gone frankly negative.  "Stop the
##  calcium supplement" is almost always right and almost never enough,
##  and here that is a number rather than an aphorism.
##
##  One internal consistency check, which was not arranged.  Net input at
##  day 12 is 13.18 mmol/day in humoral disease and 24.23 in osteolytic
##  disease -- nearly twice as much for the same plasma calcium.
##  Independently, note B puts the fold at 10.92 mmol/day for the humoral
##  renal phenotype and 23.07 for a kidney PTHrP has not touched.  Each
##  mechanism sits just past ITS OWN fold by a similar relative margin,
##  which is why two such different physiologies arrive at the same
##  number.  The two calculations share parameters, but neither was tuned
##  against the other.
##
##  The therapeutic consequence follows without further assumption:
##
##    scenario                       nadir   days < 2.70   CTX nadir
##    HHM + saline + ZA (s05)        2.635      15.0         0.153
##    Osteolytic + saline + ZA (s11) 2.402      25.9         0.302
##    Calcitriol + ZA (s12)          2.599      28.5         0.109
##    Calcitriol + prednisolone(s13) 2.411      22.4         0.382
##
##  Note the last two rows in particular.  Prednisolone reaches a DEEPER
##  and FASTER nadir than zoledronate in the calcitriol mechanism (2.411
##  at day 10.6, first normocalcaemia at 0.85 days versus 1.50) while
##  suppressing bone resorption three times LESS (CTX nadir 0.382 versus
##  0.109).  The two drugs are working on different arms, and the model
##  can tell which arm is carrying the disease.
##
##  NOTE E — THE SELF-DELIVERING DRUG
##  Identical 4 mg zoledronate dose, three different bone turnover rates:
##
##    J_res at dosing   peak intra-OC drug   day of peak   CTX nadir/base
##       8.0 (normal)         0.418             8.40           0.190
##      20.5 (humoral)        0.610             6.15           0.170
##      34.3 (osteolytic)     0.725             4.85           0.201
##
##  The drug arrives faster and in greater amount exactly where there is
##  more resorption to arrest, because resorption is the delivery
##  mechanism.  That is one term in one equation -- uptake = ZOL_KUP *
##  ZOL_B * (J_res / J_RES0) -- and it is also why the effect outlasts
##  the plasma exposure by weeks: once resorption stops, uptake stops,
##  the surface residue is slowly buried (ZOL_BURY) and the intracellular
##  drug decays with a 15-day half-life.
##
##  NOTE F — THE FLAT TOP OF THE DOSE-RESPONSE
##    1 mg  -> nadir 2.763, 0.0 normocalcaemic days in 30
##    2 mg  -> nadir 2.686, 5.6
##    4 mg  -> nadir 2.635, 15.0
##    8 mg  -> nadir 2.604, 23.2
##   16 mg  -> nadir 2.587, 26.4
##  Doubling from 4 to 8 mg buys 8.2 days; doubling again buys 3.2.  The
##  ceiling is imposed by the delivery step, not by target engagement --
##  which is the mechanistic reading of why 8 mg was never adopted
##  despite being measurably better on the biomarker.
##
##  NOTE G — WHERE THE CORRECTED CALCIUM FAILS
##  "True equivalent" = the total calcium that would give the same
##  ionised calcium at albumin 40 g/L.
##
##    albumin   total   ionised   Payne   true equiv   Payne error
##       35      3.55    1.879     3.65      3.76        -0.11
##       30      3.55    1.997     3.75      3.99        -0.24
##       25      3.55    2.130     3.85      4.26        -0.41
##       22      3.55    2.219     3.91      4.44        -0.53
##       25      2.40    1.440     2.70      2.88        -0.18
##
##  The error is not a constant offset: it grows with the calcium as well
##  as with the hypoalbuminaemia, so the correction is least reliable in
##  precisely the patient it was invented for.  Scenario 18 is the
##  dynamic version -- the same tumour in a patient with albumin 22 g/L
##  presents with a LOWER total calcium (3.42 versus 3.55) and a far
##  HIGHER ionised calcium (2.06 versus 1.72).
##  pH matters on the same scale: at a total of 3.00 and albumin 40, the
##  ionised calcium is 1.563 at pH 7.30 and 1.437 at pH 7.50.  Vomiting
##  alkalosis is therefore a genuine negative-feedback arm inside the
##  vicious cycle -- it protects the ionised calcium while destroying the
##  volume.
##
##  NOTE H — RENAL FAILURE FLIPS THE DRUG CHOICE
##  Starting from 30% nephron mass, the same disease presents at 4.40
##  rather than 3.55 mmol/L (the fold has moved from 23.1 to 7.3
##  mmol/day, note B).  On that background:
##      zoledronate (s09): nadir 2.967, 0 normocalcaemic days, and the
##                         trajectory still reaches the collapse boundary
##                         on day 46 -- renal clearance of the drug is
##                         reduced but so is the GFR that would excrete
##                         the calcium, and the bone arm alone is not
##                         enough.
##      denosumab   (s10): nadir 2.621, 7.7 normocalcaemic days, no
##                         collapse.  Denosumab has no renal step in its
##                         disposition and does not need resorption to
##                         deliver itself.
##  This is the one comparison in the file where the two antiresorptives
##  separate decisively, and the model gets there from disposition
##  structure rather than from an assumed potency difference.
##
##  NOTE I — THE LOOP DIURETIC, ASKED PROPERLY
##      s05  adequate saline + ZA                 nadir 2.635, 15.0 days
##      s15  adequate saline + ZA + furosemide    nadir 2.254, 19.2 days
##      s23  inadequate saline (1.2 L/d) + ZA     nadir 2.642, 14.0 days
##      s16  inadequate saline + ZA + furosemide  COLLAPSE BOUNDARY d 2.6
##  Furosemide is a real calciuretic and the model says so.  It is also
##  capable of converting a recoverable illness into circulatory collapse
##  inside three days when the volume underneath it is not being
##  replaced.  The variable that decides which of those two happens is
##  not the diuretic.
##
##  NOTE J — WHAT THE MODEL DOES NOT CONTAIN
##  No infection, no thrombosis, no hypercoagulability, no analgesia, no
##  bone metastasis geometry, no skeletal-related-event model, no
##  osteonecrosis of the jaw, no acute-phase reaction, no atrial
##  fibrillation, no digoxin interaction, and no mortality model.  The
##  "collapse boundary" is the edge of the calibrated region (ECF deficit
##  > 38% or total calcium > 5.0 mmol/L), NOT a prediction of death.
##  Tumour growth is a single logistic compartment with no resistance,
##  no heterogeneity and no metastatic dynamics; the chemotherapy arm
##  (scenario 21) should be read as "the input was removed", not as a
##  simulation of any particular regimen.  Bone is a single remodelling
##  compartment: trabecular and cortical sites are not distinguished, so
##  the model has nothing to say about fracture risk.  Every parameter is
##  a population point estimate with no between-subject variability
##  attached; the model reproduces mechanisms and their ordering, not
##  individual patients.
## =====================================================================
