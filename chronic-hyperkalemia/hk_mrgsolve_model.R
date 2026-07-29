## =============================================================================
##  CHRONIC HYPERKALAEMIA -- QSP / PK-PD MODEL  (mrgsolve)
##  Chronic hyperkalaemia of CKD and heart failure, and the RAASi dilemma
## =============================================================================
##
##  FILE          hk_mrgsolve_model.R
##  COMPANIONS    hk_qsp_model.dot/.svg/.png   mechanistic map
##                hk_reference_model.py        pure-Python reference + report
##                hk_model_report.txt          the numbers this file must match
##                hk_shiny_app.R               interactive dashboard
##                hk_references.md             sources for every parameter
##
##  RUN
##      library(mrgsolve); library(dplyr); library(ggplot2)
##      mod <- mread("hk_mrgsolve_model.R")
##      out <- scenario_all(mod)          # all 12 scenarios, see bottom of file
##
## -----------------------------------------------------------------------------
##  THE ONE IDEA
## -----------------------------------------------------------------------------
##  Serum potassium is not a pool.  It is the solution Ce of
##
##      K_total  =  Ce*V_ECF  +  Ci0 * LAMrel * (Ce/Ce0)^alpha * V_ICF
##
##  where LAMrel (the relative intracellular/extracellular partition) is set in
##  minutes by insulin, beta-2 tone, pH and tonicity, and K_total is moved only
##  over days by intake, kidney, colon and binders.  alpha = 0.25 is what makes
##  the intracellular space a BUFFER rather than a mirror, and it is the single
##  parameter that sets the exchange rate between "mmol/L on the report" and
##  "mmol in the patient":
##
##      fully equilibrated  224 mmol per mmol/L      (predicts, without being
##      fast pool only       66 mmol per mmol/L       told, the classical
##                                                    deficit nomogram)
##
##  Every therapy acts on exactly one term, and that -- not potency -- decides
##  when each is the right drug:
##
##      insulin/dextrose, salbutamol   -> LAMrel     0 mmol removed
##      calcium                        -> NEITHER    moves V_threshold only
##      binders, diuretics, diet, HD   -> K_total    mmol actually removed
##      alkali                         -> mostly the KIDNEY, not the cell
##
## -----------------------------------------------------------------------------
##  CALIBRATION (all fits performed in hk_reference_model.py; see section A)
## -----------------------------------------------------------------------------
##  FITTED (3 renal parameters to 3 anchors, 1 MR parameter, 1 acid parameter):
##    S0, adapt_p, adapt_max  <- steady-state serum K of 4.20 / 5.00 / 5.30
##                               at eGFR 100 / 20 / 12
##    KI_MRA                  <- RALES: spironolactone 25 mg raises K by +0.30
##    N_HCO3                  <- alkali therapy: HCO3 18->24 lowers K by 0.30
##
##  HELD OUT (predictions, mean |error| 0.06 mmol/L):
##    eGFR  60 -> K 4.47 (target 4.40)
##    eGFR  45 -> K 4.61 (target 4.55)
##    eGFR  30 -> K 4.80 (target 4.75)
##
##  EXTERNAL VALIDATION (nothing below was fitted):
##    OPAL-HK      patiromer 16.8 g/d, 4 wk    model -0.95  trial -1.01 mmol/L
##    HARMONIZE    SZC 5/10/15 g/d, 28 d       model 4.89/4.63/4.49
##                                             trial 4.8 /4.5 /4.4
##    HARMONIZE    SZC 10 g TID at 48 h        model -0.95  trial ~-1.1
##    insulin 10 U + D50                       model -0.85 at ~50 min
##                                             trial -0.6 to -1.0 at 30-60 min
##    salbutamol 20 mg neb                     model -0.74 at 2-4 h
##                                             trial -0.6 to -1.0
##    ACE inhibitor                            model +0.21 to +0.23
##                                             trial +0.1 to +0.4
##    deficit nomogram, stable serum K 3.0     model -302 mmol
##                                             literature -200 to -400 mmol
##
##  KNOWN DISCREPANCIES, carried rather than tuned away:
##    (1) FIDELIO finerenone 20 mg: assuming an MR load 35% of spironolactone
##        25 mg gives +0.10 vs the observed +0.23.  Back-solving needs 72%.
##        The model cannot separate "more occupancy" from "same occupancy,
##        different tissue distribution", which is precisely the non-steroidal
##        MRA claim.  Reported as unidentifiable, not fitted.
##    (2) The model does not reproduce the 15-20% incidence of hypoglycaemia
##        after insulin/dextrose: its insulin effect site decays with the
##        plasma level, and there is no between-patient glycogen reserve.
##    (3) The MRA potassium cost comes out FLAT across eGFR (+0.30 to +0.28
##        from eGFR 90 to 15), contradicting the prior that it should steepen.
##        What steepens is the baseline it lands on, not the increment.
##
## -----------------------------------------------------------------------------
##  UNITS.  Time = days.  Amounts = mmol (potassium), mg (drugs), g (binders).
##          Concentrations = mmol/L, mg/L.  eGFR = mL/min/1.73 m2.
## =============================================================================

$PROB
# Chronic hyperkalaemia (CKD / heart failure) -- QSP model
# 30 compartments: whole-body K distribution, renal ASDN, colon, RAAS,
# acid-base, RAASi PK, binder gut PK, acute rescue PK, ECG and outcomes.

$PARAM @annotated
// ---- anthropometry and volumes ---------------------------------------------
BW       :  70.0  : Body weight (kg)
FV_ECF   :  0.20  : ECF as fraction of body weight (-)
FV_ICF   :  0.36  : ICF as fraction of body weight (-)
F_FAST   :  0.25  : Fraction of ICF in the fast-exchanging pool (-)
CI0      : 140.0  : Reference intracellular potassium (mmol/L)
CE0      :   4.2  : Reference serum potassium (mmol/L)

// ---- transcellular partition ------------------------------------------------
ALPHA    :  0.25  : Partition exponent, Ci_target ~ Ce^ALPHA (-)
KTC_F    : 13.30  : Fast-pool exchange rate constant (1/day)
KTC_S    :  0.277 : Slow (muscle) pool exchange rate constant (1/day)
EMAX_INS :  0.120 : Max fractional rise in LAMrel from insulin (-)
EC50_INS : 60.0   : Insulin concentration for half effect (uU/mL)
EMAX_B2  :  0.070 : Max fractional rise in LAMrel from beta-2 agonism (-)
EC50_B2  :  0.012 : Beta-2 agonist effect-site EC50 (mg/L)
K_PH     :  0.190 : Transcellular pH coefficient (1/pH unit)
K_GLU    :  0.0025: Tonicity coefficient (L/mmol glucose)
EMAX_ALDC:  0.012 : Aldosterone effect on cellular uptake (-)

// ---- diet and gut -----------------------------------------------------------
INTAKE   : 80.0   : Dietary potassium intake (mmol/day)
KA_GUT   : 36.0   : Small-bowel absorption rate constant (1/day)
KT_GUT   :  4.0   : Proximal-to-colon transit rate constant (1/day)
KST_COL  :  4.0   : Colon-to-stool rate constant (1/day)
C0_COL   : 10.0   : Net colonic K secretion at reference GFR (mmol/day)
R_COL    :  0.50  : Fraction of secreted colonic K reabsorbed (-)
KC_COL   :  1.00  : Colonic up-regulation slope with nephron loss (-)
KM_COL   :  4.2   : Serum-K half-max for colonic secretion (mmol/L)
M_COL    :  0.50  : MR-occupancy exponent for colonic secretion (-)

// ---- renal handling ---------------------------------------------------------
GFR0     : 100.0  : Reference eGFR (mL/min/1.73)
FD       :  0.030 : Fraction of filtered K delivered to the ASDN (-)
S0       : 75.182 : ASDN secretory capacity at GFR0 (mmol/day) [FITTED]
ADAPT_P  :  0.879 : Per-nephron secretory up-regulation exponent (-) [FITTED]
ADAPT_MX :  6.361 : Ceiling on per-nephron up-regulation (-) [FITTED]
KM_SECR  :  3.0   : Serum-K half-max for distal secretion (mmol/L)
FLOW_EXP :  0.50  : Distal-flow (BK-channel) exponent (-)
N_HCO3   :  0.435 : Acidosis inhibition of distal K secretion (-) [FITTED]

// ---- aldosterone and the MR -------------------------------------------------
ALDO0    : 12.0   : Baseline plasma aldosterone (ng/dL)
SK_ALDO  :  1.30  : Direct potassium drive on aldosterone (L/mmol)
KALDO    : 50.0   : Aldosterone turnover rate constant (1/day)
KD_MR    :  8.0   : Aldosterone Kd at the MR (ng/dL)
OCC0     :  0.60  : Baseline MR occupancy (-)
M_ASDN   :  0.80  : Occupancy-to-transporter exponent (-)
K_ASDN   :  0.462 : ASDN transporter turnover, t1/2 1.5 d (1/day)

// ---- acid-base --------------------------------------------------------------
NEAP     : 70.0   : Net endogenous acid production (mmol/day)
NAE_EXP  :  0.60  : GFR scaling of renal net acid excretion (-)
K_BUF    :  0.70  : Bone and muscle buffering rate constant (1/day)
HCO3_0   : 24.0   : Reference serum bicarbonate (mmol/L)
PCO2_0   : 40.0   : Reference arterial pCO2 (mmHg)
F_BIC    :  0.40  : Systemic fraction of an oral alkali dose (-)
KA_BIC   : 12.0   : Oral alkali absorption rate constant (1/day)

// ---- ACE inhibitor ----------------------------------------------------------
KA_ACE   : 24.0   : ACEi absorption rate constant (1/day)
CL_ACE   : 72.0   : ACEi clearance (L/day)
V_ACE    : 40.0   : ACEi volume of distribution (L)
F_ACE    :  0.25  : ACEi bioavailability (-)
EMAX_ACE :  0.45  : Max suppression of the angiotensin II drive (-)
EC50_ACE :  0.020 : ACEi EC50 (mg/L)

// ---- mineralocorticoid receptor antagonist ----------------------------------
KA_MRA   : 28.8   : MRA absorption rate constant (1/day)
CL_MRA   : 240.0  : MRA clearance (L/day)
V_MRA    : 200.0  : MRA volume of distribution (L)
F_MRA    :  0.70  : MRA bioavailability (-)
KI_MRA   :  0.1029: MRA competitive Ki at the MR (mg/L) [FITTED to RALES]

// ---- potassium binders ------------------------------------------------------
PHIMAX_P :  0.42  : Max luminal K captured by patiromer (-)
D50_PAT  : 12.0   : Patiromer dose for half-max capture (g/day)
KT_PAT   :  4.0   : Patiromer proximal-to-colon transit (1/day)
KEL_PAT  :  4.0   : Patiromer colonic elimination (1/day)
PHIMAX_S :  0.45  : Max luminal K captured by SZC (-)
D50_SZC  :  8.0   : SZC dose for half-max capture (g/day)
KT_SZC   :  8.0   : SZC proximal transit (1/day)
KEL_SZC  :  8.0   : SZC distal elimination (1/day)

// ---- insulin and glucose ----------------------------------------------------
INS_BASE : 10.0   : Basal plasma insulin (uU/mL)
KEL_INS  : 166.0  : Insulin elimination, t1/2 ~6 min (1/day)
KE0_INS  : 50.0   : Insulin effect-site equilibration (1/day)
V_INS    : 12.0   : Insulin distribution volume (L)
GLU_BASE :  5.5   : Reference plasma glucose (mmol/L)
KEL_GLU  : 30.0   : Insulin-independent glucose disposal (1/day)
K_INSGLU :  0.010 : Insulin-dependent glucose disposal (L/uU/day)
EGP0     : 173.25 : Endogenous glucose production (mmol/L/day)
EC_EGP   : 200.0  : Insulin suppression of EGP (uU/mL)
V_GLU    : 16.0   : Glucose distribution volume (L)

// ---- beta-2 agonist ---------------------------------------------------------
CL_B2    : 720.0  : Salbutamol clearance (L/day)
V_B2     : 150.0  : Salbutamol volume of distribution (L)
F_B2     :  0.20  : Pulmonary bioavailability of a nebulised dose (-)
KE0_B2   : 12.0   : Salbutamol effect-site equilibration (1/day)

// ---- loop diuretic ----------------------------------------------------------
KA_FUR   : 48.0   : Furosemide absorption rate constant (1/day)
CL_FUR   : 216.0  : Furosemide clearance (L/day)
V_FUR    : 12.0   : Furosemide volume of distribution (L)
F_FUR    :  0.50  : Furosemide bioavailability (-)
EMAX_FUR :  1.60  : Max multiple of baseline distal flow (-)
EC50_FUR :  0.35  : Furosemide EC50 (mg/L)
SGLT2I   :  0.0   : SGLT2 inhibitor on/off (0/1)
FLOW_SGL :  1.15  : Distal-flow multiplier from an SGLT2 inhibitor (-)

// ---- calcium and electrophysiology -----------------------------------------
KEL_CA   : 16.0   : Ionised calcium increment decay (1/day)
CA_GAIN  :  0.25  : Ionised Ca rise per gram calcium gluconate (mmol/L/g)
S_CA     : 15.0   : Threshold shift per mmol/L ionised calcium (mV per mmol/L)
VTH0     : -70.0  : Threshold potential at normal ionised calcium (mV)
H_MID    : -78.0  : Sodium-channel availability midpoint (mV)
H_SLOPE  :   6.0  : Sodium-channel availability slope (mV)
QRS0     :  90.0  : Reference QRS duration (ms)

// ---- progression, outcome and the clinician --------------------------------
SLOPE0   : -4.0   : Untreated diabetic CKD eGFR slope (mL/min/1.73/yr)
SLOPE_RD :  1.40  : eGFR slope recovered at full RAASi dose (mL/min/1.73/yr)
SLOPE_MRA:  0.60  : Additional eGFR slope from an MRA (mL/min/1.73/yr)
LNHR_RD  : -0.248 : ln HR of full-dose RAASi (-)
LNHR_MRA : -0.198 : ln HR of an MRA on top (-)
BK_HI    :  0.642 : Hazard coefficient above K 5.0 (-)
PK_HI    :  1.20  : Hazard exponent above K 5.0 (-)
BK_LO    :  0.916 : Hazard coefficient below K 4.0 (-)
PK_LO    :  1.50  : Hazard exponent below K 4.0 (-)
K_DOWN   :  0.90  : Clinician down-titration gain (1/day per mmol/L)
K_UP     :  0.020 : Clinician up-titration gain (1/day per mmol/L)
K_STOP   :  5.50  : Serum K that triggers dose reduction (mmol/L)
K_SAFE   :  5.00  : Serum K below which up-titration resumes (mmol/L)

// ---- switches and infusion rates (set per scenario) ------------------------
PROGRESS :  0.0   : Allow eGFR to decline over time (0/1)
TITRATE  :  0.0   : Allow the clinician controller to move the dose (0/1)
RD_FIX   :  1.0   : Prescribed RAASi dose when TITRATE = 0 (fraction of target)
ACE_TGT  :  0.0   : Guideline-target ACE inhibitor dose (mg/day)
MRA_TGT  :  0.0   : Guideline-target MRA dose (mg/day)
PAT_RATE :  0.0   : Patiromer dose rate (g/day)
SZC_RATE :  0.0   : Sodium zirconium cyclosilicate dose rate (g/day)
FUR_RATE :  0.0   : Furosemide dose rate (mg/day)
BIC_RATE :  0.0   : Oral alkali dose rate (mmol/day)
KCL_RATE :  0.0   : Exogenous IV/oral KCl load (mmol/day)
GFR_SET  : 25.0   : Baseline eGFR (mL/min/1.73)
HCO3_SET :  0.0   : Override serum HCO3 (mmol/L; 0 = let the model set it)

$CMT @annotated
KE     : ECF potassium (mmol)
KIF    : Fast-pool intracellular potassium (mmol)
KIS    : Slow-pool (muscle) intracellular potassium (mmol)
HCO3A  : ECF bicarbonate (mmol)
ALDO   : Plasma aldosterone (ng/dL)
RASDN  : ASDN transporter abundance (relative)
GUTK   : Proximal luminal potassium (mmol)
COLK   : Colonic luminal potassium (mmol)
AACE   : ACE inhibitor absorption depot (mg)
CACE   : ACE inhibitor central concentration (mg/L)
AMRA   : MRA absorption depot (mg)
CMRA   : MRA central concentration (mg/L)
PATP   : Patiromer, proximal gut (g)
PATC   : Patiromer, colon = site of action (g)
SZCP   : SZC, stomach and proximal gut = site of action (g)
SZCC   : SZC, distal gut (g)
INS    : Plasma insulin (uU/mL)
INSE   : Insulin effect site (uU/mL)
GLU    : Plasma glucose (mmol/L)
B2C    : Beta-2 agonist plasma concentration (mg/L)
B2E    : Beta-2 agonist effect site (mg/L)
AFUR   : Furosemide depot (mg)
CFUR   : Furosemide central concentration (mg/L)
CAE    : Ionised calcium increment (mmol/L)
GFR    : Estimated GFR (mL/min/1.73)
RD     : Prescribed RAASi dose (fraction of guideline target)
BICG   : Oral alkali gut depot (mmol)
T55    : Cumulative time with serum K above 5.5 (days)
CHAZ   : Cumulative composite hazard (day)
KREM   : Cumulative potassium removed by binder (mmol)

$GLOBAL
#define V_ECF   (FV_ECF*BW)
#define V_ICF   (FV_ICF*BW)
#define V_IF    (F_FAST*V_ICF)
#define V_IS    ((1.0-F_FAST)*V_ICF)

$MAIN
// ---- initial conditions -----------------------------------------------------
// Seeded near the eGFR-appropriate steady state; run the model in for 200 days
// (see the scenarios at the bottom) before reading anything off it.
double gfr_0   = GFR_SET;
double hco3_ss = (HCO3_SET > 0) ? HCO3_SET
                 : HCO3_0 + (NEAP*pow(gfr_0/GFR0, NAE_EXP) - NEAP)/(K_BUF*V_ECF);
double ce_0    = CE0 + 0.9*(GFR0 - gfr_0)/GFR0;      // crude seed, refined by run-in

KE_0    = ce_0 * V_ECF;
KIF_0   = CI0 * pow(ce_0/CE0, ALPHA) * V_IF;
KIS_0   = CI0 * pow(ce_0/CE0, ALPHA) * V_IS;
HCO3A_0 = hco3_ss * V_ECF;
ALDO_0  = ALDO0 * exp(SK_ALDO*(ce_0 - CE0));
RASDN_0 = 1.0;
GUTK_0  = INTAKE/(KA_GUT + KT_GUT);
COLK_0  = 3.0;
INS_0   = INS_BASE;
INSE_0  = INS_BASE;
GLU_0   = GLU_BASE;
GFR_0   = gfr_0;
RD_0    = RD_FIX;

$ODE
// =============================================================================
//  ALGEBRA
// =============================================================================
double Ce   = KE/V_ECF;
double Cif  = KIF/V_IF;
double Cis  = KIS/V_IS;
double hco3 = HCO3A/V_ECF;
double gfr  = (GFR < 3.0) ? 3.0 : GFR;

// ---- acid-base: Winter's respiratory compensation ---------------------------
double pco2 = PCO2_0 - 1.2*(HCO3_0 - hco3);
if (pco2 < 15.0) pco2 = 15.0;
double pH   = 6.10 + log10((hco3 < 1e-6 ? 1e-6 : hco3)/(0.03*pco2));

// ---- mineralocorticoid receptor occupancy -----------------------------------
double kd_app = KD_MR*(1.0 + CMRA/KI_MRA);
double occ    = ALDO/(ALDO + kd_app);

// ---- transcellular partition (LAMrel) and the two exchange fluxes ----------
// Insulin acts over its WHOLE range, not just above basal: basal insulin is
// itself doing most of the work, which is why insulinopenia raises serum
// potassium with no change whatsoever in total body potassium.
double occ_i = INSE/(EC50_INS + INSE);
double occ_b = INS_BASE/(EC50_INS + INS_BASE);
double f_ins = 1.0 + EMAX_INS*(occ_i - occ_b)/(1.0 - occ_b);
double f_b2  = 1.0 + EMAX_B2*B2E/(EC50_B2 + B2E);
double f_pH  = exp(K_PH*(pH - 7.40));
double f_glu = 1.0/(1.0 + K_GLU*((GLU > GLU_BASE) ? (GLU - GLU_BASE) : 0.0));
double f_ald = 1.0 + EMAX_ALDC*(occ - OCC0)/OCC0;
double LAMrel = f_ins*f_b2*f_pH*f_glu*f_ald;

double ci_t = CI0*LAMrel*pow((Ce < 1e-6 ? 1e-6 : Ce)/CE0, ALPHA);
double J_f  = KTC_F*V_IF*(ci_t - Cif);
double J_s  = KTC_S*V_IS*(ci_t - Cis);

// ---- binder capture fraction ------------------------------------------------
// Driven by the dose rate THROUGH the site of action, so the dynamic model and
// the algebraic steady state stay on the same scale.  Patiromer acts distally
// (colon), SZC from the stomach onward -- which is why SZC is the faster of the
// two and why both work at all in a patient who is not eating.
double rate_p = KEL_PAT*PATC;
double rate_s = KT_SZC*SZCP;
double phi_p  = (rate_p > 0) ? PHIMAX_P*rate_p/(D50_PAT + rate_p) : 0.0;
double phi_s  = (rate_s > 0) ? PHIMAX_S*rate_s/(D50_SZC + rate_s) : 0.0;
double phi    = 1.0 - (1.0 - phi_p)*(1.0 - phi_s);

// ---- renal potassium excretion ----------------------------------------------
double filt  = gfr*1.44*Ce;
double deliv = FD*filt;
double adapt = pow(GFR0/gfr, ADAPT_P);
if (adapt > ADAPT_MX) adapt = ADAPT_MX;
double scap  = S0*(gfr/GFR0)*adapt;
double qrel  = (1.0 + (EMAX_FUR - 1.0)*CFUR/(EC50_FUR + CFUR))
               * (SGLT2I > 0.5 ? FLOW_SGL : 1.0);
double g_acid = pow(((hco3 < 4.0) ? 4.0 : hco3)/HCO3_0, N_HCO3);
double E_ren = deliv + scap*RASDN*(Ce/(KM_SECR + Ce))*pow(qrel, FLOW_EXP)*g_acid;

// ---- colonic secretion ------------------------------------------------------
// Secretion is GROSS; half is normally reabsorbed.  A binder in the lumen
// blocks that reabsorption, which is the whole reason a binder still works in
// a fasting patient.
double gross0 = C0_COL/(1.0 - R_COL);
double gfr_rel = (gfr/GFR0 > 1.0) ? 1.0 : gfr/GFR0;
double col_gr = gross0*(1.0 + KC_COL*(1.0 - gfr_rel))
                * ((Ce/(Ce + KM_COL))/(CE0/(CE0 + KM_COL)))
                * pow(occ/OCC0, M_COL);
double E_col  = col_gr*(1.0 - R_COL*(1.0 - phi));

// ---- gut --------------------------------------------------------------------
double abs_K = KA_GUT*GUTK*(1.0 - phi);

// ---- prescribed RAASi dose --------------------------------------------------
double rd = (RD < 0.0) ? 0.0 : ((RD > 1.0) ? 1.0 : RD);

// =============================================================================
//  DIFFERENTIAL EQUATIONS
// =============================================================================
dxdt_KE    = abs_K + KCL_RATE - E_ren - E_col - J_f - J_s;
dxdt_KIF   = J_f;
dxdt_KIS   = J_s;

dxdt_GUTK  = INTAKE - KA_GUT*GUTK - KT_GUT*GUTK;
dxdt_COLK  = KT_GUT*GUTK + E_col - KST_COL*COLK;

dxdt_HCO3A = NEAP*pow(gfr_rel, NAE_EXP) - NEAP
             + F_BIC*KA_BIC*BICG + K_BUF*V_ECF*(HCO3_0 - hco3);
dxdt_BICG  = BIC_RATE - KA_BIC*BICG;

double f_ace   = 1.0 - EMAX_ACE*CACE/(EC50_ACE + CACE);
double aldo_ss = ALDO0*f_ace*exp(SK_ALDO*(Ce - CE0));
dxdt_ALDO  = KALDO*(aldo_ss - ALDO);
dxdt_RASDN = K_ASDN*(pow(occ/OCC0, M_ASDN) - RASDN);

dxdt_AACE  = ACE_TGT*rd - KA_ACE*AACE;
dxdt_CACE  = F_ACE*KA_ACE*AACE/V_ACE - CL_ACE/V_ACE*CACE;
dxdt_AMRA  = MRA_TGT*rd - KA_MRA*AMRA;
dxdt_CMRA  = F_MRA*KA_MRA*AMRA/V_MRA - CL_MRA/V_MRA*CMRA;

dxdt_PATP  = PAT_RATE - KT_PAT*PATP;
dxdt_PATC  = KT_PAT*PATP - KEL_PAT*PATC;
dxdt_SZCP  = SZC_RATE - KT_SZC*SZCP;
dxdt_SZCC  = KT_SZC*SZCP - KEL_SZC*SZCC;

dxdt_INS   = -KEL_INS*(INS - INS_BASE);
dxdt_INSE  = KE0_INS*(INS - INSE);
dxdt_GLU   = EGP0/(1.0 + INSE/EC_EGP)
             - (KEL_GLU + K_INSGLU*((INSE > INS_BASE) ? (INSE - INS_BASE) : 0.0))*GLU;
dxdt_B2C   = -CL_B2/V_B2*B2C;
dxdt_B2E   = KE0_B2*(B2C - B2E);

dxdt_AFUR  = FUR_RATE - KA_FUR*AFUR;
dxdt_CFUR  = F_FUR*KA_FUR*AFUR/V_FUR - CL_FUR/V_FUR*CFUR;
dxdt_CAE   = -KEL_CA*CAE;

// ---- progression and the clinician in the loop ------------------------------
double mra_on = (MRA_TGT*rd > 0.0) ? 1.0 : 0.0;
double slope  = SLOPE0 + SLOPE_RD*rd + SLOPE_MRA*mra_on;
dxdt_GFR = (PROGRESS > 0.5) ? slope/365.0 : 0.0;

double dn = (Ce > K_STOP)  ? (Ce - K_STOP)  : 0.0;
double up = (Ce < K_SAFE)  ? (K_SAFE - Ce)  : 0.0;
dxdt_RD  = (TITRATE > 0.5) ? (-K_DOWN*dn*rd + K_UP*up*(1.0 - rd)) : 0.0;

// ---- trackers ---------------------------------------------------------------
double hz_K = (Ce > 5.0) ? exp(BK_HI*pow(Ce - 5.0, PK_HI))
            : ((Ce < 4.0) ? exp(BK_LO*pow(4.0 - Ce, PK_LO)) : 1.0);
dxdt_T55  = (Ce > 5.5) ? 1.0 : 0.0;
dxdt_CHAZ = hz_K*exp(LNHR_RD*rd)*exp(LNHR_MRA*mra_on);
dxdt_KREM = phi*KA_GUT*GUTK + phi*R_COL*col_gr;

$TABLE
double K       = KE/V_ECF;
double KI_F    = KIF/V_IF;
double KI_S    = KIS/V_IS;
double KTOT    = KE + KIF + KIS;
double HCO3ser = HCO3A/V_ECF;
double PCO2    = PCO2_0 - 1.2*(HCO3_0 - HCO3ser);
if (PCO2 < 15.0) PCO2 = 15.0;
double PHA     = 6.10 + log10((HCO3ser < 1e-6 ? 1e-6 : HCO3ser)/(0.03*PCO2));

double GFRo    = (GFR < 3.0) ? 3.0 : GFR;
double FILT    = GFRo*1.44*K;
double ADAPTo  = pow(GFR0/GFRo, ADAPT_P);
if (ADAPTo > ADAPT_MX) ADAPTo = ADAPT_MX;
double SCAP    = S0*(GFRo/GFR0)*ADAPTo;
double QREL    = (1.0 + (EMAX_FUR - 1.0)*CFUR/(EC50_FUR + CFUR))
                 * (SGLT2I > 0.5 ? FLOW_SGL : 1.0);
double GACID   = pow(((HCO3ser < 4.0) ? 4.0 : HCO3ser)/HCO3_0, N_HCO3);
double UK      = FD*FILT + SCAP*RASDN*(K/(KM_SECR + K))*pow(QREL, FLOW_EXP)*GACID;
double FEK     = 100.0*UK/FILT;

double KDAPP   = KD_MR*(1.0 + CMRA/KI_MRA);
double MROCC   = ALDO/(ALDO + KDAPP);

// LAMrel must be recomputed here: variables declared in $ODE are local to the
// ODE function and are not visible to $TABLE or $CAPTURE.
double OCCI    = INSE/(EC50_INS + INSE);
double OCCB    = INS_BASE/(EC50_INS + INS_BASE);
double FINS    = 1.0 + EMAX_INS*(OCCI - OCCB)/(1.0 - OCCB);
double FB2     = 1.0 + EMAX_B2*B2E/(EC50_B2 + B2E);
double FPH     = exp(K_PH*(PHA - 7.40));
double FGLU    = 1.0/(1.0 + K_GLU*((GLU > GLU_BASE) ? (GLU - GLU_BASE) : 0.0));
double FALD    = 1.0 + EMAX_ALDC*(MROCC - OCC0)/OCC0;
double LAMrel  = FINS*FB2*FPH*FGLU*FALD;

double RATEP   = KEL_PAT*PATC;
double RATES   = KT_SZC*SZCP;
double PHIP    = (RATEP > 0) ? PHIMAX_P*RATEP/(D50_PAT + RATEP) : 0.0;
double PHIS    = (RATES > 0) ? PHIMAX_S*RATES/(D50_SZC + RATES) : 0.0;
double PHI     = 1.0 - (1.0 - PHIP)*(1.0 - PHIS);

double GROSS0  = C0_COL/(1.0 - R_COL);
double GFRREL  = (GFRo/GFR0 > 1.0) ? 1.0 : GFRo/GFR0;
double COLGR   = GROSS0*(1.0 + KC_COL*(1.0 - GFRREL))
                 * ((K/(K + KM_COL))/(CE0/(CE0 + KM_COL)))*pow(MROCC/OCC0, M_COL);
double COLK_EX = COLGR*(1.0 - R_COL*(1.0 - PHI));

// ---- electrophysiology ------------------------------------------------------
double EM      = 61.5*log10(K/KI_F);
double VTH     = VTH0 + S_CA*CAE;
double GAP     = VTH - EM;
double HAVAIL  = 1.0/(1.0 + exp((EM - H_MID)/H_SLOPE));
double HREF    = 1.0/(1.0 + exp((61.5*log10(CE0/CI0) - H_MID)/H_SLOPE));
double QRS     = QRS0*sqrt(HREF/HAVAIL);
double TWAVE   = sqrt(K/CE0)*(1.0 + 0.35*((K > 5.0) ? (K - 5.0) : 0.0));

// ---- buffer capacity: the mmol-per-mmol/L exchange rate ---------------------
double DCI     = CI0*LAMrel*ALPHA/CE0*pow(K/CE0, ALPHA - 1.0);
double BUF_CHR = V_ECF + (V_IF + V_IS)*DCI;
double BUF_ACU = V_ECF + V_IF*DCI;

// ---- outcome ----------------------------------------------------------------
double HZK     = (K > 5.0) ? exp(BK_HI*pow(K - 5.0, PK_HI))
               : ((K < 4.0) ? exp(BK_LO*pow(4.0 - K, PK_LO)) : 1.0);
double RDo     = (RD < 0.0) ? 0.0 : ((RD > 1.0) ? 1.0 : RD);
double MRAON   = (MRA_TGT*RDo > 0.0) ? 1.0 : 0.0;
double HZTOT   = HZK*exp(LNHR_RD*RDo)*exp(LNHR_MRA*MRAON);

$CAPTURE @annotated
K       : Serum potassium (mmol/L)
KI_F    : Fast-pool intracellular potassium (mmol/L)
KI_S    : Muscle intracellular potassium (mmol/L)
KTOT    : Total body potassium (mmol)
LAMrel  : Relative transcellular partition (-)
BUF_CHR : Chronic buffer capacity (mmol per mmol/L)
BUF_ACU : Acute buffer capacity (mmol per mmol/L)
HCO3ser : Serum bicarbonate (mmol/L)
PHA     : Arterial pH (-)
UK      : Urine potassium excretion (mmol/day)
FEK     : Fractional excretion of potassium (%)
FILT    : Filtered potassium load (mmol/day)
SCAP    : ASDN secretory capacity (mmol/day)
COLK_EX : Net colonic potassium excretion (mmol/day)
COLGR   : Gross colonic potassium secretion (mmol/day)
PHI     : Fraction of luminal potassium captured by binder (-)
MROCC   : Mineralocorticoid receptor occupancy (-)
EM      : Resting membrane potential (mV)
VTH     : Threshold potential (mV)
GAP     : Excitability gap (mV)
HAVAIL  : Sodium-channel availability (-)
QRS     : QRS duration (ms)
TWAVE   : T-wave amplitude, relative (-)
HZK     : Potassium-attributable hazard ratio (-)
HZTOT   : Composite hazard ratio (-)
RDo     : Prescribed RAASi dose, fraction of target (-)

## =============================================================================
##  R-SIDE HELPERS AND THE TWELVE SCENARIOS
##  Everything below is ordinary R and is ignored by mread().
## =============================================================================
$ENV

## Standard regimens, expressed as parameter sets (the model takes dose RATES
## rather than dosing events for the chronic agents, which keeps a 5-year
## simulation cheap; the acute agents are given as ordinary mrgsolve boluses).
RUN_IN <- 250      # days of run-in to reach steady state before t = 0

.p <- function(...) list(...)

REGIMENS <- list(
  ckd4_naive         = .p(GFR_SET = 25, ACE_TGT = 0,  MRA_TGT = 0),
  ckd4_acei          = .p(GFR_SET = 25, ACE_TGT = 20, MRA_TGT = 0),
  ckd4_acei_mra      = .p(GFR_SET = 25, ACE_TGT = 20, MRA_TGT = 25),
  ckd4_patiromer     = .p(GFR_SET = 25, ACE_TGT = 20, MRA_TGT = 25, PAT_RATE = 16.8),
  ckd4_szc           = .p(GFR_SET = 25, ACE_TGT = 20, MRA_TGT = 25, SZC_RATE = 10),
  ckd4_lowk          = .p(GFR_SET = 25, ACE_TGT = 20, MRA_TGT = 25, INTAKE   = 50),
  ckd4_loop          = .p(GFR_SET = 25, ACE_TGT = 20, MRA_TGT = 25, FUR_RATE = 40),
  ckd4_sglt2i        = .p(GFR_SET = 25, ACE_TGT = 20, MRA_TGT = 25, SGLT2I   = 1),
  ckd4_alkali        = .p(GFR_SET = 25, ACE_TGT = 20, MRA_TGT = 25, BIC_RATE = 70),
  ckd5_acei_mra      = .p(GFR_SET = 12, ACE_TGT = 20, MRA_TGT = 25),
  hf_preserved_gfr   = .p(GFR_SET = 55, ACE_TGT = 20, MRA_TGT = 25),
  ckd4_saltsubstitute= .p(GFR_SET = 25, ACE_TGT = 20, MRA_TGT = 25, INTAKE   = 130)
)

## ---------------------------------------------------------------------------
## sim_steady(): run a regimen to steady state and return the last row.
## ---------------------------------------------------------------------------
sim_steady <- function(mod, pars = list(), tend = RUN_IN) {
  mod %>%
    param(pars) %>%
    mrgsim(end = tend, delta = 1, atol = 1e-8, rtol = 1e-8) %>%
    as.data.frame() %>%
    dplyr::filter(time == max(time))
}

## ---------------------------------------------------------------------------
## SCENARIO 1-3.  The eGFR-potassium curve and the kaliuresis reserve.
##   Reproduces section A/B of hk_model_report.txt: serum K is FLAT against
##   eGFR until the reserve is spent, and the fractional excretion of potassium
##   is the gauge that was moving the whole time.
## ---------------------------------------------------------------------------
scenario_reserve <- function(mod, gfrs = c(90, 75, 60, 45, 35, 25, 20, 15, 12)) {
  arms <- list(none      = .p(ACE_TGT = 0,  MRA_TGT = 0),
               acei      = .p(ACE_TGT = 20, MRA_TGT = 0),
               acei_mra  = .p(ACE_TGT = 20, MRA_TGT = 25),
               plus_bind = .p(ACE_TGT = 20, MRA_TGT = 25, PAT_RATE = 16.8))
  do.call(rbind, lapply(names(arms), function(a) {
    do.call(rbind, lapply(gfrs, function(g) {
      r <- sim_steady(mod, c(arms[[a]], list(GFR_SET = g)))
      data.frame(arm = a, eGFR = g, K = r$K, FEK = r$FEK, UK = r$UK,
                 COL = r$COLK_EX, SCAP = r$SCAP, ALDO = r$ALDO)
    }))
  }))
}

## ---------------------------------------------------------------------------
## SCENARIO 4.  Binder dose-response, against OPAL-HK and HARMONIZE.
## ---------------------------------------------------------------------------
scenario_binder <- function(mod, gfr = 17) {
  base <- sim_steady(mod, .p(GFR_SET = gfr, ACE_TGT = 20, MRA_TGT = 25))
  grid <- rbind(
    data.frame(drug = "patiromer", dose = c(8.4, 16.8, 25.2, 33.6)),
    data.frame(drug = "SZC",       dose = c(5, 10, 15, 30)))
  grid$K <- NA_real_; grid$PHI <- NA_real_; grid$dK <- NA_real_
  for (i in seq_len(nrow(grid))) {
    pp <- .p(GFR_SET = gfr, ACE_TGT = 20, MRA_TGT = 25)
    if (grid$drug[i] == "patiromer") pp$PAT_RATE <- grid$dose[i] else pp$SZC_RATE <- grid$dose[i]
    r <- sim_steady(mod, pp)
    grid$K[i] <- r$K; grid$PHI[i] <- r$PHI; grid$dK[i] <- r$K - base$K
  }
  attr(grid, "baseline_K") <- base$K
  grid
}

## ---------------------------------------------------------------------------
## SCENARIO 5.  Acute rescue -- the three jobs.
##   Calcium moves the threshold, insulin and salbutamol move the partition,
##   only the binder moves the pool.  Watch KREM, not K.
## ---------------------------------------------------------------------------
scenario_acute <- function(mod, gfr = 12, target_K = 6.8) {
  base <- mod %>% param(GFR_SET = gfr, ACE_TGT = 20, MRA_TGT = 25, HCO3_SET = 19)
  ss   <- sim_steady(base)
  # load the patient up to target_K by adding potassium to the body
  buf  <- ss$BUF_CHR
  dK   <- (target_K - ss$K) * buf
  init0 <- list(KE = ss$KE + dK*(0.20*70)/buf,
                KIF = ss$KIF + dK*(1 - (0.20*70)/buf)*0.25,
                KIS = ss$KIS + dK*(1 - (0.20*70)/buf)*0.75)

  run <- function(lab, ev = NULL, pars = list()) {
    m <- base %>% param(pars) %>% init(init0)
    o <- if (is.null(ev)) mrgsim(m, end = 2, delta = 1/288)
         else mrgsim(m, ev, end = 2, delta = 1/288)
    as.data.frame(o) %>% dplyr::mutate(arm = lab)
  }
  rbind(
    run("no treatment"),
    run("calcium gluconate 1 g", ev(amt = 0.25, cmt = "CAE")),
    run("insulin 10 U + D50 25 g",
        ev(amt = 10*1000/12, cmt = "INS") + ev(amt = 25/180.15*1000/16, cmt = "GLU")),
    run("salbutamol 20 mg neb", ev(amt = 20*0.20/150, cmt = "B2C")),
    run("SZC 10 g TID", NULL, .p(SZC_RATE = 30)),
    run("insulin/D50 + SZC 10 g TID",
        ev(amt = 10*1000/12, cmt = "INS") + ev(amt = 25/180.15*1000/16, cmt = "GLU"),
        .p(SZC_RATE = 30))
  )
}

## ---------------------------------------------------------------------------
## SCENARIO 6.  The RAASi-potassium dilemma with a clinician in the loop.
##   TITRATE = 1 lets the prescribed dose fall when K breaches 5.5 and creep
##   back up when it is safe.  Read the RD column, not the K column: every arm
##   ends with an acceptable potassium, and they differ only in what it cost.
## ---------------------------------------------------------------------------
scenario_dilemma <- function(mod, years = 5) {
  arms <- list(
    "ACEi+MRA, no binder"        = .p(),
    "ACEi+MRA + patiromer"       = .p(PAT_RATE = 16.8),
    "ACEi+MRA + SZC"             = .p(SZC_RATE = 10),
    "ACEi+MRA, low-K diet"       = .p(INTAKE   = 50),
    "ACEi+MRA + furosemide"      = .p(FUR_RATE = 40),
    "ACEi+MRA + SGLT2i"          = .p(SGLT2I   = 1),
    "ACEi+MRA + NaHCO3"          = .p(BIC_RATE = 70),
    "ACEi alone, no MRA"         = .p(MRA_TGT  = 0))
  do.call(rbind, lapply(names(arms), function(a) {
    pp <- c(arms[[a]], .p(GFR_SET = 25, ACE_TGT = 20, MRA_TGT = 25,
                          PROGRESS = 1, TITRATE = 1))
    if (!is.null(arms[[a]]$MRA_TGT)) pp$MRA_TGT <- arms[[a]]$MRA_TGT
    o <- mod %>% param(pp) %>%
      mrgsim(end = years*365, delta = 7, atol = 1e-8, rtol = 1e-8) %>%
      as.data.frame()
    o$arm <- a
    o
  }))
}

## ---------------------------------------------------------------------------
## SCENARIO 7.  Threshold eGFR -- at what renal function does each regimen
##   cross K 5.5?  This is the number that turns "spironolactone raises K by
##   0.3" into a clinical decision.
## ---------------------------------------------------------------------------
scenario_threshold <- function(mod, arms = NULL) {
  if (is.null(arms)) arms <- list(
    none      = .p(ACE_TGT = 0,  MRA_TGT = 0),
    acei      = .p(ACE_TGT = 20, MRA_TGT = 0),
    acei_mra  = .p(ACE_TGT = 20, MRA_TGT = 25),
    plus_bind = .p(ACE_TGT = 20, MRA_TGT = 25, PAT_RATE = 16.8))
  bisect <- function(pp) {
    lo <- 5; hi <- 120
    for (i in 1:40) {
      mid <- (lo + hi)/2
      k <- sim_steady(mod, c(pp, list(GFR_SET = mid)))$K
      if (k > 5.5) lo <- mid else hi <- mid
    }
    (lo + hi)/2
  }
  data.frame(arm = names(arms),
             threshold_eGFR = vapply(arms, bisect, numeric(1)))
}

## ---------------------------------------------------------------------------
## SCENARIO 8.  The dietary ceiling -- the maximum intake compatible with
##   K <= 5.5, and the hard bound on what any binder can do.
## ---------------------------------------------------------------------------
scenario_diet <- function(mod, gfrs = c(60, 45, 30, 20, 15)) {
  arms <- list(none = .p(ACE_TGT = 0, MRA_TGT = 0),
               acei_mra = .p(ACE_TGT = 20, MRA_TGT = 25),
               plus_szc = .p(ACE_TGT = 20, MRA_TGT = 25, SZC_RATE = 10))
  do.call(rbind, lapply(names(arms), function(a) {
    do.call(rbind, lapply(gfrs, function(g) {
      lo <- 5; hi <- 400
      for (i in 1:40) {
        mid <- (lo + hi)/2
        k <- sim_steady(mod, c(arms[[a]], list(GFR_SET = g, INTAKE = mid)))$K
        if (k < 5.5) lo <- mid else hi <- mid
      }
      data.frame(arm = a, eGFR = g, max_intake = (lo + hi)/2)
    }))
  }))
}

## ---------------------------------------------------------------------------
## SCENARIO 9.  Acid-base -- alkali is a RENAL potassium therapy.
## ---------------------------------------------------------------------------
scenario_alkali <- function(mod, gfr = 20) {
  do.call(rbind, lapply(c(16, 18, 20, 22, 24), function(h) {
    r <- sim_steady(mod, .p(GFR_SET = gfr, ACE_TGT = 20, HCO3_SET = h))
    data.frame(HCO3 = h, pH = r$PHA, LAMrel = r$LAMrel, K = r$K,
               UK = r$UK, LAM_route = "negligible", renal_route = r$UK)
  }))
}

## ---------------------------------------------------------------------------
## SCENARIO 10.  The two patients with the same number -- mass-balance vs
##   partition hyperkalaemia.  This is the scenario to run before believing a
##   potassium result.
## ---------------------------------------------------------------------------
scenario_two_patients <- function(mod) {
  normal <- sim_steady(mod, .p(GFR_SET = 100, ACE_TGT = 0, MRA_TGT = 0))
  A <- sim_steady(mod, .p(GFR_SET = 12, ACE_TGT = 20, MRA_TGT = 25, HCO3_SET = 19))
  # patient B: DKA -- collapse LAMrel via pH, glucose and insulinopenia, and
  # carry the whole-body deficit DKA always has
  B <- mod %>%
    param(GFR_SET = 60, ACE_TGT = 0, MRA_TGT = 0, HCO3_SET = 5,
          INS_BASE = 0.5, GLU_BASE = 40) %>%
    init(KE = (normal$KE + normal$KIF + normal$KIS - 400) * 0.016) %>%
    mrgsim(end = 0.02, delta = 0.001) %>% as.data.frame()
  list(normal = normal, mass_balance_patient = A, partition_patient = B)
}

## ---------------------------------------------------------------------------
## SCENARIO 11.  Virtual population -- who a binder actually rescues.
## ---------------------------------------------------------------------------
scenario_population <- function(mod, n_gfr = 10, n_diet = 10, n_hco3 = 6) {
  grid <- expand.grid(GFR_SET = seq(15, 60, length.out = n_gfr),
                      INTAKE  = seq(40, 140, length.out = n_diet),
                      HCO3_SET= seq(16, 26, length.out = n_hco3))
  grid$K_nobinder <- NA_real_; grid$K_binder <- NA_real_
  for (i in seq_len(nrow(grid))) {
    pp <- as.list(grid[i, c("GFR_SET", "INTAKE", "HCO3_SET")])
    grid$K_nobinder[i] <- sim_steady(mod, c(pp, .p(ACE_TGT = 20, MRA_TGT = 25)))$K
    grid$K_binder[i]   <- sim_steady(mod, c(pp, .p(ACE_TGT = 20, MRA_TGT = 25,
                                                   PAT_RATE = 16.8)))$K
  }
  grid
}

## ---------------------------------------------------------------------------
## SCENARIO 12.  ECG / membrane arithmetic across the potassium range, with
##   and without calcium.  Calcium restores ~29% of the lost excitability gap
##   and removes exactly nothing.
## ---------------------------------------------------------------------------
scenario_ecg <- function(mod, ks = seq(3.0, 8.0, by = 0.25)) {
  do.call(rbind, lapply(ks, function(k) {
    m <- mod %>% init(KE = k*0.20*70, KIF = 140*0.25*0.36*70, KIS = 140*0.75*0.36*70)
    o <- as.data.frame(mrgsim(m, end = 0, delta = 1))
    data.frame(K = k, EM = o$EM[1], GAP = o$GAP[1], QRS = o$QRS[1],
               TWAVE = o$TWAVE[1], HZK = o$HZK[1])
  }))
}

## ---------------------------------------------------------------------------
## Run everything.
## ---------------------------------------------------------------------------
scenario_all <- function(mod) {
  list(reserve      = scenario_reserve(mod),
       binder       = scenario_binder(mod),
       acute        = scenario_acute(mod),
       dilemma      = scenario_dilemma(mod),
       threshold    = scenario_threshold(mod),
       diet         = scenario_diet(mod),
       alkali       = scenario_alkali(mod),
       two_patients = scenario_two_patients(mod),
       ecg          = scenario_ecg(mod))
}
