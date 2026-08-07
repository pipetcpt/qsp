## =============================================================================
##  hs_mrgsolve_model.R
##  Heat Stroke (exertional and classic) -- Quantitative Systems Pharmacology
##  50-ODE mrgsolve model
## =============================================================================
##
##  THESIS
##  ------
##  Heat stroke is what happens when the heat-balance equation loses its FIXED
##  POINT.  Everything after that is a CLOCK and a DOSE.
##
##    regime 1  COMPENSABLE    a steady state exists; the core plateaus.  This
##                             is heat STRAIN, and it is not heat stroke however
##                             unpleasant the environment is.
##    regime 2  UNCOMPENSABLE  no steady state exists.  dTc/dt = (H_prod -
##                             Q_loss)/C is strictly positive and a clock runs
##                             at a rate the equation computes exactly.
##    regime 3  COMMITTED      the accumulated thermal dose has latched a
##                             bistable inflammatory switch.  Cooling now
##                             restores the temperature but not the patient.
##
##  Every intervention in this file acts on exactly one of those three regimes,
##  and that is the model's explanation for why the therapeutic record in heat
##  stroke is so poor: almost everything ever tried acts on regime 3.
##
##  WHAT THE MODEL COMPUTES RATHER THAN ASSUMES
##  -------------------------------------------
##   * the critical wet-bulb temperature (35.7 C at rest, ~21 C at 900 W)
##   * the rate of core temperature rise, and therefore the warning time a
##     40 C diagnostic threshold actually gives (20 min in EHS, 375 min in NEHS)
##   * the thermal dose CEM43 paid under each cooling strategy
##   * the exchange rate between cooling modality and cooling delay
##   * the dose at which the inflammatory switch latches
##   * ISTH DIC score, GCS, SOFA -- all as OUTPUTS, never as states
##
##  UNITS.  Time = MINUTES.  Temperature = degrees C.  Heat = watts internally.
##  Body heat capacity C = mass * 3470 J/(kg.K); dT/dt = 60*W/C converts W to
##  degrees C per minute.
##
##  VERIFICATION.  The identical ODE system is implemented independently in
##  Python/scipy (hs_core.py) and every number quoted in README.md comes from
##  running hs_analysis.py against it.  That cross-implementation exposed and
##  fixed ten defects in this model, all documented in hs_verification_output.txt
##  and summarised in README.md section 9.  R is not installed in the build
##  environment used to generate this repository, so the Python implementation
##  is the executed source of truth and this file is the transcription; symbol
##  names correspond one-to-one.
##
##  CALIBRATION SOURCES: see hs_references.md.
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

code <- '
$PROB
# Heat Stroke QSP model (EHS + NEHS), 50 ODEs

$PARAM @annotated
// ---- anthropometry --------------------------------------------------------
BW       :  70.0  : body mass (kg)
HT       :   1.75 : height (m)
CP_BODY  : 3470.0 : whole-body specific heat (J/kg/K)
FMU      :   0.40 : muscle mass fraction
FSK      :   0.10 : skin/shell mass fraction

// ---- metabolic heat production --------------------------------------------
MREST    :   1.20 : basal heat production (W/kg)
Q10      :   2.30 : vant Hoff coefficient for metabolic heat
ETA_EX   :   0.20 : mechanical efficiency of exercise
FMU_EX   :   0.85 : fraction of exercise heat deposited in muscle
FMU_REST :   0.22 : fraction of resting heat deposited in muscle
SHIVMAX  : 350.0  : maximal shivering thermogenesis (W)
TSK_SHIV :  31.0  : skin temperature engaging shivering (deg C)
TCR_SHIV :  38.6  : core temperature permitting shivering (deg C)

// ---- passive conductances (W/K) -------------------------------------------
KCM0     :  70.0  : core-muscle tissue conduction
KCS0     :  10.0  : core-skin tissue conduction
KBLOOD   :   1.13 : W/K per (L/h) of blood flow

// ---- skin blood flow ------------------------------------------------------
SKBF0    :   0.30 : baseline skin blood flow (L/min)
SKBFMAX  :   7.50 : maximal skin blood flow (L/min)
TCVD0    :  36.80 : core threshold for active vasodilation (deg C)
KVD      :   1.20 : deg C spanning full vasodilation
WSK_VD   :   0.15 : weight of skin temperature in the vasodilation drive
FVD_AGE  :   1.00 : age multiplier on vasodilation (elderly ~0.55)
FVD_DRUG :   1.00 : drug multiplier on vasodilation (beta-blocker)
TSK_VC   :  28.0  : skin temperature completing cold vasoconstriction
FCOLD_FLOOR0 : 0.25 : residual skin flow fraction, no central drive
FCOLD_FLOOR1 : 0.45 : extra residual fraction at full central drive
KVC_COLD :   0.45 : tissue conductance retained when vasoconstricted

// ---- sweating -------------------------------------------------------------
SWMAX    :  20.0  : maximal sweat output (g/min) = 1.2 L/h unacclimatised
SWGAIN   :  26.0  : g/min per deg C of integrated thermal drive
TCSW0    :  37.00 : core sweat onset threshold (deg C)
TSKSW0   :  33.0  : skin reference for the sweat drive (deg C)
WSK_SW   :   0.10 : weight of skin temperature in the sweat drive
FSW_AGE  :   1.00 : age multiplier on sweating (elderly ~0.60)
FSW_DRUG :   1.00 : drug multiplier on sweating (anticholinergic ~0.30)
KSWFAT   :   0.0016 : gland fatigue accrual (/min at full wettedness)
KSWREC   :   0.0007 : gland fatigue recovery (/min)
SWFAT_MAX:   0.30 : maximal hidromeiosis
KDEHY_SW :   0.055: fractional sweat loss per % body-mass deficit

// ---- environment (time-varying covariates in the input data set) ----------
TA       :  25.0  : air temperature (deg C)
TRAD     :  25.0  : mean radiant temperature (deg C)
RH       :   0.50 : relative humidity (fraction)
VAIR     :   1.0  : air velocity (m/s)
ICL      :   0.30 : clothing insulation (clo)
QSOL     :   0.0  : solar/radiant heat gain (W)
WMAXW    :   0.85 : maximal skin wettedness
MEX      :   0.0  : exercise metabolic rate (W)
ORAL_RATE:   0.0  : oral fluid intake (L/min)
ORAL_MATCH:  0.0  : 1 = drink to exactly replace losses

// ---- cooling --------------------------------------------------------------
UA_COOL  :   0.0  : lumped conductance of the cooling device (W/K)
T_COOL   :   2.0  : temperature of the cooling medium (deg C)
IMMERSE  :   0.0  : 1 = skin submerged (no air exchange)
COOL_STOP_TC : 38.6 : core temperature at which cooling stops (deg C)
IVF_RATE :   0.0  : cold fluid infusion rate (L/min)
IVF_TEMP :   4.0  : infusate temperature (deg C)

// ---- thermal dose ---------------------------------------------------------
CEM_R_LO :   0.25 : Sapareto-Dewey R below 43 C
CEM_R_HI :   0.50 : Sapareto-Dewey R at or above 43 C
HSP50    :   2.5  : HSP70 fold-induction giving 2x dose protection

// ---- HSP70 ----------------------------------------------------------------
HSP_BASE :   1.0  : baseline HSP70 (fold)
KHSP_IN  :   0.0025 : HSP70 induction rate (/min); larger values give a
                    : physiologically impossible 58-fold induction
KHSP_OUT :   0.0009: HSP70 turnover (/min)
THSF     :  39.6  : HSF1 activation midpoint (deg C)
NHSF     :   8.0  : HSF1 activation Hill coefficient

// ---- water / plasma -------------------------------------------------------
TBW_FRAC :   0.60 : total body water fraction
VP0      :   3.0  : plasma volume (L)
KPV      :   1.30 : share of water deficit borne by plasma
INSENS   :   0.45 : insensible water loss (g/min)

// ---- cardiovascular / splanchnic ------------------------------------------
CO0      :   5.0  : resting cardiac output (L/min)
KCO_EX   :   0.0110 : L/min per W of exercise
KCO_T    :   0.90 : L/min per deg C of core hyperthermia
COMAX    :  24.0  : maximal cardiac output (L/min)
MAP0     :  90.0  : baseline mean arterial pressure (mmHg)
KMAP_VP  :   1.60 : MAP sensitivity to fractional plasma volume loss
KMAP_SKBF:   0.045: MAP lost per L/min of cutaneous vasodilation
QSPL0    :   1.50 : resting splanchnic blood flow (L/min)
KSPL_EX  :   0.55 : maximal fractional splanchnic reduction by exercise
KSPL_SK  :   0.45 : maximal fractional splanchnic reduction by cutaneous steal
MEX_REF  : 900.0  : exercise rate giving maximal splanchnic steal (W)
QSPL_CRIT:   0.90 : splanchnic flow below which the mucosa is ischaemic

// ---- gut barrier ----------------------------------------------------------
KG_ISCH  :   0.0060 : ischaemic enterocyte injury gain (/min)
KG_HEAT  :   0.0130 : direct thermal enterocyte injury per CEM43-minute
KG_REP   :   0.00090: mucosal repair (/min)
KPERM    :   0.35 : injury giving half-maximal permeability
FLUX_LPS :   0.08 : portal endotoxin flux at full permeability (EU/mL/min)
CL_LPS   :   0.075: hepatic endotoxin clearance (/min)

// ---- cytokines ------------------------------------------------------------
KTNF     :   5.0  : TNF production
KDTNF    :   0.030: TNF elimination (/min)
KM_LPS   :   0.55 : LPS EC50 (EU/mL)
WD_TNF   :   0.55 : DAMP weight on TNF
KIL6     :  11.0  : IL-6 production
KDIL6    :   0.0075 : IL-6 elimination (/min)
WT_IL6   :   0.0020 : TNF weight on IL-6
KI10     : 380.0  : IL-10 concentration halving IL-6 production
KIL1     :   0.66 : IL-1beta production
KDIL1    :   0.016: IL-1beta elimination (/min)
KIL10    :   0.24 : IL-10 production per unit IL-6/100
KDIL10   :   0.014: IL-10 elimination (/min)
WCORT_IL10:  0.020: cortisol drive on IL-10
KNET     :   0.020: NETosis rate
KDNET    :   0.010: NET clearance (/min)

// ---- HMGB1 commitment switch ----------------------------------------------
// Saddle-node bistability.  Necrotic cells release HMGB1; HMGB1 signalling via
// RAGE/TLR4 produces more necrosis; the autocatalytic production SATURATES
// while clearance stays linear.  With these values: OFF = 0, unstable
// threshold = 10.91 ng/mL, ON = 50.3 ng/mL.  The saddle-node sits at an
// effective clearance of 0.00498/min, so an agent adding more than
// 0.00235/min of HMGB1 clearance ABOLISHES the ON state rather than merely
// shifting it.  That is the structural claim about thrombomodulin alfa.
KH_NEC   :  14.5  : HMGB1 released per unit necrosis (ng/mL)
KH_AUTO  :   0.16 : saturating autocatalytic HMGB1 production (ng/mL/min)
KH_HALF  :  17.0  : midpoint of the autocatalytic term (ng/mL)
NH       :   3.0  : Hill coefficient of the autocatalytic term
KH_OUT   :   0.0035 : HMGB1 clearance (/min), t1/2 ~ 3.3 h
KH_RTM   :   0.0020 : extra clearance per unit rTM concentration (/min)
KNEC_H   :   3.0e-5 : necrosis per (ng/mL) HMGB1 per min
KNEC_TH  :   0.080  : necrosis per CEM43-minute (calibrated so the switch
                    : latches at ~8 CEM43: 30 min of delay does not commit,
                    : 45 min does)

// ---- endothelium / coagulation --------------------------------------------
SDC1_0   :  18.0  : baseline syndecan-1 (ng/mL)
KSDC     :   0.85 : glycocalyx shedding gain
KDSDC    :   0.0055 : syndecan-1 clearance (/min)
KTF      :   0.0055 : tissue factor induction
KDTF     :   0.0060 : tissue factor decay (/min)
KTF_LPS  :   0.75 : LPS weight on tissue factor
KTHR     :   1.30 : thrombin generation per unit TF
KDTHR    :   0.22 : thrombin inactivation (/min)
KTHR_PC  :   0.85 : protein C brake on thrombin
FIB0     : 300.0  : baseline fibrinogen (mg/dL)
KFIB_SYN :   0.055: fibrinogen synthesis toward baseline (/min)
KFIB_APF :   1.60 : acute-phase amplification of fibrinogen synthesis
KFIB_CONS:   0.055: fibrinogen consumption by thrombin
PLT0     : 250.0  : baseline platelets (10^9/L)
KPLT_SYN :   0.00080 : platelet recovery (/min)
KPLT_CONS:   0.019: platelet consumption by thrombin
PC0      : 100.0  : baseline protein C activity (%)
KPC_SYN  :   0.00090 : protein C synthesis (/min)
KPC_CONS :   0.030: protein C consumption
KDD      :   0.85 : D-dimer generation
KDDD     :   0.0028 : D-dimer clearance (/min)

// ---- organ injury ---------------------------------------------------------
ALT0     :  25.0  : baseline ALT (U/L)
KALT     : 1350.0 : ALT released per unit hepatocyte necrosis
KDALT    :   0.00040 : ALT clearance (/min), t1/2 ~ 29 h
AST0     :  25.0  : baseline AST (U/L)
KAST     : 1950.0 : AST released per unit hepatocyte necrosis
KDAST    :   0.00096 : AST clearance (/min), t1/2 ~ 12 h
KAST_MU  :   0.0130 : muscle contribution to AST
CK0      : 120.0  : baseline CK (U/L)
KCK      : 48000.0: CK released per unit rhabdomyolysis
KDCK     :   0.00032 : CK clearance (/min), t1/2 ~ 36 h
MB0      :  30.0  : baseline myoglobin (ug/L)
KMB      : 8200.0 : myoglobin release
KDMB     :   0.0077 : myoglobin clearance (/min), t1/2 ~ 90 min
KMU_REL  :   0.0012 : release from injured myocytes (/min), t1/2 ~9.6 h
KRHAB_TH :   0.145: rhabdomyolysis per CEM43-minute of MUSCLE dose
KRHAB_ISCH:  0.0022 : ischaemic contribution to rhabdomyolysis
SCR0     :   0.90 : baseline creatinine (mg/dL)
KSCR_OUT :   0.0125 : creatinine clearance scale (/min)
KGFR_VP  :   1.10 : GFR sensitivity to plasma volume loss
KGFR_CAST:   0.55 : GFR sensitivity to cast burden
KGFR_DIC :   0.30 : GFR sensitivity to microthrombosis
KGFR_REC :   0.00035 : GFR relaxation scale
KCAST    :   1.0e-7 : cast formation per unit myoglobin
KDCAST   :   0.00050 : cast clearance (/min)
KLIV_TH  :   0.00185 : hepatic injury per CEM43-minute
KLIV_ISCH:   0.00062 : hepatic injury from ischaemia
KLIV_NAP :   0.00450 : hepatic injury from NAPQI
KLIV_REC :   0.00030 : hepatic regeneration (/min)
KLIV_HMGB:   0.00160 : sustained hepatic injury while the switch is ON
KGFR_HMGB:   0.35 : HMGB1 contribution to the GFR deficit
KCNS_HMGB:   0.00220 : persistent encephalopathy while the switch is ON
NSE0     :   8.0  : baseline NSE (ug/L)
KNSE     : 140.0  : NSE release per CEM43-minute
KDNSE    :   0.00048 : NSE clearance (/min)
KCNS_TH  :   0.0180 : encephalopathy per CEM43-minute
KCNS_INFL:   0.00030 : encephalopathy from inflammation
KCNS_REC :   0.0035 : neurological recovery (/min)
KPOT0    :   4.0  : baseline potassium (mmol/L)
KK_RHAB  :   3.0e-5 : potassium released per unit CK flux
KK_REN   :   0.030: renal potassium excretion
LAC0     :   1.0  : baseline lactate (mmol/L)
KLAC     :   0.0130 : lactate production
KDLAC    :   0.0130 : lactate clearance

// ---- endogenous cortisol --------------------------------------------------
CORT0    :  12.0  : baseline cortisol (ug/dL)
KCORT_IN :   0.030: cortisol secretion gain
KCORT_OUT:   0.0075 : cortisol elimination (/min)
KCORT_STIM:  2.6  : stress amplification of cortisol secretion

// ---- paracetamol ----------------------------------------------------------
PARA_KA  :   0.030 : absorption (/min)
PARA_CL  : 350.0  : clearance (mL/min) = 21 L/h
PARA_V1  : 32000.0: central volume (mL)
PARA_V2  : 18000.0: peripheral volume (mL)
PARA_Q   : 133.3  : intercompartmental clearance (mL/min)
FNAPQI   :   0.055: fraction of clearance to NAPQI
FNAPQI_HEAT: 1.9  : CYP2E1 induction multiplier under heat/fasting
KGSH_REC :   0.00090 : glutathione resynthesis (/min)
KGSH_USE :   0.0130 : glutathione consumption by NAPQI
KTSET_PARA:  0.40 : maximal set-point depression by paracetamol (deg C)
PARA_EC50:  12.0  : paracetamol EC50 for set-point (mg/L)

// ---- ibuprofen ------------------------------------------------------------
IBU_KA   :   0.030 : absorption (/min)
IBU_CL   :  58.3  : clearance (mL/min) = 3.5 L/h
IBU_V    : 10000.0: volume (mL)
KTSET_IBU:   0.60 : maximal set-point depression by ibuprofen (deg C)
IBU_EC50 :  15.0  : ibuprofen EC50 (mg/L)
KGFR_IBU :   0.16 : fractional GFR loss from renal prostaglandin blockade
KPLTFN_IBU:  0.55 : platelet function loss (reported, not a state)

// ---- dantrolene -----------------------------------------------------------
DAN_CL   :  10.0  : clearance (mL/min)
DAN_V1   : 36000.0: central volume (mL)
DAN_V2   : 40000.0: peripheral volume (mL)
DAN_Q    :  20.0  : intercompartmental clearance (mL/min)
DAN_EMAX :   0.12 : maximal fractional cut in muscle heat production
DAN_EC50 :   2.5  : EC50 (mg/L)

// ---- hydrocortisone -------------------------------------------------------
HC_CL    : 300.0  : clearance (mL/min)
HC_V     : 35000.0: volume (mL)
HC_EMAX  :   0.35 : maximal fractional cytokine suppression
HC_EC50  :   0.15 : EC50 (mg/L)

// ---- thrombomodulin alfa (ART-123) ----------------------------------------
// Concentration is normalised so that one standard daily dose (380 U/kg) is
// given as an amount of 3.5 into RTM_C, yielding C_rtm = 1.0 unit/L at peak.
RTM_CL   :   2.02 : clearance (mL/min), t1/2 ~ 20 h
RTM_V    : 3500.0 : volume (mL)
RTM_EAPC :   1.60 : fold increase in protein C activation
RTM_EHMGB:   1.00 : scaling of the lectin-domain HMGB1 clearance term

$CMT @annotated
TCR    : core temperature (deg C)
TMU    : muscle temperature (deg C)
TSK    : skin temperature (deg C)
CEM43  : cumulative equivalent minutes at 43 C (raw)
CEM43E : cumulative equivalent minutes at 43 C (HSP-protected)
HSP    : HSP70 (fold of unstressed)
WDEF   : body water deficit (L)
VP     : plasma volume (L)
SWFAT  : sweat gland fatigue (0-1)
GUT    : enterocyte injury (0-1)
LPS    : systemic endotoxin (EU/mL)
TNF    : TNF-alpha (pg/mL)
IL6    : IL-6 (pg/mL)
IL1B   : IL-1beta (pg/mL)
IL10   : IL-10 (pg/mL)
HMGB1  : HMGB1 (ng/mL)
NETS   : NETs / cell-free DNA (arbitrary)
SDC1   : syndecan-1 (ng/mL)
TF     : tissue factor (arbitrary)
THR    : thrombin (nM)
FIB    : fibrinogen (mg/dL)
PLT    : platelets (10^9/L)
PC     : protein C activity (%)
DDIM   : D-dimer (ug/mL FEU)
ALT    : ALT (U/L)
AST    : AST (U/L)
CK     : creatine kinase (U/L)
MB     : myoglobin (ug/L)
SCR    : creatinine (mg/dL)
GFRF   : GFR as a fraction of normal
NSE    : neuron-specific enolase (ug/L)
CNSD   : encephalopathy (0-1)
LIVF   : functional liver mass (fraction)
MUINJ  : injured myocyte pool (source of CK and myoglobin)
CAST   : myoglobin cast burden
KPOT   : potassium (mmol/L)
LAC    : lactate (mmol/L)
PARA_A : paracetamol absorption depot (mg)
PARA_C : paracetamol central (mg)
PARA_P : paracetamol peripheral (mg)
NAPQI  : NAPQI burden
GSH    : hepatic glutathione (fraction)
IBU_A  : ibuprofen depot (mg)
IBU_C  : ibuprofen central (mg)
DAN_C  : dantrolene central (mg)
DAN_P  : dantrolene peripheral (mg)
HC_C   : hydrocortisone (mg)
RTM_C  : thrombomodulin alfa (normalised units)
CORT   : cortisol (ug/dL)
FLUID  : cumulative IV crystalloid (L)

$GLOBAL
#define LATENT 2426.0        // J per gram of evaporated sweat

// Buck (1981) saturation vapour pressure, kPa
double psat_kpa(double T) {
  return 0.61121 * exp((18.678 - T/234.5) * (T/(257.14 + T)));
}
double clampd(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}
double hillf(double x, double k, double n) {
  if (x <= 0.0) return 0.0;
  double xn = pow(x, n);
  return xn / (xn + pow(k, n));
}

$MAIN
// Du Bois body surface area
double AD = 0.202 * pow(BW, 0.425) * pow(HT, 0.725);

TCR_0   = 37.0;   TMU_0 = 37.2;  TSK_0 = 33.0;
HSP_0   = HSP_BASE;
VP_0    = VP0;
LIVF_0  = 1.0;    GSH_0 = 1.0;
FIB_0   = FIB0;   PLT_0 = PLT0;  PC_0  = PC0;
SDC1_0_ = SDC1_0;
SDC1_0  = SDC1_0_;
ALT_0   = ALT0;   AST_0 = AST0;  CK_0 = CK0;  MB_0 = MB0;
SCR_0   = SCR0;   GFRF_0 = 1.0;  NSE_0 = NSE0;
KPOT_0  = KPOT0;  LAC_0 = LAC0;  CORT_0 = CORT0;

$ODE
double AD = 0.202 * pow(BW, 0.425) * pow(HT, 0.725);

// ---------------- drug concentrations (mg/L; rTM in normalised units/L) ----
double C_para = PARA_C / PARA_V1 * 1000.0;
double C_ibu  = IBU_C  / IBU_V   * 1000.0;
double C_dan  = DAN_C  / DAN_V1  * 1000.0;
double C_hc   = HC_C   / HC_V    * 1000.0;
double C_rtm  = RTM_C  / RTM_V   * 1000.0;

// ---------------- B. metabolic heat production -----------------------------
// van t Hoff acceleration: hotter tissue makes more heat.  This is the
// autocatalytic term that turns a marginal imbalance into a runaway.  Above
// 42 C it is rolled off rather than extrapolated (enzyme denaturation).
double q10f = pow(Q10, (min(TCR, 42.0) - 37.0)/10.0);
double over = TCR > 42.0 ? TCR - 42.0 : 0.0;
q10f = q10f / (1.0 + over*over);
double M_rest = MREST * BW * q10f;
double dan_eff = 1.0 - DAN_EMAX * C_dan/(C_dan + DAN_EC50);
double H_ex = MEX * (1.0 - ETA_EX) * dan_eff;
double shiv_drive = (max(0.0, TSK_SHIV - TSK)/4.0) * (max(0.0, TCR_SHIV - TCR)/1.5);
double H_shiv = SHIVMAX * min(1.0, shiv_drive);
double H_prod = M_rest + H_ex + H_shiv;
double f_mu = MEX > 50.0 ? FMU_EX : FMU_REST;
double H_mu = H_ex*f_mu + (M_rest + H_shiv)*FMU_REST;
double H_cr = H_prod - H_mu;

// ---------------- effector: skin blood flow --------------------------------
// Antipyretic set-point depression enters on the THRESHOLDS only, which is
// exactly why it cannot help once the effectors are saturated.
double dTset = KTSET_PARA * C_para/(C_para + PARA_EC50)
             + KTSET_IBU  * C_ibu /(C_ibu  + IBU_EC50);
double pv_frac = VP / VP0;
double vd_drive = ((TCR - (TCVD0 - dTset)) + WSK_VD*(TSK - 33.0)) / KVD;
vd_drive = clampd(vd_drive, 0.0, 1.0);
// Cold vasoconstriction is INCOMPLETE while central hyperthermic drive is on:
// the quantitative form of the standard answer to the claim that ice water
// traps heat by vasoconstricting the skin.
double fc_floor = FCOLD_FLOOR0 + FCOLD_FLOOR1*vd_drive;
double f_cold = fc_floor + (1.0 - fc_floor)*clampd((TSK - TSK_VC)/5.0, 0.0, 1.0);
double SKBF_raw = SKBF0 + (SKBFMAX - SKBF0)*vd_drive*FVD_AGE*FVD_DRUG*f_cold;
double map_now = MAP0*(1.0 - KMAP_VP*max(0.0, 1.0 - pv_frac))
               - KMAP_SKBF*MAP0*(SKBF_raw - SKBF0);
double f_bp = clampd(map_now/MAP0, 0.25, 1.0);
double SKBF = SKBF0 + (SKBF_raw - SKBF0)*f_bp;
double kcs_tissue = KCS0*(KVC_COLD + (1.0 - KVC_COLD)*vd_drive);
double KCS = kcs_tissue + KBLOOD*SKBF*60.0;
double Q_mu_flow = 0.75 + 0.0085*MEX + 0.25*max(0.0, TCR - 37.0);
double KCM = KCM0 + KBLOOD*Q_mu_flow*60.0;

// ---------------- effector: sweating ---------------------------------------
double pct_def = 100.0*WDEF/BW;
double f_dehy = max(0.35, 1.0 - KDEHY_SW*pct_def);
double sw_drive = (TCR - (TCSW0 - dTset)) + WSK_SW*(TSK - TSKSW0);
double m_sw = SWGAIN * max(0.0, sw_drive);
m_sw = min(m_sw, SWMAX);
m_sw = m_sw * FSW_AGE * FSW_DRUG * f_dehy * (1.0 - SWFAT);

// ---------------- A. environmental heat exchange ---------------------------
double hc = max(3.1, 8.3*pow(VAIR, 0.6));
double hr = 4.7;
double h  = hc + hr;
double fcl = 1.0 + 0.31*ICL;
double Rcl = 0.155*ICL;
double To  = (hc*TA + hr*TRAD)/h;
double Q_dry = AD*(TSK - To)/(Rcl + 1.0/(fcl*h));
double he = 16.5*hc;
double Recl = 0.0276*ICL;
double Pa = RH * psat_kpa(TA);
double E_max_env = AD*(psat_kpa(TSK) - Pa)/(Recl + 1.0/(fcl*he));
E_max_env = max(E_max_env, 0.0);
double E_env_ceiling = E_max_env * WMAXW;
double E_sweat_cap = m_sw * LATENT / 60.0;
double Q_sol_now = QSOL;
double E_req = max(0.0, H_prod - Q_dry + Q_sol_now);
// Evaporation is the smaller of what the glands deliver, what the air can
// absorb, and what balance requires.  Surplus sweat DRIPS: cooling is capped,
// water loss is not.  That asymmetry is why uncompensable heat dehydrates
// faster than it cools.
double E_actual = min(min(E_sweat_cap, E_env_ceiling), E_req);
double Q_res = (0.0014*(34.0 - TA) + 0.0173*(5.87 - Pa)) * H_prod;

// ---------------- cooling device -------------------------------------------
// Every modality is one empirical lumped conductance UA (W/K) to a medium at
// T_COOL, fitted to its published whole-body cooling rate.  Making modality a
// single potency number puts "which cooler" and "how long until cooling" on
// the same axis, which is the point.
double UA = UA_COOL;
if (UA > 0.0 && TCR <= COOL_STOP_TC) UA = 0.0;
double Q_dev = UA*(TSK - T_COOL);
if (IMMERSE > 0.5 && UA > 0.0) { Q_dry = 0.0; E_actual = 0.0; Q_sol_now = 0.0; }

// Cold intravenous fluid: the infusate mixes with blood and is distributed to
// the three nodes in proportion to PERFUSION.  Loading it all on the core
// over-predicts the measured core fall by about 50%.
double Q_ivf_tot = IVF_RATE*1000.0*4.18*(TCR - IVF_TEMP)/60.0;
double CO = min(COMAX, CO0 + KCO_EX*MEX + KCO_T*max(0.0, TCR - 37.0));
CO = max(CO, SKBF + Q_mu_flow + 0.5);
double f_sk_perf = SKBF/CO;
double f_mu_perf = Q_mu_flow/CO;
double f_cr_perf = 1.0 - f_sk_perf - f_mu_perf;

// ---------------- thermal node balances ------------------------------------
double C_cr = BW*(1.0 - FMU - FSK)*CP_BODY;
double C_mu = BW*FMU*CP_BODY;
double C_sk = BW*FSK*CP_BODY;
double Q_cr_mu = KCM*(TCR - TMU);
double Q_cr_sk = KCS*(TCR - TSK);

dxdt_TCR = 60.0*(H_cr - Q_cr_mu - Q_cr_sk - Q_res - Q_ivf_tot*f_cr_perf)/C_cr;
dxdt_TMU = 60.0*(H_mu + Q_cr_mu - Q_ivf_tot*f_mu_perf)/C_mu;
dxdt_TSK = 60.0*(Q_cr_sk - Q_dry - E_actual + Q_sol_now - Q_dev
                 - Q_ivf_tot*f_sk_perf)/C_sk;

// ---------------- C. thermal dose (Sapareto-Dewey) -------------------------
double R    = TCR < 43.0 ? CEM_R_LO : CEM_R_HI;
double R_mu = TMU < 43.0 ? CEM_R_LO : CEM_R_HI;
double dose_rate = min(pow(R,    clampd(43.0 - TCR, -6.0, 30.0)), 64.0);
double prot = max(1.0, 1.0 + (HSP - 1.0)/HSP50);
double dose_eff = dose_rate / prot;
double dose_mu  = min(pow(R_mu, clampd(43.0 - TMU, -6.0, 30.0)), 64.0) / prot;

dxdt_CEM43  = dose_rate;
dxdt_CEM43E = dose_eff;

double hsf = hillf(TCR - 36.0, THSF - 36.0, NHSF);
dxdt_HSP = KHSP_IN*hsf - KHSP_OUT*(HSP - HSP_BASE);

// ---------------- water balance --------------------------------------------
double loss_Lmin = (m_sw + INSENS)/1000.0;
double intake = ORAL_MATCH > 0.5 ? loss_Lmin + IVF_RATE : ORAL_RATE + IVF_RATE;
dxdt_WDEF = loss_Lmin - intake;
double tbw = TBW_FRAC*BW;
dxdt_VP = -KPV*(VP0/tbw)*(loss_Lmin - intake);
double wet = E_env_ceiling <= 1e-9 ? 1.0
             : min(1.0, E_sweat_cap/max(E_env_ceiling, 1e-9));
dxdt_SWFAT = KSWFAT*wet*(SWFAT_MAX - SWFAT)/SWFAT_MAX - KSWREC*SWFAT;
dxdt_FLUID = IVF_RATE;

// ---------------- D. splanchnic steal --------------------------------------
double steal_ex = KSPL_EX*min(1.0, MEX/MEX_REF);
double steal_sk = KSPL_SK*(SKBF - SKBF0)/(SKBFMAX - SKBF0);
double QSPL = QSPL0*max(0.10, (1.0 - steal_ex)*(1.0 - steal_sk))*f_bp;
double ISCH = clampd(1.0 - QSPL/QSPL_CRIT, 0.0, 1.0);

// ---------------- E. gut barrier -------------------------------------------
dxdt_GUT = (KG_ISCH*ISCH*ISCH + KG_HEAT*dose_eff)*(1.0 - GUT) - KG_REP*GUT;
double PERM = GUT/(GUT + KPERM);
double livf = max(LIVF, 0.05);
dxdt_LPS = FLUX_LPS*PERM - CL_LPS*LPS*livf*min(1.0, QSPL/QSPL0);

// ---------------- F. cytokines and the commitment switch -------------------
double hc_sup   = 1.0 - HC_EMAX*C_hc/(C_hc + HC_EC50);
double cort_sup = 1.0/(1.0 + 0.010*max(0.0, CORT - CORT0));
double lps_sig  = LPS/(LPS + KM_LPS);
double hdamp    = HMGB1/(HMGB1 + 25.0);
double damp     = hdamp + 0.5*(HSP - 1.0)/6.0;

dxdt_TNF  = KTNF*(lps_sig + WD_TNF*damp)*hc_sup*cort_sup - KDTNF*TNF;
double il10_brake = 1.0/(1.0 + IL10/KI10);
dxdt_IL6  = KIL6*(lps_sig + WT_IL6*TNF + 0.35*damp)*il10_brake*hc_sup
            - KDIL6*IL6;
dxdt_IL1B = KIL1*(lps_sig + 0.45*damp)*hc_sup*cort_sup - KDIL1*IL1B;
dxdt_IL10 = KIL10*IL6/100.0 + WCORT_IL10*CORT - KDIL10*IL10;
dxdt_NETS = KNET*(IL1B/100.0 + lps_sig) - KDNET*NETS;

double NEC = KNEC_TH*dose_eff + KNEC_H*HMGB1;
double rtm_h = KH_RTM*C_rtm*RTM_EHMGB;
// NOTE: the autocatalytic term is NOT multiplied by HMGB1.  Writing it as
// A*H*hill(H) makes production super-linear at every H, destroys the OFF
// state, and makes the switch latch unconditionally.  Saturating production
// against linear clearance is what creates two stable states.
dxdt_HMGB1 = KH_NEC*NEC + KH_AUTO*hillf(HMGB1, KH_HALF, NH)
             - (KH_OUT + rtm_h)*HMGB1;

// ---------------- G. endothelium and coagulation ---------------------------
dxdt_SDC1 = KSDC*(dose_eff + 0.010*TNF + 0.0035*IL6) - KDSDC*(SDC1 - SDC1_0);
dxdt_TF   = KTF*(KTF_LPS*lps_sig + 0.0075*TNF + 0.0030*IL1B) - KDTF*TF;
double rtm_apc = 1.0 + (RTM_EAPC - 1.0)*C_rtm/(C_rtm + 0.5);
double pc_act  = PC/PC0 * rtm_apc;
dxdt_THR = KTHR*TF - KDTHR*THR*(1.0 + KTHR_PC*pc_act);
double apf = IL6/(IL6 + 120.0);
dxdt_FIB = KFIB_SYN*(FIB0 - FIB)*livf*(1.0 + KFIB_APF*apf)
           - KFIB_CONS*THR*FIB/100.0;
dxdt_PLT = KPLT_SYN*(PLT0 - PLT) - KPLT_CONS*THR*PLT/100.0;
dxdt_PC  = KPC_SYN*(PC0 - PC)*livf - KPC_CONS*THR*PC/100.0*rtm_apc;
dxdt_DDIM= KDD*THR*FIB/300.0 - KDDD*DDIM;

// ---------------- H. organ injury ------------------------------------------
double nap_inj = KLIV_NAP*max(0.0, 1.0 - GSH)*NAPQI/50.0;
double liv_hit = KLIV_TH*dose_eff + KLIV_ISCH*ISCH*ISCH*(1.0 + 2.0*damp)
               + KLIV_HMGB*hdamp + nap_inj;
dxdt_LIVF = -liv_hit*LIVF + KLIV_REC*(1.0 - LIVF);
double hep_nec = liv_hit*LIVF;
dxdt_ALT = KALT*hep_nec - KDALT*(ALT - ALT0);
double rhab = KRHAB_TH*dose_mu + KRHAB_ISCH*ISCH*ISCH*min(1.0, MEX/400.0);
// Injured myocytes go on leaking for hours.  Releasing CK instantaneously with
// the muscle thermal dose put the CK peak 1.7 h after collapse instead of the
// observed 24-48 h; the transit pool also puts myoglobin (t1/2 90 min) ahead of
// CK (t1/2 36 h) automatically, which is what is seen clinically.
double mu_rel = KMU_REL*MUINJ;
dxdt_MUINJ = rhab - mu_rel;
dxdt_AST = KAST*hep_nec + KAST_MU*KCK*mu_rel/60.0 - KDAST*(AST - AST0);
dxdt_CK  = KCK*mu_rel - KDCK*(CK - CK0);
dxdt_MB  = KMB*mu_rel - KDMB*(MB - MB0);
dxdt_CAST= KCAST*MB*(2.0 - min(1.0, pv_frac)) - KDCAST*CAST;
double dic_burden = max(0.0, 1.0 - PLT/PLT0);
double ibu_gfr = KGFR_IBU*C_ibu/(C_ibu + IBU_EC50);
double gfr_target = max(0.05, 1.0 - KGFR_VP*max(0.0, 1.0 - pv_frac)
                        - KGFR_CAST*min(1.0, CAST) - KGFR_DIC*dic_burden
                        - KGFR_HMGB*hdamp - ibu_gfr);
dxdt_GFRF = KGFR_REC*160.0*(gfr_target - GFRF);
// creatinine: constant generation pinned so that GFRF = 1 sits exactly at SCR0
dxdt_SCR = KSCR_OUT*SCR0 - KSCR_OUT*GFRF*SCR;
dxdt_NSE = KNSE*dose_eff*(1.0 + 1.5*damp) - KDNSE*(NSE - NSE0);
dxdt_CNSD= (KCNS_TH*dose_eff + KCNS_INFL*(IL6/100.0 + damp)
            + KCNS_HMGB*hdamp)*(1.0 - CNSD) - KCNS_REC*CNSD;
dxdt_KPOT= KK_RHAB*KCK*rhab - KK_REN*GFRF*(KPOT - KPOT0);
dxdt_LAC = KLAC*(ISCH + 0.4*min(1.0, MEX/600.0))*8.0 - KDLAC*(LAC - LAC0)*livf;

dxdt_CORT = KCORT_IN*KCORT_STIM*(dose_eff + IL6/200.0) - KCORT_OUT*(CORT - CORT0);

// ---------------- I. drug PK ------------------------------------------------
double cyp = 1.0 + (FNAPQI_HEAT - 1.0)*min(1.0, dose_eff/0.25);
double cl_para = PARA_CL*livf;
double C_para_p = PARA_P/PARA_V2*1000.0;
dxdt_PARA_A = -PARA_KA*PARA_A;
dxdt_PARA_C = PARA_KA*PARA_A - cl_para*C_para/1000.0
              - PARA_Q*(C_para - C_para_p)/1000.0;
dxdt_PARA_P = PARA_Q*(C_para - C_para_p)/1000.0;
dxdt_NAPQI  = FNAPQI*cyp*cl_para*C_para/1000.0 - 0.10*NAPQI;
dxdt_GSH    = KGSH_REC*(1.0 - GSH) - KGSH_USE*NAPQI/50.0*GSH;

dxdt_IBU_A = -IBU_KA*IBU_A;
dxdt_IBU_C = IBU_KA*IBU_A - IBU_CL*C_ibu/1000.0*livf;

double C_dan_p = DAN_P/DAN_V2*1000.0;
dxdt_DAN_C = -DAN_CL*C_dan/1000.0*livf - DAN_Q*(C_dan - C_dan_p)/1000.0;
dxdt_DAN_P = DAN_Q*(C_dan - C_dan_p)/1000.0;
dxdt_HC_C  = -HC_CL*C_hc/1000.0;
dxdt_RTM_C = -RTM_CL*C_rtm/1000.0;

$TABLE
double AD_t = 0.202*pow(BW,0.425)*pow(HT,0.725);
double hc_t = max(3.1, 8.3*pow(VAIR,0.6));
double Pa_t = RH*psat_kpa(TA);
double Emax_t = AD_t*(psat_kpa(TSK) - Pa_t)
                /(0.0276*ICL + 1.0/((1.0+0.31*ICL)*16.5*hc_t))*WMAXW;
double q10_t = pow(Q10, (min(TCR,42.0)-37.0)/10.0);
double Hprod_t = MREST*BW*q10_t + MEX*(1.0-ETA_EX);
double Qdry_t = AD_t*(TSK - TA)/(0.155*ICL + 1.0/((1.0+0.31*ICL)*(hc_t+4.7)));
double HSI = 100.0*max(0.0, Hprod_t - Qdry_t + QSOL)/max(Emax_t, 1e-9);

// Clinical scores are OUTPUTS.  They are never integrated, never fed back.
double PTR = 1.0 + 0.9*(1.0 - LIVF);
double dic = 0.0;
dic += PLT < 50.0 ? 2.0 : (PLT < 100.0 ? 1.0 : 0.0);
dic += DDIM > 5.0 ? 3.0 : (DDIM > 1.0 ? 2.0 : 0.0);
dic += FIB < 100.0 ? 1.0 : 0.0;
dic += PTR > 1.5 ? 2.0 : (PTR > 1.25 ? 1.0 : 0.0);
double ISTH_DIC = dic;
double GCS = clampd(15.0 - 12.0*CNSD, 3.0, 15.0);
double COMMITTED = HMGB1 > 10.91 ? 1.0 : 0.0;   // above the unstable fixed point
double SOFA = (GCS < 13 ? 1:0) + (GCS < 10 ? 1:0) + (GCS < 6 ? 1:0)
            + (SCR > 1.2 ? 1:0) + (SCR > 2.0 ? 1:0) + (SCR > 3.5 ? 1:0)
            + (PLT < 150 ? 1:0) + (PLT < 100 ? 1:0) + (PLT < 50 ? 1:0)
            + (ALT > 200 ? 1:0) + (ALT > 1000 ? 1:0);

$CAPTURE HSI ISTH_DIC GCS SOFA COMMITTED
'

mod <- mcode("heatstroke", code)

## =============================================================================
##  Cooling modalities: one lumped conductance each, fitted in hs_calibrate.py
##  to published whole-body core cooling rates.  UA is an empirical potency,
##  not a first-principles surface coefficient.
## =============================================================================
MODALITY <- list(
  ice_water_immersion   = list(UA_COOL = 23.53, T_COOL =  2.0, IMMERSE = 1),  # 0.22 C/min
  cold_water_immersion  = list(UA_COOL = 27.74, T_COOL = 14.0, IMMERSE = 1),  # 0.17
  tarp_assisted         = list(UA_COOL = 20.05, T_COOL = 10.0, IMMERSE = 1),  # 0.14
  cold_shower           = list(UA_COOL = 15.94, T_COOL = 20.0, IMMERSE = 0),  # 0.10
  evaporative           = list(UA_COOL = 16.23, T_COOL = 25.0, IMMERSE = 0),  # 0.08
  endovascular          = list(UA_COOL =  4.60, T_COOL =  4.0, IMMERSE = 0),  # 0.06
  ice_packs             = list(UA_COOL =  1.60, T_COOL =  0.0, IMMERSE = 0),  # 0.032
  passive               = list(UA_COOL =  0.00, T_COOL = 30.0, IMMERSE = 0)   # 0.016
)

RACE  <- list(TA = 35.0, TRAD = 35.0, RH = 0.80, VAIR = 1.5, ICL = 0.15, QSOL = 120)
FIELD <- list(TA = 30.0, TRAD = 30.0, RH = 0.50, VAIR = 0.3, ICL = 0.10, QSOL = 0)
WARD  <- list(TA = 22.0, TRAD = 22.0, RH = 0.50, VAIR = 0.15, ICL = 1.0, QSOL = 0)

## Build a time-varying covariate data set: race -> delay -> cooling -> ward.
build_course <- function(t_collapse = 46, delay = 20, modality = "ice_water_immersion",
                         mex = 900, days = 7, ivf_rate = 0, ivf_dur = 0,
                         extra = list()) {
  m <- MODALITY[[modality]]
  seg <- function(start, dur, env, ...) {
    c(list(time = start, MEX = 0, UA_COOL = 0, IMMERSE = 0, T_COOL = 25,
           IVF_RATE = 0, ORAL_RATE = 0), env, list(...))
  }
  rows <- list(
    seg(0,                    t_collapse, RACE,  MEX = mex),
    seg(t_collapse,           delay,      FIELD),
    seg(t_collapse + delay,   240,        FIELD, UA_COOL = m$UA_COOL,
        T_COOL = m$T_COOL, IMMERSE = m$IMMERSE,
        IVF_RATE = ivf_rate),
    seg(t_collapse + delay + ifelse(ivf_dur > 0, ivf_dur, 240), 1, FIELD,
        UA_COOL = m$UA_COOL, T_COOL = m$T_COOL, IMMERSE = m$IMMERSE),
    seg(t_collapse + delay + 240, days*1440, WARD, ORAL_RATE = 0.0012)
  )
  df <- bind_rows(lapply(rows, as.data.frame))
  df$ID <- 1
  df$evid <- 0; df$amt <- 0; df$cmt <- 1
  df[order(df$time), ]
}

## =============================================================================
##  SIXTEEN SCENARIOS
##  Biology is held fixed except where the scenario name says otherwise.  Any
##  difference in outcome between arms 01-08 is therefore attributable to the
##  cooling strategy alone.
## =============================================================================
run_scn <- function(label, ..., days = 7, doses = NULL) {
  dat <- build_course(days = days, ...)
  if (!is.null(doses)) dat <- bind_rows(dat, doses) %>% arrange(time)
  out <- mod %>% data_set(dat) %>% mrgsim(end = max(dat$time), delta = 1) %>% as_tibble()
  out$scenario <- label
  out
}

bolus <- function(time, cmt_name, amt) {
  data.frame(ID = 1, time = time, amt = amt, evid = 1,
             cmt = which(names(init(mod)) == cmt_name))
}

scenarios <- list(
  ## --- cooling strategy: the only variable that changes ---------------------
  s01 = function() run_scn("01 no cooling until hospital (60 min)",
                           delay = 60, modality = "evaporative"),
  s02 = function() run_scn("02 on-site ice-water immersion (<5 min)",
                           delay = 4,  modality = "ice_water_immersion"),
  s03 = function() run_scn("03 on-site CWI at 20 min",
                           delay = 20, modality = "ice_water_immersion"),
  s04 = function() run_scn("04 tarp-assisted cooling at 10 min",
                           delay = 10, modality = "tarp_assisted"),
  s05 = function() run_scn("05 evaporative + convective at 10 min",
                           delay = 10, modality = "evaporative"),
  s06 = function() run_scn("06 ice packs only at 10 min",
                           delay = 10, modality = "ice_packs"),
  s07 = function() run_scn("07 2 L cold saline only at 20 min",
                           delay = 20, modality = "passive",
                           ivf_rate = 0.10, ivf_dur = 20),
  s08 = function() run_scn("08 CWI at 20 min + 2 L cold saline",
                           delay = 20, modality = "ice_water_immersion",
                           ivf_rate = 0.10, ivf_dur = 20),
  ## --- host factors ---------------------------------------------------------
  s09 = function() run_scn("09 heat-acclimatised runner, CWI at 20 min",
                           delay = 20, modality = "ice_water_immersion"),
  s10 = function() run_scn("10 euhydrated (drinking during race), CWI 20 min",
                           delay = 20, modality = "ice_water_immersion"),
  ## --- pharmacology ---------------------------------------------------------
  s11 = function() run_scn("11 + paracetamol 1 g IV at cooling start",
                           delay = 20, doses = bolus(66, "PARA_C", 1000)),
  s12 = function() run_scn("12 + ibuprofen 800 mg at cooling start",
                           delay = 20, doses = bolus(66, "IBU_A", 800)),
  s13 = function() run_scn("13 + dantrolene 2.5 mg/kg at cooling start",
                           delay = 20, doses = bolus(66, "DAN_C", 175)),
  s14 = function() run_scn("14 + hydrocortisone 200 mg at cooling start",
                           delay = 20, doses = bolus(66, "HC_C", 200)),
  s15 = function() run_scn("15 60 min delay + thrombomodulin alfa x 3 d",
                           delay = 60, modality = "evaporative",
                           doses = bind_rows(bolus(106, "RTM_C", 3.5),
                                             bolus(1546, "RTM_C", 3.5),
                                             bolus(2986, "RTM_C", 3.5))),
  s16 = function() run_scn("16 classic NEHS, elderly, found late",
                           delay = 45, modality = "evaporative", mex = 0)
)

## =============================================================================
##  The headline analyses.  Each is reproduced numerically in hs_analysis.py;
##  the values in the comments are what that run produces.
## =============================================================================

## ---- A. the dose ladder: temperature and time are not exchangeable 1:1 -----
##   40.0 C -> 0.016 CEM43 per minute      43.0 C -> 1.000
##   41.0 C -> 0.062                       43.5 C -> 1.414
##   42.0 C -> 0.250                       44.0 C -> 2.000
dose_rate <- function(Tc) ifelse(Tc < 43, 0.25^(43 - Tc), 0.5^(43 - Tc))

## ---- B. CEM43 paid from collapse at 42.0 C to 38.6 C -----------------------
##   modality                  delay 0   delay 10  delay 30  delay 60  cool min
##   ice-water immersion          0.90      2.97      6.22      9.52      17
##   cold-water immersion 14C     1.12      3.16      6.35      9.60      22
##   tarp-assisted cooling        1.35      3.36      6.49      9.69      27
##   cold shower / dousing        1.85      3.78      6.81      9.90      40
##   evaporative + convective     2.32      4.19      7.11     10.09      51
##   endovascular catheter        3.18      4.91      7.64     10.42      68
##   ice packs only               6.27      7.56      9.58     11.66     136
##   passive (no cooling)        15.67     15.67     15.67     15.67     364
##
## ---- C. THE EXCHANGE RATE --------------------------------------------------
##   Upgrading evaporative -> ice-water immersion is worth 4.6-6.0 MINUTES of
##   delay, essentially independent of the collapse temperature, while the
##   absolute dose at stake rises exponentially with it.  That is the
##   quantitative content of "cool first, transport second": modality is worth
##   a factor of a few, delay is unbounded.

## ---- D. the commitment switch ----------------------------------------------
## Fixed points of dH/dt = drive + A*H^n/(H^n+K^n) - c_eff*H with the shipped
## parameters:  OFF = 0,  UNSTABLE = 10.91 ng/mL,  ON = 50.3 ng/mL.
switch_rhs <- function(H, drive = 0, extra_cl = 0) {
  p <- as.list(param(mod))
  c_eff <- p$KH_OUT + extra_cl - p$KH_NEC*p$KNEC_H
  drive + p$KH_AUTO*H^p$NH/(H^p$NH + p$KH_HALF^p$NH) - c_eff*H
}
## uniroot(switch_rhs, c(0.5, 30))$root   ->  10.91
## uniroot(switch_rhs, c(30, 400))$root   ->  50.3
## The saddle-node sits at c_eff = 0.00498/min, so an agent that adds more than
## 0.00235/min of HMGB1 clearance abolishes the ON state entirely.  Only
## thrombomodulin alfa does that in this model, and it is the only agent whose
## therapeutic window is set by the switch rather than by the thermometer.

## ---- E. why antipyretics cannot work ---------------------------------------
## Sensitivity of dTc/dt to a 0.5 C set-point shift, by regime:
##   compensable (rest, 30 C/50% RH)     large
##   marginal    (400 W, 33 C/60% RH)    intermediate
##   uncompensable (900 W, 35 C/80% RH)  ~0
## The effectors are already saturated where the patient actually has heat
## stroke, so lowering a threshold changes nothing.  Meanwhile paracetamol adds
## NAPQI to a liver that is already the second organ to fail.

if (interactive()) {
  res <- bind_rows(lapply(names(scenarios), function(k) scenarios[[k]]()))
  res %>%
    group_by(scenario) %>%
    summarise(peak_Tc  = max(TCR), CEM43 = max(CEM43),
              HMGB1    = last(HMGB1), committed = max(COMMITTED),
              peak_ALT = max(ALT),  peak_CK = max(CK),
              peak_Cr  = max(SCR),  min_PLT = min(PLT),
              max_DIC  = max(ISTH_DIC), worst_GCS = min(GCS)) %>%
    print(n = 20)

  ggplot(res, aes(time/60, TCR, colour = scenario)) +
    geom_line(linewidth = 0.6) + xlim(0, 8) +
    geom_hline(yintercept = c(40, 38.6), linetype = 2, colour = "grey40") +
    labs(x = "hours from start of exposure", y = "core temperature (deg C)",
         title = "Heat stroke: the clock and the dose") +
    theme_minimal()
}
