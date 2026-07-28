##############################################################################
## tls_mrgsolve_model.R
## Tumour Lysis Syndrome — Quantitative Systems Pharmacology model
## ===========================================================================
##
## THE STRUCTURAL CLAIM
## --------------------
## TLS is normally taught as four laboratory abnormalities (K+ up, PO4 up,
## urate up, Ca down) plus acute kidney injury.  This model poses it instead as
## a RACE between two rates that are dimensionally identical and numerically
## comparable —
##
##     release flux        J_rel = Q_i * k_lys * N_lys          [mmol/h]
##     clearance capacity  C_i(GFR, urine pH, urine flow)       [mmol/h]
##
## — and then CLOSES that race into a loop, because the solute that loses it
## precipitates, and the precipitate takes away the GFR that set the clearance
## capacity:
##
##     J_rel > C_i -> concentration up -> supersaturation up -> precipitation
##                 -> GFR down -> C_i down -> ...
##
## Every therapy is then classified by WHICH TERM of that loop it touches, and
## the model's main claim is that the classes are not interchangeable:
##
##   POOL            removes solute that already exists   rasburicase, dialysis
##   FLUX            blocks new production                allopurinol, febuxostat
##   FLUX-SHAPING    lowers the release RATE              venetoclax ramp, prephase
##   DILUTION        lowers tubular concentration         hydration, furosemide
##   SPECIATION      moves the solubility curve           bicarbonate -> urine pH
##   REDISTRIBUTION  moves solute between compartments    insulin/glucose, beta2
##   SEQUESTRATION   removes solute outside the body      sevelamer, SZC
##
## 48 ODEs: 19 drug PK compartments, 3 tumour/release, 11 solute pools,
## 7 kidney/crystal, 7 safety-endpoint, 1 bookkeeping.
##
## VERIFICATION
## ------------
## The build environment for this repository has no R toolchain, so this system
## has ALSO been transcribed into `tls_reference_check.py` (numpy/scipy LSODA)
## and integrated there.  Every number quoted in README.md is produced by a
## function in that file.  The two transcriptions are meant to be the same
## system with the same parameter block; if they disagree, one of them is
## wrong.  Numbers in the comments below are from that reference integration.
##
## Usage:
##   source("tls_mrgsolve_model.R")
##   TLS_baseline()                 # calibration check
##   TLS_run_scenarios()            # the 12 shipped scenarios
##   TLS_race()                     # UA_req vs UA_crit
##   TLS_operator_decomposition()   # which term does each therapy touch
##   TLS_ph_optimum()               # the alkalinisation trade
##   TLS_leadtime()                 # allopurinol start time
##   TLS_rasburicase_capacity()     # zero-order pool operator
##   TLS_second_order_solutes()     # K+ and PO4 are consequences
##   TLS_ramp()                     # venetoclax flux-shaping
##   TLS_potassium_rescue()         # redistribution is not clearance
##   TLS_calcium_reflex()           # correcting hypocalcaemia feeds the crystal
##   TLS_bistability()              # is the loop a switch or a knee?
##   TLS_trial_ledger()             # model vs published numbers
##############################################################################

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
})

##############################################################################
## 1. THE MODEL
##############################################################################

tls_code <- '
$PROB
Tumour lysis syndrome: release flux vs clearance capacity, closed through the
kidney by intratubular and interstitial precipitation.

$PARAM @annotated
// ---------------- patient ----------------
BW      :  70.0  : body weight (kg)
BSA     :   1.73 : body surface area (m2)
VEC0    :  16.8  : baseline extracellular fluid volume (L)
VTBW    :  42.0  : total body water, creatinine/urea Vd (L)
GFR0    :   7.2  : baseline GFR (L/h) = 120 mL/min/1.73m2
HB0     :  13.5  : baseline haemoglobin (g/dL)

// ---------------- fluid / urine ----------------
FLUID   :   0.113 : total free-water input (L/h); 0.216 = 3 L/m2/day
QU_BASE :   0.083 : baseline urine flow (L/h) = 2 L/day
INSENS  :   0.030 : insensible loss (L/h)
KVOL    :   0.100 : extra urine per L of ECF excess (L/h/L)
FE_H2OM :   0.10  : maximal fractional water excretion, intact tubule
OBSTRFL :   0.90  : fraction of that lost at complete obstruction
QU_FLOOR:   0.004 : urine-flow floor (L/h)
FURO_QU :   0.15  : incremental diuresis at saturating furosemide (L/h)
EC50FUR :   1.0   : furosemide EC50 (mg/L)
KEL_FUR :   0.35  : furosemide elimination (1/h)
FURO_K  :   1.9   : multiplier on renal K clearance at saturating furosemide

// ---------------- tumour ----------------
TD      :  30.0  : tumour doubling time (h); 2000 for CLL
KLYS    :   0.1155: transit out of the committed-to-lysis pool (1/h), t1/2 6 h
KKILLSP :   0.0   : spontaneous lysis rate (1/h)
Q_K     :  38.5  : K+ released per 1e12 cells (mmol)
Q_PO4   :  60.0  : phosphate released per 1e12 cells (mmol)
Q_PUR   :  24.0  : purine released per 1e12 cells (mmol)
Q_MG    :   5.0  : Mg released per 1e12 cells (mmol)
Q_LDH   :  56.0  : LDH released per 1e12 cells (kU)
Q_CR    :   0.9  : creatinine-equivalent released per 1e12 cells (mmol)
LDH_PRD :   0.0390: basal LDH turnover (kU/h)
KEL_LDH :   0.0116: LDH elimination (1/h), t1/2 60 h

// ---------------- purine cascade ----------------
P_END   :   0.186 : endogenous purine input (mmol/h)
VXO     :   2.0   : xanthine oxidase Vmax (mmol/L/h of VEC)
KM_XO1  :   0.05  : hypoxanthine Km (mmol/L)
KM_XO2  :   0.05  : xanthine Km (mmol/L)
VSALV   :   0.075 : HGPRT salvage Vmax (mmol/L/h)
KM_SALV :   0.05  : salvage Km (mmol/L)
FE_HX   :   0.60  : fractional excretion, hypoxanthine
FE_XAN  :   0.55  : fractional excretion, xanthine
FE_UA   :   0.080 : fractional excretion, urate (basal)
FE_ALL  :   0.95  : fractional excretion, allantoin
FEUA_G  :   3.5   : gain of urate FE as URAT1 reabsorption saturates
FEUA_TH :   0.40  : urate above which FE starts to rise (mmol/L)
FEUA_C50:   0.55  : half-max of that rise (mmol/L)

// ---------------- xanthine oxidase inhibitors ----------------
KA_ALLO :   1.5   : allopurinol absorption (1/h)
V_ALLO  :  90.0   : allopurinol Vd (L)
CL_ALLO :  45.0   : allopurinol CL (L/h), t1/2 1.4 h
F_ALLO  :   0.90  : allopurinol bioavailability
FM_OXY  :   0.75  : fraction converted to oxypurinol
MW_ALLO : 136.1   : allopurinol MW
MW_OXY  : 152.1   : oxypurinol MW
V_OXY   :  40.0   : oxypurinol Vd (L)
CL_OXY  :   1.20  : oxypurinol CL (L/h) at GFR0 — SCALES WITH GFR
IC50_AL :   6.0   : allopurinol XO IC50 (mg/L)
IC50_OX :   1.10  : oxypurinol XO IC50 (mg/L)
KA_FEBU :   1.2   : febuxostat absorption (1/h)
V_FEBU  :  50.0   : febuxostat Vd (L)
CL_FEBU :   5.0   : febuxostat CL (L/h), t1/2 6.9 h, hepatic
F_FEBU  :   0.85  : febuxostat bioavailability
IC50_FE :   0.030 : febuxostat XO IC50 (mg/L)

// ---------------- rasburicase ----------------
V1_RASB :   8.0   : rasburicase central Vd (L)
V2_RASB :   6.0   : rasburicase peripheral Vd (L)
Q_RASB  :   0.6   : rasburicase intercompartmental CL (L/h)
CL_RASB :   0.30  : rasburicase CL (L/h), t1/2 18 h
VMAXRAS :   1.50  : urate oxidation, mmol/h per (mg/L) of rasburicase
KM_RASB :   0.030 : urate Km (mmol/L) — ZERO ORDER at TLS concentrations
KADA    :   0.004 : anti-drug antibody formation
ADA_CL  :   2.0   : CL multiplier at full antibody titre

// ---------------- venetoclax / cytotoxics ----------------
KA_VEN  :   0.5   : venetoclax absorption (1/h)
V_VEN   : 130.0   : venetoclax Vd/F (L)
CL_VEN  :   3.6   : venetoclax CL/F (L/h), t1/2 25 h
EMAX_VE :   0.060 : venetoclax max kill rate (1/h)
EC50_VE :   1.0   : venetoclax EC50 (mg/L)
KEL_CHE :   0.15  : generic cytotoxic elimination (1/h)
EMAX_CH :   0.055 : cytotoxic max kill rate (1/h)
EC50_CH :   1.0   : cytotoxic EC50 (relative)
KEL_PRD :   0.23  : prednisolone elimination (1/h)
EMAX_PR :   0.020 : prednisolone max kill rate (1/h)
EC50_PR :   0.10  : prednisolone EC50 (mg/L)

// ---------------- potassium ----------------
CLK0    :   0.729 : renal K clearance at GFR0 (L/h)
GFR_EXPK:   0.70  : GFR exponent on K clearance
K_ALDO  :   1.5   : aldosterone/hyperkalaemia gain
EK_COLON:   0.40  : colonic K secretion (mmol/h), GFR-independent
K_INTAKE:   3.561 : dietary K absorption (mmol/h)
KSH_INS :   0.090 : insulin-driven K shift (1/h)
EC50_INS:  60.0   : insulin EC50 (mU/L)
KEL_INS :   0.55  : insulin effect elimination (1/h)
KSH_SALB:   0.050 : beta2-driven K shift (1/h)
EC50_SAL:   8.0   : salbutamol EC50 (ng/mL)
KEL_SALB:   0.25  : salbutamol elimination (1/h)
KSH_BACK:   0.28  : return of parked K (1/h)
KSZC    :   0.55  : K bound per g of gut binder per h

// ---------------- phosphate ----------------
TMP0    :   0.9700: TmP/GFR before hormonal modulation (mmol/L)
TMP_PTH :   0.45  : fractional TmP suppression at max PTH
FGF_MAX :   0.40  : fractional TmP suppression at max FGF23
PO4_SET :   0.90  : phosphate above which FGF23 rises (mmol/L)
KFGF    :   1.20  : half-max of that rise (mmol/L)
F_ULTRA :   0.90  : ultrafilterable phosphate fraction
PO4_INT :   0.95  : net dietary phosphate absorption (mmol/h)
KSEVE   :   0.0022: phosphate bound per mg of gut binder per h

// ---------------- calcium / PTH ----------------
FE_CA   :   0.012 : fractional excretion, calcium
CA_NET0 :   0.2030: net ECF calcium input at baseline (mmol/h)
CA_BONE :   2.2   : extra bone efflux at max PTH (mmol/h)
TAU_PTH :   1.0   : PTH time constant (h)
CAIONSET:   1.15  : ionised Ca below which PTH rises (mmol/L)
CAION_HI:   1.22  : ionised Ca above which PTH is suppressed (mmol/L)
KCA_PTH :   0.12  : half-max of the PTH response (mmol/L)
PTH_MAX :   6.0   : maximal PTH (x normal)
PTH_SUPP:   0.85  : maximal PTH suppression
CAPO4_ST:   1.67  : Ca:PO4 stoichiometry in hydroxyapatite (10:6)
CAION_F :   0.50  : ionised fraction of total calcium
CAION_AL:   0.06  : loss of that fraction at full bicarbonate load
CA_HCU  :   3.0   : hypercalciuric gain on FE_CA

// ---------------- magnesium / creatinine / urea ----------------
FE_MG   :   0.030 : fractional excretion, magnesium
MG_INT  :   0.1836: net magnesium input (mmol/h)
CR_PROD :   0.55  : creatinine production (mmol/h)
UREA_PRD:  18.0   : urea production (mmol/h)
FE_UREA :   0.50  : fractional excretion, urea

// ---------------- crystal / precipitation ----------------
CF_MED  :   1.4   : medullary concentrating factor on tubular fluid
S_HU    :   0.655 : undissociated uric acid solubility (mmol/L)
PKA_UA  :   5.75  : uric acid pKa1
SS_TH_UA:   2.5   : urate metastable limit (normal urine sits at 1.9)
KN_UA   :   0.055 : urate nucleation rate constant
N_UA    :   2.0   : urate nucleation order
KDIS_UA :   0.10  : urate dissolution when SS < 1 (1/h)
S_XAN   :   0.40  : xanthine solubility (mmol/L)
PKA_XAN :   7.70  : xanthine pKa
SS_TH_XA:   1.4   : xanthine metastable limit
KN_XAN  :   0.075 : xanthine nucleation rate constant
N_XAN   :   2.0   : xanthine nucleation order
KDIS_XAN:   0.05  : xanthine dissolution (1/h)
KSP_CAP :   6.0   : tubular calcium-phosphate solubility product
PKA2_PO4:   6.80  : H2PO4-/HPO4-- pKa2
SS_TH_CAP:  4.0   : tubular CaP metastable limit
KN_CAP  :   0.004 : tubular CaP nucleation rate constant
N_CAP   :   2.0   : tubular CaP nucleation order
KDIS_CAP:   0.010 : tubular CaP dissolution (1/h)
KCLR_XT :   0.008 : clearance of deposited crystal (1/h), t1/2 3.6 d
KSP_SYS :   4.84  : systemic CaxPO4 threshold (mmol2/L2) = 60 mg2/dL2
KN_SYS  :   1.20  : systemic/interstitial deposition rate constant
N_SYS   :   2.0   : systemic deposition order
ALK_SYS :   0.90  : extra driving force at full bicarbonate load
FRAC_KID:   0.20  : share of systemic deposition that is nephrocalcinosis

// ---------------- tubular injury / nephron loss ----------------
JREF    :   0.25  : crystal flux giving unit injury drive (mmol/h)
URATETOX:   0.30  : weight of crystal-independent soluble urate toxicity
UATOX_TH:   0.48  : urate above which that term engages (mmol/L) = 8 mg/dL
UATOX_C50:  0.60  : half-max of that term (mmol/L)
KINJ    :   0.050 : injury on-rate (1/h)
KREP    :   0.030 : injury repair rate (1/h)
KNL     :   0.0015: nephron loss rate (1/h)
KNR     :   0.0040: nephron recovery rate (1/h)
XT_K50  :  12.0   : crystal mass giving half-maximal obstruction (mmol)
XT_HILL :   2.0   : obstruction Hill coefficient

// ---------------- urine pH ----------------
PH_BASE :   5.90  : baseline urine pH
PH_RISE :   1.75  : maximal rise with bicarbonate
KH_HCO3 : 140.0   : half-max bicarbonate load (mmol)
KEL_HCO3:   0.17  : bicarbonate load elimination (1/h)
TAU_PH  :   0.5   : urine pH time constant (h)

// ---------------- rasburicase oxidant limb ----------------
G6PD    :   1.0   : 1 = normal, 0 = deficient
KMETHB  :   0.0180: methaemoglobin formation per mmol H2O2
KREDMET :   0.09  : methaemoglobin reduction (1/h)
KHEM    :   0.060 : haemoglobin loss per mmol H2O2 (g/dL)
KRECHB  :   0.004 : haemoglobin recovery (1/h)

// ---------------- infusions (set by the scenario builder) ----------------
RHCO3   :   0.0   : bicarbonate infusion (mmol/h)
RCA     :   0.0   : calcium infusion (mmol/h)
RCHEMO  :   0.0   : cytotoxic infusion (relative units/h)
RINS    :   0.0   : insulin infusion (mU/L/h)
CL_HD   :   0.0   : dialysis clearance (L/h)

// ---------------- hazards ----------------
H_ARR0  :   2.0e-5: baseline arrhythmia hazard (1/h)
H_ARRMAX:   0.030 : maximal arrhythmia hazard (1/h)
ARR_K   :   5.5   : K threshold (mmol/L)
ARR_KS  :   1.50  : K half-max above threshold (mmol/L)
ARR_CA  :   1.05  : ionised Ca threshold (mmol/L)
ARR_CAS :   0.25  : ionised Ca half-max below threshold (mmol/L)
H_SZ0   :   5.0e-6: baseline seizure hazard (1/h)
H_SZMAX :   0.010 : maximal seizure hazard (1/h)
CAION_SZ:   0.90  : ionised Ca seizure threshold (mmol/L)
SZ_S    :   0.20  : half-max below that threshold (mmol/L)
H_RRT0  :   1.0e-5: baseline renal-replacement hazard (1/h)
H_RRTMAX:   0.020 : maximal renal-replacement hazard (1/h)
RRT_K   :   6.0   : K threshold (mmol/L)
RRT_KS  :   0.80  : K half-max above threshold
RRT_CR  :   2.0   : creatinine-ratio threshold
RRT_CRS :   1.50  : creatinine half-max above threshold
RRT_OLIG:   1.00  : oliguria threshold (L/day)
RRT_OLGS:   0.60  : oliguria half-max below threshold

$CMT @annotated
// ---- drug PK (19) ----
ALLO_G  : allopurinol, gut (mg)
ALLO_C  : allopurinol, central (mg)
OXY_C   : oxypurinol, central (mg)
FEBU_G  : febuxostat, gut (mg)
FEBU_C  : febuxostat, central (mg)
RASB_C  : rasburicase, central (mg)
RASB_P  : rasburicase, peripheral (mg)
ADA     : anti-rasburicase antibody (relative)
VEN_G   : venetoclax, gut (mg)
VEN_C   : venetoclax, central (mg)
CHEMO_C : cytotoxic exposure (relative)
PRED_C  : prednisolone (mg/L)
INS_C   : insulin effect compartment (mU/L)
GLU     : plasma glucose (mmol/L)
SALB_C  : salbutamol (ng/mL)
FURO_C  : furosemide (mg/L)
SEVE_G  : phosphate binder in gut (mg)
SZC_G   : potassium binder in gut (g)
HCO3_L  : bicarbonate load (mmol)
// ---- tumour / release (3) ----
N_VIA   : viable tumour (1e12 cells)
N_LYS   : committed-to-lysis pool (1e12 cells)
LDH     : serum LDH (kU in VEC)
// ---- solute pools (11) ----
K_EC    : extracellular potassium (mmol)
K_SH    : potassium parked intracellularly (mmol)
PO4_EC  : extracellular phosphate (mmol)
CA_EC   : extracellular total calcium (mmol)
MG_EC   : extracellular magnesium (mmol)
HYPOX   : hypoxanthine (mmol)
XAN     : xanthine (mmol)
URATE   : urate (mmol)
ALLANT  : allantoin (mmol)
CREAT   : creatinine (mmol)
UREA    : urea (mmol)
// ---- kidney (7) ----
NEPH    : viable nephron fraction
TUBINJ  : tubular injury (0-1)
XT_UA   : renal urate crystal (mmol)
XT_XAN  : renal xanthine crystal (mmol)
XT_CAP  : renal calcium-phosphate (mmol)
PHU     : urine pH
VEC     : extracellular fluid volume (L)
// ---- safety / endpoints (7) ----
H2O2    : cumulative peroxide from urate oxidase (mmol)
METHB   : methaemoglobin fraction
HB      : haemoglobin (g/dL)
PTH     : PTH (x normal)
CH_ARR  : cumulative arrhythmia hazard
CH_SZ   : cumulative seizure hazard
CH_RRT  : cumulative renal-replacement hazard
// ---- bookkeeping (1) ----
CUM_REL : cumulative tumour lysed (1e12 cells)

$GLOBAL
#define HILL(x, k, n) (((x) <= 0.0) ? 0.0 : pow((x), (n)) / (pow((k), (n)) + pow((x), (n))))
#define POS(x) ((x) > 0.0 ? (x) : 0.0)
// conventional-unit conversions
#define UA_MGDL   16.81
#define CR_MGDL   11.31
#define PO4_MGDL   3.10
#define CA_MGDL    4.01

$MAIN
// initial conditions are the untreated steady state reproduced by
// tls_reference_check.py::equilibrate()
N_VIA_0  = 0.05;
K_EC_0   = 4.10  * VEC0;
PO4_EC_0 = 1.15  * VEC0;
CA_EC_0  = 2.35  * VEC0;
MG_EC_0  = 0.85  * VEC0;
HYPOX_0  = 0.0004 * VEC0;
XAN_0    = 0.0004 * VEC0;
URATE_0  = 0.32  * VEC0;
CREAT_0  = 0.076 * VTBW;
UREA_0   = 5.0   * VTBW;
NEPH_0   = 1.0;
PHU_0    = PH_BASE;
VEC_0    = VEC0;
HB_0     = HB0;
PTH_0    = 1.0;
LDH_0    = 0.200 * VEC0;
GLU_0    = 6.1;

$ODE
double Vec = (VEC > 5.0) ? VEC : 5.0;

// ======================= concentrations =======================
double cK    = K_EC   / Vec;
double cPO4  = PO4_EC / Vec;
double cCaT  = CA_EC  / Vec;
double cMg   = MG_EC  / Vec;
double cHX   = HYPOX  / Vec;
double cXA   = XAN    / Vec;
double cUA   = URATE  / Vec;
double cAll  = ALLANT / Vec;
double cCr   = CREAT  / VTBW;
double cBun  = UREA   / VTBW;

double alk   = HILL(HCO3_L, KH_HCO3, 1.0);
double cCaI  = cCaT * (CAION_F - CAION_AL * alk);

// ======================= kidney =======================
double XT    = XT_UA + XT_XAN + XT_CAP;
double obstr = HILL(XT, XT_K50, XT_HILL);
double volf  = pow(Vec / VEC0, 1.5);
if (volf > 1.15) volf = 1.15;
if (volf < 0.40) volf = 0.40;
double gfr   = GFR0 * POS(NEPH) * (1.0 - obstr) * volf;
if (gfr < 0.02) gfr = 0.02;

double furo  = HILL(FURO_C, EC50FUR, 1.0);
double qu_t  = QU_BASE + POS(FLUID - QU_BASE - INSENS)
             + KVOL * (Vec - VEC0) + FURO_QU * furo;
if (qu_t < QU_FLOOR) qu_t = QU_FLOOR;
// tubular obstruction limits flow independently of GFR: this is what makes
// crystal nephropathy OLIGURIC rather than merely azotaemic
double fe_h2o = FE_H2OM * (1.0 - OBSTRFL * obstr);
double qu     = (qu_t < fe_h2o * gfr) ? qu_t : fe_h2o * gfr;

// ======================= xanthine oxidase inhibition =======================
double cAllo = ALLO_C / V_ALLO;
double cOxy  = OXY_C  / V_OXY;
double cFebu = FEBU_C / V_FEBU;
double rsum  = cAllo / IC50_AL + cOxy / IC50_OX + cFebu / IC50_FE;
double xoi   = rsum / (1.0 + rsum);
double xofree = 1.0 - xoi;

// ======================= purine fluxes =======================
double v_hx_xan = VXO   * xofree * HILL(cHX, KM_XO1, 1.0) * Vec;
double v_xan_ua = VXO   * xofree * HILL(cXA, KM_XO2, 1.0) * Vec;
double v_salv   = VSALV *          HILL(cHX, KM_SALV, 1.0) * Vec;
double cRasb    = RASB_C / V1_RASB;
double v_uox    = VMAXRAS * cRasb * HILL(cUA, KM_RASB, 1.0);

// ======================= renal solute excretion =======================
// urate reabsorption (URAT1/GLUT9) is saturable, so the FRACTIONAL EXCRETION
// rises with plasma urate.  This is the kidney's only defence, and its price
// is a higher tubular urate concentration.
double fe_ua = FE_UA * (1.0 + FEUA_G * HILL(POS(cUA - FEUA_TH), FEUA_C50, 1.0));
double e_ua  = fe_ua * gfr * cUA;
double e_hx  = FE_HX  * gfr * cHX;
double e_xan = FE_XAN * gfr * cXA;
double e_all = FE_ALL * gfr * cAll;

double aldo = 1.0 + K_ALDO * POS((cK - 4.0) / 4.0);
double clk  = CLK0 * pow(gfr / GFR0, GFR_EXPK) * aldo
            * (1.0 + (FURO_K - 1.0) * furo);
double e_k  = clk * cK + EK_COLON * (1.0 + 1.5 * POS(cK - 4.0));

double pth_f = (PTH - 1.0) / (PTH_MAX - 1.0);
if (pth_f > 1.0) pth_f = 1.0;
double fgf   = FGF_MAX * HILL(POS(cPO4 - PO4_SET), KFGF, 1.0);
double tmp   = TMP0 * (1.0 - TMP_PTH * POS(pth_f)) * (1.0 - fgf);
double filt  = gfr * cPO4 * F_ULTRA;
double reabs = (filt < tmp * gfr) ? filt : tmp * gfr;
double e_po4 = POS(filt - reabs);

double e_ca  = FE_CA * (1.0 + CA_HCU * HILL(POS(cCaT - 2.55), 0.40, 1.0))
             * gfr * cCaT;
double e_mg  = FE_MG * gfr * cMg;
double e_cr  = gfr * cCr;
double e_urea= FE_UREA * gfr * cBun;

// ======================= tubular fluid concentrations =======================
double cu_ua  = CF_MED * e_ua  / qu;
double cu_xan = CF_MED * e_xan / qu;
double cu_ca  = CF_MED * e_ca  / qu;
double cu_po4 = CF_MED * e_po4 / qu;

// ======================= solubility and supersaturation =======================
double s_ua   = S_HU  * (1.0 + pow(10.0, PHU - PKA_UA));
double ss_ua  = cu_ua / s_ua;
double s_xan  = S_XAN * (1.0 + pow(10.0, PHU - PKA_XAN));
double ss_xan = cu_xan / s_xan;
double f_hpo4 = 1.0 / (1.0 + pow(10.0, PKA2_PO4 - PHU));
double ss_cap = cu_ca * cu_po4 * f_hpo4 / KSP_CAP;
double ss_sys = cCaT * cPO4 * (1.0 + ALK_SYS * alk) / KSP_SYS;

// ======================= precipitation fluxes =======================
// nucleation rate, then the hard physical bound: a crystal cannot take up more
// solute per hour than the tubule DELIVERS per hour.  Without this cap,
// C_urine = E/qu diverges as urine flow goes to zero and an obstructed kidney
// keeps depositing crystal it is no longer being given.
double n_ua  = KN_UA  * qu * pow(POS(ss_ua  - SS_TH_UA ), N_UA );
double n_xan = KN_XAN * qu * pow(POS(ss_xan - SS_TH_XA ), N_XAN);
double n_cap = KN_CAP * qu * pow(POS(ss_cap - SS_TH_CAP), N_CAP);
double dep_ua  = (n_ua  < 0.95 * e_ua ) ? n_ua  : 0.95 * e_ua;
double dep_xan = (n_xan < 0.95 * e_xan) ? n_xan : 0.95 * e_xan;
double dep_cpt = (n_cap < 0.95 * e_po4) ? n_cap : 0.95 * e_po4;

double po4_av = HILL(cPO4, 0.20, 1.0);
double ca_av  = HILL(cCaI, 0.35, 1.0);
double j_sys  = KN_SYS * Vec * pow(POS(ss_sys - 1.0), N_SYS) * po4_av * ca_av;

double j_ua   = dep_ua  - KDIS_UA  * XT_UA  * POS(1.0 - ss_ua )
              - KCLR_XT * XT_UA;
double j_xan  = dep_xan - KDIS_XAN * XT_XAN * POS(1.0 - ss_xan)
              - KCLR_XT * XT_XAN;
double j_cptn = dep_cpt - KDIS_CAP * XT_CAP * POS(1.0 - ss_cap);
double j_cap  = j_cptn + FRAC_KID * j_sys - KCLR_XT * XT_CAP;
double j_po4_dep = POS(j_cptn) + POS(j_sys);

double dep_cap_g = POS(j_cptn) + FRAC_KID * j_sys;
double j_xtal = dep_ua + dep_xan + dep_cap_g;
double inj_urate = URATETOX * HILL(POS(cUA - UATOX_TH), UATOX_C50, 1.0);
double inj_drive = j_xtal / JREF + inj_urate;

// ======================= drug PK =======================
dxdt_ALLO_G = -KA_ALLO * ALLO_G;
dxdt_ALLO_C =  KA_ALLO * ALLO_G * F_ALLO - CL_ALLO * cAllo;
dxdt_OXY_C  =  FM_OXY * CL_ALLO * cAllo * MW_OXY / MW_ALLO
             - CL_OXY * (gfr / GFR0) * cOxy;
dxdt_FEBU_G = -KA_FEBU * FEBU_G;
dxdt_FEBU_C =  KA_FEBU * FEBU_G * F_FEBU - CL_FEBU * cFebu;

double ada_c  = (ADA < 1.0) ? ADA : 1.0;
double cl_ras = CL_RASB * (1.0 + (ADA_CL - 1.0) * ada_c);
double cRasbP = RASB_P / V2_RASB;
dxdt_RASB_C = -cl_ras * cRasb - Q_RASB * (cRasb - cRasbP);
dxdt_RASB_P =  Q_RASB * (cRasb - cRasbP);
dxdt_ADA    =  KADA * cRasb - 0.002 * ADA;

double cVen = VEN_C / V_VEN;
dxdt_VEN_G  = -KA_VEN * VEN_G;
dxdt_VEN_C  =  KA_VEN * VEN_G - CL_VEN * cVen;

dxdt_CHEMO_C = RCHEMO - KEL_CHE * CHEMO_C;
dxdt_PRED_C  = -KEL_PRD * PRED_C;
dxdt_INS_C   = RINS - KEL_INS * INS_C;
dxdt_GLU     = 0.55 - 0.09 * GLU * (1.0 + 0.030 * INS_C);
dxdt_SALB_C  = -KEL_SALB * SALB_C;
dxdt_FURO_C  = -KEL_FUR * FURO_C;
dxdt_SEVE_G  = -0.35 * SEVE_G;
dxdt_SZC_G   = -0.35 * SZC_G;
dxdt_HCO3_L  = RHCO3 - KEL_HCO3 * HCO3_L;

// ======================= tumour kill and release =======================
double kkill = KKILLSP
             + EMAX_VE * HILL(cVen,    EC50_VE, 1.0)
             + EMAX_CH * HILL(CHEMO_C, EC50_CH, 1.0)
             + EMAX_PR * HILL(PRED_C,  EC50_PR, 1.0);
dxdt_N_VIA = (log(2.0) / TD) * N_VIA - kkill * N_VIA;
dxdt_N_LYS = kkill * N_VIA - KLYS * N_LYS;
double rel = KLYS * N_LYS;
dxdt_LDH   = Q_LDH * rel + LDH_PRD - KEL_LDH * LDH;
dxdt_CUM_REL = rel;

// ======================= potassium =======================
double shift_in = (KSH_INS  * HILL(INS_C,  EC50_INS, 1.0)
                 + KSH_SALB * HILL(SALB_C, EC50_SAL, 1.0)) * K_EC;
double j_shift  = shift_in - KSH_BACK * K_SH;
double kfrac    = (cK / 4.0 < 1.0) ? cK / 4.0 : 1.0;
double j_szc    = KSZC * SZC_G * kfrac;
dxdt_K_EC = Q_K * rel + K_INTAKE - e_k - j_shift - j_szc - CL_HD * cK;
dxdt_K_SH = j_shift;

// ======================= phosphate =======================
dxdt_PO4_EC = Q_PO4 * rel + PO4_INT - e_po4 - KSEVE * SEVE_G * po4_av
            - j_po4_dep - CL_HD * cPO4 * 0.9;

// ======================= calcium / PTH =======================
dxdt_CA_EC = CA_NET0 + CA_BONE * pth_f + RCA - e_ca - CAPO4_ST * j_po4_dep;
// PTH responds only to ionised calcium OUTSIDE its set-point band, so the
// basal state is PTH = 1 with zero net PTH-driven bone efflux
double pth_ss = 1.0
              + (PTH_MAX - 1.0) * HILL(POS(CAIONSET - cCaI), KCA_PTH, 2.0)
              - PTH_SUPP        * HILL(POS(cCaI - CAION_HI), KCA_PTH, 2.0);
dxdt_PTH = (pth_ss - PTH) / TAU_PTH;

dxdt_MG_EC = Q_MG * rel + MG_INT - e_mg - CL_HD * cMg;

// ======================= purines =======================
dxdt_HYPOX  = Q_PUR * rel + P_END - v_hx_xan - v_salv - e_hx - CL_HD * cHX;
dxdt_XAN    = v_hx_xan - v_xan_ua - e_xan - POS(j_xan) - CL_HD * cXA;
dxdt_URATE  = v_xan_ua - e_ua - v_uox - POS(j_ua) - CL_HD * cUA;
dxdt_ALLANT = v_uox - e_all - CL_HD * cAll;
dxdt_CREAT  = CR_PROD + Q_CR * rel - e_cr - CL_HD * cCr;
dxdt_UREA   = UREA_PRD - e_urea - CL_HD * cBun * 1.2;

// ======================= kidney states =======================
dxdt_TUBINJ = KINJ * inj_drive * (1.0 - TUBINJ) - KREP * TUBINJ;
dxdt_NEPH   = -KNL * TUBINJ * TUBINJ * NEPH
            +  KNR * (1.0 - NEPH) * POS(1.0 - TUBINJ / 0.30);
dxdt_XT_UA  = j_ua;
dxdt_XT_XAN = j_xan;
dxdt_XT_CAP = j_cap;
double ph_ss = PH_BASE + PH_RISE * HILL(HCO3_L, KH_HCO3, 1.0);
dxdt_PHU    = (ph_ss - PHU) / TAU_PH;
dxdt_VEC    = FLUID - qu - INSENS;

// ======================= oxidant limb =======================
dxdt_H2O2 = v_uox;
double ox = v_uox * (1.0 - G6PD);
dxdt_METHB = KMETHB * ox - KREDMET * METHB;
dxdt_HB    = -KHEM * ox + KRECHB * (HB0 - HB);

// ======================= hazards =======================
double crrat = cCr * CR_MGDL / 0.86;
double quday = qu * 24.0;
dxdt_CH_ARR = H_ARR0 + H_ARRMAX * HILL(POS(cK - ARR_K), ARR_KS, 4.0)
            * (1.0 + HILL(POS(ARR_CA - cCaI), ARR_CAS, 2.0));
dxdt_CH_SZ  = H_SZ0  + H_SZMAX  * HILL(POS(CAION_SZ - cCaI), SZ_S, 3.0);
dxdt_CH_RRT = H_RRT0 + H_RRTMAX / 3.0 * (
                HILL(POS(cK - RRT_K), RRT_KS, 3.0)
              + HILL(POS(crrat - RRT_CR), RRT_CRS, 3.0)
              + HILL(POS(RRT_OLIG - quday), RRT_OLGS, 3.0));

$TABLE
double Vecx = (VEC > 5.0) ? VEC : 5.0;
double alkx = HILL(HCO3_L, KH_HCO3, 1.0);
double XTx  = XT_UA + XT_XAN + XT_CAP;
double obs  = HILL(XTx, XT_K50, XT_HILL);
double vf   = pow(Vecx / VEC0, 1.5);
if (vf > 1.15) vf = 1.15;
if (vf < 0.40) vf = 0.40;
double GFRt = GFR0 * POS(NEPH) * (1.0 - obs) * vf;
if (GFRt < 0.02) GFRt = 0.02;
double fur  = HILL(FURO_C, EC50FUR, 1.0);
double qut  = QU_BASE + POS(FLUID - QU_BASE - INSENS)
            + KVOL * (Vecx - VEC0) + FURO_QU * fur;
if (qut < QU_FLOOR) qut = QU_FLOOR;
double feh  = FE_H2OM * (1.0 - OBSTRFL * obs);
double QUt  = (qut < feh * GFRt) ? qut : feh * GFRt;

capture UA      = URATE  / Vecx;
capture UA_mgdl = URATE  / Vecx * UA_MGDL;
capture K       = K_EC   / Vecx;
capture PO4     = PO4_EC / Vecx;
capture PO4_mgdl= PO4_EC / Vecx * PO4_MGDL;
capture CaT     = CA_EC  / Vecx;
capture CaIon   = CA_EC  / Vecx * (CAION_F - CAION_AL * alkx);
capture Mg      = MG_EC  / Vecx;
capture Cr      = CREAT  / VTBW * CR_MGDL;
capture BUN     = UREA   / VTBW;
capture LDH_UL  = 1000.0 * LDH / Vecx;
capture eGFR    = GFRt / GFR0 * 120.0;
capture QU_Lday = QUt * 24.0;
capture urinepH = PHU;
capture OBSTR   = obs;
capture XTtot   = XTx;
capture CaxPO4  = (CA_EC / Vecx * CA_MGDL) * (PO4_EC / Vecx * PO4_MGDL);
capture SS_UA   = CF_MED * FE_UA
                * (1.0 + FEUA_G * HILL(POS(URATE / Vecx - FEUA_TH), FEUA_C50, 1.0))
                * GFRt * (URATE / Vecx) / QUt
                / (S_HU * (1.0 + pow(10.0, PHU - PKA_UA)));
capture SS_SYS  = (CA_EC / Vecx) * (PO4_EC / Vecx) * (1.0 + ALK_SYS * alkx) / KSP_SYS;
capture XOI     = (ALLO_C / V_ALLO / IC50_AL + OXY_C / V_OXY / IC50_OX
                 + FEBU_C / V_FEBU / IC50_FE)
                / (1.0 + ALLO_C / V_ALLO / IC50_AL + OXY_C / V_OXY / IC50_OX
                 + FEBU_C / V_FEBU / IC50_FE);
capture C_OXY   = OXY_C  / V_OXY;
capture C_FEBU  = FEBU_C / V_FEBU;
capture C_RASB  = RASB_C / V1_RASB;
capture C_VEN   = VEN_C  / V_VEN;
capture MetHb   = METHB * 100.0;
capture P_ARR   = 1.0 - exp(-CH_ARR);
capture P_SZ    = 1.0 - exp(-CH_SZ);
capture P_RRT   = 1.0 - exp(-CH_RRT);
'

mod <- mrgsolve::mcode("tls", tls_code, soloc = tempdir())

##############################################################################
## 2. HELPERS AND THE THREE CLOSED-FORM QUANTITIES
##############################################################################

pget <- function(name) as.numeric(param(mod)[[name]])

HYD_STD  <- 0.113            # 2 L/day urine
HYD_AGGR <- 3.0 * 1.73 / 24  # 3 L/m2/day = 0.216 L/h
HYD_MAX  <- 4.0 * 1.73 / 24  # 4 L/m2/day = 0.288 L/h

hillf <- function(x, k, n = 1) ifelse(x <= 0, 0, x^n / (k^n + x^n))

## fractional excretion of urate as a function of plasma urate
fe_ua_of <- function(ua) {
  pget("FE_UA") * (1 + pget("FEUA_G") *
                     hillf(pmax(0, ua - pget("FEUA_TH")), pget("FEUA_C50")))
}

## UA_req: the plasma urate at which renal excretion balances a release flux
ua_required <- function(jrel, gfr = pget("GFR0")) {
  f <- function(u) fe_ua_of(u) * gfr * u - jrel
  uniroot(f, c(1e-9, 200), tol = 1e-10)$root
}

## UA_crit: the plasma urate at which tubular fluid reaches the metastable limit
ua_critical <- function(qu, ph = pget("PH_BASE"), gfr = pget("GFR0")) {
  s_ua <- pget("S_HU") * (1 + 10^(ph - pget("PKA_UA")))
  target <- pget("SS_TH_UA") * s_ua
  f <- function(u) pget("CF_MED") * fe_ua_of(u) * gfr * u / qu - target
  uniroot(f, c(1e-9, 200), tol = 1e-10)$root
}

## Peak release flux in closed form.  With N' = -(k-g)N and L' = kN - kl*L,
##   L(t) = k*N0/(kl-a) * (exp(-a t) - exp(-kl t)),  a = k - g
## which peaks at t* = log(kl/a)/(kl-a); flux = q * kl * L(t*).
peak_release_flux <- function(N0, q = pget("Q_PUR"), kkill = pget("EMAX_CH"),
                              td = pget("TD")) {
  a  <- kkill - log(2) / td
  kl <- pget("KLYS")
  if (a <= 1e-9) return(q * kkill * N0)
  tstar <- log(kl / a) / (kl - a)
  L <- kkill * N0 / (kl - a) * (exp(-a * tstar) - exp(-kl * tstar))
  q * kl * L
}

## Potassium excretion capacity at a given serum K
k_capacity <- function(kconc, gfr = pget("GFR0")) {
  aldo <- 1 + pget("K_ALDO") * pmax(0, (kconc - 4) / 4)
  clk <- pget("CLK0") * (gfr / pget("GFR0"))^pget("GFR_EXPK") * aldo
  clk * kconc + pget("EK_COLON") * (1 + 1.5 * pmax(0, kconc - 4))
}

## Rasburicase capacity (zero order in urate at TLS concentrations)
rasburicase_capacity <- function(dose_mgkg) {
  pget("VMAXRAS") * dose_mgkg * pget("BW") / pget("V1_RASB")
}

##############################################################################
## 3. SCENARIO BUILDER
##############################################################################

## Build an event object.  Negative times are lead-time prophylaxis; the tumour
## is wound back along its own growth curve so that N0 is the burden AT t = 0
## in the absence of pre-treatment.  Without that, a 5-day steroid prephase arm
## would silently start from a tumour 2^(120/30) = 16-fold larger and the
## comparison would be rigged in its favour.
tls_events <- function(allo_start = NA, allo_dose = 300, febu_start = NA,
                       rasb_start = NA, rasb_mgkg = 0.20, rasb_days = 1,
                       szc_start = NA, seve_start = NA, furo_start = NA,
                       pred_start = NA, pred_days = 5, ven = NULL,
                       tend = 336) {
  ev <- NULL
  add <- function(e) if (is.null(ev)) e else c(ev, e)
  if (!is.na(allo_start)) {
    n <- ceiling((tend - allo_start) / 24)
    ev <- add(ev(time = allo_start, amt = allo_dose, cmt = "ALLO_G",
                 ii = 24, addl = n - 1))
  }
  if (!is.na(febu_start)) {
    n <- ceiling((tend - febu_start) / 24)
    ev <- add(ev(time = febu_start, amt = 120, cmt = "FEBU_G",
                 ii = 24, addl = n - 1))
  }
  if (!is.na(rasb_start)) {
    ev <- add(ev(time = rasb_start, amt = rasb_mgkg * pget("BW"),
                 cmt = "RASB_C", ii = 24, addl = rasb_days - 1))
  }
  if (!is.na(szc_start)) {
    n <- ceiling((tend - szc_start) / 8)
    ev <- add(ev(time = szc_start, amt = 10, cmt = "SZC_G", ii = 8, addl = n - 1))
  }
  if (!is.na(seve_start)) {
    n <- ceiling((tend - seve_start) / 8)
    ev <- add(ev(time = seve_start, amt = 1600, cmt = "SEVE_G",
                 ii = 8, addl = n - 1))
  }
  if (!is.na(furo_start)) {
    n <- ceiling((tend - furo_start) / 8)
    ev <- add(ev(time = furo_start, amt = 1.4, cmt = "FURO_C",
                 ii = 8, addl = n - 1))
  }
  if (!is.na(pred_start)) {
    ev <- add(ev(time = pred_start, amt = 1.2, cmt = "PRED_C",
                 ii = 24, addl = pred_days - 1))
  }
  if (!is.null(ven)) {
    for (i in seq_len(nrow(ven))) {
      t1 <- if (i < nrow(ven)) ven$time[i + 1] else tend
      n <- max(1, ceiling((t1 - ven$time[i]) / 24))
      ev <- add(ev(time = ven$time[i], amt = ven$dose[i], cmt = "VEN_G",
                   ii = 24, addl = n - 1))
    }
  }
  ev
}

tls_sim <- function(N0 = 3.0, hydration = HYD_AGGR, tend = 336, delta = 0.25,
                    chemo = c(0, 24, 1.0), hco3 = 0, ca_inf = 0,
                    dialysis = NULL, insulin = NULL, neph0 = 1.0,
                    xtua0 = 0, extra_param = list(), ...) {
  events <- tls_events(tend = tend, ...)
  tstart <- if (is.null(events)) 0 else min(0, min(events$time))
  ## wind the tumour back so that N0 is the burden at t = 0
  td <- if (!is.null(extra_param$TD)) extra_param$TD else pget("TD")
  n_start <- N0 * 2^(tstart / td)

  m <- mod %>%
    param(FLUID = hydration, RHCO3 = hco3, RCA = ca_inf) %>%
    param(extra_param) %>%
    init(N_VIA = n_start, NEPH = neph0, XT_UA = xtua0)

  ## piecewise-constant infusions implemented as a covariate data set
  tg <- seq(tstart, tend, by = delta)
  idata <- data.frame(ID = 1)
  dat <- data.frame(
    ID = 1, time = tg,
    RCHEMO = ifelse(!is.null(chemo) & tg >= chemo[1] & tg < chemo[2],
                    chemo[3], 0),
    RINS   = if (is.null(insulin)) 0 else
      ifelse(tg >= insulin[1] & tg < insulin[2], insulin[3], 0),
    CL_HD  = if (is.null(dialysis)) 0 else
      ifelse(tg >= dialysis[1] & tg < dialysis[2], dialysis[3], 0))

  if (is.null(events)) {
    out <- m %>% data_set(dat) %>% carry_out(RCHEMO, CL_HD) %>% mrgsim(obsonly = FALSE)
  } else {
    out <- m %>% data_set(dat) %>% ev(events) %>%
      carry_out(RCHEMO, CL_HD) %>% mrgsim(obsonly = FALSE)
  }
  as.data.frame(out)
}

##############################################################################
## 4. SUMMARY AND THE CAIRO-BISHOP CLASSIFICATION
##############################################################################

tls_summary <- function(d, label = "") {
  post <- d[d$time >= 0, ]
  cr0  <- approx(d$time, d$Cr, 0)$y
  data.frame(
    label     = label,
    UA_peak   = max(post$UA_mgdl),
    K_peak    = max(post$K),
    PO4_peak  = max(post$PO4),
    Ca_nadir  = min(post$CaIon),
    Cr_peak   = max(post$Cr),
    Cr_ratio  = max(post$Cr) / cr0,
    eGFR_nadir= min(post$eGFR),
    XT_UA     = tail(d$XT_UA, 1),
    XT_XAN    = tail(d$XT_XAN, 1),
    XT_CAP    = tail(d$XT_CAP, 1),
    XT_max    = max(post$XTtot),
    NEPH_end  = tail(d$NEPH, 1),
    LDH_peak  = max(post$LDH_UL),
    lysed     = tail(d$CUM_REL, 1),
    lysed_pre = approx(d$time, d$CUM_REL, 0)$y,
    H2O2      = tail(d$H2O2, 1),
    MetHb     = max(d$MetHb),
    Hb_nadir  = min(d$HB),
    P_arr     = tail(d$P_ARR, 1),
    P_sz      = tail(d$P_SZ, 1),
    P_rrt     = tail(d$P_RRT, 1),
    urine_pH  = max(d$urinepH),
    LTLS      = cairo_bishop(d)$LTLS,
    CTLS      = cairo_bishop(d)$CTLS,
    n_crit    = cairo_bishop(d)$n,
    stringsAsFactors = FALSE)
}

## Cairo & Bishop 2004.  Laboratory TLS: >= 2 of the four, from 3 days before to
## 7 days after therapy.  Clinical TLS: LTLS + creatinine >= 1.5 x ULN, or
## arrhythmia, or seizure.  `sample_h` evaluates the criteria only on a discrete
## blood-draw schedule, which is what actually happens on a ward.
cairo_bishop <- function(d, sample_h = NULL, window = c(-72, 168)) {
  w <- d[d$time >= window[1] & d$time <= window[2], ]
  if (!is.null(sample_h)) {
    want <- seq(max(min(w$time), window[1]), min(max(w$time), window[2]),
                by = sample_h)
    w <- w[sapply(want, function(x) which.min(abs(w$time - x))), ]
  }
  b <- d[which.min(abs(d$time - 0)), ]
  c_ua <- any(w$UA_mgdl >= 8.0 | w$UA_mgdl >= 1.25 * b$UA_mgdl)
  c_k  <- any(w$K       >= 6.0 | w$K       >= 1.25 * b$K)
  c_p  <- any(w$PO4     >= 1.45| w$PO4     >= 1.25 * b$PO4)
  c_c  <- any(w$CaT     <= 1.75| w$CaT     <= 0.75 * b$CaT)
  n <- sum(c_ua, c_k, c_p, c_c)
  ltls <- n >= 2
  crit <- any(w$Cr >= 1.5 * 0.86) ||
    tail(d$P_ARR, 1) > 0.05 || tail(d$P_SZ, 1) > 0.05
  list(LTLS = ltls, CTLS = ltls && crit, n = n,
       urate = c_ua, K = c_k, PO4 = c_p, Ca = c_c)
}

##############################################################################
## 5. ANALYSIS FUNCTIONS
## Reference values in the comments come from tls_reference_check.py.
##############################################################################

## ---------------------------------------------------------------------------
## A0.  Baseline calibration.  Expected: urate 5.17 mg/dL, K 4.10, PO4 1.15,
## Ca 2.35 (ionised 1.17), Cr 0.86 mg/dL, eGFR 120, urine 1.99 L/day,
## urate excretion 715 mg/day, K 85 mmol/day, PO4 23 mmol/day.
## ---------------------------------------------------------------------------
TLS_baseline <- function() {
  d <- tls_sim(N0 = 1e-9, chemo = NULL, tend = 600, delta = 4)
  b <- tail(d, 1)
  cat("\n=== A0. BASELINE STEADY STATE (no tumour, no drug) ===\n")
  ref <- c(UA_mgdl = "4-6", K = "3.5-5.0", PO4 = "0.8-1.45", CaT = "2.2-2.6",
           CaIon = "1.12-1.30", Cr = "0.7-1.2", BUN = "3-7", eGFR = "~120",
           QU_Lday = "1.5-2.5", urinepH = "5.5-6.5", LDH_UL = "140-280",
           SS_UA = "1.5-2.5 (metastable)", SS_SYS = "<1")
  for (k in names(ref))
    cat(sprintf("  %-10s %9.2f    target %s\n", k, b[[k]], ref[[k]]))
  cat(sprintf("\n  tubular injury %.3f, nephron mass %.3f (both must be at rest)\n",
              b$TUBINJ, b$NEPH))
  invisible(b)
}

## ---------------------------------------------------------------------------
## A1.  THE RACE.  UA_req is the plasma urate needed to EXCRETE the release
## flux; UA_crit is the plasma urate at which tubular fluid reaches the
## metastable limit.  The disease is the region where UA_req > UA_crit.
##
## Reference output: UA_crit = 6.8 / 9.0 / 10.5 mg/dL at 2 / 3 / 4 L urine,
## and 33.2 mg/dL at pH 7.0.  The Cairo-Bishop urate cut-off is 8.0 mg/dL,
## which the model therefore reproduces as a crystallisation threshold rather
## than as a chosen number.  req/crit crosses 1 near 0.6e12 cells.
##
## Solute hierarchy at intact GFR, burden 3e12: urate flux/capacity 6.2,
## potassium 0.4, phosphate 0.5.  ONLY URATE LOSES THE RACE AT INTACT GFR.
## ---------------------------------------------------------------------------
TLS_race <- function(burdens = c(0.05, 0.3, 0.5, 1, 2, 3, 5, 8)) {
  cat("\n=== A1. THE RACE: UA_req vs UA_crit ===\n")
  cat("\n  UA_crit is a property of the PRESCRIPTION, not of the tumour:\n")
  for (r in list(c(HYD_STD, 5.90, "2 L/day, pH 5.9"),
                 c(HYD_AGGR, 5.90, "3 L/m2/day, pH 5.9"),
                 c(HYD_MAX, 5.90, "4 L/m2/day, pH 5.9"),
                 c(HYD_AGGR, 7.00, "3 L/m2/day, pH 7.0"))) {
    qu <- as.numeric(r[[1]]) - pget("INSENS")
    uc <- ua_critical(qu, ph = as.numeric(r[[2]]))
    cat(sprintf("    %-22s UA_crit = %5.2f mmol/L = %5.1f mg/dL\n",
                r[[3]], uc, uc * 16.81))
  }
  qu <- HYD_AGGR - pget("INSENS")
  uc <- ua_critical(qu)
  res <- do.call(rbind, lapply(burdens, function(N0) {
    s <- tls_summary(tls_sim(N0 = N0))
    j <- peak_release_flux(N0)
    ur <- ua_required(j)
    data.frame(burden = N0, urate_g = N0 * pget("Q_PUR") * 0.168,
               Jrel = j, UA_req = ur * 16.81, UA_crit = uc * 16.81,
               req_crit = ur / uc, UA_peak = s$UA_peak, K_peak = s$K_peak,
               Cr_ratio = s$Cr_ratio, TLS = ifelse(s$CTLS, "C",
                                                   ifelse(s$LTLS, "L", "-")))
  }))
  print(res, row.names = FALSE, digits = 3)
  cat("\n  Same race for the other two solutes at 3e12 cells:\n")
  h <- data.frame(
    solute = c("urate", "potassium", "phosphate"),
    content = c(pget("Q_PUR"), pget("Q_K"), pget("Q_PO4")),
    Jrel = c(peak_release_flux(3, q = pget("Q_PUR")),
             peak_release_flux(3, q = pget("Q_K")),
             peak_release_flux(3, q = pget("Q_PO4"))),
    capacity = c(fe_ua_of(0.476) * pget("GFR0") * 0.476, k_capacity(6.0), 12.0))
  h$ratio <- h$Jrel / h$capacity
  print(h, row.names = FALSE, digits = 3)
  cat("  -> only urate loses the race at intact GFR; the other two lose it\n")
  cat("     only after urate's crystals have taken the GFR away.\n")
  invisible(res)
}

## ---------------------------------------------------------------------------
## A2.  OPERATOR DECOMPOSITION.  Which term of the loop does each therapy
## touch?  Reference (3e12 cells, peak Cr ratio):
##   nothing 2 L/day 3.06 | 3 L/m2 1.55 | 4 L/m2 1.27 | allopurinol 1.26
##   febuxostat 1.18 | rasburicase 1.15 | dialysis 1.52 | pH 7.5 1.18
##   phosphate binder 1.55 (i.e. nothing) | steroid prephase 1.04
## ---------------------------------------------------------------------------
TLS_operator_decomposition <- function(N0 = 3.0) {
  cat("\n=== A2. OPERATOR DECOMPOSITION ===\n")
  arms <- list(
    "nothing (2 L/day)"        = list(hydration = HYD_STD),
    "DILUTION 3 L/m2/day"      = list(hydration = HYD_AGGR),
    "DILUTION 4 L/m2/day"      = list(hydration = HYD_MAX),
    "DILUTION + furosemide"    = list(furo_start = 0),
    "FLUX allopurinol t=0"     = list(allo_start = 0),
    "FLUX febuxostat t=0"      = list(febu_start = 0),
    "POOL rasburicase"         = list(rasb_start = 0),
    "POOL dialysis d2-5"       = list(dialysis = c(48, 120, 1.2)),
    "SPECIATION pH 7.5"        = list(hco3 = 45),
    "SEQUESTRATION PO4 binder" = list(seve_start = 0),
    "FLUX-SHAPE prephase"      = list(pred_start = -120, pred_days = 5))
  res <- do.call(rbind, lapply(names(arms), function(nm) {
    a <- arms[[nm]]; a$N0 <- N0
    s <- tls_summary(do.call(tls_sim, a), nm)
    s[, c("label", "UA_peak", "K_peak", "PO4_peak", "Ca_nadir", "Cr_ratio",
          "XT_UA", "XT_CAP", "P_rrt")]
  }))
  print(res, row.names = FALSE, digits = 3)
  cat("\n  No single arm lowers BOTH crystal columns.\n")
  invisible(res)
}

## ---------------------------------------------------------------------------
## A3.  THE URINE pH TRADE.  Urate solubility rises 24-fold from pH 5.9 to 7.5;
## the HPO4(2-) fraction rises 7.5-fold over the same range and systemic
## alkalosis lowers ionised calcium as well.  The optimum therefore depends on
## which solute is rate-limiting.
##
## Reference: optimum urine pH 6.60 at 3e12 without rasburicase (Cr ratio 1.15
## vs 1.55 at pH 5.90); 5.90 (i.e. no alkali at all) WITH rasburicase; 7.00 at
## 6e12 without rasburicase.  This derives the guideline change rather than
## asserting it.
## ---------------------------------------------------------------------------
TLS_ph_optimum <- function(burden = 3.0, rasb = FALSE) {
  cat(sprintf("\n=== A3. URINE pH OPTIMUM (burden %.1fe12, rasburicase %s) ===\n",
              burden, ifelse(rasb, "yes", "no")))
  f0 <- 1 / (1 + 10^(pget("PKA2_PO4") - pget("PH_BASE")))
  s0 <- pget("S_HU") * (1 + 10^(pget("PH_BASE") - pget("PKA_UA")))
  curves <- do.call(rbind, lapply(c(5.5, 5.9, 6.5, 7.0, 7.5, 7.8), function(ph) {
    sua <- pget("S_HU") * (1 + 10^(ph - pget("PKA_UA")))
    fh <- 1 / (1 + 10^(pget("PKA2_PO4") - ph))
    data.frame(pH = ph, S_UA = sua,
               UA_crit = ua_critical(HYD_AGGR - pget("INSENS"), ph) * 16.81,
               f_HPO4 = fh, UA_gain = sua / s0, CaP_cost = fh / f0)
  }))
  print(curves, row.names = FALSE, digits = 3)
  res <- do.call(rbind, lapply(c(0, 16, 40, 90, 200), function(h) {
    a <- list(N0 = burden, hco3 = h,
              hydration = if (burden > 3) HYD_STD else HYD_AGGR)
    if (rasb) { a$rasb_start <- 0; a$rasb_days <- 3 }
    d <- do.call(tls_sim, a)
    s <- tls_summary(d)
    data.frame(urine_pH = tail(d$urinepH, 1), XT_UA = s$XT_UA,
               XT_CAP = s$XT_CAP, Ca_nadir = s$Ca_nadir,
               Cr_ratio = s$Cr_ratio, P_rrt = s$P_rrt, P_sz = s$P_sz)
  }))
  print(res, row.names = FALSE, digits = 3)
  cat(sprintf("  optimum urine pH %.2f (peak Cr ratio %.2f)\n",
              res$urine_pH[which.min(res$Cr_ratio)], min(res$Cr_ratio)))
  invisible(res)
}

## ---------------------------------------------------------------------------
## A4.  ALLOPURINOL LEAD TIME.  AGAINST THE HYPOTHESIS THIS ANALYSIS WAS
## WRITTEN TO EXPRESS: the lead time is worth very little, and what little it
## is worth is set by OXYPURINOL ACCUMULATION rather than by the size of the
## pre-existing urate pool.
##
## Reference: peak Cr ratio 1.26 (start at t=0), 1.23 (-12 h), 1.21 (-24 h),
## 1.19 (-72 h), 1.18 (-168 h).  Twelve hours captures most of it; a week adds
## 0.05.  The arithmetic reason is printed by the function: the pre-existing
## miscible urate pool is ~5 mmol against ~72 mmol released, so emptying it
## early cannot matter much.
## ---------------------------------------------------------------------------
TLS_leadtime <- function(N0 = 3.0) {
  cat("\n=== A4. ALLOPURINOL LEAD TIME ===\n")
  pool <- 0.32 * pget("VEC0")
  load <- N0 * pget("Q_PUR")
  cat(sprintf("  pre-existing miscible urate pool  %6.1f mmol (%.2f g)\n",
              pool, pool * 0.168))
  cat(sprintf("  purine released by lysis          %6.1f mmol (%.2f g)\n",
              load, load * 0.168))
  cat(sprintf("  ratio                             %6.1f : 1\n", load / pool))
  cat("  -> a flux operator started late has almost the same job as one\n")
  cat("     started early, because the pool it failed to empty is small.\n\n")
  res <- do.call(rbind, lapply(list(NA, 0, -12, -24, -48, -72, -120, -168),
                               function(st) {
    a <- list(N0 = N0)
    if (!is.na(st)) a$allo_start <- st
    d <- do.call(tls_sim, a)
    s <- tls_summary(d)
    data.frame(start = ifelse(is.na(st), NA, st),
               XOI_at_0 = approx(d$time, d$XOI, 0)$y,
               UA_peak = s$UA_peak, XAN_peak = max(d$XAN / pget("VEC0")) * 15.21,
               XT_UA = s$XT_UA, XT_XAN = s$XT_XAN,
               Cr_ratio = s$Cr_ratio, P_rrt = s$P_rrt)
  }))
  print(res, row.names = FALSE, digits = 3)
  cat("\n  Watch the XANTHINE column: the better the flux block, the more\n")
  cat("  xanthine there is to precipitate, and xanthine is less soluble than\n")
  cat("  the urate it replaced and is NOT rescued by alkali (pKa 7.7).\n")
  invisible(res)
}

## ---------------------------------------------------------------------------
## A5.  POOL vs FLUX at matched start time.  Reference 4-hour urate change
## (the registration-trial endpoint): nothing +26%, allopurinol 300 +18%,
## febuxostat -7%, rasburicase 0.2 mg/kg -89%.  Reported: rasburicase -86%,
## allopurinol +2% (Goldman 2001).
## ---------------------------------------------------------------------------
TLS_pool_vs_flux <- function(N0 = 3.0) {
  cat("\n=== A5. POOL vs FLUX, both started at t = 0 ===\n")
  arms <- list(
    "nothing"             = list(),
    "allopurinol 300"     = list(allo_start = 0),
    "allopurinol 600"     = list(allo_start = 0, allo_dose = 600),
    "febuxostat 120"      = list(febu_start = 0),
    "rasburicase 0.2 x1"  = list(rasb_start = 0),
    "rasburicase 0.2 x3"  = list(rasb_start = 0, rasb_days = 3),
    "rasburicase 0.15 x1" = list(rasb_start = 0, rasb_mgkg = 0.15),
    "febuxostat + rasb"   = list(febu_start = 0, rasb_start = 0))
  res <- do.call(rbind, lapply(names(arms), function(nm) {
    a <- arms[[nm]]; a$N0 <- N0
    d <- do.call(tls_sim, a); s <- tls_summary(d, nm)
    u0 <- approx(d$time, d$UA_mgdl, 0)$y
    u4 <- approx(d$time, d$UA_mgdl, 4)$y
    w <- d[d$time >= 0 & d$time <= 96, ]
    data.frame(arm = nm, UA_0h = u0, UA_4h = u4,
               pct_4h = 100 * (u4 - u0) / u0, UA_peak = s$UA_peak,
               AUC_0_96 = sum(diff(w$time) * head(w$UA_mgdl, -1)),
               Cr_ratio = s$Cr_ratio, P_rrt = s$P_rrt)
  }))
  print(res, row.names = FALSE, digits = 3)
  invisible(res)
}

## ---------------------------------------------------------------------------
## A5b. RASBURICASE IS A ZERO-ORDER POOL OPERATOR WITH A FIXED CAPACITY.
## Km of urate oxidase for urate (~25 umol/L) is two orders of magnitude below
## the urate concentrations of TLS, so the enzyme runs saturated: a dose buys a
## fixed mmol/h, not a fixed fractional reduction, and that is something a
## tumour can outgrow.
##
## Reference: 0.20 mg/kg = 2.62 mmol/h, covering ~3.2e12 cells; 0.40 mg/kg
## covers ~6.5e12.  Across burdens the 4-hour drop stays impressive (-98% at
## 0.5e12 to -43% at 10e12) while the PEAK urate goes from 5.2 to 94.7 mg/dL.
## The registration endpoint measures the pre-existing pool, not whether the
## drug keeps up with the release flux -- and the kidney outcome depends on the
## second thing.
## ---------------------------------------------------------------------------
TLS_rasburicase_capacity <- function() {
  cat("\n=== A5b. RASBURICASE CAPACITY (zero order in urate) ===\n")
  cap <- do.call(rbind, lapply(c(0.05, 0.10, 0.15, 0.20, 0.30, 0.40), function(dz) {
    cp <- rasburicase_capacity(dz)
    f <- function(n) peak_release_flux(n) - cp
    data.frame(dose_mgkg = dz, C_rasb = dz * pget("BW") / pget("V1_RASB"),
               capacity_mmol_h = cp,
               burden_covered = uniroot(f, c(1e-6, 60))$root)
  }))
  print(cap, row.names = FALSE, digits = 3)
  res <- do.call(rbind, lapply(c(0.5, 1, 3, 6, 10), function(N0) {
    d <- tls_sim(N0 = N0, rasb_start = 0); s <- tls_summary(d)
    u0 <- approx(d$time, d$UA_mgdl, 0)$y
    u4 <- approx(d$time, d$UA_mgdl, 4)$y
    data.frame(burden = N0, UA_4h = u4, pct_4h = 100 * (u4 - u0) / u0,
               UA_peak = s$UA_peak, urate_oxidised = s$H2O2,
               Cr_ratio = s$Cr_ratio)
  }))
  print(res, row.names = FALSE, digits = 3)
  invisible(list(capacity = cap, ode = res))
}

## ---------------------------------------------------------------------------
## A5c. POTASSIUM AND PHOSPHATE ARE SECOND-ORDER SOLUTES.  A1 showed that only
## urate loses the race at intact GFR.  If that is right, the potassium and
## phosphate excursions are consequences of the urate limb's kidney injury, and
## a drug with no potassium or phosphate pharmacology at all should lower them.
## That is a falsifiable prediction of the structure, and rasburicase is the
## test.  Reference (3e12, 2 L/day): K peak 6.07 -> 5.38, PO4 2.49 -> 1.92,
## on rasburicase alone.
## ---------------------------------------------------------------------------
TLS_second_order_solutes <- function() {
  cat("\n=== A5c. POTASSIUM AND PHOSPHATE AS CONSEQUENCES ===\n")
  arms <- list(
    "no prophylaxis, 2 L/day" = list(N0 = 3, hydration = HYD_STD),
    "+ rasburicase only"      = list(N0 = 3, hydration = HYD_STD, rasb_start = 0),
    "+ febuxostat only"       = list(N0 = 3, hydration = HYD_STD, febu_start = -72),
    "no ppx, 6e12 cells"      = list(N0 = 6, hydration = HYD_STD),
    "+ rasb 0.4 x5 only"      = list(N0 = 6, hydration = HYD_STD, rasb_start = 0,
                                     rasb_days = 5, rasb_mgkg = 0.4))
  res <- do.call(rbind, lapply(names(arms), function(nm) {
    s <- tls_summary(do.call(tls_sim, arms[[nm]]), nm)
    s[, c("label", "UA_peak", "K_peak", "PO4_peak", "Ca_nadir",
          "eGFR_nadir", "P_arr")]
  }))
  print(res, row.names = FALSE, digits = 3)
  cat("\n  Potassium and phosphate move in arms whose only intervention is a\n")
  cat("  urate enzyme, so in this model hyperkalaemia in TLS tracks the\n")
  cat("  creatinine rather than the LDH.\n")
  invisible(res)
}

## ---------------------------------------------------------------------------
## A6.  WHAT SURVIVES DELETION OF THE URATE LIMB?  The injury drive is a sum of
## four terms and each arm reports its decomposition, so the claim is computed
## rather than asserted.
##
## Reference: at 3e12 cells the calcium-phosphate share of the peak injury
## drive is 0%; at 6e12 it is 22%, and still 17% after the urate limb has been
## suppressed as far as the enzyme can suppress it.  The phosphate binder does
## nothing about it (Cr ratio 1.15 -> 1.15) for an arithmetic reason: lysis
## delivers 12.1 mmol/h of phosphate to the ECF against 0.95 mmol/h from the
## diet, so a gut binder can reach 7% of the load.  Only removal from the body
## (dialysis) touches the rest.
## ---------------------------------------------------------------------------
TLS_residual <- function() {
  cat("\n=== A6. WHAT SURVIVES DELETION OF THE URATE LIMB? ===\n")
  arms <- list(
    "hydration only 3 L/m2"     = list(N0 = 3),
    "hydration only 2 L/day"    = list(N0 = 3, hydration = HYD_STD),
    "+ rasburicase"             = list(N0 = 3, rasb_start = 0),
    "+ rasb + PO4 binder"       = list(N0 = 3, rasb_start = 0, seve_start = 0),
    "high burden 2 L/day"       = list(N0 = 6, hydration = HYD_STD),
    "high burden + rasb x1"     = list(N0 = 6, hydration = HYD_STD, rasb_start = 0),
    "high burden + rasb x5"     = list(N0 = 6, hydration = HYD_STD, rasb_start = 0,
                                       rasb_days = 5),
    "high burden + rasb 0.4 x5" = list(N0 = 6, hydration = HYD_STD, rasb_start = 0,
                                       rasb_days = 5, rasb_mgkg = 0.4),
    "... + max fluid"           = list(N0 = 6, hydration = HYD_MAX, rasb_start = 0,
                                       rasb_days = 5, rasb_mgkg = 0.4),
    "... + early dialysis"      = list(N0 = 6, hydration = HYD_MAX, rasb_start = 0,
                                       rasb_days = 5, rasb_mgkg = 0.4,
                                       dialysis = c(24, 120, 1.5)))
  res <- do.call(rbind, lapply(names(arms), function(nm) {
    s <- tls_summary(do.call(tls_sim, arms[[nm]]), nm)
    s[, c("label", "XT_UA", "XT_CAP", "Ca_nadir", "Cr_ratio",
          "eGFR_nadir", "P_rrt")]
  }))
  print(res, row.names = FALSE, digits = 3)
  jl <- peak_release_flux(6, q = pget("Q_PO4"))
  cat(sprintf("\n  phosphate entering ECF from lysis  %6.2f mmol/h\n", jl))
  cat(sprintf("  phosphate entering ECF from diet   %6.2f mmol/h\n",
              pget("PO4_INT")))
  cat(sprintf("  share a gut binder can reach       %6.1f %%\n",
              100 * pget("PO4_INT") / (jl + pget("PO4_INT"))))
  cat("  Sequestration is an operator on the INTAKE term, and in acute TLS\n")
  cat("  the intake term is not where the phosphate comes from.\n")
  invisible(res)
}

## ---------------------------------------------------------------------------
## A8.  FLUX-SHAPING: the venetoclax ramp.  Same drug, same target, the same
## eventual kill; only the release RATE differs, because the clearance system
## is a low-pass filter on the release flux.  CLL kinetics (TD 2000 h).
##
## Reference (5e12-cell CLL, 2 L/day): 400 mg from day 1 -> peak release flux
## 2.85 mmol/h, Cr ratio 6.85, P(RRT) 74%, clinical TLS.  The 5-week label ramp
## -> 0.84 mmol/h, Cr ratio 1.63, P(RRT) 1.6%.  All arms lyse 5.0-5.3e12 cells
## and all end with zero residual tumour, so the ramp buys safety without
## giving up kill.  The 2-step and the 5-week ramp are IDENTICAL on every
## safety column, because the peak flux is set by the FIRST dose level -- which
## is why the label puts the monitoring requirement on the 20 mg dose.
## ---------------------------------------------------------------------------
TLS_ramp <- function(N0 = 5.0, tend = 24 * 42) {
  cat("\n=== A8. FLUX-SHAPING: the venetoclax ramp ===\n")
  ramps <- list(
    "400 mg from day 1" = data.frame(time = 0, dose = 400),
    "200 mg from day 1" = data.frame(time = 0, dose = 200),
    "100 mg from day 1" = data.frame(time = 0, dose = 100),
    "50 mg from day 1"  = data.frame(time = 0, dose = 50),
    "2-step 20/50"      = data.frame(time = c(0, 168), dose = c(20, 50)),
    "5-week label ramp" = data.frame(time = c(0, 168, 336, 504, 672),
                                     dose = c(20, 50, 100, 200, 400)),
    "8-week slow ramp"  = data.frame(time = c(0, 168, 336, 504, 672, 840),
                                     dose = c(10, 20, 50, 100, 200, 400)))
  res <- do.call(rbind, lapply(names(ramps), function(nm) {
    d <- tls_sim(N0 = N0, hydration = HYD_STD, chemo = NULL, tend = tend,
                 delta = 0.5, ven = ramps[[nm]],
                 extra_param = list(TD = 2000))
    s <- tls_summary(d, nm)
    data.frame(schedule = nm,
               Jrel_peak = max(pget("Q_PUR") * pget("KLYS") * d$N_LYS),
               K_peak = s$K_peak, PO4_peak = s$PO4_peak, UA_peak = s$UA_peak,
               Cr_ratio = s$Cr_ratio, P_rrt = s$P_rrt, lysed = s$lysed,
               N_end = tail(d$N_VIA, 1),
               TLS = ifelse(s$CTLS, "C", ifelse(s$LTLS, "L", "-")))
  }))
  print(res, row.names = FALSE, digits = 3)
  invisible(res)
}

## ---------------------------------------------------------------------------
## A9.  REDISTRIBUTION IS NOT CLEARANCE.  Total body exchangeable potassium
## (extracellular + parked) is reported so that the difference between moving
## potassium and removing it is visible.
##
## Reference (starting from K 6.07 at 42 h): insulin 10 U gives K 5.47 at 2 h
## (a -0.60 shift, matching Allon 1990) but 6.06 at 12 h -- HIGHER than the
## untreated trajectory's 5.90 -- and total body exchangeable K at 24 h is
## -5 mmol against the untreated arm's -8 mmol.  Repeated insulin is worse
## still: +11 mmol.  The reason is that distal potassium secretion depends on
## the delivered concentration, so parking potassium inside cells switches off
## the only route that was removing it.  Haemodialysis gives K 1.89 at 6 h and
## -26 mmol; the gut binder -22 mmol; furosemide -15 mmol.
## ---------------------------------------------------------------------------
TLS_potassium_rescue <- function() {
  cat("\n=== A9. REDISTRIBUTION vs CLEARANCE ===\n")
  base <- tls_sim(N0 = 3, hydration = HYD_STD)
  i <- which.max(base$K)
  cat(sprintf("  starting point: t = %.0f h, K = %.2f mmol/L, eGFR = %.0f\n",
              base$time[i], base$K[i], base$eGFR[i]))
  y0 <- as.list(base[i, intersect(names(base), names(init(mod)))])
  arms <- list(
    "no rescue"              = list(),
    "insulin 10 U + glucose" = list(insulin = c(0, 1, 900)),
    "salbutamol 20 mg neb"   = list(salb = TRUE),
    "insulin + salbutamol"   = list(insulin = c(0, 1, 900), salb = TRUE),
    "SZC 10 g q8h"           = list(szc_start = 0),
    "furosemide 40 mg q8h"   = list(furo_start = 0),
    "haemodialysis 4 h"      = list(dialysis = c(2, 6, 6)))
  res <- do.call(rbind, lapply(names(arms), function(nm) {
    a <- arms[[nm]]
    salb <- isTRUE(a$salb); a$salb <- NULL
    events <- do.call(tls_events, c(a[intersect(names(a),
      c("szc_start", "furo_start"))], list(tend = 48)))
    if (salb) {
      e2 <- ev(time = 0, amt = 45, cmt = "SALB_C", ii = 6, addl = 3)
      events <- if (is.null(events)) e2 else c(events, e2)
    }
    tg <- seq(0, 48, by = 0.1)
    dat <- data.frame(ID = 1, time = tg,
      RINS = if (is.null(a$insulin)) 0 else
        ifelse(tg >= a$insulin[1] & tg < a$insulin[2], a$insulin[3], 0),
      CL_HD = if (is.null(a$dialysis)) 0 else
        ifelse(tg >= a$dialysis[1] & tg < a$dialysis[2], a$dialysis[3], 0),
      RCHEMO = 0)
    m <- mod %>% param(FLUID = HYD_STD) %>% init(y0)
    o <- if (is.null(events)) m %>% data_set(dat) %>% mrgsim(obsonly = FALSE) else
      m %>% data_set(dat) %>% ev(events) %>% mrgsim(obsonly = FALSE)
    d <- as.data.frame(o)
    tbk <- d$K_EC + d$K_SH
    at <- function(h) approx(d$time, d$K, h)$y
    data.frame(rescue = nm, K_2h = at(2), K_6h = at(6), K_12h = at(12),
               K_24h = at(24), TBK_24h = approx(d$time, tbk, 24)$y,
               dTBK = approx(d$time, tbk, 24)$y - tbk[1],
               P_arr = tail(d$P_ARR, 1))
  }))
  print(res, row.names = FALSE, digits = 3)
  cat("\n  Only the dialysis, binder and diuretic rows change total body K.\n")
  invisible(res)
}

## ---------------------------------------------------------------------------
## A10. THE CALCIUM REFLEX.  Hypocalcaemia in TLS is a CONSEQUENCE of
## calcium-phosphate precipitation, so the calcium given to correct it is a
## reactant in the reaction that caused it.
##
## Reference (3e12 cells): calcium 0 / 1.5 / 3.0 / 6.0 mmol/h gives Ca x PO4
## 59 / 73 / 80 / 89 mg2/dL2, renal CaP 0.00 / 2.32 / 11.03 / 10.75 mmol,
## peak Cr ratio 1.55 / 1.88 / 3.19 / 7.89 and P(RRT) 0.9% / 0.9% / 16% / 82%.
## ---------------------------------------------------------------------------
TLS_calcium_reflex <- function(N0 = 3.0) {
  cat("\n=== A10. THE CALCIUM REFLEX ===\n")
  res <- do.call(rbind, lapply(c(0, 0.5, 1.5, 3.0, 6.0), function(ca) {
    d <- tls_sim(N0 = N0, ca_inf = ca); s <- tls_summary(d)
    data.frame(Ca_inf_mmol_h = ca, Ca_nadir = s$Ca_nadir,
               CaxPO4_max = max(d$CaxPO4), XT_CAP = s$XT_CAP,
               Cr_ratio = s$Cr_ratio, eGFR_nadir = s$eGFR_nadir,
               P_rrt = s$P_rrt, P_sz = s$P_sz)
  }))
  print(res, row.names = FALSE, digits = 3)
  cat("\n  The seizure column is the reason the reflex exists; the renal\n")
  cat("  replacement column is its price.\n")
  invisible(res)
}

## ---------------------------------------------------------------------------
## A11. IS THE LOOP A SWITCH OR A KNEE?  AGAINST what a positive-feedback loop
## might be expected to produce, this one is NOT bistable: loop gain stays
## below 1 and the response is a steep but single-valued function with a knee
## near 3.5e12 cells.  Seeding the kidney with crystal at t = 0 shifts the
## trajectory but every seed returns to the same eGFR by day 14.
##
## Reference: nephron-reserve axis Cr ratio 1.55 (NEPH 1.00) -> 4.98 (0.30);
## insult axis 1.08 (1e12) -> 3.06 (3e12) -> 8.48 (4e12) -> 16.50 (6e12);
## crystal-seed axis 1.14 (no seed) -> 10.25 (40 mmol) but eGFR at 14 days is
## 119 / 119 / 119 / 119 / 119 / 106 -- one attractor.
## ---------------------------------------------------------------------------
TLS_bistability <- function() {
  cat("\n=== A11. SWITCH OR KNEE? ===\n")
  cat("\n  nephron reserve axis (burden fixed at 3e12):\n")
  a <- do.call(rbind, lapply(c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3), function(n) {
    d <- tls_sim(N0 = 3, neph0 = n); s <- tls_summary(d)
    data.frame(NEPH_0 = n, eGFR_0 = d$eGFR[1], Cr_ratio = s$Cr_ratio,
               eGFR_nadir = s$eGFR_nadir, NEPH_14d = s$NEPH_end,
               XT_max = s$XT_max, P_rrt = s$P_rrt)
  }))
  print(a, row.names = FALSE, digits = 3)
  cat("\n  crystal-seed axis (hysteresis test, burden fixed at 2e12):\n")
  b <- do.call(rbind, lapply(c(0, 2, 5, 10, 20, 40), function(sd) {
    d <- tls_sim(N0 = 2, xtua0 = sd); s <- tls_summary(d)
    data.frame(seed_mmol = sd, Cr_ratio = s$Cr_ratio,
               eGFR_nadir = s$eGFR_nadir, XT_14d = tail(d$XTtot, 1),
               eGFR_14d = tail(d$eGFR, 1), P_rrt = s$P_rrt)
  }))
  print(b, row.names = FALSE, digits = 3)
  cat("\n  insult axis (nephron mass intact, 2 L/day):\n")
  cc <- do.call(rbind, lapply(c(1, 1.5, 2, 2.5, 3, 3.5, 4, 6), function(N0) {
    s <- tls_summary(tls_sim(N0 = N0, hydration = HYD_STD))
    data.frame(burden = N0, Cr_ratio = s$Cr_ratio, eGFR_nadir = s$eGFR_nadir,
               XT_max = s$XT_max, NEPH_14d = s$NEPH_end, P_rrt = s$P_rrt)
  }))
  print(cc, row.names = FALSE, digits = 3)
  invisible(list(reserve = a, seed = b, insult = cc))
}

## ---------------------------------------------------------------------------
## A12. CAIRO-BISHOP AS A DETECTION FILTER.  A NULL RESULT, reported as one:
## the >=2-of-4 rule is ROBUST to blood-draw frequency in this model, because
## three of the four analytes have multi-day time constants and even the
## potassium peak is broad (a q24h draw misses only 0.06 mmol/L of it).  What
## the analysis does show is WHICH criteria fire, and that the 8 mg/dL urate
## cut-off coincides with UA_crit from A1 -- the definition's urate criterion
## is a crystallisation detector, and its other three criteria are consequence
## detectors.
## ---------------------------------------------------------------------------
TLS_detection <- function() {
  cat("\n=== A12. WHAT THE CAIRO-BISHOP DEFINITION DETECTS ===\n")
  scen <- list(
    "high burden, no ppx"  = list(N0 = 3, hydration = HYD_STD),
    "high burden + rasb"   = list(N0 = 3, rasb_start = 0),
    "high burden + bundle" = list(N0 = 3, rasb_start = 0, seve_start = 0,
                                  hydration = HYD_MAX),
    "moderate + allo -72h" = list(N0 = 1, allo_start = -72),
    "low burden"           = list(N0 = 0.2, allo_start = -72))
  uc <- ua_critical(HYD_AGGR - pget("INSENS")) * 16.81
  res <- do.call(rbind, lapply(names(scen), function(nm) {
    d <- do.call(tls_sim, scen[[nm]])
    cb <- cairo_bishop(d)
    samp <- sapply(c(4, 8, 12, 24, 48), function(h) {
      x <- cairo_bishop(d, sample_h = h)
      ifelse(x$CTLS, "C", ifelse(x$LTLS, "L", "-"))
    })
    data.frame(scenario = nm, urate = cb$urate, K = cb$K, PO4 = cb$PO4,
               Ca = cb$Ca, n = cb$n, UA_peak = max(d$UA_mgdl), UA_crit = uc,
               q4h = samp[1], q8h = samp[2], q12h = samp[3],
               q24h = samp[4], q48h = samp[5])
  }))
  print(res, row.names = FALSE, digits = 3)
  invisible(res)
}

## ---------------------------------------------------------------------------
## A13. THE OXIDANT LOAD SCALES WITH THE INDICATION.  Urate oxidase makes one
## H2O2 per urate destroyed, so the patient with the largest tumour receives
## the largest peroxide load: the toxicity is proportional to the reason for
## giving the drug.  Reference (G6PD-deficient): methaemoglobin 8.2% at 0.3e12
## rising to 18.0% at 6e12, haemoglobin nadir 12.2 -> 11.0 g/dL.
## ---------------------------------------------------------------------------
TLS_oxidant_load <- function() {
  cat("\n=== A13. OXIDANT LOAD vs TUMOUR BURDEN ===\n")
  res <- do.call(rbind, lapply(c(0.3, 1, 3, 6), function(N0)
    do.call(rbind, lapply(c(1, 0), function(g) {
      s <- tls_summary(tls_sim(N0 = N0, rasb_start = 0,
                               extra_param = list(G6PD = g)))
      data.frame(burden = N0, G6PD = ifelse(g == 1, "normal", "deficient"),
                 urate_oxidised = s$H2O2, MetHb_pct = s$MetHb,
                 Hb_nadir = s$Hb_nadir, UA_peak = s$UA_peak)
    }))))
  print(res, row.names = FALSE, digits = 3)
  invisible(res)
}

## ---------------------------------------------------------------------------
## A14. THE LOOP BITES THE DRUG TOO.  Oxypurinol is renally cleared and
## febuxostat is not, so the kidney injury allopurinol is given to prevent
## changes allopurinol's own exposure.  Reference: oxypurinol 5.8 mg/L at day 7
## in the low-burden arm and 30.3 mg/L in the high-burden arm whose eGFR has
## fallen to 16 -- a 5-fold accumulation, with XO inhibition rising from 0.84
## to 0.97.  This is the model's account of why allopurinol needs renal dose
## reduction in exactly the patients whose TLS risk is highest.
## ---------------------------------------------------------------------------
TLS_oxypurinol <- function() {
  cat("\n=== A14. OXYPURINOL ACCUMULATION IN EVOLVING AKI ===\n")
  arms <- list(
    list("allopurinol, low burden", list(N0 = 0.2, allo_start = -72), "C_OXY"),
    list("allopurinol, high burden",
         list(N0 = 6, allo_start = -72, hydration = HYD_STD), "C_OXY"),
    list("febuxostat, high burden",
         list(N0 = 6, febu_start = -72, hydration = HYD_STD), "C_FEBU"))
  res <- do.call(rbind, lapply(arms, function(a) {
    d <- do.call(tls_sim, a[[2]])
    at <- function(h, k) approx(d$time, d[[k]], h)$y
    data.frame(arm = a[[1]], C_d1 = at(24, a[[3]]), C_d3 = at(72, a[[3]]),
               C_d7 = at(168, a[[3]]), XOI_d7 = at(168, "XOI"),
               eGFR_d7 = at(168, "eGFR"))
  }))
  print(res, row.names = FALSE, digits = 3)
  invisible(res)
}

## ---------------------------------------------------------------------------
## A15. TRIAL LEDGER.  Where the model matches, and where it does not.
## ---------------------------------------------------------------------------
TLS_trial_ledger <- function() {
  cat("\n=== A15. TRIAL LEDGER ===\n")
  d_r <- tls_sim(N0 = 1, rasb_start = 0, rasb_days = 5)
  d_a <- tls_sim(N0 = 1, allo_start = 0)
  p4 <- function(d) {
    u0 <- approx(d$time, d$UA_mgdl, 0)$y
    100 * (approx(d$time, d$UA_mgdl, 4)$y - u0) / u0
  }
  auc <- function(d) { w <- d[d$time >= 0 & d$time <= 96, ]
    sum(diff(w$time) * head(w$UA_mgdl, -1)) }
  rows <- list(
    c("rasburicase urate change at 4 h", sprintf("%+.0f%%", p4(d_r)),
      "-86% (Goldman 2001 Blood)"),
    c("allopurinol urate change at 4 h", sprintf("%+.0f%%", p4(d_a)),
      "+2% (Goldman 2001)"),
    c("urate AUC0-96 rasb / allo", sprintf("%.2fx", auc(d_r) / auc(d_a)),
      "0.39x (Goldman 2001) - model OVER-separates"),
    c("urate held < 8 mg/dL on rasburicase",
      ifelse(max(d_r$UA_mgdl) < 8, "yes", "no"), "87% (Cortes 2010 JCO)"),
    c("urate held < 8 mg/dL on allopurinol",
      ifelse(max(d_a$UA_mgdl) < 8, "yes", "no"), "66% (Cortes 2010)"),
    c("dialysis, high risk, full bundle",
      sprintf("%.0f%%", 100 * tls_summary(tls_sim(N0 = 3, rasb_start = 0,
        seve_start = 0, hydration = HYD_MAX))$P_rrt),
      "1.5-5% overall (Coiffier 2008 JCO)"),
    c("dialysis, high risk, unprophylaxed",
      sprintf("%.0f%%", 100 * tls_summary(tls_sim(N0 = 3,
        hydration = HYD_STD))$P_rrt),
      "up to ~30% in historical Burkitt series"),
    c("LDH peak, high-burden Burkitt",
      sprintf("%.0f U/L", tls_summary(tls_sim(N0 = 3,
        hydration = HYD_STD))$LDH_peak),
      "commonly > 2x ULN, often > 5000"),
    c("insulin 10 U potassium shift", "see TLS_potassium_rescue()",
      "-0.6 to -1.0 mmol/L at 1-4 h (Allon 1990)"))
  for (r in rows)
    cat(sprintf("  %-38s model %10s   reported %s\n", r[1], r[2], r[3]))
  invisible(rows)
}

## ---------------------------------------------------------------------------
## A16. THE 12 SHIPPED SCENARIOS.  Reference peak Cr ratios: S1 1.00, S2 3.06,
## S3 1.55, S4 1.26, S5 1.19, S6 1.15, S7 1.15, S8 1.18, S9 1.04, S10 1.16,
## S11 1.06, S12 1.15.
## ---------------------------------------------------------------------------
TLS_run_scenarios <- function() {
  cat("\n=== A16. THE 12 SHIPPED SCENARIOS ===\n")
  scen <- list(
    "S1  low-risk solid tumour"        = list(N0 = 0.05, hydration = HYD_STD),
    "S2  high burden, fluids only"     = list(N0 = 3, hydration = HYD_STD),
    "S3  + aggressive hydration"       = list(N0 = 3),
    "S4  + allopurinol t=0"            = list(N0 = 3, allo_start = 0),
    "S5  + allopurinol -72 h"          = list(N0 = 3, allo_start = -72),
    "S6  + rasburicase 0.2 mg/kg"      = list(N0 = 3, rasb_start = 0),
    "S7  + rasb, max fluid, no alkali" = list(N0 = 3, rasb_start = 0,
                                              hydration = HYD_MAX),
    "S8  + alkalinisation pH 7.5"      = list(N0 = 3, hco3 = 45),
    "S9  steroid prephase 5 d"         = list(N0 = 3, pred_start = -120,
                                              pred_days = 5),
    "S10 full bundle"                  = list(N0 = 3, rasb_start = 0,
                                              seve_start = 0, allo_start = -72,
                                              hydration = HYD_MAX,
                                              furo_start = 24),
    "S11 CKD stage 3 host"             = list(N0 = 3, rasb_start = 0,
                                              allo_start = -72, neph0 = 0.45),
    "S12 G6PD deficient + rasb"        = list(N0 = 3, rasb_start = 0,
                                              extra_param = list(G6PD = 0)))
  res <- do.call(rbind, lapply(names(scen), function(nm) {
    s <- tls_summary(do.call(tls_sim, scen[[nm]]), nm)
    s[, c("label", "UA_peak", "K_peak", "PO4_peak", "Ca_nadir", "Cr_ratio",
          "eGFR_nadir", "LDH_peak", "lysed", "lysed_pre", "P_rrt",
          "LTLS", "CTLS")]
  }))
  print(res, row.names = FALSE, digits = 3)
  cat("\n  lysed_pre is non-zero only for the prephase arm, which is the\n")
  cat("  point of a prephase: it moves lysis earlier and slower.\n")
  invisible(res)
}

## ---------------------------------------------------------------------------
## A17. Numerical hygiene.
## ---------------------------------------------------------------------------
TLS_mass_balance <- function(N0 = 3) {
  cat("\n=== A17. NUMERICAL HYGIENE ===\n")
  d <- tls_sim(N0 = N0, rasb_start = 0)
  cat(sprintf("  purine released by lysis        %8.2f mmol\n",
              pget("Q_PUR") * N0))
  cat(sprintf("  K released by lysis             %8.2f mmol\n",
              pget("Q_K") * N0))
  cat(sprintf("  PO4 released by lysis           %8.2f mmol\n",
              pget("Q_PO4") * N0))
  cat(sprintf("  cumulative urate oxidised       %8.2f mmol\n",
              tail(d$H2O2, 1)))
  st <- names(init(mod))
  bad <- st[sapply(st, function(s) min(d[[s]]) < -1e-6)]
  cat("  negative states: ", ifelse(length(bad) == 0, "none",
                                    paste(bad, collapse = ", ")), "\n")
  invisible(d)
}

TLS_run_all <- function() {
  TLS_baseline(); TLS_race(); TLS_operator_decomposition()
  TLS_ph_optimum(3, FALSE); TLS_ph_optimum(3, TRUE); TLS_ph_optimum(6, FALSE)
  TLS_leadtime(); TLS_pool_vs_flux(); TLS_rasburicase_capacity()
  TLS_second_order_solutes(); TLS_residual(); TLS_ramp()
  TLS_potassium_rescue(); TLS_calcium_reflex(); TLS_bistability()
  TLS_detection(); TLS_oxidant_load(); TLS_oxypurinol()
  TLS_trial_ledger(); TLS_run_scenarios(); TLS_mass_balance()
  invisible(NULL)
}

message("tls_mrgsolve_model.R loaded: 48 ODEs, 12 scenarios, 18 analyses.")
message("Start with TLS_baseline(), TLS_race(), TLS_operator_decomposition().")
message("Every reference number in the comments comes from ",
        "tls_reference_check.py, which integrates the same system in scipy.")
