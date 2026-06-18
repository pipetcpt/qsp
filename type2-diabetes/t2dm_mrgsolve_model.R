## =============================================================================
##  Type 2 Diabetes Mellitus — mrgsolve QSP Model
##  Disease: T2DM  |  Version: 1.0  |  Date: 2026-06-16
##
##  Model structure:
##    - 7 drug PK compartments (Metformin, Empagliflozin, Semaglutide,
##      Sitagliptin, Glimepiride, Insulin Degludec, Pioglitazone)
##    - Glucose-insulin minimal model with tissue compartment (Bergman ext.)
##    - Glucagon dynamics
##    - β-cell mass/function trajectories (adapted from Topp et al., 2000)
##    - Incretin (GLP-1) compartment with DPP-4 degradation
##    - Hepatic & peripheral insulin resistance indices
##    - FFA / adiposity / body weight
##    - HbA1c, eGFR, UACR endpoints
##
##  Treatment scenarios (7):
##    1. No treatment (diet/exercise baseline)
##    2. Metformin 1000 mg BID
##    3. Metformin + Empagliflozin 10 mg QD
##    4. Metformin + Semaglutide 1 mg SC weekly
##    5. Triple: Metformin + Empagliflozin + Semaglutide
##    6. Metformin + Insulin Degludec 20 U QD
##    7. Metformin + Sitagliptin 100 mg QD (DPP-4i)
##
##  Parameter calibration notes:
##    - Glucose-insulin: Bergman minimal model (Bergman 1989, Diabetes)
##    - β-cell mass: Topp et al. (2000), J Theor Biol 244:501
##    - SGLT2i PD: Ferrannini et al. (2012), Diabetes Care
##    - GLP-1RA PK: Lau et al. (2015), Pharm Res (semaglutide)
##    - Metformin PK: Gong et al. (2012), Clin Pharmacokinet
##    - DECLARE-TIMI 58 (Wiviott 2019): empagliflozin → HbA1c ↓0.7%
##    - SUSTAIN-6 (Marso 2016): semaglutide → HbA1c ↓1.4%, wt ↓4.5 kg
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

# =============================================================================
# MODEL DEFINITION
# =============================================================================
t2dm_model_code <- '
$PROB T2DM QSP Model — Multi-Drug PK/PD

$PARAM
// ---- Patient baseline ----
BW      = 90,    // body weight (kg)
HbA1c0  = 9.0,   // baseline HbA1c (%)
Gp0     = 180,   // baseline plasma glucose (mg/dL)
Ip0     = 20,    // baseline plasma insulin (mU/L)
Gc0     = 150,   // baseline glucagon (pg/mL)
eGFR0   = 80,    // baseline eGFR (mL/min/1.73m2)
UACR0   = 50,    // baseline UACR (ug/mg)

// ---- Drug dosing switches (1=on, 0=off) ----
USE_MET  = 0,    // Metformin
USE_EMPA = 0,    // Empagliflozin
USE_SEMA = 0,    // Semaglutide
USE_DPP4 = 0,    // Sitagliptin
USE_SU   = 0,    // Glimepiride
USE_INS  = 0,    // Insulin Degludec
USE_PIOG = 0,    // Pioglitazone

// ---- Metformin PK (Gong 2012, Clin Pharmacokinet) ----
ka_met  = 1.2,   // absorption rate (h-1), F=55%
F_met   = 0.55,  // bioavailability
Vc_met  = 300,   // central vol (L/70kg)
k12_met = 0.15,  // central→peripheral (h-1)
k21_met = 0.08,  // peripheral→central (h-1)
CL_met  = 35,    // clearance (L/h)
MW_met  = 165.6, // molecular weight

// ---- Empagliflozin PK (Macha 2013, Clin Pharmacokinet) ----
ka_empa = 0.60,  // absorption (h-1), F=86%
F_empa  = 0.86,
Vc_empa = 73,    // (L/70kg)
CL_empa = 10.6,

// ---- Semaglutide PK (Lau 2015, Pharm Res) ----
ka_sema = 0.0085,// SC absorption (h-1), t½abs~5d
F_sema  = 0.89,
Vc_sema = 12.5,  // (L)
CL_sema = 0.053, // very low CL (t½~168h)

// ---- Sitagliptin PK (He 2009, Br J Clin Pharmacol) ----
ka_dpp4 = 0.80,
F_dpp4  = 0.87,
Vc_dpp4 = 198,
CL_dpp4 = 12.4,

// ---- Glimepiride PK (Massi-Benedetti 1996) ----
ka_su   = 0.50,
F_su    = 1.00,
Vc_su   = 12.6,
CL_su   = 3.1,

// ---- Insulin Degludec PK (Kurtzhals 2011) ----
ka_ins  = 0.048, // from SC hexamer depot (h-1)
F_ins   = 0.91,
Vc_ins  = 8.0,
CL_ins  = 0.96,

// ---- Pioglitazone PK (Eckland 2000) ----
ka_piog = 0.70,
F_piog  = 0.83,
Vc_piog = 89,
CL_piog = 5.6,

// ---- Glucose-Insulin Dynamics (Bergman + extensions) ----
SI       = 8e-4, // insulin sensitivity (dL/mU/h, normal ~1e-3)
Sg       = 0.01, // glucose effectiveness (h-1)
EGP0     = 2.4,  // basal EGP (mg/kg/min)
Rd0      = 2.4,  // basal Rd (mg/kg/min)
p2       = 0.05, // remote compartment transfer (h-1)
Vg       = 1.5,  // glucose distribution vol (dL/kg → scaled by BW)
Vi       = 0.05, // insulin vol (L/kg)

// ---- Insulin secretion ----
beta_sens = 0.6, // GSIS sensitivity (β-cell function)
beta_M0  = 1.0,  // normalized baseline β-cell mass
k_prolif = 1e-4, // β-cell proliferation constant (d-1)
k_apop   = 1e-3, // β-cell apoptosis rate (d-1) — elevated in T2DM
Gth      = 90,   // glucose threshold for GSIS (mg/dL)
phi_max  = 900,  // max insulin secretion rate (mU/L/h)

// ---- Glucagon dynamics ----
kout_Gc  = 0.3,  // glucagon turnover (h-1)
Gc_Gp50  = 100,  // Gp at 50% glucagon suppression (mg/dL)
EGP_Gc   = 0.015,// glucagon effect on EGP (mg/kg/min per pg/mL above basal)

// ---- GLP-1 dynamics ----
kin_GLP1 = 15,   // basal GLP-1 production (pmol/L/h)
kout_GLP1= 4.0,  // GLP-1 elimination (h-1)
kDPP4    = 3.5,  // DPP-4 degradation rate (h-1)
GLP1_Gp  = 0.05, // GLP-1 response to meal glucose (pmol/L per mg/dL)

// ---- DPP-4 inhibition (Sitagliptin) ----
Imax_dpp4 = 0.80,// max DPP-4 inhibition
IC50_dpp4 = 0.10,// μM (Sitagliptin)

// ---- SGLT2i PD (Empagliflozin) ----
TmG_base  = 340, // max tubular glucose reabsorption (mg/min/1.73m2)
RGT0      = 180, // renal glucose threshold (mg/dL)
Imax_empa = 0.55,// max reduction in TmG (Ferrannini 2012)
IC50_empa = 30,  // ng/mL (empagliflozin)
GFR_val   = 100, // individual GFR (mL/min) for UGE calc

// ---- SU PD ----
Emax_su   = 2.5, // fold-increase in insulin secretion (SU max)
EC50_su   = 50,  // ng/mL Glimepiride
n_su      = 1.5, // Hill coefficient

// ---- GLP-1RA / Semaglutide PD ----
Emax_sema_ins = 0.6, // max increase in GSIS (fold)
EC50_sema    = 5,    // nmol/L semaglutide
Emax_sema_wt = 4.5,  // max weight loss (kg) from GLP-1RA at 1yr
ksema_wt     = 0.004,// weight loss rate constant (h-1)

// ---- Insulin Resistance dynamics ----
IR_H0    = 1.0,  // hepatic IR index (1=normal, >1=resistant)
IR_P0    = 1.0,  // peripheral IR index
k_IR_FFA = 0.005,// FFA-driven IR increase rate
k_IR_rec = 0.001,// natural IR recovery rate

// ---- FFA / adiposity ----
FFA0     = 0.6,  // baseline plasma FFA (mmol/L)
kFFA_rel = 0.15, // lipolysis rate
kFFA_up  = 0.20, // peripheral FFA uptake/clearance
FFA_Ip50 = 15,   // Ip for 50% lipolysis inhibition (mU/L)
BW_loss_rate = 0.0,// body weight change from lifestyle (kg/day)

// ---- HbA1c kinetics ----
kHbA1c   = 0.0084,// HbA1c equilibration rate (h-1, t½≈60 days)
HbA1c_ss = 0.165, // HbA1c-to-mean glucose slope (% per mg/dL above 90)

// ---- Renal / complication dynamics ----
k_eGFR_decline = 2e-5,  // eGFR decline rate per unit HbA1c excess (mL/min/day)
k_UACR_rise    = 0.004, // UACR rise rate per unit HbA1c excess
k_eGFR_empa    = 1.5e-5,// nephroprotective effect of SGLT2i on eGFR decline
k_UACR_empa    = 0.003, // SGLT2i → UACR reduction

// ---- PPARg / Pioglitazone PD ----
Emax_piog_IR  = 0.35,  // max IR reduction (35%)
EC50_piog     = 200,   // ng/mL pioglitazone
Emax_piog_wt  = 3.0,   // max weight GAIN (kg, fluid retention)

// ---- Body weight ----
BW_rate  = 0     // net weight change rate (kg/day, baseline)

$CMT
// Drug PK
MET_GUT MET_C MET_P          // Metformin (gut, central, peripheral)
EMPA_C                        // Empagliflozin plasma
SEMA_SC SEMA_C                // Semaglutide SC depot, plasma
DPP4I_C                       // Sitagliptin plasma
SU_C                          // Glimepiride plasma
INS_SC INS_C                  // Insulin degludec SC depot, plasma
PIOG_C                        // Pioglitazone plasma

// Glucose-Insulin
Gp Gt Ip X_action             // Plasma glucose, tissue glucose, insulin, remote action

// Endocrine
Gc                             // Glucagon (pg/mL)
GLP1                           // Active GLP-1 (pmol/L)
beta_mass                      // β-cell mass (normalized)

// Intermediate
IR_H IR_P                      // Hepatic and peripheral insulin resistance
FFA                            // Plasma FFA (mmol/L)
BW_t                           // Body weight (kg)

// Endpoints
HbA1c_cmpt                     // HbA1c (%)
eGFR_cmpt                      // eGFR (mL/min/1.73m2)
UACR_cmpt                      // UACR (ug/mg creatinine)

$INIT
MET_GUT = 0, MET_C = 0, MET_P = 0,
EMPA_C = 0,
SEMA_SC = 0, SEMA_C = 0,
DPP4I_C = 0,
SU_C = 0,
INS_SC = 0, INS_C = 0,
PIOG_C = 0,
Gp = 180, Gt = 180, Ip = 20, X_action = 0.016,
Gc = 150,
GLP1 = 5,
beta_mass = 0.7,         // T2DM: ~50-70% of normal β-cell mass
IR_H = 2.5,              // elevated hepatic IR
IR_P = 2.0,              // elevated peripheral IR
FFA = 0.8,               // elevated FFA in T2DM
BW_t = 90,               // baseline body weight
HbA1c_cmpt = 9.0,
eGFR_cmpt = 80,
UACR_cmpt = 50

$MAIN
// ---- Plasma drug concentrations (ng/mL or nmol/L as noted) ----
double Cmet   = MET_C / Vc_met * 1000;   // ng/mL
double Cempa  = EMPA_C / Vc_empa * 1000; // ng/mL
double Csema  = SEMA_C / Vc_sema;        // nmol/L (MW~4114 → approx nmol/L)
double Cdpp4  = DPP4I_C / Vc_dpp4 * 1000;// ng/mL
double Csu    = SU_C / Vc_su * 1000;     // ng/mL
double Cins   = INS_C / Vc_ins * 1000;   // ng/mL
double Cpiog  = PIOG_C / Vc_piog * 1000; // ng/mL

// ---- Drug effect functions ----
// Metformin: AMPK → reduces EGP (max 30%), mild Rd increase
double Emet_EGP  = USE_MET * 0.30 * Cmet / (Cmet + 800);
double Emet_Rd   = USE_MET * 0.12 * Cmet / (Cmet + 800);

// Empagliflozin: SGLT2 inhibition → UGE increase
double E_empa = USE_EMPA * Imax_empa * Cempa / (Cempa + IC50_empa * 1000);
double TmG    = TmG_base * (1 - E_empa);   // reduced reabsorption capacity
double UGE    = fmax(0.0, (GFR_val * Gp/100 - TmG) / (24*60)); // mg/min → mg/h per dL
// Urinary glucose excretion (mg/dL plasma equivalent loss per h)
double UGE_loss = fmax(0.0, GFR_val * Gp / 100 - TmG) * 60 / BW_t; // mg/kg/h

// Semaglutide: GLP-1R agonism — amplifies insulin, suppresses glucagon, reduces weight
double E_sema_ins = USE_SEMA * Emax_sema_ins * pow(Csema,1.5) / (pow(Csema,1.5) + pow(EC50_sema,1.5));
double E_sema_Gc  = USE_SEMA * 0.40 * Csema / (Csema + EC50_sema);
double BW_sema_effect = USE_SEMA * Emax_sema_wt * (1 - exp(-ksema_wt * SOLVERTIME));

// DPP-4 inhibition: prevents GLP-1 degradation (extends GLP-1 half-life ~2x)
double DPP4_inh = USE_DPP4 * Imax_dpp4 * Cdpp4 / (Cdpp4 + IC50_dpp4 * 1000);
double keff_DPP4 = kDPP4 * (1 - DPP4_inh);

// Sulfonylurea: increases insulin secretion via SUR1 closure
double E_su = USE_SU * Emax_su * pow(Csu, n_su) / (pow(Csu, n_su) + pow(EC50_su, n_su));

// Pioglitazone: PPARg → reduces IR, increases adiponectin
double E_piog_IR = USE_PIOG * Emax_piog_IR * Cpiog / (Cpiog + EC50_piog);
double BW_piog_gain = USE_PIOG * Emax_piog_wt * (1 - exp(-0.003 * SOLVERTIME));

// ---- GLP-1 effective level (endogenous + semaglutide contribution) ----
double GLP1_total = GLP1 + USE_SEMA * Csema * 50; // sema in pmol/L equiv

// ---- Insulin secretion (GSIS) ----
double GSIS_base = phi_max * beta_mass * beta_sens
                   * pow(fmax(0.0, Gp - Gth), 1.5)
                   / (pow(fmax(0.0, Gp - Gth), 1.5) + pow(90.0, 1.5));
double GSIS_incretin = GSIS_base * (1 + 0.35 * GLP1_total/(GLP1_total + 10.0)
                                    + E_sema_ins);
double GSIS_su       = GSIS_incretin * (1 + E_su);
double dIp_sec       = GSIS_su;  // total secretion into portal

// ---- EGP: suppressed by insulin (via X_action), glucagon drives it up ----
double EGP_ins_supp  = 1.0 / (1.0 + 3.0 * X_action);
double EGP_Gc_drive  = 1.0 + EGP_Gc * fmax(0.0, Gc - Gc0);
double EGP_IR_drive  = IR_H / 2.5; // elevated with hepatic IR
double EGP_val       = EGP0 * EGP_ins_supp * EGP_Gc_drive * EGP_IR_drive
                        * (1 - Emet_EGP);
// cap to physiological range
EGP_val = fmax(0.5, fmin(EGP_val, 8.0));

// ---- Rd: insulin-stimulated glucose disposal ----
double SI_eff  = SI / IR_P;      // peripheral IR reduces SI
double Rd_ins  = SI_eff * X_action * Gt;
double Rd_val  = Rd0 + Rd_ins + Emet_Rd * Rd0;

// ---- Glucagon equation parameters ----
double Gc_ss   = Gc0 * (Gc_Gp50 / fmax(Gp, 50.0)) * (1 - E_sema_Gc)
                * (1.0 / (1.0 + 0.02 * fmax(0.0, Ip - Ip0)));

// ---- β-cell mass dynamics ----
double beta_prolif_rate = k_prolif * beta_mass * fmax(0.0, Gp - 90) / 90.0;
double beta_apop_rate   = k_apop * beta_mass
                          * (1.0 + 0.5 * fmax(0.0, FFA - 0.6))
                          * (1.0 + 0.3 * fmax(0.0, Gp - 180.0)/180.0);
// GLP-1 (incretin) protects β-cell
double GLP1_beta_prot   = 1.0 - 0.3 * GLP1_total / (GLP1_total + 10.0);
beta_apop_rate *= GLP1_beta_prot;

// ---- FFA: anti-lipolytic effect of insulin ----
double Ip_antilipol  = 1.0 / (1.0 + Ip / FFA_Ip50);
double FFA_release_r = kFFA_rel * Ip_antilipol * BW_t / 90.0;
double FFA_uptake_r  = kFFA_up * FFA;

// ---- IR dynamics ----
// Hepatic IR driven by FFA, glucotoxicity, fat accumulation
double IR_H_drive = k_IR_FFA * fmax(0.0, FFA - 0.5) + k_IR_FFA * fmax(0.0, Gp - 150)/300.0;
double IR_H_recov = k_IR_rec * (IR_H - 1.0);
// Metformin reduces hepatic IR
double IR_H_met   = USE_MET * 0.002 * Emet_EGP;
// Pioglitazone reduces both IR_H and IR_P
double dIRH_piog  = E_piog_IR * 0.003 * IR_H;
double dIRP_piog  = E_piog_IR * 0.003 * IR_P;

// ---- Body weight dynamics ----
double BW_UGE_loss   = UGE_loss / 40.0 * 0.001;   // caloric loss from glucosuria (kg/d)
double BW_sema_rate  = USE_SEMA * ksema_wt * (Emax_sema_wt - BW_sema_effect) / 24.0;
double BW_piog_rate  = USE_PIOG * 0.003 * (Emax_piog_wt - BW_piog_gain) / 24.0;
double dBW           = BW_rate - BW_UGE_loss - BW_sema_rate + BW_piog_rate;

// ---- HbA1c equilibration to current mean glucose ----
double Gp_eq         = fmax(Gp, 70.0);
double HbA1c_target  = 5.0 + HbA1c_ss * fmax(0.0, Gp_eq - 90.0);
double dHbA1c        = kHbA1c * (HbA1c_target - HbA1c_cmpt);

// ---- eGFR decline ----
double HbA1c_excess  = fmax(0.0, HbA1c_cmpt - 7.0);
double eGFR_decline  = k_eGFR_decline * HbA1c_excess * 24.0; // per day
double eGFR_protect  = USE_EMPA * k_eGFR_empa * 24.0 * E_empa;
double deGFR         = -(eGFR_decline - eGFR_protect);

// ---- UACR dynamics ----
double UACR_drive    = k_UACR_rise * HbA1c_excess * 24.0;
double UACR_protect  = USE_EMPA * k_UACR_empa * 24.0 * E_empa;
double dUACR         = UACR_drive - UACR_protect;

$ODE
// ===== Drug PK ODEs =====
// Metformin
dxdt_MET_GUT = -ka_met * MET_GUT;
dxdt_MET_C   = ka_met * MET_GUT - (CL_met/Vc_met + k12_met) * MET_C + k21_met * MET_P;
dxdt_MET_P   = k12_met * MET_C - k21_met * MET_P;

// Empagliflozin
dxdt_EMPA_C  = -CL_empa/Vc_empa * EMPA_C;

// Semaglutide
dxdt_SEMA_SC = -ka_sema * SEMA_SC;
dxdt_SEMA_C  = ka_sema * SEMA_SC - CL_sema/Vc_sema * SEMA_C;

// Sitagliptin
dxdt_DPP4I_C = -CL_dpp4/Vc_dpp4 * DPP4I_C;

// Glimepiride
dxdt_SU_C    = -CL_su/Vc_su * SU_C;

// Insulin Degludec
dxdt_INS_SC  = -ka_ins * INS_SC;
dxdt_INS_C   = ka_ins * INS_SC - CL_ins/Vc_ins * INS_C;

// Pioglitazone
dxdt_PIOG_C  = -CL_piog/Vc_piog * PIOG_C;

// ===== Glucose-Insulin ODEs (Bergman extended) =====
// Plasma glucose (mg/dL)
dxdt_Gp = (EGP_val - Rd_val) / (Vg * BW_t) * 10 - Sg * (Gp - 90) - UGE_loss / (Vg * BW_t);

// Tissue glucose
dxdt_Gt = Sg * (Gp - Gt) - Rd_ins / (Vg * BW_t) * 5;

// Plasma insulin (endogenous, mU/L)
dxdt_Ip = (dIp_sec - CL_ins * Ip / Vc_ins) / Vi / BW_t * 10;

// Remote insulin action compartment
dxdt_X_action = -p2 * X_action + p2 * SI * Ip;

// ===== Glucagon (pg/mL) =====
dxdt_Gc = kout_Gc * (Gc_ss - Gc);

// ===== Active GLP-1 (pmol/L) =====
dxdt_GLP1 = kin_GLP1 + GLP1_Gp * fmax(0.0, Gp - 90) - (kDPP4 * (1 - DPP4_inh) + kout_GLP1) * GLP1;

// ===== β-cell mass (normalized, 0–1) =====
dxdt_beta_mass = (beta_prolif_rate - beta_apop_rate) / (24.0 * 365.25);  // per year

// ===== Insulin Resistance indices =====
dxdt_IR_H = IR_H_drive - IR_H_recov - IR_H_met - dIRH_piog;
dxdt_IR_P = k_IR_FFA * fmax(0.0, FFA - 0.5)
             - k_IR_rec * (IR_P - 1.0) - dIRP_piog
             + 0.001 * fmax(0.0, Gp - 150)/300.0;

// ===== FFA (mmol/L) =====
dxdt_FFA = FFA_release_r - FFA_uptake_r;

// ===== Body weight (kg) =====
dxdt_BW_t = dBW;

// ===== Endpoints =====
dxdt_HbA1c_cmpt = dHbA1c;
dxdt_eGFR_cmpt  = deGFR;
dxdt_UACR_cmpt  = dUACR;

$CAPTURE
Cmet Cempa Csema Cdpp4 Csu Cins Cpiog
EGP_val Rd_val UGE_loss
Gc GLP1
E_empa E_sema_ins DPP4_inh E_su E_piog_IR
beta_prolif_rate beta_apop_rate
IR_H IR_P FFA BW_t
HbA1c_cmpt eGFR_cmpt UACR_cmpt
'

# ============================================================================
# COMPILE MODEL
# ============================================================================
mod <- mcode("T2DM_QSP", t2dm_model_code)

# ============================================================================
# DOSE EVENTS BUILDER
# ============================================================================
build_doses <- function(scenario, start_h = 0, days = 365) {
  events <- list()
  total_h <- days * 24

  if (scenario %in% c("metformin","combo_empa","combo_sema","triple","insulin","dpp4i")) {
    # Metformin 1000 mg BID (F=55% → 550 mg reaches gut cmpt)
    # Represent as mg bioavailable entering MET_GUT
    doses_met <- ev(amt = 1000 * 0.55, cmt = "MET_GUT", ii = 12, addl = 2*days - 1, time = start_h)
    events <- c(events, list(doses_met))
  }

  if (scenario %in% c("combo_empa","triple")) {
    # Empagliflozin 10 mg QD (F=86%)
    doses_empa <- ev(amt = 10 * 0.86 * 1e6 / 450.9, cmt = "EMPA_C",
                     ii = 24, addl = days - 1, time = start_h)
    events <- c(events, list(doses_empa))
  }

  if (scenario %in% c("combo_sema","triple")) {
    # Semaglutide 1 mg SC weekly
    doses_sema <- ev(amt = 1000 / 4113.6 * 1000, cmt = "SEMA_SC",
                     ii = 168, addl = ceiling(days/7) - 1, time = start_h)
    events <- c(events, list(doses_sema))
  }

  if (scenario == "dpp4i") {
    # Sitagliptin 100 mg QD
    doses_dpp4 <- ev(amt = 100 * 0.87 * 1e6 / 407.5, cmt = "DPP4I_C",
                     ii = 24, addl = days - 1, time = start_h)
    events <- c(events, list(doses_dpp4))
  }

  if (scenario == "insulin") {
    # Insulin degludec 20 U QD (1 U = ~0.0347 mg)
    doses_ins <- ev(amt = 20 * 0.0347 * 0.91 * 1000 / 6103, cmt = "INS_SC",
                    ii = 24, addl = days - 1, time = start_h)
    events <- c(events, list(doses_ins))
  }

  if (length(events) == 0) return(ev(amt = 0, cmt = "MET_GUT", time = 0))
  Reduce(c, events)
}

# ============================================================================
# SIMULATION FUNCTION
# ============================================================================
run_scenario <- function(scenario_name, use_flags, days = 365) {
  params_update <- list(
    USE_MET  = use_flags[["met"]],
    USE_EMPA = use_flags[["empa"]],
    USE_SEMA = use_flags[["sema"]],
    USE_DPP4 = use_flags[["dpp4"]],
    USE_SU   = use_flags[["su"]],
    USE_INS  = use_flags[["ins"]],
    USE_PIOG = use_flags[["piog"]]
  )

  dose_ev <- build_doses(scenario_name, days = days)
  sim_times <- seq(0, days * 24, by = 6)  # every 6 hours

  out <- mod %>%
    param(params_update) %>%
    ev(dose_ev) %>%
    mrgsim(end = days * 24, delta = 6) %>%
    as.data.frame()

  out$scenario <- scenario_name
  out
}

# ============================================================================
# SCENARIO DEFINITIONS (7 scenarios)
# ============================================================================
scenarios <- list(
  list(
    name  = "No treatment",
    key   = "no_tx",
    flags = list(met=0, empa=0, sema=0, dpp4=0, su=0, ins=0, piog=0)
  ),
  list(
    name  = "Metformin 1000 BID",
    key   = "metformin",
    flags = list(met=1, empa=0, sema=0, dpp4=0, su=0, ins=0, piog=0)
  ),
  list(
    name  = "Metformin + Empagliflozin",
    key   = "combo_empa",
    flags = list(met=1, empa=1, sema=0, dpp4=0, su=0, ins=0, piog=0)
  ),
  list(
    name  = "Metformin + Semaglutide",
    key   = "combo_sema",
    flags = list(met=1, empa=0, sema=1, dpp4=0, su=0, ins=0, piog=0)
  ),
  list(
    name  = "Triple (Met+Empa+Sema)",
    key   = "triple",
    flags = list(met=1, empa=1, sema=1, dpp4=0, su=0, ins=0, piog=0)
  ),
  list(
    name  = "Metformin + Insulin Degludec",
    key   = "insulin",
    flags = list(met=1, empa=0, sema=0, dpp4=0, su=0, ins=1, piog=0)
  ),
  list(
    name  = "Metformin + Sitagliptin",
    key   = "dpp4i",
    flags = list(met=1, empa=0, sema=0, dpp4=1, su=0, ins=0, piog=0)
  )
)

# ============================================================================
# RUN ALL SCENARIOS
# ============================================================================
cat("Running T2DM QSP simulations...\n")
sim_results <- map_dfr(scenarios, function(s) {
  cat(sprintf("  [%s] %s\n", s$key, s$name))
  run_scenario(s$key, s$flags, days = 365)
})

# Convert time to days
sim_results <- sim_results %>%
  mutate(time_days = time / 24,
         scenario = factor(scenario, levels = map_chr(scenarios, "key"),
                          labels = map_chr(scenarios, "name")))

cat("Simulation complete. N rows:", nrow(sim_results), "\n")

# ============================================================================
# SUMMARY TABLE AT 52 WEEKS
# ============================================================================
summary_52w <- sim_results %>%
  filter(abs(time_days - 365) < 1) %>%
  group_by(scenario) %>%
  slice_tail(n = 1) %>%
  summarise(
    HbA1c_52w     = round(mean(HbA1c_cmpt), 2),
    dHbA1c        = round(mean(HbA1c_cmpt) - 9.0, 2),
    Gp_52w        = round(mean(Gp), 1),
    BW_52w        = round(mean(BW_t), 1),
    dBW           = round(mean(BW_t) - 90.0, 1),
    eGFR_52w      = round(mean(eGFR_cmpt), 1),
    UACR_52w      = round(mean(UACR_cmpt), 0),
    beta_mass_52w = round(mean(beta_mass), 3),
    .groups = "drop"
  )

cat("\n=== 52-Week Summary ===\n")
print(summary_52w, n = 10)

# ============================================================================
# PLOTTING FUNCTIONS
# ============================================================================

# Color palette for 7 scenarios
scen_colors <- c(
  "#e41a1c","#377eb8","#4daf4a","#984ea3",
  "#ff7f00","#a65628","#f781bf"
)

# Plot 1: HbA1c trajectories
p_hba1c <- ggplot(sim_results %>% filter(time_days %% 1 < 0.3),
                  aes(x = time_days, y = HbA1c_cmpt, color = scenario)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 7.0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = scen_colors) +
  labs(title = "HbA1c Trajectory Over 1 Year",
       x = "Time (days)", y = "HbA1c (%)", color = "Treatment") +
  theme_bw() + theme(legend.position = "right")

# Plot 2: Plasma glucose
p_glucose <- ggplot(sim_results %>% filter(time_days %% 1 < 0.3),
                    aes(x = time_days, y = Gp, color = scenario)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 126, linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Plasma Glucose Over Time",
       x = "Time (days)", y = "Plasma Glucose (mg/dL)", color = "Treatment") +
  theme_bw()

# Plot 3: Body weight
p_bw <- ggplot(sim_results %>% filter(time_days %% 1 < 0.3),
               aes(x = time_days, y = BW_t, color = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Body Weight Over Time",
       x = "Time (days)", y = "Body Weight (kg)", color = "Treatment") +
  theme_bw()

# Plot 4: β-cell mass
p_beta <- ggplot(sim_results %>% filter(time_days %% 1 < 0.3),
                 aes(x = time_days, y = beta_mass, color = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = scen_colors) +
  labs(title = "β-cell Mass (Normalized)",
       x = "Time (days)", y = "β-cell Mass (relative)", color = "Treatment") +
  theme_bw()

# Plot 5: eGFR
p_egfr <- ggplot(sim_results %>% filter(time_days %% 1 < 0.3),
                 aes(x = time_days, y = eGFR_cmpt, color = scenario)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 60, linetype = "dashed", color = "orange") +
  scale_color_manual(values = scen_colors) +
  labs(title = "eGFR Over Time (Renal Protection)",
       x = "Time (days)", y = "eGFR (mL/min/1.73m²)", color = "Treatment") +
  theme_bw()

# Plot 6: Urinary Glucose Excretion (SGLT2i effect)
p_uge <- ggplot(sim_results %>%
                  filter(scenario %in% c("Metformin + Empagliflozin","Triple (Met+Empa+Sema)"),
                         time_days %% 1 < 0.3),
                aes(x = time_days, y = UGE_loss * 24, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors[c(3,5)]) +
  labs(title = "Urinary Glucose Excretion (SGLT2i)",
       x = "Time (days)", y = "UGE (mg/kg/day)", color = "Treatment") +
  theme_bw()

# Plot 7: GLP-1 with and without DPP-4 inhibition
p_glp1 <- ggplot(sim_results %>%
                   filter(scenario %in% c("No treatment","Metformin + Sitagliptin",
                                          "Metformin + Semaglutide"),
                          time_days < 7),
                 aes(x = time_days, y = GLP1, color = scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scen_colors[c(1,7,4)]) +
  labs(title = "Active GLP-1 Levels (First Week)",
       x = "Time (days)", y = "Active GLP-1 (pmol/L)", color = "Treatment") +
  theme_bw()

# Summary bar chart
p_summary <- summary_52w %>%
  ggplot(aes(x = reorder(scenario, dHbA1c), y = dHbA1c, fill = scenario)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", dHbA1c)), hjust = 1.1, color = "white") +
  scale_fill_manual(values = scen_colors) +
  coord_flip() +
  labs(title = "ΔHbA1c at 52 Weeks",
       x = NULL, y = "ΔHbA1c (%)") +
  theme_bw()

# Display plots
print(p_hba1c)
print(p_glucose)
print(p_bw)
print(p_beta)
print(p_egfr)
print(summary_52w)

# ============================================================================
# CLINICAL CALIBRATION CHECKPOINTS
# ============================================================================
cat("\n=== Clinical Calibration Validation ===\n")
calib_checks <- tibble::tribble(
  ~Drug,             ~Trial,           ~Endpoint,    ~Observed,  ~Simulated,
  "Metformin",       "UKPDS-34",       "ΔHbA1c(%)",  "-1.4",
    as.character(round(filter(summary_52w, grepl("Metformin 1000", scenario))$dHbA1c[1], 1)),
  "Empagliflozin",   "EMPA-REG",       "ΔHbA1c(%)",  "-0.7",
    as.character(round(filter(summary_52w, grepl("Empa", scenario) & !grepl("Sema", scenario))$dHbA1c[1], 1)),
  "Semaglutide",     "SUSTAIN-6",      "ΔHbA1c(%)",  "-1.4",
    as.character(round(filter(summary_52w, grepl("Sema", scenario) & !grepl("Empa", scenario))$dHbA1c[1], 1)),
  "Semaglutide",     "SUSTAIN-6",      "ΔBW(kg)",    "-4.5",
    as.character(round(filter(summary_52w, grepl("Sema", scenario) & !grepl("Empa", scenario))$dBW[1], 1)),
  "Sitagliptin",     "TECOS",          "ΔHbA1c(%)",  "-0.7",
    as.character(round(filter(summary_52w, grepl("Sitagliptin", scenario))$dHbA1c[1], 1))
)
print(calib_checks)

cat("\nModel ready. Use Shiny app for interactive exploration.\n")
