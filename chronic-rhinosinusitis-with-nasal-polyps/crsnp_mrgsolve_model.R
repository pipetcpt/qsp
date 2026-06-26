## =============================================================================
## CRSwNP QSP Model — mrgsolve Implementation
## Chronic Rhinosinusitis with Nasal Polyps
## Type 2 Inflammation / Eosinophilic Endotype
##
## Key clinical trials calibrated:
##   - LIBERTY NP SINUS-24 & SINUS-52 (dupilumab)
##   - SYNAPSE (mepolizumab 100 mg SC q4w)
##   - OSTRO (benralizumab 30 mg SC q4w/q8w)
##   - POLYP 1 & POLYP 2 (omalizumab)
##   - WAYPOINT (tezepelumab 210 mg SC q4w)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

## ---------------------------------------------------------------------------
## Model code block
## ---------------------------------------------------------------------------
code <- '
$PROB
CRSwNP QSP Model | 5 biologics + INCS | 22-ODE system
Dupilumab (anti-IL-4Ra), Mepolizumab (anti-IL-5), Benralizumab (anti-IL-5Ra),
Omalizumab (anti-IgE), Tezepelumab (anti-TSLP), INCS, Montelukast

$PARAM @annotated
// ---- Drug selection & dosing ----
DRUG       : 0     : Drug (0=none, 1=Dupilumab, 2=Mepolizumab, 3=Benralizumab, 4=Omalizumab, 5=Tezepelumab)
USE_INCS   : 1     : INCS use (0=no, 1=yes)
USE_MLK    : 0     : Montelukast (0=no, 1=yes)

// ---- Dupilumab PK (SINUS-52 popPK: 2-cmpt SC) ----
F_DUP      : 0.64  : Dupilumab SC bioavailability
KA_DUP     : 0.14  : Absorption rate constant (1/day)
V1_DUP     : 3.10  : Central volume (L)
V2_DUP     : 3.40  : Peripheral volume (L)
CL_DUP     : 0.50  : Clearance (L/day), t1/2~21d
Q_DUP      : 0.80  : Intercompartmental CL (L/day)

// ---- Mepolizumab PK (SYNAPSE popPK) ----
F_MEP      : 0.80  : Mepolizumab SC bioavailability
KA_MEP     : 0.12  : Absorption rate constant (1/day)
V1_MEP     : 3.60  : Central volume (L)
CL_MEP     : 0.22  : Clearance (L/day), t1/2~20d

// ---- Benralizumab PK (OSTRO popPK) ----
F_BEN      : 0.56  : Benralizumab SC bioavailability
KA_BEN     : 0.10  : Absorption rate constant (1/day)
V1_BEN     : 3.10  : Central volume (L)
CL_BEN     : 0.35  : Clearance (L/day), t1/2~15d

// ---- Omalizumab PK (POLYP trials popPK) ----
F_OMA      : 0.62  : Omalizumab SC bioavailability
KA_OMA     : 0.10  : Absorption rate constant (1/day)
V1_OMA     : 4.00  : Central volume (L)
CL_OMA     : 0.19  : Clearance (L/day), t1/2~26d

// ---- Tezepelumab PK (WAYPOINT popPK) ----
F_TEZ      : 0.77  : Tezepelumab SC bioavailability
KA_TEZ     : 0.09  : Absorption rate constant (1/day)
V1_TEZ     : 4.00  : Central volume (L)
CL_TEZ     : 0.18  : Clearance (L/day), t1/2~26d

// ---- PD parameters (drug effects) ----
IC50_DUP   : 1.20  : Dupilumab IC50 on IL-4Rα signal (μg/mL)
IC50_MEP   : 0.50  : Mepolizumab IC50 on free IL-5 (μg/mL)
IC50_BEN   : 0.30  : Benralizumab IC50 on IL-5Rα / ADCC Eos depletion (μg/mL)
IC50_OMA   : 2.50  : Omalizumab IC50 on free IgE (μg/mL)
IC50_TEZ   : 1.50  : Tezepelumab IC50 on free TSLP (μg/mL)
EMAX_INCS  : 0.50  : Maximum INCS Emax on cytokine production
EMAX_MLK   : 0.35  : Maximum montelukast Emax on CysLT signaling

// ---- Disease biology parameters ----
// Epithelial barrier
K_EPI_DMG  : 0.08  : Rate of epithelial damage per EosTissue (1/day/AU)
K_EPI_REP  : 0.05  : Spontaneous epithelial repair rate (1/day)
EPI_SS     : 0.40  : Disease steady-state epithelial integrity

// TSLP dynamics
K_TSLP_P   : 0.50  : TSLP production rate from damaged epithelium
K_TSLP_D   : 0.50  : TSLP first-order degradation (1/day)

// ILC2 dynamics
K_ILC2_P   : 0.30  : ILC2 activation per TSLP (AU/day)
K_ILC2_D   : 0.25  : ILC2 resolution rate (1/day)

// Th2 dynamics
K_TH2_P    : 0.20  : Th2 polarization rate per IL4 signal
K_TH2_D    : 0.18  : Th2 resolution rate (1/day)

// Cytokine dynamics
K_IL4_P    : 0.40  : IL-4 production per ILC2+Th2 (AU/day)
K_IL4_D    : 0.80  : IL-4 degradation rate (1/day), t1/2~0.9d
K_IL5_P    : 0.30  : IL-5 production rate (AU/day)
K_IL5_D    : 0.60  : IL-5 degradation rate (1/day), t1/2~1.2d
K_IL13_P   : 0.30  : IL-13 production rate (AU/day)
K_IL13_D   : 0.50  : IL-13 degradation rate (1/day)

// IgE dynamics
K_IGE_P    : 0.06  : IgE production per IL4+IL13 (kU/L/day per AU)
K_IGE_D    : 0.030 : IgE degradation (1/day), t1/2~23d
IGE_BAS    : 300   : Baseline total IgE (kU/L)

// Eosinophil dynamics
K_EOSB_P   : 0.50  : Blood eosinophil production per IL5 (AU/day)
K_EOSB_D   : 0.10  : Blood eosinophil clearance (1/day), t1/2~7d
EOSB_BAS   : 500   : Baseline blood eosinophils (cells/μL)
K_EOST_R   : 0.04  : Tissue eosinophil recruitment per EOSB*IL5
K_EOST_D   : 0.20  : Tissue eosinophil clearance (1/day)
EOST_BAS   : 80    : Baseline tissue eosinophil AU

// Goblet cells
K_GOBC_P   : 0.20  : Goblet cell expansion per IL13
K_GOBC_D   : 0.08  : Goblet cell normalization (1/day)

// TGF-beta
K_TGFB_P   : 0.15  : TGF-β production per EosTissue
K_TGFB_D   : 0.25  : TGF-β degradation (1/day)

// VEGF
K_VEGF_P   : 0.12  : VEGF production per EosTissue+TGFb
K_VEGF_D   : 0.20  : VEGF degradation (1/day)

// Nasal Polyp Score
K_NPS_G    : 0.012 : NPS growth rate per EOST*VEGF (1/day)
K_NPS_R    : 0.015 : NPS spontaneous resolution (1/day)
NPS_MAX    : 8.0   : Maximum NPS
NPS_BAS    : 5.5   : Disease-state baseline NPS (0-8)

// Symptom parameters
NPS_SNOT   : 3.5   : NPS-to-SNOT22 coefficient
EOST_OBS   : 2.5   : EosTissue to obstruction coefficient
TGFB_FIB   : 1.8   : TGF-β to fibrosis coefficient (LM CT score proxy)

$CMT @annotated
// --- Dupilumab PK ---
D_SC   : Dupilumab SC depot (mg)
D_C1   : Dupilumab central plasma (mg/L = μg/mL)
D_P1   : Dupilumab peripheral (mg)

// --- Mepolizumab PK ---
M_SC   : Mepolizumab SC depot (mg)
M_C1   : Mepolizumab central plasma (μg/mL)

// --- Benralizumab PK ---
B_SC   : Benralizumab SC depot (mg)
B_C1   : Benralizumab central plasma (μg/mL)

// --- Omalizumab PK ---
O_SC   : Omalizumab SC depot (mg)
O_C1   : Omalizumab central plasma (μg/mL)

// --- Tezepelumab PK ---
T_SC   : Tezepelumab SC depot (mg)
T_C1   : Tezepelumab central plasma (μg/mL)

// --- Disease state variables ---
EPI    : Epithelial barrier integrity (0-1)
TSLP   : TSLP (AU)
ILC2   : ILC2 cell density (AU)
TH2    : Th2 cell density (AU)
IL4    : IL-4 (AU)
IL5    : IL-5 (AU)
IL13   : IL-13 (AU)
IGE    : Total serum IgE (kU/L)
EOSB   : Blood eosinophils (cells/μL)
EOST   : Tissue eosinophils (AU)
GOBC   : Goblet cell density (AU)
TGFB   : TGF-β (AU)
VEGF   : VEGF (AU)
NPS    : Nasal polyp score (0-8)

$MAIN
// Initial conditions at disease steady state
EPI_0   = EPI_SS;
TSLP_0  = K_TSLP_P * (1 - EPI_SS) / K_TSLP_D;
ILC2_0  = K_ILC2_P * TSLP_0 / K_ILC2_D;
double th2ss = 0.8;
TH2_0   = th2ss;
double src0 = ILC2_0 + TH2_0;
IL4_0   = K_IL4_P * src0 / K_IL4_D;
IL5_0   = K_IL5_P * src0 / K_IL5_D;
IL13_0  = K_IL13_P * src0 / K_IL13_D;
IGE_0   = IGE_BAS;
EOSB_0  = EOSB_BAS;
EOST_0  = EOST_BAS;
GOBC_0  = K_GOBC_P * IL13_0 / K_GOBC_D;
TGFB_0  = K_TGFB_P * EOST_0 / K_TGFB_D;
VEGF_0  = K_VEGF_P * (EOST_0 + TGFB_0) / K_VEGF_D;
NPS_0   = NPS_BAS;

$ODE
// ============ Drug PK ODEs ============

// Dupilumab
dxdt_D_SC = -KA_DUP * D_SC;
dxdt_D_C1 =  KA_DUP * F_DUP * D_SC / V1_DUP
             - (CL_DUP/V1_DUP) * D_C1
             - (Q_DUP/V1_DUP)  * D_C1
             + (Q_DUP/V2_DUP)  * D_P1;
dxdt_D_P1 = (Q_DUP/V1_DUP) * D_C1 - (Q_DUP/V2_DUP) * D_P1;

// Mepolizumab
dxdt_M_SC = -KA_MEP * M_SC;
dxdt_M_C1 =  KA_MEP * F_MEP * M_SC / V1_MEP - (CL_MEP/V1_MEP) * M_C1;

// Benralizumab
dxdt_B_SC = -KA_BEN * B_SC;
dxdt_B_C1 =  KA_BEN * F_BEN * B_SC / V1_BEN - (CL_BEN/V1_BEN) * B_C1;

// Omalizumab
dxdt_O_SC = -KA_OMA * O_SC;
dxdt_O_C1 =  KA_OMA * F_OMA * O_SC / V1_OMA - (CL_OMA/V1_OMA) * O_C1;

// Tezepelumab
dxdt_T_SC = -KA_TEZ * T_SC;
dxdt_T_C1 =  KA_TEZ * F_TEZ * T_SC / V1_TEZ - (CL_TEZ/V1_TEZ) * T_C1;

// ============ Drug effect functions (Emax) ============
// Dupilumab: blocks IL-4Rα → reduces IL-4 AND IL-13 downstream signaling
double E_DUP = (DRUG==1) ? D_C1 / (IC50_DUP + D_C1) : 0;

// Mepolizumab: neutralizes IL-5
double E_MEP = (DRUG==2) ? M_C1 / (IC50_MEP + M_C1) : 0;

// Benralizumab: blocks IL-5Rα + ADCC → near-complete Eos depletion
double E_BEN_B = (DRUG==3) ? B_C1 / (IC50_BEN + B_C1) : 0;  // blood Eos
double E_BEN_T = (DRUG==3) ? B_C1 / (IC50_BEN * 0.5 + B_C1) : 0;  // tissue (ADCC enhanced)

// Omalizumab: neutralizes free IgE
double E_OMA = (DRUG==4) ? O_C1 / (IC50_OMA + O_C1) : 0;

// Tezepelumab: neutralizes TSLP
double E_TEZ = (DRUG==5) ? T_C1 / (IC50_TEZ + T_C1) : 0;

// INCS: broad suppression via NF-κB (Emax model, constant steady-state effect)
double E_INCS = USE_INCS * EMAX_INCS;

// Montelukast: blocks CysLT1R (symptom relief, modest effect on inflammation)
double E_MLK = USE_MLK * EMAX_MLK;

// ============ Disease Biology ODEs ============
// Bound concentrations
double EPI_c  = EPI;
double TSLP_c = TSLP;
double ILC2_c = ILC2;
double TH2_c  = TH2;
double IL4_c  = IL4;
double IL5_c  = IL5;
double IL13_c = IL13;
double EOSB_c = EOSB;
double EOST_c = EOST;
double TGFB_c = TGFB;

// --- Epithelial barrier ---
// Damage by eosinophil toxins, repair attempted; tezepelumab/dupilumab help repair
double epi_damage  = K_EPI_DMG * EOST_c * (1.0 - EPI_c);
double epi_repair  = K_EPI_REP * (1.0 - EPI_c) * (1 + E_DUP * 0.5 + E_TEZ * 0.5);
dxdt_EPI = -epi_damage + epi_repair;

// --- TSLP ---
double tslp_prod = K_TSLP_P * (1.0 - EPI_c);
double tslp_deg  = K_TSLP_D * TSLP_c * (1 + E_TEZ);  // tezepelumab neutralizes
dxdt_TSLP = tslp_prod - tslp_deg;

// --- ILC2 ---
double ilc2_act  = K_ILC2_P * TSLP_c;
double ilc2_inh  = K_ILC2_D * ILC2_c;
dxdt_ILC2 = ilc2_act - ilc2_inh;

// --- Th2 ---
// IL-4 positive feedback drives Th2 polarization; dupilumab reduces IL-4 signal
double il4_signal = IL4_c * (1 - E_DUP);
double th2_prod   = K_TH2_P * il4_signal + 0.05 * ILC2_c;  // ILC2 also supports
double th2_res    = K_TH2_D * TH2_c;
dxdt_TH2 = th2_prod - th2_res;

// --- Cytokine source (ILC2 + Th2) ---
double src = ILC2_c + TH2_c;

// --- IL-4 ---
// Production by ILC2 + Th2; suppressed by INCS and dupilumab (signal blocked)
double il4_prod = K_IL4_P * src * (1 - E_INCS);
double il4_deg  = K_IL4_D * IL4_c;
dxdt_IL4 = il4_prod - il4_deg;

// --- IL-5 ---
// Production; blocked by mepolizumab (neutralizes IL-5 itself)
double il5_prod = K_IL5_P * src * (1 - E_INCS);
double il5_deg  = K_IL5_D * IL5_c * (1 + E_MEP * 3.0);  // mepolizumab enhances clearance
dxdt_IL5 = il5_prod - il5_deg;

// --- IL-13 ---
// Production; dupilumab blocks IL-4Rα → downstream effects blocked
// INCS reduces production
double il13_prod = K_IL13_P * src * (1 - E_INCS);
double il13_deg  = K_IL13_D * IL13_c;
dxdt_IL13 = il13_prod - il13_deg;

// --- Total serum IgE ---
// Production driven by IL-4 + IL-13 class switching; omalizumab neutralizes
double ige_prod  = K_IGE_P * (IL4_c + IL13_c) * (1 - E_DUP * 0.8);
double ige_deg   = K_IGE_D * IGE;
// Omalizumab forms immune complexes → enhanced removal
double ige_oma   = E_OMA * K_IGE_D * 3.0 * IGE;
dxdt_IGE = ige_prod - ige_deg - ige_oma;

// --- Blood eosinophils (cells/μL) ---
// Produced by IL-5; benralizumab ADCC depletes; mepolizumab reduces IL-5 signal
double eosb_in  = K_EOSB_P * IL5_c * (1 - E_MEP * 0.9) * (EOSB_BAS / 100.0);
double eosb_out = K_EOSB_D * EOSB_c;
double eosb_ben = E_BEN_B * K_EOSB_D * 15.0 * EOSB_c;  // rapid ADCC depletion
dxdt_EOSB = eosb_in - eosb_out - eosb_ben;

// --- Tissue eosinophils (AU) ---
// Recruited from blood via IL-5 / eotaxins; INCS reduces recruitment
// Benralizumab ADCC depletes tissue eosinophils
double eost_in  = K_EOST_R * EOSB_c * IL5_c * (1 - E_INCS * 0.6);
double eost_out = K_EOST_D * EOST_c;
double eost_ben = E_BEN_T * K_EOST_D * 8.0 * EOST_c;
// Dupilumab reduces tissue eos indirectly (via IL-4/IL-13 → eotaxin-3 ↓)
double eost_dup = E_DUP * K_EOST_D * 0.8 * EOST_c;
dxdt_EOST = eost_in - eost_out - eost_ben - eost_dup;

// --- Goblet cell density (AU) ---
// IL-13 drives goblet cell hyperplasia; dupilumab (via IL-13Rα1/IL-4Rα) reduces
double gobc_prod = K_GOBC_P * IL13_c * (1 - E_DUP * 0.85);
double gobc_res  = K_GOBC_D * GOBC;
dxdt_GOBC = gobc_prod - gobc_res;

// --- TGF-β (AU) ---
// Produced by tissue eosinophils; contributes to fibrosis
double tgfb_prod = K_TGFB_P * EOST_c;
double tgfb_deg  = K_TGFB_D * TGFB;
dxdt_TGFB = tgfb_prod - tgfb_deg;

// --- VEGF (AU) ---
// Produced by mast cells + eosinophils; INCS reduces
double vegf_prod = K_VEGF_P * (EOST_c + TGFB) * (1 - E_INCS * 0.5);
double vegf_deg  = K_VEGF_D * VEGF;
dxdt_VEGF = vegf_prod - vegf_deg;

// --- Nasal Polyp Score (NPS, 0-8) ---
// Growth driven by tissue eosinophilia + VEGF + TGF-β (fibrosis)
// Resolution enhanced by treatment; FESS modeled via event
double nps_grow  = K_NPS_G * EOST_c * VEGF * (NPS_MAX - NPS) / NPS_MAX;
double nps_res   = K_NPS_R * NPS;
// Treatment-enhanced resolution
double nps_treat = E_DUP * 0.040 * NPS + E_MEP * 0.015 * NPS
                 + E_BEN_T * 0.015 * NPS + E_OMA * 0.018 * NPS
                 + E_TEZ * 0.022 * NPS + E_INCS * 0.008 * NPS;
dxdt_NPS = nps_grow - nps_res - nps_treat;

$TABLE
// ---- Derived outcomes ----
double Cp_DUP = D_C1;   // μg/mL
double Cp_MEP = M_C1;
double Cp_BEN = B_C1;
double Cp_OMA = O_C1;
double Cp_TEZ = T_C1;

// Nasal Obstruction VAS (0-10): driven by EOST and NPS
double OBS_VAS = fmin(10.0, fmax(0.0, 3.5 + EOST_OBS * EOST / EOST_BAS * 2.0
                   + 0.8 * NPS + E_MLK * (-2.0)));

// Olfactory score: inversely related to NPS
double OLFACT = fmax(0.0, 10.0 - 1.2 * NPS);

// SNOT-22 (0-110): composite of NPS + obstruction + olfactory loss
double SNOT22 = fmin(110.0, fmax(0.0,
    NPS * NPS_SNOT + OBS_VAS * 3.0 + (10.0 - OLFACT) * 2.0 + GOBC * 4.0));

// Lund-Mackay CT score proxy (0-24)
double LM_CT = fmin(24.0, fmax(0.0, NPS * 2.0 + TGFB * TGFB_FIB));

// Serum biomarkers
double BLD_EOS = EOSB;          // cells/μL
double SERUM_IGE = IGE;         // kU/L
double BLOOD_ECP = EOST * 12.0; // μg/L proxy

// FeNO (ppb) — paradoxically low in CRSwNP (mucosal consumption); rises with treatment
double FeNO = fmax(5.0, 25.0 - NPS * 2.0 + E_DUP * 15.0);

// NPS components
double NPS_L = NPS / 2.0;
double NPS_R = NPS / 2.0;

// Mucosal periostin
double PERIOSTIN = IL13 * 8.0 + EOST * 1.5;  // ng/mL proxy

$CAPTURE
Cp_DUP Cp_MEP Cp_BEN Cp_OMA Cp_TEZ
EPI TSLP ILC2 TH2 IL4 IL5 IL13 IGE EOSB EOST GOBC TGFB VEGF NPS
OBS_VAS OLFACT SNOT22 LM_CT BLD_EOS SERUM_IGE BLOOD_ECP FeNO
NPS_L NPS_R PERIOSTIN
'

## Load model
mod <- mread("crsnp", tempdir(), code, quiet = TRUE)

## ---------------------------------------------------------------------------
## Dosing events helper
## ---------------------------------------------------------------------------
# Study: 52-week treatment, weekly obs
obs_times <- seq(0, 365, by = 7)

make_events <- function(drug, dose_mg, interval_days, n_doses,
                        use_incs = 1, use_mlk = 0, n_benra_q4 = 3) {
  if (drug == 0) {
    return(ev(time = 0, amt = 0, cmt = 1))  # placeholder
  }
  cmt_map <- c("1" = 1, "2" = 4, "3" = 6, "4" = 8, "5" = 10)  # SC compartments
  cmt_id  <- cmt_map[as.character(drug)]

  if (drug == 3) {  # Benralizumab: q4w x 3 then q8w
    doses_q4 <- ev(time = seq(0, (n_benra_q4 - 1) * 28, by = 28),
                   amt = dose_mg, cmt = cmt_id)
    start_q8  <- n_benra_q4 * 28
    n_q8      <- ceiling((n_doses - n_benra_q4) / 1)
    doses_q8  <- ev(time = seq(start_q8, start_q8 + (n_q8 - 1) * 56, by = 56),
                    amt = dose_mg, cmt = cmt_id)
    return(c(doses_q4, doses_q8))
  }

  ev(time = seq(0, (n_doses - 1) * interval_days, by = interval_days),
     amt  = dose_mg, cmt = cmt_id)
}

## ---------------------------------------------------------------------------
## Simulation scenarios
## ---------------------------------------------------------------------------
scenarios <- list(
  list(name = "No Treatment",
       drug = 0, dose = 0, interval = 14, ndoses = 26,
       use_incs = 0, use_mlk = 0),
  list(name = "INCS Only",
       drug = 0, dose = 0, interval = 14, ndoses = 26,
       use_incs = 1, use_mlk = 0),
  list(name = "Dupilumab 300mg q2w + INCS",
       drug = 1, dose = 300, interval = 14, ndoses = 26,
       use_incs = 1, use_mlk = 0),
  list(name = "Mepolizumab 100mg q4w + INCS",
       drug = 2, dose = 100, interval = 28, ndoses = 13,
       use_incs = 1, use_mlk = 0),
  list(name = "Benralizumab 30mg q4-8w + INCS",
       drug = 3, dose = 30,  interval = 28, ndoses = 10,
       use_incs = 1, use_mlk = 0),
  list(name = "Omalizumab 300mg q4w + INCS",
       drug = 4, dose = 300, interval = 28, ndoses = 13,
       use_incs = 1, use_mlk = 0),
  list(name = "Tezepelumab 210mg q4w + INCS",
       drug = 5, dose = 210, interval = 28, ndoses = 13,
       use_incs = 1, use_mlk = 0)
)

run_scenario <- function(sc) {
  params_update <- list(
    DRUG     = sc$drug,
    USE_INCS = sc$use_incs,
    USE_MLK  = sc$use_mlk
  )
  mod_sc <- mod %>% param(params_update)
  dosing <- make_events(sc$drug, sc$dose, sc$interval, sc$ndoses,
                        sc$use_incs, sc$use_mlk)
  out <- mrgsim(mod_sc, events = dosing, obsonly = TRUE,
                tgrid = obs_times, digits = 4)
  as.data.frame(out) %>%
    mutate(Scenario = sc$name, Week = time / 7)
}

set.seed(42)
results_list <- lapply(scenarios, run_scenario)
results_all  <- bind_rows(results_list) %>%
  mutate(Scenario = factor(Scenario, levels = sapply(scenarios, `[[`, "name")))

## ---------------------------------------------------------------------------
## Clinical Endpoint Summary at Week 24 & 52
## ---------------------------------------------------------------------------
summary_table <- results_all %>%
  filter(Week %in% c(0, 24, 52)) %>%
  group_by(Scenario, Week) %>%
  summarise(
    NPS       = round(mean(NPS), 2),
    SNOT22    = round(mean(SNOT22), 1),
    OBS_VAS   = round(mean(OBS_VAS), 2),
    BLD_EOS   = round(mean(BLD_EOS), 0),
    SERUM_IGE = round(mean(SERUM_IGE), 0),
    OLFACT    = round(mean(OLFACT), 2),
    .groups   = "drop"
  )

cat("=== Clinical Endpoint Summary ===\n")
cat("Endpoint: Nasal Polyp Score (NPS, 0-8) — lower is better\n\n")
summary_table %>%
  select(Scenario, Week, NPS, SNOT22, OBS_VAS, BLD_EOS, SERUM_IGE) %>%
  arrange(Scenario, Week) %>%
  print(n = Inf)

## ---------------------------------------------------------------------------
## Figure 1: NPS Trajectories
## ---------------------------------------------------------------------------
p_nps <- ggplot(results_all,
                aes(x = Week, y = NPS, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_y_continuous(limits = c(0, 8), breaks = 0:8) +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Nasal Polyp Score (NPS) Over 52 Weeks",
       subtitle = "CRSwNP QSP Model — Primary Efficacy Endpoint",
       x = "Week", y = "NPS (0–8)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8)) +
  guides(color = guide_legend(ncol = 2), linetype = guide_legend(ncol = 2))

## ---------------------------------------------------------------------------
## Figure 2: Blood Eosinophils
## ---------------------------------------------------------------------------
p_eos <- ggplot(results_all,
                aes(x = Week, y = BLD_EOS, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  labs(title = "Blood Eosinophil Count Over 52 Weeks",
       x = "Week", y = "Eosinophils (cells/μL)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

## ---------------------------------------------------------------------------
## Figure 3: SNOT-22 PRO
## ---------------------------------------------------------------------------
p_snot <- ggplot(results_all,
                 aes(x = Week, y = SNOT22, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  labs(title = "SNOT-22 Score Over 52 Weeks (PRO)",
       subtitle = "MCID = 8.9 points",
       x = "Week", y = "SNOT-22 (0–110)") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "gray50") +
  annotate("text", x = 2, y = 21, label = "Well-controlled threshold", hjust = 0, size = 3) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

## ---------------------------------------------------------------------------
## Figure 4: Nasal Obstruction VAS
## ---------------------------------------------------------------------------
p_obs <- ggplot(results_all,
                aes(x = Week, y = OBS_VAS, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  scale_y_continuous(limits = c(0, 10)) +
  labs(title = "Nasal Obstruction VAS",
       x = "Week", y = "Obstruction VAS (0–10)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

## ---------------------------------------------------------------------------
## Figure 5: PK Profiles for Dupilumab
## ---------------------------------------------------------------------------
dup_data <- results_all %>% filter(Scenario == "Dupilumab 300mg q2w + INCS")
p_pk <- ggplot(dup_data, aes(x = Week, y = Cp_DUP)) +
  geom_line(color = "#1B7837", linewidth = 1.2) +
  labs(title = "Dupilumab PK — Cp (300 mg SC q2w)",
       subtitle = "2-compartment model, F=0.64, CL=0.50 L/d, t½≈21 days",
       x = "Week", y = "Dupilumab Cp (μg/mL)") +
  theme_bw(base_size = 12)

## ---------------------------------------------------------------------------
## Figure 6: Biomarker–response (serum IgE)
## ---------------------------------------------------------------------------
p_ige <- ggplot(results_all,
                aes(x = Week, y = SERUM_IGE, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  labs(title = "Total Serum IgE Over 52 Weeks",
       x = "Week", y = "Total IgE (kU/L)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

## ---------------------------------------------------------------------------
## Combine plots
## ---------------------------------------------------------------------------
combined_plot <- (p_nps | p_eos) / (p_snot | p_obs) / (p_pk | p_ige) +
  plot_annotation(
    title    = "CRSwNP QSP Model — Comparative Efficacy of Biologics",
    subtitle = "Simulated 52-week treatment trajectories (INCS + biologic vs. INCS alone vs. No Tx)",
    theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                     plot.subtitle = element_text(size = 11))
  )

## ---------------------------------------------------------------------------
## Treatment Discontinuation Scenario (relapse after stopping at Week 52)
## ---------------------------------------------------------------------------
relapse_scenarios <- list(
  list(name = "Dupilumab — then Stop",
       drug = 1, dose = 300, interval = 14, ndoses = 26,
       use_incs = 1, use_mlk = 0),
  list(name = "Mepolizumab — then Stop",
       drug = 2, dose = 100, interval = 28, ndoses = 13,
       use_incs = 1, use_mlk = 0)
)

obs_long <- c(obs_times, seq(366, 730, by = 7))

run_with_stop <- function(sc, stop_week = 52) {
  params_update <- list(
    DRUG     = sc$drug,
    USE_INCS = sc$use_incs,
    USE_MLK  = sc$use_mlk
  )
  mod_sc <- mod %>% param(params_update)
  dosing <- make_events(sc$drug, sc$dose, sc$interval, sc$ndoses,
                        sc$use_incs, sc$use_mlk)

  # Phase 1: with drug
  p1 <- mrgsim(mod_sc, events = dosing, obsonly = TRUE,
               tgrid = obs_times, digits = 4) %>% as.data.frame()

  # Phase 2: stop drug, INCS continues
  init_vals <- tail(p1, 1)
  init_list <- as.list(init_vals[, c("EPI","TSLP","ILC2","TH2","IL4","IL5",
                                      "IL13","IGE","EOSB","EOST","GOBC",
                                      "TGFB","VEGF","NPS")])
  mod_stop <- mod %>%
    param(DRUG = 0, USE_INCS = 1) %>%
    init(init_list)

  p2 <- mrgsim(mod_stop, obsonly = TRUE,
               tgrid = seq(0, 365, by = 7), digits = 4) %>%
    as.data.frame() %>%
    mutate(time = time + 365)

  bind_rows(p1, p2) %>%
    mutate(Scenario = sc$name, Week = time / 7, Phase = ifelse(time <= 365, "Treatment", "Post-Stop"))
}

relapse_results <- lapply(relapse_scenarios, run_with_stop) %>% bind_rows()

p_relapse <- ggplot(relapse_results,
                    aes(x = Week, y = NPS, color = Scenario,
                        linetype = Phase)) +
  geom_line(linewidth = 1.0) +
  geom_vline(xintercept = 52, linetype = "dotted", color = "red", linewidth = 1) +
  annotate("text", x = 53, y = 7, label = "Drug stopped", hjust = 0, color = "red", size = 3.5) +
  scale_color_manual(values = c("#1B7837", "#8B0000")) +
  scale_x_continuous(breaks = seq(0, 104, 8)) +
  scale_y_continuous(limits = c(0, 8), breaks = 0:8) +
  labs(title = "NPS Trajectory: Treatment & Relapse After Drug Discontinuation",
       subtitle = "INCS maintained throughout; biologic stopped at Week 52",
       x = "Week", y = "NPS (0–8)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

## ---------------------------------------------------------------------------
## Print summary
## ---------------------------------------------------------------------------
cat("\n\n=== NPS change from baseline at Week 24 ===\n")
nps_baseline <- results_all %>% filter(Week == 0) %>%
  group_by(Scenario) %>% summarise(NPS_0 = mean(NPS))
nps_wk24 <- results_all %>% filter(abs(Week - 24) < 0.5) %>%
  group_by(Scenario) %>% summarise(NPS_24 = mean(NPS))
nps_change <- left_join(nps_baseline, nps_wk24, by = "Scenario") %>%
  mutate(Delta_NPS = round(NPS_24 - NPS_0, 2),
         Pct_Change = round(100 * (NPS_24 - NPS_0) / NPS_0, 1)) %>%
  arrange(Delta_NPS)
print(nps_change)

cat("\n\n=== NPS change from baseline at Week 52 ===\n")
nps_wk52 <- results_all %>% filter(abs(Week - 52) < 0.5) %>%
  group_by(Scenario) %>% summarise(NPS_52 = mean(NPS))
nps_change52 <- left_join(nps_baseline, nps_wk52, by = "Scenario") %>%
  mutate(Delta_NPS = round(NPS_52 - NPS_0, 2),
         Pct_Change = round(100 * (NPS_52 - NPS_0) / NPS_0, 1)) %>%
  arrange(Delta_NPS)
print(nps_change52)

## ---------------------------------------------------------------------------
## Save outputs
## ---------------------------------------------------------------------------
ggsave("crsnp_qsp_results.png", plot = combined_plot,
       width = 16, height = 14, dpi = 150)
ggsave("crsnp_relapse_analysis.png", plot = p_relapse,
       width = 12, height = 6, dpi = 150)
cat("\nPlots saved: crsnp_qsp_results.png, crsnp_relapse_analysis.png\n")
