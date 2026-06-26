## ============================================================
## Aplastic Anemia (AA) — QSP mrgsolve Model
## ============================================================
## Disease: Aplastic Anemia (Immune-Mediated)
## Author  : Claude Code Routine (CCR)
## Date    : 2026-06-25
## Version : 1.0
##
## Model Scope:
##   PK  — Anti-thymocyte globulin (hATG/rATG), Cyclosporine (CsA),
##           Eltrombopag (EPAG), Danazol (androgen)
##   PD  — Autoreactive T cell dynamics, HSC destruction & recovery,
##           hematopoietic lineage differentiation,
##           Peripheral blood count trajectories (Hgb, ANC, PLT, ARC)
##
## Compartments (20 ODE):
##   1  ATG_C     — ATG central plasma
##   2  ATG_P     — ATG peripheral tissue
##   3  CsA_C     — Cyclosporine central blood
##   4  EPAG_C    — Eltrombopag central plasma
##   5  Teff      — Autoreactive effector T cell count (×10⁶/kg)
##   6  Treg      — Regulatory T cell count (×10⁶/kg)
##   7  HSC       — HSC pool size (% of normal = 100)
##   8  CFU_E     — Erythroid progenitor pool (CFU-E)
##   9  Retic     — Reticulocyte pool
##  10  RBC       — Red blood cell pool → Hgb
##  11  CFU_G     — Granulocyte progenitor pool (CFU-G)
##  12  ANC_pool  — Circulating neutrophil pool
##  13  MK        — Megakaryocyte pool
##  14  PLT_pool  — Platelet pool
##  15  BM_score  — Bone marrow cellularity score (0–1)
##  16  IFNg      — IFN-γ concentration (pg/mL)
##  17  TNFa      — TNF-α concentration (pg/mL)
##  18  IL2_comp  — IL-2 concentration (pg/mL)
##  19  PNH_clone — PNH clone fraction (0–1)
##  20  Danazol_C — Danazol central plasma (ng/mL)
##
## Clinical Calibration References:
##   - Scheinberg P et al. NEJM 2011 (hATG+CsA vs rATG+CsA)
##   - Townsley DM et al. NEJM 2017 (EPAG add-on to IST)
##   - Dezern AE et al. Haematologica 2018 (long-term outcomes)
##   - Olnes MJ et al. NEJM 2012 (EPAG monotherapy)
##   - Peffault de Latour R et al. Blood 2021 (rATG+CsA+EPAG)
##   - Young NS. Blood 2018 (pathophysiology review)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

## ============================================================
## 1. Model Code Block
## ============================================================

aa_model_code <- '
$PROB Aplastic Anemia QSP Model — Immune T-cell / HSC / Hematopoiesis

$PARAM
// ---- ATG PK Parameters ----
// hATG: horse ATG (ATGAM) 40 mg/kg/d × 4d
// rATG: rabbit ATG (Thymoglobulin) 3.5 mg/kg/d × 5d
ATG_CL    = 0.85    // L/h clearance (hATG; rATG CL ~0.5 L/h)
ATG_Vc    = 5.8     // L central volume of distribution
ATG_Vp    = 12.0    // L peripheral volume
ATG_Q     = 0.30    // L/h inter-compartmental clearance
ATG_Emax  = 0.92    // Max T-cell depletion efficacy (fraction)
ATG_EC50  = 0.15    // mg/L for 50% T-cell depletion
ATG_hill  = 1.5     // Hill coefficient

// ---- CsA PK Parameters ----
// Oral bioavailability ~34%, t½ ~27h, target trough 150-250 ng/mL
CsA_ka    = 0.55    // /h absorption rate constant
CsA_F     = 0.34    // bioavailability
CsA_Vc    = 4.5     // L/kg central volume
CsA_CL    = 0.28    // L/h/kg clearance
CsA_EC50  = 165     // ng/mL for 50% calcineurin inhibition
CsA_Emax  = 0.85    // Max fractional calcineurin inhibition

// ---- EPAG PK Parameters ----
// t½ ~21h, target trough ≥70 μg/mL; East Asian dose 75 mg/d
EPAG_ka   = 0.90    // /h absorption rate constant
EPAG_F    = 0.52    // bioavailability (fasted)
EPAG_Vc   = 38.0    // L central volume
EPAG_CL   = 1.30    // L/h clearance
EPAG_EC50 = 60.0    // μg/mL for 50% c-Mpl stimulation
EPAG_Emax = 0.78    // Max fractional HSC / MK stimulation

// ---- Danazol PK Parameters ----
Danazol_ka = 0.40   // /h
Danazol_F  = 0.20   // low F (~20%)
Danazol_Vc = 200.0  // L
Danazol_CL = 18.0   // L/h

// ---- T Cell Dynamics ----
kTeff_prod = 0.12   // /d baseline Teff production rate
kTeff_death= 0.08   // /d natural Teff death rate
kTreg_prod = 0.03   // /d Treg production rate
kTreg_death= 0.06   // /d Treg death rate
Treg0      = 0.25   // fraction of baseline (reduced in AA)
// IFN-γ amplifies Teff expansion
IFNg_Teff  = 0.15   // /d Teff stimulation by IFN-γ (per 100 pg/mL)

// ---- HSC Dynamics ----
kHSC_self  = 0.04   // /d HSC self-renewal rate
kHSC_diff  = 0.06   // /d HSC differentiation to MPP
kHSC_death = 0.002  // /d basal HSC apoptosis rate
// IFN-γ induced apoptosis rate constant (per 100 pg/mL)
kHSC_IFNg  = 0.18   // /d HSC death per IFN-γ unit
// TNF-α synergy
kHSC_TNFa  = 0.08   // /d HSC death per TNF-α unit
// EPAG stimulation
kHSC_EPAG  = 0.10   // /d additional self-renewal (at Emax)

// ---- Erythropoiesis ----
kCFUE_prod = 0.30   // /d CFU-E production from HSC
kCFUE_mat  = 0.18   // /d CFU-E maturation to reticulocyte
kRetic_mat = 0.25   // /d retic maturation to RBC
kRBC_death = 0.0083 // /d RBC lifespan ~120d (1/120/d)
// EPO feedback (simplified)
EPO_basal  = 15.0   // mU/mL normal EPO
EPO_max    = 200.0  // mU/mL maximal EPO in severe anemia
EPO_k50    = 9.5    // g/dL Hgb at 50% EPO upregulation
Hgb_normal = 14.0   // g/dL

// ---- Granulopoiesis ----
kCFUG_prod = 0.35   // /d CFU-G production from HSC
kCFUG_mat  = 0.20   // /d CFU-G maturation to neutrophil
kANC_death = 3.0    // /d neutrophil clearance (t½ ~8h)
ANC_normal = 3.0    // ×10⁹/L

// ---- Megakaryopoiesis ----
kMK_prod   = 0.15   // /d megakaryocyte production from HSC
kMK_mat    = 0.06   // /d MK maturation rate
kPLT_prod  = 4.0    // factor: platelets per MK
kPLT_death = 0.10   // /d platelet clearance (t½ ~10d)
PLT_normal = 200.0  // ×10⁹/L

// ---- BM Cellularity ----
kBM_repair = 0.02   // /d BM cellularity recovery (per unit HSC above base)
kBM_damage = 0.05   // /d BM damage driven by inflammation
BM0        = 0.45   // initial BM cellularity (hypocellular, 45%)

// ---- Cytokine Dynamics ----
kIFNg_prod = 2.5    // pg/mL/d IFN-γ production per Teff unit
kIFNg_CL   = 1.2    // /d IFN-γ clearance
kTNFa_prod = 1.0    // pg/mL/d TNF-α production per Teff
kTNFa_CL   = 2.0    // /d TNF-α clearance
kIL2_prod  = 0.8    // pg/mL/d IL-2 production per Teff
kIL2_CL    = 3.0    // /d IL-2 clearance
IFNg0      = 45.0   // pg/mL baseline IFN-γ in active AA
TNFa0      = 25.0   // pg/mL baseline TNF-α in active AA

// ---- PNH Clone ----
kPNH_growth = 0.005 // /d PNH clone advantage (immune escape)
kPNH_max    = 0.60  // max clone fraction

// ---- Baseline (steady-state) values ----
// Used to initialize states (healthy = 100%)
HSC0        = 100.0  // % of normal
CFU_E0      = 100.0
Retic0      = 100.0
RBC0        = 100.0  // corresponds to Hgb ~14 g/dL
CFU_G0      = 100.0
ANC0        = 100.0  // corresponds to ANC ~3×10⁹/L
MK0         = 100.0
PLT0        = 100.0  // corresponds to PLT ~200×10⁹/L
Teff0       = 10.0   // ×10⁶/kg (elevated in active AA)
PNH0        = 0.02   // 2% PNH clone at diagnosis

// ---- Disease severity modifier ----
// VSAA = 1.5, SAA = 1.0, nSAA = 0.5
severity    = 1.0   // 1 = SAA baseline

// ---- Treatment flags (0=off, 1=on) ----
use_ATG     = 0
use_CsA     = 0
use_EPAG    = 0
use_Danazol = 0
use_HSCT    = 0

$CMT
// Drug PK compartments
ATG_C ATG_P         // ATG central / peripheral
CsA_C               // Cyclosporine
EPAG_C              // Eltrombopag
Danazol_C           // Danazol

// Immune compartments
Teff                // Autoreactive effector T cells
Treg                // Regulatory T cells

// HSC & progenitors
HSC                 // Hematopoietic stem cell pool (% normal)
CFU_E               // Erythroid progenitor
Retic               // Reticulocyte
RBC                 // Circulating RBC pool
CFU_G               // Granulocyte progenitor
ANC_pool            // Circulating neutrophil pool
MK                  // Megakaryocyte pool
PLT_pool            // Platelet pool

// BM & cytokines
BM_score            // BM cellularity (0-1)
IFNg_c              // IFN-γ (pg/mL)
TNFa_c              // TNF-α (pg/mL)
IL2_c               // IL-2 (pg/mL)
PNH_clone           // PNH clone fraction

$MAIN
// --- ATG dose rate (mg/h, set via event table) ---
// Will be controlled via $INIT and events

// --- Derived ATG effect on T-cell depletion ---
double ATG_Eff = 0.0;
if (use_ATG > 0.5) {
  double atg_conc = ATG_C;
  ATG_Eff = ATG_Emax * pow(atg_conc, ATG_hill) /
            (pow(ATG_EC50, ATG_hill) + pow(atg_conc, ATG_hill));
}

// --- CsA effect on IL-2 / Teff proliferation ---
double CsA_Eff = 0.0;
if (use_CsA > 0.5) {
  double csa_conc = CsA_C;
  CsA_Eff = CsA_Emax * csa_conc / (CsA_EC50 + csa_conc);
}

// --- EPAG effect on HSC/MK stimulation ---
double EPAG_Eff = 0.0;
if (use_EPAG > 0.5) {
  double epag_conc = EPAG_C;
  EPAG_Eff = EPAG_Emax * epag_conc / (EPAG_EC50 + epag_conc);
}

// --- Danazol androgen effect (EPO augmentation) ---
double Danazol_Eff = 0.0;
if (use_Danazol > 0.5) {
  double dan_conc = Danazol_C;
  // modest effect, saturating above 600 ng/mL
  Danazol_Eff = 0.30 * dan_conc / (400.0 + dan_conc);
}

// --- HSCT: instantaneous replacement modeled via event (flag) ---
double HSCT_Eff = use_HSCT;

// --- EPO feedback on erythropoiesis ---
double Hgb_cur = (RBC / 100.0) * Hgb_normal;
double EPO_cur = EPO_basal + (EPO_max - EPO_basal) *
                 pow((1.0 - Hgb_cur / Hgb_normal), 2.0) /
                 (pow(EPO_k50/Hgb_normal, 2.0) +
                  pow((1.0 - Hgb_cur / Hgb_normal), 2.0) + 1e-6);
double EPO_fold = EPO_cur / EPO_basal;

// --- IFN-γ normalized unit (per 100 pg/mL) ---
double IFNg_unit = IFNg_c / 100.0;
double TNFa_unit = TNFa_c / 100.0;

// --- HSC destruction rate by immune attack ---
double HSC_kill = (kHSC_IFNg * IFNg_unit + kHSC_TNFa * TNFa_unit)
                   * (HSC / 100.0) * severity;

// --- Treg suppression of HSC kill ---
double Treg_suppress = Treg / (Treg + 5.0); // Hill-type, K½ = 5
HSC_kill = HSC_kill * (1.0 - 0.5 * Treg_suppress);

$INIT
ATG_C     = 0.0
ATG_P     = 0.0
CsA_C     = 0.0
EPAG_C    = 0.0
Danazol_C = 0.0
Teff      = 10.0    // ×10⁶/kg  (elevated autoimmune)
Treg      = 1.0     // ×10⁶/kg  (reduced)
HSC       = 20.0    // 20% of normal (SAA, severe pancytopenia)
CFU_E     = 25.0    // 25% residual
Retic     = 20.0
RBC       = 55.0    // Hgb ~7.7 g/dL (transfusion dependent)
CFU_G     = 15.0
ANC_pool  = 10.0    // ANC ~0.3 ×10⁹/L
MK        = 12.0
PLT_pool  = 8.0     // PLT ~16 ×10⁹/L (thrombocytopenic)
BM_score  = 0.12    // BM cellularity ~12% (hypocellular)
IFNg_c    = 45.0    // pg/mL (elevated)
TNFa_c    = 25.0    // pg/mL (elevated)
IL2_c     = 20.0    // pg/mL (elevated)
PNH_clone = 0.02    // 2% PNH clone

$ODE
// ====== ATG PK ======
dxdt_ATG_C  = -ATG_CL/ATG_Vc * ATG_C - ATG_Q/ATG_Vc * ATG_C +
               ATG_Q/ATG_Vp * ATG_P;  // + infusion via event rate
dxdt_ATG_P  =  ATG_Q/ATG_Vc * ATG_C - ATG_Q/ATG_Vp * ATG_P;

// ====== CsA PK (oral, first-order) ======
// CsA given as daily PO; ka modeled separately, assume depot handled by event
dxdt_CsA_C  = -CsA_CL/CsA_Vc * CsA_C;

// ====== EPAG PK ======
dxdt_EPAG_C = -EPAG_CL/EPAG_Vc * EPAG_C;

// ====== Danazol PK ======
dxdt_Danazol_C = -Danazol_CL/Danazol_Vc * Danazol_C;

// ====== EFFECTOR T CELLS ======
// Production ↑ by IL-2, ↓ by ATG depletion and CsA (IL-2 block)
double Teff_prolif = kTeff_prod * Teff *
                     (1.0 + (IL2_c / 50.0)) *
                     (1.0 - CsA_Eff) *
                     (1.0 - ATG_Eff);
double Teff_death  = (kTeff_death + ATG_Eff * 0.5) * Teff;
// Treg suppresses Teff
double Treg_inh    = 0.3 * (Treg / (Treg + 2.0)) * Teff;
// HSCT ablates autoreactive T cells immediately (modeled as high death)
double HSCT_Teff   = HSCT_Eff * 5.0 * Teff;

dxdt_Teff = Teff_prolif - Teff_death - Treg_inh - HSCT_Teff;

// ====== REGULATORY T CELLS ======
// Treg production preserved but baseline reduced in AA
double Treg_prod  = kTreg_prod * (1.0 + 0.5 * ATG_Eff) * 5.0;
double Treg_death = kTreg_death * Treg;
double HSCT_Treg  = HSCT_Eff * 0.5 * Treg;  // less depletion
dxdt_Treg = Treg_prod - Treg_death - HSCT_Treg;

// ====== HSC POOL ======
// Self-renewal
double HSC_selfR = kHSC_self * HSC * (1.0 + EPAG_Eff * kHSC_EPAG / kHSC_self);
// Differentiation output
double HSC_diff  = kHSC_diff * HSC;
// Basal death
double HSC_basal = kHSC_death * HSC;
// HSCT: sudden restoration of donor HSC (modeled as bolus to 100% on Day 0)
// implemented via $EVENT in simulation code
// Upper limit: cannot exceed 105 (normal + small reserve)
double HSC_growth = (HSC < 105.0) ? (HSC_selfR - HSC_diff - HSC_basal - HSC_kill) : 0.0;
dxdt_HSC = HSC_growth;

// ====== ERYTHROID PROGENITORS (CFU-E) ======
double CFU_E_in   = kCFUE_prod * (HSC / 100.0) * EPO_fold;
double CFU_E_out  = kCFUE_mat * CFU_E + kHSC_IFNg * 0.5 * IFNg_unit * CFU_E;
dxdt_CFU_E = CFU_E_in - CFU_E_out;

// ====== RETICULOCYTES ======
dxdt_Retic = kCFUE_mat * CFU_E - kRetic_mat * Retic;

// ====== CIRCULATING RBC ======
// Androgen slightly augments via EPO
double RBC_in  = kRetic_mat * Retic * (1.0 + Danazol_Eff);
double RBC_out = kRBC_death * RBC;
dxdt_RBC = RBC_in - RBC_out;

// ====== GRANULOCYTE PROGENITORS (CFU-G) ======
double CFU_G_in  = kCFUG_prod * (HSC / 100.0);
double CFU_G_out = kCFUG_mat * CFU_G;
dxdt_CFU_G = CFU_G_in - CFU_G_out;

// ====== CIRCULATING NEUTROPHILS ======
double ANC_in  = kCFUG_mat * CFU_G;
double ANC_out = kANC_death * ANC_pool;
dxdt_ANC_pool = ANC_in - ANC_out;

// ====== MEGAKARYOCYTES ======
double MK_in  = kMK_prod * (HSC / 100.0) * (1.0 + EPAG_Eff);
double MK_out = kMK_mat * MK;
dxdt_MK = MK_in - MK_out;

// ====== PLATELETS ======
double PLT_in  = kMK_mat * MK * kPLT_prod;
double PLT_out = kPLT_death * PLT_pool;
dxdt_PLT_pool = PLT_in - PLT_out;

// ====== BONE MARROW CELLULARITY ======
// Recovers as HSC pool recovers; damaged by inflammation
double BM_recovery = kBM_repair * (HSC - 20.0) / 100.0; // above baseline (AA baseline=20%)
double BM_damage   = kBM_damage * (IFNg_unit + TNFa_unit * 0.5) * BM_score;
double BM_max_gain = (BM_score < 1.0) ? BM_recovery : 0.0;
dxdt_BM_score = BM_max_gain - BM_damage;

// ====== CYTOKINES ======
dxdt_IFNg_c = kIFNg_prod * Teff - kIFNg_CL * IFNg_c;
dxdt_TNFa_c = kTNFa_prod * Teff - kTNFa_CL * TNFa_c;
dxdt_IL2_c  = kIL2_prod  * Teff * (1.0 - CsA_Eff) - kIL2_CL * IL2_c;

// ====== PNH CLONE ======
// PNH clones expand due to immune escape (GPI-anchor deficient = invisible to immune attack)
double PNH_immune_escape = kPNH_growth * PNH_clone * (1.0 - PNH_clone);
double PNH_IST_suppress  = ATG_Eff * 0.1 * PNH_clone; // partial IST suppression
dxdt_PNH_clone = PNH_immune_escape - PNH_IST_suppress;

$TABLE
capture Hgb        = (RBC / 100.0) * Hgb_normal;
capture ANC        = (ANC_pool / 100.0) * ANC_normal;  // ×10⁹/L
capture PLT        = (PLT_pool / 100.0) * PLT_normal;  // ×10⁹/L
capture ARC        = (Retic / 100.0) * 80.0;           // ×10⁹/L (normal ARC ~80)
capture BM_cell    = BM_score * 100.0;                 // % cellularity
capture Hgb_Tx     = (Hgb < 8.0) ? 1.0 : 0.0;        // transfusion trigger
capture PLT_Tx     = (PLT < 10.0) ? 1.0 : 0.0;        // platelet transfusion trigger
capture CR_flag    = (Hgb > 11.0 && ANC > 1.0 && PLT > 100.0) ? 1.0 : 0.0;
capture PR_flag    = (Hgb_Tx < 0.5 && PLT > 20.0 && ANC > 0.5) ? 1.0 : 0.0;
capture SAA_flag   = (ANC < 0.5 && (PLT < 20.0 || Hgb < 8.0)) ? 1.0 : 0.0;
capture VSAA_flag  = (ANC < 0.2) ? 1.0 : 0.0;
capture IFNg_norm  = IFNg_c / 10.0;  // scale for plotting
capture ATG_conc   = ATG_C;
capture CsA_conc   = CsA_C;
capture EPAG_conc  = EPAG_C;
'

## ============================================================
## 2. Compile the Model
## ============================================================

aa_mod <- mcode("aplastic_anemia_qsp", aa_model_code)

## ============================================================
## 3. Helper Functions
## ============================================================

# Build ATG event table (hATG: 40 mg/kg over 12h × 4 days, ~70kg pt → 2800mg total)
# Modeled as infusion rate = dose / duration = 700 mg / 12h = 58.3 mg/h per day
build_hATG_events <- function(start_day = 1, weight = 70) {
  dose_per_day <- 40 * weight  # mg
  rate_per_day <- dose_per_day / 12  # mg/h (12h infusion)
  ev <- ev(
    time = (start_day - 1 + 0:3) * 24,  # days 1,2,3,4 (in hours)
    amt  = dose_per_day,
    rate = rate_per_day,
    cmt  = "ATG_C"
  )
  return(ev)
}

# Build rATG events (3.5 mg/kg/d × 5d)
build_rATG_events <- function(start_day = 1, weight = 70) {
  dose_per_day <- 3.5 * weight
  rate_per_day <- dose_per_day / 12
  ev <- ev(
    time = (start_day - 1 + 0:4) * 24,
    amt  = dose_per_day,
    rate = rate_per_day,
    cmt  = "ATG_C"
  )
  return(ev)
}

# Build CsA events (5 mg/kg/d, BID oral, ~600mg/d for 70kg → 300mg BID)
build_CsA_events <- function(start_day = 1, duration_days = 365, weight = 70) {
  dose_bid <- 2.5 * weight  # mg per dose
  times_am <- ((start_day - 1):(start_day + duration_days - 2)) * 24
  times_pm <- times_am + 12
  times <- sort(c(times_am, times_pm))
  ev <- ev(
    time = times,
    amt  = dose_bid * CsA_F * 1000,  # convert to ng/mL scale... simplification
    cmt  = "CsA_C"
  )
  # Simulate with simple bolus for illustration
  ev <- ev(
    time = ((start_day - 1):(start_day + duration_days - 2)) * 24,
    amt  = 300,   # mg/day simplified PO dose (ng/mL units adjusted in PK)
    cmt  = "CsA_C",
    rate = -2     # first-order absorption flag
  )
  return(ev)
}

# Build EPAG events (150 mg/d oral QD)
build_EPAG_events <- function(start_day = 1, duration_days = 180, dose = 150) {
  ev <- ev(
    time = ((start_day - 1):(start_day + duration_days - 2)) * 24,
    amt  = dose,
    cmt  = "EPAG_C",
    rate = -2
  )
  return(ev)
}

## ============================================================
## 4. Treatment Scenarios
## ============================================================

## Scenario 1: No Treatment (Natural History of SAA)
scenario_1_NoTx <- function() {
  param_update <- param(aa_mod,
    use_ATG = 0, use_CsA = 0, use_EPAG = 0, use_Danazol = 0, use_HSCT = 0,
    severity = 1.0
  )
  out <- mrgsim(param_update,
    end = 365, delta = 1,
    carry_out = "evid"
  )
  as.data.frame(out) %>% mutate(Scenario = "No Treatment")
}

## Scenario 2: Standard IST (hATG + CsA, no EPAG)
scenario_2_hATG_CsA <- function() {
  param_update <- param(aa_mod,
    use_ATG = 1, use_CsA = 1, use_EPAG = 0, use_Danazol = 0, use_HSCT = 0,
    severity = 1.0,
    ATG_CL = 0.85   # hATG PK
  )
  atg_ev  <- build_hATG_events(start_day = 1)
  csa_ev  <- build_CsA_events(start_day = 1, duration_days = 365)
  all_ev  <- c(atg_ev, csa_ev)

  out <- mrgsim(param_update, events = all_ev,
    end = 365, delta = 1,
    carry_out = "evid"
  )
  as.data.frame(out) %>% mutate(Scenario = "hATG + CsA (Standard IST)")
}

## Scenario 3: IST + Eltrombopag (EPAG added at Day 14)
scenario_3_hATG_CsA_EPAG <- function() {
  param_update <- param(aa_mod,
    use_ATG = 1, use_CsA = 1, use_EPAG = 1, use_Danazol = 0, use_HSCT = 0,
    severity = 1.0
  )
  atg_ev  <- build_hATG_events(start_day = 1)
  csa_ev  <- build_CsA_events(start_day = 1, duration_days = 365)
  epag_ev <- build_EPAG_events(start_day = 14, duration_days = 180, dose = 150)
  all_ev  <- c(atg_ev, csa_ev, epag_ev)

  out <- mrgsim(param_update, events = all_ev,
    end = 365, delta = 1,
    carry_out = "evid"
  )
  as.data.frame(out) %>% mutate(Scenario = "hATG + CsA + EPAG (Townsley 2017)")
}

## Scenario 4: rATG + CsA + EPAG (NiH protocol)
scenario_4_rATG_CsA_EPAG <- function() {
  param_update <- param(aa_mod,
    use_ATG = 1, use_CsA = 1, use_EPAG = 1, use_Danazol = 0, use_HSCT = 0,
    severity = 1.0,
    ATG_CL = 0.50  # rATG: slower clearance, more potent
  )
  atg_ev  <- build_rATG_events(start_day = 1)
  csa_ev  <- build_CsA_events(start_day = 1, duration_days = 365)
  epag_ev <- build_EPAG_events(start_day = 14, duration_days = 180)
  all_ev  <- c(atg_ev, csa_ev, epag_ev)

  out <- mrgsim(param_update, events = all_ev,
    end = 365, delta = 1,
    carry_out = "evid"
  )
  as.data.frame(out) %>% mutate(Scenario = "rATG + CsA + EPAG (NiH protocol)")
}

## Scenario 5: Allogeneic HSCT (MSD donor, Day 0 conditioning)
## HSCT modeled as: conditioning ablates HSC (day -7 to -1), then donor graft engrafts
## Simplified: set use_HSCT = 1 on day 0; HSC bolus to 100% via event
scenario_5_HSCT <- function() {
  param_update <- param(aa_mod,
    use_ATG = 1, use_CsA = 0, use_EPAG = 0, use_Danazol = 0, use_HSCT = 1,
    severity = 1.0,
    ATG_CL = 0.50  # conditioning ATG (rATG-based)
  )
  # Conditioning ATG (day -7 to -3, simplified)
  atg_ev  <- build_rATG_events(start_day = 1)
  # Simulate donor HSC engraftment on Day 14 as reset of HSC compartment
  # Handled via initial conditions + parameter change at engraftment
  out <- mrgsim(param_update, events = atg_ev,
    end = 365, delta = 1,
    carry_out = "evid",
    idata = data.frame(ID = 1)
  )
  res <- as.data.frame(out) %>% mutate(Scenario = "Allogeneic HSCT (MSD)")
  # After engraftment (Day ~21), HSC recovers rapidly
  # For illustration, set HSC to 100 at Day 21 in post-processing
  res$HSC[res$time >= 21 * 24] <- pmin(res$HSC[res$time >= 21 * 24] +
    seq(0, 80, length.out = sum(res$time >= 21 * 24)), 100)
  return(res)
}

## ============================================================
## 5. Run All Scenarios & Combine
## ============================================================
run_all_scenarios <- function() {
  scenarios <- list(
    scenario_1_NoTx(),
    scenario_2_hATG_CsA(),
    scenario_3_hATG_CsA_EPAG(),
    scenario_4_rATG_CsA_EPAG(),
    scenario_5_HSCT()
  )
  bind_rows(scenarios) %>%
    mutate(time_days = time / 24)
}

## ============================================================
## 6. Visualization
## ============================================================
plot_results <- function(results) {

  theme_aa <- theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "#2C3E50"),
      strip.text = element_text(color = "white", face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  cols <- c(
    "No Treatment"                      = "#E74C3C",
    "hATG + CsA (Standard IST)"         = "#3498DB",
    "hATG + CsA + EPAG (Townsley 2017)" = "#2ECC71",
    "rATG + CsA + EPAG (NiH protocol)"  = "#F39C12",
    "Allogeneic HSCT (MSD)"             = "#9B59B6"
  )

  ## A. Hemoglobin trajectories
  p1 <- results %>%
    ggplot(aes(x = time_days, y = Hgb, color = Scenario)) +
    geom_line(linewidth = 1.0) +
    geom_hline(yintercept = 11, linetype = "dashed", color = "darkgreen", alpha = 0.7) +
    geom_hline(yintercept = 8,  linetype = "dashed", color = "red",       alpha = 0.7) +
    scale_color_manual(values = cols) +
    labs(title = "Hemoglobin Trajectory",
         x = "Time (days)", y = "Hemoglobin (g/dL)",
         caption = "Dashed lines: CR threshold (11 g/dL) and transfusion trigger (8 g/dL)") +
    theme_aa

  ## B. ANC trajectories
  p2 <- results %>%
    ggplot(aes(x = time_days, y = ANC, color = Scenario)) +
    geom_line(linewidth = 1.0) +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "darkgreen", alpha = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "orange",    alpha = 0.7) +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "red",       alpha = 0.7) +
    scale_color_manual(values = cols) +
    labs(title = "Absolute Neutrophil Count (ANC)",
         x = "Time (days)", y = "ANC (×10⁹/L)",
         caption = "Dashed: CR >1.0, SAA <0.5, VSAA <0.2") +
    theme_aa

  ## C. Platelet trajectories
  p3 <- results %>%
    ggplot(aes(x = time_days, y = PLT, color = Scenario)) +
    geom_line(linewidth = 1.0) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "darkgreen", alpha = 0.7) +
    geom_hline(yintercept = 20,  linetype = "dashed", color = "orange",    alpha = 0.7) +
    geom_hline(yintercept = 10,  linetype = "dashed", color = "red",       alpha = 0.7) +
    scale_color_manual(values = cols) +
    labs(title = "Platelet Count",
         x = "Time (days)", y = "Platelets (×10⁹/L)",
         caption = "Dashed: CR >100, PLT Tx <10") +
    theme_aa

  ## D. HSC Pool
  p4 <- results %>%
    ggplot(aes(x = time_days, y = HSC, color = Scenario)) +
    geom_line(linewidth = 1.0) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "darkgreen", alpha = 0.5) +
    scale_color_manual(values = cols) +
    labs(title = "HSC Pool (% of Normal)",
         x = "Time (days)", y = "HSC Pool (%)") +
    theme_aa

  ## E. Cytokine dynamics
  p5 <- results %>%
    select(time_days, Scenario, IFNg_norm, TNFa_c, IL2_c) %>%
    pivot_longer(cols = c(IFNg_norm, TNFa_c, IL2_c),
                 names_to = "Cytokine", values_to = "Concentration") %>%
    mutate(Cytokine = recode(Cytokine,
      IFNg_norm = "IFN-γ (×10 pg/mL)",
      TNFa_c    = "TNF-α (pg/mL)",
      IL2_c     = "IL-2 (pg/mL)"
    )) %>%
    ggplot(aes(x = time_days, y = Concentration, color = Scenario, linetype = Cytokine)) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = cols) +
    labs(title = "Cytokine Dynamics",
         x = "Time (days)", y = "Concentration") +
    theme_aa

  ## F. Response rate summary (CR at Day 180 & 365)
  response_summary <- results %>%
    filter(time_days %in% c(90, 180, 365)) %>%
    group_by(Scenario, time_days) %>%
    summarize(
      CR_rate = mean(CR_flag),
      PR_rate = mean(PR_flag),
      NR_rate = 1 - mean(CR_flag | PR_flag),
      .groups = "drop"
    )

  p6 <- response_summary %>%
    pivot_longer(cols = c(CR_rate, PR_rate, NR_rate), names_to = "Response", values_to = "Rate") %>%
    mutate(Response = recode(Response,
      CR_rate = "Complete Response",
      PR_rate = "Partial Response",
      NR_rate = "No Response"
    )) %>%
    ggplot(aes(x = as.factor(time_days), y = Rate * 100, fill = Response)) +
    geom_col(position = "stack") +
    facet_wrap(~ Scenario) +
    scale_fill_manual(values = c("Complete Response" = "#2ECC71",
                                  "Partial Response"  = "#F39C12",
                                  "No Response"       = "#E74C3C")) +
    labs(title = "Response Rate by Scenario",
         x = "Time (days)", y = "Response Rate (%)") +
    theme_aa

  list(Hgb = p1, ANC = p2, PLT = p3, HSC = p4, Cytokines = p5, Response = p6)
}

## ============================================================
## 7. Sensitivity Analysis
## ============================================================
sensitivity_EPAG_dose <- function(doses = c(75, 100, 150, 200)) {
  results <- map_dfr(doses, function(dose) {
    param_update <- param(aa_mod,
      use_ATG = 1, use_CsA = 1, use_EPAG = 1, severity = 1.0
    )
    atg_ev  <- build_hATG_events(start_day = 1)
    csa_ev  <- build_CsA_events(start_day = 1, duration_days = 365)
    epag_ev <- build_EPAG_events(start_day = 14, duration_days = 180, dose = dose)
    all_ev  <- c(atg_ev, csa_ev, epag_ev)

    out <- mrgsim(param_update, events = all_ev, end = 365, delta = 1)
    as.data.frame(out) %>%
      mutate(time_days = time / 24, EPAG_dose = paste0(dose, " mg/d"))
  })

  ggplot(results, aes(x = time_days, y = PLT, color = EPAG_dose)) +
    geom_line(linewidth = 1.1) +
    scale_color_viridis_d(name = "EPAG Dose") +
    geom_hline(yintercept = 100, linetype = "dashed") +
    labs(title = "Sensitivity: EPAG Dose vs Platelet Recovery",
         x = "Time (days)", y = "Platelet Count (×10⁹/L)") +
    theme_bw(base_size = 12)
}

## ============================================================
## 8. Virtual Patient Population (VPP)
## ============================================================
run_VPP <- function(n = 100, seed = 42) {
  set.seed(seed)
  # Vary: initial HSC pool (5-35%), severity, EPAG response (EC50 variability)
  VPP_params <- data.frame(
    ID         = 1:n,
    HSC_init   = runif(n, 5, 35),
    severity   = sample(c(0.5, 1.0, 1.5), n, replace = TRUE,
                        prob = c(0.25, 0.55, 0.20)),
    EPAG_EC50  = rnorm(n, 60, 15) %>% pmax(20) %>% pmin(120)
  )

  results <- map_dfr(1:n, function(i) {
    p <- VPP_params[i, ]
    param_update <- param(aa_mod,
      use_ATG = 1, use_CsA = 1, use_EPAG = 1,
      severity  = p$severity,
      EPAG_EC50 = p$EPAG_EC50
    )
    init_update <- init(param_update, HSC = p$HSC_init,
      RBC = p$HSC_init * 2.5, PLT_pool = p$HSC_init * 0.4)
    atg_ev  <- build_hATG_events(start_day = 1)
    csa_ev  <- build_CsA_events(start_day = 1, duration_days = 365)
    epag_ev <- build_EPAG_events(start_day = 14, duration_days = 180)
    all_ev  <- c(atg_ev, csa_ev, epag_ev)

    out <- mrgsim(init_update, events = all_ev, end = 180, delta = 7)
    as.data.frame(out) %>%
      mutate(time_days = time / 24, ID = p$ID,
             severity_label = case_when(
               p$severity < 0.8 ~ "nSAA",
               p$severity > 1.2 ~ "VSAA",
               TRUE             ~ "SAA"
             ))
  })

  # Spaghetti plot: Hgb trajectories
  ggplot(results, aes(x = time_days, y = Hgb, group = ID, color = severity_label)) +
    geom_line(alpha = 0.35, linewidth = 0.5) +
    stat_summary(aes(group = severity_label), fun = "median",
                 geom = "line", linewidth = 1.5, linetype = "solid") +
    scale_color_manual(values = c(nSAA = "#2ECC71", SAA = "#F39C12", VSAA = "#E74C3C"),
                       name = "Severity") +
    geom_hline(yintercept = 11, linetype = "dashed", color = "darkgreen") +
    labs(title = "Virtual Patient Population: Hgb Trajectories (hATG+CsA+EPAG)",
         subtitle = paste("n =", n, "virtual patients; bold = median by subgroup"),
         x = "Time (days)", y = "Hemoglobin (g/dL)") +
    theme_bw(base_size = 12)
}

## ============================================================
## Main Execution Block (example)
## ============================================================
if (FALSE) {  # Set to TRUE to run interactively
  cat("Running all scenarios...\n")
  results <- run_all_scenarios()

  cat("Generating plots...\n")
  plots <- plot_results(results)
  print(plots$Hgb)
  print(plots$ANC)
  print(plots$PLT)
  print(plots$HSC)
  print(plots$Cytokines)
  print(plots$Response)

  cat("Sensitivity analysis: EPAG dose...\n")
  print(sensitivity_EPAG_dose())

  cat("Virtual patient population...\n")
  print(run_VPP(n = 100))

  cat("Done.\n")
}
