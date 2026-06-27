##############################################################################
# Hypereosinophilic Syndrome (HES) — QSP Model (mrgsolve)
# Author : Claude Code Routine (CCR)
# Date   : 2026-06-27
# ─────────────────────────────────────────────────────────────────────────────
# Model structure
#   PK  : Mepolizumab (2-CMT SC + TMDD)
#         Benralizumab (2-CMT SC + ADCC)
#         Imatinib    (1-CMT PO oral)
#         Prednisolone (1-CMT PO)
#   PD  : IL-5 kinetics (TMDD-linked)
#         Bone-marrow eosinophilopoiesis (EoP → EoImm → EoMat-BM → Blood)
#         Peripheral blood AEC dynamics
#         Cardiac fibrosis progression (Löffler → EMF → RCM)
#         Pulmonary infiltration score
#
# Clinical calibration references
#   Mepolizumab PK : Roufosse 2020 NEJM; Farne & Cates 2017 Cochrane
#   Benralizumab   : Bleecker 2016 NEJM; Kolbeck 2010 J Allergy Clin Immunol
#   Imatinib       : Gleich 2002 NEJM; Cools 2003 NEJM
#   Prednisolone   : Ogbogu 2009 J Allergy Clin Immunol
#   AEC kinetics   : Ackerman 1981 Blood; Sanderson 1992 Blood
##############################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ─────────────────────────────────────────────────────────────────────────────
# MODEL BLOCK
# ─────────────────────────────────────────────────────────────────────────────
hes_model_code <- '
$PROB HES QSP Model — Eosinophil Kinetics + Multi-drug PK/PD

$PARAM @annotated
// ── Disease baseline ──────────────────────────────────────────────────────
AEC0      : 3000  : Baseline absolute eosinophil count (cells/µL)
IL5_0     :  5.0  : Baseline serum IL-5 (pg/mL)
kprod_IL5 : 0.05  : IL-5 production rate constant (pg/mL/h)
kdeg_IL5  : 0.01  : IL-5 baseline degradation rate (1/h)

// ── BM Eosinophilopoiesis ─────────────────────────────────────────────────
kprol_EoP  : 0.050  : EoP proliferation rate (1/h)
kmat_EoP   : 0.020  : EoP → Immature Eo maturation rate (1/h)
kmat_EoI   : 0.015  : Immature → Mature BM eo maturation (1/h)
krel_BM    : 0.030  : BM release rate to blood (1/h)
kdeath_EoP : 0.008  : EoP baseline apoptosis rate (1/h)

// ── Peripheral blood Eo kinetics ──────────────────────────────────────────
kel_Eo     : 0.038  : Eo elimination from blood (1/h; t½~18h)
Vd_Eo      : 5.0    : Distribution volume scaling for Eo (L)
k_tissue   : 0.010  : Eo migration to tissue (1/h)

// ── IL-5 effect on eosinophils ────────────────────────────────────────────
Emax_IL5   : 3.0    : Maximum fold-increase in EoP proliferation by IL-5
EC50_IL5   : 2.0    : EC50 of IL-5 (pg/mL) for EoP proliferation
Emax_surv  : 0.7    : Maximum fraction reduction in Eo apoptosis by IL-5
EC50_surv  : 1.5    : EC50 of IL-5 for blood Eo survival extension

// ── Cardiac fibrosis progression ─────────────────────────────────────────
k_fibrosis : 0.0002 : Rate of cardiac fibrosis accumulation (1/h per Eo unit)
krev_fibro : 0.0001 : Fibrosis reversibility rate (1/h; partial)
FIBROSIS_0 : 0.05   : Baseline fibrosis score (0–1 scale)

// ── Pulmonary infiltration ────────────────────────────────────────────────
k_pulm     : 0.0003 : Rate of pulmonary infiltration accumulation (1/h)
krev_pulm  : 0.0002 : Pulmonary infiltration reversibility rate (1/h)
PULM_0     : 0.10   : Baseline pulmonary infiltration score (0–1)

// ── Mepolizumab PK (2-CMT SC, anti-IL-5 IgG1) ───────────────────────────
F_MEPO     : 0.80   : SC bioavailability of mepolizumab
ka_MEPO    : 0.0035 : Absorption rate constant SC (1/h; ~2.8 day Tmax)
CL_MEPO    : 0.0072 : Clearance (L/h; Vss~7L, t½~22 days)
Vc_MEPO    : 3.5    : Central volume (L)
Q_MEPO     : 0.012  : Inter-compartmental clearance (L/h)
Vp_MEPO    : 3.5    : Peripheral volume (L)
kon_MEPO   : 10.0   : IL-5 binding on-rate (1/[nM·h])
koff_MEPO  : 0.001  : IL-5 binding off-rate (1/h)
kdeg_TMDD  : 0.05   : TMDD complex degradation rate (1/h)
MW_MEPO    : 149000 : Mepolizumab molecular weight (Da)

// ── Benralizumab PK (2-CMT SC, anti-IL-5Rα) ─────────────────────────────
F_BENRA    : 0.59   : SC bioavailability of benralizumab
ka_BENRA   : 0.0060 : Absorption rate constant SC (1/h)
CL_BENRA   : 0.0040 : Clearance (L/h; t½~15 days)
Vc_BENRA   : 3.0    : Central volume (L)
Q_BENRA    : 0.015  : Inter-compartmental clearance (L/h)
Vp_BENRA   : 3.8    : Peripheral volume (L)
Emax_ADCC  : 0.95   : Maximum ADCC-mediated Eo depletion
EC50_ADCC  : 0.01   : EC50 for ADCC (µg/mL)
ADCC_hill  : 1.5    : Hill coefficient for ADCC

// ── Imatinib PK (1-CMT PO) ───────────────────────────────────────────────
F_IMAT     : 0.98   : Oral bioavailability of imatinib
ka_IMAT    : 0.30   : Absorption rate constant PO (1/h)
CL_IMAT    : 12.0   : Clearance (L/h; CYP3A4)
Vc_IMAT    : 110.0  : Central volume (L; Vd ~4.5 L/kg)
IMAT_IC50  : 0.10   : IC50 for FIP1L1-PDGFRA inhibition (µg/mL)
IMAT_hill  : 1.0    : Hill coefficient for imatinib effect

// ── Prednisolone PK/PD (1-CMT PO) ────────────────────────────────────────
F_PRED     : 0.82   : Oral bioavailability of prednisolone
ka_PRED    : 0.80   : Absorption rate constant PO (1/h; rapid)
CL_PRED    : 3.5    : Clearance (L/h; CYP3A4 + hepatic)
Vc_PRED    : 38.0   : Central volume (L)
Emax_PRED  : 0.80   : Maximum IL-5 suppression by prednisolone
EC50_PRED  : 0.05   : EC50 prednisolone for IL-5 suppression (µg/mL)
Emax_APO   : 0.60   : Maximum promotion of Eo apoptosis (prednisolone)
EC50_APO   : 0.08   : EC50 prednisolone for apoptosis promotion (µg/mL)

// ── Simulation flags ─────────────────────────────────────────────────────
USE_MEPO   : 0  : 1 = use mepolizumab
USE_BENRA  : 0  : 1 = use benralizumab
USE_IMAT   : 0  : 1 = use imatinib (clonal HES)
USE_PRED   : 0  : 1 = use prednisolone
CLONAL_HES : 0  : 1 = FIP1L1-PDGFRA+ clonal HES (enhanced EoP prolif)
CLONAL_FOLD: 5.0 : Fold-increase in EoP proliferation in clonal HES

$CMT @annotated
// Mepolizumab PK
MEPO_DEPOT : Mepolizumab SC depot (µg)
MEPO_C1    : Mepolizumab central (µg/L)
MEPO_C2    : Mepolizumab peripheral (µg)
TMDD       : Mepolizumab-IL5 complex (nM)

// Benralizumab PK
BENRA_DEPOT : Benralizumab SC depot (µg)
BENRA_C1    : Benralizumab central (µg/L)
BENRA_C2    : Benralizumab peripheral (µg)

// Imatinib PK
IMAT_GUT    : Imatinib gut (mg)
IMAT_C      : Imatinib central (mg/L = µg/mL)

// Prednisolone PK
PRED_GUT    : Prednisolone gut (mg)
PRED_C      : Prednisolone central (mg/L = µg/mL)

// Disease: Eosinophil compartments
EoP         : Eosinophil progenitors in BM (relative units)
EoI         : Immature eosinophils in BM (relative units)
EoM_BM      : Mature eosinophils in BM (relative units)
EO_BLOOD    : Absolute eosinophil count in blood (cells/µL)
IL5         : Serum IL-5 (pg/mL)

// Tissue damage compartments
FIBROSIS    : Cardiac fibrosis score (0–1)
PULM_SCORE  : Pulmonary infiltration score (0–1)

$INIT @annotated
MEPO_DEPOT  : 0    : Initial mepolizumab depot
MEPO_C1     : 0    : Initial mepolizumab central
MEPO_C2     : 0    : Initial mepolizumab peripheral
TMDD        : 0    : Initial TMDD complex

BENRA_DEPOT : 0    : Initial benralizumab depot
BENRA_C1    : 0    : Initial benralizumab central
BENRA_C2    : 0    : Initial benralizumab peripheral

IMAT_GUT    : 0    : Initial imatinib gut
IMAT_C      : 0    : Initial imatinib central

PRED_GUT    : 0    : Initial prednisolone gut
PRED_C      : 0    : Initial prednisolone central

EoP         : 100  : Baseline EoP (arbitrary units; calibrated to AEC0)
EoI         : 80   : Baseline immature Eo
EoM_BM      : 150  : Baseline mature BM Eo
EO_BLOOD    : 3000 : Baseline blood AEC (cells/µL)
IL5         : 5.0  : Baseline serum IL-5 (pg/mL)

FIBROSIS    : 0.05 : Initial cardiac fibrosis
PULM_SCORE  : 0.10 : Initial pulmonary score

$ODE
// ─── IL-5 Dynamics (TMDD-linked) ─────────────────────────────────────────
// IL-5 bound by mepolizumab reduces free IL-5
double MEPO_Cp  = MEPO_C1 / Vc_MEPO;  // concentration µg/mL → nM ≈ Cp/MW*1e6
double IL5_free = IL5 > 0 ? IL5 : 0;

// Prednisolone IL-5 suppression
double PRED_Cp = PRED_C;
double pred_il5_inhib = 0.0;
if (USE_PRED > 0.5) {
    pred_il5_inhib = Emax_PRED * PRED_Cp / (EC50_PRED + PRED_Cp);
}

// Imatinib inhibition of clonal EoP expansion
double IMAT_Cp = IMAT_C;
double imat_inhib = 0.0;
if (USE_IMAT > 0.5 && CLONAL_HES > 0.5) {
    imat_inhib = IMAT_Cp / (IMAT_IC50 + IMAT_Cp);
}

// Clonal HES enhancement factor (FIP1L1-PDGFRA)
double clonal_factor = 1.0;
if (CLONAL_HES > 0.5) {
    clonal_factor = CLONAL_FOLD * (1.0 - imat_inhib);
}

// IL-5 production (modified by prednisolone and mepolizumab TMDD)
double PROD_IL5 = kprod_IL5 * (1.0 - pred_il5_inhib);
// IL-5 degradation (including TMDD sequestration by mepolizumab)
double DEG_IL5 = kdeg_IL5 * IL5;
double TMDD_on  = 0.0;
double TMDD_off = 0.0;
if (USE_MEPO > 0.5) {
    TMDD_on  = kon_MEPO * (MEPO_C1/Vc_MEPO) * IL5_free;
    TMDD_off = koff_MEPO * TMDD;
}
dxdt_IL5 = PROD_IL5 - DEG_IL5 - TMDD_on + TMDD_off;

// ─── Mepolizumab PK (2-CMT SC + TMDD) ───────────────────────────────────
double MEPO_ABS = (USE_MEPO > 0.5) ? (ka_MEPO * F_MEPO * MEPO_DEPOT) : 0.0;
dxdt_MEPO_DEPOT = -ka_MEPO * MEPO_DEPOT;
dxdt_MEPO_C1    = MEPO_ABS
                  - (CL_MEPO/Vc_MEPO) * MEPO_C1
                  - Q_MEPO * (MEPO_C1/Vc_MEPO - MEPO_C2/Vp_MEPO)
                  - TMDD_on * Vc_MEPO + TMDD_off * Vc_MEPO;
dxdt_MEPO_C2    = Q_MEPO * (MEPO_C1/Vc_MEPO - MEPO_C2/Vp_MEPO);
dxdt_TMDD       = TMDD_on - TMDD_off - kdeg_TMDD * TMDD;

// ─── Benralizumab PK (2-CMT SC) ─────────────────────────────────────────
double BENRA_ABS = (USE_BENRA > 0.5) ? (ka_BENRA * F_BENRA * BENRA_DEPOT) : 0.0;
dxdt_BENRA_DEPOT = -ka_BENRA * BENRA_DEPOT;
dxdt_BENRA_C1    = BENRA_ABS
                   - (CL_BENRA/Vc_BENRA) * BENRA_C1
                   - Q_BENRA * (BENRA_C1/Vc_BENRA - BENRA_C2/Vp_BENRA);
dxdt_BENRA_C2    = Q_BENRA * (BENRA_C1/Vc_BENRA - BENRA_C2/Vp_BENRA);

// ─── Imatinib PK (1-CMT PO) ─────────────────────────────────────────────
dxdt_IMAT_GUT = -(USE_IMAT > 0.5 ? ka_IMAT : 0.0) * IMAT_GUT;
dxdt_IMAT_C   =  (USE_IMAT > 0.5 ? ka_IMAT * F_IMAT * IMAT_GUT : 0.0)
                 - (CL_IMAT/Vc_IMAT) * IMAT_C;

// ─── Prednisolone PK (1-CMT PO) ─────────────────────────────────────────
dxdt_PRED_GUT = -(USE_PRED > 0.5 ? ka_PRED : 0.0) * PRED_GUT;
dxdt_PRED_C   =  (USE_PRED > 0.5 ? ka_PRED * F_PRED * PRED_GUT : 0.0)
                 - (CL_PRED/Vc_PRED) * PRED_C;

// ─── BM Eosinophilopoiesis ───────────────────────────────────────────────
// IL-5 Emax on EoP proliferation
double IL5_Eprol = Emax_IL5 * IL5_free / (EC50_IL5 + IL5_free);
// Prednisolone apoptosis promotion
double pred_apo = 0.0;
if (USE_PRED > 0.5) {
    pred_apo = Emax_APO * PRED_Cp / (EC50_APO + PRED_Cp);
}
// Benralizumab ADCC effect on blood Eo
double BENRA_Cp = BENRA_C1 / Vc_BENRA;
double adcc_eff = 0.0;
if (USE_BENRA > 0.5) {
    double BENRA_n = pow(BENRA_Cp, ADCC_hill);
    double EC50_n  = pow(EC50_ADCC, ADCC_hill);
    adcc_eff = Emax_ADCC * BENRA_n / (EC50_n + BENRA_n);
}

// EoP compartment
double PROL_RATE = kprol_EoP * (1.0 + IL5_Eprol) * clonal_factor * EoP;
double DEATH_EoP = kdeath_EoP * (1.0 + pred_apo) * EoP;
dxdt_EoP = PROL_RATE - kmat_EoP * EoP - DEATH_EoP;

// Immature Eo
dxdt_EoI = kmat_EoP * EoP - kmat_EoI * EoI;

// Mature BM Eo
dxdt_EoM_BM = kmat_EoI * EoI - krel_BM * EoM_BM;

// ─── Peripheral blood Eo (AEC) ───────────────────────────────────────────
// IL-5 survival extension (reduces elimination)
double IL5_surv  = Emax_surv * IL5_free / (EC50_surv + IL5_free);
double kel_eff   = kel_Eo * (1.0 - IL5_surv);
// ADCC depletion
double adcc_elim = adcc_eff * EO_BLOOD;

dxdt_EO_BLOOD = krel_BM * EoM_BM
                - kel_eff * EO_BLOOD
                - k_tissue * EO_BLOOD
                - adcc_elim;

// ─── Cardiac Fibrosis (Löffler → EMF) ────────────────────────────────────
// Driven by tissue eosinophil burden (proportional to AEC)
double EO_norm  = EO_BLOOD / 500.0;  // normalize to upper normal
double FIBRO_IN = k_fibrosis * (EO_norm > 1.0 ? EO_norm - 1.0 : 0.0) * (1.0 - FIBROSIS);
double FIBRO_REV = krev_fibro * FIBROSIS;
dxdt_FIBROSIS = FIBRO_IN - FIBRO_REV;

// ─── Pulmonary Infiltration Score ────────────────────────────────────────
double PULM_IN  = k_pulm * (EO_norm > 1.0 ? EO_norm - 1.0 : 0.0) * (1.0 - PULM_SCORE);
double PULM_REV = krev_pulm * PULM_SCORE;
dxdt_PULM_SCORE = PULM_IN - PULM_REV;

$TABLE
// ── Derived PK outputs ────────────────────────────────────────────────────
double MEPO_Cobs  = MEPO_C1 / Vc_MEPO;   // µg/mL
double BENRA_Cobs = BENRA_C1 / Vc_BENRA; // µg/mL
double IMAT_Cobs  = IMAT_C;              // µg/mL
double PRED_Cobs  = PRED_C;              // µg/mL

// ── Clinical biomarkers ───────────────────────────────────────────────────
double AEC_obs    = EO_BLOOD;             // cells/µL
double PERCENT_CHG = (AEC_obs - AEC0) / AEC0 * 100.0;  // % change from baseline
double RESP_300   = (AEC_obs < 300.0) ? 1.0 : 0.0;     // AEC < 300 response
double RESP_1500  = (AEC_obs < 1500.0) ? 1.0 : 0.0;    // AEC < 1500 control
double CARDIAC_SCORE = FIBROSIS;
double PULM_INF   = PULM_SCORE;
double IL5_obs    = IL5;

$CAPTURE
MEPO_Cobs BENRA_Cobs IMAT_Cobs PRED_Cobs
AEC_obs PERCENT_CHG RESP_300 RESP_1500
CARDIAC_SCORE PULM_INF IL5_obs
EoP EoI EoM_BM FIBROSIS PULM_SCORE
'

# ─────────────────────────────────────────────────────────────────────────────
# Compile model
# ─────────────────────────────────────────────────────────────────────────────
mod <- mcode("HES_QSP", hes_model_code)

# ─────────────────────────────────────────────────────────────────────────────
# Dosing helper functions
# ─────────────────────────────────────────────────────────────────────────────
# Mepolizumab 300 mg SC q4w (for HES; 100 mg q4w for severe asthma)
make_mepo_ev <- function(n_doses = 12, interval_wk = 4) {
  ev(cmt = "MEPO_DEPOT",
     amt = 300 * 1000,  # µg (300 mg)
     ii  = interval_wk * 7 * 24,
     addl = n_doses - 1,
     time = 0)
}

# Benralizumab 30 mg SC q4w x3 then q8w
make_benra_ev <- function(n_loading = 3, n_maint = 6) {
  loading <- ev(cmt = "BENRA_DEPOT",
                amt = 30 * 1000,  # µg
                ii  = 4 * 7 * 24,
                addl = n_loading - 1,
                time = 0)
  maint_start <- n_loading * 4 * 7 * 24
  maint <- ev(cmt = "BENRA_DEPOT",
              amt = 30 * 1000,
              ii  = 8 * 7 * 24,
              addl = n_maint - 1,
              time = maint_start)
  c(loading, maint)
}

# Imatinib 100 mg PO QD (low dose for HES)
make_imat_ev <- function(dose_mg = 100, duration_wk = 52) {
  ev(cmt = "IMAT_GUT",
     amt = dose_mg,
     ii  = 24,
     addl = duration_wk * 7 - 1,
     time = 0)
}

# Prednisolone 1 mg/kg/day PO (assume 70 kg → 70 mg/day) tapered
make_pred_ev <- function(dose_mg = 70, taper_wk = 8, maint_mg = 10, maint_wk = 44) {
  induction <- ev(cmt = "PRED_GUT",
                  amt = dose_mg,
                  ii  = 24,
                  addl = taper_wk * 7 - 1,
                  time = 0)
  maint_start <- taper_wk * 7 * 24
  maint <- ev(cmt = "PRED_GUT",
              amt = maint_mg,
              ii  = 24,
              addl = maint_wk * 7 - 1,
              time = maint_start)
  c(induction, maint)
}

# ─────────────────────────────────────────────────────────────────────────────
# 5 TREATMENT SCENARIOS
# ─────────────────────────────────────────────────────────────────────────────
sim_end <- 52 * 7 * 24  # 52 weeks in hours
delta   <- 24           # daily output

#  Scenario 1: Natural history (no treatment) ────────────────────────────────
out_untreated <- mod %>%
  param(AEC0 = 3000, USE_MEPO = 0, USE_BENRA = 0,
        USE_IMAT = 0, USE_PRED = 0, CLONAL_HES = 0) %>%
  mrgsim(end = sim_end, delta = delta) %>%
  as.data.frame() %>%
  mutate(scenario = "Untreated (Reactive HES)")

#  Scenario 2: Prednisolone monotherapy ──────────────────────────────────────
pred_dose <- make_pred_ev(dose_mg = 70, taper_wk = 8, maint_mg = 10, maint_wk = 44)
out_pred <- mod %>%
  param(AEC0 = 3000, USE_MEPO = 0, USE_BENRA = 0,
        USE_IMAT = 0, USE_PRED = 1, CLONAL_HES = 0) %>%
  mrgsim(events = pred_dose, end = sim_end, delta = delta) %>%
  as.data.frame() %>%
  mutate(scenario = "Prednisolone (1 mg/kg → taper)")

#  Scenario 3: Mepolizumab 300 mg q4w ────────────────────────────────────────
mepo_dose <- make_mepo_ev(n_doses = 13, interval_wk = 4)
out_mepo <- mod %>%
  param(AEC0 = 3000, USE_MEPO = 1, USE_BENRA = 0,
        USE_IMAT = 0, USE_PRED = 0, CLONAL_HES = 0) %>%
  mrgsim(events = mepo_dose, end = sim_end, delta = delta) %>%
  as.data.frame() %>%
  mutate(scenario = "Mepolizumab 300 mg q4w")

#  Scenario 4: Benralizumab 30 mg q4w → q8w ──────────────────────────────────
benra_dose <- make_benra_ev(n_loading = 3, n_maint = 6)
out_benra <- mod %>%
  param(AEC0 = 3000, USE_MEPO = 0, USE_BENRA = 1,
        USE_IMAT = 0, USE_PRED = 0, CLONAL_HES = 0) %>%
  mrgsim(events = benra_dose, end = sim_end, delta = delta) %>%
  as.data.frame() %>%
  mutate(scenario = "Benralizumab 30 mg q4w→q8w")

#  Scenario 5: Imatinib for clonal (FIP1L1-PDGFRA+) HES ──────────────────────
imat_dose <- make_imat_ev(dose_mg = 100, duration_wk = 52)
out_imat <- mod %>%
  param(AEC0 = 5000, USE_MEPO = 0, USE_BENRA = 0,
        USE_IMAT = 1, USE_PRED = 0,
        CLONAL_HES = 1, CLONAL_FOLD = 5.0,
        EO_BLOOD = 5000, EoP = 200) %>%
  mrgsim(events = imat_dose, end = sim_end, delta = delta) %>%
  as.data.frame() %>%
  mutate(scenario = "Imatinib 100mg/day (Clonal HES, FIP1L1-PDGFRA+)")

# Combine all scenarios
all_scenarios <- bind_rows(
  out_untreated, out_pred, out_mepo, out_benra, out_imat
) %>%
  mutate(
    time_wk = time / (7 * 24),
    scenario = factor(scenario, levels = c(
      "Untreated (Reactive HES)",
      "Prednisolone (1 mg/kg → taper)",
      "Mepolizumab 300 mg q4w",
      "Benralizumab 30 mg q4w→q8w",
      "Imatinib 100mg/day (Clonal HES, FIP1L1-PDGFRA+)"
    ))
  )

scenario_colors <- c(
  "Untreated (Reactive HES)"              = "#E74C3C",
  "Prednisolone (1 mg/kg → taper)"        = "#F39C12",
  "Mepolizumab 300 mg q4w"                = "#2ECC71",
  "Benralizumab 30 mg q4w→q8w"            = "#3498DB",
  "Imatinib 100mg/day (Clonal HES, FIP1L1-PDGFRA+)" = "#9B59B6"
)

# ─────────────────────────────────────────────────────────────────────────────
# PLOTS
# ─────────────────────────────────────────────────────────────────────────────
hes_theme <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13))

# Plot 1: Absolute eosinophil count over time
p_aec <- ggplot(all_scenarios, aes(x = time_wk, y = AEC_obs, color = scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1500, linetype = "dashed", color = "orange", linewidth = 0.8) +
  geom_hline(yintercept = 300,  linetype = "dotted", color = "green4",  linewidth = 0.8) +
  annotate("text", x = 50, y = 1600, label = "HES threshold (1500/µL)", size = 3, color = "orange") +
  annotate("text", x = 50, y = 380,  label = "Target AEC (300/µL)",     size = 3, color = "green4") +
  scale_color_manual(values = scenario_colors) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Absolute Eosinophil Count (AEC) — 5 Treatment Scenarios",
       x = "Time (weeks)", y = "AEC (cells/µL)", color = "Scenario") +
  hes_theme

# Plot 2: IL-5 serum concentration over time
p_il5 <- ggplot(all_scenarios, aes(x = time_wk, y = IL5_obs, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Serum IL-5 Dynamics",
       x = "Time (weeks)", y = "IL-5 (pg/mL)", color = "Scenario") +
  hes_theme

# Plot 3: Cardiac fibrosis progression
p_fibro <- ggplot(all_scenarios, aes(x = time_wk, y = CARDIAC_SCORE, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Cardiac Fibrosis Progression (Löffler/EMF)",
       x = "Time (weeks)", y = "Fibrosis Score (0–1)", color = "Scenario") +
  hes_theme

# Plot 4: Pulmonary infiltration
p_pulm <- ggplot(all_scenarios, aes(x = time_wk, y = PULM_INF, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Pulmonary Infiltration Score",
       x = "Time (weeks)", y = "Pulmonary Score (0–1)", color = "Scenario") +
  hes_theme

# Plot 5: Drug PK (mepolizumab & benralizumab)
mepo_data  <- all_scenarios %>% filter(scenario == "Mepolizumab 300 mg q4w")
benra_data <- all_scenarios %>% filter(scenario == "Benralizumab 30 mg q4w→q8w")

p_pk <- ggplot() +
  geom_line(data = mepo_data,  aes(x = time_wk, y = MEPO_Cobs,  color = "Mepolizumab (µg/mL)"), linewidth = 1.1) +
  geom_line(data = benra_data, aes(x = time_wk, y = BENRA_Cobs, color = "Benralizumab (µg/mL)"), linewidth = 1.1) +
  scale_color_manual(values = c("Mepolizumab (µg/mL)" = "#2ECC71",
                                "Benralizumab (µg/mL)" = "#3498DB")) +
  labs(title = "Biologic PK Profiles (Central Compartment)",
       x = "Time (weeks)", y = "Concentration (µg/mL)", color = "Drug") +
  hes_theme

# Plot 6: % Change from baseline AEC
p_pct <- ggplot(all_scenarios, aes(x = time_wk, y = PERCENT_CHG, color = scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = -50, linetype = "dashed", color = "gray40") +
  annotate("text", x = 49, y = -48, label = "−50% threshold", size = 3, color = "gray40") +
  scale_color_manual(values = scenario_colors) +
  labs(title = "% Change from Baseline AEC",
       x = "Time (weeks)", y = "% Change from Baseline", color = "Scenario") +
  hes_theme

# Composite
composite_plot <- (p_aec | p_il5) / (p_fibro | p_pulm) / (p_pk | p_pct) +
  plot_annotation(
    title = "Hypereosinophilic Syndrome — QSP Model Simulation Results",
    subtitle = "5 treatment scenarios: Untreated · Prednisolone · Mepolizumab · Benralizumab · Imatinib (Clonal HES)",
    theme = theme(plot.title    = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12))
  )

print(composite_plot)

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY TABLE (week 24 snapshot)
# ─────────────────────────────────────────────────────────────────────────────
wk24_summary <- all_scenarios %>%
  filter(abs(time_wk - 24) < 0.5) %>%
  group_by(scenario) %>%
  slice(1) %>%
  ungroup() %>%
  select(scenario, AEC_obs, PERCENT_CHG, RESP_300, RESP_1500,
         CARDIAC_SCORE, PULM_INF, IL5_obs) %>%
  rename(
    "Scenario"         = scenario,
    "AEC (cells/µL)"   = AEC_obs,
    "% Change AEC"     = PERCENT_CHG,
    "AEC<300 response" = RESP_300,
    "AEC<1500 control" = RESP_1500,
    "Cardiac Score"    = CARDIAC_SCORE,
    "Pulmonary Score"  = PULM_INF,
    "IL-5 (pg/mL)"     = IL5_obs
  )

cat("\n── Week 24 Summary Table ───────────────────────────────────────────────\n")
print(as.data.frame(wk24_summary), digits = 3)

# ─────────────────────────────────────────────────────────────────────────────
# SENSITIVITY ANALYSIS: AEC0 effect on fibrosis at week 52
# ─────────────────────────────────────────────────────────────────────────────
aec_range <- c(1500, 2000, 3000, 5000, 8000, 12000)

sens_data <- lapply(aec_range, function(aec_val) {
  out <- mod %>%
    param(AEC0 = aec_val, USE_MEPO = 0, USE_BENRA = 0,
          USE_IMAT = 0, USE_PRED = 0, CLONAL_HES = 0,
          EO_BLOOD = aec_val) %>%
    mrgsim(end = sim_end, delta = delta) %>%
    as.data.frame() %>%
    filter(abs(time - sim_end) < 25) %>%
    slice(1) %>%
    mutate(AEC_baseline = aec_val)
  out
}) %>% bind_rows()

p_sens <- ggplot(sens_data, aes(x = AEC_baseline, y = CARDIAC_SCORE)) +
  geom_line(color = "#E74C3C", linewidth = 1.2) +
  geom_point(color = "#E74C3C", size = 3) +
  scale_x_continuous(labels = scales::comma) +
  labs(title = "Sensitivity: Baseline AEC vs Cardiac Fibrosis at Week 52",
       x = "Baseline AEC (cells/µL)",
       y = "Cardiac Fibrosis Score at Wk 52") +
  hes_theme

print(p_sens)

cat("\n── Model run complete ───────────────────────────────────────────────────\n")
cat("  Compartments    : 20 ODEs\n")
cat("  Treatment scenarios : 5\n")
cat("  PK drugs modeled    : Mepolizumab (TMDD), Benralizumab (ADCC), Imatinib, Prednisolone\n")
cat("  Disease endpoints   : AEC, IL-5, Cardiac fibrosis, Pulmonary infiltration\n")
