## ============================================================
## CKD-Mineral Bone Disorder (CKD-MBD) — mrgsolve QSP Model
## FGF23 · Klotho · PTH · Vitamin D · Bone · Vascular Calcification
## 17 ODE Compartments · 7 Treatment Scenarios
## Parameter calibration refs: KDIGO 2017 CKD-MBD Guidelines,
##   Ix & Shlipak 2010 (FGF23-mortality), Tentori 2008 (DOPPS),
##   Chertow 2002 (ACHIEVE), Moe 2005 (ADVANCE)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ─── Model Code ───────────────────────────────────────────────────────────────
code <- '
$PARAM @annotated
// ── Patient / CKD baseline ──────────────────────────────────────
GFR0        : 15    : Baseline GFR (mL/min/1.73m2)  CKD G5
Pi_in       : 1200  : Dietary phosphate intake (mg/day)
Ca_in       : 800   : Dietary calcium intake (mg/day)
VitD_sun    : 0.1   : Cutaneous Vitamin D3 synthesis (nmol/day)

// ── Phosphate homeostasis ────────────────────────────────────────
kabs_Pi     : 0.60  : Fractional intestinal Pi absorption
kPi_urine   : 0.0008: Renal Pi clearance per GFR unit (L/min)
kPi_bone    : 0.005 : Pi release from bone (pool/day)
kPi_NaPi    : 0.15  : NaPi-IIa/IIc downregulation by FGF23
Pi_ss       : 4.5   : Target serum Pi (mg/dL) — normal

// ── FGF23 dynamics ──────────────────────────────────────────────
kFGF23_syn  : 0.02  : FGF23 synthesis rate (pmol/day)
kFGF23_deg  : 0.35  : FGF23 degradation rate (1/day)
EC50_Pi_FGF23: 5.5  : Pi EC50 for FGF23 stimulation (mg/dL)
Emax_Pi_FGF23: 8.0  : Maximum Pi-driven FGF23 fold increase
EC50_VitD_FGF23: 50 : VitD EC50 for FGF23 (pg/mL)
FGF23_ss    : 30    : Normal FGF23 (pg/mL)

// ── Klotho dynamics ─────────────────────────────────────────────
Klotho_ss   : 1.0   : Normal Klotho (relative units)
kKlotho_deg : 0.12  : Klotho degradation (1/day), ↑ in CKD
kKlotho_GFR : 0.04  : GFR-dependent Klotho synthesis

// ── PTH dynamics (parathyroid gland) ────────────────────────────
kPTH_syn    : 0.50  : PTH synthesis rate baseline (pmol/day)
kPTH_deg    : 1.20  : PTH degradation (1/day)
EC50_Ca_PTH : 1.10  : Ionized Ca2+ EC50 for PTH suppression (mmol/L)
IC50_VitD_PTH: 30   : VitD IC50 for PTH suppression (pg/mL)
IC50_CaSR_cin: 15   : Cinacalcet IC50 for CaSR activation (ng/mL)
Hill_PTH    : 3.5   : Hill coefficient (CaSR cooperativity)
PTH_ss      : 65    : Normal intact PTH (pg/mL)

// ── Vitamin D axis ───────────────────────────────────────────────
k25_syn     : 0.08  : 25-OH-D synthesis rate (nmol/day)
k25_deg     : 0.02  : 25-OH-D catabolism (1/day)
k125_syn    : 0.015 : 1,25-(OH)2D synthesis (CYP27B1, 1/day)
k125_deg    : 0.30  : 1,25-(OH)2D catabolism (CYP24A1, 1/day)
Imax_FGF23_CYP27B1: 0.80 : FGF23 max inhibition of CYP27B1
IC50_FGF23_CYP27B1: 200  : FGF23 IC50 for CYP27B1 inhibition (pg/mL)
VitD_ss     : 30    : Normal 25-OH-D (ng/mL → 75 nmol/L)
VitD_active_ss: 40  : Normal 1,25-OH2D (pg/mL)

// ── Calcium homeostasis ──────────────────────────────────────────
kabs_Ca     : 0.30  : Fractional intestinal Ca absorption
kCa_urine   : 0.006 : Renal Ca excretion rate constant
kCa_bone    : 0.008 : Net bone Ca exchange rate
Ca_ss       : 9.5   : Normal serum Ca (mg/dL)
Ca_dial     : 0     : Dialysate Ca contribution (0 = no dialysis)

// ── Bone remodeling ──────────────────────────────────────────────
kOB_syn     : 0.008 : Osteoblast formation rate
kOB_deg     : 0.02  : Osteoblast apoptosis (1/day)
kOC_syn     : 0.006 : Osteoclast formation rate (RANKL-driven)
kOC_deg     : 0.03  : Osteoclast apoptosis (1/day)
kRANKL_OPG  : 0.5   : OPG inhibition of RANKL (relative)
kBMD_form   : 0.0003: BMD formation rate (by OB)
kBMD_res    : 0.0004: BMD resorption rate (by OC)
BMD_ss      : 1.0   : Normal BMD (relative)

// ── Vascular calcification ───────────────────────────────────────
kVASC_calc  : 0.0001: Vascular Ca-Pi deposition rate
kMGP_syn    : 0.05  : MGP synthesis rate (carboxylated, protective)
kFetuin_prot: 0.08  : Fetuin-A protective factor
VC_threshold: 55    : Ca×Pi product threshold for calcification

// ── Drug PK parameters ──────────────────────────────────────────
// Sevelamer (phosphate binder, non-absorbed)
Sev_Emax    : 0.65  : Max Pi binding by sevelamer (fraction)
Sev_EC50    : 800   : ED50 dose (mg/day)
// Cinacalcet (calcimimetic, oral)
Cin_F       : 0.21  : Bioavailability
Cin_ka      : 1.2   : Absorption rate (1/hr)
Cin_V       : 1000  : Volume of distribution (L)
Cin_CL      : 250   : Clearance (L/hr)
// Paricalcitol (VDRa, IV for dialysis)
Par_F       : 1.0   : Bioavailability (IV)
Par_ka      : 0     : N/A (IV bolus)
Par_V       : 34    : Volume of distribution (L)
Par_CL      : 17    : Clearance (L/hr)
Par_Emax    : 0.90  : Max VDR activation (fraction)
Par_EC50    : 0.2   : EC50 (ng/mL)
// Etelcalcetide (IV calcimimetic)
Etel_ka     : 0     : N/A (IV)
Etel_V      : 11.5  : Volume of distribution (L)
Etel_CL     : 0.60  : Clearance (L/hr)
Etel_Emax   : 0.85  : Max CaSR activation
Etel_EC50   : 50    : EC50 (ng/mL)
// Denosumab (RANKL inhibitor, SC)
Den_F       : 0.62  : Bioavailability
Den_ka      : 0.006 : Absorption rate (1/hr)
Den_V       : 2.8   : Volume of distribution (L)
Den_CL      : 0.008 : Clearance (L/hr)
Den_Imax    : 0.95  : Max RANKL inhibition

// ── Dosing flags ─────────────────────────────────────────────────
DOSE_Sev    : 0     : Sevelamer daily dose (mg/day) — 0 = off
DOSE_Cin    : 0     : Cinacalcet daily dose (mg/day)
DOSE_Par    : 0     : Paricalcitol dose per session (mcg)
DOSE_Etel   : 0     : Etelcalcetide dose per session (mg)
DOSE_Den    : 0     : Denosumab dose (mg, SC q6mo)
DOSE_CaCO3  : 0     : Calcium carbonate dose (mg elemental Ca/day)

$CMT @annotated
Pi          : Serum phosphate (mg/dL)
FGF23       : Plasma FGF23 (pg/mL)
Klotho      : Soluble Klotho (rel. units)
PTH         : Intact PTH (pg/mL)
VitD25      : 25-OH-Vitamin D (nmol/L)
VitD_act    : 1,25-(OH)2D — calcitriol (pg/mL)
Ca          : Serum total calcium (mg/dL)
OB          : Osteoblast activity (rel. units)
OC          : Osteoclast activity (rel. units)
BMD         : Bone mineral density (rel. to normal)
VC          : Vascular calcification score (AU)
// Drug PK compartments
CIN_GUT     : Cinacalcet gut (mg)
CIN_PLASMA  : Cinacalcet plasma (mg/L)
PAR_PLASMA  : Paricalcitol plasma (ng/mL)
ETEL_PLASMA : Etelcalcetide plasma (ng/mL)
DEN_DEPOT   : Denosumab SC depot (mg)
DEN_PLASMA  : Denosumab plasma (mg/L)

$MAIN
// CaSR occupancy: combined cinacalcet + etelcalcetide
double Ca_ion  = Ca / 2.51;   // total → ionized Ca2+ (mmol/L)
double cin_eff = CIN_PLASMA;
double etel_eff = ETEL_PLASMA;

// CaSR activation Hill function
double CaSR_act = pow(Ca_ion, Hill_PTH) /
    (pow(EC50_Ca_PTH, Hill_PTH) + pow(Ca_ion, Hill_PTH));

// Calcimimetic boost on CaSR
double Cin_CaSR = cin_eff / (IC50_CaSR_cin + cin_eff);
double Etel_CaSR = etel_eff / (Etel_EC50 + etel_eff);
double CaSR_total = CaSR_act + (1 - CaSR_act) * (Cin_CaSR * 0.7 + Etel_CaSR * Etel_Emax);

// VDR activation: endogenous + paricalcitol
double VitD_VDR = VitD_act / (VitD_act + IC50_VitD_PTH);
double Par_VDR  = (PAR_PLASMA * Par_Emax) / (PAR_PLASMA + Par_EC50);
double VDR_act  = fmin(1.0, VitD_VDR + Par_VDR);

// Sevelamer Pi binding (GI)
double Sev_bind = (DOSE_Sev * Sev_Emax) / (DOSE_Sev + Sev_EC50);

// CaCO3 extra Ca absorption
double CaCO3_Ca = DOSE_CaCO3 * 0.40 / 10.0;   // elemental→ mmol→ mg/dL equivalent

// FGF23 stimulation by Pi
double Pi_stim_FGF23 = 1.0 + Emax_Pi_FGF23 * pow(fmax(Pi - 3.5, 0), 2) /
    (pow(EC50_Pi_FGF23 - 3.5, 2) + pow(fmax(Pi - 3.5, 0), 2));

// FGF23 stimulation by active VitD
double VitD_stim_FGF23 = 1.0 + (VitD_act / (EC50_VitD_FGF23 + VitD_act));

// GFR-dependent renal Pi clearance
double GFR_act = GFR0 * (1.0 - 0.005 * (VC + 0.001));  // vascular damage slowly ↓ GFR
double kPi_clear = kPi_urine * GFR_act;

// CYP27B1 inhibition by FGF23 and hyperparathyroidism (PTH stimulates)
double FGF23_inh_CYP27B1 = 1.0 - Imax_FGF23_CYP27B1 * FGF23 /
    (IC50_FGF23_CYP27B1 + FGF23);
double PTH_stim_CYP27B1 = 1.0 + 0.5 * PTH / (PTH + 100);

// RANKL/OPG driven osteoclast formation
double Den_RANKL_inh = (DEN_PLASMA * Den_Imax) / (DEN_PLASMA + 0.03);
double RANKL_eff = kOC_syn * OB * (1.0 - kRANKL_OPG * VDR_act) *
    (1.0 - Den_RANKL_inh);

// Ca×Pi product for vascular calcification
double CaP_prod = Ca * Pi;

$ODE
// ── [1] Serum Phosphate ──────────────────────────────────────────
dxdt_Pi = (kabs_Pi * (1.0 - Sev_bind) * Pi_in / 100.0)  // dietary absorption
        + kPi_bone * BMD                                  // bone release
        - kPi_clear * Pi                                  // renal excretion
        - 0.002 * Pi;                                     // other losses

// ── [2] FGF23 ────────────────────────────────────────────────────
dxdt_FGF23 = kFGF23_syn * Pi_stim_FGF23 * VitD_stim_FGF23
           - kFGF23_deg * FGF23;

// ── [3] Klotho ───────────────────────────────────────────────────
dxdt_Klotho = kKlotho_GFR * GFR_act
            - kKlotho_deg * Klotho;

// ── [4] PTH ──────────────────────────────────────────────────────
dxdt_PTH = kPTH_syn * (1.0 - CaSR_total) * (1.0 - VDR_act * 0.6)
         + 0.02 * Pi_stim_FGF23   // hyperphosphatemia directly stimulates
         - kPTH_deg * PTH;

// ── [5] 25-OH-Vitamin D ──────────────────────────────────────────
dxdt_VitD25 = k25_syn + VitD_sun
            - k25_deg * VitD25;

// ── [6] 1,25-(OH)2D (Calcitriol) ────────────────────────────────
dxdt_VitD_act = k125_syn * VitD25 * FGF23_inh_CYP27B1 * PTH_stim_CYP27B1
              + PAR_PLASMA * Par_Emax * 0.1   // paricalcitol contribution
              - k125_deg * VitD_act
              - 0.05 * VitD_act * VitD_act / (50 + VitD_act);  // CYP24A1 autoinduction

// ── [7] Serum Calcium ─────────────────────────────────────────────
dxdt_Ca = (kabs_Ca * (1.0 + 0.5 * VDR_act) * Ca_in / 200.0)   // intestinal
         + CaCO3_Ca                                              // supplement
         + Ca_dial                                               // dialysate
         + 0.003 * PTH * BMD                                     // PTH→bone resorption
         - kCa_urine * GFR_act * Ca                              // renal excretion
         - kCa_bone * (OB - OC * 0.7);                          // bone formation

// ── [8] Osteoblast ───────────────────────────────────────────────
dxdt_OB = kOB_syn * (1.0 + 0.4 * VDR_act) * (1.0 - 0.3 * Sclerostin_eff)
        - kOB_deg * OB;
// Sclerostin inhibition (elevated in CKD): approximated
double Sclerostin_eff = fmin(1.0, 0.3 + 0.5 * (1.0 - GFR_act / 90.0));

// ── [9] Osteoclast ───────────────────────────────────────────────
dxdt_OC = RANKL_eff * (1.0 + 0.6 * PTH / (PTH + 150))
        - kOC_deg * OC;

// ── [10] BMD ─────────────────────────────────────────────────────
dxdt_BMD = kBMD_form * OB - kBMD_res * OC - 0.0001 * (PTH / 65 - 1) * BMD;

// ── [11] Vascular Calcification ──────────────────────────────────
double VC_drive = fmax(0, CaP_prod - VC_threshold);
dxdt_VC = kVASC_calc * VC_drive * (1.0 - kMGP_syn) * (1.0 - kFetuin_prot * 0.5)
        + 0.00002 * Pi * (Pi - 3.5)   // direct VSMC osteogenic switch
        - 0.0002 * VC;                 // regression (slow)

// ── Drug PK ODEs ────────────────────────────────────────────────
// [12] Cinacalcet GI absorption
dxdt_CIN_GUT    = -Cin_ka * CIN_GUT;
// [13] Cinacalcet plasma
dxdt_CIN_PLASMA = Cin_ka * CIN_GUT * Cin_F / Cin_V - Cin_CL / Cin_V * CIN_PLASMA;
// [14] Paricalcitol plasma (IV bolus — handled via event)
dxdt_PAR_PLASMA = -Par_CL / Par_V * PAR_PLASMA;
// [15] Etelcalcetide plasma (IV bolus)
dxdt_ETEL_PLASMA = -Etel_CL / Etel_V * ETEL_PLASMA;
// [16] Denosumab SC depot
dxdt_DEN_DEPOT  = -Den_ka * DEN_DEPOT;
// [17] Denosumab plasma
dxdt_DEN_PLASMA = Den_ka * DEN_DEPOT * Den_F / Den_V - Den_CL / Den_V * DEN_PLASMA;

$TABLE
double iPTH     = PTH;
double sPi      = Pi;
double sCa      = Ca;
double s25D     = VitD25;
double s125D    = VitD_act;
double sFGF23   = FGF23;
double sKlotho  = Klotho;
double sBMD     = BMD;
double sVC      = VC;
double CaP      = Ca * Pi;
double OB_act   = OB;
double OC_act   = OC;
double sGFR     = GFR_act;

$CAPTURE iPTH sPi sCa s25D s125D sFGF23 sKlotho sBMD sVC CaP OB_act OC_act sGFR
'

## ─── Compile ──────────────────────────────────────────────────────────────────
mod <- mcode("ckdmbd_qsp", code)

## ─── Initial Conditions (CKD G4-G5, pre-dialysis) ────────────────────────────
init_ckd <- init(mod,
  Pi          = 6.2,    # hyperphosphatemia
  FGF23       = 800,    # markedly elevated FGF23
  Klotho      = 0.35,   # severely reduced Klotho
  PTH         = 420,    # secondary hyperparathyroidism
  VitD25      = 18,     # 25-OH-D insufficient
  VitD_act    = 15,     # calcitriol deficiency
  Ca          = 8.8,    # low-normal Ca
  OB          = 0.7,    # reduced osteoblast activity
  OC          = 1.4,    # increased osteoclast activity
  BMD         = 0.82,   # osteopenia
  VC          = 25,     # moderate vascular calcification
  CIN_GUT     = 0,
  CIN_PLASMA  = 0,
  PAR_PLASMA  = 0,
  ETEL_PLASMA = 0,
  DEN_DEPOT   = 0,
  DEN_PLASMA  = 0
)

## ─── Simulation Time (2 years, 730 days) ─────────────────────────────────────
tgrid <- tgrid(0, 730, 1)

## ─── Helper: Dosing Events ────────────────────────────────────────────────────
mk_cin_events <- function(dose_mg, days = seq(0, 729, 1)) {
  # Cinacalcet 30–90 mg once daily → gut bolus
  ev(cmt = "CIN_GUT", amt = dose_mg, ii = 24/24, addl = length(days) - 1, time = 0)
}
mk_par_events <- function(dose_mcg, freq_days = 2) {
  # Paricalcitol 2–4 mcg IV 3x/wk during HD
  ev(cmt = "PAR_PLASMA", amt = dose_mcg, ii = freq_days * 24, addl = floor(730 / freq_days), time = 0)
}
mk_etel_events <- function(dose_mg, freq_days = 2) {
  ev(cmt = "ETEL_PLASMA", amt = dose_mg * 1000, ii = freq_days * 24, addl = floor(730 / freq_days), time = 0)
}
mk_den_events <- function(dose_mg = 60, interval_days = 180) {
  # Denosumab 60 mg SC every 6 months
  ev(cmt = "DEN_DEPOT", amt = dose_mg, ii = interval_days * 24, addl = 3, time = 0)
}

## ─── 7 Treatment Scenarios ────────────────────────────────────────────────────

## Scenario 1: Untreated CKD G5 (natural progression)
sim_s1 <- mod %>%
  init(init_ckd) %>%
  param(DOSE_Sev = 0, DOSE_Cin = 0, DOSE_Par = 0) %>%
  mrgsim(tgrid = tgrid) %>%
  as_tibble() %>%
  mutate(Scenario = "S1: Untreated CKD G5")

## Scenario 2: Sevelamer 2400 mg/day (phosphate binder monotherapy)
sim_s2 <- mod %>%
  init(init_ckd) %>%
  param(DOSE_Sev = 2400) %>%
  mrgsim(tgrid = tgrid) %>%
  as_tibble() %>%
  mutate(Scenario = "S2: Sevelamer 2400 mg/d")

## Scenario 3: Cinacalcet 60 mg/day (calcimimetic)
e3 <- mk_cin_events(60)
sim_s3 <- mod %>%
  init(init_ckd) %>%
  mrgsim(events = e3, tgrid = tgrid) %>%
  as_tibble() %>%
  mutate(Scenario = "S3: Cinacalcet 60 mg/d")

## Scenario 4: Paricalcitol 4 mcg 3x/wk IV (VDR agonist)
e4 <- mk_par_events(4, freq_days = 2)
sim_s4 <- mod %>%
  init(init_ckd) %>%
  mrgsim(events = e4, tgrid = tgrid) %>%
  as_tibble() %>%
  mutate(Scenario = "S4: Paricalcitol 4 mcg 3x/wk")

## Scenario 5: Sevelamer + Cinacalcet (combination)
e5 <- mk_cin_events(60)
sim_s5 <- mod %>%
  init(init_ckd) %>%
  param(DOSE_Sev = 2400) %>%
  mrgsim(events = e5, tgrid = tgrid) %>%
  as_tibble() %>%
  mutate(Scenario = "S5: Sevelamer + Cinacalcet")

## Scenario 6: Etelcalcetide 5 mg IV 3x/wk (next-gen calcimimetic)
e6 <- mk_etel_events(5, freq_days = 2)
sim_s6 <- mod %>%
  init(init_ckd) %>%
  mrgsim(events = e6, tgrid = tgrid) %>%
  as_tibble() %>%
  mutate(Scenario = "S6: Etelcalcetide 5 mg 3x/wk")

## Scenario 7: Triple therapy — Sevelamer + Etelcalcetide + Denosumab
e7_etel <- mk_etel_events(5, freq_days = 2)
e7_den  <- mk_den_events(60, interval_days = 180)
e7 <- ev_seq(e7_etel, e7_den)   # stack events
sim_s7 <- mod %>%
  init(init_ckd) %>%
  param(DOSE_Sev = 2400) %>%
  mrgsim(events = e7, tgrid = tgrid) %>%
  as_tibble() %>%
  mutate(Scenario = "S7: Sevelamer + Etelcalcetide + Denosumab")

## ─── Combine Results ──────────────────────────────────────────────────────────
all_sims <- bind_rows(sim_s1, sim_s2, sim_s3, sim_s4, sim_s5, sim_s6, sim_s7)

## ─── Plotting ─────────────────────────────────────────────────────────────────
plot_var <- function(data, var, ylab, target_lo = NULL, target_hi = NULL) {
  p <- ggplot(data, aes(x = time, y = .data[[var]], color = Scenario)) +
    geom_line(linewidth = 0.8) +
    labs(x = "Day", y = ylab, title = ylab) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8))
  if (!is.null(target_lo)) p <- p + geom_hline(yintercept = target_lo, linetype = 2, color = "green4")
  if (!is.null(target_hi)) p <- p + geom_hline(yintercept = target_hi, linetype = 2, color = "red3")
  p
}

p_pth  <- plot_var(all_sims, "iPTH",   "iPTH (pg/mL)",       target_lo = 150, target_hi = 600)
p_pi   <- plot_var(all_sims, "sPi",    "Serum Pi (mg/dL)",   target_hi = 5.5)
p_ca   <- plot_var(all_sims, "sCa",    "Serum Ca (mg/dL)",   target_lo = 8.4, target_hi = 10.2)
p_fgf  <- plot_var(all_sims, "sFGF23", "FGF23 (pg/mL)")
p_vitd <- plot_var(all_sims, "s25D",   "25-OH-D (nmol/L)",   target_lo = 50)
p_bmd  <- plot_var(all_sims, "sBMD",   "BMD (rel. to normal)")
p_vc   <- plot_var(all_sims, "sVC",    "Vascular Calcification Score")
p_cap  <- plot_var(all_sims, "CaP",    "Ca×Pi Product",       target_hi = 55)

## ─── Summary Table at Day 180 and Day 365 ────────────────────────────────────
summary_tbl <- all_sims %>%
  filter(time %in% c(0, 90, 180, 365, 730)) %>%
  select(Scenario, time, iPTH, sPi, sCa, s25D, sFGF23, sBMD, sVC) %>%
  arrange(Scenario, time)

print(summary_tbl, n = 50)

## ─── Target Attainment (KDIGO 2017 Targets) ──────────────────────────────────
kdigo_check <- all_sims %>%
  filter(time == 365) %>%
  mutate(
    PTH_ok  = iPTH >= 150 & iPTH <= 600,
    Pi_ok   = sPi <= 5.5,
    Ca_ok   = sCa >= 8.4 & sCa <= 10.2,
    CaP_ok  = (sCa * sPi) < 55,
    VitD_ok = s25D >= 50
  ) %>%
  select(Scenario, iPTH, sPi, sCa, PTH_ok, Pi_ok, Ca_ok, CaP_ok, VitD_ok)

print("=== KDIGO Target Attainment at Day 365 ===")
print(kdigo_check)
