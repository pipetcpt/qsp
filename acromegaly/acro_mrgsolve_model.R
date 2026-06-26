## ============================================================
## Acromegaly QSP Model — mrgsolve Implementation
## GH-Secreting Pituitary Adenoma: PK/PD/Disease Model
##
## Disease: Acromegaly (Somatotropinoma)
## Drug targets: SSTR2/SSTR5 (SSA), GHR (pegvisomant), D2R (cabergoline)
##
## ODE States (17 compartments):
##   1.  DEPOT_SSA    — SSA microsphere/autogel depot (IM or SC)
##   2.  CENT_SSA     — SSA central plasma compartment (ng/mL)
##   3.  PERI_SSA     — SSA peripheral tissue compartment
##   4.  CENT_PEG     — Pegvisomant central compartment (ng/mL)
##   5.  PERI_PEG     — Pegvisomant tissue compartment
##   6.  GH_ADENOM    — GH production rate from adenoma (ng/mL/h driver)
##   7.  GH_PLASMA    — Plasma GH (ng/mL)
##   8.  SSTR_BOUND   — SSA-bound SSTR2/SSTR5 receptors (fraction)
##   9.  GHR_FREE     — Free GH receptors (normalized, 0–1)
##  10.  GHR_BLOCKED  — Pegvisomant-occupied GHR (fraction)
##  11.  STAT5b_ACT   — Active STAT5b signaling (normalized, 0–1)
##  12.  IGF1_LIVER   — Hepatic IGF-1 production rate (ng/mL/h)
##  13.  IGF1_PLASMA  — Plasma IGF-1 (ng/mL)
##  14.  ADENOM_VOL   — Pituitary adenoma volume (cm³)
##  15.  LVH_IDX      — LV mass index (g/m²; normal <115 M, <95 F)
##  16.  GLUCOSE      — Fasting plasma glucose (mg/dL)
##  17.  ARTH_SCORE   — Arthropathy cumulative damage score (0–100)
##
## Treatment scenarios:
##   1. Untreated acromegaly (baseline)
##   2. Successful surgery (70% GH reduction, instantaneous)
##   3. Octreotide LAR 30 mg IM q28d
##   4. Lanreotide Autogel 120 mg SC q28d
##   5. Pasireotide LAR 60 mg IM q28d
##   6. Pegvisomant 15 mg SC daily
##   7. SSA + pegvisomant combination
##   8. Failed surgery + SSA rescue
##
## Key references for parameters:
##   - Colao et al. NEJM 2009; 360:2467  (SSA primary therapy)
##   - Trainer et al. NEJM 2000; 342:1171 (pegvisomant)
##   - Chanson et al. Endocr Rev 2019 (review)
##   - Freda et al. JCEM 2005 (GH kinetics)
##   - Gatto et al. Pituitary 2015 (PK modeling)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ---- Model code -----------------------------------------------------
mod_code <- '
$PROB Acromegaly QSP Model — GH/IGF-1/Adenoma PK/PD

$PARAM
// ---- SSA PK parameters (Octreotide LAR reference) ----
KA_SSA   = 0.0060   // depot absorption rate (h-1); ~2-wk plateau
V1_SSA   = 14.0     // central volume (L)
V2_SSA   = 55.0     // peripheral volume (L)
Q_SSA    = 4.5      // intercompartmental CL (L/h)
CL_SSA   = 9.0      // systemic CL (L/h); biliary/hepatic
// Lanreotide AG modifier (same compartment, adjusted ka/CL)
KA_LAN   = 0.0040   // slower release (autogel nanotubes)
CL_LAN   = 7.5

// ---- Pasireotide LAR PK (multireceptor SSA) ----
KA_PAS   = 0.0055   // (h-1); similar to octreotide LAR
V1_PAS   = 100.0    // larger Vd (more tissue bound)
CL_PAS   = 5.0      // lower CL (longer t1/2 ~23 d)

// ---- Pegvisomant PK ----
KA_PEG   = 0.040    // SC absorption (h-1)
V1_PEG   = 7.0      // central (L)
V2_PEG   = 45.0     // tissue (L)
Q_PEG    = 2.5      // Q (L/h)
CL_PEG   = 0.45     // CL (L/h); slower than SSA (MW ~22 kDa)

// ---- SSTR2/5 receptor pharmacodynamics ----
IC50_SSTR2 = 0.4    // SSA plasma conc for 50% SSTR2 occupancy (ng/mL equiv ~0.2 nM)
IC50_SSTR5 = 8.0    // standard SSA (oct/lan); pasireotide IC50 much lower
IC50_PAS_SSTR5 = 0.25 // pasireotide IC50 for SSTR5 (ng/mL equiv)
IMAX_SSA   = 0.92   // maximum fractional GH suppression via SSTR

// ---- GHR block by pegvisomant ----
IC50_PEG_GHR = 5.0  // pegvisomant plasma conc for 50% GHR blockade (ng/mL)
IMAX_PEG     = 0.98 // max GHR blockade

// ---- Cabergoline (D2 agonist) fixed effect ----
CAB_EFFECT = 0.20   // max fractional GH reduction by cabergoline alone

// ---- GH secretion & kinetics ----
KOUT_GH    = 0.693  // GH elimination (h-1); t1/2 ~1 h plasma
GH_BASE    = 15.0   // basal adenoma GH production rate (ng/mL/h)
GH_PULSE   = 8.0    // extra GH per pulse (ng/mL added peak)
PULSE_FREQ = 0.33   // pulse frequency (pulses/h); ~1 per 3 h
// GH feedback from IGF-1
FB_IGF1    = 0.0010 // strength of IGF-1 neg feedback on GH production

// ---- STAT5b signaling ----
KACT_STAT  = 0.50   // STAT5b activation rate constant (h-1)
KDEG_STAT  = 0.35   // STAT5b deactivation/degrad (h-1)
SOCS2_FB   = 0.40   // SOCS2 negative feedback strength on STAT5b (0-1)

// ---- IGF-1 kinetics ----
KPROD_IGF1  = 1.20  // basal hepatic IGF-1 production constant (ng/mL/h per STAT5b unit)
KCLEAR_IGF1 = 0.010 // IGF-1 plasma clearance (h-1); t1/2 ~70 h (ternary complex)
IGF1_0      = 450.0 // typical active acromegaly IGF-1 at baseline (ng/mL)

// ---- Adenoma growth ----
KG_ADENOM  = 0.00005 // adenoma net growth rate (cm3/h); slow
KK_ADENOM  = 5.0     // carrying capacity (cm³)
ADENOM_0   = 1.5     // baseline macroadenoma volume (cm³)
// Surgery effect: rapid kill parameter
SURG_KILL  = 0.0     // instantaneous fraction killed by surgery (set in scenarios)

// ---- LVH dynamics ----
KIN_LVH    = 0.0020  // LVH induction rate per unit GH (g/m² per ng/mL/h)
KOUT_LVH   = 0.0004  // LVH regression rate (h-1); slow reversal
LVH_0      = 145.0   // baseline elevated LV mass index (g/m²)

// ---- Glucose dynamics ----
KINS_GH    = 0.080   // GH-induced insulin resistance → glucose increase (mg/dL per ng/mL/h)
KOUT_GLUC  = 0.0030  // glucose homeostasis return rate (h-1)
GLUC_0     = 105.0   // baseline fasting glucose in active acromegaly (mg/dL)
PAS_GLUC   = 0.0     // pasireotide-induced hyperglycemia additive effect (0 or 15)

// ---- Arthropathy ----
KARTH      = 0.000010 // IGF-1 driven arthropathy accumulation rate
KARTH_REM  = 0.000001 // minimal reversal rate (largely irreversible)

// ---- Flags ----
DRUG_TYPE  = 1   // 1=OctLAR, 2=LanAG, 3=PasLAR, 4=Pegvisomant, 5=Combo
USE_CAB    = 0   // 0=no cabergoline, 1=add cabergoline
SURG_DONE  = 0   // 0=no surgery, 1=surgery applied at t=0


$CMT DEPOT_SSA CENT_SSA PERI_SSA CENT_PEG PERI_PEG
     GH_ADENOM GH_PLASMA SSTR_BOUND GHR_FREE GHR_BLOCKED
     STAT5b_ACT IGF1_LIVER IGF1_PLASMA
     ADENOM_VOL LVH_IDX GLUCOSE ARTH_SCORE

$INIT
DEPOT_SSA  = 0
CENT_SSA   = 0
PERI_SSA   = 0
CENT_PEG   = 0
PERI_PEG   = 0
GH_ADENOM  = 15.0
GH_PLASMA  = 20.0
SSTR_BOUND = 0
GHR_FREE   = 1.0
GHR_BLOCKED = 0.0
STAT5b_ACT = 1.0
IGF1_LIVER = 5.0
IGF1_PLASMA = 450.0
ADENOM_VOL = 1.5
LVH_IDX    = 145.0
GLUCOSE    = 105.0
ARTH_SCORE = 10.0

$GLOBAL
// effective ka/CL based on drug type
double KA_eff, CL_eff;
// pulsatile GH signal (square pulse)
double PULSE_GH;

$MAIN
// Select SSA PK parameters by drug type
if(DRUG_TYPE == 2) {  // Lanreotide AG
  KA_eff = KA_LAN;
  CL_eff = CL_LAN;
} else if(DRUG_TYPE == 3) { // Pasireotide LAR
  KA_eff = KA_PAS;
  CL_eff = CL_PAS;
} else {  // Octreotide LAR (default) or combo SSA component
  KA_eff = KA_SSA;
  CL_eff = CL_SSA;
}

// Surgery: set adenoma volume fraction remaining
if(SURG_DONE == 1 && NEWIND <= 1) {
  ADENOM_VOL_0 = ADENOM_0 * (1.0 - SURG_KILL);
  GH_ADENOM_0  = GH_BASE  * (1.0 - SURG_KILL);
}

// Pulsatile GH: use sinusoidal approximation
double t_mod = fmod(TIME, 1.0 / PULSE_FREQ);
PULSE_GH = (t_mod < 0.5) ? GH_PULSE : 0.0;

$ODE
// ---- SSA Depot → Central → Peripheral PK ----
double KA_use = KA_eff;
double CL_use = CL_eff;

dxdt_DEPOT_SSA  = -KA_use * DEPOT_SSA;
dxdt_CENT_SSA   =  KA_use * DEPOT_SSA
                  - (CL_use / V1_SSA) * CENT_SSA
                  - (Q_SSA  / V1_SSA) * CENT_SSA
                  + (Q_SSA  / V2_SSA) * PERI_SSA;
dxdt_PERI_SSA   =  (Q_SSA  / V1_SSA) * CENT_SSA
                  - (Q_SSA  / V2_SSA) * PERI_SSA;

// ---- Pegvisomant PK ----
dxdt_CENT_PEG   =  KA_PEG * PERI_PEG  // note: SC given into PERI for simplicity
                  - (CL_PEG / V1_PEG) * CENT_PEG
                  - (Q_PEG  / V1_PEG) * CENT_PEG
                  + (Q_PEG  / V2_PEG) * PERI_PEG;
dxdt_PERI_PEG   = -(CL_PEG / V2_PEG) * PERI_PEG
                  - (Q_PEG  / V2_PEG) * PERI_PEG
                  + (Q_PEG  / V1_PEG) * CENT_PEG;

// ---- SSTR receptor occupancy (Emax model) ----
double SSTR2_occ, SSTR5_occ;
double IC50_s5 = (DRUG_TYPE == 3) ? IC50_PAS_SSTR5 : IC50_SSTR5;
SSTR2_occ = CENT_SSA / (CENT_SSA + IC50_SSTR2 + 1e-9);
SSTR5_occ = CENT_SSA / (CENT_SSA + IC50_s5 + 1e-9);
double SSTR_occ_total = (SSTR2_occ * 0.7 + SSTR5_occ * 0.3);  // SSTR2 dominant

dxdt_SSTR_BOUND = IMAX_SSA * SSTR_occ_total - SSTR_BOUND;  // fast equilibration

// ---- GHR blockade by pegvisomant ----
double PEG_GHR_occ = CENT_PEG / (CENT_PEG + IC50_PEG_GHR + 1e-9);
dxdt_GHR_FREE    = -IMAX_PEG * PEG_GHR_occ * GHR_FREE + 0.02 * GHR_BLOCKED;
dxdt_GHR_BLOCKED =  IMAX_PEG * PEG_GHR_occ * GHR_FREE - 0.02 * GHR_BLOCKED;

// ---- GH adenoma production ----
// Inhibition by SSA occupancy; negative feedback by IGF-1
double IGF1_FB = FB_IGF1 * IGF1_PLASMA;
double SSA_inhib = SSTR_BOUND;  // already scaled 0-1
double CAB_inhib = USE_CAB * CAB_EFFECT;
double GH_prod_rate = GH_BASE * (1.0 - SSA_inhib) * (1.0 - CAB_inhib)
                      / (1.0 + IGF1_FB);
// Proportional to adenoma volume fraction
double adenom_frac = ADENOM_VOL / ADENOM_0;
GH_prod_rate = GH_prod_rate * adenom_frac;
dxdt_GH_ADENOM = GH_prod_rate - KOUT_GH * GH_ADENOM;

// ---- Plasma GH (pulsatile secretion + clearance) ----
dxdt_GH_PLASMA = GH_ADENOM + PULSE_GH * (1.0 - SSA_inhib) * adenom_frac
                 - KOUT_GH * GH_PLASMA;

// ---- STAT5b activation (driven by GHR-available fraction) ----
double GHR_avail = GHR_FREE * (1.0 - GHR_BLOCKED);
double STAT5_input = KACT_STAT * GH_PLASMA * GHR_avail;
double STAT5_deg   = KDEG_STAT * (1.0 + SOCS2_FB * STAT5b_ACT) * STAT5b_ACT;
dxdt_STAT5b_ACT = STAT5_input - STAT5_deg;

// ---- Hepatic IGF-1 production ----
dxdt_IGF1_LIVER = KPROD_IGF1 * STAT5b_ACT - 0.50 * IGF1_LIVER;

// ---- Plasma IGF-1 dynamics ----
dxdt_IGF1_PLASMA = 50.0 * IGF1_LIVER - KCLEAR_IGF1 * IGF1_PLASMA;

// ---- Adenoma growth (logistic model with therapy) ----
double SURG_suppress = SURG_DONE * 0.0;  // handled in MAIN by volume reset
double SSA_antiprof  = 0.15 * SSTR_BOUND;   // SSA antiproliferative (weak)
double SRS_kill_rate = 0.0;                  // set nonzero for RT scenario
dxdt_ADENOM_VOL = KG_ADENOM * ADENOM_VOL
                  * (1.0 - ADENOM_VOL / KK_ADENOM)
                  * (1.0 - SSA_antiprof);

// ---- LVH dynamics ----
// GH drives LV mass directly; reverses slowly with treatment
dxdt_LVH_IDX = KIN_LVH * GH_PLASMA - KOUT_LVH * (LVH_IDX - 90.0);

// ---- Glucose homeostasis ----
// GH → insulin resistance → hyperglycemia
double PAS_glucose_add = (DRUG_TYPE == 3) ? PAS_GLUC : 0.0;
dxdt_GLUCOSE = KINS_GH * GH_PLASMA
              + PAS_glucose_add * 0.001
              - KOUT_GLUC * (GLUCOSE - 85.0);

// ---- Arthropathy accumulation ----
// Driven by IGF-1; partially reversible
dxdt_ARTH_SCORE = KARTH * IGF1_PLASMA - KARTH_REM * ARTH_SCORE;

$TABLE
// Capture output variables
double GH_mean   = GH_PLASMA;
double IGF1_out  = IGF1_PLASMA;
double SSTR_occ  = SSTR_BOUND;
double GHR_blk   = GHR_BLOCKED;
double SSA_plasma = CENT_SSA;
double PEG_plasma = CENT_PEG;
double LVmass     = LVH_IDX;
double Glucose_out= GLUCOSE;
double ArthrScore = ARTH_SCORE;
double TumorVol   = ADENOM_VOL;
// Biochemical remission flag (GH <2.5 ng/mL random, or IGF-1 normal)
double GH_ctrl    = (GH_PLASMA < 2.5) ? 1.0 : 0.0;
double IGF1_ctrl  = (IGF1_PLASMA < 250.0) ? 1.0 : 0.0;  // simplified age-norm
double Remission  = GH_ctrl * IGF1_ctrl;

$CAPTURE GH_mean IGF1_out SSTR_occ GHR_blk SSA_plasma PEG_plasma
         LVmass Glucose_out ArthrScore TumorVol Remission
         GH_ctrl IGF1_ctrl
'

## ---- Compile model --------------------------------------------------
mod <- mcode("acromegaly_qsp", mod_code)

## ====================================================================
## DOSING EVENTS — 8 Treatment Scenarios
## ====================================================================

## Time vector: 2 years (every 6 hours)
TSIM <- seq(0, 8760, by = 6)  # 8760 h = 1 year; extend to 2yr

## Helper: monthly dosing (28-day intervals × 13 doses = ~1 year)
monthly_doses <- seq(0, 8760, by = 672)  # 672 h = 28 days

## --------------------
## Scenario 1: Untreated (No treatment — natural history)
## --------------------
ev_s1 <- ev(time=0, amt=0, cmt=1)  # null dose
dat_s1 <- mod %>%
  param(DRUG_TYPE=1, SURG_DONE=0, USE_CAB=0, PAS_GLUC=0) %>%
  ev(ev_s1) %>%
  mrgsim(end=8760, delta=6) %>%
  as.data.frame() %>%
  mutate(scenario="S1: Untreated")

## --------------------
## Scenario 2: Successful Surgery (70% GH reduction at t=0)
## --------------------
ev_surg <- ev(time=0, amt=0, cmt=1)
dat_s2 <- mod %>%
  param(DRUG_TYPE=1, SURG_DONE=1, SURG_KILL=0.70, USE_CAB=0,
        GH_BASE=4.5, ADENOM_0=1.5) %>%
  init(GH_ADENOM=4.5, GH_PLASMA=6.0, ADENOM_VOL=0.45, IGF1_PLASMA=350.0) %>%
  ev(ev_surg) %>%
  mrgsim(end=8760, delta=6) %>%
  as.data.frame() %>%
  mutate(scenario="S2: Surgery (curative 70%)")

## --------------------
## Scenario 3: Octreotide LAR 30 mg IM q28d (DRUG_TYPE=1)
## --------------------
ev_oct <- lapply(monthly_doses, function(t) ev(time=t, amt=30, cmt=1, ii=0, addl=0))
ev_oct <- do.call(c, ev_oct)
dat_s3 <- mod %>%
  param(DRUG_TYPE=1, SURG_DONE=0, USE_CAB=0, PAS_GLUC=0) %>%
  ev(ev_oct) %>%
  mrgsim(tgrid=TSIM) %>%
  as.data.frame() %>%
  mutate(scenario="S3: Octreotide LAR 30 mg")

## --------------------
## Scenario 4: Lanreotide AG 120 mg SC q28d (DRUG_TYPE=2)
## --------------------
ev_lan <- lapply(monthly_doses, function(t) ev(time=t, amt=120, cmt=1, ii=0, addl=0))
ev_lan <- do.call(c, ev_lan)
dat_s4 <- mod %>%
  param(DRUG_TYPE=2, SURG_DONE=0, USE_CAB=0, PAS_GLUC=0,
        KA_SSA=0.0040, CL_SSA=7.5) %>%
  ev(ev_lan) %>%
  mrgsim(tgrid=TSIM) %>%
  as.data.frame() %>%
  mutate(scenario="S4: Lanreotide AG 120 mg")

## --------------------
## Scenario 5: Pasireotide LAR 60 mg IM q28d (DRUG_TYPE=3)
## --------------------
ev_pas <- lapply(monthly_doses, function(t) ev(time=t, amt=60, cmt=1, ii=0, addl=0))
ev_pas <- do.call(c, ev_pas)
dat_s5 <- mod %>%
  param(DRUG_TYPE=3, SURG_DONE=0, USE_CAB=0, PAS_GLUC=15.0,
        KA_SSA=0.0055, CL_SSA=5.0, V1_SSA=100.0,
        IC50_SSTR5=0.25) %>%
  ev(ev_pas) %>%
  mrgsim(tgrid=TSIM) %>%
  as.data.frame() %>%
  mutate(scenario="S5: Pasireotide LAR 60 mg")

## --------------------
## Scenario 6: Pegvisomant 15 mg SC daily (DRUG_TYPE=4)
## --------------------
daily_doses <- seq(0, 8760, by = 24)
ev_peg <- lapply(daily_doses, function(t) ev(time=t, amt=15, cmt=5, ii=0, addl=0))
ev_peg <- do.call(c, ev_peg)
dat_s6 <- mod %>%
  param(DRUG_TYPE=4, SURG_DONE=0, USE_CAB=0, PAS_GLUC=0) %>%
  ev(ev_peg) %>%
  mrgsim(tgrid=TSIM) %>%
  as.data.frame() %>%
  mutate(scenario="S6: Pegvisomant 15 mg/d")

## --------------------
## Scenario 7: Combination — Octreotide LAR + Pegvisomant (DRUG_TYPE=5)
## --------------------
ev_combo <- c(
  lapply(monthly_doses, function(t) ev(time=t, amt=20, cmt=1, ii=0, addl=0)),
  lapply(daily_doses,   function(t) ev(time=t, amt=10, cmt=5, ii=0, addl=0))
)
ev_combo <- do.call(c, ev_combo)
dat_s7 <- mod %>%
  param(DRUG_TYPE=5, SURG_DONE=0, USE_CAB=0, PAS_GLUC=0) %>%
  ev(ev_combo) %>%
  mrgsim(tgrid=TSIM) %>%
  as.data.frame() %>%
  mutate(scenario="S7: Combo Oct LAR+Pegvisomant")

## --------------------
## Scenario 8: Failed surgery + SSA rescue (post-op residual disease)
## --------------------
ev_postsurg_ssa <- c(
  ev(time=0, amt=0, cmt=1),  # surgery applied via init
  lapply(monthly_doses[monthly_doses > 730], function(t)  # start SSA at 1 month
    ev(time=t, amt=30, cmt=1, ii=0, addl=0))
)
ev_postsurg_ssa <- do.call(c, ev_postsurg_ssa)
dat_s8 <- mod %>%
  param(DRUG_TYPE=1, SURG_DONE=1, SURG_KILL=0.40, USE_CAB=0, PAS_GLUC=0,
        GH_BASE=9.0) %>
  init(GH_ADENOM=9.0, GH_PLASMA=12.0, ADENOM_VOL=0.9, IGF1_PLASMA=400.0) %>%
  ev(ev_postsurg_ssa) %>%
  mrgsim(tgrid=TSIM) %>%
  as.data.frame() %>%
  mutate(scenario="S8: Surgery+SSA rescue")

## ====================================================================
## COMBINE ALL SCENARIOS
## ====================================================================
all_scenarios <- bind_rows(dat_s1, dat_s2, dat_s3, dat_s4,
                            dat_s5, dat_s6, dat_s7, dat_s8)

## ====================================================================
## DOSE-RESPONSE ANALYSIS — Octreotide LAR doses 10/20/30 mg
## ====================================================================
doses_dr <- c(10, 20, 30)
dr_results <- lapply(doses_dr, function(d) {
  ev_dr <- lapply(monthly_doses, function(t) ev(time=t, amt=d, cmt=1, ii=0, addl=0))
  ev_dr <- do.call(c, ev_dr)
  mod %>%
    param(DRUG_TYPE=1) %>%
    ev(ev_dr) %>%
    mrgsim(tgrid=TSIM) %>%
    as.data.frame() %>%
    mutate(dose=d, label=paste0("OctLAR ", d, "mg"))
})
dr_data <- bind_rows(dr_results)

## ====================================================================
## SUMMARY TABLE
## ====================================================================
summary_tbl <- all_scenarios %>%
  filter(time %in% c(0, 2160, 4320, 8760)) %>%   # 0, 3mo, 6mo, 12mo
  group_by(scenario, time) %>%
  summarise(
    GH_mean    = mean(GH_mean, na.rm=TRUE),
    IGF1_mean  = mean(IGF1_out, na.rm=TRUE),
    LVmass     = mean(LVmass, na.rm=TRUE),
    Glucose    = mean(Glucose_out, na.rm=TRUE),
    TumorVol   = mean(TumorVol, na.rm=TRUE),
    Remission  = mean(Remission, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    TimeLabel = case_when(
      time == 0    ~ "Baseline",
      time == 2160 ~ "3 months",
      time == 4320 ~ "6 months",
      time == 8760 ~ "12 months"
    )
  )

## ====================================================================
## VISUALIZATION
## ====================================================================
theme_acro <- theme_bw(base_size=12) +
  theme(legend.position="bottom",
        strip.background=element_rect(fill="#2c3e50"),
        strip.text=element_text(color="white"))

## Figure 1: GH and IGF-1 over time
p1 <- all_scenarios %>%
  filter(time <= 8760) %>%
  pivot_longer(cols=c(GH_mean, IGF1_out),
               names_to="marker", values_to="value") %>%
  mutate(marker=recode(marker,
    GH_mean="Plasma GH (ng/mL)",
    IGF1_out="Plasma IGF-1 (ng/mL)")) %>%
  ggplot(aes(x=time/720, y=value, color=scenario)) +
  geom_line(linewidth=0.8, alpha=0.85) +
  facet_wrap(~marker, scales="free_y", nrow=2) +
  geom_hline(data=data.frame(
    marker=c("Plasma GH (ng/mL)","Plasma IGF-1 (ng/mL)"),
    yint=c(2.5, 250)),
    aes(yintercept=yint), linetype="dashed", color="red", linewidth=0.7) +
  labs(title="Acromegaly QSP: GH & IGF-1 Response by Treatment Scenario",
       x="Time (months)", y="Concentration", color="Scenario") +
  scale_x_continuous(breaks=0:12) +
  theme_acro
print(p1)

## Figure 2: LV mass index and Glucose
p2 <- all_scenarios %>%
  filter(time <= 8760) %>%
  pivot_longer(cols=c(LVmass, Glucose_out),
               names_to="param", values_to="value") %>%
  mutate(param=recode(param,
    LVmass="LV Mass Index (g/m²)",
    Glucose_out="Fasting Glucose (mg/dL)")) %>%
  ggplot(aes(x=time/720, y=value, color=scenario)) +
  geom_line(linewidth=0.8, alpha=0.85) +
  facet_wrap(~param, scales="free_y", nrow=2) +
  geom_hline(data=data.frame(
    param=c("LV Mass Index (g/m²)","Fasting Glucose (mg/dL)"),
    yint=c(115, 100)),
    aes(yintercept=yint), linetype="dashed", color="blue", linewidth=0.7) +
  labs(title="Acromegaly QSP: Cardiovascular & Metabolic Outcomes",
       x="Time (months)", y="Value", color="Scenario") +
  scale_x_continuous(breaks=0:12) +
  theme_acro
print(p2)

## Figure 3: Dose-response — Octreotide LAR
p3 <- dr_data %>%
  filter(time <= 8760) %>%
  ggplot(aes(x=time/720, y=IGF1_out, color=factor(dose), group=label)) +
  geom_line(linewidth=1.0) +
  geom_hline(yintercept=250, linetype="dashed", color="red") +
  labs(title="Dose-Response: Octreotide LAR on IGF-1 Normalization",
       x="Time (months)", y="Plasma IGF-1 (ng/mL)",
       color="Octreotide LAR\nDose (mg)") +
  scale_x_continuous(breaks=0:12) +
  theme_acro
print(p3)

## Figure 4: Tumor volume over time
p4 <- all_scenarios %>%
  filter(time <= 8760) %>%
  ggplot(aes(x=time/720, y=TumorVol, color=scenario)) +
  geom_line(linewidth=0.8, alpha=0.85) +
  labs(title="Pituitary Adenoma Volume by Treatment",
       x="Time (months)", y="Tumor Volume (cm³)", color="Scenario") +
  scale_x_continuous(breaks=0:12) +
  theme_acro
print(p4)

## ====================================================================
## PRINT SUMMARY TABLE
## ====================================================================
cat("\n=== Acromegaly QSP — 12-Month Outcome Summary ===\n")
print(summary_tbl %>%
  filter(TimeLabel %in% c("Baseline","6 months","12 months")) %>%
  select(scenario, TimeLabel, GH_mean, IGF1_mean, LVmass, Glucose, Remission) %>%
  mutate(across(where(is.numeric), ~round(.x, 1))) %>%
  arrange(scenario, TimeLabel),
  n=Inf)
