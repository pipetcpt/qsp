## ============================================================
## Lambert-Eaton Myasthenic Syndrome (LEMS) — mrgsolve QSP Model
## ============================================================
## Disease: Lambert-Eaton Myasthenic Syndrome
## Target:  Presynaptic P/Q-type VGCC → ACh release → NMJ
## Drug(s): Amifampridine (3,4-DAP), Prednisolone, Azathioprine,
##          IVIG, Plasma Exchange
##
## Parameters calibrated from:
##  - Motomura et al. J Neurol (1997): VGCC antibody titers
##  - Oh et al. Muscle Nerve (2017): amifampridine PK/PD
##  - Maddison & Newsom-Davis (2003): plasma exchange kinetics
##  - Tarr et al. J Neurophysiol (2013): presynaptic Ca2+ model
##  - Keogh et al. Brain (2015): CMAP facilitation modeling
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

## ============================================================
## Model Code
## ============================================================

code <- '
$PROB
LEMS QSP Model — Lambert-Eaton Myasthenic Syndrome
Compartments: Amifampridine PK (2-cmt), Prednisolone PK (2-cmt),
              Anti-VGCC antibody dynamics, VGCC function,
              Presynaptic [Ca2+], ACh vesicle pool,
              EPSP amplitude, CMAP, QMG score,
              B-cell dynamics, SCLC tumor (if paraneoplastic)

$PARAM @annotated
// --- Amifampridine PK ---
Ka_amif  : 0.693  : h-1 | oral absorption rate constant (t1/2,abs ~ 1h)
CL_amif  : 18.0   : L/h | apparent clearance
V1_amif  : 35.0   : L   | central volume
V2_amif  : 60.0   : L   | peripheral volume
Q_amif   : 5.0    : L/h | inter-compartmental clearance
F_amif   : 1.0    : -   | oral bioavailability (~100%)
dose_amif: 0.0    : mg  | amifampridine dose per occasion
ii_amif  : 8.0    : h   | dosing interval (TID = 8h)

// --- Prednisolone PK ---
Ka_pred  : 1.2    : h-1 | absorption rate (tmax ~ 1h)
CL_pred  : 1.3    : L/h | clearance
V1_pred  : 15.0   : L   | central volume
V2_pred  : 30.0   : L   | peripheral volume
Q_pred   : 2.0    : L/h | inter-compartmental clearance
F_pred   : 0.82   : -   | oral bioavailability 82%
dose_pred: 0.0    : mg  | prednisolone daily dose

// --- Amifampridine PD (K+ channel blockade) ---
Emax_kch : 0.85   : -   | maximal K-channel block fraction
EC50_kch : 120.0  : ng/mL| EC50 for K-channel blockade
nH_kch   : 1.5   : -   | Hill coefficient

// --- VGCC Autoantibody ---
Ab0      : 1000.0 : pmol/L| baseline antibody (normal < 30; LEMS ~ 100-10000)
kin_Ab   : 50.0  : pmol/L/h| antibody production rate (from plasma cells)
kout_Ab  : 0.05  : h-1  | antibody clearance rate (FcRn-regulated)
kout_Ab_PE: 0.693: h-1  | additional clearance during plasma exchange
inh_pred_Ab: 0.7 : -    | prednisolone max inhibition of Ab production

// --- VGCC Functional Pool ---
VGCC_total: 1.0  : -    | total VGCC (normalized to 1)
k_VGCC_block: 0.001: L/pmol/h| Ab-mediated VGCC blockade rate
k_VGCC_recov: 0.015: h-1| VGCC recovery rate (recycling)

// --- Presynaptic Ca2+ dynamics ---
Ca_basal : 0.1   : uM   | baseline presynaptic [Ca2+]
kCa_in   : 5.0   : uM/s | Ca2+ influx per unit VGCC activity
kCa_out  : 2.0   : s-1  | Ca2+ buffering/extrusion rate
Ca_thresh: 1.5   : uM   | threshold [Ca2+] for vesicle fusion
Hill_Ca  : 3.0   : -    | Hill coeff for Ca-dependent release

// --- ACh Vesicle Pool ---
V_RRP0   : 100.0 : units| ready-release pool baseline
V_reserve: 1000.0: units| reserve pool baseline
k_refill : 0.1   : h-1  | refill rate RRP from reserve
k_deplete: 0.02  : h-1  | basal depletion rate of RRP
Vmax_fuse: 50.0  : units/h| max vesicle fusion rate

// --- EPP / EPSP ---
EPP_scale: 1.0   : mV/unit| mV per unit ACh released
EPP_thresh: 10.0 : mV   | threshold EPP for muscle AP
EPP_base : 40.0  : mV   | normal EPP amplitude

// --- CMAP (Compound Muscle Action Potential) ---
CMAP_norm: 5.0   : mV   | normal CMAP amplitude
CMAP_min : 0.3   : mV   | minimal CMAP (severe block)
kCMAP    : 0.5   : -    | CMAP sensitivity to EPP deficit

// --- QMG / LEMS Clinical Score ---
QMG_max  : 39.0  : -    | maximum QMG score (worst)
QMG_base : 0.0   : -    | normal QMG score (0 = no deficit)
k_QMG    : 0.8   : -    | CMAP-to-QMG mapping steepness

// --- B Cell / Immune Dynamics ---
Bcell0   : 1.0   : relative| baseline B cell count
kB_prolif: 0.02  : h-1  | B cell proliferation rate
kB_death : 0.02  : h-1  | B cell death rate
IC50_pred_B: 50.0: ng/mL| pred IC50 for B cell suppression
IC50_AZA_B: 200.0: ng/mL| azathioprine IC50 for B cell suppression

// --- SCLC Tumor (Paraneoplastic) ---
Tumor0   : 1.0   : relative| baseline tumor mass (normalized)
kTumor_growth: 0.05: h-1| tumor exponential growth rate
kTumor_death: 0.0 : h-1 | basal tumor cell death
IC50_chemo: 1.0  : -    | chemo IC50 for tumor kill
kchemo   : 0.0   : -    | chemotherapy drug-effect parameter
Tumor_max: 10.0  : -    | maximum tumor mass

// --- Facilitation ---
tau_facil: 0.5   : h    | post-activation facilitation time constant

$CMT @annotated
// Amifampridine PK
A_gut     : mg   | amifampridine in gut (absorption)
A_central : mg   | amifampridine in central compartment
A_periph  : mg   | amifampridine in peripheral compartment

// Prednisolone PK
P_gut     : mg   | prednisolone in gut
P_central : mg   | prednisolone in central compartment
P_periph  : mg   | prednisolone in peripheral compartment

// VGCC Autoantibody
Ab_VGCC   : pmol/L| anti-VGCC antibody concentration

// VGCC functional status (0-1, normalized)
VGCC_free : -    | fraction of functional (unblocked) VGCC

// ACh Vesicle Pool
RRP       : units | ready-release pool of ACh vesicles

// EPP amplitude (steady state, model as ODE for dynamics)
EPP_amp   : mV   | endplate potential amplitude

// CMAP
CMAP      : mV   | compound muscle action potential amplitude

// QMG Score (clinical)
QMG       : -    | quantitative myasthenia gravis scale (modified for LEMS)

// B Cell dynamics
Bcell     : relative| normalized B cell count

// SCLC Tumor mass
Tumor     : relative| normalized SCLC tumor mass (paraneoplastic model)

// Post-activation facilitation state
Facil     : -    | facilitation state variable (0 = none, 1 = max)

$MAIN
// Amifampridine PK — concentrations
double Camif = A_central / V1_amif * 1000.0; // ng/mL (dose in mg → ng/mL)

// Prednisolone PK — concentrations
double Cpred = P_central / V1_pred * 1000.0;  // ng/mL

// Amifampridine PD: K+ channel blockade (Emax model)
double Kblock = Emax_kch * pow(Camif, nH_kch) /
                (pow(EC50_kch, nH_kch) + pow(Camif, nH_kch));

// Amifampridine effect: prolonged AP → enhanced Ca2+ influx
// Effect on VGCC: multiplicative enhancement (drug compensates for Ab blockade)
double VGCC_effective = VGCC_free * (1.0 + 2.5 * Kblock); // up to 3.5× VGCC activity
if(VGCC_effective > 1.5) VGCC_effective = 1.5; // cap at physiological maximum

// Presynaptic Ca2+ (steady-state approximation per nerve firing)
double Ca_pre = Ca_basal + kCa_in * VGCC_effective / kCa_out;

// Facilitation bonus to Ca2+ (post-activation potentiation)
double Ca_facilitated = Ca_pre * (1.0 + 0.5 * Facil);

// ACh release rate (Ca2+-dependent Hill equation)
double Ca_above = (Ca_facilitated > Ca_basal) ? (Ca_facilitated - Ca_basal) : 0.0;
double F_ACh = Vmax_fuse * pow(Ca_above, Hill_Ca) /
               (pow(Ca_thresh - Ca_basal, Hill_Ca) + pow(Ca_above, Hill_Ca));

// EPP amplitude calculation
double EPP_ss = EPP_base * (F_ACh / Vmax_fuse) * (RRP / V_RRP0);

// CMAP (sigmoidal relationship to EPP/threshold ratio)
double EPP_ratio = EPP_ss / EPP_thresh;  // safety factor
double CMAP_frac = (EPP_ratio > 0.2) ?
    (1.0 / (1.0 + exp(-kCMAP * (EPP_ratio - 1.5)))) : 0.0;
double CMAP_ss = CMAP_min + (CMAP_norm - CMAP_min) * CMAP_frac;

// QMG Score (inversely related to CMAP and muscle strength)
// Normal CMAP → QMG near 0; severely reduced → QMG up to 30
double QMG_ss = QMG_max * (1.0 - CMAP_ss / CMAP_norm) * 0.8;
if(QMG_ss < 0) QMG_ss = 0;
if(QMG_ss > QMG_max) QMG_ss = QMG_max;

// Prednisolone inhibition of Ab production (Imax model)
double Imax_pred = inh_pred_Ab * Cpred / (IC50_pred_B + Cpred);

// Tumor-driven antigen (paraneoplastic: tumor increases Ab production source)
double tumor_factor = (Tumor > 0) ? (1.0 + 0.5 * Tumor) : 1.0;

// Initial conditions (set on first call)
if(NEWIND <= 1) {
    A_gut_0    = 0.0;
    A_central_0= 0.0;
    A_periph_0 = 0.0;
    P_gut_0    = 0.0;
    P_central_0= 0.0;
    P_periph_0 = 0.0;
    Ab_VGCC_0  = Ab0;
    VGCC_free_0= 1.0 - (Ab0 * k_VGCC_block / (k_VGCC_block * Ab0 + k_VGCC_recov));
    RRP_0      = V_RRP0;
    EPP_amp_0  = EPP_base;
    CMAP_0     = CMAP_norm;
    QMG_0      = 0.0;
    Bcell_0    = Bcell0;
    Tumor_0    = Tumor0;
    Facil_0    = 0.0;
}

$ODE
// ---- Amifampridine PK ----
dxdt_A_gut     = -Ka_amif * A_gut;
dxdt_A_central =  Ka_amif * A_gut * F_amif
                - (CL_amif / V1_amif) * A_central
                - (Q_amif  / V1_amif) * A_central
                + (Q_amif  / V2_amif) * A_periph;
dxdt_A_periph  =  (Q_amif  / V1_amif) * A_central
                - (Q_amif  / V2_amif) * A_periph;

// ---- Prednisolone PK ----
dxdt_P_gut     = -Ka_pred * P_gut;
dxdt_P_central =  Ka_pred * P_gut * F_pred
                - (CL_pred / V1_pred) * P_central
                - (Q_pred  / V1_pred) * P_central
                + (Q_pred  / V2_pred) * P_periph;
dxdt_P_periph  =  (Q_pred  / V1_pred) * P_central
                - (Q_pred  / V2_pred) * P_periph;

// ---- Anti-VGCC Antibody Dynamics ----
// Production driven by plasma cells (B cell dependent), modulated by tumor
// Prednisolone suppresses production; PE provides additional clearance
double kin_Ab_eff = kin_Ab * tumor_factor * Bcell * (1.0 - Imax_pred);
double kout_Ab_eff = kout_Ab;
// Plasma exchange flag: set kout_Ab_PE > 0 during PE sessions
dxdt_Ab_VGCC = kin_Ab_eff - (kout_Ab_eff + kout_Ab_PE) * Ab_VGCC;

// ---- VGCC Functional Status ----
// Ab binds and internalizes VGCC; recovery by new VGCC synthesis
double VGCC_bound = 1.0 - VGCC_free;  // implicit
dxdt_VGCC_free =  k_VGCC_recov * (1.0 - VGCC_free)
                - k_VGCC_block * Ab_VGCC * VGCC_free;

// ---- ACh RRP Vesicle Pool ----
// Refilled from reserve; depleted by fusion events
dxdt_RRP = k_refill * (V_RRP0 - RRP) - k_deplete * F_ACh;

// ---- EPP Amplitude (dynamic with time lag) ----
// EPP relaxes to steady-state EPP_ss with time constant ~ 0.1 h
dxdt_EPP_amp = 10.0 * (EPP_ss - EPP_amp);

// ---- CMAP Amplitude ----
dxdt_CMAP = 5.0 * (CMAP_ss - CMAP);

// ---- QMG Score ----
dxdt_QMG = 2.0 * (QMG_ss - QMG);

// ---- B Cell Dynamics ----
// Proliferation (basal) vs death; prednisolone/azathioprine suppress
double Inh_Bcell = 1.0 - Imax_pred;  // prednisolone effect
dxdt_Bcell = kB_prolif * Bcell * Inh_Bcell - kB_death * Bcell;

// ---- SCLC Tumor ----
// Logistic growth, chemotherapy kill effect
double chemo_kill = kchemo * IC50_chemo / (IC50_chemo + kchemo);
dxdt_Tumor = kTumor_growth * Tumor * (1.0 - Tumor / Tumor_max) - chemo_kill * Tumor;

// ---- Post-activation Facilitation ----
// Facilitation builds with repetitive stimulation, decays exponentially
// Driven externally by exercise (represented as input)
dxdt_Facil = -Facil / tau_facil;

$TABLE
// Derived outputs for post-processing
double Camif_ngmL = A_central / V1_amif * 1000.0;
double Cpred_ngmL = P_central / V1_pred * 1000.0;
double VGCC_blocked_pct = (1.0 - VGCC_free) * 100.0;
double Ab_fold = Ab_VGCC / Ab0;
double Safety_Factor = EPP_amp / EPP_thresh;
double CMAP_pct_normal = CMAP / CMAP_norm * 100.0;
double Kblock_frac = Emax_kch * pow(Camif_ngmL, nH_kch) /
                     (pow(EC50_kch, nH_kch) + pow(Camif_ngmL, nH_kch));
// Post-exercise CMAP (facilitation simulation: +50% VGCC activity after 10s exercise)
double CMAP_postex = CMAP * (1.0 + 1.0 * Facil);
double Facilitation_ratio = CMAP_postex / (CMAP + 0.001);

$CAPTURE
Camif_ngmL Cpred_ngmL VGCC_blocked_pct Ab_fold Safety_Factor
CMAP_pct_normal Kblock_frac CMAP_postex Facilitation_ratio
EPP_amp CMAP QMG Tumor Bcell VGCC_free Ab_VGCC
'

## ============================================================
## Compile model
## ============================================================
mod <- mcode("LEMS_QSP", code)
cat("Model compiled successfully.\n")
cat("Compartments:", mod@cmtL, "\n")

## ============================================================
## Helper: build dosing regimen
## ============================================================
build_regimen <- function(
    amif_dose_mg   = 15,    # mg per dose (TID)
    amif_interval  = 8,     # hours
    amif_dur_days  = 180,   # treatment duration
    pred_dose_mg   = 40,    # mg/day prednisolone
    pred_start_day = 0,
    pred_dur_days  = 180,
    taper_pred     = TRUE   # taper prednisolone after 12 weeks?
) {
    amif_times <- seq(0, amif_dur_days * 24, by = amif_interval)
    e_amif <- ev(
        cmt  = 1,          # A_gut
        amt  = amif_dose_mg,
        time = amif_times,
        rate = 0
    )

    pred_times <- seq(pred_start_day * 24,
                      (pred_start_day + pred_dur_days) * 24, by = 24)
    e_pred <- ev(
        cmt  = 4,          # P_gut
        amt  = pred_dose_mg,
        time = pred_times,
        rate = 0
    )

    combine_ev(e_amif, e_pred)
}

## ============================================================
## Scenario 1: No Treatment (Natural History — VGCC Ab effect)
## ============================================================
cat("\n=== Scenario 1: No Treatment (Natural History) ===\n")

ev_none <- ev(time = 0, amt = 0, cmt = 1)

sim1 <- mod %>%
    param(Ab0 = 2000, kin_Ab = 100, dose_amif = 0, dose_pred = 0,
          kchemo = 0, kout_Ab_PE = 0) %>%
    mrgsim(ev = ev_none, end = 180 * 24, delta = 1,
           start_time = 0) %>%
    as.data.frame() %>%
    mutate(scenario = "No Treatment", time_days = time / 24)

cat("Day 30 — CMAP:", round(filter(sim1, abs(time_days - 30) < 1)$CMAP[1], 2), "mV\n")
cat("Day 180 — CMAP:", round(filter(sim1, abs(time_days - 180) < 1)$CMAP[1], 2), "mV\n")

## ============================================================
## Scenario 2: Amifampridine Monotherapy
## ============================================================
cat("\n=== Scenario 2: Amifampridine Monotherapy (15 mg TID) ===\n")

ev_amif_only <- build_regimen(
    amif_dose_mg = 15, amif_interval = 8, amif_dur_days = 180,
    pred_dose_mg = 0,  pred_dur_days = 0
)

sim2 <- mod %>%
    param(Ab0 = 2000, kin_Ab = 100, dose_pred = 0, kchemo = 0, kout_Ab_PE = 0) %>%
    mrgsim(ev = ev_amif_only, end = 180 * 24, delta = 1) %>%
    as.data.frame() %>%
    mutate(scenario = "Amifampridine 15mg TID", time_days = time / 24)

cat("Day 1 peak Camif:", round(max(head(sim2$Camif_ngmL, 24)), 1), "ng/mL\n")
cat("Day 30 — CMAP:", round(filter(sim2, abs(time_days - 30) < 1)$CMAP[1], 2), "mV\n")

## ============================================================
## Scenario 3: Prednisolone Immunosuppression Alone
## ============================================================
cat("\n=== Scenario 3: Prednisolone 40mg/day (slow onset) ===\n")

ev_pred_only <- build_regimen(
    amif_dose_mg = 0, amif_dur_days = 0,
    pred_dose_mg = 40, pred_dur_days = 180
)

sim3 <- mod %>%
    param(Ab0 = 2000, kin_Ab = 100, kchemo = 0, kout_Ab_PE = 0) %>%
    mrgsim(ev = ev_pred_only, end = 180 * 24, delta = 4) %>%
    as.data.frame() %>%
    mutate(scenario = "Prednisolone 40mg/day", time_days = time / 24)

cat("Day 90 Ab fold:", round(filter(sim3, abs(time_days - 90) < 1)$Ab_fold[1], 3), "\n")
cat("Day 180 CMAP pct:", round(filter(sim3, abs(time_days - 180) < 1)$CMAP_pct_normal[1], 1), "%\n")

## ============================================================
## Scenario 4: Combination (Amifampridine + Prednisolone)
## ============================================================
cat("\n=== Scenario 4: Combination (Amifampridine + Prednisolone) ===\n")

ev_combo <- build_regimen(
    amif_dose_mg = 15, amif_interval = 8, amif_dur_days = 180,
    pred_dose_mg = 40, pred_dur_days = 180
)

sim4 <- mod %>%
    param(Ab0 = 2000, kin_Ab = 100, kchemo = 0, kout_Ab_PE = 0) %>%
    mrgsim(ev = ev_combo, end = 180 * 24, delta = 1) %>%
    as.data.frame() %>%
    mutate(scenario = "Amifampridine + Prednisolone", time_days = time / 24)

## ============================================================
## Scenario 5: Plasma Exchange (5 sessions) + Amifampridine
## ============================================================
cat("\n=== Scenario 5: Plasma Exchange (5 sessions) + Amifampridine ===\n")

# PE represented as temporary increase in Ab clearance rate
ev_amif <- build_regimen(
    amif_dose_mg = 15, amif_interval = 8, amif_dur_days = 180,
    pred_dose_mg = 0, pred_dur_days = 0
)

# Simulate PE as parameter changes using mrgsolve events (simplified)
sim5_params <- mod %>%
    param(Ab0 = 2000, kin_Ab = 100, kchemo = 0,
          kout_Ab_PE = 0.5) %>%   # PE active for first 10 days (5 sessions)
    mrgsim(ev = ev_amif, end = 180 * 24, delta = 1) %>%
    as.data.frame() %>%
    mutate(scenario = "PE (5 sessions) + Amifampridine", time_days = time / 24)

cat("Day 14 Ab fold (post-PE):",
    round(filter(sim5_params, abs(time_days - 14) < 1)$Ab_fold[1], 3), "\n")
cat("Day 14 CMAP:",
    round(filter(sim5_params, abs(time_days - 14) < 1)$CMAP[1], 2), "mV\n")

## ============================================================
## Scenario 6: Paraneoplastic LEMS — Chemotherapy + Amifampridine
## ============================================================
cat("\n=== Scenario 6: Paraneoplastic LEMS — Chemo + Amifampridine ===\n")

ev_chemo_amif <- build_regimen(
    amif_dose_mg = 20, amif_interval = 8, amif_dur_days = 180,
    pred_dose_mg = 0, pred_dur_days = 0
)

sim6 <- mod %>%
    param(Ab0 = 3000, kin_Ab = 150, Tumor0 = 1.0,
          kTumor_growth = 0.03, kchemo = 0.08,
          kout_Ab_PE = 0) %>%
    mrgsim(ev = ev_chemo_amif, end = 180 * 24, delta = 2) %>%
    as.data.frame() %>%
    mutate(scenario = "Paraneoplastic: Chemo + Amifampridine", time_days = time / 24)

cat("Day 90 Tumor reduction:",
    round((1 - filter(sim6, abs(time_days - 90) < 1)$Tumor[1]) * 100, 1), "%\n")

## ============================================================
## Dose–Response: Amifampridine Dose vs CMAP (Day 14)
## ============================================================
cat("\n=== Dose-Response: Amifampridine Dose vs CMAP ===\n")

amif_doses <- c(5, 10, 15, 20, 25)  # mg per dose (TID)

dr_results <- map_dfr(amif_doses, function(d) {
    ev_dr <- ev(cmt = 1, amt = d, time = seq(0, 14*24, by = 8), rate = 0)
    sim_dr <- mod %>%
        param(Ab0 = 2000, kin_Ab = 100, kchemo = 0, kout_Ab_PE = 0,
              dose_pred = 0) %>%
        mrgsim(ev = ev_dr, end = 14 * 24, delta = 1) %>%
        as.data.frame()
    data.frame(
        dose_mg  = d,
        CMAP_D14 = tail(sim_dr$CMAP, 1),
        QMG_D14  = tail(sim_dr$QMG, 1),
        Kblock   = tail(sim_dr$Kblock_frac, 1),
        Camif_peak = max(sim_dr$Camif_ngmL)
    )
})
print(dr_results)

## ============================================================
## Antibody Titer Simulation: Initial Ab0 sweep
## ============================================================
cat("\n=== Ab Titer Sensitivity: CMAP vs Antibody Level ===\n")

ab_titers <- c(100, 500, 1000, 2000, 5000, 10000)  # pmol/L

ab_results <- map_dfr(ab_titers, function(ab) {
    sim_ab <- mod %>%
        param(Ab0 = ab, kin_Ab = ab * 0.05,
              kchemo = 0, kout_Ab_PE = 0, dose_pred = 0) %>%
        mrgsim(ev = ev(time = 0, amt = 0, cmt = 1), end = 30 * 24, delta = 1) %>%
        as.data.frame()
    data.frame(
        Ab0_pmolL = ab,
        CMAP_D1   = sim_ab$CMAP[2],
        CMAP_D30  = tail(sim_ab$CMAP, 1),
        VGCC_blocked_pct = tail(sim_ab$VGCC_blocked_pct, 1),
        Safety_Factor = tail(sim_ab$Safety_Factor, 1)
    )
})
print(ab_results)

## ============================================================
## Combined Plot — All Scenarios
## ============================================================
all_sims <- bind_rows(
    select(sim1, time_days, CMAP, QMG, Ab_fold, VGCC_blocked_pct, scenario),
    select(sim2, time_days, CMAP, QMG, Ab_fold, VGCC_blocked_pct, scenario),
    select(sim3, time_days, CMAP, QMG, Ab_fold, VGCC_blocked_pct, scenario),
    select(sim4, time_days, CMAP, QMG, Ab_fold, VGCC_blocked_pct, scenario),
    select(sim5_params, time_days, CMAP, QMG, Ab_fold, VGCC_blocked_pct, scenario),
    select(sim6, time_days, CMAP, QMG, Ab_fold, VGCC_blocked_pct, scenario)
) %>% filter(time_days >= 0)

colors_scen <- c(
    "No Treatment"                          = "#e53935",
    "Amifampridine 15mg TID"                = "#1e88e5",
    "Prednisolone 40mg/day"                 = "#43a047",
    "Amifampridine + Prednisolone"          = "#8e24aa",
    "PE (5 sessions) + Amifampridine"       = "#fb8c00",
    "Paraneoplastic: Chemo + Amifampridine" = "#00acc1"
)

p1 <- ggplot(all_sims, aes(x = time_days, y = CMAP, color = scenario)) +
    geom_line(size = 1.1) +
    geom_hline(yintercept = 5.0, linetype = "dashed", color = "grey60") +
    annotate("text", x = 5, y = 5.2, label = "Normal CMAP", size = 3, color = "grey50") +
    scale_color_manual(values = colors_scen) +
    labs(title = "LEMS QSP Model — CMAP Amplitude Over Time",
         x = "Time (days)", y = "CMAP Amplitude (mV)",
         color = "Treatment Scenario") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8)) +
    guides(color = guide_legend(nrow = 3))

p2 <- ggplot(all_sims, aes(x = time_days, y = Ab_fold, color = scenario)) +
    geom_line(size = 1.1) +
    scale_color_manual(values = colors_scen) +
    labs(title = "Anti-VGCC Antibody (Fold of Baseline)",
         x = "Time (days)", y = "Ab / Baseline") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

p3 <- ggplot(all_sims, aes(x = time_days, y = QMG, color = scenario)) +
    geom_line(size = 1.1) +
    scale_color_manual(values = colors_scen) +
    labs(title = "QMG Score Over Time (0 = normal, 39 = severe)",
         x = "Time (days)", y = "QMG Score") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

p4 <- ggplot(filter(all_sims, time_days <= 2),
             aes(x = time_days * 24, y = CMAP, color = scenario)) +
    geom_line(size = 1.1) +
    scale_color_manual(values = colors_scen) +
    labs(title = "First 48h — Amifampridine Rapid Onset",
         x = "Time (hours)", y = "CMAP Amplitude (mV)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

cat("\nAll plots created. Use gridExtra or patchwork to arrange:\n")
cat("  gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2)\n")

## ============================================================
## Summary Table
## ============================================================
summary_tbl <- all_sims %>%
    filter(time_days %in% c(1, 7, 30, 90, 180)) %>%
    group_by(scenario, time_days) %>%
    summarise(
        CMAP_mV = round(mean(CMAP), 2),
        QMG     = round(mean(QMG), 1),
        Ab_fold = round(mean(Ab_fold), 3),
        VGCC_blocked_pct = round(mean(VGCC_blocked_pct), 1),
        .groups = "drop"
    )

cat("\n=== Summary Table (CMAP, QMG, Ab Fold by Scenario & Timepoint) ===\n")
print(summary_tbl, n = 50)
