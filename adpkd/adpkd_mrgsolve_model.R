## ============================================================
##  ADPKD QSP Model — mrgsolve ODE Implementation
##  Autosomal Dominant Polycystic Kidney Disease
##
##  Key mechanisms:
##    • PKD1/PKD2 → PC1/PC2 loss → ↓Ca²⁺ → ↑cAMP → cyst growth
##    • cAMP → PKA/CFTR → fluid secretion → TKV enlargement
##    • mTOR/ERK → epithelial proliferation
##    • RAAS activation → hypertension
##
##  Drug targets modelled:
##    1) Tolvaptan   (V2R antagonist, 45 mg AM / 15 mg PM)
##    2) Everolimus  (mTOR inhibitor, 2.5 mg/day)
##    3) Octreotide LAR (SSTR agonist, 30 mg/28 days)
##    4) ACEi/ARB    (RAAS blocker, ramipril equivalent)
##
##  PK calibration:
##    Tolvaptan:  Torres VE et al. NEJM 2012 (TEMPO 3:4)
##    Everolimus: Serra AL et al. NEJM 2010; Walz G NEJM 2010
##    Octreotide: Meijer E et al. JAMA Intern Med 2011 (DIPAK-1)
##    ACEi:       Schrier RW et al. NEJM 2014 (HALT-PKD)
##
##  Disease progression:
##    TKV growth ~5.5%/yr (Grantham JJ, NEJM 2006; CRISP cohort)
##    eGFR decline ~3.5 mL/min/1.73m²/yr (Chapman AB, KI 2003)
##    Tolvaptan: TKV +4.4% vs +8.0% (TEMPO 3:4 Phase 3)
##    Everolimus: TKV growth ~41% reduction at 1yr (SIRENA)
##
##  Time unit: HOURS
##  Disease progression spans 3 years (26,280 h)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)

## ============================================================
##  1. MODEL DEFINITION
## ============================================================

adpkd_model <- '
$PROB
ADPKD QSP model — mrgsolve v1.0
Mechanistic PK/PD with TKV growth, eGFR decline, RAAS, cAMP, mTOR

$PARAM @annotated
// === Tolvaptan PK ===
KA_TOLV  : 1.00 : Tolvaptan absorption rate constant (/h)
CL_TOLV  : 4.00 : Tolvaptan clearance (L/h)
V1_TOLV  : 10.0 : Tolvaptan central volume (L)
Q_TOLV   : 2.00 : Tolvaptan intercompartmental CL (L/h)
V2_TOLV  : 20.0 : Tolvaptan peripheral volume (L)
F_TOLV   : 0.56 : Tolvaptan oral bioavailability

// === Tolvaptan PD (V2R inhibition) ===
EC50_TOLV : 50.0  : Tolvaptan EC50 for V2R blockade (ng/mL)
EMAX_TOLV :  1.0  : Tolvaptan maximal V2R inhibition (0-1)
HILL_TOLV :  1.5  : Hill coefficient for tolvaptan

// === Everolimus PK ===
KA_EVER  : 0.40 : Everolimus absorption rate constant (/h)
CL_EVER  : 14.0 : Everolimus clearance (L/h)
V1_EVER  : 20.0 : Everolimus central volume (L)
Q_EVER   : 5.00 : Everolimus intercompartmental CL (L/h)
V2_EVER  : 50.0 : Everolimus peripheral volume (L)

// === Everolimus PD (mTOR inhibition) ===
EC50_EVER : 5.00 : Everolimus EC50 for mTOR inhibition (ng/mL)
EMAX_EVER : 0.80 : Everolimus maximal mTOR inhibition (0-1)
HILL_EVER : 1.00 : Hill coefficient for everolimus

// === Octreotide LAR PK ===
KREL_OCT : 0.002 : Octreotide depot release rate (/h, ~21 day t1/2)
CL_OCT   : 0.50  : Octreotide clearance (L/h)
V1_OCT   : 8.00  : Octreotide central volume (L)

// === ACEi/ARB PK ===
KA_ACEI  : 0.80  : ACEi absorption rate constant (/h)
CL_ACEI  : 5.00  : ACEi clearance (L/h)
V1_ACEI  : 15.0  : ACEi central volume (L)
EC50_ACEI: 0.50  : ACEi EC50 for RAAS inhibition (mg/L)
EMAX_ACEI: 0.70  : ACEi maximal RAAS inhibition

// === Ca2+-cAMP Axis ===
KIN_CAMP  : 0.50  : Basal cAMP production rate
KOUT_CAMP : 0.50  : cAMP degradation rate (/h)
CAMP_SS   : 1.00  : cAMP steady-state baseline (normalized)
CA2_SCALE : 1.50  : cAMP amplification due to ↓Ca²⁺ in ADPKD

// === mTOR Axis ===
KIN_MTOR  : 0.30  : Basal mTOR synthesis/activation rate
KOUT_MTOR : 0.30  : mTOR deactivation rate (/h)
MTOR_SS   : 1.00  : mTOR baseline (normalized)
CAMP_MTOR : 0.30  : cAMP contribution to mTOR activation

// === RAAS / Blood Pressure ===
KIN_ANGII : 0.20  : Basal Ang II production rate
KOUT_ANGII: 0.20  : Ang II degradation (/h)
ANGII_SS  : 1.00  : Ang II baseline (normalized)
RENIN_TKV : 0.30  : RAAS amplification per unit TKV increase
KIN_BP    : 26.0  : Blood pressure synthesis term (mmHg * /h)
KOUT_BP   : 0.20  : BP regulation rate (/h)
ANGII_BP  : 5.00  : Ang II to BP scaling factor
BP_SS     : 130.0 : Basal blood pressure (mmHg)

// === Disease Progression ===
TKV0      : 1500.0 : Baseline total kidney volume (mL)
KGROW_TKV : 6.28e-6: TKV growth rate (/h = 5.5%/yr)
EGFR0     : 70.0   : Baseline eGFR (mL/min/1.73m²)
KDECL_EGFR: 4.00e-4: eGFR decline rate (mL/min/h = 3.5/yr)
TKV_EGFR  : 0.50   : TKV compression effect on eGFR decline
BP_EGFR   : 0.30   : Hypertension effect on eGFR decline

// === TKV growth modulators (fraction of growth) ===
FCAMP     : 0.50   : Fraction of TKV growth driven by cAMP
FMTOR     : 0.30   : Fraction of TKV growth driven by mTOR
FBASE     : 0.20   : Irreducible baseline growth fraction

// === Urine Osmolality PD ===
UOSM_BASE : 600.0  : Baseline urine osmolality (mOsm/kg)
UOSM_MIN  : 50.0   : Minimum urine osmolality (free water)
KOUT_UOSM : 0.50   : Urine osm equilibration rate (/h)

$CMT @annotated
AGUT    : Tolvaptan gut (mg)
ACENT   : Tolvaptan central (mg)
APERI   : Tolvaptan peripheral (mg)
EGUT    : Everolimus gut (mg)
ECENT   : Everolimus central (mg)
EPERI   : Everolimus peripheral (mg)
OCTDEP  : Octreotide LAR depot (mg)
OCTCENT : Octreotide central (ug)
ACEI_GUT  : ACEi gut (mg)
ACEI_CENT : ACEi central (mg)
AVP_ST  : AVP level (normalized, 1 = normal)
CAMP_ST : Collecting duct cAMP (normalized, 1 = baseline)
MTOR_ST : mTOR activity (normalized)
ANGII_ST: Angiotensin II level (normalized)
BP_ST   : Blood pressure (mmHg)
TKV_ST  : Total kidney volume (mL)
EGFR_ST : Estimated GFR (mL/min/1.73m²)
UOSM_ST : Urine osmolality (mOsm/kg)
NEPH_ST : Functional nephron fraction (0-1)

$MAIN
// Bioavailability for tolvaptan
F_AGUT = F_TOLV;

// Initial conditions for disease variables
AVP_ST_0   = 1.0;
CAMP_ST_0  = CA2_SCALE;   // elevated cAMP at baseline (ADPKD)
MTOR_ST_0  = 1.2;          // slightly elevated mTOR
ANGII_ST_0 = 1.0;
BP_ST_0    = BP_SS;
TKV_ST_0   = TKV0;
EGFR_ST_0  = EGFR0;
UOSM_ST_0  = UOSM_BASE;
NEPH_ST_0  = 1.0;

$ODE
// ---- Tolvaptan PK ----
double Cp_tolv_mgl  = ACENT / V1_TOLV;           // mg/L = ug/mL
double Cp_tolv_ngml = Cp_tolv_mgl * 1000.0;      // ng/mL

dxdt_AGUT  = -KA_TOLV * AGUT;
dxdt_ACENT = KA_TOLV * AGUT
             - (CL_TOLV / V1_TOLV) * ACENT
             - (Q_TOLV  / V1_TOLV) * ACENT
             + (Q_TOLV  / V2_TOLV) * APERI;
dxdt_APERI = (Q_TOLV / V1_TOLV) * ACENT
             - (Q_TOLV / V2_TOLV) * APERI;

// ---- Everolimus PK ----
double Ce_ever_mgl  = ECENT / V1_EVER;
double Ce_ever_ngml = Ce_ever_mgl * 1000.0;

dxdt_EGUT  = -KA_EVER * EGUT;
dxdt_ECENT = KA_EVER * EGUT
             - (CL_EVER / V1_EVER) * ECENT
             - (Q_EVER  / V1_EVER) * ECENT
             + (Q_EVER  / V2_EVER) * EPERI;
dxdt_EPERI = (Q_EVER / V1_EVER) * ECENT
             - (Q_EVER / V2_EVER) * EPERI;

// ---- Octreotide LAR PK ----
double Cp_oct_ngml = (OCTCENT / V1_OCT) * 1000.0;  // ug/L = ng/mL

dxdt_OCTDEP  = -KREL_OCT * OCTDEP;
dxdt_OCTCENT =  KREL_OCT * OCTDEP * 1000.0   // mg → ug conversion
                - (CL_OCT / V1_OCT) * OCTCENT;

// ---- ACEi PK ----
double Cp_acei_mgl = ACEI_CENT / V1_ACEI;

dxdt_ACEI_GUT  = -KA_ACEI * ACEI_GUT;
dxdt_ACEI_CENT = KA_ACEI * ACEI_GUT
                 - (CL_ACEI / V1_ACEI) * ACEI_CENT;

// ---- Drug PD effects (0-1 scale) ----
// Tolvaptan V2R inhibition (Hill Emax)
double INH_TOLV = EMAX_TOLV *
                  pow(Cp_tolv_ngml, HILL_TOLV) /
                  (pow(EC50_TOLV, HILL_TOLV) + pow(Cp_tolv_ngml, HILL_TOLV));

// Everolimus mTOR inhibition
double INH_EVER = EMAX_EVER * Ce_ever_ngml /
                  (EC50_EVER + Ce_ever_ngml);

// Octreotide SSTR effect on cAMP (EC50 = 1 ng/mL for SSTR2)
double INH_OCT = 0.40 * Cp_oct_ngml / (1.0 + Cp_oct_ngml);

// ACEi RAAS inhibition
double INH_ACEI = EMAX_ACEI * Cp_acei_mgl /
                  (EC50_ACEI + Cp_acei_mgl);

// ---- AVP dynamics (relatively stable in plasma) ----
dxdt_AVP_ST = KIN_CAMP - KOUT_CAMP * AVP_ST;

// ---- cAMP dynamics ----
// cAMP driven by AVP/V2R; reduced by tolvaptan and octreotide
// In ADPKD: ↓Ca²⁺ → removes Ca²⁺-mediated AC inhibition → ↑cAMP
double CAMP_IN  = KIN_CAMP * AVP_ST * CA2_SCALE
                  * (1.0 - INH_TOLV)       // tolvaptan blocks V2R
                  * (1.0 - INH_OCT);       // octreotide (Gi) reduces AC
double CAMP_OUT = KOUT_CAMP * CAMP_ST;
dxdt_CAMP_ST = CAMP_IN - CAMP_OUT;

// ---- mTOR dynamics ----
// mTOR activated by growth factors, indirectly by cAMP (PKA→B-Raf→ERK→mTOR)
double MTOR_IN  = KIN_MTOR * (1.0 + CAMP_MTOR * CAMP_ST);
double MTOR_OUT = KOUT_MTOR * MTOR_ST * (1.0 + INH_EVER);
dxdt_MTOR_ST = MTOR_IN - MTOR_OUT;

// ---- RAAS / Angiotensin II dynamics ----
// Renin release increases with TKV (intrarenal pressure)
double RENIN_FACTOR = 1.0 + RENIN_TKV * (TKV_ST / TKV0 - 1.0);
double ANGII_IN  = KIN_ANGII * RENIN_FACTOR * (1.0 - INH_ACEI);
double ANGII_OUT = KOUT_ANGII * ANGII_ST;
dxdt_ANGII_ST = ANGII_IN - ANGII_OUT;

// ---- Blood pressure dynamics ----
double BP_IN  = KIN_BP + ANGII_BP * ANGII_ST;
double BP_OUT = KOUT_BP * BP_ST;
dxdt_BP_ST = BP_IN - BP_OUT;

// ---- TKV growth model ----
// Growth driven by cAMP (fluid secretion + proliferation) and mTOR (proliferation)
// Treatment reduces cAMP and/or mTOR
double TKV_MOD = FBASE
               + FCAMP * (CAMP_ST / CAMP_SS)
               + FMTOR * (MTOR_ST / MTOR_SS);
// Normalize so that baseline TKV_MOD = 1 (FBASE + FCAMP*CA2_SCALE + FMTOR*1.2)
double TKV_MOD_NORM = FBASE + FCAMP * CA2_SCALE + FMTOR * 1.2;
double KGROW_EFF = KGROW_TKV * TKV_MOD / TKV_MOD_NORM;
dxdt_TKV_ST = KGROW_EFF * TKV_ST;

// ---- eGFR decline model ----
// Decline accelerates with TKV (compression) and hypertension
double TKV_COMPRESS = 1.0 + TKV_EGFR * (TKV_ST / TKV0 - 1.0);
double BP_DAMAGE    = 1.0 + BP_EGFR  * (BP_ST  / BP_SS  - 1.0);
// Ensure non-negative
if (TKV_COMPRESS < 0.5) TKV_COMPRESS = 0.5;
if (BP_DAMAGE    < 0.5) BP_DAMAGE    = 0.5;
double KDECL_EFF = KDECL_EGFR * TKV_COMPRESS * BP_DAMAGE;
dxdt_EGFR_ST = (EGFR_ST > 5.0) ? -KDECL_EFF : 0.0;

// ---- Functional nephron fraction ----
dxdt_NEPH_ST = -KDECL_EFF / EGFR0;

// ---- Urine osmolality (tolvaptan PD biomarker) ----
// V2R blockade → ↓AQP2 → ↓water reabsorption → dilute urine
double UOSM_TARGET = UOSM_BASE * (1.0 - 0.90 * INH_TOLV)
                   + UOSM_MIN  * 0.90 * INH_TOLV;
if (UOSM_TARGET < UOSM_MIN) UOSM_TARGET = UOSM_MIN;
dxdt_UOSM_ST = KOUT_UOSM * (UOSM_TARGET - UOSM_ST);

$TABLE
double Cp_tolv     = Cp_tolv_ngml;            // ng/mL
double Ce_ever     = Ce_ever_ngml;            // ng/mL
double Cp_oct      = Cp_oct_ngml;             // ng/mL
double V2R_OCC     = INH_TOLV * 100.0;        // % V2R occupied
double mTOR_INH    = INH_EVER * 100.0;        // % mTOR inhibition
double RAAS_INH    = INH_ACEI * 100.0;        // % RAAS inhibition
double CAMP_norm   = CAMP_ST  / CAMP_SS;       // relative cAMP
double MTOR_norm   = MTOR_ST  / MTOR_SS;       // relative mTOR
double TKV_L       = TKV_ST   / 1000.0;        // L
double TKV_pct     = (TKV_ST  / TKV0 - 1.0) * 100.0;  // % change from baseline
double eGFR        = EGFR_ST;
double BP          = BP_ST;
double Uosm        = UOSM_ST;
double ANGII_rel   = ANGII_ST;

$CAPTURE Cp_tolv Ce_ever Cp_oct V2R_OCC mTOR_INH RAAS_INH CAMP_norm MTOR_norm TKV_L TKV_pct eGFR BP Uosm ANGII_rel
'

## Compile model
mod <- mcode("adpkd_qsp", adpkd_model)


## ============================================================
##  2. DOSING EVENTS FOR 5 TREATMENT SCENARIOS
## ============================================================

## Simulation time: 3 years = 26,280 hours
SIM_HOURS <- 3 * 365.25 * 24

## Helper function: daily dosing events
make_daily_dose <- function(cmt, amt, freq_per_day, start_h = 0, end_h = SIM_HOURS) {
  # freq_per_day = 1 or 2 (AM or AM+PM)
  times <- seq(start_h, end_h, by = 24 / freq_per_day)
  ev(cmt = cmt, amt = amt, time = times)
}

## ---- Scenario 1: Placebo ----
ev_placebo <- ev(cmt = "AGUT", amt = 0, time = 0)

## ---- Scenario 2: Tolvaptan 45/15 mg (low dose, 60 mg/day) ----
## AM dose 45 mg, PM dose 15 mg (8 hours later)
ev_tolv_low <- ev(cmt = "AGUT", amt = 45, time = seq(0, SIM_HOURS, by = 24)) +
               ev(cmt = "AGUT", amt = 15, time = seq(8, SIM_HOURS + 8, by = 24))

## ---- Scenario 3: Tolvaptan 90/30 mg (high dose, 120 mg/day) ----
ev_tolv_high <- ev(cmt = "AGUT", amt = 90, time = seq(0, SIM_HOURS, by = 24)) +
                ev(cmt = "AGUT", amt = 30, time = seq(8, SIM_HOURS + 8, by = 24))

## ---- Scenario 4: Everolimus 2.5 mg/day ----
ev_ever <- ev(cmt = "EGUT", amt = 2.5, time = seq(0, SIM_HOURS, by = 24))

## ---- Scenario 5: Tolvaptan 45/15 mg + ACEi 10 mg/day ----
ev_combo <- ev_tolv_low +
            ev(cmt = "ACEI_GUT", amt = 10, time = seq(0, SIM_HOURS, by = 24))


## ============================================================
##  3. SIMULATE ALL 5 SCENARIOS
## ============================================================

# Output times: daily for disease tracking + hourly for first 48h PK
times_full  <- seq(0, SIM_HOURS, by = 24)        # daily (disease)
times_pk    <- seq(0, 48, by = 0.5)               # hourly (PK day 1-2)

sim_disease <- function(ev_dose, label) {
  out <- mod %>%
    mrgsim_e(ev_dose, end = SIM_HOURS, delta = 24,
             recover = "time") %>%
    as.data.frame() %>%
    mutate(scenario = label,
           time_yr  = time / (365.25 * 24))
  out
}

sim_pk <- function(ev_dose, label) {
  out <- mod %>%
    mrgsim_e(ev_dose, end = 72, delta = 0.5,
             recover = "time") %>%
    as.data.frame() %>%
    mutate(scenario = label)
  out
}

cat("Simulating 5 ADPKD treatment scenarios...\n")
res_placebo   <- sim_disease(ev_placebo,   "1. Placebo")
res_tolv_low  <- sim_disease(ev_tolv_low,  "2. Tolvaptan 60 mg/day")
res_tolv_high <- sim_disease(ev_tolv_high, "3. Tolvaptan 120 mg/day")
res_ever      <- sim_disease(ev_ever,      "4. Everolimus 2.5 mg/day")
res_combo     <- sim_disease(ev_combo,     "5. Tolvaptan + ACEi")

pk_placebo    <- sim_pk(ev_placebo,   "1. Placebo")
pk_tolv_low   <- sim_pk(ev_tolv_low,  "2. Tolvaptan 60 mg/day")
pk_tolv_high  <- sim_pk(ev_tolv_high, "3. Tolvaptan 120 mg/day")

res_all    <- bind_rows(res_placebo, res_tolv_low, res_tolv_high,
                        res_ever, res_combo)
pk_all     <- bind_rows(pk_placebo, pk_tolv_low, pk_tolv_high)

cat("Simulation complete.\n")


## ============================================================
##  4. PUBLICATION-STYLE PLOTS
## ============================================================

cols_scen <- c(
  "1. Placebo"           = "#888888",
  "2. Tolvaptan 60 mg/day" = "#1565C0",
  "3. Tolvaptan 120 mg/day"= "#0D47A1",
  "4. Everolimus 2.5 mg/day" = "#6A1B9A",
  "5. Tolvaptan + ACEi"  = "#00695C"
)

## ---- 4a. TKV Growth Over 3 Years ----
p_tkv <- ggplot(res_all,
                aes(x = time_yr, y = TKV_L, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = cols_scen) +
  labs(
    title    = "Total Kidney Volume (TKV) Over 3 Years",
    subtitle = "ADPKD QSP Model — Primary Disease Endpoint",
    x        = "Time (years)",
    y        = "TKV (L)",
    color    = "Scenario"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 9))

## ---- 4b. eGFR Decline Over 3 Years ----
p_egfr <- ggplot(res_all,
                 aes(x = time_yr, y = eGFR, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = cols_scen) +
  labs(
    title    = "eGFR Decline Over 3 Years",
    subtitle = "CKD Progression in ADPKD",
    x        = "Time (years)",
    y        = "eGFR (mL/min/1.73m²)",
    color    = "Scenario"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 9))

## ---- 4c. Tolvaptan PK (72 hours) ----
p_pk_tolv <- ggplot(pk_all %>% filter(scenario != "1. Placebo"),
                    aes(x = time, y = Cp_tolv, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = cols_scen) +
  labs(
    title = "Tolvaptan Plasma Concentration (0-72 h)",
    x     = "Time (h)",
    y     = "Tolvaptan Cp (ng/mL)",
    color = "Dose"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

## ---- 4d. V2R Occupancy & Urine Osmolality ----
p_v2r <- ggplot(res_tolv_high %>% filter(time_yr <= 0.1),
                aes(x = time / 24)) +
  geom_line(aes(y = V2R_OCC, color = "V2R Occupancy (%)"), linewidth = 1.2) +
  geom_line(aes(y = Uosm / 6,  color = "Urine Osm / 6 (scaled)"), linewidth = 1.2) +
  scale_color_manual(values = c("V2R Occupancy (%)" = "#1565C0",
                                "Urine Osm / 6 (scaled)" = "#E65100")) +
  labs(
    title = "Tolvaptan PD: V2R Occupancy & Urine Osmolality (First Month)",
    x     = "Time (days)",
    y     = "% or Scaled mOsm/kg",
    color = "PD Marker"
  ) +
  theme_bw(base_size = 12)

## ---- 4e. Blood Pressure Over 3 Years ----
p_bp <- ggplot(res_all,
               aes(x = time_yr, y = BP, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = cols_scen) +
  geom_hline(yintercept = 130, linetype = "dashed", color = "red") +
  annotate("text", x = 0.2, y = 131.5, label = "Target: 130 mmHg",
           color = "red", size = 3.5) +
  labs(
    title    = "Blood Pressure Over 3 Years",
    subtitle = "ACEi/ARB effect on RAAS-driven hypertension",
    x        = "Time (years)",
    y        = "Systolic BP (mmHg)",
    color    = "Scenario"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

## ---- 4f. cAMP & mTOR Response ----
p_biomarkers <- res_all %>%
  filter(time_yr %in% c(0, 0.5, 1, 2, 3)) %>%
  select(time_yr, scenario, CAMP_norm, MTOR_norm) %>%
  pivot_longer(cols = c(CAMP_norm, MTOR_norm),
               names_to = "biomarker", values_to = "value") %>%
  ggplot(aes(x = factor(time_yr), y = value, fill = scenario)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~biomarker, scales = "free_y",
             labeller = labeller(biomarker = c(CAMP_norm = "Relative cAMP",
                                               MTOR_norm = "Relative mTOR"))) +
  scale_fill_manual(values = cols_scen) +
  labs(
    title = "cAMP & mTOR Biomarkers by Scenario",
    x     = "Time (years)",
    y     = "Relative Level (1 = ADPKD baseline)",
    fill  = "Scenario"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0))

## ---- 4g. TKV % Change Summary (at 3 years) ----
tkv_summary <- res_all %>%
  group_by(scenario) %>%
  filter(time_yr == max(time_yr)) %>%
  summarise(
    TKV_change_pct = last(TKV_pct),
    eGFR_final     = last(eGFR),
    eGFR_decline   = EGFR0 - last(eGFR),
    BP_final        = last(BP),
    .groups = "drop"
  ) %>%
  rename(Scenario = scenario)

cat("\n=== 3-Year Treatment Outcomes ===\n")
print(tkv_summary %>%
  mutate(
    `TKV Growth (%)` = round(TKV_change_pct, 1),
    `eGFR (mL/min)`  = round(eGFR_final, 1),
    `eGFR Δ (mL/min)`= round(-eGFR_decline, 1),
    `BP (mmHg)`       = round(BP_final, 1)
  ) %>%
  select(Scenario, `TKV Growth (%)`, `eGFR (mL/min)`,
         `eGFR Δ (mL/min)`, `BP (mmHg)`),
  n = 10)


## ============================================================
##  5. DOSE-RESPONSE ANALYSIS: TOLVAPTAN
## ============================================================

doses_tolv <- c(15, 30, 45, 60, 90, 120)  # mg/day total

dose_response <- lapply(doses_tolv, function(d) {
  dose_am <- d * 0.75  # 75% AM dose
  dose_pm <- d * 0.25  # 25% PM dose
  ev_dose  <- ev(cmt = "AGUT", amt = dose_am, time = seq(0, SIM_HOURS, by = 24)) +
              ev(cmt = "AGUT", amt = dose_pm, time = seq(8, SIM_HOURS + 8, by = 24))
  out <- mod %>%
    mrgsim_e(ev_dose, end = SIM_HOURS, delta = 24) %>%
    as.data.frame() %>%
    filter(row_number() == n()) %>%
    mutate(total_dose_mgday = d)
  out
})

dr_df <- bind_rows(dose_response)

p_dose_response <- dr_df %>%
  select(total_dose_mgday, TKV_pct, eGFR) %>%
  pivot_longer(cols = c(TKV_pct, eGFR),
               names_to = "endpoint", values_to = "value") %>%
  ggplot(aes(x = total_dose_mgday, y = value, color = endpoint)) +
  geom_point(size = 3) + geom_line(linewidth = 1.1) +
  facet_wrap(~endpoint, scales = "free_y",
             labeller = labeller(endpoint = c(TKV_pct = "TKV Change at 3yr (%)",
                                              eGFR = "eGFR at 3yr (mL/min)"))) +
  scale_color_manual(values = c(TKV_pct = "#1565C0", eGFR = "#2E7D32")) +
  labs(
    title = "Tolvaptan Dose-Response at 3 Years",
    x     = "Total Tolvaptan Dose (mg/day)",
    y     = "Endpoint Value"
  ) +
  theme_bw(base_size = 12) + theme(legend.position = "none")


## ============================================================
##  6. VIRTUAL PATIENT ANALYSIS (Variability)
## ============================================================

set.seed(2026)
n_vp <- 50  # virtual patients

vp_params <- tibble(
  ID     = 1:n_vp,
  TKV0   = rlnorm(n_vp, meanlog = log(1500), sdlog = 0.4),  # mL
  EGFR0  = rnorm(n_vp, mean = 70, sd = 15) %>% pmax(20),
  KGROW_TKV = rlnorm(n_vp, meanlog = log(6.28e-6), sdlog = 0.3),
  KDECL_EGFR = rlnorm(n_vp, meanlog = log(4e-4), sdlog = 0.3),
  CA2_SCALE  = rnorm(n_vp, mean = 1.5, sd = 0.2) %>% pmax(1)
)

vp_tolvaptan <- lapply(seq_len(n_vp), function(i) {
  p <- as.list(vp_params[i, ])
  mod_vp <- mod %>%
    param(TKV0 = p$TKV0, EGFR0 = p$EGFR0,
          KGROW_TKV = p$KGROW_TKV, KDECL_EGFR = p$KDECL_EGFR,
          CA2_SCALE = p$CA2_SCALE)
  out_tolv <- mod_vp %>%
    mrgsim_e(ev_tolv_low, end = SIM_HOURS, delta = 24 * 30) %>%
    as.data.frame() %>%
    mutate(ID = i, scenario = "Tolvaptan", time_yr = time / (365.25 * 24))
  out_plac <- mod_vp %>%
    mrgsim_e(ev_placebo, end = SIM_HOURS, delta = 24 * 30) %>%
    as.data.frame() %>%
    mutate(ID = i, scenario = "Placebo", time_yr = time / (365.25 * 24))
  bind_rows(out_tolv, out_plac)
})

vp_df <- bind_rows(vp_tolvaptan)

# Calculate median + 5th/95th percentile bands
vp_summary <- vp_df %>%
  group_by(scenario, time_yr) %>%
  summarise(
    TKV_med  = median(TKV_L),
    TKV_lo   = quantile(TKV_L, 0.05),
    TKV_hi   = quantile(TKV_L, 0.95),
    eGFR_med = median(eGFR),
    eGFR_lo  = quantile(eGFR, 0.05),
    eGFR_hi  = quantile(eGFR, 0.95),
    .groups  = "drop"
  )

p_vp_tkv <- ggplot(vp_summary, aes(x = time_yr, color = scenario, fill = scenario)) +
  geom_ribbon(aes(ymin = TKV_lo, ymax = TKV_hi), alpha = 0.2) +
  geom_line(aes(y = TKV_med), linewidth = 1.2) +
  scale_color_manual(values = c("Placebo" = "#888888",
                                "Tolvaptan" = "#1565C0")) +
  scale_fill_manual(values = c("Placebo" = "#888888",
                               "Tolvaptan" = "#1565C0")) +
  labs(
    title    = "TKV Trajectory: Virtual Patient Analysis (n=50)",
    subtitle = "Median + 5th/95th percentile",
    x        = "Time (years)",
    y        = "TKV (L)",
    color    = NULL, fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")


## ============================================================
##  7. PRINT SUMMARY & DISPLAY PLOTS
## ============================================================

cat("\n========== ADPKD QSP MODEL SUMMARY ==========\n")
cat("Model: Tolvaptan, Everolimus, Octreotide, ACEi\n")
cat("Simulation: 3 years, 5 scenarios\n")
cat(sprintf("Tolvaptan (60mg/day) vs Placebo:\n"))
cat(sprintf("  TKV growth: %.1f%% vs %.1f%%\n",
  filter(tkv_summary, Scenario == "2. Tolvaptan 60 mg/day")$TKV_change_pct,
  filter(tkv_summary, Scenario == "1. Placebo")$TKV_change_pct))
cat(sprintf("  eGFR decline: %.1f vs %.1f mL/min\n",
  EGFR0 - filter(tkv_summary, Scenario == "2. Tolvaptan 60 mg/day")$eGFR_final,
  EGFR0 - filter(tkv_summary, Scenario == "1. Placebo")$eGFR_final))

# Arrange all plots
grid.arrange(p_tkv, p_egfr, ncol = 2, top = "ADPKD QSP Model — Treatment Comparison")
grid.arrange(p_pk_tolv, p_v2r,   ncol = 2, top = "ADPKD QSP Model — Tolvaptan PK/PD")
grid.arrange(p_bp, p_biomarkers, ncol = 2, top = "ADPKD QSP Model — Cardiovascular & Biomarkers")
grid.arrange(p_dose_response, p_vp_tkv, ncol = 2, top = "ADPKD QSP Model — Dose-Response & Variability")
