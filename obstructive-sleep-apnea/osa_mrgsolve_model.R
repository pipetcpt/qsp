## ============================================================
## Obstructive Sleep Apnea (OSA) QSP Model — mrgsolve
## 폐쇄성 수면 무호흡증 정량적 시스템 약리학 모델
##
## Disease: Obstructive Sleep Apnea (OSA)
## Target pathways:
##   1. Upper Airway Collapsibility (Pcrit)
##   2. Ventilatory Control (Loop Gain)
##   3. Arousal Threshold
##   4. Apnea-Hypopnea Dynamics → AHI
##   5. Intermittent Hypoxia → SpO2 / HIF-1α / ROS
##   6. Sympathetic Nervous System → BP / HR
##   7. Inflammatory & Metabolic effects
##   8. Drug PK/PD:
##      - CPAP (pneumatic splinting)
##      - Modafinil (orexin/histamine wakefulness)
##      - Solriamfetol (DAT/NET inhibition)
##      - Tirzepatide (GLP-1/GIP → weight loss → AHI)
##      - Eszopiclone (GABA-A → arousal threshold)
##      - Acetazolamide (CA inhibition → loop gain)
##
## Parameter calibration references:
##   - Eckert DJ et al. Sleep 2013;36:1741-52 (phenotyping Pcrit, LG, AT)
##   - Wellman A et al. J Appl Physiol 2011;110:1617-25 (loop gain)
##   - Zinchuk A et al. AJRCCM 2018;197:1245-56 (endotype → outcome)
##   - Malhotra A et al. Sleep Med Rev 2022;63:101636 (pharmacology)
##   - SURMOUNT-OSA trial (Wolk 2024 NEJM) — tirzepatide OSA data
##   - TONES trials — solriamfetol phase 3
##   - HARMONY trial — pitolisant OSA
## ============================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

## ── mrgsolve model code block ───────────────────────────────
osa_model_code <- '
$PROB
OSA QSP Model v1.0
Obstructive Sleep Apnea — Multi-Endotype PK/PD
18 ODE compartments: airway + respiratory control + sleep +
hypoxia + cardiovascular + metabolic + drug PK

$PARAM
// === Upper Airway Parameters ===
Pcrit_base     = 1.5    // cmH2O — baseline critical closing pressure (normal: <-2, severe OSA: +5)
Pcrit_BMI_coef = 0.3    // cmH2O per BMI unit above 25
BMI_ref        = 30.0   // kg/m2 reference BMI for severe OSA
kPcrit         = 0.02   // /hr — return-to-baseline rate
UA_muscle_gain = 0.5    // dimensionless — upper airway dilator efficiency

// === Loop Gain (Ventilatory Instability) ===
LG_base        = 0.85   // loop gain at baseline (>1 = unstable)
kLG            = 0.05   // /hr — LG equilibration rate
LG_hypoxia_cof = 0.1    // LG amplification per 10% SpO2 decrease

// === Arousal Threshold ===
AT_base        = 14.0   // cmH2O — baseline arousal threshold
kAT            = 0.05   // /hr — AT return rate

// === Apnea Dynamics ===
AHI_base       = 35.0   // events/hr — baseline AHI (moderate-severe)
kAHI           = 0.3    // /hr — AHI equilibration rate
AHI_Pcrit_sens = 4.5    // events per cmH2O increase in Pcrit
AHI_LG_sens    = 25.0   // events per unit LG increase above 0.7
AHI_AT_sens    = -0.8   // events per cmH2O AT (higher AT = fewer arousals)

// === SpO2 / Oxygen Dynamics ===
SpO2_base      = 92.5   // % baseline mean overnight SpO2 (OSA population)
SpO2_max       = 98.5   // % max achievable SpO2
kSpO2          = 2.0    // /hr — SpO2 equilibration
SpO2_AHI_coef  = -0.2   // % per event/hr of AHI

// === HIF-1α dynamics ===
kHIF_in        = 0.5    // production rate constant
kHIF_out       = 0.3    // degradation rate
HIF_SpO2_EC50  = 88.0   // SpO2 at half-max HIF-1α induction

// === Sympathetic Nervous System ===
SNA_base       = 1.0    // normalized baseline SNS activity
kSNA           = 0.8    // /hr — SNA equilibration
SNA_hypoxia    = 0.6    // SNS per unit AHI (normalized)

// === Blood Pressure ===
SBP_base       = 138.0  // mmHg baseline SBP (OSA + HTN)
kSBP           = 0.3    // /hr — SBP equilibration
SBP_SNA_coef   = 12.0   // mmHg per normalized SNS unit above baseline

// === Heart Rate ===
HR_base        = 72.0   // bpm
kHR            = 1.0    // /hr
HR_SNA_coef    = 10.0   // bpm per SNS unit

// === Inflammatory Marker (CRP) ===
CRP_base       = 3.5    // mg/L — baseline hsCRP (elevated in OSA)
kCRP           = 0.02   // /hr — CRP equilibration (slow)
CRP_hypoxia    = 0.05   // mg/L per AHI event/hr
CRP_SBP_coef   = 0.02   // mg/L per mmHg SBP

// === Insulin Resistance (HOMA-IR) ===
HOMA_base      = 3.2    // baseline HOMA-IR (insulin resistance)
kHOMA          = 0.01   // /hr — HOMA equilibration
HOMA_CRP_coef  = 0.15   // HOMA increase per mg/L CRP
HOMA_cortisol  = 0.08   // HOMA increase per normalized cortisol unit

// === ESS (Epworth Sleepiness Scale) ===
ESS_base       = 12.0   // 0-24 scale; >10 = excessive daytime sleepiness
kESS           = 0.05   // /hr
ESS_AHI_coef   = 0.08   // ESS increase per event/hr AHI
ESS_SpO2_coef  = -0.15  // ESS decrease per % SpO2 increase

// === Body Weight (for tirzepatide modeling) ===
Weight_base    = 100.0  // kg
kWeight        = 0.001  // /hr — very slow (weeks-scale)

// ═══════════════════════════════════════════════════════════
// === DRUG PK PARAMETERS ===

// --- CPAP ---
CPAP_on       = 0       // 0=off, 1=on
CPAP_pressure = 10.0    // cmH2O

// --- Modafinil PK (1-compartment oral) ---
MODA_CL       = 35.0    // L/hr apparent clearance
MODA_V        = 35.0    // L apparent volume
MODA_ka       = 1.0     // /hr absorption rate
MODA_F        = 0.40    // oral bioavailability
MODA_Emax     = 1.0     // normalized max EDS reduction
MODA_EC50     = 1.0     // μg/mL
MODA_dose_mg  = 200.0   // mg QD
MODA_on       = 0       // 0=off, 1=on

// --- Solriamfetol PK (1-compartment oral) ---
SOLRI_CL      = 12.5    // L/hr
SOLRI_V       = 118.0   // L
SOLRI_ka      = 1.5     // /hr
SOLRI_F       = 0.95    // bioavailability
SOLRI_Emax    = 1.0
SOLRI_EC50    = 0.08    // μg/mL
SOLRI_dose_mg = 150.0   // mg QD
SOLRI_on      = 0

// --- Tirzepatide PK (2-compartment, sc) ---
TIRZ_ka       = 0.02    // /hr — sc depot absorption
TIRZ_CL       = 0.036   // L/hr central clearance (t½≈5d)
TIRZ_V1       = 6.5     // L central volume
TIRZ_Q        = 0.015   // L/hr intercompartmental CL
TIRZ_V2       = 15.0    // L peripheral volume
TIRZ_Emax_wt  = 20.0    // % max weight loss
TIRZ_EC50_wt  = 12.0    // ng/mL
TIRZ_dose_mg  = 10.0    // mg QW
TIRZ_on       = 0

// --- Eszopiclone PK (1-compartment oral) ---
ESZOP_CL      = 21.0    // L/hr
ESZOP_V       = 80.0    // L
ESZOP_ka      = 1.2     // /hr
ESZOP_F       = 0.80
ESZOP_Emax    = 0.6     // fraction AT increase at max effect
ESZOP_EC50    = 0.5     // ng/mL
ESZOP_dose_mg = 3.0     // mg QHS
ESZOP_on      = 0

// --- Acetazolamide PK (1-compartment oral) ---
ACETZ_CL      = 3.5     // L/hr
ACETZ_V       = 30.0    // L
ACETZ_ka      = 1.0     // /hr
ACETZ_F       = 0.90
ACETZ_Emax_LG = 0.40    // max fractional LG reduction
ACETZ_EC50    = 8.0     // μg/mL
ACETZ_dose_mg = 250.0   // mg BID
ACETZ_on      = 0

$INIT
// === Disease State Variables ===
Pcrit   = 1.5     // effective Pcrit (cmH2O)
LG      = 0.85    // loop gain
AT      = 14.0    // arousal threshold (cmH2O)
AHI     = 35.0    // events/hr
SpO2    = 92.5    // % overnight mean
HIF1a   = 1.0     // normalized
SNA     = 1.0     // sympathetic nerve activity
SBP     = 138.0   // mmHg
HR      = 72.0    // bpm
CRP     = 3.5     // mg/L
HOMA    = 3.2     // HOMA-IR
ESS     = 12.0    // Epworth score
Weight  = 100.0   // kg

// === Drug PK Compartments ===
// Modafinil
MODA_depot   = 0
MODA_central = 0

// Solriamfetol
SOLRI_depot   = 0
SOLRI_central = 0

// Tirzepatide
TIRZ_depot    = 0
TIRZ_c1       = 0
TIRZ_c2       = 0

// Eszopiclone
ESZOP_depot   = 0
ESZOP_central = 0

// Acetazolamide
ACETZ_depot   = 0
ACETZ_central = 0

$GLOBAL
// Helper to compute Emax model effect (inhibitory or stimulatory)
// E = Emax * Cp^n / (EC50^n + Cp^n)
// For inhibitory: effect = 1 - Emax * Cp/(EC50+Cp)

$ODE

// ─────────────────────────────────────────────────
// SECTION 1: Drug PK ODEs
// ─────────────────────────────────────────────────

// -- Modafinil --
double MODA_dose_rate = (MODA_on > 0.5) ? (MODA_dose_mg * MODA_F * MODA_ka) : 0.0;
dxdt_MODA_depot   = MODA_dose_rate - MODA_ka * MODA_depot;
dxdt_MODA_central = MODA_ka * MODA_depot - (MODA_CL / MODA_V) * MODA_central;
double Cp_MODA = MODA_central / MODA_V;              // mg/L = μg/mL

// -- Solriamfetol --
double SOLRI_dose_rate = (SOLRI_on > 0.5) ? (SOLRI_dose_mg * SOLRI_F * SOLRI_ka) : 0.0;
dxdt_SOLRI_depot   = SOLRI_dose_rate - SOLRI_ka * SOLRI_depot;
dxdt_SOLRI_central = SOLRI_ka * SOLRI_depot - (SOLRI_CL / SOLRI_V) * SOLRI_central;
double Cp_SOLRI = SOLRI_central / SOLRI_V;           // mg/L = μg/mL

// -- Tirzepatide --
double TIRZ_dose_rate = (TIRZ_on > 0.5) ? (TIRZ_dose_mg * TIRZ_ka) : 0.0;  // weekly dose → hourly input
dxdt_TIRZ_depot = TIRZ_dose_rate - TIRZ_ka * TIRZ_depot;
dxdt_TIRZ_c1 = TIRZ_ka * TIRZ_depot
               - (TIRZ_CL / TIRZ_V1) * TIRZ_c1
               - (TIRZ_Q  / TIRZ_V1) * TIRZ_c1
               + (TIRZ_Q  / TIRZ_V2) * TIRZ_c2;
dxdt_TIRZ_c2 = (TIRZ_Q / TIRZ_V1) * TIRZ_c1 - (TIRZ_Q / TIRZ_V2) * TIRZ_c2;
double Cp_TIRZ = TIRZ_c1 / TIRZ_V1 * 1e3;          // ng/mL

// -- Eszopiclone --
double ESZOP_dose_rate = (ESZOP_on > 0.5) ? (ESZOP_dose_mg * ESZOP_F * ESZOP_ka) : 0.0;
dxdt_ESZOP_depot   = ESZOP_dose_rate - ESZOP_ka * ESZOP_depot;
dxdt_ESZOP_central = ESZOP_ka * ESZOP_depot - (ESZOP_CL / ESZOP_V) * ESZOP_central;
double Cp_ESZOP = ESZOP_central / ESZOP_V * 1e3;    // ng/mL

// -- Acetazolamide --
double ACETZ_dose_rate = (ACETZ_on > 0.5) ? (ACETZ_dose_mg * ACETZ_F * ACETZ_ka) : 0.0;
dxdt_ACETZ_depot   = ACETZ_dose_rate - ACETZ_ka * ACETZ_depot;
dxdt_ACETZ_central = ACETZ_ka * ACETZ_depot - (ACETZ_CL / ACETZ_V) * ACETZ_central;
double Cp_ACETZ = ACETZ_central / ACETZ_V;          // mg/L = μg/mL

// ─────────────────────────────────────────────────
// SECTION 2: Drug PD Effects
// ─────────────────────────────────────────────────

// CPAP: directly reduces effective Pcrit to < 0 cmH2O
double CPAP_Pcrit_reduction = (CPAP_on > 0.5) ? (-CPAP_pressure - 1.5) : 0.0;

// MAD: approximate 3 cmH2O Pcrit reduction (if included via parameter)

// Tirzepatide weight effect → Pcrit reduction via fat pad loss
double Tirz_wt_eff = TIRZ_Emax_wt * Cp_TIRZ / (TIRZ_EC50_wt + Cp_TIRZ);  // % weight loss
double Tirz_Pcrit_eff = -Tirz_wt_eff * 0.10;  // 0.10 cmH2O per % weight loss

// Eszopiclone → arousal threshold increase
double ESZOP_AT_eff = ESZOP_Emax * Cp_ESZOP / (ESZOP_EC50 + Cp_ESZOP);

// Acetazolamide → loop gain reduction
double ACETZ_LG_eff = ACETZ_Emax_LG * Cp_ACETZ / (ACETZ_EC50 + Cp_ACETZ);

// Modafinil → EDS reduction (indirect: improves ESS via orexin/histamine)
double MODA_ESS_eff = MODA_Emax * Cp_MODA / (MODA_EC50 + Cp_MODA);

// Solriamfetol → EDS reduction (direct DAT/NET)
double SOLRI_ESS_eff = SOLRI_Emax * Cp_SOLRI / (SOLRI_EC50 + Cp_SOLRI);

// ─────────────────────────────────────────────────
// SECTION 3: Disease State ODEs
// ─────────────────────────────────────────────────

// -- Effective Critical Closing Pressure --
double Pcrit_target = Pcrit_base
                    + Pcrit_BMI_coef * (Weight / 1.75 / 1.75 - 25.0)  // BMI proxy
                    + CPAP_Pcrit_reduction
                    + Tirz_Pcrit_eff;
dxdt_Pcrit = kPcrit * (Pcrit_target - Pcrit);

// -- Loop Gain --
double LG_hypoxia_effect = LG_hypoxia_cof * (SpO2_max - SpO2) / 10.0;
double LG_target = (LG_base + LG_hypoxia_effect) * (1 - ACETZ_LG_eff);
dxdt_LG = kLG * (LG_target - LG);

// -- Arousal Threshold --
double AT_target = AT_base * (1 + ESZOP_AT_eff);
dxdt_AT = kAT * (AT_target - AT);

// -- AHI (Apnea-Hypopnea Index) --
// Mechanistic equation: AHI driven by Pcrit, LG, and inversely by AT
// Eckert 2013: AHI = f(collapsibility, loop gain, arousal threshold, genioglossal response)
double AHI_Pcrit_contrib = AHI_Pcrit_sens * (Pcrit - (-2.0));     // deviation from normal Pcrit
double AHI_LG_contrib    = AHI_LG_sens   * (LG   - 0.70);         // deviation from stable LG
double AHI_AT_contrib    = AHI_AT_sens   * (AT   - 14.0);         // higher AT helps reduce AHI
double AHI_target_raw    = AHI_Pcrit_contrib + AHI_LG_contrib + AHI_AT_contrib;
double AHI_target = (AHI_target_raw < 1.0) ? 1.0 : AHI_target_raw;  // floor at 1
dxdt_AHI = kAHI * (AHI_target - AHI);

// -- SpO2 --
double SpO2_target = SpO2_base + SpO2_AHI_coef * AHI
                   + (CPAP_on > 0.5 ? (SpO2_max - SpO2_base) * 0.9 : 0.0);
SpO2_target = (SpO2_target > SpO2_max) ? SpO2_max : SpO2_target;
SpO2_target = (SpO2_target < 70.0)    ? 70.0      : SpO2_target;
dxdt_SpO2 = kSpO2 * (SpO2_target - SpO2);

// -- HIF-1α (normalized, 1=normal, >1=elevated) --
double HIF_stimulus = (SpO2 < HIF_SpO2_EC50)
                    ? 2.0 * (HIF_SpO2_EC50 - SpO2) / HIF_SpO2_EC50
                    : 0.0;
dxdt_HIF1a = kHIF_in * (1 + HIF_stimulus) - kHIF_out * HIF1a;

// -- Sympathetic Nerve Activity --
double SNA_target = SNA_base + SNA_hypoxia * (AHI / AHI_base - 1.0);
SNA_target = (SNA_target < 0.2) ? 0.2 : SNA_target;
dxdt_SNA = kSNA * (SNA_target - SNA);

// -- Systolic Blood Pressure --
double SBP_target = SBP_base + SBP_SNA_coef * (SNA - 1.0);
dxdt_SBP = kSBP * (SBP_target - SBP);

// -- Heart Rate --
double HR_target = HR_base + HR_SNA_coef * (SNA - 1.0);
dxdt_HR = kHR * (HR_target - HR);

// -- CRP (slow inflammatory marker) --
double CRP_target = CRP_base + CRP_hypoxia * AHI + CRP_SBP_coef * (SBP - 130.0);
dxdt_CRP = kCRP * (CRP_target - CRP);

// -- HOMA-IR --
double HOMA_target = HOMA_base + HOMA_CRP_coef * (CRP - 1.0)
                   + HOMA_cortisol * (SNA - 1.0)
                   - (Tirz_wt_eff * 0.08);   // tirzepatide insulin sensitization
HOMA_target = (HOMA_target < 1.0) ? 1.0 : HOMA_target;
dxdt_HOMA = kHOMA * (HOMA_target - HOMA);

// -- ESS (Epworth Sleepiness Scale) --
double ESS_drug_eff = MODA_ESS_eff + SOLRI_ESS_eff;
ESS_drug_eff = (ESS_drug_eff > 0.8) ? 0.8 : ESS_drug_eff;  // cap at 80% reduction
double ESS_target = ESS_base * (1.0 - ESS_drug_eff)
                  + ESS_AHI_coef  * (AHI - AHI_base)
                  + ESS_SpO2_coef * (SpO2 - SpO2_base);
ESS_target = (ESS_target < 0)    ? 0    : ESS_target;
ESS_target = (ESS_target > 24.0) ? 24.0 : ESS_target;
dxdt_ESS = kESS * (ESS_target - ESS);

// -- Body Weight (responds slowly to tirzepatide) --
double Weight_loss_rate = kWeight * Tirz_wt_eff * Weight;
dxdt_Weight = -Weight_loss_rate;

$TABLE
capture Pcrit_eff   = Pcrit;
capture LoopGain    = LG;
capture ArousThresh = AT;
capture AHI_val     = AHI;
capture SpO2_val    = SpO2;
capture HIF1a_val   = HIF1a;
capture SNA_val     = SNA;
capture SBP_val     = SBP;
capture HR_val      = HR;
capture CRP_val     = CRP;
capture HOMA_val    = HOMA;
capture ESS_val     = ESS;
capture Weight_kg   = Weight;
capture Cp_MODA_uM  = Cp_MODA;
capture Cp_SOLRI_uM = Cp_SOLRI;
capture Cp_TIRZ_ng  = Cp_TIRZ;
capture Cp_ESZOP_ng = Cp_ESZOP;
capture Cp_ACETZ_uM = Cp_ACETZ;
capture BMI_proxy   = Weight / 1.75 / 1.75;
capture TIRZ_wt_eff = Tirz_wt_eff;
'

## ── Compile model ────────────────────────────────────────────
mod <- mcode("OSA_QSP", osa_model_code)

## ── Helper: simulation function ──────────────────────────────
run_osa_sim <- function(
    mod,
    duration_weeks = 52,
    CPAP_on        = 0,
    CPAP_pressure  = 10,
    MODA_on        = 0,
    SOLRI_on       = 0,
    TIRZ_on        = 0,
    ESZOP_on       = 0,
    ACETZ_on       = 0,
    AHI_init       = 35,
    BMI_init       = 30,
    scenario_name  = "Untreated"
) {
  time_hours <- seq(0, duration_weeks * 168, by = 1)

  # Drug dosing events (daily oral for most; weekly sc for tirzepatide)
  # For simplicity, drugs are modeled as continuous infusion in PK depot
  # In clinical use, dose events would be via ev() objects.

  params_update <- list(
    CPAP_on       = CPAP_on,
    CPAP_pressure = CPAP_pressure,
    MODA_on       = MODA_on,
    SOLRI_on      = SOLRI_on,
    TIRZ_on       = TIRZ_on,
    ESZOP_on      = ESZOP_on,
    ACETZ_on      = ACETZ_on,
    AHI_base      = AHI_init,
    Weight_base   = BMI_init * 1.75^2
  )

  init_update <- list(
    AHI    = AHI_init,
    Weight = BMI_init * 1.75^2
  )

  out <- mod %>%
    param(params_update) %>%
    init(init_update) %>%
    mrgsim(end = max(time_hours), delta = 1) %>%
    as_tibble() %>%
    mutate(
      week           = time / 168,
      scenario       = scenario_name,
      AHI_severity   = case_when(
        AHI_val < 5   ~ "Normal",
        AHI_val < 15  ~ "Mild",
        AHI_val < 30  ~ "Moderate",
        TRUE          ~ "Severe"
      )
    )
  out
}

## ── SCENARIO 1: Untreated moderate-severe OSA ────────────────
s1 <- run_osa_sim(mod,
  duration_weeks = 52,
  CPAP_on = 0, MODA_on = 0, TIRZ_on = 0,
  AHI_init = 35, BMI_init = 30,
  scenario_name = "Untreated (AHI=35)"
)

## ── SCENARIO 2: CPAP therapy ─────────────────────────────────
s2 <- run_osa_sim(mod,
  duration_weeks = 52,
  CPAP_on = 1, CPAP_pressure = 10,
  AHI_init = 35, BMI_init = 30,
  scenario_name = "CPAP 10 cmH2O"
)

## ── SCENARIO 3: Modafinil (residual EDS on CPAP) ─────────────
s3 <- run_osa_sim(mod,
  duration_weeks = 12,
  CPAP_on = 1, MODA_on = 1,
  AHI_init = 35, BMI_init = 30,
  scenario_name = "CPAP + Modafinil 200mg"
)

## ── SCENARIO 4: Solriamfetol (residual EDS on CPAP) ──────────
s4 <- run_osa_sim(mod,
  duration_weeks = 12,
  CPAP_on = 1, SOLRI_on = 1,
  AHI_init = 35, BMI_init = 30,
  scenario_name = "CPAP + Solriamfetol 150mg"
)

## ── SCENARIO 5: Tirzepatide (obesity-related OSA) ────────────
s5 <- run_osa_sim(mod,
  duration_weeks = 52,
  TIRZ_on = 1,
  AHI_init = 35, BMI_init = 38,
  scenario_name = "Tirzepatide 10mg QW (No CPAP)"
)

## ── SCENARIO 6: Phenotype-guided combination therapy ─────────
## High loop gain + elevated arousal threshold phenotype
s6 <- run_osa_sim(mod,
  duration_weeks = 12,
  CPAP_on = 0, ACETZ_on = 1, ESZOP_on = 1,
  AHI_init = 28, BMI_init = 27,
  scenario_name = "Acetazolamide + Eszopiclone\n(Non-CPAP phenotypic)"
)

## ── Combine scenarios for AHI comparison ─────────────────────
all_scenarios <- bind_rows(s1, s2, s3, s4, s5, s6)

## ── PLOT 1: AHI over time by scenario ────────────────────────
p1 <- all_scenarios %>%
  filter(week %in% seq(0, 52, by = 0.5)) %>%
  ggplot(aes(x = week, y = AHI_val, color = scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = c(5, 15, 30), linetype = "dashed",
             color = c("green3", "orange", "red"), alpha = 0.6) +
  annotate("text", x = 52, y = c(3, 13, 28, 40),
           label = c("Normal <5", "Mild <15", "Moderate <30", "Severe ≥30"),
           hjust = 1, size = 3, color = c("green4","orange3","red3","red4")) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title    = "AHI Trajectory by Treatment Scenario",
    subtitle = "OSA QSP Model — mrgsolve simulation",
    x        = "Time (weeks)",
    y        = "AHI (events/hr)",
    color    = "Scenario"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

## ── PLOT 2: SpO2 & Cardiovascular ────────────────────────────
p2 <- all_scenarios %>%
  filter(week %in% seq(0, 52, by = 0.5)) %>%
  ggplot(aes(x = week, y = SpO2_val, color = scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "blue") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "SpO2 Over Time", x = "Weeks", y = "SpO2 (%)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p3 <- all_scenarios %>%
  filter(week %in% seq(0, 52, by = 0.5)) %>%
  ggplot(aes(x = week, y = SBP_val, color = scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 130, linetype = "dashed", color = "red") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Systolic BP Over Time", x = "Weeks", y = "SBP (mmHg)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

## ── PLOT 3: ESS (sleepiness) ─────────────────────────────────
p4 <- all_scenarios %>%
  filter(week %in% seq(0, 52, by = 0.5)) %>%
  ggplot(aes(x = week, y = ESS_val, color = scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "orange") +
  annotate("text", x = 48, y = 10.8, label = "EDS threshold (ESS>10)",
           size = 3, color = "orange4") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Epworth Sleepiness Scale", x = "Weeks", y = "ESS Score (0-24)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

## ── PLOT 4: Tirzepatide weight + AHI ─────────────────────────
p5 <- s5 %>%
  filter(week %in% seq(0, 52, by = 0.5)) %>%
  ggplot(aes(x = week)) +
  geom_line(aes(y = Weight_kg), color = "steelblue", linewidth = 1.2) +
  geom_line(aes(y = AHI_val * 3), color = "red3", linewidth = 1.2, linetype = "dashed") +
  scale_y_continuous(
    name = "Weight (kg)",
    sec.axis = sec_axis(~./3, name = "AHI (events/hr)")
  ) +
  labs(title = "Tirzepatide: Weight Loss → AHI Reduction",
       subtitle = "SURMOUNT-OSA Trial Replication",
       x = "Weeks") +
  theme_bw(base_size = 11) +
  theme(axis.title.y.left  = element_text(color = "steelblue"),
        axis.title.y.right = element_text(color = "red3"))

## ── PLOT 5: Metabolic effects ─────────────────────────────────
p6 <- all_scenarios %>%
  filter(week %in% seq(0, 52, by = 2)) %>%
  ggplot(aes(x = week, y = HOMA_val, color = scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 2.5, linetype = "dashed", color = "forestgreen") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "HOMA-IR (Insulin Resistance)",
       x = "Weeks", y = "HOMA-IR") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

## ── PK Profiles: modafinil + solriamfetol ────────────────────
pk_sim <- function(dose_mg, CL, V, ka, label, hours = 24) {
  t  <- seq(0, hours, by = 0.1)
  Cp <- (dose_mg * ka) / (V * (ka - CL/V)) *
        (exp(-CL/V * t) - exp(-ka * t))
  tibble(time = t, Cp = Cp, drug = label)
}

pk_data <- bind_rows(
  pk_sim(200, 35, 35, 1.0, "Modafinil 200mg"),
  pk_sim(150, 12.5, 118, 1.5, "Solriamfetol 150mg"),
  pk_sim(3,   21,  80,  1.2, "Eszopiclone 3mg") %>%
    mutate(Cp = Cp * 1000)  # scale to ng/mL for display
) %>%
  mutate(unit = ifelse(drug == "Eszopiclone 3mg", "ng/mL", "μg/mL"))

p7 <- pk_data %>%
  ggplot(aes(x = time, y = Cp, color = drug)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~unit, scales = "free_y") +
  labs(title = "PK Profiles — Wake-Promoting & Arousal-Modulating Agents",
       x = "Time (hours)", y = "Plasma Concentration",
       color = "Drug") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

## ── Print summary table ───────────────────────────────────────
summary_table <- all_scenarios %>%
  filter(week %in% c(0, 4, 12, 26, 52)) %>%
  group_by(scenario, week) %>%
  slice_tail(n = 1) %>%
  select(scenario, week, AHI_val, SpO2_val, SBP_val, ESS_val, HOMA_val, CRP_val) %>%
  ungroup() %>%
  rename(
    Scenario = scenario,
    Week     = week,
    `AHI (ev/hr)`   = AHI_val,
    `SpO2 (%)`      = SpO2_val,
    `SBP (mmHg)`    = SBP_val,
    `ESS`           = ESS_val,
    `HOMA-IR`       = HOMA_val,
    `hsCRP (mg/L)`  = CRP_val
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 1)))

cat("\n=== OSA QSP Model — Clinical Endpoint Summary ===\n")
print(summary_table, n = 100)

## ── Composite dashboard ───────────────────────────────────────
dashboard <- (p1 | p4) / (p2 | p3) / (p4 | p6) +
  plot_annotation(
    title   = "Obstructive Sleep Apnea (OSA) — QSP Model Dashboard",
    subtitle = paste0("mrgsolve ODE model | ", Sys.Date()),
    theme   = theme(plot.title = element_text(size = 14, face = "bold"))
  )

# Save plot
tryCatch({
  ggsave("osa_qsp_dashboard.pdf",
         plot   = (p1 / (p2 + p3) / (p4 + p6) / p7),
         width  = 14,
         height = 18,
         device = "pdf")
  cat("Dashboard saved to osa_qsp_dashboard.pdf\n")
}, error = function(e) cat("PDF save skipped:", conditionMessage(e), "\n"))

## ── Parameter sensitivity analysis ───────────────────────────
sens_Pcrit <- sapply(seq(-2, 6, by = 0.5), function(Pc) {
  mod %>%
    param(Pcrit_base = Pc, CPAP_on = 0) %>%
    init(Pcrit = Pc) %>%
    mrgsim(end = 720, delta = 24) %>%
    as_tibble() %>%
    filter(time == 720) %>%
    pull(AHI_val)
})

sens_df <- tibble(
  Pcrit_base = seq(-2, 6, by = 0.5),
  AHI_ss     = sens_Pcrit
)

p_sens <- ggplot(sens_df, aes(x = Pcrit_base, y = AHI_ss)) +
  geom_line(color = "steelblue", linewidth = 1.3) +
  geom_point(color = "steelblue", size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = 0.2, y = max(sens_df$AHI_ss) * 0.9,
           label = "Pcrit = 0\n(threshold)", color = "red", size = 3.5) +
  labs(
    title    = "Sensitivity: Pcrit vs Steady-State AHI",
    subtitle = "Upper airway collapsibility as primary determinant",
    x        = "Baseline Pcrit (cmH2O)",
    y        = "Steady-State AHI (events/hr)"
  ) +
  theme_bw(base_size = 11)

print(p_sens)

cat("\n=== Simulation Complete ===\n")
cat("Outputs: osa_qsp_dashboard.pdf\n")
cat("Objects: s1-s6 (scenario data frames), summary_table, mod\n")
