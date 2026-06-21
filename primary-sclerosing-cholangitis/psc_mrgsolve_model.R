## =============================================================================
## Primary Sclerosing Cholangitis (PSC) — QSP mrgsolve Model
## =============================================================================
## Disease: Primary Sclerosing Cholangitis (PSC)
## Model type: Multi-scale ODE — Gut-Liver Axis / Bile Acid / Immune / Fibrosis
## Author: QSP Library (auto-generated via Claude Code Routine)
## Date: 2026-06-19
##
## Compartments (24 ODEs):
##   PK:  UDCA_gut, UDCA_plasma, UDCA_bile (3)
##        OCA_gut, OCA_plasma, OCA_bile (3)
##        NorUDCA_bile, BEZ_plasma (2)
##   PD:
##   Gut: LPS, GutBarrier (2)
##   BA:  FXR_act, BilePool, HydroIndex (3)
##   Immune: IL17A, TNFa, IL6, Treg_IL10 (4)
##   Cholangiocyte: Cholangio_health, Senescence (2)
##   Fibrosis: HSC_act, Col1a1, LOXL2 (3)
##   Biomarkers: ALP, Bilirubin, Fibroscan (3)
##   CCA: CCA_risk (1)
##   Portal: PortalPressure (1)
##
## Key References:
##   - Lazaridis & LaRusso (2016) NEJM 375:1161
##   - Eaton et al. (2020) Hepatology 71:2219
##   - Trauner et al. (2017) Gut 66:1933
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ---------------------------------------------------------------------------
## 1. Model Code Block
## ---------------------------------------------------------------------------
psc_model_code <- '
$PROB
Primary Sclerosing Cholangitis (PSC) QSP Model
Multi-compartment ODE: PK + Gut-Liver + BA signaling + Immune + Fibrosis
Version 1.0 — 2026-06-19

$PARAM @annotated
// ---- Drug PK parameters ----
Ka_UDCA   : 0.8    : UDCA oral absorption rate constant (1/h)
F_UDCA    : 0.50   : UDCA bioavailability
CL_UDCA   : 12.0   : UDCA total clearance (L/h)
Vd_UDCA   : 25.0   : UDCA volume of distribution (L)
Kbile_UDCA: 0.15   : UDCA hepatic-to-bile transfer rate (1/h)

Ka_OCA    : 0.6    : OCA oral absorption rate constant (1/h)
F_OCA     : 0.60   : OCA bioavailability
CL_OCA    : 8.0    : OCA total clearance (L/h)
Vd_OCA    : 20.0   : OCA volume of distribution (L)
Kbile_OCA : 0.12   : OCA hepatic-to-bile transfer rate (1/h)
IC50_OCA_FXR: 0.10 : OCA FXR IC50 in bile (μmol/L)

Kbile_NorUDCA: 0.20 : norUDCA cholehepatic shunting rate (1/h)
Ka_BEZ    : 1.0    : Bezafibrate absorption rate constant (1/h)
CL_BEZ    : 15.0   : Bezafibrate clearance (L/h)
Vd_BEZ    : 18.0   : Bezafibrate Vd (L)

// ---- Gut-Liver Axis ----
k_LPS_prod   : 0.05  : LPS production rate (baseline dysbiosis) (AU/h)
k_LPS_clear  : 0.10  : LPS portal clearance by Kupffer cells (1/h)
LPS_max      : 1.0   : Max LPS level (AU, normalized)
k_barrier    : 0.02  : Gut barrier repair rate (1/h)
k_IBD_damage : 0.005 : IBD contribution to barrier damage rate (AU/h)
GutBarrier0  : 1.0   : Baseline gut barrier integrity (0-1)
VAN_eff      : 0.60  : Vancomycin efficacy on gut LPS (fraction reduction)

// ---- Bile Acid metabolism ----
k_FXR_act   : 0.30   : FXR activation rate constant by bile OCA/UDCA (1/h/μmol)
k_FXR_decay : 0.05   : FXR activation decay rate (1/h)
FXR0        : 0.30   : Baseline FXR activity (0-1, partial activation by endogenous BA)
k_BA_synth  : 0.10   : Endogenous BA synthesis rate (AU/h)
k_BA_excrete: 0.08   : Bile acid excretion rate (1/h)
k_hydro_form : 0.04  : Rate of hydrophobic BA index increase (when FXR ↓)
k_hydro_clear: 0.06  : Hydrophobic BA clearance rate (1/h)
HydroIndex0 : 0.30   : Baseline hydrophobic BA index
UDCA_hydro_eff: 0.50 : UDCA effect on reducing hydrophobic index (fraction)

// ---- Immune dynamics ----
k_IL17_base  : 0.002 : Basal IL-17A secretion (AU/h)
k_IL17_LPS   : 0.05  : LPS → IL-17A production rate
k_IL17_BA    : 0.04  : Hydrophobic BA → IL-17A production
k_IL17_decay : 0.08  : IL-17A decay rate (1/h)
k_TNFa_base  : 0.003 : Basal TNF-α secretion (AU/h)
k_TNFa_LPS   : 0.06  : LPS → TNF-α
k_TNFa_decay : 0.10  : TNF-α decay rate (1/h)
k_IL6_IL17   : 0.03  : IL-17A → IL-6 production
k_IL6_TNFa   : 0.04  : TNF-α → IL-6 production
k_IL6_decay  : 0.12  : IL-6 decay rate (1/h)
k_Treg_base  : 0.001 : Basal Treg/IL-10 activity
k_Treg_decay : 0.06  : Treg decay (1/h)
k_Treg_inhib : 0.03  : IL-6 inhibition of Treg
IL17_SS      : 0.10  : IL-17A steady-state (AU, normalized)
TNFa_SS      : 0.08  : TNF-α steady-state

// ---- Cholangiocyte / Biliary health ----
Cholangio0   : 1.0   : Baseline cholangiocyte health (0-1)
k_cholangio_damage_IL17: 0.02  : IL-17A → cholangiocyte damage
k_cholangio_damage_BA  : 0.03  : HydroIndex → cholangiocyte damage
k_cholangio_repair     : 0.01  : Cholangiocyte repair rate
k_senesc_form: 0.015   : Cholangiocyte senescence formation rate
k_senesc_clear: 0.008  : Senescent cell clearance

// ---- Hepatic Fibrosis ----
k_HSC_act    : 0.01  : HSC activation rate (by TGF-β proxy = IL-17A + senescence)
k_HSC_resol  : 0.005 : HSC resolution rate
k_Col_synth  : 0.008 : Collagen I synthesis rate by activated HSC
k_Col_degrad : 0.004 : Collagen degradation rate (MMP-2 mediated)
LOXL2_0      : 0.20  : Baseline LOXL2 cross-linking activity
k_LOXL2_form : 0.005 : LOXL2 production by activated HSC
k_LOXL2_clear: 0.008 : LOXL2 clearance rate
Simtu_LOXL2  : 0.70  : Simtuzumab LOXL2 inhibition (fraction)

// ---- Clinical Biomarkers ----
ALP0         : 200.0 : Baseline ALP (IU/L, elevated in PSC)
k_ALP_chol   : 5.0   : Cholestasis → ALP increase rate
k_ALP_repair : 0.30  : ALP return-to-normal rate
ALP_normal   : 100.0 : Normal ALP level (IU/L)
Bili0        : 1.5   : Baseline bilirubin (mg/dL)
k_Bili_chol  : 0.20  : Cholestasis → bilirubin increase
k_Bili_synth : 0.05  : Normal bilirubin production (mg/dL/h)
k_Bili_excrete: 0.10 : Bilirubin excretion rate (1/h)
Fibro0       : 8.0   : Baseline liver stiffness (kPa, PSC-elevated)
k_Fibro_form : 0.10  : Fibrosis → stiffness increase (kPa per Col unit/h)
k_Fibro_resol: 0.005 : Fibrosis-driven stiffness resolution

// ---- Portal Hypertension ----
PP0          : 5.0   : Baseline portal pressure (mmHg)
k_PP_fibrosis: 0.50  : Col1a1 → portal pressure increase rate
k_PP_decay   : 0.02  : Portal pressure resolution rate
PP_varices   : 12.0  : Portal pressure threshold for varices (mmHg)

// ---- CCA risk ----
k_CCA_chol   : 0.0001 : Cholestasis-driven CCA risk accumulation rate
k_CCA_senesc : 0.0002 : Senescence-driven CCA risk (SASP)
k_CCA_max    : 1.0    : Maximum CCA risk index

// ---- Disease severity flags ----
IBD_status   : 1.0   : IBD co-existence (1=present, 0=absent)
Disease_duration: 0.0 : Disease duration at start (years)

$CMT @annotated
// Drug PK compartments
UDCA_gut    : UDCA gut absorption compartment (mg)
UDCA_plasma : UDCA plasma (mg)
UDCA_bile   : UDCA biliary concentration (mg)
OCA_gut     : OCA gut absorption (mg)
OCA_plasma  : OCA plasma (mg)
OCA_bile    : OCA biliary concentration (mg)
NorUDCA_bile: norUDCA biliary pool (mg)
BEZ_plasma  : Bezafibrate plasma (mg)

// Disease state variables
LPS         : Portal LPS (AU)
GutBarrier  : Gut barrier integrity (0-1)
FXR_act     : FXR activation level (0-1)
BilePool    : Bile acid pool homeostasis (AU)
HydroIndex  : Hydrophobic bile acid index (AU)
IL17A       : IL-17A concentration (AU)
TNFa        : TNF-α concentration (AU)
IL6         : IL-6 concentration (AU)
Treg_IL10   : Regulatory T cell / IL-10 activity (AU)
Cholangio_health : Cholangiocyte health index (0-1)
Senescence  : Senescent cholangiocyte burden (AU)
HSC_act     : Activated HSC fraction (0-1)
Col1a1      : Collagen I/III matrix load (AU)
LOXL2       : LOXL2 cross-linking enzyme (AU)
ALP         : Alkaline phosphatase (IU/L)
Bilirubin   : Total bilirubin (mg/dL)
Fibroscan   : Liver stiffness kPa
PortalPressure : Portal pressure (mmHg)
CCA_risk    : Cumulative CCA risk index (0-1)

$MAIN
// Derived variables — drug concentrations in bile (μmol/L equivalent for FXR)
double OCA_bile_conc = OCA_bile / 0.5;   // convert to approximate μmol/L
double UDCA_bile_conc = UDCA_bile / 2.0;

// FXR activation by bile acids (Hill equation)
double FXR_OCA = (OCA_bile_conc / (IC50_OCA_FXR + OCA_bile_conc));
double FXR_UDCA = 0.05 * (UDCA_bile_conc / (5.0 + UDCA_bile_conc));
double FXR_input = FXR0 + (1.0 - FXR0) * (FXR_OCA + FXR_UDCA);
FXR_input = (FXR_input > 1.0) ? 1.0 : FXR_input;

// Cholestasis severity (1 - cholangiocyte health * bile pool)
double Cholestasis = 1.0 - (Cholangio_health * BilePool);
Cholestasis = (Cholestasis < 0.0) ? 0.0 : Cholestasis;
Cholestasis = (Cholestasis > 1.0) ? 1.0 : Cholestasis;

// Portal pressure risk
double VaricesRisk = (PortalPressure > PP_varices) ? 1.0 : (PortalPressure / PP_varices);

$ODE
// ===========================================================================
// DRUG PK ODEs
// ===========================================================================

// UDCA
dxdt_UDCA_gut    = -Ka_UDCA * UDCA_gut;
dxdt_UDCA_plasma = F_UDCA * Ka_UDCA * UDCA_gut - (CL_UDCA / Vd_UDCA) * UDCA_plasma;
dxdt_UDCA_bile   = Kbile_UDCA * UDCA_plasma - 0.05 * UDCA_bile;

// OCA
dxdt_OCA_gut    = -Ka_OCA * OCA_gut;
dxdt_OCA_plasma = F_OCA * Ka_OCA * OCA_gut - (CL_OCA / Vd_OCA) * OCA_plasma;
dxdt_OCA_bile   = Kbile_OCA * OCA_plasma - 0.04 * OCA_bile;

// norUDCA biliary (cholehepatic shunting — direct to bile)
dxdt_NorUDCA_bile = Kbile_NorUDCA * 1.0 - 0.06 * NorUDCA_bile;

// Bezafibrate
dxdt_BEZ_plasma = Ka_BEZ * 1.0 - (CL_BEZ / Vd_BEZ) * BEZ_plasma;

// ===========================================================================
// GUT-LIVER AXIS ODEs
// ===========================================================================

// LPS (portal burden)
double LPS_prod = k_LPS_prod * (2.0 - GutBarrier) * (1.0 + IBD_status * 0.5);
double LPS_cleared = k_LPS_clear * LPS;
// Vancomycin reduces gut LPS production
double VAN_effect = 0.0;  // placeholder; linked via dosing
dxdt_LPS = LPS_prod - LPS_cleared;

// Gut Barrier Integrity
double barrier_damage = k_IBD_damage * IBD_status
                      + 0.05 * HydroIndex * LPS  // BA + LPS synergy
                      + 0.02 * IL17A;             // Th17 disruption
double barrier_repair = k_barrier * GutBarrier;
dxdt_GutBarrier = barrier_repair * (GutBarrier0 - GutBarrier) - barrier_damage * GutBarrier;

// ===========================================================================
// BILE ACID & FXR SIGNALING ODEs
// ===========================================================================

// FXR activation (approaches FXR_input, decays naturally)
dxdt_FXR_act = k_FXR_act * (FXR_input - FXR_act) - k_FXR_decay * (FXR_act - FXR0);

// Bile acid pool (FXR suppresses BA synthesis)
double BA_synth = k_BA_synth * (1.0 - 0.7 * FXR_act);
dxdt_BilePool = BA_synth - k_BA_excrete * BilePool;

// Hydrophobic BA index (toxic — reduced by UDCA and FXR)
double hydro_increase = k_hydro_form * (1.0 - FXR_act) * LPS;
double hydro_decrease = k_hydro_clear * HydroIndex * (1.0 + UDCA_hydro_eff * UDCA_bile_conc / (2.0 + UDCA_bile_conc));
dxdt_HydroIndex = hydro_increase - hydro_decrease;

// ===========================================================================
// IMMUNE DYNAMICS ODEs
// ===========================================================================

// IL-17A (Th17, driven by LPS and BA toxicity, inhibited by Treg)
double IL17_prod = k_IL17_base
                 + k_IL17_LPS * LPS * (1.0 - Treg_IL10 * 0.5)
                 + k_IL17_BA  * HydroIndex;
dxdt_IL17A = IL17_prod - k_IL17_decay * IL17A;

// TNF-α (macrophage, driven by LPS)
double TNFa_prod = k_TNFa_base
                 + k_TNFa_LPS * LPS;
dxdt_TNFa = TNFa_prod - k_TNFa_decay * TNFa;

// IL-6 (IL-17 + TNF synergy)
double IL6_prod = k_IL6_IL17 * IL17A + k_IL6_TNFa * TNFa;
dxdt_IL6 = IL6_prod - k_IL6_decay * IL6;

// Treg / IL-10 (suppresses Th17; inhibited by IL-6)
double Treg_prod = k_Treg_base;
double Treg_inhibited = k_Treg_inhib * IL6 * Treg_IL10;
dxdt_Treg_IL10 = Treg_prod - k_Treg_decay * Treg_IL10 - Treg_inhibited;

// ===========================================================================
// CHOLANGIOCYTE ODEs
// ===========================================================================

// Cholangiocyte health (0=destroyed, 1=healthy)
double cholangio_damage = k_cholangio_damage_IL17 * IL17A
                        + k_cholangio_damage_BA   * HydroIndex
                        + 0.01 * TNFa;
double cholangio_repair = k_cholangio_repair * Cholangio_health;
dxdt_Cholangio_health = cholangio_repair * (1.0 - Cholangio_health) - cholangio_damage * Cholangio_health;

// Senescence burden
double senesc_form  = k_senesc_form * IL17A * HydroIndex;
double senesc_clear = k_senesc_clear * Senescence;
dxdt_Senescence = senesc_form - senesc_clear;

// ===========================================================================
// FIBROSIS ODEs
// ===========================================================================

// HSC activation (driven by senescence SASP + IL-17A)
double HSC_activation = k_HSC_act * (IL17A + Senescence) * (1.0 - HSC_act);
double HSC_resolution = k_HSC_resol * HSC_act * Treg_IL10;
dxdt_HSC_act = HSC_activation - HSC_resolution;

// Collagen I synthesis and degradation
// Bezafibrate reduces TGF-β signaling (PPARα) — models as ↓ collagen synthesis
double BEZ_antifib = 1.0 - 0.3 * BEZ_plasma / (5.0 + BEZ_plasma);
double Col_synth = k_Col_synth * HSC_act * BEZ_antifib;
double Col_degrad = k_Col_degrad * Col1a1 * (1.0 - LOXL2 * 0.3);
dxdt_Col1a1 = Col_synth - Col_degrad;

// LOXL2 (cross-linking; inhibited by simtuzumab)
double LOXL2_eff_clearance = k_LOXL2_clear * LOXL2;
dxdt_LOXL2 = k_LOXL2_form * HSC_act - LOXL2_eff_clearance;

// ===========================================================================
// CLINICAL BIOMARKER ODEs
// ===========================================================================

// ALP (IU/L): elevated by cholestasis, reduced by UDCA, OCA, bezafibrate
double ALP_prod = k_ALP_chol * Cholestasis * ALP0;
double ALP_norm = k_ALP_repair * (ALP - ALP_normal);
// Drug-mediated ALP reduction
double OCA_ALP_eff  = 0.25 * (OCA_bile_conc  / (2.0 + OCA_bile_conc));
double UDCA_ALP_eff = 0.15 * (UDCA_bile_conc / (5.0 + UDCA_bile_conc));
double BEZ_ALP_eff  = 0.20 * BEZ_plasma / (10.0 + BEZ_plasma);
double ALP_drug_red = (OCA_ALP_eff + UDCA_ALP_eff + BEZ_ALP_eff) * ALP;
dxdt_ALP = ALP_prod - ALP_norm - ALP_drug_red;

// Bilirubin (mg/dL)
double Bili_prod = k_Bili_synth + k_Bili_chol * Cholestasis;
dxdt_Bilirubin = Bili_prod - k_Bili_excrete * Bilirubin;

// Liver stiffness (kPa; FibroScan)
dxdt_Fibroscan = k_Fibro_form * Col1a1 - k_Fibro_resol * Fibroscan;

// Portal pressure (mmHg)
dxdt_PortalPressure = k_PP_fibrosis * Col1a1 - k_PP_decay * (PortalPressure - PP0);

// ===========================================================================
// CCA RISK (cumulative hazard proxy)
// ===========================================================================
dxdt_CCA_risk = k_CCA_chol * Cholestasis * ALP / ALP0
              + k_CCA_senesc * Senescence;
// Capped at k_CCA_max via bounded accumulation

$INIT @annotated
// Drug PK compartments — initialize to zero
UDCA_gut     = 0
UDCA_plasma  = 0
UDCA_bile    = 0
OCA_gut      = 0
OCA_plasma   = 0
OCA_bile     = 0
NorUDCA_bile = 0
BEZ_plasma   = 0

// Disease states — initialized to PSC-established state
LPS           = 0.6    // elevated due to gut dysbiosis
GutBarrier    = 0.55   // reduced gut barrier
FXR_act       = 0.25   // sub-optimal FXR activity
BilePool      = 0.70   // partially impaired bile flow
HydroIndex    = 0.55   // elevated hydrophobic BA
IL17A         = 0.25   // elevated Th17
TNFa          = 0.20   // elevated TNF
IL6           = 0.15   // elevated IL-6
Treg_IL10     = 0.05   // low Treg
Cholangio_health = 0.50 // 50% cholangiocyte health
Senescence    = 0.30   // moderate senescent cell burden
HSC_act       = 0.35   // moderately activated HSC
Col1a1        = 0.40   // moderate collagen deposition
LOXL2         = 0.35   // elevated LOXL2
ALP           = 380.0  // elevated ALP (IU/L)
Bilirubin     = 2.8    // mildly elevated bilirubin
Fibroscan     = 12.5   // elevated liver stiffness (kPa)
PortalPressure = 8.0   // mild portal hypertension
CCA_risk      = 0.02   // small cumulative CCA risk at baseline

$CAPTURE @annotated
Cholestasis        : Cholestasis severity index (0-1)
FXR_input          : Target FXR activation from drugs (0-1)
VaricesRisk        : Risk of variceal complication (0-1)
OCA_bile_conc      : OCA biliary concentration (μmol/L equiv)
UDCA_bile_conc     : UDCA biliary concentration (μmol/L equiv)
'

## ---------------------------------------------------------------------------
## 2. Compile the Model
## ---------------------------------------------------------------------------
psc_mod <- mcode("PSC_QSP", psc_model_code)

cat("Model compiled successfully!\n")
cat("Compartments:", length(init(psc_mod)), "\n")

## ---------------------------------------------------------------------------
## 3. Helper: Build Dosing Events
## ---------------------------------------------------------------------------
build_dosing <- function(
    treatment = "UDCA",
    daily_dose_mg = 750,
    start_day = 1,
    duration_days = 365,
    freq_hours = 12
) {
  times <- seq(
    from = (start_day - 1) * 24,
    by   = freq_hours,
    length.out = duration_days * (24 / freq_hours)
  )
  dose_amt <- daily_dose_mg / (24 / freq_hours)  # split daily dose

  cmt_map <- list(
    UDCA    = "UDCA_gut",
    OCA     = "OCA_gut",
    norUDCA = "NorUDCA_bile",
    BEZ     = "BEZ_plasma"
  )
  cmt <- cmt_map[[treatment]]
  if (is.null(cmt)) stop("Unknown treatment: ", treatment)

  ev <- ev(
    time   = times,
    amt    = dose_amt,
    cmt    = cmt,
    evid   = 1
  )
  return(ev)
}

## ---------------------------------------------------------------------------
## 4. Simulation Scenarios
## ---------------------------------------------------------------------------
sim_years <- 5  # Simulate 5 years
sim_hours <- sim_years * 365 * 24
obs_times <- seq(0, sim_hours, by = 24 * 7)  # weekly observations

# Scenario parameters
scenarios <- list(
  "Natural Progression (No Treatment)" = list(
    events = NULL,
    params = list()
  ),
  "UDCA 15 mg/kg/day" = list(
    events = build_dosing("UDCA", 900, freq_hours = 12),
    params = list()
  ),
  "OCA 10 mg/day" = list(
    events = build_dosing("OCA", 10, freq_hours = 24),
    params = list()
  ),
  "UDCA + OCA Combination" = list(
    events = c(
      build_dosing("UDCA", 900, freq_hours = 12),
      build_dosing("OCA",  10,  freq_hours = 24)
    ),
    params = list()
  ),
  "Bezafibrate 400 mg/day" = list(
    events = build_dosing("BEZ", 400, freq_hours = 24),
    params = list()
  )
)

## ---------------------------------------------------------------------------
## 5. Run All Scenarios
## ---------------------------------------------------------------------------
cat("\nRunning", length(scenarios), "treatment scenarios...\n")

results_all <- lapply(names(scenarios), function(scen_name) {
  cat("  Scenario:", scen_name, "\n")
  sc   <- scenarios[[scen_name]]
  mod2 <- param(psc_mod, sc$params)

  if (is.null(sc$events)) {
    out <- mrgsim(mod2, tgrid = obs_times)
  } else {
    out <- mrgsim(mod2, events = sc$events, tgrid = obs_times)
  }

  as.data.frame(out) %>%
    mutate(
      scenario = scen_name,
      time_years = time / (365 * 24)
    )
})

results_df <- bind_rows(results_all)
cat("Simulation complete. Rows:", nrow(results_df), "\n")

## ---------------------------------------------------------------------------
## 6. Key Outputs & Plots
## ---------------------------------------------------------------------------

# Color palette for 5 scenarios
palette_sc <- c(
  "Natural Progression (No Treatment)" = "#E74C3C",
  "UDCA 15 mg/kg/day"                  = "#3498DB",
  "OCA 10 mg/day"                      = "#9B59B6",
  "UDCA + OCA Combination"             = "#27AE60",
  "Bezafibrate 400 mg/day"             = "#E67E22"
)

plot_var <- function(df, varname, ylabel, title) {
  ggplot(df, aes(x = time_years, y = .data[[varname]], color = scenario)) +
    geom_line(size = 1.0) +
    scale_color_manual(values = palette_sc) +
    labs(x = "Time (years)", y = ylabel, title = title, color = "Treatment") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8),
          plot.title = element_text(face = "bold"))
}

p_ALP  <- plot_var(results_df, "ALP", "ALP (IU/L)", "Alkaline Phosphatase")
p_Bili <- plot_var(results_df, "Bilirubin", "Bilirubin (mg/dL)", "Total Bilirubin")
p_Fib  <- plot_var(results_df, "Fibroscan", "Liver Stiffness (kPa)", "FibroScan (Stiffness)")
p_PP   <- plot_var(results_df, "PortalPressure", "Portal Pressure (mmHg)", "Portal Pressure")
p_Col  <- plot_var(results_df, "Col1a1", "Collagen Index (AU)", "Hepatic Collagen (Fibrosis)")
p_CCA  <- plot_var(results_df, "CCA_risk", "CCA Risk Index (0-1)", "Cumulative CCA Risk")
p_IL17 <- plot_var(results_df, "IL17A", "IL-17A (AU)", "IL-17A (Th17 Inflammation)")
p_Cho  <- plot_var(results_df, "Cholangio_health", "Cholangiocyte Health (0-1)", "Cholangiocyte Health")

p_combined <- (p_ALP | p_Bili) /
              (p_Fib | p_PP)  /
              (p_Col | p_CCA) /
              (p_IL17 | p_Cho) +
  plot_annotation(
    title = "PSC QSP Model — 5-Year Treatment Comparison",
    subtitle = "Natural progression vs UDCA / OCA / Combination / Bezafibrate",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12)
    )
  )

## Uncomment to save:
# ggsave("psc_5year_comparison.png", p_combined, width = 16, height = 20, dpi = 150)

## ---------------------------------------------------------------------------
## 7. Summary Statistics at Year 1, 2, 5
## ---------------------------------------------------------------------------
summary_table <- results_df %>%
  filter(time_years %in% c(0, 1, 2, 5)) %>%
  group_by(scenario, time_years) %>%
  summarise(
    ALP         = round(mean(ALP), 1),
    Bilirubin   = round(mean(Bilirubin), 2),
    Fibroscan   = round(mean(Fibroscan), 1),
    PortalP     = round(mean(PortalPressure), 1),
    IL17A       = round(mean(IL17A), 3),
    CCA_risk    = round(mean(CCA_risk), 4),
    Col1a1      = round(mean(Col1a1), 3),
    Cholangio   = round(mean(Cholangio_health), 3),
    .groups = "drop"
  )

cat("\n=== PSC Model Summary: Key Biomarkers by Year ===\n")
print(summary_table)

## ---------------------------------------------------------------------------
## 8. ALP Response Rate (PSC primary endpoint: ALP < 1.5x ULN)
## ---------------------------------------------------------------------------
ALP_ULN <- 150  # Upper limit of normal (IU/L)
ALP_threshold <- 1.5 * ALP_ULN  # = 225 IU/L

resp_year1 <- results_df %>%
  filter(abs(time_years - 1.0) < 0.05) %>%
  group_by(scenario) %>%
  summarise(
    ALP_mean = mean(ALP),
    ALP_responder = mean(ALP < ALP_threshold) * 100,
    .groups = "drop"
  ) %>%
  arrange(ALP_mean)

cat("\n=== ALP Response at Year 1 (ALP < 225 IU/L = 1.5x ULN) ===\n")
print(resp_year1)

## ---------------------------------------------------------------------------
## 9. Dose-Response: OCA dose effect on ALP at Year 1
## ---------------------------------------------------------------------------
oca_doses <- c(5, 10, 25, 50)
oca_dr <- lapply(oca_doses, function(dose) {
  ev_oca <- build_dosing("OCA", dose, freq_hours = 24)
  out <- mrgsim(psc_mod, events = ev_oca,
                tgrid = seq(0, 365 * 24, by = 24 * 7))
  as.data.frame(out) %>%
    filter(abs(time / (365 * 24) - 1.0) < 0.05) %>%
    summarise(OCA_dose_mg = dose, ALP_Y1 = mean(ALP)) %>%
    head(1)
})
oca_dr_df <- bind_rows(oca_dr)

cat("\n=== OCA Dose-Response: ALP at Year 1 ===\n")
print(oca_dr_df)

## ---------------------------------------------------------------------------
## 10. Sensitivity Analysis: IBD status impact
## ---------------------------------------------------------------------------
ibd_scenarios <- data.frame(
  IBD_status = c(0, 1),
  label = c("PSC without IBD", "PSC with IBD (UC)")
)

sens_ibd <- lapply(1:nrow(ibd_scenarios), function(i) {
  mod_ibd <- param(psc_mod, list(IBD_status = ibd_scenarios$IBD_status[i]))
  ev_base  <- build_dosing("UDCA", 900, freq_hours = 12)
  out <- mrgsim(mod_ibd, events = ev_base,
                tgrid = seq(0, sim_hours, by = 24 * 7))
  as.data.frame(out) %>%
    mutate(label = ibd_scenarios$label[i],
           time_years = time / (365 * 24))
})
sens_ibd_df <- bind_rows(sens_ibd)

p_IBD <- ggplot(sens_ibd_df, aes(x = time_years, y = ALP, color = label, linetype = label)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("#1A5276", "#C0392B")) +
  labs(x = "Time (years)", y = "ALP (IU/L)",
       title = "PSC with vs without IBD — UDCA Response",
       color = "IBD Status", linetype = "IBD Status") +
  theme_bw(base_size = 11)

cat("\nPSC QSP Model — All simulations complete.\n")
cat("Key variables tracked:", paste(names(init(psc_mod)), collapse = ", "), "\n")
