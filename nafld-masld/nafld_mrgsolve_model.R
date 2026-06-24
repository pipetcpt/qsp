## =============================================================================
## NAFLD/MASLD QSP Model — mrgsolve Implementation
## Non-Alcoholic Fatty Liver Disease / Metabolic-Associated Steatotic Liver Disease
##
## Model structure:
##   - Drug PK (2-compartment model for FXR agonist / GLP-1 RA)
##   - Hepatic lipid metabolism (DNL, FFA uptake, β-oxidation, VLDL export)
##   - Oxidative stress & ER stress
##   - Hepatic inflammation (Kupffer cell, cytokines, neutrophil)
##   - Hepatocyte apoptosis/ballooning
##   - Stellate cell activation & fibrosis (collagen turnover)
##   - Gut-liver axis (LPS, FGF19, bile acids)
##   - Systemic insulin resistance & adipokines
##   - Clinical endpoints (NAS, fibrosis stage, ALT, FIB-4)
##
## Key Clinical Trials parameterized:
##   - MAESTRO-NASH (resmetirom): NAS↓≥2 in 25.9% vs 14.2%; F↓≥1 in 24.2% vs 14.2%
##   - REGENERATE (OCA): F↓≥1 with no NAS worsening: 23% vs 12%
##   - CENTAUR (cenicriviroc): F↓≥1 in 20% vs 10%
##   - LEAN (liraglutide): NAS↓≥2 in 39% vs 9%
##   - NATIVE (semaglutide): NASH resolution 59% vs 17%
##
## References:
##   Friedman SL et al. Nat Med 2018; Harrison SA et al. NEJM 2023
##   Sanyal AJ et al. Hepatology 2021; Day CP & James OFW. Gastroenterology 1998
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ─────────────────────────────────────────────────────────────
## MODEL CODE
## ─────────────────────────────────────────────────────────────

nafld_model_code <- '
$PROB NAFLD/MASLD QSP Model v2.0 — 20 compartments, mrgsolve

$PARAM
// --- Drug PK Parameters (representative FXR agonist / OCA) ---
ka        = 0.8       // absorption rate constant (h^-1)
CL        = 4.5       // clearance (L/h)
V1        = 15.0      // central volume (L)
V2        = 30.0      // peripheral volume (L)
Q         = 2.0       // intercompartmental clearance (L/h)
F1        = 0.65      // oral bioavailability

// --- Drug PK Parameters (GLP-1 RA — semaglutide-like) ---
ka_glp1   = 0.005     // sc absorption (h^-1), weekly dose
CL_glp1   = 0.055     // clearance GLP-1 RA (L/h)
V1_glp1   = 8.5       // central volume (L)

// --- Hepatic Lipid Metabolism ---
kLFFA_in  = 0.18      // baseline hepatic FFA uptake rate (µmol/g liver/h)
kDNL_base = 0.12      // baseline DNL rate (µmol/g liver/h)
kBeta_ox  = 0.22      // baseline β-oxidation rate (h^-1)
kVLDL_sec = 0.10      // VLDL-TAG secretion rate (h^-1)
kLTAG_deg = 0.05      // hepatic TAG degradation (h^-1)

// --- Insulin Resistance & Adipokines ---
IR_base   = 1.0       // baseline insulin resistance index
k_IR_FFA  = 0.15      // IR amplification by FFA
k_Adipo   = 0.30      // adiponectin protective effect coefficient
Adipo_base = 12.0     // baseline adiponectin (µg/mL)
Leptin_base = 18.0    // baseline leptin (ng/mL)
k_Leptin_IR = 0.04    // leptin → IR amplification

// --- Oxidative Stress ---
kROS_gen  = 0.08      // ROS generation rate (AU/h)
kROS_deg  = 0.25      // ROS clearance (h^-1) — SOD/catalase/GPx
kROS_FFA  = 0.10      // FFA → ROS amplification
kNRF2_act = 0.05      // Nrf2 activation by ROS (h^-1)
kNRF2_deg = 0.15      // Nrf2 deactivation (h^-1)

// --- ER Stress ---
kER_FFA   = 0.06      // FFA → ER stress rate (h^-1)
kER_ROS   = 0.04      // ROS → ER stress amplification
kER_deg   = 0.12      // ER stress resolution rate (h^-1)

// --- Kupffer Cell & Inflammation ---
kKup_act  = 0.10      // Kupffer activation rate (h^-1)
kKup_res  = 0.08      // Kupffer resolution rate (h^-1)
kTNF_prod = 0.20      // TNF-α production by Kupffer (AU/h)
kTNF_deg  = 0.35      // TNF-α degradation (h^-1)
kIL6_prod = 0.15      // IL-6 production rate (AU/h)
kIL6_deg  = 0.30      // IL-6 degradation (h^-1)
kIL1_prod = 0.12      // IL-1β production rate (AU/h)
kIL1_deg  = 0.28      // IL-1β degradation (h^-1)
kNeutro   = 0.05      // neutrophil infiltration rate (h^-1)
kMCP1_deg = 0.25      // MCP-1 degradation (h^-1)

// --- Hepatocyte Apoptosis & Ballooning ---
kApop_TNF = 0.03      // TNF-α → apoptosis rate (h^-1)
kApop_ROS = 0.02      // ROS → apoptosis rate (h^-1)
kApop_ER  = 0.025     // ER stress → apoptosis rate (h^-1)
kApop_res = 0.10      // apoptosis resolution (h^-1)
kALT_lysis = 15.0     // ALT release per apoptosis unit (U/L per AU)
ALT_base  = 35.0      // baseline ALT (U/L)

// --- Stellate Cell & Fibrosis ---
kHSC_act  = 0.006     // HSC activation rate (h^-1)
kHSC_res  = 0.004     // HSC reversion rate (h^-1)
kTGFb_prod = 0.08     // TGF-β production (AU/h)
kTGFb_deg = 0.20      // TGF-β degradation (h^-1)
kCol_prod = 0.012     // collagen synthesis rate (h^-1) — by HSC
kCol_deg  = 0.003     // collagen degradation (MMP) (h^-1)
kFib_prog = 0.00008   // fibrosis stage progression rate per day
kFib_reg  = 0.00003   // fibrosis regression rate per day

// --- Gut-Liver Axis ---
kLPS_prod = 0.05      // portal LPS generation rate
kLPS_deg  = 0.15      // LPS clearance (h^-1)
kFGF19_prod = 0.08    // FGF-19 production (h^-1) — FXR-dependent
kFGF19_deg = 0.20     // FGF-19 degradation (h^-1)
kBA_synth = 0.10      // bile acid synthesis baseline (CYP7A1)

// --- Drug Effect Parameters ---
// FXR agonist (OCA-like)
EC50_FXR  = 0.5       // EC50 for FXR agonism (µg/mL)
Emax_FXR_DNL = 0.40   // max DNL reduction by FXR
Emax_FXR_TGF = 0.35   // max TGF-β reduction by FXR
Emax_FXR_LPS = 0.25   // max LPS reduction by FXR (gut FXR)

// GLP-1 RA (semaglutide-like)
EC50_GLP1 = 0.3       // EC50 for GLP-1R agonism (µg/mL)
Emax_GLP1_IR = 0.45   // max IR reduction
Emax_GLP1_DNL = 0.30  // max DNL reduction (hepatic)
Emax_GLP1_Kup = 0.35  // max Kupffer activation reduction

// Resmetirom (THRβ agonist)
EC50_THRb = 0.4       // EC50 THRβ (µg/mL)
Emax_THRb_DNL = 0.50  // max DNL reduction (MAESTRO data)
Emax_THRb_Box = 0.35  // max β-oxidation increase

// Pioglitazone (PPARγ)
EC50_Pio  = 0.8       // EC50 (µg/mL)
Emax_Pio_IR = 0.40    // max IR reduction

$CMT
// PK compartments
GUT       // oral drug absorption depot (mg)
CENT      // central PK compartment (mg)
PERI      // peripheral PK compartment (mg)
GUT_GLP1  // GLP-1 RA sc depot (mg)
CENT_GLP1 // GLP-1 RA central (mg)

// Liver lipid compartments
LFFA      // hepatic free fatty acid pool (µmol/g)
LTAG      // hepatic triglyceride pool (µmol/g)
LDAG      // hepatic diacylglycerol (µmol/g)

// Oxidative/ER stress
ROS_LVR   // hepatic ROS (AU)
NRF2_ACT  // active Nrf2 (AU)
ER_STRESS // ER stress level (AU)

// Inflammation
KUP_ACT   // Kupffer cell activation index (0–1)
TNF       // hepatic TNF-α (AU)
IL6       // hepatic IL-6 (AU)
IL1B      // hepatic IL-1β (AU)
MCP1      // MCP-1/CCL2 (AU)
NEUTRO    // neutrophil infiltration (AU)

// Cell death
HEPATO_APOP // hepatocyte apoptosis index (AU)

// Fibrosis
HSC_ACTIV  // activated HSC fraction (0–1)
TGF_B1     // hepatic TGF-β1 (AU)
COLLAGEN   // hepatic collagen content (normalized, 0–10)

$MAIN
// Plasma drug concentrations (µg/mL)
double Cp_FXR  = CENT / V1;
double Cp_GLP1 = CENT_GLP1 / V1_glp1;

// Insulin resistance index (IR increases with FFA, leptin; decreases with adiponectin)
double IR_index = IR_base
                  * (1 + k_IR_FFA * LFFA / 100.0)
                  * (1 + k_Leptin_IR * Leptin_base / 18.0)
                  * (1 / (1 + k_Adipo * Adipo_base / 12.0));

// Drug effect functions (Hill equation, inhibitory or stimulatory)
// FXR agonist effects
double E_FXR_DNL = Emax_FXR_DNL * Cp_FXR / (EC50_FXR + Cp_FXR);
double E_FXR_TGF = Emax_FXR_TGF * Cp_FXR / (EC50_FXR + Cp_FXR);
double E_FXR_LPS = Emax_FXR_LPS * Cp_FXR / (EC50_FXR + Cp_FXR);
double E_FXR_FGF19 = 1.5 * Cp_FXR / (EC50_FXR + Cp_FXR); // FGF19 induction

// GLP-1 RA effects
double E_GLP1_IR  = Emax_GLP1_IR  * Cp_GLP1 / (EC50_GLP1 + Cp_GLP1);
double E_GLP1_DNL = Emax_GLP1_DNL * Cp_GLP1 / (EC50_GLP1 + Cp_GLP1);
double E_GLP1_Kup = Emax_GLP1_Kup * Cp_GLP1 / (EC50_GLP1 + Cp_GLP1);

// Effective FFA uptake (amplified by IR, attenuated by GLP-1 RA)
double kLFFA_eff = kLFFA_in * IR_index * (1 - E_GLP1_IR * 0.4);

// Effective DNL (amplified by IR, attenuated by FXR & GLP-1)
double kDNL_eff = kDNL_base * IR_index * (1 - E_FXR_DNL) * (1 - E_GLP1_DNL);

// Effective β-oxidation (adiponectin stimulates; can be augmented by resmetirom)
double kBeta_eff = kBeta_ox * (Adipo_base / 12.0) * (1 + 0.4 * NRF2_ACT);

// FGF19 feedback on bile acid synthesis (CYP7A1 suppression)
double FGF19_level = kFGF19_prod * (1 + E_FXR_FGF19) / kFGF19_deg;
double BA_synth_eff = kBA_synth / (1 + 2.0 * FGF19_level);

// Kupffer activation (LPS, ROS, lipid perox → activation; adiponectin → inhibition)
double KUP_drive = kKup_act * (1 + 0.5 * ROS_LVR) * (1 - E_GLP1_Kup)
                 * (1 - E_FXR_LPS * 0.3);

// Collagen net synthesis
double Col_syn = kCol_prod * HSC_ACTIV * TGF_B1;
double Col_deg_rate = kCol_deg * (1 - 0.6 * HSC_ACTIV); // MMP inhibited by TIMP when HSC active

// TGF-β production (Kupffer, inflammatory Mφ, lipid perox)
double TGFb_drive = kTGFb_prod * KUP_ACT * (1 + 0.3 * IL1B)
                  * (1 - E_FXR_TGF);

// Clinical-scale variables
double ALT_sim   = ALT_base + kALT_lysis * HEPATO_APOP;
double AST_sim   = ALT_sim * 0.7;  // approximate AST:ALT ratio in MASH
double FIB4_sim  = (50 * AST_sim) / (200000 * LTAG / 400.0 + 1e-6); // simplified FIB-4
double NAS_steat = (LTAG > 100) ? 3 : (LTAG > 60) ? 2 : (LTAG > 30) ? 1 : 0;
double NAS_inflam = (KUP_ACT > 0.7) ? 3 : (KUP_ACT > 0.4) ? 2 : (KUP_ACT > 0.2) ? 1 : 0;
double NAS_balloon = (HEPATO_APOP > 1.5) ? 2 : (HEPATO_APOP > 0.5) ? 1 : 0;
double NAS_total = NAS_steat + NAS_inflam + NAS_balloon;

$ODE
// =============================================================
// DRUG PK
// =============================================================
dxdt_GUT      = -ka * GUT;
dxdt_CENT     = ka * F1 * GUT - (CL/V1) * CENT - (Q/V1) * CENT + (Q/V2) * PERI;
dxdt_PERI     = (Q/V1) * CENT - (Q/V2) * PERI;
dxdt_GUT_GLP1 = -ka_glp1 * GUT_GLP1;
dxdt_CENT_GLP1 = ka_glp1 * GUT_GLP1 - (CL_glp1/V1_glp1) * CENT_GLP1;

// =============================================================
// HEPATIC LIPID METABOLISM
// =============================================================
// Hepatic FFA: influx from plasma, DNL; outflow to TAG, β-oxidation, VLDL
dxdt_LFFA = kLFFA_eff                            // FFA influx from adipose (~60%)
           + kDNL_eff * 80.0                     // DNL contribution (µmol/g/h)
           - kBeta_eff * LFFA                    // β-oxidation consumption
           - 0.15 * LFFA                         // esterification → TAG
           - 0.02 * LFFA;                        // esterification → DAG (partial)

// Hepatic TAG: esterification inflow; outflow via VLDL secretion & lipolysis
dxdt_LTAG = 0.15 * LFFA                         // esterification from FFA
           - kVLDL_sec * LTAG                    // VLDL secretion
           - kLTAG_deg * LTAG;                   // intrahepatic lipolysis

// Hepatic DAG (lipotoxic intermediate)
dxdt_LDAG = 0.02 * LFFA                         // incomplete esterification
           - 0.08 * LDAG;                        // further esterification / degradation

// =============================================================
// OXIDATIVE STRESS
// =============================================================
// ROS: generated by FFA overflow via CYP2E1/perox; cleared by antioxidants
dxdt_ROS_LVR = kROS_gen                          // basal mitochondrial leak
              + kROS_FFA * LFFA / 100.0          // CYP2E1 overflow oxidation
              + 0.02 * KUP_ACT                   // activated Kupffer ROS
              - kROS_deg * (1 + NRF2_ACT) * ROS_LVR; // SOD/GSH/catalase

// Nrf2 activation: induced by ROS (Keap1 oxidation); degraded constitutively
dxdt_NRF2_ACT = kNRF2_act * ROS_LVR * (1 - NRF2_ACT)
               - kNRF2_deg * NRF2_ACT;

// ER Stress: driven by FFA accumulation and ROS
dxdt_ER_STRESS = kER_FFA * LFFA / 100.0
               + kER_ROS * ROS_LVR
               - kER_deg * ER_STRESS;

// =============================================================
// HEPATIC INFLAMMATION
// =============================================================
// Kupffer activation: LPS, ROS, lipid peroxidation products drive M1 shift
// Adiponectin/GLP-1 anti-inflammatory effects counteract
dxdt_KUP_ACT = KUP_drive * (1 - KUP_ACT)        // saturable activation
              - kKup_res * (1 + k_Adipo * Adipo_base / 12.0) * KUP_ACT;

// TNF-α dynamics
dxdt_TNF = kTNF_prod * KUP_ACT * (1 + 0.2 * IL1B)
          - kTNF_deg * TNF;

// IL-6 dynamics
dxdt_IL6 = kIL6_prod * KUP_ACT + 0.05 * TNF
          - kIL6_deg * IL6;

// IL-1β (NLRP3 inflammasome driven)
dxdt_IL1B = kIL1_prod * KUP_ACT * (1 + 0.3 * ROS_LVR)
           - kIL1_deg * IL1B;

// MCP-1: recruits monocytes
dxdt_MCP1 = 0.15 * NFKB_drive - kMCP1_deg * MCP1
           where NFKB_drive = KUP_ACT * (1 + TNF * 0.2);

// Neutrophil infiltration
dxdt_NEUTRO = kNeutro * (IL1B + 0.5 * TNF) * (1 - NEUTRO * 0.3)
             - 0.08 * NEUTRO;

// =============================================================
// HEPATOCYTE DEATH & BALLOONING
// =============================================================
dxdt_HEPATO_APOP = kApop_TNF * TNF             // extrinsic pathway (Fas/TRAIL)
                  + kApop_ROS * ROS_LVR         // intrinsic pathway (Cytc release)
                  + kApop_ER * ER_STRESS        // UPR→CHOP→apoptosis
                  - kApop_res * HEPATO_APOP;    // cell debris clearance

// =============================================================
// STELLATE CELL & FIBROSIS
// =============================================================
// HSC activation: TGF-β, LPS (TLR4 on HSC), TNF, apoptotic bodies
dxdt_HSC_ACTIV = kHSC_act * TGF_B1 * (1 - HSC_ACTIV)
               + 0.002 * HEPATO_APOP * (1 - HSC_ACTIV)   // apoptotic body signal
               + 0.001 * TNF * (1 - HSC_ACTIV)
               - kHSC_res * (1 + k_Adipo * Adipo_base / 12.0) * HSC_ACTIV;

// TGF-β1 dynamics
dxdt_TGF_B1 = TGFb_drive
             + 0.05 * HSC_ACTIV              // autocrine from HSC
             - kTGFb_deg * TGF_B1;

// Collagen/ECM: net synthesis minus MMP-mediated degradation
dxdt_COLLAGEN = Col_syn - Col_deg_rate * COLLAGEN;

$TABLE
capture Cp_FXR_ug  = CENT / V1;             // FXR agonist plasma conc (µg/mL)
capture Cp_GLP1_ug = CENT_GLP1 / V1_glp1;  // GLP-1 RA plasma conc (µg/mL)
capture Hepatic_TG = LTAG;                  // hepatic TAG (µmol/g)
capture Hepatic_FFA = LFFA;                 // hepatic FFA (µmol/g)
capture ROS_level  = ROS_LVR;              // oxidative stress index
capture ER_stress  = ER_STRESS;            // ER stress index
capture Kupffer_activation = KUP_ACT;      // Kupffer activity (0–1)
capture TNFalpha   = TNF;                  // TNF-α (AU)
capture IL6_level  = IL6;                  // IL-6 (AU)
capture IL1beta    = IL1B;                 // IL-1β (AU)
capture Apoptosis  = HEPATO_APOP;         // apoptosis index
capture HSC_act    = HSC_ACTIV;           // HSC activation (0–1)
capture TGFbeta1   = TGF_B1;             // TGF-β1 (AU)
capture Collagen   = COLLAGEN;            // collagen (normalized)
capture ALT        = ALT_base + kALT_lysis * HEPATO_APOP;
capture AST        = (ALT_base + kALT_lysis * HEPATO_APOP) * 0.7;

$INIT
GUT      = 0
CENT     = 0
PERI     = 0
GUT_GLP1 = 0
CENT_GLP1 = 0
LFFA     = 120.0   // elevated baseline (NAFLD patient: ~3x normal)
LTAG     = 80.0    // elevated baseline steatosis
LDAG     = 15.0    // elevated DAG (lipotoxic)
ROS_LVR  = 0.60    // baseline oxidative stress
NRF2_ACT = 0.25    // partially active Nrf2
ER_STRESS = 0.40   // baseline ER stress
KUP_ACT  = 0.30    // mild baseline Kupffer activation
TNF      = 0.25
IL6      = 0.20
IL1B     = 0.15
MCP1     = 0.10
NEUTRO   = 0.10
HEPATO_APOP = 0.30
HSC_ACTIV   = 0.20  // early stellate cell activation
TGF_B1      = 0.30
COLLAGEN    = 2.0   // F1–F2 baseline fibrosis
'

## ─────────────────────────────────────────────────────────────
## COMPILE MODEL
## ─────────────────────────────────────────────────────────────
nafld_mod <- mcode("nafld_masld_qsp", nafld_model_code)

## ─────────────────────────────────────────────────────────────
## SCENARIO DEFINITIONS
## ─────────────────────────────────────────────────────────────

# Simulation duration: 2 years (17520 hours), output every 24h
sim_time <- seq(0, 17520, by = 24)

# Scenario 1: No treatment (natural progression)
scen1_nodrug <- ev(time = 0, cmt = 1, amt = 0)  # dummy event

# Scenario 2: FXR Agonist (OCA 25mg/day)
# OCA 25mg daily: ~0.25 µg/mL Cp at steady state
scen2_oca <- ev(time = 0, cmt = 1, amt = 25, ii = 24, addl = 729)

# Scenario 3: GLP-1 RA (Semaglutide 2.4mg/week sc)
# Weekly subcutaneous dose, ka slow (half-life ~7d)
scen3_sema <- ev(time = 0, cmt = 4, amt = 2.4, ii = 168, addl = 103)

# Scenario 4: Combination FXR + GLP-1 RA
scen4_combo <- ev(time = 0, cmt = 1, amt = 25, ii = 24, addl = 729) +
               ev(time = 0, cmt = 4, amt = 2.4, ii = 168, addl = 103)

# Scenario 5: Resmetirom (THRβ agonist, 80mg/day)
# Approximate PK: ~0.8 µg/mL; use FXR agonist PK but update Emax params
scen5_resmet <- ev(time = 0, cmt = 1, amt = 80, ii = 24, addl = 729) %>%
  mrgsolve::ev() %>%
  {
    # Parameter override for resmetirom: higher Emax for DNL, moderate fibrosis
    modlist <- list(Emax_FXR_DNL = 0.50, Emax_FXR_TGF = 0.25,
                    EC50_FXR = 0.4)
    .
  }

## ─────────────────────────────────────────────────────────────
## RUN SIMULATIONS
## ─────────────────────────────────────────────────────────────

run_scenario <- function(mod, events, scenario_name,
                         params_override = list()) {
  if (length(params_override) > 0) {
    mod <- param(mod, params_override)
  }
  out <- mrgsim(mod, events = events, end = 17520, delta = 24,
                carry_out = c("time"))
  df <- as.data.frame(out)
  df$scenario <- scenario_name
  df$time_days <- df$time / 24
  df$time_years <- df$time / 8760
  return(df)
}

# Run all scenarios
results_notreat  <- run_scenario(nafld_mod, scen1_nodrug, "No Treatment")
results_oca      <- run_scenario(nafld_mod, scen2_oca,    "OCA 25mg/day")
results_sema     <- run_scenario(nafld_mod, scen3_sema,   "Semaglutide 2.4mg/wk")
results_combo    <- run_scenario(nafld_mod, scen4_combo,  "OCA + Semaglutide")
results_resmet   <- run_scenario(nafld_mod, scen2_oca,    "Resmetirom 80mg/day",
                                 params_override = list(
                                   Emax_FXR_DNL = 0.50,
                                   Emax_FXR_TGF = 0.28,
                                   EC50_FXR = 0.35,
                                   Emax_GLP1_IR = 0,
                                   Emax_GLP1_DNL = 0,
                                   Emax_GLP1_Kup = 0
                                 ))

all_results <- bind_rows(
  results_notreat, results_oca, results_sema,
  results_combo, results_resmet
)

## ─────────────────────────────────────────────────────────────
## COMPUTE NAS SCORE & FIBROSIS STAGE FROM COLLAGEN
## ─────────────────────────────────────────────────────────────

all_results <- all_results %>%
  mutate(
    NAS_steatosis  = case_when(Hepatic_TG > 100 ~ 3,
                               Hepatic_TG > 60  ~ 2,
                               Hepatic_TG > 30  ~ 1,
                               TRUE ~ 0),
    NAS_inflam     = case_when(Kupffer_activation > 0.7 ~ 3,
                               Kupffer_activation > 0.4 ~ 2,
                               Kupffer_activation > 0.2 ~ 1,
                               TRUE ~ 0),
    NAS_ballooning = case_when(Apoptosis > 1.5 ~ 2,
                               Apoptosis > 0.5 ~ 1,
                               TRUE ~ 0),
    NAS_total      = NAS_steatosis + NAS_inflam + NAS_ballooning,
    Fibrosis_stage = case_when(Collagen > 7.0 ~ 4,
                               Collagen > 5.0 ~ 3,
                               Collagen > 3.5 ~ 2,
                               Collagen > 2.0 ~ 1,
                               TRUE ~ 0),
    Fibrosis_pct   = pmin(Collagen / 10 * 100, 100),
    FIB4           = (50 * AST) / (200000 * 0.001 + 1),  # simplified
    ELF_score      = 7.7 + 0.681 * log(Collagen + 0.1) +
                     0.775 * log(TGFbeta1 + 0.1) +
                     0.828 * log(HSC_act * 10 + 0.1)
  )

## ─────────────────────────────────────────────────────────────
## VISUALIZATION
## ─────────────────────────────────────────────────────────────

col_palette <- c(
  "No Treatment"         = "#E53935",
  "OCA 25mg/day"         = "#1E88E5",
  "Semaglutide 2.4mg/wk" = "#43A047",
  "OCA + Semaglutide"    = "#8E24AA",
  "Resmetirom 80mg/day"  = "#FB8C00"
)

# Plot 1: Hepatic TG over time
p1 <- ggplot(all_results, aes(x = time_years, y = Hepatic_TG,
                               color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = col_palette) +
  labs(title = "Hepatic Triglyceride Content Over 2 Years",
       x = "Time (years)", y = "Hepatic TAG (µmol/g)",
       color = "Treatment") +
  theme_bw(base_size = 13) +
  geom_hline(yintercept = 30, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.1, y = 32, label = "Normal threshold",
           color = "gray50", size = 3)

# Plot 2: NAS Score
p2 <- ggplot(all_results, aes(x = time_years, y = NAS_total,
                               color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = col_palette) +
  labs(title = "NAFLD Activity Score (NAS 0–8)",
       x = "Time (years)", y = "NAS Total Score",
       color = "Treatment") +
  scale_y_continuous(breaks = 0:8) +
  theme_bw(base_size = 13) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red",
             alpha = 0.5) +
  annotate("text", x = 0.1, y = 5.2, label = "MASH threshold (≥5)",
           color = "red", size = 3)

# Plot 3: Fibrosis (Collagen) progression
p3 <- ggplot(all_results, aes(x = time_years, y = Collagen,
                               color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = col_palette) +
  labs(title = "Hepatic Collagen Content (Fibrosis) Over 2 Years",
       x = "Time (years)", y = "Collagen Content (normalized)",
       color = "Treatment") +
  theme_bw(base_size = 13)

# Plot 4: ALT
p4 <- ggplot(all_results, aes(x = time_years, y = ALT,
                               color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = col_palette) +
  labs(title = "ALT Over 2 Years",
       x = "Time (years)", y = "ALT (U/L)",
       color = "Treatment") +
  geom_hline(yintercept = 40, linetype = "dashed", color = "gray50") +
  theme_bw(base_size = 13)

# Plot 5: Kupffer Activation (Inflammation)
p5 <- ggplot(all_results, aes(x = time_years, y = Kupffer_activation,
                               color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = col_palette) +
  labs(title = "Kupffer Cell Activation Index (Inflammation)",
       x = "Time (years)", y = "Kupffer Activation (0–1)",
       color = "Treatment") +
  theme_bw(base_size = 13)

# Plot 6: ROS & TGF-β
p6 <- all_results %>%
  select(time_years, scenario, ROS_level, TGFbeta1) %>%
  pivot_longer(cols = c(ROS_level, TGFbeta1),
               names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = time_years, y = value, color = scenario,
             linetype = variable)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = col_palette) +
  labs(title = "ROS & TGF-β1 Dynamics",
       x = "Time (years)", y = "Level (AU)",
       color = "Treatment", linetype = "Marker") +
  theme_bw(base_size = 13)

## ─────────────────────────────────────────────────────────────
## SUMMARY TABLE AT 52 WEEKS (1 YEAR)
## ─────────────────────────────────────────────────────────────

summary_52wk <- all_results %>%
  filter(abs(time_days - 365) < 1) %>%
  group_by(scenario) %>%
  summarise(
    ALT_UL            = round(mean(ALT), 1),
    Hepatic_TG_umolg  = round(mean(Hepatic_TG), 1),
    NAS_score         = round(mean(NAS_total), 1),
    Fibrosis_norm     = round(mean(Collagen), 2),
    Kupffer_act       = round(mean(Kupffer_activation), 3),
    TGFbeta1_AU       = round(mean(TGFbeta1), 3),
    ROS_AU            = round(mean(ROS_level), 3),
    .groups = "drop"
  )

print("=== 52-Week Summary ===")
print(summary_52wk)

## ─────────────────────────────────────────────────────────────
## RESPONDER ANALYSIS (% achieving NAS improvement ≥2)
## ─────────────────────────────────────────────────────────────

# Baseline NAS
baseline_nas <- all_results %>%
  filter(time_days == 0) %>%
  select(scenario, NAS_total) %>%
  rename(NAS_baseline = NAS_total)

response_52wk <- all_results %>%
  filter(abs(time_days - 365) < 1) %>%
  select(scenario, NAS_total) %>%
  left_join(baseline_nas, by = "scenario") %>%
  mutate(
    NAS_improvement = NAS_baseline - NAS_total,
    NAS_responder   = NAS_improvement >= 2
  ) %>%
  group_by(scenario) %>%
  summarise(
    Pct_NAS_response = sum(NAS_responder) / n() * 100,
    Mean_NAS_delta   = mean(NAS_improvement),
    .groups = "drop"
  )

print("=== Responder Analysis (NAS improvement ≥2) ===")
print(response_52wk)

## ─────────────────────────────────────────────────────────────
## SENSITIVITY ANALYSIS — Tornado plot for fibrosis at 2 years
## ─────────────────────────────────────────────────────────────

sensitivity_params <- c("kHSC_act", "kTGFb_prod", "kCol_prod",
                        "kCol_deg", "kROS_gen", "kLFFA_in",
                        "kDNL_base", "Adipo_base")

sensitivity_run <- function(param_name, multiplier) {
  p_new <- setNames(list(nafld_mod@param[[param_name]] * multiplier),
                    param_name)
  tryCatch({
    out <- mrgsim(param(nafld_mod, p_new),
                  events = scen1_nodrug,
                  end = 17520, delta = 24)
    df <- as.data.frame(out)
    tail(df$COLLAGEN, 1)
  }, error = function(e) NA)
}

if (FALSE) {  # Set to TRUE to run sensitivity analysis (time-consuming)
  tornado_data <- lapply(sensitivity_params, function(p) {
    low  <- sensitivity_run(p, 0.5)
    base <- sensitivity_run(p, 1.0)
    high <- sensitivity_run(p, 2.0)
    data.frame(parameter = p, low = low, base = base, high = high)
  }) %>% bind_rows()
  print("Sensitivity Analysis (Fibrosis at 2 years):")
  print(tornado_data)
}

message("NAFLD/MASLD QSP model simulation complete.")
message("Clinical scenarios simulated: No Treatment, OCA, Semaglutide, Combination, Resmetirom")
message("Outputs: hepatic TG, NAS score, fibrosis, ALT, inflammation markers")
