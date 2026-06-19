################################################################################
# Chronic Recurrent Urolithiasis — QSP Model (mrgsolve)
#
# Model structure:
#   - Drug PK: Thiazide (2-cpt), Allopurinol/Oxypurinol (2-cpt),
#              Potassium Citrate (1-cpt), Tamsulosin (1-cpt)
#   - Calcium dynamics: gut absorption, plasma pool, renal excretion
#   - Oxalate dynamics: hepatic synthesis, gut absorption, plasma, urine
#   - Uric acid dynamics: production (XO), plasma, renal excretion
#   - Citrate dynamics: plasma and urine
#   - Supersaturation (CaOx, CaP, UA) surrogate models
#   - Stone growth kinetics (CaOx dominant)
#   - Renal inflammatory injury index
#   - GFR trajectory
#
# Calibrated to published clinical trial data:
#   - Thiazide for idiopathic hypercalciuria (Laerum & Larsen, 1984)
#   - Allopurinol for UA/CaOx stones (Ettinger et al. NEJM 1986)
#   - Potassium citrate for CaOx/UA stones (Pak et al. 1985; NEJM 2014)
#   - Lumasiran PH1 (ILLUMINATE-A, NEJM 2021)
#   - ESWL recurrence data (Pearle et al., AUA Guidelines 2014)
#
# References: see uri_references.md
################################################################################

library(mrgsolve)
library(tidyverse)
library(ggplot2)

# ─────────────────────────────────────────────────────────────────────────────
# MODEL CODE
# ─────────────────────────────────────────────────────────────────────────────
uri_model_code <- '
$PROB
Chronic Recurrent Urolithiasis QSP Model
15 ODE compartments + algebraic outputs

$PARAM @annotated
// Patient demographics
WT      : 80   : Body weight (kg)
SEX     : 1    : Sex (1=male 0=female)
AGE     : 45   : Age (years)

// Dietary inputs (baseline daily amounts)
Ca_diet  : 1000  : Dietary calcium (mg/day)
Ox_diet  : 200   : Dietary oxalate (mg/day)
Pro_diet : 96    : Animal protein intake (g/day)
Na_diet  : 150   : Dietary sodium (mEq/day)
FluIn    : 2.0   : Fluid intake (L/day)
VitC     : 500   : Vitamin C intake (mg/day)

// Comorbidity flags
HPT     : 0    : Primary hyperparathyroidism (1=yes)
MetSyn  : 0    : Metabolic syndrome (1=yes)
IBD     : 0    : IBD/malabsorption (1=yes)
PH1     : 0    : Primary hyperoxaluria type 1 (1=yes)
Cystin  : 0    : Cystinuria (1=yes)

// Baseline physiological values
GFR_0   : 90   : Baseline GFR (mL/min/1.73m2)
PTH_0   : 40   : Baseline PTH (pg/mL)
VitD3_0 : 35   : Baseline 1,25(OH)2D (pg/mL)

// ── Calcium physiology parameters ──
k_Ca_abs   : 0.30  : Fractional intestinal Ca absorption (baseline)
k_Ca_bone  : 50    : Ca from bone resorption (mg/day)
Emax_PTH_Ca : 0.20 : PTH effect on Ca absorption (Emax)
EC50_PTH_Ca : 60   : PTH EC50 for Ca absorption (pg/mL)
CL_Ca_renal : 0.02 : Renal Ca clearance (fraction of filtered)
k_Ca_out    : 1.5  : Ca plasma elimination rate (1/day)

// ── Oxalate physiology parameters ──
k_Ox_hepatic : 15    : Baseline hepatic oxalate synthesis (mg/day)
k_Ox_gut_abs : 0.10  : Fractional dietary oxalate absorption
k_Ox_PH1     : 5.0   : Multiplier for PH1 hepatic OX synthesis
k_Ox_IBD     : 2.0   : Multiplier for IBD gut OX absorption
k_Ox_elim    : 2.0   : Oxalate plasma elimination rate (1/day)
k_Ox_VitC    : 0.005 : Conversion rate VitC→oxalate (mg_OX/mg_VitC/day)

// ── Uric acid physiology parameters ──
k_UA_prod0  : 800   : Baseline UA production rate (mg/day)
k_UA_diet   : 0.003 : Diet contribution per g protein/day
k_UA_MetSyn : 1.3   : Metabolic syndrome multiplier
k_UA_elim   : 0.8   : UA plasma elimination (1/day)
f_UA_renal  : 0.70  : Fraction of UA excreted in urine
Emax_XO     : 0.90  : Allopurinol max XO inhibition
EC50_XO_Oxp : 2.0   : Oxypurinol EC50 for XO (mg/L)
Hill_XO     : 1.2   : Hill coefficient XO inhibition

// ── Citrate physiology parameters ──
k_Cit_base   : 600  : Baseline urinary citrate excretion (mg/day)
k_Cit_PTH    : 0.8  : PTH reduces citrate (fraction decrease per 2x PTH)
k_Cit_acid   : 0.5  : Metabolic acidosis reduces citrate
k_Cit_K      : 1.3  : Alkalinization / K load increases citrate
k_Cit_elim   : 3.0  : Citrate plasma elimination (1/day)

// ── Urine volume physiology ──
k_Urine_base : 1.5  : Baseline urine volume (L/day)
k_Urine_fluid: 0.85 : Fraction of excess fluid -> urine

// ── Supersaturation parameters ──
// CaOx SS = [(uCa/V) * (uOx/V)] / Ksp_CaOx_eff
Ksp_CaOx  : 2.57e-9 : CaOx solubility product (mol^2/L^2)
Ksp_CaP   : 2.35e-7 : Brushite Ksp
Ksp_UA    : 9.8e-5  : UA solubility (mol/L)
Mg_frac   : 0.25    : Mg inhibition factor on CaOx SS
Cit_frac  : 0.35    : Citrate chelation fraction of free Ca

// ── Stone growth parameters ──
k_stone_grow : 0.02  : CaOx stone growth rate constant (mm/yr per SS unit)
n_stone      : 1.5   : Stone growth SS exponent
k_stone_pass : 0.08  : Spontaneous passage rate constant (1/yr)
Stone_0      : 2.0   : Baseline stone size (mm)

// ── Renal injury parameters ──
k_infl_crystal : 0.10 : Crystal→inflammation induction rate
k_infl_decay   : 0.30 : Inflammation resolution rate
k_GFR_loss     : 0.005: GFR loss per unit inflammation-time (mL/min/yr)

// ── Thiazide PK parameters (HCTZ) ──
ka_HCTZ  : 1.5  : Absorption rate constant (1/h)
CL_HCTZ  : 18.0 : Total clearance (L/h)
V1_HCTZ  : 40.0 : Central volume (L)
V2_HCTZ  : 80.0 : Peripheral volume (L)
Q_HCTZ   : 10.0 : Intercompartmental clearance (L/h)
F_HCTZ   : 0.65 : Oral bioavailability
Emax_HCTZ: 0.45 : HCTZ maximum Ca reduction (Emax)
EC50_HCTZ: 0.05 : HCTZ EC50 on uCa reduction (mg/L)

// ── Allopurinol/Oxypurinol PK parameters ──
ka_Allo  : 2.0   : Allo absorption rate (1/h)
CL_Allo  : 10.0  : Allo clearance (L/h)
V_Allo   : 35.0  : Allo central volume (L)
F_Allo   : 0.90  : Allo bioavailability
k_Allo_Oxp: 0.30 : Allo → Oxypurinol conversion rate (1/h)
CL_Oxp   : 2.5   : Oxypurinol renal clearance (L/h)
V_Oxp    : 100.0 : Oxypurinol distribution volume (L)

// ── Potassium Citrate PK parameters ──
ka_KCit  : 1.2   : KCit absorption rate (1/h)
CL_KCit  : 15.0  : KCit clearance (L/h)
V_KCit   : 25.0  : KCit volume (L)
F_KCit   : 0.95  : KCit bioavailability
Emax_KCit: 0.60  : Max urinary citrate increase (fraction)
EC50_KCit: 0.3   : KCit EC50 (mg/L)

// ── Tamsulosin PK parameters ──
ka_Tam   : 0.8   : Tamsulosin absorption (1/h)
CL_Tam   : 3.5   : Tamsulosin clearance (L/h)
V_Tam    : 65.0  : Tamsulosin volume (L)
F_Tam    : 0.90  : Bioavailability
Emax_Tam : 0.28  : Max stone passage rate increase
EC50_Tam : 0.001 : Tamsulosin EC50 (mg/L)

// Dosing (mg/day for drugs, 0=off)
Dose_HCTZ : 0   : HCTZ daily dose (mg/day)
Dose_Allo  : 0   : Allopurinol daily dose (mg/day)
Dose_KCit  : 0   : K-citrate daily dose (mEq/day, as mEq)
Dose_Tam   : 0   : Tamsulosin daily dose (mg/day)
Lumasiran  : 0   : Lumasiran (1=on, PH1 treatment)

$CMT @annotated
// Drug PK compartments
HCTZ1   : HCTZ gut depot (mg)
HCTZ2   : HCTZ central (mg)
HCTZ3   : HCTZ peripheral (mg)
ALLO1   : Allopurinol gut depot (mg)
ALLO2   : Allopurinol central (mg)
OXP     : Oxypurinol central (mg)
KCIT1   : K-Citrate gut depot (mEq)
KCIT2   : K-Citrate central (mEq)
TAM1    : Tamsulosin gut depot (mg)
TAM2    : Tamsulosin central (mg)
// Physiology compartments
CA_PL   : Plasma Ca pool (mg)
OX_PL   : Plasma Oxalate pool (mg)
UA_PL   : Plasma Uric Acid pool (mg)
CIT_PL  : Plasma Citrate pool (mg)
// Stone compartment
STONE   : Stone size (mm)
// Inflammation/injury
INFL    : Renal inflammation index (AU)
// GFR trajectory
GFR_CUR : Current GFR (mL/min)

$INIT
HCTZ1  = 0, HCTZ2  = 0, HCTZ3  = 0
ALLO1  = 0, ALLO2  = 0, OXP    = 0
KCIT1  = 0, KCIT2  = 0
TAM1   = 0, TAM2   = 0
CA_PL  = 840    // ~10.5 mg/dL × 80L plasma → actual in mg total
OX_PL  = 1.2    // ~15 µmol/L plasma oxalate pool (mg)
UA_PL  = 400    // ~5 mg/dL × 80L
CIT_PL = 12     // ~0.15 mmol/L
STONE  = 2.0    // starting stone mm
INFL   = 0.05   // baseline low inflammation
GFR_CUR= 90     // baseline GFR

$OMEGA @block 0
@labels ETA_CA ETA_OX ETA_UA ETA_Stone
0.04
0.02 0.09
0.02 0.03 0.04
0.00 0.00 0.00 0.00

$SIGMA @block 0
@labels EPS_uCa EPS_uOx EPS_uUA
0.04 0.02 0.04

$MAIN
// ── Urine volume (L/day → L/h) ──
double UV_day = k_Urine_base + k_Urine_fluid * (FluIn - 1.5);
if (UV_day < 0.5) UV_day = 0.5;
double UV_h = UV_day / 24.0;  // L/h

// ── PTH dynamics (simplified algebraic) ──
double PTH = PTH_0 * (1.0 + HPT * 2.0);   // HPT doubles PTH

// ── VitD3 (1,25(OH)2D, pg/mL) ──
double VitD3 = VitD3_0 * (PTH / PTH_0) * 0.8;

// ── Drug concentrations ──
double C_HCTZ = HCTZ2 / V1_HCTZ;    // mg/L
double C_ALLO = ALLO2 / V_Allo;      // mg/L
double C_OXP  = OXP   / V_Oxp;      // mg/L
double C_KCIT = KCIT2 / V_KCit;     // mEq/L (approx as mg/L equivalent)
double C_TAM  = TAM2  / V_Tam;      // mg/L

// ── Daily dose rates (mg/h or mEq/h) for infusion-style ──
// Doses given once daily (approximated as continuous for simplicity)
double Rate_HCTZ = Dose_HCTZ / 24.0 * F_HCTZ;
double Rate_Allo = Dose_Allo / 24.0 * F_Allo;
double Rate_KCit = Dose_KCit / 24.0 * F_KCit;
double Rate_Tam  = Dose_Tam  / 24.0 * F_Tam;

// ── XO inhibition by oxypurinol ──
double XO_inhib = Emax_XO * pow(C_OXP, Hill_XO) /
                  (pow(EC50_XO_Oxp, Hill_XO) + pow(C_OXP, Hill_XO));

// ── HCTZ effect on urinary Ca ──
double HCTZ_Ca_effect = Emax_HCTZ * C_HCTZ / (EC50_HCTZ + C_HCTZ);

// ── KCit effect on urinary citrate ──
double KCit_Cit_effect = Emax_KCit * C_KCIT / (EC50_KCit + C_KCIT);

// ── Lumasiran effect on hepatic oxalate synthesis ──
double Luma_OX_inhib = Lumasiran * 0.53;  // 53% reduction (ILLUMINATE-A)

// ── Calcium: intestinal absorption ──
double Ca_abs_frac = k_Ca_abs * (1.0 + Emax_PTH_Ca * PTH / (EC50_PTH_Ca + PTH));
double Ca_input    = Ca_diet * Ca_abs_frac / 24.0 +   // absorbed Ca (mg/h)
                     k_Ca_bone / 24.0;                  // bone resorption

// ── Oxalate input ──
double Ox_hepatic = k_Ox_hepatic * (1.0 + PH1 * (k_Ox_PH1 - 1.0)) *
                    (1.0 - Luma_OX_inhib) / 24.0;
double Ox_gut_abs = (Ox_diet * k_Ox_gut_abs * (1.0 + IBD * (k_Ox_IBD - 1.0)) +
                     VitC * k_Ox_VitC) / 24.0;

// ── UA production ──
double UA_prod = k_UA_prod0 * (1.0 + k_UA_diet * Pro_diet / 24.0 *
                 (k_UA_MetSyn - 1.0) * MetSyn) * (1.0 - XO_inhib) / 24.0;

// ── Citrate input (renal excretion modeled via pool) ──
double Cit_input = k_Cit_base * (1.0 - k_Cit_PTH * (PTH/PTH_0 - 1.0) * 0.3) /24.0;

// IIV on key parameters
double Ca_PL_mg_dL  = CA_PL / 80.0 * (1.0 + ETA(1));
double Ox_PL_umolL  = OX_PL / 0.088 * (1.0 + ETA(2));  // MW oxalate=88
double UA_PL_mg_dL  = UA_PL / 80.0 * (1.0 + ETA(3));

// ── Urinary excretion rates (mg/h) ──
double uCa_rate = (CA_PL / 80.0 / 100.0) * GFR_CUR * 60.0 * CL_Ca_renal *
                  (1.0 - HCTZ_Ca_effect) * (1.0 + Na_diet / 150.0 * 0.1);
// Convert: [Ca] mg/dL × GFR L/h × CL fraction → mg/h

double uOx_rate  = OX_PL * k_Ox_elim * 0.6;  // ~60% renally excreted
double uUA_rate  = UA_PL * k_UA_elim * f_UA_renal / 24.0;
double uCit_rate = CIT_PL * k_Cit_elim * (1.0 + KCit_Cit_effect) / 24.0;

// ── 24h urinary excretions (mg/day equivalent for SS calculation) ──
double uCa_24  = uCa_rate  * 24.0;
double uOx_24  = uOx_rate  * 24.0;
double uUA_24  = uUA_rate  * 24.0;
double uCit_24 = uCit_rate * 24.0;

// ── Supersaturation (CaOx) ──
// [Ca] mmol/L = uCa_24/400.0 / UV_day; [Ox] mmol/L = uOx_24/88.0 / UV_day
// CaOx SS = ([Ca] * [Ox]) / Ksp_CaOx (dimensionless)
double conc_Ca_mM = uCa_24 / 400.0 / UV_day;   // mmol/L in urine
double conc_Ox_mM = uOx_24 / 88.0  / UV_day;
double conc_UA_mM = uUA_24 / 168.0 / UV_day;
double conc_Cit_mM= uCit_24/ 192.0 / UV_day;

// Effective free Ca (reduced by citrate chelation)
double freeCa_mM  = conc_Ca_mM * (1.0 - Cit_frac * conc_Cit_mM /
                    (conc_Cit_mM + 2.0)) * (1.0 - Mg_frac * 0.3);

double SS_CaOx = (freeCa_mM * 1e-3) * (conc_Ox_mM * 1e-3) / Ksp_CaOx;
double SS_CaP  = (freeCa_mM * 1e-3) * (0.8e-3) / Ksp_CaP;  // simplified PO4
double urine_pH = 5.8 - MetSyn * 0.5 + KCit_Cit_effect * 1.2;
double SS_UA   = conc_UA_mM * 1e-3 / (Ksp_UA * pow(10.0, urine_pH - 5.35));

// Clamp SS ≥ 0
if(SS_CaOx < 0) SS_CaOx = 0;
if(SS_CaP  < 0) SS_CaP  = 0;
if(SS_UA   < 0) SS_UA   = 0;

// ── Tamsulosin effect on passage ──
double Tam_pass = Emax_Tam * C_TAM / (EC50_Tam + C_TAM);
double k_passage = k_stone_pass * (1.0 + Tam_pass);

// ── Stone growth driver ──
double SS_eff = (SS_CaOx > 1.0) ? pow(SS_CaOx - 1.0, n_stone) : 0.0;
double dStone_grow = k_stone_grow * SS_eff * 24.0;  // mm/yr → mm/day (×24 for h)
double dStone_pass = k_passage / 365.0 * STONE;     // passage per day

// ── Crystal load for inflammation ──
double crystal_load = (SS_CaOx > 1.0) ? (SS_CaOx - 1.0) : 0.0;

$ODE
// ── Drug PK ODEs ──
// HCTZ (2-compartment, first-order absorption)
dxdt_HCTZ1 = Rate_HCTZ - ka_HCTZ * HCTZ1;
dxdt_HCTZ2 = ka_HCTZ * HCTZ1 - (CL_HCTZ/V1_HCTZ) * HCTZ2 -
              (Q_HCTZ/V1_HCTZ) * HCTZ2 + (Q_HCTZ/V2_HCTZ) * HCTZ3;
dxdt_HCTZ3 = (Q_HCTZ/V1_HCTZ) * HCTZ2 - (Q_HCTZ/V2_HCTZ) * HCTZ3;

// Allopurinol → Oxypurinol
dxdt_ALLO1 = Rate_Allo - ka_Allo * ALLO1;
dxdt_ALLO2 = ka_Allo * ALLO1 - (CL_Allo/V_Allo) * ALLO2 -
              k_Allo_Oxp * ALLO2;
dxdt_OXP   = k_Allo_Oxp * ALLO2 - (CL_Oxp/V_Oxp) * OXP;

// K-Citrate
dxdt_KCIT1 = Rate_KCit - ka_KCit * KCIT1;
dxdt_KCIT2 = ka_KCit * KCIT1 - (CL_KCit/V_KCit) * KCIT2;

// Tamsulosin
dxdt_TAM1 = Rate_Tam - ka_Tam * TAM1;
dxdt_TAM2 = ka_Tam * TAM1 - (CL_Tam/V_Tam) * TAM2;

// ── Physiological ODEs ──
// Calcium plasma pool (mg; pool approximates total body Ca flux)
dxdt_CA_PL = Ca_input - k_Ca_out * CA_PL / 80.0 -
              uCa_rate;

// Oxalate plasma pool (mg)
dxdt_OX_PL = Ox_hepatic + Ox_gut_abs - k_Ox_elim * OX_PL;

// Uric Acid plasma pool (mg)
dxdt_UA_PL = UA_prod - k_UA_elim * UA_PL / 24.0;

// Citrate plasma pool (mg)
dxdt_CIT_PL = Cit_input - k_Cit_elim * CIT_PL / 24.0;

// Stone growth (mm): driven by CaOx SS
dxdt_STONE = (dStone_grow - dStone_pass) / (24.0 * 365.0);
// [k_stone_grow mm/yr per SS unit → /h conversion]
// Simplified: net stone growth per hour
if(STONE < 0) STONE = 0;

// Renal inflammation (AU)
dxdt_INFL = k_infl_crystal * crystal_load - k_infl_decay * INFL;

// GFR trajectory (mL/min)
dxdt_GFR_CUR = -k_GFR_loss * INFL;
if(GFR_CUR < 15) GFR_CUR = 15;  // lower bound ESRD

$TABLE
// Observed quantities with residual error
capture uCa_mgday  = uCa_24  * (1.0 + EPS(1));
capture uOx_mgday  = uOx_24  * (1.0 + EPS(2));
capture uUA_mgday  = uUA_24  * (1.0 + EPS(3));
capture uCit_mgday = uCit_24;
capture uVol_Lday  = UV_day;
capture urinepH    = urine_pH;

capture SS_CaOx_out = SS_CaOx;
capture SS_CaP_out  = SS_CaP;
capture SS_UA_out   = SS_UA;

capture StoneSize_mm = STONE;
capture GFR_out      = GFR_CUR;
capture Inflammation = INFL;

capture C_HCTZ_mgL   = C_HCTZ;
capture C_OXP_mgL    = C_OXP;
capture C_KCIT_mEqL  = C_KCIT;
capture C_TAM_mgL    = C_TAM;

// Flag risk thresholds
capture Hypercalciuria = (uCa_24 > (SEX ? 300 : 250)) ? 1.0 : 0.0;
capture Hyperoxaluria  = (uOx_24 > 45)  ? 1.0 : 0.0;
capture Hyperuricosuria= (uUA_24 > (SEX ? 800 : 750)) ? 1.0 : 0.0;
capture Hypocitraturia = (uCit_24 < 320) ? 1.0 : 0.0;
capture LowUV         = (UV_day < 2.0)   ? 1.0 : 0.0;

capture Stone_recur_risk = SS_CaOx * 0.4 + Hypercalciuria * 0.2 +
                           Hyperoxaluria * 0.2 + Hypocitraturia * 0.15 +
                           LowUV * 0.05;
'

# ─────────────────────────────────────────────────────────────────────────────
# Compile model
# ─────────────────────────────────────────────────────────────────────────────
mod <- mcode("urolithiasis_qsp", uri_model_code)

# ─────────────────────────────────────────────────────────────────────────────
# Helper: run simulation
# ─────────────────────────────────────────────────────────────────────────────
run_sim <- function(model, params = list(), duration_yr = 5, n = 1) {
  tmax_h <- duration_yr * 365 * 24
  times   <- seq(0, tmax_h, length.out = 500)

  if (length(params) > 0) model <- param(model, params)

  mrgsim(model, end = tmax_h, delta = tmax_h / 500, nid = n) %>%
    as_tibble() %>%
    mutate(time_yr = time / (365 * 24))
}

# ─────────────────────────────────────────────────────────────────────────────
# SCENARIO 1: Untreated idiopathic CaOx stone former
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== Scenario 1: Untreated CaOx stone former ===\n")
scen1 <- run_sim(mod, params = list(
  Ca_diet = 1000, Ox_diet = 250, FluIn = 1.8, Pro_diet = 120,
  Na_diet = 200, GFR_0 = 90, PTH_0 = 40
))

cat(sprintf("  5-yr final stone size: %.2f mm\n",
            tail(scen1$StoneSize_mm, 1)))
cat(sprintf("  Mean CaOx SS: %.2f\n", mean(scen1$SS_CaOx_out)))
cat(sprintf("  Mean uCa: %.0f mg/day\n", mean(scen1$uCa_mgday)))

# ─────────────────────────────────────────────────────────────────────────────
# SCENARIO 2: Thiazide treatment (HCTZ 25mg/day)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== Scenario 2: HCTZ 25mg/day ===\n")
scen2 <- run_sim(mod, params = list(
  Ca_diet = 1000, Ox_diet = 250, FluIn = 1.8, Pro_diet = 120,
  Na_diet = 200, GFR_0 = 90, PTH_0 = 40,
  Dose_HCTZ = 25
))

cat(sprintf("  5-yr final stone size: %.2f mm\n",
            tail(scen2$StoneSize_mm, 1)))
cat(sprintf("  uCa reduction: %.0f → %.0f mg/day\n",
            mean(scen1$uCa_mgday), mean(scen2$uCa_mgday)))

# ─────────────────────────────────────────────────────────────────────────────
# SCENARIO 3: Potassium citrate (60 mEq/day)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== Scenario 3: K-Citrate 60 mEq/day ===\n")
scen3 <- run_sim(mod, params = list(
  Ca_diet = 1000, Ox_diet = 250, FluIn = 1.8, Pro_diet = 120,
  Na_diet = 200, GFR_0 = 90, PTH_0 = 40,
  Dose_KCit = 60
))

cat(sprintf("  5-yr final stone size: %.2f mm\n",
            tail(scen3$StoneSize_mm, 1)))
cat(sprintf("  uCit increase: %.0f → %.0f mg/day\n",
            mean(scen1$uCit_mgday), mean(scen3$uCit_mgday)))

# ─────────────────────────────────────────────────────────────────────────────
# SCENARIO 4: Allopurinol 300mg/day (UA stone / CaOx with hyperuricosuria)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== Scenario 4: Allopurinol 300mg/day ===\n")
scen4 <- run_sim(mod, params = list(
  Ca_diet = 1000, Ox_diet = 150, FluIn = 2.0,
  Pro_diet = 150, Na_diet = 160, MetSyn = 1, GFR_0 = 85,
  Dose_Allo = 300
))

cat(sprintf("  5-yr uUA: %.0f mg/day (baseline ~%.0f)\n",
            mean(scen4$uUA_mgday),
            mean(run_sim(mod, params=list(MetSyn=1,Pro_diet=150))$uUA_mgday)))
cat(sprintf("  UA SS reduction: %.2f → %.2f\n",
            mean(run_sim(mod, params=list(MetSyn=1))$SS_UA_out),
            mean(scen4$SS_UA_out)))

# ─────────────────────────────────────────────────────────────────────────────
# SCENARIO 5: Lumasiran for Primary Hyperoxaluria Type 1
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== Scenario 5: PH1 — Lumasiran ===\n")
scen5_no_Rx <- run_sim(mod, params = list(PH1 = 1, FluIn = 3.0))
scen5_Luma  <- run_sim(mod, params = list(PH1 = 1, FluIn = 3.0,
                                           Lumasiran = 1))

cat(sprintf("  PH1 untreated uOx: %.0f mg/day\n",  mean(scen5_no_Rx$uOx_mgday)))
cat(sprintf("  PH1 + Lumasiran uOx: %.0f mg/day\n", mean(scen5_Luma$uOx_mgday)))
cat(sprintf("  CaOx SS reduction: %.2f → %.2f\n",
            mean(scen5_no_Rx$SS_CaOx_out),
            mean(scen5_Luma$SS_CaOx_out)))

# ─────────────────────────────────────────────────────────────────────────────
# SCENARIO 6: Comprehensive lifestyle + triple therapy
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== Scenario 6: Lifestyle + HCTZ + K-Citrate ===\n")
scen6 <- run_sim(mod, params = list(
  Ca_diet = 1000, Ox_diet = 100, FluIn = 3.0, Pro_diet = 80,
  Na_diet = 80, GFR_0 = 90,
  Dose_HCTZ = 25, Dose_KCit = 60
))

cat(sprintf("  5-yr stone size: %.2f mm (vs. %.2f untreated)\n",
            tail(scen6$StoneSize_mm, 1), tail(scen1$StoneSize_mm, 1)))
cat(sprintf("  CaOx SS: %.2f (vs. %.2f untreated)\n",
            mean(scen6$SS_CaOx_out), mean(scen1$SS_CaOx_out)))

# ─────────────────────────────────────────────────────────────────────────────
# SCENARIO 7: Metabolic Syndrome UA stone + Allopurinol + K-Citrate
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== Scenario 7: MetSyn UA stone — Combo Rx ===\n")
scen7 <- run_sim(mod, params = list(
  MetSyn = 1, FluIn = 2.5, Dose_Allo = 300, Dose_KCit = 60, Pro_diet = 120
))
cat(sprintf("  UA SS under combo Rx: %.2f\n", mean(scen7$SS_UA_out)))
cat(sprintf("  Urine pH: %.1f\n", mean(scen7$urinepH)))

# ─────────────────────────────────────────────────────────────────────────────
# VISUALIZATION
# ─────────────────────────────────────────────────────────────────────────────
scenarios_list <- list(
  "1. Untreated"          = scen1,
  "2. HCTZ 25mg"          = scen2,
  "3. K-Citrate 60mEq"   = scen3,
  "4. Allopurinol 300mg" = scen4,
  "5. PH1+Lumasiran"     = scen5_Luma,
  "6. Lifestyle+Combo"   = scen6
)

df_all <- bind_rows(lapply(names(scenarios_list), function(nm) {
  scenarios_list[[nm]] %>%
    select(time_yr, StoneSize_mm, SS_CaOx_out, uCa_mgday, uCit_mgday,
           uOx_mgday, GFR_out, Inflammation) %>%
    mutate(Scenario = nm)
}))

# Plot 1: Stone size trajectories
p1 <- ggplot(df_all, aes(time_yr, StoneSize_mm, color = Scenario)) +
  geom_line(size = 1.1) +
  labs(title = "Urolithiasis QSP: Stone Size Trajectories",
       x = "Time (years)", y = "Stone Size (mm)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 13) +
  scale_color_brewer(palette = "Set2")

# Plot 2: CaOx Supersaturation
p2 <- ggplot(df_all, aes(time_yr, SS_CaOx_out, color = Scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  annotate("text", x = 0.5, y = 1.1, label = "SS = 1 (crystallization threshold)",
           color = "red", size = 3) +
  labs(title = "CaOx Supersaturation by Treatment Scenario",
       x = "Time (years)", y = "CaOx SS (relative)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 13) +
  scale_color_brewer(palette = "Set2")

# Plot 3: Urinary Ca over time
p3 <- ggplot(df_all %>% filter(Scenario %in% c("1. Untreated", "2. HCTZ 25mg",
                                                "6. Lifestyle+Combo")),
             aes(time_yr, uCa_mgday, color = Scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 300, linetype = "dashed", color = "orange",
             alpha = 0.7) +
  annotate("text", x = 0.5, y = 310, label = "Hypercalciuria threshold (M)",
           color = "darkorange", size = 3) +
  labs(title = "24h Urinary Calcium Under Thiazide Therapy",
       x = "Time (years)", y = "Urinary Ca (mg/day)") +
  theme_bw(base_size = 13) +
  scale_color_brewer(palette = "Dark2")

# Plot 4: GFR trajectory
p4 <- ggplot(df_all, aes(time_yr, GFR_out, color = Scenario)) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = c(60, 30, 15), linetype = "dashed",
             color = c("yellow3", "orange", "red")) +
  annotate("text", x = 4.5, y = c(62, 32, 17),
           label = c("CKD G3a", "CKD G4", "CKD G5"),
           color = c("darkgoldenrod", "darkorange", "red"), size = 3) +
  labs(title = "GFR Trajectory by Treatment Scenario",
       x = "Time (years)", y = "eGFR (mL/min/1.73m²)") +
  theme_bw(base_size = 13) +
  scale_color_brewer(palette = "Set2")

# Arrange plots
if (requireNamespace("gridExtra", quietly = TRUE)) {
  library(gridExtra)
  gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
} else {
  print(p1); print(p2); print(p3); print(p4)
}

# ─────────────────────────────────────────────────────────────────────────────
# SENSITIVITY ANALYSIS: Fluid intake vs. stone risk
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== Sensitivity: Fluid Intake vs. CaOx SS ===\n")
fluid_vals <- seq(1.0, 4.0, by = 0.5)
ss_vals <- sapply(fluid_vals, function(fl) {
  s <- run_sim(mod, params = list(FluIn = fl))
  mean(s$SS_CaOx_out)
})

df_sens <- data.frame(Fluid_L = fluid_vals, SS_CaOx = ss_vals)
print(df_sens)

p_sens <- ggplot(df_sens, aes(Fluid_L, SS_CaOx)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(size = 3, color = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = "Sensitivity Analysis: Fluid Intake vs. CaOx Supersaturation",
       x = "Fluid Intake (L/day)", y = "Mean CaOx SS") +
  theme_bw(base_size = 13)

print(p_sens)

cat("\n=== QSP Model Simulation Complete ===\n")
