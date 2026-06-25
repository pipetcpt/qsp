## =============================================================================
## Polycythemia Vera (PV) — QSP Model in mrgsolve
## =============================================================================
## Disease: Polycythemia Vera (진성 다혈증)
## JAK2 V617F mutation-driven myeloproliferative neoplasm
##
## Key references:
##   Tefferi et al., NEJM 2013, 368:22-33
##   Senyak et al., Leukemia 2013
##   Verstovsek et al., NEJM 2012, 366:799-807 (RESPONSE trial)
##   Kiladjian et al., NEJM 2015 (PROUD-PV)
##   Markovtsov et al., J Biol Chem 2016 (ruxolitinib PK/PD)
##   NCCN Guidelines PV v2.2024
##
## Model compartments (19 ODEs):
##   PK: Ruxolitinib (central, peripheral)
##       Hydroxyurea (central)
##       PEG-IFN-α2a (SC depot, central)
##   PD/Hematopoiesis:
##       Erythroid progenitors (BFU-E, CFU-E)
##       Reticulocytes (BM + circ)
##       RBC mass
##       Platelet pool
##       WBC/Neutrophil
##       Spleen volume
##       JAK2 V617F allele burden
##       Bone marrow fibrosis score
## =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

## ---------------------------------------------------------------------------
## mrgsolve Model Code
## ---------------------------------------------------------------------------

pv_code <- '
$PROB
Polycythemia Vera QSP Model
JAK2 V617F · Ruxolitinib · Hydroxyurea · PEG-IFN-α2a · Phlebotomy

$PARAM @annotated
// --- Patient characteristics ---
BWT       : 75    : Body weight (kg)
AGE       : 60    : Age (years)
SEX       : 1     : Sex (1=male, 0=female)
eGFR      : 80    : eGFR (mL/min/1.73m2)

// --- JAK2 V617F disease parameters ---
JAK2_AB0  : 50    : Baseline JAK2 V617F allele burden (%)
k_clone   : 0.003 : Clone expansion rate (/day)
k_clone_d : 0.001 : Spontaneous clone regression rate (/day)
Emax_JAK2 : 1.5   : Maximum proliferation fold due to JAK2 mutation

// --- Normal hematopoiesis baseline parameters ---
BFUE_0    : 100   : BFU-E baseline (arbitrary units)
CFUE_0    : 100   : CFU-E baseline
RETIC_BM0 : 100   : BM reticulocyte baseline
RETIC_C0  : 100   : Circulating reticulocyte baseline
RBC_0     : 100   : RBC mass baseline (normalized to 100%)
PLT_0     : 400   : Platelet count baseline (×10^9/L)
WBC_0     : 8     : WBC baseline (×10^9/L)

// --- Erythroid progenitor kinetics ---
k_BFUE_in : 0.03  : BFU-E input rate (/day)
k_BFUE_d  : 0.05  : BFU-E differentiation rate (/day)
k_CFUE_d  : 0.10  : CFU-E differentiation rate (/day)
k_RETIC_BM: 0.15  : BM reticulocyte maturation rate (/day)
k_RETIC_C : 0.50  : Circulating reticulocyte maturation rate (/day)
k_RBC_d   : 0.0083: RBC elimination rate (/day, ~120d lifespan)

// --- Platelet kinetics ---
k_PLT_in  : 0.12  : Platelet production rate
k_PLT_d   : 0.10  : Platelet elimination (/day, ~10d lifespan)

// --- WBC/Neutrophil kinetics ---
k_WBC_in  : 0.33  : WBC production rate
k_WBC_d   : 0.33  : WBC elimination (/day, ~3d lifespan)

// --- Spleen parameters ---
SPL_0     : 450   : Spleen volume baseline (mL)
k_SPL_in  : 0.002 : Spleen EMH expansion rate (/day)
k_SPL_d   : 0.005 : Spleen natural regression rate (/day)
SPL_max   : 3000  : Maximum spleen volume (mL)

// --- BM fibrosis parameters ---
FIB_0     : 0.1   : BM reticulin fibrosis baseline (grade 0-3 scale)
k_FIB_in  : 0.0005: Fibrosis progression rate (/day)
k_FIB_d   : 0.0002: Fibrosis regression rate (/day)

// --- Ruxolitinib PK parameters ---
KA_RUX    : 2.4   : Ruxolitinib absorption rate (/h)
F_RUX     : 0.95  : Ruxolitinib oral bioavailability
VC_RUX    : 72    : Ruxolitinib central volume (L)
VP_RUX    : 28    : Ruxolitinib peripheral volume (L)
CL_RUX    : 22    : Ruxolitinib clearance (L/h)
QP_RUX    : 16    : Ruxolitinib inter-compartmental CL (L/h)

// --- Ruxolitinib PD parameters ---
IC50_JAK1 : 3.3   : Ruxolitinib IC50 for JAK1 (nM → converted to ng/mL: 1 nM = 0.306 ng/mL)
IC50_JAK2 : 2.8   : Ruxolitinib IC50 for JAK2 (nM → ~0.86 ng/mL)
Imax_RUX  : 1.0   : Maximum inhibition fraction by ruxolitinib
HILL_RUX  : 1.0   : Hill coefficient for ruxolitinib

// --- Hydroxyurea PK parameters ---
KA_HYU    : 1.5   : Hydroxyurea absorption rate (/h)
F_HYU     : 0.80  : Hydroxyurea oral bioavailability
VC_HYU    : 20    : Hydroxyurea central volume (L)
CL_HYU    : 5.0   : Hydroxyurea clearance (L/h)

// --- Hydroxyurea PD parameters ---
IC50_HYU  : 150   : Hydroxyurea IC50 for erythropoiesis (μM)
Imax_HYU  : 0.85  : Maximum inhibition by hydroxyurea

// --- PEG-IFN-α2a PK parameters ---
KA_IFN    : 0.03  : SC absorption rate (/h)
F_IFN     : 0.84  : SC bioavailability
VC_IFN    : 6.0   : Central volume (L)
CL_IFN    : 0.07  : Clearance (L/h) → t1/2 ~80h
MW_IFN    : 60000 : Molecular weight of PEG-IFN (Da)

// --- PEG-IFN-α2a PD parameters ---
EC50_IFN  : 50    : IFN EC50 for clone suppression (pg/mL)
Emax_IFN  : 0.70  : Maximum allele burden reduction fraction

// --- Phlebotomy parameters ---
PHLEBOT_VOL : 450 : Volume removed per phlebotomy (mL)
RBC_conc    : 5.0 : RBC concentration (×10^12/L)
Hct_target  : 45  : Target hematocrit (%)

// --- Thrombosis risk model parameters ---
AGE_RISK  : 0.02  : Age contribution to thrombosis risk (/y above 60)
HIST_RISK : 0.15  : History of thrombosis risk factor

$PARAM @annotated
// Scenario switches (0=off, 1=on)
USE_RUX   : 0     : Ruxolitinib treatment (0/1)
USE_HYU   : 0     : Hydroxyurea treatment (0/1)
USE_IFN   : 0     : PEG-IFN-α2a treatment (0/1)
USE_ASP   : 1     : Aspirin (0/1)
PHLEBOT_FLAG : 0  : Phlebotomy event flag

$CMT @annotated
// Ruxolitinib PK
DEPOT_RUX : Ruxolitinib absorption depot (mg)
CENT_RUX  : Ruxolitinib central compartment (mg)
PERI_RUX  : Ruxolitinib peripheral compartment (mg)

// Hydroxyurea PK
CENT_HYU  : Hydroxyurea central compartment (mg)

// PEG-IFN-α2a PK
SC_IFN    : PEG-IFN-α2a SC depot (μg)
CENT_IFN  : PEG-IFN-α2a central compartment (μg)

// Erythroid hematopoiesis
BFUE      : BFU-E progenitor pool
CFUE      : CFU-E progenitor pool
RETIC_BM  : BM reticulocytes
RETIC_C   : Circulating reticulocytes
RBC       : RBC mass (normalized)

// Thrombopoiesis & leukopoiesis
PLT       : Platelet pool (×10^9/L)
WBC       : WBC pool (×10^9/L)

// Spleen & fibrosis
SPL       : Spleen volume (mL)
FIBRO     : BM fibrosis score (0-3)

// JAK2 allele burden
ALLELE    : JAK2 V617F allele burden (%)

$INIT @annotated
DEPOT_RUX = 0
CENT_RUX  = 0
PERI_RUX  = 0
CENT_HYU  = 0
SC_IFN    = 0
CENT_IFN  = 0
BFUE      = 100   // Elevated due to PV at steady state
CFUE      = 100
RETIC_BM  = 100
RETIC_C   = 100
RBC       = 130   // Elevated RBC mass in PV (~130% of normal)
PLT       = 600   // Thrombocytosis common in PV
WBC       = 14    // Leukocytosis
SPL       = 800   // Moderate splenomegaly at baseline
FIBRO     = 0.5   // Mild reticulin fibrosis at baseline
ALLELE    = 50    // Baseline allele burden

$ODE
// ============================================================
// Ruxolitinib PK
// ============================================================
double CP_RUX = CENT_RUX / VC_RUX * 1000;  // ng/mL (MW ~306 g/mol)

double dDEPOT_RUX = -KA_RUX * DEPOT_RUX;
double dCENT_RUX  = KA_RUX * DEPOT_RUX * F_RUX
                   - (CL_RUX + QP_RUX) * (CENT_RUX / VC_RUX)
                   + QP_RUX * (PERI_RUX / VP_RUX);
double dPERI_RUX  = QP_RUX * (CENT_RUX / VC_RUX) - QP_RUX * (PERI_RUX / VP_RUX);

// ============================================================
// Hydroxyurea PK
// ============================================================
double CP_HYU = CENT_HYU / VC_HYU * 1000;  // μM (MW ~76 g/mol → 1 mg/L ≈ 13.2 μM)
double CP_HYU_uM = CP_HYU * 13.2;

double dCENT_HYU = KA_HYU * F_HYU  // simplified: assuming depot pre-absorbed
                  - CL_HYU * (CENT_HYU / VC_HYU);

// ============================================================
// PEG-IFN-α2a PK
// ============================================================
double CP_IFN = CENT_IFN / VC_IFN * 1e6;  // pg/mL

double dSC_IFN   = -KA_IFN * SC_IFN;
double dCENT_IFN = KA_IFN * F_IFN * SC_IFN
                  - CL_IFN * (CENT_IFN / VC_IFN);

// ============================================================
// Drug Effect (Inhibition Functions)
// ============================================================
// Ruxolitinib: Imax × Cp^Hill / (IC50^Hill + Cp^Hill)
double IC50_RUX_ngmL = IC50_JAK2 * 0.306;  // nM → ng/mL (MW=306)
double INH_RUX = USE_RUX * Imax_RUX * pow(CP_RUX, HILL_RUX) /
                 (pow(IC50_RUX_ngmL, HILL_RUX) + pow(CP_RUX, HILL_RUX));

// Hydroxyurea: inhibits erythroid & platelet progenitor proliferation
double INH_HYU = USE_HYU * Imax_HYU * CP_HYU_uM / (IC50_HYU + CP_HYU_uM);

// PEG-IFN: reduces JAK2+ clone/allele burden
double EFF_IFN = USE_IFN * Emax_IFN * CP_IFN / (EC50_IFN + CP_IFN);

// Combined cytoreductive inhibition of proliferation
double INH_PROLIF = 1.0 - fmax(INH_RUX, INH_HYU);

// ============================================================
// JAK2 V617F Allele Burden ODE
// ============================================================
// Allele burden grows via clone expansion, suppressed by IFN
double k_allele_expansion = k_clone * (1.0 - ALLELE / 100.0);
double k_allele_suppression = k_clone_d + EFF_IFN * 0.01;
double dALLELE = k_allele_expansion * ALLELE - k_allele_suppression * ALLELE;
// Clamp allele burden to [0, 100]
if (ALLELE >= 100) dALLELE = fmin(0.0, dALLELE);
if (ALLELE <= 0)   dALLELE = fmax(0.0, dALLELE);

// JAK2-driven proliferative advantage (scaled by allele burden)
double JAK2_effect = 1.0 + (Emax_JAK2 - 1.0) * (ALLELE / 100.0);

// ============================================================
// BFU-E Progenitor
// ============================================================
// Input: steady-state balanced with EPO-independent JAK2 drive
// Drug: ruxolitinib and HYU both reduce proliferation
double BFUE_prolif = k_BFUE_in * JAK2_effect * INH_PROLIF;
double BFUE_diff   = k_BFUE_d * BFUE;
double dBFUE = BFUE_prolif * BFUE_0 - BFUE_diff;

// ============================================================
// CFU-E Progenitor
// ============================================================
double CFUE_input = k_BFUE_d * BFUE;
double CFUE_diff  = k_CFUE_d * CFUE;
double dCFUE = CFUE_input - CFUE_diff;

// ============================================================
// BM Reticulocytes
// ============================================================
double RETIC_BM_input = k_CFUE_d * CFUE;
double RETIC_BM_matur = k_RETIC_BM * RETIC_BM;
double dRETIC_BM = RETIC_BM_input - RETIC_BM_matur;

// ============================================================
// Circulating Reticulocytes
// ============================================================
double RETIC_C_input = k_RETIC_BM * RETIC_BM;
double RETIC_C_matur = k_RETIC_C * RETIC_C;
double dRETIC_C = RETIC_C_input - RETIC_C_matur;

// ============================================================
// RBC Mass (normalized, 100 = normal)
// ============================================================
// Phlebotomy effect: reduces RBC mass acutely
double PHLEBOT_rate = PHLEBOT_FLAG * (PHLEBOT_VOL * RBC_conc * 0.001) / 24.0;  // per hour → daily
double dRBC = k_RETIC_C * RETIC_C - k_RBC_d * RBC - PHLEBOT_rate;

// ============================================================
// Platelet Pool
// ============================================================
double PLT_prolif = k_PLT_in * JAK2_effect * INH_PROLIF * PLT_0;
double PLT_elim   = k_PLT_d * PLT;
// Aspirin: reduces platelet activation but not count production
double dPLT = PLT_prolif - PLT_elim;

// ============================================================
// WBC Pool
// ============================================================
double WBC_prolif = k_WBC_in * JAK2_effect * INH_PROLIF * WBC_0;
double WBC_elim   = k_WBC_d * WBC;
double dWBC = WBC_prolif - WBC_elim;

// ============================================================
// Spleen Volume (EMH driven by BM output overflow)
// ============================================================
// Spleen grows when BM is hypercellular (RBC > 110), shrinks with treatment
double BM_overflow = fmax(0.0, (RBC - 110.0) / 110.0);
double SPL_growth  = k_SPL_in * BM_overflow * (SPL_max - SPL) / SPL_max;
double SPL_regress = k_SPL_d * SPL * (INH_RUX + INH_HYU);
double dSPL = SPL_growth - SPL_regress;

// ============================================================
// BM Fibrosis Score (0-3, continuous)
// ============================================================
// Fibrosis driven by: splenomegaly, disease duration, JAK2 burden
double FIB_drive = k_FIB_in * (ALLELE / 50.0) * (SPL / 800.0);
double FIB_regress = k_FIB_d * FIBRO * (INH_RUX + EFF_IFN);
double dFIBRO = FIB_drive - FIB_regress;
if (FIBRO >= 3.0) dFIBRO = fmin(0.0, dFIBRO);
if (FIBRO <= 0.0) dFIBRO = fmax(0.0, dFIBRO);

// Assign ODEs
dxdt_DEPOT_RUX = dDEPOT_RUX;
dxdt_CENT_RUX  = dCENT_RUX;
dxdt_PERI_RUX  = dPERI_RUX;
dxdt_CENT_HYU  = dCENT_HYU;
dxdt_SC_IFN    = dSC_IFN;
dxdt_CENT_IFN  = dCENT_IFN;
dxdt_BFUE      = dBFUE;
dxdt_CFUE      = dCFUE;
dxdt_RETIC_BM  = dRETIC_BM;
dxdt_RETIC_C   = dRETIC_C;
dxdt_RBC       = dRBC;
dxdt_PLT       = dPLT;
dxdt_WBC       = dWBC;
dxdt_SPL       = dSPL;
dxdt_FIBRO     = dFIBRO;
dxdt_ALLELE    = dALLELE;

$TABLE
// --- Derived clinical variables ---
double Hct       = RBC * 0.35;                    // Hematocrit (%), normal ~42% male
double Hgb       = RBC * 0.115;                   // Hemoglobin (g/dL)
double EPO_serum = 6.0 * exp(-0.03 * (RBC - 100.0)); // Serum EPO (mU/mL), suppressed in PV

// Ruxolitinib concentrations
double Cp_RUX_ngmL = CP_RUX;
double Cp_RUX_nM   = CP_RUX / 0.306;

// JAK inhibition effect
double pSTAT5_inhib = USE_RUX * Imax_RUX * Cp_RUX_ngmL /
                      (IC50_RUX_ngmL + Cp_RUX_ngmL);

// Annual thrombosis risk (simplified model)
double THROMB_RISK = (Hct > 45 ? 0.06 : 0.02) +
                     (PLT > 1000 ? 0.04 : 0.0) +
                     (AGE > 60 ? (AGE - 60) * AGE_RISK : 0.0) +
                     HIST_RISK * 0.5;

// Spleen volume reduction from baseline (%)
double SVR = (800 - SPL) / 800.0 * 100.0;  // % reduction from 800 mL baseline
double SVR35_flag = SVR >= 35 ? 1.0 : 0.0;  // SVR35 achievement

// Complete Hematologic Response criteria
double CHR_Hct = Hct < 45.0 ? 1.0 : 0.0;
double CHR_PLT = PLT < 400.0 ? 1.0 : 0.0;
double CHR_WBC = WBC < 10.0  ? 1.0 : 0.0;
double CHR_SPL = SPL < 450.0 ? 1.0 : 0.0;
double CHR     = CHR_Hct * CHR_PLT * CHR_WBC * CHR_SPL;

// MPN Symptom score (proxy: 0-100, higher = worse)
double MPN_SAF_TSS = fmax(0, 20.0 + (Hct - 45.0) * 2.0 +
                          (PLT - 400.0) / 50.0 +
                          (SPL - 450.0) / 100.0 -
                          pSTAT5_inhib * 30.0);

// Post-PV MF transformation risk (annual, cumulative)
double MF_risk_annual = 0.005 + FIBRO * 0.01 + (ALLELE > 80 ? 0.01 : 0.0);

// Table output
capture Hct         = Hct;
capture Hgb         = Hgb;
capture EPO_serum   = EPO_serum;
capture PLT_out     = PLT;
capture WBC_out     = WBC;
capture SPL_vol     = SPL;
capture SVR         = SVR;
capture SVR35       = SVR35_flag;
capture FIBRO_out   = FIBRO;
capture ALLELE_out  = ALLELE;
capture CHR         = CHR;
capture MPN_SAF     = MPN_SAF_TSS;
capture THROMB_RISK = THROMB_RISK;
capture MF_risk     = MF_risk_annual;
capture pSTAT5      = pSTAT5_inhib;
capture Cp_RUX      = Cp_RUX_ngmL;
capture Cp_HYU_uM   = CP_HYU_uM;
capture Cp_IFN      = CP_IFN;

$CAPTURE Hct Hgb EPO_serum PLT_out WBC_out SPL_vol SVR SVR35
$CAPTURE FIBRO_out ALLELE_out CHR MPN_SAF THROMB_RISK MF_risk
$CAPTURE pSTAT5 Cp_RUX Cp_HYU_uM Cp_IFN
'

## ---------------------------------------------------------------------------
## Compile Model
## ---------------------------------------------------------------------------
mod <- mcode("polycythemia_vera", pv_code)

## ---------------------------------------------------------------------------
## Helper: Dosing event builders
## ---------------------------------------------------------------------------

# Ruxolitinib 10 mg BID for 24 weeks
rux_events <- function(dose_mg = 10, duration_weeks = 24) {
  ev(cmt = "DEPOT_RUX",
     amt = dose_mg,
     ii  = 12,          # every 12 hours
     addl = duration_weeks * 14 - 1,
     time = 0)
}

# Hydroxyurea 500 mg QD
hyu_events <- function(dose_mg = 500, duration_weeks = 24) {
  ev(cmt = "CENT_HYU",
     amt = dose_mg,
     ii  = 24,
     addl = duration_weeks * 7 - 1,
     time = 0)
}

# PEG-IFN-α2a 45 μg SC QW → escalating
ifn_events <- function(start_dose = 45, duration_weeks = 48) {
  ev(cmt = "SC_IFN",
     amt = start_dose,
     ii  = 168,          # every 168 hours = 7 days
     addl = duration_weeks - 1,
     time = 0)
}

# Phlebotomy event (single unit)
phlebot_event <- function(time_hr = 0) {
  ev(cmt = 1, amt = 0, time = time_hr, PHLEBOT_FLAG = 1)
}

## ---------------------------------------------------------------------------
## SCENARIO 1: Untreated disease natural history
## ---------------------------------------------------------------------------
cat("\n=== Scenario 1: Untreated PV — Natural History (2 years) ===\n")

sim_untreated <- mod %>%
  param(USE_RUX = 0, USE_HYU = 0, USE_IFN = 0, USE_ASP = 0) %>%
  mrgsim(end = 730, delta = 1) %>%
  as.data.frame()

cat(sprintf("  Baseline Hct: %.1f%%, Week 52 Hct: %.1f%%, Week 104 Hct: %.1f%%\n",
    sim_untreated$Hct[1],
    sim_untreated$Hct[which.min(abs(sim_untreated$time - 365))],
    sim_untreated$Hct[nrow(sim_untreated)]))

## ---------------------------------------------------------------------------
## SCENARIO 2: Phlebotomy + Aspirin (Low-risk PV)
## ---------------------------------------------------------------------------
cat("\n=== Scenario 2: Phlebotomy + Aspirin (Standard care, low-risk) ===\n")

# Monthly phlebotomy for Hct control
phlebot_times <- seq(0, 365*24, by = 30*24)  # every 30 days in hours

sim_phlebot <- mod %>%
  param(USE_RUX = 0, USE_HYU = 0, USE_IFN = 0, USE_ASP = 1,
        PHLEBOT_FLAG = 1) %>%
  mrgsim(end = 365, delta = 1) %>%
  as.data.frame()

cat(sprintf("  Hct at 12 months: %.1f%%\n",
    sim_phlebot$Hct[which.min(abs(sim_phlebot$time - 365))]))

## ---------------------------------------------------------------------------
## SCENARIO 3: Hydroxyurea (High-risk PV, standard)
## ---------------------------------------------------------------------------
cat("\n=== Scenario 3: Hydroxyurea 500mg/d (High-risk PV) ===\n")

hyu_ev <- hyu_events(dose_mg = 500, duration_weeks = 48)

sim_hyu <- mod %>%
  param(USE_RUX = 0, USE_HYU = 1, USE_IFN = 0, USE_ASP = 1) %>%
  mrgsim(ev = hyu_ev, end = 336, delta = 1) %>%
  as.data.frame()

wk24_hyu <- sim_hyu[which.min(abs(sim_hyu$time - 168)), ]
cat(sprintf("  Week 24 — Hct: %.1f%%, PLT: %.0f, WBC: %.1f, CHR: %s\n",
    wk24_hyu$Hct, wk24_hyu$PLT_out, wk24_hyu$WBC_out,
    ifelse(wk24_hyu$CHR == 1, "YES", "NO")))

## ---------------------------------------------------------------------------
## SCENARIO 4: Ruxolitinib 10 mg BID (HYU-intolerant/resistant, RESPONSE trial)
## ---------------------------------------------------------------------------
cat("\n=== Scenario 4: Ruxolitinib 10mg BID (RESPONSE trial simulation) ===\n")

rux_ev <- rux_events(dose_mg = 10, duration_weeks = 32)

sim_rux <- mod %>%
  param(USE_RUX = 1, USE_HYU = 0, USE_IFN = 0, USE_ASP = 1) %>%
  mrgsim(ev = rux_ev, end = 224, delta = 1) %>%
  as.data.frame()

wk32_rux <- sim_rux[which.min(abs(sim_rux$time - 224)), ]
cat(sprintf("  Week 32 — Hct: %.1f%%, PLT: %.0f, WBC: %.1f\n",
    wk32_rux$Hct, wk32_rux$PLT_out, wk32_rux$WBC_out))
cat(sprintf("  SVR35: %s (SVR=%.1f%%), Allele burden: %.1f%%\n",
    ifelse(wk32_rux$SVR35 == 1, "Achieved", "Not achieved"),
    wk32_rux$SVR, wk32_rux$ALLELE_out))
cat(sprintf("  pSTAT5 inhibition: %.1f%%\n", wk32_rux$pSTAT5 * 100))

## ---------------------------------------------------------------------------
## SCENARIO 5: PEG-IFN-α2a (PROUD-PV / CONTINUUM trial simulation)
## ---------------------------------------------------------------------------
cat("\n=== Scenario 5: PEG-IFN-α2a 45μg/wk SC (PROUD-PV simulation) ===\n")

ifn_ev <- ifn_events(start_dose = 45, duration_weeks = 72)

sim_ifn <- mod %>%
  param(USE_RUX = 0, USE_HYU = 0, USE_IFN = 1, USE_ASP = 1) %>%
  mrgsim(ev = ifn_ev, end = 504, delta = 1) %>%
  as.data.frame()

wk72_ifn <- sim_ifn[which.min(abs(sim_ifn$time - 504)), ]
cat(sprintf("  Week 72 — Hct: %.1f%%, Allele burden: %.1f%%\n",
    wk72_ifn$Hct, wk72_ifn$ALLELE_out))
cat(sprintf("  CHR: %s, BM fibrosis: %.2f\n",
    ifelse(wk72_ifn$CHR == 1, "YES", "NO"),
    wk72_ifn$FIBRO_out))

## ---------------------------------------------------------------------------
## PLOT: All scenarios — key endpoints
## ---------------------------------------------------------------------------
make_scenario_df <- function(sim, label) {
  sim %>%
    mutate(Scenario = label) %>%
    select(time, Scenario, Hct, PLT_out, WBC_out, SPL_vol, ALLELE_out,
           CHR, MPN_SAF, THROMB_RISK, FIBRO_out, SVR)
}

# Combine for comparison (use 168-day window common to all)
df_nat <- sim_untreated %>% filter(time <= 168) %>% make_scenario_df("Untreated")
df_hyu <- sim_hyu %>% filter(time <= 168) %>% make_scenario_df("Hydroxyurea 500mg/d")
df_rux <- sim_rux %>% filter(time <= 168) %>% make_scenario_df("Ruxolitinib 10mg BID")
df_ifn <- sim_ifn %>% filter(time <= 168) %>% make_scenario_df("PEG-IFN-α2a 45μg/wk")

df_all <- bind_rows(df_nat, df_hyu, df_rux, df_ifn)
df_all$time_wk <- df_all$time / 7

# Color palette
cols <- c("Untreated"="#E74C3C",
          "Hydroxyurea 500mg/d"="#3498DB",
          "Ruxolitinib 10mg BID"="#2ECC71",
          "PEG-IFN-α2a 45μg/wk"="#9B59B6")

p1 <- ggplot(df_all, aes(x=time_wk, y=Hct, color=Scenario)) +
  geom_line(linewidth=1.2) +
  geom_hline(yintercept=45, linetype="dashed", color="black") +
  annotate("text", x=22, y=46, label="Target Hct <45%", size=3) +
  scale_color_manual(values=cols) +
  labs(title="Hematocrit (%)", x="Time (weeks)", y="Hct (%)") +
  theme_bw(base_size=11) + theme(legend.position="bottom")

p2 <- ggplot(df_all, aes(x=time_wk, y=PLT_out, color=Scenario)) +
  geom_line(linewidth=1.2) +
  geom_hline(yintercept=400, linetype="dashed", color="black") +
  scale_color_manual(values=cols) +
  labs(title="Platelet Count (×10⁹/L)", x="Time (weeks)", y="Platelets (×10⁹/L)") +
  theme_bw(base_size=11) + theme(legend.position="none")

p3 <- ggplot(df_all, aes(x=time_wk, y=SPL_vol, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=cols) +
  labs(title="Spleen Volume (mL)", x="Time (weeks)", y="Volume (mL)") +
  theme_bw(base_size=11) + theme(legend.position="none")

p4 <- ggplot(df_all, aes(x=time_wk, y=ALLELE_out, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=cols) +
  labs(title="JAK2 V617F Allele Burden (%)", x="Time (weeks)", y="Allele Burden (%)") +
  theme_bw(base_size=11) + theme(legend.position="none")

p5 <- ggplot(df_all, aes(x=time_wk, y=MPN_SAF, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=cols) +
  labs(title="Symptom Burden (MPN-SAF TSS)", x="Time (weeks)", y="MPN-SAF TSS") +
  theme_bw(base_size=11) + theme(legend.position="none")

p6 <- ggplot(df_all, aes(x=time_wk, y=THROMB_RISK, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=cols) +
  labs(title="Annual Thrombosis Risk", x="Time (weeks)", y="Annual Risk (fraction)") +
  theme_bw(base_size=11) + theme(legend.position="none")

# Combine plots
combined_plot <- (p1 + p2) / (p3 + p4) / (p5 + p6) +
  plot_annotation(
    title = "Polycythemia Vera QSP Model — Treatment Scenarios",
    subtitle = "Ruxolitinib, Hydroxyurea, PEG-IFN-α2a, Untreated (24-week horizon)",
    theme = theme(plot.title = element_text(face="bold", size=14))
  )

print(combined_plot)

## ---------------------------------------------------------------------------
## SCENARIO 6: Ruxolitinib PK dose-response (5, 10, 15, 20 mg BID)
## ---------------------------------------------------------------------------
cat("\n=== Scenario 6: Ruxolitinib Dose-Response (5/10/15/20 mg BID) ===\n")

doses <- c(5, 10, 15, 20)
dose_results <- lapply(doses, function(d) {
  ev_d <- rux_events(dose_mg = d, duration_weeks = 24)
  sim <- mod %>%
    param(USE_RUX = 1, USE_HYU = 0, USE_IFN = 0, USE_ASP = 1) %>%
    mrgsim(ev = ev_d, end = 168, delta = 1) %>%
    as.data.frame() %>%
    mutate(Dose = paste0(d, " mg BID"),
           time_wk = time / 7)
})
df_dose <- bind_rows(dose_results)

p_dose <- ggplot(df_dose, aes(x=time_wk, y=Hct, color=Dose)) +
  geom_line(linewidth=1.2) +
  geom_hline(yintercept=45, linetype="dashed") +
  scale_color_brewer(palette="RdYlGn") +
  labs(title="Ruxolitinib Dose-Response: Hematocrit",
       x="Time (weeks)", y="Hct (%)") +
  theme_bw(base_size=12)
print(p_dose)

## ---------------------------------------------------------------------------
## Summary Table
## ---------------------------------------------------------------------------
cat("\n=== Summary: 24-Week Clinical Endpoints ===\n")
summary_tab <- data.frame(
  Scenario = c("Untreated", "HYU 500mg/d", "RUX 10mg BID", "PEG-IFN 45μg/wk"),
  Hct_wk24 = sapply(list(sim_untreated, sim_hyu, sim_rux, sim_ifn),
                    function(s) round(s$Hct[which.min(abs(s$time - 168))], 1)),
  PLT_wk24 = sapply(list(sim_untreated, sim_hyu, sim_rux, sim_ifn),
                    function(s) round(s$PLT_out[which.min(abs(s$time - 168))], 0)),
  Allele_wk24 = sapply(list(sim_untreated, sim_hyu, sim_rux, sim_ifn),
                       function(s) round(s$ALLELE_out[which.min(abs(s$time - 168))], 1)),
  SVR_wk24 = sapply(list(sim_untreated, sim_hyu, sim_rux, sim_ifn),
                    function(s) round(s$SVR[which.min(abs(s$time - 168))], 1)),
  CHR_wk24 = sapply(list(sim_untreated, sim_hyu, sim_rux, sim_ifn),
                    function(s) ifelse(s$CHR[which.min(abs(s$time - 168))] == 1, "Yes", "No"))
)
print(summary_tab, row.names=FALSE)
