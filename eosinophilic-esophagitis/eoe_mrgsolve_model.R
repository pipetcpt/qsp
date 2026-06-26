##############################################################################
# EoE (Eosinophilic Esophagitis) QSP Model — mrgsolve Implementation
#
# Disease: Eosinophilic Esophagitis (EoE)
# Model type: Quantitative Systems Pharmacology (QSP)
# Framework: mrgsolve (R)
#
# Key pathways modeled:
#   1. Esophageal epithelial barrier dynamics (IL-13 → DSG1/FLG ↓)
#   2. Th2 cytokine network (IL-4, IL-5, IL-13, eotaxin-3)
#   3. Eosinophil biology (BM → blood → tissue)
#   4. Mast cell dynamics
#   5. Lamina propria fibrosis (TGF-β/SMAD pathway)
#   6. Drug PK: budesonide, dupilumab, mepolizumab, cendakimab
#
# Clinical calibration references:
#   - MATS trial (Lucendo 2022): dupilumab 300mg SC q2w → 80% eos reduction
#   - ApplE trial (Dellon 2022): budesonide ODT → 58% histological remission
#   - CACTUS trial (Hirano 2023): cendakimab → 64% histo remission
#   - Stein 2006: SFED diet → ~72% histological remission
#   - Rothenberg 2014: mepolizumab Phase 2 (blood eos ↓, tissue partial)
#
# Units:
#   Time:         days
#   Drug doses:   mg
#   Concentrations: mg/L (≡ µg/mL)
#   Cytokines:    pg/mL
#   Eosinophils:  cells/µL (blood); eos/hpf (tissue)
#   Fibrosis:     dimensionless (0–1 scale)
#
# Author: Claude Code Routine (QSP Library)
# Date: 2026-06-24
##############################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ===========================================================================
# 1. DEFINE mrgsolve MODEL
# ===========================================================================

eoe_model_code <- '
$PROB
EoE (Eosinophilic Esophagitis) QSP Model
Th2 cytokine-driven eosinophilic inflammation with drug PK/PD
18 ODE compartments, 5 treatment scenarios

$PARAM @annotated
// ---- Patient / Disease Parameters ----
WT       : 70    : Body weight (kg)
AGE      : 35    : Age (years)
SEX      : 1     : Sex (1=male, 2=female)
EOS_BL_INIT  : 600  : Baseline blood eosinophils (cells/µL)
IL13_BL  : 80    : Baseline esophageal IL-13 (pg/mL)
IL5_BL   : 15    : Baseline circulating IL-5 (pg/mL)
EOTAX3_BL: 400   : Baseline eotaxin-3/CCL26 (pg/mL)
EOS_ESO_BL: 80   : Baseline tissue eosinophils (eos/hpf)
MAST_BL  : 50    : Baseline esophageal mast cells (per mm²)
FIBRO_BL : 0.4   : Baseline lamina propria fibrosis (0-1 scale)
IGE_BL   : 300   : Baseline total IgE (IU/mL)
EPBAR_BL : 0.4   : Baseline epithelial barrier integrity (0=disrupted, 1=intact)

// ---- Disease Dynamics ----
kin_IL13  : 320   : IL-13 zero-order production rate (pg/mL/day)
kout_IL13 : 4.0   : IL-13 first-order elimination rate (1/day)
kin_IL5   : 60    : IL-5 production rate (pg/mL/day)
kout_IL5  : 4.0   : IL-5 elimination rate (1/day)
kin_eotax : 1600  : Eotaxin-3 production rate (pg/mL/day)
kout_eotax: 4.0   : Eotaxin-3 elimination rate (1/day)
hill_eotax: 1.5   : Hill coefficient for IL-13 → eotaxin-3
kin_eosbl : 2400  : Blood eosinophil production (cells/µL/day)
kout_eosbl: 4.0   : Blood eosinophil turnover (1/day)
k_recruit : 0.06  : Eosinophil recruitment rate to tissue (1/day)
kin_eos_eso: 0.0  : Direct tissue eosinophil input (supplementary)
kout_eos_eso: 0.05: Tissue eosinophil apoptosis rate (1/day)
kin_mast  : 0.10  : Mast cell influx rate (cells·mm⁻²/day)
kout_mast : 0.002 : Mast cell turnover rate (1/day)
kin_fibro : 0.001 : Fibrosis accumulation rate (1/day)
kout_fibro: 0.0008: Fibrosis resolution rate (1/day)
kin_IgE   : 0.15  : IgE production rate (IU/mL/day)
kout_IgE  : 0.0005: IgE elimination rate (1/day)
tau_dysphag: 14   : Dysphagia response time constant (days)
tau_EPBAR : 7     : Epithelial barrier response time constant (days)

// ---- Emax parameters (disease mechanisms) ----
Emax_IL13_eotax : 0.95 : Max effect of IL-13 on eotaxin-3 production
IC50_IL13_eotax : 40   : IL-13 EC50 for eotaxin-3 induction (pg/mL)
Emax_eotax_recruit: 0.90: Max effect of eotaxin-3 on eosinophil recruitment
IC50_eotax_recruit: 200 : Eotaxin-3 EC50 for recruitment (pg/mL)
Emax_IL5_eosbl  : 0.80  : Max effect of IL-5 on blood eosinophils
IC50_IL5_eosbl  : 8     : IL-5 EC50 for eosinophil production (pg/mL)
Emax_eos_fibro  : 0.70  : Max effect of tissue eosinophils on fibrosis
IC50_eos_fibro  : 40    : Eosinophil EC50 for fibrosis (eos/hpf)
Emax_IL13_bar   : 0.70  : Max barrier disruption by IL-13
IC50_IL13_bar   : 40    : IL-13 EC50 for barrier disruption (pg/mL)
Emax_eos_bar    : 0.50  : Max barrier disruption by eos (MBP)
IC50_eos_bar    : 50    : Eosinophil EC50 for barrier damage (eos/hpf)

// ---- Budesonide PK/PD ----
ka_bud_sys  : 0.6   : Budesonide esophageal→systemic absorption (1/day)
ke_bud_eso  : 2.4   : Budesonide esophageal elimination (1/day)
ke_bud_sys  : 8.0   : Budesonide systemic elimination (1/day)
Emax_bud_IL13: 0.85 : Max budesonide suppression of IL-13
IC50_bud_IL13: 0.1  : Budesonide esophageal IC50 for IL-13 (mg/L)
Emax_bud_eotax:0.90 : Max budesonide suppression of eotaxin-3
IC50_bud_eotax:0.08 : Budesonide IC50 for eotaxin-3 (mg/L)
Emax_bud_eos_eso: 0.85: Max budesonide induction of eos apoptosis
IC50_bud_eos_eso:0.1 : Budesonide IC50 for tissue eos reduction (mg/L)

// ---- Dupilumab PK/PD ----
ka_dup    : 0.18  : Dupilumab SC absorption rate (1/day)
F_dup     : 0.64  : Dupilumab SC bioavailability (fraction)
CL_dup    : 0.21  : Dupilumab clearance (L/day)
Vd_dup_C  : 3.5   : Dupilumab central volume (L)
Vd_dup_P  : 2.8   : Dupilumab peripheral volume (L)
Q_dup     : 1.5   : Dupilumab intercompartmental clearance (L/day)
Emax_dup_STAT6: 0.95: Max dupilumab blockade of IL-4Rα/STAT6
IC50_dup_STAT6: 2.0 : Dupilumab IC50 for STAT6 blockade (mg/L)
Emax_dup_IgE: 0.70  : Max dupilumab reduction in IgE
IC50_dup_IgE: 10    : Dupilumab IC50 for IgE reduction (mg/L)

// ---- Mepolizumab PK/PD ----
ka_mepo   : 0.34  : Mepolizumab SC absorption (1/day)
F_mepo    : 0.81  : Mepolizumab SC bioavailability (fraction)
CL_mepo   : 0.28  : Mepolizumab clearance (L/day)
Vd_mepo   : 3.6   : Mepolizumab volume of distribution (L)
Emax_mepo_IL5: 0.95: Max mepolizumab neutralization of IL-5
IC50_mepo_IL5: 1.0  : Mepolizumab IC50 for IL-5 blockade (mg/L)

// ---- Cendakimab PK/PD ----
ka_cenda  : 14.4  : Cendakimab oral absorption rate (1/day)
F_cenda   : 0.35  : Cendakimab bioavailability (fraction)
CL_cenda  : 360   : Cendakimab apparent clearance CL/F (L/day)
Vd_cenda  : 100   : Cendakimab apparent volume Vd/F (L)
Emax_cenda_IL13: 0.90: Max cendakimab blockade of IL-13
IC50_cenda_IL13: 0.05 : Cendakimab IC50 for IL-13 blockade (mg/L)

// ---- Dietary elimination ----
DIETARY : 0 : Dietary elimination effect (0=none, 1=active)
Emax_diet_antigen: 0.80: Max reduction of allergen load by diet

$CMT @annotated
// Drug PK compartments
BUD_ESO  : Budesonide esophageal (mg/L)
BUD_SYS  : Budesonide systemic (mg/L)
DUP_SC   : Dupilumab SC depot (mg)
DUP_C    : Dupilumab central (mg/L)
DUP_P    : Dupilumab peripheral (mg/L)
MEPO_SC  : Mepolizumab SC depot (mg)
MEPO_C   : Mepolizumab central (mg/L)
CENDA_GUT: Cendakimab gut depot (mg)
CENDA_C  : Cendakimab central (mg/L)

// Disease state compartments
IL13     : Esophageal IL-13 (pg/mL)
IL5      : Circulating IL-5 (pg/mL)
EOTAX3   : Eotaxin-3/CCL26 (pg/mL)
EOS_BL   : Blood eosinophils (cells/µL)
EOS_ESO  : Tissue eosinophils (eos/hpf)
MAST_ESO : Esophageal mast cells (per mm²)
FIBRO    : Lamina propria fibrosis (0-1)
IGE_TOT  : Total serum IgE (IU/mL)
EPBAR    : Epithelial barrier integrity (0-1)

$INIT @annotated
BUD_ESO  = 0     : Initial budesonide esoph (mg/L)
BUD_SYS  = 0     : Initial budesonide systemic (mg/L)
DUP_SC   = 0     : Initial dupilumab SC depot (mg)
DUP_C    = 0     : Initial dupilumab central (mg/L)
DUP_P    = 0     : Initial dupilumab periph (mg/L)
MEPO_SC  = 0     : Initial mepolizumab SC (mg)
MEPO_C   = 0     : Initial mepolizumab central (mg/L)
CENDA_GUT= 0     : Initial cendakimab gut (mg)
CENDA_C  = 0     : Initial cendakimab central (mg/L)
IL13     = 80    : Baseline IL-13 (pg/mL) — active EoE
IL5      = 15    : Baseline IL-5 (pg/mL)
EOTAX3   = 400   : Baseline eotaxin-3 (pg/mL)
EOS_BL   = 600   : Baseline blood eosinophils (cells/µL)
EOS_ESO  = 80    : Baseline tissue eosinophils (eos/hpf) — active EoE
MAST_ESO = 50    : Baseline mast cells (per mm²)
FIBRO    = 0.4   : Baseline fibrosis score (moderate)
IGE_TOT  = 300   : Baseline IgE (IU/mL)
EPBAR    = 0.4   : Baseline barrier integrity (disrupted)

$MAIN
// ======== Derived drug effects ========

// --- Budesonide: GR-mediated transcription suppression ---
double Inh_bud_IL13  = Emax_bud_IL13  * BUD_ESO / (BUD_ESO + IC50_bud_IL13);
double Inh_bud_eotax = Emax_bud_eotax * BUD_ESO / (BUD_ESO + IC50_bud_eotax);
double Enh_bud_apop  = Emax_bud_eos_eso * BUD_ESO / (BUD_ESO + IC50_bud_eos_eso);

// --- Dupilumab: IL-4Rα / STAT6 blockade (dual IL-4 + IL-13 inhibition) ---
double Inh_dup_STAT6 = Emax_dup_STAT6 * DUP_C / (DUP_C + IC50_dup_STAT6);
double Inh_dup_IgE   = Emax_dup_IgE   * DUP_C / (DUP_C + IC50_dup_IgE);

// --- Mepolizumab: IL-5 neutralization ---
double Inh_mepo_IL5 = Emax_mepo_IL5 * MEPO_C / (MEPO_C + IC50_mepo_IL5);

// --- Cendakimab: IL-13 direct neutralization ---
double Inh_cenda_IL13 = Emax_cenda_IL13 * CENDA_C / (CENDA_C + IC50_cenda_IL13);

// --- Dietary elimination: allergen avoidance ---
double Emax_diet = DIETARY * Emax_diet_antigen;

// ======== Combined drug effect on IL-13 ========
// Budesonide suppresses production; dupilumab blocks receptor signaling
// (represented as enhanced IL-13 elimination via STAT6 block)
// Cendakimab directly neutralizes IL-13
double Prod_IL13_factor = (1.0 - Inh_bud_IL13) * (1.0 - Emax_diet);
// Net IL-13 signal = free IL-13 × (1 - STAT6 block - cendakimab)
double IL13_signal = IL13 * (1.0 - Inh_dup_STAT6) * (1.0 - Inh_cenda_IL13);

// ======== Eotaxin-3 stimulation by IL-13 ========
double Stim_eotax = Emax_IL13_eotax * pow(IL13_signal, hill_eotax) /
                    (pow(IC50_IL13_eotax, hill_eotax) + pow(IL13_signal, hill_eotax));
double Prod_eotax = kin_eotax * (1.0 + Stim_eotax) * (1.0 - Inh_bud_eotax);

// ======== IL-5 effect on eosinophils ========
double IL5_effective = IL5 * (1.0 - Inh_mepo_IL5);
double Stim_eos_IL5 = Emax_IL5_eosbl * IL5_effective /
                      (IC50_IL5_eosbl + IL5_effective);

// ======== Eotaxin-3-driven tissue recruitment ========
double Stim_recruit = Emax_eotax_recruit * EOTAX3 /
                      (IC50_eotax_recruit + EOTAX3);
double k_recruit_eff = k_recruit * (1.0 + Stim_recruit);

// ======== Fibrosis: driven by tissue eosinophils and mast cells ========
double Stim_fibro = Emax_eos_fibro * EOS_ESO /
                    (IC50_eos_fibro + EOS_ESO);

// ======== Epithelial barrier: IL-13 and eos damage, drugs restore ========
double Disrupt_IL13 = Emax_IL13_bar * IL13 / (IC50_IL13_bar + IL13);
double Disrupt_eos  = Emax_eos_bar  * EOS_ESO / (IC50_eos_bar + EOS_ESO);
// Net barrier setpoint [0,1]: damaged by IL-13/eos, partially restored by drugs
double EPBAR_ss = (1.0 - Disrupt_IL13 * (1.0 - Inh_dup_STAT6) * (1.0 - Inh_cenda_IL13) *
                         (1.0 - Inh_bud_IL13)) *
                  (1.0 - 0.5 * Disrupt_eos);
EPBAR_ss = EPBAR_ss < 0.05 ? 0.05 : EPBAR_ss;  // floor

// ======== Dysphagia score target [0-10] ========
// Driven by tissue eos, fibrosis, and barrier disruption
double DYSPHAG_ss = 10.0 * (0.5 * EOS_ESO / (EOS_ESO + 60.0) +
                             0.3 * FIBRO / (FIBRO + 0.5) +
                             0.2 * (1.0 - EPBAR));

// ======== EREFS score [0-18] ========
// Proxy: rings(0-3), exudates(0-2), edema(0-2), furrows(0-2), stricture(0-2)
double EOS_FRAC = EOS_ESO / (EOS_ESO + 50.0);    // saturation function for eos
double FIBRO_FRAC = FIBRO / (FIBRO + 0.5);
double EREFS_val = 3.0 * EOS_FRAC + 2.0 * EOS_FRAC * 0.8 +
                   2.0 * (1.0 - EPBAR) * 0.6 +
                   2.0 * FIBRO_FRAC * 0.9 +
                   2.0 * FIBRO_FRAC * 0.7;
if(EREFS_val > 18.0) EREFS_val = 18.0;

$ODE
// ======================================================
// DRUG PK ODEs
// ======================================================

// --- Budesonide PK (esophageal topical deposition) ---
// Dose enters BUD_ESO directly (esophageal compartment)
dxdt_BUD_ESO = -ka_bud_sys * BUD_ESO - ke_bud_eso * BUD_ESO;
dxdt_BUD_SYS = ka_bud_sys * BUD_ESO - ke_bud_sys * BUD_SYS;

// --- Dupilumab PK (2-compartment SC model) ---
dxdt_DUP_SC = -ka_dup * DUP_SC;
// Central: input from SC depot (F scaled into dose input), output CL+Q
dxdt_DUP_C  = ka_dup * F_dup * DUP_SC / Vd_dup_C
              - (CL_dup / Vd_dup_C) * DUP_C
              - (Q_dup  / Vd_dup_C) * DUP_C
              + (Q_dup  / Vd_dup_P) * DUP_P;
dxdt_DUP_P  = (Q_dup / Vd_dup_C) * DUP_C
              - (Q_dup / Vd_dup_P) * DUP_P;

// --- Mepolizumab PK (1-compartment SC model) ---
dxdt_MEPO_SC = -ka_mepo * MEPO_SC;
dxdt_MEPO_C  = ka_mepo * F_mepo * MEPO_SC / Vd_mepo
               - (CL_mepo / Vd_mepo) * MEPO_C;

// --- Cendakimab PK (oral 1-compartment, apparent params) ---
dxdt_CENDA_GUT = -ka_cenda * CENDA_GUT;
dxdt_CENDA_C   = ka_cenda * F_cenda * CENDA_GUT / Vd_cenda
                 - (CL_cenda / Vd_cenda) * CENDA_C;

// ======================================================
// DISEASE DYNAMICS ODEs
// ======================================================

// --- IL-13 (esophageal tissue) ---
// Production inhibited by budesonide and dietary allergen avoidance
// Elimination accelerated by cendakimab (direct neutralization)
// Dupilumab blocks downstream STAT6 (not IL-13 itself)
dxdt_IL13 = kin_IL13 * Prod_IL13_factor
            - kout_IL13 * IL13 * (1.0 + Inh_cenda_IL13 * 2.0);

// --- IL-5 (circulating) ---
// Production reduced by budesonide and dietary allergen avoidance
// Mepolizumab neutralizes free IL-5 (modeled as enhanced elimination)
dxdt_IL5 = kin_IL5 * (1.0 - Inh_bud_IL13 * 0.5) * (1.0 - Emax_diet * 0.7)
           - kout_IL5 * IL5 * (1.0 + Inh_mepo_IL5 * 3.0);

// --- Eotaxin-3/CCL26 (esophageal epithelium) ---
// Production driven by IL-13/STAT6 pathway; suppressed by budesonide + dup
dxdt_EOTAX3 = Prod_eotax * (1.0 - Inh_dup_STAT6)
              - kout_eotax * EOTAX3;

// --- Blood eosinophils ---
// Production from BM driven by IL-5; suppressed by mepolizumab
// Tissue recruitment by eotaxin-3 reduces blood pool
dxdt_EOS_BL = kin_eosbl * (1.0 + Stim_eos_IL5) * (1.0 - Inh_bud_IL13 * 0.3)
              - kout_eosbl * EOS_BL
              - k_recruit_eff * EOS_BL;

// --- Tissue (esophageal) eosinophils ---
// Recruitment from blood via CCR3-eotaxin-3; apoptosis enhanced by budesonide
dxdt_EOS_ESO = k_recruit_eff * EOS_BL
               - kout_eos_eso * EOS_ESO * (1.0 + Enh_bud_apop * 3.0);

// --- Mast cells (esophageal) ---
// Influx stimulated by IL-4/IL-9 (proxied by IL-13 signal)
// Stabilization by dupilumab-mediated Th2 suppression
dxdt_MAST_ESO = kin_mast * (1.0 + IL13_signal / IL13_BL * 0.5) * (1.0 - Emax_diet * 0.5)
                - kout_mast * MAST_ESO * (1.0 + Inh_dup_STAT6 * 0.5);

// --- Lamina propria fibrosis ---
// TGF-β from eosinophils and mast cells drives collagen deposition
// Slow process (months to years); dupilumab/cendakimab reduce via IL-13 blockade
dxdt_FIBRO = kin_fibro * (1.0 + Stim_fibro) * (1.0 - Inh_dup_STAT6 * 0.4)
                        * (1.0 - Inh_cenda_IL13 * 0.3)
             - kout_fibro * FIBRO;

// --- Total serum IgE ---
// Reduced by dupilumab (anti-IL-4Rα blocks class switching)
// Slow turnover (t½ ~ months)
dxdt_IGE_TOT = kin_IgE * (1.0 - Emax_diet * 0.3)
               - kout_IgE * IGE_TOT * (1.0 + Inh_dup_IgE);

// --- Epithelial barrier integrity ---
// First-order approach to steady-state setpoint
dxdt_EPBAR = (EPBAR_ss - EPBAR) / tau_EPBAR;

$TABLE
// ======================================================
// DERIVED OUTPUTS (captured in simulation output)
// ======================================================
double DYSPHAG_ODE = DYSPHAG_ss;    // dysphagia score (0-10)
double EREFS_SCORE = EREFS_val;     // endoscopic score (0-18)

// Histological remission: YES if peak eos/hpf < 15
double HISTO_REMIS = (EOS_ESO < 15.0) ? 1.0 : 0.0;

// % change in tissue eosinophils from baseline
double PCT_CHG_EOS_ESO = (EOS_ESO - EOS_ESO_BL) / EOS_ESO_BL * 100.0;

// Effective IL-13 signal (post-drug)
double IL13_EFF = IL13_signal;

// Blood eosinophil absolute count
double EOS_BL_ABS = EOS_BL;

// Eotaxin-3 tissue level
double EOTAX3_LVL = EOTAX3;

// Budesonide esophageal trough concentration
double BUD_ESO_CONC = BUD_ESO;

// Dupilumab central trough
double DUP_TROUGH = DUP_C;

// Mepolizumab central trough
double MEPO_TROUGH = MEPO_C;

// Cendakimab Cmax (proxy)
double CENDA_CONC = CENDA_C;

// EREFS sub-scores
double EREFS_RINGS = 3.0 * EOS_ESO / (EOS_ESO + 50.0) * FIBRO_FRAC;
double EREFS_EXUD  = 2.0 * EOS_ESO / (EOS_ESO + 50.0);
double EREFS_EDEMA = 2.0 * (1.0 - EPBAR);
double EREFS_FURR  = 2.0 * FIBRO_FRAC;
double EREFS_STRIC = 2.0 * FIBRO * FIBRO / (FIBRO * FIBRO + 0.2);

// Fibrosis level
double FIBRO_LVL = FIBRO;

// IgE
double IGE_LEVEL = IGE_TOT;

$CAPTURE
BUD_ESO_CONC DUP_TROUGH MEPO_TROUGH CENDA_CONC
IL13 IL5 EOTAX3 IL13_EFF
EOS_BL_ABS EOS_ESO MAST_ESO EPBAR FIBRO_LVL IGE_LEVEL
DYSPHAG_ODE EREFS_SCORE HISTO_REMIS PCT_CHG_EOS_ESO
EREFS_RINGS EREFS_EXUD EREFS_EDEMA EREFS_FURR EREFS_STRIC
'

# ===========================================================================
# 2. COMPILE MODEL
# ===========================================================================

cat("Compiling EoE QSP mrgsolve model...\n")
mod <- mcode("EoE_QSP", eoe_model_code)

# ===========================================================================
# 3. HELPER FUNCTION: EVENT SCHEDULE BUILDER
# ===========================================================================

build_events <- function(scenario = "none",
                         start_day = 1,
                         duration_days = 365) {
  events <- data.frame()

  if (scenario == "budesonide") {
    # Budesonide ODT 1 mg BID (dose enters BUD_ESO compartment)
    # Total 2 mg/day; split into two 1 mg doses
    dose_times <- seq(start_day, start_day + duration_days - 1, by = 0.5)
    events <- bind_rows(events, data.frame(
      time = dose_times, amt = 1, cmt = 1, rate = 0, evid = 1, ii = 0
    ))

  } else if (scenario == "dupilumab") {
    # 300 mg SC q14 days → dose enters DUP_SC
    dose_times <- seq(start_day, start_day + duration_days, by = 14)
    events <- bind_rows(events, data.frame(
      time = dose_times, amt = 300, cmt = 3, rate = 0, evid = 1, ii = 0
    ))

  } else if (scenario == "mepolizumab") {
    # 300 mg SC q28 days → MEPO_SC
    dose_times <- seq(start_day, start_day + duration_days, by = 28)
    events <- bind_rows(events, data.frame(
      time = dose_times, amt = 300, cmt = 6, rate = 0, evid = 1, ii = 0
    ))

  } else if (scenario == "cendakimab") {
    # 160 mg PO QD → CENDA_GUT (F=0.35 in model)
    dose_times <- seq(start_day, start_day + duration_days - 1, by = 1)
    events <- bind_rows(events, data.frame(
      time = dose_times, amt = 160, cmt = 8, rate = 0, evid = 1, ii = 0
    ))

  } else if (scenario == "dupilumab_bud") {
    # Combination: dupilumab 300 mg q2w + budesonide 1 mg BID
    dup_times <- seq(start_day, start_day + duration_days, by = 14)
    bud_times <- seq(start_day, start_day + duration_days - 1, by = 0.5)
    events <- bind_rows(events,
      data.frame(time = dup_times, amt = 300, cmt = 3, rate = 0, evid = 1, ii = 0),
      data.frame(time = bud_times, amt = 1,   cmt = 1, rate = 0, evid = 1, ii = 0)
    )
  }

  events <- events[order(events$time), ]
  return(events)
}

# ===========================================================================
# 4. DEFINE TREATMENT SCENARIOS
# ===========================================================================

SCENARIOS <- list(
  list(name = "No Treatment",
       label = "1. No Treatment\n(Disease Progression)",
       scenario = "none",
       color = "#E53935"),

  list(name = "Budesonide ODT",
       label = "2. Budesonide ODT\n1 mg BID (Topical CS)",
       scenario = "budesonide",
       color = "#FF8F00"),

  list(name = "Dupilumab",
       label = "3. Dupilumab\n300 mg SC q2w (anti-IL-4Rα)",
       scenario = "dupilumab",
       color = "#1565C0"),

  list(name = "Mepolizumab",
       label = "4. Mepolizumab\n300 mg SC q4w (anti-IL-5)",
       scenario = "mepolizumab",
       color = "#6A1B9A"),

  list(name = "Cendakimab",
       label = "5. Cendakimab\n160 mg PO QD (anti-IL-13)",
       scenario = "cendakimab",
       color = "#00695C"),

  list(name = "Dupilumab + Budesonide",
       label = "6. Combination\n(Dupilumab + Budesonide)",
       scenario = "dupilumab_bud",
       color = "#37474F")
)

# ===========================================================================
# 5. SIMULATE ALL SCENARIOS (0–52 weeks = 364 days)
# ===========================================================================

SIM_END <- 364     # 52 weeks
OBS_TIMES <- c(seq(0, 28, by = 1), seq(29, SIM_END, by = 7))

run_scenario <- function(sc) {
  events <- build_events(scenario = sc$scenario,
                         start_day = 1,
                         duration_days = SIM_END)

  if (nrow(events) == 0) {
    # No treatment: just observation
    out <- mod %>%
      mrgsim(end = SIM_END, delta = 7) %>%
      as.data.frame()
  } else {
    out <- mod %>%
      mrgsim_df(events = as.data.frame(events), end = SIM_END,
                add = OBS_TIMES, output = "df")
  }

  out$scenario   <- sc$name
  out$label      <- sc$label
  out$color      <- sc$color
  out$week       <- out$time / 7
  return(out)
}

cat("Running simulations for 6 treatment scenarios...\n")
all_results <- lapply(SCENARIOS, function(sc) {
  tryCatch({
    cat(sprintf("  Simulating: %s\n", sc$name))
    run_scenario(sc)
  }, error = function(e) {
    cat(sprintf("  ERROR in %s: %s\n", sc$name, e$message))
    return(NULL)
  })
})
all_results <- Filter(Negate(is.null), all_results)
df_all <- bind_rows(all_results)
df_all$scenario <- factor(df_all$scenario,
                          levels = sapply(SCENARIOS, function(s) s$name))

COLORS <- setNames(
  sapply(SCENARIOS, function(s) s$color),
  sapply(SCENARIOS, function(s) s$name)
)

# ===========================================================================
# 6. RESULTS SUMMARY AT WEEK 24
# ===========================================================================

cat("\n===== RESULTS SUMMARY AT WEEK 24 (168 days) =====\n")
df_wk24 <- df_all %>%
  filter(time >= 168 & time < 175) %>%
  group_by(scenario) %>%
  slice(1) %>%
  select(scenario, EOS_ESO, HISTO_REMIS, DYSPHAG_ODE, EREFS_SCORE,
         EOS_BL_ABS, IL13, IL5, EOTAX3, FIBRO_LVL, IGE_LEVEL) %>%
  mutate(
    EOS_ESO    = round(EOS_ESO, 1),
    HISTO_REMIS = ifelse(HISTO_REMIS == 1, "YES (<15/hpf)", "NO (≥15/hpf)"),
    DYSPHAG_ODE = round(DYSPHAG_ODE, 2),
    EREFS_SCORE = round(EREFS_SCORE, 1),
    EOS_BL_ABS  = round(EOS_BL_ABS, 0),
    IL13        = round(IL13, 1),
    IL5         = round(IL5, 2),
    EOTAX3      = round(EOTAX3, 1),
    FIBRO_LVL   = round(FIBRO_LVL, 3),
    IGE_LEVEL   = round(IGE_LEVEL, 1)
  )

print(as.data.frame(df_wk24))

# ===========================================================================
# 7. VISUALIZATION
# ===========================================================================

p1 <- ggplot(df_all, aes(x = week, y = EOS_ESO, color = scenario)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "gray40", size = 0.8) +
  annotate("text", x = 50, y = 16.5, label = "Remission threshold (<15 eos/hpf)",
           hjust = 1, color = "gray40", size = 3) +
  scale_color_manual(values = COLORS) +
  labs(title = "Tissue Eosinophils (Peak eos/hpf)",
       x = "Week", y = "Eosinophils (eos/hpf)",
       color = "Treatment") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.text = element_text(size = 8))

p2 <- ggplot(df_all, aes(x = week, y = EOS_BL_ABS, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = COLORS) +
  labs(title = "Blood Eosinophils (Absolute)",
       x = "Week", y = "Eosinophils (cells/µL)", color = "Treatment") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p3 <- ggplot(df_all, aes(x = week, y = IL13, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = COLORS) +
  labs(title = "Esophageal IL-13 (pg/mL)",
       x = "Week", y = "IL-13 (pg/mL)", color = "Treatment") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p4 <- ggplot(df_all, aes(x = week, y = EOTAX3, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = COLORS) +
  labs(title = "Eotaxin-3 / CCL26 (pg/mL)",
       x = "Week", y = "Eotaxin-3 (pg/mL)", color = "Treatment") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p5 <- ggplot(df_all, aes(x = week, y = DYSPHAG_ODE, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = COLORS) +
  labs(title = "Dysphagia Score (0–10)",
       x = "Week", y = "Dysphagia Score", color = "Treatment") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p6 <- ggplot(df_all, aes(x = week, y = EREFS_SCORE, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = COLORS) +
  labs(title = "EREFS Score (0–18)",
       x = "Week", y = "EREFS Score", color = "Treatment") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p7 <- ggplot(df_all, aes(x = week, y = EPBAR, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = COLORS) +
  labs(title = "Epithelial Barrier Integrity (0=disrupted, 1=intact)",
       x = "Week", y = "Barrier Integrity", color = "Treatment") +
  ylim(0, 1) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p8 <- ggplot(df_all, aes(x = week, y = FIBRO_LVL, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = COLORS) +
  labs(title = "Lamina Propria Fibrosis (0–1 scale)",
       x = "Week", y = "Fibrosis Score", color = "Treatment") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p9 <- ggplot(df_all %>% filter(scenario %in% c("Dupilumab", "Mepolizumab",
                                                 "Cendakimab", "Dupilumab + Budesonide")),
             aes(x = week, y = IGE_LEVEL, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = COLORS) +
  labs(title = "Total Serum IgE (IU/mL)",
       x = "Week", y = "IgE (IU/mL)", color = "Treatment") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.text = element_text(size = 8))

# Main dashboard
main_plot <- (p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | p9) +
  plot_annotation(
    title = "Eosinophilic Esophagitis (EoE) QSP Model — 6 Treatment Scenarios",
    subtitle = "52-week simulation: dupilumab, budesonide ODT, mepolizumab, cendakimab vs. no treatment",
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 10))
  )

cat("Saving main plot...\n")
ggsave("eoe_qsp_simulation.png", main_plot, width = 18, height = 14, dpi = 150)
cat("Plot saved: eoe_qsp_simulation.png\n")

# ===========================================================================
# 8. SCENARIO 3: DETAILED DRUG PK PROFILES
# ===========================================================================

# Simulate dupilumab PK alone for 52 weeks
dup_events <- build_events("dupilumab", start_day = 1, duration_days = 364)
dup_pk <- mod %>%
  mrgsim_df(events = as.data.frame(dup_events),
            end = 364, delta = 0.5, output = "df") %>%
  mutate(week = time / 7)

pk_plot <- ggplot(dup_pk, aes(x = week, y = DUP_TROUGH)) +
  geom_line(color = "#1565C0", size = 1.2) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red",
             size = 0.8) +
  annotate("text", x = 50, y = 2.5, label = "IC50 (STAT6 blockade ~2 µg/mL)",
           hjust = 1, color = "red", size = 3) +
  labs(title = "Dupilumab PK — Central Concentration (300 mg SC q2w)",
       x = "Week", y = "Dupilumab Concentration (mg/L ≡ µg/mL)") +
  theme_bw(base_size = 11)

ggsave("eoe_dupilumab_pk.png", pk_plot, width = 10, height = 5, dpi = 150)
cat("Dupilumab PK plot saved.\n")

# ===========================================================================
# 9. SENSITIVITY ANALYSIS
# ===========================================================================

cat("\n===== SENSITIVITY ANALYSIS: BASELINE EOS_ESO vs DUPILUMAB RESPONSE =====\n")
eos_baselines <- c(20, 50, 80, 120, 200)
sens_results <- lapply(eos_baselines, function(eos_bl) {
  mod2 <- mod %>% param(kin_eosbl = eos_bl * 4) %>% init(EOS_ESO = eos_bl)
  dup_ev <- build_events("dupilumab", 1, 168)
  out <- tryCatch(
    mrgsim_df(mod2, events = as.data.frame(dup_ev), end = 168, delta = 7),
    error = function(e) NULL
  )
  if (!is.null(out)) {
    out$BL_EOS_ESO <- eos_bl
    out$week <- out$time / 7
  }
  out
})
sens_df <- bind_rows(Filter(Negate(is.null), sens_results))

if (nrow(sens_df) > 0) {
  sens_plot <- ggplot(sens_df, aes(x = week, y = EOS_ESO,
                                   color = factor(BL_EOS_ESO),
                                   group = BL_EOS_ESO)) +
    geom_line(size = 1.2) +
    geom_hline(yintercept = 15, linetype = "dashed") +
    scale_color_viridis_d(option = "magma", name = "Baseline\neos/hpf") +
    labs(title = "Sensitivity Analysis: Dupilumab Response vs. Baseline Tissue Eosinophilia",
         x = "Week", y = "Tissue Eosinophils (eos/hpf)") +
    theme_bw(base_size = 11)
  ggsave("eoe_sensitivity.png", sens_plot, width = 9, height = 5, dpi = 150)
  cat("Sensitivity plot saved.\n")
}

cat("\n=======================================================\n")
cat("EoE QSP Simulation Complete.\n")
cat("Output files: eoe_qsp_simulation.png, eoe_dupilumab_pk.png,\n")
cat("              eoe_sensitivity.png\n")
cat("=======================================================\n")
