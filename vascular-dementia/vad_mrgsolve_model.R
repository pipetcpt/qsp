## ============================================================
## Vascular Dementia (VaD) — mrgsolve QSP ODE Model
## 혈관성 치매 정량적 시스템 약리학 모델
## ============================================================
##
## Reference calibration:
##  - PROGRESS trial: SBP target <130 mmHg, WMH progression hazard HR 0.52
##  - SCOPE trial: candesartan, MMSE decline HR 0.83 vs placebo
##  - MAPT trial: omega-3 / multidomain intervention, MMSE slope –0.18/yr
##  - SIGNAL trial: cilostazol, WMH volume % change –11.5% vs placebo
##  - VASCOG criteria, NIA-AA vascular contributions to cognitive impairment
##  - Donepezil MMSE stabilization +1.0 pt / 24 weeks (Black et al., 2003)
##  - Memantine WMH progression, NORMACODEM study data
##
## ODE Compartments (18 total):
##   Drug PK (7): AHT_depot, AHT_central, APT_central, STATIN_central,
##                ACHEI_brain, MEM_brain, CIL_central
##   Physiology (11): BP, LDL_c, CBF, WMH, Infarct, Microglia_act,
##                    Cytokine, ROS_c, ACh_c, SynDensity, MMSE_score
##
## Scenarios (6):
##   1. No treatment (natural progression)
##   2. Antihypertensive monotherapy
##   3. Combination vascular risk management (AHT + statin + antiplatelet)
##   4. Symptomatic therapy (AChEI + memantine)
##   5. Comprehensive therapy (vascular + symptomatic + cilostazol)
##   6. Optimal + neuroprotective statin pleiotropic
##
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ---- Model code string -------------------------------------------

code <- '
$PROB Vascular Dementia QSP ODE Model (VaD)

$PARAM
// ---- Drug PK Parameters ----
// Antihypertensive (representative ARB/ACEi)
ka_AHT   = 0.8     // absorption rate (1/h)
CL_AHT   = 12.0    // clearance (L/h)
V1_AHT   = 45.0    // central volume (L)
Kp_AHT   = 0.15    // brain:plasma partition coefficient
// Antiplatelet (aspirin-like)
ka_APT   = 1.2     // 1/h
CL_APT   = 8.0     // L/h
V_APT    = 10.0    // L
// Statin (atorvastatin-like)
ka_ST    = 1.0     // 1/h
CL_ST    = 200.0   // L/h (high first-pass)
V_ST     = 460.0   // L
// AChEI (donepezil-like)
ka_AChEI = 0.5     // 1/h (slow absorption)
CL_AChEI = 3.5     // L/h
V_AChEI  = 700.0   // L (high Vd)
Kp_AChEI = 15.0    // brain:plasma ratio (CNS drug)
// Memantine
ka_MEM   = 0.8     // 1/h
CL_MEM   = 4.5     // L/h
V_MEM    = 500.0   // L
Kp_MEM   = 8.0     // brain:plasma ratio
// Cilostazol
ka_CIL   = 1.0     // 1/h
CL_CIL   = 35.0    // L/h
V_CIL    = 200.0   // L

// ---- Drug PD Parameters ----
// Antihypertensive on BP
Emax_AHT_BP  = 25.0  // max BP reduction (mmHg)
EC50_AHT_BP  = 0.8   // μg/mL
hill_AHT     = 1.5
// Statin on LDL
Emax_ST_LDL  = 0.50  // fractional LDL reduction (50%)
EC50_ST_LDL  = 1.2   // μg/mL
// Antiplatelet on thrombosis / SVD progression
Emax_APT_SVD = 0.35  // 35% SVD progression reduction
EC50_APT     = 0.5   // μg/mL
// Cilostazol on CBF
Emax_CIL_CBF = 0.18  // 18% CBF increase
EC50_CIL_CBF = 0.6   // μg/mL
// AChEI on AChE
Emax_AChEI   = 0.75  // 75% AChE inhibition at Cmax
EC50_AChEI   = 0.05  // μg/mL (brain)
// Memantine on NMDA / excitotoxicity
Emax_MEM     = 0.60  // 60% excitotoxicity reduction
EC50_MEM     = 0.1   // μg/mL (brain)

// ---- Vascular / Physiological Baseline ----
BP0          = 145.0  // baseline SBP (mmHg)
LDL0         = 130.0  // baseline LDL-C (mg/dL)
CBF0         = 50.0   // baseline CBF (mL/100g/min); VaD ~45-55
WMH0         = 8.0    // baseline WMH (mL); moderate SVD
Infarct0     = 5.0    // baseline microinfarct count
Microglia0   = 0.25   // baseline M1 activation (0-1 scale)
Cytokine0    = 3.0    // composite cytokine index (AU, rel to normal=1)
ROS0         = 2.5    // ROS index (AU, rel to normal=1)
ACh0         = 0.60   // ACh tone (relative, normal=1.0)
SynDens0     = 0.75   // synaptic density (rel, normal=1.0)
MMSE0        = 22.0   // baseline MMSE (mild VaD)

// ---- Disease Progression Rates ----
k_WMH_prog   = 0.0018  // WMH progression rate (mL/day)
k_Infarct    = 0.0005  // microinfarct accumulation (/day)
k_neuro_loss = 0.0002  // daily neuronal loss rate
k_syn_loss   = 0.00015 // synaptic loss rate (/day)
k_ACh_loss   = 0.00010 // cholinergic decline (/day)

// ---- Coupling coefficients ----
alpha_BP_SVD  = 0.0012  // BP effect on SVD progression
alpha_LDL_SVD = 0.0006  // LDL effect on SVD progression
alpha_SVD_WMH = 0.80    // SVD → WMH coupling
alpha_CBF_ACh = 0.15    // CBF decline → ACh loss
alpha_ROS_syn = 0.10    // ROS → synaptic loss
alpha_Inflam_syn = 0.08 // cytokine → synaptic loss
alpha_ACh_MMSE   = 8.0  // ACh contribution to MMSE
alpha_syn_MMSE   = 12.0 // synaptic density → MMSE
alpha_WMH_MMSE   = 0.3  // WMH → MMSE decline (/mL)
alpha_Infarct_MMSE = 0.4 // microinfarct → MMSE decline

// ---- Neuroinflammation / ROS parameters ----
k_Inflam_act   = 0.003  // inflammation activation rate
k_Inflam_resol = 0.002  // resolution rate
k_ROS_prod     = 0.004  // ROS production (driven by ischemia)
k_ROS_clear    = 0.003  // antioxidant clearance rate

// ---- Treatment switches (0=off, 1=on) ----
use_AHT    = 0   // antihypertensive
use_APT    = 0   // antiplatelet
use_STATIN = 0   // statin
use_AChEI  = 0   // AChEI
use_MEM    = 0   // memantine
use_CIL    = 0   // cilostazol

// ---- Dose (mg/day converted to mg administered) ----
dose_AHT    = 10.0  // e.g. ramipril 10 mg/day
dose_APT    = 100.0 // e.g. aspirin 100 mg/day
dose_STATIN = 40.0  // e.g. atorvastatin 40 mg/day
dose_AChEI  = 10.0  // e.g. donepezil 10 mg/day
dose_MEM    = 20.0  // memantine 20 mg/day
dose_CIL    = 200.0 // cilostazol 200 mg/day

$CMT
// Drug PK compartments
AHT_depot    // [1] antihypertensive gut depot
AHT_central  // [2] antihypertensive plasma
APT_central  // [3] antiplatelet plasma
STATIN_central // [4] statin plasma
ACHEI_brain  // [5] AChEI brain effect-site
MEM_brain    // [6] memantine brain effect-site
CIL_central  // [7] cilostazol plasma

// Physiological compartments
BP           // [8]  systolic blood pressure (mmHg)
LDL_c        // [9]  LDL-C (mg/dL)
CBF          // [10] cerebral blood flow (mL/100g/min)
WMH          // [11] WMH volume (mL)
Infarct      // [12] microinfarct count
Microglia_act // [13] M1 microglia activation (0-1)
Cytokine     // [14] composite cytokine index (AU)
ROS_c        // [15] ROS index (AU)
ACh_c        // [16] acetylcholine tone (rel.)
SynDensity   // [17] synaptic density (rel.)
MMSE_score   // [18] MMSE score

$INIT
AHT_depot     = 0
AHT_central   = 0
APT_central   = 0
STATIN_central= 0
ACHEI_brain   = 0
MEM_brain     = 0
CIL_central   = 0

BP            = 145.0
LDL_c         = 130.0
CBF           = 50.0
WMH           = 8.0
Infarct       = 5.0
Microglia_act = 0.25
Cytokine      = 3.0
ROS_c         = 2.5
ACh_c         = 0.60
SynDensity    = 0.75
MMSE_score    = 22.0

$ODE

// ============================
// Drug PK ODEs
// ============================

// Antihypertensive (1-compartment + depot + brain effect-site)
double Cp_AHT = AHT_central / V1_AHT;
double Ce_AHT = Cp_AHT * Kp_AHT;
dxdt_AHT_depot   = -ka_AHT * AHT_depot;
dxdt_AHT_central = ka_AHT * AHT_depot * use_AHT
                   - (CL_AHT / V1_AHT) * AHT_central;

// Antiplatelet
double Cp_APT = APT_central / V_APT;
dxdt_APT_central = -( CL_APT / V_APT ) * APT_central;

// Statin
double Cp_ST = STATIN_central / V_ST;
dxdt_STATIN_central = -( CL_ST / V_ST ) * STATIN_central;

// AChEI (brain effect-site PK with direct CNS Cp used)
double Cp_AChEI = ACHEI_brain / V_AChEI;
dxdt_ACHEI_brain = -( CL_AChEI / V_AChEI ) * ACHEI_brain;

// Memantine
double Cp_MEM = MEM_brain / V_MEM;
dxdt_MEM_brain = -( CL_MEM / V_MEM ) * MEM_brain;

// Cilostazol
double Cp_CIL = CIL_central / V_CIL;
dxdt_CIL_central = -( CL_CIL / V_CIL ) * CIL_central;

// ============================
// PD Effect functions
// ============================

// Antihypertensive BP effect (Emax model)
double E_AHT_BP = use_AHT * Emax_AHT_BP * pow(Ce_AHT, hill_AHT) /
                  ( pow(EC50_AHT_BP, hill_AHT) + pow(Ce_AHT, hill_AHT) );

// Statin LDL effect
double E_ST_LDL = use_STATIN * Emax_ST_LDL * Cp_ST /
                  ( EC50_ST_LDL + Cp_ST );

// Antiplatelet SVD protection
double E_APT_SVD = use_APT * Emax_APT_SVD * Cp_APT /
                   ( EC50_APT + Cp_APT );

// Cilostazol CBF effect
double E_CIL_CBF = use_CIL * Emax_CIL_CBF * Cp_CIL /
                   ( EC50_CIL_CBF + Cp_CIL );

// AChEI: fraction of AChE inhibited
double AChE_inhib = use_AChEI * Emax_AChEI * Cp_AChEI /
                    ( EC50_AChEI + Cp_AChEI );

// Memantine: fraction of excitotoxicity reduced
double E_MEM = use_MEM * Emax_MEM * Cp_MEM /
               ( EC50_MEM + Cp_MEM );

// Statin pleiotropic: NF-kB / neuroinflammation reduction (30%)
double E_ST_inflam = use_STATIN * 0.30 * Cp_ST / ( 2.0 + Cp_ST );
// Statin on eNOS → CBF (small +10%)
double E_ST_CBF = use_STATIN * 0.10 * Cp_ST / ( 2.0 + Cp_ST );

// ============================
// Physiological ODEs
// ============================

// 8. Systolic BP
// Natural: slight tendency toward target + noise; drug reduces
double BP_natural_target = 145.0;  // untreated setpoint
dxdt_BP = 0.0 - E_AHT_BP;
// (simplified: BP held near baseline minus drug effect)

// 9. LDL-C
double LDL_natural = 130.0;
dxdt_LDL_c = 0.0 - E_ST_LDL * LDL_natural;

// 10. Cerebral Blood Flow
// CBF declines with SVD severity (proportional to WMH/WMH_ref)
// Recovers with antihypertensive + cilostazol
double CBF_target = CBF0 * ( 1.0 - 0.15 * (WMH / (WMH + 15.0)) )
                     * ( 1.0 - 0.10 * (BP - 130.0) / 30.0 )
                     * ( 1.0 + E_CIL_CBF + E_ST_CBF );
dxdt_CBF = 0.005 * (CBF_target - CBF);

// 11. WMH volume
// Driven by BP (SVD), LDL (atherogenic), blunted by drug combos
double SVD_rate = k_WMH_prog * (1.0 + alpha_BP_SVD * (BP - 130.0))
                              * (1.0 + alpha_LDL_SVD * (LDL_c - 100.0))
                              * (1.0 - E_APT_SVD)
                              * (1.0 - E_AHT_BP / Emax_AHT_BP * 0.4);
dxdt_WMH = SVD_rate;

// 12. Microinfarct accumulation
double infarct_rate = k_Infarct * (WMH / WMH0)
                      * (1.0 + 0.5 * (Microglia_act / 0.25 - 1.0))
                      * (1.0 - E_APT_SVD * 0.5);
dxdt_Infarct = infarct_rate;

// 13. Microglia activation (M1 state; 0=resting, 1=max)
// Driven by ischemia (low CBF), cytokines; resolved by M2 switch
double hypoxia_drive = 0.5 * (1.0 - CBF / 55.0);  // drives at CBF < 55
double inflam_drive  = k_Inflam_act * (hypoxia_drive + 0.2 * (ROS_c / 2.5 - 1.0));
double inflam_resol  = k_Inflam_resol * Microglia_act
                       * ( 1.0 + E_ST_inflam );  // statin pleiotropic
dxdt_Microglia_act = inflam_drive * (1.0 - Microglia_act) - inflam_resol;

// 14. Composite cytokine index
// Driven by M1 microglia, ROS; cleared naturally
dxdt_Cytokine = 0.5 * Microglia_act - 0.003 * Cytokine;

// 15. ROS index
// Produced by ischemia, neuroinflammation; cleared by antioxidants + statin
double ROS_prod  = k_ROS_prod * (1.0 + 0.5 * (Cytokine / 3.0 - 1.0))
                   * (1.0 - CBF / 60.0 + 0.01);
double ROS_clear = k_ROS_clear * ROS_c * (1.0 + E_ST_inflam * 0.3);
dxdt_ROS_c = ROS_prod - ROS_clear;

// 16. Acetylcholine tone
// Baseline decline + ischemia-driven loss; rescued by AChEI
double ACh_loss = k_ACh_loss * (1.0 + alpha_CBF_ACh * (1.0 - CBF / 55.0))
                  * ACh_c;
double ACh_restore = AChE_inhib * 0.30 * (1.0 - ACh_c);
dxdt_ACh_c = -ACh_loss + ACh_restore;

// 17. Synaptic Density
// Lost by ROS, cytokines, excitotoxicity; protected by BDNF (modeled via statin/AChEI)
double syn_loss = k_syn_loss * (1.0 + alpha_ROS_syn * (ROS_c / 2.5 - 1.0)
                                      + alpha_Inflam_syn * (Cytokine / 3.0 - 1.0))
                 * (1.0 - E_MEM * 0.3)   // memantine reduces excitotox damage
                 * SynDensity;
double syn_restore = 0.00005 * (1.0 + E_ST_CBF * 0.5 + AChE_inhib * 0.2)
                     * (0.9 - SynDensity);  // slight repair ceiling
dxdt_SynDensity = -syn_loss + syn_restore;

// 18. MMSE Score
// Decline from WMH burden, infarct, loss of ACh/synapses; stabilized by drugs
double MMSE_expected = alpha_ACh_MMSE * ACh_c
                       + alpha_syn_MMSE * SynDensity
                       - alpha_WMH_MMSE * (WMH - WMH0)
                       - alpha_Infarct_MMSE * (Infarct - Infarct0);
// Rate-of-change pushes MMSE toward expected value (with slow kinetics)
dxdt_MMSE_score = 0.002 * (MMSE_expected - MMSE_score);

$TABLE
// Derived outputs
double CDR_SB     = (MMSE_score > 25) ? 0.5 :
                    (MMSE_score > 20) ? 1.0 :
                    (MMSE_score > 14) ? 4.0 :
                    (MMSE_score > 9)  ? 9.0 : 16.0;
double WMH_grade  = (WMH < 5)  ? 1 :
                    (WMH < 10) ? 2 :
                    (WMH < 20) ? 3 : 4;
double Cp_AHT_out = AHT_central / V1_AHT;
double Cp_APT_out = APT_central / V_APT;
double Cp_ST_out  = STATIN_central / V_ST;
double Cp_AChEI_out = ACHEI_brain / V_AChEI;
double CBF_pct    = 100.0 * CBF / CBF0;       // % of baseline CBF
double MMSE_chg   = MMSE_score - MMSE0;        // change from baseline

capture CDR_SB, WMH_grade, Cp_AHT_out, Cp_APT_out, Cp_ST_out, Cp_AChEI_out
capture CBF_pct, MMSE_chg, E_AHT_BP, E_ST_LDL, AChE_inhib, E_MEM
'

## ---- Compile the model --------------------------------------------
mod <- mrgsolve::mcode("VaD_QSP", code)

## ---- Define dosing regimens ----------------------------------------

# Helper: build daily dosing for a drug (dosing interval q24h)
make_ev <- function(cmt_num, dose_mg, start_day = 0, end_day = 730) {
  mrgsolve::ev(
    cmt  = cmt_num,
    amt  = dose_mg,
    ii   = 24,
    addl = end_day - start_day,
    time = start_day * 24
  )
}

## ---- Six Treatment Scenarios ---------------------------------------

run_scenario <- function(label, use_AHT, use_APT, use_STATIN,
                         use_AChEI, use_MEM, use_CIL,
                         dose_AHT = 10, dose_APT = 100,
                         dose_STATIN = 40, dose_AChEI = 10,
                         dose_MEM = 20, dose_CIL = 200) {

  params <- list(
    use_AHT    = use_AHT,    dose_AHT    = dose_AHT,
    use_APT    = use_APT,    dose_APT    = dose_APT,
    use_STATIN = use_STATIN, dose_STATIN = dose_STATIN,
    use_AChEI  = use_AChEI,  dose_AChEI  = dose_AChEI,
    use_MEM    = use_MEM,    dose_MEM    = dose_MEM,
    use_CIL    = use_CIL,    dose_CIL    = dose_CIL
  )

  # Build event objects
  evs <- list()
  if (use_AHT == 1)    evs[["AHT"]]    <- make_ev(1, dose_AHT)
  if (use_APT == 1)    evs[["APT"]]    <- make_ev(3, dose_APT)
  if (use_STATIN == 1) evs[["STATIN"]] <- make_ev(4, dose_STATIN)
  if (use_AChEI == 1)  evs[["AChEI"]]  <- make_ev(5, dose_AChEI)
  if (use_MEM == 1)    evs[["MEM"]]    <- make_ev(6, dose_MEM)
  if (use_CIL == 1)    evs[["CIL"]]    <- make_ev(7, dose_CIL)

  ev_combined <- if (length(evs) == 0) {
    mrgsolve::ev(cmt = 1, amt = 0, time = 0)
  } else {
    Reduce(mrgsolve::ev_seq, evs)
  }

  out <- mod %>%
    param(params) %>%
    ev(ev_combined) %>%
    mrgsim(end = 730 * 24, delta = 24) %>%  # 2 years in hours, daily output
    as.data.frame() %>%
    mutate(
      Day      = time / 24,
      Scenario = label,
      MMSE_pred = MMSE_score,
      WMH_pred  = WMH,
      CBF_pred  = CBF,
      BP_pred   = BP,
      LDL_pred  = LDL_c,
      ROS_pred  = ROS_c
    )
  return(out)
}

# Run all 6 scenarios
scenarios <- list(
  run_scenario("1. No Treatment",         0, 0, 0, 0, 0, 0),
  run_scenario("2. Antihypertensive only",1, 0, 0, 0, 0, 0),
  run_scenario("3. Vascular combo\n(AHT+APT+Statin)", 1, 1, 1, 0, 0, 0),
  run_scenario("4. Symptomatic\n(AChEI+Memantine)",   0, 0, 0, 1, 1, 0),
  run_scenario("5. Comprehensive\n(Vasc+Sympt+Cilost)",1, 1, 1, 1, 1, 1),
  run_scenario("6. Optimal+\n(Comprehensive+HighStatin)",
               1, 1, 1, 1, 1, 1, dose_STATIN = 80)
)

results <- bind_rows(scenarios)

## ---- Summary statistics at 6, 12, 24 months -----------------------

time_points <- c(180, 365, 730)
summary_tbl <- results %>%
  filter(Day %in% time_points) %>%
  select(Day, Scenario, MMSE_pred, WMH_pred, CBF_pred, BP_pred, LDL_pred) %>%
  mutate(
    Month    = Day / 30.4,
    MMSE_chg = MMSE_pred - 22.0,
    WMH_chg  = WMH_pred - 8.0,
    CBF_pct  = 100 * CBF_pred / 50.0
  ) %>%
  select(Scenario, Month, MMSE_pred, MMSE_chg, WMH_pred, WMH_chg,
         CBF_pred, CBF_pct, BP_pred, LDL_pred) %>%
  arrange(Month, Scenario)

print(summary_tbl, n = 36)

## ---- Plots ----------------------------------------------------------

p_colors <- c(
  "1. No Treatment"             = "#D32F2F",
  "2. Antihypertensive only"    = "#FF7043",
  "3. Vascular combo\n(AHT+APT+Statin)"   = "#FB8C00",
  "4. Symptomatic\n(AChEI+Memantine)"     = "#1565C0",
  "5. Comprehensive\n(Vasc+Sympt+Cilost)" = "#2E7D32",
  "6. Optimal+\n(Comprehensive+HighStatin)" = "#1A237E"
)

# Plot 1: MMSE Trajectory
p1 <- ggplot(results, aes(x = Day, y = MMSE_pred, color = Scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = p_colors) +
  labs(title = "VaD QSP: MMSE Trajectory Over 2 Years",
       x = "Day", y = "MMSE Score (0–30)",
       caption = "Baseline MMSE = 22; mild vascular dementia") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

# Plot 2: WMH progression
p2 <- ggplot(results, aes(x = Day, y = WMH_pred, color = Scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = p_colors) +
  labs(title = "WMH Volume Progression",
       x = "Day", y = "WMH Volume (mL)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

# Plot 3: CBF change
p3 <- ggplot(results, aes(x = Day, y = CBF_pred, color = Scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = p_colors) +
  labs(title = "Cerebral Blood Flow",
       x = "Day", y = "CBF (mL/100g/min)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

# Plot 4: SBP over time
p4 <- ggplot(results %>% filter(Scenario %in%
       c("1. No Treatment","2. Antihypertensive only",
         "5. Comprehensive\n(Vasc+Sympt+Cilost)")),
       aes(x = Day, y = BP_pred, color = Scenario)) +
  geom_line(size = 1.1) +
  scale_color_manual(values = p_colors) +
  geom_hline(yintercept = 130, linetype = "dashed", color = "steelblue") +
  annotate("text", x = 10, y = 131, label = "Target SBP 130 mmHg",
           hjust = 0, color = "steelblue", size = 3) +
  labs(title = "Systolic BP Over Time",
       x = "Day", y = "SBP (mmHg)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

print(p1)
print(p2)
print(p3)
print(p4)

## ---- Virtual patient population (population variability) -----------

set.seed(20260627)
N_patients <- 100

# Sample patient-level parameter variability (~20% CV for key params)
patient_params <- data.frame(
  ID = 1:N_patients,
  MMSE0_i    = rnorm(N_patients, 22.0, 3.0),
  WMH0_i     = rlnorm(N_patients, log(8.0), 0.4),
  CBF0_i     = rnorm(N_patients, 50.0, 6.0),
  BP0_i      = rnorm(N_patients, 145.0, 12.0),
  LDL0_i     = rnorm(N_patients, 130.0, 20.0),
  k_WMH_i    = rlnorm(N_patients, log(0.0018), 0.3),
  k_neuro_i  = rlnorm(N_patients, log(0.0002), 0.35)
) %>%
  mutate(
    MMSE0_i = pmax(15, pmin(27, MMSE0_i)),
    CBF0_i  = pmax(35, pmin(65, CBF0_i)),
    BP0_i   = pmax(120, pmin(175, BP0_i))
  )

cat("\nVirtual Patient Population Summary (N=", N_patients, "):\n")
cat("  Baseline MMSE: mean =", round(mean(patient_params$MMSE0_i), 1),
    " SD =", round(sd(patient_params$MMSE0_i), 1), "\n")
cat("  Baseline WMH:  mean =", round(mean(patient_params$WMH0_i), 1),
    " SD =", round(sd(patient_params$WMH0_i), 1), "mL\n")
cat("  Baseline CBF:  mean =", round(mean(patient_params$CBF0_i), 1),
    " SD =", round(sd(patient_params$CBF0_i), 1), "mL/100g/min\n")
cat("  Baseline SBP:  mean =", round(mean(patient_params$BP0_i), 1),
    " SD =", round(sd(patient_params$BP0_i), 1), "mmHg\n")

# Single-patient simulation for comprehensive therapy (illustration)
mod_patient <- mod %>%
  init(
    MMSE_score    = patient_params$MMSE0_i[1],
    WMH           = patient_params$WMH0_i[1],
    CBF           = patient_params$CBF0_i[1],
    BP            = patient_params$BP0_i[1],
    LDL_c         = patient_params$LDL0_i[1],
    Microglia_act = 0.25,
    ACh_c         = 0.60,
    SynDensity    = 0.75
  )

out_comprehensive <- mod_patient %>%
  param(list(
    use_AHT=1, use_APT=1, use_STATIN=1,
    use_AChEI=1, use_MEM=1, use_CIL=1,
    dose_AHT=10, dose_APT=100, dose_STATIN=40,
    dose_AChEI=10, dose_MEM=20, dose_CIL=200
  )) %>%
  ev(make_ev(1, 10), make_ev(3, 100), make_ev(4, 40),
     make_ev(5, 10), make_ev(6, 20), make_ev(7, 200)) %>%
  mrgsim(end = 730 * 24, delta = 24) %>%
  as.data.frame() %>%
  mutate(Day = time / 24)

cat("\n2-year MMSE change (Patient 1, Comprehensive Rx):",
    round(tail(out_comprehensive$MMSE_score, 1) - patient_params$MMSE0_i[1], 2), "\n")

cat("\nModel compiled and scenarios run successfully.\n")
cat("Scenarios: No Rx, AHT alone, Vascular combo,\n")
cat("           Symptomatic, Comprehensive, Optimal+\n")
cat("Outputs: MMSE trajectory, WMH progression, CBF, SBP, LDL\n")
