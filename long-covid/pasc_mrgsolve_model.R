################################################################################
## Long COVID (PASC) — mrgsolve QSP Model
## Post-Acute Sequelae of SARS-CoV-2 Infection
##
## Model Architecture:
##   22 ODE compartments spanning:
##     - Viral kinetics (acute → reservoir → antigen)
##     - Immune dysregulation (T/B cell, IFN, cytokines)
##     - Endothelial/coagulation (microthrombi, D-dimer)
##     - Neuroinflammation (BBB, microglia, serotonin)
##     - Autonomic dysfunction (POTS)
##     - Mitochondrial/energy (ATP, ROS, lactate)
##     - PK/PD: nirmatrelvir, metformin, SSRIs, LDN, anticoagulants
##
## Parameter Calibration Sources:
##   - Viral kinetics: Ke et al. 2021 Nature (NHM); Kissler et al. 2021 Science
##   - Immune: Phetsouphanh et al. 2022 Nature Immunology
##   - Microthrombi: Pretorius et al. 2021 Cardiovasc Diabetol
##   - Neuro: Fernández-Castañeda et al. 2022 Cell
##   - POTS: Dani et al. 2021 Clinical Medicine
##   - Metformin: COVID-OUT RCT (Bramante 2023 NEJM Evid)
##   - Nirmatrelvir: Hammond et al. 2022 NEJM; STOP-PASC trial
##
## Treatment Scenarios (7 total):
##   S1 = No treatment (natural history)
##   S2 = Extended nirmatrelvir/ritonavir (15-day course)
##   S3 = Metformin 500mg BID
##   S4 = Low-dose naltrexone (LDN) 4.5mg QD
##   S5 = Sertraline 50mg QD
##   S6 = Nirmatrelvir + Metformin (combination)
##   S7 = Full combination: Nirmatrelvir + Metformin + LDN + Sertraline
##
## Units: time = days, concentrations = ng/mL or normalized (0–1 scale)
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

pasc_code <- '
$PROB Long COVID (PASC) QSP Model — mrgsolve

$PARAM @annotated
// ---- Viral Kinetics ----
kViral    : 0.80  : /day : viral replication rate constant
kClear    : 1.20  : /day : viral clearance rate (innate immune)
kReservoir: 0.05  : /day : rate of viral seeding to tissue reservoirs
kActiv    : 0.02  : /day : reservoir reactivation rate (antigen persistence)
kAntigen  : 0.30  : /day : antigen clearance rate
IC50_nirm : 0.003 : ug/mL: nirmatrelvir IC50 for viral replication (3 nM)

// ---- Immune Dysregulation ----
kIFN      : 0.80  : /day : IFN-I induction by viral antigen
kdIFN     : 0.25  : /day : IFN-I decay
kCD8exh   : 0.10  : /day : CD8 T cell exhaustion induction
kdCD8     : 0.05  : /day : CD8 T cell recovery rate
kAutoAb   : 0.04  : /day : autoantibody production rate
kdAutoAb  : 0.008 : /day : autoantibody decay
kIL6      : 0.50  : /day : IL-6 production (from cytokine cascade)
kdIL6     : 0.35  : /day : IL-6 degradation
kTNF      : 0.30  : /day : TNF-alpha production
kdTNF     : 0.45  : /day : TNF-alpha degradation
IL6_basal : 0.10  : ng/mL: IL-6 baseline (healthy = ~1.5 pg/mL normalized)
EC50_met_IL6 : 500 : ng/mL: metformin EC50 for IL-6 suppression

// ---- Endothelial / Coagulation ----
kFibrin   : 0.15  : /day : fibrin microclot formation rate
kdFibrin  : 0.08  : /day : fibrin microclot lysis rate
kDdimer   : 0.20  : /day : D-dimer generation from fibrin
kdDdimer  : 0.12  : /day : D-dimer clearance
kAnticoag : 2.00  : /day : anticoagulant effect on fibrin (PD)

// ---- Neuroinflammation / BBB ----
kBBB      : 0.20  : /day : BBB disruption rate (from IL-6, TNF)
kdBBB     : 0.05  : /day : BBB restoration rate
kMicroglia: 0.40  : /day : microglial activation (from BBB breach + IL-6)
kdMicroglia:0.10  : /day : microglial deactivation
k5HT      : 0.15  : /day : serotonin depletion rate (neuroinflam-driven)
kd5HT     : 0.08  : /day : serotonin recovery rate
kSSRI_5HT : 0.25  : /day : SSRI effect on 5-HT restoration
kLDN_MG   : 1.50  : /day : LDN EC50 for microglial suppression (PD)
EC50_LDN  : 2.00  : ng/mL: LDN EC50 (plasma)

// ---- Autonomic / POTS ----
kAutoNom  : 0.10  : /day : autonomic dysfunction onset (from autoAb + SFN)
kdAutoNom : 0.03  : /day : autonomic function recovery
kPOTS     : 0.25  : /day : POTS severity development

// ---- Mitochondrial / Energy ----
kROS      : 0.30  : /day : ROS accumulation (driven by mito dysfunction)
kdROS     : 0.20  : /day : ROS scavenging
kMitoDmg  : 0.12  : /day : mitochondrial damage rate (ROS-driven)
kdMitoDmg : 0.04  : /day : mitochondrial recovery (biogenesis)
kLactate  : 0.20  : /day : lactate accumulation (anaerobic shift)
kdLactate : 0.30  : /day : lactate clearance
kMet_mito : 0.60  : /day : metformin AMPK-mediated mito protection

// ---- Drug PK Parameters ----
// Nirmatrelvir/Ritonavir
F_nirm    : 0.74  :      : nirmatrelvir oral bioavailability
ka_nirm   : 1.50  : /day : nirmatrelvir absorption rate
CL_nirm   : 8.00  : L/h  : nirmatrelvir total clearance (boosted by RTV)
Vd_nirm   : 105.0 : L    : volume of distribution
// Metformin
F_met     : 0.55  :      : metformin bioavailability
ka_met    : 0.70  : /day : absorption rate constant
CL_met    : 30.0  : L/h  : renal clearance
Vd_met    : 654.0 : L    : volume of distribution (extensive tissue)
// Sertraline
F_sert    : 0.44  :      : sertraline bioavailability
ka_sert   : 0.50  : /day : absorption rate
CL_sert   : 2.14  : L/h  : CL (hepatic, CYP2D6)
Vd_sert   : 2052.0: L    : Vd (lipophilic)
// Low-Dose Naltrexone
F_LDN     : 0.96  :      : naltrexone oral bioavailability
ka_LDN    : 2.40  : /day : absorption rate
CL_LDN    : 95.0  : L/h  : hepatic CL (high extraction)
Vd_LDN    : 1340.0: L    : Vd

// ---- Dosing switches (0=off, 1=on) ----
use_nirm  : 0  : : administer nirmatrelvir
use_met   : 0  : : administer metformin
use_sert  : 0  : : administer sertraline
use_LDN   : 0  : : administer LDN
use_anticoag : 0 :: administer anticoagulation

$CMT @annotated
// Viral
V_GUT    : viral depot (gut absorption)
V_PLASMA : viral load in plasma (RNA equiv)
V_RES    : tissue viral reservoir
V_AG     : persistent viral antigen
// Immune
IFN      : type I interferon (normalized)
CD8_exh  : CD8 T cell exhaustion index
Auto_Ab  : autoantibody level (normalized AU)
IL6      : IL-6 (pg/mL scale)
TNF      : TNF-alpha (normalized)
// Vascular
Fibrin   : fibrin microclot burden (normalized)
Ddimer   : D-dimer level (normalized)
// Neuro
BBB      : BBB disruption index (0=intact, 1=max disruption)
Microglia: microglial activation index
Serotonin: CNS serotonin (normalized, 1=normal)
// Autonomic
AutNom   : autonomic dysfunction index
// Mitochondria
ROS      : reactive oxygen species
MitoDmg  : mitochondrial damage index
Lactate  : blood lactate (normalized)
// Drug PK
A_nirm   : nirmatrelvir gut depot (ug)
C_nirm   : nirmatrelvir central (ug/L = ng/mL)
A_met    : metformin gut depot (ug)
C_met    : metformin central (ug/L = ng/mL)
A_sert   : sertraline gut depot (ug)
C_sert   : sertraline central (ug/L = ng/mL)
A_LDN    : LDN gut depot (ug)
C_LDN    : LDN central (ug/L = ng/mL)

$MAIN
// ---- PK-derived quantities ----
double ke_nirm = (CL_nirm * 1000) / (24.0 * Vd_nirm); // /day
double ke_met  = (CL_met  * 1000) / (24.0 * Vd_met);
double ke_sert = (CL_sert * 1000) / (24.0 * Vd_sert);
double ke_LDN  = (CL_LDN  * 1000) / (24.0 * Vd_LDN);

$ODE
// -------- PK --------
// Nirmatrelvir
dxdt_A_nirm = -ka_nirm * A_nirm;
dxdt_C_nirm = use_nirm * (ka_nirm * A_nirm * F_nirm / Vd_nirm) - ke_nirm * C_nirm;

// Metformin
dxdt_A_met  = -ka_met * A_met;
dxdt_C_met  = use_met * (ka_met * A_met * F_met / Vd_met) - ke_met * C_met;

// Sertraline
dxdt_A_sert = -ka_sert * A_sert;
dxdt_C_sert = use_sert * (ka_sert * A_sert * F_sert / Vd_sert) - ke_sert * C_sert;

// LDN
dxdt_A_LDN  = -ka_LDN * A_LDN;
dxdt_C_LDN  = use_LDN * (ka_LDN * A_LDN * F_LDN / Vd_LDN) - ke_LDN * C_LDN;

// -------- PD / Disease --------

// Nirmatrelvir viral suppression (Hill equation, n=1)
double nirm_eff  = C_nirm / (IC50_nirm + C_nirm);  // 0→1

// Metformin effects
double met_IL6   = 1.0 - C_met / (EC50_met_IL6 + C_met); // IL-6 reduction factor
double met_mito  = C_met / (EC50_met_IL6 + C_met) * kMet_mito;

// SSRI effect on 5-HT
double sert_5HT_eff = C_sert / (50.0 + C_sert);  // EC50 ~50 ng/mL for SERT occupancy

// LDN microglial suppression
double LDN_eff   = C_LDN / (EC50_LDN + C_LDN);

// ---- Viral ----
// Assumes patient post-acute: starts with low viral load from reservoir
dxdt_V_PLASMA = kViral * (1 - nirm_eff) * V_PLASMA * (1 - V_PLASMA/100)
              - kClear * IFN * V_PLASMA
              - kReservoir * V_PLASMA;
dxdt_V_RES    = kReservoir * V_PLASMA - kActiv * V_RES;
dxdt_V_AG     = kActiv * V_RES - kAntigen * V_AG;

// ---- Immune ----
dxdt_IFN      = kIFN * V_AG * (1 - IFN) - kdIFN * IFN;
dxdt_CD8_exh  = kCD8exh * IFN * (1 - CD8_exh) - kdCD8 * (1 - IFN) * CD8_exh;
dxdt_Auto_Ab  = kAutoAb * V_AG * (1 - Auto_Ab) - kdAutoAb * Auto_Ab;
dxdt_IL6      = (kIL6 * (V_AG + 0.5*Auto_Ab) * met_IL6) - kdIL6 * (IL6 - IL6_basal);
if (IL6 < IL6_basal) dxdt_IL6 = 0;
dxdt_TNF      = kTNF * (V_AG + 0.3*IL6) * (1 - TNF) - kdTNF * TNF;

// ---- Vascular ----
dxdt_Fibrin   = kFibrin * (IL6 * 0.5 + Auto_Ab * 0.5) * (1 - Fibrin)
              - kdFibrin * Fibrin
              - use_anticoag * kAnticoag * Fibrin;
dxdt_Ddimer   = kDdimer * Fibrin - kdDdimer * Ddimer;

// ---- Neuroinflammation ----
dxdt_BBB      = kBBB * (IL6 + TNF + Fibrin) / 3.0 * (1 - BBB)
              - kdBBB * BBB;
dxdt_Microglia= kMicroglia * BBB * (1 - Microglia) * (1 - LDN_eff)
              - kdMicroglia * (1 - BBB) * Microglia;
dxdt_Serotonin= kd5HT * (1 - Serotonin) + kSSRI_5HT * sert_5HT_eff * (1 - Serotonin)
              - k5HT * Microglia * Serotonin;

// ---- Autonomic ----
dxdt_AutNom   = kAutoNom * (Auto_Ab + 0.5 * CD8_exh) * (1 - AutNom)
              - kdAutoNom * AutNom;

// ---- Mitochondria ----
dxdt_ROS      = kROS * (Microglia + BBB + MitoDmg) / 3.0 - kdROS * ROS;
dxdt_MitoDmg  = kMitoDmg * ROS * (1 - MitoDmg) - kdMitoDmg * MitoDmg
              - met_mito * MitoDmg * 0.1;
dxdt_Lactate  = kLactate * MitoDmg - kdLactate * Lactate;

$TABLE
// ---- Composite Endpoints / Biomarkers ----

// Fatigue Severity Score (FSS, 1-7 scale; higher=worse)
double FSS = 1.0 + 6.0 * (0.4*MitoDmg + 0.3*Lactate + 0.2*AutNom + 0.1*Microglia);
if (FSS > 7.0) FSS = 7.0;

// VO2max (% predicted; normal=100)
double VO2max_pct = 100.0 * (1.0 - 0.5*MitoDmg - 0.3*Fibrin - 0.2*Lactate);
if (VO2max_pct < 10.0) VO2max_pct = 10.0;

// Cognitive score (MoCA proxy, 0-30; higher=better)
double MoCA = 30.0 * (1.0 - 0.5*(1-Serotonin) - 0.3*Microglia - 0.2*BBB);
if (MoCA < 0) MoCA = 0;

// POTS score (orthostatic HR change in bpm)
double POTS_HR_delta = 30.0 * AutNom + 15.0 * (1.0 - Serotonin);
if (POTS_HR_delta < 0) POTS_HR_delta = 0;

// SF-36 Physical Component Score (PCS, 0-100)
double SF36_PCS = 100.0 - 30.0*FSS/7.0 - 20.0*MitoDmg - 15.0*AutNom
                        - 15.0*(1-VO2max_pct/100.0) - 10.0*Fibrin - 10.0*Ddimer;
if (SF36_PCS < 0) SF36_PCS = 0;

// Neurofilament light (NfL, pg/mL; normal<10)
double NfL = 5.0 + 45.0 * (0.5*BBB + 0.3*Microglia + 0.2*(1-Serotonin));

// Dyspnea score (mMRC 0-4)
double dyspnea = 4.0 * (0.4*Fibrin + 0.3*MitoDmg + 0.3*(1-VO2max_pct/100.0));
if (dyspnea > 4.0) dyspnea = 4.0;

capture FSS       = FSS;
capture VO2max    = VO2max_pct;
capture MoCA      = MoCA;
capture POTS_HR   = POTS_HR_delta;
capture SF36_PCS  = SF36_PCS;
capture NfL_pg    = NfL;
capture mMRC      = dyspnea;
capture CRP_proxy = IL6 * 5.0;  // CRP ~ 5x IL-6 (normalized)
capture Ddimer_out= Ddimer;
capture Viral     = V_PLASMA;
capture Antigen   = V_AG;
capture Reservoir = V_RES;

$CAPTURE C_nirm C_met C_sert C_LDN

$INIT
// Starting state: 3-4 months post-acute COVID (established PASC)
V_GUT     = 0.0
V_PLASMA  = 0.5    // low-level viral RNA persisting
V_RES     = 3.0    // established tissue reservoirs
V_AG      = 0.8    // significant persistent antigen
IFN       = 0.4    // partially dysfunctional IFN response
CD8_exh   = 0.55   // significant T cell exhaustion
Auto_Ab   = 0.35   // moderate autoantibodies
IL6       = 0.50   // elevated IL-6 (norm ~2x baseline = 0.5 normalized)
TNF       = 0.30
Fibrin    = 0.45   // substantial microthrombus burden
Ddimer    = 0.40
BBB       = 0.40
Microglia = 0.50
Serotonin = 0.55   // depleted
AutNom    = 0.45   // significant autonomic dysfunction
ROS       = 0.40
MitoDmg   = 0.45
Lactate   = 0.35
A_nirm    = 0.0
C_nirm    = 0.0
A_met     = 0.0
C_met     = 0.0
A_sert    = 0.0
C_sert    = 0.0
A_LDN     = 0.0
C_LDN     = 0.0
'

mod <- mcode("pasc_qsp", pasc_code, quiet = TRUE)

################################################################################
## DOSING REGIMENS
################################################################################

# Base events shared by all scenarios
# Nirmatrelvir 300mg BID (q12h), 15-day course
nirm_dose <- 300 * 1000  # ug (per dose)
ev_nirm <- ev(cmt="A_nirm", amt=nirm_dose, ii=0.5, addl=29, time=0)  # 15 days BID

# Metformin 500mg BID (maintenance)
met_dose <- 500 * 1000   # ug
ev_met  <- ev(cmt="A_met",  amt=met_dose, ii=0.5, addl=364*2, time=0)  # 1 year

# Sertraline 50mg QD
sert_dose <- 50 * 1000   # ug
ev_sert <- ev(cmt="A_sert", amt=sert_dose, ii=1, addl=364, time=0)

# LDN 4.5mg QD
LDN_dose <- 4.5 * 1000   # ug
ev_LDN  <- ev(cmt="A_LDN",  amt=LDN_dose, ii=1, addl=364, time=0)

################################################################################
## SIMULATION FUNCTION
################################################################################

simulate_scenario <- function(scenario_name, use_nirm_=0, use_met_=0,
                                use_sert_=0, use_LDN_=0, use_ac_=0,
                                events=NULL) {
  p <- list(use_nirm=use_nirm_, use_met=use_met_,
            use_sert=use_sert_, use_LDN=use_LDN_,
            use_anticoag=use_ac_)
  if (is.null(events)) {
    out <- mod %>% param(p) %>%
      mrgsim(end=365, delta=1) %>%
      as.data.frame()
  } else {
    out <- mod %>% param(p) %>%
      mrgsim(events=events, end=365, delta=1) %>%
      as.data.frame()
  }
  out$scenario <- scenario_name
  out
}

################################################################################
## SCENARIO DEFINITIONS
################################################################################

# S1: Natural history (no treatment)
s1 <- simulate_scenario("S1: No Treatment", 0, 0, 0, 0, 0)

# S2: Extended Nirmatrelvir 15-day course
s2 <- simulate_scenario("S2: Nirmatrelvir 15d", 1, 0, 0, 0, 0,
                         events=ev_nirm)

# S3: Metformin 500mg BID (1 year)
s3 <- simulate_scenario("S3: Metformin", 0, 1, 0, 0, 0,
                         events=ev_met)

# S4: LDN 4.5mg QD
s4 <- simulate_scenario("S4: LDN 4.5mg", 0, 0, 0, 1, 0,
                         events=ev_LDN)

# S5: Sertraline 50mg QD
s5 <- simulate_scenario("S5: Sertraline 50mg", 0, 0, 1, 0, 0,
                         events=ev_sert)

# S6: Nirmatrelvir + Metformin
s6_ev <- c(ev_nirm, ev_met)
s6 <- simulate_scenario("S6: Nirm + Metformin", 1, 1, 0, 0, 0,
                         events=s6_ev)

# S7: Full combination (Nirm + Met + LDN + SSRI)
s7_ev <- c(ev_nirm, ev_met, ev_LDN, ev_sert)
s7 <- simulate_scenario("S7: Full Combo", 1, 1, 1, 1, 0,
                         events=s7_ev)

all_scenarios <- bind_rows(s1, s2, s3, s4, s5, s6, s7)

################################################################################
## VIRTUAL PATIENT POPULATION (n=200)
################################################################################

set.seed(2026)
n_vp <- 200

vp_params <- tibble(
  ID       = 1:n_vp,
  kViral   = rlnorm(n_vp, log(0.80), 0.30),
  kClear   = rlnorm(n_vp, log(1.20), 0.25),
  kMitoDmg = rlnorm(n_vp, log(0.12), 0.35),
  kIL6     = rlnorm(n_vp, log(0.50), 0.30),
  kAutoAb  = rlnorm(n_vp, log(0.04), 0.40),
  kAutoNom = rlnorm(n_vp, log(0.10), 0.35)
)

# Initial states variability (e.g., severity of PASC)
init_var <- tibble(
  ID         = 1:n_vp,
  V_RES_init = rlnorm(n_vp, log(3.0), 0.5),
  CD8_init   = pmin(0.95, rlnorm(n_vp, log(0.55), 0.25)),
  Mito_init  = pmin(0.95, rlnorm(n_vp, log(0.45), 0.30)),
  AutNom_init= pmin(0.95, rlnorm(n_vp, log(0.45), 0.35))
)

simulate_vp <- function(trt_name, use_nirm_, use_met_, use_sert_, use_LDN_,
                         events_=NULL) {
  purrr::map_dfr(1:n_vp, function(i) {
    p_i <- list(
      kViral   = vp_params$kViral[i],
      kClear   = vp_params$kClear[i],
      kMitoDmg = vp_params$kMitoDmg[i],
      kIL6     = vp_params$kIL6[i],
      kAutoAb  = vp_params$kAutoAb[i],
      kAutoNom = vp_params$kAutoNom[i],
      use_nirm = use_nirm_, use_met = use_met_,
      use_sert = use_sert_, use_LDN  = use_LDN_
    )
    init_i <- list(
      V_RES    = init_var$V_RES_init[i],
      CD8_exh  = init_var$CD8_init[i],
      MitoDmg  = init_var$Mito_init[i],
      AutNom   = init_var$AutNom_init[i]
    )
    if (is.null(events_)) {
      out <- mod %>% param(p_i) %>% init(init_i) %>%
        mrgsim(end=365, delta=7) %>% as.data.frame()
    } else {
      out <- mod %>% param(p_i) %>% init(init_i) %>%
        mrgsim(events=events_, end=365, delta=7) %>% as.data.frame()
    }
    out$ID <- i; out$trt <- trt_name; out
  })
}

cat("Simulating virtual patients for S1 and S7...\n")
vp_s1 <- simulate_vp("S1: No Treatment", 0, 0, 0, 0)
vp_s7 <- simulate_vp("S7: Full Combo",   1, 1, 1, 1, c(ev_nirm, ev_met, ev_LDN, ev_sert))

################################################################################
## OUTCOME SUMMARY TABLES
################################################################################

week52_endpoints <- all_scenarios %>%
  filter(time == 364) %>%
  select(scenario, FSS, VO2max, MoCA, POTS_HR, SF36_PCS, NfL_pg, mMRC,
         CRP_proxy, Ddimer_out) %>%
  mutate(across(where(is.numeric), ~round(., 2)))

cat("\n===== Week 52 Endpoint Summary =====\n")
print(week52_endpoints)

# Response rates (FSS ≤ 4 = responder, improvement ≥ 2 points)
baseline_FSS <- s1 %>% filter(time == 0) %>% pull(FSS) %>% mean()
response_table <- all_scenarios %>%
  filter(time == 364) %>%
  group_by(scenario) %>%
  summarise(
    FSS_mean        = mean(FSS),
    FSS_improve     = baseline_FSS - mean(FSS),
    pct_FSS_responder = mean(FSS <= 4) * 100,
    VO2max_mean     = mean(VO2max),
    MoCA_mean       = mean(MoCA),
    POTS_HR_mean    = mean(POTS_HR),
    .groups="drop"
  )

cat("\n===== Response Rates at Week 52 =====\n")
print(response_table)

################################################################################
## VISUALIZATIONS
################################################################################

colors_scenarios <- c(
  "S1: No Treatment"   = "#E74C3C",
  "S2: Nirmatrelvir 15d" = "#E67E22",
  "S3: Metformin"      = "#F1C40F",
  "S4: LDN 4.5mg"      = "#2ECC71",
  "S5: Sertraline 50mg"= "#3498DB",
  "S6: Nirm + Metformin" = "#9B59B6",
  "S7: Full Combo"     = "#1ABC9C"
)

# Plot 1: Fatigue severity over time
p1 <- ggplot(all_scenarios, aes(x=time, y=FSS, color=scenario)) +
  geom_line(size=1.1) +
  geom_hline(yintercept=4, linetype="dashed", color="gray50") +
  annotate("text", x=5, y=3.85, label="Responder threshold", size=3, color="gray50") +
  scale_color_manual(values=colors_scenarios) +
  labs(title="Fatigue Severity Score (FSS)", x="Day", y="FSS (1-7 scale)",
       color="Treatment") +
  theme_bw(base_size=12) + theme(legend.position="bottom")

# Plot 2: VO2max recovery
p2 <- ggplot(all_scenarios, aes(x=time, y=VO2max, color=scenario)) +
  geom_line(size=1.1) +
  scale_color_manual(values=colors_scenarios) +
  labs(title="VO2max (% Predicted)", x="Day", y="VO2max %",
       color="Treatment") +
  theme_bw(base_size=12) + theme(legend.position="none")

# Plot 3: Cognitive function
p3 <- ggplot(all_scenarios, aes(x=time, y=MoCA, color=scenario)) +
  geom_line(size=1.1) +
  geom_hline(yintercept=26, linetype="dashed", color="gray50") +
  annotate("text", x=5, y=25.7, label="MoCI (≥26 = normal)", size=3, color="gray50") +
  scale_color_manual(values=colors_scenarios) +
  labs(title="MoCA Cognitive Score", x="Day", y="MoCA (0-30)",
       color="Treatment") +
  theme_bw(base_size=12) + theme(legend.position="none")

# Plot 4: POTS HR delta
p4 <- ggplot(all_scenarios, aes(x=time, y=POTS_HR, color=scenario)) +
  geom_line(size=1.1) +
  geom_hline(yintercept=30, linetype="dashed", color="red", alpha=0.6) +
  annotate("text", x=5, y=29, label="POTS threshold (≥30 bpm)", size=3, color="red") +
  scale_color_manual(values=colors_scenarios) +
  labs(title="POTS: Orthostatic HR Change (bpm)", x="Day", y="ΔHR (bpm)",
       color="Treatment") +
  theme_bw(base_size=12) + theme(legend.position="none")

# Plot 5: SF-36 PCS
p5 <- ggplot(all_scenarios, aes(x=time, y=SF36_PCS, color=scenario)) +
  geom_line(size=1.1) +
  scale_color_manual(values=colors_scenarios) +
  labs(title="SF-36 Physical Component Score", x="Day", y="SF-36 PCS (0-100)",
       color="Treatment") +
  theme_bw(base_size=12) + theme(legend.position="none")

# Plot 6: PK profiles (day 1-15, nirmatrelvir)
pk_data <- all_scenarios %>%
  filter(scenario == "S2: Nirmatrelvir 15d", time <= 15) %>%
  select(time, C_nirm, C_met, C_sert, C_LDN) %>%
  pivot_longer(-time, names_to="drug", values_to="conc")

p6 <- ggplot(filter(all_scenarios, scenario == "S6: Nirm + Metformin", time <= 30),
             aes(x=time)) +
  geom_line(aes(y=C_nirm, color="Nirmatrelvir (ng/mL)"), size=1) +
  geom_line(aes(y=C_met/50, color="Metformin/50 (ng/mL)"), size=1, linetype=2) +
  scale_y_continuous("Nirmatrelvir (ng/mL)",
                     sec.axis = sec_axis(~.*50, name="Metformin (ng/mL)")) +
  scale_color_manual(values=c("Nirmatrelvir (ng/mL)"="#9B59B6",
                               "Metformin/50 (ng/mL)"="#F1C40F")) +
  labs(title="PK: Nirmatrelvir & Metformin (S6)", x="Day", color="Drug") +
  theme_bw(base_size=12) + theme(legend.position="bottom")

# Plot 7: VP distribution - FSS at week 52
vp_week52 <- bind_rows(vp_s1, vp_s7) %>%
  filter(time == 357) %>%
  select(ID, trt, FSS, VO2max, MoCA)

p7 <- ggplot(vp_week52, aes(x=FSS, fill=trt)) +
  geom_density(alpha=0.5, bw=0.3) +
  geom_vline(xintercept=4, linetype="dashed") +
  scale_fill_manual(values=c("S1: No Treatment"="#E74C3C", "S7: Full Combo"="#1ABC9C")) +
  labs(title="Virtual Population: FSS Distribution (Week 52)",
       x="Fatigue Severity Score", y="Density", fill="Treatment") +
  theme_bw(base_size=12)

# Combine plots
combined <- (p1 + p2) / (p3 + p4) / (p5 + p6)

cat("\nAll simulations complete. Generating plots...\n")

# Print key outputs
cat("\n=== KEY FINDINGS ===\n")
cat("Baseline FSS:", round(s1$FSS[1], 2), "\n")
cat("Week 52 FSS by scenario:\n")
print(week52_endpoints[, c("scenario", "FSS", "MoCA", "VO2max", "SF36_PCS")])

cat("\n=== CLINICAL TRIAL BENCHMARKS ===\n")
cat("COVID-OUT (Bramante 2023): Metformin reduced long COVID incidence by 41% (HR=0.59)\n")
cat("STOP-PASC: Extended nirmatrelvir pending results (NCT05595369)\n")
cat("RECOVER-VITAL: Nirmatrelvir x15d vs x5d for PASC (NCT05595369)\n")
cat("LDN PASC pilot: Ongoing trials at Stanford/Yale 2023-2025\n")

################################################################################
## SENSITIVITY ANALYSIS: Tornado plot
################################################################################

param_ranges <- list(
  kViral    = c(0.4, 1.6),
  kMitoDmg  = c(0.06, 0.24),
  kIL6      = c(0.25, 1.0),
  kAutoAb   = c(0.02, 0.08),
  kAutoNom  = c(0.05, 0.20),
  kROS      = c(0.15, 0.60)
)

sensitivity_results <- purrr::map_dfr(names(param_ranges), function(pname) {
  p_low  <- list(); p_low[[pname]]  <- param_ranges[[pname]][1]
  p_high <- list(); p_high[[pname]] <- param_ranges[[pname]][2]
  p_low[["use_nirm"]] <- 0; p_high[["use_nirm"]] <- 0
  p_low[["use_met"]]  <- 0; p_high[["use_met"]]  <- 0
  p_low[["use_sert"]] <- 0; p_high[["use_sert"]] <- 0
  p_low[["use_LDN"]]  <- 0; p_high[["use_LDN"]]  <- 0
  p_low[["use_anticoag"]] <- 0; p_high[["use_anticoag"]] <- 0

  fss_low  <- mod %>% param(p_low)  %>% mrgsim(end=364, delta=364) %>%
    as.data.frame() %>% tail(1) %>% pull(FSS)
  fss_high <- mod %>% param(p_high) %>% mrgsim(end=364, delta=364) %>%
    as.data.frame() %>% tail(1) %>% pull(FSS)
  tibble(param=pname, FSS_low=fss_low, FSS_high=fss_high,
         delta=abs(fss_high - fss_low))
})

cat("\n=== Sensitivity Analysis (impact on Week 52 FSS) ===\n")
print(sensitivity_results %>% arrange(desc(delta)))
