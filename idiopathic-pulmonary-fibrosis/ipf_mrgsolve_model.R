## =============================================================================
## IPF QSP Model — mrgsolve ODE Implementation
## Disease: Idiopathic Pulmonary Fibrosis (IPF)
## Drugs  : Pirfenidone (ESBRIET) · Nintedanib (OFEV) · Combination
## Author : Claude Code Routine (CCR) · 2026-06-17
##
## Key References (calibration anchors):
##  - ASCEND (King et al. NEJM 2014): pirfenidone –47.9% FVC decline vs placebo
##  - INPULSIS-1/2 (Richeldi et al. NEJM 2014): nintedanib –50.1% FVC decline
##  - Pirfenidone PK: Rubino (2009) Clin Pharmacokinet; F=81%, t½=2.4h, ka=1.74/h
##  - Nintedanib PK: Stopfer (2011) Clin Pharmacokinet; F=4.7%, t½=10h, ka=0.8/h
##  - Natural history FVC decline ~200 mL/yr (Raghu 2011 ATS Guidelines)
## =============================================================================

library(mrgsolve)
library(tidyverse)
library(ggplot2)
library(gridExtra)

## ─── 1. MODEL DEFINITION ────────────────────────────────────────────────────

ipf_model_code <- '
$PROB IPF QSP Model — Pirfenidone & Nintedanib PK/PD

$PARAM @annotated
// ── Pirfenidone PK ──────────────────────────────────────────────────────────
ka_P   : 1.74   : Absorption rate constant pirfenidone (1/h)
F_P    : 0.81   : Oral bioavailability pirfenidone
CL_P   : 8.4    : Clearance pirfenidone (L/h)
V1_P   : 20.4   : Central volume pirfenidone (L)
V2_P   : 14.0   : Peripheral volume pirfenidone (L)
Q_P    : 3.6    : Inter-compartment clearance pirfenidone (L/h)

// ── Nintedanib PK ────────────────────────────────────────────────────────────
ka_N   : 0.80   : Absorption rate constant nintedanib (1/h)
F_N    : 0.047  : Oral bioavailability nintedanib
CL_N   : 22.0   : Clearance nintedanib (L/h)
V1_N   : 730.0  : Central volume nintedanib (L)
V2_N   : 900.0  : Peripheral volume nintedanib (L)
Q_N    : 15.0   : Inter-compartment clearance nintedanib (L/h)

// ── Disease PD — AEC & TGF-β ─────────────────────────────────────────────────
AEC2_ss : 1.0   : AEC-II steady-state (normalized)
k_AEC   : 0.005 : AEC-II damage rate (per hour, from ROS/injury)
k_rep   : 0.003 : AEC-II repair rate (per hour)
kprod_TGFb : 0.08 : TGF-β1 basal production rate
kdeg_TGFb  : 0.04 : TGF-β1 degradation rate (1/h)
kact_TGFb  : 0.12 : TGF-β1 activation from damaged AEC (proportionality)
EC50_P_TGFb: 30.0 : Pirfenidone EC50 for TGF-β inhibition (µg/mL)
Emax_P_TGFb: 0.65 : Pirfenidone Emax for TGF-β inhibition

// ── Disease PD — Fibroblast/Myofibroblast ─────────────────────────────────────
kact_F  : 0.015 : Fibroblast activation rate by TGF-β
kdiff_M : 0.025 : Myofibroblast differentiation rate from activated fibroblast
kapop_F : 0.008 : Fibroblast apoptosis (reduced in IPF)
F_ss    : 1.0   : Baseline fibroblast level
EC50_N_F: 15.0  : Nintedanib EC50 for fibroblast proliferation inhibition (nM)
Emax_N_F: 0.70  : Nintedanib Emax for fibroblast inhibition

// ── Disease PD — ECM / Collagen ───────────────────────────────────────────────
kprod_Col : 0.010 : Collagen production rate by myofibroblasts (per unit time)
kdeg_Col  : 0.002 : Collagen degradation (MMP-mediated)
kprod_MMP : 0.020 : MMP production rate
kdeg_MMP  : 0.030 : MMP degradation rate
kprod_TIMP: 0.018 : TIMP production rate (TGF-β driven)
kdeg_TIMP : 0.020 : TIMP degradation rate
Col_ss    : 1.0   : Baseline collagen (normalized)

// ── Disease PD — Macrophage ───────────────────────────────────────────────────
kM2_act   : 0.012 : M2 macrophage activation rate
kdeg_M2   : 0.008 : M2 macrophage clearance rate
EC50_P_M2 : 25.0  : Pirfenidone EC50 for M2/cytokine inhibition (µg/mL)

// ── Oxidative stress ──────────────────────────────────────────────────────────
kprod_ROS  : 0.04 : Basal ROS production rate
kdeg_ROS   : 0.05 : ROS clearance rate (GSH/antioxidant)
kfb_ROS    : 0.02 : ROS positive feedback from AEC damage
Emax_P_ROS : 0.45 : Pirfenidone Emax antioxidant effect

// ── Lung Function (FVC) ───────────────────────────────────────────────────────
FVC_base   : 80.0 : Baseline FVC % predicted (typical mild-moderate IPF)
k_FVC_loss : 0.0026 : FVC loss rate per unit collagen excess (per hour)
// ~200 mL/yr ≈ 2.5% predicted/yr ≈ 0.000285%/h → calibrated

// ── DLCO ─────────────────────────────────────────────────────────────────────
DLCO_base  : 65.0  : Baseline DLCO % predicted
k_DLCO_loss: 0.0018 : DLCO loss rate from AEC damage

// ── Molecular weights for unit conversion ────────────────────────────────────
MW_P  : 185.22 : Molecular weight pirfenidone (g/mol)
MW_N  : 539.63 : Molecular weight nintedanib (g/mol)

$CMT @annotated
// Pirfenidone PK
DEPOT_P  : Pirfenidone gut depot (µg)
CENT_P   : Pirfenidone central compartment (µg)
PERI_P   : Pirfenidone peripheral compartment (µg)
// Nintedanib PK
DEPOT_N  : Nintedanib gut depot (ng)
CENT_N   : Nintedanib central compartment (ng)
PERI_N   : Nintedanib peripheral compartment (ng)
// Disease PD
AEC2     : AEC-II cell population (normalized)
TGFb     : Active TGF-β1 level (normalized)
M2       : M2 macrophage activation (normalized)
ROS      : Reactive oxygen species (normalized)
FIBRO    : Activated fibroblast pool (normalized)
MYOFIB   : Myofibroblast pool (normalized)
COLLAGEN : ECM collagen accumulation (normalized)
MMP      : Matrix metalloproteinase activity (normalized)
TIMP     : TIMP activity (normalized)
FVC_st   : FVC state variable (% predicted)
DLCO_st  : DLCO state variable (% predicted)

$MAIN
// ── Derived PK concentrations ─────────────────────────────────────────────────
double Cp_P = CENT_P / V1_P;          // pirfenidone µg/mL
double Cn_N = CENT_N / V1_N * 1000.0; // nintedanib ng/mL → nM: /MW_N*1e6 simplified to nM

// nM conversion for nintedanib: (ng/mL) / MW_N * 1000 = nM
double Cn_nM = (CENT_N / V1_N) / MW_N * 1e3;  // nM

// ── Pharmacodynamic inhibition (Hill equation) ────────────────────────────────
double inh_P_TGFb = Emax_P_TGFb * Cp_P / (EC50_P_TGFb + Cp_P);
double inh_P_M2   = Emax_P_TGFb * Cp_P / (EC50_P_M2   + Cp_P);
double inh_P_ROS  = Emax_P_ROS  * Cp_P / (EC50_P_TGFb + Cp_P);
double inh_N_F    = Emax_N_F    * Cn_nM / (EC50_N_F    + Cn_nM);

// ── MMP:TIMP ratio → net ECM remodeling ──────────────────────────────────────
double net_ECM = (TIMP > 0) ? MMP / TIMP : 1.0;

// ── Initial conditions (set in $INIT) ────────────────────────────────────────
// All PD states start at baseline (normalized = 1 except collagen)

$INIT
DEPOT_P  = 0
CENT_P   = 0
PERI_P   = 0
DEPOT_N  = 0
CENT_N   = 0
PERI_N   = 0
AEC2     = 1.0
TGFb     = 1.0
M2       = 1.0
ROS      = 1.0
FIBRO    = 1.0
MYOFIB   = 1.0
COLLAGEN = 1.0
MMP      = 1.0
TIMP     = 1.0
FVC_st   = 80.0
DLCO_st  = 65.0

$ODE
// ── Pirfenidone PK ───────────────────────────────────────────────────────────
dxdt_DEPOT_P = -ka_P * DEPOT_P;
dxdt_CENT_P  =  ka_P * F_P * DEPOT_P
                - (CL_P + Q_P) / V1_P * CENT_P
                + Q_P / V2_P * PERI_P;
dxdt_PERI_P  =  Q_P / V1_P * CENT_P - Q_P / V2_P * PERI_P;

// ── Nintedanib PK ─────────────────────────────────────────────────────────────
dxdt_DEPOT_N = -ka_N * DEPOT_N;
dxdt_CENT_N  =  ka_N * F_N * DEPOT_N
                - (CL_N + Q_N) / V1_N * CENT_N
                + Q_N / V2_N * PERI_N;
dxdt_PERI_N  =  Q_N / V1_N * CENT_N - Q_N / V2_N * PERI_N;

// ── AEC-II dynamics ───────────────────────────────────────────────────────────
// Loss: ROS-driven damage, senescence by TGF-β feedback
// Gain: repair (EGF/HGF), proportional to remaining cells
double AEC2_damage_rate = k_AEC * ROS * AEC2;
double AEC2_repair_rate = k_rep * AEC2_ss * (1.0 - AEC2);
dxdt_AEC2 = AEC2_repair_rate - AEC2_damage_rate;

// ── TGF-β1 dynamics ──────────────────────────────────────────────────────────
// Production: M2 macrophages, damaged AEC, myofibroblasts
// Inhibition: pirfenidone
double TGFb_prod = kprod_TGFb
                   + kact_TGFb * (AEC2_ss - AEC2)  // more damage → more TGF-β
                   + 0.06 * M2
                   + 0.04 * MYOFIB;
double TGFb_deg  = kdeg_TGFb * TGFb;
dxdt_TGFb = TGFb_prod * (1.0 - inh_P_TGFb) - TGFb_deg;

// ── M2 macrophage dynamics ────────────────────────────────────────────────────
double M2_prod = kM2_act * TGFb;  // TGF-β promotes M2 polarization
dxdt_M2 = M2_prod * (1.0 - inh_P_M2) - kdeg_M2 * M2;

// ── ROS dynamics ─────────────────────────────────────────────────────────────
double ROS_prod = kprod_ROS
                  + kfb_ROS * (AEC2_ss - AEC2) * TGFb;
double ROS_clear = kdeg_ROS * ROS;
dxdt_ROS = ROS_prod * (1.0 - inh_P_ROS) - ROS_clear;

// ── Fibroblast activation ─────────────────────────────────────────────────────
// Activated by TGF-β, PDGF (proxied by nintedanib block)
double F_activ = kact_F * TGFb * F_ss;
double F_apop  = kapop_F * FIBRO;
dxdt_FIBRO = F_activ * (1.0 - inh_N_F) - F_apop;

// ── Myofibroblast differentiation ────────────────────────────────────────────
double M_diff   = kdiff_M * FIBRO * TGFb;
double M_apop   = 0.006 * MYOFIB;  // very low apoptosis in IPF (bcl-2 upregulated)
dxdt_MYOFIB = M_diff * (1.0 - inh_N_F * 0.5) - M_apop;

// ── ECM — Collagen ────────────────────────────────────────────────────────────
double Col_prod = kprod_Col * MYOFIB;
double Col_deg  = kdeg_Col  * MMP / (TIMP + 0.1) * COLLAGEN;
dxdt_COLLAGEN = Col_prod - Col_deg;

// ── MMP dynamics (MMP-1, -7, -9) ─────────────────────────────────────────────
double MMP_prod = kprod_MMP * FIBRO * 0.5 + kprod_MMP * M2 * 0.5;
double MMP_deg  = kdeg_MMP * MMP;
dxdt_MMP = MMP_prod - MMP_deg;

// ── TIMP dynamics (TGF-β upregulates TIMP-1) ─────────────────────────────────
double TIMP_prod = kprod_TIMP * TGFb * MYOFIB;
double TIMP_deg  = kdeg_TIMP * TIMP;
dxdt_TIMP = TIMP_prod - TIMP_deg;

// ── FVC (% predicted) — primary endpoint ─────────────────────────────────────
// Collagen excess above baseline drives FVC decline
double col_excess = (COLLAGEN > 1.0) ? (COLLAGEN - 1.0) : 0.0;
double fvc_loss = k_FVC_loss * col_excess * FVC_st;
dxdt_FVC_st = -fvc_loss;

// ── DLCO (% predicted) ────────────────────────────────────────────────────────
double dlco_loss = k_DLCO_loss * (AEC2_ss - AEC2 + 0.5 * col_excess) * DLCO_st;
dxdt_DLCO_st = -dlco_loss;

$TABLE
double Cp_pirf  = CENT_P / V1_P;                        // µg/mL
double Cn_nint  = (CENT_N / V1_N) / MW_N * 1e3;         // nM
double FVC_pct  = FVC_st;
double DLCO_pct = DLCO_st;
double FVC_decline_yr = k_FVC_loss * ((COLLAGEN > 1.0) ? (COLLAGEN - 1.0) : 0.0) * FVC_st * 8760.0;
double Col_norm = COLLAGEN;
double AEC2_norm = AEC2;
double TGFb_norm = TGFb;
double MMP_norm  = MMP;
double TIMP_norm = TIMP;
double MMP_TIMP_ratio = (TIMP > 0) ? MMP / TIMP : 1.0;
double Myofib_norm = MYOFIB;
double ROS_norm  = ROS;
// Periostin proxy (correlates with collagen deposition)
double Periostin_proxy = COLLAGEN * 1.2 * MYOFIB;
// MMP-7 serum proxy
double MMP7_proxy = MMP * TGFb * 0.8;
// KL-6 proxy
double KL6_proxy = (AEC2_ss - AEC2 + 0.5) * 800.0;  // U/mL scale

$CAPTURE
Cp_pirf Cn_nint FVC_pct DLCO_pct FVC_decline_yr Col_norm
AEC2_norm TGFb_norm MMP_norm TIMP_norm MMP_TIMP_ratio Myofib_norm ROS_norm
Periostin_proxy MMP7_proxy KL6_proxy
'

## ─── 2. COMPILE MODEL ───────────────────────────────────────────────────────

mod <- mcode("IPF_QSP", ipf_model_code)
cat("Model compiled successfully.\n")
cat("Compartments:", nrow(init(mod)), "\n")

## ─── 3. DOSING REGIMENS ─────────────────────────────────────────────────────

# Pirfenidone 801 mg TID (every 8h) × MW=185.22 g/mol
dose_pirf_mg  <- 801     # mg
dose_pirf_ug  <- dose_pirf_mg * 1e3  # µg

# Nintedanib 150 mg BID (every 12h) × MW=539.63 g/mol
dose_nint_mg  <- 150
dose_nint_ng  <- dose_nint_mg * 1e6  # ng

# Simulation duration (52 weeks = 364 days)
sim_duration <- 52 * 7 * 24  # hours
obs_times    <- seq(0, sim_duration, by=24)  # daily observations

## Treatment scenarios
mk_pirf_ev <- function() {
  ev(cmt=1, amt=dose_pirf_ug, ii=8, addl=sim_duration/8, time=0)
}
mk_nint_ev <- function() {
  ev(cmt=4, amt=dose_nint_ng, ii=12, addl=sim_duration/12, time=0)
}
mk_combo_ev <- function() {
  c(ev(cmt=1, amt=dose_pirf_ug, ii=8,  addl=sim_duration/8,  time=0),
    ev(cmt=4, amt=dose_nint_ng,  ii=12, addl=sim_duration/12, time=0))
}

## ─── 4. SIMULATION — 5 TREATMENT SCENARIOS ──────────────────────────────────

cat("\nRunning 5 treatment scenarios...\n")

scenarios <- list(
  "Placebo (Natural History)"         = ev(time=0, amt=0),
  "Pirfenidone 801 mg TID"            = mk_pirf_ev(),
  "Nintedanib 150 mg BID"             = mk_nint_ev(),
  "Combination (Pirf + Nint)"         = mk_combo_ev(),
  "Pirfenidone Low Dose (267 mg TID)" = ev(cmt=1, amt=267e3, ii=8, addl=sim_duration/8)
)

results <- lapply(names(scenarios), function(sc) {
  ev_obj <- scenarios[[sc]]
  out <- mod %>%
    mrgsim(events=ev_obj, end=sim_duration, delta=24, digits=4) %>%
    as.data.frame()
  out$Scenario <- sc
  out
}) %>% bind_rows()

# Convert time to weeks
results$Week <- results$time / (24 * 7)

cat("Simulation complete. Rows:", nrow(results), "\n")

## ─── 5. PK PROFILE (first 72h) ──────────────────────────────────────────────

pk_data <- mod %>%
  mrgsim(
    events = c(ev(cmt=1, amt=dose_pirf_ug, time=0),
               ev(cmt=4, amt=dose_nint_ng, time=0)),
    end = 72, delta = 0.5
  ) %>%
  as.data.frame()

## ─── 6. DOSE-RESPONSE ANALYSIS ──────────────────────────────────────────────

cat("\nDose-response analysis...\n")

pirf_doses <- c(267, 534, 801, 1068) * 1e3  # µg
nint_doses  <- c(50, 100, 150, 200) * 1e6   # ng

dr_pirf <- lapply(pirf_doses, function(d) {
  ev_d <- ev(cmt=1, amt=d, ii=8, addl=sim_duration/8)
  out  <- mod %>%
    mrgsim(events=ev_d, end=sim_duration, delta=24) %>%
    as.data.frame() %>%
    filter(time == max(time))
  data.frame(Drug="Pirfenidone", Dose_mg=d/1e3,
             FVC_final=out$FVC_pct, DLCO_final=out$DLCO_pct,
             Col_final=out$Col_norm)
}) %>% bind_rows()

dr_nint <- lapply(nint_doses, function(d) {
  ev_d <- ev(cmt=4, amt=d, ii=12, addl=sim_duration/12)
  out  <- mod %>%
    mrgsim(events=ev_d, end=sim_duration, delta=24) %>%
    as.data.frame() %>%
    filter(time == max(time))
  data.frame(Drug="Nintedanib", Dose_mg=d/1e6,
             FVC_final=out$FVC_pct, DLCO_final=out$DLCO_pct,
             Col_final=out$Col_norm)
}) %>% bind_rows()

dose_response <- bind_rows(dr_pirf, dr_nint)

## ─── 7. BIOMARKER TRAJECTORIES ──────────────────────────────────────────────

bm_data <- results %>%
  select(Week, Scenario, TGFb_norm, Myofib_norm, Col_norm,
         Periostin_proxy, MMP7_proxy, KL6_proxy, ROS_norm,
         MMP_TIMP_ratio, AEC2_norm)

## ─── 8. PLOT FUNCTIONS ──────────────────────────────────────────────────────

scenario_colors <- c(
  "Placebo (Natural History)"         = "#E74C3C",
  "Pirfenidone 801 mg TID"            = "#2471A3",
  "Nintedanib 150 mg BID"             = "#27AE60",
  "Combination (Pirf + Nint)"         = "#8E44AD",
  "Pirfenidone Low Dose (267 mg TID)" = "#E67E22"
)

## FVC over time
p_fvc <- ggplot(results, aes(x=Week, y=FVC_pct, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=scenario_colors) +
  labs(title="FVC % Predicted Over 52 Weeks (IPF QSP Model)",
       x="Time (weeks)", y="FVC (% predicted)",
       caption="ASCEND: pirfenidone −47.9% decline reduction; INPULSIS: nintedanib −50.1%") +
  theme_bw(base_size=13) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  geom_hline(yintercept=c(50, 70, 80), linetype="dashed", color="gray60", alpha=0.5)

## DLCO over time
p_dlco <- ggplot(results, aes(x=Week, y=DLCO_pct, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=scenario_colors) +
  labs(title="DLCO % Predicted Over 52 Weeks",
       x="Time (weeks)", y="DLCO (% predicted)") +
  theme_bw(base_size=13) +
  theme(legend.position="bottom", legend.title=element_blank())

## PK profiles
p_pk_pirf <- ggplot(pk_data, aes(x=time, y=Cp_pirf)) +
  geom_line(color="#2471A3", linewidth=1.5) +
  labs(title="Pirfenidone PK — Single Dose (801 mg)",
       x="Time (h)", y="Plasma Concentration (µg/mL)") +
  geom_hline(yintercept=30, linetype="dashed", color="red", label="EC50") +
  annotate("text", x=60, y=32, label="EC50 ~30 µg/mL", color="red", size=3.5) +
  theme_bw(base_size=13)

p_pk_nint <- ggplot(pk_data, aes(x=time, y=Cn_nint)) +
  geom_line(color="#27AE60", linewidth=1.5) +
  labs(title="Nintedanib PK — Single Dose (150 mg)",
       x="Time (h)", y="Plasma Concentration (nM)") +
  geom_hline(yintercept=20, linetype="dashed", color="red") +
  annotate("text", x=60, y=22, label="IC50 ~20 nM (FGFR1)", color="red", size=3.5) +
  theme_bw(base_size=13)

## Biomarker dynamics
p_tgfb <- ggplot(results, aes(x=Week, y=TGFb_norm, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=scenario_colors) +
  labs(title="TGF-β1 Level (Normalized)", x="Time (weeks)", y="TGF-β1 (normalized)") +
  theme_bw(base_size=13) + theme(legend.position="none")

p_col <- ggplot(results, aes(x=Week, y=Col_norm, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=scenario_colors) +
  labs(title="Collagen Deposition (Normalized)", x="Time (weeks)", y="Collagen") +
  theme_bw(base_size=13) + theme(legend.position="none")

p_mmp7 <- ggplot(results, aes(x=Week, y=MMP7_proxy, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=scenario_colors) +
  labs(title="Serum MMP-7 Proxy", x="Time (weeks)", y="MMP-7 (proxy)") +
  theme_bw(base_size=13) + theme(legend.position="none")

p_kl6 <- ggplot(results, aes(x=Week, y=KL6_proxy, color=Scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=scenario_colors) +
  labs(title="KL-6 Biomarker Proxy (U/mL)", x="Time (weeks)", y="KL-6 (U/mL equiv.)") +
  theme_bw(base_size=13) + theme(legend.position="none")

## Dose-response
p_dr <- ggplot(dose_response, aes(x=Dose_mg, y=FVC_final, color=Drug)) +
  geom_line(linewidth=1.4) + geom_point(size=3) +
  scale_color_manual(values=c("Pirfenidone"="#2471A3", "Nintedanib"="#27AE60")) +
  labs(title="Dose-Response: Final FVC at 52 Weeks",
       x="Dose (mg)", y="FVC % predicted at 52 weeks") +
  theme_bw(base_size=13)

## ─── 9. SUMMARY STATISTICS TABLE ────────────────────────────────────────────

fvc_summary <- results %>%
  filter(Week %in% c(0, 13, 26, 39, 52)) %>%
  group_by(Scenario, Week) %>%
  summarise(FVC=round(mean(FVC_pct), 2),
            DLCO=round(mean(DLCO_pct), 2),
            TGFb=round(mean(TGFb_norm), 3),
            Col=round(mean(Col_norm), 3),
            .groups="drop")

cat("\n=== FVC Summary at Key Timepoints ===\n")
print(as.data.frame(fvc_summary))

## ─── 10. CLINICAL TRIAL CALIBRATION CHECK ───────────────────────────────────

placebo_52  <- results %>% filter(Scenario=="Placebo (Natural History)", Week==52)
pirf_52     <- results %>% filter(Scenario=="Pirfenidone 801 mg TID", Week==52)
nint_52     <- results %>% filter(Scenario=="Nintedanib 150 mg BID", Week==52)
combo_52    <- results %>% filter(Scenario=="Combination (Pirf + Nint)", Week==52)

fvc_0       <- 80.0
pirf_0      <- results %>% filter(Scenario=="Pirfenidone 801 mg TID", Week==0)
placebo_decl <- fvc_0 - mean(placebo_52$FVC_pct)
pirf_decl    <- fvc_0 - mean(pirf_52$FVC_pct)
nint_decl    <- fvc_0 - mean(nint_52$FVC_pct)
combo_decl   <- fvc_0 - mean(combo_52$FVC_pct)

cat("\n=== Clinical Calibration (52-week FVC decline in % predicted) ===\n")
cat(sprintf("  Placebo:     −%.2f%%  (Literature target: ~2.5-3%% predicted/yr)\n", placebo_decl))
cat(sprintf("  Pirfenidone: −%.2f%%  (Target: ~47-50%% reduction vs placebo)\n", pirf_decl))
cat(sprintf("  Nintedanib:  −%.2f%%  (Target: ~50%% reduction vs placebo)\n", nint_decl))
cat(sprintf("  Combination: −%.2f%%  (Target: ≥50%% reduction vs placebo)\n", combo_decl))
if (placebo_decl > 0) {
  cat(sprintf("  Pirf reduction vs placebo: %.1f%%\n", (1 - pirf_decl/placebo_decl)*100))
  cat(sprintf("  Nint reduction vs placebo: %.1f%%\n", (1 - nint_decl/placebo_decl)*100))
}

## ─── 11. SAVE PLOTS ──────────────────────────────────────────────────────────

cat("\nPlots generated (use gridExtra or patchwork to display):\n")
cat("  p_fvc, p_dlco, p_pk_pirf, p_pk_nint\n")
cat("  p_tgfb, p_col, p_mmp7, p_kl6, p_dr\n")

## Combined plot (4×2 grid)
if (requireNamespace("gridExtra", quietly=TRUE)) {
  grid_plot <- gridExtra::grid.arrange(
    p_fvc, p_dlco, p_pk_pirf, p_pk_nint,
    p_tgfb, p_col, p_dr, p_kl6,
    ncol=2, nrow=4
  )
}

cat("\nIPF QSP model simulation complete.\n")
cat("Results in 'results' data frame, summaries in 'fvc_summary'.\n")
