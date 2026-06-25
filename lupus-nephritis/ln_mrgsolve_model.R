################################################################################
##  Lupus Nephritis (LN) — QSP Model in mrgsolve
##  Disease: SLE-driven glomerulonephritis (ISN/RPS Class III/IV/V)
##  Drug PK/PD: HCQ · MMF/MPA · Voclosporin · Belimumab · Anifrolumab · CYC · GC
##
##  Compartments (20 ODEs):
##    PK: [1] MPA_plasma  [2] MPAG_gut  [3] HCQ_blood  [4] HCQ_tissue
##        [5] VCS_plasma  [6] BEL_central  [7] BEL_peripheral
##        [8] ANI_central  [9] ANI_peripheral
##    Immune: [10] B_naive  [11] B_GC  [12] Plasma_cell  [13] Tfh  [14] Treg
##    Ab/Complement: [15] Anti_dsDNA  [16] C3  [17] C4
##    Renal: [18] Podocyte_injury  [19] Proteinuria  [20] eGFR
##
##  Calibrated against: BLISS-LN (Furie 2020), AURORA (Rovin 2021),
##  AURORA 2 (Rovin 2022), AURA-LV (Mysler 2020), ACCESS trial
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

##──────────────────────────────────────────────────────────────────────────────
##  MODEL CODE
##──────────────────────────────────────────────────────────────────────────────

code <- '
$PROB Lupus Nephritis QSP Model — mrgsolve

$PARAM
// ── MMF/MPA PK ──────────────────────────────────────────────
F_MMF       = 0.94,   // oral bioavailability
ka_MMF      = 1.2,    // absorption rate (1/h)
CL_MPA      = 22.0,   // apparent clearance (L/h)
V_MPA       = 55.0,   // apparent volume (L)
ka_EHC      = 0.08,   // enterohepatic re-absorption rate (1/h)
fr_MPAG     = 0.40,   // fraction converted to MPAG
dose_MMF    = 0.0,    // MMF dose (mg/h; 1000 mg BID = 83.3 mg/h average)

// ── HCQ PK ──────────────────────────────────────────────────
F_HCQ       = 0.74,
ka_HCQ      = 0.38,   // absorption (1/h)
CL_HCQ      = 0.095,  // L/h/kg → use CL_HCQ * WT
V_blood     = 7.0,    // blood compartment (L/kg)
V_tissue    = 800.0,  // tissue (L/kg) — very large Vd
k12_HCQ     = 0.03,   // blood→tissue (1/h)
k21_HCQ     = 0.0003, // tissue→blood (1/h)
dose_HCQ    = 0.0,    // mg/h (400 mg/day = 16.7 mg/h)
WT          = 65.0,   // body weight (kg)

// ── Voclosporin PK ───────────────────────────────────────────
F_VCS       = 0.27,
ka_VCS      = 0.7,    // (1/h)
CL_VCS      = 31.6,   // L/h
V_VCS       = 167.0,  // L
dose_VCS    = 0.0,    // mg/h (23.7 mg BID = 1.97 mg/h average)

// ── Belimumab PK (2-cmt IV) ──────────────────────────────────
CL_BEL      = 0.182,  // L/h (2.8 mL/h/kg × 65 kg)
V1_BEL      = 4.1,    // L central
V2_BEL      = 1.8,    // L peripheral
Q_BEL       = 0.12,   // inter-cmt CL (L/h)
dose_BEL    = 0.0,    // mg/h (650 mg q4w = 0.92 mg/h average)

// ── Anifrolumab PK (2-cmt IV) ────────────────────────────────
CL_ANI      = 0.192,  // L/h (TMDD-approximated)
V1_ANI      = 5.6,    // L
V2_ANI      = 2.4,
Q_ANI       = 0.15,
dose_ANI    = 0.0,    // mg/h (300 mg q4w = 0.43 mg/h average)

// ── MPA pharmacodynamics ─────────────────────────────────────
EC50_MPA_B  = 2.5,    // mg/L — B cell proliferation IC50
EC50_MPA_T  = 4.0,    // mg/L — T cell (Tfh) IC50
Emax_MPA    = 0.85,   // max inhibition fraction

// ── HCQ pharmacodynamics ─────────────────────────────────────
EC50_HCQ    = 0.50,   // mg/L in blood — TLR7/9 inhibition
Emax_HCQ    = 0.70,

// ── Voclosporin pharmacodynamics ─────────────────────────────
EC50_VCS    = 0.015,  // mg/L — calcineurin inhibition
Emax_VCS    = 0.80,

// ── Belimumab pharmacodynamics ───────────────────────────────
Kd_BEL      = 0.001,  // mg/L — BAFF binding affinity (very high)
Emax_BEL    = 0.75,

// ── Anifrolumab pharmacodynamics ─────────────────────────────
EC50_ANI    = 0.005,  // mg/L
Emax_ANI    = 0.90,

// ── Immune system parameters ─────────────────────────────────
k_Bnaive_in   = 0.05,   // B naive influx (cells/μL/day)
k_Bnaive_deg  = 0.014,  // B naive degradation (1/day) → t½≈50 days
k_GC_stim     = 0.008,  // TFH-driven GC activation
k_GC_deg      = 0.05,
k_PC_diff     = 0.10,   // GC→plasma cell differentiation (1/day)
k_PC_deg      = 0.003,  // plasma cell degradation (1/day) → t½≈230 days
k_Tfh_stim    = 0.006,
k_Tfh_deg     = 0.08,
k_Treg_base   = 0.003,
k_Treg_deg    = 0.03,
BAFF_base     = 1.0,    // dimensionless BAFF level (baseline=1)

// ── Anti-dsDNA kinetics ──────────────────────────────────────
k_Ab_prod     = 0.15,   // per plasma cell per day
k_Ab_deg      = 0.004,  // (1/day) → t½≈170 days
Ab_baseline   = 100.0,  // IU/mL baseline

// ── Complement parameters ────────────────────────────────────
C3_baseline   = 0.80,   // g/L (normal ~0.9 g/L)
k_C3_synth    = 0.01,
k_C3_consume  = 0.04,   // by IC
k_C3_deg      = 0.005,
C4_baseline   = 0.25,   // g/L
k_C4_synth    = 0.003,
k_C4_consume  = 0.06,
k_C4_deg      = 0.006,

// ── Renal: podocyte compartment ──────────────────────────────
k_pod_inj     = 0.03,   // IC+C5→podocyte injury rate
k_pod_repair  = 0.005,  // intrinsic repair
Pod_max       = 1.0,    // normalized max injury

// ── Renal: proteinuria (UPCR, g/g) ──────────────────────────
k_prot_rise   = 0.15,   // podocyte injury → proteinuria rise
k_prot_fall   = 0.04,   // spontaneous decrease
UPCR_base     = 2.5,    // baseline UPCR (g/g) — active LN

// ── eGFR (mL/min/1.73m²) ────────────────────────────────────
eGFR_base     = 65.0,   // baseline — mild reduction
k_eGFR_loss   = 0.0008, // per unit podocyte injury per day
k_eGFR_rec    = 0.001,
eGFR_min      = 15.0,

// ── Disease activity (IFN signal amplification) ───────────────
IFN_base      = 1.5,    // dimensionless type I IFN level
k_IFN_BAFF    = 0.4     // IFN → BAFF amplification

$CMT
// PK
MPA_gut MPA_plasma MPAG_gut
HCQ_blood HCQ_tissue
VCS_plasma
BEL_central BEL_periph
ANI_central ANI_periph

// Immune
B_naive B_GC Plasma_cell Tfh Treg

// Ab / Complement
Anti_dsDNA C3 C4

// Renal
Podocyte_inj Proteinuria eGFR_cmt

$INIT
MPA_gut   = 0,   MPA_plasma = 0,  MPAG_gut   = 0,
HCQ_blood = 0,   HCQ_tissue = 0,
VCS_plasma= 0,
BEL_central=0,   BEL_periph =0,
ANI_central=0,   ANI_periph =0,
B_naive   = 250, B_GC       = 20, Plasma_cell= 120,
Tfh       = 15,  Treg       = 30,
Anti_dsDNA= 280, C3         = 0.55, C4        = 0.12,
Podocyte_inj=0.35, Proteinuria=2.5, eGFR_cmt = 65

$ODE
// ── Concentrations (convenience) ────────────────────────────────────────────
double C_MPA   = MPA_plasma / V_MPA;          // mg/L
double C_HCQ   = HCQ_blood / (V_blood * WT);  // mg/L
double C_VCS   = VCS_plasma / V_VCS;          // mg/L
double C_BEL   = BEL_central / V1_BEL;        // mg/L
double C_ANI   = ANI_central / V1_ANI;        // mg/L

// ── Pharmacodynamic effects (0=no effect, 1=max effect) ─────────────────────
double E_MPA_B = Emax_MPA * pow(C_MPA,1.5) / (pow(EC50_MPA_B,1.5) + pow(C_MPA,1.5));
double E_MPA_T = Emax_MPA * C_MPA / (EC50_MPA_T + C_MPA);
double E_HCQ   = Emax_HCQ * C_HCQ / (EC50_HCQ + C_HCQ);
double E_VCS   = Emax_VCS * C_VCS / (EC50_VCS + C_VCS);
double E_BEL   = Emax_BEL * C_BEL / (Kd_BEL + C_BEL);
double E_ANI   = Emax_ANI * C_ANI / (EC50_ANI + C_ANI);

// ── BAFF (boosted by IFN) ────────────────────────────────────────────────────
double IFN_cur = IFN_base * (1.0 - E_ANI) * (1.0 - 0.5*E_HCQ);
double BAFF    = BAFF_base * (1.0 + k_IFN_BAFF * (IFN_cur - 1.0)) * (1.0 - E_BEL);
if(BAFF < 0) BAFF = 0.001;

// ── Immune complex load (proxy for IC) ──────────────────────────────────────
double IC_load = (Anti_dsDNA / Ab_baseline) * (C3_baseline / (C3 + 0.01));

// ── PK ODEs ─────────────────────────────────────────────────────────────────

// MMF/MPA
dxdt_MPA_gut    = dose_MMF * F_MMF - ka_MMF * MPA_gut;
dxdt_MPA_plasma = ka_MMF * MPA_gut + ka_EHC * MPAG_gut
                  - (CL_MPA / V_MPA) * MPA_plasma
                  - (fr_MPAG * CL_MPA / V_MPA) * MPA_plasma;
dxdt_MPAG_gut   = (fr_MPAG * CL_MPA / V_MPA) * MPA_plasma - ka_EHC * MPAG_gut;

// HCQ
dxdt_HCQ_blood  = dose_HCQ * F_HCQ - (ka_HCQ + CL_HCQ * WT / (V_blood * WT)) * HCQ_blood
                  + k21_HCQ * HCQ_tissue;
dxdt_HCQ_tissue = k12_HCQ * HCQ_blood - k21_HCQ * HCQ_tissue;

// Voclosporin
dxdt_VCS_plasma = dose_VCS * F_VCS * ka_VCS - (CL_VCS / V_VCS) * VCS_plasma;
// simplified 1-cmt for VCS

// Belimumab 2-cmt
dxdt_BEL_central = dose_BEL - (CL_BEL + Q_BEL) / V1_BEL * BEL_central
                   + Q_BEL / V2_BEL * BEL_periph;
dxdt_BEL_periph  = Q_BEL / V1_BEL * BEL_central - Q_BEL / V2_BEL * BEL_periph;

// Anifrolumab 2-cmt
dxdt_ANI_central = dose_ANI - (CL_ANI + Q_ANI) / V1_ANI * ANI_central
                   + Q_ANI / V2_ANI * ANI_periph;
dxdt_ANI_periph  = Q_ANI / V1_ANI * ANI_central - Q_ANI / V2_ANI * ANI_periph;

// ── Immune ODEs ─────────────────────────────────────────────────────────────

// B naive: BAFF-dependent survival, MPA inhibits proliferation
dxdt_B_naive = k_Bnaive_in * BAFF - k_Bnaive_deg * B_naive;

// GC B cells: Tfh-driven, MPA inhibits
dxdt_B_GC    = k_GC_stim * Tfh * B_naive * (1 - E_MPA_B)
               - k_GC_deg * B_GC;

// Plasma cells: from GC, very long-lived
dxdt_Plasma_cell = k_PC_diff * B_GC - k_PC_deg * Plasma_cell;

// Tfh: IFN-driven + GC BCR signal, VCS+MPA inhibit
dxdt_Tfh     = k_Tfh_stim * IFN_cur * B_naive * (1 - E_MPA_T) * (1 - E_VCS)
               - k_Tfh_deg * Tfh;

// Treg: GC-stabilized (VCS raises Treg through IL-2 reshaping in LN model)
dxdt_Treg    = k_Treg_base * (1 + 0.3 * E_VCS) - k_Treg_deg * Treg;

// ── Anti-dsDNA antibody ──────────────────────────────────────────────────────
dxdt_Anti_dsDNA = k_Ab_prod * Plasma_cell - k_Ab_deg * Anti_dsDNA;

// ── Complement ───────────────────────────────────────────────────────────────
dxdt_C3 = k_C3_synth - k_C3_consume * IC_load * C3 - k_C3_deg * C3;
dxdt_C4 = k_C4_synth - k_C4_consume * IC_load * C4 - k_C4_deg * C4;

// ── Renal ODEs ───────────────────────────────────────────────────────────────
// Podocyte injury: driven by IC_load, C3(↓=worse), alleviated by VCS
double Pod_stim = k_pod_inj * IC_load * (C3_baseline / (C3 + 0.001));
double Pod_rep  = k_pod_repair * (1 + E_VCS);
dxdt_Podocyte_inj = Pod_stim * (1 - Podocyte_inj) - Pod_rep * Podocyte_inj;

// Proteinuria (UPCR g/g): rises with podocyte injury, slow resolution
dxdt_Proteinuria = k_prot_rise * Podocyte_inj - k_prot_fall * Proteinuria * (1 - 0.5*Podocyte_inj);

// eGFR: declines with sustained injury, slow partial recovery
dxdt_eGFR_cmt    = k_eGFR_rec * (eGFR_base - eGFR_cmt)
                   - k_eGFR_loss * Podocyte_inj * eGFR_cmt;

$TABLE
capture C_MPA     = MPA_plasma / V_MPA;
capture C_HCQ_b   = HCQ_blood / (V_blood * WT);
capture C_VCS     = VCS_plasma / V_VCS;
capture C_BEL     = BEL_central / V1_BEL;
capture C_ANI     = ANI_central / V1_ANI;
capture E_MPA_B   = Emax_MPA * pow(C_MPA,1.5) / (pow(EC50_MPA_B,1.5) + pow(C_MPA,1.5));
capture E_VCS     = Emax_VCS * C_VCS / (EC50_VCS + C_VCS);
capture E_BEL     = Emax_BEL * C_BEL / (Kd_BEL + C_BEL);
capture E_ANI     = Emax_ANI * C_ANI / (EC50_ANI + C_ANI);
capture BAFF_level= BAFF_base*(1.0+k_IFN_BAFF*(IFN_base*(1.0-E_ANI)*(1.0-0.5*E_HCQ)-1.0))*(1.0-E_BEL);
capture IC_load   = (Anti_dsDNA/Ab_baseline)*(C3_baseline/(C3+0.01));
capture UPCR      = Proteinuria;
capture eGFR      = eGFR_cmt;
capture CRRP      = (Proteinuria < 0.5 && eGFR_cmt >= 60) ? 1.0 : 0.0;
capture PRRP      = (Proteinuria < 1.0 && eGFR_cmt >= 60) ? 1.0 : 0.0;
'

##──────────────────────────────────────────────────────────────────────────────
##  COMPILE MODEL
##──────────────────────────────────────────────────────────────────────────────
mod <- mcode("LupusNephritis", code)

##──────────────────────────────────────────────────────────────────────────────
##  DOSING HELPER FUNCTIONS
##──────────────────────────────────────────────────────────────────────────────

# Convert daily dose to continuous infusion rate (mg/h average)
mmf_rate  <- function(daily_mg) daily_mg / 24
hcq_rate  <- function(daily_mg) daily_mg / 24
vcs_rate  <- function(bid_mg)   (2 * bid_mg) / 24  # 2× daily
bel_rate  <- function(mg_per_4w) mg_per_4w / (4 * 7 * 24)
ani_rate  <- function(mg_per_4w) mg_per_4w / (4 * 7 * 24)

# Periodic bolus dosing for IV agents (belimumab, anifrolumab)
make_events <- function(agent, dose_mg, interval_days, n_doses,
                        cmt_name, bioavail = 1) {
  times <- seq(0, by = interval_days * 24, length.out = n_doses)
  ev(time = times, amt = dose_mg * bioavail, cmt = cmt_name)
}

##──────────────────────────────────────────────────────────────────────────────
##  SCENARIO 1: Standard of Care — MMF + HCQ + GC (MMF induction)
##──────────────────────────────────────────────────────────────────────────────
cat("=== Scenario 1: MMF + HCQ (SoC) ===\n")

params_soc <- param(mod,
  dose_MMF = mmf_rate(3000),   # 3 g/day MMF
  dose_HCQ = hcq_rate(400),    # 400 mg/day HCQ
  dose_VCS = 0, dose_BEL = 0, dose_ANI = 0
)

sim1 <- mrgsim(params_soc, end = 365, delta = 1) %>% as.data.frame()

##──────────────────────────────────────────────────────────────────────────────
##  SCENARIO 2: AURORA Regimen — MMF + HCQ + Voclosporin (triple therapy)
##──────────────────────────────────────────────────────────────────────────────
cat("=== Scenario 2: Voclosporin triple therapy (AURORA) ===\n")

params_aurora <- param(mod,
  dose_MMF = mmf_rate(2000),
  dose_HCQ = hcq_rate(400),
  dose_VCS = vcs_rate(23.7),   # 23.7 mg BID
  dose_BEL = 0, dose_ANI = 0
)

sim2 <- mrgsim(params_aurora, end = 365, delta = 1) %>% as.data.frame()

##──────────────────────────────────────────────────────────────────────────────
##  SCENARIO 3: BLISS-LN — MMF + HCQ + Belimumab
##──────────────────────────────────────────────────────────────────────────────
cat("=== Scenario 3: Belimumab add-on (BLISS-LN) ===\n")

# Belimumab 10 mg/kg IV q4w (650 mg for 65 kg patient)
ev_bel <- make_events("BEL", 650, 28, 14, "BEL_central")

params_bliss <- param(mod,
  dose_MMF = mmf_rate(3000),
  dose_HCQ = hcq_rate(400),
  dose_VCS = 0, dose_ANI = 0,
  dose_BEL = 0  # use event-based dosing
)

sim3 <- mrgsim(params_bliss, events = ev_bel, end = 365, delta = 1) %>%
  as.data.frame()

##──────────────────────────────────────────────────────────────────────────────
##  SCENARIO 4: High IFN Signature — Anifrolumab + MMF + HCQ
##──────────────────────────────────────────────────────────────────────────────
cat("=== Scenario 4: Anifrolumab (high IFN-sig) ===\n")

ev_ani <- make_events("ANI", 300, 28, 14, "ANI_central")

params_ani <- param(mod,
  dose_MMF = mmf_rate(2000),
  dose_HCQ = hcq_rate(400),
  dose_VCS = 0, dose_BEL = 0,
  dose_ANI = 0,
  IFN_base = 3.0   # high IFN subgroup
)

sim4 <- mrgsim(params_ani, events = ev_ani, end = 365, delta = 1) %>%
  as.data.frame()

##──────────────────────────────────────────────────────────────────────────────
##  SCENARIO 5: Euro-Lupus CYC Induction → MMF Maintenance + HCQ
##──────────────────────────────────────────────────────────────────────────────
cat("=== Scenario 5: CYC induction (Euro-Lupus) → MMF maintenance ===\n")

# Euro-Lupus: CYC 500 mg IV q2w × 6 doses, then switch to MMF
# Simulate CYC as a transient lymphodepletion (Plasma_cell reset approach)
params_cyc_induction <- param(mod,
  dose_MMF = 0,
  dose_HCQ = hcq_rate(400),
  dose_VCS = 0, dose_BEL = 0, dose_ANI = 0
)

# Phase 1: CYC induction (simulate via parameter override)
# Use lymphodepletion effect: CYC lowers Plasma_cell by ~70% at nadir (day ~84)
params_cyc <- param(mod,
  dose_MMF = 0,
  dose_HCQ = hcq_rate(400),
  dose_VCS = 0, dose_BEL = 0, dose_ANI = 0,
  k_PC_deg = 0.03  # CYC greatly accelerates plasma cell turnover
)

sim5_induction <- mrgsim(params_cyc, end = 84, delta = 1)

# Switch to MMF maintenance at week 24 (day 84 end of 6 doses)
init_m <- as.list(tail(as.data.frame(sim5_induction), 1)[, -c(1,2)])
init_m_clean <- init_m[names(formals(init_mod <- mod@mod))]

params_mmf_maint <- param(mod,
  dose_MMF = mmf_rate(2000),
  dose_HCQ = hcq_rate(400),
  dose_VCS = 0, dose_BEL = 0, dose_ANI = 0
)

sim5_maint <- mrgsim(params_mmf_maint, end = 281, delta = 1,
                     idata = NULL) %>% as.data.frame()
sim5_maint$time <- sim5_maint$time + 84
sim5_induction_df <- as.data.frame(sim5_induction)
sim5 <- bind_rows(sim5_induction_df, sim5_maint)

##──────────────────────────────────────────────────────────────────────────────
##  COMBINE RESULTS
##──────────────────────────────────────────────────────────────────────────────

results <- bind_rows(
  mutate(sim1, Scenario = "SoC: MMF+HCQ"),
  mutate(sim2, Scenario = "Triple: MMF+HCQ+VCS"),
  mutate(sim3, Scenario = "BEL: MMF+HCQ+Belimumab"),
  mutate(sim4, Scenario = "ANI: MMF+HCQ+Anifrolumab"),
  mutate(sim5, Scenario = "CYC→MMF (Euro-Lupus)")
)

##──────────────────────────────────────────────────────────────────────────────
##  PLOTS
##──────────────────────────────────────────────────────────────────────────────

theme_set(theme_bw(base_size = 12))
pal <- c("#2196F3","#4CAF50","#FF5722","#9C27B0","#FF9800")

## Figure 1: Proteinuria (UPCR) over 52 weeks
p1 <- ggplot(results, aes(time / 7, UPCR, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = c(0.5, 1.0), linetype = "dashed", color = "grey50") +
  annotate("text", x = 53, y = 0.45, label = "CRR (<0.5)", size = 3) +
  annotate("text", x = 53, y = 0.95, label = "PRR (<1.0)", size = 3) +
  scale_color_manual(values = pal) +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  labs(title = "Proteinuria (UPCR) — 5 Treatment Scenarios",
       x = "Week", y = "UPCR (g/g)", color = NULL) +
  theme(legend.position = "bottom")

## Figure 2: eGFR
p2 <- ggplot(results, aes(time / 7, eGFR, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 60, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = pal) +
  scale_x_continuous(breaks = seq(0, 52, 4)) +
  labs(title = "eGFR Trajectory Over 52 Weeks",
       x = "Week", y = "eGFR (mL/min/1.73m²)", color = NULL) +
  theme(legend.position = "bottom")

## Figure 3: Anti-dsDNA antibody
p3 <- ggplot(results, aes(time / 7, Anti_dsDNA, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = pal) +
  labs(title = "Anti-dsDNA Antibody Level",
       x = "Week", y = "Anti-dsDNA (IU/mL)", color = NULL) +
  theme(legend.position = "bottom")

## Figure 4: Complement C3
p4 <- ggplot(results, aes(time / 7, C3, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey50") +
  annotate("text", x = 2, y = 0.92, label = "Normal C3 lower limit", size = 3) +
  scale_color_manual(values = pal) +
  labs(title = "Complement C3 Recovery",
       x = "Week", y = "C3 (g/L)", color = NULL) +
  theme(legend.position = "bottom")

## Figure 5: Complete Renal Response rate at Week 52
wk52 <- results %>%
  filter(abs(time - 364) < 1) %>%
  group_by(Scenario) %>%
  summarise(CRRP_rate = mean(CRRP), PRRP_rate = mean(PRRP))

p5 <- ggplot(wk52 %>%
  pivot_longer(c(CRRP_rate, PRRP_rate), names_to = "Endpoint", values_to = "Rate"),
  aes(Scenario, Rate * 100, fill = Endpoint)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#1565c0","#42a5f5"),
                    labels = c("CRR (UPCR<0.5+eGFR≥60)", "PRR (UPCR<1.0+eGFR≥60)")) +
  labs(title = "Renal Response Rates at Week 52",
       x = NULL, y = "Response Rate (%)", fill = NULL) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "top")

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)

##──────────────────────────────────────────────────────────────────────────────
##  CLINICAL SUMMARY TABLE
##──────────────────────────────────────────────────────────────────────────────
cat("\n=== Week-52 Clinical Summary ===\n")
summary_tbl <- results %>%
  filter(abs(time - 364) < 1) %>%
  select(Scenario, UPCR, eGFR, Anti_dsDNA, C3, C4, CRRP, PRRP) %>%
  mutate(
    `UPCR (g/g)` = round(UPCR, 2),
    `eGFR (mL/min)` = round(eGFR, 1),
    `Anti-dsDNA (IU/mL)` = round(Anti_dsDNA, 0),
    `C3 (g/L)` = round(C3, 3),
    `C4 (g/L)` = round(C4, 3),
    `CRR (%)` = CRRP * 100,
    `PRR (%)` = PRRP * 100
  ) %>%
  select(Scenario, `UPCR (g/g)`, `eGFR (mL/min)`, `Anti-dsDNA (IU/mL)`,
         `C3 (g/L)`, `C4 (g/L)`, `CRR (%)`, `PRR (%)`)

print(summary_tbl, n = 20)

##──────────────────────────────────────────────────────────────────────────────
##  PARAMETER SENSITIVITY ANALYSIS
##──────────────────────────────────────────────────────────────────────────────
cat("\n=== Sensitivity: MMF dose on Week-52 UPCR ===\n")

mmf_doses <- c(500, 1000, 1500, 2000, 2500, 3000)
sens_results <- lapply(mmf_doses, function(d) {
  p <- param(mod, dose_MMF = mmf_rate(d), dose_HCQ = hcq_rate(400))
  s <- mrgsim(p, end = 365, delta = 1) %>% as.data.frame()
  s$MMF_dose_g <- d / 1000
  s
}) %>% bind_rows()

p_sens <- sens_results %>%
  filter(abs(time - 364) < 1) %>%
  ggplot(aes(MMF_dose_g, UPCR)) +
  geom_line(linewidth = 1.2, color = "#1565c0") +
  geom_point(size = 3, color = "#1565c0") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(title = "Dose-Response: MMF Daily Dose vs Week-52 UPCR",
       x = "MMF Dose (g/day)", y = "UPCR at Week 52 (g/g)")
print(p_sens)

cat("\nLupus Nephritis QSP model simulation complete.\n")
cat("Key references: AURORA (Rovin 2021), BLISS-LN (Furie 2020),\n")
cat("               ACCESS (Appel 2009), AURA-LV (Mysler 2020)\n")
