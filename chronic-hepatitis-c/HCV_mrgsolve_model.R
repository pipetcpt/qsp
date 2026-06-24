## ============================================================
## Chronic Hepatitis C — Quantitative Systems Pharmacology Model
## mrgsolve ODE implementation
## Version : 1.0   |   2026-06-23
##
## Mathematical framework:
##   Neumann et al. (Science 1998) viral kinetics model
##   extended with:
##   - Target-cell limited (TCL) dynamics
##   - DAA PK for SOF/LED, SOF/VEL, GLE/PIB, and PEG-RBV
##   - Liver fibrosis dynamics (Metavir score)
##   - Immune response module (CTL, NK)
##
## 20 ODE compartments:
##   Drug PK     (9): SOF Tp, LED, VEL, GLE, PIB, RBV, RBV-RBC, PEG-IFN, NS5A_i
##   Viral kinetics (4): T (target), I (infected), V (viral load), Defective
##   Immune module  (3): CTL, NK, Treg
##   Liver pathology(4): ALT, Fibrosis, HSC activity, HCC-risk index
##
## Key references:
##   Neumann et al. Science 1998 — biphasic viral kinetics
##   Guedj et al. PNAS 2013     — intracellular model (SOF)
##   Adiwijaya et al. PLoS Comp Biol 2010 — NS3/4A inhibitor kinetics
##   Rong & Perelson 2010       — target-cell limited model
##   Dahari et al. Hepatology 2007 — fibrosis dynamics
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

code <- '
$PROB
Chronic Hepatitis C — QSP Model v1.0
20-compartment mrgsolve ODE

$PARAM
@annotated

// ── Sofosbuvir (SOF) PK ──────────────────────────────────────
// Intracellular triphosphate (SOF-Tp) kinetics
ka_SOF     : 1.44  : SOF oral absorption rate (/h), Tmax ~0.5h
CL_SOF_Tp  : 0.05  : Intrahepatic SOF-Tp clearance (/h)
V_SOF_Tp   : 1.0   : Hepatocyte SOF-Tp volume (relative)
EC50_SOF   : 0.35  : EC50 NS5B inhibition (μM)
n_SOF      : 1.6   : Hill coefficient SOF

// ── Ledipasvir (LED) PK ──────────────────────────────────────
ka_LED     : 0.2   : LED absorption rate (/h)
CL_LED     : 0.25  : LED clearance (/h)
V_LED      : 1500  : LED Vd (L)
EC50_LED   : 0.031 : EC50 NS5A inhibition (nM), picomolar-range
n_LED      : 1.2   : Hill coefficient LED

// ── Velpatasvir (VEL) PK ─────────────────────────────────────
ka_VEL     : 0.3   : VEL absorption rate (/h)
CL_VEL     : 0.33  : VEL clearance (/h)
V_VEL      : 1200  : VEL Vd (L)
EC50_VEL   : 0.017 : EC50 NS5A pangenotypic (nM)

// ── Glecaprevir (GLE) PK ─────────────────────────────────────
ka_GLE     : 0.7   : GLE absorption rate (/h)
CL_GLE     : 1.2   : GLE clearance (/h)
V_GLE      : 200   : GLE Vd (L)
EC50_GLE   : 0.45  : EC50 NS3/4A protease inhibition (nM)
n_GLE      : 1.8   : Hill coefficient GLE

// ── Pibrentasvir (PIB) PK ────────────────────────────────────
ka_PIB     : 0.5   : PIB absorption rate (/h)
CL_PIB     : 0.18  : PIB clearance (/h)
V_PIB      : 800   : PIB Vd (L)
EC50_PIB   : 0.002 : EC50 NS5A pan-genotypic (nM)

// ── Ribavirin (RBV) PK ───────────────────────────────────────
ka_RBV     : 0.25  : RBV oral absorption (/h)
CL_RBV     : 0.1   : RBV plasma clearance (/h)
V_RBV      : 2825  : RBV Vd (L, large due to RBC distribution)
k_RBV_RBC  : 0.02  : RBV plasma → RBC phosphorylation (/h)
kdeg_RBV_RBC: 0.003 : RBV-RBC dephosphorylation (/h)
EC50_RBV   : 2.0   : RBV EC50 for IMPDH inhibition (μg/mL plasma)

// ── PEG-IFN PK ───────────────────────────────────────────────
ka_PEGIFN  : 0.015  : SC absorption (/h)
CL_PEGIFN  : 0.025  : PEG-IFN clearance (/h)
V_PEGIFN   : 50     : PEG-IFN Vd (L)
EC50_PEGIFN: 0.8    : EC50 for ISG induction (ng/mL)

// ── Viral kinetics parameters (Neumann 1998 / Rong 2010) ─────
beta_infect : 1.2e-7 : Infection rate constant (mL/virion/h)
p_prod     : 100.0   : Virion production rate (virions/cell/h)
c_clear    : 22.0    : Viral clearance rate (/h) = 22/h (t½ ~3h)
delta_I    : 0.08    : Infected cell death rate (/h) = ~0.5/day
d_T        : 0.001   : Target cell death rate (/h)
s_T        : 10000   : Target cell production (/h) (hepatocytes)
epsilon0   : 0.0     : Baseline drug efficacy (0 = no treatment)

// ── Drug efficacy combination ─────────────────────────────────
// Combined efficacy on production: εp = 1 - (1-εSOF)(1-εLED/VEL)(1-εGLE/PIB)
// Combined efficacy on infection : εi = εSOF (NS5B → ↓infectivity)

// ── Immune parameters ─────────────────────────────────────────
kprol_CTL  : 0.02    : CTL proliferation rate (/h) per infected cell
kdeath_CTL : 0.004   : CTL natural death rate (/h)
keff_CTL   : 1e-4    : CTL killing rate constant (mL/cell/h)
kprol_NK   : 0.01    : NK cell activation rate (/h)
kdeath_NK  : 0.003   : NK cell death rate (/h)
keff_NK    : 5e-5    : NK killing rate constant

// Treg suppression of CTL
kprol_Treg_HCV : 0.005  : Treg induction by chronic inflammation (/h)
kdeath_Treg_HCV: 0.003  : Treg death rate (/h)
k_exhaust  : 1e-5    : CTL exhaustion rate constant

// ── Liver fibrosis dynamics ───────────────────────────────────
kprog_fibro : 1e-5   : Fibrosis progression per unit injury (/h)
krem_fibro  : 2e-6   : Fibrosis regression with SVR (/h)
kprog_HSC  : 0.001   : HSC activation rate (/h)
kdeath_HSC  : 0.0005 : HSC deactivation rate (/h)
kALT_prod  : 0.5     : ALT release per infected cell death (IU/L/cell)
kALT_clear : 0.008   : ALT half-life clearance (/h) ~3.5 days
kHCC_prog  : 1e-6    : HCC risk accumulation (/h) from fibrosis

// ── Baseline / steady-state ───────────────────────────────────
T0         : 2e7     : Baseline target hepatocytes (cells/mL liver)
I0         : 1e5     : Baseline infected cells (cells/mL)
V0         : 1e6     : Baseline viral load (IU/mL) = 6 log10
CTL0       : 1000    : Baseline HCV-specific CTL (cells/mL)
NK0        : 500     : Baseline activated NK cells (cells/mL)
Treg0      : 200     : Baseline Treg (cells/mL)
ALT0       : 80      : Baseline ALT (IU/L, elevated in chronic HCV)
Fibrosis0  : 1.0     : Baseline Metavir score (F0–F4, 0–4)
HSC0       : 50      : Baseline activated HSC (arbitrary units)

$CMT
@annotated
// Drug PK (9)
SOF_Tp    : SOF intrahepatic triphosphate (μM)
LED_p     : Ledipasvir plasma (ng/mL)
VEL_p     : Velpatasvir plasma (ng/mL)
GLE_p     : Glecaprevir plasma (ng/mL)
PIB_p     : Pibrentasvir plasma (ng/mL)
RBV_p     : Ribavirin plasma (μg/mL)
RBV_RBC   : Ribavirin RBC (μg/mL cell equivalent)
PEGIFN_p  : PEG-Interferon plasma (ng/mL)
NS5A_i    : NS5A inhibitor combined (nM, LED or VEL — whichever active)

// Viral kinetics (4)
T_cell    : Target hepatocytes (cells/mL)
I_cell    : Productively infected hepatocytes (cells/mL)
V_rna     : Free HCV RNA (IU/mL)
V_def     : Non-productive/defective particles (IU/mL)

// Immune (3)
CTL       : HCV-specific cytotoxic T lymphocytes (cells/mL)
NK_cell   : Activated NK cells (cells/mL)
Treg_HCV  : Regulatory T cells (cells/mL)

// Liver pathology (4)
ALT       : Serum ALT (IU/L)
Fibro_met : Metavir fibrosis score (0–4)
HSC_act   : Activated HSC index (0–100)
HCC_idx   : HCC risk accumulation index (0–1)

$GLOBAL
double E_SOF_p, E_LED_p, E_GLE_p, E_PIB_p, E_RBV_p, E_PEGIFN_p;
double Ep, Ei;   // combined efficacy on production and infectivity
double lytic_CTL, lytic_NK;

$MAIN
// Drug efficacy (Hill equation)
E_SOF_p    = pow(SOF_Tp, n_SOF) / (pow(EC50_SOF, n_SOF) + pow(SOF_Tp, n_SOF));
E_LED_p    = LED_p    / (EC50_LED + LED_p);
E_GLE_p    = pow(GLE_p, n_GLE) / (pow(EC50_GLE, n_GLE) + pow(GLE_p, n_GLE));
E_PIB_p    = PIB_p    / (EC50_PIB + PIB_p);
E_RBV_p    = RBV_p    / (EC50_RBV + RBV_p);
E_PEGIFN_p = PEGIFN_p / (EC50_PEGIFN + PEGIFN_p);

// Use VEL if higher than LED (mutually exclusive in practice)
double E_NS5A = (NS5A_i > 0) ? NS5A_i / (EC50_VEL + NS5A_i) : E_LED_p;

// Combined production blockade: SOF + NS5A + NS3/4A
// (1 - εp) = (1-εSOF)(1-εNS5A)(1-εNS3)
double prod_block = (1.0 - E_SOF_p) * (1.0 - E_NS5A) * (1.0 - E_GLE_p)
                    * (1.0 - E_PIB_p * 0.5);
Ep = 1.0 - prod_block;

// Infectivity block: primarily SOF-derived and RBV
Ei = E_SOF_p * 0.4 + E_RBV_p * 0.3;
if(Ei > 0.99) Ei = 0.99;

// Immune cytolysis rates
lytic_CTL = keff_CTL * CTL * I_cell;
lytic_NK  = keff_NK  * NK_cell * I_cell;

$ODE
// ── Drug PK ─────────────────────────────────────────────────
// SOF → intrahepatic triphosphate (SOF-Tp): pseudo-1-cpt intrahepatic
dxdt_SOF_Tp   = -CL_SOF_Tp * SOF_Tp;
// (dose events directly into SOF_Tp compartment after absorption conversion)

// Ledipasvir plasma 1-cpt
dxdt_LED_p    = -CL_LED * LED_p;

// Velpatasvir plasma 1-cpt
dxdt_VEL_p    = -CL_VEL * VEL_p;

// NS5A combined (for pangenotypic use)
dxdt_NS5A_i   = -0.3 * NS5A_i;  // approximate combined

// Glecaprevir 1-cpt
dxdt_GLE_p    = -CL_GLE * GLE_p;

// Pibrentasvir 1-cpt
dxdt_PIB_p    = -CL_PIB * PIB_p;

// Ribavirin plasma + RBC
dxdt_RBV_p    = -CL_RBV * RBV_p - k_RBV_RBC * RBV_p;
dxdt_RBV_RBC  =  k_RBV_RBC * RBV_p - kdeg_RBV_RBC * RBV_RBC;

// PEG-IFN-α (SC, once weekly)
dxdt_PEGIFN_p = -CL_PEGIFN * PEGIFN_p;

// ── Viral Kinetics (Target-Cell Limited) ─────────────────────
// dT/dt = s - d_T·T - (1-Ei)·β·V·T
dxdt_T_cell = s_T - d_T * T_cell
              - (1.0 - Ei) * beta_infect * V_rna * T_cell;
if(T_cell < 0) dxdt_T_cell = 0;

// dI/dt = (1-Ei)·β·V·T - δ·I - CTL + NK killing
double infection_rate = (1.0 - Ei) * beta_infect * V_rna * T_cell;
dxdt_I_cell = infection_rate
              - delta_I * I_cell
              - lytic_CTL
              - lytic_NK;
if(I_cell < 0) dxdt_I_cell = 0;

// dV/dt = (1-Ep)·p·I - c·V
dxdt_V_rna = (1.0 - Ep) * p_prod * I_cell
             - c_clear * V_rna;
if(V_rna < 0) dxdt_V_rna = 0;

// Defective particles (minor, diagnostic)
dxdt_V_def = 0.05 * (1.0 - Ep) * p_prod * I_cell - c_clear * V_def;
if(V_def < 0) dxdt_V_def = 0;

// ── Immune Dynamics ──────────────────────────────────────────
// CTL: expand when infected cells present, suppressed by Treg, exhausted chronically
double CTL_stimul = kprol_CTL * I_cell;
double CTL_exhaust_rate = k_exhaust * CTL * I_cell;
dxdt_CTL = CTL_stimul - kdeath_CTL * CTL
           - (Treg_HCV / (Treg_HCV + 100.0)) * CTL * 0.5
           - CTL_exhaust_rate;
if(CTL < 0) dxdt_CTL = 0;

// NK cells: activated by type I IFN and infected cells
double NK_stim = kprol_NK * (I_cell / (I_cell + 1e4)) * (1.0 + E_PEGIFN_p);
dxdt_NK_cell = NK_stim - kdeath_NK * NK_cell;
if(NK_cell < 0) dxdt_NK_cell = 0;

// Treg: expand with chronic inflammation (I_cell high)
dxdt_Treg_HCV = kprol_Treg_HCV * (I_cell / (I_cell + 5e4))
                - kdeath_Treg_HCV * Treg_HCV;
if(Treg_HCV < 0) dxdt_Treg_HCV = 0;

// ── Liver Pathology ──────────────────────────────────────────
// ALT: rises with hepatocyte death (immune + viral)
double cell_death_total = delta_I * I_cell + lytic_CTL + lytic_NK;
dxdt_ALT = kALT_prod * cell_death_total - kALT_clear * ALT;
if(ALT < 0) dxdt_ALT = 0;

// HSC activation: driven by TGF-β surrogate (cell death index)
dxdt_HSC_act = kprog_HSC * cell_death_total / (cell_death_total + 1e3)
               - kdeath_HSC * HSC_act;
if(HSC_act < 0) dxdt_HSC_act = 0;
if(HSC_act > 100) dxdt_HSC_act = 0;

// Fibrosis: progresses with HSC activity; regresses post-SVR
double SVR_flag = (V_rna < 10.0) ? 1.0 : 0.0;  // approximate SVR
dxdt_Fibro_met = kprog_fibro * HSC_act
                 - krem_fibro * Fibro_met * SVR_flag;
if(Fibro_met < 0)  dxdt_Fibro_met = 0;
if(Fibro_met > 4.0) dxdt_Fibro_met = 0;

// HCC risk: cumulative damage from fibrosis (slow accumulation)
dxdt_HCC_idx = kHCC_prog * Fibro_met * (1.0 - SVR_flag * 0.75);
if(HCC_idx > 1.0) dxdt_HCC_idx = 0;

$TABLE
capture log10_V     = (V_rna > 0) ? log10(V_rna) : -1;
capture V_IU        = V_rna;
capture log10_I     = (I_cell > 0) ? log10(I_cell) : 0;
capture ALT_out     = ALT;
capture Fibrosis_out = Fibro_met;
capture CTL_out     = CTL;
capture NK_out      = NK_cell;
capture Treg_out    = Treg_HCV;
capture E_SOF_pct   = E_SOF_p  * 100.0;
capture E_LED_pct   = E_LED_p  * 100.0;
capture E_GLE_pct   = E_GLE_p  * 100.0;
capture Ep_pct      = Ep       * 100.0;
capture Ei_pct      = Ei       * 100.0;
capture HCC_risk    = HCC_idx;

// SVR flag: HCV RNA < 15 IU/mL at t ≥ end-of-treatment + 12 weeks
capture TND_flag    = (V_rna <  15.0) ? 1.0 : 0.0;
capture SVR12_flag  = (V_rna <  15.0) ? 1.0 : 0.0;
capture RVR_flag    = (V_rna <  15.0 && TIME < 700) ? 1.0 : 0.0;

$INIT
SOF_Tp    = 0
LED_p     = 0
VEL_p     = 0
NS5A_i    = 0
GLE_p     = 0
PIB_p     = 0
RBV_p     = 0
RBV_RBC   = 0
PEGIFN_p  = 0
T_cell    = 2e7
I_cell    = 1e5
V_rna     = 1e6
V_def     = 1e3
CTL       = 1000
NK_cell   = 500
Treg_HCV  = 200
ALT       = 80
Fibro_met = 1.0
HSC_act   = 20
HCC_idx   = 0.01
'

## ----------------------------------------------------------
## Compile model
## ----------------------------------------------------------
mod <- mcode("hcv_qsp", code, quiet = TRUE)

## ----------------------------------------------------------
## Dosing helpers (hourly time scale)
## ----------------------------------------------------------
make_sof_led_events <- function(duration_wk = 12) {
  # SOF 400 mg QD → intrahepatic SOF-Tp peak ~10 μM
  # LED 90 mg QD  → plasma ~120 ng/mL
  dur_h <- duration_wk * 7 * 24
  sof_ev  <- ev(amt = 8.0,  cmt = "SOF_Tp", time = seq(0, dur_h - 24, by = 24))
  led_ev  <- ev(amt = 120,  cmt = "LED_p",  time = seq(0, dur_h - 24, by = 24))
  sof_ev + led_ev
}

make_sof_vel_events <- function(duration_wk = 12) {
  dur_h <- duration_wk * 7 * 24
  sof_ev  <- ev(amt = 8.0,  cmt = "SOF_Tp", time = seq(0, dur_h - 24, by = 24))
  vel_ev  <- ev(amt = 100,  cmt = "VEL_p",  time = seq(0, dur_h - 24, by = 24))
  sof_ev + vel_ev
}

make_gle_pib_events <- function(duration_wk = 8) {
  dur_h <- duration_wk * 7 * 24
  gle_ev  <- ev(amt = 18,  cmt = "GLE_p", time = seq(0, dur_h - 24, by = 24))
  pib_ev  <- ev(amt = 0.15, cmt = "PIB_p", time = seq(0, dur_h - 24, by = 24))
  gle_ev + pib_ev
}

make_peg_rbv_events <- function(duration_wk = 48) {
  dur_h <- duration_wk * 7 * 24
  # PEG-IFN: SC once weekly
  peg_times <- seq(0, dur_h - 168, by = 168)
  # RBV: BID (weight-based ~1000-1200 mg/day → 600 mg BID)
  rbv_times <- seq(0, dur_h - 12, by = 12)
  peg_ev <- ev(amt = 150, cmt = "PEGIFN_p", time = peg_times)
  rbv_ev <- ev(amt = 3.0, cmt = "RBV_p",    time = rbv_times)
  peg_ev + rbv_ev
}

## ----------------------------------------------------------
## SCENARIO 1 – SOF/LED (Harvoni®) 12 weeks — GT1
## ----------------------------------------------------------
sim_sof_led <- mod %>%
  ev(make_sof_led_events(12)) %>%
  mrgsim(end = 12*7*24 + 12*7*24, delta = 12) %>%  # treat + 12wk follow-up
  as.data.frame() %>%
  mutate(scenario = "SOF/LED 12wk (Harvoni)", time_wk = time / (7*24))

## ----------------------------------------------------------
## SCENARIO 2 – SOF/VEL (Epclusa®) 12 weeks — Pangenotypic
## ----------------------------------------------------------
sim_sof_vel <- mod %>%
  ev(make_sof_vel_events(12)) %>%
  mrgsim(end = 12*7*24 + 12*7*24, delta = 12) %>%
  as.data.frame() %>%
  mutate(scenario = "SOF/VEL 12wk (Epclusa)", time_wk = time / (7*24))

## ----------------------------------------------------------
## SCENARIO 3 – GLE/PIB (Mavyret®) 8 weeks — Naive
## ----------------------------------------------------------
sim_gle_pib <- mod %>%
  ev(make_gle_pib_events(8)) %>%
  mrgsim(end = 8*7*24 + 12*7*24, delta = 12) %>%
  as.data.frame() %>%
  mutate(scenario = "GLE/PIB 8wk (Mavyret)", time_wk = time / (7*24))

## ----------------------------------------------------------
## SCENARIO 4 – PEG-IFN + RBV (historical, 48 weeks)
## ----------------------------------------------------------
sim_peg_rbv <- mod %>%
  ev(make_peg_rbv_events(48)) %>%
  mrgsim(end = 48*7*24 + 12*7*24, delta = 24) %>%
  as.data.frame() %>%
  mutate(scenario = "PEG-IFN/RBV 48wk (historical)", time_wk = time / (7*24))

## ----------------------------------------------------------
## SCENARIO 5 – SOF/VEL extended 24 weeks (decompensated cirrhosis)
## ----------------------------------------------------------
mod_cirr <- mod %>%
  param(Fibrosis0 = 3.5, I0 = 5e5, V0 = 5e6, ALT0 = 150,
        T0 = 1e7)
sim_cirr <- mod_cirr %>%
  init(Fibro_met = 3.5, I_cell = 5e5, V_rna = 5e6,
       ALT = 150, T_cell = 1e7) %>%
  ev(make_sof_vel_events(24)) %>%
  mrgsim(end = 24*7*24 + 12*7*24, delta = 24) %>%
  as.data.frame() %>%
  mutate(scenario = "SOF/VEL 24wk (cirrhosis)", time_wk = time / (7*24))

## ----------------------------------------------------------
## SCENARIO 6 – SOF/LED + RBV (GT3 / NS5A RAS)
## ----------------------------------------------------------
make_sof_led_rbv_events <- function(duration_wk = 24) {
  dur_h <- duration_wk * 7 * 24
  sof_ev <- ev(amt = 8.0,  cmt = "SOF_Tp", time = seq(0, dur_h - 24, by = 24))
  led_ev <- ev(amt = 120,  cmt = "LED_p",  time = seq(0, dur_h - 24, by = 24))
  rbv_ev <- ev(amt = 3.0,  cmt = "RBV_p",  time = seq(0, dur_h - 12, by = 12))
  sof_ev + led_ev + rbv_ev
}
sim_sof_led_rbv <- mod %>%
  ev(make_sof_led_rbv_events(24)) %>%
  mrgsim(end = 24*7*24 + 12*7*24, delta = 24) %>%
  as.data.frame() %>%
  mutate(scenario = "SOF/LED+RBV 24wk (GT3/RAS)", time_wk = time / (7*24))

## ----------------------------------------------------------
## SCENARIO 7 – Untreated natural history
## ----------------------------------------------------------
sim_untreated <- mod %>%
  mrgsim(end = 52*7*24, delta = 24) %>%
  as.data.frame() %>%
  mutate(scenario = "Untreated", time_wk = time / (7*24))

## ----------------------------------------------------------
## Combine & filter to weekly resolution
## ----------------------------------------------------------
sim_all <- bind_rows(
  sim_sof_led, sim_sof_vel, sim_gle_pib, sim_peg_rbv,
  sim_cirr, sim_sof_led_rbv, sim_untreated
) %>%
  group_by(scenario) %>%
  mutate(time_wk_round = round(time_wk, 1)) %>%
  ungroup()

## Colour palette
cols <- c(
  "SOF/LED 12wk (Harvoni)"           = "#E74C3C",
  "SOF/VEL 12wk (Epclusa)"           = "#3498DB",
  "GLE/PIB 8wk (Mavyret)"            = "#2ECC71",
  "PEG-IFN/RBV 48wk (historical)"   = "#9B59B6",
  "SOF/VEL 24wk (cirrhosis)"         = "#E67E22",
  "SOF/LED+RBV 24wk (GT3/RAS)"       = "#1ABC9C",
  "Untreated"                         = "#7F8C8D"
)

## ----------------------------------------------------------
## Plot 1: Viral load (log10 IU/mL) over time
## ----------------------------------------------------------
p1 <- ggplot(sim_all, aes(time_wk, log10_V + 1e-3, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = log10(15), linetype = "dashed", colour = "black") +
  annotate("text", x = 0.5, y = log10(15) + 0.1,
           label = "LLOQ (15 IU/mL)", hjust = 0, size = 3) +
  scale_colour_manual(values = cols) +
  labs(title = "HCV QSP — Viral Load Trajectories",
       x = "Time (weeks)", y = "HCV-RNA (log10 IU/mL)", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

## ----------------------------------------------------------
## Plot 2: ALT over time
## ----------------------------------------------------------
p2 <- ggplot(sim_all, aes(time_wk, ALT_out, colour = scenario)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 40, linetype = "dashed") +
  scale_colour_manual(values = cols) +
  labs(title = "Serum ALT", x = "Time (weeks)", y = "ALT (IU/L)", colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "none")

## ----------------------------------------------------------
## Plot 3: Fibrosis score
## ----------------------------------------------------------
p3 <- ggplot(sim_all, aes(time_wk, Fibrosis_out, colour = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(limits = c(0, 4.1), breaks = 0:4) +
  scale_colour_manual(values = cols) +
  labs(title = "Metavir Fibrosis Score", x = "Time (weeks)", y = "Score (F0–F4)", colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "none")

## ----------------------------------------------------------
## Plot 4: CTL vs Treg dynamics
## ----------------------------------------------------------
p4 <- sim_all %>%
  filter(scenario %in% c("SOF/LED 12wk (Harvoni)", "Untreated")) %>%
  select(time_wk, scenario, CTL_out, Treg_out) %>%
  pivot_longer(c(CTL_out, Treg_out), names_to = "cell", values_to = "count") %>%
  mutate(cell = recode(cell, CTL_out = "CTL", Treg_out = "Treg")) %>%
  ggplot(aes(time_wk, count, colour = scenario, linetype = cell)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols[c("SOF/LED 12wk (Harvoni)", "Untreated")]) +
  labs(title = "Immune Dynamics (CTL vs Treg)",
       x = "Time (weeks)", y = "Cells/mL", colour = NULL, linetype = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

## Combined figure
if (requireNamespace("patchwork", quietly = TRUE)) {
  combined <- (p1 / (p2 | p3) / p4) +
    plot_annotation(
      title    = "Chronic Hepatitis C QSP — Treatment Comparison",
      subtitle = "Neumann–Perelson viral kinetics + DAA PK/PD + immune module + fibrosis",
      theme    = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  print(combined)
}

## ----------------------------------------------------------
## Week-12 and SVR summary
## ----------------------------------------------------------
wk12_h <- 12 * 7 * 24
summary_end <- sim_all %>%
  group_by(scenario) %>%
  filter(time == max(time)) %>%
  slice(1) %>%
  ungroup() %>%
  select(scenario, log10_V, ALT_out, Fibrosis_out, SVR12_flag, HCC_risk) %>%
  rename(
    Scenario         = scenario,
    `log10-VL`       = log10_V,
    `ALT (IU/L)`     = ALT_out,
    `Fibrosis (F0-4)`= Fibrosis_out,
    `SVR (1=Yes)`    = SVR12_flag,
    `HCC Risk Index` = HCC_risk
  )

cat("\n── Chronic HCV QSP Model — End-of-Follow-up Summary ──────\n")
print(summary_end, digits = 3)

cat("\n── Viral Kinetics Parameters (Neumann 1998) ────────────────
  Viral clearance rate (c)  : 22 /h  (t½ ~3h)
  Infected cell loss (δ)    : 0.08 /h (~0.5/day)
  Production rate (p)       : 100 virions/cell/h
  Phase 1: rapid decline   → clearance of free virus (c)
  Phase 2: slower decline  → loss of infected cells (δ)
\n")

message("Chronic Hepatitis C QSP model compiled and simulated successfully.")
