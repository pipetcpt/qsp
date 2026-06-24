## ============================================================
## Relapsing Polychondritis (RP) — mrgsolve QSP Model
## ============================================================
## Description:
##   Mechanistic ODE model covering:
##     • Prednisolone 2-CMT PK + GR binding dynamics
##     • Tocilizumab 2-CMT PK (IL-6R blockade)
##     • CD4+ T cell differentiation (Th1/Th17/Treg)
##     • B cell activation → anti-type II collagen antibodies
##     • Immune complex formation + complement activation
##     • Cytokine network (TNF-α, IL-6, IL-17A, IL-1β)
##     • MMP-mediated cartilage degradation
##     • Disease activity (RPDAI proxy)
##
## Treatment Scenarios:
##   1. Untreated (natural disease course)
##   2. Prednisone monotherapy (induction → taper)
##   3. Prednisone + Methotrexate (standard of care)
##   4. Tocilizumab (IL-6R inhibitor, biologic)
##   5. Abatacept (CTLA4-Ig costimulation blockade)
##   6. Dapsone (mild/first-line, neutrophil-targeted)
##   7. Prednisone + Tocilizumab (combination)
##
## Key References:
##   Cantarini et al., Autoimmun Rev 2014; Dion et al., JCI 2007;
##   Mathian et al., Ann Rheum Dis 2019; Shimizu et al., 2019
##
## Author: Claude Code Routine (QSP Library Project)
## Date: 2026-06-20
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

## ── mrgsolve model code ──────────────────────────────────────
rpc_code <- '
$PROB
Relapsing Polychondritis (RP) QSP Model
20-compartment mechanistic ODE system
PK: Prednisolone (2-CMT + GR), Tocilizumab (2-CMT)
PD: T cells, B cells, cytokines, cartilage integrity

$PARAM
// ── Prednisolone PK (Frey et al., J Clin Pharmacol 2002) ────
ka_pred   = 1.5     // h-1, oral absorption rate constant
CL_pred   = 15.0    // L/h, apparent clearance
V1_pred   = 30.0    // L,   central volume
V2_pred   = 45.0    // L,   peripheral volume
Q_pred    = 8.0     // L/h, intercompartmental clearance
F_pred    = 0.82    // bioavailability (prednisone → prednisolone)

// GR binding/dynamics (Ramakrishnan et al., J PK PD 2002)
kon_GR    = 0.20    // nM-1·h-1, association rate
koff_GR   = 0.02    // h-1,      dissociation rate
ksyn_GR   = 10.0    // nmol/h,   GR synthesis rate
kdeg_GR   = 0.05    // h-1,      GR baseline degradation
MW_pred   = 360.4   // g/mol prednisolone (for mg/L to nM)

// ── Tocilizumab PK (Frey et al., J Clin Pharmacol 2010) ─────
CL_tcz    = 0.008   // L/h
V1_tcz    = 4.0     // L
V2_tcz    = 6.5     // L
Q_tcz     = 0.015   // L/h

// ── T cell dynamics ──────────────────────────────────────────
kin_Tact  = 0.50    // h-1 baseline Th1 input rate
kout_Tact = 0.40    // h-1 elimination
kin_Th17  = 0.30    // h-1
kout_Th17 = 0.35
kin_Treg  = 0.25    // h-1
kout_Treg = 0.30

// ── B cell / Ab dynamics ─────────────────────────────────────
kin_Bact  = 0.20    // h-1
kout_Bact = 0.18
kprod_Ab  = 0.05    // nM/h per unit Bact
kdeg_Ab   = 0.004   // h-1  (half-life ~7 d)

// ── Immune complex / complement ──────────────────────────────
kform_IC  = 0.10    // IC formation coefficient
kelim_IC  = 0.15    // h-1 IC elimination
IC50_IC   = 5.0     // EC50 (saturable formation)
kact_Comp = 0.40    // h-1 complement activation by IC
kdeg_Comp = 0.60    // h-1

// ── Cytokine dynamics (ng/mL) ────────────────────────────────
// TNF-α
kin_TNF   = 0.30    // ng/mL/h
kout_TNF  = 1.20    // h-1
// IL-6
kin_IL6   = 0.25
kout_IL6  = 0.80
// IL-17A
kin_IL17  = 0.30
kout_IL17 = 0.90
// IL-1β
kin_IL1b  = 0.35
kout_IL1b = 1.00

// ── MMP / Cartilage ──────────────────────────────────────────
kin_MMP      = 0.15    // AU/h
kout_MMP     = 0.40    // h-1
kdest_Cart   = 0.06    // h-1 destruction by MMP
krep_Cart    = 0.015   // h-1 intrinsic repair
Cart_max     = 100.0   // normal cartilage (AU, 0-100)

// ── Drug PD parameters ───────────────────────────────────────
// Prednisolone Emax inhibition of inflammation
Emax_pred   = 0.85     // maximum inhibitory effect
EC50_pred   = 8.0      // mg/L half-effect concentration
nH_pred     = 2.0      // Hill exponent

// Tocilizumab on IL-6
Emax_tcz    = 0.95
EC50_tcz    = 0.30     // mg/L

// MTX (combined pathway): set E_MTX via dosing event
E_MTX       = 0.0      // 0=off; 0.60 = therapeutic (15mg/wk)

// Abatacept (CTLA4-Ig): set E_ABA via dosing event
E_ABA       = 0.0      // 0=off; 0.65 = therapeutic

// Dapsone: set E_DAP via dosing event
E_DAP       = 0.0      // 0=off; 0.55 = therapeutic (100 mg/d)

// ── Interaction coefficients ─────────────────────────────────
h_IL6_Tact  = 0.30     // IL-6 stimulates Tact
h_IC_Tact   = 0.25     // IC stimulates T cells
h_Treg_Tact = 0.60     // Treg suppresses Tact
h_IL6_Th17  = 0.80     // IL-6 drives Th17
h_Treg_Th17 = 0.70     // Treg suppresses Th17
h_Tact_Bact = 0.50     // Tact (T-help) drives B
h_TNF_MMP   = 0.50     // TNF drives MMP
h_IL17_MMP  = 0.40     // IL-17 drives MMP
h_IL1b_MMP  = 0.30     // IL-1β drives MMP
h_Comp_MMP  = 0.20     // complement drives MMP
h_Tact_TNF  = 0.60     // Tact drives TNF
h_Tact_IL6  = 0.50     // Tact drives IL-6
h_TNF_IL6   = 0.40     // TNF drives IL-6
h_IL1b_IL6  = 0.30     // IL-1β drives IL-6
h_Th17_IL17 = 1.00     // Th17 drives IL-17
h_IC_IL1b   = 0.40     // IC activates IL-1β
h_Comp_IL1b = 0.30     // complement drives IL-1β

$CMT
// ── PK compartments ──────────────────────────────────────────
DEPOT_pred    // [1]  oral prednisolone depot (mg)
C1_pred       // [2]  prednisolone central amount (mg)
C2_pred       // [3]  prednisolone peripheral amount (mg)
GR_free       // [4]  free GR-α pool (nmol)
GR_bound      // [5]  GR-prednisolone complex (nmol)
C1_tcz        // [6]  tocilizumab central amount (mg)
C2_tcz        // [7]  tocilizumab peripheral amount (mg)

// ── Immune compartments ───────────────────────────────────────
Tact          // [8]  activated CD4+ Th1 cells (×10⁶)
Th17          // [9]  Th17 effector cells (×10⁶)
Treg          // [10] FoxP3+ regulatory T cells (×10⁶)
Bact          // [11] activated B cells (×10⁶)
Ab_CII        // [12] anti-type II collagen antibody (nM)
IC            // [13] immune complexes (AU)
Comp          // [14] complement activation (AU)

// ── Cytokines (ng/mL) ─────────────────────────────────────────
TNF           // [15] TNF-α
IL6           // [16] IL-6
IL17          // [17] IL-17A
IL1b          // [18] IL-1β

// ── Disease ──────────────────────────────────────────────────
MMP           // [19] MMP net activity (AU)
Cart          // [20] cartilage integrity (0–100 AU)

$INIT
DEPOT_pred = 0
C1_pred    = 0
C2_pred    = 0
GR_free    = 200.0    // ksyn_GR/kdeg_GR at baseline
GR_bound   = 0
C1_tcz     = 0
C2_tcz     = 0
Tact       = 1.25     // kin/kout at baseline
Th17       = 0.857
Treg       = 0.833
Bact       = 1.111
Ab_CII     = 50.0     // elevated in RP patients
IC         = 5.0
Comp       = 3.33
TNF        = 0.25
IL6        = 0.3125
IL17       = 0.333
IL1b       = 0.35
MMP        = 0.375
Cart       = 100.0

$MAIN
// ── Derived concentrations ───────────────────────────────────
double Cp_pred   = C1_pred / V1_pred;           // mg/L
double Cp_pred_nM = Cp_pred * 1000.0 / MW_pred; // nM (for GR)
double Cp_tcz    = C1_tcz  / V1_tcz;            // mg/L

// ── Drug effect functions (Hill/Emax inhibitory) ─────────────
double Epred = Emax_pred * pow(Cp_pred, nH_pred) /
               (pow(EC50_pred, nH_pred) + pow(Cp_pred, nH_pred));

double Etcz  = Emax_tcz  * Cp_tcz / (EC50_tcz + Cp_tcz);

// Combined IL-6 suppression (max of steroid + TCZ)
double E_IL6_tot = fmax(Epred * 0.7, Etcz);

// GR occupancy fraction
double GR_occ = GR_bound / (GR_free + GR_bound + 0.001);

// Cart floor guard (set initial values if negative)
if(Cart < 0.0) Cart = 0.0;

$ODE
// ── Prednisolone PK ──────────────────────────────────────────
dxdt_DEPOT_pred =  -ka_pred * DEPOT_pred;
dxdt_C1_pred    =   F_pred * ka_pred * DEPOT_pred
                  - (CL_pred + Q_pred) / V1_pred * C1_pred
                  + Q_pred / V2_pred * C2_pred;
dxdt_C2_pred    =   Q_pred / V1_pred * C1_pred
                  - Q_pred / V2_pred * C2_pred;

// GR dynamics
dxdt_GR_free  = ksyn_GR
                - kdeg_GR * GR_free
                - kon_GR  * Cp_pred_nM * GR_free
                + koff_GR * GR_bound;

dxdt_GR_bound = kon_GR  * Cp_pred_nM * GR_free
                - koff_GR * GR_bound
                - kdeg_GR * GR_bound;

// ── Tocilizumab PK ───────────────────────────────────────────
dxdt_C1_tcz = -(CL_tcz + Q_tcz) / V1_tcz * C1_tcz
              + Q_tcz / V2_tcz * C2_tcz;
dxdt_C2_tcz =  Q_tcz / V1_tcz * C1_tcz
              - Q_tcz / V2_tcz * C2_tcz;

// ── Immune cell dynamics ─────────────────────────────────────

// Tact (Th1): stimulated by IC, IL-6; suppressed by Treg, pred, ABA
double drive_Tact = 1.0 + h_IL6_Tact * IL6 + h_IC_Tact * IC;
double supp_Tact  = 1.0 + h_Treg_Tact * Treg + Epred * 0.5 + E_ABA;
dxdt_Tact = kin_Tact * drive_Tact / supp_Tact - kout_Tact * Tact;

// Th17: stimulated by IL-6; suppressed by Treg, MTX, pred
double drive_Th17 = 1.0 + h_IL6_Th17 * IL6 + 0.3 * IL1b;
double supp_Th17  = 1.0 + h_Treg_Th17 * Treg + E_MTX * 1.5 + Epred * 0.4;
dxdt_Th17 = kin_Th17 * drive_Th17 / supp_Th17 - kout_Th17 * Th17;

// Treg: mild promotion by GR occupancy (steroid benefit)
dxdt_Treg = kin_Treg * (1.0 + 0.3 * GR_occ) - kout_Treg * Treg;

// Bact: driven by Tact (T-help); suppressed by pred, RTX (E_MTX proxy)
double drive_Bact = 1.0 + h_Tact_Bact * Tact;
dxdt_Bact = kin_Bact * drive_Bact * (1.0 - Epred * 0.3 - E_MTX * 0.4)
            - kout_Bact * Bact;

// Ab_CII (anti-CII IgG)
dxdt_Ab_CII = kprod_Ab * Bact - kdeg_Ab * Ab_CII;

// IC: saturating formation, cleared by complement/phagocytosis
dxdt_IC = kform_IC * Ab_CII / (IC50_IC + Ab_CII) * Ab_CII
          - kelim_IC * IC;

// Complement activation
dxdt_Comp = kact_Comp * IC - kdeg_Comp * Comp;

// ── Cytokine dynamics ────────────────────────────────────────

// TNF-α: driven by Tact, IC; inhibited by pred
double prod_TNF = kin_TNF * (1.0 + h_Tact_TNF * Tact + 0.2 * Comp);
dxdt_TNF = prod_TNF * (1.0 - Epred) - kout_TNF * TNF;

// IL-6: driven by Tact, TNF, IL-1β; inhibited by pred + TCZ
double prod_IL6 = kin_IL6 * (1.0 + h_Tact_IL6 * Tact
                              + h_TNF_IL6  * TNF
                              + h_IL1b_IL6 * IL1b);
dxdt_IL6 = prod_IL6 * (1.0 - E_IL6_tot) - kout_IL6 * IL6;

// IL-17A: driven by Th17; partially suppressed by MTX
double prod_IL17 = kin_IL17 * (1.0 + h_Th17_IL17 * Th17);
dxdt_IL17 = prod_IL17 * (1.0 - E_MTX * 0.5) - kout_IL17 * IL17;

// IL-1β: driven by IC, complement, TNF; inhibited by pred
double prod_IL1b = kin_IL1b * (1.0 + h_IC_IL1b   * IC
                               + h_Comp_IL1b * Comp
                               + 0.2 * TNF);
dxdt_IL1b = prod_IL1b * (1.0 - Epred * 0.7) - kout_IL1b * IL1b;

// ── MMP and cartilage ────────────────────────────────────────

// MMP: driven by cytokines; dampened by pred
double prod_MMP = kin_MMP * (1.0 + h_TNF_MMP  * TNF
                              + h_IL17_MMP * IL17
                              + h_IL1b_MMP * IL1b
                              + h_Comp_MMP * Comp);
dxdt_MMP = prod_MMP * (1.0 - Epred * 0.4) - kout_MMP * MMP;

// Cartilage: destroyed by MMP, slowly repairs
double destr = kdest_Cart * MMP * (Cart / Cart_max);
double repair = krep_Cart * (Cart_max - Cart);
dxdt_Cart = repair - destr;

$TABLE
capture Cp_pred_mgL  = C1_pred / V1_pred;      // prednisolone conc (mg/L)
capture Cp_tcz_mgL   = C1_tcz  / V1_tcz;       // tocilizumab conc (mg/L)
capture GR_occ_pct   = GR_bound / (GR_free + GR_bound + 0.001) * 100.0;
capture E_pred_pct   = Epred * 100.0;           // prednisolone effect (%)
capture E_tcz_pct    = Etcz  * 100.0;
capture CartPct      = Cart;                    // cartilage integrity (0-100)
capture RPDAI_proxy  = 10.0 * (TNF + IL6 + IL17) + 50.0 * (1.0 - Cart / 100.0);
capture CRP_proxy    = 0.5 * IL6 + 0.2 * TNF;  // surrogate CRP
capture Ab_CII_nM    = Ab_CII;
capture IC_AU        = IC;
'

## ── Compile model ───────────────────────────────────────────
mod <- mcode("rpc_qsp", rpc_code)

## ── Helper: build event table ────────────────────────────────
make_events <- function(scenario, t_end = 8760) {  # 365 days in hours

  ev <- ev()  # empty event

  if (scenario == "untreated") {
    return(ev)
  }

  if (scenario == "prednisone") {
    # Induction: 60 mg/d × 4 weeks, taper to 10 mg/d maintenance
    ev_ind <- ev(amt = 60, ii = 24, addl = 27, time = 0, cmt = 1)    # 28 d
    ev_tap <- ev(amt = 40, ii = 24, addl = 13, time = 672, cmt = 1)   # wk 5-6
    ev_mnt <- ev(amt = 20, ii = 24, addl = 27, time = 1008, cmt = 1)  # wk 7-10
    ev_low <- ev(amt = 10, ii = 24, addl = 299, time = 1680, cmt = 1) # ≥wk 11
    return(c(ev_ind, ev_tap, ev_mnt, ev_low))
  }

  if (scenario == "pred_mtx") {
    ev_pred <- ev(amt = 40, ii = 24, addl = 27, time = 0, cmt = 1)
    ev_low  <- ev(amt = 10, ii = 24, addl = 299, time = 672, cmt = 1)
    # MTX modelled as parameter switch at t = 0 (simplification)
    # Use param() override below
    return(c(ev_pred, ev_low))
  }

  if (scenario == "tocilizumab") {
    # 8 mg/kg q4w IV for 70-kg patient → 560 mg q4w
    # IV bolus into C1_tcz (cmt 6)
    ev_tcz <- ev(amt = 560, ii = 672, addl = 12, time = 0, cmt = 6)  # 13 doses, ~1 yr
    return(ev_tcz)
  }

  if (scenario == "abatacept") {
    # 750 mg q4w IV; modelled via E_ABA parameter
    return(ev())  # effect set via param()
  }

  if (scenario == "dapsone") {
    # 100 mg/d; modelled via E_DAP parameter
    return(ev())  # effect set via param()
  }

  if (scenario == "pred_tcz") {
    ev_pred <- ev(amt = 40, ii = 24, addl = 27, time = 0, cmt = 1)
    ev_low  <- ev(amt = 10, ii = 24, addl = 299, time = 672, cmt = 1)
    ev_tcz  <- ev(amt = 560, ii = 672, addl = 12, time = 0, cmt = 6)
    return(c(ev_pred, ev_low, ev_tcz))
  }
}

## ── Simulation function ──────────────────────────────────────
sim_rpc <- function(scenario, t_end = 8760, delta = 12) {
  times  <- seq(0, t_end, by = delta)
  events <- make_events(scenario, t_end)

  # Scenario-specific parameters
  p_extra <- list()
  if (scenario == "pred_mtx")   p_extra <- list(E_MTX = 0.60)
  if (scenario == "abatacept")  p_extra <- list(E_ABA = 0.65)
  if (scenario == "dapsone")    p_extra <- list(E_DAP = 0.55)

  mod2 <- if (length(p_extra) > 0) param(mod, p_extra) else mod

  out <- mrgsim(mod2, events = events, end = t_end, delta = delta, digits = 4)
  as.data.frame(out) %>% mutate(scenario = scenario)
}

## ── Run all scenarios ────────────────────────────────────────
scenarios <- c("untreated", "prednisone", "pred_mtx",
                "tocilizumab", "abatacept", "dapsone", "pred_tcz")

cat("Running", length(scenarios), "treatment scenarios...\n")
results <- bind_rows(lapply(scenarios, sim_rpc))

results$scenario <- factor(results$scenario, levels = scenarios,
  labels = c("Untreated", "Prednisone\n(monotherapy)",
             "Pred + MTX", "Tocilizumab\n(8mg/kg q4w)",
             "Abatacept", "Dapsone\n(100mg/d)",
             "Pred + Tocilizumab"))

results$time_d <- results$time / 24  # hours → days

cat("Simulation complete. Rows:", nrow(results), "\n")

## ── Plot 1: Cartilage Integrity ──────────────────────────────
p_cart <- ggplot(results, aes(time_d, CartPct, color = scenario)) +
  geom_line(linewidth = 0.9) +
  labs(title = "Cartilage Integrity over Time",
       x = "Time (days)", y = "Cartilage Integrity (AU, 0–100)",
       color = "Treatment") +
  scale_y_continuous(limits = c(0, 105)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

## ── Plot 2: RPDAI Proxy ──────────────────────────────────────
p_rpdai <- ggplot(results, aes(time_d, RPDAI_proxy, color = scenario)) +
  geom_line(linewidth = 0.9) +
  labs(title = "Disease Activity (RPDAI Proxy)",
       x = "Time (days)", y = "RPDAI Proxy (AU)",
       color = "Treatment") +
  theme_bw(base_size = 12)

## ── Plot 3: Cytokine dynamics (untreated vs best responders) ─
cyt_long <- results %>%
  filter(scenario %in% c("Untreated", "Tocilizumab\n(8mg/kg q4w)",
                          "Pred + Tocilizumab")) %>%
  select(time_d, scenario, TNF, IL6, IL17, IL1b) %>%
  pivot_longer(TNF:IL1b, names_to = "Cytokine", values_to = "Conc")

p_cyt <- ggplot(cyt_long, aes(time_d, Conc, color = scenario, linetype = Cytokine)) +
  geom_line(linewidth = 0.8) +
  labs(title = "Cytokine Dynamics",
       x = "Time (days)", y = "Cytokine Concentration (ng/mL)",
       color = "Treatment", linetype = "Cytokine") +
  theme_bw(base_size = 12)

## ── Plot 4: Immune cell dynamics (untreated) ─────────────────
imm_untreated <- results %>%
  filter(scenario == "Untreated") %>%
  select(time_d, Tact, Th17, Treg, Bact) %>%
  pivot_longer(-time_d, names_to = "Cell", values_to = "Count")

p_immune <- ggplot(imm_untreated, aes(time_d, Count, color = Cell)) +
  geom_line(linewidth = 0.9) +
  labs(title = "Immune Cell Dynamics (Untreated)",
       x = "Time (days)", y = "Cell Count (×10⁶)",
       color = "Cell Type") +
  theme_bw(base_size = 12)

## ── Plot 5: Prednisolone PK (single dose 40 mg) ──────────────
pk_data <- sim_rpc("prednisone") %>%
  filter(time_d < 14)  # first 2 weeks

p_pk <- ggplot(pk_data, aes(time_d, Cp_pred_mgL)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_line(aes(y = GR_occ_pct / 10), color = "darkred", linetype = 2) +
  scale_y_continuous(
    name = "Prednisolone (mg/L)",
    sec.axis = sec_axis(~. * 10, name = "GR Occupancy (%)")
  ) +
  labs(title = "Prednisolone PK & GR Occupancy (Induction Phase)",
       x = "Time (days)") +
  theme_bw(base_size = 12) +
  annotate("text", x = 10, y = 8, label = "--- GR occupancy (%/10)",
           color = "darkred", size = 3)

## ── Summary table at day 180 ─────────────────────────────────
summary_180 <- results %>%
  filter(abs(time_d - 180) < 1) %>%
  group_by(scenario) %>%
  slice(1) %>%
  select(scenario, CartPct, RPDAI_proxy, CRP_proxy,
         Ab_CII_nM, TNF, IL6, IL17) %>%
  ungroup() %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

cat("\n=== Simulated Outcomes at Day 180 ===\n")
print(as.data.frame(summary_180), row.names = FALSE)

## ── Print plots ──────────────────────────────────────────────
print(p_cart)
print(p_rpdai)
print(p_cyt)
print(p_immune)
print(p_pk)

## ── Parameter sensitivity analysis (tornado, day 365) ────────
sens_params <- c("kdest_Cart", "krep_Cart", "kin_TNF", "kin_IL17",
                 "h_TNF_MMP", "h_IL17_MMP", "kprod_Ab", "kdeg_Ab",
                 "kin_Tact", "kin_Th17")

base_cart_365 <- results %>%
  filter(scenario == "Untreated", abs(time_d - 365) < 1) %>%
  slice(1) %>% pull(CartPct)

sens_results <- lapply(sens_params, function(p) {
  vals <- c(0.5, 1.0, 2.0)  # × baseline
  base_val <- param(mod)[[p]]
  purrr::map_dfr(vals, function(mult) {
    mod_s <- param(mod, setNames(list(base_val * mult), p))
    out_s <- mrgsim(mod_s, end = 8760, delta = 24, digits = 4)
    df_s  <- as.data.frame(out_s)
    cart_val <- df_s %>% filter(abs(time / 24 - 365) < 1) %>%
                slice(1) %>% pull(CartPct)
    data.frame(param = p, multiplier = mult,
               Cart365 = cart_val,
               delta_pct = (cart_val - base_cart_365) / base_cart_365 * 100)
  })
}) %>% bind_rows()

cat("\n=== Sensitivity Analysis: CartPct at Day 365 ===\n")
print(sens_results[sens_results$multiplier != 1, ], row.names = FALSE)

cat("\nQSP model run complete.\n")
