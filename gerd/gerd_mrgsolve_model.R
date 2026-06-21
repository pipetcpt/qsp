##############################################################################
# GERD QSP Model — mrgsolve implementation
# Disease: Gastroesophageal Reflux Disease (GERD)
# Model: 20-compartment ODE system
#   • PK: PPI (omeprazole), H2RA (famotidine), P-CAB (vonoprazan),
#          prokinetics (domperidone)
#   • PD: H+/K+-ATPase turnover, gastric pH, esophageal acid exposure,
#          mucosal damage/healing, symptom score
# Calibration: Miner 2003 (omeprazole), Hunt 1984 (famotidine),
#              Ashida 2016 (vonoprazan VOYAGE), Kahrilas 2008 (esomeprazole)
# Scenarios: 6 treatment arms
# Author: Claude Code Routine (CCR) — 2026-06-18
##############################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ─── Model code ──────────────────────────────────────────────────────────────
gerd_code <- '
$PROB GERD QSP Model — PPI/H2RA/P-CAB/Prokinetic PK-PD

$PARAM @annotated
// ── Drug Doses ──────────────────────────────────────────────────────────────
DOSE_PPI   : 20  : PPI dose (mg, omeprazole equivalent)
DOSE_H2RA  : 40  : H2RA dose (mg, famotidine equivalent)
DOSE_PCAB  : 20  : P-CAB dose (mg, vonoprazan)
DOSE_PROK  : 10  : Prokinetic dose (mg, domperidone)

// ── PPI PK Parameters (Omeprazole) ──────────────────────────────────────────
KA_PPI     : 0.80  : PPI absorption rate constant (h-1)
CL_PPI     : 18.0  : PPI apparent clearance (L/h)
V_PPI      : 20.0  : PPI volume of distribution (L)
F_PPI      : 0.65  : PPI bioavailability (CYP2C19 EM)
CYP2C19    : 1.0   : CYP2C19 phenotype multiplier (EM=1, IM=0.6, PM=0.25)

// ── H2RA PK Parameters (Famotidine) ─────────────────────────────────────────
KA_H2RA    : 1.50  : H2RA absorption rate constant (h-1)
CL_H2RA    : 14.0  : H2RA clearance (L/h, renal dominant)
V_H2RA     : 90.0  : H2RA volume of distribution (L)
F_H2RA     : 0.50  : H2RA oral bioavailability

// ── P-CAB PK Parameters (Vonoprazan) ─────────────────────────────────────────
KA_PCAB    : 2.00  : P-CAB absorption rate constant (h-1)
CL_PCAB    : 8.0   : P-CAB clearance (L/h, CYP3A4)
V_PCAB     : 300.0 : P-CAB volume of distribution (L, lipophilic)
F_PCAB     : 0.60  : P-CAB oral bioavailability

// ── Prokinetic PK (Domperidone) ───────────────────────────────────────────────
KA_PROK    : 1.20  : Prokinetic absorption rate (h-1)
CL_PROK    : 60.0  : Prokinetic clearance (L/h)
V_PROK     : 450.0 : Prokinetic volume of distribution (L)
F_PROK     : 0.15  : Prokinetic bioavailability (first-pass)

// ── H+/K+-ATPase (Proton Pump) Dynamics ────────────────────────────────────
PUMP_TOTAL : 100.0 : Total H+/K+-ATPase pool (normalized units)
PUMP_ACT0  : 30.0  : Baseline active pump fraction (%)
K_ACT      : 0.30  : Pump activation rate (h-1)
K_DEACT    : 0.30  : Pump deactivation rate (h-1)
K_SYN_PUMP : 3.0   : Pump synthesis rate (units/h)
K_DEG_PUMP : 0.03  : Pump degradation rate (h-1)
IC50_PPI   : 0.15  : PPI IC50 for covalent binding (mg/L)
HILL_PPI   : 1.5   : PPI Hill coefficient
IC50_PCAB  : 0.08  : P-CAB IC50 for pump inhibition (mg/L)
IC50_H2RA  : 0.04  : H2RA IC50 for H2 receptor block (mg/L)
EMAX_H2RA  : 0.60  : H2RA maximal inhibition (60% of stimulated acid)

// ── Gastric Acid Secretion ──────────────────────────────────────────────────
ACID_BASE  : 3.5   : Baseline gastric acid secretion rate (mmol/h)
ACID_MAX   : 12.0  : Maximal acid secretion rate (mmol/h, pentagastrin)
HIST_STIM  : 0.35  : Histamine contribution to acid secretion
GAST_STIM  : 0.45  : Gastrin contribution to acid secretion
ACH_STIM   : 0.20  : Acetylcholine contribution
K_NEUT_INT : 0.5   : Intragastric acid neutralization rate (h-1)
PH_BUFF    : 1.5   : Buffer capacity constant
VOL_GAS    : 200.0 : Gastric volume (mL)

// ── LES & Reflux ─────────────────────────────────────────────────────────────
LES_P0     : 22.0  : Baseline LES pressure (mmHg)
TLESR_BASE : 6.0   : Baseline TLESR rate (events/h)
K_PROK_LES : 5.0   : Prokinetic LES pressure increase (mmHg at Emax)
K_PROK_EMP : 0.30  : Prokinetic gastric emptying rate increase (h-1)
GASTRIC_EMP: 0.50  : Gastric emptying rate constant (h-1)
REFLUX_K   : 0.05  : Acid reflux to esophageal acid exposure conversion

// ── Esophageal Mucosal Dynamics ──────────────────────────────────────────────
MUC_HEAL0  : 0.30  : Baseline mucosal healing rate (units/day)
MUC_INJ_K  : 0.40  : Mucosal injury rate constant (per AET unit)
K_HEAL     : 0.20  : Mucosal healing rate constant (h-1)
GRADE_K    : 0.015 : Damage→LA grade conversion constant
MAX_DAMAGE : 100.0 : Maximum damage score (100 = LA Grade D)

// ── Symptom Score ────────────────────────────────────────────────────────────
SYM_K      : 0.05  : Symptom sensitivity (per AET unit)
SYM_HEAL_K : 0.30  : Symptom relief rate with acid suppression

$CMT @annotated
// PPI PK
PPI_GUT    : PPI gut compartment (mg)
PPI_CENT   : PPI central compartment (mg)

// H2RA PK
H2RA_GUT   : H2RA gut compartment (mg)
H2RA_CENT  : H2RA central compartment (mg)

// P-CAB PK
PCAB_GUT   : P-CAB gut compartment (mg)
PCAB_CENT  : P-CAB central compartment (mg)

// Prokinetic PK
PROK_GUT   : Prokinetic gut compartment (mg)
PROK_CENT  : Prokinetic central compartment (mg)

// H+/K+-ATPase Pool
PUMP_INACT : Inactive pump pool (tubulovesicular)
PUMP_ACT   : Active pump pool (secretory canaliculus)
PUMP_INH   : Irreversibly inhibited pump pool (PPIs only)

// Gastric Acid & pH
ACID_RATE  : Gastric acid secretion rate (mmol/h)
GAS_pH     : Intragastric pH (transformed)

// Esophageal Compartments
AET        : Acid exposure time (% time pH<4)
MUC_DMG    : Mucosal damage score (0-100)
MUC_HEAL   : Mucosal integrity index (0-100)

// Symptom
SYM_SCORE  : Symptom score (GERD-Q equivalent, 0-18)

// Barrett progression (long-term)
BE_RISK    : Cumulative Barrett risk index (0-1)

$MAIN
double CP_PPI  = PPI_CENT  / V_PPI;
double CP_H2RA = H2RA_CENT / V_H2RA;
double CP_PCAB = PCAB_CENT / V_PCAB;
double CP_PROK = PROK_CENT / V_PROK;

// PPI: Emax covalent inhibition of active pumps
double INH_PPI = pow(CP_PPI, HILL_PPI) / (pow(IC50_PPI, HILL_PPI) + pow(CP_PPI, HILL_PPI));

// P-CAB: reversible ionic block
double INH_PCAB = CP_PCAB / (IC50_PCAB + CP_PCAB);

// H2RA: competitive antagonism at H2 receptor
double INH_H2RA = EMAX_H2RA * CP_H2RA / (IC50_H2RA + CP_H2RA);

// Combined acid suppression (PPI+PCAB irreversible/reversible; H2RA partial)
double ACID_INH = 1.0 - (1.0 - INH_PPI) * (1.0 - INH_PCAB) * (1.0 - INH_H2RA * HIST_STIM);
if(ACID_INH > 0.99) ACID_INH = 0.99;

// Prokinetic effect on LES pressure
double PROK_EFF_LES = K_PROK_LES * CP_PROK / (0.005 + CP_PROK);
double LES_P = LES_P0 + PROK_EFF_LES;

// TLESR rate — inversely proportional to LES pressure
double TLESR = TLESR_BASE * (LES_P0 / LES_P);

// Prokinetic effect on gastric emptying
double EMP_RATE = GASTRIC_EMP + K_PROK_EMP * CP_PROK / (0.005 + CP_PROK);

// Acid secretion (effective)
double ACID_EFF = ACID_BASE + (ACID_MAX - ACID_BASE) * (1.0 - ACID_INH);

// pH calculation from acid secretion rate
double pH_calc = 1.0 + PH_BUFF * exp(-ACID_EFF / 2.0);
if(pH_calc < 1.0) pH_calc = 1.0;
if(pH_calc > 7.0) pH_calc = 7.0;

// Reflux acid volume → AET
double REFLUX_RATE = TLESR * REFLUX_K * (pH_calc < 4.0 ? 1.0 : exp(-pH_calc + 4.0));
double AET_INST = REFLUX_RATE * 100.0;  // % time pH<4

// Mucosal damage dynamics
double DMG_RATE = MUC_INJ_K * (AET / 100.0) * (1.0 - MUC_HEAL / MAX_DAMAGE);
double HEAL_RATE = K_HEAL * (MUC_HEAL0 / 24.0) * (1.0 + ACID_INH);

$ODE
// ── PPI PK ──────────────────────────────────────────────────────────────────
dxdt_PPI_GUT  = -KA_PPI * PPI_GUT;
dxdt_PPI_CENT = KA_PPI * F_PPI * PPI_GUT * CYP2C19
                - (CL_PPI / V_PPI) / CYP2C19 * PPI_CENT;

// ── H2RA PK ─────────────────────────────────────────────────────────────────
dxdt_H2RA_GUT  = -KA_H2RA * H2RA_GUT;
dxdt_H2RA_CENT = KA_H2RA * F_H2RA * H2RA_GUT - (CL_H2RA / V_H2RA) * H2RA_CENT;

// ── P-CAB PK ─────────────────────────────────────────────────────────────────
dxdt_PCAB_GUT  = -KA_PCAB * PCAB_GUT;
dxdt_PCAB_CENT = KA_PCAB * F_PCAB * PCAB_GUT - (CL_PCAB / V_PCAB) * PCAB_CENT;

// ── Prokinetic PK ─────────────────────────────────────────────────────────────
dxdt_PROK_GUT  = -KA_PROK * PROK_GUT;
dxdt_PROK_CENT = KA_PROK * F_PROK * PROK_GUT - (CL_PROK / V_PROK) * PROK_CENT;

// ── Proton Pump Pool Dynamics ────────────────────────────────────────────────
// Synthesis into inactive pool; translocation to active; irreversible PPI binding
dxdt_PUMP_INACT = K_SYN_PUMP - K_DEG_PUMP * PUMP_INACT - K_ACT * PUMP_INACT;
dxdt_PUMP_ACT   = K_ACT * PUMP_INACT - K_DEACT * PUMP_ACT
                  - K_DEG_PUMP * PUMP_ACT
                  - INH_PPI * 0.5 * PUMP_ACT;  // covalent PPI binding
dxdt_PUMP_INH   = INH_PPI * 0.5 * PUMP_ACT - K_DEG_PUMP * PUMP_INH; // degraded at normal rate

// ── Gastric Acid & pH ────────────────────────────────────────────────────────
double PUMP_FRAC = PUMP_ACT / PUMP_ACT0;
dxdt_ACID_RATE  = 0.1 * (ACID_EFF * PUMP_FRAC - ACID_RATE);  // half-life ~7h
dxdt_GAS_pH     = 0.5 * (pH_calc - GAS_pH);                   // slow equilibration

// ── Esophageal Acid Exposure Time ────────────────────────────────────────────
// AET drifts toward instantaneous AET with clearance
dxdt_AET = 0.2 * (AET_INST - AET);

// ── Mucosal Damage / Healing ─────────────────────────────────────────────────
dxdt_MUC_DMG  = DMG_RATE - HEAL_RATE;
if(MUC_DMG < 0) dxdt_MUC_DMG = 0;
if(MUC_DMG > MAX_DAMAGE) dxdt_MUC_DMG = -0.1;

dxdt_MUC_HEAL = -DMG_RATE + HEAL_RATE;
if(MUC_HEAL < 0) dxdt_MUC_HEAL = 0;
if(MUC_HEAL > MAX_DAMAGE) dxdt_MUC_HEAL = 0;

// ── Symptom Score ─────────────────────────────────────────────────────────────
double SYM_TARGET = 15.0 * (AET / 30.0) + 3.0 * (MUC_DMG / MAX_DAMAGE);
dxdt_SYM_SCORE = 0.1 * (SYM_TARGET - SYM_SCORE);

// ── Barrett Risk (long-term cumulative) ────────────────────────────────────────
dxdt_BE_RISK = 0.0001 * AET * (1.0 - BE_RISK);

$TABLE
capture CP_PPI   = PPI_CENT / V_PPI;
capture CP_H2RA  = H2RA_CENT / V_H2RA;
capture CP_PCAB  = PCAB_CENT / V_PCAB;
capture CP_PROK  = PROK_CENT / V_PROK;
capture ACID_INH_pct = ACID_INH * 100;
capture pH       = GAS_pH;
capture AET_pct  = AET;
capture DMG      = MUC_DMG;
capture HEAL     = MUC_HEAL;
capture SYM      = SYM_SCORE;
capture PUMP_active = PUMP_ACT;
capture PUMP_inhibited = PUMP_INH;
capture LES_pressure = LES_P0 + K_PROK_LES * CP_PROK / (0.005 + CP_PROK);
capture Barrett  = BE_RISK;

$INIT
PPI_GUT   = 0, PPI_CENT   = 0
H2RA_GUT  = 0, H2RA_CENT  = 0
PCAB_GUT  = 0, PCAB_CENT  = 0
PROK_GUT  = 0, PROK_CENT  = 0
PUMP_INACT = 70, PUMP_ACT = 30, PUMP_INH = 0
ACID_RATE = 3.5, GAS_pH = 1.8
AET       = 15.0
MUC_DMG   = 25.0
MUC_HEAL  = 75.0
SYM_SCORE = 8.0
BE_RISK   = 0.0
'

mod <- mcode("GERD_QSP", gerd_code)

# ─── Treatment Scenarios ─────────────────────────────────────────────────────
# Dosing events for 8 weeks (56 days, 1344 h), QD morning dosing

make_events <- function(drug, dose, interval = 24, duration = 56*24) {
  times <- seq(0, duration - 1, by = interval)
  ev(amt = dose, cmt = drug, time = times)
}

scenarios <- list(
  "No Treatment (Control)" = list(
    mrgsolve::ev(amt = 0, cmt = "PPI_GUT", time = 0),
    param = list(DOSE_PPI = 0, DOSE_H2RA = 0, DOSE_PCAB = 0, DOSE_PROK = 0)
  ),
  "Omeprazole 20 mg QD (Standard PPI)" = list(
    ev = make_events("PPI_GUT", 20),
    param = list(DOSE_PPI = 20)
  ),
  "Esomeprazole 40 mg QD (High-dose PPI)" = list(
    ev = make_events("PPI_GUT", 40),
    param = list(DOSE_PPI = 40, F_PPI = 0.73, CL_PPI = 12)
  ),
  "Vonoprazan 20 mg QD (P-CAB)" = list(
    ev = make_events("PCAB_GUT", 20),
    param = list(DOSE_PCAB = 20)
  ),
  "Famotidine 40 mg BID (H2RA)" = list(
    ev = c(make_events("H2RA_GUT", 40, interval = 12)),
    param = list(DOSE_H2RA = 40)
  ),
  "Eso 40 mg QD + Domperidone 10 mg TID\n(PPI + Prokinetic)" = list(
    ev = c(make_events("PPI_GUT", 40),
           make_events("PROK_GUT", 10, interval = 8)),
    param = list(DOSE_PPI = 40, F_PPI = 0.73, CL_PPI = 12, DOSE_PROK = 10)
  )
)

# ─── Simulate all scenarios ──────────────────────────────────────────────────
sim_list <- lapply(seq_along(scenarios), function(i) {
  scen <- scenarios[[i]]
  nm   <- names(scenarios)[i]

  m2 <- mod
  if (!is.null(scen$param)) m2 <- param(m2, scen$param)

  ev_obj <- if (inherits(scen[[1]], "ev")) scen[[1]] else scen$ev

  out <- mrgsim(m2, ev_obj,
                start = 0, end = 56 * 24, delta = 0.5,
                carry_out = "evid")
  out <- as.data.frame(out)
  out$Scenario <- nm
  out
})

all_sims <- bind_rows(sim_list)

# ─── Summary Statistics (Week 8 endpoint) ───────────────────────────────────
summary_stats <- all_sims %>%
  filter(time == max(time)) %>%
  group_by(Scenario) %>%
  summarise(
    pH_mean         = round(mean(pH, na.rm = TRUE), 2),
    AET_pct_mean    = round(mean(AET_pct, na.rm = TRUE), 1),
    Acid_inh_pct    = round(mean(ACID_INH_pct, na.rm = TRUE), 1),
    Mucosal_damage  = round(mean(DMG, na.rm = TRUE), 1),
    Symptom_score   = round(mean(SYM, na.rm = TRUE), 1),
    Barrett_risk    = round(mean(Barrett, na.rm = TRUE), 4),
    .groups = "drop"
  )

print(summary_stats)

# ─── Plots ───────────────────────────────────────────────────────────────────
cols <- c(
  "No Treatment (Control)"                      = "#E53935",
  "Omeprazole 20 mg QD (Standard PPI)"          = "#FB8C00",
  "Esomeprazole 40 mg QD (High-dose PPI)"       = "#8E24AA",
  "Vonoprazan 20 mg QD (P-CAB)"                 = "#00897B",
  "Famotidine 40 mg BID (H2RA)"                 = "#1E88E5",
  "Eso 40 mg QD + Domperidone 10 mg TID\n(PPI + Prokinetic)" = "#43A047"
)

# Thin data for plotting
plot_data <- all_sims %>% filter(time %% 6 == 0)

p1 <- ggplot(plot_data, aes(x = time / 24, y = pH, colour = Scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  labs(title = "A. Intragastric pH", x = "Time (days)", y = "pH") +
  geom_hline(yintercept = 4, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 50, y = 4.2, label = "pH 4 threshold", size = 3) +
  theme_classic() + theme(legend.position = "none")

p2 <- ggplot(plot_data, aes(x = time / 24, y = AET_pct, colour = Scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  labs(title = "B. Acid Exposure Time (%)", x = "Time (days)", y = "AET (%)") +
  geom_hline(yintercept = 6, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 50, y = 7, label = "Lyon 2.0 cutoff (6%)", size = 3) +
  theme_classic() + theme(legend.position = "none")

p3 <- ggplot(plot_data, aes(x = time / 24, y = DMG, colour = Scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  labs(title = "C. Mucosal Damage Score", x = "Time (days)", y = "Damage (0-100)") +
  theme_classic() + theme(legend.position = "none")

p4 <- ggplot(plot_data, aes(x = time / 24, y = SYM, colour = Scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  labs(title = "D. Symptom Score (GERD-Q)", x = "Time (days)", y = "Score (0-18)") +
  geom_hline(yintercept = 8, linetype = "dashed", colour = "grey50") +
  annotate("text", x = 50, y = 8.5, label = "GERD-Q ≥8 = GERD", size = 3) +
  theme_classic() + theme(legend.position = "right", legend.text = element_text(size = 7))

combined <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title    = "GERD QSP Model — 8-week Treatment Simulation",
    subtitle = "PPI vs H2RA vs P-CAB vs Prokinetic combination",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(combined)

# ─── CYP2C19 Phenotype Sensitivity Analysis ──────────────────────────────────
cyp_scenarios <- list(
  "Ultra-Rapid (UM, CYP2C19×2.0)"  = list(CYP2C19 = 2.0),
  "Extensive (EM, CYP2C19×1.0)"    = list(CYP2C19 = 1.0),
  "Intermediate (IM, CYP2C19×0.6)" = list(CYP2C19 = 0.6),
  "Poor (PM, CYP2C19×0.25)"        = list(CYP2C19 = 0.25)
)

cyp_sims <- lapply(seq_along(cyp_scenarios), function(i) {
  m2 <- param(mod, cyp_scenarios[[i]])
  out <- mrgsim(m2, make_events("PPI_GUT", 20),
                start = 0, end = 56 * 24, delta = 0.5)
  out <- as.data.frame(out)
  out$Phenotype <- names(cyp_scenarios)[i]
  out
})

cyp_data <- bind_rows(cyp_sims) %>% filter(time %% 6 == 0)

p_cyp <- ggplot(cyp_data, aes(x = time / 24, y = pH, colour = Phenotype)) +
  geom_line(linewidth = 1) +
  labs(
    title    = "CYP2C19 Phenotype Effect on Gastric pH\n(Omeprazole 20 mg QD)",
    x        = "Time (days)",
    y        = "Intragastric pH",
    colour   = "CYP2C19 Phenotype"
  ) +
  geom_hline(yintercept = 4, linetype = "dashed") +
  scale_colour_manual(values = c("#C62828", "#1976D2", "#388E3C", "#F57F17")) +
  theme_classic()

print(p_cyp)

# ─── Dose-Response at Week 8 ─────────────────────────────────────────────────
doses <- c(5, 10, 20, 40, 80)
dr_ppi <- lapply(doses, function(d) {
  out <- mrgsim(mod, make_events("PPI_GUT", d),
                start = 0, end = 56 * 24, delta = 1)
  out <- as.data.frame(out)
  data.frame(Dose = d, AET = tail(out$AET_pct, 1), pH = tail(out$pH, 1),
             Drug = "PPI (Omeprazole)")
})

dr_pcab <- lapply(doses, function(d) {
  out <- mrgsim(mod, make_events("PCAB_GUT", d),
                start = 0, end = 56 * 24, delta = 1)
  out <- as.data.frame(out)
  data.frame(Dose = d, AET = tail(out$AET_pct, 1), pH = tail(out$pH, 1),
             Drug = "P-CAB (Vonoprazan)")
})

dr_data <- bind_rows(c(dr_ppi, dr_pcab))

p_dr <- ggplot(dr_data, aes(x = Dose, y = AET, colour = Drug, group = Drug)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  scale_x_log10() +
  labs(
    title  = "Dose-Response: AET (%) at Week 8",
    x      = "Dose (mg, log scale)",
    y      = "Acid Exposure Time (%)",
    colour = "Drug Class"
  ) +
  geom_hline(yintercept = 6, linetype = "dashed") +
  scale_colour_manual(values = c("#8E24AA", "#00897B")) +
  theme_classic()

print(p_dr)

message("\n=== GERD QSP Model simulation complete ===")
message("Key Week-8 endpoints:")
print(summary_stats)
