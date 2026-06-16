##############################################################################
#  QSP Model — Pulmonary Arterial Hypertension (PAH)
#  Framework : mrgsolve (R)
#  Author    : Auto-generated QSP scaffold
#  Date      : 2026-06-16
#  Reference : See pah_references.md for full citation list
##############################################################################
#
#  Model structure (compartments):
#   PK  — ERA (endothelin receptor antagonist)
#          PDE5i (sildenafil/tadalafil)
#          PGI2 (prostacyclin analogue)
#   PD  — ET-1 dynamics
#          cGMP / NO signalling
#          cAMP / PGI2 signalling
#          PASMC proliferation / apoptosis index (remodelling score)
#          Right-ventricular function: PVR, mPAP, CO, BNP
#          Clinical endpoints: 6MWD, WHO-FC logistic
#
##############################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ── mrgsolve model code string ──────────────────────────────────────────────
pah_model_code <- '
$PROB PAH QSP Model — mrgsolve implementation

$PLUGIN Rcpp

##############################################################################
# PARAMETERS (all have literature-derived defaults; can be overridden)
##############################################################################
$PARAM
// ── ERA PK (bosentan-like; 2-compartment oral) ───────────────────────────
ERA_dose   = 125,   // mg, twice daily
ERA_F      = 0.50,  // oral bioavailability
ERA_ka     = 0.693, // h-1  absorption rate  (t½abs ~ 1 h)
ERA_CL     = 15,    // L/h  apparent clearance
ERA_Vc     = 35,    // L    central volume
ERA_Q      = 9,     // L/h  inter-compartment clearance
ERA_Vp     = 50,    // L    peripheral volume
ERA_ke0    = 0.3,   // h-1  effect-site equilibration
ERA_IC50   = 1.0,   // ng/mL  inhibitory concentration (ETA)
ERA_Emax   = 0.85,  // max fractional ETA block
ERA_n_hill = 1.5,   // Hill coefficient

// ── PDE5i PK (tadalafil-like; 1-compartment oral) ────────────────────────
PDE5_dose  = 40,    // mg, once daily
PDE5_F     = 0.80,
PDE5_ka    = 0.35,  // h-1
PDE5_CL    = 3.5,   // L/h
PDE5_Vc    = 63,    // L
PDE5_ke0   = 0.2,   // h-1
PDE5_IC50  = 0.5,   // ng/mL
PDE5_Emax  = 0.90,
PDE5_n     = 1.2,

// ── PGI2 PK (epoprostenol-like; IV infusion, t½ ~3 min) ─────────────────
PGI2_rate  = 2.0,   // ng/kg/min (IV infusion rate, continuous)
PGI2_BW    = 70,    // kg body weight
PGI2_CL    = 800,   // L/h  (very high clearance)
PGI2_Vc    = 30,    // L
PGI2_ke0   = 2.0,   // h-1
PGI2_EC50  = 0.1,   // ng/mL
PGI2_Emax  = 0.95,
PGI2_n     = 1.0,

// ── ET-1 PD dynamics ─────────────────────────────────────────────────────
ET1_baseline = 3.0,   // pg/mL  healthy reference
ET1_PAH      = 9.0,   // pg/mL  typical PAH baseline
ET1_kprod    = 0.15,  // h-1    ET-1 production rate constant
ET1_kdeg     = 0.15,  // h-1    ET-1 clearance rate constant
ET1_stimHyp  = 1.8,   // fold increase in production due to hypoxia/PAH state
Kd_ETA       = 1.0,   // dimensionless ETA fractional activation at baseline

// ── NO-cGMP dynamics ─────────────────────────────────────────────────────
cGMP_baseline = 5.0,  // nM (arbitrary units)
cGMP_kprod    = 1.0,  // production rate (stimulated by NO)
cGMP_kdeg_0   = 0.8,  // h-1 basal PDE5 degradation rate
NO_synth_eff  = 0.6,  // fractional eNOS efficiency in PAH (reduced)
sGC_baseline  = 1.0,  // sGC activity (normalized)

// ── cAMP dynamics ────────────────────────────────────────────────────────
cAMP_baseline = 4.0,  // nM
cAMP_kprod    = 0.9,
cAMP_kdeg     = 0.7,
IP_r_eff      = 1.0,  // IP receptor coupling efficiency

// ── Vascular Remodelling Index (VRI) ─────────────────────────────────────
// Lumped ODE for PASMC proliferation vs apoptosis balance
// VRI = 1 at disease onset, 0 = normal, increases over time without Rx
VRI_0         = 1.0,   // initial remodelling severity (PAH onset)
VRI_kgrowth   = 0.015, // h-1  net remodelling progression rate
VRI_kdrug_ERA = 0.30,  // h-1  ERA contribution to remodelling reversal
VRI_kdrug_PDE = 0.20,
VRI_kdrug_PGI = 0.45,
VRI_max       = 3.0,   // maximum remodelling index (plexiform lesion)

// ── Haemodynamic model ────────────────────────────────────────────────────
// PVR (dyn·s·cm-5) modeled as function of tone + remodelling
PVR_normal    = 100,   // normal PVR dyn·s·cm-5
PVR_PAH0      = 900,   // baseline PAH PVR
PVR_tone_frac = 0.45,  // fraction of PVR from vascular tone (reversible)
PVR_remod_frac= 0.55,  // fraction from remodelling (partially reversible)

// mPAP = CO * PVR/80 + PAWP  (simplified Ohm's law)
PAWP_baseline = 10,    // mmHg  pulmonary arterial wedge pressure
CO_baseline   = 3.5,   // L/min  baseline cardiac output in PAH
CO_normal     = 5.5,

// ── RV function model ─────────────────────────────────────────────────────
Ees_0         = 0.8,   // mmHg/mL  RV end-systolic elastance (PAH onset)
Ees_max       = 2.5,   // max Ees (compensatory hypertrophy)
Ees_kHyp      = 0.003, // h-1  rate of RV hypertrophy
Ees_kFail     = 0.006, // h-1  rate of RV decompensation
RV_fail_thresh= 1.5,   // Ea threshold for decompensation (mmHg/mL)

// ── BNP dynamics ──────────────────────────────────────────────────────────
BNP_normal    = 20,    // pg/mL
BNP_kprod     = 0.1,   // h-1  production proportional to wall stress
BNP_kdeg      = 0.08,  // h-1

// ── 6MWD model ────────────────────────────────────────────────────────────
sixMWD_normal = 580,   // m    healthy reference
sixMWD_PAH0   = 330,   // m    typical PAH baseline
// 6MWD decreases with ↑PVR and ↑BNP, increases with treatment

$CMT
// PK compartments
ERA_gut ERA_central ERA_periph ERA_effect
PDE5_gut PDE5_central PDE5_effect
PGI2_iv PGI2_central PGI2_effect

// PD / Biomarker compartments
ET1       // pg/mL
cGMP      // nM
cAMP      // nM
VRI       // vascular remodelling index (dimensionless)
Ees_RV    // RV end-systolic elastance (mmHg/mL)
BNP_conc  // pg/mL

$MAIN
// ── Steady-state initial conditions ──
ET1_0    = ET1_PAH;
cGMP_0   = cGMP_baseline;
cAMP_0   = cAMP_baseline;
VRI_0_ic = VRI_0;
Ees_RV_0 = Ees_0;
BNP_0    = 200;     // typical PAH BNP at onset

if(NEWIND <=1) {
    ERA_gut_0    = 0;
    ERA_central_0= 0;
    ERA_periph_0 = 0;
    ERA_effect_0 = 0;
    PDE5_gut_0   = 0;
    PDE5_central_0=0;
    PDE5_effect_0= 0;
    PGI2_iv_0    = 0;
    PGI2_central_0=0;
    PGI2_effect_0= 0;
}

$ODE
// ─────────────────────────────────────────────────────────────────────────
//  SECTION A: ERA PK (2-compartment, oral)
// ─────────────────────────────────────────────────────────────────────────
double ERA_abs   = ERA_ka * ERA_gut;
double ERA_distr = ERA_Q  * (ERA_central/ERA_Vc - ERA_periph/ERA_Vp);
double ERA_elim  = (ERA_CL/ERA_Vc) * ERA_central;
double ERA_eq    = ERA_ke0 * (ERA_central/ERA_Vc - ERA_effect);

dxdt_ERA_gut     = -ERA_abs;
dxdt_ERA_central =  ERA_abs * ERA_F - ERA_distr - ERA_elim;
dxdt_ERA_periph  =  ERA_distr;
dxdt_ERA_effect  =  ERA_eq;

// ─────────────────────────────────────────────────────────────────────────
//  SECTION B: PDE5i PK (1-compartment, oral)
// ─────────────────────────────────────────────────────────────────────────
double PDE5_abs  = PDE5_ka * PDE5_gut;
double PDE5_elim = (PDE5_CL/PDE5_Vc) * PDE5_central;
double PDE5_eq   = PDE5_ke0 * (PDE5_central/PDE5_Vc - PDE5_effect);

dxdt_PDE5_gut     = -PDE5_abs;
dxdt_PDE5_central =  PDE5_abs * PDE5_F - PDE5_elim;
dxdt_PDE5_effect  =  PDE5_eq;

// ─────────────────────────────────────────────────────────────────────────
//  SECTION C: PGI2 PK (1-compartment IV, very rapid t½)
// ─────────────────────────────────────────────────────────────────────────
// PGI2_iv tracks administered dose (infusion handled via $EVENT)
double PGI2_abs  = PGI2_ke0 * PGI2_central;   // rapid equilibration
double PGI2_elim = (PGI2_CL/PGI2_Vc) * PGI2_central;
double PGI2_eq   = PGI2_ke0 * (PGI2_central/PGI2_Vc - PGI2_effect);

dxdt_PGI2_iv      = 0;   // controlled externally (RATE infusion)
dxdt_PGI2_central = PGI2_iv - PGI2_elim;
dxdt_PGI2_effect  = PGI2_eq;

// ─────────────────────────────────────────────────────────────────────────
//  SECTION D: PD DRUG EFFECTS (Emax models)
// ─────────────────────────────────────────────────────────────────────────
// ERA — block ETA receptor (Hill equation)
double ERA_Ce   = ERA_effect;
double ERA_Inh  = ERA_Emax * pow(ERA_Ce, ERA_n_hill) /
                  (pow(ERA_IC50, ERA_n_hill) + pow(ERA_Ce, ERA_n_hill));

// PDE5i — inhibit cGMP degradation
double PDE5_Ce  = PDE5_effect;
double PDE5_Inh = PDE5_Emax * pow(PDE5_Ce, PDE5_n) /
                  (pow(PDE5_IC50, PDE5_n) + pow(PDE5_Ce, PDE5_n));

// PGI2 — activate IP receptor → cAMP
double PGI2_Ce  = PGI2_effect;
double PGI2_Act = PGI2_Emax * PGI2_Ce / (PGI2_EC50 + PGI2_Ce);

// ─────────────────────────────────────────────────────────────────────────
//  SECTION E: ET-1 DYNAMICS
// ─────────────────────────────────────────────────────────────────────────
// ET-1 production: elevated by PAH state (hypoxia, inflammation)
// ERA blocks ETA receptor (does not directly reduce ET-1 production,
// but via ETB clearance effect; simplification: ERA slightly ↓ ET-1)
double ET1_prod = ET1_kprod * ET1_stimHyp * (1.0 - 0.15 * ERA_Inh);
double ET1_clear= ET1_kdeg;
dxdt_ET1 = ET1_prod * (ET1_PAH - ET1) - ET1_clear * ET1;
// Note: This is a turnover model. At steady state ET1 = ET1_PAH in absence of Rx

// ─────────────────────────────────────────────────────────────────────────
//  SECTION F: cGMP DYNAMICS (NO-sGC-PDE5 axis)
// ─────────────────────────────────────────────────────────────────────────
// cGMP production: reduced eNOS in PAH, stimulated by PDE5i/riociguat
// cGMP degradation: accelerated by PDE5 (countered by PDE5i)
double cGMP_prod = cGMP_kprod * NO_synth_eff * sGC_baseline;
double PDE5_eff_rate = cGMP_kdeg_0 * (1.0 - PDE5_Inh);  // PDE5 inhibited
dxdt_cGMP = cGMP_prod - PDE5_eff_rate * cGMP;

// ─────────────────────────────────────────────────────────────────────────
//  SECTION G: cAMP DYNAMICS (PGI2-IP-AC axis)
// ─────────────────────────────────────────────────────────────────────────
double cAMP_prod = cAMP_kprod * IP_r_eff * (1.0 + PGI2_Act);
dxdt_cAMP = cAMP_prod - cAMP_kdeg * cAMP;

// ─────────────────────────────────────────────────────────────────────────
//  SECTION H: VASCULAR REMODELLING INDEX (VRI)
// ─────────────────────────────────────────────────────────────────────────
// VRI grows over time (disease progression) and is reduced by treatments
// VRI drives structural PVR increase
// Drug effects on remodelling (slower time-scale, weeks-months)
double VRI_drug_effect = VRI_kdrug_ERA * ERA_Inh +
                         VRI_kdrug_PDE * (cGMP/cGMP_baseline - 1.0) +
                         VRI_kdrug_PGI * PGI2_Act;
double VRI_net = VRI_kgrowth * VRI * (1.0 - VRI/VRI_max) - VRI_drug_effect * VRI;
dxdt_VRI = VRI_net;

// ─────────────────────────────────────────────────────────────────────────
//  SECTION I: RIGHT VENTRICULAR ELASTANCE (adaptive → maladaptive)
// ─────────────────────────────────────────────────────────────────────────
// Ea = effective arterial elastance = mPAP_sys / SV
// Ees initially increases (hypertrophy), then declines (failure)
// Simplified: Ees rate depends on Ea relative to threshold
double PVR_curr = PVR_normal + (PVR_PAH0 - PVR_normal) *
                  (PVR_tone_frac * (ET1/ET1_PAH) * (1.0 - ERA_Inh) *
                   (cGMP_baseline/cGMP) * (cAMP_baseline/cAMP) +
                   PVR_remod_frac * (VRI/VRI_0_ic));
double CO_curr  = CO_baseline * (PVR_PAH0 / PVR_curr);  // simplified inverse
double mPAP_curr= CO_curr * PVR_curr / 80.0 + PAWP_baseline;
double Ea_curr  = mPAP_curr / (CO_curr * 1000.0 / 60.0);  // mmHg/mL

// Ees changes: hypertrophy if Ea < threshold, failure if Ea > threshold
double dEes;
if(Ea_curr < RV_fail_thresh) {
    dEes = Ees_kHyp * (Ees_max - Ees_RV);  // adaptive hypertrophy
} else {
    dEes = -Ees_kFail * (Ees_RV - Ees_0 * 0.5);  // decompensation
}
dxdt_Ees_RV = dEes;

// ─────────────────────────────────────────────────────────────────────────
//  SECTION J: BNP DYNAMICS
// ─────────────────────────────────────────────────────────────────────────
// BNP production proportional to RV wall stress ~ mPAP × RV volume
double BNP_stim = BNP_kprod * (mPAP_curr / 25.0);   // normalised to mPAP
dxdt_BNP_conc = BNP_stim * BNP_normal - BNP_kdeg * BNP_conc;

$TABLE
// ─────────────────────────────────────────────────────────────────────────
//  DERIVED VARIABLES (for output)
// ─────────────────────────────────────────────────────────────────────────
// Drug effect fractions
double ERA_effect_frac = ERA_Emax * pow(ERA_effect, ERA_n_hill) /
    (pow(ERA_IC50, ERA_n_hill) + pow(ERA_effect, ERA_n_hill));
double PDE5_effect_frac = PDE5_Emax * pow(PDE5_effect, PDE5_n) /
    (pow(PDE5_IC50, PDE5_n) + pow(PDE5_effect, PDE5_n));
double PGI2_effect_frac = PGI2_Emax * PGI2_effect / (PGI2_EC50 + PGI2_effect);

// PVR (dyn·s·cm-5)
double PVR_sim = PVR_normal + (PVR_PAH0 - PVR_normal) *
    (PVR_tone_frac * (ET1/ET1_PAH) * (1.0 - ERA_effect_frac) *
     (cGMP_baseline/cGMP) * (cAMP_baseline/cAMP) +
     PVR_remod_frac * (VRI/VRI_0));

// mPAP (mmHg)
double CO_sim   = CO_baseline * (PVR_PAH0 / PVR_sim);
double mPAP_sim = CO_sim * PVR_sim / 80.0 + PAWP_baseline;

// 6-minute walk distance (m) — linear function of PVR & BNP
double sixMWD_sim = sixMWD_normal -
    (PVR_sim - PVR_normal) * (sixMWD_normal - sixMWD_PAH0) /
    (PVR_PAH0 - PVR_normal);
if(sixMWD_sim < 50)  sixMWD_sim = 50;
if(sixMWD_sim > sixMWD_normal) sixMWD_sim = sixMWD_normal;

// WHO Functional Class (1-4, approximated as continuous)
// FC driven by 6MWD: FC IV < 150m, FC III 150-300, FC II 300-440, FC I > 440
double WHO_FC_sim;
if(sixMWD_sim > 440)      WHO_FC_sim = 1.0;
else if(sixMWD_sim > 300) WHO_FC_sim = 1.0 + (440-sixMWD_sim)/140.0;
else if(sixMWD_sim > 150) WHO_FC_sim = 2.0 + (300-sixMWD_sim)/150.0;
else                       WHO_FC_sim = 3.0 + (150-sixMWD_sim)/150.0;
if(WHO_FC_sim > 4.0) WHO_FC_sim = 4.0;

// RV-PA coupling (Ees/Ea)
double Ea_sim = mPAP_sim / (CO_sim * 1000.0/60.0);
double coupling_ratio = Ees_RV / Ea_sim;

capture PVR_dyn  = PVR_sim;
capture mPAP_mmHg= mPAP_sim;
capture CO_Lmin  = CO_sim;
capture sixMWD_m = sixMWD_sim;
capture WHO_FC   = WHO_FC_sim;
capture BNP_pg   = BNP_conc;
capture ET1_pg   = ET1;
capture cGMP_nM  = cGMP;
capture cAMP_nM  = cAMP;
capture VRI_idx  = VRI;
capture Ees_mmHg = Ees_RV;
capture Ea_mmHg  = Ea_sim;
capture RV_PA_coupling = coupling_ratio;
capture ERA_Ce_ng   = ERA_effect;
capture PDE5_Ce_ng  = PDE5_effect;
capture PGI2_Ce_ng  = PGI2_effect;
capture ERA_Inh_frac= ERA_effect_frac;
capture PDE5_Inh_frac=PDE5_effect_frac;
capture PGI2_Act_frac=PGI2_effect_frac;
'

# ── Compile the model ────────────────────────────────────────────────────────
mod <- mread_cache("pah_qsp", tempdir(), pah_model_code)

cat("Model compiled successfully. Compartments:", mod@cmtL, "\n")

# ─────────────────────────────────────────────────────────────────────────────
#  DOSING EVENTS
# ─────────────────────────────────────────────────────────────────────────────
make_dosing <- function(
    era_dose   = 125,   # mg BID
    pde5_dose  = 40,    # mg QD
    pgi2_rate  = 2.0,   # ng/kg/min (0 = no PGI2)
    bw         = 70,    # kg
    duration   = 12 * 7 * 24  # 12 weeks in hours
) {
  ev_list <- list()

  # ERA: twice daily oral
  if (era_dose > 0) {
    ev_list$era <- ev(
      cmt  = "ERA_gut",
      amt  = era_dose * 1000,  # convert mg → μg (adjust for unit consistency)
      ii   = 12,
      addl = as.integer(duration / 12) - 1
    )
  }

  # PDE5i: once daily oral
  if (pde5_dose > 0) {
    ev_list$pde5 <- ev(
      cmt  = "PDE5_gut",
      amt  = pde5_dose * 1000,
      ii   = 24,
      addl = as.integer(duration / 24) - 1
    )
  }

  # PGI2: continuous IV infusion (rate in total ng/h)
  if (pgi2_rate > 0) {
    total_rate_ng_h <- pgi2_rate * bw * 60  # ng/kg/min → ng/h
    ev_list$pgi2 <- ev(
      cmt  = "PGI2_central",
      amt  = total_rate_ng_h * duration,
      rate = total_rate_ng_h,
      time = 0
    )
  }

  if (length(ev_list) == 0) return(ev(time = 0, amt = 0, cmt = 1))
  Reduce(c, ev_list)
}

# ─────────────────────────────────────────────────────────────────────────────
#  SIMULATION SCENARIOS
# ─────────────────────────────────────────────────────────────────────────────
sim_duration <- 12 * 7 * 24  # 12 weeks = 2016 hours
sim_times    <- c(seq(0, 168, by = 1),        # first week: hourly
                  seq(168, sim_duration, by = 24))  # rest: daily

scenarios <- list(
  "No Treatment"              = list(era=0,   pde5=0,  pgi2=0),
  "ERA monotherapy"           = list(era=125, pde5=0,  pgi2=0),
  "PDE5i monotherapy"         = list(era=0,   pde5=40, pgi2=0),
  "ERA + PDE5i combination"   = list(era=125, pde5=40, pgi2=0),
  "PGI2 monotherapy"          = list(era=0,   pde5=0,  pgi2=2),
  "Triple therapy"            = list(era=125, pde5=40, pgi2=2)
)

results_list <- lapply(names(scenarios), function(scen_name) {
  sc <- scenarios[[scen_name]]
  dose_ev <- make_dosing(era_dose  = sc$era,
                         pde5_dose = sc$pde5,
                         pgi2_rate = sc$pgi2,
                         duration  = sim_duration)
  out <- mod %>%
    init(ET1      = 9.0,
         cGMP     = 5.0,
         cAMP     = 4.0,
         VRI      = 1.0,
         Ees_RV   = 0.8,
         BNP_conc = 200) %>%
    ev(dose_ev) %>%
    mrgsim(end = sim_duration, delta = 1, obsonly = TRUE) %>%
    as.data.frame() %>%
    mutate(scenario = scen_name,
           time_days = time / 24)
  out
})

results <- bind_rows(results_list)

# ─────────────────────────────────────────────────────────────────────────────
#  PLOTTING
# ─────────────────────────────────────────────────────────────────────────────
scenarios_ordered <- names(scenarios)
results$scenario <- factor(results$scenario, levels = scenarios_ordered)

pal <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628")

# Keep daily snapshots for plotting
plot_data <- results %>% filter(time %% 24 < 1)

p1 <- ggplot(plot_data, aes(time_days, PVR_dyn, colour = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = pal) +
  labs(title = "A — Pulmonary Vascular Resistance over 12 weeks",
       x = "Time (days)", y = "PVR (dyn·s·cm⁻⁵)", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

p2 <- ggplot(plot_data, aes(time_days, mPAP_mmHg, colour = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = pal) +
  geom_hline(yintercept = 25, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 5, y = 26, label = "PAH threshold (25 mmHg)", size = 3) +
  labs(title = "B — Mean Pulmonary Artery Pressure (mPAP)",
       x = "Time (days)", y = "mPAP (mmHg)", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p3 <- ggplot(plot_data, aes(time_days, sixMWD_m, colour = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = pal) +
  labs(title = "C — 6-Minute Walk Distance (6MWD)",
       x = "Time (days)", y = "6MWD (m)", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p4 <- ggplot(plot_data, aes(time_days, BNP_pg, colour = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = pal) +
  labs(title = "D — BNP / NT-proBNP",
       x = "Time (days)", y = "BNP (pg/mL)", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p5 <- ggplot(plot_data, aes(time_days, RV_PA_coupling, colour = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = pal) +
  geom_hline(yintercept = 1.0, linetype = "dashed", colour = "red") +
  annotate("text", x = 5, y = 1.05, label = "Ees/Ea = 1\n(decompensation threshold)",
           size = 3, colour = "red") +
  labs(title = "E — RV-PA Coupling (Ees/Ea)",
       x = "Time (days)", y = "Ees/Ea (dimensionless)", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p6 <- ggplot(plot_data, aes(time_days, VRI_idx, colour = scenario)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = pal) +
  labs(title = "F — Vascular Remodelling Index",
       x = "Time (days)", y = "VRI (dimensionless)", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# Arrange all panels
library(patchwork)
combined_plot <- (p1 + p2) / (p3 + p4) / (p5 + p6) +
  plot_annotation(
    title   = "PAH QSP Model — 12-Week Simulation Across Treatment Scenarios",
    subtitle= "Baseline: mPAP 45 mmHg, PVR 900 dyn·s·cm⁻⁵, 6MWD 330 m",
    caption = "ERA = endothelin receptor antagonist | PDE5i = PDE-5 inhibitor | PGI2 = prostacyclin analogue"
  )

print(combined_plot)

# ─────────────────────────────────────────────────────────────────────────────
#  DOSE-RESPONSE ANALYSIS (ERA IC50 sensitivity)
# ─────────────────────────────────────────────────────────────────────────────
dose_levels_era  <- c(0, 31.25, 62.5, 125, 250)
dose_levels_pde5 <- c(0, 10, 20, 40, 80)

dr_results <- lapply(dose_levels_era, function(d) {
  dose_ev <- make_dosing(era_dose = d, pde5_dose = 0, pgi2_rate = 0,
                         duration = sim_duration)
  mod %>%
    init(ET1 = 9, cGMP = 5, cAMP = 4, VRI = 1, Ees_RV = 0.8, BNP_conc = 200) %>%
    ev(dose_ev) %>%
    mrgsim(end = sim_duration, delta = 24, obsonly = TRUE) %>%
    as.data.frame() %>%
    filter(time == max(time)) %>%
    mutate(ERA_dose_mg = d)
})

dr_df <- bind_rows(dr_results)

p_dr <- ggplot(dr_df, aes(ERA_dose_mg, sixMWD_m)) +
  geom_point(size = 3, colour = "#377EB8") +
  geom_line(colour = "#377EB8") +
  scale_x_continuous(breaks = dose_levels_era) +
  labs(title = "ERA Dose-Response: 6MWD at 12 Weeks",
       x = "ERA dose (mg BID)",
       y = "6MWD at week 12 (m)") +
  theme_bw(base_size = 11)

print(p_dr)

# ─────────────────────────────────────────────────────────────────────────────
#  SUMMARY TABLE AT WEEK 12
# ─────────────────────────────────────────────────────────────────────────────
summary_tbl <- results %>%
  filter(time == sim_duration) %>%
  select(scenario, PVR_dyn, mPAP_mmHg, sixMWD_m, BNP_pg, WHO_FC,
         RV_PA_coupling, VRI_idx) %>%
  mutate(across(where(is.numeric), ~round(.x, 1)))

cat("\n=== PAH QSP Model — Week 12 Summary ===\n")
print(summary_tbl, n = Inf)

# Return compiled model for Shiny use
invisible(list(model = mod, results = results, summary = summary_tbl))
