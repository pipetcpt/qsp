## ============================================================
## Raynaud's Phenomenon (Primary & Secondary) — QSP mrgsolve Model
## Mechanistic ODE: Adrenergic · ET-1 · NO/cGMP · RhoA/ROCK
##                  PGI2/cAMP · Platelet · Neuropeptide Pathways
## Drugs: Nifedipine (CCB), Sildenafil (PDE5i), Bosentan (ERA),
##         Iloprost (PGI2 analog), Prazosin (α1-blocker)
## Compartments: 18 ODE states
## Reference calibration:
##   - Wigley FM et al. N Engl J Med 2002 (Raynaud epidemiology & pathophysiology)
##   - Matucci-Cerinic M et al. Arthritis Rheum 2011 (bosentan DU trial)
##   - Fries R et al. J Am Coll Cardiol 2005 (sildenafil Raynaud)
##   - Thompson AE & Pope JE, Rheumatology 2005 (CCB meta-analysis)
##   - Belch JJ et al. Ann Rheum Dis 2008 (iloprost RCT)
## ============================================================

[PROB]
Raynaud's Phenomenon QSP Model
Simulates vasospasm frequency, digital blood flow, RCS over 12 weeks
Primary (idiopathic) and secondary (SSc-associated) subtypes

[PARAM]
// --- Sympathetic / Adrenergic ---
k_NE_base  = 0.05    // basal NE tone (nmol/L/h)
k_cold_NE  = 0.80    // cold stimulus amplitude on NE
alpha2_sens = 1.20   // α2-AR sensitivity multiplier (1.0=normal, 1.5=Raynaud)
k_alpha2_constr = 0.30  // α2→vasoconstriction rate constant
k_NE_clear = 0.50    // NE clearance (h-1)
kRho_NE    = 0.15    // NE → RhoA activation

// --- VSMC Calcium / Contraction ---
Cai_base   = 0.10    // basal [Ca2+]i (μM)
k_VGCC_open = 0.25   // VGCC opening rate
k_Ca_clear  = 0.40   // Ca2+ clearance rate
EC50_Ca_contr = 0.25  // EC50 for Ca2+-driven contraction
Emax_Ca_contr = 1.0   // maximal contraction from Ca

// --- RhoA / ROCK ---
RhoA_base  = 0.10    // basal RhoA-GTP (AU)
k_ROCK_act = 0.50    // RhoA → ROCK activation
k_ROCK_clear = 0.30  // ROCK clearance
EC50_ROCK  = 0.30    // ROCK → MLC phosphorylation EC50

// --- NO / cGMP pathway ---
eNOS_base  = 1.0     // basal eNOS activity (AU; ↓0.6 in secondary)
k_eNOS_NO  = 0.40    // eNOS → NO production
k_NO_decay = 0.60    // NO half-life (h-1; short: ~2 min → ~0.6 h-1)
k_sGC_cGMP = 0.35    // NO → cGMP via sGC
k_cGMP_clear = 0.25  // PDE5 cGMP clearance
EC50_cGMP  = 0.20    // cGMP → MLCP EC50

// --- PGI2 / cAMP ---
PGI2_base  = 0.50    // basal PGI2 (pg/mL; ↓ in Raynaud)
k_IP_cAMP  = 0.30    // PGI2 → cAMP
k_cAMP_clear = 0.40  // PDE cAMP clearance
EC50_cAMP  = 0.25    // cAMP → VSMC relaxation EC50

// --- Endothelin-1 ---
ET1_base   = 0.80    // basal plasma ET-1 (pg/mL; ↑ in secondary)
k_ET1_synth = 0.05   // basal synthesis
k_ET1_clear = 0.10   // ET-1 clearance (h-1)
k_ETA_constr = 0.40  // ETA → vasoconstriction
k_ETA_Rho  = 0.20    // ETA → RhoA activation

// --- Digital blood flow (outcome) ---
DBF_max    = 10.0    // max digital blood flow (mL/min/100g)
DBF_base   = 5.0     // resting DBF
k_DBF_contr = 0.60   // contraction → DBF decrease
k_DBF_relax = 0.40   // relaxation → DBF recovery

// --- Vasospasm Episode threshold ---
VasoThresh  = 0.55   // VSMC contraction threshold for attack
AttackDecay = 0.80   // attack resolution rate (h-1)
EpiFreq_base = 3.0   // basal episodes/week (untreated primary)

// --- ROS / Inflammation ---
ROS_base   = 0.20    // basal ROS (AU)
k_ROS_prod = 0.15    // ischemia/NE → ROS
k_ROS_decay = 0.35   // antioxidant clearance
k_NF_kB    = 0.12    // ROS → NF-κB
k_ET1_inf  = 0.08    // NF-κB → ET-1 (inflammation amplifies)

// --- PK parameters: Nifedipine (CCB extended-release) ---
ka_NIF = 0.50        // absorption rate (h-1)
F_NIF  = 0.85        // bioavailability
Vd_NIF = 120.0       // volume of distribution (L)
CL_NIF = 60.0        // clearance (L/h)

// --- PK parameters: Sildenafil ---
ka_SIL = 0.80
F_SIL  = 0.40
Vd_SIL = 105.0
CL_SIL = 41.0

// --- PK parameters: Bosentan ---
ka_BOS = 0.60
F_BOS  = 0.50
Vd_BOS = 18.0
CL_BOS = 4.0

// --- PK parameters: Iloprost (IV infusion—direct input) ---
Vd_ILO = 25.0
CL_ILO = 15.0

// --- PK parameters: Prazosin ---
ka_PRA = 1.20
F_PRA  = 0.68
Vd_PRA = 97.0
CL_PRA = 40.0

// --- PD: Drug effect parameters ---
Emax_CCB  = 0.75     // max fractional VGCC block by CCB
EC50_CCB  = 15.0     // nifedipine EC50 (ng/mL)
Emax_PDE5 = 0.80     // max cGMP elevation by PDE5i
EC50_PDE5 = 50.0     // sildenafil EC50 (ng/mL)
Emax_ERA  = 0.70     // max ET-1 block by ERA
EC50_ERA  = 600.0    // bosentan EC50 (ng/mL)
Emax_PGI2a = 0.90    // max IP-R activation by iloprost
EC50_PGI2a = 0.5     // iloprost EC50 (ng/mL)
Emax_Praz = 0.60     // max α1-block by prazosin
EC50_Praz = 2.0      // prazosin EC50 (ng/mL)

// --- Disease subtype switch (0=primary, 1=secondary SSc) ---
secondary = 0
// Secondary modifiers applied below in code

// --- Cold challenge (0=basal, 1=active) ---
cold_challenge = 0

[CMT]
// Drug PK (10 compartments)
GUT_NIF PLASMA_NIF
GUT_SIL PLASMA_SIL
GUT_BOS PLASMA_BOS
PLASMA_ILO
GUT_PRA PLASMA_PRA

// Disease PD (8 compartments)
NE_lvl        // sympathetic NE (nmol/L)
RhoA_GTP      // active RhoA (AU)
Cai_VSMC      // [Ca2+]i in VSMC (μM)
cGMP_VSMC     // cGMP in VSMC (nM)
cAMP_VSMC     // cAMP in VSMC (nM)
ET1_plasma    // plasma ET-1 (pg/mL)
ROS_lvl       // ROS/oxidative stress (AU)
DBF           // digital blood flow (mL/min/100g)

[INIT]
GUT_NIF = 0
PLASMA_NIF = 0
GUT_SIL = 0
PLASMA_SIL = 0
GUT_BOS = 0
PLASMA_BOS = 0
PLASMA_ILO = 0
GUT_PRA = 0
PLASMA_PRA = 0
NE_lvl    = 0.05
RhoA_GTP  = 0.10
Cai_VSMC  = 0.12
cGMP_VSMC = 0.80
cAMP_VSMC = 0.50
ET1_plasma = 0.80
ROS_lvl   = 0.20
DBF       = 5.0

[ODE]
// ── Secondary modifier ──────────────────────────────────────
double sec_mod = 1.0 + secondary * 0.5;  // ET-1 and eNOS mod
double eNOS_eff = eNOS_base * (1.0 - secondary * 0.35);

// ── Cold/stress stimulus ─────────────────────────────────────
double cold_stim = cold_challenge * k_cold_NE;

// ── Drug concentrations ──────────────────────────────────────
double Cp_NIF = PLASMA_NIF / Vd_NIF;   // ng/mL
double Cp_SIL = PLASMA_SIL / Vd_SIL;
double Cp_BOS = PLASMA_BOS / Vd_BOS;
double Cp_ILO = PLASMA_ILO / Vd_ILO;
double Cp_PRA = PLASMA_PRA / Vd_PRA;

// ── Drug PD effects ───────────────────────────────────────────
double E_CCB  = Emax_CCB  * Cp_NIF / (EC50_CCB  + Cp_NIF);
double E_PDE5 = Emax_PDE5 * Cp_SIL / (EC50_PDE5 + Cp_SIL);
double E_ERA  = Emax_ERA  * Cp_BOS / (EC50_ERA  + Cp_BOS);
double E_PGI2 = Emax_PGI2a * Cp_ILO / (EC50_PGI2a + Cp_ILO);
double E_Praz = Emax_Praz * Cp_PRA / (EC50_Praz + Cp_PRA);

// ── NE dynamics ──────────────────────────────────────────────
double NE_synth = k_NE_base * alpha2_sens + cold_stim;
dxdt_NE_lvl = NE_synth - k_NE_clear * NE_lvl;

// ── RhoA activation ──────────────────────────────────────────
double RhoA_input = kRho_NE * NE_lvl
                  + k_ETA_Rho * ET1_plasma / (0.5 + ET1_plasma);
dxdt_RhoA_GTP = RhoA_input - k_ROCK_clear * RhoA_GTP;

// ── VSMC [Ca2+]i ─────────────────────────────────────────────
// VGCC opening driven by NE/ET1, blocked by CCB and cGMP/cAMP
double VGCC_open = k_VGCC_open * NE_lvl * (1.0 - E_CCB)
                 * (1.0 - 0.4 * cGMP_VSMC / (EC50_cGMP + cGMP_VSMC))
                 * (1.0 - 0.3 * cAMP_VSMC / (EC50_cAMP + cAMP_VSMC));
double Ca_from_SR = k_alpha2_constr * NE_lvl;
double Ca_efflux  = k_Ca_clear * Cai_VSMC;
dxdt_Cai_VSMC = VGCC_open + Ca_from_SR - Ca_efflux;

// ── cGMP dynamics ────────────────────────────────────────────
double NO_conc = eNOS_eff * k_eNOS_NO / (1.0 + 0.5 * ROS_lvl);
double cGMP_synth = k_sGC_cGMP * NO_conc;
// PDE5 clearance reduced by sildenafil
double PDE5_activity = 1.0 - E_PDE5;
dxdt_cGMP_VSMC = cGMP_synth - k_cGMP_clear * PDE5_activity * cGMP_VSMC;

// ── cAMP dynamics ────────────────────────────────────────────
double PGI2_eff = PGI2_base * (1.0 - secondary * 0.30) + E_PGI2 * 3.0;
double cAMP_synth = k_IP_cAMP * PGI2_eff;
dxdt_cAMP_VSMC = cAMP_synth - k_cAMP_clear * cAMP_VSMC;

// ── ET-1 dynamics ─────────────────────────────────────────────
double ET1_synth_rate = (k_ET1_synth + k_ET1_inf * ROS_lvl) * sec_mod
                       * (1.0 - E_ERA);
dxdt_ET1_plasma = ET1_synth_rate - k_ET1_clear * ET1_plasma;

// ── ROS dynamics ─────────────────────────────────────────────
double ROS_from_NE   = k_ROS_prod * NE_lvl;
double ROS_ischemia  = 0.05 * (1.0 - DBF / DBF_max);
dxdt_ROS_lvl = ROS_base * 0.1 + ROS_from_NE + ROS_ischemia
             - k_ROS_decay * ROS_lvl;

// ── Digital blood flow ───────────────────────────────────────
// VSMC contraction index
double VSMC_idx = (Cai_VSMC / (EC50_Ca_contr + Cai_VSMC))
                + (RhoA_GTP / (EC50_ROCK + RhoA_GTP))
                + (ET1_plasma * k_ETA_constr / (0.5 + ET1_plasma));
// α1-block by prazosin reduces sympathetic vasoconstriction
double alpha1_block = E_Praz * 0.40;
double DBF_target = DBF_max / (1.0 + 2.5 * VSMC_idx * (1.0 - alpha1_block));
double DBF_relax_effect = cGMP_VSMC / (EC50_cGMP + cGMP_VSMC) * 0.5
                        + cAMP_VSMC / (EC50_cAMP + cAMP_VSMC) * 0.3;
DBF_target = DBF_target * (1.0 + DBF_relax_effect);
if(DBF_target > DBF_max) DBF_target = DBF_max;
dxdt_DBF = 5.0 * (DBF_target - DBF);

// ── PK ODEs: Nifedipine ─────────────────────────────────────
dxdt_GUT_NIF    = -ka_NIF * GUT_NIF;
dxdt_PLASMA_NIF =  ka_NIF * F_NIF * GUT_NIF - (CL_NIF / Vd_NIF) * PLASMA_NIF;

// ── PK ODEs: Sildenafil ─────────────────────────────────────
dxdt_GUT_SIL    = -ka_SIL * GUT_SIL;
dxdt_PLASMA_SIL =  ka_SIL * F_SIL * GUT_SIL - (CL_SIL / Vd_SIL) * PLASMA_SIL;

// ── PK ODEs: Bosentan ────────────────────────────────────────
dxdt_GUT_BOS    = -ka_BOS * GUT_BOS;
dxdt_PLASMA_BOS =  ka_BOS * F_BOS * GUT_BOS - (CL_BOS / Vd_BOS) * PLASMA_BOS;

// ── PK ODEs: Iloprost (IV) ───────────────────────────────────
dxdt_PLASMA_ILO = -(CL_ILO / Vd_ILO) * PLASMA_ILO;

// ── PK ODEs: Prazosin ────────────────────────────────────────
dxdt_GUT_PRA    = -ka_PRA * GUT_PRA;
dxdt_PLASMA_PRA =  ka_PRA * F_PRA * GUT_PRA - (CL_PRA / Vd_PRA) * PLASMA_PRA;

[TABLE]
// ── Derived clinical endpoints ───────────────────────────────
double VSMC_contraction = (Cai_VSMC / (EC50_Ca_contr + Cai_VSMC))
                        + (RhoA_GTP / (EC50_ROCK + RhoA_GTP));

// Vasospasm episodes per week (driven by DBF reduction)
double VasoEp_rate = EpiFreq_base * exp(-2.0 * DBF / DBF_base)
                   * alpha2_sens * (1.0 + secondary * 0.8);
double VasoEp_wk   = VasoEp_rate;  // model-predicted weekly episodes

// Attack duration (min): inversely related to DBF recovery
double AttackDur_min = 15.0 * exp(-DBF / DBF_base) + 5.0;

// VAS pain (0-10)
double VAS_pain = 8.0 * (1.0 - DBF / DBF_max);

// Raynaud Condition Score (0-10, composite)
double RCS = 0.4 * VAS_pain + 0.3 * VasoEp_wk + 0.3 * (10.0 - DBF);
if(RCS > 10) RCS = 10.0;
if(RCS < 0)  RCS = 0.0;

// Digital ulcer risk (% probability per episode in secondary)
double DU_risk = secondary * 0.15 * (1.0 - DBF / DBF_max);

// Plasma drug concentrations (ng/mL)
double Cp_NIF_out = PLASMA_NIF / Vd_NIF;
double Cp_SIL_out = PLASMA_SIL / Vd_SIL;
double Cp_BOS_out = PLASMA_BOS / Vd_BOS;
double Cp_ILO_out = PLASMA_ILO / Vd_ILO;
double Cp_PRA_out = PLASMA_PRA / Vd_PRA;

[CAPTURE]
Cp_NIF_out Cp_SIL_out Cp_BOS_out Cp_ILO_out Cp_PRA_out
VSMC_contraction DBF VasoEp_wk AttackDur_min
VAS_pain RCS DU_risk
NE_lvl RhoA_GTP Cai_VSMC cGMP_VSMC cAMP_VSMC
ET1_plasma ROS_lvl

// ================================================================
// R SIMULATION SCENARIOS (run from this file after mread)
// ================================================================

/*** R
library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

mod <- mread("raynaud_mrgsolve_model.R")

# ── Helper: weekly dosing regimen ──────────────────────────────
make_ev <- function(cmt, amt, ii, addl, start=0) {
  ev(cmt=cmt, amt=amt, ii=ii, addl=addl, time=start)
}

# Simulation period: 84 days = 2016 hours
sim_end <- 2016

# ──────────────────────────────────────────────────────────────
# SCENARIO 1: Untreated Primary Raynaud's
# ──────────────────────────────────────────────────────────────
e1 <- ev(time=0, cmt="GUT_NIF", amt=0)  # no drug
out1 <- mod %>%
  param(secondary=0, cold_challenge=0) %>%
  mrgsim(ev=e1, end=sim_end, delta=1) %>%
  as.data.frame() %>%
  mutate(scenario="Untreated Primary")

# ──────────────────────────────────────────────────────────────
# SCENARIO 2: Nifedipine 30mg QD (CCB, standard first-line)
# ──────────────────────────────────────────────────────────────
e2 <- make_ev("GUT_NIF", amt=30*1000, ii=24, addl=83)  # 84 days
out2 <- mod %>%
  param(secondary=0, cold_challenge=0) %>%
  mrgsim(ev=e2, end=sim_end, delta=1) %>%
  as.data.frame() %>%
  mutate(scenario="Nifedipine 30mg QD")

# ──────────────────────────────────────────────────────────────
# SCENARIO 3: Sildenafil 50mg BID (PDE5i)
# ──────────────────────────────────────────────────────────────
e3a <- make_ev("GUT_SIL", amt=50*1000, ii=12, addl=167)  # BID ×84d
out3 <- mod %>%
  param(secondary=0, cold_challenge=0) %>%
  mrgsim(ev=e3a, end=sim_end, delta=1) %>%
  as.data.frame() %>%
  mutate(scenario="Sildenafil 50mg BID")

# ──────────────────────────────────────────────────────────────
# SCENARIO 4: Bosentan 125mg BID (ERA — secondary/SSc)
# ──────────────────────────────────────────────────────────────
e4 <- make_ev("GUT_BOS", amt=125*1000, ii=12, addl=167)
out4 <- mod %>%
  param(secondary=1, cold_challenge=0) %>%
  mrgsim(ev=e4, end=sim_end, delta=1) %>%
  as.data.frame() %>%
  mutate(scenario="Bosentan 125mg BID (2°)")

# ──────────────────────────────────────────────────────────────
# SCENARIO 5: Iloprost IV infusion (0.5-2 ng/kg/min × 6h/d × 5d)
#   Modeled as bolus dose cycles
# ──────────────────────────────────────────────────────────────
ilo_doses <- do.call(rbind, lapply(0:4, function(d) {
  ev(time=d*24, cmt="PLASMA_ILO", amt=50, ii=0, addl=0)
}))
out5 <- mod %>%
  param(secondary=1, cold_challenge=0) %>%
  mrgsim(ev=ilo_doses, end=sim_end, delta=1) %>%
  as.data.frame() %>%
  mutate(scenario="Iloprost IV ×5d (2°)")

# ──────────────────────────────────────────────────────────────
# SCENARIO 6: Prazosin 1mg BID (α1-blocker)
# ──────────────────────────────────────────────────────────────
e6 <- make_ev("GUT_PRA", amt=1*1000, ii=12, addl=167)
out6 <- mod %>%
  param(secondary=0, cold_challenge=0) %>%
  mrgsim(ev=e6, end=sim_end, delta=1) %>%
  as.data.frame() %>%
  mutate(scenario="Prazosin 1mg BID")

# ──────────────────────────────────────────────────────────────
# SCENARIO 7: Nifedipine + Sildenafil combination
# ──────────────────────────────────────────────────────────────
e7_NIF <- make_ev("GUT_NIF", amt=30*1000, ii=24, addl=83)
e7_SIL <- make_ev("GUT_SIL", amt=25*1000, ii=12, addl=167)
e7 <- c(e7_NIF, e7_SIL)
out7 <- mod %>%
  param(secondary=0, cold_challenge=0) %>%
  mrgsim(ev=e7, end=sim_end, delta=1) %>%
  as.data.frame() %>%
  mutate(scenario="Nifedipine + Sildenafil")

# ──────────────────────────────────────────────────────────────
# SCENARIO 8: Untreated Secondary Raynaud's (SSc)
# ──────────────────────────────────────────────────────────────
out8 <- mod %>%
  param(secondary=1, cold_challenge=0) %>%
  mrgsim(ev=e1, end=sim_end, delta=1) %>%
  as.data.frame() %>%
  mutate(scenario="Untreated Secondary (SSc)")

# ──────────────────────────────────────────────────────────────
# SCENARIO 9: Cold challenge simulation (acute, 1h exposure)
# ──────────────────────────────────────────────────────────────
cold_ev <- ev(time=0, cmt="GUT_NIF", amt=0)
out9 <- mod %>%
  param(secondary=0, cold_challenge=1) %>%
  mrgsim(ev=cold_ev, end=24, delta=0.1) %>%
  as.data.frame() %>%
  mutate(scenario="Cold Challenge")

# ──────────────────────────────────────────────────────────────
# Combine & plot key endpoints
all_out <- bind_rows(out1, out2, out3, out4, out5, out6, out7, out8)
all_out$time_wk <- all_out$time / 168  # hours → weeks

# Weekly vasospasm episodes plot
ggplot(all_out %>% filter(time_wk <= 12), aes(x=time_wk, y=VasoEp_wk, color=scenario)) +
  geom_line(size=1) +
  labs(title="Predicted Vasospasm Episodes/Week",
       x="Time (weeks)", y="Episodes/week", color="Scenario") +
  theme_bw()

# RCS over time
ggplot(all_out %>% filter(time_wk <= 12), aes(x=time_wk, y=RCS, color=scenario)) +
  geom_line(size=1) +
  labs(title="Raynaud Condition Score (0-10)",
       x="Time (weeks)", y="RCS", color="Scenario") +
  theme_bw()

# Digital blood flow
ggplot(all_out %>% filter(time_wk <= 12), aes(x=time_wk, y=DBF, color=scenario)) +
  geom_line(size=1) +
  labs(title="Digital Blood Flow (mL/min/100g)",
       x="Time (weeks)", y="DBF", color="Scenario") +
  theme_bw()

# ET-1 in secondary scenarios
et1_scen <- bind_rows(out4, out5, out8)
ggplot(et1_scen %>% filter(time_wk <= 12), aes(x=time_wk, y=ET1_plasma, color=scenario)) +
  geom_line(size=1) +
  labs(title="Plasma ET-1 (pg/mL) — Secondary Raynaud",
       x="Time (weeks)", y="ET-1 (pg/mL)") +
  theme_bw()

cat("Raynaud's QSP simulation complete. Key endpoint summary:\n")
summary_df <- all_out %>%
  group_by(scenario) %>%
  filter(time_wk >= 11) %>%
  summarise(
    Mean_VasoEp = mean(VasoEp_wk, na.rm=TRUE),
    Mean_DBF    = mean(DBF, na.rm=TRUE),
    Mean_RCS    = mean(RCS, na.rm=TRUE),
    Mean_ET1    = mean(ET1_plasma, na.rm=TRUE)
  )
print(summary_df)
***/
