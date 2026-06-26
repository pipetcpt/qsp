## ============================================================
## CTEPH QSP mrgsolve Model
## Chronic Thromboembolic Pulmonary Hypertension
## ============================================================
## Compartments (19 ODEs):
##   Drug PK (7): Riociguat C/P/Met, Macitentan C/P/Met, Treprostinil C
##   Disease PD (12): TB, PVR_fixed, PVR_var, ET1, cGMP, cAMP,
##                    RV_work, mPAP, CO, SaO2, BNP, 6MWD
## Treatment scenarios (6):
##   1. No treatment (natural history)
##   2. Riociguat monotherapy (2.5 mg TID)
##   3. Macitentan monotherapy (10 mg QD)
##   4. Riociguat + Macitentan combination
##   5. BPA procedure + Riociguat
##   6. Post-PEA surgery
## References: CHEST-1 (2013), SERAPHIN (2013), MERIT (2022),
##             CTREPH (2018), Galie et al. ESC/ERS Guidelines 2022
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ---- Model code block ----------------------------------------
cteph_code <- '
$PROB
CTEPH QSP Model
Riociguat + Macitentan + Treprostinil + PEA/BPA
Version 1.0 | 2026-06-26

$PARAM
// ---- Riociguat PK (2-compartment, TID dosing) ----
// Ref: Bayer data, Frey et al. 2014 J Clin Pharmacol
CL_RIO   = 2.4    // L/h clearance
V1_RIO   = 30.0   // L central volume
V2_RIO   = 20.0   // L peripheral volume
Q_RIO    = 1.8    // L/h inter-compartmental
ka_RIO   = 0.9    // 1/h absorption rate
F_RIO    = 0.94   // bioavailability
fm_RIO   = 0.15   // fraction metabolized to M1

// ---- Riociguat M1 Metabolite PK ----
CL_M1    = 3.6    // L/h (active metabolite)
V_M1     = 25.0   // L

// ---- Macitentan PK (1-compartment + metabolite) ----
// Ref: Sidharta et al. 2012 Br J Clin Pharmacol
CL_MAC   = 1.1    // L/h
V1_MAC   = 50.0   // L central
V2_MAC   = 40.0   // L peripheral
Q_MAC    = 0.8    // L/h
ka_MAC   = 0.35   // 1/h (Tmax ~8h)
F_MAC    = 0.75
fm_MAC   = 0.70   // ~70% converted to ACT-132577

// ---- ACT-132577 (active metabolite of macitentan) ----
CL_ACT   = 0.55   // L/h
V_ACT    = 60.0   // L

// ---- Treprostinil PK (SC infusion, 1-compartment) ----
// Ref: McNeil et al. 2013, Wade et al. 2004
CL_TREP  = 4.0    // L/h
V_TREP   = 14.0   // L
ka_TREP  = 2.0    // 1/h (SC absorption)
F_TREP   = 0.79   // SC bioavailability

// ---- Riociguat PD parameters ----
// sGC stimulation -> cGMP
Emax_RIO    = 0.70   // maximal cGMP increase fraction
EC50_RIO    = 2.5    // ng/mL (Cp that gives 50% Emax)
kin_cGMP    = 0.15   // pmol/mL/h cGMP production baseline
kout_cGMP   = 0.10   // 1/h cGMP degradation (PDE5)
cGMP0       = 1.5    // pmol/mL baseline

// ---- Macitentan/ERA PD parameters ----
// ET-1 receptor blockade -> PVR variable reduction
Emax_ERA    = 0.55   // maximal PVR_var reduction
EC50_ERA    = 0.8    // ng/mL (unbound macitentan)
ku_ERA      = 0.01   // dissociation constant correction

// ---- Treprostinil PD parameters ----
// IP receptor -> cAMP -> vasodilation
Emax_TREP   = 0.45   // maximal PVR_var reduction
EC50_TREP   = 3.0    // ng/mL
kin_cAMP    = 0.12   // pmol/mL/h baseline cAMP production
kout_cAMP   = 0.08   // 1/h cAMP degradation
cAMP0       = 1.5    // pmol/mL baseline

// ---- Disease Parameters ----
// Thrombotic burden (TB): normalized 0=none, 1=max burden
TB0         = 0.75   // baseline (75% typical at diagnosis)
kdeg_TB     = 0.002  // 1/h very slow natural thrombus resolution
kform_TB    = 0.003  // 1/h thrombus formation/organization
TB_PE_effect= 0.15   // increment in TB per recurrent PE

// PVR_fixed (dyn.s/cm5): fixed mechanical obstruction
PVRf0       = 550    // dyn.s/cm5 baseline fixed component
// Note: PEA reduces PVRf by ~70%; BPA by ~10-15% per session

// PVR_variable (dyn.s/cm5): vasomotor/remodeling component
PVRv0       = 350    // dyn.s/cm5 baseline variable component
kin_PVRv    = 0.010  // 1/h PVR_var progression rate
kout_PVRv   = 0.008  // 1/h natural decay
PVRv_ET1fx  = 0.25   // ET-1 contribution fraction
PVRv_HVfx   = 0.15   // hypoxic vasoconstriction fraction

// ET-1 kinetics
ET1_0       = 3.5    // pg/mL baseline plasma ET-1
kin_ET1     = 0.50   // pg/mL/h production
kout_ET1    = 0.14   // 1/h clearance

// Cardiac parameters
mPAP0       = 46     // mmHg baseline mean PAP (CHEST-1 criteria >25)
CO0         = 3.6    // L/min baseline cardiac output (reduced)
RVwork0     = 18     // g.m/m2 baseline RV stroke work index (elevated)
SaO2_0      = 90     // % baseline arterial O2 saturation
BNP0        = 180    // pg/mL baseline BNP

// 6MWD
sixMWD0     = 342    // m baseline (CHEST-1 mean)
k6MWD_CO    = 45     // m per L/min CO improvement
k6MWD_SaO2  = 2.5    // m per % SaO2 improvement

// Procedural/Surgical intervention flags
PEA_done    = 0      // 1 = PEA performed (reduces PVRf by 70%)
BPA_done    = 0      // 1 = BPA performed (reduces PVRf by ~35%)
PEA_time    = 0      // hours post-PEA
BPA_sessions= 0      // number of BPA sessions (each ~10% PVRf reduction)

// Anticoagulation effect
AC_effect   = 1      // 1=on anticoagulation (reduces TB formation)

$CMT
// Drug PK compartments
DEPOT_RIO     // Riociguat absorption depot
C1_RIO        // Riociguat central
C2_RIO        // Riociguat peripheral
MET_RIO       // M1 metabolite central

DEPOT_MAC     // Macitentan depot
C1_MAC        // Macitentan central
C2_MAC        // Macitentan peripheral
MET_MAC       // ACT-132577 metabolite

DEPOT_TREP    // Treprostinil depot
C1_TREP       // Treprostinil central

// Disease state variables
TB            // Thrombotic burden (normalized 0-1)
PVR_fixed     // Fixed PVR (dyn.s/cm5) — thrombus-related
PVR_var       // Variable PVR (dyn.s/cm5) — vasomotor/remodeling
ET1           // Plasma endothelin-1 (pg/mL)
cGMP          // cGMP levels (pmol/mL)
cAMP          // cAMP levels (pmol/mL)
RV_work       // RV stroke work index (g.m/m2)
mPAP          // Mean PAP (mmHg)
CO            // Cardiac output (L/min)
SaO2          // Arterial O2 saturation (%)
BNP           // BNP (pg/mL)
sixMWD        // 6-minute walk distance (m)

$INIT
DEPOT_RIO = 0, C1_RIO = 0, C2_RIO = 0, MET_RIO = 0,
DEPOT_MAC = 0, C1_MAC = 0, C2_MAC = 0, MET_MAC = 0,
DEPOT_TREP = 0, C1_TREP = 0,
TB      = 0.75,
PVR_fixed = 550,
PVR_var   = 350,
ET1     = 3.5,
cGMP    = 1.5,
cAMP    = 1.5,
RV_work = 18,
mPAP    = 46,
CO      = 3.6,
SaO2    = 90,
BNP     = 180,
sixMWD  = 342

$ODE
// ==============================================================
// RIOCIGUAT PK
// ==============================================================
double k10_RIO = CL_RIO / V1_RIO;
double k12_RIO = Q_RIO / V1_RIO;
double k21_RIO = Q_RIO / V2_RIO;
double kmet_RIO = fm_RIO * k10_RIO;

dxdt_DEPOT_RIO = -ka_RIO * DEPOT_RIO;
dxdt_C1_RIO    =  ka_RIO * F_RIO * DEPOT_RIO
                 - (k10_RIO + k12_RIO) * C1_RIO
                 + k21_RIO * C2_RIO;
dxdt_C2_RIO    =  k12_RIO * C1_RIO - k21_RIO * C2_RIO;
dxdt_MET_RIO   =  kmet_RIO * C1_RIO * V1_RIO / V_M1
                 - (CL_M1/V_M1) * MET_RIO;

// Riociguat plasma concentrations (ng/mL)
double Cp_RIO  = C1_RIO / V1_RIO;   // ng/mL
double Cm_RIO  = MET_RIO / V_M1;    // ng/mL M1

// ==============================================================
// MACITENTAN PK
// ==============================================================
double k10_MAC = CL_MAC / V1_MAC;
double k12_MAC = Q_MAC  / V1_MAC;
double k21_MAC = Q_MAC  / V2_MAC;
double kmet_MAC = fm_MAC * k10_MAC;

dxdt_DEPOT_MAC = -ka_MAC * DEPOT_MAC;
dxdt_C1_MAC    =  ka_MAC * F_MAC * DEPOT_MAC
                 - (k10_MAC + k12_MAC) * C1_MAC
                 + k21_MAC * C2_MAC;
dxdt_C2_MAC    =  k12_MAC * C1_MAC - k21_MAC * C2_MAC;
dxdt_MET_MAC   =  kmet_MAC * C1_MAC * V1_MAC / V_ACT
                 - (CL_ACT/V_ACT) * MET_MAC;

double Cp_MAC  = C1_MAC / V1_MAC;   // ng/mL macitentan
double Cm_MAC  = MET_MAC / V_ACT;   // ng/mL ACT-132577

// ==============================================================
// TREPROSTINIL PK (SC infusion modeled as bolus/depot)
// ==============================================================
dxdt_DEPOT_TREP = -ka_TREP * DEPOT_TREP;
dxdt_C1_TREP    =  ka_TREP * F_TREP * DEPOT_TREP
                  - (CL_TREP/V_TREP) * C1_TREP;

double Cp_TREP = C1_TREP / V_TREP; // ng/mL

// ==============================================================
// PHARMACODYNAMICS
// ==============================================================

// ---- cGMP: Riociguat sGC stimulation ----
double E_RIO = (Emax_RIO * Cp_RIO) / (EC50_RIO + Cp_RIO);
double E_M1  = (Emax_RIO * 0.4 * Cm_RIO) / (EC50_RIO * 0.8 + Cm_RIO); // M1 ~40% activity
double E_total_RIO = 1.0 + E_RIO + E_M1;  // stimulation fold
dxdt_cGMP = kin_cGMP * E_total_RIO - kout_cGMP * cGMP;

// ---- cAMP: Treprostinil IP receptor ----
double E_TREP = (Emax_TREP * Cp_TREP) / (EC50_TREP + Cp_TREP);
double E_total_TREP = 1.0 + E_TREP;
dxdt_cAMP = kin_cAMP * E_total_TREP - kout_cAMP * cAMP;

// ---- ET-1 kinetics ----
// ERA blocks ET receptor => plasma ET-1 rises (pharmacokinetic paradox)
double ERA_total = (Cp_MAC + 0.5 * Cm_MAC);  // total ERA activity
double ET1_feedback = 1.0 + 0.3 * ERA_total / (0.5 + ERA_total); // ET-1 plasma rise
dxdt_ET1 = kin_ET1 * ET1_feedback - kout_ET1 * ET1;

// ---- Thrombotic Burden ----
double AC_factor = AC_effect * 0.6;   // anticoag reduces formation by 60%
double kform_eff = kform_TB * (1.0 - AC_factor);
double kdeg_eff  = kdeg_TB * (1.0 + 0.5 * E_TREP); // prostacyclin aids fibrinolysis slightly
double TB_capped = (TB < 0) ? 0 : (TB > 1) ? 1 : TB;
dxdt_TB = kform_eff * (1.0 - TB_capped) - kdeg_eff * TB_capped;

// ---- PVR_fixed: thrombus-related obstruction ----
// PEA/BPA interventions applied as step reduction
double PEA_reduction = PEA_done * 0.70;   // PEA removes ~70% of fixed PVR
double BPA_reduction = BPA_done * (1.0 - pow(0.88, BPA_sessions)); // each session ~12% reduction
double PVRf_floor    = 80.0; // minimum residual PVRf (dyn.s/cm5)
double PVRf_target   = PVR_fixed * (1.0 - PEA_reduction) * (1.0 - BPA_reduction);
if (PVRf_target < PVRf_floor) PVRf_target = PVRf_floor;
// slow re-organization pressure drives toward target
double kPVRf = 0.0015; // 1/h rate of structural change
dxdt_PVR_fixed = -kPVRf * (PVR_fixed - PVRf_target);

// ---- PVR_variable: vasomotor/remodeling component ----
// Drug effects on variable PVR
double cGMP_ratio = cGMP / cGMP0;
double cAMP_ratio = cAMP / cAMP0;
double ERA_eff    = (Emax_ERA * ERA_total) / (EC50_ERA + ERA_total);

// Contributions to PVR_var change
double ET1_ratio  = ET1 / ET1_0;                // ET-1 drives PVR up
double PVRv_drive = kin_PVRv * ET1_ratio;       // ET-1 driven increase
double PVRv_inhibit = kout_PVRv * (cGMP_ratio + cAMP_ratio + ERA_eff); // drugs reduce
double PVRv_floor = 100.0;
dxdt_PVR_var = PVRv_drive - PVRv_inhibit * PVR_var;
// clamp to floor in capture
if (PVR_var < PVRv_floor) dxdt_PVR_var = (dxdt_PVR_var < 0) ? 0 : dxdt_PVR_var;

// ==============================================================
// HEMODYNAMICS - Algebraic-style approximation via slow ODEs
// ==============================================================

// Total PVR
double PVR_total = PVR_fixed + PVR_var;  // dyn.s/cm5

// mPAP target: Wood formula mPAP = CO * PVR/80 + PCWP_est
// PCWP ~10 mmHg in CTEPH (pre-capillary)
double PCWP_est  = 10.0;
double mPAP_target = (CO * PVR_total / 80.0) + PCWP_est;
dxdt_mPAP = 0.02 * (mPAP_target - mPAP); // 1/h adaptation

// Cardiac output: depends on RV function
// As PVR rises, CO falls (afterload-sensitive RV)
double CO_target;
double PVR_norm  = PVR_total / 900.0; // normalized to severe value
if (PVR_norm > 1.0) PVR_norm = 1.0;
CO_target = 5.2 * (1.0 - 0.45 * PVR_norm); // max 5.2 L/min, falls with PVR
// RV function effect
double RV_adj    = (RV_work < 15.0) ? (0.7 + 0.02 * RV_work) : 1.0; // low RVwork -> reduced CO
CO_target = CO_target * RV_adj;
if (CO_target < 1.5) CO_target = 1.5;
dxdt_CO = 0.015 * (CO_target - CO);

// RV stroke work index (proportional to mPAP and CO)
// RVSWI = (mPAP - mRAP) * SV * 0.0136; SV = CO/HR * 1000
// Simplified: rises with pressure load, falls with treatment
double HR_est    = 80.0;  // bpm assumed
double SV_est    = (CO * 1000.0) / HR_est;  // mL
double mRAP_est  = 8.0 + 0.15 * (mPAP - 20.0);  // estimated mRAP rises with mPAP
double RVwork_target = (mPAP - mRAP_est) * SV_est * 0.0136;
dxdt_RV_work = 0.01 * (RVwork_target - RV_work);

// SaO2: inversely related to V/Q mismatch (PVR_fixed drives V/Q mismatch)
// Simplified: SaO2 falls with fixed PVR (dead space effect)
double SaO2_target = 99.0 - 0.012 * PVR_fixed - 0.008 * (CO0 - CO) * 10.0;
if (SaO2_target < 70.0) SaO2_target = 70.0;
if (SaO2_target > 99.0) SaO2_target = 99.0;
dxdt_SaO2 = 0.05 * (SaO2_target - SaO2);

// BNP: rises with RV wall stress (mPAP * RV_volume proxy)
// Simplified: BNP proportional to mPAP and inversely to CO
double BNP_target = 20.0 * exp(0.08 * (mPAP - 20.0)) * (5.0 / (CO + 0.01));
if (BNP_target > 2000.0) BNP_target = 2000.0;
dxdt_BNP = 0.02 * (BNP_target - BNP);

// 6MWD: clinical composite endpoint
double sixMWD_target = sixMWD0
  + k6MWD_CO * (CO - CO0)
  + k6MWD_SaO2 * (SaO2 - SaO2_0)
  - 0.05 * (mPAP - mPAP0);
if (sixMWD_target < 50)  sixMWD_target = 50;
if (sixMWD_target > 600) sixMWD_target = 600;
dxdt_sixMWD = 0.008 * (sixMWD_target - sixMWD);

$TABLE
// Derived PK outputs
double Cp_RIO_out  = C1_RIO / V1_RIO;
double Cp_MAC_out  = C1_MAC / V1_MAC;
double Cp_TREP_out = C1_TREP / V_TREP;
double Cm_RIO_out  = MET_RIO / V_M1;
double Cm_MAC_out  = MET_MAC / V_ACT;

// Total PVR
double PVR_total_out = PVR_fixed + PVR_var;

// WHO FC approximation (1-4 scale)
double WHO_FC;
if (sixMWD > 440)       WHO_FC = 1.5;
else if (sixMWD > 300)  WHO_FC = 2.0;
else if (sixMWD > 150)  WHO_FC = 3.0;
else                    WHO_FC = 4.0;

// Haemodynamic response flag
double HAEMO_RESP = (PVR_total_out < 480 && mPAP < 38) ? 1 : 0;

$CAPTURE
Cp_RIO_out Cm_RIO_out Cp_MAC_out Cm_MAC_out Cp_TREP_out
ET1 cGMP cAMP
TB PVR_fixed PVR_var PVR_total_out
mPAP CO SaO2 RV_work BNP sixMWD
WHO_FC HAEMO_RESP
'

## ---- Compile model -------------------------------------------
mod <- mcode("cteph", cteph_code)

## ============================================================
## TREATMENT SCENARIOS
## ============================================================

# Simulation parameters
TEND   <- 52 * 7 * 24  # 52 weeks in hours
DELTA  <- 12            # output every 12 hours

# Create dosing events
# Riociguat: 2.5 mg TID (every 8h)
ev_rio <- ev(amt = 2.5, cmt = "DEPOT_RIO", ii = 8, addl = round(TEND/8))

# Macitentan: 10 mg QD (every 24h)
ev_mac <- ev(amt = 10,  cmt = "DEPOT_MAC", ii = 24, addl = round(TEND/24))

# Treprostinil SC: modeled as bolus doses q6h (approximation of infusion)
# Dose escalated: start 1.25 ng/kg/min, target ~40-60 ng/kg/min
# For 70 kg: ~3 µg/h start → ~170 µg/h target; here using mg units simplified
ev_trep <- ev(amt = 0.5, cmt = "DEPOT_TREP", ii = 6, addl = round(TEND/6))

# BPA: modeled as step change in parameters at week 4
# PEA: modeled as step change at day 0 (post-PEA simulation)

## ---- Helper: run scenario ------------------------------------
run_scenario <- function(model, events, params_override = list(),
                         tend = TEND, delta = DELTA) {
  mod2 <- param(model, params_override)
  if (is.null(events)) {
    mrgsim(mod2, end = tend, delta = delta) %>% as.data.frame()
  } else {
    mrgsim(mod2, events = events, end = tend, delta = delta) %>% as.data.frame()
  }
}

## ---- Scenario 1: Natural history (no treatment) --------------
message("Scenario 1: Natural history...")
scen1 <- run_scenario(mod, NULL,
  params_override = list(AC_effect = 0, PEA_done = 0, BPA_done = 0))
scen1$scenario <- "1. No Treatment"

## ---- Scenario 2: Anticoagulation only ------------------------
message("Scenario 2: Anticoagulation alone...")
scen2 <- run_scenario(mod, NULL,
  params_override = list(AC_effect = 1, PEA_done = 0, BPA_done = 0))
scen2$scenario <- "2. Anticoagulation Only"

## ---- Scenario 3: Riociguat monotherapy (2.5 mg TID) ---------
message("Scenario 3: Riociguat monotherapy...")
scen3 <- run_scenario(mod, ev_rio,
  params_override = list(AC_effect = 1, PEA_done = 0, BPA_done = 0))
scen3$scenario <- "3. Riociguat 2.5 mg TID"

## ---- Scenario 4: Macitentan monotherapy (10 mg QD) ----------
message("Scenario 4: Macitentan monotherapy...")
scen4 <- run_scenario(mod, ev_mac,
  params_override = list(AC_effect = 1, PEA_done = 0, BPA_done = 0))
scen4$scenario <- "4. Macitentan 10 mg QD"

## ---- Scenario 5: Riociguat + Macitentan combination ---------
message("Scenario 5: Riociguat + Macitentan combination...")
ev_combo <- ev_rio + ev_mac
scen5 <- run_scenario(mod, ev_combo,
  params_override = list(AC_effect = 1, PEA_done = 0, BPA_done = 0))
scen5$scenario <- "5. Riociguat + Macitentan"

## ---- Scenario 6: BPA (5 sessions) + Riociguat ---------------
message("Scenario 6: BPA + Riociguat...")
scen6 <- run_scenario(mod, ev_rio,
  params_override = list(AC_effect = 1, PEA_done = 0, BPA_done = 1,
                         BPA_sessions = 5))
scen6$scenario <- "6. BPA (5x) + Riociguat"

## ---- Scenario 7 (bonus): Post-PEA + medical therapy ---------
message("Scenario 7: Post-PEA + Riociguat + Macitentan...")
scen7 <- run_scenario(mod, ev_combo,
  params_override = list(AC_effect = 1, PEA_done = 1, BPA_done = 0,
                         PVR_fixed = 165,  # post-PEA residual
                         mPAP = 30, CO = 4.5, SaO2 = 95, BNP = 60,
                         sixMWD = 440))
scen7$scenario <- "7. Post-PEA + Riociguat+MAC"

## ---- Combine results -----------------------------------------
results <- bind_rows(scen1, scen2, scen3, scen4, scen5, scen6, scen7)
results$time_weeks <- results$time / (7 * 24)

## ============================================================
## VISUALIZATION
## ============================================================

theme_cteph <- theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "bottom",
    legend.title  = element_text(face = "bold"),
    strip.background = element_rect(fill = "#2C3E50"),
    strip.text    = element_text(color = "white", face = "bold")
  )

colors_scen <- c(
  "1. No Treatment"           = "#95A5A6",
  "2. Anticoagulation Only"   = "#E67E22",
  "3. Riociguat 2.5 mg TID"  = "#2980B9",
  "4. Macitentan 10 mg QD"   = "#8E44AD",
  "5. Riociguat + Macitentan" = "#1ABC9C",
  "6. BPA (5x) + Riociguat"  = "#E74C3C",
  "7. Post-PEA + Riociguat+MAC" = "#27AE60"
)

## Plot 1: mPAP over time
p1 <- ggplot(results, aes(x = time_weeks, y = mPAP, color = scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 25, linetype = "dashed", color = "red", alpha = 0.7) +
  scale_color_manual(values = colors_scen) +
  labs(title = "Mean Pulmonary Artery Pressure (mPAP) Over Time",
       x = "Time (weeks)", y = "mPAP (mmHg)", color = "Scenario") +
  annotate("text", x = 5, y = 25.5, label = "PH threshold (25 mmHg)",
           hjust = 0, size = 3, color = "red") +
  theme_cteph

## Plot 2: PVR over time
p2 <- results %>%
  mutate(PVR_total = PVR_fixed + PVR_var) %>%
  ggplot(aes(x = time_weeks, y = PVR_total, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colors_scen) +
  labs(title = "Total Pulmonary Vascular Resistance",
       x = "Time (weeks)", y = "PVR (dyn·s/cm⁵)", color = "Scenario") +
  theme_cteph

## Plot 3: Cardiac output
p3 <- ggplot(results, aes(x = time_weeks, y = CO, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colors_scen) +
  labs(title = "Cardiac Output",
       x = "Time (weeks)", y = "CO (L/min)", color = "Scenario") +
  theme_cteph

## Plot 4: 6MWD
p4 <- ggplot(results, aes(x = time_weeks, y = sixMWD, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colors_scen) +
  labs(title = "6-Minute Walk Distance (6MWD)",
       x = "Time (weeks)", y = "6MWD (m)", color = "Scenario") +
  theme_cteph

## Plot 5: BNP
p5 <- ggplot(results, aes(x = time_weeks, y = BNP, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colors_scen) +
  labs(title = "BNP (Right Ventricular Stress Biomarker)",
       x = "Time (weeks)", y = "BNP (pg/mL)", color = "Scenario") +
  theme_cteph

## Plot 6: PK profiles (Scenario 3: Riociguat, first 48h)
pk_data <- scen3 %>%
  filter(time_weeks <= 2) %>%
  pivot_longer(cols = c(Cp_RIO_out, Cm_RIO_out),
               names_to = "compound", values_to = "concentration")
pk_data$compound <- factor(pk_data$compound,
  levels = c("Cp_RIO_out", "Cm_RIO_out"),
  labels = c("Riociguat Parent", "M1 Metabolite"))

p6 <- ggplot(pk_data, aes(x = time_weeks * 7 * 24, y = concentration,
                           color = compound)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("Riociguat Parent" = "#2980B9",
                                "M1 Metabolite" = "#85C1E9")) +
  labs(title = "Riociguat PK Profile (First 2 Weeks)",
       x = "Time (hours)", y = "Concentration (ng/mL)", color = "") +
  theme_cteph

## Plot 7: cGMP vs cAMP (PD signal)
pd_data <- scen5 %>%
  filter(time_weeks <= 52) %>%
  pivot_longer(cols = c(cGMP, cAMP),
               names_to = "second_messenger", values_to = "level")

p7 <- ggplot(pd_data, aes(x = time_weeks, y = level,
                           color = second_messenger)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = c("cGMP" = "#2980B9", "cAMP" = "#148F77")) +
  labs(title = "Second Messenger Dynamics (Scenario 5: Combination)",
       x = "Time (weeks)", y = "Level (pmol/mL)", color = "Molecule") +
  theme_cteph

## Plot 8: PVR decomposition (fixed vs variable, Scenario 5)
pvr_data <- scen5 %>%
  select(time_weeks, PVR_fixed, PVR_var) %>%
  pivot_longer(cols = c(PVR_fixed, PVR_var),
               names_to = "component", values_to = "value")

p8 <- ggplot(pvr_data, aes(x = time_weeks, y = value, fill = component)) +
  geom_area(alpha = 0.7, position = "stack") +
  scale_fill_manual(values = c("PVR_fixed" = "#C0392B", "PVR_var" = "#E74C3C"),
                    labels = c("Fixed (Thrombus)", "Variable (Vasomotor)")) +
  labs(title = "PVR Decomposition: Fixed vs Variable (Combination Therapy)",
       x = "Time (weeks)", y = "PVR (dyn·s/cm⁵)", fill = "Component") +
  theme_cteph

## Print plots
print(p1); print(p2); print(p3); print(p4)
print(p5); print(p6); print(p7); print(p8)

## ============================================================
## SUMMARY TABLE: Week 52 endpoints
## ============================================================
summary_table <- results %>%
  filter(abs(time_weeks - 52) == min(abs(time_weeks - 52))) %>%
  group_by(scenario) %>%
  slice(1) %>%
  mutate(
    PVR_total = PVR_fixed + PVR_var,
    delta_6MWD = sixMWD - 342,
    delta_mPAP = mPAP - 46
  ) %>%
  select(scenario, mPAP, PVR_total, CO, SaO2, BNP, sixMWD,
         delta_6MWD, delta_mPAP, WHO_FC) %>%
  ungroup()

message("\n=== Week 52 Summary ===")
print(as.data.frame(summary_table), digits = 3)

## ============================================================
## CLINICAL VALIDATION NOTE
## ============================================================
message("
Clinical Calibration Notes:
- CHEST-1 trial: Riociguat 2.5mg TID -> +46m 6MWD at week 16 (Ghofrani 2013)
- SERAPHIN: Macitentan 10mg -> -47% risk of clinical worsening (Pulido 2013)
- MERIT: Combination ERA + sGC in CTEPH inoperable (Meyer 2022)
- CTREPH: Riociguat reduced PVR by 31% vs placebo (Ghofrani 2014)
- PEA: Mean PVR reduction ~65-75% in eligible patients (Madani 2012)
- BPA: PVR reduction ~50% after series of sessions (Wiedenroth 2020)
- mPAP threshold for CTEPH diagnosis: >20 mmHg (ESC/ERS 2022 updated)
")
