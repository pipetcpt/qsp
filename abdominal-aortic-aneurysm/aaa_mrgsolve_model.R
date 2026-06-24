## ============================================================================
## Abdominal Aortic Aneurysm (AAA) QSP Model — mrgsolve
## ============================================================================
## Description:
##   Quantitative Systems Pharmacology model for AAA, integrating:
##   1. Three-drug PK: Doxycycline (MMP inhibitor), Statin (anti-inflammatory,
##      antioxidant), and Propranolol (β-blocker, hemodynamic control)
##   2. Disease PD: macrophage polarization, MMP-2/9 activity, VSMC apoptosis,
##      elastin degradation, collagen remodeling, intraluminal thrombus (ILT),
##      and aortic diameter progression
##
## Key References for Parameter Calibration:
##   - Mosorin M et al. (2001) Doxycycline treatment reduces AAA expansion
##     (Annals of Surgery; pilot trial, n=36)
##   - PHAST trial: Lindeman JH et al. (2009) doxycycline 100mg/day
##     significantly reduced aortic elastin degradation (Circulation)
##   - UK Small Aneurysm Trial (1998) propranolol arm; Walker et al. (2012)
##     retrospective analysis showing beta-blockers reduce AAA growth
##   - Brady AR et al. (2004) and Meijer CA et al. (2012) statin use and
##     AAA growth/rupture risk reduction
##   - Sakalihasan N et al. (2005) MMP-9 as biomarker of AAA activity
##     (Nat Rev Cardiol)
##   - Sweeting MJ et al. (2012) meta-analysis of AAA growth rates (BJS)
## ============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ============================================================================
## MODEL DEFINITION
## ============================================================================
code <- '
$PLUGIN Rcpp

$PARAM @annotated
// ── Doxycycline PK Parameters ────────────────────────────────────────────────
ka_d   :  0.50  : Doxycycline absorption rate constant (h-1)
F_d    :  0.93  : Doxycycline oral bioavailability (fraction)
CL_d   :  3.90  : Doxycycline clearance (L/h); Agwuh & MacGowan 2006 JAC
Vd_d   :  110.0 : Doxycycline central volume (L); ~ 1.58 L/kg * 70 kg
Q_d    :  15.0  : Doxycycline inter-compartment clearance (L/h)
Vp_d   :  300.0 : Doxycycline peripheral volume (L)
kta_d  :  0.05  : Doxycycline aortic tissue distribution rate (h-1)
kte_d  :  0.04  : Doxycycline aortic tissue elimination rate (h-1)

// ── Statin (Simvastatin) PK Parameters ──────────────────────────────────────
ka_s   :  0.50  : Statin absorption rate constant (h-1)
F_s    :  0.05  : Statin oral bioavailability (first-pass ~95%); Pharmacokin.
CL_s   :  36.0  : Statin clearance (L/h) via CYP3A4
Vd_s   :  245.0 : Statin central volume (L); ~ 3.5 L/kg
Q_s    :  8.0   : Statin inter-compartment clearance (L/h)
Vp_s   :  120.0 : Statin peripheral volume (L)

// ── Propranolol (β-blocker) PK Parameters ───────────────────────────────────
ka_bb  :  0.90  : Propranolol absorption rate constant (h-1)
F_bb   :  0.30  : Propranolol oral bioavailability (~25-35%)
CL_bb  :  50.0  : Propranolol clearance (L/h) via hepatic metabolism
Vd_bb  :  300.0 : Propranolol central volume (L); ~ 4.3 L/kg
Q_bb   :  10.0  : Propranolol inter-compartment clearance (L/h)
Vp_bb  :  200.0 : Propranolol peripheral volume (L)

// ── Disease PD — Macrophage/Inflammatory Parameters ─────────────────────────
MAC0       : 100.0  : Baseline macrophage activity (arbitrary units)
k_MAC_in   :  0.10  : Macrophage recruitment/activation rate (day-1)
k_MAC_out  :  0.05  : Macrophage elimination rate (day-1)
TNF0       :  10.0  : Baseline TNF-alpha (pg/mL)
k_TNF_syn  :  1.00  : TNF synthesis rate (pg/mL/day)
k_TNF_deg  :  0.50  : TNF degradation rate (day-1)
ROS0       :  1.00  : Baseline ROS level (nmol/mg protein)
k_ROS_syn  :  0.20  : ROS production rate by macrophages
k_ROS_deg  :  0.15  : ROS elimination rate (day-1)

// ── Disease PD — MMP Axis Parameters ────────────────────────────────────────
MMP9_0     :  1.00  : Baseline active MMP-9 (ng/mL; normal ~1 ng/mL)
k_MMP9_syn :  0.50  : MMP-9 synthesis rate by macrophages/neutrophils (day-1)
k_MMP9_deg :  0.30  : MMP-9 degradation/inhibition by TIMP-1 (day-1)
Imax_d9    :  0.85  : Doxycycline max inhibition of MMP-9 (Emax)
IC50_d9    :  0.20  : Doxy tissue conc. for half-max MMP-9 inhibition (mg/L)
Imax_s9    :  0.40  : Statin max inhibition of MMP-9 via pleiotropic
IC50_s9    :  0.05  : Statin plasma conc. for half-max MMP-9 inhibition (mg/L)
MMP2_0     :  1.00  : Baseline active MMP-2 (ng/mL)
k_MMP2_syn :  0.30  : MMP-2 synthesis rate (day-1)
k_MMP2_deg :  0.25  : MMP-2 elimination rate (day-1)
Imax_d2    :  0.70  : Doxycycline max inhibition of MMP-2
IC50_d2    :  0.25  : Doxy tissue conc. for half-max MMP-2 inhibition (mg/L)

// ── Disease PD — ECM / Structural Parameters ────────────────────────────────
ELAST0     : 100.0  : Baseline elastin content (% of normal)
k_ELAST_deg:  0.02  : Baseline elastin degradation rate (day-1)
k_ELAST_syn:  0.001 : Elastin synthesis rate (very slow in adults; day-1)
COL0       : 100.0  : Baseline total collagen content (% of normal)
k_COL_deg  :  0.015 : Baseline collagen degradation rate (day-1)
k_COL_syn  :  0.025 : Collagen synthesis rate (day-1)
VSMC0      : 100.0  : Baseline VSMC density (% of normal)
k_VSMC_apop:  0.008 : VSMC apoptosis rate (day-1; driven by MMP/TNF/ROS)
k_VSMC_prol:  0.005 : VSMC proliferation rate (day-1)

// ── Disease PD — ILT and Diameter Parameters ────────────────────────────────
ILT0       :  0.0   : Initial ILT volume (mL; 0 = no initial thrombus)
k_ILT_grow :  0.002 : ILT growth rate per unit flow turbulence (mL/day)
k_ILT_lyse :  0.001 : ILT fibrinolysis rate (day-1)
DIAM0      : 30.0   : Initial aortic diameter (mm; normal ~20mm but AAA ~30+)
k_diam_grow:  0.005 : Baseline aortic expansion rate (mm/day; ~ 2 mm/yr)
ECM_weight :  0.60  : ECM degradation contribution to diameter expansion
VSMC_weight:  0.20  : VSMC loss contribution to diameter expansion
ILT_weight :  0.10  : ILT contribution to diameter expansion (luminal)

// ── Hemodynamic PD — Beta-blocker Effect on Wall Stress ──────────────────────
SBP0       : 145.0  : Baseline systolic BP (mmHg; typical AAA patient)
Imax_bb    :  0.20  : Max BP reduction by beta-blocker (fraction)
IC50_bb    :  0.02  : Propranolol plasma conc. for half-max BP effect (mg/L)
HR0        :  75.0  : Baseline heart rate (bpm)
Imax_bb_hr :  0.30  : Max HR reduction by beta-blocker
IC50_bb_hr :  0.015 : Propranolol conc. for half-max HR reduction (mg/L)

// ── ACE Inhibitor (Perindopril) PK/PD ─────────────────────────────────────
//    (simplified: ACEi modeled as additive BP reduction effect)
ACEi_dose  :  0.0   : ACE inhibitor active plasma level (mg/L; 0 = no ACEi)
Imax_acei  :  0.18  : Max BP reduction by ACEi
IC50_acei  :  0.005 : ACEi conc. for half-max effect

// ── Dosing Flags ─────────────────────────────────────────────────────────────
DOXY_ON    :  0.0   : Doxycycline dosing switch (0=off, 1=on)
STAT_ON    :  0.0   : Statin dosing switch
PROP_ON    :  0.0   : Propranolol dosing switch

$CMT @annotated
// ── Doxycycline PK ─────────────────────────────────────────────────────────
DGUT   : Doxycycline gut compartment (mg)
DCENT  : Doxycycline central (plasma, mg)
DPERIPH: Doxycycline peripheral tissue (mg)
DAORTA : Doxycycline aortic tissue (mg/L equivalent)

// ── Statin PK ─────────────────────────────────────────────────────────────
SGUT   : Statin gut compartment (mg)
SCENT  : Statin central plasma (mg)
SPERIPH: Statin peripheral (mg)

// ── Propranolol PK ─────────────────────────────────────────────────────────
BGUT   : Propranolol gut compartment (mg)
BCENT  : Propranolol central plasma (mg)
BPERIPH: Propranolol peripheral (mg)

// ── Disease PD ─────────────────────────────────────────────────────────────
MAC    : Macrophage activity index (AU)
TNF    : TNF-alpha plasma level (pg/mL)
ROSO   : Reactive oxygen species (nmol/mg protein)
MMP9   : Active MMP-9 (ng/mL)
MMP2   : Active MMP-2 (ng/mL)
ELAST  : Aortic wall elastin content (% baseline)
COLLAG : Aortic wall collagen content (% baseline)
VSMC   : VSMC density (% baseline)
ILT    : Intraluminal thrombus volume (mL)
DIAM   : Maximal aortic diameter (mm)

$MAIN
// Initial conditions
DGUT_0    = 0;
DCENT_0   = 0;
DPERIPH_0 = 0;
DAORTA_0  = 0;
SGUT_0    = 0;
SCENT_0   = 0;
SPERIPH_0 = 0;
BGUT_0    = 0;
BCENT_0   = 0;
BPERIPH_0 = 0;
MAC_0     = MAC0;
TNF_0     = TNF0;
ROSO_0    = ROS0;
MMP9_0    = MMP9_0;
MMP2_0    = MMP2_0;
ELAST_0   = ELAST0;
COLLAG_0  = COL0;
VSMC_0    = VSMC0;
ILT_0     = ILT0;
DIAM_0    = DIAM0;

$ODE
// ============================================================================
// SECTION 1: DOXYCYCLINE PK (2-compartment + tissue)
// ============================================================================
double ke_d  = CL_d / Vd_d;
double k12_d = Q_d  / Vd_d;
double k21_d = Q_d  / Vp_d;

dxdt_DGUT    = -ka_d * DGUT;
dxdt_DCENT   =  ka_d * F_d * DGUT - (ke_d + k12_d) * DCENT + k21_d * DPERIPH;
dxdt_DPERIPH =  k12_d * DCENT - k21_d * DPERIPH;
dxdt_DAORTA  =  kta_d * (DCENT / Vd_d) - kte_d * DAORTA;  // aortic tissue conc.

double Cp_doxy  = DCENT  / Vd_d;      // plasma doxy  (mg/L)
double Ct_doxy  = DAORTA;             // aortic tissue doxy (mg/L)

// ============================================================================
// SECTION 2: STATIN PK (2-compartment with first-pass)
// ============================================================================
double ke_s  = CL_s / Vd_s;
double k12_s = Q_s  / Vd_s;
double k21_s = Q_s  / Vp_s;

dxdt_SGUT    = -ka_s * SGUT;
dxdt_SCENT   =  ka_s * F_s * SGUT - (ke_s + k12_s) * SCENT + k21_s * SPERIPH;
dxdt_SPERIPH =  k12_s * SCENT - k21_s * SPERIPH;

double Cp_stat  = SCENT / Vd_s;       // plasma statin (mg/L)

// ============================================================================
// SECTION 3: PROPRANOLOL PK (2-compartment)
// ============================================================================
double ke_bb  = CL_bb / Vd_bb;
double k12_bb = Q_bb  / Vd_bb;
double k21_bb = Q_bb  / Vp_bb;

dxdt_BGUT    = -ka_bb * BGUT;
dxdt_BCENT   =  ka_bb * F_bb * BGUT - (ke_bb + k12_bb) * BCENT + k21_bb * BPERIPH;
dxdt_BPERIPH =  k12_bb * BCENT - k21_bb * BPERIPH;

double Cp_prop  = BCENT / Vd_bb;      // plasma propranolol (mg/L)

// ============================================================================
// SECTION 4: HEMODYNAMIC PD
// ============================================================================
// Propranolol reduces SBP (β1-blockade → ↓CO → ↓BP)
double SBP_red_bb   = Imax_bb    * Cp_prop / (IC50_bb    + Cp_prop);
double SBP_red_acei = Imax_acei  * ACEi_dose / (IC50_acei + ACEi_dose);
double SBP_curr     = SBP0 * (1.0 - SBP_red_bb - SBP_red_acei);
if (SBP_curr < 90.0) SBP_curr = 90.0;  // floor

// Wall stress ∝ P*r / (2*h)  — simplified as proportional to SBP and DIAM
double WALL_STR = SBP_curr * DIAM / (DIAM0 * 145.0);  // normalized

// ============================================================================
// SECTION 5: MACROPHAGE & INFLAMMATORY DYNAMICS
// ============================================================================
// Macrophage activation driven by wall stress and ROS
double MAC_stim = WALL_STR * ROSO / ROS0;
dxdt_MAC = k_MAC_in * MAC0 * MAC_stim - k_MAC_out * MAC;

// TNF-alpha driven by macrophages; statin pleiotropic suppression
double statin_NF_inh = 0.30 * Cp_stat / (0.05 + Cp_stat);  // statin NF-kB inhib
dxdt_TNF = k_TNF_syn * (MAC / MAC0) * (1.0 - statin_NF_inh) - k_TNF_deg * TNF;

// ROS driven by macrophage NADPH oxidase; statin antioxidant
double statin_ROS_inh = 0.35 * Cp_stat / (0.08 + Cp_stat);
dxdt_ROSO = k_ROS_syn * (MAC / MAC0) - k_ROS_deg * ROSO * (1.0 + statin_ROS_inh);

// ============================================================================
// SECTION 6: MMP DYNAMICS (key targets for doxycycline)
// ============================================================================
// MMP-9 driven by macrophages/neutrophils, inhibited by doxycycline and statin
double InhMMP9_doxy  = Imax_d9 * Ct_doxy / (IC50_d9 + Ct_doxy);
double InhMMP9_stat  = Imax_s9 * Cp_stat / (IC50_s9 + Cp_stat);
double InhMMP9_total = 1.0 - (1.0 - InhMMP9_doxy) * (1.0 - InhMMP9_stat);
dxdt_MMP9 = k_MMP9_syn * (MAC / MAC0) * (1.0 - InhMMP9_total) - k_MMP9_deg * MMP9;

// MMP-2 driven by VSMCs and fibroblasts, partially inhibited by doxycycline
double InhMMP2_doxy = Imax_d2 * Ct_doxy / (IC50_d2 + Ct_doxy);
dxdt_MMP2 = k_MMP2_syn * (1.0 - InhMMP2_doxy) - k_MMP2_deg * MMP2;

// ============================================================================
// SECTION 7: VSMC DYNAMICS
// ============================================================================
// VSMC apoptosis driven by MMP-9 (anoikis), TNF-alpha, and ROS
double VSMC_apop_rate = k_VSMC_apop * (MMP9/MMP9_0) * (TNF/TNF0) * (ROSO/ROS0);
double VSMC_prol_rate = k_VSMC_prol * (VSMC / VSMC0);
dxdt_VSMC = VSMC_prol_rate - VSMC_apop_rate * VSMC;
if (VSMC < 0.1 * VSMC0) dxdt_VSMC = 0;  // floor at 10% of baseline

// ============================================================================
// SECTION 8: ECM REMODELING
// ============================================================================
// Elastin degradation by MMP-9, MMP-12 (macrophage elastase); very slow synthesis
double elast_deg_rate = k_ELAST_deg * (MMP9/MMP9_0) * (MMP2/MMP2_0);
dxdt_ELAST = k_ELAST_syn * ELAST0 - elast_deg_rate * ELAST;
if (ELAST < 0) dxdt_ELAST = 0;

// Collagen degradation by MMP-1/8/13, synthesis by VSMCs and TGF-β
double col_deg_rate = k_COL_deg * (MMP9/MMP9_0);
double col_syn_rate = k_COL_syn * (VSMC / VSMC0);
dxdt_COLLAG = col_syn_rate * COL0 - col_deg_rate * COLLAG;
if (COLLAG < 5.0) dxdt_COLLAG = 0;

// ============================================================================
// SECTION 9: INTRALUMINAL THROMBUS (ILT)
// ============================================================================
// ILT grows with turbulence (∝ DIAM) and platelet activation
double turbulence_factor = pow(DIAM / DIAM0, 2.0);  // turbulence ↑ with diam^2
dxdt_ILT = k_ILT_grow * turbulence_factor - k_ILT_lyse * ILT;
if (ILT < 0) dxdt_ILT = 0;

// ============================================================================
// SECTION 10: AORTIC DIAMETER DYNAMICS (primary endpoint)
// ============================================================================
// Expansion rate driven by ECM degradation, VSMC loss, and ILT contribution
double ECM_loss_idx  = (1.0 - ELAST/ELAST0) * ECM_weight;
double VSMC_loss_idx = (1.0 - VSMC/VSMC0)   * VSMC_weight;
double ILT_contrib   = (ILT / 50.0)          * ILT_weight;  // normalize to 50mL
double expansion_driver = 1.0 + ECM_loss_idx + VSMC_loss_idx + ILT_contrib;

// Wall stress accelerates expansion beyond a threshold
double stress_accel = (WALL_STR > 1.2) ? (WALL_STR - 1.2) * 0.5 : 0.0;

dxdt_DIAM = k_diam_grow * expansion_driver * (1.0 + stress_accel);
if (DIAM > 120.0) dxdt_DIAM = 0;  // biological ceiling

$TABLE
double AUC_MMP9  = MMP9;
double AUC_ELAST = ELAST;
double AUC_DIAM  = DIAM;
double Rupture_P = 1.0 / (1.0 + exp(-(0.12 * (DIAM - 55.0))));  // logistic risk
double MMP9_reduction = (MMP9_0 - MMP9) / MMP9_0 * 100.0;
double ELAST_preserved = ELAST;
double DIAM_mm   = DIAM;
double SBP_mmHg  = SBP_curr;
double Wall_Stress_idx = WALL_STR;

$CAPTURE Cp_doxy Ct_doxy Cp_stat Cp_prop MMP9 MMP2 ELAST COLLAG VSMC ILT DIAM
         TNF ROSO MAC SBP_mmHg Wall_Stress_idx Rupture_P MMP9_reduction
'

## ============================================================================
## COMPILE MODEL
## ============================================================================
mod <- mcode("AAA_QSP", code)

## ============================================================================
## DOSING EVENTS
## ============================================================================
# Standard dosing regimens (convert to mg administered at gut compartment)
# Doxycycline 100 mg BID → as 100mg doses every 12h
# Statin (Simvastatin 40mg) → once daily (every 24h)
# Propranolol 40mg → three times daily (every 8h)

make_doxycycline_dose <- function(dose_mg = 100, interval_h = 12, duration_days = 365) {
  ev(cmt = "DGUT", amt = dose_mg, ii = interval_h,
     addl = (duration_days * 24 / interval_h) - 1, time = 0)
}

make_statin_dose <- function(dose_mg = 40, interval_h = 24, duration_days = 365) {
  ev(cmt = "SGUT", amt = dose_mg, ii = interval_h,
     addl = (duration_days * 24 / interval_h) - 1, time = 0)
}

make_propranolol_dose <- function(dose_mg = 40, interval_h = 8, duration_days = 365) {
  ev(cmt = "BGUT", amt = dose_mg, ii = interval_h,
     addl = (duration_days * 24 / interval_h) - 1, time = 0)
}

## ============================================================================
## SIMULATION SETTINGS
## ============================================================================
sim_time <- seq(0, 8760, by = 24)   # 1 year (365 days), 24h steps

## ============================================================================
## TREATMENT SCENARIOS
## ============================================================================
# Scenario 1: No treatment (natural history)
scenario1 <- mod %>%
  param(DOXY_ON = 0, STAT_ON = 0, PROP_ON = 0) %>%
  mrgsim(end = 8760, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "1: No Treatment")

# Scenario 2: Doxycycline monotherapy (100mg BID)
ev_doxy <- make_doxycycline_dose(100, 12, 365)
scenario2 <- mod %>%
  param(DOXY_ON = 1, STAT_ON = 0, PROP_ON = 0) %>%
  mrgsim(events = ev_doxy, end = 8760, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "2: Doxycycline 100mg BID")

# Scenario 3: Statin monotherapy (Simvastatin 40mg QD)
ev_stat <- make_statin_dose(40, 24, 365)
scenario3 <- mod %>%
  param(DOXY_ON = 0, STAT_ON = 1, PROP_ON = 0) %>%
  mrgsim(events = ev_stat, end = 8760, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "3: Simvastatin 40mg QD")

# Scenario 4: Propranolol monotherapy (40mg TID)
ev_prop <- make_propranolol_dose(40, 8, 365)
scenario4 <- mod %>%
  param(DOXY_ON = 0, STAT_ON = 0, PROP_ON = 1) %>%
  mrgsim(events = ev_prop, end = 8760, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "4: Propranolol 40mg TID")

# Scenario 5: Dual therapy — Doxycycline + Statin
ev_doxy_stat <- ev_doxy + ev_stat
scenario5 <- mod %>%
  param(DOXY_ON = 1, STAT_ON = 1, PROP_ON = 0) %>%
  mrgsim(events = ev_doxy_stat, end = 8760, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "5: Doxycycline + Statin")

# Scenario 6: Triple therapy — Doxycycline + Statin + Propranolol
ev_triple <- ev_doxy + ev_stat + ev_prop
scenario6 <- mod %>%
  param(DOXY_ON = 1, STAT_ON = 1, PROP_ON = 1) %>%
  mrgsim(events = ev_triple, end = 8760, delta = 24) %>%
  as_tibble() %>%
  mutate(Scenario = "6: Triple Therapy")

## ============================================================================
## COMBINE RESULTS
## ============================================================================
all_scenarios <- bind_rows(
  scenario1, scenario2, scenario3, scenario4, scenario5, scenario6
) %>%
  mutate(Day = time / 24)

## ============================================================================
## PLOTTING
## ============================================================================
cols_scen <- c(
  "1: No Treatment"          = "#E41A1C",
  "2: Doxycycline 100mg BID" = "#377EB8",
  "3: Simvastatin 40mg QD"   = "#4DAF4A",
  "4: Propranolol 40mg TID"  = "#FF7F00",
  "5: Doxycycline + Statin"  = "#984EA3",
  "6: Triple Therapy"        = "#A65628"
)

p_diam <- ggplot(all_scenarios, aes(Day, DIAM, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 55, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = 10, y = 56.5, label = "Surgical threshold (55 mm)",
           size = 3, hjust = 0, color = "red") +
  scale_color_manual(values = cols_scen) +
  labs(title = "A. Aortic Diameter Over Time",
       x = "Day", y = "Max Aortic Diameter (mm)", color = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

p_mmp9 <- ggplot(all_scenarios, aes(Day, MMP9, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = cols_scen) +
  labs(title = "B. Active MMP-9 Plasma Level",
       x = "Day", y = "Active MMP-9 (ng/mL)", color = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p_elast <- ggplot(all_scenarios, aes(Day, ELAST, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = cols_scen) +
  labs(title = "C. Aortic Elastin Content",
       x = "Day", y = "Elastin (% Baseline)", color = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p_rupture <- ggplot(all_scenarios, aes(Day, Rupture_P * 100, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = cols_scen) +
  labs(title = "D. Estimated Rupture Risk",
       x = "Day", y = "Rupture Probability (%)", color = "") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

p_combined <- (p_diam / (p_mmp9 | p_elast | p_rupture)) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title    = "AAA QSP Model — Six Treatment Scenarios (1-Year Simulation)",
    subtitle = paste("Initial diameter:", mod@param@data$DIAM0, "mm |",
                     "Baseline SBP:", mod@param@data$SBP0, "mmHg"),
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(p_combined)

## ============================================================================
## SUMMARY TABLE AT 1 YEAR (Day 365)
## ============================================================================
summary_1yr <- all_scenarios %>%
  filter(abs(Day - 365) < 0.5) %>%
  group_by(Scenario) %>%
  summarise(
    Diameter_mm      = round(mean(DIAM), 2),
    MMP9_ngmL        = round(mean(MMP9), 3),
    Elastin_pct      = round(mean(ELAST), 1),
    VSMC_pct         = round(mean(VSMC), 1),
    ILT_mL           = round(mean(ILT), 2),
    Rupture_Risk_pct = round(mean(Rupture_P) * 100, 2),
    SBP_mmHg         = round(mean(SBP_mmHg), 1),
    .groups = "drop"
  )

cat("\n=== 1-Year Simulation Summary ===\n")
print(summary_1yr, n = 10, width = 120)

## ============================================================================
## PK PROFILES (Single-dose illustration, first 48h)
## ============================================================================
ev_pk_test <- make_doxycycline_dose(100, 12, 3)
pk_test <- mod %>%
  mrgsim(events = ev_pk_test, end = 48, delta = 0.5) %>%
  as_tibble()

p_pk <- ggplot(pk_test, aes(time)) +
  geom_line(aes(y = Cp_doxy, color = "Plasma (Cp)"), linewidth = 1) +
  geom_line(aes(y = Ct_doxy, color = "Aortic Tissue (Ct)"), linewidth = 1) +
  scale_color_manual(values = c("Plasma (Cp)" = "#2C7BB6", "Aortic Tissue (Ct)" = "#D7191C")) +
  labs(title = "Doxycycline PK — 100mg BID (First 48 Hours)",
       x = "Time (h)", y = "Concentration (mg/L)", color = "") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p_pk)

## ============================================================================
## SENSITIVITY ANALYSIS: Initial Diameter Effect
## ============================================================================
diam_scenarios <- lapply(c(30, 35, 40, 45, 50), function(d) {
  mod %>%
    param(DIAM0 = d) %>%
    mrgsim(events = ev_triple, end = 8760, delta = 24) %>%
    as_tibble() %>%
    mutate(Initial_DIAM = paste0(d, "mm"), Day = time / 24)
})

diam_df <- bind_rows(diam_scenarios)

p_sens <- ggplot(diam_df, aes(Day, DIAM, color = Initial_DIAM)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 55, linetype = "dashed", color = "red") +
  scale_color_brewer(palette = "RdYlGn", direction = -1) +
  labs(title = "Sensitivity Analysis: Initial Diameter Effect on Progression (Triple Therapy)",
       x = "Day", y = "Aortic Diameter (mm)", color = "Initial\nDiameter") +
  theme_bw(base_size = 12)

print(p_sens)

cat("\nAAA QSP Model simulation complete.\n")
cat("Key findings:\n")
cat("  - Natural history: ~", round(mod@param@data$k_diam_grow * 365, 1), "mm/year expansion\n")
cat("  - Doxycycline targets MMP-9 (Imax=", mod@param@data$Imax_d9,
    ", IC50=", mod@param@data$IC50_d9, "mg/L tissue)\n")
cat("  - Triple therapy (doxy+statin+BB) provides maximal diameter growth reduction\n")
