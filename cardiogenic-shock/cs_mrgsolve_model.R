# =============================================================================
# Cardiogenic Shock (CS) — QSP / mrgsolve model
# -----------------------------------------------------------------------------
# Mechanistic PK/PD model integrating:
#   * Hemodynamic compartments (LVEF, CO, MAP, SVR, PCWP, lactate, SvO2)
#   * Sympathetic / RAAS / inflammation / iNOS-NO surge (vasoplegia)
#   * Multi-organ failure markers (creatinine, ALT, urine output)
#   * PK/PD for: norepinephrine, epinephrine, dobutamine, milrinone,
#     levosimendan (+ OR-1896 active metabolite), vasopressin
#   * Mechanical circulatory support (MCS) effects: IABP, Impella, VA-ECMO
#   * Composite mortality hazard driven by SCAI stage proxy
#
# Calibration references (qualitative):
#   - SHOCK trial (Hochman NEJM 1999/JAMA 2006)
#   - IABP-SHOCK II (Thiele NEJM 2012; Lancet 2019)
#   - SOAP-II dopamine vs NE (De Backer NEJM 2010)
#   - DanGer Shock — Impella CP in STEMI-CS (Møller NEJM 2024)
#   - ECLS-SHOCK — VA-ECMO in AMI-CS (Thiele NEJM 2023)
#   - DOREMI dobutamine vs milrinone (Mathew NEJM 2021)
#   - SURVIVE / LIDO levosimendan trials
#
# Five+ therapeutic scenarios are provided (`run_cs_scenarios()` at the end).
# =============================================================================

library(mrgsolve)
library(dplyr)
library(tibble)
library(ggplot2)

cs_code <- '
$PROB
# Cardiogenic Shock QSP model (LV/RV failure + neurohormonal +
# inflammation + microcirculation + drug PK/PD + MCS).

$PARAM @annotated
// ---- Patient baseline ----
BW       :  75   : body weight (kg)
BSA      :   1.8 : body surface area (m^2)
LVEF0    :   0.22: baseline LVEF (fraction)
CO0      :   2.0 : baseline CO (L/min, shock state)
MAP0     :  62   : baseline MAP (mmHg)
SVR0     : 1500  : baseline SVR (dyn*s*cm^-5)
PCWP0    :  24   : baseline PCWP (mmHg)
HR0      : 110   : baseline HR (bpm)
LAC0     :   4.0 : baseline lactate (mmol/L)
SVO20    :   0.55: baseline SvO2 (fraction)
CR0      :   1.2 : baseline creatinine (mg/dL)
ALT0     :  40   : baseline ALT (U/L)
UO0      :   0.4 : baseline urine output (mL/kg/h)

// ---- Drug PK / receptor binding ----
// Norepinephrine
NE_CL    :  30   : NE clearance (L/h)
NE_V     :   8   : NE Vd (L)
NE_KD_a1 :   0.04: NE EC50 alpha1 (µg/min)  (effect on SVR)
NE_KD_b1 :   0.3 : NE EC50 beta1 (µg/min)
NE_EMAX_a1:  1.0
NE_EMAX_b1:  0.4

// Epinephrine
EPI_CL   :  90
EPI_V    :  35
EPI_KD_a1:   0.06
EPI_KD_b1:   0.05
EPI_KD_b2:   0.04

// Dobutamine
DOBU_CL  : 120
DOBU_V   :  20
DOBU_KD_b1: 2.5      // µg/kg/min for half-max
DOBU_KD_b2: 6.0

// Dopamine
DA_CL    :  60
DA_V     :  15
DA_EC50_d1: 3        // µg/kg/min (renal vasodil)
DA_EC50_a1: 12

// Milrinone (PDE3 inhibitor)
MIL_CL   :   8      // L/h
MIL_V    :  35
MIL_IC50 :   0.15   // µg/mL (PDE3 inhibition)
MIL_EMAX :   0.9

// Levosimendan + OR-1896
LEVO_CL  :  17
LEVO_V   :  12
LEVO_KM  :   0.05    // OR-1896 conversion rate (1/h)
OR_CL    :   0.5
OR_V     :  20
LEVO_TNC :   0.45    // Ca-sensitization Emax
OR_TNC   :   0.30

// Vasopressin
VAS_CL   :  45
VAS_V    :   7
VAS_EC50 :   0.030   // U/min
VAS_EMAX :   0.4     // fractional SVR rise

// ---- Pathophysiology rate constants ----
k_stun_recov : 0.005 // stunning recovery per hour (long)
k_isch_dmg   : 0.01  // ischemic damage rate (per hr if perfusion low)
k_remod      : 0.002 // remodeling per day

k_SVR_synth  : 0.5
k_SVR_deg    : 0.4
k_NO_synth   : 0.08
k_NO_deg     : 0.15
k_NO_indEX   : 0.6   // factor amplifying NO under TNFa
k_TNF_synth  : 0.05  // TNF synthesis (DAMP driven)
k_TNF_deg    : 0.2
k_DAMP_synth : 0.2
k_DAMP_deg   : 0.1

k_lac_form   : 0.6   // lactate production from low DO2
k_lac_clear  : 0.3
k_SvO2_norm  : 0.7   // target SvO2 normalized
k_Cr_rise    : 0.04  // creatinine rise under AKI
k_Cr_clear   : 0.025
k_ALT_rise   : 0.10
k_ALT_clear  : 0.05
k_UO_baseline: 0.4

// ---- MCS support flags & magnitudes ----
IABP_ON    : 0
IABP_CO    : 0.5     // L/min augment
Impella_ON : 0
Impella_CO : 3.5     // L/min direct LV unload
Impella_unload : 0.25
ECMO_ON    : 0
ECMO_CO    : 4.0
ECMO_unload: -0.10   // negative = INCREASES LV afterload risk

// ---- Drug infusions ----
NE_inf     : 0       // µg/min
EPI_inf    : 0       // µg/min
DOBU_inf   : 0       // µg/kg/min
DA_inf     : 0       // µg/kg/min
MIL_inf    : 0       // µg/kg/min
LEVO_inf   : 0       // µg/kg/min
VAS_inf    : 0       // U/min

// Mortality hazard coefficients
h_base     : 0.001   // /h baseline
b_lac      : 0.05
b_MAP      : 0.04
b_CO       : 0.20
b_NOex     : 0.30

$CMT @annotated
// Drug compartments (8)
NE_amt   : NE amount  (µg)
EPI_amt  : Epi amount
DOBU_amt : Dobutamine
DA_amt   : Dopamine
MIL_amt  : Milrinone
LEVO_amt : Levosimendan parent
OR_amt   : OR-1896 active metabolite
VAS_amt  : Vasopressin

// Hemodynamic / physiologic state (12)
EF       : LV ejection fraction
COx      : Cardiac output (L/min)
MAP      : mean arterial pressure (mmHg)
SVR      : systemic vascular resistance (state)
PCWP     : pulmonary cap wedge
HR       : heart rate (bpm)
LAC      : lactate (mmol/L)
SVO2     : mixed venous SaO2 (fraction)
CR       : creatinine (mg/dL)
ALT      : ALT (U/L)
UO       : urine output (mL/kg/h)
STUN     : stunned myocardium fraction (0-1)

// Pathophysiology mediators (6)
NOex     : NO excess index
TNFa     : TNF-alpha
DAMP     : DAMPs index
RAAS     : RAAS activity index
SNS      : sympathetic activity
CUMHAZ   : cumulative mortality hazard

$MAIN
// Initial conditions
EF_0       = LVEF0;
COx_0      = CO0;
MAP_0      = MAP0;
SVR_0      = SVR0;
PCWP_0     = PCWP0;
HR_0       = HR0;
LAC_0      = LAC0;
SVO2_0     = SVO20;
CR_0       = CR0;
ALT_0      = ALT0;
UO_0       = UO0;
STUN_0     = 0.4;
NOex_0     = 0.3;
TNFa_0     = 1.0;
DAMP_0     = 1.0;
RAAS_0     = 1.5;
SNS_0      = 2.0;
CUMHAZ_0   = 0;

$ODE
// ---- Drug PK ----
double NE_conc   = NE_amt   / NE_V;       // µg/L plasma proxy
double EPI_conc  = EPI_amt  / EPI_V;
double DOBU_conc = DOBU_amt / DOBU_V;
double DA_conc   = DA_amt   / DA_V;
double MIL_conc  = MIL_amt  / MIL_V;
double LEVO_conc = LEVO_amt / LEVO_V;
double OR_conc   = OR_amt   / OR_V;
double VAS_conc  = VAS_amt  / VAS_V;

// First-order elimination (CL*C)
dxdt_NE_amt   = NE_inf*60.0   - NE_CL   * NE_conc;
dxdt_EPI_amt  = EPI_inf*60.0  - EPI_CL  * EPI_conc;
dxdt_DOBU_amt = DOBU_inf*BW*60 - DOBU_CL* DOBU_conc;
dxdt_DA_amt   = DA_inf*BW*60   - DA_CL  * DA_conc;
dxdt_MIL_amt  = MIL_inf*BW*60  - MIL_CL * MIL_conc;
dxdt_LEVO_amt = LEVO_inf*BW*60 - LEVO_CL* LEVO_conc - LEVO_KM*LEVO_amt;
dxdt_OR_amt   = LEVO_KM*LEVO_amt - OR_CL*OR_conc;
dxdt_VAS_amt  = VAS_inf*60    - VAS_CL  * VAS_conc;

// ---- Receptor occupancy / drug effects ----
double E_NE_a1   = NE_EMAX_a1  * NE_inf  / (NE_KD_a1  + NE_inf  + 1e-9);
double E_NE_b1   = NE_EMAX_b1  * NE_inf  / (NE_KD_b1  + NE_inf  + 1e-9);
double E_EPI_a1  = 1.0         * EPI_inf / (EPI_KD_a1 + EPI_inf + 1e-9);
double E_EPI_b1  = 0.8         * EPI_inf / (EPI_KD_b1 + EPI_inf + 1e-9);
double E_EPI_b2  = 0.6         * EPI_inf / (EPI_KD_b2 + EPI_inf + 1e-9);
double E_DOBU_b1 = 0.9         * DOBU_inf/ (DOBU_KD_b1+ DOBU_inf+ 1e-9);
double E_DOBU_b2 = 0.4         * DOBU_inf/ (DOBU_KD_b2+ DOBU_inf+ 1e-9);
double E_DA_d1   = 0.4         * DA_inf  / (DA_EC50_d1+ DA_inf  + 1e-9);
double E_DA_a1   = 0.9         * pow(DA_inf,2) / (pow(DA_EC50_a1,2) + pow(DA_inf,2) + 1e-9);
double E_MIL     = MIL_EMAX    * MIL_conc/ (MIL_IC50 + MIL_conc + 1e-9);
double E_LEVO    = LEVO_TNC    * LEVO_conc /(0.05 + LEVO_conc + 1e-9);
double E_OR      = OR_TNC      * OR_conc /(0.02 + OR_conc + 1e-9);
double E_VAS     = VAS_EMAX    * VAS_inf  /(VAS_EC50 + VAS_inf + 1e-9);

// Combined inotropy (Ca-sensitization + cAMP-PDE3 + beta1 agonism)
double Inotropy  = E_NE_b1 + E_EPI_b1 + E_DOBU_b1 + E_MIL + E_LEVO + E_OR;

// Combined chronotropy
double Chrono    = 0.6*(E_NE_b1 + E_EPI_b1 + E_DOBU_b1) + 0.2*E_MIL;

// Combined SVR effects (alpha = ↑ ; beta2/PDE3/Levo = ↓)
double Vaso_up   = E_NE_a1 + E_EPI_a1 + E_DA_a1 + E_VAS;
double Vaso_down = 0.4*E_EPI_b2 + 0.3*E_DOBU_b2 + 0.4*E_MIL + 0.3*E_LEVO + 0.2*E_DA_d1;

// ---- Inflammation / pathophysiology dynamics ----
double damp_drv  = (LAC > 2.0 ? (LAC-2.0)*0.1 : 0.0) + STUN*0.05;
dxdt_DAMP   = k_DAMP_synth*(1.0 + damp_drv) - k_DAMP_deg*DAMP;
dxdt_TNFa   = k_TNF_synth*DAMP - k_TNF_deg*TNFa;
dxdt_NOex   = k_NO_synth*(1.0 + k_NO_indEX*TNFa) - k_NO_deg*NOex;

// Sympathetic and RAAS: driven by low MAP and perfusion
double MAP_dev = (75.0 - MAP)/75.0;
dxdt_SNS    = 0.05*fmax(0.0, MAP_dev) - 0.04*(SNS - 1.0);
dxdt_RAAS   = 0.04*fmax(0.0, MAP_dev) - 0.04*(RAAS - 1.0);

// ---- Myocardial recovery / damage ----
// Stunning recovers slowly if perfusion adequate
double perf_adeq = (MAP >= 65 ? 1.0 : MAP/65.0);
dxdt_STUN   = -k_stun_recov*STUN*perf_adeq +
              k_isch_dmg*(LAC>4 ? 1.0 : 0.0) -
              0.005*(E_NE_b1 + E_DOBU_b1)*STUN +    // adrenergic stunning recovery
              0.008*Impella_unload*Impella_ON*STUN  // mechanical unloading aids recovery
              -0.0;
// EF dynamics: improved by inotropy, hurt by stunning & NO myocardial depression
double EF_target = LVEF0 + 0.10*Inotropy - 0.15*(NOex - 0.3) - 0.05*STUN + 0.04*Impella_ON;
if (EF_target < 0.10) EF_target = 0.10;
if (EF_target > 0.65) EF_target = 0.65;
dxdt_EF = 0.05*(EF_target - EF);

// HR dynamics
double HR_target = HR0 + 25.0*Chrono - 5.0*E_VAS - 5.0*E_LEVO;
dxdt_HR = 0.10*(HR_target - HR);

// SVR dynamics
double SVR_target = 1200.0 + 1200.0*Vaso_up - 600.0*Vaso_down - 800.0*(NOex - 0.3);
if (SVR_target < 400) SVR_target = 400;
if (SVR_target > 3000) SVR_target = 3000;
dxdt_SVR = k_SVR_synth*(SVR_target - SVR) - k_SVR_deg*0;   // first-order toward target

// CO and MAP
double CO_intrinsic = (EF*100.0)*HR / 1000.0 * (BSA/1.8);   // simplified L/min
double CO_MCS = IABP_ON*IABP_CO + Impella_ON*Impella_CO + ECMO_ON*ECMO_CO;
double CO_target = CO_intrinsic + CO_MCS;
if (CO_target < 1.0) CO_target = 1.0;
dxdt_COx = 0.20*(CO_target - COx);

double MAP_target = COx * SVR / 80.0;     // simplified Ohm relation
dxdt_MAP = 0.20*(MAP_target - MAP);

// PCWP — improved by inotropy, vasodilation, MCS unload
double PCWP_target = PCWP0 - 8.0*E_MIL - 6.0*E_LEVO - 4.0*Impella_ON*Impella_unload + 2.0*ECMO_ON*ECMO_unload;
if (PCWP_target < 6) PCWP_target = 6;
dxdt_PCWP = 0.10*(PCWP_target - PCWP);

// SvO2 and lactate (oxygen delivery vs demand)
double DO2 = COx * 13.4 * 0.95;   // assume Hb=10 g/dL, SaO2 0.95
double VO2_baseline = 250.0;       // mL/min
double O2_extr = VO2_baseline / fmax(DO2, 1e-3);
double SvO2_target = 0.95 - O2_extr;
if (SvO2_target < 0.20) SvO2_target = 0.20;
if (SvO2_target > 0.85) SvO2_target = 0.85;
dxdt_SVO2 = 0.20*(SvO2_target - SVO2);

double lac_prod = (DO2 < 600 ? (600 - DO2)*0.005 : 0.0);
dxdt_LAC = lac_prod - k_lac_clear*LAC + 0.05*Epi_or_isch_demand_marker(EPI_inf);

// Multi-organ damage proxies
double aki_drv = (MAP < 65 ? (65-MAP)*0.03 : 0.0) + (LAC>3 ? (LAC-3)*0.04 : 0.0);
dxdt_CR = k_Cr_rise*aki_drv - k_Cr_clear*(CR - CR0);
dxdt_ALT = k_ALT_rise*aki_drv*0.5 - k_ALT_clear*(ALT - ALT0);
double UO_target = (MAP >= 70 ? UO0 : UO0*MAP/70.0);
dxdt_UO = 0.2*(UO_target - UO);

// Mortality hazard accumulation
double inst_h = h_base + b_lac*fmax(0.0, LAC - 2.0) +
                b_MAP*fmax(0.0, (65.0 - MAP)/65.0) +
                b_CO*fmax(0.0, (2.2 - COx)/2.2) +
                b_NOex*fmax(0.0, NOex - 0.6);
dxdt_CUMHAZ = inst_h;

$TABLE
double SaO2 = 0.95;
double DO2_calc = COx * 13.4 * SaO2;
double SurvP = exp(-CUMHAZ);
double SCAI = 0;
if (LAC > 8.0 || MAP < 50.0 || COx < 1.2) SCAI = 5;   // E
else if (LAC > 5.0 || MAP < 60.0) SCAI = 4;            // D
else if (LAC > 3.0 || MAP < 65.0) SCAI = 3;            // C
else if (LAC > 2.0)              SCAI = 2;             // B
else                             SCAI = 1;             // A

$CAPTURE EF COx MAP SVR PCWP HR LAC SVO2 CR ALT UO NOex TNFa CUMHAZ SurvP SCAI DO2_calc Inotropy Vaso_up Vaso_down

$GLOBAL
// Helper inline; keep mathematical interpretation simple
inline double Epi_or_isch_demand_marker(double epi_inf) {
  return 0.02 * epi_inf;   // epinephrine drives lactate (well-known clinical observation)
}
'

# ---- Compile -----------------------------------------------------------------
# Note: To run, you may comment out the build at file load and call build_cs()
build_cs <- function() {
  mread("cs_qsp_model", code = cs_code)
}

# =============================================================================
# Therapeutic Scenarios -------------------------------------------------------
# 1. Untreated baseline (no support, no drugs) → severe SCAI E
# 2. Norepinephrine + Dobutamine (standard guideline-directed)
# 3. Norepinephrine + Milrinone (DOREMI alternative)
# 4. Levosimendan rescue (SURVIVE-style)
# 5. NE + Dobutamine + IABP (IABP-SHOCK II analogue)
# 6. NE + Impella CP (DanGer Shock arm)
# 7. NE + VA-ECMO (ECLS-SHOCK ECMO arm)
# =============================================================================
run_cs_scenarios <- function(end = 72, delta = 0.1) {
  mod <- build_cs()
  ev_baseline <- ev(amt = 0, ID = 1) %>% mutate(scen = "Baseline")
  scen_list <- list(
    Baseline = list(),
    NE_Dobu  = list(NE_inf = 0.10, DOBU_inf = 5),
    NE_Mil   = list(NE_inf = 0.10, MIL_inf  = 0.25),
    Levosim  = list(NE_inf = 0.05, LEVO_inf = 0.10),
    NE_IABP  = list(NE_inf = 0.10, DOBU_inf = 5, IABP_ON = 1),
    NE_Impella = list(NE_inf = 0.10, Impella_ON = 1),
    NE_ECMO  = list(NE_inf = 0.10, ECMO_ON = 1)
  )
  out_all <- bind_rows(lapply(names(scen_list), function(nm) {
    pars <- scen_list[[nm]]
    sim <- mod %>% param(pars) %>% mrgsim(end = end, delta = delta) %>% as_tibble()
    sim$Scenario <- nm
    sim
  }))
  out_all
}

# Helper: quick plot
plot_cs <- function(sim) {
  ggplot(sim, aes(x = time)) +
    geom_line(aes(y = MAP, colour = Scenario), linewidth = 0.6) +
    facet_wrap(~"MAP (mmHg)") + theme_bw()
}

# Run when sourced interactively:
#   sim <- run_cs_scenarios(); plot_cs(sim)
#   ggplot(sim, aes(x=time, y=LAC, colour=Scenario)) + geom_line()
#   ggplot(sim, aes(x=time, y=COx,  colour=Scenario)) + geom_line()
#   ggplot(sim, aes(x=time, y=SurvP, colour=Scenario)) + geom_line()
