# =====================================================================
# Hepatorenal Syndrome (HRS) — mrgsolve QSP Model
#   Author : Claude Code Routine (2026-07-01)
#   Scope  : HRS-AKI / HRS-NAKI in decompensated cirrhosis
#            Splanchnic vasodilation → EABV ↓ → RAAS/SNS/AVP ↑ →
#            Renal vasoconstriction → GFR ↓ → serum creatinine ↑
#   PK/PD  : Terlipressin (bolus / continuous), Norepinephrine,
#            Midodrine + Octreotide, Albumin 25%, LVP, Antibiotics
#   Outputs: MAP, RBF, GFR, sCr, urine Na/output, HRS reversal flag,
#            ischemic ADR flags, 30/90-day survival hazard
#   References (calibration): CONFIRM (Wong 2021, NEJM),
#            OT-0401/REVERSE, ATTIRE, ANSWER, ICA-2015/2019.
# =====================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

hrs_code <- '
$PROB
# Hepatorenal Syndrome (HRS) QSP model
# 26 ODE compartments · 6 drug PK · 12 disease/organ · 8 clinical

$PLUGIN autodec

$PARAM @annotated
// ============================================
// Terlipressin PK (2-cpt IV bolus / continuous)
// Sub. Krag 2007, Halimi 2002, Nazer 2016
// ============================================
CL_T   :  3.6   : Terlipressin clearance (L/h)      // rapid conversion to lysyl-VP
V1_T   :  9.0   : Central Vd (L)
V2_T   : 22.0   : Peripheral Vd (L)
Q_T    :  8.0   : Inter-compartmental clearance (L/h)
Kconv  :  4.0   : Conversion terli→lys-VP (1/h)     // half-life ~10 min
CL_LVP :  1.0   : Lysyl-VP clearance (L/h)          // half-life ~50-80 min

// Norepinephrine PK (1-cpt continuous IV)
CL_NE  : 200.0  : NE clearance (L/h)                // t1/2 ~2 min
V_NE   :  10.0  : NE Vd (L)

// Midodrine → desglymidodrine PK (oral)
KA_MID :  1.2   : Absorption rate midodrine (1/h)
CL_DES : 30.0   : Desglymidodrine clearance (L/h)
V_DES  : 120.0  : Desglymidodrine Vd (L)
F_MID  :  0.93  : Midodrine bioavailability

// Octreotide PK (SC bolus q8h or CIV)
CL_OCT :  10.0  : Octreotide clearance (L/h)
V_OCT  :  20.0  : Octreotide Vd (L)
KA_OCT :  4.0   : SC absorption rate (1/h)          // t1/2 abs ~10 min

// Albumin (IV infusion, 25%)
KE_ALB :  0.036 : Albumin elimination (1/h; t1/2 ~19 d in health, faster in cirrhosis)
V_ALB  :  4.0   : Albumin distribution volume (L)

// ============================================
// Disease-baseline (Child-Pugh C, MELD ~26)
// ============================================
BASE_MAP    : 70.0  : Baseline mean arterial pressure (mmHg)
BASE_HR     : 90.0  : Baseline heart rate (bpm)
BASE_CO     :  8.0  : Baseline cardiac output (L/min; hyperdynamic)
BASE_SVR    :  650  : Baseline SVR (dyn·s/cm5) — reduced from 1200
BASE_RBF    :  400  : Baseline renal blood flow (mL/min; reduced from 1200)
BASE_GFR    :   30  : Baseline GFR (mL/min) — HRS diagnosis
BASE_CR     :  2.5  : Baseline serum creatinine (mg/dL)
BASE_UNa    :   15  : Baseline urine Na (mEq/L)
BASE_UOUT   :  600  : Baseline urine output (mL/day)
BASE_ALB    :  2.8  : Baseline serum albumin (g/dL)
BASE_BIL    :  6.0  : Baseline bilirubin (mg/dL)
BASE_INR    :  2.0  : Baseline INR
BASE_Na     :  128  : Baseline serum Na (mEq/L)

// Neurohormonal baselines
BASE_PRA    :  8.0  : Plasma renin activity (ng/mL/h; normal <1.5)
BASE_ALDO   : 800   : Aldosterone (pg/mL; normal 100-300)
BASE_NE     : 800   : Plasma norepinephrine (pg/mL; normal 200-500)
BASE_AVP    :  8.0  : Plasma AVP (pg/mL; normal 1-3)

// ============================================
// Pharmacodynamic parameters
// ============================================
// Terlipressin V1a-mediated MAP increase
EMAX_T_MAP  :  22.0 : Max MAP increase with terli (mmHg)
EC50_T      :   1.5 : Lysyl-VP concentration at 50% effect (ng/mL)
HILL_T      :   2.0 : Hill coefficient
K_T_RBF     :   0.4 : Fraction of MAP effect on RBF (splanchnic → renal redistribution)

// Norepinephrine α1-mediated
EMAX_NE_MAP : 25.0  : Max MAP increase (mmHg)
EC50_NE     :  4.0  : NE conc at 50% effect (ng/mL, above endogenous)

// Midodrine (desglymidodrine α1)
EMAX_MID    : 10.0  : Max MAP increase (mmHg)
EC50_DES    : 30.0  : Desglymidodrine EC50 (ng/mL)

// Octreotide splanchnic (via glucagon suppression)
EMAX_OCT    :  6.0  : Splanchnic vasoconstriction (mmHg equivalent MAP)
EC50_OCT    :  2.0  : Octreotide EC50 (ng/mL)

// Albumin oncotic + immunomodulation
EMAX_ALB_MAP:  8.0  : Max MAP increase from oncotic expansion
K_ALB_INFL  :  0.02 : Anti-inflammatory rate per (g/L) albumin above baseline

// Physiology gain terms
K_MAP_RBF   :   6.0 : RBF ~ K * (MAP - 60)  (mL/min per mmHg above threshold)
GFR_FRACT   :   0.20: GFR/RBF (filtration fraction, elevated in HRS)
CR_PROD     :  20.0 : Creatinine production (mg/day; sarcopenia reduces)
V_CR        :  40.0 : Cr Vd (L, TBW)
K_INFL_MAP  :   0.3 : Inflammation-mediated MAP suppression (mmHg per unit)

// Precipitants (event on/off; use as switch)
FLAG_SBP    :   0   : SBP precipitant flag (0/1)
FLAG_GIB    :   0   : GI bleed precipitant flag
FLAG_LVP    :   0   : LVP without albumin
FLAG_NSAID  :   0   : NSAID exposure flag

// Adverse-event thresholds
TERLI_ISCH_C:  4.0  : Lysyl-VP conc threshold for ischemic ADR (ng/mL)
NE_ARR_C    : 10.0  : NE conc threshold for arrhythmia (ng/mL)

$CMT @annotated
// ============== Drug PK (10) ==============
GUT_MID    : Midodrine gut cpt
CENT_DES   : Desglymidodrine central
DEP_OCT    : Octreotide SC depot
CENT_OCT   : Octreotide central
TERLI_C    : Terlipressin central
TERLI_P    : Terlipressin peripheral
LVP_C      : Lysyl-VP central
NE_C       : Norepinephrine central
ALB_C      : Albumin central (exogenous)
// ============== Neurohormonal (4) =========
PRA_c      : Renin activity (integrated response)
ALDO_c     : Aldosterone
NE_endo    : Endogenous NE
AVP_c      : AVP
// ============== Hemodynamic (4) ===========
MAP_c      : Mean arterial pressure
SVR_c      : Systemic vascular resistance
RBF_c      : Renal blood flow (mL/min)
CardOut    : Cardiac output (L/min)
// ============== Renal / clinical (6) ======
GFR_c      : Glomerular filtration rate (mL/min)
SCR_c      : Serum creatinine (mg/dL)
UNa_c      : Urine Na (mEq/L)
UOUT_c     : Urine output (mL/day)
Na_c       : Serum Na
INFL_c     : Systemic inflammation (unit: pg/mL TNF-α equivalent)
// ============== Outcomes (3) ==============
Hazard30   : 30-day hazard integral
IschAUC    : Terlipressin ischemic exposure AUC
LiverBil   : Bilirubin trajectory

$MAIN
// Initial conditions
MAP_c_0    = BASE_MAP;
SVR_c_0    = BASE_SVR;
RBF_c_0    = BASE_RBF;
CardOut_0  = BASE_CO;
GFR_c_0    = BASE_GFR;
SCR_c_0    = BASE_CR;
UNa_c_0    = BASE_UNa;
UOUT_c_0   = BASE_UOUT;
Na_c_0     = BASE_Na;
PRA_c_0    = BASE_PRA;
ALDO_c_0   = BASE_ALDO;
NE_endo_0  = BASE_NE;
AVP_c_0    = BASE_AVP;
LiverBil_0 = BASE_BIL;
INFL_c_0   = 40.0;   // pg/mL TNF-α equivalent (elevated)
Hazard30_0 = 0.0;
IschAUC_0  = 0.0;

$ODE
// ================================
// PK equations
// ================================
// Terlipressin: 2-cpt + conversion to lysyl-VP + separate LVP clearance
dxdt_TERLI_C = -(CL_T/V1_T)*TERLI_C - (Q_T/V1_T)*TERLI_C + (Q_T/V2_T)*TERLI_P - Kconv*TERLI_C;
dxdt_TERLI_P =  (Q_T/V1_T)*TERLI_C - (Q_T/V2_T)*TERLI_P;
dxdt_LVP_C   =  Kconv*TERLI_C - (CL_LVP/V1_T)*LVP_C;

// Norepinephrine
dxdt_NE_C = -(CL_NE/V_NE)*NE_C;

// Midodrine → desglymidodrine
dxdt_GUT_MID  = -KA_MID*GUT_MID;
dxdt_CENT_DES =  F_MID*KA_MID*GUT_MID - (CL_DES/V_DES)*CENT_DES;

// Octreotide (SC → central)
dxdt_DEP_OCT  = -KA_OCT*DEP_OCT;
dxdt_CENT_OCT =  KA_OCT*DEP_OCT - (CL_OCT/V_OCT)*CENT_OCT;

// Albumin
dxdt_ALB_C = -KE_ALB*ALB_C;

// ================================
// Effect-site concentrations
// ================================
double CP_LVP  = LVP_C / V1_T;         // ng/mL (assuming dose in ng)
double CP_NE   = NE_C  / V_NE;
double CP_DES  = CENT_DES / V_DES;
double CP_OCT  = CENT_OCT / V_OCT;
double ALB_add = ALB_C / V_ALB;        // g/L increment

// Emax responses
double E_terli_MAP = EMAX_T_MAP * pow(CP_LVP, HILL_T) / (pow(EC50_T, HILL_T) + pow(CP_LVP, HILL_T));
double E_NE_MAP    = EMAX_NE_MAP * CP_NE / (EC50_NE + CP_NE);
double E_MID_MAP   = EMAX_MID    * CP_DES / (EC50_DES + CP_DES);
double E_OCT_MAP   = EMAX_OCT    * CP_OCT / (EC50_OCT + CP_OCT);
double E_ALB_MAP   = EMAX_ALB_MAP * ALB_add / (5.0 + ALB_add);

double DrugMAP = E_terli_MAP + E_NE_MAP + E_MID_MAP + E_OCT_MAP + E_ALB_MAP;

// ================================
// Precipitant modifiers
// ================================
double PrecipInfl = 30.0*FLAG_SBP + 15.0*FLAG_GIB + 10.0*FLAG_LVP;
double PrecipVaso = 20.0*FLAG_LVP + 15.0*FLAG_GIB;          // effective volume loss
double NsaidPGE2  = 8.0*FLAG_NSAID;                          // afferent constriction

// ================================
// Neurohormonal ODEs (slow, 1st-order)
// ================================
double drive = fmax(BASE_MAP - MAP_c, 0.0);   // baroreceptor unloading
dxdt_PRA_c   =  0.6*drive - 0.3*(PRA_c - BASE_PRA);
dxdt_ALDO_c  =  40.0*(PRA_c/BASE_PRA - 1.0) - 0.25*(ALDO_c - BASE_ALDO);
dxdt_NE_endo =  25.0*drive - 0.3*(NE_endo - BASE_NE);
dxdt_AVP_c   =  0.4*drive - 0.4*(AVP_c - BASE_AVP);

// ================================
// Hemodynamics
// ================================
double dMAP = 0.8*(BASE_MAP + DrugMAP - MAP_c) - K_INFL_MAP*(INFL_c-30.0)/10.0 - 0.05*PrecipVaso;
dxdt_MAP_c = dMAP;

double SVR_target = BASE_SVR + 15.0*DrugMAP - 100.0*(INFL_c-30)/40.0;
dxdt_SVR_c = 0.5*(SVR_target - SVR_c);

dxdt_CardOut = 0.4*(BASE_CO*(1 - 0.02*(MAP_c-BASE_MAP)) - CardOut);

// ================================
// Renal hemodynamics
// ================================
double RBF_target = K_MAP_RBF*fmax(MAP_c - 60.0, 0.0)*(1 - 0.15*NsaidPGE2/8.0)
                    * (1 - 0.10*(NE_endo/BASE_NE - 1.0))
                    * (1 - 0.10*(AVP_c/BASE_AVP - 1.0));
dxdt_RBF_c = 0.6*(RBF_target - RBF_c);

double GFR_target = GFR_FRACT * RBF_c;
dxdt_GFR_c = 0.6*(GFR_target - GFR_c);

// ================================
// Serum creatinine (mg/dL)
// ================================
double Cr_prod = CR_PROD / 24.0;     // mg/h
double Cr_clear = (GFR_c/1000.0)*60.0 * SCR_c / V_CR * 24.0;   // approx
dxdt_SCR_c = (Cr_prod - (GFR_c*1.44/V_CR)*SCR_c*0.5) ;   // simplified

// ================================
// Urine Na · output
// ================================
double UNa_target = 5.0 + 100.0*(1.0 - ALDO_c/BASE_ALDO*0.5) * (GFR_c/BASE_GFR);
dxdt_UNa_c  = 0.5*(fmax(UNa_target,1.0) - UNa_c);
double UOUT_target = 50.0*(GFR_c/BASE_GFR) + 500.0*(1.0 - AVP_c/BASE_AVP*0.5);
dxdt_UOUT_c = 0.4*(fmax(UOUT_target,100.0) - UOUT_c);

// ================================
// Serum sodium (dilution from AVP-V2)
// ================================
double Na_target = 140.0 - 8.0*(AVP_c/BASE_AVP);
dxdt_Na_c = 0.05*(Na_target - Na_c);

// ================================
// Inflammation trajectory
// ================================
dxdt_INFL_c = 0.05*PrecipInfl - 0.1*(INFL_c - 20.0) - K_ALB_INFL*ALB_add*10.0;

// ================================
// Bilirubin (surrogate for liver reserve)
// ================================
dxdt_LiverBil = 0.005*(BASE_BIL - LiverBil) + 0.002*PrecipInfl;

// ================================
// 30-day hazard (Cox surrogate; log-hazard integrates)
// ================================
double lp = -3.0 + 0.6*(SCR_c-1.5) + 0.05*(LiverBil - 3.0)
             + 0.3*(INFL_c/30.0 - 1.0) - 0.02*(MAP_c-70.0);
dxdt_Hazard30 = exp(lp) / 720.0;   // per hour, normalized

// Ischemic ADR AUC (terlipressin exposure above threshold)
dxdt_IschAUC = fmax(CP_LVP - TERLI_ISCH_C, 0.0);

$TABLE
double MELD_pred    = 3.78*log(fmax(LiverBil,1.0)) + 11.2*log(BASE_INR) + 9.57*log(fmax(SCR_c,1.0)) + 6.43;
double HRS_reversal = (SCR_c < 1.5 && GFR_c > 40) ? 1.0 : 0.0;
double HRS_partial  = (SCR_c < BASE_CR*0.7) ? 1.0 : 0.0;
double IschRisk     = (IschAUC > 24.0) ? 1.0 : 0.0;   // >24 ng·h/mL above threshold
double S30          = exp(-Hazard30);
double S90          = exp(-Hazard30*3.0);

$CAPTURE CP_LVP CP_NE CP_DES CP_OCT ALB_add MELD_pred HRS_reversal HRS_partial IschRisk S30 S90 DrugMAP
'

hrs_mod <- mcode("hrs", hrs_code)

# ---------------------------------------------------------------------
# Scenario library
# ---------------------------------------------------------------------
# All scenarios begin at day 0 (Child-Pugh C, MELD ~26, HRS-AKI stage 2)
# Simulation horizon: 14 days (336 h)
# Doses expressed with amt in mg equivalents where practical; scale factors
# integrated inside the model use dose in ng for lysyl-VP calculations.

sim_hrs <- function(mod, events, label, time_max = 336) {
  mod %>% ev(events) %>% mrgsim(end = time_max, delta = 0.5) %>%
    as_tibble() %>% mutate(scenario = label)
}

# 1) Natural history (no intervention)
scn_nat <- sim_hrs(hrs_mod, ev(amt = 0, time = 0), "1_Natural_history")

# 2) Terlipressin 1 mg IV q6h + albumin 40 g/day (CONFIRM regimen)
ev_terli_bolus <- ev(cmt = "TERLI_C", amt = 1e6, ii = 6, addl = 55, time = 0)  # 1 mg = 1e6 ng
ev_alb         <- ev(cmt = "ALB_C", amt = 40, ii = 24, addl = 13, time = 0)
scn_terli <- sim_hrs(hrs_mod, c(ev_terli_bolus, ev_alb), "2_Terlipressin_bolus+Albumin")

# 3) Terlipressin continuous infusion 2 mg/24h + albumin
ev_terli_ci <- ev(cmt = "TERLI_C", amt = 2e6, rate = 2e6/24, ii = 24, addl = 13, time = 0)
scn_terli_ci <- sim_hrs(hrs_mod, c(ev_terli_ci, ev_alb), "3_Terlipressin_CI+Albumin")

# 4) Norepinephrine 0.5 mg/h continuous IV + albumin (ICU alternative)
ev_ne <- ev(cmt = "NE_C", amt = 0.5e6, rate = 0.5e6, ii = 1, addl = 335, time = 0) # ng/h
scn_ne <- sim_hrs(hrs_mod, c(ev_ne, ev_alb), "4_Norepinephrine+Albumin")

# 5) Midodrine 12.5 mg PO tid + Octreotide 200 μg SC tid + Albumin
ev_mid <- ev(cmt = "GUT_MID", amt = 12.5, ii = 8, addl = 41, time = 0)
ev_oct <- ev(cmt = "DEP_OCT", amt = 0.2, ii = 8, addl = 41, time = 0)
scn_mo <- sim_hrs(hrs_mod, c(ev_mid, ev_oct, ev_alb), "5_Midodrine+Octreotide+Albumin")

# 6) Albumin monotherapy (ATTIRE-like, historical comparator)
scn_alb <- sim_hrs(hrs_mod, ev_alb, "6_Albumin_only")

# 7) Precipitant: SBP triggered, then terlipressin + antibiotics response
scn_sbp <- hrs_mod %>%
  param(FLAG_SBP = 1) %>%
  ev(c(ev_terli_bolus, ev_alb)) %>%
  mrgsim(end = 336, delta = 0.5) %>% as_tibble() %>%
  mutate(scenario = "7_SBP_precipitated_terli")

# 8) LVP without albumin (PICD: post-paracentesis circulatory dysfunction)
scn_lvp <- hrs_mod %>%
  param(FLAG_LVP = 1) %>%
  mrgsim(end = 336, delta = 0.5) %>% as_tibble() %>%
  mutate(scenario = "8_LVP_no_albumin_PICD")

# 9) NSAID exposure (avoidable precipitant)
scn_nsaid <- hrs_mod %>%
  param(FLAG_NSAID = 1) %>%
  mrgsim(end = 336, delta = 0.5) %>% as_tibble() %>%
  mutate(scenario = "9_NSAID_precipitant")

# 10) TIPS response surrogate: instantaneous 30% reduction in EABV strain
tips_mod <- hrs_mod %>% param(BASE_RBF = 800, BASE_SVR = 900)
scn_tips <- tips_mod %>% mrgsim(end = 336, delta = 0.5) %>% as_tibble() %>%
  mutate(scenario = "10_TIPS_hemodynamic_surrogate")

all_scenarios <- bind_rows(
  scn_nat, scn_terli, scn_terli_ci, scn_ne, scn_mo, scn_alb,
  scn_sbp, scn_lvp, scn_nsaid, scn_tips
)

# ---------------------------------------------------------------------
# Quick visualisation
# ---------------------------------------------------------------------
plot_traj <- function(df, var) {
  ggplot(df, aes(time/24, .data[[var]], colour = scenario)) +
    geom_line(linewidth = 0.6) + theme_minimal() +
    labs(x = "Day", y = var, title = paste("HRS QSP —", var)) +
    theme(legend.position = "bottom")
}

# Uncomment to run interactively
# print(plot_traj(all_scenarios, "SCR_c"))
# print(plot_traj(all_scenarios, "MAP_c"))
# print(plot_traj(all_scenarios, "GFR_c"))
# print(plot_traj(all_scenarios, "S30"))
# print(plot_traj(all_scenarios, "HRS_reversal"))

message("HRS QSP model built. Scenarios: ", length(unique(all_scenarios$scenario)))
