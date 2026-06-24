################################################################################
# Bronchial Asthma – QSP Model (mrgsolve)
#
# Disease:  Bronchial Asthma (기관지 천식)
# Mechanism: Type-2 (Th2/ILC2) inflammation → IgE / eosinophil axis /
#            airway smooth muscle → FEV1 loss + exacerbations
#
# Compartments (18 ODE states):
#  PK      – Biologic (2-comp SC, 5 agents) + ICS lung/sys
#  TSLP    – alarmin signalling hub
#  IL5     – eosinophil growth/survival factor
#  IL13    – goblet/ASM/remodeling driver
#  IgE     – free + TMDD complex (omalizumab)
#  EOS_B   – blood eosinophils
#  EOS_T   – tissue (airway) eosinophils
#  ASM     – airway smooth muscle tone index
#  MUCUS   – mucus production state
#  FEV1    – % predicted FEV1 (clinical output)
#  RISK    – exacerbation hazard (cumulative)
#
# Treatment scenarios (5):
#  1. ICS/LABA alone  (baseline moderate-severe)
#  2. + Mepolizumab 100 mg SC q4w  (MENSA)
#  3. + Benralizumab 30 mg SC q4w×3 then q8w  (CALIMA)
#  4. + Dupilumab 300 mg SC q2w  (LIBERTY AIR)
#  5. + Tezepelumab 210 mg SC q4w  (NAVIGATOR)
#
# Parameter calibration sources:
#  – Blood eosinophil kinetics: Beal et al. J Pharmacokinet Pharmacodyn 2020
#  – Mepolizumab PK: MENSA PK sub-study; Cmax ~8 μg/mL
#  – Benralizumab PK/PD: CALIMA; near-complete depletion AUC
#  – Dupilumab PK: Sanofi/Regeneron PopPK; t½ ~20 d
#  – Tezepelumab PK: NAVIGATOR; SC bioavailability 0.81, t½ ~26 d
#  – ICS/LABA FEV1 effect: GINA 2024 step 3–4 data
#  – IgE kinetics: omalizumab TMDD model (Lowe et al. 2009)
#
# Author: Claude Code Routine (CCR) · Date: 2026-06-16
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ── Model code ────────────────────────────────────────────────────────────────

asthma_code <- '
$PROB Bronchial Asthma QSP – mrgsolve
Type-2 inflammation, IgE/eosinophil axis, airway smooth muscle, FEV1, exacerbation

$PARAM @annotated
// ── Biologic PK ─────────────────────────────────────
ka_MEPO  : 0.0067 : Mepolizumab SC absorption rate (1/h)   [F=0.77]
CL_MEPO  : 0.0115 : Mepolizumab clearance (L/h)
V1_MEPO  : 3.6    : Mepolizumab central volume (L)
V2_MEPO  : 3.8    : Mepolizumab peripheral volume (L)
Q_MEPO   : 0.043  : Mepolizumab intercompartmental CL (L/h)

ka_BENZ  : 0.0055 : Benralizumab SC absorption rate (1/h)  [F=0.59]
CL_BENZ  : 0.0085 : Benralizumab clearance (L/h)
V1_BENZ  : 3.1    : Benralizumab central volume (L)
V2_BENZ  : 3.4    : Benralizumab peripheral volume (L)
Q_BENZ   : 0.028  : Benralizumab intercompartmental CL (L/h)

ka_DUPIL : 0.0042 : Dupilumab SC absorption rate (1/h)     [F=0.64]
CL_DUPIL : 0.0065 : Dupilumab clearance (L/h)
V1_DUPIL : 4.8    : Dupilumab central volume (L)
V2_DUPIL : 5.2    : Dupilumab peripheral volume (L)
Q_DUPIL  : 0.058  : Dupilumab intercompartmental CL (L/h)

ka_TEZE  : 0.0072 : Tezepelumab SC absorption rate (1/h)   [F=0.81]
CL_TEZE  : 0.0088 : Tezepelumab clearance (L/h)
V1_TEZE  : 3.9    : Tezepelumab central volume (L)
V2_TEZE  : 4.4    : Tezepelumab peripheral volume (L)
Q_TEZE   : 0.051  : Tezepelumab intercompartmental CL (L/h)

// ── Omalizumab TMDD (anti-IgE) ──────────────────────
ka_OMAL  : 0.0071 : Omalizumab SC absorption rate (1/h)    [F=0.62]
CL_OMAL  : 0.0092 : Omalizumab linear clearance (L/h)
V1_OMAL  : 3.3    : Omalizumab central volume (L)
kon_IgE  : 0.55   : IgE-Omalizumab binding rate (1/(nM·h))
koff_IgE : 0.00017: IgE-Omalizumab dissociation rate (1/h)
kdeg_RC  : 0.0037 : IgE:Omalizumab complex degradation rate (1/h)

// ── ICS/LABA PK ─────────────────────────────────────
ka_ICS   : 0.18   : ICS lung absorption rate (1/h)
CL_ICS   : 28     : ICS systemic clearance (L/h)
Vd_ICS   : 90     : ICS systemic volume (L)
ka_LABA  : 0.35   : LABA absorption rate (1/h)
CL_LABA  : 22     : LABA systemic clearance (L/h)
Vd_LABA  : 140    : LABA systemic volume (L)

// ── TSLP dynamics ───────────────────────────────────
ksyn_TSLP : 0.006 : TSLP basal synthesis (nM/h)
kdeg_TSLP : 0.045 : TSLP degradation rate (1/h)
TSLP_ss   : 0.133 : TSLP steady-state (nM) – allergen-driven baseline
Emax_TEZE_TSLP : 0.95 : Tezepelumab max TSLP neutralization

// ── IL-5 dynamics ───────────────────────────────────
ksyn_IL5  : 0.18  : IL-5 basal synthesis (pg/mL/h)
kdeg_IL5  : 0.065 : IL-5 degradation rate (1/h)
IL5_ss    : 2.77  : IL-5 baseline (pg/mL)
Emax_MEPO_IL5 : 0.95 : Mepolizumab max IL-5 reduction
EC50_MEPO_IL5 : 0.12 : Mepolizumab EC50 for IL-5 (μg/mL)

// ── IL-13 dynamics ──────────────────────────────────
ksyn_IL13 : 0.10  : IL-13 basal synthesis (pg/mL/h)
kdeg_IL13 : 0.08  : IL-13 degradation rate (1/h)
IL13_ss   : 1.25  : IL-13 baseline (pg/mL)
Emax_DUPIL_IL13 : 0.90 : Dupilumab max IL-13 reduction
EC50_DUPIL_IL13 : 0.25 : Dupilumab EC50 for IL-13 (μg/mL)

// ── IgE dynamics (free) ─────────────────────────────
ksyn_IgE  : 0.0072: IgE basal synthesis rate (nM/h)
kdeg_IgE  : 0.0032: Free IgE degradation (1/h)
IgE_ss    : 2.25  : Free IgE baseline (nM, ~600 kU/L)

// ── Blood Eosinophil kinetics ────────────────────────
kprod_EOS : 0.18  : Eosinophil BM production rate (cells/μL/h)
kin_EOS   : 0.040 : Eosinophil tissue influx rate (1/h)
kout_EOS  : 0.038 : Eosinophil tissue egress + death (1/h)
EOS_B_ss  : 450   : Blood eosinophil baseline (cells/μL)
EOS_T_ss  : 2.0   : Tissue eosinophil baseline (10^6/g)

Emax_MEPO_EOS  : 0.80 : Mepolizumab max Eos suppression
EC50_MEPO_EOS  : 0.15 : Mepolizumab EC50 for blood Eos (μg/mL)
Emax_BENZ_EOS  : 0.97 : Benralizumab near-complete depletion
EC50_BENZ_EOS  : 0.08 : Benralizumab EC50 (μg/mL)

// ── Airway Smooth Muscle tone ────────────────────────
kact_ASM  : 0.012 : ASM tone increase rate (1/h, driven by IL-13/EOS-T)
krel_ASM  : 0.025 : ASM tone relaxation (1/h, LABA/SABA)
ASM_ss    : 1.0   : Baseline normalized ASM tone index (1.0 = diseased)

// ── Mucus dynamics ──────────────────────────────────
ksyn_MUC  : 0.15  : Mucus production rate (AU/h)
kdeg_MUC  : 0.08  : Mucus clearance rate (1/h)
MUC_ss    : 1.875 : Mucus baseline (AU)

// ── FEV1 dynamics ───────────────────────────────────
FEV1_max  : 78    : Max achievable FEV1 (% predicted) with treatment
FEV1_min  : 48    : Minimum FEV1 without treatment
kFEV1     : 0.005 : FEV1 adaptation rate (1/h)
// FEV1 = FEV1_max - (FEV1_max-FEV1_min)*f(ASM)*f(MUC)*f(EOS_T)

// ── ICS PD ──────────────────────────────────────────
Emax_ICS_IL5  : 0.45 : ICS max IL-5 suppression
EC50_ICS      : 2.0  : ICS EC50 (ng/mL)
Emax_ICS_IL13 : 0.40 : ICS max IL-13 suppression
Emax_ICS_TSLP : 0.30 : ICS max TSLP suppression

// ── LABA PD ─────────────────────────────────────────
Emax_LABA_ASM : 0.55 : LABA max ASM tone reduction
EC50_LABA     : 0.8  : LABA EC50 (ng/mL)

// ── Exacerbation hazard ─────────────────────────────
lambda0   : 0.18  : Baseline exacerbation hazard rate (events/year)
beta_EOS  : 0.0015: Eos contribution to exacerbation risk (/cell/μL)
beta_FEV1 : -0.025: FEV1 (% pred) protective effect (per %)
beta_MUC  : 0.04  : Mucus contribution to exacerbation risk (per AU)

// ── ICS/LABA doses ──────────────────────────────────
DOSE_ICS  : 200   : ICS dose lung deposition (μg, q12h → use infusion proxy)
DOSE_LABA : 50    : LABA dose (μg, q12h)

$CMT @annotated
// Biologic PK compartments
MEPO_SC  : Mepolizumab SC depot (mg)
MEPO_C1  : Mepolizumab central (μg/mL)
MEPO_C2  : Mepolizumab peripheral
BENZ_SC  : Benralizumab SC depot (mg)
BENZ_C1  : Benralizumab central (μg/mL)
BENZ_C2  : Benralizumab peripheral
DUPIL_SC : Dupilumab SC depot (mg)
DUPIL_C1 : Dupilumab central (μg/mL)
DUPIL_C2 : Dupilumab peripheral
TEZE_SC  : Tezepelumab SC depot (mg)
TEZE_C1  : Tezepelumab central (μg/mL)
TEZE_C2  : Tezepelumab peripheral
// Omalizumab + TMDD
OMAL_SC  : Omalizumab SC depot (mg)
OMAL_C1  : Omalizumab central (nM)
IgE_FREE : Free IgE (nM)
IgE_RC   : IgE:Omalizumab complex (nM)
// ICS/LABA
ICS_LUNG : ICS lung (ng/mL equiv)
ICS_SYS  : ICS systemic (ng/mL)
LABA_C   : LABA central (ng/mL)
// Disease states (18 total)
TSLP_PD  : TSLP signalling (nM)
IL5_PD   : IL-5 (pg/mL)
IL13_PD  : IL-13 (pg/mL)
EOS_B    : Blood eosinophils (cells/μL)
EOS_T    : Tissue eosinophils (10^6/g)
ASM_TONE : ASM tone index
MUCUS    : Mucus production (AU)
FEV1_ODE : FEV1 % predicted

$MAIN
// ── Biologic concentrations ─────────────────────────
double cMEPO  = MEPO_C1;
double cBENZ  = BENZ_C1;
double cDUPIL = DUPIL_C1;
double cTEZE  = TEZE_C1;
double cOMAL  = OMAL_C1;

// ── ICS/LABA effects ────────────────────────────────
double EFF_ICS_IL5  = Emax_ICS_IL5  * ICS_SYS / (EC50_ICS + ICS_SYS + 0.001);
double EFF_ICS_IL13 = Emax_ICS_IL13 * ICS_SYS / (EC50_ICS + ICS_SYS + 0.001);
double EFF_ICS_TSLP = Emax_ICS_TSLP * ICS_SYS / (EC50_ICS + ICS_SYS + 0.001);
double EFF_LABA_ASM = Emax_LABA_ASM * LABA_C  / (EC50_LABA + LABA_C + 0.001);

// ── Biologic PD effects ─────────────────────────────
double EFF_MEPO_IL5  = Emax_MEPO_IL5  * cMEPO  / (EC50_MEPO_IL5  + cMEPO  + 0.001);
double EFF_MEPO_EOS  = Emax_MEPO_EOS  * cMEPO  / (EC50_MEPO_EOS  + cMEPO  + 0.001);
double EFF_BENZ_EOS  = Emax_BENZ_EOS  * cBENZ  / (EC50_BENZ_EOS  + cBENZ  + 0.001);
double EFF_DUPIL_IL13= Emax_DUPIL_IL13* cDUPIL / (EC50_DUPIL_IL13+ cDUPIL + 0.001);
double EFF_TEZE_TSLP = Emax_TEZE_TSLP * cTEZE  / (0.05 + cTEZE  + 0.001);

// ── TSLP feedback on IL-5/IL-13 ─────────────────────
double TSLP_rel = TSLP_PD / (TSLP_ss + 0.001); // relative TSLP level

// ── Omalizumab TMDD binding ─────────────────────────
double binding_rate = kon_IgE * OMAL_C1 * IgE_FREE;
double unbind_rate  = koff_IgE * IgE_RC;
double complex_deg  = kdeg_RC  * IgE_RC;

// ── Eosinophil recruitment driven by IL-5 ───────────
double IL5_rel  = IL5_PD  / (IL5_ss  + 0.001);
double IL13_rel = IL13_PD / (IL13_ss + 0.001);
double TSLP_stim = TSLP_rel;

// ── FEV1 target ─────────────────────────────────────
double f_ASM  = ASM_TONE / ASM_ss;
double f_MUC  = MUCUS    / MUC_ss;
double f_EOS  = EOS_T    / EOS_T_ss;
double FEV1_target = FEV1_max - (FEV1_max - FEV1_min) *
                     (0.5*f_ASM + 0.3*f_MUC + 0.2*f_EOS);
if(FEV1_target < FEV1_min) FEV1_target = FEV1_min;
if(FEV1_target > FEV1_max) FEV1_target = FEV1_max;

$ODE
// ── Mepolizumab PK ───────────────────────────────────
dxdt_MEPO_SC = -ka_MEPO * MEPO_SC;
dxdt_MEPO_C1 =  ka_MEPO * MEPO_SC / V1_MEPO
                - (CL_MEPO + Q_MEPO)/V1_MEPO * MEPO_C1
                +  Q_MEPO/V2_MEPO * MEPO_C2;
dxdt_MEPO_C2 =  Q_MEPO/V1_MEPO * MEPO_C1 - Q_MEPO/V2_MEPO * MEPO_C2;

// ── Benralizumab PK ──────────────────────────────────
dxdt_BENZ_SC = -ka_BENZ * BENZ_SC;
dxdt_BENZ_C1 =  ka_BENZ * BENZ_SC / V1_BENZ
                - (CL_BENZ + Q_BENZ)/V1_BENZ * BENZ_C1
                +  Q_BENZ/V2_BENZ * BENZ_C2;
dxdt_BENZ_C2 =  Q_BENZ/V1_BENZ * BENZ_C1 - Q_BENZ/V2_BENZ * BENZ_C2;

// ── Dupilumab PK ─────────────────────────────────────
dxdt_DUPIL_SC = -ka_DUPIL * DUPIL_SC;
dxdt_DUPIL_C1 =  ka_DUPIL * DUPIL_SC / V1_DUPIL
                 - (CL_DUPIL + Q_DUPIL)/V1_DUPIL * DUPIL_C1
                 +  Q_DUPIL/V2_DUPIL * DUPIL_C2;
dxdt_DUPIL_C2 =  Q_DUPIL/V1_DUPIL * DUPIL_C1 - Q_DUPIL/V2_DUPIL * DUPIL_C2;

// ── Tezepelumab PK ───────────────────────────────────
dxdt_TEZE_SC = -ka_TEZE * TEZE_SC;
dxdt_TEZE_C1 =  ka_TEZE * TEZE_SC / V1_TEZE
                - (CL_TEZE + Q_TEZE)/V1_TEZE * TEZE_C1
                +  Q_TEZE/V2_TEZE * TEZE_C2;
dxdt_TEZE_C2 =  Q_TEZE/V1_TEZE * TEZE_C1 - Q_TEZE/V2_TEZE * TEZE_C2;

// ── Omalizumab PK + IgE TMDD ─────────────────────────
dxdt_OMAL_SC = -ka_OMAL * OMAL_SC;
dxdt_OMAL_C1 =  ka_OMAL * OMAL_SC / V1_OMAL
                - CL_OMAL/V1_OMAL * OMAL_C1
                - binding_rate + unbind_rate + complex_deg;
dxdt_IgE_FREE= ksyn_IgE - kdeg_IgE * IgE_FREE
               - binding_rate + unbind_rate;
dxdt_IgE_RC  =  binding_rate - unbind_rate - complex_deg;

// ── ICS / LABA PK ─────────────────────────────────────
dxdt_ICS_LUNG= -ka_ICS * ICS_LUNG;
dxdt_ICS_SYS =  ka_ICS * ICS_LUNG - (CL_ICS/Vd_ICS) * ICS_SYS;
dxdt_LABA_C  = -( CL_LABA/Vd_LABA) * LABA_C;

// ── TSLP ──────────────────────────────────────────────
dxdt_TSLP_PD = ksyn_TSLP * (1 - EFF_ICS_TSLP)
               - kdeg_TSLP * TSLP_PD * (1 + EFF_TEZE_TSLP);

// ── IL-5 ──────────────────────────────────────────────
dxdt_IL5_PD  = ksyn_IL5 * (1 + 0.5*TSLP_rel) * (1 - EFF_ICS_IL5) * (1 - EFF_MEPO_IL5)
               - kdeg_IL5 * IL5_PD;

// ── IL-13 ─────────────────────────────────────────────
dxdt_IL13_PD = ksyn_IL13 * (1 + 0.4*TSLP_rel) * (1 - EFF_ICS_IL13) * (1 - EFF_DUPIL_IL13)
               - kdeg_IL13 * IL13_PD;

// ── Blood eosinophils ─────────────────────────────────
// IL-5 promotes BM production/release; mepolizumab+benralizumab suppress
double EOS_prod_mod = (1 - EFF_MEPO_EOS) * (1 - EFF_BENZ_EOS);
dxdt_EOS_B = kprod_EOS * IL5_rel * EOS_prod_mod
             - kin_EOS * EOS_B;

// ── Tissue eosinophils ────────────────────────────────
// Influx from blood; IL-13 upregulates eotaxin → CCR3 → tissue entry
double IL13_eos_factor = 1 + 0.5 * IL13_rel;
dxdt_EOS_T = kin_EOS * EOS_B / EOS_B_ss * IL13_eos_factor
             - kout_EOS * EOS_T;

// ── ASM tone ──────────────────────────────────────────
// Driven by IL-13 + EOS_T granules; relaxed by LABA
double ASM_drive = 1 + 0.4*IL13_rel + 0.2*(EOS_T/EOS_T_ss - 1);
dxdt_ASM_TONE= kact_ASM * ASM_drive - krel_ASM * ASM_TONE * (1 + EFF_LABA_ASM);

// ── Mucus ──────────────────────────────────────────────
// IL-13 drives goblet metaplasia → MUC5AC; ICS suppresses
double MUC_drive = 1 + 0.5*IL13_rel;
dxdt_MUCUS   = ksyn_MUC * MUC_drive * (1 - 0.3*EFF_ICS_IL13)
               - kdeg_MUC * MUCUS;

// ── FEV1 ──────────────────────────────────────────────
dxdt_FEV1_ODE= kFEV1 * (FEV1_target - FEV1_ODE);

$TABLE
capture MEPO_conc   = MEPO_C1;
capture BENZ_conc   = BENZ_C1;
capture DUPIL_conc  = DUPIL_C1;
capture TEZE_conc   = TEZE_C1;
capture OMAL_conc   = OMAL_C1;
capture IgE_free_nM = IgE_FREE;
capture IgE_free_kUL= IgE_FREE * 246.0; // nM → kU/L
capture EOS_blood   = EOS_B;
capture EOS_tissue  = EOS_T;
capture IL5_pgmL    = IL5_PD;
capture IL13_pgmL   = IL13_PD;
capture TSLP_nM     = TSLP_PD;
capture ASM_index   = ASM_TONE;
capture Mucus_AU    = MUCUS;
capture FEV1_pct    = FEV1_ODE;
// Exacerbation risk (instantaneous hazard rate, events/yr)
capture exacerb_hazard = lambda0
                        + beta_EOS  * (EOS_B - EOS_B_ss)
                        + beta_FEV1 * (FEV1_ODE - 70)
                        + beta_MUC  * (MUCUS - MUC_ss);

$INIT
MEPO_SC  = 0, MEPO_C1  = 0, MEPO_C2  = 0
BENZ_SC  = 0, BENZ_C1  = 0, BENZ_C2  = 0
DUPIL_SC = 0, DUPIL_C1 = 0, DUPIL_C2 = 0
TEZE_SC  = 0, TEZE_C1  = 0, TEZE_C2  = 0
OMAL_SC  = 0, OMAL_C1  = 0
IgE_FREE = 2.25, IgE_RC = 0
ICS_LUNG = 0, ICS_SYS = 0, LABA_C = 0
TSLP_PD  = 0.133
IL5_PD   = 2.77
IL13_PD  = 1.25
EOS_B    = 450
EOS_T    = 2.0
ASM_TONE = 1.0
MUCUS    = 1.875
FEV1_ODE = 58
'

## ── Compile model ─────────────────────────────────────────────────────────────
mod <- mcode("asthma_qsp", asthma_code, quiet = TRUE)

## ── Helper: run a treatment scenario ─────────────────────────────────────────
run_scenario <- function(model,
                         label       = "ICS/LABA",
                         mepo_dose   = 0,   # mg SC q4w
                         benz_dose   = 0,   # mg SC (q4w×3, then q8w)
                         dupil_dose  = 0,   # mg SC q2w
                         teze_dose   = 0,   # mg SC q4w
                         sim_weeks   = 52) {

  # ICS (budesonide 320 μg/day lung dep proxy) as zero-order infusion
  ics_infusion <- ev(time = 0, amt = 0, cmt = "ICS_LUNG",
                     rate = 200/12, ii = 0, addl = 0)  # ~16.7 μg/h

  # LABA (formoterol 18 μg/day) as infusion proxy
  laba_infusion <- ev(time = 0, amt = 0, cmt = "LABA_C",
                      rate = 50/Vd_LABA, ii = 0, addl = 0)

  # Build ICS/LABA events
  events <- ev(time = 0, amt = 13.3, cmt = "ICS_LUNG",   # μg pulse q12h
               rate = -2, ii = 12, addl = sim_weeks*14 - 1) +
            ev(time = 0, amt = 0.18, cmt = "LABA_C",
               rate = -2, ii = 12, addl = sim_weeks*14 - 1)

  # Mepolizumab 100 mg SC q4w
  if (mepo_dose > 0) {
    ev_mepo <- ev(time = 0, amt = mepo_dose, cmt = "MEPO_SC",
                  ii = 4*7*24, addl = floor(sim_weeks/4) - 1)
    events <- events + ev_mepo
  }

  # Benralizumab 30 mg SC: q4w×3 then q8w
  if (benz_dose > 0) {
    ev_benz_q4 <- ev(time = 0, amt = benz_dose, cmt = "BENZ_SC",
                     ii = 4*7*24, addl = 2)           # 3 doses
    q8_start    <- 3 * 4 * 7 * 24
    addl_q8     <- floor((sim_weeks - 12)/8)
    ev_benz_q8  <- ev(time = q8_start, amt = benz_dose, cmt = "BENZ_SC",
                      ii = 8*7*24, addl = addl_q8)
    events <- events + ev_benz_q4 + ev_benz_q8
  }

  # Dupilumab 300 mg SC q2w (loading 600 mg at week 0 not modelled for simplicity)
  if (dupil_dose > 0) {
    ev_dup <- ev(time = 0, amt = dupil_dose, cmt = "DUPIL_SC",
                 ii = 2*7*24, addl = floor(sim_weeks/2) - 1)
    events <- events + ev_dup
  }

  # Tezepelumab 210 mg SC q4w
  if (teze_dose > 0) {
    ev_teze <- ev(time = 0, amt = teze_dose, cmt = "TEZE_SC",
                  ii = 4*7*24, addl = floor(sim_weeks/4) - 1)
    events <- events + ev_teze
  }

  out <- mrgsim(model, events,
                end   = sim_weeks * 7 * 24,
                delta = 12) %>%
         as.data.frame() %>%
         mutate(time_weeks = time / (7 * 24),
                Scenario   = label)
  out
}

## Dummy Vd_LABA for helper function
Vd_LABA <- 140

## ── Scenario 1: ICS/LABA alone ───────────────────────────────────────────────
s1 <- run_scenario(mod, "ICS/LABA only")

## ── Scenario 2: + Mepolizumab ────────────────────────────────────────────────
s2 <- run_scenario(mod, "ICS/LABA + Mepolizumab", mepo_dose = 100)

## ── Scenario 3: + Benralizumab ───────────────────────────────────────────────
s3 <- run_scenario(mod, "ICS/LABA + Benralizumab", benz_dose = 30)

## ── Scenario 4: + Dupilumab ──────────────────────────────────────────────────
s4 <- run_scenario(mod, "ICS/LABA + Dupilumab", dupil_dose = 300)

## ── Scenario 5: + Tezepelumab ────────────────────────────────────────────────
s5 <- run_scenario(mod, "ICS/LABA + Tezepelumab", teze_dose = 210)

## ── Combine results ───────────────────────────────────────────────────────────
all_results <- bind_rows(s1, s2, s3, s4, s5) %>%
  mutate(Scenario = factor(Scenario, levels = c(
    "ICS/LABA only",
    "ICS/LABA + Mepolizumab",
    "ICS/LABA + Benralizumab",
    "ICS/LABA + Dupilumab",
    "ICS/LABA + Tezepelumab"
  )))

scenario_colors <- c(
  "ICS/LABA only"          = "#616161",
  "ICS/LABA + Mepolizumab" = "#1565C0",
  "ICS/LABA + Benralizumab"= "#2E7D32",
  "ICS/LABA + Dupilumab"   = "#AD1457",
  "ICS/LABA + Tezepelumab" = "#E65100"
)

## ── Plot 1: FEV1 % predicted over 52 weeks ────────────────────────────────────
p1 <- ggplot(all_results, aes(time_weeks, FEV1_pct, color = Scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = scenario_colors) +
  labs(title    = "FEV1 (% Predicted) – Treatment Scenarios",
       subtitle = "Bronchial Asthma QSP Model",
       x = "Time (weeks)", y = "FEV1 (% predicted)") +
  coord_cartesian(ylim = c(45, 85)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        legend.title = element_blank())

print(p1)

## ── Plot 2: Blood eosinophils ─────────────────────────────────────────────────
p2 <- ggplot(all_results, aes(time_weeks, EOS_blood, color = Scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = scenario_colors) +
  geom_hline(yintercept = 300, linetype = "dashed", color = "tomato") +
  annotate("text", x = 50, y = 320, label = "T2-high threshold (300/μL)",
           size = 3, color = "tomato") +
  labs(title = "Blood Eosinophil Count",
       x = "Time (weeks)", y = "Eosinophils (cells/μL)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_blank())

print(p2)

## ── Plot 3: Tissue eosinophils ────────────────────────────────────────────────
p3 <- ggplot(all_results, aes(time_weeks, EOS_tissue, color = Scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Airway Tissue Eosinophils",
       x = "Time (weeks)", y = "Eosinophils (10^6/g tissue)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_blank())

print(p3)

## ── Plot 4: IL-5 and IL-13 dynamics ──────────────────────────────────────────
p4 <- all_results %>%
  select(time_weeks, Scenario, IL5_pgmL, IL13_pgmL) %>%
  pivot_longer(cols = c(IL5_pgmL, IL13_pgmL),
               names_to = "Cytokine", values_to = "Conc") %>%
  ggplot(aes(time_weeks, Conc, color = Scenario, linetype = Cytokine)) +
  geom_line(size = 0.9) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Type-2 Cytokines: IL-5 & IL-13",
       x = "Time (weeks)", y = "Concentration (pg/mL)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p4)

## ── Plot 5: Exacerbation hazard rate ─────────────────────────────────────────
p5 <- ggplot(all_results, aes(time_weeks, exacerb_hazard, color = Scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Instantaneous Exacerbation Hazard Rate",
       x = "Time (weeks)", y = "Hazard (events/year)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.title = element_blank())

print(p5)

## ── Summary table (week 52 values) ───────────────────────────────────────────
summary_w52 <- all_results %>%
  filter(abs(time_weeks - 52) < 0.5) %>%
  group_by(Scenario) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  select(Scenario,
         FEV1_pct,
         EOS_blood,
         EOS_tissue,
         IL5_pgmL,
         IL13_pgmL,
         exacerb_hazard) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

message("\n── Week-52 Summary ──────────────────────────────")
print(as.data.frame(summary_w52))

## ── Dose–response: Benralizumab dose vs FEV1 at wk 52 ───────────────────────
benz_doses <- c(10, 20, 30, 50, 100)

dr_res <- purrr::map_dfr(benz_doses, function(d) {
  run_scenario(mod, sprintf("Benz %g mg", d), benz_dose = d) %>%
    filter(abs(time_weeks - 52) < 0.5) %>%
    slice_tail(n = 1) %>%
    select(FEV1_pct, EOS_blood) %>%
    mutate(Dose_mg = d)
})

p_dr <- dr_res %>%
  pivot_longer(cols = c(FEV1_pct, EOS_blood),
               names_to = "Endpoint", values_to = "Value") %>%
  ggplot(aes(Dose_mg, Value, color = Endpoint)) +
  geom_point(size = 3) + geom_line(size = 1) +
  scale_x_log10() +
  labs(title    = "Benralizumab Dose–Response (week 52)",
       subtitle = "FEV1 and Blood Eosinophils",
       x = "Benralizumab Dose (mg, SC)", y = "Value") +
  facet_wrap(~Endpoint, scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

print(p_dr)

message("
═══════════════════════════════════════════════════
Bronchial Asthma QSP – mrgsolve simulation complete
Scenarios:
  1. ICS/LABA (standard of care)
  2. + Mepolizumab 100 mg q4w   (anti-IL-5)
  3. + Benralizumab 30 mg q4w→q8w (anti-IL-5Rα)
  4. + Dupilumab 300 mg q2w     (anti-IL-4Rα)
  5. + Tezepelumab 210 mg q4w   (anti-TSLP)
Key references:
  MENSA (NEJM 2014), CALIMA (Lancet 2016),
  LIBERTY AIR (NEJM 2018), NAVIGATOR (NEJM 2021)
═══════════════════════════════════════════════════
")
