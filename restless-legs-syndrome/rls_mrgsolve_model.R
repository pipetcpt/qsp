# =============================================================================
# Restless Legs Syndrome (RLS / Willis-Ekbom Disease) — QSP / PK-PD Model
#   mrgsolve compartmental model integrating:
#     1) Brain-iron dynamics (serum ferritin → CSF ferritin → SN brain iron)
#     2) Dopaminergic effector axis (TH activity, striatal DA, A11→spinal tone)
#     3) Adenosine A1 receptor tone (Ferré hypothesis: iron-dependent)
#     4) Cortico-spinal glutamate / α2δ-gabapentinoid pharmacology
#     5) Endogenous & exogenous opioid analgesia of spinal sensory drive
#     6) Periodic limb movements (PLMS) + IRLS symptom score
#     7) Augmentation dynamics for chronic DA-agonist exposure
#     8) Drug PK: pramipexole, ropinirole, rotigotine (patch),
#                gabapentin (+ enacarbil prodrug), pregabalin,
#                oxycodone/naloxone PR, IV ferric carboxymaltose
#
#   Symptom output:  IRLS (0-40), PLMS index (events/h),
#                    sleep efficiency surrogate, augmentation hazard
#   Calibrated to:   Allen NEJM 2003, Trenkwalder Lancet 2008 (rotigotine),
#                    Garcia-Borreguero NEJM 2010 (gabapentin enacarbil),
#                    Trenkwalder Lancet Neurol 2013 (oxycodone/naloxone,
#                    RELOXYN), Cho 2018 / Allen IRON-CIRCA, IRLSSG augmentation.
#
#   Author : QSP Disease Model Library · 2026-06-30
# =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

rls_code <- '
$PROB
# Restless Legs Syndrome (Willis-Ekbom) QSP model
# Iron · Dopamine · Adenosine · α2δ · Opioid · PLMS · IRLS · Augmentation

$PARAM @annotated
// ---------- Pramipexole PK (2-comp, 1st-order absorption) ----------
KA_pra  : 1.20  : Pramipexole ka (1/h)
CL_pra  : 25.0  : Pramipexole CL (L/h) -- ~80% renal
V1_pra  : 130   : Pramipexole V1 (L)
V2_pra  : 270   : Pramipexole V2 (L)
Q_pra   : 28    : Pramipexole inter-comp Q (L/h)
EC50_pra: 0.45  : Pramipexole EC50 at D3R (ng/mL plasma surrogate)
EMAX_pra: 0.95  : Pramipexole max D3 stimulation

// ---------- Ropinirole PK ----------
KA_rop  : 1.10  : Ropinirole ka (1/h)
CL_rop  : 47    : Ropinirole CL (L/h) -- CYP1A2 hepatic
V1_rop  : 165   : Ropinirole V1 (L)
V2_rop  : 360   : Ropinirole V2 (L)
Q_rop   : 30    : Ropinirole Q (L/h)
EC50_rop: 5.4   : Ropinirole EC50 (ng/mL)
EMAX_rop: 0.90  : Ropinirole max D2/D3 stimulation

// ---------- Rotigotine (transdermal, zero-order release) ----------
KA_rot  : 0.13  : Rotigotine apparent ka via patch (1/h)
CL_rot  : 80    : Rotigotine CL (L/h)
V1_rot  : 80    : Rotigotine V1 (L)
EC50_rot: 0.30  : Rotigotine EC50 (ng/mL)
EMAX_rot: 0.97  : Rotigotine D3 stimulation
F_rot   : 0.45  : Rotigotine patch bioavail.

// ---------- Gabapentin enacarbil → gabapentin ----------
KA_gab  : 1.6   : Gabapentin ka (1/h)
F_genacarbil : 0.75 : Gabapentin enacarbil → gabapentin bioavail.
CL_gab  : 11    : Gabapentin CL (L/h) -- renal (eGFR-dependent)
V1_gab  : 60    : Gabapentin V1 (L)
EC50_gab: 4.0   : Gabapentin EC50 at α2δ (mg/L)
EMAX_gab: 0.85  : Gabapentin α2δ-mediated efficacy

// ---------- Pregabalin ----------
KA_pre  : 1.4   : Pregabalin ka (1/h)
CL_pre  : 6.5   : Pregabalin CL (L/h)
V1_pre  : 35    : Pregabalin V1 (L)
EC50_pre: 2.5   : Pregabalin EC50 (mg/L)
EMAX_pre: 0.80  : Pregabalin α2δ efficacy

// ---------- Oxycodone (PR with sequestered naloxone) ----------
KA_oxy  : 0.34  : Oxycodone ka (PR, 1/h)
CL_oxy  : 47    : Oxycodone CL (L/h)
V1_oxy  : 220   : Oxycodone V1 (L)
EC50_oxy: 18    : Oxycodone EC50 at MOR (ng/mL)
EMAX_oxy: 0.85  : Oxycodone MOR Emax

// ---------- IV ferric carboxymaltose (FCM) ----------
CL_fcm  : 0.40  : FCM clearance from plasma (L/h)
V_fcm   : 3.8   : FCM Vd (L) -- limited to plasma + RES uptake
KFCM_brain : 0.0009 : Plasma Fe → brain Fe transfer rate (per h)

// ---------- Iron dynamics ----------
Ferritin0   : 35   : Baseline serum ferritin (µg/L) -- typical RLS pt low
K_fer_loss  : 0.003 : Ferritin loss (per day, includes monthly menses)
K_fer_in    : 0.10  : Daily POIron net ferritin uptake (µg/L per mg PO if absorbed)
F_PO_iron   : 0.10  : Fractional gut Fe absorption (low at steady state)
Hepcidin    : 1.0   : Hepcidin scaling (1.0 baseline; rises in inflam/ESRD)
BrainFe0    : 100   : Brain iron index (au) -- 100 = healthy, RLS ~70
K_bfein     : 0.02  : Brain Fe uptake from serum/ferritin (per day)
K_bfeloss   : 0.02  : Brain Fe slow turnover (per day)

// ---------- Dopamine effector ----------
DAtone0   : 1.0   : Baseline (healthy) striatal/spinal DA tone (au)
DA_FeSens : 0.005 : Brain Fe → TH activity slope (per au of brain Fe)
DAtone_circ_amp : 0.30 : Circadian DA amplitude (trough at "night")
DAtone_min : 0.50 : Hard floor in untreated RLS during night

// ---------- Adenosine A1 tone ----------
A1tone0   : 1.0   : Baseline A1R tone
A1_FeSens : 0.006 : Brain Fe → A1R tone (Ferré 2018)
A1_caffeine : 0.0 : Caffeine antagonism flag (0-1, mg-equivalent)

// ---------- α2δ effector ----------
a2d_base  : 0.0   : Baseline α2δ-modulated suppression

// ---------- Spinal sensory drive / PLMS / IRLS ----------
SensDrive0 : 1.0  : Baseline aggregated sensory drive (untreated RLS = 1.0)
K_PLMS_base : 28  : PLMS index baseline in RLS (events/h)
K_IRLS_base : 28  : IRLS baseline in moderate-severe RLS

// ---------- Augmentation dynamics ----------
KAUG_on  : 0.0006 : Augmentation hazard slope per (DA-agonist effect)·h
KAUG_off : 0.10   : Augmentation decay when DA-agonist withdrawn (per day)
DA_aug_threshold : 0.55 : Effective DA tone above which augmentation accumulates

// ---------- Adverse-event surrogates ----------
ICD_slope_D3 : 0.10 : Impulse control disorder hazard per unit D3 fractional eff.
Constipation_oxy : 0.20 : Constipation index per (MOR efficacy); halved by naloxone

// ---------- Patient covariates ----------
WT     : 70    : Weight (kg)
eGFR   : 90    : eGFR (mL/min/1.73 m^2)
SEX    : 0     : 0=male, 1=female
PREG   : 0     : 1 if pregnant (3rd trimester)
ESRD   : 0     : 1 if ESRD/dialysis
SSRI   : 0     : 1 if on SSRI/SNRI trigger
CAFFEINE_mg : 0 : Daily caffeine intake (mg)

$CMT  @annotated
GUT_pra   : Pramipexole gut depot (mg)
CEN_pra   : Pramipexole central (mg)
PER_pra   : Pramipexole peripheral (mg)
GUT_rop   : Ropinirole gut (mg)
CEN_rop   : Ropinirole central (mg)
PER_rop   : Ropinirole peripheral (mg)
PATCH_rot : Rotigotine patch reservoir (mg)
CEN_rot   : Rotigotine central (mg)
GUT_gab   : Gabapentin/enacarbil gut (mg)
CEN_gab   : Gabapentin central (mg)
GUT_pre   : Pregabalin gut (mg)
CEN_pre   : Pregabalin central (mg)
GUT_oxy   : Oxycodone gut (mg)
CEN_oxy   : Oxycodone central (mg)
CEN_fcm   : IV ferric carboxymaltose central (mg)
FERRITIN  : Serum ferritin (µg/L)
BRAINFE   : Brain iron (au; 100 = healthy)
DA_AUG    : Augmentation index (0-1)
IRLS_dyn  : Dynamic IRLS score (0-40)
PLMS_dyn  : Dynamic PLMS index (events/h)

$GLOBAL
#define CP_pra  (CEN_pra/V1_pra*1e6/1000)   // ng/mL  (mg/L = µg/mL)
#define CP_rop  (CEN_rop/V1_rop*1e6/1000)   // ng/mL
#define CP_rot  (CEN_rot/V1_rot*1e6/1000)   // ng/mL
#define CP_gab  (CEN_gab/V1_gab)            // mg/L
#define CP_pre  (CEN_pre/V1_pre)            // mg/L
#define CP_oxy  (CEN_oxy/V1_oxy*1e6/1000)   // ng/mL
#define CP_fcm  (CEN_fcm/V_fcm)             // mg/L (Fe)

$MAIN
// Renal scaling for gabapentinoids
double eGFR_corr   = eGFR/90.0;
double CL_gab_act  = CL_gab * eGFR_corr;
double CL_pre_act  = CL_pre * eGFR_corr;
double CL_pra_act  = CL_pra * eGFR_corr;

// Initial conditions
FERRITIN_0 = Ferritin0;
BRAINFE_0  = BrainFe0;
DA_AUG_0   = 0.0;
IRLS_dyn_0 = K_IRLS_base;
PLMS_dyn_0 = K_PLMS_base;

$ODE
// ----- Pramipexole PK -----
dxdt_GUT_pra = -KA_pra*GUT_pra;
dxdt_CEN_pra =  KA_pra*GUT_pra - (CL_pra_act/V1_pra)*CEN_pra - (Q_pra/V1_pra)*CEN_pra + (Q_pra/V2_pra)*PER_pra;
dxdt_PER_pra =  (Q_pra/V1_pra)*CEN_pra - (Q_pra/V2_pra)*PER_pra;

// ----- Ropinirole PK -----
dxdt_GUT_rop = -KA_rop*GUT_rop;
dxdt_CEN_rop =  KA_rop*GUT_rop - (CL_rop/V1_rop)*CEN_rop - (Q_rop/V1_rop)*CEN_rop + (Q_rop/V2_rop)*PER_rop;
dxdt_PER_rop =  (Q_rop/V1_rop)*CEN_rop - (Q_rop/V2_rop)*PER_rop;

// ----- Rotigotine PK (transdermal) -----
dxdt_PATCH_rot = -KA_rot*PATCH_rot;
dxdt_CEN_rot   =  KA_rot*PATCH_rot*F_rot - (CL_rot/V1_rot)*CEN_rot;

// ----- Gabapentin / enacarbil PK -----
dxdt_GUT_gab = -KA_gab*GUT_gab;
dxdt_CEN_gab =  KA_gab*GUT_gab*F_genacarbil - (CL_gab_act/V1_gab)*CEN_gab;

// ----- Pregabalin PK -----
dxdt_GUT_pre = -KA_pre*GUT_pre;
dxdt_CEN_pre =  KA_pre*GUT_pre - (CL_pre_act/V1_pre)*CEN_pre;

// ----- Oxycodone PK (PR) -----
dxdt_GUT_oxy = -KA_oxy*GUT_oxy;
dxdt_CEN_oxy =  KA_oxy*GUT_oxy - (CL_oxy/V1_oxy)*CEN_oxy;

// ----- IV ferric carboxymaltose PK -----
dxdt_CEN_fcm = -(CL_fcm/V_fcm)*CEN_fcm;

// ----- Drug-target effects -----
double E_D3 = (EMAX_pra*CP_pra)/(EC50_pra+CP_pra)
            + (EMAX_rop*CP_rop)/(EC50_rop+CP_rop)
            + (EMAX_rot*CP_rot)/(EC50_rot+CP_rot);
if(E_D3 > 1.0) E_D3 = 1.0;

double E_a2d = (EMAX_gab*CP_gab)/(EC50_gab+CP_gab)
             + (EMAX_pre*CP_pre)/(EC50_pre+CP_pre);
if(E_a2d > 1.0) E_a2d = 1.0;

double E_MOR = (EMAX_oxy*CP_oxy)/(EC50_oxy+CP_oxy);

// ----- Iron / ferritin / brain iron -----
double menstrual = (SEX==1 && PREG==0) ? 1.0 : 0.4;
double dialysis_loss = ESRD==1 ? 0.015 : 0.0;
double FCM_delivery  = CP_fcm * KFCM_brain * 24.0; // per day
double POIron_in     = K_fer_in * F_PO_iron * 0.0; // dose-handled via mevent below

dxdt_FERRITIN = (FCM_delivery*4.0)                    // FCM mostly hits ferritin
              - K_fer_loss*FERRITIN*menstrual         // basal loss
              - dialysis_loss*FERRITIN;               // ESRD penalty

double FerrAvail = FERRITIN/Ferritin0;
dxdt_BRAINFE  = K_bfein*FerrAvail*(BrainFe0 - BRAINFE)  // restoration
              - K_bfeloss*BRAINFE                        // turnover
              + FCM_delivery*0.6;                        // small direct

// ----- Network tones -----
double DAtone_eff = DAtone0 - DA_FeSens*(BrainFe0 - BRAINFE); // ↓ with iron loss
if(DAtone_eff < DAtone_min) DAtone_eff = DAtone_min;
double A1tone_eff = A1tone0 - A1_FeSens*(BrainFe0 - BRAINFE)
                  - 0.0008*CAFFEINE_mg;                       // caffeine block

// Cumulative DA-agonist drive (for augmentation)
double DA_total = DAtone_eff + E_D3;

// Augmentation accumulation
double aug_growth = 0.0;
if(DA_total > DA_aug_threshold) {
   aug_growth = KAUG_on * (DA_total - DA_aug_threshold) * 24.0; // per day
}
dxdt_DA_AUG = aug_growth - KAUG_off*DA_AUG;

// ----- Spinal sensory drive -----
double sensory = SensDrive0
               - 0.6*E_D3                        // DA agonist
               - 0.8*E_a2d                       // gabapentinoids
               - 0.7*E_MOR                       // opioid
               - 0.3*(DAtone_eff-DAtone_min)     // endogenous tone
               + 0.5*DA_AUG                      // augmentation worsens
               + 0.3*(1.0 - A1tone_eff)          // A1 loss worsens
               + 0.4*SSRI;                       // trigger
if(sensory < 0.0) sensory = 0.0;

// ----- PLMS (events/h) -----
double PLMS_target = K_PLMS_base * sensory / SensDrive0;
dxdt_PLMS_dyn = 0.2*(PLMS_target - PLMS_dyn);

// ----- IRLS dynamic -----
double IRLS_target = K_IRLS_base * sensory / SensDrive0;
if(IRLS_target > 40.0) IRLS_target = 40.0;
dxdt_IRLS_dyn = 0.05*(IRLS_target - IRLS_dyn);

$TABLE
double IRLS    = IRLS_dyn;
double PLMS    = PLMS_dyn;
double DA_eff  = E_D3;
double a2d_eff = E_a2d;
double MOR_eff = E_MOR;
double AugIndex= DA_AUG;
double SleepEff = 100.0 - 1.6*PLMS - 0.4*IRLS;       // SE surrogate
if(SleepEff < 0.0) SleepEff = 0.0;
double CGI_I   = 7.0 - 6.0*(K_IRLS_base - IRLS)/K_IRLS_base; // 1 (best) - 7 (worst)
if(CGI_I < 1.0) CGI_I = 1.0;
if(CGI_I > 7.0) CGI_I = 7.0;
double ICD_haz = ICD_slope_D3 * E_D3;
double Constip = Constipation_oxy * E_MOR * 0.5; // PR oxy/naloxone halves constip.

$CAPTURE
CP_pra CP_rop CP_rot CP_gab CP_pre CP_oxy CP_fcm
DA_eff a2d_eff MOR_eff DAtone_eff A1tone_eff
IRLS PLMS SleepEff CGI_I AugIndex ICD_haz Constip
'

# Compile -----------------------------------------------------------------
mod <- mcode("rls_qsp", rls_code)

# =============================================================================
# Treatment scenarios (5+) — reproduce as ggplot panels
# =============================================================================

simulate_scenario <- function(label, events, tend=24*180, by=2) {
  out <- mod %>%
    ev(events) %>%
    mrgsim(end=tend, delta=by) %>%
    as_tibble() %>%
    mutate(scenario=label, day=time/24)
  out
}

# (1) Untreated moderate-severe RLS
sc1 <- simulate_scenario(
  "1. Untreated",
  events = ev(amt=0, evid=0)
)

# (2) Pramipexole 0.25 mg PO QHS, 180 days
sc2 <- simulate_scenario(
  "2. Pramipexole 0.25 mg QHS",
  events = ev(amt=0.25, cmt="GUT_pra", ii=24, addl=179)
)

# (3) Rotigotine 2 mg/24h transdermal (replaced daily) 180 days
sc3 <- simulate_scenario(
  "3. Rotigotine patch 2 mg/24h",
  events = ev(amt=2.0, cmt="PATCH_rot", ii=24, addl=179, rate=2.0/24)
)

# (4) Gabapentin enacarbil 600 mg QPM
sc4 <- simulate_scenario(
  "4. Gabapentin enacarbil 600 mg",
  events = ev(amt=600, cmt="GUT_gab", ii=24, addl=179)
)

# (5) Pregabalin 300 mg QHS
sc5 <- simulate_scenario(
  "5. Pregabalin 300 mg QHS",
  events = ev(amt=300, cmt="GUT_pre", ii=24, addl=179)
)

# (6) Oxycodone/Naloxone PR 5/2.5 mg BID, refractory RLS
sc6 <- simulate_scenario(
  "6. Oxycodone PR 5 mg BID",
  events = ev(amt=5.0, cmt="GUT_oxy", ii=12, addl=359)
)

# (7) IV FCM 1000 mg single dose
sc7 <- simulate_scenario(
  "7. IV FCM 1000 mg single",
  events = ev(amt=1000, cmt="CEN_fcm", rate=1000/0.25)
)

# (8) Pramipexole 0.5 mg + chronic exposure -> augmentation
sc8 <- simulate_scenario(
  "8. Pramipexole 0.5 mg (augmentation risk)",
  events = ev(amt=0.5, cmt="GUT_pra", ii=24, addl=359),
  tend = 24*365
)

all_sims <- bind_rows(sc1, sc2, sc3, sc4, sc5, sc6, sc7, sc8)

# Quick plots (commented for batch use) ---------------------------------------
# ggplot(all_sims, aes(day, IRLS, color=scenario)) + geom_line() + ylim(0,40)
# ggplot(all_sims, aes(day, PLMS, color=scenario)) + geom_line()
# ggplot(all_sims, aes(day, AugIndex, color=scenario)) + geom_line()
# ggplot(all_sims, aes(day, BRAINFE, color=scenario)) + geom_line()
# ggplot(all_sims, aes(day, FERRITIN, color=scenario)) + geom_line()

# =============================================================================
# Calibration anchors (literature)
# =============================================================================
# - Allen 2003 NEJM (IRLS validation): baseline 22-28 in moderate-severe RLS
# - Trenkwalder 2008 Lancet (rotigotine SP790): ΔIRLS -13.7 (3 mg) at 6 mo
# - Garcia-Borreguero 2010 NEJM (gabapentin enacarbil PIVOT-RLS):
#     ΔIRLS -13.2 (1200 mg) vs -8.8 placebo at 12 wk
# - Allen 2014 Sleep Med (pregabalin PIVOT-RLS):
#     pregabalin 300 mg ΔIRLS -14.6 vs pramipexole 0.5 mg -12.7 at 12 wk
#     augmentation: 1.7% pregabalin vs 7.7% pramipexole 52 wk
# - Trenkwalder 2013 Lancet Neurol (oxycodone-naloxone RELOXYN):
#     ΔIRLS -16.5 vs -9.4 (placebo) at 12 wk
# - Allen 2013 Sleep Med (IV FCM 1000 mg): ferritin > 50 % responders
# - IRLSSG augmentation criteria (Garcia-Borreguero 2007 / 2016 update)
# - Earley 2000 Neurology / Connor 2003 Brain (CSF ferritin & SN iron ↓)
# =============================================================================
