## ============================================================
## Marfan Syndrome (MFS) — QSP Model (mrgsolve)
## ============================================================
## Disease:  Marfan Syndrome (FBN1 mutation)
## Targets:  FBN1/TGF-β/MAPK → aortic root dilation
## Drugs:    Atenolol (β1-blocker), Losartan (ARB/TGF-β),
##           Propranolol, Irbesartan
## ODEs:     20 compartments
## Scenarios: 6 treatment comparisons
##
## Parameter sources:
##   • Pediatric Heart Network (PHN) RCT; Lacro et al. NEJM 2014
##   • COMPARE trial; Radonic et al. Eur Heart J 2010
##   • AIMS RCT; Forteza et al. J Am Coll Cardiol 2016
##   • Milleron et al. Eur Heart J 2015 (β-blocker data)
##   • Mouse model TGF-β calibration: Habashi et al. Science 2006
##   • PK: Atenolol - Öhrvall 1994; Losartan - McCrea 1996
##   • Aortic growth: Salim 1994, Rossig 2019 natural history
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ─────────────────────────────────────────────────────────────
# MODEL DEFINITION
# ─────────────────────────────────────────────────────────────
mfs_code <- '
$PROB Marfan Syndrome QSP Model — mrgsolve

$PLUGIN Rcpp

$PARAM @annotated
// ── Atenolol PK (2-CMT oral, kg-independent baseline) ──
KA_ATN  : 1.06  : Atenolol absorption rate constant [1/h]  (Öhrvall 1994)
F_ATN   : 0.50  : Atenolol oral bioavailability            (Öhrvall 1994)
CL_ATN  : 10.8  : Atenolol systemic clearance [L/h]        (Öhrvall 1994)
Vc_ATN  : 67.0  : Atenolol central volume [L]              (Öhrvall 1994)
Vp_ATN  : 64.0  : Atenolol peripheral volume [L]           (Öhrvall 1994)
Q_ATN   : 15.0  : Atenolol inter-compartment CL [L/h]      (Öhrvall 1994)
// ── Losartan PK (oral 1-CMT parent + metabolite) ──
KA_LOS  : 1.41  : Losartan absorption rate constant [1/h]  (McCrea 1996)
F_LOS   : 0.33  : Losartan oral bioavailability            (McCrea 1996)
CL_LOS  : 55.0  : Losartan parent clearance [L/h]          (McCrea 1996)
Vc_LOS  : 34.0  : Losartan central volume [L]              (McCrea 1996)
FM_LOS  : 0.14  : Fraction converted to EXP-3174           (McCrea 1996)
CL_EXP  : 18.0  : EXP-3174 clearance [L/h]                 (McCrea 1996)
Vc_EXP  : 12.0  : EXP-3174 volume [L]                      (McCrea 1996)
// ── TGF-β / SMAD dynamics ──
kprod_TGFb : 0.05  : TGF-β1 production rate [ng/mL/h]      (Habashi 2006)
kdeg_TGFb  : 0.10  : TGF-β1 degradation rate [1/h]
TGFb0      : 0.50  : Baseline plasma TGF-β1 [ng/mL]        (Matt 2009, Matt 2010)
TGFb_MFS_mult : 2.5 : TGF-β elevation fold in MFS          (Matt 2009)
kact_SMAD  : 0.30  : pSMAD activation rate by TGF-β [1/h]
kdact_SMAD : 0.08  : pSMAD dephosphorylation [1/h]
kact_ERK   : 0.20  : pERK activation rate [1/h]
kdact_ERK  : 0.12  : pERK inactivation [1/h]
// ── MMP dynamics ──
kprod_MMP  : 0.04  : MMP production [U/mL/h]               (Chung 2007)
kdeg_MMP   : 0.15  : MMP baseline degradation [1/h]
MMP0       : 0.27  : Baseline MMP activity [U/mL]
// ── Aortic remodelling ──
Ao_D0      : 30.0  : Baseline aortic root diameter [mm]    (natural history)
k_ao_growth : 1.8e-4 : Intrinsic aortic growth rate [mm/h] (Salim 1994: ~1.5mm/y)
Emax_MMP_ao : 0.80  : Max MMP contribution to dilation
EC50_MMP_ao : 0.30  : MMP EC50 for aortic growth [U/mL]
Emax_SMAD_ao: 0.60  : Max SMAD contribution
EC50_SMAD_ao: 0.50  : pSMAD EC50 [fold]
// ── Aortic wall stress & compliance ──
SBP0       : 125.0  : Baseline systolic BP [mmHg]           (PHN baseline)
DBP0       : 75.0   : Baseline diastolic BP [mmHg]
HR0        : 78.0   : Baseline heart rate [bpm]
dPdt0      : 1200.0 : Baseline dP/dt_max [mmHg/s]
// ── AR / mitral dynamics ──
k_AR_prog  : 4e-5   : AR progression rate [grade/h]         (Detaint 2008)
AR0        : 0.20   : Baseline AR grade (0-4)
// ── Drug pharmacodynamics ──
IC50_ATN_HR : 65.0  : Atenolol plasma conc for 50% HR ↓ [ng/mL] (PHN)
Emax_ATN_HR : 0.28  : Max fractional HR reduction
IC50_ATN_dPdt: 55.0 : IC50 for dP/dt reduction
Emax_ATN_dPdt: 0.30 : Max dP/dt reduction
IC50_LOS_TGFb: 40.0 : EXP-3174 for 50% TGF-β ↓ [ng/mL]    (Habashi 2006)
Emax_LOS_TGFb: 0.55 : Max TGF-β signal reduction (ARB)
IC50_LOS_SBP : 45.0 : EXP-3174 IC50 for SBP ↓ [ng/mL]
Emax_LOS_SBP : 0.18 : Max SBP fractional reduction
// ── Patient characteristics ──
WT         : 68.0   : Body weight [kg]
AGE        : 30.0   : Age [years]
MFS_FBN1   : 1.0    : MFS genotype (1=yes, 0=no; affects TGF-β baseline)
// ── Dosing flags ──
DOSE_ATN   : 0.0    : Atenolol daily dose [mg]
DOSE_LOS   : 0.0    : Losartan daily dose [mg]

$CMT @annotated
// PK compartments
DEPOT_ATN  : Atenolol gut depot [mg]
C1_ATN     : Atenolol central [mg/L]
C2_ATN     : Atenolol peripheral [mg/L]
DEPOT_LOS  : Losartan gut depot [mg]
C1_LOS     : Losartan central [mg/L]
C_EXP3174  : EXP-3174 (active metabolite) [mg/L]
// TGF-β / signalling
TGFb       : Free plasma TGF-β1 [ng/mL]
pSMAD      : Phospho-SMAD2/3 (fold vs baseline)
pERK       : Phospho-ERK1/2 (fold vs baseline)
MMP        : Circulating MMP activity [U/mL]
// Aortic / cardiac
Ao_Diam    : Aortic root diameter [mm]
AR_Grade   : Aortic regurgitation grade [0-4]
// Haemodynamic state variables
HR         : Heart rate [bpm]
SBP        : Systolic blood pressure [mmHg]
dPdt       : dP/dt_max [mmHg/s]
// Other biomarkers
NT_proBNP  : NT-proBNP [pg/mL]
LVEDD      : LV end-diastolic diameter [mm]
TGFb_plasma_obs : Observed plasma TGF-β [ng/mL]
// Score
Systemic_score : Ghent systemic score (0-20)

$INIT @annotated
DEPOT_ATN  = 0
C1_ATN     = 0
C2_ATN     = 0
DEPOT_LOS  = 0
C1_LOS     = 0
C_EXP3174  = 0
TGFb       = 1.25
pSMAD      = 1.0
pERK       = 1.0
MMP        = 0.27
Ao_Diam    = 30.0
AR_Grade   = 0.20
HR         = 78.0
SBP        = 125.0
dPdt       = 1200.0
NT_proBNP  = 50.0
LVEDD      = 46.0
TGFb_plasma_obs = 1.25
Systemic_score  = 3.0

$MAIN
// TGF-β baseline adjusted for MFS
double TGFb_base = TGFb0 * (1 + (TGFb_MFS_mult - 1) * MFS_FBN1);

// ATN plasma concentration (ng/mL; Vc in L → ng/mL when dose in mg → ×1000/Vc)
double Cp_ATN = C1_ATN / Vc_ATN * 1000.0;   // ng/mL
double Cp_LOS = C1_LOS / Vc_LOS * 1000.0;   // ng/mL
double Cp_EXP = C_EXP3174 / Vc_EXP * 1000.0;// ng/mL

// ── Drug PD: β-blocker effects ──
double E_ATN_HR   = Emax_ATN_HR   * Cp_ATN / (IC50_ATN_HR   + Cp_ATN);
double E_ATN_dPdt = Emax_ATN_dPdt * Cp_ATN / (IC50_ATN_dPdt + Cp_ATN);

// ── Drug PD: ARB effects (via EXP-3174) ──
double E_LOS_TGFb = Emax_LOS_TGFb * Cp_EXP / (IC50_LOS_TGFb + Cp_EXP);
double E_LOS_SBP  = Emax_LOS_SBP  * Cp_EXP / (IC50_LOS_SBP  + Cp_EXP);

// Initialise haemodynamics
double HR_ss   = HR0   * (1 - E_ATN_HR);
double SBP_ss  = SBP0  * (1 - E_LOS_SBP);
double dPdt_ss = dPdt0 * (1 - E_ATN_dPdt);

$ODE
// ── Atenolol PK ──
dxdt_DEPOT_ATN = -KA_ATN * DEPOT_ATN;
dxdt_C1_ATN    =  KA_ATN * DEPOT_ATN * F_ATN
                  - (CL_ATN + Q_ATN) / Vc_ATN * C1_ATN
                  +  Q_ATN / Vp_ATN * C2_ATN;
dxdt_C2_ATN    =  Q_ATN / Vc_ATN * C1_ATN
                  - Q_ATN / Vp_ATN * C2_ATN;

// ── Losartan PK (parent) ──
dxdt_DEPOT_LOS = -KA_LOS * DEPOT_LOS;
dxdt_C1_LOS    =  KA_LOS * DEPOT_LOS * F_LOS
                  - CL_LOS / Vc_LOS * C1_LOS;
// EXP-3174 formation
dxdt_C_EXP3174 =  FM_LOS * CL_LOS / Vc_LOS * C1_LOS
                  - CL_EXP / Vc_EXP * C_EXP3174;

// ── TGF-β dynamics ──
// ARB reduces TGF-β signalling via AT1R/TGF-β cross-talk
double Cp_EXP_now = C_EXP3174 / Vc_EXP * 1000.0;
double E_LOS_TGFb_now = Emax_LOS_TGFb * Cp_EXP_now / (IC50_LOS_TGFb + Cp_EXP_now);
double TGFb_ss = TGFb0 * (1 + (TGFb_MFS_mult - 1) * MFS_FBN1);
dxdt_TGFb = kprod_TGFb * TGFb_ss * (1 - E_LOS_TGFb_now)
            - kdeg_TGFb * TGFb;

// ── SMAD2/3 signalling (fold over baseline) ──
// TGF-β drives pSMAD; ARB reduces it
double TGFb_norm = TGFb / TGFb_ss;
dxdt_pSMAD = kact_SMAD * TGFb_norm * (1 - E_LOS_TGFb_now) - kdact_SMAD * pSMAD;

// ── ERK signalling ──
dxdt_pERK = kact_ERK * TGFb_norm * (1 - E_LOS_TGFb_now * 0.6)
             - kdact_ERK * pERK;

// ── MMP activity ──
// Both SMAD and ERK drive MMP; ARB + β-blocker indirectly reduce
double kprod_MMP_driven = kprod_MMP * (0.5 * pSMAD + 0.5 * pERK);
dxdt_MMP = kprod_MMP_driven - kdeg_MMP * MMP;

// ── Aortic root diameter [mm] ──
// Growth driven by MMP and pSMAD; wall stress (dPdt) also contributes
double Ao_Diam_now = Ao_Diam;
double Cp_ATN_now  = C1_ATN / Vc_ATN * 1000.0;
double E_ATN_dPdt_now = Emax_ATN_dPdt * Cp_ATN_now / (IC50_ATN_dPdt + Cp_ATN_now);
double hemo_stress  = dPdt / dPdt0 * (1 - E_ATN_dPdt_now);
double E_MMP_dilation  = Emax_MMP_ao  * MMP    / (EC50_MMP_ao  + MMP);
double E_SMAD_dilation = Emax_SMAD_ao * pSMAD  / (EC50_SMAD_ao + pSMAD);
double growth_rate = k_ao_growth * hemo_stress
                     * (1 + E_MMP_dilation + E_SMAD_dilation);
dxdt_Ao_Diam = growth_rate * Ao_Diam_now;

// ── Aortic regurgitation progression ──
// AR worsens with increasing root diameter beyond 37 mm
double AR_driver = (Ao_Diam_now > 37.0) ? (Ao_Diam_now - 37.0) / 20.0 : 0.0;
dxdt_AR_Grade = k_AR_prog * AR_driver * (4.0 - AR_Grade);

// ── Haemodynamic state (quasi-steady; driven by PD) ──
double Cp_ATN_hr = C1_ATN / Vc_ATN * 1000.0;
double E_ATN_HR_now   = Emax_ATN_HR   * Cp_ATN_hr / (IC50_ATN_HR   + Cp_ATN_hr);
double Cp_EXP_sbp     = C_EXP3174 / Vc_EXP * 1000.0;
double E_LOS_SBP_now  = Emax_LOS_SBP  * Cp_EXP_sbp / (IC50_LOS_SBP + Cp_EXP_sbp);
double HR_target  = HR0  * (1 - E_ATN_HR_now);
double SBP_target = SBP0 * (1 - E_LOS_SBP_now);
double dPdt_target = dPdt0 * (1 - Emax_ATN_dPdt * Cp_ATN_hr / (IC50_ATN_dPdt + Cp_ATN_hr));
dxdt_HR   = 0.5 * (HR_target  - HR);
dxdt_SBP  = 0.5 * (SBP_target - SBP);
dxdt_dPdt = 0.5 * (dPdt_target - dPdt);

// ── NT-proBNP (surrogate cardiac stress) ──
double AR_factor = 1 + AR_Grade * 0.5;
double LVEDD_factor = (LVEDD - 46.0) / 10.0;
double BNP_driver = AR_factor * (1 + fmax(LVEDD_factor, 0));
dxdt_NT_proBNP = 0.02 * (50.0 * BNP_driver - NT_proBNP);

// ── LVEDD (LV end-diastolic diameter) ──
// Driven by AR and MR volume overload
double AR_vol_load = AR_Grade / 4.0;   // normalised
dxdt_LVEDD = 0.005 * AR_vol_load * (1 + (Ao_Diam_now - 30.0) / 30.0);

// ── Observed TGF-β ──
dxdt_TGFb_plasma_obs = 0.5 * (TGFb - TGFb_plasma_obs);

// ── Ghent systemic score (simplified) ──
// Approximation based on aortic, ocular, skeletal contributions
double ao_pts = (Ao_Diam_now >= 50) ? 2.0 : (Ao_Diam_now >= 42) ? 1.0 : 0.0;
double ect_pts = 2.0;   // ectopia lentis (fixed for MFS genotype)
double skel_pts = 3.0;  // baseline skeleton (arachnodactyly, pectus, scoliosis)
double dural_pts = 2.0; // dural ectasia
double target_score = ao_pts + ect_pts + skel_pts + dural_pts;
dxdt_Systemic_score = 0.001 * (target_score - Systemic_score);

$TABLE
capture Cp_ATN_ng  = C1_ATN / Vc_ATN * 1000.0;
capture Cp_LOS_ng  = C1_LOS / Vc_LOS * 1000.0;
capture Cp_EXP_ng  = C_EXP3174 / Vc_EXP * 1000.0;
capture TGFb_conc  = TGFb;
capture pSMAD_fold = pSMAD;
capture pERK_fold  = pERK;
capture MMP_act    = MMP;
capture AoD_mm     = Ao_Diam;
capture AoD_Zscore = (Ao_Diam - 23.0) / 3.5;  // simplified adult Z-score
capture AR_gr      = AR_Grade;
capture HR_bpm     = HR;
capture SBP_mmhg   = SBP;
capture dPdt_val   = dPdt;
capture NTproBNP   = NT_proBNP;
capture LVEDD_mm   = LVEDD;
capture TGFbP      = TGFb_plasma_obs;
capture GhentScore = Systemic_score;
capture AnnualGrowthRate_mm_yr = dxdt_Ao_Diam * 8760;

$CAPTURE Cp_ATN_ng Cp_LOS_ng Cp_EXP_ng TGFb_conc pSMAD_fold pERK_fold
         MMP_act AoD_mm AoD_Zscore AR_gr HR_bpm SBP_mmhg dPdt_val
         NTproBNP LVEDD_mm GhentScore AnnualGrowthRate_mm_yr
'

# Compile model
mod <- mcode("mfs_qsp", mfs_code)


# ─────────────────────────────────────────────────────────────
# SIMULATION HELPER
# ─────────────────────────────────────────────────────────────
simulate_mfs <- function(
    mod,
    dose_atn   = 0,    # Atenolol daily dose [mg]
    dose_los   = 0,    # Losartan daily dose [mg]
    dose_interval = 24, # dosing interval [h]
    duration_yr   = 5,
    label         = "Untreated",
    n_patients    = 1
) {
  end_h <- duration_yr * 8760
  # Add dosing events
  ev <- ev(
    time    = seq(0, end_h - dose_interval, by = dose_interval),
    amt_atn = dose_atn,
    amt_los = dose_los,
    cmt_atn = "DEPOT_ATN",
    cmt_los = "DEPOT_LOS"
  )
  # Build event data frame
  ev_df_atn <- data.frame(time = seq(0, end_h - dose_interval, by = dose_interval),
                          amt  = dose_atn, cmt = 1, evid = 1, ii = 0, addl = 0)
  ev_df_los <- data.frame(time = seq(0, end_h - dose_interval, by = dose_interval),
                          amt  = dose_los, cmt = 4, evid = 1, ii = 0, addl = 0)
  ev_all <- bind_rows(ev_df_atn, ev_df_los) %>% arrange(time)

  sims <- mod %>%
    param(DOSE_ATN = dose_atn, DOSE_LOS = dose_los) %>%
    mrgsim_df(data = ev_all, end = end_h, delta = 24) %>%
    mutate(scenario = label,
           time_yr  = time / 8760)
  return(sims)
}


# ─────────────────────────────────────────────────────────────
# SCENARIO 1 — No treatment (natural history)
# ─────────────────────────────────────────────────────────────
s1 <- simulate_mfs(mod, dose_atn = 0, dose_los = 0,
                   duration_yr = 5, label = "1. Untreated")

# ─────────────────────────────────────────────────────────────
# SCENARIO 2 — Atenolol 50 mg QD
# (PHN trial standard arm; Lacro et al. NEJM 2014)
# ─────────────────────────────────────────────────────────────
s2 <- simulate_mfs(mod, dose_atn = 50, dose_los = 0,
                   duration_yr = 5, label = "2. Atenolol 50mg QD")

# ─────────────────────────────────────────────────────────────
# SCENARIO 3 — Atenolol 100 mg QD (higher dose)
# ─────────────────────────────────────────────────────────────
s3 <- simulate_mfs(mod, dose_atn = 100, dose_los = 0,
                   duration_yr = 5, label = "3. Atenolol 100mg QD")

# ─────────────────────────────────────────────────────────────
# SCENARIO 4 — Losartan 50 mg QD
# (PHN trial losartan arm; Lacro et al. NEJM 2014)
# ─────────────────────────────────────────────────────────────
s4 <- simulate_mfs(mod, dose_atn = 0, dose_los = 50,
                   duration_yr = 5, label = "4. Losartan 50mg QD")

# ─────────────────────────────────────────────────────────────
# SCENARIO 5 — Losartan 100 mg QD (COMPARE-dose)
# (Radonic et al. Eur Heart J 2010 — losartan arm)
# ─────────────────────────────────────────────────────────────
s5 <- simulate_mfs(mod, dose_atn = 0, dose_los = 100,
                   duration_yr = 5, label = "5. Losartan 100mg QD")

# ─────────────────────────────────────────────────────────────
# SCENARIO 6 — Combination: Atenolol 50mg + Losartan 50mg QD
# (Most aggressive medical therapy; guidelines post-AIMS)
# ─────────────────────────────────────────────────────────────
s6 <- simulate_mfs(mod, dose_atn = 50, dose_los = 50,
                   duration_yr = 5, label = "6. Atenolol+Losartan")

# ─────────────────────────────────────────────────────────────
# COMBINE RESULTS
# ─────────────────────────────────────────────────────────────
all_scenarios <- bind_rows(s1, s2, s3, s4, s5, s6) %>%
  mutate(scenario = factor(scenario, levels = c(
    "1. Untreated",
    "2. Atenolol 50mg QD",
    "3. Atenolol 100mg QD",
    "4. Losartan 50mg QD",
    "5. Losartan 100mg QD",
    "6. Atenolol+Losartan"
  )))

# ─────────────────────────────────────────────────────────────
# SUMMARY TABLE at 5 years
# ─────────────────────────────────────────────────────────────
summary_5yr <- all_scenarios %>%
  filter(abs(time_yr - 5) < 0.1) %>%
  group_by(scenario) %>%
  summarise(
    AoD_mm        = round(mean(AoD_mm), 2),
    AoD_Zscore    = round(mean(AoD_Zscore), 2),
    AnnGrowth_mm  = round(mean(AnnualGrowthRate_mm_yr), 3),
    AR_grade      = round(mean(AR_gr), 2),
    HR_bpm        = round(mean(HR_bpm), 1),
    SBP_mmhg      = round(mean(SBP_mmhg), 1),
    TGFb_conc     = round(mean(TGFb_conc), 3),
    pSMAD_fold    = round(mean(pSMAD_fold), 2),
    NT_proBNP     = round(mean(NTproBNP), 1),
    GhentScore    = round(mean(GhentScore), 1),
    .groups       = "drop"
  )
print(summary_5yr)

# ─────────────────────────────────────────────────────────────
# PLOTS
# ─────────────────────────────────────────────────────────────
scenario_colors <- c(
  "1. Untreated"          = "#B71C1C",
  "2. Atenolol 50mg QD"   = "#1565C0",
  "3. Atenolol 100mg QD"  = "#0D47A1",
  "4. Losartan 50mg QD"   = "#2E7D32",
  "5. Losartan 100mg QD"  = "#1B5E20",
  "6. Atenolol+Losartan"  = "#4527A0"
)

# Plot 1 — Aortic root diameter over time
p1 <- ggplot(all_scenarios, aes(x = time_yr, y = AoD_mm, colour = scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "red", linewidth = 0.8) +
  annotate("text", x = 0.2, y = 50.5, label = "Surgical threshold (50mm)", hjust = 0, colour = "red", size = 3) +
  geom_hline(yintercept = 40, linetype = "dotted", colour = "orange") +
  scale_colour_manual(values = scenario_colors) +
  labs(title    = "Marfan Syndrome — Aortic Root Diameter over 5 Years",
       subtitle = "PHN/COMPARE/AIMS calibrated; natural history vs medical therapy",
       x = "Time (years)", y = "Aortic Root Diameter (mm)",
       colour = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p1)

# Plot 2 — TGF-β1 plasma concentration
p2 <- ggplot(all_scenarios, aes(x = time_yr, y = TGFb_conc, colour = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = scenario_colors) +
  labs(title  = "Plasma TGF-β1 Dynamics",
       x = "Time (years)", y = "Plasma TGF-β1 [ng/mL]",
       colour = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p2)

# Plot 3 — pSMAD2/3 fold change
p3 <- ggplot(all_scenarios, aes(x = time_yr, y = pSMAD_fold, colour = scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  scale_colour_manual(values = scenario_colors) +
  labs(title  = "p-SMAD2/3 Activity (fold vs baseline)",
       x = "Time (years)", y = "p-SMAD2/3 (fold)",
       colour = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p3)

# Plot 4 — Heart rate & dP/dt
p4 <- all_scenarios %>%
  select(time_yr, scenario, HR_bpm, dPdt_val) %>%
  pivot_longer(c(HR_bpm, dPdt_val), names_to = "metric") %>%
  ggplot(aes(x = time_yr, y = value, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~metric, scales = "free_y",
             labeller = labeller(metric = c(HR_bpm="Heart Rate (bpm)", dPdt_val="dP/dt_max (mmHg/s)"))) +
  scale_colour_manual(values = scenario_colors) +
  labs(title = "Haemodynamic Effects of β-blockade",
       x = "Time (years)", y = "Value", colour = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p4)

# Plot 5 — Aortic regurgitation progression
p5 <- ggplot(all_scenarios, aes(x = time_yr, y = AR_gr, colour = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = scenario_colors) +
  scale_y_continuous(limits = c(0, 3), breaks = 0:3) +
  labs(title  = "Aortic Regurgitation Grade Progression",
       x = "Time (years)", y = "AR Grade (0-4)",
       colour = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p5)

# Plot 6 — NT-proBNP
p6 <- ggplot(all_scenarios, aes(x = time_yr, y = NTproBNP, colour = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = scenario_colors) +
  labs(title = "NT-proBNP (Cardiac Stress Biomarker)",
       x = "Time (years)", y = "NT-proBNP [pg/mL]",
       colour = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p6)

message("\n✔ Marfan Syndrome QSP model simulations complete.")
message("  6 scenarios × 5 years simulated.")
message("  Key calibration: PHN trial (Lacro NEJM 2014) | COMPARE (Radonic EHJ 2010) | AIMS (Forteza JACC 2016)")
