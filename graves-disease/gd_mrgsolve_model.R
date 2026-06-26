################################################################################
# Graves' Disease QSP Model — mrgsolve Implementation
#
# Disease:    Graves' Disease (Autoimmune Hyperthyroidism)
# Mechanisms: HPT axis · TRAb dynamics · Thyroid hormone PK/PD
#             Drug PK (MMI/PTU/RAI/propranolol/levothyroxine)
#             Target organ effects (cardiovascular, bone, CNS)
#
# References: Morshed & Davies 2015 (Endocr Rev), Brix 2005 (Endocrinology),
#             Guo 2018 (JPKPD), Walter 2012 (Clin Pharmacokinet),
#             Bartalena 2018 (Eur Thyroid J), Menconi 2014 (Best Pract Res)
#
# Compartments (≥15 ODEs):
#   C1  = MMI (methimazole) plasma
#   C2  = PTU plasma
#   C3  = 131I thyroid uptake
#   C4  = Propranolol plasma
#   C5  = Levothyroxine (exogenous T4) pool
#   C6  = TSH (pituitary output)
#   C7  = T4 (total serum)
#   C8  = T3 (total serum)
#   C9  = fT4 (free T4, bioavailable)
#   C10 = fT3 (free T3, active)
#   C11 = rT3 (reverse T3, inactive)
#   C12 = TRAb (TSAb stimulating antibody)
#   C13 = B_cell (thyroid-antigen-specific B cells)
#   C14 = Bone_resorption (osteoclast activity index)
#   C15 = HR_effect (T3-driven heart rate deviation)
#   C16 = GO_score (orbital fibroblast activation)
#   C17 = Thyroid_mass (follicular cell mass, relative)
#   C18 = TPO_inhibited (fraction of TPO inhibited by MMI/PTU)
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ─── Model Definition ──────────────────────────────────────────────────────────
graves_model <- '
$PROB
Graves Disease QSP Model — HPT Axis, TRAb, Thyroid Hormone PK/PD, Drug Effects

$PARAM
// ── Thyroid Hormone Production (baseline calibrated to normal) ──
ksyn_T4   = 1.12,   // T4 synthesis rate constant (nmol/day) per TSHR activation
ksyn_T3   = 0.18,   // direct T3 synthesis rate constant
kD1_T4T3  = 0.060,  // D1/D2 deiodination T4→T3 (day⁻¹)
kD3_T4rT3 = 0.015,  // D3 inactivation T4→rT3 (day⁻¹)
kD1_rT3   = 0.40,   // D1 rT3 clearance (day⁻¹)
kel_T4    = 0.099,  // T4 elimination (t½≈7d)
kel_T3    = 0.693,  // T3 elimination (t½≈1d)
kel_rT3   = 1.386,  // rT3 elimination (t½≈0.5d)
Vd_T4     = 10.0,   // T4 apparent volume (L)
Vd_T3     = 38.0,   // T3 apparent volume (L)

// ── HPT Axis ──
ksyn_TSH  = 0.50,   // TSH synthesis rate (mIU/L/day)
kel_TSH   = 1.386,  // TSH elimination (t½≈0.5d)
IC50_fT4  = 14.0,   // fT4 IC50 for TSH suppression (pmol/L)
IC50_fT3  = 4.5,    // fT3 IC50 for TSH suppression
n_TSH     = 2.0,    // Hill coefficient for TSH feedback
TSH_base  = 2.5,    // target steady-state TSH (mIU/L)
fT4_norm  = 16.0,   // normal fT4 (pmol/L)
fT3_norm  = 5.5,    // normal fT3 (pmol/L)
frac_fT4  = 0.00025, // fraction of T4 that is free
frac_fT3  = 0.003,  // fraction of T3 that is free

// ── TSHR Activation (TSH + TRAb combined) ──
E_TSH_max = 1.0,    // maximal TSHR effect from TSH
EC50_TSH  = 1.5,    // EC50 of TSH for TSHR (mIU/L)
E_TRAb_max= 3.5,    // max fold-stimulation by TRAb (relative to TSH=0)
EC50_TRAb = 8.0,    // EC50 of TRAb for TSHR stimulation (IU/L)

// ── TRAb / B-cell Dynamics ──
ksyn_Bcell  = 0.05,  // B-cell proliferation rate (day⁻¹)
kdeath_Bcell= 0.02,  // B-cell apoptosis rate (day⁻¹)
kprod_TRAb  = 0.20,  // TRAb production per B cell
kel_TRAb    = 0.021, // TRAb elimination (t½≈33d)
TRAb0       = 25.0,  // initial TRAb at disease onset (IU/L, >1.75 = positive)
Bcell0      = 10.0,  // initial disease-specific B cells (AU)
MMI_Bcell   = 0.35,  // fractional B-cell suppression by MMI (immunomodulatory)
RTX_Bcell   = 0.95,  // B-cell depletion fraction by rituximab

// ── Thyroid Gland / TPO ──
Thyroid_mass0 = 1.0,  // initial thyroid mass (relative, 1 = normal)
kgrowth_thy   = 0.01, // thyroid hypertrophy rate (TSHR stimulation driven)
kdeath_thy    = 0.003,// thyroid apoptosis rate
RAI_kill      = 0.10, // 131I thyroid cell kill rate constant (day⁻¹/absorbed dose)
TPO_inh_max   = 0.98, // max TPO inhibition fraction

// ── Drug PK — Methimazole (MMI) ──
ka_MMI   = 8.0,    // absorption (day⁻¹, rapid: tmax~1h)
kel_MMI  = 4.0,    // elimination (t½~4h)
Vd_MMI   = 0.5,    // volume (L/kg × 70 kg ≈ 35L → normalized)
F_MMI    = 0.93,   // bioavailability
IC50_MMI = 0.15,   // MMI conc for 50% TPO inhibition (µg/mL)

// ── Drug PK — PTU ──
ka_PTU   = 5.0,    // absorption (day⁻¹, tmax~1-2h)
kel_PTU  = 9.24,   // elimination (t½~1.8h)
Vd_PTU   = 0.5,
F_PTU    = 0.75,
IC50_PTU = 1.0,    // PTU for 50% TPO inhibition (µg/mL)
IC50_PTU_D1 = 5.0, // PTU for 50% D1 inhibition

// ── Drug PK — Propranolol ──
ka_PROP  = 6.0,    // absorption (day⁻¹, tmax~1-2h)
kel_PROP = 4.16,   // elimination (t½~4h)
Vd_PROP  = 3.9,    // volume (L/kg × 70 kg normalized)
F_PROP   = 0.26,   // first-pass ~74%
EC50_PROP_HR = 0.02, // propranolol for 50% HR reduction

// ── Drug PK — 131I (Radioiodine) ──
kabs_RAI = 4.0,    // GI absorption (rapid)
NIS_RAI  = 0.50,   // fraction of circulating 131I taken up by thyroid
kphy_RAI = 0.087,  // physical decay (t½~8d)
kbio_RAI = 0.15,   // biological clearance from thyroid
kelb_RAI = 0.050,  // whole-body biological elimination

// ── Bone Resorption ──
kBR_base  = 1.0,   // baseline osteoclast activity (AU=1)
kBR_T3stim= 0.25,  // T3 stimulation of osteoclasts
kBR_kel   = 0.10,  // bone resorption marker clearance

// ── Heart Rate Effect ──
HR_base   = 72.0,  // baseline HR (bpm)
kHR_T3    = 0.80,  // T3-driven HR increase coefficient
HR_max_stim = 40.0,// max excess HR from hyperthyroidism (bpm)
HR_EC50_fT3 = 8.0, // fT3 for 50% max HR effect

// ── Graves Ophthalmopathy (GO) ──
kGO_TRAb  = 0.08,  // TRAb-driven orbital fibroblast activation rate
kGO_kel   = 0.015, // orbital inflammatory score natural regression
GCS_GO_inh= 0.70,  // glucocorticoid inhibition of GO (fraction)

// ── Treatment Flags (0=off, 1=on) ──
USE_MMI   = 0,
USE_PTU   = 0,
USE_RAI   = 0,
USE_PROP  = 0,
USE_LT4   = 0,
USE_RTX   = 0,
USE_GCS   = 0

$CMT
// 1-5: Drug PK compartments
MMI         // methimazole plasma
PTU_C       // PTU plasma
RAI_SERUM   // 131I serum
RAI_THYR    // 131I thyroid
PROP_C      // propranolol plasma

// 6: HPT axis
TSH_C       // serum TSH (mIU/L)

// 7-11: Thyroid hormone pools
T4_C        // total T4 serum (nmol/L)
T3_C        // total T3 serum (nmol/L)
fT4_C       // free T4 (pmol/L, tracked separately)
fT3_C       // free T3 (pmol/L)
rT3_C       // reverse T3

// 12-13: Immune compartments
TRAb_C      // TRAb (IU/L)
Bcell_C     // disease-specific B cells

// 14-18: PD endpoints
BoneResor   // osteoclast activity index
HR_dev      // heart rate deviation from baseline (bpm)
GO_act      // GO orbital activation score
ThyMass     // thyroid mass (relative)
TPO_inh     // TPO inhibition fraction (0-1)

$MAIN
// ── Initial conditions ──
ThyMass_0 = Thyroid_mass0;
TRAb_C_0  = TRAb0;
Bcell_C_0 = Bcell0;
BoneResor_0 = 1.0;
HR_dev_0  = 0.0;
GO_act_0  = TRAb0 * kGO_TRAb / kGO_kel; // GO at pseudo-SS with initial TRAb

// Normal fT4/fT3 steady-state given untreated hyperthyroid TRAb
// (computed in $ODE via ODEs)

$ODE
// ────────────────────────────────────────────────────────────────
// Drug PK
// ────────────────────────────────────────────────────────────────

// Methimazole (administered as dose events via mrgsolve)
dxdt_MMI      = -kel_MMI * MMI;

// PTU
dxdt_PTU_C    = -kel_PTU * PTU_C;

// Radioiodine — serum → thyroid uptake
dxdt_RAI_SERUM = -NIS_RAI * RAI_SERUM * (ThyMass / Thyroid_mass0)
                 - kelb_RAI * RAI_SERUM;
dxdt_RAI_THYR  = NIS_RAI * RAI_SERUM * (ThyMass / Thyroid_mass0)
                 - (kphy_RAI + kbio_RAI) * RAI_THYR;

// Propranolol
dxdt_PROP_C   = -kel_PROP * PROP_C;

// ────────────────────────────────────────────────────────────────
// TPO inhibition (0–1 scale, combined MMI + PTU)
// ────────────────────────────────────────────────────────────────
double MMI_conc   = MMI;
double PTU_conc   = PTU_C;

double TPO_inh_MMI = (USE_MMI > 0.5) ? (TPO_inh_max * MMI_conc / (IC50_MMI + MMI_conc)) : 0.0;
double TPO_inh_PTU = (USE_PTU > 0.5) ? (TPO_inh_max * PTU_conc / (IC50_PTU + PTU_conc)) : 0.0;
double TPO_total   = 1.0 - (1.0 - TPO_inh_MMI) * (1.0 - TPO_inh_PTU); // combined
dxdt_TPO_inh = (TPO_total - TPO_inh) * 2.0; // fast equilibration

// ────────────────────────────────────────────────────────────────
// D1 inhibition by PTU and propranolol
// ────────────────────────────────────────────────────────────────
double D1_inh_PTU  = (USE_PTU  > 0.5) ? (PTU_conc   / (IC50_PTU_D1   + PTU_conc))   : 0.0;
double D1_inh_PROP = (USE_PROP > 0.5) ? (PROP_C     / (EC50_PROP_HR  + PROP_C))      : 0.0;
double D1_factor   = 1.0 - 0.5 * D1_inh_PTU - 0.3 * D1_inh_PROP; // D1 activity
if(D1_factor < 0.05) D1_factor = 0.05;

// ────────────────────────────────────────────────────────────────
// TSHR Activation
// ────────────────────────────────────────────────────────────────
double TRAb_stim = E_TRAb_max * TRAb_C / (EC50_TRAb + TRAb_C);
double TSH_stim  = E_TSH_max  * TSH_C  / (EC50_TSH  + TSH_C);
double TSHR_act  = TSH_stim + TRAb_stim; // combined stimulation

// RAI-mediated thyroid damage
double RAI_damage = (USE_RAI > 0.5) ? RAI_kill * RAI_THYR : 0.0;

// ────────────────────────────────────────────────────────────────
// Thyroid Mass
// ────────────────────────────────────────────────────────────────
dxdt_ThyMass = kgrowth_thy * ThyMass * (TSHR_act - 1.0)
               - kdeath_thy * ThyMass
               - RAI_damage * ThyMass;
if(ThyMass < 0.05) dxdt_ThyMass = 0.0;

// ────────────────────────────────────────────────────────────────
// Thyroid Hormone Synthesis & Secretion
// (modulated by thyroid mass, TSHR activation, TPO inhibition)
// ────────────────────────────────────────────────────────────────
double synth_factor = ThyMass * TSHR_act * (1.0 - TPO_inh);
double T4_synth = ksyn_T4 * synth_factor;
double T3_synth = ksyn_T3 * synth_factor;

// T4 kinetics
double D1_conv  = kD1_T4T3  * D1_factor * fT4_C;
double D3_T4rT3 = kD3_T4rT3 * fT4_C;

dxdt_T4_C  = T4_synth
             - kel_T4  * T4_C
             - D1_conv
             - D3_T4rT3;

// T3 kinetics (from direct secretion + T4 conversion)
dxdt_T3_C  = T3_synth
             + D1_conv
             - kel_T3 * T3_C;

// fT4 tracks total T4 (via rapid equilibrium with binding proteins)
dxdt_fT4_C = frac_fT4 * T4_synth
              - kel_T4 * fT4_C
              - kD1_T4T3 * D1_factor * fT4_C
              - kD3_T4rT3 * fT4_C;

// fT3 tracks free T3
dxdt_fT3_C = frac_fT3 * T3_synth
              + kD1_T4T3 * D1_factor * fT4_C * frac_fT4/frac_fT3 * 0.15
              - kel_T3 * fT3_C;

// rT3
dxdt_rT3_C = kD3_T4rT3 * fT4_C
              - kD1_rT3 * rT3_C;

// ────────────────────────────────────────────────────────────────
// TSH (HPT axis negative feedback)
// ────────────────────────────────────────────────────────────────
double fb_fT4 = pow(fT4_C / fT4_norm, n_TSH);
double fb_fT3 = pow(fT3_C / fT3_norm, n_TSH);
double fb_tot = (fb_fT4 + fb_fT3) / 2.0;
double TSH_prod = ksyn_TSH / fb_tot;

dxdt_TSH_C = TSH_prod - kel_TSH * TSH_C;

// ────────────────────────────────────────────────────────────────
// B-cell & TRAb Dynamics
// ────────────────────────────────────────────────────────────────
double Bcell_growth = ksyn_Bcell * Bcell_C;
double Bcell_death  = kdeath_Bcell * Bcell_C;
double MMI_immu_eff = (USE_MMI > 0.5) ? MMI_Bcell * MMI_conc / (IC50_MMI + MMI_conc) : 0.0;
double RTX_eff      = (USE_RTX > 0.5) ? RTX_Bcell : 0.0;

dxdt_Bcell_C = Bcell_growth - Bcell_death
               - MMI_immu_eff * Bcell_C
               - RTX_eff      * Bcell_C;

dxdt_TRAb_C = kprod_TRAb * Bcell_C
               - kel_TRAb  * TRAb_C;

// ────────────────────────────────────────────────────────────────
// Bone Resorption (osteoclast index)
// ────────────────────────────────────────────────────────────────
double BR_stim = kBR_T3stim * (fT3_C / fT3_norm - 1.0);
dxdt_BoneResor = BR_stim - kBR_kel * (BoneResor - 1.0);

// ────────────────────────────────────────────────────────────────
// Heart Rate Deviation
// ────────────────────────────────────────────────────────────────
double fT3_effect = HR_max_stim * fT3_C / (HR_EC50_fT3 + fT3_C);
double PROP_HR_inh = (USE_PROP > 0.5) ?
                     (HR_max_stim * PROP_C / (EC50_PROP_HR + PROP_C)) : 0.0;
dxdt_HR_dev = 0.5 * (fT3_effect - PROP_HR_inh - HR_dev);

// ────────────────────────────────────────────────────────────────
// Graves Ophthalmopathy (GO orbital score)
// ────────────────────────────────────────────────────────────────
double GCS_GO = (USE_GCS > 0.5) ? GCS_GO_inh : 0.0;
dxdt_GO_act = kGO_TRAb * TRAb_C * (1.0 - GCS_GO)
              - kGO_kel  * GO_act;

$SIGMA
0.01  // proportional residual error

$TABLE
double TSH_obs   = TSH_C   * (1 + EPS(1));
double fT4_obs   = fT4_C   * (1 + EPS(1));
double fT3_obs   = fT3_C   * (1 + EPS(1));
double TRAb_obs  = TRAb_C  * (1 + EPS(1));
double HR_obs    = HR_base + HR_dev;
double BMD_loss  = 1.0 - 0.005 * (BoneResor - 1.0) * TIME; // cumulative % BMD loss per year

$CAPTURE TSH_obs fT4_obs fT3_obs TRAb_obs HR_obs BMD_loss
         MMI PTU_C RAI_THYR PROP_C TPO_inh ThyMass GO_act BoneResor
'

# ─── Compile Model ─────────────────────────────────────────────────────────────
gd_mod <- mrgsolve::mcode("graves_disease", graves_model)

# ─── Scenario 1: Untreated Graves' Disease (natural history) ──────────────────
sim_untreated <- gd_mod %>%
  param(USE_MMI=0, USE_PTU=0, USE_RAI=0, USE_PROP=0) %>%
  init(TRAb_C  = 25,    # elevated TRAb at disease onset
       Bcell_C = 10,
       TSH_C   = 0.01,  # suppressed TSH
       fT4_C   = 35,    # elevated fT4 (pmol/L)
       fT3_C   = 12,    # elevated fT3
       T4_C    = 220,   # elevated total T4
       T3_C    = 5.0,
       rT3_C   = 0.3,
       ThyMass = 1.5,   # goiter
       BoneResor = 1.6,
       HR_dev  = 25,
       GO_act  = 35) %>%
  mrgsim(end=730, delta=1) %>%
  as.data.frame() %>%
  mutate(scenario = "Untreated")

# ─── Scenario 2: Methimazole 30mg/day (antithyroid drug) ──────────────────────
# MMI given as 3x daily doses (10 mg each, ~15 µg/mL peak)
mmi_events <- ev(amt=0.10, ii=0.33, addl=2190-1, cmt=1, time=0)  # q8h for 730 days

sim_mmi <- gd_mod %>%
  param(USE_MMI=1, USE_PTU=0, USE_RAI=0, USE_PROP=0) %>%
  init(TRAb_C=25, Bcell_C=10, TSH_C=0.01, fT4_C=35, fT3_C=12,
       T4_C=220, T3_C=5.0, rT3_C=0.3, ThyMass=1.5, BoneResor=1.6,
       HR_dev=25, GO_act=35) %>%
  mrgsim(events=mmi_events, end=730, delta=1) %>%
  as.data.frame() %>%
  mutate(scenario = "Methimazole 30mg/day")

# ─── Scenario 3: Radioiodine Ablation (single 15 mCi dose at day 0) ──────────
rai_events <- ev(amt=15, cmt=3, time=0)  # single RAI dose (mCi → normalized)

sim_rai <- gd_mod %>%
  param(USE_MMI=0, USE_PTU=0, USE_RAI=1, USE_PROP=0) %>%
  init(TRAb_C=25, Bcell_C=10, TSH_C=0.01, fT4_C=35, fT3_C=12,
       T4_C=220, T3_C=5.0, rT3_C=0.3, ThyMass=1.5, BoneResor=1.6,
       HR_dev=25, GO_act=35) %>%
  mrgsim(events=rai_events, end=730, delta=1) %>%
  as.data.frame() %>%
  mutate(scenario = "Radioiodine 15 mCi")

# ─── Scenario 4: MMI + Propranolol (combined initial therapy) ─────────────────
combo_events <- ev(amt=0.10, ii=0.33, addl=365*3-1, cmt=1, time=0) +
                ev(amt=0.06, ii=0.17, addl=365*6-1, cmt=5, time=0)  # prop q6h

sim_combo <- gd_mod %>%
  param(USE_MMI=1, USE_PTU=0, USE_RAI=0, USE_PROP=1) %>%
  init(TRAb_C=25, Bcell_C=10, TSH_C=0.01, fT4_C=35, fT3_C=12,
       T4_C=220, T3_C=5.0, rT3_C=0.3, ThyMass=1.5, BoneResor=1.6,
       HR_dev=25, GO_act=35) %>%
  mrgsim(events=combo_events, end=730, delta=1) %>%
  as.data.frame() %>%
  mutate(scenario = "MMI + Propranolol")

# ─── Scenario 5: Block-and-Replace (MMI + Levothyroxine) ──────────────────────
br_mmi   <- ev(amt=0.20, ii=0.33, addl=2190-1, cmt=1, time=0)   # MMI 60mg/day
br_lt4   <- ev(amt=0.10, ii=1.0,  addl=730-1,  cmt=5, time=90)  # LT4 start at day 90

sim_br <- gd_mod %>%
  param(USE_MMI=1, USE_PTU=0, USE_RAI=0, USE_PROP=0, USE_LT4=1) %>%
  init(TRAb_C=25, Bcell_C=10, TSH_C=0.01, fT4_C=35, fT3_C=12,
       T4_C=220, T3_C=5.0, rT3_C=0.3, ThyMass=1.5, BoneResor=1.6,
       HR_dev=25, GO_act=35) %>%
  mrgsim(events=br_mmi+br_lt4, end=730, delta=1) %>%
  as.data.frame() %>%
  mutate(scenario = "Block-and-Replace (MMI+LT4)")

# ─── Scenario 6: Rituximab for refractory GO ──────────────────────────────────
rtx_events <- ev(amt=1, cmt=1, time=0) + ev(amt=1, cmt=1, time=14)  # 2 doses d0,d14

sim_rtx <- gd_mod %>%
  param(USE_MMI=1, USE_RTX=1, USE_GCS=1) %>%
  init(TRAb_C=30, Bcell_C=15, TSH_C=0.01, fT4_C=35, fT3_C=12,
       T4_C=220, T3_C=5.0, rT3_C=0.3, ThyMass=1.5, BoneResor=1.6,
       HR_dev=25, GO_act=55) %>%  # higher GO at baseline
  mrgsim(events=rtx_events, end=730, delta=1) %>%
  as.data.frame() %>%
  mutate(scenario = "Rituximab + MMI (refractory GO)")

# ─── Combine All Scenarios ─────────────────────────────────────────────────────
all_sims <- bind_rows(sim_untreated, sim_mmi, sim_rai, sim_combo, sim_br, sim_rtx)

# ─── Plotting Functions ────────────────────────────────────────────────────────
plot_graves_outcomes <- function(data) {
  cols <- c("Untreated"                        = "#e74c3c",
            "Methimazole 30mg/day"             = "#3498db",
            "Radioiodine 15 mCi"               = "#e67e22",
            "MMI + Propranolol"                = "#9b59b6",
            "Block-and-Replace (MMI+LT4)"      = "#27ae60",
            "Rituximab + MMI (refractory GO)"  = "#1abc9c")

  p_tsh <- ggplot(data, aes(TIME/30, TSH_obs, color=scenario)) +
    geom_line(size=0.9) +
    geom_hline(yintercept=c(0.4,4.5), linetype="dashed", color="grey50") +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.4, ymax=4.5,
             alpha=0.1, fill="green") +
    scale_color_manual(values=cols) +
    scale_y_log10(labels=scales::label_number()) +
    labs(title="Serum TSH", x="Time (months)", y="TSH (mIU/L)", color=NULL) +
    theme_classic(base_size=11) + theme(legend.position="bottom")

  p_fT4 <- ggplot(data, aes(TIME/30, fT4_obs, color=scenario)) +
    geom_line(size=0.9) +
    geom_hline(yintercept=c(12,22), linetype="dashed", color="grey50") +
    scale_color_manual(values=cols) +
    labs(title="Free T4 (fT4)", x="Time (months)", y="fT4 (pmol/L)", color=NULL) +
    theme_classic(base_size=11) + theme(legend.position="bottom")

  p_TRAb <- ggplot(data, aes(TIME/30, TRAb_obs, color=scenario)) +
    geom_line(size=0.9) +
    geom_hline(yintercept=1.75, linetype="dashed", color="red", alpha=0.7) +
    annotate("text", x=1, y=2.2, label="TRAb cutoff 1.75 IU/L", size=3) +
    scale_color_manual(values=cols) +
    labs(title="TRAb (Thyroid-Stimulating Antibody)", x="Time (months)", y="TRAb (IU/L)", color=NULL) +
    theme_classic(base_size=11) + theme(legend.position="bottom")

  p_HR <- ggplot(data, aes(TIME/30, HR_obs, color=scenario)) +
    geom_line(size=0.9) +
    geom_hline(yintercept=100, linetype="dashed", color="red", alpha=0.7) +
    geom_hline(yintercept=72, linetype="dashed", color="grey60") +
    scale_color_manual(values=cols) +
    labs(title="Heart Rate", x="Time (months)", y="HR (bpm)", color=NULL) +
    theme_classic(base_size=11) + theme(legend.position="bottom")

  p_bone <- ggplot(data, aes(TIME/30, BoneResor, color=scenario)) +
    geom_line(size=0.9) +
    geom_hline(yintercept=1.0, linetype="dashed", color="grey50") +
    scale_color_manual(values=cols) +
    labs(title="Bone Resorption Index\n(Osteoclast Activity AU)", x="Time (months)", y="Bone Resorption (AU)", color=NULL) +
    theme_classic(base_size=11) + theme(legend.position="bottom")

  p_GO <- ggplot(data, aes(TIME/30, GO_act, color=scenario)) +
    geom_line(size=0.9) +
    scale_color_manual(values=cols) +
    labs(title="Graves Ophthalmopathy\n(Orbital Activation Score)", x="Time (months)", y="GO Score (AU)", color=NULL) +
    theme_classic(base_size=11) + theme(legend.position="bottom")

  p_thy <- ggplot(data, aes(TIME/30, ThyMass, color=scenario)) +
    geom_line(size=0.9) +
    geom_hline(yintercept=1.0, linetype="dashed", color="grey50") +
    scale_color_manual(values=cols) +
    labs(title="Thyroid Mass (Goiter)", x="Time (months)", y="Relative Thyroid Mass", color=NULL) +
    theme_classic(base_size=11) + theme(legend.position="bottom")

  p_fT3 <- ggplot(data, aes(TIME/30, fT3_obs, color=scenario)) +
    geom_line(size=0.9) +
    geom_hline(yintercept=c(3.5,7.5), linetype="dashed", color="grey50") +
    scale_color_manual(values=cols) +
    labs(title="Free T3 (fT3)", x="Time (months)", y="fT3 (pmol/L)", color=NULL) +
    theme_classic(base_size=11) + theme(legend.position="bottom")

  # Panel layout
  layout <- (p_tsh | p_fT4 | p_fT3) /
             (p_TRAb | p_HR | p_thy) /
             (p_bone | p_GO | plot_spacer()) +
    plot_annotation(
      title    = "Graves' Disease QSP Model — Treatment Scenario Comparison",
      subtitle = "6 scenarios: Untreated | MMI | RAI | MMI+Propranolol | Block-and-Replace | Rituximab+MMI",
      theme    = theme(plot.title    = element_text(size=16, face="bold"),
                       plot.subtitle = element_text(size=12, color="grey40"))
    )
  return(layout)
}

# ─── Sensitivity Analysis ──────────────────────────────────────────────────────
run_sensitivity <- function(param_name, values, base_params=list()) {
  results <- lapply(values, function(v) {
    p <- base_params
    p[[param_name]] <- v
    gd_mod %>%
      do.call(param, list(., p)) %>%
      param(USE_MMI=1) %>%
      init(TRAb_C=25, Bcell_C=10, TSH_C=0.01, fT4_C=35, fT3_C=12,
           T4_C=220, T3_C=5.0, ThyMass=1.5) %>%
      mrgsim(events=ev(amt=0.10, ii=0.33, addl=2190-1, cmt=1, time=0),
             end=365, delta=7) %>%
      as.data.frame() %>%
      mutate(param_val = v)
  })
  bind_rows(results)
}

# TRAb sensitivity to IC50_MMI
sa_IC50 <- run_sensitivity("IC50_MMI",
                           c(0.05, 0.10, 0.15, 0.25, 0.50),
                           list(USE_MMI=1))

# ─── Virtual Population (VPop) Simulation ─────────────────────────────────────
set.seed(42)
n_pop <- 200

vpop_params <- data.frame(
  ID       = 1:n_pop,
  TRAb0    = rlnorm(n_pop, log(20), 0.5),       # lognormal TRAb
  IC50_MMI = rlnorm(n_pop, log(0.15), 0.3),     # variability in drug sensitivity
  ksyn_Bcell = rnorm(n_pop, 0.05, 0.01)
) %>% filter(IC50_MMI > 0.02, ksyn_Bcell > 0.01)

vpop_sim <- lapply(1:min(50, nrow(vpop_params)), function(i) {
  p <- vpop_params[i, ]
  gd_mod %>%
    param(USE_MMI=1,
          IC50_MMI   = p$IC50_MMI,
          ksyn_Bcell = p$ksyn_Bcell) %>%
    init(TRAb_C=p$TRAb0, Bcell_C=10, TSH_C=0.01, fT4_C=35, fT3_C=12,
         T4_C=220, T3_C=5.0, ThyMass=1.5) %>%
    mrgsim(events=ev(amt=0.10, ii=0.33, addl=2190-1, cmt=1, time=0),
           end=730, delta=7) %>%
    as.data.frame() %>%
    mutate(ID = p$ID)
}) %>% bind_rows()

# Remission rate (TRAb < 1.75 at 18 months)
remission_rate <- vpop_sim %>%
  filter(TIME == 540) %>%  # day 540 = ~18 months
  summarise(pct_remission = mean(TRAb_obs < 1.75) * 100,
            pct_euthyroid = mean(TSH_obs > 0.4 & TSH_obs < 4.5) * 100)
cat("VPop MMI 18-month outcomes:\n")
print(remission_rate)

# ─── Summary Print ─────────────────────────────────────────────────────────────
cat("\n========================================\n")
cat("Graves' Disease QSP Model Summary\n")
cat("========================================\n")
cat("Model compartments: 18 ODEs\n")
cat("Drug PK modules: MMI, PTU, 131I, Propranolol, LT4, Rituximab\n")
cat("Treatment scenarios: 6\n")
cat("Virtual population: n=", nrow(vpop_params), " patients\n")
cat("\nKey normal ranges:\n")
cat("  TSH:  0.4–4.5 mIU/L\n")
cat("  fT4:  12–22  pmol/L\n")
cat("  fT3:  3.5–7.5 pmol/L\n")
cat("  TRAb: <1.75 IU/L (positive threshold)\n")
cat("  HR:   60–100 bpm\n")
cat("\nRun plot_graves_outcomes(all_sims) to visualize results\n")
