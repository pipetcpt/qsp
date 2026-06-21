## ============================================================
## Sickle Cell Disease (SCD) — QSP Model (mrgsolve)
## ============================================================
## Compartments (24 ODEs):
##   Drug PK  : HU_gut, HU_plasma, VOX_plasma, VOX_RBC,
##              CRIZ_C, CRIZ_P, LG
##   RBC/Hb   : CFU_E, RET, RBC_S, RBC_N, HbF_frac, Hgb,
##              free_Hb, Haptoglobin
##   Biomarkers: LDH, Bilirubin, NO, P_selectin,
##              VOC, NADH, Iron, TRV, eGFR
##
## Treatments:
##   1. No treatment (baseline)
##   2. Hydroxyurea (HU) monotherapy
##   3. Voxelotor (VOX) monotherapy
##   4. Crizanlizumab (CRIZ) monotherapy
##   5. L-Glutamine monotherapy
##   6. HU + VOX combination
##   7. HU + VOX + CRIZ triple therapy
##
## Key References:
##   Charache 1995 (NEJM, hydroxyurea RCT MSH trial)
##   Vichinsky 2019 (NEJM, voxelotor HOPE trial)
##   Ataga 2017 (NEJM, crizanlizumab SUSTAIN trial)
##   Niihara 2018 (NEJM, L-glutamine phase III)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ── Model code ───────────────────────────────────────────────
scd_model_code <- '
$PROB Sickle Cell Disease QSP Model — mrgsolve
  24-compartment ODE system: Drug PK + RBC/Hb dynamics + Biomarkers

$PARAM
  // ─── Hydroxyurea PK ───
  ka_HU   = 2.0,    // Absorption rate constant (h⁻¹)
  CL_HU   = 3.5,    // Clearance (L/h)
  Vd_HU   = 28.0,   // Volume of distribution (L)
  F_HU    = 0.80,   // Bioavailability
  MW_HU   = 76.1,   // Molecular weight (g/mol)

  // ─── Voxelotor PK ───
  ka_VOX  = 0.45,   // Absorption rate constant (h⁻¹)
  CL_VOX  = 2.1,    // Apparent clearance (L/h)
  Vd_VOX  = 97.0,   // Apparent volume of distribution (L)
  kon_VOX = 0.15,   // RBC binding on-rate (h⁻¹·(μM)⁻¹)
  koff_VOX= 0.005,  // RBC binding off-rate (h⁻¹)
  RBC_capacity = 50, // RBC binding capacity (μM equivalent)

  // ─── Crizanlizumab PK (2-compartment mAb) ───
  CL_CRIZ = 0.008,  // Clearance (L/h)
  Vc_CRIZ = 3.5,    // Central volume (L)
  Vp_CRIZ = 2.8,    // Peripheral volume (L)
  k12_CRIZ= 0.015,  // Distribution (central→peripheral, h⁻¹)
  k21_CRIZ= 0.012,  // Return (peripheral→central, h⁻¹)

  // ─── L-Glutamine PK ───
  ka_LG   = 1.2,    // Absorption rate constant (h⁻¹)
  CL_LG   = 45.0,   // Rapid clearance (L/h)
  Vd_LG   = 18.0,   // Volume of distribution (L)

  // ─── Erythropoiesis baseline ───
  kprod_CFU = 0.015, // CFU-E production rate (h⁻¹)
  kdiff_CFU = 0.08,  // CFU-E → RET differentiation (h⁻¹)
  kdiff_RET = 0.035, // RET → RBC maturation (h⁻¹)
  kdeath_S  = 0.0055,// Sickle RBC elimination (h⁻¹, t½≈5.3 d)
  kdeath_N  = 0.00014,//Normal RBC elimination (h⁻¹, t½≈30 d)
  RET_0     = 45.0,  // Baseline reticulocyte (10⁹ cells/L)
  RBC_S_0   = 1800,  // Baseline sickle RBC (10⁹/L)
  RBC_N_0   = 200,   // Baseline normal RBC (10⁹/L, HbF/HbA)
  CFU_E_0   = 12.0,  // Baseline CFU-E (10⁹/L)

  // ─── HbF dynamics ───
  HbF0      = 0.07,  // Baseline HbF fraction (7%)
  kHbF_deg  = 0.002, // HbF fraction degradation rate (h⁻¹)
  EC50_HU_HbF = 15.0,// EC50 of HU for HbF induction (μM)
  Emax_HU_HbF = 0.18, // Maximum HU-induced HbF increase

  // ─── Hemoglobin dynamics ───
  Hgb_0     = 8.5,   // Baseline Hgb (g/dL)
  kprod_Hgb = 0.003, // Hgb production rate (g/dL/h)
  kdeg_Hgb  = 0.00035,// Baseline Hgb loss rate (h⁻¹)

  // ─── Hemolysis markers ───
  LDH_0     = 480.0, // Baseline LDH (U/L, SCD typical)
  kprod_LDH = 0.15,  // LDH production (h⁻¹ × baseline)
  kdeg_LDH  = 0.08,  // LDH degradation (h⁻¹)
  Bili_0    = 42.0,  // Baseline bilirubin (μmol/L)
  kdeg_Bili = 0.02,  // Bilirubin clearance (h⁻¹)
  Hp_0      = 0.25,  // Baseline haptoglobin (g/L, low in SCD)
  kdeg_Hp   = 0.03,  // Haptoglobin consumption rate (h⁻¹)
  kprod_Hp  = 0.008, // Haptoglobin synthesis (g/L/h)
  kfree_Hb  = 0.12,  // Free Hb clearance (h⁻¹)

  // ─── NO & Vascular ───
  NO_0      = 0.40,  // Baseline NO index (arbitrary units, 0–1)
  kNO_prod  = 0.05,  // NO production rate (h⁻¹)
  kNO_scav  = 0.8,   // NO scavenging by free Hb (per μg/mL·h)

  // ─── P-Selectin dynamics ───
  Psel_0    = 1.0,   // Baseline P-selectin (relative, 1=normal SCD)
  kPsel_deg = 0.04,  // P-selectin turnover (h⁻¹)
  EC50_CRIZ_Psel = 0.5, // CRIZ EC50 for P-sel inhibition (μg/mL)
  Imax_CRIZ_Psel = 0.80,// Maximum P-sel inhibition by CRIZ

  // ─── VOC dynamics ───
  VOC_0     = 0.0004,// Baseline VOC rate (crises/h → ~3.5/year)
  kVOC_Psel = 0.8,   // P-selectin contribution to VOC
  kVOC_NO   = 0.6,   // NO deficiency contribution to VOC
  kVOC_HbF  = 0.5,   // HbF fraction modulates VOC

  // ─── NADH/Redox ───
  NADH_0    = 0.55,  // Baseline NADH (relative units)
  kNADH_prod= 0.04,  // NADH production
  kNADH_ox  = 0.06,  // NADH oxidation rate

  // ─── Iron / Ferritin ───
  Iron_0    = 900.0, // Baseline ferritin (μg/L, elevated in SCD)
  kIron_acc = 0.002, // Iron accumulation from hemolysis (h⁻¹)
  kIron_deg = 0.001, // Iron mobilization/excretion (h⁻¹)

  // ─── TRV (Pulmonary HTN proxy) ───
  TRV_0     = 2.65,  // Baseline TRV (m/s, borderline elevated)
  kTRV_prod = 0.0001,// TRV increase rate
  kTRV_deg  = 0.0002,// TRV decrease rate

  // ─── eGFR (Renal function) ───
  eGFR_0    = 105.0, // Baseline eGFR (mL/min/1.73m², SCD often hyperfiltrates)
  keGFR_dec = 0.000005,// Slow CKD progression rate (h⁻¹)

  // ─── Drug doses (on/off flags) ───
  DOSE_HU   = 0,     // HU dose (mg/kg/day)
  BW        = 65,    // Body weight (kg)
  DOSE_VOX  = 0,     // Voxelotor dose (1500 mg/day = 1 if on)
  DOSE_CRIZ = 0,     // Crizanlizumab (1 if on, q4w)
  DOSE_LG   = 0      // L-Glutamine (1 if on)

$CMT
  // Drug PK (7 compartments)
  HU_gut HU_plasma VOX_plasma VOX_RBC CRIZ_C CRIZ_P LG

  // Erythropoiesis & RBC (7 compartments)
  CFU_E RET RBC_S RBC_N HbF_frac Hgb free_Hb

  // Biomarkers (10 compartments)
  Haptoglobin LDH Bilirubin NO P_selectin VOC NADH Iron TRV eGFR

$INIT
  HU_gut    = 0,
  HU_plasma = 0,
  VOX_plasma= 0,
  VOX_RBC   = 0,
  CRIZ_C    = 0,
  CRIZ_P    = 0,
  LG        = 0,

  CFU_E     = 12.0,
  RET       = 45.0,
  RBC_S     = 1800,
  RBC_N     = 200,
  HbF_frac  = 0.07,
  Hgb       = 8.5,
  free_Hb   = 8.0,

  Haptoglobin= 0.25,
  LDH       = 480.0,
  Bilirubin = 42.0,
  NO        = 0.40,
  P_selectin= 1.0,
  VOC       = 0.0004,
  NADH      = 0.55,
  Iron      = 900.0,
  TRV       = 2.65,
  eGFR      = 105.0

$ODE

  // ──────────────────────────────────────────────────────
  // Derived quantities
  // ──────────────────────────────────────────────────────
  double HU_mg   = DOSE_HU * BW;         // Total HU dose (mg/day)
  double HU_umol = HU_mg * 1000 / MW_HU; // Convert to μmol/day → per h below
  double HU_rate = HU_umol / 24.0;       // μmol/h input rate

  // VOX: 1500 mg/day oral (MW = 392.4 g/mol)
  double VOX_mg   = DOSE_VOX * 1500.0;
  double VOX_umol = VOX_mg * 1000 / 392.4;
  double VOX_rate = VOX_umol / 24.0;

  // LG: 5000 mg BID = 10g/day (MW = 146.1)
  double LG_mg   = DOSE_LG * 10000.0;
  double LG_umol = LG_mg * 1000 / 146.1;
  double LG_rate = LG_umol / 24.0;

  // ──────────────────────────────────────────────────────
  // Drug PK ODEs
  // ──────────────────────────────────────────────────────

  // Hydroxyurea (1-compartment oral)
  dxdt_HU_gut    = F_HU * HU_rate - ka_HU * HU_gut;
  dxdt_HU_plasma = ka_HU * HU_gut / Vd_HU
                   - (CL_HU / Vd_HU) * HU_plasma;
  double HU_uM   = HU_plasma;  // μM in plasma

  // Voxelotor (1-cmpt oral + RBC binding)
  double VOX_free = VOX_plasma;
  dxdt_VOX_plasma = VOX_rate / Vd_VOX
                   - (CL_VOX / Vd_VOX) * VOX_plasma
                   - kon_VOX * VOX_free * (RBC_capacity - VOX_RBC)
                   + koff_VOX * VOX_RBC;
  dxdt_VOX_RBC   = kon_VOX * VOX_free * (RBC_capacity - VOX_RBC)
                   - koff_VOX * VOX_RBC;

  // Crizanlizumab (2-compartment IV mAb, q4w dosing handled in events)
  dxdt_CRIZ_C    = - (CL_CRIZ / Vc_CRIZ) * CRIZ_C
                   - k12_CRIZ * CRIZ_C
                   + k21_CRIZ * CRIZ_P;
  dxdt_CRIZ_P    = k12_CRIZ * CRIZ_C - k21_CRIZ * CRIZ_P;
  double CRIZ_ugml = CRIZ_C * 148000 / 1e6;  // Convert nM to μg/mL (MW~148 kDa)

  // L-Glutamine (rapid 1-compartment)
  dxdt_LG        = LG_rate / Vd_LG - (CL_LG / Vd_LG) * LG;

  // ──────────────────────────────────────────────────────
  // Drug PD effects
  // ──────────────────────────────────────────────────────

  // HU → HbF induction (Emax Hill model, Hill=1.5)
  double E_HU_HbF = Emax_HU_HbF * pow(HU_uM, 1.5) /
                    (pow(EC50_HU_HbF, 1.5) + pow(HU_uM, 1.5));

  // VOX → increases O₂ affinity → reduces deoxy-HbS fraction
  // Proportional to RBC-bound VOX
  double VOX_occ  = VOX_RBC / (RBC_capacity + 1e-6);  // 0–1
  double VOX_polymerInhib = 0.45 * VOX_occ;           // Max 45% sickling reduction

  // CRIZ → P-Selectin inhibition (Imax model)
  double I_CRIZ_Psel = Imax_CRIZ_Psel * CRIZ_ugml /
                       (EC50_CRIZ_Psel + CRIZ_ugml);

  // LG → NADH support
  double E_LG_NADH = 0.15 * DOSE_LG * (LG / (LG + 100));

  // ──────────────────────────────────────────────────────
  // HbF fraction (key protective factor)
  // ──────────────────────────────────────────────────────
  dxdt_HbF_frac = (HbF0 + E_HU_HbF) * kHbF_deg
                  - kHbF_deg * HbF_frac;

  // ──────────────────────────────────────────────────────
  // Erythropoiesis
  // ──────────────────────────────────────────────────────
  // EPO feedback: anemia drives erythroid output
  double EPO_stim = pow(Hgb_0 / (Hgb + 0.01), 1.5);  // Compensatory

  dxdt_CFU_E = kprod_CFU * EPO_stim - kdiff_CFU * CFU_E;
  dxdt_RET   = kdiff_CFU * CFU_E - kdiff_RET * RET;

  // RBC populations (sickle vs normal/HbF)
  double frac_sickle = 1 - HbF_frac;
  double frac_normal = HbF_frac;

  // Voxelotor reduces effective sickling via polymerization inhibition
  double eff_sickle = frac_sickle * (1 - VOX_polymerInhib);

  dxdt_RBC_S = kdiff_RET * RET * eff_sickle
               - kdeath_S * RBC_S;

  dxdt_RBC_N = kdiff_RET * RET * (1 - eff_sickle)
               - kdeath_N * RBC_N;

  // ──────────────────────────────────────────────────────
  // Hemoglobin dynamics
  // ──────────────────────────────────────────────────────
  double RBC_total = RBC_S + RBC_N;
  double Hgb_prod  = kprod_Hgb * (RBC_total / (RBC_S_0 + RBC_N_0));
  double Hgb_loss  = kdeg_Hgb * Hgb;

  dxdt_Hgb = Hgb_prod - Hgb_loss;

  // ──────────────────────────────────────────────────────
  // Hemolysis: Free Hb
  // ──────────────────────────────────────────────────────
  double hemolysis_rate = kdeath_S * RBC_S * 0.015;  // 1.5% per lysis event → free Hb
  dxdt_free_Hb = hemolysis_rate - kfree_Hb * free_Hb
                 - kdeg_Hp * Haptoglobin * free_Hb;

  // ──────────────────────────────────────────────────────
  // Haptoglobin (consumed by free Hb)
  // ──────────────────────────────────────────────────────
  dxdt_Haptoglobin = kprod_Hp - kdeg_Hp * free_Hb * Haptoglobin;

  // ──────────────────────────────────────────────────────
  // LDH (rises with hemolysis)
  // ──────────────────────────────────────────────────────
  double LDH_prod = kprod_LDH * (hemolysis_rate / (kdeath_S * RBC_S_0 * 0.015 + 1e-6))
                    * LDH_0;
  dxdt_LDH = LDH_prod - kdeg_LDH * LDH;

  // ──────────────────────────────────────────────────────
  // Bilirubin (unconjugated, from heme catabolism)
  // ──────────────────────────────────────────────────────
  double bili_prod = 0.006 * hemolysis_rate * Bili_0;
  dxdt_Bilirubin = bili_prod - kdeg_Bili * Bilirubin;

  // ──────────────────────────────────────────────────────
  // Nitric Oxide (scavenged by free Hb)
  // ──────────────────────────────────────────────────────
  dxdt_NO = kNO_prod - kNO_scav * free_Hb * NO
            - 0.01 * NO;  // baseline turnover

  // ──────────────────────────────────────────────────────
  // P-Selectin (upregulated by inflammation/ischemia)
  // ──────────────────────────────────────────────────────
  double Psel_stim = 1.0 + 0.3 * (1 - NO / NO_0);  // NO deficiency → Psel↑
  double Psel_inhib = 1 - I_CRIZ_Psel;

  dxdt_P_selectin = kPsel_deg * Psel_stim * Psel_inhib
                    - kPsel_deg * P_selectin;

  // ──────────────────────────────────────────────────────
  // VOC rate (vaso-occlusive crisis rate)
  // ──────────────────────────────────────────────────────
  double VOC_driver = P_selectin * kVOC_Psel
                      * (1.0 / (NO + 0.1)) * kVOC_NO
                      * (1 - HbF_frac * kVOC_HbF);
  double VOC_norm   = Psel_0 * kVOC_Psel
                      * (1.0 / (NO_0 + 0.1)) * kVOC_NO
                      * (1 - HbF0 * kVOC_HbF);

  dxdt_VOC = VOC_0 * (VOC_driver / VOC_norm + 1e-9) - 0.05 * VOC;

  // ──────────────────────────────────────────────────────
  // NADH/Redox (L-Glutamine supports NAD⁺ regeneration)
  // ──────────────────────────────────────────────────────
  dxdt_NADH = kNADH_prod * (1 + E_LG_NADH) - kNADH_ox * NADH;

  // ──────────────────────────────────────────────────────
  // Iron/Ferritin (accumulates from chronic hemolysis)
  // ──────────────────────────────────────────────────────
  double iron_accum = kIron_acc * hemolysis_rate * 100;
  dxdt_Iron = iron_accum - kIron_deg * Iron;

  // ──────────────────────────────────────────────────────
  // TRV (Tricuspid Regurgitant Velocity — pulm HTN proxy)
  // ──────────────────────────────────────────────────────
  double TRV_drive = 1 + 0.002 * (1 - NO / NO_0) + 0.001 * (LDH / LDH_0 - 1);
  dxdt_TRV = kTRV_prod * TRV_drive - kTRV_deg * TRV;

  // ──────────────────────────────────────────────────────
  // eGFR (slow CKD progression driven by sickle nephropathy)
  // ──────────────────────────────────────────────────────
  double GFR_loss = keGFR_dec * (1 + 0.5 * (RBC_S / RBC_S_0 - 1));
  dxdt_eGFR = - GFR_loss * eGFR;

$TABLE
  double VOC_annual = VOC * 8760;   // Convert h⁻¹ → per year
  double HbF_pct    = HbF_frac * 100;
  double Hb_resp    = Hgb;
  double Psel_rel   = P_selectin;
  double ret_pct    = RET / (RBC_S + RBC_N + RET) * 100;
  double freeHb_mgL = free_Hb * 0.64;  // μM → mg/dL approx
  double CRIZ_conc  = CRIZ_C * 148000 / 1e6;  // nM → μg/mL

  capture(VOC_annual, HbF_pct, Hb_resp, Psel_rel, ret_pct,
          freeHb_mgL, CRIZ_conc, LDH, Bilirubin, NO, NADH,
          Iron, TRV, eGFR, Haptoglobin, HU_plasma,
          VOX_RBC, CRIZ_C)
'

## ── Compile model ─────────────────────────────────────────────
mod <- mcode("SCD_QSP", scd_model_code)
mod <- param(mod, DOSE_HU=0, DOSE_VOX=0, DOSE_CRIZ=0, DOSE_LG=0)

## ── Simulation time: 2 years (17520 h) ───────────────────────
sim_times <- c(0, seq(24, 17520, by=24))  # Daily samples

## ── Steady-state dosing events ────────────────────────────────

make_dosing <- function(hu=FALSE, vox=FALSE, criz=FALSE, lg=FALSE,
                        bw=65, sim_days=730) {
  ev_list <- list()

  if (hu) {
    # HU: oral daily, 20 mg/kg/day
    # Handled via DOSE_HU parameter (rate-based input)
  }
  if (vox) {
    # VOX: oral daily (rate-based)
  }
  if (criz) {
    # CRIZ IV q4w: bolus events into CRIZ_C compartment
    # 5 mg/kg → 325 mg → ~2.2 μmol → into Vc=3.5 L → ~630 nM at t=0
    criz_dose_nmol <- (5 * bw * 1e6) / 148000  # nmol
    criz_times <- seq(0, sim_days * 24, by = 28*24)
    ev_list[["criz"]] <- ev(amt = criz_dose_nmol, cmt = "CRIZ_C",
                            time = criz_times, rate = -2)
  }
  do.call(c, ev_list[sapply(ev_list, function(x) !is.null(x))])
}

## ── Treatment scenarios ───────────────────────────────────────

scenarios <- list(
  `1. No Treatment (Baseline)` = list(
    DOSE_HU=0, DOSE_VOX=0, DOSE_CRIZ=0, DOSE_LG=0, criz=FALSE),
  `2. Hydroxyurea (20 mg/kg/d)` = list(
    DOSE_HU=20, DOSE_VOX=0, DOSE_CRIZ=0, DOSE_LG=0, criz=FALSE),
  `3. Voxelotor (1500 mg/d)` = list(
    DOSE_HU=0, DOSE_VOX=1, DOSE_CRIZ=0, DOSE_LG=0, criz=FALSE),
  `4. Crizanlizumab (5 mg/kg q4w)` = list(
    DOSE_HU=0, DOSE_VOX=0, DOSE_CRIZ=1, DOSE_LG=0, criz=TRUE),
  `5. L-Glutamine (5g BID)` = list(
    DOSE_HU=0, DOSE_VOX=0, DOSE_CRIZ=0, DOSE_LG=1, criz=FALSE),
  `6. HU + Voxelotor` = list(
    DOSE_HU=20, DOSE_VOX=1, DOSE_CRIZ=0, DOSE_LG=0, criz=FALSE),
  `7. HU + VOX + CRIZ (Triple)` = list(
    DOSE_HU=20, DOSE_VOX=1, DOSE_CRIZ=1, DOSE_LG=0, criz=TRUE)
)

## ── Run all scenarios ─────────────────────────────────────────

run_scenario <- function(mod, scen_params, bw=65) {
  p <- modifyList(list(DOSE_HU=0, DOSE_VOX=0, DOSE_CRIZ=0,
                        DOSE_LG=0, BW=bw), scen_params)
  criz_on <- isTRUE(p$criz)

  mod_p <- param(mod, DOSE_HU=p$DOSE_HU, DOSE_VOX=p$DOSE_VOX,
                 DOSE_CRIZ=p$DOSE_CRIZ, DOSE_LG=p$DOSE_LG, BW=bw)

  if (criz_on) {
    ev_criz <- make_dosing(criz=TRUE, bw=bw, sim_days=730)
    out <- mrgsim(mod_p, events=ev_criz, times=sim_times,
                  delta=0.5, end=17520)
  } else {
    out <- mrgsim(mod_p, times=sim_times, delta=0.5, end=17520)
  }
  as.data.frame(out)
}

results <- lapply(names(scenarios), function(nm) {
  cat("Running:", nm, "\n")
  df <- run_scenario(mod, scenarios[[nm]])
  df$scenario <- nm
  df
})

results_df <- bind_rows(results) %>%
  mutate(time_days = time / 24,
         time_weeks = time / 168)

## ── Plotting functions ────────────────────────────────────────

theme_scd <- function() {
  theme_bw(base_size=13) +
    theme(
      legend.position="bottom",
      legend.key.size=unit(0.5,"cm"),
      plot.title=element_text(face="bold", size=14),
      strip.background=element_rect(fill="#3A5F8A"),
      strip.text=element_text(colour="white", face="bold")
    )
}

scenario_colors <- c(
  "1. No Treatment (Baseline)"       = "#7F8C8D",
  "2. Hydroxyurea (20 mg/kg/d)"      = "#2980B9",
  "3. Voxelotor (1500 mg/d)"         = "#27AE60",
  "4. Crizanlizumab (5 mg/kg q4w)"   = "#E67E22",
  "5. L-Glutamine (5g BID)"          = "#9B59B6",
  "6. HU + Voxelotor"                = "#C0392B",
  "7. HU + VOX + CRIZ (Triple)"      = "#1ABC9C"
)

## Plot 1: Hemoglobin response
p_hgb <- ggplot(results_df %>% filter(time_days %% 7 < 1),
                aes(time_days, Hb_resp, colour=scenario)) +
  geom_line(linewidth=0.9) +
  scale_colour_manual(values=scenario_colors, name="") +
  labs(title="Hemoglobin Response Over 2 Years",
       x="Time (days)", y="Hemoglobin (g/dL)") +
  geom_hline(yintercept=9, linetype="dashed", colour="grey40") +
  annotate("text", x=50, y=9.15, label="Target Hgb ≥9 g/dL",
           size=3, colour="grey40") +
  theme_scd()

## Plot 2: HbF Fraction
p_hbf <- ggplot(results_df %>% filter(time_days %% 7 < 1),
                aes(time_days, HbF_pct, colour=scenario)) +
  geom_line(linewidth=0.9) +
  scale_colour_manual(values=scenario_colors, name="") +
  labs(title="HbF Percentage Over 2 Years",
       x="Time (days)", y="HbF (%)") +
  geom_hline(yintercept=20, linetype="dashed", colour="#27AE60") +
  annotate("text", x=50, y=20.5, label="Target HbF ≥20%",
           size=3, colour="#27AE60") +
  theme_scd()

## Plot 3: VOC Annual Rate
p_voc <- ggplot(results_df %>% filter(time_days %% 7 < 1),
                aes(time_days, VOC_annual, colour=scenario)) +
  geom_line(linewidth=0.9) +
  scale_colour_manual(values=scenario_colors, name="") +
  labs(title="Estimated Annual VOC Rate",
       x="Time (days)", y="VOC Rate (crises/year)") +
  theme_scd()

## Plot 4: LDH (hemolysis marker)
p_ldh <- ggplot(results_df %>% filter(time_days %% 7 < 1),
                aes(time_days, LDH, colour=scenario)) +
  geom_line(linewidth=0.9) +
  scale_colour_manual(values=scenario_colors, name="") +
  labs(title="LDH (Hemolysis Marker)",
       x="Time (days)", y="LDH (U/L)") +
  geom_hline(yintercept=250, linetype="dashed", colour="grey40") +
  annotate("text", x=50, y=260, label="ULN = 250 U/L",
           size=3, colour="grey40") +
  theme_scd()

## Plot 5: NO Bioavailability
p_no <- ggplot(results_df %>% filter(time_days %% 7 < 1),
               aes(time_days, NO, colour=scenario)) +
  geom_line(linewidth=0.9) +
  scale_colour_manual(values=scenario_colors, name="") +
  labs(title="Nitric Oxide Bioavailability",
       x="Time (days)", y="NO Index (a.u.)") +
  theme_scd()

## Plot 6: P-Selectin (adhesion)
p_psel <- ggplot(results_df %>% filter(time_days %% 7 < 1),
                 aes(time_days, Psel_rel, colour=scenario)) +
  geom_line(linewidth=0.9) +
  scale_colour_manual(values=scenario_colors, name="") +
  labs(title="P-Selectin Expression (Relative to Baseline SCD)",
       x="Time (days)", y="P-Selectin (relative units)") +
  theme_scd()

## Plot 7: TRV (Pulmonary HTN)
p_trv <- ggplot(results_df %>% filter(time_days %% 7 < 1),
                aes(time_days, TRV, colour=scenario)) +
  geom_line(linewidth=0.9) +
  scale_colour_manual(values=scenario_colors, name="") +
  labs(title="Tricuspid Regurgitant Velocity (TRV)",
       x="Time (days)", y="TRV (m/s)") +
  geom_hline(yintercept=2.5, linetype="dashed", colour="#E74C3C") +
  annotate("text", x=50, y=2.52, label="PH threshold: >2.5 m/s",
           size=3, colour="#E74C3C") +
  theme_scd()

## Plot 8: eGFR trend
p_gfr <- ggplot(results_df %>% filter(time_days %% 7 < 1),
                aes(time_days, eGFR, colour=scenario)) +
  geom_line(linewidth=0.9) +
  scale_colour_manual(values=scenario_colors, name="") +
  labs(title="eGFR Trajectory (CKD Progression)",
       x="Time (days)", y="eGFR (mL/min/1.73m²)") +
  theme_scd()

## ── Summary table at 6 months ─────────────────────────────────
summary_6mo <- results_df %>%
  filter(time_days >= 179, time_days <= 181) %>%
  group_by(scenario) %>%
  summarise(
    `Hgb (g/dL)`      = round(mean(Hb_resp), 2),
    `HbF (%)`         = round(mean(HbF_pct), 1),
    `VOC (crises/yr)` = round(mean(VOC_annual), 2),
    `LDH (U/L)`       = round(mean(LDH), 0),
    `P-Selectin`      = round(mean(Psel_rel), 2),
    `NO Index`        = round(mean(NO), 3),
    `TRV (m/s)`       = round(mean(TRV), 3),
    `eGFR`            = round(mean(eGFR), 1),
    .groups="drop"
  )

## ── Summary table at 1 year ───────────────────────────────────
summary_1yr <- results_df %>%
  filter(time_days >= 364, time_days <= 366) %>%
  group_by(scenario) %>%
  summarise(
    `Hgb (g/dL)`      = round(mean(Hb_resp), 2),
    `ΔHgb vs baseline`= round(mean(Hb_resp) - 8.5, 2),
    `HbF (%)`         = round(mean(HbF_pct), 1),
    `VOC (crises/yr)` = round(mean(VOC_annual), 2),
    `VOC reduction %` = round((1 - mean(VOC_annual)/
                                 filter(results_df, scenario=="1. No Treatment (Baseline)",
                                        time_days>=364, time_days<=366)$VOC_annual[1]) * 100, 1),
    `LDH (U/L)`       = round(mean(LDH), 0),
    `eGFR`            = round(mean(eGFR), 1),
    .groups="drop"
  )

## ── Print results ─────────────────────────────────────────────
cat("\n=== SCD QSP Model: 6-Month Summary ===\n")
print(summary_6mo, n=Inf)
cat("\n=== SCD QSP Model: 1-Year Summary ===\n")
print(summary_1yr, n=Inf)

## ── Sensitivity analysis: HU dose-response ───────────────────
hu_doses <- c(5, 10, 15, 20, 25, 30, 35)  # mg/kg/day

sa_results <- lapply(hu_doses, function(d) {
  p <- list(DOSE_HU=d, DOSE_VOX=0, DOSE_CRIZ=0, DOSE_LG=0, criz=FALSE)
  df <- run_scenario(mod, p)
  df$HU_dose <- paste0(d, " mg/kg/d")
  df$dose_num <- d
  df
})
sa_df <- bind_rows(sa_results) %>%
  mutate(time_days = time / 24)

p_sa_hbf <- ggplot(sa_df %>% filter(time_days %% 7 < 1),
                   aes(time_days, HbF_pct, colour=factor(dose_num))) +
  geom_line(linewidth=0.9) +
  scale_colour_brewer(palette="Blues", name="HU Dose (mg/kg/d)") +
  labs(title="Hydroxyurea Dose-Response: HbF Induction",
       x="Time (days)", y="HbF (%)") +
  theme_scd()

p_sa_hgb <- ggplot(sa_df %>% filter(time_days %% 7 < 1),
                   aes(time_days, Hb_resp, colour=factor(dose_num))) +
  geom_line(linewidth=0.9) +
  scale_colour_brewer(palette="Blues", name="HU Dose (mg/kg/d)") +
  labs(title="Hydroxyurea Dose-Response: Hemoglobin",
       x="Time (days)", y="Hemoglobin (g/dL)") +
  theme_scd()

cat("\nModel compilation and simulation complete.\n")
cat("Key outputs: p_hgb, p_hbf, p_voc, p_ldh, p_no, p_psel, p_trv, p_gfr\n")
cat("Sensitivity: p_sa_hbf, p_sa_hgb\n")
cat("Summary tables: summary_6mo, summary_1yr\n")
