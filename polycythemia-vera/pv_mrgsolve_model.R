## ============================================================
## Polycythemia Vera (PV) — mrgsolve QSP ODE Model
## ============================================================
## Disease: Polycythemia Vera (WHO 2022: JAK2V617F / exon12 mutation
##          + HCT ≥45% M / ≥42% F, or elevated RBC mass)
## Key pathways: JAK2V617F → JAK-STAT5 → EPO-independent erythropoiesis
##               Thrombosis, splenomegaly, MF progression
## Drug classes: Phlebotomy · Hydroxyurea · Ruxolitinib ·
##               Ropeginterferon-α2b · Low-dose aspirin ·
##               Fedratinib · Anagrelide
## Calibration references (key trials):
##   CYTOREDUCE (phlebotomy vs HCT<45%)
##   RESPONSE/RESPONSE-2 (ruxolitinib vs BAT, NEJM 2015/Blood 2017)
##   PROUD-PV / CONTINUATION-PV (ropeg-IFN, Blood 2019/Leukemia 2020)
##   MAJIC-PV (ruxolitinib vs HU, Lancet Haematol 2017)
##   EV-Rux (fedratinib, ASH 2021)
## Author: Claude Code Routine — Catholic Univ. Seoul, 2026-06-24
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ============================================================
## Model code block
## ============================================================
pv_model_code <- '
$PROB
Polycythemia Vera QSP Model
Compartments: 22 ODE states
Treatments: Phlebotomy · HU · Ruxolitinib · Ropeg-IFN · Aspirin · Fedratinib · Anagrelide

$PARAM @annotated
// ---- Clonal Hematopoiesis ----
k_clone_grow  : 0.012  : Mutant HSC net growth rate (/day)
k_clone_death : 0.003  : Mutant HSC basal death rate (/day)
k_wt_grow     : 0.008  : Wild-type HSC growth rate (/day)
k_wt_death    : 0.005  : Wild-type HSC death rate (/day)
tau_allele    : 30     : Time constant for allele burden equilibration (days)

// ---- JAK-STAT Signaling ----
kin_stat5     : 0.80   : STAT5 activation rate constant (/day)
kout_stat5    : 0.75   : STAT5 deactivation rate (/day, t1/2~0.9d)
STAT5_base    : 1.0    : Baseline STAT5 activity (normalized)

// ---- Erythropoiesis ----
k_bfu_prod    : 18.0   : Mutant BFU-E production rate (cells/µL/day, per STAT5 unit)
k_bfu_prod_wt : 6.0    : Wild-type BFU-E production (/day)
k_bfu_mat     : 0.15   : BFU-E maturation rate to CFU-E (/day)
k_rbc_prod    : 0.028  : RBC production from progenitors (/day)
k_rbc_dest    : 0.00833: RBC destruction rate (/day, t1/2=120d)
RBC_base      : 5.5    : Baseline RBC (×10^12/L)
HCT_per_RBC   : 0.082  : HCT fraction per RBC unit (L/L per 10^12/L)
phlebotomy_eff : 0.045 : RBC removed per phlebotomy session (×10^12/L)

// ---- Platelet / WBC ----
k_plt_prod    : 2.5    : Platelet production (×10^9/L/day per STAT5 unit)
k_plt_dest    : 0.10   : Platelet destruction (/day, t1/2=10d)
PLT_base      : 350    : Baseline PLT (×10^9/L)
k_wbc_prod    : 0.30   : WBC production (×10^9/L/day per STAT5 unit)
k_wbc_dest    : 0.14   : WBC destruction (/day, t1/2=7d)
WBC_base      : 8.5    : Baseline WBC (×10^9/L)

// ---- Spleen & EMH ----
k_spl_grow    : 0.0015 : Spleen growth rate (/day per STAT5 unit above base)
k_spl_base    : 0.0004 : Baseline spleen volume regulation (/day)
Spleen_norm   : 400    : Normal spleen volume (mL)

// ---- MPN Symptom Score ----
k_tss_rise    : 0.04   : MPN-SAF-TSS increase rate per STAT5 unit (/day)
k_tss_fall    : 0.035  : MPN-SAF-TSS natural remission (/day)
TSS_base      : 8.0    : Baseline symptom score (0–100)

// ---- Fibrosis Progression ----
k_mf          : 0.00008: Marrow fibrosis progression rate (/day)
TGFb_driver   : 1.5    : TGF-β fibrosis driver (proportional to STAT5)
MF_max        : 3.0    : Maximum fibrosis grade

// ---- Thrombosis Hazard ----
k_thromb_base : 0.005  : Baseline thrombosis hazard rate (/day)
HCT_thr_threshold : 45.0: HCT threshold for thrombosis risk increase (%)
k_thromb_HCT  : 0.0008 : HCT excess → thrombosis hazard (/day per % above threshold)
k_thromb_PLT  : 0.0001 : PLT → thrombosis hazard (/day per 100×10^9/L above 400)

// ---- Ruxolitinib PK (2-compartment oral) ----
// RESPONSE trial: 10mg BID → Cmax ~1200 ng/mL, t1/2 ~3h
// MAJIC-PV: 10–25mg BID dosing
ka_rux        : 2.16   : Ruxolitinib absorption rate (/h)
CL_rux        : 17.7   : Ruxolitinib clearance (L/h)
Vd1_rux       : 50.0   : Central volume (L)
Vd2_rux       : 22.0   : Peripheral volume (L)
CLD_rux       : 12.0   : Inter-compartmental clearance (L/h)
Emax_rux      : 0.92   : Maximum JAK2 inhibition (%)
IC50_rux      : 450    : IC50 ruxolitinib JAK2 (ng/mL)
Hill_rux      : 1.2    : Hill coefficient (ruxolitinib)

// ---- Hydroxyurea PK (1-compartment oral) ----
// t1/2 ~3-4h; Cmax ~300-700 µmol/L at 1500mg dose
ka_hu         : 2.88   : HU absorption rate (/h)
CL_hu         : 18.5   : HU clearance (L/h)
Vd_hu         : 50.0   : HU distribution volume (L)
Emax_hu_rbc   : 0.65   : Maximum RBC suppression (HU)
Emax_hu_plt   : 0.70   : Maximum PLT suppression (HU)
Emax_hu_wbc   : 0.75   : Maximum WBC suppression (HU)
EC50_hu       : 8.5    : HU EC50 (mg/L, ~100 µmol/L)
Hill_hu       : 1.1    : Hill coefficient (HU)
MW_hu         : 76.05  : Molecular weight HU (g/mol)

// ---- Ropeginterferon-α2b PK (SC, 1-compartment) ----
// PROUD-PV: 100mcg q2w → 250mcg/mL Cmax; t1/2 80-130h
ka_ifn        : 0.15   : IFN absorption rate from SC (/day)
CL_ifn        : 1.2    : IFN clearance (L/day)
Vd_ifn        : 12.0   : IFN volume of distribution (L)
Emax_ifn_clone : 0.90  : Maximum clone suppression by IFN
EC50_ifn      : 250    : IFN EC50 for clone suppression (IU/mL)
Hill_ifn      : 1.4    : Hill coefficient (IFN)

// ---- Low-Dose Aspirin PK ----
ka_asp        : 4.32   : Aspirin absorption (/h)
CL_asp        : 25.0   : Aspirin clearance (L/h)
Vd_asp        : 8.0    : Aspirin volume (L)
Emax_asp_plt  : 0.85   : Max platelet activation inhibition
EC50_asp      : 0.05   : Aspirin EC50 (mg/L, irreversible COX-1)

// ---- Fedratinib PK (1-compartment oral) ----
// JAKARTA-2: 400mg/day → Cmax ~2500 ng/mL, t1/2 ~41h
ka_fed        : 0.48   : Fedratinib absorption (/h)
CL_fed        : 4.2    : Fedratinib clearance (L/h)
Vd_fed        : 212    : Fedratinib volume (L)
Emax_fed      : 0.88   : Max JAK2 inhibition (fedratinib)
IC50_fed      : 600    : Fedratinib IC50 (ng/mL)
Hill_fed      : 1.2    : Hill coefficient (fedratinib)

// ---- Anagrelide PK (1-compartment oral) ----
ka_ana        : 2.4    : Anagrelide absorption (/h)
CL_ana        : 21.0   : Anagrelide clearance (L/h)
Vd_ana        : 45.0   : Anagrelide volume (L)
Emax_ana_plt  : 0.80   : Max platelet reduction (anagrelide)
EC50_ana      : 2.0    : Anagrelide EC50 (ng/mL)

// ---- Dosing flags (0=off, 1=on) ----
use_phlebotomy  : 0    : Phlebotomy (1=active)
phlebotomy_freq : 0    : Phlebotomy sessions per year
use_HU          : 0    : Hydroxyurea (1=active)
HU_daily_dose   : 0    : HU daily dose (mg)
use_RUX         : 0    : Ruxolitinib (1=active)
RUX_BID_dose    : 0    : Ruxolitinib BID dose (mg)
use_IFN         : 0    : Ropeg-IFN (1=active)
IFN_q2w_dose    : 0    : Ropeg-IFN q2w dose (mcg)
use_ASP         : 0    : Aspirin (1=active)
ASP_daily_dose  : 0    : Aspirin dose (mg/day)
use_FED         : 0    : Fedratinib (1=active)
FED_daily_dose  : 0    : Fedratinib daily dose (mg)
use_ANA         : 0    : Anagrelide (1=active)
ANA_daily_dose  : 0    : Anagrelide daily dose (mg)

$INIT @annotated
// ---- Clonal State (normalized 0–1) ----
mut_clone   = 0.30  : Mutant HSC clone fraction (initial 30% = early disease)
wt_clone    = 0.70  : Wild-type HSC fraction

// ---- Signaling ----
STAT5       = 1.40  : Active STAT5 (elevated in PV)
allele_burden = 45.0 : JAK2 allele burden (%, initial)

// ---- Erythroid Progenitors ----
BFU_E_mut   = 80.0  : Mutant BFU-E (cells/µL)
BFU_E_wt    = 40.0  : Wild-type BFU-E (cells/µL)

// ---- Mature Blood Cells ----
RBC         = 6.8   : Total RBC (×10^12/L, elevated in PV)
HCT         = 55.0  : Hematocrit (%, elevated)
PLT         = 550   : Platelet count (×10^9/L, elevated)
WBC         = 12.5  : WBC count (×10^9/L, elevated)

// ---- Organ & Symptom State ----
Spleen_vol  = 850   : Spleen volume (mL, enlarged in PV)
MPN_SAF     = 20.0  : MPN-SAF TSS (0–100, symptomatic baseline)
MF_score    = 0.30  : Fibrosis score (MF-0~MF-1 at baseline)

// ---- Risk/Events ----
Thromb_hazard = 0.0 : Cumulative thrombosis hazard (dimensionless)

// ---- Drug PK States ----
RUX_gut     = 0.0   : Ruxolitinib in GI (mg)
RUX_cent    = 0.0   : Ruxolitinib central (ng/mL)
RUX_periph  = 0.0   : Ruxolitinib peripheral (ng/mL)
HU_cent     = 0.0   : Hydroxyurea central (mg/L)
IFN_sc      = 0.0   : IFN-α SC depot (IU)
IFN_cent    = 0.0   : IFN-α central (IU/mL)
ASP_cent    = 0.0   : Aspirin central (mg/L)
FED_cent    = 0.0   : Fedratinib central (ng/mL)
ANA_cent    = 0.0   : Anagrelide central (ng/mL)

$MAIN
// ---- PD effect calculations ----
// Ruxolitinib JAK2 inhibition (Emax model)
double RUX_C_ngmL = RUX_cent;
double rux_inhib = use_RUX * Emax_rux * pow(RUX_C_ngmL, Hill_rux) /
                   (pow(IC50_rux, Hill_rux) + pow(RUX_C_ngmL, Hill_rux));

// Fedratinib JAK2 inhibition
double FED_C_ngmL = FED_cent;
double fed_inhib = use_FED * Emax_fed * pow(FED_C_ngmL, Hill_fed) /
                   (pow(IC50_fed, Hill_fed) + pow(FED_C_ngmL, Hill_fed));

// Combined JAK2 inhibition (max of independent inhibitors, simplified)
double jak2_inhib = 1.0 - (1.0 - rux_inhib) * (1.0 - fed_inhib);

// Hydroxyurea cytoreduction
double HU_C = HU_cent;
double hu_eff_rbc = use_HU * Emax_hu_rbc * pow(HU_C, Hill_hu) /
                    (pow(EC50_hu, Hill_hu) + pow(HU_C, Hill_hu));
double hu_eff_plt = use_HU * Emax_hu_plt * pow(HU_C, Hill_hu) /
                    (pow(EC50_hu, Hill_hu) + pow(HU_C, Hill_hu));
double hu_eff_wbc = use_HU * Emax_hu_wbc * pow(HU_C, Hill_hu) /
                    (pow(EC50_hu, Hill_hu) + pow(HU_C, Hill_hu));

// Ropeginterferon clone suppression
double IFN_C = IFN_cent;
double ifn_clone_eff = use_IFN * Emax_ifn_clone * pow(IFN_C, Hill_ifn) /
                       (pow(EC50_ifn, Hill_ifn) + pow(IFN_C, Hill_ifn));

// Aspirin platelet inhibition
double asp_plt_inhib = use_ASP * Emax_asp_plt * ASP_cent /
                       (EC50_asp + ASP_cent);

// Anagrelide platelet reduction
double ana_plt_eff = use_ANA * Emax_ana_plt * ANA_cent /
                     (EC50_ana + ANA_cent);

// ---- Phlebotomy effect (pulse, modeled as continuous drain) ----
// Approx: phlebotomy_freq sessions/year × 0.045 ×10^12/L per session
double phlebotomy_drain = use_phlebotomy * phlebotomy_freq * phlebotomy_eff / 365.0;

// ---- Derived disease quantities ----
double STAT5_eff = STAT5 * (1.0 - jak2_inhib);  // effective STAT5 after JAK2 inhibition
double HCT_excess = (HCT > HCT_thr_threshold) ? (HCT - HCT_thr_threshold) : 0.0;
double PLT_excess = (PLT > 400.0) ? (PLT - 400.0) / 100.0 : 0.0;

$ODE
// ============================================================
// 1. CLONAL HEMATOPOIESIS
// ============================================================
// Mutant clone logistic growth with IFN-mediated elimination
dxdt_mut_clone = k_clone_grow * mut_clone * (1.0 - mut_clone - wt_clone)
               - k_clone_death * mut_clone
               - ifn_clone_eff * mut_clone;

// Wild-type HSC (stabilized by loss of mutant competition)
dxdt_wt_clone = k_wt_grow * wt_clone * (1.0 - mut_clone - wt_clone)
              - k_wt_death * wt_clone;

// Allele burden equilibrates to mut_clone fraction (×100 for %)
dxdt_allele_burden = (mut_clone / (mut_clone + wt_clone + 1e-6) * 100.0 - allele_burden) / tau_allele;

// ============================================================
// 2. JAK-STAT5 SIGNALING
// ============================================================
// STAT5 driven by allele burden (proportional), inhibited by JAK inhibitors
dxdt_STAT5 = kin_stat5 * (allele_burden / 50.0) * (1.0 - jak2_inhib)
           - kout_stat5 * STAT5;

// ============================================================
// 3. ERYTHROID PROGENITORS
// ============================================================
// Mutant BFU-E: EPO-independent, driven by STAT5
dxdt_BFU_E_mut = k_bfu_prod * STAT5_eff * (1.0 - hu_eff_rbc)
               - k_bfu_mat * BFU_E_mut
               - ifn_clone_eff * BFU_E_mut;

// Wild-type BFU-E: suppressed by HU
dxdt_BFU_E_wt = k_bfu_prod_wt * (1.0 - hu_eff_rbc)
              - k_bfu_mat * BFU_E_wt;

// ============================================================
// 4. MATURE BLOOD CELLS
// ============================================================
// RBC: produced from progenitors, destroyed at 1/120d, removed by phlebotomy
dxdt_RBC = k_rbc_prod * (BFU_E_mut + BFU_E_wt)
          - k_rbc_dest * RBC
          - phlebotomy_drain;

// HCT: proportional to RBC (equilibration)
dxdt_HCT = (RBC * HCT_per_RBC * 100.0 - HCT) / 5.0;  // equilibrate over ~5 days

// Platelets: STAT5-driven production, removed by HU and anagrelide
dxdt_PLT = k_plt_prod * STAT5_eff * (1.0 - hu_eff_plt) * (1.0 - ana_plt_eff)
          - k_plt_dest * PLT;

// WBC (neutrophils dominant): STAT5-driven, removed by HU
dxdt_WBC = k_wbc_prod * STAT5_eff * (1.0 - hu_eff_wbc)
          - k_wbc_dest * WBC;

// ============================================================
// 5. SPLEEN VOLUME (mL)
// ============================================================
// Grows when STAT5-driven EMH exceeds normal; ruxolitinib suppresses EMH
dxdt_Spleen_vol = k_spl_grow * (STAT5_eff - STAT5_base) * Spleen_vol * (1.0 - jak2_inhib)
                - k_spl_base * (Spleen_vol - Spleen_norm);

// ============================================================
// 6. MPN SYMPTOM SCORE (MPN-SAF-TSS, 0–100)
// ============================================================
dxdt_MPN_SAF = k_tss_rise * (STAT5_eff - STAT5_base) * 20.0
             - k_tss_fall * MPN_SAF * (1.0 + jak2_inhib * 1.5);
// Clamp to 0–100 handled by ALAG/scale below conceptually

// ============================================================
// 7. BONE MARROW FIBROSIS SCORE (0–3)
// ============================================================
dxdt_MF_score = k_mf * TGFb_driver * STAT5_eff * (MF_max - MF_score);

// ============================================================
// 8. THROMBOSIS HAZARD (cumulative, dimensionless)
// ============================================================
dxdt_Thromb_hazard = k_thromb_base
                   + k_thromb_HCT * HCT_excess
                   + k_thromb_PLT * PLT_excess
                   - asp_plt_inhib * k_thromb_base * 0.40;

// ============================================================
// 9. RUXOLITINIB PK (2-compartment oral)
// ============================================================
double RUX_dose_h = use_RUX * RUX_BID_dose;  // mg per dose (BID handled by events)
dxdt_RUX_gut   = -ka_rux * RUX_gut;  // units: mg
dxdt_RUX_cent  = ka_rux * RUX_gut / Vd1_rux * 1000.0  // convert mg/L to ng/mL (×1000)
               - (CL_rux + CLD_rux) / Vd1_rux * RUX_cent
               + CLD_rux / Vd1_rux * RUX_periph;
dxdt_RUX_periph = CLD_rux / Vd1_rux * RUX_cent
                - CLD_rux / Vd2_rux * RUX_periph;

// ============================================================
// 10. HYDROXYUREA PK (1-compartment)
// ============================================================
double HU_dose_h = use_HU * HU_daily_dose / 24.0;  // mg/h
dxdt_HU_cent   = HU_dose_h / Vd_hu         // input (mg/L/h)
               - CL_hu / Vd_hu * HU_cent;  // elimination

// ============================================================
// 11. ROPEGINTERFERON α-2b PK (SC → central)
// ============================================================
// IFN_sc in IU; dosed q2w (handled by events in IU)
dxdt_IFN_sc    = -ka_ifn * IFN_sc;
dxdt_IFN_cent  = ka_ifn * IFN_sc / Vd_ifn
               - CL_ifn / Vd_ifn * IFN_cent;

// ============================================================
// 12. ASPIRIN PK (1-compartment)
// ============================================================
double ASP_dose_h = use_ASP * ASP_daily_dose / 24.0;  // mg/h
dxdt_ASP_cent  = ASP_dose_h / Vd_asp
               - CL_asp / Vd_asp * ASP_cent;

// ============================================================
// 13. FEDRATINIB PK (1-compartment)
// ============================================================
double FED_dose_h = use_FED * FED_daily_dose / 24.0;  // mg/h
dxdt_FED_cent  = FED_dose_h / Vd_fed * 1000.0  // mg/h → ng/mL/h
               - CL_fed / Vd_fed * FED_cent;

// ============================================================
// 14. ANAGRELIDE PK (1-compartment)
// ============================================================
double ANA_dose_h = use_ANA * ANA_daily_dose / 24.0;  // mg/h
dxdt_ANA_cent  = ANA_dose_h / Vd_ana * 1e6  // mg/h → ng/mL/h
               - CL_ana / Vd_ana * ANA_cent;

$TABLE
// Derived outputs for plotting
double HCT_pct       = HCT;
double PLT_out       = PLT;
double WBC_out       = WBC;
double allele_out    = allele_burden;
double spleen_out    = Spleen_vol;
double tss_out       = (MPN_SAF < 0) ? 0 : ((MPN_SAF > 100) ? 100 : MPN_SAF);
double mf_out        = MF_score;
double thromb_out    = Thromb_hazard;
double rux_cp        = RUX_cent;
double hu_cp         = HU_cent;
double ifn_cp        = IFN_cent;
double jak2_inhib_out = 1.0 - (1.0 - rux_inhib) * (1.0 - fed_inhib);
double svr35         = (Spleen_vol <= 260) ? 1 : 0;  // SVR35 from 400mL baseline
double hct_ctrl      = (HCT < 45.0) ? 1 : 0;         // HCT control endpoint
double tss50         = (tss_out <= TSS_base * 0.5) ? 1 : 0;  // TSS50 response

$CAPTURE
HCT_pct allele_out PLT_out WBC_out spleen_out tss_out mf_out thromb_out
rux_cp hu_cp ifn_cp jak2_inhib_out svr35 hct_ctrl tss50 RBC
'

## ============================================================
## Compile model
## ============================================================
pv_mod <- mcode("pv_qsp", pv_model_code)
cat("PV QSP model compiled successfully.\n")

## ============================================================
## Treatment Scenarios
## ============================================================
## S0: Natural history (no treatment, 5 years)
## S1: Phlebotomy + Low-dose aspirin (low-risk management)
## S2: Hydroxyurea 1500mg/day + Phlebotomy + Aspirin (standard HU)
## S3: Ruxolitinib 10mg BID (HU-resistant/intolerant, RESPONSE)
## S4: Ruxolitinib 20mg BID (high-risk, aggressive)
## S5: Ropeginterferon-α2b 100mcg q2w (molecular remission aim)
## S6: Fedratinib 400mg/day (JAK2/FLT3, HU-failed)

sim_time_days <- 0:(5 * 365)  # 5-year simulation

## ---- Helper: Ruxolitinib BID dosing event table ----
make_rux_events <- function(dose_mg, start_day = 1, end_day = 5*365) {
  ev <- ev(
    time  = c(seq(start_day, end_day, by = 1),
              seq(start_day + 0.5, end_day + 0.5, by = 1)),  # BID ~12h apart
    cmt   = "RUX_gut",
    amt   = dose_mg,
    evid  = 1
  )
  return(ev)
}

## ---- Helper: IFN q2w dosing ----
make_ifn_events <- function(dose_mcg, start_day = 1, end_day = 5*365) {
  # Convert mcg to IU: ~3×10^6 IU per 100mcg for ropeginterferon
  dose_IU <- dose_mcg * 3e4  # simplified conversion
  ev_ifn <- ev(
    time = seq(start_day, end_day, by = 14),  # q2w
    cmt  = "IFN_sc",
    amt  = dose_IU,
    evid = 1
  )
  return(ev_ifn)
}

## ============================================================
## Run all scenarios
## ============================================================

## S0: Natural history
pv_S0 <- pv_mod %>%
  param(use_phlebotomy = 0, use_HU = 0, use_RUX = 0,
        use_IFN = 0, use_ASP = 0, use_FED = 0, use_ANA = 0) %>%
  mrgsim(end = max(sim_time_days), delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "S0: Natural History")

## S1: Phlebotomy + Aspirin
pv_S1 <- pv_mod %>%
  param(use_phlebotomy = 1, phlebotomy_freq = 8,   # ~8 phlebotomies/year initially
        use_ASP = 1, ASP_daily_dose = 100,
        use_HU = 0, use_RUX = 0, use_IFN = 0) %>%
  mrgsim(end = max(sim_time_days), delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "S1: Phlebotomy + Aspirin")

## S2: Hydroxyurea + Phlebotomy + Aspirin (standard cytoreduction)
pv_S2 <- pv_mod %>%
  param(use_phlebotomy = 1, phlebotomy_freq = 3,
        use_HU = 1, HU_daily_dose = 1500,
        use_ASP = 1, ASP_daily_dose = 100,
        use_RUX = 0, use_IFN = 0) %>%
  mrgsim(end = max(sim_time_days), delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "S2: HU + Phlebotomy + Aspirin")

## S3: Ruxolitinib 10mg BID (RESPONSE trial: HU-resistant/intolerant)
ev_rux10 <- make_rux_events(dose_mg = 10)
pv_S3 <- pv_mod %>%
  param(use_RUX = 1, RUX_BID_dose = 10,
        use_phlebotomy = 1, phlebotomy_freq = 2,
        use_ASP = 1, ASP_daily_dose = 100,
        use_HU = 0, use_IFN = 0) %>%
  mrgsim(events = ev_rux10, end = max(sim_time_days), delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "S3: Ruxolitinib 10mg BID")

## S4: Ruxolitinib 20mg BID (high-risk, aggressive disease)
ev_rux20 <- make_rux_events(dose_mg = 20)
pv_S4 <- pv_mod %>%
  param(use_RUX = 1, RUX_BID_dose = 20,
        use_phlebotomy = 1, phlebotomy_freq = 1,
        use_ASP = 1, ASP_daily_dose = 100,
        use_HU = 0, use_IFN = 0) %>%
  mrgsim(events = ev_rux20, end = max(sim_time_days), delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "S4: Ruxolitinib 20mg BID")

## S5: Ropeginterferon-α2b 100mcg q2w → molecular remission
ev_ifn100 <- make_ifn_events(dose_mcg = 100)
pv_S5 <- pv_mod %>%
  param(use_IFN = 1, IFN_q2w_dose = 100,
        use_ASP = 1, ASP_daily_dose = 100,
        use_phlebotomy = 1, phlebotomy_freq = 4,
        use_HU = 0, use_RUX = 0) %>%
  mrgsim(events = ev_ifn100, end = max(sim_time_days), delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "S5: Ropeg-IFN 100mcg q2w")

## S6: Fedratinib 400mg/day (HU + ruxolitinib failed)
pv_S6 <- pv_mod %>%
  param(use_FED = 1, FED_daily_dose = 400,
        use_ASP = 1, ASP_daily_dose = 100,
        use_phlebotomy = 1, phlebotomy_freq = 2,
        use_HU = 0, use_RUX = 0, use_IFN = 0) %>%
  mrgsim(end = max(sim_time_days), delta = 1) %>%
  as.data.frame() %>%
  mutate(scenario = "S6: Fedratinib 400mg/day")

## ============================================================
## Combine and summarize results
## ============================================================
all_results <- bind_rows(pv_S0, pv_S1, pv_S2, pv_S3, pv_S4, pv_S5, pv_S6) %>%
  mutate(time_yr = time / 365)

## Summary at key time points (Week 32 ≈ Day 224, Year 1, Year 3, Year 5)
key_timepoints <- c(0, 112, 224, 365, 730, 1095, 1825)  # days

summary_table <- all_results %>%
  filter(time %in% key_timepoints) %>%
  select(scenario, time, HCT_pct, allele_out, PLT_out, WBC_out,
         spleen_out, tss_out, mf_out, jak2_inhib_out) %>%
  mutate(
    time_label = case_when(
      time == 0   ~ "Baseline",
      time == 112 ~ "Week 16",
      time == 224 ~ "Week 32 (RESPONSE endpoint)",
      time == 365 ~ "Year 1",
      time == 730 ~ "Year 2",
      time == 1095 ~ "Year 3",
      time == 1825 ~ "Year 5"
    )
  )

cat("\n=== PV QSP Model — 5-Year Simulation Summary ===\n")
print(summary_table %>% filter(time %in% c(224, 1095)) %>%
  select(scenario, time_label, HCT_pct, PLT_out, spleen_out, tss_out, allele_out) %>%
  arrange(time_label, scenario))

## ============================================================
## Clinical Response Endpoints at Week 32 (RESPONSE trial design)
## ============================================================
response_wk32 <- all_results %>%
  filter(abs(time - 224) < 1) %>%
  mutate(
    HCT_ctrl         = HCT_pct < 45,
    Spleen_response  = spleen_out < (850 * 0.65),  # SVR35
    TSS_response     = tss_out < (20.0 * 0.50),    # TSS50
    CHR              = HCT_ctrl & PLT_out < 400 & WBC_out < 10
  ) %>%
  select(scenario, HCT_pct, PLT_out, WBC_out, spleen_out, tss_out,
         allele_out, HCT_ctrl, Spleen_response, TSS_response, CHR)

cat("\n=== Response at Week 32 (RESPONSE Trial Primary Endpoint) ===\n")
print(response_wk32 %>%
  select(scenario, HCT_pct, PLT_out, spleen_out, tss_out, allele_out,
         HCT_ctrl, Spleen_response, CHR))

## ============================================================
## Visualization
## ============================================================
scenario_colors <- c(
  "S0: Natural History"           = "#E74C3C",
  "S1: Phlebotomy + Aspirin"     = "#F39C12",
  "S2: HU + Phlebotomy + Aspirin" = "#27AE60",
  "S3: Ruxolitinib 10mg BID"     = "#2980B9",
  "S4: Ruxolitinib 20mg BID"     = "#8E44AD",
  "S5: Ropeg-IFN 100mcg q2w"     = "#16A085",
  "S6: Fedratinib 400mg/day"     = "#D35400"
)

## Plot 1: Hematocrit over 5 years
p_hct <- ggplot(all_results, aes(x = time_yr, y = HCT_pct, color = scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 45, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = 4.8, y = 46.5, label = "Target <45%", color = "red", size = 3.5) +
  scale_color_manual(values = scenario_colors) +
  labs(title    = "Polycythemia Vera — Hematocrit Control",
       subtitle = "Red dashed line: HCT <45% target (CYTOREDUCE trial)",
       x        = "Time (years)",
       y        = "Hematocrit (%)",
       color    = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", legend.key.size = unit(0.5, "cm"))

## Plot 2: JAK2 Allele Burden
p_allele <- ggplot(all_results, aes(x = time_yr, y = allele_out, color = scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 50, linetype = "dotted", color = "#E74C3C", linewidth = 0.8) +
  geom_hline(yintercept = 10, linetype = "dotted", color = "#27AE60", linewidth = 0.8) +
  annotate("text", x = 4.7, y = 52, label = "50% (high burden)", color = "#E74C3C", size = 3.3) +
  annotate("text", x = 4.7, y = 12, label = "10% (mol. remission)", color = "#27AE60", size = 3.3) +
  scale_color_manual(values = scenario_colors) +
  labs(title    = "Polycythemia Vera — JAK2 V617F Allele Burden",
       subtitle = "IFN targets clonal elimination; JAK inhibitors reduce signaling without allele elimination",
       x        = "Time (years)",
       y        = "JAK2 V617F Allele Burden (%)",
       color    = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

## Plot 3: Spleen Volume
p_spleen <- ggplot(all_results, aes(x = time_yr, y = spleen_out, color = scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 850 * 0.65, linetype = "dashed", color = "#2980B9", linewidth = 0.8) +
  annotate("text", x = 4.7, y = 580, label = "SVR35 threshold", color = "#2980B9", size = 3.3) +
  scale_color_manual(values = scenario_colors) +
  labs(title    = "Polycythemia Vera — Spleen Volume",
       subtitle = "SVR35 (≥35% reduction): primary endpoint RESPONSE trial",
       x        = "Time (years)",
       y        = "Spleen Volume (mL)",
       color    = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

## Plot 4: MPN-SAF Symptom Score
p_tss <- ggplot(all_results, aes(x = time_yr, y = tss_out, color = scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "#8E44AD", linewidth = 0.8) +
  annotate("text", x = 4.7, y = 11.5, label = "TSS50 response", color = "#8E44AD", size = 3.3) +
  scale_color_manual(values = scenario_colors) +
  labs(title    = "Polycythemia Vera — MPN-SAF Symptom Score",
       subtitle = "TSS50 response = ≥50% reduction from baseline",
       x        = "Time (years)",
       y        = "MPN-SAF Total Symptom Score (0–100)",
       color    = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

## Plot 5: Platelet count
p_plt <- ggplot(all_results, aes(x = time_yr, y = PLT_out, color = scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 400, linetype = "dashed", color = "#E67E22", linewidth = 0.8) +
  geom_hline(yintercept = 1500, linetype = "dotted", color = "#E74C3C", linewidth = 0.8) +
  annotate("text", x = 4.7, y = 415, label = "Normal upper limit", color = "#E67E22", size = 3.0) +
  annotate("text", x = 4.7, y = 1520, label = "1500 (bleeding risk)", color = "#E74C3C", size = 3.0) +
  scale_color_manual(values = scenario_colors) +
  labs(title    = "Polycythemia Vera — Platelet Count",
       subtitle = "PLT >1500 = acquired von Willebrand disease / bleeding risk",
       x        = "Time (years)",
       y        = "Platelet Count (×10⁹/L)",
       color    = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

## Plot 6: Fibrosis progression
p_mf <- ggplot(all_results, aes(x = time_yr, y = mf_out, color = scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  geom_hline(yintercept = 2.0, linetype = "dotted", color = "#922B21", linewidth = 0.8) +
  annotate("text", x = 4.7, y = 1.1, label = "MF-1", color = "#E74C3C", size = 3.5) +
  annotate("text", x = 4.7, y = 2.1, label = "MF-2", color = "#922B21", size = 3.5) +
  scale_color_manual(values = scenario_colors) +
  labs(title    = "Polycythemia Vera — Marrow Fibrosis Progression",
       subtitle = "Post-PV MF: fibrosis grade ≥2 (Silver stain)",
       x        = "Time (years)",
       y        = "WHO Fibrosis Grade (0–3)",
       color    = "Treatment") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

## Print plots
print(p_hct)
print(p_allele)
print(p_spleen)
print(p_tss)
print(p_plt)
print(p_mf)

## ============================================================
## Ruxolitinib PK Simulation (10mg BID, steady-state)
## ============================================================
rux_pk_mod <- pv_mod %>%
  param(use_RUX = 1, RUX_BID_dose = 10)

ev_rux_ss <- ev(
  time = c(seq(0, 96, by = 24), seq(12, 108, by = 24)),  # BID for 5 days
  cmt  = "RUX_gut",
  amt  = 10,
  evid = 1
) %>% arrange(time)

rux_pk_sim <- rux_pk_mod %>%
  mrgsim(events = ev_rux_ss, end = 120, delta = 0.5) %>%
  as.data.frame()

p_pk_rux <- ggplot(rux_pk_sim, aes(x = time, y = rux_cp)) +
  geom_line(color = "#2980B9", linewidth = 1.2) +
  geom_hline(yintercept = 450, linetype = "dashed", color = "#E74C3C", linewidth = 0.8) +
  annotate("text", x = 110, y = 490, label = "IC50 JAK2 (~450 ng/mL)", color = "#E74C3C", size = 3.5) +
  labs(title    = "Ruxolitinib PK — 10mg BID, Multiple Dose",
       subtitle = "2-compartment model; t½ ~3h; Cmax ~1200 ng/mL",
       x        = "Time (hours)",
       y        = "Ruxolitinib Plasma Conc. (ng/mL)") +
  theme_bw(base_size = 12)

print(p_pk_rux)

## ============================================================
## Summary statistics
## ============================================================
cat("\n=== 5-Year Model Predictions Summary ===\n")
yr5_summary <- all_results %>%
  filter(abs(time - 1825) < 1) %>%
  select(scenario, HCT_pct, allele_out, PLT_out, WBC_out,
         spleen_out, tss_out, mf_out, thromb_out) %>%
  mutate(across(where(is.numeric), round, 1))

print(yr5_summary)

cat("\n=== Clinical Parameter Reference Ranges ===\n")
cat("HCT target: <45% (men), <42% (women) — CYTOREDUCE trial primary endpoint\n")
cat("SVR35: spleen volume ≥35% reduction from baseline at Week 32\n")
cat("TSS50: ≥50% reduction in MPN-SAF TSS — RESPONSE trial secondary endpoint\n")
cat("JAK2 allele burden: <50% (response), <10% (molecular remission) — PROUD-PV\n")
cat("Thrombosis risk: annual rate ~2.5% (low-risk PV with HCT control)\n")
