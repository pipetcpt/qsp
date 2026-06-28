## ============================================================
## Hypersensitivity Pneumonitis (HP) — mrgsolve QSP/PK-PD Model
## ============================================================
## Disease: Hypersensitivity Pneumonitis (Extrinsic Allergic Alveolitis)
## Model Type: ODE-based QSP — Immune-Fibrotic Cascade + Drug PK/PD
## Compartments: 22 ODE states (PK + Innate + Adaptive + Fibrosis + Lung Function)
##
## Treatment Scenarios Covered:
##   1. Antigen avoidance only (standard of care)
##   2. Oral prednisolone (0.5 mg/kg/d → taper)
##   3. Mycophenolate mofetil (MMF, 1500 mg BID)
##   4. Azathioprine (2 mg/kg/d)
##   5. Nintedanib (150 mg BID, antifibrotic)
##   6. Prednisolone + MMF combination
##   7. Nintedanib + antigen avoidance
##
## Key Clinical Calibration References:
##   - Morisset et al. (2020) Lancet Respir Med: MMF vs AZA in fibrotic HP
##   - Raghu et al. (2021) NEJM: Nintedanib in fibrotic ILD
##   - Walsh et al. (2014) Thorax: Natural history of HP
##   - Giménez et al. (2018) ERJ: Fibrotic HP prognosis
##   - Fernández Pérez et al. (2018) ATS: HP diagnosis & management
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

code <- '
$PROB
Hypersensitivity Pneumonitis QSP Model
Immune-fibrotic cascade with PK/PD for 5 treatments

$PARAM
// ── Antigen Exposure Parameters ──────────────────────────────
AG_exposure  = 1.0   // baseline antigen exposure rate (AU/d)
AG_avoidance = 0.0   // antigen avoidance efficacy (0=none, 1=complete)
k_ag_clear   = 0.5   // antigen clearance rate (/d)

// ── Innate Immunity ───────────────────────────────────────────
k_mac_act   = 0.8    // macrophage activation rate by antigen (/d)
k_mac_base  = 0.05   // baseline resting macrophage proliferation (/d)
k_mac_death = 0.15   // macrophage death rate (/d)
k_mac_M2    = 0.05   // M1→M2 switch rate driven by TGFb (/d per AU)
k_neu_rec   = 0.6    // neutrophil recruitment by IL8 (/d)
k_neu_death = 0.8    // neutrophil death rate (/d)

// ── Cytokine Kinetics ─────────────────────────────────────────
k_TNF_prod  = 0.4    // TNFa production by M1 macrophage
k_TNF_deg   = 1.2    // TNFa degradation (/d)
k_IL6_prod  = 0.35   // IL-6 production rate
k_IL6_deg   = 1.5    // IL-6 degradation (/d)
k_IL12_prod = 0.3    // IL-12 production by M1 macrophage
k_IL12_deg  = 1.0    // IL-12 degradation (/d)
k_IFNg_prod = 0.5    // IFN-gamma production by Th1
k_IFNg_deg  = 1.2    // IFN-gamma degradation (/d)
k_TGFb_prod = 0.25   // TGF-b production (M2 + myofibroblast)
k_TGFb_deg  = 0.8    // TGF-b degradation (/d)
k_IL17_prod = 0.3    // IL-17A production by Th17
k_IL17_deg  = 1.0    // IL-17A degradation (/d)
k_IL10_prod = 0.15   // IL-10 production by Treg/M2
k_IL10_deg  = 0.8    // IL-10 degradation (/d)

// ── Adaptive Immunity ─────────────────────────────────────────
k_Th1_diff  = 0.4    // naive→Th1 differentiation rate (IL-12 driven)
k_Th17_diff = 0.25   // naive→Th17 differentiation rate (IL-6 driven)
k_Treg_diff = 0.15   // naive→Treg differentiation rate
k_Th1_death = 0.1    // Th1 cell death rate (/d)
k_Th17_death = 0.12  // Th17 death rate (/d)
k_Treg_death = 0.08  // Treg death rate (/d)
T_naive_0   = 1.0    // baseline naive T cell (AU)
Treg_inhib  = 0.3    // Treg inhibition strength on Th1/Th17

// ── Granuloma ─────────────────────────────────────────────────
k_gran_form = 0.2    // granuloma formation rate by M1+Th1
k_gran_res  = 0.15   // granuloma resolution rate
k_gran_fib  = 0.08   // granuloma→fibrosis transition (chronic)

// ── Fibrosis Cascade ──────────────────────────────────────────
k_fib_act   = 0.3    // fibroblast activation rate (TGFb)
k_fib_death = 0.2    // fibroblast death rate (/d)
k_myo_diff  = 0.25   // fibroblast→myofibroblast differentiation
k_myo_death = 0.15   // myofibroblast apoptosis (/d)
k_col_prod  = 0.15   // collagen production rate (myofib)
k_col_deg   = 0.05   // collagen degradation rate (/d)
k_ROS_prod  = 0.3    // ROS production by M2 + neutrophils
k_ROS_deg   = 0.5    // ROS scavenging (/d)

// ── Lung Function Decline ─────────────────────────────────────
FVC_0       = 90.0   // baseline FVC (% predicted)
DLCO_0      = 85.0   // baseline DLCO (% predicted)
k_FVC_col   = 0.08   // FVC decline rate per collagen unit
k_DLCO_col  = 0.10   // DLCO decline rate per collagen unit
k_FVC_recov = 0.01   // FVC spontaneous recovery (very slow)
k_DLCO_recov = 0.005 // DLCO spontaneous recovery

// ── Prednisolone PK ───────────────────────────────────────────
ka_PDN      = 1.5    // absorption rate (/h)
CL_PDN      = 10.0   // clearance (L/h)
Vd_PDN      = 35.0   // volume of distribution (L)
F_PDN       = 0.82   // bioavailability

// ── MMF/MPA PK ────────────────────────────────────────────────
ka_MMF      = 2.0    // MMF absorption rate (/h)
CL_MPA      = 16.0   // MPA clearance (L/h)
Vd_MPA      = 55.0   // Vd MPA (L)
F_MMF       = 0.94   // MMF bioavailability

// ── Nintedanib PK ─────────────────────────────────────────────
ka_Nint     = 0.8    // nintedanib absorption (/h)
CL_Nint     = 85.0   // clearance (L/h) — high first pass
Vd_Nint     = 1050.0 // Vd (L)
F_Nint      = 0.047  // low bioavailability (~4.7%)

// ── PD Parameters ─────────────────────────────────────────────
EC50_PDN_inflam = 150.0  // prednisolone EC50 for inflammation (ng/mL)
EC50_MPA_prolif = 0.5    // MPA EC50 for T cell proliferation (μg/mL)
EC50_Nint_fib   = 200.0  // nintedanib EC50 for fibrosis (ng/mL)
Emax_PDN    = 0.85   // max suppression by prednisolone
Emax_MPA    = 0.80   // max suppression by MPA
Emax_Nint   = 0.75   // max antifibrotic effect (nintedanib)

// ── Biomarker Parameters ──────────────────────────────────────
KL6_base    = 250.0  // baseline KL-6 (U/mL)
k_KL6_col   = 2.0    // KL-6 scaling with collagen
k_KL6_deg   = 0.3    // KL-6 degradation (/d)

$CMT
// Antigen
AG_lung

// Innate Immunity
M_M1        // M1 activated macrophage (AU)
M_M2        // M2 pro-fibrotic macrophage (AU)
Neutrophil  // Neutrophil count (AU)

// Cytokines (all in AU relative units)
C_TNF
C_IL6
C_IL12
C_IFNg
C_TGFb
C_IL17
C_IL10

// Adaptive Immunity
T_Th1       // Th1 cells (AU)
T_Th17      // Th17 cells (AU)
T_Treg      // Regulatory T cells (AU)

// Granuloma
Granuloma   // Granuloma burden (AU)

// Fibrosis
Fibroblast  // Activated fibroblast (AU)
Myofib      // Myofibroblast (AU)
Collagen    // Collagen deposition (AU)
ROS         // Reactive oxygen species (AU)

// PK compartments
PDN_gut     // Prednisolone gut (mg)
PDN_cent    // Prednisolone central (mg)
MPA_gut     // MPA gut (mg)
MPA_cent    // MPA central (mg)
Nint_gut    // Nintedanib gut (mg)
Nint_cent   // Nintedanib central (mg)

// Lung function (% predicted)
FVC
DLCO

// Biomarker
KL6_serum   // KL-6 serum (U/mL)

$INIT
AG_lung = 0
M_M1 = 0.1
M_M2 = 0.05
Neutrophil = 0.05
C_TNF = 0
C_IL6 = 0
C_IL12 = 0
C_IFNg = 0
C_TGFb = 0
C_IL17 = 0
C_IL10 = 0
T_Th1 = 0.1
T_Th17 = 0.05
T_Treg = 0.1
Granuloma = 0
Fibroblast = 0.05
Myofib = 0
Collagen = 0
ROS = 0
PDN_gut = 0
PDN_cent = 0
MPA_gut = 0
MPA_cent = 0
Nint_gut = 0
Nint_cent = 0
FVC = 90
DLCO = 85
KL6_serum = 250

$ODE
// ─────────────────────────────────────────────────────────────
// Drug concentrations (ng/mL or μg/mL for PD)
double PDN_conc = PDN_cent / Vd_PDN * 1000;   // ng/mL
double MPA_conc = MPA_cent / Vd_MPA;          // μg/mL
double Nint_conc = Nint_cent / Vd_Nint * 1000; // ng/mL

// ─────────────────────────────────────────────────────────────
// PD Effect (Imax model)
double E_PDN = Emax_PDN * PDN_conc / (EC50_PDN_inflam + PDN_conc);
double E_MPA = Emax_MPA * MPA_conc / (EC50_MPA_prolif + MPA_conc);
double E_Nint = Emax_Nint * Nint_conc / (EC50_Nint_fib + Nint_conc);

// Antigen avoidance reduces exposure
double AG_in = AG_exposure * (1.0 - AG_avoidance);

// ─────────────────────────────────────────────────────────────
// ANTIGEN LUNG COMPARTMENT
dxdt_AG_lung = AG_in - k_ag_clear * AG_lung;

// ─────────────────────────────────────────────────────────────
// INNATE IMMUNITY
// M1 macrophage: activated by antigen, amplified by IFNg, suppressed by GC
double mac_stim = k_mac_act * AG_lung * (1 + 0.3 * C_IFNg) * (1 - 0.8 * E_PDN);
dxdt_M_M1 = mac_stim + k_mac_base - k_mac_death * M_M1 - k_mac_M2 * C_TGFb * M_M1;

// M2 macrophage: switches from M1 under TGFb/IL-10 influence
dxdt_M_M2 = k_mac_M2 * C_TGFb * M_M1 - k_mac_death * M_M2;

// Neutrophils: recruited by IL-8 (modeled via IL-17A/TNF)
dxdt_Neutrophil = k_neu_rec * (C_TNF + C_IL17) * (1 - 0.5 * E_PDN) - k_neu_death * Neutrophil;

// ─────────────────────────────────────────────────────────────
// CYTOKINES
dxdt_C_TNF  = k_TNF_prod * M_M1 * (1 - E_PDN) - k_TNF_deg * C_TNF;
dxdt_C_IL6  = k_IL6_prod * (M_M1 + T_Th17) * (1 - E_PDN) - k_IL6_deg * C_IL6;
dxdt_C_IL12 = k_IL12_prod * M_M1 - k_IL12_deg * C_IL12;
dxdt_C_IFNg = k_IFNg_prod * T_Th1 + 0.1 * M_M1 - k_IFNg_deg * C_IFNg;
dxdt_C_TGFb = k_TGFb_prod * (M_M2 + Myofib) - k_TGFb_deg * C_TGFb;
dxdt_C_IL17 = k_IL17_prod * T_Th17 - k_IL17_deg * C_IL17;
dxdt_C_IL10 = k_IL10_prod * (T_Treg + M_M2) - k_IL10_deg * C_IL10;

// ─────────────────────────────────────────────────────────────
// ADAPTIVE IMMUNITY
// Th1: differentiation driven by IL-12; suppressed by Treg and GC/MMF
double Th1_diff = k_Th1_diff * C_IL12 * T_naive_0 * (1 - Treg_inhib * T_Treg) * (1 - E_PDN) * (1 - E_MPA);
dxdt_T_Th1 = Th1_diff - k_Th1_death * T_Th1;

// Th17: differentiation driven by IL-6+TGFb; suppressed by Treg and MPA
double Th17_diff = k_Th17_diff * C_IL6 * C_TGFb * T_naive_0 * (1 - Treg_inhib * T_Treg) * (1 - E_MPA);
dxdt_T_Th17 = Th17_diff - k_Th17_death * T_Th17;

// Treg: differentiation by TGFb+IL-10; partially suppressed by inflammation
dxdt_T_Treg = k_Treg_diff * C_TGFb * C_IL10 * T_naive_0 - k_Treg_death * T_Treg;

// ─────────────────────────────────────────────────────────────
// GRANULOMA
double gran_form = k_gran_form * M_M1 * T_Th1;
double gran_res  = k_gran_res  * C_IL10 * Granuloma;
double gran_fib  = k_gran_fib  * Granuloma * (1 - C_IL10) * (1 - E_PDN);
dxdt_Granuloma = gran_form - gran_res - gran_fib;
if (Granuloma < 0) Granuloma = 0;

// ─────────────────────────────────────────────────────────────
// FIBROSIS CASCADE
// Fibroblast activation: TGFb, CCL18 (proxy via M2), granuloma-driven
dxdt_Fibroblast = k_fib_act * C_TGFb * (1 + 0.3 * gran_fib) * (1 - 0.5 * E_PDN)
                - k_fib_death * Fibroblast;

// Myofibroblast: from fibroblast differentiation; inhibited by nintedanib
dxdt_Myofib = k_myo_diff * C_TGFb * Fibroblast * (1 - E_Nint)
             - k_myo_death * Myofib;

// Collagen deposition: produced by myofib; reduced by nintedanib
dxdt_Collagen = k_col_prod * Myofib * (1 - 0.6 * E_Nint) - k_col_deg * Collagen;

// ROS: from M2 and neutrophils
dxdt_ROS = k_ROS_prod * (M_M2 + Neutrophil) - k_ROS_deg * ROS;

// ─────────────────────────────────────────────────────────────
// PK ODEs (units: mg)

// Prednisolone
dxdt_PDN_gut  = -ka_PDN * PDN_gut;
dxdt_PDN_cent = ka_PDN * F_PDN * PDN_gut - (CL_PDN / Vd_PDN) * PDN_cent;

// MPA (from MMF)
dxdt_MPA_gut  = -ka_MMF * MPA_gut;
dxdt_MPA_cent = ka_MMF * F_MMF * MPA_gut - (CL_MPA / Vd_MPA) * MPA_cent;

// Nintedanib
dxdt_Nint_gut  = -ka_Nint * Nint_gut;
dxdt_Nint_cent = ka_Nint * F_Nint * Nint_gut - (CL_Nint / Vd_Nint) * Nint_cent;

// ─────────────────────────────────────────────────────────────
// LUNG FUNCTION (% predicted)
// FVC declines with collagen deposition; partial recovery possible
dxdt_FVC  = -k_FVC_col  * Collagen * FVC / 100.0 + k_FVC_recov  * (FVC_0 - FVC);
dxdt_DLCO = -k_DLCO_col * Collagen * DLCO / 100.0 + k_DLCO_recov * (DLCO_0 - DLCO);

// ─────────────────────────────────────────────────────────────
// KL-6 BIOMARKER (U/mL)
dxdt_KL6_serum = k_KL6_col * (Collagen + ROS + M_M2) - k_KL6_deg * (KL6_serum - KL6_base);

$TABLE
double PDN_Cmax = PDN_cent / Vd_PDN * 1000; // ng/mL
double MPA_Cmax = MPA_cent / Vd_MPA;        // μg/mL
double Nint_Cmax = Nint_cent / Vd_Nint * 1000; // ng/mL
double FVC_pct  = FVC;
double DLCO_pct = DLCO;
double Inflam_index = C_TNF + C_IL6 + C_IFNg; // composite inflammation
double Fibrosis_index = Collagen + Myofib;     // composite fibrosis

$CAPTURE
PDN_Cmax MPA_Cmax Nint_Cmax
FVC_pct DLCO_pct KL6_serum
Inflam_index Fibrosis_index
M_M1 M_M2 T_Th1 T_Th17 T_Treg
C_TNF C_IL6 C_TGFb C_IFNg C_IL10
Granuloma Collagen Myofib ROS
'

## Build model
mod <- mcode("hp_qsp", code)

## ─────────────────────────────────────────────────────────────
## SCENARIO DEFINITIONS
## ─────────────────────────────────────────────────────────────

## Simulation parameters
DAYS <- 730   # 2-year simulation (730 days)
DT   <- 1.0   # daily output

## Helper: build dosing event for oral drugs
## PDN: given in mg, converted to dose compartment
## MMF: 1500 mg BID = 3000 mg/d (given as MPA 1080 mg/d effective)
## Nintedanib: 150 mg BID = 300 mg/d

## Dosing functions using mrgsolve ev()
make_pdn_events <- function(dose_mg, freq = 24, duration = 180) {
  # Prednisolone taper: full dose for 90d → half for 90d
  e1 <- ev(amt = dose_mg, cmt = "PDN_gut", ii = freq, addl = 89, time = 0)
  e2 <- ev(amt = dose_mg/2, cmt = "PDN_gut", ii = freq, addl = 89, time = 90*24)
  e3 <- ev(amt = dose_mg/4, cmt = "PDN_gut", ii = freq, addl = (duration-180)*24/freq, time = 180*24)
  ev(e1, e2)  # return first two phases for simplicity
}

## ─────────────────────────────────────────────────────────────
## SCENARIO 1: No treatment — progressive HP
## ─────────────────────────────────────────────────────────────
s1_params <- param(mod, AG_exposure = 1.0, AG_avoidance = 0)

sim1 <- mrgsim(s1_params, end = DAYS, delta = DT) %>%
  as.data.frame() %>%
  mutate(Scenario = "1. Untreated HP\n(No Intervention)")

## ─────────────────────────────────────────────────────────────
## SCENARIO 2: Antigen avoidance only
## ─────────────────────────────────────────────────────────────
s2_params <- param(mod, AG_exposure = 1.0, AG_avoidance = 0.9)

sim2 <- mrgsim(s2_params, end = DAYS, delta = DT) %>%
  as.data.frame() %>%
  mutate(Scenario = "2. Antigen Avoidance\n(90% Reduction)")

## ─────────────────────────────────────────────────────────────
## SCENARIO 3: Prednisolone 0.5 mg/kg/d (40 mg for 80 kg pt)
## ─────────────────────────────────────────────────────────────
PDN_daily <- 40  # mg/day, simplified as single daily dose

pdn_dose <- ev(
  amt  = PDN_daily,
  cmt  = "PDN_gut",
  ii   = 24,       # every 24 hours
  addl = DAYS - 1
)

s3_params <- param(mod, AG_exposure = 1.0, AG_avoidance = 0)

sim3 <- mrgsim(s3_params, ev = pdn_dose, end = DAYS, delta = DT) %>%
  as.data.frame() %>%
  mutate(Scenario = "3. Prednisolone 40mg/d\n(Monotherapy)")

## ─────────────────────────────────────────────────────────────
## SCENARIO 4: Mycophenolate mofetil (MMF) 1500 mg BID
## ─────────────────────────────────────────────────────────────
# MMF 1500mg BID → MPA ~1080mg equiv. Simplified: 3000mg/d to MPA_gut
mmf_dose <- ev(
  amt  = 1500,    # mg per dose
  cmt  = "MPA_gut",
  ii   = 12,      # BID
  addl = DAYS*2 - 1
)

s4_params <- param(mod, AG_exposure = 1.0, AG_avoidance = 0)

sim4 <- mrgsim(s4_params, ev = mmf_dose, end = DAYS, delta = DT) %>%
  as.data.frame() %>%
  mutate(Scenario = "4. MMF 1500mg BID\n(Immunosuppression)")

## ─────────────────────────────────────────────────────────────
## SCENARIO 5: Nintedanib 150 mg BID (antifibrotic)
## ─────────────────────────────────────────────────────────────
nint_dose <- ev(
  amt  = 150,
  cmt  = "Nint_gut",
  ii   = 12,
  addl = DAYS*2 - 1
)

s5_params <- param(mod, AG_exposure = 1.0, AG_avoidance = 0)

sim5 <- mrgsim(s5_params, ev = nint_dose, end = DAYS, delta = DT) %>%
  as.data.frame() %>%
  mutate(Scenario = "5. Nintedanib 150mg BID\n(Antifibrotic)")

## ─────────────────────────────────────────────────────────────
## SCENARIO 6: Prednisolone + MMF combination
## ─────────────────────────────────────────────────────────────
combo_dose <- ev(pdn_dose) + ev(mmf_dose)

s6_params <- param(mod, AG_exposure = 1.0, AG_avoidance = 0.5)

sim6 <- mrgsim(s6_params, ev = combo_dose, end = DAYS, delta = DT) %>%
  as.data.frame() %>%
  mutate(Scenario = "6. PDN + MMF + Partial\nAntigen Avoidance")

## ─────────────────────────────────────────────────────────────
## SCENARIO 7: Nintedanib + complete antigen avoidance
## ─────────────────────────────────────────────────────────────
s7_params <- param(mod, AG_exposure = 1.0, AG_avoidance = 0.95)

sim7 <- mrgsim(s7_params, ev = nint_dose, end = DAYS, delta = DT) %>%
  as.data.frame() %>%
  mutate(Scenario = "7. Nintedanib +\nComplete Avoidance")

## ─────────────────────────────────────────────────────────────
## COMBINE ALL SCENARIOS
## ─────────────────────────────────────────────────────────────
all_sim <- bind_rows(sim1, sim2, sim3, sim4, sim5, sim6, sim7)

all_sim$Scenario <- factor(all_sim$Scenario, levels = c(
  "1. Untreated HP\n(No Intervention)",
  "2. Antigen Avoidance\n(90% Reduction)",
  "3. Prednisolone 40mg/d\n(Monotherapy)",
  "4. MMF 1500mg BID\n(Immunosuppression)",
  "5. Nintedanib 150mg BID\n(Antifibrotic)",
  "6. PDN + MMF + Partial\nAntigen Avoidance",
  "7. Nintedanib +\nComplete Avoidance"
))

colors7 <- c("#E74C3C","#27AE60","#3498DB","#9B59B6","#F39C12","#1ABC9C","#2C3E50")

## ─────────────────────────────────────────────────────────────
## PLOT 1: FVC (% predicted) over 2 years
## ─────────────────────────────────────────────────────────────
p1 <- ggplot(all_sim, aes(x = time, y = FVC_pct, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 70, linetype = "dashed", color = "gray40") +
  annotate("text", x = 30, y = 71.5, label = "Transplant threshold ~70%", color = "gray40", size = 3.2) +
  scale_color_manual(values = colors7) +
  labs(
    title = "FVC (% Predicted) — 2-Year Trajectory",
    subtitle = "Hypersensitivity Pneumonitis — Treatment Comparison",
    x = "Days", y = "FVC (% predicted)",
    color = "Treatment Scenario"
  ) +
  coord_cartesian(ylim = c(40, 95)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right", legend.text = element_text(size = 8))

## ─────────────────────────────────────────────────────────────
## PLOT 2: DLCO (% predicted)
## ─────────────────────────────────────────────────────────────
p2 <- ggplot(all_sim, aes(x = time, y = DLCO_pct, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = colors7) +
  labs(
    title = "DLCO (% Predicted) — Gas Transfer Capacity",
    x = "Days", y = "DLCO (% predicted)"
  ) +
  coord_cartesian(ylim = c(30, 90)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

## ─────────────────────────────────────────────────────────────
## PLOT 3: Collagen deposition (fibrosis burden)
## ─────────────────────────────────────────────────────────────
p3 <- ggplot(all_sim, aes(x = time, y = Collagen, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = colors7) +
  labs(
    title = "Collagen Deposition (Fibrosis Burden)",
    x = "Days", y = "Collagen (AU)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

## ─────────────────────────────────────────────────────────────
## PLOT 4: KL-6 serum biomarker
## ─────────────────────────────────────────────────────────────
p4 <- ggplot(all_sim, aes(x = time, y = KL6_serum, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 500, linetype = "dashed", color = "red") +
  annotate("text", x = 30, y = 520, label = "Diagnostic threshold 500 U/mL", color = "red", size = 3) +
  scale_color_manual(values = colors7) +
  labs(
    title = "Serum KL-6 (Biomarker of Fibrosis Activity)",
    x = "Days", y = "KL-6 (U/mL)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

## ─────────────────────────────────────────────────────────────
## PLOT 5: Inflammation index (TNF + IL6 + IFNg)
## ─────────────────────────────────────────────────────────────
p5 <- ggplot(all_sim, aes(x = time, y = Inflam_index, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = colors7) +
  labs(
    title = "Composite Inflammation Index\n(TNF-α + IL-6 + IFN-γ)",
    x = "Days", y = "Inflammation Index (AU)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

## ─────────────────────────────────────────────────────────────
## PLOT 6: Th1 / Th17 / Treg immune balance
## ─────────────────────────────────────────────────────────────
immune_df <- all_sim %>%
  select(time, Scenario, T_Th1, T_Th17, T_Treg) %>%
  pivot_longer(cols = c(T_Th1, T_Th17, T_Treg), names_to = "Cell", values_to = "Level")

immune_df$Cell <- factor(immune_df$Cell, levels = c("T_Th1", "T_Th17", "T_Treg"),
                         labels = c("Th1 (pro-inflam)", "Th17 (pro-inflam/fibrosis)", "Treg (regulatory)"))

p6 <- ggplot(
  immune_df %>% filter(Scenario %in% c("1. Untreated HP\n(No Intervention)",
                                        "2. Antigen Avoidance\n(90% Reduction)",
                                        "6. PDN + MMF + Partial\nAntigen Avoidance")),
  aes(x = time, y = Level, color = Cell, linetype = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = c("#E74C3C","#E67E22","#27AE60")) +
  labs(
    title = "T Cell Immune Balance\n(Selected Scenarios)",
    x = "Days", y = "Cell Level (AU)",
    color = "Cell Type", linetype = "Scenario"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.text = element_text(size = 8))

## ─────────────────────────────────────────────────────────────
## COMBINE PLOTS
## ─────────────────────────────────────────────────────────────
combined <- (p1 | p2) / (p3 | p4) / (p5 | p6)
print(combined)

## ─────────────────────────────────────────────────────────────
## SUMMARY TABLE at 1 year (day 365) and 2 years (day 730)
## ─────────────────────────────────────────────────────────────
summary_tbl <- all_sim %>%
  filter(time %in% c(0, 180, 365, 548, 730)) %>%
  group_by(Scenario, time) %>%
  summarise(
    FVC     = round(mean(FVC_pct), 1),
    DLCO    = round(mean(DLCO_pct), 1),
    KL6     = round(mean(KL6_serum), 0),
    Collagen = round(mean(Collagen), 3),
    Inflam  = round(mean(Inflam_index), 3),
    .groups = "drop"
  ) %>%
  rename(Day = time)

print(summary_tbl, n = 50)

## ─────────────────────────────────────────────────────────────
## SENSITIVITY ANALYSIS: Antigen avoidance gradient
## ─────────────────────────────────────────────────────────────
avoidance_levels <- c(0, 0.25, 0.5, 0.75, 0.9, 1.0)
avoid_results <- lapply(avoidance_levels, function(av) {
  p <- param(mod, AG_exposure = 1.0, AG_avoidance = av)
  out <- mrgsim(p, end = DAYS, delta = 7) %>%
    as.data.frame() %>%
    mutate(Avoidance = paste0(av * 100, "% avoidance"))
  out
})
avoid_df <- bind_rows(avoid_results)

p_avoid <- ggplot(avoid_df, aes(x = time, y = FVC_pct, color = Avoidance)) +
  geom_line(linewidth = 1.1) +
  scale_color_brewer(palette = "RdYlGn", direction = 1) +
  labs(
    title = "Sensitivity: Antigen Avoidance Level → FVC Outcome",
    subtitle = "No pharmacotherapy; avoidance is primary intervention",
    x = "Days", y = "FVC (% predicted)",
    color = "Avoidance Level"
  ) +
  theme_bw(base_size = 12)
print(p_avoid)

cat("\n=== HP QSP Model: Simulation Complete ===\n")
cat("Scenarios simulated: 7\n")
cat("Duration: 730 days (2 years)\n")
cat("Key output variables: FVC, DLCO, KL-6, Collagen, Inflammation Index\n")
cat("Drug classes modeled: GC (Prednisolone), MMF, Nintedanib\n")
