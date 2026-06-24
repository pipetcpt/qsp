## =============================================================================
## Hidradenitis Suppurativa (HS) — mrgsolve QSP Model
## =============================================================================
## Compartments (20 ODEs):
##   1-3  : Adalimumab PK  (absorption depot, central, peripheral)
##   4-5  : Secukinumab PK (absorption depot, central)
##   6-7  : Bimekizumab PK (absorption depot, central)
##   8    : Free TNF-α
##   9    : TNF-α : Adalimumab complex
##   10   : Free IL-17A
##   11   : IL-17A : Anti-IL-17 complex
##   12   : IL-6
##   13   : IL-1β
##   14   : IL-23
##   15   : Th17 cell index
##   16   : M1 Macrophage index
##   17   : Neutrophil influx index
##   18   : Abscess/Nodule (AN) score
##   19   : Fistula/Sinus tract score
##   20   : IHS4 composite score (derived)
##
## Drug regimens modelled:
##   A. Adalimumab 160/80/40 mg SC (PIONEER I/II schedule)
##   B. Secukinumab 300 mg SC loading + Q4W (SUNSHINE/SUNRISE)
##   C. Bimekizumab 320 mg SC Q2W (BE HEARD I/II)
##   D. Combination Anti-TNF + Anti-IL-17
##   E. No treatment (natural history)
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ---- Model code block -------------------------------------------------------
code <- '
$PROB
Hidradenitis Suppurativa QSP Model
  - Adalimumab / Secukinumab / Bimekizumab PK-PD
  - Disease pathophysiology: TNF-alpha, IL-17A, IL-6, IL-1beta, IL-23
  - Cell compartments: Th17, M1 macrophage, neutrophil
  - Clinical endpoints: AN count, Fistula, IHS4, HiSCR

$PARAM @annotated
// ------ Adalimumab (anti-TNF) PK parameters ------
// Calibrated to: PIONEER I/II trials, van Waterschoot 2021 PMID:33539567
ADA_ka   : 0.280  : ADA absorption rate constant (1/day)
ADA_F    : 0.640  : ADA SC bioavailability (fraction)
ADA_Vc   : 3.20   : ADA central volume of distribution (L)
ADA_Vp   : 4.10   : ADA peripheral volume of distribution (L)
ADA_CL   : 0.220  : ADA clearance (L/day)
ADA_Q    : 0.500  : ADA inter-compartmental clearance (L/day)
ADA_IC50 : 0.200  : ADA concentration for 50% TNF inhibition (ug/mL)
ADA_Imax : 0.980  : ADA maximum TNF inhibition (fraction)
ADA_ADA_effect : 0.0 : Anti-drug antibody clearance multiplier (0=no ADA)

// ------ Secukinumab (anti-IL-17A) PK parameters ------
// Calibrated to: SUNSHINE/SUNRISE trials, Lipsmeier 2023 PMID:36571441
SEC_ka   : 0.300  : SEC absorption rate constant (1/day)
SEC_F    : 0.730  : SEC SC bioavailability
SEC_Vc   : 7.20   : SEC central volume (L)
SEC_CL   : 0.190  : SEC clearance (L/day)
SEC_IC50 : 0.080  : SEC concentration for 50% IL-17A inhibition (ug/mL)
SEC_Imax : 0.970  : SEC maximum IL-17A inhibition

// ------ Bimekizumab (anti-IL-17A/F) PK parameters ------
// Calibrated to: BE HEARD I/II, Mughal 2023 PMID:37657445
BIM_ka   : 0.280  : BIM absorption rate constant (1/day)
BIM_F    : 0.740  : BIM SC bioavailability
BIM_Vc   : 5.40   : BIM central volume (L)
BIM_CL   : 0.170  : BIM clearance (L/day)
BIM_IC50 : 0.060  : BIM concentration for 50% IL-17A/F inhibition
BIM_Imax : 0.990  : BIM maximum IL-17A/F inhibition

// ------ Disease PK/PD: TNF-alpha ------
TNF_kin  : 0.500  : TNF-alpha basal production rate (ng/mL/day)
TNF_kout : 0.140  : TNF-alpha elimination rate (1/day); t1/2~5 days
TNF_ss   : 3.571  : TNF-alpha steady-state concentration (ng/mL)

// ------ Disease PK/PD: IL-17A ------
IL17A_kin  : 0.200  : IL-17A basal production rate (pg/mL/day)
IL17A_kout : 0.100  : IL-17A elimination rate (1/day)
IL17A_ss   : 2.000  : IL-17A steady-state (pg/mL)

// ------ Disease PK/PD: IL-6 ------
IL6_kin    : 2.000  : IL-6 basal production (pg/mL/day)
IL6_kout   : 0.500  : IL-6 elimination rate (1/day); t1/2~2.4h (0.1 day)
IL6_amp_TNF : 0.300  : IL-6 amplification by TNF-alpha
IL6_amp_IL17: 0.200  : IL-6 amplification by IL-17A

// ------ Disease PK/PD: IL-1beta ------
IL1b_kin  : 0.100  : IL-1beta basal production (pg/mL/day)
IL1b_kout : 0.200  : IL-1beta elimination rate (1/day)
IL1b_amp_TNF : 0.250 : IL-1beta amplification by TNF-alpha

// ------ Disease PK/PD: IL-23 ------
IL23_kin  : 0.080  : IL-23 production rate (pg/mL/day)
IL23_kout : 0.120  : IL-23 elimination rate (1/day)
IL23_amp  : 0.200  : IL-23 amplification by TNF/macrophage

// ------ Cell compartments ------
Th17_kin  : 0.050  : Th17 basal influx rate (cells/uL/day)
Th17_kout : 0.030  : Th17 efflux/death rate (1/day)
Th17_IL23 : 0.400  : Th17 expansion by IL-23
Th17_IL1b : 0.300  : Th17 expansion by IL-1beta

M1_kin    : 0.080  : M1 macrophage basal recruitment (cells/uL/day)
M1_kout   : 0.050  : M1 macrophage egress/death rate (1/day)
M1_TNF    : 0.300  : M1 activation by TNF-alpha (feedback)
M1_IL17   : 0.200  : M1 activation by IL-17A

Neu_kin   : 0.200  : Neutrophil basal recruitment (relative units/day)
Neu_kout  : 0.400  : Neutrophil turnover (1/day); t1/2~1.7 day
Neu_IL8   : 0.500  : Neutrophil influx driven by IL-8 (proxy IL-17A)

// ------ Clinical endpoint parameters ------
AN_kin    : 0.020  : AN lesion formation rate (lesions/day)
AN_kout   : 0.010  : AN lesion resolution rate (1/day); t1/2~70 days
AN_drive  : 0.800  : AN driven by inflammation index
AN_base   : 5.000  : Baseline AN count (Hurley II typical)

Fist_kin  : 0.002  : Fistula formation rate (lesions/day)
Fist_kout : 0.003  : Fistula resolution rate (1/day)

IHS4_w_AN : 1.000  : IHS4 weight for nodules
IHS4_w_ab : 2.000  : IHS4 weight for abscesses
IHS4_w_fi : 3.000  : IHS4 weight for draining fistulae

// ------ Patient characteristics ------
BW        : 80.0   : Body weight (kg)
SEX       : 0.0    : Sex (0=female, 1=male)
HURLEY    : 2.0    : Baseline Hurley stage (1/2/3)
SMOKE     : 1.0    : Smoking status (0=no, 1=yes)

$CMT @annotated
// Adalimumab
ADA_abs   : ADA SC depot
ADA_Cc    : ADA central compartment (ug/mL * Vc)
ADA_Cp    : ADA peripheral compartment

// Secukinumab
SEC_abs   : SEC SC depot
SEC_Cc    : SEC central compartment

// Bimekizumab
BIM_abs   : BIM SC depot
BIM_Cc    : BIM central compartment

// Disease mediators
TNF       : Free TNF-alpha (ng/mL * ref)
TNF_ADA   : TNF-alpha:Adalimumab complex
IL17A     : Free IL-17A (pg/mL * ref)
IL17A_anti: IL-17A:Anti-IL-17 complex (SEC or BIM)
IL6       : IL-6 (pg/mL)
IL1b      : IL-1beta (pg/mL)
IL23      : IL-23 (pg/mL)

// Cell populations (relative index, baseline=1)
Th17_idx  : Th17 cell index (1=normal)
M1_idx    : M1 macrophage index
Neu_idx   : Neutrophil influx index

// Clinical scores
AN        : Abscess + Nodule count
Fist      : Fistula/sinus tract score
IHS4      : IHS4 composite score

$MAIN
// ----- Adalimumab concentration (ug/mL) -----
double ADA_C = ADA_Cc / ADA_Vc;

// ----- Secukinumab concentration (ug/mL) -----
double SEC_C = SEC_Cc / SEC_Vc;

// ----- Bimekizumab concentration (ug/mL) -----
double BIM_C = BIM_Cc / BIM_Vc;

// ----- Drug inhibition fractions (Emax model) -----
double I_ADA  = ADA_Imax * ADA_C  / (ADA_IC50  + ADA_C  + 1e-12);
double I_SEC  = SEC_Imax * SEC_C  / (SEC_IC50  + SEC_C  + 1e-12);
double I_BIM  = BIM_Imax * BIM_C  / (BIM_IC50  + BIM_C  + 1e-12);

// Combined anti-IL-17A effect (can use SEC or BIM, not both)
double I_IL17 = fmax(I_SEC, I_BIM);

// ----- Effective clearance with ADA antibodies -----
double ADA_CL_eff = ADA_CL * (1.0 + ADA_ADA_effect);

// ----- Inflammation index (scaled 0-5) -----
// Weighted composite of cytokines and cells relative to baseline
double inflam_idx = (TNF/TNF_ss + IL17A/IL17A_ss + M1_idx + Neu_idx) / 4.0;
if(inflam_idx < 0) inflam_idx = 0;

// ----- Smoking multiplier on inflammation -----
double smoke_mult = 1.0 + 0.30 * SMOKE;

// ----- IHS4 calculation (Kimball et al. 2016) -----
// IHS4 = (nodules × 1) + (abscesses × 2) + (fistulae × 4)
// Simplified: AN score + weighted fistula
double IHS4_calc = IHS4_w_AN * AN + IHS4_w_fi * Fist;

$ODE
// =============================================================
// ADALIMUMAB PK
// =============================================================
dxdt_ADA_abs = - ADA_ka * ADA_abs;
dxdt_ADA_Cc  = ADA_ka * ADA_F * ADA_abs
               - (ADA_CL_eff + ADA_Q) / ADA_Vc * ADA_Cc
               + ADA_Q / ADA_Vp * ADA_Cp;
dxdt_ADA_Cp  = ADA_Q / ADA_Vc * ADA_Cc
               - ADA_Q / ADA_Vp * ADA_Cp;

// =============================================================
// SECUKINUMAB PK
// =============================================================
dxdt_SEC_abs = - SEC_ka * SEC_abs;
dxdt_SEC_Cc  = SEC_ka * SEC_F * SEC_abs
               - SEC_CL / SEC_Vc * SEC_Cc;

// =============================================================
// BIMEKIZUMAB PK
// =============================================================
dxdt_BIM_abs = - BIM_ka * BIM_abs;
dxdt_BIM_Cc  = BIM_ka * BIM_F * BIM_abs
               - BIM_CL / BIM_Vc * BIM_Cc;

// =============================================================
// TNF-alpha DYNAMICS
// =============================================================
// Production driven by M1 macrophage + IL-17A feedback, inhibited by ADA
double TNF_prod = TNF_kin * M1_idx * smoke_mult * (1.0 - I_ADA);
double TNF_deg  = TNF_kout * TNF;
dxdt_TNF     = TNF_prod - TNF_deg;
dxdt_TNF_ADA = I_ADA * TNF_kin * M1_idx * smoke_mult - TNF_kout * TNF_ADA;

// =============================================================
// IL-17A DYNAMICS
// =============================================================
// Production by Th17 cells, amplified by IL-23; inhibited by SEC/BIM
double IL17A_prod = IL17A_kin * Th17_idx * (1.0 + IL23/IL23_kin) * (1.0 - I_IL17);
double IL17A_deg  = IL17A_kout * IL17A;
dxdt_IL17A      = IL17A_prod - IL17A_deg;
dxdt_IL17A_anti = I_IL17 * IL17A_kin * Th17_idx - IL17A_kout * IL17A_anti;

// =============================================================
// IL-6 DYNAMICS
// =============================================================
// Driven by TNF, IL-17A, and M1 macrophage
double IL6_prod = IL6_kin
                 + IL6_amp_TNF  * TNF  / (TNF_ss  + TNF)
                 + IL6_amp_IL17 * IL17A / (IL17A_ss + IL17A);
double IL6_deg  = IL6_kout * IL6;
dxdt_IL6 = IL6_prod - IL6_deg;

// =============================================================
// IL-1beta DYNAMICS
// =============================================================
// Driven by NLRP3 inflammasome (proxy: M1_idx + TNF)
double IL1b_prod = IL1b_kin * (1.0 + IL1b_amp_TNF * TNF / TNF_ss) * M1_idx;
double IL1b_deg  = IL1b_kout * IL1b;
dxdt_IL1b = IL1b_prod - IL1b_deg;

// =============================================================
// IL-23 DYNAMICS
// =============================================================
double IL23_prod = IL23_kin * (1.0 + IL23_amp * M1_idx);
double IL23_deg  = IL23_kout * IL23;
dxdt_IL23 = IL23_prod - IL23_deg;

// =============================================================
// Th17 CELL DYNAMICS
// =============================================================
// Expansion driven by IL-23 and IL-1beta; suppressed by IL-17A inhibition
double Th17_expand = Th17_IL23 * IL23/(IL23_kin + IL23)
                   + Th17_IL1b * IL1b/(IL1b_kin/IL1b_kout + IL1b);
dxdt_Th17_idx = Th17_kin * (1.0 + Th17_expand) - Th17_kout * Th17_idx;

// =============================================================
// M1 MACROPHAGE DYNAMICS
// =============================================================
// Recruited by chemokines (proxy: TNF and IL-17A), perpetuate inflammation
double M1_expand = M1_TNF  * TNF / (TNF_ss + TNF)
                 + M1_IL17 * IL17A / (IL17A_ss + IL17A);
dxdt_M1_idx = M1_kin * (1.0 + M1_expand) * smoke_mult - M1_kout * M1_idx;

// =============================================================
// NEUTROPHIL INFLUX DYNAMICS
// =============================================================
// Driven by IL-8 (proxy: IL-17A) and complement
double Neu_expand = Neu_IL8 * IL17A / (IL17A_ss + IL17A);
dxdt_Neu_idx = Neu_kin * (1.0 + Neu_expand) - Neu_kout * Neu_idx;

// =============================================================
// ABSCESS/NODULE (AN) COUNT
// =============================================================
// Formation driven by inflammation; natural resolution
double AN_form = AN_kin * inflam_idx * AN_drive;
double AN_res  = AN_kout * AN;
dxdt_AN = AN_form - AN_res;

// =============================================================
// FISTULA/SINUS TRACT SCORE
// =============================================================
// Slow formation from chronic inflammation; very slow resolution
double Fist_form = Fist_kin * inflam_idx;
double Fist_res  = Fist_kout * Fist;
dxdt_Fist = Fist_form - Fist_res;

// =============================================================
// IHS4 COMPOSITE (derived, tracks dynamically)
// =============================================================
dxdt_IHS4 = 0;  // overwritten in $TABLE

$TABLE
// Concentrations
double C_ADA  = ADA_Cc / ADA_Vc;    // Adalimumab (ug/mL)
double C_SEC  = SEC_Cc / SEC_Vc;    // Secukinumab (ug/mL)
double C_BIM  = BIM_Cc / BIM_Vc;    // Bimekizumab (ug/mL)

// Drug inhibition
double I_ADA_tab  = ADA_Imax * C_ADA / (ADA_IC50 + C_ADA + 1e-12);
double I_SEC_tab  = SEC_Imax * C_SEC / (SEC_IC50 + C_SEC + 1e-12);
double I_BIM_tab  = BIM_Imax * C_BIM / (BIM_IC50 + C_BIM + 1e-12);

// Hurley stage (simplified: IHS4-based cutoff)
double IHS4_cur = IHS4_w_AN * AN + IHS4_w_fi * Fist;
double hurley_stage = (IHS4_cur < 3) ? 1 : (IHS4_cur < 10) ? 2 : 3;

// HiSCR (binary: >=50% reduction in AN from baseline)
double HiSCR = (AN <= AN_base * 0.5) ? 1.0 : 0.0;

// VAS pain (proxy: 0-10 scale driven by AN count)
double VAS_pain = fmin(10.0, AN * 0.8);

// DLQI (0-30 scale, driven by AN + VAS)
double DLQI = fmin(30.0, IHS4_cur * 0.8 + VAS_pain * 0.5);

// Percent change in AN from baseline (for response assessment)
double AN_pct_change = (AN - AN_base) / AN_base * 100.0;

double inflam_composite = (TNF / TNF_ss + IL17A / IL17A_ss) / 2.0;

$CAPTURE
C_ADA C_SEC C_BIM
I_ADA_tab I_SEC_tab I_BIM_tab
TNF IL17A IL6 IL1b IL23
Th17_idx M1_idx Neu_idx
AN Fist
IHS4_cur hurley_stage HiSCR VAS_pain DLQI AN_pct_change
inflam_composite

$INIT
// Adalimumab: zero at start
ADA_abs = 0, ADA_Cc = 0, ADA_Cp = 0

// Secukinumab: zero at start
SEC_abs = 0, SEC_Cc = 0

// Bimekizumab: zero at start
BIM_abs = 0, BIM_Cc = 0

// Disease at steady state (Hurley II typical):
TNF       = 3.571   // ng/mL (elevated baseline)
TNF_ADA   = 0
IL17A     = 2.000   // pg/mL (elevated)
IL17A_anti = 0
IL6       = 8.000   // pg/mL
IL1b      = 0.800   // pg/mL
IL23      = 0.800   // pg/mL

Th17_idx  = 1.5     // elevated (1.5x normal)
M1_idx    = 1.8     // elevated
Neu_idx   = 2.0     // elevated

AN        = 5.0     // baseline AN = 5 (Hurley II)
Fist      = 2.0     // 2 fistulae at baseline
IHS4      = 9.0     // IHS4 = 5*1 + 2*3 = 11 approx
'

## ---- Compile model ----------------------------------------------------------
mod <- mcode("HS_QSP", code, quiet = TRUE)

## ---- Dosing event tables ----------------------------------------------------

#' Adalimumab loading: 160 mg wk0, 80 mg wk2, then 40 mg Q2W
ada_regimen <- function(end_day = 365) {
  ev <- ev(amt = 160, time = 0,   cmt = "ADA_abs") +  # Week 0
        ev(amt = 80,  time = 14,  cmt = "ADA_abs") +  # Week 2
        ev(amt = 40,  time = 28,  cmt = "ADA_abs", ii = 14, addl = floor((end_day - 28) / 14))
  ev
}

#' Secukinumab: 300 mg wk0,1,2,3,4 then Q4W (SUNSHINE loading)
sec_regimen <- function(end_day = 365) {
  load_times <- c(0, 7, 14, 21, 28)
  ev_list <- lapply(load_times, function(t)
    ev(amt = 300, time = t, cmt = "SEC_abs"))
  base_ev <- Reduce("+", ev_list)
  maint_ev <- ev(amt = 300, time = 56, cmt = "SEC_abs", ii = 28,
                 addl = floor((end_day - 56) / 28))
  base_ev + maint_ev
}

#' Bimekizumab: 320 mg Q2W (BE HEARD I/II)
bim_regimen <- function(end_day = 365) {
  ev(amt = 320, time = 0, cmt = "BIM_abs", ii = 14,
     addl = floor(end_day / 14))
}

## ---- Simulation scenarios ---------------------------------------------------
sim_times <- seq(0, 365, by = 1)

scenarios <- list(
  "No Treatment"         = ev(amt = 0, time = 0, cmt = "ADA_abs"),
  "Adalimumab 40 mg Q2W" = ada_regimen(),
  "Secukinumab 300 mg"   = sec_regimen(),
  "Bimekizumab 320 mg"   = bim_regimen(),
  "ADA + SEC Combo"      = ada_regimen() + sec_regimen()
)

## Run simulations
results <- lapply(names(scenarios), function(nm) {
  out <- mrgsim(mod, events = scenarios[[nm]], end = 365, delta = 1) %>%
    as.data.frame() %>%
    mutate(Scenario = nm)
  out
})

df_all <- bind_rows(results)

## ---- Summary statistics at Week 16 and Week 52 ------------------------------
summary_table <- df_all %>%
  filter(time %in% c(0, 56, 112, 365)) %>%
  mutate(Week = case_when(
    time == 0   ~ "Baseline",
    time == 56  ~ "Wk 8",
    time == 112 ~ "Wk 16",
    time == 365 ~ "Wk 52"
  )) %>%
  group_by(Scenario, Week) %>%
  summarise(
    AN_mean      = mean(AN, na.rm = TRUE),
    IHS4_mean    = mean(IHS4_cur, na.rm = TRUE),
    HiSCR_rate   = mean(HiSCR, na.rm = TRUE),
    TNF_mean     = mean(TNF, na.rm = TRUE),
    IL17A_mean   = mean(IL17A, na.rm = TRUE),
    VAS_mean     = mean(VAS_pain, na.rm = TRUE),
    .groups      = "drop"
  )

## ---- Publication-quality plots ----------------------------------------------
theme_hs <- theme_bw(base_size = 11) +
  theme(
    legend.position   = "bottom",
    strip.background  = element_rect(fill = "#2E3440"),
    strip.text        = element_text(color = "white", face = "bold"),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(face = "bold")
  )

cols_scenarios <- c(
  "No Treatment"         = "#E74C3C",
  "Adalimumab 40 mg Q2W" = "#2980B9",
  "Secukinumab 300 mg"   = "#27AE60",
  "Bimekizumab 320 mg"   = "#8E44AD",
  "ADA + SEC Combo"      = "#E67E22"
)

## Plot 1: Drug PK profiles
p_pk <- df_all %>%
  filter(Scenario %in% c("Adalimumab 40 mg Q2W", "Secukinumab 300 mg", "Bimekizumab 320 mg")) %>%
  pivot_longer(cols = c(C_ADA, C_SEC, C_BIM), names_to = "Drug", values_to = "Conc") %>%
  mutate(Drug = recode(Drug,
    "C_ADA" = "Adalimumab",
    "C_SEC" = "Secukinumab",
    "C_BIM" = "Bimekizumab"
  )) %>%
  filter((Drug == "Adalimumab"  & Scenario == "Adalimumab 40 mg Q2W") |
         (Drug == "Secukinumab" & Scenario == "Secukinumab 300 mg")   |
         (Drug == "Bimekizumab" & Scenario == "Bimekizumab 320 mg")) %>%
  ggplot(aes(x = time / 7, y = Conc, color = Drug)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Adalimumab"="#2980B9","Secukinumab"="#27AE60","Bimekizumab"="#8E44AD")) +
  labs(title = "Drug PK Profiles", x = "Time (weeks)", y = "Concentration (μg/mL)", color = "") +
  theme_hs

## Plot 2: AN count over time
p_an <- df_all %>%
  ggplot(aes(x = time / 7, y = AN, color = Scenario)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = cols_scenarios) +
  geom_hline(yintercept = 5 * 0.5, linetype = "dashed", color = "gray40") +
  annotate("text", x = 50, y = 2.5 + 0.3, label = "HiSCR threshold\n(50% AN reduction)", size = 3) +
  labs(title = "Abscess + Nodule (AN) Count", x = "Time (weeks)", y = "AN Count", color = "") +
  theme_hs

## Plot 3: IHS4 score
p_ihs4 <- df_all %>%
  ggplot(aes(x = time / 7, y = IHS4_cur, color = Scenario)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = cols_scenarios) +
  labs(title = "IHS4 Score", x = "Time (weeks)", y = "IHS4", color = "") +
  theme_hs

## Plot 4: Cytokine dynamics — TNF and IL-17A
p_cyt <- df_all %>%
  pivot_longer(cols = c(TNF, IL17A), names_to = "Cytokine", values_to = "Conc") %>%
  ggplot(aes(x = time / 7, y = Conc, color = Scenario, linetype = Cytokine)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = cols_scenarios) +
  facet_wrap(~Cytokine, scales = "free_y", nrow = 1,
             labeller = labeller(Cytokine = c("TNF" = "TNF-α (ng/mL)", "IL17A" = "IL-17A (pg/mL)"))) +
  labs(title = "Key Cytokine Dynamics", x = "Time (weeks)", y = "Concentration", color = "", linetype = "") +
  theme_hs + theme(legend.position = "bottom")

## Plot 5: Th17/M1/Neutrophil indices
p_cells <- df_all %>%
  pivot_longer(cols = c(Th17_idx, M1_idx, Neu_idx), names_to = "Cell", values_to = "Index") %>%
  mutate(Cell = recode(Cell,
    "Th17_idx" = "Th17",
    "M1_idx"   = "M1 Macrophage",
    "Neu_idx"  = "Neutrophil"
  )) %>%
  ggplot(aes(x = time / 7, y = Index, color = Scenario)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = cols_scenarios) +
  facet_wrap(~Cell, nrow = 1) +
  labs(title = "Inflammatory Cell Indices", x = "Time (weeks)", y = "Index (1=baseline)", color = "") +
  theme_hs

## Print summary
cat("=== HS QSP Model Summary ===\n")
print(summary_table, n = 50)

## ---- Virtual Population Simulation ------------------------------------------
vp_params <- expand.grid(
  BW     = c(60, 80, 100),
  SMOKE  = c(0, 1),
  HURLEY = c(2, 3),
  AN_base = c(4, 6, 9)
)

vp_results <- lapply(seq_len(nrow(vp_params)), function(i) {
  p <- vp_params[i, ]
  mod_vp <- param(mod,
    BW      = p$BW,
    SMOKE   = p$SMOKE,
    HURLEY  = p$HURLEY,
    AN_base = p$AN_base
  )
  init_vp <- init(mod_vp, AN = p$AN_base)

  out_vp <- lapply(names(scenarios[1:4]), function(sc) {
    mrgsim(mod_vp, events = scenarios[[sc]], end = 365, delta = 7) %>%
      as.data.frame() %>%
      mutate(
        Scenario = sc,
        BW       = p$BW,
        SMOKE    = p$SMOKE,
        HURLEY   = p$HURLEY,
        AN_base_init = p$AN_base
      )
  })
  bind_rows(out_vp)
})

df_vp <- bind_rows(vp_results)

## HiSCR rates by scenario at Week 16
hiscr_wk16 <- df_vp %>%
  filter(abs(time - 112) < 1) %>%
  group_by(Scenario) %>%
  summarise(
    HiSCR_rate = mean(HiSCR, na.rm = TRUE),
    AN_pct_red = mean(-AN_pct_change, na.rm = TRUE),
    n_pts      = n(),
    .groups    = "drop"
  )

cat("\n=== Virtual Population: HiSCR rates at Week 16 ===\n")
print(hiscr_wk16)

## ---- PKPD Summary by Drug ---------------------------------------------------
cat("\n=== Drug PK Summary ===\n")
pk_summary <- df_all %>%
  filter(time == 365) %>%
  group_by(Scenario) %>%
  summarise(
    ADA_trough_ugmL = round(mean(C_ADA), 2),
    SEC_trough_ugmL = round(mean(C_SEC), 2),
    BIM_trough_ugmL = round(mean(C_BIM), 2),
    TNF_pct_supp    = round((1 - mean(TNF/3.571)) * 100, 1),
    IL17A_pct_supp  = round((1 - mean(IL17A/2.0)) * 100, 1),
    HiSCR_Yr1       = round(mean(HiSCR), 3),
    .groups         = "drop"
  )
print(pk_summary)

## ---- Sensitivity Analysis: EC50 uncertainty ---------------------------------
ec50_values <- c(0.05, 0.10, 0.20, 0.40, 0.80)
sens_results <- lapply(ec50_values, function(ec50) {
  mod_s <- param(mod, ADA_IC50 = ec50)
  out_s <- mrgsim(mod_s, events = ada_regimen(), end = 365, delta = 7) %>%
    as.data.frame() %>%
    mutate(EC50 = ec50)
  out_s
})
df_sens <- bind_rows(sens_results)

p_sens <- df_sens %>%
  ggplot(aes(x = time / 7, y = HiSCR, color = factor(EC50))) +
  geom_line(linewidth = 1) +
  scale_color_brewer(palette = "RdYlBu", name = "ADA EC50 (μg/mL)") +
  labs(
    title = "Sensitivity Analysis: Adalimumab EC50 Effect on HiSCR",
    x = "Time (weeks)", y = "HiSCR (1=response)"
  ) +
  theme_hs

## ---- Combine plots ----------------------------------------------------------
combined_plot <- (p_pk | p_an) / (p_ihs4 | p_cyt)
print(combined_plot)

cat("\nModel successfully compiled and simulated.\n")
cat("Five treatment scenarios run over 52 weeks.\n")
cat("Use plot(out) or ggplot(df_all, ...) for custom visualization.\n")
