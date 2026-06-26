## ============================================================
## Melanoma QSP Model — mrgsolve ODE Implementation
## BRAF-mutant cutaneous melanoma: targeted therapy + immunotherapy
##
## Key Biology:
##   BRAF V600E/K → constitutive MEK-ERK → proliferation
##   BRAFi/MEKi → ERK inhibition → tumor apoptosis (but resistance emerges)
##   PD-1/CTLA-4 blockade → T cell reactivation → immune-mediated killing
##
## Compartments (16 ODEs):
##   PK: GUT_BRAF, CENT_BRAF, GUT_MEK, CENT_MEK, CENT_ICI, PERI_ICI
##   PD: ERK_act, RESIST, TUMOR, CD8_TIL, PD1_RO, CTLA4_RO,
##       Treg_frac, IFNg_TME, LDH_ser, S100B_ser
##
## Treatment Scenarios (6):
##   1. Untreated (BRAF V600E metastatic melanoma)
##   2. Vemurafenib 960mg BID (BRAF inhibitor mono)
##   3. Dabrafenib 150mg BID + Trametinib 2mg QD (BRAFi+MEKi combo)
##   4. Pembrolizumab 200mg IV q3w (PD-1 inhibitor)
##   5. Nivolumab 1mg/kg + Ipilimumab 3mg/kg q3w x4 (dual checkpoint)
##   6. Dabrafenib+Trametinib (wk 0-24) → Pembrolizumab (wk 24+)
##
## Clinical calibration references:
##   BRIM-3 (Chapman 2011 NEJM): Vmfnb ORR 48%, mPFS 5.3mo vs 1.6mo
##   COMBI-d (Long 2014 NEJM): Dbfnb+Tram mPFS 9.3mo, 5yr OS 28%
##   CheckMate 067 (Larkin 2015 NEJM): Nivo+Ipi mPFS 11.5mo, 5yr OS 52%
##   KEYNOTE-006 (Robert 2015 NEJM): Pembro mPFS 5.6mo, 3yr OS 50%
## ============================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

## ============================================================
## Model Code Block
## ============================================================
melanoma_code <- '
$PARAM
// ---- BRAF Inhibitor PK (Vemurafenib/Dabrafenib-like) ----
ka_B    = 0.42,   // /h  absorption rate
F_B     = 0.64,   // oral bioavailability
Vd_B    = 100.0,  // L   apparent volume
CL_B    = 13.0,   // L/h apparent clearance (t1/2 ~5.3h PK-proxy)
DOSE_B  = 0.0,    // mg  daily dose (oral, BID = dose/2 q12h)

// ---- MEK Inhibitor PK (Trametinib-like) ----
ka_M    = 0.35,   // /h
F_M     = 0.72,   // oral bioavailability
Vd_M    = 214.0,  // L
CL_M    = 4.9,    // L/h (t1/2 ~30h trametinib, long)
DOSE_M  = 0.0,    // mg  daily dose (oral, QD)

// ---- Immune Checkpoint Inhibitor PK (Pembrolizumab/Nivolumab-like IgG4) ----
CL_I    = 0.244,  // L/day central clearance
Vc_I    = 3.61,   // L   central volume
Vp_I    = 2.75,   // L   peripheral volume
Q_I     = 0.448,  // L/day inter-compartmental
DOSE_I  = 0.0,    // mg  flat dose per infusion (IV bolus approximation)
DOSE_C4 = 0.0,    // mg  CTLA-4 inhibitor dose (ipilimumab-like)

// ---- ERK PD: tumor MAPK activity ----
kout_ERK    = 0.12,  // /h  ERK turnover
EC50_BRAFi  = 0.28,  // µg/mL BRAFi IC50 (on ERK suppression)
Emax_BRAFi  = 0.90,  // max ERK suppression by BRAFi
EC50_MEKi   = 12.0,  // ng/mL MEKi IC50
Emax_MEKi   = 0.85,  // max ERK suppression by MEKi
hill_B      = 1.8,   // Hill coefficient BRAFi
hill_M      = 1.5,   // Hill coefficient MEKi

// ---- Tumor Growth PD (modified Gompertz) ----
kg0         = 0.0055, // /h  intrinsic tumor growth rate
kmax_tum    = 10.0,   // fold max tumor mass (Gompertz ceiling)
kd_ERK      = 0.0040, // /h  BRAFi/MEKi kill rate (ERK-dep)
kd_immune   = 0.0050, // /h  immune-mediated kill rate (CD8 TIL dep)
BRAF_V600   = 1.0,    // 1=V600E/K (BRAFi sens), 0=WT/NRAS (no BRAFi)

// ---- Acquired Resistance ----
kR_on       = 0.00005, // /h  resistance accumulation (under BRAFi)
kR_off      = 0.0,     // /h  partial reversibility (simplified 0)
Rmax        = 1.0,     // max resistance = full resistance

// ---- CD8+ TIL Dynamics ----
kin_CD8     = 0.018,  // /h  basal CD8 TIL recruitment
kout_CD8    = 0.010,  // /h  CD8 TIL decay/emigration
CD8_ICI_fac = 1.8,    // fold PD-1 blockade boost on CD8 TIL (per RO unit)
CD8_CTLA4   = 1.4,    // fold CTLA-4 blockade boost
Treg_sup    = 0.60,   // Treg suppression factor on CD8 (max 60%)

// ---- PD-1 / CTLA-4 Receptor Occupancy ----
EC50_PD1    = 12.0,   // µg/mL half-effect conc PD-1 occupancy
EC50_CTLA4  = 20.0,   // µg/mL CTLA-4 occupancy

// ---- Treg Fraction in TME ----
Treg0       = 0.15,   // baseline Treg fraction (15%)
Treg_ICI    = 0.40,   // CTLA-4 blockade reduces Treg by 40%

// ---- IFN-gamma in TME ----
kin_IFNg    = 0.002,  // /h  basal IFN-gamma production (norm)
kout_IFNg   = 0.05,   // /h  IFN-gamma decay
IFNg_CD8    = 0.02,   // IFN-gamma production per unit CD8 TIL
kPDL1_IFNg  = 0.30,   // IFN-gamma → PD-L1 induction (secondary resistance)

// ---- LDH (serum) ----
LDH0        = 200.0,  // U/L baseline LDH
LDH_tum_fac = 1.8,    // max LDH fold-increase at peak tumor burden
kout_LDH    = 0.03,   // /h  LDH half-life

// ---- S100B (serum) ----
S100B_0     = 0.10,   // µg/L baseline S100B
S100B_tum   = 2.0,    // fold tumor-dependent S100B

$CMT
GUT_BRAF    // oral BRAF inhibitor (µg)
CENT_BRAF   // plasma BRAF inhibitor (µg/mL * L = µg)
GUT_MEK     // oral MEK inhibitor (ng equiv)
CENT_MEK    // plasma MEK inhibitor (ng/mL * L = ng)
CENT_ICI    // checkpoint inhibitor central (mg)
PERI_ICI    // checkpoint inhibitor peripheral (mg)
ERK_act     // ERK relative activity (1 = uninhibited)
RESIST      // resistance factor (0 = sensitive, 1 = fully resistant)
TUMOR       // tumor burden (normalized: 1 = baseline)
CD8_TIL     // CD8+ TIL (normalized: 1 = baseline)
PD1_RO      // PD-1 receptor occupancy (0-1)
CTLA4_RO    // CTLA-4 receptor occupancy (0-1)
Treg_frac   // Treg fraction in TME (0-1)
IFNg_TME    // IFN-gamma in tumor microenvironment (normalized)
LDH_ser     // serum LDH (U/L)
S100B_ser   // S100B (µg/L)

$MAIN
// BRAFi/MEKi plasma concentrations (µg/mL and ng/mL)
double Cb = CENT_BRAF / Vd_B;    // µg/mL BRAF inhibitor
double Cm = CENT_MEK  / Vd_M;    // ng/mL MEK inhibitor
double Ci = CENT_ICI  / Vc_I;    // µg/mL checkpoint inhibitor

// ERK inhibition (Emax models; combined BRAFi and MEKi effect)
double inh_B = BRAF_V600 * Emax_BRAFi * pow(Cb, hill_B) /
               (pow(EC50_BRAFi, hill_B) + pow(Cb, hill_B));
double inh_M = Emax_MEKi * pow(Cm, hill_M) /
               (pow(EC50_MEKi, hill_M) + pow(Cm, hill_M));
// Combination: additive-to-synergistic ERK suppression (Loewe)
double ERK_inh = inh_B + inh_M - inh_B * inh_M;
double ERK_target = 1.0 - ERK_inh;  // target ERK level (0 = full inhib)

// ERK reactivation by acquired resistance (R shifts ERK back toward 1)
ERK_target = ERK_target + RESIST * (1.0 - ERK_target);
if (ERK_target < 0.0) ERK_target = 0.0;
if (ERK_target > 1.0) ERK_target = 1.0;

// PD-1 receptor occupancy
double PD1_eq  = Ci / (Ci + EC50_PD1);
// CTLA-4 occupancy (separate ICI or dual)
double CTA4_eq = Ci / (Ci + EC50_CTLA4);

// IFN-gamma-induced PD-L1 (adaptive resistance to immunotherapy)
double PDL1_ind = 1.0 + kPDL1_IFNg * IFNg_TME;

// Effective CD8 TIL activity (modulated by Treg, PD-1, PD-L1)
double PD1_sup   = 1.0 - PD1_RO;          // residual PD-1 suppression
double CD8_eff   = CD8_TIL * (1.0 - Treg_frac * Treg_sup) *
                   (1.0 - 0.70 * PD1_sup / PDL1_ind);
if (CD8_eff < 0.0) CD8_eff = 0.0;

// Tumor kill rate
double kill_ERK    = kd_ERK    * (1.0 - ERK_act) * BRAF_V600;
double kill_immune = kd_immune * CD8_eff;

// Initial conditions (set once at simulation start)
TUMOR_0    = 1.0;
CD8_TIL_0  = 1.0;
ERK_act_0  = 1.0;
RESIST_0   = 0.0;
PD1_RO_0   = 0.0;
CTLA4_RO_0 = 0.0;
Treg_frac_0= Treg0;
IFNg_TME_0 = 1.0;
LDH_ser_0  = LDH0;
S100B_ser_0= S100B_0;

$ODE
// ---- BRAF Inhibitor PK ----
dxdt_GUT_BRAF  = -ka_B * GUT_BRAF;
dxdt_CENT_BRAF = F_B * ka_B * GUT_BRAF - (CL_B / Vd_B) * CENT_BRAF;

// ---- MEK Inhibitor PK ----
dxdt_GUT_MEK   = -ka_M * GUT_MEK;
dxdt_CENT_MEK  = F_M * ka_M * GUT_MEK - (CL_M / Vd_M) * CENT_MEK;

// ---- ICI PK (2-compartment, IV infusion/bolus) ----
double Ci_cur = CENT_ICI / Vc_I;
double Cp_cur = PERI_ICI / Vp_I;
dxdt_CENT_ICI  = -(CL_I / Vc_I) * CENT_ICI - Q_I * (Ci_cur - Cp_cur);
dxdt_PERI_ICI  =  Q_I * (Ci_cur - Cp_cur);

// ---- ERK Activity (indirect response: production offset by inhibition) ----
// kin = kout*(1) so at SS ERK_act = 1; inhibition lowers production
double kin_ERK = kout_ERK * ERK_target;
dxdt_ERK_act   = kin_ERK - kout_ERK * ERK_act;

// ---- Acquired Resistance (grows under BRAFi/MEKi exposure) ----
double drug_pressure = (CENT_BRAF > 0 || CENT_MEK > 0) ? 1.0 : 0.0;
dxdt_RESIST    = kR_on * drug_pressure * (Rmax - RESIST) - kR_off * RESIST;

// ---- Tumor Burden (modified Gompertz-Simeoni) ----
// Gompertz growth: kg0 * TUMOR * ln(kmax_tum / TUMOR)
double tumor_grow  = kg0 * TUMOR * log(kmax_tum / TUMOR);
double tumor_kill  = (kill_ERK + kill_immune) * TUMOR;
dxdt_TUMOR     = tumor_grow - tumor_kill;

// ---- CD8+ TIL ----
double cd8_recruit = kin_CD8 * (1.0 + (CD8_ICI_fac - 1.0) * PD1_RO +
                                      (CD8_CTLA4 - 1.0) * CTLA4_RO);
dxdt_CD8_TIL   = cd8_recruit - kout_CD8 * CD8_TIL;

// ---- PD-1 Receptor Occupancy ----
dxdt_PD1_RO    = 0.5 * (PD1_eq  - PD1_RO);   // first-order approach to Emax

// ---- CTLA-4 Receptor Occupancy ----
dxdt_CTLA4_RO  = 0.5 * (CTA4_eq - CTLA4_RO);

// ---- Treg Fraction (CTLA-4 blockade depletes Tregs) ----
double Treg_eq  = Treg0 * (1.0 - Treg_ICI * CTLA4_RO);
dxdt_Treg_frac  = 0.1 * (Treg_eq - Treg_frac);  // slow equilibration

// ---- IFN-gamma in TME ----
double kin_IFNg_eff = kin_IFNg + IFNg_CD8 * CD8_TIL * PD1_RO;
dxdt_IFNg_TME  = kin_IFNg_eff - kout_IFNg * IFNg_TME;

// ---- LDH (tumor necrosis surrogate) ----
double LDH_eq   = LDH0 * (1.0 + (LDH_tum_fac - 1.0) * (TUMOR - 1.0));
if (LDH_eq < LDH0 * 0.5) LDH_eq = LDH0 * 0.5;
dxdt_LDH_ser   = kout_LDH * (LDH_eq - LDH_ser);

// ---- S100B (melanoma burden marker) ----
double S100B_eq = S100B_0 * (1.0 + (S100B_tum - 1.0) * (TUMOR - 1.0));
if (S100B_eq < S100B_0 * 0.3) S100B_eq = S100B_0 * 0.3;
dxdt_S100B_ser = 0.05 * (S100B_eq - S100B_ser);

$CAPTURE
// Derived outputs for plotting
double Cb_um    = CENT_BRAF / Vd_B;     // BRAFi plasma (µg/mL)
double Cm_nm    = CENT_MEK  / Vd_M;     // MEKi plasma (ng/mL)
double Ci_um    = CENT_ICI  / Vc_I;     // ICI plasma (µg/mL)
double ERK_pct  = ERK_act * 100.0;      // ERK activity %
double Tumor_pct = TUMOR * 100.0;       // tumor % of baseline
double LDH_val  = LDH_ser;
double S100B_val = S100B_ser;
double CD8_val  = CD8_TIL;
double Resist_pct = RESIST * 100.0;
double PD1_RO_pct = PD1_RO * 100.0;
double IFNg_val   = IFNg_TME;
double Treg_pct   = Treg_frac * 100.0;
'

## Build the model
mel_mod <- mcode("melanoma_qsp", melanoma_code)

## ============================================================
## Helper: Build Dosing Events
## ============================================================
build_ev_oral <- function(dose_mg, tau_h, duration_wk) {
  n_doses <- floor(duration_wk * 7 * 24 / tau_h)
  ev(amt = dose_mg, ii = tau_h, addl = n_doses - 1, cmt = 1, rate = 0)
}

build_ev_iv_q3w <- function(dose_mg, duration_wk, cmt_num = 5) {
  n_doses <- floor(duration_wk / 3)
  ev(amt = dose_mg, ii = 21 * 24, addl = max(0, n_doses - 1),
     cmt = cmt_num, rate = dose_mg / 0.5)   # 30-min infusion
}

## ============================================================
## Simulation Parameters
## ============================================================
duration_wk <- 52    # 1-year simulation
delta_t     <- 1     # hourly output
sim_times   <- seq(0, duration_wk * 7 * 24, by = delta_t)

## ============================================================
## Scenario 1: No treatment (untreated BRAF V600E HetFH baseline)
## ============================================================
scen1 <- mel_mod %>%
  param(BRAF_V600 = 1.0) %>%
  mrgsim(ev = ev(), end = max(sim_times), delta = delta_t) %>%
  as.data.frame() %>%
  mutate(Scenario = "1. Untreated (BRAF V600E)")

## ============================================================
## Scenario 2: Vemurafenib 960mg BID (BRAFi monotherapy)
##   BRIM-3: ORR 48%, mPFS 5.3mo; calibrate kg0, kd_ERK
## ============================================================
ev2 <- ev(amt = 960, ii = 12, addl = 364 * 2 - 1, cmt = 1)   # BID oral

scen2 <- mel_mod %>%
  param(DOSE_B = 960, BRAF_V600 = 1.0) %>%
  mrgsim(ev = ev2, end = max(sim_times), delta = delta_t) %>%
  as.data.frame() %>%
  mutate(Scenario = "2. Vemurafenib 960mg BID")

## ============================================================
## Scenario 3: Dabrafenib 150mg BID + Trametinib 2mg QD (COMBI-d)
##   mPFS 9.3mo; ORR 67%; 5yr OS 28%
## ============================================================
ev3_braf <- ev(amt = 150, ii = 12, addl = 364 * 2 - 1, cmt = 1)
ev3_mek  <- ev(amt = 2,   ii = 24, addl = 364 - 1,     cmt = 3)
ev3 <- ev3_braf + ev3_mek

scen3 <- mel_mod %>%
  param(DOSE_B = 150, DOSE_M = 2, BRAF_V600 = 1.0,
        kd_ERK = 0.0045) %>%
  mrgsim(ev = ev3, end = max(sim_times), delta = delta_t) %>%
  as.data.frame() %>%
  mutate(Scenario = "3. Dabrafenib+Trametinib")

## ============================================================
## Scenario 4: Pembrolizumab 200mg IV q3w (KEYNOTE-006)
##   ORR 33%, mPFS 5.6mo, 3yr OS 50%
## ============================================================
ev4 <- ev(amt = 200, ii = 21 * 24, addl = 17, cmt = 5,
          rate = 200 / 0.5)  # q3w x 18 cycles = 54 weeks

scen4 <- mel_mod %>%
  param(BRAF_V600 = 1.0, DOSE_I = 200,
        kd_immune = 0.0060) %>%
  mrgsim(ev = ev4, end = max(sim_times), delta = delta_t) %>%
  as.data.frame() %>%
  mutate(Scenario = "4. Pembrolizumab 200mg q3w")

## ============================================================
## Scenario 5: Nivolumab 1mg/kg + Ipilimumab 3mg/kg q3w x4,
##             then Nivo 3mg/kg q2w (CheckMate 067)
##   mPFS 11.5mo, 5yr OS 52%, ORR 58%
## ============================================================
ev5_combo <- ev(amt = 200, ii = 21 * 24, addl = 3, cmt = 5,
                rate = 200 / 0.5)  # first 4 doses combo
ev5_ipi   <- ev(amt = 300, ii = 21 * 24, addl = 3, cmt = 5,
                rate = 300 / 0.5)
ev5_maint <- ev(amt = 240, ii = 14 * 24, addl = 25, time = 4 * 21 * 24,
                cmt = 5, rate = 240 / 0.5)   # nivo maintenance
ev5 <- ev5_combo + ev5_ipi + ev5_maint

scen5 <- mel_mod %>%
  param(BRAF_V600 = 1.0, DOSE_I = 200, DOSE_C4 = 300,
        kd_immune = 0.0080, CD8_ICI_fac = 2.2, CD8_CTLA4 = 1.8,
        Treg_ICI = 0.60) %>%
  mrgsim(ev = ev5, end = max(sim_times), delta = delta_t) %>%
  as.data.frame() %>%
  mutate(Scenario = "5. Nivolumab+Ipilimumab")

## ============================================================
## Scenario 6: Sequential BRAFi/MEKi → Immunotherapy
##   Dabrafenib+Trametinib wk 0-24, then Pembrolizumab wk 24+
## ============================================================
dur_braf_mek_h <- 24 * 7 * 24   # 24 weeks in hours
ev6_braf <- ev(amt = 150, ii = 12, addl = 24*7*2 - 1, cmt = 1)
ev6_mek  <- ev(amt = 2,   ii = 24, addl = 24*7 - 1,   cmt = 3)
ev6_ici  <- ev(amt = 200, ii = 21 * 24, addl = 12,
               time = dur_braf_mek_h, cmt = 5, rate = 200/0.5)
ev6 <- ev6_braf + ev6_mek + ev6_ici

scen6 <- mel_mod %>%
  param(DOSE_B = 150, DOSE_M = 2, DOSE_I = 200, BRAF_V600 = 1.0,
        kd_ERK = 0.0045, kd_immune = 0.0060) %>%
  mrgsim(ev = ev6, end = max(sim_times), delta = delta_t) %>%
  as.data.frame() %>%
  mutate(Scenario = "6. BRAFi/MEKi → Pembrolizumab")

## ============================================================
## Combine all scenarios and convert time to weeks
## ============================================================
all_scen <- bind_rows(scen1, scen2, scen3, scen4, scen5, scen6) %>%
  mutate(time_wk = time / (7 * 24))

## ============================================================
## Summary Statistics (at Week 12, 24, 52)
## ============================================================
summary_tbl <- all_scen %>%
  filter(time_wk %in% c(0, 12, 24, 52)) %>%
  group_by(Scenario, time_wk) %>%
  summarise(
    Tumor_pct   = mean(Tumor_pct, na.rm = TRUE),
    LDH_val     = mean(LDH_val, na.rm = TRUE),
    CD8_val     = mean(CD8_val, na.rm = TRUE),
    ERK_pct     = mean(ERK_pct, na.rm = TRUE),
    Resist_pct  = mean(Resist_pct, na.rm = TRUE),
    PD1_RO_pct  = mean(PD1_RO_pct, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_tbl)

## ============================================================
## Plots
## ============================================================
theme_mel <- theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(face = "bold"))

scen_colors <- c(
  "1. Untreated (BRAF V600E)"     = "#E53935",
  "2. Vemurafenib 960mg BID"      = "#FF8F00",
  "3. Dabrafenib+Trametinib"      = "#F9A825",
  "4. Pembrolizumab 200mg q3w"    = "#1E88E5",
  "5. Nivolumab+Ipilimumab"       = "#00897B",
  "6. BRAFi/MEKi → Pembrolizumab" = "#8E24AA"
)

## Plot 1: Tumor Burden Over Time
p1 <- ggplot(all_scen %>% filter(time_wk <= 52),
             aes(x = time_wk, y = Tumor_pct, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 130, linetype = "dotted", color = "red",
             alpha = 0.6) +
  geom_hline(yintercept = 70, linetype = "dotted", color = "green3",
             alpha = 0.7) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Tumor Burden Over Time",
       x = "Time (weeks)", y = "Tumor Burden (% Baseline)") +
  annotate("text", x = 50, y = 132, label = "Progressive Disease (>120%)",
           size = 3, color = "red") +
  annotate("text", x = 50, y = 62, label = "Partial Response (<70%)",
           size = 3, color = "green4") +
  theme_mel
print(p1)

## Plot 2: ERK Activity (BRAFi/MEKi effect)
p2 <- ggplot(all_scen %>% filter(time_wk <= 12),
             aes(x = time_wk, y = ERK_pct, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = scen_colors) +
  labs(title = "ERK Activity (% of Baseline)",
       x = "Time (weeks)", y = "ERK Relative Activity (%)") +
  theme_mel
print(p2)

## Plot 3: Acquired Resistance Accumulation
p3 <- ggplot(all_scen %>% filter(time_wk <= 52, Scenario %in% c(
  "2. Vemurafenib 960mg BID",
  "3. Dabrafenib+Trametinib",
  "6. BRAFi/MEKi → Pembrolizumab")),
             aes(x = time_wk, y = Resist_pct, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Acquired Resistance Accumulation (BRAFi/MEKi)",
       x = "Time (weeks)", y = "Resistance Score (%)") +
  theme_mel
print(p3)

## Plot 4: CD8+ TIL Dynamics
p4 <- ggplot(all_scen %>% filter(time_wk <= 52),
             aes(x = time_wk, y = CD8_val, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = scen_colors) +
  labs(title = "CD8+ TIL Dynamics",
       x = "Time (weeks)", y = "CD8+ TIL (Relative to Baseline)") +
  theme_mel
print(p4)

## Plot 5: Serum LDH Over Time
p5 <- ggplot(all_scen %>% filter(time_wk <= 52),
             aes(x = time_wk, y = LDH_val, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 250, linetype = "dashed", color = "red",
             alpha = 0.7) +
  scale_color_manual(values = scen_colors) +
  labs(title = "Serum LDH Over Time",
       x = "Time (weeks)", y = "LDH (U/L)") +
  annotate("text", x = 45, y = 260, label = "ULN (250 U/L)",
           size = 3, color = "red") +
  theme_mel
print(p5)

## ============================================================
## Waterfall Plot: Best % Change in Tumor Burden (Week 24)
## ============================================================
wf_data <- all_scen %>%
  filter(time_wk == 24) %>%
  group_by(Scenario) %>%
  summarise(BestChange = mean(Tumor_pct - 100), .groups = "drop") %>%
  arrange(BestChange)

p6 <- ggplot(wf_data, aes(x = reorder(Scenario, BestChange),
                           y = BestChange, fill = Scenario)) +
  geom_col() +
  geom_hline(yintercept = -30, linetype = "dashed", color = "green4") +
  geom_hline(yintercept = 20,  linetype = "dashed", color = "red") +
  scale_fill_manual(values = scen_colors) +
  coord_flip() +
  labs(title = "Waterfall Plot: Best % Change from Baseline (Week 24)",
       x = NULL, y = "% Change from Baseline") +
  annotate("text", x = 0.5, y = -32, label = "PR threshold (-30%)",
           size = 3, color = "green4", hjust = 0) +
  annotate("text", x = 0.5, y = 22, label = "PD threshold (+20%)",
           size = 3, color = "red", hjust = 0) +
  theme_mel + theme(legend.position = "none")
print(p6)

message("Melanoma QSP Model complete. 6 scenarios simulated.")
message("Key readouts: Tumor burden, ERK activity, CD8+ TIL, resistance, LDH, S100B")
