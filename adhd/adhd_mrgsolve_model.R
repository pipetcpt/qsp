## ============================================================
## ADHD QSP Model — mrgsolve ODE
## Attention Deficit Hyperactivity Disorder
## Dopamine · Norepinephrine · PFC Executive Function · Drug PK/PD
##
## Compartments (25 total):
##   MPH PK  : GUT1, CENT1, PER1                          [1-3]
##   AMP PK  : GUT2, CENT2, PER2                          [4-6]
##   ATX PK  : GUT3, CENT3, PER3                          [7-9]
##   GFN PK  : GUT4, CENT4, PER4                          [10-12]
##   VLX PK  : GUT5, CENT5                                [13-14]
##   Disease Biology:
##     DA_syn  : Synaptic dopamine                         [15]
##     NE_syn  : Synaptic norepinephrine                   [16]
##     DAT_occ : DAT occupancy (fraction)                  [17]
##     NET_occ : NET occupancy (fraction)                  [18]
##     PFC_DA  : PFC dopamine tone (index 0-1)             [19]
##     PFC_NE  : PFC norepinephrine tone (index 0-1)       [20]
##     WM_idx  : Working memory index (0-1, higher=better) [21]
##     ExecFun : Executive function index (0-1)            [22]
##   Clinical Endpoints:
##     ADHD_RS : ADHD-RS-5 total score (0-54)              [23]
##     CGI_S   : CGI-Severity (1-7)                        [24]
##     QoL_idx : Quality of life index (0-1)               [25]
##
## Calibration: Multiple clinical trials (CADDRA, MTA, AHRQ, NICE)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ---- Model specification ----------------------------------------------------
code <- '
$PROB ADHD QSP Model — DA/NE/PFC/Clinical Endpoints

$PARAM
// ── MPH PK parameters (Methylphenidate IR/ER)
ka1   = 1.2    // h-1, absorption (IR); use 0.4 for ER
CL1   = 31.5   // L/h, apparent clearance (70 kg adult)
V2c   = 448.0  // L, central volume
Q1    = 75.6   // L/h, inter-compartmental CL
V2p   = 560.0  // L, peripheral volume
F1    = 0.22   // bioavailability (IR)
Dose_MPH = 0   // mg, dose per admin

// ── AMP PK parameters (Mixed amphetamine salts)
ka2   = 1.0    // h-1
CL2   = 39.2   // L/h (EM; CYP2D6)
V3c   = 245.0  // L
Q2    = 42.0   // L/h
V3p   = 350.0  // L
F2    = 0.75   // bioavailability
Dose_AMP = 0   // mg

// ── ATX PK parameters (Atomoxetine; CYP2D6 EM)
ka3   = 2.0    // h-1
CL3   = 24.5   // L/h (EM); 5.5 L/h for PM
V4c   = 59.5   // L
Q3    = 28.0   // L/h
V4p   = 168.0  // L
F3    = 0.63   // EM; 0.94 for PM
Dose_ATX = 0   // mg

// ── GFN PK parameters (Guanfacine ER)
ka4   = 0.6    // h-1
CL4   = 3.5    // L/h
V5c   = 196.0  // L
Q4    = 5.6    // L/h
V5p   = 1120.0 // L (extensive distribution)
F4    = 0.80   // bioavailability
Dose_GFN = 0   // mg

// ── VLX PK parameters (Viloxazine ER)
ka5   = 0.5    // h-1
CL5   = 19.6   // L/h
V6c   = 105.0  // L
F5    = 0.88   // bioavailability
Dose_VLX = 0   // mg

// ── DA / NE synaptic kinetics
DA_base = 1.0  // nM, baseline synaptic DA
krel_DA = 0.5  // h-1, DA release rate
kel_DA  = 2.0  // h-1, DA elimination rate (DAT + MAO)
NE_base = 0.8  // nM, baseline synaptic NE
krel_NE = 0.4  // h-1, NE release rate
kel_NE  = 1.8  // h-1, NE elimination rate (NET + MAO)

// ── DAT/NET occupancy parameters
Ki_MPH_DAT = 0.05  // µM, Ki for MPH at DAT
Ki_MPH_NET = 0.34  // µM, Ki for MPH at NET
Ki_AMP_DAT = 0.10  // µM, Ki for AMP at DAT (substrate)
Ki_AMP_NET = 0.04  // µM, Ki for AMP at NET
Ki_ATX_NET = 0.002 // µM, Ki for ATX at NET
Ki_GFN_A2A = 0.001 // µM, Ki for GFN at α2A-AR
Ki_VLX_NET = 0.042 // µM, Ki for VLX at NET

// ── Emax model parameters for PFC tone
Emax_DA_WM = 0.8   // max DA effect on WM
EC50_DA_WM = 1.5   // nM, EC50 DA for WM
Emax_NE_WM = 0.9   // max NE effect on WM
EC50_NE_WM = 1.2   // nM, EC50 NE for WM
hill        = 1.5  // Hill coefficient for dose-response

// ── Inverted-U DA function parameters (PFC optimal tone)
DA_opt  = 2.0  // nM, optimal DA for PFC function
DA_bw   = 1.5  // bandwidth (sigma) of Gaussian
NE_opt  = 1.5  // nM, optimal NE
NE_bw   = 1.2  // bandwidth

// ── Symptom dynamics
ADHD_RS_base = 32.0  // baseline ADHD-RS total score (moderate ADHD)
k_symp  = 0.1   // h-1, rate of symptom change
k_WM_effect = 0.6  // relative weight of WM on symptom score

// ── Population variability flags
WT      = 70    // kg body weight
AGE     = 10    // years
CYP2D6  = 1     // 1=EM, 0.5=IM, 0.1=PM

$CMT
// Drug PK
GUT1 CENT1 PER1     // MPH
GUT2 CENT2 PER2     // AMP
GUT3 CENT3 PER3     // ATX
GUT4 CENT4 PER4     // GFN
GUT5 CENT5          // VLX
// Disease biology
DA_syn NE_syn DAT_occ NET_occ
PFC_DA PFC_NE WM_idx ExecFun
// Clinical endpoints
ADHD_RS CGI_S QoL_idx

$INIT
// PK: all zero at start
GUT1=0, CENT1=0, PER1=0
GUT2=0, CENT2=0, PER2=0
GUT3=0, CENT3=0, PER3=0
GUT4=0, CENT4=0, PER4=0
GUT5=0, CENT5=0
// Disease: baseline values
DA_syn=1.0, NE_syn=0.8
DAT_occ=0, NET_occ=0
PFC_DA=0.5, PFC_NE=0.5
WM_idx=0.35, ExecFun=0.30
ADHD_RS=32.0, CGI_S=4.0, QoL_idx=0.40

$ODE
// ─── MPH PK (2-compartment) ────────────────────────────────
double C1 = CENT1 / (V2c * WT/70.0);   // µg/mL = µM for MW=233.7
double CL1_adj = CL1 * pow(WT/70.0, 0.75);
dxdt_GUT1  = -ka1 * GUT1;
dxdt_CENT1 =  ka1 * GUT1 - (CL1_adj/V2c)*CENT1 - (Q1/V2c)*CENT1 + (Q1/V2p)*PER1;
dxdt_PER1  =  (Q1/V2c)*CENT1 - (Q1/V2p)*PER1;

// ─── AMP PK (2-compartment) ────────────────────────────────
double C2 = CENT2 / (V3c * WT/70.0);
double CL2_adj = CL2 * pow(WT/70.0, 0.75);
dxdt_GUT2  = -ka2 * GUT2;
dxdt_CENT2 =  ka2 * GUT2 - (CL2_adj/V3c)*CENT2 - (Q2/V3c)*CENT2 + (Q2/V3p)*PER2;
dxdt_PER2  =  (Q2/V3c)*CENT2 - (Q2/V3p)*PER2;

// ─── ATX PK (2-compartment, CYP2D6-dependent) ─────────────
double C3 = CENT3 / (V4c * WT/70.0);
double CL3_adj = CL3 * CYP2D6 * pow(WT/70.0, 0.75);
dxdt_GUT3  = -ka3 * GUT3;
dxdt_CENT3 =  ka3 * GUT3 - (CL3_adj/V4c)*CENT3 - (Q3/V4c)*CENT3 + (Q3/V4p)*PER3;
dxdt_PER3  =  (Q3/V4c)*CENT3 - (Q3/V4p)*PER3;

// ─── GFN PK (2-compartment) ────────────────────────────────
double C4 = CENT4 / (V5c * WT/70.0);
double CL4_adj = CL4 * pow(WT/70.0, 0.75);
dxdt_GUT4  = -ka4 * GUT4;
dxdt_CENT4 =  ka4 * GUT4 - (CL4_adj/V5c)*CENT4 - (Q4/V5c)*CENT4 + (Q4/V5p)*PER4;
dxdt_PER4  =  (Q4/V5c)*CENT4 - (Q4/V5p)*PER4;

// ─── VLX PK (1-compartment) ────────────────────────────────
double C5 = CENT5 / (V6c * WT/70.0);
double CL5_adj = CL5 * pow(WT/70.0, 0.75);
dxdt_GUT5  = -ka5 * GUT5;
dxdt_CENT5 =  ka5 * GUT5 - (CL5_adj/V6c)*CENT5;

// ─── DAT/NET occupancy (competitive inhibition) ───────────
// Convert µg/mL to µM:
//  MPH: MW 233.7 → µM = C1*1000/233.7
//  AMP: MW 135.2 → µM = C2*1000/135.2
//  ATX: MW 291.8 → µM = C3*1000/291.8
//  GFN: MW 246.7 → µM = C4*1000/246.7
//  VLX: MW 237.3 → µM = C5*1000/237.3
double MPH_uM = C1 * 1000.0 / 233.7;
double AMP_uM = C2 * 1000.0 / 135.2;
double ATX_uM = C3 * 1000.0 / 291.8;
double GFN_uM = C4 * 1000.0 / 246.7;
double VLX_uM = C5 * 1000.0 / 237.3;

// DAT occupancy (fraction 0-1)
double DAT_occ_ss = (MPH_uM/Ki_MPH_DAT + AMP_uM/Ki_AMP_DAT) /
                    (1.0 + MPH_uM/Ki_MPH_DAT + AMP_uM/Ki_AMP_DAT);
// NET occupancy (fraction 0-1)
double NET_occ_ss = (MPH_uM/Ki_MPH_NET + AMP_uM/Ki_AMP_NET +
                     ATX_uM/Ki_ATX_NET + VLX_uM/Ki_VLX_NET) /
                    (1.0 + MPH_uM/Ki_MPH_NET + AMP_uM/Ki_AMP_NET +
                     ATX_uM/Ki_ATX_NET + VLX_uM/Ki_VLX_NET);

dxdt_DAT_occ = 2.0 * (DAT_occ_ss - DAT_occ);  // fast equilibrium
dxdt_NET_occ = 2.0 * (NET_occ_ss - NET_occ);

// ─── Synaptic DA dynamics ─────────────────────────────────
// Drug effects:
//   MPH/AMP block DAT → reduce kel_DA proportionally to occupancy
//   AMP also causes active efflux → extra DA release
double AMP_efflux_factor = 1.0 + 2.0 * AMP_uM / (AMP_uM + Ki_AMP_DAT);
double kel_DA_eff = kel_DA * (1.0 - 0.85 * DAT_occ);  // max 85% block
double krel_DA_eff = krel_DA * AMP_efflux_factor;
dxdt_DA_syn = krel_DA_eff * DA_base - kel_DA_eff * DA_syn;

// ─── Synaptic NE dynamics ─────────────────────────────────
// α2A-AR agonism (GFN) reduces presynaptic NE release
double GFN_A2A_occ = GFN_uM / (GFN_uM + Ki_GFN_A2A);
double krel_NE_eff = krel_NE * (1.0 - 0.6 * GFN_A2A_occ);  // autoreceptor
double kel_NE_eff  = kel_NE * (1.0 - 0.85 * NET_occ);
dxdt_NE_syn = krel_NE_eff * NE_base - kel_NE_eff * NE_syn;

// ─── PFC DA/NE tone (Inverted-U) ─────────────────────────
// Gaussian optimal: PFC_DA_tone = exp(-(DA-DA_opt)^2/(2*DA_bw^2))
double PFC_DA_ss = exp(-pow(DA_syn - DA_opt, 2.0) / (2.0 * DA_bw * DA_bw));
double PFC_NE_ss = exp(-pow(NE_syn - NE_opt, 2.0) / (2.0 * NE_bw * NE_bw));
// Guanfacine directly activates α2A post-synaptic → boosts PFC NE effect
double GFN_post = 1.0 + 0.5 * GFN_A2A_occ;
double PFC_NE_adj = fmin(PFC_NE_ss * GFN_post, 1.0);
dxdt_PFC_DA = 0.5 * (PFC_DA_ss - PFC_DA);
dxdt_PFC_NE = 0.5 * (PFC_NE_adj - PFC_NE);

// ─── Working Memory & Executive Function ─────────────────
double WM_ss  = 0.35 + 0.35 * PFC_DA + 0.30 * PFC_NE;  // baseline 0.35 ADHD
double EF_ss  = 0.30 + 0.30 * PFC_DA + 0.40 * PFC_NE;
dxdt_WM_idx  = 0.3 * (WM_ss  - WM_idx);
dxdt_ExecFun = 0.3 * (EF_ss  - ExecFun);

// ─── ADHD-RS-5 score dynamics (0-54; lower=better) ───────
// Score driven by WM + ExecFun improvement from baseline
double symptom_relief = (WM_idx - 0.35)/0.65 * k_WM_effect +
                        (ExecFun - 0.30)/0.70 * (1.0 - k_WM_effect);
double ADHD_RS_ss = ADHD_RS_base * (1.0 - symptom_relief);
dxdt_ADHD_RS = k_symp * (ADHD_RS_ss - ADHD_RS);

// ─── CGI-Severity (1-7; lower=better) ────────────────────
double CGI_S_ss = 4.0 - 2.5 * symptom_relief;  // 4=moderate, 1.5=min
if(CGI_S_ss < 1.0) CGI_S_ss = 1.0;
dxdt_CGI_S = 0.05 * (CGI_S_ss - CGI_S);

// ─── Quality of Life (0-1; higher=better) ─────────────────
double QoL_ss = 0.40 + 0.55 * symptom_relief;
if(QoL_ss > 1.0) QoL_ss = 1.0;
dxdt_QoL_idx = 0.05 * (QoL_ss - QoL_idx);

$TABLE
double Cp_MPH = CENT1 / (V2c * WT/70.0);  // µg/mL
double Cp_AMP = CENT2 / (V3c * WT/70.0);
double Cp_ATX = CENT3 / (V4c * WT/70.0);
double Cp_GFN = CENT4 / (V5c * WT/70.0);
double Cp_VLX = CENT5 / (V6c * WT/70.0);
double Cp_MPH_nM = Cp_MPH * 1000.0 / 233.7;
double Cp_AMP_nM = Cp_AMP * 1000.0 / 135.2;
double Cp_ATX_nM = Cp_ATX * 1000.0 / 291.8;
double Cp_GFN_nM = Cp_GFN * 1000.0 / 246.7;
double Cp_VLX_nM = Cp_VLX * 1000.0 / 237.3;
double response_pct = 100.0 * (1.0 - ADHD_RS / ADHD_RS_base);

$CAPTURE Cp_MPH Cp_AMP Cp_ATX Cp_GFN Cp_VLX
         Cp_MPH_nM Cp_AMP_nM Cp_ATX_nM Cp_GFN_nM Cp_VLX_nM
         DA_syn NE_syn DAT_occ NET_occ
         PFC_DA PFC_NE WM_idx ExecFun
         ADHD_RS CGI_S QoL_idx response_pct
'

## ---- Compile ----------------------------------------------------------------
mod <- mcode("adhd_qsp", code)

## ---- Dosing Events ----------------------------------------------------------

## Helper: multiple-dose event table
make_ev <- function(dose_mg, cmt, ii = 24, addl = 83) {
  ev(amt = dose_mg, cmt = cmt, ii = ii, addl = addl, time = 0)
}

## ---- Treatment Scenarios (7 scenarios) --------------------------------------

## Scenario 1: Untreated ADHD (natural history, no drug)
ev_untreated <- ev(amt = 0, cmt = "GUT1", time = 0)

## Scenario 2: Methylphenidate IR 10 mg TID (q8h)
ev_MPH_low <- ev(amt = 10, cmt = "GUT1", ii = 8, addl = 251, time = 0,
                 F1 = 0.22)

## Scenario 3: Methylphenidate IR 20 mg BID → ER 36 mg QD
## (Concerta paradigm: ka=0.4 for ER, dose=36 mg)
ev_MPH_ER <- ev(amt = 36, cmt = "GUT1", ii = 24, addl = 83, time = 0) %>%
  param(ka1 = 0.4, F1 = 0.22)

## Scenario 4: Mixed Amphetamine Salts (Adderall XR) 20 mg QD
ev_AMP <- ev(amt = 20, cmt = "GUT2", ii = 24, addl = 83, time = 0)

## Scenario 5: Atomoxetine 80 mg QD (EM phenotype)
ev_ATX <- ev(amt = 80, cmt = "GUT3", ii = 24, addl = 83, time = 0)

## Scenario 6: Guanfacine ER 4 mg QD (ADHD >6yr)
ev_GFN <- ev(amt = 4, cmt = "GUT4", ii = 24, addl = 83, time = 0)

## Scenario 7: Viloxazine ER 400 mg QD
ev_VLX <- ev(amt = 400, cmt = "GUT5", ii = 24, addl = 83, time = 0)

## ---- Simulation settings ----------------------------------------------------
sim_time  <- seq(0, 84 * 24, by = 1)   # 84 days (12 weeks)
acute_time <- seq(0, 24, by = 0.25)    # 24h PK profile

## ---- PK simulation (first-day) -----------------------------------------------
pk_scenarios <- list(
  MPH_IR_10mg  = list(ev = make_ev(10,  "GUT1", ii = 8,  addl = 2),   label = "MPH IR 10mg TID"),
  MPH_ER_36mg  = list(ev = make_ev(36,  "GUT1", ii = 24, addl = 0) %>%
                         param(ka1 = 0.4),                              label = "MPH ER 36mg QD"),
  AMP_XR_20mg  = list(ev = make_ev(20,  "GUT2", ii = 24, addl = 0),    label = "AMP XR 20mg QD"),
  ATX_80mg     = list(ev = make_ev(80,  "GUT3", ii = 24, addl = 0),    label = "ATX 80mg QD"),
  GFN_ER_4mg   = list(ev = make_ev(4,   "GUT4", ii = 24, addl = 0),    label = "GFN ER 4mg QD"),
  VLX_ER_400mg = list(ev = make_ev(400, "GUT5", ii = 24, addl = 0),    label = "VLX ER 400mg QD")
)

run_pk_scenario <- function(sc) {
  mrgsim(mod, ev = sc$ev, delta = 0.25, end = 36) %>%
    as_tibble() %>%
    mutate(scenario = sc$label)
}

pk_results <- bind_rows(lapply(pk_scenarios, run_pk_scenario))

## ── Plot PK profiles ----------------------------------------------------------
plot_pk <- function(drug_col, ylab) {
  ggplot(pk_results, aes(x = time, y = .data[[drug_col]], color = scenario)) +
    geom_line(linewidth = 1) +
    labs(x = "Time (h)", y = ylab,
         title = paste("Day 1 PK:", ylab),
         color = "Treatment") +
    theme_bw() +
    theme(legend.position = "bottom")
}

## ---- Long-term simulation (12 weeks) ----------------------------------------
scenarios_12wk <- list(
  list(ev = ev_untreated,  lbl = "Untreated"),
  list(ev = ev_MPH_low,    lbl = "MPH IR 10mg TID"),
  list(ev = ev_MPH_ER,     lbl = "MPH ER 36mg QD"),
  list(ev = ev_AMP,        lbl = "AMP XR 20mg QD"),
  list(ev = ev_ATX,        lbl = "ATX 80mg QD"),
  list(ev = ev_GFN,        lbl = "GFN ER 4mg QD"),
  list(ev = ev_VLX,        lbl = "VLX ER 400mg QD")
)

run_12wk <- function(sc) {
  mrgsim(mod, ev = sc$ev, delta = 24, end = 84 * 24) %>%
    as_tibble() %>%
    mutate(scenario = sc$lbl, week = time / 168)
}

results_12wk <- bind_rows(lapply(scenarios_12wk, run_12wk))

## ── Plot 12-week ADHD-RS trajectories ----------------------------------------
p_adhd_rs <- ggplot(results_12wk, aes(x = week, y = ADHD_RS, color = scenario)) +
  geom_line(linewidth = 1.2) +
  scale_y_reverse(limits = c(35, 5)) +  # lower = better
  labs(x = "Week", y = "ADHD-RS-5 Total Score",
       title = "ADHD-RS-5 over 12 Weeks by Treatment",
       caption = "Lower score = fewer symptoms") +
  theme_bw() + theme(legend.position = "bottom")

## ── Response rates at 12 weeks ------------------------------------------------
response_12wk <- results_12wk %>%
  filter(near(week, 12, tol = 0.5)) %>%
  group_by(scenario) %>%
  summarise(
    ADHD_RS_mean  = mean(ADHD_RS),
    Response_pct  = mean(response_pct),
    CGI_S_mean    = mean(CGI_S),
    WM_pct_chg    = 100 * (mean(WM_idx) - 0.35) / 0.35,
    EF_pct_chg    = 100 * (mean(ExecFun) - 0.30) / 0.30,
    QoL_final     = mean(QoL_idx)
  )

cat("\n=== 12-Week Response Summary ===\n")
print(as.data.frame(response_12wk), digits = 3)

## ── Plot DA/NE trajectories ---------------------------------------------------
p_DA_NE <- results_12wk %>%
  select(week, scenario, DA_syn, NE_syn) %>%
  pivot_longer(c(DA_syn, NE_syn), names_to = "transmitter", values_to = "conc") %>%
  ggplot(aes(x = week, y = conc, color = scenario, linetype = transmitter)) +
  geom_line(linewidth = 1) +
  labs(x = "Week", y = "Synaptic Concentration (nM)",
       title = "Synaptic DA/NE Levels over 12 Weeks") +
  theme_bw() + theme(legend.position = "bottom")

## ── Inverted-U function illustration -----------------------------------------
DA_seq <- seq(0, 5, by = 0.05)
inv_U_DA <- data.frame(
  DA = DA_seq,
  PFC_function = exp(-(DA_seq - 2.0)^2 / (2 * 1.5^2))
)

p_invU <- ggplot(inv_U_DA, aes(x = DA, y = PFC_function)) +
  geom_line(color = "steelblue", linewidth = 1.5) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 2.0, linetype = "dashed", color = "green") +
  annotate("text", x = 1.0, y = 0.85, label = "ADHD\nbaseline", color = "red", size = 3) +
  annotate("text", x = 2.0, y = 0.95, label = "Optimal\ntone", color = "green", size = 3) +
  labs(x = "[DA] PFC (nM)", y = "PFC Executive Function Index",
       title = "Inverted-U: DA Tone vs PFC Function\n(Arnsten model)") +
  theme_bw()

## ── Pediatric vs Adult PK comparison ----------------------------------------
compare_age <- function(wt, age_label, dose = 18) {
  mrgsim(mod, ev = make_ev(dose, "GUT1", ii = 24, addl = 0),
         param(WT = wt, ka1 = 0.4),
         delta = 0.5, end = 48) %>%
    as_tibble() %>%
    mutate(group = age_label)
}

age_results <- bind_rows(
  compare_age(25, "Child (25 kg)", 18),
  compare_age(50, "Adolescent (50 kg)", 36),
  compare_age(70, "Adult (70 kg)", 36)
)

p_age <- ggplot(age_results, aes(x = time, y = Cp_MPH_nM, color = group)) +
  geom_line(linewidth = 1.2) +
  labs(x = "Time (h)", y = "MPH Plasma (nM)",
       title = "Age/Weight Effect on MPH PK\n(Allometric scaling)") +
  theme_bw()

## ---- Biomarker correlation output ------------------------------------------
cat("\n=== DA/NET Occupancy at 12 Weeks (trough) ===\n")
results_12wk %>%
  filter(near(week, 12, tol = 0.5)) %>%
  select(scenario, DAT_occ, NET_occ, PFC_DA, PFC_NE) %>%
  group_by(scenario) %>%
  summarise(across(everything(), mean)) %>%
  print()

cat("\n=== Clinical Benchmarks ===\n")
cat("MPH IR 10mg TID: ADHD-RS reduction ~10 pts (MTA trial)\n")
cat("AMP XR 20mg QD: ADHD-RS reduction ~12 pts (Biederman 2002)\n")
cat("ATX 80mg QD:    ADHD-RS reduction ~8 pts  (Michelson 2001)\n")
cat("GFN ER 4mg QD:  ADHD-RS reduction ~7 pts  (Sallee 2009)\n")
cat("VLX ER 400mg:   ADHD-RS reduction ~8 pts  (Nasser 2021)\n")
cat("\nModel output should approximate above benchmarks.\n")

## ---- Save plots (if ggplot2 available) --------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggsave("adhd_rs_trajectories.png", p_adhd_rs, width = 10, height = 6, dpi = 150)
  ggsave("adhd_inverted_u.png",      p_invU,    width = 7,  height = 5, dpi = 150)
  ggsave("adhd_da_ne_levels.png",    p_DA_NE,   width = 10, height = 6, dpi = 150)
  cat("\nPlots saved.\n")
}

cat("\nADHD QSP mrgsolve model ready.\n")
