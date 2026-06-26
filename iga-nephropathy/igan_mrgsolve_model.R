##########################################################################
# IgA Nephropathy (IgAN) – QSP mrgsolve Model
# Four-Hit Hypothesis + Drug PK/PD
#
# Disease State Variables:
#   GdIgA1   – Circulating Gd-IgA1 (normalized; healthy = 1, IgAN ≈ 1.5)
#   AutoAb   – Anti-Gd-IgA1 IgG autoantibodies (normalized)
#   IC_mes   – Mesangial immune-complex deposit index
#   CompAP   – Alternative-pathway complement activity (normalized)
#   Mesangial– Mesangial activation index
#   Podocyte – Podocyte integrity (1 = normal, 0 = fully depleted)
#   TIF      – Tubulointerstitial fibrosis index (0–1)
#   UPCR     – Urine protein:creatinine ratio (g/g)
#   eGFR     – Estimated GFR (mL/min/1.73 m²)
#   BP_sys   – Systolic blood pressure (mmHg)
#
# Drug compartments:
#   Budesonide TRF (Nefecon) · Sparsentan · Iptacopan · Sibeprenlimab
#
# References:
#   Barratt J et al. Kidney Int 2023 (NefIgArd)
#   Heerspink HJL et al. Lancet 2023 (PROTECT)
#   Rovin BH et al. NEJM 2024 (APPLAUSE-IgAN)
#   Barratt J et al. NEJM 2023 (AFFINITY)
#   Suzuki H et al. JASN 2011 (four-hit review)
##########################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)

# =========================================================================
# MODEL CODE
# =========================================================================
igan_code <- '
$PROB
IgA Nephropathy (IgAN) QSP Model – Four-Hit Hypothesis
Drug classes: Budesonide TRF | Sparsentan | Iptacopan | Sibeprenlimab | RAAS inh.

$PARAM @annotated
// ---------- Gd-IgA1 turnover ----------
k_syn_IgA1  : 0.05780  : /day | Gd-IgA1 synthesis (t1/2 ~ 12 d)
k_deg_IgA1  : 0.05780  : /day | Gd-IgA1 degradation

// ---------- Autoantibody (IgG) ----------
k_syn_AB    : 0.01980  : /day | anti-Gd-IgA1 IgG synthesis (GdIgA1-driven)
k_deg_AB    : 0.01980  : /day | anti-Gd-IgA1 IgG degradation (t1/2 ~ 35 d)

// ---------- IC mesangial deposition ----------
k_form_IC   : 0.10000  : /day | IC formation rate (GdIgA1 × AutoAb product)
k_clear_IC  : 0.12000  : /day | mesangial IC clearance

// ---------- Complement alternative pathway ----------
k_syn_CP    : 0.09000  : /day | AP complement activation (IC-driven)
k_deg_CP    : 0.18000  : /day | complement inactivation / regulators

// ---------- Mesangial activation ----------
k_act_MES   : 0.07000  : /day | mesangial activation (CompAP- and IC-driven)
k_res_MES   : 0.04000  : /day | mesangial return-to-baseline rate

// ---------- Podocyte integrity ----------
k_inj_Pod   : 0.00450  : /day | podocyte injury rate (MAC + inflammatory)
k_rep_Pod   : 0.00180  : /day | podocyte repair rate

// ---------- Tubulointerstitial fibrosis ----------
k_syn_TIF   : 0.00280  : /day | TIF progression (proteinuria and hypoxia)
k_deg_TIF   : 0.00040  : /day | TIF regression (very slow)

// ---------- UPCR dynamics ----------
k_syn_UPCR  : 2.00000  : g/g/day | proteinuria synthesis scaled to pod injury
k_deg_UPCR  : 0.28000  : /day    | proteinuria clearance dynamics

// ---------- eGFR ----------
k_loss_GFR  : 0.00018  : /day    | eGFR loss per TIF unit
k_RAAS_GFR  : 0.00014  : /day    | eGFR loss from intraglom. hypertension
eGFR_0      : 75       : mL/min  | initial eGFR at presentation

// ---------- BP homeostasis ----------
BP_0        : 135      : mmHg    | baseline systolic BP
k_BP_on     : 0.04500  : /day    | BP activation by AngII-like drive
k_BP_off    : 0.04500  : /day    | BP return-to-set-point rate

// ---------- Budesonide TRF PK (Nefecon 16 mg/day oral) ----------
ka_BUD      : 0.50000  : /h   | gut absorption rate constant
CL_BUD      : 90.00    : L/h  | apparent oral clearance (high hepatic extraction)
V_BUD       : 200.00   : L    | apparent distribution volume
F_BUD       : 0.15     : frac | oral bioavailability
Emax_BUD    : 0.55     : frac | max Gd-IgA1 reduction (budesonide)
EC50_BUD    : 8.00     : ng/mL| EC50 for Gd-IgA1 suppression
Hill_BUD    : 1.50     :      | Hill coefficient

// ---------- Sparsentan PK (400 mg/day oral) ----------
ka_SPA      : 0.40000  : /h   | gut absorption rate
CL_SPA      : 14.00    : L/h  | apparent oral clearance
V_SPA       : 120.00   : L    | apparent Vd
F_SPA       : 0.85     : frac | oral bioavailability
Emax_SPA    : 0.62     : frac | max UPCR reduction (dual AT1R+ETB)
EC50_SPA    : 200.00   : ng/mL| EC50 for proteinuria reduction
Hill_SPA    : 1.20     :      | Hill coefficient

// ---------- Iptacopan PK (200 mg BID oral) ----------
ka_IPT      : 0.80000  : /h   | gut absorption
CL_IPT      : 8.00     : L/h  | apparent oral clearance
V_IPT       : 60.00    : L    | apparent Vd
F_IPT       : 0.90     : frac | oral bioavailability
Emax_IPT    : 0.88     : frac | max AP complement reduction
EC50_IPT    : 50.00    : ng/mL| EC50 for factor B inhibition
Hill_IPT    : 2.00     :      | Hill coefficient

// ---------- Sibeprenlimab PK (500 mg SC Q4W) ----------
ka_SIB      : 0.01200  : /h   | SC absorption rate
CL_SIB      : 0.00800  : L/h  | mAb clearance
V1_SIB      : 3.50     : L    | central volume
V2_SIB      : 2.50     : L    | peripheral volume
Q_SIB       : 0.02000  : L/h  | intercompartmental clearance
F_SIB       : 0.80     : frac | SC bioavailability
Emax_SIB    : 0.68     : frac | max Gd-IgA1 reduction (anti-APRIL)
EC50_SIB    : 5.00     : ug/mL| EC50 (APRIL neutralization)
Hill_SIB    : 1.00     :      | Hill coefficient

// ---------- RAAS inhibitor (binary switch) ----------
E_RAAS      : 0.0      :      | RAAS inhibitor on/off (0=off, 1=full ACEi/ARB)
Emax_RAAS_P : 0.35     : frac | max UPCR reduction by RAAS blockade
Emax_RAAS_G : 0.22     : frac | max eGFR protection by RAAS

$CMT @annotated
BUD_gut     : Budesonide gut (mg)
BUD_central : Budesonide central (mg)
SPA_gut     : Sparsentan gut (mg)
SPA_central : Sparsentan central (mg)
IPT_gut     : Iptacopan gut (mg)
IPT_central : Iptacopan central (mg)
SIB_depot   : Sibeprenlimab SC depot (mg)
SIB_central : Sibeprenlimab central (mg)
SIB_periph  : Sibeprenlimab peripheral (mg)
GdIgA1      : Gd-IgA1 level (normalized)
AutoAb      : Anti-Gd-IgA1 IgG (normalized)
IC_mes      : Mesangial IC deposit (normalized)
CompAP      : Complement AP activity (normalized)
Mesangial   : Mesangial activation index
Podocyte    : Podocyte integrity (0–1)
TIF         : Tubulointerstitial fibrosis (0–1)
UPCR        : Urine P:Cr ratio (g/g)
eGFR        : eGFR (mL/min/1.73 m2)
BP_sys      : Systolic BP (mmHg)

$MAIN
// Concentrations
double C_BUD = BUD_central / V_BUD;      // ng/mL (dose mg → conc ng/mL after unit alignment)
double C_SPA = SPA_central / V_SPA;      // ng/mL
double C_IPT = IPT_central / V_IPT;      // ng/mL
double C_SIB = SIB_central / V1_SIB;    // ug/mL (mg / L = mg/L; 1 mg/L = 1 ug/mL)

// Hill-function drug effects
double E_bud = Emax_BUD * pow(C_BUD, Hill_BUD) /
               (pow(EC50_BUD, Hill_BUD) + pow(C_BUD, Hill_BUD) + 1e-12);
double E_spa = Emax_SPA * pow(C_SPA, Hill_SPA) /
               (pow(EC50_SPA, Hill_SPA) + pow(C_SPA, Hill_SPA) + 1e-12);
double E_ipt = Emax_IPT * pow(C_IPT, Hill_IPT) /
               (pow(EC50_IPT, Hill_IPT) + pow(C_IPT, Hill_IPT) + 1e-12);
double E_sib = Emax_SIB * pow(C_SIB, Hill_SIB) /
               (pow(EC50_SIB, Hill_SIB) + pow(C_SIB, Hill_SIB) + 1e-12);

// Combined mucosal effect on Gd-IgA1 (budesonide + sibeprenlimab independent)
double E_mucosal = 1.0 - (1.0 - E_bud) * (1.0 - E_sib);

// RAAS-derived effects
double E_raas_p = Emax_RAAS_P * E_RAAS;
double E_raas_g = Emax_RAAS_G * E_RAAS;

// Bounded states (defensive clipping)
double pod_ok  = (Podocyte > 0.0) ? Podocyte : 0.0;
double tif_cap = (TIF < 1.0) ? TIF : 1.0;
double egfr_ok = (eGFR > 5.0) ? eGFR : 5.0;

// Injury driver: complement MAC + inflammatory mesangial signals
double inj_drv = CompAP * 0.55 + Mesangial * 0.45;

$ODE
// ---- Budesonide TRF PK (oral, 1-compartment) ----
dxdt_BUD_gut     = -ka_BUD * BUD_gut;
dxdt_BUD_central = ka_BUD * F_BUD * BUD_gut - (CL_BUD / V_BUD) * BUD_central;

// ---- Sparsentan PK (oral, 1-compartment) ----
dxdt_SPA_gut     = -ka_SPA * SPA_gut;
dxdt_SPA_central = ka_SPA * F_SPA * SPA_gut - (CL_SPA / V_SPA) * SPA_central;

// ---- Iptacopan PK (oral, 1-compartment) ----
dxdt_IPT_gut     = -ka_IPT * IPT_gut;
dxdt_IPT_central = ka_IPT * F_IPT * IPT_gut - (CL_IPT / V_IPT) * IPT_central;

// ---- Sibeprenlimab PK (SC, 2-compartment) ----
dxdt_SIB_depot   = -ka_SIB * SIB_depot;
dxdt_SIB_central = ka_SIB * F_SIB * SIB_depot
                   - (CL_SIB / V1_SIB) * SIB_central
                   - (Q_SIB  / V1_SIB) * SIB_central
                   + (Q_SIB  / V2_SIB) * SIB_periph;
dxdt_SIB_periph  = (Q_SIB / V1_SIB) * SIB_central
                   - (Q_SIB / V2_SIB) * SIB_periph;

// ================================================================
// DISEASE PHARMACODYNAMIC ODEs
// ================================================================

// Hit 1 – Gd-IgA1 turnover
//   Synthesis suppressed by budesonide TRF (mucosal) and sibeprenlimab (APRIL)
dxdt_GdIgA1 = k_syn_IgA1 * (1.0 - E_mucosal) - k_deg_IgA1 * GdIgA1;

// Hit 2 – Autoantibody (GdIgA1-driven BCR stimulation)
dxdt_AutoAb = k_syn_AB * GdIgA1 - k_deg_AB * AutoAb;

// Hit 3 – IC mesangial deposition
dxdt_IC_mes = k_form_IC * GdIgA1 * AutoAb - k_clear_IC * IC_mes;

// Complement alternative pathway (IC-driven; blocked by iptacopan)
dxdt_CompAP = k_syn_CP * IC_mes * (1.0 - E_ipt) - k_deg_CP * CompAP;

// Mesangial activation (CompAP + direct IC via FcαRI)
dxdt_Mesangial = k_act_MES * (CompAP + IC_mes * 0.30) - k_res_MES * Mesangial;

// Podocyte integrity (bounded 0→1)
dxdt_Podocyte = -k_inj_Pod * inj_drv * pod_ok
                + k_rep_Pod * (1.0 - pod_ok);

// Tubulointerstitial fibrosis (driven by proteinuria + pod depletion)
dxdt_TIF = k_syn_TIF * (1.0 - pod_ok) * (1.0 + UPCR / 3.0) * (1.0 - tif_cap)
           - k_deg_TIF * tif_cap;

// UPCR dynamics
//   Production: glomerular injury × intraglomerular-pressure factor
//   Suppression: sparsentan (AT1R+ETB) and RAAS inhibitors
double UPCR_syn = k_syn_UPCR
                  * (1.0 - pod_ok)
                  * (1.0 + Mesangial * 0.40)
                  * (1.0 - E_spa)
                  * (1.0 - E_raas_p);
dxdt_UPCR = UPCR_syn - k_deg_UPCR * UPCR;

// eGFR decline
//   Loss from TIF accumulation + podocyte-driven glomerulosclerosis
//   RAAS inhibition provides partial eGFR protection
double GFR_loss = (k_loss_GFR * tif_cap
                   + k_RAAS_GFR * (1.0 - E_raas_g) * (1.0 - pod_ok));
dxdt_eGFR = -GFR_loss * egfr_ok;

// Systolic BP (AngII-driven equilibrium model)
double BP_driver = 1.0 + Mesangial * 0.12 - E_raas_p * 0.45;
dxdt_BP_sys = k_BP_on * (BP_driver * BP_0 - BP_sys)
              - k_BP_off * (BP_sys - BP_0);

$TABLE
double C_BUD_obs  = BUD_central / V_BUD;
double C_SPA_obs  = SPA_central / V_SPA;
double C_IPT_obs  = IPT_central / V_IPT;
double C_SIB_obs  = SIB_central / V1_SIB;
double UPCR_obs   = UPCR;
double eGFR_obs   = eGFR;
double GdIgA1_obs = GdIgA1;
double CompAP_obs = CompAP;
double Pod_obs    = Podocyte;
double TIF_obs    = TIF;
double IC_obs     = IC_mes;
double BP_obs     = BP_sys;
double AutoAb_obs = AutoAb;
// Derived: % UPCR reduction from baseline (assuming UPCR_base = 2.5)
double UPCR_pct_chg = 100.0 * (UPCR_obs - 2.5) / 2.5;
// CR50 flag: ≥50% reduction from 2.5 g/g baseline
double CR50 = (UPCR_obs <= 1.25) ? 1.0 : 0.0;
// eGFR slope (annualized)
double eGFR_slope = -(k_loss_GFR * TIF_obs
                       + k_RAAS_GFR * (1.0 - Emax_RAAS_G * E_RAAS)
                       * (1.0 - Podocyte)) * eGFR_obs * 365.0;

$CAPTURE C_BUD_obs C_SPA_obs C_IPT_obs C_SIB_obs
         UPCR_obs eGFR_obs GdIgA1_obs CompAP_obs
         Pod_obs TIF_obs IC_obs BP_obs AutoAb_obs
         UPCR_pct_chg CR50 eGFR_slope
'

# =========================================================================
# Compile
# =========================================================================
mod <- mcode("IgAN_QSP", igan_code)

# =========================================================================
# Initial conditions (typical high-risk IgAN patient at biopsy)
# =========================================================================
init_vals <- list(
  BUD_gut = 0, BUD_central = 0,
  SPA_gut = 0, SPA_central = 0,
  IPT_gut = 0, IPT_central = 0,
  SIB_depot = 0, SIB_central = 0, SIB_periph = 0,
  GdIgA1   = 1.55,   # ~55% above healthy (Moldoveanu 2007)
  AutoAb   = 1.30,   # elevated at presentation
  IC_mes   = 0.85,   # moderate mesangial IC deposits
  CompAP   = 0.65,   # active AP complement
  Mesangial = 0.75,  # moderate mesangial activation
  Podocyte  = 0.82,  # 82% integrity – some foot-process effacement
  TIF       = 0.12,  # early TIF (Oxford T0)
  UPCR      = 2.50,  # g/g – high-risk criterion (PROTECT/NefIgArd entry)
  eGFR      = 75.0,  # CKD stage G2
  BP_sys    = 136.0  # mild hypertension
)

mod <- mod %>% init(init_vals)

SIM_DAYS <- 730   # 2 years

# =========================================================================
# Dosing events
# =========================================================================
# Budesonide TRF: 16 mg once daily
e_BUD <- ev(amt = 16,  ii = 24, addl = SIM_DAYS - 1, cmt = "BUD_gut")
# Sparsentan: 400 mg once daily
e_SPA <- ev(amt = 400, ii = 24, addl = SIM_DAYS - 1, cmt = "SPA_gut")
# Iptacopan: 200 mg twice daily
e_IPT <- ev(amt = 200, ii = 12, addl = 2 * SIM_DAYS - 1, cmt = "IPT_gut")
# Sibeprenlimab: 500 mg SC every 28 days (Q4W)
e_SIB <- ev(amt = 500, ii = 28 * 24, addl = 25, cmt = "SIB_depot")

# =========================================================================
# 7 Treatment Scenarios
# =========================================================================
run_scenario <- function(mod, evs = NULL, raas = 0, label) {
  m <- mod %>% param(E_RAAS = raas)
  if (is.null(evs)) {
    out <- m %>% mrgsim(end = SIM_DAYS, delta = 1)
  } else {
    out <- m %>% mrgsim(events = evs, end = SIM_DAYS, delta = 1)
  }
  as_tibble(out) %>% mutate(Scenario = label)
}

scenarios <- bind_rows(
  run_scenario(mod, NULL,                  raas = 0, label = "1. Untreated"),
  run_scenario(mod, NULL,                  raas = 1, label = "2. RAAS inhibitor"),
  run_scenario(mod, e_BUD,                 raas = 1, label = "3. RAAS + Budesonide TRF"),
  run_scenario(mod, e_SPA,                 raas = 0, label = "4. Sparsentan"),
  run_scenario(mod, e_IPT,                 raas = 1, label = "5. RAAS + Iptacopan"),
  run_scenario(mod, e_SIB,                 raas = 1, label = "6. RAAS + Sibeprenlimab"),
  run_scenario(mod, c(e_BUD, e_SPA, e_IPT),raas = 0, label = "7. Triple Combo (BUD+SPA+IPT)")
) %>%
  mutate(
    time_wk  = time / 7,
    time_mo  = time / 30.44
  )

pal <- c(
  "1. Untreated"                    = "#B71C1C",
  "2. RAAS inhibitor"               = "#E65100",
  "3. RAAS + Budesonide TRF"        = "#1B5E20",
  "4. Sparsentan"                   = "#0D47A1",
  "5. RAAS + Iptacopan"             = "#4A148C",
  "6. RAAS + Sibeprenlimab"         = "#006064",
  "7. Triple Combo (BUD+SPA+IPT)"   = "#212121"
)

# =========================================================================
# Plots
# =========================================================================
theme_qsp <- theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 8),
    plot.title      = element_text(face = "bold", size = 13)
  )

p_upcr <- ggplot(scenarios, aes(x = time_mo, y = UPCR_obs, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = pal) +
  geom_hline(yintercept = 1.25, linetype = "dashed", color = "grey40") +
  annotate("text", x = 1, y = 1.15, label = "CR50 threshold (1.25 g/g)", size = 3, hjust = 0) +
  labs(title    = "IgAN: Proteinuria (UPCR) Trajectory",
       subtitle = "Primary efficacy endpoint; 2-year simulation",
       x = "Time (months)", y = "UPCR (g/g)") +
  guides(color = guide_legend(nrow = 3)) +
  theme_qsp

p_egfr <- ggplot(scenarios, aes(x = time_mo, y = eGFR_obs, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = pal) +
  geom_hline(yintercept = c(60, 30), linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.5, y = 61.5, label = "CKD G3a (60)", size = 3, hjust = 0) +
  annotate("text", x = 0.5, y = 31.5, label = "CKD G4 (30)", size = 3, hjust = 0) +
  labs(title = "IgAN: eGFR Trajectory",
       x = "Time (months)", y = "eGFR (mL/min/1.73 m²)") +
  guides(color = guide_legend(nrow = 3)) +
  theme_qsp

p_gdiga1 <- ggplot(scenarios, aes(x = time_mo, y = GdIgA1_obs, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = pal) +
  labs(title = "IgAN: Serum Gd-IgA1 (Hit 1 Biomarker)",
       x = "Time (months)", y = "Gd-IgA1 (normalized)") +
  guides(color = guide_legend(nrow = 3)) +
  theme_qsp

p_comp <- ggplot(scenarios, aes(x = time_mo, y = CompAP_obs, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = pal) +
  labs(title = "IgAN: Complement AP Activity",
       x = "Time (months)", y = "CompAP (normalized)") +
  guides(color = guide_legend(nrow = 3)) +
  theme_qsp

p_pod <- ggplot(scenarios, aes(x = time_mo, y = Pod_obs, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = pal) +
  labs(title = "IgAN: Podocyte Integrity",
       x = "Time (months)", y = "Podocyte integrity (0–1)") +
  guides(color = guide_legend(nrow = 3)) +
  theme_qsp

p_tif <- ggplot(scenarios, aes(x = time_mo, y = TIF_obs, color = Scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = pal) +
  labs(title = "IgAN: Tubulointerstitial Fibrosis Index",
       x = "Time (months)", y = "TIF index (0–1)") +
  guides(color = guide_legend(nrow = 3)) +
  theme_qsp

# PK profiles (first 14 days)
pk_day14 <- scenarios %>% filter(time <= 14)

p_pk_bud <- ggplot(
  filter(pk_day14, Scenario == "3. RAAS + Budesonide TRF"),
  aes(x = time, y = C_BUD_obs)) +
  geom_line(color = "#1B5E20", linewidth = 1.3) +
  labs(title = "Budesonide TRF PK (first 14 days)",
       x = "Time (days)", y = "Plasma conc. (ng/mL)") +
  theme_qsp

p_pk_spa <- ggplot(
  filter(pk_day14, Scenario == "4. Sparsentan"),
  aes(x = time, y = C_SPA_obs)) +
  geom_line(color = "#0D47A1", linewidth = 1.3) +
  labs(title = "Sparsentan PK (first 14 days)",
       x = "Time (days)", y = "Plasma conc. (ng/mL)") +
  theme_qsp

p_pk_sib <- ggplot(
  filter(scenarios, Scenario == "6. RAAS + Sibeprenlimab", time <= 120),
  aes(x = time, y = C_SIB_obs)) +
  geom_line(color = "#006064", linewidth = 1.3) +
  labs(title = "Sibeprenlimab PK (first 120 days, SC Q4W)",
       x = "Time (days)", y = "Serum conc. (μg/mL)") +
  theme_qsp

# =========================================================================
# Summary table at Week 36 (primary endpoint) and Week 104 (2 year)
# =========================================================================
summary_tbl <- scenarios %>%
  filter(time %in% c(252, 728)) %>%
  group_by(Scenario, time) %>%
  summarise(
    UPCR      = round(mean(UPCR_obs), 2),
    pct_UPCR  = round(mean(UPCR_pct_chg), 1),
    eGFR      = round(mean(eGFR_obs), 1),
    GdIgA1    = round(mean(GdIgA1_obs), 2),
    CompAP    = round(mean(CompAP_obs), 2),
    Podocyte  = round(mean(Pod_obs), 3),
    TIF       = round(mean(TIF_obs), 3),
    CR50_flag = round(mean(CR50), 0),
    .groups   = "drop"
  ) %>%
  mutate(Timepoint = ifelse(time == 252, "Wk 36", "Wk 104")) %>%
  select(Scenario, Timepoint, UPCR, pct_UPCR, eGFR, GdIgA1, CompAP, Podocyte, TIF, CR50_flag)

# =========================================================================
# Print & display
# =========================================================================
message("=== IgA Nephropathy QSP Model — Simulation Results ===")
print(summary_tbl, n = 50)

message("\nGenerating multi-panel figure...")
fig_main <- grid.arrange(p_upcr, p_egfr, p_gdiga1, p_comp, nrow = 2,
                         top = "IgA Nephropathy QSP Model — 2-Year Simulation")

fig_detail <- grid.arrange(p_pod, p_tif, p_pk_bud, p_pk_sib, nrow = 2,
                            top = "IgAN QSP — Podocyte, Fibrosis & Drug PK")

message("Done. Plots stored in fig_main and fig_detail.")
message("Objects: mod (compiled), scenarios (tidy tibble), summary_tbl (results)")

# =========================================================================
# Dose–response analysis for sparsentan (Week 36 UPCR)
# =========================================================================
spa_doses <- c(50, 100, 200, 400, 800)
dr_results <- lapply(spa_doses, function(d) {
  e_spa_dr <- ev(amt = d, ii = 24, addl = 251, cmt = "SPA_gut")
  mod %>%
    mrgsim(events = e_spa_dr, end = 252, delta = 1) %>%
    as_tibble() %>%
    filter(time == 252) %>%
    mutate(Dose_mg = d)
}) %>% bind_rows()

p_dr <- ggplot(dr_results, aes(x = Dose_mg, y = UPCR_obs)) +
  geom_line(color = "#0D47A1", linewidth = 1.3) +
  geom_point(size = 3, color = "#0D47A1") +
  scale_x_log10() +
  labs(title = "Sparsentan Dose–Response (UPCR at Week 36)",
       x = "Sparsentan dose (mg, log scale)", y = "UPCR at Wk 36 (g/g)") +
  theme_qsp

print(p_dr)
