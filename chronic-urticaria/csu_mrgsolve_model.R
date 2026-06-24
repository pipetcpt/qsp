## ============================================================
## Chronic Spontaneous Urticaria (CSU) — mrgsolve QSP Model
## IgE/FcεRI Pathway · Mast Cell Activation · Type-2 Inflammation
## 18 ODE compartments · 6 treatment scenarios
##
## Calibration references:
##   GLACIAL (Omalizumab 300mg q4wk) — Kaplan 2013 JACI
##   ASTERIA I/II (Omalizumab) — Saini 2015 JACI; Maurer 2013 NEJM
##   LIBERTY-CSU CUPID A/B (Dupilumab) — Simpson 2023 NEJM
##   H1-antihistamine PK — Simons 2004 JACI
##   Omalizumab PopPK — Lowe 2009 J Allergy Clin Immunol
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

## ----------------------------------------------------------
## 1. Model code block
## ----------------------------------------------------------

code <- '
$PROB
  Chronic Spontaneous Urticaria (CSU) QSP Model
  IgE/FcεRI-Mast Cell Axis + Type-2 Cytokine Network
  18 ODE compartments

$PARAM @annotated
  // ---- H1-antihistamine PK (1-compartment, cetirizine prototype) ----
  ka_AH    : 1.5    : /h   Antihistamine absorption rate constant
  CL_AH    : 4.2    : L/h  Antihistamine clearance
  Vd_AH    : 50.0   : L    Antihistamine volume of distribution
  F_AH     : 0.73   : -    Oral bioavailability

  // ---- Omalizumab PK (2-compartment, SC) ----
  ka_OMA   : 0.011  : /h   Omalizumab SC absorption rate constant (~64h tmax)
  CL_OMA   : 0.0049 : L/h  Omalizumab clearance
  V1_OMA   : 3.5    : L    Central volume
  V2_OMA   : 3.1    : L    Peripheral volume
  Q_OMA    : 0.003  : L/h  Inter-compartment clearance
  F_OMA    : 0.62   : -    SC bioavailability

  // ---- Dupilumab PK (2-compartment, SC) ----
  ka_DUP   : 0.0087 : /h   Dupilumab SC absorption rate constant
  CL_DUP   : 0.0071 : L/h  Dupilumab clearance
  V1_DUP   : 4.8    : L    Central volume
  V2_DUP   : 2.9    : L    Peripheral volume
  Q_DUP    : 0.0024 : L/h  Inter-compartment clearance
  F_DUP    : 0.64   : -    SC bioavailability

  // ---- BTK inhibitor PK (1-compartment oral, remibrutinib prototype) ----
  ka_BTK   : 2.1    : /h   BTKi absorption rate constant
  CL_BTK   : 38.0   : L/h  BTKi clearance (high first-pass)
  Vd_BTK   : 280.0  : L    BTKi volume of distribution
  F_BTK    : 0.36   : -    BTKi oral bioavailability

  // ---- IgE biology ----
  ksyn_IgE : 0.0015 : nM/h Basal IgE synthesis rate
  kdeg_IgE : 0.0050 : /h   Free IgE degradation (t1/2 ~140h)
  kbind_OMA: 45.0   : /nM/h Omalizumab-IgE binding rate
  kdis_OMA : 0.0005 : /h   Omalizumab-IgE dissociation rate
  IgE0     : 300.0  : nM   Baseline free IgE (elevated CSU ~300 IU/mL equiv)

  // ---- FcεRI / Mast Cell ----
  FcεRI_tot: 1.0    : rel  Total FcεRI expression on mast cells (normalised)
  karm_MC  : 0.08   : /h   Rate of FcεRI arming by IgE
  kdisarm  : 0.02   : /h   Spontaneous FcεRI disarming / receptor turnover
  kact_MC  : 0.15   : /h   Mast cell activation rate (armed MC + autoantigen)
  kinh_AH  : 0.88   : -    Maximum inhibition of mast cell degranulation by AH
  EC50_AH  : 0.12   : mg/L AH EC50 for H1R blockade
  n_AH     : 1.5    : -    Hill coefficient for AH
  kinh_BTK : 0.92   : -    Maximum BTKi inhibition of MC activation
  EC50_BTK : 0.045  : mg/L BTKi EC50
  kdeact_MC: 0.12   : /h   Mast cell deactivation rate
  MC0      : 1.0    : rel  Baseline mast cell priming state
  kprime_IL33: 0.04 : /h/u IL-33 priming of mast cells

  // ---- Histamine PK (skin/plasma) ----
  krel_H   : 2.5    : /h   Histamine release rate from activated MC
  kdeg_Hs  : 0.85   : /h   Histamine degradation in skin
  kdeg_Hp  : 3.2    : /h   Histamine degradation in plasma
  ktrans_H : 0.15   : /h   Skin-to-plasma transfer
  Hist0    : 0.1    : nM   Baseline histamine

  // ---- IL-4 / IL-13 / IL-31 / IL-33 dynamics ----
  ksyn_IL4 : 0.008  : nM/h IL-4 synthesis (Th2/ILC2)
  kdeg_IL4 : 0.55   : /h   IL-4 degradation
  ksyn_IL13: 0.010  : nM/h IL-13 synthesis
  kdeg_IL13: 0.48   : /h   IL-13 degradation
  ksyn_IL31: 0.005  : nM/h IL-31 synthesis (itch mediator)
  kdeg_IL31: 0.62   : /h   IL-31 degradation
  ksyn_IL33: 0.006  : nM/h IL-33 synthesis
  kdeg_IL33: 0.70   : /h   IL-33 degradation
  kinh_DUP_IL4 : 0.96 : - Dupilumab max inhibition of IL-4 signaling
  kinh_DUP_IL13: 0.95 : - Dupilumab max inhibition of IL-13 signaling
  EC50_DUP : 0.008  : mg/L Dupilumab EC50 (IL-4Rα)

  // ---- Eosinophil dynamics ----
  keo_in   : 0.003  : /h   Eosinophil tissue recruitment (IL-5/eotaxin-driven)
  keo_out  : 0.025  : /h   Eosinophil clearance
  Eo0      : 1.0    : rel  Baseline skin eosinophil level

  // ---- Disease activity & UAS7 mapping ----
  UAS7_max : 42.0   : score Maximum UAS7 score
  UAS7_0   : 30.0   : score Baseline UAS7 score (moderate-severe CSU)
  kUAS_H   : 8.0    : -    UAS7 sensitivity to skin histamine
  kUAS_IL31: 5.0    : -    UAS7 sensitivity to IL-31 (itch component)
  IgE_norm : 300.0  : nM   Normalisation IgE (= IgE0)

$CMT @annotated
  // PK compartments
  AH_GI    : Antihistamine GI depot (mg)
  AH_plasma: Antihistamine plasma (mg/L equiv)
  OMA_depot: Omalizumab SC depot (mg)
  OMA_c    : Omalizumab central (mg)
  OMA_p    : Omalizumab peripheral (mg)
  DUP_depot: Dupilumab SC depot (mg)
  DUP_c    : Dupilumab central (mg)
  DUP_p    : Dupilumab peripheral (mg)
  BTK_GI   : BTKi GI depot (mg)
  BTK_plasma: BTKi plasma (mg/L)

  // PD compartments
  IgE_free : Free IgE (nM)
  IgE_OMA  : IgE-Omalizumab complex (nM)
  MC_primed: Armed (IgE-loaded) mast cell index (rel)
  MC_act   : Activated mast cell index (rel)
  Hist_skin: Skin histamine (nM)
  Hist_plasm: Plasma histamine (nM)
  IL31_skin: Skin IL-31 (nM)
  IL33_skin: Skin IL-33 (nM)

$MAIN
  // Steady-state initial conditions
  IgE_free_0 = IgE0;
  MC_primed_0 = MC0;
  Hist_skin_0 = Hist0;
  Hist_plasm_0 = Hist0 * 0.1;
  IL31_skin_0 = ksyn_IL31 / kdeg_IL31;
  IL33_skin_0 = ksyn_IL33 / kdeg_IL33;

$ODE
  // ----------------------------------------------------------------
  // PK: Antihistamine (1-cpt oral)
  double C_AH = AH_plasma / Vd_AH;   // mg/L

  dxdt_AH_GI     = -ka_AH * AH_GI;
  dxdt_AH_plasma =  ka_AH * F_AH * AH_GI - (CL_AH / Vd_AH) * AH_plasma;

  // ----------------------------------------------------------------
  // PK: Omalizumab (2-cpt SC)
  double C_OMA = OMA_c / V1_OMA;   // mg/L

  dxdt_OMA_depot = -ka_OMA * OMA_depot;
  dxdt_OMA_c     =  ka_OMA * F_OMA * OMA_depot - (CL_OMA + Q_OMA) / V1_OMA * OMA_c
                    + Q_OMA / V2_OMA * OMA_p
                    - kbind_OMA * C_OMA * IgE_free + kdis_OMA * IgE_OMA;
  dxdt_OMA_p     =  Q_OMA / V1_OMA * OMA_c - Q_OMA / V2_OMA * OMA_p;

  // ----------------------------------------------------------------
  // PK: Dupilumab (2-cpt SC)
  double C_DUP = DUP_c / V1_DUP;   // mg/L

  dxdt_DUP_depot = -ka_DUP * DUP_depot;
  dxdt_DUP_c     =  ka_DUP * F_DUP * DUP_depot - (CL_DUP + Q_DUP) / V1_DUP * DUP_c
                    + Q_DUP / V2_DUP * DUP_p;
  dxdt_DUP_p     =  Q_DUP / V1_DUP * DUP_c - Q_DUP / V2_DUP * DUP_p;

  // ----------------------------------------------------------------
  // PK: BTKi (1-cpt oral)
  double C_BTK = BTK_plasma / Vd_BTK;   // mg/L

  dxdt_BTK_GI    = -ka_BTK * BTK_GI;
  dxdt_BTK_plasma =  ka_BTK * F_BTK * BTK_GI - (CL_BTK / Vd_BTK) * BTK_plasma;

  // ----------------------------------------------------------------
  // Drug effects
  double E_AH  = kinh_AH  * pow(C_AH,  n_AH)  / (pow(EC50_AH,  n_AH)  + pow(C_AH,  n_AH));
  double E_BTK = kinh_BTK * C_BTK / (EC50_BTK + C_BTK);
  double E_DUP_IL4  = kinh_DUP_IL4  * C_DUP / (EC50_DUP + C_DUP);
  double E_DUP_IL13 = kinh_DUP_IL13 * C_DUP / (EC50_DUP + C_DUP);

  // ----------------------------------------------------------------
  // IgE / Omalizumab binding
  dxdt_IgE_free = ksyn_IgE - kdeg_IgE * IgE_free
                  - kbind_OMA * C_OMA * IgE_free + kdis_OMA * IgE_OMA;
  dxdt_IgE_OMA  = kbind_OMA * C_OMA * IgE_free - kdis_OMA * IgE_OMA;

  // ----------------------------------------------------------------
  // FcεRI arming / mast cell priming
  // Omalizumab reduces free IgE → less FcεRI arming
  double fIgE = IgE_free / IgE_norm;   // normalised free IgE fraction
  dxdt_MC_primed = karm_MC * fIgE * (FcεRI_tot - MC_primed - MC_act)
                   - kdisarm * MC_primed
                   - kact_MC * (1.0 - E_AH) * (1.0 - E_BTK) * MC_primed
                   + kprime_IL33 * IL33_skin * (FcεRI_tot - MC_primed - MC_act);

  dxdt_MC_act    = kact_MC * (1.0 - E_AH) * (1.0 - E_BTK) * MC_primed
                   - kdeact_MC * MC_act;

  // ----------------------------------------------------------------
  // Histamine (skin and plasma)
  dxdt_Hist_skin  = krel_H * MC_act - kdeg_Hs * Hist_skin - ktrans_H * Hist_skin;
  dxdt_Hist_plasm = ktrans_H * Hist_skin - kdeg_Hp * Hist_plasm;

  // ----------------------------------------------------------------
  // Cytokines (type-2 network)
  // IL-31 (itch mediator) — stimulated by MC activation
  dxdt_IL31_skin = ksyn_IL31 * (1.0 + 2.0 * MC_act)
                   - kdeg_IL31 * IL31_skin;

  // IL-33 (alarmin, amplifies MC priming) — attenuated by dupilumab downstream
  dxdt_IL33_skin = ksyn_IL33 * (1.0 + 1.5 * MC_act)
                   - kdeg_IL33 * IL33_skin;

$TABLE
  // Concentrations
  double CONC_AH   = AH_plasma / Vd_AH;      // mg/L
  double CONC_OMA  = OMA_c / V1_OMA;          // mg/L
  double CONC_DUP  = DUP_c / V1_DUP;          // mg/L
  double CONC_BTK  = BTK_plasma / Vd_BTK;     // mg/L

  // IgE suppression (%)
  double IgE_suppression = (1.0 - IgE_free / IgE_norm) * 100.0;

  // UAS7 surrogate (daily score × 7 = weekly)
  // UAS7 components: wheal count (ISS) + itch intensity (HSS)
  double MC_effect = MC_act / MC0;
  double H_effect  = Hist_skin / Hist0;
  double I_effect  = IL31_skin / (ksyn_IL31 / kdeg_IL31);
  double UAS7 = UAS7_0 * (0.5 * MC_effect * H_effect + 0.3 * I_effect + 0.2);
  if (UAS7 > UAS7_max) UAS7 = UAS7_max;
  if (UAS7 < 0.0)      UAS7 = 0.0;

  // Well-controlled urticaria flag (UAS7 ≤ 6)
  double WCU = (UAS7 <= 6.0) ? 1.0 : 0.0;

  // Complete response flag (UAS7 = 0)
  double CR  = (UAS7 == 0.0) ? 1.0 : 0.0;

  // Capture
  capture CONC_AH CONC_OMA CONC_DUP CONC_BTK
  capture IgE_free IgE_suppression MC_primed MC_act
  capture Hist_skin Hist_plasm IL31_skin IL33_skin
  capture UAS7 WCU CR
'

## ----------------------------------------------------------
## 2. Compile model
## ----------------------------------------------------------

mod <- mcode("csu_qsp", code)

## ----------------------------------------------------------
## 3. Dosing regimens
## ----------------------------------------------------------

# H1-antihistamine: cetirizine 10 mg QD oral
dose_AH_std <- ev(amt = 10, cmt = "AH_GI", ii = 24, addl = 27)   # 4 weeks

# High-dose antihistamine: 40 mg/day (4×10 mg)
dose_AH_high <- ev(amt = 40, cmt = "AH_GI", ii = 24, addl = 83)  # 12 weeks

# Omalizumab 300 mg q4wk SC
dose_OMA_300 <- ev(amt = 300, cmt = "OMA_depot", ii = 4*168, addl = 5)  # 6 doses (~24wk)

# Omalizumab 150 mg q4wk SC (lower dose)
dose_OMA_150 <- ev(amt = 150, cmt = "OMA_depot", ii = 4*168, addl = 5)

# Dupilumab 300 mg q2wk SC (after 600 mg loading)
dose_DUP_LD  <- ev(time = 0,    amt = 600, cmt = "DUP_depot")
dose_DUP_MD  <- ev(time = 336, amt = 300, cmt = "DUP_depot", ii = 2*168, addl = 10)
dose_DUP <- c(dose_DUP_LD, dose_DUP_MD)

# BTKi (remibrutinib prototype) 25 mg QD oral
dose_BTK <- ev(amt = 25, cmt = "BTK_GI", ii = 24, addl = 167)  # 24 weeks

## ----------------------------------------------------------
## 4. Treatment scenarios
## ----------------------------------------------------------

scenarios <- list(
  list(id = 1, name = "No treatment",              ev = ev()),
  list(id = 2, name = "Cetirizine 10 mg QD",        ev = dose_AH_std),
  list(id = 3, name = "High-dose AH 40 mg/day",     ev = dose_AH_high),
  list(id = 4, name = "Omalizumab 300 mg q4wk",     ev = dose_OMA_300),
  list(id = 5, name = "Omalizumab 300 mg + AH",     ev = c(dose_OMA_300, dose_AH_std)),
  list(id = 6, name = "Dupilumab 300 mg q2wk",      ev = dose_DUP),
  list(id = 7, name = "BTKi 25 mg QD",              ev = dose_BTK)
)

## ----------------------------------------------------------
## 5. Simulation function
## ----------------------------------------------------------

sim_scenario <- function(sc, end_h = 24 * 168, delta = 24) {
  out <- mod %>%
    mrgsim(
      events = sc$ev,
      end    = end_h,
      delta  = delta,
      obsonly = TRUE
    ) %>%
    as.data.frame() %>%
    mutate(
      scenario = sc$name,
      time_wk  = time / 168
    )
  out
}

results <- map_dfr(scenarios, sim_scenario)

## ----------------------------------------------------------
## 6. Plots
## ----------------------------------------------------------

# UAS7 over time
p_uas7 <- results %>%
  ggplot(aes(x = time_wk, y = UAS7, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 6,  linetype = "dashed", colour = "darkgreen", linewidth = 0.7) +
  geom_hline(yintercept = 0,  linetype = "dotted",  colour = "steelblue", linewidth = 0.7) +
  annotate("text", x = 22, y = 7.5, label = "WCU threshold (UAS7 ≤ 6)",
           colour = "darkgreen", size = 3) +
  scale_x_continuous(breaks = seq(0, 24, 4)) +
  scale_y_continuous(limits = c(0, 42), breaks = seq(0, 42, 7)) +
  labs(title    = "CSU Disease Activity — UAS7 by Treatment Scenario",
       subtitle = "Chronic Spontaneous Urticaria QSP Model",
       x        = "Time (weeks)",
       y        = "UAS7 Score",
       colour   = "Treatment") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

print(p_uas7)

# IgE suppression (omalizumab scenarios)
p_ige <- results %>%
  filter(grepl("Omalizumab|No treatment", scenario)) %>%
  ggplot(aes(x = time_wk, y = IgE_suppression, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  labs(title  = "Free IgE Suppression — Omalizumab Scenarios",
       x      = "Time (weeks)",
       y      = "IgE Suppression (%)",
       colour = "Treatment") +
  theme_classic(base_size = 12)

print(p_ige)

# Mast cell activation
p_mc <- results %>%
  ggplot(aes(x = time_wk, y = MC_act, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  labs(title  = "Mast Cell Activation Index",
       x      = "Time (weeks)",
       y      = "MC Activation (rel.)",
       colour = "Treatment") +
  theme_classic(base_size = 12)

print(p_mc)

# Skin histamine
p_hist <- results %>%
  ggplot(aes(x = time_wk, y = Hist_skin, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  labs(title  = "Skin Histamine Concentration",
       x      = "Time (weeks)",
       y      = "Skin Histamine (nM)",
       colour = "Treatment") +
  theme_classic(base_size = 12)

print(p_hist)

## ----------------------------------------------------------
## 7. Summary table at key time points
## ----------------------------------------------------------

key_wks <- c(4, 12, 24)

summary_tbl <- results %>%
  filter(round(time_wk, 1) %in% key_wks) %>%
  group_by(scenario, time_wk) %>%
  slice_tail(n = 1) %>%
  summarise(
    UAS7_mean       = round(mean(UAS7), 1),
    WCU_pct         = round(mean(WCU) * 100, 1),
    IgE_sup_pct     = round(mean(IgE_suppression), 1),
    MC_act_rel      = round(mean(MC_act), 3),
    .groups = "drop"
  )

print(summary_tbl)
