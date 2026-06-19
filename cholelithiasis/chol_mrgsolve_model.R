## =============================================================================
## Cholelithiasis (Gallstone Disease) — QSP Model
## mrgsolve ODE-based PK/PD Simulation
## =============================================================================
## Reference parameters calibrated from:
##   - Bachrach WH & Hofmann AF (1982) Ursodeoxycholic acid in the treatment
##     of cholesterol cholelithiasis. Dig Dis Sci 27:737-761
##   - Paumgartner G & Beuers U (2002) Ursodeoxycholic acid in cholestatic
##     liver disease. Hepatology 36:525-531
##   - Jazrawi RP et al (1992) Gut 33:381-386
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ---- mrgsolve model code block ----------------------------------------------
code <- '
$PROB Cholelithiasis QSP Model - UDCA/Statin PK/PD + Gallstone Dynamics

$PARAM
// ---- UDCA PK Parameters ----
DOSE_UDCA   = 750,   // total daily UDCA dose (mg/day), split TID
BWT         = 70,    // body weight (kg)
ka_UDCA     = 0.80,  // absorption rate constant UDCA (h^-1)
F_UDCA      = 0.50,  // oral bioavailability UDCA (fraction)
Vd_UDCA     = 12.0,  // volume of distribution UDCA (L/kg) → 840 L for 70 kg
CL_UDCA     = 40.0,  // total clearance UDCA (L/h)
kEHC_UDCA   = 0.35,  // enterohepatic cycling rate constant (h^-1)
f_hep       = 0.70,  // hepatic first-pass extraction fraction
f_bile      = 0.80,  // fraction of hepatic UDCA secreted into bile

// ---- Statin PK Parameters (Simvastatin) ----
DOSE_STAT   = 0,     // statin daily dose (mg), 0 = off; 40 = typical
ka_STAT     = 0.60,  // absorption rate constant statin (h^-1)
F_STAT      = 0.05,  // oral bioavailability statin (first-pass ~95% extraction)
Vd_STAT     = 3.2,   // volume of distribution statin (L/kg)
CL_STAT     = 70.0,  // total clearance statin (L/h)
Km_STAT     = 0.008, // Km for HMGCR inhibition (mg/L)
Emax_STAT   = 0.85,  // maximum HMGCR inhibition by statin (fraction)

// ---- Ezetimibe PK Parameters ----
DOSE_EZET   = 0,     // ezetimibe daily dose (mg), 0 = off; 10 = typical
ka_EZET     = 0.50,
F_EZET      = 0.35,
Vd_EZET     = 4.0,
CL_EZET     = 8.0,
Km_EZET     = 0.006, // Km for NPC1L1 inhibition (mg/L)
Emax_EZET   = 0.80,  // maximum intestinal cholesterol absorption inhibition

// ---- Bile Acid Pool Dynamics ----
BA_synth0   = 0.52,  // baseline bile acid synthesis rate (g/day ≈ 0.022 g/h)
BA_pool0    = 3.5,   // baseline total BA pool (g)
kBA_EHC     = 0.42,  // BA enterohepatic cycling rate (h^-1; ~6-10 cycles/day)
kBA_fecal   = 0.014, // BA fecal loss rate constant (h^-1)
E_UDCA_BA   = 0.25,  // fractional increase in BA pool by UDCA (at full dose)
KD_FXR      = 8.0,   // FXR activation concentration for BA feedback (µmol/L)

// ---- Hepatic Cholesterol ----
CHOL_h0     = 15.0,  // baseline hepatic free cholesterol pool (mmol)
k_CHOL_syn  = 1.50,  // hepatic cholesterol synthesis rate (mmol/h)
k_CHOL_deg  = 0.10,  // hepatic cholesterol utilization/degradation (h^-1)
E_STAT_CHOL = 0.60,  // fractional reduction in cholesterol synthesis by statin
E_FXR_CHOL  = 0.20,  // FXR-mediated reduction in biliary cholesterol secretion

// ---- Biliary Composition & Saturation ----
k_CHOL_bil  = 0.080, // rate of biliary cholesterol secretion (mmol/h)
k_PL_bil    = 0.18,  // rate of biliary phospholipid secretion (mmol/h)
PL_bil0     = 12.0,  // baseline biliary PL (mmol/L in GB bile)
BA_bil0     = 35.0,  // baseline biliary BA concentration (mmol/L in GB bile)
CHOL_bil0   = 4.2,   // baseline biliary cholesterol (mmol/L)
// CSI = CHOL_bil / (0.1875*BA_bil + 0.1429*PL_bil)

// ---- Gallbladder Volume & Motility ----
GB_vol0     = 30.0,  // fasting gallbladder volume (mL)
GB_vol_min  = 5.0,   // minimum GB volume (emptied, mL)
k_GB_fill   = 0.025, // GB filling rate constant (h^-1)
CCK_peak    = 1.0,   // normalized peak CCK (post-prandial)
k_GB_empty  = 0.30,  // GB emptying rate constant (h^-1)

// ---- Gallstone Dynamics ----
CSI_thresh  = 1.05,  // CSI threshold for nucleation
k_nucleat   = 0.0005,// crystal nucleation rate constant (mL/h per CSI unit above threshold)
k_grow      = 0.012, // stone growth rate constant (mL/h)
k_dissol    = 0.0,   // stone dissolution rate constant (mL/h) — set by UDCA effect
E_UDCA_dis  = 0.025, // UDCA dissolution effect coefficient (/mmol/L UDCA in bile)
Stone_vol0  = 0.0,   // initial stone volume (mL), 0=prevention, >0=dissolution
Stone_max   = 5.0,   // maximum stone volume (mL, for Emax model of growth)

// ---- Inflammatory Markers ----
IL6_base    = 2.0,   // baseline IL-6 (pg/mL)
CRP_base    = 0.5,   // baseline CRP (mg/L)
k_IL6_prod  = 0.05,  // IL-6 production rate per stone volume unit
k_IL6_elim  = 0.15,  // IL-6 elimination rate (h^-1)
k_CRP_prod  = 0.30,  // CRP production stimulated by IL-6
k_CRP_elim  = 0.035, // CRP elimination rate (h^-1)

// ---- Simulation Flags ----
WLOSS       = 0,     // weight loss intervention (1=yes)
k_WL        = 0.002  // weight loss rate (fraction/h, ~1.5 kg/week)

$CMT
// PK compartments - UDCA
A_gut_UDCA   // gut (absorption) compartment for UDCA [mg]
A_plas_UDCA  // plasma UDCA [mg]
A_hep_UDCA   // hepatic UDCA [mg]
A_bile_UDCA  // biliary UDCA [mg]
A_gb_UDCA    // gallbladder UDCA [mg]

// PK compartments - Statin
A_gut_STAT   // gut statin [mg]
A_plas_STAT  // plasma statin [mg]

// PK compartments - Ezetimibe
A_gut_EZET   // gut ezetimibe [mg]
A_plas_EZET  // plasma ezetimibe [mg]

// Bile acid & cholesterol dynamics
BA_pool      // total bile acid pool [g]
CHOL_h       // hepatic free cholesterol [mmol]
CHOL_bil     // biliary cholesterol [mmol/L equivalent]
PL_bil       // biliary phospholipid [mmol/L equivalent]

// Gallbladder & stone
GB_vol       // gallbladder volume [mL]
Crystal_mass // cholesterol crystal mass [mg]
Stone_V      // gallstone volume [mL]

// Inflammatory markers
IL6          // IL-6 [pg/mL]
CRP_plas     // CRP [mg/L]

$MAIN
// Derived PK volumes
double Vd_UDCA_L = Vd_UDCA * BWT;
double Vd_STAT_L = Vd_STAT * BWT;
double Vd_EZET_L = Vd_EZET * BWT;

// UDCA concentration in plasma (mg/L)
double C_UDCA_plas = A_plas_UDCA / Vd_UDCA_L;

// UDCA concentration in bile (µmol/L, MW UDCA = 392.6)
double C_UDCA_bile = (A_bile_UDCA / 392.6) * 1000.0;  // µmol/L

// Statin plasma concentration (mg/L)
double C_STAT_plas = A_plas_STAT / Vd_STAT_L;

// Ezetimibe plasma concentration (mg/L)
double C_EZET_plas = A_plas_EZET / Vd_EZET_L;

// Statin effect on HMGCR (Hill function, inhibition)
double E_STAT = (Emax_STAT * C_STAT_plas) / (Km_STAT + C_STAT_plas);

// Ezetimibe effect on intestinal cholesterol absorption
double E_EZET = (Emax_EZET * C_EZET_plas) / (Km_EZET + C_EZET_plas);

// FXR activation fraction (by BA pool via enterohepatic-delivered BAs)
double BA_conc_portal = (BA_pool / 3.5) * 40.0;  // µmol/L proxy
double FXR_act = BA_conc_portal / (KD_FXR + BA_conc_portal);

// Weight during weight loss (if enabled)
double BWT_t = BWT * (1.0 - WLOSS * k_WL * SOLVERTIME);
if(BWT_t < BWT * 0.75) BWT_t = BWT * 0.75;

// Cholesterol synthesis rate (inhibited by statin and SREBP2 feedback)
double k_CHOL_syn_eff = k_CHOL_syn * (1.0 - E_STAT) * (1.0 + 0.3*(1.0 - CHOL_h/CHOL_h0));

// Biliary cholesterol secretion rate (reduced by FXR, UDCA effect in bile)
double C_UDCA_bile_norm = C_UDCA_bile / 500.0;  // normalize to typical biliary UDCA
double E_UDCA_CSI = E_FXR_CHOL + E_UDCA_dis * C_UDCA_bile_norm;
if(E_UDCA_CSI > 0.70) E_UDCA_CSI = 0.70;
double k_CHOL_bil_eff = k_CHOL_bil * (1.0 - E_UDCA_CSI) * (1.0 - E_STAT * 0.3);

// Compute CSI (Cholesterol Saturation Index)
// CSI = CHOL_bil / (0.1875 * BA_bil + 0.1429 * PL_bil) from Admirand-Small diagram
double BA_bil = BA_pool / 0.10;  // rough [BA]_bile from pool size
double CSI = CHOL_bil / (0.1875 * BA_bil + 0.1429 * PL_bil + 1e-6);

// Stone dissolution rate (UDCA-enhanced)
double k_dissol_eff = E_UDCA_dis * C_UDCA_bile_norm * 0.5;

// Crystal nucleation (occurs above CSI threshold)
double delta_CSI = (CSI > CSI_thresh) ? (CSI - CSI_thresh) : 0.0;
double nucleat_rate = k_nucleat * delta_CSI * GB_vol;

// Stone growth (sigmoidal inhibition by bile capacity)
double growth_rate = k_grow * delta_CSI * Stone_V * (1.0 - Stone_V / Stone_max);

// IL-6 production from stone-induced inflammation
double k_IL6_stim = k_IL6_prod * Stone_V * (Stone_V > 0.1 ? 1.0 : 0.0);

// Initial conditions
double Stone_V_IC = Stone_vol0;
if(NEWIND <= 1) {
    _INIT(Stone_V)      = Stone_vol0;
    _INIT(Crystal_mass) = Stone_vol0 * 100.0;  // mg per mL stone
    _INIT(BA_pool)      = BA_pool0;
    _INIT(CHOL_h)       = CHOL_h0;
    _INIT(CHOL_bil)     = CHOL_bil0;
    _INIT(PL_bil)       = PL_bil0;
    _INIT(GB_vol)       = GB_vol0;
    _INIT(IL6)          = IL6_base;
    _INIT(CRP_plas)     = CRP_base;
}

$ODE
// ---- UDCA PK ----
double dose_UDCA_per = DOSE_UDCA / 3.0;  // TID dosing handled via event table

// Gut → Plasma (with first-pass hepatic extraction)
dxdt_A_gut_UDCA  = -ka_UDCA * A_gut_UDCA;
double UDCA_absorbed = ka_UDCA * A_gut_UDCA;
double UDCA_systemic = UDCA_absorbed * (1.0 - f_hep);
double UDCA_liver_in = UDCA_absorbed * f_hep;

// Plasma ↔ Peripheral (simplified 1-cpt for plasma)
dxdt_A_plas_UDCA = UDCA_systemic - (CL_UDCA / Vd_UDCA_L) * A_plas_UDCA;

// Hepatic UDCA
double UDCA_bile_out = f_bile * kEHC_UDCA * A_hep_UDCA;
dxdt_A_hep_UDCA  = UDCA_liver_in - kEHC_UDCA * A_hep_UDCA;

// Biliary UDCA (in bile ducts + GB)
double UDCA_gb_in = 0.40 * UDCA_bile_out;  // fraction going to GB during fasting
dxdt_A_bile_UDCA = UDCA_bile_out - 0.50 * kEHC_UDCA * A_bile_UDCA;

// GB UDCA (concentration during fasting)
dxdt_A_gb_UDCA   = UDCA_gb_in - kEHC_UDCA * A_gb_UDCA;

// ---- Statin PK ----
dxdt_A_gut_STAT  = -ka_STAT * A_gut_STAT;
dxdt_A_plas_STAT = ka_STAT * A_gut_STAT * F_STAT - (CL_STAT / Vd_STAT_L) * A_plas_STAT;

// ---- Ezetimibe PK ----
dxdt_A_gut_EZET  = -ka_EZET * A_gut_EZET;
dxdt_A_plas_EZET = ka_EZET * A_gut_EZET * F_EZET - (CL_EZET / Vd_EZET_L) * A_plas_EZET;

// ---- Bile Acid Pool ----
// BA synthesis (suppressed by FXR feedback from UDCA/BAs)
double BA_syn_rate = (BA_synth0 / 24.0) * (1.0 - 0.50 * FXR_act) * (1.0 + E_UDCA_BA * C_UDCA_bile_norm);
// BA fecal loss
double BA_fecal = kBA_fecal * BA_pool;
// Net BA pool dynamics
dxdt_BA_pool = BA_syn_rate - BA_fecal;

// ---- Hepatic Cholesterol ----
// Synthesis inhibited by statin, LDLR-mediated uptake adds
double CHOL_uptake = 0.20 * (1.0 + E_STAT * 0.8);  // LDLR upregulated by statin
dxdt_CHOL_h = k_CHOL_syn_eff + CHOL_uptake - k_CHOL_deg * CHOL_h - k_CHOL_bil_eff;

// ---- Biliary Cholesterol & Phospholipids ----
dxdt_CHOL_bil = k_CHOL_bil_eff - 0.08 * CHOL_bil;  // secretion minus removal/dilution
dxdt_PL_bil   = k_PL_bil - 0.06 * PL_bil;           // PL secretion vs removal

// ---- Gallbladder Volume ----
// Fasting: fills slowly; post-prandial: CCK → empties
dxdt_GB_vol = k_GB_fill * (GB_vol0 - GB_vol) - k_GB_empty * CCK_peak * GB_vol;
if(GB_vol < GB_vol_min) dxdt_GB_vol = 0.0;

// ---- Crystal & Stone Dynamics ----
// Crystal nucleation and growth
dxdt_Crystal_mass = nucleat_rate * 50.0 + growth_rate * 80.0 - k_dissol_eff * Crystal_mass;
if(Crystal_mass < 0.0) dxdt_Crystal_mass = 0.0;

// Stone volume
dxdt_Stone_V = growth_rate - k_dissol_eff * Stone_V;
if(Stone_V < 0.0) dxdt_Stone_V = 0.0;

// ---- Inflammatory Markers ----
dxdt_IL6     = IL6_base * k_IL6_elim + k_IL6_stim - k_IL6_elim * IL6;
dxdt_CRP_plas = k_CRP_prod * IL6 / (IL6_base + 1.0) * CRP_base - k_CRP_elim * CRP_plas;

$TABLE
double CSI_out = CHOL_bil / (0.1875 * (BA_pool / 0.10) + 0.1429 * PL_bil + 1e-6);
double UDCA_plas_conc = A_plas_UDCA / (Vd_UDCA * BWT);   // mg/L
double UDCA_bile_conc_umol = (A_bile_UDCA / 392.6) * 1000.0; // µmol/L
double STAT_plas_conc = A_plas_STAT / (Vd_STAT * BWT);   // mg/L
double EZET_plas_conc = A_plas_EZET / (Vd_EZET * BWT);   // mg/L
double Stone_mm = pow(Stone_V * 6.0 / 3.14159, 1.0/3.0) * 10.0; // approx diameter mm (sphere)
double BA_pool_g = BA_pool;
double CHOL_sat_pct = CSI_out * 100.0;

capture CSI_out UDCA_plas_conc UDCA_bile_conc_umol STAT_plas_conc Stone_V
capture Stone_mm BA_pool_g CHOL_sat_pct IL6 CRP_plas CHOL_bil PL_bil CHOL_h
'

## ---- Compile model ----------------------------------------------------------
mod <- mcode("cholelithiasis_qsp", code)

## ---- Helper: create dosing events -------------------------------------------
make_events <- function(dose_UDCA = 750, dose_STAT = 0, dose_EZET = 0,
                        dur_days = 365, freq_UDCA = "TID") {
  per_UDCA <- if(freq_UDCA == "TID") dose_UDCA / 3 else dose_UDCA
  int_UDCA <- if(freq_UDCA == "TID") 8 else 24   # hours between doses

  ev_list <- list()

  if(dose_UDCA > 0)
    ev_list$udca  <- ev(amt = per_UDCA, cmt = "A_gut_UDCA",
                        ii = int_UDCA, addl = dur_days * (24 / int_UDCA) - 1)
  if(dose_STAT > 0)
    ev_list$stat  <- ev(amt = dose_STAT,   cmt = "A_gut_STAT",
                        ii = 24, addl = dur_days - 1)
  if(dose_EZET > 0)
    ev_list$ezet  <- ev(amt = dose_EZET,   cmt = "A_gut_EZET",
                        ii = 24, addl = dur_days - 1)

  if(length(ev_list) == 0) return(ev(time = 0, amt = 0, cmt = 1))
  Reduce(ev_seq, ev_list)
}

## ---- Treatment Scenarios -----------------------------------------------------
scenarios <- list(
  "Scenario 1: No Treatment (Natural History)" = list(
    DOSE_UDCA = 0, DOSE_STAT = 0, DOSE_EZET = 0,
    Stone_vol0 = 0.5, WLOSS = 0
  ),
  "Scenario 2: UDCA 750 mg/day (Standard)" = list(
    DOSE_UDCA = 750, DOSE_STAT = 0, DOSE_EZET = 0,
    Stone_vol0 = 0.5, WLOSS = 0
  ),
  "Scenario 3: UDCA 1050 mg/day (High Dose)" = list(
    DOSE_UDCA = 1050, DOSE_STAT = 0, DOSE_EZET = 0,
    Stone_vol0 = 0.5, WLOSS = 0
  ),
  "Scenario 4: UDCA + Simvastatin 40 mg" = list(
    DOSE_UDCA = 750, DOSE_STAT = 40, DOSE_EZET = 0,
    Stone_vol0 = 0.5, WLOSS = 0
  ),
  "Scenario 5: Ezetimibe 10 mg Prevention" = list(
    DOSE_UDCA = 0, DOSE_STAT = 0, DOSE_EZET = 10,
    Stone_vol0 = 0.0, WLOSS = 0
  ),
  "Scenario 6: Lifestyle (Weight Loss) + UDCA" = list(
    DOSE_UDCA = 750, DOSE_STAT = 0, DOSE_EZET = 0,
    Stone_vol0 = 0.5, WLOSS = 1
  )
)

## ---- Run simulations ---------------------------------------------------------
sim_duration  <- 365     # days
sim_delta     <- 1       # hourly resolution? No, daily output: every 24h
sim_end_h     <- sim_duration * 24
sim_times     <- seq(0, sim_end_h, by = 24)  # daily output

run_scenario <- function(sc_name, params) {
  ev_data <- make_events(
    dose_UDCA = params$DOSE_UDCA,
    dose_STAT = params$DOSE_STAT,
    dose_EZET = params$DOSE_EZET,
    dur_days  = sim_duration
  )

  out <- mod %>%
    param(DOSE_UDCA = params$DOSE_UDCA,
          DOSE_STAT  = params$DOSE_STAT,
          DOSE_EZET  = params$DOSE_EZET,
          Stone_vol0 = params$Stone_vol0,
          WLOSS      = params$WLOSS) %>%
    mrgsim(events = ev_data,
           end    = sim_end_h,
           delta  = 24,
           obsonly = TRUE) %>%
    as_tibble() %>%
    mutate(Scenario = sc_name,
           Day = time / 24)

  return(out)
}

results <- bind_rows(
  mapply(run_scenario, names(scenarios), scenarios, SIMPLIFY = FALSE)
)

## ---- Plot 1: Stone Volume Dissolution / Progression -------------------------
p1 <- results %>%
  ggplot(aes(x = Day, y = Stone_V, color = Scenario, linetype = Scenario)) +
  geom_line(linewidth = 1.1) +
  labs(title    = "Gallstone Volume Over Time",
       subtitle = "Scenarios 1–4: Stone Volume (mL) — Dissolution vs. Growth",
       x = "Time (days)", y = "Stone Volume (mL)",
       color = "Treatment", linetype = "Treatment") +
  scale_color_manual(values = c("#D32F2F","#1565C0","#2E7D32","#6A1B9A",
                                "#F57F17","#00695C")) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(nrow = 3))

## ---- Plot 2: Cholesterol Saturation Index -----------------------------------
p2 <- results %>%
  ggplot(aes(x = Day, y = CHOL_sat_pct, color = Scenario)) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = 20, y = 102, label = "CSI = 1.0 (Saturation threshold)",
           color = "red", size = 3.5) +
  geom_line(linewidth = 1.0) +
  labs(title    = "Biliary Cholesterol Saturation Index (%)",
       subtitle = "CSI > 100% = lithogenic bile",
       x = "Time (days)", y = "CSI (%)") +
  scale_color_manual(values = c("#D32F2F","#1565C0","#2E7D32","#6A1B9A",
                                "#F57F17","#00695C")) +
  theme_classic(base_size = 13)

## ---- Plot 3: Bile Acid Pool Dynamics ----------------------------------------
p3 <- results %>%
  ggplot(aes(x = Day, y = BA_pool_g, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  labs(title = "Total Bile Acid Pool (g)",
       x = "Time (days)", y = "BA Pool (g)") +
  scale_color_manual(values = c("#D32F2F","#1565C0","#2E7D32","#6A1B9A",
                                "#F57F17","#00695C")) +
  theme_classic(base_size = 13)

## ---- Plot 4: UDCA Biliary Concentration (Scenarios with UDCA) ---------------
p4 <- results %>%
  filter(grepl("UDCA", Scenario)) %>%
  ggplot(aes(x = Day, y = UDCA_bile_conc_umol, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  labs(title    = "UDCA Biliary Concentration",
       subtitle = "Steady-state biliary UDCA (µmol/L)",
       x = "Time (days)", y = "UDCA in Bile (µmol/L)") +
  theme_classic(base_size = 13)

## ---- Plot 5: Inflammatory Markers -------------------------------------------
p5 <- results %>%
  ggplot(aes(x = Day, y = CRP_plas, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  labs(title = "Plasma CRP Over Time",
       x = "Time (days)", y = "CRP (mg/L)") +
  scale_color_manual(values = c("#D32F2F","#1565C0","#2E7D32","#6A1B9A",
                                "#F57F17","#00695C")) +
  theme_classic(base_size = 13)

## ---- Plot 6: Hepatic Cholesterol Over Time ----------------------------------
p6 <- results %>%
  ggplot(aes(x = Day, y = CHOL_h, color = Scenario)) +
  geom_line(linewidth = 1.0) +
  labs(title = "Hepatic Cholesterol Pool (mmol)",
       x = "Time (days)", y = "Hepatic Cholesterol (mmol)") +
  scale_color_manual(values = c("#D32F2F","#1565C0","#2E7D32","#6A1B9A",
                                "#F57F17","#00695C")) +
  theme_classic(base_size = 13)

## ---- Summary Table: Key Endpoints at 6 months and 12 months ----------------
endpoint_summary <- results %>%
  filter(Day %in% c(0, 180, 365)) %>%
  group_by(Scenario, Day) %>%
  summarise(
    Stone_Vol_mL    = round(mean(Stone_V), 3),
    Stone_Diam_mm   = round(mean(Stone_mm), 1),
    CSI_pct         = round(mean(CHOL_sat_pct), 1),
    BA_pool_g       = round(mean(BA_pool_g), 2),
    UDCA_bile_uM    = round(mean(UDCA_bile_conc_umol), 0),
    CRP_mgL         = round(mean(CRP_plas), 2),
    .groups = "drop"
  ) %>%
  arrange(Scenario, Day)

print(endpoint_summary)

## ---- Sensitivity Analysis: UDCA dose vs dissolution rate at Day 365 --------
udca_doses <- seq(250, 1500, by = 250)
dose_response <- lapply(udca_doses, function(d) {
  ev_d <- make_events(dose_UDCA = d, dur_days = 365)
  out  <- mod %>%
    param(DOSE_UDCA = d, Stone_vol0 = 0.5) %>%
    mrgsim(events = ev_d, end = 365*24, delta = 24, obsonly = TRUE) %>%
    as_tibble() %>%
    filter(time == 365*24) %>%
    mutate(Dose_UDCA = d,
           Pct_dissolv = (1 - Stone_V / 0.5) * 100)
  out
}) %>% bind_rows()

p_dose_resp <- dose_response %>%
  ggplot(aes(x = Dose_UDCA, y = Pct_dissolv)) +
  geom_point(size = 3, color = "#1565C0") +
  geom_line(color = "#1565C0", linewidth = 1.1) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "gray50") +
  labs(title    = "UDCA Dose-Response: Stone Dissolution at 12 Months",
       subtitle = "% stone dissolved from initial volume (0.5 mL)",
       x = "UDCA Daily Dose (mg/day)", y = "Stone Dissolution (%)") +
  theme_classic(base_size = 13)

cat("\n=== Cholelithiasis QSP Model — Simulation Complete ===\n")
cat("Scenarios simulated:", length(scenarios), "\n")
cat("Simulation duration:", sim_duration, "days\n")
cat("\nPlots generated: p1 (stone volume), p2 (CSI), p3 (BA pool),\n")
cat("  p4 (UDCA bile conc), p5 (CRP), p6 (hepatic CHOL), p_dose_resp\n")
cat("\nEndpoint summary table printed above.\n")
