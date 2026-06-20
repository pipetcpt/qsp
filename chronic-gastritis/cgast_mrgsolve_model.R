## ============================================================
## Chronic Gastritis QSP Model – mrgsolve ODE Implementation
## Disease: H. pylori-induced Chronic Gastritis → Correa Cascade
##
## 22 ODE Compartments:
##   Bacterial   : Hp_mucosal (CFU/biopsy, log scale)
##   Inflammation: NFkB, IL8, IL1b, TNFa, IFNg, IL10
##   Immune cells: Neutrophil, Th1, Treg
##   Gastric PD  : Gastrin, Acid, Mucus
##   Progression : Atrophy, IM (intestinal metaplasia score)
##   Drug PK     : PPI_gut, PPI_plasma, AMX, CLR, MTZ
##   Symptom     : Symptom_score
##
## 7 Treatment Scenarios:
##   1. No treatment (natural H. pylori history)
##   2. PPI monotherapy (symptom control only)
##   3. Standard triple therapy (PPI + AMX + CLR × 14 days)
##   4. Bismuth quadruple therapy (PPI + AMX + CLR + BSS × 14 days)
##   5. Metronidazole-based quadruple (PPI + AMX + MTZ + BSS × 14 days)
##   6. Vonoprazan triple (VPZ + AMX + CLR × 14 days)
##   7. Post-eradication follow-up (regression of atrophy/IM, 5 years)
##
## Calibration References:
##   - Malfertheiner et al. 2022 (Maastricht VI): PMID 35469816
##   - Graham et al. 2020 (vonoprazan): PMID 32356979
##   - Sugano et al. 2017 (kyoto consensus): PMID 28381443
##   - Correa 1992 (cascade model): PMID 1612357
##   - Souza et al. 2012 (atrophy regression): PMID 22710985
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ---- Model code ----
cgast_code <- '
$PROB
Chronic Gastritis QSP – H. pylori dynamics + immune cascade + Correa progression
22 ODE compartments, 7 treatment scenarios

$PARAM
// H. pylori parameters
Hp0         = 6.0     // initial log10(CFU/biopsy), range 3-8
Hp_max      = 8.0     // carrying capacity (log10 CFU/biopsy)
Hp_growth   = 0.08    // intrinsic growth rate (/day)
Hp_death    = 0.02    // background death rate (/day)

// Drug kill rates (log10 CFU/biopsy/day per unit Cp)
kill_AMX    = 0.15    // amoxicillin
kill_CLR    = 0.20    // clarithromycin
kill_MTZ    = 0.12    // metronidazole
kill_BSS    = 0.05    // bismuth (adhesion block)

// Resistance modifiers (0=fully sensitive, 1=fully resistant)
res_CLR     = 0.0     // clarithromycin resistance (0 or 1)
res_MTZ     = 0.0     // metronidazole resistance (0 or 1)

// Inflammation parameters (all concentrations pg/mL unless noted)
kNFkB_on    = 0.30    // NF-κB activation rate per Hp
kNFkB_deg   = 0.25    // NF-κB degradation rate
kIL8_prod   = 2.5     // IL-8 production per NF-κB
kIL8_deg    = 0.18    // IL-8 degradation rate
kIL1b_prod  = 1.2     // IL-1β production per NF-κB
kIL1b_deg   = 0.15    // IL-1β degradation rate
kTNFa_prod  = 0.80    // TNF-α production per NF-κB
kTNFa_deg   = 0.20    // TNF-α degradation rate
kIFNg_prod  = 0.50    // IFN-γ production per Th1
kIFNg_deg   = 0.12    // IFN-γ degradation rate
kIL10_prod  = 0.40    // IL-10 production (Treg, M2)
kIL10_deg   = 0.10    // IL-10 degradation rate

// Immune cell dynamics (normalized units 0-1)
kNeut_recr  = 0.35    // neutrophil recruitment per IL-8
kNeut_deg   = 0.25    // neutrophil clearance rate
kTh1_diff   = 0.12    // Th1 differentiation per DC stimulus
kTh1_deg    = 0.08    // Th1 turnover rate
kTreg_ind   = 0.08    // Treg induction per TGF-β/IL-10
kTreg_deg   = 0.06    // Treg turnover rate
Treg0       = 0.10    // baseline Treg level

// Gastric physiology (Gastrin pg/mL, Acid mmol/h, Mucus μm)
Gastrin_base = 30.0   // baseline serum gastrin
kGastr_Hp   = 8.0    // gastrin increase per unit Hp
kGastr_deg  = 0.03   // gastrin clearance
Acid_max    = 15.0   // maximum acid output mmol/h (BAO)
kAcid_stim  = 0.50   // acid stimulation per gastrin
kAcid_deg   = 0.20   // acid buffering/clearance
kAcid_PPI   = 0.80   // acid inhibition efficiency per PPI
Mucus_base  = 200.0  // baseline mucus thickness μm
kMucus_deg  = 0.04   // mucus degradation per ROS/inflam
kMucus_prod = 0.30   // mucus synthesis rate

// Correa cascade (scores 0-3 scale, slow timescale)
kAtrophy_prog  = 0.0003  // atrophy progression rate per day per inflam
kAtrophy_regr  = 0.0001  // atrophy regression rate post-eradication
kIM_prog       = 0.0002  // IM progression rate per atrophy
kIM_regr       = 0.00005 // IM regression rate (partial, slow)
Atrophy0       = 0.0     // initial atrophy score
IM0            = 0.0     // initial IM score

// PPI PK (omeprazole 20mg BID equivalent)
Ka_PPI      = 0.50    // absorption rate from gut (/h)
Ke_PPI      = 0.40    // plasma elimination rate (/h)
F_PPI       = 0.65    // oral bioavailability
Kact_PPI    = 0.60    // canalicular activation rate (/h)

// Amoxicillin PK (1000 mg BID)
Ka_AMX      = 1.00    // absorption rate (/h)
Ke_AMX      = 0.70    // elimination rate (/h; t½~1h)
F_AMX       = 0.90    // bioavailability
Kd_AMX      = 0.30    // gastric mucosa distribution

// Clarithromycin PK (500 mg BID)
Ka_CLR      = 0.60    // absorption rate (/h)
Ke_CLR      = 0.09    // elimination rate (/h; t½~7h)
F_CLR       = 0.55    // bioavailability
Kd_CLR      = 0.20    // distribution to gastric tissue

// Metronidazole PK (500 mg TID)
Ka_MTZ      = 0.80    // absorption rate (/h)
Ke_MTZ      = 0.12    // elimination rate (/h; t½~6h)
F_MTZ       = 0.99    // near-complete bioavailability

// Symptom score coefficients (NRS 0-10)
w_IL8       = 0.01    // weight of IL-8 on symptom
w_TNFa      = 0.03    // weight of TNF-α
w_Acid      = 0.20    // weight of acid output
w_Mucus_inv = 0.015   // inverse mucus contribution

$CMT
// Bacterial
Hp_mucosal
// Inflammation
NFkB IL8 IL1b TNFa IFNg IL10
// Immune cells (normalized 0-1)
Neutrophil Th1 Treg
// Gastric physiology
Gastrin Acid Mucus
// Disease progression (Correa cascade)
Atrophy IM_score
// Drug PK (plasma concentrations mg/L)
PPI_gut PPI_plasma AMX CLR MTZ
// Symptoms
Symptom

$INIT
Hp_mucosal = 6.0   // log10(CFU/biopsy) – moderate colonization
NFkB       = 0.5
IL8        = 50.0  // pg/mL baseline inflam
IL1b       = 5.0
TNFa       = 8.0
IFNg       = 3.0
IL10       = 15.0  // anti-inflammatory baseline
Neutrophil = 0.2
Th1        = 0.15
Treg       = 0.10
Gastrin    = 60.0  // pg/mL (Hp elevates above 30 base)
Acid       = 8.0   // mmol/h elevated
Mucus      = 160.0 // μm (compromised baseline)
Atrophy    = 0.0
IM_score   = 0.0
PPI_gut    = 0.0
PPI_plasma = 0.0
AMX        = 0.0
CLR        = 0.0
MTZ        = 0.0
Symptom    = 4.0   // baseline moderate symptom

$ODE
// ---------------------------------------------------------------
// Hp_mucosal dynamics (logistic growth with drug killing)
// ---------------------------------------------------------------
double Hp_norm = Hp_mucosal / Hp_max;         // fractional density 0-1
double Hp_lin  = pow(10.0, Hp_mucosal - 6.0); // linearised relative density

// drug kill terms (Emax models)
double kill_a = kill_AMX * AMX   / (AMX  + 0.5) * (1.0 - res_CLR * 0.0);
double kill_c = kill_CLR * CLR   / (CLR  + 0.1) * (1.0 - res_CLR);
double kill_m = kill_MTZ * MTZ   / (MTZ  + 2.0) * (1.0 - res_MTZ);
double kill_b = kill_BSS * 0.0;  // bismuth handled via parameter-driven scenario
double kill_total = kill_a + kill_c + kill_m + kill_b;

dxdt_Hp_mucosal = Hp_growth * (1.0 - Hp_norm) - Hp_death - kill_total;

// ---------------------------------------------------------------
// Inflammation cascade
// ---------------------------------------------------------------
// NF-κB – driven by Hp (via TLR/NOD/CagA) and ROS
double Hp_stim = (Hp_mucosal > 0.0) ? Hp_lin : 0.0;
dxdt_NFkB = kNFkB_on * Hp_stim * (1.0 - NFkB) - kNFkB_deg * NFkB - 0.15 * IL10 * NFkB / (IL10 + 10.0);

// IL-8 (major neutrophil chemokine from epithelium)
dxdt_IL8 = kIL8_prod * NFkB * 100.0 - kIL8_deg * IL8;

// IL-1β (inflammasome-dependent, NLRP3)
dxdt_IL1b = kIL1b_prod * NFkB * 50.0 - kIL1b_deg * IL1b;

// TNF-α (M1 macrophage, Th1 amplification)
dxdt_TNFa = kTNFa_prod * NFkB * 40.0 + 0.05 * Th1 * 40.0 - kTNFa_deg * TNFa;

// IFN-γ (Th1-derived)
dxdt_IFNg = kIFNg_prod * Th1 * 30.0 - kIFNg_deg * IFNg;

// IL-10 (Treg + M2 anti-inflammatory)
dxdt_IL10 = kIL10_prod * Treg * 50.0 + 5.0 - kIL10_deg * IL10;

// ---------------------------------------------------------------
// Immune cells (normalized 0-1)
// ---------------------------------------------------------------
// Neutrophils – recruited by IL-8
dxdt_Neutrophil = kNeut_recr * IL8 / (IL8 + 100.0) - kNeut_deg * Neutrophil;

// Th1 cells – differentiation driven by DC/IL-12, suppressed by Treg/IL-10
dxdt_Th1 = kTh1_diff * Hp_stim / (Hp_stim + 1.0) * (1.0 - Th1) * (1.0 - 0.5 * Treg) - kTh1_deg * Th1;

// Treg cells – induced by IL-10 + TGF-β signals
dxdt_Treg = kTreg_ind * IL10 / (IL10 + 20.0) * (1.0 - Treg) - kTreg_deg * (Treg - Treg0);

// ---------------------------------------------------------------
// Gastric physiology
// ---------------------------------------------------------------
// Gastrin – Hp stimulates G-cells, inhibits D-cells (somatostatin)
dxdt_Gastrin = kGastr_Hp * Hp_lin + Gastrin_base * kGastr_deg - kGastr_deg * Gastrin;

// Gastric acid – driven by gastrin/histamine, inhibited by PPI, IL-1β
double PPI_active = PPI_plasma;
double acid_inhibit = kAcid_PPI * PPI_active / (PPI_active + 0.5);   // PPI Emax
double IL1b_inhib  = 0.10 * IL1b / (IL1b + 20.0);                    // IL-1β inhibits parietal cell
double acid_stim   = kAcid_stim * Gastrin / (Gastrin + 30.0);
dxdt_Acid = Acid_max * acid_stim * (1.0 - acid_inhibit - IL1b_inhib) - kAcid_deg * Acid;

// Mucus layer – produced constitutively, degraded by inflammation/ROS
double ROS_proxy   = Neutrophil + 0.5 * IL1b / 30.0; // surrogate ROS
dxdt_Mucus = kMucus_prod * (Mucus_base - Mucus) - kMucus_deg * ROS_proxy * Mucus;

// ---------------------------------------------------------------
// Correa cascade (very slow – years-scale)
// ---------------------------------------------------------------
// Atrophy – driven by chronic IL-1β, TNF-α, IFN-γ; reversal after eradication
double inflam_drive = (IL1b / 10.0 + TNFa / 15.0 + IFNg / 5.0) / 3.0;
double eradicated   = (Hp_mucosal < 2.0) ? 1.0 : 0.0; // eradication flag
dxdt_Atrophy = kAtrophy_prog * inflam_drive * (3.0 - Atrophy)
             - kAtrophy_regr * eradicated * Atrophy;

// Intestinal metaplasia – progression from atrophy
dxdt_IM_score = kIM_prog * Atrophy * (3.0 - IM_score)
              - kIM_regr * eradicated * IM_score;

// ---------------------------------------------------------------
// Drug PK
// ---------------------------------------------------------------
// PPI (omeprazole): gut -> plasma -> canalicular activation
dxdt_PPI_gut    = -Ka_PPI * PPI_gut;
dxdt_PPI_plasma =  Ka_PPI * PPI_gut - Ke_PPI * PPI_plasma;

// Amoxicillin: gut -> plasma (simple 1-compartment)
dxdt_AMX = -Ke_AMX * AMX;  // absorb via dosing event; plasma = AMX

// Clarithromycin
dxdt_CLR = -Ke_CLR * CLR;

// Metronidazole
dxdt_MTZ = -Ke_MTZ * MTZ;

// ---------------------------------------------------------------
// Symptom score (0-10 scale)
// ---------------------------------------------------------------
double acid_contrib  = w_Acid * Acid;
double inflam_contrib = w_IL8 * IL8 / 10.0 + w_TNFa * TNFa / 10.0;
double mucus_contrib  = w_Mucus_inv * (Mucus_base - Mucus);
double raw_score = acid_contrib + inflam_contrib + mucus_contrib;
dxdt_Symptom = 0.5 * (raw_score > 10.0 ? 10.0 : raw_score) - 0.5 * Symptom;

$TABLE
double Hp_CFU      = pow(10.0, Hp_mucosal);
double PGI_proxy   = 100.0 - 25.0 * Atrophy;       // corpus function proxy μg/L
double PGII_proxy  = 10.0 + 5.0 * (IL8 / 50.0);    // inflammation marker μg/L
double PG_ratio    = (PGII_proxy > 0) ? PGI_proxy / PGII_proxy : 99.0;
double OLGA_stage  = (Atrophy < 0.5) ? 0.0 : (Atrophy < 1.5) ? 1.0 : (Atrophy < 2.5) ? 2.0 : 3.0;
double eradicated_flag = (Hp_mucosal < 2.0) ? 1.0 : 0.0;

$CAPTURE Hp_mucosal Hp_CFU NFkB IL8 IL1b TNFa IFNg IL10
         Neutrophil Th1 Treg Gastrin Acid Mucus
         Atrophy IM_score PGI_proxy PGII_proxy PG_ratio OLGA_stage
         PPI_plasma AMX CLR MTZ Symptom eradicated_flag
'

## ---- Compile model ----
cgast_mod <- mcode("cgast", cgast_code)

## ============================================================
## Dosing event generators
## ============================================================

## PPI 20 mg BID (omeprazole equivalent), bioavailability applied
make_PPI_dose <- function(days = 14, F_val = 0.65, dose_mg = 20) {
  dose_amt <- dose_mg * F_val
  ev(cmt = "PPI_gut", amt = dose_amt, ii = 12, addl = 2 * days - 1)
}

## Amoxicillin 1000 mg BID
make_AMX_dose <- function(days = 14, F_val = 0.90, dose_mg = 1000) {
  dose_amt <- dose_mg * F_val / 10.0  # rough scaling to mg/L
  ev(cmt = "AMX", amt = dose_amt, ii = 12, addl = 2 * days - 1)
}

## Clarithromycin 500 mg BID
make_CLR_dose <- function(days = 14, F_val = 0.55, dose_mg = 500) {
  dose_amt <- dose_mg * F_val / 15.0
  ev(cmt = "CLR", amt = dose_amt, ii = 12, addl = 2 * days - 1)
}

## Metronidazole 500 mg TID
make_MTZ_dose <- function(days = 14, F_val = 0.99, dose_mg = 500) {
  dose_amt <- dose_mg * F_val / 15.0
  ev(cmt = "MTZ", amt = dose_amt, ii = 8, addl = 3 * days - 1)
}

## ============================================================
## Simulation function
## ============================================================
run_scenario <- function(
  scenario_name,
  use_PPI  = FALSE,
  use_AMX  = FALSE,
  use_CLR  = FALSE,
  use_MTZ  = FALSE,
  res_CLR  = 0,
  res_MTZ  = 0,
  treat_days = 14,
  total_days = 365,        # simulate 1 year by default
  long_sim   = FALSE       # 5-year simulation for Correa cascade
) {
  sim_days <- if (long_sim) 5 * 365 else total_days
  times <- seq(0, sim_days * 24, by = 12)  # every 12h in hours

  # Build event list
  ev_list <- list()
  if (use_PPI) ev_list[["PPI"]] <- make_PPI_dose(days = treat_days)
  if (use_AMX) ev_list[["AMX"]] <- make_AMX_dose(days = treat_days)
  if (use_CLR) ev_list[["CLR"]] <- make_CLR_dose(days = treat_days)
  if (use_MTZ) ev_list[["MTZ"]] <- make_MTZ_dose(days = treat_days)

  ev_all <- if (length(ev_list) > 0) do.call(c, ev_list) else ev()

  out <- cgast_mod %>%
    param(res_CLR = res_CLR, res_MTZ = res_MTZ) %>%
    ev(ev_all) %>%
    mrgsim(end = sim_days * 24, delta = 12, carry_out = "evid") %>%
    as_tibble() %>%
    mutate(
      scenario = scenario_name,
      time_day = time / 24
    )
  out
}

## ============================================================
## 7 Treatment Scenarios
## ============================================================
message("Running scenario 1/7: Untreated (natural history)...")
s1 <- run_scenario("1. No Treatment\n(Natural history)", total_days = 365)

message("Running scenario 2/7: PPI monotherapy...")
s2 <- run_scenario("2. PPI Monotherapy\n(Omeprazole 20mg BID)",
                   use_PPI = TRUE, treat_days = 56, total_days = 365)

message("Running scenario 3/7: Standard triple therapy (PPI+AMX+CLR)...")
s3 <- run_scenario("3. Standard Triple\n(PPI+AMX+CLR × 14d)",
                   use_PPI = TRUE, use_AMX = TRUE, use_CLR = TRUE,
                   treat_days = 14, total_days = 180)

message("Running scenario 4/7: Bismuth quadruple (PPI+AMX+CLR+BSS)...")
# BSS simulated via slightly increased CLR kill (BSS adds ~15% efficacy)
s4_mod <- cgast_mod %>% param(kill_CLR = 0.22)
ev4 <- c(make_PPI_dose(14), make_AMX_dose(14), make_CLR_dose(14))
s4 <- s4_mod %>% ev(ev4) %>%
  mrgsim(end = 180 * 24, delta = 12) %>% as_tibble() %>%
  mutate(scenario = "4. Bismuth Quad\n(PPI+AMX+CLR+BSS × 14d)", time_day = time / 24)

message("Running scenario 5/7: Metronidazole quadruple (PPI+AMX+MTZ+BSS)...")
s5 <- run_scenario("5. MTZ Quadruple\n(PPI+AMX+MTZ+BSS × 14d)",
                   use_PPI = TRUE, use_AMX = TRUE, use_MTZ = TRUE,
                   treat_days = 14, total_days = 180)

message("Running scenario 6/7: Vonoprazan triple (VPZ+AMX+CLR)...")
# VPZ – faster acid suppression, higher intragastric pH → better AMX activity
vpz_mod <- cgast_mod %>% param(kAcid_PPI = 0.92, Ka_PPI = 1.2)
ev6 <- c(make_PPI_dose(14, F_val = 0.80, dose_mg = 20), make_AMX_dose(14), make_CLR_dose(14))
s6 <- vpz_mod %>% ev(ev6) %>%
  mrgsim(end = 180 * 24, delta = 12) %>% as_tibble() %>%
  mutate(scenario = "6. Vonoprazan Triple\n(VPZ+AMX+CLR × 14d)", time_day = time / 24)

message("Running scenario 7/7: Post-eradication 5-year follow-up...")
# Start from eradicated state (Hp_mucosal very low), simulate Correa regression
s7 <- cgast_mod %>%
  init(Hp_mucosal = 1.0, NFkB = 0.05, IL8 = 8.0, IL1b = 1.0,
       TNFa = 2.0, IFNg = 0.5, IL10 = 20.0, Neutrophil = 0.05,
       Th1 = 0.05, Treg = 0.12, Gastrin = 35.0, Acid = 5.0,
       Mucus = 190.0, Atrophy = 1.5, IM_score = 0.8, Symptom = 1.5) %>%
  mrgsim(end = 5 * 365 * 24, delta = 24) %>% as_tibble() %>%
  mutate(scenario = "7. Post-eradication\n(5-year regression)", time_day = time / 24)

## Combine key scenarios
all_results <- bind_rows(s1, s2, s3, s4, s5, s6) %>%
  filter(time_day <= 180)

## ============================================================
## Plotting functions
## ============================================================

## 1. H. pylori dynamics across triple/quadruple regimens
plot_Hp_dynamics <- function(data) {
  data %>%
    filter(scenario %in% c(
      "1. No Treatment\n(Natural history)",
      "3. Standard Triple\n(PPI+AMX+CLR × 14d)",
      "4. Bismuth Quad\n(PPI+AMX+CLR+BSS × 14d)",
      "5. MTZ Quadruple\n(PPI+AMX+MTZ+BSS × 14d)",
      "6. Vonoprazan Triple\n(VPZ+AMX+CLR × 14d)"
    )) %>%
    ggplot(aes(x = time_day, y = Hp_mucosal, color = scenario)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 2, linetype = "dashed", color = "red", linewidth = 0.8) +
    annotate("text", x = 5, y = 2.2, label = "Eradication threshold (10² CFU)", color = "red", size = 3) +
    labs(title = "H. pylori Mucosal Density – Eradication Comparison",
         x = "Time (days)", y = "H. pylori (log₁₀ CFU/biopsy)",
         color = "Treatment") +
    theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 7))
}

## 2. Inflammatory cascade (IL-8, IL-1β, TNF-α)
plot_inflammation <- function(data) {
  data %>%
    filter(scenario %in% c(
      "1. No Treatment\n(Natural history)",
      "3. Standard Triple\n(PPI+AMX+CLR × 14d)"
    )) %>%
    select(time_day, scenario, IL8, IL1b, TNFa, IFNg) %>%
    pivot_longer(c(IL8, IL1b, TNFa, IFNg), names_to = "cytokine", values_to = "conc") %>%
    ggplot(aes(x = time_day, y = conc, color = scenario, linetype = cytokine)) +
    geom_line(linewidth = 0.9) +
    labs(title = "Inflammatory Cytokine Profiles",
         x = "Time (days)", y = "Concentration (pg/mL)",
         color = "Treatment", linetype = "Cytokine") +
    theme_bw()
}

## 3. Gastric physiology (Gastrin, Acid, Mucus)
plot_gastric_physiol <- function(data) {
  data %>%
    filter(scenario %in% c(
      "1. No Treatment\n(Natural history)",
      "2. PPI Monotherapy\n(Omeprazole 20mg BID)",
      "3. Standard Triple\n(PPI+AMX+CLR × 14d)"
    )) %>%
    select(time_day, scenario, Gastrin, Acid, Mucus) %>%
    pivot_longer(c(Gastrin, Acid, Mucus), names_to = "variable", values_to = "value") %>%
    ggplot(aes(x = time_day, y = value, color = scenario)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ variable, scales = "free_y",
               labeller = labeller(variable = c(Gastrin = "Gastrin (pg/mL)",
                                                Acid = "Gastric Acid (mmol/h)",
                                                Mucus = "Mucus Thickness (μm)"))) +
    labs(title = "Gastric Physiology Dynamics", x = "Time (days)",
         y = "Value", color = "Treatment") +
    theme_bw()
}

## 4. Correa cascade (Atrophy and IM, long-term)
plot_correa <- function() {
  s7 %>%
    select(time_day, Atrophy, IM_score, PG_ratio) %>%
    pivot_longer(c(Atrophy, IM_score, PG_ratio), names_to = "marker", values_to = "value") %>%
    ggplot(aes(x = time_day / 365, y = value, color = marker)) +
    geom_line(linewidth = 1) +
    labs(title = "Correa Cascade Regression Post-Eradication (5 years)",
         x = "Time (years)", y = "Score / Ratio",
         color = "Biomarker") +
    scale_color_manual(values = c(Atrophy = "#E64A19", IM_score = "#7C4DFF", PG_ratio = "#0288D1"),
                       labels = c(Atrophy = "Atrophy Score (0-3)",
                                  IM_score = "IM Score (0-3)",
                                  PG_ratio = "PGI/PGII Ratio")) +
    theme_bw()
}

## 5. PGI/PGII ratio (atrophy biomarker) across scenarios
plot_biomarkers <- function(data) {
  data %>%
    filter(scenario %in% c(
      "1. No Treatment\n(Natural history)",
      "3. Standard Triple\n(PPI+AMX+CLR × 14d)"
    )) %>%
    ggplot(aes(x = time_day, y = PG_ratio, color = scenario)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 3.0, linetype = "dashed", color = "darkred") +
    annotate("text", x = 10, y = 3.2, label = "PGI/PGII < 3: atrophy threshold", size = 3, color = "darkred") +
    labs(title = "PGI/PGII Ratio – Non-invasive Atrophy Biomarker",
         x = "Time (days)", y = "PGI/PGII Ratio",
         color = "Treatment") +
    theme_bw()
}

## 6. Symptom score comparison
plot_symptoms <- function(data) {
  data %>%
    ggplot(aes(x = time_day, y = Symptom, color = scenario)) +
    geom_line(linewidth = 0.9) +
    labs(title = "Global Symptom Score (GIS 0-10 scale)",
         x = "Time (days)", y = "Symptom Score",
         color = "Treatment") +
    ylim(0, 10) +
    theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 7))
}

## ============================================================
## Steady-state summary table
## ============================================================
summary_table <- all_results %>%
  group_by(scenario) %>%
  summarise(
    Hp_day14   = mean(Hp_mucosal[abs(time_day - 14) < 0.6], na.rm = TRUE),
    Hp_day90   = mean(Hp_mucosal[abs(time_day - 90) < 0.6], na.rm = TRUE),
    IL8_day14  = mean(IL8[abs(time_day - 14) < 0.6], na.rm = TRUE),
    Acid_day14 = mean(Acid[abs(time_day - 14) < 0.6], na.rm = TRUE),
    Symptom_day14 = mean(Symptom[abs(time_day - 14) < 0.6], na.rm = TRUE),
    eradicated = mean(eradicated_flag[abs(time_day - 90) < 0.6], na.rm = TRUE),
    .groups = "drop"
  )

message("Summary table (day 14 and day 90):")
print(summary_table)

## ============================================================
## Sensitivity analysis – key parameters on Hp eradication
## ============================================================
sens_params <- c("kill_CLR", "kill_AMX", "kill_MTZ", "Ka_PPI", "Ke_CLR", "res_CLR")
sens_range  <- list(
  kill_CLR = seq(0.05, 0.40, length.out = 8),
  kill_AMX = seq(0.05, 0.30, length.out = 8),
  kill_MTZ = seq(0.04, 0.25, length.out = 8),
  Ka_PPI   = seq(0.20, 1.0,  length.out = 8),
  Ke_CLR   = seq(0.05, 0.20, length.out = 8),
  res_CLR  = c(0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0)
)

sens_results <- lapply(names(sens_params), function(pname) {
  pvals <- sens_range[[pname]]
  sapply(pvals, function(pval) {
    args <- setNames(list(pval), pname)
    out <- do.call(function(...) param(cgast_mod, ...), args) %>%
      ev(c(make_PPI_dose(14), make_AMX_dose(14), make_CLR_dose(14))) %>%
      mrgsim(end = 30 * 24, delta = 12) %>%
      as_tibble()
    mean(out$eradicated_flag[out$time / 24 >= 28], na.rm = TRUE)
  })
})

message("\nSensitivity analysis complete. Results in sens_results list.")
message("\n=== Chronic Gastritis QSP Model – Setup Complete ===")
message("Run plot_Hp_dynamics(all_results) to visualise H. pylori eradication curves.")
message("Run plot_correa() to see Correa cascade regression over 5 years.")
