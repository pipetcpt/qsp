# ==============================================================================
# Celiac Disease (CD) — Quantitative Systems Pharmacology Model
# Implemented in mrgsolve (R)
# ==============================================================================
# Disease:  Celiac Disease (Gluten-sensitive enteropathy)
# Version:  1.0  (2026-06-18)
# Author:   QSP Library (Claude Code Routine)
#
# Key References:
#   Shan et al. (2002) Science 297:2275-2279 [33-mer gliadin]
#   Sollid et al. (2012) Immunity 36:907-919 [HLA-DQ2/DQ8]
#   Jabri & Sollid (2006) Nat Clin Pract Gastroenterol 3:516-525
#   Mention et al. (2003) Gastroenterology 125:730-745 [IL-15]
#   Leffler et al. (2015) Gastroenterology 148:1311-1319 [Larazotide]
#   Schuppan et al. (2017) Aliment Pharmacol Ther 45:1145-1155 [ZED1227]
#   Mukherjee et al. (2017) PLoS ONE 12:e0172518 [IL-15/AMG714]
# ==============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ==============================================================================
# mrgsolve model code block
# ==============================================================================
cd_model_code <- '
$PROB Celiac Disease QSP Model — 21-Compartment ODE System
  Disease: Celiac Disease (CD)
  Pathways: Gluten→tTG2 deamidation, innate IL-15, Th1/Th17 adaptive immunity,
            B cell serology, intestinal histopathology (Marsh), nutritional sequelae
  Treatments: GFD, Larazotide (AT-1001), ZED1227 (TG2i), AMG714 (anti-IL-15)

$PARAM @annotated
// --- Gluten intake ---
GlutenIntake_g_day : 12  : Daily gluten intake on unrestricted Western diet (g/day)
GFD               : 0   : Gluten-free diet compliance (0=none, 1=strict)
GFD_leak          : 0   : GFD leakage fraction (0=no leak; 0.05=5% trace contamination)

// --- Gluten processing & deamidation ---
kabs_gliadin   : 0.25   : h-1, gliadin peptide absorption into lamina propria
kdeg_lumen     : 0.40   : h-1, luminal degradation of gliadin by proteases
k_deamid       : 0.30   : h-1, tTG2-mediated deamidation of gliadin peptides
kdeg_DGP       : 0.35   : h-1, clearance of deamidated gliadin peptides (DGP)

// --- tTG2 enzyme kinetics ---
tTG2_basal     : 1.0    : Baseline tTG2 activity (normalized, 1=normal)
k_tTG2_TGFb    : 0.05   : Rate of tTG2 upregulation by TGF-β / IFN-γ
k_tTG2_deg     : 0.04   : h-1, tTG2 turnover rate

// --- Intestinal permeability ---
IP_basal       : 1.0    : Baseline permeability index (normalized)
k_IP_IL15      : 0.08   : IP increase rate driven by IL-15 (zonulin pathway)
k_IP_restore   : 0.03   : h-1, permeability restoration rate (tight junction repair)
EC50_IP_IL15   : 2.5    : ng/mL, EC50 of IL-15 for zonulin/IP induction

// --- Innate IL-15 ---
k_IL15_prod    : 0.15   : ng/mL/h, IL-15 production by epithelial cells
k_IL15_deg     : 0.12   : h-1, IL-15 degradation rate (half-life ~5.8 h)
EC50_IL15_glia : 0.8    : AU, EC50 of gliadin for IL-15 production

// --- IEL dynamics ---
IEL_normal     : 20     : Normal IEL count (per 100 enterocytes; threshold=25)
k_IEL_IL15     : 0.06   : IL-15-driven IEL expansion rate
k_IEL_death    : 0.025  : h-1, IEL apoptosis / turnover rate
IEL_max        : 120    : Maximum IEL count per 100 enterocytes
EC50_IEL_IL15  : 3.0    : ng/mL, EC50 of IL-15 for IEL expansion

// --- CD4+ T cell dynamics ---
CD4T_basal     : 100    : Baseline antigen-specific CD4+ T cells (AU)
k_CD4_prol     : 0.07   : DGP-driven CD4+ T cell proliferation rate
k_CD4_death    : 0.03   : h-1, antigen-specific T cell death rate
EC50_CD4_DGP   : 0.5    : AU, EC50 of DGP for CD4+ T cell activation

// --- Th1/IFN-γ axis ---
k_IFNg_prod    : 0.18   : ng/mL/h, IFN-γ production rate per unit CD4+ T cell
k_IFNg_deg     : 0.35   : h-1, IFN-γ degradation (half-life ~2 h)

// --- Th17/IL-17A axis ---
k_IL17_prod    : 0.09   : ng/mL/h, IL-17A production (Th17)
k_IL17_deg     : 0.30   : h-1

// --- IL-21 (Th17/Tfh) ---
k_IL21_prod    : 0.10   : ng/mL/h, IL-21 production (Th17/Tfh)
k_IL21_deg     : 0.25   : h-1

// --- B cell & plasma cell dynamics ---
Bcell_basal    : 50     : Baseline B cells (AU, lamina propria)
k_Bprol        : 0.05   : DGP + IL-21 driven B cell proliferation
k_Bdeath       : 0.02   : h-1, B cell apoptosis

// --- Anti-tTG IgA serology ---
k_AntiTTG_prod : 0.06   : U/mL/h, anti-tTG IgA production per B cell
k_AntiTTG_deg  : 0.005  : h-1, IgA turnover (half-life ~5.8 days)
AntiTTG_ULN   : 10     : Upper limit of normal (U/mL); positive if >10 ULN

// --- Intestinal histopathology (Marsh model) ---
VH0            : 1.0    : Normalized baseline villous height (1=normal)
VH_min         : 0.05   : Minimum villous height (complete atrophy)
k_VH_damage    : 0.004  : h-1, IFN-γ/TNF-α-mediated villous damage rate
k_VH_repair    : 0.002  : h-1, epithelial repair rate on antigen withdrawal
EC50_VH_IFNg   : 4.0    : ng/mL, EC50 of IFN-γ for villous damage
CD0            : 1.0    : Normalized baseline crypt depth (1=normal)
k_CD_hyp       : 0.003  : h-1, IL-15-driven crypt hyperplasia rate
k_CD_norm      : 0.001  : h-1, crypt normalization rate

// --- Nutritional compartments ---
Iron0          : 1000   : Baseline iron stores (mg)
k_Iron_abs_hr  : 0.083  : mg/h, baseline iron absorption (= 2 mg/day)
k_Iron_loss_hr : 0.042  : mg/h, iron utilization/loss (= 1 mg/day)
BMD0           : 1.0    : Baseline BMD (normalized; 1=normal)
k_BMD_loss     : 0.0001 : Fractional BMD loss per unit absorption deficit (per h)
k_BMD_repair   : 0.00005 : Fractional BMD repair rate on adequate Ca absorption (per h)

// --- Drug PK parameters (Larazotide example) ---
F_oral         : 0.5    : Bioavailability (gut-local drug; systemic F~0.01)
ka_drug        : 0.6    : h-1, first-order absorption rate
CL_drug_L_h   : 8.0    : L/h, total drug clearance
Vd_drug_L     : 6.0    : L, volume of distribution
MW_drug        : 406    : g/mol, molecular weight (larazotide)

// --- Drug PD parameters ---
Drug_type      : 0      : 0=none, 1=Larazotide (IP), 2=ZED1227 (tTG2), 3=AMG714 (IL-15)
EC50_lara_IP   : 50     : ng/mL, EC50 of larazotide for IP reduction
EC50_ZED_tTG   : 80     : ng/mL, EC50 of ZED1227 for tTG2 inhibition
EC50_AMG_IL15  : 15     : ng/mL, EC50 of AMG714 for IL-15 neutralization
Emax_drug      : 0.85   : Maximum drug efficacy (85% of pathway inhibition)

$INIT @annotated
GlutenLumen  : 0    : Gluten in gut lumen (g)
GlutenLP     : 0    : Gliadin peptides in lamina propria (AU)
DGP          : 0    : Deamidated gliadin peptides in LP (AU)
tTG2_act     : 1.0  : tTG2 enzyme activity (normalized)
IP           : 1.0  : Intestinal permeability index (1=baseline)
IL15         : 0.1  : Epithelial IL-15 concentration (ng/mL)
IEL          : 20   : Intraepithelial lymphocyte count (per 100 enterocytes)
CD4T         : 100  : CD4+ antigen-specific T cells (AU)
IFNg         : 0.1  : IFN-γ concentration (ng/mL)
IL17         : 0.05 : IL-17A concentration (ng/mL)
IL21         : 0.05 : IL-21 concentration (ng/mL)
Bcell        : 50   : B cells in lamina propria (AU)
AntiTTG      : 2.0  : Anti-tTG IgA level (U/mL)
VH           : 1.0  : Villous height (normalized; 1=normal)
CrD          : 1.0  : Crypt depth (normalized; 1=normal)
AbsArea      : 1.0  : Intestinal absorption surface area (normalized)
IronStores   : 1000 : Body iron stores (mg)
BMD          : 1.0  : Bone mineral density (normalized)
DrugGut      : 0    : Drug in GI tract (mg)
DrugPlasma   : 0    : Drug in plasma (ng/mL)

$ODE
// ---- Effective gluten input (accounting for GFD) ----
double gluten_rate_h = GlutenIntake_g_day *
    ((1.0 - GFD) + GFD * GFD_leak) / 24.0;

// ---- Drug effect functions ----
double Cp = DrugPlasma;  // ng/mL
double E_lara  = (Drug_type == 1) ? Emax_drug * Cp / (Cp + EC50_lara_IP)  : 0.0;
double E_ZED   = (Drug_type == 2) ? Emax_drug * Cp / (Cp + EC50_ZED_tTG)  : 0.0;
double E_AMG   = (Drug_type == 3) ? Emax_drug * Cp / (Cp + EC50_AMG_IL15) : 0.0;

// ---- Hill activation helpers ----
double h_IL15_IP  = IL15 / (IL15 + EC50_IP_IL15);
double h_IL15_IEL = IL15 / (IL15 + EC50_IEL_IL15);
double h_glia_IL15 = GlutenLP / (GlutenLP + EC50_IL15_glia);
double h_DGP_CD4   = DGP / (DGP + EC50_CD4_DGP);
double h_IFNg_VH   = IFNg / (IFNg + EC50_VH_IFNg);

// -------- ODEs --------

// [1] Gluten in lumen (g): input - degradation - absorption × permeability
dxdt_GlutenLumen = gluten_rate_h
    - kdeg_lumen * GlutenLumen
    - kabs_gliadin * GlutenLumen * IP;

// [2] Gliadin peptides in lamina propria (AU)
dxdt_GlutenLP = kabs_gliadin * GlutenLumen * IP
    - k_deamid * tTG2_act * GlutenLP * (1.0 - E_ZED)  // tTG2 deamidates → DGP
    - kdeg_DGP * GlutenLP;                               // non-deamidated clearance

// [3] Deamidated gliadin peptides (DGP)
dxdt_DGP = k_deamid * tTG2_act * GlutenLP * (1.0 - E_ZED)
    - kdeg_DGP * DGP;

// [4] tTG2 activity: upregulated by TGF-β / IFN-γ; inhibited by ZED1227
dxdt_tTG2_act = k_tTG2_TGFb * (1.0 + 0.3 * IFNg / (IFNg + 2.0))
    - k_tTG2_deg * tTG2_act;

// [5] Intestinal permeability: IL-15 zonulin pathway; restored by Larazotide
dxdt_IP = k_IP_IL15 * h_IL15_IP * (1.0 - E_lara)
    - k_IP_restore * (IP - IP_basal);

// [6] Epithelial IL-15: induced by gliadin; suppressed by AMG714
dxdt_IL15 = k_IL15_prod * h_glia_IL15 * (1.0 - E_AMG)
    - k_IL15_deg * IL15;

// [7] IEL: expanded by IL-15 signaling; max capacity constraint
double IEL_prol = k_IEL_IL15 * h_IL15_IEL * IEL * (1.0 - IEL / IEL_max);
dxdt_IEL = IEL_prol - k_IEL_death * (IEL - IEL_normal);

// [8] Antigen-specific CD4+ T cells: expanded by DGP-MHC signaling
dxdt_CD4T = k_CD4_prol * h_DGP_CD4 * CD4T
    - k_CD4_death * (CD4T - CD4T_basal);

// [9] IFN-γ (Th1 signature cytokine)
dxdt_IFNg = k_IFNg_prod * (CD4T / 100.0) * 0.4
    - k_IFNg_deg * IFNg;

// [10] IL-17A (Th17 cytokine)
dxdt_IL17 = k_IL17_prod * (CD4T / 100.0) * 0.2 * h_IL15_IP
    - k_IL17_deg * IL17;

// [11] IL-21 (Th17/Tfh; B cell helper and IL-15 amplifier)
dxdt_IL21 = k_IL21_prod * (CD4T / 100.0) * 0.3
    - k_IL21_deg * IL21;

// [12] B cells: DGP antigen + IL-21 T cell help
double Bcell_prol = k_Bprol * h_DGP_CD4 * (IL21 / (IL21 + 0.1)) * Bcell;
dxdt_Bcell = Bcell_prol - k_Bdeath * Bcell;

// [13] Anti-tTG IgA serology (U/mL; diagnostic marker)
dxdt_AntiTTG = k_AntiTTG_prod * (Bcell / 50.0)
    - k_AntiTTG_deg * AntiTTG;

// [14] Villous height (normalized; damage by IFN-γ/IL-17; repair on withdrawal)
double VH_damage = k_VH_damage * h_IFNg_VH * (IL17 / (IL17 + 1.0)) * VH;
double VH_repair = k_VH_repair * (1.0 - h_DGP_CD4) * (VH0 - VH);
dxdt_VH = fmax(VH_repair - VH_damage, -VH + VH_min);

// [15] Crypt depth (normalized; hyperplasia by IL-15; normalization on healing)
dxdt_CrD = k_CD_hyp * h_IL15_IP * (2.5 - CrD)   // expand toward 2.5× max
    - k_CD_norm * (CrD - CD0) * (GlutenLP < 0.1 ? 1.0 : 0.0);

// [16] Absorption surface area (tracks villous height with delay)
double Abs_target = fmax(0.05, VH / (VH + 0.2));  // sigmoid of VH
dxdt_AbsArea = 0.015 * (Abs_target - AbsArea);

// [17] Iron stores (mg): abs proportional to surface area
double Iron_abs = k_Iron_abs_hr * AbsArea;
dxdt_IronStores = Iron_abs - k_Iron_loss_hr;

// [18] BMD (normalized)
dxdt_BMD = k_BMD_repair * AbsArea
    - k_BMD_loss * (1.0 - AbsArea);

// [19] Drug in GI tract (mg): cleared by absorption
dxdt_DrugGut = -ka_drug * DrugGut;

// [20] Drug in plasma (ng/mL)
dxdt_DrugPlasma = F_oral * ka_drug * DrugGut * 1000.0 / Vd_drug_L
    - (CL_drug_L_h / Vd_drug_L) * DrugPlasma;

$TABLE
// Derived clinical outputs
capture VH_CD_ratio  = VH / CrD;            // V:C ratio (normal ~3:1; atrophy <1)
capture Marsh_score  = (VH_CD_ratio < 0.3) ? 3 :
                       (VH_CD_ratio < 0.7) ? 2 :
                       (VH_CD_ratio < 1.0) ? 1 : 0;
capture Serology_pos = (AntiTTG > AntiTTG_ULN) ? 1 : 0;  // positive serology
capture Hgb_g_dL     = 8.0 + 6.0 * (IronStores / 1000.0); // simplified Hgb proxy
capture Ferritin_ug  = 15.0 * (IronStores / 1000.0);       // ferritin proxy
capture BMD_Tscore   = (BMD - 1.0) / 0.1 * (-1);           // T-score proxy
capture IEL_elevated = (IEL > 25) ? 1 : 0;                  // Marsh 1+ flag
capture GFD_flag     = GFD;
capture Drug_Cp      = DrugPlasma;

$CAPTURE VH_CD_ratio Marsh_score Serology_pos Hgb_g_dL Ferritin_ug
         BMD_Tscore IEL_elevated GFD_flag Drug_Cp AntiTTG VH CrD AbsArea
'

# ==============================================================================
# Build the model
# ==============================================================================
cd_mod <- mcode("CeliacDisease_QSP", cd_model_code)

# ==============================================================================
# Treatment Scenarios
# ==============================================================================
scenarios <- list(
  "1_Untreated_Normal_Diet" = list(
    GFD = 0, GFD_leak = 0, Drug_type = 0, Drug_dose = 0,
    label = "Untreated — Normal Diet (10-12 g gluten/day)",
    color = "#E53935"
  ),
  "2_Strict_GFD" = list(
    GFD = 1, GFD_leak = 0, Drug_type = 0, Drug_dose = 0,
    label = "Strict GFD (0 g gluten/day)",
    color = "#43A047"
  ),
  "3_Partial_GFD_5pct_leak" = list(
    GFD = 1, GFD_leak = 0.05, Drug_type = 0, Drug_dose = 0,
    label = "Partial GFD (5% trace contamination)",
    color = "#FB8C00"
  ),
  "4_GFD_plus_Larazotide" = list(
    GFD = 1, GFD_leak = 0.10, Drug_type = 1, Drug_dose = 2,
    label = "GFD + Larazotide (2 mg TID, 10% leakage)",
    color = "#1E88E5"
  ),
  "5_GFD_plus_ZED1227" = list(
    GFD = 1, GFD_leak = 0.10, Drug_type = 2, Drug_dose = 300,
    label = "GFD + ZED1227 TG2i (300 mg QD, 10% leakage)",
    color = "#8E24AA"
  ),
  "6_GFD_plus_AMG714_RCD" = list(
    GFD = 1, GFD_leak = 0.20, Drug_type = 3, Drug_dose = 150,
    label = "GFD + AMG714 anti-IL-15 (150 mg SC, RCD)",
    color = "#00897B"
  )
)

# Simulation duration: 2 years (17520 h) with burn-in
sim_time <- seq(0, 8760, by = 24)  # 1 year, daily sampling

run_scenario <- function(scen_name, scen_params) {
  # Create dosing event (for drug scenarios)
  if (scen_params$Drug_dose > 0) {
    dose_interval_h <- switch(as.character(scen_params$Drug_type),
      "1" = 8,    # Larazotide TID
      "2" = 24,   # ZED1227 QD
      "3" = 168,  # AMG714 weekly SC
      24
    )
    ev <- ev(
      time = seq(0, max(sim_time) - dose_interval_h, by = dose_interval_h),
      amt  = scen_params$Drug_dose,
      cmt  = "DrugGut",
      rate = -2  # bolus
    )
  } else {
    ev <- ev(time = 0, amt = 0, cmt = "DrugGut")
  }

  # Set parameters
  params <- param(cd_mod,
    GFD       = scen_params$GFD,
    GFD_leak  = scen_params$GFD_leak,
    Drug_type = scen_params$Drug_type
  )

  # Simulate
  out <- mrgsim(
    x    = params,
    events = ev,
    end  = max(sim_time),
    delta = 24,
    obsonly = TRUE
  ) %>%
    as.data.frame() %>%
    mutate(
      Scenario  = scen_name,
      Label     = scen_params$label,
      Time_days = time / 24
    )
  out
}

# Run all scenarios
results <- lapply(names(scenarios), function(nm) {
  run_scenario(nm, scenarios[[nm]])
}) %>% bind_rows()

# Add color map
color_map <- sapply(scenarios, `[[`, "color")
names(color_map) <- sapply(scenarios, `[[`, "label")

# ==============================================================================
# PLOT 1: Disease Biomarker Dynamics
# ==============================================================================
p1 <- ggplot(results, aes(Time_days, AntiTTG, color = Label)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "gray40") +
  annotate("text", x = 300, y = 11.5, label = "Positive threshold (ULN=10 U/mL)",
           color = "gray40", size = 3.5) +
  scale_color_manual(values = color_map) +
  labs(title = "Anti-tTG IgA Serology (U/mL)", x = "Time (days)", y = "Anti-tTG IgA (U/mL)",
       color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

p2 <- ggplot(results, aes(Time_days, VH, color = Label)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = color_map) +
  labs(title = "Villous Height (normalized)", x = "Time (days)", y = "VH (1=normal)",
       color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p3 <- ggplot(results, aes(Time_days, VH_CD_ratio, color = Label)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  annotate("text", x = 50, y = 3.15, label = "Normal V:C ratio ~3:1",
           color = "steelblue", size = 3.5) +
  annotate("text", x = 50, y = 0.85, label = "Atrophy threshold V:C <1",
           color = "red", size = 3.5) +
  scale_color_manual(values = color_map) +
  labs(title = "Villi:Crypt (V:C) Ratio", x = "Time (days)", y = "V:C Ratio",
       color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p4 <- ggplot(results, aes(Time_days, IFNg, color = Label)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = color_map) +
  labs(title = "IFN-γ (Th1 Cytokine, ng/mL)", x = "Time (days)", y = "IFN-γ (ng/mL)",
       color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p5 <- ggplot(results, aes(Time_days, Hgb_g_dL, color = Label)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 12, linetype = "dashed", color = "tomato") +
  scale_color_manual(values = color_map) +
  labs(title = "Hemoglobin (proxy, g/dL)", x = "Time (days)", y = "Hgb (g/dL)",
       color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

p6 <- ggplot(results, aes(Time_days, BMD, color = Label)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = color_map) +
  labs(title = "Bone Mineral Density (normalized)", x = "Time (days)",
       y = "BMD (1=normal)", color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# Combine
combined_plot <- (p1 / (p2 | p3) / (p4 | p5 | p6)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ==============================================================================
# PLOT 2: Marsh Score Trajectory
# ==============================================================================
marsh_labels <- c("0" = "Marsh 0 (Normal)", "1" = "Marsh 1 (↑IEL)",
                  "2" = "Marsh 2 (+Crypt)", "3" = "Marsh 3 (Atrophy)")

p_marsh <- ggplot(results, aes(Time_days, Marsh_score, color = Label)) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(breaks = 0:3, labels = marsh_labels) +
  scale_color_manual(values = color_map) +
  labs(title = "Marsh Histopathology Score Over Time",
       x = "Time (days)", y = "Marsh Score", color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 9))

# ==============================================================================
# PLOT 3: Cytokine Panel
# ==============================================================================
cyto_long <- results %>%
  select(Time_days, Label, IFNg, IL17, IL21, IL15, IEL) %>%
  pivot_longer(cols = c(IFNg, IL17, IL21, IL15, IEL),
               names_to = "Analyte", values_to = "Value")

p_cyto <- ggplot(cyto_long, aes(Time_days, Value, color = Label)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Analyte, scales = "free_y", ncol = 3) +
  scale_color_manual(values = color_map) +
  labs(title = "Cytokine & Immune Cell Dynamics",
       x = "Time (days)", y = "Concentration (ng/mL) / Count", color = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

# ==============================================================================
# PLOT 4: Endpoint Summary at 1 Year
# ==============================================================================
summary_1yr <- results %>%
  filter(Time_days == max(Time_days)) %>%
  group_by(Label) %>%
  summarize(
    Anti_tTG_IgA  = round(mean(AntiTTG), 1),
    VH_ratio      = round(mean(VH), 3),
    VC_ratio      = round(mean(VH_CD_ratio), 2),
    Marsh_score   = round(mean(Marsh_score), 0),
    IFNg_ng_mL    = round(mean(IFNg), 2),
    Hgb_g_dL     = round(mean(Hgb_g_dL), 1),
    BMD_norm      = round(mean(BMD), 3),
    IEL_count     = round(mean(IEL), 1),
    Serology_pos  = max(Serology_pos),
    .groups = "drop"
  )

print("=== 1-Year Clinical Endpoint Summary ===")
print(summary_1yr)

# ==============================================================================
# Output
# ==============================================================================
message("\nCeliac Disease QSP Model — Simulation Complete")
message("Scenarios simulated: ", length(unique(results$Scenario)))
message("Time points per scenario: ", length(unique(results$Time_days)))
message("Key model compartments: 21 ODEs")
message("Key outputs: Anti-tTG IgA, Marsh score, V:C ratio, IFN-γ, Hgb, BMD")

list(
  model    = cd_mod,
  results  = results,
  summary  = summary_1yr,
  plots    = list(
    biomarkers = combined_plot,
    marsh      = p_marsh,
    cytokines  = p_cyto
  )
)
