# =============================================================================
# Bladder Cancer (Urothelial Carcinoma) — QSP mrgsolve Model
# =============================================================================
# 20 ODE compartments, 7 treatment scenarios
# Calibrated to: SWOG S8507 (BCG), von der Maase 2000 (GC),
#   KEYNOTE-045 (pembrolizumab), IMvigor210/211 (atezolizumab),
#   BLC2001 (erdafitinib), EV-301 (enfortumab vedotin)
# =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

blca_model_code <- '
$PROB Bladder Cancer (BLCA) QSP Model — 20 compartments, 7 scenarios

$PARAM
// ---- Patient flags ----
FGFR3_flag   = 0      // 1 = FGFR3-altered patient (erdafitinib eligible)
Nectin4_flag = 1      // 1 = Nectin-4+ (EV eligible; near-universal)
TROP2_flag   = 1      // 1 = TROP2+ (SG eligible; near-universal)
PDL1_CPS     = 0      // 1 = PDL1 CPS>=10 (preferred pembrolizumab)

// ---- PK parameters: BCG intravesical ----
kabs_BCG  = 0.5       // h-1, absorption into bladder wall effect compartment
kelim_BCG = 0.12      // h-1, BCG elimination from wall

// ---- PK parameters: Cisplatin (70 mg/m²) ----
kelim_Cis = 0.35      // h-1, renal + nonrenal clearance (t1/2 ~2h)

// ---- PK parameters: Gemcitabine (1000 mg/m²) ----
kelim_Gem = 1.2       // h-1, rapid metabolic clearance (t1/2 ~0.6h)

// ---- PK parameters: Pembrolizumab 200 mg q3w (2-cmpt) ----
CL_Pembro  = 0.22     // L/h
V1_Pembro  = 3.4      // L (central)
Q_Pembro   = 0.55     // L/h intercompartmental
V2_Pembro  = 4.8      // L (peripheral)
IC50_Pembro = 0.40    // ug/mL, PD-1 half-maximal binding
hill_Pembro = 1.5

// ---- PK parameters: Atezolizumab 1200 mg q3w (2-cmpt) ----
CL_Atezo   = 0.20     // L/h
V1_Atezo   = 3.2      // L
Q_Atezo    = 0.48     // L/h
V2_Atezo   = 5.1      // L
IC50_Atezo = 0.60     // ug/mL, PD-L1 binding
hill_Atezo = 1.4

// ---- PK parameters: Erdafitinib 8 mg PO QD ----
ka_Erda    = 0.8      // h-1, oral absorption
F_Erda     = 0.70     // bioavailability 70%
kelim_Erda = 0.058    // h-1, t1/2 ~12h, CYP2C9/3A4
EC50_Erda  = 0.12     // ug/mL, FGFR1-4 inhibition
hill_Erda  = 1.8

// ---- PK parameters: Enfortumab Vedotin 1.25 mg/kg (q28d d1,8,15) ----
kelim_EV   = 0.045    // h-1, ADC elimination (~t1/2 3-4d)
EC50_EV    = 0.08     // ug/mL, MMAE-mediated kill
hill_EV    = 1.6

// ---- PD: BCG immune effects ----
Emax_BCG_CD8  = 3.0   // max fold stimulation of CD8 by BCG
EC50_BCG_CD8  = 0.5   // BCG_eff at half-max (normalized)
Emax_BCG_kill = 0.65  // max direct tumor kill fraction
EC50_BCG_kill = 0.6
hill_BCG      = 1.5

// ---- PD: Cisplatin tumor kill ----
Emax_Cis  = 0.55      // max kill rate coefficient h-1
EC50_Cis  = 8.0       // ug/mL
hill_Cis  = 1.2

// ---- PD: Gemcitabine tumor kill ----
Emax_Gem  = 0.45
EC50_Gem  = 2.5       // ug/mL
hill_Gem  = 1.3

// ---- PD: IO-enhanced CTL kill ----
Emax_IO_kill  = 0.30  // max additional kill rate via IO-rescued CD8
CD8_EC50_kill = 0.8   // CD8 concentration for half-max CTL kill

// ---- PD: Erdafitinib FGFR kill ----
Emax_Erda_kill = 0.40 // (FGFR3-altered, ORR 40.4% BLC2001)

// ---- PD: Enfortumab Vedotin ADC kill ----
Emax_EV_kill  = 0.42  // (ORR 40.6% EV-301)

// ---- Tumor biology ----
kg_tumor   = 0.0025   // h-1, net tumor growth rate (doubling ~12 days NMIBC)
Kmax_tumor = 1.0e10   // cells, carrying capacity
T0_cells   = 1.0e8    // baseline tumor burden (NMIBC)
SLD0       = 55.0     // mm, baseline sum of longest diameters (MIBC/mUBC)

// ---- Immune compartment steady-state ----
kin_CD8    = 0.02     // cells/h reference production
kout_CD8   = 0.008    // h-1 CD8 natural decay
kin_Treg   = 0.005
kout_Treg  = 0.003
kin_MDSC   = 0.004
kout_MDSC  = 0.005
suppTreg   = 0.30     // max suppression coefficient (Treg on CD8)
suppMDSC   = 0.25     // max suppression coefficient (MDSC on CD8)
Treg_EC50  = 1.0      // Treg level at half-max suppression
MDSC_EC50  = 1.0

// ---- IFN-gamma ----
kprod_IFNg = 0.015    // h-1 produced per CD8 unit
kdeg_IFNg  = 0.10     // h-1

// ---- PD-L1 expression ----
kprod_PDL1 = 0.008
kdeg_PDL1  = 0.04
IFNg_stim  = 2.5      // fold IFN-gamma-mediated PD-L1 upregulation

// ---- FGFR3 pathway activity ----
kprod_FGFR3 = 0.05
kdeg_FGFR3  = 0.02

// ---- NMP22 biomarker ----
kprod_NMP22 = 0.003
kdeg_NMP22  = 0.015

// ---- SLD (Sum Longest Diameters) dynamics ----
kg_SLD  = 0.00015     // h-1 SLD growth proportional to tumor burden

$CMT
BCG_depot BCG_eff
Cis_plasm Gem_plasm
Pembro_c Pembro_p
Atezo_c Atezo_p
Erda_dep Erda_plasm
EnFV_c
CD8_eff Treg MDSC_cmt
TumBurd
FGFR3act PDL1_lvl IFNg_cmt
NMP22_cmt SLD_cmt

$MAIN
// Initialize steady-state immune compartments
if (NEWIND <= 1) {
  double CD8_ss  = kin_CD8 / kout_CD8;
  double Treg_ss = kin_Treg / kout_Treg;
  double MDSC_ss = kin_MDSC / kout_MDSC;

  CD8_eff_0   = CD8_ss;
  Treg_0      = Treg_ss;
  MDSC_cmt_0  = MDSC_ss;
  TumBurd_0   = T0_cells;
  FGFR3act_0  = (FGFR3_flag > 0.5) ? kprod_FGFR3 / kdeg_FGFR3 : 0.1;
  PDL1_lvl_0  = kprod_PDL1 / kdeg_PDL1;
  IFNg_cmt_0  = kprod_IFNg * CD8_ss / kdeg_IFNg;
  NMP22_cmt_0 = kprod_NMP22 * T0_cells / kdeg_NMP22;
  SLD_cmt_0   = SLD0;
}

$ODE
// ---------- BCG PK ----------
dxdt_BCG_depot = -kabs_BCG * BCG_depot;
dxdt_BCG_eff   =  kabs_BCG * BCG_depot - kelim_BCG * BCG_eff;

// ---------- Chemotherapy PK ----------
dxdt_Cis_plasm = -kelim_Cis * Cis_plasm;
dxdt_Gem_plasm = -kelim_Gem * Gem_plasm;

// ---------- Pembrolizumab 2-cmpt ----------
double k10_P  = CL_Pembro / V1_Pembro;
double k12_P  = Q_Pembro  / V1_Pembro;
double k21_P  = Q_Pembro  / V2_Pembro;
dxdt_Pembro_c = -(k10_P + k12_P) * Pembro_c + k21_P * Pembro_p;
dxdt_Pembro_p =  k12_P * Pembro_c - k21_P * Pembro_p;

// ---------- Atezolizumab 2-cmpt ----------
double k10_A  = CL_Atezo / V1_Atezo;
double k12_A  = Q_Atezo  / V1_Atezo;
double k21_A  = Q_Atezo  / V2_Atezo;
dxdt_Atezo_c  = -(k10_A + k12_A) * Atezo_c + k21_A * Atezo_p;
dxdt_Atezo_p  =  k12_A * Atezo_c - k21_A * Atezo_p;

// ---------- Erdafitinib oral PK ----------
dxdt_Erda_dep   = -ka_Erda * Erda_dep;
dxdt_Erda_plasm =  F_Erda * ka_Erda * Erda_dep - kelim_Erda * Erda_plasm;

// ---------- Enfortumab Vedotin ----------
dxdt_EnFV_c = -kelim_EV * EnFV_c;

// ---------- Drug effects (Hill equations) ----------
double BCG_E_norm = BCG_eff / (1.0 + BCG_eff); // simplified 0-1

double E_BCG_CD8  = Emax_BCG_CD8 * pow(BCG_eff, hill_BCG) /
                    (pow(EC50_BCG_CD8, hill_BCG) + pow(BCG_eff, hill_BCG));
double E_BCG_kill = Emax_BCG_kill * pow(BCG_eff, hill_BCG) /
                    (pow(EC50_BCG_kill, hill_BCG) + pow(BCG_eff, hill_BCG));

double E_Cis  = Emax_Cis * pow(Cis_plasm, hill_Cis) /
                (pow(EC50_Cis, hill_Cis) + pow(Cis_plasm, hill_Cis));
double E_Gem  = Emax_Gem * pow(Gem_plasm, hill_Gem) /
                (pow(EC50_Gem, hill_Gem) + pow(Gem_plasm, hill_Gem));

double Pembro_conc_mgL = Pembro_c / V1_Pembro;
double Atezo_conc_mgL  = Atezo_c  / V1_Atezo;

double Pembro_RO = pow(Pembro_conc_mgL, hill_Pembro) /
                   (pow(IC50_Pembro, hill_Pembro) + pow(Pembro_conc_mgL, hill_Pembro));
double Atezo_RO  = pow(Atezo_conc_mgL, hill_Atezo) /
                   (pow(IC50_Atezo, hill_Atezo) + pow(Atezo_conc_mgL, hill_Atezo));
double IO_RO_total = (Pembro_RO > Atezo_RO) ? Pembro_RO : Atezo_RO;

double E_Erda = FGFR3_flag * Emax_Erda_kill * pow(Erda_plasm, hill_Erda) /
                (pow(EC50_Erda, hill_Erda) + pow(Erda_plasm, hill_Erda));

double E_EV   = Nectin4_flag * Emax_EV_kill * pow(EnFV_c, hill_EV) /
                (pow(EC50_EV, hill_EV) + pow(EnFV_c, hill_EV));

// ---------- CD8 dynamics ----------
double suppression = 1.0 + suppTreg * Treg / (Treg_EC50 + Treg) +
                           suppMDSC * MDSC_cmt / (MDSC_EC50 + MDSC_cmt);
double CD8_input   = kin_CD8 * (1.0 + E_BCG_CD8 + 3.0 * IO_RO_total);
dxdt_CD8_eff = CD8_input - kout_CD8 * suppression * CD8_eff;

// ---------- Treg dynamics ----------
double TGFb_drive = TumBurd / (1.0e9 + TumBurd);
dxdt_Treg = kin_Treg * (1.0 + 2.0 * TGFb_drive) - kout_Treg * Treg;

// ---------- MDSC dynamics ----------
double tumor_drive = TumBurd / (5.0e9 + TumBurd);
dxdt_MDSC_cmt = kin_MDSC * (1.0 + 1.5 * tumor_drive) - kout_MDSC * MDSC_cmt;

// ---------- Tumor burden (logistic + kills) ----------
double CTL_kill = Emax_IO_kill * CD8_eff * IO_RO_total /
                  (CD8_EC50_kill + CD8_eff * IO_RO_total);
double BCG_CTL  = E_BCG_kill * CD8_eff / (0.5 + CD8_eff);

double kdeath_total = BCG_CTL + (E_Cis + E_Gem) * 0.012 +
                      CTL_kill * 0.01 +
                      E_Erda * 0.015 + E_EV * 0.014;
double Tgrow = kg_tumor * TumBurd * (1.0 - TumBurd / Kmax_tumor);
dxdt_TumBurd = Tgrow - kdeath_total * TumBurd;
if (TumBurd < 1.0) dxdt_TumBurd = 0.0;

// ---------- FGFR3 pathway activity ----------
double kinh_FGFR3 = E_Erda / (0.01 + E_Erda);
dxdt_FGFR3act = kprod_FGFR3 * FGFR3_flag - (kdeg_FGFR3 + 0.05 * kinh_FGFR3) * FGFR3act;

// ---------- PD-L1 expression ----------
double IFNg_norm = IFNg_cmt / (0.1 + IFNg_cmt);
dxdt_PDL1_lvl = kprod_PDL1 * (1.0 + IFNg_stim * IFNg_norm) - kdeg_PDL1 * PDL1_lvl;

// ---------- IFN-gamma ----------
double BCG_IFNg = 1.5 * E_BCG_CD8;
dxdt_IFNg_cmt  = kprod_IFNg * CD8_eff * (1.0 + BCG_IFNg) - kdeg_IFNg * IFNg_cmt;

// ---------- NMP22 biomarker ----------
dxdt_NMP22_cmt = kprod_NMP22 * TumBurd - kdeg_NMP22 * NMP22_cmt;

// ---------- SLD (Sum Longest Diameters) ----------
double T_norm  = TumBurd / T0_cells;
dxdt_SLD_cmt   = kg_SLD * SLD_cmt * (T_norm - 1.0);

$TABLE
// Capture variables
double Pembro_conc_ugmL = Pembro_c / V1_Pembro;
double Atezo_conc_ugmL  = Atezo_c  / V1_Atezo;
double PD1_RO  = pow(Pembro_conc_ugmL, hill_Pembro) /
                 (pow(IC50_Pembro, hill_Pembro) + pow(Pembro_conc_ugmL, hill_Pembro));
double PDL1_RO = pow(Atezo_conc_ugmL, hill_Atezo) /
                 (pow(IC50_Atezo, hill_Atezo) + pow(Atezo_conc_ugmL, hill_Atezo));

double TumBurd_log = log10(TumBurd + 1.0);
double TumorRed_pct = 100.0 * (T0_cells - TumBurd) / T0_cells;
if (TumorRed_pct < -100.0) TumorRed_pct = -100.0;

double SLD_change_pct = 100.0 * (SLD_cmt - SLD0) / SLD0;

// Simplified ARR (annualized recurrence rate for NMIBC, proxy)
double ARR_proxy = 1.0 - exp(-TumBurd / T0_cells * 0.8);

double E_BCG_kill_out  = Emax_BCG_kill * pow(BCG_eff, hill_BCG) /
                         (pow(EC50_BCG_kill, hill_BCG) + pow(BCG_eff, hill_BCG));
double E_Cis_out = Emax_Cis * pow(Cis_plasm, hill_Cis) /
                   (pow(EC50_Cis, hill_Cis) + pow(Cis_plasm, hill_Cis));
double E_Gem_out = Emax_Gem * pow(Gem_plasm, hill_Gem) /
                   (pow(EC50_Gem, hill_Gem) + pow(Gem_plasm, hill_Gem));

$CAPTURE
BCG_eff Cis_plasm Gem_plasm
Pembro_conc_ugmL Atezo_conc_ugmL PD1_RO PDL1_RO
Erda_plasm EnFV_c
CD8_eff Treg MDSC_cmt
TumBurd TumBurd_log TumorRed_pct SLD_cmt SLD_change_pct ARR_proxy
FGFR3act PDL1_lvl IFNg_cmt NMP22_cmt
E_BCG_kill_out E_Cis_out E_Gem_out
'

# =============================================================================
# Model compilation
# =============================================================================
mod <- mcode("blca_qsp", blca_model_code)

# =============================================================================
# Dosing event builders
# =============================================================================

# BCG intravesical: 81 mg weekly x6 then q3w maintenance (SWOG S8507)
# Dose in model units (mg) — BCG_depot receives the dose
make_BCG_events <- function(induction_weeks = 6, maint_cycles = 3, dur_days = 365) {
  induction <- ev(cmt = "BCG_depot", amt = 81,
                  time = seq(0, by = 7*24, length.out = induction_weeks))
  maint_start <- induction_weeks * 7 * 24
  maint_times <- maint_start + seq(0, by = 3*7*24, length.out = maint_cycles)
  maint <- ev(cmt = "BCG_depot", amt = 81, time = maint_times)
  c(induction, maint)
}

# GC: Cisplatin 70 mg/m² + Gemcitabine 1000 mg/m² q3w x6 cycles
make_GC_events <- function(cycles = 6) {
  times_h <- seq(0, by = 21*24, length.out = cycles)
  cis <- ev(cmt = "Cis_plasm", amt = 70 * 1.73, time = times_h)
  gem <- ev(cmt = "Gem_plasm", amt = 1000 * 1.73, time = times_h)
  c(cis, gem)
}

# Pembrolizumab 200 mg IV q3w
make_Pembro_events <- function(cycles = 18, dose_mg = 200) {
  times_h <- seq(0, by = 21*24, length.out = cycles)
  ev(cmt = "Pembro_c", amt = dose_mg, time = times_h)
}

# Atezolizumab 1200 mg IV q3w
make_Atezo_events <- function(cycles = 18, dose_mg = 1200) {
  times_h <- seq(0, by = 21*24, length.out = cycles)
  ev(cmt = "Atezo_c", amt = dose_mg, time = times_h)
}

# Erdafitinib 8 mg PO QD (continuous)
make_Erda_events <- function(days = 365, dose_mg = 8) {
  times_h <- seq(0, by = 24, length.out = days)
  ev(cmt = "Erda_dep", amt = dose_mg, time = times_h)
}

# Enfortumab vedotin 1.25 mg/kg IV days 1, 8, 15 q28d
make_EV_events <- function(cycles = 8, dose_mgkg = 1.25, weight_kg = 70) {
  dose_mg <- dose_mgkg * weight_kg
  times_h <- unlist(lapply(0:(cycles - 1), function(c) {
    c * 28 * 24 + c(0, 7*24, 14*24)
  }))
  ev(cmt = "EnFV_c", amt = dose_mg, time = times_h)
}

# =============================================================================
# Run all 7 scenarios
# =============================================================================
run_scenario <- function(drug, FGFR3 = 0, Nectin4 = 1, duration_d = 365) {
  end_h <- duration_d * 24
  delta <- 6  # hours step

  base_par <- list(FGFR3_flag = FGFR3, Nectin4_flag = Nectin4)

  if (drug == "untreated") {
    dose <- ev(time = 0, amt = 0, cmt = 1)
  } else if (drug == "BCG") {
    dose <- make_BCG_events()
  } else if (drug == "GC") {
    dose <- make_GC_events()
  } else if (drug == "pembrolizumab") {
    dose <- make_Pembro_events()
  } else if (drug == "atezolizumab") {
    dose <- make_Atezo_events()
  } else if (drug == "erdafitinib") {
    base_par$FGFR3_flag <- 1
    dose <- make_Erda_events(days = duration_d)
  } else if (drug == "enfortumab-vedotin") {
    dose <- make_EV_events()
  } else {
    stop("Unknown drug: ", drug)
  }

  mrgsim(mod, ev = dose, param = base_par,
         end = end_h, delta = delta, obsonly = TRUE) |>
    as.data.frame() |>
    mutate(Drug = drug, time_d = time / 24)
}

drugs_list <- c("untreated", "BCG", "GC",
                "pembrolizumab", "atezolizumab",
                "erdafitinib", "enfortumab-vedotin")

set.seed(42)
results <- lapply(drugs_list, function(d) run_scenario(d))
all_results <- bind_rows(results)

# =============================================================================
# Calibration reference table
# =============================================================================
calibration_table <- data.frame(
  Scenario = c("BCG (SWOG S8507)",
               "GC (von der Maase 2000)",
               "Pembrolizumab (KEYNOTE-045)",
               "Atezolizumab (IMvigor210)",
               "Erdafitinib (BLC2001)",
               "Enfortumab Vedotin (EV-301)"),
  Endpoint = c("CRR (CIS) 55-70%",
               "ORR 49%, mOS 13.8 mo",
               "ORR 21.1%, OS HR 0.73",
               "ORR 15% (all), 26% (IC2/3)",
               "ORR 40.4% (FGFR-altered)",
               "ORR 40.6%, OS HR 0.70"),
  Source = c("Lamm et al. J Urol 2000",
             "von der Maase et al. JCO 2000",
             "Bellmunt et al. NEJM 2017",
             "Balar et al. Lancet 2017",
             "Loriot et al. NEJM 2019",
             "Powles et al. NEJM 2021")
)

# =============================================================================
# Quick plots
# =============================================================================
p1 <- ggplot(all_results, aes(x = time_d, y = TumorRed_pct, color = Drug)) +
  geom_line(size = 1) +
  geom_hline(yintercept = c(-30, 30), linetype = "dashed", alpha = 0.4) +
  labs(title = "Bladder Cancer QSP — Tumor Reduction by Treatment",
       x = "Time (days)", y = "Tumor Reduction (%)",
       caption = "BLCA QSP Model | SWOG/KEYNOTE/IMvigor/BLC2001/EV-301 calibrated") +
  theme_bw()

p2 <- ggplot(all_results, aes(x = time_d, y = CD8_eff, color = Drug)) +
  geom_line(size = 1) +
  labs(title = "CD8+ Effector T Cell Dynamics",
       x = "Time (days)", y = "CD8+ CTL (relative units)") +
  theme_bw()

p3 <- ggplot(all_results, aes(x = time_d, y = SLD_change_pct, color = Drug)) +
  geom_line(size = 1) +
  geom_hline(yintercept = c(-30, 20), linetype = "dashed", alpha = 0.5) +
  labs(title = "Sum of Longest Diameters (SLD) Change",
       x = "Time (days)", y = "SLD Change from Baseline (%)") +
  theme_bw()

p4 <- ggplot(all_results, aes(x = time_d, y = NMP22_cmt, color = Drug)) +
  geom_line(size = 1) +
  labs(title = "NMP22 Urinary Biomarker Dynamics",
       x = "Time (days)", y = "NMP22 (relative units)") +
  theme_bw()

print(p1); print(p2); print(p3); print(p4)

cat("\n=== BLCA QSP Model Summary ===\n")
cat("Compartments: 20 ODEs (BCG PK, Cis/Gem PK, Pembro 2-cmpt,\n")
cat("  Atezo 2-cmpt, Erdafitinib 1-cmpt oral, EV 1-cmpt,\n")
cat("  CD8/Treg/MDSC immune, TumBurd, FGFR3act, PDL1, IFNg, NMP22, SLD)\n")
cat("Treatment scenarios (7):\n")
for (d in drugs_list) cat(" ", d, "\n")
cat("\nCalibration references:\n")
print(calibration_table, row.names = FALSE)
