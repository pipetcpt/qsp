###############################################################################
# Chronic Pyelonephritis (CPN) — QSP mrgsolve Model
# 만성 신우신염 정량적 시스템 약리학 모델
#
# Compartments (18 ODEs):
#   PK  : Cipro_gut, Cipro_C, Cipro_P, TMP_gut, TMP_C, NIT_gut, NIT_urine
#   Disease : Bacteria, Biofilm, Neutrophil, Macrophage, IL6, TGFb1,
#             Collagen, RenalScar, GFR
#   (+ Creatinine from TABLE)
#
# Key references:
#   - Rybak MJ (2006) Pharmacodynamics: Relation to Antimicrobial Resistance
#     Am J Infect Control 34(5 Suppl):S38–45
#   - Craig WA (1998) Pharmacokinetic/pharmacodynamic parameters – Clin Infect Dis
#   - Roberts KB (2011) UTI clinical practice guideline – Pediatrics 128:595
#   - Eddy AA (2014) Progression in CKD: the road less traveled – Am J Physiol
#   - Weiss R (2014) VUR and renal scarring – Pediatr Nephrol
###############################################################################

library(mrgsolve)
library(tidyverse)
library(ggplot2)

# --------------------------------------------------------------------------- #
# MODEL CODE
# --------------------------------------------------------------------------- #
cpn_code <- '
$PARAM
// -------- Ciprofloxacin PK (500 mg PO BID) --------------------------------
ka_cipro    = 1.50,   // absorption rate constant (h-1)  Craig 1998
F_cipro     = 0.70,   // oral bioavailability
Vc_cipro    = 140,    // central volume (L)  Vd=3-4 L/kg × 40 kg
Vp_cipro    = 100,    // peripheral volume (L)
k12_cipro   = 0.12,   // central→periph (h-1)
k21_cipro   = 0.08,   // periph→central (h-1)
CL_cipro    = 25,     // total clearance (L/h)  t1/2≈4h
fu_cipro    = 0.65,   // unbound fraction in plasma/renal tissue

// -------- TMP-SMX PK (160/800 mg PO BID) -----------------------------------
ka_tmp      = 0.80,   // h-1
F_tmp       = 0.95,   // bioavailability
Vc_tmp      = 85,     // L  (Vd≈1.6 L/kg)
CL_tmp      = 3.50,   // L/h  (t1/2≈10h)
ka_smx      = 0.80,
F_smx       = 0.85,
Vc_smx      = 20,     // L  (Vd≈0.2 L/kg)
CL_smx      = 1.50,   // L/h  (t1/2≈10h)
dose_tmp    = 160,    // mg
dose_smx    = 800,    // mg

// -------- Nitrofurantoin PK (100 mg PO QD, prophylaxis) --------------------
ka_nit      = 2.00,   // h-1  rapid absorption
F_nit       = 0.75,   // bioavailability
ke_nit      = 0.40,   // urinary elimination (h-1)
fuUr_nit    = 0.40,   // fraction excreted unchanged in urine
dose_nit    = 100,    // mg

// -------- PD: Ciprofloxacin (AUC/MIC-driven kill) --------------------------
MIC_cipro   = 0.25,   // μg/mL (susceptible EUCAST breakpoint for E. coli)
Emax_cipro  = 3.50,   // max kill (log10 CFU per h)
EC50_cipro  = 2.00,   // fAUC/MIC at 50% Emax  (target > 125)
Hill_cipro  = 1.50,   // Hill slope

// -------- PD: TMP-SMX (T>MIC-dependent) ------------------------------------
MIC_tmp     = 2.00,   // mg/L  TMP component  (susceptible EUCAST ≤2)
Emax_tmp    = 2.50,   // max kill (log10 CFU/h)
EC50_tmp    = 1.50,   // T>MIC (mg/L) for 50% effect
Hill_tmp    = 2.00,

// -------- PD: Nitrofurantoin (urinary Cmax-dependent) ----------------------
MIC_nit     = 32,     // μg/mL urinary (susceptible ≤32)
Emax_nit    = 2.00,   // max kill (log10 CFU/h)
EC50_nit    = 64,     // μg/mL urinary
Hill_nit    = 2.00,

// -------- Bacterial dynamics -----------------------------------------------
kgrow       = 0.15,   // intrinsic growth rate (h-1)  ~doubling time 4.6h
Bmax        = 9.00,   // carrying capacity (log10 CFU/g)
kdeath      = 0.02,   // spontaneous death rate (h-1)
kbiofilm    = 0.04,   // biofilm formation (/h)
kdisbio     = 0.01,   // biofilm dispersion (/h)
phi_biofilm = 0.30,   // fraction of biofilm kill vs. planktonic
krecurr     = 0.0005, // spontaneous re-seeding rate (h-1)  reinfection
Resist_F    = 1.00,   // resistance factor (1=susceptible, >1=resistant)

// -------- Inflammation dynamics --------------------------------------------
kNin        = 2.00,   // max neutrophil inflow stimulus
kNout       = 0.30,   // neutrophil efflux (h-1)
kMin        = 0.50,   // max macrophage inflow
kMout       = 0.05,   // macrophage efflux (h-1)
kIL6        = 0.80,   // IL-6 production coefficient
kIL6_deg    = 0.20,   // IL-6 degradation (h-1)
kTGFb1      = 0.30,   // TGF-β1 production (h-1 per macrophage unit)
kTGFb1_deg  = 0.10,   // TGF-β1 degradation (h-1)
neutKill    = 0.50,   // neutrophil bacterial kill contribution (log10/h per unit)
macKill     = 0.30,   // macrophage bacterial kill (log10/h per unit)

// -------- Fibrosis dynamics ------------------------------------------------
kCol        = 0.001,  // collagen synthesis coefficient
kCol_deg    = 0.0002, // collagen degradation (h-1)
kScar       = 0.0005, // scar formation from collagen (h-1)
ScarMax     = 0.90,   // maximum scar fraction
ColMax      = 5.00,   // max normalized collagen

// -------- GFR dynamics -----------------------------------------------------
GFR0        = 100,    // baseline GFR (mL/min/1.73m²)  patient-specific
kGFR_scar   = 5.00,   // GFR lost per unit scar (mL/min)
kGFR_inf    = 0.10,   // GFR reduction from active inflammation
kGFR_return = 0.01,   // rate of GFR recovery when scar stops growing (h-1)
GFR_floor   = 5.00,   // ESKD threshold

// -------- Simulation flags (treatment on=1 / off=0) ------------------------
use_cipro   = 0,
use_tmp     = 0,
use_nit     = 0

$INIT
Cipro_gut  = 0,
Cipro_C    = 0,
Cipro_P    = 0,
TMP_gut    = 0,
TMP_C      = 0,
NIT_gut    = 0,
NIT_urine  = 0,
Bacteria   = 6,    // log10 CFU/g (acute pyelonephritis baseline ~10^6)
Biofilm    = 0.10,
Neutrophil = 1,    // normalised (1 = resting)
Macrophage = 1,
IL6        = 1,    // normalised (1 = baseline ~1 pg/mL)
TGFb1      = 1,    // normalised
Collagen   = 1,    // normalised (1 = normal ECM)
RenalScar  = 0,    // 0 = no scar, 1 = fully scarred
GFR        = 100   // mL/min/1.73m²

$OMEGA @name PK_IIV @block
0.04 0.02 0.09  // CL_cipro, Vc_cipro covariance

$OMEGA @name PD_IIV
0.0625  // MIC_cipro log-normal η
0.16    // GFR0 IIV

$SIGMA
0.04  // proportional residual (PK)
0.01  // additive residual (disease markers)

$ODE
// =========================================================
// I.  CIPROFLOXACIN PK  (2-compartment oral)
// =========================================================
dxdt_Cipro_gut = -ka_cipro * Cipro_gut;
double Cp_cipro = Cipro_C;                        // concentration (μg/mL) = amount/Vc
dxdt_Cipro_C   =  ka_cipro * Cipro_gut / Vc_cipro
                 - (CL_cipro / Vc_cipro + k12_cipro) * Cipro_C
                 + k21_cipro * (Cipro_P / Vp_cipro);
dxdt_Cipro_P   =  k12_cipro * Vc_cipro * Cipro_C
                 - k21_cipro * Cipro_P;
double CiproFree = fu_cipro * Cp_cipro;            // free renal concentration

// =========================================================
// II. TMP-SMX PK  (1-compartment oral, independent fitting)
// =========================================================
dxdt_TMP_gut = -ka_tmp * TMP_gut;
dxdt_TMP_C   =  ka_tmp * TMP_gut / Vc_tmp  - (CL_tmp / Vc_tmp) * TMP_C;

// SMX: 1-compartment
double SMX_C_est = 0.0;  // placeholder (no ODE state for SMX; TMP drives PD)

// =========================================================
// III. NITROFURANTOIN PK  (urinary compartment model)
// =========================================================
dxdt_NIT_gut   = -ka_nit * NIT_gut;
dxdt_NIT_urine =  ka_nit * NIT_gut * fuUr_nit - ke_nit * NIT_urine;

// =========================================================
// IV.  PHARMACODYNAMIC KILL RATES
// =========================================================
// -- Ciprofloxacin: AUC/MIC index (time-averaged AUC ≈ Cp at SS) -------------
double MIC_c_eff = MIC_cipro * Resist_F;
double fAUC_MIC  = (CiproFree > 0) ? (CiproFree / MIC_c_eff) : 0.0;
double Kill_cip  = (use_cipro > 0.5) ?
    Emax_cipro * pow(fAUC_MIC, Hill_cipro) /
    (pow(EC50_cipro, Hill_cipro) + pow(fAUC_MIC, Hill_cipro)) : 0.0;

// -- TMP-SMX: T>MIC index ----------------------------------------------------
double MIC_t_eff = MIC_tmp * Resist_F;
double Kill_tmp_v = (use_tmp > 0.5) ?
    Emax_tmp * pow(TMP_C, Hill_tmp) /
    (pow(EC50_tmp * MIC_t_eff, Hill_tmp) + pow(TMP_C, Hill_tmp)) : 0.0;

// -- Nitrofurantoin: urinary Cmax-dependent ----------------------------------
double Kill_nit_v = (use_nit > 0.5) ?
    Emax_nit * pow(NIT_urine, Hill_nit) /
    (pow(EC50_nit, Hill_nit) + pow(NIT_urine, Hill_nit)) : 0.0;

// -- Combined kill (biofilm reduces penetration by phi_biofilm fraction) -----
double Kill_drug = (Kill_cip + Kill_tmp_v + Kill_nit_v)
                    * (1.0 - Biofilm * (1.0 - phi_biofilm));

// =========================================================
// V.   BACTERIAL DYNAMICS  (log10 CFU/g)
// =========================================================
double BactNorm = (Bacteria > 0) ? Bacteria / Bmax : 0.0;
double GrowthRate = kgrow * (1.0 - BactNorm);

// Immune-mediated kill (above resting level)
double ImmuneKill = neutKill * (Neutrophil - 1.0) + macKill * (Macrophage - 1.0);
if (ImmuneKill < 0.0) ImmuneKill = 0.0;

double dBact = GrowthRate - Kill_drug - kdeath - ImmuneKill + krecurr;
// Enforce lower bound
dxdt_Bacteria = (Bacteria <= 0.0 && dBact < 0.0) ? 0.0 : dBact;

// =========================================================
// VI.  BIOFILM DYNAMICS
// =========================================================
double BioStim  = Biofilm * (1.0 - Biofilm) * BactNorm;
double BioDis   = kdisbio * Biofilm * Kill_drug;
double dBiofilm = kbiofilm * BioStim - BioDis;
// clamp 0–1
if (Biofilm <= 0.0 && dBiofilm < 0.0) dBiofilm = 0.0;
if (Biofilm >= 1.0 && dBiofilm > 0.0) dBiofilm = 0.0;
dxdt_Biofilm = dBiofilm;

// =========================================================
// VII. INFLAMMATION
// =========================================================
double BactSignal = (Bacteria > 0.0) ? Bacteria / 6.0 : 0.0;

dxdt_Neutrophil = kNin  * BactSignal - kNout * Neutrophil;
dxdt_Macrophage = kMin  * BactSignal + 0.10 * (Neutrophil - 1.0)
                  - kMout * Macrophage;
dxdt_IL6        = kIL6  * BactSignal * (1.0 + 0.5 * (Neutrophil - 1.0))
                  - kIL6_deg * IL6;

// TGF-β1: M2 macrophage switch in chronic phase
double M2signal  = (Macrophage > 2.0) ? (Macrophage - 2.0) : 0.0;
dxdt_TGFb1 = kTGFb1 * (M2signal + 0.5) + 0.05 * (IL6 - 1.0)
              - kTGFb1_deg * TGFb1;

// =========================================================
// VIII. FIBROSIS
// =========================================================
dxdt_Collagen = kCol * TGFb1 * (1.0 - Collagen / ColMax)
                - kCol_deg * Collagen;

double ScarDrive = (Collagen > 1.0) ? kScar * (Collagen - 1.0) : 0.0;
dxdt_RenalScar  = ScarDrive * (1.0 - RenalScar / ScarMax);
if (dxdt_RenalScar < 0.0) dxdt_RenalScar = 0.0;

// =========================================================
// IX.  GFR DECLINE
// =========================================================
double GFR_target = GFR0 * (1.0 - kGFR_scar / GFR0 * RenalScar)
                    - kGFR_inf * (IL6 - 1.0);
if (GFR_target < GFR_floor) GFR_target = GFR_floor;
dxdt_GFR = -kGFR_return * (GFR - GFR_target);
if (GFR <= GFR_floor && dxdt_GFR < 0.0) dxdt_GFR = 0.0;

$TABLE
double Creatinine  = 100.0 / (GFR > 5.0 ? GFR : 5.0);
double CKD_Stage   = (GFR >= 90) ? 1 : (GFR >= 60) ? 2 :
                     (GFR >= 30) ? 3 : (GFR >= 15) ? 4 : 5;
double ScarPct     = RenalScar * 100.0;
double Bacteremia  = (Bacteria > 7.5) ? 1.0 : 0.0;  // flag for urosepsis risk
double fAUC_MIC_OB = CiproFree / (MIC_cipro * Resist_F);  // PK/PD index

$CAPTURE
Cipro_C, TMP_C, NIT_urine, Bacteria, Biofilm, Neutrophil, Macrophage,
IL6, TGFb1, Collagen, RenalScar, GFR, Creatinine, CKD_Stage,
ScarPct, Bacteremia, fAUC_MIC_OB, Kill_drug
'

# --------------------------------------------------------------------------- #
# Compile model
# --------------------------------------------------------------------------- #
mod <- mcode("cpn_qsp", cpn_code, quiet = TRUE)

# --------------------------------------------------------------------------- #
# Dosing helpers
# --------------------------------------------------------------------------- #
make_cipro_events <- function(dose = 500, start = 0, ndays = 14) {
  ev(amt = dose * 0.70, ii = 12, addl = ndays * 2 - 1, cmt = "Cipro_gut",
     time = start)
}
make_tmp_events <- function(dose = 160, start = 0, ndays = 14) {
  ev(amt = dose * 0.95, ii = 12, addl = ndays * 2 - 1, cmt = "TMP_gut",
     time = start)
}
make_nit_events <- function(dose = 100, start = 0, ndays = 180) {
  ev(amt = dose * 0.75, ii = 24, addl = ndays - 1, cmt = "NIT_gut",
     time = start)
}

# --------------------------------------------------------------------------- #
# SCENARIO DEFINITIONS
# --------------------------------------------------------------------------- #
scenarios <- list(

  # 1. Untreated — baseline disease progression
  S1_untreated = list(
    label   = "S1: Untreated (No Antibiotics)",
    params  = list(use_cipro = 0, use_tmp = 0, use_nit = 0),
    events  = NULL,
    color   = "#cc0000"
  ),

  # 2. Ciprofloxacin 14-day course (acute pyelonephritis)
  S2_cipro14 = list(
    label   = "S2: Ciprofloxacin 500 mg BID × 14 days",
    params  = list(use_cipro = 1, use_tmp = 0, use_nit = 0),
    events  = make_cipro_events(500, 0, 14),
    color   = "#0066cc"
  ),

  # 3. TMP-SMX 14-day course
  S3_tmp14 = list(
    label   = "S3: TMP-SMX 160/800 mg BID × 14 days",
    params  = list(use_cipro = 0, use_tmp = 1, use_nit = 0),
    events  = make_tmp_events(160, 0, 14),
    color   = "#009900"
  ),

  # 4. Ciprofloxacin 7-day (guideline-based shorter course)
  S4_cipro7 = list(
    label   = "S4: Ciprofloxacin 500 mg BID × 7 days",
    params  = list(use_cipro = 1, use_tmp = 0, use_nit = 0),
    events  = make_cipro_events(500, 0, 7),
    color   = "#6699ff"
  ),

  # 5. Nitrofurantoin 100 mg QD prophylaxis (6 months) — recurrence prevention
  S5_nit_ppx = list(
    label   = "S5: Nitrofurantoin 100 mg QD × 6-month prophylaxis",
    params  = list(use_cipro = 0, use_tmp = 0, use_nit = 1),
    events  = make_nit_events(100, 0, 180),
    color   = "#ff9900"
  ),

  # 6. Cipro 14 days + Nitrofurantoin prophylaxis (aggressive strategy)
  S6_cipro_then_ppx = list(
    label   = "S6: Cipro × 14 d then Nitrofurantoin × 6 mo",
    params  = list(use_cipro = 1, use_tmp = 0, use_nit = 1),
    events  = c(make_cipro_events(500, 0, 14),
                make_nit_events(100, 14 * 24, 180)),
    color   = "#990099"
  ),

  # 7. TMP-SMX-resistant strain (Resist_F = 4) — clinical failure
  S7_resistant = list(
    label   = "S7: TMP-SMX-Resistant Strain (MIC shift × 4)",
    params  = list(use_cipro = 0, use_tmp = 1, use_nit = 0, Resist_F = 4.0),
    events  = make_tmp_events(160, 0, 14),
    color   = "#888800"
  )
)

# --------------------------------------------------------------------------- #
# Run all scenarios (1-year simulation, hourly)
# --------------------------------------------------------------------------- #
run_scenario <- function(sc) {
  p_update <- sc$params
  e        <- sc$events

  m <- param(mod, p_update)
  if (!is.null(e)) {
    mrgsim(m, events = e, end = 8760, delta = 1, carry_out = "evid")
  } else {
    mrgsim(m, end = 8760, delta = 1)
  }
}

results <- lapply(scenarios, function(sc) {
  message("Running: ", sc$label)
  out <- run_scenario(sc)
  as.data.frame(out) |> mutate(Scenario = sc$label)
})
all_results <- bind_rows(results)

# --------------------------------------------------------------------------- #
# PLOTS
# --------------------------------------------------------------------------- #

theme_qsp <- theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 9),
    strip.background = element_rect(fill = "#e8f0fe"),
    panel.grid.minor = element_blank()
  )

col_map <- setNames(
  sapply(scenarios, `[[`, "color"),
  sapply(scenarios, `[[`, "label")
)

# -- Panel 1: Bacterial Load -------------------------------------------------
p1 <- ggplot(all_results, aes(x = time / 24, y = Bacteria, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = col_map) +
  labs(title = "Renal Bacterial Burden",
       x = "Time (days)", y = "log₁₀ CFU/g kidney tissue",
       color = NULL) +
  geom_hline(yintercept = c(3, 5), linetype = "dashed", color = "grey50") +
  annotate("text", x = 200, y = 3.2, label = "clearance threshold") +
  theme_qsp

# -- Panel 2: Ciprofloxacin PK -----------------------------------------------
pk_cipro <- all_results |>
  filter(grepl("Cipro.*14", Scenario)) |>
  filter(time <= 336)  # first 14 days

p2 <- ggplot(pk_cipro, aes(x = time, y = Cipro_C)) +
  geom_line(color = "#0066cc", linewidth = 1) +
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "red",
             linewidth = 0.8) +
  annotate("text", x = 50, y = 0.40, label = "MIC = 0.25 μg/mL",
           color = "red", size = 3.5) +
  labs(title = "Ciprofloxacin PK (500 mg BID, first 14 days)",
       x = "Time (h)", y = "Plasma concentration (μg/mL)") +
  theme_qsp

# -- Panel 3: GFR over time --------------------------------------------------
p3 <- ggplot(all_results, aes(x = time / 24, y = GFR, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = col_map) +
  scale_y_continuous(limits = c(0, 110)) +
  geom_hline(yintercept = c(60, 30, 15), linetype = "dashed",
             color = c("#ffa500", "#ff6600", "#cc0000")) +
  annotate("text", x = 300, y = 62, label = "CKD G3a", size = 3, color = "#ffa500") +
  annotate("text", x = 300, y = 32, label = "CKD G4",  size = 3, color = "#ff6600") +
  labs(title = "GFR Trajectory (1-year simulation)",
       x = "Time (days)", y = "GFR (mL/min/1.73m²)",
       color = NULL) +
  theme_qsp

# -- Panel 4: Renal Scar -----------------------------------------------------
p4 <- ggplot(all_results, aes(x = time / 24, y = ScarPct, color = Scenario)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = col_map) +
  labs(title = "Cumulative Renal Scar Formation",
       x = "Time (days)", y = "Scar fraction (%)",
       color = NULL) +
  theme_qsp

# -- Panel 5: Inflammatory Markers (IL-6, TGF-β1) ----------------------------
p5 <- all_results |>
  pivot_longer(cols = c(IL6, TGFb1), names_to = "Marker", values_to = "Level") |>
  ggplot(aes(x = time / 24, y = Level, color = Scenario, linetype = Marker)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = col_map) +
  labs(title = "Inflammation & Fibrotic Mediators",
       x = "Time (days)", y = "Normalised level (1 = baseline)",
       color = NULL, linetype = "Marker") +
  theme_qsp

# -- Panel 6: PK/PD Target Attainment ----------------------------------------
fAUC_summary <- all_results |>
  filter(grepl("Cipro.*14", Scenario)) |>
  summarise(
    fAUC_peak  = max(fAUC_MIC_OB, na.rm = TRUE),
    fAUC_trough = min(fAUC_MIC_OB[fAUC_MIC_OB > 0], na.rm = TRUE),
    fAUC_mean  = mean(fAUC_MIC_OB, na.rm = TRUE)
  )

# --------------------------------------------------------------------------- #
# PRINT SUMMARY TABLE
# --------------------------------------------------------------------------- #
summary_tbl <- all_results |>
  group_by(Scenario) |>
  summarise(
    GFR_start  = GFR[which.min(time)],
    GFR_end    = GFR[which.max(time)],
    GFR_delta  = GFR_end - GFR_start,
    Scar_final = ScarPct[which.max(time)],
    Bact_nadir = min(Bacteria, na.rm = TRUE),
    Bact_final = Bacteria[which.max(time)],
    CKD_final  = CKD_Stage[which.max(time)]
  ) |>
  mutate(across(where(is.numeric), round, 2))

cat("\n=== Chronic Pyelonephritis QSP Model — 1-Year Simulation Summary ===\n\n")
print(summary_tbl, n = 20, width = 120)

cat("\n\nKey PD indices (Ciprofloxacin 14-day course):\n")
print(fAUC_summary)

message("\nAll scenarios complete. Access results via `all_results` data frame.")
message("Plots: p1 (bacterial), p2 (PK), p3 (GFR), p4 (scar), p5 (inflammation)")

# Print plots
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
