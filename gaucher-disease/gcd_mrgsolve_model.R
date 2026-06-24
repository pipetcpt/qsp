################################################################################
## Gaucher Disease (GCD) — QSP mrgsolve ODE Model
## Glucocerebrosidase Deficiency · Lysosomal Storage · ERT / SRT PK-PD
##
## Compartments (26 ODEs):
##   Drug PK (8): ERT plasma/tissue, miglustat gut/plasma,
##                eliglustat gut/plasma, venglustat gut/plasma
##   Enzyme & Substrate (4): GBA activity (macrophage), GC macrophage,
##                           GC spleen, GC liver, GC bone marrow
##   Biomarkers (5): GL-1 plasma, lyso-GL1, chitotriosidase,
##                   serum ferritin, angiopoietin-2
##   Organ volumes (2): spleen volume, liver volume
##   Hematology (2): hemoglobin, platelets
##   Bone (3): BMD, osteoclast activity, osteoblast activity
##   Inflammation (2): cytokine burden (IL-6/TNF composite), NF-kB
##
## Scenarios:
##   S1: Natural history (no treatment)
##   S2: Imiglucerase ERT 60 U/kg Q2W IV
##   S3: Velaglucerase alfa ERT 60 U/kg Q2W IV
##   S4: Eliglustat SRT (CYP2D6 extensive metabolizer)
##   S5: Eliglustat SRT (CYP2D6 poor metabolizer — higher exposure)
##   S6: Combination: Low-dose ERT + eliglustat SRT
##
## Parameter calibration notes:
##   • ERT PK: Aerts 2003 NEJM (imiglucerase); Zimran 2010 Blood (velaglucerase)
##   • GC turnover half-life ~4-8 weeks (Charrow 2004 J Inherit Metab Dis)
##   • GL-1 response to ERT: Giraldo 2012 Mol Genet Metab
##   • Eliglustat PK: Lukina 2014 Orphanet J Rare Dis; Mistry 2015 NEJM
##   • Organ volume responses: Weinreb 2007 Blood; Cox 2008 J Inherit Metab Dis
##   • Hematologic responses: Grabowski 2004 ANN Int Med
##   • Bone responses (BMD): Wenstrup 2007 Clin Orthop Relat Res
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ============================================================
## MODEL DEFINITION
## ============================================================
gcd_model <- '
$PROB Gaucher Disease QSP Model — GBA Deficiency, ERT/SRT PK-PD

$PARAM
// ─── Patient characteristics ─────────────────────────────
BW       = 70      // body weight [kg]
TYPE     = 1       // Gaucher type: 1=non-neuro, 3=chronic neuro

// ─── Baseline disease state ──────────────────────────────
GBA0     = 5       // baseline GBA activity [nmol/h/mg] (normal ~30)
GC0      = 100     // baseline GC burden macrophage [% normal accumulation]
GL1_0    = 8.5     // baseline plasma GL-1 [μg/L] (normal <1)
LYSO_0   = 18      // baseline lyso-GL1 [ng/mL]
CHITR0   = 4800    // baseline chitotriosidase [nmol/h/mL] (normal <100)
FERRIT0  = 220     // baseline serum ferritin [μg/L]
SV0      = 2200    // baseline spleen volume [mL] (normal ~300)
LV0      = 1.65    // baseline liver volume [×normal]
HGB0     = 10.5    // baseline hemoglobin [g/dL]
PLT0     = 65      // baseline platelets [×10⁹/L]
BMD0     = -1.8    // baseline lumbar spine T-score
OC0      = 2.5     // baseline osteoclast activity [AU]
OB0      = 0.9     // baseline osteoblast activity [AU]
IL6_0    = 12      // baseline IL-6 composite cytokine [pg/mL]
NFKB0    = 1.0     // baseline NF-κB activity [AU] (1=normal)

// ─── GBA natural history ────────────────────────────────────
k_gba_deg   = 0.0035  // GBA degradation rate [1/h]
k_gba_syn   = 0.0175  // GBA synthesis rate (gives GBA0 at steady state)
k_gc_syn    = 0.028   // GC synthesis rate via GCS [AU/h]
k_gc_deg    = 0.014   // GC degradation by baseline GBA [1/h]
k_gc_prog   = 0.0002  // GC disease progression rate [1/h per AU excess]

// ─── Biomarker kinetics ────────────────────────────────────
k_gl1_out   = 0.055   // GL-1 plasma clearance [1/h]
k_lyso_out  = 0.045   // lyso-GL1 clearance [1/h]
k_chitr_out = 0.008   // chitotriosidase clearance [1/h]
k_ferr_out  = 0.012   // ferritin clearance [1/h]

// ─── Organ volume kinetics ─────────────────────────────────
k_sv_in     = 0.0005  // spleen accumulation rate driven by GC_spleen [1/h]
k_sv_out    = 0.0003  // spleen regression [1/h]
k_lv_in     = 0.0004  // liver accumulation rate [1/h]
k_lv_out    = 0.0002  // liver regression [1/h]

// ─── Hematology kinetics ───────────────────────────────────
k_hgb_prod  = 0.0018  // Hb production [g/dL/h]
k_hgb_deg   = 0.0017  // Hb degradation [1/h]
k_plt_prod  = 0.018   // platelet production [10⁹/L/h]
k_plt_deg   = 0.013   // platelet degradation [1/h]

// ─── Bone kinetics ─────────────────────────────────────────
k_bmd_form  = 0.0005  // BMD formation by OB [1/h]
k_bmd_res   = 0.0008  // BMD resorption by OC [1/h]
k_oc_stim   = 0.0012  // OC stimulation by IL-6/RANKL [1/h]
k_oc_deg    = 0.0010  // OC degradation [1/h]
k_ob_stim   = 0.0004  // OB stimulation
k_ob_deg    = 0.0008  // OB degradation [1/h]

// ─── Inflammation kinetics ─────────────────────────────────
k_il6_prod  = 0.025   // IL-6 production from NF-κB [pg/mL/h]
k_il6_deg   = 0.030   // IL-6 clearance [1/h]
k_nfkb_on   = 0.018   // NF-κB activation by GC [1/h]
k_nfkb_off  = 0.015   // NF-κB deactivation [1/h]

// ─── ERT PK (imiglucerase/velaglucerase as class) ─────────
CL_ERT   = 1.4        // clearance [L/h/kg]
V1_ERT   = 0.18       // central volume [L/kg]
Q_ERT    = 0.35       // intercompartmental clearance [L/h/kg]
V2_ERT   = 0.55       // peripheral volume [L/kg]
KM6P     = 0.006      // M6P receptor Km [U/mL]
UPTK_ERT = 0.45       // tissue uptake rate [1/h]
EMAX_ERT = 0.85       // maximal GBA activity restoration [fraction]
EC50_ERT = 0.6        // EC50 for GBA restoration [U/mL in macrophage]

// ─── Velaglucerase modifier ─────────────────────────────────
VELA_MOD = 1.0        // =1 for imiglucerase, =1.05 for velaglucerase (slightly higher Emax)

// ─── Eliglustat PK (extensive metabolizer default) ─────────
KA_ELIS  = 0.80       // absorption rate [1/h]
F_ELIS   = 0.20       // bioavailability (EM)
CL_ELIS  = 38.0       // clearance [L/h]
V_ELIS   = 106        // volume of distribution [L]
IC50_ELIS = 0.010     // GCS IC50 [μg/mL]
IMAX_ELIS = 0.95      // max GCS inhibition

// ─── Miglustat PK ──────────────────────────────────────────
KA_MIGS  = 0.60       // absorption [1/h]
F_MIGS   = 0.97       // bioavailability
CL_MIGS  = 4.5        // clearance [L/h]
V_MIGS   = 28         // volume [L]
IC50_MIGS = 50.0      // GCS IC50 [μg/mL] (low selectivity)
IMAX_MIGS = 0.80      // max inhibition

// ─── Venglustat PK ─────────────────────────────────────────
KA_VENG  = 0.50       // absorption [1/h]
F_VENG   = 0.70       // bioavailability
CL_VENG  = 22.0       // clearance [L/h]
V_VENG   = 480        // high Vd (CNS penetration) [L]
Kp_CNS   = 0.30       // brain/plasma ratio
IC50_VENG = 0.015     // GCS IC50 (lower than miglustat)
IMAX_VENG = 0.92

// ─── Scenario flags (0=off, 1=on) ─────────────────────────
USE_ERT   = 0         // imiglucerase/velaglucerase ERT
USE_VELA  = 0         // velaglucerase modifier
USE_ELIS  = 0         // eliglustat SRT
USE_MIGS  = 0         // miglustat SRT
USE_VENG  = 0         // venglustat SRT (CNS)
DOSE_ERT  = 60        // ERT dose [U/kg]
DOSE_ELIS = 84        // eliglustat dose [mg BID]
DOSE_MIGS = 100       // miglustat dose [mg TID]
DOSE_VENG = 15        // venglustat dose [mg QD]

$CMT
// ERT
ERT_C      // ERT central plasma [U/kg]
ERT_T      // ERT tissue / macrophage

// Eliglustat
ELIS_GUT   // eliglustat gut
ELIS_C     // eliglustat plasma [μg/mL]

// Miglustat
MIGS_GUT   // miglustat gut
MIGS_C     // miglustat plasma [μg/mL]

// Venglustat
VENG_GUT   // venglustat gut
VENG_C     // venglustat plasma [μg/mL]

// Disease compartments
GBA        // lysosomal GBA activity [nmol/h/mg protein]
GC_MAC     // glucocerebroside macrophage [AU]
GC_SP      // GC spleen burden [AU]
GC_LV      // GC liver burden [AU]
GC_BM      // GC bone marrow burden [AU]

// Biomarkers
GL1        // plasma GL-1 [μg/L]
LYSOGL1    // plasma lyso-GL1 [ng/mL]
CHITR      // chitotriosidase [nmol/h/mL]
FERRIT     // serum ferritin [μg/L]

// Organ volumes
SV         // spleen volume [mL]
LV         // liver volume [×normal]

// Hematology
HGB        // hemoglobin [g/dL]
PLT        // platelets [×10⁹/L]

// Bone
BMD        // lumbar spine T-score
OC         // osteoclast activity [AU]
OB         // osteoblast activity [AU]

// Inflammation
IL6        // cytokine composite (IL-6, TNF) [pg/mL]
NFKB       // NF-κB activity [AU]

$INIT
ERT_C    = 0,    ERT_T    = 0
ELIS_GUT = 0,    ELIS_C   = 0
MIGS_GUT = 0,    MIGS_C   = 0
VENG_GUT = 0,    VENG_C   = 0
GBA      = 5
GC_MAC   = 100
GC_SP    = 100
GC_LV    = 100
GC_BM    = 100
GL1      = 8.5
LYSOGL1  = 18
CHITR    = 4800
FERRIT   = 220
SV       = 2200
LV       = 1.65
HGB      = 10.5
PLT      = 65
BMD      = -1.8
OC       = 2.5
OB       = 0.9
IL6      = 12
NFKB     = 1.0

$ODE
// ────────────────────────────────────────────────────────────
// ERT PK (2-compartment IV)
// ────────────────────────────────────────────────────────────
double k10_ert = (CL_ERT * BW) / (V1_ERT * BW);
double k12_ert = Q_ERT / (V1_ERT * BW);
double k21_ert = Q_ERT / (V2_ERT * BW);

dxdt_ERT_C = -(k10_ert + k12_ert) * ERT_C + k21_ert * ERT_T;
dxdt_ERT_T = k12_ert * ERT_C - k21_ert * ERT_T;

// ────────────────────────────────────────────────────────────
// Eliglustat PK (1-compartment oral)
// ────────────────────────────────────────────────────────────
dxdt_ELIS_GUT = -KA_ELIS * ELIS_GUT;
dxdt_ELIS_C   = KA_ELIS * ELIS_GUT * F_ELIS / V_ELIS - (CL_ELIS / V_ELIS) * ELIS_C;

// ────────────────────────────────────────────────────────────
// Miglustat PK
// ────────────────────────────────────────────────────────────
dxdt_MIGS_GUT = -KA_MIGS * MIGS_GUT;
dxdt_MIGS_C   = KA_MIGS * MIGS_GUT * F_MIGS / V_MIGS - (CL_MIGS / V_MIGS) * MIGS_C;

// ────────────────────────────────────────────────────────────
// Venglustat PK
// ────────────────────────────────────────────────────────────
dxdt_VENG_GUT = -KA_VENG * VENG_GUT;
dxdt_VENG_C   = KA_VENG * VENG_GUT * F_VENG / V_VENG - (CL_VENG / V_VENG) * VENG_C;

// ────────────────────────────────────────────────────────────
// Drug effects on GBA and GCS
// ────────────────────────────────────────────────────────────
// ERT effect: restores GBA activity via M6P-receptor uptake
double ERT_effect = USE_ERT * EMAX_ERT * VELA_MOD * ERT_T / (EC50_ERT + ERT_T);

// SRT effects: inhibit GCS → reduce substrate synthesis
double ELIS_inh = USE_ELIS * IMAX_ELIS * ELIS_C / (IC50_ELIS + ELIS_C);
double MIGS_inh = USE_MIGS * IMAX_MIGS * MIGS_C / (IC50_MIGS + MIGS_C);
double VENG_inh = USE_VENG * IMAX_VENG * VENG_C / (IC50_VENG + VENG_C);
double SRT_inh  = 1.0 - (1.0 - ELIS_inh) * (1.0 - MIGS_inh) * (1.0 - VENG_inh);

// Effective GBA: baseline mutant activity + ERT restoration
double GBA_eff = GBA + ERT_effect * (30.0 - GBA); // 30 = normal GBA

// ────────────────────────────────────────────────────────────
// GBA enzyme dynamics
// ────────────────────────────────────────────────────────────
dxdt_GBA = k_gba_syn - k_gba_deg * GBA + ERT_effect * 0.005;

// ────────────────────────────────────────────────────────────
// Glucocerebroside dynamics
// ────────────────────────────────────────────────────────────
double GCS_activity = k_gc_syn * (1.0 - SRT_inh);
double GC_clearance = (k_gc_deg * GBA_eff / GBA0 + ERT_effect * 0.025);

dxdt_GC_MAC = GCS_activity - GC_clearance * GC_MAC;
dxdt_GC_SP  = 0.15 * GC_MAC - 0.018 * (GBA_eff / GBA0) * GC_SP;
dxdt_GC_LV  = 0.12 * GC_MAC - 0.016 * (GBA_eff / GBA0) * GC_LV;
dxdt_GC_BM  = 0.18 * GC_MAC - 0.014 * (GBA_eff / GBA0) * GC_BM;

// ────────────────────────────────────────────────────────────
// NF-κB and inflammation
// ────────────────────────────────────────────────────────────
double GC_excess = (GC_MAC > GC0) ? (GC_MAC - GC0) / GC0 : 0.0;
dxdt_NFKB = k_nfkb_on * GC_excess - k_nfkb_off * NFKB;
dxdt_IL6  = k_il6_prod * NFKB - k_il6_deg * IL6;

// ────────────────────────────────────────────────────────────
// Biomarkers
// ────────────────────────────────────────────────────────────
double gl1_drive  = 0.045 * GC_MAC;   // GC → GL-1 production
double lyso_drive = 0.025 * GC_MAC;
double chitr_drive = 3.5  * IL6;      // macrophage activation
double ferr_drive  = 0.55 * IL6;

dxdt_GL1    = gl1_drive   - k_gl1_out  * GL1;
dxdt_LYSOGL1= lyso_drive  - k_lyso_out * LYSOGL1;
dxdt_CHITR  = chitr_drive - k_chitr_out * CHITR;
dxdt_FERRIT = ferr_drive  - k_ferr_out * FERRIT;

// ────────────────────────────────────────────────────────────
// Organ volumes
// ────────────────────────────────────────────────────────────
double sv_drive = k_sv_in * GC_SP;
double sv_regress = k_sv_out * (SV - 300.0) * (GBA_eff / GBA0);  // 300 mL = normal
dxdt_SV = sv_drive - sv_regress;

double lv_drive = k_lv_in * GC_LV;
double lv_regress = k_lv_out * (LV - 1.0) * (GBA_eff / GBA0);
dxdt_LV = lv_drive - lv_regress;

// ────────────────────────────────────────────────────────────
// Hematology — suppressed by hypersplenism & bone marrow
// ────────────────────────────────────────────────────────────
double spleen_factor = (SV > 300) ? 300 / SV : 1.0;  // sequestration
double bm_factor = (GC_BM > 50) ? exp(-0.008 * (GC_BM - 50)) : 1.0;

dxdt_HGB = k_hgb_prod * bm_factor * spleen_factor - k_hgb_deg * HGB;
dxdt_PLT = k_plt_prod * bm_factor * spleen_factor - k_plt_deg * PLT;

// ────────────────────────────────────────────────────────────
// Bone — OB/OC balance, IL-6/RANKL driven
// ────────────────────────────────────────────────────────────
dxdt_OC = k_oc_stim * IL6 - k_oc_deg * OC;
dxdt_OB = k_ob_stim / (1.0 + 0.3 * IL6) - k_ob_deg * OB;
dxdt_BMD = k_bmd_form * OB - k_bmd_res * OC;

$TABLE
double GC_TOTAL = GC_MAC + GC_SP + GC_LV + GC_BM;
double GBA_PCT_NORMAL = 100.0 * GBA / 30.0;  // % of normal GBA activity
double ERT_Cplasma = ERT_C * BW;              // U total
double SRT_GCS_inh = SRT_inh * 100.0;        // % GCS inhibition

$CAPTURE
GBA GBA_PCT_NORMAL GC_MAC GC_SP GC_LV GC_BM GC_TOTAL
GL1 LYSOGL1 CHITR FERRIT
SV LV HGB PLT BMD OC OB IL6 NFKB
ERT_C ERT_T ELIS_C MIGS_C VENG_C
ERT_Cplasma SRT_GCS_inh
'

mod <- mcode("gcd_qsp", gcd_model)

## ============================================================
## DOSING HELPERS
## ============================================================
make_ert_doses <- function(bw = 70, dose_ukg = 60, n_doses = 52, interval_h = 336) {
  # Q2W = 336 h; n_doses = 52 → ~2 years
  ev(cmt = "ERT_C", amt = dose_ukg * bw / (0.18 * bw), # normalised to conc
     ii = interval_h, addl = n_doses - 1, rate = -2) # infuse over 1 h (rate=-2)
}

make_oral_doses <- function(cmt, dose_mg, freq_h, n_days = 730) {
  n_doses <- floor(n_days * 24 / freq_h)
  ev(cmt = cmt, amt = dose_mg, ii = freq_h, addl = n_doses - 1)
}

## ============================================================
## SCENARIO DEFINITIONS
## ============================================================
scenarios <- list(
  S1_NaturalHistory = list(
    label = "S1: Natural History",
    param = list(USE_ERT=0, USE_ELIS=0, USE_MIGS=0, USE_VENG=0),
    events = NULL
  ),
  S2_Imiglucerase = list(
    label = "S2: Imiglucerase 60 U/kg Q2W",
    param = list(USE_ERT=1, USE_VELA=0, DOSE_ERT=60),
    events = make_ert_doses(70, 60, 52, 336)
  ),
  S3_Velaglucerase = list(
    label = "S3: Velaglucerase α 60 U/kg Q2W",
    param = list(USE_ERT=1, USE_VELA=1, VELA_MOD=1.05, DOSE_ERT=60),
    events = make_ert_doses(70, 60, 52, 336)
  ),
  S4_Eliglustat_EM = list(
    label = "S4: Eliglustat 84 mg BID (CYP2D6 EM)",
    param = list(USE_ELIS=1, F_ELIS=0.20, DOSE_ELIS=84),
    events = make_oral_doses("ELIS_GUT", 84, 12, 730)
  ),
  S5_Eliglustat_PM = list(
    label = "S5: Eliglustat 84 mg QD (CYP2D6 PM — higher AUC)",
    param = list(USE_ELIS=1, F_ELIS=0.35, CL_ELIS=8.0, DOSE_ELIS=84),
    events = make_oral_doses("ELIS_GUT", 84, 24, 730)
  ),
  S6_CombinationLowERT_ELIS = list(
    label = "S6: Low-dose ERT (30 U/kg) + Eliglustat (84 mg BID)",
    param = list(USE_ERT=1, USE_ELIS=1, DOSE_ERT=30, DOSE_ELIS=84),
    events = c(make_ert_doses(70, 30, 52, 336),
               make_oral_doses("ELIS_GUT", 84, 12, 730))
  )
)

## ============================================================
## RUN ALL SCENARIOS
## ============================================================
sim_end <- 730 * 24   # 2 years in hours
sim_delta <- 24       # sample every 24 h

run_scenario <- function(sc) {
  m <- param(mod, sc$param)
  if (is.null(sc$events)) {
    out <- mrgsim(m, end = sim_end, delta = sim_delta)
  } else {
    out <- mrgsim(m, events = sc$events, end = sim_end, delta = sim_delta)
  }
  as.data.frame(out) %>% mutate(scenario = sc$label, time_days = time / 24)
}

results <- bind_rows(lapply(scenarios, run_scenario))

## ============================================================
## VISUALISATION HELPERS
## ============================================================
plot_outcomes <- function(data, var, ylab, title, yline = NULL) {
  p <- ggplot(data, aes(x = time_days, y = .data[[var]], color = scenario)) +
    geom_line(linewidth = 1.1) +
    labs(x = "Time (days)", y = ylab, title = title,
         color = "Scenario") +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom", legend.text = element_text(size = 8)) +
    guides(color = guide_legend(ncol = 2))
  if (!is.null(yline)) p <- p + geom_hline(yintercept = yline, linetype = "dashed", color = "black")
  p
}

## ============================================================
## KEY PLOTS
## ============================================================
p1 <- plot_outcomes(results, "GL1",    "Plasma GL-1 (μg/L)",          "Plasma Glucosylceramide (GL-1)", yline = 1.0)
p2 <- plot_outcomes(results, "LYSOGL1","Lyso-GL1 (ng/mL)",            "Lyso-Glucosylceramide", yline = 1.5)
p3 <- plot_outcomes(results, "CHITR",  "Chitotriosidase (nmol/h/mL)", "Chitotriosidase Activity", yline = 100)
p4 <- plot_outcomes(results, "SV",     "Spleen Volume (mL)",          "Spleen Volume", yline = 500)
p5 <- plot_outcomes(results, "LV",     "Liver Volume (×normal)",      "Liver Volume", yline = 1.1)
p6 <- plot_outcomes(results, "HGB",    "Hemoglobin (g/dL)",           "Hemoglobin", yline = 12)
p7 <- plot_outcomes(results, "PLT",    "Platelets (×10⁹/L)",          "Platelets", yline = 100)
p8 <- plot_outcomes(results, "BMD",    "Lumbar Spine T-score",        "Bone Mineral Density", yline = -2.5)
p9 <- plot_outcomes(results, "GBA_PCT_NORMAL", "GBA Activity (% Normal)", "GBA Enzyme Activity")
p10 <- plot_outcomes(results, "SRT_GCS_inh",   "GCS Inhibition (%)",      "Glucosylceramide Synthase Inhibition")

## Print a summary at 6, 12, 24 months
summary_times <- c(182, 365, 548, 730)
summary_df <- results %>%
  filter(time_days %in% summary_times) %>%
  select(scenario, time_days, GL1, LYSOGL1, CHITR, SV, LV, HGB, PLT, BMD, GBA_PCT_NORMAL) %>%
  arrange(scenario, time_days)

cat("\n==== Gaucher Disease QSP Model — Scenario Summary ====\n")
print(summary_df, n = Inf)

## Percent change from baseline at 12 months
baseline_df <- results %>% filter(time_days == 0) %>%
  select(scenario, GL1, LYSOGL1, CHITR, SV, HGB, PLT) %>%
  rename_with(~paste0(.x, "_BL"), -scenario)

response_12m <- results %>%
  filter(time_days == 365) %>%
  left_join(baseline_df, by = "scenario") %>%
  transmute(
    scenario,
    GL1_pct     = (GL1 - GL1_BL)    / GL1_BL    * 100,
    LYSOGL1_pct = (LYSOGL1 - LYSOGL1_BL) / LYSOGL1_BL * 100,
    CHITR_pct   = (CHITR - CHITR_BL) / CHITR_BL * 100,
    SV_pct      = (SV - SV_BL)       / SV_BL    * 100,
    HGB_delta   = HGB - HGB_BL,
    PLT_pct     = (PLT - PLT_BL)     / PLT_BL   * 100
  )

cat("\n==== 12-Month Treatment Response ====\n")
print(response_12m)

## ============================================================
## VIRTUAL PATIENT POPULATION (Monte Carlo)
## ============================================================
set.seed(42)
n_vp <- 200

vp_params <- tibble(
  ID       = 1:n_vp,
  BW       = rnorm(n_vp, 70, 12),
  GBA0     = pmax(rnorm(n_vp, 5, 1.5), 0.5),   # severe: lower GBA
  GC0      = pmax(rnorm(n_vp, 100, 25), 30),
  GL1_0    = pmax(rnorm(n_vp, 8.5, 2.5), 1),
  SV0      = pmax(rnorm(n_vp, 2200, 400), 400),
  HGB0     = pmax(rnorm(n_vp, 10.5, 1.2), 7),
  PLT0     = pmax(rnorm(n_vp, 65, 18), 20),
  F_ELIS   = pmin(pmax(rnorm(n_vp, 0.20, 0.06), 0.05), 0.50)  # inter-individual CYP2D6
)

run_vp <- function(id, vp_row, sc_param, events) {
  m <- param(mod, c(as.list(vp_row[-1]), sc_param))
  init(m, list(
    GBA = vp_row$GBA0, GC_MAC = vp_row$GC0, GC_SP = vp_row$GC0,
    GC_LV = vp_row$GC0 * 0.8, GC_BM = vp_row$GC0 * 1.1,
    GL1 = vp_row$GL1_0, SV = vp_row$SV0,
    HGB = vp_row$HGB0, PLT = vp_row$PLT0
  ))
  if (is.null(events)) {
    out <- mrgsim(m, end = sim_end, delta = sim_delta)
  } else {
    out <- mrgsim(m, events = events, end = sim_end, delta = sim_delta)
  }
  as.data.frame(out) %>%
    mutate(ID = id, time_days = time / 24) %>%
    filter(time_days %in% c(0, 182, 365, 730))
}

## Run VP population for imiglucerase scenario
vp_ert_results <- bind_rows(
  mapply(run_vp,
         id     = vp_params$ID,
         vp_row = split(vp_params, seq_len(nrow(vp_params))),
         MoreArgs = list(
           sc_param = list(USE_ERT=1),
           events   = make_ert_doses(70, 60, 52, 336)
         ),
         SIMPLIFY = FALSE)
)

cat("\n==== Virtual Population Hb Response at 12 months (ERT) ====\n")
vp_resp <- vp_ert_results %>%
  filter(time_days %in% c(0, 365)) %>%
  select(ID, time_days, HGB) %>%
  pivot_wider(names_from = time_days, values_from = HGB, names_prefix = "t") %>%
  mutate(delta_HGB = t365 - t0, responder = delta_HGB >= 1.0)

cat(sprintf("  Hb response rate (≥1 g/dL): %.1f%%\n",
            mean(vp_resp$responder, na.rm = TRUE) * 100))
cat(sprintf("  Median delta Hb: %.2f g/dL (IQR: %.2f–%.2f)\n",
            median(vp_resp$delta_HGB, na.rm = TRUE),
            quantile(vp_resp$delta_HGB, 0.25, na.rm = TRUE),
            quantile(vp_resp$delta_HGB, 0.75, na.rm = TRUE)))

cat("\n==== Model complete. Use Shiny app for interactive exploration. ====\n")
