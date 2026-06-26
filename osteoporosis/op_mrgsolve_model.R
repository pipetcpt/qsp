## ============================================================
## Osteoporosis QSP Model — mrgsolve ODE
## Disease: Postmenopausal Osteoporosis / GIOP
## Model states: 22 ODE compartments
##
## Compartments:
##   PK:
##     [1]  ALN_BONE   — Alendronate bone depot (nmol/g)
##     [2]  ZOL_BONE   — Zoledronate bone depot (nmol/g)
##     [3]  DMAB_SC    — Denosumab SC depot (µg)
##     [4]  DMAB_C     — Denosumab central (µg/mL)
##     [5]  TPTD_C     — Teriparatide central (pg/mL)
##     [6]  ROMO_SC    — Romosozumab SC depot (µg)
##     [7]  ROMO_C     — Romosozumab central (µg/mL)
##     [8]  ROMO_P     — Romosozumab peripheral (µg/mL)
##   Bone Cells:
##     [9]  OB         — Active osteoblast (cells/mm³) [relative to baseline]
##     [10] OC         — Active osteoclast (cells/mm³) [relative to baseline]
##     [11] PREOB      — Pre-osteoblast pool
##     [12] PREOC      — Pre-osteoclast pool
##   Bone Matrix:
##     [13] BMD        — Bone mineral density (g/cm², lumbar spine proxy)
##     [14] CTX        — CTX bone resorption marker (ng/mL)
##     [15] P1NP       — P1NP bone formation marker (µg/L)
##   Mediators:
##     [16] RANKL      — Soluble RANKL (pmol/L)
##     [17] OPG        — Osteoprotegerin (pmol/L)
##     [18] SCLER      — Sclerostin (pmol/L)
##     [19] PTH_sys    — PTH 1-84 systemic (pg/mL)
##     [20] E2         — Estradiol (pg/mL)
##     [21] Ca_sys     — Serum calcium (mg/dL)
##     [22] FRACT_RISK — Cumulative fracture probability (0–1)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ---- Model code block ----
code <- '
$PROB Osteoporosis QSP Model — mrgsolve
Bone remodeling with RANK/RANKL/OPG axis, PTH, estrogen, sclerostin
Drugs: Alendronate, Zoledronic acid, Denosumab, Teriparatide, Romosozumab

$PARAM
// --- Drug dosing (switches, 0=off, 1=on) ---
USE_ALN  = 0   // Alendronate 70mg/week PO
USE_ZOL  = 0   // Zoledronic acid 5mg/year IV
USE_DMAB = 0   // Denosumab 60mg/6mo SC
USE_TPTD = 0   // Teriparatide 20µg/day SC
USE_ROMO = 0   // Romosozumab 210mg/month SC
MENO     = 1   // Menopause (1=yes, 0=no)
GIO      = 0   // Glucocorticoid-induced OP (1=yes)

// --- Drug PK parameters ---
// Alendronate
ALN_dose   = 70000    // nmol/week oral
ALN_F      = 0.007    // bioavailability 0.7%
ALN_ka_bone= 0.05     // /h bone uptake from plasma
ALN_kout   = 0.0000004// /h bone release (t½ >10yr)
ALN_FPP_IC50 = 50     // nmol/g half-maximal FPP inhibition

// Zoledronic acid
ZOL_dose   = 5000     // nmol/year IV (5mg ≈ 10.4µmol; model in relative units)
ZOL_ka_bone= 0.10     // /h bone uptake
ZOL_kout   = 0.0001   // /h bone release (t½ ~1yr)
ZOL_FPP_IC50 = 10     // nmol/g (more potent)

// Denosumab
DMAB_dose  = 60       // µg SC every 6 months
DMAB_ka    = 0.005    // /h SC absorption
DMAB_CL    = 0.0019   // L/h central clearance
DMAB_Vc    = 2.46     // L central volume
DMAB_RANKL_IC50 = 0.003  // µg/mL (Kd ~3 pM → ~0.44ng/mL; model units)
DMAB_Emax  = 0.98

// Teriparatide
TPTD_dose  = 20000    // pg/day SC = 20µg/day
TPTD_ka    = 5.0      // /h SC absorption
TPTD_ke    = 1.4      // /h central elimination (t½≈0.5h)
TPTD_Vc    = 1500     // mL (1.5L)
TPTD_PTH1R_EC50 = 100 // pg/mL
TPTD_Emax_OB = 0.8    // max fold-increase OB
TPTD_Emax_OC = 0.3    // transient OC increase

// Romosozumab
ROMO_dose  = 210      // µg SC monthly
ROMO_ka    = 0.01     // /h SC absorption
ROMO_CL    = 0.015    // L/h
ROMO_Vc    = 3.5      // L central
ROMO_Vp    = 2.5      // L peripheral
ROMO_Q     = 0.04     // L/h inter-compartmental
ROMO_SOST_IC50 = 0.15 // µg/mL
ROMO_Emax_OB   = 1.2  // max fold-OB increase (anabolic)
ROMO_Emax_AC   = 0.6  // max fold-OC decrease (anti-resorptive)

// --- Bone cell biology ---
OB_ss     = 1.0    // normalized osteoblast steady state
OC_ss     = 1.0    // normalized osteoclast steady state
PREOB_ss  = 2.0    // pre-OB steady state
PREOC_ss  = 2.0    // pre-OC steady state

kOB_in    = 0.02   // /h OB formation from PREOB
kOB_out   = 0.004  // /h OB apoptosis (t½~7d)
kOC_in    = 0.015  // /h OC formation from PREOC
kOC_out   = 0.01   // /h OC apoptosis (t½~3d)
kPREOB_in = 0.015  // /h PREOB production (MSC)
kPREOC_in = 0.012  // /h PREOC production (HSC)

// --- RANKL/OPG/Sclerostin ---
RANKL_ss  = 10.0   // pmol/L baseline sRANKL
OPG_ss    = 40.0   // pmol/L baseline OPG
SCLER_ss  = 35.0   // pmol/L baseline sclerostin
kRANKL_in = 0.10   // /h production rate
kRANKL_out= 0.10   // /h degradation
kOPG_in   = 0.08   // /h production
kOPG_out  = 0.08   // /h degradation
kSCLER_in = 0.05   // /h sclerostin production
kSCLER_out= 0.05   // /h sclerostin degradation

// --- Hormones ---
PTH_ss    = 35.0   // pg/mL baseline PTH
E2_ss_pre = 80.0   // pg/mL premenopausal
E2_ss_post= 12.0   // pg/mL postmenopausal
Ca_ss     = 9.4    // mg/dL normal serum calcium
kPTH_in   = 3.5    // pg/mL/h
kPTH_out  = 0.10   // /h

// --- Menopause effect coefficients ---
MENO_RANKL_fold = 1.6    // RANKL increases 60% post-menopause
MENO_OPG_fold   = 0.7    // OPG decreases 30%
MENO_OB_fold    = 0.85   // OB activity decreases 15%

// --- BMD dynamics ---
BMD_ss    = 1.0    // normalized baseline (=0.85 g/cm² LS)
BMD_form  = 0.002  // /h bone formation rate (OB-driven)
BMD_resor = 0.002  // /h bone resorption rate (OC-driven)

// --- CTX / P1NP dynamics ---
CTX_ss    = 0.35   // ng/mL
P1NP_ss   = 50     // µg/L
kCTX_in   = 0.07   // /h
kCTX_out  = 0.20   // /h
kP1NP_in  = 1.5    // µg/L/h
kP1NP_out = 0.03   // /h

// --- GIO parameters ---
GIO_OB_supp  = 0.5  // 50% OB suppression
GIO_OC_stim  = 1.3  // 30% OC stimulation
GIO_RANKL_up = 1.4

// --- Fracture risk ---
FxRisk_k  = 0.001   // /h base rate (at BMD=1)
FxRisk_BMD_slope = 2.0  // BMD sensitivity: 1 SD ↓ ≈ 2x risk

$CMT ALN_BONE ZOL_BONE DMAB_SC DMAB_C TPTD_C ROMO_SC ROMO_C ROMO_P
     OB OC PREOB PREOC BMD CTX P1NP RANKL OPG SCLER PTH_sys E2 Ca_sys FRACT_RISK

$INIT
ALN_BONE = 0, ZOL_BONE = 0, DMAB_SC = 0, DMAB_C = 0, TPTD_C = 0,
ROMO_SC  = 0, ROMO_C   = 0, ROMO_P  = 0,
OB       = 1.0,  OC    = 1.0, PREOB = 2.0, PREOC = 2.0,
BMD      = 1.0,
CTX      = 0.35, P1NP  = 50.0,
RANKL    = 10.0, OPG   = 40.0, SCLER = 35.0,
PTH_sys  = 35.0, E2    = 80.0, Ca_sys = 9.4,
FRACT_RISK = 0

$ODE

// ================================================================
// DRUG PK
// ================================================================

// --- Alendronate: bone depot (simplified — PO weekly)
// Plasma transient, mostly excreted; bone uptake modeled directly
double ALN_EFF = (USE_ALN==1) ? ALN_BONE/(ALN_BONE + ALN_FPP_IC50) : 0;
dxdt_ALN_BONE = - ALN_kout*ALN_BONE;

// --- Zoledronic acid: bone depot
double ZOL_EFF = (USE_ZOL==1) ? ZOL_BONE/(ZOL_BONE + ZOL_FPP_IC50) : 0;
dxdt_ZOL_BONE = - ZOL_kout*ZOL_BONE;

// --- Denosumab: SC→Central (2-cmt simplified)
double DMAB_abs = (USE_DMAB==1) ? DMAB_ka*DMAB_SC : 0;
dxdt_DMAB_SC = -DMAB_abs;
dxdt_DMAB_C  = DMAB_abs/DMAB_Vc - (DMAB_CL/DMAB_Vc)*DMAB_C;
double DMAB_inh_RANKL = (USE_DMAB==1) ? DMAB_Emax*DMAB_C/(DMAB_C + DMAB_RANKL_IC50) : 0;

// --- Teriparatide: SC→Central (1-cmt with ka)
double TPTD_abs = (USE_TPTD==1) ? (TPTD_dose/(24.0))*TPTD_ka*exp(-TPTD_ka*SOLVERTIME) : 0;
dxdt_TPTD_C  = TPTD_abs/TPTD_Vc - TPTD_ke*TPTD_C;
double TPTD_stim = (USE_TPTD==1) ? TPTD_C/(TPTD_C + TPTD_PTH1R_EC50) : 0;

// --- Romosozumab: SC→Central→Peripheral (2-cmt)
double ROMO_abs = (USE_ROMO==1) ? ROMO_ka*ROMO_SC : 0;
dxdt_ROMO_SC = -ROMO_abs;
dxdt_ROMO_C  = ROMO_abs/ROMO_Vc - (ROMO_CL/ROMO_Vc)*ROMO_C
               - (ROMO_Q/ROMO_Vc)*ROMO_C + (ROMO_Q/ROMO_Vp)*ROMO_P;
dxdt_ROMO_P  = (ROMO_Q/ROMO_Vc)*ROMO_C - (ROMO_Q/ROMO_Vp)*ROMO_P;
double ROMO_inh_SOST = (USE_ROMO==1) ? ROMO_Emax_OB*ROMO_C/(ROMO_C + ROMO_SOST_IC50) : 0;
double ROMO_inh_OC   = (USE_ROMO==1) ? ROMO_Emax_AC*ROMO_C/(ROMO_C + ROMO_SOST_IC50) : 0;

// ================================================================
// MEDIATORS
// ================================================================

// Estrogen: determined by menopause status
double E2_target = (MENO==1) ? E2_ss_post : E2_ss_pre;
dxdt_E2 = 0.1*(E2_target - E2);

// PTH: regulated by Ca²⁺ (simplified negative feedback)
double PTH_target = PTH_ss * (Ca_ss/Ca_sys);
dxdt_PTH_sys = 0.5*(PTH_target - PTH_sys);

// Calcium homeostasis (simplified)
double Ca_OC_release  = 0.05*(OC - 1.0);  // net release from bone resorption
double Ca_PTH_kidney  = 0.01*(PTH_sys - PTH_ss);
dxdt_Ca_sys = -0.2*(Ca_sys - Ca_ss) + Ca_OC_release + Ca_PTH_kidney;

// RANKL: produced by OB, osteocytes; regulated by estrogen, GIO, denosumab
double RANKL_meno_factor = (MENO==1) ? MENO_RANKL_fold : 1.0;
double RANKL_GIO_factor  = (GIO==1)  ? GIO_RANKL_up    : 1.0;
double PTH_RANKL_up = 1 + 0.3*(PTH_sys/PTH_ss - 1);  // PTH upregulates RANKL
double E2_RANKL_down= (E2/E2_ss_pre);                 // estrogen suppresses
double RANKL_in_eff = kRANKL_in*RANKL_ss*RANKL_meno_factor*RANKL_GIO_factor*PTH_RANKL_up/E2_RANKL_down;
dxdt_RANKL = RANKL_in_eff - kRANKL_out*RANKL*(1 - DMAB_inh_RANKL);

// OPG: produced by OB, T-regs; estrogen-dependent
double OPG_meno_factor  = (MENO==1) ? MENO_OPG_fold : 1.0;
double OPG_GIO_factor   = (GIO==1)  ? 0.7 : 1.0;
double OPG_E2_factor    = (E2/E2_ss_pre);
double OPG_ROMO_up      = 1.0 + 0.5*ROMO_inh_OC;  // romosozumab↑OPG
dxdt_OPG = kOPG_in*OPG_ss*OPG_meno_factor*OPG_GIO_factor*OPG_E2_factor*OPG_ROMO_up - kOPG_out*OPG;

// Sclerostin: produced by osteocytes; inhibited by PTH (intermittent), romosozumab
double SOST_PTH_inh = (USE_TPTD==1) ? (1 - 0.3*TPTD_stim) : 1.0; // TPTD transiently ↓SOST
double SOST_ROMO_inh = (1 - ROMO_inh_SOST);
dxdt_SCLER = kSCLER_in*SCLER_ss*SOST_PTH_inh - kSCLER_out*SCLER*SOST_ROMO_inh;

// ================================================================
// BONE CELLS — coupled RANKL/OPG model
// ================================================================

double RANKL_OPG_ratio = RANKL/(OPG + 0.001);
double RANKL_OPG_ss    = RANKL_ss/OPG_ss;

// Wnt signaling (drives OB): inhibited by sclerostin
double Wnt_activity = 1 / (1 + SCLER/SCLER_ss);

// Pre-OB production: Wnt-dependent, PTH anabolic (intermittent)
double TPTD_OB_stim = (USE_TPTD==1) ? (1 + TPTD_Emax_OB*TPTD_stim) : 1.0;
double ROMO_OB_stim = (1 + ROMO_inh_SOST);
double GIO_OB_factor= (GIO==1) ? GIO_OB_supp : 1.0;
double MENO_OB_fact = (MENO==1) ? MENO_OB_fold : 1.0;
dxdt_PREOB = kPREOB_in*PREOB_ss*Wnt_activity*TPTD_OB_stim*ROMO_OB_stim*GIO_OB_factor*MENO_OB_fact - kOB_in*PREOB;
dxdt_OB    = kOB_in*PREOB - kOB_out*OB;

// Pre-OC production: RANKL/OPG-dependent; inhibited by E2
double RANKL_ratio_norm = RANKL_OPG_ratio / RANKL_OPG_ss;
double GIO_OC_factor    = (GIO==1) ? GIO_OC_stim : 1.0;
double TPTD_OC_stim     = (USE_TPTD==1) ? (1 + TPTD_Emax_OC*TPTD_stim) : 1.0; // coupling
dxdt_PREOC = kPREOC_in*PREOC_ss*RANKL_ratio_norm*GIO_OC_factor*TPTD_OC_stim - kOC_in*PREOC;

// OC: formed from PREOC; inhibited by bisphosphonate, denosumab, romosozumab
double BP_EFF_OC    = (ALN_EFF + ZOL_EFF)*0.8; // combined BP effect on OC apoptosis
double DMAB_OC_inh  = DMAB_inh_RANKL;           // denosumab: OC ↓ via RANKL block
double ROMO_OC_eff  = ROMO_inh_OC;
double kOC_out_eff  = kOC_out*(1 + BP_EFF_OC + DMAB_OC_inh + ROMO_OC_eff);
dxdt_OC    = kOC_in*PREOC - kOC_out_eff*OC;

// ================================================================
// BONE MINERAL DENSITY
// ================================================================
// BMD changes driven by OB-OC balance
double BMD_form_rate  = BMD_form * OB;
double BMD_resor_rate = BMD_resor * OC;
dxdt_BMD = BMD_form_rate - BMD_resor_rate;

// ================================================================
// BONE TURNOVER MARKERS
// ================================================================
// CTX (resorption): driven by OC activity
double CTX_in_eff = kCTX_in * CTX_ss * OC;
dxdt_CTX = CTX_in_eff - kCTX_out*CTX;

// P1NP (formation): driven by OB activity
double P1NP_in_eff = kP1NP_in * OB;
dxdt_P1NP = P1NP_in_eff - kP1NP_out*P1NP;

// ================================================================
// FRACTURE RISK (simplified Poisson hazard)
// ================================================================
// Risk increases exponentially as BMD falls below 1.0
double BMD_SD = (1.0 - BMD)/0.1;   // every 0.1 unit ≈ 1 SD at LS
double FxRisk_h = FxRisk_k * exp(FxRisk_BMD_slope * BMD_SD);
dxdt_FRACT_RISK = FxRisk_h * (1.0 - FRACT_RISK);

$TABLE
capture BMD_gcm2    = BMD * 0.85;     // absolute LS BMD (g/cm²)
capture T_score_val = (BMD_gcm2 - 0.955) / 0.120; // T-score (norm: 0.955±0.120)
capture CTX_ngl     = CTX;
capture P1NP_ugl    = P1NP;
capture RANKL_pm    = RANKL;
capture OPG_pm      = OPG;
capture SCLER_pm    = SCLER;
capture OB_rel      = OB;
capture OC_rel      = OC;
capture E2_pg       = E2;
capture PTH_pg      = PTH_sys;
capture Ca_mgdl     = Ca_sys;
capture FxRisk10yr  = 1 - exp(-FxRisk_k * exp(FxRisk_BMD_slope*(1-BMD)/0.1)*10*365*24);
capture DMAB_C_ug   = DMAB_C;
capture TPTD_C_pg   = TPTD_C;
capture ROMO_C_ug   = ROMO_C;
'

## ---- Compile model ----
mod <- mrgsolve::mcode("osteoporosis_qsp", code)

## ---- Helper: dosing event builder ----
build_events <- function(
    use_aln  = FALSE, use_zol = FALSE,
    use_dmab = FALSE, use_tptd = FALSE, use_romo = FALSE,
    dur_years = 3
) {
  ev_list <- list()
  total_h  <- dur_years * 365 * 24

  # Alendronate 70mg/week oral → directly load ALN_BONE (simplified)
  if (use_aln) {
    wk_times <- seq(0, total_h - 168, by = 168)
    ev_list$aln <- ev(cmt = "ALN_BONE", amt = 1500, time = wk_times)
  }
  # Zoledronic acid 5mg/year IV → directly load ZOL_BONE
  if (use_zol) {
    yr_times <- seq(0, total_h - 8760, by = 8760)
    ev_list$zol <- ev(cmt = "ZOL_BONE", amt = 5000, time = yr_times)
  }
  # Denosumab 60mg SC every 6 months
  if (use_dmab) {
    mo6_times <- seq(0, total_h - 4380, by = 4380)
    ev_list$dmab <- ev(cmt = "DMAB_SC", amt = 60, time = mo6_times)
  }
  # Teriparatide 20µg/day SC — continuous infusion model via rate
  if (use_tptd) {
    ev_list$tptd <- ev(cmt = "TPTD_C", amt = 20000*dur_years*365, rate = 20000/24,
                       time = 0, tinf = dur_years*365*24)
  }
  # Romosozumab 210mg SC monthly (12 months, then switch)
  if (use_romo) {
    mo_times <- seq(0, min(total_h - 720, 11*720), by = 720)
    ev_list$romo <- ev(cmt = "ROMO_SC", amt = 210, time = mo_times)
  }

  if (length(ev_list) == 0) return(ev(time = 0, amt = 0, cmt = 1))
  Reduce(c, ev_list)
}

## ---- Simulation scenarios ----
run_scenario <- function(
    scenario_name,
    use_aln  = FALSE, use_zol  = FALSE,
    use_dmab = FALSE, use_tptd = FALSE, use_romo = FALSE,
    meno = 1, gio = 0, dur_years = 3
) {
  params <- list(
    USE_ALN  = as.integer(use_aln),
    USE_ZOL  = as.integer(use_zol),
    USE_DMAB = as.integer(use_dmab),
    USE_TPTD = as.integer(use_tptd),
    USE_ROMO = as.integer(use_romo),
    MENO     = meno,
    GIO      = gio
  )
  # Set initial E2 based on menopause status
  if (meno == 1) {
    init_vals <- init(mod, E2 = 12.0, RANKL = 16.0, OPG = 28.0, OC = 1.4, OB = 0.9)
  } else {
    init_vals <- init(mod)
  }

  evs <- build_events(use_aln, use_zol, use_dmab, use_tptd, use_romo, dur_years)

  out <- mod %>%
    param(params) %>%
    init(init_vals) %>%
    ev(evs) %>%
    mrgsim(end = dur_years*365*24, delta = 12) %>%
    as_tibble() %>%
    mutate(
      time_yr = time / (365 * 24),
      scenario = scenario_name
    )
  out
}

## ---- Define 6 scenarios ----
scenarios <- list(
  S1_Untreated = run_scenario("S1: Postmenopausal (untreated)", meno = 1),
  S2_ALN = run_scenario("S2: Alendronate 70mg/wk", use_aln = TRUE, meno = 1),
  S3_ZOL = run_scenario("S3: Zoledronate 5mg/yr", use_zol = TRUE, meno = 1),
  S4_DMAB = run_scenario("S4: Denosumab 60mg/6mo", use_dmab = TRUE, meno = 1),
  S5_TPTD = run_scenario("S5: Teriparatide 20µg/day", use_tptd = TRUE, meno = 1),
  S6_ROMO_DMAB = run_scenario("S6: Romosozumab→Denosumab", use_romo = TRUE, use_dmab = TRUE, meno = 1)
)

results <- bind_rows(scenarios)

## ---- Plots ----
theme_qsp <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom", strip.background = element_blank())

# BMD over time
p_bmd <- ggplot(results, aes(time_yr, BMD_gcm2, color = scenario)) +
  geom_line(size = 0.9) +
  geom_hline(yintercept = 0.85*0.75, linetype = "dashed", color = "red", alpha = 0.5) +
  annotate("text", x = 0.1, y = 0.85*0.74, label = "Osteoporosis threshold", hjust = 0, size = 3) +
  labs(title = "Lumbar Spine BMD Over Time",
       x = "Years", y = "BMD (g/cm²)", color = NULL) +
  scale_color_brewer(palette = "Set1") + theme_qsp

# CTX (resorption marker)
p_ctx <- ggplot(results, aes(time_yr, CTX_ngl, color = scenario)) +
  geom_line(size = 0.9) +
  geom_hline(yintercept = 0.35, linetype = "dashed", alpha = 0.5) +
  labs(title = "CTX (Bone Resorption Marker)",
       x = "Years", y = "CTX (ng/mL)", color = NULL) +
  scale_color_brewer(palette = "Set1") + theme_qsp

# P1NP (formation marker)
p_p1np <- ggplot(results, aes(time_yr, P1NP_ugl, color = scenario)) +
  geom_line(size = 0.9) +
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.5) +
  labs(title = "P1NP (Bone Formation Marker)",
       x = "Years", y = "P1NP (µg/L)", color = NULL) +
  scale_color_brewer(palette = "Set1") + theme_qsp

# OB/OC dynamics
p_cells <- results %>%
  pivot_longer(c(OB_rel, OC_rel), names_to = "cell", values_to = "activity") %>%
  mutate(cell = ifelse(cell == "OB_rel", "Osteoblast", "Osteoclast")) %>%
  ggplot(aes(time_yr, activity, color = scenario, linetype = cell)) +
  geom_line(size = 0.8) +
  facet_wrap(~cell, scales = "free_y") +
  labs(title = "Bone Cell Dynamics",
       x = "Years", y = "Relative Activity", color = NULL, linetype = NULL) +
  scale_color_brewer(palette = "Set1") + theme_qsp

# RANKL/OPG
p_rankl_opg <- results %>%
  filter(scenario %in% c("S1: Postmenopausal (untreated)", "S4: Denosumab 60mg/6mo")) %>%
  pivot_longer(c(RANKL_pm, OPG_pm), names_to = "marker", values_to = "conc") %>%
  ggplot(aes(time_yr, conc, color = scenario, linetype = marker)) +
  geom_line(size = 0.9) +
  labs(title = "RANKL & OPG Dynamics",
       x = "Years", y = "Concentration (pmol/L)", color = NULL, linetype = NULL) +
  theme_qsp

# Fracture risk accumulation
p_fx <- ggplot(results, aes(time_yr, FxRisk10yr*100, color = scenario)) +
  geom_line(size = 0.9) +
  labs(title = "Estimated 10-Year Fracture Risk",
       x = "Years", y = "10-yr Fracture Probability (%)", color = NULL) +
  scale_color_brewer(palette = "Set1") + theme_qsp

## ---- Summary table at 3 years ----
summary_tbl <- results %>%
  filter(abs(time_yr - 3) == min(abs(time_yr - 3)), .by = scenario) %>%
  select(scenario, BMD_gcm2, T_score_val, CTX_ngl, P1NP_ugl,
         OB_rel, OC_rel, FxRisk10yr) %>%
  mutate(across(where(is.numeric), ~round(.x, 3))) %>%
  rename(
    Scenario       = scenario,
    `BMD (g/cm²)`  = BMD_gcm2,
    `T-score`      = T_score_val,
    `CTX (ng/mL)`  = CTX_ngl,
    `P1NP (µg/L)`  = P1NP_ugl,
    `OB (rel.)`    = OB_rel,
    `OC (rel.)`    = OC_rel,
    `Fx Risk 10yr` = FxRisk10yr
  )

print(summary_tbl)

## ---- GIOP scenario ----
gio_scen <- run_scenario("GIOP untreated", meno = 0, gio = 1)
gio_tptd <- run_scenario("GIOP + Teriparatide", use_tptd = TRUE, meno = 0, gio = 1)
gio_zol  <- run_scenario("GIOP + Zoledronate",  use_zol  = TRUE, meno = 0, gio = 1)
gio_results <- bind_rows(gio_scen, gio_tptd, gio_zol)

p_gio <- ggplot(gio_results, aes(time_yr, BMD_gcm2, color = scenario)) +
  geom_line(size = 0.9) +
  labs(title = "Glucocorticoid-Induced Osteoporosis: BMD Response",
       x = "Years", y = "BMD (g/cm²)", color = NULL) +
  scale_color_brewer(palette = "Dark2") + theme_qsp

## ---- Sequential therapy (ROMO → DMAB → ALN) ----
# Romosozumab 12 months, then denosumab 12 months, then alendronate
romo_seq <- run_scenario("ROMO (yr1) → DMAB (yr2–3)",
                         use_romo = TRUE, use_dmab = TRUE, meno = 1, dur_years = 3)

cat("\n=== Osteoporosis QSP Model Summary ===\n")
cat("States: 22 ODE compartments\n")
cat("Drugs: Alendronate, Zoledronate, Denosumab, Teriparatide, Romosozumab\n")
cat("Scenarios: 6 postmenopausal + 3 GIOP\n\n")

gridExtra::grid.arrange(p_bmd, p_ctx, p_p1np, p_cells, nrow = 2)
