## ============================================================
## Familial Hypercholesterolemia (FH) QSP Model
## mrgsolve ODE-based PK/PD model
##
## Disease:   Familial Hypercholesterolemia (LDLR/APOB/PCSK9)
## Drugs:     Statins, PCSK9 inhibitors (evolocumab), Ezetimibe,
##            Inclisiran, Bempedoic acid, Lomitapide
## Scenarios: HetFH baseline, statin mono, statin+EZE,
##            PCSK9i mono, statin+PCSK9i, lomitapide (HomFH)
##
## References: Stecula et al. 2019, Harrington et al. 2023,
##             FOURIER trial 2017, ORION trials 2020-2022,
##             Watts et al. 2020 (FH global guidelines)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ----------------------------------------------------------
## Model code
## ----------------------------------------------------------
fh_code <- '
$PLUGIN autodiff

$PARAM
// === Statin PK (Rosuvastatin-like, oral, mg → μg/mL) ===
ka_S    = 0.45,    // GI absorption rate constant (h^-1)
F_S     = 0.20,    // oral bioavailability (hepatic first-pass ~80%)
Vd_S    = 2.0,     // apparent volume statin (L/kg body weight equiv)
CL_S    = 0.60,    // systemic clearance statin (L/h)
CL_H    = 0.45,    // hepatic clearance (L/h) — predominant route
DOSE_S  = 0.0,     // statin daily dose (mg), 0 = no statin

// === PCSK9 inhibitor PK (Evolocumab-like, subcutaneous mAb) ===
ka_P    = 0.0245,  // SC absorption rate (h^-1) ~ bioavail/t_abs
F_P     = 0.72,    // SC bioavailability evolocumab
Vc_P    = 3.1,     // central volume (L)
Vp_P    = 2.6,     // peripheral volume (L)
CL_P    = 0.0060,  // linear clearance evolocumab (L/h)
Q_P     = 0.0043,  // intercompartment clearance (L/h)
kTMDD   = 0.00015, // target-mediated (PCSK9-bound) clearance (L/h/nmol)
DOSE_P  = 0.0,     // PCSK9i dose (mg), 0 = none  [420 mg q4w]
MW_P    = 144000,  // molecular weight evolocumab (Da)

// === Ezetimibe PK (oral, mg) ===
ka_EZE  = 0.35,    // absorption rate (h^-1)
F_EZE   = 0.35,    // oral bioavailability (active + inactive EH recycling)
Vd_EZE  = 2.5,     // volume of distribution (L/kg)
CL_EZE  = 0.15,    // systemic clearance (L/h)
DOSE_EZE= 0.0,     // ezetimibe dose (mg, typically 10 mg/day)

// === HMGCR/LDLR turnover ===
kout_HMG = 0.28,   // HMGCR turnover rate (h^-1) ~2.5 h half-life
EC50_SI  = 0.08,   // statin hepatic conc for 50% HMGCR inhibition (μg/mL)
Emax_SI  = 0.92,   // maximum fraction HMGCR inhibition by statin
kin_LR   = 0.025,  // LDLR synthesis rate (h^-1, relative to baseline)
kout_LR  = 0.020,  // LDLR degradation baseline (h^-1)
EC50_LD  = 0.30,   // HMGCR inhib fraction for 50% LDLR upregulation
Emax_LD  = 2.5,    // max fold-increase in LDLR synthesis due to SREBP2

// === PCSK9 dynamics ===
kin_PK9  = 10.0,   // PCSK9 synthesis rate (ng/mL/h) → steady-state ~400 ng/mL
kout_PK9 = 0.024,  // PCSK9 clearance (h^-1), t1/2 ~29 h
kon_PK9  = 0.012,  // PCSK9i binding to PCSK9 (1/nM/h) = kon
koff_PK9 = 1e-5,   // PCSK9i dissociation (h^-1), Kd ~1 pM
EC50_PK9_LR = 200, // PCSK9 conc for 50% LDLR degradation (ng/mL)
kP9_LDLR = 0.010,  // PCSK9-mediated extra LDLR degradation (per unit PCSK9)

// === Lipoprotein dynamics (mg/dL) ===
LDL0_het = 280,    // HetFH baseline LDL-C (mg/dL)
LDL0_hom = 680,    // HomFH baseline LDL-C (mg/dL)
VLDL0    = 32,     // baseline VLDL-C (mg/dL)
IDL0     = 15,     // baseline IDL-C (mg/dL)
HDL0     = 46,     // baseline HDL-C (mg/dL)
TG0      = 155,    // baseline TG (mg/dL)
PCSK9_0  = 400,    // baseline plasma PCSK9 (ng/mL)
LDLR_0   = 1.0,    // baseline LDLR (normalized, =1)

// VLDL kinetics
kVLDL_s  = 0.040,  // VLDL → IDL conversion (h^-1)
kIDL_s   = 0.085,  // IDL → LDL conversion (h^-1)
kLDL_cl  = 0.0065, // LDL clearance baseline (h^-1)
kLDL_cl_LR = 0.028,// Additional LDL clearance per LDLR unit (h^-1)
kHDL_eq  = 0.012,  // HDL equilibration rate
fEZE_abs = 0.55,   // fraction reduction in intestinal cholesterol abs by EZE
fLomit   = 0.80,   // max MTP inhibition by lomitapide (→ ↓VLDL)
DOSE_LOMT= 0.0,    // Lomitapide dose flag (1 = on)
DOSE_BEMP= 0.0,    // Bempedoic acid flag (1 = on, adds ~18% LDL reduction)
fBemp    = 0.18,   // max additional LDL-C reduction by bempedoic acid

// === Genetic parameters ===
LDLR_fxn = 0.50,   // LDLR residual function (0.5=HetFH, 0.05=HomFH, 1=normal)

// === Body parameters ===
BW       = 75.0,   // body weight (kg)
Vmix     = 4.0     // plasma mixing volume for LP (L)

$CMT
// PK compartments
GUT_S              // Statin GI (mg)
CENT_S             // Statin central plasma (μg)
LIV_S              // Statin liver (μg) — primary site of action
SC_PCSK9I          // PCSK9i subcutaneous depot (mg)
CENT_PCSK9I        // PCSK9i central plasma (mg)
PERI_PCSK9I        // PCSK9i peripheral compartment (mg)
COMP_PK9           // PCSK9i–PCSK9 complex (nmol/L equiv)
GUT_EZE            // Ezetimibe GI (mg)
CENT_EZE           // Ezetimibe plasma (mg)

// PD compartments
HMGCR_rel          // HMGCR relative activity (1 = normal)
LDLR_rel           // LDLR surface level (1 = baseline for genotype)
PCSK9_pl           // Plasma PCSK9 (ng/mL)

// Lipoprotein compartments (mg/dL)
VLDL_C             // VLDL-C
IDL_C              // IDL-C
LDL_C              // LDL-C
HDL_C              // HDL-C
TG_C               // Triglycerides

$INIT
GUT_S       = 0,
CENT_S      = 0,
LIV_S       = 0,
SC_PCSK9I   = 0,
CENT_PCSK9I = 0,
PERI_PCSK9I = 0,
COMP_PK9    = 0,
GUT_EZE     = 0,
CENT_EZE    = 0,
HMGCR_rel   = 1.0,
LDLR_rel    = 1.0,
PCSK9_pl    = PCSK9_0,
VLDL_C      = VLDL0,
IDL_C       = IDL0,
LDL_C       = LDL0_het,
HDL_C       = HDL0,
TG_C        = TG0

$ODE

// -------------------------------------------------------
// [1-3] STATIN PK
// -------------------------------------------------------
double dose_rate_S = 0;  // oral bolus handled via events

double C_liv_S = LIV_S / (Vd_S * BW * 0.26); // hepatic conc (μg/mL), liver ~26% Vd
double C_sys_S = CENT_S / (Vd_S * BW);         // systemic conc (μg/mL)

dxdt_GUT_S  = -ka_S * GUT_S;
dxdt_CENT_S =  ka_S * F_S * GUT_S - (CL_S / (Vd_S * BW)) * CENT_S
               - (CL_H / (Vd_S * BW * 0.26)) * CENT_S;
dxdt_LIV_S  =  ka_S * F_S * GUT_S * 0.80   // first-pass hepatic extraction ~80%
               - (CL_H / (Vd_S * BW * 0.26)) * LIV_S;

// -------------------------------------------------------
// [4-7] PCSK9 INHIBITOR PK (2-cmpt + TMDD)
// -------------------------------------------------------
double PCSK9i_nM = CENT_PCSK9I / (MW_P * 1e-6 * Vc_P); // approx nM

dxdt_SC_PCSK9I   = -ka_P * SC_PCSK9I;
dxdt_CENT_PCSK9I =  ka_P * F_P * SC_PCSK9I
                    - (CL_P / Vc_P) * CENT_PCSK9I
                    - (Q_P  / Vc_P) * CENT_PCSK9I
                    + (Q_P  / Vp_P) * PERI_PCSK9I
                    - kTMDD * PCSK9_pl * CENT_PCSK9I;
dxdt_PERI_PCSK9I = (Q_P / Vc_P) * CENT_PCSK9I
                   - (Q_P / Vp_P) * PERI_PCSK9I;
dxdt_COMP_PK9    =  kon_PK9 * PCSK9i_nM * (PCSK9_pl / 1000.0)
                   - koff_PK9 * COMP_PK9
                   - kout_PK9 * COMP_PK9;  // complex is cleared

// -------------------------------------------------------
// [8-9] EZETIMIBE PK
// -------------------------------------------------------
double C_EZE = CENT_EZE / (Vd_EZE * BW);  // mg/L

dxdt_GUT_EZE  = -ka_EZE * GUT_EZE;
dxdt_CENT_EZE =  ka_EZE * F_EZE * GUT_EZE
                 - (CL_EZE / (Vd_EZE * BW)) * CENT_EZE;

// -------------------------------------------------------
// [10] HMGCR ACTIVITY (inhibited by hepatic statin)
// -------------------------------------------------------
double Inh_S   = Emax_SI * C_liv_S / (EC50_SI + C_liv_S);   // 0→1
double kin_HMG = kout_HMG;   // at steady state HMGCR_rel = 1
double Bemp_eff = DOSE_BEMP * fBemp * 0.55;  // bempedoic acid ACLY→HMGCR effect

dxdt_HMGCR_rel = kin_HMG * (1.0 - Inh_S - Bemp_eff) - kout_HMG * HMGCR_rel;

// -------------------------------------------------------
// [11] LDLR SURFACE EXPRESSION
// -------------------------------------------------------
// SREBP2 activated when HMGCR is inhibited (↓intracellular cholesterol)
double HMGCR_inh_frac = 1.0 - HMGCR_rel;   // fraction inhibited 0→1
double LDLR_up = Emax_LD * HMGCR_inh_frac / (EC50_LD + HMGCR_inh_frac); // fold-increase

// PCSK9 degrades LDLR (concentration-dependent)
double PCSK9_free = PCSK9_pl - COMP_PK9 * 1000.0; // free PCSK9 (ng/mL)
if(PCSK9_free < 0) PCSK9_free = 0;
double PCSK9_effect_LDLR = kP9_LDLR * PCSK9_free / (EC50_PK9_LR + PCSK9_free); // extra degradation

// Genotype scaling
double LDLR_max = LDLR_fxn; // peak possible = genotype-scaled
double kin_LR_eff = kin_LR * LDLR_max * (1.0 + LDLR_up);
double kout_LR_eff = kout_LR + PCSK9_effect_LDLR;

dxdt_LDLR_rel = kin_LR_eff - kout_LR_eff * LDLR_rel;

// -------------------------------------------------------
// [12] PLASMA PCSK9
// -------------------------------------------------------
// Statin increases PCSK9 transcription (SREBP2 feedback)
double PCSK9_statin_feedback = 1.0 + 0.40 * HMGCR_inh_frac; // statins ↑PCSK9 ~40%
double PCSK9_inclisiran_red  = 1.0;  // flag for inclisiran (simplified)

// PCSK9i binds PCSK9
double PCSK9_pl_use = PCSK9_pl;
if(PCSK9_pl_use < 0) PCSK9_pl_use = 0;

dxdt_PCSK9_pl = kin_PK9 * PCSK9_statin_feedback
                - kout_PK9 * PCSK9_pl_use
                - kon_PK9 * PCSK9i_nM * PCSK9_pl_use
                + koff_PK9 * COMP_PK9 * 1000.0;

// -------------------------------------------------------
// [13-17] LIPOPROTEIN DYNAMICS
// -------------------------------------------------------
// Ezetimibe reduces intestinal cholesterol input → ↓VLDL substrate
double EZE_effect_VLDL = 1.0 - fEZE_abs * C_EZE / (0.05 + C_EZE);
// Lomitapide reduces VLDL assembly
double LOMT_eff = 1.0 - DOSE_LOMT * fLomit;
// VLDL production rate
double k_VLDL_prod = kVLDL_s * VLDL0 * EZE_effect_VLDL * LOMT_eff;

dxdt_VLDL_C = k_VLDL_prod
              - kVLDL_s * VLDL_C   // VLDL → IDL lipolysis
              - (HMGCR_inh_frac * 0.20) * VLDL_C; // statins modestly ↓VLDL

// IDL
dxdt_IDL_C  = kVLDL_s * VLDL_C
              - kIDL_s  * IDL_C;   // IDL → LDL

// LDL clearance: baseline + LDLR-mediated + PCSK9i enhanced
double kLDL_total = kLDL_cl + kLDL_cl_LR * LDLR_rel;
// Bempedoic acid additional LDL-C reduction through HMGCR
double bemp_LDL_add = DOSE_BEMP * fBemp;

dxdt_LDL_C  = kIDL_s * IDL_C
              - kLDL_total * LDL_C
              - bemp_LDL_add * kLDL_cl * LDL_C;

// HDL: statins ↑HDL ~6%, CETP inhibition not modelled separately
double HDL_target = HDL0 * (1.0 + 0.06 * HMGCR_inh_frac);
dxdt_HDL_C  = kHDL_eq * (HDL_target - HDL_C);

// TG: statins reduce TG, fibrates not modelled
double TG_target = TG0 * (1.0 - 0.15 * HMGCR_inh_frac);
dxdt_TG_C   = kHDL_eq * (TG_target - TG_C);

$TABLE
// Derived PK
double C_statin_sys  = CENT_S  / (Vd_S  * BW);   // μg/mL plasma statin
double C_statin_liv  = LIV_S   / (Vd_S * BW * 0.26); // μg/mL liver statin
double C_pcsk9i_mgL  = CENT_PCSK9I / Vc_P;        // mg/L evolocumab central
double C_eze_mgL     = CENT_EZE / (Vd_EZE * BW);  // mg/L ezetimibe

// Derived PD
double LDL_reduction_pct  = (LDL0_het > 0) ? (LDL0_het - LDL_C) / LDL0_het * 100.0 : 0;
double NonHDL_C = LDL_C + VLDL_C + IDL_C;
double TC       = LDL_C + VLDL_C + HDL_C + IDL_C;
double LDLR_pct = LDLR_rel * 100.0;

// CVD risk surrogates (simplified Framingham-style)
double CVD_risk_10yr = 0.01 * exp(0.0035 * LDL_C - 0.01 * HDL_C);

// LDL goals (ESC 2019 very-high risk: <55 mg/dL; high risk: <70 mg/dL)
double LDL_goal_55  = (LDL_C <= 55) ? 1 : 0;
double LDL_goal_70  = (LDL_C <= 70) ? 1 : 0;
double PCSK9_free   = (PCSK9_pl - COMP_PK9 * 1000.0 < 0) ? 0 : PCSK9_pl - COMP_PK9 * 1000.0;

$CAPTURE
C_statin_sys, C_statin_liv, C_pcsk9i_mgL, C_eze_mgL,
HMGCR_rel, LDLR_rel, LDLR_pct,
PCSK9_pl, PCSK9_free, COMP_PK9,
VLDL_C, IDL_C, LDL_C, HDL_C, TG_C,
NonHDL_C, TC, LDL_reduction_pct, CVD_risk_10yr,
LDL_goal_55, LDL_goal_70
'

## ----------------------------------------------------------
## Compile model
## ----------------------------------------------------------
mod <- mcode("FH_QSP", fh_code)
cat("Model compiled successfully. Compartments:", length(Init(mod)), "\n")

## ----------------------------------------------------------
## Helper: build event data for daily oral + q4w SC dosing
## ----------------------------------------------------------
build_ev <- function(dose_statin_mg  = 0,
                     dose_pcsk9i_mg  = 0,    # 0 or 420 mg q4w
                     dose_eze_mg     = 0,
                     n_days          = 365) {
  evlist <- list()

  # Statin: once daily oral (add to GUT_S)
  if (dose_statin_mg > 0) {
    ev_s <- ev(amt = dose_statin_mg, cmt = "GUT_S", ii = 24, addl = n_days - 1)
    evlist[["statin"]] <- ev_s
  }

  # PCSK9i: q4w subcutaneous (every 28 days)
  if (dose_pcsk9i_mg > 0) {
    n_inj <- floor(n_days / 28)
    ev_p <- ev(amt = dose_pcsk9i_mg, cmt = "SC_PCSK9I", ii = 28 * 24, addl = n_inj - 1)
    evlist[["pcsk9i"]] <- ev_p
  }

  # Ezetimibe: once daily oral
  if (dose_eze_mg > 0) {
    ev_e <- ev(amt = dose_eze_mg, cmt = "GUT_EZE", ii = 24, addl = n_days - 1)
    evlist[["eze"]] <- ev_e
  }

  if (length(evlist) == 0) return(NULL)
  Reduce(c, evlist)
}

## ----------------------------------------------------------
## Treatment Scenarios
## ----------------------------------------------------------
scenarios <- list(
  list(
    label    = "1. HetFH — No Treatment",
    LDLR_fxn = 0.50,
    LDL0     = 280,
    ev_      = NULL,
    DOSE_P   = 0, DOSE_S = 0, DOSE_EZE = 0,
    DOSE_LOMT= 0, DOSE_BEMP = 0,
    color    = "#E74C3C"
  ),
  list(
    label    = "2. HetFH — Rosuvastatin 40 mg/d",
    LDLR_fxn = 0.50,
    LDL0     = 280,
    ev_      = build_ev(dose_statin_mg = 40),
    DOSE_P   = 0, DOSE_S = 40, DOSE_EZE = 0,
    DOSE_LOMT= 0, DOSE_BEMP = 0,
    color    = "#E67E22"
  ),
  list(
    label    = "3. HetFH — Rosuvastatin 40 mg + Ezetimibe 10 mg",
    LDLR_fxn = 0.50,
    LDL0     = 280,
    ev_      = build_ev(dose_statin_mg = 40, dose_eze_mg = 10),
    DOSE_P   = 0, DOSE_S = 40, DOSE_EZE = 10,
    DOSE_LOMT= 0, DOSE_BEMP = 0,
    color    = "#F39C12"
  ),
  list(
    label    = "4. HetFH — Evolocumab 420 mg q4w (PCSK9i mono)",
    LDLR_fxn = 0.50,
    LDL0     = 280,
    ev_      = build_ev(dose_pcsk9i_mg = 420),
    DOSE_P   = 420, DOSE_S = 0, DOSE_EZE = 0,
    DOSE_LOMT= 0, DOSE_BEMP = 0,
    color    = "#27AE60"
  ),
  list(
    label    = "5. HetFH — Rosuvastatin 40 mg + Evolocumab 420 mg q4w",
    LDLR_fxn = 0.50,
    LDL0     = 280,
    ev_      = build_ev(dose_statin_mg = 40, dose_pcsk9i_mg = 420),
    DOSE_P   = 420, DOSE_S = 40, DOSE_EZE = 0,
    DOSE_LOMT= 0, DOSE_BEMP = 0,
    color    = "#2E86C1"
  ),
  list(
    label    = "6. HomFH — Lomitapide + Statin + Evolocumab",
    LDLR_fxn = 0.05,
    LDL0     = 680,
    ev_      = build_ev(dose_statin_mg = 40, dose_pcsk9i_mg = 420),
    DOSE_P   = 420, DOSE_S = 40, DOSE_EZE = 0,
    DOSE_LOMT= 1, DOSE_BEMP = 0,
    color    = "#8E44AD"
  )
)

## ----------------------------------------------------------
## Run simulations
## ----------------------------------------------------------
sim_results <- lapply(scenarios, function(sc) {
  # Update initial LDL for genotype
  init_vals <- Init(mod)
  init_vals["LDL_C"]    <- sc$LDL0
  init_vals["LDLR_rel"] <- sc$LDLR_fxn
  init_vals["PCSK9_pl"] <- 400

  params_update <- list(
    LDLR_fxn  = sc$LDLR_fxn,
    LDL0_het  = sc$LDL0,
    DOSE_S    = sc$DOSE_S,
    DOSE_P    = sc$DOSE_P,
    DOSE_EZE  = sc$DOSE_EZE,
    DOSE_LOMT = sc$DOSE_LOMT,
    DOSE_BEMP = sc$DOSE_BEMP
  )

  m2 <- mod %>%
    param(params_update) %>%
    init(init_vals)

  out <- if (!is.null(sc$ev_)) {
    mrgsim(m2, sc$ev_, end = 365 * 24, delta = 12, carry_out = "evid")
  } else {
    mrgsim(m2, end = 365 * 24, delta = 12)
  }

  as.data.frame(out) %>%
    mutate(Scenario = sc$label, Color = sc$color,
           time_days = time / 24)
})

results_df <- bind_rows(sim_results)

## ----------------------------------------------------------
## Plot 1: LDL-C over 52 weeks
## ----------------------------------------------------------
p1 <- ggplot(results_df, aes(time_days, LDL_C, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 70, linetype = "dashed", color = "navy", alpha = 0.7) +
  geom_hline(yintercept = 55, linetype = "dotted", color = "darkred", alpha = 0.7) +
  annotate("text", x = 350, y = 75,  label = "70 mg/dL (ESC high risk)",
           color = "navy", size = 3, hjust = 1) +
  annotate("text", x = 350, y = 60,  label = "55 mg/dL (ESC very high risk)",
           color = "darkred", size = 3, hjust = 1) +
  scale_color_manual(values = setNames(sapply(scenarios, `[[`, "color"),
                                       sapply(scenarios, `[[`, "label"))) +
  labs(title = "Familial Hypercholesterolemia — LDL-C Time Course",
       x = "Time (days)", y = "LDL-C (mg/dL)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", legend.text = element_text(size = 8)) +
  guides(color = guide_legend(ncol = 2))

## ----------------------------------------------------------
## Plot 2: LDLR Surface Expression
## ----------------------------------------------------------
p2 <- ggplot(results_df, aes(time_days, LDLR_pct, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = setNames(sapply(scenarios, `[[`, "color"),
                                       sapply(scenarios, `[[`, "label"))) +
  labs(title = "Hepatic LDLR Surface Expression",
       x = "Time (days)", y = "LDLR (%  of normal)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

## ----------------------------------------------------------
## Plot 3: Plasma PCSK9 levels
## ----------------------------------------------------------
p3 <- ggplot(results_df, aes(time_days, PCSK9_pl, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = setNames(sapply(scenarios, `[[`, "color"),
                                       sapply(scenarios, `[[`, "label"))) +
  labs(title = "Plasma PCSK9 Concentration",
       x = "Time (days)", y = "PCSK9 (ng/mL)",
       color = "Treatment Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

## ----------------------------------------------------------
## Plot 4: Lipid Panel Summary at Week 52
## ----------------------------------------------------------
summary_52w <- results_df %>%
  filter(abs(time_days - 364) < 1) %>%
  group_by(Scenario) %>%
  summarise(
    LDL_C    = mean(LDL_C),
    HDL_C    = mean(HDL_C),
    TG       = mean(TG_C),
    NonHDL_C = mean(NonHDL_C),
    .groups  = "drop"
  ) %>%
  pivot_longer(cols = c(LDL_C, HDL_C, TG, NonHDL_C),
               names_to = "Lipid", values_to = "mgdL")

p4 <- ggplot(summary_52w, aes(Scenario, mgdL, fill = Lipid)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Lipid Panel at Week 52",
       x = "", y = "Concentration (mg/dL)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

## ----------------------------------------------------------
## Plot 5: % LDL Reduction waterfall
## ----------------------------------------------------------
ldl_reduc <- results_df %>%
  filter(abs(time_days - 364) < 1) %>%
  group_by(Scenario, Color) %>%
  summarise(ldl_pct = mean(LDL_reduction_pct), .groups = "drop") %>%
  arrange(desc(ldl_pct))

p5 <- ggplot(ldl_reduc, aes(reorder(Scenario, ldl_pct), ldl_pct, fill = Scenario)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", ldl_pct)), hjust = -0.1, size = 4) +
  scale_fill_manual(values = setNames(ldl_reduc$Color, ldl_reduc$Scenario)) +
  coord_flip(ylim = c(-50, 90)) +
  labs(title = "LDL-C Reduction at Week 52\n(vs. HetFH untreated baseline 280 mg/dL)",
       x = "", y = "% LDL-C Reduction") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

## ----------------------------------------------------------
## Plot 6: HMGCR Activity over time
## ----------------------------------------------------------
p6 <- ggplot(results_df, aes(time_days, HMGCR_rel * 100, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = setNames(sapply(scenarios, `[[`, "color"),
                                       sapply(scenarios, `[[`, "label"))) +
  labs(title = "HMGCR Relative Activity",
       x = "Time (days)", y = "HMGCR Activity (% baseline)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

## ----------------------------------------------------------
## Combine & save
## ----------------------------------------------------------
combined_plot <- (p1 + p5) / (p2 + p3) / (p4 + p6) +
  plot_annotation(
    title    = "Familial Hypercholesterolemia (FH) — QSP Simulation Dashboard",
    subtitle = "mrgsolve ODE model: PK/PD of statin, PCSK9 inhibitor, ezetimibe, lomitapide",
    theme    = theme(plot.title = element_text(size = 16, face = "bold"),
                     plot.subtitle = element_text(size = 12))
  )

ggsave("fh_simulation_dashboard.png", combined_plot,
       width = 18, height = 16, dpi = 150)
cat("Dashboard saved: fh_simulation_dashboard.png\n")

## ----------------------------------------------------------
## Steady-state sensitivity: LDL reduction vs statin dose
## ----------------------------------------------------------
statin_doses <- c(5, 10, 20, 40, 80)

ss_results <- lapply(statin_doses, function(d) {
  m2 <- mod %>%
    param(LDLR_fxn = 0.50, LDL0_het = 280, DOSE_S = d) %>%
    init(LDL_C = 280, LDLR_rel = 0.50)
  ev_ <- build_ev(dose_statin_mg = d)
  out <- mrgsim(m2, ev_, end = 180 * 24, delta = 24)
  tail(as.data.frame(out), 1) %>%
    mutate(Dose = d, LDL_reduction = LDL_reduction_pct)
})
ss_df <- bind_rows(ss_results)

cat("\nSteady-state LDL-C reduction by statin dose (HetFH):\n")
cat(sprintf("  Dose %3d mg/d → LDL-C = %5.1f mg/dL (↓%.0f%%)\n",
            ss_df$Dose, ss_df$LDL_C, ss_df$LDL_reduction))

## ----------------------------------------------------------
## Print summary table
## ----------------------------------------------------------
cat("\n=== LDL-C Summary at Week 52 ===\n")
summ <- results_df %>%
  filter(abs(time_days - 364) < 1) %>%
  group_by(Scenario) %>%
  summarise(
    LDL_C_mgdL = round(mean(LDL_C), 1),
    HDL_C_mgdL = round(mean(HDL_C), 1),
    TG_mgdL    = round(mean(TG_C), 1),
    PCSK9_ngmL = round(mean(PCSK9_pl), 0),
    LDLR_pct   = round(mean(LDLR_pct), 0),
    LDL_red_pct= round(mean(LDL_reduction_pct), 1),
    LDL_goal55 = round(mean(LDL_goal_55) * 100),
    LDL_goal70 = round(mean(LDL_goal_70) * 100),
    .groups    = "drop"
  )
print(as.data.frame(summ))
