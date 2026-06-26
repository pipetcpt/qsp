##############################################################################
# Hashimoto's Thyroiditis — Quantitative Systems Pharmacology Model
# mrgsolve ODE-based PK/PD Model
#
# Model states (18 ODEs):
#   HPT axis    : TRH, TSH
#   Thyroid     : T4thy, T3thy (synthesis compartments)
#   Plasma/tissue: T4p, T3p, T4t, T3t
#   Autoimmune  : Th1, Treg, Bcell, AntiTPOAb, AntiTgAb
#   Disease     : DmgThy (thyroid damage 0-1)
#   Drug PK     : LT4g (gut), LT4c (central), LT4p (peripheral)
#   Selenium    : Se
#
# Treatment scenarios:
#   1. Untreated Hashimoto's (natural progression)
#   2. Standard levothyroxine (LT4 1.6 μg/kg/day)
#   3. Selenium supplementation (200 μg/day)
#   4. Combination LT4 + selenium
#   5. LT4 + liothyronine (T3) combination
#   6. High-dose LT4 (TSH suppressive)
#   7. Early intervention + selenium (pre-overt hypothyroid)
#
# Key references:
#   - Bianco AC et al. Endocr Rev 2019 (Deiodinase kinetics)
#   - Jonklaas J et al. Thyroid 2014 (LT4 PK parameters)
#   - Gärtner R et al. J Clin Endocrinol Metab 2002 (Selenium RCT)
#   - Eisenberg M et al. Thyroid 2008 (TSH-T4 model)
#   - Pandiyan B et al. Thyroid 2011 (Mathematical HPT model)
##############################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ---- Model definition -------------------------------------------------------
ht_model_code <- '
$PROB
Hashimoto Thyroiditis QSP Model
18-state ODE: HPT axis, thyroid synthesis, peripheral metabolism,
autoimmune cell populations, thyroid damage, LT4 PK, selenium

$PARAM
// ---- HPT Axis Parameters ----
k_TRH_prod  = 0.5      // TRH synthesis rate (nmol/L/day)
k_TRH_deg   = 1.0      // TRH degradation rate (1/day)
k_TSH_stim  = 4.0      // TRH-stimulated TSH production (mIU/L/day per normalized TRH)
k_TSH_deg   = 6.0      // TSH degradation (1/day, t1/2~2.8h)
EC50_T3_TSH = 0.005    // T3 conc. for 50% TSH suppression (nmol/L relative scale)
hill_TSH    = 2.5      // Hill coefficient, TSH-T3 feedback
TSH0        = 2.0      // Baseline TSH (mIU/L)
TRH0        = 0.5      // Baseline TRH (normalized)

// ---- Thyroid Synthesis ----
k_T4_syn    = 0.12     // Basal T4 synthesis rate (nmol/L/day per TSH unit)
k_T3_syn    = 0.008    // Basal T3 synthesis rate (nmol/L/day per TSH unit)
k_T4sec     = 1.5      // T4 secretion rate from gland (1/day)
k_T3sec     = 1.5      // T3 secretion from gland (1/day)
EC50_TSH_T4 = 2.0      // TSH for half-max T4 stimulation (mIU/L)
hill_T4     = 1.5      // Hill coefficient for TSH-T4 stimulation
T4thy0      = 0.12     // Baseline thyroidal T4 (nmol/L equiv.)
T3thy0      = 0.008

// ---- Peripheral T4/T3 Metabolism ----
k_T4_deg    = 0.099    // T4 plasma elimination (1/day, t1/2~7 days)
k_T3_deg    = 0.693    // T3 plasma elimination (1/day, t1/2~1 day)
k_D1        = 0.3      // DIO1-mediated T4→T3 conversion (1/day)
k_D2        = 0.15     // DIO2-mediated T4→T3 (1/day)
k_T4_tissue = 0.2      // T4 tissue uptake (1/day)
k_T4_tissue_ret = 0.1  // T4 return from tissue
k_T3_tissue = 0.35     // T3 tissue uptake
k_T3_tissue_ret = 0.15 // T3 return
T4p0        = 0.12     // Steady-state plasma free T4 (nmol/L ~ 12 pmol/L)
T3p0        = 0.005    // Steady-state plasma free T3 (nmol/L ~ 5 pmol/L)

// ---- Autoimmune Dynamics ----
k_Th1_base  = 0.01     // Basal Th1 activation rate
k_Th1_stim  = 2.5      // Th1 stimulation by antigen/APC (1/day)
k_Th1_treg  = 3.0      // Treg suppression of Th1 (1/day per Treg unit)
k_Th1_deg   = 0.5      // Th1 turnover (1/day)
k_Treg_prod = 0.3      // Treg homeostatic production
k_Treg_deg  = 0.4      // Treg turnover
k_Treg_inh  = 0.8      // Th1-mediated Treg suppression
Treg_base   = 0.75     // Baseline Treg level (relative, Hashimoto reduced)
k_Bcell_stim= 1.5      // B cell activation by Th1
k_Bcell_deg = 0.6      // B cell turnover
k_Ab_prod   = 0.3      // Anti-TPO Ab production by Bcell
k_Ab_prod2  = 0.2      // Anti-Tg Ab production
k_Ab_deg    = 0.05     // Antibody clearance (t1/2~14 days)
Th1_0       = 0.15     // Baseline Th1 (elevated in Hashimoto)
Treg_0      = 0.65
Bcell_0     = 0.20
Ab_TPO_0    = 150.0    // Baseline anti-TPO Ab (IU/mL, Hashimoto patient)
Ab_Tg_0     = 120.0    // Baseline anti-Tg Ab

// ---- Thyroid Damage ----
k_dmg       = 0.8      // Damage rate (per unit Th1 x AbTPO)
k_repair    = 0.2      // Repair rate by Treg
DmgThy_0    = 0.25     // Initial damage (moderate, established Hashimoto)
k_dmg_scale = 0.001    // Scale factor

// ---- Levothyroxine PK (2-compartment) ----
// Per-dose parameters; dose given as μg/day (will be converted to nmol)
MW_T4       = 776.87   // Molecular weight LT4 (g/mol)
F_LT4       = 0.75     // Bioavailability (~75%)
ka_LT4      = 0.48     // Absorption rate (1/h = 11.5/day)
CL_LT4      = 1.3      // Clearance (L/day)
Vc_LT4      = 10.0     // Central volume (L)
Vp_LT4      = 22.0     // Peripheral volume (L)
k12_LT4     = 0.25     // Distribution (1/day)
k21_LT4     = 0.11     // Return (1/day)
kel_LT4     = 0.13     // Elimination from central (1/day, t1/2~5-7 days)
f_T4conv    = 0.35     // Fraction LT4 → T3 via deiodinase

// ---- Liothyronine PK ----
ka_LiT3     = 5.0      // T3 rapid absorption (1/day)
CL_LiT3     = 25.0     // T3 clearance (L/day)
Vc_LiT3     = 40.0     // T3 volume of distribution
kel_LiT3    = 0.625    // T3 elimination (t1/2~1.1 days)

// ---- Selenium PK/PD ----
F_Se        = 0.85     // Selenium bioavailability
k_Se_abs    = 2.0      // Absorption rate (1/day)
k_Se_elim   = 0.15     // Selenium elimination (1/day)
Vd_Se       = 50.0     // Se distribution volume (L)
Se_0        = 80.0     // Baseline plasma Se (μg/L, normal ~90-120, low in HT)
Emax_Se_GPx = 0.6      // Max GPx enhancement by Se
EC50_Se_GPx = 90.0     // Se conc. for 50% GPx effect (μg/L)
Emax_Se_Ab  = 0.55     // Max anti-TPO Ab reduction by Se
EC50_Se_Ab  = 95.0     // Se for 50% Ab reduction
Se_dose     = 0.0      // Selenium dose (μg/day, 0=off, 200=therapeutic)

// ---- Drug doses (μg/day) ----
LT4_dose    = 0.0      // Levothyroxine dose (μg/day, typical 75-150 μg)
LiT3_dose   = 0.0      // Liothyronine dose (μg/day, typical 5-10 μg)
MMI_block   = 0.0      // Methimazole effect (0=none, 0-1 scale)

$INIT
// HPT Axis
TRH    = 0.5       // normalized TRH
TSH    = 2.0       // mIU/L

// Thyroid compartments
T4thy  = 0.12      // nmol/L equiv.
T3thy  = 0.008

// Plasma hormones
T4p    = 0.10      // Free T4 (nmol/L); slightly low (Hashimoto)
T3p    = 0.0045    // Free T3 (nmol/L); slightly low

// Tissue
T4t    = 0.15
T3t    = 0.008

// Autoimmune
Th1    = 0.15
Treg   = 0.65
Bcell  = 0.20
AntiTPOAb = 150.0  // IU/mL
AntiTgAb  = 120.0  // IU/mL

// Thyroid damage (0-1 scale)
DmgThy = 0.25

// LT4 PK
LT4g   = 0.0
LT4c   = 0.0
LT4p   = 0.0

// Liothyronine
LiT3c  = 0.0

// Selenium
Se     = 80.0      // μg/L

$OMEGA
0.09   // ETA1: BSV on LT4 CL
0.09   // ETA2: BSV on Thyroid damage
0.09   // ETA3: BSV on Anti-TPO baseline

$SIGMA
0.04   // Proportional residual error T4
0.04   // Proportional residual error TSH
0.09   // Proportional residual error Ab

$MAIN
// Individual parameters with BSV
double CLi    = CL_LT4  * exp(ETA(1));
double DmgSS  = DmgThy_0 * exp(ETA(2));
double AbSS   = Ab_TPO_0 * exp(ETA(3));

// Selenium effect on GPx (normalizes ROS)
double Se_eff_GPx = Emax_Se_GPx * Se / (EC50_Se_GPx + Se);

// Selenium effect on Anti-TPO reduction
double Se_eff_Ab  = 1.0 - Emax_Se_Ab * Se / (EC50_Se_Ab + Se);

// TSH-stimulated thyroid production (sigmoidal Emax)
double TSH_stim   = pow(TSH, hill_T4) / (pow(EC50_TSH_T4, hill_T4) + pow(TSH, hill_T4));

// Thyrocyte integrity (inversely related to damage)
double thyroid_func = 1.0 - DmgThy;

// Methimazole blocking of TPO
double TPO_activity = 1.0 - MMI_block;

// LT4 contribution to plasma T4 (convert LT4c from μg to nmol/L)
// LT4c is in μg, divide by MW_T4*Vc to get nmol/L equivalent
double LT4c_nmol = LT4c / (MW_T4 * 1e-3 * Vc_LT4);  // nmol/L

// LiT3 contribution to plasma T3
double LiT3c_nmol = LiT3c / (65.5 * 1e-3 * Vc_LiT3);  // T3 MW=651 Da

$ODE
// ================================================================
// HPT AXIS
// ================================================================
// TRH: produced by hypothalamus, inhibited by T3/T4 feedback
double T3_fb = T3p + 0.5 * T4p * k_D2;  // effective T3 at hypothalamus
double TRH_inhib = 1.0 / (1.0 + pow(T3_fb / EC50_T3_TSH, hill_TSH));
dxdt_TRH = k_TRH_prod * TRH_inhib - k_TRH_deg * TRH;

// TSH: stimulated by TRH, inhibited by pituitary T3 (via DIO2)
double T3_pit = T3p + k_D2 * T4p * 2.0;  // pituitary T3 enriched by DIO2
double TSH_inhib = pow(EC50_T3_TSH, hill_TSH) /
                   (pow(EC50_T3_TSH, hill_TSH) + pow(T3_pit * 0.5, hill_TSH));
dxdt_TSH = k_TSH_stim * TRH * TSH_inhib - k_TSH_deg * TSH;

// ================================================================
// THYROID GLAND
// ================================================================
// T4 in thyroid: synthesis (TSH-driven, TPO-dependent) - secretion
double T4_synthesis = k_T4_syn * TSH_stim * thyroid_func * TPO_activity;
dxdt_T4thy = T4_synthesis - k_T4sec * T4thy;

// T3 in thyroid: direct synthesis + T4 conversion
double T3_synthesis = k_T3_syn * TSH_stim * thyroid_func * TPO_activity;
dxdt_T3thy = T3_synthesis - k_T3sec * T3thy;

// ================================================================
// PLASMA T4/T3 (free hormone compartments)
// ================================================================
// T4 plasma: input from thyroid + LT4, output via DIO1, DIO2, DIO3, elimination
// Se effect on deiodinase activity
double DIO_Se_eff = 1.0 + Se_eff_GPx * 0.3;  // Se boosts DIO via selenoprotein

dxdt_T4p = k_T4sec * T4thy                         // from thyroid
           + LT4c_nmol * kel_LT4                    // from LT4 PK (approximate)
           - k_D1 * DIO_Se_eff * T4p                // DIO1 conversion
           - k_D2 * DIO_Se_eff * T4p                // DIO2 conversion
           - k_T4_deg * T4p                          // degradation
           - k_T4_tissue * T4p                       // tissue uptake
           + k_T4_tissue_ret * T4t;                  // return

// T3 plasma: from thyroid + DIO conversion + LiT3
dxdt_T3p = k_T3sec * T3thy                          // thyroidal T3
           + (k_D1 + k_D2) * DIO_Se_eff * T4p       // T4→T3 conversion
           + LiT3c_nmol * kel_LiT3                   // from LiT3 PK
           - k_T3_deg * T3p
           - k_T3_tissue * T3p
           + k_T3_tissue_ret * T3t;

// ================================================================
// TISSUE COMPARTMENTS
// ================================================================
dxdt_T4t = k_T4_tissue * T4p - k_T4_tissue_ret * T4t - k_D2 * T4t * 0.5;
dxdt_T3t = k_T3_tissue * T3p + k_D2 * T4t * 0.5
           - k_T3_tissue_ret * T3t - k_T3_deg * T3t * 0.5;

// ================================================================
// AUTOIMMUNE COMPARTMENTS
// ================================================================
// Antigenic drive proportional to thyroid damage (damaged thyrocytes release Ag)
double Ag_drive = 0.1 + DmgThy * 0.9;

// Th1: activated by Ag, suppressed by Treg
dxdt_Th1 = k_Th1_base + k_Th1_stim * Ag_drive * (1.0 - Th1)
           - k_Th1_treg * Treg * Th1
           - k_Th1_deg * Th1;

// Treg: homeostatic production, suppressed by Th1 (Hashimoto: reduced FOXP3)
dxdt_Treg = k_Treg_prod * Treg_base
            - k_Treg_inh * Th1 * Treg
            - k_Treg_deg * Treg;

// B cells: activated by Th1/Tfh signals
dxdt_Bcell = k_Bcell_stim * Th1 * (1.0 - Bcell)
             - k_Bcell_deg * Bcell;

// Anti-TPO Ab: produced by B/plasma cells, reduced by selenium
dxdt_AntiTPOAb = k_Ab_prod * Bcell * 1000.0        // production (IU/mL scale)
                 * Se_eff_Ab                         // selenium reduces
                 - k_Ab_deg * AntiTPOAb;

// Anti-Tg Ab
dxdt_AntiTgAb  = k_Ab_prod2 * Bcell * 800.0
                 * Se_eff_Ab
                 - k_Ab_deg * AntiTgAb;

// ================================================================
// THYROID DAMAGE (0-1 scale)
// ================================================================
// Damage driven by CTL/Ab/complement; repair by Treg
double dmg_drive = k_dmg * Th1 * (AntiTPOAb / 500.0) * (1.0 + Se_eff_GPx * (-0.4));
dxdt_DmgThy = dmg_drive * k_dmg_scale * (1.0 - DmgThy)
              - k_repair * Treg * DmgThy;

// ================================================================
// LEVOTHYROXINE PK (2-compartment, oral)
// ================================================================
// Dose enters gut compartment (μg/day, continuous infusion model)
double LT4_input = LT4_dose * F_LT4;  // effective daily LT4 absorbed (μg)

// ka per hour (11.5/day), here using daily rates
dxdt_LT4g = LT4_input - ka_LT4 * LT4g;  // gut absorption (day units)
dxdt_LT4c = ka_LT4 * LT4g
             - (kel_LT4 + k12_LT4) * LT4c
             + k21_LT4 * LT4p;
dxdt_LT4p = k12_LT4 * LT4c - k21_LT4 * LT4p;

// ================================================================
// LIOTHYRONINE PK (1-compartment, rapid)
// ================================================================
double LiT3_input = LiT3_dose * 0.95;  // ~95% bioavailability
dxdt_LiT3c = LiT3_input - (kel_LiT3 + k12_LT4 * 0.5) * LiT3c;

// ================================================================
// SELENIUM KINETICS
// ================================================================
double Se_input = Se_dose * F_Se;
dxdt_Se = Se_input / Vd_Se - k_Se_elim * Se;

$TABLE
// ================================================================
// DERIVED OUTPUTS & CLINICAL BIOMARKERS
// ================================================================
// Convert relative T4p to clinical units (nmol/L → pmol/L for fT4)
double fT4_pmolL = T4p * 1000.0;        // pmol/L (normal ~12 pmol/L)
double fT3_pmolL = T3p * 1000.0;        // pmol/L (normal ~5 pmol/L)

// Thyroid function status
double TSH_clinical = TSH;              // mIU/L

// Thyroid volume (relative, inversely related to damage long-term)
// Initially volume may increase (goiter), then decrease
double ThyVol_rel = 1.0 + 0.3 * (0.5 - DmgThy) * 2.0 + 0.1 * (AntiTPOAb / 300.0);

// Clinical hypothyroid score (0-10 scale, based on T3t deficiency)
double T3_norm = T3t / 0.008;  // relative to baseline
double HypoScore = 10.0 * (1.0 - fmin(T3_norm, 1.0));

// Cardiovascular: heart rate effect (T3 on cardiac)
double HR_effect = 60.0 + 20.0 * fmin(T3_norm, 1.5);  // bpm

// LDL-C surrogate (rises with hypothyroid; LDLr expression ↓)
double LDL_rel = 4.0 / fmax(T3_norm, 0.3);  // mmol/L equiv.

// Fatigue (inversely related to T3t)
double Fatigue_score = 5.0 * (1.0 - T3_norm) + 2.0;

// BMD risk (T-score decline with excess LT4)
double BMD_risk = -0.5 * (LT4c / 100.0);  // T-score units

// Add residual error
double fT4_obs  = fT4_pmolL  * (1.0 + EPS(1));
double TSH_obs  = TSH_clinical * (1.0 + EPS(2));
double Ab_obs   = AntiTPOAb  * (1.0 + EPS(3));

$CAPTURE
fT4_pmolL fT3_pmolL TSH_clinical AntiTPOAb AntiTgAb DmgThy
Th1 Treg Bcell ThyVol_rel HypoScore HR_effect LDL_rel Fatigue_score BMD_risk
T4p T3p T4t T3t LT4c LiT3c Se
fT4_obs TSH_obs Ab_obs
'

# Compile model
ht_mod <- mrgsolve::mcode("hashimoto_qsp", ht_model_code)

cat("Model compiled successfully!\n")
cat("Number of ODEs:", length(init(ht_mod)), "\n")
cat("Parameters:", length(param(ht_mod)), "\n")

##############################################################################
# SCENARIO SIMULATIONS
##############################################################################

# Simulation time: 5 years with weekly resolution
sim_time <- seq(0, 5*365, by = 7)  # days

# ---- Helper: simulate one scenario ----------------------------------------
run_scenario <- function(mod, dose_LT4=0, dose_Se=0, dose_LiT3=0, mmI=0,
                         label="Untreated") {
  mod %>%
    param(LT4_dose=dose_LT4, Se_dose=dose_Se,
          LiT3_dose=dose_LiT3, MMI_block=mmI) %>%
    mrgsim(end=5*365, delta=7, obsonly=TRUE) %>%
    as_tibble() %>%
    mutate(Scenario=label, Time_years = time / 365)
}

# ---- 7 Treatment Scenarios -------------------------------------------------
cat("\nRunning 7 treatment scenarios...\n")

scen1 <- run_scenario(ht_mod, 0, 0, 0, 0,
                       "1. Untreated Hashimoto's")

scen2 <- run_scenario(ht_mod, 100, 0, 0, 0,
                       "2. LT4 monotherapy (100 μg/day)")

scen3 <- run_scenario(ht_mod, 0, 200, 0, 0,
                       "3. Selenium (200 μg/day)")

scen4 <- run_scenario(ht_mod, 100, 200, 0, 0,
                       "4. LT4 + Selenium combination")

scen5 <- run_scenario(ht_mod, 100, 0, 7.5, 0,
                       "5. LT4 + Liothyronine (100+7.5 μg/day)")

scen6 <- run_scenario(ht_mod, 175, 0, 0, 0,
                       "6. High-dose LT4 (175 μg/day, TSH-suppressive)")

scen7 <- run_scenario(ht_mod, 75, 200, 0, 0,
                       "7. Early LT4 + Selenium (75 μg/day)")

all_scen <- bind_rows(scen1, scen2, scen3, scen4, scen5, scen6, scen7)

##############################################################################
# VISUALIZATION
##############################################################################

scenario_colors <- c(
  "1. Untreated Hashimoto's"            = "#E74C3C",
  "2. LT4 monotherapy (100 μg/day)"     = "#3498DB",
  "3. Selenium (200 μg/day)"            = "#2ECC71",
  "4. LT4 + Selenium combination"       = "#9B59B6",
  "5. LT4 + Liothyronine (100+7.5 μg/day)" = "#E67E22",
  "6. High-dose LT4 (175 μg/day, TSH-suppressive)" = "#1ABC9C",
  "7. Early LT4 + Selenium (75 μg/day)" = "#F39C12"
)

# ---- Panel 1: TSH over time ------------------------------------------------
p1 <- all_scen %>%
  ggplot(aes(x=Time_years, y=TSH_clinical, color=Scenario)) +
  geom_line(linewidth=0.9, alpha=0.9) +
  geom_hline(yintercept=c(0.4, 4.0), linetype="dashed", color="grey40", linewidth=0.6) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=0.4, ymax=4.0, alpha=0.05, fill="#2ECC71") +
  annotate("text", x=4.8, y=2.0, label="Normal TSH\n(0.4–4.0)", color="grey40", size=3) +
  scale_color_manual(values=scenario_colors) +
  scale_y_log10() +
  labs(title="A  TSH Dynamics (mIU/L)", x=NULL, y="TSH (mIU/L, log)") +
  theme_bw(base_size=10) +
  theme(legend.position="none")

# ---- Panel 2: Free T4 -------------------------------------------------------
p2 <- all_scen %>%
  ggplot(aes(x=Time_years, y=fT4_pmolL, color=Scenario)) +
  geom_line(linewidth=0.9, alpha=0.9) +
  geom_hline(yintercept=c(9, 23), linetype="dashed", color="grey40", linewidth=0.6) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=9, ymax=23, alpha=0.05, fill="#2ECC71") +
  scale_color_manual(values=scenario_colors) +
  labs(title="B  Free T4 (pmol/L)", x=NULL, y="fT4 (pmol/L)") +
  theme_bw(base_size=10) +
  theme(legend.position="none")

# ---- Panel 3: Anti-TPO Antibody --------------------------------------------
p3 <- all_scen %>%
  ggplot(aes(x=Time_years, y=AntiTPOAb, color=Scenario)) +
  geom_line(linewidth=0.9, alpha=0.9) +
  geom_hline(yintercept=34, linetype="dashed", color="grey40", linewidth=0.6) +
  annotate("text", x=4.8, y=50, label="Upper normal\n34 IU/mL", color="grey40", size=3) +
  scale_color_manual(values=scenario_colors) +
  labs(title="C  Anti-TPO Antibody (IU/mL)", x=NULL, y="Anti-TPO Ab (IU/mL)") +
  theme_bw(base_size=10) +
  theme(legend.position="none")

# ---- Panel 4: Thyroid Damage -----------------------------------------------
p4 <- all_scen %>%
  ggplot(aes(x=Time_years, y=DmgThy*100, color=Scenario)) +
  geom_line(linewidth=0.9, alpha=0.9) +
  scale_color_manual(values=scenario_colors) +
  labs(title="D  Thyroid Damage (% loss)", x=NULL, y="Damage (%)") +
  theme_bw(base_size=10) +
  theme(legend.position="none")

# ---- Panel 5: Hypothyroid Score --------------------------------------------
p5 <- all_scen %>%
  ggplot(aes(x=Time_years, y=HypoScore, color=Scenario)) +
  geom_line(linewidth=0.9, alpha=0.9) +
  scale_color_manual(values=scenario_colors) +
  labs(title="E  Hypothyroid Symptom Score (0-10)", x="Time (years)", y="Score") +
  theme_bw(base_size=10) +
  theme(legend.position="none")

# ---- Panel 6: Th1/Treg balance --------------------------------------------
p6 <- all_scen %>%
  ggplot(aes(x=Time_years, y=Th1/Treg, color=Scenario)) +
  geom_line(linewidth=0.9, alpha=0.9) +
  scale_color_manual(values=scenario_colors) +
  labs(title="F  Th1/Treg Ratio (immune balance)", x="Time (years)", y="Th1/Treg") +
  theme_bw(base_size=10) +
  theme(legend.position="bottom", legend.text=element_text(size=7),
        legend.key.size=unit(0.4, "cm")) +
  guides(color=guide_legend(ncol=2))

# ---- Combine panels --------------------------------------------------------
fig_main <- (p1 + p2) / (p3 + p4) / (p5 + p6) +
  plot_annotation(
    title = "Hashimoto's Thyroiditis QSP Model — 7-Scenario Comparison (5-Year Simulation)",
    subtitle = "Levothyroxine · Selenium · Liothyronine · High-dose LT4 · Combination Regimens",
    theme = theme(plot.title = element_text(face="bold", size=13),
                  plot.subtitle = element_text(size=10, color="grey40"))
  )

cat("Main figure created.\n")

##############################################################################
# SENSITIVITY ANALYSIS: Selenium Dose-Response
##############################################################################

Se_doses <- c(0, 50, 100, 150, 200, 300)
Se_dose_response <- lapply(Se_doses, function(sd) {
  ht_mod %>%
    param(LT4_dose=100, Se_dose=sd) %>%
    mrgsim(end=365, delta=30, obsonly=TRUE) %>%
    as_tibble() %>%
    filter(time == 365) %>%
    mutate(Se_dose_ugday = sd)
}) %>% bind_rows()

p_se_dr <- Se_dose_response %>%
  ggplot(aes(x=Se_dose_ugday)) +
  geom_line(aes(y=AntiTPOAb/max(AntiTPOAb)*100, color="Anti-TPO Ab (%)"), linewidth=1.2) +
  geom_point(aes(y=AntiTPOAb/max(AntiTPOAb)*100, color="Anti-TPO Ab (%)"), size=3) +
  scale_color_manual(values=c("Anti-TPO Ab (%)" = "#E74C3C")) +
  labs(title="Selenium Dose-Response on Anti-TPO Antibody (at 1 year, with LT4 100 μg/day)",
       x="Selenium Dose (μg/day)", y="Anti-TPO Ab (% of no-Se baseline)",
       color=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom")

##############################################################################
# POPULATION VARIABILITY (n=50 virtual patients)
##############################################################################

set.seed(42)

pop_sim <- ht_mod %>%
  param(LT4_dose=100, Se_dose=200) %>%
  idata_set(data.frame(ID=1:50)) %>%
  mrgsim(end=365, delta=30, obsonly=TRUE) %>%
  as_tibble()

p_pop <- pop_sim %>%
  ggplot(aes(x=time/30, y=TSH_clinical, group=ID)) +
  geom_line(alpha=0.25, color="#3498DB") +
  geom_smooth(aes(group=NULL), method="loess", se=FALSE, color="#1A237E", linewidth=1.5) +
  geom_hline(yintercept=c(0.4, 4.0), linetype="dashed", color="red") +
  labs(title="Population Variability: TSH Under LT4 100 μg + Se 200 μg (n=50 virtual patients)",
       x="Time (months)", y="TSH (mIU/L)") +
  theme_bw(base_size=11)

##############################################################################
# LT4 DOSE-RESPONSE TABLE
##############################################################################

LT4_doses <- c(0, 25, 50, 75, 100, 125, 150, 175)

lt4_response <- lapply(LT4_doses, function(ld) {
  ht_mod %>%
    param(LT4_dose=ld, Se_dose=0) %>%
    mrgsim(end=365, delta=365, obsonly=TRUE) %>%
    as_tibble() %>%
    filter(time == 365) %>%
    transmute(
      `LT4 Dose (μg/day)` = ld,
      `TSH (mIU/L)` = round(TSH_clinical, 2),
      `fT4 (pmol/L)` = round(fT4_pmolL, 1),
      `fT3 (pmol/L)` = round(fT3_pmolL, 2),
      `HypoScore` = round(HypoScore, 1),
      `HR (bpm)` = round(HR_effect, 0),
      `LDL (mmol/L)` = round(LDL_rel, 2),
      `BMD risk` = round(BMD_risk, 3)
    )
}) %>% bind_rows()

cat("\n=== LT4 Dose-Response Summary (at 1 year) ===\n")
print(lt4_response)

##############################################################################
# SELENIUM MECHANISTIC VALIDATION
# (Replication of Gärtner 2002 finding: ~50% Ab reduction at 200 μg/day)
##############################################################################

validation_sim <- ht_mod %>%
  param(LT4_dose=0, Se_dose=200) %>%
  mrgsim(end=270, delta=30, obsonly=TRUE) %>%  # 9-month RCT
  as_tibble()

baseline_Ab <- validation_sim$AntiTPOAb[1]
final_Ab    <- validation_sim$AntiTPOAb[nrow(validation_sim)]
reduction_pct <- (1 - final_Ab / baseline_Ab) * 100

cat("\n=== Selenium Validation (Gärtner 2002 RCT replication) ===\n")
cat("Baseline Anti-TPO Ab:", round(baseline_Ab, 1), "IU/mL\n")
cat("9-month Anti-TPO Ab :", round(final_Ab, 1), "IU/mL\n")
cat("Reduction:           ", round(reduction_pct, 1), "%\n")
cat("(Expected from RCT:  ~50% reduction)\n")

##############################################################################
# SUMMARY FUNCTION
##############################################################################

summarize_endpoints <- function(sim_data) {
  sim_data %>%
    group_by(Scenario) %>%
    summarise(
      TSH_1yr    = round(mean(TSH_clinical[Time_years >= 0.9 & Time_years <= 1.1]), 2),
      fT4_1yr    = round(mean(fT4_pmolL[Time_years >= 0.9 & Time_years <= 1.1]), 1),
      AntiTPO_1yr = round(mean(AntiTPOAb[Time_years >= 0.9 & Time_years <= 1.1]), 0),
      DmgThy_5yr = round(mean(DmgThy[Time_years >= 4.9]),  3),
      HypoScore_1yr = round(mean(HypoScore[Time_years >= 0.9 & Time_years <= 1.1]), 1),
      .groups = "drop"
    )
}

summary_table <- summarize_endpoints(all_scen)
cat("\n=== 7-Scenario Clinical Endpoint Summary ===\n")
print(summary_table)

cat("\n✔ Hashimoto's Thyroiditis QSP model simulation complete.\n")
cat("Output objects: ht_mod, all_scen, fig_main, p_se_dr, p_pop, lt4_response, summary_table\n")
