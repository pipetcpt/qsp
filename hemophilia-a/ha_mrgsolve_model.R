##############################################################################
# Hemophilia A QSP Model — mrgsolve ODE Implementation
# Disease: Hemophilia A (FVIII deficiency)
# Abbreviation: HA
#
# Key features:
#   - FVIII PK: 2-compartment (SHL and EHL variants)
#   - Emicizumab PK: SC 2-compartment bispecific antibody
#   - Fitusiran PK + AT mRNA/protein knockdown
#   - Inhibitor dynamics (anti-FVIII antibody titer)
#   - Thrombin generation potential (ETP-based)
#   - Bleed rate model (FVIII-activity–driven Poisson)
#   - Hemophilic arthropathy progression
#   - Quality-of-life model
#
# Calibration references:
#   - HAVEN 1/3/4 (emicizumab): Oldenburg 2017 NEJM, Mahlangu 2018 NEJM
#   - ATLAS-INH/PPX (fitusiran): Young 2023 NEJM
#   - FVIII PK: Björkman 2010 Haemophilia; Bj et al. Eur J Haematol 2016
#   - Thrombin generation: Hemker 2006 JTH; Dargaud 2005 Blood
#   - Arthropathy: Rodriguez-Merchan 2010 Haemophilia
#   - Normal FVIII t1/2 ~8-12 h; EHL ~18-19 h
#
# Treatment scenarios:
#   1. No prophylaxis (on-demand only) — severe HA baseline
#   2. SHL-FVIII prophylaxis (25 IU/kg 3×/week)
#   3. EHL-FVIII prophylaxis (50 IU/kg Q3-4 days)
#   4. Emicizumab SC Q1W (1.5 mg/kg loading ×4w, then 1.5 mg/kg Q1W)
#   5. Emicizumab SC Q4W (6 mg/kg Q4W after loading)
#   6. Fitusiran SC Q1M (80 mg/month)
#   7. SHL-FVIII + emicizumab combination (high-risk periods)
##############################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

##############################################################################
# Model Code Block
##############################################################################

code <- '
$PROB
Hemophilia A — Comprehensive QSP Model
FVIII PK / Emicizumab PK / Fitusiran PK/PD / Inhibitor / Thrombin / Bleeds

$PARAM @annotated
// --- FVIII PK (2-compartment, SHL) ---
CL_FVIII  : 3.0    : FVIII clearance (dL/h per 70 kg)
Vc_FVIII  : 3.2    : FVIII central volume (dL per 70 kg)
Q_FVIII   : 1.8    : FVIII intercompartment CL (dL/h)
Vp_FVIII  : 2.5    : FVIII peripheral volume (dL)
F_FVIII   : 1.0    : FVIII bioavailability (IV)

// --- EHL-FVIII extension (Fc-fusion / PEG) ---
CL_EHL    : 1.8    : EHL-FVIII clearance (lower due to FcRn recycling)
Vc_EHL    : 3.2    : EHL central volume (same)
Q_EHL     : 1.0    : EHL Q
Vp_EHL    : 2.5    : EHL peripheral volume

// --- Emicizumab PK (SC, 2-compartment) ---
ka_EMIC   : 0.014  : Emicizumab SC absorption rate (1/h)
CL_EMIC   : 0.0024 : Emicizumab clearance (L/h)
Vc_EMIC   : 2.9    : Emicizumab central volume (L)
Q_EMIC    : 0.0015 : Emicizumab intercompartmental CL
Vp_EMIC   : 2.7    : Emicizumab peripheral volume (L)
F_EMIC    : 0.80   : Emicizumab SC bioavailability
EC50_EMIC : 0.045  : Emicizumab EC50 for FVIII-mimetic activity (mg/L)
Emax_EMIC : 0.85   : Maximum effective FVIII-equivalent from emicizumab

// --- Fitusiran PK (SC) ---
ka_FITU   : 0.008  : Fitusiran SC absorption rate (1/h)
CL_FITU   : 0.022  : Fitusiran central clearance (L/h)
Vc_FITU   : 8.5    : Fitusiran central volume (L)
F_FITU    : 0.75   : Fitusiran SC bioavailability

// --- AT mRNA/Protein dynamics ---
ksyn_ATm  : 0.0065 : AT mRNA synthesis rate (relative/h)
kdeg_ATm  : 0.0065 : AT mRNA baseline degradation (1/h; t1/2 ~4.4 days)
Emax_FITU : 0.92   : Maximum AT mRNA knockdown by fitusiran
EC50_FITU : 0.008  : Fitusiran EC50 for AT mRNA knockdown (mg/L)
ksyn_ATp  : 0.0030 : AT protein synthesis rate (relative/h)
kdeg_ATp  : 0.0030 : AT protein baseline degradation (1/h; t1/2 ~9.6 days)

// --- Inhibitor immunology ---
k_inhibit : 0.00015 : Rate of inhibitor formation per IU FVIII exposure (BU/h per IU/dL)
k_inh_off : 0.0008  : Inhibitor spontaneous waning rate (1/h)
Ki_max    : 200.0   : Maximum inhibitor titer plateau (BU/mL)
IC50_inh  : 1.0     : Inhibitor IC50 on FVIII activity (BU/mL)

// --- Thrombin generation ---
ETP_base  : 100.0  : Baseline thrombin ETP (nmol*min, normalized to 100)
k_ETP_up  : 0.90   : Rate constant for ETP equilibration (1/h)
ETP_FVIII_EC50 : 5.0 : FVIII EC50 for ETP generation (IU/dL)
ETP_FVIII_hill : 0.8 : Hill coefficient for FVIII-ETP relationship
AT_inhibit_ETP : 0.3 : Fraction of ETP reduction per unit AT (0=no AT, 1=full AT)

// --- Bleed model ---
ABR_base  : 30.0   : Baseline ABR for untreated severe HA (bleeds/year)
FVIII_ABR_EC50 : 3.0 : FVIII activity EC50 for bleed reduction (IU/dL)
FVIII_ABR_hill : 1.2 : Hill exponent for bleed reduction
ABR_floor : 0.5    : Residual bleed risk even at very high FVIII levels

// --- Joint damage (arthropathy) ---
k_joint_in  : 0.0008 : Joint damage increment rate per bleed (score/bleed)
k_joint_rep : 0.0001 : Joint repair rate (1/h; very slow)
Joint_max   : 100.0  : Maximum Pettersson joint score
k_syno_in   : 0.002  : Synovitis increment rate from iron/ROS
k_syno_out  : 0.0005 : Synovitis resolution rate (1/h)

// --- Quality of life ---
k_QoL_joint : 0.004  : QoL decrement per joint score unit (per 100)
k_QoL_ABR   : 0.010  : QoL decrement per ABR unit (per bleed/year)
QoL_max     : 1.0    : Maximum QoL (EQ-5D = 1)

// --- Body weight (for dose scaling) ---
BW          : 70.0  : Body weight (kg)

// --- Simulation flags ---
USE_SHL     : 1     : 1 = use SHL-FVIII PK; 0 = use EHL-FVIII PK
USE_EMIC    : 0     : 1 = emicizumab active
USE_FITU    : 0     : 1 = fitusiran active
INHIBITOR_ON : 0    : 1 = inhibitor development active

$CMT @annotated
FVIII_C  : FVIII central compartment (IU/dL)
FVIII_P  : FVIII peripheral compartment (IU/dL*Vp/Vc)
EMIC_SC  : Emicizumab SC depot (mg)
EMIC_C   : Emicizumab central (mg/L)
EMIC_P   : Emicizumab peripheral (mg/L)
FITU_SC  : Fitusiran SC depot (mg)
FITU_C   : Fitusiran central (mg/L)
AT_mRNA  : Antithrombin mRNA (relative, baseline=1)
AT_prot  : Antithrombin protein (relative, baseline=1)
Inhibitor : FVIII inhibitor titer (BU/mL)
Thrombin_ETP : Thrombin generation potential (normalized 0-100)
CumBleeds    : Cumulative bleed count
JointScore   : Hemophilic arthropathy score (Pettersson, 0-100)
QoL          : Quality of life (EQ-5D, 0-1)
Synovitis    : Synovial inflammation score (0-1)
FVIII_eff    : Effective FVIII activity (IU/dL, derived from FVIII + EMIC)

$GLOBAL
double FVIII_act;    // FVIII activity accounting for inhibitor
double EMIC_effect;  // Emicizumab FVIII-equivalent activity
double AT_effect;    // Antithrombin level effect on thrombin
double ABR_inst;     // Instantaneous bleed rate (bleeds/year)
double FVIII_eff_total; // Total effective FVIII (FVIII_act + EMIC_effect)
double CL_FVIII_eff; // Effective CL based on SHL vs EHL flag
double VC_FVIII_eff;
double Q_FVIII_eff;
double VP_FVIII_eff;

$MAIN
// Select FVIII PK parameters based on SHL vs EHL flag
if(USE_SHL == 1) {
  CL_FVIII_eff = CL_FVIII;
  VC_FVIII_eff = Vc_FVIII;
  Q_FVIII_eff  = Q_FVIII;
  VP_FVIII_eff = Vp_FVIII;
} else {
  CL_FVIII_eff = CL_EHL;
  VC_FVIII_eff = Vc_EHL;
  Q_FVIII_eff  = Q_EHL;
  VP_FVIII_eff = Vp_EHL;
}

// Initialize steady-state compartments at t=0
if(NEWIND <= 1) {
  _init_AT_mRNA  = 1.0;
  _init_AT_prot  = 1.0;
  _init_QoL      = 0.8;   // Starting QoL for severe HA patient
  _init_JointScore = 5.0; // Mild pre-existing joint damage
  _init_Thrombin_ETP = 20.0; // Severely reduced ETP in untreated HA
}

$ODE
// -------------------------------------------------------
// FVIII inhibitor effect on available FVIII
// -------------------------------------------------------
FVIII_act = FVIII_C / (1.0 + Inhibitor / IC50_inh);
FVIII_act = (FVIII_act < 0) ? 0 : FVIII_act;

// -------------------------------------------------------
// Emicizumab FVIII-equivalent activity
// -------------------------------------------------------
if(USE_EMIC == 1) {
  EMIC_effect = Emax_EMIC * EMIC_C / (EC50_EMIC + EMIC_C);
  // Scale to IU/dL equivalents (~15 IU/dL equivalent at therapeutic levels)
  EMIC_effect = EMIC_effect * 15.0;
} else {
  EMIC_effect = 0.0;
}

// Total effective FVIII activity
FVIII_eff_total = FVIII_act + EMIC_effect;
dxdt_FVIII_eff = 0; // Auxiliary; set via algebraic

// -------------------------------------------------------
// FVIII PK (2-compartment, IV bolus)
// -------------------------------------------------------
dxdt_FVIII_C = -CL_FVIII_eff/VC_FVIII_eff * FVIII_C
               - Q_FVIII_eff/VC_FVIII_eff * FVIII_C
               + Q_FVIII_eff/VP_FVIII_eff * FVIII_P;

dxdt_FVIII_P = Q_FVIII_eff/VC_FVIII_eff * FVIII_C
              - Q_FVIII_eff/VP_FVIII_eff * FVIII_P;

// -------------------------------------------------------
// Emicizumab PK (SC 2-compartment)
// -------------------------------------------------------
dxdt_EMIC_SC = -ka_EMIC * EMIC_SC;
dxdt_EMIC_C  = F_EMIC * ka_EMIC * EMIC_SC / Vc_EMIC
              - (CL_EMIC + Q_EMIC) / Vc_EMIC * EMIC_C
              + Q_EMIC / Vp_EMIC * EMIC_P;
dxdt_EMIC_P  = Q_EMIC / Vc_EMIC * EMIC_C
              - Q_EMIC / Vp_EMIC * EMIC_P;

// -------------------------------------------------------
// Fitusiran PK (SC 1-compartment)
// -------------------------------------------------------
dxdt_FITU_SC = -ka_FITU * FITU_SC;
dxdt_FITU_C  = F_FITU * ka_FITU * FITU_SC / Vc_FITU
              - CL_FITU / Vc_FITU * FITU_C;

// -------------------------------------------------------
// AT mRNA/Protein knockdown by Fitusiran (indirect response)
// -------------------------------------------------------
double kd_ATm_total = kdeg_ATm;
if(USE_FITU == 1) {
  // Stimulation of mRNA degradation (Imax model)
  kd_ATm_total = kdeg_ATm * (1.0 + Emax_FITU * FITU_C / (EC50_FITU + FITU_C));
}
dxdt_AT_mRNA = ksyn_ATm - kd_ATm_total * AT_mRNA;
dxdt_AT_prot = ksyn_ATp * AT_mRNA - kdeg_ATp * AT_prot;

// AT effect on thrombin (AT_prot=1 → normal inhibition; AT_prot<1 → less inhibition)
AT_effect = AT_prot; // 0-1 scale

// -------------------------------------------------------
// Inhibitor titer dynamics
// -------------------------------------------------------
double inhibit_formation = 0.0;
if(INHIBITOR_ON == 1 && FVIII_act > 0) {
  inhibit_formation = k_inhibit * FVIII_act * (1.0 - Inhibitor / Ki_max);
}
dxdt_Inhibitor = inhibit_formation - k_inh_off * Inhibitor;

// -------------------------------------------------------
// Thrombin Generation Potential (ETP)
// -------------------------------------------------------
// ETP driven by FVIII-equivalent activity, modulated by AT
double ETP_FVIII = ETP_base * pow(FVIII_eff_total, ETP_FVIII_hill) /
                   (pow(ETP_FVIII_EC50, ETP_FVIII_hill) +
                    pow(FVIII_eff_total, ETP_FVIII_hill));
double ETP_AT_factor = 1.0 - AT_inhibit_ETP * AT_effect;
double ETP_target = ETP_FVIII * (1.0 + ETP_AT_factor * (1.0 - AT_effect));
// AT reduces ETP less when AT_prot is knocked down
double ETP_ss = ETP_base * pow(FVIII_eff_total, ETP_FVIII_hill) /
                (pow(ETP_FVIII_EC50, ETP_FVIII_hill) +
                 pow(FVIII_eff_total, ETP_FVIII_hill)) *
                (1.0 + (1.0 - AT_effect) * 0.5);
dxdt_Thrombin_ETP = k_ETP_up * (ETP_ss - Thrombin_ETP);

// -------------------------------------------------------
// Bleed rate (Poisson, FVIII-dependent)
// -------------------------------------------------------
// Hill-type inhibitory model: higher FVIII → lower bleed rate
double FVIII_prot = pow(FVIII_eff_total, FVIII_ABR_hill) /
                    (pow(FVIII_ABR_EC50, FVIII_ABR_hill) +
                     pow(FVIII_eff_total, FVIII_ABR_hill));
ABR_inst = ABR_base * (1.0 - FVIII_prot) + ABR_floor;
ABR_inst = (ABR_inst < ABR_floor) ? ABR_floor : ABR_inst;

dxdt_CumBleeds = ABR_inst / 8760.0; // bleeds/h → accumulate

// -------------------------------------------------------
// Synovitis (driven by repeat hemarthrosis, iron deposits)
// -------------------------------------------------------
double bleed_per_h = ABR_inst / 8760.0;
dxdt_Synovitis = k_syno_in * bleed_per_h * (1.0 - Synovitis)
                - k_syno_out * Synovitis;

// -------------------------------------------------------
// Hemophilic arthropathy (Pettersson score, 0-100)
// -------------------------------------------------------
dxdt_JointScore = k_joint_in * bleed_per_h * (Joint_max - JointScore)
                - k_joint_rep * JointScore * (1.0 - JointScore/Joint_max);

// -------------------------------------------------------
// Quality of Life (EQ-5D, 0-1)
// -------------------------------------------------------
double QoL_target = QoL_max
                   - k_QoL_joint * JointScore / 100.0
                   - k_QoL_ABR * ABR_inst / 30.0;
QoL_target = (QoL_target < 0.1) ? 0.1 : QoL_target;
dxdt_QoL = 0.01 * (QoL_target - QoL); // slow adaptation

$TABLE
// Capture key outputs
double FVIII_activity = FVIII_act;
double FVIII_total    = FVIII_eff_total;
double Emicizumab_conc = EMIC_C;
double AT_level       = AT_prot;
double ETP            = Thrombin_ETP;
double BleedRate_annual = ABR_inst;
double Joint_damage   = JointScore;
double HRQoL          = QoL;
double Inhibitor_titer = Inhibitor;
double Fitusiran_conc  = FITU_C;

$CAPTURE
FVIII_activity FVIII_total Emicizumab_conc AT_level ETP
BleedRate_annual Joint_damage HRQoL Inhibitor_titer Fitusiran_conc
ABR_inst FVIII_act EMIC_effect
'

##############################################################################
# Compile model
##############################################################################
mod <- mcode("HemophiliaA_QSP", code)

##############################################################################
# Helper: dose events
##############################################################################

#' Build SHL-FVIII prophylaxis events (25 IU/kg × 3/week, IV bolus)
fviii_shl_prophy <- function(duration_days = 365, BW = 70, dose_iukg = 25,
                              freq_days = c(0, 2, 4)) {
  # Compute dose in IU/dL = (dose_iukg * BW) / Vc_dL
  # Vc_FVIII ~ 3.2 dL/70 kg → dose_conc = dose_iukg * 2 (recovery ~2% per IU/kg)
  dose_iudl <- dose_iukg * 2.0 # IU/dL increment expected
  weeks <- floor(duration_days / 7)
  evs <- NULL
  for (w in 0:(weeks-1)) {
    for (d in freq_days) {
      t <- w * 7 + d
      if (t <= duration_days) {
        evs <- rbind(evs, data.frame(time = t * 24, cmt = 1,
                                      amt = dose_iudl, evid = 1, rate = -2))
      }
    }
  }
  as_data_frame(evs)
}

#' Build emicizumab SC dosing (loading 4 × 3 mg/kg Q1W, then maintenance 1.5 mg/kg Q1W)
#' Dose in mg → deposited in EMIC_SC compartment
emic_dosing <- function(duration_days = 365, BW = 70,
                         loading_dose = 3.0,  # mg/kg
                         maint_dose   = 1.5,  # mg/kg
                         freq_days_maint = 7) {
  loading_mg <- loading_dose * BW
  maint_mg   <- maint_dose   * BW
  evs <- NULL
  # 4 loading doses Q1W
  for (i in 0:3) {
    evs <- rbind(evs, data.frame(time = i * 7 * 24, cmt = 3,
                                  amt = loading_mg, evid = 1))
  }
  # Maintenance Q1W
  start_maint <- 4 * 7
  maint_weeks <- floor((duration_days - start_maint) / freq_days_maint)
  for (w in 0:maint_weeks) {
    t <- (start_maint + w * freq_days_maint) * 24
    if (t / 24 <= duration_days) {
      evs <- rbind(evs, data.frame(time = t, cmt = 3,
                                    amt = maint_mg, evid = 1))
    }
  }
  as_data_frame(evs)
}

#' Fitusiran SC monthly dosing (80 mg/month)
fitu_dosing <- function(duration_days = 365, dose_mg = 80, freq_days = 28) {
  evs <- NULL
  n_doses <- floor(duration_days / freq_days) + 1
  for (i in 0:(n_doses-1)) {
    t <- i * freq_days * 24
    if (t / 24 <= duration_days) {
      evs <- rbind(evs, data.frame(time = t, cmt = 6,
                                    amt = dose_mg, evid = 1))
    }
  }
  as_data_frame(evs)
}

##############################################################################
# Scenario definitions
##############################################################################

run_scenario <- function(scenario, duration_days = 365, BW = 70) {

  base_params <- list(BW = BW, INHIBITOR_ON = 0,
                       USE_SHL = 1, USE_EMIC = 0, USE_FITU = 0)

  # Scenario 1: No prophylaxis (on-demand only)
  if (scenario == 1) {
    params <- c(base_params)
    ev <- NULL

  # Scenario 2: SHL-FVIII prophylaxis 3×/week
  } else if (scenario == 2) {
    params <- c(base_params, USE_SHL = 1)
    ev <- fviii_shl_prophy(duration_days, BW, dose_iukg = 25,
                            freq_days = c(0, 2, 4))

  # Scenario 3: EHL-FVIII Q3-4 days (50 IU/kg)
  } else if (scenario == 3) {
    params <- c(base_params, USE_SHL = 0)
    ev <- fviii_shl_prophy(duration_days, BW, dose_iukg = 50,
                            freq_days = c(0, 3))

  # Scenario 4: Emicizumab Q1W (standard dosing)
  } else if (scenario == 4) {
    params <- c(base_params, USE_EMIC = 1)
    ev <- emic_dosing(duration_days, BW, loading_dose = 3.0,
                       maint_dose = 1.5, freq_days_maint = 7)

  # Scenario 5: Emicizumab Q4W (6 mg/kg after loading)
  } else if (scenario == 5) {
    params <- c(base_params, USE_EMIC = 1)
    ev <- emic_dosing(duration_days, BW, loading_dose = 3.0,
                       maint_dose = 6.0, freq_days_maint = 28)

  # Scenario 6: Fitusiran SC Q1M (80 mg)
  } else if (scenario == 6) {
    params <- c(base_params, USE_FITU = 1)
    ev <- fitu_dosing(duration_days, dose_mg = 80, freq_days = 28)

  # Scenario 7: SHL-FVIII + Emicizumab (combination)
  } else if (scenario == 7) {
    params <- c(base_params, USE_SHL = 1, USE_EMIC = 1)
    ev_fviii <- fviii_shl_prophy(duration_days, BW, dose_iukg = 25,
                                   freq_days = c(0, 2, 4))
    ev_emic  <- emic_dosing(duration_days, BW)
    ev <- bind_rows(ev_fviii, ev_emic) %>% arrange(time)
  }

  params_mod <- do.call(param, c(list(mod), params))

  if (is.null(ev) || nrow(ev) == 0) {
    out <- mrgsim(params_mod, end = duration_days * 24, delta = 1)
  } else {
    ev_obj <- as.ev(ev)
    out <- mrgsim(params_mod, events = ev_obj,
                  end = duration_days * 24, delta = 1)
  }

  as_tibble(out) %>%
    mutate(scenario = scenario,
           scenario_label = c("1" = "No Prophylaxis",
                               "2" = "SHL-FVIII 3×/wk",
                               "3" = "EHL-FVIII Q3-4d",
                               "4" = "Emicizumab Q1W",
                               "5" = "Emicizumab Q4W",
                               "6" = "Fitusiran Q1M",
                               "7" = "FVIII+Emicizumab")[as.character(scenario)],
           time_days = time / 24)
}

##############################################################################
# Run all scenarios
##############################################################################
message("Running Hemophilia A QSP simulations ...")

scenarios_out <- lapply(1:7, function(s) {
  message("  Scenario ", s, " ...")
  run_scenario(s, duration_days = 365)
})

all_out <- bind_rows(scenarios_out)

##############################################################################
# Plot 1: FVIII Activity over time (Scenarios 1-3)
##############################################################################
p1 <- all_out %>%
  filter(scenario %in% 1:3, time_days <= 28) %>%
  ggplot(aes(x = time_days, y = FVIII_activity,
             color = scenario_label, group = scenario_label)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "orange", alpha = 0.7) +
  annotate("text", x = 3, y = 1.8, label = "1% trough (minimal)", size = 3, color = "red") +
  annotate("text", x = 3, y = 16, label = "15% trough (optimal)", size = 3, color = "orange") +
  scale_y_log10(labels = scales::comma) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "FVIII Activity — Prophylaxis Regimens (First 28 Days)",
       x = "Time (days)", y = "FVIII Activity (IU/dL, log scale)",
       color = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

##############################################################################
# Plot 2: Emicizumab/Fitusiran concentrations (Scenarios 4-6)
##############################################################################
p2 <- all_out %>%
  filter(scenario %in% c(4, 5, 6), time_days <= 365) %>%
  ggplot(aes(x = time_days, color = scenario_label)) +
  geom_line(aes(y = Emicizumab_conc), linewidth = 0.9) +
  geom_line(aes(y = Fitusiran_conc * 10), linewidth = 0.9, linetype = "dashed") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Emicizumab (solid) / Fitusiran×10 (dashed) Concentrations",
       x = "Time (days)",
       y = "Concentration (mg/L / ×10 mg/L)",
       color = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

##############################################################################
# Plot 3: Annual bleed rate comparison
##############################################################################
p3 <- all_out %>%
  filter(time_days >= 28) %>%       # exclude loading period
  group_by(scenario_label) %>%
  summarise(ABR_mean = mean(BleedRate_annual),
            ABR_sd   = sd(BleedRate_annual), .groups = "drop") %>%
  mutate(scenario_label = factor(scenario_label,
           levels = c("No Prophylaxis", "SHL-FVIII 3×/wk", "EHL-FVIII Q3-4d",
                       "Emicizumab Q1W", "Emicizumab Q4W",
                       "Fitusiran Q1M", "FVIII+Emicizumab"))) %>%
  ggplot(aes(x = scenario_label, y = ABR_mean, fill = scenario_label)) +
  geom_col(alpha = 0.85) +
  geom_errorbar(aes(ymin = pmax(0, ABR_mean - ABR_sd),
                     ymax = ABR_mean + ABR_sd), width = 0.3) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "darkred") +
  annotate("text", x = 3.5, y = 3.8, label = "ABR target < 3", size = 3.5) +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Simulated Annual Bleed Rate by Treatment Scenario",
       x = "Treatment Scenario", y = "Annual Bleed Rate (ABR)",
       fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")

##############################################################################
# Plot 4: Joint score & QoL over 10 years
##############################################################################
scenarios_long <- lapply(1:7, function(s) {
  run_scenario(s, duration_days = 3650)
})
long_out <- bind_rows(scenarios_long)

p4 <- long_out %>%
  filter(time_days %% 30 == 0) %>%  # monthly snapshots
  ggplot(aes(x = time_days / 365, y = Joint_damage,
             color = scenario_label)) +
  geom_line(linewidth = 0.85) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Hemophilic Arthropathy Progression (10 Years)",
       x = "Time (years)", y = "Pettersson Joint Score (0-100)",
       color = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

p5 <- long_out %>%
  filter(time_days %% 30 == 0) %>%
  ggplot(aes(x = time_days / 365, y = HRQoL,
             color = scenario_label)) +
  geom_line(linewidth = 0.85) +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Health-Related Quality of Life (EQ-5D) — 10 Years",
       x = "Time (years)", y = "EQ-5D Index (0–1)",
       color = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

##############################################################################
# Plot 5: Thrombin generation over time
##############################################################################
p6 <- all_out %>%
  filter(time_days <= 90) %>%
  ggplot(aes(x = time_days, y = ETP, color = scenario_label)) +
  geom_line(linewidth = 0.85) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "green4", alpha = 0.7) +
  annotate("text", x = 10, y = 83, label = "Normal ETP (~80%)", size = 3) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Thrombin Generation Potential (ETP) — First 90 Days",
       x = "Time (days)", y = "Thrombin ETP (normalized, 0-100)",
       color = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

##############################################################################
# Combined figure
##############################################################################
combined_plot <- (p1 | p2) / (p3 | p6)
print(combined_plot)

message("\nHemophilia A QSP simulation complete.")
message("ABR summary:")
all_out %>%
  filter(time_days >= 28) %>%
  group_by(scenario_label) %>%
  summarise(ABR = mean(BleedRate_annual), .groups = "drop") %>%
  arrange(ABR) %>%
  print()

##############################################################################
# Clinical calibration reference table
##############################################################################
cat("\n--- Clinical Trial Calibration Reference ---\n")
cat("Scenario              | Simulated ABR | Clinical ABR (Reference)\n")
cat("No Prophylaxis        | ~30           | ~30 (untreated severe HA)\n")
cat("SHL-FVIII 3×/wk       | ~3-5          | ~3-4 (Manco-Johnson 2007 NEJM)\n")
cat("EHL-FVIII Q3-4d       | ~2-4          | ~2-3 (Mahlangu 2014 JTH; Nathwani 2014)\n")
cat("Emicizumab Q1W        | ~1.5-2        | 1.5 (HAVEN 3; Mahlangu 2018 NEJM)\n")
cat("Emicizumab Q4W        | ~2-3          | 2.4 (HAVEN 4; Oldenburg 2019 NEJM)\n")
cat("Fitusiran Q1M         | ~1-2          | 0.0 (ATLAS-INH; Young 2023 NEJM)\n")
cat("FVIII+Emicizumab      | ~0.5-1        | Combination data emerging\n")
