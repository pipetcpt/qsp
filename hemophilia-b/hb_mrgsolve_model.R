##############################################################################
# Hemophilia B QSP Model — mrgsolve ODE Implementation
# Disease: Hemophilia B (FIX deficiency, Christmas disease)
# Abbreviation: HB
#
# Key features:
#   - FIX PK: 2-compartment for SHL (BeneFIX/Rixubis) and EHL variants
#       (rFIX-Fc/Alprolix, rFIX-albumin/Idelvion, glycoPEGylated/Refixia)
#   - AAV gene therapy compartment: vector transduction -> episomal FIX-Padua
#       transgene expression -> capsid immune response -> transaminitis ->
#       corticosteroid-modulated expression decline
#   - Concizumab / Marstacimab PK + TFPI-neutralization -> thrombin rebalancing
#   - Fitusiran PK + antithrombin (AT) mRNA/protein knockdown
#   - Inhibitor immunology (anti-FIX antibody titer, Bethesda-scaled)
#   - Thrombin generation potential (ETP-based)
#   - Bleed rate model (FIX-activity-driven annualized bleed rate, ABR)
#   - Hemophilic arthropathy (joint score) & synovitis progression
#   - Quality-of-life (Haem-A-QoL-like, 0-1 scale)
#
# Calibration references (see hb_references.md for full list):
#   - FIX SHL PK: Björkman 2012 Haemophilia (population PK), BeneFIX/Rixubis label
#   - EHL rFIX-Fc: B-LONG (Powell 2013 NEJM) — t1/2 ~82h
#   - EHL rFIX-albumin: PROLONG-9FP (Santagostino 2016 Blood) — t1/2 ~102-104h
#   - GlycoPEGylated rFIX: pathfinder2 (Collins 2014 JTH) — t1/2 ~93h
#   - Concizumab: explorer7/8 (Shapiro 2019 JTH; Chowdary 2023)
#   - Marstacimab: BASIS (Pipe 2023 NEJM) — weekly SC
#   - Fitusiran: ATLAS-A/B, ATLAS-INH (Young 2023 NEJM)
#   - AAV gene therapy: HOPE-B / Hemgenix (Pipe 2023 NEJM); BENEGENE-2 / Beqvez
#       (Konkle 2024); FIX-Padua (R338L): Simioni 2009 NEJM (~8x specific activity)
#   - Arthropathy/ABR modeling framework adapted from Hemophilia A QSP precedent
#
# Treatment scenarios (10):
#   1. No prophylaxis (on-demand only) — severe HB baseline
#   2. SHL-rFIX prophylaxis (40 IU/kg 2x/week IV)
#   3. EHL rFIX-Fc prophylaxis (Alprolix, 50 IU/kg Q7-10 days IV)
#   4. EHL rFIX-albumin prophylaxis (Idelvion, 75 IU/kg Q14 days IV)
#   5. GlycoPEGylated rFIX prophylaxis (Refixia, 40 IU/kg Q7 days IV)
#   6. Concizumab SC daily prophylaxis (non-factor rebalancing)
#   7. Marstacimab SC weekly prophylaxis (non-factor rebalancing)
#   8. Fitusiran SC monthly prophylaxis (AT-lowering rebalancing)
#   9. AAV gene therapy single IV dose (etranacogene dezaparvovec-like) + steroid taper
#  10. Inhibitor-positive patient on ITI protocol + bypassing-agent bleed management
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
Hemophilia B — Comprehensive QSP Model
FIX PK (SHL/EHL) / AAV Gene Therapy / Concizumab-Marstacimab-Fitusiran PK-PD /
Inhibitor Immunology / Thrombin Generation / Bleeds / Arthropathy / QoL

$PARAM @annotated
// --- FIX PK (2-compartment, SHL: BeneFIX/Rixubis) ---
CL_FIX    : 0.30   : FIX clearance, SHL (dL/h per kg-normalized)
Vc_FIX    : 6.0    : FIX central volume, SHL (dL per 70 kg; recovery ~1 IU/dL per IU/kg)
Q_FIX     : 0.15   : FIX intercompartment CL (dL/h)
Vp_FIX    : 4.0    : FIX peripheral volume (dL)

// --- EHL-FIX variants (Fc-fusion / albumin-fusion / glycoPEGylated) ---
CL_FcFusion   : 0.075 : rFIX-Fc (Alprolix) clearance — t1/2 ~82h
Vc_FcFusion   : 6.0   : rFIX-Fc central volume
Q_FcFusion    : 0.05  : rFIX-Fc intercompartment CL
Vp_FcFusion   : 4.0   : rFIX-Fc peripheral volume

CL_AlbFusion  : 0.048 : rFIX-albumin (Idelvion) clearance — t1/2 ~102-104h
Vc_AlbFusion  : 5.8   : rFIX-albumin central volume
Q_AlbFusion   : 0.03  : rFIX-albumin intercompartment CL
Vp_AlbFusion  : 3.8   : rFIX-albumin peripheral volume

CL_GlycoPEG   : 0.065 : GlycoPEGylated rFIX (Refixia) clearance — t1/2 ~93h
Vc_GlycoPEG   : 6.0   : GlycoPEGylated rFIX central volume
Q_GlycoPEG    : 0.04  : GlycoPEGylated rFIX intercompartment CL
Vp_GlycoPEG   : 4.0   : GlycoPEGylated rFIX peripheral volume

// --- AAV gene therapy (etranacogene dezaparvovec-like, AAV5-FIX-Padua) ---
ka_AAV        : 0.35   : Vector clearance from blood into liver (1/h, rapid hepatic uptake)
k_transduce   : 1.00   : Hepatocyte transduction efficiency rate constant (1/h; calibrated so a standard dose transduces ~95% of achievable hepatocyte pool)
NAb_block     : 0.0    : Fractional block of transduction by pre-existing anti-AAV NAb (0-1)
k_expr_ramp   : 0.0015 : Transgene expression ramp-up rate (1/h; time constant ~28d, plateau by ~8-12 wk, HOPE-B-like)
FIXPadua_fold : 8.0    : FIX-Padua (R338L) specific activity fold-increase vs wild-type
Expr_plateau_max : 40.0 : Maximum steady-state endogenous FIX:C from AAV (IU/dL, HOPE-B-like)
k_capsid_immune : 0.006 : Rate of capsid-specific T-cell immune activation (1/h)
k_immune_decay  : 0.004 : Immune response resolution rate (1/h)
k_ALT_rise      : 0.08  : ALT rise rate constant per unit immune activation
k_ALT_fall      : 0.03  : ALT resolution rate constant (1/h)
ALT_thresh      : 1.5   : ALT fold-elevation threshold triggering steroid response
steroid_suppress : 0.75 : Fractional suppression of capsid immune response by corticosteroid
k_vector_dilution : 0.0000099 : Hepatocyte turnover -> episomal vector genome dilution (1/h; ~8-year half-life, gradual multi-year decline)
k_antigen_clear : 0.000825 : Intracellular capsid antigen clearance rate (1/h; ~5-week half-life; transient, distinct from durable episomal transgene)
k_immune_erosion : 0.0015 : Rate of immune-mediated transduced-hepatocyte loss per unit Capsid_Immune (1/h; active only during transient antigen-presentation window)

// --- Concizumab / Marstacimab PK (anti-TFPI mAb, SC) ---
ka_CONC   : 0.020  : Concizumab SC absorption rate (1/h)
CL_CONC   : 0.0021 : Concizumab clearance (L/h)
Vc_CONC   : 3.4    : Concizumab central volume (L)
F_CONC    : 0.65   : Concizumab SC bioavailability
EC50_CONC : 45.0   : Concizumab EC50 for TFPI neutralization (ng/mL)

ka_MARS   : 0.018  : Marstacimab SC absorption rate (1/h)
CL_MARS   : 0.0032 : Marstacimab clearance (L/h)
Vc_MARS   : 4.6    : Marstacimab central volume (L)
F_MARS    : 0.70   : Marstacimab SC bioavailability
EC50_MARS : 60.0   : Marstacimab EC50 for TFPI neutralization (ng/mL)

Emax_TFPI  : 0.90  : Maximum fractional TFPI neutralization achievable
FIXeq_TFPI : 18.0  : Max FIX-equivalent thrombin-generation boost from full TFPI blockade (IU/dL-equiv)
FIXeq_AT   : 12.0  : Max FIX-equivalent hemostatic boost from full AT knockdown (IU/dL-equiv, fitusiran)

// --- Fitusiran PK (SC) + AT mRNA/protein dynamics ---
ka_FITU   : 0.008  : Fitusiran SC absorption rate (1/h)
CL_FITU   : 0.022  : Fitusiran central clearance (L/h)
Vc_FITU   : 8.5    : Fitusiran central volume (L)
F_FITU    : 0.75   : Fitusiran SC bioavailability
ksyn_ATm  : 0.0065 : AT mRNA synthesis rate (relative/h)
kdeg_ATm  : 0.0065 : AT mRNA baseline degradation (1/h; t1/2 ~4.4 days)
Emax_FITU : 0.92   : Maximum AT mRNA knockdown by fitusiran
EC50_FITU : 0.008  : Fitusiran EC50 for AT mRNA knockdown (mg/L)
ksyn_ATp  : 0.0030 : AT protein synthesis rate (relative/h)
kdeg_ATp  : 0.0030 : AT protein baseline degradation (1/h; t1/2 ~9.6 days)

// --- Inhibitor immunology (anti-FIX) ---
k_inhibit : 0.00010 : Rate of inhibitor formation per IU/dL FIX exposure (BU/h per IU/dL)
k_inh_off : 0.0006  : Inhibitor spontaneous waning rate (1/h)
Ki_max    : 150.0   : Maximum inhibitor titer plateau (BU/mL)
IC50_inh  : 1.0     : Inhibitor IC50 on FIX activity (BU/mL)
NULL_MUT_MULT : 1.0 : Multiplier on inhibitor formation rate for null-mutation genotype

// --- Thrombin generation ---
ETP_base  : 100.0  : Baseline thrombin ETP (nmol*min, normalized to 100)
k_ETP_up  : 0.85   : Rate constant for ETP equilibration (1/h)
ETP_FIX_EC50 : 4.0 : FIX EC50 for ETP generation (IU/dL)
ETP_FIX_hill : 0.85 : Hill coefficient for FIX-ETP relationship
AT_inhibit_ETP : 0.3 : Fractional ETP modulation by AT level

// --- Bleed model ---
ABR_base  : 28.0   : Baseline ABR for untreated severe HB (bleeds/year)
FIX_ABR_EC50 : 3.5  : FIX activity EC50 for bleed reduction (IU/dL)
FIX_ABR_hill : 1.15 : Hill exponent for bleed reduction
ABR_floor : 0.4    : Residual bleed risk even at very high FIX levels

// --- Joint damage (arthropathy) ---
k_joint_in  : 0.0007 : Joint damage increment rate per bleed (score/bleed)
k_joint_rep : 0.0001 : Joint repair rate (1/h; very slow)
Joint_max   : 100.0  : Maximum Pettersson-like joint score
k_syno_in   : 0.002  : Synovitis increment rate from iron/ROS
k_syno_out  : 0.0005 : Synovitis resolution rate (1/h)

// --- Quality of life ---
k_QoL_joint : 0.004  : QoL decrement per joint score unit (per 100)
k_QoL_ABR   : 0.010  : QoL decrement per ABR unit (per bleed/year)
QoL_max     : 1.0    : Maximum QoL (Haem-A-QoL-derived utility = 1)

// --- Body weight (for dose scaling) ---
BW          : 70.0  : Body weight (kg)

// --- Simulation flags ---
FIX_TYPE    : 1     : 1=SHL, 2=Fc-fusion(EHL), 3=Albumin-fusion(EHL), 4=GlycoPEG(EHL)
USE_TFPI_mAb : 0    : 1 = concizumab/marstacimab active (choose via TFPI_DRUG)
TFPI_DRUG   : 1     : 1=concizumab, 2=marstacimab (only used if USE_TFPI_mAb=1)
USE_FITU    : 0     : 1 = fitusiran active
USE_AAV     : 0     : 1 = AAV gene therapy administered
USE_STEROID : 0     : 1 = corticosteroid taper active (for AAV transaminitis)
INHIBITOR_ON : 0    : 1 = inhibitor development active

$CMT @annotated
FIX_C     : FIX central compartment, SHL (IU/dL)
FIX_P     : FIX peripheral compartment, SHL
FIXe_C    : FIX central compartment, EHL variant (IU/dL)
FIXe_P    : FIX peripheral compartment, EHL variant
AAV_Vector : Circulating AAV vector genome (relative units)
Transduced_Hep : Transduced hepatocyte fraction (0-1)
Capsid_Antigen : Transient intracellular capsid antigen pool (relative, drives immune response only)
Transgene_Expr : Endogenous FIX-Padua transgene expression (IU/dL)
Capsid_Immune  : Capsid-specific immune activation (0-1)
ALT_level      : ALT fold-elevation (baseline=1)
CONC_SC   : Concizumab SC depot (mg)
CONC_C    : Concizumab central (ng/mL)
MARS_SC   : Marstacimab SC depot (mg)
MARS_C    : Marstacimab central (ng/mL)
FITU_SC   : Fitusiran SC depot (mg)
FITU_C    : Fitusiran central (mg/L)
AT_mRNA   : Antithrombin mRNA (relative, baseline=1)
AT_prot   : Antithrombin protein (relative, baseline=1)
Inhibitor : Anti-FIX inhibitor titer (BU/mL)
Thrombin_ETP : Thrombin generation potential (normalized 0-100)
CumBleeds    : Cumulative bleed count
JointScore   : Hemophilic arthropathy score (0-100)
Synovitis    : Synovial inflammation score (0-1)
QoL          : Quality of life (utility, 0-1)

$GLOBAL
double FIX_act;        // FIX activity accounting for inhibitor (replacement-derived)
double FIX_replace_raw; // Raw replacement FIX before inhibitor neutralization
double TFPI_effect;    // TFPI-mAb FIX-equivalent activity
double AT_effect;      // Antithrombin level effect on thrombin
double ABR_inst;       // Instantaneous bleed rate (bleeds/year)
double FIX_eff_total;  // Total effective FIX (replacement + transgene + TFPI-equiv)
double CL_eff, VC_eff, Q_eff, VP_eff; // Active EHL PK parameter set

$MAIN
// Select EHL-FIX PK parameter set based on FIX_TYPE flag (2/3/4)
if(FIX_TYPE == 2) {
  CL_eff = CL_FcFusion; VC_eff = Vc_FcFusion; Q_eff = Q_FcFusion; VP_eff = Vp_FcFusion;
} else if(FIX_TYPE == 3) {
  CL_eff = CL_AlbFusion; VC_eff = Vc_AlbFusion; Q_eff = Q_AlbFusion; VP_eff = Vp_AlbFusion;
} else if(FIX_TYPE == 4) {
  CL_eff = CL_GlycoPEG; VC_eff = Vc_GlycoPEG; Q_eff = Q_GlycoPEG; VP_eff = Vp_GlycoPEG;
} else {
  CL_eff = CL_FIX; VC_eff = Vc_FIX; Q_eff = Q_FIX; VP_eff = Vp_FIX;
}

if(NEWIND <= 1) {
  _init_AT_mRNA  = 1.0;
  _init_AT_prot  = 1.0;
  _init_QoL      = 0.78;   // Starting QoL for severe HB patient
  _init_JointScore = 6.0;  // Mild pre-existing joint damage
  _init_Thrombin_ETP = 18.0; // Severely reduced ETP in untreated HB
}

$ODE
// -------------------------------------------------------
// Replacement FIX pooled across SHL/EHL depots, inhibitor-neutralized
// -------------------------------------------------------
FIX_replace_raw = FIX_C + FIXe_C;
FIX_act = FIX_replace_raw / (1.0 + Inhibitor / IC50_inh);
FIX_act = (FIX_act < 0) ? 0 : FIX_act;

// -------------------------------------------------------
// TFPI-mAb (concizumab/marstacimab) FIX-equivalent activity
// -------------------------------------------------------
if(USE_TFPI_mAb == 1) {
  if(TFPI_DRUG == 1) {
    TFPI_effect = Emax_TFPI * CONC_C / (EC50_CONC + CONC_C);
  } else {
    TFPI_effect = Emax_TFPI * MARS_C / (EC50_MARS + MARS_C);
  }
  TFPI_effect = TFPI_effect * FIXeq_TFPI;
} else {
  TFPI_effect = 0.0;
}

// -------------------------------------------------------
// FIX SHL PK (2-compartment, IV bolus)
// -------------------------------------------------------
dxdt_FIX_C = -CL_FIX/Vc_FIX * FIX_C - Q_FIX/Vc_FIX * FIX_C + Q_FIX/Vp_FIX * FIX_P;
dxdt_FIX_P =  Q_FIX/Vc_FIX * FIX_C - Q_FIX/Vp_FIX * FIX_P;

// -------------------------------------------------------
// FIX EHL PK (2-compartment, IV bolus) — active variant selected by FIX_TYPE
// -------------------------------------------------------
dxdt_FIXe_C = -CL_eff/VC_eff * FIXe_C - Q_eff/VC_eff * FIXe_C + Q_eff/VP_eff * FIXe_P;
dxdt_FIXe_P =  Q_eff/VC_eff * FIXe_C - Q_eff/VP_eff * FIXe_P;

// -------------------------------------------------------
// AAV gene therapy: vector clearance -> hepatocyte transduction ->
// transgene expression ramp-up -> capsid immune response -> ALT rise ->
// steroid-modulated immune suppression -> long-term vector dilution
// -------------------------------------------------------
double transduce_rate = 0.0;
if(USE_AAV == 1) {
  transduce_rate = k_transduce * (1.0 - NAb_block);
}
double transduction_flux = transduce_rate * AAV_Vector/100.0 * (1.0 - Transduced_Hep);
dxdt_AAV_Vector = -ka_AAV * AAV_Vector;

// Capsid antigen is presented transiently within transduced hepatocytes (~5-week
// intracellular clearance) — distinct from the durable episomal transgene, so the
// adaptive immune response resolves even though FIX-Padua expression persists.
dxdt_Capsid_Antigen = transduction_flux - k_antigen_clear * Capsid_Antigen;

double immune_drive = k_capsid_immune * Capsid_Antigen;
if(USE_STEROID == 1) {
  immune_drive = immune_drive * (1.0 - steroid_suppress);
}
dxdt_Capsid_Immune = immune_drive * (1.0 - Capsid_Immune) - k_immune_decay * Capsid_Immune;
dxdt_ALT_level = k_ALT_rise * Capsid_Immune * (1.0 - ALT_level/5.0) - k_ALT_fall * (ALT_level - 1.0);

// Transduced hepatocyte pool: gains from transduction, slow loss from hepatocyte
// turnover (vector dilution) plus transient immune-mediated erosion while the
// capsid antigen window is active (self-limiting; long-term durability dominated
// by k_vector_dilution once Capsid_Immune resolves).
dxdt_Transduced_Hep = transduction_flux - k_vector_dilution * Transduced_Hep
                      - k_immune_erosion * Capsid_Immune * Transduced_Hep;

double Expr_target = Expr_plateau_max * Transduced_Hep;
dxdt_Transgene_Expr = k_expr_ramp * (Expr_target - Transgene_Expr);

// -------------------------------------------------------
// Concizumab SC PK (1-compartment)
// -------------------------------------------------------
dxdt_CONC_SC = -ka_CONC * CONC_SC;
dxdt_CONC_C  = F_CONC * ka_CONC * CONC_SC * 1000.0 / Vc_CONC - CL_CONC/Vc_CONC * CONC_C;

// -------------------------------------------------------
// Marstacimab SC PK (1-compartment)
// -------------------------------------------------------
dxdt_MARS_SC = -ka_MARS * MARS_SC;
dxdt_MARS_C  = F_MARS * ka_MARS * MARS_SC * 1000.0 / Vc_MARS - CL_MARS/Vc_MARS * MARS_C;

// -------------------------------------------------------
// Fitusiran PK (SC 1-compartment)
// -------------------------------------------------------
dxdt_FITU_SC = -ka_FITU * FITU_SC;
dxdt_FITU_C  = F_FITU * ka_FITU * FITU_SC / Vc_FITU - CL_FITU/Vc_FITU * FITU_C;

// -------------------------------------------------------
// AT mRNA/Protein knockdown by Fitusiran (indirect response)
// -------------------------------------------------------
double kd_ATm_total = kdeg_ATm;
if(USE_FITU == 1) {
  kd_ATm_total = kdeg_ATm * (1.0 + Emax_FITU * FITU_C / (EC50_FITU + FITU_C));
}
dxdt_AT_mRNA = ksyn_ATm - kd_ATm_total * AT_mRNA;
dxdt_AT_prot = ksyn_ATp * AT_mRNA - kdeg_ATp * AT_prot;
AT_effect = AT_prot; // 0-1 scale

// -------------------------------------------------------
// Inhibitor titer dynamics (anti-FIX IgG)
// -------------------------------------------------------
double inhibit_formation = 0.0;
if(INHIBITOR_ON == 1 && FIX_act > 0) {
  inhibit_formation = k_inhibit * NULL_MUT_MULT * FIX_act * (1.0 - Inhibitor/Ki_max);
}
dxdt_Inhibitor = inhibit_formation - k_inh_off * Inhibitor;

// -------------------------------------------------------
// Total effective FIX activity (replacement + transgene + TFPI-mAb-equivalent)
// -------------------------------------------------------
double AT_rebalance_effect = FIXeq_AT * (1.0 - AT_effect);
FIX_eff_total = FIX_act + Transgene_Expr + TFPI_effect + AT_rebalance_effect;

// -------------------------------------------------------
// Thrombin Generation Potential (ETP)
// -------------------------------------------------------
double ETP_FIX = ETP_base * pow(FIX_eff_total, ETP_FIX_hill) /
                 (pow(ETP_FIX_EC50, ETP_FIX_hill) + pow(FIX_eff_total, ETP_FIX_hill));
double ETP_ss = ETP_FIX * (1.0 + (1.0 - AT_effect) * 0.5);
dxdt_Thrombin_ETP = k_ETP_up * (ETP_ss - Thrombin_ETP);

// -------------------------------------------------------
// Bleed rate (FIX-activity-dependent, Hill inhibitory model)
// -------------------------------------------------------
double FIX_prot = pow(FIX_eff_total, FIX_ABR_hill) /
                   (pow(FIX_ABR_EC50, FIX_ABR_hill) + pow(FIX_eff_total, FIX_ABR_hill));
ABR_inst = ABR_base * (1.0 - FIX_prot) + ABR_floor;
ABR_inst = (ABR_inst < ABR_floor) ? ABR_floor : ABR_inst;
dxdt_CumBleeds = ABR_inst / 8760.0;

// -------------------------------------------------------
// Synovitis & hemophilic arthropathy
// -------------------------------------------------------
double bleed_per_h = ABR_inst / 8760.0;
dxdt_Synovitis = k_syno_in * bleed_per_h * (1.0 - Synovitis) - k_syno_out * Synovitis;
dxdt_JointScore = k_joint_in * bleed_per_h * (Joint_max - JointScore) -
                  k_joint_rep * JointScore * (1.0 - JointScore/Joint_max);

// -------------------------------------------------------
// Quality of Life
// -------------------------------------------------------
double QoL_target = QoL_max - k_QoL_joint * JointScore/100.0 - k_QoL_ABR * ABR_inst/28.0;
QoL_target = (QoL_target < 0.1) ? 0.1 : QoL_target;
dxdt_QoL = 0.01 * (QoL_target - QoL);

$TABLE
double FIX_activity   = FIX_act;
double FIX_total      = FIX_eff_total;
double Transgene_Level = Transgene_Expr;
double Concizumab_conc = CONC_C;
double Marstacimab_conc = MARS_C;
double Fitusiran_conc  = FITU_C;
double AT_level        = AT_prot;
double ETP             = Thrombin_ETP;
double BleedRate_annual = ABR_inst;
double Joint_damage    = JointScore;
double HRQoL           = QoL;
double Inhibitor_titer = Inhibitor;
double ALT_fold        = ALT_level;
double ImmuneActivation = Capsid_Immune;

$CAPTURE
FIX_activity FIX_total Transgene_Level Concizumab_conc Marstacimab_conc
Fitusiran_conc AT_level ETP BleedRate_annual Joint_damage HRQoL
Inhibitor_titer ALT_fold ImmuneActivation TFPI_effect
'

##############################################################################
# Compile model
##############################################################################
mod <- mcode("HemophiliaB_QSP", code)

##############################################################################
# Helper: dose events
##############################################################################

#' SHL/EHL FIX prophylaxis events (IV bolus into cmt 1 [SHL] or 3 [EHL])
#' recovery ~1 IU/dL per IU/kg for SHL/EHL FIX (lower than FVIII's ~2%)
fix_prophy <- function(duration_days = 365, BW = 70, dose_iukg = 40,
                        freq_days = c(0, 3.5), cmt = 1, recovery = 1.0) {
  dose_iudl <- dose_iukg * recovery
  n_cycles <- ceiling(duration_days / 7) + 1
  evs <- NULL
  for (w in 0:(n_cycles - 1)) {
    for (d in freq_days) {
      t <- w * 7 + d
      if (t <= duration_days) {
        evs <- rbind(evs, data.frame(time = t * 24, cmt = cmt,
                                      amt = dose_iudl, evid = 1, rate = -2))
      }
    }
  }
  as_data_frame(evs)
}

#' Simple fixed-interval IV bolus dosing (e.g. EHL-FIX Q7/Q10/Q14 days)
fix_prophy_interval <- function(duration_days = 365, BW = 70, dose_iukg = 50,
                                 interval_days = 7, cmt = 3, recovery = 1.0) {
  dose_iudl <- dose_iukg * recovery
  times <- seq(0, duration_days, by = interval_days)
  data.frame(time = times * 24, cmt = cmt, amt = dose_iudl, evid = 1, rate = -2)
}

#' Concizumab/Marstacimab SC dosing (loading + maintenance)
tfpi_dosing <- function(duration_days = 365, BW = 70, cmt = 11,
                         loading_mg = 210, maint_mg = 15, freq_days_maint = 1) {
  evs <- NULL
  evs <- rbind(evs, data.frame(time = 0, cmt = cmt, amt = loading_mg, evid = 1))
  n_doses <- floor((duration_days - 1) / freq_days_maint)
  for (i in 0:n_doses) {
    t <- (1 + i * freq_days_maint) * 24
    if (t / 24 <= duration_days) {
      evs <- rbind(evs, data.frame(time = t, cmt = cmt, amt = maint_mg, evid = 1))
    }
  }
  as_data_frame(evs)
}

#' Fitusiran SC monthly dosing (50 mg/month)
fitu_dosing <- function(duration_days = 365, dose_mg = 50, freq_days = 28) {
  evs <- NULL
  n_doses <- floor(duration_days / freq_days) + 1
  for (i in 0:(n_doses - 1)) {
    t <- i * freq_days * 24
    if (t / 24 <= duration_days) {
      evs <- rbind(evs, data.frame(time = t, cmt = 15, amt = dose_mg, evid = 1))
    }
  }
  as_data_frame(evs)
}

#' AAV gene therapy single IV dose (relative vector genome units)
aav_dosing <- function(vector_dose = 100) {
  data.frame(time = 0, cmt = 5, amt = vector_dose, evid = 1)
}

##############################################################################
# Scenario definitions
##############################################################################

run_scenario <- function(scenario, duration_days = 365, BW = 70) {

  base_params <- list(BW = BW, INHIBITOR_ON = 0, FIX_TYPE = 1,
                       USE_TFPI_mAb = 0, TFPI_DRUG = 1, USE_FITU = 0,
                       USE_AAV = 0, USE_STEROID = 0)

  # Scenario 1: No prophylaxis (on-demand only)
  if (scenario == 1) {
    params <- base_params
    ev <- NULL

  # Scenario 2: SHL-rFIX prophylaxis 2x/week (40 IU/kg)
  } else if (scenario == 2) {
    params <- c(base_params, FIX_TYPE = 1)
    ev <- fix_prophy(duration_days, BW, dose_iukg = 40, freq_days = c(0, 3.5),
                      cmt = 1, recovery = 1.0)

  # Scenario 3: EHL rFIX-Fc (Alprolix-like) Q7-10 days (50 IU/kg)
  } else if (scenario == 3) {
    params <- c(base_params, FIX_TYPE = 2)
    ev <- fix_prophy_interval(duration_days, BW, dose_iukg = 50,
                               interval_days = 7, cmt = 3, recovery = 1.0)

  # Scenario 4: EHL rFIX-albumin (Idelvion-like) Q14 days (75 IU/kg)
  } else if (scenario == 4) {
    params <- c(base_params, FIX_TYPE = 3)
    ev <- fix_prophy_interval(duration_days, BW, dose_iukg = 75,
                               interval_days = 14, cmt = 3, recovery = 1.0)

  # Scenario 5: GlycoPEGylated rFIX (Refixia-like) Q7 days (40 IU/kg)
  } else if (scenario == 5) {
    params <- c(base_params, FIX_TYPE = 4)
    ev <- fix_prophy_interval(duration_days, BW, dose_iukg = 40,
                               interval_days = 7, cmt = 3, recovery = 1.0)

  # Scenario 6: Concizumab SC daily prophylaxis
  } else if (scenario == 6) {
    params <- c(base_params, USE_TFPI_mAb = 1, TFPI_DRUG = 1)
    ev <- tfpi_dosing(duration_days, BW, cmt = 11, loading_mg = 210,
                       maint_mg = 15, freq_days_maint = 1)

  # Scenario 7: Marstacimab SC weekly prophylaxis
  } else if (scenario == 7) {
    params <- c(base_params, USE_TFPI_mAb = 1, TFPI_DRUG = 2)
    ev <- tfpi_dosing(duration_days, BW, cmt = 13, loading_mg = 300,
                       maint_mg = 150, freq_days_maint = 7)

  # Scenario 8: Fitusiran SC monthly prophylaxis
  } else if (scenario == 8) {
    params <- c(base_params, USE_FITU = 1)
    ev <- fitu_dosing(duration_days, dose_mg = 50, freq_days = 28)

  # Scenario 9: AAV gene therapy single dose + steroid taper
  } else if (scenario == 9) {
    params <- c(base_params, USE_AAV = 1, USE_STEROID = 1)
    ev <- aav_dosing(vector_dose = 100)

  # Scenario 10: Inhibitor-positive patient, ITI + on-demand bypass management
  } else if (scenario == 10) {
    params <- c(base_params, INHIBITOR_ON = 1, NULL_MUT_MULT = 3.0, FIX_TYPE = 1)
    ev <- fix_prophy(duration_days, BW, dose_iukg = 50, freq_days = c(0, 2, 4),
                      cmt = 1, recovery = 1.0)
  }

  params_mod <- do.call(param, c(list(mod), params))

  if (is.null(ev) || nrow(ev) == 0) {
    out <- mrgsim(params_mod, end = duration_days * 24, delta = 1)
  } else {
    ev_obj <- as.ev(ev)
    out <- mrgsim(params_mod, events = ev_obj, end = duration_days * 24, delta = 1)
  }

  as_tibble(out) %>%
    mutate(scenario = scenario,
           scenario_label = c("1" = "No Prophylaxis",
                               "2" = "SHL-rFIX 2x/wk",
                               "3" = "EHL-FIX-Fc Q7-10d",
                               "4" = "EHL-FIX-Alb Q14d",
                               "5" = "GlycoPEG-FIX Q7d",
                               "6" = "Concizumab Daily",
                               "7" = "Marstacimab Q1W",
                               "8" = "Fitusiran Q1M",
                               "9" = "AAV Gene Therapy",
                               "10" = "Inhibitor + ITI/Bypass")[as.character(scenario)],
           time_days = time / 24)
}

##############################################################################
# Run all scenarios
##############################################################################
message("Running Hemophilia B QSP simulations ...")

scenarios_out <- lapply(1:10, function(s) {
  message("  Scenario ", s, " ...")
  run_scenario(s, duration_days = 365 * 2)
})

all_out <- bind_rows(scenarios_out)

##############################################################################
# Plot 1: FIX activity — replacement regimens (first 28 days)
##############################################################################
p1 <- all_out %>%
  filter(scenario %in% 2:5, time_days <= 28) %>%
  ggplot(aes(x = time_days, y = FIX_activity, color = scenario_label)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "orange", alpha = 0.7) +
  scale_y_log10() +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "FIX Activity — SHL vs EHL Prophylaxis (First 28 Days)",
       x = "Time (days)", y = "FIX Activity (IU/dL, log scale)", color = "Scenario") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

##############################################################################
# Plot 2: Non-factor rebalancing agents concentration/effect
##############################################################################
p2 <- all_out %>%
  filter(scenario %in% c(6, 7, 8), time_days <= 180) %>%
  ggplot(aes(x = time_days, y = FIX_total, color = scenario_label)) +
  geom_line(linewidth = 0.9) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Effective Hemostatic FIX-Equivalent — Non-Factor Rebalancing Agents",
       x = "Time (days)", y = "FIX-Equivalent Activity (IU/dL)", color = "Scenario") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

##############################################################################
# Plot 3: AAV gene therapy trajectory
##############################################################################
p3 <- all_out %>%
  filter(scenario == 9, time_days <= 365) %>%
  ggplot(aes(x = time_days)) +
  geom_line(aes(y = Transgene_Level, color = "Transgene FIX Expression"), linewidth = 1) +
  geom_line(aes(y = ALT_fold * 10, color = "ALT (x10 fold)"), linewidth = 0.8, linetype = "dashed") +
  scale_color_manual(values = c("Transgene FIX Expression" = "#2b9348", "ALT (x10 fold)" = "#d90429")) +
  labs(title = "AAV Gene Therapy: Endogenous FIX Expression & Transaminitis",
       x = "Time (days)", y = "Level", color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

##############################################################################
# Plot 4: Annual bleed rate comparison across all scenarios
##############################################################################
p4 <- all_out %>%
  filter(time_days >= 60, time_days <= 365) %>%
  group_by(scenario_label) %>%
  summarise(ABR_mean = mean(BleedRate_annual), ABR_sd = sd(BleedRate_annual), .groups = "drop") %>%
  mutate(scenario_label = factor(scenario_label, levels = c(
    "No Prophylaxis", "SHL-rFIX 2x/wk", "EHL-FIX-Fc Q7-10d", "EHL-FIX-Alb Q14d",
    "GlycoPEG-FIX Q7d", "Concizumab Daily", "Marstacimab Q1W", "Fitusiran Q1M",
    "AAV Gene Therapy", "Inhibitor + ITI/Bypass"))) %>%
  ggplot(aes(x = scenario_label, y = ABR_mean, fill = scenario_label)) +
  geom_col(alpha = 0.85) +
  geom_errorbar(aes(ymin = pmax(0, ABR_mean - ABR_sd), ymax = ABR_mean + ABR_sd), width = 0.3) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "darkred") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Simulated Annual Bleed Rate (ABR) by Treatment Scenario",
       x = NULL, y = "ABR (bleeds/year)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")

##############################################################################
# Plot 5: Joint score & QoL over 2 years (prophylaxis vs none)
##############################################################################
p5 <- all_out %>%
  filter(scenario %in% c(1, 3, 9)) %>%
  ggplot(aes(x = time_days / 365, y = JointScore, color = scenario_label)) +
  geom_line(linewidth = 0.9) +
  labs(title = "Hemophilic Arthropathy Progression (2-Year Horizon)",
       x = "Time (years)", y = "Joint Damage Score (0-100)", color = "Scenario") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

combined <- (p1 + p2) / (p3 + p4) / p5
ggsave("hb_qsp_simulation_plots.png", combined, width = 14, height = 16, dpi = 150)

message("Hemophilia B QSP simulation complete. Plots saved to hb_qsp_simulation_plots.png")
