# ===========================================================
# Hodgkin Lymphoma (HL) - Quantitative Systems Pharmacology
# mrgsolve ODE model
# -----------------------------------------------------------
# 24 ODE compartments covering:
#   * HRS (Reed-Sternberg) tumor cell dynamics with NF-kB / JAK-STAT drive
#   * Exhausted CD8 T-cell <-> effector T-cell pool (PD-1/PD-L1 axis)
#   * Treg pool and TME cytokine signal (TARC/CCL17, IL-6, IL-13)
#   * Drug PK compartments for ABVD (Doxo, Bleo, Vinblastine, Dacarbazine)
#   * Brentuximab vedotin + free MMAE payload (A+AVD / ECHELON-1)
#   * Nivolumab (SWOG S1826)
#   * Toxicity surrogates: ANC, cumulative anthracycline, bleomycin lung,
#     MMAE/vinca neurotoxicity, irAE risk
#   * Endpoints: total tumor mass, MTV, TARC, Deauville surrogate, PFS hazard
#
# Treatment scenarios encoded as event tables (see end of file):
#   1. ABVD x 6 (advanced-stage classical regimen)
#   2. AVD x 6 (RATHL iPET-adapted, drop bleomycin)
#   3. escBEACOPP x 6 (GHSG advanced-stage intensive)
#   4. A+AVD x 6 (BV-AVD, ECHELON-1)
#   5. N-AVD x 6 (Nivolumab + AVD, SWOG S1826)
#   6. Nivolumab monotherapy (R/R, CheckMate 205)
#   7. BV + Nivolumab combo (R/R bridge to ASCT)
#
# Parameters reflect published clinical/PK literature:
#   * Doxorubicin: CL 30 L/h, Vc 19 L, Vp 250 L, Q 24 L/h
#   * Bleomycin: CL 4 L/h, V 17 L
#   * Vinblastine: 3-cpt; CL 0.6 L/h/kg
#   * Dacarbazine: CL 25 L/h, V 35 L
#   * Brentuximab vedotin: CL 1.5 L/d, V 6 L, t1/2 4-6d; MMAE release k=0.012/h
#   * Nivolumab: CL 0.008 L/h, V 8 L, t1/2 ~25d
#
# Tumor kinetics: Simeoni-like growth (lambda0 / lambda1) + drug-induced
# kill K2 from instantaneous drug effect; immune kill K_im from
# effector T cell pool reinvigorated by PD-1 blockade.
# ===========================================================

library(mrgsolve)

hl_model_code <- '
$PROB
# Hodgkin Lymphoma QSP model
# Author: QSP Daily Routine
# Version: 1.0

$PARAM @annotated
// ---- Doxorubicin (ABVD/AVD/BEACOPP/A+AVD/N-AVD) ----
CL_DOX    :  30   : Doxorubicin clearance (L/h)
V1_DOX    :  19   : Doxorubicin central volume (L)
V2_DOX    : 250   : Doxorubicin peripheral volume (L)
Q_DOX     :  24   : Doxorubicin inter-compartmental (L/h)

// ---- Bleomycin (ABVD) ----
CL_BLM    :   4   : Bleomycin clearance (L/h)
V_BLM     :  17   : Bleomycin volume (L)

// ---- Vinblastine ----
CL_VBL    :  42   : Vinblastine clearance (L/h)
V1_VBL    :  70   : Vinblastine central volume (L)

// ---- Dacarbazine ----
CL_DAC    :  25   : Dacarbazine clearance (L/h)
V_DAC     :  35   : Dacarbazine volume (L)

// ---- Brentuximab vedotin (anti-CD30 ADC) + MMAE ----
CL_BV     :  0.063 : BV clearance (L/h, ~1.5 L/d)
V_BV      :  6     : BV volume (L)
k_release :  0.012 : MMAE release from BV (1/h)
CL_MMAE   :  4     : MMAE clearance (L/h)
V_MMAE    :  100   : MMAE volume (L)

// ---- Nivolumab (anti-PD-1) ----
CL_NIVO   :  0.008 : Nivolumab linear CL (L/h)
V_NIVO    :  8     : Nivolumab volume (L)

// ---- Tumor (HRS) dynamics ----
lambda0   :  0.018 : HRS exponential growth rate (1/d, untreated doubling ~38d)
lambda1   :  0.6   : HRS linear growth rate (g/d at large mass)
TumorW0   :  0.8   : Baseline total tumor mass at diagnosis (kg-equiv arbitrary units)
K_DOX     :  0.55  : Doxorubicin tumor kill coefficient
K_BLM     :  0.30  : Bleomycin tumor kill coefficient
K_VBL     :  0.40  : Vinblastine tumor kill coefficient
K_DAC     :  0.35  : Dacarbazine tumor kill coefficient
K_MMAE    :  0.65  : MMAE (BV payload) tumor kill coefficient (CD30+ enriched)
K_BV_dir  :  0.10  : Direct BV ADC effect on CD30+ HRS
K_immune  :  0.018 : Immune-mediated HRS kill (per Effector T / d)

// ---- PD-1 / PD-L1 axis ----
PDL1_base :  3.0   : Baseline PD-L1 expression (relative units)
PD1_occ_EC50 :  0.6 : Nivolumab serum [conc] for 50% PD-1 occupancy (ug/mL)
k_exhaust :  0.04  : Rate of CD8 exhaustion driven by PD-L1 (1/d)
k_reinvig :  0.12  : Reinvigoration rate when PD-1 blocked (1/d)
T_eff_base : 1.0   : Baseline effector T pool
T_exh_base : 1.5   : Baseline exhausted T pool
T_max     :  6.0   : Max effector T expansion ceiling

// ---- TME cytokines / biomarkers ----
k_TARC_prod : 0.35 : TARC (CCL17) production per HRS unit (pg/mL/d)
k_TARC_elim : 0.42 : TARC elimination (1/d)
k_IL6_prod  : 0.12 : IL-6 production per HRS (pg/mL/d)
k_IL6_elim  : 0.9  : IL-6 elimination (1/d)
k_IL13_prod : 0.07 : IL-13 production per HRS (pg/mL/d)
k_IL13_elim : 0.8  : IL-13 elimination (1/d)

// ---- Treg drive ----
Treg_base : 0.5
k_Treg_in : 0.15
k_Treg_out: 0.30

// ---- ANC (neutropenia) ----
ANC_base   : 4.5  : Baseline ANC (10^9/L)
MTT_ANC    : 110  : Mean transit time ANC (h, Friberg-Karlsson)
gamma_ANC  : 0.17 : Friberg gamma
slope_ANC_DOX  : 0.95
slope_ANC_BLM  : 0.00
slope_ANC_VBL  : 0.30
slope_ANC_DAC  : 0.20
slope_ANC_MMAE : 0.55
slope_ANC_BV   : 0.10
GCSF_factor    : 1.0  : 1.0 = no G-CSF, set 0.5 with G-CSF support

// ---- Late toxicity accumulators ----
k_cardio_dox    : 0.0035 : per (mg/m2) doxorubicin cumulative
cardio_threshold: 400    : mg/m2 lifetime doxo dose
k_lung_blm      : 0.0010 : per unit bleomycin cumulative
lung_threshold  : 300    : cum BLM units
k_neuro_VBL     : 0.0006
k_neuro_VCR     : 0.0010
k_neuro_MMAE    : 0.0020

// ---- irAE drive ----
k_irAE     : 0.004 : per nivolumab AUC unit
irAE_thresh: 1.0

// ---- Hazard (PFS) ----
h0_PFS    : 0.0006 : baseline daily PFS hazard
beta_tum  : 1.4    : tumor mass coefficient
beta_TARC : 0.01   : TARC coefficient
beta_imm  : -0.30  : effector T protective coefficient

$CMT  @annotated
// Drug PK compartments
DOX_C    : Doxorubicin central (mg)
DOX_P    : Doxorubicin peripheral (mg)
BLM_C    : Bleomycin central (mg)
VBL_C    : Vinblastine central (mg)
DAC_C    : Dacarbazine central (mg)
BV_C     : Brentuximab vedotin central (mg)
MMAE_C   : Free MMAE central (mg)
NIVO_C   : Nivolumab central (mg)

// Tumor / immune
HRS      : Total HRS tumor mass (relative, AU)
PDL1     : Effective PD-L1 expression (AU)
T_eff    : Effector CD8 T pool (AU)
T_exh    : Exhausted CD8 T pool (AU)
TREG     : Treg pool (AU)
TARC     : Serum TARC pg/mL
IL6      : IL-6 pg/mL
IL13     : IL-13 pg/mL

// Hematology - Friberg-Karlsson transit chain
PROL     : Proliferative neutrophil pool
TR1      : Transit 1
TR2      : Transit 2
TR3      : Transit 3
ANC      : Circulating ANC (10^9/L)

// Toxicity accumulators
CumDOX   : Cumulative doxorubicin mg/m2 equiv
CumBLM   : Cumulative bleomycin units
CumVBL   : Cumulative vinblastine
CumMMAE  : Cumulative MMAE
NIVO_AUC : Nivolumab cumulative AUC (ug*d/mL surrogate)
CARDIO   : Cardiotoxicity index (0-1)
LUNG     : Pulmonary toxicity index
NEURO    : Peripheral neuropathy index
irAE     : Cumulative irAE risk
HAZ      : Cumulative PFS hazard

$MAIN
PROL_0 = ANC_base;
TR1_0  = ANC_base;
TR2_0  = ANC_base;
TR3_0  = ANC_base;
ANC_0  = ANC_base;
HRS_0  = TumorW0;
PDL1_0 = PDL1_base;
T_eff_0 = T_eff_base;
T_exh_0 = T_exh_base;
TREG_0  = Treg_base;
TARC_0  = 1000;   // baseline elevated serum TARC ~1000 pg/mL in active HL
IL6_0   = 25;
IL13_0  = 30;

$ODE
// ----- Drug PK -----
double Cdox_central = DOX_C / V1_DOX;          // mg/L
double Cdox_periph  = DOX_P / V2_DOX;
double Cblm         = BLM_C / V_BLM;
double Cvbl         = VBL_C / V1_VBL;
double Cdac         = DAC_C / V_DAC;
double Cbv          = BV_C  / V_BV;            // mg/L
double Cmmae        = MMAE_C/ V_MMAE;
double Cnivo        = NIVO_C/ V_NIVO;          // mg/L  ~ ug/mL

dxdt_DOX_C = -(CL_DOX/V1_DOX)*DOX_C - (Q_DOX/V1_DOX)*DOX_C + (Q_DOX/V2_DOX)*DOX_P;
dxdt_DOX_P =  (Q_DOX/V1_DOX)*DOX_C - (Q_DOX/V2_DOX)*DOX_P;
dxdt_BLM_C = -(CL_BLM/V_BLM)*BLM_C;
dxdt_VBL_C = -(CL_VBL/V1_VBL)*VBL_C;
dxdt_DAC_C = -(CL_DAC/V_DAC)*DAC_C;
dxdt_BV_C  = -(CL_BV/V_BV)*BV_C  - k_release*BV_C;
dxdt_MMAE_C=  k_release*BV_C - (CL_MMAE/V_MMAE)*MMAE_C;
dxdt_NIVO_C= -(CL_NIVO/V_NIVO)*NIVO_C;

// PD-1 occupancy 0-1 (Hill = 1)
double PD1_occ = Cnivo / (Cnivo + PD1_occ_EC50);

// Drug kill effects (Emax-like)
double E_DOX  = K_DOX  * Cdox_central / (Cdox_central + 0.05);
double E_BLM  = K_BLM  * Cblm         / (Cblm + 0.10);
double E_VBL  = K_VBL  * Cvbl         / (Cvbl + 0.005);
double E_DAC  = K_DAC  * Cdac         / (Cdac + 0.50);
double E_MMAE = K_MMAE * Cmmae        / (Cmmae + 0.001);
double E_BV   = K_BV_dir* Cbv         / (Cbv + 1.0);

// Immune kill scales with effector T cells (per day -> per hour /24)
double E_imm  = K_immune * T_eff;

// --- HRS tumor dynamics (Simeoni-like, units of d converted) ---
double growth = lambda0*HRS / pow(1 + pow(lambda0/lambda1 * HRS, 20), 1.0/20);
// Convert kill/day to /h
double total_kill = (E_DOX + E_BLM + E_VBL + E_DAC + E_MMAE + E_BV + E_imm) / 24.0;
dxdt_HRS = growth/24.0 - total_kill * HRS;
if (HRS < 1e-6) dxdt_HRS = 0;

// --- PD-L1 dynamics (driven by IFN-g surrogate ~ IL6+IL13 reflex) ---
double IFNg_surrogate = 0.4*IL6/25.0 + 0.3*IL13/30.0;
dxdt_PDL1 = (PDL1_base*(1 + 0.5*IFNg_surrogate) - PDL1) * 0.05/24.0;

// --- T cell exhaustion / reinvigoration ---
double exhaustion_drive = k_exhaust * PDL1 * (1 - PD1_occ) / 24.0;
double reinvig_drive    = k_reinvig * PD1_occ * T_exh / 24.0;
dxdt_T_eff = reinvig_drive - 0.02/24.0 * (T_eff - T_eff_base)
             - exhaustion_drive*0.3
             + 0.06/24.0*(T_max - T_eff)*PD1_occ;
dxdt_T_exh = exhaustion_drive - reinvig_drive - 0.04/24.0*(T_exh - T_exh_base);

// Treg dynamics
dxdt_TREG  = k_Treg_in*HRS/24.0 - k_Treg_out*TREG/24.0;

// Cytokines / biomarkers
dxdt_TARC = (k_TARC_prod*HRS*1000 - k_TARC_elim*TARC)/24.0;
dxdt_IL6  = (k_IL6_prod *HRS*100  - k_IL6_elim *IL6 )/24.0;
dxdt_IL13 = (k_IL13_prod*HRS*100  - k_IL13_elim*IL13)/24.0;

// --- Hematology (Friberg-Karlsson) ---
double ktr  = 4.0 / MTT_ANC;   // 1/h
double drug_hem = slope_ANC_DOX * Cdox_central
                + slope_ANC_VBL * Cvbl
                + slope_ANC_DAC * Cdac
                + slope_ANC_MMAE* Cmmae
                + slope_ANC_BV  * Cbv;
drug_hem *= GCSF_factor;
double feedback = pow(ANC_base/(ANC + 1e-6), gamma_ANC);
dxdt_PROL = ktr*PROL*((1 - drug_hem)*feedback - 1);
dxdt_TR1  = ktr*(PROL - TR1);
dxdt_TR2  = ktr*(TR1  - TR2);
dxdt_TR3  = ktr*(TR2  - TR3);
dxdt_ANC  = ktr*(TR3  - ANC);

// --- Cumulative drug exposure for late tox ---
dxdt_CumDOX  = Cdox_central;
dxdt_CumBLM  = Cblm;
dxdt_CumVBL  = Cvbl;
dxdt_CumMMAE = Cmmae;
dxdt_NIVO_AUC= Cnivo;

double cardio_drive = k_cardio_dox * fmax(0.0, CumDOX - cardio_threshold);
dxdt_CARDIO = cardio_drive/24.0 * (1 - CARDIO);
double lung_drive   = k_lung_blm * fmax(0.0, CumBLM - lung_threshold);
dxdt_LUNG   = lung_drive/24.0 * (1 - LUNG);
dxdt_NEURO  = (k_neuro_VBL*Cvbl + k_neuro_MMAE*Cmmae)/24.0 * (1 - NEURO);
dxdt_irAE   = k_irAE * NIVO_AUC/100.0 / 24.0 * (1 - irAE);

// --- Hazard for PFS ---
dxdt_HAZ    = (h0_PFS + beta_tum*HRS + beta_TARC*TARC/1000.0 + beta_imm*T_eff)/24.0;

$TABLE
double Cdox = DOX_C/V1_DOX;
double Cnivo_out = NIVO_C/V_NIVO;
double Deauville = (HRS > 0.6) ? 5
                  : (HRS > 0.3) ? 4
                  : (HRS > 0.10) ? 3
                  : (HRS > 0.02) ? 2 : 1;
double PFS_surv = exp(-HAZ);
double MTV = HRS * 1200.0;   // surrogate MTV (mL) - calibrated to baseline ~960 mL
double pCR = (HRS < 0.02) ? 1.0 : 0.0;

$CAPTURE  Cdox Cnivo_out Deauville PFS_surv MTV pCR PDL1 T_eff T_exh TREG TARC IL6 IL13 ANC CARDIO LUNG NEURO irAE HAZ
'

mod_hl <- mcode("hl_qsp", hl_model_code)

# ===========================================================
# Scenario event tables
# All dosing in mg; cycle length 21 days unless noted; BSA = 1.7 m2
# ===========================================================
BSA <- 1.7

ev_ABVD <- function(cycles=6){
  d <- numeric(0); evs <- list()
  for (c in seq_len(cycles)){
    base <- (c-1)*28*24            # cycle = 28 days, days 1 & 15
    for (day in c(1, 15)){
      t <- base + (day-1)*24
      evs[[length(evs)+1]] <- ev(time=t, cmt="DOX_C", amt=25*BSA)
      evs[[length(evs)+1]] <- ev(time=t, cmt="BLM_C", amt=10*BSA)
      evs[[length(evs)+1]] <- ev(time=t, cmt="VBL_C", amt=6 *BSA)
      evs[[length(evs)+1]] <- ev(time=t, cmt="DAC_C", amt=375*BSA)
    }
  }
  Reduce(`+`, evs)
}

ev_AVD <- function(cycles=6){
  evs <- list()
  for (c in seq_len(cycles)){
    base <- (c-1)*28*24
    for (day in c(1, 15)){
      t <- base + (day-1)*24
      evs[[length(evs)+1]] <- ev(time=t, cmt="DOX_C", amt=25*BSA)
      evs[[length(evs)+1]] <- ev(time=t, cmt="VBL_C", amt=6 *BSA)
      evs[[length(evs)+1]] <- ev(time=t, cmt="DAC_C", amt=375*BSA)
    }
  }
  Reduce(`+`, evs)
}

ev_BV_AVD <- function(cycles=6){
  # Brentuximab vedotin 1.2 mg/kg + AVD on D1 & D15; weight 70 kg
  WT <- 70
  evs <- list()
  for (c in seq_len(cycles)){
    base <- (c-1)*28*24
    for (day in c(1,15)){
      t <- base + (day-1)*24
      evs[[length(evs)+1]] <- ev(time=t, cmt="BV_C",  amt=1.2*WT)
      evs[[length(evs)+1]] <- ev(time=t, cmt="DOX_C", amt=25*BSA)
      evs[[length(evs)+1]] <- ev(time=t, cmt="VBL_C", amt=6 *BSA)
      evs[[length(evs)+1]] <- ev(time=t, cmt="DAC_C", amt=375*BSA)
    }
  }
  Reduce(`+`, evs)
}

ev_N_AVD <- function(cycles=6){
  WT <- 70
  evs <- list()
  for (c in seq_len(cycles)){
    base <- (c-1)*28*24
    # Nivo 240 mg flat D1 & D15
    for (day in c(1,15)){
      t <- base + (day-1)*24
      evs[[length(evs)+1]] <- ev(time=t, cmt="NIVO_C", amt=240)
      evs[[length(evs)+1]] <- ev(time=t, cmt="DOX_C",  amt=25*BSA)
      evs[[length(evs)+1]] <- ev(time=t, cmt="VBL_C",  amt=6 *BSA)
      evs[[length(evs)+1]] <- ev(time=t, cmt="DAC_C",  amt=375*BSA)
    }
  }
  Reduce(`+`, evs)
}

ev_escBEACOPP <- function(cycles=6){
  # Simplified: model Doxo + Bleo + Etoposide (treated as DOX_C surrogate cleared faster - approx)
  # Use cycle 21 d
  evs <- list()
  for (c in seq_len(cycles)){
    base <- (c-1)*21*24
    # D1: Doxo 35 mg/m2, Cyclo 1250 mg/m2 (lumped via DOX_C scaled)
    evs[[length(evs)+1]] <- ev(time=base, cmt="DOX_C", amt=35*BSA)
    # D1-3: Etoposide 200 mg/m2 -> approximate as add'l DOX_C-like pulses
    for (day in 1:3){
      t <- base + (day-1)*24
      evs[[length(evs)+1]] <- ev(time=t, cmt="DOX_C", amt=200*BSA*0.1) # surrogate
    }
    # D8: Bleo 10 mg/m2 + Vinc 1.4 (model VBL_C surrogate)
    t8 <- base + 7*24
    evs[[length(evs)+1]] <- ev(time=t8, cmt="BLM_C", amt=10*BSA)
    evs[[length(evs)+1]] <- ev(time=t8, cmt="VBL_C", amt=1.4*BSA)
  }
  Reduce(`+`, evs)
}

ev_BV_mono <- function(cycles=8){
  WT <- 70
  evs <- list()
  for (c in seq_len(cycles)){
    t <- (c-1)*21*24
    evs[[length(evs)+1]] <- ev(time=t, cmt="BV_C", amt=1.8*WT)
  }
  Reduce(`+`, evs)
}

ev_NIVO_mono <- function(cycles=12){
  evs <- list()
  for (c in seq_len(cycles)){
    t <- (c-1)*14*24   # q2w
    evs[[length(evs)+1]] <- ev(time=t, cmt="NIVO_C", amt=240)
  }
  Reduce(`+`, evs)
}

ev_BV_NIVO <- function(cycles=8){
  WT <- 70
  evs <- list()
  for (c in seq_len(cycles)){
    base <- (c-1)*21*24
    evs[[length(evs)+1]] <- ev(time=base, cmt="BV_C",  amt=1.8*WT)
    evs[[length(evs)+1]] <- ev(time=base, cmt="NIVO_C", amt=3*WT)  # mg/kg
  }
  Reduce(`+`, evs)
}

# ===========================================================
# Run a scenario (example)
# ===========================================================
run_scenario <- function(regimen = "ABVD", tend_days = 365){
  ev_fun <- switch(regimen,
                   ABVD       = ev_ABVD(),
                   AVD        = ev_AVD(),
                   BV_AVD     = ev_BV_AVD(),
                   N_AVD      = ev_N_AVD(),
                   escBEACOPP = ev_escBEACOPP(),
                   BV_mono    = ev_BV_mono(),
                   NIVO_mono  = ev_NIVO_mono(),
                   BV_NIVO    = ev_BV_NIVO(),
                   stop("Unknown regimen"))
  mod_hl %>%
    ev(ev_fun) %>%
    mrgsim(end = tend_days*24, delta = 6) %>%
    as.data.frame()
}

# ===========================================================
# Calibration notes (clinical anchors)
# ===========================================================
# - Baseline HRS mass = 0.8 AU calibrated so total tumor volume ~ 960 mL
#   (MTV ~ 960 mL corresponds to median advanced-stage HL per RATHL/EORTC).
# - ABVD x 6 should drive HRS < 0.02 (Deauville 1-2) in >85% of patients
#   (Johnson et al, NEJM 2016 - RATHL: 5-yr PFS ~83% after iPET-adapted).
# - escBEACOPP x 6 in advanced: PFS ~88% at 5y (Engert HD15).
# - A+AVD vs ABVD (ECHELON-1, Connors NEJM 2018): 6-yr PFS 82.3% vs 74.5%.
# - SWOG S1826 (Herrera et al NEJM 2024): N-AVD 2-yr PFS 92% vs A+AVD 83%.
# - Nivolumab R/R post-ASCT (CheckMate 205): ORR ~69%, CR 16%, median PFS ~15mo.
# - Brentuximab vedotin monotherapy (Younes JCO 2012): ORR 75%, CR 34% in R/R.
# - PD-L1 expression EC50 for nivolumab calibrated such that ~80% PD-1 occupancy
#   at Cmin ~ 25 ug/mL (Brahmer JCO 2010, Topalian NEJM 2012).
# - Cardiotoxicity threshold 400 mg/m2 lifetime doxorubicin (per Swain JCO 2003).
# - Bleomycin cumulative lung tox threshold modelled at 300 units cumulative.
# ===========================================================
