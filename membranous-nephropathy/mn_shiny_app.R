# =============================================================================
# Membranous Nephropathy (MN) QSP Model - Shiny Application
# =============================================================================
# Disease: Primary Membranous Nephropathy (autoimmune glomerular disease)
# Key drugs: Rituximab, Tacrolimus, Cyclophosphamide, ACE inhibitors
# Key biomarkers: anti-PLA2R1, proteinuria, eGFR, serum albumin,
#                 CD20+ B cells, complement (sC5b-9)
# =============================================================================

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(DT)

# Try to load mrgsolve; fall back to deSolve-based ODE solver
has_mrgsolve <- requireNamespace("mrgsolve", quietly = TRUE)
has_deSolve  <- requireNamespace("deSolve",  quietly = TRUE)

if (has_mrgsolve) {
  library(mrgsolve)
  message("[MN-App] Using mrgsolve for ODE simulation")
} else if (has_deSolve) {
  library(deSolve)
  message("[MN-App] mrgsolve not available; using deSolve")
} else {
  message("[MN-App] Neither mrgsolve nor deSolve available; using analytical approximations")
}

# =============================================================================
# INLINE mrgsolve MODEL CODE
# =============================================================================
MN_MODEL_CODE <- "
$PARAM
// --- Patient covariates ---
WT   = 70,      // body weight (kg)
AGE  = 45,      // age (years)
BEGFR = 60,     // baseline eGFR (mL/min/1.73m2)
BPROT = 6,      // baseline proteinuria (g/day)
BPLA2R = 200,   // baseline anti-PLA2R1 titer (RU/mL)
BALB  = 2.8,    // baseline serum albumin (g/dL)
STAGE = 2,      // disease stage (1/2/3)

// --- Rituximab PK (2-compartment) ---
RTX_CL  = 0.35,  // clearance (L/day) - pop estimate
RTX_V1  = 3.1,   // central volume (L)
RTX_Q   = 0.57,  // intercompartmental CL (L/day)
RTX_V2  = 3.7,   // peripheral volume (L)
RTX_KA  = 0,     // absorption (IV -> 0)
RTX_F   = 1,     // bioavailability

// --- Tacrolimus PK (1-compartment oral) ---
TAC_CL  = 2.25,  // apparent clearance (L/h) -> /24 for /day
TAC_V   = 101,   // apparent volume (L)
TAC_KA  = 1.4,   // absorption rate (1/h)
TAC_F   = 0.25,  // bioavailability

// --- Cyclophosphamide PK (1-compartment IV) ---
CYC_CL  = 5.1,   // clearance (L/h)
CYC_V   = 38,    // volume (L)

// --- PD parameters ---
EC50_RTX   = 0.5,   // RTX Cmin for 50% B-cell depletion (ug/mL)
EMAX_RTX   = 0.98,  // max B-cell depletion fraction
kout_Bcell = 0.02,  // B-cell recovery rate (1/day)
kin_Bcell  = 1,     // B-cell production rate (relative)

kd_PLA2R   = 0.015, // anti-PLA2R1 decline rate with B-cell depletion (1/day)
kr_PLA2R   = 0.003, // anti-PLA2R1 spontaneous decline (1/day)
ks_PLA2R   = 0.001, // spontaneous anti-PLA2R1 production (RU/mL/day)

kd_prot    = 0.02,  // proteinuria decline rate (1/day) at max effect
kr_prot    = 0.005, // proteinuria recovery rate (1/day)

kd_alb     = 0.08,  // albumin recovery rate (per unit proteinuria reduction)
kr_alb     = 0.02,  // albumin decline rate

kd_eGFR    = 0.0012,// eGFR decline per unit prot per day
kr_eGFR    = 0.008, // eGFR partial recovery per day

// --- Tacrolimus PD ---
EC50_TAC   = 5,     // trough ng/mL for 50% immune suppression
EMAX_TAC   = 0.7,

// --- Cyclophosphamide PD ---
EC50_CYC   = 1,     // ug/mL threshold
EMAX_CYC   = 0.85,

// --- Complement ---
kd_C3    = 0.01,
kd_sC5b9 = 0.015,

// --- Drug dosing flags ---
USE_RTX  = 0,
USE_TAC  = 0,
USE_CYC  = 0,
USE_ACE  = 0,
ACE_EFF  = 0.4,      // ACE inhibitor proteinuria reduction fraction

// --- Spontaneous remission ---
P_SPONT  = 0.03      // spontaneous remission probability/month (stage-dependent override in $MAIN)

$CMT
// Rituximab (2-compartment IV)
RTX_C  RTX_P
// Tacrolimus (oral 1-compartment)
TAC_GUT TAC_C
// Cyclophosphamide (1-compartment)
CYC_C
// PD compartments
BCELL     // CD20+ B cells (fraction of baseline)
PLA2R     // anti-PLA2R1 titer (RU/mL)
COMPL     // complement MAC activity (fraction)
ALB       // serum albumin (g/dL)
EGFR      // eGFR (mL/min/1.73m2)
PROT      // proteinuria (g/day)
IGG4      // serum IgG4 (relative)
LIPID     // LDL cholesterol (relative)

$MAIN
// Weight-based PK adjustments
double CL_RTX  = RTX_CL * pow(WT/70.0, 0.75);
double V1_RTX  = RTX_V1 * (WT/70.0);
double Q_RTX   = RTX_Q  * pow(WT/70.0, 0.75);
double V2_RTX  = RTX_V2 * (WT/70.0);

double CL_TAC  = (TAC_CL/24.0) * pow(WT/70.0, 0.75);
double V_TAC   = TAC_V * (WT/70.0);
double KA_TAC  = TAC_KA/24.0;

double CL_CYC  = (CYC_CL/24.0) * pow(WT/70.0, 0.75);
double V_CYC   = CYC_V * (WT/70.0);

// Initials
if(NEWIND <= 1) {
  RTX_C_0  = 0;
  RTX_P_0  = 0;
  TAC_GUT_0 = 0;
  TAC_C_0  = 0;
  CYC_C_0  = 0;
  BCELL_0  = 1.0;
  PLA2R_0  = BPLA2R;
  COMPL_0  = 1.0;
  ALB_0    = BALB;
  EGFR_0   = BEGFR;
  PROT_0   = BPROT;
  IGG4_0   = 1.0;
  LIPID_0  = 1.0 + (4.0 - BALB) * 0.3;
}

$ODE
// --- Rituximab PK ---
double k10_rtx = CL_RTX / V1_RTX;
double k12_rtx = Q_RTX  / V1_RTX;
double k21_rtx = Q_RTX  / V2_RTX;

dxdt_RTX_C = -k10_rtx * RTX_C - k12_rtx * RTX_C + k21_rtx * RTX_P;
dxdt_RTX_P =  k12_rtx * RTX_C - k21_rtx * RTX_P;

// RTX concentration in ug/mL (V in L, amount in mg -> ug/mL = mg/L * 1000/V)
double Crtx = RTX_C / V1_RTX * 1000.0;  // ug/mL
double Ertx = USE_RTX * EMAX_RTX * Crtx / (EC50_RTX + Crtx + 1e-6);

// --- Tacrolimus PK ---
dxdt_TAC_GUT = -KA_TAC * TAC_GUT;
dxdt_TAC_C   =  KA_TAC * TAC_GUT * TAC_F - (CL_TAC/V_TAC) * TAC_C;

double Ctac = TAC_C / V_TAC * 1000.0; // ng/mL (approx)
double Etac = USE_TAC * EMAX_TAC * Ctac / (EC50_TAC + Ctac + 1e-6);

// --- Cyclophosphamide PK ---
dxdt_CYC_C = -(CL_CYC/V_CYC) * CYC_C;
double Ccyc = CYC_C / V_CYC * 1000.0;
double Ecyc = USE_CYC * EMAX_CYC * Ccyc / (EC50_CYC + Ccyc + 1e-6);

// Combined immunosuppression effect (max principle)
double Eimm = fmax(fmax(Ertx, Etac), Ecyc);

// --- B cells ---
double Ebc = Eimm;
dxdt_BCELL = kin_Bcell * (1.0 - Ebc) - kout_Bcell * BCELL;
double Bcell_cur = fmax(BCELL, 0.001);

// --- anti-PLA2R1 ---
// B-cell depletion drives antibody clearance; spontaneous slow decline
double PLA2R_drain = (kd_PLA2R * (1.0 - Bcell_cur) + kr_PLA2R) * PLA2R;
double PLA2R_prod  = ks_PLA2R * Bcell_cur * PLA2R;
dxdt_PLA2R = -PLA2R_drain + PLA2R_prod;

// --- Complement MAC (sC5b-9) driven by PLA2R-immune complex ---
double norm_PLA2R = PLA2R / (BPLA2R + 1e-6);
dxdt_COMPL = kd_C3 * (norm_PLA2R - COMPL);

// --- IgG4 (parallels B cells and antibody production) ---
dxdt_IGG4 = -kd_sC5b9 * (1.0 - Bcell_cur) * IGG4 + 0.002 * Bcell_cur;

// --- Proteinuria (driven by complement MAC on podocytes) ---
double ACE_prot_red = USE_ACE * ACE_EFF;
double prot_drive  = COMPL * norm_PLA2R;  // immune injury driver
double prot_target = BPROT * prot_drive * (1.0 - ACE_prot_red) * (1.0 - 0.6 * Eimm);
dxdt_PROT  = 0.03 * (prot_target - PROT);

// --- Serum albumin ---
double prot_norm = PROT / (BPROT + 1e-6);
double alb_target = 4.2 - 1.5 * prot_norm;
dxdt_ALB   = 0.04 * (alb_target - ALB);

// --- eGFR ---
double egfr_target = BEGFR * (1.0 - 0.003 * STAGE) - kd_eGFR * PROT * t + kr_eGFR * (1.0 - prot_norm) * (BEGFR - EGFR);
dxdt_EGFR  = 0.01 * (egfr_target - EGFR);

// --- Lipids (rise with hypoalbuminemia) ---
double lipid_target = 1.0 + (4.2 - ALB) * 0.35;
dxdt_LIPID = 0.02 * (lipid_target - LIPID);

$TABLE
double RTX_CONC = RTX_C / V1_RTX * 1000.0;  // ug/mL
double TAC_CONC = TAC_C / V_TAC * 1000.0;   // ng/mL
double CYC_CONC = CYC_C / V_CYC;            // ug/mL (mg/L)
double BCELL_PCT = BCELL * 100.0;
double CR = (PROT < 0.3) ? 1.0 : 0.0;
double PR = (PROT < 3.5 && PROT < BPROT * 0.5) ? 1.0 : 0.0;

$CAPTURE RTX_CONC TAC_CONC CYC_CONC BCELL_PCT CR PR ALB EGFR PROT PLA2R COMPL IGG4 LIPID
"

# =============================================================================
# ODE SIMULATION FUNCTIONS
# =============================================================================

# Build mrgsolve model (if available)
build_mrgsolve_model <- function() {
  if (!has_mrgsolve) return(NULL)
  tryCatch(
    mrgsolve::mcode("mn_model", MN_MODEL_CODE, quiet = TRUE),
    error = function(e) { message("mrgsolve compile error: ", e$message); NULL }
  )
}

MN_MOD <- if (has_mrgsolve) build_mrgsolve_model() else NULL

# deSolve ODE system (fallback)
mn_ode_desolve <- function(t, state, params) {
  with(as.list(c(state, params)), {
    # Weight-based PK
    CL_RTX <- RTX_CL * (WT/70)^0.75
    V1_RTX <- RTX_V1 * (WT/70)
    Q_RTX  <- RTX_Q  * (WT/70)^0.75
    V2_RTX <- RTX_V2 * (WT/70)
    CL_TAC <- (TAC_CL/24) * (WT/70)^0.75
    V_TAC  <- TAC_V * (WT/70)
    KA_TAC <- TAC_KA/24
    CL_CYC <- (CYC_CL/24) * (WT/70)^0.75
    V_CYC  <- CYC_V * (WT/70)

    k10_rtx <- CL_RTX/V1_RTX
    k12_rtx <- Q_RTX/V1_RTX
    k21_rtx <- Q_RTX/V2_RTX

    dRTX_C  <- -k10_rtx*RTX_C - k12_rtx*RTX_C + k21_rtx*RTX_P
    dRTX_P  <-  k12_rtx*RTX_C - k21_rtx*RTX_P
    dTAC_GUT <- -KA_TAC*TAC_GUT
    dTAC_C   <-  KA_TAC*TAC_GUT*TAC_F - (CL_TAC/V_TAC)*TAC_C
    dCYC_C   <- -(CL_CYC/V_CYC)*CYC_C

    Crtx <- max(RTX_C/V1_RTX*1000, 0)
    Ctac <- max(TAC_C/V_TAC*1000, 0)
    Ccyc <- max(CYC_C/V_CYC*1000, 0)

    Ertx <- USE_RTX * EMAX_RTX * Crtx / (EC50_RTX + Crtx + 1e-6)
    Etac <- USE_TAC * EMAX_TAC * Ctac / (EC50_TAC + Ctac + 1e-6)
    Ecyc <- USE_CYC * EMAX_CYC * Ccyc / (EC50_CYC + Ccyc + 1e-6)
    Eimm <- max(Ertx, Etac, Ecyc)

    Bcell_cur <- max(BCELL, 0.001)
    dBCELL  <- kin_Bcell*(1-Eimm) - kout_Bcell*Bcell_cur
    norm_PLA2R <- PLA2R / (BPLA2R + 1e-6)
    dPLA2R  <- -(kd_PLA2R*(1-Bcell_cur) + kr_PLA2R)*PLA2R + ks_PLA2R*Bcell_cur*PLA2R
    dCOMPL  <- kd_C3*(norm_PLA2R - COMPL)
    dIGG4   <- -kd_sC5b9*(1-Bcell_cur)*IGG4 + 0.002*Bcell_cur
    ACE_red <- USE_ACE * ACE_EFF
    prot_drive  <- COMPL * norm_PLA2R
    prot_target <- BPROT * prot_drive * (1 - ACE_red) * (1 - 0.6*Eimm)
    dPROT   <- 0.03*(prot_target - PROT)
    prot_norm <- PROT / (BPROT + 1e-6)
    alb_target <- 4.2 - 1.5*prot_norm
    dALB    <- 0.04*(alb_target - ALB)
    egfr_target <- max(BEGFR*(1-0.003*STAGE) - 0.0012*PROT*t + 0.008*(1-prot_norm)*(BEGFR-EGFR), 10)
    dEGFR   <- 0.01*(egfr_target - EGFR)
    lipid_target <- 1 + (4.2 - ALB)*0.35
    dLIPID  <- 0.02*(lipid_target - LIPID)

    list(c(dRTX_C, dRTX_P, dTAC_GUT, dTAC_C, dCYC_C,
           dBCELL, dPLA2R, dCOMPL, dALB, dEGFR, dPROT, dIGG4, dLIPID))
  })
}

# Analytical / approximate fallback (no ODE solver)
simulate_analytical <- function(params, times_days) {
  with(params, {
    n <- length(times_days)
    t <- times_days

    # B-cell depletion kinetics
    RTX_t50   <- ifelse(USE_RTX == 1, 7, Inf)
    Bcell_nadir <- ifelse(USE_RTX == 1, 0.02, ifelse(USE_TAC == 1, 0.5, ifelse(USE_CYC == 1, 0.1, 1.0)))
    Bcell_rec <- 180
    BCELL_pct <- ifelse(t < RTX_t50,
                        1 - (1 - Bcell_nadir) * (t / RTX_t50),
                        Bcell_nadir + (1 - Bcell_nadir) * (1 - exp(-(t - RTX_t50) / Bcell_rec)))
    BCELL_pct <- pmax(pmin(BCELL_pct, 1), 0)

    # anti-PLA2R1
    PLA2R_t <- BPLA2R * exp(-(kd_PLA2R * (1 - BCELL_pct) + kr_PLA2R) * t) +
               ks_PLA2R * BCELL_pct * t
    PLA2R_t <- pmax(PLA2R_t, 0)

    # Complement
    norm_pla2r <- PLA2R_t / (BPLA2R + 1e-6)
    COMPL_t <- 1 - (1 - pmin(norm_pla2r, 1)) * (1 - exp(-kd_C3 * t))

    # Proteinuria
    ACE_red <- USE_ACE * ACE_EFF
    Eimm    <- ifelse(USE_RTX == 1, EMAX_RTX * (1 - BCELL_pct),
                      ifelse(USE_TAC == 1, 0.5,
                             ifelse(USE_CYC == 1, 0.7, 0)))
    prot_target <- BPROT * COMPL_t * norm_pla2r * (1 - ACE_red) * (1 - 0.6 * Eimm)
    prot_target <- pmax(prot_target, 0.05)
    PROT_t <- BPROT + (prot_target - BPROT) * (1 - exp(-0.03 * t))

    # Albumin
    prot_norm <- PROT_t / (BPROT + 1e-6)
    ALB_t <- 4.2 - 1.5 * prot_norm
    ALB_t <- pmax(pmin(ALB_t, 4.5), 1.0)

    # eGFR
    egfr_loss <- 0.0012 * PROT_t * t * (STAGE / 2)
    EGFR_t <- pmax(BEGFR - egfr_loss + 0.008 * (1 - prot_norm) * pmax(BEGFR - 20, 0), 5)

    # RTX PK (2-compartment, analytical bi-exponential)
    RTX_C <- rep(0, n)
    TAC_C <- rep(0, n)
    CYC_C <- rep(0, n)

    data.frame(time = t,
               RTX_CONC = RTX_C,
               TAC_CONC = TAC_C,
               CYC_CONC = CYC_C,
               BCELL_PCT = BCELL_pct * 100,
               PLA2R = PLA2R_t,
               COMPL = COMPL_t,
               ALB = ALB_t,
               EGFR = EGFR_t,
               PROT = PROT_t,
               IGG4 = pmax(1 - 0.6 * (1 - BCELL_pct), 0.1),
               LIPID = 1 + (4.2 - ALB_t) * 0.35)
  })
}

# Main simulation dispatcher
run_simulation <- function(params, times_days,
                           rtx_doses = NULL, rtx_times = NULL,
                           tac_doses = NULL, tac_interval = NULL,
                           cyc_doses = NULL, cyc_times = NULL) {
  if (!is.null(MN_MOD) && has_mrgsolve) {
    # --- mrgsolve path ---
    param_list <- as.list(params)

    ev <- mrgsolve::ev(time = 0, amt = 0, cmt = 1)  # dummy

    if (!is.null(rtx_doses) && !is.null(rtx_times)) {
      ev_rtx <- mrgsolve::ev(time = rtx_times, amt = rtx_doses, cmt = 1, rate = -2)
      ev <- mrgsolve::ev_seq(ev, ev_rtx)
    }
    if (!is.null(tac_doses) && !is.null(tac_interval)) {
      tac_t  <- seq(0, max(times_days), by = tac_interval)
      ev_tac <- mrgsolve::ev(time = tac_t, amt = tac_doses, cmt = 3)
      ev <- mrgsolve::ev_seq(ev, ev_tac)
    }
    if (!is.null(cyc_doses) && !is.null(cyc_times)) {
      ev_cyc <- mrgsolve::ev(time = cyc_times, amt = cyc_doses, cmt = 5, rate = -2)
      ev <- mrgsolve::ev_seq(ev, ev_cyc)
    }

    out <- tryCatch({
      MN_MOD %>%
        mrgsolve::param(param_list) %>%
        mrgsolve::mrgsim(events = ev, tgrid = times_days, obsonly = TRUE) %>%
        as.data.frame()
    }, error = function(e) NULL)

    if (!is.null(out)) return(out)
  }

  if (has_deSolve) {
    # --- deSolve fallback ---
    state <- c(RTX_C = 0, RTX_P = 0, TAC_GUT = 0, TAC_C = 0, CYC_C = 0,
               BCELL = 1, PLA2R = params$BPLA2R, COMPL = 1,
               ALB = params$BALB, EGFR = params$BEGFR, PROT = params$BPROT,
               IGG4 = 1, LIPID = 1 + (4.2 - params$BALB) * 0.3)

    p <- c(params,
           list(kin_Bcell = 1, kout_Bcell = 0.02,
                kd_PLA2R = 0.015, kr_PLA2R = 0.003, ks_PLA2R = 0.001,
                kd_C3 = 0.01, kd_sC5b9 = 0.015,
                RTX_CL = 0.35, RTX_V1 = 3.1, RTX_Q = 0.57, RTX_V2 = 3.7,
                TAC_CL = 2.25, TAC_V = 101, TAC_KA = 1.4, TAC_F = 0.25,
                CYC_CL = 5.1, CYC_V = 38,
                EC50_RTX = 0.5, EMAX_RTX = 0.98,
                EC50_TAC = 5, EMAX_TAC = 0.7,
                EC50_CYC = 1, EMAX_CYC = 0.85,
                ACE_EFF = 0.4))

    out <- tryCatch({
      res <- deSolve::ode(y = state, times = times_days, func = mn_ode_desolve, parms = p)
      df <- as.data.frame(res)
      names(df)[1] <- "time"
      df$RTX_CONC <- pmax(df$RTX_C / (p$RTX_V1 * (p$WT/70)) * 1000, 0)
      df$TAC_CONC <- pmax(df$TAC_C / (p$TAC_V * (p$WT/70)) * 1000, 0)
      df$CYC_CONC <- pmax(df$CYC_C / (p$CYC_V * (p$WT/70)) * 1000, 0)
      df$BCELL_PCT <- df$BCELL * 100
      df$PLA2R  <- df$PLA2R
      df$COMPL  <- df$COMPL
      df
    }, error = function(e) NULL)

    if (!is.null(out)) return(out)
  }

  # --- Pure analytical fallback ---
  simulate_analytical(params, times_days)
}

# PK only simulation for Tab 2
simulate_pk_only <- function(drug, dose, weight, tau_h = NULL, n_dose = 1) {
  if (drug == "Rituximab") {
    CL <- 0.35 * (weight/70)^0.75; V1 <- 3.1*(weight/70); Q <- 0.57*(weight/70)^0.75; V2 <- 3.7*(weight/70)
    t  <- seq(0, 21, by = 0.1)
    # Biexponential for IV bolus
    alpha <- ((Q/V1 + Q/V2 + CL/V1) + sqrt((Q/V1 + Q/V2 + CL/V1)^2 - 4*(CL/V1)*(Q/V2)))/2
    beta  <- ((Q/V1 + Q/V2 + CL/V1) - sqrt((Q/V1 + Q/V2 + CL/V1)^2 - 4*(CL/V1)*(Q/V2)))/2
    A <- (alpha - Q/V2) / ((alpha - beta)*V1)
    B <- (Q/V2 - beta)  / ((alpha - beta)*V1)
    Ct <- dose * (A*exp(-alpha*t) + B*exp(-beta*t)) * 1000  # ug/mL
    data.frame(time = t * 24, conc = pmax(Ct, 0), unit = "ug/mL", drug = "Rituximab")

  } else if (drug == "Tacrolimus") {
    CL_app <- (2.25/24) * (weight/70)^0.75; V_app <- 101*(weight/70); ka <- 1.4/24; F <- 0.25
    tau <- if (is.null(tau_h)) 12 else tau_h
    tau_d <- tau / 24
    t  <- seq(0, tau_d * max(n_dose, 7), by = 0.01)
    Ct <- rep(0, length(t))
    for (i in seq_len(n_dose)) {
      t0 <- (i-1)*tau_d
      tlag <- t - t0
      idx <- tlag > 0
      Ct[idx] <- Ct[idx] + (dose * ka * F / (V_app * (ka - CL_app/V_app))) *
        (exp(-(CL_app/V_app)*tlag[idx]) - exp(-ka*tlag[idx]))
    }
    data.frame(time = t * 24, conc = pmax(Ct, 0) * 1000, unit = "ng/mL", drug = "Tacrolimus")

  } else {  # Cyclophosphamide
    CL <- (5.1/24) * (weight/70)^0.75; V <- 38*(weight/70)
    t  <- seq(0, 2, by = 0.01)
    Ct <- (dose / V) * exp(-(CL/V)*t) * 1000  # ug/mL
    data.frame(time = t * 24, conc = pmax(Ct, 0), unit = "ug/mL", drug = "Cyclophosphamide")
  }
}

# Compute PK metrics
pk_metrics <- function(pk_df) {
  cmax <- max(pk_df$conc, na.rm = TRUE)
  t_cmax <- pk_df$time[which.max(pk_df$conc)]
  # AUC by trapezoidal
  auc <- sum(diff(pk_df$time) * (head(pk_df$conc, -1) + tail(pk_df$conc, -1)) / 2)
  trough <- min(tail(pk_df$conc, 10), na.rm = TRUE)
  # Half-life from terminal slope
  idx_term <- which(pk_df$time > t_cmax)
  if (length(idx_term) > 5) {
    fit <- tryCatch(lm(log(pmax(pk_df$conc[idx_term], 1e-10)) ~ pk_df$time[idx_term]), error = function(e) NULL)
    t12 <- if (!is.null(fit)) -log(2)/coef(fit)[2] else NA
  } else t12 <- NA
  data.frame(Metric = c("Cmax", "Tmax (h)", "AUC", "Trough", "Half-life (h)"),
             Value  = round(c(cmax, t_cmax, auc, trough, t12), 2),
             Unit   = c(pk_df$unit[1], "h", paste0(pk_df$unit[1], "*h"), pk_df$unit[1], "h"))
}

# =============================================================================
# SCENARIO DEFINITIONS
# =============================================================================
SCENARIOS <- list(
  "No Treatment"                    = list(USE_RTX=0, USE_TAC=0, USE_CYC=0, USE_ACE=0, color="#E41A1C"),
  "Rituximab Monotherapy"           = list(USE_RTX=1, USE_TAC=0, USE_CYC=0, USE_ACE=0, color="#377EB8"),
  "Tacrolimus Monotherapy"          = list(USE_RTX=0, USE_TAC=1, USE_CYC=0, USE_ACE=0, color="#4DAF4A"),
  "Cyclophosphamide + Steroids"     = list(USE_RTX=0, USE_TAC=0, USE_CYC=1, USE_ACE=0, color="#FF7F00"),
  "Rituximab + Tacrolimus"          = list(USE_RTX=1, USE_TAC=1, USE_CYC=0, USE_ACE=0, color="#984EA3"),
  "ACE Inhibitor (Conservative)"    = list(USE_RTX=0, USE_TAC=0, USE_CYC=0, USE_ACE=1, color="#A65628")
)

# =============================================================================
# COLOR SCHEME
# =============================================================================
COLS <- list(
  healthy   = "#2166AC",  # blue
  worsening = "#D6604D",  # red
  remission = "#1A9850",  # green
  warning   = "#F4A736",  # orange
  neutral   = "#636363",  # grey
  bg_light  = "#F7FBFF"
)

THEME_MN <- theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#2166AC"),
        strip.text = element_text(color = "white", face = "bold"))

# =============================================================================
# UI
# =============================================================================
ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Membranous Nephropathy QSP"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Profile",     tabName = "tab_patient",   icon = icon("user")),
      menuItem("PK / Pharmacokinetics", tabName = "tab_pk",      icon = icon("chart-line")),
      menuItem("PD — Immune Markers", tabName = "tab_pd_immune", icon = icon("dna")),
      menuItem("Clinical Endpoints",  tabName = "tab_clinical",  icon = icon("heartbeat")),
      menuItem("Scenario Comparison", tabName = "tab_scenario",  icon = icon("balance-scale")),
      menuItem("Biomarker Dashboard", tabName = "tab_biomarker", icon = icon("tachometer-alt"))
    ),
    hr(),
    div(style = "padding:10px; font-size:11px; color:#aaa;",
        "MN QSP Model v1.0",
        br(), "ODE via mrgsolve / deSolve",
        br(), "For research use only")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .skin-blue .main-sidebar { background-color: #1a3a5c; }
      .skin-blue .sidebar-menu > li.active > a { background-color: #2166AC; border-left-color: #5AADE0; }
      .info-box-icon { background-color: #2166AC !important; }
      .box.box-primary { border-top-color: #2166AC; }
      .traffic-green  { background-color: #1A9850; color: white; border-radius: 50%; width:30px; height:30px;
                        display:inline-block; text-align:center; line-height:30px; font-weight:bold; }
      .traffic-yellow { background-color: #F4A736; color: white; border-radius: 50%; width:30px; height:30px;
                        display:inline-block; text-align:center; line-height:30px; font-weight:bold; }
      .traffic-red    { background-color: #D6604D; color: white; border-radius: 50%; width:30px; height:30px;
                        display:inline-block; text-align:center; line-height:30px; font-weight:bold; }
      .remission-badge { padding:5px 10px; border-radius:5px; font-weight:bold; font-size:14px; display:inline-block; }
    "))),

    tabItems(
      # -----------------------------------------------------------------------
      # TAB 1: PATIENT PROFILE
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_patient",
        fluidRow(
          box(title = "Patient Demographics", width = 4, status = "primary", solidHeader = TRUE,
              sliderInput("age",    "Age (years)",              min=18, max=80, value=45, step=1),
              sliderInput("weight","Weight (kg)",               min=40, max=120,value=70, step=1),
              sliderInput("begfr", "Baseline eGFR (mL/min/1.73m²)", min=20,max=120,value=60,step=1),
              sliderInput("bprot", "Baseline Proteinuria (g/day)",   min=1, max=15, value=6, step=0.5),
              sliderInput("bpla2r","Anti-PLA2R1 Titer (RU/mL)",      min=0, max=500,value=200,step=10),
              sliderInput("balb",  "Serum Albumin (g/dL)",            min=1.5,max=4.5,value=2.8,step=0.1),
              selectInput("stage","Disease Stage",
                          choices = c("Stage I (low risk)"=1,
                                      "Stage II (moderate risk)"=2,
                                      "Stage III (high risk)"=3), selected=2)
          ),
          box(title = "Risk Assessment & Clinical Summary", width = 8, status = "primary", solidHeader = TRUE,
              fluidRow(
                valueBoxOutput("vbox_prot",   width=4),
                valueBoxOutput("vbox_alb",    width=4),
                valueBoxOutput("vbox_egfr",   width=4)
              ),
              fluidRow(
                valueBoxOutput("vbox_pla2r",  width=4),
                valueBoxOutput("vbox_spont",  width=4),
                valueBoxOutput("vbox_esrd",   width=4)
              ),
              hr(),
              h4("Toronto Risk Score Components", style="font-weight:bold;"),
              tableOutput("risk_table"),
              hr(),
              h4("Disease Interpretation", style="font-weight:bold;"),
              uiOutput("disease_interp")
          )
        ),
        fluidRow(
          box(title = "Baseline Clinical Profile Radar", width = 12, status = "info",
              plotOutput("patient_radar", height = "300px"))
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 2: PK / PHARMACOKINETICS
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_pk",
        fluidRow(
          box(title = "Drug & Dosing Parameters", width = 3, status = "primary", solidHeader = TRUE,
              selectInput("pk_drug", "Select Drug",
                          choices = c("Rituximab", "Tacrolimus", "Cyclophosphamide")),
              conditionalPanel("input.pk_drug == 'Rituximab'",
                numericInput("rtx_dose", "Dose (mg/infusion)", value = 1000, min = 375, max = 2000, step = 125),
                numericInput("rtx_ndose","Number of Infusions", value = 2, min=1, max=4),
                numericInput("rtx_interval","Interval (days)", value = 14, min=7, max=28),
                helpText("Standard: 1000 mg x2 (days 1, 15) or 375 mg/m² x4 weekly")
              ),
              conditionalPanel("input.pk_drug == 'Tacrolimus'",
                numericInput("tac_dose","Dose per administration (mg)", value=2, min=0.5, max=10, step=0.5),
                numericInput("tac_tau", "Dosing interval (hours)", value=12, min=8, max=24),
                numericInput("tac_ndose","Number of doses to simulate", value=14, min=1, max=60),
                helpText("Target trough: 3–8 ng/mL for MN")
              ),
              conditionalPanel("input.pk_drug == 'Cyclophosphamide'",
                numericInput("cyc_dose","Dose (mg/pulse)", value=1000, min=500, max=2500, step=100),
                helpText("IV pulse: 0.5–1 g/m² monthly x6 (Ponticelli protocol)")
              ),
              actionButton("run_pk", "Simulate PK", class="btn-primary", icon=icon("play"))
          ),
          box(title = "PK Concentration-Time Profile", width = 9, status = "primary", solidHeader = TRUE,
              plotOutput("pk_plot", height="350px"),
              hr(),
              h4("Key PK Metrics", style="font-weight:bold;"),
              tableOutput("pk_metrics_table")
          )
        ),
        fluidRow(
          box(title = "PK / Clinical Notes", width = 12, status = "info",
              uiOutput("pk_notes"))
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 3: PD - IMMUNE MARKERS
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_pd_immune",
        fluidRow(
          box(title = "Treatment Selection", width = 3, status = "primary", solidHeader = TRUE,
              checkboxInput("pd_rtx", "Rituximab", value = TRUE),
              checkboxInput("pd_tac", "Tacrolimus", value = FALSE),
              checkboxInput("pd_cyc", "Cyclophosphamide", value = FALSE),
              checkboxInput("pd_ace", "ACE Inhibitor", value = FALSE),
              sliderInput("pd_months","Simulation Duration (months)", min=6, max=36, value=24, step=3),
              actionButton("run_pd", "Simulate", class="btn-primary", icon=icon("play")),
              hr(),
              helpText("Rituximab depletes CD20+ B cells, which produce anti-PLA2R1 IgG4 that activates complement, injuring podocytes and causing proteinuria.")
          ),
          box(title = "B Cell Kinetics & Antibody Response", width = 9, status = "primary", solidHeader = TRUE,
              plotOutput("pd_bcell_plot", height="280px"),
              hr(),
              plotOutput("pd_pla2r_plot", height="280px")
          )
        ),
        fluidRow(
          box(title = "Complement & Immunological Remission", width = 6, status = "info",
              plotOutput("pd_compl_plot", height="280px")),
          box(title = "Serum IgG4 Levels", width = 6, status = "info",
              plotOutput("pd_igg4_plot", height="280px"))
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 4: CLINICAL ENDPOINTS
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_clinical",
        fluidRow(
          box(title = "Treatment Settings", width = 3, status = "primary", solidHeader = TRUE,
              checkboxInput("clin_rtx", "Rituximab", value = TRUE),
              checkboxInput("clin_tac", "Tacrolimus", value = FALSE),
              checkboxInput("clin_cyc", "Cyclophosphamide", value = FALSE),
              checkboxInput("clin_ace", "ACE Inhibitor", value = TRUE),
              sliderInput("clin_months","Observation Period (months)", min=6, max=36, value=24),
              actionButton("run_clin", "Simulate", class="btn-primary", icon=icon("play")),
              hr(),
              h5("Remission Criteria:"),
              tags$ul(
                tags$li(tags$b("Complete:"), " Proteinuria < 0.3 g/day"),
                tags$li(tags$b("Partial:"), " > 50% reduction AND < 3.5 g/day")
              )
          ),
          box(title = "Remission Status", width = 9, status = "primary", solidHeader = TRUE,
              fluidRow(
                column(4, uiOutput("remission_badge")),
                column(4, uiOutput("time_to_pr_box")),
                column(4, uiOutput("time_to_cr_box"))
              ),
              hr(),
              plotOutput("clin_prot_plot", height="250px")
          )
        ),
        fluidRow(
          box(title = "Serum Albumin", width = 4, status = "info",
              plotOutput("clin_alb_plot", height="220px")),
          box(title = "eGFR Trajectory", width = 4, status = "info",
              plotOutput("clin_egfr_plot", height="220px")),
          box(title = "Serum Creatinine (estimated)", width = 4, status = "info",
              plotOutput("clin_scr_plot", height="220px"))
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 5: SCENARIO COMPARISON
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_scenario",
        fluidRow(
          box(title = "Select Treatment Scenarios", width = 3, status = "primary", solidHeader = TRUE,
              lapply(names(SCENARIOS), function(s) {
                checkboxInput(paste0("sc_", gsub("[^A-Za-z0-9]", "_", s)), s,
                              value = s %in% c("No Treatment","Rituximab Monotherapy","ACE Inhibitor (Conservative)"))
              }),
              sliderInput("sc_months","Comparison Period (months)", min=12, max=36, value=24),
              actionButton("run_scenario","Compare Scenarios", class="btn-primary btn-lg btn-block", icon=icon("chart-bar"))
          ),
          box(title = "Proteinuria Over Time by Scenario", width = 9, status = "primary", solidHeader = TRUE,
              plotOutput("sc_prot_plot", height="280px"),
              plotOutput("sc_egfr_plot", height="280px"))
        ),
        fluidRow(
          box(title = "Anti-PLA2R1 Titer Comparison", width = 6, status = "info",
              plotOutput("sc_pla2r_plot", height="260px")),
          box(title = "Scenario Outcome Summary Table", width = 6, status = "info",
              DTOutput("sc_summary_table"))
        )
      ),

      # -----------------------------------------------------------------------
      # TAB 6: BIOMARKER DASHBOARD
      # -----------------------------------------------------------------------
      tabItem(tabName = "tab_biomarker",
        fluidRow(
          box(title = "Biomarker Status Overview", width = 12, status = "primary", solidHeader = TRUE,
              fluidRow(
                column(2, div(class="text-center", h5("Anti-PLA2R1"), uiOutput("bm_pla2r_light"))),
                column(2, div(class="text-center", h5("Proteinuria"), uiOutput("bm_prot_light"))),
                column(2, div(class="text-center", h5("Albumin"),     uiOutput("bm_alb_light"))),
                column(2, div(class="text-center", h5("eGFR"),        uiOutput("bm_egfr_light"))),
                column(2, div(class="text-center", h5("sC5b-9"),      uiOutput("bm_compl_light"))),
                column(2, div(class="text-center", h5("Thrombosis"),  uiOutput("bm_thrombo_light")))
              )
          )
        ),
        fluidRow(
          box(title = "Anti-PLA2R1 & Complement Trajectory (with selected Tx)", width = 6, status = "primary",
              checkboxInput("bm_rtx","Rituximab",  TRUE),
              checkboxInput("bm_tac","Tacrolimus",FALSE),
              checkboxInput("bm_cyc","Cyclophosphamide",FALSE),
              actionButton("run_bm","Update",class="btn-info"),
              plotOutput("bm_pla2r_plot",  height="220px"),
              plotOutput("bm_compl_plot",  height="220px")
          ),
          box(title = "Renal & Metabolic Biomarkers", width = 6, status = "info",
              plotOutput("bm_renal_plot",  height="220px"),
              plotOutput("bm_lipid_plot",  height="220px")
          )
        ),
        fluidRow(
          box(title = "Thrombosis Risk Score", width = 12, status = "warning",
              uiOutput("thrombo_score_ui"),
              plotOutput("thrombo_gauge", height="200px"))
        )
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  # ---- Reactive: base patient params ----
  base_params <- reactive({
    list(
      WT     = input$weight,
      AGE    = input$age,
      BEGFR  = input$begfr,
      BPROT  = input$bprot,
      BPLA2R = input$bpla2r,
      BALB   = input$balb,
      STAGE  = as.numeric(input$stage)
    )
  })

  # ==========================================================================
  # TAB 1: PATIENT PROFILE
  # ==========================================================================
  output$vbox_prot <- renderValueBox({
    col <- if (input$bprot < 3.5) "yellow" else "red"
    valueBox(paste0(input$bprot, " g/day"), "Baseline Proteinuria", icon=icon("tint"), color=col)
  })
  output$vbox_alb <- renderValueBox({
    col <- if (input$balb >= 3.0) "green" else if (input$balb >= 2.0) "yellow" else "red"
    valueBox(paste0(input$balb, " g/dL"), "Serum Albumin", icon=icon("flask"), color=col)
  })
  output$vbox_egfr <- renderValueBox({
    col <- if (input$begfr >= 60) "green" else if (input$begfr >= 30) "yellow" else "red"
    valueBox(paste0(input$begfr), "eGFR (mL/min/1.73m²)", icon=icon("filter"), color=col)
  })
  output$vbox_pla2r <- renderValueBox({
    col <- if (input$bpla2r < 50) "green" else if (input$bpla2r < 150) "yellow" else "red"
    valueBox(paste0(input$bpla2r, " RU/mL"), "Anti-PLA2R1", icon=icon("microscope"), color=col)
  })
  output$vbox_spont <- renderValueBox({
    # Toronto criteria: ~30% spontaneous remission at 24 months, higher for low titer
    stage_val <- as.numeric(input$stage)
    p_spont <- round(max(0.05,
                         0.35 - 0.08 * stage_val -
                           0.0003 * input$bpla2r +
                           0.01 * pmax(input$begfr - 60, 0) / 10), 2)
    p_spont <- min(p_spont, 0.50)
    valueBox(paste0(round(p_spont * 100), "%"), "Spontaneous Remission (24 mo)",
             icon=icon("chart-pie"), color="light-blue")
  })
  output$vbox_esrd <- renderValueBox({
    stage_val <- as.numeric(input$stage)
    p_esrd <- round(pmin(0.05 + 0.12*(stage_val-1) +
                           0.002 * pmax(60 - input$begfr, 0) +
                           0.0008 * input$bpla2r, 0.6), 2)
    col <- if (p_esrd < 0.1) "green" else if (p_esrd < 0.25) "yellow" else "red"
    valueBox(paste0(round(p_esrd*100), "%"), "ESRD Risk (5-year)", icon=icon("exclamation-triangle"), color=col)
  })

  output$risk_table <- renderTable({
    stage_val <- as.numeric(input$stage)
    data.frame(
      Factor = c("Proteinuria", "eGFR", "Anti-PLA2R1", "Albumin", "Stage"),
      Value  = c(paste0(input$bprot, " g/day"), paste0(input$begfr, " mL/min/1.73m²"),
                 paste0(input$bpla2r, " RU/mL"), paste0(input$balb, " g/dL"),
                 paste0("Stage ", stage_val)),
      `Risk Level` = c(
        ifelse(input$bprot < 3.5, "Low", ifelse(input$bprot < 8, "Moderate", "High")),
        ifelse(input$begfr >= 60, "Low", ifelse(input$begfr >= 30, "Moderate", "High")),
        ifelse(input$bpla2r < 50, "Low", ifelse(input$bpla2r < 200, "Moderate", "High")),
        ifelse(input$balb >= 3.0, "Normal", ifelse(input$balb >= 2.0, "Low", "Very Low")),
        c("I — Low","II — Moderate","III — High")[stage_val]
      )
    )
  }, striped = TRUE, hover = TRUE)

  output$disease_interp <- renderUI({
    prot <- input$bprot; alb <- input$balb; egfr <- input$begfr; pla2r <- input$bpla2r
    msgs <- list()
    if (prot > 8)   msgs[[length(msgs)+1]] <- tags$li(tags$b("Severe nephrotic syndrome"), ": High proteinuria warrants immediate immunosuppression consideration.")
    if (alb < 2.5)  msgs[[length(msgs)+1]] <- tags$li(tags$b("Hypoalbuminemia"), ": Increased risk of thrombotic complications. Anticoagulation may be indicated.")
    if (egfr < 45)  msgs[[length(msgs)+1]] <- tags$li(tags$b("Reduced renal function"), ": CKD stage ≥3. Closer monitoring required; Tacrolimus dose adjustment needed.")
    if (pla2r > 300)msgs[[length(msgs)+1]] <- tags$li(tags$b("High anti-PLA2R1"), ": Unlikely to spontaneously remit. Early rituximab therapy recommended (MENTOR trial evidence).")
    if (length(msgs) == 0) msgs[[1]] <- tags$li("Moderate disease activity. Consider watchful waiting with supportive therapy initially.")
    tagList(tags$ul(msgs))
  })

  output$patient_radar <- renderPlot({
    df <- data.frame(
      Variable = c("Proteinuria\n(inverse)", "eGFR", "Albumin", "Anti-PLA2R1\n(inverse)", "Renal\nReserve"),
      Score    = c(
        1 - pmin(input$bprot / 15, 1),
        pmin(input$begfr / 120, 1),
        pmin((input$balb - 1.5) / 3, 1),
        1 - pmin(input$bpla2r / 500, 1),
        pmin((input$begfr - 20) / 100, 1)
      )
    )
    ggplot(df, aes(x = Variable, y = Score, fill = Score)) +
      geom_col(alpha = 0.8, color = "white", width = 0.6) +
      geom_hline(yintercept = 0.6, linetype = "dashed", color = COLS$remission, size = 0.8) +
      scale_fill_gradient(low = COLS$worsening, high = COLS$healthy, limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      labs(title = "Patient Profile: Normalized Clinical Indicators",
           subtitle = "Green dashed line = target/healthy threshold",
           x = NULL, y = "Relative Score (0 = worst, 1 = best)") +
      THEME_MN + theme(legend.position = "none")
  })

  # ==========================================================================
  # TAB 2: PHARMACOKINETICS
  # ==========================================================================
  pk_sim <- eventReactive(input$run_pk, {
    drug <- input$pk_drug
    wt   <- input$weight
    if (drug == "Rituximab") {
      dose <- input$rtx_dose * input$rtx_ndose  # total
      simulate_pk_only("Rituximab", input$rtx_dose, wt)
    } else if (drug == "Tacrolimus") {
      simulate_pk_only("Tacrolimus", input$tac_dose, wt,
                       tau_h = input$tac_tau, n_dose = input$tac_ndose)
    } else {
      simulate_pk_only("Cyclophosphamide", input$cyc_dose, wt)
    }
  }, ignoreNULL = FALSE)

  output$pk_plot <- renderPlot({
    df <- pk_sim()
    ylab <- paste0("Concentration (", df$unit[1], ")")
    gg <- ggplot(df, aes(x = time, y = conc)) +
      geom_line(color = COLS$healthy, size = 1.2) +
      geom_area(fill = COLS$healthy, alpha = 0.15) +
      labs(title = paste(df$drug[1], "PK Profile"),
           x = "Time (hours)", y = ylab) +
      THEME_MN

    if (df$drug[1] == "Rituximab") {
      gg <- gg + geom_hline(yintercept = 0.5, linetype = "dashed",
                            color = COLS$remission, size = 0.8) +
        annotate("text", x = max(df$time)*0.7, y = 0.7, label = "EC50 threshold (0.5 ug/mL)",
                 color = COLS$remission, size = 3.5)
    }
    if (df$drug[1] == "Tacrolimus") {
      gg <- gg +
        geom_hline(yintercept = 3, linetype = "dashed", color = COLS$remission, size = 0.7) +
        geom_hline(yintercept = 8, linetype = "dashed", color = COLS$worsening, size = 0.7) +
        annotate("text", x = 5, y = 4, label = "Therapeutic window (3–8 ng/mL)",
                 color = COLS$neutral, size = 3.5)
    }
    gg
  })

  output$pk_metrics_table <- renderTable({
    df <- pk_sim()
    pk_metrics(df)
  }, striped = TRUE, hover = TRUE, width = "100%")

  output$pk_notes <- renderUI({
    drug <- input$pk_drug
    notes <- switch(drug,
      "Rituximab" = tagList(
        tags$h5("Rituximab PK Notes"),
        tags$ul(
          tags$li("Two-compartment IV model: CL ~0.35 L/day (WT-adjusted), t½β ~20 days"),
          tags$li("MENTOR trial (2019): 1000 mg x2 (D1, D15) — superior to Cyclosporine at 24 months"),
          tags$li("B-cell depletion maintained for 6–12 months; re-dosing when B cells > 1% or anti-PLA2R1 rises"),
          tags$li("Reference: Fervenza FC et al., NEJM 2019; 381:36-46")
        )
      ),
      "Tacrolimus" = tagList(
        tags$h5("Tacrolimus PK Notes"),
        tags$ul(
          tags$li("1-compartment oral model: high variability (CV ~40%), F ~25%"),
          tags$li("Target trough 3–8 ng/mL for MN; CYP3A4/P-gp interactions important"),
          tags$li("Inhibits calcineurin → reduces IL-2 → T-cell suppression"),
          tags$li("Renal toxicity risk (vasoconstriction) at high troughs; eGFR monitoring essential"),
          tags$li("Reference: Chen M et al., JASN 2019; Praga M et al., JASN 2007")
        )
      ),
      "Cyclophosphamide" = tagList(
        tags$h5("Cyclophosphamide PK Notes"),
        tags$ul(
          tags$li("IV pulse (Ponticelli protocol): 0.5–1 g/m² monthly x 6 cycles"),
          tags$li("Pro-drug activated by CYP2B6; active metabolites cause DNA crosslinks"),
          tags$li("Combined with alternate-month steroids in Ponticelli regimen"),
          tags$li("Reference: Ponticelli C et al., JASN 1998; Jha V et al., NEJM 2007")
        )
      )
    )
    notes
  })

  # ==========================================================================
  # TAB 3: PD - IMMUNE MARKERS
  # ==========================================================================
  pd_sim <- eventReactive(input$run_pd, {
    t_days <- seq(0, input$pd_months * 30.4, by = 1)
    p <- c(base_params(),
           list(USE_RTX = as.numeric(input$pd_rtx),
                USE_TAC = as.numeric(input$pd_tac),
                USE_CYC = as.numeric(input$pd_cyc),
                USE_ACE = as.numeric(input$pd_ace)))
    run_simulation(p, t_days)
  }, ignoreNULL = FALSE)

  output$pd_bcell_plot <- renderPlot({
    df <- pd_sim()
    bcell_col <- if ("BCELL_PCT" %in% names(df)) "BCELL_PCT" else "BCELL"
    ggplot(df, aes(x = time/30.4, y = .data[[bcell_col]])) +
      geom_ribbon(aes(ymin = 0, ymax = .data[[bcell_col]]), fill = COLS$healthy, alpha = 0.2) +
      geom_line(color = COLS$healthy, size = 1.2) +
      geom_hline(yintercept = 1, linetype = "dashed", color = COLS$neutral, size = 0.7) +
      labs(title = "CD20+ B Cell Depletion and Reconstitution",
           x = "Time (months)", y = "B Cells (% of baseline)") +
      THEME_MN
  })

  output$pd_pla2r_plot <- renderPlot({
    df <- pd_sim()
    ggplot(df, aes(x = time/30.4, y = PLA2R)) +
      geom_ribbon(aes(ymin = 0, ymax = PLA2R), fill = COLS$worsening, alpha = 0.15) +
      geom_line(color = COLS$worsening, size = 1.2) +
      geom_hline(yintercept = 14, linetype = "dashed", color = COLS$remission, size = 0.8) +
      annotate("text", x = 1, y = 20, label = "Serological remission threshold (14 RU/mL)",
               hjust = 0, color = COLS$remission, size = 3.5) +
      labs(title = "Anti-PLA2R1 IgG4 Titer",
           x = "Time (months)", y = "Anti-PLA2R1 (RU/mL)") +
      THEME_MN
  })

  output$pd_compl_plot <- renderPlot({
    df <- pd_sim()
    ggplot(df, aes(x = time/30.4, y = COMPL)) +
      geom_line(color = COLS$warning, size = 1.2) +
      geom_hline(yintercept = 0.2, linetype = "dashed", color = COLS$remission, size = 0.7) +
      labs(title = "Complement MAC Activity (sC5b-9)",
           x = "Time (months)", y = "MAC Activity (relative, 1 = baseline)") +
      THEME_MN
  })

  output$pd_igg4_plot <- renderPlot({
    df <- pd_sim()
    ggplot(df, aes(x = time/30.4, y = IGG4)) +
      geom_line(color = "#984EA3", size = 1.2) +
      labs(title = "Serum IgG4 Levels",
           x = "Time (months)", y = "IgG4 (relative to baseline)") +
      THEME_MN
  })

  # ==========================================================================
  # TAB 4: CLINICAL ENDPOINTS
  # ==========================================================================
  clin_sim <- eventReactive(input$run_clin, {
    t_days <- seq(0, input$clin_months * 30.4, by = 1)
    p <- c(base_params(),
           list(USE_RTX = as.numeric(input$clin_rtx),
                USE_TAC = as.numeric(input$clin_tac),
                USE_CYC = as.numeric(input$clin_cyc),
                USE_ACE = as.numeric(input$clin_ace)))
    run_simulation(p, t_days)
  }, ignoreNULL = FALSE)

  output$remission_badge <- renderUI({
    df <- clin_sim()
    latest_prot <- tail(df$PROT, 1)
    baseline_prot <- input$bprot
    status <- if (latest_prot < 0.3) {
      div(class="remission-badge", style="background:#1A9850;color:white;", "COMPLETE REMISSION")
    } else if (latest_prot < 3.5 && latest_prot < baseline_prot * 0.5) {
      div(class="remission-badge", style="background:#F4A736;color:white;", "PARTIAL REMISSION")
    } else {
      div(class="remission-badge", style="background:#D6604D;color:white;", "NO REMISSION")
    }
    tagList(h5("Final Remission Status"), status)
  })

  output$time_to_pr_box <- renderUI({
    df <- clin_sim()
    baseline <- input$bprot
    pr_idx <- which(df$PROT < 3.5 & df$PROT < baseline * 0.5)
    t_pr <- if (length(pr_idx) > 0) round(df$time[min(pr_idx)] / 30.4, 1) else NA
    tagList(h5("Time to Partial Remission"),
            div(class="remission-badge",
                style = if (!is.na(t_pr)) "background:#4292C6;color:white;" else "background:#636363;color:white;",
                if (!is.na(t_pr)) paste0(t_pr, " months") else "Not achieved"))
  })

  output$time_to_cr_box <- renderUI({
    df <- clin_sim()
    cr_idx <- which(df$PROT < 0.3)
    t_cr <- if (length(cr_idx) > 0) round(df$time[min(cr_idx)] / 30.4, 1) else NA
    tagList(h5("Time to Complete Remission"),
            div(class="remission-badge",
                style = if (!is.na(t_cr)) "background:#1A9850;color:white;" else "background:#636363;color:white;",
                if (!is.na(t_cr)) paste0(t_cr, " months") else "Not achieved"))
  })

  output$clin_prot_plot <- renderPlot({
    df <- clin_sim()
    ggplot(df, aes(x = time/30.4, y = PROT)) +
      geom_ribbon(aes(ymin = 0, ymax = PROT), fill = COLS$worsening, alpha = 0.15) +
      geom_line(color = COLS$worsening, size = 1.2) +
      geom_hline(yintercept = 0.3, linetype = "dashed", color = COLS$remission,  size = 0.9) +
      geom_hline(yintercept = 3.5, linetype = "dashed", color = COLS$warning,    size = 0.9) +
      annotate("text", x = 0.5, y = 0.15, label = "CR threshold (<0.3)", hjust=0, color=COLS$remission, size=3.5) +
      annotate("text", x = 0.5, y = 3.8,  label = "PR threshold (<3.5)", hjust=0, color=COLS$warning,   size=3.5) +
      labs(title = "Proteinuria Time Course", x = "Time (months)", y = "Proteinuria (g/day)") +
      THEME_MN
  })

  output$clin_alb_plot <- renderPlot({
    df <- clin_sim()
    ggplot(df, aes(x = time/30.4, y = ALB)) +
      geom_line(color = COLS$healthy, size = 1.2) +
      geom_hline(yintercept = 3.5, linetype = "dashed", color = COLS$remission, size = 0.8) +
      labs(title = "Serum Albumin", x = "Time (months)", y = "Albumin (g/dL)") +
      THEME_MN
  })

  output$clin_egfr_plot <- renderPlot({
    df <- clin_sim()
    ggplot(df, aes(x = time/30.4, y = EGFR)) +
      geom_line(color = COLS$healthy, size = 1.2) +
      geom_hline(yintercept = 60, linetype = "dashed", color = COLS$warning, size = 0.8) +
      labs(title = "eGFR", x = "Time (months)", y = "eGFR (mL/min/1.73m²)") +
      THEME_MN
  })

  output$clin_scr_plot <- renderPlot({
    df <- clin_sim()
    # Cockcroft-Gault inverse: SCr ~ (140-age)*wt/(72*eGFR) (male) *0.85 female
    df$SCr <- (140 - input$age) * input$weight / (72 * pmax(df$EGFR, 5))
    ggplot(df, aes(x = time/30.4, y = SCr)) +
      geom_line(color = COLS$worsening, size = 1.2) +
      geom_hline(yintercept = 1.2, linetype = "dashed", color = COLS$warning, size = 0.8) +
      labs(title = "Serum Creatinine (estimated)", x = "Time (months)", y = "Creatinine (mg/dL)") +
      THEME_MN
  })

  # ==========================================================================
  # TAB 5: SCENARIO COMPARISON
  # ==========================================================================
  scenario_sims <- eventReactive(input$run_scenario, {
    t_days <- seq(0, input$sc_months * 30.4, by = 2)
    bp <- base_params()

    results <- list()
    for (sc_name in names(SCENARIOS)) {
      sc_id <- paste0("sc_", gsub("[^A-Za-z0-9]", "_", sc_name))
      if (isTRUE(input[[sc_id]])) {
        sc_params <- SCENARIOS[[sc_name]]
        p <- c(bp, list(
          USE_RTX = sc_params$USE_RTX,
          USE_TAC = sc_params$USE_TAC,
          USE_CYC = sc_params$USE_CYC,
          USE_ACE = sc_params$USE_ACE
        ))
        sim <- run_simulation(p, t_days)
        sim$Scenario <- sc_name
        sim$color    <- sc_params$color
        results[[sc_name]] <- sim
      }
    }
    if (length(results) == 0) return(NULL)
    bind_rows(results)
  }, ignoreNULL = FALSE)

  output$sc_prot_plot <- renderPlot({
    df <- scenario_sims()
    if (is.null(df)) return(NULL)
    ggplot(df, aes(x = time/30.4, y = PROT, color = Scenario)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = setNames(
        sapply(unique(df$Scenario), function(s) SCENARIOS[[s]]$color),
        unique(df$Scenario))) +
      geom_hline(yintercept = 0.3, linetype="dashed", color=COLS$remission,  size=0.8) +
      geom_hline(yintercept = 3.5, linetype="dashed", color=COLS$warning, size=0.8) +
      labs(title = "Proteinuria by Treatment Scenario",
           x = "Time (months)", y = "Proteinuria (g/day)") +
      THEME_MN
  })

  output$sc_egfr_plot <- renderPlot({
    df <- scenario_sims()
    if (is.null(df)) return(NULL)
    ggplot(df, aes(x = time/30.4, y = EGFR, color = Scenario)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = setNames(
        sapply(unique(df$Scenario), function(s) SCENARIOS[[s]]$color),
        unique(df$Scenario))) +
      labs(title = "eGFR by Treatment Scenario",
           x = "Time (months)", y = "eGFR (mL/min/1.73m²)") +
      THEME_MN
  })

  output$sc_pla2r_plot <- renderPlot({
    df <- scenario_sims()
    if (is.null(df)) return(NULL)
    ggplot(df, aes(x = time/30.4, y = PLA2R, color = Scenario)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = setNames(
        sapply(unique(df$Scenario), function(s) SCENARIOS[[s]]$color),
        unique(df$Scenario))) +
      geom_hline(yintercept = 14, linetype="dashed", color=COLS$remission, size=0.8) +
      labs(title = "Anti-PLA2R1 Titer by Scenario",
           x = "Time (months)", y = "Anti-PLA2R1 (RU/mL)") +
      THEME_MN
  })

  output$sc_summary_table <- renderDT({
    df <- scenario_sims()
    if (is.null(df)) return(datatable(data.frame(Message="No scenarios selected")))
    bp <- input$bprot
    summary <- df %>%
      group_by(Scenario) %>%
      summarise(
        `Final Proteinuria (g/day)` = round(tail(PROT, 1), 2),
        `CR Rate (%)` = round(mean(PROT < 0.3) * 100, 1),
        `PR Rate (%)` = round(mean(PROT < 3.5 & PROT < bp*0.5) * 100, 1),
        `eGFR at endpoint` = round(tail(EGFR, 1), 1),
        `Anti-PLA2R1 final (RU/mL)` = round(tail(PLA2R, 1), 1),
        .groups = "drop"
      ) %>%
      arrange(`Final Proteinuria (g/day)`)

    datatable(summary, options = list(dom = "t", pageLength = 10),
              rownames = FALSE) %>%
      formatStyle("CR Rate (%)",
                  background = styleColorBar(c(0,100), COLS$remission),
                  backgroundSize = "80% 70%", backgroundRepeat = "no-repeat",
                  backgroundPosition = "center")
  })

  # ==========================================================================
  # TAB 6: BIOMARKER DASHBOARD
  # ==========================================================================
  bm_sim <- eventReactive(input$run_bm, {
    t_days <- seq(0, 24 * 30.4, by = 1)
    p <- c(base_params(),
           list(USE_RTX = as.numeric(input$bm_rtx),
                USE_TAC = as.numeric(input$bm_tac),
                USE_CYC = as.numeric(input$bm_cyc),
                USE_ACE = 0))
    run_simulation(p, t_days)
  }, ignoreNULL = FALSE)

  # Traffic lights
  traffic_light <- function(val, g_thresh, y_thresh, higher_worse = TRUE) {
    if (higher_worse) {
      cls <- if (val <= g_thresh) "traffic-green" else if (val <= y_thresh) "traffic-yellow" else "traffic-red"
    } else {
      cls <- if (val >= g_thresh) "traffic-green" else if (val >= y_thresh) "traffic-yellow" else "traffic-red"
    }
    div(class = cls, "●")
  }

  output$bm_pla2r_light  <- renderUI({ traffic_light(input$bpla2r, 14, 150, higher_worse=TRUE) })
  output$bm_prot_light   <- renderUI({ traffic_light(input$bprot,   0.3, 3.5, higher_worse=TRUE) })
  output$bm_alb_light    <- renderUI({ traffic_light(input$balb,    3.5, 2.5, higher_worse=FALSE) })
  output$bm_egfr_light   <- renderUI({ traffic_light(input$begfr,   60,  30,  higher_worse=FALSE) })
  output$bm_compl_light  <- renderUI({
    # Complement activity: 1 = baseline (bad), <0.2 = remission (good)
    traffic_light(1 - input$bpla2r/500, 0.6, 0.3, higher_worse=FALSE)
  })
  output$bm_thrombo_light <- renderUI({
    score <- (input$bprot/15)*2 + (4.5-input$balb)/3 + (1 - pmin(input$begfr/120,1))*1
    traffic_light(score, 1.5, 2.5, higher_worse=TRUE)
  })

  output$bm_pla2r_plot <- renderPlot({
    df <- bm_sim()
    ggplot(df, aes(x = time/30.4, y = PLA2R)) +
      geom_line(color = COLS$worsening, size = 1.2) +
      geom_ribbon(aes(ymin=0, ymax=PLA2R), fill=COLS$worsening, alpha=0.1) +
      geom_hline(yintercept=14,  linetype="dashed", color=COLS$remission, size=0.8) +
      geom_hline(yintercept=150, linetype="dashed", color=COLS$warning, size=0.7) +
      labs(title="Anti-PLA2R1 Titer (24-month trajectory)",
           x="Time (months)", y="Anti-PLA2R1 (RU/mL)") +
      THEME_MN
  })

  output$bm_compl_plot <- renderPlot({
    df <- bm_sim()
    df_long <- df %>%
      transmute(time, `MAC Activity` = COMPL, `IgG4 Level` = IGG4) %>%
      pivot_longer(-time, names_to="Marker", values_to="Value")
    ggplot(df_long, aes(x=time/30.4, y=Value, color=Marker)) +
      geom_line(size=1.1) +
      scale_color_manual(values=c("MAC Activity"=COLS$warning, "IgG4 Level"="#984EA3")) +
      labs(title="Complement MAC & IgG4", x="Time (months)", y="Relative level") +
      THEME_MN
  })

  output$bm_renal_plot <- renderPlot({
    df <- bm_sim()
    # NGAL and KIM-1 as proxies driven by proteinuria
    df <- df %>%
      mutate(NGAL_rel = 1 + (PROT / input$bprot - 0.5) * 2,
             KIM1_rel = 1 + (PROT / input$bprot - 0.5) * 1.5)
    df_long <- df %>%
      transmute(time,
                `NGAL (relative)` = pmax(NGAL_rel, 0.1),
                `KIM-1 (relative)` = pmax(KIM1_rel, 0.1)) %>%
      pivot_longer(-time, names_to="Marker", values_to="Value")
    ggplot(df_long, aes(x=time/30.4, y=Value, color=Marker)) +
      geom_line(size=1.1) +
      scale_color_manual(values=c("NGAL (relative)"="#1F78B4","KIM-1 (relative)"="#33A02C")) +
      geom_hline(yintercept=1.0, linetype="dashed", color=COLS$neutral, size=0.7) +
      labs(title="Tubular Injury Biomarkers (NGAL, KIM-1)",
           x="Time (months)", y="Relative level (1 = normal)") +
      THEME_MN
  })

  output$bm_lipid_plot <- renderPlot({
    df <- bm_sim()
    df <- df %>%
      mutate(LDL_rel   = LIPID,
             TChol_rel = LIPID * 1.15,
             TG_rel    = 1 + (4.2 - ALB) * 0.25)
    df_long <- df %>%
      transmute(time,
                `LDL (relative)` = LDL_rel,
                `Total Cholesterol` = TChol_rel,
                `Triglycerides` = TG_rel) %>%
      pivot_longer(-time, names_to="Lipid", values_to="Value")
    ggplot(df_long, aes(x=time/30.4, y=Value, color=Lipid)) +
      geom_line(size=1.1) +
      scale_color_manual(values=c("LDL (relative)"="#E41A1C",
                                  "Total Cholesterol"="#FF7F00",
                                  "Triglycerides"="#FFFF33")) +
      geom_hline(yintercept=1.0, linetype="dashed", color=COLS$neutral, size=0.7) +
      labs(title="Serum Lipids (nephrotic syndrome-related dyslipidemia)",
           x="Time (months)", y="Relative level (1 = normal upper limit)") +
      THEME_MN
  })

  output$thrombo_score_ui <- renderUI({
    alb <- input$balb; prot <- input$bprot; egfr <- input$begfr
    score <- 0
    factors <- character(0)
    if (alb < 2.5)  { score <- score + 2; factors <- c(factors, "Serum albumin < 2.5 g/dL (+2)") }
    if (prot > 8)   { score <- score + 2; factors <- c(factors, "Proteinuria > 8 g/day (+2)") }
    if (egfr < 30)  { score <- score + 1; factors <- c(factors, "eGFR < 30 (+1)") }
    if (input$bpla2r > 200) { score <- score + 1; factors <- c(factors, "High anti-PLA2R1 (+1)") }
    risk_cat <- if (score >= 4) "HIGH" else if (score >= 2) "MODERATE" else "LOW"
    risk_col <- if (score >= 4) "#D6604D" else if (score >= 2) "#F4A736" else "#1A9850"
    tagList(
      fluidRow(
        column(4,
          h4("Thrombosis Risk Score"),
          div(style = paste0("font-size:36px; font-weight:bold; color:", risk_col, ";"), score, "/6"),
          div(class="remission-badge", style=paste0("background:", risk_col, "; color:white;"), risk_cat),
          br(),
          helpText("Score ≥ 4: Consider prophylactic anticoagulation (heparin/warfarin)")
        ),
        column(8,
          h5("Risk Factors Present:"),
          if (length(factors) > 0) tags$ul(lapply(factors, tags$li))
          else tags$p("No major thrombosis risk factors identified.", style="color:green;")
        )
      )
    )
  })

  output$thrombo_gauge <- renderPlot({
    alb <- input$balb; prot <- input$bprot; egfr <- input$begfr
    score <- min(6, (alb < 2.5)*2 + (prot > 8)*2 + (egfr < 30)*1 + (input$bpla2r > 200)*1)
    gauge_df <- data.frame(
      Component = c("Score", "Remaining"),
      Value     = c(score, 6 - score)
    )
    fill_col <- if (score >= 4) COLS$worsening else if (score >= 2) COLS$warning else COLS$remission
    ggplot(gauge_df, aes(x = "", y = Value, fill = Component)) +
      geom_col(width = 0.5, color = "white") +
      scale_fill_manual(values = c("Score" = fill_col, "Remaining" = "#EEEEEE")) +
      coord_polar(theta = "y", start = pi) +
      labs(title = paste("Thrombosis Risk:", round(score/6*100), "%")) +
      theme_void() +
      theme(legend.position = "none",
            plot.title = element_text(hjust=0.5, face="bold", size=16, color=fill_col))
  })
}

# =============================================================================
# LAUNCH
# =============================================================================
shinyApp(ui = ui, server = server)
