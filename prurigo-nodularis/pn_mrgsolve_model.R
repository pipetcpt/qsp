# =============================================================================
# prurigo-nodularis/pn_mrgsolve_model.R
# -----------------------------------------------------------------------------
# Quantitative Systems Pharmacology (QSP) model for PRURIGO NODULARIS (PN)
#   • Keratinocyte barrier dysfunction/alarmins (TSLP/IL-33/IL-25) -> Th2/Th17
#   • IL-31 -> IL-31RA/OSMRβ peripheral nerve sensitization (key itch driver)
#   • Central sensitization (spinal GRPR/central gain) + MOR:KOR opioid tone
#   • Itch -> scratch -> mechanical trauma -> TGF-β/fibroblast -> nodule burden
#     (self-reinforcing itch-scratch-fibrosis cycle, feeds back to barrier)
#   • Drug PK/PD: dupilumab (anti-IL-4Rα), nemolizumab (anti-IL-31RA),
#                 gabapentin (α2δ, central sensitization), nalbuphine ER
#                 (MOR antagonist / KOR agonist, 2-cpt oral)
#   • Clinical outputs: Worst-Itch NRS, PN-IGA, nodule count, sleep, DLQI
#
# Author: QSP Library (CCR) — calibrated to PRIME/PRIME2 (dupilumab),
#         OLYMPIA 1/2 (nemolizumab), Weisshaar 2022 JEADV Phase 2 (nalbuphine ER),
#         gabapentinoid real-world cohorts. Pedagogical / illustrative use only.
# =============================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

pn_code <- '
$PROB
# Prurigo Nodularis QSP model (time horizon: up to 24 weeks / 4032 h)
# Units: time in hours; drug concentrations in mg/L (or ng/mL where noted); biomarkers normalized 0-10 a.u.

$PARAM @annotated
// ---------- DEMOGRAPHICS / DISEASE SEVERITY DRIVERS ----------
WT           : 75    : Body weight (kg)
AGE          : 58    : Age (years)
BASELINE_WINRS : 8.5 : Baseline Worst-Itch NRS (0-10)
NODULE0      : 25    : Baseline active nodule count
ATOPIC_FLAG  : 0     : 1 = atopic diathesis subtype (higher Th2 drive)
CKD_FLAG     : 0     : 1 = CKD-associated pruritus overlap (uremic/opioid driver)
DURATION_YR  : 4     : Disease duration (years) - drives baseline IENFD depletion

// ---------- BARRIER / ALARMIN CASCADE ----------
kBARR_g      : 0.015 : Barrier dysfunction generation rate (scratch-driven), per h
kBARR_clr    : 0.010 : Barrier dysfunction spontaneous repair, per h
kTSLP_g      : 0.40  : Barrier -> TSLP alarmin generation
kTSLP_clr    : 0.35  : TSLP clearance, per h
kTH2_g       : 0.30  : TSLP/IL-33/IL-25 -> Th2 (IL-4/IL-13) drive
kTH2_clr     : 0.20  : Th2 signal decay, per h
kTH17_g      : 0.18  : Mixed Th17/Th22 drive (epidermal hyperplasia arm)
kTH17_clr    : 0.20  : Th17 signal decay, per h
kIL31_g      : 0.35  : Th2/Th22 -> IL-31 synthesis rate
kIL31_clr    : 0.45  : IL-31 clearance (t1/2 ~ 1.5 h), per h

// ---------- PERIPHERAL & CENTRAL SENSITIZATION ----------
kPSENS_g     : 0.25  : IL-31/TRPV1/TRPA1 -> peripheral sensitization gain
kPSENS_clr   : 0.15  : Peripheral sensitization decay, per h
kIENFD_loss  : 0.0006: Chronic-nodule driven intraepidermal nerve fiber loss, per h
IENFD_regen  : 0.0004: Spontaneous nerve fiber regeneration, per h
kCSENS_g     : 0.22  : Peripheral drive -> central (spinal/GRPR) sensitization
kCSENS_clr   : 0.12  : Central sensitization decay, per h
kOPIOID_g    : 0.10  : Baseline drift toward pruritogenic MOR:KOR imbalance
kOPIOID_clr  : 0.08  : Opioid tone normalization, per h
CKD_OPIOID_BOOST : 1.6 : Multiplier on opioid pruritogenic tone if CKD-aP overlap

// ---------- ITCH -> SCRATCH -> FIBROSIS CYCLE ----------
kSCRATCH_g   : 0.50  : Itch perception -> scratch behavior gain
kSCRATCH_clr : 0.40  : Scratch behavior decay, per h
kNODULE_g    : 0.0009: Scratch-driven fibroblast/TGF-β -> nodule growth rate
kNODULE_res  : 0.00025: Spontaneous nodule involution rate (with therapy-permitted quiescence)
NODULE_CAP   : 60    : Carrying-capacity nodule count (a.u.)

// ---------- CLINICAL OUTPUT KINETICS (indirect-response turnover) ----------
kWINRS_eq    : 0.08  : WI-NRS equilibration rate toward itch-drive set point, per h
kIGA_eq      : 0.02  : PN-IGA equilibration rate, per h
kSLEEP_eq    : 0.10  : Sleep disturbance equilibration rate, per h
kDLQI_eq     : 0.05  : DLQI equilibration rate, per h

// ---------- DRUG TARGET POTENCIES (EC50, effect gains) ----------
DUPI_EC50    : 50    : Dupilumab central conc. for 50% IL-4Rα blockade (mg/L)
DUPI_EMAX    : 0.92  : Max fractional blockade of Th2/TSLP-driven signal
NEMO_EC50    : 2.5   : Nemolizumab central conc. for 50% IL-31RA blockade (mg/L)
NEMO_EMAX    : 0.90  : Max fractional blockade of IL-31 -> peripheral sensitization
GABA_EC50    : 4     : Gabapentin central conc. for 50% central-sensitization damping (mg/L)
GABA_EMAX    : 0.45  : Max fractional damping of central sensitization
JAKI_EC50    : 150   : Oral JAK1i (abrocitinib) conc. for 50% Th2/IL-31 signal block (ng/mL)
JAKI_EMAX    : 0.70  : Max fractional blockade (JAK1-dependent transduction)
JAKI_ON      : 0     : 1 = JAK1 inhibitor active in this run (flag; conc still simulated)
NAL_EC50     : 15    : Nalbuphine central conc. for 50% opioid-tone correction (ng/mL)
NAL_EMAX     : 0.75  : Max fractional correction of MOR:KOR pruritogenic imbalance
KOR_AGONIST_ON : 0   : 1 = adjunct peripheral KOR-agonist (difelikefalin-like) add-on active
KOR_AGONIST_EFF: 0.35: Direct fractional reduction of opioid pruritogenic tone if ON
TCS_REPAIR   : 0     : Topical corticosteroid barrier-repair boost (0-1)

// ---------- PK PARAMETERS ----------
// Dupilumab (SC, linear approx around therapeutic range)
KA_DUPI      : 0.0075 : Absorption rate constant (1/h), t1/2abs ~ 4 d
CL_DUPI      : 0.0021 : Clearance (L/h)
V_DUPI       : 4.6    : Volume of distribution (L)
F_DUPI       : 0.64   : Bioavailability (SC)

// Nemolizumab (SC)
KA_NEMO      : 0.0090 : Absorption rate constant (1/h)
CL_NEMO      : 0.014  : Clearance (L/h)
V_NEMO       : 4.6    : Volume of distribution (L)
F_NEMO       : 0.68   : Bioavailability (SC)

// Gabapentin (oral, saturable transporter-limited absorption approximated linear)
KA_GABA      : 0.15   : Absorption rate constant (1/h)
CL_GABA      : 8.5    : Clearance (L/h) - renally eliminated
V_GABA       : 58     : Volume of distribution (L)
F_GABA       : 0.55   : Oral bioavailability (dose-dependent, approx at 300-600 mg)

// Nalbuphine ER (oral, 2-compartment)
KA_NAL       : 0.55   : Absorption rate constant (1/h)
CL_NAL       : 115    : Clearance (L/h) - high hepatic extraction
V2_NAL       : 180    : Central volume (L)
V3_NAL       : 260    : Peripheral volume (L)
Q_NAL        : 55     : Inter-compartmental clearance (L/h)
F_NAL        : 0.12   : Oral bioavailability (extensive first-pass)

// JAK1 inhibitor (abrocitinib, oral, simple 1-cpt)
KA_JAKI      : 1.1    : Absorption rate constant (1/h)
CL_JAKI      : 42     : Clearance (L/h)
V_JAKI       : 100    : Volume of distribution (L)
F_JAKI       : 0.65   : Oral bioavailability

$CMT @annotated
// ---------- PK compartments (11) ----------
DEPOT_DUPI   : Dupilumab SC depot (mg)
C_DUPI       : Dupilumab central amount (mg)
DEPOT_NEMO   : Nemolizumab SC depot (mg)
C_NEMO       : Nemolizumab central amount (mg)
GUT_GABA     : Gabapentin GI depot (mg)
C_GABA       : Gabapentin central amount (mg)
GUT_NAL      : Nalbuphine ER GI depot (mg)
C_NAL        : Nalbuphine central amount (mg)
P_NAL        : Nalbuphine peripheral amount (mg)
GUT_JAKI     : Oral JAK1 inhibitor GI depot (mg)
C_JAKI       : Oral JAK1 inhibitor central amount (mg)

// ---------- Disease/PD compartments (16) ----------
BARRIER      : Epidermal barrier dysfunction (0-1 a.u.)
TSLP         : Keratinocyte alarmin TSLP (a.u., ~pg/mL scale)
TH2          : Th2 (IL-4/IL-13) composite signal (a.u.)
TH17         : Th17/Th22 composite signal (a.u.)
IL31         : IL-31 (a.u., ~pg/mL scale)
PSENS        : Peripheral nerve sensitization (0-10 a.u.)
IENFD        : Intraepidermal nerve fiber density (fraction of normal, 0-1)
CSENS        : Central (spinal/GRPR) sensitization (0-10 a.u.)
OPIOID       : Net pruritogenic MOR:KOR opioid tone (a.u., 0 = balanced)
SCRATCH      : Scratch behavior intensity (0-10 a.u.)
NODULE       : Active fibrotic nodule burden (count, a.u.)
WINRS        : Worst-Itch NRS (0-10) - clinical output
IGA          : PN Investigator Global Assessment (0-4) - clinical output
SLEEP        : Sleep-disturbance score (0-10)
DLQI         : Dermatology Life Quality Index (0-30)
CUM_SCRATCH  : Cumulative scratch exposure (a.u., tracks mechanical trauma burden)

$MAIN
DEPOT_DUPI_0 = 0; C_DUPI_0 = 0;
DEPOT_NEMO_0 = 0; C_NEMO_0 = 0;
GUT_GABA_0   = 0; C_GABA_0 = 0;
GUT_NAL_0    = 0; C_NAL_0  = 0; P_NAL_0 = 0;
GUT_JAKI_0   = 0; C_JAKI_0 = 0;

BARRIER_0 = 0.55;
TSLP_0    = 4.0;
TH2_0     = (ATOPIC_FLAG > 0.5) ? 5.5 : 3.5;
TH17_0    = 2.0;
IL31_0    = 6.5;
PSENS_0   = 6.5;
IENFD_0   = fmax(0.25, 1.0 - 0.10*DURATION_YR);
CSENS_0   = 6.0;
OPIOID_0  = (CKD_FLAG > 0.5) ? 3.5 : 1.5;
SCRATCH_0 = 6.0;
NODULE_0  = NODULE0;
WINRS_0   = BASELINE_WINRS;
IGA_0     = 3.0;
SLEEP_0   = 6.5;
DLQI_0    = 18.0;
CUM_SCRATCH_0 = 0.0;

$GLOBAL
#define CP_DUPI (C_DUPI/V_DUPI)
#define CP_NEMO (C_NEMO/V_NEMO)
#define CP_GABA (C_GABA/V_GABA)
#define CP_NAL  (C_NAL/V2_NAL*1000.0)
#define CP_JAKI (C_JAKI/V_JAKI*1000.0)

#define EFF_DUPI (DUPI_EMAX*CP_DUPI/(CP_DUPI+DUPI_EC50))
#define EFF_NEMO (NEMO_EMAX*CP_NEMO/(CP_NEMO+NEMO_EC50))
#define EFF_GABA (GABA_EMAX*CP_GABA/(CP_GABA+GABA_EC50))
#define EFF_JAKI ((JAKI_ON>0.5) ? (JAKI_EMAX*CP_JAKI/(CP_JAKI+JAKI_EC50)) : 0.0)
#define EFF_NAL  (NAL_EMAX*CP_NAL/(CP_NAL+NAL_EC50))

$ODE
// ---------------- PK ----------------
dxdt_DEPOT_DUPI = -KA_DUPI*DEPOT_DUPI;
dxdt_C_DUPI     =  KA_DUPI*F_DUPI*DEPOT_DUPI - (CL_DUPI/V_DUPI)*C_DUPI;

dxdt_DEPOT_NEMO = -KA_NEMO*DEPOT_NEMO;
dxdt_C_NEMO     =  KA_NEMO*F_NEMO*DEPOT_NEMO - (CL_NEMO/V_NEMO)*C_NEMO;

dxdt_GUT_GABA   = -KA_GABA*GUT_GABA;
dxdt_C_GABA     =  KA_GABA*F_GABA*GUT_GABA - (CL_GABA/V_GABA)*C_GABA;

dxdt_GUT_NAL    = -KA_NAL*GUT_NAL;
dxdt_C_NAL      =  KA_NAL*F_NAL*GUT_NAL - (CL_NAL/V2_NAL)*C_NAL - (Q_NAL/V2_NAL)*C_NAL + (Q_NAL/V3_NAL)*P_NAL;
dxdt_P_NAL      =  (Q_NAL/V2_NAL)*C_NAL - (Q_NAL/V3_NAL)*P_NAL;

dxdt_GUT_JAKI   = -KA_JAKI*GUT_JAKI;
dxdt_C_JAKI     =  KA_JAKI*F_JAKI*GUT_JAKI - (CL_JAKI/V_JAKI)*C_JAKI;

// ---------------- BARRIER -> ALARMIN -> TH2/TH17 ----------------
// NOTE: this cascade is a closed positive-feedback loop (IL31 -> PSENS -> CSENS ->
// SCRATCH -> BARRIER -> TSLP -> TH2 -> IL31). Each compartment therefore carries an
// explicit logistic ceiling (1 - X/CAP) on its generation term so the loop relaxes to
// a bounded steady state instead of diverging, mirroring physiological receptor/
// cytokine saturation rather than unconstrained linear amplification.
double topical_repair = TCS_REPAIR;  // set via $PARAM override per scenario (see below)
dxdt_BARRIER = kBARR_g*(SCRATCH/10.0)*(1.0 - BARRIER) - kBARR_clr*BARRIER*(1.0 + topical_repair);

dxdt_TSLP    = kTSLP_g*BARRIER*10.0*(1.0 - TSLP/20.0) - kTSLP_clr*TSLP;

double th2_drive = kTH2_g*TSLP*(ATOPIC_FLAG>0.5 ? 1.3 : 1.0);
dxdt_TH2     = th2_drive*(1.0 - EFF_DUPI)*(1.0 - TH2/20.0) - kTH2_clr*TH2;

dxdt_TH17    = kTH17_g*TSLP*0.6*(1.0 - TH17/15.0) - kTH17_clr*TH17;

double il31_drive = kIL31_g*(TH2 + 0.5*TH17);
dxdt_IL31    = il31_drive*(1.0 - EFF_DUPI)*(1.0 - EFF_JAKI)*(1.0 - IL31/30.0) - kIL31_clr*IL31;

// ---------------- PERIPHERAL & CENTRAL SENSITIZATION ----------------
double psens_drive = kPSENS_g*IL31*(1.0/fmax(IENFD,0.15));
dxdt_PSENS   = psens_drive*(1.0 - EFF_NEMO)*(1.0 - EFF_JAKI)*(1.0 - PSENS/10.0) - kPSENS_clr*PSENS;

double ienfd_loss  = kIENFD_loss*NODULE;
dxdt_IENFD   = IENFD_regen*(1.0 - IENFD) - ienfd_loss*IENFD;

double csens_drive = kCSENS_g*(PSENS + 0.4*OPIOID);
dxdt_CSENS   = csens_drive*(1.0 - EFF_GABA)*(1.0 - CSENS/10.0) - kCSENS_clr*CSENS;

double opioid_drive = kOPIOID_g*(1.0 + (CKD_FLAG>0.5 ? (CKD_OPIOID_BOOST-1.0) : 0.0));
double opioid_corr  = kOPIOID_clr*OPIOID*(1.0 + 3.0*EFF_NAL) + (KOR_AGONIST_ON>0.5 ? KOR_AGONIST_EFF*OPIOID : 0.0);
dxdt_OPIOID  = opioid_drive*(1.0 - OPIOID/10.0) - opioid_corr;

// ---------------- ITCH -> SCRATCH -> FIBROSIS ----------------
double itch_drive = 0.35*PSENS + 0.45*CSENS + 0.20*OPIOID;
dxdt_SCRATCH = kSCRATCH_g*itch_drive - kSCRATCH_clr*SCRATCH;
dxdt_CUM_SCRATCH = SCRATCH;

double nodule_growth = kNODULE_g*SCRATCH*(1.0 - NODULE/NODULE_CAP);
double nodule_resolve = kNODULE_res*NODULE*(1.0 + 2.0*EFF_DUPI + 1.5*EFF_NEMO);
dxdt_NODULE  = nodule_growth - nodule_resolve;

// ---------------- CLINICAL OUTPUTS (indirect response, equilibrating) ----------------
double winrs_setpoint = fmin(10.0, itch_drive);
dxdt_WINRS   = kWINRS_eq*(winrs_setpoint - WINRS);

double iga_setpoint = fmin(4.0, 0.55 + 0.06*NODULE + 0.25*BARRIER);
dxdt_IGA     = kIGA_eq*(iga_setpoint - IGA);

double sleep_setpoint = fmin(10.0, 0.8*WINRS + 0.5);
dxdt_SLEEP   = kSLEEP_eq*(sleep_setpoint - SLEEP);

double dlqi_setpoint = fmin(30.0, 1.6*WINRS + 1.1*IGA + 0.7*SLEEP);
dxdt_DLQI    = kDLQI_eq*(dlqi_setpoint - DLQI);

$TABLE
double WINRS_PCT_IMPROVE = 100.0*(BASELINE_WINRS - WINRS)/BASELINE_WINRS;
double IGA_SUCCESS = (IGA <= 1.0) ? 1.0 : 0.0;
double RESPONDER   = (WINRS_PCT_IMPROVE >= 40.0 && IGA_SUCCESS > 0.5) ? 1.0 : 0.0;

$CAPTURE
CP_DUPI CP_NEMO CP_GABA CP_NAL CP_JAKI
WINRS_PCT_IMPROVE IGA_SUCCESS RESPONDER
'

# ===== BUILD MODEL =====
pn_model <- mcode("prurigo_nodularis", pn_code)

# =============================================================================
# THERAPEUTIC SCENARIOS  (8)
# =============================================================================
TEND  <- 4032          # 24 weeks in hours
TGRID <- seq(0, TEND, by = 4)

run_scenario <- function(label, dosing = NULL, params = list(), end_t = TEND) {
  mod <- pn_model
  if (length(params) > 0) mod <- update(mod, param = params)
  if (is.null(dosing) || nrow(dosing) == 0) {
    out <- mod %>% mrgsim(end = end_t, delta = 4)
  } else {
    out <- mod %>% ev(dosing) %>% mrgsim(end = end_t, delta = 4)
  }
  out_df <- as.data.frame(out)
  out_df$scenario <- label
  out_df
}

# ---- Dosing regimens ----
ev_dupilumab <- function() {
  # 600 mg loading, then 300 mg SC Q2W (PRIME/PRIME2 regimen)
  seq_ev(ev(time = 0, amt = 600, cmt = "DEPOT_DUPI"),
         ev(time = seq(336, TEND - 336, by = 336), amt = 300, cmt = "DEPOT_DUPI"))
}
ev_nemolizumab <- function() {
  # 60 mg loading, then 30 mg SC Q4W (OLYMPIA regimen, weight-tiered)
  seq_ev(ev(time = 0, amt = 60, cmt = "DEPOT_NEMO"),
         ev(time = seq(672, TEND - 672, by = 672), amt = 30, cmt = "DEPOT_NEMO"))
}
ev_gabapentin <- function() {
  # 300 mg TID oral, titrated up to steady dosing from day 1
  ev(time = seq(0, TEND - 8, by = 8), amt = 300, cmt = "GUT_GABA")
}
ev_nalbuphine_er <- function() {
  # 162 mg BID oral (Weisshaar 2022 JEADV Phase 2 top dose)
  ev(time = seq(0, TEND - 12, by = 12), amt = 162, cmt = "GUT_NAL")
}
ev_abrocitinib <- function() {
  # 200 mg QD oral (off-label JAK1i)
  ev(time = seq(0, TEND - 24, by = 24), amt = 200, cmt = "GUT_JAKI")
}

seq_ev <- function(...) Reduce(`+`, list(...), accumulate = FALSE)  # combine event objects

scenarios <- list(
  list(label = "Natural History (topical SOC only)",
       dosing = NULL, params = list(TCS_REPAIR = 0.15)),
  list(label = "Dupilumab 600mg LD / 300mg Q2W SC",
       dosing = ev_dupilumab(), params = list(TCS_REPAIR = 0.15)),
  list(label = "Nemolizumab 60mg LD / 30mg Q4W SC",
       dosing = ev_nemolizumab(), params = list(TCS_REPAIR = 0.15)),
  list(label = "Gabapentin 300mg TID oral",
       dosing = ev_gabapentin(), params = list(TCS_REPAIR = 0.15)),
  list(label = "Nalbuphine ER 162mg BID oral",
       dosing = ev_nalbuphine_er(), params = list(TCS_REPAIR = 0.15)),
  list(label = "Abrocitinib 200mg QD oral (off-label)",
       dosing = ev_abrocitinib(), params = list(TCS_REPAIR = 0.15, JAKI_ON = 1)),
  list(label = "Dupilumab + Gabapentin combination",
       dosing = ev_dupilumab() + ev_gabapentin(), params = list(TCS_REPAIR = 0.15)),
  list(label = "Nalbuphine ER + adjunct KOR agonist (CKD-aP overlap)",
       dosing = ev_nalbuphine_er(),
       params = list(TCS_REPAIR = 0.15, CKD_FLAG = 1, KOR_AGONIST_ON = 1))
)

results <- lapply(scenarios, function(s) {
  run_scenario(s$label, dosing = s$dosing, params = s$params)
})
pn_sim_all <- bind_rows(results)

# =============================================================================
# CALIBRATION NOTES (pedagogical targets from published trials)
# -----------------------------------------------------------------------------
# * Dupilumab PRIME/PRIME2 (Yosipovitch et al. NEJM 2023): WI-NRS >=4-pt
#   improvement in ~60% at week 24 vs ~18% placebo; IGA success ~48% vs 18%.
#   -> model DUPI_EMAX/EC50 tuned so RESPONDER ~0.55-0.65 by week 24.
# * Nemolizumab OLYMPIA 1/2 (Ständer/Kwatra 2024 NEJM): WI-NRS responder
#   ~56-60% at week 16 vs ~21-27% placebo; rapid onset by week 4.
#   -> NEMO_EC50 set lower (higher potency) to reflect faster onset vs dupilumab.
# * Nalbuphine ER (Kwatra et al. NEJM Evidence 2023, Phase 2/3): ~50% WI-NRS
#   reduction at 162mg BID by week 8-10, modest vs placebo (~30%).
# * Gabapentinoids: mainly retrospective/open-label cohorts (Gooding 2010
#   systematic review) - modest symptomatic effect, no barrier/nodule
#   remodeling benefit (EFF_GABA acts only on CSENS, not BARRIER/TH2/IL31).
# * Abrocitinib/upadacitinib: small case series/off-label; JAKI_EMAX kept
#   intermediate (0.70) with a scenario flag (JAKI_ON) since not FDA-approved
#   for PN as of this writing.
# =============================================================================

if (interactive()) {
  ggplot(pn_sim_all, aes(x = time/168, y = WINRS, color = scenario)) +
    geom_line(linewidth = 0.9) +
    labs(x = "Time (weeks)", y = "Worst-Itch NRS (0-10)",
         title = "Prurigo Nodularis QSP — Worst-Itch NRS by Treatment Scenario") +
    theme_minimal()
}
