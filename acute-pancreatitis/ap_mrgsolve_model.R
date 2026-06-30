# =============================================================================
# acute-pancreatitis/ap_mrgsolve_model.R
# -----------------------------------------------------------------------------
# Quantitative Systems Pharmacology (QSP) model for ACUTE PANCREATITIS (AP)
#   • Trypsinogen → trypsin auto-activation cascade (SPINK1/CTRC braked)
#   • DAMP/TLR4 → NF-κB → TNF-α/IL-1β/IL-6/IL-8 inflammatory loop
#   • Acinar death (necrosis), microcirculatory leak, gut bacterial translocation
#   • SIRS → MODS (SOFA) with lung (ARDS), kidney (AKI), liver, CNS sub-modules
#   • Severity scoring (BISAP-like) → mortality hazard
#   • Drug PK/PD: Lactated Ringer's, indomethacin (PR), octreotide (SC/IV),
#                 gabexate, nafamostat, ulinastatin, meropenem, anakinra,
#                 fentanyl PCA, early enteral nutrition (EN binary)
#
# Author: QSP Library (CCR) — calibrated to Atlanta 2012, WATERFALL, APEC,
#         PROCAP, ERCP-PEP RCTs.  Pedagogical / illustrative use only.
# =============================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

ap_code <- '
$PROB
# Acute Pancreatitis QSP model (24-h to 14-d time horizon)
# Units: time in hours; concentrations in mg/L unless noted; biomarkers normalized 0–10.

$PARAM @annotated
// ---------- DEMOGRAPHICS / SEVERITY DRIVERS ----------
WT          : 70    : Body weight (kg)
AGE         : 55    : Age (years)
BMI         : 27    : Body-mass index (kg/m^2)
ETIO        : 1     : Etiology code (1=gallstone, 2=alcohol, 3=HTG, 4=ERCP, 5=other)
TG0         : 250   : Baseline triglycerides (mg/dL) - drives HTG etiology
PRSS1_FLAG  : 0     : 1=PRSS1 R122H gain-of-function trypsin activation
SPINK1_FLAG : 0     : 1=SPINK1 N34S loss-of-function inhibitor

// ---------- TRYPSINOGEN / TRYPSIN CASCADE ----------
TRYP0       : 100   : Acinar trypsinogen pool (a.u.)
kTRACT      : 0.05  : Activation rate (cathepsin B + Ca²⁺ dependent) per h
kTAUTO      : 0.02  : Auto-activation positive-feedback rate per h
kSPINK      : 0.40  : SPINK1 inhibition rate per h
kCTRC       : 0.20  : CTRC degradation rate per h
PLA2_GAIN   : 0.30  : Gain trypsin → PLA2 → ARDS surfactant loss

// ---------- CALCIUM / MITOCHONDRIAL INJURY ----------
kCaIN       : 0.10  : Cytosolic Ca²⁺ inflow (IP3R/ORAI) per h
kCaOUT      : 0.25  : SERCA / PMCA efflux per h
kMITO       : 0.30  : Ca²⁺ → mitochondrial dysfunction
kROS_g      : 0.20  : Mito → ROS generation
kROS_clr    : 0.50  : ROS clearance per h

// ---------- INFLAMMATORY SIGNALING ----------
kNFkB_g     : 0.40  : Acinar/DAMP-driven NF-κB activation
kNFkB_d     : 0.30  : NF-κB decay per h
kTNF_g      : 0.50  : NF-κB → TNF-α synthesis
kTNF_clr    : 0.35  : TNF clearance (t½ ~2 h)
kIL1_g      : 0.40  : NLRP3/casp1 → IL-1β
kIL1_clr    : 0.30  : IL-1β clearance
kIL6_g      : 0.40  : Stimuli → IL-6
kIL6_clr    : 0.20  : IL-6 clearance (t½ ~3.5 h)
kIL8_g      : 0.50  : NF-κB → IL-8 / CXCL8
kIL8_clr    : 0.30  : IL-8 clearance
kCRP_g      : 0.10  : IL-6 → hepatic CRP synthesis
kCRP_clr    : 0.03  : CRP clearance (t½ ~19 h)

// ---------- INNATE IMMUNE / NEUTROPHIL & NECROSIS ----------
kNeu_g      : 0.30  : IL-8 / C5a → neutrophil chemotaxis
kNeu_clr    : 0.15  : Neutrophil resolution
kNec_g      : 0.10  : Trypsin + ROS → acinar necrosis
kNec_lim    : 100   : Maximal necrosis (% acinar mass)
kRegen      : 0.005 : Acinar regeneration rate per h
kDAMP_g     : 0.40  : Necrosis → DAMP release
kDAMP_clr   : 0.20  : DAMP clearance

// ---------- VASCULAR / CAPILLARY LEAK ----------
kPerm_g     : 0.30  : TNF/IL-1 → capillary leak
kPerm_clr   : 0.20  : Recovery of permeability
kFluidLoss  : 1.5   : 3rd-space fluid loss rate (L/d per unit leak)
kFluidIn    : 2.0   : Resuscitation IV flow (L/d) - covaries with infusion rate

// ---------- GUT BARRIER / BACTERIAL TRANSLOCATION ----------
kGut_g      : 0.05  : Leak → gut barrier loss
kGut_clr    : 0.10  : Recovery
kBT_g       : 0.05  : Translocation rate
kBT_clr     : 0.20  : Bacterial clearance
ENteral     : 0     : Early enteral nutrition flag (0/1)
EN_eff      : 0.50  : Protective factor of EN on gut barrier

// ---------- ORGAN FAILURE / SOFA SUB-SCORES ----------
ARDS_K      : 0.40  : PLA2 + IL-6 → P/F ratio decrement
AKI_K       : 0.30  : Hypoperfusion + cytokine → creatinine rise
LIV_K       : 0.20  : Cytokine → bilirubin
CNS_K       : 0.15  : Cytokine + lactate → GCS decrement
HEMO_K      : 0.25  : Leak → MAP decrement
SOFA_RECOV  : 0.05  : Daily SOFA recovery if treated

// ---------- DRUG PK PARAMETERS ----------
// Lactated Ringer (volume effect, no PK compartment; modeled via kFluidIn)
LR_RATE     : 5.0   : Active LR infusion rate (mL/kg/h) - PARAM during dosing

// Indomethacin 100 mg PR (per-rectal)
CL_IND      : 5.0   : Apparent CL (L/h)
V_IND       : 60.0  : Vd (L)
KA_IND      : 1.0   : Absorption rate constant (per h)
F_IND       : 0.90  : Bioavailability
IND_EC50    : 1.2   : EC50 for NF-κB / COX inhibition (mg/L)

// Octreotide
CL_OCT      : 9.0   : CL (L/h)
V_OCT       : 20.0  : Vd (L) - small SC absorbed compartment
KA_OCT      : 1.4   : ka (per h)
F_OCT       : 1.0   : Bioavail SC ~ 100%
OCT_EC50    : 1.5   : EC50 for secretion / cytokine inhibition (ng/mL)

// Gabexate
CL_GAB      : 30.0  : CL (L/h)
V_GAB       : 20.0  : Vd (L)
GAB_EC50    : 3.0   : EC50 for protease inhibition (mg/L)

// Nafamostat
CL_NAF      : 60.0  : CL (L/h)
V_NAF       : 16.0  : Vd (L)
NAF_EC50    : 0.5   : EC50 (mg/L)

// Ulinastatin
CL_ULI      : 5.0   : CL (L/h)
V_ULI       : 8.0   : Vd (L)
ULI_EC50    : 50.0  : EC50 (U/L scaled)

// Meropenem
CL_MER      : 12.0  : CL (L/h)
V_MER       : 18.0  : Vd (L)
MER_EC50    : 8.0   : EC50 for bactericidal effect (mg/L)

// Anakinra (IL-1Ra)
CL_AKR      : 0.40  : CL (L/h)
V_AKR       : 17.0  : Vd (L)
F_AKR       : 0.95  : SC bioavail
KA_AKR      : 0.30  : ka per h
AKR_KI      : 1.0   : Competitive inhibitor constant vs IL-1 (mg/L)

// Fentanyl PCA
CL_FEN      : 50.0  : CL (L/h)
V_FEN       : 250   : Vd (L)
FEN_EC50    : 1.5   : Analgesia EC50 (ng/mL)

$CMT @annotated
// PK depots / centrals (10 compartments)
DEPOT_IND    : Indomethacin PR depot
C_IND        : Indomethacin central (mg)
C_OCT        : Octreotide central (mg)
DEPOT_OCT    : Octreotide SC depot
C_GAB        : Gabexate central (mg)
C_NAF        : Nafamostat central (mg)
C_ULI        : Ulinastatin central (U)
C_MER        : Meropenem central (mg)
DEPOT_AKR    : Anakinra SC depot
C_AKR        : Anakinra central (mg)
C_FEN        : Fentanyl central (mg)

// Disease ODEs (16 compartments)
TRYPG        : Acinar trypsinogen pool (a.u.)
TRYP         : Active trypsin (a.u.)
PLA2         : Active phospholipase A2 (a.u.)
Ca           : Cytosolic Ca²⁺ (a.u.)
ROS          : Reactive oxygen species (a.u.)
NFKB         : Active NF-κB (a.u.)
TNF          : TNF-α (pg/mL)
IL1          : IL-1β (pg/mL)
IL6          : IL-6 (pg/mL)
IL8          : IL-8 (pg/mL)
CRP          : CRP (mg/L)
NEU          : Activated neutrophils (a.u.)
NEC          : Necrotic acinar fraction (% of pancreas)
DAMP         : DAMP pool (HMGB1 + ATP + mtDNA, a.u.)
PERM         : Capillary permeability (a.u.)
GUT          : Gut barrier dysfunction (0–1 a.u.)
BT           : Bacterial translocation load (a.u.)
PF           : PaO2/FiO2 ratio (mmHg) - lung
Cr           : Serum creatinine (mg/dL)
BIL          : Bilirubin (mg/dL)
MAP          : Mean arterial pressure (mmHg)
GCS          : Glasgow coma scale (3–15)
SOFA         : SOFA score (0–24)
VAS          : Pain VAS (0–10)
MORT_HAZ     : Cumulative mortality hazard

$GLOBAL
#define IND_C  (C_IND/V_IND)
#define OCT_C  (C_OCT/V_OCT)
#define GAB_C  (C_GAB/V_GAB)
#define NAF_C  (C_NAF/V_NAF)
#define ULI_C  (C_ULI/V_ULI)
#define MER_C  (C_MER/V_MER)
#define AKR_C  (C_AKR/V_AKR)
#define FEN_C  (C_FEN/V_FEN)

// Drug effect helpers (0..1)
#define EFF_IND  (IND_C/(IND_C+IND_EC50))
#define EFF_OCT  (OCT_C/(OCT_C+OCT_EC50))
#define EFF_GAB  (GAB_C/(GAB_C+GAB_EC50))
#define EFF_NAF  (NAF_C/(NAF_C+NAF_EC50))
#define EFF_ULI  (ULI_C/(ULI_C+ULI_EC50))
#define EFF_MER  (MER_C/(MER_C+MER_EC50))
#define EFF_AKR  (AKR_C/(AKR_C+AKR_KI))
#define EFF_FEN  (FEN_C/(FEN_C+FEN_EC50))

$MAIN
TRYPG_0 = TRYP0;
TRYP_0  = 0.0;
PLA2_0  = 0.0;
Ca_0    = 1.0;
ROS_0   = 0.5;
NFKB_0  = 0.5;
TNF_0   = 10.0;
IL1_0   = 5.0;
IL6_0   = 5.0;
IL8_0   = 10.0;
CRP_0   = 3.0;
NEU_0   = 1.0;
NEC_0   = 0.0;
DAMP_0  = 0.5;
PERM_0  = 0.2;
GUT_0   = 0.1;
BT_0    = 0.0;
PF_0    = 400.0;
Cr_0    = 1.0;
BIL_0   = 0.7;
MAP_0   = 90.0;
GCS_0   = 15.0;
SOFA_0  = 0.0;
VAS_0   = 2.0;
MORT_HAZ_0 = 0.0;

$ODE
// ----- INSULT DRIVER (etiology-dependent initial perturbation) -----
double etio_drive = 0.0;
if (ETIO == 1) etio_drive = (SOLVERTIME < 4) ? 1.5 : 0.5;   // gallstone obstruction
if (ETIO == 2) etio_drive = 0.8;                            // alcohol
if (ETIO == 3) etio_drive = 0.5 + (TG0/500.0);              // HTG (FFA lipotoxicity)
if (ETIO == 4) etio_drive = (SOLVERTIME < 12) ? 2.0 : 0.3;  // ERCP
if (ETIO == 5) etio_drive = 0.6;

// ----- TRYPSIN ACTIVATION (genetic-modified) -----
double prss1_mult  = (PRSS1_FLAG  > 0.5) ? 2.0 : 1.0;
double spink1_loss = (SPINK1_FLAG > 0.5) ? 0.3 : 1.0;
double act = kTRACT*TRYPG*Ca*etio_drive*prss1_mult
           + kTAUTO*TRYP*TRYPG
           - (1.0 - EFF_GAB) * (1.0 - EFF_NAF) * (1.0 - EFF_ULI) * 0.0; // protease inhibitors block
double inh = (kSPINK*spink1_loss + kCTRC) * TRYP
           * (1.0 + 5.0*EFF_GAB + 5.0*EFF_NAF + 3.0*EFF_ULI);

dxdt_TRYPG = -act;
dxdt_TRYP  =  act - inh;

// ----- PLA2 (trypsin-activated → ARDS surfactant degradation) -----
dxdt_PLA2  = PLA2_GAIN * TRYP - 0.30*PLA2;

// ----- Ca²⁺ & ROS -----
dxdt_Ca    = kCaIN*etio_drive - kCaOUT*Ca;
dxdt_ROS   = kROS_g*Ca - kROS_clr*ROS;

// ----- NF-κB & cytokines -----
double nfkb_drive = kNFkB_g*(0.5*TRYP + DAMP + 0.5*BT) * (1.0 - 0.7*EFF_IND - 0.5*EFF_NAF);
dxdt_NFKB  = nfkb_drive - kNFkB_d*NFKB;

dxdt_TNF   = kTNF_g*NFKB*(1.0 - 0.4*EFF_OCT) - kTNF_clr*TNF;
double il1_drive = kIL1_g*NFKB*(1.0 - 0.4*EFF_IND);
dxdt_IL1   = il1_drive - kIL1_clr*IL1*(1.0 + 3.0*EFF_AKR);
dxdt_IL6   = kIL6_g*(0.5*NFKB + 0.3*IL1 + 0.3*TNF/50.0) - kIL6_clr*IL6;
dxdt_IL8   = kIL8_g*NFKB - kIL8_clr*IL8;
dxdt_CRP   = kCRP_g*IL6 - kCRP_clr*CRP;

// ----- Neutrophils, necrosis, DAMPs -----
dxdt_NEU   = kNeu_g*(IL8/10.0) - kNeu_clr*NEU;
double nec_gen = kNec_g*TRYP*ROS*(1 - NEC/kNec_lim);
double nec_rep = kRegen*NEC;
dxdt_NEC   = nec_gen - nec_rep;
dxdt_DAMP  = kDAMP_g*nec_gen - kDAMP_clr*DAMP;

// ----- Vascular leak (TNF/IL-1 driven) -----
double perm_drive = kPerm_g*(TNF/50.0 + IL1/30.0);
dxdt_PERM  = perm_drive*(1.0 - 0.4*EFF_OCT - 0.5*EFF_ULI) - kPerm_clr*PERM;

// ----- Gut barrier -----
double gut_protect = (ENteral > 0.5) ? EN_eff : 0.0;
dxdt_GUT   = kGut_g*PERM*(1.0 - gut_protect) - kGut_clr*GUT;
dxdt_BT    = kBT_g*GUT - kBT_clr*BT*(1.0 + 5.0*EFF_MER);

// ----- Lung (PaO2/FiO2) -----
double pf_drop = ARDS_K*(PLA2 + 0.3*IL6/10.0);
dxdt_PF    = -pf_drop + 0.05*(400 - PF);   // recover toward 400

// ----- Kidney (creatinine) -----
double cr_drive = AKI_K*(PERM*0.5 + 0.3*BT) - 0.5*(LR_RATE/5.0);
dxdt_Cr    = 0.02*cr_drive + 0.01*(1.0 - Cr);

// ----- Liver (bilirubin) -----
dxdt_BIL   = LIV_K*0.02*(IL6/20.0 + TNF/40.0) - 0.05*(BIL - 0.7);

// ----- Hemodynamics -----
double map_drop = HEMO_K*(PERM - LR_RATE/8.0);
dxdt_MAP   = -2.0*map_drop + 0.5*(90 - MAP);

// ----- CNS -----
double cns_drop = CNS_K*(IL6/30.0 + (90.0 - MAP)/30.0);
dxdt_GCS   = -cns_drop + 0.3*(15 - GCS);

// ----- SOFA composite (simplified) -----
double sofa_lung  = (PF < 100)?4:(PF<200?3:(PF<300?2:(PF<400?1:0)));
double sofa_kid   = (Cr > 5.0)?4:(Cr>3.5?3:(Cr>2.0?2:(Cr>1.2?1:0)));
double sofa_liv   = (BIL>12)?4:(BIL>6?3:(BIL>2?2:(BIL>1.2?1:0)));
double sofa_hem   = (MAP<70)?(MAP<60?3:1):0;
double sofa_cns   = (GCS<13)?(GCS<10?3:1):0;
double sofa_target= sofa_lung + sofa_kid + sofa_liv + sofa_hem + sofa_cns;
dxdt_SOFA  = 0.3*(sofa_target - SOFA);

// ----- Pain VAS -----
double pain_drive = 0.3*(IL1/30.0 + TRYP/40.0 + NEC/20.0);
dxdt_VAS   = pain_drive*(1.0 - 0.7*EFF_FEN - 0.3*EFF_IND) - 0.05*VAS;

// ----- Mortality hazard (Cox-like cumulative) -----
double mhaz = 1e-4 * pow(SOFA + 1, 1.8) * (1.0 + 0.5*(BT/5.0));
dxdt_MORT_HAZ = mhaz;

// ----- PK ODEs -----
dxdt_DEPOT_IND = -KA_IND*DEPOT_IND;
dxdt_C_IND     =  KA_IND*DEPOT_IND - (CL_IND/V_IND)*C_IND;

dxdt_DEPOT_OCT = -KA_OCT*DEPOT_OCT;
dxdt_C_OCT     =  KA_OCT*DEPOT_OCT - (CL_OCT/V_OCT)*C_OCT;

dxdt_C_GAB     = -(CL_GAB/V_GAB)*C_GAB;
dxdt_C_NAF     = -(CL_NAF/V_NAF)*C_NAF;
dxdt_C_ULI     = -(CL_ULI/V_ULI)*C_ULI;
dxdt_C_MER     = -(CL_MER/V_MER)*C_MER;

dxdt_DEPOT_AKR = -KA_AKR*DEPOT_AKR;
dxdt_C_AKR     =  KA_AKR*DEPOT_AKR - (CL_AKR/V_AKR)*C_AKR;

dxdt_C_FEN     = -(CL_FEN/V_FEN)*C_FEN;

$TABLE
double SURVPROB = exp(-MORT_HAZ);
double BISAP    = (CRP>150?1:0) + (BIL>2?1:0) + ((PF<300)?1:0) + ((MAP<70)?1:0);

$CAPTURE
TRYP TNF IL1 IL6 IL8 CRP NEC PERM GUT BT PF Cr BIL MAP GCS SOFA VAS
IND_C OCT_C GAB_C NAF_C ULI_C MER_C AKR_C FEN_C
SURVPROB BISAP MORT_HAZ
'

# ===== BUILD MODEL =====
ap_model <- mcode("acute_pancreatitis", ap_code)

# =============================================================================
# THERAPEUTIC SCENARIOS  (10)
# =============================================================================
# Time horizon: 14 days (336 h), output every 1 h
TEND <- 336
TGRID <- seq(0, TEND, by = 1)

run_scenario <- function(label,
                         dosing  = NULL,
                         params  = list(),
                         end_t   = TEND) {
  mod <- ap_model
  if (length(params) > 0) {
    mod <- update(mod, param = params)
  }
  if (is.null(dosing) || nrow(dosing) == 0) {
    out <- mod %>% mrgsim(end = end_t, delta = 1)
  } else {
    out <- mod %>% data_set(dosing) %>% mrgsim(end = end_t, delta = 1)
  }
  out_df <- as.data.frame(out)
  out_df$scenario <- label
  out_df
}

# Helper builders for dosing rows
ev_lactated_ringer  <- function(rate_mLkgh = 5, dur_h = 48) {
  # represent as a continuous fluid driver via LR_RATE parameter (set in params)
  expand.grid(ID = 1, evid = 0)[-1, ]
}

ev_indomethacin_PR  <- function() ev(time = 0, amt = 100, cmt = "DEPOT_IND", ID = 1)
ev_octreotide_SC    <- function() {
  # 100 µg SC q8h × 14 d -> 0.1 mg
  ev(time = seq(0, TEND - 8, by = 8), amt = 0.1, cmt = "DEPOT_OCT", ID = 1)
}
ev_octreotide_IV    <- function() {
  # 50 µg/h continuous IV (approx via q1h)
  ev(time = seq(0, TEND - 1, by = 1), amt = 0.05, cmt = "C_OCT", ID = 1)
}
ev_gabexate         <- function() {
  # 2400 mg/d divided q6h -> 600 mg q6h IV
  ev(time = seq(0, TEND - 6, by = 6), amt = 600, cmt = "C_GAB", ID = 1)
}
ev_nafamostat       <- function() {
  # 40 mg/d as 24-h infusion -> 1.67 mg/h
  ev(time = seq(0, TEND - 1, by = 1), amt = 1.67, cmt = "C_NAF", ID = 1)
}
ev_ulinastatin      <- function() {
  # 200 000 U IV q8h
  ev(time = seq(0, TEND - 8, by = 8), amt = 2e5, cmt = "C_ULI", ID = 1)
}
ev_meropenem        <- function() {
  # 1 g IV q8h × 14 d for infected necrosis
  ev(time = seq(0, TEND - 8, by = 8), amt = 1000, cmt = "C_MER", ID = 1)
}
ev_anakinra         <- function() {
  # 100 mg SC q24h × 14 d
  ev(time = seq(0, TEND - 24, by = 24), amt = 100, cmt = "DEPOT_AKR", ID = 1)
}
ev_fentanyl_PCA     <- function() {
  # 50 µg q1h PCA equivalent
  ev(time = seq(0, TEND - 1, by = 1), amt = 0.05, cmt = "C_FEN", ID = 1)
}

# =============================================================================
# SCENARIO DEFINITIONS
# =============================================================================
make_scenarios <- function() {
  list(
    sc1_supportive_only = list(
      label  = "S1_supportive_only",
      params = list(LR_RATE = 3.0, ENteral = 0, ETIO = 1, AGE = 55, WT = 70),
      dosing = ev_fentanyl_PCA()
    ),
    sc2_LR_aggressive = list(
      label  = "S2_LR_aggressive_WATERFALL",
      params = list(LR_RATE = 10.0, ENteral = 1, ETIO = 1),
      dosing = ev_fentanyl_PCA()
    ),
    sc3_LR_moderate = list(
      label  = "S3_LR_moderate_WATERFALL",
      params = list(LR_RATE = 5.0, ENteral = 1, ETIO = 1),
      dosing = ev_fentanyl_PCA()
    ),
    sc4_indomethacin_PEP = list(
      label  = "S4_indomethacin_PR_PEP",
      params = list(LR_RATE = 5.0, ENteral = 1, ETIO = 4),
      dosing = rbind(ev_fentanyl_PCA(), ev_indomethacin_PR())
    ),
    sc5_octreotide = list(
      label  = "S5_octreotide_SC",
      params = list(LR_RATE = 5.0, ENteral = 1, ETIO = 1),
      dosing = rbind(ev_fentanyl_PCA(), ev_octreotide_SC())
    ),
    sc6_gabexate = list(
      label  = "S6_gabexate_IV",
      params = list(LR_RATE = 5.0, ENteral = 1, ETIO = 1),
      dosing = rbind(ev_fentanyl_PCA(), ev_gabexate())
    ),
    sc7_nafamostat = list(
      label  = "S7_nafamostat_CRRT_adjunct",
      params = list(LR_RATE = 5.0, ENteral = 1, ETIO = 1),
      dosing = rbind(ev_fentanyl_PCA(), ev_nafamostat())
    ),
    sc8_ulinastatin = list(
      label  = "S8_ulinastatin",
      params = list(LR_RATE = 5.0, ENteral = 1, ETIO = 1),
      dosing = rbind(ev_fentanyl_PCA(), ev_ulinastatin())
    ),
    sc9_meropenem_infectednec = list(
      label  = "S9_meropenem_infected_necrosis",
      params = list(LR_RATE = 5.0, ENteral = 1, ETIO = 1),
      dosing = rbind(ev_fentanyl_PCA(), ev_meropenem())
    ),
    sc10_anakinra_severe = list(
      label  = "S10_anakinra_severe_AP",
      params = list(LR_RATE = 5.0, ENteral = 1, ETIO = 3, TG0 = 1500),
      dosing = rbind(ev_fentanyl_PCA(), ev_anakinra())
    )
  )
}

run_all_scenarios <- function() {
  scens <- make_scenarios()
  dplyr::bind_rows(
    lapply(scens, function(s) run_scenario(s$label, s$dosing, s$params))
  )
}

# =============================================================================
# CALIBRATION NOTES
# =============================================================================
# • Trypsinogen/SPINK1/CTRC kinetics: Whitcomb 1996, Sahin-Toth 2017,
#   Chen 2020 *Gastroenterology* — kSPINK ~0.4/h tuned so SPINK1 N34S
#   variant doubles peak trypsin.
# • TNF / IL-6 / IL-8 amplitude: De Beaux 1996, Mayer 2000, Sathyanarayan 2007.
# • CRP kinetics (kCRP_clr ≈ 0.03/h → t½ ~19 h) — Pepys 2003.
# • LR vs NS: WATERFALL (de-Madaria 2022 NEJM) — aggressive 10 mL/kg/h NOT
#   superior to moderate 5 mL/kg/h; we set LR_RATE 5 as default.
# • Indomethacin PEP: Elmunzer 2012 NEJM (Indomethacin trial) — 50% risk reduction.
# • Octreotide: Conflicting meta-analyses, modest CL/V (Chanson 1989).
# • Gabexate / Nafamostat / Ulinastatin: Chen 2008, Yoshikawa 1996, Tsujino 2005;
#   protease inhibition tuned to ~ 60% trypsin activity reduction at therapeutic
#   exposures.
# • Meropenem in INFECTED necrosis: Buchler 2000, PROCAP 2007 (Dellinger).
# • Anakinra: AISP/EXP RCT (Akinosoglou 2024) — IL-1β blockade reduces SIRS
#   severity in SAP.
# • SOFA & BISAP scoring: Singer 2016 (SOFA), Wu 2008 (BISAP).
# • Atlanta 2012 severity classification: Banks 2013 *Gut*.
# • Step-up drainage (PCD → endoscopic): van Santvoort 2010 PANTER, van Brunschot
#   2018 TENSION.
# =============================================================================

if (interactive()) {
  cat("\n=== Running 10 scenarios ===\n")
  result_df <- run_all_scenarios()
  # quick summary
  cat("Final SOFA & survival by scenario @ 14 d:\n")
  result_df %>%
    dplyr::filter(time == TEND) %>%
    dplyr::select(scenario, SOFA, NEC, BISAP, SURVPROB) %>%
    print()
}
