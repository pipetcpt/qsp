# ============================================================================
# Postherpetic Neuralgia (PHN) — Quantitative Systems Pharmacology Model
# mrgsolve PK/PD platform
#
# Compartments (24 ODEs)
#   PK (per drug, simplified one/two-compartment models)
#     1.  GBP_GUT, GBP_CENT     gabapentin (with LAT1-saturable absorption)
#     2.  PGB_CENT              pregabalin (1-cpt linear)
#     3.  AMI_CENT, AMI_NOR     amitriptyline → nortriptyline (CYP2D6 metabolite)
#     4.  DLX_CENT              duloxetine (1-cpt)
#     5.  LIDO_SKIN             topical lidocaine 5% (dermal depot)
#     6.  CAP_SKIN              capsaicin 8% (single-application, residual depot)
#     7.  VAL_CENT              valaciclovir / aciclovir (active species)
#     8.  RZV_VAC               recombinant zoster vaccine antigen dose (gE+AS01B)
#     9.  TRA_CENT              tramadol  (μ + SNRI, simplified)
#   PD / disease biology
#    10.  VZV_LOAD              viral DNA copies in skin/ganglion (zoster phase)
#    11.  GANG_INJ              dorsal-root-ganglion injury index (0-1)
#    12.  IENF                  intra-epidermal nerve fiber density (% of healthy)
#    13.  NAV_ACT               peripheral NaV1.7/1.8 ectopic activity
#    14.  CSEN                  central sensitization / spinal LTP index (0-1)
#    15.  MICROG                spinal microglial activation (0-1)
#    16.  KCC2                  KCC2 chloride extrusion (% of baseline)
#    17.  NGF_SKIN              skin NGF concentration (relative)
#    18.  NMDA_TONE             NMDA-driven WDR tone (0-1)
#    19.  CMI                   anti-VZV cell-mediated immunity (IFN-γ ELISpot units)
#    20.  PAIN                  NRS 0-10 daily pain score
#    21.  ALLO                  dynamic mechanical allodynia (0-10)
#    22.  SLEEP                 sleep disturbance (0-10, PSQI-like)
#    23.  MOOD                  PHQ-9 / HADS-D anxiety-depression composite
#    24.  AE_SED                cumulative sedation / dizziness AE index
#
# Calibration anchors (literature targets):
#   * RZV (Shingrix) HZ efficacy ≥97% (ZOE-50/70 — Lal 2015, Cunningham 2016)
#   * Valaciclovir 1 g tid ↓ acute pain time-to-resolution by ~30% (Tyring 2000)
#   * Pregabalin 300-600 mg/d → ~50% pain reduction in ~35-50% (Dworkin 2003, 2018 AAN)
#   * Gabapentin 1800-3600 mg/d → ~30% NRS reduction (Rice 2001; Wallace 2010)
#   * Amitriptyline 25-100 mg HS → NNT ~2.7 for 50% pain (Watson 1992)
#   * Duloxetine 60 mg → ~25-30% pain reduction (off-label, NeuPSIG)
#   * Lidocaine 5% patch → 30-50% allodynia reduction, low systemic Css (Argoff 2000)
#   * Capsaicin 8% single 60-min application → 30% reduction sustained 12 wk
#     (PACE: Backonja 2008; STRIDE: Webster 2017)
#   * Tramadol 200-400 mg → 30% NRS reduction (Boureau 2003)
#
# Notes
#   - Doses encoded in CMT 1 (GBP_GUT) etc via amt=, cmt= in events
#   - PK parameters from labels & published popPK studies (citations in references)
#   - Disease ODEs are deliberately phenomenological/semimechanistic —
#     intended for in-silico scenario comparison and trial design, not regulatory.
# ============================================================================

library(mrgsolve)

code <- '
$PROB
  PHN QSP model v1.0 — PK/PD of analgesics + RZV prophylaxis
  Disease compartments: VZV → ganglion injury → peripheral & central sensitization
                       → multimodal pain phenotype → sleep/mood/QoL endpoints.

$PARAM @annotated
// ===== Patient covariates =====
AGE     : 72   : age (yr)
WT      : 70   : weight (kg)
CRCL    : 60   : creatinine clearance (mL/min)
CYP2D6  : 1    : CYP2D6 activity score (0.5 PM, 1 IM, 2 EM, 3 UM)
SEXF    : 0    : female (1) / male (0)
DM      : 0    : diabetes (0/1)
IS      : 0    : immunosuppressed (0/1)
HZ_OPHTH: 0    : ophthalmic HZ (0/1)

// ===== Gabapentin PK (saturable LAT1 absorption, renal CL) =====
GBP_KA   : 0.40 : /h
GBP_VMAX : 230  : LAT1 Vmax (mg/h)  — produces ~60% F at 900 mg, ~30% at 3600 mg/d
GBP_KM   : 110  : LAT1 Km (mg)
GBP_V    : 70   : L
GBP_CLr0 : 13   : L/h reference at CRCL 100
GBP_F    : 1    : flag for absorption

// ===== Pregabalin PK (linear, mostly renal) =====
PGB_KA   : 1.4  : /h
PGB_V    : 60   : L
PGB_F    : 0.90 : oral bioavailability
PGB_CLr0 : 7.0  : L/h at CRCL 100

// ===== Amitriptyline → Nortriptyline (CYP2D6) =====
AMI_KA   : 1.2  : /h
AMI_V    : 700  : L (high Vss)
AMI_F    : 0.50 :
AMI_CL   : 50   : L/h hepatic (parent disappearance)
NOR_FM   : 0.40 : fraction converted to nortriptyline
NOR_V    : 1300 : L
NOR_CL   : 30   : L/h

// ===== Duloxetine =====
DLX_KA   : 0.75 : /h
DLX_V    : 1620 : L
DLX_F    : 0.50 :
DLX_CL   : 100  : L/h

// ===== Topical lidocaine 5% patch (12-h on/off cycle) =====
LIDO_KSKIN : 0.05 : /h (release into local pool)
LIDO_KEL   : 0.40 : /h (local clearance)
LIDO_DOSE  : 350  : mg (typical 1-patch up-to-3 daily, but kept low Css)

// ===== Capsaicin 8% (single 60-min application; persistent receptor block) =====
CAP_KDEC   : 0.005 : /day  — slow defunctionalization decay
CAP_E50    : 0.5   :        normalized topical dose for 50% TRPV1 block

// ===== Valaciclovir (acyclovir) =====
VAL_KA   : 1.5  : /h
VAL_V    : 50   : L
VAL_F    : 0.55 :
VAL_CL   : 19   : L/h (renal)

// ===== Tramadol (combined μ + SNRI) =====
TRA_KA   : 1.0  : /h
TRA_V    : 200  : L
TRA_F    : 0.70 :
TRA_CL   : 40   : L/h

// ===== Recombinant Zoster Vaccine (RZV) =====
RZV_KDEG : 0.005 : /day (gE antigen depot decay)
RZV_EFF  : 0.97  : maximum efficacy ceiling on HZ probability

// ===== Disease dynamics =====
KVZ_IN     : 1.0   : VZV reactivation rate constant (zoster onset)
KVZ_OUT    : 0.10  : intrinsic VZV clearance per day
ACV_EMAX   : 0.85  : aciclovir Emax on VZV
ACV_EC50   : 1.2   : mg/L active

KGANG_IN   : 0.4   : ganglion injury accrual from VZV load
KGANG_REC  : 0.005 : intrinsic ganglion recovery
NGF_IN     : 0.6   : NGF response to ganglion injury
NGF_OUT    : 0.05  : daily NGF decay
KNAV_IN    : 1.0   : NaV ectopic activity from NGF + injury
KNAV_OUT   : 0.10  : intrinsic NaV decay
KCS_IN     : 1.0   : central sensitization input from NaV
KCS_OUT    : 0.05  : intrinsic CSEN decay
KMG_IN     : 0.7   : microglia activation by CSEN
KMG_OUT    : 0.05  :
KCC2_BASE  : 1.0   :
KKCC2_LOSS : 0.5   : loss per unit microglia/d
KKCC2_REC  : 0.05  :
KNMDA_IN   : 0.8   :
KNMDA_OUT  : 0.10  :
KIENF_LOSS : 0.005 : IENF/day loss from injury (slow)
KIENF_REC  : 0.0015 : IENF/day intrinsic recovery
CMI_BASE   : 1.0   : ELISpot units baseline
KCMI_LOSS  : 0.02  : senescence loss /yr (folded into IC age effect)
KCMI_REC   : 0.01  : recovery rate

// ===== Drug-target effects (Emax) =====
PGB_EMAX   : 0.55 : pregabalin maximal CSEN reduction
PGB_EC50   : 3.5  : mg/L  (free PGB)
GBP_EMAX   : 0.45 : gabapentin maximal CSEN reduction
GBP_EC50   : 5.0  : mg/L
AMI_EMAX   : 0.50 : amitriptyline+NOR descending facilitation
AMI_EC50   : 0.05 : mg/L (sum AMI+NOR)
DLX_EMAX   : 0.40 : duloxetine
DLX_EC50   : 0.05 :
TRA_EMAX   : 0.45 : tramadol composite
TRA_EC50   : 0.15 : mg/L
LIDO_EMAX  : 0.50 : local NaV block (on skin pool)
LIDO_EC50  : 1.5  : mg/g skin (very rough)
CAP_EMAX   : 0.70 : TRPV1 defunctionalization (Pain reduction)
KMICRO_INH : 0.20 : microglial dampening by anti-inflammatories / steroids

// ===== Pain composite weights =====
W_NRS_NAV  : 4.0  : peripheral component (allodynia, sharp)
W_NRS_CSEN : 4.0  : central component (burning, hyperalgesia)
W_NRS_NMDA : 2.0  : NMDA wind-up component
B_PAIN_REC : 0.05 : intrinsic relaxation toward homeostasis
W_ALLO     : 6.0  : allodynia weight on NAV/IENF
W_SLEEP    : 0.6  :
W_MOOD     : 0.4  :
KAE_SED    : 0.10 : sedation AE accrual per (PGB+GBP+AMI+TRA)
KAE_DECAY  : 0.20 : daily sedation decay

$CMT  @annotated
GBP_GUT  : gabapentin gut
GBP_CENT : gabapentin central
PGB_CENT : pregabalin central
AMI_CENT : amitriptyline central
NOR_CENT : nortriptyline central
DLX_CENT : duloxetine central
LIDO_SKIN: lidocaine skin pool
CAP_SKIN : capsaicin skin pool
VAL_CENT : aciclovir central
TRA_CENT : tramadol central
RZV_VAC  : RZV antigen depot

VZV_LOAD : VZV load (relative)
GANG_INJ : ganglion injury (0-1)
IENF     : IENF density (%)
NAV_ACT  : peripheral NaV ectopic activity
CSEN     : central sensitization
MICROG   : microglial activation
KCC2     : KCC2 function (% baseline)
NGF_SKIN : skin NGF (rel)
NMDA_TONE: NMDA tone
CMI      : anti-VZV CMI

PAIN     : NRS 0-10
ALLO     : allodynia 0-10
SLEEP    : sleep disturbance 0-10
MOOD     : mood 0-10
AE_SED   : sedation/AE index

$GLOBAL
  // CrCl-dependent CL scaling helpers
  double crcl_factor(double crcl){ return (crcl/100.0); }

$MAIN
  // Initialize disease state at simulation start
  // For acute herpes-zoster scenario:
  //   VZV_LOAD_0 = 1, GANG_INJ_0 = 0, IENF_0 = 100, etc.
  GBP_GUT_0  = 0;
  GBP_CENT_0 = 0;
  PGB_CENT_0 = 0;
  AMI_CENT_0 = 0;
  NOR_CENT_0 = 0;
  DLX_CENT_0 = 0;
  LIDO_SKIN_0= 0;
  CAP_SKIN_0 = 0;
  VAL_CENT_0 = 0;
  TRA_CENT_0 = 0;
  RZV_VAC_0  = 0;

  VZV_LOAD_0 = 1.0;
  GANG_INJ_0 = 0.0;
  IENF_0     = 100.0;
  NAV_ACT_0  = 0.05;
  CSEN_0     = 0.05;
  MICROG_0   = 0.05;
  KCC2_0     = 1.0;
  NGF_SKIN_0 = 0.1;
  NMDA_TONE_0= 0.05;
  CMI_0      = 1.0 * (1.0 - 0.005*(AGE-50));   // immunosenescence

  PAIN_0     = 0;
  ALLO_0     = 0;
  SLEEP_0    = 1;
  MOOD_0     = 1;
  AE_SED_0   = 0;

$ODE
  // ---- PK ----
  // Gabapentin saturable absorption: dGBP_GUT/dt = -Vmax*A/(Km+A)
  double absRate_GBP = GBP_VMAX*GBP_GUT/(GBP_KM + GBP_GUT + 1e-9);
  dxdt_GBP_GUT  = -absRate_GBP;
  double GBP_CL = GBP_CLr0 * crcl_factor(CRCL);
  dxdt_GBP_CENT =  absRate_GBP - (GBP_CL/GBP_V)*GBP_CENT;

  double PGB_CL = PGB_CLr0 * crcl_factor(CRCL);
  dxdt_PGB_CENT = -(PGB_CL/PGB_V)*PGB_CENT;

  // Amitriptyline → Nortriptyline (CYP2D6 scaling)
  double cypf = CYP2D6/1.0;
  double ami_k10 = (AMI_CL*cypf)/AMI_V;
  dxdt_AMI_CENT = - ami_k10*AMI_CENT;
  dxdt_NOR_CENT =   NOR_FM*ami_k10*AMI_CENT*(AMI_V/NOR_V) - (NOR_CL/NOR_V)*NOR_CENT;

  dxdt_DLX_CENT = -(DLX_CL/DLX_V)*DLX_CENT;

  // Lidocaine skin depot (local; minimal systemic spillover modeled implicit)
  dxdt_LIDO_SKIN = - LIDO_KEL*LIDO_SKIN;

  // Capsaicin depot — slow decay
  dxdt_CAP_SKIN  = - CAP_KDEC*CAP_SKIN;

  // Valaciclovir → aciclovir composite
  dxdt_VAL_CENT = -(VAL_CL/VAL_V)*VAL_CENT;

  // Tramadol
  dxdt_TRA_CENT = -(TRA_CL/TRA_V)*TRA_CENT;

  // RZV antigen depot
  dxdt_RZV_VAC  = - RZV_KDEG * RZV_VAC;

  // ---- Concentrations ----
  double Cp_GBP = GBP_CENT/GBP_V;
  double Cp_PGB = PGB_CENT/PGB_V;
  double Cp_AMI = AMI_CENT/AMI_V;
  double Cp_NOR = NOR_CENT/NOR_V;
  double Cp_DLX = DLX_CENT/DLX_V;
  double Cp_VAL = VAL_CENT/VAL_V;
  double Cp_TRA = TRA_CENT/TRA_V;
  double Cs_LIDO = LIDO_SKIN/100.0;  // arbitrary normalized skin concentration
  double Cs_CAP  = CAP_SKIN /1.0;

  double Edrug_VZV  = ACV_EMAX*Cp_VAL/(ACV_EC50 + Cp_VAL);
  double Edrug_PGB  = PGB_EMAX*Cp_PGB/(PGB_EC50 + Cp_PGB);
  double Edrug_GBP  = GBP_EMAX*Cp_GBP/(GBP_EC50 + Cp_GBP);
  double Edrug_AMI  = AMI_EMAX*(Cp_AMI + Cp_NOR)/(AMI_EC50 + Cp_AMI + Cp_NOR);
  double Edrug_DLX  = DLX_EMAX*Cp_DLX/(DLX_EC50 + Cp_DLX);
  double Edrug_TRA  = TRA_EMAX*Cp_TRA/(TRA_EC50 + Cp_TRA);
  double Edrug_LIDO = LIDO_EMAX*Cs_LIDO/(LIDO_EC50 + Cs_LIDO);
  double Edrug_CAP  = CAP_EMAX*Cs_CAP/(CAP_E50  + Cs_CAP);

  // ---- Disease ODEs ----
  // Vaccine effect on probability of zoster (interpreted as reduced reactivation drive)
  double vac_protection = RZV_EFF*RZV_VAC/(1.0 + RZV_VAC);
  // Force VZV_LOAD growth gating by (1 - vac_protection)
  double vzv_growth = KVZ_IN*(1.0 - vac_protection)*(1.0 + 0.5*IS) -
                      KVZ_OUT*VZV_LOAD*(1.0 + 4*Edrug_VZV);
  dxdt_VZV_LOAD = vzv_growth;

  // Ganglion injury accrues from VZV load; recovers slowly
  dxdt_GANG_INJ = KGANG_IN*VZV_LOAD*(1.0 - GANG_INJ) - KGANG_REC*GANG_INJ;

  // IENF: loss proportional to active ganglion injury and ongoing CAP defunctionalization;
  //       slow recovery
  dxdt_IENF = - KIENF_LOSS*100.0*GANG_INJ
              - 0.4*Edrug_CAP*IENF*0.01
              + KIENF_REC*(100.0 - IENF);

  // NGF response to injury (kept positive)
  dxdt_NGF_SKIN = NGF_IN*GANG_INJ - NGF_OUT*NGF_SKIN;

  // NaV peripheral ectopic activity — driven by NGF and IENF loss
  double iennorm = (100.0 - IENF)/100.0;
  dxdt_NAV_ACT = KNAV_IN*(NGF_SKIN + iennorm) - KNAV_OUT*NAV_ACT;

  // Microglial activation: amplified by CSEN, dampened by IL10/steroid surrogate (not modeled drug)
  dxdt_MICROG = KMG_IN*CSEN*(1.0 - MICROG) - KMG_OUT*MICROG;

  // KCC2 loss from microglial BDNF
  dxdt_KCC2 = -KKCC2_LOSS*MICROG*KCC2 + KKCC2_REC*(KCC2_BASE - KCC2);

  // NMDA tone — driven by glutamatergic transmission, gated by (1/KCC2)
  dxdt_NMDA_TONE = KNMDA_IN*CSEN/KCC2 - KNMDA_OUT*NMDA_TONE;

  // Central sensitization — input from NaV; dampened by α2δ ligands (GBP/PGB), descending NE (AMI, DLX, TRA),
  // microglial inhibitors, and topical lidocaine (peripheral input reduction).
  double central_drug = (1.0 - Edrug_PGB)*(1.0 - Edrug_GBP)*(1.0 - Edrug_AMI)*
                        (1.0 - Edrug_DLX)*(1.0 - Edrug_TRA);
  double periph_drug  = (1.0 - Edrug_LIDO)*(1.0 - Edrug_CAP);
  dxdt_CSEN = KCS_IN*NAV_ACT*periph_drug*central_drug*(1.0 - CSEN) - KCS_OUT*CSEN;

  // CMI immune surveillance — RZV boosts it, age erodes it
  double rzv_boost = 1.0 * RZV_VAC/(0.5 + RZV_VAC);
  dxdt_CMI = KCMI_REC*(CMI_BASE - CMI) + rzv_boost*0.05 - KCMI_LOSS*CMI*((AGE-50)>0?(AGE-50)/50.0:0);

  // ---- Pain endpoints ----
  // raw pain build-up
  double pain_raw =  W_NRS_NAV*NAV_ACT*periph_drug
                   + W_NRS_CSEN*CSEN*central_drug
                   + W_NRS_NMDA*NMDA_TONE*(1.0 - Edrug_CAP*0.3);
  dxdt_PAIN = (pain_raw - PAIN)*B_PAIN_REC*10.0;

  double allo_raw = W_ALLO*NAV_ACT*(100.0 - IENF)/100.0*(1.0 - Edrug_LIDO)*(1.0 - Edrug_CAP);
  dxdt_ALLO = (allo_raw - ALLO)*B_PAIN_REC*10.0;

  dxdt_SLEEP = (W_SLEEP*PAIN - SLEEP)*B_PAIN_REC*10.0;
  dxdt_MOOD  = (W_MOOD *PAIN - MOOD )*B_PAIN_REC*10.0;

  // Sedation/AE — driven by total CNS-active analgesic exposure
  double cns_exp = Cp_PGB/PGB_EC50 + Cp_GBP/GBP_EC50 +
                   (Cp_AMI+Cp_NOR)/AMI_EC50 + Cp_TRA/TRA_EC50;
  dxdt_AE_SED = KAE_SED*cns_exp - KAE_DECAY*AE_SED;

$TABLE
  capture nrs_pain      = (PAIN  > 10 ? 10 : PAIN);
  capture allo_score    = (ALLO  > 10 ? 10 : ALLO);
  capture sleep_score   = (SLEEP > 10 ? 10 : SLEEP);
  capture mood_score    = (MOOD  > 10 ? 10 : MOOD);
  capture nav_drive     = NAV_ACT;
  capture central_sens  = CSEN;
  capture vzv_log       = log10(VZV_LOAD + 1e-6);
  capture cmi_score     = CMI;
  capture ae_sed_score  = AE_SED;
'

mod <- mcode("phn_qsp", code)

# ============================================================================
# Scenario library (10 scenarios)
# ============================================================================
# Time grid: simulation typically 0–180 d (zoster) or 0–365 d (PHN management)

library(dplyr)

# --- 1. Untreated / placebo PHN ---
e_placebo  <- ev(amt = 0)

# --- 2. Acute zoster: valaciclovir 1 g q8h × 7 d ---
e_val <- ev(amt = 1000, ii = 8, addl = 20, cmt = "VAL_CENT")

# --- 3. RZV (Shingrix) prophylaxis: 2 doses 0 and 60 d ---
e_rzv <- ev(amt = 50, cmt = "RZV_VAC", time = 0) +
         ev(amt = 50, cmt = "RZV_VAC", time = 60*24)

# --- 4. Gabapentin titration: 300→600→900 tid (3600/d) ---
e_gbp <- ev(amt = 300, ii = 8, addl = 2, cmt = "GBP_GUT", time = 0) +
         ev(amt = 600, ii = 8, addl = 2, cmt = "GBP_GUT", time = 24) +
         ev(amt = 900, ii = 8, addl = 30*3 - 1, cmt = "GBP_GUT", time = 48)

# --- 5. Pregabalin 75→150 bid (300/d) ---
e_pgb <- ev(amt = 75,  ii = 12, addl = 5, cmt = "PGB_CENT", time = 0) +
         ev(amt = 150, ii = 12, addl = 80, cmt = "PGB_CENT", time = 72)

# --- 6. Amitriptyline 25 mg HS, titrate to 75 mg ---
e_ami <- ev(amt = 25, ii = 24, addl = 6,  cmt = "AMI_CENT") +
         ev(amt = 50, ii = 24, addl = 6,  cmt = "AMI_CENT", time = 7*24) +
         ev(amt = 75, ii = 24, addl = 60, cmt = "AMI_CENT", time = 14*24)

# --- 7. Duloxetine 30→60 mg qd ---
e_dlx <- ev(amt = 30, ii = 24, addl = 6,  cmt = "DLX_CENT") +
         ev(amt = 60, ii = 24, addl = 60, cmt = "DLX_CENT", time = 7*24)

# --- 8. Lidocaine 5% patch — apply daily 12 h ---
# Approximated as bolus into LIDO_SKIN at 0, 24, 48 ... (skin pool)
e_lido <- ev(amt = 70, ii = 24, addl = 60, cmt = "LIDO_SKIN")

# --- 9. Capsaicin 8% patch — single 60-min, modeled as bolus q90 d ---
e_cap  <- ev(amt = 100, ii = 90*24, addl = 3, cmt = "CAP_SKIN")

# --- 10. Combination (acute zoster + RZV-prophy + PGB + LIDO patch + AMI) ---
e_combo <- e_val + e_rzv + e_pgb + e_lido + e_ami

# ============================================================================
# Helper runner
# ============================================================================
run_scenario <- function(model, events, end_d = 120, label = "scenario") {
  out <- mrgsim_e(model, events, end = end_d*24, delta = 6) %>% as.data.frame()
  out$scenario <- label
  out$day      <- out$time/24
  out
}

if (interactive()) {
  scens <- list(
    placebo  = e_placebo,
    val      = e_val,
    rzv      = e_rzv,
    gbp      = e_gbp,
    pgb      = e_pgb,
    ami      = e_ami,
    dlx      = e_dlx,
    lido     = e_lido,
    cap      = e_cap,
    combo    = e_combo
  )
  all_out <- do.call(rbind, lapply(names(scens), function(nm)
    run_scenario(mod, scens[[nm]], 180, nm)))
  message("Scenarios simulated:\n", paste(names(scens), collapse=", "))
}
