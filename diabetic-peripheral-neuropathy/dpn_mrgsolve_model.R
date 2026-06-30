## =====================================================================
## Diabetic Peripheral Neuropathy (DPN) — mrgsolve QSP Model
## ---------------------------------------------------------------------
## Author: Claude Code Routine — QSP Disease Model Library
## File:   dpn_mrgsolve_model.R
##
## Scope (26 ODE compartments):
##   • Drug PK   — pregabalin, duloxetine, alpha-lipoic acid (ALA),
##                  epalrestat (aldose-reductase inhibitor), capsaicin
##                  8% patch effect compartment, lidocaine 5% patch.
##   • Disease PD — glycaemic burden (HbA1c, intracellular flux),
##                  polyol-pathway sorbitol, AGE accumulation, ROS,
##                  endoneurial blood flow, NGF, IENFD, NCV,
##                  C-fibre integrity, peripheral & central
##                  sensitisation, pain NRS (0-10), MNSI, BPI,
##                  sleep, QoL, foot-ulcer hazard.
##   • 7 treatment scenarios prepared in companion examples (below).
##
## Calibration anchors:
##   - DCCT/EDIC: tight glucose control → 60% reduction in
##                neuropathy incidence at 5y (NEJM 1993, Ann Intern
##                Med 2010).
##   - SENZA-PDN HF10 SCS: 76% pain responders at 3m vs CMM 5%
##                (Petersen JAMA Neurol 2021).
##   - Pregabalin 300-600 mg/d: mean NRS Δ −1.3 vs placebo
##                (Freeman Diabetes Care 2008).
##   - Duloxetine 60-120 mg/d: NRS Δ −1.4 (Goldstein Pain 2005).
##   - α-Lipoic acid 600 mg IV ×3w: TSS Δ −2.7 (ALADIN, Ziegler
##                Diabetologia 1995); oral NATHAN trial neutral at
##                4y but improved subgroups.
##   - Capsaicin 8%: NRS Δ −1.0 sustained 12 wk (STEP).
##   - Epalrestat: stabilises NCV over 3 y in mild DPN
##                (Hotta Diabetes Care 2006).
##
## Usage:
##   library(mrgsolve)
##   mod <- mread("dpn_mrgsolve_model.R")
##   ev_preg  <- ev(amt=300, ii=24, addl=364, cmt="GUT_PG", time=0)  # mg
##   out      <- mrgsim(mod, events=ev_preg, end=365, delta=1)
##   plot(out, NRS~time)
##
## Notes:
##   • Time is in DAYS. PK absorption uses fast 1st-order kinetics
##     averaged to daily; for hour-resolution PK switch delta=1/24.
##   • The pain score is bounded [0,10] via a logistic transform.
##   • Foot-ulcer hazard is a Cox-like cumulative incidence.
## =====================================================================

[ PROB ]
Diabetic Peripheral Neuropathy QSP Model — v1.0 (2026-06-30)
Domains: glycaemic drive, polyol-AGE-ROS, microvascular, nerve fibre
loss, pain sensitisation, clinical outcomes; six drug PK + non-pharm.

[ PLUGIN ] Rcpp

[ PARAM ] @annotated
// === Patient phenotype ===
HbA1c0      :  8.5 : Baseline HbA1c (%)
HbA1c_tgt   :  7.0 : Target HbA1c with treatment (%)
k_HbA1c     :  0.030 : HbA1c relaxation rate (1/day)
DM_DURATION : 10   : Diabetes duration at t0 (years)
AGE_PT      : 60   : Patient age (years)
BMI         : 30   : BMI (kg/m^2)
HTN_FLAG    :  1   : Hypertension (0/1)
SMOKE_FLAG  :  0   : Current smoker (0/1)
EGFR        : 75   : eGFR mL/min/1.73m^2 (affects pregabalin CL)

// === Glycaemia → intracellular flux ===
KFLUX       : 0.18 : Intracellular glucose flux scaling (per %HbA1c)
GLUC_VAR    : 1.0  : Glycaemic variability multiplier (1=normal)

// === Polyol pathway ===
kAR         : 0.40 : Aldose reductase rate (Sorbitol prod/day)
kSDH        : 0.30 : Sorbitol dehydrogenase (Sorbitol→Fruc)
kSORB_OUT   : 0.05 : Osmotic clearance
Epalrestat_IC50 : 5  : Epalrestat IC50 (mg/L) on AR

// === AGE / RAGE ===
kAGE_form   : 0.005 : AGE formation rate (per (Glu·Fruc))
kAGE_clr    : 0.0008 : AGE clearance (1/day)
kRAGE_sig   : 0.5   : RAGE→NFkB transduction
Amino_IC50  : 50    : Aminoguanidine IC50 on AGE formation

// === ROS / mitochondria ===
kROS_prod   : 1.0  : ROS production (per polyol + AGE)
kROS_clr    : 0.8  : ROS clearance (GSH-dependent)
GSH0        : 1.0  : Baseline GSH
kGSH_loss   : 0.05 : GSH loss per ROS
ALA_kbase   : 0.25 : ALA antioxidant rate constant
ALA_EC50    : 0.5  : ALA EC50 (mg/L plasma)
Nrf2_max    : 1.5  : Maximum Nrf2 induction

// === Microvascular ===
kVASA_dmg   : 0.02 : Vasa nervorum damage per ROS·AGE
kVASA_heal  : 0.01 : Vasa nervorum repair
ACEi_protect: 0.30 : ACEi/ARB protection fraction

// === Inflammation ===
kNFkB_act   : 0.5  : NFkB activation (ROS+RAGE driven)
kNFkB_dec   : 0.4  : NFkB decay
kTNF_prod   : 0.6  : TNFa production
kTNF_clr    : 0.7  : TNFa clearance

// === Nerve compartments ===
kNGF_prod   : 1.0  : NGF baseline production
kNGF_dec    : 0.05 : NGF degradation
kNGF_loss   : 0.4  : NGF loss per ROS/hypoxia
IENFD0      : 12   : Baseline IENFD (fibres/mm)
IENFD_min   : 1    : Minimum viable IENFD
kIENFD_loss : 0.0008 : Daily IENFD loss per damage unit
kIENFD_reg  : 0.002 : Regrowth rate (NGF dependent)
NCV0        : 50   : Baseline NCV (m/s, sural)
NCV_floor   : 25   : Floor for severe DPN
kNCV_loss   : 0.0010 : NCV decline rate
kNCV_rec    : 0.0005 : NCV recovery rate
kSchwann    : 0.001 : Schwann cell injury rate

// === Pain pathways ===
kSensP_on   : 0.3  : Peripheral sensitisation gain
kSensP_off  : 0.2  : Resolution
kSensC_on   : 0.25 : Central sensitisation gain
kSensC_off  : 0.18 : Central resolution
DESC_TONE0  : 1.0  : Baseline descending modulation
kDESC_dec   : 0.10 : Descending tone decline with disease
Capsaicin_def : 0.6 : Cfiber defunctionalisation amplitude per patch
Capsaicin_tau : 90  : Time-constant of effect (days)
Lidocaine_PMax : 1.0 : Lidocaine peripheral block efficacy (relative)
Lidocaine_tau  : 12  : Lidocaine residence (h equivalent in days)

// === Pregabalin PK (oral, 1-cmt) ===
PG_KA       : 12   : Absorption rate (1/day)
PG_VD       : 35   : Volume (L)  (~0.5 L/kg)
PG_CL_base  : 80   : Baseline CL (L/day) at eGFR 90  (~3.3 L/h)
PG_CLrenal_frac : 0.95 : Renal elimination fraction
PG_EC50_pain    : 4 : Plasma EC50 (mg/L) for pain effect
PG_Emax_pain    : 0.45 : Max fractional pain reduction
PG_EC50_sleep   : 3 : EC50 sleep
PG_Emax_sleep   : 0.30 : Max sleep effect

// === Duloxetine PK (oral, 1-cmt simplified) ===
DLX_KA      : 6    : Absorption (1/day)
DLX_VD      : 1640 : Apparent V/F (L) (~23 L/kg)
DLX_CL      : 2160 : CL/F (L/day) ~90 L/h
DLX_EC50    : 0.08 : EC50 (mg/L) for pain & descending tone
DLX_Emax    : 0.35 : Max pain reduction
DLX_CYP2D6  : 1.0  : Phenotype multiplier (PM:0.4 / EM:1.0 / UM:2.0)

// === Alpha-lipoic acid PK (oral) ===
ALA_KA      : 24   : Absorption (1/day) (rapid)
ALA_VD      : 350  : Apparent V (L) (~5 L/kg)
ALA_CL      : 3000 : CL (L/day, hepatic)
ALA_F_oral  : 0.30 : Oral bioavailability
ALA_F_IV    : 1.00 : IV bioavailability

// === Epalrestat PK ===
EP_KA       : 24   : Absorption (1/day)
EP_VD       : 12   : V (L)
EP_CL       : 65   : CL (L/day)

// === Capsaicin 8% patch — effect compartment ===
CAP_DOSE_FRAC : 0.05 : Fraction systemic (negligible);
                       drive via patch_on event
CAP_ON      : 0    : Switch (=1 day of patch application)

// === Hazard / outcomes ===
DFU_BASE_HAZ : 0.00015 : Baseline daily foot-ulcer hazard
DFU_BETA_IENFD : 0.10  : Hazard increase per 1 fibre/mm loss
DFU_BETA_PAIN  : 0.005 : Hazard increase per NRS unit
QoL_MAX     : 100  : Maximum QoL (Norfolk-QoL inv. scale)

[ CMT ] @annotated
// PK gut + central compartments per drug
GUT_PG     : Pregabalin gut amount (mg)
PLAS_PG    : Pregabalin plasma amount (mg)
GUT_DLX    : Duloxetine gut amount (mg)
PLAS_DLX   : Duloxetine plasma amount (mg)
GUT_ALA    : ALA gut amount (mg)
PLAS_ALA   : ALA plasma amount (mg)
GUT_EP     : Epalrestat gut amount (mg)
PLAS_EP    : Epalrestat plasma amount (mg)
CAP_EFF    : Capsaicin effect compartment (0-1)
LID_EFF    : Lidocaine effect compartment (0-1)
// Disease state compartments
HbA1c      : Glycated haemoglobin (%)
SORB       : Polyol/sorbitol burden (a.u.)
AGE        : AGE burden (a.u.)
ROS        : ROS level (a.u.)
GSH        : Glutathione (a.u.)
NFKB       : NFkB activity (a.u.)
TNFa       : TNFa (a.u.)
VASA       : Vasa nervorum integrity (1=intact, 0=lost)
NGF        : Nerve growth factor (a.u., 1=baseline)
IENFD      : Intra-epidermal nerve fibre density (fibres/mm)
NCV        : Sural nerve conduction velocity (m/s)
PERIPH_S   : Peripheral sensitisation (a.u.)
CENTRAL_S  : Central sensitisation (a.u.)
DESC_TONE  : Descending inhibition tone (a.u., 1=baseline)
PAIN_DRIVE : Latent pain drive (a.u.)
DFU_HAZ_C  : Cumulative foot-ulcer hazard

[ MAIN ]
// Renal adjustment for pregabalin
double PG_CL = PG_CL_base * pow(EGFR/90.0, PG_CLrenal_frac);
// CYP2D6-modulated duloxetine CL
double DLX_CL_eff = DLX_CL / DLX_CYP2D6;

[ ODE ]
// ===== PK: Pregabalin =====
dxdt_GUT_PG   = -PG_KA * GUT_PG;
dxdt_PLAS_PG  =  PG_KA * GUT_PG  - (PG_CL / PG_VD) * PLAS_PG;
double Cp_PG  = PLAS_PG / PG_VD;                 // mg/L

// ===== PK: Duloxetine =====
dxdt_GUT_DLX  = -DLX_KA * GUT_DLX;
dxdt_PLAS_DLX =  DLX_KA * GUT_DLX - (DLX_CL_eff / DLX_VD) * PLAS_DLX;
double Cp_DLX = PLAS_DLX / DLX_VD;               // mg/L

// ===== PK: α-Lipoic acid =====
dxdt_GUT_ALA  = -ALA_KA * GUT_ALA;
dxdt_PLAS_ALA =  ALA_KA * GUT_ALA - (ALA_CL / ALA_VD) * PLAS_ALA;
double Cp_ALA = PLAS_ALA / ALA_VD;

// ===== PK: Epalrestat =====
dxdt_GUT_EP   = -EP_KA * GUT_EP;
dxdt_PLAS_EP  =  EP_KA * GUT_EP - (EP_CL / EP_VD) * PLAS_EP;
double Cp_EP  = PLAS_EP / EP_VD;

// ===== Topical effect compartments =====
dxdt_CAP_EFF  = -CAP_EFF / Capsaicin_tau;       // wash-out only;
                                                // CAP_EFF set by event
dxdt_LID_EFF  = -LID_EFF / Lidocaine_tau;

// ===== Glycaemia =====
double targ = HbA1c_tgt;        // treatment-targeted HbA1c
dxdt_HbA1c  = -k_HbA1c * (HbA1c - targ);

// ===== Polyol / sorbitol =====
double AR_inh = 1.0 / (1.0 + Cp_EP / Epalrestat_IC50);
double bena   = 1.0; // benfotiamine flag could be added here
dxdt_SORB     =  kAR * (HbA1c/HbA1c0) * AR_inh
               - kSDH * SORB
               - kSORB_OUT * SORB;

// ===== AGE / RAGE =====
double amino_inh = 1.0; // not in default dosing
dxdt_AGE      =  kAGE_form * (HbA1c/HbA1c0) * SORB * amino_inh
               - kAGE_clr * AGE;

// ===== ROS / GSH =====
double ALA_eff = (Cp_ALA / (ALA_EC50 + Cp_ALA));
double GSH_eff = GSH;
dxdt_ROS      =  kROS_prod * (SORB + AGE) * GLUC_VAR
               - kROS_clr  * ROS * GSH_eff * (1.0 + ALA_kbase*ALA_eff);
dxdt_GSH      = -kGSH_loss * ROS * GSH
               + 0.05 * (GSH0 - GSH)
               + 0.10 * ALA_eff;

// ===== Inflammation =====
double rage_signal = kRAGE_sig * AGE;
dxdt_NFKB     =  kNFkB_act * (ROS + rage_signal) - kNFkB_dec * NFKB;
dxdt_TNFa     =  kTNF_prod * NFKB - kTNF_clr * TNFa;

// ===== Microvascular: vasa nervorum =====
double vasa_dmg = kVASA_dmg * ROS * AGE * (1.0 - ACEi_protect*0.0);
dxdt_VASA     = -vasa_dmg * VASA + kVASA_heal * (1.0 - VASA);

// ===== NGF =====
double hypoxia = (1.0 - VASA);
dxdt_NGF      =  kNGF_prod - kNGF_dec * NGF - kNGF_loss * ROS * hypoxia;

// ===== IENFD =====
double dmg_unit = ROS + AGE + 2.0*TNFa + 1.5*hypoxia;
dxdt_IENFD    = -kIENFD_loss * dmg_unit * (IENFD - IENFD_min)
               + kIENFD_reg  * NGF * (IENFD0 - IENFD);

// ===== NCV =====
dxdt_NCV      = -kNCV_loss * (AGE + 0.5*ROS + hypoxia) * (NCV - NCV_floor)
               + kNCV_rec  * NGF * (NCV0 - NCV);

// ===== Sensitisation =====
double cap_block = Capsaicin_def * CAP_EFF;
double lid_block = Lidocaine_PMax * LID_EFF;
double drive_pos = (ROS + TNFa + (1.0 - VASA)) * (1.0 - cap_block - 0.5*lid_block);
dxdt_PERIPH_S =  kSensP_on  * drive_pos - kSensP_off * PERIPH_S;
dxdt_CENTRAL_S=  kSensC_on  * (PERIPH_S + 0.5*TNFa) - kSensC_off * CENTRAL_S;

// ===== Descending modulation =====
double DLX_eff = Cp_DLX / (DLX_EC50 + Cp_DLX);
dxdt_DESC_TONE = -kDESC_dec * (DESC_TONE - DESC_TONE0*(1.0 - 0.3*hypoxia))
                + 0.4 * DLX_Emax * DLX_eff;

// ===== Pain drive (latent → NRS) =====
double PG_eff   = PG_Emax_pain * (Cp_PG / (PG_EC50_pain + Cp_PG));
double pain_in  = PERIPH_S + 1.2*CENTRAL_S - 0.8*DESC_TONE - PG_eff;
dxdt_PAIN_DRIVE = 0.2 * (pain_in - PAIN_DRIVE);

// ===== Foot-ulcer hazard =====
double IENFD_loss = (IENFD0 - IENFD);
double NRS_now    = 10.0 * (1.0 / (1.0 + exp(-(PAIN_DRIVE - 1.0))));
double haz        = DFU_BASE_HAZ
                  * exp(DFU_BETA_IENFD*IENFD_loss + DFU_BETA_PAIN*NRS_now);
dxdt_DFU_HAZ_C    = haz;

[ TABLE ]
// Bounded pain NRS (0..10) via logistic on PAIN_DRIVE
double NRS = 10.0 * (1.0 / (1.0 + exp(-(PAIN_DRIVE - 1.0))));
// Symptom & clinical composites
double MNSI    = 13.0
               - 0.9 * (IENFD / IENFD0) * 4.0
               - 0.7 * (NCV   / NCV0)   * 4.0
               - 0.4 * VASA            * 2.0
               + 0.5 * NRS;
if (MNSI < 0)  MNSI = 0;
if (MNSI > 13) MNSI = 13;

double TCNS    = 2.0 * (1.0 - IENFD/IENFD0) * 10.0
               + 2.0 * (1.0 - NCV/NCV0)     * 9.0;
if (TCNS > 19) TCNS = 19;

double BPI_INT = 0.7 * NRS + 0.3 * (1.0 - VASA) * 10.0;
double SLEEP_INTERF = 0.5 * NRS - 0.4 * PG_Emax_sleep * (Cp_PG/(PG_EC50_sleep+Cp_PG));
if (SLEEP_INTERF < 0) SLEEP_INTERF = 0;

double Norfolk_QoL = QoL_MAX
                   - 4.0 * NRS
                   - 3.0 * (IENFD0 - IENFD)
                   - 2.0 * (NCV0 - NCV)
                   - 5.0 * SLEEP_INTERF;
if (Norfolk_QoL < 0) Norfolk_QoL = 0;

// Cumulative foot-ulcer incidence (Cox-like)
double DFU_INC = 1.0 - exp(-DFU_HAZ_C);

[ CAPTURE ] @annotated
Cp_PG    : Pregabalin plasma (mg/L)
Cp_DLX   : Duloxetine plasma (mg/L)
Cp_ALA   : ALA plasma (mg/L)
Cp_EP    : Epalrestat plasma (mg/L)
NRS      : Pain numerical rating scale (0-10)
MNSI     : Michigan Neuropathy Symptom Score
TCNS     : Toronto Clinical Neuropathy Score
BPI_INT  : Brief Pain Inventory — interference
SLEEP_INTERF : Sleep interference (0-10)
Norfolk_QoL  : Norfolk QoL (0-100)
DFU_INC  : Cumulative foot-ulcer incidence
hypoxia  : Endoneurial hypoxia index
IENFD_loss : IENFD loss vs baseline
NRS_now  : Real-time NRS used in hazard

/*
=========================================================================
EXAMPLE SCENARIOS (run in R after `mod <- mread(...)`):

# 1) Untreated control — 1 year
ev0 <- ev(amt=0, time=0)
s0  <- mrgsim(mod, events=ev0, end=365, delta=1)

# 2) Pregabalin 300 mg PO BID
ev1 <- ev(amt=150, ii=12, addl=730, cmt="GUT_PG")
s1  <- mrgsim(mod, events=ev1, end=365, delta=1)

# 3) Duloxetine 60 mg PO QD
ev2 <- ev(amt=60, ii=24, addl=364, cmt="GUT_DLX")
s2  <- mrgsim(mod, events=ev2, end=365, delta=1)

# 4) α-Lipoic acid 600 mg IV daily ×3 wk then oral 600 mg/d
ev3a <- ev(amt=600, ii=24, addl=20, cmt="PLAS_ALA")    # IV bolus
ev3b <- ev(amt=600, ii=24, addl=340, cmt="GUT_ALA", time=21)
s3   <- mrgsim(mod, events=c(ev3a, ev3b), end=365, delta=1)

# 5) Epalrestat 150 mg PO QD (Japan)
ev4 <- ev(amt=150, ii=24, addl=364, cmt="GUT_EP")
s4  <- mrgsim(mod, events=ev4, end=365, delta=1)

# 6) Capsaicin 8% patch q90d — apply by pulsing CAP_EFF=1.0
ev5 <- ev(amt=1, time=c(0,90,180,270), cmt="CAP_EFF")
s5  <- mrgsim(mod, events=ev5, end=365, delta=1)

# 7) Combination: pregabalin + duloxetine + ALA + intensive glucose
mod2 <- update(mod, param=list(HbA1c_tgt=6.5))
ev6  <- c(ev(amt=150, ii=12, addl=730, cmt="GUT_PG"),
          ev(amt=60,  ii=24, addl=364, cmt="GUT_DLX"),
          ev(amt=600, ii=24, addl=364, cmt="GUT_ALA"))
s6   <- mrgsim(mod2, events=ev6, end=365, delta=1)
=========================================================================
*/
