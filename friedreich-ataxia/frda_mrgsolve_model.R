# =============================================================================
# Friedreich Ataxia (FRDA) — QSP / mrgsolve model
# -----------------------------------------------------------------------------
# Captures:
#   * GAA-repeat → frataxin (FXN) protein dynamics
#   * Mitochondrial Fe-S cluster pool, mitochondrial labile-iron, ETC capacity,
#     ATP, ROS, NRF2-ARE antioxidant pool
#   * Tissue compartments: dorsal root ganglion + cerebellar (CNS),
#     cardiomyocyte (LV mass), pancreatic β-cell (insulin → glucose)
#   * Clinical scales: mFARS, T25FW (proxy ambulation)
#   * Drugs:
#       - Omaveloxolone (Skyclarys 150 mg QD)   : NRF2 stabilizer
#       - Idebenone PO                          : antioxidant / ETC bypass
#       - Deferiprone PO                        : mitochondrial Fe chelator
#       - Nomlabofusp SC (CTI-1601, TAT-FXN)    : protein replacement
#       - AAVrh10-FXN gene therapy (single IV)  : durable FXN restoration
#       - ACEi (cardiac)                        : LV mass reduction
# -----------------------------------------------------------------------------
# Author : QSP library
# Use    : mrgsolve::mread("frda_mrgsolve_model.R") |> mrgsim(events = ev(...))
# =============================================================================

library(mrgsolve)

code <- '
$PROB
# Friedreich Ataxia QSP model (15+ compartments, multi-drug)

$PARAM @annotated
// ----- Genotype / disease severity -----
GAA1     :   650 : Short allele GAA repeats
GAA2     :   850 : Long allele GAA repeats
AAO      :    12 : Age at onset (years)
SEX      :     1 : Sex (1=F,0=M)

// ----- Frataxin synthesis / turnover -----
FXN_max  : 100   : Maximum FXN expression (% of WT)
kFXN_in  :  0.06 : FXN synthesis rate (1/day)
kFXN_out :  0.06 : FXN degradation rate (1/day)
GAA_IC50 :  500  : GAA repeat IC50 for FXN suppression

// ----- Fe-S biogenesis / mito iron -----
kFeS_in  :  1.0  : Fe-S synthesis rate constant (1/day)
kFeS_out :  0.10 : Fe-S turnover (1/day)
KmFXN_FeS:  20   : FXN level for half-max Fe-S synthesis
kFe_in   :  0.05 : Mito Fe import (1/day, baseline)
kFe_out  :  0.04 : Fe export (1/day)
Fe_basal :  1.0  : Baseline mito labile iron (a.u.)

// ----- ETC / ATP -----
ETC_max  :  100  : Max ETC capacity (% normal)
kETC_FeS :  0.10 : Fe-S to ETC coupling (1/day)
kETC_loss:  0.05 : ETC turnover (1/day)
kATP     :  1.5  : ATP regen rate (1/day per unit ETC)
kATPuse  :  1.5  : ATP utilization (1/day)

// ----- ROS / antioxidant balance -----
kROSgen_b:  0.20 : Baseline ROS generation (a.u./day)
kROSgen_F:  0.50 : Iron-driven ROS scaling (a.u./day per Fe unit above 1)
kROS_AOX :  1.0  : ROS quenching by AOX (1/day per unit AOX)
AOX_b    :  1.0  : Baseline AOX (a.u.)
kAOX_in  :  0.10 : AOX synthesis (1/day)
kAOX_out :  0.10 : AOX turnover (1/day)
EmaxNRF2 :  3.0  : Max NRF2 fold-induction of AOX

// ----- DRG / sensory neuron pool -----
DRG_b    :  100  : Baseline DRG neuron count (%)
kDRG_loss:  0.0015: DRG loss rate (1/day) at unit ROS above baseline
kDRG_min :  20   : Asymptotic floor (%)

// ----- Cerebellar neuronal function (dentate) -----
CB_b     :  100  : Baseline cerebellar function (%)
kCB_loss :  0.0012: Cerebellar loss/day at ROS=1
kCB_rec  :  0.0002: Recovery rate (PT/OT)

// ----- Cardiac LV mass (HCM-like) -----
LVMI_b   :  60   : Baseline LVMI g/m²
kLVMI_up :  0.015: LV mass growth (1/day at low ATP)
LVMI_max :  200  : Maximum LVMI g/m²
ATP_LV50 :  0.6  : ATP level at half-max LV hypertrophy
kLVMI_dn :  0.005: LV mass regression (1/day, with ACEi)

// ----- β-cell / glucose -----
Bcell_b  :  100  : Baseline β-cell mass (%)
kBcell_loss: 0.0010: β-cell loss/day (ROS-driven)
kGlu_in  :  6.0  : Endogenous glucose input (mmol/L/day)
kGlu_out :  0.5  : Glucose clearance (1/day per insulin unit)
Insulin_max: 100 : Max insulin secretion (μU/mL eq.)

// ----- mFARS clinical scale -----
mFARS_max :  93  : Maximum mFARS score
kmFARS_DRG: 0.50 : Slope: mFARS per %DRG loss
kmFARS_CB : 0.30 : Slope: mFARS per %CB loss

// ----- T25FW (timed 25-foot walk, sec) -----
T25FW_b   :  5.0 : Baseline T25FW seconds
kT25FW_dis:  0.10: Slope per unit mFARS above 30

// ----- Omaveloxolone PK (150 mg PO QD, fed) -----
OMAV_ka  : 0.45  : Absorption rate (1/h)
OMAV_CL  : 1.5   : Apparent clearance (L/h)
OMAV_V   : 250   : Apparent Vd (L)
OMAV_F   : 0.40  : Bioavailability (food state)
OMAV_EC50: 30    : NRF2 activation EC50 (ng/mL)
OMAV_Emax: 2.5   : Max NRF2 fold-induction
OMAV_AE  : 0.001 : ALT elevation per ng/mL/day

// ----- Idebenone PK (PO TID) -----
IDB_ka   : 1.2
IDB_CL   : 250
IDB_V    : 90
IDB_EC50 : 500   : ETC bypass EC50 (ng/mL)
IDB_Emax : 0.30  : Max +30% ETC capacity

// ----- Deferiprone PK (PO TID) -----
DFP_ka   : 1.5
DFP_CL   : 14
DFP_V    : 110
DFP_EC50 : 50    : Mito Fe chelation EC50 (μmol/L)
DFP_Emax : 0.60  : Max 60% reduction Fe_in

// ----- Nomlabofusp SC (CTI-1601, TAT-FXN) -----
NOM_ka   : 0.020 : SC absorption (1/h)
NOM_CL   : 0.30
NOM_V    : 5
NOM_FXN_delta: 25 : FXN restoration at Css (% WT)
NOM_EC50 : 100   : Css for half-max FXN restoration

// ----- AAV gene therapy (single IV dose) -----
AAV_dose_flag: 0 : Flag for AAV dose
AAV_t_on :  21   : Days to onset of expression
AAV_FXN_delta: 50: Steady FXN restoration (% WT)
AAV_decay: 0     : Decay (effectively durable)

// ----- ACEi -----
ACEi_eff : 0.0   : ACEi effect flag (0–1)

$CMT @annotated
// PK depots (oral / SC)
OMAV_GUT : Omaveloxolone gut depot (mg)
OMAV_CEN : Omaveloxolone central (mg)
IDB_GUT  : Idebenone gut (mg)
IDB_CEN  : Idebenone central (mg)
DFP_GUT  : Deferiprone gut (mg)
DFP_CEN  : Deferiprone central (mg)
NOM_SC   : Nomlabofusp SC depot (mg)
NOM_CEN  : Nomlabofusp central (mg)

// Disease / physiology
FXN      : Frataxin protein (% WT)
FeS      : Fe-S cluster pool (% normal)
Fe_mito  : Mitochondrial labile iron (a.u.)
ETC      : ETC capacity (% normal)
ATP      : Cellular ATP (a.u.)
ROS      : ROS pool (a.u.)
AOX      : Antioxidant capacity (a.u.)
DRG_pool : DRG neuron count (%)
CB_func  : Cerebellar function (%)
LVMI     : LV mass index (g/m²)
Bcell    : β-cell mass (%)
Glucose  : Plasma glucose (mmol/L)
AAV_FXN  : AAV-delivered FXN (% WT) (added to endogenous)

// Safety / PD trackers
ALT_OMAV : ALT elevation tracker (U/L)

$GLOBAL
#define OMAV_CP (OMAV_CEN/OMAV_V*1000.0)    // ng/mL
#define IDB_CP  (IDB_CEN/IDB_V*1000.0)
#define DFP_CP  (DFP_CEN/DFP_V*1.0)         // μmol/L proxy
#define NOM_CP  (NOM_CEN/NOM_V*1000.0)

$MAIN
// -- Initial conditions (steady-state of untreated FRDA, by GAA) --
double GAA_eff = (GAA1+GAA2)/2.0;
double FXN_ss  = FXN_max / (1.0 + pow(GAA_eff/GAA_IC50, 1.8));
FXN_0   = FXN_ss;
FeS_0   = 100.0 * FXN_ss / (FXN_ss + KmFXN_FeS);
Fe_mito_0 = Fe_basal * (1.0 + 0.8*(1.0 - FXN_ss/100.0));
ETC_0   = 100.0 * FeS_0/100.0;
ATP_0   = ETC_0/100.0;
ROS_0   = kROSgen_b + kROSgen_F*(Fe_mito_0 - 1.0);
AOX_0   = AOX_b;
DRG_pool_0 = DRG_b;
CB_func_0  = CB_b;
LVMI_0  = LVMI_b * (1.0 + 0.5*(1.0 - FXN_ss/100.0));
Bcell_0 = Bcell_b;
Glucose_0 = 5.5;
AAV_FXN_0 = 0;
ALT_OMAV_0 = 25;  // baseline U/L

$ODE
// ------- Drug PK -------
dxdt_OMAV_GUT = -OMAV_ka * OMAV_GUT;
dxdt_OMAV_CEN =  OMAV_ka * OMAV_GUT * OMAV_F - (OMAV_CL/OMAV_V) * OMAV_CEN;

dxdt_IDB_GUT  = -IDB_ka * IDB_GUT;
dxdt_IDB_CEN  =  IDB_ka * IDB_GUT - (IDB_CL/IDB_V) * IDB_CEN;

dxdt_DFP_GUT  = -DFP_ka * DFP_GUT;
dxdt_DFP_CEN  =  DFP_ka * DFP_GUT - (DFP_CL/DFP_V) * DFP_CEN;

dxdt_NOM_SC   = -NOM_ka * NOM_SC;
dxdt_NOM_CEN  =  NOM_ka * NOM_SC - (NOM_CL/NOM_V) * NOM_CEN;

// ------- Frataxin -------
double FXN_target = FXN_max / (1.0 + pow(GAA_eff/GAA_IC50, 1.8));
double nom_eff    = NOM_FXN_delta * NOM_CP / (NOM_CP + NOM_EC50);
double FXN_total  = FXN + AAV_FXN + nom_eff;        // effective FXN
dxdt_FXN    = kFXN_in * FXN_target - kFXN_out * FXN;

// AAV expression (simple onset, durable)
double AAV_drive = AAV_dose_flag * (TIME > AAV_t_on ? 1.0 : 0.0);
dxdt_AAV_FXN = (AAV_drive ? (AAV_FXN_delta - AAV_FXN)/14.0 : 0.0) - AAV_decay*AAV_FXN;

// ------- Fe-S -------
dxdt_FeS = kFeS_in * 100.0 * FXN_total/(FXN_total + KmFXN_FeS) - kFeS_out * FeS;

// ------- Mito Fe -------
double dfp_chel = DFP_Emax * DFP_CP/(DFP_CP + DFP_EC50);
dxdt_Fe_mito = kFe_in * (1.0 - dfp_chel) * (1.0 + 0.8*(1.0 - FXN_total/100.0))
              - kFe_out * (Fe_mito - 1.0);

// ------- ETC -------
double ide_bonus = IDB_Emax * IDB_CP/(IDB_CP + IDB_EC50);
dxdt_ETC = kETC_FeS * (FeS - ETC) + kETC_loss * (100.0*(1.0 + ide_bonus) - ETC) * 0.10;

// ------- ATP -------
dxdt_ATP = kATP * (ETC/100.0) - kATPuse * ATP;

// ------- NRF2 / AOX -------
double NRF2_drive  = 1.0 + (OMAV_Emax - 1.0) * OMAV_CP/(OMAV_CP + OMAV_EC50);
double AOX_target  = AOX_b * NRF2_drive;
dxdt_AOX = kAOX_in * AOX_target - kAOX_out * AOX;

// ------- ROS -------
double ROS_gen  = kROSgen_b + kROSgen_F * (Fe_mito - 1.0);
double ROS_quench = kROS_AOX * AOX;
dxdt_ROS = ROS_gen - ROS_quench * ROS;

// ------- DRG sensory pool -------
double drg_drive = kDRG_loss * (ROS - kROSgen_b);
if (drg_drive < 0) drg_drive = 0;
dxdt_DRG_pool = -drg_drive * (DRG_pool - kDRG_min);

// ------- Cerebellar function -------
double cb_drive = kCB_loss * (ROS - kROSgen_b);
if (cb_drive < 0) cb_drive = 0;
dxdt_CB_func = -cb_drive * CB_func + kCB_rec * (100.0 - CB_func);

// ------- LV mass -------
double hyp_drive = kLVMI_up * (1.0 - ATP/ATP_LV50);
if (hyp_drive < 0) hyp_drive = 0;
dxdt_LVMI = hyp_drive * (LVMI_max - LVMI) - kLVMI_dn * ACEi_eff * (LVMI - LVMI_b);

// ------- β-cell / glucose -------
double bcell_drive = kBcell_loss * (ROS - kROSgen_b);
if (bcell_drive < 0) bcell_drive = 0;
dxdt_Bcell = -bcell_drive * Bcell;
double Ins_sec = Insulin_max * (Bcell/100.0);
dxdt_Glucose = kGlu_in - kGlu_out * (Ins_sec/100.0) * Glucose;

// ------- Safety: ALT (omav class effect) -------
dxdt_ALT_OMAV = OMAV_AE * OMAV_CP - 0.05*(ALT_OMAV - 25.0);

$TABLE
double mFARS = kmFARS_DRG * (DRG_b - DRG_pool) + kmFARS_CB * (CB_b - CB_func);
if (mFARS < 0)         mFARS = 0;
if (mFARS > mFARS_max) mFARS = mFARS_max;

double T25FW = T25FW_b * (1.0 + kT25FW_dis * (mFARS > 30 ? (mFARS - 30) : 0)/30.0);
double HbA1c = 4.5 + 0.6 * (Glucose - 5.0);
if (HbA1c < 4.5) HbA1c = 4.5;

double LVH_flag = (LVMI > 95) ? 1.0 : 0.0;
double DM_flag  = (HbA1c >= 6.5) ? 1.0 : 0.0;

$CAPTURE @annotated
FXN     : Frataxin (%WT)
FeS     : Fe-S pool (%)
ETC     : ETC capacity (%)
ATP     : ATP (a.u.)
ROS     : ROS (a.u.)
AOX     : Antioxidant capacity
DRG_pool: DRG (%)
CB_func : Cerebellar (%)
LVMI    : LV mass index (g/m²)
Bcell   : β-cell (%)
Glucose : Glucose (mmol/L)
mFARS   : mFARS clinical score
T25FW   : Timed 25-foot walk (s)
HbA1c   : HbA1c (%)
LVH_flag: LV hypertrophy flag
DM_flag : Diabetes flag
ALT_OMAV: ALT (U/L)
OMAV_CP : Omav plasma (ng/mL)
'

mod <- mcode("frda_qsp", code)

# =============================================================================
# Scenario library
# Usage example:
#   library(mrgsolve); library(dplyr); library(ggplot2)
#   out <- mod %>% mrgsim(events = ev_scenarios()$omav_chronic, end = 365, delta = 1)
#   plot(out, mFARS + LVMI + Glucose + FXN ~ time)
# =============================================================================

ev_scenarios <- function() {
  list(
    # 1. Natural history (no drug)
    natural   = ev(amt = 0, cmt = "OMAV_GUT", time = 0),

    # 2. Omaveloxolone 150 mg QD chronic
    omav_chronic = ev(amt = 150, cmt = "OMAV_GUT", ii = 24, addl = 364),

    # 3. Idebenone 450 mg TID (legacy)
    idebenone = ev(amt = 450, cmt = "IDB_GUT", ii = 8, addl = 365*3-1),

    # 4. Deferiprone 25 mg/kg/day BID (Fe chelation)
    deferiprone = ev(amt = 800, cmt = "DFP_GUT", ii = 12, addl = 365*2-1),

    # 5. Nomlabofusp 50 mg SC QD (FXN replacement)
    nomlabofusp = ev(amt = 50, cmt = "NOM_SC", ii = 24, addl = 364),

    # 6. AAVrh10-FXN single IV dose (set AAV_dose_flag externally)
    aav_gene = ev(amt = 0, cmt = "OMAV_GUT"),

    # 7. Combination: omaveloxolone + ACEi + insulin (real-world)
    combo = c(
      ev(amt = 150, cmt = "OMAV_GUT", ii = 24, addl = 364)
    )
  )
}

# =============================================================================
# CALIBRATION NOTES
# -----------------------------------------------------------------------------
# 1. FXN-vs-GAA: Saveliev 2003 / Punga 2006: FXN ≈ 5–35% WT at GAA ~600–1000
# 2. Omaveloxolone MOXIe Part 2 (NCT02255435): mFARS Δ = -2.41 (placebo-corrected)
#    at 48w; 150 mg QD; ALT elevation in ~30%
# 3. Cardiac LV: 60–80% of FRDA patients have LV hypertrophy (LVMI > 95 g/m²)
# 4. Diabetes: ~30% of FRDA patients develop DM by adulthood
# 5. Deferiprone Boddaert 2007: 20 mg/kg/d for 6 mo → ↓cardiac iron T2*; mixed
#    neurologic results (some worsening of ataxia at higher doses, Pandolfo 2014)
# 6. Idebenone: equivocal in IONIA, declined approval; CoQ analog used historically
# 7. Nomlabofusp (CTI-1601): Phase 1 data shows dose-dependent FXN restoration
# 8. AAVrh10-FXN (LX2006): Preclinical durable expression; ongoing trials
# =============================================================================
