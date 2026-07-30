# =============================================================================
#  cca_mrgsolve_model.R
#  Cholangiocarcinoma (biliary tract cancer) — QSP model for mrgsolve
# =============================================================================
#
#  54 ODE compartments.  TIME IS DAYS FROM THE START OF SYSTEMIC THERAPY.
#
#  ---------------------------------------------------------------------------
#  WHAT THIS MODEL IS FOR
#  ---------------------------------------------------------------------------
#  Cholangiocarcinoma is usually written as a linear cascade — risk factor,
#  biliary mutation, tumour grows, chemotherapy shrinks it, patient lives
#  longer or does not.  That cascade cannot survive seven observations:
#
#   (A) Response and survival are only loosely coupled.  ABC-02 gemcitabine/
#       cisplatin: ORR ~26%, mPFS 8.0 mo, mOS 11.7 mo.  Three quarters of the
#       patients who benefit never meet a RECIST response.
#       [Valle 2010 PMID 20375404]
#   (B) The commonest terminal event is biliary, not bulk.  Cholangitis and
#       obstructive liver failure arrive long before tumour volume alone would
#       be lethal, especially in perihilar and distal disease.
#   (C) Biliary drainage has ZERO antitumour activity and is nonetheless the
#       intervention without which nothing else can be given at dose — both
#       gemcitabine and cisplatin are held or reduced on bilirubin.
#   (D) FGFR2-fusion tumours nearly always progress on pemigatinib through
#       POLYCLONAL kinase-domain mutations (V564F, N550K, E566A, K660M) that
#       are detectable in ctDNA before the imaging changes.
#       [Goyal 2017 PMID 28034880; Goyal 2019 PMID 31575540]
#   (E) TOPAZ-1 moved median OS from 11.5 to 12.8 months — barely — while
#       moving 24-month OS from 10.4% to 24.9%.  A model that tracks a mean
#       tumour size cannot separate tails.  [Oh 2022 PMID 36406756]
#   (F) Ivosidenib in IDH1-mutant iCCA: ORR ~2%, yet mPFS 2.7 vs 1.4 mo,
#       HR 0.37.  Benefit with essentially no shrinkage.
#       [Abou-Alfa 2020 PMID 32416072; Zhu 2021 PMID 34554208]
#   (G) Anatomic subtype reorders the failure mode.  FGFR2 fusions and IDH1
#       mutations are essentially confined to intrahepatic disease; hilar
#       disease fails by obstruction.
#
#  ---------------------------------------------------------------------------
#  THREE STRUCTURAL COMMITMENTS
#  ---------------------------------------------------------------------------
#  1. THE DELIVERED DOSE IS AN OUTPUT, NOT AN INPUT.
#     Every trial protocol writes gemcitabine 1000 mg/m2 on days 1 and 8, and
#     every trial then reports a relative dose intensity well below 100%.  In
#     this model the prescribed dose enters as an intention and is multiplied,
#     at the moment of administration, by four gates computed from the state:
#
#        RULE 1  bilirubin      GATE = 1 / (1 + (BILI/BILI50)^6),  BILI50 2.6
#        RULE 2  neutrophils    GATE = 1 / (1 + (1.0/ANC)^8)
#        RULE 3  creat. clear.  GATE = 1 / (1 + (45/CRCL)^8)   [cisplatin only]
#        RULE 4  performance    GATE = 1 / (1 + (PS/2.6)^8)
#
#     Because tumour at the hilum produces obstruction, obstruction produces
#     bilirubin, and bilirubin closes rule 1, the model contains a CLOSED
#     POSITIVE-FEEDBACK LOOP: growth -> obstruction -> held dose -> growth.
#     Biliary drainage is the only edge that cuts it, and drainage kills no
#     tumour cells anywhere in this file.  That is the point.
#
#     THE FALSIFIER IS ONE PARAMETER.  GATE_ON = 0 removes all four rules.
#     Then delaying drainage by 60 days must stop mattering.
#
#  2. RESISTANCE IS SELECTED, NOT INDUCED.
#     The FGFR2 kinase-domain-mutant clone T_R is seeded at t = 0 at
#     frequency MU_RES * ln(N_cells) — Goldie-Coldman — and nothing in this
#     file ever converts a sensitive cell into a mutant one.  There is no
#     "time to resistance" parameter.
#     THIS COMMITMENT PRODUCED A NEGATIVE RESULT AND IT IS REPORTED, NOT
#     HIDDEN: at a defensible seeding frequency (~2.5e-5) the clone cannot
#     reach clinical dominance inside a 7-month PFS — the arithmetic forbids
#     it, ln(0.2/2.5e-5)/0.026 per day is about 11 months.  So the model
#     assigns median PFS on an FGFR inhibitor to DRUG TOLERANCE (T_P) and
#     assigns the clone to the later events it can actually explain: loss of
#     an established response at around a year, ctDNA rising before imaging,
#     and the cross-resistance pattern between a reversible and a covalent
#     inhibitor.  Setting MU_RES = 0 therefore changes almost nothing in the
#     first year, and the model says so.
#
#  3. SURVIVAL IS TWO COMPETING HAZARDS ON DIFFERENT CLOCKS.
#        h_tumour  slow, burden-driven:  volume/liver + reserve loss + cachexia
#        h_biliary fast, recurrent:      cholangitis burden + ALBI
#     S(t) = exp(-INT(h_tumour + h_biliary)).  Population survival is the MEAN
#     of individual S(t), which is why the model can report a median and a
#     24-month tail from the same run — and why immunotherapy, which acts on
#     a MIXTURE fraction PI_IMMUNE of tumours rather than on everybody a
#     little, moves the tail without moving the median.
#
#  ---------------------------------------------------------------------------
#  WHAT WAS FITTED
#  ---------------------------------------------------------------------------
#  Eight parameters, marked [FIT] below:
#     K_GEM, K_CIS      cytotoxic potency
#     KSP0, ALPHA_P     drug-tolerant persister entry and residual killing
#     KFGFR_KILL        FGFR-inhibitor cytotoxic component
#     KIMM              immune kill scale
#     HT_T, HB_CH       the two hazard slopes
#  Everything else is a published PK parameter, a label-derived dose, a
#  physiological constant, or a structural choice fixed a priori (PSI = 20,
#  the Friberg transit count = 3, the Goldie-Coldman seeding rule).
#
#  ---------------------------------------------------------------------------
#  RELATIONSHIP TO cca_reference_impl.py
#  ---------------------------------------------------------------------------
#  The Python file in this directory is an exact twin of the equations below
#  (same 54 states, same names, same numbers) and it RUNS WITHOUT ANY
#  DEPENDENCIES.  Every number in README.md was produced by that file.  The
#  two differ in exactly two mechanical respects, both documented here:
#    * stent exchange.  Python re-stents when patency falls below 0.25 (a
#      state-triggered event).  mrgsolve has no state-triggered events, so the
#      R driver below schedules exchanges at fixed intervals instead.
#    * integrator.  Python uses fixed-step RK4 with substepping around the
#      infusions because gemcitabine has a 6-minute plasma half-life; mrgsolve
#      uses LSODA and needs none of that.
#
#  Render / run:
#     library(mrgsolve); mod <- mread("cca_mrgsolve_model.R")
#     source("cca_mrgsolve_model.R")   # then e.g. run_topaz1()
# =============================================================================

library(mrgsolve)
library(dplyr)

code <- '
$PROB
# Cholangiocarcinoma QSP model
# 54 ODEs. Delivered dose is an output; resistance is a pre-existing clone;
# survival is two competing hazards on different clocks.

$PARAM @annotated
// ---------------- patient / anatomy ----------------
BSA      : 1.72   : Body surface area (m2)
LIVVOL   : 1500   : Liver parenchymal volume (mL)
FHILAR   : 0.72   : Hilar-position fraction (pCCA/dCCA high, iCCA low)
IMMENG   : 0      : Immune-engaged tumour indicator (mixture member, 0/1)
FGFR2    : 0      : FGFR2 fusion present (0/1)
IDH1     : 0      : IDH1 R132 mutation present (0/1)
TT0      : 100    : Baseline tumour volume for RECIST reference (cm3)

// ---------------- gemcitabine PK ----------------
CL_GEM   : 3800   : Gemcitabine clearance (L/d)
V1_GEM   : 22     : Gemcitabine central volume (L)
Q_GEM    : 40     : Gemcitabine intercompartmental clearance (L/d)
V2_GEM   : 35     : Gemcitabine peripheral volume (L)
VMAX_TP  : 95     : Maximal dFdCTP formation rate (AU/d)
KM_TP    : 6      : dCK half-saturating gemcitabine concentration (mg/L)
KOUT_TP  : 2.08   : dFdCTP elimination (1/d)

// ---------------- cisplatin PK (free platinum) ----------------
CL_CIS   : 430    : Free platinum clearance (L/d)
V1_CIS   : 18     : Free platinum central volume (L)
Q_CIS    : 25     : Platinum intercompartmental clearance (L/d)
V2_CIS   : 45     : Platinum peripheral volume (L)
KFORM_PT : 9      : Pt-DNA adduct formation (AU per mg/L per d)
KREP_PT  : 0.85   : ERCC1 nucleotide excision repair (1/d)

// ---------------- durvalumab PK ----------------
CL_DUR   : 0.232  : Durvalumab clearance (L/d)
V1_DUR   : 5.6    : Durvalumab central volume (L)
Q_DUR    : 0.68   : Durvalumab intercompartmental clearance (L/d)
V2_DUR   : 3.9    : Durvalumab peripheral volume (L)
IC50_DUR : 0.35   : PD-L1 occupancy EC50 (mg/L)

// ---------------- FGFR inhibitor PK ----------------
KA_FGI   : 12     : Oral absorption rate (1/d)
CL_FGI   : 254    : Apparent clearance CL/F (L/d)
V_FGI    : 235    : Apparent volume V/F (L)
IC50_FGI : 0.008  : FGFR pathway IC50 in plasma units (mg/L)
KCOV     : 2.2    : Covalent inactivation rate, futibatinib (1/d)
KCOV_OFF : 0.10   : FGFR2 receptor resynthesis (1/d)
RHO      : 0.05   : Residual activity vs gatekeeper (0.05 reversible 0.45 covalent)
COVAL    : 0      : Inhibitor is covalent (0/1)

// ---------------- ivosidenib ----------------
CL_IVO   : 6.5    : Ivosidenib apparent clearance (L/d)
V_IVO    : 180    : Ivosidenib apparent volume (L)
IC50_IVO : 0.9    : 2-HG suppression EC50 (mg/L)
EMAX_IDH : 0.62   : Fractional growth-rate reduction at full 2-HG block

// ---------------- capecitabine / 5-FU ----------------
CL_FU    : 600    : 5-FU surrogate clearance (L/d)
V_FU     : 60     : 5-FU surrogate volume (L)
FCONV    : 0.030  : Capecitabine to systemic 5-FU conversion fraction

// ---------------- tumour ----------------
LAM0     : 0.0100 : Exponential growth rate (1/d)
LAM1     : 2.4    : Linear-phase growth ceiling (cm3/d)
PSI      : 20     : Simeoni exponential-to-linear switch sharpness (fixed a priori)
GP       : 0.55   : Persister growth rate relative to sensitive
KSP0     : 3.50   : Drug-driven entry into persistence [FIT]
KPS0     : 0.008  : Reversion out of persistence (1/d)
ALPHA_P  : 0.45   : Persister fraction of cytotoxic kill [FIT]
K_GEM    : 0.168  : Kill per AU dFdCTP per d [FIT]
K_CIS    : 0.192  : Kill per AU Pt-DNA adduct per d [FIT]
K_FU     : 0.100  : Kill per mg/L 5-FU per d
K_SYN    : 0.55   : Gemcitabine interference with platinum adduct repair
EFGFR    : 0.88   : Growth suppression at full FGFR block
KFGF_K   : 0.018  : Direct kill at full FGFR block (1/d) [FIT]
KIMM     : 0.055  : Immune kill scale (1/d) [FIT]
KIMM50   : 180    : Immune kill saturation volume (cm3)
KTR_D    : 0.55   : Damaged-cell transit rate (1/d)
MU_RES   : 1e-6   : Per-division rate to ANY FGFR2 kinase-domain mutation
KMET     : 8e-5   : Extrahepatic seeding rate (1/d per cm3)
LAMMET   : 0.020  : Metastatic growth rate (1/d)
K_2L     : 0      : Second-line cytotoxic kill after progression (1/d)
SL_LAG   : 21     : Days from progression to second line
SL_DUR   : 168    : Duration of second line (d)
T_PROG   : -1     : Progression time supplied by the driver (d; -1 = none yet)

// ---------------- stroma ----------------
KCAF     : 0.05   : CAF activation (1/d)
KCAF_OFF : 0.03   : CAF turnover (1/d)
KECM     : 0.030  : Matrix deposition (1/d)
KECM_OFF : 0.012  : Matrix turnover (1/d)
FPEN_MIN : 0.30   : Drug penetration floor at maximal desmoplasia
ECM50    : 1.0    : Matrix density halving penetration

// ---------------- immune ----------------
KREC     : 0.13   : CD8 recruitment (1/d)
KDEC     : 0.075  : CD8 decay (1/d)
KSUP     : 0.030  : Suppressor-mediated CD8 loss (1/d)
KICD     : 0.55   : Immunogenic-cell-death amplification
KSUPP_ON : 0.05   : Suppressor accumulation (1/d)
KSUPP_OF : 0.04   : Suppressor turnover (1/d)
PDL1EXP  : 1.0    : Tumour PD-L1 expression scale
KPD      : 0.06   : PD-L1 suppression half-effect
KIL6     : 1.0    : IL-6 half-effect on albumin
KIL6_ON  : 0.020  : IL-6 production (1/d)
KIL6_OFF : 0.09   : IL-6 clearance (1/d)

// ---------------- biliary ----------------
KOBS_EQ  : 0.30   : Obstruction equilibration (1/d)
KOBS     : 110    : Tumour volume for half-maximal duct occlusion (cm3)
DR_EFF   : 0.88   : Fraction of obstruction relieved by a patent stent
KOCC     : 0.0029 : Stent patency loss (1/d; 0.0029 SEMS 0.0077 plastic)
KING     : 0.9    : Tumour ingrowth acceleration of stent occlusion
KIN_BIL  : 0.15   : Bilirubin equilibration (1/d)
BILI0    : 0.6    : Unobstructed bilirubin (mg/dL)
EMAX_BIL : 24     : Maximal obstruction-driven bilirubin (mg/dL)
OBS50    : 0.55   : Obstruction fraction for half-maximal bilirubin
HBIL     : 3      : Hill coefficient obstruction to bilirubin
BIL_LIV  : 8      : Hepatocellular contribution to bilirubin (mg/dL)
KALP     : 0.10   : ALP equilibration (1/d)
ALP0     : 90     : Baseline ALP (U/L)
EMAX_ALP : 520    : Maximal obstruction-driven ALP (U/L)
KCH_ON   : 0.055  : Cholangitis generation (1/d)
KCH_OFF  : 0.16   : Cholangitis resolution (1/d)
INF_STEN : 1.4    : Occluded-stent amplification of cholangitis

// ---------------- liver ----------------
KFLR     : 0.06   : Functional reserve equilibration (1/d)
FOBS_FLR : 0.55   : Obstruction penalty on functional reserve
FLR_MIN  : 0.03   : Functional reserve floor
KALB     : 0.035  : Albumin turnover (1/d)
ALB_MAX  : 4.3    : Maximal albumin (g/dL)
ALB_FLR  : 0.35   : Obstruction-independent albumin floor fraction

// ---------------- myelosuppression (Friberg) ----------------
CIRC0    : 4.2    : Baseline ANC (10^9/L)
MTT      : 125    : Mean maturation time (h)
GAMMA    : 0.16   : Rebound feedback exponent
SLOPE_GE : 1.40   : Gemcitabine myelosuppression slope (per AU dFdCTP)
SLOPE_CI : 0.90   : Cisplatin myelosuppression slope (per AU Pt-DNA)
PLT0     : 250    : Baseline platelets (10^9/L)
KPLT     : 0.09   : Platelet turnover (1/d)
SLOPE_PL : 0.90   : Gemcitabine thrombocytopenia slope

// ---------------- renal / neuro / host ----------------
CRCL0    : 88     : Baseline creatinine clearance (mL/min)
KNEPH    : 0.055  : Platinum nephrotoxicity (per mg/L per d)
KCRCL_RE : 0.0016 : Renal recovery (1/d)
KNEURO   : 0.075  : Platinum neuropathy accumulation (per mg/L per d)
KNEURO_O : 0.0022 : Neuropathy resolution (1/d)
LBM0     : 42     : Baseline lean body mass (kg)
KCACH    : 0.0060 : Cachexia rate (1/d)
KLBM_REC : 0.0016 : Lean mass recovery (1/d)
KIL6C    : 1.4    : IL-6 half-effect on cachexia
KPS      : 0.10   : Performance-status equilibration (1/d)

// ---------------- biomarkers ----------------
KCA      : 0.20   : CA 19-9 equilibration (1/d)
CA_TUM   : 6.2    : CA 19-9 per cm3 tumour (U/mL)
CA_BILI  : 140    : CA 19-9 from cholestasis (U/mL) — the confound
CA0      : 22     : Baseline CA 19-9 (U/mL)
KCT      : 1.6    : ctDNA equilibration (1/d)
CT_SHED  : 0.020  : ctDNA shedding per cm3
KHG      : 1.1    : 2-HG equilibration (1/d)
HG0      : 1.0    : Baseline 2-HG (AU)
PHOS0    : 3.4    : Baseline serum phosphate (mg/dL)
KPHOS    : 0.28   : Phosphate equilibration (1/d)
EMAX_PHO : 3.0    : Maximal FGFR-inhibitor phosphate rise (mg/dL)
KIRAE_ON : 0.0016 : Immune-related AE accumulation (1/d)
KIRAE_OF : 0.020  : Immune-related AE resolution (1/d)

// ---------------- hazards ----------------
HT0      : 9.6e-5 : Baseline tumour hazard (1/d)
HT_T     : 0.0202 : Hazard per unit tumour/liver volume fraction (1/d) [FIT]
HT_FLR   : 0.0053 : Hazard per unit reserve deficit (1/d)
FLRCRIT  : 0.55   : Functional reserve below which hazard rises
HT_CACH  : 0.0038 : Hazard per unit fractional lean-mass loss (1/d)
HB0      : 6.4e-5 : Baseline biliary hazard (1/d)
HB_CH    : 0.0132 : Hazard per unit cholangitis burden (1/d) [FIT]
HB_ALBI  : 0.00053: Hazard per unit ALBI excess (1/d)
ALBI_REF : -2.60  : ALBI grade 1 upper bound

// ---------------- THE GATE ----------------
GATE_ON  : 1      : Master switch for dose-modification rules 1-4 (FALSIFIER)
BILI50   : 2.6    : Bilirubin at which half the dose is held (mg/dL)
HB_GATE  : 6      : Steepness of the bilirubin rule
ANC50    : 1.0    : ANC at which half the dose is held (10^9/L)
HA_GATE  : 8      : Steepness of the neutrophil rule
PS50     : 2.60   : ECOG at which half the dose is held
HP_GATE  : 8      : Steepness of the performance rule
CRCL50   : 45     : CrCl at which cisplatin is half omitted (mL/min)
HC_GATE  : 8      : Steepness of the renal rule

// ---------------- prescribed regimen (INTENTION, not delivery) ----------
GEM_MGM2 : 0      : Prescribed gemcitabine (mg/m2, days 1 and 8 of 21)
CIS_MGM2 : 0      : Prescribed cisplatin (mg/m2, days 1 and 8 of 21)
NCYCLE   : 0      : Number of 21-day gemcitabine/cisplatin cycles
DUR_MG   : 0      : Durvalumab dose (mg; q21d during chemo then q28d)
DUR_DAYS : 1460   : Durvalumab treatment duration (d)
FGI_MG   : 0      : Oral FGFR inhibitor daily dose (mg)
FGI_ON   : 14     : FGFR inhibitor days on
FGI_OFF  : 7      : FGFR inhibitor days off
FGI_ST   : 0      : FGFR inhibitor start day
IVO_MG   : 0      : Ivosidenib daily dose (mg)
CAP_MGM2 : 0      : Capecitabine dose (mg/m2 BID, 14 of 21)
CAP_DAYS : 0      : Capecitabine treatment duration (d)

$CMT @annotated
// --- PK ---
GEM_C   : Gemcitabine central (mg)
GEM_P   : Gemcitabine peripheral (mg)
DFDCTP  : Intracellular dFdCTP (AU)
CIS_C   : Free platinum central (mg)
CIS_P   : Free platinum peripheral (mg)
PTDNA   : Pt-DNA adducts (AU)
DUR_C   : Durvalumab central (mg)
DUR_P   : Durvalumab peripheral (mg)
FGI_A   : FGFR inhibitor depot (mg)
FGI_C   : FGFR inhibitor central (mg)
COV     : Covalently inactivated FGFR2 fraction
IVO_C   : Ivosidenib central (mg)
FU_C    : 5-FU surrogate central (mg)
// --- tumour ---
TS      : Drug-sensitive tumour clone (cm3)
TP      : Drug-tolerant persister clone (cm3)
TR      : FGFR2 kinase-domain-mutant clone (cm3)
TD1     : Damaged cells transit 1 (cm3)
TD2     : Damaged cells transit 2 (cm3)
TD3     : Damaged cells transit 3 (cm3)
TMET    : Extrahepatic metastatic burden (cm3)
NADIR   : RECIST nadir tracker (relative diameter)
// --- microenvironment ---
CAF     : Activated cancer-associated fibroblasts (AU)
ECM     : Desmoplastic matrix density (AU)
TEFF    : Intratumoural CD8 effectors (AU)
SUPP    : Suppressor compartment Treg/MDSC/TAM (AU)
IL6     : Tumour and stromal IL-6 (AU)
// --- biliary / hepatic ---
OBSR    : Biliary obstruction fraction
PATN    : Stent patency fraction
BILI    : Total bilirubin (mg/dL)
ALP     : Alkaline phosphatase (U/L)
ALB     : Serum albumin (g/dL)
FLR     : Functional liver reserve fraction
CHOLI   : Cholangitis burden (AU)
// --- host / toxicity ---
PROL    : Proliferating myeloid pool (10^9/L)
TRN1    : Neutrophil transit 1 (10^9/L)
TRN2    : Neutrophil transit 2 (10^9/L)
TRN3    : Neutrophil transit 3 (10^9/L)
ANC     : Circulating neutrophils (10^9/L)
PLT     : Platelets (10^9/L)
CRCL    : Creatinine clearance (mL/min)
NEURO   : Cumulative sensory neuropathy (AU)
LBM     : Lean body mass (kg)
PS      : ECOG performance status (continuous)
// --- biomarkers ---
CA199   : CA 19-9 (U/mL)
CTDNA   : ctDNA variant allele fraction (AU)
HG2     : Plasma 2-hydroxyglutarate (AU)
PHOS    : Serum phosphate (mg/dL)
IRAE    : Immune-related adverse event burden (AU)
// --- survival / bookkeeping ---
CUMHT   : Cumulative tumour hazard
CUMHB   : Cumulative biliary hazard
SURV    : Survival probability
CUMCIS  : Cumulative cisplatin dose (mg)
RDI     : Delivered dose accumulator
NDOSE   : Prescribed dose accumulator

$GLOBAL
#define INFW  (0.5/24.0)      // 30-minute chemotherapy infusion
#define DURW  (1.0/24.0)      // 1-hour durvalumab infusion
double _pos(double x){ return x > 0.0 ? x : 0.0; }
double _flr(double x, double lo){ return x < lo ? lo : x; }

$MAIN
if (NEWIND < 2) {
  // Goldie-Coldman seeding.  The resistant clone is PRESENT AT DIAGNOSIS.
  double frac_R = MU_RES * log(TT0 * 1e9);
  double tr0    = TT0 * frac_R;
  TS_0    = TT0 - tr0;
  TR_0    = tr0;
  NADIR_0 = 1.0;
  CAF_0   = 0.9;
  ECM_0   = 2.2;
  // A patient enrolled on a systemic-therapy trial has ALREADY been drained.
  // The driver sets PATN_0 = 0 for the delayed-drainage scenario.
  double obs0 = FHILAR * TT0 / (KOBS + TT0) * (1.0 - DR_EFF);
  PATN_0  = 1.0;
  OBSR_0  = obs0;
  BILI_0  = BILI0 + EMAX_BIL * pow(obs0, HBIL) / (pow(OBS50, HBIL) + pow(obs0, HBIL));
  ALP_0   = ALP0 + EMAX_ALP * obs0;
  FLR_0   = _flr((1.0 - TT0 / LIVVOL) * (1.0 - FOBS_FLR * obs0), FLR_MIN);
  ALB_0   = ALB_MAX * (ALB_FLR + (1.0 - ALB_FLR) * FLR_0);
  PROL_0  = CIRC0;  TRN1_0 = CIRC0;  TRN2_0 = CIRC0;  TRN3_0 = CIRC0;
  ANC_0   = CIRC0;  PLT_0  = PLT0;   CRCL_0 = CRCL0;  LBM_0  = LBM0;
  PS_0    = 1.0;
  CA199_0 = CA0 + CA_TUM * TT0 + CA_BILI * obs0 * obs0;
  CTDNA_0 = CT_SHED * TT0;
  HG2_0   = HG0 * (1.0 + 9.0 * IDH1);
  PHOS_0  = PHOS0;
  SURV_0  = 1.0;
}

$ODE
double t   = SOLVERTIME;
double vTS = _pos(TS), vTP = _pos(TP), vTR = _pos(TR);
double TT  = vTS + vTP + vTR;
double TTs = TT > 1e-9 ? TT : 1e-9;
double vTM = _pos(TMET);

double vECM  = _pos(ECM);
double fpen  = FPEN_MIN + (1.0 - FPEN_MIN) / (1.0 + vECM / ECM50);
double vBILI = _flr(BILI, 0.10);
double vALB  = _flr(ALB, 1.0);
double vANC  = _flr(ANC, 0.02);
double vPS   = _pos(PS);
double vCRCL = _flr(CRCL, 3.0);
double vFLR  = _flr(FLR, FLR_MIN);
double vPATN = _pos(PATN);
double vCHOL = _pos(CHOLI);
double vCAF  = _pos(CAF);
double vIL6  = _pos(IL6);
double vSUPP = _pos(SUPP);
double vTEFF = _pos(TEFF);
double albi  = 0.66 * log10(vBILI * 17.1) - 0.085 * (vALB * 10.0);

// ================= THE GATE — rules 1 to 4 ==================================
double g_chem = 1.0, g_cr = 1.0, g_ps = 1.0;
if (GATE_ON > 0.5) {
  double g_bil = 1.0 / (1.0 + pow(vBILI / BILI50, HB_GATE));
  double g_anc = 1.0 / (1.0 + pow(ANC50 / vANC,  HA_GATE));
  g_ps         = 1.0 / (1.0 + pow(vPS  / PS50,   HP_GATE));
  g_cr         = 1.0 / (1.0 + pow(CRCL50 / vCRCL, HC_GATE));
  g_chem = g_bil * g_anc * g_ps;
}

// ================= scheduled INTENTIONS =====================================
double pulse = 0.0, gem_rate = 0.0, cis_rate = 0.0;
if (NCYCLE > 0 && t < NCYCLE * 21.0) {
  double ph = fmod(t, 21.0);
  if (ph < INFW || (ph >= 8.0 && ph < 8.0 + INFW)) {
    pulse    = 1.0 / INFW;
    gem_rate = GEM_MGM2 * BSA * pulse * g_chem;          // DELIVERED, not prescribed
    cis_rate = CIS_MGM2 * BSA * pulse * g_chem * g_cr;
  }
}
double dur_rate = 0.0;
if (DUR_MG > 0.0 && t < DUR_DAYS) {
  double chend = NCYCLE * 21.0;
  double ph = (t < chend || NCYCLE == 0) ? fmod(t, 21.0) : fmod(t - chend, 28.0);
  if (ph < DURW) dur_rate = DUR_MG * (1.0 / DURW) * g_ps;
}
double fgi_in = 0.0;
if (FGI_MG > 0.0 && t >= FGI_ST) {
  if (fmod(t - FGI_ST, FGI_ON + FGI_OFF) < FGI_ON) fgi_in = FGI_MG;
}
double fu_in = 0.0;
if (CAP_MGM2 > 0.0 && t < CAP_DAYS && fmod(t, 21.0) < 14.0)
  fu_in = CAP_MGM2 * BSA * 2.0 * FCONV;

// ================= PK =======================================================
double Cg  = GEM_C / V1_GEM,  Cgp = GEM_P / V2_GEM;
dxdt_GEM_C = gem_rate - (CL_GEM / V1_GEM) * GEM_C - Q_GEM * (Cg - Cgp);
dxdt_GEM_P = Q_GEM * (Cg - Cgp);
dxdt_DFDCTP = VMAX_TP * fpen * Cg / (KM_TP + Cg) - KOUT_TP * DFDCTP;

double Cc  = CIS_C / V1_CIS,  Ccp = CIS_P / V2_CIS;
dxdt_CIS_C = cis_rate - (CL_CIS / V1_CIS) * CIS_C - Q_CIS * (Cc - Ccp);
dxdt_CIS_P = Q_CIS * (Cc - Ccp);
double tp_c = _pos(DFDCTP);
double krep = KREP_PT / (1.0 + K_SYN * tp_c / (1.0 + tp_c));   // gem blocks NER
dxdt_PTDNA = KFORM_PT * fpen * Cc - krep * PTDNA;

double Cd = DUR_C / V1_DUR, Cdp = DUR_P / V2_DUR;
dxdt_DUR_C = dur_rate - (CL_DUR / V1_DUR) * DUR_C - Q_DUR * (Cd - Cdp);
dxdt_DUR_P = Q_DUR * (Cd - Cdp);
double durocc = Cd / (Cd + IC50_DUR);

dxdt_FGI_A = fgi_in - KA_FGI * FGI_A;
double Cf = FGI_C / V_FGI;
dxdt_FGI_C = KA_FGI * FGI_A - (CL_FGI / V_FGI) * FGI_C;
double occ_rev = Cf / (Cf + IC50_FGI);
dxdt_COV = KCOV * occ_rev * (1.0 - COV) - KCOV_OFF * COV;
double occ_S = (COVAL > 0.5) ? _pos(COV) : occ_rev;
double occ_R = occ_S * RHO;                 // the gatekeeper mutation lives here

double Ci = IVO_C / V_IVO;
dxdt_IVO_C = IVO_MG - (CL_IVO / V_IVO) * IVO_C;
double occ_ivo = Ci / (Ci + IC50_IVO);

double Cu = FU_C / V_FU;
dxdt_FU_C = fu_in - (CL_FU / V_FU) * FU_C;

// ================= tumour ===================================================
double lam_mod = 1.0, lam_modR = 1.0;
if (FGFR2 > 0.5) { lam_mod *= (1.0 - EFGFR * occ_S); lam_modR *= (1.0 - EFGFR * occ_R); }
if (IDH1  > 0.5) { lam_mod *= (1.0 - EMAX_IDH * occ_ivo);
                   lam_modR *= (1.0 - EMAX_IDH * occ_ivo); }

double rr    = LAM0 * TTs / LAM1;
double denom = (rr > 4.0) ? rr : pow(1.0 + pow(rr, PSI), 1.0 / PSI);
double gS = LAM0 * lam_mod  * vTS / denom;
double gP = LAM0 * GP       * vTP / denom;
double gR = LAM0 * lam_modR * vTR / denom;

double chemo_hit = K_GEM * tp_c + K_CIS * _pos(PTDNA) + K_FU * Cu;
// Second line is not a scheduled input either: it starts a fixed lag after the
// progression the driver detected, and it passes through the same gate.
if (K_2L > 0.0 && T_PROG >= 0.0 && t >= T_PROG + SL_LAG && t < T_PROG + SL_LAG + SL_DUR)
  chemo_hit += K_2L * g_chem;

double fg_S = (FGFR2 > 0.5) ? KFGF_K * occ_S : 0.0;
double fg_R = (FGFR2 > 0.5) ? KFGF_K * occ_R : 0.0;
double imm_hit = KIMM * vTEFF * KIMM50 / (KIMM50 + TT);

double killS = (chemo_hit + fg_S + imm_hit) * vTS;
double killP = (ALPHA_P * (chemo_hit + fg_S) + imm_hit) * vTP;
double killR = (chemo_hit + fg_R + imm_hit) * vTR;
double e_sp  = KSP0 * (chemo_hit + fg_S);

dxdt_TS = gS - killS - e_sp * vTS + KPS0 * vTP;
dxdt_TP = gP - killP + e_sp * vTS - KPS0 * vTP;
dxdt_TR = gR - killR;
dxdt_TD1 = killS + killP + killR - KTR_D * TD1;
dxdt_TD2 = KTR_D * (TD1 - TD2);
dxdt_TD3 = KTR_D * (TD2 - TD3);
dxdt_TMET = KMET * TT + LAMMET * vTM - (chemo_hit + imm_hit) * vTM;

double sld = pow(TTs / TT0, 1.0 / 3.0);
dxdt_NADIR = (NADIR > sld) ? -4.0 * (NADIR - sld) : 0.0;

// ================= stroma ===================================================
dxdt_CAF = KCAF * TT / (TT + 60.0) - KCAF_OFF * vCAF;
dxdt_ECM = KECM * vCAF - KECM_OFF * vECM;
dxdt_IL6 = KIL6_ON * (vCAF + TT / 60.0 + 3.0 * vCHOL) - KIL6_OFF * vIL6;

// ================= immune ===================================================
double icd = KICD * (killS + killP + killR) / (5.0 + TT);
double supp_frac = (1.0 - durocc) * PDL1EXP / (PDL1EXP + KPD);
dxdt_TEFF = KREC * IMMENG * (1.0 + icd) * (1.0 - supp_frac)
            - KDEC * vTEFF - KSUP * vTEFF * vSUPP;
dxdt_SUPP = KSUPP_ON * (vCAF + vIL6) - KSUPP_OF * vSUPP;

// ================= biliary ==================================================
double obs_raw = FHILAR * TT / (KOBS + TT);
double obs_tgt = obs_raw * (1.0 - DR_EFF * vPATN);
double obsr = OBSR; if (obsr < 0.0) obsr = 0.0; if (obsr > 0.999) obsr = 0.999;
dxdt_OBSR = KOBS_EQ * (obs_tgt - OBSR);
dxdt_PATN = -KOCC * vPATN * (1.0 + KING * TT / (KOBS + TT));

double ob3 = pow(obsr, HBIL);
double bil_tgt = BILI0 + EMAX_BIL * ob3 / (pow(OBS50, HBIL) + ob3)
                 + BIL_LIV * ((vFLR > 0.45) ? 0.0 : (0.45 - vFLR) / 0.45);
dxdt_BILI = KIN_BIL * (bil_tgt - BILI);
dxdt_ALP  = KALP * (ALP0 + EMAX_ALP * obsr - ALP);

double neutropenic = 1.0 / (1.0 + pow(vANC / 0.5, 6.0));
dxdt_CHOLI = KCH_ON * obsr * obsr * (1.0 + INF_STEN * (1.0 - vPATN))
             * (1.0 + 0.6 * neutropenic) - KCH_OFF * vCHOL;

// ================= liver reserve (a FRACTION that is chased, not integrated) =
double flr_tgt = (1.0 - TT / LIVVOL) * (1.0 - FOBS_FLR * obsr);
if (flr_tgt < FLR_MIN) flr_tgt = FLR_MIN;
dxdt_FLR = KFLR * (flr_tgt - FLR);
dxdt_ALB = KALB * (ALB_MAX * (ALB_FLR + (1.0 - ALB_FLR) * vFLR)
                   * (1.0 - 0.32 * vIL6 / (vIL6 + KIL6)) - ALB);

// ================= myelosuppression (Friberg) ===============================
double ktr = 4.0 / (MTT / 24.0);
double edrug = SLOPE_GE * tp_c + SLOPE_CI * _pos(PTDNA);
if (edrug > 0.95) edrug = 0.95;
double fb = pow(CIRC0 / vANC, GAMMA);
if (fb > 3.0) fb = 3.0;
dxdt_PROL = ktr * PROL * ((1.0 - edrug) * fb - 1.0);
dxdt_TRN1 = ktr * (PROL - TRN1);
dxdt_TRN2 = ktr * (TRN1 - TRN2);
dxdt_TRN3 = ktr * (TRN2 - TRN3);
dxdt_ANC  = ktr * (TRN3 - ANC);
double plt_sup = SLOPE_PL * tp_c; if (plt_sup > 0.9) plt_sup = 0.9;
dxdt_PLT = KPLT * (PLT0 * (1.0 - plt_sup) - PLT);

// ================= renal / neuro / host =====================================
dxdt_CRCL   = -KNEPH * Cc * CRCL0 + KCRCL_RE * (CRCL0 - CRCL);
dxdt_NEURO  = KNEURO * Cc - KNEURO_O * NEURO;
dxdt_CUMCIS = cis_rate;
dxdt_LBM    = -KCACH * vIL6 / (vIL6 + KIL6C) * LBM + KLBM_REC * (LBM0 - LBM);

double bx = _pos(vBILI - 1.2);
double lbm_loss = _pos(1.0 - LBM / LBM0);
double ps_tgt = 0.6 + 1.7 * bx / (bx + 5.0) + 1.8 * vCHOL / (vCHOL + 0.8)
                + 3.3 * lbm_loss + 1.2 * (TT + vTM) / (TT + vTM + 900.0)
                + 0.9 * _pos(IRAE) + 0.5 * NEURO / (NEURO + 2.0);
if (ps_tgt > 4.0) ps_tgt = 4.0;
dxdt_PS = KPS * (ps_tgt - PS);

// ================= biomarkers ===============================================
dxdt_CA199 = KCA * (CA0 + CA_TUM * (TT + vTM) + CA_BILI * obsr * obsr - CA199);
dxdt_CTDNA = KCT * (CT_SHED * (TT + vTM) * (1.0 + 3.0 * vTR / TTs) - CTDNA);
dxdt_HG2   = KHG * (HG0 * (1.0 + 9.0 * IDH1 * (1.0 - 0.96 * occ_ivo)) - HG2);
dxdt_PHOS  = KPHOS * (PHOS0 + ((FGI_MG > 0.0) ? EMAX_PHO * occ_S : 0.0) - PHOS);
dxdt_IRAE  = ((DUR_MG > 0.0) ? KIRAE_ON * durocc : 0.0) - KIRAE_OF * IRAE;

// ================= the two hazards ==========================================
double ht = HT0 + HT_T * (TT + vTM) / LIVVOL + HT_FLR * _pos(FLRCRIT - vFLR)
            + HT_CACH * lbm_loss;
double hb = HB0 + HB_CH * vCHOL + HB_ALBI * _pos(albi - ALBI_REF);
dxdt_CUMHT = ht;
dxdt_CUMHB = hb;
dxdt_SURV  = -(ht + hb) * SURV;

dxdt_RDI   = pulse * g_chem * ((CIS_MGM2 > 0.0) ? g_cr : 1.0);
dxdt_NDOSE = pulse;

$TABLE
double TTOT   = _pos(TS) + _pos(TP) + _pos(TR);
double SLD    = pow((TTOT > 1e-9 ? TTOT : 1e-9) / TT0, 1.0 / 3.0);
double PERSFR = TTOT > 1e-9 ? _pos(TP) / TTOT : 0.0;
double RESFR  = TTOT > 1e-9 ? _pos(TR) / TTOT : 0.0;
double ALBI   = 0.66 * log10(_flr(BILI, 0.10) * 17.1) - 0.085 * (_flr(ALB, 1.0) * 10.0);
double ALBIG  = (ALBI <= -2.60) ? 1.0 : ((ALBI <= -1.39) ? 2.0 : 3.0);
double HAZT   = HT0 + HT_T * (TTOT + _pos(TMET)) / LIVVOL
                + HT_FLR * _pos(FLRCRIT - _flr(FLR, FLR_MIN))
                + HT_CACH * _pos(1.0 - LBM / LBM0);
double HAZB   = HB0 + HB_CH * _pos(CHOLI) + HB_ALBI * _pos(ALBI - ALBI_REF);
double HAZFR  = (HAZT + HAZB) > 0 ? HAZB / (HAZT + HAZB) : 0.0;
double GBILI  = (GATE_ON > 0.5) ? 1.0 / (1.0 + pow(_flr(BILI,0.10) / BILI50, HB_GATE)) : 1.0;
double GANC   = (GATE_ON > 0.5) ? 1.0 / (1.0 + pow(ANC50 / _flr(ANC,0.02), HA_GATE)) : 1.0;
double GPS    = (GATE_ON > 0.5) ? 1.0 / (1.0 + pow(_pos(PS) / PS50, HP_GATE)) : 1.0;
double GCRCL  = (GATE_ON > 0.5) ? 1.0 / (1.0 + pow(CRCL50 / _flr(CRCL,3.0), HC_GATE)) : 1.0;
double GATE   = GBILI * GANC * GPS;
double RDI_PC = (NDOSE > 1e-9) ? 100.0 * RDI / NDOSE : 100.0;

$CAPTURE @annotated
TTOT   : Total tumour volume (cm3)
SLD    : RECIST relative sum of diameters
PERSFR : Persister fraction of the tumour
RESFR  : FGFR2-mutant clone fraction of the tumour
ALBI   : ALBI score
ALBIG  : ALBI grade
HAZT   : Instantaneous tumour hazard (1/d)
HAZB   : Instantaneous biliary hazard (1/d)
HAZFR  : Fraction of the total hazard that is biliary
GBILI  : Rule 1 gate (bilirubin)
GANC   : Rule 2 gate (neutrophils)
GPS    : Rule 4 gate (performance status)
GCRCL  : Rule 3 gate (renal, cisplatin only)
GATE   : Combined chemotherapy gate
RDI_PC : Delivered relative dose intensity (%)
'

mod <- mcode_cache("cca_qsp", code)

# =============================================================================
#  SCENARIO DRIVERS
#  Each driver builds an event object for the DISCRETE interventions only
#  (stent placement and exchange).  The cytotoxic, oral and antibody schedules
#  live inside $ODE because they must be multiplied by the gate at the exact
#  moment of administration — which is the whole argument of this model.
# =============================================================================

# --- stent events -----------------------------------------------------------
# mrgsolve has no state-triggered events, so exchanges are scheduled here at
# the interval at which the corresponding stent type actually occludes.
# The Python twin instead re-stents when patency crosses 0.25; the two agree
# to within a few percent on every endpoint in README.md.
stent_events <- function(kind = c("sems", "plastic", "none"),
                         first_day = 0, tend = 1100) {
  kind <- match.arg(kind)
  if (kind == "none") return(NULL)
  every <- if (kind == "sems") 240 else 90
  days  <- seq(first_day, tend, by = every)
  ev(time = days, cmt = "PATN", amt = 1, evid = 8)   # evid 8 = replace state
}

sim_arm <- function(par = list(), stent = "sems", stent_day = 0,
                    tend = 1100, delta = 1, ...) {
  m <- mod
  if (length(par)) m <- param(m, par)
  if (stent_day > 0) m <- init(m, PATN = 0)   # start jaundiced (S4)
  e <- stent_events(stent, stent_day, tend)
  if (is.null(e)) mrgsim(m, end = tend, delta = delta, ...)
  else            mrgsim(m, events = e, end = tend, delta = delta, ...)
}

# --- S1  best supportive care + drainage ------------------------------------
run_bsc <- function(tend = 1100)
  sim_arm(list(FHILAR = 0.72), stent = "sems", tend = tend)

# --- S2  gemcitabine/cisplatin (ABC-02) -------------------------------------
run_gemcis <- function(tend = 1100)
  sim_arm(list(FHILAR = 0.72, GEM_MGM2 = 1000, CIS_MGM2 = 25, NCYCLE = 8,
               K_2L = 0.010), stent = "sems", tend = tend)

# --- S3  gem/cis + durvalumab (TOPAZ-1) -------------------------------------
# IMMENG is a MIXTURE INDICATOR, not a dose.  Run the arm twice — engaged and
# not engaged — and weight by PI_IMMUNE.  Averaging the two S(t) curves is
# what produces a tail separation with almost no median separation.
run_topaz1 <- function(pi_immune = 0.22, tend = 1100) {
  p <- list(FHILAR = 0.72, GEM_MGM2 = 1000, CIS_MGM2 = 25, NCYCLE = 8,
            DUR_MG = 1500, K_2L = 0.010)
  hot  <- sim_arm(c(p, list(IMMENG = 1)), tend = tend)
  cold <- sim_arm(c(p, list(IMMENG = 0)), tend = tend)
  list(hot = hot, cold = cold, pi = pi_immune,
       surv = pi_immune * hot$SURV + (1 - pi_immune) * cold$SURV)
}

# --- S4  the same regimen with drainage delayed 60 days ---------------------
#     The ONLY difference from S3 is when the stent goes in.  No antitumour
#     parameter changes anywhere.
run_delayed_drainage <- function(delay = 60, tend = 1100)
  sim_arm(list(FHILAR = 0.72, GEM_MGM2 = 1000, CIS_MGM2 = 25, NCYCLE = 8,
               DUR_MG = 1500, K_2L = 0.010),
          stent = "sems", stent_day = delay, tend = tend)

# --- S6/S7  FGFR inhibitors --------------------------------------------------
run_pemigatinib <- function(tend = 1100)
  sim_arm(list(FHILAR = 0.16, FGFR2 = 1, FGI_MG = 13.5, FGI_ON = 14,
               FGI_OFF = 7, COVAL = 0, RHO = 0.05), tend = tend)

run_futibatinib <- function(tend = 1100)
  sim_arm(list(FHILAR = 0.16, FGFR2 = 1, FGI_MG = 20, FGI_ON = 1, FGI_OFF = 0,
               COVAL = 1, RHO = 0.45), tend = tend)

# --- S11  a PREDICTION: what the 1-week pemigatinib holiday costs ------------
run_pemigatinib_continuous <- function(tend = 1100)
  sim_arm(list(FHILAR = 0.16, FGFR2 = 1, FGI_MG = 13.5, FGI_ON = 1, FGI_OFF = 0,
               COVAL = 0, RHO = 0.05), tend = tend)

# --- S8  ivosidenib ---------------------------------------------------------
run_ivosidenib <- function(tend = 1100)
  sim_arm(list(FHILAR = 0.16, IDH1 = 1, IVO_MG = 500), tend = tend)

# --- S9/S10  adjuvant --------------------------------------------------------
run_bilcap <- function(cape = TRUE, tend = 2400) {
  p <- list(FHILAR = 0.72, TT0 = 100)
  s <- sim_arm(if (cape) c(p, list(CAP_MGM2 = 1250, CAP_DAYS = 180)) else p,
               stent = "sems", tend = tend)
  # resection is an initial condition, not a rate: 99.99% of the burden is gone
  s
}
run_resected <- function(cape = TRUE, tend = 2400) {
  m <- param(mod, list(FHILAR = 0.72,
                       CAP_MGM2 = if (cape) 1250 else 0,
                       CAP_DAYS = if (cape) 180 else 0))
  m <- init(m, TS = 0.01, TR = 0)     # microscopic residual disease
  mrgsim(m, events = stent_events("sems", 0, tend), end = tend, delta = 7)
}

# =============================================================================
#  FALSIFIERS — one parameter each, nothing else refitted
# =============================================================================
# F1  the gate is the mechanism by which drainage matters.  Switch it off and
#     delaying drainage by 60 days must stop costing survival.
run_F1 <- function(tend = 1100)
  sim_arm(list(FHILAR = 0.72, GEM_MGM2 = 1000, CIS_MGM2 = 25, NCYCLE = 8,
               DUR_MG = 1500, K_2L = 0.010, GATE_ON = 0),
          stent = "sems", stent_day = 60, tend = tend)

# F2  remove the pre-existing resistant clone entirely.
run_F2 <- function(tend = 1100)
  sim_arm(list(FHILAR = 0.16, FGFR2 = 1, FGI_MG = 13.5, MU_RES = 0,
               COVAL = 0, RHO = 0.05), tend = tend)

# F3  make every tumour immune-engaged: durvalumab must then move the MEDIAN,
#     which TOPAZ-1 says it does not.
run_F3 <- function(tend = 1100)
  sim_arm(list(FHILAR = 0.72, GEM_MGM2 = 1000, CIS_MGM2 = 25, NCYCLE = 8,
               DUR_MG = 1500, IMMENG = 1, K_2L = 0.010), tend = tend)

# F4  delete the stromal penetration barrier: gem/cis ORR must overshoot.
run_F4 <- function(tend = 1100)
  sim_arm(list(FHILAR = 0.72, GEM_MGM2 = 1000, CIS_MGM2 = 25, NCYCLE = 8,
               FPEN_MIN = 1.0, K_2L = 0.010), tend = tend)

# =============================================================================
#  VIRTUAL POPULATION
#  Population survival is the MEAN of individual S(t).  Median OS is where the
#  mean curve crosses 0.5; the 24-month tail is that same curve at t = 730.
#  Both come out of one simulation, which is the point of commitment 3.
# =============================================================================
run_population <- function(base = list(), n = 120, pi_immune = 0.22,
                           tend = 1100, seed = 1104) {
  set.seed(seed)
  ln <- function(cv) exp(rnorm(n, -0.5 * log(1 + cv^2), sqrt(log(1 + cv^2))))
  idata <- data.frame(
    ID     = seq_len(n),
    LAM0   = mod@param$LAM0   * ln(0.30),
    FHILAR = pmin(0.97, pmax(0.04, (base$FHILAR %||% mod@param$FHILAR) * ln(0.22))),
    KOBS   = mod@param$KOBS   * ln(0.30),
    HT_T   = mod@param$HT_T   * ln(0.35),
    HB_CH  = mod@param$HB_CH  * ln(0.40),
    K_GEM  = mod@param$K_GEM  * ln(0.42),
    K_CIS  = mod@param$K_CIS  * ln(0.42),
    KFGF_K = mod@param$KFGF_K * ln(0.42),
    CIRC0  = mod@param$CIRC0  * ln(0.20),
    CRCL0  = mod@param$CRCL0  * ln(0.22),
    IMMENG = as.numeric(runif(n) < pi_immune),
    TT0    = 100 * ln(0.55)
  )
  m <- param(mod, base[setdiff(names(base), "FHILAR")])
  out <- mrgsim(m, idata = idata, events = stent_events("sems", 0, tend),
                end = tend, delta = 7)
  as.data.frame(out) %>%
    group_by(time) %>%
    summarise(S = mean(SURV), TT = median(TTOT), RDI = mean(RDI_PC),
              .groups = "drop")
}
`%||%` <- function(a, b) if (is.null(a)) b else a

median_os <- function(pop) {
  i <- which(pop$S <= 0.5)[1]
  if (is.na(i)) return(NA_real_)
  x0 <- pop$time[i - 1]; y0 <- pop$S[i - 1]; y1 <- pop$S[i]
  (x0 + (y0 - 0.5) / (y0 - y1) * (pop$time[i] - x0)) / 30.44
}
os_at <- function(pop, days) approx(pop$time, pop$S, days)$y * 100

if (interactive()) {
  message("Cholangiocarcinoma QSP model — 54 ODEs")
  message("Virtual trials: source cca_reference_impl.py for the calibrated numbers.")
  ctrl <- run_population(list(GEM_MGM2 = 1000, CIS_MGM2 = 25, NCYCLE = 8, K_2L = 0.010))
  message(sprintf("gem/cis  median OS %.1f mo   24-mo OS %.1f%%",
                  median_os(ctrl), os_at(ctrl, 730)))
}
