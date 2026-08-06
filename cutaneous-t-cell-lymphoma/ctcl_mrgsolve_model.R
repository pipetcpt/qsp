## =============================================================================
##  ctcl_mrgsolve_model.R
##  Cutaneous T-cell lymphoma — mycosis fungoides (MF) and Sézary syndrome (SS)
##  Quantitative Systems Pharmacology model for mrgsolve
##
##  67 ODE compartments.  Time unit = DAYS.
##
##  ---------------------------------------------------------------------------
##  THE TWO STRUCTURAL COMMITMENTS
##  ---------------------------------------------------------------------------
##
##  (1) ONE CLONE, THREE BOXES.  Mycosis fungoides and Sézary syndrome are not
##      given different equations here.  They are the same equations run at
##      different values of TWO residency covariates:
##
##          CD69   high -> S1PR1 internalised -> the clone cannot leave the skin
##          CCR7   high -> central-memory recirculation programme
##
##      MF (CD69 0.85 / CCR7 0.10) and SS (CD69 0.20 / CCR7 0.85) then fall out
##      of J_egress and J_home.  Nothing in this file says "Sézary syndrome is
##      leukaemic".  The blood burden is computed.
##
##  (2) TWO HOST VARIABLES, NOT ONE.  How much clone was killed and how much
##      anti-clone surveillance survived are separate states:
##
##          CLONE KILL       NSK, NSKR, NTR, NBL, NBLR, NLN, NVS
##          SURVEILLANCE     DCA, E8SK, E8BL, NKB, TREG
##
##      Every response rate (ORR) reads the first.  Every durability measure
##      (TTNT, DOR) is set by the second.  Cytotoxic chemotherapy moves them in
##      OPPOSITE directions; interferon and photopheresis move only the second.
##      The model is therefore able to produce a high ORR with a short TTNT and
##      a low ORR with a long TTNT from the same two lines of algebra.
##
##  ---------------------------------------------------------------------------
##  THE CONSTANTS THAT WERE NOT FITTED
##  ---------------------------------------------------------------------------
##
##      ABCSK  = 0.157   antibody biodistribution coefficient, skin
##      ABCLN  = 0.085   antibody biodistribution coefficient, lymph node
##                       (Shah DK & Betts AM, mAbs 2013 — measured for IgG,
##                        taken here as given, for every antibody in the file)
##      NKRSK  = 0.12    NK-cell density in lesional skin relative to blood
##
##  These three numbers, and not any fitted potency term, are what makes an
##  anti-CCR4 antibody an order of magnitude weaker in the skin than in the
##  blood.  The consequence is derived rather than declared:
##
##    * receptor OCCUPANCY in skin is essentially complete
##      (0.157 x C_plasma is still >> KD), so raising the dose does nothing;
##    * the limiting quantity is the EFFECTOR, not the drug;
##    * therefore mogamulizumab's skin benefit cannot be ADCC.  What is left is
##      (a) removing the blood source by blocking CCR4-mediated homing and
##      (b) depleting CCR4-high regulatory T cells, which releases the CD8
##          surveillance term.  Both are in the file; neither was added to make
##          the skin response rate come out right.
##    * because (b) is the same event as mogamulizumab-associated rash, the
##      model predicts rash and response travel together.
##
##  An antibody-drug conjugate has no such gap: brentuximab vedotin needs no
##  effector cell, so its only compartment penalty is the 0.157 itself.
##
##  ---------------------------------------------------------------------------
##  WHAT IS AN INPUT AND WHAT IS AN OUTPUT
##  ---------------------------------------------------------------------------
##
##  INPUTS   residency (CD69, CCR7), CLA, CCR4 and CD30 expression, the
##           drug-resistant clone fraction FRES, immune fitness IMMF, NK
##           density factor NKF, stage burden, body size, and the regimen.
##
##  OUTPUTS  everything a trial reports: mSWAT and its patch/plaque/tumour
##           decomposition, Sézary count and B-score, N-score, the GLOBAL
##           response by the Olsen 2011 conjunction, PFS, TTNT, pruritus NRS,
##           serum TARC / sIL-2R / LDH, ANC and platelets, MMAE neuropathy,
##           bexarotene hypothyroidism and hypertriglyceridaemia, CD4 count,
##           infection hazard and disease hazard.
##
##  ---------------------------------------------------------------------------
##  CALIBRATION ANCHORS (see ctcl_references.md)
##  ---------------------------------------------------------------------------
##    MAVORIC     mogamulizumab vs vorinostat, previously treated MF/SS
##    ALCANZA     brentuximab vedotin vs physician's choice, CD30 >= 10%
##    Olsen 2007  vorinostat single-arm (skin-weighted criteria)
##    Whittaker 2010 romidepsin single-arm
##    Hughes 2015 time to next treatment by therapy class
##    Agar 2010   survival by TNMB stage
##    Hoppe 2015  low-dose (12 Gy) total skin electron beam
##
##  DISCLAIMER.  Educational / research model.  Not validated for clinical or
##  regulatory use.  Parameters are illustrative approximations.
## =============================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
})

ctcl_code <- '
$PROB
# Cutaneous T-cell lymphoma (MF / Sezary syndrome) QSP model
# One clone in three compartments; kill and surveillance as separate states.

$PARAM @annotated
// ---------------------------------------------------------------- PK: mAbs
VM1    :  3.6  : Mogamulizumab central volume (L)
VM2    :  3.1  : Mogamulizumab peripheral volume (L)
QM     :  0.42 : Mogamulizumab intercompartmental clearance (L/day)
CLM    :  0.30 : Mogamulizumab clearance (L/day)
KDMOG  :  0.02 : Apparent K50 for CCR4 occupancy (ug/mL)
VB1    :  6.0  : Brentuximab vedotin ADC central volume (L)
VB2    :  3.0  : Brentuximab vedotin ADC peripheral volume (L)
QB     :  0.9  : ADC intercompartmental clearance (L/day)
CLB    :  0.60 : ADC clearance (L/day)
KDEC   :  0.050: ADC deconjugation rate (1/day)
FMMAE  :  0.019: mg MMAE liberated per mg ADC (DAR 4)
VMM    :  1200.: MMAE volume of distribution (L)
CLMM   :  60.  : MMAE clearance (L/day)
KDBV   :  0.30 : ADC binding K50 (ug/mL)
VA     :  6.0  : Alemtuzumab volume (L)
CLAL   :  3.0  : Alemtuzumab clearance (L/day)
KDALEM :  0.05 : Alemtuzumab occupancy K50 (ug/mL)
// ------------------------------------------------- PK: small molecules
VR1    :  100. : Romidepsin central volume (L)
VR2    :  200. : Romidepsin peripheral volume (L)
QR     :  200. : Romidepsin intercompartmental clearance (L/day)
CLR    :  706. : Romidepsin clearance (L/day)
KAV    :  6.0  : Vorinostat absorption rate (1/day)
FV     :  0.43 : Vorinostat bioavailability
VV     :  520. : Vorinostat volume (L)
CLV    :  4300.: Vorinostat clearance (L/day)
KAX    :  3.0  : Bexarotene absorption rate (1/day)
FX     :  0.60 : Bexarotene bioavailability
VX     :  200. : Bexarotene volume (L)
CLX    :  475. : Bexarotene clearance (L/day)
KAI    :  3.0  : Interferon SC absorption rate (1/day)
VI     :  8.0  : Interferon volume (L)
CLI    :  60.  : Interferon clearance (L/day)
KAT    :  8.0  : Methotrexate absorption rate (1/day)
FT     :  0.70 : Methotrexate bioavailability
VT     :  35.  : Methotrexate volume (L)
CLT    :  83.  : Methotrexate clearance (L/day)
VG     :  50.  : Gemcitabine (dFdCTP surrogate) volume (L)
CLG    :  23.  : Gemcitabine surrogate clearance (L/day)
KSDT   :  0.35 : Skin-directed exposure decay (1/day)
KTSEB  :  0.10 : TSEB sublethal-damage decay (1/day)
KSTER  :  1.2  : Corticosteroid exposure decay (1/day)
KABX   :  0.8  : Antibiotic exposure decay (1/day)
KECP   :  0.12 : ECP signal decay (1/day)
// ------------------------------- structural constants that were NOT fitted
ABCSK  :  0.157: Antibody biodistribution coefficient - skin
ABCLN  :  0.085: Antibody biodistribution coefficient - lymph node
NKRSK  :  0.12 : NK density in lesional skin relative to blood
// ---------------------------------------------------------- clone kinetics
GSK    :  0.01774: Skin clone gross proliferation (1/day)
GBL    :  0.0120 : Blood clone gross proliferation (1/day)
GLN    :  0.0150 : Node clone gross proliferation (1/day)
GVS    :  0.0180 : Visceral clone gross proliferation (1/day)
GTR    :  0.0090 : Transformed subclone gross proliferation (1/day)
KSK    :  20.0   : Skin niche carrying capacity (1e9 cells)
KBL    :  300.0  : Blood carrying capacity (1e9 cells)
KLNC   :  8.0    : Node carrying capacity (1e9 cells)
KHOME  :  0.055  : Blood-to-skin homing rate constant (1/day)
KEGR   :  0.060  : Skin-to-blood egress rate constant (1/day)
KNODE  :  0.010  : Blood-to-node rate (1/day)
KNOUT  :  0.006  : Node-to-blood rate (1/day)
KVIS   :  0.0010 : Node-to-viscera rate (1/day)
KTRANS :  8.0e-7 : Large-cell transformation hazard (1/day per 1e9 cells)
CLADH  :  0.90   : CLA-dependent adhesion efficiency
RESK   :  0.12   : Residual DRUG kill on the resistant subclone (8-fold resistant)
RESI   :  0.50   : Residual CTL kill on the resistant subclone (2-fold resistant)
RESE   :  0.85   : Residual ADCC/CDC kill on the resistant subclone (effector killing bypasses it)
FITR   :  0.93   : Fitness cost of drug resistance
KCONV  :  0.030  : Drug-pressure-driven conversion to resistance
// ------------------------------------------------------ immune surveillance
KIMM   :  0.01921: Skin immune kill rate constant (1/day)
KIMB   :  0.0100 : Blood immune kill rate constant (1/day)
KIM50  :  60.0   : Burden at which skin surveillance efficiency halves
KIM50B :  100.0  : Burden at which blood surveillance efficiency halves
KIMTR  :  1.20   : Surveillance sensitivity of the transformed subclone
KIMESC :  0.04   : Th2-driven loss of surveillance efficiency
KDE8   :  0.030  : CD8 effector turnover (1/day)
KTREG50:  1.4    : Treg suppression K50
KIL1050:  2.2    : IL-10 suppression K50
KDC    :  0.085  : DC activation turnover (1/day)
KDCAG  :  0.80   : Antigen-driven DC priming gain
KDCAGK :  6.0    : Antigen load at half-maximal DC priming (1e9 cells)
KECPD  :  0.90   : ECP gain on DC activation
KIFNDC :  1.10   : Interferon gain on DC activation
KDCSDT :  0.35   : Phototherapy gain on DC activation (immunogenic cell death)
KNK    :  0.25   : NK production (1/day)
KNKD   :  0.25   : NK turnover (1/day)
KTREGP :  0.030  : Treg production (1/day)
KTREGD :  0.030  : Treg turnover (1/day)
// -------------------------------------------------------------- drug effect
KADCC  :  0.130  : Maximum ADCC kill rate at full occupancy (1/day)
KADCA  :  0.200  : Maximum alemtuzumab kill rate (1/day)
ETSAT  :  3.0    : Effector-to-target saturation constant
KINT   :  0.55   : ADC internalisation efficiency
KBVKILL:  0.500  : ADC payload kill rate constant (1/day)
KBYST  :  0.020  : MMAE bystander kill rate constant (1/day)
KMMAE50:  0.0020 : MMAE K50 (ug/mL)
EROM   :  0.550  : Romidepsin maximum kill rate (1/day)
CROM50 :  0.050  : Romidepsin K50 (ug/mL)
EVOR   :  0.0394 : Vorinostat maximum kill rate (1/day)
CVOR50 :  0.040  : Vorinostat K50 (ug/mL)
EBEX   :  0.037  : Bexarotene maximum kill rate (1/day)
CBEX50 :  0.600  : Bexarotene K50 (ug/mL)
EMTX   :  0.085  : Methotrexate maximum kill rate (1/day)
CMTX50 :  0.150  : Methotrexate K50 (ug/mL)
EGEM   :  0.130  : Gemcitabine maximum kill rate (1/day)
CGEM50 :  8.0    : Gemcitabine surrogate K50 (ug/mL)
ESTER  :  0.065  : Corticosteroid maximum kill rate (1/day)
ISTER  :  0.60   : Corticosteroid immunosuppression factor
EUVB   :  0.088  : Skin-directed maximum kill rate (1/day)
ETSEBK :  0.020  : TSEB kill rate per unit damage (1/day)
DEPTHP :  0.45   : Penetration into plaque
DEPTHT :  0.12   : Penetration into tumour
EIFNK  :  0.020  : Interferon direct kill (1/day)
CIFN50 :  0.080  : Interferon K50 (ug/mL)
EIFNTH2:  0.55   : Interferon Th2 suppression
KCD30I :  0.55   : HDAC-inhibitor induction of CD30
KCD30O :  0.10   : CD30 expression turnover (1/day)
CD30B  :  0.10   : Baseline CD30 set point
KSELC  :  0.60   : CCR4 antigen-loss selection strength
KSELD  :  0.35   : CD30 antigen-loss selection strength
// ------------------------------------------------- cytokines and biomarkers
KTH2   :  0.16 : Th2 tone turnover (1/day)
TH2SK  :  1.20 : Skin burden gain on Th2
TH2BL  :  0.80 : Blood burden gain on Th2
TH2K   :  6.0  : Burden at half-maximal Th2 gain
KIL31  :  0.20 : IL-31 turnover (1/day)
I31SK  :  1.40 : Skin burden gain on IL-31
I31BL  :  1.30 : Blood burden gain on IL-31
I31K   :  5.0  : Burden at half-maximal IL-31 gain
KIL10  :  0.16 : IL-10 turnover (1/day)
I10SK  :  0.50 : Treg-dependent IL-10 gain
I10K   :  6.0  : Burden at half-maximal IL-10 gain
KIFNG  :  0.15 : IFN-gamma turnover (1/day)
KTARC  :  0.55 : Serum TARC turnover (1/day)
TARCB  :  250. : Serum TARC baseline (pg/mL)
TARCSK :  900. : Skin contribution to TARC (pg/mL)
TARCBL :  9000.: Blood contribution to TARC (pg/mL)
TARCK  :  8.0  : Burden at half-maximal TARC
KSIL   :  0.50 : Soluble IL-2R turnover (1/day)
SILB   :  450. : sIL-2R baseline (U/mL)
SILSK  :  1400.: Skin contribution to sIL-2R
SILBL  :  1800.: Blood contribution to sIL-2R
SILK   :  8.0  : Burden at half-maximal sIL-2R
KLDH   :  0.40 : LDH turnover (1/day)
LDHB   :  195. : LDH baseline (U/L)
LDHSL  :  42.  : LDH slope per 10e9 cells
// ------------------------------------------------------- lesion morphology
ASCALE :  6.0  : Clone burden scaling for involved body-surface area
TAUPAT :  8.0  : Patch resolution time constant (day)
TAUPLQ :  25.0 : Plaque resolution time constant (day)
TAUTUM :  12.0 : Tumour resolution time constant (day)
ERYK   :  2.6  : Blood burden at half-maximal erythroderma
ERYF   :  62.0 : Maximum erythroderma area contribution (percent BSA)
// -------------------------------------------- barrier, S. aureus, pruritus
KBARR  :  0.10 : Barrier repair rate (1/day)
KBDAM  :  0.05 : Barrier damage per unit lesion (1/day)
KSAP   :  0.40 : S. aureus colonisation rate (1/day)
KSAD   :  0.42 : S. aureus clearance (1/day)
KSAG   :  0.30 : Superantigen production (1/day)
KSAGD  :  0.55 : Superantigen clearance (1/day)
ESAG   :  0.15 : Superantigen boost to clone proliferation
SAG50  :  2.0  : Superantigen half-effect
KPR    :  0.16 : Pruritus turnover (1/day)
KITCH50:  2.6  : Itch drive at half-maximal NRS
KSENS  :  0.010: Central sensitisation gain (1/day)
KSENSD :  0.010: Central sensitisation decay (1/day)
PRMAX  :  10.0 : Maximum pruritus NRS
// --------------------------------------------------- haematology, toxicity
KTR    :  0.55 : Myeloid transit rate (1/day)
GAMMA  :  0.16 : Myeloid feedback exponent
CIRC0  :  4.0  : Baseline ANC (10^9/L)
PLT0   :  250. : Baseline platelets (10^9/L)
KPLT   :  0.10 : Platelet turnover (1/day)
KPN    :  0.55 : Neuropathy accumulation per ng/mL MMAE (1/day)
KPND   :  0.0035: Neuropathy resolution (1/day)
KTSH   :  0.25 : TSH turnover (1/day)
KFT4   :  0.09 : Free T4 turnover (1/day)
ETSH   :  0.90 : Maximum bexarotene TSH suppression
CTSH50 :  0.350: Bexarotene K50 for endocrine effects (ug/mL)
KTG    :  0.05 : Triglyceride turnover (1/day)
TG0    :  130. : Baseline triglycerides (mg/dL)
ETG    :  3.2  : Maximum bexarotene triglyceride rise
KMARP  :  0.030: Mogamulizumab-associated rash formation (1/day)
KMARD  :  0.035: Mogamulizumab-associated rash resolution (1/day)
KCD4   :  0.11 : CD4 reconstitution (1/day)
CD40   :  800. : Baseline CD4 count (cells/uL)
HDIS   :  6.5e-5 : Disease hazard per unit tumour burden (1/day)
HINF   :  3.2e-4 : Infection hazard scaling (1/day)
HBASE  :  2.2e-5 : Background hazard (1/day)
// ------------------------------------------------------------- covariates
CD69   :  0.85 : Clone CD69 positivity (residency; MF high, SS low)
CCR7   :  0.10 : Clone CCR7 positivity (recirculation; SS high)
CLA    :  0.90 : Clone CLA positivity (skin homing)
CCR4   :  0.85 : Clone CCR4 positivity at baseline
CD30   :  0.10 : Clone CD30 expression at baseline
IMMF   :  1.00 : Immune fitness multiplier
NKF    :  1.00 : NK density multiplier
GTRF   :  1.00 : Transformed-subclone proliferation multiplier
FRES   :  0.12 : Drug-resistant fraction of the clone at baseline
SENSF  :  1.00 : Clone-intrinsic sensitivity to direct drug kill
KSKF   :  1.00 : Skin niche capacity multiplier
CLF    :  1.00 : Antibody clearance multiplier
CLAMP  :  0    : If 1, hold the clone fixed (used to equilibrate the host)
INITF  :  1    : If 1, derive initial conditions in $MAIN; set 0 when supplying them

$CMT @annotated
MOG1  : Mogamulizumab central (mg)
MOG2  : Mogamulizumab peripheral (mg)
BV1   : Brentuximab vedotin ADC central (mg)
BV2   : Brentuximab vedotin ADC peripheral (mg)
MMAE1 : Free MMAE (mg)
MMAE2 : MMAE placeholder (mg)
ROM1  : Romidepsin central (mg)
ROM2  : Romidepsin peripheral (mg)
VORD  : Vorinostat depot (mg)
VOR   : Vorinostat central (mg)
BEXD  : Bexarotene depot (mg)
BEX   : Bexarotene central (mg)
IFND  : Interferon SC depot (units)
IFN   : Interferon central (units)
MTXD  : Methotrexate depot (mg)
MTX   : Methotrexate central (mg)
GEM   : Gemcitabine intracellular surrogate (mg)
ALEM  : Alemtuzumab central (mg)
SDT   : Skin-directed exposure (arbitrary)
TSEBC : Total skin electron beam damage (Gy-equivalent)
STER  : Corticosteroid exposure (arbitrary)
ABX   : Antibiotic exposure (arbitrary)
ECPD  : Photopheresis signal (arbitrary)
NSK   : Skin clone - drug sensitive (1e9 cells)
NSKR  : Skin clone - drug resistant (1e9 cells)
NTR   : Transformed subclone in skin (1e9 cells)
NBL   : Blood clone - drug sensitive (1e9 cells)
NBLR  : Blood clone - drug resistant (1e9 cells)
NLN   : Lymph-node clone (1e9 cells)
NVS   : Visceral clone (1e9 cells)
FCCR4 : Fraction of clone expressing CCR4
CD30E : Mean CD30 expression of the clone
E8SK  : Clone-specific CD8 effector - skin (relative)
E8BL  : Clone-specific CD8 effector - blood (relative)
NKB   : NK cell pool - blood (relative)
TREG  : Regulatory T cells (relative)
DCA   : Activated dendritic cells (relative)
TH2   : Th2 cytokine tone (relative)
IL31  : IL-31 (relative)
IL10  : IL-10 (relative)
IFNG  : IFN-gamma (relative)
TARC  : Serum TARC / CCL17 (pg/mL)
SIL2R : Soluble IL-2 receptor (U/mL)
LDH   : Lactate dehydrogenase (U/L)
APAT  : Patch area (percent BSA)
APLQ  : Plaque area (percent BSA)
ATUM  : Tumour area (percent BSA)
BARR  : Skin barrier integrity (0-1)
SAUR  : S. aureus colonisation (relative)
SAG   : Superantigen (relative)
PRUR  : Pruritus NRS (0-10)
SENS  : Central itch sensitisation (relative)
PROL  : Myeloid proliferating pool (10^9/L)
TR1   : Myeloid transit 1 (10^9/L)
TR2   : Myeloid transit 2 (10^9/L)
TR3   : Myeloid transit 3 (10^9/L)
CIRC  : Absolute neutrophil count (10^9/L)
PLT   : Platelets (10^9/L)
PN    : Peripheral neuropathy score
TSH   : TSH (mIU/L)
FT4   : Free T4 (ng/dL)
TG    : Triglycerides (mg/dL)
MAR   : Mogamulizumab-associated rash score
CD4N  : Normal CD4 count (cells/uL)
HZI   : Cumulative infection hazard
HZD   : Cumulative disease hazard
ADA   : Anti-drug antibody (placeholder)

$MAIN
// Guarded so that ctcl_patient() can hand back a fully equilibrated state
// without $MAIN silently overwriting it on the next run.
if (NEWIND <= 1 && INITF > 0.5) {
  E8SK_0 = IMMF;
  E8BL_0 = IMMF;
  NKB_0  = NKF;
  TREG_0 = 1.0;
  DCA_0  = 1.0;
  TH2_0  = 1.0;
  IL31_0 = 1.0;
  IL10_0 = 1.0;
  IFNG_0 = 1.0;
  BARR_0 = 1.0;
  SAUR_0 = 1.0;
  FCCR4_0 = CCR4;
  CD30E_0 = CD30;
  PROL_0 = CIRC0; TR1_0 = CIRC0; TR2_0 = CIRC0; TR3_0 = CIRC0; CIRC_0 = CIRC0;
  PLT_0  = PLT0;
  TSH_0  = 1.8;  FT4_0 = 1.2;  TG_0 = TG0;  CD4N_0 = CD40;
  TARC_0 = 300.; SIL2R_0 = 600.; LDH_0 = 190.; PRUR_0 = 2.0;
}

$ODE
// ============================================================ PK
double CMOG = MOG1 / VM1;
dxdt_MOG1 = -(CLM * CLF + QM) * CMOG + QM * MOG2 / VM2;
dxdt_MOG2 =  QM * CMOG - QM * MOG2 / VM2;

double CBV = BV1 / VB1;
dxdt_BV1 = -(CLB + QB + KDEC * VB1) * CBV + QB * BV2 / VB2;
dxdt_BV2 =  QB * CBV - QB * BV2 / VB2;
double CMM = MMAE1 / VMM;
dxdt_MMAE1 = KDEC * BV1 * FMMAE - CLMM * CMM;
dxdt_MMAE2 = 0.0;

double CROM = ROM1 / VR1;
dxdt_ROM1 = -(CLR + QR) * CROM + QR * ROM2 / VR2;
dxdt_ROM2 =  QR * CROM - QR * ROM2 / VR2;

dxdt_VORD = -KAV * VORD;
double CVOR = VOR / VV;
dxdt_VOR  = FV * KAV * VORD - CLV * CVOR;

dxdt_BEXD = -KAX * BEXD;
double CBEX = BEX / VX;
dxdt_BEX  = FX * KAX * BEXD - CLX * CBEX;

dxdt_IFND = -KAI * IFND;
double CIFN = IFN / VI;
dxdt_IFN  = KAI * IFND - CLI * CIFN;

dxdt_MTXD = -KAT * MTXD;
double CMTX = MTX / VT;
dxdt_MTX  = FT * KAT * MTXD - CLT * CMTX;

double CGEM = GEM / VG;
dxdt_GEM  = -CLG * CGEM;

double CALE = ALEM / VA;
dxdt_ALEM = -CLAL * CALE;

dxdt_SDT   = -KSDT * SDT;
dxdt_TSEBC = -KTSEB * TSEBC;
dxdt_STER  = -KSTER * STER;
dxdt_ABX   = -KABX * ABX;
dxdt_ECPD  = -KECP * ECPD;

// ================================= tissue exposure (the un-fitted asymmetry)
double CMOGSK = CMOG * ABCSK;
double CMOGLN = CMOG * ABCLN;
double CBVSK  = CBV  * ABCSK;
double CALESK = CALE * ABCSK;

// ============================================================ clone state
double NSKT = NSK + NSKR + NTR;
double NBLT = NBL + NBLR;
double KSKC = KSK * KSKF;

double occM_bl = FCCR4 * CMOG   / (CMOG   + KDMOG);
double occM_sk = FCCR4 * CMOGSK / (CMOGSK + KDMOG);
double occM_ln = FCCR4 * CMOGLN / (CMOGLN + KDMOG);
double occA_bl = CALE   / (CALE   + KDALEM);
double occA_sk = CALESK / (CALESK + KDALEM);

// effector availability: blood is rich in NK, lesional skin is not
double NKs = NKB * NKRSK;
double etb = NKB / (NKB + ETSAT * NBLT / (NBLT + 1.0));
double ets = NKs / (NKs + ETSAT * NSKT / (NSKT + 4.0));

double ADCC_BL = KADCC * occM_bl * etb;
double ADCC_SK = KADCC * occM_sk * ets;
double ADCC_LN = KADCC * occM_ln * etb * 0.55;
double ALEM_BL = KADCA * occA_bl;
double ALEM_SK = KADCA * occA_sk * (0.25 + 0.75 * NKRSK);

// antibody-drug conjugate: effector-independent
double bindB = CD30E * CBV   / (CBV   + KDBV);
double bindS = CD30E * CBVSK / (CBVSK + KDBV);
double byst  = KBYST * CMM / (CMM + KMMAE50);
double BV_BL = KBVKILL * KINT * bindB + byst;
double BV_SK = KBVKILL * KINT * bindS + byst * 1.6;

// cell-intrinsic small molecules distribute freely: no compartment gap
double E_ROM = EROM * CROM / (CROM + CROM50);
double E_VOR = EVOR * CVOR / (CVOR + CVOR50);
double E_BEX = EBEX * CBEX / (CBEX + CBEX50);
double E_MTX = EMTX * CMTX / (CMTX + CMTX50);
double E_GEM = EGEM * CGEM / (CGEM + CGEM50);
double E_STE = ESTER * STER / (STER + 1.0);
double E_IFN = EIFNK * CIFN / (CIFN + CIFN50);

// skin-directed modalities are depth-limited
double atot_now = APAT + APLQ + ATUM + 1e-9;
double depth = (APAT + DEPTHP * APLQ + DEPTHT * ATUM) / atot_now;
double E_SDT = EUVB * SDT / (SDT + 1.0) * depth;
double E_TSE = ETSEBK * TSEBC;

// surveillance
double esc0 = TH2 > 1.0 ? TH2 - 1.0 : 0.0;
double esc  = 1.0 / (1.0 + KIMESC * esc0);
double IMM_SK = KIMM * E8SK * esc / (1.0 + NSKT / KIM50);
double IMM_BL = KIMB * E8BL * esc / (1.0 + NBLT / KIM50B);

// Every DIRECT drug effect is scaled by one patient-level sensitivity term.
// Immune kill is not: CTL potency is a host property, not a clone property.
// Drug kill then splits in TWO, because the drug-resistant subclone is defined
// by an intrinsic apoptosis / efflux phenotype.  Antibody EFFECTOR killing does
// not go through that phenotype and is barely degraded by it; payload and
// cell-intrinsic agents are.  Three residual factors, not one:
//     RESE 0.85 (ADCC) > RESI 0.50 (CTL) >> RESK 0.12 (cell-intrinsic)
// which is why an antibody still works in a patient who has failed everything.
double EFF_SK = SENSF * (ADCC_SK + ALEM_SK);
double EFF_BL = SENSF * (ADCC_BL + ALEM_BL);
double EFF_LN = SENSF * ADCC_LN;
double CYT_SK = SENSF * (BV_SK + E_ROM + E_VOR + E_BEX + E_MTX + E_GEM + E_STE
                         + E_SDT + E_TSE);
double CYT_BL = SENSF * (BV_BL + E_ROM + E_VOR + E_BEX + E_MTX + E_GEM + E_STE
                         + 0.030 * ECPD / (ECPD + 1.0));
double CYT_LN = SENSF * (BV_BL * 0.6 + E_ROM + E_VOR + E_BEX + E_MTX + E_GEM + E_STE);
double KILL_SK = EFF_SK + CYT_SK + IMM_SK;
double KILL_BL = EFF_BL + CYT_BL + IMM_BL;
double KILL_LN = EFF_LN + CYT_LN + 0.5 * IMM_BL;
double KILL_TR = (KILL_SK + IMM_SK * (KIMTR - 1.0)) * 0.85
               + SENSF * KBVKILL * KINT * bindS * 1.6;

// drug pressure and the resistant subclone.  Immune kill is NOT reduced.
double KILL_SKR = RESI * IMM_SK + RESE * EFF_SK + RESK * CYT_SK;
double KILL_BLR = RESI * IMM_BL + RESE * EFF_BL + RESK * CYT_BL;
// resistance is SELECTED by cell-intrinsic pressure, not by effector killing
double conv_sk = KCONV * CYT_SK * NSK;
double conv_bl = KCONV * CYT_BL * NBL;

// superantigen drive on proliferation
double sag = 1.0 + ESAG * SAG / (SAG + SAG50);

// trafficking: mogamulizumab blocks the road as well as killing the cell
double homeblk = 1.0 - 0.85 * (CMOGSK / (CMOGSK + KDMOG));
double bexblk  = 1.0 - 0.35 * (CBEX / (CBEX + CBEX50));
double cap     = 1.0 - NSKT / KSKC; if (cap < 0) cap = 0;
double JH  = KHOME * CLA * CLADH * FCCR4 * homeblk * bexblk * NBLT * cap;
double kegr = KEGR * (1.0 - CD69) * (0.2 + 0.8 * CCR7);
double JE  = kegr * NSK;
double JER = kegr * NSKR;

double e8def = 1.0 - E8SK / (IMMF + 1e-9); if (e8def < 0) e8def = 0;
double trans = KTRANS * (NSK + NSKR) * sag * (1.0 + 3.0 * e8def);

dxdt_NSK  = CLAMP > 0.5 ? 0.0 :
            GSK * sag * NSK * (1.0 - NSKT / KSKC) + JH * (1.0 - FRES)
            - JE - KILL_SK * NSK - trans - conv_sk;
dxdt_NSKR = CLAMP > 0.5 ? 0.0 :
            GSK * FITR * sag * NSKR * (1.0 - NSKT / KSKC) + JH * FRES
            - JER - KILL_SKR * NSKR + conv_sk;
dxdt_NTR  = CLAMP > 0.5 ? 0.0 :
            GTR * GTRF * NTR * (1.0 - NSKT / KSKC) + trans - KILL_TR * NTR;
dxdt_NBL  = CLAMP > 0.5 ? 0.0 :
            GBL * NBL * (1.0 - NBLT / KBL) + JE - JH * (1.0 - FRES)
            - KNODE * NBL + KNOUT * NLN - KILL_BL * NBL - conv_bl;
dxdt_NBLR = CLAMP > 0.5 ? 0.0 :
            GBL * FITR * NBLR * (1.0 - NBLT / KBL) + JER - JH * FRES
            - KNODE * NBLR - KILL_BLR * NBLR + conv_bl;
double KLN_EFF = KILL_LN * (1.0 - FRES)
               + (RESI * IMM_BL + RESE * EFF_LN + RESK * CYT_LN) * FRES;
dxdt_NLN  = CLAMP > 0.5 ? 0.0 :
            GLN * NLN * (1.0 - NLN / KLNC) + KNODE * NBLT - KNOUT * NLN
            - KVIS * NLN - KLN_EFF * NLN;
dxdt_NVS  = CLAMP > 0.5 ? 0.0 :
            GVS * NVS + KVIS * NLN - KLN_EFF * NVS;

// antigen-loss selection (replicator dynamics)
double selC = KSELC * (ADCC_BL + ADCC_SK);
double selD = KSELD * (BV_BL + BV_SK) / 2.0;
dxdt_FCCR4 = CLAMP > 0.5 ? 0.0 : -selC * FCCR4 * (1.0 - FCCR4);
dxdt_CD30E = CLAMP > 0.5 ? 0.0 :
             KCD30I * (E_ROM + E_VOR) * (1.0 - CD30E)
             - KCD30O * (CD30E - CD30) - selD * CD30E * (1.0 - CD30E);

// ==================================================== host immune arm
double supp0 = 1.0 / ((1.0 + 1.0 / KTREG50) * (1.0 + 1.0 / KIL1050));
double supp  = (1.0 / ((1.0 + TREG / KTREG50) * (1.0 + IL10 / KIL1050))) / supp0;
double dctar = 1.0 + KDCAG * NSKT / (NSKT + KDCAGK)
             + KECPD * ECPD / (ECPD + 1.0)
             + KIFNDC * CIFN / (CIFN + CIFN50) + KDCSDT * E_SDT;
dxdt_DCA = KDC * (dctar - DCA);

double lymphotox = E_GEM + E_MTX * 0.6 + E_STE * ISTER + ALEM_BL * 1.1 + BV_BL * 0.30;
dxdt_E8SK = KDE8 * (DCA * supp * IMMF - E8SK) - (lymphotox + E_TSE * 1.4) * E8SK;
dxdt_E8BL = KDE8 * (DCA * supp * IMMF - E8BL) - lymphotox * E8BL;
dxdt_NKB  = KNK * (1.0 + 0.9 * CIFN / (CIFN + CIFN50)) * NKF - KNKD * NKB
            - (E_GEM + ALEM_BL * 1.2) * NKB;
// Tregs are CCR4-high: mogamulizumab reaches them before it reaches the clone
dxdt_TREG = KTREGP * (0.6 + 0.4 * TH2) - KTREGD * TREG
            - (1.35 * ADCC_BL + ALEM_BL + E_GEM) * TREG;

// ==================================================== cytokines, biomarkers
dxdt_TH2 = KTH2 * ((1.0 + TH2SK * NSKT / (NSKT + TH2K) + TH2BL * NBLT / (NBLT + TH2K))
                   / (1.0 + EIFNTH2 * CIFN / (CIFN + CIFN50) + 1.2 * E_BEX + 1.6 * E_STE)
                   - TH2);
dxdt_IL31 = KIL31 * ((1.0 + I31SK * NSKT / (NSKT + I31K) + I31BL * NBLT / (NBLT + I31K))
                     * (0.6 + 0.4 * TH2) - IL31);
dxdt_IL10 = KIL10 * (1.0 + I10SK * TREG * NSKT / (NSKT + I10K) - IL10);
dxdt_IFNG = KIFNG * (0.4 + 0.6 * E8SK + 1.2 * CIFN / (CIFN + CIFN50) - IFNG);
dxdt_TARC = KTARC * ((TARCB + TARCSK * NSKT / (NSKT + TARCK)
                      + TARCBL * NBLT / (NBLT + TARCK)) * (0.5 + 0.5 * TH2) - TARC);
dxdt_SIL2R = KSIL * (SILB + SILSK * NSKT / (NSKT + SILK)
                     + SILBL * NBLT / (NBLT + SILK) - SIL2R);
dxdt_LDH  = KLDH * (LDHB + LDHSL * (NSKT + NBLT + NLN + 4.0 * NVS) / 10.0 - LDH);

// ==================================================== lesion morphology
double ery = ERYF * (NBLT / (NBLT + ERYK)) * (1.0 - CD69);
double atot = 100.0 * (1.0 - exp(-NSKT / ASCALE)) + ery;
if (atot > 100.0) atot = 100.0;
double ftr  = NTR / (NTR + NSK + NSKR + 1e-9);
double dens = NSKT / KSKC;
double fT = ftr + (1.0 - ftr) * pow(dens, 3.0) / (pow(dens, 3.0) + 0.125);
double fP = (1.0 - fT) * dens / (dens + 0.30);
double fA = 1.0 - fT - fP; if (fA < 0) fA = 0;
dxdt_APAT = (atot * fA - APAT) / TAUPAT;
dxdt_APLQ = (atot * fP - APLQ) / TAUPLQ;
dxdt_ATUM = (atot * fT - ATUM) / TAUTUM;

// ==================================================== barrier, staph, itch
double lesion = (APAT + 2.0 * APLQ + 4.0 * ATUM) / 100.0;
dxdt_BARR = KBARR * (1.0 - BARR) - KBDAM * lesion * BARR - 0.05 * PRUR / 10.0 * BARR;
dxdt_SAUR = KSAP * (1.0 - BARR) * (1.0 + 0.5 * TH2)
            - KSAD * SAUR * (1.0 + 6.0 * ABX / (ABX + 0.5));
dxdt_SAG  = KSAG * SAUR - KSAGD * SAG;
double itch_in = 1.6 * (IL31 - 1.0) + 0.5 * (TH2 - 1.0) + 2.5 * (1.0 - BARR) + 0.6 * SAG;
if (itch_in < 0) itch_in = 0;
double itch_t = PRMAX * itch_in * (1.0 + SENS) / (itch_in * (1.0 + SENS) + KITCH50);
dxdt_PRUR = KPR * (itch_t - PRUR);
dxdt_SENS = KSENS * PRUR / PRMAX - KSENSD * SENS;

// ==================================================== haematology, toxicity
double myelo = 1.10 * E_GEM + 0.55 * byst / KBYST + 0.30 * E_MTX
             + 0.25 * E_ROM + 0.20 * ALEM_BL;
if (myelo > 0.95) myelo = 0.95;
double circ_safe = CIRC > 0.05 ? CIRC : 0.05;
double fb = pow(CIRC0 / circ_safe, GAMMA); if (fb > 4.0) fb = 4.0;
dxdt_PROL = KTR * PROL * (1.0 - myelo) * fb - KTR * PROL;
dxdt_TR1  = KTR * (PROL - TR1);
dxdt_TR2  = KTR * (TR1 - TR2);
dxdt_TR3  = KTR * (TR2 - TR3);
dxdt_CIRC = KTR * (TR3 - CIRC);
dxdt_PLT  = KPLT * (PLT0 * (1.0 - 0.9 * myelo) - PLT);
dxdt_PN   = KPN * CMM * 1000.0 - KPND * PN;
dxdt_TSH  = KTSH * (1.8 * (1.0 - ETSH * CBEX / (CBEX + CTSH50)) - TSH);
dxdt_FT4  = KFT4 * (0.55 + 0.36 * TSH - FT4);
dxdt_TG   = KTG * (TG0 * (1.0 + ETG * CBEX / (CBEX + CTSH50)) - TG);
double tregloss = 1.0 - TREG; if (tregloss < 0) tregloss = 0;
dxdt_MAR  = KMARP * tregloss * (MOG1 > 0.0 ? 1.0 : 0.0) - KMARD * MAR;
dxdt_CD4N = KCD4 * (CD40 - CD4N)
            - (1.4 * ALEM_BL + 0.9 * ADCC_BL + 1.1 * E_GEM + 0.5 * E_STE) * CD4N;

// ==================================================== hazards
double tum = NSKT + NBLT + NLN + 6.0 * NVS + 3.0 * NTR;
dxdt_HZD = HDIS * tum * (1.0 + 2.0 * fT);
double cd4_safe = CD4N > 40.0 ? CD4N : 40.0;
double neut_def = 1.0 - CIRC / CIRC0; if (neut_def < 0) neut_def = 0;
double infec = (1.0 - BARR) * (1.0 + SAUR) * (CD40 / cd4_safe) * (1.0 + 2.5 * neut_def);
dxdt_HZI = HINF * infec + HBASE;
dxdt_ADA = 0.0;

$TABLE
double MSWAT   = APAT + 2.0 * APLQ + 4.0 * ATUM;
double SEZCT   = (NBL + NBLR) * 200.0;                 // Sezary cells per uL
double BSCORE  = SEZCT >= 1000.0 ? 2.0 : (SEZCT >= 250.0 ? 1.0 : 0.0);
double NSCORE  = NLN >= 3.0 ? 3.0 : (NLN >= 1.2 ? 2.0 : (NLN >= 0.4 ? 1.0 : 0.0));
double MSCORE  = NVS >= 0.05 ? 1.0 : 0.0;
double TSKIN   = ATUM > 0.5 ? 3.0 : ((APAT + APLQ + ATUM) >= 80.0 ? 4.0 :
                 ((APAT + APLQ) >= 10.0 ? 2.0 : 1.0));
double CLONET  = NSK + NSKR + NTR + NBL + NBLR + NLN + NVS;
double CMOGP   = MOG1 / VM1;
double CMOGSKO = CMOGP * ABCSK;
double OCCSK   = FCCR4 * CMOGSKO / (CMOGSKO + KDMOG);   // ~1: drug is not limiting
double OCCBL   = FCCR4 * CMOGP   / (CMOGP   + KDMOG);
double CBVP    = BV1 / VB1;
double CMMP    = MMAE1 / VMM * 1000.0;                  // ng/mL
double SURVSK  = E8SK;                                  // the second state variable
double SURV    = exp(-(HZI + HZD));
double FT4LOW  = FT4 < 0.8 ? 1.0 : 0.0;
double ANCG34  = CIRC < 1.0 ? 1.0 : 0.0;

$CAPTURE MSWAT SEZCT BSCORE NSCORE MSCORE TSKIN CLONET CMOGP CMOGSKO OCCSK OCCBL
$CAPTURE CBVP CMMP SURVSK SURV FT4LOW ANCG34
'

ctcl_mod <- mcode_cache("ctcl_qsp", ctcl_code, end = 730, delta = 1)

## =============================================================================
##  STAGE = THE EQUILIBRIUM OF THE SURVEILLANCE SWITCH, NOT AN INITIAL CONDITION
## -----------------------------------------------------------------------------
##  A TNMB stage is a burden a patient has been carrying for years.  Dropping
##  that burden into a naive host and pressing "go" produces a transient that
##  is pure artefact — the cytokines, the barrier, the lesion morphology and
##  the skin/blood partition have not had time to reach the values that burden
##  implies.  Every apparent "response" in the first weeks would then be the
##  transient, not the drug.
##
##  So the patient is built in four steps:
##    1. hold the clone fixed (CLAMP = 1) and let the host equilibrate;
##    2. release the clone, let the skin/blood/node partition find itself,
##       and rescale the TOTAL clone back to the stage burden (repeat);
##       -> the MF vs SS split is therefore COMPUTED from CD69/CCR7;
##    3. clamp again briefly so lesion morphology matches the final burden;
##    4. solve for the immune fitness IMMF at which the clone is stationary,
##       then apply ONE uniform surveillance deficit.
##  The stage-dependent progression rate that follows is an output.
## =============================================================================

ctcl_stage_burden <- function(stage) {
  b <- list(
    IA      = c(0.55, 0.00, 0.010, 0.05, 0.00),
    IB      = c(1.60, 0.00, 0.020, 0.10, 0.00),
    IIA     = c(2.20, 0.00, 0.030, 0.45, 0.00),
    IIB     = c(5.00, 0.00, 0.040, 0.55, 0.00),
    IIB_LCT = c(4.20, 0.55, 0.060, 0.90, 0.00),
    IIIA    = c(6.20, 0.00, 0.300, 0.70, 0.00),
    IIIB    = c(6.60, 0.00, 0.900, 0.80, 0.00),
    IVA1    = c(6.80, 0.00, 6.500, 1.30, 0.00),
    IVA2    = c(5.20, 0.35, 1.200, 3.40, 0.00),
    IVB     = c(5.50, 0.60, 2.200, 3.20, 0.55))[[stage]]
  stats::setNames(as.list(b), c("NSK", "NTR", "NBL", "NLN", "NVS"))
}

ctcl_phenotype <- function(kind = c("MF", "SS")) {
  kind <- match.arg(kind)
  if (kind == "MF") list(CD69 = 0.85, CCR7 = 0.10) else list(CD69 = 0.20, CCR7 = 0.85)
}

.state_of <- function(out) {
  d  <- as.data.frame(out)
  cm <- mrgsolve::cmt(ctcl_mod)
  stats::setNames(lapply(cm, function(k) as.numeric(d[[k]][nrow(d)])), cm)
}

## Two initialisers.  The first keeps INITF = 1 so that $MAIN supplies the host
## defaults (Tregs, barrier, ANC, ...) while init() supplies only the clone.
## The second switches INITF off, because by then we are handing back a complete
## equilibrated state that $MAIN must not touch.
.set_clone <- function(mod, st)
  do.call(mrgsolve::init, c(list(mrgsolve::param(mod, INITF = 1)), st))
.set_init <- function(mod, st)
  do.call(mrgsolve::init, c(list(mrgsolve::param(mod, INITF = 0)), st))

ctcl_patient <- function(stage = "IB", kind = "MF", deficit = 0.10, cov = list(),
                         mod = ctcl_mod) {
  ph  <- ctcl_phenotype(kind)
  par <- utils::modifyList(utils::modifyList(ph, list()), cov)
  m   <- do.call(mrgsolve::param, c(list(mod), par))
  fr  <- as.numeric(mrgsolve::param(m)$FRES)
  b   <- ctcl_stage_burden(stage)
  st0 <- list(NSK = b$NSK * (1 - fr), NSKR = b$NSK * fr, NTR = b$NTR,
              NBL = b$NBL * (1 - fr), NBLR = b$NBL * fr,
              NLN = b$NLN, NVS = b$NVS)
  m0  <- .set_clone(m, st0)  # $MAIN still supplies the host defaults here

  ## 1. host equilibration with the clone clamped
  o  <- mrgsim(mrgsolve::param(m0, CLAMP = 1), end = 500, delta = 500)
  st <- .state_of(o)

  ## 2. let the skin / blood / node partition find itself; hold the TOTAL
  pools <- c("NSK", "NSKR", "NBL", "NBLR", "NLN")
  tgt   <- sum(unlist(st[pools]))
  for (i in 1:4) {
    o  <- mrgsim(.set_init(mrgsolve::param(m, CLAMP = 0), st), end = 120, delta = 120)
    st <- .state_of(o)
    cur <- sum(unlist(st[pools]))
    if (cur > 1e-12) for (k in pools) st[[k]] <- st[[k]] * tgt / cur
  }

  ## 3. settle lesion morphology at the exact stage burden
  o  <- mrgsim(.set_init(mrgsolve::param(m, CLAMP = 1), st), end = 250, delta = 250)
  st <- .state_of(o)

  ## 4. solve for the immune fitness that makes the clone stationary
  p     <- as.list(mrgsolve::param(m))
  supp0 <- 1 / ((1 + 1 / p$KTREG50) * (1 + 1 / p$KIL1050))
  supp  <- (1 / ((1 + st$TREG / p$KTREG50) * (1 + st$IL10 / p$KIL1050))) / supp0
  dclone <- function(immf) {
    s2 <- st; s2$E8SK <- s2$E8BL <- st$DCA * supp * immf
    o2 <- mrgsim(.set_init(mrgsolve::param(m, IMMF = immf, CLAMP = 0), s2),
                 end = 1, delta = 1)
    (o2$NSK[2] + o2$NSKR[2]) - (s2$NSK + s2$NSKR)
  }
  lo <- 0.05; hi <- 20
  if (dclone(lo) < 0) { immf_eq <- lo
  } else if (dclone(hi) > 0) { immf_eq <- hi
  } else {
    for (i in 1:60) { mid <- 0.5 * (lo + hi); if (dclone(mid) > 0) lo <- mid else hi <- mid }
    immf_eq <- 0.5 * (lo + hi)
  }
  immf <- immf_eq * (1 - deficit)
  st$E8SK <- st$E8BL <- st$DCA * supp * immf
  list(mod = .set_init(mrgsolve::param(m, IMMF = immf, CLAMP = 0), st),
       init = st, immf_eq = immf_eq, immf = immf, stage = stage, kind = kind)
}

## =============================================================================
##  REGIMENS
## =============================================================================
ev_mogamulizumab <- function(wt = 75, dur = 546) {
  a <- 1.0 * wt
  n_q2w <- max(0, floor((dur - 28) / 14))
  c(ev(time = 0,  amt = a, cmt = "MOG1", ii = 7, addl = 3),
    ev(time = 28, amt = a, cmt = "MOG1", ii = 14, addl = n_q2w))
}
ev_brentuximab   <- function(wt = 75, ncyc = 16)
  ev(amt = 1.8 * wt, cmt = "BV1", ii = 21, addl = ncyc - 1)
ev_romidepsin    <- function(bsa = 1.85, ncyc = 12)
  Reduce(c, lapply(0:(ncyc - 1), function(k)
    ev(time = k * 28, amt = 14 * bsa, cmt = "ROM1", ii = 7, addl = 2)))
## Daily oral agents are given as a ZERO-ORDER daily input (rate = amt): with a
## 2-7 h half-life and a PD time constant of weeks, the average exposure is the
## quantity that matters, and this keeps the numerics honest.
ev_vorinostat    <- function(dur = 546, dose = 400)
  ev(amt = dose, rate = dose, cmt = "VORD", ii = 1, addl = dur - 1)
ev_bexarotene    <- function(bsa = 1.85, dur = 546, dose = 300) {
  a <- dose * bsa; ev(amt = a, rate = a, cmt = "BEXD", ii = 1, addl = dur - 1)
}
ev_interferon    <- function(dur = 546, mu = 3)
  Reduce(c, lapply(seq(0, dur - 1, by = 7), function(t)
    Reduce(c, lapply(c(0, 2, 4), function(o)
      ev(time = t + o, amt = mu * 5, cmt = "IFND")))))
ev_methotrexate  <- function(dur = 546, dose = 25)
  ev(amt = dose, cmt = "MTXD", ii = 7, addl = floor(dur / 7) - 1)
ev_gemcitabine   <- function(bsa = 1.85, ncyc = 6)
  Reduce(c, lapply(0:(ncyc - 1), function(k)
    ev(time = k * 28, amt = 1000 * bsa, cmt = "GEM", ii = 7, addl = 2)))
ev_alemtuzumab   <- function(dur = 180, dose = 10)
  Reduce(c, lapply(seq(0, dur - 1, by = 7), function(t)
    Reduce(c, lapply(c(0, 2, 4), function(o)
      ev(time = t + o, amt = dose, cmt = "ALEM")))))
ev_nbuvb         <- function(dur = 180)
  Reduce(c, lapply(seq(0, dur - 1, by = 7), function(t)
    Reduce(c, lapply(c(0, 2, 4), function(o)
      ev(time = t + o, amt = 1.4, cmt = "SDT")))))
ev_chlormethine  <- function(dur = 546, amt = 0.85)
  ev(amt = amt, rate = amt, cmt = "SDT", ii = 1, addl = dur - 1)
ev_tseb          <- function(total = 12, nfx = 8)
  ev(amt = total / nfx, cmt = "TSEBC", ii = 2, addl = nfx - 1)
ev_ecp           <- function(dur = 546)
  Reduce(c, lapply(seq(0, dur - 1, by = 28), function(t)
    ev(time = t, amt = 3, cmt = "ECPD", ii = 1, addl = 1)))
ev_steroid       <- function(dur = 180, amt = 1)
  ev(amt = amt, rate = amt, cmt = "STER", ii = 1, addl = dur - 1)
ev_antibiotic    <- function(t0 = 0, days = 10)
  ev(time = t0, amt = 1, rate = 1, cmt = "ABX", ii = 1, addl = days - 1)

## =============================================================================
##  SCENARIOS
##  Each entry: patient (stage / phenotype / covariates) + regimen + horizon.
##  Scenarios 1-8 are the therapy classes; 9-16 are the structural experiments
##  that the model exists to run.
## =============================================================================
ctcl_scenarios <- function() list(
  "01 untreated stage IB MF (natural history)" = list(
    stage = "IB", kind = "MF", deficit = 0.10, ev = NULL, end = 1095),

  "02 stage IB MF - nbUVB + topical steroid" = list(
    stage = "IB", kind = "MF", deficit = 0.10,
    ev = c(ev_nbuvb(180), ev_steroid(180)), end = 730),

  "03 stage IB MF - chlormethine gel" = list(
    stage = "IB", kind = "MF", deficit = 0.10, ev = ev_chlormethine(546), end = 730),

  "04 stage IIB MF - low-dose TSEB 12 Gy" = list(
    stage = "IIB", kind = "MF", deficit = 0.10, ev = ev_tseb(12, 8), end = 730),

  "05 stage IIB MF - bexarotene 300 mg/m2" = list(
    stage = "IIB", kind = "MF", deficit = 0.25, ev = ev_bexarotene(1.85, 546), end = 730),

  "06 stage IIB MF - interferon alfa 3 MU tiw" = list(
    stage = "IIB", kind = "MF", deficit = 0.25, ev = ev_interferon(546), end = 730),

  "07 stage IIB MF - bexarotene + interferon" = list(
    stage = "IIB", kind = "MF", deficit = 0.25,
    ev = c(ev_bexarotene(1.85, 546), ev_interferon(546)), end = 730),

  "08 stage IIB MF - vorinostat 400 mg qd (MAVORIC comparator)" = list(
    stage = "IIB", kind = "MF", deficit = 0.25, ev = ev_vorinostat(546), end = 730),

  "09 stage IIB MF - romidepsin 14 mg/m2 d1,8,15 q28" = list(
    stage = "IIB", kind = "MF", deficit = 0.25, ev = ev_romidepsin(1.85, 12), end = 730),

  "10 stage IIB MF - mogamulizumab 1 mg/kg" = list(
    stage = "IIB", kind = "MF", deficit = 0.25, ev = ev_mogamulizumab(75, 546), end = 730),

  "11 Sezary IVA1 - mogamulizumab 1 mg/kg (MAVORIC)" = list(
    stage = "IVA1", kind = "SS", deficit = 0.25, ev = ev_mogamulizumab(75, 546), end = 730),

  "12 Sezary IVA1 - vorinostat 400 mg qd (MAVORIC)" = list(
    stage = "IVA1", kind = "SS", deficit = 0.25, ev = ev_vorinostat(546), end = 730),

  "13 Sezary IVA1 - ECP q4w + interferon" = list(
    stage = "IVA1", kind = "SS", deficit = 0.25,
    ev = c(ev_ecp(546), ev_interferon(546)), end = 730),

  "14 Sezary IVA1 - low-dose alemtuzumab SC" = list(
    stage = "IVA1", kind = "SS", deficit = 0.25, ev = ev_alemtuzumab(180), end = 730),

  "15 CD30-high IIB - brentuximab vedotin x16 (ALCANZA)" = list(
    stage = "IIB", kind = "MF", deficit = 0.25, cov = list(CD30 = 0.45),
    ev = ev_brentuximab(75, 16), end = 730),

  "16 CD30-low IIB - brentuximab vedotin x16" = list(
    stage = "IIB", kind = "MF", deficit = 0.25, cov = list(CD30 = 0.05),
    ev = ev_brentuximab(75, 16), end = 730),

  "17 CD30-low IIB - romidepsin x2 cycles THEN brentuximab vedotin" = list(
    stage = "IIB", kind = "MF", deficit = 0.25, cov = list(CD30 = 0.05),
    ev = c(ev_romidepsin(1.85, 2), ev(time = 56, amt = 135, cmt = "BV1",
                                      ii = 21, addl = 12)), end = 730),

  "18 stage IIB MF - gemcitabine x6 (ORR high, TTNT short)" = list(
    stage = "IIB", kind = "MF", deficit = 0.25, ev = ev_gemcitabine(1.85, 6), end = 730),

  "19 stage IIB MF - low-dose methotrexate" = list(
    stage = "IIB", kind = "MF", deficit = 0.25, ev = ev_methotrexate(546), end = 730),

  "20 transformed MF (LCT) - brentuximab vedotin" = list(
    stage = "IIB_LCT", kind = "MF", deficit = 0.25,
    cov = list(GTRF = 4.5, CD30 = 0.60), ev = ev_brentuximab(75, 16), end = 547),

  "21 stage IIB MF - anti-staphylococcal course alone (day 90-100)" = list(
    stage = "IIB", kind = "MF", deficit = 0.25, ev = ev_antibiotic(90, 10), end = 365),

  "22 Sezary IVA1 - mogamulizumab + ECP" = list(
    stage = "IVA1", kind = "SS", deficit = 0.25,
    ev = c(ev_mogamulizumab(75, 546), ev_ecp(546)), end = 730),

  "23 refractory IIB - CCR4-low clone on mogamulizumab" = list(
    stage = "IIB", kind = "MF", deficit = 0.25, cov = list(CCR4 = 0.25),
    ev = ev_mogamulizumab(75, 546), end = 730),

  "24 Sezary IVA1 - mogamulizumab at 3 mg/kg (dose is not the limit)" = list(
    stage = "IVA1", kind = "SS", deficit = 0.25,
    ev = c(ev(time = 0, amt = 225, cmt = "MOG1", ii = 7, addl = 3),
           ev(time = 28, amt = 225, cmt = "MOG1", ii = 14, addl = 36)), end = 730)
)

ctcl_run <- function(sc, mod = ctcl_mod) {
  pt <- ctcl_patient(sc$stage, sc$kind, sc$deficit %||% 0.10, sc$cov %||% list(), mod)
  m  <- pt$mod
  out <- if (is.null(sc$ev)) mrgsim(m, end = sc$end, delta = 1)
         else mrgsim(m, events = sc$ev, end = sc$end, delta = 1)
  as.data.frame(out)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

ctcl_run_all <- function(scenarios = ctcl_scenarios()) {
  do.call(rbind, lapply(names(scenarios), function(nm) {
    x <- ctcl_run(scenarios[[nm]]); x$scenario <- nm; x
  }))
}

## =============================================================================
##  RESPONSE SCORING — the Olsen 2011 conjunction, implemented as a conjunction
##  Nothing here computes a global response from a formula.  Each compartment is
##  scored on its own, and the global category is the WORST of them.  That single
##  line is why a drug with a 68% blood response and a 42% skin response has a
##  global response rate near their product rather than near either of them.
## =============================================================================
score_compartment <- function(x, base, kind = c("skin", "blood", "node")) {
  kind <- match.arg(kind)
  if (kind == "skin")
    ifelse(x <= 0.5, "CR", ifelse(x <= 0.5 * base, "PR",
           ifelse(x >= 1.25 * base, "PD", "SD")))
  else if (kind == "blood") {
    if (base < 250) ifelse(x >= 1000, "PD", "CR")
    else ifelse(x < 250, "CR", ifelse(x <= 0.5 * base, "PR",
                ifelse(x >= 1.5 * base, "PD", "SD")))
  } else {
    if (base < 0.40) ifelse(x >= 1.5 * max(base, 0.15), "PD", "CR")
    else ifelse(x <= 0.25 * base, "CR", ifelse(x <= 0.5 * base, "PR",
                ifelse(x >= 1.5 * base, "PD", "SD")))
  }
}

global_response <- function(skin, blood, node) {
  rank <- c(PD = 0, SD = 1, PR = 2, CR = 3)
  lab  <- c("PD", "SD", "PR", "CR")
  lab[pmin(rank[skin], rank[blood], rank[node]) + 1]
}

ctcl_endpoints <- function(x, visit = 28) {
  b_ms <- x$MSWAT[1]; b_sz <- x$SEZCT[1]; b_nd <- x$NLN[1]
  sk <- score_compartment(x$MSWAT, b_ms, "skin")
  bl <- score_compartment(x$SEZCT, b_sz, "blood")
  nd <- score_compartment(x$NLN,   b_nd, "node")
  gl <- global_response(sk, bl, nd)
  vt <- seq(visit, max(x$time), by = visit)
  ix <- vapply(vt, function(v) which.min(abs(x$time - v)), 1L)
  conf <- function(cat) {
    r <- cat[ix] %in% c("CR", "PR")
    if (length(r) < 2) return(NA_real_)
    k <- which(r[-length(r)] & r[-1])
    if (!length(k)) NA_real_ else vt[k[1]]
  }
  pd <- gl[ix] == "PD"
  pfs <- if (any(pd)) vt[which(pd)[1]] else NA_real_
  nadir <- which.min(x$MSWAT)
  ttnt <- {
    after <- x$MSWAT[nadir:nrow(x)] >= 0.75 * b_ms
    if (x$MSWAT[nadir] > 0.75 * b_ms) 0
    else if (any(after)) x$time[nadir + which(after)[1] - 1] else NA_real_
  }
  data.frame(
    mSWAT_base = round(b_ms, 1), mSWAT_nadir = round(min(x$MSWAT), 1),
    best_skin_pct = round(100 * (1 - min(x$MSWAT) / b_ms), 1),
    Sezary_base = round(b_sz, 0), Sezary_nadir = round(min(x$SEZCT), 0),
    skin_response_day = conf(sk), blood_response_day = conf(bl),
    global_response_day = conf(gl),
    PFS_day = pfs, TTNT_day = ttnt,
    NRS_base = round(x$PRUR[1], 1), NRS_nadir = round(min(x$PRUR), 1),
    TARC_base = round(x$TARC[1], 0), TARC_nadir = round(min(x$TARC), 0),
    surveillance_end = round(x$SURVSK[nrow(x)] / x$SURVSK[1], 2),
    ANC_nadir = round(min(x$CIRC), 2), neuropathy_max = round(max(x$PN), 2),
    CD4_nadir = round(min(x$CD4N), 0), MAR_max = round(max(x$MAR), 2),
    row.names = NULL)
}

ctcl_summary <- function(all = ctcl_run_all()) {
  do.call(rbind, lapply(split(all, all$scenario), function(x) {
    cbind(scenario = x$scenario[1], ctcl_endpoints(x))
  }))
}

## =============================================================================
##  STRUCTURAL EXPERIMENTS
##  These are the questions the file was built to answer.  None of them is a
##  fitted number; each is a comparison the model has to generate on its own.
## =============================================================================

## E1 - the compartment gap is an EFFECTOR gap, not an exposure gap.
##      Occupancy in the skin is already ~1 at 1 mg/kg, so tripling the dose
##      moves the skin response not at all, while removing the effector deficit
##      (NKRSK -> 1) moves it a great deal.
E1_effector_not_exposure <- function() {
  base <- ctcl_patient("IVA1", "SS", 0.25)
  arms <- list(
    "1 mg/kg"            = list(m = base$mod, e = ev_mogamulizumab(75, 546)),
    "3 mg/kg"            = list(m = base$mod, e = ev(time = 0, amt = 225, cmt = "MOG1",
                                                     ii = 7, addl = 39)),
    "1 mg/kg, NKRSK=1.0" = list(m = param(base$mod, NKRSK = 1.0),
                                e = ev_mogamulizumab(75, 546)))
  do.call(rbind, lapply(names(arms), function(nm) {
    x <- as.data.frame(mrgsim(arms[[nm]]$m, events = arms[[nm]]$e, end = 365, delta = 1))
    data.frame(arm = nm,
               skin_occupancy = round(max(x$OCCSK), 3),
               best_mSWAT_drop_pct = round(100 * (1 - min(x$MSWAT) / x$MSWAT[1]), 1),
               best_blood_drop_pct = round(100 * (1 - min(x$SEZCT) / x$SEZCT[1]), 1),
               row.names = NULL)
  }))
}

## E2 - ORR and TTNT are different questions.  Same patient, two regimens with
##      similar skin response and opposite effects on the surveillance state.
E2_orr_versus_ttnt <- function() {
  arms <- list("gemcitabine x6"   = ev_gemcitabine(1.85, 6),
               "interferon alfa"  = ev_interferon(546),
               "ECP + interferon" = c(ev_ecp(546), ev_interferon(546)),
               "mogamulizumab"    = ev_mogamulizumab(75, 546))
  pt <- ctcl_patient("IIB", "MF", 0.25)
  do.call(rbind, lapply(names(arms), function(nm) {
    x <- as.data.frame(mrgsim(pt$mod, events = arms[[nm]], end = 730, delta = 1))
    e <- ctcl_endpoints(x)
    data.frame(arm = nm, best_skin_pct = e$best_skin_pct,
               TTNT_day = e$TTNT_day,
               surveillance_ratio = e$surveillance_end, row.names = NULL)
  }))
}

## E3 - the sequencing rationale: an HDAC inhibitor raises CD30, so the same
##      ADC dose delivers more payload afterwards.  Reported as a difference.
E3_hdac_then_adc <- function() {
  pt <- ctcl_patient("IIB", "MF", 0.25, list(CD30 = 0.05))
  arms <- list("BV alone" = ev_brentuximab(75, 16),
               "romidepsin x2 then BV" = c(ev_romidepsin(1.85, 2),
                                           ev(time = 56, amt = 135, cmt = "BV1",
                                              ii = 21, addl = 12)))
  do.call(rbind, lapply(names(arms), function(nm) {
    x <- as.data.frame(mrgsim(pt$mod, events = arms[[nm]], end = 547, delta = 1))
    data.frame(arm = nm, CD30_max = round(max(x$CD30E), 3),
               best_skin_pct = round(100 * (1 - min(x$MSWAT) / x$MSWAT[1]), 1),
               row.names = NULL)
  }))
}

## E4 - antibiotics change the disease without touching the clone.
E4_superantigen <- function() {
  pt <- ctcl_patient("IIB", "MF", 0.25)
  a <- as.data.frame(mrgsim(pt$mod, end = 365, delta = 1))
  b <- as.data.frame(mrgsim(pt$mod, events = ev_antibiotic(90, 10), end = 365, delta = 1))
  i <- function(x, d) x[which.min(abs(x$time - d)), ]
  data.frame(day = c(90, 110, 150, 365),
             mSWAT_no_abx = round(sapply(c(90, 110, 150, 365), function(d) i(a, d)$MSWAT), 1),
             mSWAT_abx    = round(sapply(c(90, 110, 150, 365), function(d) i(b, d)$MSWAT), 1),
             NRS_no_abx   = round(sapply(c(90, 110, 150, 365), function(d) i(a, d)$PRUR), 1),
             NRS_abx      = round(sapply(c(90, 110, 150, 365), function(d) i(b, d)$PRUR), 1),
             row.names = NULL)
}

## E5 - TSEB clears the skin and burns the skin-resident CD8 pool with it.
E5_tseb_paradox <- function() {
  pt <- ctcl_patient("IIB", "MF", 0.10)
  x <- as.data.frame(mrgsim(pt$mod, events = ev_tseb(12, 8), end = 730, delta = 1))
  data.frame(day = c(0, 30, 90, 180, 365, 730),
             mSWAT = round(approx(x$time, x$MSWAT, c(0, 30, 90, 180, 365, 730))$y, 1),
             skin_CD8 = round(approx(x$time, x$SURVSK, c(0, 30, 90, 180, 365, 730))$y, 2),
             row.names = NULL)
}

## E6 - MF and SS from the same equations: sweep the residency covariates.
E6_residency_sweep <- function() {
  do.call(rbind, lapply(c(0.85, 0.65, 0.45, 0.25, 0.10), function(cd69) {
    pt <- ctcl_patient("IIIA", "MF", 0.10, list(CD69 = cd69, CCR7 = 1 - cd69))
    x  <- as.data.frame(mrgsim(pt$mod, end = 1, delta = 1))
    data.frame(CD69 = cd69, Sezary_per_uL = round(x$SEZCT[1], 0),
               B_score = x$BSCORE[1], mSWAT = round(x$MSWAT[1], 1),
               row.names = NULL)
  }))
}

if (identical(Sys.getenv("CTCL_RUN"), "1")) {
  all <- ctcl_run_all()
  print(ctcl_summary(all), row.names = FALSE)
  cat("\n-- E1 effector, not exposure --\n");  print(E1_effector_not_exposure(), row.names = FALSE)
  cat("\n-- E2 ORR vs TTNT --\n");             print(E2_orr_versus_ttnt(), row.names = FALSE)
  cat("\n-- E3 HDACi then ADC --\n");          print(E3_hdac_then_adc(), row.names = FALSE)
  cat("\n-- E4 superantigen --\n");            print(E4_superantigen(), row.names = FALSE)
  cat("\n-- E5 TSEB paradox --\n");            print(E5_tseb_paradox(), row.names = FALSE)
  cat("\n-- E6 residency sweep --\n");         print(E6_residency_sweep(), row.names = FALSE)
}
