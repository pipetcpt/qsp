## =============================================================================
##  nb_mrgsolve_model.R
##  High-Risk Neuroblastoma — Quantitative Systems Pharmacology model
##  46 ODEs · anti-GD2 antibody (two Fc arms) · isotretinoin · cytotoxics + ASCT
##             · 131I-MIBG dosimetry · ALK inhibition · IL-2/GM-CSF · endpoints
##
##  WHAT THIS MODEL IS ORGANISED AROUND
##  ----------------------------------------------------------------------------
##  1.  GD2 is expressed at ~8e6 copies/cell, i.e. 13.3 nmol of antigen per gram
##      of tumour.  A whole course of dinutuximab is 373 nmol.  Antibody is NOT
##      the scarce resource and fractional site occupancy is NOT the driver:
##      capillary permeability is.  The PD driver in this model is therefore
##      BOUND IgG PER CELL (molecules/cell), not occupancy.
##
##  2.  The antibody reaches THREE target compartments with three different
##      permeability-surface products per gram:
##          solid tumour   PSG_TU  0.020 L/(d*kg)  + an interstitial-pressure
##                                                   penalty that grows with size
##          bone marrow    PSG_BM  0.300 L/(d*kg)  sinusoidal, no barrier
##          DRG / nerve    PSG_NRV 4.000 L/(d*kg)  fenestrated, NO blood-nerve
##                                                   barrier at the ganglion
##      Everything clinically distinctive about anti-GD2 therapy follows from
##      that 200-fold spread: marrow disease is saturated at a tenth of the
##      clinical dose, solid tumour is still linear at ten times it, and the
##      nociceptor is the best-perfused target of the three.
##
##  3.  The two Fc effector functions have DIFFERENT HILL COEFFICIENTS in bound
##      IgG per cell, because they have different stoichiometries:
##          ADCC : n = 1   one Fc engages one FcgammaR       (ADCC50 = 2e4/cell)
##          CDC  : n = 2   C1q must BRIDGE A PAIR of Fc      (CDC50  = 5e4/cell)
##      The therapeutic arm is linear in surface density; the dose-limiting toxic
##      arm is quadratic.  This is why the K322A Fc mutation — which removes C1q
##      binding while leaving FcgammaR binding intact — is the only lever in the
##      model that moves the therapeutic index without paying efficacy for it.
##
##  4.  The four modalities differ in PROLIFERATION-DEPENDENCE, and that alone
##      forces the observed treatment order:
##          cytotoxics    kill TP strongly, TQ weakly (FQ_CT), TD almost not at
##                        all (FD_CT) -> antagonised by differentiation therapy
##          ADCC          proliferation-INDEPENDENT -> compatible with retinoids
##          MIBG beta-    proliferation-INDEPENDENT -> hits the quiescent pool
##          retinoids     move cells OUT of the cytotoxic-sensitive pool
##
##  UNITS
##    time          days
##    tumour        1e9 cells   (1e9 cells ~ 1 g ~ 1 cm3)
##    antibody      nmol (amounts) / nM (concentrations)
##    isotretinoin  umol / uM
##    cytotoxic     mg (amount) / mg/L
##    MIBG          MBq (activity), umol (molecular mass), Gy (absorbed dose)
##    ANC, PLT      1e9/L
##
##  CALIBRATION TARGETS (see nb_references.md for sources)
##    dinutuximab   Cmax ~11.5 ug/mL, terminal t1/2 ~10 d at 17.5 mg/m2/d x 4 d
##    isotretinoin  peak ~2-4 uM at 160 mg/m2/d divided BID; exposure FALLS
##                  across the 14-day course via CYP26A1 autoinduction
##    131I-MIBG     whole-body absorbed dose ~0.20 mGy/MBq; 2 Gy whole body is
##                  the stem-cell-rescue threshold; tumour dose median ~15-30 Gy
##    ANBL0032      2-year EFS 66% (immunotherapy + isotretinoin) vs 46%
##    HR-NBL1       adding IL-2 to dinutuximab beta: no EFS gain, more toxicity
##    ADVL0912      crizotinib active in only 1/11 ALK-mutant neuroblastoma
##
##  Requires: mrgsolve (>= 1.0), dplyr, tidyr, ggplot2
## =============================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(tidyr)
})

nb_code <- '
$PROB High-Risk Neuroblastoma QSP — two Fc arms, three permeabilities

$PARAM @annotated
// ------------------------------------------------------------------ patient
BSA    : 0.80   : Body surface area (m2)
WT     : 20.0   : Body weight (kg)

// ------------------------------------------------------------ tumour biology
KPROL      : 0.0300 : Proliferating-pool net growth rate (1/d)
KPROL_MYCN : 1.70   : KPROL multiplier when MYCN amplified (-)
MYCN       : 1      : MYCN amplification flag (1 = amplified)
KDEATH     : 0.0100 : Spontaneous tumour cell loss (1/d)
KQ         : 0.0300 : TP -> TQ entry into quiescence (1/d)
KRE        : 0.0200 : TQ -> TP re-entry into cycle (1/d)
KRE_MYCN   : 1.80   : KRE multiplier when MYCN amplified (-)
KDIS       : 0.0040 : TP -> TM marrow dissemination (1/d)
KSEED      : 0.0030 : TM -> TP re-seeding (1/d)
KDDEATH    : 0.0200 : Differentiated-cell loss (1/d)
KPROL_M    : 0.85   : Marrow-pool growth relative to primary (-)
KVTU       : 0.20   : Tumour-volume relaxation toward cell number (1/d)
FINT       : 0.35   : Interstitial volume fraction of tissue (-)
NEXT       : 1e-9   : Extinction floor = one cell (1e9 cells)

// ---------------------------------------------------------------------- GD2
GD2DENS     : 8.0e6  : GD2 surface density (molecules/cell)
GD2_MES     : 0.15   : GD2 density multiplier in the MES cell state (-)
FMES        : 0.05   : MES-state fraction of the tumour (-)
FGD2_RA     : 1.00   : Retinoid effect on GD2 density (1 = neutral; sign disputed)
GD2DENS_NRV : 0.5e6  : GD2 density on nerve/DRG (molecules/cell-equivalent)
NRVCELLS    : 2.0    : GD2-positive nerve tissue (1e9 cell-equivalents)
MNRV        : 0.0020 : Mass of GD2-positive DRG/nerve tissue (kg)
MBM         : 0.300  : Red-marrow mass (kg)

// -------------------------------------------------- anti-GD2 antibody PK/PD
MWAB    : 150000 : Antibody molecular weight (g/mol)
V1AB    : 3.30   : Central volume (L)
V2AB    : 1.90   : Peripheral volume (L)
CLAB    : 0.41   : Clearance (L/d)
QAB     : 0.40   : Intercompartmental clearance (L/d)
KON     : 5.00   : Association rate constant (1/(nM*d)) ~ 6e4 /M/s
KDAB    : 20.0   : GD2 dissociation constant (nM)
NVAL    : 2.0    : GD2 sites occupied per bound antibody (-)
KINT    : 0.030  : Internalisation of bound complex (1/d) — GD2 is a glycolipid
PSG_TU  : 0.020  : Permeability-surface per kg, solid tumour (L/(d*kg))
PSG_BM  : 0.300  : Permeability-surface per kg, bone marrow (L/(d*kg))
PSG_NRV : 4.000  : Permeability-surface per kg, DRG/nerve (L/(d*kg))
VIFP    : 0.033  : Tumour volume at which tumour PS halves (L) — IFP/ECM barrier
C1QEFF  : 1.00   : C1q-binding efficiency of the Fc (1.00 ch14.18; 0.10 K322A)
FCGRX   : 1.00   : FcgammaR-binding multiplier of the Fc variant (-)

// ---------------------------------------------------- Fc ARM A : ADCC (n=1)
KADCC    : 0.90   : Maximum ADCC kill rate (1/d)
ADCC50   : 2.0e4  : Bound IgG per cell for half-maximal ADCC (molecules/cell)
KE       : 180.0  : Effector half-saturation (cells/uL equivalents)
WG       : 0.045  : Granulocyte weight relative to NK in ADCC (-)
FCG      : 1.00   : FcgammaR genotype multiplier (1.00 F/F, 1.35 V/F, 1.80 V/V)
FTM_ADCC : 1.25   : ADCC efficiency in marrow relative to solid tumour (-)

// ----------------------------------------------------- Fc ARM B : CDC (n=2)
KCDC   : 1.00   : Complement activation flux scale (-)
CDC50  : 5.0e4  : Bound IgG per cell for half-maximal CDC (molecules/cell)
HCDC   : 2.0    : CDC Hill coefficient — C1q bridges a PAIR of Fc domains (-)
KPAIN  : 9.00   : Pain gain from complement flux (-)
KOFFP  : 1.20   : Pain resolution rate (1/d)
KSYNC  : 0.90   : Complement pool resynthesis (1/d)
KCONSC : 1.80   : Complement consumption per unit flux (1/d)

// ------------------------------------------------- NK / Treg / exhaustion
NKB0        : 200.0 : Baseline blood NK cells (cells/uL)
KOUTNK      : 0.10  : NK turnover (1/d)
KTRAF       : 0.020 : NK trafficking out of blood (1/d)
FINF        : 0.10  : Effector infiltration efficiency into tumour (-)
KDEGNKT     : 0.15  : Intratumoural NK loss (1/d)
EMAX_IL2_NK : 3.50  : IL-2 Emax on NK expansion (-)
EC50_IL2_NK : 40.0  : IL-2 EC50 for NK expansion (IU-equivalent)
TREG0       : 40.0  : Baseline Treg (cells/uL)
KOUTTR      : 0.08  : Treg turnover (1/d)
EMAX_IL2_TR : 4.00  : IL-2 Emax on Treg expansion (-)
EC50_IL2_TR : 5.00  : IL-2 EC50 for Treg — CD25 is the HIGH-affinity receptor
TREG50      : 110.0 : Treg for half-maximal ADCC suppression (cells/uL)
KEXH        : 0.30  : NK exhaustion induction (1/d)
KREC        : 0.10  : NK exhaustion recovery (1/d)
EC50_EXH    : 30.0  : IL-2 EC50 for NK exhaustion (IU-equivalent)
EMAX_GM_NK  : 0.60  : GM-CSF Emax on NK/CD16 (-)
EC50_GM_NK  : 3.00  : GM-CSF EC50 on NK (ug/L)

// ------------------------------------------- myelosuppression (Friberg 3-tr)
CIRC0        : 4.00  : Baseline ANC (1e9/L)
MTT          : 4.583 : Mean transit time (d) = 110 h
GAM          : 0.16  : Feedback exponent (-)
SLOPE_CT     : 0.550 : Cytotoxic slope on progenitors (L/mg)
SLOPE_RAD    : 1.60  : Whole-body dose-rate slope on progenitors (1/(Gy/d))
EMAX_GM_ANC  : 2.20  : GM-CSF Emax on ANC setpoint (-)
EC50_GM_ANC  : 3.00  : GM-CSF EC50 on ANC (ug/L)
PLT0         : 250.0 : Baseline platelets (1e9/L)
KPLT         : 0.42  : Platelet relaxation rate (1/d)
SLOPE_PLT    : 0.300 : Cytotoxic slope on platelets (L/mg)

// -------------------------------------------------------------- isotretinoin
MWRA     : 300.4 : Isotretinoin molecular weight (g/mol)
KARA     : 1.80  : Absorption rate constant (1/d)
VRA      : 70.0  : Apparent volume V/F (L)
CLRA     : 73.0  : Apparent clearance CL/F (L/d)
EMAXI    : 1.50  : CYP26A1 autoinduction Emax (CL x 2.5 at plateau) (-)
EC50I    : 1.20  : Autoinduction EC50 (uM)
KENZ     : 0.35  : CYP26A1 turnover (1/d)
KDIFFMAX : 0.080 : Maximum retinoid-driven differentiation rate (1/d)
EC50RA   : 1.50  : Differentiation EC50 (uM)

// ---------------------------------------------------- composite cytotoxic PK
V1CT   : 12.0  : Central volume (L)
V2CT   : 20.0  : Peripheral volume (L)
CLCT   : 30.0  : Clearance (L/d)
QCT    : 8.00  : Intercompartmental clearance (L/d)
EMAXCT : 1.95  : Maximum cytotoxic kill rate on TP (1/d)
EC50CT : 1.10  : Cytotoxic EC50 (mg/L)
FQ_CT  : 0.12  : Relative cytotoxic kill of the QUIESCENT pool (-)
FD_CT  : 0.02  : Relative cytotoxic kill of the DIFFERENTIATED pool (-)
FM_CT  : 0.55  : Relative cytotoxic kill in the marrow sanctuary (-)
FPLAT  : 0.55  : Platinum fraction of the composite exposure (-)
KOTO   : 0.0016: Ototoxicity accumulation (per mg/L/d)
ESTS   : 0.00  : Sodium-thiosulfate cochlear protection (0-0.75)
FTUMSTS: 0.00  : Sodium-thiosulfate TUMOUR protection (0-0.30) — the trade-off

// ------------------------------------------------------- external-beam RT
RTON   : 1     : External-beam radiotherapy flag
TRT    : 170.0 : RT start (d)
TRTDUR : 16.0  : RT duration (d)
RTGY   : 21.6  : RT total dose (Gy)
ARAD_EB: 0.177 : LQ effective kill per Gy, 1.8 Gy fractions (alpha 0.15, a/b 10)

// ------------------------------------------------------------------ surgery
SURGON : 1     : Surgical resection flag
TSURG  : 98.0  : Resection time (d)
FRESID : 0.10  : Residual fraction after resection (-)

// ---------------------------------------------------------------- 131I-MIBG
LAMPHYS  : 0.0864 : 131I physical decay constant (1/d) = t1/2 8.02 d
KEXCR    : 0.190  : Whole-body biological clearance (1/d)
CLMIBG   : 9.00   : MIBG renal clearance (L/d)
VMIBG    : 14.0   : MIBG distribution volume (L)
VMAXNET  : 1.05   : Max NET-mediated uptake (umol/(d*1e9 cells))
KMNET    : 0.30   : NET Km for MIBG (uM)
NETX     : 1.00   : NET expression / blockade multiplier (labetalol etc. < 1)
KWASH    : 0.28   : Tumour washout, set by VMAT vesicular retention (1/d)
SWB      : 5.54e-5: Whole-body S factor (Gy/(MBq*d)) -> 0.20 mGy/MBq
STU      : 2.63e-3: Tumour S factor (Gy*kg/(MBq*d)), 131I mean 0.19 MeV/decay
ALPHARAD : 0.060  : Effective LQ alpha at MIBG dose rates (1/Gy)
FTHY     : 0.020  : Fraction of blood activity trapped by thyroid (1/d)
KIBLOCK  : 1.00   : Potassium-iodide thyroid blockade (0.10 = blocked)

// ------------------------------------------------------------ ALK inhibitor
KAALK   : 6.00  : Absorption rate constant (1/d)
VALK    : 60.0  : Volume (L)
CLALK   : 45.0  : Clearance (L/d)
FUALK   : 0.09  : FREE fraction (crizotinib 0.09; lorlatinib 0.34)
IC50ALK : 0.30  : FREE-drug IC50 (uM) — R1275Q 0.30, F1174L 2.50, lorla 0.05
EMAXALK : 0.55  : Maximum fractional reduction of KPROL (-)

// ---------------------------------------------------------------- cytokines
KAIL2 : 8.00 : IL-2 absorption (1/d)
VIL2  : 4.00 : IL-2 volume (L)
CLIL2 : 60.0 : IL-2 clearance (L/d)
KAGM  : 6.00 : GM-CSF absorption (1/d)
VGM   : 4.50 : GM-CSF volume (L)
CLGM  : 40.0 : GM-CSF clearance (L/d)

// --------------------------------------------------------------- biomarkers
KSYNH : 1.20 : Urinary HVA production per 1e9 tumour cells (1/d)
KELH  : 2.20 : HVA elimination (1/d)

// ------------------------------------------------------------------ switches
IMMUNO : 1 : Anti-GD2 immunotherapy active (-)
USE_RA : 1 : Isotretinoin differentiation effect active (-)

$CMT @annotated
TP    : Proliferating tumour cells, solid disease (1e9 cells)
TQ    : Quiescent (G0) tumour cells (1e9 cells)
TD    : Retinoid-differentiated post-mitotic cells (1e9 cells)
TM    : Marrow minimal-residual-disease pool (1e9 cells)
VTU   : Tumour volume (L)
ABC   : Antibody, central (nmol)
ABP   : Antibody, peripheral (nmol)
ABT   : Antibody free in tumour interstitium (nmol)
BNDT  : Antibody bound to tumour GD2 (nmol)
ABM   : Antibody free in marrow interstitium (nmol)
BNDM  : Antibody bound to marrow GD2 (nmol)
ABN   : Antibody free in DRG/nerve interstitium (nmol)
BNDN  : Antibody bound to nerve GD2 (nmol)
NKB   : Blood NK cells (cells/uL)
NKT   : Tumour-infiltrating NK cells (cells/uL equivalent)
TREG  : Regulatory T cells (cells/uL)
NKEXH : NK exhaustion fraction (-)
CPL   : Functional complement pool (fraction of normal)
PROL  : Marrow proliferative progenitors (1e9/L)
TR1   : Transit 1 (1e9/L)
TR2   : Transit 2 (1e9/L)
TR3   : Transit 3 (1e9/L)
ANC   : Circulating absolute neutrophil count (1e9/L)
PLT   : Platelets (1e9/L)
RAG   : Isotretinoin gut depot (umol)
RAA   : Isotretinoin central amount (umol)
FIND  : CYP26A1 induction state (-)
CTC   : Composite cytotoxic, central (mg)
CTP   : Composite cytotoxic, peripheral (mg)
AWB   : 131I whole-body retained activity (MBq)
MBC   : 131I-MIBG blood activity (MBq)
MASS  : Total MIBG molecular mass, hot + cold (umol)
MBT   : 131I-MIBG tumour-retained activity (MBq)
DWB   : Cumulative whole-body absorbed dose (Gy)
DTU   : Cumulative tumour absorbed dose (Gy)
THY   : Thyroid retained activity (MBq)
AKG   : ALK inhibitor gut depot (umol)
AKC   : ALK inhibitor central (umol)
IL2D  : IL-2 subcutaneous depot (IU-equivalent)
IL2C  : IL-2 central (IU-equivalent)
GMD   : GM-CSF subcutaneous depot (ug)
GMC   : GM-CSF central (ug)
PAIN  : Allodynia / neuropathic pain intensity (-)
HVA   : Urinary homovanillic acid (relative)
OTO   : Cumulative ototoxicity (relative)
PLATAUC : Cumulative platinum AUC (mg*d/L)

$GLOBAL
#define NA_NMOL 6.022e14    // molecules per nmol

$MAIN
CPL_0   = 1.0;
NKB_0   = NKB0;
TREG_0  = TREG0;
PROL_0  = CIRC0;
TR1_0   = CIRC0;
TR2_0   = CIRC0;
TR3_0   = CIRC0;
ANC_0   = CIRC0;
PLT_0   = PLT0;

$ODE
// ------------------------------------------------------------ concentrations
double CP   = ABC / V1AB;
double CPP  = ABP / V2AB;
double NSOL  = TP + TQ + TD;
double NSOLe = (NSOL > NEXT) ? NSOL : NEXT;     // extinction floor = 1 cell
double TMe   = (TM   > NEXT) ? TM   : NEXT;
double VTUf  = (VTU > NEXT/1000.0) ? VTU : NEXT/1000.0;
double VTUi  = FINT * VTUf;
double VBMi  = FINT * MBM;
double VNRi  = FINT * MNRV;
double CT = ABT / VTUi;
double CM = ABM / VBMi;
double CN = ABN / VNRi;
double CRA  = RAA / VRA;
double CCT  = CTC / V1CT;
double CALK = FUALK * AKC / VALK;
double CIL2 = IL2C / VIL2;
double CGM  = GMC  / VGM;
double CMIB = MASS / VMIBG;

// ------------------------------------------------------- GD2 site pools (nmol)
double dens = GD2DENS * ((1.0 - FMES) + FMES * GD2_MES)
              * (1.0 + (FGD2_RA - 1.0) * CRA / (CRA + EC50RA));
double AGT = NSOLe * 1e9 * dens / NA_NMOL;
double AGM = TMe   * 1e9 * dens / NA_NMOL;
double AGN = NRVCELLS * 1e9 * GD2DENS_NRV / NA_NMOL;

// -------------------------------- binding: TMDD in each target compartment
double KOFF = KON * KDAB;
double fsT = AGT - NVAL * BNDT;  if (fsT < 0.0) fsT = 0.0;
double fsM = AGM - NVAL * BNDM;  if (fsM < 0.0) fsM = 0.0;
double fsN = AGN - NVAL * BNDN;  if (fsN < 0.0) fsN = 0.0;
double bT = KON * CT * fsT / NVAL - KOFF * BNDT;
double bM = KON * CM * fsM / NVAL - KOFF * BNDM;
double bN = KON * CN * fsN / NVAL - KOFF * BNDN;

// ---- permeability-surface products. THE 200-FOLD SPREAD IS THE WHOLE MODEL.
double PS_T = PSG_TU  * (NSOLe / 1000.0) / (1.0 + VTUf / VIFP);
double PS_M = PSG_BM  * MBM;
double PS_N = PSG_NRV * MNRV;
double fT = PS_T * (CP - CT);
double fM = PS_M * (CP - CM);
double fN = PS_N * (CP - CN);

dxdt_ABC  = - CLAB * CP - QAB * (CP - CPP) - fT - fM - fN;
dxdt_ABP  =   QAB * (CP - CPP);
dxdt_ABT  =   fT - bT;
dxdt_BNDT =   bT - KINT * BNDT;
dxdt_ABM  =   fM - bM;
dxdt_BNDM =   bM - KINT * BNDM;
dxdt_ABN  =   fN - bN;
dxdt_BNDN =   bN - KINT * BNDN;

// ------------------ BOUND IgG PER CELL — the PD driver, not site occupancy
double BPCT = BNDT * NA_NMOL / (NSOLe * 1e9);
double BPCM = BNDM * NA_NMOL / (TMe   * 1e9);
double BPCN = BNDN * NA_NMOL / (NRVCELLS * 1e9);

// ------------------------------------------------------------------ cytokines
dxdt_IL2D = - KAIL2 * IL2D;
dxdt_IL2C =   KAIL2 * IL2D - CLIL2 * CIL2;
dxdt_GMD  = - KAGM * GMD;
dxdt_GMC  =   KAGM * GMD - CLGM * CGM;

// -------------------------------------------------------------- effector cells
double eIL2NK = EMAX_IL2_NK * CIL2 / (CIL2 + EC50_IL2_NK);
double eGMNK  = EMAX_GM_NK  * CGM  / (CGM  + EC50_GM_NK);
double eIL2TR = EMAX_IL2_TR * CIL2 / (CIL2 + EC50_IL2_TR);
dxdt_NKB   = KOUTNK * NKB0 * (1.0 + eIL2NK + eGMNK) - (KOUTNK + KTRAF) * NKB;
dxdt_NKT   = KTRAF * NKB * FINF - KDEGNKT * NKT;
dxdt_TREG  = KOUTTR * TREG0 * (1.0 + eIL2TR) - KOUTTR * TREG;
dxdt_NKEXH = KEXH * CIL2 / (CIL2 + EC50_EXH) * (1.0 - NKEXH) - KREC * NKEXH;

// ---------------------------- Fc ARM A : ADCC, Hill n = 1 in bound IgG/cell
double EFF  = NKT + WG * ANC * 1000.0 * FINF;
double FSUP = 1.0 / (1.0 + TREG / TREG50);
double comm = KADCC * EFF / (EFF + KE) * FCG * FCGRX * (1.0 - NKEXH) * FSUP * IMMUNO;
double ADCC  = comm * BPCT / (BPCT + ADCC50);
double ADCCM = comm * FTM_ADCC * BPCM / (BPCM + ADCC50);

// ------------------- Fc ARM B : CDC, Hill n = 2 (C1q bridges a PAIR of Fc)
double rN  = pow(BPCN / CDC50, HCDC);
double CDCFLUX = KCDC * rN / (1.0 + rN) * CPL * C1QEFF * IMMUNO;
dxdt_CPL  = KSYNC * (1.0 - CPL) - KCONSC * CDCFLUX;
dxdt_PAIN = KPAIN * CDCFLUX - KOFFP * PAIN;

// -------------------------------------------------------- cytotoxic PK + kill
dxdt_CTC = - CLCT * CCT - QCT * (CCT - CTP / V2CT);
dxdt_CTP =   QCT * (CCT - CTP / V2CT);
double KILLCT = EMAXCT * CCT / (CCT + EC50CT) * (1.0 - FTUMSTS);
dxdt_PLATAUC = FPLAT * CCT;
dxdt_OTO     = KOTO * FPLAT * CCT * (1.0 - ESTS);

// --------------------------------------- isotretinoin + CYP26A1 autoinduction
dxdt_RAG  = - KARA * RAG;
dxdt_RAA  =   KARA * RAG - CLRA * (1.0 + FIND) * CRA;
dxdt_FIND =   KENZ * (EMAXI * CRA / (CRA + EC50I) - FIND);
double KDIFF = KDIFFMAX * CRA / (CRA + EC50RA) * USE_RA;

// ----------------------------------------------------------------- 131I-MIBG
dxdt_AWB = - (KEXCR + LAMPHYS) * AWB;
double upt_mass = VMAXNET * NETX * NSOLe * CMIB / (KMNET + CMIB);
// specific activity of the circulating pool falls with time: 131I decays, the
// carrier molecule does not. This is why no-carrier-added MIBG delivers more
// tumour dose per MBq than carrier-added MIBG.
double SA_now = MBC / ((MASS > 1e-9) ? MASS : 1e-9);
dxdt_MASS = - CLMIBG * CMIB - upt_mass;
dxdt_MBC  = - CLMIBG * MBC / VMIBG - LAMPHYS * MBC - upt_mass * SA_now;
dxdt_MBT  =   upt_mass * SA_now - (KWASH + LAMPHYS) * MBT;
double DR_WB = SWB * AWB;
double tum_kg = (NSOLe / 1000.0 > 1e-5) ? NSOLe / 1000.0 : 1e-5;
double DR_TU = STU * MBT / tum_kg;
dxdt_DWB = DR_WB;
dxdt_DTU = DR_TU;
dxdt_THY = FTHY * KIBLOCK * MBC - (LAMPHYS + 0.10) * THY;
double KILLMIBG = ALPHARAD * DR_TU;

// ------------------------------------------------------- external-beam RT
double RTRATE = (RTON > 0.5 && SOLVERTIME >= TRT && SOLVERTIME < TRT + TRTDUR)
                ? RTGY / TRTDUR : 0.0;
double KILLEB = ARAD_EB * RTRATE;

// ------------------------------------------------------------------- surgery
double RESECT = (SURGON > 0.5 && SOLVERTIME >= TSURG && SOLVERTIME < TSURG + 0.25)
                ? -log(FRESID) / 0.25 : 0.0;

// --------------------------------------------------------------- ALK inhibitor
dxdt_AKG = - KAALK * AKG;
dxdt_AKC =   KAALK * AKG - CLALK * AKC / VALK;
double ALKINH = EMAXALK * CALK / (CALK + IC50ALK);

// ----------------------------------------------------------- tumour dynamics
double kp  = KPROL * ((MYCN > 0.5) ? KPROL_MYCN : 1.0) * (1.0 - ALKINH);
double kre = KRE   * ((MYCN > 0.5) ? KRE_MYCN   : 1.0);
// proliferation-INDEPENDENT kill reaches every pool; cytotoxics do not
double KRADALL = KILLMIBG + KILLEB;

dxdt_TP = kp * TP - KDEATH * TP - KQ * TP + kre * TQ - KDIFF * TP
          - KDIS * TP + KSEED * TM
          - (KILLCT + ADCC + KRADALL + RESECT) * TP;
dxdt_TQ = KQ * TP - kre * TQ - KDEATH * TQ
          - (FQ_CT * KILLCT + ADCC + KRADALL + RESECT) * TQ;
dxdt_TD = KDIFF * TP - KDDEATH * TD
          - (FD_CT * KILLCT + ADCC + KRADALL + RESECT) * TD;
dxdt_TM = kp * KPROL_M * TM + KDIS * TP - KSEED * TM - KDEATH * TM
          - (FM_CT * KILLCT + ADCCM + KILLMIBG) * TM;
dxdt_VTU = KVTU * (NSOLe / 1000.0 - VTUf);

// ---- bound antibody is CARRIED BY CELLS: when a cell is destroyed its bound
// IgG leaves with it. Omitting this makes bound-IgG-per-cell diverge as burden
// falls and ADCC becomes self-reinforcing (a defect found by running the model).
double gT = (dxdt_TP + dxdt_TQ + dxdt_TD) / NSOLe;
double gM = dxdt_TM / TMe;
if (gT < 0.0) { dxdt_BNDT += gT * BNDT; dxdt_ABT += gT * ABT; }
if (gM < 0.0) { dxdt_BNDM += gM * BNDM; }

// --------------------------------------------------------- myelosuppression
double EDRUG = SLOPE_CT * CCT + SLOPE_RAD * DR_WB;
if (EDRUG > 0.98) EDRUG = 0.98;
double eGM  = EMAX_GM_ANC * CGM / (CGM + EC50_GM_ANC);
double KTR  = (4.0 / MTT) * (1.0 + 0.35 * eGM);
double ANCf = (ANC > 0.02) ? ANC : 0.02;
double fb   = pow(CIRC0 * (1.0 + eGM) / ANCf, GAM);
dxdt_PROL = KTR * PROL * ((1.0 - EDRUG) * fb - 1.0);
dxdt_TR1  = KTR * (PROL - TR1);
dxdt_TR2  = KTR * (TR1 - TR2);
dxdt_TR3  = KTR * (TR2 - TR3);
dxdt_ANC  = KTR * (TR3 - ANC);
double pdep = SLOPE_PLT * CCT;  if (pdep > 0.95) pdep = 0.95;
dxdt_PLT  = KPLT * (PLT0 * (1.0 - pdep) - PLT);

dxdt_HVA = KSYNH * (TP + TQ + TM) - KELH * HVA;

$TABLE
double BURDEN = TP + TQ + TD + TM;
double CPn    = ABC / V1AB;
double CPug   = CPn * MWAB / 1e6;           // ug/mL
double CRAo   = RAA / VRA;
double CCTo   = CTC / V1CT;
double NSOLo  = TP + TQ + TD;
double NSOLc  = (NSOLo > NEXT) ? NSOLo : NEXT;
double TMc    = (TM > NEXT) ? TM : NEXT;
double BPCTo  = BNDT * NA_NMOL / (NSOLc * 1e9);
double BPCMo  = BNDM * NA_NMOL / (TMc * 1e9);
double BPCNo  = BNDN * NA_NMOL / (NRVCELLS * 1e9);
double densO  = GD2DENS * ((1.0 - FMES) + FMES * GD2_MES);
double OCCT   = NVAL * BPCTo / densO;
double OCCN   = NVAL * BPCNo / GD2DENS_NRV;
double FREEALK= FUALK * AKC / VALK;
double SAo    = MBC / ((MASS > 1e-9) ? MASS : 1e-9);
double LOGB   = log10((BURDEN > 1e-12) ? BURDEN : 1e-12);
double MRD    = TM;

$CAPTURE BURDEN LOGB CPn CPug BPCTo BPCMo BPCNo OCCT OCCN CRAo CCTo FREEALK SAo MRD
'

nb_mod <- mcode("nb_qsp", nb_code, atol = 1e-10, rtol = 1e-8, maxsteps = 1e6)

## =============================================================================
##  SECTION 2 — dosing helpers
## =============================================================================

MW_AB <- 150000; MW_RA <- 300.4

## dinutuximab: 17.5 mg/m2/d over `hours`, x `days`
ev_ab <- function(t0, BSA = 0.80, mgm2 = 17.5, days = 4, hours = 10) {
  amt <- mgm2 * BSA / 1000 / MW_AB * 1e9                # nmol per day
  ev(time = t0, cmt = "ABC", amt = amt, rate = amt / (hours / 24),
     ii = 1, addl = days - 1)
}

## isotretinoin 160 mg/m2/d divided BID x 14 d
ev_ra <- function(t0, BSA = 0.80, mgm2 = 160, days = 14, FREL = 1.0) {
  amt <- mgm2 * BSA / 2 / MW_RA * 1000 * FREL           # umol per dose
  ev(time = t0, cmt = "RAG", amt = amt, ii = 0.5, addl = 2 * days - 1)
}

## GM-CSF 250 ug/m2/d SC
ev_gm <- function(t0, BSA = 0.80, days = 14, ugm2 = 250) {
  ev(time = t0, cmt = "GMD", amt = ugm2 * BSA, ii = 1, addl = days - 1)
}

## IL-2 4.5 MIU/m2/d as a 24-h infusion, two 4-day blocks
ev_il2 <- function(t0, BSA = 0.80, mium2 = 4.5, days = 4) {
  amt <- mium2 * BSA * 10
  ev(time = t0, cmt = "IL2D", amt = amt, rate = amt, ii = 1, addl = days - 1)
}

## composite cytotoxic cycle, specified as a target AUC (mg*d/L)
ev_chemo <- function(t0, auc = 6.0, days = 3, CLCT = 30) {
  amt <- auc * CLCT / days
  ev(time = t0, cmt = "CTC", amt = amt, rate = amt, ii = 1, addl = days - 1)
}

## autologous stem-cell rescue: reinfused graft repopulates the progenitor pool
ev_asct <- function(t0, cells = 3.0) ev(time = t0, cmt = "PROL", amt = cells)

## 131I-MIBG: activity and molecular mass must be dosed TOGETHER, because the
## ratio (specific activity) is what determines transporter competition.
ev_mibg <- function(t0, MBq, SA = 1.10e5) {
  bind_rows(
    ev(time = t0, cmt = "MBC",  amt = MBq,      rate = MBq / 0.08),
    ev(time = t0, cmt = "AWB",  amt = MBq,      rate = MBq / 0.08),
    ev(time = t0, cmt = "MASS", amt = MBq / SA, rate = MBq / SA / 0.08))
}

## ALK inhibitor, BID oral
ev_alk <- function(t0, t1, umol) {
  ev(time = t0, cmt = "AKG", amt = umol, ii = 0.5, addl = 2 * (t1 - t0) - 1)
}

## ---------------------------------------------------------------- initial state
nb_init <- function(mod, tp0 = 100, tm0 = 5) {
  init(mod, TP = tp0, TQ = 0.25 * tp0, TM = tm0, VTU = 1.25 * tp0 / 1000,
       HVA = 1.20 * (1.25 * tp0 + tm0) / 2.20)
}

## =============================================================================
##  SECTION 3 — the COG/SIOPEN high-risk protocol
## =============================================================================
IND_T   <- c(0, 21, 42, 63, 84, 105)
T_SURG  <- 98
T_HDCT  <- 133
T_ASCT  <- 140
T_RT    <- 170
IMMU_T  <- c(200, 228, 256, 284, 312, 340)

##' Build the full high-risk regimen.
##' @param immuno       give anti-GD2 antibody
##' @param use_il2      give IL-2 in cycles 2 and 4 (COG) vs omit it (SIOPEN)
##' @param use_gm       give GM-CSF in cycles 1, 3, 5
##' @param use_ra       give isotretinoin
##' @param hdct         "CEM" or "BuMel" (BuMel = deeper log-kill)
##' @param tandem       second consolidation + ASCT (ANBL0532)
##' @param ra_concurrent  give isotretinoin DURING induction (the antagonism test)
##' @param ab_mgm2      dinutuximab daily dose
##' @param ab_hours     infusion duration (naxitamab-style short infusion = 0.75)
nb_regimen <- function(immuno = TRUE, use_il2 = TRUE, use_gm = TRUE,
                       use_ra = TRUE, hdct = "CEM", tandem = FALSE,
                       ra_concurrent = FALSE, ab_mgm2 = 17.5, ab_hours = 10,
                       ab_cycles = c(1, 2, 3, 4, 5), BSA = 0.80, FREL = 1.0,
                       induction = TRUE, mibg = NULL, alk = NULL) {
  e <- NULL
  add <- function(x) if (is.null(e)) x else bind_rows(e, x)
  if (induction) {
    for (t0 in IND_T) e <- add(ev_chemo(t0, auc = 6.0, days = 3))
    if (ra_concurrent) for (t0 in IND_T) e <- add(ev_ra(t0 + 3, BSA, FREL = FREL))
    e <- add(ev_chemo(T_HDCT, auc = if (hdct == "BuMel") 38 else 32, days = 4))
    e <- add(ev_asct(T_ASCT))
    if (tandem) {
      e <- add(ev_chemo(T_HDCT + 42, auc = 22, days = 4))
      e <- add(ev_asct(T_ASCT + 42))
    }
  }
  if (!is.null(mibg)) e <- add(ev_mibg(mibg$t, mibg$MBq, mibg$SA %||% 1.10e5))
  if (!is.null(alk))  e <- add(ev_alk(alk$t0, alk$t1, alk$umol))
  for (i in seq_along(IMMU_T)) {
    t0 <- IMMU_T[i]
    if (use_ra && !ra_concurrent) e <- add(ev_ra(t0 + 11, BSA, FREL = FREL))
    if (!immuno || !(i %in% ab_cycles)) next
    if (i %in% c(1, 3, 5)) {
      if (use_gm) e <- add(ev_gm(t0, BSA))
      e <- add(ev_ab(t0 + 3, BSA, mgm2 = ab_mgm2, hours = ab_hours))
    } else {
      if (use_il2) {
        e <- add(ev_il2(t0, BSA)); e <- add(ev_il2(t0 + 7, BSA))
      }
      e <- add(ev_ab(t0 + 7, BSA, mgm2 = ab_mgm2, hours = ab_hours))
    }
  }
  e
}
`%||%` <- function(a, b) if (is.null(a)) b else a

##' Post-consolidation immunotherapy phase alone, starting from a chosen burden.
nb_immuno_only <- function(tstart = 0, ...) {
  IMMU_T <<- tstart + 28 * (0:5)
  on.exit(IMMU_T <<- c(200, 228, 256, 284, 312, 340))
  nb_regimen(induction = FALSE, ...)
}

nb_run <- function(mod = nb_mod, e, tp0 = 100, tm0 = 5, end = 1200, delta = 0.25,
                   param = list()) {
  m <- nb_init(mod, tp0, tm0)
  if (length(param)) m <- do.call(mrgsolve::param, c(list(m), param))
  m %>% mrgsim(events = e, end = end, delta = delta) %>% as_tibble()
}

## time at which burden re-crosses the detection threshold after nadir
nb_relapse <- function(out, from = 210, thr = 0.5) {
  d <- dplyr::filter(out, time >= from)
  if (!nrow(d)) return(Inf)
  i <- which.min(d$BURDEN)
  d2 <- d[i:nrow(d), ]
  j <- which(d2$BURDEN >= thr)
  if (!length(j)) Inf else d2$time[j[1]]
}

## =============================================================================
##  SECTION 4 — the nine scenarios
## =============================================================================

## ---- SCENARIO 1 -------------------------------------------------------------
## The two Fc arms against dose. ADCC is linear in bound IgG/cell (n=1), CDC is
## quadratic (n=2) but its compartment is already saturated at the clinical dose,
## so the dose-response geometry differs between the marrow and the solid tumour.
nb_scen_dose <- function(doses = c(0.175, 0.875, 1.75, 4.375, 8.75, 17.5, 35, 70, 175)) {
  purrr::map_dfr(doses, function(d) {
    o <- nb_run(e = nb_immuno_only(0, ab_mgm2 = d), tp0 = 0.002, tm0 = 0.001,
                end = 190, delta = 0.05)
    tibble(mgm2 = d, Cmax_nM = max(o$CPn), BPC_tumour = max(o$BPCTo),
           BPC_marrow = max(o$BPCMo), BPC_nerve = max(o$BPCNo),
           pain_AUC = sum(diff(o$time) * head(o$PAIN, -1)),
           pain_peak = max(o$PAIN),
           logkill = log10(o$BURDEN[1] / max(min(o$BURDEN), 1e-12)))
  })
}

## ---- SCENARIO 2 -------------------------------------------------------------
## Permeability, not dose. Sweeping tumour PS moves the tumour ADCC arm; sweeping
## dose barely does, because the interstitium equilibrates.
nb_scen_perm <- function(folds = c(0.1, 0.25, 0.5, 1, 2, 5, 10, 25)) {
  purrr::map_dfr(folds, function(f) {
    o <- nb_run(e = nb_immuno_only(0), tp0 = 0.002, tm0 = 0.001, end = 190,
                delta = 0.05, param = list(PSG_TU = 0.020 * f))
    tibble(PS_fold = f, BPC_tumour = max(o$BPCTo),
           logkill = log10(o$BURDEN[1] / max(min(o$BURDEN), 1e-12)))
  })
}

## ---- SCENARIO 3 -------------------------------------------------------------
## The K322A Fc mutation: C1q binding x0.1, FcgammaR binding untouched. Compare
## with the dose reduction of wild-type Fc that produces the same pain exposure.
nb_scen_k322a <- function() {
  f <- function(lbl, ...) {
    o <- nb_run(e = nb_immuno_only(0, ...), tp0 = 0.002, tm0 = 0.001,
                end = 190, delta = 0.05,
                param = list(...)[intersect(names(list(...)), "C1QEFF")])
    tibble(arm = lbl, pain_AUC = sum(diff(o$time) * head(o$PAIN, -1)),
           pain_peak = max(o$PAIN), CPL_min = min(o$CPL),
           logkill = log10(o$BURDEN[1] / max(min(o$BURDEN), 1e-12)))
  }
  bind_rows(
    {o <- nb_run(e = nb_immuno_only(0), tp0 = 0.002, tm0 = 0.001, end = 190,
                 delta = 0.05)
     tibble(arm = "ch14.18 (dinutuximab)",
            pain_AUC = sum(diff(o$time) * head(o$PAIN, -1)),
            pain_peak = max(o$PAIN), CPL_min = min(o$CPL),
            logkill = log10(o$BURDEN[1] / max(min(o$BURDEN), 1e-12)))},
    {o <- nb_run(e = nb_immuno_only(0), tp0 = 0.002, tm0 = 0.001, end = 190,
                 delta = 0.05, param = list(C1QEFF = 0.10))
     tibble(arm = "hu14.18K322A",
            pain_AUC = sum(diff(o$time) * head(o$PAIN, -1)),
            pain_peak = max(o$PAIN), CPL_min = min(o$CPL),
            logkill = log10(o$BURDEN[1] / max(min(o$BURDEN), 1e-12)))})
}

## ---- SCENARIO 4 -------------------------------------------------------------
## Burden dependence: the interstitial-pressure term is what makes anti-GD2 a
## minimal-residual-disease therapy rather than a debulking agent.
nb_scen_burden <- function(grams = c(0.001, 0.01, 0.1, 0.5, 1, 3, 7, 15, 30, 60, 125, 300)) {
  purrr::map_dfr(grams, function(g) {
    o <- nb_run(e = nb_immuno_only(0), tp0 = g, tm0 = g * 0.05, end = 100,
                delta = 0.1)
    tibble(burden_g = g, BPC_tumour = max(o$BPCTo),
           logkill = log10(o$BURDEN[1] / max(min(o$BURDEN), 1e-12)))
  })
}

## ---- SCENARIO 5 -------------------------------------------------------------
## IL-2 and GM-CSF. IL-2 expands NK cells but expands Treg at a LOWER EC50
## (CD25 is the high-affinity receptor) and induces exhaustion; GM-CSF raises
## the ANC, and the ANC IS the granulocyte ADCC effector.
nb_scen_cytokine <- function() {
  arms <- list("Ab alone" = list(use_il2 = FALSE, use_gm = FALSE),
               "Ab + GM-CSF" = list(use_il2 = FALSE, use_gm = TRUE),
               "Ab + IL-2" = list(use_il2 = TRUE, use_gm = FALSE),
               "Ab + GM-CSF + IL-2 (COG)" = list(use_il2 = TRUE, use_gm = TRUE))
  purrr::imap_dfr(arms, function(a, nm) {
    o <- nb_run(e = do.call(nb_immuno_only, c(list(0), a)), tp0 = 0.002,
                tm0 = 0.001, end = 190, delta = 0.05)
    tibble(arm = nm, NK_peak = max(o$NKB), Treg_peak = max(o$TREG),
           exhaustion = max(o$NKEXH), ANC_peak = max(o$ANC),
           pain_AUC = sum(diff(o$time) * head(o$PAIN, -1)),
           logkill = log10(o$BURDEN[1] / max(min(o$BURDEN), 1e-12)))
  })
}

## ---- SCENARIO 6 -------------------------------------------------------------
## Isotretinoin: CYP26A1 autoinduction destroys its own exposure across the
## 14-day course; opening the capsule into food costs ~40% of it; and giving it
## concurrently with cytotoxics is ANTAGONISTIC because differentiated cells are
## post-mitotic and therefore invisible to S/M-phase drugs.
nb_scen_retinoid <- function() {
  o1 <- nb_run(e = ev_ra(0), tp0 = 1, tm0 = 0.05, end = 20, delta = 0.02)
  o2 <- nb_run(e = ev_ra(0, FREL = 0.60), tp0 = 1, tm0 = 0.05, end = 20,
               delta = 0.02)
  peak_day <- function(o, d) max(o$CRAo[o$time >= d & o$time < d + 1])
  seq_arm  <- nb_run(e = nb_regimen(), end = 1200)
  conc_arm <- nb_run(e = nb_regimen(ra_concurrent = TRUE), end = 1200)
  list(exposure = tibble(
         day = c(1, 4, 8, 14),
         swallowed_uM = sapply(c(1, 4, 8, 14), function(d) peak_day(o1, d)),
         opened_uM    = sapply(c(1, 4, 8, 14), function(d) peak_day(o2, d))),
       schedule = tibble(
         arm = c("isotretinoin AFTER consolidation", "isotretinoin DURING induction"),
         post_induction_burden = c(seq_arm$BURDEN[which.min(abs(seq_arm$time - 132))],
                                  conc_arm$BURDEN[which.min(abs(conc_arm$time - 132))]),
         relapse_day = c(nb_relapse(seq_arm), nb_relapse(conc_arm))))
}

## ---- SCENARIO 7 -------------------------------------------------------------
## 131I-MIBG. Three things steal tumour dose: carrier MIBG competing for NET,
## a NET-blocking home medication, and a fixed mCi/kg prescription that ignores
## whole-body clearance.
nb_scen_mibg <- function(MBq = 666 * 20) {
  grid <- tibble(
    arm = c("no-carrier-added", "carrier-added (10x cold)",
            "NCA + labetalol not stopped", "NCA, half activity"),
    SA  = c(1.10e5, 1.10e4, 1.10e5, 1.10e5),
    NETX = c(1, 1, 0.35, 1),
    act = c(MBq, MBq, MBq, MBq / 2))
  purrr::pmap_dfr(grid, function(arm, SA, NETX, act) {
    o <- nb_run(e = ev_mibg(0, act, SA), tp0 = 20, tm0 = 1, end = 60,
                delta = 0.05, param = list(NETX = NETX))
    tibble(arm = arm, MBq = act, whole_body_Gy = max(o$DWB),
           tumour_Gy = max(o$DTU), mGy_per_MBq = 1000 * max(o$DWB) / act,
           ANC_nadir = min(o$ANC), thyroid_MBq = max(o$THY),
           logkill = log10(o$BURDEN[1] / max(min(o$BURDEN), 1e-12)))
  })
}

## ---- SCENARIO 8 -------------------------------------------------------------
## ALK inhibition is a FREE-FRACTION story. Crizotinib is 91% protein bound;
## lorlatinib ~66%. The same "total Cmax" means very different free exposure.
nb_scen_alk <- function() {
  grid <- tibble(
    drug = c("crizotinib / R1275Q", "crizotinib / F1174L",
             "lorlatinib / R1275Q", "lorlatinib / F1174L"),
    FUALK = c(0.09, 0.09, 0.34, 0.34),
    IC50 = c(0.30, 2.50, 0.05, 0.05),
    umol = c(1.10, 1.10, 0.42, 0.42))         # BID oral dose, umol
  purrr::pmap_dfr(grid, function(drug, FUALK, IC50, umol) {
    o <- nb_run(e = ev_alk(0, 28, umol), tp0 = 20, tm0 = 1, end = 40,
                delta = 0.05, param = list(FUALK = FUALK, IC50ALK = IC50))
    tibble(arm = drug, free_Cmax_uM = max(o$FREEALK),
           free_Ctrough_uM = min(o$FREEALK[o$time > 20 & o$time < 28]),
           IC50_uM = IC50,
           trough_over_IC50 = min(o$FREEALK[o$time > 20 & o$time < 28]) / IC50,
           logkill = log10(o$BURDEN[1] / max(min(o$BURDEN), 1e-12)))
  })
}

## ---- SCENARIO 9 -------------------------------------------------------------
## Sodium thiosulfate protects the cochlea. If it also protects the tumour, the
## trade-off is computable rather than rhetorical.
nb_scen_sts <- function(prot = c(0, 0.10, 0.20, 0.30)) {
  purrr::map_dfr(prot, function(fp) {
    o <- nb_run(e = nb_regimen(), end = 1200,
                param = list(ESTS = 0.75, FTUMSTS = fp))
    o0 <- nb_run(e = nb_regimen(), end = 1200)
    tibble(tumour_protection = fp, oto_STS = max(o$OTO), oto_noSTS = max(o0$OTO),
           relapse_STS = nb_relapse(o), relapse_noSTS = nb_relapse(o0))
  })
}

## =============================================================================
##  SECTION 5 — virtual cohort and EFS
## =============================================================================
nb_cohort <- function(n = 200, seed = 20260730) {
  set.seed(seed)
  lnv <- function(cv) exp(rnorm(n, 0, sqrt(log(1 + cv^2))))
  tibble(
    ID       = seq_len(n),
    KPROL    = 0.0300 * lnv(0.28),
    EMAXCT   = 1.95   * lnv(0.26),
    PSG_TU   = 0.020  * lnv(0.45),
    PSG_BM   = 0.300  * lnv(0.35),
    NKB0     = 200    * lnv(0.32),
    KADCC    = 0.90   * lnv(0.25),
    GD2DENS  = 8.0e6  * lnv(0.30),
    FMES     = pmin(rbeta(n, 1.4, 10) * 1.4, 0.60),
    FCG      = sample(c(1.00, 1.35, 1.80), n, TRUE, c(0.35, 0.45, 0.20)),
    MYCN     = rbinom(n, 1, 0.45),
    tp0      = 100 * lnv(0.45),
    tm0      = 5   * lnv(0.90))
}

nb_efs <- function(cohort, end = 1150, ...) {
  e <- nb_regimen(...)
  purrr::pmap_dbl(cohort, function(ID, tp0, tm0, ...) {
    pp <- list(...)
    o <- nb_run(e = e, tp0 = tp0, tm0 = tm0, end = end, delta = 1,
                param = pp)
    nb_relapse(o)
  })
}

nb_efs_rate <- function(times, day) mean(times > day)

## =============================================================================
##  SECTION 6 — example session
## =============================================================================
if (interactive()) {
  ## the headline contrast
  a <- nb_run(e = nb_regimen(immuno = FALSE), end = 1200)
  b <- nb_run(e = nb_regimen(immuno = TRUE),  end = 1200)
  cat(sprintf("isotretinoin only : relapse day %.0f\n", nb_relapse(a)))
  cat(sprintf("+ anti-GD2        : relapse day %.0f\n", nb_relapse(b)))

  print(nb_scen_dose())        # the two Fc arms against dose
  print(nb_scen_perm())        # permeability, not dose
  print(nb_scen_k322a())       # the Fc mutation
  print(nb_scen_burden())      # why it is an MRD therapy
  print(nb_scen_cytokine())    # IL-2 and GM-CSF
  print(nb_scen_retinoid())    # autoinduction, bioavailability, scheduling
  print(nb_scen_mibg())        # transporter competition and dosimetry
  print(nb_scen_alk())         # free fraction
  print(nb_scen_sts())         # ear versus tumour

  co <- nb_cohort(200)
  t_ctl <- nb_efs(co, immuno = FALSE)
  t_imm <- nb_efs(co, immuno = TRUE)
  cat(sprintf("2-year EFS  control %.1f%%   immunotherapy %.1f%%\n",
              100 * nb_efs_rate(t_ctl, 730), 100 * nb_efs_rate(t_imm, 730)))
}
