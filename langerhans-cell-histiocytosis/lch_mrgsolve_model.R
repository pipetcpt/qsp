# =============================================================================
# lch_mrgsolve_model.R
# Langerhans Cell Histiocytosis (LCH) — Quantitative Systems Pharmacology model
# =============================================================================
#
# 62 ODEs. Time unit = DAYS. Concentrations = mg/L. Lesional burden in
# "burden units" (1 unit ~ 1e9 LCH cells). ANC in 1e9/L.
#
# -----------------------------------------------------------------------------
# WHAT THIS MODEL COMMITS TO
# -----------------------------------------------------------------------------
# Textbook LCH models are "a clonal cell proliferates, drugs kill it". That
# picture cannot simultaneously reproduce six things the literature insists on:
#
#   (i)   the same driver (BRAF V600E) gives a self-healing skull lesion in one
#         child and fatal risk-organ disease in another;
#   (ii)  lesional Ki-67 is only ~5-10% (Brabencova 1998, PMID 9588881) yet the
#         lesion behaves like an inflammatory storm;
#   (iii) MAPK inhibitors produce clinical responses within DAYS but plasma
#         BRAF V600E cfDNA rarely clears (Eckstein 2019, PMID 30718231;
#         Evseev 2021, PMID 34383272);
#   (iv)  >75% reactivate within weeks-months of stopping a MAPK inhibitor
#         (Donadieu 2019, PMID 31513482; Cohen Aubart 2017, PMID 28667012);
#   (v)   cladribine works, even though the target cells barely divide;
#   (vi)  central DI and neurodegeneration are irreversible and correlate with
#         how long disease was active, not with how bad it got at its worst
#         (Grois 2006, PMID 16047354; Haupt 2004, PMID 15049016).
#
# So the model makes three structural commitments, each with a KILL SWITCH so
# it can be falsified:
#
# (1) THE CELL-OF-ORIGIN PARTITION, NOT ANY GROWTH RATE, SETS THE PHENOTYPE.
#     One mutant precursor pool is described by origin descriptors: FSR
#     (self-renewal competence of the stage at which the driver arose), PRECM0
#     (its size), and the seeding partition THB/THS/THR/THP/THC/THL. Every
#     LESIONAL kinetic constant (KPROL, KDL, KIMM, KDP, KSR, LMAX, KSEED, KDC)
#     is IDENTICAL across phenotypes. SS-b vs MS RO+ falls out of the origin
#     descriptors alone. Kill switch: give the SS-b phenotype the MS RO+
#     descriptors and it becomes risk-organ-positive with no rate changed.
#
# (2) LESIONS ARE MAINTAINED BY RECRUITMENT, AND MAPK INHIBITION IS CYTOSTATIC.
#     KPROL (0.022/day) is deliberately set BELOW the minimum immune clearance
#     rate KDL + KIMM*(1-TREGMAX) = 0.0288/day, so no lesion can sustain itself
#     by local division — it needs continuous supply from the precursor pool.
#     That is the quantitative content of the low Ki-67. MAPK inhibitors lower
#     proliferation (CCND), survival signal (BCL2A1) and the secretory program
#     (SASP) but add NO kill term (SL_MAPKI_KILL = 0); nucleoside analogues,
#     vinblastine and glucocorticoid do. Two consequences FALL OUT rather than
#     being coded: DAS collapses in days (secretome) while mass and cfDNA fall
#     over weeks and plateau; and withdrawal rebounds at the reservoir's own
#     regrowth rate. Kill switch: set SL_MAPKI_KILL > 0 and the rebound
#     disappears, contradicting (iv).
#
#     Two supporting mechanisms make the plateau real rather than fitted:
#       * FN_APO — the marrow niche supplies ERK-independent survival signal,
#         so losing BCL2A1 sensitises lesional cells much more than reservoir
#         cells. The reservoir therefore plateaus instead of decaying away.
#       * PNICHE — a minimum clone size below which self-renewal collapses (a
#         deterministic surrogate for stochastic clonal extinction). MAPK
#         inhibition never reaches it; six cycles of 2-CdA/Ara-C do, three do
#         not. Dose intensity therefore has a threshold, not a gradient.
#
# (3) PERMANENT SEQUELAE ARE TIME-INTEGRALS OF LOCAL ACTIVITY WITH THRESHOLDS.
#     AVPN (AVP-secreting neurons), ANTPIT, NEUR (cerebellar/pontine neurons),
#     BILF and LUNGC are MONOTONE — they never recover. CDI is declared when
#     AVPN < AVP_CRIT, clinical ND when NEUR < NEUR_CRIT. Therefore
#     time-to-effective-therapy (accumulator TTET) can dominate drug potency,
#     and the model COMPUTES whether a slower-but-deeper regimen loses to a
#     faster-but-shallower one rather than assuming it.
#
# -----------------------------------------------------------------------------
# VERIFICATION
# -----------------------------------------------------------------------------
# `lch_python_twin.py` in this directory re-implements this exact right-hand
# side and parameter block in dependency-free Python (RK4) and asserts 49
# quantitative checks. Run `python3 lch_python_twin.py`. If you edit this file,
# edit that one.
#
# Requires: mrgsolve, dplyr, ggplot2, tidyr
# =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)

lch_code <- '
$PARAM @annotated
// ---- 1. Patient / genotype / phenotype -------------------------------------
WT       :  12.0 : Body weight (kg) - drives allometric PK scaling
ALLO_CL  :  0.75 : Allometric exponent for clearances
ALLO_V   :  1.00 : Allometric exponent for volumes
GENO     :  1    : 1 = BRAF V600E, 2 = MAP2K1/ARAF/driver-negative
SMOKE    :  0    : Cigarette smoke exposure (0/1), pulmonary LCH driver
IL17ON   :  1    : Keep the contested IL-17A / DC-fusion arm (0/1)

// ---- 2. Cell-of-origin descriptors (STRUCTURAL COMMITMENT 1) ---------------
PRECM0   :  0.45 : Initial mutant marrow precursor pool (burden units)
FSR      :  1.00 : Self-renewal competence of the cell of origin (0-1)
PRMAX    :  3.00 : Niche carrying capacity of the marrow reservoir
PNICHE   :  0.0030 : Minimum clone size still supporting self-renewal
THB      :  0.28 : Seeding fraction to bone
THS      :  0.16 : Seeding fraction to skin
THR      :  0.34 : Seeding fraction to risk organs (liver/spleen/marrow)
THP      :  0.08 : Seeding fraction to pituitary / hypothalamus
THC      :  0.06 : Seeding fraction to CNS parenchyma
THL      :  0.08 : Seeding fraction to lung
LBONE0   :  0.05 : Locally seeded initial bone lesion
LSKIN0   :  0.02 : Locally seeded initial skin lesion
LRO0     :  0.02 : Locally seeded initial risk-organ lesion
LPIT0    :  0.01 : Locally seeded initial pituitary lesion
LCNS0    :  0.005 : Locally seeded initial CNS lesion
LLUNG0   :  0.005 : Locally seeded initial lung lesion

// ---- 3. Precursor and lesion kinetics (IDENTICAL ACROSS PHENOTYPES) --------
KSR      :  0.070 : Max self-renewal rate of the mutant precursor (1/day)
KDP      :  0.015 : Precursor baseline death (1/day)
KEX      :  0.080 : Progeny output to blood (1/day, proliferation-coupled)
KDC      :  1.00  : Circulating precursor clearance (1/day)
KSEED    :  0.50  : Circulating to tissue seeding (1/day)
KPROL    :  0.022 : Intralesional proliferation at CCND=1 (1/day), set below the minimum immune clearance rate
KDL      :  0.020 : Lesional baseline death (1/day)
KIMM     :  0.035 : Immune-mediated lesional clearance at TREG=0 (1/day)
KAPO     :  1.60  : Apoptosis amplification when BCL2A1 is lost
FN_APO   :  0.25  : Fraction of that amplification the reservoir feels
LMAX     :  8.0   : Per-site carrying capacity
FREC     :  0.45  : IL-1beta-driven CCR6 recruitment amplification
RECMAX   :  3.0   : Cap on recruitment amplification
LREF     :  2.0   : Burden units used to normalise lesion mass

// ---- 4. MAPK signalling ---------------------------------------------------
KERK     : 24.0   : pERK equilibration (1/day)
ERKB     :  0.05  : ERK activity floor (non-driver signalling)
KCC      :  2.00  : Cyclin D1 turnover (1/day)
KBCL     :  1.50  : BCL2A1 turnover (1/day)
KSASP    :  3.00  : Secretory program turnover (1/day) - FAST
EC50C    :  0.35  : ERK -> cyclin D1, EC50
HC       :  2.0   : ERK -> cyclin D1, Hill coefficient
EC50B    :  0.30  : ERK -> BCL2A1, EC50
HB       :  1.5   : ERK -> BCL2A1, Hill coefficient
EC50S    :  0.25  : ERK -> SASP, EC50
HS       :  1.2   : ERK -> SASP, Hill coefficient
PMAX     :  0.60  : Max paradoxical MAPK activation in BRAF-WT cells

// ---- 5. Vemurafenib PK/PD -------------------------------------------------
KA_VEM   :  8.0   : Absorption rate (1/day)
CL_VEM   : 32.0   : Apparent clearance at 70 kg (L/day)
V2_VEM   : 90.0   : Central volume at 70 kg (L)
Q_VEM    : 12.0   : Intercompartmental clearance at 70 kg (L/day)
V3_VEM   : 25.0   : Peripheral volume at 70 kg (L)
FU_VEM   :  0.005 : Unbound fraction (>99% protein bound)
IC50_VEM :  0.018 : Free conc for 50% BRAF inhibition (mg/L)
KENZ     :  0.10  : Autoinduction turnover (1/day)
EMAX_IND :  0.35  : Max fractional clearance autoinduction
EC50_IND : 30.0   : Concentration for half-maximal autoinduction (mg/L)

// ---- 6. Dabrafenib PK/PD --------------------------------------------------
KA_DAB   : 20.0   : Absorption rate (1/day)
CL_DAB   : 840.0  : Apparent clearance at 70 kg (L/day)
V_DAB    : 400.0  : Apparent volume at 70 kg (L)
FM_DAB   :  0.40  : Molar fraction converted to hydroxy-dabrafenib
CL_DABM  : 520.0  : Metabolite clearance at 70 kg (L/day)
V_DABM   : 180.0  : Metabolite volume at 70 kg (L)
RPOT_M   :  1.00  : Metabolite potency relative to parent
FU_DAB   :  0.005 : Unbound fraction
IC50_DAB :  0.00050 : Free conc for 50% BRAF inhibition (mg/L)
FCNS_DAB :  0.35  : CNS penetration factor (better than vemurafenib)

// ---- 7. Trametinib PK/PD --------------------------------------------------
KA_TRA   :  6.0   : Absorption rate (1/day)
CL_TRA   : 118.0  : Apparent clearance at 70 kg (L/day)
V2_TRA   : 250.0  : Central volume at 70 kg (L)
Q_TRA    : 25.0   : Intercompartmental clearance at 70 kg (L/day)
V3_TRA   : 900.0  : Deep peripheral volume at 70 kg (L) - gives t1/2 ~4 days
FU_TRA   :  0.026 : Unbound fraction
IC50_TRA :  0.00015 : Free conc for 50% MEK inhibition (mg/L)

// ---- 8. Cobimetinib PK/PD -------------------------------------------------
KA_COB   :  8.0    : Absorption rate (1/day)
CL_COB   : 322.0   : Apparent clearance at 70 kg (L/day)
V_COB    : 1000.0  : Apparent volume at 70 kg (L)
FU_COB   :  0.050  : Unbound fraction
IC50_COB :  0.0012 : Free conc for 50% MEK inhibition (mg/L)

// ---- 9. Vinblastine PK/PD -------------------------------------------------
CL_VBL   : 400.0  : Clearance at 70 kg (L/day)
V1_VBL   : 10.0   : Central volume at 70 kg (L)
Q_VBL    : 300.0  : Intercompartmental clearance at 70 kg (L/day)
V2_VBL   : 700.0  : Peripheral volume at 70 kg (L)
EC50_VBL :  0.0050 : Concentration for half-maximal tubulin effect (mg/L)
EMAX_VBL :  0.35  : Max mitotic-arrest kill rate (1/day)

// ---- 10. Prednisolone PK/PD ----------------------------------------------
KA_PRE   : 12.0   : Absorption rate (1/day)
CL_PRE   : 190.0  : Apparent clearance at 70 kg (L/day)
V_PRE    : 40.0   : Apparent volume at 70 kg (L)
KEO_P    :  4.0   : Effect-compartment equilibration (1/day)
EMAX_GR  :  0.70  : Max NF-kB transrepression of the secretome
EC50_GR  :  0.050 : Effect-site conc for half-maximal transrepression (mg/L)
EMAX_GRK :  0.28  : Max direct apoptotic kill rate (1/day)
EC50_GRK :  0.060 : Effect-site conc for half-maximal apoptosis (mg/L)

// ---- 11. Cytarabine PK/PD ------------------------------------------------
CL_ARAC  : 1915.0 : Clearance at 70 kg (L/day) - t1/2 ~20 min
V_ARAC   : 40.0   : Apparent volume at 70 kg (L)
KIN_ARAC :  1.0   : Intracellular ara-CTP formation rate
KOUT_ARATP : 4.5  : Intracellular ara-CTP elimination (1/day)
FCNS_ARAC : 0.15  : CNS exposure factor
EMAX_ARAC : 0.25  : Max ara-CTP kill rate (1/day)
EC50_ARAC : 0.020 : Intracellular ara-CTP for half-maximal kill

// ---- 12. Cladribine PK/PD ------------------------------------------------
CL_CLAD  : 150.0  : Clearance at 70 kg (L/day)
V_CLAD   : 55.0   : Apparent volume at 70 kg (L)
KIN_CLAD :  1.0   : Intracellular Cd-ATP formation rate
KOUT_CLATP : 0.70 : Intracellular Cd-ATP elimination (1/day)
FCNS_CLAD : 0.25  : CNS exposure factor (CSF ~25% of plasma)
EMAX_CLAD : 0.20  : Max Cd-ATP kill rate (1/day)
EC50_CLAD : 0.015 : Intracellular Cd-ATP for half-maximal kill

// ---- 13. Maintenance (6-mercaptopurine / methotrexate) -------------------
KA_MNT   :  2.0   : Absorption rate (1/day)
KE_MNT   :  1.0   : Elimination rate (1/day)
EMAX_MNT :  0.35  : Max fractional suppression of precursor self-renewal
EC50_MNT :  0.30  : Exposure for half-maximal suppression

// ---- 14. The cytostatic / cytotoxic split (STRUCTURAL COMMITMENT 2) ------
PFMAX    :  0.09  : Proliferative fraction of LCH cells at CCND=1 (Ki-67)
SL_MAPKI_KILL : 0.0 : KILL SWITCH - must stay 0 for commitment 2
W_ARAC_P :  0.55  : ara-CTP kill weight on the marrow precursor pool
W_CLAD_P :  1.40  : Cd-ATP kill weight on the precursor (cycling-independent)
W_VBL_P  :  0.30  : Vinblastine kill weight on the precursor
W_GR_P   :  0.35  : Glucocorticoid kill weight on the precursor

// ---- 15. Microenvironment -----------------------------------------------
KTR      :  0.25  : Treg equilibration (1/day)
TREGMAX  :  0.75  : Max lesional Treg index
EC50_TREG : 1.0   : Normalised burden for half-maximal Treg expansion
KOC      :  0.50  : Osteoclast pool equilibration (1/day)
KR_RANKL :  1.20  : RANKL concentration for half-maximal RANK signalling
AOSM     :  0.45  : Oncostatin-M amplification of osteoclastogenesis

// ---- 16. Secretome ------------------------------------------------------
KP_IL1   :  1.00  : IL-1beta production coefficient
KD_IL1   :  6.0   : IL-1beta elimination (1/day)
IL1B0    :  0.05  : IL-1beta baseline
KP_TNF   :  1.00  : TNF-alpha production coefficient
KD_TNF   :  6.0   : TNF-alpha elimination (1/day)
TNFA0    :  0.08  : TNF-alpha baseline
KP_IL6   :  1.20  : IL-6 production coefficient
KD_IL6   :  5.0   : IL-6 elimination (1/day)
IL60     :  0.10  : IL-6 baseline
KP_OSM   :  0.80  : Oncostatin-M production coefficient
KD_OSM   :  4.0   : Oncostatin-M elimination (1/day)
OSM0     :  0.05  : Oncostatin-M baseline
KP_MMP   :  1.00  : MMP-9 production coefficient
KD_MMP   :  2.0   : MMP-9 elimination (1/day)
MMP90    :  0.20  : MMP-9 baseline
KP_RNK   :  1.00  : RANKL production coefficient
KD_RNK   :  2.0   : RANKL elimination (1/day)
RANKL0   :  0.30  : RANKL baseline
W_TCELL_RNK : 0.35 : Lesional T-cell contribution to RANKL

// ---- 17. Biomarkers -----------------------------------------------------
KP_CF    : 14.0   : cfDNA production per unit mutant-cell death flux
KD_CF    :  8.0   : cfDNA elimination (1/day, t1/2 ~2 h)
CF_LOD   :  0.005 : ddPCR limit of detection (model units)
KP_163   :  1.0   : Soluble CD163 production coefficient
KD_163   :  1.0   : Soluble CD163 elimination (1/day)
SCD1630  :  1.0   : Soluble CD163 baseline
KP_CRP   :  1.6   : CRP production coefficient
KD_CRP   :  2.0   : CRP elimination (1/day)
CRP0     :  0.8   : CRP baseline (mg/L)
KP_FER   :  1.0   : Ferritin production coefficient
KD_FER   :  0.7   : Ferritin elimination (1/day)
FERR0    : 100.0  : Ferritin baseline (ng/mL)
W_HLH    :  2.5   : Risk-organ weight on the HLH-overlap ferritin term

// ---- 18. Irreversible organ pools (STRUCTURAL COMMITMENT 3) -------------
KAVP     :  0.010 : AVP-neuron loss rate per unit pituitary activity (1/day)
AVP_CRIT :  0.15  : AVP-neuron fraction below which central DI is manifest
KANT     :  0.008 : Anterior pituitary reserve loss rate (1/day)
KND      :  0.010 : Cerebellar/pontine neuron loss rate (1/day)
NEUR_CRIT : 0.70  : Neuron fraction below which clinical ND is manifest
W_CDI_ND :  0.35  : CDI as a neurodegeneration risk multiplier
KBIL     :  0.020 : Sclerosing-cholangitis fibrosis rate (1/day)
KCYST    :  0.0012 : Pulmonary cystic destruction rate (1/day)
ASMK     :  1.60  : Smoke amplification of pulmonary lesion growth and cysts
KRES     :  0.60  : Osteoclastic resorption rate into lesion volume
KHEAL    :  0.030 : Reossification rate (1/day)
KINH     :  3.0   : Inhibition of healing by residual LCH cells
BVMAX    : 30.0   : Maximum lytic lesion volume (cm3)

// ---- 19. Myelosuppression (Friberg) ------------------------------------
KTR_F    :  0.96  : Transit rate (1/day), MTT ~100 h
ANC0     :  4.0   : Baseline ANC (1e9/L)
GAM      :  0.16  : Feedback exponent
KCIRC    :  2.30  : Circulating neutrophil elimination (1/day)
SL_ARAC  :  9.0   : ara-CTP slope on proliferating marrow
SL_CLAD  :  3.0   : Cd-ATP slope on proliferating marrow
SL_VBL   : 25.0   : Vinblastine slope on proliferating marrow
MSUP     :  0.45  : Marrow suppression by risk-organ disease itself

// ---- 20. Targeted-therapy toxicity ------------------------------------
KSK      :  0.35  : Cutaneous paradox toxicity accrual (1/day)
KSKR     :  0.10  : Cutaneous toxicity resolution (1/day)
LVEF0    : 62.0   : Baseline LVEF (%)
KLV      :  0.65  : MEK-inhibitor LVEF decrement rate
KREC_LV  :  0.06  : LVEF recovery rate (1/day)

// ---- 21. Endpoint definitions ----------------------------------------
DAS_ACTIVE : 3.0  : DAS above which disease counts as active
DAS_DX     : 3.5  : DAS at which the patient presents
W_DAS_BONE : 1.4  : DAS weight, bone
W_DAS_SKIN : 0.9  : DAS weight, skin
W_DAS_RO   : 3.2  : DAS weight, risk organs
W_DAS_PIT  : 1.6  : DAS weight, pituitary
W_DAS_CNS  : 1.8  : DAS weight, CNS
W_DAS_LUNG : 1.5  : DAS weight, lung
W_DAS_INFL : 4.0  : DAS weight, systemic inflammation

$CMT @annotated
// ---- drug PK (22) ---------------------------------------------------------
VEMG   : Vemurafenib depot (mg)
VEMC   : Vemurafenib central (mg)
VEMP   : Vemurafenib peripheral (mg)
VEMI   : CYP3A4 autoinduction state (relative, 1 = uninduced)
DABG   : Dabrafenib depot (mg)
DABC   : Dabrafenib central (mg)
DABM   : Hydroxy-dabrafenib (mg)
TRAG   : Trametinib depot (mg)
TRAC   : Trametinib central (mg)
TRAP   : Trametinib deep peripheral (mg)
COBG   : Cobimetinib depot (mg)
COBC   : Cobimetinib central (mg)
VBLC   : Vinblastine central (mg)
VBLP   : Vinblastine peripheral (mg)
PREG   : Prednisolone depot (mg)
PREC   : Prednisolone central (mg)
ARAC   : Cytarabine central (mg)
ARATP  : Intracellular ara-CTP (relative)
CLAC   : Cladribine central (mg)
CLATP  : Intracellular Cd-ATP (relative)
MNTG   : 6-MP/MTX maintenance depot (relative)
MNTC   : 6-MP/MTX maintenance exposure (relative)
// ---- signal transduction (5) ---------------------------------------------
ERK    : pERK activity (1 = untreated LCH lesion)
CCND   : Cyclin D1 / proliferative drive (relative)
BCL    : BCL2A1 / MCL-1 survival signal (relative)
SASP   : Senescence-associated secretory program (relative)
GRE    : Glucocorticoid effect site (mg/L)
// ---- cell pools (10) -----------------------------------------------------
PRECM  : Mutant marrow precursor reservoir (burden units)
CIRC   : Circulating mutant DC precursors (burden units)
LBONE  : Bone lesional burden (burden units)
LSKIN  : Skin lesional burden
LRO    : Risk-organ lesional burden (liver/spleen/marrow)
LPIT   : Pituitary / hypothalamic lesional burden
LCNS   : CNS parenchymal mutant myeloid burden
LLUNG  : Pulmonary lesional burden
TREG   : Lesional regulatory T-cell index (0-1)
OCL    : Activated osteoclast pool (relative)
// ---- secretome (6) -------------------------------------------------------
IL1B   : Interleukin-1 beta (relative)
TNFA   : Tumour necrosis factor alpha (relative)
IL6    : Interleukin-6 (relative)
OSM    : Oncostatin M (relative)
MMP9   : MMP-9 / MMP-12 (relative)
RANKL  : Effective RANKL to OPG ratio (relative)
// ---- biomarkers (4) ------------------------------------------------------
CFDNA  : Plasma cell-free BRAF V600E signal
SCD163 : Soluble CD163 (relative)
CRP    : C-reactive protein (mg/L)
FERR   : Ferritin (ng/mL)
// ---- irreversible organ pools (6) ----------------------------------------
AVPN   : AVP-secreting magnocellular neuron pool (fraction, 1 -> 0)
ANTPIT : Anterior pituitary reserve (fraction)
NEUR   : Cerebellar / pontine neuron pool (fraction)
BILF   : Biliary fibrosis / sclerosing cholangitis (0 -> 1)
LUNGC  : Pulmonary cystic destruction (0 -> 1)
BVOL   : Osteolytic lesion volume (cm3)
// ---- myelosuppression, Friberg (4) ---------------------------------------
PROL   : Proliferating marrow progenitors (1e9/L equivalent)
TR1    : Maturation transit compartment 1
TR2    : Maturation transit compartment 2
ANC    : Absolute neutrophil count (1e9/L)
// ---- targeted-therapy toxicity (2) ---------------------------------------
SKTOX  : Cutaneous paradox toxicity index
LVEF   : Left ventricular ejection fraction (%)
// ---- endpoint accumulators (3) -------------------------------------------
AUCERK : Days spent with >80% pERK suppression
TTET   : Days of uncontrolled active disease (time to effective therapy)
CUMDAS : Cumulative disease-activity exposure (DAS-days)

$MAIN
  VEMI_0   = 1.0;
  ERK_0    = 1.0;
  CCND_0   = 1.0;
  BCL_0    = 1.0;
  SASP_0   = 1.0;
  PRECM_0  = PRECM0;
  // circulating precursors start at their quasi-steady state
  CIRC_0   = KEX * PRECM0 / (KSEED + KDC);
  LBONE_0  = LBONE0;
  LSKIN_0  = LSKIN0;
  LRO_0    = LRO0;
  LPIT_0   = LPIT0;
  LCNS_0   = LCNS0;
  LLUNG_0  = LLUNG0;
  TREG_0   = 0.05;
  OCL_0    = 0.10;
  IL1B_0   = IL1B0;
  TNFA_0   = TNFA0;
  IL6_0    = IL60;
  OSM_0    = OSM0;
  MMP9_0   = MMP90;
  RANKL_0  = RANKL0;
  SCD163_0 = SCD1630;
  CRP_0    = CRP0;
  FERR_0   = FERR0;
  AVPN_0   = 1.0;
  ANTPIT_0 = 1.0;
  NEUR_0   = 1.0;
  PROL_0   = ANC0;
  TR1_0    = ANC0;
  TR2_0    = ANC0;
  ANC_0    = ANC0;
  LVEF_0   = LVEF0;
  // cfDNA starts at the quasi-steady state implied by the baseline
  // mutant-cell death flux (cfDNA t1/2 ~2 h, so it tracks flux closely)
  double lt0  = LBONE0 + LSKIN0 + LRO0 + LPIT0 + LCNS0 + LLUNG0;
  double clr0 = KDL + KIMM * (1.0 - 0.05);
  CFDNA_0 = KP_CF * (KDP * PRECM0 + clr0 * lt0 + KDC * CIRC_0) / KD_CF;

$GLOBAL
// f(x)/f(1): a Hill transducer normalised so that ERK = 1 maps to 1
#define HILLN(x, ec50, h) ( ((x) > 0 ? pow((x),(h)) / (pow((x),(h)) + pow((ec50),(h))) : 0.0) \\
                            * (1.0 + pow((ec50),(h))) )
// differentiable indicator of (x > thr)
#define SSTEP(x, thr, w) ( 1.0 / (1.0 + exp(-(((x)-(thr))/(w)))) )

$ODE
  // =========================================================================
  // Allometric scaling. All PK parameters are quoted at 70 kg; without this a
  // 12 kg infant on 20 mg/kg/day vemurafenib would be predicted to sit at
  // ~7 mg/L instead of the observed 20-60 mg/L.
  // =========================================================================
  double fcl = pow(WT / 70.0, ALLO_CL);
  double fv  = pow(WT / 70.0, ALLO_V);

  // ---- vemurafenib --------------------------------------------------------
  double CVEM  = VEMC / (V2_VEM * fv);
  double CVEMP = VEMP / (V3_VEM * fv);
  double QV    = Q_VEM * fcl;
  double CLVEM = CL_VEM * fcl * VEMI;
  dxdt_VEMG = -KA_VEM * VEMG;
  dxdt_VEMC = KA_VEM * VEMG - CLVEM * CVEM - QV * CVEM + QV * CVEMP;
  dxdt_VEMP = QV * CVEM - QV * CVEMP;
  dxdt_VEMI = KENZ * ((1.0 + EMAX_IND * CVEM / (EC50_IND + CVEM)) - VEMI);

  // ---- dabrafenib + active hydroxy metabolite -----------------------------
  double CDAB  = DABC / (V_DAB * fv);
  double CDABM = DABM / (V_DABM * fv);
  dxdt_DABG = -KA_DAB * DABG;
  dxdt_DABC = KA_DAB * DABG - CL_DAB * fcl * CDAB;
  dxdt_DABM = FM_DAB * CL_DAB * fcl * CDAB - CL_DABM * fcl * CDABM;

  // ---- trametinib ---------------------------------------------------------
  double CTRA  = TRAC / (V2_TRA * fv);
  double CTRAP = TRAP / (V3_TRA * fv);
  double QT    = Q_TRA * fcl;
  dxdt_TRAG = -KA_TRA * TRAG;
  dxdt_TRAC = KA_TRA * TRAG - CL_TRA * fcl * CTRA - QT * CTRA + QT * CTRAP;
  dxdt_TRAP = QT * CTRA - QT * CTRAP;

  // ---- cobimetinib --------------------------------------------------------
  double CCOB = COBC / (V_COB * fv);
  dxdt_COBG = -KA_COB * COBG;
  dxdt_COBC = KA_COB * COBG - CL_COB * fcl * CCOB;

  // ---- vinblastine --------------------------------------------------------
  double CVBL  = VBLC / (V1_VBL * fv);
  double CVBLP = VBLP / (V2_VBL * fv);
  double QB    = Q_VBL * fcl;
  dxdt_VBLC = -CL_VBL * fcl * CVBL - QB * CVBL + QB * CVBLP;
  dxdt_VBLP = QB * CVBL - QB * CVBLP;

  // ---- prednisolone -------------------------------------------------------
  double CPRE = PREC / (V_PRE * fv);
  dxdt_PREG = -KA_PRE * PREG;
  dxdt_PREC = KA_PRE * PREG - CL_PRE * fcl * CPRE;
  dxdt_GRE  = KEO_P * (CPRE - GRE);

  // ---- cytarabine / cladribine and their intracellular triphosphates -----
  double CARAC = ARAC / (V_ARAC * fv);
  dxdt_ARAC  = -CL_ARAC * fcl * CARAC;
  dxdt_ARATP = KIN_ARAC * CARAC - KOUT_ARATP * ARATP;
  double CCLAD = CLAC / (V_CLAD * fv);
  dxdt_CLAC  = -CL_CLAD * fcl * CCLAD;
  dxdt_CLATP = KIN_CLAD * CCLAD - KOUT_CLATP * CLATP;

  // ---- maintenance --------------------------------------------------------
  dxdt_MNTG = -KA_MNT * MNTG;
  dxdt_MNTC = KA_MNT * MNTG - KE_MNT * MNTC;

  // =========================================================================
  // Target engagement. Occupancies are additive on their shared node, then
  // saturated - so combining two BRAF inhibitors is not double-counted.
  // =========================================================================
  double RV = FU_VEM * CVEM / IC50_VEM;
  double RD = FU_DAB * (CDAB + RPOT_M * CDABM) / IC50_DAB;
  double RB = RV + RD;
  double RT = FU_TRA * CTRA / IC50_TRA;
  double RC = FU_COB * CCOB / IC50_COB;
  double RM = RT + RC;

  double IMEK = RM / (1.0 + RM);
  double braf_active = (GENO < 1.5) ? 1.0 : 0.0;
  double IBRAF = braf_active * RB / (1.0 + RB);
  // Paradoxical MAPK activation: a BRAF inhibitor bound to a cell WITHOUT the
  // V600E monomer transactivates CRAF dimers. In BRAF-WT keratinocytes this
  // drives keratoacanthoma/cuSCC; in MAP2K1-driven LCH it drives the disease.
  // A MEK inhibitor sits downstream and suppresses both.
  double PARADOX      = (1.0 - braf_active) * PMAX * (RB / (1.0 + RB)) * (1.0 - IMEK);
  double PARADOX_SKIN = PMAX * (RB / (1.0 + RB)) * (1.0 - IMEK);

  double ERKt = (1.0 - IBRAF) * (1.0 - IMEK) * (1.0 + PARADOX) + ERKB;
  if (ERKt > 1.0 + PMAX) ERKt = 1.0 + PMAX;
  dxdt_ERK = KERK * (ERKt - ERK);

  double ERKs = (ERK > 1e-9) ? ERK : 1e-9;
  dxdt_CCND = KCC  * (HILLN(ERKs, EC50C, HC) - CCND);
  dxdt_BCL  = KBCL * (HILLN(ERKs, EC50B, HB) - BCL);
  dxdt_SASP = KSASP * (HILLN(ERKs, EC50S, HS) - SASP);

  // =========================================================================
  // Drug effects. NOTE the asymmetry that is the whole point of commitment 2:
  // the MAPK-inhibitor terms appear ONLY through CCND / BCL / SASP, never in
  // KILL_LES or KILL_PRE (SL_MAPKI_KILL is a kill switch pinned at 0).
  // =========================================================================
  double FGR      = 1.0 - EMAX_GR * GRE / (EC50_GR + GRE);
  double EGR_KILL = EMAX_GRK * GRE / (EC50_GRK + GRE);
  double E_ARAC   = EMAX_ARAC * ARATP / (EC50_ARAC + ARATP);
  double E_CLAD   = EMAX_CLAD * CLATP / (EC50_CLAD + CLATP);
  double E_VBL    = EMAX_VBL * CVBL / (EC50_VBL + CVBL);
  double E_MNT    = EMAX_MNT * MNTC / (EC50_MNT + MNTC);

  double PFN = CCND;   // cycling fraction normalised to CCND = 1

  // Cytarabine and vinblastine need cells in S phase or mitosis, so they are
  // gated by PFN. Cladribine is NOT: Cd-ATP kills non-dividing monocytoid
  // cells. With a Ki-67 of only ~5-10% that difference is the pharmacological
  // reason 2-CdA works in LCH at all.
  double KILL_LES = E_ARAC * PFN + E_VBL * PFN + E_CLAD + EGR_KILL
                  + SL_MAPKI_KILL * (IBRAF + IMEK);
  double KILL_PRE = W_ARAC_P * E_ARAC * PFN + W_VBL_P * E_VBL * PFN
                  + W_CLAD_P * E_CLAD + W_GR_P * EGR_KILL
                  + SL_MAPKI_KILL * (IBRAF + IMEK);
  double KILL_CNS = E_ARAC * FCNS_ARAC * PFN + E_CLAD * FCNS_CLAD + EGR_KILL;

  double APO = 1.0 + KAPO * (1.0 - BCL);
  double CLR = KDL * APO + KIMM * (1.0 - TREG);

  // =========================================================================
  // Cell pools
  // =========================================================================
  // The mutant precursor divides asymmetrically: one daughter self-renews in
  // the marrow niche, the other is exported. Export therefore does not debit
  // the reservoir - which is why a purely cytostatic drug cannot empty it.
  // The niche also shields the reservoir from BCL2A1-loss apoptosis (FN_APO),
  // and self-renewal collapses below a minimum clone size (PNICHE).
  double SR = KSR * FSR * CCND * (1.0 - PRECM / PRMAX) * (1.0 - E_MNT)
              * PRECM / (PRECM + PNICHE);
  double APO_PRE = 1.0 + FN_APO * KAPO * (1.0 - BCL);
  dxdt_PRECM = (SR - KDP * APO_PRE - KILL_PRE) * PRECM;
  dxdt_CIRC  = KEX * CCND * PRECM - (KSEED + KDC + KILL_PRE) * CIRC;

  double LTOT  = LBONE + LSKIN + LRO + LPIT + LCNS + LLUNG;
  double LTOTn = LTOT / LREF;
  double REC   = 1.0 + FREC * (IL1B / IL1B0 - 1.0);
  if (REC < 0.2)    REC = 0.2;
  if (REC > RECMAX) REC = RECMAX;
  double seed = KSEED * CIRC * REC;
  double smoke_boost = ASMK * SMOKE;

  dxdt_LBONE = THB * seed + KPROL * CCND * LBONE * (1.0 - LBONE / LMAX)
               - CLR * LBONE - KILL_LES * LBONE;
  dxdt_LSKIN = THS * seed + KPROL * CCND * LSKIN * (1.0 - LSKIN / LMAX)
               - CLR * LSKIN - KILL_LES * LSKIN;
  dxdt_LRO   = THR * seed + KPROL * CCND * LRO * (1.0 - LRO / LMAX)
               - CLR * LRO - KILL_LES * LRO;
  dxdt_LPIT  = THP * seed + KPROL * CCND * LPIT * (1.0 - LPIT / LMAX)
               - CLR * LPIT - KILL_LES * LPIT;
  dxdt_LCNS  = THC * seed + KPROL * CCND * LCNS * (1.0 - LCNS / LMAX)
               - CLR * LCNS - KILL_CNS * LCNS;
  dxdt_LLUNG = THL * seed
               + KPROL * CCND * LLUNG * (1.0 - LLUNG / LMAX) * (1.0 + smoke_boost)
               - CLR * LLUNG - KILL_LES * LLUNG;

  dxdt_TREG = KTR * (TREGMAX * LTOTn / (EC50_TREG + LTOTn) - TREG);

  double rnk = RANKL / (RANKL + KR_RANKL);
  double oc_target = rnk * (1.0 + AOSM * OSM / OSM0 * 0.25) * (1.0 + 0.30 * IL17ON);
  dxdt_OCL = KOC * (oc_target - OCL);

  // =========================================================================
  // Secretome: ERK/SASP-driven, proportional to active lesional mass, and
  // suppressible by glucocorticoid NF-kB transrepression (FGR).
  // =========================================================================
  double SRC = SASP * LTOTn * FGR;
  dxdt_IL1B  = KP_IL1 * SRC - KD_IL1 * (IL1B - IL1B0);
  dxdt_TNFA  = KP_TNF * SRC - KD_TNF * (TNFA - TNFA0);
  dxdt_IL6   = KP_IL6 * SRC - KD_IL6 * (IL6 - IL60);
  dxdt_OSM   = KP_OSM * SRC - KD_OSM * (OSM - OSM0);
  dxdt_MMP9  = KP_MMP * SRC - KD_MMP * (MMP9 - MMP90);
  dxdt_RANKL = KP_RNK * SRC * (1.0 + W_TCELL_RNK) - KD_RNK * (RANKL - RANKL0);

  // =========================================================================
  // Biomarkers. cfDNA tracks the DEATH FLUX of mutant cells, not their mass -
  // which is exactly why a cytostatic drug produces a shallow plateau and a
  // nucleoside analogue produces a deep nadir.
  // =========================================================================
  double death_flux = KDP * APO_PRE * PRECM + KILL_PRE * PRECM
                    + (CLR + KILL_LES) * (LTOT - LCNS)
                    + (CLR + KILL_CNS) * LCNS
                    + KDC * CIRC;
  dxdt_CFDNA  = KP_CF * death_flux - KD_CF * CFDNA;
  dxdt_SCD163 = KP_163 * LTOTn * FGR - KD_163 * (SCD163 - SCD1630);
  dxdt_CRP    = KP_CRP * IL6 / IL60 - KD_CRP * CRP;
  double hlh  = W_HLH * (LRO / LREF);
  dxdt_FERR   = KP_FER * (100.0 * LTOTn + 400.0 * hlh) - KD_FER * (FERR - FERR0);

  // =========================================================================
  // Irreversible organ pools: monotone integrals of LOCAL activity. No term
  // in any of these five equations can be positive for a recovering patient.
  // =========================================================================
  double infl = IL1B / IL1B0;
  dxdt_AVPN   = -KAVP * (LPIT / LREF) * sqrt(infl) * AVPN;
  dxdt_ANTPIT = -KANT * (LPIT / LREF) * ANTPIT;
  double cdi_flag = SSTEP(AVP_CRIT, AVPN, 0.02);
  dxdt_NEUR   = -KND * ((LCNS / LREF)
                        + W_CDI_ND * cdi_flag * (LCNS / LREF + 0.15)) * NEUR;
  dxdt_BILF   = KBIL * (LRO / LREF) * (1.0 - BILF);
  dxdt_LUNGC  = KCYST * (LLUNG / LREF) * (1.0 + smoke_boost) * (1.0 - LUNGC);
  dxdt_BVOL   = KRES * OCL * (LBONE / LREF) * (1.0 - BVOL / BVMAX)
                - KHEAL * BVOL / (1.0 + KINH * LBONE / LREF);

  // =========================================================================
  // Myelosuppression: Friberg semi-mechanistic model, with the extra fact
  // that risk-organ LCH suppresses the marrow by itself.
  // =========================================================================
  double EDRUG = SL_ARAC * ARATP + SL_CLAD * CLATP + SL_VBL * CVBL;
  if (EDRUG > 0.95) EDRUG = 0.95;
  double dis_sup = MSUP * LRO / LREF;
  if (dis_sup > 0.8) dis_sup = 0.8;
  double ancs = (ANC > 1e-3) ? ANC : 1e-3;
  double fb = pow(ANC0 / ancs, GAM);
  dxdt_PROL = KTR_F * PROL * (1.0 - EDRUG) * (1.0 - dis_sup) * fb
              - KTR_F * PROL;
  dxdt_TR1  = KTR_F * (PROL - TR1);
  dxdt_TR2  = KTR_F * (TR1 - TR2);
  dxdt_ANC  = KTR_F * TR2 - KCIRC * ANC;

  // ---- targeted-therapy toxicity -----------------------------------------
  dxdt_SKTOX = KSK * PARADOX_SKIN - KSKR * SKTOX;
  dxdt_LVEF  = KREC_LV * (LVEF0 - LVEF) - KLV * IMEK;

  // ---- endpoint accumulators ---------------------------------------------
  double dasx = W_DAS_BONE * (LBONE / LREF) / (1.0 + LBONE / LREF)
              + W_DAS_SKIN * (LSKIN / LREF) / (1.0 + LSKIN / LREF)
              + W_DAS_RO   * (LRO   / LREF) / (1.0 + LRO   / LREF)
              + W_DAS_PIT  * (LPIT  / LREF) / (1.0 + LPIT  / LREF)
              + W_DAS_CNS  * (LCNS  / LREF) / (1.0 + LCNS  / LREF)
              + W_DAS_LUNG * (LLUNG / LREF) / (1.0 + LLUNG / LREF);
  double inflx = infl - 1.0;
  if (inflx < 0.0) inflx = 0.0;
  dasx = dasx + W_DAS_INFL * inflx / (1.0 + inflx);

  dxdt_AUCERK = SSTEP(0.20, ERK, 0.02);
  dxdt_TTET   = SSTEP(dasx, DAS_ACTIVE, 0.3);
  dxdt_CUMDAS = dasx;

$TABLE
  double fclT = pow(WT / 70.0, ALLO_CL);
  double fvT  = pow(WT / 70.0, ALLO_V);
  capture CVEM_mgL  = VEMC / (V2_VEM * fvT);
  capture CDAB_mgL  = DABC / (V_DAB * fvT);
  capture CTRA_ngmL = TRAC / (V2_TRA * fvT) * 1000.0;
  capture CCOB_mgL  = COBC / (V_COB * fvT);
  capture CVBL_mgL  = VBLC / (V1_VBL * fvT);
  capture CPRE_mgL  = PREC / (V_PRE * fvT);
  capture pERK_pct  = 100.0 * ERK;
  capture LTOT_units = LBONE + LSKIN + LRO + LPIT + LCNS + LLUNG;

  double inflT = IL1B / IL1B0 - 1.0;
  if (inflT < 0.0) inflT = 0.0;
  capture DAS = W_DAS_BONE * (LBONE / LREF) / (1.0 + LBONE / LREF)
              + W_DAS_SKIN * (LSKIN / LREF) / (1.0 + LSKIN / LREF)
              + W_DAS_RO   * (LRO   / LREF) / (1.0 + LRO   / LREF)
              + W_DAS_PIT  * (LPIT  / LREF) / (1.0 + LPIT  / LREF)
              + W_DAS_CNS  * (LCNS  / LREF) / (1.0 + LCNS  / LREF)
              + W_DAS_LUNG * (LLUNG / LREF) / (1.0 + LLUNG / LREF)
              + W_DAS_INFL * inflT / (1.0 + inflT);

  capture CDI       = (AVPN < AVP_CRIT)  ? 1.0 : 0.0;   // central diabetes insipidus
  capture ND        = (NEUR < NEUR_CRIT) ? 1.0 : 0.0;   // clinical neurodegeneration
  capture GHD       = (ANTPIT < 0.55)    ? 1.0 : 0.0;   // GH deficiency
  capture CHOLANG   = (BILF > 0.50)      ? 1.0 : 0.0;   // sclerosing cholangitis
  capture G4NEUT    = (ANC < 0.50)       ? 1.0 : 0.0;   // grade 4 neutropenia
  capture CFDETECT  = (CFDNA > CF_LOD)   ? 1.0 : 0.0;   // cfDNA above assay LOD
  capture RESERVOIR = (PRECM > PNICHE)   ? 1.0 : 0.0;   // clone still viable
  capture ACTIVE    = (DAS > DAS_ACTIVE) ? 1.0 : 0.0;
'

mod <- mcode("lch_qsp", lch_code, atol = 1e-8, rtol = 1e-8, maxsteps = 200000)

# =============================================================================
# PATIENT PHENOTYPES
# -----------------------------------------------------------------------------
# NOTE WHAT DOES AND DOES NOT CHANGE. Only the cell-of-origin descriptors
# (FSR, PRECM0, the seeding thetas, the locally seeded lesions), body weight
# and the smoking flag differ. Every lesional kinetic constant is untouched.
# That is structural commitment (1) in executable form.
# =============================================================================
pheno <- list(
  SSb = list(   # single-system bone; driver arose in a committed tissue
                # precursor, so the clone cannot self-renew (FSR = 0.10)
    PRECM0 = 0.020, FSR = 0.10, THB = 0.70, THS = 0.10, THR = 0.04,
    THP = 0.06, THC = 0.02, THL = 0.08, LBONE0 = 0.80, LSKIN0 = 0,
    LRO0 = 0, LPIT0 = 0, LCNS0 = 0, LLUNG0 = 0
  ),
  MSROneg = list(
    PRECM0 = 0.12, FSR = 0.72, THB = 0.42, THS = 0.22, THR = 0.06,
    THP = 0.14, THC = 0.06, THL = 0.10, LRO0 = 0
  ),
  MSROpos = list(   # driver arose at the HSC stage: full self-renewal, heavy
                    # risk-organ seeding, age < 2 years
    PRECM0 = 0.45, FSR = 1.00, THB = 0.28, THS = 0.16, THR = 0.34,
    THP = 0.08, THC = 0.06, THL = 0.08
  ),
  CNSrisk = list(   # craniofacial CNS-risk lesion with pituitary involvement
    PRECM0 = 0.25, FSR = 0.85, THB = 0.30, THS = 0.10, THR = 0.10,
    THP = 0.30, THC = 0.14, THL = 0.06, LPIT0 = 0.06
  ),
  PLCH = list(      # adult pulmonary LCH in a smoker
    PRECM0 = 0.06, FSR = 0.24, THB = 0.10, THS = 0.04, THR = 0.02,
    THP = 0.02, THC = 0.02, THL = 0.80, SMOKE = 1, LLUNG0 = 0.20,
    LBONE0 = 0, LSKIN0 = 0, LRO0 = 0, LPIT0 = 0, LCNS0 = 0, WT = 65
  )
)

# Run-in to presentation. Response and time-to-control are only meaningful
# against the burden the patient actually presents with, so every scenario
# starts from an established disease state rather than from t = 0 seeds.
runin <- list(SSb = 30, MSROneg = 500, MSROpos = 500, CNSrisk = 500,
              PLCH = 365)

establish <- function(kind, das_dx = 3.5) {
  p <- pheno[[kind]]
  m <- param(mod, p)
  cap <- runin[[kind]]
  if (kind %in% c("SSb", "PLCH")) {
    out <- m %>% mrgsim(end = cap, delta = 1, recsort = 3) %>% as.data.frame()
    y <- out[nrow(out), ]
  } else {
    out <- m %>% mrgsim(end = cap, delta = 1, recsort = 3) %>% as.data.frame()
    hit <- which(out$DAS >= das_dx)
    y <- if (length(hit)) out[hit[1], ] else out[nrow(out), ]
  }
  # organ damage accrued before diagnosis is KEPT; the day counters are zeroed
  y$AUCERK <- 0; y$TTET <- 0; y$CUMDAS <- 0
  list(p = p, init = y, t_dx = y$time)
}

init_from <- function(m, y) {
  cmts <- names(init(m))
  vals <- as.list(y[, intersect(cmts, names(y)), drop = FALSE])
  do.call(init, c(list(m), vals))
}

# =============================================================================
# DOSING REGIMENS
# =============================================================================
BSA <- function(wt) 0.024265 * (wt ^ 0.5378) * (85 ^ 0.3964)   # Haycock-like

# LCH-III style front line: vinblastine 6 mg/m2 IV weekly + prednisolone
# 40 mg/m2/day x 4 weeks, taper, then 3-weekly pulses to 12 months, with
# 6-MP/MTX maintenance during continuation.
r_vblpred <- function(wt, start = 0, weeks_ind = 6, cont_to = 365) {
  bsa <- BSA(wt); vbl <- 6 * bsa; pre <- 40 * bsa
  ev_vbl <- ev(time = start, amt = vbl, cmt = "VBLC",
               ii = 7, addl = weeks_ind - 1)
  ev_pre <- ev(time = start, amt = pre, cmt = "PREG", ii = 0.5, addl = 55)
  ev_tap <- ev(time = start + 28, amt = pre * 0.5, cmt = "PREG",
               ii = 1, addl = 13)
  out <- c(ev_vbl, ev_pre, ev_tap)
  t <- start + 42
  while (t < cont_to) {
    out <- c(out,
             ev(time = t, amt = vbl, cmt = "VBLC"),
             ev(time = t, amt = pre, cmt = "PREG", ii = 0.5, addl = 9))
    t <- t + 21
  }
  c(out, ev(time = start + 42, amt = 1, cmt = "MNTG",
            ii = 1, addl = max(0, cont_to - start - 43)))
}

# Salvage: cladribine 9 mg/m2/day + cytarabine 500 mg/m2/day, 5-day CIV, q28d
# (Donadieu 2015, PMID 26194764).
r_cladarac <- function(wt, start = 0, cycles = 6, clad = 9, arac = 500,
                       cyc = 28) {
  bsa <- BSA(wt); out <- NULL
  for (i in seq_len(cycles)) {
    s <- start + (i - 1) * cyc
    out <- c(out,
             ev(time = s, amt = clad * bsa * 5, cmt = "CLAC",
                rate = clad * bsa),
             ev(time = s, amt = arac * bsa * 5, cmt = "ARAC",
                rate = arac * bsa))
  }
  out
}

r_maint <- function(start, stop) {
  ev(time = start, amt = 1, cmt = "MNTG", ii = 1,
     addl = max(0, stop - start - 1))
}

# Vemurafenib 20 mg/kg/day divided BID (Donadieu 2019, PMID 31513482)
r_vem <- function(wt, start = 0, stop = 730, mgkg = 20) {
  ev(time = start, amt = mgkg * wt / 2, cmt = "VEMG",
     ii = 0.5, addl = max(0, (stop - start) * 2 - 1))
}

# Dabrafenib + trametinib (Whitlock 2023, PMID 36884302)
r_dabtram <- function(wt, start = 0, stop = 730, dab = 4.5, tram = 0.032) {
  c(ev(time = start, amt = dab * wt / 2, cmt = "DABG",
       ii = 0.5, addl = max(0, (stop - start) * 2 - 1)),
    ev(time = start, amt = tram * wt, cmt = "TRAG",
       ii = 1, addl = max(0, stop - start - 1)))
}

# =============================================================================
# THE TEN SCENARIOS
# =============================================================================
scenarios <- list(
  S1_observation = list(
    kind = "SSb", end = 730, dose = function(wt) NULL,
    label = "1. SS-b, observation only (self-limited)"),
  S2_LCH3_ROneg = list(
    kind = "MSROneg", end = 730,
    dose = function(wt) r_vblpred(wt, 0, 6, 365),
    label = "2. MS RO-negative, LCH-III VBL/prednisolone 12 months"),
  S3_frontline_fail = list(
    kind = "MSROpos", end = 730,
    dose = function(wt) c(r_vblpred(wt, 0, 6, 42),
                          r_cladarac(wt, 42, 5), r_maint(180, 545)),
    label = "3. MS RO+ week-6 non-responder, salvage 2-CdA/Ara-C"),
  S4_cladarac_upfront = list(
    kind = "MSROpos", end = 730,
    dose = function(wt) c(r_cladarac(wt, 0, 6), r_maint(175, 545)),
    label = "4. MS RO+ refractory, 2-CdA/Ara-C x6 up front"),
  S5_vem_continuous = list(
    kind = "MSROpos", end = 730, dose = function(wt) r_vem(wt, 0, 730),
    label = "5. BRAF V600E MS RO+, vemurafenib continuous"),
  S6_vem_stop = list(
    kind = "MSROpos", end = 730, dose = function(wt) r_vem(wt, 0, 365),
    label = "6. Vemurafenib 12 months then STOP (rebound)"),
  S7_dabtram = list(
    kind = "MSROpos", end = 730, dose = function(wt) r_dabtram(wt, 0, 730),
    label = "7. Dabrafenib + trametinib continuous"),
  S8_bridge_consolidate = list(
    kind = "MSROpos", end = 730,
    dose = function(wt) c(r_vem(wt, 0, 112), r_cladarac(wt, 56, 3),
                          r_maint(145, 420)),
    label = "8. MAPKi bridge 8 wk -> 2-CdA/Ara-C consolidation -> stop all"),
  S9_delayed_dx = list(
    kind = "CNSrisk", end = 730,
    dose = function(wt) r_vblpred(wt, 180, 6, 545),
    label = "9a. CNS-risk lesion, therapy delayed 6 months"),
  S9b_early_dx = list(
    kind = "CNSrisk", end = 730,
    dose = function(wt) r_vblpred(wt, 14, 6, 379),
    label = "9b. Same patient, therapy at day 14 (comparator)"),
  S10_plch_smoke = list(
    kind = "PLCH", end = 1095, dose = function(wt) NULL,
    label = "10a. Adult pulmonary LCH, continued smoking"),
  S10b_plch_quit = list(
    kind = "PLCH", end = 1095, dose = function(wt) NULL, quit = 90,
    label = "10b. Adult pulmonary LCH, smoking cessation at month 3")
)

run_scenario <- function(name, extra = list()) {
  sc <- scenarios[[name]]
  est <- establish(sc$kind)
  wt <- if (!is.null(est$p$WT)) est$p$WT else 12
  m <- mod %>% param(est$p)
  if (length(extra)) m <- m %>% param(extra)
  m <- m %>% init_from(est$init)
  d <- sc$dose(wt)
  if (is.null(sc$quit)) {
    out <- m %>% (function(x) if (is.null(d)) x else x %>% ev(d)) %>%
      mrgsim(end = sc$end, delta = 0.5, recsort = 3) %>% as.data.frame()
  } else {
    # smoking cessation: two parameter segments spliced at `quit`
    o1 <- m %>% (function(x) if (is.null(d)) x else x %>% ev(d)) %>%
      mrgsim(end = sc$quit, delta = 0.5, recsort = 3) %>% as.data.frame()
    m2 <- mod %>% param(est$p)
    if (length(extra)) m2 <- m2 %>% param(extra)
    m2 <- m2 %>% param(SMOKE = 0) %>% init_from(o1[nrow(o1), ])
    o2 <- m2 %>% mrgsim(end = sc$end - sc$quit, delta = 0.5,
                        recsort = 3) %>% as.data.frame()
    o2$time <- o2$time + sc$quit
    out <- rbind(o1, o2[-1, ])
  }
  out$scenario <- name
  out$label <- sc$label
  out
}

run_all <- function() bind_rows(lapply(names(scenarios), run_scenario))

# =============================================================================
# CALIBRATION NOTES — where each number came from
# -----------------------------------------------------------------------------
# PK
#   Vemurafenib CL/F 32 L/day, V 90+25 L, t1/2 ~2.5 d, >99% bound, mild CYP3A4
#     autoinduction: Zhang 2017 (PMID 28255850), Grippo 2014 (PMID 24178368).
#     On 20 mg/kg/day in a 12 kg child the model gives Css ~22 mg/L, inside the
#     range reported in the paediatric vemurafenib series (PMID 31513482).
#   Dabrafenib CL/F 840 L/day, t1/2 ~8 h, active hydroxy metabolite at ~40%
#     molar formation: Ouellet 2014 (PMID 24408395), Falchook 2012
#     (PMID 22608338); combination exposure Balakirouchenane 2020
#     (PMID 32283865).
#   Trametinib CL/F 118 L/day with a deep peripheral compartment giving a
#     terminal t1/2 of ~4 days and a trough of ~9 ng/mL on 0.032 mg/kg/day:
#     Ouellet 2016 (PMID 26940938), Infante 2012 (PMID 22805291).
#   Cladribine CL 150 L/day, CSF ~25% of plasma, intracellular Cd-ATP t1/2
#     ~24 h: Liliemark 1997 (PMID 9068927), Albertioni 1998 (PMID 9533533).
#   Cytarabine CL 1915 L/day (t1/2 ~20 min), ara-CTP t1/2 ~3.7 h, CSF ~15%:
#     Slevin 1983 (PMID 6583325), Heinemann 1988 (PMID 3383195).
#   Vinblastine 3-compartment behaviour collapsed to 2 with t1/2z ~30 h:
#     Balis 1983 (PMID 6189661).
#   Prednisolone CL/F 190 L/day, V 40 L: Petersen 2003 (PMID 12698270).
#   Myelosuppression MTT ~100 h, gamma 0.16, ANC0 4.0: Friberg 2002
#     (PMID 12488418); the grade-4 nadir on 2-CdA/Ara-C matches Donadieu 2015.
#
# PD / disease
#   IC50s are set so that vemurafenib monotherapy suppresses pERK to ~19% of
#     baseline, i.e. just past the >80% blockade that Bollag 2010
#     (PMID 20823850) showed is required for response - and the Hill
#     coefficients (HC = 2 on cyclin D1) put a THRESHOLD there rather than a
#     gradient.
#   PFMAX = 0.09 is the lesional Ki-67 of Brabencova 1998 (PMID 9588881).
#   KPROL = 0.022 < KDL + KIMM*(1-TREGMAX) = 0.0288: lesions cannot sustain
#     themselves locally. This is the quantitative content of that Ki-67.
#   SASP with KSASP = 3/day reproduces resolution of fever/rash/pain within
#     days of starting a MAPK inhibitor: Bigenwald 2021 (PMID 33958797) for the
#     senescence-secretome mechanism, Donadieu 2019 for the clinical tempo.
#   FN_APO = 0.25 and PNICHE = 0.003 are calibrated so that 12 months of MAPK
#     inhibition leaves cfDNA above the ddPCR LOD with a viable reservoir
#     (Eckstein 2019, PMID 30718231; Evseev 2021, PMID 34383272) while six
#     cycles of 2-CdA/Ara-C do not (Donadieu 2015, PMID 26194764).
#   KAVP / AVP_CRIT are set so an untreated CNS-risk patient crosses into
#     permanent central DI over ~6 months, matching the risk-factor structure
#     of Grois 2006 (PMID 16047354) and the permanent-consequence rates of
#     Haupt 2004 (PMID 15049016).
#   RANKL0, KR_RANKL from Makras 2012 (PMID 22278426).
#   ASMK / KCYST from the natural history of adult pulmonary LCH and the
#     nodule-resolves / cyst-persists asymmetry: Tazi 2012 (PMID 22441752),
#     Benattia 2022 (PMID 34675043).
#
# WHAT THE MODEL DOES NOT CLAIM
#   At these parameters the cytotoxic arm is NOT slower to control disease than
#   the MAPK inhibitor, so there is no crossover on the CNS endpoint when both
#   start at presentation. The crossover appears only once the cytotoxic arm is
#   delayed (~15 days in the twin's scan). The model computes that boundary; it
#   does not assume it.
# =============================================================================

# =============================================================================
# EXAMPLE USE
# =============================================================================
if (interactive()) {

  ## --- Commitment 2: response tempo vs mass vs cfDNA -----------------------
  s5 <- run_scenario("S5_vem_continuous")
  s5 %>%
    select(time, DAS, LTOT_units, CFDNA, pERK_pct, SASP) %>%
    filter(time <= 120) %>%
    pivot_longer(-time) %>%
    ggplot(aes(time, value)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~name, scales = "free_y") +
    labs(title = "Vemurafenib: pERK and secretome collapse in days, mass and cfDNA in weeks",
         x = "Days since diagnosis") +
    theme_bw()

  ## --- Commitment 2: rebound on withdrawal --------------------------------
  bind_rows(run_scenario("S5_vem_continuous"),
            run_scenario("S6_vem_stop"),
            run_scenario("S8_bridge_consolidate")) %>%
    ggplot(aes(time, DAS, colour = scenario)) +
    geom_line(linewidth = 0.7) +
    geom_vline(xintercept = 365, linetype = 2) +
    labs(title = "Stopping a cytostatic drug rebounds; consolidating with a cytotoxic one does not",
         x = "Days since diagnosis", y = "Disease Activity Score") +
    theme_bw()

  ## --- Commitment 2 kill switch: rebound must disappear -------------------
  bind_rows(
    run_scenario("S6_vem_stop") %>% mutate(arm = "SL_MAPKI_KILL = 0"),
    run_scenario("S6_vem_stop", list(SL_MAPKI_KILL = 0.25)) %>%
      mutate(arm = "SL_MAPKI_KILL = 0.25 (falsifier)")
  ) %>%
    ggplot(aes(time, CFDNA, colour = arm)) + geom_line() +
    labs(title = "Kill switch: giving the MAPK inhibitor a kill term abolishes the rebound") +
    theme_bw()

  ## --- Commitment 3: sequelae are integrals -------------------------------
  bind_rows(run_scenario("S9_delayed_dx"), run_scenario("S9b_early_dx")) %>%
    select(time, scenario, AVPN, NEUR, TTET, DAS) %>%
    pivot_longer(c(AVPN, NEUR, TTET, DAS)) %>%
    ggplot(aes(time, value, colour = scenario)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~name, scales = "free_y") +
    labs(title = "A six-month diagnostic delay costs AVP neurons that no therapy returns",
         x = "Days since diagnosis") +
    theme_bw()

  ## --- Commitment 1: partition alone flips the phenotype ------------------
  swap <- pheno$MSROpos
  bind_rows(
    run_scenario("S1_observation") %>% mutate(arm = "SS-b origin descriptors"),
    run_scenario("S1_observation", swap) %>%
      mutate(arm = "same rates, MS RO+ origin descriptors")
  ) %>%
    ggplot(aes(time, LRO, colour = arm)) + geom_line(linewidth = 0.7) +
    labs(title = "Only the cell-of-origin descriptors differ",
         y = "Risk-organ burden") +
    theme_bw()

  ## --- All scenarios, key endpoints at end of follow-up -------------------
  run_all() %>%
    group_by(scenario, label) %>%
    slice_tail(n = 1) %>%
    select(scenario, label, DAS, CFDNA, AVPN, NEUR, BILF, LUNGC, BVOL,
           CDI, ND, TTET) %>%
    as.data.frame() %>%
    print()
}
