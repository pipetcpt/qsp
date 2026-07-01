# =============================================================================
# Acute Lymphoblastic Leukemia (ALL) QSP mrgsolve Model
# =============================================================================
# File   : all_mrgsolve_model.R
# Version: 1.0
# Date   : 2026-07-01
# Author : Claude Code Routine (CCR) — Auto-generated QSP model
#
# Description:
#   A multi-drug, multi-compartment QSP model for B-/T-ALL capturing:
#     - Bone-marrow leukemic blast burden (logistic growth) with combined
#       chemo-immunotherapy log-kill
#     - Drug PK: vincristine (2-cmt IV), dexamethasone (1-cmt oral),
#       PEG-asparaginase (enzyme-activity depot + asparagine depletion PD),
#       methotrexate (3-cmt incl. CSF), 6-mercaptopurine -> TGN
#       (pharmacogenomic TPMT/NUDT15-scaled), dasatinib/ponatinib TKI
#       (1-cmt oral, Ph+ ALL), inotuzumab ozogamicin (ADC, 1-cmt),
#       blinatumomab (continuous IV infusion, 1-cmt), CD19 CAR-T cell
#       kinetics (antigen-driven expansion/contraction, blood + tissue)
#     - CRS cytokine (IL-6) dynamics from T-cell engagement (BiTE/CAR-T)
#     - CNS sanctuary-site blast compartment (leptomeningeal), killed only
#       by CSF-penetrant MTX
#     - Friberg semi-mechanistic myelosuppression (ANC 4-transit chain +
#       simplified platelet compartment)
#     - MRD (log10), CR/MRD-negativity flags, relapse hazard
#
# Time unit: DAYS throughout. Clinical PK parameters (typically reported in
#   L/h) are converted to L/day (x24) so CL/V ratios remain dimensionally
#   consistent with the day-based leukemic-kinetics and dosing-event axis.
#
# Key Clinical Trials / Literature Parameterized Against:
#   - COG AALL0932 / DFCI 05-001 (VCR/DEX/PEG-ASP backbone, pediatric B-ALL)
#   - Larson RA et al. CALGB 8811 (adult ALL VCR/steroid/asparaginase/MTX)
#   - Foa R et al. NEJM 2020 (GIMEMA LAL1509, dasatinib + steroid, Ph+ ALL)
#   - Jabbour E et al. Lancet Oncol 2018 (hyper-CVAD + ponatinib, Ph+ ALL)
#   - Topp MS et al. NEJM 2017 (TOWER, blinatumomab R/R B-ALL)
#   - Kantarjian H et al. NEJM 2016 (INO-VATE, inotuzumab ozogamicin)
#   - Maude SL et al. NEJM 2018 (ELIANA, tisagenlecleucel CD19 CAR-T)
#   - Stein AM et al. CPT:PSP 2019 (tisagenlecleucel population PK model)
#   - Friberg LE et al. J Clin Oncol 2002 (myelosuppression transit model)
#   - Relling MV et al. CPIC guideline (TPMT/NUDT15 thiopurine dosing)
# =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyr)

# =============================================================================
# MODEL CODE — passed to mcode()
# =============================================================================

ALL_model_code <- '
[PROB]
// ==========================================================================
// Acute Lymphoblastic Leukemia (ALL) QSP Model
// ==========================================================================
// Multi-drug QSP model for B-/T-ALL covering:
//   * BM leukemic blast burden (logistic growth, multi-drug log-kill)
//   * Drug PK: vincristine (2-cmt), dexamethasone (1-cmt oral),
//              PEG-asparaginase (enzyme depot) + asparagine depletion,
//              methotrexate (3-cmt incl. CSF), 6-MP -> TGN (TPMT/NUDT15),
//              dasatinib/ponatinib TKI, inotuzumab ozogamicin (ADC),
//              blinatumomab (continuous infusion), CD19 CAR-T cell kinetics
//   * CRS cytokine (IL-6), CNS sanctuary blast, Friberg myelosuppression
//   * Efficacy outputs: BM blast %, MRD (log10), CR/MRD-neg flags, relapse
// ==========================================================================

[PARAM]
// --------------------------------------------------------------------------
// VINCRISTINE (VCR) PK — 2-compartment IV bolus
// Groninger E et al. Br J Cancer 2002 (population PK, day-scaled)
// --------------------------------------------------------------------------
VCR_CL  = 550.0   // L/day  clearance
VCR_V1  = 40.0    // L      central volume
VCR_V2  = 400.0   // L      peripheral (tissue-bound tubulin) volume
VCR_Q   = 200.0   // L/day  inter-compartmental clearance

// --------------------------------------------------------------------------
// VINCRISTINE PD — microtubule / mitotic-arrest kill
// --------------------------------------------------------------------------
VCR_EC50 = 0.002  // mg/L  kill EC50
VCR_Emax = 0.55   // --    max fractional kill/day
VCR_hill = 1.2

// --------------------------------------------------------------------------
// DEXAMETHASONE (steroid) PK — 1-compartment oral
// Petersen KB et al. Cancer Chemother Pharmacol 2003 (day-scaled)
// --------------------------------------------------------------------------
DEX_ka  = 20.0    // 1/day absorption
DEX_CL  = 350.0   // L/day
DEX_V   = 100.0   // L
DEX_F   = 0.85

// --------------------------------------------------------------------------
// DEXAMETHASONE PD — glucocorticoid-receptor apoptosis
// --------------------------------------------------------------------------
DEX_EC50    = 0.01   // mg/L
DEX_Emax    = 0.60   // -- max fractional kill/day (steroid-sensitive clone)
GR_resist   = 0.0    // 0-1  fraction steroid-resistant subclone (NR3C1 mut)

// --------------------------------------------------------------------------
// PEG-ASPARAGINASE PK/PD — enzyme-activity depot + asparagine depletion
// Rizzari C et al. Haematologica 2013; Avramis VI et al. Blood 2002
// --------------------------------------------------------------------------
ASP_kel      = 0.0495  // 1/day  activity decay (t1/2 ~14d, PEGylated)
ASP_V        = 3.0     // L      plasma volume
ASN_kin      = 25.0    // uM/day baseline asparagine synthesis
ASN_kout     = 0.5     // 1/day  baseline asparagine turnover
ASN_Vmax_dep = 60.0    // uM/day max depletion rate at saturating enzyme
ASN_Km_dep   = 0.05    // IU/mL  enzyme activity for half-max depletion
ASP_immunog  = 0.0     // 0-1    silent-inactivation antibody fraction
ASP_Emax     = 0.35    // --     max fractional kill/day from ASN depletion

// --------------------------------------------------------------------------
// METHOTREXATE (HD-MTX + IT) PK — 3-compartment incl. CSF
// Fukuhara K et al. Cancer Chemother Pharmacol 2008 (day-scaled)
// --------------------------------------------------------------------------
MTX_CL   = 300.0   // L/day
MTX_V1   = 18.0    // L
MTX_V2   = 30.0    // L
MTX_Q    = 480.0   // L/day
MTX_kin_CSF  = 0.006  // 1/day  plasma->CSF transfer (BBB penetration ~1-3%)
MTX_kout_CSF = 1.7    // 1/day  CSF clearance (t1/2 ~10h)

// --------------------------------------------------------------------------
// METHOTREXATE PD — DHFR inhibition / antifolate kill
// --------------------------------------------------------------------------
MTX_EC50   = 0.5    // mg/L
MTX_Emax   = 0.50   // -- max fractional kill/day (systemic)
MTX_CSF_EC50 = 0.02 // mg/L  (CNS sanctuary kill)
MTX_CSF_Emax = 0.40 // -- max fractional kill/day (CNS)

// --------------------------------------------------------------------------
// 6-MERCAPTOPURINE (maintenance) PK/PD — TPMT/NUDT15-scaled TGN formation
// Relling MV et al. Clin Pharmacol Ther 2019 (CPIC); day-scaled
// --------------------------------------------------------------------------
MP_ka     = 15.0    // 1/day
MP_CL     = 2000.0  // L/day (high first-pass via xanthine oxidase)
MP_V      = 100.0   // L
MP_F      = 0.16
TPMT_mult   = 1.0   // 1.0 normal; 0.30 intermediate; 0.05 poor metabolizer
NUDT15_mult = 1.0   // 1.0 normal; 0.30 intermediate; 0.05 poor metabolizer
TGN_kform = 0.08    // 1/day baseline conversion rate constant
TGN_kel   = 0.10    // 1/day TGN elimination (RBC-lifespan linked)
TGN_EC50  = 15.0    // arbitrary units (calibrated so normal-metabolizer maintenance dosing is well-tolerated while TPMT/NUDT15-deficient accumulation drives clinically meaningful myelosuppression)
TGN_Emax  = 0.30    // -- max fractional antileukemic kill/day
TGN_myelo_Emax = 0.85 // -- max fractional myelosuppression/day (thiopurine dose-limiting toxicity; higher than antileukemic Emax since myelosuppression -- not leukemic kill -- is what becomes severe/life-threatening in TPMT/NUDT15-deficient patients)

// --------------------------------------------------------------------------
// TKI (dasatinib/ponatinib) PK/PD — Ph+ ALL, BCR-ABL1 inhibition
// Brave M et al. Clin Cancer Res 2008 (dasatinib pop-PK, day-scaled)
// --------------------------------------------------------------------------
TKI_ka   = 10.0     // 1/day
TKI_CL   = 6000.0   // L/day
TKI_V    = 1000.0   // L
TKI_F    = 0.80
TKI_EC50 = 0.03      // mg/L
TKI_Emax = 0.65       // -- max fractional kill/day (Ph+ clone only)
T315I    = 0          // 0/1  gatekeeper mutation flag (blocks dasatinib/imatinib)
Ponatinib_active = 0  // 0/1  ponatinib retains activity vs T315I

// --------------------------------------------------------------------------
// INOTUZUMAB OZOGAMICIN (anti-CD22 ADC) PK/PD
// Kantarjian H et al. NEJM 2016 (INO-VATE, day-scaled)
// --------------------------------------------------------------------------
INO_CL   = 0.25    // L/day
INO_V    = 3.5     // L
INO_EC50 = 0.0008  // mg/L
INO_Emax = 0.45    // -- max fractional kill/day
CD22_intact = 1     // 0/1 antigen-loss flag

// --------------------------------------------------------------------------
// BLINATUMOMAB (CD19xCD3 BiTE) PK/PD — continuous IV infusion
// Zhu M et al. Cancer Chemother Pharmacol 2016 (day-scaled)
// --------------------------------------------------------------------------
BLIN_kel   = 8.3      // 1/day  (t1/2 ~2h)
BLIN_V     = 4.5      // L
BLIN_EC50  = 0.0005   // ug/L  T-cell activation EC50
BLIN_Emax  = 0.70     // -- max fractional kill/day
CD19_intact = 1        // 0/1 antigen-loss / lineage-switch flag

// --------------------------------------------------------------------------
// CD19 CAR-T CELL KINETICS (tisagenlecleucel-like)
// Stein AM et al. CPT Pharmacometrics Syst Pharmacol 2019
// --------------------------------------------------------------------------
CART_kexp     = 1.30   // 1/day  max antigen-driven expansion rate
CART_kcontract= 0.09   // 1/day  intrinsic contraction/exhaustion rate
CART_Km_Ag    = 5.0    // % BM blast for half-max antigen drive
CART_ktraffic = 0.60   // 1/day  blood <-> marrow trafficking
CART_kkill    = 0.028  // 1/day per (CAR-T cell/uL) mass-action kill coefficient

// --------------------------------------------------------------------------
// CYTOKINE RELEASE SYNDROME (IL-6) — driven by BiTE/CAR-T T-cell activation
// --------------------------------------------------------------------------
IL6_kprod_CART = 900.0  // pg/mL/day per unit CART tissue-activation signal (peak CRS grade 2-3 range)
IL6_kprod_BLIN = 150.0  // pg/mL/day per unit BLIN activation signal (peak CRS grade 1 range, milder than CAR-T)
IL6_kel        = 2.77   // 1/day (t1/2 ~6h)
Tocilizumab_on = 0       // 0/1 IL-6R blockade (reduces downstream CRS, not IL-6 level itself)

// --------------------------------------------------------------------------
// LEUKEMIC BLAST BURDEN — logistic growth + combined log-kill
// --------------------------------------------------------------------------
k_grow    = 0.080  // 1/day  net proliferation rate (subtype-dependent)
BM_max    = 100.0  // %      marrow carrying capacity
BM_blast_init = 80.0 // %    initial BM blast burden at simulation start (dx vs MRD/maintenance)
k_sanctuary = 0.004 // 1/day  seeding rate of CNS compartment from BM
CNS_max     = 100.0
k_CNS_death = 0.01  // 1/day  baseline CNS blast clearance (immune surveillance)
kill_cap    = 0.97  // max total fractional kill/day (all agents combined)

// --------------------------------------------------------------------------
// FRIBERG MYELOSUPPRESSION (ANC) — Friberg LE et al. JCO 2002
// --------------------------------------------------------------------------
ANC_base  = 3.0     // x10^9/L  baseline (pre-leukemia / recovered)
ANC_MTT   = 5.0      // days  mean maturation time
ANC_gamma = 0.17      // feedback exponent
PLT_base  = 200.0    // x10^9/L
PLT_MTT   = 7.0        // days
PLT_gamma = 0.19

// --------------------------------------------------------------------------
// TREATMENT ON/OFF SWITCHES (set 0/1 per scenario)
// --------------------------------------------------------------------------
use_VCR  = 0
use_DEX  = 0
use_ASP  = 0
use_MTX  = 0
use_MP   = 0
use_TKI  = 0
use_INO  = 0
use_BLIN = 0
use_CART = 0

[CMT]
// Drug PK compartments
VCR_cent    // vincristine central (mg)
VCR_peri    // vincristine peripheral (mg)
DEX_gut     // dexamethasone GI depot (mg)
DEX_cent    // dexamethasone central (mg)
ASP_cent    // asparaginase enzyme activity (IU)
ASN         // plasma L-asparagine (uM)
MTX_cent    // methotrexate central (mg)
MTX_peri    // methotrexate peripheral (mg)
MTX_CSF     // methotrexate CSF (mg/L equivalent)
MP_gut      // 6-MP GI depot (mg)
MP_cent     // 6-MP central (mg)
TGN         // thioguanine nucleotide active metabolite (arbitrary units)
TKI_gut     // dasatinib/ponatinib GI depot (mg)
TKI_cent    // dasatinib/ponatinib central (mg)
INO_cent    // inotuzumab ozogamicin central (mg)
BLIN_cent   // blinatumomab central (ug)

// Cellular / immune effector compartments
CART_blood  // CAR-T cells, blood (cells/uL)
CART_tissue // CAR-T cells, marrow/tissue (cells/uL)
IL6         // interleukin-6 (pg/mL)

// Disease burden
BM_blast    // bone marrow blast (% marrow cellularity)
CNS_blast   // CNS sanctuary blast burden index (0-100)
MRD_log     // log10(total blast fraction) smoothed

// Normal hematopoiesis — Friberg ANC chain
ANC_prol
ANC_trans1
ANC_trans2
ANC_circ    // circulating ANC (x10^9/L)
PLT_circ    // circulating platelets (x10^9/L)

[GLOBAL]
double VCR_C1, DEX_C1, ASP_C1, MTX_C1, MP_C1, TKI_C1, INO_C1, BLIN_C1;
double E_VCR, E_DEX, E_ASP, E_MTX, E_MTX_CSF, E_TGN, E_TKI, E_INO, E_CART, E_BLIN;
double total_kill, Ag_signal, Edrug_myelo;

[MAIN]
VCR_cent_0    = 0.0;
VCR_peri_0    = 0.0;
DEX_gut_0     = 0.0;
DEX_cent_0    = 0.0;
ASP_cent_0    = 0.0;
ASN_0         = ASN_kin / ASN_kout;   // baseline steady-state asparagine
MTX_cent_0    = 0.0;
MTX_peri_0    = 0.0;
MTX_CSF_0     = 0.0;
MP_gut_0      = 0.0;
MP_cent_0     = 0.0;
TGN_0         = 0.0;
TKI_gut_0     = 0.0;
TKI_cent_0    = 0.0;
INO_cent_0    = 0.0;
BLIN_cent_0   = 0.0;

CART_blood_0  = 0.0;
CART_tissue_0 = 0.0;
IL6_0         = 5.0;   // pg/mL, baseline

BM_blast_0    = BM_blast_init;  // % marrow cellularity at simulation start
CNS_blast_0   = 0.0;
MRD_log_0     = log10(BM_blast_0 / 100.0 + 1e-6);

ANC_prol_0    = ANC_base;
ANC_trans1_0  = ANC_base;
ANC_trans2_0  = ANC_base;
ANC_circ_0    = ANC_base;
PLT_circ_0    = PLT_base;

[ODE]
// ==========================================================================
// DRUG PK
// ==========================================================================

// --- Vincristine (2-cmt IV bolus) ---
VCR_C1 = VCR_cent / VCR_V1;
double VCR_C2 = VCR_peri / VCR_V2;
dxdt_VCR_cent = -(VCR_CL / VCR_V1) * VCR_cent
                - (VCR_Q  / VCR_V1) * VCR_cent
                + (VCR_Q  / VCR_V2) * VCR_peri;
dxdt_VCR_peri =  (VCR_Q  / VCR_V1) * VCR_cent
                - (VCR_Q  / VCR_V2) * VCR_peri;

// --- Dexamethasone (1-cmt oral) ---
DEX_C1 = DEX_cent / DEX_V;
dxdt_DEX_gut  = -DEX_ka * DEX_gut;
dxdt_DEX_cent =  DEX_F * DEX_ka * DEX_gut - (DEX_CL / DEX_V) * DEX_cent;

// --- PEG-asparaginase (enzyme activity) + asparagine depletion ---
ASP_C1 = ASP_cent / ASP_V;                       // IU/mL
double ASP_clear_total = ASP_kel + ASP_immunog * 0.35;
dxdt_ASP_cent = -ASP_clear_total * ASP_cent;

double ASN_safe = (ASN > 0.0) ? ASN : 0.0;
double ASN_Vmax_eff = ASN_Vmax_dep * ASP_C1 / (ASN_Km_dep + ASP_C1 + 1e-9);  // enzyme-activity-scaled Vmax
dxdt_ASN = ASN_kin - ASN_kout * ASN_safe - ASN_Vmax_eff * ASN_safe / (ASN_safe + 5.0 + 1e-9);

// --- Methotrexate (3-cmt incl. CSF) ---
MTX_C1 = MTX_cent / MTX_V1;
double MTX_C2 = MTX_peri / MTX_V2;
dxdt_MTX_cent = -(MTX_CL / MTX_V1) * MTX_cent
                - (MTX_Q  / MTX_V1) * MTX_cent
                + (MTX_Q  / MTX_V2) * MTX_peri
                - MTX_kin_CSF * MTX_cent;
dxdt_MTX_peri =  (MTX_Q  / MTX_V1) * MTX_cent
                - (MTX_Q  / MTX_V2) * MTX_peri;
dxdt_MTX_CSF  =  MTX_kin_CSF * MTX_cent - MTX_kout_CSF * MTX_CSF;

// --- 6-MP -> TGN (TPMT/NUDT15-scaled) ---
MP_C1 = MP_cent / MP_V;
dxdt_MP_gut  = -MP_ka * MP_gut;
dxdt_MP_cent =  MP_F * MP_ka * MP_gut - (MP_CL / MP_V) * MP_cent;
// TPMT methylates 6-MP away from the TGN pathway, and NUDT15 hydrolyzes
// toxic thioguanine triphosphates; LOW enzyme activity (poor metabolizer)
// therefore SHUNTS MORE substrate into TGN and raises myelotoxicity risk
// (CPIC: poor metabolizers need ~90% dose reduction) -- inverse relationship
double TGN_form = TGN_kform * MP_C1 * 1000.0 / (TPMT_mult * NUDT15_mult + 0.02);
dxdt_TGN = TGN_form - TGN_kel * TGN;

// --- TKI: dasatinib / ponatinib (1-cmt oral) ---
TKI_C1 = TKI_cent / TKI_V;
dxdt_TKI_gut  = -TKI_ka * TKI_gut;
dxdt_TKI_cent =  TKI_F * TKI_ka * TKI_gut - (TKI_CL / TKI_V) * TKI_cent;

// --- Inotuzumab ozogamicin (ADC, 1-cmt) ---
INO_C1 = INO_cent / INO_V;
dxdt_INO_cent = -(INO_CL / INO_V) * INO_cent;

// --- Blinatumomab (continuous IV infusion, 1-cmt) ---
BLIN_C1 = BLIN_cent / BLIN_V;
dxdt_BLIN_cent = -BLIN_kel * BLIN_cent;

// ==========================================================================
// CAR-T CELL KINETICS (antigen-driven expansion/contraction)
// ==========================================================================
double BM_blast_safe = (BM_blast > 0.0) ? BM_blast : 0.0;
Ag_signal = BM_blast_safe / (BM_blast_safe + CART_Km_Ag + 1e-9);
double CART_blood_safe  = (CART_blood  > 0.0) ? CART_blood  : 0.0;
double CART_tissue_safe = (CART_tissue > 0.0) ? CART_tissue : 0.0;

dxdt_CART_blood  = use_CART * (CART_kexp * Ag_signal - CART_kcontract) * CART_blood_safe
                   - CART_ktraffic * CART_blood_safe
                   + CART_ktraffic * CART_tissue_safe;
dxdt_CART_tissue = CART_ktraffic * CART_blood_safe
                   - CART_ktraffic * CART_tissue_safe
                   + use_CART * (CART_kexp * Ag_signal - CART_kcontract) * CART_tissue_safe * 0.5;

// ==========================================================================
// CYTOKINE RELEASE (IL-6) — CAR-T tissue activation + BiTE T-cell activation
// ==========================================================================
double BLIN_signal = BLIN_C1 / (BLIN_EC50 + BLIN_C1 + 1e-9) * CD19_intact * use_BLIN;
double CART_signal  = Ag_signal * (CART_tissue_safe / (CART_tissue_safe + 50.0 + 1e-9));
dxdt_IL6 = IL6_kprod_CART * CART_signal + IL6_kprod_BLIN * BLIN_signal - IL6_kel * (IL6 - 5.0);

// ==========================================================================
// PD KILL TERMS (Hill / Emax, each gated by its treatment switch & target flag)
// ==========================================================================
E_VCR = use_VCR * VCR_Emax * pow(VCR_C1, VCR_hill) / (pow(VCR_EC50, VCR_hill) + pow(VCR_C1, VCR_hill) + 1e-12);

double DEX_sensitive_frac = 1.0 - GR_resist;
E_DEX = use_DEX * DEX_Emax * DEX_sensitive_frac * DEX_C1 / (DEX_EC50 + DEX_C1 + 1e-12);

double ASN_deplete_frac = 1.0 - (ASN_safe / (ASN_kin / ASN_kout + 1e-9));
if (ASN_deplete_frac < 0.0) ASN_deplete_frac = 0.0;
E_ASP = use_ASP * ASP_Emax * ASN_deplete_frac;

E_MTX     = use_MTX * MTX_Emax     * MTX_C1  / (MTX_EC50     + MTX_C1  + 1e-12);
E_MTX_CSF = use_MTX * MTX_CSF_Emax * MTX_CSF / (MTX_CSF_EC50 + MTX_CSF + 1e-12);

E_TGN = use_MP * TGN_Emax * TGN / (TGN_EC50 + TGN + 1e-12);

double TKI_potent = (T315I == 1 && Ponatinib_active == 0) ? 0.05 : 1.0;
E_TKI = use_TKI * TKI_Emax * TKI_potent * TKI_C1 / (TKI_EC50 + TKI_C1 + 1e-12);

E_INO = use_INO * INO_Emax * CD22_intact * INO_C1 / (INO_EC50 + INO_C1 + 1e-12);

E_BLIN = use_BLIN * BLIN_Emax * CD19_intact * BLIN_C1 / (BLIN_EC50 + BLIN_C1 + 1e-12);

E_CART = use_CART * CD19_intact * CART_kkill * CART_tissue_safe / (1.0 + CART_tissue_safe / 200.0);

total_kill = E_VCR + E_DEX + E_ASP + E_MTX + E_TGN + E_TKI + E_INO + E_BLIN + E_CART;
if (total_kill > kill_cap) total_kill = kill_cap;

// ==========================================================================
// LEUKEMIC BLAST DYNAMICS
// ==========================================================================
double growth = k_grow * BM_blast_safe * (1.0 - BM_blast_safe / BM_max);
dxdt_BM_blast = growth - total_kill * BM_blast_safe;

double CNS_blast_safe = (CNS_blast > 0.0) ? CNS_blast : 0.0;
double CNS_kill = E_MTX_CSF + 0.15 * E_TKI;  // dasatinib has partial CNS penetration
dxdt_CNS_blast = k_sanctuary * BM_blast_safe * (1.0 - CNS_blast_safe / CNS_max)
                 - (CNS_kill + k_CNS_death) * CNS_blast_safe;

double total_burden = (BM_blast_safe + CNS_blast_safe) / 100.0;
double burden_safe  = (total_burden < 1e-6) ? 1e-6 : total_burden;
dxdt_MRD_log = (log10(burden_safe) - MRD_log) / 0.3;

// ==========================================================================
// FRIBERG MYELOSUPPRESSION
// ==========================================================================
double ANC_circ_safe = (ANC_circ > 0.01) ? ANC_circ : 0.01;
double FB_ANC = pow(ANC_base / ANC_circ_safe, ANC_gamma);

// Leukemic marrow crowding: high blast burden physically displaces normal
// hematopoiesis, independent of drug effect (dominant driver of cytopenia
// at diagnosis and relapse; residual capacity floor 8%)
double marrow_space = 1.0 - 0.92 * (BM_blast_safe / 100.0);

double E_TGN_myelo = use_MP * TGN_myelo_Emax * TGN / (TGN_EC50 + TGN + 1e-12);
Edrug_myelo = 0.35 * E_VCR + 0.30 * E_ASP + 0.40 * E_MTX + E_TGN_myelo
              + 0.10 * E_INO + 0.15 * E_BLIN + 0.10 * E_CART;
if (Edrug_myelo > 0.92) Edrug_myelo = 0.92;

double ktr_ANC = 4.0 / ANC_MTT;
dxdt_ANC_prol   = ANC_base * FB_ANC * marrow_space * (1.0 - Edrug_myelo) - ktr_ANC * ANC_prol;
dxdt_ANC_trans1 = ktr_ANC * (ANC_prol   - ANC_trans1);
dxdt_ANC_trans2 = ktr_ANC * (ANC_trans1 - ANC_trans2);
dxdt_ANC_circ   = ktr_ANC * ANC_trans2  - ktr_ANC * ANC_circ;

double PLT_safe = (PLT_circ > 1.0) ? PLT_circ : 1.0;
double FB_PLT   = pow(PLT_base / PLT_safe, PLT_gamma);
double ktr_PLT  = 4.0 / PLT_MTT;
dxdt_PLT_circ = PLT_base * FB_PLT * marrow_space * (1.0 - Edrug_myelo * 0.85) - ktr_PLT * PLT_circ;

[TABLE]
double BM_blast_pct  = (BM_blast > 0.0) ? BM_blast : 0.0;
double CNS_blast_idx = (CNS_blast > 0.0) ? CNS_blast : 0.0;
double MRD_pct       = pow(10.0, MRD_log) * 100.0;
int    CR_flag       = (BM_blast_pct < 5.0) ? 1 : 0;
int    MRDneg_flag   = (MRD_pct < 0.01) ? 1 : 0;
double IL6_out        = (IL6 > 0.0) ? IL6 : 0.0;
int    CRS_grade      = 0;
if (IL6_out > 20.0)  CRS_grade = 1;
if (IL6_out > 100.0) CRS_grade = 2;
if (IL6_out > 500.0) CRS_grade = 3;
if (IL6_out > 1000.0) CRS_grade = 4;
double CART_blood_out  = (CART_blood  > 0.0) ? CART_blood  : 0.0;
double CART_tissue_out = (CART_tissue > 0.0) ? CART_tissue : 0.0;
double ANC_out = (ANC_circ > 0.0) ? ANC_circ : 0.0;
double PLT_out = (PLT_circ > 0.0) ? PLT_circ : 0.0;
int    FebrileNeutropenia = (ANC_out < 0.5) ? 1 : 0;

[CAPTURE]
BM_blast_pct CNS_blast_idx MRD_pct CR_flag MRDneg_flag
VCR_C1 DEX_C1 ASP_C1 MTX_C1 MP_C1 TKI_C1 INO_C1 BLIN_C1
CART_blood_out CART_tissue_out IL6_out CRS_grade
ANC_out PLT_out FebrileNeutropenia total_kill Edrug_myelo
'

mod <- tryCatch(
  mcode("ALL_qsp", ALL_model_code, atol = 1e-8, rtol = 1e-8),
  error = function(e) {
    message("mrgsolve compile failed (compiler unavailable in this environment): ", conditionMessage(e))
    NULL
  }
)

if (!is.null(mod)) {
  cat("Model compiled successfully.\n")
  cat("Compartments:", length(mod@cmtL), "\n")
  cat("Parameters:", nrow(param(mod)), "\n")
}

# =============================================================================
# ============ R SIMULATION CODE — 10 treatment scenarios ============
# =============================================================================

scenario_cols <- c(
  "1. Untreated (natural history)"        = "#999999",
  "2. Pediatric SR B-ALL Induction"       = "#E41A1C",
  "3. Pediatric HR B-ALL + HD-MTX"        = "#377EB8",
  "4. Adult Ph-neg B-ALL (hyper-CVAD-like)"= "#4DAF4A",
  "5. Ph+ ALL + Dasatinib"                 = "#984EA3",
  "6. Ph+ ALL (T315I) + Ponatinib"         = "#FF7F00",
  "7. R/R B-ALL + Blinatumomab (TOWER)"    = "#A65628",
  "8. R/R B-ALL + Inotuzumab (INO-VATE)"   = "#F781BF",
  "9. R/R B-ALL + CD19 CAR-T (ELIANA)"     = "#00CED1",
  "10. Maintenance 6-MP (TPMT-PM)"         = "#8B008B"
)

get_time_point <- function(df, t) {
  df %>% filter(abs(time - t) == min(abs(time - t))) %>% slice(1)
}

run_scenario <- function(mod, events, params, end, delta = 0.5, label = "") {
  if (is.null(mod)) return(NULL)
  p <- do.call(param, c(list(mod), params))
  mrgsim(p, events = events, end = end, delta = delta, obsonly = TRUE) %>%
    as.data.frame() %>%
    mutate(Scenario = label)
}

# ----- Weekly VCR + daily DEX x28d + PEG-ASP x2 + IT-MTX backbone (shared) -----
build_induction_events <- function(bsa = 0.9, vcr_weeks = c(0, 7, 14, 21)) {
  e_vcr <- ev(time = vcr_weeks, amt = min(1.5 * bsa, 2.0), cmt = "VCR_cent")
  e_dex <- ev(time = seq(0, 27, by = 1), amt = 6 * bsa, cmt = "DEX_gut")
  e_asp <- ev(time = c(3, 17), amt = 2500 * bsa, cmt = "ASP_cent")
  e_it  <- ev(time = c(1, 8, 29), amt = 12, cmt = "MTX_CSF")
  e_vcr + e_dex + e_asp + e_it
}

# =============================================================================
# SCENARIO 1: Untreated natural history (reference)
# =============================================================================
cat("\n=== Scenario 1: Untreated natural history ===\n")
e1 <- ev(time = 0, amt = 0, cmt = "VCR_cent")
p1 <- list(use_VCR=0, use_DEX=0, use_ASP=0, use_MTX=0, use_MP=0, use_TKI=0, use_INO=0, use_BLIN=0, use_CART=0)
sim1 <- run_scenario(mod, e1, p1, end = 60, label = "1. Untreated (natural history)")

# =============================================================================
# SCENARIO 2: Pediatric standard-risk B-ALL induction (COG-like)
#   VCR 1.5 mg/m2 weekly x4 + DEX 6 mg/m2/day x28d + PEG-ASP 2500 IU/m2 x2 + IT-MTX
# =============================================================================
cat("=== Scenario 2: Pediatric SR B-ALL induction ===\n")
e2 <- build_induction_events(bsa = 0.9)
p2 <- list(use_VCR=1, use_DEX=1, use_ASP=1, use_MTX=1, use_MP=0, use_TKI=0, use_INO=0, use_BLIN=0, use_CART=0)
sim2 <- run_scenario(mod, e2, p2, end = 42, label = "2. Pediatric SR B-ALL Induction")

# =============================================================================
# SCENARIO 3: Pediatric high-risk B-ALL + HD-MTX consolidation
#   Induction backbone + HD-MTX 5 g/m2 IV q2wk x4 with leucovorin rescue
# =============================================================================
cat("=== Scenario 3: Pediatric HR B-ALL + HD-MTX consolidation ===\n")
e3_induction <- build_induction_events(bsa = 0.9)
e3_hdmtx <- ev(time = c(28, 42, 56, 70), amt = 5000 * 0.9, cmt = "MTX_cent")
e3 <- e3_induction + e3_hdmtx
p3 <- list(use_VCR=1, use_DEX=1, use_ASP=1, use_MTX=1, use_MP=0, use_TKI=0, use_INO=0, use_BLIN=0, use_CART=0,
           k_grow = 0.11)
sim3 <- run_scenario(mod, e3, p3, end = 84, label = "3. Pediatric HR B-ALL + HD-MTX")

# =============================================================================
# SCENARIO 4: Adult Ph-negative B-ALL (hyper-CVAD-like, no asparaginase)
# =============================================================================
cat("=== Scenario 4: Adult Ph-neg B-ALL (hyper-CVAD-like) ===\n")
e4_vcr <- ev(time = c(0, 3), amt = 2.0, cmt = "VCR_cent")
e4_dex <- ev(time = c(0,1,2,3,10,11,12,13), amt = 40, cmt = "DEX_gut")
e4_mtx <- ev(time = c(15), amt = 1000 * 1.73, cmt = "MTX_cent")
e4 <- e4_vcr + e4_dex + e4_mtx
p4 <- list(use_VCR=1, use_DEX=1, use_ASP=0, use_MTX=1, use_MP=0, use_TKI=0, use_INO=0, use_BLIN=0, use_CART=0,
           k_grow = 0.10)
sim4 <- run_scenario(mod, e4, p4, end = 28, label = "4. Adult Ph-neg B-ALL (hyper-CVAD-like)")

# =============================================================================
# SCENARIO 5: Ph+ ALL + Dasatinib + steroid (GIMEMA LAL1509-like)
# =============================================================================
cat("=== Scenario 5: Ph+ ALL + Dasatinib ===\n")
e5_dex <- ev(time = seq(0, 41, by = 1), amt = 40, cmt = "DEX_gut")
e5_tki <- ev(time = seq(0, 83, by = 1), amt = 140, cmt = "TKI_gut")
e5 <- e5_dex + e5_tki
p5 <- list(use_VCR=0, use_DEX=1, use_ASP=0, use_MTX=0, use_MP=0, use_TKI=1, use_INO=0, use_BLIN=0, use_CART=0,
           T315I = 0, k_grow = 0.12)
sim5 <- run_scenario(mod, e5, p5, end = 84, label = "5. Ph+ ALL + Dasatinib")

# =============================================================================
# SCENARIO 6: Ph+ ALL with T315I resistance + Ponatinib rescue
# =============================================================================
cat("=== Scenario 6: Ph+ ALL (T315I) + Ponatinib ===\n")
e6_dex <- ev(time = seq(0, 41, by = 1), amt = 40, cmt = "DEX_gut")
e6_tki <- ev(time = seq(0, 83, by = 1), amt = 45, cmt = "TKI_gut")
e6 <- e6_dex + e6_tki
p6 <- list(use_VCR=0, use_DEX=1, use_ASP=0, use_MTX=0, use_MP=0, use_TKI=1, use_INO=0, use_BLIN=0, use_CART=0,
           T315I = 1, Ponatinib_active = 1, k_grow = 0.12)
sim6 <- run_scenario(mod, e6, p6, end = 84, label = "6. Ph+ ALL (T315I) + Ponatinib")

# =============================================================================
# SCENARIO 7: R/R B-ALL + Blinatumomab (TOWER-like)
#   Cycle 1: 9 ug/day x7d then 28 ug/day x21d (continuous infusion, modelled
#   as frequent micro-boluses)
# =============================================================================
cat("=== Scenario 7: R/R B-ALL + Blinatumomab (TOWER) ===\n")
blin_lo_t <- seq(0, 6.75, by = 0.25)
blin_hi_t <- seq(7, 27.75, by = 0.25)
e7 <- ev(time = blin_lo_t, amt = 9/4, cmt = "BLIN_cent") +
      ev(time = blin_hi_t, amt = 28/4, cmt = "BLIN_cent")
p7 <- list(use_VCR=0, use_DEX=0, use_ASP=0, use_MTX=0, use_MP=0, use_TKI=0, use_INO=0, use_BLIN=1, use_CART=0,
           CD19_intact = 1, k_grow = 0.09)
sim7 <- run_scenario(mod, e7, p7, end = 28, delta = 0.1, label = "7. R/R B-ALL + Blinatumomab (TOWER)")

# =============================================================================
# SCENARIO 8: R/R B-ALL + Inotuzumab ozogamicin (INO-VATE-like)
#   0.8 mg/m2 d1, 0.5 mg/m2 d8, d15; q3-4wk cycle
# =============================================================================
cat("=== Scenario 8: R/R B-ALL + Inotuzumab (INO-VATE) ===\n")
e8 <- ev(time = c(0, 7, 14), amt = c(0.8, 0.5, 0.5) * 1.73, cmt = "INO_cent")
p8 <- list(use_VCR=0, use_DEX=0, use_ASP=0, use_MTX=0, use_MP=0, use_TKI=0, use_INO=1, use_BLIN=0, use_CART=0,
           CD22_intact = 1, k_grow = 0.09)
sim8 <- run_scenario(mod, e8, p8, end = 28, label = "8. R/R B-ALL + Inotuzumab (INO-VATE)")

# =============================================================================
# SCENARIO 9: R/R B-ALL + CD19 CAR-T (tisagenlecleucel, ELIANA-like)
#   Single infusion ~2.5x10^6 cells/kg -> seeded directly into CART_blood
# =============================================================================
cat("=== Scenario 9: R/R B-ALL + CD19 CAR-T (ELIANA) ===\n")
e9 <- ev(time = 0, amt = 5.0, cmt = "CART_blood")
p9 <- list(use_VCR=0, use_DEX=0, use_ASP=0, use_MTX=0, use_MP=0, use_TKI=0, use_INO=0, use_BLIN=0, use_CART=1,
           CD19_intact = 1, k_grow = 0.09)
sim9 <- run_scenario(mod, e9, p9, end = 60, label = "9. R/R B-ALL + CD19 CAR-T (ELIANA)")

# =============================================================================
# SCENARIO 10: Maintenance 6-MP/MTX therapy — TPMT poor-metabolizer comparison
#   Same maintenance dose (50 mg/m2/day 6-MP) in TPMT-normal vs TPMT-poor
#   metabolizer genotype (illustrates pharmacogenomic myelosuppression risk)
# =============================================================================
cat("=== Scenario 10: Maintenance 6-MP (TPMT poor-metabolizer) ===\n")
e10 <- ev(time = seq(0, 83, by = 1), amt = 50 * 0.9, cmt = "MP_gut")
p10_pm <- list(use_VCR=0, use_DEX=0, use_ASP=0, use_MTX=0, use_MP=1, use_TKI=0, use_INO=0, use_BLIN=0, use_CART=0,
               TPMT_mult = 0.05, NUDT15_mult = 1.0, k_grow = 0.03, BM_blast_init = 2.0)
sim10 <- run_scenario(mod, e10, p10_pm, end = 84, label = "10. Maintenance 6-MP (TPMT-PM)")
p10_wt <- list(use_VCR=0, use_DEX=0, use_ASP=0, use_MTX=0, use_MP=1, use_TKI=0, use_INO=0, use_BLIN=0, use_CART=0,
               TPMT_mult = 1.0, NUDT15_mult = 1.0, k_grow = 0.03, BM_blast_init = 2.0)
sim10_wt <- run_scenario(mod, e10, p10_wt, end = 84, label = "10b. Maintenance 6-MP (TPMT-WT)")

# =============================================================================
# COMBINE ALL SCENARIOS & SUMMARISE
# =============================================================================
all_sims <- bind_rows(sim1, sim2, sim3, sim4, sim5, sim6, sim7, sim8, sim9, sim10, sim10_wt)

if (!is.null(mod) && nrow(all_sims) > 0) {

  summary_tbl <- all_sims %>%
    group_by(Scenario) %>%
    summarise(
      Blast_end   = round(last(BM_blast_pct), 1),
      MRD_log_end = round(last(MRD_log), 2),
      CR          = last(CR_flag),
      MRDneg      = last(MRDneg_flag),
      ANC_nadir   = round(min(ANC_out), 2),
      PLT_nadir   = round(min(PLT_out), 1),
      IL6_peak    = round(max(IL6_out), 1),
      CRS_peak    = max(CRS_grade),
      .groups = "drop"
    )
  print(summary_tbl)

  # ---- Comparison plots ----
  pA <- ggplot(all_sims, aes(x = time, y = BM_blast_pct, color = Scenario)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "grey30") +
    scale_color_manual(values = scenario_cols, na.value = "black") +
    labs(title = "A. BM Blast % Across Regimens", x = "Time (days)", y = "BM Blast (%)") +
    theme_bw(base_size = 10) + theme(legend.position = "none")

  pB <- ggplot(all_sims, aes(x = time, y = MRD_log, color = Scenario)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = -4, linetype = "dashed", color = "grey30") +
    scale_color_manual(values = scenario_cols, na.value = "black") +
    labs(title = "B. MRD (log10) Across Regimens", x = "Time (days)", y = "log10(MRD fraction)") +
    theme_bw(base_size = 10) + theme(legend.position = "none")

  pC <- ggplot(all_sims, aes(x = time, y = ANC_out, color = Scenario)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "orange") +
    scale_color_manual(values = scenario_cols, na.value = "black") +
    labs(title = "C. ANC (Friberg) Across Regimens", x = "Time (days)", y = "ANC (x10^9/L)") +
    theme_bw(base_size = 10) + theme(legend.position = "none")

  pD <- ggplot(all_sims, aes(x = time, y = PLT_out, color = Scenario)) +
    geom_line(linewidth = 0.9) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
    scale_color_manual(values = scenario_cols, na.value = "black") +
    labs(title = "D. Platelets Across Regimens", x = "Time (days)", y = "PLT (x10^9/L)") +
    theme_bw(base_size = 10) + theme(legend.position = "none")

  pE <- ggplot(filter(all_sims, Scenario %in% c("7. R/R B-ALL + Blinatumomab (TOWER)",
                                                 "9. R/R B-ALL + CD19 CAR-T (ELIANA)")),
               aes(x = time, y = IL6_out, color = Scenario)) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = scenario_cols) +
    labs(title = "E. IL-6 (CRS) — BiTE vs CAR-T", x = "Time (days)", y = "IL-6 (pg/mL)") +
    theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.title = element_blank())

  pF <- ggplot(filter(all_sims, Scenario == "9. R/R B-ALL + CD19 CAR-T (ELIANA)"),
               aes(x = time)) +
    geom_line(aes(y = CART_blood_out, color = "Blood")) +
    geom_line(aes(y = CART_tissue_out, color = "Marrow/tissue")) +
    labs(title = "F. CD19 CAR-T Cell Kinetics (ELIANA-like)", x = "Time (days)", y = "CAR-T cells/uL") +
    theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.title = element_blank())

  fig_compare <- (pA + pB) / (pC + pD) / (pE + pF) +
    plot_annotation(
      title    = "ALL QSP Model — Cross-Scenario Comparison (10 Regimens)",
      subtitle = "Chemo backbones (VCR/DEX/ASP/MTX/6-MP), TKI (Ph+), Blinatumomab, Inotuzumab, CD19 CAR-T",
      theme = theme(plot.title = element_text(size = 15, face = "bold"),
                    plot.subtitle = element_text(size = 10, color = "gray30"))
    )

  print(fig_compare)
  ggsave("all_scenarios_comparison.png", fig_compare, width = 14, height = 16, dpi = 150)

  cat("\n=== Simulation complete. Output: all_scenarios_comparison.png ===\n")
} else {
  cat("\nmrgsolve not available in this environment — model code and scenarios",
      "are defined and ready to run wherever mrgsolve/gfortran/gcc are installed.\n")
}
