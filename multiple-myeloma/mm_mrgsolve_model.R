## ============================================================
## Multiple Myeloma (MM) QSP mrgsolve Model
## ============================================================
## Components:
##   - Tumor cell dynamics (sensitive + resistant clones)
##   - M-protein / serum free light chain (sFLC) PD
##   - Bone remodeling (OB/OC coupling + RANKL/OPG axis)
##   - Drug PK: Bortezomib (2-cmt), Carfilzomib (2-cmt),
##              Lenalidomide (1-cmt), Daratumumab (2-cmt + TMDD),
##              Dexamethasone (1-cmt), Venetoclax (1-cmt),
##              Zoledronic acid (2-cmt bone)
##   - Cytokines: IL-6, VEGF (simplified)
##   - Anemia marker: Hemoglobin
## ============================================================
## Key References:
##   Dingli et al. Blood 2007 (myeloma growth model)
##   Lemaire et al. PLoS ONE 2004 (bone remodeling ODE)
##   Chari et al. NEJM 2019 (daratumumab TMDD)
##   Lacy et al. JCO 2010 (lenalidomide PK/PD)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ============================================================
## MODEL DEFINITION
## ============================================================
mm_model_code <- '
$PROB
Multiple Myeloma QSP Model
Bone marrow microenvironment, myeloma cell dynamics,
drug PK/PD (BTZ, CFZ, LEN, DARA, DEX, VEN, ZOL)

$PARAM
// ---- Tumor Cell Biology ----
kprol  = 0.020    // MM cell net growth rate (1/day); baseline doubling ~35d
kapop  = 0.005    // baseline apoptosis rate (1/day)
kgrow  = 0.015    // net proliferation = kprol - kapop
Kmax   = 5e9      // carrying capacity (cells equivalent)
kres   = 1e-4     // de-novo resistance emergence rate (1/day/cell)
kgrow_res = 0.018 // resistant clone growth rate

// ---- M-Protein / sFLC Kinetics ----
kMP_prod = 1.5e-9  // M-protein production per MM cell (g/dL/cell/day)
kMP_elim = 0.08    // M-protein elimination rate (1/day); t½ ~ 9 days
kFLC_prod = 2.0e-10 // sFLC production per MM cell (mg/dL/cell/day)
kFLC_elim = 0.15   // sFLC elimination rate (1/day); t½ ~ 4.6 days

// ---- IL-6 / VEGF (simplified) ----
kIL6_prod = 0.05   // IL-6 production rate (arbitrary units/day)
kIL6_elim = 1.2    // IL-6 turnover (1/day)
kIL6_stim  = 0.3   // IL-6 stimulation of MM cell growth (per unit IL6)
kVEGF_prod = 0.02
kVEGF_elim = 0.8

// ---- Bone Remodeling (Lemaire-type) ----
k_OB_form  = 0.0016  // OB formation rate (cells/day)
k_OB_apop  = 0.0020  // OB apoptosis rate (1/day)
k_OC_form  = 0.0006  // OC formation rate
k_OC_apop  = 0.0025  // OC apoptosis rate
k_RANKL    = 0.005   // RANKL production by MM cells (AU/cell/day) ×1e-9
k_OPG      = 0.002   // OPG production by OB
k_OPG_inh  = 0.3     // OPG suppression by MM cells
BoneV0     = 100.0   // Initial bone volume (%)
k_bone_res = 0.0003  // OC-mediated bone resorption rate
k_bone_form= 0.0002  // OB-mediated bone formation rate
DKK1_base  = 1.0     // DKK1 baseline (AU)
kDKK1_MM   = 0.002   // DKK1 production by MM cells

// ---- Hemoglobin ----
Hgb0       = 14.0   // normal Hgb (g/dL)
kHgb_recover= 0.02  // Hgb recovery rate (1/day)
kHgb_supp  = 0.3    // MM-mediated Hgb suppression (per tumor burden)
IL6_Hgb_inh= 0.5    // IL-6 mediated EPO blunting

// ---- Beta2-Microglobulin ----
kB2M_prod  = 1.0e-10  // β₂M production per MM cell
kB2M_elim  = 0.25     // β₂M clearance (1/day)

// ---- Drug PK: Bortezomib (BTZ, IV/SC) ----
// 3-compartment model based on Moreau et al. Blood 2011
BTZ_CL   = 9.0     // clearance (L/hr)
BTZ_V1   = 4.7     // central volume (L)
BTZ_V2   = 34.5    // peripheral 1 (L)
BTZ_V3   = 298.0   // peripheral 2 (L)
BTZ_Q2   = 22.0    // inter-compartmental CL 1
BTZ_Q3   = 5.6     // inter-compartmental CL 2
BTZ_F    = 0.83    // SC bioavailability

// ---- Drug PK: Carfilzomib (CFZ, IV) ----
// Based on Wang et al. Clin Pharmacokinet 2013
CFZ_CL   = 245.0   // clearance (L/hr) - high first-pass
CFZ_V1   = 8.0     // central volume (L)
CFZ_V2   = 12.0    // peripheral (L)
CFZ_Q    = 25.0

// ---- Drug PK: Lenalidomide (LEN, oral) ----
// Based on Chen et al. J Clin Pharmacol 2012
LEN_Ka   = 1.2     // absorption rate (1/hr)
LEN_CL   = 3.2     // clearance (L/hr)
LEN_V    = 65.0    // volume (L)
LEN_F    = 0.90    // oral bioavailability

// ---- Drug PK: Daratumumab (DARA, IV/SC) ----
// 2-cmt + target-mediated disposition (TMDD)
DARA_CL  = 0.0029  // linear clearance (L/hr/kg) - 70kg
DARA_V1  = 0.043   // central volume (L/kg)
DARA_V2  = 0.028   // peripheral (L/kg)
DARA_Q   = 0.0019  // inter-cmt CL
DARA_kon = 1.85    // CD38 binding on-rate (1/nM/day)
DARA_koff= 0.37    // off-rate (1/day)
DARA_kdeg= 0.83    // DARA-CD38 complex deg (1/day)
CD38_0   = 50.0    // baseline CD38 expression (nM)

// ---- Drug PK: Dexamethasone (DEX, oral) ----
DEX_Ka   = 2.5     // (1/hr)
DEX_CL   = 20.0    // (L/hr)
DEX_V    = 130.0   // (L)
DEX_F    = 0.78

// ---- Drug PK: Venetoclax (VEN, oral) ----
VEN_Ka   = 0.8     // (1/hr)
VEN_CL   = 14.0    // (L/hr)
VEN_V    = 256.0   // (L)
VEN_F    = 0.55    // food-dependent (with fatty meal)

// ---- Drug PK: Zoledronic Acid (ZOL, IV) ----
ZOL_CL   = 3.8     // (L/hr)
ZOL_V1   = 4.0     // (L)
ZOL_V2   = 22.0    // (L)
ZOL_Q    = 0.7
ZOL_kbone= 0.04    // bone uptake rate (1/hr)
ZOL_krel = 0.0001  // bone release rate (1/hr)

// ---- PD Parameters: Bortezomib (Proteasome inhibition) ----
// Proteasome inhibition → ER stress → apoptosis
BTZ_Emax = 0.90    // max additional apoptosis effect
BTZ_EC50 = 2.5     // EC50 (nM)
BTZ_hill = 2.0     // Hill coefficient

// ---- PD Parameters: Carfilzomib ----
CFZ_Emax = 0.95
CFZ_EC50 = 1.8     // nM (more potent than BTZ)
CFZ_hill = 2.5

// ---- PD Parameters: Lenalidomide (CRBN/IKZF1/3 → IRF4 → apoptosis) ----
LEN_Emax = 0.65    // max effect on proliferation suppression
LEN_EC50 = 0.8     // µg/mL
LEN_hill = 1.5

// ---- PD Parameters: Daratumumab (ADCC/CDC) ----
DARA_Emax = 0.85
DARA_EC50 = 0.15   // µg/mL (free DARA)
DARA_hill = 2.0

// ---- PD Parameters: Dexamethasone (GC-mediated apoptosis) ----
DEX_Emax = 0.55
DEX_EC50 = 0.10    // µg/mL
DEX_hill = 1.5

// ---- PD Parameters: Venetoclax (BCL-2 BH3 mimetic) ----
VEN_Emax = 0.80    // effective in t(11;14) / BCL-2 high
VEN_EC50 = 0.50    // µg/mL
VEN_hill = 1.8

// ---- PD: Bisphosphonate (OC apoptosis) ----
ZOL_OC_Emax = 0.75
ZOL_OC_EC50 = 0.05  // µg/mL bone-bound

$INIT
// Tumor compartments (cells)
MM_S    = 1e7    // sensitive MM cells at diagnosis
MM_R    = 1e4    // pre-existing resistant cells (rare)

// M-protein (g/dL) and sFLC (mg/dL)
MP      = 3.0    // baseline M-protein
FLC     = 150.0  // baseline sFLC

// IL-6, VEGF (AU)
IL6     = 0.1
VEGF    = 0.05

// Bone remodeling
OB      = 100.0  // osteoblasts (AU, normalized)
OC      = 100.0  // osteoclasts (AU)
BV      = 100.0  // bone volume (%)
NTX     = 1.0    // N-telopeptide (AU, bone resorption marker)
PINP    = 1.0    // P1NP (AU, bone formation marker)

// Biomarkers
Hgb     = 14.0   // hemoglobin (g/dL)
B2M     = 2.5    // β₂-microglobulin (mg/L)

// Drug PK compartments (amounts in µg for mAbs, ng for small molecules)
// Bortezomib (nM·L equivalent stored as ng)
BTZ1    = 0.0   // central
BTZ2    = 0.0   // peripheral 1
BTZ3    = 0.0   // peripheral 2

// Carfilzomib
CFZ1    = 0.0
CFZ2    = 0.0

// Lenalidomide
LEN_gut = 0.0
LEN1    = 0.0

// Daratumumab (µg)
DARA1   = 0.0   // central (µg)
DARA2   = 0.0   // peripheral
DARA_CD38 = 0.0 // bound complex

// Dexamethasone
DEX_gut = 0.0
DEX1    = 0.0

// Venetoclax
VEN_gut = 0.0
VEN1    = 0.0

// Zoledronic acid
ZOL1    = 0.0
ZOL2    = 0.0
ZOL_bone= 0.0

$ODE
// ==== DERIVED CONCENTRATIONS ====
double BTZ_Cp   = BTZ1 / BTZ_V1;          // nM equivalent
double CFZ_Cp   = CFZ1 / CFZ_V1;          // nM
double LEN_Cp   = LEN1 / LEN_V;           // µg/mL
double DARA_Cp  = DARA1 / (DARA_V1*70.0); // µg/mL (70 kg BW)
double DEX_Cp   = DEX1 / DEX_V;           // µg/mL
double VEN_Cp   = VEN1 / VEN_V;           // µg/mL
double ZOL_bone_Cp = ZOL_bone / ZOL_V2;  // µg/mL bone

// ==== DRUG EFFECTS (Emax models) ====
double E_BTZ = BTZ_Emax * pow(BTZ_Cp, BTZ_hill) /
               (pow(BTZ_EC50, BTZ_hill) + pow(BTZ_Cp, BTZ_hill));
double E_CFZ = CFZ_Emax * pow(CFZ_Cp, CFZ_hill) /
               (pow(CFZ_EC50, CFZ_hill) + pow(CFZ_Cp, CFZ_hill));
double E_LEN = LEN_Emax * pow(LEN_Cp, LEN_hill) /
               (pow(LEN_EC50, LEN_hill) + pow(LEN_Cp, LEN_hill));
double E_DARA= DARA_Emax * pow(DARA_Cp, DARA_hill) /
               (pow(DARA_EC50, DARA_hill) + pow(DARA_Cp, DARA_hill));
double E_DEX = DEX_Emax * pow(DEX_Cp, DEX_hill) /
               (pow(DEX_EC50, DEX_hill) + pow(DEX_Cp, DEX_hill));
double E_VEN = VEN_Emax * pow(VEN_Cp, VEN_hill) /
               (pow(VEN_EC50, VEN_hill) + pow(VEN_Cp, VEN_hill));
double E_ZOL = ZOL_OC_Emax * ZOL_bone_Cp /
               (ZOL_OC_EC50 + ZOL_bone_Cp);

double total_kill = E_BTZ + E_CFZ + E_LEN + E_DARA + E_DEX + E_VEN;
if(total_kill > 0.99) total_kill = 0.99;

// ==== TUMOR CELLS ====
double N_total = MM_S + MM_R;
double IL6_effect = 1.0 + kIL6_stim * IL6;

dxdt_MM_S = kgrow * MM_S * (1.0 - N_total/Kmax) * IL6_effect
            - kapop * MM_S
            - (kapop * total_kill / (1.0 - total_kill)) * MM_S
            - kres * MM_S;

dxdt_MM_R = kgrow_res * MM_R * (1.0 - N_total/Kmax)
            - kapop * MM_R
            + kres * MM_S;

// ==== M-PROTEIN ====
dxdt_MP  = kMP_prod * MM_S - kMP_elim * MP;

// ==== sFLC ====
dxdt_FLC = kFLC_prod * MM_S - kFLC_elim * FLC;

// ==== CYTOKINES ====
double MM_burden = N_total / 1e9; // normalized
dxdt_IL6  = kIL6_prod + 0.5*MM_burden - kIL6_elim * IL6;
dxdt_VEGF = kVEGF_prod + 0.3*MM_burden - kVEGF_elim * VEGF;

// ==== BONE REMODELING (Lemaire-type) ====
double RANKL_eff = 1.0 + k_RANKL * MM_burden * 1e3;
double OPG_eff   = k_OPG * OB / (1.0 + k_OPG_inh * MM_burden);
double DKK1      = DKK1_base + kDKK1_MM * MM_burden;
double Wnt_inh   = 1.0 / (1.0 + DKK1);  // reduced OB activity

dxdt_OB = k_OB_form * Wnt_inh - k_OB_apop * OB
          - 0.0001 * (1.0 + E_DEX*0.3) * OB; // GC suppresses OB

dxdt_OC = k_OC_form * RANKL_eff / (1.0 + OPG_eff)
          - k_OC_apop * OC * (1.0 + E_ZOL); // ZOL → OC apoptosis

dxdt_BV = k_bone_form * OB - k_bone_res * OC;
if(BV < 0) BV = 0;

dxdt_NTX  = k_bone_res * OC - 0.5 * NTX;  // bone resorption marker
dxdt_PINP = k_bone_form * OB - 0.5 * PINP; // bone formation marker

// ==== HEMOGLOBIN ====
double IL6_Hgb_eff = 1.0 + IL6_Hgb_inh * IL6;
dxdt_Hgb = kHgb_recover * (Hgb0 - Hgb) / IL6_Hgb_eff
           - kHgb_supp * MM_burden;

// ==== BETA2-MICROGLOBULIN ====
dxdt_B2M = kB2M_prod * N_total - kB2M_elim * B2M;

// ==== DRUG PK: BORTEZOMIB (3-cmt, IV/SC) ====
double BTZ_dose_rate = 0.0; // driven by event
dxdt_BTZ1 = BTZ_dose_rate
            - (BTZ_CL + BTZ_Q2 + BTZ_Q3)/BTZ_V1 * BTZ1
            + BTZ_Q2/BTZ_V2 * BTZ2
            + BTZ_Q3/BTZ_V3 * BTZ3;
dxdt_BTZ2 = BTZ_Q2/BTZ_V1 * BTZ1 - BTZ_Q2/BTZ_V2 * BTZ2;
dxdt_BTZ3 = BTZ_Q3/BTZ_V1 * BTZ1 - BTZ_Q3/BTZ_V3 * BTZ3;

// ==== DRUG PK: CARFILZOMIB (2-cmt, IV) ====
dxdt_CFZ1 = -(CFZ_CL + CFZ_Q)/CFZ_V1 * CFZ1 + CFZ_Q/CFZ_V2 * CFZ2;
dxdt_CFZ2 =  CFZ_Q/CFZ_V1 * CFZ1 - CFZ_Q/CFZ_V2 * CFZ2;

// ==== DRUG PK: LENALIDOMIDE (1-cmt oral) ====
dxdt_LEN_gut = -LEN_Ka * LEN_gut;
dxdt_LEN1    =  LEN_Ka * LEN_gut * LEN_F - (LEN_CL/LEN_V) * LEN1;

// ==== DRUG PK: DARATUMUMAB (2-cmt + TMDD) ====
double CD38_free = CD38_0 - DARA_CD38;
if(CD38_free < 0) CD38_free = 0;
dxdt_DARA1    = -(DARA_CL + DARA_Q)*70.0/DARA_V1/70.0 * DARA1
               + DARA_Q*70.0/DARA_V2/70.0 * DARA2
               - DARA_kon * DARA1 * CD38_free
               + DARA_koff * DARA_CD38;
dxdt_DARA2    =  DARA_Q*70.0/DARA_V1/70.0 * DARA1
               - DARA_Q*70.0/DARA_V2/70.0 * DARA2;
dxdt_DARA_CD38 = DARA_kon * DARA1 * CD38_free
                - (DARA_koff + DARA_kdeg) * DARA_CD38;

// ==== DRUG PK: DEXAMETHASONE (1-cmt oral) ====
dxdt_DEX_gut = -DEX_Ka * DEX_gut;
dxdt_DEX1    =  DEX_Ka * DEX_gut * DEX_F - (DEX_CL/DEX_V) * DEX1;

// ==== DRUG PK: VENETOCLAX (1-cmt oral) ====
dxdt_VEN_gut = -VEN_Ka * VEN_gut;
dxdt_VEN1    =  VEN_Ka * VEN_gut * VEN_F - (VEN_CL/VEN_V) * VEN1;

// ==== DRUG PK: ZOLEDRONIC ACID (2-cmt + bone) ====
dxdt_ZOL1    = -(ZOL_CL + ZOL_Q)/ZOL_V1 * ZOL1 + ZOL_Q/ZOL_V2 * ZOL2
               - ZOL_kbone * ZOL1;
dxdt_ZOL2    =  ZOL_Q/ZOL_V1 * ZOL1 - ZOL_Q/ZOL_V2 * ZOL2;
dxdt_ZOL_bone = ZOL_kbone * ZOL1 - ZOL_krel * ZOL_bone;

$TABLE
// Response criteria
double M_protein_pct = (MP / 3.0) * 100.0;  // % of baseline
double MM_total = MM_S + MM_R;
double BTZ_Cp_out  = BTZ1 / BTZ_V1;
double CFZ_Cp_out  = CFZ1 / CFZ_V1;
double LEN_Cp_out  = LEN1 / LEN_V;
double DARA_Cp_out = DARA1 / (DARA_V1*70.0);
double DEX_Cp_out  = DEX1 / DEX_V;
double VEN_Cp_out  = VEN1 / VEN_V;
double E_drug_total = E_BTZ + E_CFZ + E_LEN + E_DARA + E_DEX + E_VEN;
double ISS_score = (B2M > 3.5 && Albumin_sim < 3.5) ? 3 : (B2M > 3.5 ? 2 : 1);

// Pseudo-albumin (declines with disease burden)
double Albumin_sim = 4.0 - 0.3 * (MM_total / 1e9);

$CAPTURE
MM_S MM_R MM_total MP FLC IL6 VEGF
OB OC BV NTX PINP Hgb B2M
BTZ_Cp_out CFZ_Cp_out LEN_Cp_out DARA_Cp_out DEX_Cp_out VEN_Cp_out
E_BTZ E_CFZ E_LEN E_DARA E_DEX E_VEN E_drug_total
M_protein_pct Albumin_sim ISS_score
'

mm_mod <- mread_cache("mm_qsp", tempdir(), mm_model_code)

## ============================================================
## DOSING REGIMENS
## ============================================================

# Helper to build standard regimen event tables
make_BTZ_events <- function(start=0, n_cycles=8) {
  # Bortezomib VD: 1.3 mg/m² SC days 1,4,8,11 of 21-day cycle
  # Assume BSA = 1.8 m² → dose = 2.34 mg → ~2340 ng IV bolus into BTZ1
  dose_ng <- 2.34e3  # ng (equivalent for nM in V1=4.7L)
  days <- unlist(lapply(0:(n_cycles-1), function(c) {
    c*21 + c(1, 4, 8, 11)
  }))
  ev(amt=dose_ng, cmt="BTZ1", time=days*24, ii=0, addl=0)
}

make_LEN_events <- function(start=0, n_cycles=8, dose_mg=25) {
  # Lenalidomide 25 mg/day PO days 1-21 of 28-day cycle
  dose_ug <- dose_mg * 1e3
  days <- unlist(lapply(0:(n_cycles-1), function(c) {
    (c*28 + 0:20)
  }))
  ev(amt=dose_ug * LEN_F, cmt="LEN_gut", time=days*24, ii=0, addl=0)
}

# Use explicit parameters for simulation
LEN_F <- 0.90

make_DARA_events <- function(n_weeks_q1w=8, n_q2w=8, n_q4w=12) {
  # Daratumumab 16 mg/kg IV: weekly x8, q2w x8, q4w thereafter
  dose_ug <- 16e3 * 70  # µg (16 mg/kg × 70 kg)
  t_q1w <- seq(0, (n_weeks_q1w-1)*7, by=7) * 24
  t_q2w <- (n_weeks_q1w*7 + seq(0, (n_q2w-1)*14, by=14)) * 24
  t_q4w <- (n_weeks_q1w*7 + n_q2w*14 +
              seq(0, (n_q4w-1)*28, by=28)) * 24
  ev(amt=dose_ug, cmt="DARA1", time=c(t_q1w, t_q2w, t_q4w))
}

make_DEX_events <- function(n_cycles=8, dose_mg=40) {
  # Dexamethasone 40 mg PO days 1,8,15,22 of 28-day cycle
  dose_ug <- dose_mg * 1e3
  days <- unlist(lapply(0:(n_cycles-1), function(c) {
    c*28 + c(1, 8, 15, 22)
  }))
  ev(amt=dose_ug * DEX_F, cmt="DEX_gut", time=days*24, ii=0, addl=0)
}

DEX_F <- 0.78

make_ZOL_events <- function(n_infusions=12) {
  # Zoledronic acid 4 mg IV q4 weeks x 12 cycles = 1 year
  dose_ug <- 4e3
  times <- seq(0, (n_infusions-1)*28, by=28) * 24
  ev(amt=dose_ug, cmt="ZOL1", time=times)
}

## ============================================================
## SCENARIO 1: Untreated disease progression
## ============================================================
sim_time <- seq(0, 720*24, by=24)  # 720 days

s1_obs <- mm_mod %>%
  mrgsim(end=720*24, delta=24) %>%
  as_tibble() %>%
  mutate(scenario = "No Treatment", day = time/24)

## ============================================================
## SCENARIO 2: VRd (Bortezomib + Lenalidomide + Dexamethasone)
## Standard frontline regimen (SWOG S0777)
## ============================================================
ev_VRd <- as_data_frame(
  bind_rows(
    make_BTZ_events(n_cycles=8),
    ev(amt=25e3*LEN_F, cmt="LEN_gut",
       time=seq(0, 7*28, by=28)*24 + seq(0,20)*24,
       ii=0, addl=0) |> slice(rep(1:21, 8)) |>
      mutate(time = unlist(lapply(0:7, function(c) (c*28 + 0:20)*24))),
    make_DEX_events(n_cycles=8, dose_mg=20)
  )
) %>% distinct()

s2_obs <- mm_mod %>%
  data_set(ev(amt=25e3*LEN_F, cmt="LEN_gut",
              time=seq(0, 8*28-1)*24, ii=24*7, addl=3) %>%
             as_data_frame() %>%
             bind_rows(make_BTZ_events(n_cycles=8) %>% as_data_frame()) %>%
             bind_rows(make_DEX_events(n_cycles=8, dose_mg=20) %>%
                         as_data_frame())) %>%
  mrgsim(end=720*24, delta=24) %>%
  as_tibble() %>%
  mutate(scenario = "VRd (Bortezomib+Len+Dex)", day = time/24)

## ============================================================
## SCENARIO 3: Daratumumab + Rd (DRd) — MAIA trial regimen
## For transplant-ineligible NDMM
## ============================================================
s3_obs <- mm_mod %>%
  data_set(
    bind_rows(
      make_DARA_events(n_weeks_q1w=8, n_q2w=8, n_q4w=8) %>%
        as_data_frame(),
      ev(amt=25e3*LEN_F, cmt="LEN_gut",
         time=seq(0, 8*28-1)*24, ii=24*7, addl=3) %>%
        as_data_frame(),
      make_DEX_events(n_cycles=8, dose_mg=40) %>% as_data_frame()
    )
  ) %>%
  mrgsim(end=720*24, delta=24) %>%
  as_tibble() %>%
  mutate(scenario = "DRd (Dara+Len+Dex)", day = time/24)

## ============================================================
## SCENARIO 4: Carfilzomib + Lenalidomide + Dexamethasone (KRd)
## ASPIRE trial — relapsed/refractory MM
## ============================================================
s4_obs <- mm_mod %>%
  data_set(
    bind_rows(
      ev(amt=36e3, cmt="CFZ1",
         time=c(1,2,8,9,15,16)*24, ii=21*24, addl=7) %>%
        as_data_frame(),
      ev(amt=25e3*LEN_F, cmt="LEN_gut",
         time=seq(0, 8*28-1)*24, ii=24*7, addl=3) %>%
        as_data_frame(),
      make_DEX_events(n_cycles=8) %>% as_data_frame()
    )
  ) %>%
  mrgsim(end=720*24, delta=24) %>%
  as_tibble() %>%
  mutate(scenario = "KRd (Carfilzomib+Len+Dex)", day = time/24)

## ============================================================
## SCENARIO 5: Venetoclax + Dexamethasone (t(11;14) subtype)
## BCL-2 high expression — Venetoclax + Dex
## ============================================================
mm_mod_bcl2 <- mm_mod %>%
  param(VEN_Emax = 0.92,   # t(11;14) high BCL-2 — very sensitive
        kgrow    = 0.018)

s5_obs <- mm_mod_bcl2 %>%
  data_set(
    bind_rows(
      ev(amt=800e3*VEN_F, cmt="VEN_gut",
         time=seq(0, 720)*24, ii=0, addl=0) %>%
        as_data_frame(),
      make_DEX_events(n_cycles=24, dose_mg=40) %>% as_data_frame()
    )
  ) %>%
  mrgsim(end=720*24, delta=24) %>%
  as_tibble() %>%
  mutate(scenario = "VenDex (t(11;14)/BCL-2 high)", day = time/24)

## ============================================================
## SCENARIO 6: High-risk disease — Daratumumab + VRd (DVRd)
## PERSEUS trial for NDMM
## ============================================================
s6_obs <- mm_mod %>%
  param(kgrow = 0.022, del17p_adj = 1.3,  # faster growth high-risk
        BTZ_Emax = 0.75) %>%  # slightly reduced BTZ effect del17p
  data_set(
    bind_rows(
      make_DARA_events(n_weeks_q1w=8, n_q2w=8, n_q4w=12) %>%
        as_data_frame(),
      make_BTZ_events(n_cycles=8) %>% as_data_frame(),
      ev(amt=25e3*LEN_F, cmt="LEN_gut",
         time=seq(0, 8*28-1)*24, ii=24*7, addl=3) %>%
        as_data_frame(),
      make_DEX_events(n_cycles=8) %>% as_data_frame()
    )
  ) %>%
  mrgsim(end=720*24, delta=24) %>%
  as_tibble() %>%
  mutate(scenario = "DVRd High-risk (PERSEUS)", day = time/24)

## ============================================================
## COMBINE ALL SCENARIOS
## ============================================================
all_scenarios <- bind_rows(s1_obs, s2_obs, s3_obs, s4_obs, s5_obs, s6_obs)

## ============================================================
## VISUALIZATION
## ============================================================

p1 <- ggplot(all_scenarios, aes(x=day, y=MP, color=scenario)) +
  geom_line(size=1) +
  labs(title="M-Protein Response by Treatment Regimen",
       x="Day", y="M-Protein (g/dL)",
       color="Regimen") +
  geom_hline(yintercept=c(1.5, 0.6, 0.3), linetype="dashed",
             color=c("orange","blue","green"),
             alpha=0.7) +
  annotate("text", x=700, y=c(1.55, 0.65, 0.35),
           label=c("PR", "VGPR", "CR"), size=3) +
  scale_color_brewer(palette="Dark2") +
  theme_bw() +
  theme(legend.position="bottom",
        legend.text=element_text(size=7))

p2 <- ggplot(all_scenarios, aes(x=day, y=log10(MM_total+1), color=scenario)) +
  geom_line(size=1) +
  labs(title="Total MM Cell Burden (log10)",
       x="Day", y="log10(MM cells + 1)",
       color="Regimen") +
  scale_color_brewer(palette="Dark2") +
  theme_bw() +
  theme(legend.position="bottom",
        legend.text=element_text(size=7))

p3 <- ggplot(all_scenarios %>% filter(scenario %in%
               c("No Treatment","VRd (Bortezomib+Len+Dex)")),
             aes(x=day, y=BV, color=scenario)) +
  geom_line(size=1) +
  labs(title="Bone Volume (% baseline)", x="Day", y="Bone Volume (%)") +
  scale_color_brewer(palette="Set1") +
  theme_bw()

p4 <- ggplot(all_scenarios, aes(x=day, y=Hgb, color=scenario)) +
  geom_line(size=1) +
  geom_hline(yintercept=10, linetype="dashed", color="red") +
  annotate("text", x=50, y=10.2, label="Anemia threshold (10 g/dL)",
           size=3, color="red") +
  labs(title="Hemoglobin over Time",
       x="Day", y="Hemoglobin (g/dL)") +
  scale_color_brewer(palette="Dark2") +
  theme_bw() +
  theme(legend.position="bottom",
        legend.text=element_text(size=7))

## PK profiles — day 1 BTZ cycle 1
btz_pk <- mm_mod %>%
  ev(amt=2340, cmt="BTZ1", time=0) %>%
  mrgsim(end=72, delta=0.5) %>%
  as_tibble() %>%
  mutate(BTZ_Cp = BTZ1 / 4.7)

p5 <- ggplot(btz_pk, aes(x=time, y=BTZ_Cp)) +
  geom_line(color="#7e22ce", size=1.2) +
  geom_hline(yintercept=2.5, linetype="dashed", color="red") +
  annotate("text", x=50, y=2.7, label="EC50 BTZ", size=3, color="red") +
  labs(title="Bortezomib PK (Single IV Dose 1.3 mg/m²)",
       x="Time (hr)", y="Plasma Concentration (nM)") +
  theme_bw()

## Response depth comparison at day 360
resp_d360 <- all_scenarios %>%
  filter(day == 360) %>%
  mutate(
    resp_cat = case_when(
      MP < 0.1 * 3.0 ~ "CR",      # <10% of baseline
      MP < 0.3 * 3.0 ~ "VGPR",    # <30% of baseline
      MP < 1.5 * 3.0 ~ "PR",      # <50% reduction from 3.0
      MP >= 1.5 * 3.0 & MP < 3.75 ~ "SD",
      TRUE ~ "PD"
    )
  )

print(resp_d360 %>% select(scenario, MP, FLC, Hgb, BV, B2M, resp_cat))

## ============================================================
## POPULATION SIMULATION (variability in response)
## ============================================================
pop_params <- expand.grid(
  kgrow   = rnorm(50, 0.015, 0.003),
  BTZ_EC50 = rlnorm(50, log(2.5), 0.4),
  LEN_EC50 = rlnorm(50, log(0.8), 0.4)
) %>%
  slice_sample(n=50) %>%
  mutate(ID = row_number())

pop_sim <- lapply(1:nrow(pop_params), function(i) {
  mm_mod %>%
    param(kgrow   = pop_params$kgrow[i],
          BTZ_EC50 = pop_params$BTZ_EC50[i],
          LEN_EC50 = pop_params$LEN_EC50[i]) %>%
    data_set(
      bind_rows(
        make_BTZ_events(n_cycles=8) %>% as_data_frame(),
        ev(amt=25e3*LEN_F, cmt="LEN_gut",
           time=seq(0, 8*28-1)*24, ii=24*7, addl=3) %>%
          as_data_frame(),
        make_DEX_events(n_cycles=8, dose_mg=20) %>%
          as_data_frame()
      )
    ) %>%
    mrgsim(end=720*24, delta=24) %>%
    as_tibble() %>%
    mutate(ID = pop_params$ID[i], day=time/24)
}) %>% bind_rows()

p6 <- ggplot(pop_sim, aes(x=day, y=MP, group=ID)) +
  geom_line(alpha=0.2, color="#0077b6") +
  stat_summary(aes(group=1), fun=median,
               geom="line", color="red", size=1.2) +
  stat_summary(aes(group=1), fun=function(x) quantile(x, 0.25),
               geom="line", color="red", size=0.8, linetype=2) +
  stat_summary(aes(group=1), fun=function(x) quantile(x, 0.75),
               geom="line", color="red", size=0.8, linetype=2) +
  geom_hline(yintercept=c(1.5, 0.6), linetype="dashed",
             color=c("orange","blue")) +
  labs(title="Population PK/PD Simulation — VRd (n=50)\nM-Protein Response with Variability",
       subtitle="Red = median; dashed = 25th/75th percentile",
       x="Day", y="M-Protein (g/dL)") +
  theme_bw()

## Print all plots
print(p1); print(p2); print(p3); print(p4); print(p5); print(p6)

cat("\n=== Simulation Summary ===\n")
cat("Day 360 Response Summary:\n")
print(resp_d360 %>%
        select(scenario, day, MP, FLC, Hgb, BV, resp_cat) %>%
        arrange(MP))
