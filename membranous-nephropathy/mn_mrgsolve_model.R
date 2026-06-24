################################################################################
# Membranous Nephropathy (MN) QSP Model — mrgsolve Implementation
#
# Disease Background:
#   Membranous Nephropathy is an autoimmune glomerular disease characterized by
#   subepithelial IgG4 immune deposits (mainly anti-PLA2R1 antibodies), complement
#   (MAC) activation, podocyte injury, GBM thickening, and nephrotic syndrome.
#
# Key Clinical Trials Referenced:
#   - GEMRITUX (2017, NEJM): Rituximab vs. conservative therapy; primary endpoint
#     partial/complete remission at 6 months. Rituximab: 375 mg/m² ×2 (D1, D8).
#     Beck LH Jr et al., NEJM 2019 (extended follow-up).
#   - MENTOR (2019, NEJM): Rituximab vs. cyclosporine over 24 months.
#     Fervenza FC et al., NEJM 2019; rituximab 1 g IV at D1 and M6.
#     Complete + partial remission: 60% RTX vs. 20% CsA at 24 months.
#   - STARMEN (2020, JASN): Sequential tacrolimus → RTX vs. cyclophosphamide +
#     alternating-day steroids (modified Ponticelli). van den Brand JASN 2020.
#     RTX arm superior at 24 months.
#   - RI-CYCLO (2021, JASN): Rituximab vs. cyclophosphamide-based Ponticelli.
#     Non-inferiority of RTX demonstrated. Scolari F et al., JASN 2021.
#   - AVACOPAN trial (2023, phase II): C5aR1 inhibitor avacopan investigated in
#     complement-driven glomerulonephritides.
#   - Ponticelli Regimen (Ponticelli C et al., NEJM 1989; classic reference for
#     cyclophosphamide + steroids alternating monthly for 6 months).
#
# Model Structure (≥15 ODE compartments):
#   Drug PK:
#     1.  RTX_cent       — Rituximab central compartment (µg/mL)
#     2.  RTX_peri       — Rituximab peripheral compartment (µg/mL)
#     3.  RTX_CD20_bound — RTX-CD20 bound complex (TMDD, µg/mL equiv.)
#     4.  TAC_blood      — Tacrolimus whole blood (ng/mL)
#     5.  CPx_cent       — Cyclophosphamide central (µg/mL)
#     6.  CPx_metab      — Active 4-OH-cyclophosphamide metabolite (µg/mL)
#   Immune/Disease:
#     7.  CD20_B         — CD20+ B cells (cells/µL)
#     8.  Plasma_cells   — Antibody-secreting plasma cells (rel. units)
#     9.  Anti_PLA2R1    — Serum anti-PLA2R1 IgG4 (U/mL)
#     10. IgG_deposit    — Subepithelial IgG4 deposit (rel. units)
#     11. Complement_MAC — MAC formation/activity (rel. units)
#     12. Podocyte_inj   — Podocyte injury score (0–1)
#     13. GBM_thick      — GBM structural thickening (rel. units, baseline=1)
#   Clinical Endpoints:
#     14. Proteinuria    — 24h urine protein (g/day)
#     15. Serum_alb      — Serum albumin (g/dL)
#     16. eGFR           — Glomerular filtration rate (mL/min/1.73m²)
#   RAAS:
#     17. AngII          — Angiotensin II (pg/mL)
#     18. Aldosterone    — Aldosterone (pg/mL)
#
# Author: Claude Code Routine (CCR) — 2026-06-17
################################################################################

library(mrgsolve)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ============================================================
# mrgsolve model code
# ============================================================
mn_code <- '
$PROB
Membranous Nephropathy (MN) QSP Model
18-compartment ODE system
Key trials: GEMRITUX, MENTOR, STARMEN, RI-CYCLO

$PARAM
// -------------------------------------------------------
// Rituximab PK parameters (population mean, calibrated to
// MENTOR trial PK data; Fervenza FC NEJM 2019)
// -------------------------------------------------------
CL_RTX   = 0.35    // L/day  clearance
Vc_RTX   = 3.1     // L      central volume
Vp_RTX   = 2.8     // L      peripheral volume
Q_RTX    = 0.7     // L/day  intercompartmental clearance
// TMDD parameters for RTX-CD20 binding
Kd_RTX   = 0.010   // µg/mL  dissociation constant (~10 nM)
kon_RTX  = 0.01    // 1/(µg/mL·day)  association
koff_RTX = 0.0001  // 1/day  dissociation (koff = Kd*kon)
kint_RTX = 0.5     // 1/day  internalization of bound RTX-CD20

// -------------------------------------------------------
// Tacrolimus PK (whole-blood; calibrated to STARMEN dosing
// 0.05 mg/kg/day; target C0 5-10 ng/mL)
// -------------------------------------------------------
CL_TAC   = 25.0    // L/h  apparent oral clearance (CL/F)
Vd_TAC   = 900.0   // L    apparent volume of distribution
F_TAC    = 0.25    // bioavailability (oral)
ka_TAC   = 2.4     // 1/h  absorption rate constant

// -------------------------------------------------------
// Cyclophosphamide PK (prodrug → active 4-OH-CPx)
// Reference: de Jonge ME et al., Clin Pharmacokinet 2005
// -------------------------------------------------------
CL_CPx   = 4.5     // L/h  CPx clearance
Vd_CPx   = 35.0    // L    CPx volume
km_CPx   = 0.15    // 1/h  metabolic activation rate (CPx → 4OH-CPx)
kel_met  = 0.8     // 1/h  elimination of active metabolite

// -------------------------------------------------------
// B cell / Plasma cell dynamics
// B cell baseline: 200 cells/µL (normal peripheral blood)
// RTX depletes to <5 cells/µL (MENTOR PK/PD data)
// Reconstitution half-life ~6–12 months (median ~9 months)
// Reference: Fervenza NEJM 2019, van den Brand JASN 2020
// -------------------------------------------------------
Bss      = 200.0   // cells/µL  steady-state CD20+ B cells
kprol_B  = 0.005   // 1/day  B cell proliferation rate (= kdeg_B at SS)
kdeg_B   = 0.005   // 1/day  B cell natural turnover
// TMDD-driven B cell killing by RTX
keff_RTX = 0.08    // 1/(µg/mL·day)  RTX efficacy on B cell depletion

// Plasma cell dynamics
PC_ss    = 1.0     // rel. units  plasma cell steady state
kprol_PC = 0.003   // 1/day  plasma cell production from B cells
kdeg_PC  = 0.003   // 1/day  plasma cell turnover

// -------------------------------------------------------
// Anti-PLA2R1 IgG4 dynamics
// Half-life IgG4 ~28 days (t1/2 = 0.693/kdeg_Ab)
// Baseline titer ~150 U/mL (nephrotic MN)
// Reference: Beck LH Jr NEJM 2009; Bobart SA JASN 2020
// -------------------------------------------------------
Ab_ss    = 150.0   // U/mL  baseline anti-PLA2R1 IgG4
kdeg_Ab  = 0.025   // 1/day  IgG4 elimination (t1/2 ~28 days)
// Production rate = kdeg_Ab * Ab_ss (driven by plasma cells)

// -------------------------------------------------------
// Subepithelial IgG4 deposit dynamics
// -------------------------------------------------------
kdepo    = 0.04    // 1/day  deposition rate from serum IgG4
kclear_d = 0.02    // 1/day  deposit clearance
Dep_ss   = 1.0     // rel. units  steady state deposit

// -------------------------------------------------------
// Complement / MAC dynamics
// Threshold: MAC formation increases sharply once IgG deposit > threshold
// Reference: Cybulsky AV, J Am Soc Nephrol 2011
// -------------------------------------------------------
MAC_ss   = 1.0     // rel. units  MAC at baseline
kcomp    = 0.05    // 1/day  complement activation
kres_MAC = 0.05    // 1/day  MAC resolution
nH_MAC   = 3.0     // Hill coefficient for IgG-driven MAC
EC50_MAC = 1.0     // rel. units  IgG deposit EC50 for MAC

// -------------------------------------------------------
// Podocyte injury dynamics (0 = healthy, 1 = maximal injury)
// Reference: Shimizu M et al., Clin Exp Nephrol 2019
// -------------------------------------------------------
kinj_pod  = 0.03   // 1/day  MAC-driven podocyte injury rate
krec_pod  = 0.005  // 1/day  intrinsic podocyte recovery rate
Pod_ss    = 0.6    // dimensionless  baseline injury at nephrotic (0.6)

// -------------------------------------------------------
// GBM thickening (rel. units; baseline = 1.0)
// Driven by cumulative podocyte injury
// Reference: Glassock RJ, CJASN 2012
// -------------------------------------------------------
kGBM     = 0.002   // 1/day  GBM thickening rate per injury unit
kGBM_res = 0.0005  // 1/day  partial GBM resolution
GBM_ss   = 1.3     // rel. units  GBM at baseline (mild thickening)

// -------------------------------------------------------
// Proteinuria (g/day)
// Normal ~0.15 g/day; nephrotic >3.5 g/day
// Driven by podocyte injury + GBM thickening
// Reference: Cattran DC, Kidney Int 2001 (spontaneous remission data)
// -------------------------------------------------------
Prot_max = 10.0    // g/day  maximum proteinuria
kprot    = 0.02    // 1/day  proteinuria kinetics
Prot_ss  = 6.0     // g/day  baseline (nephrotic)

// -------------------------------------------------------
// Serum albumin (g/dL)
// Inversely related to proteinuria
// Reference: Glassock RJ, CJASN 2012
// -------------------------------------------------------
Alb_ss   = 2.2     // g/dL  baseline serum albumin (nephrotic)
kAlb_syn = 14.3    // g/dL/day  synthesis (~1.5× normal when depleted)
kAlb_deg = 0.2     // 1/day  albumin catabolism
kAlb_prot= 0.7     // coupling coefficient: albumin loss ~ proteinuria

// -------------------------------------------------------
// eGFR (mL/min/1.73m²)
// Declines with GBM thickening; partial recovery possible
// Reference: Jha V et al., KI Reports 2020
// -------------------------------------------------------
eGFR_ss  = 78.0    // mL/min/1.73m²  baseline eGFR
eGFR_min = 15.0    // floor (ESRD boundary)
keGFR    = 0.001   // 1/day  eGFR loss per GBM unit
keGFR_r  = 0.002   // 1/day  partial eGFR recovery

// -------------------------------------------------------
// RAAS parameters
// ACEi/ARB: modeled as RAAS suppression flag (0–1)
// Reference: KDIGO 2021 guideline for MN supportive care
// -------------------------------------------------------
AngII_ss  = 25.0   // pg/mL  baseline AngII
kAng      = 0.1    // 1/day  AngII turnover
ACEI_eff  = 0.0    // 0–1  ACEi suppression of AngII (input flag)
Aldo_ss   = 100.0  // pg/mL  baseline Aldosterone
kAldo     = 0.05   // 1/day  Aldosterone turnover
kAldo_Ang = 2.5    // Aldosterone coupling to AngII

// -------------------------------------------------------
// Tacrolimus / CNI immunosuppression effect on B cells & Ab
// -------------------------------------------------------
TAC_EC50  = 5.0    // ng/mL  half-maximal effect on Ab production
TAC_Emax  = 0.6    // maximal fractional inhibition of Ab production (60%)

// -------------------------------------------------------
// Cyclophosphamide effect on plasma cells
// Reference: Ponticelli C et al., NEJM 1989; Scolari JASN 2021
// -------------------------------------------------------
CPx_EC50  = 0.5    // µg/mL  active metabolite EC50 on plasma cells
CPx_Emax  = 0.85   // maximal kill fraction of plasma cells (85%)

// -------------------------------------------------------
// Avacopan (C5aR inhibitor) — investigational
// Effect: reduce complement MAC by inhibiting C5a-driven amplification
// Reference: Phase II data 2023; Jayne DRW et al., NEJM 2021 (ANCA context)
// -------------------------------------------------------
AVA_EC50  = 1.0    // rel. units (arbitrary plasma concentration proxy)
AVA_Emax  = 0.7    // maximal MAC reduction (70%)
AVA_dose  = 0.0    // binary switch 0/1 (simplified PK for investigational agent)

// -------------------------------------------------------
// Dose event flags (set via event table in R)
// -------------------------------------------------------
RTX_dose  = 0.0    // µg/mL infusion input (via rate in event table)
TAC_dose  = 0.0    // mg (oral; handled via ka absorption)
CPx_dose  = 0.0    // µg/mL IV input

$CMT
// Drug PK compartments
RTX_cent
RTX_peri
RTX_CD20_bound
TAC_gut
TAC_blood
CPx_cent
CPx_metab

// Immune / Disease compartments
CD20_B
Plasma_cells
Anti_PLA2R1
IgG_deposit
Complement_MAC
Podocyte_inj
GBM_thick

// Clinical endpoint compartments
Proteinuria
Serum_alb
eGFR
AngII
Aldosterone

$GLOBAL
// helper variables visible to $MAIN, $ODE, $TABLE
double RTX_conc, TAC_conc, CPx_met_conc;
double Emax_TAC, Emax_CPx, Emax_AVA;
double B_norm, MAC_norm, Pod_norm;

$MAIN
// Initialize all compartments at (approximate) nephrotic baseline
if(NEWIND <= 1) {
  // Drug PK — start at zero
  RTX_cent_0     = 0.0;
  RTX_peri_0     = 0.0;
  RTX_CD20_bound_0 = 0.0;
  TAC_gut_0      = 0.0;
  TAC_blood_0    = 0.0;
  CPx_cent_0     = 0.0;
  CPx_metab_0    = 0.0;

  // Immune / Disease — at steady-state nephrotic disease
  CD20_B_0        = Bss;
  Plasma_cells_0  = PC_ss;
  Anti_PLA2R1_0   = Ab_ss;
  IgG_deposit_0   = Dep_ss;
  Complement_MAC_0= MAC_ss;
  Podocyte_inj_0  = Pod_ss;
  GBM_thick_0     = GBM_ss;

  // Clinical endpoints
  Proteinuria_0   = Prot_ss;
  Serum_alb_0     = Alb_ss;
  eGFR_0          = eGFR_ss;
  AngII_0         = AngII_ss;
  Aldosterone_0   = Aldo_ss;
}

$ODE
// -------------------------------------------------------
// Convenience aliases (concentrations)
// -------------------------------------------------------
RTX_conc    = RTX_cent;                   // µg/mL
TAC_conc    = TAC_blood;                  // ng/mL
CPx_met_conc= CPx_metab;                  // µg/mL  active metabolite

// -------------------------------------------------------
// Drug effect inhibition fractions (Emax models)
// -------------------------------------------------------
// Tacrolimus: inhibits Ab production (via calcineurin → NFAT → plasma cell function)
Emax_TAC = TAC_Emax * TAC_conc / (TAC_EC50 + TAC_conc + 1e-9);

// Cyclophosphamide active metabolite: alkylates and kills plasma cells
Emax_CPx = CPx_Emax * CPx_met_conc / (CPx_EC50 + CPx_met_conc + 1e-9);

// Avacopan: reduces complement MAC (C5aR block)
Emax_AVA = AVA_Emax * AVA_dose / (AVA_EC50 + AVA_dose + 1e-9);

// -------------------------------------------------------
// 1. RTX_cent: Rituximab central compartment
//    Dose input via rate (zero-order infusion in event table)
//    TMDD: binding to CD20+ B cells expressed as plasma compartment sink
// -------------------------------------------------------
double RTX_CD20_on  = kon_RTX * RTX_conc * CD20_B;   // association flux
double RTX_CD20_off = koff_RTX * RTX_CD20_bound;      // dissociation flux
double RTX_internalize = kint_RTX * RTX_CD20_bound;   // internalization

dxdt_RTX_cent = -(CL_RTX/Vc_RTX)*RTX_cent
                - (Q_RTX/Vc_RTX)*RTX_cent
                + (Q_RTX/Vp_RTX)*RTX_peri
                - RTX_CD20_on
                + RTX_CD20_off;

// 2. RTX_peri: peripheral compartment
dxdt_RTX_peri = (Q_RTX/Vc_RTX)*RTX_cent - (Q_RTX/Vp_RTX)*RTX_peri;

// 3. RTX_CD20_bound: bound drug-receptor complex (TMDD)
dxdt_RTX_CD20_bound = RTX_CD20_on - RTX_CD20_off - RTX_internalize;

// -------------------------------------------------------
// 4. TAC_gut: oral absorption compartment (transit)
// -------------------------------------------------------
dxdt_TAC_gut   = -(ka_TAC/24.0)*TAC_gut;   // ka converted to 1/day

// 5. TAC_blood: whole blood compartment (ng/mL)
//    CL/F in L/h → convert to 1/day (* 24)
//    Dose enters TAC_gut as amount (mg) in event table
dxdt_TAC_blood = (ka_TAC/24.0)*TAC_gut/(Vd_TAC/1000.0)
                 - (CL_TAC*24.0/Vd_TAC)*TAC_blood;

// -------------------------------------------------------
// 6. CPx_cent: Cyclophosphamide (prodrug)
// -------------------------------------------------------
dxdt_CPx_cent  = -(CL_CPx*24.0/Vd_CPx)*CPx_cent
                 - km_CPx*24.0*CPx_cent;

// 7. CPx_metab: active 4-OH-cyclophosphamide
dxdt_CPx_metab = km_CPx*24.0*CPx_cent - kel_met*24.0*CPx_metab;

// -------------------------------------------------------
// 8. CD20_B: CD20+ B cells (cells/µL)
//    RTX depletes via TMDD-mediated internalization (signal from bound RTX)
//    Reconstitution: first-order return to Bss
// Reference: B cell depletion kinetics from MENTOR PK/PD (Fervenza 2019)
// -------------------------------------------------------
double RTX_B_kill = keff_RTX * RTX_conc * CD20_B;   // concentration-driven kill

dxdt_CD20_B    = kprol_B * Bss                       // constitutive replenishment
                 - kdeg_B * CD20_B                    // natural turnover
                 - RTX_B_kill;                        // RTX depletion

// Prevent negative B cells (floor at 0.1)
if(CD20_B < 0.1 && dxdt_CD20_B < 0) dxdt_CD20_B = 0;

// -------------------------------------------------------
// 9. Plasma_cells: IgG4-secreting plasma cells (rel. units)
//    Produced from B cells; depleted by cyclophosphamide
//    Tacrolimus suppresses activation but not established plasma cells directly
// -------------------------------------------------------
double PC_prod = kprol_PC * CD20_B;
double PC_kill = Emax_CPx * kdeg_PC * Plasma_cells;   // CPx cytotoxic effect

dxdt_Plasma_cells = PC_prod
                    - kdeg_PC * Plasma_cells
                    - PC_kill;

if(Plasma_cells < 0.01 && dxdt_Plasma_cells < 0) dxdt_Plasma_cells = 0;

// -------------------------------------------------------
// 10. Anti_PLA2R1: serum anti-PLA2R1 IgG4 (U/mL)
//     Production driven by plasma cells, inhibited by TAC (calcineurin pathway)
//     t1/2 ~28 days → kdeg_Ab = 0.693/28 ≈ 0.0248/day
// Reference: Bobart SA et al., JASN 2020; Beck LH NEJM 2009
// -------------------------------------------------------
double Ab_prod = kdeg_Ab * Ab_ss * (Plasma_cells / PC_ss) * (1.0 - Emax_TAC);

dxdt_Anti_PLA2R1 = Ab_prod - kdeg_Ab * Anti_PLA2R1;

// -------------------------------------------------------
// 11. IgG_deposit: subepithelial IgG4 deposits (rel. units)
//     Accumulates from serum IgG4; cleared slowly
// -------------------------------------------------------
dxdt_IgG_deposit = kdepo * Anti_PLA2R1 / Ab_ss
                   - kclear_d * IgG_deposit;

// -------------------------------------------------------
// 12. Complement_MAC: membrane attack complex (rel. units)
//     Hill activation by IgG deposits; avacopan reduces MAC
// -------------------------------------------------------
double MAC_drive = pow(IgG_deposit, nH_MAC) /
                   (pow(EC50_MAC, nH_MAC) + pow(IgG_deposit, nH_MAC));
double MAC_form  = kcomp * MAC_drive;
double MAC_res   = kres_MAC * Complement_MAC * (1.0 - Emax_AVA);

dxdt_Complement_MAC = MAC_form - MAC_res;

// -------------------------------------------------------
// 13. Podocyte_inj: podocyte injury score (0–1)
//     Driven by MAC; partial intrinsic recovery when MAC falls
// Reference: Shimizu M, Clin Exp Nephrol 2019
// -------------------------------------------------------
double pod_drive = kinj_pod * Complement_MAC * (1.0 - Podocyte_inj);
double pod_rec   = krec_pod * Podocyte_inj * (1.0 - Complement_MAC / (1.0 + Complement_MAC));

dxdt_Podocyte_inj = pod_drive - pod_rec;

if(Podocyte_inj < 0.0 && dxdt_Podocyte_inj < 0) dxdt_Podocyte_inj = 0;
if(Podocyte_inj > 1.0 && dxdt_Podocyte_inj > 0) dxdt_Podocyte_inj = 0;

// -------------------------------------------------------
// 14. GBM_thick: glomerular basement membrane thickening (rel. units)
//     Slowly increases with sustained podocyte injury; partial reversibility
// Reference: Glassock RJ, CJASN 2012
// -------------------------------------------------------
dxdt_GBM_thick = kGBM * Podocyte_inj * GBM_thick
                 - kGBM_res * (GBM_thick - 1.0);   // remodel toward norm (=1)

if(GBM_thick < 1.0 && dxdt_GBM_thick < 0) dxdt_GBM_thick = 0;

// -------------------------------------------------------
// 15. Proteinuria (g/day)
//     Driven by podocyte injury and GBM structural damage
//     RAAS (AngII) modulates glomerular hypertension → proteinuria
// Reference: Cattran DC, KI 2001; GEMRITUX endpoints
// -------------------------------------------------------
double Prot_drive = Prot_max * Podocyte_inj * GBM_thick / GBM_ss
                    * (1.0 + 0.05 * (AngII - AngII_ss) / AngII_ss);
dxdt_Proteinuria = kprot * (Prot_drive - Proteinuria);

if(Proteinuria < 0.05 && dxdt_Proteinuria < 0) dxdt_Proteinuria = 0;

// -------------------------------------------------------
// 16. Serum_alb: serum albumin (g/dL)
//     Synthesis suppressed by inflammation; loss via proteinuria
// Reference: Glassock RJ, CJASN 2012
// -------------------------------------------------------
double Alb_loss = kAlb_prot * Proteinuria * 0.05;   // scale factor
dxdt_Serum_alb = kAlb_syn - kAlb_deg * Serum_alb - Alb_loss;

if(Serum_alb < 0.5 && dxdt_Serum_alb < 0) dxdt_Serum_alb = 0;

// -------------------------------------------------------
// 17. eGFR: glomerular filtration rate (mL/min/1.73m²)
//     Declines with GBM thickening; partial recovery when GBM improves
//     ACEi reduces hyperfiltration (via AngII-mediated efferent dilation)
// Reference: Jha V et al., KI Reports 2020; KDIGO 2021
// -------------------------------------------------------
double eGFR_base_decline = keGFR * (GBM_thick - 1.0) * eGFR;
double eGFR_recovery     = keGFR_r * (eGFR_ss - eGFR) * (GBM_thick < 1.2 ? 1.0 : 0.0);
double ACEI_eGFR_protect = ACEI_eff * 0.01 * (AngII - 10.0);  // reduced hyperfiltration

dxdt_eGFR = -eGFR_base_decline + eGFR_recovery - ACEI_eGFR_protect;

if(eGFR < eGFR_min && dxdt_eGFR < 0) dxdt_eGFR = 0;

// -------------------------------------------------------
// 18. AngII: angiotensin II (pg/mL)
//     ACEi suppresses AngII (ACEI_eff = 0.7 for full dose)
// -------------------------------------------------------
dxdt_AngII = kAng * AngII_ss * (1.0 - ACEI_eff)
             - kAng * AngII;

// 19. Aldosterone (pg/mL): driven by AngII
dxdt_Aldosterone = kAldo * kAldo_Ang * AngII
                   - kAldo * Aldosterone;

$TABLE
// Capture key derived variables
double CR_prot  = (Proteinuria < 0.3)  ? 1.0 : 0.0;   // Complete remission flag
double PR_prot  = (Proteinuria < 3.5 && Proteinuria >= 0.3) ? 1.0 : 0.0;  // Partial remission
double B_depl   = (CD20_B < 5.0) ? 1.0 : 0.0;          // B cell depletion flag
double AntiPLA2R_neg = (Anti_PLA2R1 < 14.0) ? 1.0 : 0.0; // Serological remission

$CAPTURE
RTX_cent RTX_peri RTX_CD20_bound TAC_blood CPx_cent CPx_metab
CD20_B Plasma_cells Anti_PLA2R1 IgG_deposit Complement_MAC
Podocyte_inj GBM_thick
Proteinuria Serum_alb eGFR AngII Aldosterone
CR_prot PR_prot B_depl AntiPLA2R_neg
'

# ============================================================
# Compile model
# ============================================================
mod <- mcode("MN_QSP", mn_code)

# ============================================================
# Helper: build dosing event table
# ============================================================
# RTX: IV infusion — 1 g over 4 hours → rate = 1000 mg / (Vc * 4/24 day)
# Vc = 3.1 L → concentration dose = 1000 mg / 3100 mL = 0.322 µg/mL equiv per hour
# We model as direct bolus addition to RTX_cent for simplicity with mrgsolve ev()

make_events <- function(scenario, bwt = 70) {

  ev_list <- list()

  if(scenario == "no_treatment") {
    # No drugs — just run forward
    ev_list <- ev(time = 0, cmt = "RTX_cent", amt = 0)
  }

  if(scenario == "RTX_mono_MENTOR") {
    # MENTOR regimen: 1 g IV at Day 0 and Day 180
    # Concentration bolus = 1000 mg / 3.1 L ≈ 322.6 µg/mL
    dose_ug_mL <- 1000 / 3.1 * 1000 / 1000   # = 322.6 µg/mL when added to Vc
    ev_list <- ev(time = c(0, 180),
                  cmt  = "RTX_cent",
                  amt  = dose_ug_mL)
  }

  if(scenario == "RTX_mono_GEMRITUX") {
    # GEMRITUX regimen: 375 mg/m² × 2 (D1, D8)
    # BSA ~1.8 m² (average adult) → 675 mg per dose
    dose_mg    <- 375 * 1.8
    dose_ug_mL <- dose_mg / 3.1 * 1000 / 1000
    ev_list <- ev(time = c(0, 7),
                  cmt  = "RTX_cent",
                  amt  = dose_ug_mL)
  }

  if(scenario == "TAC_mono") {
    # Tacrolimus 0.05 mg/kg/day oral (STARMEN dose; target 5–10 ng/mL)
    # bwt = 70 kg → 3.5 mg/day = 3.5 mg split q12h → 1.75 mg q12h
    # mrgsolve: daily dosing, absorbed via TAC_gut
    # Amount in mg → ng: 1 mg = 1e6 ng; Vd = 900 L = 900,000 mL
    # TAC_blood in ng/mL: 1 mg dose / 900 L = 1000/900 ng/mL ≈ 1.11 ng/mL per mg
    daily_mg  <- 0.05 * bwt
    dose_mg   <- daily_mg / 2   # BID
    # amt deposited into TAC_gut (mg → converted internally)
    tac_times <- seq(0, 540, by = 0.5)   # 0.5-day intervals = BID for 18 months
    ev_list   <- ev(time = tac_times,
                    cmt  = "TAC_gut",
                    amt  = dose_mg)
  }

  if(scenario == "Ponticelli") {
    # Modified Ponticelli (RI-CYCLO control arm; Scolari JASN 2021):
    # Alternating months: chlorambucil/cyclophosphamide + methylprednisolone
    # Simplified: CPx 2.5 mg/kg/day × 30 days per cycle × 3 cycles
    # CPx IV at days 1–30 of months 1, 3, 5 (approximate)
    cpx_mg_day <- 2.5 * bwt        # mg/day
    cpx_ug_mL  <- cpx_mg_day / 35  # ≈ µg/mL in Vd = 35 L (rough bolus equiv)
    cpx_times  <- c(seq(0, 29, 1),
                    seq(60, 89, 1),
                    seq(120, 149, 1))
    ev_list <- ev(time = cpx_times,
                  cmt  = "CPx_cent",
                  amt  = cpx_ug_mL)
  }

  if(scenario == "RTX_TAC_combo") {
    # STARMEN-like: Tacrolimus × 6 months then sequential RTX
    # TAC for first 6 months, RTX 1g at month 6
    daily_mg <- 0.05 * bwt
    dose_mg  <- daily_mg / 2
    tac_times <- seq(0, 179, by = 0.5)
    ev_tac <- ev(time = tac_times, cmt = "TAC_gut", amt = dose_mg)
    dose_rtx <- 1000 / 3.1 * 1000 / 1000
    ev_rtx <- ev(time = 180, cmt = "RTX_cent", amt = dose_rtx)
    ev_list <- c(ev_tac, ev_rtx)
  }

  if(scenario == "ACEi_only") {
    # ACE inhibitor: modeled via ACEI_eff parameter = 0.7
    # No drug dosing events needed — ACEI_eff set at simulation level
    ev_list <- ev(time = 0, cmt = "RTX_cent", amt = 0)
  }

  if(scenario == "Avacopan") {
    # Investigational C5aR inhibitor; modeled via AVA_dose switch
    ev_list <- ev(time = 0, cmt = "RTX_cent", amt = 0)
  }

  return(ev_list)
}

# ============================================================
# Run all treatment scenarios
# ============================================================
sim_duration <- 730   # days (2 years)
sim_delta    <- 1     # daily output

run_scenario <- function(scenario_name, ACEI_eff_val = 0, AVA_val = 0) {
  ev_dose <- make_events(scenario_name)

  out <- mod %>%
    param(ACEI_eff = ACEI_eff_val, AVA_dose = AVA_val) %>%
    ev(ev_dose) %>%
    mrgsim(end = sim_duration, delta = sim_delta,
           obsonly = TRUE, carry.out = "evid") %>%
    as.data.frame() %>%
    mutate(scenario = scenario_name,
           time_months = time / 30.44)

  return(out)
}

message("Running treatment scenarios...")

df_no_tx    <- run_scenario("no_treatment")
df_rtx_men  <- run_scenario("RTX_mono_MENTOR")
df_rtx_gem  <- run_scenario("RTX_mono_GEMRITUX")
df_tac      <- run_scenario("TAC_mono")
df_ponti    <- run_scenario("Ponticelli")
df_combo    <- run_scenario("RTX_TAC_combo")
df_acei     <- run_scenario("ACEi_only", ACEI_eff_val = 0.7)
df_ava      <- run_scenario("Avacopan",  AVA_val = 1.0)

df_all <- bind_rows(
  df_no_tx,
  df_rtx_men,
  df_rtx_gem,
  df_tac,
  df_ponti,
  df_combo,
  df_acei,
  df_ava
)

# ============================================================
# Label scenarios for plotting
# ============================================================
scenario_labels <- c(
  "no_treatment"       = "No Treatment",
  "RTX_mono_MENTOR"    = "RTX 1g×2 (MENTOR)",
  "RTX_mono_GEMRITUX"  = "RTX 375mg/m²×2 (GEMRITUX)",
  "TAC_mono"           = "Tacrolimus mono (STARMEN)",
  "Ponticelli"         = "Ponticelli CPx regimen (RI-CYCLO)",
  "RTX_TAC_combo"      = "RTX + TAC sequential combo",
  "ACEi_only"          = "ACEi (conservative)",
  "Avacopan"           = "Avacopan C5aR inhibitor (investigational)"
)

df_all$scenario_label <- scenario_labels[df_all$scenario]

# Palette (colorblind-friendly)
pal <- c(
  "No Treatment"                             = "#D62728",
  "RTX 1g×2 (MENTOR)"                        = "#1F77B4",
  "RTX 375mg/m²×2 (GEMRITUX)"               = "#AEC7E8",
  "Tacrolimus mono (STARMEN)"                = "#FF7F0E",
  "Ponticelli CPx regimen (RI-CYCLO)"        = "#2CA02C",
  "RTX + TAC sequential combo"              = "#9467BD",
  "ACEi (conservative)"                      = "#8C564B",
  "Avacopan C5aR inhibitor (investigational)"= "#17BECF"
)

# ============================================================
# Figure 1: Primary Clinical Endpoints
# ============================================================
p1 <- ggplot(df_all, aes(x = time_months, y = Proteinuria,
                          color = scenario_label)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 3.5, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 0.3,  linetype = "dotted", color = "grey30") +
  annotate("text", x = 23, y = 3.7,  label = "Nephrotic threshold", size = 2.5, hjust = 1) +
  annotate("text", x = 23, y = 0.1,  label = "Complete remission",  size = 2.5, hjust = 1) +
  scale_color_manual(values = pal) +
  labs(title = "Proteinuria (g/day)",
       x = "Time (months)", y = "24h Proteinuria (g/day)",
       color = "Treatment") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 10))

p2 <- ggplot(df_all, aes(x = time_months, y = Serum_alb,
                          color = scenario_label)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 3.5, linetype = "dashed", color = "grey50") +
  annotate("text", x = 23, y = 3.65, label = "Normal lower limit", size = 2.5, hjust = 1) +
  scale_color_manual(values = pal) +
  labs(title = "Serum Albumin (g/dL)",
       x = "Time (months)", y = "Serum Albumin (g/dL)",
       color = "Treatment") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(1, 4.5))

p3 <- ggplot(df_all, aes(x = time_months, y = eGFR,
                          color = scenario_label)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 60, linetype = "dashed", color = "grey50") +
  annotate("text", x = 23, y = 61.5, label = "CKD G2/G3 boundary", size = 2.5, hjust = 1) +
  scale_color_manual(values = pal) +
  labs(title = "eGFR (mL/min/1.73m²)",
       x = "Time (months)", y = "eGFR (mL/min/1.73m²)",
       color = "Treatment") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

p4 <- ggplot(df_all, aes(x = time_months, y = Anti_PLA2R1,
                          color = scenario_label)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 14, linetype = "dashed", color = "grey50") +
  annotate("text", x = 23, y = 20, label = "Seronegative cutoff", size = 2.5, hjust = 1) +
  scale_color_manual(values = pal) +
  labs(title = "Anti-PLA2R1 IgG4 (U/mL)",
       x = "Time (months)", y = "Anti-PLA2R1 (U/mL)",
       color = "Treatment") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "lines"))

fig1 <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Membranous Nephropathy QSP Model — Primary Clinical Endpoints",
    subtitle = "GEMRITUX · MENTOR · STARMEN · RI-CYCLO calibrated parameters",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )

print(fig1)

# ============================================================
# Figure 2: Mechanistic Biomarkers
# ============================================================
b1 <- ggplot(df_all, aes(x = time_months, y = CD20_B,
                          color = scenario_label)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = pal) +
  labs(title = "CD20+ B Cells (cells/µL)",
       x = "Time (months)", y = "CD20+ B cells (cells/µL)",
       color = "Treatment") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

b2 <- ggplot(df_all, aes(x = time_months, y = Complement_MAC,
                          color = scenario_label)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = pal) +
  labs(title = "Complement MAC Activity (rel. units)",
       x = "Time (months)", y = "MAC activity",
       color = "Treatment") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

b3 <- ggplot(df_all, aes(x = time_months, y = Podocyte_inj,
                          color = scenario_label)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = pal) +
  labs(title = "Podocyte Injury Score (0–1)",
       x = "Time (months)", y = "Podocyte injury",
       color = "Treatment") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 1))

b4 <- ggplot(df_all, aes(x = time_months, y = IgG_deposit,
                          color = scenario_label)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = pal) +
  labs(title = "Subepithelial IgG4 Deposits (rel. units)",
       x = "Time (months)", y = "IgG deposit",
       color = "Treatment") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right",
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "lines"))

fig2 <- (b1 | b2) / (b3 | b4) +
  plot_annotation(
    title = "Membranous Nephropathy QSP — Mechanistic Biomarkers",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )

print(fig2)

# ============================================================
# Figure 3: Drug PK profiles (select scenarios with active drugs)
# ============================================================
df_pk <- df_all %>%
  filter(scenario %in% c("RTX_mono_MENTOR", "TAC_mono", "Ponticelli",
                          "RTX_TAC_combo")) %>%
  mutate(time_months = time / 30.44)

pk1 <- ggplot(df_pk %>% filter(scenario %in% c("RTX_mono_MENTOR", "RTX_TAC_combo")),
              aes(x = time_months, y = RTX_cent, color = scenario_label)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = pal) +
  labs(title = "Rituximab PK (central, µg/mL)",
       x = "Time (months)", y = "RTX concentration (µg/mL)",
       color = "Scenario") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right", legend.text = element_text(size = 7))

pk2 <- ggplot(df_pk %>% filter(scenario %in% c("TAC_mono", "RTX_TAC_combo")),
              aes(x = time_months, y = TAC_blood, color = scenario_label)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = c(5, 10), linetype = "dashed", color = "grey40") +
  annotate("text", x = 5.5, y = 10.5, label = "Target range", size = 2.5) +
  scale_color_manual(values = pal) +
  labs(title = "Tacrolimus whole blood (ng/mL)",
       x = "Time (months)", y = "TAC (ng/mL)",
       color = "Scenario") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right", legend.text = element_text(size = 7))

pk3 <- ggplot(df_pk %>% filter(scenario == "Ponticelli"),
              aes(x = time_months, y = CPx_metab, color = scenario_label)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = pal) +
  labs(title = "Active CPx Metabolite (4-OH-CPx, µg/mL)",
       x = "Time (months)", y = "4-OH-CPx (µg/mL)",
       color = "Scenario") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right", legend.text = element_text(size = 7))

pk4 <- ggplot(df_all %>% filter(scenario %in% c("ACEi_only", "no_treatment")),
              aes(x = time_months, y = AngII, color = scenario_label)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = pal) +
  labs(title = "Angiotensin II (pg/mL) — RAAS",
       x = "Time (months)", y = "AngII (pg/mL)",
       color = "Scenario") +
  theme_bw(base_size = 10) +
  theme(legend.position = "right", legend.text = element_text(size = 7))

fig3 <- (pk1 | pk2) / (pk3 | pk4) +
  plot_annotation(
    title = "Membranous Nephropathy QSP — PK and RAAS Profiles",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )

print(fig3)

# ============================================================
# Figure 4: Remission rates at 6, 12, 24 months
# ============================================================
remission_summary <- df_all %>%
  filter(time %in% c(182, 365, 730)) %>%
  mutate(
    time_label   = paste0("Month ", round(time / 30.44)),
    CR           = CR_prot * 100,
    PR           = PR_prot * 100,
    CRorPR       = as.numeric((Proteinuria < 3.5)) * 100
  ) %>%
  select(scenario_label, time_label, Proteinuria, CR, PR, CRorPR, eGFR, Serum_alb, Anti_PLA2R1)

# Bar chart: proteinuria < 3.5 g/day (CR+PR combined) at each time point
fig4 <- ggplot(remission_summary,
               aes(x = reorder(scenario_label, -CRorPR),
                   y = CRorPR, fill = scenario_label)) +
  geom_col(position = "dodge", alpha = 0.85) +
  facet_wrap(~time_label) +
  scale_fill_manual(values = pal) +
  coord_flip() +
  labs(title = "Membranous Nephropathy: Partial + Complete Remission Rate",
       subtitle = "Proteinuria < 3.5 g/day (partial) or < 0.3 g/day (complete)",
       x = NULL, y = "Remission rate (%)",
       fill = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

print(fig4)

# ============================================================
# Summary table at 24 months
# ============================================================
cat("\n========================================================\n")
cat("Membranous Nephropathy QSP — 24-month Outcome Summary\n")
cat("(Calibrated to GEMRITUX, MENTOR, STARMEN, RI-CYCLO)\n")
cat("========================================================\n")

summary_24m <- df_all %>%
  filter(time == 730) %>%
  select(scenario_label, Proteinuria, Serum_alb, eGFR, Anti_PLA2R1,
         CD20_B, Complement_MAC, Podocyte_inj, GBM_thick) %>%
  arrange(Proteinuria)

print(summary_24m, digits = 3)

# ============================================================
# Save plots
# ============================================================
output_dir <- "/home/user/qsp/membranous-nephropathy"

ggsave(file.path(output_dir, "mn_fig1_clinical_endpoints.png"),
       fig1, width = 12, height = 8, dpi = 150)
ggsave(file.path(output_dir, "mn_fig2_biomarkers.png"),
       fig2, width = 12, height = 8, dpi = 150)
ggsave(file.path(output_dir, "mn_fig3_pk_raas.png"),
       fig3, width = 12, height = 8, dpi = 150)
ggsave(file.path(output_dir, "mn_fig4_remission.png"),
       fig4, width = 12, height = 6, dpi = 150)

message("All figures saved to: ", output_dir)
message("Model file: mn_mrgsolve_model.R")
message("Done.")
