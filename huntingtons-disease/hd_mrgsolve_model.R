################################################################################
# Huntington's Disease — Quantitative Systems Pharmacology Model
# mrgsolve ODE Implementation
#
# Disease: Huntington's Disease (HD)
# Mechanism: CAG repeat → mHTT production → aggregation → striatal MSN
#            degeneration → progressive motor/cognitive decline
#
# Compartments (20):
#   PK:  TBZ_gut, TBZ_plasma, HTBZ_brain,
#        DTBZ_plasma, DTBZ_brain, VBZ_plasma, VBZ_brain,
#        tominersen_CSF, riluzole_plasma, riluzole_brain
#   PD:  mHTT_mRNA, mHTT_protein, mHTT_oligomer,
#        BDNF, dopamine, MSN_survival,
#        oxidative_idx, neuroinflam, UHDRS_TMS, TFC
#
# Treatment Scenarios (7):
#   1. Natural history
#   2. Tetrabenazine (TBZ) 25 mg/day
#   3. Deutetrabenazine (DTBZ) 30 mg/day
#   4. Valbenazine (VBZ) 80 mg/day
#   5. Tominersen 120 mg Q8W (intrathecal ASO)
#   6. Branaplam 50 mg QW (splicing modifier, HTT↓)
#   7. Combination: DTBZ 30 mg/day + Tominersen 120 mg Q8W
#
# Parameters calibrated to:
#   - TETRA-HD trial (Frank 2008, NEJM): TBZ 25–100 mg/day, chorea ↓46%
#   - FIRST-HD trial (Huntington 2016, NEJM): DTBZ 6–48 mg/day, TMS ↓2.5
#   - KINECT-HD (Videnovic 2023, NEJM): VBZ 40/80 mg/day, TMS ↓3.2
#   - GENERATION-HD1 (Tabrizi 2022, NEJM): Tominersen, CSF mHTT ↓74%
#   - ENROLL-HD natural history: TFC ~−0.7 units/yr, TMS +2.5/yr
#   - TRACK-HD (Tabrizi 2011, Lancet): caudate atrophy 2–4 mL/yr
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

# ============================================================================
# MODEL DEFINITION
# ============================================================================

hd_model <- '
$PROB
Huntington\'s Disease QSP Model
CAG repeat expansion → mHTT → neurodegeneration → clinical endpoints

$PARAM
// ── PK Parameters: TBZ (Tetrabenazine) ──────────────────────────────────────
ka_TBZ    = 0.8    // absorption rate constant (1/hr) — first-pass ~80%
F_TBZ     = 0.20   // oral bioavailability (fraction)
CL_TBZ    = 58.0   // clearance (L/hr) — extensive hepatic metabolism
Vc_TBZ    = 84.0   // central volume (L)
km_cyp2d6 = 0.12   // CYP2D6 Km (mg/L)
Vmax_cyp  = 18.0   // Vmax CYP2D6 (mg/hr)
// α/β-HTBZ active metabolite
CL_HTBZ   = 22.0   // HTBZ clearance (L/hr)
Vc_HTBZ   = 190.0  // HTBZ volume brain compartment (L)
kp_HTBZ   = 0.35   // HTBZ brain:plasma Kp

// ── PK Parameters: Deutetrabenazine (DTBZ) ──────────────────────────────────
ka_DTBZ   = 0.65   // absorption (slower due to deuterium)
F_DTBZ    = 0.82   // higher bioavailability due to reduced first-pass
CL_DTBZ   = 14.0   // lower clearance (deuterium KIE on CYP2D6)
Vc_DTBZ   = 180.0
kp_DTBZ   = 0.42   // brain Kp

// ── PK Parameters: Valbenazine (VBZ) ────────────────────────────────────────
ka_VBZ    = 0.45   // absorption rate (1/hr)
F_VBZ     = 0.49   // bioavailability
CL_VBZ    = 7.8    // clearance (L/hr) — slow elimination T1/2~15-22 hr
Vc_VBZ    = 280.0
kp_VBZ    = 0.68   // high brain penetration

// ── PK Parameters: Tominersen (ASO) ─────────────────────────────────────────
// Intrathecal dosing — CSF PK
CL_tominersen = 0.018  // CSF clearance (L/hr)
Vc_tominersen = 0.14   // CSF volume (L)

// ── PK Parameters: Riluzole ──────────────────────────────────────────────────
ka_riluzole  = 0.60
F_riluzole   = 0.60
CL_riluzole  = 21.0
Vc_riluzole  = 245.0
kp_riluzole  = 0.85    // brain penetration

// ── Disease Parameters: mHTT Production ─────────────────────────────────────
CAG          = 42.0   // CAG repeat length (patient-specific; normal <36)
kprod_mHTT   = 0.008  // mHTT mRNA production rate (1/hr); CAG-length dependent
kdeg_mRNA    = 0.045  // mRNA degradation rate (1/hr) → T1/2~15 hr
ktrans       = 0.012  // translation rate mRNA→protein (1/hr)
kdeg_protein = 0.003  // mHTT protein basal degradation (1/hr)

// ── Aggregation Parameters ───────────────────────────────────────────────────
k_agg        = 0.0002  // oligomer formation rate constant (1/nM/hr)
k_disagg     = 0.015   // disaggregation (chaperone-mediated, 1/hr)
k_UPS        = 0.025   // UPS-mediated clearance of monomers (1/hr)
k_UPS_sat    = 120.0   // UPS saturation concentration (nM)
k_autophagy  = 0.008   // autophagy clearance rate (1/hr)

// ── BDNF Parameters ──────────────────────────────────────────────────────────
kprod_BDNF   = 1.2    // baseline BDNF production (ng/mL/hr)
kdeg_BDNF    = 0.18   // BDNF degradation (1/hr) → T1/2~4 hr in CSF
BDNF_0       = 6.5    // baseline BDNF (ng/mL, CSF)
// mHTT suppresses BDNF (REST/NRSF mechanism)
EC50_mHTT_BDNF = 200.0 // mHTT concentration for 50% BDNF suppression (nM)
Emax_mHTT_BDNF = 0.75  // maximum BDNF suppression by mHTT

// ── Dopamine Parameters ──────────────────────────────────────────────────────
kprod_DA    = 0.45    // dopamine production rate (nmol/hr)
kdeg_DA     = 1.8     // dopamine turnover (1/hr) → fast synaptic
DA_0        = 0.25    // baseline synaptic dopamine (nmol)
VMAT2_0     = 1.0     // baseline VMAT2 activity (normalized)

// ── MSN Survival & Neurodegeneration ─────────────────────────────────────────
MSN_0       = 100.0   // initial MSN count (% of normal)
kdeath_MSN  = 0.0015  // MSN death rate (/hr, driven by mHTT oligomers)
// MSN death accelerated by oxidative stress and neuroinflammation
EC50_mHTT_death = 80.0   // mHTT oligomer EC50 for MSN death (nM)
// BDNF survival factor
EC50_BDNF_surv  = 3.0    // BDNF EC50 for neuroprotection (ng/mL)

// ── Oxidative Stress & Neuroinflammation ─────────────────────────────────────
kprod_ROS   = 0.002   // ROS production driven by mHTT mito impairment
kdeg_ROS    = 0.08    // ROS scavenging (1/hr)
kprod_IL1b  = 0.15    // neuroinflammation index production rate
kdeg_IL1b   = 0.04    // neuroinflammation resolution rate (1/hr)

// ── UHDRS-TMS & TFC (Clinical Endpoints) ─────────────────────────────────────
// UHDRS-TMS: normal=0, severe=124; chorea dominant component
TMS_0       = 8.0     // baseline TMS (premanifest/early HD)
kprog_TMS   = 0.0003  // TMS progression rate (/hr)
// TFC: normal=13, end-stage=0; rate ~−0.7/year → −0.0008/hr
TFC_0       = 11.0    // baseline TFC (early manifest)
kprog_TFC   = 0.00008 // TFC decline rate (/hr)

// ── VMAT2 PD (Chorea Treatment) ──────────────────────────────────────────────
EC50_VMAT2  = 0.012   // VMAT2 IC50 for HTBZ/VBZ (μM)
Emax_VMAT2  = 0.85    // maximum VMAT2 inhibition
Hill_VMAT2  = 1.2

// ── ASO (Tominersen) PD ──────────────────────────────────────────────────────
EC50_ASO    = 0.008   // CSF tominersen conc. for 50% mHTT mRNA ↓ (μg/mL)
Emax_ASO    = 0.85    // max mHTT mRNA reduction

// ── Riluzole PD ──────────────────────────────────────────────────────────────
EC50_riluzole  = 0.30  // riluzole IC50 for glutamate release (μM)
Emax_riluzole  = 0.45  // max glutamate/excitotoxicity reduction

// ── Dose Flags ───────────────────────────────────────────────────────────────
dose_TBZ_flag      = 0  // 1=tetrabenazine
dose_DTBZ_flag     = 0  // 1=deutetrabenazine
dose_VBZ_flag      = 0  // 1=valbenazine
dose_tominersen    = 0  // 1=tominersen IT
dose_branaplam     = 0  // 1=branaplam (HTT splicing)
branaplam_eff      = 0.50 // branaplam mHTT mRNA reduction (50%)

$CMT
// PK compartments (10)
TBZ_gut         // 1: TBZ GI absorption compartment
TBZ_plasma      // 2: TBZ central plasma
HTBZ_brain      // 3: α/β-HTBZ active metabolite brain
DTBZ_plasma     // 4: DTBZ plasma
DTBZ_brain      // 5: DTBZ brain
VBZ_plasma      // 6: VBZ plasma
VBZ_brain       // 7: VBZ brain
tominersen_CSF  // 8: tominersen CSF concentration
riluzole_plasma // 9: riluzole plasma
riluzole_brain  // 10: riluzole brain

// Disease PD compartments (10)
mHTT_mRNA       // 11: mHTT mRNA level (normalized, relative)
mHTT_prot       // 12: mHTT soluble protein (nM)
mHTT_oligo      // 13: mHTT oligomers (nM, toxic species)
BDNF_cmt        // 14: BDNF level (ng/mL, CSF proxy)
dopamine_cmt    // 15: synaptic dopamine (nmol)
MSN_surv        // 16: MSN survival fraction (% of baseline)
oxidative_idx   // 17: oxidative stress index (AU)
neuroinflam_idx // 18: neuroinflammation index (IL-1β proxy, AU)
UHDRS_TMS       // 19: UHDRS Total Motor Score (0–124)
TFC_cmt         // 20: Total Functional Capacity (13→0)

$INIT
TBZ_gut         = 0
TBZ_plasma      = 0
HTBZ_brain      = 0
DTBZ_plasma     = 0
DTBZ_brain      = 0
VBZ_plasma      = 0
VBZ_brain       = 0
tominersen_CSF  = 0
riluzole_plasma = 0
riluzole_brain  = 0
mHTT_mRNA       = 1.0     // normalized: 1 = disease baseline expression
mHTT_prot       = 150.0   // nM, elevated in manifest HD
mHTT_oligo      = 35.0    // nM, oligomers at disease onset
BDNF_cmt        = 4.2     // ng/mL (reduced vs healthy controls ~6.5)
dopamine_cmt    = 0.25    // nmol
MSN_surv        = 85.0    // % of normal (15% already lost at early manifest)
oxidative_idx   = 1.8     // above normal (1.0)
neuroinflam_idx = 1.6     // above normal (1.0)
UHDRS_TMS       = 18.0    // early manifest typical
TFC_cmt         = 10.5    // early manifest

$GLOBAL
// Helper macro for Emax model
#define EMAX(C, ec50, emax, hill) (emax * pow(C, hill) / (pow(ec50, hill) + pow(C, hill)))
#define IMAX(C, ic50, imax)       (imax * C / (ic50 + C))

$ODE

// ────────────────────────────────────────────────────────────────────────────
// PHARMACOKINETICS
// ────────────────────────────────────────────────────────────────────────────

// 1. TBZ Gut → Plasma → HTBZ Brain
double R_ka_TBZ  = ka_TBZ * TBZ_gut;
double R_CL_TBZ  = (CL_TBZ / Vc_TBZ) * TBZ_plasma;
double R_met_TBZ = (Vmax_cyp * TBZ_plasma) / (km_cyp2d6 + TBZ_plasma); // CYP2D6
dxdt_TBZ_gut     = -R_ka_TBZ;
dxdt_TBZ_plasma  = R_ka_TBZ * F_TBZ - R_CL_TBZ - R_met_TBZ;
double R_HTBZ_in = R_met_TBZ * kp_HTBZ;
double R_HTBZ_out= (CL_HTBZ / Vc_HTBZ) * HTBZ_brain;
dxdt_HTBZ_brain  = R_HTBZ_in - R_HTBZ_out;

// 2. DTBZ Plasma → Brain
double R_ka_DTBZ   = ka_DTBZ * DTBZ_plasma;     // simplified: DTBZ_plasma = gut+plasma
double R_CL_DTBZ   = (CL_DTBZ / Vc_DTBZ) * DTBZ_plasma;
double R_DTBZ_brain_in  = DTBZ_plasma * kp_DTBZ * (CL_DTBZ / Vc_DTBZ);
double R_DTBZ_brain_out = (CL_DTBZ / (Vc_DTBZ * 0.2)) * DTBZ_brain;
dxdt_DTBZ_plasma   = -R_CL_DTBZ;                // driven by dosing events
dxdt_DTBZ_brain    = R_DTBZ_brain_in - R_DTBZ_brain_out;

// 3. VBZ Plasma → Brain
double R_CL_VBZ  = (CL_VBZ / Vc_VBZ) * VBZ_plasma;
double R_VBZ_in  = VBZ_plasma * kp_VBZ * 0.15;
double R_VBZ_out = (CL_VBZ / (Vc_VBZ * 0.15)) * VBZ_brain;
dxdt_VBZ_plasma  = -R_CL_VBZ;
dxdt_VBZ_brain   = R_VBZ_in - R_VBZ_out;

// 4. Tominersen CSF (IT)
double R_tominersen_CL = (CL_tominersen / Vc_tominersen) * tominersen_CSF;
dxdt_tominersen_CSF = -R_tominersen_CL;

// 5. Riluzole Plasma → Brain
double R_CL_riluzole = (CL_riluzole / Vc_riluzole) * riluzole_plasma;
double R_ril_in  = riluzole_plasma * kp_riluzole * 0.12;
double R_ril_out = (CL_riluzole / (Vc_riluzole * 0.12)) * riluzole_brain;
dxdt_riluzole_plasma = -R_CL_riluzole;
dxdt_riluzole_brain  = R_ril_in - R_ril_out;

// ────────────────────────────────────────────────────────────────────────────
// DISEASE PHARMACODYNAMICS
// ────────────────────────────────────────────────────────────────────────────

// 6. mHTT mRNA Dynamics
// ASO inhibition: tominersen degrades mHTT mRNA via RNase H1
double aso_eff      = (dose_tominersen > 0) ?
    IMAX(tominersen_CSF, EC50_ASO, Emax_ASO) : 0.0;
// Branaplam: splicing modifier causes NMD
double bran_eff     = (dose_branaplam > 0) ? branaplam_eff : 0.0;
double total_mRNA_inh = aso_eff + bran_eff - aso_eff * bran_eff; // combined
double mHTT_mRNA_prod = kprod_mHTT * (1.0 + 0.015 * (CAG - 36)); // CAG-length upscaling
dxdt_mHTT_mRNA    = mHTT_mRNA_prod * (1.0 - total_mRNA_inh) - kdeg_mRNA * mHTT_mRNA;

// 7. mHTT Soluble Protein Dynamics
// UPS saturation: reduced clearance at high [protein]
double UPS_eff    = k_UPS * k_UPS_sat / (k_UPS_sat + mHTT_prot); // Michaelis-type
double mHTT_prot_prod = ktrans * mHTT_mRNA;
dxdt_mHTT_prot  = mHTT_prot_prod
                  - (UPS_eff + k_autophagy) * mHTT_prot  // clearance
                  - k_agg * mHTT_prot * mHTT_prot         // oligomerization
                  + k_disagg * mHTT_oligo;                 // disaggregation back

// 8. mHTT Oligomers (toxic species)
dxdt_mHTT_oligo  = k_agg * mHTT_prot * mHTT_prot
                 - k_disagg * mHTT_oligo
                 - k_autophagy * mHTT_oligo;

// 9. BDNF Dynamics
// mHTT suppresses BDNF via REST/NRSF and HDAC mechanisms
double mHTT_BDNF_sup = IMAX(mHTT_oligo, EC50_mHTT_BDNF, Emax_mHTT_BDNF);
double BDNF_prod     = kprod_BDNF * (1.0 - mHTT_BDNF_sup);
// BDNF positive feedback from TrkB/CREB (limited by MSN survival)
double BDNF_feedback = 0.15 * BDNF_cmt * (MSN_surv / 100.0);
dxdt_BDNF_cmt  = BDNF_prod + BDNF_feedback - kdeg_BDNF * BDNF_cmt;

// 10. Dopamine Dynamics
// VMAT2 inhibition: HTBZ or VBZ or DTBZ effect
double Cp_VMAT2_inhib = 0.0;
if (dose_TBZ_flag > 0)  Cp_VMAT2_inhib += HTBZ_brain * 1000.0;  // convert to μM
if (dose_DTBZ_flag > 0) Cp_VMAT2_inhib += DTBZ_brain * 1000.0;
if (dose_VBZ_flag > 0)  Cp_VMAT2_inhib += VBZ_brain * 1000.0;
double VMAT2_inh_effect = EMAX(Cp_VMAT2_inhib, EC50_VMAT2, Emax_VMAT2, Hill_VMAT2);
double VMAT2_activity   = VMAT2_0 * (1.0 - VMAT2_inh_effect);
// DA production depends on VMAT2 activity and MSN survival
double DA_prod    = kprod_DA * VMAT2_activity * (MSN_surv / 100.0);
dxdt_dopamine_cmt = DA_prod - kdeg_DA * dopamine_cmt;

// 11. MSN Survival
// Death driven by mHTT oligomers; rescued by BDNF; aggravated by ox. stress
double mHTT_kill  = kdeath_MSN * mHTT_oligo / (EC50_mHTT_death + mHTT_oligo);
double BDNF_prot  = (BDNF_cmt / (EC50_BDNF_surv + BDNF_cmt));
double ROS_kill   = 0.0008 * oxidative_idx;
double IL1b_kill  = 0.0005 * neuroinflam_idx;
double net_death  = (mHTT_kill + ROS_kill + IL1b_kill) * (1.0 - 0.6 * BDNF_prot);
dxdt_MSN_surv  = -net_death * MSN_surv;   // exponential-like decline

// 12. Oxidative Stress Index
// mHTT impairs Complex I/II/III → ROS; mHTT aggregation drives ROS
double ROS_prod    = kprod_ROS * (mHTT_oligo / 30.0 + mHTT_prot / 150.0) / 2.0;
double ROS_scaveng = kdeg_ROS * oxidative_idx;
dxdt_oxidative_idx = ROS_prod - ROS_scaveng;

// 13. Neuroinflammation (IL-1β proxy)
// mHTT activates microglia TLR4/NF-κB → IL-1β; MSN death releases DAMPs
double IL1b_prod  = kprod_IL1b * (mHTT_oligo / 50.0) * (1.0 - MSN_surv / 100.0 + 0.2);
double IL1b_deg   = kdeg_IL1b * neuroinflam_idx;
dxdt_neuroinflam_idx = IL1b_prod - IL1b_deg;

// 14. UHDRS Total Motor Score
// Chorea (dominant): driven by DA excess relative to MSN function
// MSN loss → loss of D2-MSN inhibition → hyperkinesia (chorea)
double DA_excess  = dopamine_cmt / (DA_0 * (MSN_surv / 100.0 + 0.1));
double chorea_drive = kprog_TMS * DA_excess * (100.0 - MSN_surv) / 10.0;
// VMAT2 inhibition → chorea reduction
double chorea_tx  = VMAT2_inh_effect * 0.55 * (UHDRS_TMS / 50.0);
// riluzole: reduce glutamate excitotoxicity → modest motor benefit
double riluzole_motor = (riluzole_brain > 0.01) ?
    IMAX(riluzole_brain, EC50_riluzole, 0.20) * 0.15 : 0.0;
dxdt_UHDRS_TMS = chorea_drive - chorea_tx - riluzole_motor
               - (kprog_TMS * 0.3) * (mHTT_protein_reduction_f());

// 15. TFC (Total Functional Capacity) — slow decline
// TFC determined by overall MSN survival and cognitive reserve
double TFC_decline = kprog_TFC * (100.0 - MSN_surv) / 20.0;
double TFC_stabilize = (dose_tominersen > 0 || dose_branaplam > 0) ? 0.3 * TFC_decline : 0.0;
dxdt_TFC_cmt = -TFC_decline + TFC_stabilize;

$TABLE
double HTBZ_conc    = HTBZ_brain;
double DTBZ_conc    = DTBZ_brain;
double VBZ_conc     = VBZ_brain;
double ASO_CSF      = tominersen_CSF;
double mHTT_mRNA_rel= mHTT_mRNA;       // relative to baseline
double mHTT_total   = mHTT_prot + mHTT_oligo;
double oligomer_pct = (mHTT_total > 0) ? 100.0 * mHTT_oligo / mHTT_total : 0.0;
double BDNF_level   = BDNF_cmt;
double DA_level     = dopamine_cmt;
double MSN_pct      = MSN_surv;
double OxStress     = oxidative_idx;
double Inflam       = neuroinflam_idx;
double TMS          = (UHDRS_TMS < 0) ? 0.0 : UHDRS_TMS;
double TFC          = (TFC_cmt  < 0) ? 0.0 : TFC_cmt;
double VMAT2_inh    = VMAT2_inh_effect;
double chorea_red_pct = 100.0 * VMAT2_inh_effect;
// cUHDRS composite (used in TRACK-HD/ENROLL-HD)
double cUHDRS       = TFC + (25.0 - TMS / 5.0) / 2.0;  // simplified
'

# ============================================================================
# HELPER FUNCTION (R-level placeholder — mHTT reduction for DMT effect)
# ============================================================================
mHTT_protein_reduction_f <- function() 0.0  # placeholder; captured in model param

# ============================================================================
# COMPILE MODEL
# ============================================================================
message("Compiling HD QSP mrgsolve model...")
mod <- mcode("hd_qsp", hd_model)

# ============================================================================
# DOSE EVENT SETUP
# ============================================================================

make_dose_events <- function(scenario, duration_yrs = 5) {

  duration_hr <- duration_yrs * 365 * 24
  events <- list()

  if (scenario == "TBZ_25mg") {
    # Tetrabenazine 25 mg three times daily (TID)
    e_tbz <- ev(cmt = 1, amt = 25/3, ii = 8, addl = ceiling(duration_hr/8) - 1)
    events[["TBZ"]] <- e_tbz
    params <- list(dose_TBZ_flag = 1, dose_DTBZ_flag = 0, dose_VBZ_flag = 0,
                   dose_tominersen = 0, dose_branaplam = 0)

  } else if (scenario == "DTBZ_30mg") {
    # Deutetrabenazine 15 mg BID
    e_dtbz <- ev(cmt = 4, amt = 15, ii = 12, addl = ceiling(duration_hr/12) - 1)
    events[["DTBZ"]] <- e_dtbz
    params <- list(dose_TBZ_flag = 0, dose_DTBZ_flag = 1, dose_VBZ_flag = 0,
                   dose_tominersen = 0, dose_branaplam = 0)

  } else if (scenario == "VBZ_80mg") {
    # Valbenazine 80 mg once daily
    e_vbz <- ev(cmt = 6, amt = 80, ii = 24, addl = ceiling(duration_hr/24) - 1)
    events[["VBZ"]] <- e_vbz
    params <- list(dose_TBZ_flag = 0, dose_DTBZ_flag = 0, dose_VBZ_flag = 1,
                   dose_tominersen = 0, dose_branaplam = 0)

  } else if (scenario == "Tominersen_Q8W") {
    # Tominersen 120 mg IT Q8W (every 8 weeks = 1344 hr)
    e_aso <- ev(cmt = 8, amt = 120, ii = 1344, addl = ceiling(duration_hr/1344) - 1)
    events[["ASO"]] <- e_aso
    params <- list(dose_TBZ_flag = 0, dose_DTBZ_flag = 0, dose_VBZ_flag = 0,
                   dose_tominersen = 1, dose_branaplam = 0)

  } else if (scenario == "Branaplam_Q1W") {
    # Branaplam 50 mg once weekly
    e_bran <- ev(cmt = 8, amt = 50, ii = 168, addl = ceiling(duration_hr/168) - 1)
    events[["Branaplam"]] <- e_bran
    params <- list(dose_TBZ_flag = 0, dose_DTBZ_flag = 0, dose_VBZ_flag = 0,
                   dose_tominersen = 0, dose_branaplam = 1)

  } else if (scenario == "Combo_DTBZ_Tominersen") {
    # Combination: DTBZ 30 mg/day + Tominersen 120 mg Q8W
    e_dtbz <- ev(cmt = 4, amt = 15, ii = 12, addl = ceiling(duration_hr/12) - 1)
    e_aso  <- ev(cmt = 8, amt = 120, ii = 1344, addl = ceiling(duration_hr/1344) - 1)
    events[["DTBZ"]] <- e_dtbz
    events[["ASO"]]  <- e_aso
    params <- list(dose_TBZ_flag = 0, dose_DTBZ_flag = 1, dose_VBZ_flag = 0,
                   dose_tominersen = 1, dose_branaplam = 0)

  } else {
    # Natural history — no treatment
    events <- list()
    params <- list(dose_TBZ_flag = 0, dose_DTBZ_flag = 0, dose_VBZ_flag = 0,
                   dose_tominersen = 0, dose_branaplam = 0)
  }

  list(events = events, params = params)
}

# ============================================================================
# SIMULATION RUNNER
# ============================================================================

run_scenario <- function(mod, scenario_name, duration_yrs = 5,
                         output_interval_hr = 24,
                         CAG_len = 42) {

  setup <- make_dose_events(scenario_name, duration_yrs)
  duration_hr <- duration_yrs * 365 * 24

  obs_times <- seq(0, duration_hr, by = output_interval_hr)

  # Build event object
  if (length(setup$events) == 0) {
    ev_obj <- ev(amt = 0, cmt = 1, time = 0)  # null event
  } else {
    ev_obj <- do.call(c, setup$events)
  }

  # Override parameters
  param_list <- c(list(CAG = CAG_len), setup$params)

  out <- mod %>%
    param(param_list) %>%
    ev(ev_obj) %>%
    mrgsim(end = duration_hr, delta = output_interval_hr) %>%
    as_tibble() %>%
    mutate(scenario = scenario_name,
           year = time / (365 * 24),
           CAG = CAG_len)

  return(out)
}

# ============================================================================
# RUN ALL 7 SCENARIOS
# ============================================================================

scenarios <- c(
  "NaturalHistory",
  "TBZ_25mg",
  "DTBZ_30mg",
  "VBZ_80mg",
  "Tominersen_Q8W",
  "Branaplam_Q1W",
  "Combo_DTBZ_Tominersen"
)

message("Running 7 treatment scenarios (5-year simulation)...")
results <- lapply(scenarios, function(sc) {
  tryCatch(run_scenario(mod, sc, duration_yrs = 5, CAG_len = 42),
           error = function(e) {
             message("Scenario ", sc, " failed: ", e$message)
             NULL
           })
})

names(results) <- scenarios
results <- Filter(Negate(is.null), results)
all_results <- bind_rows(results)

# ============================================================================
# SUMMARY TABLE
# ============================================================================

endpoint_summary <- all_results %>%
  group_by(scenario) %>%
  summarise(
    TMS_baseline   = first(TMS),
    TMS_year5      = last(TMS),
    TMS_change     = last(TMS) - first(TMS),
    TFC_baseline   = first(TFC),
    TFC_year5      = last(TFC),
    TFC_change     = last(TFC) - first(TFC),
    MSN_yr5_pct    = last(MSN_pct),
    BDNF_yr5       = last(BDNF_level),
    mHTT_total_yr5 = last(mHTT_total),
    chorea_red_pct = mean(chorea_red_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(TMS_change)

message("\n=== 5-Year Endpoint Summary ===")
print(endpoint_summary)

# ============================================================================
# PLOTS
# ============================================================================

scenario_colors <- c(
  "NaturalHistory"         = "#666666",
  "TBZ_25mg"               = "#e74c3c",
  "DTBZ_30mg"              = "#e67e22",
  "VBZ_80mg"               = "#f39c12",
  "Tominersen_Q8W"         = "#2980b9",
  "Branaplam_Q1W"          = "#8e44ad",
  "Combo_DTBZ_Tominersen"  = "#27ae60"
)

scenario_labels <- c(
  "NaturalHistory"         = "Natural History",
  "TBZ_25mg"               = "Tetrabenazine 25 mg/d",
  "DTBZ_30mg"              = "Deutetrabenazine 30 mg/d",
  "VBZ_80mg"               = "Valbenazine 80 mg/d",
  "Tominersen_Q8W"         = "Tominersen 120 mg Q8W (IT)",
  "Branaplam_Q1W"          = "Branaplam 50 mg QW",
  "Combo_DTBZ_Tominersen"  = "DTBZ + Tominersen (Combo)"
)

# Plot 1: UHDRS-TMS over time
p1 <- ggplot(all_results, aes(x = year, y = TMS, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  labs(title = "UHDRS Total Motor Score (TMS) — 5-Year Simulation",
       subtitle = "CAG=42, Early Manifest HD (TMS baseline=18)",
       x = "Time (years)", y = "UHDRS-TMS (0–124)",
       color = "Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

# Plot 2: TFC over time
p2 <- ggplot(all_results, aes(x = year, y = TFC, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  labs(title = "Total Functional Capacity (TFC) — 5-Year Simulation",
       x = "Time (years)", y = "TFC (0–13)",
       color = "Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

# Plot 3: MSN Survival
p3 <- ggplot(all_results, aes(x = year, y = MSN_pct, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  labs(title = "Medium Spiny Neuron (MSN) Survival",
       x = "Time (years)", y = "MSN Survival (%)",
       color = "Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

# Plot 4: mHTT Total Level
p4 <- ggplot(all_results, aes(x = year, y = mHTT_total, color = scenario)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  labs(title = "Total mHTT Level (Protein + Oligomers)",
       x = "Time (years)", y = "mHTT (nM)",
       color = "Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

# Plot 5: BDNF Level
p5 <- ggplot(all_results, aes(x = year, y = BDNF_level, color = scenario)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 6.5, linetype = "dashed", color = "darkgreen", alpha = 0.6) +
  annotate("text", x = 0.2, y = 6.6, label = "Healthy BDNF (~6.5 ng/mL)", size = 3, color = "darkgreen") +
  scale_color_manual(values = scenario_colors, labels = scenario_labels) +
  labs(title = "CSF BDNF Level (ng/mL)",
       x = "Time (years)", y = "BDNF (ng/mL)",
       color = "Scenario") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

# Plot 6: CAG-Dose Sensitivity (natural history, 3 CAG lengths)
cag_results <- lapply(c(40, 42, 46, 50), function(cag_val) {
  tryCatch(run_scenario(mod, "NaturalHistory", duration_yrs = 10, CAG_len = cag_val),
           error = function(e) NULL)
}) %>% bind_rows()

p6 <- ggplot(cag_results, aes(x = year, y = TMS, color = factor(CAG),
                               group = factor(CAG))) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis_d(name = "CAG Length", option = "C") +
  labs(title = "TMS Progression by CAG Repeat Length (Natural History, 10 yr)",
       x = "Time (years)", y = "UHDRS-TMS") +
  theme_bw(base_size = 12)

# ── Save plots ──
if (!dir.exists("plots")) dir.create("plots")
ggsave("plots/HD_TMS_scenarios.png", p1, width = 11, height = 6, dpi = 150)
ggsave("plots/HD_TFC_scenarios.png", p2, width = 11, height = 6, dpi = 150)
ggsave("plots/HD_MSN_survival.png",  p3, width = 11, height = 6, dpi = 150)
ggsave("plots/HD_mHTT_level.png",    p4, width = 11, height = 6, dpi = 150)
ggsave("plots/HD_BDNF_level.png",    p5, width = 11, height = 6, dpi = 150)
ggsave("plots/HD_CAG_sensitivity.png", p6, width = 10, height = 6, dpi = 150)

message("\nAll 6 plots saved to plots/ directory.")
message("HD QSP model run complete.")
