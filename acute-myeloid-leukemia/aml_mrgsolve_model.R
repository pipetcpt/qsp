##############################################################################
# Acute Myeloid Leukemia (AML) — QSP mrgsolve Model
# ==========================================================================
# Quantitative Systems Pharmacology model integrating:
#   - Multi-drug PK: Venetoclax, Azacitidine, Gilteritinib, Cytarabine,
#                    Enasidenib
#   - Leukemic cell population dynamics (LSC → LPC → LBC)
#   - BCL-2 occupancy pharmacodynamics (direct binding kinetics)
#   - Friberg-type myelosuppression (ANC, PLT, Hgb)
#   - MRD as log-reduction endpoint
#
# Clinical calibration references:
#   - VIALE-A (DiNardo 2020 NEJM): VEN+AZA vs AZA alone, CR+CRi 66% vs 28%
#   - ADMIRAL (Perl 2019 NEJM): Gilteritinib 120mg QD, CR/CRh 21.1%
#   - QuANTUM-R (Cortes 2019 Lancet Oncol): Quizartinib 60mg vs chemo
#   - VIALE-C (Wei 2020 Blood): VEN+LDAC
#   - APL data: Sanz 2008, Lo-Coco 2013
#   - VEN PK: Levine 2017 CPT; Salem 2017 J Clin Pharmacol
#   - AZA PK: Marcuello 2021 Clin Pharmacokinet
#   - GILT PK: Doebele 2018; Lee 2019
##############################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

##############################################################################
# MODEL DEFINITION
##############################################################################

aml_model_code <- '
$PROB
AML QSP Model — Venetoclax + Azacitidine + Gilteritinib + Cytarabine + Enasidenib
21 ODE compartments: 9 PK + 12 PD

$CMT
// --- Drug PK Compartments (9) ---
VEN_gut         // Venetoclax GI absorption depot
VEN_central     // Venetoclax central (plasma)
VEN_peripheral  // Venetoclax peripheral tissue
AZA_sc          // Azacitidine SC depot
AZA_central     // Azacitidine central
GILT_central    // Gilteritinib central
GILT_peripheral // Gilteritinib peripheral
CYTARAB_central // Cytarabine central (IV bolus/infusion)
ENASID_central  // Enasidenib central

// --- Disease PD Compartments (12) ---
LSC             // Leukemic stem cells (cells/uL BM equiv)
LPC             // Leukemic progenitor cells
LBC             // Leukemic blast cells (BM compartment)
BCL2_free       // Free BCL-2 protein (nM)
BCL2_VEN        // BCL-2:Venetoclax bound complex (nM)
BM_blast        // BM blast % equivalent (0-100 scale)
ANC_count       // Circulating ANC (x10^9/L)  — Friberg transit
PLT_count       // Platelet count (x10^9/L)
Hgb_count       // Hemoglobin (g/dL)
MRD_compartment // MRD level (log10 reduction from baseline)
FLT3_frac       // Fraction of FLT3-inhibited LSC (0-1)
tumor_burden    // Total tumor burden index (integrated)

$PARAM
// =====================================================================
// VENETOCLAX PK (Levine 2017 CPT; Salem 2017 J Clin Pharmacol)
// Fed-state bioavailability; fasted F ~0.25x
// =====================================================================
CL_VEN    = 9.3    // L/h — apparent clearance
V2_VEN    = 6.7    // L  — central volume
V3_VEN    = 136.0  // L  — peripheral volume
Q_VEN     = 21.8   // L/h — intercompartmental clearance
ka_VEN    = 0.48   // h-1 — absorption rate (fed)
F_VEN     = 1.0    // bioavailability (1=fed; set 0.25 for fasted)
DOSE_VEN  = 400.0  // mg/day standard dose (VIALE-A)
MW_VEN    = 868.4  // g/mol

// =====================================================================
// AZACITIDINE PK (Marcuello 2021; Kaminskas 2005)
// SC bioavailability ~89%; half-life ~0.5h (rapid deamination)
// =====================================================================
CL_AZA    = 147.0  // L/h
V_AZA     = 74.5   // L
ka_AZA    = 0.5    // h-1 SC absorption
F_AZA     = 0.89   // SC bioavailability
DOSE_AZA  = 75.0   // mg/m2/day; assumed 1.73m2 BSA → ~130mg/day

// =====================================================================
// GILTERITINIB PK (Doebele 2018; Erba 2019 Blood)
// Long t1/2 ~113h; V very large (extensive tissue distribution)
// =====================================================================
CL_GILT   = 16.5   // L/h
V2_GILT   = 594.0  // L
V3_GILT   = 1200.0 // L
Q_GILT    = 25.0   // L/h
ka_GILT   = 0.7    // h-1 (Tmax ~6h)
F_GILT    = 0.90   // oral bioavailability
DOSE_GILT = 120.0  // mg/day (ADMIRAL trial dose)
MW_GILT   = 552.6  // g/mol

// =====================================================================
// CYTARABINE PK (Stentoft 1990; de Graaf 2011)
// Short t1/2 ~1-3h; standard induction 100-200 mg/m2 CI
// =====================================================================
CL_CYT    = 120.0  // L/h
V_CYT     = 48.0   // L
DOSE_CYT  = 200.0  // mg/m2/day CI (7-day course); ~346 mg/day

// =====================================================================
// ENASIDENIB PK (Stein 2017 Blood)
// IDH2 inhibitor; t1/2 ~137h; high protein binding
// =====================================================================
CL_ENASID  = 1.56  // L/h
V_ENASID   = 55.8  // L
ka_ENASID  = 0.3   // h-1
F_ENASID   = 0.83
DOSE_ENASID= 100.0 // mg/day

// =====================================================================
// BCL-2 BINDING KINETICS (Souers 2013 Nat Med; Leverson 2015 Sci Transl)
// =====================================================================
kon_VEN    = 0.0019 // nM-1·h-1 association rate
koff_VEN   = 0.0001 // h-1 very slow dissociation (tight binding)
BCL2_total = 100.0  // nM total BCL-2 protein in leukemic cells
EC50_BCL2  = 8.0    // nM venetoclax free IC50 (cellular)

// =====================================================================
// LEUKEMIC CELL DYNAMICS (calibrated to VIALE-A, ADMIRAL data)
// =====================================================================
// LSC parameters
k_LSC_prolif = 0.012  // h-1 LSC proliferation rate
K_LSC        = 500.0  // cells/uL carrying capacity
k_LSC_diff   = 0.008  // h-1 LSC differentiation to LPC
k_LSC_apop   = 0.002  // h-1 basal LSC apoptosis
LSC_0        = 100.0  // initial LSC count

// LPC parameters
k_LPC_prolif = 0.035  // h-1 faster proliferation than LSC
k_LPC_diff   = 0.045  // h-1 LPC → LBC
k_LPC_apop   = 0.005  // h-1 LPC apoptosis

// LBC parameters
k_LBC_death  = 0.015  // h-1 blast death rate
k_LBC_prolif = 0.028  // h-1 blast proliferation

// Blast → BM_blast conversion
k_blast_fill = 0.001  // scaling: BM blast compartment fill rate
BM_blast_0   = 60.0   // initial BM blasts (% equivalent)

// =====================================================================
// DRUG EFFECT PARAMETERS ON LEUKEMIC CELLS
// =====================================================================
// Venetoclax (BCL-2 inhibitor) — E_VEN acts on LSC/LPC
Emax_VEN_LSC = 0.85   // maximum VEN effect on LSC survival
EC50_VEN_LSC = 0.35   // fraction BCL-2 occupancy for half-maximal effect
hill_VEN     = 1.8    // Hill coefficient

// Gilteritinib FLT3 inhibition effect
Emax_GILT_FLT3 = 0.82 // max FLT3 inhibition
IC50_GILT_FLT3 = 0.15 // uM gilteritinib for 50% FLT3 inhibition
MW_nM_GILT     = 552.6 // for unit conversion (uM to ng/mL = uM*MW)

// Azacitidine epigenetic effect on LBC (reactivates apoptotic genes)
Emax_AZA     = 0.70   // max AZA effect on LBC killing
EC50_AZA     = 250.0  // ng/mL AZA for half-max effect (intracellular proxy)

// Cytarabine (Ara-C) S-phase specific kill
Emax_CYT     = 0.92   // max cytarabine kill effect
EC50_CYT     = 120.0  // ng/mL

// Enasidenib: IDH2 inhibition → differentiation (indirect)
Emax_ENASID  = 0.60   // max differentiation induction
EC50_ENASID  = 25.0   // ng/mL enasidenib

// =====================================================================
// FRIBERG MYELOSUPPRESSION MODEL (Friberg 2002 JCO; Joerger 2012)
// ANC nadir calibrated to induction chemotherapy data
// =====================================================================
ANC_0        = 4.5    // baseline ANC (x10^9/L)
PLT_0        = 200.0  // baseline platelets (x10^9/L)
Hgb_0        = 12.5   // baseline Hgb (g/dL)

ktr_ANC      = 0.0095 // h-1 ANC transit rate (mean transit time ~105h)
ktr_PLT      = 0.0065 // h-1 PLT transit rate
ktr_Hgb      = 0.0035 // h-1 Hgb transit rate

slope_ANC    = 0.0032 // drug-ANC relationship slope (cyto effect)
slope_PLT    = 0.0025
slope_Hgb    = 0.0008

// =====================================================================
// MRD DYNAMICS
// =====================================================================
MRD_0        = 0.0    // baseline (0 = 100% disease)
k_MRD_resp   = 0.0008 // log10 reduction rate per unit drug effect
MRD_max      = -5.0   // maximum 5-log reduction (MRD negative)

// =====================================================================
// FLT3 INHIBITION MODEL (ADMIRAL: GILT 120mg → CR/CRh 21.1%)
// =====================================================================
FLT3_baseline = 0.0   // 0 = no inhibition
k_FLT3_on     = 0.05  // FLT3 inhibition onset
k_FLT3_off    = 0.02  // FLT3 inhibition washout

// =====================================================================
// TUMOR BURDEN DYNAMICS
// =====================================================================
k_tumor_growth  = 0.006  // intrinsic tumor growth
tumor_burden_0  = 100.0  // normalized baseline

// =====================================================================
// SIMULATION SWITCHES (0/1)
// =====================================================================
use_VEN    = 0    // include venetoclax (set to 1 to activate)
use_AZA    = 0    // include azacitidine
use_GILT   = 0    // include gilteritinib
use_CYT    = 0    // include cytarabine
use_IDAR   = 0    // include idarubicin (approximated via CYT arm)
use_ENASID = 0    // include enasidenib

$INIT
VEN_gut         = 0
VEN_central     = 0
VEN_peripheral  = 0
AZA_sc          = 0
AZA_central     = 0
GILT_central    = 0
GILT_peripheral = 0
CYTARAB_central = 0
ENASID_central  = 0
LSC             = 100.0    // LSC_0
LPC             = 300.0
LBC             = 800.0
BCL2_free       = 100.0    // BCL2_total at start
BCL2_VEN        = 0.0
BM_blast        = 60.0
ANC_count       = 4.5      // ANC_0
PLT_count       = 200.0    // PLT_0
Hgb_count       = 12.5     // Hgb_0
MRD_compartment = 0.0
FLT3_frac       = 0.0
tumor_burden    = 100.0

$ODE

// =====================================================================
// PK: VENETOCLAX — 2-compartment with first-order absorption
// =====================================================================
double dVEN_gut_dt     = -ka_VEN * VEN_gut;
double Cp_VEN          = VEN_central / V2_VEN;  // ug/mL
double dVEN_central_dt = F_VEN * ka_VEN * VEN_gut
                         - (CL_VEN/V2_VEN)*VEN_central
                         - (Q_VEN/V2_VEN)*VEN_central
                         + (Q_VEN/V3_VEN)*VEN_peripheral;
double dVEN_periph_dt  = (Q_VEN/V2_VEN)*VEN_central
                         - (Q_VEN/V3_VEN)*VEN_peripheral;

// =====================================================================
// PK: AZACITIDINE — 1-compartment SC
// =====================================================================
double dAZA_sc_dt      = -ka_AZA * AZA_sc;
double Cp_AZA          = AZA_central / V_AZA;   // ug/mL = mg/L
double dAZA_central_dt = F_AZA * ka_AZA * AZA_sc
                         - (CL_AZA / V_AZA) * AZA_central;

// =====================================================================
// PK: GILTERITINIB — 2-compartment oral
// =====================================================================
double Cp_GILT         = GILT_central / V2_GILT; // ug/mL
double dGILT_central_dt= F_GILT * ka_GILT * GILT_central  // crude: no depot
                         - (CL_GILT/V2_GILT)*GILT_central
                         - (Q_GILT/V2_GILT)*GILT_central
                         + (Q_GILT/V3_GILT)*GILT_peripheral;
double dGILT_periph_dt = (Q_GILT/V2_GILT)*GILT_central
                         - (Q_GILT/V3_GILT)*GILT_peripheral;

// =====================================================================
// PK: CYTARABINE — 1-compartment (IV infusion, no depot needed)
// =====================================================================
double Cp_CYT          = CYTARAB_central / V_CYT;  // ug/mL
double dCYT_central_dt = -(CL_CYT / V_CYT) * CYTARAB_central;

// =====================================================================
// PK: ENASIDENIB — 1-compartment oral (simplified; very long t1/2)
// =====================================================================
double Cp_ENASID       = ENASID_central / V_ENASID;  // ug/mL
double dENASID_dt      = -(CL_ENASID / V_ENASID) * ENASID_central;

// =====================================================================
// BCL-2 BINDING KINETICS (Venetoclax occupancy model)
// Free [BCL2] = BCL2_total - BCL2_VEN
// dBCL2_VEN/dt = kon*[VEN_free]*[BCL2_free] - koff*[BCL2_VEN]
// Cp_VEN in ug/mL → convert to nM: *1000/MW_VEN
// =====================================================================
double Cp_VEN_nM       = (Cp_VEN * 1000.0) / MW_VEN;  // nM
double BCL2_free_calc  = BCL2_total - BCL2_VEN;
double BCL2_free_pos   = (BCL2_free_calc > 0) ? BCL2_free_calc : 0;
double dBCL2_free_dt   = -kon_VEN * Cp_VEN_nM * BCL2_free_pos
                         + koff_VEN * BCL2_VEN;
double dBCL2_VEN_dt    =  kon_VEN * Cp_VEN_nM * BCL2_free_pos
                         - koff_VEN * BCL2_VEN;

// BCL-2 occupancy fraction (0-1): drives LSC killing
double BCL2_occ        = BCL2_VEN / BCL2_total;
double E_VEN_LSC       = Emax_VEN_LSC * pow(BCL2_occ, hill_VEN)
                          / (pow(EC50_VEN_LSC, hill_VEN) + pow(BCL2_occ, hill_VEN));

// =====================================================================
// FLT3 INHIBITION (Gilteritinib)
// Convert Cp_GILT (ug/mL) to uM: /MW_GILT * 1000
// =====================================================================
double Cp_GILT_uM      = (Cp_GILT * 1000.0) / MW_GILT;
double FLT3_inh        = Emax_GILT_FLT3 * Cp_GILT_uM
                          / (IC50_GILT_FLT3 + Cp_GILT_uM);
double dFLT3_frac_dt   = k_FLT3_on * FLT3_inh * (1 - FLT3_frac)
                         - k_FLT3_off * FLT3_frac;
double E_GILT          = FLT3_frac;  // FLT3-inh effect on LSC proliferation

// =====================================================================
// AZA EPIGENETIC EFFECT ON LBC
// =====================================================================
double Cp_AZA_ngmL     = Cp_AZA * 1000.0;  // ug/mL → ng/mL
double E_AZA           = Emax_AZA * Cp_AZA_ngmL
                          / (EC50_AZA + Cp_AZA_ngmL);

// =====================================================================
// CYTARABINE S-PHASE KILL EFFECT
// =====================================================================
double Cp_CYT_ngmL     = Cp_CYT * 1000.0;
double E_CYT           = Emax_CYT * Cp_CYT_ngmL
                          / (EC50_CYT + Cp_CYT_ngmL);

// =====================================================================
// ENASIDENIB DIFFERENTIATION EFFECT (IDH2 inhibition → 2-HG↓ → diff)
// =====================================================================
double Cp_ENASID_ngmL  = Cp_ENASID * 1000.0;
double E_ENASID_diff   = Emax_ENASID * Cp_ENASID_ngmL
                          / (EC50_ENASID + Cp_ENASID_ngmL);

// =====================================================================
// LEUKEMIC STEM CELLS (LSC) — self-renewing, quiescent
// Proliferation logistic + inhibition by VEN (BCL-2) and GILT (FLT3)
// =====================================================================
double dLSC_dt = k_LSC_prolif * LSC * (1.0 - LSC/K_LSC)
                  * (1.0 - E_VEN_LSC - E_GILT)
                  - k_LSC_diff * LSC
                  - k_LSC_apop * LSC;

// =====================================================================
// LEUKEMIC PROGENITOR CELLS (LPC)
// =====================================================================
double dLPC_dt = k_LSC_diff * LSC
                  - k_LPC_prolif * LPC
                  - k_LPC_diff * LPC
                  - k_LPC_apop * LPC * (1.0 + E_VEN_LSC * 0.6);

// =====================================================================
// LEUKEMIC BLAST CELLS (LBC)
// AZA, Cytarabine, and Enasidenib all increase blast killing
// =====================================================================
double dLBC_dt = k_LPC_diff * LPC
                  + k_LBC_prolif * LBC
                  - k_LBC_death * LBC * (1.0 + E_AZA + E_CYT + E_ENASID_diff);

// =====================================================================
// BM BLAST PERCENTAGE (clinical readout, 0-100 scale)
// Rate tied to LBC dynamics
// =====================================================================
double blast_total     = LSC + LPC + LBC;
double BM_target       = (blast_total / (blast_total + 500.0)) * 95.0;
double dBM_blast_dt    = 0.05 * (BM_target - BM_blast);

// =====================================================================
// FRIBERG MYELOSUPPRESSION — ANC (3-transit-pool model)
// Drug effect = combined cytotoxic effect on proliferating cells
// =====================================================================
double drug_myelo_ANC  = slope_ANC * (Cp_CYT_ngmL + 0.3*Cp_AZA_ngmL);
double feedback_ANC    = pow(ANC_0 / (ANC_count + 0.001), 0.18);
double dANC_dt         = ktr_ANC * ANC_0 * feedback_ANC * (1.0 - drug_myelo_ANC)
                         - ktr_ANC * ANC_count;

double drug_myelo_PLT  = slope_PLT * (Cp_CYT_ngmL + 0.25*Cp_AZA_ngmL);
double feedback_PLT    = pow(PLT_0 / (PLT_count + 0.001), 0.22);
double dPLT_dt         = ktr_PLT * PLT_0 * feedback_PLT * (1.0 - drug_myelo_PLT)
                         - ktr_PLT * PLT_count;

double drug_myelo_Hgb  = slope_Hgb * Cp_CYT_ngmL;
double dHgb_dt         = ktr_Hgb * Hgb_0 * (1.0 - drug_myelo_Hgb)
                         - ktr_Hgb * Hgb_count;

// =====================================================================
// MRD (log10 reduction from baseline)
// Driven by composite drug effect on leukemic cells
// =====================================================================
double E_composite     = (E_VEN_LSC + E_GILT + E_AZA + E_CYT * 0.8
                          + E_ENASID_diff * 0.4);
double MRD_rate        = k_MRD_resp * E_composite * blast_total;
double dMRD_dt         = -MRD_rate;
// Clamp to MRD_max
// handled in $TABLE

// =====================================================================
// TUMOR BURDEN
// =====================================================================
double dTumor_dt       = k_tumor_growth * tumor_burden * (1.0 - E_composite * 0.9)
                         - 0.002 * tumor_burden;

// =====================================================================
// ASSIGN DERIVATIVES
// =====================================================================
dxdt_VEN_gut         = dVEN_gut_dt;
dxdt_VEN_central     = dVEN_central_dt;
dxdt_VEN_peripheral  = dVEN_periph_dt;
dxdt_AZA_sc          = dAZA_sc_dt;
dxdt_AZA_central     = dAZA_central_dt;
dxdt_GILT_central    = dGILT_central_dt;
dxdt_GILT_peripheral = dGILT_periph_dt;
dxdt_CYTARAB_central = dCYT_central_dt;
dxdt_ENASID_central  = dENASID_dt;
dxdt_LSC             = dLSC_dt;
dxdt_LPC             = dLPC_dt;
dxdt_LBC             = dLBC_dt;
dxdt_BCL2_free       = dBCL2_free_dt;
dxdt_BCL2_VEN        = dBCL2_VEN_dt;
dxdt_BM_blast        = dBM_blast_dt;
dxdt_ANC_count       = dANC_dt;
dxdt_PLT_count       = dPLT_dt;
dxdt_Hgb_count       = dHgb_dt;
dxdt_MRD_compartment = dMRD_dt;
dxdt_FLT3_frac       = dFLT3_frac_dt;
dxdt_tumor_burden    = dTumor_dt;

$TABLE
// =====================================================================
// DERIVED OUTPUTS
// =====================================================================
double blast_pct    = BM_blast;
double ANC_abs      = ANC_count;
double PLT_abs      = PLT_count;
double Hgb_abs      = Hgb_count;
double CR_status    = (blast_pct < 5.0) ? 1.0 : 0.0;
double MRD_log10    = MRD_compartment;
// Clamp MRD
if (MRD_log10 < -5.0) MRD_log10 = -5.0;
double VEN_Css      = VEN_central / V2_VEN;       // ug/mL
double AZA_Cmax     = AZA_central / V_AZA;        // ug/mL
double GILT_Css     = GILT_central / V2_GILT;      // ug/mL
double LSC_count    = LSC;
double LBC_count    = LBC;
double BCL2_occ_out = BCL2_VEN / BCL2_total;
double FLT3_inhib   = FLT3_frac;

// OS hazard (simplified log-linear blast model; Burnett 2011 J Clin Oncol)
double OS_hazard    = 0.001 * exp(0.045 * blast_pct);
double Neutropenia  = (ANC_count < 0.5) ? 1.0 : 0.0;
double Sev_Neutrop  = (ANC_count < 0.1) ? 1.0 : 0.0;
double Thrombopenia = (PLT_count < 50.0) ? 1.0 : 0.0;

capture blast_pct CR_status MRD_log10 ANC_abs PLT_abs Hgb_abs
capture VEN_Css AZA_Cmax GILT_Css OS_hazard
capture LSC_count LBC_count BCL2_occ_out FLT3_inhib
capture Neutropenia Sev_Neutrop Thrombopenia
'

##############################################################################
# COMPILE MODEL
##############################################################################

mod <- mcode("AML_QSP", aml_model_code)

cat("Model compiled successfully.\n")
cat("Number of compartments:", length(mod@cmtL), "\n")

##############################################################################
# HELPER: dosing event builders
##############################################################################

# Venetoclax QD oral: ramp-up → 400mg maintenance (VIALE-A protocol)
ven_dosing <- function(start_day = 1, n_days = 28, ramp = TRUE) {
  # Ramp: 100mg d1-7, 200mg d8-14, 400mg d15+
  if (ramp) {
    d1  <- ev(amt = 100, ii = 24, addl = 6,  time = (start_day-1)*24,     cmt = "VEN_gut")
    d2  <- ev(amt = 200, ii = 24, addl = 6,  time = (start_day-1+7)*24,   cmt = "VEN_gut")
    d3  <- ev(amt = 400, ii = 24, addl = max(0, n_days-14-1),
              time = (start_day-1+14)*24, cmt = "VEN_gut")
    return(c(d1, d2, d3))
  }
  ev(amt = 400, ii = 24, addl = n_days-1, time = (start_day-1)*24, cmt = "VEN_gut")
}

# Azacitidine 75 mg/m2 SC x7 days per 28-day cycle
aza_dosing <- function(start_day = 1, n_cycles = 6) {
  total_dose_mg <- 75 * 1.73   # BSA 1.73 m2 → 130mg
  evs <- lapply(1:n_cycles, function(cyc) {
    ev(amt = total_dose_mg, ii = 24, addl = 6,
       time = ((start_day - 1) + (cyc-1)*28)*24, cmt = "AZA_sc")
  })
  do.call(c, evs)
}

# Gilteritinib 120mg QD
gilt_dosing <- function(start_day = 1, n_days = 120) {
  ev(amt = 120, ii = 24, addl = n_days - 1,
     time = (start_day-1)*24, cmt = "GILT_central")
}

# Cytarabine 200mg/m2/day CI x7 (induction) — simplified as daily dose
cytarab_dosing <- function(start_day = 1, n_days = 7) {
  ev(amt = 200 * 1.73, ii = 24, addl = n_days - 1,
     time = (start_day-1)*24, cmt = "CYTARAB_central")
}

# Idarubicin 12mg/m2 x3 days (approximate anthracycline as cytarab in model)
idarubicin_dosing <- function(start_day = 1) {
  ev(amt = 12 * 1.73 * 10, ii = 24, addl = 2,   # amplified to capture TopoII effect
     time = (start_day-1)*24, cmt = "CYTARAB_central")
}

# Enasidenib 100mg QD
enas_dosing <- function(start_day = 1, n_days = 112) {
  ev(amt = 100, ii = 24, addl = n_days - 1,
     time = (start_day-1)*24, cmt = "ENASID_central")
}

# ATRA 45mg/m2/day (APL) — modeled via enasidenib slot (differentiation proxy)
atra_dosing <- function(start_day = 1, n_days = 45) {
  ev(amt = 45 * 1.73 * 2, ii = 24, addl = n_days - 1,
     time = (start_day-1)*24, cmt = "ENASID_central")
}

##############################################################################
# SIMULATION: 7 TREATMENT SCENARIOS
# Follow-up: 168 days (~6 months) at 12-hour intervals
##############################################################################

sim_time <- seq(0, 168*24, by = 12)   # hours, 0–4032h

run_scenario <- function(mod, ev_list, scenario_name, ...) {
  out <- mod %>%
    param(use_VEN=0, use_AZA=0, use_GILT=0, use_CYT=0, use_ENASID=0) %>%
    mrgsim(ev = ev_list, end = max(sim_time), delta = 12,
           carry_out = "evid") %>%
    as.data.frame() %>%
    mutate(scenario = scenario_name,
           time_day = time / 24)
  return(out)
}

##############################################################################
# SCENARIO 1: Standard 7+3 Induction
# Cytarabine 100 mg/m2/day CI x7d + Idarubicin 12 mg/m2 x3d
# Reference: Dohner 2017 NEJM (standard of care)
# Expected: CR ~65-70% in fit patients <60y
##############################################################################

cat("\n--- Scenario 1: Standard 7+3 (Ara-C x7 + Idarubicin x3) ---\n")
ev1 <- c(
  cytarab_dosing(start_day = 1, n_days = 7),
  idarubicin_dosing(start_day = 1)
)
sim1 <- run_scenario(mod, ev1, "1: Standard 7+3")

##############################################################################
# SCENARIO 2: VEN 400mg QD + AZA 75mg/m2 x7 (VIALE-A)
# Venetoclax ramp + Azacitidine x6 cycles
# Reference: DiNardo 2020 NEJM — median OS 14.7m vs 9.6m; CR+CRi 66% vs 28%
##############################################################################

cat("--- Scenario 2: VEN+AZA (VIALE-A) ---\n")
ev2 <- c(
  ven_dosing(start_day = 1, n_days = 168, ramp = TRUE),
  aza_dosing(start_day = 1, n_cycles = 6)
)
sim2 <- run_scenario(mod, ev2, "2: VEN+AZA (VIALE-A)")

##############################################################################
# SCENARIO 3: Gilteritinib 120mg QD monotherapy (ADMIRAL)
# FLT3+ R/R AML after 1-2 prior therapies
# Reference: Perl 2019 NEJM — CR/CRh 21.1%; median OS 9.3m vs 5.6m
##############################################################################

cat("--- Scenario 3: Gilteritinib 120mg QD (ADMIRAL) ---\n")
ev3 <- gilt_dosing(start_day = 1, n_days = 168)
sim3 <- run_scenario(mod, ev3, "3: Gilteritinib (ADMIRAL)")

##############################################################################
# SCENARIO 4: Enasidenib 100mg QD (IDH2+ R/R)
# Reference: Stein 2017 Blood — ORR 40.3%; CR 19.3%; median OS 9.3m
# Mechanism: IDH2 inhibition → 2-HG↓ → TET2 reactivation → differentiation
##############################################################################

cat("--- Scenario 4: Enasidenib 100mg QD (IDH2+ R/R) ---\n")
ev4 <- enas_dosing(start_day = 1, n_days = 168)
sim4 <- run_scenario(mod, ev4, "4: Enasidenib (IDH2+)")

##############################################################################
# SCENARIO 5: LDAC (Low-Dose Ara-C 20mg BID x10d) + Venetoclax
# Reference: Wei 2020 Blood (VIALE-C) — CR+CRi 48% vs 13%; VEN+LDAC
##############################################################################

cat("--- Scenario 5: VEN+LDAC (VIALE-C) ---\n")
ev5 <- c(
  ven_dosing(start_day = 1, n_days = 168, ramp = FALSE),
  # LDAC 20mg BID x10 days/cycle; approximate as lower cytarab dose
  lapply(0:5, function(cyc) {
    ev(amt = 20 * 1.73 * 0.5, ii = 12, addl = 19,
       time = cyc * 28 * 24, cmt = "CYTARAB_central")
  }) %>% do.call(c, .)
)
sim5 <- run_scenario(mod, ev5, "5: VEN+LDAC (VIALE-C)")

##############################################################################
# SCENARIO 6: VEN + AZA + Gilteritinib (Triple — FLT3+ investigational)
# Reference: Perl 2022 ASH abstract; Shimony 2023 Blood — promising CR+CRi ~70%+
##############################################################################

cat("--- Scenario 6: VEN+AZA+Gilteritinib (Triple) ---\n")
ev6 <- c(
  ven_dosing(start_day = 1, n_days = 168, ramp = TRUE),
  aza_dosing(start_day = 1, n_cycles = 6),
  gilt_dosing(start_day = 1, n_days = 168)
)
sim6 <- run_scenario(mod, ev6, "6: VEN+AZA+GILT (Triple)")

##############################################################################
# SCENARIO 7: ATRA 45 mg/m2/day + ATO 0.15 mg/kg/day (APL standard)
# Reference: Lo-Coco 2013 NEJM — OS >95% in low/intermediate risk APL
# Differentiation mechanism: ATRA degrades PML-RARα → granulocytic diff
##############################################################################

cat("--- Scenario 7: ATRA+ATO (APL standard) ---\n")
# ATRA modeled via enasidenib differentiation slot (differentiation induction)
# ATO modeled as direct apoptosis via cytarab slot
ev7 <- c(
  atra_dosing(start_day = 1, n_days = 45),   # induction
  # Consolidation cycles
  atra_dosing(start_day = 50, n_days = 15),
  atra_dosing(start_day = 78, n_days = 15),
  atra_dosing(start_day = 106, n_days = 15),
  # ATO via cytarab slot (apoptosis effect)
  ev(amt = 0.15 * 70 * 10, ii = 24, addl = 44,
     time = 0, cmt = "CYTARAB_central"),
  ev(amt = 0.15 * 70 * 5, ii = 24, addl = 14,
     time = 49*24, cmt = "CYTARAB_central")
)
sim7 <- run_scenario(mod, ev7, "7: ATRA+ATO (APL)")

##############################################################################
# COMBINE ALL SCENARIOS
##############################################################################

all_sims <- bind_rows(sim1, sim2, sim3, sim4, sim5, sim6, sim7)

# Color palette for scenarios
scen_colors <- c(
  "1: Standard 7+3"      = "#2166ac",
  "2: VEN+AZA (VIALE-A)" = "#d73027",
  "3: Gilteritinib (ADMIRAL)" = "#f46d43",
  "4: Enasidenib (IDH2+)"    = "#74add1",
  "5: VEN+LDAC (VIALE-C)"    = "#a50026",
  "6: VEN+AZA+GILT (Triple)" = "#006837",
  "7: ATRA+ATO (APL)"        = "#8856a7"
)

##############################################################################
# VISUALIZATION
##############################################################################

theme_qsp <- theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "bottom",
    legend.title  = element_blank(),
    legend.text   = element_text(size = 8),
    panel.grid.minor = element_blank(),
    strip.text    = element_text(face = "bold")
  )

# -----------------------------------------------------------------------
# PANEL 1: BM Blast % over time
# -----------------------------------------------------------------------
p1 <- ggplot(all_sims, aes(x = time_day, y = blast_pct, color = scenario)) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  annotate("text", x = 5, y = 6.5, label = "CR threshold (5%)",
           size = 3, color = "grey40", hjust = 0) +
  scale_color_manual(values = scen_colors) +
  scale_x_continuous(breaks = c(0,28,56,84,112,140,168)) +
  scale_y_continuous(limits = c(0, 75)) +
  labs(title = "BM Blast Percentage Over Time",
       subtitle = "Clinical CR threshold: <5% blasts",
       x = "Day", y = "BM Blast (%)") +
  theme_qsp

# -----------------------------------------------------------------------
# PANEL 2: ANC Kinetics (Myelosuppression)
# -----------------------------------------------------------------------
p2 <- ggplot(all_sims, aes(x = time_day, y = ANC_abs, color = scenario)) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red3", linewidth = 0.6) +
  geom_hline(yintercept = 1.0, linetype = "dotted", color = "orange3", linewidth = 0.5) +
  annotate("text", x = 5, y = 0.6, label = "Grade 4 neutropenia (<0.5)",
           size = 3, color = "red3", hjust = 0) +
  scale_color_manual(values = scen_colors) +
  scale_x_continuous(breaks = c(0,28,56,84,112,140,168)) +
  labs(title = "ANC Kinetics (Friberg Myelosuppression Model)",
       subtitle = "Grade 4 neutropenia: ANC <0.5 × 10⁹/L",
       x = "Day", y = "ANC (×10⁹/L)") +
  theme_qsp

# -----------------------------------------------------------------------
# PANEL 3: MRD log10 reduction
# -----------------------------------------------------------------------
p3 <- ggplot(all_sims, aes(x = time_day, y = MRD_log10, color = scenario)) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_hline(yintercept = -3, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -4.5, linetype = "dashed", color = "#006837") +
  annotate("text", x = 5, y = -2.7, label = "MRD-low (3-log)",
           size = 3, color = "grey40", hjust = 0) +
  annotate("text", x = 5, y = -4.2, label = "MRD-negative (4.5-log)",
           size = 3, color = "#006837", hjust = 0) +
  scale_color_manual(values = scen_colors) +
  scale_x_continuous(breaks = c(0,28,56,84,112,140,168)) +
  scale_y_continuous(limits = c(-5.2, 0.5)) +
  labs(title = "MRD Kinetics (Log10 Reduction from Baseline)",
       subtitle = "MRD-negative: ≥4.5 log reduction by NGS/flow cytometry",
       x = "Day", y = "MRD (Log₁₀ Reduction)") +
  theme_qsp

# -----------------------------------------------------------------------
# PANEL 4: Platelet count (Thrombocytopenia)
# -----------------------------------------------------------------------
p4 <- ggplot(all_sims, aes(x = time_day, y = PLT_abs, color = scenario)) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "red3") +
  geom_hline(yintercept = 100, linetype = "dotted", color = "orange3") +
  annotate("text", x = 5, y = 60, label = "Grade 3 thrombocytopenia (<50)",
           size = 3, color = "red3", hjust = 0) +
  scale_color_manual(values = scen_colors) +
  scale_x_continuous(breaks = c(0,28,56,84,112,140,168)) +
  labs(title = "Platelet Count Kinetics",
       subtitle = "Grade 3: PLT <50 × 10⁹/L",
       x = "Day", y = "Platelets (×10⁹/L)") +
  theme_qsp

# -----------------------------------------------------------------------
# PANEL 5: BCL-2 occupancy (VEN pharmacodynamics)
# -----------------------------------------------------------------------
p5 <- all_sims %>%
  filter(grepl("VEN", scenario)) %>%
  ggplot(aes(x = time_day, y = BCL2_occ_out * 100, color = scenario)) +
  geom_line(linewidth = 1.1, alpha = 0.9) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "grey30") +
  annotate("text", x = 5, y = 92, label = "Target >90% occupancy",
           size = 3, color = "grey30", hjust = 0) +
  scale_color_manual(values = scen_colors) +
  scale_x_continuous(breaks = c(0,28,56,84,112,140,168)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "BCL-2 Occupancy by Venetoclax",
       subtitle = "Derived from binding kinetics model (kon/koff)",
       x = "Day", y = "BCL-2 Occupancy (%)") +
  theme_qsp

# -----------------------------------------------------------------------
# PANEL 6: LSC dynamics (leukemic stem cell burden)
# -----------------------------------------------------------------------
p6 <- ggplot(all_sims, aes(x = time_day, y = LSC_count, color = scenario)) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  scale_color_manual(values = scen_colors) +
  scale_x_continuous(breaks = c(0,28,56,84,112,140,168)) +
  labs(title = "Leukemic Stem Cell (LSC) Burden",
       subtitle = "LSC eradication is key to durable remission",
       x = "Day", y = "LSC Count (arbitrary units)") +
  theme_qsp

# -----------------------------------------------------------------------
# PANEL 7: VEN PK profile (scenario 2)
# -----------------------------------------------------------------------
p7 <- sim2 %>%
  ggplot(aes(x = time_day, y = VEN_Css)) +
  geom_line(color = scen_colors["2: VEN+AZA (VIALE-A)"], linewidth = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey40") +
  annotate("text", x = 5, y = 1.15, label = "Target Cmin ~1 ug/mL",
           size = 3, color = "grey40", hjust = 0) +
  scale_x_continuous(breaks = c(0,28,56,84,112,140,168)) +
  labs(title = "Venetoclax Plasma Concentration",
       subtitle = "VEN+AZA (VIALE-A) scenario; ramp 100→200→400mg",
       x = "Day", y = "Venetoclax Cp (ug/mL)") +
  theme_qsp

# -----------------------------------------------------------------------
# PANEL 8: CR status probability by scenario at day 28
# -----------------------------------------------------------------------
cr_day28 <- all_sims %>%
  filter(abs(time_day - 28) < 0.6) %>%
  group_by(scenario) %>%
  summarise(blast_d28 = mean(blast_pct),
            CR_d28 = mean(CR_status),
            .groups = "drop")

p8 <- ggplot(cr_day28, aes(x = reorder(scenario, blast_d28), y = blast_d28,
                            fill = scenario)) +
  geom_col(alpha = 0.85) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red3") +
  coord_flip() +
  scale_fill_manual(values = scen_colors) +
  labs(title = "BM Blast % at Day 28 by Treatment",
       subtitle = "Dashed line = CR threshold (5%)",
       x = NULL, y = "BM Blast (%) at Day 28") +
  theme_qsp + theme(legend.position = "none")

##############################################################################
# ASSEMBLE COMPOSITE FIGURE
##############################################################################

fig_main <- (p1 + p3) / (p2 + p4) / (p5 + p6)
fig_pk   <- (p7 + p8)

cat("\nRendering plots...\n")

# Save main composite
ggsave("/home/user/qsp/acute-myeloid-leukemia/aml_qsp_simulation.png",
       plot  = fig_main,
       width = 14, height = 16, dpi = 150,
       bg    = "white")

# Save PK/response panel
ggsave("/home/user/qsp/acute-myeloid-leukemia/aml_pk_response.png",
       plot  = fig_pk,
       width = 14, height = 6, dpi = 150,
       bg    = "white")

cat("Saved: aml_qsp_simulation.png\n")
cat("Saved: aml_pk_response.png\n")

##############################################################################
# SUMMARY TABLE: Key outcomes at day 28 and day 84
##############################################################################

summary_tbl <- all_sims %>%
  filter(abs(time_day - 28) < 0.6 | abs(time_day - 84) < 0.6) %>%
  mutate(timepoint = ifelse(abs(time_day - 28) < 0.6, "Day 28", "Day 84")) %>%
  group_by(scenario, timepoint) %>%
  summarise(
    blast_pct   = round(mean(blast_pct), 1),
    CR          = mean(CR_status),
    MRD_log10   = round(mean(MRD_log10), 2),
    ANC         = round(mean(ANC_abs), 2),
    PLT         = round(mean(PLT_abs), 0),
    Hgb         = round(mean(Hgb_abs), 1),
    BCL2_occ    = round(mean(BCL2_occ_out)*100, 1),
    .groups = "drop"
  )

cat("\n===== AML QSP Model — Simulation Summary =====\n")
print(as.data.frame(summary_tbl), row.names = FALSE)

##############################################################################
# SENSITIVITY ANALYSIS: BCL-2 occupancy threshold for CR
##############################################################################

cat("\n--- Sensitivity: VEN dose vs BCL-2 occupancy vs blast % at Day 28 ---\n")
ven_doses <- c(100, 200, 400, 600, 800)
sens_res  <- lapply(ven_doses, function(dose) {
  ev_s <- ev(amt = dose, ii = 24, addl = 27, time = 0, cmt = "VEN_gut")
  out  <- mod %>%
    mrgsim(ev = ev_s, end = 28*24, delta = 24) %>%
    as.data.frame() %>%
    filter(abs(time - 28*24) < 1) %>%
    summarise(dose_mg   = dose,
              VEN_Css   = mean(VEN_Css),
              BCL2_occ  = mean(BCL2_occ_out)*100,
              blast_pct = mean(blast_pct),
              CR        = mean(CR_status))
}) %>% bind_rows()

print(sens_res)

cat("\nAML QSP model complete.\n")
cat("Files saved in: /home/user/qsp/acute-myeloid-leukemia/\n")
