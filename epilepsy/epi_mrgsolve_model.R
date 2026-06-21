################################################################################
## Epilepsy (뇌전증) QSP Model — mrgsolve ODE Implementation
## Version 1.0 | 2026-06-21
##
## Model structure:
##   16 ODE compartments:
##     VPA PK: gut (AGUT) · central (ACENT) · peripheral (APER)
##     LEV PK: gut (BGUT) · central (BCENT)
##     CBZ PK: gut (CGUT) · central (CCENT) · epoxide metabolite (CMETA)
##     LTG PK: gut (DGUT) · central (DCENT)
##     PD:     GABA brain level (GABA) · synaptic glutamate (SYNAP)
##             SV2A occupancy (SV2A_OCC) · Na-channel blockade (NAV_BLOCK)
##             Seizure threshold (STHRES) · P-gp expression (PGP)
##
## Treatment scenarios:
##   1. Untreated — baseline seizure frequency
##   2. VPA monotherapy (1,000 mg/day)
##   3. LEV monotherapy (3,000 mg/day)
##   4. CBZ monotherapy (600 mg/day)
##   5. LTG monotherapy (200 mg/day)
##   6. VPA + LTG polytherapy (DDI: VPA inhibits LTG glucuronidation)
##   7. CBZ + LTG polytherapy (DDI: CBZ induces LTG clearance)
##   8. Drug-resistant epilepsy (P-gp overexpression x3)
##   9. Status epilepticus rescue — IV lorazepam + phenobarbital
##  10. TSC / mTOR pathway (everolimus adjunct simulation)
##
## Key references:
##   Brodie MJ et al. Lancet Neurol 2020; Perucca E et al. Lancet 2019;
##   Loscher W et al. Nat Rev Drug Discov 2019; Deckers LP et al. Epilepsia 2000
################################################################################

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

# ─── Model Definition ────────────────────────────────────────────────────────

epi_code <- '
$PROB
Epilepsy QSP Model
Ion Channel Biology · AED PK/PD · Seizure Threshold · Drug Resistance
Compounds modelled: VPA, LEV, CBZ (+CBZ-E), LTG
Authors: QSP Routine | 2026-06-21

$PARAM @annotated
// === VPA Pharmacokinetics ===
ka_VPA    : 1.20  : VPA absorption rate constant (1/h)
Vc_VPA    : 9.10  : VPA central volume of distribution (L)
Vp_VPA    : 12.0  : VPA peripheral volume (L)
CL_VPA    : 0.47  : VPA total clearance (L/h) [t1/2~14h]
Q_VPA     : 0.23  : VPA inter-compartment clearance (L/h)
F_VPA     : 0.90  : VPA oral bioavailability fraction
fu_VPA    : 0.10  : VPA unbound fraction (90% protein bound)

// === LEV Pharmacokinetics ===
ka_LEV    : 1.50  : LEV absorption rate constant (1/h)
Vc_LEV    : 42.0  : LEV central volume (L) [0.6 L/kg, 70 kg]
CL_LEV    : 3.80  : LEV renal clearance (L/h) [t1/2~7h]
F_LEV     : 1.00  : LEV oral bioavailability fraction

// === CBZ Pharmacokinetics (with autoinduction) ===
ka_CBZ    : 0.50  : CBZ absorption rate constant (1/h)
Vc_CBZ    : 51.0  : CBZ central volume (L)
CL_CBZ0   : 3.00  : CBZ initial clearance (L/h, before autoinduction)
CL_CBZ_max: 6.50  : CBZ maximal clearance (L/h, fully autoinduced)
EC50_auto : 3.00  : CBZ concentration for 50% autoinduction (mcg/mL)
F_CBZ     : 0.75  : CBZ oral bioavailability
fm_CBZ    : 0.40  : Fraction of CBZ converted to CBZ-10,11-epoxide
CL_CBZE   : 2.00  : CBZ-epoxide elimination clearance (L/h)
Vc_CBZE   : 20.0  : CBZ-epoxide distribution volume (L)

// === LTG Pharmacokinetics (UGT1A4 metabolism, DDI affected) ===
ka_LTG    : 0.80  : LTG absorption rate constant (1/h)
Vc_LTG    : 77.0  : LTG central volume (L) [1.1 L/kg]
CL_LTG0   : 1.50  : LTG baseline clearance (L/h) [t1/2~36h, monotherapy]
F_LTG     : 0.98  : LTG oral bioavailability

// === Drug-Drug Interaction Flags ===
DDI_VPA   : 0.0   : 1 = VPA co-administered (halves LTG clearance)
DDI_CBZ   : 0.0   : 1 = CBZ co-administered (doubles LTG clearance)

// === PD: GABA Indirect Response Model ===
kin_GABA  : 0.100 : GABA zero-order synthesis rate (nmol/mg/h, normalized=1)
kout_GABA : 0.100 : GABA first-order elimination rate (1/h)
GABA0     : 1.000 : Baseline brain GABA (normalized)
IC50_VPA_GABA : 50.0 : VPA IC50 for GABA-T inhibition (mcg/mL)
Imax_VPA_GABA : 0.55 : VPA maximum GABA-T inhibition fraction
n_VPA_GABA    : 1.20 : Hill coefficient (VPA-GABA)

// === PD: Synaptic Glutamate (indirect response) ===
kin_SYNAP  : 0.100 : Synaptic Glu synthesis rate (normalized)
kout_SYNAP : 0.100 : Synaptic Glu clearance rate (1/h)
SYNAP0     : 1.000 : Baseline synaptic glutamate (normalized)
IC50_LEV_SYNAP : 5.0 : LEV IC50 for vesicle-Glu release suppression (mcg/mL)
Imax_LEV_SYNAP : 0.45: LEV maximum Glu release inhibition

// === PD: SV2A Occupancy (LEV) ===
Kd_SV2A   : 4.0   : LEV Kd for SV2A binding (mcg/mL)
kon_SV2A  : 0.50  : SV2A association rate constant (1/h per mcg/mL)
koff_SV2A : 2.00  : SV2A dissociation rate constant (1/h)

// === PD: Sodium Channel Blockade (CBZ + LTG) ===
IC50_CBZ_Nav : 3.5  : CBZ IC50 for Nav slow inactivation (mcg/mL)
IC50_LTG_Nav : 2.0  : LTG IC50 for Nav state-dep block (mcg/mL)
IC50_PHT_Nav : 2.5  : PHT IC50 for Nav inactivation (mcg/mL)
Emax_Nav     : 0.90 : Maximum Nav channel blocked fraction
n_Nav        : 1.50 : Hill coefficient for Nav blockade

// === Seizure Threshold & Frequency ===
STHRES0     : 1.000 : Baseline seizure threshold (normalized)
kthres_rec  : 0.020 : Seizure threshold recovery rate constant (1/h)
alpha_GABA  : 0.350 : GABA contribution to threshold elevation
alpha_Nav   : 0.500 : Nav blockade contribution to threshold
alpha_SV2A  : 0.250 : SV2A occupancy contribution to threshold
SeizBasal   : 8.00  : Baseline seizure frequency (per 28 days)
k_seiz      : 12.0  : Seizure frequency-threshold relationship (steepness)

// === P-gp Drug Resistance ===
PGP0       : 1.000 : Baseline P-gp expression (normalized to 1)
k_PGP_ind  : 0.003 : P-gp induction rate constant (1/h per seizure burden)
k_PGP_deg  : 0.050 : P-gp degradation rate constant (1/h)
pgp_CBZ    : 0.55  : P-gp effect coefficient on CBZ CNS concentration
pgp_VPA    : 0.40  : P-gp effect coefficient on VPA CNS concentration
pgp_LTG    : 0.35  : P-gp effect coefficient on LTG CNS concentration
pgp_LEV    : 0.15  : P-gp effect coefficient on LEV CNS concentration

// === PHT/GBP/PER dosing flags (binary, simplified) ===
PHT_DOSE   : 0.0   : PHT plasma concentration if used (mcg/mL, user sets)
GBP_DOSE   : 0.0   : GBP concentration proxy (mcg/mL)
BZD_BOLUS  : 0.0   : BZD effect flag: 1 = IV bolus administered

// === mTOR pathway (TSC/FCD) ===
mTOR_activ : 1.0   : mTOR pathway activity (1=normal, >1=hyperactivated in TSC)
Emax_mTOR  : 0.50  : Max seizure threshold increase with everolimus
EC50_mTOR  : 0.5   : Everolimus EC50 for mTOR inhibition (normalized dose)
ever_dose  : 0.0   : Everolimus relative dose (0–1, 1=full inhibition)

$CMT @annotated
AGUT    : VPA gut compartment (mg)
ACENT   : VPA central compartment (mg)
APER    : VPA peripheral compartment (mg)
BGUT    : LEV gut compartment (mg)
BCENT   : LEV central compartment (mg)
CGUT    : CBZ gut compartment (mg)
CCENT   : CBZ central compartment (mg)
CMETA   : CBZ-10,11-epoxide metabolite (mg)
DGUT    : LTG gut compartment (mg)
DCENT   : LTG central compartment (mg)
GABA    : Brain GABA level — indirect response (normalized)
SYNAP   : Synaptic glutamate — indirect response (normalized)
SV2A_OCC : SV2A occupancy fraction by LEV (0–1)
NAV_BLOCK : Sodium channel blocked fraction combined (0–1)
STHRES  : Seizure threshold (normalized)
PGP     : P-glycoprotein expression (normalized)

$GLOBAL
// Convenience macros for plasma concentrations (mcg/mL)
#define C_VPA   (ACENT / Vc_VPA)
#define C_LEV   (BCENT / Vc_LEV)
#define C_CBZ   (CCENT / Vc_CBZ)
#define C_CBZE  (CMETA / Vc_CBZE)
#define C_LTG   (DCENT / Vc_LTG)
#define fmaxd(a,b) ((a)>(b)?(a):(b))

$MAIN
// ── LTG clearance with DDI adjustment ──
double CL_LTG = CL_LTG0;
if (DDI_VPA > 0.5) CL_LTG = CL_LTG0 * 0.50;  // VPA inhibits UGT1A4: t1/2 doubled
if (DDI_CBZ > 0.5) CL_LTG = CL_LTG  * 2.00;  // CBZ induces CYP3A4: CL doubled

// ── CBZ autoinduction: Michaelis-Menten saturable induction ──
double CL_CBZ = CL_CBZ0 + (CL_CBZ_max - CL_CBZ0) * C_CBZ / (EC50_auto + C_CBZ);

// ── P-gp corrected CNS drug concentrations ──
double pgp_factor = fmaxd(1.0, PGP);  // PGP >= 1 always
double CNS_VPA = (fu_VPA * C_VPA) / (1.0 + pgp_VPA * (pgp_factor - 1.0));
double CNS_LEV = C_LEV               / (1.0 + pgp_LEV * (pgp_factor - 1.0));
double CNS_CBZ = C_CBZ               / (1.0 + pgp_CBZ * (pgp_factor - 1.0));
double CNS_LTG = C_LTG               / (1.0 + pgp_LTG * (pgp_factor - 1.0));
double CNS_CBZE= C_CBZE              / (1.0 + pgp_CBZ * (pgp_factor - 1.0));

// ── GABA-T inhibition by VPA (Hill equation) ──
double VPA_pow = pow(CNS_VPA, n_VPA_GABA);
double IC50_pow= pow(IC50_VPA_GABA, n_VPA_GABA);
double Imax_GABA_eff = Imax_VPA_GABA * VPA_pow / (IC50_pow + VPA_pow);

// ── Synaptic Glu suppression by LEV (SV2A mechanism) ──
double Imax_SYNAP_eff = Imax_LEV_SYNAP * CNS_LEV / (IC50_LEV_SYNAP + CNS_LEV);

// ── Nav channel blockade (CBZ + CBZ-E + LTG combined, non-redundant) ──
double Nav_CBZ_eff = Emax_Nav * pow(CNS_CBZ + 0.5*CNS_CBZE, n_Nav) /
                     (pow(IC50_CBZ_Nav, n_Nav) + pow(CNS_CBZ + 0.5*CNS_CBZE, n_Nav));
double Nav_LTG_eff = Emax_Nav * pow(CNS_LTG, n_Nav) /
                     (pow(IC50_LTG_Nav, n_Nav) + pow(CNS_LTG, n_Nav));
double Nav_PHT_eff = Emax_Nav * PHT_DOSE / (IC50_PHT_Nav + PHT_DOSE);
double Nav_combined = 1.0 - (1.0 - Nav_CBZ_eff) * (1.0 - Nav_LTG_eff) * (1.0 - Nav_PHT_eff);

// ── BZD acute GABA-A potentiation (SE rescue) ──
double BZD_GABA_boost = (BZD_BOLUS > 0.5) ? 0.60 : 0.0;  // 60% GABA-A activity boost

// ── mTOR pathway: everolimus reduces threshold in TSC/FCD ──
double mTOR_effect = Emax_mTOR * ever_dose / (EC50_mTOR + ever_dose);
double mTOR_thresh_adj = (mTOR_activ > 1.0) ? mTOR_effect * (mTOR_activ - 1.0) : 0.0;

// ── Seizure frequency (episodes/28 days) ──
double SeizFreq = SeizBasal * exp(-k_seiz * (STHRES - STHRES0));
if (SeizFreq < 0.0) SeizFreq = 0.0;

$ODE
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// VPA Pharmacokinetics (2-compartment model)
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
dxdt_AGUT  = -ka_VPA * AGUT;
dxdt_ACENT =  ka_VPA * AGUT
              - (CL_VPA / Vc_VPA) * ACENT
              - (Q_VPA  / Vc_VPA) * ACENT
              + (Q_VPA  / Vp_VPA) * APER;
dxdt_APER  =  (Q_VPA  / Vc_VPA) * ACENT
              - (Q_VPA  / Vp_VPA) * APER;

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// LEV Pharmacokinetics (1-compartment, renal)
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
dxdt_BGUT  = -ka_LEV * BGUT;
dxdt_BCENT =  ka_LEV * BGUT - (CL_LEV / Vc_LEV) * BCENT;

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// CBZ PK (1-comp + autoinduction + CBZ-E)
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
dxdt_CGUT  = -ka_CBZ * CGUT;
dxdt_CCENT =  ka_CBZ * CGUT - (CL_CBZ / Vc_CBZ) * CCENT;
dxdt_CMETA =  fm_CBZ * (CL_CBZ / Vc_CBZ) * CCENT
              - (CL_CBZE / Vc_CBZE) * CMETA;

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// LTG PK (1-comp, UGT1A4, DDI-sensitive)
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
dxdt_DGUT  = -ka_LTG * DGUT;
dxdt_DCENT =  ka_LTG * DGUT - (CL_LTG / Vc_LTG) * DCENT;

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// GABA indirect response (VPA → GABA-T inhib → GABA↑)
// dGABA/dt = kin × (1 + Imax_GABA_eff + BZD boost) - kout × GABA
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
dxdt_GABA  =  kin_GABA * (1.0 + Imax_GABA_eff + BZD_GABA_boost)
              - kout_GABA * GABA;

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// Synaptic glutamate indirect response (LEV → SV2A → Glu release↓)
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
dxdt_SYNAP =  kin_SYNAP * (1.0 - Imax_SYNAP_eff)
              - kout_SYNAP * SYNAP;

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// SV2A occupancy (fast equilibrium approximation)
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
double SV2A_eq = CNS_LEV / (Kd_SV2A + CNS_LEV);
dxdt_SV2A_OCC = kon_SV2A * CNS_LEV * (1.0 - SV2A_OCC)
                - koff_SV2A * SV2A_OCC;

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// Nav channel blocked fraction (fast equilibrium)
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
dxdt_NAV_BLOCK = 5.0 * (Nav_combined - NAV_BLOCK);  // fast equil, t1/2~0.14h

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// Seizure threshold dynamic model
// dSThres/dt = kthres_rec × [STHRES_target - STHRES]
// STHRES_target = STHRES0 + contributions from GABA, Nav, SV2A, mTOR
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
double STHRES_target = STHRES0
                       + alpha_GABA * (GABA - GABA0)
                       + alpha_Nav  * NAV_BLOCK
                       + alpha_SV2A * SV2A_OCC
                       + mTOR_thresh_adj;
dxdt_STHRES = kthres_rec * (STHRES_target - STHRES);

// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
// P-glycoprotein expression (drug resistance)
// Driven by ongoing seizure burden; degrades toward baseline
// ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
double PGP_drive = k_PGP_ind * SeizFreq;
dxdt_PGP = PGP_drive - k_PGP_deg * (PGP - PGP0);

$TABLE
capture C_VPA_mcg   = C_VPA;
capture C_LEV_mcg   = C_LEV;
capture C_CBZ_mcg   = C_CBZ;
capture C_CBZE_mcg  = C_CBZE;
capture C_LTG_mcg   = C_LTG;
capture CNS_VPA_f   = CNS_VPA;
capture CNS_CBZ_f   = CNS_CBZ;
capture GABA_norm   = GABA;
capture SYNAP_norm  = SYNAP;
capture SV2A_frac   = SV2A_OCC;
capture Nav_frac    = NAV_BLOCK;
capture Thresh      = STHRES;
capture PGP_exp     = PGP;
capture SeizFreq_obs = SeizFreq;
capture Responder   = (SeizFreq < SeizBasal * 0.50) ? 1.0 : 0.0;
capture SeizFree    = (SeizFreq < 0.1) ? 1.0 : 0.0;
capture CL_LTG_eff  = CL_LTG;
capture CL_CBZ_eff  = CL_CBZ;
'

mod <- mcode("epilepsy_qsp", epi_code)

# ─── Helper: build event table for chronic dosing ─────────────────────────────
make_events <- function(dose_VPA=0, dose_LEV=0, dose_CBZ=0, dose_LTG=0,
                        freq_VPA=2, freq_LEV=2, freq_CBZ=2, freq_LTG=2,
                        duration_days=180) {
  e <- ev_seq()
  if (dose_VPA > 0) {
    e <- e + ev(cmt="AGUT", amt=dose_VPA/freq_VPA, ii=24/freq_VPA, addl=duration_days*freq_VPA-1)
  }
  if (dose_LEV > 0) {
    e <- e + ev(cmt="BGUT", amt=dose_LEV/freq_LEV, ii=24/freq_LEV, addl=duration_days*freq_LEV-1)
  }
  if (dose_CBZ > 0) {
    e <- e + ev(cmt="CGUT", amt=dose_CBZ/freq_CBZ, ii=24/freq_CBZ, addl=duration_days*freq_CBZ-1)
  }
  if (dose_LTG > 0) {
    e <- e + ev(cmt="DGUT", amt=dose_LTG/freq_LTG, ii=24/freq_LTG, addl=duration_days*freq_LTG-1)
  }
  e
}

# ─── Initial conditions ────────────────────────────────────────────────────────
init_epi <- init(mod, GABA=1.0, SYNAP=1.0, STHRES=1.0, PGP=1.0,
                 SV2A_OCC=0, NAV_BLOCK=0)

sim_times <- seq(0, 4320, by=1)  # 180 days in hours

# ═══════════════════════════════════════════════════════════════════════════════
# SCENARIO 1: Untreated (Baseline)
# ═══════════════════════════════════════════════════════════════════════════════
out_s1 <- mrgsim(init_epi, events=ev(), end=4320, delta=1,
                 param(mod, SeizBasal=8)) %>% as_tibble() %>% mutate(scenario="Untreated")

# ═══════════════════════════════════════════════════════════════════════════════
# SCENARIO 2: VPA Monotherapy 1,000 mg/day BID
# ═══════════════════════════════════════════════════════════════════════════════
ev_s2 <- ev(cmt="AGUT", amt=500, ii=12, addl=359)  # 500 mg BID × 360 doses
out_s2 <- mrgsim(init_epi, events=ev_s2, end=4320, delta=1) %>%
  as_tibble() %>% mutate(scenario="VPA 1,000 mg/day")

# ═══════════════════════════════════════════════════════════════════════════════
# SCENARIO 3: LEV Monotherapy 3,000 mg/day BID
# ═══════════════════════════════════════════════════════════════════════════════
ev_s3 <- ev(cmt="BGUT", amt=1500, ii=12, addl=359)
out_s3 <- mrgsim(init_epi, events=ev_s3, end=4320, delta=1) %>%
  as_tibble() %>% mutate(scenario="LEV 3,000 mg/day")

# ═══════════════════════════════════════════════════════════════════════════════
# SCENARIO 4: CBZ Monotherapy 600 mg/day BID
# ═══════════════════════════════════════════════════════════════════════════════
ev_s4 <- ev(cmt="CGUT", amt=300, ii=12, addl=359)
out_s4 <- mrgsim(init_epi, events=ev_s4, end=4320, delta=1) %>%
  as_tibble() %>% mutate(scenario="CBZ 600 mg/day")

# ═══════════════════════════════════════════════════════════════════════════════
# SCENARIO 5: LTG Monotherapy 200 mg/day BID
# ═══════════════════════════════════════════════════════════════════════════════
ev_s5 <- ev(cmt="DGUT", amt=100, ii=12, addl=359)
out_s5 <- mrgsim(init_epi, events=ev_s5, end=4320, delta=1) %>%
  as_tibble() %>% mutate(scenario="LTG 200 mg/day")

# ═══════════════════════════════════════════════════════════════════════════════
# SCENARIO 6: VPA 500 mg/day + LTG 100 mg/day (DDI: VPA halves LTG CL)
# ═══════════════════════════════════════════════════════════════════════════════
ev_s6 <- ev(cmt="AGUT", amt=250, ii=12, addl=359) +
         ev(cmt="DGUT", amt=50,  ii=12, addl=359)
out_s6 <- mrgsim(init_epi, events=ev_s6, end=4320, delta=1,
                 param(mod, DDI_VPA=1)) %>%
  as_tibble() %>% mutate(scenario="VPA+LTG (DDI: LTG t½×2)")

# ═══════════════════════════════════════════════════════════════════════════════
# SCENARIO 7: CBZ 600 mg/day + LTG 400 mg/day (DDI: CBZ doubles LTG CL)
# ═══════════════════════════════════════════════════════════════════════════════
ev_s7 <- ev(cmt="CGUT", amt=300, ii=12, addl=359) +
         ev(cmt="DGUT", amt=200, ii=12, addl=359)
out_s7 <- mrgsim(init_epi, events=ev_s7, end=4320, delta=1,
                 param(mod, DDI_CBZ=1)) %>%
  as_tibble() %>% mutate(scenario="CBZ+LTG (DDI: LTG CL×2)")

# ═══════════════════════════════════════════════════════════════════════════════
# SCENARIO 8: Drug-Resistant Epilepsy (P-gp 3× baseline)
# ═══════════════════════════════════════════════════════════════════════════════
init_dre <- init(mod, GABA=1.0, SYNAP=1.0, STHRES=1.0, PGP=3.0,
                 SV2A_OCC=0, NAV_BLOCK=0)
out_s8 <- mrgsim(init_dre, events=ev_s2, end=4320, delta=1) %>%
  as_tibble() %>% mutate(scenario="DRE: VPA (P-gp 3×)")

# ═══════════════════════════════════════════════════════════════════════════════
# SCENARIO 9: Status Epilepticus Rescue (IV BZD at t=24h)
# ═══════════════════════════════════════════════════════════════════════════════
# Simulate 24h of SE then BZD bolus (modeled as BZD_BOLUS=1 parameter change)
init_se <- init(mod, GABA=0.5, SYNAP=2.0, STHRES=0.3, PGP=1.0,
                SV2A_OCC=0, NAV_BLOCK=0)
ev_se_bzd <- ev(time=24, cmt="BGUT", amt=0)  # placeholder; BZD via param
# Two phases: pre-BZD (low threshold) and post-BZD
out_pre  <- mrgsim(init_se, end=24, delta=0.5) %>% as_tibble()
init_post <- as.list(out_pre[nrow(out_pre), c("AGUT","ACENT","APER","BGUT","BCENT",
                                               "CGUT","CCENT","CMETA","DGUT","DCENT",
                                               "GABA","SYNAP","SV2A_OCC","NAV_BLOCK",
                                               "STHRES","PGP")])
out_post <- mrgsim(do.call(init, c(list(mod), init_post)),
                   param(mod, BZD_BOLUS=1),
                   end=120, delta=0.5) %>% as_tibble()
out_pre$time  <- out_pre$time
out_post$time <- out_post$time + 24
out_s9 <- bind_rows(
  out_pre  %>% mutate(scenario="SE: Pre-BZD rescue"),
  out_post %>% mutate(scenario="SE: Post-IV BZD (t+24h)")
)

# ═══════════════════════════════════════════════════════════════════════════════
# SCENARIO 10: TSC + Everolimus (mTOR hyperactivation)
# ═══════════════════════════════════════════════════════════════════════════════
out_s10_pre <- mrgsim(init_epi, events=ev_s2, end=4320, delta=1,
                      param(mod, mTOR_activ=2.0, ever_dose=0.0)) %>%
  as_tibble() %>% mutate(scenario="TSC without everolimus")
out_s10_post <- mrgsim(init_epi, events=ev_s2, end=4320, delta=1,
                       param(mod, mTOR_activ=2.0, ever_dose=0.8)) %>%
  as_tibble() %>% mutate(scenario="TSC + Everolimus (mTOR↓)")

# ─── Combine Core Scenarios 1–8 ───────────────────────────────────────────────
results_all <- bind_rows(out_s1, out_s2, out_s3, out_s4,
                         out_s5, out_s6, out_s7, out_s8) %>%
  mutate(day = time / 24,
         scenario = factor(scenario, levels=c(
           "Untreated",
           "VPA 1,000 mg/day",
           "LEV 3,000 mg/day",
           "CBZ 600 mg/day",
           "LTG 200 mg/day",
           "VPA+LTG (DDI: LTG t½×2)",
           "CBZ+LTG (DDI: LTG CL×2)",
           "DRE: VPA (P-gp 3×)"
         )))

# ─── Plot 1: Seizure Frequency over time ──────────────────────────────────────
p1 <- results_all %>%
  filter(day <= 180) %>%
  ggplot(aes(x=day, y=SeizFreq_obs, color=scenario)) +
  geom_line(linewidth=0.8) +
  geom_hline(yintercept=4, linetype="dashed", color="gray50", alpha=0.6) +  # 50% reduction line
  annotate("text", x=170, y=4.3, label="50% reduction threshold", size=2.8, color="gray40") +
  scale_color_brewer(palette="Dark2") +
  labs(title="Seizure Frequency Over 6 Months",
       subtitle="Epilepsy QSP Model — 8 Treatment Scenarios",
       x="Time (days)", y="Seizure frequency (episodes/28 days)",
       color="Scenario") +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", legend.text=element_text(size=7)) +
  guides(color=guide_legend(nrow=3))

# ─── Plot 2: Seizure Threshold Dynamics ───────────────────────────────────────
p2 <- results_all %>%
  filter(day <= 90) %>%
  ggplot(aes(x=day, y=Thresh, color=scenario)) +
  geom_line(linewidth=0.8) +
  geom_hline(yintercept=1.0, linetype="dotted", color="gray50") +
  scale_color_brewer(palette="Dark2") +
  labs(title="Seizure Threshold (STHRES) — First 90 Days",
       x="Time (days)", y="Threshold (normalized)",
       color="Scenario") +
  theme_bw(base_size=11) +
  theme(legend.position="none")

# ─── Plot 3: Brain GABA Level ─────────────────────────────────────────────────
p3 <- results_all %>%
  filter(day <= 90) %>%
  ggplot(aes(x=day, y=GABA_norm, color=scenario)) +
  geom_line(linewidth=0.8) +
  scale_color_brewer(palette="Dark2") +
  labs(title="Brain GABA Level", x="Time (days)",
       y="GABA (normalized)", color="Scenario") +
  theme_bw(base_size=11) +
  theme(legend.position="none")

# ─── Plot 4: Nav Channel Blockade ─────────────────────────────────────────────
p4 <- results_all %>%
  filter(day <= 90) %>%
  ggplot(aes(x=day, y=Nav_frac, color=scenario)) +
  geom_line(linewidth=0.8) +
  scale_y_continuous(labels=scales::percent_format()) +
  scale_color_brewer(palette="Dark2") +
  labs(title="Sodium Channel Blockade (CBZ+LTG+PHT)",
       x="Time (days)", y="Nav blocked fraction (%)",
       color="Scenario") +
  theme_bw(base_size=11) +
  theme(legend.position="none")

# ─── Plot 5: P-gp Expression (Drug Resistance) ────────────────────────────────
p5 <- results_all %>%
  filter(day <= 180) %>%
  ggplot(aes(x=day, y=PGP_exp, color=scenario)) +
  geom_line(linewidth=0.8) +
  geom_hline(yintercept=3, linetype="dashed", color="red", alpha=0.5) +
  annotate("text", x=150, y=3.2, label="DRE threshold (P-gp 3×)",
           size=2.8, color="red") +
  scale_color_brewer(palette="Dark2") +
  labs(title="P-glycoprotein Expression (Drug Resistance Marker)",
       x="Time (days)", y="P-gp expression (normalized)",
       color="Scenario") +
  theme_bw(base_size=11) +
  theme(legend.position="none")

# ─── Plot 6: PK — VPA steady-state ────────────────────────────────────────────
p6 <- results_all %>%
  filter(scenario %in% c("VPA 1,000 mg/day","VPA+LTG (DDI: LTG t½×2)")) %>%
  filter(day >= 14 & day <= 16) %>%
  ggplot(aes(x=(time - 336), y=C_VPA_mcg, color=scenario)) +
  geom_line(linewidth=1.2) +
  labs(title="VPA Plasma Concentration — Day 15 (Steady-State)",
       x="Time from last dose (h)", y="VPA plasma (mcg/mL)",
       color="Scenario") +
  theme_bw(base_size=11) +
  theme(legend.position="bottom")

# ─── Combined Dashboard ────────────────────────────────────────────────────────
dashboard <- (p1 + p2) / (p3 + p4) / (p5 + p6)

print(dashboard)

# ─── Summary Table ─────────────────────────────────────────────────────────────
summary_tbl <- results_all %>%
  filter(day >= 150 & day <= 180) %>%
  group_by(scenario) %>%
  summarise(
    SeizFreq_SS    = mean(SeizFreq_obs, na.rm=TRUE),
    Responder_pct  = mean(Responder, na.rm=TRUE) * 100,
    SeizFree_pct   = mean(SeizFree,  na.rm=TRUE) * 100,
    STHRES_SS      = mean(Thresh,    na.rm=TRUE),
    GABA_SS        = mean(GABA_norm, na.rm=TRUE),
    PGP_SS         = mean(PGP_exp,   na.rm=TRUE),
    Nav_block_pct  = mean(Nav_frac,  na.rm=TRUE) * 100,
    .groups="drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 2)))

cat("\n====== EPILEPSY QSP MODEL — STEADY-STATE SUMMARY (Day 150–180) ======\n")
print(as.data.frame(summary_tbl))
cat("\nNote: SeizFreq = seizures/28 days; Responder = >50% reduction from baseline (8/mo)\n")
cat("      SeizFree = <0.1 seizures/28 days; Nav_block = combined CBZ+LTG+PHT fraction\n\n")

# ─── SE Rescue Visualization ──────────────────────────────────────────────────
p_se <- out_s9 %>%
  mutate(phase = scenario) %>%
  ggplot(aes(x=time, y=Thresh, color=phase)) +
  geom_line(linewidth=1.2) +
  geom_vline(xintercept=24, linetype="dashed", color="red") +
  annotate("text", x=26, y=0.7, label="IV BZD admin", size=3, color="red") +
  labs(title="Status Epilepticus: Seizure Threshold Before & After IV BZD",
       x="Time (h)", y="Seizure threshold (normalized)",
       color="Phase") +
  theme_bw(base_size=12)

print(p_se)

# ─── TSC/Everolimus Scenario ──────────────────────────────────────────────────
p_tsc <- bind_rows(out_s10_pre, out_s10_post) %>%
  mutate(day=time/24) %>%
  ggplot(aes(x=day, y=SeizFreq_obs, color=scenario)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=c("darkorange","steelblue")) +
  labs(title="TSC: Seizure Frequency With & Without Everolimus (mTOR Inhibitor)",
       x="Time (days)", y="Seizure frequency (episodes/28 days)",
       color="Scenario") +
  theme_bw(base_size=12)

print(p_tsc)
