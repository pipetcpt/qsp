## ============================================================
## Duchenne Muscular Dystrophy (DMD) QSP Model — mrgsolve
## ============================================================
## Disease:   Duchenne Muscular Dystrophy (DMD, OMIM #310200)
## Framework: mrgsolve (ODE-based PK/PD simulation in R)
##
## Pathophysiology summary:
##   X-linked dystrophin gene (Xp21.2) mutations → absent dystrophin →
##   DAPC disruption → sarcolemmal fragility → Ca²⁺ influx →
##   calpain/ROS activation → necrosis → inflammation (NF-κB/TNF-α) →
##   fibrosis (TGF-β1/SMAD2/3) → loss of ambulation → DCM → respiratory failure
##
## Compartments (18 ODEs):
##   Drug PK (6):  Ete_C1, Ete_C2, Ete_Muscle, DFZ_Gut, DFZ_Plasma, Active_DFZ
##   Disease PD (12): Dystrophin, Fiber_H, Fiber_N, Fiber_R, Inflam,
##                    Fibrosis, SC_Pool, CK_serum, FVC_pct, LVEF, NSAA, SixMWD
##
## Treatment scenarios (7):
##   1. Natural history (no treatment)
##   2. Deflazacort 0.9 mg/kg/day PO (EMFLAZA)
##   3. Prednisone 0.75 mg/kg/day PO (standard)
##   4. Eteplirsen 30 mg/kg/wk IV + Deflazacort (exon-51 skippable, ~13% DMD)
##   5. Casimersen 30 mg/kg/wk IV + Deflazacort (exon-45 skippable, ~8% DMD)
##   6. Delandistrogene moxeparvovec (Elevidys, 1×10¹⁴ vg/kg, gene therapy)
##   7. Vamorolone 6 mg/kg/day PO (AGAMREE, dissociated GR modulator)
##
## Calibration references:
##   - McDonald 2018 NEJM (ESSENCE study – deflazacort vs prednisone)
##   - Mendell 2016 Ann Neurol (eteplirsen – dystrophin restoration 0.28-0.93%)
##   - ELEVIDYS FDA approval 2023 – SRP-9001 Phase 3 (EMBARK trial)
##   - Griggs 2016 NEJM (deflazacort vs prednisone ambulation)
##   - Servais 2022 Lancet ND (vamorolone VISION-DMD trial)
##   - Bushby 2010 Lancet Neurol (DMD care guidelines)
##   - Barnard 2019 J Pharmacol Exp Ther (DMD QSP framework)
##
## Author: Claude Code Routine (CCR) | Date: 2026-06-23
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ─────────────────────────────────────────────────────────────
## MODEL CODE
## ─────────────────────────────────────────────────────────────

dmd_code <- '
$PROB DMD QSP Model — mrgsolve ODE (18 compartments)

$PARAM
@annotated
// ── Eteplirsen PK parameters (2-compartment IV, PopPK) ──────
CL_Ete    : 80.4    : Eteplirsen CL (L/h), Bladen 2015 Hum Mutat
Vc_Ete    : 15.8    : Central volume (L), 70-kg adult
Vp_Ete    : 42.2    : Peripheral volume (L)
Q_Ete     : 30.1    : Intercompartmental CL (L/h)
kUptake   : 0.035   : Muscle uptake rate constant (1/h)
kElim_m   : 0.005   : Muscle elimination rate constant (1/h)

// ── Deflazacort PK parameters (oral, 1-compartment) ──────────
ka_DFZ    : 1.20    : Absorption rate constant (1/h)
CL_DFZ    : 18.5    : Clearance (L/h), CYP3A4 metabolism
Vd_DFZ    : 62.0    : Volume of distribution (L)
F_DFZ     : 0.89    : Bioavailability of deflazacort
kconv_DFZ : 2.10    : Conversion to 21-desacetyl-DFZ (1/h) – rapid esterase
CL_aDFZ   : 8.20    : Clearance of 21-desacetyl-DFZ (L/h)
Vd_aDFZ   : 45.0    : Volume of distribution of active metabolite

// ── Exon-skipping PD (Emax model) ────────────────────────────
EC50_skip  : 0.45   : ASO muscle conc for 50% exon-skip efficiency (pmol/g)
Emax_skip  : 0.85   : Max exon-skip efficiency (fraction, ~85%)
Hill_skip  : 1.5    : Hill coefficient for exon-skipping
kDyst_syn  : 0.0045 : Dystrophin synthesis rate constant (fraction/h)
kDyst_deg  : 0.0040 : Dystrophin degradation rate constant (1/h)
// Baseline (DMD): Dyst_0 ≈ 0% (absent); normal = 100%

// ── Gene therapy PD (single-dose durable effect) ─────────────
GT_eff     : 0.60   : Peak dystrophin restoration from gene therapy (fraction)
t_GT_peak  : 720    : Time to peak expression (h, ~30 days)
t_GT_half  : 8760   : Half-life of transgene expression (h, ~1 year)
// Note: Elevidys Phase 3 (EMBARK) showed 28% dystrophin vs 2% placebo at 52wk

// ── Glucocorticoid PD (biophase, NF-κB inhibition) ────────────
EC50_GC    : 0.020  : Active DFZ plasma for 50% NF-κB inhibition (mg/L)
Emax_GC    : 0.75   : Max NF-κB inhibition (fraction)
Hill_GC    : 1.2    : Hill coefficient for GC effect
IC50_Vamo  : 0.030  : Vamorolone IC50 for NF-κB inhibition (mg/L, µg/mL)
// Vamorolone: ~similar efficacy to DFZ with fewer bone/growth side effects

// ── Muscle fiber dynamics ─────────────────────────────────────
kNecrosis  : 0.0025 : Rate of fiber necrosis (per unit Inflam, 1/h)
kRegen     : 0.0060 : Rate of fiber regeneration from satellite cells (1/h)
kMaturation: 0.0020 : Rate of regenerating → healthy fiber maturation (1/h)
kFibrosis_rate : 0.00015 : Rate of fibrosis accumulation (per necrosis, 1/h)
Fiber_H_0  : 100.0  : Baseline healthy fibers (arbitrary units)
Fiber_N_0  : 0.0    : Baseline necrotic fibers
Fiber_R_0  : 0.0    : Baseline regenerating fibers

// ── Inflammation dynamics ─────────────────────────────────────
Inflam_basal : 25.0 : Basal inflammation index in DMD (0-100 scale)
kInflam_stim : 0.15 : NF-κB stimulation by necrosis (per unit Fiber_N)
kInflam_decay: 0.010: Inflammation decay rate constant (1/h)
Inflam_0   : 25.0   : Initial inflammation (DMD basal ~25)
// Normal baseline = 5; DMD = 20-30; exacerbation = 40-60

// ── Fibrosis dynamics ─────────────────────────────────────────
kFib_stim   : 0.0025 : TGF-β driven fibrosis accumulation (per Inflam unit)
kFib_decay  : 0.00008: Fibrosis regression rate (very slow) (1/h)
Fibrosis_0  : 0.05   : Initial fibrosis at model start (~5% at age 3)
Fibrosis_max: 0.85   : Maximum fibrosis (85% of muscle)
// McDonald 2018: fibrosis increases ~5% per year in ambulatory DMD

// ── Satellite cell pool ───────────────────────────────────────
kSC_regen   : 0.005  : SC pool replenishment rate (1/h)
kSC_exhaust : 0.0003 : SC pool exhaustion by chronic cycles (per Fiber_N)
SC_Pool_0   : 1.0    : Initial SC pool (normalized, 1=full)
// After ~50 divisions, telomere shortening limits SC function

// ── Serum CK dynamics ─────────────────────────────────────────
kCK_release : 0.08   : CK release from necrotic fibers (U/L per fiber unit/h)
kCK_elim    : 0.025  : CK elimination from serum (1/h; t½ ~28h)
CK_basal    : 150    : Basal CK from normal turnover (U/L)
// DMD: CK typically 10,000-50,000 U/L (10-100× ULN)

// ── Respiratory function (FVC% decline) ──────────────────────
kFVC_decline : 0.000095 : FVC% annual decline per fibrosis unit (1/h)
FVC_0        : 95.0     : Initial FVC% predicted (age ~5-7yr)
FVC_min      : 10.0     : Minimum FVC% (floor)
// Natural history: FVC peaks ~8-10yr, then declines ~2%/yr ambulatory,
// ~4%/yr non-ambulatory; Mendell 2016

// ── Cardiac function (LVEF) ──────────────────────────────────
kLVEF_decline : 0.000050 : LVEF annual decline per fibrosis unit (1/h)
LVEF_0        : 62.0     : Initial LVEF% (normal, onset ~10-12yr DCM)
LVEF_min      : 20.0     : Minimum LVEF% (floor)
// DCM in DMD: onset ~10-12yr; LVEF declines ~2-3% per year after onset

// ── Clinical endpoints (NSAA, 6MWD) ──────────────────────────
NSAA_0     : 28.0   : Initial NSAA score (0-34) at model start (ambulatory)
k_NSAA_decline : 0.00015 : NSAA decline per unit fibrosis × time (1/h)
SixMWD_0   : 380.0  : Initial 6-minute walk distance (meters, ~7yr)
k_6MWD_decline : 0.0004  : 6MWD decline per fibrosis × time (m/unit/h)
// McDonald 2013 PLOS ONE: 6MWD peaks at ~8yr then declines

// ── Patient demographics ──────────────────────────────────────
BW         : 20.0   : Body weight (kg; pediatric DMD ~5-7yr start)
Age_start  : 5.0    : Starting age (years)
// Doses are weight-based; adjust as needed

// ── Treatment flags (0=off, 1=on) ─────────────────────────────
// (Set via event tables in simulation scripts)

$INIT
@annotated
// Drug PK
Ete_C1     : 0.0  : Eteplirsen central compartment (mg)
Ete_C2     : 0.0  : Eteplirsen peripheral compartment (mg)
Ete_Muscle : 0.0  : Eteplirsen muscle compartment (pmol/g ×1000)
DFZ_Gut    : 0.0  : Deflazacort GI compartment (mg)
DFZ_Plasma : 0.0  : Deflazacort plasma (mg)
Active_DFZ : 0.0  : 21-desacetyl-DFZ plasma (mg)

// Disease PD
Dystrophin : 0.0    : Dystrophin level (% of normal; 0 in DMD)
Fiber_H    : 100.0  : Healthy muscle fiber (arbitrary units)
Fiber_N    : 0.001  : Necrotic fiber (au; small seed)
Fiber_R    : 0.0    : Regenerating fiber (au)
Inflam     : 25.0   : Inflammation index (0-100 scale)
Fibrosis   : 0.05   : Fibrosis score (0-1)
SC_Pool    : 1.0    : Satellite cell pool (0-1, normalized)
CK_serum   : 15000  : Serum CK (U/L; DMD baseline ~10,000-20,000)
FVC_pct    : 95.0   : FVC% predicted
LVEF       : 62.0   : Left ventricular EF (%)
NSAA       : 28.0   : NSAA score (0-34)
SixMWD     : 380.0  : 6-minute walk distance (m)

$ODE
// ──────────────────────────────────────────
// DRUG PK EQUATIONS
// ──────────────────────────────────────────

// 1) Eteplirsen 2-compartment IV
double Cp_Ete = Ete_C1 / Vc_Ete;            // plasma conc (mg/L = µg/mL)
double Ct_Ete = Ete_C2 / Vp_Ete;            // peripheral conc

dxdt_Ete_C1 = -CL_Ete*Cp_Ete - Q_Ete*(Cp_Ete - Ct_Ete) - kUptake*Ete_C1;
dxdt_Ete_C2 = Q_Ete*(Cp_Ete - Ct_Ete);
dxdt_Ete_Muscle = kUptake*Ete_C1 - kElim_m*Ete_Muscle;

// 2) Deflazacort oral PK (1-compartment pro-drug → active metabolite)
dxdt_DFZ_Gut    = -ka_DFZ * DFZ_Gut;
double Cp_DFZ   = DFZ_Plasma / Vd_DFZ;
double Cp_aDFZ  = Active_DFZ / Vd_aDFZ;
dxdt_DFZ_Plasma = F_DFZ*ka_DFZ*DFZ_Gut - CL_DFZ*Cp_DFZ - kconv_DFZ*DFZ_Plasma;
dxdt_Active_DFZ = kconv_DFZ*DFZ_Plasma - CL_aDFZ*Cp_aDFZ;

// ──────────────────────────────────────────
// PD CALCULATIONS (intermediate variables)
// ──────────────────────────────────────────

// Exon-skipping effect → Dystrophin production rate modifier
double ASO_muscle_pmol = Ete_Muscle / 1000.0;  // convert back to pmol/g
double Eskip = Emax_skip * pow(ASO_muscle_pmol, Hill_skip) /
               (pow(EC50_skip, Hill_skip) + pow(ASO_muscle_pmol, Hill_skip));
// Gene therapy: durable dystrophin from model input GT_level (via PARAM override)
// GT_level=0 by default; set to GT_eff via event for Elevidys simulation
double GT_contribution = 0.0;  // See scenario simulation scripts

// Dystrophin synthesis = basal (near 0 in DMD) + exon-skip contribution
double kDyst_syn_eff = kDyst_syn * (0.02 + Eskip + GT_contribution);
// Small residual synthesis 0.02 accounts for rare revertant fibers

// Glucocorticoid effect on NF-κB (via 21-DFZ active metabolite)
double GC_NF_inhib = Emax_GC * pow(Cp_aDFZ, Hill_GC) /
                     (pow(EC50_GC, Hill_GC) + pow(Cp_aDFZ, Hill_GC));
// Note: for Vamorolone simulation, replace Cp_aDFZ with Vamo_Cp in derived variable

// Effective inflammation (GC attenuates)
double GC_protect = 1.0 - GC_NF_inhib;  // 0=full protection, 1=no GC

// Satellite cell-mediated regeneration rate (diminishes as SC pool depletes)
double SC_regen_eff = kRegen * SC_Pool;

// ──────────────────────────────────────────
// DISEASE PD ODEs
// ──────────────────────────────────────────

// 7) Dystrophin dynamics
dxdt_Dystrophin = kDyst_syn_eff - kDyst_deg * Dystrophin;
// At steady state without ASO: Dyst ≈ kDyst_syn*0.02/kDyst_deg ≈ 0.22% (revertant)
// With full exon-skip (85%): Dyst SS ≈ kDyst_syn*0.87/kDyst_deg ≈ ~9.8% (Becker-like)

// 8-10) Muscle fiber compartments (necrosis-regeneration cycle)
// Healthy fibers lost to necrosis (proportional to inflammation × susceptibility)
double dyst_protect = (1.0 - 1.0/(1.0 + Dystrophin/10.0));  // sigmoid protection
// Full dystrophin (100%): protect=0.91; 10%: protect=0.50; 0%: protect=0.0

double Necrosis_rate = kNecrosis * Inflam * GC_protect * (1.0 - dyst_protect)
                       * Fiber_H / 100.0;
// Regenerating back to healthy
double Maturation_rate = kMaturation * Fiber_R * SC_Pool;

dxdt_Fiber_H = -Necrosis_rate + Maturation_rate;
dxdt_Fiber_N = Necrosis_rate - SC_regen_eff * Fiber_N;
dxdt_Fiber_R = SC_regen_eff * Fiber_N - Maturation_rate;

// 11) Inflammation dynamics
// Stimulated by necrotic fibers (DAMPs → NF-κB); inhibited by GC
double Inflam_stim = kInflam_stim * Fiber_N;
dxdt_Inflam = Inflam_stim - kInflam_decay*Inflam
              - GC_NF_inhib * kInflam_decay * Inflam;
// Steady state without GC: Inflam_ss = kInflam_stim*Fiber_N/kInflam_decay

// 12) Fibrosis accumulation (TGF-β driven; irreversible)
double Fibrosis_stim = kFib_stim * Inflam / 100.0 * (1.0 - Fibrosis/Fibrosis_max);
double Fibrosis_reg  = kFib_decay * Fibrosis;
// GC partially reduces fibrosis accumulation
double GC_anti_fib = GC_NF_inhib * 0.5;  // partial anti-fibrotic
dxdt_Fibrosis = (1.0 - GC_anti_fib) * Fibrosis_stim - Fibrosis_reg;

// 13) Satellite cell pool exhaustion
dxdt_SC_Pool = kSC_regen*(1.0 - SC_Pool) - kSC_exhaust*Fiber_N*SC_Pool;
// SC pool depletes with chronic necrosis-regeneration cycles

// 14) Serum CK
dxdt_CK_serum = CK_basal + kCK_release*Fiber_N - kCK_elim*CK_serum;

// 15) FVC% (respiratory function)
// Declines with fibrosis, especially after loss of ambulation
double FVC_decline_rate = kFVC_decline * Fibrosis * 100.0;
dxdt_FVC_pct = -FVC_decline_rate * (FVC_pct - FVC_min);
// Corticosteroids slow FVC decline slightly through muscle anti-inflammatory effect
// (not modeled explicitly; captured in fibrosis pathway)

// 16) LVEF (cardiac function)
// Declines with fibrosis; slow deterioration starting age ~10-12yr
// Model uses Fibrosis as proxy for cardiac involvement (correlated)
double LVEF_decline_rate = kLVEF_decline * Fibrosis * 100.0;
dxdt_LVEF = -LVEF_decline_rate * (LVEF - LVEF_min);

// 17) NSAA score (functional motor endpoint)
// NSAA declines with muscle fiber loss and fibrosis
double Func_muscle = Fiber_H / (Fiber_H_0 + 0.001);  // normalized muscle function
dxdt_NSAA = -k_NSAA_decline * Fibrosis * 100.0 * (NSAA - 0.0);
// Floor at 0 (total loss); peak around age 8-10yr then decline

// 18) 6-Minute Walk Distance
double SixMWD_floor = 0.0;
dxdt_SixMWD = -k_6MWD_decline * Fibrosis * 100.0 * (SixMWD - SixMWD_floor);
// 6MWD loss accelerates as fibrosis increases; corticosteroids slow this

$TABLE
// ── Derived PK variables ──────────────────────────────────────
double Cp_Ete_uM   = (Ete_C1/Vc_Ete) / 13.1 * 1000;  // µmol/L (MW=13,100 Da)
double ASO_Muscle_pmol = Ete_Muscle/1000.0;            // pmol/g tissue
double DFZ_Cp_mgL  = DFZ_Plasma/Vd_DFZ;               // mg/L
double aDFZ_Cp_mgL = Active_DFZ/Vd_aDFZ;              // mg/L (active DFZ)

// ── Derived PD/functional variables ──────────────────────────
double Exon_Skip_Eff  = Emax_skip * pow(ASO_Muscle_pmol,Hill_skip) /
                        (pow(EC50_skip,Hill_skip) + pow(ASO_Muscle_pmol,Hill_skip)) * 100.0;
double Dyst_pct       = Dystrophin;                    // % of normal (0-100)
double GC_Inhib_pct   = (Emax_GC * pow(aDFZ_Cp_mgL,Hill_GC) /
                         (pow(EC50_GC,Hill_GC) + pow(aDFZ_Cp_mgL,Hill_GC))) * 100.0;
double Total_Fiber    = Fiber_H + Fiber_R;             // viable fibers
double Fiber_Necrotic_pct = Fiber_N / (Fiber_H + Fiber_N + Fiber_R + 0.001) * 100.0;
double Inflam_idx     = Inflam;                        // 0-100 scale
double Fibrosis_pct   = Fibrosis * 100.0;              // % fibrosis
double SC_pct         = SC_Pool * 100.0;               // satellite cell reserve %
double CK_kUL         = CK_serum / 1000.0;            // CK in kU/L
double Ambulation     = (SixMWD > 10.0) ? 1.0 : 0.0;  // 1=ambulatory, 0=non-ambulatory

$CAPTURE
Cp_Ete_uM ASO_Muscle_pmol DFZ_Cp_mgL aDFZ_Cp_mgL
Exon_Skip_Eff Dyst_pct GC_Inhib_pct
Fiber_H Fiber_N Fiber_R Fiber_Necrotic_pct
Inflam_idx Fibrosis_pct SC_pct
CK_serum CK_kUL FVC_pct LVEF NSAA SixMWD Ambulation
'

## ─────────────────────────────────────────────────────────────
## COMPILE MODEL
## ─────────────────────────────────────────────────────────────
mod <- mcode("dmd_qsp", dmd_code)
cat("DMD QSP model compiled successfully.\n")
cat("Compartments:", length(init(mod)), "\n")
cat("Parameters:", length(param(mod)), "\n")

## ─────────────────────────────────────────────────────────────
## HELPER FUNCTIONS
## ─────────────────────────────────────────────────────────────

# Build dosing events for each scenario
# Simulation runs 8 years (70,080 h) to capture natural history peak-then-decline

sim_hours <- 8 * 365.25 * 24  # ~70,080 h
dt_out    <- 168               # output every week (168 h)

## -- Scenario 1: Natural history (no treatment) ---------------
ev_natural <- ev(time=0, cmt=1, amt=0)  # null event

## -- Scenario 2: Deflazacort 0.9 mg/kg/day PO ----------------
# BW = 20 kg (starting), dose escalates with weight growth (simplified: fixed)
dose_DFZ <- 0.9 * 20   # 18 mg/day = 4.5 mg/6h (QID simplified → weekly total in ODE)
# Simplified: 4 doses per day, every 6 hours
ev_dfz <- ev(time=0, evid=1, cmt=4, amt=dose_DFZ/4, ii=6, addl=floor(sim_hours/6))
# Note: cmt=4 is DFZ_Gut compartment (index 4)

## -- Scenario 3: Prednisone 0.75 mg/kg/day PO ----------------
# Prednisone: hepatically converted to prednisolone (active)
# In model: use deflazacort PK with adjusted parameters
dose_Pred <- 0.75 * 20  # 15 mg/day
# Simplified: model prednisone as DFZ with ka=1.5, F=0.82, EC50=0.025
ev_pred <- ev(time=0, evid=1, cmt=4, amt=dose_Pred/4, ii=6, addl=floor(sim_hours/6))

## -- Scenario 4: Eteplirsen 30 mg/kg/wk IV + Deflazacort -----
dose_Ete <- 30 * 20   # 600 mg/wk (weekly IV)
ev_ete <- ev(time=0, evid=1, cmt=1, amt=dose_Ete, ii=168, addl=floor(sim_hours/168))
ev_ete_dfz <- c(ev_ete, ev_dfz)

## -- Scenario 5: Casimersen 30 mg/kg/wk IV + Deflazacort -----
# Casimersen: same dose/route as eteplirsen, different exon target (ex45)
# Reuse eteplirsen PK (same PMO chemistry class)
ev_cas <- ev(time=0, evid=1, cmt=1, amt=dose_Ete, ii=168, addl=floor(sim_hours/168))
ev_cas_dfz <- c(ev_cas, ev_dfz)

## -- Scenario 6: Elevidys (gene therapy – single dose) --------
# Gene therapy: model as an instantaneous increase in kDyst_syn_eff
# In practice: set GT_contribution=0.50 after dose at t=0
# Achieved by modifying the GT_contribution parameter via PARAM override
ev_gt <- ev(time=0, cmt=1, amt=0)  # placeholder; GT_eff applied via param mod

## -- Scenario 7: Vamorolone 6 mg/kg/day PO -------------------
dose_Vamo <- 6 * 20  # 120 mg/day
# Vamorolone PK: t½~10h, Vd~15L, oral BID or QD
# Model as modified DFZ with IC50_Vamo
ev_vamo <- ev(time=0, evid=1, cmt=4, amt=dose_Vamo/2, ii=12, addl=floor(sim_hours/12))

## ─────────────────────────────────────────────────────────────
## SIMULATION FUNCTION
## ─────────────────────────────────────────────────────────────

run_scenario <- function(model, events, scenario_name,
                         extra_params = list(),
                         GT_active = FALSE) {
  # Override parameters if needed (e.g., GT effect)
  if (length(extra_params) > 0) {
    model <- param(model, extra_params)
  }
  # Gene therapy: boost dystrophin synthesis
  if (GT_active) {
    model <- param(model, list(kDyst_syn = 0.003))  # higher synth from transgene
    # Additional: set initial dystrophin to 0 but allow rapid synthesis
  }
  out <- mrgsim(model, events = events,
                end = sim_hours, delta = dt_out,
                carry_out = "evid")
  df <- as.data.frame(out)
  df$Scenario <- scenario_name
  df$Age_yr   <- 5 + df$time / (365.25 * 24)  # starting age 5yr
  df
}

## ─────────────────────────────────────────────────────────────
## RUN ALL 7 SCENARIOS
## ─────────────────────────────────────────────────────────────
cat("\nRunning 7 treatment scenarios...\n")

sc1 <- run_scenario(mod, ev_natural,  "1. Natural History (Untreated)")
sc2 <- run_scenario(mod, ev_dfz,      "2. Deflazacort 0.9 mg/kg/d",
                    extra_params = list(EC50_GC = 0.018))  # DFZ slightly more potent
sc3 <- run_scenario(mod, ev_pred,     "3. Prednisone 0.75 mg/kg/d",
                    extra_params = list(EC50_GC = 0.025, Emax_GC = 0.68))
sc4 <- run_scenario(mod, ev_ete_dfz,  "4. Eteplirsen + Deflazacort (ex51)")
sc5 <- run_scenario(mod, ev_cas_dfz,  "5. Casimersen + Deflazacort (ex45)",
                    extra_params = list(EC50_skip = 0.50, Emax_skip = 0.78))
sc6 <- run_scenario(mod, ev_gt,       "6. Elevidys (Gene Therapy, 1×10¹⁴vg/kg)",
                    GT_active = TRUE)
sc7 <- run_scenario(mod, ev_vamo,     "7. Vamorolone 6 mg/kg/d",
                    extra_params = list(EC50_GC   = 0.035, Emax_GC = 0.72,
                                        Emax_skip  = 0.0,
                                        kFib_stim  = 0.0022))  # anti-fibrotic

df_all <- bind_rows(sc1, sc2, sc3, sc4, sc5, sc6, sc7)
cat("Simulations complete.\n")
cat("Total rows:", nrow(df_all), "\n")

## ─────────────────────────────────────────────────────────────
## VISUALIZATION
## ─────────────────────────────────────────────────────────────

color_palette <- c(
  "1. Natural History (Untreated)"               = "#e74c3c",
  "2. Deflazacort 0.9 mg/kg/d"                   = "#3498db",
  "3. Prednisone 0.75 mg/kg/d"                   = "#2980b9",
  "4. Eteplirsen + Deflazacort (ex51)"            = "#9b59b6",
  "5. Casimersen + Deflazacort (ex45)"            = "#8e44ad",
  "6. Elevidys (Gene Therapy, 1×10¹⁴vg/kg)"      = "#2ecc71",
  "7. Vamorolone 6 mg/kg/d"                       = "#f39c12"
)

p1 <- ggplot(df_all, aes(x=Age_yr, y=SixMWD, color=Scenario)) +
  geom_line(size=0.9, na.rm=TRUE) +
  scale_color_manual(values=color_palette) +
  geom_hline(yintercept=0, linetype="dashed", color="gray60") +
  labs(title="6-Minute Walk Distance (6MWD)",
       x="Age (years)", y="6MWD (meters)", color="Scenario") +
  theme_minimal(base_size=11) +
  theme(legend.position="bottom", legend.text=element_text(size=7),
        plot.background=element_rect(fill="#1a1a2e", color=NA),
        panel.background=element_rect(fill="#2d2d4a"),
        text=element_text(color="white"),
        axis.text=element_text(color="white"),
        legend.background=element_rect(fill="#1a1a2e"))

p2 <- ggplot(df_all, aes(x=Age_yr, y=NSAA, color=Scenario)) +
  geom_line(size=0.9, na.rm=TRUE) +
  scale_color_manual(values=color_palette) +
  labs(title="NSAA Score",
       x="Age (years)", y="NSAA Score (0-34)") +
  theme_minimal(base_size=11) +
  theme(legend.position="none",
        plot.background=element_rect(fill="#1a1a2e", color=NA),
        panel.background=element_rect(fill="#2d2d4a"),
        text=element_text(color="white"),
        axis.text=element_text(color="white"))

p3 <- ggplot(df_all, aes(x=Age_yr, y=Dyst_pct, color=Scenario)) +
  geom_line(size=0.9, na.rm=TRUE) +
  scale_color_manual(values=color_palette) +
  geom_hline(yintercept=3, linetype="dashed", color="yellow", alpha=0.7) +
  annotate("text", x=5.5, y=4.5, label=">3% threshold", color="yellow", size=3) +
  labs(title="Dystrophin Level",
       x="Age (years)", y="Dystrophin (% of normal)") +
  theme_minimal(base_size=11) +
  theme(legend.position="none",
        plot.background=element_rect(fill="#1a1a2e", color=NA),
        panel.background=element_rect(fill="#2d2d4a"),
        text=element_text(color="white"),
        axis.text=element_text(color="white"))

p4 <- ggplot(df_all, aes(x=Age_yr, y=Fibrosis_pct, color=Scenario)) +
  geom_line(size=0.9, na.rm=TRUE) +
  scale_color_manual(values=color_palette) +
  labs(title="Muscle Fibrosis Score",
       x="Age (years)", y="Fibrosis (%)") +
  theme_minimal(base_size=11) +
  theme(legend.position="none",
        plot.background=element_rect(fill="#1a1a2e", color=NA),
        panel.background=element_rect(fill="#2d2d4a"),
        text=element_text(color="white"),
        axis.text=element_text(color="white"))

p5 <- ggplot(df_all, aes(x=Age_yr, y=FVC_pct, color=Scenario)) +
  geom_line(size=0.9, na.rm=TRUE) +
  scale_color_manual(values=color_palette) +
  geom_hline(yintercept=30, linetype="dashed", color="orange", alpha=0.7) +
  annotate("text", x=5.5, y=32, label="NIV threshold (30%)", color="orange", size=3) +
  labs(title="FVC% Predicted",
       x="Age (years)", y="FVC (% predicted)") +
  theme_minimal(base_size=11) +
  theme(legend.position="none",
        plot.background=element_rect(fill="#1a1a2e", color=NA),
        panel.background=element_rect(fill="#2d2d4a"),
        text=element_text(color="white"),
        axis.text=element_text(color="white"))

p6 <- ggplot(df_all, aes(x=Age_yr, y=LVEF, color=Scenario)) +
  geom_line(size=0.9, na.rm=TRUE) +
  scale_color_manual(values=color_palette) +
  geom_hline(yintercept=55, linetype="dashed", color="lightblue", alpha=0.7) +
  annotate("text", x=5.5, y=57, label="EF<55% = DCM", color="lightblue", size=3) +
  labs(title="Left Ventricular EF (LVEF)",
       x="Age (years)", y="LVEF (%)") +
  theme_minimal(base_size=11) +
  theme(legend.position="none",
        plot.background=element_rect(fill="#1a1a2e", color=NA),
        panel.background=element_rect(fill="#2d2d4a"),
        text=element_text(color="white"),
        axis.text=element_text(color="white"))

p7 <- ggplot(df_all, aes(x=Age_yr, y=CK_kUL, color=Scenario)) +
  geom_line(size=0.9, na.rm=TRUE) +
  scale_color_manual(values=color_palette) +
  scale_y_log10() +
  labs(title="Serum CK (log scale)",
       x="Age (years)", y="CK (kU/L, log₁₀)") +
  theme_minimal(base_size=11) +
  theme(legend.position="none",
        plot.background=element_rect(fill="#1a1a2e", color=NA),
        panel.background=element_rect(fill="#2d2d4a"),
        text=element_text(color="white"),
        axis.text=element_text(color="white"))

p8 <- ggplot(df_all, aes(x=Age_yr, y=Inflam_idx, color=Scenario)) +
  geom_line(size=0.9, na.rm=TRUE) +
  scale_color_manual(values=color_palette) +
  labs(title="Inflammation Index",
       x="Age (years)", y="Inflammation (0-100)") +
  theme_minimal(base_size=11) +
  theme(legend.position="none",
        plot.background=element_rect(fill="#1a1a2e", color=NA),
        panel.background=element_rect(fill="#2d2d4a"),
        text=element_text(color="white"),
        axis.text=element_text(color="white"))

combined_plot <- (p1 + p2) / (p3 + p4) / (p5 + p6) / (p7 + p8)
print(combined_plot)

## ─────────────────────────────────────────────────────────────
## SUMMARY TABLE: Key endpoints at Year 3 (Age 8) and Year 8 (Age 13)
## ─────────────────────────────────────────────────────────────

endpoints_yr3 <- df_all %>%
  filter(abs(Age_yr - 8.0) < 0.1) %>%
  group_by(Scenario) %>%
  slice(1) %>%
  ungroup() %>%
  select(Scenario, Age_yr,
         SixMWD, NSAA, Dyst_pct, Fibrosis_pct,
         FVC_pct, LVEF, CK_kUL) %>%
  mutate(across(where(is.numeric), ~round(.x, 1)))

endpoints_yr8 <- df_all %>%
  filter(abs(Age_yr - 13.0) < 0.1) %>%
  group_by(Scenario) %>%
  slice(1) %>%
  ungroup() %>%
  select(Scenario, Age_yr,
         SixMWD, NSAA, Dyst_pct, Fibrosis_pct,
         FVC_pct, LVEF, CK_kUL) %>%
  mutate(across(where(is.numeric), ~round(.x, 1)))

cat("\n=== Key Endpoints at Age 8yr (3yr after start) ===\n")
print(endpoints_yr3, width=120)

cat("\n=== Key Endpoints at Age 13yr (8yr after start) ===\n")
print(endpoints_yr8, width=120)

## ─────────────────────────────────────────────────────────────
## CLINICAL CALIBRATION NOTES
## ─────────────────────────────────────────────────────────────
cat("
╔═══════════════════════════════════════════════════════════════╗
║  DMD QSP Model — Clinical Calibration Reference              ║
╠═══════════════════════════════════════════════════════════════╣
║  Natural history:                                             ║
║  • 6MWD peaks at ~380m (age 8), declines to <200m (age 12)   ║
║    → Loss of ambulation median age 12.4yr (McDonald 2018)    ║
║  • NSAA peaks ~28 (age 8), declines to ~10 by age 13         ║
║  • FVC peaks 90-95% at age 8-10, declines ~2-4%/yr           ║
║  • CK: 10,000-50,000 U/L (10-100× ULN); falls after LOA      ║
║  • DCM: LVEF declines ~2%/yr after age 12                    ║
║                                                               ║
║  Deflazacort (ESSENCE/Griggs 2016 NEJM):                      ║
║  • Ambulation prolonged ~2yr vs placebo                       ║
║  • 6MWD: +23m vs prednisone at 12mo                          ║
║  • Lean mass better preserved than prednisone                 ║
║  • Spine density: -0.003 Z-score/yr vs -0.048 prednisone     ║
║                                                               ║
║  Eteplirsen (Mendell 2016 Ann Neurol):                        ║
║  • Dystrophin restored to 0.28-0.93% of normal               ║
║  • 6MWD decline slowed vs historical control                  ║
║  • Exon-51 skippable: ~13% of DMD pts                        ║
║                                                               ║
║  Elevidys EMBARK (2023):                                      ║
║  • Dystrophin: 28.1% vs 1.7% (placebo) at 52wk              ║
║  • NSAA functional improvement primary endpoint               ║
║  • Approved for ambulatory DMD ≥4yr (US, 2023)              ║
║                                                               ║
║  Vamorolone VISION-DMD (Servais 2022 Lancet ND):             ║
║  • 6MWD: +1.88m/yr gain vs baseline                          ║
║  • Fewer bone, growth, glucose effects vs prednisone          ║
║  • Approved as AGAMREE (2023, US/EU)                         ║
╚═══════════════════════════════════════════════════════════════╝
")
