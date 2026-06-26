## ============================================================
##  FSGS QSP mrgsolve ODE Model
##  Focal Segmental Glomerulosclerosis — 22-Compartment PK/PD
##  Author  : Claude Code Routine (CCR)
##  Date    : 2026-06-18
##  Version : 1.0
## ============================================================
##
##  DISEASE BIOLOGY OVERVIEW
##  ─────────────────────────────────────────────────────────
##  FSGS results from podocyte injury across three etiologic
##  axes that converge on common end-organ damage:
##
##  1. PRIMARY (circulating permeability factor axis)
##     CLCF1 / suPAR / anti-CD40 → αvβ3 integrin activation
##     → actin cytoskeleton disruption → foot process effacement
##     → massive proteinuria
##
##  2. SECONDARY (haemodynamic / nephron-loss axis)
##     RAAS activation → Ang II → efferent arteriolar constriction
##     → intraglomerular hypertension → podocyte mechanical shear
##     → progressive loss; drugs / obesity / hyperfiltration
##
##  3. GENETIC (structural / scaffold axis)
##     NPHS1/NPHS2 (nephrin/podocin), APOL1 (pore-forming),
##     INF2/TRPC6 (actin/calcium) mutations weaken the slit
##     diaphragm; modelled as reduced baseline POD and lower
##     repair capacity
##
##  FINAL COMMON PATHWAY
##     Podocyte loss → glomerular barrier failure → sclerosis
##     TGF-β / SMAD2/3 → irreversible ECM deposition (SCAR)
##     mTOR overactivation → impaired autophagy → apoptosis
##     Complement (C5b-9 MAC) → sublytic podocyte injury
##     B-cell axis → circulating permeability factor production
##
## ============================================================
##
##  CLINICAL CALIBRATION TARGETS (literature)
##  ─────────────────────────────────────────────────────────
##  Steroid monotherapy:
##    ~25% CR, ~30% PR at 6 months (Gipson et al. 2011 NEJM)
##  Cyclosporine / Tacrolimus:
##    ~70% PR at 6 months (Cattran et al. 1999 Kidney Int)
##  Rituximab (refractory):
##    ~50–60% response (Kronbichler 2014 Am J Nephrol)
##  Sparsentan DUPLEX trial (FDA approved 2023):
##    42% CR/PR vs 26% irbesartan at 1 year (Heerspink 2023 NEJM)
##  PLANET-I: UPCR –46% sparsentan vs –19% irbesartan @ 36 wk
##  eGFR decline untreated primary FSGS: 8–15 mL/min/year
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ─────────────────────────────────────────────────────────
##  mrgsolve model code
## ─────────────────────────────────────────────────────────

fsgs_code <- '
$PROB
FSGS QSP Model v1.0 — 22-compartment PK/PD
Focal Segmental Glomerulosclerosis
Includes: primary/secondary/genetic FSGS pathophysiology,
prednisolone + tacrolimus + rituximab + sparsentan PK/PD

$PARAM @annotated
// ── DISEASE SUBTYPE FLAGS ────────────────────────────────
FSGS_type : 1   : 1=primary, 2=secondary, 3=genetic

// ── BASELINE DISEASE STATE ───────────────────────────────
CLCF_base : 1.0 : Circulating permeability factor baseline (normalized)
POD_base  : 1.0 : Podocyte fraction baseline
FPE_base  : 0.02: Foot process effacement baseline (healthy)
PROT_base : 0.15: Proteinuria baseline g/day (healthy)
GFR_base  : 100 : eGFR baseline mL/min/1.73m2
SCAR_base : 0.01: Glomerular sclerosis fraction baseline
TGFb_base : 1.0 : TGF-beta signaling index baseline
COMP_base : 1.0 : Complement activity baseline
INFLAM_base: 1.0: Glomerular inflammation baseline
RAAS_base : 1.0 : RAAS activity index baseline
BCELL_base: 1.0 : B-cell fraction baseline

// ── DISEASE KINETIC PARAMETERS ───────────────────────────
k_CLCF_prod  : 0.06 : CLCF production rate constant /h
k_CLCF_elim  : 0.06 : CLCF elimination rate constant /h (t1/2~11.5h)
k_FPE_on     : 0.15 : FPE induction rate by CLCF (/h per unit CLCF)
k_FPE_repair : 0.03 : FPE repair rate when CLCF normalised /h
k_POD_loss   : 0.005: Podocyte loss rate (/h) driven by FPE + shear
k_POD_replen : 0.001: Podocyte replenishment rate /h (limited)
k_PROT_on    : 4.0  : Proteinuria drive rate g/day/h
k_PROT_off   : 0.08 : Proteinuria resolution rate /h
k_GFR_loss   : 0.0005: GFR loss rate driven by sclerosis /h
k_SCAR_form  : 0.004 : Sclerosis formation rate from TGFb + POD loss /h
k_SCAR_max   : 0.95 : Maximum sclerosis fraction (near-irreversible)
k_TGFb_prod  : 0.05 : TGF-beta production rate /h
k_TGFb_elim  : 0.05 : TGF-beta elimination /h
k_COMP_act   : 0.04 : Complement activation by inflammation /h
k_COMP_elim  : 0.04 : Complement turnover /h
k_INFLAM_prod: 0.06 : Inflammation production /h
k_INFLAM_elim: 0.06 : Inflammation resolution /h
k_RAAS_act   : 0.03 : RAAS activation by hyperfiltration /h
k_RAAS_elim  : 0.03 : RAAS turnover /h
k_BCELL_replen:0.004 : B-cell replenishment rate /h
k_POD_apop   : 0.008: mTOR-driven podocyte apoptosis rate /h

// ── PRIMARY FSGS — elevated circulating factor ───────────
CLCF_primary : 2.5  : CLCF set-point multiplier for primary FSGS
RAAS_secondary: 1.5 : RAAS set-point for secondary FSGS
POD_genetic   : 0.75 : Reduced podocyte baseline for genetic FSGS

// ── INTERACTION COEFFICIENTS ─────────────────────────────
alpha_CLCF_FPE : 0.6  : Sensitivity of FPE to CLCF elevation
alpha_RAAS_shear: 0.3 : RAAS-driven mechanical shear on POD
alpha_COMP_POD : 0.2  : Complement-mediated sublytic POD injury
alpha_TGFb_SCAR: 0.5  : TGFb contribution to sclerosis formation
alpha_FPE_POD  : 0.4  : FPE-driven podocyte detachment
alpha_INFLAM_TGFb: 0.35: Inflammation → TGFb cross-talk
alpha_mTOR_apop: 0.25 : mTOR overactivation → autophagy impairment

// ── PREDNISOLONE PK (2-CMT oral) ────────────────────────
CL_PRED  : 12    : Prednisolone clearance L/h
V1_PRED  : 30    : Prednisolone central volume L
V2_PRED  : 45    : Prednisolone peripheral volume L
Q_PRED   : 8     : Prednisolone intercompartmental CL L/h
ka_PRED  : 2.0   : Prednisolone absorption rate /h
F_PRED   : 0.82  : Prednisolone oral bioavailability
EC50_GR  : 0.150 : Prednisolone GR EC50 ug/mL (150 ng/mL)
hill_GR  : 1.5   : GR Hill coefficient
Imax_GR  : 0.80  : Maximum immunosuppression by prednisolone

// ── TACROLIMUS PK (2-CMT oral, whole blood) ──────────────
CL_TAC   : 2.5   : Tacrolimus clearance L/h
V1_TAC   : 85    : Tacrolimus central volume L
V2_TAC   : 150   : Tacrolimus peripheral volume L
Q_TAC    : 15    : Tacrolimus intercompartmental CL L/h
ka_TAC   : 0.4   : Tacrolimus absorption rate /h
EC50_CaN : 8.0   : Calcineurin inhibition EC50 ng/mL
hill_CaN : 2.0   : Calcineurin Hill coefficient
Imax_CaN : 0.85  : Maximum calcineurin inhibition

// ── RITUXIMAB PK (IV 2-CMT) ──────────────────────────────
CL_RTX   : 0.014 : Rituximab clearance L/h
V1_RTX   : 3.1   : Rituximab central volume L
V2_RTX   : 1.7   : Rituximab peripheral volume L
Q_RTX    : 0.012 : Rituximab intercompartmental CL L/h
EC50_Bkill:0.025 : RTX saturable B-cell kill EC50 mg/L (~25 ug/mL)
kmax_Bkill: 0.30 : Maximum B-cell kill rate /h by RTX

// ── SPARSENTAN PK (oral 2-CMT, dual AT1R+ETR antagonist) ─
CL_SPARS : 8.5   : Sparsentan clearance L/h
V1_SPARS : 60    : Sparsentan central volume L
V2_SPARS : 80    : Sparsentan peripheral volume L
Q_SPARS  : 5.0   : Sparsentan intercompartmental CL L/h
ka_SPARS : 1.2   : Sparsentan absorption rate /h
EC50_dual: 0.200 : Sparsentan dual AT1R/ETR EC50 ug/mL (200 ng/mL)
hill_dual: 1.2   : Sparsentan Hill coefficient
Imax_RAAS: 0.70  : Maximum RAAS/ET axis blockade by sparsentan

// ── BODY-SIZE / DOSE SCALING ──────────────────────────────
BWT      : 70    : Body weight kg (for mg/kg dosing)
BSA      : 1.73  : Body surface area m2 (for RTX dosing mg/m2)

$CMT @annotated
// ── DISEASE STATE (11 CMTs) ──────────────────────────────
CLCF    : Circulating permeability factor (normalized units)
POD     : Podocyte fraction (0-1)
FPE     : Foot process effacement index (0-1)
PROT    : Proteinuria (g/day)
GFR_c   : eGFR (mL/min/1.73m2)
SCAR    : Glomerular sclerosis fraction (0-1)
TGFb    : TGF-beta signaling index (normalized)
COMP    : Complement activity (normalized)
INFLAM  : Glomerular inflammation (normalized)
RAAS    : RAAS activity index (normalized)
BCELL   : B-cell fraction (0-1)

// ── PREDNISOLONE PK (3 CMTs) ─────────────────────────────
PRED_DEP : Prednisolone depot/GI absorption (mg)
PRED1    : Prednisolone central compartment (mg)
PRED2    : Prednisolone peripheral compartment (mg)

// ── TACROLIMUS PK (3 CMTs) ───────────────────────────────
TAC_DEP  : Tacrolimus depot/GI absorption (ug)
TAC1     : Tacrolimus central compartment (ug)
TAC2     : Tacrolimus peripheral compartment (ug)

// ── RITUXIMAB PK (2 CMTs) ────────────────────────────────
RTX1     : Rituximab central compartment (mg)
RTX2     : Rituximab peripheral compartment (mg)

// ── SPARSENTAN PK (3 CMTs) ───────────────────────────────
SPARS_DEP: Sparsentan depot/GI absorption (mg)
SPARS1   : Sparsentan central compartment (mg)
SPARS2   : Sparsentan peripheral compartment (mg)

$INIT
// ── DISEASE STATE INIT (subtype-dependent, set in $MAIN) ─
CLCF    = 1.0
POD     = 1.0
FPE     = 0.02
PROT    = 0.15
GFR_c   = 100.0
SCAR    = 0.01
TGFb    = 1.0
COMP    = 1.0
INFLAM  = 1.0
RAAS    = 1.0
BCELL   = 1.0

// ── PK INIT ───────────────────────────────────────────────
PRED_DEP = 0.0
PRED1    = 0.0
PRED2    = 0.0
TAC_DEP  = 0.0
TAC1     = 0.0
TAC2     = 0.0
RTX1     = 0.0
RTX2     = 0.0
SPARS_DEP= 0.0
SPARS1   = 0.0
SPARS2   = 0.0

$MAIN
// ─────────────────────────────────────────────────────────
//  FSGS SUBTYPE INITIALISATION
//  Called at t=0 (NEWIND>0) to set disease-specific ICs
// ─────────────────────────────────────────────────────────
if (NEWIND <= 1) {

  // Primary FSGS: high circulating permeability factor
  if (FSGS_type == 1) {
    CLCF_0   = CLCF_primary;  // ~2.5 at presentation
    POD_0    = 1.0;
    FPE_0    = 0.35;           // moderate FPE already present
    PROT_0   = 6.0;            // nephrotic-range proteinuria
    GFR_0    = 80.0;           // mild reduction
    RAAS_0   = 1.1;
    BCELL_0  = 1.2;            // B-cell driven permeability
  }

  // Secondary FSGS: haemodynamic + RAAS axis dominant
  else if (FSGS_type == 2) {
    CLCF_0  = 1.2;
    POD_0   = 0.85;            // early podocyte loss
    FPE_0   = 0.20;
    PROT_0  = 3.0;             // sub-nephrotic to nephrotic
    GFR_0   = 65.0;            // reduced (e.g. solitary kidney)
    RAAS_0  = RAAS_secondary;  // 1.5 — elevated
    BCELL_0 = 1.0;
  }

  // Genetic FSGS: structurally impaired podocytes
  else {
    CLCF_0  = 1.1;
    POD_0   = POD_genetic;     // 0.75 — genetically reduced
    FPE_0   = 0.40;
    PROT_0  = 7.0;
    GFR_0   = 70.0;
    RAAS_0  = 1.2;
    BCELL_0 = 1.0;
  }

  CLCF_0    = CLCF_0;
  POD_0     = POD_0;
  FPE_0     = FPE_0;
  PROT_0    = PROT_0;
  GFR_0     = GFR_0;
  RAAS_0    = RAAS_0;
  BCELL_0   = BCELL_0;
}

// ─────────────────────────────────────────────────────────
//  DERIVED PK CONCENTRATIONS
// ─────────────────────────────────────────────────────────

// Prednisolone: convert mg → ug/mL using V1 (L)
double Cp_PRED = PRED1 / V1_PRED;   // ug/mL (mg/L)

// Tacrolimus: whole blood ng/mL = ug / V1 (L) * 1000
double Cp_TAC  = (TAC1 / V1_TAC) * 1000.0;  // ng/mL

// Rituximab: mg/L
double Cp_RTX  = RTX1 / V1_RTX;     // mg/L (= ug/mL)

// Sparsentan: ug/mL
double Cp_SPARS = SPARS1 / V1_SPARS; // mg/L ~ ug/mL

// ─────────────────────────────────────────────────────────
//  PHARMACODYNAMIC EFFECT FUNCTIONS  (0 = no drug, 1 = max)
// ─────────────────────────────────────────────────────────

// Prednisolone: glucocorticoid receptor occupancy
double EFF_GR = Imax_GR * pow(Cp_PRED, hill_GR) /
                (pow(EC50_GR, hill_GR) + pow(Cp_PRED, hill_GR));

// Tacrolimus: calcineurin inhibition
double EFF_CaN = Imax_CaN * pow(Cp_TAC, hill_CaN) /
                 (pow(EC50_CaN, hill_CaN) + pow(Cp_TAC, hill_CaN));

// Sparsentan: dual AT1R + ETR blockade → reduces RAAS axis
double EFF_dual = Imax_RAAS * pow(Cp_SPARS, hill_dual) /
                  (pow(EC50_dual, hill_dual) + pow(Cp_SPARS, hill_dual));

// Combined immunosuppression (prednisolone + tacrolimus, not additive beyond 1)
double IS_combined = fmin(1.0, EFF_GR + EFF_CaN - EFF_GR * EFF_CaN);

// ─────────────────────────────────────────────────────────
//  GUARD CURRENT STATE (protect against negatives)
// ─────────────────────────────────────────────────────────

double clcf   = fmax(0.01, CLCF);
double pod    = fmax(0.001, fmin(1.0, POD));
double fpe    = fmax(0.0,   fmin(1.0, FPE));
double prot   = fmax(0.0,   PROT);
double gfr    = fmax(1.0,   GFR_c);
double scar   = fmax(0.0,   fmin(k_SCAR_max, SCAR));
double tgfb   = fmax(0.01,  TGFb);
double comp   = fmax(0.01,  COMP);
double inflam = fmax(0.01,  INFLAM);
double raas   = fmax(0.01,  RAAS);
double bcell  = fmax(0.0,   fmin(1.0, BCELL));

$ODE
// =========================================================
//  DISEASE STATE ODES  (compartments 1–11)
// =========================================================

// ── 1. CLCF — Circulating permeability factor ─────────────
//    Production stimulated by B-cells; suppressed by IS
//    Subtype set-point drives chronic elevation in primary FSGS
double CLCF_setpt = (FSGS_type == 1) ? CLCF_primary : CLCF_base;
dxdt_CLCF = k_CLCF_prod * bcell * (1.0 - IS_combined) * CLCF_setpt
           - k_CLCF_elim * clcf;

// ── 2. POD — Podocyte fraction ─────────────────────────────
//    Loss: FPE-mediated detachment + RAAS mechanical shear
//          + complement sublytic injury + mTOR apoptosis
//    Gain: limited regeneration from parietal epithelial cells
//    EFF_GR suppresses loss (steroids are cytoprotective)
double POD_loss = k_POD_loss   * fpe  * alpha_FPE_POD
                + k_POD_loss   * (raas - 1.0) * alpha_RAAS_shear
                + k_POD_loss   * (comp - 1.0) * alpha_COMP_POD
                + k_POD_apop   * alpha_mTOR_apop;
double POD_gain = k_POD_replen * (1.0 - pod) * (1.0 + 0.5 * EFF_GR);
dxdt_POD = POD_gain - POD_loss * pod * (1.0 - 0.4 * EFF_GR);

// ── 3. FPE — Foot process effacement ──────────────────────
//    Driven by CLCF via αvβ3 integrin signalling
//    Suppressed by IS_combined (steroid + calcineurin inhibitor)
//    Repair limited; partial reversal with treatment
double FPE_drive = k_FPE_on * (clcf - 1.0) * alpha_CLCF_FPE * (1.0 - scar);
double FPE_repair= k_FPE_repair * fpe * (1.0 + IS_combined);
dxdt_FPE = FPE_drive - FPE_repair;

// ── 4. PROT — Proteinuria (g/day) ─────────────────────────
//    Driven by FPE and reduced podocyte coverage (1-pod)
//    RAAS haemodynamic component (sparsentan acts here)
//    Resolved by treatment; floor at 0.15 healthy
double PROT_drive  = k_PROT_on * fpe * (2.0 - pod)
                   + 2.0 * (raas - 1.0) * (1.0 - pod);
double PROT_resolve= k_PROT_off * prot * (1.0 + IS_combined + EFF_dual);
dxdt_PROT = PROT_drive - PROT_resolve;

// ── 5. GFR_c — eGFR (mL/min/1.73m²) ──────────────────────
//    Decline driven by sclerosis (irreversible) + acute proteinuric
//    haemodynamic injury; RAAS blockade is renoprotective
double GFR_decline = k_GFR_loss * gfr
                   * (scar + 0.3 * (prot / 10.0) + 0.2 * (raas - 1.0));
double GFR_protect = k_GFR_loss * gfr * (EFF_dual * 0.6 + EFF_GR * 0.2);
dxdt_GFR_c = -(GFR_decline - GFR_protect);

// ── 6. SCAR — Glomerular sclerosis (irreversible) ─────────
//    Driven by TGFb × (1-pod); near-irreversible once formed
//    mTOR (alpha_mTOR_apop) accelerates fibrosis via autophagy block
//    IS_combined modestly slows (not reverses) progression
double SCAR_rate = k_SCAR_form
                 * tgfb * alpha_TGFb_SCAR
                 * (1.0 - pod)
                 * (k_SCAR_max - scar)
                 * (1.0 - 0.25 * IS_combined);
dxdt_SCAR = fmax(0.0, SCAR_rate);  // unidirectional — sclerosis is permanent

// ── 7. TGFb — TGF-β signalling index ──────────────────────
//    Produced by inflamed mesangium and podocytes under mechanical stress
//    Suppressed by IS; fuels ECM fibrosis (SCAR)
dxdt_TGFb = k_TGFb_prod * (inflam * alpha_INFLAM_TGFb + (raas - 1.0) * 0.3)
           - k_TGFb_elim * tgfb * (1.0 + 0.3 * IS_combined);

// ── 8. COMP — Complement activity ─────────────────────────
//    C3/C5 activation driven by immune complex deposition (inflam)
//    C5b-9 MAC causes sublytic podocyte injury
//    Steroids suppress complement activation modestly
dxdt_COMP = k_COMP_act * inflam
           - k_COMP_elim * comp * (1.0 + 0.2 * EFF_GR);

// ── 9. INFLAM — Glomerular inflammation ───────────────────
//    Macrophage/lymphocyte infiltration driven by CLCF + proteinuria
//    Primary therapeutic target for steroids and calcineurin inhibitors
double INFLAM_drive = k_INFLAM_prod * clcf * (1.0 + 0.5 * (prot / 8.0));
double INFLAM_suppress = k_INFLAM_elim * inflam * (1.0 + IS_combined * 2.0);
dxdt_INFLAM = INFLAM_drive - INFLAM_suppress;

// ── 10. RAAS — RAAS activity index ────────────────────────
//    Secondary FSGS set-point = 1.5; feedback from proteinuria
//    and reduced GFR (macula densa signal)
//    Dual blockade by sparsentan (EFF_dual) — main renoprotective axis
double RAAS_setpt = (FSGS_type == 2) ? RAAS_secondary : RAAS_base;
double RAAS_stim  = k_RAAS_act * (1.0 + 0.3 * (prot / 5.0)
                  + 0.2 * (1.0 - gfr / GFR_base));
dxdt_RAAS = RAAS_stim * RAAS_setpt
           - k_RAAS_elim * raas * (1.0 + EFF_dual * 2.5);

// ── 11. BCELL — B-cell fraction ───────────────────────────
//    Drives production of circulating permeability factor (primary FSGS)
//    Depleted by rituximab via saturable CD20 killing
//    Replenishment from bone marrow (slow, ~6 months for full recovery)
double RTX_kill = kmax_Bkill * Cp_RTX / (EC50_Bkill + Cp_RTX);
dxdt_BCELL = k_BCELL_replen * (1.0 - bcell) - RTX_kill * bcell;

// =========================================================
//  PREDNISOLONE PK ODES  (compartments 12–14)
// =========================================================
//  Dose administered to PRED_DEP depot; F_PRED fraction absorbed
//  2-CMT with central (PRED1) and peripheral (PRED2) volumes
// ─────────────────────────────────────────────────────────
dxdt_PRED_DEP = -ka_PRED * PRED_DEP;
dxdt_PRED1    = ka_PRED * F_PRED * PRED_DEP
               + (Q_PRED / V2_PRED) * PRED2
               - (CL_PRED / V1_PRED + Q_PRED / V1_PRED) * PRED1;
dxdt_PRED2    = (Q_PRED / V1_PRED) * PRED1
               - (Q_PRED / V2_PRED) * PRED2;

// =========================================================
//  TACROLIMUS PK ODES  (compartments 15–17)
// =========================================================
//  Whole-blood concentration drives calcineurin inhibition
//  High Vd reflects extensive erythrocyte/lymphocyte binding
// ─────────────────────────────────────────────────────────
dxdt_TAC_DEP = -ka_TAC * TAC_DEP;
dxdt_TAC1    = ka_TAC * TAC_DEP
              + (Q_TAC / V2_TAC) * TAC2
              - (CL_TAC / V1_TAC + Q_TAC / V1_TAC) * TAC1;
dxdt_TAC2    = (Q_TAC / V1_TAC) * TAC1
              - (Q_TAC / V2_TAC) * TAC2;

// =========================================================
//  RITUXIMAB PK ODES  (compartments 18–19)
// =========================================================
//  IV 2-CMT; slow clearance (~50-day half-life central)
//  Additional target-mediated elimination via B-cell binding
//  Simplified: CL_RTX absorbs TMDD into apparent CL
// ─────────────────────────────────────────────────────────
dxdt_RTX1    = (Q_RTX / V2_RTX) * RTX2
              - (CL_RTX / V1_RTX + Q_RTX / V1_RTX) * RTX1;
dxdt_RTX2    = (Q_RTX / V1_RTX) * RTX1
              - (Q_RTX / V2_RTX) * RTX2;

// =========================================================
//  SPARSENTAN PK ODES  (compartments 20–22)
// =========================================================
//  Oral 2-CMT; dual AT1R + ETR antagonist (FDA approved 2023)
//  DUPLEX trial: 800 mg QD; rapid absorption ka=1.2/h
// ─────────────────────────────────────────────────────────
dxdt_SPARS_DEP = -ka_SPARS * SPARS_DEP;
dxdt_SPARS1    = ka_SPARS * SPARS_DEP
                + (Q_SPARS / V2_SPARS) * SPARS2
                - (CL_SPARS / V1_SPARS + Q_SPARS / V1_SPARS) * SPARS1;
dxdt_SPARS2    = (Q_SPARS / V1_SPARS) * SPARS1
                - (Q_SPARS / V2_SPARS) * SPARS2;

$TABLE
// ─────────────────────────────────────────────────────────
//  DERIVED CLINICAL ENDPOINTS
// ─────────────────────────────────────────────────────────

// Prednisolone concentration (ug/mL)
double PRED_Cp  = PRED1 / V1_PRED;

// Tacrolimus trough (ng/mL whole blood)
double TAC_trough = (TAC1 / V1_TAC) * 1000.0;

// Rituximab serum (mg/L)
double RTX_serum  = RTX1 / V1_RTX;

// Sparsentan plasma (ug/mL)
double SPARS_Cp   = SPARS1 / V1_SPARS;

// Remission classification
// Complete  : PROT < 0.3 g/day  AND GFR_c stable (< 25% decline from GFR_base)
// Partial   : 50% reduction in PROT from baseline (>= 3.5 if initial nephrotic)
// No response: <50% PROT reduction

// UPCR proxy (g/g; approximating g/day as equivalent for 1.73m2 BSA standard)
double UPCR      = PROT / 0.15;   // normalized ratio relative to healthy

// Nephrotic syndrome flag: PROT > 3.5 g/day
double NS_flag   = (PROT > 3.5) ? 1.0 : 0.0;

// Chronic kidney disease stage (KDIGO)
double CKD_stage = (GFR_c >= 90) ? 1.0 :
                   (GFR_c >= 60) ? 2.0 :
                   (GFR_c >= 45) ? 3.0 :
                   (GFR_c >= 30) ? 4.0 : 5.0;

// Drug effect summaries
double DRUG_IS    = IS_combined;     // combined immunosuppression
double DRUG_RAAS  = EFF_dual;        // RAAS/ET blockade
double DRUG_GR    = EFF_GR;          // glucocorticoid receptor

$CAPTURE
PRED_Cp TAC_trough RTX_serum SPARS_Cp
EFF_GR EFF_CaN EFF_dual IS_combined
UPCR NS_flag CKD_stage
DRUG_IS DRUG_RAAS DRUG_GR
'

## ─────────────────────────────────────────────────────────
##  Compile model
## ─────────────────────────────────────────────────────────

mod <- mcode("fsgs_qsp", fsgs_code, quiet = TRUE)

## ─────────────────────────────────────────────────────────
##  HELPER: build a dosing event table
##    tstart  : start time (h)
##    tend    : end time (h)
##    amt     : dose amount (model units)
##    ii      : dosing interval (h)
##    cmt     : depot compartment number
##    addl    : additional doses (= total - 1)
## ─────────────────────────────────────────────────────────

## ─────────────────────────────────────────────────────────
##  SCENARIO DEFINITIONS
##  All simulations: primary FSGS (type=1), t = 0 → 8760 h (365 d)
##  PRED dose : 70 mg QD oral (1 mg/kg × 70 kg)
##  TAC  dose : 3 mg BID oral (typical induction for FSGS)
##  RTX  dose : 375 mg/m2 IV × 4 wk = 649 mg × 4 doses (BSA 1.73)
##  SPARS dose: 800 mg QD oral
## ─────────────────────────────────────────────────────────

SIM_END  <- 8760   # 365 days in hours
BWT_val  <- 70     # kg
BSA_val  <- 1.73   # m2

PRED_dose   <- BWT_val * 1.0   # 70 mg QD
TAC_dose    <- 3.0             # 3 mg BID
RTX_dose    <- round(375 * BSA_val, 0)  # ~649 mg per infusion
SPARS_dose  <- 800             # mg QD

# Tapering prednisolone: full dose 12 wk, then halve every 4 wk
PRED_taper_wk12 <- PRED_dose / 2
PRED_taper_wk16 <- PRED_dose / 4
PRED_taper_wk20 <- PRED_dose / 8   # maintenance or discontinue at wk24

## ─────────────────────────────────────────────────────────
##  SCENARIO 1: Natural History (no treatment)
## ─────────────────────────────────────────────────────────

ev_NatHist <- ev(time = 0, amt = 0, cmt = 1)   # dummy event

## ─────────────────────────────────────────────────────────
##  SCENARIO 2: PRED Monotherapy
##  1 mg/kg/d oral × 12 wk then taper to wk 24
## ─────────────────────────────────────────────────────────

ev_pred_full  <- ev(amt = PRED_dose,        cmt = "PRED_DEP",
                    ii  = 24, addl = 83,    time = 0)      # wk 0–12  (84 d)
ev_pred_half  <- ev(amt = PRED_taper_wk12,  cmt = "PRED_DEP",
                    ii  = 24, addl = 27,    time = 84*24)  # wk 12–16
ev_pred_qtr   <- ev(amt = PRED_taper_wk16,  cmt = "PRED_DEP",
                    ii  = 24, addl = 27,    time = 112*24) # wk 16–20
ev_pred_maint <- ev(amt = PRED_taper_wk20,  cmt = "PRED_DEP",
                    ii  = 24, addl = 27,    time = 140*24) # wk 20–24

ev_PRED_mono <- ev_pred_full + ev_pred_half + ev_pred_qtr + ev_pred_maint

## ─────────────────────────────────────────────────────────
##  SCENARIO 3: PRED + TAC (steroid-resistant protocol)
##  PRED same taper as above; TAC 3 mg BID throughout 52 wk
## ─────────────────────────────────────────────────────────

ev_tac <- ev(amt = TAC_dose, cmt = "TAC_DEP",
             ii  = 12, addl = (365*2 - 1), time = 0)   # BID × 52 wk

ev_PRED_TAC <- ev_PRED_mono + ev_tac

## ─────────────────────────────────────────────────────────
##  SCENARIO 4: PRED + TAC + RTX (refractory FSGS protocol)
##  RTX 375 mg/m2 IV × 4 weekly doses (wk 0,1,2,3)
##  Calibration: Kronbichler 2014 ~50-60% response
## ─────────────────────────────────────────────────────────

ev_rtx <- ev(amt = RTX_dose, cmt = "RTX1",
             ii  = 168, addl = 3, time = 0)   # 4 weekly IV infusions

ev_PRED_TAC_RTX <- ev_PRED_TAC + ev_rtx

## ─────────────────────────────────────────────────────────
##  SCENARIO 5: Sparsentan Monotherapy (DUPLEX trial protocol)
##  800 mg QD oral × 52 wk
##  Calibration: 42% CR/PR (Heerspink 2023 NEJM); UPCR –46% vs –19% irbesartan
## ─────────────────────────────────────────────────────────

ev_SPARS_mono <- ev(amt = SPARS_dose, cmt = "SPARS_DEP",
                    ii  = 24, addl = 364, time = 0)

## ─────────────────────────────────────────────────────────
##  SCENARIO 6: FULL Combination
##  PRED (taper) + TAC (BID) + SPARS (QD) — no RTX
##  Rationale: maximal immunosuppression + haemodynamic protection
## ─────────────────────────────────────────────────────────

ev_FULL <- ev_PRED_mono + ev_tac + ev_SPARS_mono

## ─────────────────────────────────────────────────────────
##  SIMULATION FUNCTION
## ─────────────────────────────────────────────────────────

run_scenario <- function(model, events, label,
                         fsgs_type = 1,
                         sim_end   = SIM_END,
                         delta     = 24) {

  out <- model %>%
    param(FSGS_type = fsgs_type) %>%
    mrgsim(events = events,
           end    = sim_end,
           delta  = delta,
           obsonly= TRUE,
           carry_out = "evid") %>%
    as.data.frame() %>%
    mutate(scenario = label,
           day      = time / 24)
  return(out)
}

## ─────────────────────────────────────────────────────────
##  RUN ALL 6 SCENARIOS
## ─────────────────────────────────────────────────────────

message("Running 6 FSGS treatment scenarios...")

results <- bind_rows(
  run_scenario(mod, ev_NatHist,       "1_NatHist"),
  run_scenario(mod, ev_PRED_mono,     "2_PRED_mono"),
  run_scenario(mod, ev_PRED_TAC,      "3_PRED_TAC"),
  run_scenario(mod, ev_PRED_TAC_RTX,  "4_PRED_TAC_RTX"),
  run_scenario(mod, ev_SPARS_mono,    "5_SPARS_mono"),
  run_scenario(mod, ev_FULL,          "6_FULL")
)

## Rename for clarity
results <- results %>%
  mutate(
    scenario_label = recode(scenario,
      "1_NatHist"      = "Natural History",
      "2_PRED_mono"    = "Prednisolone",
      "3_PRED_TAC"     = "Pred + Tacrolimus",
      "4_PRED_TAC_RTX" = "Pred + TAC + RTX",
      "5_SPARS_mono"   = "Sparsentan",
      "6_FULL"         = "Pred + TAC + Sparsentan"
    )
  )

## ─────────────────────────────────────────────────────────
##  SUMMARY TABLE AT KEY TIMEPOINTS
## ─────────────────────────────────────────────────────────

summary_tbl <- results %>%
  filter(day %in% c(0, 90, 180, 270, 365)) %>%
  group_by(scenario_label, day) %>%
  summarise(
    PROT_gday   = round(mean(PROT),   2),
    GFR         = round(mean(GFR_c),  1),
    POD_frac    = round(mean(POD),    3),
    SCAR_frac   = round(mean(SCAR),   3),
    BCELL_frac  = round(mean(BCELL),  3),
    CLCF_idx    = round(mean(CLCF),   2),
    .groups     = "drop"
  ) %>%
  arrange(day, scenario_label)

message("\n===== FSGS MODEL SUMMARY =====")
print(summary_tbl, n = 40)

## ─────────────────────────────────────────────────────────
##  REMISSION CLASSIFICATION AT 12 MONTHS
## ─────────────────────────────────────────────────────────

remission_d365 <- results %>%
  filter(day == 365) %>%
  mutate(
    baseline_prot = 6.0,   # primary FSGS initial PROT
    prot_change   = (baseline_prot - PROT) / baseline_prot * 100,
    remission = case_when(
      PROT < 0.3                      ~ "Complete",
      prot_change >= 50               ~ "Partial",
      TRUE                            ~ "No Response"
    )
  ) %>%
  select(scenario_label, PROT, GFR_c, prot_change, remission)

message("\n===== REMISSION STATUS AT 12 MONTHS =====")
print(remission_d365)

## ─────────────────────────────────────────────────────────
##  PLOTTING
## ─────────────────────────────────────────────────────────

scenario_colors <- c(
  "Natural History"        = "#E74C3C",
  "Prednisolone"           = "#F39C12",
  "Pred + Tacrolimus"      = "#3498DB",
  "Pred + TAC + RTX"       = "#9B59B6",
  "Sparsentan"             = "#27AE60",
  "Pred + TAC + Sparsentan"= "#2C3E50"
)

## Panel A: Proteinuria over time
p_prot <- ggplot(results, aes(day, PROT,
                              color = scenario_label,
                              linetype = scenario_label)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 3.5, linetype = "dashed", color = "gray50", alpha = 0.7) +
  geom_hline(yintercept = 0.3, linetype = "dotted", color = "steelblue", alpha = 0.7) +
  annotate("text", x = 370, y = 3.7, label = "Nephrotic (3.5 g/d)",
           size = 2.8, hjust = 1, color = "gray50") +
  annotate("text", x = 370, y = 0.1, label = "CR (<0.3 g/d)",
           size = 2.8, hjust = 1, color = "steelblue") +
  scale_color_manual(values = scenario_colors) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(title    = "A. Proteinuria (g/day)",
       x        = "Day",
       y        = "Proteinuria (g/day)",
       color    = "Scenario",
       linetype = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "right",
        legend.text     = element_text(size = 8))

## Panel B: eGFR over time
p_gfr <- ggplot(results, aes(day, GFR_c,
                             color = scenario_label,
                             linetype = scenario_label)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = c(60, 30), linetype = "dashed",
             color = "gray60", alpha = 0.6) +
  scale_color_manual(values = scenario_colors) +
  labs(title    = "B. eGFR (mL/min/1.73m²)",
       x        = "Day",
       y        = "eGFR",
       color    = "Scenario",
       linetype = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

## Panel C: Podocyte fraction
p_pod <- ggplot(results, aes(day, POD,
                             color = scenario_label,
                             linetype = scenario_label)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scenario_colors) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(title    = "C. Podocyte Fraction",
       x        = "Day",
       y        = "Podocyte fraction (0–1)",
       color    = "Scenario",
       linetype = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

## Panel D: Glomerular Sclerosis
p_scar <- ggplot(results, aes(day, SCAR,
                              color = scenario_label,
                              linetype = scenario_label)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scenario_colors) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(title    = "D. Glomerular Sclerosis",
       x        = "Day",
       y        = "Sclerosis fraction (0–1)",
       color    = "Scenario",
       linetype = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

## Panel E: B-cell fraction (Rituximab effect)
p_bcell <- ggplot(results, aes(day, BCELL,
                               color = scenario_label,
                               linetype = scenario_label)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scenario_colors) +
  scale_y_continuous(limits = c(0, 1.2)) +
  labs(title    = "E. B-Cell Fraction (CD20+)",
       x        = "Day",
       y        = "B-cell fraction (0–1)",
       color    = "Scenario",
       linetype = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

## Panel F: CLCF circulating factor
p_clcf <- ggplot(results, aes(day, CLCF,
                              color = scenario_label,
                              linetype = scenario_label)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 1.0, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = scenario_colors) +
  labs(title    = "F. Circulating Permeability Factor (CLCF)",
       x        = "Day",
       y        = "CLCF index (normalized)",
       color    = "Scenario",
       linetype = "Scenario") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

## Combine panels (requires patchwork)
combined_plot <- (p_prot | p_gfr) /
                 (p_pod  | p_scar) /
                 (p_bcell| p_clcf) +
  plot_annotation(
    title    = "FSGS QSP Model — 6 Treatment Scenarios",
    subtitle = "Primary FSGS; 22-compartment mrgsolve ODE; 365-day simulation",
    caption  = paste0(
      "Calibration targets: CR ~25% steroid mono (Gipson 2011 NEJM); ",
      "TAC ~70% PR (Cattran 1999 Kidney Int);\n",
      "RTX ~55% response (Kronbichler 2014 Am J Nephrol); ",
      "Sparsentan 42% CR/PR (Heerspink 2023 NEJM DUPLEX)"
    )
  )

## ─────────────────────────────────────────────────────────
##  PK CONCENTRATION PROFILES (subset first 7 days for clarity)
## ─────────────────────────────────────────────────────────

pk_early <- results %>%
  filter(day <= 7) %>%
  select(day, scenario_label, PRED_Cp, TAC_trough, RTX_serum, SPARS_Cp)

p_pk_pred <- ggplot(
    filter(pk_early, scenario_label %in%
             c("Prednisolone", "Pred + Tacrolimus",
               "Pred + TAC + RTX", "Pred + TAC + Sparsentan")),
    aes(day, PRED_Cp, color = scenario_label)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Prednisolone Concentration",
       x = "Day", y = "Cp (ug/mL)", color = NULL) +
  theme_bw(base_size = 10)

p_pk_tac <- ggplot(
    filter(pk_early, scenario_label %in%
             c("Pred + Tacrolimus", "Pred + TAC + RTX",
               "Pred + TAC + Sparsentan")),
    aes(day, TAC_trough, color = scenario_label)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = c(5, 10), linetype = "dashed", color = "gray50") +
  annotate("text", x = 7, y = 5.3, label = "Target 5 ng/mL",
           size = 2.5, hjust = 1, color = "gray40") +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Tacrolimus Whole-Blood",
       x = "Day", y = "Conc (ng/mL)", color = NULL) +
  theme_bw(base_size = 10)

p_pk_rtx <- ggplot(
    filter(results, scenario_label == "Pred + TAC + RTX", day <= 30),
    aes(day, RTX_serum)) +
  geom_line(color = "#9B59B6", linewidth = 1) +
  labs(title = "Rituximab Serum (first 30 d)",
       x = "Day", y = "Conc (mg/L)") +
  theme_bw(base_size = 10)

p_pk_spars <- ggplot(
    filter(pk_early, scenario_label %in%
             c("Sparsentan", "Pred + TAC + Sparsentan")),
    aes(day, SPARS_Cp, color = scenario_label)) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = scenario_colors) +
  labs(title = "Sparsentan Plasma",
       x = "Day", y = "Cp (ug/mL)", color = NULL) +
  theme_bw(base_size = 10)

pk_plot <- (p_pk_pred | p_pk_tac) / (p_pk_rtx | p_pk_spars) +
  plot_annotation(title = "FSGS QSP — Drug PK Profiles (first 7 days)")

## ─────────────────────────────────────────────────────────
##  PRINT & SAVE
## ─────────────────────────────────────────────────────────

print(combined_plot)
print(pk_plot)

## Save to files (PDF for vector quality)
tryCatch({
  ggsave("fsgs_qsp_scenarios.pdf",
         plot = combined_plot,
         width = 14, height = 16, device = "pdf")
  ggsave("fsgs_qsp_pk.pdf",
         plot = pk_plot,
         width = 12, height = 8, device = "pdf")
  message("Plots saved: fsgs_qsp_scenarios.pdf, fsgs_qsp_pk.pdf")
}, error = function(e) {
  message("Plot save skipped (run interactively to save): ", e$message)
})

## ─────────────────────────────────────────────────────────
##  RETURN OBJECTS (for interactive use)
## ─────────────────────────────────────────────────────────

message("\n===== MODEL COMPILATION COMPLETE =====")
message("Objects available: mod, results, summary_tbl, remission_d365")
message("Plots: combined_plot, pk_plot")

invisible(list(
  model          = mod,
  results        = results,
  summary        = summary_tbl,
  remission_1yr  = remission_d365,
  plot_scenarios = combined_plot,
  plot_pk        = pk_plot
))
