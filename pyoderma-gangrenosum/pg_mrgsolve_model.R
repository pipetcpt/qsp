# =============================================================================
# Pyoderma Gangrenosum (PG) — mrgsolve QSP Model
# =============================================================================
# Compartments (≥ 15 ODEs):
#   Drug PK:        ADA_SC, ADA_CENT, ADA_PERI,
#                   IFX_CENT, IFX_PERI,
#                   ANA_SC, ANA_CENT,
#                   CsA_GUT, CsA_CENT,
#                   UST_SC, UST_CENT,
#                   PRED_GUT, PRED_CENT
#   Inflammation:   TNFa, IL1b, IL17A, IL6, IL23, IL8
#                   Neutrophil, NET, Th17, Treg, Macro_M1
#   Tissue:         Keratinocyte, MMP9_act, ROS_lesion,
#                   UlcerArea (cm²), HealedFraction
#   Endpoints:      PARACELSUS, Pain_VAS, DLQI
#   PD bridge:      CRP, Calprotectin
# -----------------------------------------------------------------------------
# Reference parameters drawn from:
#   - Adalimumab PIONEER-PG case series; Kimball AB 2016
#   - Infliximab Brooklyn trial (Brooklyn JE, Gut 2006 — PG/UC)
#   - Anakinra in PAPA syndrome (Brenner M 2009)
#   - Cyclosporine STOP-GAP (Ormerod AD 2015, BMJ)
#   - Ustekinumab Goldminz AM 2012 case series
#   - Prednisone STOP-GAP comparator (Ormerod AD 2015)
# =============================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

pg_qsp_code <- '
$PROB
# Pyoderma Gangrenosum QSP model
# v1.0 — Drug PK/PD, neutrophilic inflammation, ulcer dynamics, healing
# Time unit: days

$PARAM @annotated
// -------- Adalimumab (ADA) --------
KA_ADA       :  0.30   : ADA absorption rate (1/d)
F1_ADA       :  0.64   : ADA SC bioavailability
CL_ADA       :  0.31   : ADA clearance (L/d)
V1_ADA       :  4.7    : ADA central volume (L)
Q_ADA        :  0.45   : ADA inter-compartmental Q (L/d)
V2_ADA       :  2.6    : ADA peripheral volume (L)

// -------- Infliximab (IFX) --------
CL_IFX       :  0.32   : IFX clearance (L/d)
V1_IFX       :  3.5    : IFX central volume (L)
Q_IFX        :  0.40   : IFX Q (L/d)
V2_IFX       :  2.0    : IFX peripheral volume (L)

// -------- Anakinra (ANA) --------
KA_ANA       :  6.0    : ANA absorption rate (1/d)
F1_ANA       :  0.95   : ANA SC bioavailability
CL_ANA       :  17.0   : ANA clearance (L/d)
V_ANA        :  20.0   : ANA volume of distribution (L)

// -------- Cyclosporine (CsA) --------
KA_CSA       :  1.4    : CsA absorption rate (1/d)
F1_CSA       :  0.30   : CsA bioavailability
CL_CSA       :  27.0   : CsA clearance (L/d)
V_CSA        :  85.0   : CsA volume of distribution (L)

// -------- Ustekinumab (UST) --------
KA_UST       :  0.20   : UST absorption rate (1/d)
F1_UST       :  0.57   : UST bioavailability
CL_UST       :  0.45   : UST clearance (L/d)
V_UST        :  4.6    : UST volume of distribution (L)

// -------- Prednisone (PRED) --------
KA_PRED      :  4.0    : PRED absorption rate (1/d)
CL_PRED      :  90.0   : PRED clearance (L/d)
V_PRED       :  35.0   : PRED volume of distribution (L)

// -------- PD: cytokine homeostasis (turnover) --------
kin_TNF      :  3.0    : TNF synthesis (pg/mL/d)
kout_TNF     :  3.0    : TNF degradation (1/d)
kin_IL1      :  4.0    : IL-1b synthesis baseline
kout_IL1     :  4.0    : IL-1b degradation (1/d)
kin_IL17     :  2.0    : IL-17A synthesis baseline
kout_IL17    :  2.0    : IL-17A degradation (1/d)
kin_IL6      :  3.0    : IL-6 baseline
kout_IL6     :  3.0    : IL-6 degradation
kin_IL23     :  1.5    : IL-23 baseline
kout_IL23    :  1.5    : IL-23 degradation
kin_IL8      :  6.0    : IL-8 baseline
kout_IL8     :  6.0    : IL-8 degradation

// -------- Disease forcing (PG-specific) --------
DSEV         :  3.0    : Disease severity multiplier (1=remission, 3-5=active)
F_pathergy   :  0.5    : Pathergy event toggle (0/1)
F_comorb_IBD :  0.6    : IBD comorbidity coefficient

// -------- Inflammatory drivers (cross-cytokine) --------
k_TNF_drive   : 0.20   : disease → TNF
k_IL1_drive   : 0.25   : disease → IL-1
k_IL17_drive  : 0.20   : Th17/IL-23 → IL-17
k_IL6_drive   : 0.15   : TNF/IL-1 → IL-6
k_IL23_drive  : 0.10   : DC → IL-23
k_IL8_drive   : 0.40   : TNF+IL-1+IL-17 → IL-8

// -------- Cell dynamics --------
kin_Neu      : 0.6     : neutrophil influx baseline (1e9/L/d)
kout_Neu     : 0.5     : neutrophil clearance (1/d)
k_NETform    : 0.5     : NET formation per neutrophil unit
kout_NET     : 0.8     : NET degradation/clearance
kin_Th17     : 0.05    : Th17 baseline (au/d)
kout_Th17    : 0.05
kin_Treg     : 0.06
kout_Treg    : 0.06
kin_M1       : 0.20
kout_M1      : 0.20

// -------- Tissue / ulcer kinetics --------
k_ulcer_grow : 0.30    : ulcer growth from inflammation (cm²/d)
k_heal_max   : 0.18    : maximum healing rate (cm²/d)
EC50_heal_drug : 0.20  : drug effect on healing rate
UlcerArea_0  : 12.0    : baseline ulcer area cm²
k_MMP        : 0.5     : MMP-9 synthesis rate
kout_MMP     : 0.4     : MMP-9 degradation
k_ROS        : 0.6     : ROS gen
kout_ROS     : 0.5

// -------- Endpoints --------
k_CRP_syn    : 0.6     : CRP from IL-6 (mg/L per IL-6 unit)
kout_CRP     : 0.3
k_Calp_syn   : 0.3
kout_Calp    : 0.3

// -------- Drug effect potencies (EC50; suppression Emax = 0.9 by default) --------
IC50_ADA_TNF  : 0.20    : ug/mL on TNF
IC50_IFX_TNF  : 0.10
IC50_ANA_IL1  : 0.50    : ug/mL on IL-1 axis
IC50_UST_IL23 : 0.50
IC50_CSA_Th17 : 0.10    : ng/mL (Th17 inhibition)
IC50_PRED     : 5.0     : ng/mL on NF-kB-dependent cytokines (TNF/IL-1/IL-6)
Emax_drug     : 0.92

$CMT @annotated
ADA_SC    : adalimumab SC depot (mg)
ADA_CENT  : adalimumab central (mg)
ADA_PERI  : adalimumab peripheral (mg)
IFX_CENT  : infliximab central (mg)
IFX_PERI  : infliximab peripheral (mg)
ANA_SC    : anakinra SC depot (mg)
ANA_CENT  : anakinra central (mg)
CSA_GUT   : cyclosporine gut (mg)
CSA_CENT  : cyclosporine central (mg)
UST_SC    : ustekinumab SC depot (mg)
UST_CENT  : ustekinumab central (mg)
PRED_GUT  : prednisone gut (mg)
PRED_CENT : prednisone central (mg)
TNFa      : TNF-alpha (pg/mL)
IL1b      : IL-1beta (pg/mL)
IL17A     : IL-17A (pg/mL)
IL6       : IL-6 (pg/mL)
IL23      : IL-23 (pg/mL)
IL8       : IL-8 (pg/mL)
Neutroph  : tissue neutrophil index (au)
NET       : NET burden (au)
Th17      : Th17 index (au)
Treg      : Treg index (au)
M1        : macrophage M1 index (au)
MMP9_act  : active MMP-9 (au)
ROS_les   : lesional ROS (au)
Ulcer     : ulcer area (cm²)
Healed    : cumulative healed fraction (0-1)
CRP       : C-reactive protein (mg/L)
Calprot   : serum calprotectin (ng/mL)

$MAIN
ADA_CENT_0 = 0;
TNFa_0     = 3.0;
IL1b_0     = 4.0;
IL17A_0    = 2.0;
IL6_0      = 3.0;
IL23_0     = 1.5;
IL8_0      = 6.0;
Neutroph_0 = 1.0;
NET_0      = 0.5;
Th17_0     = 1.0;
Treg_0     = 1.0;
M1_0       = 1.0;
MMP9_act_0 = 1.0;
ROS_les_0  = 1.0;
Ulcer_0    = UlcerArea_0;
Healed_0   = 0.0;
CRP_0      = 6.0;
Calprot_0  = 100.0;

$ODE
// ===== Drug PK =====
double CONC_ADA  = ADA_CENT / V1_ADA;     // mg/L = ug/mL
double CONC_IFX  = IFX_CENT / V1_IFX;
double CONC_ANA  = ANA_CENT / V_ANA;
double CONC_CSA  = CSA_CENT / V_CSA * 1000.0; // ng/mL (mg/L*1000)
double CONC_UST  = UST_CENT / V_UST;
double CONC_PRED = PRED_CENT / V_PRED * 1000.0; // ng/mL

dxdt_ADA_SC    = -KA_ADA  * ADA_SC;
dxdt_ADA_CENT  =  KA_ADA  * F1_ADA * ADA_SC
                 - (CL_ADA/V1_ADA) * ADA_CENT
                 - (Q_ADA/V1_ADA) * ADA_CENT + (Q_ADA/V2_ADA) * ADA_PERI;
dxdt_ADA_PERI  =  (Q_ADA/V1_ADA) * ADA_CENT - (Q_ADA/V2_ADA) * ADA_PERI;

dxdt_IFX_CENT  = -(CL_IFX/V1_IFX) * IFX_CENT
                 - (Q_IFX/V1_IFX) * IFX_CENT + (Q_IFX/V2_IFX) * IFX_PERI;
dxdt_IFX_PERI  =  (Q_IFX/V1_IFX) * IFX_CENT - (Q_IFX/V2_IFX) * IFX_PERI;

dxdt_ANA_SC    = -KA_ANA  * ANA_SC;
dxdt_ANA_CENT  =  KA_ANA  * F1_ANA * ANA_SC - (CL_ANA/V_ANA) * ANA_CENT;

dxdt_CSA_GUT   = -KA_CSA  * CSA_GUT;
dxdt_CSA_CENT  =  KA_CSA  * F1_CSA * CSA_GUT - (CL_CSA/V_CSA) * CSA_CENT;

dxdt_UST_SC    = -KA_UST  * UST_SC;
dxdt_UST_CENT  =  KA_UST  * F1_UST * UST_SC - (CL_UST/V_UST) * UST_CENT;

dxdt_PRED_GUT  = -KA_PRED * PRED_GUT;
dxdt_PRED_CENT =  KA_PRED * PRED_GUT - (CL_PRED/V_PRED) * PRED_CENT;

// ===== PD: drug effects (Emax sigmoidal suppression) =====
double EFF_ADA_TNF  = Emax_drug * CONC_ADA  / (IC50_ADA_TNF + CONC_ADA);
double EFF_IFX_TNF  = Emax_drug * CONC_IFX  / (IC50_IFX_TNF + CONC_IFX);
double EFF_ANA_IL1  = Emax_drug * CONC_ANA  / (IC50_ANA_IL1 + CONC_ANA);
double EFF_UST_IL23 = Emax_drug * CONC_UST  / (IC50_UST_IL23 + CONC_UST);
double EFF_CSA_Th17 = Emax_drug * CONC_CSA  / (IC50_CSA_Th17 + CONC_CSA);
double EFF_PRED     = Emax_drug * CONC_PRED / (IC50_PRED + CONC_PRED);

// Combined TNF suppression (multiplicative escape)
double SUP_TNF  = 1.0 - (1.0 - (1-EFF_ADA_TNF)*(1-EFF_IFX_TNF)*(1-EFF_PRED)) ;
double SUP_IL1  = (1 - EFF_ANA_IL1) * (1 - EFF_PRED);
double SUP_IL17 = (1 - EFF_UST_IL23) * (1 - EFF_CSA_Th17);
double SUP_IL23 = (1 - EFF_UST_IL23);
double SUP_IL6  = (1 - EFF_PRED) * (1 - EFF_IFX_TNF * 0.6);

// Disease-driven forcing
double F_dis = DSEV * (1.0 + 0.5*F_pathergy + 0.4*F_comorb_IBD);

// ===== Cytokine ODEs (synthesis driven; drug-suppressed) =====
dxdt_TNFa  = (kin_TNF  + k_TNF_drive  * F_dis * (Neutroph + Th17 + M1))
               * (1 - EFF_ADA_TNF)*(1 - EFF_IFX_TNF)*(1 - EFF_PRED)
             - kout_TNF * TNFa;

dxdt_IL1b  = (kin_IL1  + k_IL1_drive  * F_dis * (Neutroph + M1)) * SUP_IL1
             - kout_IL1 * IL1b;

dxdt_IL17A = (kin_IL17 + k_IL17_drive * F_dis * Th17) * SUP_IL17
             - kout_IL17 * IL17A;

dxdt_IL6   = (kin_IL6  + k_IL6_drive  * F_dis * (TNFa + IL1b)/10.0) * SUP_IL6
             - kout_IL6 * IL6;

dxdt_IL23  = (kin_IL23 + k_IL23_drive * F_dis * M1) * SUP_IL23
             - kout_IL23 * IL23;

dxdt_IL8   = (kin_IL8  + k_IL8_drive  * F_dis * (TNFa + IL1b + IL17A)/15.0)
               * (1 - EFF_PRED*0.6)
             - kout_IL8 * IL8;

// ===== Cell ODEs =====
double rec_neu = IL8 / 6.0;            // recruitment scaled by IL-8
dxdt_Neutroph = kin_Neu + k_IL8_drive * 0.1 * rec_neu * F_dis - kout_Neu * Neutroph;
dxdt_NET      = k_NETform * Neutroph * (1 + ROS_les) - kout_NET * NET;
dxdt_Th17     = kin_Th17 + 0.1 * IL23/1.5 * F_dis * SUP_IL17 - kout_Th17 * Th17;
dxdt_Treg     = kin_Treg - kout_Treg * Treg - 0.02 * Th17;   // Th17 suppresses Treg
dxdt_M1       = kin_M1   + 0.05 * (TNFa + IL1b)/7.0 * F_dis - kout_M1 * M1;

// ===== Tissue: MMP-9, ROS, ulcer =====
dxdt_MMP9_act = k_MMP * (Neutroph + 0.3 * M1) - kout_MMP * MMP9_act;
dxdt_ROS_les  = k_ROS * Neutroph - kout_ROS * ROS_les;

// Ulcer growth from MMP-9 / ROS / NET ; healing accelerated by drugs (cytokine suppression)
double inflam_drive = (MMP9_act + 0.6*ROS_les + 0.5*NET);
double drug_heal_eff = EFF_ADA_TNF + EFF_IFX_TNF + EFF_ANA_IL1 + EFF_UST_IL23 + EFF_CSA_Th17 + 0.5*EFF_PRED;
double heal_rate = k_heal_max * drug_heal_eff / (EC50_heal_drug + drug_heal_eff);

dxdt_Ulcer = k_ulcer_grow * inflam_drive * (Ulcer / (Ulcer + 5)) - heal_rate * Ulcer;
dxdt_Healed = heal_rate * Ulcer / UlcerArea_0;   // fraction healed

// ===== Endpoints =====
dxdt_CRP     = k_CRP_syn * IL6 - kout_CRP * CRP;
dxdt_Calprot = k_Calp_syn * (Neutroph + 0.5*NET) - kout_Calp * Calprot;

$TABLE
double Pain_VAS   = std::min(10.0, 1.5 + 0.4*Ulcer + 0.1*TNFa);
double PARACELSUS = std::min(60.0, 5 + 1.5*Ulcer + 0.5*Pain_VAS + 2.0*(Neutroph>2 ? 1 : 0)*5
                                  + 0.8*M1 + 0.4*Th17);
double DLQI       = std::min(30.0, 4 + 0.6*Ulcer + 0.3*Pain_VAS);
double HiSCRpseudo = (Ulcer < 0.5 * UlcerArea_0) ? 1.0 : 0.0;  // ≥50% reduction
double CompleteHeal = (Ulcer < 0.05 * UlcerArea_0) ? 1.0 : 0.0;

$CAPTURE @annotated
CONC_ADA  : Adalimumab conc (ug/mL)
CONC_IFX  : Infliximab conc (ug/mL)
CONC_ANA  : Anakinra conc (ug/mL)
CONC_CSA  : Cyclosporine conc (ng/mL)
CONC_UST  : Ustekinumab conc (ug/mL)
CONC_PRED : Prednisone conc (ng/mL)
Pain_VAS  : VAS pain 0-10
PARACELSUS : PARACELSUS score
DLQI      : DLQI score
HiSCRpseudo : >=50% ulcer reduction
CompleteHeal : ulcer cleared
'

# -----------------------------------------------------------------------------
# Compile model
# -----------------------------------------------------------------------------
pg_qsp_model <- mcode("pg_qsp", pg_qsp_code)

# -----------------------------------------------------------------------------
# Helper: simulate a scenario over 168 days (24 weeks) with a defined dose regimen
# -----------------------------------------------------------------------------
simulate_pg <- function(model, scenario, tmax = 168, ...) {
  ev <- switch(scenario,
    # 1) Standard of care: Prednisone 60 mg/d ×14d → taper × 12 weeks
    "Prednisone_SOC" = ev(amt = 60, cmt = "PRED_GUT", ii = 1, addl = 13) +
                       ev(amt = 40, cmt = "PRED_GUT", time = 14, ii = 1, addl = 13) +
                       ev(amt = 20, cmt = "PRED_GUT", time = 28, ii = 1, addl = 27) +
                       ev(amt = 10, cmt = "PRED_GUT", time = 56, ii = 1, addl = 55),
    # 2) Cyclosporine: 4 mg/kg/d (≈ 280 mg) ×6 mo
    "Cyclosporine"   = ev(amt = 140, cmt = "CSA_GUT", ii = 0.5, addl = 335),  # BID
    # 3) Infliximab induction 5 mg/kg @0,2,6 wk then q8w (assume 70 kg = 350 mg)
    "Infliximab"     = ev(amt = 350, cmt = "IFX_CENT", time = 0) +
                       ev(amt = 350, cmt = "IFX_CENT", time = 14) +
                       ev(amt = 350, cmt = "IFX_CENT", time = 42) +
                       ev(amt = 350, cmt = "IFX_CENT", time = 98) +
                       ev(amt = 350, cmt = "IFX_CENT", time = 154),
    # 4) Adalimumab: 80 mg loading → 40 mg q1w
    "Adalimumab"     = ev(amt = 80, cmt = "ADA_SC", time = 0) +
                       ev(amt = 40, cmt = "ADA_SC", time = 7, ii = 7, addl = 22),
    # 5) Anakinra 100 mg/d SC
    "Anakinra"       = ev(amt = 100, cmt = "ANA_SC", ii = 1, addl = 167),
    # 6) Ustekinumab 90 mg SC q12w (after 0 wk + 4 wk loading)
    "Ustekinumab"    = ev(amt = 90, cmt = "UST_SC", time = 0) +
                       ev(amt = 90, cmt = "UST_SC", time = 28) +
                       ev(amt = 90, cmt = "UST_SC", time = 112),
    # 7) Combo: Prednisone short course + Cyclosporine maintenance
    "Combo_PRED_CSA" = ev(amt = 60, cmt = "PRED_GUT", ii = 1, addl = 13) +
                       ev(amt = 30, cmt = "PRED_GUT", time = 14, ii = 1, addl = 27) +
                       ev(amt = 140, cmt = "CSA_GUT", time = 0, ii = 0.5, addl = 335),
    # 8) Combo: Infliximab + Methotrexate-like Th17 suppression simulated via CsA low-dose
    "Combo_IFX_low_CsA" = ev(amt = 350, cmt = "IFX_CENT", time = 0) +
                          ev(amt = 350, cmt = "IFX_CENT", time = 14) +
                          ev(amt = 350, cmt = "IFX_CENT", time = 42) +
                          ev(amt = 350, cmt = "IFX_CENT", time = 98) +
                          ev(amt = 70, cmt = "CSA_GUT", ii = 0.5, addl = 335),
    # 9) No treatment (natural history)
    "NoTreatment"    = ev(amt = 0, cmt = "PRED_GUT", time = 0)
  )

  out <- model %>%
    ev(ev) %>%
    mrgsim(end = tmax, delta = 0.5)

  as.data.frame(out) %>% mutate(scenario = scenario)
}

# -----------------------------------------------------------------------------
# Run 9 scenarios (use any subset as needed)
# -----------------------------------------------------------------------------
scenarios <- c("NoTreatment", "Prednisone_SOC", "Cyclosporine",
               "Infliximab", "Adalimumab", "Anakinra",
               "Ustekinumab", "Combo_PRED_CSA", "Combo_IFX_low_CsA")

if (interactive() || identical(Sys.getenv("RUN_PG_SIM"), "1")) {
  res <- bind_rows(lapply(scenarios, function(s)
    simulate_pg(pg_qsp_model, s, tmax = 168)))

  message("Endpoint summary @ Week 24:")
  res %>% filter(abs(time - 168) < 0.6) %>%
    select(scenario, Ulcer, Healed, PARACELSUS, DLQI, Pain_VAS, CRP) %>%
    print()
}

# -----------------------------------------------------------------------------
# Virtual population helper — sample CL_ADA, IC50, baseline ulcer area, comorbidities
# -----------------------------------------------------------------------------
make_vpop <- function(n = 200, seed = 42) {
  set.seed(seed)
  data.frame(
    ID            = seq_len(n),
    CL_ADA        = rlnorm(n, log(0.31), 0.25),
    CL_IFX        = rlnorm(n, log(0.32), 0.30),
    CL_CSA        = rlnorm(n, log(27),   0.25),
    UlcerArea_0   = pmax(2, rlnorm(n, log(10), 0.5)),
    DSEV          = pmax(1, rnorm(n, 3.0, 0.7)),
    F_pathergy    = rbinom(n, 1, 0.45),
    F_comorb_IBD  = rbinom(n, 1, 0.40)
  )
}

# Example: run virtual population on adalimumab
run_vpop_scenario <- function(model, scenario = "Adalimumab",
                              n = 50, tmax = 168) {
  vp <- make_vpop(n = n)
  bind_rows(lapply(seq_len(n), function(i) {
    m <- model %>% param(as.list(vp[i, c("CL_ADA","CL_IFX","CL_CSA","UlcerArea_0",
                                         "DSEV","F_pathergy","F_comorb_IBD")]))
    simulate_pg(m, scenario, tmax = tmax) %>% mutate(ID = i)
  }))
}

cat("[PG QSP] Model compiled.  Scenarios available:\n  ",
    paste(scenarios, collapse = "\n   "),
    "\nCall simulate_pg(pg_qsp_model, '<scenario>') to generate trajectories.\n")
