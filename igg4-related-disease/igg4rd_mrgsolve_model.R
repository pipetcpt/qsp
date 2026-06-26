## ============================================================
## IgG4-Related Disease (IgG4-RD) QSP Model
## mrgsolve ODE Implementation
##
## Disease: IgG4-RD — systemic fibroinflammatory condition
## Mechanism: Tfh2-driven IgG4+ plasmablasts, M2-TGF-β fibrosis
##
## Drugs modelled:
##   1. Rituximab (anti-CD20, IV)
##   2. Prednisone (glucocorticoid, oral)
##   3. Dupilumab (anti-IL-4Rα, SC)  [investigational]
##
## ODE Compartments (23 total):
##   PK: CENT_RTX, PERI_RTX, CD20_FREE, RTX_CD20   (4)
##       GUT_PRED, CENT_PRED                          (2)
##       SC_DUP, CENT_DUP, IL4RA_FREE, DUP_IL4RA     (4)
##   Immunology: BNV, GCB, PB, PC, TFH2, CTL4       (6)
##   Cytokines:  IgG4_SER, IL4, IL10, TGFB           (4)
##   Fibrosis/Damage: MYOFIB, ECM, IRI               (3)
##
## Key Calibration Sources:
##   - Khosroshahi 2012 Ann Rheum Dis: RTX IgG4-RD (n=10)
##   - Carruthers 2015 Ann Rheum Dis: RTX responder (n=30)
##   - Lanzillotta 2020 Lancet Rheum: RTX vs GC
##   - Hart 2021 N Engl J Med: RTX in IgG4-RD
##   - Perugino 2021 N Engl J Med: IgG4-RD pathogenesis
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## ── Model code ──────────────────────────────────────────────────────────────
code <- '
$PROB
IgG4-Related Disease (IgG4-RD) QSP Model
23-compartment ODE: Drug PK + Immunology + Fibrosis

$PARAM
// ── Rituximab PK (TMDD 2-compartment) ─────────────────────────────────────
CL_RTX    = 0.35    // L/day, nonspecific clearance (TMDD offset)
V1_RTX    = 3.0     // L, central volume
Q_RTX     = 0.8     // L/day, intercompartmental CL
V2_RTX    = 4.5     // L, peripheral volume
kon_RTX   = 0.55    // 1/(nM·day), RTX-CD20 on-rate
koff_RTX  = 0.002   // 1/day, RTX-CD20 off-rate
CD20_ss   = 180.0   // nM, baseline CD20 on B cells
ksyn_CD20 = 90.0    // nM/day, CD20 synthesis rate (= kd × CD20_ss)
kdeg_CD20 = 0.5     // 1/day, CD20 baseline degradation
kint_RTX  = 0.3     // 1/day, RTX-CD20 complex internalisation rate

// ── Prednisone PK (1-compartment oral) ────────────────────────────────────
ka_PRED   = 3.0     // 1/day, absorption
CL_PRED   = 18.0    // L/day, hepatic clearance (bioactivation → prednisolone)
V_PRED    = 40.0    // L, Vd
F_PRED    = 0.80    // oral bioavailability

// ── Dupilumab PK (SC, 2-compartment) ──────────────────────────────────────
ka_DUP    = 0.08    // 1/day, SC absorption
CL_DUP    = 0.25    // L/day, nonspecific CL
V1_DUP    = 4.5     // L
kon_DUP   = 0.30    // 1/(nM·day)
koff_DUP  = 0.001   // 1/day
IL4RA_ss  = 10.0    // nM, baseline IL-4Rα
ksyn_IL4RA = 0.5    // nM/day
kdeg_IL4RA = 0.05   // 1/day
kint_DUP  = 0.15    // 1/day, DUP-IL4Rα complex internalisation

// ── Immunology: B cells (relative units, 1.0 = normal) ────────────────────
kprolif_BNV = 0.05  // 1/day, naïve B proliferation
kdeath_BNV  = 0.05  // 1/day, naïve B turnover
kgc_BNV     = 0.02  // 1/day, naïve B → GC B (antigen-driven)
kprolif_GCB = 0.20  // 1/day, GC B proliferation (Tfh2-driven)
kdeath_GCB  = 0.18  // 1/day, GC B apoptosis (selection)
kpb_GCB     = 0.10  // 1/day, GC B → plasmablast
kmat_PB     = 0.30  // 1/day, plasmablast → plasma cell maturation
kdeath_PB   = 0.40  // 1/day, plasmablast turnover (short-lived, t½~2d)
kdeath_PC   = 0.005 // 1/day, plasma cell turnover (long-lived, t½~4-6mo)

// ── Tfh2 & CTL4 dynamics ──────────────────────────────────────────────────
kprolif_TFH2 = 0.08 // 1/day, Tfh2 expansion (IL-4-driven)
kdeath_TFH2  = 0.05 // 1/day, Tfh2 turnover
IL4_EC50_TFH2 = 0.5 // relative IL-4 EC50 for Tfh2 expansion
kprolif_CTL4 = 0.04 // 1/day, CD4 CTL expansion
kdeath_CTL4  = 0.04 // 1/day

// ── Cytokine dynamics ─────────────────────────────────────────────────────
kprod_IL4   = 0.10  // relative units/day, IL-4 production (Tfh2)
kdeg_IL4    = 1.0   // 1/day, IL-4 degradation (t½~16h)
kprod_IL10  = 0.06  // relative units/day, IL-10 production (Treg/Tfh2)
kdeg_IL10   = 0.8   // 1/day
kprod_TGFB  = 0.08  // relative units/day, TGF-β from CTL4/M2
kdeg_TGFB   = 0.5   // 1/day

// ── IgG4 dynamics ─────────────────────────────────────────────────────────
ksec_IgG4   = 0.15  // mg/dL/day per PC unit, IgG4 secretion
kdeg_IgG4   = 0.015 // 1/day, IgG4 catabolism (t½~21-28d)
IgG4_norm   = 60.0  // mg/dL, upper normal limit
IgG4_0      = 450.0 // mg/dL, typical active IgG4-RD at baseline

// ── Fibrosis dynamics ─────────────────────────────────────────────────────
kact_MYOFIB = 0.03  // 1/day, fibroblast activation by TGF-β
kinact_MYOFIB = 0.02 // 1/day, myofibroblast spontaneous inactivation
TGFB_EC50_FIB = 1.5  // EC50 for TGF-β → myofibroblast
kcol_ECM    = 0.04  // 1/day, ECM synthesis by myofibroblast
kdeg_ECM    = 0.005 // 1/day, spontaneous ECM remodelling

// ── Drug effect parameters ────────────────────────────────────────────────
Emax_RTX    = 0.98   // max B cell depletion by RTX
EC50_RTX    = 50.0   // nM RTX (central) for 50% B cell kill
Emax_PRED   = 0.85   // max immunosuppression by GC
EC50_PRED   = 80.0   // nM prednisolone
Emax_DUP    = 0.90   // max IL-4/IL-13 blockade
// Dupilumab PD already mechanistic via IL4RA complex

// ── IRI calculation weights ───────────────────────────────────────────────
w_IgG4  = 0.30
w_ECM   = 0.35
w_PC    = 0.20
w_TFH2  = 0.15

$INIT
// Rituximab PK
CENT_RTX  = 0
PERI_RTX  = 0
CD20_FREE = 180.0   // nM (= CD20_ss)
RTX_CD20  = 0

// Prednisone PK
GUT_PRED  = 0
CENT_PRED = 0

// Dupilumab PK
SC_DUP    = 0
CENT_DUP  = 0
IL4RA_FREE = 10.0   // nM (= IL4RA_ss)
DUP_IL4RA = 0

// Immunology (relative units; 1.0 = normal steady-state)
BNV   = 1.0         // naïve B cells
GCB   = 2.5         // GC B cells (elevated in IgG4-RD)
PB    = 3.0         // plasmablasts (elevated)
PC    = 4.0         // long-lived plasma cells (elevated)
TFH2  = 3.5         // Tfh2 cells (markedly elevated)
CTL4  = 2.0         // cytotoxic CD4+ T cells (elevated)

// Cytokines (relative; 1.0 = normal)
IgG4_SER = 450.0    // mg/dL serum IgG4 (active disease)
IL4       = 3.0     // relative units (elevated)
IL10      = 2.5     // relative units (elevated)
TGFB      = 3.0     // relative units (elevated)

// Fibrosis/Damage
MYOFIB = 2.0        // myofibroblast activation (elevated)
ECM    = 2.5        // ECM fibrosis index (elevated)
IRI    = 8.0        // IgG4-RD Responder Index (0-24, typically 8-16 at diagnosis)

$ODE
// ─── Rituximab PK (TMDD) ─────────────────────────────────────────────────
double kon_r = kon_RTX * CENT_RTX * CD20_FREE;
double koff_r = koff_RTX * RTX_CD20;
double kint_r = kint_RTX * RTX_CD20;

dxdt_CENT_RTX = - CL_RTX/V1_RTX * CENT_RTX
                - Q_RTX/V1_RTX * CENT_RTX
                + Q_RTX/V2_RTX * PERI_RTX
                - kon_r + koff_r;
dxdt_PERI_RTX = Q_RTX/V1_RTX * CENT_RTX - Q_RTX/V2_RTX * PERI_RTX;
dxdt_CD20_FREE = ksyn_CD20 - kdeg_CD20 * CD20_FREE - kon_r + koff_r;
dxdt_RTX_CD20  = kon_r - koff_r - kint_r;

// ─── Prednisone PK ────────────────────────────────────────────────────────
dxdt_GUT_PRED  = - ka_PRED * GUT_PRED;
dxdt_CENT_PRED = ka_PRED * GUT_PRED * F_PRED
                 - CL_PRED/V_PRED * CENT_PRED;

// ─── Dupilumab PK (TMDD) ──────────────────────────────────────────────────
double kon_d = kon_DUP * CENT_DUP * IL4RA_FREE;
double koff_d = koff_DUP * DUP_IL4RA;
double kint_d = kint_DUP * DUP_IL4RA;

dxdt_SC_DUP    = - ka_DUP * SC_DUP;
dxdt_CENT_DUP  = ka_DUP * SC_DUP
                 - CL_DUP/V1_DUP * CENT_DUP
                 - kon_d + koff_d;
dxdt_IL4RA_FREE = ksyn_IL4RA - kdeg_IL4RA * IL4RA_FREE - kon_d + koff_d;
dxdt_DUP_IL4RA  = kon_d - koff_d - kint_d;

// ─── Drug effect calculations ─────────────────────────────────────────────
// RTX: fraction of B cells killed (TMDD-driven; CD20 occupancy)
double CD20_occ = RTX_CD20 / (RTX_CD20 + CD20_FREE + 0.001);
double Ekill_RTX = Emax_RTX * CD20_occ;

// Prednisone: immunosuppressive effect (0-1)
double Eimmu_PRED = Emax_PRED * CENT_PRED / (EC50_PRED + CENT_PRED);

// Dupilumab: IL-4Rα occupancy → IL-4/13 blockade fraction
double IL4RA_occ = DUP_IL4RA / (DUP_IL4RA + IL4RA_FREE + 0.001);

// ─── B Cell Dynamics ──────────────────────────────────────────────────────
// Naïve B cells
double BNV_death = kdeath_BNV * BNV * (1 + Ekill_RTX + Eimmu_PRED);
dxdt_BNV = kprolif_BNV * BNV
           - kdeath_BNV * BNV
           - kgc_BNV * BNV
           - Ekill_RTX * kprolif_BNV * BNV;

// GC B cells: expanded by Tfh2, killed by RTX
double TFH2_stim = TFH2 / (1.0 + TFH2);  // saturating Tfh2 stimulation
dxdt_GCB = kgc_BNV * BNV
           + kprolif_GCB * GCB * TFH2_stim * (1 - IL4RA_occ * 0.5)
           - kdeath_GCB * GCB
           - kpb_GCB * GCB
           - Ekill_RTX * kprolif_GCB * GCB
           - Eimmu_PRED * 0.8 * GCB;

// Plasmablasts (short-lived)
dxdt_PB = kpb_GCB * GCB
          - kmat_PB * PB
          - kdeath_PB * PB
          - Ekill_RTX * (kmat_PB + kdeath_PB) * PB
          - Eimmu_PRED * 0.9 * PB;

// Long-lived plasma cells (RTX-resistant once CD20-)
double PC_RTX_kill = Ekill_RTX * 0.15;  // partial kill (few CD20+)
dxdt_PC = kmat_PB * PB
          - kdeath_PC * PC
          - PC_RTX_kill * PC
          - Eimmu_PRED * 0.2 * PC;

// ─── T Cell Dynamics ──────────────────────────────────────────────────────
// Tfh2: driven by IL-4 (autocrine loop), suppressed by GC
double IL4_drive = IL4 / (IL4_EC50_TFH2 + IL4);
dxdt_TFH2 = kprolif_TFH2 * TFH2 * IL4_drive * (1 - IL4RA_occ)
            - kdeath_TFH2 * TFH2
            - Eimmu_PRED * 0.85 * TFH2;

// CTL4 (SLAMF7+ CD4 CTL)
dxdt_CTL4 = kprolif_CTL4 * CTL4
            - kdeath_CTL4 * CTL4
            - Eimmu_PRED * 0.7 * CTL4;

// ─── Cytokine Dynamics ────────────────────────────────────────────────────
// IL-4: from Tfh2 + Th2 cells, blocked by dupilumab
dxdt_IL4 = kprod_IL4 * TFH2 * (1 - IL4RA_occ * 0.6)
           - kdeg_IL4 * IL4
           - Eimmu_PRED * 0.6 * IL4;

// IL-10: from Treg + Tfh2, immunosuppressive
dxdt_IL10 = kprod_IL10 * (TFH2 + 0.5)
            - kdeg_IL10 * IL10
            - Eimmu_PRED * 0.3 * IL10;

// TGF-β: from CTL4 + M2 macrophages (modelled as function of CTL4)
dxdt_TGFB = kprod_TGFB * CTL4 * (1 - IL4RA_occ * 0.3)
             - kdeg_TGFB * TGFB
             - Eimmu_PRED * 0.4 * TGFB;

// ─── IgG4 Serum Level ────────────────────────────────────────────────────
// Produced by PC, switch driven by IL-4/IL-10, catabolic t½ ~21d
double IgG4_switch = IL4 * IL10 / (1 + IL4 * IL10);  // synergistic switch
dxdt_IgG4_SER = ksec_IgG4 * PC * IgG4_switch
                - kdeg_IgG4 * IgG4_SER;

// ─── Fibrosis Dynamics ────────────────────────────────────────────────────
// Myofibroblast activation: TGF-β drives, GC suppresses
double TGFB_eff = TGFB^2 / (TGFB_EC50_FIB^2 + TGFB^2);
dxdt_MYOFIB = kact_MYOFIB * TGFB_eff
              - kinact_MYOFIB * MYOFIB
              - Eimmu_PRED * 0.35 * MYOFIB;

// ECM: synthesized by myofibroblast
dxdt_ECM = kcol_ECM * MYOFIB
           - kdeg_ECM * ECM
           - Eimmu_PRED * 0.2 * ECM;

// ─── IRI (IgG4-RD Responder Index) ──────────────────────────────────────
// Composite score 0-24 scale (0=complete response, 24=severe active)
double IgG4_comp = (IgG4_SER / IgG4_0) * 10.0 * w_IgG4;
double ECM_comp  = (ECM / 2.5)          * 10.0 * w_ECM;
double PC_comp   = (PC  / 4.0)          * 10.0 * w_PC;
double TFH2_comp = (TFH2/ 3.5)          * 10.0 * w_TFH2;
double IRI_new = IgG4_comp + ECM_comp + PC_comp + TFH2_comp;
// Clamp 0-24
if(IRI_new < 0) IRI_new = 0;
if(IRI_new > 24) IRI_new = 24;
dxdt_IRI = (IRI_new - IRI) * 0.5;  // exponential moving average

$TABLE
double RTX_conc_ug_mL = CENT_RTX * 148000 / (1e6) * 1000; // nM → μg/mL (MW~148kDa)
double PRED_conc_nM   = CENT_PRED;
double DUP_conc_ug_mL = CENT_DUP * 146000 / (1e6) * 1000;
double Bcell_pct      = BNV * 100;  // % of baseline
double IgG4_mgdL      = IgG4_SER;
double Fibrosis_idx   = ECM;
double Activity_IRI   = IRI;
double CD20_occ_pct   = RTX_CD20 / (RTX_CD20 + CD20_FREE + 0.001) * 100;
double IL4RA_occ_pct  = DUP_IL4RA / (DUP_IL4RA + IL4RA_FREE + 0.001) * 100;
double TGFB_rel       = TGFB;
double PC_rel         = PC;
double TFH2_rel       = TFH2;
double CR_flag        = (IRI < 1.0) ? 1.0 : 0.0; // complete response
double PR_flag        = (IRI >= 1.0 && IRI < 4.0) ? 1.0 : 0.0;

$CAPTURE RTX_conc_ug_mL PRED_conc_nM DUP_conc_ug_mL Bcell_pct
         IgG4_mgdL Fibrosis_idx Activity_IRI CD20_occ_pct
         IL4RA_occ_pct TGFB_rel PC_rel TFH2_rel CR_flag PR_flag
'

## ── Compile model ──────────────────────────────────────────────────────────
mod <- mcode("IgG4RD_QSP", code)

## ── Helper: mg → nM concentration ─────────────────────────────────────────
mg_to_nM <- function(dose_mg, volume_L, MW_kDa) {
  (dose_mg / (MW_kDa * 1000)) * 1e9 / volume_L  # nM
}

## ── Dosing regimens ────────────────────────────────────────────────────────

# Rituximab IV: Khosroshahi protocol (1g IV D1, D15)
# MW ~148 kDa; dose expressed as μg bolus into CENT_RTX (nM)
RTX_1g_nM <- mg_to_nM(1000, 3.0, 148)  # ~2252 nM
RTX_375_nM <- function(BSA_m2) mg_to_nM(375 * BSA_m2, 3.0, 148)

# Prednisone 40mg/d oral → mg in GUT_PRED
# 40 mg × F=0.8 → bioavailability handled inside ODE
# dose in μg for nM calculation: prednisone MW = 358.4 Da
PRED_40mg_dose <- 40 * 1e6 / 358.4 / 40.0  # nM after absorption to Vd

# Dupilumab 300mg SC q2w → SC_DUP (nM)
DUP_300mg_nM <- mg_to_nM(300, 4.5, 146)  # ~457 nM

## ── Treatment scenarios ────────────────────────────────────────────────────

## Scenario 1: Untreated natural history (24 months)
e1 <- ev(time = 0, cmt = "CENT_RTX", amt = 0)  # no drug
out1 <- mod %>%
  mrgsim(events = e1, end = 730, delta = 1) %>%
  as.data.frame() %>%
  mutate(Scenario = "S1: Untreated Natural History")

## Scenario 2: Prednisone 40mg/d × 4wk taper over 6mo
# Approximate as fixed-dose then taper
pred_doses <- bind_rows(
  # 40mg/d × 4wk = 28 days
  lapply(0:27, function(d) ev(time = d, cmt = "GUT_PRED",
                               amt = 40 * 1e6 / 358.4 / 40.0,
                               rate = -2)),
  # 30mg/d wk 5-8
  lapply(28:55, function(d) ev(time = d, cmt = "GUT_PRED",
                                amt = 30 * 1e6 / 358.4 / 40.0,
                                rate = -2)),
  # taper 20→10→5mg/d by mo 3-6
  lapply(56:167, function(d) ev(time = d, cmt = "GUT_PRED",
                                 amt = max(5, 20 - floor((d-56)/28)*5) *
                                   1e6 / 358.4 / 40.0,
                                 rate = -2))
)
out2 <- mod %>%
  mrgsim(events = do.call(c, pred_doses), end = 730, delta = 1) %>%
  as.data.frame() %>%
  mutate(Scenario = "S2: Prednisone 40mg/d Taper")

## Scenario 3: Rituximab 1g IV × 2 doses (D1, D15) — Khosroshahi 2012
e3 <- ev(time = c(0, 14), cmt = "CENT_RTX", amt = RTX_1g_nM)
out3 <- mod %>%
  mrgsim(events = e3, end = 730, delta = 1) %>%
  as.data.frame() %>%
  mutate(Scenario = "S3: RTX 1g IV D1+D15 (Khosroshahi)")

## Scenario 4: Rituximab 375 mg/m² × 4 weekly (oncology protocol)
BSA <- 1.7
e4 <- ev(time = c(0, 7, 14, 21), cmt = "CENT_RTX",
         amt = RTX_375_nM(BSA))
out4 <- mod %>%
  mrgsim(events = e4, end = 730, delta = 1) %>%
  as.data.frame() %>%
  mutate(Scenario = "S4: RTX 375mg/m² ×4 Weekly")

## Scenario 5: Rituximab + maintenance 500mg q6m × 2 years
e5_induct <- ev(time = c(0, 14), cmt = "CENT_RTX",
                amt = mg_to_nM(1000, 3.0, 148))
e5_maint  <- ev(time = c(180, 360), cmt = "CENT_RTX",
                amt = mg_to_nM(500, 3.0, 148))
e5 <- c(e5_induct, e5_maint)
out5 <- mod %>%
  mrgsim(events = e5, end = 730, delta = 1) %>%
  as.data.frame() %>%
  mutate(Scenario = "S5: RTX 1g D1+D15 + Maintenance 500mg q6m")

## Scenario 6: Dupilumab 300mg SC q2w (investigational)
dup_times <- seq(0, 700, by = 14)
e6 <- ev(time = dup_times, cmt = "SC_DUP", amt = DUP_300mg_nM)
out6 <- mod %>%
  mrgsim(events = e6, end = 730, delta = 1) %>%
  as.data.frame() %>%
  mutate(Scenario = "S6: Dupilumab 300mg SC q2w (Investigational)")

## ── Combine results ────────────────────────────────────────────────────────
all_out <- bind_rows(out1, out2, out3, out4, out5, out6) %>%
  mutate(Scenario = factor(Scenario))

## ── Summary Statistics ─────────────────────────────────────────────────────
summary_stats <- all_out %>%
  group_by(Scenario) %>%
  summarise(
    IgG4_baseline   = first(IgG4_mgdL),
    IgG4_wk12       = IgG4_mgdL[which.min(abs(time - 84))],
    IgG4_wk24       = IgG4_mgdL[which.min(abs(time - 168))],
    IgG4_nadir      = min(IgG4_mgdL),
    IRI_baseline    = first(Activity_IRI),
    IRI_wk12        = Activity_IRI[which.min(abs(time - 84))],
    IRI_wk24        = Activity_IRI[which.min(abs(time - 168))],
    CR_pct_wk24     = mean(CR_flag[time >= 155 & time <= 181]) * 100,
    Bcell_nadir_pct = min(Bcell_pct),
    .groups = "drop"
  )

cat("\n=== IgG4-RD QSP Model — Treatment Summary (2-Year Simulation) ===\n")
print(as.data.frame(summary_stats), digits = 3)

## ── Calibration Validation ────────────────────────────────────────────────
cat("\n=== Calibration Check (Published Clinical Trial Data) ===\n")
cat("Khosroshahi 2012 (RTX 1g ×2): IgG4 fall ~75-80% at 3mo\n")
rtx_result <- all_out %>%
  filter(Scenario == "S3: RTX 1g IV D1+D15 (Khosroshahi)", time == 84) %>%
  mutate(IgG4_pct_change = (IgG4_mgdL - 450) / 450 * 100)
cat(sprintf("  Model: IgG4 at 3mo = %.1f mg/dL (%.1f%% change)\n",
            rtx_result$IgG4_mgdL[1], rtx_result$IgG4_pct_change[1]))

cat("Carruthers 2015 (RTX): B cell nadir <5/μL (95% depletion)\n")
rtx_bcell <- all_out %>%
  filter(Scenario == "S3: RTX 1g IV D1+D15 (Khosroshahi)") %>%
  summarise(nadir = min(Bcell_pct))
cat(sprintf("  Model: B cell nadir = %.1f%% of baseline\n", rtx_bcell$nadir))

cat("Lanzillotta 2020 (GC): IRI response, 75% respond by wk12\n")
gc_result <- all_out %>%
  filter(Scenario == "S2: Prednisone 40mg/d Taper", time == 84)
cat(sprintf("  Model: IRI at 3mo = %.1f (baseline 8.0)\n",
            gc_result$Activity_IRI[1]))

## ── Plots ─────────────────────────────────────────────────────────────────
p1 <- ggplot(all_out, aes(x = time/30.4, y = IgG4_mgdL, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 135, linetype = "dashed", color = "black") +
  annotate("text", x = 0.5, y = 140, label = "ULN 135 mg/dL", size = 3) +
  labs(title = "IgG4-RD: Serum IgG4 Over Time by Treatment",
       x = "Time (months)", y = "Serum IgG4 (mg/dL)") +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7)) +
  scale_color_brewer(palette = "Dark2")

p2 <- ggplot(all_out, aes(x = time/30.4, y = Activity_IRI, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 3, linetype = "dashed", color = "orange") +
  annotate("text", x = 0.5, y = 3.3, label = "PR threshold (IRI<3)", size = 3) +
  labs(title = "IgG4-RD: Disease Activity (IRI) Over Time",
       x = "Time (months)", y = "IgG4-RD Responder Index (0-24)") +
  ylim(0, 10) +
  theme_bw() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7)) +
  scale_color_brewer(palette = "Dark2")

p3 <- ggplot(all_out %>% filter(grepl("RTX", Scenario)),
             aes(x = time/30.4, y = Bcell_pct, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 5, linetype = "dashed") +
  labs(title = "B Cell Depletion by Rituximab Regimen",
       x = "Time (months)", y = "B Cells (% of baseline)") +
  theme_bw() +
  theme(legend.position = "bottom")

p4 <- ggplot(all_out, aes(x = time/30.4, y = Fibrosis_idx, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  labs(title = "Fibrosis Index (ECM) Over Time",
       x = "Time (months)", y = "ECM Fibrosis Index (relative)") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2")

p5 <- ggplot(all_out %>% filter(Scenario == "S6: Dupilumab 300mg SC q2w (Investigational)"),
             aes(x = time/30.4)) +
  geom_line(aes(y = IL4RA_occ_pct, color = "IL-4Rα Occupancy (%)"), linewidth=1.2) +
  geom_line(aes(y = IgG4_mgdL/10, color = "IgG4 (÷10, mg/dL)"), linewidth=1.2) +
  labs(title = "Dupilumab: IL-4Rα Occupancy & IgG4 Response",
       x = "Time (months)", y = "Value") +
  theme_bw() +
  scale_color_manual(values = c("IL-4Rα Occupancy (%)" = "#2196F3",
                                 "IgG4 (÷10, mg/dL)" = "#FF5722"))

gridExtra_ok <- requireNamespace("gridExtra", quietly = TRUE)
if (gridExtra_ok) {
  library(gridExtra)
  grid.arrange(p1, p2, p3, p4, ncol = 2)
} else {
  print(p1); print(p2); print(p3); print(p4); print(p5)
}

cat("\nModel simulation complete. 23-compartment IgG4-RD QSP model.\n")
cat("Scenarios 1-6 simulated over 730 days (24 months).\n")
