## =============================================================================
## Juvenile Idiopathic Arthritis (JIA) — QSP Model (mrgsolve)
## =============================================================================
## Disease: Juvenile Idiopathic Arthritis (JIA) — polyarticular & systemic subtypes
## Key drugs: MTX, Etanercept, Tocilizumab, Canakinumab, Abatacept, Baricitinib
## Author: QSP Disease Model Library (automated CCR session)
## Date: 2026-06-25
## Parameters calibrated from:
##   - Lovell et al. (NEJM 1998) — Etanercept pediatric trial
##   - De Benedetti et al. (NEJM 2012) — Tocilizumab sJIA trial
##   - Ruperto et al. (NEJM 2012) — Canakinumab sJIA trial
##   - Consolaro et al. (Arthritis Rheum 2009) — JADAS-27 validation
##   - Giannini et al. (Arthritis Rheum 1997) — ACR Pediatric criteria
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(patchwork)

## ─────────────────────────────────────────────────────────────────────────────
## mrgsolve model code
## ─────────────────────────────────────────────────────────────────────────────

code <- '
$PROB
JIA QSP Model — PK/PD for pediatric inflammatory arthritis
Compartments (21 ODE): MTX PK (4) | Biologic PK ETN (3) | TCZ (3) |
                       Canakinumab PK (2) | GC (1) | Cytokines (4) |
                       Biomarkers (2) | Tissue (2)

$PARAM @annotated
// ---- MTX PK (weekly oral/SC dosing) ----
ka_mtx   : 0.50   : MTX GI absorption rate constant (1/h)
F_mtx    : 0.72   : MTX oral bioavailability (fraction)
V1_mtx   : 11.0   : MTX central volume (L)
V2_mtx   : 28.0   : MTX peripheral volume (L)
CL_mtx   : 4.2    : MTX total clearance (L/h) [renal dominant]
Q_mtx    : 7.5    : MTX inter-compartmental clearance (L/h)
kpg      : 0.007  : MTX polyglutamation rate constant (1/h per mg/L)
kdepg    : 0.0015 : MTX polyglutamate de-polyglutamation (1/h)

// ---- Etanercept PK (biweekly/weekly SC 25-50 mg) ----
ka_eta   : 0.019  : Etanercept SC absorption rate (1/h; t_peak ~48h)
F_eta    : 0.76   : Etanercept SC bioavailability
V1_eta   : 6.5    : Etanercept central volume (L)
V2_eta   : 4.1    : Etanercept peripheral volume (L)
CL_eta   : 0.060  : Etanercept linear clearance (L/h; t½~70h)
Q_eta    : 0.090  : Etanercept inter-compartmental CL (L/h)
Vmax_eta : 0.045  : Etanercept saturable clearance Vmax (mg/h, TMDD)
Km_eta   : 2.5    : Etanercept TMDD Michaelis constant (mg/L)

// ---- Tocilizumab PK (q2w or q4w SC/IV) ----
ka_tcz   : 0.018  : Tocilizumab SC absorption rate (1/h)
F_tcz    : 0.80   : Tocilizumab SC bioavailability
V1_tcz   : 4.8    : Tocilizumab central volume (L)
V2_tcz   : 2.2    : Tocilizumab peripheral volume (L)
CL_tcz   : 0.021  : Tocilizumab linear clearance (L/h; t½~11d)
Q_tcz    : 0.035  : Tocilizumab inter-compartmental CL (L/h)
Vmax_tcz : 0.048  : Tocilizumab saturable CL Vmax (IL-6R-mediated)
Km_tcz   : 1.8    : Tocilizumab Km for saturable CL (mg/L)

// ---- Canakinumab PK (q4w SC 4 mg/kg) ----
ka_can   : 0.0138 : Canakinumab SC absorption rate (1/h; t½ SC ~7d)
F_can    : 0.70   : Canakinumab SC bioavailability
V1_can   : 6.8    : Canakinumab central volume (L)
CL_can   : 0.0065 : Canakinumab total clearance (L/h; t½~26d)

// ---- Prednisolone PK (oral, variable dose) ----
ka_gc    : 1.80   : Prednisolone oral absorption rate (1/h)
V_gc     : 42.0   : Prednisolone volume of distribution (L)
CL_gc    : 15.0   : Prednisolone clearance (L/h; t½~3h)

// ---- Baricitinib PK (daily oral, pediatric weight-based) ----
ka_bar   : 0.95   : Baricitinib absorption rate (1/h)
V_bar    : 76.0   : Baricitinib Vd (L)
CL_bar   : 6.4    : Baricitinib clearance (L/h; t½~12h)

// ---- PD: Cytokine baseline production rates ----
kprod_tnf : 6.0   : TNF-α basal production (pg/mL/h)
kdeg_tnf  : 0.55  : TNF-α degradation rate (1/h; t½~1.3h)
kprod_il6 : 3.5   : IL-6 basal production (pg/mL/h)
kdeg_il6  : 0.32  : IL-6 degradation rate (1/h; t½~2.2h)
kprod_il1 : 2.5   : IL-1β basal production (pg/mL/h)
kdeg_il1  : 0.45  : IL-1β degradation rate (1/h; t½~1.5h)
kprod_il18: 1.8   : IL-18 basal production (pg/mL/h; sJIA)
kdeg_il18 : 0.22  : IL-18 degradation rate (1/h; t½~3.2h)

// ---- PD: Cytokine cross-amplification ----
kTNF_IL6  : 0.055 : TNF-α stimulation of IL-6 (pg/mL normalized)
kIL1_TNF  : 0.048 : IL-1β stimulation of TNF-α
kIL6_auto : 0.030 : IL-6 positive feedback on IL-6 production
kTNF_feed : 0.085 : TNF-α autocrine NF-kB amplification
kIL1_IL18 : 0.060 : IL-1β stimulation of IL-18 (NLRP3 loop)

// ---- PD: Biomarker production ----
kprod_crp : 0.080 : CRP production rate constant (per IL-6 pg/mL/h)
kdeg_crp  : 0.028 : CRP elimination rate (1/h; t½~24h)
kprod_esr : 0.040 : ESR increase rate (per combined cytokines/h)
kdeg_esr  : 0.008 : ESR normalization rate (1/h; t½~3.6d)

// ---- PD: Joint damage ----
kprog_cart: 0.0012 : Cartilage degradation rate (1/h per TNF+IL1 score)
krep_cart : 0.0020 : Cartilage repair potential (1/h)
kprog_bone: 0.0006 : Periarticular bone erosion rate (1/h)
kloss_bmd : 0.0004 : BMD loss per GC unit (1/h per mg/L GC)
krec_bmd  : 0.0002 : BMD recovery rate (1/h)

// ---- PD Drug Emax models ----
Emax_eta  : 0.87   : Etanercept maximal TNF inhibition
EC50_eta  : 0.50   : Etanercept EC50 (mg/L)
Emax_tcz  : 0.91   : Tocilizumab maximal IL-6 signaling inhibition
EC50_tcz  : 0.60   : Tocilizumab EC50 (mg/L)
Emax_can  : 0.88   : Canakinumab maximal IL-1β inhibition
EC50_can  : 0.90   : Canakinumab EC50 (mg/L)
Emax_gc   : 0.68   : GC maximal anti-inflammatory effect
EC50_gc   : 0.12   : GC EC50 (mg/L)
Emax_mtx  : 0.52   : MTX polyglut. maximal anti-inflam
EC50_mtx  : 0.006  : MTX EC50 (umol/L poly)
Emax_bar  : 0.78   : Baricitinib maximal JAK inhibition
EC50_bar  : 0.042  : Baricitinib EC50 (mg/L; IC50 JAK1 ~5.9 nM)
Hill_bar  : 1.3    : Baricitinib Hill coefficient

// ---- AJC and JADAS baseline ----
AJC_base  : 13.0   : Untreated active joint count (polyarticular ~13)
kAJC_IL6  : 140.0  : IL-6 half-effect concentration for AJC (pg/mL)
kAJC_TNF  : 60.0   : TNF-α half-effect for AJC contribution
JADAS_base: 20.0   : Untreated JADAS-27 (high activity)

$CMT @annotated
// MTX compartments
MTX_GI     : MTX gut lumen (mg)
MTX_C      : MTX central plasma (mg)
MTX_P      : MTX peripheral tissue (mg)
MTX_Poly   : MTX polyglutamates intracellular (umol)

// Etanercept compartments
ETA_SC     : Etanercept SC depot (mg)
ETA_C      : Etanercept central plasma (mg)
ETA_P      : Etanercept peripheral (mg)

// Tocilizumab compartments
TCZ_SC     : Tocilizumab SC depot (mg)
TCZ_C      : Tocilizumab central (mg)
TCZ_P      : Tocilizumab peripheral (mg)

// Canakinumab compartments
CAN_SC     : Canakinumab SC depot (mg)
CAN_C      : Canakinumab central (mg)

// Glucocorticoid
GC_C       : Prednisolone central (mg)

// Baricitinib
BAR_C      : Baricitinib central (mg)

// Inflammatory mediators
TNF        : TNF-alpha synovial/serum (pg/mL)
IL6        : IL-6 synovial/serum (pg/mL)
IL1        : IL-1beta synovial (pg/mL)
IL18       : IL-18 systemic (pg/mL) [sJIA]

// Biomarkers
CRP        : C-reactive protein (mg/L)
ESR        : ESR (mm/h)

// Tissue
Cartilage  : Articular cartilage integrity (% normal, 0-100)
BMD        : Periarticular bone mineral density (% normal)

$MAIN
// ── Plasma concentrations ──────────────────────────────────────────────────
double C_MTX  = MTX_C  / V1_mtx;    // mg/L
double C_ETA  = ETA_C  / V1_eta;    // mg/L
double C_TCZ  = TCZ_C  / V1_tcz;    // mg/L
double C_CAN  = CAN_C  / V1_can;    // mg/L
double C_GC   = GC_C   / V_gc;      // mg/L
double C_BAR  = BAR_C  / V_bar;     // mg/L
double C_poly = MTX_Poly;            // umol/L (intracellular proxy)

// ── Drug inhibition fractions (Emax / Hill) ───────────────────────────────
double Ieta  = Emax_eta * C_ETA  / (EC50_eta  + C_ETA);
double Itcz  = Emax_tcz * C_TCZ  / (EC50_tcz  + C_TCZ);
double Ican  = Emax_can * C_CAN  / (EC50_can  + C_CAN);
double Igc   = Emax_gc  * C_GC   / (EC50_gc   + C_GC);
double Imtx  = Emax_mtx * C_poly / (EC50_mtx  + C_poly);
double Ibar  = Emax_bar * pow(C_BAR, Hill_bar) /
               (pow(EC50_bar, Hill_bar) + pow(C_BAR, Hill_bar));

// Combined small-molecule effect (GC + MTX + baricitinib)
// Additive inhibition on cytokine production (simplified Bliss independence)
double I_sm  = 1.0 - (1.0 - Igc) * (1.0 - Imtx) * (1.0 - Ibar * 0.5);

// ── Active Joint Count (derived) ──────────────────────────────────────────
double AJC = AJC_base * (IL6 / (kAJC_IL6 + IL6) * 0.6 +
                         TNF / (kAJC_TNF + TNF) * 0.4);
double AJC_norm = AJC / AJC_base;  // 0-1

// ── JADAS-27 composite (simplified continuous approximation) ──────────────
// JADAS-27 = AJC (0-27) + CRP-normalized (0-10) + PGA (0-10) + PtGA (0-10)
// Here modeled as function of AJC and CRP with assumed PGA scaling
double JADAS = AJC * (27.0 / AJC_base) +
               log1p(CRP) / log1p(200.0) * 10.0 +
               5.0 * AJC_norm;   // PGA proxy

$ODE
// ═══════════════════════════════════════════════════════════════════
// MTX PHARMACOKINETICS
// ═══════════════════════════════════════════════════════════════════
dxdt_MTX_GI   = -ka_mtx * MTX_GI;

dxdt_MTX_C    =  F_mtx * ka_mtx * MTX_GI
                - (CL_mtx / V1_mtx + Q_mtx / V1_mtx) * MTX_C
                + (Q_mtx / V2_mtx) * MTX_P;

dxdt_MTX_P    =  (Q_mtx / V1_mtx) * MTX_C
                - (Q_mtx / V2_mtx) * MTX_P;

// MTX polyglutamation (slow intracellular retention)
dxdt_MTX_Poly =  kpg * C_MTX - kdepg * MTX_Poly;

// ═══════════════════════════════════════════════════════════════════
// ETANERCEPT PHARMACOKINETICS
// ═══════════════════════════════════════════════════════════════════
dxdt_ETA_SC   = -ka_eta * ETA_SC;

dxdt_ETA_C    =  F_eta * ka_eta * ETA_SC
                - (CL_eta / V1_eta + Q_eta / V1_eta) * ETA_C
                + (Q_eta / V2_eta) * ETA_P
                - Vmax_eta * C_ETA / (Km_eta + C_ETA);  // TMDD

dxdt_ETA_P    =  (Q_eta / V1_eta) * ETA_C
                - (Q_eta / V2_eta) * ETA_P;

// ═══════════════════════════════════════════════════════════════════
// TOCILIZUMAB PHARMACOKINETICS
// ═══════════════════════════════════════════════════════════════════
dxdt_TCZ_SC   = -ka_tcz * TCZ_SC;

dxdt_TCZ_C    =  F_tcz * ka_tcz * TCZ_SC
                - (CL_tcz / V1_tcz + Q_tcz / V1_tcz) * TCZ_C
                + (Q_tcz / V2_tcz) * TCZ_P
                - Vmax_tcz * C_TCZ / (Km_tcz + C_TCZ);  // IL-6R-mediated CL

dxdt_TCZ_P    =  (Q_tcz / V1_tcz) * TCZ_C
                - (Q_tcz / V2_tcz) * TCZ_P;

// ═══════════════════════════════════════════════════════════════════
// CANAKINUMAB PHARMACOKINETICS
// ═══════════════════════════════════════════════════════════════════
dxdt_CAN_SC   = -ka_can * CAN_SC;

dxdt_CAN_C    =  F_can * ka_can * CAN_SC
                - (CL_can / V1_can) * CAN_C;

// ═══════════════════════════════════════════════════════════════════
// PREDNISOLONE PHARMACOKINETICS
// ═══════════════════════════════════════════════════════════════════
dxdt_GC_C     = -(CL_gc / V_gc) * GC_C;  // IV or oral bolus events

// ═══════════════════════════════════════════════════════════════════
// BARICITINIB PHARMACOKINETICS
// ═══════════════════════════════════════════════════════════════════
dxdt_BAR_C    = -(CL_bar / V_bar) * BAR_C;

// ═══════════════════════════════════════════════════════════════════
// INFLAMMATORY MEDIATOR DYNAMICS
// ═══════════════════════════════════════════════════════════════════

// TNF-α: produced by macrophages, amplified by IL-1 and NF-kB;
//        inhibited by Etanercept (primary) and small molecules
dxdt_TNF = kprod_tnf
           * (1.0 + kTNF_feed * TNF / (100.0 + TNF))
           * (1.0 + kIL1_TNF  * IL1 / (15.0  + IL1))
           * (1.0 - Ieta)
           * (1.0 - I_sm)
           - kdeg_tnf * TNF;

// IL-6: stimulated by TNF-α, self-amplified; inhibited by Tocilizumab (IL-6R)
//       and small molecules
dxdt_IL6 = kprod_il6
           * (1.0 + kTNF_IL6  * TNF  / (50.0 + TNF))
           * (1.0 + kIL6_auto * IL6  / (20.0 + IL6))
           * (1.0 - Itcz)
           * (1.0 - I_sm)
           - kdeg_il6 * IL6;

// IL-1β: NLRP3 inflammasome-driven; inhibited by Canakinumab and small molecules
dxdt_IL1 = kprod_il1
           * (1.0 + 0.04 * TNF / (50.0 + TNF))
           * (1.0 - Ican)
           * (1.0 - I_sm)
           - kdeg_il1 * IL1;

// IL-18: key in sJIA/MAS; amplified by IL-1β via NLRP3
dxdt_IL18 = kprod_il18
            * (1.0 + kIL1_IL18 * IL1 / (10.0 + IL1))
            * (1.0 - Ican * 0.6)    // partial canakinumab effect on IL-18
            * (1.0 - I_sm * 0.4)
            - kdeg_il18 * IL18;

// ═══════════════════════════════════════════════════════════════════
// BIOMARKER DYNAMICS
// ═══════════════════════════════════════════════════════════════════

// CRP: produced by hepatocytes in response to IL-6 (primary) and IL-1
dxdt_CRP = kprod_crp * (IL6 + 0.4 * IL1)
           - kdeg_crp * CRP;

// ESR: driven by fibrinogen (IL-6→fibrinogen), slow kinetics
dxdt_ESR = kprod_esr * (IL6 + 0.3 * IL1 + 0.1 * IL18)
           - kdeg_esr * ESR;

// ═══════════════════════════════════════════════════════════════════
// JOINT DAMAGE (irreversible component)
// ═══════════════════════════════════════════════════════════════════

// Cartilage: degraded by TNF/IL-1 via MMP/ADAMTS; partial repair possible
dxdt_Cartilage = -kprog_cart * (TNF / 50.0 + IL1 / 15.0) * Cartilage
                 + krep_cart * (100.0 - Cartilage)
                   * (1.0 - 0.3 * (1.0 - Ieta))  // TNFi helps repair
                 * (Cartilage > 5.0 ? 1.0 : 0.0); // floor at 5%

// Periarticular BMD: eroded by RANKL-driven osteoclasts (TNF/IL-17/IL-1)
//                    and GC; partial recovery with disease control
dxdt_BMD = -kprog_bone * (TNF / 50.0 + IL1 / 15.0) * BMD / 100.0
           - kloss_bmd  * C_GC
           + krec_bmd   * (100.0 - BMD) * Ieta; // recovery when TNFi active

$TABLE
// Plasma concentrations (mg/L)
double Cp_MTX  = MTX_C  / V1_mtx;
double Cp_ETA  = ETA_C  / V1_eta;
double Cp_TCZ  = TCZ_C  / V1_tcz;
double Cp_CAN  = CAN_C  / V1_can;
double Cp_GC   = GC_C   / V_gc;
double Cp_BAR  = BAR_C  / V_bar;

// MTX polyglutamates (umol/L proxy)
double Cp_poly = MTX_Poly;

// Drug inhibition fractions (%)
double pct_Ieta  = 100.0 * Emax_eta * Cp_ETA / (EC50_eta + Cp_ETA);
double pct_Itcz  = 100.0 * Emax_tcz * Cp_TCZ / (EC50_tcz + Cp_TCZ);
double pct_Ican  = 100.0 * Emax_can * Cp_CAN / (EC50_can + Cp_CAN);
double pct_Igc   = 100.0 * Emax_gc  * Cp_GC  / (EC50_gc  + Cp_GC);
double pct_Imtx  = 100.0 * Emax_mtx * Cp_poly/ (EC50_mtx + Cp_poly);
double pct_Ibar  = 100.0 * Emax_bar * pow(Cp_BAR, Hill_bar) /
                  (pow(EC50_bar, Hill_bar) + pow(Cp_BAR, Hill_bar));

// Clinical outcome measures
double AJC_out   = AJC_base * (IL6 / (kAJC_IL6 + IL6) * 0.6 +
                               TNF / (kAJC_TNF + TNF) * 0.4);
double JADAS_out = AJC_out * (27.0 / AJC_base) +
                  log1p(CRP) / log1p(200.0) * 10.0 +
                  5.0 * AJC_out / AJC_base;

// % improvement from untreated baseline
double pct_AJC_improve   = 100.0 * (1.0 - AJC_out / AJC_base);
double pct_JADAS_improve = 100.0 * (1.0 - JADAS_out / JADAS_base);

// ACR Pediatric criteria (binary flags)
double ACR30 = (pct_AJC_improve >= 30.0) ? 1.0 : 0.0;
double ACR50 = (pct_AJC_improve >= 50.0) ? 1.0 : 0.0;
double ACR70 = (pct_AJC_improve >= 70.0) ? 1.0 : 0.0;
double ACR90 = (pct_AJC_improve >= 90.0) ? 1.0 : 0.0;

// Joint space narrowing (proxy for cartilage damage)
double JSN_pct = 100.0 - Cartilage;

// Remission (Wallace criteria: JADAS ≤ 1.0 for ≥6 months)
double Remission = (JADAS_out <= 1.0) ? 1.0 : 0.0;

// CID (clinically inactive disease) proxy
double CID = (AJC_out < 0.5 && CRP < 1.0) ? 1.0 : 0.0;

$CAPTURE Cp_MTX Cp_ETA Cp_TCZ Cp_CAN Cp_GC Cp_BAR Cp_poly
         pct_Ieta pct_Itcz pct_Ican pct_Igc pct_Imtx pct_Ibar
         TNF IL6 IL1 IL18 CRP ESR
         AJC_out JADAS_out pct_AJC_improve pct_JADAS_improve
         ACR30 ACR50 ACR70 ACR90
         Cartilage BMD JSN_pct Remission CID
'

## ─────────────────────────────────────────────────────────────────────────────
## Compile model
## ─────────────────────────────────────────────────────────────────────────────

mod <- mcode("JIA_QSP", code)

## ─────────────────────────────────────────────────────────────────────────────
## Initial conditions (steady-state untreated polyarticular JIA)
## ─────────────────────────────────────────────────────────────────────────────

# Cytokine steady states (untreated, high disease activity)
# TNF_ss: kprod_tnf * (1 + feed) / kdeg_tnf ≈ 6*2/0.55 ≈ 22 pg/mL
# IL6_ss: ≈ 3.5*(1+0.055*22/72)*(1+self)/0.32 ≈ 12 pg/mL
# IL1_ss: ≈ 2.5*(1+0.04*22/72)/0.45 ≈ 6 pg/mL
# CRP_ss: 0.08*(12+0.4*6)/0.028 ≈ 40 mg/L
# ESR_ss: 0.04*(12+0.3*6)/0.008 ≈ 69 mm/h

init_vals <- c(
  MTX_GI=0, MTX_C=0, MTX_P=0, MTX_Poly=0,
  ETA_SC=0, ETA_C=0, ETA_P=0,
  TCZ_SC=0, TCZ_C=0, TCZ_P=0,
  CAN_SC=0, CAN_C=0,
  GC_C=0, BAR_C=0,
  TNF=22, IL6=12, IL1=6, IL18=25,
  CRP=40, ESR=65,
  Cartilage=92, BMD=95
)

mod <- mod %>% init(init_vals)

## ─────────────────────────────────────────────────────────────────────────────
## Helper: build event dosing tables
## ─────────────────────────────────────────────────────────────────────────────

make_events <- function(type="mtx", dose=15, start=0, end=52*7, interval=168) {
  # interval in hours
  times <- seq(start, end, by=interval)
  cmt_map <- c(mtx="MTX_GI", etanercept="ETA_SC", tocilizumab="TCZ_SC",
               canakinumab="CAN_SC", prednisolone="GC_C", baricitinib="BAR_C")
  cmt <- cmt_map[type]
  ev(amt=dose, cmt=cmt, ii=interval, addl=length(times)-1, time=times[1])
}

## ─────────────────────────────────────────────────────────────────────────────
## Scenario 1: Natural history (no treatment)
## ─────────────────────────────────────────────────────────────────────────────

sim_notrx <- mod %>%
  mrgsim(end=52*7, delta=24) %>%
  as.data.frame() %>%
  mutate(scenario="No Treatment", week=time/168)

cat("\n=== Scenario 1: Natural history ===\n")
cat("Week 24 JADAS-27:", round(sim_notrx$JADAS_out[sim_notrx$week==24][1], 1), "\n")
cat("Week 24 AJC:    ", round(sim_notrx$AJC_out[sim_notrx$week==24][1], 1), "\n")

## ─────────────────────────────────────────────────────────────────────────────
## Scenario 2: MTX Monotherapy (15 mg/week, oral)
## ─────────────────────────────────────────────────────────────────────────────

ev_mtx <- ev(cmt="MTX_GI", amt=15, ii=168, addl=51, time=0)
sim_mtx <- mod %>%
  ev(ev_mtx) %>%
  mrgsim(end=52*7, delta=24) %>%
  as.data.frame() %>%
  mutate(scenario="MTX Mono", week=time/168)

cat("\n=== Scenario 2: MTX 15 mg/week ===\n")
cat("Week 24 ACR30:", round(mean(sim_mtx$ACR30[sim_mtx$week==24]), 2), "\n")
cat("Week 24 JADAS:", round(sim_mtx$JADAS_out[sim_mtx$week==24][1], 1), "\n")

## ─────────────────────────────────────────────────────────────────────────────
## Scenario 3: MTX + Etanercept (25 mg biweekly SC)
## ─────────────────────────────────────────────────────────────────────────────

ev_eta <- ev(cmt="ETA_SC", amt=25, ii=336, addl=25, time=0)  # biweekly
ev_combo1 <- c(ev_mtx, ev_eta)

sim_combo1 <- mod %>%
  ev(ev_combo1) %>%
  mrgsim(end=52*7, delta=24) %>%
  as.data.frame() %>%
  mutate(scenario="MTX + Etanercept", week=time/168)

cat("\n=== Scenario 3: MTX + Etanercept ===\n")
w24 <- sim_combo1$week == 24
cat("Week 24 ACR30:", round(mean(sim_combo1$ACR30[w24]), 2), "\n")
cat("Week 24 ACR70:", round(mean(sim_combo1$ACR70[w24]), 2), "\n")
cat("Week 24 JADAS:", round(sim_combo1$JADAS_out[w24][1], 1), "\n")

## ─────────────────────────────────────────────────────────────────────────────
## Scenario 4: Tocilizumab (sJIA: 8 mg/kg IV q2w or 162 mg SC q2w)
## ─────────────────────────────────────────────────────────────────────────────

# SC 162 mg every 2 weeks for sJIA
ev_tcz <- ev(cmt="TCZ_SC", amt=162, ii=336, addl=25, time=0)

sim_tcz <- mod %>%
  ev(ev_tcz) %>%
  mrgsim(end=52*7, delta=24) %>%
  as.data.frame() %>%
  mutate(scenario="Tocilizumab (sJIA)", week=time/168)

cat("\n=== Scenario 4: Tocilizumab (sJIA subtype) ===\n")
w24 <- sim_tcz$week == 24
cat("Week 24 IL-6:", round(sim_tcz$IL6[w24][1], 1), "pg/mL\n")
cat("Week 24 CRP: ", round(sim_tcz$CRP[w24][1], 1), "mg/L\n")
cat("Week 24 JADAS:", round(sim_tcz$JADAS_out[w24][1], 1), "\n")

## ─────────────────────────────────────────────────────────────────────────────
## Scenario 5: Canakinumab (sJIA: 4 mg/kg q4w SC, max 300 mg)
## ─────────────────────────────────────────────────────────────────────────────

# 150 mg q4w (representative 25-30 kg child at ~5 mg/kg → 150 mg)
ev_can <- ev(cmt="CAN_SC", amt=150, ii=672, addl=12, time=0)  # q4w

sim_can <- mod %>%
  ev(ev_can) %>%
  mrgsim(end=52*7, delta=24) %>%
  as.data.frame() %>%
  mutate(scenario="Canakinumab (sJIA)", week=time/168)

cat("\n=== Scenario 5: Canakinumab (sJIA) ===\n")
w24 <- sim_can$week == 24
cat("Week 24 IL-1β:", round(sim_can$IL1[w24][1], 1), "pg/mL\n")
cat("Week 24 IL-18:", round(sim_can$IL18[w24][1], 1), "pg/mL\n")
cat("Week 24 Remission rate:", round(mean(sim_can$Remission[w24]), 2), "\n")

## ─────────────────────────────────────────────────────────────────────────────
## Scenario 6: MTX + Prednisolone bridge (short course) + step up to Etanercept
## ─────────────────────────────────────────────────────────────────────────────

# Step-up strategy: MTX + GC bridge (wk 0-12) then add ETN (wk 12+)
ev_gc_bridge <- ev(cmt="GC_C", amt=10, ii=24, addl=11*7-1, time=0) # 0.3 mg/kg for 12wk
ev_eta_late  <- ev(cmt="ETA_SC", amt=25, ii=336, addl=19, time=12*7*24)
ev_stepup    <- c(ev_mtx, ev_gc_bridge, ev_eta_late)

sim_stepup <- mod %>%
  ev(ev_stepup) %>%
  mrgsim(end=52*7, delta=24) %>%
  as.data.frame() %>%
  mutate(scenario="Step-up Strategy\n(MTX→GC→+ETN)", week=time/168)

cat("\n=== Scenario 6: Step-up (MTX + GC bridge → ETN) ===\n")
w52 <- sim_stepup$week == 52
cat("Week 52 JADAS:", round(sim_stepup$JADAS_out[w52][1], 1), "\n")
cat("Week 52 BMD:  ", round(sim_stepup$BMD[w52][1], 1), "% normal\n")

## ─────────────────────────────────────────────────────────────────────────────
## Scenario 7: Baricitinib (JAK1/2, 4 mg daily; age ≥2 yr, ≥10 kg)
## ─────────────────────────────────────────────────────────────────────────────

ev_bar <- ev(cmt="BAR_C", amt=4, ii=24, addl=52*7-1, time=0)

sim_bar <- mod %>%
  ev(ev_bar) %>%
  mrgsim(end=52*7, delta=24) %>%
  as.data.frame() %>%
  mutate(scenario="Baricitinib", week=time/168)

cat("\n=== Scenario 7: Baricitinib 4 mg QD ===\n")
w12 <- sim_bar$week == 12
cat("Week 12 AJC:  ", round(sim_bar$AJC_out[w12][1], 1), "\n")
cat("Week 12 JADAS:", round(sim_bar$JADAS_out[w12][1], 1), "\n")

## ─────────────────────────────────────────────────────────────────────────────
## Combine all scenarios
## ─────────────────────────────────────────────────────────────────────────────

all_sim <- bind_rows(
  sim_notrx, sim_mtx, sim_combo1, sim_tcz, sim_can, sim_stepup, sim_bar
)

## ─────────────────────────────────────────────────────────────────────────────
## Visualization
## ─────────────────────────────────────────────────────────────────────────────

colors7 <- c(
  "No Treatment"       = "#E53935",
  "MTX Mono"           = "#FB8C00",
  "MTX + Etanercept"   = "#43A047",
  "Tocilizumab (sJIA)" = "#1E88E5",
  "Canakinumab (sJIA)" = "#8E24AA",
  "Step-up Strategy\n(MTX→GC→+ETN)" = "#00ACC1",
  "Baricitinib"        = "#6D4C41"
)

p1 <- all_sim %>%
  filter(week <= 52) %>%
  ggplot(aes(week, JADAS_out, color=scenario)) +
  geom_line(linewidth=1.0) +
  geom_hline(yintercept=1.0, linetype="dashed", color="grey40",
             linewidth=0.7) +
  annotate("text", x=2, y=2.2, label="Remission threshold\n(JADAS ≤ 1)", size=2.8) +
  labs(title="JADAS-27 Over Time", x="Week", y="JADAS-27",
       color="Treatment") +
  scale_color_manual(values=colors7) +
  theme_bw(base_size=11) +
  theme(legend.position="right", legend.text=element_text(size=8))

p2 <- all_sim %>%
  filter(week <= 52) %>%
  ggplot(aes(week, AJC_out, color=scenario)) +
  geom_line(linewidth=1.0) +
  labs(title="Active Joint Count", x="Week", y="AJC (joints)",
       color="Treatment") +
  scale_color_manual(values=colors7) +
  theme_bw(base_size=11) +
  theme(legend.position="none")

p3 <- all_sim %>%
  filter(week <= 52) %>%
  ggplot(aes(week, CRP, color=scenario)) +
  geom_line(linewidth=1.0) +
  geom_hline(yintercept=5, linetype="dashed", color="grey40") +
  labs(title="CRP (mg/L)", x="Week", y="CRP (mg/L)",
       color="Treatment") +
  scale_color_manual(values=colors7) +
  theme_bw(base_size=11) +
  theme(legend.position="none")

p4 <- all_sim %>%
  filter(week <= 52) %>%
  ggplot(aes(week, IL6, color=scenario)) +
  geom_line(linewidth=1.0) +
  labs(title="IL-6 (pg/mL)", x="Week", y="IL-6 (pg/mL)",
       color="Treatment") +
  scale_color_manual(values=colors7) +
  theme_bw(base_size=11) +
  theme(legend.position="none")

p5 <- all_sim %>%
  filter(week <= 52) %>%
  ggplot(aes(week, Cartilage, color=scenario)) +
  geom_line(linewidth=1.0) +
  labs(title="Articular Cartilage Integrity (%)",
       x="Week", y="Cartilage Integrity (%)", color="Treatment") +
  scale_color_manual(values=colors7) +
  scale_y_continuous(limits=c(0, 100)) +
  theme_bw(base_size=11) +
  theme(legend.position="none")

p6 <- all_sim %>%
  filter(week <= 52) %>%
  ggplot(aes(week, IL18, color=scenario)) +
  geom_line(linewidth=1.0) +
  labs(title="IL-18 (pg/mL) — sJIA/MAS marker",
       x="Week", y="IL-18 (pg/mL)", color="Treatment") +
  scale_color_manual(values=colors7) +
  theme_bw(base_size=11) +
  theme(legend.position="none")

# Combined figure
(p1 | p2) / (p3 | p4) / (p5 | p6) +
  plot_annotation(
    title="JIA QSP Model — Treatment Scenario Comparison",
    subtitle="Juvenile Idiopathic Arthritis (Polyarticular + sJIA Subtypes)",
    theme=theme(plot.title=element_text(size=14, face="bold"),
                plot.subtitle=element_text(size=11))
  )

## ─────────────────────────────────────────────────────────────────────────────
## ACR Pedi response summary at Week 24 (key regulatory endpoint)
## ─────────────────────────────────────────────────────────────────────────────

week24_summary <- all_sim %>%
  filter(abs(week - 24) < 0.5) %>%
  group_by(scenario) %>%
  summarise(
    JADAS27     = round(mean(JADAS_out), 1),
    AJC         = round(mean(AJC_out), 1),
    CRP         = round(mean(CRP), 1),
    IL6         = round(mean(IL6), 1),
    ACR30_pct   = round(mean(ACR30) * 100, 0),
    ACR50_pct   = round(mean(ACR50) * 100, 0),
    ACR70_pct   = round(mean(ACR70) * 100, 0),
    Cartilage   = round(mean(Cartilage), 1),
    Remission   = round(mean(Remission) * 100, 0)
  )

print(week24_summary)

## ─────────────────────────────────────────────────────────────────────────────
## PK profile: Single dose etanercept (25 mg SC)
## ─────────────────────────────────────────────────────────────────────────────

mod_pk <- mod %>% init(
  MTX_GI=0, MTX_C=0, MTX_P=0, MTX_Poly=0,
  ETA_SC=0, ETA_C=0, ETA_P=0,
  TCZ_SC=0, TCZ_C=0, TCZ_P=0,
  CAN_SC=0, CAN_C=0,
  GC_C=0, BAR_C=0,
  TNF=0, IL6=0, IL1=0, IL18=0,
  CRP=0, ESR=0, Cartilage=100, BMD=100
)

pk_single <- mod_pk %>%
  ev(ev(cmt="ETA_SC", amt=25, time=0)) %>%
  mrgsim(end=336, delta=2) %>%
  as.data.frame()

ggplot(pk_single, aes(time/24, Cp_ETA)) +
  geom_line(color="#1565C0", linewidth=1.2) +
  labs(title="Etanercept PK — Single SC Dose (25 mg)",
       x="Time (days)", y="Etanercept Cp (mg/L)") +
  theme_bw(base_size=12)
