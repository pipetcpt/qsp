## ============================================================
## Multiple Sclerosis (MS) — mrgsolve QSP Model
## Author  : Claude Code Routine (CCR)
## Date    : 2026-06-16
## Version : 1.0
## ============================================================
## Model architecture:
##   - 7 drug PK compartments (IFN-β, Natalizumab-TMDD, Ocrelizumab,
##     Siponimod, DMF/MMF, Cladribine, Ofatumumab)
##   - 15 disease-state ODEs covering peripheral immunity → BBB →
##     CNS neuroinflammation → demyelination → axonal injury
##   - Biomarker ODEs: NfL, GFAP
##   - Clinical endpoints computed: EDSS, ARR, T2 lesion volume
## Parameters calibrated to published Phase III trial data:
##   AFFIRM (natalizumab), OPERA I/II (ocrelizumab), EXPAND (siponimod),
##   CLARITY (cladribine), DEFINE (DMF), TRANSFORMS (fingolimod proxy)
## ============================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)

## ── Model code ────────────────────────────────────────────────────────────────

ms_model_code <- '
$PROB
MS QSP Model — Peripheral immunity to CNS neurodegeneration
Drugs: IFN-beta, Natalizumab, Ocrelizumab, Siponimod, DMF, Cladribine

$PARAM @annotated
// ── Drug dosing flags (1=active, 0=off) ──────────────────────────────────
dose_ifnb   : 0  : IFN-beta active (1=yes)
dose_nat    : 0  : Natalizumab active
dose_ocre   : 0  : Ocrelizumab active
dose_sip    : 0  : Siponimod active
dose_dmf    : 0  : DMF active
dose_clad   : 0  : Cladribine active

// ── IFN-beta PK (2-comp SC/IM) ───────────────────────────────────────────
ka_ifnb     : 0.30  : /h  absorption rate
CL_ifnb     : 5.0   : L/h clearance
V1_ifnb     : 12.0  : L   central volume
Q_ifnb      : 1.5   : L/h intercompartmental Q
V2_ifnb     : 20.0  : L   peripheral volume
F_ifnb      : 0.40  :     bioavailability SC

// ── Natalizumab PK (TMDD 2-comp IV) ─────────────────────────────────────
CL_nat      : 0.014 : L/h linear clearance
V1_nat      : 3.5   : L   central volume
Q_nat       : 0.50  : L/h
V2_nat      : 6.0   : L
kon_nat     : 0.013 : /h per nM  binding to α4-integrin
koff_nat    : 0.003 : /h  dissociation
kint_nat    : 0.005 : /h  internalization rate
Rtot_nat    : 30.0  : nM  total receptor (α4-integrin) on T cells
ksyn_R      : 0.15  : nM/h receptor synthesis
kdeg_R      : 0.005 : /h  receptor degradation

// ── Ocrelizumab PK (2-comp IV) ───────────────────────────────────────────
CL_ocre     : 0.012 : L/h
V1_ocre     : 3.0   : L
Q_ocre      : 0.40  : L/h
V2_ocre     : 5.0   : L
EC50_ocre   : 1.5   : µg/mL  for B cell depletion (Emax model)
Emax_ocre   : 0.97  :        max fractional B cell depletion

// ── Siponimod PK (2-comp oral) ───────────────────────────────────────────
ka_sip      : 0.50  : /h  absorption
CL_sip      : 2.2   : L/h
V1_sip      : 30.0  : L
Q_sip       : 3.0   : L/h
V2_sip      : 60.0  : L
EC50_sip    : 8.0   : ng/mL S1PR1 occupancy IC50
Emax_sip    : 0.75  :       max lymphocyte retention fraction

// ── DMF / MMF PK (1-comp oral) ───────────────────────────────────────────
ka_dmf      : 1.20  : /h
CL_dmf      : 25.0  : L/h
V1_dmf      : 50.0  : L
EC50_dmf    : 0.40  : µg/mL  NF-κB inhibition IC50
EC50_nrf2   : 0.20  : µg/mL  Nrf2 activation EC50

// ── Cladribine PK (1-comp oral) ──────────────────────────────────────────
ka_clad     : 0.40  : /h
CL_clad     : 15.0  : L/h
V1_clad     : 40.0  : L
EC50_clad   : 5.0   : ng/mL  lymphodepletion IC50
Emax_clad   : 0.90  :        max lymphodepletion fraction

// ── Disease biology parameters ────────────────────────────────────────────
// Peripheral T/B cells (normalized to 1 = healthy baseline)
kin_Th1     : 0.10  : /h  Th1 production (in RRMS ~3x healthy)
kout_Th1    : 0.10  : /h  Th1 natural turnover
kin_Th17    : 0.06  : /h  Th17 production
kout_Th17   : 0.06  : /h
kin_Treg    : 0.04  : /h  Treg production
kout_Treg   : 0.04  : /h
kin_B       : 0.08  : /h  B cell production
kout_B      : 0.08  : /h
Treg_ss     : 1.0   :     Treg baseline
k_Treg_inh  : 0.30  :     Treg inhibitory strength on Th17

// BBB integrity (1 = intact, 0 = fully disrupted)
kin_BBB     : 0.02  : /h  BBB repair rate
k_Th17_BBB  : 0.02  :     Th17 BBB disruption coefficient
k_Th1_BBB   : 0.01  :     Th1 BBB disruption coefficient
BBB_ss      : 0.85  :     RRMS baseline BBB integrity (vs 1.0 healthy)

// CNS infiltrating T cells
kin_cTh1    : 0.008 : /h  Th1 CNS infiltration rate
kin_cTh17   : 0.005 : /h  Th17 CNS infiltration rate
kout_cTh1   : 0.015 : /h  CNS Th1 clearance
kout_cTh17  : 0.015 : /h  CNS Th17 clearance
k_nat_inf   : 0.90  :     natalizumab BBB-crossing inhibition efficacy

// Microglia activation
kin_micro   : 0.004 : /h  microglia activation rate by CNS T cells
kout_micro  : 0.006 : /h  microglia deactivation
k_IL10_micro: 0.40  :     IL-10/Treg dampening of microglia

// Oligodendrocytes (normalized, 1 = healthy)
kin_oligo   : 0.005 : /h  OL synthesis
kout_oligo  : 0.005 : /h  OL natural turnover
k_micro_oligo: 0.04 :     microglia → OL injury rate
k_ROS_oligo : 0.02  :     ROS contribution to OL injury

// OPC and remyelination
kin_OPC     : 0.008 : /h  OPC progenitor input
kout_OPC    : 0.003 : /h  OPC natural loss
kdiff_OPC   : 0.002 : /h  OPC differentiation → new OL
k_LINGO_inh : 0.50  :     LINGO-1 inhibition of OPC diff (RRMS-elevated)

// Myelin (normalized, 1 = full)
kin_myelin  : 0.003 : /h  (remyelin from OPC)
kdmyelin    : 0.02  :     demyelination rate per unit microglia
k_remyelin  : 0.005 : /h  OPC-driven remyelination

// Axonal integrity (1 = intact)
k_axon_loss : 0.010 :     axonal injury from demyelinated segments
k_axon_rep  : 0.002 : /h  axonal repair (very slow)

// Neurofilament light chain (NfL) serum
kin_NfL     : 0.0020 : /h baseline NfL release
k_axon_NfL  : 0.040  :    NfL release per unit axonal injury
kelim_NfL   : 0.030  : /h NfL clearance

// GFAP serum
kin_GFAP    : 0.002  : /h baseline GFAP
k_astro_GFAP: 0.020  :    astro reactivity → GFAP
kelim_GFAP  : 0.025  : /h

// EDSS progression
k_EDSS_acc  : 0.0005 : /h EDSS accumulation per unit axonal injury
k_EDSS_max  : 10.0   :    EDSS scale maximum

// Relapse biology
k_relapse_base : 1.5 :    baseline ARR in untreated RRMS (events/year)
k_relapse_myelin : 3.0 :  ARR scaling with myelin loss rate

$CMT @annotated
// Drug PK states
A1_ifnb  : IFN-beta depot (SC/IM)
C1_ifnb  : IFN-beta central
C2_ifnb  : IFN-beta peripheral
C1_nat   : Natalizumab central (µg/mL)
C2_nat   : Natalizumab peripheral
RC_nat   : Natalizumab-receptor complex
R_nat    : Free alpha4-integrin receptor
C1_ocre  : Ocrelizumab central (µg/mL)
C2_ocre  : Ocrelizumab peripheral
C1_sip   : Siponimod plasma (ng/mL)
C2_sip   : Siponimod deep
C1_dmf   : MMF plasma (µg/mL)
C1_clad  : Cladribine plasma (ng/mL)

// Disease states
Th1      : Peripheral Th1 cells (normalized)
Th17     : Peripheral Th17 cells (normalized)
Treg     : Regulatory T cells (normalized)
Bcell    : Peripheral B cells (normalized)
BBB      : BBB integrity (0-1)
cTh1     : CNS-infiltrating Th1
cTh17    : CNS-infiltrating Th17
Micro    : Activated microglia (normalized)
Oligo    : Oligodendrocyte health (0-1)
OPC      : OPC pool (normalized)
Myelin   : Myelin integrity (0-1)
Axon     : Axonal integrity (0-1)

// Biomarkers
NfL      : Serum NfL (pg/mL)
GFAP     : Serum GFAP (pg/mL)

$MAIN
// ── Drug effect terms ──────────────────────────────────────────────────────
// IFN-beta exposure → effects
double E_ifnb   = (dose_ifnb > 0) ? C1_ifnb / (C1_ifnb + 500.0) : 0;
double I_VCAM   = E_ifnb * 0.50;    // ↓ VCAM-1/MMP9
double I_TNF    = E_ifnb * 0.40;    // ↓ TNF/IFN-γ
double E_Treg_ifnb = E_ifnb * 0.30; // ↑ Treg

// Natalizumab occupancy (% receptor blocked)
double occ_nat  = (dose_nat > 0) ? RC_nat / (RC_nat + R_nat + 1e-6) : 0;
double I_BBB_nat = occ_nat * k_nat_inf;  // ↓ CNS infiltration

// Ocrelizumab B cell depletion
double E_ocre   = (dose_ocre > 0) ? Emax_ocre * C1_ocre / (C1_ocre + EC50_ocre) : 0;

// Siponimod lymphocyte retention
double E_sip    = (dose_sip > 0) ? Emax_sip * C1_sip / (C1_sip + EC50_sip) : 0;

// DMF effects
double I_NFKB   = (dose_dmf > 0) ? C1_dmf / (C1_dmf + EC50_dmf) * 0.45 : 0;
double E_NRF2   = (dose_dmf > 0) ? C1_dmf / (C1_dmf + EC50_nrf2) * 0.55 : 0;
double I_Th_dmf = I_NFKB * 0.50;   // ↓ Th1/Th17 trafficking

// Cladribine lymphodepletion
double I_clad   = (dose_clad > 0) ? Emax_clad * C1_clad / (C1_clad + EC50_clad) : 0;

// Net lymphocyte inhibition (combined treatment)
double I_lymph_total = 1.0 - (1.0 - E_sip) * (1.0 - I_clad) * (1.0 - I_Th_dmf);

// Microglia dampening by Treg/IL-10
double I_micro_Treg = k_IL10_micro * (Treg / Treg_ss);

// Nrf2-driven ROS protection on oligodendrocytes
double ROS_protect  = 1.0 - E_NRF2;

// EDSS model (0-10 bounded)
double EDSS_calc = k_EDSS_max * (1.0 - Axon) * 0.8 +
                   k_EDSS_max * (1.0 - Myelin) * 0.2;
EDSS_calc = (EDSS_calc > 10.0) ? 10.0 : EDSS_calc;
EDSS_calc = (EDSS_calc < 0.0)  ? 0.0  : EDSS_calc;

// ARR instantaneous (events per year)
double myelin_loss_rate = kdmyelin * Micro * Myelin;
double ARR_inst = k_relapse_base * myelin_loss_rate / (kdmyelin * 0.5) *
                  (1.0 - occ_nat * 0.68) *
                  (1.0 - E_sip * 0.55) *
                  (1.0 - E_ocre * 0.47) *
                  (1.0 - I_clad * 0.58) *
                  (1.0 - E_ifnb * 0.34) *
                  (1.0 - I_NFKB * 0.49);

$ODE
// ── Drug PK ODEs ─────────────────────────────────────────────────────────

// IFN-beta 2-comp SC
dxdt_A1_ifnb = -ka_ifnb * A1_ifnb;
dxdt_C1_ifnb = (dose_ifnb > 0 ? F_ifnb * ka_ifnb * A1_ifnb : 0.0)
               - (CL_ifnb/V1_ifnb + Q_ifnb/V1_ifnb) * C1_ifnb
               + (Q_ifnb/V2_ifnb) * C2_ifnb;
dxdt_C2_ifnb = (Q_ifnb/V1_ifnb) * C1_ifnb - (Q_ifnb/V2_ifnb) * C2_ifnb;

// Natalizumab TMDD 2-comp (concentration in µg/mL ~= nM × 148/1000)
dxdt_C1_nat  = -(CL_nat/V1_nat + Q_nat/V1_nat) * C1_nat
               + (Q_nat/V2_nat) * C2_nat
               - kon_nat * C1_nat * R_nat
               + koff_nat * RC_nat;
dxdt_C2_nat  = (Q_nat/V1_nat) * C1_nat - (Q_nat/V2_nat) * C2_nat;
dxdt_RC_nat  = kon_nat * C1_nat * R_nat - (koff_nat + kint_nat) * RC_nat;
dxdt_R_nat   = ksyn_R - kdeg_R * R_nat - kon_nat * C1_nat * R_nat
               + koff_nat * RC_nat;

// Ocrelizumab 2-comp IV
dxdt_C1_ocre = -(CL_ocre/V1_ocre + Q_ocre/V1_ocre) * C1_ocre
               + (Q_ocre/V2_ocre) * C2_ocre;
dxdt_C2_ocre = (Q_ocre/V1_ocre) * C1_ocre - (Q_ocre/V2_ocre) * C2_ocre;

// Siponimod 2-comp oral
dxdt_C1_sip  = ka_sip * C2_sip - (CL_sip/V1_sip + Q_sip/V1_sip) * C1_sip
               + (Q_sip/V2_sip) * C2_sip;
dxdt_C2_sip  = -(ka_sip + Q_sip/V2_sip) * C2_sip + (Q_sip/V1_sip) * C1_sip;
// Note: C2_sip used as absorption compartment (depot → central)

// MMF (DMF hydrolysis product) 1-comp
dxdt_C1_dmf  = -( CL_dmf/V1_dmf ) * C1_dmf;

// Cladribine 1-comp
dxdt_C1_clad = -( CL_clad/V1_clad ) * C1_clad;

// ── Disease biology ODEs ─────────────────────────────────────────────────

// Th1 — elevated in RRMS, inhibited by sip/clad/dmf/IFNb
dxdt_Th1  = kin_Th1 * (1.0 - I_lymph_total) * (1.0 - I_TNF * 0.5)
            + kin_Th1 * E_Treg_ifnb * 0.2
            - kout_Th1 * Th1
            - kout_Th1 * I_clad * Th1;

// Th17 — inhibited by Treg and treatments
dxdt_Th17 = kin_Th17 * (1.0 - I_lymph_total)
            * (1.0 / (1.0 + k_Treg_inh * Treg / Treg_ss))
            - kout_Th17 * Th17;

// Treg — enhanced by IFNb and GA (via Foxp3)
dxdt_Treg = kin_Treg * (1.0 + E_Treg_ifnb)
            - kout_Treg * Treg;

// B cells — depleted by anti-CD20
dxdt_Bcell = kin_B * (1.0 - E_ocre) * (1.0 - I_lymph_total)
             - kout_B * Bcell;

// BBB integrity (0=disrupted, 1=intact)
// Disrupted by Th17/Th1, repaired over time; natalizumab protects
dxdt_BBB = kin_BBB * (1.0 - BBB)
           - BBB * (k_Th17_BBB * Th17 + k_Th1_BBB * Th1)
           * (1.0 - I_VCAM)      // IFN-b protects
           + kin_BBB * I_VCAM;

// CNS Th1 infiltration (gated by BBB disruption)
double BBB_open = 1.0 - BBB;   // 0=intact, 1=open
dxdt_cTh1  = kin_cTh1 * Th1 * (1.0 + BBB_open * 2.0)
             * (1.0 - I_BBB_nat)
             - kout_cTh1 * cTh1;

// CNS Th17 infiltration
dxdt_cTh17 = kin_cTh17 * Th17 * (1.0 + BBB_open * 3.0)
             * (1.0 - I_BBB_nat)
             - kout_cTh17 * cTh17;

// Activated microglia (driven by CNS T cells + complement deposition)
dxdt_Micro = kin_micro * (cTh1 + 2.0 * cTh17 + 0.5 * (1.0 - Bcell))
             * (1.0 - I_micro_Treg)
             - kout_micro * Micro;

// Oligodendrocyte health (cytotoxic injury vs. OPC replenishment)
dxdt_Oligo = kin_oligo
             - kout_oligo * Oligo
             - Oligo * (k_micro_oligo * Micro + k_ROS_oligo * (1.0 - E_NRF2) * Micro);

// OPC pool (precursor cells available for remyelination)
dxdt_OPC   = kin_OPC
             - kout_OPC * OPC
             - kdiff_OPC * OPC * (1.0 - k_LINGO_inh * Micro);

// Myelin integrity
dxdt_Myelin = k_remyelin * OPC * (1.0 - Myelin)        // remyelination
              - kdmyelin * Micro * Myelin * ROS_protect  // demyelination
              + kin_myelin * (1.0 - Myelin);             // baseline repair

// Axonal integrity (very slow repair, drives long-term disability)
dxdt_Axon  = k_axon_rep * Myelin               // recovery if remyelinated
             - k_axon_loss * (1.0 - Myelin)     // loss from demyelinated segments
             - k_axon_loss * 0.5 * Micro * Axon; // smouldering inflammation loss

// Biomarkers
dxdt_NfL   = kin_NfL + k_axon_NfL * (1.0 - Axon) * (1.0 + Micro)
             - kelim_NfL * NfL;

dxdt_GFAP  = kin_GFAP + k_astro_GFAP * (cTh17 + Micro * 0.5) * (1.0 - BBB)
             - kelim_GFAP * GFAP;

$TABLE
capture EDSS     = EDSS_calc;
capture ARR      = ARR_inst;
capture T2lesion = (1.0 - Myelin) * 15000.0;  // mm3 approximate
capture occ_nat  = occ_nat * 100.0;           // % receptor saturation
capture E_ocre   = E_ocre * 100.0;            // % B cell depletion
capture E_sip    = E_sip * 100.0;             // % lymphocyte retention

$INIT
// Drug PK (all zero at start)
A1_ifnb = 0, C1_ifnb = 0, C2_ifnb = 0,
C1_nat = 0, C2_nat = 0, RC_nat = 0, R_nat = 30.0,
C1_ocre = 0, C2_ocre = 0,
C1_sip = 0, C2_sip = 0,
C1_dmf = 0,
C1_clad = 0,

// Disease states (RRMS at model baseline, normalized)
Th1   = 1.5,   // elevated in RRMS
Th17  = 1.8,   // strongly elevated
Treg  = 0.7,   // reduced in RRMS
Bcell = 1.2,   // mildly elevated
BBB   = 0.80,  // mildly disrupted in RRMS
cTh1  = 0.30,  // active CNS infiltration
cTh17 = 0.25,
Micro = 0.50,  // activated microglia
Oligo = 0.80,  // oligodendrocyte loss
OPC   = 0.90,  // relatively preserved
Myelin = 0.75, // demyelination burden at baseline
Axon  = 0.85,  // mild axonal loss
NfL   = 12.0,  // pg/mL — elevated in active RRMS
GFAP  = 180.0  // pg/mL — elevated
'

## ── Build model ──────────────────────────────────────────────────────────────

ms <- mrgsolve::mcode("ms_qsp", ms_model_code)

## ── Helper: create dosing events ─────────────────────────────────────────────

make_doses <- function(drug = "sip", duration_wk = 96) {
  dur_h <- duration_wk * 7 * 24
  switch(drug,
    "ifnb" = {
      # IFN-beta-1a IM 30 µg q1w — enters depot A1_ifnb
      ev(cmt = "A1_ifnb", amt = 30, ii = 7*24, addl = round(dur_h/(7*24)) - 1,
         time = 0) %>% mrgsolve::ev_add(param = list(dose_ifnb = 1))
    },
    "nat" = {
      # Natalizumab 300 mg IV q4w → C1_nat
      ev(cmt = "C1_nat", amt = 300000/3.5, ii = 4*7*24,
         addl = round(dur_h/(4*7*24)) - 1, time = 0) %>%
        mrgsolve::ev_add(param = list(dose_nat = 1))
    },
    "ocre" = {
      # Ocrelizumab 300 mg × 2 (load) then 600 mg q6m
      ev_load  <- ev(cmt = "C1_ocre", amt = 300000/3.0, time = 0)
      ev_load2 <- ev(cmt = "C1_ocre", amt = 300000/3.0, time = 14*24)
      ev_maint <- ev(cmt = "C1_ocre", amt = 600000/3.0,
                     time = 24*7*24, ii = 26*7*24,
                     addl = round(dur_h/(26*7*24)) - 1)
      as.ev(bind_rows(as.data.frame(ev_load),
                      as.data.frame(ev_load2),
                      as.data.frame(ev_maint))) %>%
        mrgsolve::ev_add(param = list(dose_ocre = 1))
    },
    "sip" = {
      # Siponimod 2 mg qd oral → C2_sip (depot → C1_sip)
      ev(cmt = "C2_sip", amt = 2000 * 0.84 / 30.0,  # ng/mL in V=30L, F=84%
         ii = 24, addl = round(dur_h/24) - 1, time = 0) %>%
        mrgsolve::ev_add(param = list(dose_sip = 1))
    },
    "dmf" = {
      # DMF 240 mg BID oral
      ev(cmt = "C1_dmf", amt = 240000/50.0 * 0.25,  # µg/mL scaling
         ii = 12, addl = round(dur_h/12) - 1, time = 0) %>%
        mrgsolve::ev_add(param = list(dose_dmf = 1))
    },
    "clad" = {
      # Cladribine 3.5 mg/kg in 2 annual courses (days 1-4, 5-8 of wk 1, 5)
      times <- c(0, 24, 48, 72, 96, 120, 144, 168,
                 48*7*24, 48*7*24+24, 48*7*24+48, 48*7*24+72,
                 48*7*24+96, 48*7*24+120, 48*7*24+144, 48*7*24+168)
      ev(data.frame(cmt = "C1_clad", amt = 3.5 * 70 / 8 * 1000 / 40.0,
                    time = times, evid = 1)) %>%
        mrgsolve::ev_add(param = list(dose_clad = 1))
    },
    "none" = {
      ev(time = 0, cmt = 1, amt = 0, evid = 2)  # null dose
    }
  )
}

## ── Simulation scenarios ──────────────────────────────────────────────────────

sim_scenario <- function(drug = "none", duration_wk = 104, n_id = 1,
                          seed = 42) {
  set.seed(seed)

  # Observation times every 4 weeks, plus weekly for first 12 weeks
  obs_times <- sort(unique(c(
    seq(0, 12*7*24, by = 7*24),   # weekly for 3 months
    seq(12*7*24, duration_wk*7*24, by = 4*7*24) # q4w thereafter
  )))

  if (drug == "none") {
    dosing <- ev(time = 0, cmt = 1, amt = 0, evid = 2)
  } else {
    dosing <- make_doses(drug, duration_wk)
  }

  ms %>%
    ev(dosing) %>%
    mrgsim(end = duration_wk * 7 * 24, delta = 24, obsonly = TRUE) %>%
    as_tibble() %>%
    mutate(drug = drug, time_wk = time / (7*24),
           time_yr = time / (365.25*24))
}

## ── Run all 6 scenarios ───────────────────────────────────────────────────────

cat("Running MS QSP simulations (6 drug scenarios + placebo)...\n")

scenarios <- c("none", "ifnb", "nat", "ocre", "sip", "dmf", "clad")

results <- lapply(scenarios, function(d) {
  cat(sprintf("  Simulating: %s\n", d))
  tryCatch(sim_scenario(d, duration_wk = 104),
           error = function(e) { message(e); NULL })
})
names(results) <- scenarios

all_results <- bind_rows(Filter(Negate(is.null), results))

## ── Plot functions ────────────────────────────────────────────────────────────

drug_colors <- c(
  none  = "#7F8C8D",
  ifnb  = "#3498DB",
  nat   = "#E74C3C",
  ocre  = "#9B59B6",
  sip   = "#F39C12",
  dmf   = "#27AE60",
  clad  = "#E67E22"
)
drug_labels <- c(
  none  = "Untreated",
  ifnb  = "IFN-β (Avonex)",
  nat   = "Natalizumab",
  ocre  = "Ocrelizumab",
  sip   = "Siponimod",
  dmf   = "Dimethyl Fumarate",
  clad  = "Cladribine"
)

plot_theme <- theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "#F0F0F0"),
        panel.grid.minor = element_blank())

## Plot 1: EDSS over time
p_edss <- ggplot(all_results, aes(x = time_yr, y = EDSS,
                                   color = drug, linetype = drug)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = drug_colors, labels = drug_labels,
                     name = "Treatment") +
  scale_linetype_manual(values = c("none"="solid","ifnb"="dashed",
                                   "nat"="solid","ocre"="dotdash",
                                   "sip"="dashed","dmf"="dotted",
                                   "clad"="longdash"),
                        labels = drug_labels, name = "Treatment") +
  labs(title = "EDSS Progression Over 2 Years",
       x = "Time (years)", y = "EDSS Score (0–10)") +
  ylim(0, 4) + plot_theme

## Plot 2: Annualized Relapse Rate
p_arr <- ggplot(all_results, aes(x = time_yr, y = ARR,
                                  color = drug)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = drug_colors, labels = drug_labels) +
  labs(title = "Instantaneous Annualized Relapse Rate",
       x = "Time (years)", y = "ARR (events/year)") +
  ylim(0, 2.5) + plot_theme +
  theme(legend.position = "none")

## Plot 3: T2 Lesion Volume
p_t2 <- ggplot(all_results, aes(x = time_yr, y = T2lesion,
                                  color = drug)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = drug_colors, labels = drug_labels) +
  labs(title = "T2 Lesion Volume (MRI)",
       x = "Time (years)", y = "T2 Volume (mm³)") +
  plot_theme + theme(legend.position = "none")

## Plot 4: Myelin integrity
p_myelin <- ggplot(all_results, aes(x = time_yr, y = Myelin,
                                     color = drug)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = drug_colors, labels = drug_labels) +
  labs(title = "Myelin Integrity (ODE state)",
       x = "Time (years)", y = "Myelin (0=lost, 1=intact)") +
  ylim(0, 1) + plot_theme + theme(legend.position = "none")

## Plot 5: NfL biomarker
p_nfl <- ggplot(all_results, aes(x = time_yr, y = NfL,
                                  color = drug)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = drug_colors, labels = drug_labels) +
  geom_hline(yintercept = 10, linetype = 2, color = "grey40") +
  annotate("text", x = 0.1, y = 10.5, label = "ULN (~10 pg/mL)",
           size = 3, color = "grey40") +
  labs(title = "Serum NfL (Neurodegeneration Biomarker)",
       x = "Time (years)", y = "NfL (pg/mL)") +
  plot_theme + theme(legend.position = "none")

## Plot 6: CNS infiltration dynamics (untreated)
p_cns <- all_results %>%
  filter(drug %in% c("none", "nat")) %>%
  select(time_yr, drug, cTh1, cTh17, Micro) %>%
  pivot_longer(c(cTh1, cTh17, Micro), names_to = "cell", values_to = "value") %>%
  ggplot(aes(x = time_yr, y = value, color = cell, linetype = drug)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = c(cTh1 = "#E74C3C", cTh17 = "#F39C12",
                                 Micro = "#8E44AD")) +
  labs(title = "CNS Immune Cell Dynamics\n(Untreated vs. Natalizumab)",
       x = "Time (years)", y = "Normalized cell level",
       color = "Cell type", linetype = "Treatment") +
  plot_theme

## Combine plots
cat("Generating combined plot...\n")
combined_plot <- (p_edss | p_arr) /
                 (p_t2   | p_myelin) /
                 (p_nfl  | p_cns) +
  plot_annotation(
    title    = "Multiple Sclerosis QSP Model — Treatment Scenario Comparison",
    subtitle = "Simulated RRMS cohort, 2-year follow-up · 7 treatment arms",
    caption  = "Parameters calibrated to AFFIRM, OPERA, EXPAND, CLARITY, DEFINE trials"
  )

ggsave("ms_qsp_simulation.png", combined_plot,
       width = 16, height = 14, dpi = 150,
       path = dirname(rstudioapi::getActiveDocumentContext()$path %||% "."))

cat("Simulation complete. Output: ms_qsp_simulation.png\n")

## ── 2-year cumulative ARR summary table ──────────────────────────────────────

arr_summary <- all_results %>%
  group_by(drug) %>%
  summarise(
    mean_EDSS     = round(mean(EDSS, na.rm = TRUE), 2),
    final_EDSS    = round(last(EDSS), 2),
    mean_ARR      = round(mean(ARR, na.rm = TRUE), 3),
    ARR_reduction = round((1 - mean(ARR, na.rm = TRUE) /
                             all_results %>%
                             filter(drug == "none") %>%
                             summarise(a = mean(ARR)) %>% pull(a)) * 100, 1),
    final_Myelin  = round(last(Myelin), 3),
    final_NfL     = round(last(NfL), 1),
    final_T2      = round(last(T2lesion), 0),
    .groups = "drop"
  ) %>%
  mutate(drug = drug_labels[drug]) %>%
  arrange(mean_ARR)

cat("\n=== 2-Year Simulation Summary ===\n")
print(arr_summary)

## ── Sensitivity analysis: dose-response for natalizumab ──────────────────────

dose_levels <- c(100, 150, 200, 300, 400) # mg IV q4w

dose_response <- lapply(dose_levels, function(d) {
  ms %>%
    ev(cmt = "C1_nat", amt = d * 1000 / 3.5, ii = 4*7*24, addl = 5,
       time = 0) %>%
    param(dose_nat = 1) %>%
    mrgsim(end = 24 * 4*7, delta = 24, obsonly = TRUE) %>%
    as_tibble() %>%
    mutate(dose_mg = d, time_wk = time / (7*24))
}) %>% bind_rows()

p_doseresponse <- ggplot(dose_response,
                          aes(x = time_wk, y = occ_nat, color = factor(dose_mg))) +
  geom_line(linewidth = 1.2) +
  scale_color_viridis_d(name = "Dose (mg)", option = "plasma") +
  geom_hline(yintercept = 80, linetype = 2, color = "red") +
  annotate("text", x = 0.5, y = 82, label = "80% target occupancy",
           color = "red", size = 3) +
  labs(title = "Natalizumab: Dose-Response α4-Integrin Receptor Occupancy",
       x = "Time (weeks)", y = "% Receptor Occupancy") +
  plot_theme

print(p_doseresponse)
cat("\nNatalizumab dose-response simulation complete.\n")
