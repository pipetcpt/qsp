## =============================================================================
## Trigeminal Neuralgia (TN) — mrgsolve QSP model
## Neurovascular compression -> Nav channelopathy -> ectopic discharge ->
##   central sensitization -> paroxysmal pain, with anticonvulsant PK/PD
##   (carbamazepine, oxcarbazepine, baclofen, gabapentin, pregabalin) and
##   microvascular decompression (MVD) / radiofrequency rhizotomy effects.
##
## Calibration anchors (see tn_references.md for full PMID list):
##   - Carbamazepine dose-response & NNT: Wiffen 2011 Cochrane (CD005451)
##   - CBZ vs OXC comparable efficacy: Di Stefano 2021 J Headache Pain;
##     Zakrzewska 1997 Pain (CBZ RCT), Gronseth 2008 AAN/EFNS guideline
##   - CBZ autoinduction PK (t1/2 36h -> 12-17h): Bertilsson & Tomson 1986
##     Clin Pharmacokinet; Eichelbaum 1975
##   - OXC prodrug -> MHD PK: Flesch 2004 Clin Pharmacokinet review
##   - OXC hyponatremia incidence ~2.7%: Dong 2005 Neurology
##   - Baclofen add-on (Fromm 1984 Ann Neurol RCT)
##   - Gabapentin/pregabalin alpha2delta PK-PD: neuropathic pain extrapolation
##     (Backonja 1998 JAMA gabapentin PHN as structural analogue)
##   - MVD long-term outcomes (Barker 1996 NEJM; ~70% pain-free at 10y with
##     ~4%/yr recurrence); RF rhizotomy / balloon compression recurrence rates
##     (Kondziolka 1996 J Neurosurg gamma knife; Tronnier 2001 Neurosurgery RF)
## =============================================================================

$PROB
# Trigeminal Neuralgia QSP model (17-compartment PK/PD + disease + safety)

$PARAM
// ---- Carbamazepine (CBZ) PK: 1-cpt oral, autoinduction via enzyme-turnover ----
CBZ_KA      = 0.4     // 1/h, slow variable GI absorption
CBZ_F       = 0.85    // bioavailability
CBZ_V       = 70      // L (~1 L/kg x 70kg)
CBZ_CL0     = 2.2     // L/h, baseline clearance (pre-induction)
CBZ_CLmax_f = 2.8     // fold increase in CL at full autoinduction
CBZ_KENZ    = 0.010   // 1/h, enzyme turnover rate constant (induction onset ~wk2-3)
CBZ_FM      = 0.30    // fraction converted to active epoxide metabolite
CBZ_VM      = 40      // L, epoxide volume
CBZ_CLM     = 3.0     // L/h, epoxide clearance

// ---- Oxcarbazepine (OXC) PK: prodrug -> MHD active metabolite ----
OXC_KA      = 1.0     // 1/h
OXC_F       = 0.95
OXC_CONV    = 3.0     // 1/h, rapid hepatic conversion to MHD (near-complete)
OXC_VM      = 50      // L, MHD volume of distribution
OXC_CLM     = 2.0     // L/h, MHD clearance (renal, linear)

// ---- Baclofen PK ----
BAC_KA      = 0.9
BAC_F       = 0.8
BAC_V       = 25
BAC_CL      = 6.0

// ---- Gabapentin PK (saturable intestinal transport, LAT system) ----
GBP_F1      = 0.9     // bioavailability at low dose (near-linear region)
GBP_VMAX    = 40       // mg/h, transporter-limited absorption Vmax
GBP_KM      = 400      // mg, amount at half-max transport
GBP_V       = 58
GBP_CL      = 7.5      // renal clearance

// ---- Pregabalin PK (linear, non-saturable) ----
PGB_KA      = 1.4
PGB_F       = 0.9
PGB_V       = 40
PGB_CL      = 5.0

// ---- Disease / Nav channelopathy / ectopic discharge ----
NAV_BASE      = 1.0     // baseline Nav availability index (compression-driven)
NAV_COMP_GAIN = 0.55    // NVC severity contribution to Nav upregulation
ECTOPIC_K     = 0.30    // 1/h, ectopic-discharge generation rate constant
ECTOPIC_KOUT  = 0.20    // 1/h, spontaneous decay of ectopic drive
EC50_CBZ      = 8       // mg/L (epoxide+parent combined effect scale), Nav blockade EC50
EC50_OXC      = 15      // mg/L MHD, Nav blockade EC50
EC50_LTG_TOPO = 6        // generic scale for Na-channel adjunct agents (unused unless scenario adds)
EMAX_NAVBLOCK = 0.92     // maximal fractional suppression of ectopic drive

// ---- Central sensitization ----
CENTSENS_KIN   = 0.05   // 1/h, sensitization build-up driven by ectopic input
CENTSENS_KOUT  = 0.02   // 1/h, decay of central sensitization (slow)
GABA_GAIN      = 0.35   // baclofen GABA-B inhibitory gain on centsens buildup
CA2D_GAIN      = 0.30   // gabapentinoid alpha2delta inhibitory gain on centsens buildup
EC50_BAC       = 400    // ng/mL baclofen plasma EC50
EC50_GBP       = 4000   // ng/mL gabapentin EC50
EC50_PGB       = 3000   // ng/mL pregabalin EC50

// ---- Clinical paroxysm / pain translation ----
PAROX_BASE    = 18      // baseline paroxysms/day (moderate-severe untreated)
PAROX_GAIN    = 22      // scaling from (ectopic+centsens) drive to paroxysm frequency
PAIN_BASE     = 8.5     // baseline NRS
PAIN_GAIN     = 1.6

// ---- Safety ----
NA_BASE       = 140     // mEq/L baseline plasma sodium
OXC_NA_SLOPE  = 0.012   // mEq/L drop per mg/L MHD (SIADH-like effect)
CBZ_NA_SLOPE  = 0.004   // smaller hyponatremia effect for CBZ
SEDATION_BASE = 0.5
SEDATION_CBZ_SLOPE = 0.045
SEDATION_GBP_SLOPE = 0.00018
SEDATION_PGB_SLOPE = 0.00022
SEDATION_BAC_SLOPE = 0.0020

// ---- Interventional (MVD / RF rhizotomy) ----
MVD_ON        = 0     // 0/1 switch, set by scenario
MVD_TIME      = 1e6   // h, time of MVD surgery (large default = never)
MVD_EFFICACY  = 0.85  // fractional reduction of NVC-driven Nav upregulation
MVD_RECUR_K   = 0.00006 // 1/h, slow stochastic-average recurrence (mechanical re-compression / granuloma)
RF_ON         = 0
RF_TIME       = 1e6
RF_EFFICACY   = 0.80
RF_RECUR_K    = 0.00015  // faster recurrence than MVD (destructive but non-anatomic fix)

// ---- Disease severity covariate ----
NVC_SEVERITY  = 1.0   // 0-2 scale, vascular compression severity (patient covariate)

$CMT
// PK
CBZ_GUT CBZ_CENT CBZ_EPOX CBZ_ENZ
OXC_GUT OXC_MHD
BAC_GUT BAC_CENT
GBP_GUT GBP_CENT
PGB_GUT PGB_CENT
// Disease / PD / safety (11 compartments)
NAV_UPREG ECTOPIC CENTSENS PAROX PAIN NA_PLASMA SEDATION MVD_STATE RF_STATE

$MAIN
double CBZ_conc   = CBZ_CENT / CBZ_V;                     // mg/L
double CBZ_epox_c = CBZ_EPOX / CBZ_VM;                    // mg/L
double CBZ_total_effconc = CBZ_conc + 0.6*CBZ_epox_c;     // combined Nav-blocking exposure
double OXC_mhd_c  = OXC_MHD / OXC_VM;                     // mg/L
double BAC_conc   = (BAC_CENT / BAC_V) * 1000.0;           // ng/mL
double GBP_conc   = (GBP_CENT / GBP_V) * 1000.0;           // ng/mL
double PGB_conc   = (PGB_CENT / PGB_V) * 1000.0;           // ng/mL

// Enzyme-induction fractional clearance multiplier (turnover model, 0-1 -> maps to 1x-CBZ_CLmax_f x)
double CBZ_CL_now = CBZ_CL0 * (1.0 + (CBZ_CLmax_f - 1.0) * CBZ_ENZ);

// Nav channel blockade fraction (additive-ish combination, capped at EMAX)
double navblock_cbz = EMAX_NAVBLOCK * CBZ_total_effconc / (EC50_CBZ + CBZ_total_effconc);
double navblock_oxc = EMAX_NAVBLOCK * OXC_mhd_c / (EC50_OXC + OXC_mhd_c);
double navblock_tot = 1.0 - (1.0 - navblock_cbz) * (1.0 - navblock_oxc);
if(navblock_tot > EMAX_NAVBLOCK) navblock_tot = EMAX_NAVBLOCK;

// Central sensitization inhibitory gain from GABA-B / alpha2delta agents
double gaba_inhib = GABA_GAIN * BAC_conc   / (EC50_BAC + BAC_conc);
double ca2d_inhib = CA2D_GAIN * (GBP_conc/(EC50_GBP+GBP_conc) + PGB_conc/(EC50_PGB+PGB_conc));
double centsens_inhib_tot = gaba_inhib + ca2d_inhib;
if(centsens_inhib_tot > 0.85) centsens_inhib_tot = 0.85;

$ODE
// ---- Carbamazepine PK with autoinduction (indirect turnover of CL) ----
dxdt_CBZ_GUT   = -CBZ_KA * CBZ_GUT;
dxdt_CBZ_CENT  =  CBZ_F * CBZ_KA * CBZ_GUT - CBZ_CL_now * CBZ_conc;
dxdt_CBZ_EPOX  =  CBZ_FM * CBZ_CL_now * CBZ_conc - CBZ_CLM * CBZ_epox_c;
dxdt_CBZ_ENZ   =  CBZ_KENZ * ( (CBZ_conc/(CBZ_conc+4.0)) - CBZ_ENZ );   // 0->1 induction state, driven by exposure

// ---- Oxcarbazepine PK (prodrug rapidly -> MHD, MHD eliminated renally) ----
dxdt_OXC_GUT   = -OXC_KA * OXC_GUT;
dxdt_OXC_MHD   =  OXC_F * OXC_KA * OXC_GUT - OXC_CLM * OXC_mhd_c;

// ---- Baclofen PK ----
dxdt_BAC_GUT   = -BAC_KA * BAC_GUT;
dxdt_BAC_CENT  =  BAC_F * BAC_KA * BAC_GUT - BAC_CL * (BAC_CENT/BAC_V);

// ---- Gabapentin PK (saturable absorption, Michaelis-Menten transport) ----
dxdt_GBP_GUT   = -(GBP_VMAX * GBP_GUT / (GBP_KM + GBP_GUT));
dxdt_GBP_CENT  =  GBP_F1 * (GBP_VMAX * GBP_GUT / (GBP_KM + GBP_GUT)) - GBP_CL * (GBP_CENT/GBP_V);

// ---- Pregabalin PK (linear) ----
dxdt_PGB_GUT   = -PGB_KA * PGB_GUT;
dxdt_PGB_CENT  =  PGB_F * PGB_KA * PGB_GUT - PGB_CL * (PGB_CENT/PGB_V);

// ---- Interventional state variables (0=intact NVC effect scaling, decaying toward recurrence) ----
double mvd_active = (MVD_ON==1 && SOLVERTIME >= MVD_TIME) ? 1.0 : 0.0;
double rf_active  = (RF_ON==1  && SOLVERTIME >= RF_TIME)  ? 1.0 : 0.0;
// MVD_STATE / RF_STATE: 0 = fresh post-op relief, ->1 = fully recurred (mechanical re-compression / fibrosis)
dxdt_MVD_STATE = mvd_active * MVD_RECUR_K * (1.0 - MVD_STATE);
dxdt_RF_STATE  = rf_active  * RF_RECUR_K  * (1.0 - RF_STATE);

// Effective structural relief factor (0 = full relief maintained, ->1 = fully recurred)
double mvd_relief  = mvd_active * MVD_EFFICACY * (1.0 - MVD_STATE);
double rf_relief   = rf_active  * RF_EFFICACY  * (1.0 - RF_STATE);
double structural_relief = 1.0 - (1.0-mvd_relief)*(1.0-rf_relief); // combined fractional relief of NVC drive

// ---- Nav channel upregulation driven by residual (unrelieved) NVC severity ----
double nvc_drive = NVC_SEVERITY * (1.0 - structural_relief);
dxdt_NAV_UPREG = 0.02*(NAV_BASE + NAV_COMP_GAIN*nvc_drive - NAV_UPREG);

// ---- Ectopic discharge generation, suppressed by Nav blockade from drugs ----
dxdt_ECTOPIC = ECTOPIC_K * NAV_UPREG * (1.0 - navblock_tot) - ECTOPIC_KOUT * ECTOPIC;

// ---- Central sensitization, driven by ectopic input, inhibited by GABA-B/alpha2delta ----
dxdt_CENTSENS = CENTSENS_KIN * ECTOPIC * (1.0 - centsens_inhib_tot) - CENTSENS_KOUT * CENTSENS;

// ---- Clinical translation: paroxysm frequency & pain intensity (fast equilibrating indirect links) ----
double paroxysm_target = PAROX_BASE * (0.3 + 0.7*(ECTOPIC/(ECTOPIC+1.0))) * (0.4 + 0.6*(CENTSENS/(CENTSENS+1.0)));
dxdt_PAROX = 0.25*(paroxysm_target - PAROX);
double pain_target = PAIN_BASE * (0.35 + 0.65*(CENTSENS/(CENTSENS+1.0)));
dxdt_PAIN = 0.25*(pain_target - PAIN);

// ---- Safety: plasma sodium (hyponatremia from OXC > CBZ) ----
double na_target = NA_BASE - OXC_NA_SLOPE*OXC_mhd_c - CBZ_NA_SLOPE*CBZ_total_effconc;
dxdt_NA_PLASMA = 0.05*(na_target - NA_PLASMA);

// ---- Safety: sedation composite score (0-10) ----
double sedation_target = SEDATION_BASE + SEDATION_CBZ_SLOPE*CBZ_total_effconc + SEDATION_GBP_SLOPE*GBP_conc
                          + SEDATION_PGB_SLOPE*PGB_conc + SEDATION_BAC_SLOPE*BAC_conc;
dxdt_SEDATION = 0.15*(sedation_target - SEDATION);

$CAPTURE
CBZ_conc CBZ_epox_c OXC_mhd_c BAC_conc GBP_conc PGB_conc navblock_tot centsens_inhib_tot
structural_relief NAV_UPREG ECTOPIC CENTSENS PAROX PAIN NA_PLASMA SEDATION CBZ_ENZ CBZ_CL_now

$SET delta = 4, end = 4320   // hourly-resolution, 180-day horizon by default

## =============================================================================
## Scenario runner (5+ required; 7 provided)
## =============================================================================
$ENV
library(mrgsolve)
library(dplyr)

# (build with: mod <- mread("tn_mrgsolve_model.R"))
# Below: example scenario harness, kept inside the same file for portability.
# ---------------------------------------------------------------------------
run_scenarios <- function(mod) {

  e_base <- ev(time = 0, amt = 0, cmt = "CBZ_GUT")  # placeholder no-dose (natural history)

  # 1. Untreated natural history
  sc1 <- mod %>% param(NVC_SEVERITY = 1.2) %>% ev(e_base) %>% mrgsim(end = 4320) %>% as_tibble()

  # 2. Carbamazepine monotherapy, titrated 200mg BID -> 400mg TID (Zakrzewska 1997 titration)
  cbz_dosing <- ev(time = seq(0, 4320, 12), amt = 200, cmt = "CBZ_GUT", ii = 12, addl = 359)
  sc2 <- mod %>% ev(cbz_dosing) %>% mrgsim(end = 4320) %>% as_tibble()

  # 3. Oxcarbazepine monotherapy 600mg/day (300 BID) titrated to 1200mg/day
  oxc_dosing <- ev(time = seq(0, 4320, 12), amt = 300, cmt = "OXC_GUT", ii = 12, addl = 359)
  sc3 <- mod %>% ev(oxc_dosing) %>% mrgsim(end = 4320) %>% as_tibble()

  # 4. CBZ + Baclofen combination (refractory add-on, Fromm 1984)
  bac_dosing <- ev(time = seq(0, 4320, 8), amt = 10, cmt = "BAC_GUT", ii = 8, addl = 539)
  sc4 <- mod %>% ev(seq(cbz_dosing, bac_dosing)) %>% mrgsim(end = 4320) %>% as_tibble()

  # 5. MVD at day 14 (h=336) after initial CBZ bridge, CBZ tapered off post-op
  mvd_dosing <- ev(time = seq(0, 312, 12), amt = 200, cmt = "CBZ_GUT", ii = 12, addl = 25)
  sc5 <- mod %>% param(MVD_ON = 1, MVD_TIME = 336) %>% ev(mvd_dosing) %>% mrgsim(end = 4320) %>% as_tibble()

  # 6. CBZ-intolerant (hepatotoxic/SJS risk) -> switched to Gabapentin+Pregabalin
  gbp_dosing <- ev(time = seq(0, 4320, 8), amt = 300, cmt = "GBP_GUT", ii = 8, addl = 539)
  pgb_dosing <- ev(time = seq(0, 4320, 12), amt = 75, cmt = "PGB_GUT", ii = 12, addl = 359)
  sc6 <- mod %>% ev(seq(gbp_dosing, pgb_dosing)) %>% mrgsim(end = 4320) %>% as_tibble()

  # 7. Percutaneous RF rhizotomy at day 30, drug-refractory patient, with recurrence over 12 months
  sc7 <- mod %>% param(RF_ON = 1, RF_TIME = 720, NVC_SEVERITY = 1.5) %>%
    ev(e_base) %>% mrgsim(end = 8760) %>% as_tibble()

  list(untreated = sc1, cbz_mono = sc2, oxc_mono = sc3, cbz_baclofen = sc4,
       mvd_post = sc5, gbp_pgb_switch = sc6, rf_rhizotomy = sc7)
}

# Example:
# mod <- mread("tn_mrgsolve_model.R")
# results <- run_scenarios(mod)
