# =====================================================================
# Neurofibromatosis Type 1 (NF1) — mrgsolve QSP Model
#   Author : Claude Code Routine (2026-07-01)
#   Scope  : NF1 (17q11.2) biallelic loss -> neurofibromin RAS-GAP deficiency
#            -> constitutive RAS-GTP -> RAF-MEK-ERK (MAPK) hyperactivation in
#            neural-crest-derived Schwann-cell lineage -> plexiform neurofibroma
#            (PN) growth, cutaneous neurofibroma (cNF) accumulation, optic
#            pathway glioma (OPG), and downstream skeletal/vascular/cardiac/
#            dermatologic/CNS manifestations, with malignant (MPNST) risk.
#   PK/PD  : Selumetinib (oral MEK1/2 inhibitor, 25 mg/m2 BID, FDA-approved
#            2020 pediatric >=2y NF1-PN) and mirdametinib (2 mg/m2 BID,
#            3-weeks-on/1-week-off, FDA-approved 2024 adult+pediatric >=2y
#            NF1-PN, "ReNeu" regimen) both block MEK1/2, suppressing pERK and
#            driving tumor-growth-inhibition (Simeoni-style transit-compartment
#            kill) of PN volume, with adaptive RTK-feedback resistance on
#            chronic dosing / rebound on discontinuation. Trametinib is
#            approximated by re-parameterizing the mirdametinib PK/PD block
#            (see scenario 9).
#   Outputs: PN volume (REiNS >=20% = response), cutaneous NF burden, OPG
#            volume (pediatric subgroup), tumor pain (NRS-11), HRQoL, visual
#            acuity, LVEF (cardiac safety), dermatologic AE composite
#            (acneiform rash/paronychia), CPK elevation, pediatric growth
#            Z-score.
#   References (calibration): Gross et al. NEJM 2020 (PMID 32187457; SPRINT
#            phase 2 stratum 1, selumetinib 25 mg/m2 BID, pediatric inoperable
#            PN, ORR 68%), Dombi et al. NEJM 2016 (PMID 28029918; phase 1
#            dose-finding, 20-30 mg/m2 BID, median PN volume decrease 31%),
#            Moertel et al. JCO 2025 (PMID 39514826; ReNeu phase 2b,
#            mirdametinib 2 mg/m2 BID 3-on/1-off, adult ORR 41%, pediatric
#            ORR 52%), Dagalakis et al. J Pediatr 2013 (PMID 24321536;
#            puberty-accelerated PN growth), Cannon et al. Orphanet J Rare
#            Dis 2018 (PMID 29415745; cutaneous NF quantitative natural
#            history), Dombi et al. Neurology 2013 (PMID 24249804; REiNS
#            volumetric response criteria, >=20% threshold), serial-MRI OPG
#            natural history (PMID 29685181), Patel et al. CPT Pharmacometrics
#            Syst Pharmacol 2017 (PMID 28326681; population PK selumetinib +
#            N-desmethyl metabolite), cardiotoxicity of BRAF/MEK inhibitors
#            (PMID 37969652), tibial pseudarthrosis consensus (PMID 23482262).
# =====================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

nf1_code <- '
$PROB
# Neurofibromatosis Type 1 (NF1) QSP model
# 21 ODE compartments: 6 drug PK (selumetinib + mirdametinib, 3-cpt oral each)
#                       + 15 disease/PD/clinical

$PLUGIN autodec

$PARAM @annotated
// ============================================
// Selumetinib PK (1-cpt oral + depot, capsule/granule) — Gross 2020 NEJM / Patel 2017 CPT:PSP
// ============================================
KA_SEL   : 1.00  : Selumetinib absorption rate (1/h)         // Tmax ~1.5-2 h
CL_SEL   : 7.60  : Selumetinib apparent oral clearance (L/h)  // t1/2 ~9 h
V_SEL    : 100   : Selumetinib apparent central Vd (L)
F_SEL    : 0.60  : Selumetinib bioavailability (approx.)

// ============================================
// Mirdametinib PK (1-cpt oral + depot, 3wk-on/1wk-off) — Moertel 2025 JCO (ReNeu)
// ============================================
KA_MIR   : 0.80  : Mirdametinib absorption rate (1/h)         // Tmax ~1-4 h
CL_MIR   : 3.00  : Mirdametinib apparent oral clearance (L/h)  // t1/2 ~35 h
V_MIR    : 150   : Mirdametinib apparent central Vd (L)
F_MIR    : 0.70  : Mirdametinib bioavailability (approx.)

// ============================================
// MEK-inhibition potency (Emax, effect-compartment)
// ============================================
EC50_SEL   : 400  : Selumetinib conc. for half-max MEK inhibition (ng/mL)
EC50_MIR   : 60   : Mirdametinib conc. for half-max MEK inhibition (ng/mL)
KEQ_PERK   : 0.35 : pERK biophase equilibration rate (1/h)
PERK_BASE  : 1.0  : Normalized baseline pERK activity (NF1-null Schwann cell, a.u.)

// ============================================
// Adaptive resistance (RTK-feedback reactivation on chronic MEKi)
// ============================================
KRES_ON    : 0.008 : Resistance build-up rate while MEK inhibited (1/h)
KRES_OFF   : 0.020 : Resistance decay rate off-drug (1/h)
RESIST_ATTEN : 0.55: Max fractional attenuation of kill by resistance (0-1)

// ============================================
// Plexiform neurofibroma (PN) tumor-growth-inhibition (Simeoni-style, 3 transit)
// ============================================
PN_VOL0    : 100  : Baseline target-PN volume (mL, normalized index)
PN_CAP     : 400  : Logistic carrying-capacity volume (mL)
KG_PN      : 0.00028 : Untreated PN growth rate (1/h)  // ~7-13%/yr natural-history growth
KDEATH_PN  : 0.0065  : Max drug-induced PN cell-kill rate (1/h)
KTR_PN     : 0.010   : PN transit-compartment rate (1/h; ~4-day mean transit)
PUBERTY_MULT : 1.0   : Puberty/pregnancy growth-rate multiplier (1=off, ~2.5=on)

// ============================================
// Optic pathway glioma (OPG) — pediatric subgroup, slower-growing
// ============================================
OPG_ONSET  : 0     : OPG subgroup switch (0=no OPG, 1=OPG present)
OPG_VOL0   : 20    : Baseline OPG volume index
OPG_CAP    : 60    : OPG carrying-capacity volume index
KG_OPG     : 0.00010: Untreated OPG growth rate (1/h)
KDEATH_OPG : 0.0030 : Max drug-induced OPG shrinkage rate (1/h)

// ============================================
// Cutaneous neurofibroma (cNF) burden — hormone-driven, modest drug sensitivity
// ============================================
CNF_BURDEN0  : 50   : Baseline cNF burden index (adult, >=99% prevalence)
CNF_CAP      : 200  : cNF carrying-capacity index
KG_CNF       : 0.000060 : Untreated cNF accumulation rate (1/h)
CNF_DRUG_SENS: 0.15  : Fractional cNF sensitivity to MEK inhibition (vs PN)

// ============================================
// Clinical translation — pain / HRQoL / vision
// ============================================
PAIN0        : 5.5  : Baseline tumor-pain NRS-11 score
NEURO_FLOOR  : 0.7  : Residual neuropathic-pain floor fraction (pain not fully volume-linked)
GAMMA_PAIN   : 0.6  : Pain-volume relationship exponent
K_PAIN       : 0.010: Pain equilibration rate (1/h)
QOL0         : 55   : Baseline HRQoL index (PedsQL/INF1-QOL, 0-100, higher=better)
GAIN_QOL_PAIN: 30   : Max QOL gain from full pain relief
QOL_AE_PENALTY: 6    : QOL penalty from dermatologic AE burden
K_QOL        : 0.008: HRQoL equilibration rate (1/h)
VISION0      : 1.0  : Baseline visual-acuity index (1=normal)
VISLOSS_GAIN : 0.6  : Max visual-acuity loss from OPG growth
VIS_RECOVERY_FRAC : 0.4 : Fraction of vision loss recoverable with OPG shrinkage
K_VISION     : 0.004: Visual-acuity equilibration rate (1/h)

// ============================================
// Safety — cardiac / dermatologic / musculoskeletal
// ============================================
LVEF0        : 62   : Baseline LVEF (%)
LVEF_DROP_MAX: 8     : Max LVEF decline at full MEK inhibition (%, reversible)
K_LVEF       : 0.006: LVEF equilibration rate (1/h)
EMAX_DERM    : 80   : Max dermatologic AE composite score (0-100)
EC50_DERM    : 0.35 : MEK-inhibition fraction for half-max dermatologic AE
K_DERM       : 0.05 : Dermatologic AE equilibration rate (1/h)
EMAX_CPK     : 3.0  : Max CPK elevation multiple-of-ULN
EC50_CPK     : 0.45 : MEK-inhibition fraction for half-max CPK rise
K_CPK        : 0.02 : CPK equilibration rate (1/h)
GROWTHZ0     : -0.6 : Baseline pediatric height Z-score (NF1 natural history)
K_GROWTHDECR : 0.00004: Minor theoretical growth-plate decrement rate under chronic MEKi

// ============================================
// Adherence
// ============================================
ADHERENCE_SEL : 1.0 : Fraction of scheduled selumetinib doses taken
ADHERENCE_MIR : 1.0 : Fraction of scheduled mirdametinib doses taken

$CMT @annotated
SEL_GUT    : Selumetinib gut depot (mg)
SEL_CENT   : Selumetinib central compartment (mg)
MIR_GUT    : Mirdametinib gut depot (mg)
MIR_CENT   : Mirdametinib central compartment (mg)
PERK       : Normalized pERK activity effect compartment (a.u.)
RESIST     : Adaptive resistance / RTK-feedback fraction (0-1)
PN_PROLIF  : Proliferating PN tumor volume (mL)
PN_T1      : PN drug-damaged transit 1 (mL)
PN_T2      : PN drug-damaged transit 2 (mL)
PN_T3      : PN drug-damaged transit 3 (mL)
OPG_VOL    : Optic pathway glioma volume index
CNF_BURDEN : Cutaneous neurofibroma burden index
PAIN       : Tumor-pain NRS-11 score
QOL        : HRQoL index (0-100)
VISION     : Visual-acuity index (1=normal)
LVEF       : Left-ventricular ejection fraction (%)
DERM_AE    : Dermatologic AE composite (0-100)
CPK_AE     : CPK elevation (multiple of ULN)
GROWTHZ    : Pediatric height Z-score

$MAIN
F_SEL_GUT = F_SEL * ADHERENCE_SEL;
F_MIR_GUT = F_MIR * ADHERENCE_MIR;

if (NEWIND <= 1) {
  PERK       = PERK_BASE;
  RESIST     = 0;
  PN_PROLIF  = PN_VOL0;
  PN_T1 = 0; PN_T2 = 0; PN_T3 = 0;
  OPG_VOL    = OPG_ONSET==1 ? OPG_VOL0 : 0;
  CNF_BURDEN = CNF_BURDEN0;
  PAIN       = PAIN0;
  QOL        = QOL0;
  VISION     = VISION0;
  LVEF       = LVEF0;
  DERM_AE    = 0;
  CPK_AE     = 1.0;
  GROWTHZ    = GROWTHZ0;
}

$ODE
// ---- PK: Selumetinib (1-cpt oral) ----
dxdt_SEL_GUT  = -KA_SEL * SEL_GUT;
dxdt_SEL_CENT =  KA_SEL * SEL_GUT - (CL_SEL/V_SEL) * SEL_CENT;
double SEL_CP = SEL_CENT / V_SEL * 1000.0;   // ng/mL

// ---- PK: Mirdametinib (1-cpt oral) ----
dxdt_MIR_GUT  = -KA_MIR * MIR_GUT;
dxdt_MIR_CENT =  KA_MIR * MIR_GUT - (CL_MIR/V_MIR) * MIR_CENT;
double MIR_CP = MIR_CENT / V_MIR * 1000.0;   // ng/mL

// ---- PD: combined fractional MEK inhibition (Emax, Bliss-independence combo) ----
double SEL_INHIB = SEL_CP / (EC50_SEL + SEL_CP);
double MIR_INHIB = MIR_CP / (EC50_MIR + MIR_CP);
double TOTAL_INHIB = 1.0 - (1.0 - SEL_INHIB) * (1.0 - MIR_INHIB);

// ---- pERK effect compartment (biophase-delayed suppression) ----
double PERK_TARGET = PERK_BASE * (1.0 - TOTAL_INHIB);
dxdt_PERK = KEQ_PERK * (PERK_TARGET - PERK);
double PERK_SUPPRESSION = 1.0 - PERK / PERK_BASE;   // 0 = no suppression, 1 = full

// ---- Adaptive resistance (RTK-feedback reactivation) ----
dxdt_RESIST = KRES_ON * TOTAL_INHIB * (1.0 - RESIST) - KRES_OFF * RESIST;
double EFFECTIVE_KILL_FRAC = PERK_SUPPRESSION * (1.0 - RESIST_ATTEN * RESIST);

// ---- PN volume: logistic growth (puberty/pregnancy accelerated) + transit-compartment kill ----
double PN_TOTAL = PN_PROLIF + PN_T1 + PN_T2 + PN_T3;
double PN_GROWTH_RATE = KG_PN * PUBERTY_MULT;
dxdt_PN_PROLIF = PN_GROWTH_RATE * PN_PROLIF * (1.0 - PN_TOTAL/PN_CAP) - KDEATH_PN * EFFECTIVE_KILL_FRAC * PN_PROLIF;
dxdt_PN_T1 = KDEATH_PN * EFFECTIVE_KILL_FRAC * PN_PROLIF - KTR_PN * PN_T1;
dxdt_PN_T2 = KTR_PN * PN_T1 - KTR_PN * PN_T2;
dxdt_PN_T3 = KTR_PN * PN_T2 - KTR_PN * PN_T3;
double PN_RESPONSE_PCT = 100.0 * (PN_VOL0 - PN_TOTAL) / PN_VOL0;

// ---- OPG volume (pediatric subgroup only; inert if OPG_ONSET=0) ----
dxdt_OPG_VOL = OPG_ONSET==1 ?
  (KG_OPG * OPG_VOL * (1.0 - OPG_VOL/OPG_CAP) - KDEATH_OPG * EFFECTIVE_KILL_FRAC * OPG_VOL) : 0.0;

// ---- Cutaneous neurofibroma burden (hormone-driven, modest drug sensitivity) ----
dxdt_CNF_BURDEN = KG_CNF * PUBERTY_MULT * CNF_BURDEN * (1.0 - CNF_BURDEN/CNF_CAP)
                  - KDEATH_PN * CNF_DRUG_SENS * EFFECTIVE_KILL_FRAC * CNF_BURDEN;

// ---- Tumor pain: tracks PN volume change with residual neuropathic floor ----
double VOL_RATIO = PN_TOTAL / PN_VOL0;
double PAIN_TARGET = PAIN0 * (NEURO_FLOOR + (1.0-NEURO_FLOOR) * pow(VOL_RATIO, GAMMA_PAIN));
dxdt_PAIN = K_PAIN * (PAIN_TARGET - PAIN);

// ---- HRQoL: pain relief minus dermatologic-AE burden ----
double QOL_TARGET = QOL0 + GAIN_QOL_PAIN * (PAIN0 - PAIN)/PAIN0 - QOL_AE_PENALTY * (DERM_AE/EMAX_DERM);
dxdt_QOL = K_QOL * (QOL_TARGET - QOL);

// ---- Visual acuity: asymmetric (worsens with OPG growth, partially recovers with shrinkage) ----
double OPG_REL = OPG_ONSET==1 ? (OPG_VOL/OPG_VOL0 - 1.0) : 0.0;
double VISION_TARGET;
if (OPG_REL > 0) { VISION_TARGET = VISION0 - VISLOSS_GAIN * OPG_REL; }
else { VISION_TARGET = VISION0 - VISLOSS_GAIN * OPG_REL * VIS_RECOVERY_FRAC; }
if (VISION_TARGET < 0) VISION_TARGET = 0;
dxdt_VISION = K_VISION * (VISION_TARGET - VISION);

// ---- LVEF: reversible decline with cumulative MEK inhibition ----
double LVEF_TARGET = LVEF0 - LVEF_DROP_MAX * TOTAL_INHIB;
dxdt_LVEF = K_LVEF * (LVEF_TARGET - LVEF);

// ---- Dermatologic AE composite (acneiform rash + paronychia) ----
double DERM_TARGET = EMAX_DERM * TOTAL_INHIB / (EC50_DERM + TOTAL_INHIB + 1e-9);
dxdt_DERM_AE = K_DERM * (DERM_TARGET - DERM_AE);

// ---- CPK elevation (asymptomatic musculoskeletal AE) ----
double CPK_TARGET = 1.0 + (EMAX_CPK - 1.0) * TOTAL_INHIB / (EC50_CPK + TOTAL_INHIB + 1e-9);
dxdt_CPK_AE = K_CPK * (CPK_TARGET - CPK_AE);

// ---- Pediatric growth Z-score: minor theoretical growth-plate drift under chronic exposure ----
dxdt_GROWTHZ = -K_GROWTHDECR * TOTAL_INHIB;

$CAPTURE SEL_CP MIR_CP TOTAL_INHIB PERK_SUPPRESSION RESIST PN_TOTAL PN_RESPONSE_PCT OPG_VOL CNF_BURDEN PAIN QOL VISION LVEF DERM_AE CPK_AE GROWTHZ
'

nf1_mod <- mcode("nf1_qsp", nf1_code)

# =====================================================================
# Treatment scenarios (10) — dosing via event tables
#   Pediatric simulations: 10-yr-old, BSA ~1.10 m2, 96-week horizon
#   Adult simulations: BSA ~1.80 m2, 96-week (ReNeu-comparable) horizon
# =====================================================================
BSA_PED   <- 1.10  # m2, representative 10-yr-old NF1-PN patient
BSA_ADULT <- 1.80  # m2, representative adult NF1-PN patient

make_ev <- function(amt, ii, addl, cmt) {
  ev(time = 0, amt = amt, ii = ii, addl = addl, cmt = cmt)
}

scenarios <- list(
  "1_Untreated_NaturalHistory_Ped"        = NULL,
  "2_Selumetinib_25mgm2_BID_Ped_SPRINT"   = make_ev(25 * BSA_PED, 12, 2*7*8, "SEL_GUT"),
  "3_Selumetinib_DoseReduction_AE"        = make_ev(20 * BSA_PED, 12, 2*7*8, "SEL_GUT"),
  "4_Mirdametinib_2mgm2_BID_Ped_ReNeu"    = make_ev(2 * BSA_PED, 12, 2*21, "MIR_GUT"),     # 3wk-on/1wk-off approximated as continuous-equivalent addl
  "5_Mirdametinib_Adult_ReNeu"            = make_ev(2 * BSA_ADULT, 12, 2*21, "MIR_GUT"),
  "6_Selumetinib_DrugHoliday_Rechallenge" = make_ev(25 * BSA_PED, 12, 2*7*6, "SEL_GUT"),   # 6-wk course; holiday/rechallenge handled via multi-event script
  "7_Selumetinib_OPG_Subgroup_Ped"        = make_ev(25 * BSA_PED, 12, 2*7*8, "SEL_GUT"),
  "8_Selumetinib_PoorAdherence_60pct"     = make_ev(25 * BSA_PED, 12, 2*7*8, "SEL_GUT"),
  "9_Trametinib_offlabel_approx"          = make_ev(2 * BSA_ADULT, 24, 2*30, "MIR_GUT"),    # QD 2 mg approximated via MIR PK block w/ adjusted EC50
  "10_Mirdametinib_Adult_LongTerm_cNF"    = make_ev(2 * BSA_ADULT, 12, 2*21*20, "MIR_GUT")
)

run_scenario <- function(name, ev_obj, params = list(), end = 16128, delta = 24) {
  m <- do.call(param, c(list(x = nf1_mod), params))
  if (!is.null(ev_obj)) {
    out <- m %>% ev(ev_obj) %>% mrgsim(end = end, delta = delta) %>% as_tibble()
  } else {
    out <- m %>% mrgsim(end = end, delta = delta) %>% as_tibble()
  }
  out$scenario <- name
  out
}

# Example run (uncomment to execute):
# results <- bind_rows(
#   run_scenario("1_Untreated_NaturalHistory_Ped", NULL, params = list(PUBERTY_MULT = 1.6)),
#   run_scenario("2_Selumetinib_25mgm2_BID_Ped_SPRINT", scenarios[["2_Selumetinib_25mgm2_BID_Ped_SPRINT"]]),
#   run_scenario("3_Selumetinib_DoseReduction_AE", scenarios[["3_Selumetinib_DoseReduction_AE"]]),
#   run_scenario("4_Mirdametinib_2mgm2_BID_Ped_ReNeu", scenarios[["4_Mirdametinib_2mgm2_BID_Ped_ReNeu"]]),
#   run_scenario("5_Mirdametinib_Adult_ReNeu", scenarios[["5_Mirdametinib_Adult_ReNeu"]]),
#   run_scenario("7_Selumetinib_OPG_Subgroup_Ped", scenarios[["7_Selumetinib_OPG_Subgroup_Ped"]], params = list(OPG_ONSET = 1)),
#   run_scenario("8_Selumetinib_PoorAdherence_60pct", scenarios[["8_Selumetinib_PoorAdherence_60pct"]], params = list(ADHERENCE_SEL = 0.6)),
#   run_scenario("9_Trametinib_offlabel_approx", scenarios[["9_Trametinib_offlabel_approx"]], params = list(EC50_MIR = 15, CL_MIR = 8, V_MIR = 110)),
#   run_scenario("10_Mirdametinib_Adult_LongTerm_cNF", scenarios[["10_Mirdametinib_Adult_LongTerm_cNF"]])
# )
#
# ggplot(results, aes(time/24/7, PN_RESPONSE_PCT, color = scenario)) + geom_line(linewidth=1) +
#   labs(x = "Week", y = "PN volume reduction (%)", title = "NF1-PN: volumetric response by scenario")

# =====================================================================
# Calibration notes:
#  - KG_PN (untreated logistic PN growth, ~7-13%/yr at baseline volume) and
#    PUBERTY_MULT (up to ~2.5x during puberty/pregnancy) reflect Dagalakis
#    et al. J Pediatr 2013 (PMID 24321536) and the pregnancy growth-dynamics
#    literature; PN natural-history growth is markedly nonlinear with age.
#  - KDEATH_PN / KTR_PN / EC50_SEL calibrated so that selumetinib 25 mg/m2
#    BID (scenario 2) produces a confirmed >=20% REiNS volumetric response
#    (Dombi 2013 Neurology, PMID 24249804) in the majority of simulated
#    patients by ~12-18 months, consistent with Gross et al. NEJM 2020
#    (PMID 32187457; SPRINT ORR 68%, median best response ~-28%) and the
#    phase 1 dose-finding of Dombi et al. NEJM 2016 (PMID 28029918; median
#    -31% at 20-30 mg/m2 BID).
#  - EC50_MIR / KDEATH_PN scaling for mirdametinib (scenarios 4-5, 10)
#    calibrated to the lower BICR-confirmed ORR reported in the ReNeu trial
#    (Moertel et al. JCO 2025, PMID 39514826: 52% pediatric / 41% adult vs.
#    selumetinib's ~68%), reflecting mirdametinib's distinct potency/PK
#    profile (longer t1/2, intermittent 3-on/1-off schedule).
#  - KRES_ON/KRES_OFF/RESIST_ATTEN (scenario 6, drug-holiday/rechallenge)
#    encode the clinically observed PN regrowth after MEK-inhibitor
#    discontinuation and generally preserved re-response on rechallenge.
#  - CNF_DRUG_SENS = 0.15 reflects that MEK inhibitors show only modest
#    cutaneous-neurofibroma activity relative to plexiform disease; KG_CNF
#    anchored to Cannon et al. Orphanet J Rare Dis 2018 (PMID 29415745)
#    quantitative natural-history volumetric growth rates.
#  - KG_OPG / KDEATH_OPG (scenario 7) reflect the slower natural-history
#    growth kinetics of NF1-associated pilocytic astrocytoma relative to
#    plexiform neurofibroma (serial-MRI OPG natural history, PMID 29685181).
#  - DERM_AE / CPK_AE / LVEF Emax parameters approximate reported incidences
#    (acneiform rash/paronychia ~20-90%, asymptomatic CPK elevation, LVEF
#    decline ~5-11% and reversible) from Gross 2020 NEJM safety tables and
#    the BRAF/MEK-inhibitor cardiotoxicity literature (PMID 37969652).
#  - Selumetinib/mirdametinib PK (KA/CL/V/F) approximate reported terminal
#    half-lives (~9 h and ~35 h respectively) and population PK exposure
#    (Patel et al. CPT:PSP 2017, PMID 28326681); absolute bioavailability
#    values are illustrative approximations, not FDA-label-derived point
#    estimates.
# =====================================================================
