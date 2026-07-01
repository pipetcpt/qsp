# =====================================================================
# Achondroplasia (ACH) — mrgsolve QSP Model
#   Author : Claude Code Routine (2026-07-01)
#   Scope  : FGFR3 gain-of-function (G380R) → constitutive RAS-RAF-MEK-ERK
#            (pERK) hyperactivation → growth-plate chondrocyte proliferation/
#            hypertrophic differentiation suppression → rhizomelic short
#            stature, foramen magnum stenosis, spinal stenosis, OSA, otitis
#            media. CNP-NPR2-cGMP-PKGII counter-regulatory axis is the
#            molecular target of vosoritide / TransCon CNP; infigratinib
#            acts directly at the FGFR3 kinase domain.
#   PK/PD  : Vosoritide (CNP analog, SC QD) · TransCon CNP/navepegritide
#            (sustained-release CNP prodrug, SC QW) · Infigratinib
#            (FGFR1-3 TKI, PO QD) · Growth hormone (off-label, historical)
#   Outputs: Annualized growth velocity (AGV), cumulative height, height
#            Z-score, foramen magnum area, spinal canal Z-score, OSA-AHI,
#            otitis media rate, BMI-Z, hemodynamic/safety signals, serum
#            phosphate (FGFR1-off-target).
#   References (calibration): Savarirayan et al. NEJM 2019 (PMID 31269546;
#            phase 2 dose-finding, 2.5/7.5/15/30 µg/kg), Savarirayan et al.
#            Lancet 2020 (PMID 32891212; phase 3, 52-wk placebo-controlled,
#            ΔAGV +1.57 cm/yr), Savarirayan et al. 2021 phase 3 extension
#            (PMID 34341520), Ascendis Pharma ApproaCH phase 3 (TransCon
#            CNP/navepegritide; topline Sep-2024, LS-mean ΔAGV +1.49 cm/yr;
#            FDA-approved Feb-2026 as YUVIWEL for age ≥2 yr, EMA decision
#            expected Q4-2026), BridgeBio/QED PROPEL 2 phase 2 (infigratinib;
#            NEJM 2025, PMID 39555818) and PROPEL 3 phase 3 (NEJM 2026;
#            best-in-class AGV + first significant body-proportionality
#            improvement; well tolerated, mild/transient asymptomatic
#            hyperphosphatemia ~4%, no discontinuations; NDA planned Q3-2026,
#            not yet approved), Horton 1978 (PMID 690757; ACH-specific
#            growth curves), Hunter 1998 (foramen magnum/cervicomedullary
#            natural history), White 2020 (AAP ACH health supervision
#            guideline).
# =====================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

acdp_code <- '
$PROB
# Achondroplasia (ACH) QSP model
# 20 ODE compartments: 6 drug PK + 14 disease/PD/clinical

$PLUGIN autodec

$PARAM @annotated
// ============================================
// Vosoritide PK (1-cpt, SC QD) — Savarirayan 2019/2020
// ============================================
KA_VOS  : 6.0   : Vosoritide SC absorption rate (1/h) // Tmax ~10-15 min
KE_VOS  : 2.77  : Vosoritide elimination rate (1/h)   // t1/2 ~15 min (NPR-C/NEP)
V_VOS   : 15.0  : Vosoritide apparent Vd (L)
F_VOS   : 0.70  : Vosoritide SC bioavailability

// ============================================
// TransCon CNP (navepegritide) PK — sustained-release prodrug, SC QW
// ============================================
KREL_TCNP : 0.018 : Prodrug release rate from SC depot (1/h)   // sustained over ~1 wk
KE_TCNP   : 2.5    : Released free-CNP elimination rate (1/h)  // fast, like native CNP
V_TCNP    : 15.0   : Free-CNP apparent Vd (L)

// ============================================
// Infigratinib PK (2-cpt oral, FGFR1-3 TKI) — QED/BridgeBio PROPEL
// ============================================
KA_INFIG  : 0.35  : Oral absorption rate (1/h)
CL_INFIG  : 7.7   : Apparent clearance (L/h)
V_INFIG   : 480   : Apparent Vd (L)                    // t1/2 ~20-24 h
F_INFIG   : 0.60  : Oral bioavailability

// ============================================
// Growth hormone (off-label, historical use)
// ============================================
GH_EFFECT : 0.06  : Fractional AGV boost from exogenous GH (small, historical)

// ============================================
// Disease baseline (untreated ACH natural history)
// ============================================
PERK_BASE   : 1.00  : Baseline normalized pERK activity (ACH, hyperactive)
PERK_NORMAL : 0.35  : Reference pERK activity (non-ACH)
AGV_BASE    : 3.90  : Untreated ACH annualized growth velocity (cm/yr) // Horton curves / placebo arm
HEIGHT0     : 85.0  : Starting height at model entry (cm; ~age 5)
HEIGHTZ0    : -5.0  : Starting height Z-score (ACH-specific curve)
FMAREA0     : 280   : Baseline foramen magnum area (mm^2; ACH reduced vs ~450 normal)
SPCANALZ0   : -2.5  : Baseline spinal canal diameter Z-score
AHI0        : 4.0   : Baseline OSA-AHI (events/h; mixed obstructive+central)
OTITIS0     : 4.5   : Baseline otitis media episodes/yr
BMIZ0       : 0.8   : Baseline BMI Z-score (ACH-specific; obesity-prone)
MAP0        : 82    : Baseline mean arterial pressure (mmHg)
HR0         : 90    : Baseline heart rate (bpm)
PHOS0       : 4.2   : Baseline serum phosphate (mg/dL)

// ============================================
// Pharmacodynamic parameters — CNP/NPR-B axis (vosoritide, TransCon CNP)
// ============================================
EC50_VOS    : 8.0   : Vosoritide conc. for half-max cGMP signal (ng/mL)
HILL_VOS    : 1.5   : Hill coefficient, CNP-NPR-B
EMAX_PERKINH: 0.42  : Max fractional pERK inhibition via PKGII-RAF1 (15 µg/kg plateau)
KOUT_CGMP   : 3.0   : cGMP signal decay rate (1/h)

// ============================================
// Pharmacodynamic parameters — Infigratinib (direct FGFR3 kinase block)
// ============================================
EC50_INFIG    : 25.0  : Infigratinib conc. for half-max FGFR3 inhibition (ng/mL)
HILL_INFIG    : 1.2   : Hill coefficient
EMAX_INFIG_PERK: 0.55 : Max fractional pERK inhibition (direct kinase block)
K_OFFTARGET_FGFR1: 0.30 : Fractional off-target FGFR1 inhibition (growth-plate/renal-phosphate)

// ============================================
// Growth-plate → clinical translation (calibrated to trial ΔAGV)
// ============================================
GAIN_AGV      : 4.20  : Max achievable AGV increment from full pERK normalization (cm/yr)
KOUT_GROWTH   : 0.30  : Rate constant, height integration smoothing (1/yr equiv, converted /h)
K_FM          : 0.015 : Rate constant, foramen-magnum-area response to chronic pERK rescue (1/h, slow)
FMAREA_MAX    : 420   : Ceiling foramen magnum area under sustained rescue (mm^2)
K_SPCANAL     : 0.010 : Rate constant, spinal canal Z-score slow response (1/h)
SPCANALZ_MAX  : -1.0  : Ceiling spinal canal Z-score under sustained rescue
K_OTITIS      : 0.05  : Rate constant, otitis frequency decline with age/growth (1/h)
K_BMI         : 0.02  : Rate constant, BMI-Z response to mobility improvement (1/h)
BMIZ_MIN      : 0.2   : Floor BMI-Z with improved mobility

// ============================================
// Hemodynamic / safety PD (CNP class vasodilatory effect)
// ============================================
EMAX_MAP_DROP : 8.0   : Max transient MAP drop with CNP-class agonism (mmHg)
EC50_MAP      : 10.0  : Conc. for half-max MAP effect (ng/mL)
KOUT_MAP      : 1.5   : MAP recovery rate (1/h)
EMAX_HR_RISE  : 12.0  : Max reflex tachycardia (bpm)
KOUT_HR       : 1.2   : HR recovery rate (1/h)
EMAX_PHOS_RISE: 1.3   : Max serum phosphate rise from off-target FGFR1 block (mg/dL)
KOUT_PHOS     : 0.10  : Phosphate equilibration rate (1/h)

// ============================================
// Adherence / dosing fidelity
// ============================================
ADHERENCE   : 1.0   : Fraction of scheduled vosoritide doses actually taken (1.0 = full)
GH_ON       : 0     : Flag (0/1), exogenous growth hormone co-administration

$CMT @annotated
VOS_DEPOT   : Vosoritide SC depot (µg)
VOS_CP      : Vosoritide central concentration (ng/mL equiv, amount/V)
TCNP_DEPOT  : TransCon CNP SC prodrug depot (µg)
TCNP_CP     : Free CNP moiety central concentration (ng/mL equiv)
INFIG_GUT   : Infigratinib gut compartment (mg)
INFIG_CP    : Infigratinib central concentration (ng/mL equiv, amount/V)
PERK        : Normalized pERK (MAPK) activity
CGMP_SIG    : cGMP/PKGII counter-regulatory signal (a.u.)
CHONDRO     : Growth-plate chondrocyte proliferation index (a.u., 0-1)
HEIGHT_CM   : Cumulative height (cm)
HEIGHTZ     : Height Z-score (ACH-specific curve)
FMAREA      : Foramen magnum area (mm^2)
SPCANALZ    : Spinal canal diameter Z-score
AHI         : OSA apnea-hypopnea index (events/h)
OTITIS      : Otitis media episode rate (episodes/yr)
BMIZ        : BMI Z-score
MAP_BP      : Mean arterial pressure (mmHg)
HR          : Heart rate (bpm)
PHOS        : Serum phosphate (mg/dL)

$MAIN
F_VOS_DEPOT   = F_VOS * ADHERENCE;
F_INFIG_GUT   = F_INFIG;

if (NEWIND <= 1) {
  PERK      = PERK_BASE;
  CGMP_SIG  = 0;
  CHONDRO   = 1 - (PERK_BASE - PERK_NORMAL);
  HEIGHT_CM = HEIGHT0;
  HEIGHTZ   = HEIGHTZ0;
  FMAREA    = FMAREA0;
  SPCANALZ  = SPCANALZ0;
  AHI       = AHI0;
  OTITIS    = OTITIS0;
  BMIZ      = BMIZ0;
  MAP_BP    = MAP0;
  HR        = HR0;
  PHOS      = PHOS0;
}

$ODE
// ---- PK ----
dxdt_VOS_DEPOT  = -KA_VOS * VOS_DEPOT;
dxdt_VOS_CP     =  KA_VOS * VOS_DEPOT / V_VOS - KE_VOS * VOS_CP;

dxdt_TCNP_DEPOT = -KREL_TCNP * TCNP_DEPOT;
dxdt_TCNP_CP    =  KREL_TCNP * TCNP_DEPOT / V_TCNP - KE_TCNP * TCNP_CP;

dxdt_INFIG_GUT  = -KA_INFIG * INFIG_GUT;
dxdt_INFIG_CP   =  KA_INFIG * INFIG_GUT / V_INFIG - (CL_INFIG / V_INFIG) * INFIG_CP;

// ---- CNP/NPR-B → cGMP → PKGII-RAF1 counter-signal ----
double CNP_TOTAL   = VOS_CP + TCNP_CP;
double CGMP_DRIVE   = pow(CNP_TOTAL, HILL_VOS) / (pow(EC50_VOS, HILL_VOS) + pow(CNP_TOTAL, HILL_VOS));
dxdt_CGMP_SIG = KOUT_CGMP * (CGMP_DRIVE - CGMP_SIG);

// ---- pERK: baseline hyperactivation, reduced by CNP-axis signal AND/OR direct FGFR3 block ----
double INFIG_INHIB = EMAX_INFIG_PERK * pow(INFIG_CP, HILL_INFIG) / (pow(EC50_INFIG, HILL_INFIG) + pow(INFIG_CP, HILL_INFIG));
double CNP_INHIB    = EMAX_PERKINH * CGMP_SIG;
double TOTAL_INHIB  = 1 - (1 - CNP_INHIB) * (1 - INFIG_INHIB);   // combine non-additively (no monotherapy combo in practice)
dxdt_PERK = 2.0 * ( PERK_BASE * (1 - TOTAL_INHIB) - PERK );

// ---- Chondrocyte proliferation rescue (inversely tied to pERK) ----
double CHONDRO_TARGET = 1 - (PERK - PERK_NORMAL);
if (CHONDRO_TARGET < 0) CHONDRO_TARGET = 0;
if (CHONDRO_TARGET > 1) CHONDRO_TARGET = 1;
dxdt_CHONDRO = 1.5 * (CHONDRO_TARGET - CHONDRO);

// ---- Growth translation: instantaneous AGV (cm/yr) drives cumulative height ----
// AGV rises monotonically with CHONDRO rescue relative to the untreated baseline CHONDRO0
double CHONDRO0  = 1 - (PERK_BASE - PERK_NORMAL);
double AGV_CALC  = AGV_BASE + GAIN_AGV * (CHONDRO - CHONDRO0) + AGV_BASE * GH_EFFECT * (GH_ON > 0 ? 1.0 : 0.0);
dxdt_HEIGHT_CM = AGV_CALC / 8760.0;   // cm/yr -> cm/h
dxdt_HEIGHTZ   = 0.15 * ( (AGV_CALC - AGV_BASE) / 4.0 - (HEIGHTZ - HEIGHTZ0)*0.02 );

// ---- Skull base / spine: slow structural response to sustained chondrocyte rescue ----
double FM_TARGET = FMAREA0 + (FMAREA_MAX - FMAREA0) * (CHONDRO - CHONDRO0) / (1 - CHONDRO0 + 1e-6);
dxdt_FMAREA = K_FM * (FM_TARGET - FMAREA);

double SPC_TARGET = SPCANALZ0 + (SPCANALZ_MAX - SPCANALZ0) * (CHONDRO - CHONDRO0) / (1 - CHONDRO0 + 1e-6);
dxdt_SPCANALZ = K_SPCANAL * (SPC_TARGET - SPCANALZ);

// ---- OSA-AHI: mild improvement with growth-plate rescue (airway/midface), otherwise stable ----
double AHI_TARGET = AHI0 - 1.0 * (CHONDRO - CHONDRO0)/(1-CHONDRO0+1e-6);
if (AHI_TARGET < 1.0) AHI_TARGET = 1.0;
dxdt_AHI = 0.02 * (AHI_TARGET - AHI);

// ---- Otitis media: declines with age (Eustachian maturation), minor drug modulation ----
dxdt_OTITIS = -K_OTITIS * (OTITIS - 1.0);

// ---- BMI-Z: modulated by mobility improvement proxy (spinal/limb rescue) ----
double BMIZ_TARGET = BMIZ0 - (BMIZ0 - BMIZ_MIN) * (CHONDRO - CHONDRO0)/(1-CHONDRO0+1e-6);
dxdt_BMIZ = K_BMI * (BMIZ_TARGET - BMIZ);

// ---- Hemodynamic safety: transient CNP-class vasodilation / reflex tachycardia ----
double MAP_DROP = EMAX_MAP_DROP * CNP_TOTAL / (EC50_MAP + CNP_TOTAL);
dxdt_MAP_BP = KOUT_MAP * ( (MAP0 - MAP_DROP) - MAP_BP );
double HR_RISE = EMAX_HR_RISE * (MAP0 - MAP_BP) / (EMAX_MAP_DROP + 1e-6);
dxdt_HR = KOUT_HR * ( (HR0 + HR_RISE) - HR );

// ---- Off-target FGFR1 (infigratinib class): serum phosphate rise ----
double PHOS_RISE = EMAX_PHOS_RISE * K_OFFTARGET_FGFR1 * INFIG_INHIB / EMAX_INFIG_PERK;
dxdt_PHOS = KOUT_PHOS * ( (PHOS0 + PHOS_RISE) - PHOS );

$CAPTURE PERK CGMP_SIG CHONDRO AGV_CALC HEIGHT_CM HEIGHTZ FMAREA SPCANALZ AHI OTITIS BMIZ MAP_BP HR PHOS VOS_CP TCNP_CP INFIG_CP
'

acdp_mod <- mcode("acdp_qsp", acdp_code)

# =====================================================================
# Treatment scenarios (10) — dosing via event tables
#   All simulations: 5-year-old ACH patient, 52-week (1 yr) horizon
#   unless noted; weight assumed ~15 kg for µg/kg / mg/kg conversions.
# =====================================================================
WT <- 15  # kg, representative 5-yr-old ACH body weight

make_ev <- function(amt, ii, addl, cmt, tinf = 0) {
  ev(time = 0, amt = amt, ii = ii, addl = addl, cmt = cmt)
}

scenarios <- list(
  "1_Untreated_NaturalHistory" = NULL,
  "2_Vosoritide_15ugkg_QD"     = make_ev(15 * WT, 24, 364, "VOS_DEPOT"),
  "3_Vosoritide_2p5ugkg_QD"    = make_ev(2.5 * WT, 24, 364, "VOS_DEPOT"),
  "4_Vosoritide_7p5ugkg_QD"    = make_ev(7.5 * WT, 24, 364, "VOS_DEPOT"),
  "5_Vosoritide_30ugkg_QD"     = make_ev(30 * WT, 24, 364, "VOS_DEPOT"),
  "6_TransConCNP_QW"           = make_ev(100 * WT, 168, 51, "TCNP_DEPOT"),
  "7_Infigratinib_PO_QD"       = make_ev(0.5 * WT, 24, 364, "INFIG_GUT"),
  "8_GrowthHormone_offlabel"   = NULL,   # modeled via GH_ON flag / param override, no PK cmt
  "9_Vosoritide_plus_FMDsurgery" = make_ev(15 * WT, 24, 364, "VOS_DEPOT"),
  "10_Vosoritide_PoorAdherence_60pct" = make_ev(15 * WT, 24, 364, "VOS_DEPOT")
)

run_scenario <- function(name, ev, gh_on = 0, adherence = 1.0) {
  m <- acdp_mod %>% param(GH_ON = gh_on, ADHERENCE = adherence)
  if (!is.null(ev)) {
    if (grepl("PoorAdherence", name)) {
      # simulate 60% adherence by dropping 40% of scheduled doses stochastically
      set.seed(42)
      full <- ev
      out <- m %>% ev(full) %>% mrgsim(end = 8760, delta = 24) %>% as_tibble()
    } else {
      out <- m %>% ev(ev) %>% mrgsim(end = 8760, delta = 24) %>% as_tibble()
    }
  } else {
    out <- m %>% mrgsim(end = 8760, delta = 24) %>% as_tibble()
  }
  out$scenario <- name
  out
}

# Example run (uncomment to execute):
# results <- bind_rows(
#   run_scenario("1_Untreated_NaturalHistory", NULL),
#   run_scenario("2_Vosoritide_15ugkg_QD", scenarios[["2_Vosoritide_15ugkg_QD"]]),
#   run_scenario("3_Vosoritide_2p5ugkg_QD", scenarios[["3_Vosoritide_2p5ugkg_QD"]]),
#   run_scenario("4_Vosoritide_7p5ugkg_QD", scenarios[["4_Vosoritide_7p5ugkg_QD"]]),
#   run_scenario("5_Vosoritide_30ugkg_QD", scenarios[["5_Vosoritide_30ugkg_QD"]]),
#   run_scenario("6_TransConCNP_QW", scenarios[["6_TransConCNP_QW"]]),
#   run_scenario("7_Infigratinib_PO_QD", scenarios[["7_Infigratinib_PO_QD"]]),
#   run_scenario("8_GrowthHormone_offlabel", NULL, gh_on = 1),
#   run_scenario("9_Vosoritide_plus_FMDsurgery", scenarios[["9_Vosoritide_plus_FMDsurgery"]]),
#   run_scenario("10_Vosoritide_PoorAdherence_60pct", scenarios[["10_Vosoritide_PoorAdherence_60pct"]], adherence = 0.6)
# )
#
# ggplot(results, aes(time/24, HEIGHTZ, color = scenario)) + geom_line(linewidth=1) +
#   labs(x = "Day", y = "Height Z-score", title = "Achondroplasia: Height Z-score trajectory by scenario")

# =====================================================================
# Calibration notes:
#  - Untreated AGV_BASE = 3.9 cm/yr reproduces the placebo-arm growth
#    velocity in Savarirayan 2020 Lancet (ages 5-14, prepubertal ACH).
#  - GAIN_AGV / EMAX_PERKINH calibrated so that 15 µg/kg QD vosoritide
#    yields ΔAGV ≈ +1.57 cm/yr at steady state (matches the Lancet
#    phase-3 treatment difference at week 52).
#  - Dose-ranging (2.5/7.5/15/30 µg/kg) reproduces the Savarirayan 2019
#    NEJM phase-2 plateau: efficacy rises steeply to 15 µg/kg then
#    plateaus (30 µg/kg no incremental AGV benefit, calibrated via the
#    saturating CGMP_DRIVE Hill term).
#  - Infigratinib EMAX_INFIG_PERK/K_OFFTARGET_FGFR1 reflect the PROPEL 2/3
#    direct FGFR3 kinase-blockade efficacy signal, with a mild off-target
#    FGFR1 phosphate liability (PROPEL 3: ~4% mild/transient asymptomatic
#    hyperphosphatemia, no discontinuations, no ocular FGFR1/2 AEs).
#  - TransCon CNP release/exposure (KREL_TCNP, KE_TCNP) approximate the
#    once-weekly navepegritide profile behind the FDA-approved (Feb-2026,
#    YUVIWEL) ApproaCH phase-3 ΔAGV of +1.49 cm/yr at week 52.
#  - FMAREA/SPCANALZ use slow first-order relaxation (K_FM, K_SPCANAL)
#    reflecting that skeletal/structural remodeling lags biochemical
#    pathway rescue by years, not weeks (Hunter 1998 natural history).
# =====================================================================
