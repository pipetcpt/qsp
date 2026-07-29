## =============================================================================
## cipn_mrgsolve_model.R
## Chemotherapy-Induced Peripheral Neuropathy (CIPN) — QSP model
## 항암화학요법 유발 말초신경병증 정량적 시스템 약리학 모델
##
##   34 ODEs · 4 drug PK models · 18 mechanism/biomarker states
##   18 treatment scenarios · virtual population · therapeutic-index optimiser
##
## MODEL OF RECORD.  A stdlib-only Python re-implementation of the identical
## ODE system and parameter set lives in `cipn_reference_model.py`; running
## that script prints every number quoted in README.md.  Use it to verify this
## file if you do not have R/mrgsolve available.
##
## -----------------------------------------------------------------------------
## THE THREE STRUCTURAL COMMITMENTS THAT GENERATE THE INTERESTING BEHAVIOUR
## -----------------------------------------------------------------------------
##  1. COASTING IS EMERGENT, NOT IMPOSED.  The SARM1 axon-death program has a
##     slow off-rate (t1/2 = 23 d) and regeneration is GATED by SARM1 activity
##     (an axon cannot regrow while its own degeneration program is running).
##     Severity therefore keeps worsening for 4-12 weeks after the last dose
##     without any post-treatment forcing term.
##
##  2. Cmax vs AUC SEPARATE TOXICITY FROM EFFICACY.  Nerve injury is a CONVEX
##     (Hill, h = 1.51) function of DRG proteasome occupancy, while tumour kill
##     is a SATURATING function of blood occupancy.  Two routes with identical
##     AUC but an 11-fold Cmax difference therefore have very different
##     neurotoxicity and identical efficacy.  The DRG occupancy AUC differs by
##     only ~0.85x between SC and IV; it is the convexity, not the exposure,
##     that produces the toxicity difference.
##
##  3. THE ANTI-TUMOUR EXPOSURE-RESPONSE SATURATES; NEUROTOXICITY DOES NOT.
##     Oxaliplatin's log-kill is Emax*CUM/(CUM + 250 mg/m2) while axonal injury
##     keeps integrating.  A therapeutic-index optimum in cumulative dose
##     therefore exists and is DERIVED, not assumed — reproducing the IDEA
##     3-vs-6-month result.
##
## -----------------------------------------------------------------------------
## CALIBRATION: SEVEN fitted parameters against SEVEN trial observations
## -----------------------------------------------------------------------------
##   KDAM_PT    <- IDEA FOLFOX 6 months, grade >=2 CIPN 47.7%
##   SIGMA_S    <- MOSAIC FOLFOX 6 months, grade >=3 CIPN 12.4%
##   PHI_RELIEF <- IDEA FOLFOX 3 months, grade >=2 CIPN 16.6%
##   KDAM_TAX   <- ECOG 1199 weekly paclitaxel 80 mg/m2 x12, grade >=2 27%
##   KDAM_BTZ   <- MMY-3021 bortezomib 1.3 mg/m2 IV, grade >=2 PN 41%
##   BTZ_JH     <- MMY-3021 bortezomib 1.3 mg/m2 SC, grade >=2 PN 24%
##   KCS_DUL    <- Smith 2013 JAMA, duloxetine net -0.73 BPI over 5 weeks
##
## Every other number the model produces is an OUT-OF-SAMPLE PREDICTION:
##   CAPOX 3 / 6 months grade >=2   (obs ~15% / ~45%; model 18.4% / 50.0%)
##   paclitaxel q3w 175 x4          (obs 20%;   model 14.9%)
##   FOLFOX 3 months grade >=3      (obs 2.7%;  model 0.8%)
##   MOSAIC grade 3 recovery at 12 and 48 months (obs 1.1% / 0.7%)
##   IDEA 3-year DFS in both risk strata
##   coasting: peak severity 4-12 weeks after the last dose
##   the BPI-vs-IENFD dissociation under duloxetine
##
## KNOWN MISSES, stated up front rather than buried:
##   - grade >=3 is systematically under-predicted (a single lognormal
##     susceptibility cannot widen the upper tail without inflating grade >=2,
##     which is pinned by the fit).
##   - once-weekly bortezomib is over-rewarded: the direction of the schedule
##     benefit is right but the magnitude is exaggerated, because inter-dose
##     recovery of mitochondrial and transport capacity is nearly complete
##     over 7 days in this model.
##
## References: see cipn_references.md (102 entries).
## =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)
library(tidyr)

## =============================================================================
## 1.  THE MODEL
## =============================================================================

cipn_code <- '
$PROB
CIPN QSP model: neurotoxic exposure -> DRG/axon injury -> committed
degeneration & coasting -> hyperexcitability -> CIPN20 / CTCAE / dose intensity

$GLOBAL
#define HILL(x, k, n) (((x) <= 0.0) ? 0.0 : (pow((x),(n)) / (pow((x),(n)) + pow((k),(n)))))

$PARAM @annotated
// ---------------------------------------------------------------- demographics
BSA      : 1.8   : Body surface area (m2)

// ------------------------------------------- oxaliplatin PK (ultrafilterable Pt)
// Graham 2000 Clin Cancer Res 6:1205 - Cmax ~0.81 ug/mL after 85 mg/m2 / 2 h,
// CL 13.3 L/h, t1/2 alpha 0.43 h, t1/2 gamma 391 h. The deep compartment IS the
// long-lived tissue-bound platinum that keeps driving injury after treatment.
OXA_V1   : 100   : Oxaliplatin central volume (L)
OXA_CL   : 320   : Oxaliplatin clearance (L/day)
OXA_Q2   : 1920  : Oxaliplatin shallow intercompartmental clearance (L/day)
OXA_V2   : 350   : Oxaliplatin shallow peripheral volume (L)
OXA_Q3   : 24    : Oxaliplatin deep intercompartmental clearance (L/day)
OXA_V3   : 500   : Oxaliplatin deep tissue-bound volume (L)

// ---------------------------------- paclitaxel PK (3-cmt, saturable elimination)
// Henningsson 2001 JCO 19:4065 reports Vmax ~30 mg/h, Km ~0.73 mg/L. We use a
// MILDER saturation here (Km 8 mg/L) because the strongly saturable version
// makes total AUC superproportional to dose, and the model then predicts q3w
// 175 mg/m2 to be more neurotoxic than weekly 80 mg/m2 - the opposite of
// ECOG 1199. The cost is that the q3w Cmax is under-predicted (~2.5 vs ~4.3
// umol/L observed); AUC and time-above-threshold, which are the established
// taxane PD drivers, are reproduced within ~10%.
PAC_V1   : 11    : Paclitaxel central volume (L)
PAC_VMAX : 3920  : Paclitaxel maximum elimination rate (mg/day)
PAC_KM   : 8.0   : Paclitaxel Michaelis constant (mg/L) - mild saturation,
         :       : chosen so total AUC stays dose-proportional (the taxane PD driver)
PAC_Q2   : 700   : Paclitaxel shallow intercompartmental clearance (L/day)
PAC_V2   : 60    : Paclitaxel shallow peripheral volume (L)
PAC_Q3   : 150   : Paclitaxel deep intercompartmental clearance (L/day)
PAC_V3   : 250   : Paclitaxel deep peripheral volume (L)
PAC_THR  : 0.0427: Threshold for time-above-concentration (mg/L = 0.05 umol/L)

// ------------------------------------------ bortezomib PK (SC depot + 2-cmt)
// Moreau 2011 Lancet Oncol 12:431 (MMY-3021) - IV bolus Cmax 223 ng/mL,
// SC Cmax 20.4 ng/mL, AUC equivalent (151 vs 155 ng.h/mL).
BTZ_V1   : 10.5  : Bortezomib central volume (L)
BTZ_CL   : 2448  : Bortezomib clearance (L/day)
BTZ_Q    : 600   : Bortezomib intercompartmental clearance (L/day)
BTZ_V2   : 500   : Bortezomib peripheral volume (L)
BTZ_KA   : 24    : Bortezomib subcutaneous absorption rate (1/day)

// ------------------------------------------------- duloxetine / pregabalin PK
DUL_V    : 1640  : Duloxetine apparent volume (L)
DUL_CL   : 2424  : Duloxetine apparent clearance (L/day)
DUL_KA   : 57.6  : Duloxetine absorption rate (1/day)
DUL_F    : 0.5   : Duloxetine bioavailability
PGB_V    : 39    : Pregabalin apparent volume (L)
PGB_CL   : 115   : Pregabalin apparent clearance (L/day)
PGB_F    : 0.9   : Pregabalin bioavailability

// --------------------------------------------------------- biophase / effect
// KOUT_PT: DRG platinum washes out with t1/2 = 8 d, so exposure outlives the
// infusion - one of the two clocks behind coasting.
KIN_PT   : 1.00  : DRG platinum uptake rate constant (1/day per mg/L)
KOUT_PT  : 0.0866: DRG platinum elimination rate (1/day)
KIN_TAX  : 1.00  : Taxane nerve-burden uptake rate (1/day per mg/L)
KOUT_TAX : 0.0578: Taxane nerve-burden clearance (1/day) - t1/2 = 12 d. Taxanes
         :       : are retained in peripheral nerve far longer than in plasma,
         :       : so the burden integrates CUMULATIVE dose. Modelling this as a
         :       : fast effect-site equilibration instead makes injury track peak
         :       : concentration, and the model then ranks q3w 175 mg/m2 as MORE
         :       : neurotoxic than weekly 80 mg/m2 - the opposite of ECOG 1199.
KD_TAX   : 0.60  : Taxane tubulin-occupancy Kd (mg/L)
KEO_BTZ  : 200   : Bortezomib DRG equilibration rate (1/day) - FAST because the
                 : DRG capillaries are fenestrated, so nerve exposure tracks
                 : PLASMA Cmax. This is what makes SC vs IV separate.
KD_BTZ   : 150   : DRG 20S proteasome Kd (ng/mL) - unsaturated in nerve
BTZ_J50  : 0.35  : DRG occupancy giving half-maximal injury flux
BTZ_JH   : 1.5113: Hill exponent of the injury flux (>1 => CONVEX => Cmax matters)
KEO_ACU  : 24    : Fast oxaliplatin channel-site equilibration rate (1/day)

// --------------------------------------------------- platinum genotoxic arm
KADD     : 0.050 : Pt-DNA adduct formation rate (1/day)
KREP     : 0.050 : Nucleotide-excision repair rate (1/day)
KAD50    : 0.45  : Adduct level giving half-maximal nucleolar stress
KAD_H    : 1.5   : Hill exponent of nucleolar stress
KPSYN_R  : 0.030 : Somal protein-synthesis recovery rate (1/day)
KPSYN_D  : 0.120 : Somal protein-synthesis damage rate (1/day)
KNEURON  : 0.0060: DRG neuron apoptosis rate (1/day) - IRREVERSIBLE
KN50     : 0.62  : Adduct level giving half-maximal neuron apoptosis

// -------------------------------------------- transport / mitochondria / ROS
KAT_REC  : 0.100 : Axonal transport recovery rate (1/day)
KAT_DAM  : 0.45  : Axonal transport damage rate (1/day)
KAT_PSYN : 0.060 : Transport damage from loss of somal supply (1/day)
KM_REC   : 0.150 : Mitochondrial recovery rate (1/day)
KM_DAM   : 0.20  : Mitochondrial damage rate (1/day)
KM_ROS   : 0.020 : ROS contribution to mitochondrial damage
KROS     : 0.55  : ROS generation from mitochondrial dysfunction (1/day)
KROS_MAC : 0.25  : ROS generation from activated macrophages (1/day)
KROS_CL  : 1.10  : ROS clearance rate (1/day)

// ------------------------------------------------- SARM1 axon-death program
// RISK is the FRACTION of the sensory axon population whose metabolic demand
// exceeds supply. Axons differ in length and therefore demand, so this is a
// smooth sigmoid in ENERGY/E_THR: a single hard threshold makes the virtual
// population bimodal (untouched or devastated) and the observed graded
// distribution of CTCAE grades cannot be reproduced.
E_THR    : 0.62  : Axonal energy adequacy threshold
RISK_OFF : 0.75  : Offset making RISK exactly zero in a healthy axon
RISK_K   : 0.80  : Half-maximal risk deficit
RISK_H   : 1.0   : Hill exponent of the risk function (Michaelis - deliberately
         :       : non-cooperative, so the transition band is broad)
KS_ON    : 0.050 : SARM1 activation rate (1/day)
KS_OFF   : 0.0301: SARM1 deactivation rate (1/day) - t1/2 = 23 d => COASTING

// --------------------------------------------------- axon loss / regeneration
// PHI_RELIEF: dying-back is SELF-LIMITING - losing distal axon reduces the
// transport and bioenergetic burden on the survivors, so ENERGY recovers as
// AXON falls. Without it, severity is all-or-nothing rather than graded.
KDEG     : 0.030 : Axon degeneration rate (1/day per unit SARM1)
KREGEN   : 0.0140: Axon regeneration rate (1/day)
PHI_RELIEF: 0.5009: Energy relief per unit fractional axon loss (FITTED)
KREGEN_G : 0.080 : SARM1 activity giving half-maximal block of regeneration

// ---------------------------------------------- neuroinflammation, biomarkers
// The ROS -> macrophage -> ROS loop gain must stay below 1 or the model
// diverges: gain = (KROS_MAC/KROS_CL) x (0.30 x KMAC/KMACD) = 0.12 here.
KMAC     : 0.60  : DRG macrophage recruitment rate (1/day)
KMACD    : 0.35  : DRG macrophage resolution rate (1/day)
KIL      : 0.45  : IL-1beta production rate (1/day)
KILD     : 0.35  : IL-1beta elimination rate (1/day)
KNFL     : 580   : NfL release per unit degeneration flux (pg/mL per day) -
         :       : scaled so the peak rise matches the 3-5x increase reported in
         :       : CIPN (Karteri 2022, Huehnchen 2022). NfL feeds nothing, so this
         :       : is a pure output scaling and changes no other prediction.
NFL_BASE : 0.495 : Baseline NfL production (pg/mL/day) - gives 10 pg/mL at rest
KELNFL   : 0.0495: NfL elimination rate (1/day) - t1/2 = 14 d

// --------------------------------------------- excitability, central gain
KEX      : 0.075 : Chronic hyperexcitability generation rate (1/day)
KEXD     : 0.050 : Chronic hyperexcitability decay rate (1/day)
KEX_IL   : 0.40  : IL-1beta contribution to hyperexcitability
KCOLD    : 1.30  : Acute cold-allodynia generation rate (1/day per mg/L)
KCOLDD   : 0.693 : Acute cold-allodynia decay rate (1/day) - t1/2 = 1 d
KCOLD_S  : 1.60  : Sensitization of the acute response by accumulated axon loss
KCS      : 0.55  : Central sensitization generation rate (1/day)
KCSD     : 0.16  : Central sensitization decay rate (1/day)
KCS_DUL  : 0.3535: Duloxetine-driven suppression of central gain (1/day) (FITTED)
KCS_PGB  : 0.40  : Pregabalin-driven suppression of central gain (1/day)
KNAT     : 1.10  : Spinal noradrenergic tone build-up rate (1/day)
KNATD    : 0.50  : Spinal noradrenergic tone decay rate (1/day)
EC50_DUL : 25    : Duloxetine EC50 for noradrenergic tone (ng/mL)
EC50_PGB : 3.0   : Pregabalin EC50 for alpha2-delta occupancy (mg/L)

// -------------------------------------------------- clinical score mapping
W_NEURO  : 0.80  : Weight of irreversible DRG neuron loss in the deficit
W_STRUCT : 0.70  : Structural weight in the CTCAE severity index
W_SYMPT  : 0.30  : Symptom weight in the CTCAE severity index
// Grade thresholds are fixed A PRIORI on fractional sensory-fibre loss
// (10 / 25 / 45%), not fitted: sural-nerve and IENFD series place the
// symptomatic threshold near 20-30% fibre loss. Keeping them fixed is what
// makes the non-fitted arms genuine out-of-sample predictions.
G1       : 0.100 : CTCAE grade 1 threshold
G2       : 0.250 : CTCAE grade 2 threshold
G3       : 0.450 : CTCAE grade 3 threshold
PAIN_WC  : 0.55  : Central-sensitization weight in perceived pain
PAIN_WE  : 0.30  : Peripheral-excitability weight in perceived pain
PAIN_WD  : 0.15  : Acute cold-allodynia weight in perceived pain
PAIN_SC  : 11.5  : Scale factor mapping to BPI 0-10

// ----------------------------- FITTED susceptibility scalers (see header)
KDAM_PT  : 0.3799: Platinum damage scaler   (fitted: IDEA FOLFOX 6 mo 47.7%)
KDAM_TAX : 0.1626: Taxane damage scaler     (fitted: ECOG 1199 weekly 27%)
KDAM_BTZ : 147.602: Bortezomib damage scaler (fitted: MMY-3021 IV 41%)
SIGMA_S  : 0.3664: Population log-SD of susceptibility (fitted: MOSAIC g>=3 12.4%)

// ----------------------------------- host covariates (1.0 = reference patient)
S        : 1.0   : Individual susceptibility multiplier (lognormal, median 1)
AXON0    : 100   : Baseline distal axon density (% of healthy-young norm)
RESERVE  : 1.0   : Bioenergetic reserve (age, diabetes, fitness)
REGEN    : 1.0   : Regenerative capacity (age, diabetes)
REPAIR   : 1.0   : Nucleotide-excision repair capacity (ERCC1)
OCT2     : 1.0   : DRG platinum uptake capacity (SLC22A2)
CRYO     : 0.0   : Fraction of distal delivery blocked by cryotherapy/compression

// ---------------------------------------------------------- oncology mapping
// Saturating log-kill jointly calibrated to IDEA 3-year DFS in BOTH risk
// strata. ONC_N0 is the residual clonogen number with NO oxaliplatin.
ONC_ALPHA: 0.449 : Oxaliplatin log-kill coefficient
ONC_C50  : 250   : Cumulative dose giving half-maximal log-kill (mg/m2)
ONC_N0   : 0.6315: Residual clonogens without oxaliplatin (high-risk stage III)

$CMT @annotated
A1_OXA  : Oxaliplatin ultrafilterable Pt, central (mg)
A2_OXA  : Oxaliplatin shallow peripheral (mg)
A3_OXA  : Oxaliplatin deep tissue-bound Pt (mg)
A1_PAC  : Paclitaxel central (mg)
A2_PAC  : Paclitaxel shallow peripheral (mg)
A3_PAC  : Paclitaxel deep peripheral (mg)
AB_BTZ  : Bortezomib subcutaneous depot (mg)
A1_BTZ  : Bortezomib central (mg)
A2_BTZ  : Bortezomib peripheral (mg)
AB_DUL  : Duloxetine gut depot (mg)
A1_DUL  : Duloxetine central (mg)
A1_PGB  : Pregabalin central (mg)
CE_PT   : DRG platinum biophase (accumulating, a.u.)
CE_TAX  : Taxane effect site (mg/L)
CE_BTZ  : Bortezomib DRG effect site (ng/mL)
CE_ACU  : Fast oxaliplatin channel effect site (mg/L)
ADDUCT  : Pt-DNA adduct burden (a.u.)
PSYN    : Somal protein-synthesis capacity (0-1)
MITO    : Mitochondrial functional capacity (0-1)
ROS     : Axonal oxidative stress (a.u.)
ATRANS  : Axonal transport capacity (0-1)
SARM    : SARM1 axon-death program activity (0-1)
AXON    : Distal sensory axon density (% of healthy norm)
NEURON  : Surviving DRG sensory neurons (% of baseline)
MAC     : DRG macrophage / neuroinflammation index (a.u.)
IL1B    : DRG IL-1beta (a.u.)
NFL     : Plasma neurofilament light chain (pg/mL)
EXCITC  : Chronic axonal hyperexcitability index (a.u.)
COLDA   : Acute cold-allodynia state (0-1)
CENTS   : Central sensitization gain (a.u.)
NATONE  : Spinal noradrenergic tone fraction (0-1)
CUMPT   : Cumulative oxaliplatin dose (mg/m2)
CUMTAX  : Cumulative taxane dose (mg/m2)
TCTHR   : Time paclitaxel above threshold (days)

$MAIN
PSYN_0   = 1.0;
MITO_0   = 1.0;
ATRANS_0 = 1.0;
AXON_0   = AXON0;
NEURON_0 = 100.0;
NFL_0    = NFL_BASE / KELNFL;   // 10 pg/mL

$ODE
// ============================== PK ==========================================
double c1   = A1_OXA / OXA_V1;                 // oxaliplatin UF Pt (mg/L)
double c2   = A2_OXA / OXA_V2;
double c3   = A3_OXA / OXA_V3;
double f12  = OXA_Q2 * (c1 - c2);
double f13  = OXA_Q3 * (c1 - c3);
dxdt_A1_OXA = - OXA_CL * c1 - f12 - f13;
dxdt_A2_OXA = f12;
dxdt_A3_OXA = f13;

double p1   = A1_PAC / PAC_V1;                 // paclitaxel (mg/L)
double q2p  = A2_PAC / PAC_V2;
double q3p  = A3_PAC / PAC_V3;
double g12  = PAC_Q2 * (p1 - q2p);
double g13  = PAC_Q3 * (p1 - q3p);
double elim = PAC_VMAX * p1 / (PAC_KM + p1);   // Michaelis-Menten
dxdt_A1_PAC = - elim - g12 - g13;
dxdt_A2_PAC = g12;
dxdt_A3_PAC = g13;

double bc1  = A1_BTZ / BTZ_V1 * 1000.0;        // bortezomib (ng/mL)
double ka_b = BTZ_KA * AB_BTZ;
double h12  = BTZ_Q * (A1_BTZ / BTZ_V1 - A2_BTZ / BTZ_V2);
dxdt_AB_BTZ = - ka_b;
dxdt_A1_BTZ = ka_b - BTZ_CL * A1_BTZ / BTZ_V1 - h12;
dxdt_A2_BTZ = h12;

dxdt_AB_DUL = - DUL_KA * AB_DUL;
dxdt_A1_DUL = DUL_KA * AB_DUL - DUL_CL * A1_DUL / DUL_V;
double cdul = A1_DUL / DUL_V * 1000.0;         // ng/mL
dxdt_A1_PGB = - PGB_CL * A1_PGB / PGB_V;
double cpgb = A1_PGB / PGB_V;                  // mg/L

// ========================= BIOPHASE / EFFECT SITES ==========================
// CRYO acts here, not on plasma: limb cooling / compression reduces the drug
// delivered to the DISTAL nerve without changing systemic exposure - which is
// exactly why it can be neuroprotective without being tumour-protective.
double delivery = 1.0 - CRYO;
dxdt_CE_PT  = KIN_PT * c1 * OCT2 * delivery - KOUT_PT * CE_PT;
dxdt_CE_TAX = KIN_TAX * p1 * delivery - KOUT_TAX * CE_TAX;
dxdt_CE_BTZ = KEO_BTZ  * (bc1 * delivery - CE_BTZ);
dxdt_CE_ACU = KEO_ACU  * (c1 - CE_ACU);

double TOCC   = CE_TAX / (CE_TAX + KD_TAX);          // tubulin occupancy
double PI_DRG = CE_BTZ / (CE_BTZ + KD_BTZ);          // DRG 20S occupancy
double J_BTZ  = HILL(PI_DRG, BTZ_J50, BTZ_JH);       // CONVEX injury flux

// The susceptibility multiplier S enters each drug chain EXACTLY ONCE:
// platinum at adduct formation, taxane/bortezomib at the insult signal.
double kpt   = KDAM_PT  * S;
double I_TAX = KDAM_TAX * S * TOCC;
double I_BTZ = KDAM_BTZ * S * J_BTZ;

// ======================= PLATINUM GENOTOXIC ARM =============================
dxdt_ADDUCT = KADD * kpt * CE_PT - KREP * REPAIR * ADDUCT;
double I_PT = HILL(ADDUCT, KAD50, KAD_H);            // nucleolar stress
dxdt_PSYN   = KPSYN_R * (1.0 - PSYN) - KPSYN_D * I_PT * PSYN;

// ========================= AXONAL TRANSPORT =================================
dxdt_ATRANS = KAT_REC * (1.0 - ATRANS)
            - KAT_DAM * (I_TAX + 0.45 * I_BTZ) * ATRANS
            - KAT_PSYN * (1.0 - PSYN) * ATRANS;

// ====================== MITOCHONDRIA & OXIDATIVE STRESS =====================
double mito_ins = 0.85 * I_PT + 0.55 * I_TAX + 0.45 * I_BTZ + KM_ROS * ROS;
dxdt_MITO = KM_REC * (1.0 - MITO) - KM_DAM * mito_ins * MITO;
dxdt_ROS  = KROS * (1.0 - MITO) / (RESERVE > 0.1 ? RESERVE : 0.1)
          + KROS_MAC * MAC - KROS_CL * ROS;

// ================= ENERGY, RISK AND THE SARM1 COMMITMENT ====================
double axl = (AXON0 - AXON) / AXON0;
if (axl < 0.0) axl = 0.0;
if (axl > 1.0) axl = 1.0;
double ENERGY = MITO * (0.25 + 0.75 * ATRANS) * RESERVE
              * (1.0 + PHI_RELIEF * axl);
double ratio  = E_THR / (ENERGY > 0.02 ? ENERGY : 0.02);
double RISK   = HILL((ratio - RISK_OFF > 0.0 ? ratio - RISK_OFF : 0.0),
                     RISK_K, RISK_H);
dxdt_SARM = KS_ON * RISK * (1.0 - SARM) - KS_OFF * SARM;

// ========================== AXON & NEURON ===================================
double axmax   = NEURON * AXON0 / 100.0;
double degflux = KDEG * SARM * AXON;
// An axon cannot regenerate while its own degeneration program is executing.
double regen_gate = 1.0 - HILL(SARM, KREGEN_G, 2.0);
dxdt_AXON   = - degflux
            + KREGEN * REGEN * regen_gate * (axmax > AXON ? axmax - AXON : 0.0);
dxdt_NEURON = - KNEURON * HILL(ADDUCT, KN50, 3.0) * NEURON;

// ==================== NEUROINFLAMMATION & BIOMARKERS ========================
dxdt_MAC  = KMAC * (degflux / 100.0 + 0.30 * ROS) - KMACD * MAC;
dxdt_IL1B = KIL * MAC / (MAC + 1.0) - KILD * IL1B;
dxdt_NFL  = KNFL * degflux / 100.0 + NFL_BASE - KELNFL * NFL;

// ===================== EXCITABILITY (chronic and acute) =====================
dxdt_EXCITC = KEX * (axl + KEX_IL * IL1B) - KEXD * EXCITC;
dxdt_COLDA  = KCOLD * CE_ACU * (1.0 + KCOLD_S * axl) * (1.0 - COLDA)
            - KCOLDD * COLDA;

// ========================= CENTRAL SENSITIZATION ============================
double A2D     = cpgb / (cpgb + EC50_PGB);
double barrage = EXCITC + 2.0 * COLDA;
dxdt_CENTS  = KCS * barrage
            - (KCSD + KCS_DUL * NATONE + KCS_PGB * A2D) * CENTS;
dxdt_NATONE = KNAT * (cdul / (cdul + EC50_DUL)) * (1.0 - NATONE)
            - KNATD * NATONE;

// ============================ CUMULATIVE TRACKERS ===========================
dxdt_CUMPT  = 0.0;    // incremented by dosing records via a companion cmt
dxdt_CUMTAX = 0.0;
dxdt_TCTHR  = (p1 > PAC_THR) ? 1.0 : 0.0;

$TABLE
// ----------------------------- structural deficit ---------------------------
// Structural deficit is LINEAR in fibre loss, measured FROM THE PATIENT'S OWN
// BASELINE (AXON0), plus an explicit neuronopathy term (a lost soma can never be
// reinnervated). CTCAE grades TREATMENT-EMERGENT neuropathy, so a diabetic's
// pre-existing subclinical deficit is not already graded CIPN - reduced reserve
// must act through bioenergetics (RESERVE) and regeneration (REGEN), not by
// moving the patient closer to a grading threshold.
double SDEF = (AXON0 - AXON) / 100.0 + W_NEURO * (100.0 - NEURON) / 100.0;
if (SDEF < 0.0) SDEF = 0.0;
if (SDEF > 1.0) SDEF = 1.0;

// --------------------------------- symptoms ---------------------------------
double cn = CENTS  / (CENTS  + 1.6);
double en = EXCITC / (EXCITC + 1.6);
double painraw = PAIN_WC * cn + PAIN_WE * en + PAIN_WD * COLDA;
double BPI = PAIN_SC * painraw;  if (BPI > 10.0) BPI = 10.0;

// CTCAE is graded at pre-dose visits, i.e. AFTER the previous cycle acute
// oxaliplatin syndrome has resolved, so the severity index is CHRONIC only.
double SYMPT = 0.62 * cn + 0.28 * en;  if (SYMPT > 1.0) SYMPT = 1.0;
double CS = W_STRUCT * SDEF + W_SYMPT * SYMPT;

// CIPN20 is patient-reported and DOES capture the acute syndrome.
double CIPN20 = 100.0 * (0.62 * SDEF + 0.26 * painraw / 0.60 + 0.12 * COLDA);
if (CIPN20 > 100.0) CIPN20 = 100.0;

double GRADE = 0.0;
if (CS >= G1) GRADE = 1.0;
if (CS >= G2) GRADE = 2.0;
if (CS >= G3) GRADE = 3.0;

double IENFD = 7.0 * AXON / 100.0;              // fibres/mm (7.0 = normal)
double TNSc  = 28.0 * (0.60 * SDEF + 0.25 * en + 0.15 * (1.0 - AXON / 100.0));

// --------- 3-year DFS from cumulative oxaliplatin (saturating log-kill) ------
double kfrac = CUMPT / (CUMPT + ONC_C50);
double NRES  = ONC_N0 * exp(-ONC_ALPHA * kfrac);
double DFS3  = exp(-NRES);

// ---------------------------- exposed PK/PD signals -------------------------
double C_OXA  = A1_OXA / OXA_V1;
double C_PAC  = A1_PAC / PAC_V1;
double C_BTZ  = A1_BTZ / BTZ_V1 * 1000.0;
double C_DUL  = A1_DUL / DUL_V * 1000.0;
double TOCCo   = CE_TAX / (CE_TAX + KD_TAX);
double PI_DRGo = CE_BTZ / (CE_BTZ + KD_BTZ);
double ENERGYo = MITO * (0.25 + 0.75 * ATRANS) * RESERVE
               * (1.0 + PHI_RELIEF * ((AXON0 - AXON) / AXON0));

$CAPTURE @annotated
CS      : CTCAE severity index (0-1)
GRADE   : CTCAE sensory neuropathy grade (0-3)
CIPN20  : EORTC QLQ-CIPN20 sensory subscale (0-100)
BPI     : Average pain, Brief Pain Inventory (0-10)
SDEF    : Structural sensory deficit (0-1)
IENFD   : Intraepidermal nerve fibre density (fibres/mm)
TNSc    : Total Neuropathy Score clinical (0-28)
DFS3    : Predicted 3-year disease-free survival
C_OXA   : Oxaliplatin ultrafilterable Pt (mg/L)
C_PAC   : Paclitaxel plasma concentration (mg/L)
C_BTZ   : Bortezomib plasma concentration (ng/mL)
C_DUL   : Duloxetine plasma concentration (ng/mL)
TOCCo   : Tubulin occupancy (0-1)
PI_DRGo : DRG 20S proteasome occupancy (0-1)
ENERGYo : Axonal energy adequacy
'

mod <- mcode("cipn", cipn_code, end = 400, delta = 1)

## =============================================================================
## 2.  DOSING HELPERS
##
## NOTE ON CUMPT/CUMTAX: mrgsolve cannot increment a compartment from inside
## $ODE using the dosing rate, so cumulative dose is supplied as an explicit
## bolus into the CUMPT / CUMTAX compartments alongside each drug record. That
## keeps the tracker exactly equal to the delivered mg/m2.
## =============================================================================

BSA_DEF <- 1.8

ev_oxa <- function(n, dose = 85, ii = 14, bsa = BSA_DEF, start = 0) {
  drug <- ev(amt = dose * bsa, cmt = "A1_OXA", time = start,
             ii = ii, addl = n - 1, tinf = 2 / 24)
  trk  <- ev(amt = dose, cmt = "CUMPT", time = start, ii = ii, addl = n - 1)
  c(drug, trk)
}

ev_pac <- function(n, dose = 80, ii = 7, hours = 1, bsa = BSA_DEF, start = 0) {
  drug <- ev(amt = dose * bsa, cmt = "A1_PAC", time = start,
             ii = ii, addl = n - 1, tinf = hours / 24)
  trk  <- ev(amt = dose, cmt = "CUMTAX", time = start, ii = ii, addl = n - 1)
  c(drug, trk)
}

## Bortezomib d1,4,8,11 of a 21-day cycle (or weekly d1,8,15,22 of 35 days)
ev_btz <- function(cycles = 8, route = c("IV", "SC"), dose = 1.3,
                   bsa = BSA_DEF, weekly = FALSE) {
  route <- match.arg(route)
  cmt   <- if (route == "IV") "A1_BTZ" else "AB_BTZ"
  days  <- if (weekly) c(0, 7, 14, 21) else c(0, 3, 7, 10)
  cyc   <- if (weekly) 35 else 21
  out   <- NULL
  for (c_i in seq_len(cycles)) for (d in days) {
    e <- ev(amt = dose * bsa, cmt = cmt, time = (c_i - 1) * cyc + d)
    out <- if (is.null(out)) e else c(out, e)
  }
  out
}

## Oral daily dosing of a symptomatic drug over [t0, t1]
ev_oral <- function(cmt, t0, t1, mg) {
  n <- max(1, floor(t1 - t0))
  ev(amt = mg, cmt = cmt, time = t0, ii = 1, addl = n - 1)
}

ev_duloxetine <- function(t0, t1, mg = 60) ev_oral("AB_DUL", t0, t1, mg)
ev_pregabalin <- function(t0, t1, mg = 300) ev_oral("A1_PGB", t0, t1, mg)

## =============================================================================
## 3.  EIGHTEEN TREATMENT SCENARIOS
## =============================================================================

scenarios <- list(

  ## ---- platinum: the IDEA question -----------------------------------------
  S01 = list(label = "FOLFOX 6 months (12 x 85 mg/m2 q2w, 1020 mg/m2)",
             ev = ev_oxa(12), end = 500, par = list()),
  S02 = list(label = "FOLFOX 3 months (6 x 85 mg/m2 q2w, 510 mg/m2)",
             ev = ev_oxa(6), end = 400, par = list()),
  S03 = list(label = "CAPOX 3 months (4 x 130 mg/m2 q3w, 520 mg/m2)",
             ev = ev_oxa(4, dose = 130, ii = 21), end = 400, par = list()),
  S04 = list(label = "CAPOX 6 months (8 x 130 mg/m2 q3w, 1040 mg/m2)",
             ev = ev_oxa(8, dose = 130, ii = 21), end = 500, par = list()),

  ## ---- schedule at MATCHED cumulative dose ---------------------------------
  S05 = list(label = "FOLFOX 1020 mg/m2 spread q3w (12 x 85, 9 months)",
             ev = ev_oxa(12, ii = 21), end = 600, par = list()),
  S06 = list(label = "FOLFOX stop-and-go / OPTIMOX (6 on, 8 wk off, 6 on)",
             ev = c(ev_oxa(6), ev_oxa(6, start = 140)), end = 600,
             par = list()),
  S07 = list(label = "FOLFOX 6 months, 20% dose reduction (12 x 68)",
             ev = ev_oxa(12, dose = 68), end = 500, par = list()),

  ## ---- taxanes -------------------------------------------------------------
  S08 = list(label = "Paclitaxel weekly 80 mg/m2 x 12 (960 mg/m2)",
             ev = ev_pac(12), end = 400, par = list()),
  S09 = list(label = "Paclitaxel q3w 175 mg/m2 x 4 (700 mg/m2)",
             ev = ev_pac(4, dose = 175, ii = 21, hours = 3), end = 400,
             par = list()),

  ## ---- bortezomib route: same AUC, 1/11 the Cmax ---------------------------
  S10 = list(label = "Bortezomib 1.3 mg/m2 IV bolus, d1/4/8/11 x 8 cycles",
             ev = ev_btz(8, "IV"), end = 500, par = list()),
  S11 = list(label = "Bortezomib 1.3 mg/m2 SUBCUTANEOUS, d1/4/8/11 x 8",
             ev = ev_btz(8, "SC"), end = 500, par = list()),
  S12 = list(label = "Bortezomib 1.3 mg/m2 IV WEEKLY d1/8/15/22 x 8",
             ev = ev_btz(8, "IV", weekly = TRUE), end = 600, par = list()),

  ## ---- prevention ----------------------------------------------------------
  S13 = list(label = "FOLFOX 6 months + cryotherapy every cycle (45% block)",
             ev = ev_oxa(12), end = 500, par = list(CRYO = 0.45)),
  S14 = list(label = "FOLFOX 6 months + compression therapy (25% block)",
             ev = ev_oxa(12), end = 500, par = list(CRYO = 0.25)),

  ## ---- host phenotypes at identical exposure -------------------------------
  S15 = list(label = "FOLFOX 6 months, DIABETIC host (reduced reserve)",
             ev = ev_oxa(12), end = 500,
             par = list(AXON0 = 82, RESERVE = 0.90, REGEN = 0.70)),
  S16 = list(label = "FOLFOX 6 months, age 75 (reduced regeneration)",
             ev = ev_oxa(12), end = 500,
             par = list(REGEN = 0.70, RESERVE = 0.94)),

  ## ---- symptomatic treatment of established CIPN ---------------------------
  S17 = list(label = "FOLFOX 6 mo, then duloxetine 60 mg d250-285",
             ev = c(ev_oxa(12), ev_duloxetine(250, 285)), end = 400,
             par = list()),
  S18 = list(label = "FOLFOX 6 mo, then pregabalin 300 mg d250-285",
             ev = c(ev_oxa(12), ev_pregabalin(250, 285)), end = 400,
             par = list())
)

run_scenario <- function(s, ...) {
  m <- mod
  if (length(s$par)) m <- param(m, s$par)
  mrgsim(m, events = s$ev, end = s$end, delta = 1, ...) %>%
    as_tibble() %>%
    mutate(scenario = s$label)
}

run_all <- function() bind_rows(lapply(scenarios, run_scenario))

## Summary table: peak severity, coasting lag, nadir IENFD, predicted DFS
summarise_scenarios <- function(sims = run_all()) {
  sims %>%
    group_by(scenario) %>%
    summarise(
      peak_CS      = max(CS),
      day_peak_CS  = time[which.max(CS)],
      peak_CIPN20  = max(CIPN20),
      peak_BPI     = max(BPI),
      nadir_IENFD  = min(IENFD),
      day_nadir    = time[which.min(IENFD)],
      cum_Pt       = max(CUMPT),
      cum_taxane   = max(CUMTAX),
      DFS3         = last(DFS3),
      CIPN20_d365  = CIPN20[which.min(abs(time - 365))],
      .groups = "drop"
    ) %>%
    arrange(desc(peak_CS))
}

## =============================================================================
## 4.  POPULATION INCIDENCE
##
## Peak severity is MONOTONE in the susceptibility multiplier S, so
##     P(peak grade >= g) = P(S >= s*_g)
## can be obtained exactly by bisecting for s*_g and evaluating the lognormal
## survival function - no Monte-Carlo noise, ~30 simulations per estimate.
## =============================================================================

peak_CS_for_S <- function(S, evnt, end, extra = list()) {
  m <- param(mod, c(list(S = S), extra))
  max(mrgsim(m, events = evnt, end = end, delta = 1)$CS)
}

incidence <- function(evnt, end, grade = 2, extra = list()) {
  thr <- unname(param(mod)[[paste0("G", grade)]])
  lo <- 1e-4; hi <- 1e4
  if (peak_CS_for_S(lo, evnt, end, extra) > thr) return(1)
  if (peak_CS_for_S(hi, evnt, end, extra) < thr) return(0)
  for (i in 1:32) {
    mid <- sqrt(lo * hi)
    if (peak_CS_for_S(mid, evnt, end, extra) < thr) lo <- mid else hi <- mid
  }
  sigma <- unname(param(mod)$SIGMA_S)
  1 - pnorm(log(sqrt(lo * hi)) / sigma)
}

## Reproduce the calibration + validation table
validation_table <- function() {
  rows <- list(
    list("FOLFOX 6 months", ev_oxa(12), 500, c(47.7, 12.4), "FIT"),
    list("FOLFOX 3 months", ev_oxa(6), 400, c(16.6, 2.7), "PREDICTION"),
    list("Paclitaxel weekly x12", ev_pac(12), 400, c(27, NA), "FIT"),
    list("Paclitaxel q3w x4", ev_pac(4, dose = 175, ii = 21, hours = 3),
         400, c(20, NA), "PREDICTION"),
    list("Bortezomib IV", ev_btz(8, "IV"), 500, c(41, NA), "FIT"),
    list("Bortezomib SC", ev_btz(8, "SC"), 500, c(24, NA), "PREDICTION")
  )
  bind_rows(lapply(rows, function(r) tibble(
    arm        = r[[1]],
    model_g2   = 100 * incidence(r[[2]], r[[3]], 2),
    model_g3   = 100 * incidence(r[[2]], r[[3]], 3),
    observed_g2 = r[[4]][1],
    observed_g3 = r[[4]][2],
    status     = r[[5]]
  )))
}

## Multivariate virtual population (for the covariate analysis only)
virtual_population <- function(n = 300, evnt = ev_oxa(12), end = 300,
                               seed = 20260729) {
  set.seed(seed)
  sig  <- unname(param(mod)$SIGMA_S)
  diab <- runif(n) < 0.20
  age  <- runif(n, 45, 75)
  idata <- tibble(
    ID      = seq_len(n),
    S       = exp(rnorm(n, 0, sig)),
    AXON0   = pmax(60, 100 - ifelse(diab, 18, 0) - 0.18 * pmax(0, age - 50)),
    RESERVE = pmax(0.6, 1 - ifelse(diab, 0.18, 0) - 0.003 * pmax(0, age - 50)),
    REGEN   = pmax(0.35, 1 - ifelse(diab, 0.30, 0) - 0.010 * pmax(0, age - 50))
  )
  sim <- mrgsim(mod, events = evnt, idata = idata, end = end, delta = 2) %>%
    as_tibble()
  peak <- sim %>% group_by(ID) %>%
    summarise(peak_grade = max(GRADE), peak_CIPN20 = max(CIPN20),
              nadir_IENFD = min(IENFD), .groups = "drop") %>%
    left_join(idata %>% mutate(diabetes = diab, age = age), by = "ID")
  list(sim = sim, subjects = peak,
       incidence = peak %>% summarise(
         g2 = 100 * mean(peak_grade >= 2), g3 = 100 * mean(peak_grade >= 3)))
}

## =============================================================================
## 5.  THE THERAPEUTIC-INDEX OPTIMISER
##
## Net utility = (3-year DFS gain over no oxaliplatin, percentage points)
##               - w * (grade >=2 CIPN incidence, percentage points)
## The CIPN cost curve is identical in every risk stratum; only the DFS gain
## available differs. The optimum number of cycles therefore depends on RISK,
## not on neurotoxicity - which is the IDEA recommendation, derived rather
## than assumed.
## =============================================================================

therapeutic_index <- function(cycles = c(2, 4, 6, 8, 10, 12, 14), w = 0.15,
                              N0 = c(high_risk = 0.6315, low_risk = 0.2623)) {
  dfs <- function(cum, n0) {
    a <- unname(param(mod)$ONC_ALPHA); c50 <- unname(param(mod)$ONC_C50)
    100 * exp(-n0 * exp(-a * cum / (cum + c50)))
  }
  bind_rows(lapply(cycles, function(nc) {
    i2 <- 100 * incidence(ev_oxa(nc), 14 * nc + 250, 2)
    cum <- nc * 85
    tibble(cycles = nc, cum_Pt = cum, grade2_pct = i2,
           DFS_high = dfs(cum, N0[["high_risk"]]),
           DFS_low  = dfs(cum, N0[["low_risk"]]),
           utility_high = dfs(cum, N0[["high_risk"]]) - dfs(0, N0[["high_risk"]])
                          - w * i2,
           utility_low  = dfs(cum, N0[["low_risk"]]) - dfs(0, N0[["low_risk"]])
                          - w * i2)
  }))
}

## =============================================================================
## 6.  PLOTS
## =============================================================================

plot_coasting <- function() {
  bind_rows(run_scenario(scenarios$S01), run_scenario(scenarios$S02)) %>%
    select(time, scenario, CIPN20, IENFD, SARM, AXON) %>%
    pivot_longer(c(CIPN20, IENFD, SARM, AXON)) %>%
    ggplot(aes(time, value, colour = scenario)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = c(154, 70), linetype = 3) +
    facet_wrap(~name, scales = "free_y") +
    labs(title = "Coasting: severity peaks weeks AFTER the last dose",
         subtitle = "dotted lines mark the last oxaliplatin infusion",
         x = "Day", y = NULL) +
    theme_bw() + theme(legend.position = "bottom")
}

plot_route <- function() {
  bind_rows(run_scenario(scenarios$S10), run_scenario(scenarios$S11)) %>%
    filter(time <= 21) %>%
    select(time, scenario, C_BTZ, PI_DRGo, CIPN20) %>%
    pivot_longer(c(C_BTZ, PI_DRGo, CIPN20)) %>%
    ggplot(aes(time, value, colour = scenario)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~name, scales = "free_y", ncol = 1) +
    labs(title = "Bortezomib SC vs IV: equal AUC, 1/11 the Cmax",
         x = "Day", y = NULL) +
    theme_bw() + theme(legend.position = "bottom")
}

plot_therapeutic_index <- function(ti = therapeutic_index()) {
  ti %>%
    select(cum_Pt, grade2_pct, utility_high, utility_low) %>%
    pivot_longer(-cum_Pt) %>%
    ggplot(aes(cum_Pt, value, colour = name)) +
    geom_line(linewidth = 0.9) + geom_point() +
    labs(title = "Therapeutic-index optimum in cumulative oxaliplatin dose",
         subtitle = paste("neurotoxicity integrates linearly;",
                          "anti-tumour benefit saturates by ~500 mg/m2"),
         x = "Cumulative oxaliplatin (mg/m2)", y = "percentage points") +
    theme_bw() + theme(legend.position = "bottom")
}

## =============================================================================
## 7.  EXAMPLE SESSION
## =============================================================================
if (interactive()) {
  sims <- run_all()
  print(summarise_scenarios(sims), n = 20)
  print(validation_table())
  print(therapeutic_index())
  vp <- virtual_population(300)
  print(vp$incidence)
  print(plot_coasting())
  print(plot_route())
  print(plot_therapeutic_index())
}
