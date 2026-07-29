# =============================================================================
#  hg_mrgsolve_model.R
#  Hyperemesis Gravidarum (HG) — QSP model for mrgsolve
# =============================================================================
#
#  42 ODE compartments.  Time is DAYS OF GESTATIONAL AGE by LMP; negative times
#  are pre-conception, which the model needs because its central claim is that
#  the therapeutic window for the emetic set-point closes at conception.
#
#  ---------------------------------------------------------------------------
#  WHAT THIS MODEL IS FOR
#  ---------------------------------------------------------------------------
#  It exists to explain five observations that are awkward together:
#
#   (A) Maternal GDF15 rises through the first trimester and stays high (indeed
#       keeps climbing) into the third, yet nausea and vomiting of pregnancy
#       peaks at GA 9-11 weeks and has essentially gone by 16-20 weeks.
#       [Marjono 2003 PMID 12495665; Fejzo 2019 PMID 31515515]
#   (B) Women with beta-thalassaemia have chronically HIGH GDF15 and report
#       very low levels of NVP.  [Fejzo 2024 PMID 38092039]
#   (C) LOW pre-pregnancy GDF15 is a RISK factor.  [same]
#   (D) Pre-pregnancy metformin -- which raises GDF15 and is stopped before the
#       nausea would start -- is associated with aRR 0.29 for HG; pre-pregnancy
#       tobacco, which raises GDF15 less, gives aRR 0.51.
#       [Sharma 2025 PMID 40588059]
#   (E) In the first placebo-controlled HG trial, ondansetron -- the
#       guideline-recommended second-line drug, at ~98% 5-HT3 occupancy --
#       moved PUQE-24 by only -0.51 (95% CI -2.32 to 1.30), while mirtazapine
#       moved it -1.86 (-3.61 to -0.12), separating further after day 4.
#       [Ostenfeld 2026 PMID 41478546]
#
#  ---------------------------------------------------------------------------
#  THREE STRUCTURAL COMMITMENTS
#  ---------------------------------------------------------------------------
#  1. THE DETECTOR IS A RATIO, NOT A LEVEL.
#     The GDF15 -> GFRAL -> area postrema axis carries an adapting set-point SP
#     that tracks log(GDF15) with a time constant TAU_SP of about four weeks.
#     Emetic drive is GDF15 / (adapted set-point), not GDF15.
#       - the placental ramp (doubling every ~4 days at GA 5-9 wk) outruns the
#         set-point, so the ratio spikes;
#       - once the ramp flattens the set-point catches up and the ratio decays
#         to ~1, so symptoms remit WHILE THE HORMONE STAYS HIGH  -> (A);
#       - anything that raises GDF15 BEFORE conception raises SP before the
#         placenta arrives, so the same absolute pregnancy level produces a
#         much smaller ratio -> (B), (C), (D).
#     ALPHA is "adaptation completeness".  ALPHA = 1 is a pure fold-change
#     detector; ALPHA = 0 is a pure level detector.  Setting ALPHA = 0 turns
#     this file into the naive comparator model, and every one of (A)-(D)
#     inverts.  That single switch is how the claim is shown to be load-bearing.
#
#  2. EFFICACY IS SET BY NODE POSITION.
#     Drive reaches the emetic pattern generator through the NTS, which sums a
#     LARGE GFRAL/area-postrema term (w_AP) and a SMALL peripheral vagal/5-HT3
#     term (W_VAG).  A drug cannot outperform the weight of the node it
#     occupies.  Drugs are therefore attached to SPECIFIC receptor classes at
#     specific nodes, with occupancy computed from published Ki, unbound
#     fraction, brain penetration and PK -- not to a generic "antiemetic
#     effect".  This predicts the ordering
#        ligand > alpha-2-delta > multi-receptor > single H1 > peripheral 5-HT3
#        > nothing at all (corticosteroids, a published null) -> (E).
#
#  3. THE NUTRITIONAL CASCADE RUNS ON A DIFFERENT CLOCK.
#     Vomiting drives volume/Cl-/K+ loss into a chloride-responsive alkalosis
#     that will not correct until chloride is replaced, and empties a ~28 mg
#     thiamine store with a ~15 day half-life.  So an antiemetic can normalise
#     PUQE while Wernicke risk is still climbing, and IV dextrose given before
#     thiamine makes it worse.  [Oudman 2019 PMID 30889425]
#
#  ---------------------------------------------------------------------------
#  CALIBRATION — WHAT IS FITTED AND WHAT IS NOT
#  ---------------------------------------------------------------------------
#  FITTED (4 parameters):
#    ALPHA, TAU_SP   adaptation completeness and time constant, to the natural
#                    history (peak GA week; remission by 16 weeks)
#    W_VAG           peripheral 5-HT3 limb share, to ondansetron's -0.51
#    E0              central node-authority scale, to mirtazapine's -1.86
#
#  PREDICTED (not fitted):
#    beta-thalassaemia protection; pre-conception metformin protection AND its
#    failure when started at GA 6 wk; the metformin < tobacco ordering;
#    doxylamine/pyridoxine effect size; gabapentin superiority over an active
#    comparator; clonidine effect size; the corticosteroid null; the post-day-4
#    widening of mirtazapine vs ondansetron (from PK alone, no onset
#    parameter); thyrotoxicosis on the fetal-production axis only; and the
#    harm of IV dextrose given without thiamine.
#
#  A STATED MISS: metoclopramide comes out essentially inactive, whereas trials
#  find it comparable to promethazine.  This is the cleanest place to try to
#  falsify commitment 2.  See hg_reference_impl.py section 5.
#
#  ---------------------------------------------------------------------------
#  VERIFICATION
#  ---------------------------------------------------------------------------
#  Every number quoted above is regenerated by the dependency-free Python twin
#  of this same ODE system:
#
#      python3 hg_reference_impl.py            # all scenarios
#      python3 hg_reference_impl.py --check    # PASS/FAIL (currently 15/15)
#
#  The twin needs only python3, so the calibration can be checked without R,
#  mrgsolve or a compiler.
#
#  UNITS
#    time          days of gestational age (LMP); negative = pre-conception
#    GDF15         pg/mL          hCG            IU/L
#    drug amounts  nmol           concentrations nM
#    volume        L              electrolytes   mmol (amounts), mmol/L (conc)
#    weight        kg             thiamine       mg
#    PUQE-24       3..15 (the real instrument's range)
#
#  DISCLAIMER: educational / research QSP model.  Semi-quantitative, not
#  independently validated, and NOT for clinical decision-making.  In
#  particular the GDF15-targeting scenarios are untested in pregnancy and the
#  model itself flags a fetal-safety concern for them (GDF15 supports
#  trophoblast invasion: PMID 37272232, 40157640).
# =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

hg_code <- '
$PARAM @annotated
// ---- gestational timing -----------------------------------------------
T_IMP     :  22.0  : Implantation, day GA by LMP (placental secretion starts)

// ---- placental trophoblast --------------------------------------------
KG_PLAC   :  0.100 : Logistic growth rate of trophoblast mass (1/day)
PLAC0     :  0.020 : Relative trophoblast mass seeded at implantation
PLAC_CAPG :  0.0025: Growth of carrying capacity, so mass rises to term (1/day)
TROPH_GAIN:  1.000 : Trophoblast secretory gain -- the FETAL axis of HG risk

// ---- trophoblast integrated stress response ---------------------------
ISR_HI    :  1.75  : ISR index while the placenta is hypoxic (pre-perfusion)
ISR_LO    :  1.00  : ISR index once intervillous flow is established
T_PERF1   :  66.0  : Start of the perfusion transition (day GA)
T_PERF2   : 102.0  : End of the perfusion transition (day GA)
TAU_ISR   :  4.0   : Lag of the ISR state behind placental pO2 (day)

// ---- maternal GDF15 ---------------------------------------------------
KEL_GDF   :  5.545 : GDF15 elimination rate (1/day; t1/2 ~3 h)
GDF_BASE  : 250.0  : Typical non-pregnant plasma GDF15 (pg/mL)
KS_GDF_PL : 48000  : Placental GDF15 synthesis per unit (mass x ISR)
GREF      : 250.0  : Population reference GDF15 for the level detector

// ---- THE FOLD-CHANGE DETECTOR (commitment 1) --------------------------
TAU_SP    : 30.0   : Adaptation time constant of the set-point (day) [FITTED]
ALPHA     :  0.92  : Adaptation completeness; 1 = ratio, 0 = level  [FITTED]

// ---- GFRAL receptor pool ----------------------------------------------
KSYN_GF   :  0.35  : GFRAL synthesis (1/day)
KDEG_GF   :  0.35  : GFRAL constitutive degradation (1/day)
KINT_GF   :  0.20  : Ligand-driven GFRAL internalisation (1/day)
KD_GF     : 5000   : GDF15 concentration half-maximal for internalisation

// ---- area postrema / NTS / emesis -------------------------------------
SENS      :  1.00  : Constitutional GFRAL gain -- the MATERNAL axis of risk
EC50_AP   :  4.10  : Effective drive at half-maximal AP output
HILL_AP   :  4.0   : AP is a threshold detector, not a proportional one
TAU_AP    :  0.30  : AP activation time constant (day)
TAU_NTS   :  0.25  : NTS integrator time constant (day)
TAU_NAUS  :  0.35  : Nausea state time constant (day)
TAU_CPG   :  0.20  : Emetic pattern generator time constant (day)
W_TOT     :  1.000 : Total NTS input weight (w_AP is the remainder)
W_VAG     :  0.034 : Peripheral 5-HT3 / vagal limb share          [FITTED]
W_VEST    :  0.055 : Vestibular / H1 limb share
W_HCG     :  0.030 : hCG-associated limb share
E0        :  0.044 : Central node-authority scale                 [FITTED]
R_H1      :  1.0   : Relative NTS authority of H1
R_5HT2    :  1.4   : Relative NTS authority of 5-HT2A/2C
R_5HT3C   :  0.6   : Relative NTS authority of central 5-HT3
R_M1      :  0.4   : Relative NTS authority of M1
R_D2      :  0.5   : Relative AP authority of D2
R_A2      :  5.8   : Relative NTS authority of presynaptic alpha-2
R_A2D     :  7.5   : Relative NTS authority of Cav alpha-2-delta
EC50_NAUS :  0.60  : NTS activity at half-maximal nausea hours
HILL_NAUS :  6.0   : Nausea Hill coefficient
EC50_VOM  :  0.66  : NTS activity at half-maximal vomiting (higher threshold)
HILL_VOM  :  7.0   : Vomiting Hill coefficient
VOM_MAX   : 12.0   : Emetic episodes per day at full CPG output
RETCH_RAT :  1.45  : Retches per vomit

// ---- gastric emptying / peripheral limb -------------------------------
TAU_GAS   :  1.0   : Gastric delay index time constant (day)
GAS_P4    :  0.35  : Progesterone-driven gastric delay
GAS_GDF   :  0.30  : GDF15-driven gastric delay
EC50_EC   :  0.55  : Gastric delay at half-maximal EC-cell 5-HT release

// ---- hCG --------------------------------------------------------------
KS_HCG    : 28000  : hCG synthesis per unit (mass x cytotrophoblast index)
KEL_HCG   :  0.28  : hCG elimination (1/day)
T_CTB1    : 60.0   : Cytotrophoblast index begins to fall (day GA)
T_CTB2    : 130.0  : Cytotrophoblast index reaches plateau (day GA)
CTB_PLAT  :  0.16  : Cytotrophoblast plateau fraction

// ---- fluid / electrolytes ---------------------------------------------
VOL0      : 14.0   : Extracellular fluid volume (L)
ORAL_FLUID:  2.30  : Oral fluid intake at normal eating (L/day)
WATER_FLR :  0.55  : Water intake retained at zero food intake (fraction)
THIRST_G  :  0.70  : Thirst amplification per unit RAAS activation
INSENS    :  0.90  : Insensible losses (L/day)
VOM_VOL   :  0.055 : Gastric juice lost per emetic episode (L)
KREN_VOL  :  0.55  : Renal volume regulation gain (1/day)
URINE_MIN :  0.40  : Minimum urine output (L/day)
GJ_CL     : 140.0  : Chloride in gastric juice (mmol/L)
GJ_H      :  75.0  : Hydrogen ion in gastric juice (mmol/L)
GJ_NA     :  42.0  : Sodium in gastric juice (mmol/L)
GJ_K      :  11.0  : Potassium in gastric juice (mmol/L)
NA_TOT0   : 1960   : Initial exchangeable sodium (mmol)
CL_TOT0   : 1442   : Initial exchangeable chloride (mmol)
HCO30     : 24.0   : Baseline plasma bicarbonate (mmol/L)
K_TOT0    : 3500   : Initial total body potassium (mmol)
K_INTAKE  : 72.0   : Dietary potassium at full intake (mmol/day)
K_URINE0  : 72.0   : Baseline renal potassium excretion (mmol/day)
K_ALDO    :  1.2   : Aldosterone amplification of renal K+ loss
K_RETAIN  :  0.86  : Fractional body K+ at which the kidney stops excreting
KREN_HCO3 :  0.28  : Renal bicarbonate excretion gain (1/day)
HCO3_ESC  :  0.30  : Chloride-independent fraction of HCO3- excretion
NA_HCO3_C :  0.45  : Fraction of excreted HCO3- obligated to carry Na+
NA_FLOOR  : 126.0  : Plasma Na+ at which the Na+/HCO3- coupling gives way
NA_INTAKE : 150.0  : Dietary sodium at full intake (mmol/day)
CL_INTAKE : 150.0  : Dietary chloride at full intake (mmol/day)

// ---- energy / weight --------------------------------------------------
WT0       : 65.0   : Pre-pregnancy weight (kg)
KCAL_NEED : 2100   : Energy requirement (kcal/day)
KCAL_FULL : 2100   : Energy intake at full oral intake (kcal/day)
KCAL_PERKG: 7700   : Energy per kg of tissue (kcal)
WT_FETAL  :  0.052 : Conceptus + uterus + plasma gain, second half (kg/day)

// ---- ketones / liver --------------------------------------------------
KET_MAX   :  5.2   : Beta-hydroxybutyrate at total starvation (mmol/L)
TAU_KET   :  0.5   : Ketone state time constant (day)
ALT_BASE  : 20.0   : Baseline ALT (U/L)
ALT_MAX   : 190.0  : ALT ceiling in severe starvation (U/L)
TAU_ALT   :  3.0   : ALT time constant (day)

// ---- thiamine ---------------------------------------------------------
THI0      : 28.0   : Whole-body thiamine store (mg)
THI_DIET  :  1.45  : Thiamine absorbed at full oral intake (mg/day)
KEL_THI   :  0.046 : Thiamine elimination (1/day; t1/2 ~15 d)
THI_CRIT  : 12.0   : Store below which Wernicke hazard starts (mg)
HZ_WE     :  0.060 : Wernicke hazard scale (1/day)
THI_PERCHO:  0.0022: Thiamine consumed per gram of glucose infused (mg/g)
THI_REFEED:  2.6   : Amplification of that when the store is already low

// ---- thyroid ----------------------------------------------------------
FT40      : 14.0   : Baseline free T4 (pmol/L)
TSH0      :  1.6   : Baseline TSH (mIU/L)
KD_T4     :  0.10  : Free T4 turnover (1/day)
KD_TSH    :  9.0   : TSH turnover (1/day)
AH_HCG    :  1.8e-5: TSH-equivalents per IU/L of hCG
FT4_50    : 14.0   : Free T4 at half-maximal TSH suppression
HILL_TSH  :  5.0   : TSH feedback Hill coefficient

// ---- drug PK (CL in L/day, V in L, KA in 1/day; Ki/IC50 in nM) --------
OND_KA    :  6.0   : Ondansetron absorption
OND_CL    : 660.0  : Ondansetron clearance
OND_V     : 160.0  : Ondansetron volume
OND_F     :  0.60  : Ondansetron bioavailability
OND_FU    :  0.27  : Ondansetron unbound fraction
OND_KP    :  0.60  : Ondansetron brain-to-plasma (unbound)
OND_KI    :  0.50  : Ondansetron 5-HT3 Ki
DOX_KA    :  1.1   : Doxylamine absorption
DOX_CL    : 84.0   : Doxylamine clearance
DOX_V     : 180.0  : Doxylamine volume
DOX_F     :  0.90  : Doxylamine bioavailability
DOX_FU    :  0.20  : Doxylamine unbound fraction
DOX_KP    :  0.80  : Doxylamine brain-to-plasma
DOX_KI_H1 :  5.0   : Doxylamine H1 Ki
MCP_KA    :  3.0   : Metoclopramide absorption
MCP_CL    : 720.0  : Metoclopramide clearance
MCP_V     : 200.0  : Metoclopramide volume
MCP_F     :  0.75  : Metoclopramide bioavailability
MCP_FU    :  0.70  : Metoclopramide unbound fraction
MCP_KP    :  0.25  : Metoclopramide brain-to-plasma
MCP_KI_D2 : 200.0  : Metoclopramide D2 Ki
PMZ_CL    : 1440   : Promethazine clearance
PMZ_V     : 970.0  : Promethazine volume
PMZ_FU    :  0.07  : Promethazine unbound fraction
PMZ_KP    :  2.0   : Promethazine brain-to-plasma
PMZ_KI_H1 :  2.0   : Promethazine H1 Ki
PMZ_KI_M1 : 60.0   : Promethazine M1 Ki
MIR_KA    :  1.4   : Mirtazapine absorption
MIR_CL    : 300.0  : Mirtazapine clearance
MIR_V     : 530.0  : Mirtazapine volume (t1/2 ~29 h -- generates the day-4 gap)
MIR_F     :  0.50  : Mirtazapine bioavailability
MIR_FU    :  0.15  : Mirtazapine unbound fraction
MIR_KP    :  1.00  : Mirtazapine brain-to-plasma
MIR_KI_H1 :  0.14  : Mirtazapine H1 Ki
MIR_KI_5T2: 50.0   : Mirtazapine 5-HT2A/2C Ki
MIR_KI_5T3:  8.1   : Mirtazapine 5-HT3 Ki
STE_KA    :  3.0   : Prednisolone absorption
STE_CL    : 192.0  : Prednisolone clearance
STE_V     : 45.0   : Prednisolone volume
STE_F     :  0.80  : Prednisolone bioavailability
STE_APP   :  0.10  : Prednisolone appetite gain (its ONLY action here)
STE_EC50  : 400.0  : Prednisolone EC50 for that appetite gain
GBP_KA    :  1.0   : Gabapentin absorption
GBP_CL    : 178.0  : Gabapentin clearance
GBP_V     : 60.0   : Gabapentin volume
GBP_F     :  0.45  : Gabapentin bioavailability
GBP_FU    :  0.97  : Gabapentin unbound fraction
GBP_KP    :  0.10  : Gabapentin CSF-to-plasma
GBP_IC50  : 3000   : Gabapentin alpha-2-delta IC50
CLO_KA    :  0.25  : Clonidine transdermal input rate
CLO_CL    : 302.0  : Clonidine clearance
CLO_V     : 190.0  : Clonidine volume
CLO_F     :  1.00  : Clonidine transdermal availability
CLO_FU    :  0.80  : Clonidine unbound fraction
CLO_KP    :  1.50  : Clonidine brain-to-plasma
CLO_KI_A2 :  3.0   : Clonidine alpha-2A Ki
MET_CL    : 1150   : Metformin clearance
MET_V     : 650.0  : Metformin volume
MET_EMAX  :  2.20  : Metformin fractional increase in basal GDF15 synthesis
MET_EC50  : 8000   : Metformin EC50 for that increase
MAB_CL    :  0.072 : Anti-GDF15 mAb clearance (t1/2 ~31 d)
MAB_V     :  3.2   : Anti-GDF15 mAb volume
MAB_KON   :  0.42  : Anti-GDF15 mAb ligand capture (1/day per nM)
RGD_CL    :  0.55  : Recombinant long-acting GDF15 clearance (t1/2 ~6 d)
RGD_V     :  5.0   : Recombinant GDF15 volume
RGD_POT   :  2.9e4 : Plasma GDF15 (pg/mL) contributed per nM of rGDF15

// ---- non-drug modifiers of basal GDF15 --------------------------------
THAL_ON   :  0     : Beta-thalassaemia flag (0/1)
THAL_FOLD : 12.0   : Beta-thalassaemia basal GDF15 fold-increase
TOBAC_ON  :  0     : Pre-conception tobacco flag (0/1)
TOBACCO_F :  0.30  : Tobacco fractional increase in basal GDF15
LOWGDF_ON :  0     : Low-GDF15 risk-allele flag (0/1)
LOWGDF_F  :  0.60  : Low-GDF15 allele multiplier on basal GDF15

// ---- supportive-care infusions (set by the scenario, not by dosing) ---
IV_FLUID  :  0     : Crystalloid infusion (L/day)
IV_NA_C   : 154.0  : Sodium in the crystalloid (mmol/L)
IV_CL_C   : 154.0  : Chloride in the crystalloid (mmol/L)
IV_K_RATE :  0     : Potassium infusion (mmol/day)
IV_THI    :  0     : Thiamine infusion (mg/day)
IV_DEX_G  :  0     : Dextrose infusion (g/day)

$CMT @annotated
PLAC   : Relative trophoblast mass
ISRS   : Trophoblast integrated stress response index
GDF    : Maternal plasma GDF15 (pg/mL)
SP     : Adapted log set-point of the GFRAL axis
GFRAL  : Relative GFRAL surface receptor pool
AP     : Area postrema activation
VAG    : Vagal afferent firing index
GAS    : Gastric emptying delay index
NTSC   : NTS integrator activity
CPG    : Emetic central pattern generator output
NAUS   : Nausea intensity (hours per 24 h)
HCG    : Maternal plasma hCG (IU/L)
VOL    : Extracellular fluid volume (L)
NATOT  : Exchangeable sodium (mmol)
CLTOT  : Exchangeable chloride (mmol)
HCO3   : Plasma bicarbonate (mmol/L)
KTOT   : Total body potassium (mmol)
WT     : Maternal weight (kg)
KET    : Beta-hydroxybutyrate (mmol/L)
THI    : Whole-body thiamine store (mg)
WERISK : Cumulative Wernicke hazard
ALT    : Alanine aminotransferase (U/L)
FT4    : Free T4 (pmol/L)
TSH    : TSH (mIU/L)
OND_D  : Ondansetron depot (nmol)
OND_C  : Ondansetron central (nmol)
DOX_D  : Doxylamine depot (nmol)
DOX_C  : Doxylamine central (nmol)
MCP_D  : Metoclopramide depot (nmol)
MCP_C  : Metoclopramide central (nmol)
PMZ_C  : Promethazine central (nmol)
MIR_D  : Mirtazapine depot (nmol)
MIR_C  : Mirtazapine central (nmol)
STE_D  : Prednisolone depot (nmol)
STE_C  : Prednisolone central (nmol)
GBP_D  : Gabapentin depot (nmol)
GBP_C  : Gabapentin central (nmol)
CLO_D  : Clonidine patch depot (nmol)
CLO_C  : Clonidine central (nmol)
MET_C  : Metformin central (nmol)
MAB_C  : Anti-GDF15 mAb central (nmol)
RGD_C  : Recombinant GDF15 central (nmol)

$GLOBAL
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define BASFOLD ( (THAL_ON > 0.5 ? THAL_FOLD : 1.0)                       \\
                * (TOBAC_ON > 0.5 ? (1.0 + TOBACCO_F) : 1.0)              \\
                * (LOWGDF_ON > 0.5 ? LOWGDF_F : 1.0) )
#define GDF0 (GDF_BASE * BASFOLD)

// smooth (cosine) transition, used for the perfusion and cytotrophoblast ramps
double hgramp(double t, double t0, double t1, double y0, double y1) {
  if (t <= t0) return y0;
  if (t >= t1) return y1;
  double f = 0.5 * (1.0 - cos(M_PI * (t - t0) / (t1 - t0)));
  return y0 + (y1 - y0) * f;
}
double hgocc(double c, double ki) {          // fractional receptor occupancy
  return (c > 0.0) ? c / (c + ki) : 0.0;
}
double hghill(double x, double ec50, double n) {
  if (x <= 0.0) return 0.0;
  double xn = pow(x, n);
  return xn / (pow(ec50, n) + xn);
}
// Continuous form of a PUQE-24 item.  The instrument scores 1-5 in bins; a
// step function would make the output jump, so we interpolate between edges.
double hgpuqe_item(double x, double b1, double b2, double b3, double b4) {
  if (x <= b1) return 1.0 + (x - 0.0) / (b1 > 0 ? b1 : 1e-9);
  if (x <= b2) return 2.0 + (x - b1) / (b2 - b1);
  if (x <= b3) return 3.0 + (x - b2) / (b3 - b2);
  if (x <= b4) return 4.0 + (x - b3) / (b4 - b3);
  return 5.0;
}

$MAIN
PLAC_0  = 0.0;
ISRS_0  = ISR_HI;
GDF_0   = GDF0;
SP_0    = log(GDF0);
GFRAL_0 = KSYN_GF / (KDEG_GF + KINT_GF * GDF0 / (GDF0 + KD_GF));
AP_0    = 0.0;
VAG_0   = 0.0;
GAS_0   = 0.0;
NTSC_0  = 0.0;
CPG_0   = 0.0;
NAUS_0  = 0.0;
HCG_0   = 0.0;
VOL_0   = VOL0;
NATOT_0 = NA_TOT0;
CLTOT_0 = CL_TOT0;
HCO3_0  = HCO30;
KTOT_0  = K_TOT0;
WT_0    = WT0;
KET_0   = 0.1;
THI_0   = THI0;
WERISK_0= 0.0;
ALT_0   = ALT_BASE;
FT4_0   = FT40;
TSH_0   = TSH0;

$ODE
// ===================== 1. placenta ======================================
double plac = PLAC > 0.0 ? PLAC : 0.0;
double cap  = 1.0 + PLAC_CAPG * (SOLVERTIME - T_IMP);
// implantation seeds the compartment over a 1-day window; without the seed the
// logistic term can never leave zero
dxdt_PLAC = (SOLVERTIME >= T_IMP)
          ? KG_PLAC * plac * (1.0 - plac / cap)
            + ((SOLVERTIME < T_IMP + 1.0) ? PLAC0 : 0.0)
          : 0.0;

// GDF15 is an ISR/ATF4-CHOP target gene.  Placental pO2 stays low until the
// intervillous circulation opens at ~10-12 weeks, so the SPECIFIC secretion
// rate is front-loaded -- highest while placental mass is still small.  That is
// what makes the GDF15 ramp so steep at GA 5-9 weeks and why it flattens early.
double isr_target = hgramp(SOLVERTIME, T_PERF1, T_PERF2, ISR_HI, ISR_LO);
dxdt_ISRS = (isr_target - ISRS) / TAU_ISR;

// ===================== 2. drug concentrations ===========================
double c_ond = (OND_C > 0 ? OND_C : 0) / OND_V;
double c_dox = (DOX_C > 0 ? DOX_C : 0) / DOX_V;
double c_mcp = (MCP_C > 0 ? MCP_C : 0) / MCP_V;
double c_pmz = (PMZ_C > 0 ? PMZ_C : 0) / PMZ_V;
double c_mir = (MIR_C > 0 ? MIR_C : 0) / MIR_V;
double c_ste = (STE_C > 0 ? STE_C : 0) / STE_V;
double c_gbp = (GBP_C > 0 ? GBP_C : 0) / GBP_V;
double c_clo = (CLO_C > 0 ? CLO_C : 0) / CLO_V;
double c_met = (MET_C > 0 ? MET_C : 0) / MET_V;
double c_mab = (MAB_C > 0 ? MAB_C : 0) / MAB_V;
double c_rgd = (RGD_C > 0 ? RGD_C : 0) / RGD_V;

double f_ond_p = c_ond * OND_FU;                 // peripheral (vagal terminals)
double f_ond_c = f_ond_p * OND_KP;               // NTS
double f_dox_c = c_dox * DOX_FU * DOX_KP;
double f_mcp_p = c_mcp * MCP_FU;
double f_mcp_c = f_mcp_p * MCP_KP;
double f_pmz_c = c_pmz * PMZ_FU * PMZ_KP;
double f_mir_c = c_mir * MIR_FU * MIR_KP;
double f_gbp_c = c_gbp * GBP_FU * GBP_KP;
double f_clo_c = c_clo * CLO_FU * CLO_KP;

// ===================== 3. maternal GDF15 ================================
double gdf = GDF > 1.0 ? GDF : 1.0;
// metformin raises circulating GDF15 (Coll 2020 PMID 31875646) -- this is its
// ONLY action in the model, and it is why timing decides its sign
double basfold = BASFOLD * (1.0 + MET_EMAX * hgocc(c_met, MET_EC50));
double syn_basal = GDF_BASE * KEL_GDF * basfold;
double syn_plac  = (SOLVERTIME >= T_IMP)
                 ? KS_GDF_PL * TROPH_GAIN * plac * ISRS : 0.0;
dxdt_GDF = syn_basal + syn_plac + RGD_POT * c_rgd
         - KEL_GDF * gdf - MAB_KON * c_mab * gdf;

// ===================== 4. THE FOLD-CHANGE DETECTOR ======================
// SP is the adapted log set-point of the GFRAL axis.  ALPHA mixes a pure
// fold-change detector (1) with a pure level detector (0).
dxdt_SP = (log(gdf) - SP) / TAU_SP;
double refc = exp(ALPHA * SP + (1.0 - ALPHA) * log(GREF));
double fold = gdf / refc;

dxdt_GFRAL = KSYN_GF - KDEG_GF * GFRAL
           - KINT_GF * GFRAL * gdf / (gdf + KD_GF);
double gdf_eff = fold * GFRAL / (KSYN_GF / KDEG_GF);

// ===================== 5. hCG ===========================================
double ctb = hgramp(SOLVERTIME, T_CTB1, T_CTB2, 1.0, CTB_PLAT);
double hcg = HCG > 0.0 ? HCG : 0.0;
dxdt_HCG = ((SOLVERTIME >= T_IMP) ? KS_HCG * TROPH_GAIN * plac * ctb : 0.0)
         - KEL_HCG * hcg;

// ===================== 6. area postrema =================================
double ap_target = SENS * hghill(gdf_eff, EC50_AP, HILL_AP)
                 * (1.0 - E0 * R_D2 * hgocc(f_mcp_c, MCP_KI_D2));
dxdt_AP = (ap_target - AP) / TAU_AP;

// ===================== 7. gastric emptying, peripheral 5-HT3 limb =======
double p4_idx = (SOLVERTIME - T_IMP) / 50.0;
if (p4_idx < 0.0) p4_idx = 0.0;
if (p4_idx > 1.0) p4_idx = 1.0;
double gas_target = (GAS_P4 * p4_idx + GAS_GDF * hghill(gdf_eff, EC50_AP, 1.4))
                  * (1.0 - 0.45 * hgocc(f_mcp_p, MCP_KI_D2));
dxdt_GAS = (gas_target - GAS) / TAU_GAS;
double gas = GAS > 0.0 ? GAS : 0.0;

double ec = hghill(gas, EC50_EC, 2.0);
double vag_target = ec * (1.0 - hgocc(f_ond_p, OND_KI));
dxdt_VAG = (vag_target - VAG) / TAU_AP;

double vest = (0.55 + 0.45 * hghill(gdf_eff, EC50_AP, 1.2))
            * (1.0 - hgocc(f_dox_c, DOX_KI_H1))
            * (1.0 - hgocc(f_pmz_c, PMZ_KI_H1))
            * (1.0 - hgocc(f_mir_c, MIR_KI_H1));
double hcg_limb = hcg / 60000.0;

// ===================== 8. NTS integrator ================================
// w_AP is the REMAINDER, so changing W_VAG re-apportions drive between the
// central and peripheral limbs without altering total drive.  That is what
// makes W_VAG identifiable from ondansetron alone while leaving the untreated
// natural history untouched.
double w_ap = W_TOT - W_VAG - W_VEST - W_HCG;
double nts_in = w_ap * AP + W_VAG * (VAG > 0 ? VAG : 0)
              + W_VEST * vest + W_HCG * hcg_limb;

double o_h1  = 1.0 - (1.0 - hgocc(f_dox_c, DOX_KI_H1))
                   * (1.0 - hgocc(f_pmz_c, PMZ_KI_H1))
                   * (1.0 - hgocc(f_mir_c, MIR_KI_H1));
double o_5t2 = hgocc(f_mir_c, MIR_KI_5T2);
double o_5t3 = 1.0 - (1.0 - hgocc(f_ond_c, OND_KI))
                   * (1.0 - hgocc(f_mir_c, MIR_KI_5T3));
double o_a2  = hgocc(f_clo_c, CLO_KI_A2);    // clonidine: presynaptic AGONIST.
                                             // Mirtazapine's alpha-2 ANTAGONISM
                                             // is deliberately NOT credited
                                             // here -- opposite sign, so its
                                             // antiemetic action is assigned to
                                             // H1/5-HT2/5-HT3 instead.
double o_a2d = hgocc(f_gbp_c, GBP_IC50);
double o_m1  = hgocc(f_pmz_c, PMZ_KI_M1);

double inh = (1.0 - E0 * R_H1    * o_h1)
           * (1.0 - E0 * R_5HT2  * o_5t2)
           * (1.0 - E0 * R_5HT3C * o_5t3)
           * (1.0 - E0 * R_A2    * o_a2)
           * (1.0 - E0 * R_A2D   * o_a2d)
           * (1.0 - E0 * R_M1    * o_m1);

dxdt_NTSC = (nts_in * inh - NTSC) / TAU_NTS;
double nts = NTSC > 0.0 ? NTSC : 0.0;

// ===================== 9. symptoms ======================================
dxdt_NAUS = (24.0 * hghill(nts, EC50_NAUS, HILL_NAUS) - NAUS) / TAU_NAUS;
dxdt_CPG  = (hghill(nts, EC50_VOM, HILL_VOM) - CPG) / TAU_CPG;
double vom    = VOM_MAX * (CPG > 0 ? CPG : 0);
double naus_h = NAUS > 0 ? NAUS : 0;

// ===================== 10. oral intake ==================================
double app_boost = STE_APP * hgocc(c_ste, STE_EC50)
                 + 0.12 * hgocc(f_mir_c, MIR_KI_H1);
// intake collapses steeply, not linearly: a woman nauseated 16 h a day is not
// eating two thirds of normal
double oral_base = 1.0 - 0.985 * naus_h / 24.0;
if (oral_base < 0.0) oral_base = 0.0;
double oral = pow(oral_base, 1.5);
if (oral < 0.04) oral = 0.04;
oral = oral * (1.0 + app_boost);
if (oral > 1.0) oral = 1.0;

// ===================== 11. fluid and electrolytes =======================
double vol  = VOL > 5.0 ? VOL : 5.0;
double aldo = (VOL0 - vol) / VOL0 * 6.0;
if (aldo < 0.0) aldo = 0.0;
if (aldo > 1.0) aldo = 1.0;
double l_vom = VOM_VOL * vom;

// A nauseated woman stops EATING long before she stops SIPPING, and thirst
// rises with depletion.  That is why HG produces hyponatraemia rather than the
// hypernatraemia that loss of (hypotonic) gastric juice alone would predict.
double thirst   = 1.0 + THIRST_G * aldo;
double water_in = ORAL_FLUID * (WATER_FLR + (1.0 - WATER_FLR) * oral) * thirst;
double urine    = water_in + IV_FLUID - INSENS + KREN_VOL * (vol - VOL0);
if (urine < URINE_MIN) urine = URINE_MIN;
dxdt_VOL = water_in + IV_FLUID - urine - INSENS - l_vom;

double cl_conc  = (CLTOT > 1.0 ? CLTOT : 1.0) / vol;
double na_conc  = (NATOT > 1.0 ? NATOT : 1.0) / vol;
double cl_avail = (cl_conc - 84.0) / (103.0 - 84.0);
if (cl_avail < 0.0) cl_avail = 0.0;
if (cl_avail > 1.0) cl_avail = 1.0;

// Gastric HCl secretion needs plasma chloride to draw on, so severe chloride
// depletion is self-limiting: the vomitus becomes less acid and less
// chloride-rich as the deficit grows.
double clscale = cl_conc / 103.0;
if (clscale < 0.35) clscale = 0.35;
if (clscale > 1.20) clscale = 1.20;

double cl_in = CL_INTAKE * oral + IV_CL_C * IV_FLUID;
double cl_gi = GJ_CL * clscale * l_vom;
double cl_net = cl_in - cl_gi;
dxdt_CLTOT = cl_net - (cl_net > 0 ? cl_net : 0) * (1.0 - 0.85 * aldo);

// Classic chloride-responsive metabolic alkalosis: the kidney cannot excrete
// the excess bicarbonate without chloride to accompany it, so the alkalosis
// persists until Cl- is replaced -- however much of anything else is given.
double hco3_excr = KREN_HCO3 * ((HCO3 - HCO30) > 0 ? (HCO3 - HCO30) : 0.0)
                 * (HCO3_ESC + (1.0 - HCO3_ESC) * cl_avail);
dxdt_HCO3 = GJ_H * clscale * l_vom / vol - hco3_excr;

// Much of the excreted HCO3- is obligated to carry Na+ with it.  This is why a
// volume-depleted woman whose kidney is retaining sodium as hard as it can
// still becomes hyponatraemic.  As Na+ falls the coupling gives way to
// K+/NH4+ excretion, which puts a floor under it.
double f_na = (na_conc - NA_FLOOR) / (140.0 - NA_FLOOR);
if (f_na < 0.0) f_na = 0.0;
if (f_na > 1.0) f_na = 1.0;
double na_in = NA_INTAKE * oral + IV_NA_C * IV_FLUID;
double na_net = na_in - GJ_NA * l_vom;
dxdt_NATOT = na_net - (na_net > 0 ? na_net : 0) * (1.0 - 0.85 * aldo)
           - hco3_excr * vol * NA_HCO3_C * f_na;

// renal K+ handling must be able to RETAIN as well as waste, or a depleted
// store never refills once vomiting stops
double ktot  = KTOT > 1.0 ? KTOT : 1.0;
double k_ren = (ktot / K_TOT0 - K_RETAIN) / (1.0 - K_RETAIN);
if (k_ren < 0.05) k_ren = 0.05;
if (k_ren > 1.00) k_ren = 1.00;
dxdt_KTOT = K_INTAKE * oral + IV_K_RATE - GJ_K * l_vom
          - K_URINE0 * (1.0 + K_ALDO * aldo) * (0.15 + 0.85 * oral) * k_ren;

// ===================== 12. energy, weight, ketones, liver ===============
double kcal = KCAL_FULL * oral + IV_DEX_G * 3.4;
// adaptive thermogenesis: resting expenditure falls with a sustained deficit
double kcal_out = KCAL_NEED * (0.80 + 0.20 * oral);
double fetal = (SOLVERTIME - 84.0) / 40.0;
if (fetal < 0.0) fetal = 0.0;
if (fetal > 1.0) fetal = 1.0;
dxdt_WT = (kcal - kcal_out) / KCAL_PERKG + dxdt_VOL + WT_FETAL * fetal;

double ket_ss = 0.1 + KET_MAX * ((1.0 - kcal / 1250.0) > 0
                                 ? (1.0 - kcal / 1250.0) : 0.0);
dxdt_KET = (ket_ss - KET) / TAU_KET;

double starv = (1.0 - kcal / 1600.0) > 0 ? (1.0 - kcal / 1600.0) : 0.0;
dxdt_ALT = (ALT_BASE + (ALT_MAX - ALT_BASE) * pow(starv, 1.6) - ALT) / TAU_ALT;

// ===================== 13. thiamine and Wernicke hazard =================
double thi = THI > 0.0 ? THI : 0.0;
double refeed = 1.0 + THI_REFEED * ((1.0 - thi / THI_CRIT) > 0
                                    ? (1.0 - thi / THI_CRIT) : 0.0);
dxdt_THI = THI_DIET * oral + IV_THI - KEL_THI * thi
         - THI_PERCHO * IV_DEX_G * refeed;
double thidef = (1.0 - thi / THI_CRIT) > 0 ? (1.0 - thi / THI_CRIT) : 0.0;
dxdt_WERISK = HZ_WE * thidef * thidef;

// ===================== 14. thyroid ======================================
// hCG cross-activates the TSH receptor (Rodien 2004 PMID 15073140).  Because
// hCG and GDF15 are both syncytiotrophoblast products, biochemical
// thyrotoxicosis tracks TROPH_GAIN and NOT SENS.
double tshr = (TSH > 0 ? TSH : 0) + AH_HCG * hcg;
dxdt_FT4 = KD_T4 * FT40 * (tshr / TSH0) - KD_T4 * FT4;
dxdt_TSH = KD_TSH * TSH0 * 2.0 / (1.0 + pow(FT4 / FT4_50, HILL_TSH))
         - KD_TSH * TSH;

// ===================== 15. drug PK ======================================
dxdt_OND_D = -OND_KA * OND_D;
dxdt_OND_C =  OND_F * OND_KA * OND_D - OND_CL * (OND_C > 0 ? OND_C : 0) / OND_V;
dxdt_DOX_D = -DOX_KA * DOX_D;
dxdt_DOX_C =  DOX_F * DOX_KA * DOX_D - DOX_CL * (DOX_C > 0 ? DOX_C : 0) / DOX_V;
dxdt_MCP_D = -MCP_KA * MCP_D;
dxdt_MCP_C =  MCP_F * MCP_KA * MCP_D - MCP_CL * (MCP_C > 0 ? MCP_C : 0) / MCP_V;
dxdt_PMZ_C = -PMZ_CL * (PMZ_C > 0 ? PMZ_C : 0) / PMZ_V;
dxdt_MIR_D = -MIR_KA * MIR_D;
dxdt_MIR_C =  MIR_F * MIR_KA * MIR_D - MIR_CL * (MIR_C > 0 ? MIR_C : 0) / MIR_V;
dxdt_STE_D = -STE_KA * STE_D;
dxdt_STE_C =  STE_F * STE_KA * STE_D - STE_CL * (STE_C > 0 ? STE_C : 0) / STE_V;
dxdt_GBP_D = -GBP_KA * GBP_D;
dxdt_GBP_C =  GBP_F * GBP_KA * GBP_D - GBP_CL * (GBP_C > 0 ? GBP_C : 0) / GBP_V;
dxdt_CLO_D = -CLO_KA * CLO_D;
dxdt_CLO_C =  CLO_F * CLO_KA * CLO_D - CLO_CL * (CLO_C > 0 ? CLO_C : 0) / CLO_V;
dxdt_MET_C = -MET_CL * (MET_C > 0 ? MET_C : 0) / MET_V;
dxdt_MAB_C = -MAB_CL * (MAB_C > 0 ? MAB_C : 0) / MAB_V;
dxdt_RGD_C = -RGD_CL * (RGD_C > 0 ? RGD_C : 0) / RGD_V;

$TABLE
double GDF15    = GDF;
double SETPOINT = exp(SP);
double FOLD     = GDF / exp(ALPHA * SP + (1.0 - ALPHA) * log(GREF));
double DRIVE    = FOLD * GFRAL / (KSYN_GF / KDEG_GF);
double NAUSEA_H = NAUS > 0 ? NAUS : 0;
double VOMITS   = VOM_MAX * (CPG > 0 ? CPG : 0);
double RETCHES  = VOMITS * RETCH_RAT;
double PUQE     = hgpuqe_item(NAUSEA_H, 1.0, 2.0, 4.0, 6.0)
                + hgpuqe_item(VOMITS,   1.0, 3.0, 5.0, 7.0)
                + hgpuqe_item(RETCHES,  1.0, 3.0, 5.0, 7.0);
if (PUQE > 15.0) PUQE = 15.0;
double K_PLASMA  = 4.0 + (KTOT - K_TOT0) / 300.0;
double NA_PLASMA = NATOT / (VOL > 1.0 ? VOL : 1.0);
double CL_PLASMA = CLTOT / (VOL > 1.0 ? VOL : 1.0);
double HCO3_PL   = HCO3;
double THIAMINE  = THI;
double P_WE      = 1.0 - exp(-WERISK);
double WT_LOSS_P = (WT0 - WT) / WT0 * 100.0;
double GA_WEEKS  = SOLVERTIME / 7.0;
double C_OND     = OND_C / OND_V;
double C_MIR     = MIR_C / MIR_V;
double C_GBP     = GBP_C / GBP_V;
double C_CLO     = CLO_C / CLO_V;
double C_MET     = MET_C / MET_V;

$CAPTURE GA_WEEKS GDF15 SETPOINT FOLD DRIVE NAUSEA_H VOMITS RETCHES PUQE
$CAPTURE K_PLASMA NA_PLASMA CL_PLASMA HCO3_PL THIAMINE P_WE WT_LOSS_P
$CAPTURE C_OND C_MIR C_GBP C_CLO C_MET
'

mod <- mcode("hg_qsp", hg_code, atol = 1e-8, rtol = 1e-8)

# =============================================================================
#  DOSING HELPERS
#  Molecular weights, so every dose can be written in mg and converted to nmol.
# =============================================================================
MW <- c(OND = 293.4, DOX = 270.4, MCP = 299.8, PMZ = 284.4, MIR = 265.4,
        STE = 360.4, GBP = 171.2, CLO = 230.1, MET = 129.2)
nmol <- function(mg, drug) mg * 1e6 / MW[[drug]]

CMTN <- function(name) which(names(mrgsolve::init(mod)) == name)

dose <- function(cmt, amt, time, ii = 0, addl = 0, rate = 0) {
  data.frame(ID = 1, time = time, cmt = CMTN(cmt), amt = amt,
             evid = 1, ii = ii, addl = addl, rate = rate)
}

# The model starts pre-conception (t < 0), so every scenario runs from -70 days.
T_START <- -70
T_END   <- 200

# =============================================================================
#  SCENARIOS
#  Sixteen, grouped so the model's three commitments are each exercised.
# =============================================================================

# ---- A. natural history / cohorts -------------------------------------------
sc_normal <- list(
  name = "1. Normal pregnancy",
  par  = list(SENS = 0.62), ev = NULL)

sc_hg_sens <- list(
  name = "2. HG, maternal-sensitivity phenotype",
  par  = list(SENS = 1.00), ev = NULL)

sc_hg_prod <- list(
  name = "3. HG, fetal-production phenotype (hCG rises too)",
  par  = list(SENS = 0.62, TROPH_GAIN = 2.00), ev = NULL)

sc_lowgdf <- list(
  name = "4. Low-GDF15 risk allele (low pre-pregnancy set-point)",
  par  = list(SENS = 1.00, LOWGDF_ON = 1), ev = NULL)

sc_thal <- list(
  name = "5. Beta-thalassaemia (lifelong high GDF15) -- predicted protected",
  par  = list(SENS = 1.00, THAL_ON = 1), ev = NULL)

sc_level_model <- list(
  name = "6. COMPARATOR: level detector (ALPHA = 0), same equations",
  par  = list(SENS = 1.00, ALPHA = 0.0), ev = NULL)

# ---- B. the VOMIT trial (Ostenfeld 2026) -----------------------------------
TR <- 56   # randomisation at GA 8 weeks
sc_placebo <- list(
  name = "7. VOMIT: placebo", par = list(SENS = 1.00), ev = NULL)

sc_ondansetron <- list(
  name = "8. VOMIT: ondansetron 8 mg PO q8h x 14 d (observed dPUQE -0.51)",
  par  = list(SENS = 1.00),
  ev   = dose("OND_D", nmol(8, "OND"), TR, ii = 1/3, addl = 41))

sc_mirtazapine <- list(
  name = "9. VOMIT: mirtazapine 30 mg PO qHS x 14 d (observed dPUQE -1.86)",
  par  = list(SENS = 1.00),
  ev   = dose("MIR_D", nmol(30, "MIR"), TR, ii = 1, addl = 13))

# ---- C. the node-position law ----------------------------------------------
sc_doxylamine <- list(
  name = "10. Doxylamine/pyridoxine 10-10-20 mg daily (H1 + vestibular)",
  par  = list(SENS = 1.00),
  ev   = rbind(dose("DOX_D", nmol(10, "DOX"), TR + 0.30, ii = 1, addl = 13),
               dose("DOX_D", nmol(10, "DOX"), TR + 0.65, ii = 1, addl = 13),
               dose("DOX_D", nmol(20, "DOX"), TR + 0.92, ii = 1, addl = 13)))

sc_gabapentin <- list(
  name = "11. Gabapentin 600 mg q8h (NTS alpha-2-delta -- highest authority)",
  par  = list(SENS = 1.00),
  ev   = dose("GBP_D", nmol(600, "GBP"), TR, ii = 1/3, addl = 41))

sc_clonidine <- list(
  name = "12. Clonidine 5 mg patch, ~0.15 mg/day (NTS presynaptic alpha-2)",
  par  = list(SENS = 1.00),
  ev   = dose("CLO_D", nmol(0.15, "CLO"), TR, ii = 1, addl = 13))

sc_steroid <- list(
  name = "13. Methylprednisolone 125 mg IV + prednisone taper (predicted null)",
  par  = list(SENS = 1.00),
  ev   = rbind(dose("STE_C", nmol(125, "STE"), TR),
               dose("STE_D", nmol(40,  "STE"), TR + 1),
               dose("STE_D", nmol(20,  "STE"), TR + 2, ii = 1, addl = 2),
               dose("STE_D", nmol(10,  "STE"), TR + 5, ii = 1, addl = 2),
               dose("STE_D", nmol(5,   "STE"), TR + 8, ii = 1, addl = 6)))

# ---- D. prevention: the window closes at conception ------------------------
sc_metformin_pre <- list(
  name = "14. Metformin 1 g BID PRE-CONCEPTION only, stopped at positive test",
  par  = list(SENS = 1.00),
  ev   = dose("MET_C", nmol(1000, "MET") * 0.55, T_START, ii = 0.5, addl = 167))

sc_metformin_late <- list(
  name = "15. The SAME metformin started at GA 6 wk -- predicted no benefit",
  par  = list(SENS = 1.00),
  ev   = dose("MET_C", nmol(1000, "MET") * 0.55, 42, ii = 0.5, addl = 195))

sc_antigdf15 <- list(
  name = "16. Anti-GDF15 mAb 90 nmol at GA 8 wk (largest predicted effect; UNTESTED in pregnancy)",
  par  = list(SENS = 1.00),
  ev   = dose("MAB_C", 90, TR))

SCENARIOS <- list(sc_normal, sc_hg_sens, sc_hg_prod, sc_lowgdf, sc_thal,
                  sc_level_model, sc_placebo, sc_ondansetron, sc_mirtazapine,
                  sc_doxylamine, sc_gabapentin, sc_clonidine, sc_steroid,
                  sc_metformin_pre, sc_metformin_late, sc_antigdf15)

# =============================================================================
#  SUPPORTIVE CARE
#  Given as PARAMETER changes over a window rather than as dosing events,
#  because IV fluid/thiamine/dextrose are continuous rates. Splitting the run
#  at the window edges keeps them piecewise-constant.
# =============================================================================
run_with_support <- function(par = list(), ev = NULL,
                             t_on = 56, t_off = 84,
                             fluid = 3.0, k_rate = 40, thi = 100, dex = 150,
                             end = T_END, delta = 0.25) {
  base_par <- par
  p1 <- mod %>% param(base_par) %>%
    param(IV_FLUID = 0, IV_K_RATE = 0, IV_THI = 0, IV_DEX_G = 0)
  seg1 <- p1 %>% mrgsim(events = if (is.null(ev)) NULL else ev,
                        start = T_START, end = t_on, delta = delta)
  y1 <- unlist(tail(as.data.frame(seg1)[, names(init(mod))], 1))

  p2 <- mod %>% param(base_par) %>%
    param(IV_FLUID = fluid, IV_K_RATE = k_rate, IV_THI = thi, IV_DEX_G = dex) %>%
    init(as.list(y1))
  seg2 <- p2 %>% mrgsim(events = if (is.null(ev)) NULL else ev,
                        start = t_on, end = t_off, delta = delta)
  y2 <- unlist(tail(as.data.frame(seg2)[, names(init(mod))], 1))

  p3 <- mod %>% param(base_par) %>%
    param(IV_FLUID = 0, IV_K_RATE = 0, IV_THI = 0, IV_DEX_G = 0) %>%
    init(as.list(y2))
  seg3 <- p3 %>% mrgsim(events = if (is.null(ev)) NULL else ev,
                        start = t_off, end = end, delta = delta)

  bind_rows(as.data.frame(seg1), as.data.frame(seg2), as.data.frame(seg3)) %>%
    distinct(time, .keep_all = TRUE) %>% arrange(time)
}

# =============================================================================
#  RUN
# =============================================================================
run_scenario <- function(sc, end = T_END, delta = 0.25) {
  m <- mod
  if (length(sc$par)) m <- m %>% param(sc$par)
  out <- m %>% mrgsim(events = sc$ev, start = T_START, end = end, delta = delta)
  as.data.frame(out) %>% mutate(scenario = sc$name)
}

message("Running ", length(SCENARIOS), " HG scenarios ...")
sims <- bind_rows(lapply(SCENARIOS, run_scenario))

# ---- summary ---------------------------------------------------------------
summ <- sims %>%
  filter(time >= 28, time <= 160) %>%
  group_by(scenario) %>%
  summarise(
    PUQE_max     = max(PUQE),
    peak_GA_wk   = GA_WEEKS[which.max(PUQE)],
    PUQE_at_16wk = PUQE[which.min(abs(time - 112))],
    GDF15_max    = max(GDF15),
    fold_at_9wk  = FOLD[which.min(abs(time - 63))],
    wt_loss_pct  = max(WT_LOSS_P),
    K_min        = min(K_PLASMA),
    HCO3_max     = max(HCO3_PL),
    Cl_min       = min(CL_PLASMA),
    THI_min      = min(THIAMINE),
    P_WE_pct     = max(P_WE) * 100,
    TSH_min      = min(TSH),
    .groups      = "drop")
print(as.data.frame(summ), digits = 3)

# ---- the VOMIT comparison --------------------------------------------------
# NOTE: a single deterministic trajectory OVERSTATES a drug's mean trial effect,
# because PUQE-24 has a hard ceiling of 15 and the transductions are steep Hill
# functions.  The published effect sizes are reproduced by averaging over a
# heterogeneous cohort -- see hg_reference_impl.py, which sweeps SENS over
# c(0.82, 0.88, 0.94, 1.00, 1.08, 1.18, 1.30) and recovers -0.52 and -1.92 for
# ondansetron and mirtazapine.  Use trial_population() below to do the same here.
vomit_delta <- function(sims, arm, ref = sc_placebo$name, day) {
  g <- function(nm, t) {
    s <- sims[sims$scenario == nm, ]
    s$PUQE[which.min(abs(s$time - t))]
  }
  (g(arm, TR + day) - g(arm, TR)) - (g(ref, TR + day) - g(ref, TR))
}
cat("\nVOMIT trial, single typical subject (see caveat above):\n")
for (d in c(2, 4, 7, 14)) {
  cat(sprintf("  day %2d   ondansetron %+6.2f   mirtazapine %+6.2f\n", d,
              vomit_delta(sims, sc_ondansetron$name, day = d),
              vomit_delta(sims, sc_mirtazapine$name, day = d)))
}

trial_population <- function(arm_ev, sens = c(0.82, 0.88, 0.94, 1.00,
                                              1.08, 1.18, 1.30),
                             days = c(2, 4, 7, 14)) {
  res <- lapply(sens, function(s) {
    a <- run_scenario(list(name = "arm", par = list(SENS = s), ev = arm_ev),
                      end = TR + 22)
    p <- run_scenario(list(name = "pbo", par = list(SENS = s), ev = NULL),
                      end = TR + 22)
    gp <- function(df, t) df$PUQE[which.min(abs(df$time - t))]
    sapply(days, function(d)
      (gp(a, TR + d) - gp(a, TR)) - (gp(p, TR + d) - gp(p, TR)))
  })
  setNames(rowMeans(do.call(cbind, res)), paste0("day", days))
}
cat("\nVOMIT trial, 7-subject cohort mean (this is the published comparison):\n")
cat("  ondansetron: ",
    paste(sprintf("%+.2f", trial_population(sc_ondansetron$ev)), collapse = "  "),
    "\n")
cat("  mirtazapine: ",
    paste(sprintf("%+.2f", trial_population(sc_mirtazapine$ev)), collapse = "  "),
    "\n")
cat("  observed:     -0.51 (ondansetron) and -1.86 (mirtazapine) at day 2\n")

# ---- the thiamine / dextrose ordering result -------------------------------
cat("\nSupportive care, severe protracted course (SENS = 1.18):\n")
support_arms <- list(
  "IV saline + KCl, no dextrose" = list(dex = 0,   thi = 0),
  "IV dextrose, NO thiamine"     = list(dex = 150, thi = 0),
  "IV dextrose + thiamine"       = list(dex = 150, thi = 100))
for (nm in names(support_arms)) {
  a <- support_arms[[nm]]
  df <- run_with_support(par = list(SENS = 1.18), dex = a$dex, thi = a$thi,
                         end = 150)
  w <- df[df$time >= 42 & df$time <= 145, ]
  cat(sprintf("  %-30s THI_min %5.1f mg   P(WE) %5.1f%%   HCO3_max %4.1f\n",
              nm, min(w$THIAMINE), max(w$P_WE) * 100, max(w$HCO3_PL)))
}
cat("  -> dextrose without thiamine is better hydrated and WORSE for Wernicke\n")
cat("     risk than saline alone. Thiamine first, then glucose.\n")

# =============================================================================
#  PLOTS
# =============================================================================
theme_set(theme_bw(base_size = 10))

# The single most important figure: the hormone keeps rising, the drive does not.
p_core <- sims %>%
  filter(scenario == sc_hg_sens$name, time >= 14, time <= 200) %>%
  select(GA_WEEKS, GDF15, SETPOINT) %>%
  tidyr::pivot_longer(-GA_WEEKS) %>%
  ggplot(aes(GA_WEEKS, value, colour = name)) +
  geom_line(linewidth = 0.9) +
  scale_y_log10() +
  labs(title = "GDF15 and its adapted set-point",
       subtitle = "the set-point chases the hormone with a ~4 week lag",
       x = "gestational age (weeks)", y = "pg/mL (log)", colour = NULL)

p_drive <- sims %>%
  filter(scenario %in% c(sc_hg_sens$name, sc_level_model$name,
                         sc_thal$name, sc_normal$name),
         time >= 14, time <= 200) %>%
  ggplot(aes(GA_WEEKS, PUQE, colour = scenario)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 13, linetype = 2, colour = "grey40") +
  labs(title = "PUQE-24: fold-change detector vs level detector",
       subtitle = "ALPHA = 0 never remits and makes beta-thalassaemia the worst case",
       x = "gestational age (weeks)", y = "PUQE-24 (3-15)", colour = NULL) +
  theme(legend.position = "bottom", legend.direction = "vertical")

p_elec <- sims %>%
  filter(scenario == sc_hg_sens$name, time >= 28, time <= 160) %>%
  select(GA_WEEKS, K_PLASMA, CL_PLASMA, HCO3_PL, NA_PLASMA) %>%
  tidyr::pivot_longer(-GA_WEEKS) %>%
  ggplot(aes(GA_WEEKS, value)) +
  geom_line(linewidth = 0.8, colour = "#00838f") +
  facet_wrap(~name, scales = "free_y") +
  labs(title = "Hypochloraemic metabolic alkalosis with hypokalaemia",
       x = "gestational age (weeks)", y = "mmol/L")

p_thi <- sims %>%
  filter(scenario == sc_hg_sens$name, time >= 28, time <= 160) %>%
  ggplot(aes(GA_WEEKS)) +
  geom_line(aes(y = THIAMINE), linewidth = 0.8, colour = "#c62828") +
  geom_hline(yintercept = 12, linetype = 2) +
  labs(title = "Thiamine store -- the slowest clock in the illness",
       subtitle = "PUQE can normalise while this is still falling",
       x = "gestational age (weeks)", y = "whole-body thiamine (mg)")

print(p_core); print(p_drive); print(p_elec); print(p_thi)

# =============================================================================
#  END.  Cross-check every quoted number with:
#      python3 hg_reference_impl.py --check
# =============================================================================
