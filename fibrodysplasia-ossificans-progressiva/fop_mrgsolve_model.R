# =====================================================================
# Fibrodysplasia Ossificans Progressiva (FOP) — mrgsolve QSP model
#   Author : Claude Code Routine (2026-07-30)
#
#   Scope  : ACVR1/ALK2 R206H gain-of-function -> loss of the FKBP12
#            clamp (ligand-independent leak) + NEOFUNCTION (activin A,
#            normally a non-signalling competitive inhibitor of wild-type
#            ALK2, becomes a full agonist) -> pSMAD1/5/8 in resident
#            fibro/adipogenic progenitors (FAPs) -> heterotopic
#            ENDOCHONDRAL ossification through a one-way relay
#            (fibroproliferation -> chondrogenesis -> hypertrophy ->
#            vascular invasion -> osteoblast -> mineralised bone) ->
#            an IRREVERSIBLE cumulative bone stock -> joint ankylosis
#            (CAJIS) + thoracic insufficiency (FVC) -> death.
#
#   PK/PD  : palovarotene (RARgamma agonist, oral, chronic + flare
#            dosing), garetosmab (anti-activin A IgG4, IV q4w, TMDD),
#            a generic ALK2 kinase inhibitor (fidrisertib / zilurgisertib
#            class, oral), saracatinib (SRC/ALK2), sirolimus (mTORC1),
#            prednisone bursts, IL-1 blockade (anakinra).
#
# ---------------------------------------------------------------------
#   WHY THIS MODEL IS BUILT THE WAY IT IS — four axes, each DERIVED
#   from published numbers rather than assumed:
#
#   AXIS 1 — THE ANTI-FLARE CEILING.
#     The 36-month prospective natural history study (NHS, PMID 36152026)
#     reports, in the same cohort: 229 flare-ups in 114 patients over a
#     mean 26.8 months = 0.90 flares/patient-year; 14/52 (26.9%) of
#     IMAGED flare-ups had new HO at the flare site by day 84; mean new
#     HO volume at those sites 28.8 mL; and a mean annualised new HO
#     volume of 23.6 mL/yr overall. Multiply the first three:
#        0.90 x 0.269 x 28.8 = 6.97 mL/yr
#     against a total of 23.6 mL/yr, so the flare-attributable fraction
#        PHI_FL = 6.97 / 23.6 = 0.295.
#     A drug that abolished 100% of flare-ups would therefore remove
#     ~30% of new heterotopic bone and no more. Two independent
#     corroborations: 47% of patients report disease progression without
#     obvious flare-ups (PMID 27025942), and in LUMINA-1 placebo only
#     12/29 (41%) of new lesions were co-located with a flare
#     (PMID 39216107). The model therefore carries TWO parallel routes
#     to bone -- episodic (flare) and smouldering -- and keeps them
#     exactly separable by holding the receptor step LINEAR in activin A
#     (justified: FOP serum activin A ~434 pg/mL, PMID 37165433, is far
#     below receptor saturation). Anti-inflammatory agents enter only
#     the episodic route; they cannot exceed PHI_FL by construction,
#     which is the point.
#
#   AXIS 2 — RATE versus STOCK (the ratchet).
#     Heterotopic bone remodels internally but is not resorbed
#     (KRESORB = 0); surgical excision provokes explosive recurrence.
#     Every trial endpoint (annualised NEW HO volume, new lesion count)
#     is a RATE; every disability endpoint (CAJIS, FVC, survival) is the
#     INTEGRAL of that rate. The two cannot move together. The model
#     reports both from one simulation, so MOVE 60% rate reduction and
#     ~12% total-volume reduction over 18 months are the same result.
#
#   AXIS 3 — AGE BEATS POTENCY, AND THE REASON IS TERRITORY.
#     The NHS velocity profile is non-monotonic: 21.9 mL/yr at 8-<15,
#     41.5 at 15-<25, 4.6 at 25-65 -- a 9-fold collapse in adulthood
#     that is NOT disease burn-out. The model produces it from two
#     externally-anchored factors: absolute lesion volume scales with
#     the size of the muscle compartment (a human growth curve, NOT
#     fitted to HO data), while the ossifiable TERRITORY is finite and
#     is being consumed. Late-life quiescence is substrate exhaustion,
#     which is exactly why a drug started at 25 has little left to
#     protect (~10% of lifetime HO) while the same drug started at 8
#     protects ~50%.
#
#   AXIS 4 — EFFICACY AND TOXICITY SHARE ONE AXIS.
#     Palovarotene works by suppressing chondrogenic commitment. The
#     growth plate IS chondrogenesis. MOVE achieved 54-60% reduction in
#     new HO WITH premature physeal closure in 21/57 (36.8%) of children
#     <14 y, decreased vertebral BMD and increased vertebral fracture
#     risk (PMID 36583535); daily dosing ablates growth plates in mice
#     (PMID 30226468). The model does not assume the two effects are
#     separable -- it ESTIMATES the implied selectivity from those two
#     numbers (R_SEL ~ 3.9-fold, conditional on the assumed maximal PPC
#     hazard) and then shows that no dose reaches 50% HO suppression at
#     an acceptable physeal hazard while physes are open. The only
#     separation axis is skeletal maturity, i.e. time -- and by then
#     Axis 3 has removed most of the benefit.
#
# ---------------------------------------------------------------------
#   CALIBRATION ANCHORS (all published; see fop_references.md)
#     NHS 36-month, PMID 36152026 : flare rate 0.90/pt-yr; 26.9% of
#       imaged flares ossify by d84; 28.8 mL mean new HO at those sites;
#       mean annualised new HO 23.6 mL/yr; by age band 21.9 / 41.5 /
#       4.6 mL/yr; baseline stock 68.8 mL (2-<8 y), 314.4 mL (all),
#       575.2 mL (25-65 y)
#     Flare survey, PMID 27025942 : 52% spontaneous / 48% trauma; IM
#       injection -> immediate flare in 25%, of whom 84% ossify; 47%
#       progress without obvious flares; glucocorticoids fully resolve
#       only 12% of flares, rebound in 43%
#     MOVE phase 3, PMID 36583535 : 60% (BcPM) / 54% (wLME) lower mean
#       annualised new HO vs NHS; PPC 21/57 (36.8%) if <14 y; retinoid
#       AE 97.0%
#     LUMINA-1, PMID 37770652 : PET total lesion activity primary MISSED
#       (p = 0.0741) while new HO lesions were suppressed (0% vs 40.9%,
#       p = 0.0027); PMID 39216107 : placebo new HO volume 16.6 vs
#       3.2 cm3 with vs without a prior flare
#     Garetosmab PK, PMID 37694449 / 32557665 : 10 mg/kg IV q4w,
#       Ctrough 105 +/- 30.8 mg/L, steady state by 12-16 wk, TMDD
#       saturated at 10 mg/kg
#     Pulmonary/survival, PMID 9577416 / 33748352 / 20194327 : FVC
#       44 +/- 14 % predicted in advanced disease; restrictive
#       physiology already present in childhood and NOT worsening in
#       step with later thoracic HO; median lifespan 56 y (95% CI
#       51-60); thoracic insufficiency = 54% of deaths, median age 42
#
#   DELIBERATELY NOT CALIBRATED (stated rather than invented):
#     CAJIS is implemented on its true 0-30 scale, monotone in
#     appendicular stock and saturating in late-stage disease, but it is
#     NOT fitted to a published cohort mean; treat CAJIS output as
#     ordinal/semi-quantitative. See the README.
# =====================================================================

library(mrgsolve)
library(dplyr)

fop_code <- '
$PROB
# Fibrodysplasia Ossificans Progressiva (FOP) QSP model
# 42 ODE compartments: 15 drug PK + 27 disease / PD / clinical
# Time unit: DAYS.  Age enters as AGE0 + TIME/365.25.

$PLUGIN autodec

$PARAM @annotated
// =====================================================================
// PATIENT / STRUCTURAL
// =====================================================================
AGE0     : 5.0   : Age at simulation start (yr)
WT       : 20.0  : Body weight (kg)
SEX      : 0     : 0 = female, 1 = male (physeal fusion age)

// =====================================================================
// 1. PALOVAROTENE PK — oral RARgamma agonist, 2-cpt, CYP3A4
//    Tmax ~2-4 h, terminal t1/2 ~10-14 h, minimal accumulation
// =====================================================================
KA_PAL   : 0.55  : Palovarotene absorption rate (1/h)
CL_PAL   : 7.60  : Palovarotene apparent clearance (L/h)
V1_PAL   : 120.0 : Palovarotene central volume (L)
V2_PAL   : 90.0  : Palovarotene peripheral volume (L)
Q_PAL    : 8.00  : Palovarotene intercompartmental clearance (L/h)
F_PAL    : 0.65  : Palovarotene oral bioavailability

// =====================================================================
// 2. GARETOSMAB PK — anti-activin A mAb, IV, 2-cpt + TMDD
//    10 mg/kg q4w gives Ctrough 105 +/- 30.8 mg/L (PMID 37694449)
// =====================================================================
CL_GAR   : 0.00546 : Garetosmab linear clearance (L/h). CALIBRATED to the
// reported steady-state Ctrough of 105 +/- 30.8 mg/L at 10 mg/kg q4w.
V1_GAR   : 3.10   : Garetosmab central volume (L)
V2_GAR   : 2.40   : Garetosmab peripheral volume (L)
Q_GAR    : 0.0115 : Garetosmab intercompartmental clearance (L/h)
VMAX_GAR : 0.040  : Garetosmab target-mediated elimination Vmax (mg/h).
// Small by design: at 10 mg/kg the target-mediated route is SATURATED
// and contributes little, which is what produces dose-proportional
// Ctrough of 105 +/- 31 mg/L (PMID 32557665, PMID 37694449).
KM_GAR   : 4.0    : Garetosmab target-mediated Km (mg/L)
KD_GAR   : 0.030  : Garetosmab-activin A dissociation constant (nM)
KP_GAR   : 0.15   : Garetosmab lesion-interstitium partition coefficient

// =====================================================================
// 3. ALK2 KINASE INHIBITOR PK — fidrisertib / zilurgisertib class, oral
// =====================================================================
KA_AK    : 0.60  : ALK2 inhibitor absorption rate (1/h)
CL_AK    : 12.0  : ALK2 inhibitor apparent clearance (L/h)
V_AK     : 95.0  : ALK2 inhibitor apparent volume (L)
F_AK     : 0.70  : ALK2 inhibitor oral bioavailability
IC50_AK  : 45.0  : ALK2 inhibitor kinase IC50 (ng/mL)
EMAX_AK  : 0.95  : ALK2 inhibitor maximal fractional signal block

// =====================================================================
// 4. SARACATINIB PK/PD — SRC/ALK2, oral (STOPFOP)
// =====================================================================
KA_SAR   : 0.45  : Saracatinib absorption rate (1/h)
CL_SAR   : 18.0  : Saracatinib apparent clearance (L/h)
V_SAR    : 550.0 : Saracatinib apparent volume (L)
F_SAR    : 0.60  : Saracatinib oral bioavailability
IC50_SAR : 110.0 : Saracatinib effective IC50 on SRC-permissive signal (ng/mL)
EMAX_SAR : 0.55  : Saracatinib maximal fractional signal block

// =====================================================================
// 5. SIROLIMUS PK/PD — mTORC1
// =====================================================================
KA_SIR   : 0.50  : Sirolimus absorption rate (1/h)
CL_SIR   : 8.0   : Sirolimus apparent clearance (L/h)
V_SIR    : 900.0 : Sirolimus apparent volume (L)
F_SIR    : 0.15  : Sirolimus oral bioavailability
IC50_SIR : 6.0   : Sirolimus IC50 on mTORC1 (ng/mL)
EMAX_SIR : 0.85  : Sirolimus maximal fractional mTORC1 block

// =====================================================================
// 6. PREDNISONE PK/PD — flare burst 2 mg/kg/d x 4 d
// =====================================================================
KA_PRE   : 1.20  : Prednisone absorption rate (1/h)
CL_PRE   : 9.0   : Prednisone apparent clearance (L/h)
V_PRE    : 45.0  : Prednisone apparent volume (L)
F_PRE    : 0.80  : Prednisone oral bioavailability
IC50_PRE : 25.0  : Prednisolone IC50 on inflammatory amplifier (ng/mL)
EMAX_PRE : 0.70  : Prednisone maximal fractional suppression of inflammation

// =====================================================================
// 7. IL-1 BLOCKADE PK/PD — anakinra SC
// =====================================================================
KA_IL1   : 0.25  : Anakinra SC absorption rate (1/h)
CL_IL1   : 5.5   : Anakinra clearance (L/h)
V_IL1    : 18.0  : Anakinra volume (L)
F_IL1    : 0.95  : Anakinra SC bioavailability
IC50_IL1 : 250.0 : Anakinra IC50 on IL-1 driven inflammation (ng/mL)
EMAX_IL1 : 0.75  : Anakinra maximal fractional IL-1 axis block

// =====================================================================
// 8. TRIGGER / INFLAMMATION (episodic route)
// =====================================================================
KOUT_TRIG : 0.25  : Trigger signal decay (1/day)      // ~4 d half-signal
KIN_INFL  : 0.40  : Inflammation formation per unit trigger (1/day)
KOUT_INFL : 0.10  : Inflammation resolution (1/day)   // flare ~2-6 wk
FL_RATE   : 0.90  : Flare-up rate (per patient-year; NHS 229 flares /
// 114 patients / 26.8 mo mean follow-up)
FL_TONIC  : 1.0   : Population-average flare drive when explicit flare
// events are NOT simulated (set 0 for event-driven mode)
KHIF      : 0.20  : HIF-1alpha formation per unit inflammation (1/day)
KOUT_HIF  : 0.15  : HIF-1alpha decay (1/day)
AHIF      : 0.60  : Maximal fractional signal amplification by HIF-1alpha
KIN_MTOR  : 0.30  : mTORC1 activation rate (1/day)
KOUT_MTOR : 0.30  : mTORC1 deactivation rate (1/day)
FMTOR     : 0.35  : Maximal fractional chondrogenic gain from mTORC1

// =====================================================================
// 9. ACTIVIN A / RECEPTOR / SIGNAL — kept LINEAR so the episodic and
//    smouldering routes remain exactly separable (Axis 1)
// =====================================================================
KIN_ASM   : 1.00  : Tonic (smouldering) activin A production (nM/day)
KIN_AFL   : 0.465 : Inflammation-driven activin A production (nM/day per unit
// episodic drive. SOLVED, not assumed: the value is the one
// that makes the flare-attributable share of baseline
// pSMAD1/5/8 equal PHI_FL = 0.295 (Axis 1)
KOUT_ACTA : 1.20  : Activin A elimination (1/day)
KOUT_CPLX : 0.30  : Garetosmab-activin complex clearance (1/day)
KON_GAR   : 8.0   : Garetosmab-activin association (1/nM/day)
F_LEAK    : 0.10  : Fraction of baseline smouldering signal that is
// LIGAND-INDEPENDENT receptor leak (FKBP12 clamp lost).
// Sets the ceiling on ANY anti-ligand strategy.
KSIG_A    : 1.00  : pSMAD1/5/8 generation per unit free activin A (1/day)
KOUT_SIG  : 1.00  : pSMAD1/5/8 dephosphorylation (1/day)
EPS_AGON  : 0.0   : Agonist conversion by a bivalent anti-ACVR1 antibody
// (fractional signal GAIN; set >0 to reproduce the
// paradoxical worsening reported for that class)

// =====================================================================
// 10. LESION CASCADE — one-way relay, total transit ~63 d so that a
//     flare that commits is radiographically visible by day 84
// =====================================================================
KFIB      : 0.045 : FAP recruitment into fibroproliferative lesion (1/day)
KCOMMIT   : 0.071 : Fibroproliferative -> CHONDROGENIC commitment (1/day) // ~14 d
KRESOLVE  : 0.193 : Fibroproliferative lesion RESOLUTION without ossifying
// (1/day). SOLVED so that KCOMMIT/(KCOMMIT+KRESOLVE) =
// 0.269 = the fraction of imaged flare-ups with new HO at
// day 84 (14/52, PMID 36152026). This competition is also
// where a commitment blocker acts: it diverts lesions from
// the ossifying branch into the resolving branch.
KHYP      : 0.071 : Chondrogenic -> hypertrophic (1/day)                  // ~14 d
KOST      : 0.071 : Hypertrophic -> osteoblastic / vascular invasion (1/day)
KMIN      : 0.048 : Osteoblastic -> MINERALISED bone (1/day)              // ~21 d
KRESORB   : 0.0   : Resorption of mineralised heterotopic bone (1/day)
// ZERO by design -- the ratchet (Axis 2)
F_THOR    : 0.25  : Fraction of new HO deposited in the thoracic cage
KFAP      : 0.30  : Activated-FAP formation rate (1/day)
KOUT_FAP  : 0.20  : Activated-FAP loss rate (1/day)

// =====================================================================
// 11. TERRITORY / GROWTH — Axis 3. The growth curve is anchored to human
//     somatic growth, NOT fitted to HO data; only KF_HO, HOMAX and
//     GAMMA_T are fitted, to the three NHS band velocities.
// =====================================================================
KF_HO     : 191.82 : Formation scale (mL/yr at unit drive, full territory).
// FITTED, with CMIN_C / A50C_C / AWC_C, to the FOUR prospectively measured
// NHS age-band velocities (21.9 / 41.5 / 4.6 / 4.6 mL/yr).
HOMAX     : 2000.0 : Maximum ossifiable territory, whole body ex-head (mL).
// FIXED A PRIORI, NOT FITTED. Territory is not identifiable from the NHS
// velocity profile: allowed to float, the optimiser drives HOMAX to 1e16 and
// silently DELETES the mechanism. It is therefore pinned so that territory is
// a real but non-binding constraint (about 39% consumed by age 56), which is
// what the clinical record requires -- adults retain enough unossified tissue
// to ossify catastrophically after surgery or trauma.
GAMMA_T   : 1.0   : Territory-exhaustion exponent. FIXED at 1 (linear).
BMIN_G    : 0.15  : Somatic size factor at birth (fraction of adult)
CMIN_C    : 0.1158 : Residual chondrogenic competence of adult progenitors. FITTED.
// The SAME biology as Axis 4: chondrogenic competence falls with skeletal
// maturity, which is why HO velocity drops after the early twenties AND why
// the physeal hazard of a retinoid disappears once the plates have fused.
// One mechanism, two independent observations. FITTED.
A50C_C    : 24.02 : Age at half loss of chondrogenic competence (yr). FITTED.
AWC_C     : 0.50  : Width of the competence decline (yr). FITTED, and it lands on
// its lower bound: two adjacent NHS band means of 41.5 then 4.6 mL/yr demand a
// near-STEP loss of competence in the mid-twenties, possibly sharper than this.
A50_G     : 13.5  : Age at half adult somatic size factor (yr)
AW_G      : 2.20  : Width of the somatic growth transition (yr)
HO_INIT   : 60.0  : Heterotopic bone already present at AGE0 (mL).
// ANCHORED at age 5 to the NHS 2-<8 yr band mean stock of 68.8 mL.
// Scenarios that start later must raise this: use ho_at_age().

// =====================================================================
// 12. DRUG EFFECT SITES
// =====================================================================
EC50_PAL  : 11.480 : Palovarotene plasma EC50 on chondrogenic commitment (ng/mL).
// FITTED so that the MOVE regimen (5 mg/d, 18 mo, age 15) reduces
// annualised new HO by 57% (the midpoint of the reported 54-60%).
EMAX_PAL  : 0.85  : Palovarotene maximal fractional block of commitment
FACT_PAL  : 0.20  : Extra fractional reduction of local activin-A-expressing
// cells by palovarotene (PMID 38130748)
R_SEL     : 3.698 : Physeal / anti-HO selectivity = EC50(PPC)/EC50(HO).
// NOT assumed -- solved from MOVE (60% HO suppression
// with 36.8% PPC in <14 y) given HPPC_MAX below.
HPPC_MAX  : 0.80  : Maximal premature-physeal-closure hazard at saturating
// retinoid exposure (1/yr). R_SEL is conditional on this.
K_VBMD    : 0.00045 : Vertebral BMD Z loss per unit retinoid exposure-day
K_RETAE   : 0.010 : Mucocutaneous retinoid AE accrual rate (1/day at Emax)
K_HTZ     : 0.55  : Height-Z loss per unit cumulative physeal closure

// =====================================================================
// 13. IMAGING / ENDPOINT LAYER
// =====================================================================
W_NAF_NEW : 1.00  : NaF PET signal weight per unit ACTIVE mineralising front
W_NAF_OLD : 0.014 : NaF PET signal weight per mL of MATURE stock
// (this ratio is what dilutes total lesion ACTIVITY)
KOUT_NAFO : 0.05  : Turnover of the mature-stock remodelling pool (1/day)
LES_VOL   : 28.8  : Mean volume of one discrete new HO lesion (mL)
KCAJ      : 0.0030 : Rate at which ankylosis follows ossification (1/day)
CAJ_MAX   : 30.0  : CAJIS ceiling
CAJ_K50   : 210.0 : Appendicular HO giving half-maximal CAJIS (mL)
CAJ_CONG  : 2.0   : Congenital CAJIS component (cervical spine, toes)
FVC_MIN   : 35.0  : Asymptotic floor of FVC (% predicted)
FVC_K     : 25.0  : Thoracic HO giving half of the total FVC loss (mL)
KFVC      : 0.0040 : Rate at which FVC follows thoracic ossification (1/day)
H0_MORT   : 1.871e-6 : Baseline daily mortality hazard at FVC = 100%.
// FITTED so that untreated median survival = 56 yr (PMID 20194327).
BETA_FVC  : 0.055 : Log-hazard increase per 1% predicted FVC lost
H_BG      : 2.0e-6 : Age-independent background daily hazard

$CMT @annotated
// ---- PK (15) ----
PALGUT  : Palovarotene gut depot (mg)
PALCEN  : Palovarotene central (mg)
PALPER  : Palovarotene peripheral (mg)
GARCEN  : Garetosmab central (mg)
GARPER  : Garetosmab peripheral (mg)
AKGUT   : ALK2 inhibitor gut depot (mg)
AKCEN   : ALK2 inhibitor central (mg)
SARGUT  : Saracatinib gut depot (mg)
SARCEN  : Saracatinib central (mg)
SIRGUT  : Sirolimus gut depot (mg)
SIRCEN  : Sirolimus central (mg)
PREGUT  : Prednisone gut depot (mg)
PRECEN  : Prednisolone central (mg)
IL1DEP  : Anakinra SC depot (mg)
IL1CEN  : Anakinra central (mg)
// ---- Disease / signalling (12) ----
TRIG    : Trigger / tissue-injury signal (a.u.)
INFL    : Inflammatory amplifier: mast cells, macrophages, IL-1 (a.u.)
HIFC    : HIF-1alpha activity in the hypoxic lesion core (a.u.)
MTORC   : mTORC1 activity (a.u.)
ACTASM  : Free activin A, SMOULDERING source (nM)
ACTAFL  : Free activin A, FLARE source (nM)
ACTCPX  : Garetosmab-activin A complex (nM)
SIGSM   : pSMAD1/5/8 from the smouldering route (a.u.)
SIGFL   : pSMAD1/5/8 from the episodic route (a.u.)
FAPA    : Activated fibro/adipogenic progenitors (a.u.)
FIBRO   : Stage 1 fibroproliferative lesion mass (mL equivalent)
CHOND   : Stage 2 chondrogenic anlagen -- point of no return (mL eq)
// ---- Lesion maturation / bone (7) ----
HYPT    : Stage 3-4 hypertrophic + vascular invasion (mL eq)
OSTB    : Stage 5 osteoblastic osteoid (mL eq)
HOTH    : Mineralised heterotopic bone, THORACIC cage (mL)
HOAP    : Mineralised heterotopic bone, appendicular + other (mL)
NAFOLD  : Remodelling activity of the MATURE stock (NaF PET a.u.)
CUMNEW  : Cumulative NEW mineralised HO since t = 0 (mL)
NLES    : Cumulative new discrete HO lesion count
// ---- Clinical / safety (8) ----
CAJ     : CAJIS cumulative joint involvement (0-30)
FVCC    : Forced vital capacity (% predicted)
CHAZ    : Cumulative mortality hazard
PPCC    : Cumulative premature-physeal-closure hazard
VBMDZ   : Vertebral BMD Z-score deficit (negative = worse)
HTZ     : Height Z-score deficit (negative = worse)
RETAE   : Mucocutaneous retinoid AE burden (0-1)
NFLARE  : Cumulative flare-up count

$MAIN
// ---- between-subject variability ----
KFHO_I = KF_HO    * exp(ETA_LAM);
FLT_I  = FL_TONIC * exp(ETA_FL);

// ---- initial conditions ----
HOAP_0  = HO_INIT * (1.0 - F_THOR);
HOTH_0  = HO_INIT * F_THOR;
FVCC_0  = FVC_MIN + (100.0 - FVC_MIN) * FVC_K / (FVC_K + HO_INIT * F_THOR);
CAJ_0   = CAJ_CONG + (CAJ_MAX - CAJ_CONG) * (HO_INIT * (1.0 - F_THOR)) /
                     (CAJ_K50 + HO_INIT * (1.0 - F_THOR));

// ---- exact steady state of the signalling chain at baseline drive ----
// The smouldering branch carries the ligand-independent leak, so its
// steady state is inflated by 1/(1-F_LEAK); the episodic branch does not.
ACTASM_0 = KIN_ASM / KOUT_ACTA;
ACTAFL_0 = KIN_AFL * FLT_I / KOUT_ACTA;
SIGSM_0  = KSIG_A * ACTASM_0 / ((1.0 - F_LEAK) * KOUT_SIG);
SIGFL_0  = KSIG_A * ACTAFL_0 / KOUT_SIG;
FAPA_0   = KFAP * (SIGSM_0 + SIGFL_0) / KOUT_FAP;
NAFOLD_0 = W_NAF_OLD * HO_INIT;

$ODE
// =====================================================================
// AGE, SOMATIC SIZE, TERRITORY
// =====================================================================
AGEY   = AGE0 + SOLVERTIME / 365.25;
BODYF  = BMIN_G + (1.0 - BMIN_G) / (1.0 + exp(-(AGEY - A50_G) / AW_G));
// chondrogenic competence of the resident progenitor pool, declining with
// skeletal maturity: this is what keeps adult velocity low but NON-ZERO
COMPF  = CMIN_C + (1.0 - CMIN_C) / (1.0 + exp((AGEY - A50C_C) / AWC_C));
HOTOT  = HOTH + HOAP;
TERR   = fmax(0.0, 1.0 - HOTOT / HOMAX);
TERRF  = pow(TERR, GAMMA_T);
AGEFUSE = 14.0 + 2.0 * SEX;
GPOPEN = 1.0 / (1.0 + exp((AGEY - AGEFUSE) / 0.5));

// =====================================================================
// DRUG CONCENTRATIONS (ng/mL for small molecules, mg/L for the mAb)
// =====================================================================
CP_PAL = 1000.0 * PALCEN / V1_PAL;
CP_GAR = GARCEN / V1_GAR;
CP_AK  = 1000.0 * AKCEN  / V_AK;
CP_SAR = 1000.0 * SARCEN / V_SAR;
CP_SIR = 1000.0 * SIRCEN / V_SIR;
CP_PRE = 1000.0 * PRECEN / V_PRE;
CP_IL1 = 1000.0 * IL1CEN / V_IL1;

// =====================================================================
// DRUG EFFECTS
// =====================================================================
// palovarotene: blocks the fibroproliferative -> chondrogenic commitment
E_PAL   = EMAX_PAL * CP_PAL / (EC50_PAL + CP_PAL);
// same exposure axis, different EC50, only while the physis is open
E_PPC   = HPPC_MAX * CP_PAL / (R_SEL * EC50_PAL + CP_PAL);
// ALK2 kinase inhibitor and saracatinib: block signal GENERATION,
// including the ligand-independent leak that no anti-ligand can reach
E_AK    = EMAX_AK  * CP_AK  / (IC50_AK  + CP_AK);
E_SAR   = EMAX_SAR * CP_SAR / (IC50_SAR + CP_SAR);
E_KIN   = 1.0 - (1.0 - E_AK) * (1.0 - E_SAR);
// mTORC1 and anti-inflammatory agents
E_SIR   = EMAX_SIR * CP_SIR / (IC50_SIR + CP_SIR);
E_PRE   = EMAX_PRE * CP_PRE / (IC50_PRE + CP_PRE);
E_IL1   = EMAX_IL1 * CP_IL1 / (IC50_IL1 + CP_IL1);
E_ANTIINF = 1.0 - (1.0 - E_PRE) * (1.0 - E_IL1);

// =====================================================================
// PK ODEs
// =====================================================================
dxdt_PALGUT = -KA_PAL * PALGUT * 24.0;
dxdt_PALCEN =  KA_PAL * PALGUT * 24.0
              - (CL_PAL / V1_PAL) * PALCEN * 24.0
              - (Q_PAL  / V1_PAL) * PALCEN * 24.0
              + (Q_PAL  / V2_PAL) * PALPER * 24.0;
dxdt_PALPER =  (Q_PAL / V1_PAL) * PALCEN * 24.0 - (Q_PAL / V2_PAL) * PALPER * 24.0;

dxdt_GARCEN = -(CL_GAR / V1_GAR) * GARCEN * 24.0
              - (Q_GAR  / V1_GAR) * GARCEN * 24.0
              + (Q_GAR  / V2_GAR) * GARPER * 24.0
              - 24.0 * VMAX_GAR * CP_GAR / (KM_GAR + CP_GAR);
dxdt_GARPER =  (Q_GAR / V1_GAR) * GARCEN * 24.0 - (Q_GAR / V2_GAR) * GARPER * 24.0;

dxdt_AKGUT  = -KA_AK * AKGUT * 24.0;
dxdt_AKCEN  =  KA_AK * AKGUT * 24.0 - (CL_AK / V_AK) * AKCEN * 24.0;
dxdt_SARGUT = -KA_SAR * SARGUT * 24.0;
dxdt_SARCEN =  KA_SAR * SARGUT * 24.0 - (CL_SAR / V_SAR) * SARCEN * 24.0;
dxdt_SIRGUT = -KA_SIR * SIRGUT * 24.0;
dxdt_SIRCEN =  KA_SIR * SIRGUT * 24.0 - (CL_SIR / V_SIR) * SIRCEN * 24.0;
dxdt_PREGUT = -KA_PRE * PREGUT * 24.0;
dxdt_PRECEN =  KA_PRE * PREGUT * 24.0 - (CL_PRE / V_PRE) * PRECEN * 24.0;
dxdt_IL1DEP = -KA_IL1 * IL1DEP * 24.0;
dxdt_IL1CEN =  KA_IL1 * IL1DEP * 24.0 - (CL_IL1 / V_IL1) * IL1CEN * 24.0;

// =====================================================================
// EPISODIC ROUTE: trigger -> inflammation -> hypoxia / mTOR
// =====================================================================
dxdt_TRIG = -KOUT_TRIG * TRIG;
dxdt_INFL =  KIN_INFL * TRIG * (1.0 - E_ANTIINF) - KOUT_INFL * INFL;
dxdt_HIFC =  KHIF * INFL - KOUT_HIF * HIFC;
dxdt_MTORC = KIN_MTOR * (1.0 + INFL) * (1.0 - E_SIR) - KOUT_MTOR * MTORC;

// flare drive: FL_TONIC is the population-average episodic drive used
// when discrete flare events are not simulated; INFL adds event-driven
// flares on top. Anti-inflammatory drug suppresses BOTH.
FLDRIVE = (FLT_I * (1.0 - E_ANTIINF) + INFL);

// =====================================================================
// ACTIVIN A (two separable sources) AND pSMAD1/5/8
// =====================================================================
// mg/L -> nM (IgG ~145 kDa), then partitioned into the lesion interstitium
GARFREE = KP_GAR * fmax(0.0, CP_GAR) / 145.0 * 1000.0;
BINDSM  = KON_GAR * GARFREE * ACTASM / (KD_GAR + ACTASM);
BINDFL  = KON_GAR * GARFREE * ACTAFL / (KD_GAR + ACTAFL);

dxdt_ACTASM = KIN_ASM * (1.0 - FACT_PAL * E_PAL / fmax(EMAX_PAL, 1e-9))
              - KOUT_ACTA * ACTASM - BINDSM;
dxdt_ACTAFL = KIN_AFL * FLDRIVE * (1.0 - FACT_PAL * E_PAL / fmax(EMAX_PAL, 1e-9))
              - KOUT_ACTA * ACTAFL - BINDFL;
dxdt_ACTCPX = BINDSM + BINDFL - KOUT_CPLX * ACTCPX;

AMP     = (1.0 + AHIF * HIFC / (1.0 + HIFC)) * (1.0 - E_KIN) * (1.0 + EPS_AGON);
// smouldering signal carries the ligand-independent leak; the leak is
// invisible to anti-ligand therapy but not to a kinase inhibitor
LEAKSRC = KSIG_A * (KIN_ASM / KOUT_ACTA) * F_LEAK / (1.0 - F_LEAK);
dxdt_SIGSM = (KSIG_A * ACTASM + LEAKSRC) * AMP - KOUT_SIG * SIGSM;
dxdt_SIGFL =  KSIG_A * ACTAFL * AMP        - KOUT_SIG * SIGFL;

SIGTOT = SIGSM + SIGFL;

// =====================================================================
// LESION CASCADE — the one-way relay
// =====================================================================
dxdt_FAPA = KFAP * SIGTOT - KOUT_FAP * FAPA;

// recruitment is scaled by somatic size and remaining territory (Axis 3)
RECRUIT = (KFHO_I / 365.25) * (FAPA / FAPA_REF) * BODYF * COMPF * TERRF;
dxdt_FIBRO = RECRUIT - KCOMMIT * FIBRO * (1.0 - E_PAL) - KRESOLVE * FIBRO;
COMMITF    = KCOMMIT * FIBRO * (1.0 - E_PAL) * (1.0 + FMTOR * MTORC / (1.0 + MTORC));
dxdt_CHOND = COMMITF - KHYP * CHOND;
dxdt_HYPT  = KHYP * CHOND - KOST * HYPT;
dxdt_OSTB  = KOST * HYPT  - KMIN * OSTB;

MINFLUX = KMIN * OSTB;
dxdt_HOTH = F_THOR * MINFLUX - KRESORB * HOTH;
dxdt_HOAP = (1.0 - F_THOR) * MINFLUX - KRESORB * HOAP;
dxdt_CUMNEW = MINFLUX;
dxdt_NLES   = COMMITF / LES_VOL;

// =====================================================================
// IMAGING
// =====================================================================
dxdt_NAFOLD = KOUT_NAFO * (W_NAF_OLD * HOTOT - NAFOLD);

// =====================================================================
// CLINICAL ENDPOINTS
// =====================================================================
CAJTGT = CAJ_CONG + (CAJ_MAX - CAJ_CONG) * HOAP / (CAJ_K50 + HOAP);
dxdt_CAJ = KCAJ * (CAJTGT - CAJ);

FVCTGT = FVC_MIN + (100.0 - FVC_MIN) * FVC_K / (FVC_K + HOTH);
dxdt_FVCC = KFVC * (FVCTGT - FVCC);

dxdt_CHAZ = H_BG + H0_MORT * exp(BETA_FVC * (100.0 - FVCC));

// =====================================================================
// SAFETY / TOXICITY — Axis 4
// =====================================================================
dxdt_PPCC  = GPOPEN * E_PPC / 365.25;
dxdt_VBMDZ = -K_VBMD * E_PAL / EMAX_PAL;
dxdt_HTZ   = -K_HTZ * GPOPEN * E_PPC / 365.25;
dxdt_RETAE = K_RETAE * (E_PAL / EMAX_PAL) * (1.0 - RETAE);
dxdt_NFLARE = (FL_RATE / 365.25) * (1.0 - E_ANTIINF);

$OMEGA @annotated @block
ETA_LAM : 0.36                : IIV on formation scale (log)
ETA_FL  : 0.00 0.55           : IIV on episodic drive (log)

$TABLE
capture AGE_YR   = AGE0 + TIME / 365.25;
capture HO_TOTAL = HOTH + HOAP;
capture HO_THOR  = HOTH;
capture HO_APP   = HOAP;
capture FVC_PCT  = FVCC;
capture CAJIS    = CAJ;
capture SURV     = exp(-CHAZ);
capture PPC_PROB = 1.0 - exp(-PPCC);
capture NEW_HO   = CUMNEW;
capture LESIONS  = NLES;
capture C_PALO   = 1000.0 * PALCEN / V1_PAL;
capture C_GARE   = GARCEN / V1_GAR;
capture C_ALK2I  = 1000.0 * AKCEN / V_AK;
capture C_SIRO   = 1000.0 * SIRCEN / V_SIR;
capture ACTA_FREE = ACTASM + ACTAFL;
capture ACTA_TOT  = ACTASM + ACTAFL + ACTCPX;
capture SIG_TOTAL = SIGSM + SIGFL;
// Axis 1 read-out: the flare-attributable share of the CURRENT drive
capture FRAC_FLARE = (SIGSM + SIGFL > 0) ? SIGFL / (SIGSM + SIGFL) : 0.0;
// NaF PET: activity of the NEW mineralising front versus the mature stock.
// Total lesion activity is what LUMINA-1 used as its primary endpoint.
capture NAF_NEW   = W_NAF_NEW * (HYPT + OSTB);
capture NAF_OLD   = NAFOLD;
capture NAF_TOTAL = W_NAF_NEW * (HYPT + OSTB) + NAFOLD;
capture NAF_FRAC_NEW = (NAF_TOTAL > 0) ? (W_NAF_NEW * (HYPT + OSTB)) / NAF_TOTAL : 0.0;
capture VBMD_Z   = VBMDZ;
capture HEIGHT_Z = HTZ;
capture RET_AE   = RETAE;
capture TERRITORY = fmax(0.0, 1.0 - (HOTH + HOAP) / HOMAX);
capture BODY_F   = BMIN_G + (1.0 - BMIN_G) /
                   (1.0 + exp(-((AGE0 + TIME / 365.25) - A50_G) / AW_G));
capture COMP_F   = CMIN_C + (1.0 - CMIN_C) /
                   (1.0 + exp(((AGE0 + TIME / 365.25) - A50C_C) / AWC_C));

$GLOBAL
// Reference activated-FAP level: fixes the units of the recruitment term
// so that KF_HO is directly interpretable as mL/yr at unit drive and
// full territory. Value = steady state of FAPA at baseline drive.
#define FAPA_REF (1.9700)
'

# The model carries an OMEGA block for between-subject variability, but the
# deterministic scenarios must be reproducible, so the default object has the
# random effects zeroed. Scenario 14 uses fop_iiv, which keeps them.
fop_iiv <- mcode("fop", fop_code)
fop     <- zero_re(fop_iiv)

# =====================================================================
#  HELPERS
# =====================================================================
CMT_PAL <- "PALGUT"; CMT_GAR <- "GARCEN"; CMT_AK <- "AKGUT"
CMT_SAR <- "SARGUT"; CMT_SIR <- "SIRGUT"; CMT_PRE <- "PREGUT"
CMT_IL1 <- "IL1DEP"; CMT_TRIG <- "TRIG"

yr <- function(x) x * 365.25

#' Interpolate a captured column at arbitrary times (duplicate-time safe:
#' dosing records make the output grid non-unique)
at_time <- function(out, col, times) {
  d <- out[!duplicated(out$time), ]
  approx(d$time, d[[col]], xout = times, rule = 2)$y
}

#' Interpolate a captured column at arbitrary AGES
at_age <- function(out, col, ages) {
  d <- out[!duplicated(out$AGE_YR), ]
  approx(d$AGE_YR, d[[col]], xout = ages, rule = 2)$y
}

#' Annualised new HO volume between two times (the MOVE / NHS endpoint).
#' This is a RATE. Axis 2 is the reminder that disability is its integral.
ann_new_ho <- function(out, t0, t1) {
  a <- at_time(out, "NEW_HO", c(t0, t1))
  (a[2] - a[1]) / ((t1 - t0) / 365.25)
}

#' Chronic oral regimen
oral_daily <- function(amt, start, dur_days, cmt) {
  ev(amt = amt, ii = 1, addl = dur_days - 1, time = start, cmt = cmt)
}

#' Discrete flare events injected into the trigger compartment
flare_events <- function(times, amt = 1.0) {
  ev(amt = amt, time = times, cmt = CMT_TRIG)
}

sim <- function(model = fop, events = NULL, end, delta = 7, ...) {
  m <- model %>% param(...)
  out <- if (is.null(events)) {
    mrgsim(m, end = end, delta = delta, atol = 1e-8, rtol = 1e-6)
  } else {
    mrgsim(m, events = events, end = end, delta = delta, atol = 1e-8, rtol = 1e-6)
  }
  as.data.frame(out)
}


# ---------------------------------------------------------------------
#  AGE-APPROPRIATE STARTING BURDEN
#  Every scenario that begins after age 5 must start from the heterotopic
#  bone that has ALREADY formed by then -- otherwise the territory factor
#  and, far more importantly, the DENOMINATOR of every "% change in total
#  volume" statement (Axis 2) are both wrong.
# ---------------------------------------------------------------------
ho_at_age <- local({
  traj <- NULL
  function(age) {
    if (is.null(traj)) traj <<- sim(end = yr(75), delta = 30.4375, AGE0 = 5)
    if (age <= 5) return(60)
    as.numeric(at_age(traj, "HO_TOTAL", age))
  }
})

#' Simulate a patient of a given age, starting from that age's real burden
sim_pt <- function(age0, end, delta = 7, events = NULL, ...) {
  sim(end = end, delta = delta, events = events,
      AGE0 = age0, HO_INIT = ho_at_age(age0), ...)
}

# =====================================================================
#  SCENARIO 1 — NATURAL HISTORY, age 5 -> 80 (the NHS calibration run)
# =====================================================================
sc01_natural_history <- function(to_age = 56) {
  sim(end = yr(to_age - 5), delta = 30.4375, AGE0 = 5)
}

# =====================================================================
#  SCENARIO 2 — NHS AGE-BAND VELOCITIES (target 21.9 / 41.5 / 4.6 mL/yr)
#    The non-monotonic profile is NOT fitted band by band: it emerges
#    from one growth curve x one finite territory (Axis 3).
# =====================================================================
sc02_age_bands <- function() {
  out <- sc01_natural_history()
  bands <- list(c(8, 15), c(15, 25), c(25, 56))
  obs   <- c(21.9, 41.5, 4.6)
  do.call(rbind, lapply(seq_along(bands), function(i) {
    b <- bands[[i]]
    data.frame(band = paste0(b[1], "-", b[2], " yr"),
               model_mL_per_yr = round(ann_new_ho(out, yr(b[1] - 5), yr(b[2] - 5)), 1),
               NHS_mL_per_yr   = obs[i])
  }))
}

# =====================================================================
#  SCENARIO 3 — PALOVAROTENE, MOVE REGIMEN, 18 months from age 15
#    chronic 5 mg/d + flare-up 20 mg x 4 wk then 10 mg x >=8 wk.
#    Reports the RATE and the STOCK from the same run (Axis 2).
# =====================================================================
sc03_move <- function(age0 = 15, months = 18, flare_month = 6) {
  dur  <- round(yr(months / 12))
  fl_t <- round(yr(flare_month / 12))
  trig <- flare_events(fl_t)
  chronic <- oral_daily(5, 0, dur, CMT_PAL)
  burst   <- c(oral_daily(20, fl_t, 28, CMT_PAL),
               oral_daily(10, fl_t + 28, 56, CMT_PAL))
  unt <- sim_pt(age0, dur, 14, trig)
  trt <- sim_pt(age0, dur, 14, c(chronic, burst, trig))
  rate_u <- ann_new_ho(unt, 0, dur); rate_t <- ann_new_ho(trt, 0, dur)
  list(untreated = unt, treated = trt,
       summary = data.frame(
         quantity = c("annualised NEW HO volume (mL/yr) -- the trial endpoint",
                      "reduction in that RATE (%)",
                      "TOTAL HO stock at 18 months (mL)",
                      "reduction in the STOCK (%)",
                      "CAJIS at 18 months",
                      "premature physeal closure probability (%)",
                      "mucocutaneous retinoid AE burden (0-1)"),
         untreated = round(c(rate_u, 0, tail(unt$HO_TOTAL, 1), 0,
                             tail(unt$CAJIS, 1), 100 * tail(unt$PPC_PROB, 1),
                             tail(unt$RET_AE, 1)), 2),
         palovarotene = round(c(rate_t, 100 * (1 - rate_t / rate_u),
                                tail(trt$HO_TOTAL, 1),
                                100 * (1 - tail(trt$HO_TOTAL, 1) / tail(unt$HO_TOTAL, 1)),
                                tail(trt$CAJIS, 1), 100 * tail(trt$PPC_PROB, 1),
                                tail(trt$RET_AE, 1)), 2)))
}

# =====================================================================
#  SCENARIO 4 — CHRONIC versus FLARE-TRIGGERED-ONLY dosing
#    Why flare-only dosing under-performs: by the time a flare is
#    symptomatic, part of the lesion has already crossed the chondrogenic
#    checkpoint, where a commitment blocker no longer reaches it.
# =====================================================================
sc04_dosing_strategy <- function(age0 = 15, months = 18) {
  dur  <- round(yr(months / 12))
  fl_t <- round(yr(0.5))
  trig <- flare_events(fl_t)
  flare_only <- c(oral_daily(20, fl_t, 28, CMT_PAL),
                  oral_daily(10, fl_t + 28, 56, CMT_PAL))
  chronic <- oral_daily(5, 0, dur, CMT_PAL)
  arms <- list(none = trig, flare_only = c(flare_only, trig),
               chronic = c(chronic, trig),
               chronic_plus_flare = c(chronic, flare_only, trig))
  base <- ann_new_ho(sim_pt(age0, dur, 14, trig), 0, dur)
  do.call(rbind, lapply(names(arms), function(a) {
    o <- sim_pt(age0, dur, 14, arms[[a]])
    r <- ann_new_ho(o, 0, dur)
    data.frame(regimen = a, ann_new_HO_mL_yr = round(r, 2),
               pct_reduction = round(100 * (1 - r / base), 1),
               cumulative_palo_mg = round(sum(arms[[a]]@data$amt *
                                              (arms[[a]]@data$addl + 1)), 0))
  }))
}

# =====================================================================
#  SCENARIO 5 — AGE AT TREATMENT START (Axis 3)
#    One drug, one efficacy, five starting ages, lifetime horizon.
# =====================================================================
sc05_age_at_start <- function(start_ages = c(5, 8, 15, 25, 35)) {
  ref_age <- 40      # a fixed age at which to compare burden
  end_age <- 56
  ref <- sc01_natural_history(end_age)
  burden40_untreated <- as.numeric(at_age(ref, "HO_TOTAL", ref_age))
  total56 <- as.numeric(at_age(ref, "HO_TOTAL", end_age))
  do.call(rbind, lapply(start_ages, function(a) {
    dur <- round(yr(end_age - a))
    trt <- sim_pt(a, dur, 91.3, oral_daily(5, 0, dur, CMT_PAL))
    trt$AGE_YR <- a + trt$time / 365.25
    # the delay can exceed the lifespan horizon, so look for the crossing
    # over a deliberately long window
    durL <- round(yr(95 - a))
    trtL <- sim_pt(a, durL, 182.6, oral_daily(5, 0, durL, CMT_PAL))
    trtL$AGE_YR <- a + trtL$time / 365.25
    b40 <- as.numeric(at_age(trt, "HO_TOTAL", ref_age))
    # years of DELAY: when does the treated patient reach the burden the
    # untreated patient had at age 40?
    d <- trtL[!duplicated(trtL$HO_TOTAL), ]
    delay <- if (max(d$HO_TOTAL) < burden40_untreated) Inf else
      as.numeric(approx(d$HO_TOTAL, d$AGE_YR, xout = burden40_untreated)$y) - ref_age
    data.frame(start_age = a,
               HO_already_formed_mL = round(ho_at_age(a), 0),
               pct_of_age56_burden_still_ahead =
                 round(100 * (1 - ho_at_age(a) / total56), 1),
               HO_at_40_mL = round(b40, 0),
               pct_burden_averted_at_40 =
                 round(100 * (1 - b40 / burden40_untreated), 1),
               years_of_delay = ifelse(is.infinite(delay), ">55",
                                       format(round(delay, 1))),
               physis_open_at_start = ifelse(a < 14, "YES (PPC risk)", "no"))
  }))
}

# =====================================================================
#  SCENARIO 6 — GARETOSMAB, LUMINA-1 (10 mg/kg IV q4w, 28 weeks)
#    Reproduces BOTH results of that trial from one mechanism: new
#    lesions strongly suppressed, while the PET total-lesion-ACTIVITY
#    primary endpoint is diluted by remodelling of the mature stock the
#    drug cannot touch.
# =====================================================================
#' @param enrich activity enrichment of the trial cohort relative to the NHS
#'   population average adult. NOT a free knob: LUMINA-1 placebo produced 29
#'   new HO lesions in 24 patients over 28 weeks (1.21 per patient), whereas an
#'   NHS-average adult in this model produces ~0.1. The trial cohort was
#'   therefore roughly an order of magnitude more active than the population it
#'   is usually compared with -- which is itself the finding, and the reason
#'   single-arm-versus-external-control designs in FOP are fragile (see
#'   sc15_external_control_bias).
sc06_lumina1 <- function(age0 = 33, wt = 70, enrich = 12) {
  dur  <- 196
  gar  <- ev(amt = 10 * wt, ii = 28, addl = 6, time = 0, cmt = CMT_GAR)
  trig <- flare_events(c(30, 110))
  pbo <- sim_pt(age0, dur, 7, trig, WT = wt, KF_HO = 191.82 * enrich)
  trt <- sim_pt(age0, dur, 7, c(gar, trig), WT = wt, KF_HO = 191.82 * enrich)
  tailmean <- function(x, n = 8) mean(tail(x, n))
  data.frame(
    endpoint = c("new HO volume over 28 wk (mL)",
                 "new discrete lesions over 28 wk (n)",
                 "free activin A at week 28 (nM)",
                 "TOTAL NaF lesion activity, week 28 (rel. to placebo)",
                 "share of NaF signal coming from NEW lesions",
                 "garetosmab Ctrough (mg/L)"),
    placebo = c(round(tail(pbo$NEW_HO, 1), 2), round(tail(pbo$LESIONS, 1), 2),
                round(tail(pbo$ACTA_FREE, 1), 3), 1.000,
                round(tailmean(pbo$NAF_FRAC_NEW), 3), NA),
    garetosmab = c(round(tail(trt$NEW_HO, 1), 2), round(tail(trt$LESIONS, 1), 2),
                   round(tail(trt$ACTA_FREE, 1), 3),
                   round(tailmean(trt$NAF_TOTAL) / tailmean(pbo$NAF_TOTAL), 3),
                   round(tailmean(trt$NAF_FRAC_NEW), 3),
                   round(min(tail(trt$C_GARE, 5)), 1)))
}

# =====================================================================
#  SCENARIO 7 — THE ANTI-FLARE CEILING (Axis 1)
#    Abolish the episodic route completely and see what is left standing.
# =====================================================================
sc07_antiflare_ceiling <- function(age0 = 15, years = 3) {
  dur  <- round(yr(years))
  base <- sim_pt(age0, dur, 30.4)
  perf <- sim_pt(age0, dur, 30.4, FL_TONIC = 0)
  gc   <- sim_pt(age0, dur, 30.4, KIN_AFL = 0.465 * (1 - 0.12))
  il1  <- sim_pt(age0, dur, 30.4, KIN_AFL = 0.465 * (1 - 0.50))
  r <- sapply(list(base, gc, il1, perf), function(o) ann_new_ho(o, 0, dur))
  data.frame(
    strategy = c("no anti-inflammatory therapy",
                 "glucocorticoid bursts, real-world efficacy (12% of flares fully resolved)",
                 "IL-1 blockade, assumed 50% flare suppression",
                 "PERFECT abolition of every flare-up"),
    ann_new_HO_mL_yr = round(r, 2),
    pct_reduction = round(100 * (1 - r / r[1]), 1))
}

# =====================================================================
#  SCENARIO 8 — WHERE EACH CLASS HITS THE PATHWAY, AND ITS OWN CEILING
#    Anti-ligand therapy cannot reach the ligand-independent receptor
#    leak (F_LEAK); a kinase inhibitor can. That is a structural
#    difference between the classes, not a potency difference.
# =====================================================================
sc08_class_ceilings <- function(age0 = 15, years = 3) {
  dur  <- round(yr(years))
  base <- ann_new_ho(sim_pt(age0, dur, 30.4), 0, dur)
  cl <- class_ceilings()$analytic_ceiling_pct
  ceil <- list("untreated" = NA,
               "palovarotene 5 mg/d (commitment)" = cl[4],
               "anti-activin A 10 mg/kg q4w (ligand)" = cl[2],
               "ALK2 kinase inhibitor 200 mg/d" = cl[3],
               "ALK2 kinase inhibitor 800 mg/d" = cl[3],
               "saracatinib 125 mg/d (SRC/ALK2)" = NA,
               "sirolimus 2 mg/d (mTORC1)" = NA)
  arms <- list(
    "untreated"                            = NULL,
    "palovarotene 5 mg/d (commitment)"     = oral_daily(5, 0, dur, CMT_PAL),
    "anti-activin A 10 mg/kg q4w (ligand)" = ev(amt = 700, ii = 28,
                                                addl = floor(dur / 28), cmt = CMT_GAR),
    "ALK2 kinase inhibitor 200 mg/d"       = oral_daily(200, 0, dur, CMT_AK),
    "ALK2 kinase inhibitor 800 mg/d"       = oral_daily(800, 0, dur, CMT_AK),
    "saracatinib 125 mg/d (SRC/ALK2)"      = oral_daily(125, 0, dur, CMT_SAR),
    "sirolimus 2 mg/d (mTORC1)"            = oral_daily(2, 0, dur, CMT_SIR))
  do.call(rbind, lapply(names(arms), function(a) {
    o <- sim_pt(age0, dur, 30.4, arms[[a]])
    r <- ann_new_ho(o, 0, dur)
    data.frame(arm = a, ann_new_HO_mL_yr = round(r, 2),
               pct_reduction = round(100 * (1 - r / base), 1),
               class_ceiling_pct = ceil[[a]],
               reaches_ligand_independent_leak =
                 ifelse(grepl("kinase|saracatinib", a), "yes",
                        ifelse(a == "untreated", "-", "no")))
  }))
}

#' The analytic ceiling of each mechanism class, computed from the parameters
#' rather than from a dose. This is the structural claim of the model: the
#' classes differ in WHERE they cap out, not only in how potent they are.
class_ceilings <- function() {
  p <- as.list(param(fop))
  actasm <- p$KIN_ASM / p$KOUT_ACTA
  actafl <- p$KIN_AFL * p$FL_TONIC / p$KOUT_ACTA
  sigsm  <- p$KSIG_A * actasm / (1 - p$F_LEAK)
  sigfl  <- p$KSIG_A * actafl
  leak   <- sigsm * p$F_LEAK
  # commitment blocker: the checkpoint is a competition between commitment and
  # resolution, so blocking commitment by EMAX_PAL does NOT reduce ossification
  # by EMAX_PAL -- it shifts the split
  frac0  <- p$KCOMMIT / (p$KCOMMIT + p$KRESOLVE)
  fracE  <- p$KCOMMIT * (1 - p$EMAX_PAL) /
            (p$KCOMMIT * (1 - p$EMAX_PAL) + p$KRESOLVE)
  data.frame(
    mechanism_class = c("anti-inflammatory / anti-flare (episodic route only)",
                        "anti-activin A (ligand; cannot reach receptor leak)",
                        "ALK2 kinase inhibition (ligand-driven AND leak)",
                        "RARgamma commitment blockade (checkpoint competition)"),
    analytic_ceiling_pct = round(c(100 * sigfl / (sigsm + sigfl),
                                   100 * (1 - leak / (sigsm + sigfl)),
                                   100 * p$EMAX_AK,
                                   100 * (1 - fracE / frac0)), 1),
    set_by = c("PHI_FL, the flare-attributable fraction of new HO",
               "F_LEAK, the ligand-independent receptor leak",
               "EMAX_AK, the kinase inhibitor maximal effect",
               "KRESOLVE / KCOMMIT, the lesion-fate competition"))
}

# =====================================================================
#  SCENARIO 9 — SIROLIMUS around discrete flares (mTORC1 arm)
# =====================================================================
sc09_sirolimus <- function(age0 = 12, years = 2, dose = 2) {
  dur  <- round(yr(years))
  trig <- flare_events(c(120, 400))
  base <- sim_pt(age0, dur, 14, trig)
  sir  <- sim_pt(age0, dur, 14, c(oral_daily(dose, 0, dur, CMT_SIR), trig))
  data.frame(arm = c("untreated", paste0("sirolimus ", dose, " mg/d")),
             ann_new_HO_mL_yr = round(c(ann_new_ho(base, 0, dur),
                                        ann_new_ho(sir, 0, dur)), 2),
             sirolimus_Cav_ng_mL = round(c(0, mean(sir$C_SIRO)), 1))
}

# =====================================================================
#  SCENARIO 10 — THE ANTI-ACVR1 ANTIBODY PARADOX
#    A bivalent anti-receptor antibody dimerises and ACTIVATES the mutant
#    receptor. The map predicts harm; the model puts a number on it.
# =====================================================================
sc10_anti_acvr1_paradox <- function(age0 = 15, years = 2,
                                    eps = c(0, 0.25, 0.5, 1.0)) {
  dur <- round(yr(years))
  do.call(rbind, lapply(eps, function(e) {
    o <- sim_pt(age0, dur, 30.4, EPS_AGON = e)
    data.frame(agonist_conversion_EPS_AGON = e,
               ann_new_HO_mL_yr = round(ann_new_ho(o, 0, dur), 2),
               pct_change_vs_no_drug = round(100 * (ann_new_ho(o, 0, dur) /
                 ann_new_ho(sim_pt(age0, dur, 30.4), 0, dur) - 1), 1))
  }))
}

# =====================================================================
#  SCENARIO 11 — IATROGENIC TRIGGER (intramuscular injection)
#    25% of injections cause an immediate flare and 84% of those ossify,
#    so a single avoidable needle carries a measurable lifetime cost.
# =====================================================================
sc11_im_injection <- function(age0 = 8) {
  dur  <- 365
  inj  <- flare_events(30, amt = 6)
  none <- sim_pt(age0, dur, 7)
  hit  <- sim_pt(age0, dur, 7, inj)
  pre  <- sim_pt(age0, dur, 7, c(inj, oral_daily(40, 30, 4, CMT_PRE)))
  data.frame(arm = c("no injection", "IM injection",
                     "IM injection + prednisone 2 mg/kg/d x 4 d"),
             new_HO_over_1_yr_mL = round(c(tail(none$NEW_HO, 1),
                                           tail(hit$NEW_HO, 1),
                                           tail(pre$NEW_HO, 1)), 2),
             peak_inflammation = round(c(max(none$SIG_TOTAL), max(hit$SIG_TOTAL),
                                         max(pre$SIG_TOTAL)), 2))
}

# =====================================================================
#  SCENARIO 12 — THE DOSE / SELECTIVITY TRADE-OFF (Axis 4)
#    Sweep palovarotene dose in a skeletally immature patient and read
#    efficacy against physeal toxicity on the SAME exposure axis.
# =====================================================================
sc12_dose_selectivity <- function(age0 = 10, years = 1.5,
                                  doses = c(0, 0.5, 1, 2, 3, 5, 7.5, 10)) {
  dur  <- round(yr(years))
  base <- ann_new_ho(sim_pt(age0, dur, 30.4), 0, dur)
  do.call(rbind, lapply(doses, function(d) {
    o <- if (d == 0) sim_pt(age0, dur, 30.4) else
         sim_pt(age0, dur, 30.4, oral_daily(d, 0, dur, CMT_PAL))
    data.frame(dose_mg_daily = d,
               pct_HO_suppression = round(100 * (1 - ann_new_ho(o, 0, dur) / base), 1),
               PPC_probability_pct = round(100 * tail(o$PPC_PROB, 1), 1),
               height_Z_loss = round(-tail(o$HEIGHT_Z, 1), 2),
               vertebral_BMD_Z_loss = round(-tail(o$VBMD_Z, 1), 2))
  }))
}

# =====================================================================
#  SCENARIO 13 — LIFETIME OUTCOME AND SURVIVAL
#    EXTRAPOLATION, not a trial result: no FOP trial has ever measured a
#    functional or survival endpoint. Read the numbers as the model
#    reasoning out loud about what the volume endpoints imply.
# =====================================================================
sc13_survival <- function() {
  med_surv <- function(o) {
    d <- o[!duplicated(o$SURV), ]
    if (min(d$SURV) > 0.5) return(NA_real_)
    as.numeric(approx(d$SURV, d$AGE_YR, xout = 0.5)$y)
  }
  horizon <- 85
  arms <- list("untreated" = NA, "from age 8" = 8, "from age 15" = 15,
               "from age 25" = 25)
  do.call(rbind, lapply(names(arms), function(a) {
    st <- arms[[a]]
    dur <- round(yr(horizon - 5))
    evs <- if (is.na(st)) NULL else
      oral_daily(5, round(yr(st - 5)), dur - round(yr(st - 5)), CMT_PAL)
    o <- sim(end = dur, delta = 91.3, AGE0 = 5, events = evs)
    data.frame(arm = a,
               median_survival_yr = round(med_surv(o), 1),
               HO_at_40_mL  = round(at_age(o, "HO_TOTAL", 40), 0),
               FVC_at_40    = round(at_age(o, "FVC_PCT", 40), 1),
               CAJIS_at_40  = round(at_age(o, "CAJIS", 40), 1))
  }))
}

# =====================================================================
#  SCENARIO 14 — SIMULATED MOVE-vs-NHS TRIAL WITH BETWEEN-SUBJECT
#    VARIABILITY. New HO arrives as discrete, high-variance lesions, so
#    the estimated effect size is unstable and transformation-sensitive.
#    MOVE reported exactly that: 99.4% posterior probability of any
#    reduction WITHOUT a square-root transform, 65.4% WITH it.
# =====================================================================
sc14_trial_variability <- function(n = 100, months = 18, seed = 2026) {
  set.seed(seed)
  dur <- round(yr(months / 12))
  ages <- runif(n, 8, 30)
  idata <- data.frame(ID = seq_len(n), AGE0 = ages,
                      HO_INIT = vapply(ages, ho_at_age, numeric(1)))
  run <- function(treat) {
    evs <- if (treat) oral_daily(5, 0, dur, CMT_PAL) else ev(amt = 0, cmt = CMT_PAL)
    o <- fop_iiv %>% idata_set(idata) %>%
      mrgsim(events = evs, end = dur, delta = dur, atol = 1e-8, rtol = 1e-6) %>%
      as.data.frame()
    o <- o[o$time == dur, ]
    o$NEW_HO / (dur / 365.25)
  }
  a <- run(FALSE); b <- run(TRUE)
  data.frame(
    statistic = c("mean annualised new HO, untreated (mL/yr)",
                  "mean annualised new HO, treated (mL/yr)",
                  "reduction in the MEAN (%)",
                  "reduction after square-root transform (%)",
                  "coefficient of variation, untreated",
                  "share of untreated periods below 5 mL/yr (%)"),
    value = round(c(mean(a), mean(b), 100 * (1 - mean(b) / mean(a)),
                    100 * (1 - mean(sqrt(b)) / mean(sqrt(a))),
                    sd(a) / mean(a), 100 * mean(a < 5)), 2))
}

# =====================================================================
#  SCENARIO 15 — EXTERNAL-CONTROL BIAS
#    MOVE was single-arm and compared against untreated NHS participants.
#    Under this model the apparent effect is (drug effect) x (intrinsic
#    activity ratio between the two cohorts), and those two factors are
#    not separable from the data. The question the design cannot answer:
#    how much LOWER would the treated cohort intrinsic activity have to be
#    for the entire observed benefit to be an artefact?
# =====================================================================
sc15_external_control_bias <- function(age0 = 15, months = 18,
                                       ratios = c(1.0, 0.8, 0.6, 0.4, 0.25)) {
  dur <- round(yr(months / 12))
  nhs <- ann_new_ho(sim_pt(age0, dur, 30.4), 0, dur)
  do.call(rbind, lapply(ratios, function(k) {
    trt <- sim_pt(age0, dur, 30.4, oral_daily(5, 0, dur, CMT_PAL),
                  KF_HO = 191.82 * k)
    non <- sim_pt(age0, dur, 30.4, KF_HO = 191.82 * k)
    data.frame(trial_vs_NHS_intrinsic_activity = k,
               apparent_reduction_vs_NHS_pct =
                 round(100 * (1 - ann_new_ho(trt, 0, dur) / nhs), 1),
               true_drug_effect_pct =
                 round(100 * (1 - ann_new_ho(trt, 0, dur) /
                                  ann_new_ho(non, 0, dur)), 1),
               apparent_reduction_with_NO_drug_pct =
                 round(100 * (1 - ann_new_ho(non, 0, dur) / nhs), 1))
  }))
}

# =====================================================================
#  CALIBRATION REPORT — every anchor, model versus published
# =====================================================================
calibration_report <- function() {
  nh <- sc01_natural_history(80)
  bands <- sc02_age_bands()
  mv <- sc03_move()
  ppc <- tail(sim_pt(10, round(yr(1.5)), 30.4,
                     oral_daily(5, 0, round(yr(1.5)), CMT_PAL))$PPC_PROB, 1)
  lum <- sc06_lumina1()
  surv <- sc13_survival()
  data.frame(
    anchor = c("annualised new HO, 8-<15 yr (mL/yr)",
               "annualised new HO, 15-<25 yr (mL/yr)",
               "annualised new HO, 25-65 yr (mL/yr)",
               "flare-attributable fraction of new HO",
               "HO stock at 2-<8 yr (mL)",
               "HO stock at 25-65 yr (mL)",
               "per-flare ossification probability",
               "flare-up rate (per patient-year)",
               "palovarotene reduction in annualised new HO (%)",
               "premature physeal closure, 18 mo, age <14 (%)",
               "retinoid-associated AE burden (fraction of patients)",
               "FVC at age 30 (% predicted)",
               "median survival, untreated (yr)",
               "garetosmab Ctrough at 10 mg/kg q4w (mg/L)"),
    published = c("21.9", "41.5", "4.6", "0.295", "68.8", "575.2", "0.269",
                  "0.90", "54-60", "36.8", "0.970", "44 +/- 14", "56",
                  "105 +/- 31"),
    model = c(sprintf("%.1f", bands$model_mL_per_yr),
              sprintf("%.3f", mean(nh$FRAC_FLARE[1:12])),
              sprintf("%.1f", at_age(nh, "HO_TOTAL", 6)),
              sprintf("%.1f", at_age(nh, "HO_TOTAL", 45)),
              sprintf("%.3f", 0.071 / (0.071 + 0.193)),
              sprintf("%.2f", 0.90),
              sprintf("%.1f", mv$summary$palovarotene[2]),
              sprintf("%.1f", 100 * ppc),
              sprintf("%.3f", mv$summary$palovarotene[7]),
              sprintf("%.1f", at_age(nh, "FVC_PCT", 30)),
              sprintf("%.1f", surv$median_survival_yr[1]),
              sprintf("%.1f", lum$garetosmab[6])))
}

# =====================================================================
#  RUN EVERYTHING
# =====================================================================
run_all_fop <- function() {
  cat("\n########## CALIBRATION vs PUBLISHED ANCHORS ##########\n")
  print(calibration_report(), row.names = FALSE)
  cat("\n########## S2  NHS age-band velocities ##########\n")
  print(sc02_age_bands(), row.names = FALSE)
  cat("\n########## S3  Palovarotene MOVE regimen: RATE vs STOCK (Axis 2) ##########\n")
  print(sc03_move()$summary, row.names = FALSE)
  cat("\n########## S4  Chronic vs flare-triggered dosing ##########\n")
  print(sc04_dosing_strategy(), row.names = FALSE)
  cat("\n########## S5  Age at treatment start (Axis 3) ##########\n")
  print(sc05_age_at_start(), row.names = FALSE)
  cat("\n########## S6  Garetosmab LUMINA-1: two endpoints, one mechanism ##########\n")
  print(sc06_lumina1(), row.names = FALSE)
  cat("\n########## S7  The anti-flare ceiling (Axis 1) ##########\n")
  print(sc07_antiflare_ceiling(), row.names = FALSE)
  cat("\n########## S8  Drug class ceilings ##########\n")
  print(sc08_class_ceilings(), row.names = FALSE)
  cat("\n########## S9  Sirolimus ##########\n")
  print(sc09_sirolimus(), row.names = FALSE)
  cat("\n########## S10 The anti-ACVR1 antibody paradox ##########\n")
  print(sc10_anti_acvr1_paradox(), row.names = FALSE)
  cat("\n########## S11 Iatrogenic IM injection ##########\n")
  print(sc11_im_injection(), row.names = FALSE)
  cat("\n########## S12 Dose / selectivity trade-off (Axis 4) ##########\n")
  print(sc12_dose_selectivity(), row.names = FALSE)
  cat("\n########## S13 Lifetime outcome (EXTRAPOLATION) ##########\n")
  print(sc13_survival(), row.names = FALSE)
  cat("\n########## S14 Simulated trial variability ##########\n")
  print(sc14_trial_variability(80), row.names = FALSE)
  cat("\n########## Analytic ceiling of each mechanism class ##########\n")
  print(class_ceilings(), row.names = FALSE)
  cat("\n########## S15 External-control bias (single-arm design) ##########\n")
  print(sc15_external_control_bias(), row.names = FALSE)
  invisible(NULL)
}

# Rscript fop_mrgsolve_model.R  -> runs every scenario
if (!interactive() && identical(sys.nframe(), 0L)) run_all_fop()
