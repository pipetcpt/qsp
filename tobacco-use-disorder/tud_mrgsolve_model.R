## ===========================================================================
##  Tobacco Use Disorder (Nicotine Dependence) — QSP model for mrgsolve
##  ---------------------------------------------------------------------------
##  File     : tud_mrgsolve_model.R
##  Model    : 31 ODE compartments · 10 therapy scenarios
##  Requires : mrgsolve (>= 1.0), dplyr, tidyr, ggplot2
##
##  ---------------------------------------------------------------------------
##  CENTRAL IDEA OF THIS MODEL
##  ---------------------------------------------------------------------------
##  Nicotine dependence is NOT modelled here as "nicotine leaves and the
##  smoker feels bad".  It is modelled as a mismatch between TWO PROCESSES
##  WITH VERY DIFFERENT TIME CONSTANTS acting on the same receptor pool:
##
##     (a) alpha4beta2* nAChR DESENSITIZATION   —  t1/2 ~ 7 min on / ~2 h off
##     (b) alpha4beta2* nAChR UPREGULATION      —  t1/2 ~ 14 days
##
##  DOPAMINERGIC SUPPLY is fast: DA tracks the activated, non-desensitized
##  receptor fraction and re-equilibrates in minutes.  The HEDONIC SET-POINT
##  is slow: SETP is a ~14-day moving average of DA (Koob's opponent-process
##  allostasis), on the same time constant as the receptor pool RUP.
##
##     dSETP/dt = KSET*(DA - SETP)              (slow, t1/2 = 14 d)
##     DEFICIT  = max(0, SETP - DA)             (the withdrawal driver)
##
##  In a smoker at steady state DEFICIT ~ 0 BY CONSTRUCTION — the set-point
##  has adapted to whatever tone the smoking pattern delivers.  On quitting,
##  SUPPLY collapses within hours while the SET-POINT decays over weeks, so
##  the withdrawal syndrome is the TRANSIENT BETWEEN TWO MATCHED STATES, not
##  a state of its own.  Several clinical observations then fall out of the
##  model rather than being asserted:
##
##  WHAT THE MODEL ACTUALLY PRODUCES (computed, not asserted — see section 7)
##
##   1. WITHDRAWAL IS A TIME-CONSTANT PHENOMENON. With no symptom time course
##      supplied anywhere, an unaided quit gives MNWS peaking on day 5.2 and
##      staying above half its peak for 35 d, above a quarter for 54 d, while
##      beta2* upregulation decays from +70% with t50 = 15 d. HONEST CAVEAT:
##      the clinical course is faster (peak day 1-3, largely resolved by
##      2-4 weeks), so the model runs roughly 1.5-2x slow. The shape and the
##      mechanism are right; the constants KSET/KOUTA are not fully calibrated.
##
##   2. PARTIAL AGONISM IS NOT AN EFFICACY OPTIMUM — AND THE MODEL SAYS WHY.
##      Sweeping EMAXV from 0 to 1 (experiment 5c) raises CAR MONOTONICALLY
##      (32.8% -> 39.4%) while the lapse-blockade term stays pinned at 0.227.
##      The reason is arithmetic: blockade depends on OCCUPANCY, and occupancy
##      is independent of intrinsic activity, so a full agonist keeps all the
##      blockade and gains all the relief. The therapeutic argument for PARTIAL
##      agonism is therefore NOT superior efficacy at the receptor; it is that
##      partial intrinsic activity plus non-pulsatile oral delivery limits the
##      drug own reinforcing value. Consistent with this, adding a patch
##      (full-agonist activity) to varenicline lowers peak MNWS 1.43 -> 1.31.
##
##   3. NICOTINE TITRATION + A FIXED NRT DOSE GENERATES THE NMR RESULT.
##      One CYP2A6 scalar sets clearance; fast metabolizers therefore smoke
##      more (14 -> 27 cig/day), accumulate more reinforcement episodes and a
##      HIGHER set-point, yet a fixed 21 mg patch delivers them LESS nicotine
##      (replacement ratio 0.85 -> 0.42). Patch CAR falls 31.0% -> 23.3% across
##      the CYP2A6 range while varenicline (renally cleared) moves only
##      37.7% -> 34.2%, and that residual is the shared dependence gradient,
##      not a drug effect. No interaction term is fitted anywhere.
##
##   4. beta2* OCCUPANCY SATURATES, SO IT CANNOT DISCRIMINATE DELIVERY MODE.
##      With EC50 = 0.87 ng/mL, a patch and a cigarette both sit near 0.9
##      occupancy. The model therefore carries an explicit route-resolved
##      PHASIC term (FPHCIG 1.0 / FPHORAL 0.45 / FPHPATCH 0.05) scaled by
##      reinforcement-episode frequency. The same saturation predicts that
##      reduced-nicotine cigarettes do nothing until nicotine per cigarette
##      falls BELOW ~0.10 mg: CAR 12.5 -> 13.1 -> 14.3% at 1.1/0.55/0.25 mg,
##      then 17.2 -> 23.8 -> 28.2% at 0.10/0.03/0.015 mg (experiment 5e).
##
##  ---------------------------------------------------------------------------
##  SCOPE / ASSUMPTIONS — read before using any number from this file
##  ---------------------------------------------------------------------------
##  * PK/PD trajectories are simulated CONDITIONAL ON REMAINING ABSTINENT
##    after the target quit date (TQD). PABST is the probability of still
##    being continuously abstinent (a survival state), reported separately;
##    it does not feed back into the nicotine input. This is the standard
##    conditional formulation and it means the model UNDER-estimates relief
##    in relapsers (who re-dose themselves).
##  * Smoking input is a diurnal continuous rate over a 16 h waking day, not
##    20 discrete arterial boluses. Rate-of-rise reinforcement is therefore
##    NOT resolved kinetically; it is represented as an explicit per-product
##    attribute (FPHCIG / FPHORAL / FPHPATCH) feeding the phasic dopamine
##    term and the lapse hazard. This is a modelling choice, not a derivation.
##  * The BEHAVIOURAL parameters (HAZ0, B1QSU, B2WD, B3HAB, GBLOCK, GBLOCKV,
##    GREINF) are FITTED to the 12-week continuous-abstinence rates of 7 trial
##    arms and are marked [FITTED] in $PARAM. Everything upstream of the
##    hazard — PK, metabolism, occupancy, desensitization, upregulation,
##    dopamine, set-point — is fixed from mechanism and never fitted to CAR.
##  * Reported CARs are POPULATION AVERAGES over a 5-point CYP2A6 activity
##    distribution (0.35/0.65/1.0/1.5/2.0 weighted 0.10/0.20/0.40/0.20/0.10),
##    because a fixed-dose NRT arm is heterogeneous by construction.
##  * Effective in vivo Kd values (KDNIC, KDVAR, KDCYT) are ~10x the in vitro
##    Ki because they absorb brain partitioning, non-specific binding and the
##    plasma-to-biophase step. They are calibrated to PET occupancy, NOT to
##    binding assays — see the calibration table below.
##  * Educational / research model. Not validated for clinical or regulatory
##    use. Do not use to choose therapy for a patient.
##
##  ---------------------------------------------------------------------------
##  PARAMETER CALIBRATION TABLE  (target -> source)
##  ---------------------------------------------------------------------------
##  Nicotine CL 72 L/h, Vss ~182 L, t1/2 ~1.8 h      Benowitz 2009 Handb Exp Pharmacol
##  20 cig/day -> plasma nicotine ~12-17 ng/mL       Benowitz 1997 / Russell 1980
##  Patch 21 mg/24 h -> Css ~12-17 ng/mL             Fant 2000 / Benowitz 2008
##  Cotinine t1/2 16 h; 3HC/COT (NMR) ~0.20-0.45     Dempsey 2004; Lerman 2015
##  beta2* occupancy EC50 = 0.87 ng/mL nicotine      Brody 2006 Arch Gen Psychiatry
##  1 cigarette -> 88% beta2* occupancy               Brody 2006
##  Varenicline 1 mg BID -> ~90-95% beta2* occupancy  Lotfipour 2012 / Rollema 2007
##  Varenicline Css ~6-9 ng/mL, t1/2 24 h, Vd 415 L   Faessel 2006 (Chantix label)
##  beta2* upregulation normalizes in 3-4 weeks       Cosgrove 2009 Arch Gen Psychiatry
##  MNWS peaks day 1-3, resolves 2-4 weeks            Hughes 2007 Nicotine Tob Res
##  CAR wk 9-12: VAR 33.5 / BUP 22.6 / PATCH 23.4 /
##               PBO 12.5 (%)                         EAGLES, Anthenelli 2016 Lancet
##  Cytisinicline 3 mg TID CAR ~32% vs 7% pbo         ORCA-2, Rigotti 2023 JAMA
##  Combination NRT > patch alone                     Theodoulou 2023 Cochrane
##  Varenicline + patch: MIXED evidence               Koegelenberg 2014 JAMA (+)
##    -- the model PREDICTS non-additivity; see the note in section 5f.
##  Post-cessation weight gain ~4-5 kg at 12 mo       Aubin 2012 BMJ
##  FEV1 decline 60 vs 30 mL/yr, slope recovery       Fletcher & Peto 1977; Anthonisen 1994
## ===========================================================================

library(mrgsolve)
suppressMessages({
  library(dplyr)
  library(tidyr)
})

# ---------------------------------------------------------------------------
# 1) MODEL CODE
# ---------------------------------------------------------------------------
tud_code <- '
$PROB
# Tobacco Use Disorder QSP model
- 31 ODE compartments
- nicotine PK + CYP2A6 pharmacogenetics + a4b2* receptor dynamics
- mesolimbic dopamine, allostatic withdrawal, craving, lapse hazard
- varenicline / cytisinicline / bupropion / NRT pharmacology

$PARAM @annotated
// ---------------- study design ------------------------------------------
TQD      : 2160 : Target quit date, h from t=0 (90 d smoking run-in)
TSTOPRX  : 4176 : Medication stop time, h (default TQD + 12 wk)
TSTARTRX : 1992 : Medication start time, h (default TQD - 7 d)
COUNSEL  : 1    : Behavioural support flag (0/1)

// ---------------- baseline smoker phenotype ------------------------------
CPD0     : 20   : Baseline cigarettes per day
NICCIG   : 1.1  : Nicotine absorbed per cigarette, mg
YRSMOKE  : 25   : Years of smoking (sets habit strength)
TITEXP   : 0.5  : Nicotine-titration exponent (intake vs clearance)
CLREF    : 72   : Reference nicotine clearance for titration, L/h
CPDREF   : 20   : Reference cigarettes/day for episode-frequency scaling
EPEXP    : 0.7  : Exponent, reinforcement episodes/day -> phasic gain
F2A6     : 1.0  : CYP2A6 activity multiplier (0.35 slow / 1.0 normal / 1.7 fast)
A5RISK   : 0    : CHRNA5 rs16969968 risk allele count (0/1/2)
PSYHX    : 0    : Psychiatric history flag (0/1)

// ---------------- nicotine PK (nmol, L, h) -------------------------------
CLNP     : 18   : Non-CYP2A6 nicotine clearance, L/h
CLN2A6   : 54   : CYP2A6-dependent nicotine clearance, L/h
V2N      : 60   : Nicotine central volume, L
V3N      : 122  : Nicotine peripheral volume, L
QN       : 140  : Nicotine intercompartmental clearance, L/h
KEON     : 13.9 : Nicotine brain biophase rate constant, 1/h
KASKIN   : 0.06 : Transdermal absorption rate constant, 1/h (t1/2 ~12 h)
KAMOU    : 1.4  : Buccal absorption rate constant, 1/h
FCOT     : 0.70 : Fraction of nicotine clearance to cotinine
CLCOT    : 3.0  : Cotinine clearance, L/h
VCOT     : 70   : Cotinine volume, L
FHC      : 0.45 : Fraction of cotinine clearance to 3HC (CYP2A6)
CLHC     : 3.5  : 3-hydroxycotinine clearance, L/h
VHC      : 45   : 3-hydroxycotinine volume, L

// ---------------- varenicline PK ----------------------------------------
KAV      : 0.5  : Varenicline absorption rate constant, 1/h
CLV      : 12.0 : Varenicline CL/F, L/h
VV       : 415  : Varenicline V/F, L
KEOV     : 0.9  : Varenicline biophase rate constant, 1/h
RFV      : 1.0  : Renal function multiplier on varenicline CL (CrCl<30 -> 0.5)

// ---------------- cytisinicline PK --------------------------------------
KAC      : 1.2  : Cytisinicline absorption rate constant, 1/h
CLC      : 9.1  : Cytisinicline CL/F, L/h
VC       : 63   : Cytisinicline V/F, L

// ---------------- bupropion PK ------------------------------------------
KAB      : 0.5  : Bupropion SR absorption rate constant, 1/h
CLB      : 110  : Bupropion CL/F, L/h
VB       : 1500 : Bupropion V/F, L
FOH      : 0.50 : Fraction of bupropion CL to hydroxybupropion
CLOH     : 8.0  : Hydroxybupropion clearance, L/h
VOH      : 200  : Hydroxybupropion volume, L
F2B6     : 1.0  : CYP2B6 activity multiplier (0.45 for *6/*6)

// ---------------- a4b2* receptor pharmacology (nM) ----------------------
KDNIC    : 5.4  : Effective in vivo Kd, nicotine at a4b2* (Brody 2006)
KDVAR    : 3.6  : Effective in vivo Kd, varenicline at a4b2*
KDCYT    : 35   : Effective in vivo Kd, cytisinicline (low brain penetration)
EMAXN    : 1.00 : Intrinsic activity, nicotine (reference full agonist)
EMAXV    : 0.45 : Intrinsic activity, varenicline (partial agonist)
EMAXC    : 0.35 : Intrinsic activity, cytisinicline (partial agonist)
KINACHR  : 8000 : Hydroxybupropion non-competitive nAChR block Ki, nM
KIDAT    : 2500 : Hydroxybupropion DAT/NET inhibition IC50, nM
GDAT     : 1.30 : Gain of DAT/NET inhibition on dopamine tone
FPATCH   : 0.68 : Transdermal bioavailability of the nominal patch dose
FPHCIG   : 1.00 : Phasic (rate-of-rise) reinforcement fraction, cigarette
FPHORAL  : 0.45 : Phasic reinforcement fraction, oral PRN NRT
FPHPATCH : 0.05 : Phasic reinforcement fraction, transdermal patch
EMAXPH   : 0.85 : Maximal phasic contribution to dopamine tone
PHID     : 0.78 : Fraction of signalling lost per unit desensitization
KOND     : 6.0  : Desensitization on-rate, 1/h
KOFFD    : 0.35 : Resensitization rate, 1/h
KOUTR    : 0.00206 : Receptor pool turnover, 1/h (t1/2 = 14 d)
EMAXR    : 0.75 : Maximal desensitization-driven upregulation

// ---------------- dopamine / allostasis / symptoms -----------------------
KDA      : 6.0  : Dopamine tone turnover, 1/h
EMAXDA   : 0.55 : Maximal tonic dopamine increase from nAChR activation
EC50DA   : 0.32 : Activation producing half-maximal dopamine response
HDA      : 1.6  : Hill coefficient, activation -> dopamine
KSET     : 0.00206 : Hedonic set-point adaptation rate, 1/h (t1/2 = 14 d)
KDEF     : 0.55 : Half-saturation of the normalized dopaminergic deficit
KINA     : 0.014: Allostatic load accrual rate, 1/h
KOUTA    : 0.0096 : Allostatic load resolution rate, 1/h (t1/2 = 3 d)
KWD      : 0.058: Withdrawal score turnover, 1/h (t1/2 = 12 h)
WDMAX    : 3.2  : Maximal MNWS composite score (0-4 scale)
WWDEF    : 0.55 : Weight of instantaneous deficit in withdrawal
WWALLO   : 0.45 : Weight of allostatic load in withdrawal
KQ       : 0.30 : Craving score turnover, 1/h
QSUMAX   : 6.5  : Maximal QSU-brief craving score (1-7 scale)
WQWD     : 0.45 : Weight of withdrawal on craving
WQDEF    : 0.30 : Weight of deficit on craving
WQCUE    : 0.25 : Weight of cue/habit on craving
GREL     : 0.22 : Craving relief per unit total a4b2* occupancy
GRELV    : 0.18 : Additional craving relief specific to partial agonists
GRELPH   : 0.30 : Craving relief from phasic (pulsatile) nicotine delivery
CUEAMP   : 1.0  : Waking-hours cue pressure amplitude
KHOFF    : 0.00048 : Habit extinction rate, 1/h (t1/2 = 60 d)
GREINF   : 1.377: Phasic-reinforcement gain on lapse hazard [FITTED]

// ---------------- lapse hazard / behaviour ------------------------------
HAZ0     : 0.000422 : Baseline lapse hazard, 1/h [FITTED to placebo CAR 12.5%]
B1QSU    : 0.089: Log-hazard coefficient on craving above threshold [FITTED]
B2WD     : 0.783: Log-hazard coefficient on withdrawal [FITTED]
B3HAB    : 0.604: Log-hazard coefficient on habit strength [FITTED]
BCOUNS   : -0.22: Log-hazard effect of behavioural support
BPSY     : 0.18 : Log-hazard effect of psychiatric history
BA5      : 0.12 : Log-hazard effect per CHRNA5 risk allele
QSUTHR   : 2.0  : Craving score below which hazard is not amplified
GBLOCK   : 0.540: Lapse-reward blockade per unit total occupancy [FITTED]
GBLOCKV  : 0.326: Additional blockade specific to partial agonists [FITTED]

// ---------------- downstream / safety -----------------------------------
KINCO    : 0.055: Carboxyhaemoglobin formation per cig/day, %/h
KOUTCO   : 0.173: COHb elimination, 1/h (t1/2 = 4 h)
KWT      : 0.00032 : Weight-gain rate constant, 1/h (t1/2 = 90 d)
WTSS     : 5.0  : Steady-state post-cessation weight gain, kg
FEVSMOKE : 60   : FEV1 decline while smoking, mL/yr
FEVQUIT  : 30   : FEV1 decline after quitting, mL/yr
FEV1BASE : 3000 : Baseline FEV1, mL
KNAUS    : 0.20 : Nausea index turnover, 1/h
GNAUSV   : 0.85 : Nausea gain, varenicline occupancy
GNAUSC   : 0.55 : Nausea gain, cytisinicline occupancy
GNAUSN   : 0.30 : Nausea gain, exogenous nicotine
KTOL     : 0.0072 : Nausea tolerance development, 1/h (t1/2 = 4 d)
KSLP     : 0.12 : Sleep disturbance turnover, 1/h
GSLPV    : 0.45 : Sleep disturbance gain, varenicline
KDISC    : 0.00018 : Discontinuation rate per unit AE burden, 1/h

$CMT @annotated
NICC   : Nicotine central, nmol
NICP   : Nicotine peripheral, nmol
NICB   : Nicotine brain biophase, nmol-equivalent
ASKIN  : Transdermal nicotine depot, nmol
AMOUTH : Buccal nicotine depot, nmol
COT    : Cotinine, nmol
HCOT   : 3-hydroxycotinine, nmol
AGV    : Varenicline gut, nmol
VARC   : Varenicline central, nmol
VARB   : Varenicline biophase, nmol-equivalent
AGC    : Cytisinicline gut, nmol
CYTC   : Cytisinicline central, nmol
AGB    : Bupropion gut, nmol
BUPC   : Bupropion central, nmol
OHB    : Hydroxybupropion, nmol
DES    : Desensitized a4b2* fraction (0-1)
RUP    : a4b2* receptor pool, relative to never-smoker
DA     : NAc dopamine tone, relative to never-smoker
ALLO   : Allostatic load (0-1)
WD     : MNWS withdrawal composite (0-4)
QSU    : QSU-brief craving (1-7)
HABIT  : Cue/habit strength (0-1)
PABST  : Probability of continuous abstinence (0-1)
COHB   : Carboxyhaemoglobin, %
WT     : Weight change from baseline, kg
FEV    : FEV1, mL
NAUS   : Nausea index (0-1)
TOL    : Nausea tolerance (0-1)
SLP    : Sleep disturbance index (0-1)
PADH   : Probability still on medication (0-1)
SETP   : Hedonic set-point (14-day moving average of dopamine tone)

$GLOBAL
#define CNIC   (NICC/V2N)                 // plasma nicotine, nM
#define CNB    (NICB/V2N)                 // brain biophase nicotine, nM
#define CCOT   (COT/VCOT)                 // cotinine, nM
#define CHC    (HCOT/VHC)                 // 3HC, nM
#define CVAR   (VARC/VV)                  // varenicline plasma, nM
#define CVB    (VARB/VV)                  // varenicline biophase, nM
#define CCYT   (CYTC/VC)                  // cytisinicline, nM
#define CBUP   (BUPC/VB)                  // bupropion, nM
#define COHBUP (OHB/VOH)                  // hydroxybupropion, nM

// unit helpers -----------------------------------------------------------
#define NMOL_PER_MG_NIC 6165.0            // 1 mg nicotine  = 6165 nmol
#define NGML_PER_NM_NIC 0.16223           // nM -> ng/mL (MW 162.23)
#define NGML_PER_NM_COT 0.17622
#define NGML_PER_NM_HC  0.19222
#define NGML_PER_NM_VAR 0.21128
#define NGML_PER_NM_CYT 0.19025
#define NGML_PER_NM_BUP 0.23972
#define NGML_PER_NM_OHB 0.25572

$PREAMBLE
// nothing to precompute

$MAIN
// ---- initial conditions: chronic smoker at quasi-steady state ----------
// Receptor pool and habit start already adapted, so that the model does not
// spend the run-in period climbing out of a never-smoker state.
// The default 90-day smoking run-in lets DES / RUP / DA / SETP settle onto
// the true smoker steady state, so these initials only have to be close.
double DES_SS  = 0.62;
double RUP_SS  = 1.0 + EMAXR*DES_SS;
double CL_SS   = CLNP + CLN2A6*F2A6;
double CPD_SS  = CPD0*pow(CL_SS/CLREF, TITEXP);      // nicotine titration

RUP_0   = RUP_SS;
DES_0   = DES_SS;
HABIT_0 = YRSMOKE/(YRSMOKE + 8.0);                   // 25 yr -> 0.76
DA_0    = 1.70;                                      // settles during run-in
SETP_0  = 1.70;                                      // matched: DEFICIT ~ 0
QSU_0   = 1.4;
WD_0    = 0.15;
PABST_0 = 1.0;
PADH_0  = 1.0;
FEV_0   = FEV1BASE;
COHB_0  = KINCO*CPD_SS/KOUTCO;
WT_0    = 0.0;
ALLO_0  = 0.0;
NAUS_0  = 0.0;
TOL_0   = 0.0;
SLP_0   = 0.0;

$ODE
// ======================= 0) time-varying context ========================
double HOD    = fmod(SOLVERTIME, 24.0);
double AWAKE  = (HOD >= 7.0 && HOD < 23.0) ? 1.0 : 0.0;
double SMOKE  = (SOLVERTIME < TQD) ? 1.0 : 0.0;
double ONRX   = (SOLVERTIME >= TSTARTRX && SOLVERTIME < TSTOPRX) ? 1.0 : 0.0;

// ======================= 1) nicotine PK =================================
double CLNIC  = CLNP + CLN2A6*F2A6;                       // L/h
// NICOTINE TITRATION: to defend a brain nicotine level, intake must scale
// with clearance. Partial compensation (exponent TITEXP ~0.5) reproduces the
// observed NMR-vs-cigarettes/day relationship: F2A6 0.35 -> 14 cig/day,
// F2A6 1.0 -> 20, F2A6 2.0 -> 27. This is what makes a FIXED patch dose
// under-replace fast metabolizers (see experiment 5b).
double CPDT   = CPD0*pow(CLNIC/CLREF, TITEXP);
// smoking input spread over the 16 h waking day (nmol/h)
double RSMOKE = SMOKE*AWAKE*CPDT*NICCIG*NMOL_PER_MG_NIC/16.0;

double kelN   = CLNIC/V2N;
double k23    = QN/V2N;
double k32    = QN/V3N;

dxdt_NICC   = RSMOKE + KASKIN*ASKIN + KAMOU*AMOUTH
              - kelN*NICC - k23*NICC + k32*NICP;
dxdt_NICP   = k23*NICC - k32*NICP;
dxdt_NICB   = KEON*(NICC - NICB);
dxdt_ASKIN  = -KASKIN*ASKIN;
dxdt_AMOUTH = -KAMOU*AMOUTH;

// metabolites (molar; CYP2A6 scales BOTH oxidation steps -> sets NMR)
double RCOT = FCOT*(CLN2A6*F2A6/CLNIC)*kelN*NICC;         // nicotine -> cotinine
double RHC  = FHC*F2A6*(CLCOT/VCOT)*COT;                  // cotinine -> 3HC
dxdt_COT    = RCOT - (CLCOT/VCOT)*COT;
dxdt_HCOT   = RHC  - (CLHC/VHC)*HCOT;

// ======================= 2) drug PK =====================================
dxdt_AGV  = -KAV*AGV;
dxdt_VARC =  KAV*AGV - (CLV*RFV/VV)*VARC;
dxdt_VARB =  KEOV*(VARC - VARB);

dxdt_AGC  = -KAC*AGC;
dxdt_CYTC =  KAC*AGC - (CLC/VC)*CYTC;

dxdt_AGB  = -KAB*AGB;
dxdt_BUPC =  KAB*AGB - (CLB/VB)*BUPC;
dxdt_OHB  =  FOH*F2B6*(CLB/VB)*BUPC - (CLOH/VOH)*OHB;

// ======================= 3) a4b2* competitive occupancy =================
double BN   = CNB/KDNIC;
double BV   = CVB/KDVAR;
double BC   = CCYT/KDCYT;
double DEN  = 1.0 + BN + BV + BC;
double FN   = BN/DEN;                                     // fractional occ, nicotine
double FV   = BV/DEN;                                     // varenicline
double FC   = BC/DEN;                                     // cytisinicline
double OCCT = FN + FV + FC;                               // total (PET-comparable)

// non-competitive block by hydroxybupropion
double BLK  = 1.0/(1.0 + COHBUP/KINACHR);

// numeric guards: the stiff solver can overshoot into small negatives, and
// pow(negative, 1.6) is NaN. Clamp before any fractional power is taken.
double DESc = DES; if(DESc < 0.0) DESc = 0.0; if(DESc > 1.0) DESc = 1.0;
double RUPc = RUP; if(RUPc < 0.0) RUPc = 0.0;

// intrinsic-activity-weighted activation of the AVAILABLE receptor pool.
// PHID < 1 because desensitization is neither complete nor uniform across
// cell types (alpha6beta2beta3 on DA terminals desensitizes far less than
// alpha4beta2 on VTA somata) — Nashmi 2007.
double IA    = EMAXN*FN + EMAXV*FV + EMAXC*FC;
double AVAIL = 1.0 - PHID*DESc;
double ACT   = IA*AVAIL*RUPc*BLK; if(ACT < 0.0) ACT = 0.0;

// ROUTE-RESOLVED PHASIC REINFORCEMENT.
// beta2* occupancy saturates (EC50 0.87 ng/mL), so occupancy cannot tell a
// cigarette from a patch. What distinguishes them is the RATE of rise. The
// current nicotine input is decomposed by route and weighted by that route
// phasic fraction; PHAS then feeds both dopamine and the lapse hazard.
double RSM  = RSMOKE;                                     // cigarettes
double RPA  = KASKIN*ASKIN;                               // transdermal
double RBU  = KAMOU*AMOUTH;                               // buccal/oral
double RTOT = RSM + RPA + RBU + 1e-9;
// Each cigarette is a discrete reinforcement EPISODE. Smoothing the smoking
// input over the waking day throws episode frequency away, so restore it
// explicitly (sub-linearly). This is what makes a fast CYP2A6 metabolizer
// smoking 27/day MORE dependent than a slow one on 14/day, even though the
// fast metabolizer plasma nicotine is the LOWER of the two.
double EPFREQ = pow(CPDT/CPDREF, EPEXP);
double FPH  = (FPHCIG*EPFREQ*RSM + FPHORAL*RBU + FPHPATCH*RPA)/RTOT;
double PHAS = FPH*FN;
double RISE = GREINF*PHAS;

// ======================= 4) receptor dynamics ===========================
dxdt_DES = KOND*OCCT*(1.0 - DESc) - KOFFD*DESc;
dxdt_RUP = KOUTR*(1.0 + EMAXR*DESc) - KOUTR*RUPc;

// ======================= 5) dopamine, need, deficit =====================
double IDAT  = COHBUP/(COHBUP + KIDAT);                   // DAT/NET inhibition
double DRIVE = 1.0 + EMAXDA*pow(ACT,HDA)/(pow(EC50DA,HDA) + pow(ACT,HDA))
                   + EMAXPH*PHAS;
dxdt_DA = KDA*(DRIVE*(1.0 + GDAT*IDAT) - DA);

// ALLOSTATIC SET-POINT: a ~14-day moving average of dopaminergic tone.
// This is the single structural commitment of the model (Koob & Le Moal 2001).
dxdt_SETP = KSET*(DA - SETP);
double DEF  = SETP - DA;  if(DEF < 0.0) DEF = 0.0;
double DEFn = DEF/(KDEF + DEF);                           // 0-1 saturating

// ======================= 6) allostasis, withdrawal, craving =============
dxdt_ALLO = KINA*DEFn*(1.0 - ALLO) - KOUTA*ALLO;

double ALLOc = ALLO; if(ALLOc < 0.0) ALLOc = 0.0;
double WDtarget = WDMAX*(WWDEF*DEFn + WWALLO*ALLOc);
dxdt_WD = KWD*(WDtarget - WD);

double CUE    = CUEAMP*AWAKE*HABIT;
// Craving relief has three additive sources: any occupant of the receptor,
// an extra term specific to slowly-dissociating partial agonists, and
// pulsatile delivery (why PRN gum relieves peak urges that a patch cannot).
double RELIEF = 1.0 - GREL*OCCT - GRELV*(FV + FC) - GRELPH*PHAS;
if(RELIEF < 0.05) RELIEF = 0.05;
double QSUt   = 1.0 + (QSUMAX - 1.0)*RELIEF*
                (WQWD*(WD/WDMAX) + WQDEF*DEFn + WQCUE*CUE);
dxdt_QSU = KQ*(QSUt - QSU);

dxdt_HABIT = -KHOFF*(1.0 - SMOKE)*HABIT + KHOFF*SMOKE*(1.0 - HABIT);

// ======================= 7) lapse hazard / abstinence ===================
double QEX  = QSU - QSUTHR;  if(QEX < 0.0) QEX = 0.0;
double LHAZ = HAZ0*exp(B1QSU*QEX + B2WD*WD + B3HAB*HABIT
                       + BCOUNS*COUNSEL + BPSY*PSYHX + BA5*A5RISK);
// A lapse only becomes a relapse if the lapse cigarette is still rewarding.
// Any receptor occupant blunts it (GBLOCK); partial agonists blunt it further
// because they dissociate slowly and are not displaced by an arterial spike
// as readily as steady-state nicotine is (GBLOCKV).
double LAPREW = 1.0 - GBLOCK*OCCT - GBLOCKV*(FV + FC);
if(LAPREW < 0.05) LAPREW = 0.05;
double HAZ = (SOLVERTIME < TQD) ? 0.0 : LHAZ*(0.45 + 0.55*LAPREW)*(1.0 + RISE);
dxdt_PABST = -HAZ*PABST;

// ======================= 8) biomarkers & organ systems ==================
dxdt_COHB = KINCO*CPDT*SMOKE - KOUTCO*COHB;

double ABSTF = (SOLVERTIME < TQD) ? 0.0 : PABST;          // abstinent fraction
dxdt_WT   = KWT*(WTSS*ABSTF - WT);

double SMOKEF = (SOLVERTIME < TQD) ? 1.0 : (1.0 - PABST);
double dFEV   = -(FEVSMOKE*SMOKEF + FEVQUIT*(1.0 - SMOKEF))/8760.0;
dxdt_FEV = dFEV;

// ======================= 9) tolerability ================================
double NAUSt = (GNAUSV*FV + GNAUSC*FC + GNAUSN*FN*(1.0 - SMOKE))*(1.0 - TOL);
dxdt_NAUS = KNAUS*(NAUSt - NAUS);
dxdt_TOL  = KTOL*(1.0 - TOL)*(FV + FC + FN) - 0.25*KTOL*TOL;

double SLPt = 0.55*(WD/WDMAX) + GSLPV*FV*ONRX;
dxdt_SLP  = KSLP*(SLPt - SLP);

double AEB = 0.6*NAUS + 0.4*SLP;
dxdt_PADH = -KDISC*AEB*PADH*ONRX;

$TABLE
double HOD_    = fmod(TIME, 24.0);
double SMOKE_  = (TIME < TQD) ? 1.0 : 0.0;
double CLNIC_  = CLNP + CLN2A6*F2A6;
double CPDT_   = CPD0*pow(CLNIC_/CLREF, TITEXP);

// concentrations in clinical units
double NIC_NGML = (NICC/V2N)*NGML_PER_NM_NIC;
double COT_NGML = (COT/VCOT)*NGML_PER_NM_COT;
double HC_NGML  = (HCOT/VHC)*NGML_PER_NM_HC;
double VAR_NGML = (VARC/VV)*NGML_PER_NM_VAR;
double CYT_NGML = (CYTC/VC)*NGML_PER_NM_CYT;
double BUP_NGML = (BUPC/VB)*NGML_PER_NM_BUP;
double OHB_NGML = (OHB/VOH)*NGML_PER_NM_OHB;

// nicotine metabolite ratio (the clinical biomarker)
double NMR_     = (COT > 1e-9) ? ((HCOT/VHC)/(COT/VCOT)) : 0.0;

// occupancy re-derived for reporting
double BN_ = (NICB/V2N)/KDNIC;
double BV_ = (VARB/VV)/KDVAR;
double BC_ = (CYTC/VC)/KDCYT;
double DEN_= 1.0 + BN_ + BV_ + BC_;
double OCC_NIC = BN_/DEN_;
double OCC_VAR = BV_/DEN_;
double OCC_CYT = BC_/DEN_;
double OCC_TOT = OCC_NIC + OCC_VAR + OCC_CYT;

double BLK_ = 1.0/(1.0 + (OHB/VOH)/KINACHR);
double IA_  = EMAXN*OCC_NIC + EMAXV*OCC_VAR + EMAXC*OCC_CYT;
double ACT_ = IA_*(1.0 - PHID*DES)*RUP*BLK_; if(ACT_ < 0.0) ACT_ = 0.0;

double DEF_  = SETP - DA; if(DEF_ < 0.0) DEF_ = 0.0;
double REPL_ = 0.0;  // nicotine replacement ratio vs the smoker set-point

// clinical readouts
double CO_PPM   = COHB*4.0;                    // exhaled CO ~ 4 x COHb%
double CPD_OBS  = CPDT_*(SMOKE_ + (1.0 - SMOKE_)*(1.0 - PABST));
double CAR_PCT  = 100.0*PABST;                 // continuous abstinence, %
double FEV_PCT  = 100.0*FEV/FEV1BASE;
double ADH_PCT  = 100.0*PADH;

$CAPTURE @annotated
NIC_NGML : Plasma nicotine, ng/mL
COT_NGML : Cotinine, ng/mL
HC_NGML  : 3-hydroxycotinine, ng/mL
NMR_     : Nicotine metabolite ratio (3HC/cotinine)
VAR_NGML : Varenicline, ng/mL
CYT_NGML : Cytisinicline, ng/mL
BUP_NGML : Bupropion, ng/mL
OHB_NGML : Hydroxybupropion, ng/mL
OCC_NIC  : Fractional a4b2* occupancy by nicotine
OCC_VAR  : Fractional a4b2* occupancy by varenicline
OCC_CYT  : Fractional a4b2* occupancy by cytisinicline
OCC_TOT  : Total a4b2* occupancy (PET-comparable)
ACT_     : Effective a4b2* activation signal
DEF_     : Dopaminergic deficit (SETP - DA)
CPDT_    : Titrated baseline cigarettes per day
CO_PPM   : Exhaled carbon monoxide, ppm
CPD_OBS  : Observed cigarettes per day
CAR_PCT  : Continuous abstinence probability, %
FEV_PCT  : FEV1, % of baseline
ADH_PCT  : Medication persistence, %
CLNIC_   : Nicotine clearance, L/h
'

mod <- mcode("tud", tud_code)

# ---------------------------------------------------------------------------
# 2) DOSING BUILDERS
#    All amounts are converted to nmol so the whole model is molar-consistent.
# ---------------------------------------------------------------------------
NMOL_MG_NIC <- 6165.0            # nicotine   MW 162.23
NMOL_MG_VAR <- 4733.0            # varenicline MW 211.28  (1 mg = 4733 nmol)
NMOL_MG_CYT <- 5256.0            # cytisinicline MW 190.25
NMOL_MG_BUP <- 4172.0            # bupropion HCl free base MW 239.7

TQD_H   <- 2160                  # quit at day 90 (after the smoking run-in)
FPATCH  <- 0.68                  # transdermal bioavailability of nominal dose
WK12    <- TQD_H + 12*7*24       # end of 12-week treatment
TEND    <- TQD_H + 52*7*24       # 1-year follow-up

# nicotine patch: 24 h "infusion" into the skin depot, given as a daily bolus
patch_ev <- function(mg, start, weeks, taper = TRUE) {
  n  <- weeks*7
  df <- data.frame(
    ID   = 1,
    time = start + (0:(n - 1))*24,
    cmt  = 4,                                  # ASKIN
    amt  = mg*NMOL_MG_NIC*FPATCH,
    evid = 1
  )
  if (taper && weeks >= 10) {                  # 21 -> 14 -> 7 mg standard taper
    df$amt[df$time >= start + 6*7*24] <- 14*NMOL_MG_NIC*FPATCH
    df$amt[df$time >= start + 8*7*24] <-  7*NMOL_MG_NIC*FPATCH
  }
  df
}

# PRN oral NRT (gum/lozenge) — 8 pieces/day during waking hours
prn_nrt_ev <- function(mg, start, weeks, per_day = 8) {
  days  <- 0:(weeks*7 - 1)
  hours <- seq(8, 22, length.out = per_day)
  tt    <- as.vector(outer(days*24, hours, "+")) + start
  data.frame(ID = 1, time = sort(tt), cmt = 5,     # AMOUTH
             amt = mg*NMOL_MG_NIC*0.55,            # ~55% swallowed/lost
             evid = 1)
}

# varenicline with the standard 1-week titration
var_ev <- function(start, weeks = 12) {
  sched <- function(d) {
    if (d < 3)      list(mg = 0.5, times = 9)                # d1-3  0.5 mg QD
    else if (d < 7) list(mg = 0.5, times = c(9, 21))         # d4-7  0.5 mg BID
    else            list(mg = 1.0, times = c(9, 21))         # d8+   1 mg BID
  }
  out <- lapply(0:(weeks*7 - 1), function(d) {
    s <- sched(d)
    data.frame(ID = 1, time = start + d*24 + s$times, cmt = 8,
               amt = s$mg*NMOL_MG_VAR, evid = 1)
  })
  do.call(rbind, out)
}

cyt_ev <- function(start, weeks = 12, mg = 3) {
  days  <- 0:(weeks*7 - 1)
  tt    <- as.vector(outer(days*24, c(8, 14, 20), "+")) + start
  data.frame(ID = 1, time = sort(tt), cmt = 11, amt = mg*NMOL_MG_CYT, evid = 1)
}

bup_ev <- function(start, weeks = 12) {
  out <- lapply(0:(weeks*7 - 1), function(d) {
    tt <- if (d < 3) start + d*24 + 8 else start + d*24 + c(8, 20)   # 150 QD -> BID
    data.frame(ID = 1, time = tt, cmt = 13, amt = 150*NMOL_MG_BUP, evid = 1)
  })
  do.call(rbind, out)
}

# ---------------------------------------------------------------------------
# 3) TEN THERAPY SCENARIOS
# ---------------------------------------------------------------------------
scenarios <- list(

  `1. Continued smoking (no quit attempt)` = list(
    par = list(TQD = 1e6, TSTARTRX = 1e6, TSTOPRX = 1e6, COUNSEL = 0),
    ev  = NULL),

  `2. Unaided abrupt quit (placebo + brief advice)` = list(
    par = list(TQD = TQD_H, TSTARTRX = 1e6, TSTOPRX = 1e6, COUNSEL = 1),
    ev  = NULL),

  `3. Nicotine patch 21 mg (10 wk + taper)` = list(
    par = list(TQD = TQD_H, TSTARTRX = TQD_H, TSTOPRX = WK12, COUNSEL = 1),
    ev  = patch_ev(21, TQD_H, 12)),

  `4. Combination NRT (patch 21 mg + 2 mg gum PRN)` = list(
    par = list(TQD = TQD_H, TSTARTRX = TQD_H, TSTOPRX = WK12, COUNSEL = 1),
    ev  = rbind(patch_ev(21, TQD_H, 12), prn_nrt_ev(2, TQD_H, 12))),

  `5. Varenicline 1 mg BID (standard, TQD wk 1)` = list(
    par = list(TQD = TQD_H, TSTARTRX = TQD_H - 168, TSTOPRX = WK12, COUNSEL = 1),
    ev  = var_ev(TQD_H - 168, 13)),

  `6. Varenicline preloading (4 wk before TQD)` = list(
    par = list(TQD = TQD_H + 504, TSTARTRX = TQD_H - 168, TSTOPRX = WK12 + 504,
               COUNSEL = 1),
    ev  = var_ev(TQD_H - 168, 16)),

  `7. Bupropion SR 150 mg BID` = list(
    par = list(TQD = TQD_H, TSTARTRX = TQD_H - 168, TSTOPRX = WK12, COUNSEL = 1),
    ev  = bup_ev(TQD_H - 168, 13)),

  `8. Cytisinicline 3 mg TID (12 wk)` = list(
    par = list(TQD = TQD_H, TSTARTRX = TQD_H - 168, TSTOPRX = WK12, COUNSEL = 1),
    ev  = cyt_ev(TQD_H - 168, 13)),

  `9. Varenicline + nicotine patch (combination)` = list(
    par = list(TQD = TQD_H, TSTARTRX = TQD_H - 168, TSTOPRX = WK12, COUNSEL = 1),
    ev  = rbind(var_ev(TQD_H - 168, 13), patch_ev(21, TQD_H, 12))),

  `10. Patch 21 mg in a SLOW CYP2A6 metabolizer (F2A6 = 0.35)` = list(
    par = list(TQD = TQD_H, TSTARTRX = TQD_H, TSTOPRX = WK12, COUNSEL = 1,
               F2A6 = 0.35),
    ev  = patch_ev(21, TQD_H, 12))
)

run_scenario <- function(name, s, end = TEND, delta = 1) {
  m <- mod
  if (length(s$par)) m <- param(m, s$par)
  ev <- s$ev
  out <- if (is.null(ev)) {
    mrgsim(m, end = end, delta = delta, hmax = 0.5)
  } else {
    mrgsim(m, data = ev[order(ev$time), ], end = end, delta = delta,
           hmax = 0.5, recover = TRUE)
  }
  as.data.frame(out) |> mutate(scenario = name)
}

# ---------------------------------------------------------------------------
# 4) RUN + SUMMARISE
# ---------------------------------------------------------------------------
run_all <- function() {
  res <- bind_rows(lapply(names(scenarios),
                          function(n) run_scenario(n, scenarios[[n]])))
  res
}

## IMPORTANT: report 24-h MEANS, not point samples. Patch and oral NRT have
## strong diurnal swings, and sampling at an exact multiple of 24 h lands in
## the trough — which understates exposure by an order of magnitude.
mean24 <- function(time, x, centre) {
  k <- abs(time - centre) <= 12
  if (!any(k)) return(NA_real_)
  mean(x[k], na.rm = TRUE)
}

summarise_endpoints <- function(res) {
  res |>
    group_by(scenario) |>
    summarise(
      # PK / exposure  (24-h means)
      nic_ngml_wk4   = mean24(time, NIC_NGML, TQD_H + 4*7*24),
      nic_ngml_pre   = mean24(time, NIC_NGML, TQD_H - 24),
      cot_ngml_pre   = mean24(time, COT_NGML, TQD_H - 24),
      NMR            = mean24(time, NMR_,     TQD_H - 24),
      cpd_titrated   = mean24(time, CPDT_,    TQD_H - 24),
      # receptor pharmacology
      occ_tot_wk4    = mean24(time, OCC_TOT,  TQD_H + 4*7*24),
      des_wk4        = mean24(time, DES,      TQD_H + 4*7*24),
      rup_wk4        = mean24(time, RUP,      TQD_H + 4*7*24),
      rup_wk12       = mean24(time, RUP,      TQD_H + 12*7*24),
      da_wk4         = mean24(time, DA,       TQD_H + 4*7*24),
      def_wk4        = mean24(time, DEF_,     TQD_H + 4*7*24),
      # symptoms
      mnws_peak      = max(WD[time >= TQD_H & time <= TQD_H + 28*24]),
      mnws_peak_day  = (time[which.max(replace(WD, time < TQD_H |
                          time > TQD_H + 28*24, NA))] - TQD_H)/24,
      qsu_peak       = max(QSU[time >= TQD_H & time <= TQD_H + 28*24]),
      qsu_wk4        = mean24(time, QSU,      TQD_H + 4*7*24),
      # primary + secondary endpoints
      CAR_wk9_12     = approx(time, CAR_PCT,  TQD_H + 12*7*24)$y,
      PP_wk24        = approx(time, CAR_PCT,  TQD_H + 24*7*24)$y,
      PP_wk52        = approx(time, CAR_PCT,  TQD_H + 52*7*24)$y,
      co_ppm_wk12    = approx(time, CO_PPM,   TQD_H + 12*7*24)$y,
      wt_kg_wk52     = approx(time, WT,       TQD_H + 52*7*24)$y,
      fev_pct_wk52   = approx(time, FEV_PCT,  TQD_H + 52*7*24)$y,
      adherence_wk12 = approx(time, ADH_PCT,  TQD_H + 12*7*24)$y,
      nausea_peak    = max(NAUS),
      .groups = "drop"
    ) |>
    mutate(across(where(is.numeric), \(x) round(x, 3)))
}

# ---------------------------------------------------------------------------
# 4b) POPULATION AVERAGING OVER CYP2A6 ACTIVITY
#     A fixed-dose NRT arm is heterogeneous by construction, so the trial-
#     comparable CAR is a population average, not a single-subject value.
#     Weights approximate the NMR quartile spread of Lerman 2015.
# ---------------------------------------------------------------------------
F2A6_DIST <- data.frame(
  F2A6   = c(0.35, 0.65, 1.00, 1.50, 2.00),
  weight = c(0.10, 0.20, 0.40, 0.20, 0.10),
  label  = c("slow", "intermediate", "normal", "fast", "very fast")
)

car_population <- function(name, at_week = 12, delta = 6) {
  s0  <- scenarios[[name]]
  tq  <- if (is.finite(s0$par$TQD) && s0$par$TQD < 1e5) s0$par$TQD else TQD_H
  end <- tq + at_week*7*24
  per <- vapply(F2A6_DIST$F2A6, function(f) {
    s <- s0; s$par$F2A6 <- f
    r <- run_scenario(name, s, end = end, delta = delta)
    tail(r$CAR_PCT, 1)
  }, numeric(1))
  list(population = sum(F2A6_DIST$weight*per),
       per_stratum = setNames(round(per, 2), F2A6_DIST$label))
}

## Trial anchors used for the behavioural-parameter fit (EAGLES + ORCA-2 +
## Cochrane). Reproduced here so `validate_against_trials()` is self-contained.
TRIAL_TARGETS <- c(
  `2. Unaided abrupt quit (placebo + brief advice)`     = 12.5,
  `3. Nicotine patch 21 mg (10 wk + taper)`             = 23.4,
  `4. Combination NRT (patch 21 mg + 2 mg gum PRN)`     = 28.0,
  `5. Varenicline 1 mg BID (standard, TQD wk 1)`        = 33.5,
  `7. Bupropion SR 150 mg BID`                          = 22.6,
  `8. Cytisinicline 3 mg TID (12 wk)`                   = 32.0,
  `9. Varenicline + nicotine patch (combination)`       = 38.0
)

validate_against_trials <- function() {
  out <- lapply(names(TRIAL_TARGETS), function(n) {
    cp <- car_population(n, 12)
    data.frame(arm = n, model_CAR = round(cp$population, 1),
               trial_CAR = TRIAL_TARGETS[[n]],
               diff = round(cp$population - TRIAL_TARGETS[[n]], 1))
  })
  res <- bind_rows(out)
  res$RMSE <- round(sqrt(mean(res$diff^2)), 2)
  res
}

# ---------------------------------------------------------------------------
# 5) FOCUSED EXPERIMENTS
# ---------------------------------------------------------------------------

## 5a) The Brody 2006 occupancy calibration check ---------------------------
## Analytic: beta2* occupancy = C/(C + KDNIC) with KDNIC = 5.4 nM.
## 0.87 ng/mL = 5.36 nM  -> 50% occupancy (Brody's reported EC50)
## 1 cigarette venous peak ~20 ng/mL = 123 nM -> 95.8% (Brody: 88%, arterial
## kinetics and the 3-min PET frame make the in vivo value lower).
occupancy_curve <- function() {
  ngml <- 10^seq(-1.5, 2, length.out = 60)
  nM   <- ngml/0.16223
  data.frame(nic_ngml = ngml, occupancy = nM/(nM + 5.4))
}

## 5b) NMR x treatment interaction (Lerman 2015 pharmacogenetic trial) ------
## ONE parameter (F2A6) is varied. It sets nicotine clearance, hence the NMR
## biomarker, hence the TITRATED cigarette consumption (fast metabolizers
## smoke more and so adapt to a higher set-point), hence the fraction of that
## set-point a FIXED 21 mg patch can replace. Varenicline is renally cleared
## and untouched. No interaction term is fitted anywhere in the model.
##
## The key emergent quantity is the REPLACEMENT RATIO: patch nicotine at week
## 4 divided by the subject own pre-quit nicotine level.
nmr_experiment <- function() {
  arms <- list(patch = scenarios[[3]], varenicline = scenarios[[5]])
  out  <- list()
  for (a in names(arms)) for (i in seq_len(nrow(F2A6_DIST))) {
    f <- F2A6_DIST$F2A6[i]
    s <- arms[[a]]; s$par$F2A6 <- f
    r <- run_scenario(sprintf("%s | %s", a, F2A6_DIST$label[i]), s,
                      end = WK12, delta = 4)
    pre <- mean24(r$time, r$NIC_NGML, TQD_H - 24)
    wk4 <- mean24(r$time, r$NIC_NGML, TQD_H + 4*7*24)
    out[[length(out) + 1]] <- data.frame(
      arm = a, metabolizer = F2A6_DIST$label[i], F2A6 = f,
      cpd_titrated   = mean24(r$time, r$CPDT_,   TQD_H - 24),
      NMR            = mean24(r$time, r$NMR_,    TQD_H - 24),
      nic_pre        = pre,
      nic_wk4        = wk4,
      replacement    = wk4/pre,
      occ_wk4        = mean24(r$time, r$OCC_TOT, TQD_H + 4*7*24),
      def_wk4        = mean24(r$time, r$DEF_,    TQD_H + 4*7*24),
      mnws_peak      = max(r$WD[r$time >= TQD_H]),
      CAR_wk12       = tail(r$CAR_PCT, 1)
    )
  }
  bind_rows(out) |> mutate(across(where(is.numeric), \(x) round(x, 3)))
}

## 5c) Intrinsic activity sweep: is PARTIAL agonism an efficacy optimum? -----
## Two terms are at stake as EMAXV goes 0 -> 1:
##   * RELIEF   (via DA -> DEF -> MNWS/QSU) scales WITH intrinsic activity.
##   * BLOCKADE (LAPREW = 1 - GBLOCK*OCCT - GBLOCKV*(FV+FC)) depends only on
##     OCCUPANCY, which is INDEPENDENT of intrinsic activity.
## RESULT (computed, 2026): there is NO interior optimum. CAR rises
## monotonically 32.8% -> 39.4% while lapse_reward stays pinned at 0.227.
## The model therefore REFUTES the common claim that partial agonism is
## optimal *for efficacy* and relocates the argument: partial activity plus
## non-pulsatile oral delivery is what limits the drug own abuse liability,
## which this model does not score. Report both columns and read them
## together — the flat lapse_reward column is the substantive finding.
partial_agonism_sweep <- function() {
  ia <- seq(0, 1, by = 0.1)
  out <- lapply(ia, function(e) {
    s <- scenarios[[5]]; s$par$EMAXV <- e
    r <- run_scenario(sprintf("EMAXV = %.2f", e), s, end = WK12, delta = 4)
    occv <- mean24(r$time, r$OCC_VAR, TQD_H + 4*7*24)
    occt <- mean24(r$time, r$OCC_TOT, TQD_H + 4*7*24)
    pp   <- as.list(param(mod))
    data.frame(EMAXV        = e,
               occ_var_wk4  = occv,
               da_wk4       = mean24(r$time, r$DA,   TQD_H + 4*7*24),
               def_wk4      = mean24(r$time, r$DEF_, TQD_H + 4*7*24),
               mnws_peak    = max(r$WD[r$time >= TQD_H]),
               qsu_wk4      = mean24(r$time, r$QSU,  TQD_H + 4*7*24),
               lapse_reward = max(0.05, 1 - pp$GBLOCK*occt - pp$GBLOCKV*occv),
               CAR_wk12     = tail(r$CAR_PCT, 1))
  })
  bind_rows(out) |> mutate(across(where(is.numeric), \(x) round(x, 3)))
}

## 5d) Receptor-normalization and withdrawal-window timescales ---------------
## The model is never told that withdrawal lasts 2-4 weeks or that beta2* PET
## normalizes in 3-4 weeks. Both come out of KOUTR / KSET (t1/2 = 14 d).
receptor_normalization <- function() {
  r <- run_scenario("unaided", scenarios[[2]], end = TQD_H + 12*7*24, delta = 6)
  post <- r[r$time >= TQD_H, ]
  r0   <- post$RUP[1]
  frac <- (post$RUP - 1)/(r0 - 1)
  pk   <- max(post$WD)
  above <- function(fr) {
    idx <- which(post$WD >= fr*pk)
    round((post$time[max(idx)] - TQD_H)/24, 1)
  }
  data.frame(
    RUP_at_quit      = round(r0, 3),
    pct_upregulation = round(100*(r0 - 1), 1),
    RUP_t50_days     = round((post$time[which(frac <= 0.50)[1]] - TQD_H)/24, 1),
    RUP_t90_days     = round((post$time[which(frac <= 0.10)[1]] - TQD_H)/24, 1),
    RUP_wk4          = round(mean24(post$time, post$RUP, TQD_H + 4*7*24), 3),
    MNWS_peak        = round(pk, 3),
    MNWS_peak_day    = round((post$time[which.max(post$WD)] - TQD_H)/24, 1),
    MNWS_above50_d   = above(0.50),
    MNWS_above25_d   = above(0.25),
    MNWS_above10_d   = above(0.10)
  )
}

## 5e) Very-low-nicotine cigarettes (reduced-nicotine product standard) -----
## Because beta2* occupancy saturates at an EC50 of 0.87 ng/mL, nicotine
## content has to fall very far before any dependence-relevant variable
## moves. The sweep locates that threshold instead of assuming it.
## CAVEAT: compensatory smoking (deeper puffs, more cigarettes) is NOT
## modelled here, so these are upper bounds on the benefit of a VLNC standard.
vlnc_experiment <- function() {
  yields <- c(1.1, 0.55, 0.25, 0.10, 0.05, 0.03, 0.015)
  out <- lapply(yields, function(mg) {
    s <- scenarios[[2]]; s$par$NICCIG <- mg
    r <- run_scenario(sprintf("%.3f mg/cig", mg), s, end = WK12, delta = 4)
    data.frame(nic_per_cig_mg = mg,
               nic_ngml_pre = mean24(r$time, r$NIC_NGML, TQD_H - 24),
               occ_pre      = mean24(r$time, r$OCC_TOT,  TQD_H - 24),
               des_pre      = mean24(r$time, r$DES,      TQD_H - 24),
               RUP_pre      = mean24(r$time, r$RUP,      TQD_H - 24),
               SETP_pre     = mean24(r$time, r$SETP,     TQD_H - 24),
               mnws_peak    = max(r$WD[r$time >= TQD_H]),
               CAR_wk12     = tail(r$CAR_PCT, 1))
  })
  bind_rows(out) |> mutate(across(where(is.numeric), \(x) round(x, 3)))
}

## 5f) A FALSIFIABLE PREDICTION: varenicline + patch is NOT additive --------
## Nicotine and varenicline compete for the same site. Adding a 21 mg patch to
## varenicline raises TOTAL occupancy only slightly (it is already ~0.90) while
## substantially reducing varenicline OWN occupancy — and it is varenicline
## own occupancy that carries the partial-agonist-specific blockade term
## GBLOCKV. The model therefore predicts the combination is roughly equal to,
## not clearly better than, varenicline alone.
##
## This DISAGREES with Koegelenberg 2014 (JAMA), which reported a benefit, and
## AGREES with the later neutral trials. The model names the mechanism that
## would decide it: measure varenicline beta2* occupancy with and without a
## concurrent patch. If occupancy is not displaced, the model is wrong here.
combination_displacement <- function() {
  arms <- c(`varenicline alone` = 5, `varenicline + patch` = 9,
            `patch alone` = 3)
  out <- lapply(names(arms), function(n) {
    r <- run_scenario(n, scenarios[[arms[[n]]]], end = WK12, delta = 4)
    data.frame(arm = n,
               occ_total_wk4 = mean24(r$time, r$OCC_TOT, TQD_H + 4*7*24),
               occ_var_wk4   = mean24(r$time, r$OCC_VAR, TQD_H + 4*7*24),
               occ_nic_wk4   = mean24(r$time, r$OCC_NIC, TQD_H + 4*7*24),
               da_wk4        = mean24(r$time, r$DA,      TQD_H + 4*7*24),
               mnws_peak     = max(r$WD[r$time >= TQD_H]),
               CAR_wk12      = tail(r$CAR_PCT, 1))
  })
  bind_rows(out) |> mutate(across(where(is.numeric), \(x) round(x, 3)))
}

# ---------------------------------------------------------------------------
# 6) PLOTS (optional; requires ggplot2)
# ---------------------------------------------------------------------------
plot_panel <- function(res) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  library(ggplot2)
  long <- res |>
    mutate(day = (time - TQD_H)/24) |>
    filter(day >= -14, day <= 120) |>
    select(day, scenario, NIC_NGML, OCC_TOT, RUP, DA, WD, QSU, CAR_PCT, WT) |>
    pivot_longer(-c(day, scenario))
  ggplot(long, aes(day, value, colour = scenario)) +
    geom_line(linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
    facet_wrap(~name, scales = "free_y", ncol = 2) +
    labs(x = "Days from target quit date", y = NULL,
         title = "Tobacco Use Disorder QSP model — 10 therapy scenarios") +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom", legend.text = element_text(size = 6)) +
    guides(colour = guide_legend(ncol = 2))
}

# ---------------------------------------------------------------------------
# 7) REFERENCE RESULTS
#    Values produced by this parameter set (independent LSODA integration of
#    the identical ODE system, 2026-07-29). Re-running section 8 should
#    reproduce them to within solver tolerance; a large deviation means a
#    parameter or an equation has drifted.
# ---------------------------------------------------------------------------
## Baseline smoker steady state (F2A6 = 1.0, 20 cig/day, 90-day run-in)
##   plasma nicotine 12.7 ng/mL (24-h mean) · cotinine 194 ng/mL · NMR 0.39
##   a4b2* occupancy 0.94 (daytime) · desensitized fraction 0.94
##   receptor pool RUP 1.70 (= +70% beta2* upregulation) · DA 1.81 · SETP 1.86
##   dopaminergic deficit ~0 by construction · exhaled CO 25 ppm
##
## CYP2A6 phenotypes (titration: intake scales with clearance^0.5)
##   F2A6 0.35 -> 14.3 cig/day · nicotine 17.6 ng/mL · NMR 0.14
##   F2A6 1.00 -> 20.0 cig/day · nicotine 12.7 ng/mL · NMR 0.39
##   F2A6 2.00 -> 26.5 cig/day · nicotine  9.6 ng/mL · NMR 0.78
##
## Population CAR weeks 9-12 (averaged over the CYP2A6 distribution)
##   arm                    model   trial anchor
##   unaided quit            12.5    12.5   (fit target)
##   nicotine patch 21 mg    27.0    23.4
##   combination NRT         25.1    28.0
##   varenicline             35.8    33.5
##   bupropion SR            22.2    22.6
##   cytisinicline           33.8    32.0
##   varenicline + patch     33.0    38.0   <-- model predicts non-additivity
##   RMSE = 2.80 percentage points over the 7 arms
##
## Withdrawal / receptor timescales after an unaided quit
##   MNWS peak 1.58 on day 5.2; >50% of peak for 35 d; >25% for 54 d
##   RUP 1.701 at quit -> t50 15.0 d -> 1.183 at week 4 -> t90 47.5 d
##   SETP 1.862 at quit -> 1.218 at week 4 -> 1.014 at week 12
##
## Varenicline + patch displacement (week 4, 24-h means)
##   varenicline alone   OCCtot 0.897  OCCvar 0.897  OCCnic 0.000  MNWS 1.43
##   varenicline + patch OCCtot 0.944  OCCvar 0.501  OCCnic 0.443  MNWS 1.31
##   -> total occupancy +5% but varenicline own occupancy nearly HALVED
##
## Reduced-nicotine cigarettes (unaided quit, no compensation modelled)
##   mg/cig   1.100  0.550  0.250  0.100  0.050  0.030  0.015
##   occ      0.942  0.890  0.786  0.597  0.427  0.309  0.184
##   RUP      1.701  1.696  1.685  1.660  1.626  1.588  1.518
##   MNWS pk  1.58   1.54   1.48   1.33   1.15   0.97   0.71
##   CAR (%)  12.5   13.1   14.3   17.2   20.6   23.8   28.2
##
## KNOWN DISCREPANCIES — stated rather than tuned away
##   (a) Withdrawal resolves ~1.5-2x slower than the clinical 2-4 weeks.
##   (b) Combination NRT lands just BELOW patch alone (25.1 vs 27.0) whereas
##       Cochrane finds it clearly better. Cause: the phasic-reinforcement
##       gain GREINF makes PRN oral nicotine mildly hazard-increasing, which
##       here outweighs its craving relief. GREINF is a fitted parameter.
##   (c) Varenicline + patch is predicted non-additive (see 5f). This
##       disagrees with Koegelenberg 2014 and agrees with later neutral
##       trials; the model names the measurement that would settle it.
##   (d) Bupropion relieves withdrawal more than the literature supports
##       (peak MNWS 1.06 vs 1.58 unaided) because GDAT enters DRIVE
##       multiplicatively. Its CAR is nonetheless on target.
##   (e) WT is a POPULATION mean over abstainers and relapsers, so the
##       52-week value (+2.0 kg) is far below the ~5 kg seen in sustained
##       abstainers; a sustained abstainer in this model reaches WTSS = 5 kg.

# ---------------------------------------------------------------------------
# 8) MAIN
# ---------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  message("== A) 10 therapy scenarios ==")
  res <- run_all()
  print(as.data.frame(summarise_endpoints(res)), width = 250)

  message("\n== B) Validation against trial CARs (population-averaged) ==")
  print(validate_against_trials())

  message("\n== C) Receptor + withdrawal timescales after an unaided quit ==")
  print(t(receptor_normalization()))

  message("\n== D) NMR x treatment interaction (Lerman 2015) ==")
  print(nmr_experiment())

  message("\n== E) Why partial agonism (EMAXV sweep) ==")
  print(partial_agonism_sweep())

  message("\n== F) Very-low-nicotine cigarettes ==")
  print(vlnc_experiment())

  message("\n== G) Varenicline + patch: competitive displacement prediction ==")
  print(combination_displacement())

  p <- plot_panel(res)
  if (!is.null(p)) ggplot2::ggsave("tud_scenarios.png", p,
                                   width = 10, height = 12, dpi = 150)
}
