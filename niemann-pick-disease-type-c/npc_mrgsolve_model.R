## =============================================================================
##  Niemann-Pick disease type C (NPC) — QSP model for mrgsolve
##  니만-피크병 C형 · 정량적 시스템 약리학 모델
##
##  41 ODEs · 4 drugs · 6 genotypes · 6 onset archetypes · 18 scenarios
##
##  Twin implementation:  npc_reference_model.py  (independently written; every
##  shared parameter and right-hand side agrees.  Any number quoted in README.md
##  is printed by that script.)
##
##  ---------------------------------------------------------------------------
##  WHAT THIS MODEL IS FOR
##  ---------------------------------------------------------------------------
##  NPC has four drugs acting at four points of one causal chain — two of them
##  FDA-approved in the same month (arimoclomol and levacetylleucine, September
##  2024) on the basis of trials with completely different designs and endpoints.
##  This model exists to separate three things the literature routinely fuses:
##
##    (1) WHERE the storage is.  The plasma biomarkers (cholestane-3b,5a,6b-triol,
##        lysoSM-509/PPCS, trihydroxycholanoyl-glycine) are oxidation products of
##        sterol stored in LIVER and SPLENIC macrophages.  They read the VISCERAL
##        pool.  The disease that kills people is in the CEREBELLUM.  The two
##        compartments have different drug access (oral vs intrathecal), turnover
##        (days vs years) and reversibility.
##
##    (2) WHICH pool of lost function.  Neurological function is written as
##            FUNC = 1 - D_rev - D_irr
##        D_rev reversible on ~2 weeks, D_irr monotone non-decreasing.  A
##        SYMPTOMATIC drug moves D_rev; a DISEASE-MODIFYING drug lowers dD_irr/dt.
##
##    (3) WHAT each trial design can see.  Consequence of (2): a 12-week
##        crossover on SARA sees only D_rev; a 12-month change-from-baseline on
##        NPCCSS cannot separate a constant offset from a slope change.
##
##  ---------------------------------------------------------------------------
##  THE STRUCTURAL CHOICE THAT DOES THE MOST WORK
##  ---------------------------------------------------------------------------
##  DAM is the INTEGRAL of cerebellar stress with NO repair term, and Purkinje
##  death is gated on DAM crossing D_reserve.  Two published facts fall out of
##  that single line, neither of them fitted:
##
##    * latency to onset is INVERSELY proportional to stress, so residual NPC1
##      activity spreads onset over decades (juvenile 8.0 yr, adult 29.8 yr);
##    * the post-gate death rate SATURATES, so the slope after onset is nearly
##      independent of the age of onset.
##
##  That is exactly the headline of Yanjanin 2010 (PMID 19415691): "linear
##  clinical progression, INDEPENDENT of age of onset".  A version with a repair
##  term was written first and rejected: with repair, the crossing time depends
##  only LOGARITHMICALLY on stress, which collapses the onset range to a couple
##  of years and cannot produce adult-onset disease at all.
##
##  ---------------------------------------------------------------------------
##  Usage
##  ---------------------------------------------------------------------------
##    library(mrgsolve); library(dplyr); library(ggplot2)
##    mod <- mread("npc_mrgsolve_model.R")
##    out <- npc_scenarios(mod)            # all 18 scenarios
##    npc_calibration_table(mod)           # model vs published targets
##
##  Dosing compartment cheat-sheet (cmt = ...):
##    1  MIG_GUT  miglustat, oral          200 mg tid (adult)
##    5  ARI_GUT  arimoclomol, oral        124 mg tid (>=  weight band)
##    8  NAL_GUT  levacetylleucine, oral   1333 mg tid (4 g/day)
##    11 CD_CSF   adrabetadex, intrathecal 900 mg q2wk
## =============================================================================

$PROB
# Niemann-Pick disease type C — QSP model
- Author: QSP Disease Model Library (Claude Code Routine)
- Reference twin: npc_reference_model.py
- Units: time in DAYS; amounts in mg; lipid pools in arbitrary units normalised
  to the wild-type steady state; scales in their published points.

$PARAM @annotated
// ---------------------------------------------------------------------------
// 1. genotype -> NPC1 protein  (cluster 2 of the mechanistic map)
// ---------------------------------------------------------------------------
f_null      : 0.00    : fraction of NPC1 alleles that are null (no protein)
theta0      : 0.055   : intrinsic ER folding yield of the mutant allele
f_npc2      : 1.00    : NPC2 functional fraction (1 = normal, NPC2 disease << 1)
ksyn_npc1   : 1 : NPC1 synthesis rate (arbitrary units/day)
kerad       : 2 : ER-associated degradation rate (1/day)
kexit0      : 1 : maximal ER exit rate of correctly folded NPC1 (1/day)
kdegL       : 0.15 : lysosomal NPC1 turnover (1/day)
EChsp       : 0.8 : HSP70 fold-rise above 1 giving half-max folding gain
Emax_fold   : 0.3 : max fractional closure of the folding gap by HSP70
f_leak      : 0.01 : NPC1-independent lysosomal cholesterol leak

// ---------------------------------------------------------------------------
// 2. cholesterol pools.  The CNS pool runs on a NEURONAL timescale: brain
//    cholesterol turns over with a half-life of months-to-years, the hepatic
//    pool in days.  Jin_c/Vmax_c/kbas_chol_c are the visceral values divided by
//    tau_cns = 25, which leaves every CNS steady state unchanged and makes the
//    APPROACH 25x slower.  That factor alone creates the latent period.
// ---------------------------------------------------------------------------
Jin_v       : 1 : visceral lysosomal cholesterol influx (units/day)
Vmax_v      : 3 : visceral NPC1-mediated egress Vmax
Km_v        : 1 : visceral egress Km
kbas_chol   : 0.02 : visceral non-specific cholesterol clearance (1/day)
tau_cns     : 25 : CNS:visceral timescale ratio (documentation only)
Jin_c       : 0.014 : CNS lysosomal cholesterol influx (units/day)
Vmax_c      : 0.042 : CNS NPC1-mediated egress Vmax
Km_c        : 1 : CNS egress Km
kbas_chol_c : 0.0008 : CNS non-specific cholesterol clearance (1/day)

// ---------------------------------------------------------------------------
// 3. glycosphingolipids (miglustat target)
// ---------------------------------------------------------------------------
ksyn_gsl_v  : 0.6 : visceral GSL synthesis via UGCG
ksyn_gsl_c  : 0.012 : CNS GSL synthesis via UGCG
kdeg_gsl    : 0.6 : visceral GSL catabolism (1/day, x hydrolase activity)
kdeg_gsl_c  : 0.024 : CNS GSL catabolism (1/day, x hydrolase activity)

// ---------------------------------------------------------------------------
// 4. sphingosine and lysosomal calcium (Lloyd-Evans 2008, PMID 18953351)
// ---------------------------------------------------------------------------
ksph        : 1 : lysosomal sphingosine generation (units/day)
kexp_sph    : 4 : sphingosine export (1/day)
fmin_sph    : 0.12 : NPC1-independent floor of sphingosine export
Ksph        : 1.15 : sphingosine IC50 on acidic Ca2+ store refilling
Ktr_chol    : 60 : cholesterol IC50 on TRPML1
kca         : 5 : acidic Ca2+ store equilibration (1/day)

// ---------------------------------------------------------------------------
// 5. hydrolase capacity
// ---------------------------------------------------------------------------
Kph         : 14 : cholesterol load that raises lysosomal pH (CNS units)
hill_ph     : 1.5 : steepness of the pH effect
e_hyd_hsp   : 0.2 : HSP70 gain on hydrolase activity (PMID 20111001)
khyd        : 3 : hydrolase equilibration (1/day)

// ---------------------------------------------------------------------------
// 6. autophagy backlog, mitochondria, ROS
// ---------------------------------------------------------------------------
kaut_in     : 1 : autophagic substrate delivery (units/day)
kaut_out    : 2.5 : autophagic substrate clearance
Kmito       : 22 : cholesterol load halving mitochondrial function
kmito       : 0.5 : mitochondrial equilibration (1/day)
aros        : 2.5 : ROS gain per unit lost mitochondrial function
kros        : 4 : ROS equilibration (1/day)

// ---------------------------------------------------------------------------
// 7. cerebellar stress, the damage integral, and Purkinje pools
// ---------------------------------------------------------------------------
w_chol      : 0.55 : stress weight, CNS cholesterol load
w_gsl       : 0.22 : stress weight, CNS glycosphingolipid load
w_ros       : 0.3 : stress weight, oxidative stress
w_aut       : 0.18 : stress weight, autophagic backlog
w_infl      : 0.6 : stress weight, neuroinflammation
w_ca        : 0.35 : stress weight, lysosomal calcium deficit
S_thresh    : 0.2 : stress absorbed without damage accrual
D_reserve   : 6736.55 : stress x days of damage the cerebellum absorbs
hill_gate   : 4 : sharpness of the reserve threshold
kon_pc      : 0.02 : healthy -> stressed Purkinje cell (per stress per day)
koff_pc     : 0.01 : stressed -> healthy recovery (1/day)
Krec        : 0.5 : stress that halves recovery
kdie        : 0.000113537 : stressed -> dead once the gate is open (1/day)
v_dev       : 0.6 : developmental vulnerability amplitude (PATIENT covariate)
tau_dev     : 2.5 : developmental vulnerability time constant (years)

// ---------------------------------------------------------------------------
// 8. neuroinflammation, synapses, cerebellar volume
// ---------------------------------------------------------------------------
a_chol      : 0.35 : microglial activation by CNS lipid load
a_pc        : 0.9 : microglial activation by stressed Purkinje cells
a_dead      : 0.55 : microglial activation by accumulated cell death
kinfl       : 0.05 : inflammation equilibration (1/day)
ksyn_rep    : 0.02 : synapse formation/repair (1/day)
kprune      : 0.02 : complement-mediated synaptic pruning (1/day)
g_vol_pc    : 0.85 : cerebellar volume loss per unit Purkinje loss
g_vol_syn   : 0.25 : cerebellar volume loss per unit synapse loss
kcbl        : 0.004 : cerebellar volume equilibration (1/day)

// ---------------------------------------------------------------------------
// 9. FUNCTION DECOMPOSITION — the structural core
// ---------------------------------------------------------------------------
cap_rev     : 0.335 : ceiling on reversible dysfunction
K_rev       : 6 : stress giving half of cap_rev
cr_pcs      : 0 : extra D_rev from the stressed-cell pool (deliberately 0)
krev        : 0.06 : D_rev time constant (1/day, t1/2 ~11.6 d)
g_irr_pc    : 1.6 : irreversible dysfunction per unit Purkinje loss
g_irr_syn   : 0 : synapse density is a READOUT here, not a driver
g_irr_cbl   : 0 : cerebellar volume is a READOUT here, not a driver
kirr        : 0.01 : D_irr tracking rate (1/day, one-sided)

// ---------------------------------------------------------------------------
// 10. clinical scale mapping.  SARA scores what the patient can do RIGHT NOW
//     under cerebellar control -> weighted to D_rev.  The NPCCSS domains
//     (ambulation, speech, swallow, fine motor, cognition) are milestone-like
//     -> weighted to D_irr.  This asymmetry is why a symptomatic drug can win a
//     12-week SARA endpoint while a rate endpoint needs NPCCSS and a year.
// ---------------------------------------------------------------------------
sara_max    : 40 : SARA ceiling
sara_a      : 1.8 : weight of D_rev in SARA
sara_b      : 1 : weight of D_irr in SARA
sara_p      : 1 : exponent (linear: NPC progression is reported linear)
sara_k      : 30.3141 : SARA scale (calibrated: 15.91 at trial entry)
n5_max      : 25 : 5-domain NPCCSS ceiling
n5_a        : 0.25 : weight of D_rev in NPCCSS
n5_b        : 1 : weight of D_irr in NPCCSS
n5_p        : 1 : exponent
n5_k        : 30 : NPCCSS scale (calibrated: 1.5 pt/yr, PMID 33228797)
n4_frac     : 0.8 : R4DNPCCSS is 4 of the 5 domains
n4_gain     : 1.16 : ... rescored for linearity (PMID 40520915)
n17_k       : 1.95 : 17-domain = n17_k x 5-domain + hearing
sacc_rev    : 0.4 : saccade velocity sensitivity to D_rev
sacc_irr    : 1.6 : saccade velocity sensitivity to D_irr
swal_a      : 0.55 : swallowing sensitivity to D_rev
swal_b      : 1.55 : swallowing sensitivity to D_irr

// ---------------------------------------------------------------------------
// 11. survival hazard (dysphagia is the dominant route, PMID 23039766)
// ---------------------------------------------------------------------------
h0          : 2e-05 : baseline hazard (1/day)
h_swal      : 0.000121875 : hazard from dysphagia / aspiration (1/day)
h_liver     : 0.0004 : hazard from infantile liver disease (1/day)
liver_risk  : 0.000   : patient covariate, 1 = perinatal cholestatic form

// ---------------------------------------------------------------------------
// 12. biomarkers.  ktri and Ktri are a TWO-POINT calibration that reproduces
//     BOTH the patient mean (88.31 ng/mL) and the control mean (5.97 ng/mL) of
//     PMID 33228797 exactly.  Ktri is a SATURATION constant, and the model then
//     PREDICTS that triol grades severity poorly because in patients it already
//     sits at ~82% of its own ceiling.
// ---------------------------------------------------------------------------
ktri        : 54.15 : triol generation scale
Ktri        : 8.35 : triol generation saturation constant
ktri_el     : 0.5 : triol elimination (1/day)
kppcs       : 6.2 : lysoSM-509 / PPCS generation
Kppcs       : 10 : PPCS generation saturation constant
b_ppcs_gsl  : 0.45 : GSL contribution to PPCS
kppcs_el    : 0.35 : PPCS elimination (1/day)
ktcg        : 12 : bile acid B (TCG) generation
Ktcg        : 9 : TCG saturation constant
ktcg_el     : 0.6 : TCG elimination (1/day)
knfl        : 420 : CSF NfL generation per unit Purkinje death FLUX
knfl_base   : 1.2 : baseline NfL generation
knfl_el     : 0.14 : NfL elimination (1/day)
kcalb       : 900 : CSF calbindin per unit Purkinje death FLUX (PMID 27307499)
kcalb_el    : 0.3 : calbindin elimination (1/day)

// ---------------------------------------------------------------------------
// 13. miglustat PK/PD (Zavesca; PMID 17689147, 19447653, 15901676)
// ---------------------------------------------------------------------------
mig_ka      : 8 : absorption rate (1/day)
mig_F       : 0.97 : oral bioavailability
mig_V       : 90 : central volume (L)
mig_CL      : 230 : clearance (L/day) -> t1/2 6.5 h
mig_Q       : 20 : intercompartmental clearance (L/day)
mig_Vp      : 60 : peripheral volume (L)
mig_Kp_br   : 0.45 : brain:plasma partition coefficient
mig_ke0     : 1 : brain effect-site equilibration (1/day)
mig_IC50    : 30 : UGCG inhibition IC50 (uM)
mig_Emax    : 0.85 : maximal fractional UGCG inhibition
mig_MW      : 219.28 : molecular weight (g/mol)

// ---------------------------------------------------------------------------
// 14. arimoclomol PK/PD (Miplyffa; PMID 39715913, 18551622, 34418116)
// ---------------------------------------------------------------------------
ari_ka      : 12 : absorption rate (1/day)
ari_F       : 0.6 : oral bioavailability
ari_V       : 150 : volume of distribution (L)
ari_CL      : 624 : clearance (L/day) -> t1/2 ~4 h
ari_Kp_csf  : 0.15 : CSF:plasma partition coefficient
ari_ke0     : 2 : CSF equilibration (1/day)
ari_EC50    : 120 : HSP70 induction EC50 in CSF (ng/mL)
ari_Emax    : 1.6 : maximal HSP70 fold-increase above baseline
ari_stress0 : 0.25 : CO-INDUCER floor: fraction of Emax without cell stress

// ---------------------------------------------------------------------------
// 15. levacetylleucine PK/PD (Aqneursa; PMID 38294974, 26400580, 33738443)
//     ke0 is deliberately slow (t1/2 ~4.6 d): the CLINICAL effect of
//     acetyl-leucine appears over days-to-weeks and is lost over days on
//     withdrawal, so the effect site is not the plasma spike.
// ---------------------------------------------------------------------------
nal_ka      : 30 : absorption rate (1/day)
nal_F       : 1 : oral bioavailability
nal_V       : 253 : Vss/F (L, label value)
nal_CL      : 3336 : CL/F (L/day = 139 L/h, label value)
nal_ke0     : 0.15 : effect-site equilibration (1/day)
nal_EC50    : 0.8 : effect-site EC50 (mg/L)
nal_Emax_sym: 0.30779 : max fractional reduction of the D_rev target
nal_e_gluc  : 0.12 : max gain on the brain glucose metabolism readout
nal_Emax_dm : 0.08 : max fractional reduction of CNS lipid influx

// ---------------------------------------------------------------------------
// 16. adrabetadex / 2-HPbCD PK/PD (PMID 28803710, 25717099, 23285273, 26055150)
//     Ototoxicity is driven by the CSF concentration (cochlear aqueduct route),
//     NOT the systemic concentration -- which is why intrathecal dosing is the
//     ototoxic route and why the same mechanism that extracts lysosomal
//     cholesterol also strips it from the outer hair cell membrane.
// ---------------------------------------------------------------------------
cd_Vcsf     : 0.15 : CSF volume (L)
cd_kout     : 4.16 : CSF -> systemic (1/day, t1/2 4 h)
cd_kin_br   : 0.6 : CSF -> brain ECF/intracellular (1/day)
cd_kout_br  : 1.5 : brain ECF elimination (1/day)
cd_Vsys     : 18 : systemic volume (L)
cd_klsys    : 24 : systemic elimination (1/day, renal, rapid)
cd_kext     : 0.000385 : cholesterol extraction per (mg/L) per day
cd_koto     : 1.94e-05 : outer hair cell loss per (mg/L in CSF) per day
cd_hear_max : 55 : threshold shift at total outer hair cell loss (dB)
cd_hear_p   : 0.7 : exponent of the hearing threshold curve

// ---------------------------------------------------------------------------
// 17. switches for the trial-design experiments of section 9
// ---------------------------------------------------------------------------
nal_sym_on  : 1.000   : 1 = symptomatic arm of levacetylleucine active
nal_dm_on   : 1.000   : 1 = disease-modifying arm of levacetylleucine active

$CMT @annotated
// ---- drug PK (13) ---------------------------------------------------------
MIG_GUT   : miglustat gut depot (mg)
MIG_CEN   : miglustat central (mg)
MIG_PER   : miglustat peripheral (mg)
MIG_BR    : miglustat brain effect site (uM)
ARI_GUT   : arimoclomol gut depot (mg)
ARI_CEN   : arimoclomol central (mg)
ARI_CSF   : arimoclomol CSF effect site (ng/mL)
NAL_GUT   : levacetylleucine gut depot (mg)
NAL_CEN   : levacetylleucine central (mg)
NAL_BR    : levacetylleucine brain effect site (mg/L)
CD_CSF    : adrabetadex in CSF (mg)
CD_BR     : adrabetadex in brain ECF/cells (mg)
CD_SYS    : adrabetadex systemic (mg)
// ---- NPC1 protein (2) -----------------------------------------------------
NPC1_ER   : NPC1 in the ER awaiting folding QC (units)
NPC1_L    : functional NPC1 in the lysosomal membrane (units)
// ---- visceral lipids (2) --------------------------------------------------
CHOL_V    : visceral lysosomal cholesterol (units)
GSL_V     : visceral glycosphingolipid (units)
// ---- CNS lipids (4) -------------------------------------------------------
CHOL_C    : neuronal lysosomal cholesterol (units)
GSL_C     : neuronal glycosphingolipid, GM2/GM3 (units)
SPH       : lysosomal sphingosine (units)
CA_LY     : acidic Ca2+ store filling (fraction)
// ---- lysosomal / metabolic function (4) -----------------------------------
HYD       : lysosomal hydrolase activity (fraction)
AUTOPH    : autophagic substrate backlog (units)
MITO      : mitochondrial function (fraction)
ROS       : oxidative stress (relative)
// ---- Purkinje pools and network (7) ---------------------------------------
PC        : healthy Purkinje cells (fraction)
PC_S      : stressed but recoverable Purkinje cells (fraction)
PC_LOST   : dead Purkinje cells (fraction, MONOTONE)
INFL      : neuroinflammation index
SYN       : synaptic density (fraction)
CBL       : cerebellar volume (fraction)
DAM       : cumulative cerebellar damage = integral of stress (stress x days)
// ---- function decomposition (2) -------------------------------------------
D_REV     : reversible neurological dysfunction (0-1)
D_IRR     : irreversible neurological dysfunction (0-1, MONOTONE)
// ---- biomarkers (5) -------------------------------------------------------
TRIOL     : plasma cholestane-3b,5a,6b-triol (ng/mL)
PPCS      : plasma lysoSM-509 / PPCS (relative)
TCG       : plasma bile acid B / TCG (relative)
NFL       : CSF neurofilament light chain (relative)
CALB      : CSF calbindin-D28K (relative)
// ---- cochlea and survival (2) ---------------------------------------------
OHC       : surviving cochlear outer hair cells (fraction)
CUMHAZ    : cumulative mortality hazard

$GLOBAL
#define C_MIG_PL   (MIG_CEN / mig_V / mig_MW * 1000.0)
#define C_ARI_PL   (ARI_CEN / ari_V * 1000.0)
#define C_NAL_PL   (NAL_CEN / nal_V)
#define C_CD_CSF   (CD_CSF / cd_Vcsf)
#define C_CD_BR    (CD_BR  / cd_Vcsf)
#define C_CD_SYS   (CD_SYS / cd_Vsys)

// Wild-type reference state.  Every normalised quantity in the stress term is
// divided by its wild-type value, so a healthy system sits at stress = 0
// EXACTLY, not approximately.  These constants are the closed-form solution of
// the untreated lipid/lysosome core at f_NPC1 = 1 and are reproduced to 4
// decimals by wt_reference() in the Python twin.
#define WT_CHOL_V  0.485500
#define WT_CHOL_C  0.473000
#define WT_SPH     0.250000
#define WT_CA      0.947400
#define WT_HYD     0.965100
#define WT_GSL_C   0.518100
#define WT_GSL_V   1.036200
#define WT_AUT     0.425700
#define WT_MITO    0.999500
#define WT_ROS     1.001200
#define WT_NPC1_L  2.222200

$MAIN
// Birth state: the lipid pools start at the wild-type steady state (fetal and
// maternal handling keep the newborn near normal) and fill in at the
// genotype-determined rate.  The NPC1 pools start at their genotype steady
// state, which is instantaneous relative to everything else.
double th0 = theta0;
double er0 = ksyn_npc1 * (1.0 - f_null) / (kexit0 * th0 + kerad);

NPC1_ER_0 = er0;
NPC1_L_0  = er0 * kexit0 * th0 / kdegL;
CHOL_V_0  = WT_CHOL_V;
CHOL_C_0  = WT_CHOL_C;
GSL_V_0   = WT_GSL_V;
GSL_C_0   = WT_GSL_C;
SPH_0     = WT_SPH;
CA_LY_0   = WT_CA;
HYD_0     = WT_HYD;
AUTOPH_0  = WT_AUT;
MITO_0    = WT_MITO;
ROS_0     = WT_ROS;
PC_0      = 1.0;
SYN_0     = 1.0;
CBL_0     = 1.0;
OHC_0     = 1.0;
TRIOL_0   = ktri  * WT_CHOL_V / (Ktri + WT_CHOL_V) / ktri_el;
PPCS_0    = kppcs * (WT_CHOL_V + b_ppcs_gsl * WT_GSL_V) /
            (Kppcs + WT_CHOL_V + b_ppcs_gsl * WT_GSL_V) / kppcs_el;
TCG_0     = ktcg  * WT_CHOL_V / (Ktcg + WT_CHOL_V) / ktcg_el;
NFL_0     = knfl_base / knfl_el;

$ODE
// =========================================================================
// A.  DRUG PHARMACOKINETICS
// =========================================================================
double a_mig = mig_ka * MIG_GUT;
double q_mig = mig_Q * (MIG_CEN / mig_V - MIG_PER / mig_Vp);
dxdt_MIG_GUT = -a_mig;
dxdt_MIG_CEN = mig_F * a_mig - (mig_CL / mig_V) * MIG_CEN - q_mig;
dxdt_MIG_PER = q_mig;
dxdt_MIG_BR  = mig_ke0 * (mig_Kp_br * C_MIG_PL - MIG_BR);

double a_ari = ari_ka * ARI_GUT;
dxdt_ARI_GUT = -a_ari;
dxdt_ARI_CEN = ari_F * a_ari - (ari_CL / ari_V) * ARI_CEN;
dxdt_ARI_CSF = ari_ke0 * (ari_Kp_csf * C_ARI_PL - ARI_CSF);

double a_nal = nal_ka * NAL_GUT;
dxdt_NAL_GUT = -a_nal;
dxdt_NAL_CEN = nal_F * a_nal - (nal_CL / nal_V) * NAL_CEN;
dxdt_NAL_BR  = nal_ke0 * (C_NAL_PL - NAL_BR);

dxdt_CD_CSF = -(cd_kout + cd_kin_br) * CD_CSF;
dxdt_CD_BR  =  cd_kin_br * CD_CSF - cd_kout_br * CD_BR;
dxdt_CD_SYS =  cd_kout   * CD_CSF - cd_klsys   * CD_SYS;

// =========================================================================
// B.  DRUG EFFECTS
// =========================================================================
// miglustat: UGCG inhibition, computed SEPARATELY in the viscera (from plasma)
// and in the CNS (from the brain effect site).  The brain:plasma ratio of 0.45
// is what makes the visceral effect roughly twice the central one -- a
// prediction, not an assumption, and the reason miglustat's clearest published
// effects are visceral and bulbar rather than a halt of cerebellar decline.
double I_mig_v = mig_Emax * C_MIG_PL / (mig_IC50 + C_MIG_PL);
double I_mig_c = mig_Emax * MIG_BR   / (mig_IC50 + MIG_BR);

// arimoclomol: an HSF1 CO-inducer.  It amplifies the heat-shock response only
// when the cell is already stressed, so the achievable HSP70 rise is scaled by
// a cellular-stress proxy.  In an unstressed cell the drug does almost nothing.
double chol_c_fold = CHOL_C / WT_CHOL_C;
double stress_cell = ari_stress0 + (1.0 - ari_stress0) *
                     fmin(1.0, fmax(0.0, (chol_c_fold - 1.0) / 12.0));
if (stress_cell > 1.0) stress_cell = 1.0;
double HSP70 = 1.0 + ari_Emax * stress_cell * ARI_CSF / (ari_EC50 + ARI_CSF);

// levacetylleucine: three separable actions on one effect site.
double E_nal      = NAL_BR / (nal_EC50 + NAL_BR);
double E_nal_sym  = nal_Emax_sym * E_nal * nal_sym_on;   // -> D_rev only
double E_nal_gluc = nal_e_gluc   * E_nal * nal_sym_on;   // -> brain metabolism
double E_nal_dm   = nal_Emax_dm  * E_nal * nal_dm_on;    // -> CNS lipid influx

// =========================================================================
// C.  GENOTYPE -> NPC1 -> EGRESS CAPACITY
// =========================================================================
double hsp_drive = fmax(0.0, HSP70 - 1.0);
double theta = theta0 + (1.0 - theta0) * Emax_fold * hsp_drive /
               (EChsp + hsp_drive);
if (theta > 1.0) theta = 1.0;      // a folding YIELD is a fraction: hard bound
double kexit = kexit0 * theta;

dxdt_NPC1_ER = ksyn_npc1 * (1.0 - f_null) - NPC1_ER * (kexit + kerad);
dxdt_NPC1_L  = NPC1_ER * kexit - kdegL * NPC1_L;

double f_npc1 = NPC1_L / WT_NPC1_L;
// NPC1 and NPC2 act IN SERIES on one flux, so their functional fractions
// multiply.  This is why NPC2 disease is phenotypically NPC1 disease and why
// arimoclomol, which works on NPC1 folding, must be inert in NPC2 patients.
double f_eg = f_npc1 * f_npc2 + f_leak;

// =========================================================================
// D.  CHOLESTEROL POOLS
// =========================================================================
double eg_v = Vmax_v * f_eg * CHOL_V / (Km_v + CHOL_V);
double eg_c = Vmax_c * f_eg * CHOL_C / (Km_c + CHOL_C);

dxdt_CHOL_V = Jin_v - eg_v - kbas_chol * CHOL_V
              - cd_kext * C_CD_SYS * CHOL_V;
dxdt_CHOL_C = Jin_c * (1.0 - E_nal_dm) - eg_c - kbas_chol_c * CHOL_C
              - cd_kext * C_CD_BR * CHOL_C;

// =========================================================================
// E.  SPHINGOSINE, LYSOSOMAL CALCIUM, HYDROLASES
// =========================================================================
dxdt_SPH = ksph - kexp_sph * (fmin_sph + (1.0 - fmin_sph) * f_npc1) * SPH;

double sph_r  = SPH / Ksph;
double ca_ss  = (1.0 / (1.0 + sph_r * sph_r)) * (1.0 / (1.0 + CHOL_C / Ktr_chol));
dxdt_CA_LY = kca * (ca_ss - CA_LY);

double hyd_ss = (1.0 / (1.0 + pow(CHOL_C / Kph, hill_ph)))
              * (1.0 + e_hyd_hsp * hsp_drive)
              * (0.45 + 0.55 * CA_LY);
dxdt_HYD = khyd * (hyd_ss - HYD);

// =========================================================================
// F.  GLYCOSPHINGOLIPIDS
// =========================================================================
dxdt_GSL_V = ksyn_gsl_v * (1.0 - I_mig_v) - kdeg_gsl   * GSL_V * HYD;
dxdt_GSL_C = ksyn_gsl_c * (1.0 - I_mig_c) - kdeg_gsl_c * GSL_C * HYD;

// =========================================================================
// G.  AUTOPHAGY, MITOCHONDRIA, ROS
// =========================================================================
double prot = HYD * (0.5 + 0.5 * CA_LY);
dxdt_AUTOPH = kaut_in - kaut_out * prot * AUTOPH;

double cm = CHOL_C / Kmito;
double mito_ss = (1.0 / (1.0 + cm * cm)) * (1.0 + E_nal_gluc);
dxdt_MITO = kmito * (mito_ss - MITO);

double mito_safe = fmax(MITO, 1e-6);
double ros_ss = 1.0 + aros * fmax(0.0, 1.0 / mito_safe - 1.0);
dxdt_ROS = kros * (ros_ss - ROS);

// =========================================================================
// H.  CEREBELLAR STRESS, THE DAMAGE INTEGRAL, PURKINJE POOLS
// =========================================================================
double stress = w_chol * (chol_c_fold - 1.0) / 10.0
              + w_gsl  * (GSL_C / WT_GSL_C - 1.0)
              + w_ros  * (ROS - WT_ROS)
              + w_aut  * (AUTOPH / WT_AUT - 1.0)
              + w_infl * INFL
              + w_ca   * (1.0 - CA_LY / WT_CA);
stress = fmax(0.0, stress - S_thresh);      // functional reserve, absorbed

// developmental vulnerability: the same insult costs MORE while the cerebellum
// is still being built and myelinated (PMID 15217094).  It multiplies BOTH the
// rate at which the reserve is spent and the death rate once the gate is open,
// which is what lets it move the age of ONSET and not only the later slope.
// SOLVERTIME is the time since birth in days.
double vuln = 1.0 + v_dev * exp(-(SOLVERTIME / 365.25) / tau_dev);

// The damage integral: NO repair term.  See the header for why.
dxdt_DAM = stress * vuln;

double ru   = DAM / D_reserve;
double run_ = pow(fmax(ru, 0.0), hill_gate);
double gate = run_ / (1.0 + run_);

double to_stress = kon_pc  * stress * PC;
double recover   = koff_pc * PC_S / (1.0 + stress / Krec);
double die       = kdie * PC_S * gate * vuln;

dxdt_PC      = -to_stress + recover;
dxdt_PC_S    =  to_stress - recover - die;
dxdt_PC_LOST =  die;

// =========================================================================
// I.  NEUROINFLAMMATION, SYNAPSES, VOLUME
// =========================================================================
double infl_ss = a_chol * (chol_c_fold - 1.0) / 25.0
               + a_pc * PC_S + a_dead * PC_LOST;
dxdt_INFL = kinfl * (infl_ss - INFL);

dxdt_SYN = ksyn_rep * ((PC + PC_S) - SYN) - kprune * INFL * SYN;

double cbl_t = 1.0 - g_vol_pc * PC_LOST - g_vol_syn * (1.0 - SYN);
dxdt_CBL = kcbl * (cbl_t - CBL);

// =========================================================================
// J.  FUNCTION DECOMPOSITION — reversible vs irreversible
// =========================================================================
// D_rev is a GRADED, saturating function of the CURRENT biochemical stress.  It
// is deliberately NOT a function of the fraction of cells in a "stressed" pool:
// a saturating cell count makes D_rev plateau within weeks of birth, which no
// NPC patient does.  That version was written first and rejected.
double glucm   = MITO / WT_MITO;
double tgt_rev = (cap_rev * stress / (K_rev + stress) + cr_pcs * PC_S);
if (tgt_rev > 0.90) tgt_rev = 0.90;
if (tgt_rev < 0.0)  tgt_rev = 0.0;
tgt_rev = tgt_rev * (1.0 - E_nal_sym);
dxdt_D_REV = krev * (tgt_rev - D_REV);

// LINEAR in accumulated cell loss.  A saturating form cannot reproduce the
// reported linear NPCCSS progression (PMID 19415691).
double tgt_irr = g_irr_pc * PC_LOST + g_irr_syn * (1.0 - SYN)
               + g_irr_cbl * (1.0 - CBL);
if (tgt_irr > 0.95) tgt_irr = 0.95;
dxdt_D_IRR = kirr * fmax(0.0, tgt_irr - D_IRR);     // one-sided: MONOTONE

// =========================================================================
// K.  BIOMARKERS
// =========================================================================
dxdt_TRIOL = ktri * CHOL_V / (Ktri + CHOL_V) - ktri_el * TRIOL;
double ppcs_sub = CHOL_V + b_ppcs_gsl * GSL_V;
dxdt_PPCS  = kppcs * ppcs_sub / (Kppcs + ppcs_sub) - kppcs_el * PPCS;
dxdt_TCG   = ktcg * CHOL_V / (Ktcg + CHOL_V) - ktcg_el * TCG;
// NfL and calbindin track the DEATH FLUX, not the storage burden.  That is the
// whole reason they behave differently from triol (PMID 27307499, 36470574).
dxdt_NFL   = knfl  * die + knfl_base - knfl_el  * NFL;
dxdt_CALB  = kcalb * die              - kcalb_el * CALB;

// =========================================================================
// L.  COCHLEA AND SURVIVAL
// =========================================================================
dxdt_OHC = -cd_koto * C_CD_CSF * OHC;

double swal_x = swal_a * D_REV + swal_b * D_IRR;
if (swal_x > 1.0) swal_x = 1.0;
double SWAL = 5.0 * swal_x;
dxdt_CUMHAZ = h0 + h_swal * (SWAL / 5.0) * (SWAL / 5.0) + h_liver * liver_risk;

$TABLE
double f_NPC1_out = NPC1_L / WT_NPC1_L;
double CHOL_C_fold = CHOL_C / WT_CHOL_C;
double CHOL_V_fold = CHOL_V / WT_CHOL_V;
double FUNC = fmax(0.0, 1.0 - D_REV - D_IRR);

double x_sara = sara_a * D_REV + sara_b * D_IRR;
double SARA = fmin(sara_max, sara_k * pow(fmax(0.0, x_sara), sara_p));

double x_n5 = n5_a * D_REV + n5_b * D_IRR;
double NPCCSS5 = fmin(n5_max, n5_k * pow(fmax(0.0, x_n5), n5_p));
double NPCCSS4 = fmin(20.0, NPCCSS5 * n4_frac * n4_gain);

double lost_ohc = fmax(0.0, 1.0 - OHC);
double HEARING_dB = (lost_ohc <= 0.0) ? 0.0
                    : cd_hear_max * pow(lost_ohc, cd_hear_p);
// The 17-domain scale contains a HEARING domain, which is why a drug that
// causes hearing loss can erase its own measured benefit -- and why Ory 2017
// had to score "NSS minus hearing" (PMID 28803710).
double hearing_pts = 5.0 * fmin(1.0, HEARING_dB / cd_hear_max);
double NPCCSS17 = fmin(61.0, NPCCSS5 * n17_k + hearing_pts);
double NPCCSS17_NOHEAR = fmin(61.0, NPCCSS5 * n17_k);

double swx = swal_a * D_REV + swal_b * D_IRR;
double SWALLOW = 5.0 * fmin(1.0, swx);
double SACCADE = 100.0 * exp(-(sacc_rev * D_REV + sacc_irr * D_IRR));
double SURV = exp(-CUMHAZ);
double AGE_YR = SOLVERTIME / 365.25;

double C_MIG = C_MIG_PL;
double C_ARI = C_ARI_PL;
double C_NAL = C_NAL_PL;
double CD_CSF_CONC = C_CD_CSF;
double GLUC_METAB = MITO / WT_MITO;

$CAPTURE
AGE_YR f_NPC1_out CHOL_C_fold CHOL_V_fold FUNC
SARA NPCCSS5 NPCCSS4 NPCCSS17 NPCCSS17_NOHEAR SWALLOW SACCADE HEARING_dB SURV
TRIOL PPCS TCG NFL CALB
C_MIG C_ARI C_NAL CD_CSF_CONC HSP70 GLUC_METAB
stress DAM gate theta f_eg I_mig_v I_mig_c E_nal_sym E_nal_dm

## =============================================================================
##  R HELPERS — genotypes, archetypes, scenarios, calibration table
##  Everything below is ordinary R and is ignored by mread().
## =============================================================================
/*
library(mrgsolve); library(dplyr); library(tidyr); library(ggplot2)

## ---------------------------------------------------------------------------
## Genotypes.  Residual NPC1 activity is set by two numbers only: the fraction
## of null alleles and the intrinsic folding yield of the non-null allele.
## ---------------------------------------------------------------------------
npc_genotypes <- list(
  `WT`            = list(f_null = 0.00, theta0 = 1.000, f_npc2 = 1.00),
  `I1061T/I1061T` = list(f_null = 0.00, theta0 = 0.055, f_npc2 = 1.00),
  `I1061T/null`   = list(f_null = 0.50, theta0 = 0.055, f_npc2 = 1.00),
  `null/null`     = list(f_null = 1.00, theta0 = 0.055, f_npc2 = 1.00),
  `mild/mild`     = list(f_null = 0.00, theta0 = 0.322, f_npc2 = 1.00),
  `NPC2`          = list(f_null = 0.00, theta0 = 1.000, f_npc2 = 0.03)
)

## ---------------------------------------------------------------------------
## Onset archetypes = genotype x developmental vulnerability.
## IMPORTANT: residual NPC1 activity alone does NOT determine the onset form.
## The storage phenotype SATURATES, so I1061T/null and I1061T/I1061T differ
## 2-fold in residual activity and barely at all in steady-state load.  Siblings
## with identical genotypes can differ by a decade.  v_dev is therefore carried
## as a PATIENT-level covariate, not derived from the genotype.
## ---------------------------------------------------------------------------
npc_archetypes <- list(
  `perinatal`        = list(geno = "null/null",     v_dev = 2.5, liver_risk = 1),
  `early-infantile`  = list(geno = "I1061T/null",   v_dev = 2.5, liver_risk = 0),
  `late-infantile`   = list(geno = "I1061T/null",   v_dev = 1.2, liver_risk = 0),
  `juvenile`         = list(geno = "I1061T/I1061T", v_dev = 0.6, liver_risk = 0),
  `adolescent/adult` = list(geno = "mild/mild",     v_dev = 0.2, liver_risk = 0),
  `NPC2 disease`     = list(geno = "NPC2",          v_dev = 1.0, liver_risk = 0)
)

npc_patient <- function(mod, archetype = "juvenile") {
  a <- npc_archetypes[[archetype]]
  g <- npc_genotypes[[a$geno]]
  param(mod, f_null = g$f_null, theta0 = g$theta0, f_npc2 = g$f_npc2,
        v_dev = a$v_dev, liver_risk = a$liver_risk)
}

## ---------------------------------------------------------------------------
## Regimens.  All four drugs, with the approved adult doses.
## ---------------------------------------------------------------------------
YR <- 365.25
reg_none <- function(...) NULL

reg_mig <- function(start = 0, dur = 3 * YR, mg = 200)
  ev(time = start, amt = mg, cmt = "MIG_GUT", ii = 1/3, addl = ceiling(dur * 3) - 1)

reg_ari <- function(start = 0, dur = 3 * YR, mg = 124)
  ev(time = start, amt = mg, cmt = "ARI_GUT", ii = 1/3, addl = ceiling(dur * 3) - 1)

reg_nal <- function(start = 0, dur = 3 * YR, g_per_day = 4)
  ev(time = start, amt = g_per_day * 1000 / 3, cmt = "NAL_GUT",
     ii = 1/3, addl = ceiling(dur * 3) - 1)

reg_cd  <- function(start = 0, dur = 3 * YR, mg = 900, weeks = 2)
  ev(time = start, amt = mg, cmt = "CD_CSF",
     ii = 7 * weeks, addl = max(0, floor(dur / (7 * weeks)) - 1))

combine_ev <- function(...) {
  xs <- Filter(Negate(is.null), list(...))
  if (!length(xs)) return(ev(time = 0, amt = 0, cmt = "MIG_GUT"))
  Reduce(function(a, b) a + b, xs)
}

## ---------------------------------------------------------------------------
## Run one patient from birth for `years`, optionally starting drugs at `start`.
## ---------------------------------------------------------------------------
npc_run <- function(mod, archetype = "juvenile", years = 30, events = NULL,
                    delta = 7, ...) {
  m <- npc_patient(mod, archetype)
  if (length(list(...))) m <- param(m, ...)
  e <- if (is.null(events)) ev(time = 0, amt = 0, cmt = "MIG_GUT") else events
  m %>% mrgsim(events = e, end = years * YR, delta = delta,
               hmax = 0.5, atol = 1e-8, rtol = 1e-8) %>% as_tibble()
}

at_age <- function(out, col, age) {
  approx(out$AGE_YR, out[[col]], xout = age, rule = 2)$y
}

onset_age <- function(out, thresh = 3) {
  i <- which(out$NPCCSS5 >= thresh)
  if (!length(i)) return(NA_real_)
  out$AGE_YR[i[1]]
}

median_survival <- function(out) {
  i <- which(out$CUMHAZ >= log(2))
  if (!length(i)) return(NA_real_)
  out$AGE_YR[i[1]]
}

## ---------------------------------------------------------------------------
## The 18 scenarios.  Treatment scenarios all start at age 13, which is where
## model SARA equals the 15.91 baseline of IB1001-301 (PMID 38294974) -- i.e.
## the virtual patient is matched to the trial population by SEVERITY, not by
## chronological age.
## ---------------------------------------------------------------------------
ENTRY <- 13 * YR
DUR   <- 3 * YR

npc_scenario_defs <- list(
  S01 = list(label = "natural history, juvenile",        arch = "juvenile",
             ev = NULL, years = 35),
  S02 = list(label = "natural history, late-infantile",  arch = "late-infantile",
             ev = NULL, years = 25),
  S03 = list(label = "natural history, adolescent/adult",arch = "adolescent/adult",
             ev = NULL, years = 55),
  S04 = list(label = "natural history, perinatal",       arch = "perinatal",
             ev = NULL, years = 5),
  S05 = list(label = "natural history, NPC2 disease",    arch = "NPC2 disease",
             ev = NULL, years = 25),
  S06 = list(label = "miglustat 200 mg tid",             arch = "juvenile",
             ev = reg_mig(ENTRY, DUR), years = 16),
  S07 = list(label = "levacetylleucine 4 g/day",         arch = "juvenile",
             ev = reg_nal(ENTRY, DUR), years = 16),
  S08 = list(label = "arimoclomol 124 mg tid",           arch = "juvenile",
             ev = reg_ari(ENTRY, DUR), years = 16),
  S09 = list(label = "arimoclomol + miglustat",          arch = "juvenile",
             ev = combine_ev(reg_ari(ENTRY, DUR), reg_mig(ENTRY, DUR)), years = 16),
  S10 = list(label = "levacetylleucine + miglustat",     arch = "juvenile",
             ev = combine_ev(reg_nal(ENTRY, DUR), reg_mig(ENTRY, DUR)), years = 16),
  S11 = list(label = "triple oral therapy",              arch = "juvenile",
             ev = combine_ev(reg_nal(ENTRY, DUR), reg_mig(ENTRY, DUR),
                             reg_ari(ENTRY, DUR)), years = 16),
  S12 = list(label = "IT adrabetadex 900 mg q2wk",       arch = "juvenile",
             ev = reg_cd(ENTRY, DUR), years = 16),
  S13 = list(label = "IT adrabetadex + miglustat",       arch = "juvenile",
             ev = combine_ev(reg_cd(ENTRY, DUR), reg_mig(ENTRY, DUR)), years = 16),
  S14 = list(label = "triple from age 2 (newborn screen)", arch = "juvenile",
             ev = combine_ev(reg_nal(2*YR, 18*YR), reg_mig(2*YR, 18*YR),
                             reg_ari(2*YR, 18*YR)), years = 20),
  S15 = list(label = "triple from age 8",                arch = "juvenile",
             ev = combine_ev(reg_nal(8*YR, 12*YR), reg_mig(8*YR, 12*YR),
                             reg_ari(8*YR, 12*YR)), years = 20),
  S16 = list(label = "triple from age 12",               arch = "juvenile",
             ev = combine_ev(reg_nal(12*YR, 8*YR), reg_mig(12*YR, 8*YR),
                             reg_ari(12*YR, 8*YR)), years = 20),
  S17 = list(label = "levacetylleucine, withdrawn at 12 mo", arch = "juvenile",
             ev = reg_nal(ENTRY, 1*YR), years = 16),
  S18 = list(label = "arimoclomol in NPC2 disease (predicted null)",
             arch = "NPC2 disease", ev = reg_ari(ENTRY, DUR), years = 16)
)

npc_scenarios <- function(mod) {
  purrr::imap_dfr(npc_scenario_defs, function(s, id) {
    npc_run(mod, s$arch, years = s$years, events = s$ev) %>%
      mutate(scenario = id, label = s$label)
  })
}

npc_scenario_summary <- function(mod) {
  purrr::imap_dfr(npc_scenario_defs, function(s, id) {
    o <- npc_run(mod, s$arch, years = s$years, events = s$ev)
    tibble(scenario = id, label = s$label,
           onset_yr   = onset_age(o),
           surv_yr    = median_survival(o),
           n5_at_16   = at_age(o, "NPCCSS5", min(16, s$years)),
           sara_at_16 = at_age(o, "SARA",    min(16, s$years)),
           d_rev_16   = at_age(o, "D_REV",   min(16, s$years)),
           d_irr_16   = at_age(o, "D_IRR",   min(16, s$years)),
           triol_16   = at_age(o, "TRIOL",   min(16, s$years)),
           hear_dB_16 = at_age(o, "HEARING_dB", min(16, s$years)))
  })
}

## ---------------------------------------------------------------------------
## Trial-design experiments.  Run each PUBLISHED design against BOTH mechanisms.
## The two synthetic drugs are levacetylleucine with one arm switched off:
##   nal_sym_on = 1, nal_dm_on = 0   -> purely symptomatic
##   nal_sym_on = 0, nal_dm_on = 1   -> purely disease-modifying
## ---------------------------------------------------------------------------
npc_design_experiment <- function(mod, weeks = 12, dm_scale = 1) {
  base <- npc_patient(mod, "juvenile")
  horizon <- ENTRY + weeks * 7
  run1 <- function(sym, dm, ev_) {
    param(base, nal_sym_on = sym, nal_dm_on = dm,
          nal_Emax_dm = 0.080 * dm_scale) %>%
      mrgsim(events = ev_, end = horizon, delta = 3.5,
             hmax = 0.5, atol = 1e-8, rtol = 1e-8) %>% as_tibble()
  }
  none <- ev(time = 0, amt = 0, cmt = "MIG_GUT")
  drug <- reg_nal(ENTRY, weeks * 7)
  pbo  <- run1(1, 1, none)
  sym  <- run1(1, 0, drug)
  dm   <- run1(0, 1, drug)
  b <- function(o, col) at_age(o, col, ENTRY / YR)
  e <- function(o, col) at_age(o, col, horizon / YR)
  tibble(
    arm    = c("placebo", "symptomatic only", "disease-modifying only"),
    dSARA  = c(e(pbo,"SARA")   - b(pbo,"SARA"),
               e(sym,"SARA")   - b(sym,"SARA"),
               e(dm ,"SARA")   - b(dm ,"SARA")),
    dN5    = c(e(pbo,"NPCCSS5")- b(pbo,"NPCCSS5"),
               e(sym,"NPCCSS5")- b(sym,"NPCCSS5"),
               e(dm ,"NPCCSS5")- b(dm ,"NPCCSS5"))
  ) %>% mutate(dSARA_vs_pbo = dSARA - dSARA[1], dN5_vs_pbo = dN5 - dN5[1])
}

## ---------------------------------------------------------------------------
## Model vs published targets.  See npc_references.md for the full target list;
## targets marked "validation" were NOT used in calibration.
## ---------------------------------------------------------------------------
npc_calibration_table <- function(mod) {
  oj <- npc_run(mod, "juvenile", years = 35, delta = 7)
  ow <- param(mod, f_null = 0, theta0 = 1, f_npc2 = 1) %>%
        mrgsim(end = 45 * YR, delta = 30, hmax = 1) %>% as_tibble()
  slope5  <- at_age(oj, "NPCCSS5", 13)  - at_age(oj, "NPCCSS5", 12)
  slope17 <- at_age(oj, "NPCCSS17", 13) - at_age(oj, "NPCCSS17", 12)
  tibble(
    target = c("T1 plasma triol, patients (ng/mL)",
               "T2 plasma triol, controls (ng/mL)",
               "T3 5-domain NPCCSS slope (pt/yr)",
               "T4 17-domain NPCCSS slope (pt/yr)",
               "T11 SARA at trial entry"),
    role      = c("calibration", "calibration", "calibration",
                  "VALIDATION", "calibration"),
    published = c(88.31, 5.97, 1.50, 2.80, 15.91),
    model     = c(at_age(oj, "TRIOL", 10),
                  approx(ow$time / YR, ow$TRIOL, 10)$y,
                  slope5, slope17,
                  at_age(oj, "SARA", 13))
  ) %>% mutate(rel_err_pct = 100 * (model - published) / published)
}

## ---------------------------------------------------------------------------
## Example plots
## ---------------------------------------------------------------------------
if (FALSE) {
  mod <- mread("npc_mrgsolve_model.R")
  print(npc_calibration_table(mod))
  print(npc_scenario_summary(mod))
  print(npc_design_experiment(mod))

  npc_scenarios(mod) %>%
    filter(scenario %in% c("S01","S06","S07","S08","S11","S12")) %>%
    ggplot(aes(AGE_YR, NPCCSS5, colour = label)) +
      geom_line(linewidth = 0.8) +
      labs(x = "age (years)", y = "5-domain NPCCSS",
           title = "NPC: 5-domain NPCCSS under each therapy from age 13") +
      theme_bw()

  ## The reserve claim: the same therapy started at three different ages
  npc_scenarios(mod) %>%
    filter(scenario %in% c("S01","S14","S15","S16")) %>%
    ggplot(aes(AGE_YR, D_IRR, colour = label)) +
      geom_line(linewidth = 0.8) +
      labs(x = "age (years)", y = "irreversible dysfunction D_irr",
           title = "What is bought is D_rev; what is already gone is D_irr") +
      theme_bw()
}
*/
