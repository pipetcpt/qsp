## =============================================================================
## Aneurysmal Subarachnoid Haemorrhage -> Delayed Cerebral Ischaemia (aSAH-DCI)
## 동맥류성 지주막하출혈 후 지연성 뇌허혈 — mrgsolve QSP 모델
##
## 39 ODEs (8 drug PK/exposure + 31 disease / haemodynamic / injury states) ·
## 15 treatment scenarios · 12 analysis functions.  Time unit = DAYS.
##
## -----------------------------------------------------------------------------
## WHAT THIS MODEL ASSERTS
## -----------------------------------------------------------------------------
## Delayed cerebral ischaemia is almost always described as a chain:
##
##      cisternal blood -> angiographic vasospasm -> ischaemia -> infarct
##
## That chain cannot explain the four facts that dominate the aSAH literature:
##
##   (i)   clazosentan cuts moderate-severe angiographic vasospasm by ~65%
##         (CONSCIOUS-1) and moves 3-month outcome by nothing (CONSCIOUS-2/3);
##   (ii)  nimodipine barely changes angiographic calibre and is the ONLY drug
##         that improves outcome;
##   (iii) a large minority of DCI happens with no angiographic vasospasm, and
##         most patients with severe vasospasm never develop DCI;
##   (iv)  induced hypertension and milrinone help some patients and not others,
##         with no consistent trial signal (HIMALAIA).
##
## This model answers all four with THREE structural claims.
##
## CLAIM 1 — ONE SHARED ARTERIOLAR RESERVE, FOUR COMPETING CONSUMERS.
##   Cerebral vascular resistance is split into three series elements:
##       R_large  conducting arteries, of which only the SPASTIC SEGMENT
##                (FSEG = 0.33 of R_large) obeys Poiseuille r^-4;
##       RA       arterioles — the ONLY adjustable element, range
##                RA0 0.85 down to RAMIN 0.20, i.e. a reserve of 0.65;
##       R_micro  microvessels and capillaries.
##   Autoregulation sets RA to whatever holds CBF at target:
##       RA_desired = CPP/CBF_target - R_large - R_micro
##       RA         = AREG*clamp(RA_desired) + (1-AREG)*RA0
##       CBF        = CPP / (R_large + RA + R_micro)
##   Four things consume that ONE reserve: large-artery spasm, microvascular
##   tone, capillary microthrombosis, and CPP loss — plus spreading
##   depolarisation, which CONSTRICTS the reserve element itself (inverse
##   neurovascular coupling, the KINV term).
##
##   Consequence, and the whole point: the mapping from any single consumer to
##   the endpoint is nonlinear and CONTEXT-DEPENDENT.  The same 33% calibre loss
##   is harmless in a patient with an untouched reserve and infarcting in one
##   whose reserve is already spent on microthrombosis or whose autoregulation
##   is gone.  Blocking one consumer helps only the subpopulation in which that
##   consumer happened to be the MARGINAL one.  This is why a drug can win the
##   angiogram and lose the trial.
##
## CLAIM 2 — THE MICROCIRCULATION IS NOT A RESISTANCE FOOTNOTE; IT DEGRADES
##           OXYGEN EXTRACTION.
##   Following Ostergaard/Jespersen, capillary transit-time heterogeneity
##       CTH      = (RMIC - 1) + 2.5*MTHR
##       OEFmax   = 0.85 / (1 + KCTH*CTH)
##   so the maximum extractable oxygen fraction FALLS with capillary
##   dysfunction.  Ischaemia is therefore
##       ISCH = 1 - DO2*OEFmax / CMRO2,      DO2 = CBF * CaO2
##   and the same CBF that is adequate in a homogeneous bed is not adequate in a
##   plugged one.  This is what makes microthrombosis and pericyte tone
##   first-class consumers, and it is why anaemia (CaO2) and fever (CMRO2) act
##   on exactly the same node as vasospasm.
##
## CLAIM 3 — THE CLOCK IS HAEMOGLOBIN HANDLING, NOT A CALENDAR.
##   There is NO day-4 switch anywhere in this model.  The window emerges from
##   a transit chain plus an inducible sink:
##       CLOT --k_mat--> RBCL --k_hem--> OXYHB --(HP, CD163, HO-1, bulk CSF)-->
##   Erythrocytes must first enter haemolysis (that transit is the delay), and
##   the dominant clearance term (HO-1) is INDUCIBLE with t1/2 1.9 d (that is
##   the lag).  Free oxyHb then drives ET-1, scavenges NO, and feeds ROS.
##   Because the input is a transit and the sink is inducible, the effector
##   curve peaks in the second half of the first week and resolves by ~day 14.
##   Two predictions follow directly: heavier clot moves the peak EARLIER and
##   HIGHER; and CSF drainage works by TRUNCATING THE INPUT FUNCTION, which is
##   structurally different from — and not interchangeable with — blocking an
##   effector downstream.
##
## Together the three claims say the therapeutic levers are NOT substitutes:
##   input-side levers (drainage, irrigation) shrink every consumer at once;
##   effector-side levers (clazosentan, nicardipine, fasudil) shrink ONE;
##   reserve-side levers (induced hypertension, milrinone) only help when the
##   reserve is pressure-passive (AREG low); demand/content-side levers
##   (transfusion, normothermia, SD suppression) never touch resistance at all.
##   And every lever has a harm channel into the SAME endpoint.
##
## -----------------------------------------------------------------------------
## TWELVE COMPUTED RESULTS (functions, not prose)
## -----------------------------------------------------------------------------
##  1  sah_scenarios()       15 arms: incidence, calibre, DCI, infarct, outcome
##  2  sah_risk_ratios()     RR vs standard of care, with trial comparators
##  3  sah_dissociation()    r^2 between angiographic calibre and infarct;
##                           DCI without spasm; spasm without DCI; TCD vs PbtO2
##                           as diagnostic tests
##  4  sah_reserve_decomp()  which consumer is eating the reserve, by day
##  5  sah_window()          the emergent time window by mFisher and Hp genotype
##  6  sah_knockouts()       perfect single-consumer blocks and their pairs;
##                           the superadditivity that single-target trials miss
##  7  sah_pressure_test()   induced hypertension / milrinone stratified by
##                           autoregulatory state
##  8  sah_clazosentan()     where the angiographic win is lost: redundancy
##                           first, harm channel second
##  9  sah_monitoring()      lead times of TCD, PbtO2 and infarct onset
## 10  sah_systemic()        sodium, ICP, SD burden, autoregulation nadirs
## 11  sah_anaemia()         CaO2 as a consumer of the same node; transfusion
## 12  sah_route()           oral vs IV nimodipine — same exposure, different
##                           CPP, opposite conclusion
##
## -----------------------------------------------------------------------------
## CALIBRATION NOTES (targets from the trial literature)
## -----------------------------------------------------------------------------
##   moderate-severe angiographic vasospasm, placebo   66%   (CONSCIOUS-1)
##   DCI incidence, standard care                      ~30%  (pooled)
##   poor outcome mRS 4-6 at 90 d                      30-35%
##   hyponatraemia < 135 mmol/L                        30-50%
##   clazosentan RR moderate-severe vasospasm          0.35  (CONSCIOUS-1)
##   clazosentan RR poor outcome                       ~1.05 (CONSCIOUS-2)
##   nimodipine  RR poor outcome                       0.67  (Cochrane)
##   cilostazol  RR DCI                                ~0.47 (meta-analysis)
##   simvastatin RR poor outcome                       ~1.00 (STASH, null)
##   lumbar drainage RR unfavourable outcome           0.76  (EARLYDRAIN)
##
##   Every number quoted in README.md is produced by sah_reference_check.py,
##   an independent vectorised numpy/RK4 transcription of THIS model, whose
##   full output is committed as sah_reference_check_output.txt.  Two places
##   where the model is knowingly shallower than the literature are recorded
##   there and in README.md rather than tuned away:
##     - nimodipine RR poor outcome 0.79 vs Cochrane 0.67
##     - DCI without moderate-severe angiographic spasm 10% vs 20-30% reported
##
## -----------------------------------------------------------------------------
## Requires: mrgsolve, dplyr, tidyr, ggplot2
## =============================================================================

library(mrgsolve)
library(dplyr)
library(tidyr)

sah_code <- '
$PARAM @annotated
// ---------------- haemorrhage clearance -------------------------------------
KMAT    : 0.250  : clot -> erythrocytes entering haemolysis (1/d)
KMATV   : 0.185  : same for intraventricular compartment (1/d)
KHEM    : 0.230  : haemolysis, RBC pool -> free oxyHb (1/d)
KDRAIN  : 0.036  : extra physical clearance per unit drainage intensity (1/d)
KREL    : 9.20   : oxyHb released per unit haemolysed (RU)
KHPCL   : 1.40   : haptoglobin-mediated oxyHb scavenging (1/d)
KMHP    : 0.45   : Km for haptoglobin saturation (RU)
KMAC    : 0.55   : CD163 / macrophage uptake, baseline (1/d)
KCSF    : 0.22   : bulk CSF clearance of free oxyHb (1/d)
KHPIN   : 0.60   : haptoglobin synthesis (1/d)
KHPOUT  : 0.60   : haptoglobin turnover (1/d)
KHPUSE  : 0.85   : haptoglobin consumption by complex formation (1/d)
KHEME   : 0.85   : oxyHb -> free haem (1/d)
KHOX    : 0.95   : HO-1 dependent haem degradation (1/d)
KHO1IN  : 1.30   : HO-1 induction rate (1/d)
KHO1OUT : 0.36   : HO-1 turnover, t1/2 1.9 d = THE LAG (1/d)
KBOXIN  : 0.55   : bilirubin oxidation product formation (1/d)
KBOXOUT : 0.42   : BOX elimination (1/d)
// ---------------- effectors --------------------------------------------------
KETIN   : 0.80   : endothelin-1 formation (1/d)
KETOUT  : 0.80   : endothelin-1 clearance (1/d)
TAUNO   : 0.30   : NO bioavailability time constant (d)
KNOHB   : 0.85   : oxyHb weight in NO scavenging (1/RU)
KNOROS  : 0.35   : ROS weight in NO loss (1/RU)
KROSIN  : 0.60   : ROS formation (1/d)
KROSOUT : 1.05   : ROS scavenging (1/d)
TAURHO  : 1.00   : Rho-kinase Ca-sensitisation time constant (d)
KINFIN  : 0.40   : neuroinflammation formation (1/d)
KINFOUT : 0.55   : neuroinflammation resolution (1/d)
KPAIIN  : 0.45   : PAI-1 / prothrombotic tone formation (1/d)
KPAIOUT : 0.60   : PAI-1 turnover (1/d)
// ---------------- vascular ---------------------------------------------------
KSPIN   : 0.550  : large-artery constriction on-rate (1/d)
KSPOUT  : 0.620  : large-artery relaxation rate (1/d)
SPMAX   : 0.820  : maximum attainable constriction fraction (-)
KSTIN   : 0.045  : structural remodelling formation (1/d)
KSTOUT  : 0.045  : structural remodelling regression (1/d)
TAURM   : 0.500  : microvascular tone time constant (d)
KRMIC   : 1.150  : microvascular tone gain (-)
KMTIN   : 0.0480 : microthrombosis formation (1/d)
KMTOUT  : 0.1300 : microthrombus lysis (1/d)
KAREGR  : 0.160  : autoregulation recovery (1/d)
KAREGL  : 0.350  : autoregulation loss (1/d)
// ---------------- spreading depolarisation -----------------------------------
TAUSDS  : 0.40   : SD susceptibility time constant (d)
SDMAX   : 7.00   : maximum SD rate (events/d)
KSDB    : 0.60   : SD burden decay (1/d)
KINV    : 0.55   : SD inverse-coupling gain on arteriolar resistance (-)
// ---------------- perfusion / injury -----------------------------------------
KEDIN   : 0.30   : oedema formation (1/d)
KEDOUT  : 0.35   : oedema resolution (1/d)
KHYIN   : 0.30   : ventriculomegaly formation (1/d)
KHYOUT  : 0.22   : ventriculomegaly resolution (1/d)
KEVD    : 2.40   : EVD multiplier on CSF egress (-)
TAUICP  : 0.30   : ICP equilibration (d)
TAUMAP  : 0.20   : MAP equilibration (d)
KREP    : 1.55   : oxygen-debt repair (1/d)
KINFV   : 0.42   : infarct conversion rate (1/d)
OGD50   : 0.300  : oxygen debt producing half-maximal conversion (-)
OGDN    : 4.00   : Hill coefficient for conversion (-)
KEBIREC : 0.030  : early-brain-injury resolution (1/d)
KNAD    : 1.55   : natriuresis / SIADH sodium loss (mmol/L/d)
KNAR    : 0.42   : sodium regulation toward 140 (1/d)
KFLIN   : 0.56   : fluid-retention accumulation (1/d)
KFLOUT  : 0.50   : fluid-retention resolution (1/d)
// ---------------- haemodynamic constants -------------------------------------
CBF0    : 50.0   : reference CBF (mL/100g/min)
CPP0    : 80.0   : reference CPP (mmHg)
RL0     : 0.250  : baseline large-artery resistance (mmHg per mL/100g/min)
RA0     : 0.850  : baseline arteriolar resistance
RM0     : 0.500  : baseline microvascular resistance
FSEG    : 0.330  : fraction of large-artery resistance in the spastic segment
RAMIN   : 0.200  : minimum arteriolar resistance (vasodilatory limit)
RAMAX   : 2.200  : maximum arteriolar resistance
OEFMAX  : 0.850  : maximum O2 extraction fraction, homogeneous bed
OEF0    : 0.350  : baseline O2 extraction fraction
KCTH    : 1.750  : penalty on OEFmax per unit transit heterogeneity
// ---------------- drug PK ----------------------------------------------------
KANIM   : 36.0   : nimodipine absorption (1/d)
KELNIM  : 9.79   : nimodipine elimination, t1/2 1.7 h (1/d)
VNIM    : 100.0  : nimodipine apparent volume (L)
FNIM    : 0.13   : nimodipine oral bioavailability, CYP3A4 first pass (-)
KELCLZ  : 8.32   : clazosentan elimination, t1/2 2 h (1/d)
VCLZ    : 20.0   : clazosentan volume (L)
KELCIL  : 1.512  : cilostazol elimination, t1/2 11 h (1/d)
VCIL    : 90.0   : cilostazol apparent volume (L)
FCIL    : 0.90   : cilostazol bioavailability (-)
KELMIL  : 6.65   : milrinone elimination, t1/2 2.5 h (1/d)
VMIL    : 30.0   : milrinone volume (L)
KELNIC  : 1.00   : intrathecal nicardipine CSF elimination (1/d)
VNIC    : 0.50   : effective CSF distribution volume (L)
KELKET  : 1.00   : ketamine effect-site turnover (1/d)
TAUSTAT : 1.00   : statin effect-site time constant (d)
// ---------------- drug PD ----------------------------------------------------
EC50NIM   : 30.0   : nimodipine EC50 (ng/mL)
EMXNIM_SP : 0.16   : nimodipine Emax, large-artery relaxation (-)
EMXNIM_RM : 0.45   : nimodipine Emax, microvascular relaxation (-)
EMXNIM_SD : 0.55   : nimodipine Emax, SD suppression (-)
EMXNIM_NP : 0.50   : nimodipine Emax, infarct-conversion block (-)
DMAP_NIM_PO : 5.0  : MAP reduction, oral nimodipine (mmHg)
DMAP_NIM_IV : 15.0 : MAP reduction, IV nimodipine (mmHg)
EC50CLZ   : 400.0  : clazosentan EC50 (ng/mL)
EMXCLZ_ET : 0.70   : clazosentan Emax on ET-1 signalling (-)
EMXCLZ_RM : 0.08   : clazosentan Emax, microvascular (-)
DMAP_CLZ  : 6.0    : MAP reduction, clazosentan (mmHg)
DHGB_CLZ  : 0.80   : haemoglobin reduction, clazosentan (g/dL)
EC50CIL   : 600.0  : cilostazol EC50 (ng/mL)
EMXCIL_MT : 0.70   : cilostazol Emax, microthrombosis block (-)
EMXCIL_RM : 0.28   : cilostazol Emax, microvascular (-)
EMXCIL_SP : 0.10   : cilostazol Emax, large-artery (-)
EC50MIL   : 150.0  : milrinone EC50 (ng/mL)
EMXMIL_RM : 0.34   : milrinone Emax, microvascular (-)
EMXMIL_SP : 0.14   : milrinone Emax, large-artery (-)
DMAP_MIL  : -4.0   : MAP change with milrinone (mmHg, negative = rise)
EC50NIC   : 2000.0 : intrathecal nicardipine CSF EC50 (ng/mL)
EMXNIC_SP : 0.55   : nicardipine Emax, large-artery (-)
EMXNIC_RM : 0.22   : nicardipine Emax, microvascular (-)
EC50KET   : 1000.0 : ketamine EC50 (ng/mL)
EMXKET_SD : 0.65   : ketamine Emax, SD suppression (-)
EMXSTAT_NO  : 0.03 : statin Emax on NO bioavailability (-)
EMXSTAT_INF : 0.04 : statin Emax on inflammation resolution (-)
// ---------------- subject covariates ----------------------------------------
AGE     : 56.0   : age (years)
WFNS    : 3.0    : WFNS grade 1-5 (-)
MFI     : 3.0    : modified Fisher grade 1-4 (-)
HP22    : 0.0    : Hp2-2 genotype indicator (0/1)
MAP0    : 96.0   : baseline mean arterial pressure (mmHg)
HGB0    : 12.6   : baseline haemoglobin (g/dL)
THRP    : 1.0    : microthrombosis propensity multiplier (-)
SPSN    : 1.0    : large-artery vasoreactivity multiplier (-)
RMSN    : 1.0    : microvascular reactivity multiplier (-)
AREGB   : 0.86   : baseline autoregulatory gain (-)
EBIB    : 0.40   : baseline early-brain-injury burden (-)
RA0i    : 0.850  : subject baseline arteriolar resistance
RAMINi  : 0.200  : subject minimum arteriolar resistance
CLF     : 1.0    : clearance multiplier (-)
VRISKi  : 42.0   : tissue volume at risk (mL)
EVDF    : 1.0    : external ventricular drain in place (0/1)
// ---------------- treatment switches ----------------------------------------
NIMPO   : 0.0    : oral nimodipine dose per administration (mg)
NIMIV   : 0.0    : IV nimodipine rate (mg/h)
CLZR    : 0.0    : clazosentan rate (mg/h)
CILD    : 0.0    : cilostazol dose (mg bid)
STATON  : 0.0    : simvastatin 40 mg on/off (0/1)
MILR    : 0.0    : milrinone rate (ug/kg/min)
NICD    : 0.0    : intrathecal nicardipine release (mg/d)
KETON   : 0.0    : ketamine infusion on/off (0/1)
DRAIN   : 1.0    : drainage intensity multiplier (1 = none)
DRAINHB : 1.0    : CSF free-Hb clearance multiplier (-)
MAPH    : 0.0    : induced-hypertension MAP increment, d4-14 (mmHg)
HGBTGT  : 0.0    : transfusion threshold (g/dL, 0 = none)
TSTART  : 0.5    : treatment start (d)
TSTOP   : 21.0   : treatment stop (d)
KOLARGE : 0.0    : knockout, large-artery consumer (0/1)
KOMICRO : 0.0    : knockout, microvascular tone (0/1)
KOTHROMB: 0.0    : knockout, microthrombosis (0/1)
KOSD    : 0.0    : knockout, spreading depolarisation (0/1)
KOAREG  : 0.0    : knockout, autoregulatory failure (0/1)

$CMT @annotated
// ---- drug PK / exposure (8) ----
NIMG   : nimodipine gut depot (mg)
NIMC   : nimodipine plasma concentration (ng/mL)
CLAZ   : clazosentan plasma concentration (ng/mL)
CILO   : cilostazol plasma concentration (ng/mL)
STAT   : statin effect-site (fraction)
MILRC  : milrinone plasma concentration (ng/mL)
NICA   : intrathecal nicardipine CSF concentration (ng/mL)
KETA   : ketamine effect-site concentration (ng/mL)
// ---- clot / haemoglobin (8) ----
CLOT   : cisternal clot mass (relative)
IVH    : intraventricular clot mass (relative)
RBCL   : erythrocytes entering haemolysis (relative)
OXYHB  : CSF free oxyhaemoglobin (RU)
HP     : CSF haptoglobin scavenging capacity (RU)
HEME   : free haem (RU)
HO1    : haem oxygenase-1 induction (relative)
BOX    : bilirubin oxidation products (RU)
// ---- effectors (6) ----
ET1    : endothelin-1 signal (RU)
NOB    : NO bioavailability (fraction of normal)
ROS    : oxidative stress burden (RU)
RHOK   : Rho-kinase / Ca-sensitisation (RU)
INFL   : neuroinflammatory burden (RU)
PAI    : prothrombotic tone, PAI-1 dominant (RU)
// ---- vascular mechanics (5) ----
SPASM  : functional large-artery constriction (fraction)
STRUCT : structural arterial remodelling (fraction)
RMIC   : microvascular resistance multiplier (-)
MTHR   : capillary occlusion fraction (-)
AREG   : autoregulatory gain (-)
// ---- spreading depolarisation (3) ----
SDSUS  : SD susceptibility (-)
SDBUR  : SD burden trace (events/d)
SDCUM  : cumulative SD count (events)
// ---- perfusion / injury (7) ----
EDEMA  : global oedema (relative)
HYDRO  : ventriculomegaly (relative)
ICP    : intracranial pressure (mmHg)
MAP    : mean arterial pressure (mmHg)
OGD    : oxygen debt (-)
INFVOL : DCI infarct volume (mL)
EBI    : early brain injury burden (-)
// ---- systemic (2) ----
NAS    : serum sodium (mmol/L)
FLUID  : fluid-retention burden (-)

$MAIN
if (NEWIND < 2) {
  // modified Fisher grade sets the two clot compartments
  double c0 = (MFI < 1.5) ? 0.35 : ((MFI < 2.5) ? 0.42 : ((MFI < 3.5) ? 0.95 : 1.00));
  double i0 = (MFI < 1.5) ? 0.00 : ((MFI < 2.5) ? 0.60 : ((MFI < 3.5) ? 0.00 : 0.80));
  CLOT_0  = c0;
  IVH_0   = i0;
  HP_0    = 1.0;
  NOB_0   = 1.0;
  RMIC_0  = 1.0;
  AREG_0  = AREGB;
  ICP_0   = 8.0 + 1.3 + 5.0 * i0;
  MAP_0   = MAP0;
  NAS_0   = 139.0;
  EBI_0   = EBIB;
}

$ODE
// ===========================================================================
// 0.  ON/OFF and drug effects
// ===========================================================================
double ON   = (SOLVERTIME >= TSTART && SOLVERTIME <= TSTOP) ? 1.0 : 0.0;
double ONH  = (SOLVERTIME >= 4.0    && SOLVERTIME <= 14.0)  ? 1.0 : 0.0;

double ENIM  = NIMC  / (NIMC  + EC50NIM);
double ECLZ  = CLAZ  / (CLAZ  + EC50CLZ);
double ECIL  = CILO  / (CILO  + EC50CIL);
double EMIL  = MILRC / (MILRC + EC50MIL);
double ENIC  = NICA  / (NICA  + EC50NIC);
double EKET  = KETA  / (KETA  + EC50KET);
double ESTAT = (STAT < 0) ? 0.0 : ((STAT > 1) ? 1.0 : STAT);

// ===========================================================================
// 1.  HAEMODYNAMICS — the shared reserve (algebraic, evaluated every step)
// ===========================================================================
double SPT = SPASM + STRUCT;
if (SPT > 0.85) SPT = 0.85;
if (KOLARGE > 0.5) SPT = 0.0;

double RMICv = (RMIC < 1.0) ? 1.0 : RMIC;
if (KOMICRO > 0.5) RMICv = 1.0;

double MT = MTHR;
if (MT > 0.75) MT = 0.75;
if (MT < 0.0)  MT = 0.0;
if (KOTHROMB > 0.5) MT = 0.0;

// only the spastic SEGMENT obeys r^-4
double oneMinusS = 1.0 - SPT;
if (oneMinusS < 0.15) oneMinusS = 0.15;
double RL = RL0 * ((1.0 - FSEG) + FSEG / pow(oneMinusS, 4.0));

double RMtone = RM0 * (RMICv - 1.0);
double oneMinusMT = 1.0 - MT;
if (oneMinusMT < 0.25) oneMinusMT = 0.25;
double RMthr = RM0 * RMICv * (1.0 / oneMinusMT - 1.0);
double RM = RM0 + RMtone + RMthr;

double CPP = MAP - ICP;
if (CPP < 20.0)  CPP = 20.0;
if (CPP > 200.0) CPP = 200.0;

double SDfrac = SDBUR / SDMAX;
if (SDfrac < 0.0) SDfrac = 0.0;
if (SDfrac > 1.0) SDfrac = 1.0;
if (KOSD > 0.5) SDfrac = 0.0;
double demand = 1.0 + 0.30 * SDfrac;

double AREGv = (KOAREG > 0.5) ? 1.0 : AREG;
if (AREGv < 0.0) AREGv = 0.0;
if (AREGv > 1.0) AREGv = 1.0;

double SDS = SDSUS;
if (SDS < 0.0) SDS = 0.0;
if (SDS > 1.0) SDS = 1.0;

double CBFt  = CBF0 * demand;
double RAdes = CPP / CBFt - RL - RM;
if (RAdes < RAMINi) RAdes = RAMINi;
if (RAdes > RAMAX)  RAdes = RAMAX;
double RA = AREGv * RAdes + (1.0 - AREGv) * RA0i;
// SD inverse neurovascular coupling: in compromised tissue SD CONSTRICTS
RA = RA * (1.0 + KINV * SDfrac * SDS);
if (RA < RAMINi) RA = RAMINi;
if (RA > RAMAX)  RA = RAMAX;

double CBF = CPP / (RL + RA + RM);

// arterial oxygen content, relative to the reference subject
double hgb = HGB0;
if (HGBTGT > 0.0 && hgb < HGBTGT) hgb = HGBTGT;
if (CLZR > 0.0) hgb = hgb - DHGB_CLZ;
double CaO2 = (hgb * 1.34 * 0.97 + 0.285) / (12.6 * 1.34 * 0.97 + 0.285);

double DO2   = CBF * CaO2;
double CMRO2 = CBF0 * OEF0 * demand;

// capillary transit-time heterogeneity degrades maximum extraction
double CTH    = (RMICv - 1.0) + 2.5 * MT;
double oefmax = OEFMAX / (1.0 + KCTH * CTH);
double OEF    = CMRO2 / ((DO2 > 1e-6) ? DO2 : 1e-6);
if (OEF > oefmax) OEF = oefmax;

double ISCH = 1.0 - DO2 * oefmax / ((CMRO2 > 1e-6) ? CMRO2 : 1e-6);
if (ISCH < 0.0) ISCH = 0.0;
if (ISCH > 1.0) ISCH = 1.0;

// ===========================================================================
// 2.  DRUG PK
// ===========================================================================
dxdt_NIMG = -KANIM * NIMG;
dxdt_NIMC = KANIM * FNIM * NIMG / VNIM * 1000.0 - KELNIM * CLF * NIMC
            + (NIMIV > 0.0 ? NIMIV * 24.0 * 1000.0 / VNIM * ON : 0.0);
dxdt_CLAZ = -KELCLZ * CLF * CLAZ
            + (CLZR > 0.0 ? CLZR * 24.0 * 1000.0 / VCLZ * ON : 0.0);
dxdt_CILO = -KELCIL * CLF * CILO
            + (CILD > 0.0 ? CILD * 2.0 * FCIL * 1000.0 / VCIL * ON : 0.0);
dxdt_STAT = (STATON * ON - STAT) / TAUSTAT;
dxdt_MILRC = -KELMIL * CLF * MILRC
             + (MILR > 0.0 ? MILR * 70.0 * 1440.0 / 1000.0 * 1000.0 / VMIL * ONH : 0.0);
dxdt_NICA = -KELNIC * NICA + (NICD > 0.0 ? NICD * 1000.0 / VNIC * ON : 0.0);
dxdt_KETA = -KELKET * KETA + (KETON > 0.0 ? 1480.0 * KELKET * ON : 0.0);

// ===========================================================================
// 3.  CLOT -> ERYTHROCYTES -> OXYHAEMOGLOBIN  (the clock)
// ===========================================================================
double krem = KDRAIN * ((DRAIN > 1.0) ? (DRAIN - 1.0) : 0.0);
double mat  = KMAT  * CLOT;
double matv = KMATV * IVH;
double hemo = KHEM  * RBCL;

dxdt_CLOT = -mat  - krem * CLOT;
dxdt_IVH  = -matv - krem * IVH;
dxdt_RBCL = mat + 0.70 * matv - hemo - krem * RBCL;

double effhp = 1.0 - 0.40 * HP22;          // Hp2-2 scavenges less well
double hpav  = HP / (HP + KMHP);
double clrhb = KHPCL * effhp * hpav + KMAC * (1.0 + 3.0 * HO1) + KCSF * DRAINHB;

dxdt_OXYHB = KREL * hemo - clrhb * OXYHB;
dxdt_HP    = KHPIN * (1.0 + 0.8 * INFL) - KHPOUT * HP - KHPUSE * hpav * OXYHB;
dxdt_HEME  = KHEME * OXYHB - KHOX * (0.20 + HO1) * HEME;
double drvho1 = HEME + 0.5 * OXYHB;
dxdt_HO1   = KHO1IN * drvho1 / (1.0 + drvho1) - KHO1OUT * HO1;
dxdt_BOX   = KBOXIN * KHOX * (0.20 + HO1) * HEME - KBOXOUT * BOX;

// ===========================================================================
// 4.  EFFECTORS
// ===========================================================================
dxdt_ET1 = KETIN * (0.10 + OXYHB + 0.30 * ROS + 0.25 * INFL) - KETOUT * ET1;
double ETeff = ET1 * (1.0 - EMXCLZ_ET * ECLZ);

double nobc = (NOB < 0.0) ? 0.0 : ((NOB > 1.0) ? 1.0 : NOB);
double nob_t = (1.0 + EMXSTAT_NO * ESTAT) / (1.0 + KNOHB * OXYHB + KNOROS * ROS);
dxdt_NOB = (nob_t - NOB) / TAUNO;

dxdt_ROS = KROSIN * (0.10 + 0.90 * (HEME + 0.40 * BOX + 0.50 * INFL))
           - KROSOUT * ROS;

double rhok_t = 0.55 * ETeff + 0.30 * ROS + 0.35 * (1.0 - nobc);
dxdt_RHOK = (rhok_t - RHOK) / TAURHO;

dxdt_INFL = KINFIN * (OXYHB + 0.50 * IVH + 0.30 * EBI)
            - KINFOUT * INFL * (1.0 + EMXSTAT_INF * ESTAT);

dxdt_PAI = KPAIIN * (0.10 + INFL + 0.50 * OXYHB) - KPAIOUT * PAI;

// ===========================================================================
// 5.  VASCULAR MECHANICS — the four consumers
// ===========================================================================
double vasod_sp = EMXNIM_SP * ENIM + EMXNIC_SP * ENIC
                + EMXMIL_SP * EMIL + EMXCIL_SP * ECIL;
if (vasod_sp > 0.85) vasod_sp = 0.85;
if (vasod_sp < 0.0)  vasod_sp = 0.0;

double drive_sp = SPSN * (ETeff + 0.55 * RHOK + 0.35 * (1.0 - nobc) + 0.25 * BOX)
                  * (1.0 - vasod_sp);
dxdt_SPASM  = KSPIN * drive_sp * (SPMAX - SPASM) - KSPOUT * SPASM;
dxdt_STRUCT = KSTIN * SPASM * RHOK - KSTOUT * STRUCT;

double vasod_rm = EMXNIM_RM * ENIM + EMXCLZ_RM * ECLZ + EMXCIL_RM * ECIL
                + EMXMIL_RM * EMIL + EMXNIC_RM * ENIC;
if (vasod_rm > 0.80) vasod_rm = 0.80;
if (vasod_rm < 0.0)  vasod_rm = 0.0;

double rmic_t = 1.0 + RMSN * KRMIC * (0.24 * ETeff + 0.14 * (1.0 - nobc)
                + 0.10 * ROS + 0.08 * INFL) * (1.0 - vasod_rm);
dxdt_RMIC = (rmic_t - RMIC) / TAURM;

double anti = EMXCIL_MT * ECIL;
if (anti > 0.85) anti = 0.85;
dxdt_MTHR = KMTIN * PAI * (1.0 + 0.50 * INFL) * THRP * (1.0 - MTHR) * (1.0 - anti)
            - KMTOUT * MTHR;

dxdt_AREG = KAREGR * (AREGB - AREG)
            - KAREGL * AREG * (0.90 * ISCH + 0.10 * ROS + 0.07 * INFL + 0.25 * EBI);

// ===========================================================================
// 6.  SPREADING DEPOLARISATION
// ===========================================================================
double sds_t = 0.10 + 1.25 * ISCH + 0.25 * (1.0 - nobc) + 0.35 * EBI
             + 0.40 * INFVOL / ((VRISKi > 1.0) ? VRISKi : 1.0);
if (sds_t < 0.0) sds_t = 0.0;
if (sds_t > 1.0) sds_t = 1.0;
dxdt_SDSUS = (sds_t - SDSUS) / TAUSDS;

double sdrate = SDMAX * SDS * SDS * (1.0 - EMXNIM_SD * ENIM - EMXKET_SD * EKET);
if (sdrate < 0.0) sdrate = 0.0;
dxdt_SDBUR = sdrate - KSDB * SDBUR;
dxdt_SDCUM = sdrate;

// ===========================================================================
// 7.  PERFUSION, ICP, INJURY
// ===========================================================================
dxdt_EDEMA = KEDIN * (ISCH + 0.40 * INFL
             + 0.50 * INFVOL / ((VRISKi > 1.0) ? VRISKi : 1.0)) - KEDOUT * EDEMA;
dxdt_HYDRO = KHYIN * (IVH + 0.30 * INFL) - KHYOUT * HYDRO * (1.0 + KEVD * EVDF);

double icp_t = 8.0 + 13.0 * HYDRO + 11.0 * EDEMA + 5.0 * IVH;
dxdt_ICP = (icp_t - ICP) / TAUICP;

double map_t = MAP0 + MAPH * ONH
             - (NIMPO > 0.0 ? DMAP_NIM_PO * ENIM : 0.0)
             - (NIMIV > 0.0 ? DMAP_NIM_IV * ENIM : 0.0)
             - DMAP_CLZ * ECLZ - DMAP_MIL * EMIL;
dxdt_MAP = (map_t - MAP) / TAUMAP;

dxdt_OGD = ISCH - KREP * OGD;

double ogdp = pow((OGD > 0.0 ? OGD : 0.0), OGDN);
double conv = ogdp / (ogdp + pow(OGD50, OGDN));
double npro = 1.0 - EMXNIM_NP * ENIM;
double atrisk = VRISKi - INFVOL;
if (atrisk < 0.0) atrisk = 0.0;
dxdt_INFVOL = KINFV * atrisk * conv * npro;

dxdt_EBI = -KEBIREC * EBI;

dxdt_NAS = -KNAD * (0.50 * IVH + 0.40 * INFL + 0.30 * EBI) + KNAR * (140.0 - NAS);
dxdt_FLUID = KFLIN * ECLZ + 0.25 * ((MAPH > 0.0) ? ONH : 0.0) - KFLOUT * FLUID;

$TABLE
// --- recompute observables for output -------------------------------------
double oSPT = SPASM + STRUCT;
if (oSPT > 0.85) oSPT = 0.85;
if (KOLARGE > 0.5) oSPT = 0.0;
double oRMIC = (RMIC < 1.0) ? 1.0 : RMIC;
if (KOMICRO > 0.5) oRMIC = 1.0;
double oMT = MTHR; if (oMT > 0.75) oMT = 0.75; if (oMT < 0) oMT = 0;
if (KOTHROMB > 0.5) oMT = 0.0;
double oms = 1.0 - oSPT; if (oms < 0.15) oms = 0.15;
double oRL = RL0 * ((1.0 - FSEG) + FSEG / pow(oms, 4.0));
double oRMt = RM0 * (oRMIC - 1.0);
double omt2 = 1.0 - oMT; if (omt2 < 0.25) omt2 = 0.25;
double oRMc = RM0 * oRMIC * (1.0 / omt2 - 1.0);
double oRM = RM0 + oRMt + oRMc;
double oCPP = MAP - ICP; if (oCPP < 20) oCPP = 20; if (oCPP > 200) oCPP = 200;
double oSDf = SDBUR / SDMAX; if (oSDf < 0) oSDf = 0; if (oSDf > 1) oSDf = 1;
if (KOSD > 0.5) oSDf = 0.0;
double odem = 1.0 + 0.30 * oSDf;
double oAREG = (KOAREG > 0.5) ? 1.0 : AREG;
if (oAREG < 0) oAREG = 0; if (oAREG > 1) oAREG = 1;
double oSDS = SDSUS; if (oSDS < 0) oSDS = 0; if (oSDS > 1) oSDS = 1;
double oRAdes = oCPP / (CBF0 * odem) - oRL - oRM;
if (oRAdes < RAMINi) oRAdes = RAMINi;
if (oRAdes > RAMAX)  oRAdes = RAMAX;
double oRA = oAREG * oRAdes + (1.0 - oAREG) * RA0i;
oRA = oRA * (1.0 + KINV * oSDf * oSDS);
if (oRA < RAMINi) oRA = RAMINi;
if (oRA > RAMAX)  oRA = RAMAX;
double oCBF = oCPP / (oRL + oRA + oRM);
double ohgb = HGB0;
if (HGBTGT > 0.0 && ohgb < HGBTGT) ohgb = HGBTGT;
if (CLZR > 0.0) ohgb = ohgb - DHGB_CLZ;
double oCaO2 = (ohgb * 1.34 * 0.97 + 0.285) / (12.6 * 1.34 * 0.97 + 0.285);
double oDO2 = oCBF * oCaO2;
double oCMRO2 = CBF0 * OEF0 * odem;
double oCTH = (oRMIC - 1.0) + 2.5 * oMT;
double ooefmax = OEFMAX / (1.0 + KCTH * oCTH);
double oOEF = oCMRO2 / ((oDO2 > 1e-6) ? oDO2 : 1e-6);
if (oOEF > ooefmax) oOEF = ooefmax;
double oISCH = 1.0 - oDO2 * ooefmax / ((oCMRO2 > 1e-6) ? oCMRO2 : 1e-6);
if (oISCH < 0) oISCH = 0; if (oISCH > 1) oISCH = 1;

capture CALIBER  = 100.0 * (1.0 - oSPT);        // % of baseline diameter
capture SPTOT    = oSPT;                        // fractional calibre loss
capture CBFo     = oCBF;                        // mL/100g/min
capture CPPo     = oCPP;                        // mmHg
capture ISCHo    = oISCH;                        // 0-1
capture PBTO2    = 25.0 * (oDO2 * (1.0 - oOEF)) / (CBF0 * (1.0 - OEF0));
capture TCD      = 60.0 * (oCBF / CBF0) / pow((oms > 0.20 ? oms : 0.20), 2.0);
capture OEFo     = oOEF;
capture OEFMAXo  = ooefmax;
capture CTHo     = oCTH;
capture RLo      = oRL;
capture RAo      = oRA;
capture RMo      = oRM;
capture RESERVE  = RA0i - RAMINi;
capture DLARGE   = oRL - RL0;                   // reserve demand, large artery
capture DMICRO   = oRMt;                        // reserve demand, micro tone
capture DTHROMB  = oRMc;                        // reserve demand, thrombosis
capture DCPPd    = (CPP0 - oCPP > 0 ? CPP0 - oCPP : 0) / CBF0;
capture DSDd     = oRA * KINV * oSDf * oSDS;
capture CBFMAX   = oCPP / (oRL + RAMINi + oRM);
capture CBFREQ   = oCMRO2 / ((oCaO2 * ooefmax > 1e-6) ? oCaO2 * ooefmax : 1e-6);
capture MARGIN   = CBFMAX / ((CBFREQ > 1e-6) ? CBFREQ : 1e-6);
capture LINDEG   = TCD / 40.0;
'

sah_mod <- mcode("sah_dci", sah_code)

## =============================================================================
## VIRTUAL POPULATION
## =============================================================================
sah_population <- function(n = 900, seed = 20260728) {
  set.seed(seed)
  wfns <- sample(1:5, n, TRUE, prob = c(.34, .16, .10, .24, .16))
  mfi  <- sample(1:4, n, TRUE, prob = c(.10, .20, .34, .36))
  ivh0 <- c(0, .60, 0, .80)[mfi]
  data.frame(
    ID      = seq_len(n),
    AGE     = pmin(pmax(rnorm(n, 56, 12), 25), 88),
    WFNS    = wfns,
    MFI     = mfi,
    HP22    = rbinom(n, 1, 0.36),
    MAP0    = pmin(pmax(rnorm(n, 96, 11), 70), 132),
    HGB0    = pmin(pmax(rnorm(n, 12.6, 1.7), 7.5), 17),
    THRP    = rlnorm(n, 0, 0.90),
    SPSN    = rlnorm(n, 0, 0.45),
    RMSN    = rlnorm(n, 0, 0.30),
    AREGB   = pmin(pmax(rnorm(n, 0.86, 0.11) - 0.055 * (wfns - 1), 0.10), 0.99),
    EBIB    = pmin(pmax(0.08 + 0.155 * (wfns - 1) + rnorm(n, 0, 0.07), 0), 1),
    RA0i    = 0.850 * rlnorm(n, 0, 0.12),
    RAMINi  = 0.200 * rlnorm(n, 0, 0.15),
    CLF     = rlnorm(n, 0, 0.30),
    VRISKi  = 42 * rlnorm(n, 0, 0.65),
    EVDF    = as.numeric(ivh0 > 0.30)
  )
}

## =============================================================================
## FIFTEEN TREATMENT SCENARIOS
## =============================================================================
SOC <- list(NIMPO = 60, TSTART = 0.5, TSTOP = 21)

sah_scenario_defs <- list(
  S1  = list(label = "supportive only (no nimodipine)"),
  S2  = c(list(label = "SoC: oral nimodipine 60 mg q4h x21d"), SOC),
  S3  = list(label = "IV nimodipine 2 mg/h x14d", NIMIV = 2, TSTART = 0.5, TSTOP = 14),
  S4  = c(list(label = "SoC + clazosentan 15 mg/h d1-14", CLZR = 15), SOC),
  S5  = c(list(label = "SoC + clazosentan 5 mg/h d1-14",  CLZR = 5),  SOC),
  S6  = c(list(label = "SoC + cilostazol 100 mg bid", CILD = 100), SOC),
  S7  = c(list(label = "SoC + simvastatin 40 mg qd", STATON = 1), SOC),
  S8  = c(list(label = "SoC + early lumbar drainage", DRAIN = 2.20, DRAINHB = 1.85), SOC),
  S9  = c(list(label = "SoC + intrathecal nicardipine implant", NICD = 4), SOC),
  S10 = c(list(label = "SoC + induced hypertension (MAP +20)", MAPH = 20), SOC),
  S11 = c(list(label = "SoC + milrinone 0.5 ug/kg/min d4-14", MILR = 0.5), SOC),
  S12 = c(list(label = "SoC + ketamine (SD suppression)", KETON = 1), SOC),
  S13 = c(list(label = "SoC + transfusion to Hb 10 g/dL", HGBTGT = 10), SOC),
  S14 = c(list(label = "SoC + clazosentan + cilostazol + drainage",
               CLZR = 15, CILD = 100, DRAIN = 2.20, DRAINHB = 1.85), SOC),
  S15 = c(list(label = "SoC + full multi-path stack",
               CLZR = 15, CILD = 100, DRAIN = 2.20, DRAINHB = 1.85,
               NICD = 4, KETON = 1, MAPH = 15), SOC)
)

sah_knockout_defs <- list(
  K0    = c(list(label = "SoC (reference)"), SOC),
  K1    = c(list(label = "perfect large-artery block", KOLARGE = 1), SOC),
  K2    = c(list(label = "perfect microvascular-tone block", KOMICRO = 1), SOC),
  K3    = c(list(label = "perfect microthrombosis block", KOTHROMB = 1), SOC),
  K4    = c(list(label = "perfect SD block", KOSD = 1), SOC),
  K5    = c(list(label = "autoregulation preserved", KOAREG = 1), SOC),
  K12   = c(list(label = "large + micro", KOLARGE = 1, KOMICRO = 1), SOC),
  K13   = c(list(label = "large + thrombosis", KOLARGE = 1, KOTHROMB = 1), SOC),
  K23   = c(list(label = "micro + thrombosis", KOMICRO = 1, KOTHROMB = 1), SOC),
  K1234 = c(list(label = "all four consumers blocked", KOLARGE = 1, KOMICRO = 1,
                 KOTHROMB = 1, KOSD = 1), SOC)
)

## Oral nimodipine is dosed q4h into the gut depot; everything else is either a
## continuous infusion or a switch handled inside $ODE.
sah_events <- function(pop, sc) {
  if (is.null(sc$NIMPO) || sc$NIMPO <= 0) return(NULL)
  t0 <- if (is.null(sc$TSTART)) 0.5 else sc$TSTART
  t1 <- if (is.null(sc$TSTOP))  21   else sc$TSTOP
  n  <- floor((t1 - t0) * 6) + 1
  ev(amt = sc$NIMPO, cmt = "NIMG", time = t0, ii = 1 / 6, addl = n - 1)
}

sah_simulate <- function(pop, sc, end = 21, delta = 0.1) {
  pars <- sc[setdiff(names(sc), "label")]
  idata <- pop
  for (nm in names(pars)) idata[[nm]] <- pars[[nm]]
  e <- sah_events(pop, sc)
  if (is.null(e)) {
    sah_mod %>% idata_set(idata) %>% mrgsim(end = end, delta = delta) %>% as_tibble()
  } else {
    sah_mod %>% idata_set(idata) %>% ev(e) %>%
      mrgsim(end = end, delta = delta) %>% as_tibble()
  }
}

## =============================================================================
## ENDPOINTS
## =============================================================================
##  angiographic vasospasm : peak calibre loss over days 4-11
##                           any >= 25%, moderate-severe >= 33%, severe >= 50%
##  DCI                    : new infarct >= 3 mL after day 3, OR cumulative
##                           time with ISCH > 0.15 after day 3 >= 0.25 d
##  poor outcome           : logistic in EBI, log(1+infarct), ICP burden, age,
##                           SD burden, treatment harm, hyponatraemia, rebleed
## =============================================================================
sah_endpoints <- function(sim, pop) {
  b <- list(B0 = -4.10, B_EBI = 3.15, B_INF = 1.08, B_HYD = 1.05,
            B_AGE = 0.036, B_HARM = 1.55, B_NA = 0.45, B_SD = 0.85)
  sim %>%
    group_by(ID) %>%
    summarise(
      spt_max = max(SPTOT[time >= 4 & time <= 11]),
      tcd_max = max(TCD[time >= 4 & time <= 11]),
      pbt_min = min(PBTO2[time >= 3]),
      icp_max = max(ICP[time >= 2]),
      na_min  = min(NAS),
      fluid_max = max(FLUID),
      areg_min  = min(AREG),
      sdcum   = last(SDCUM),
      infvol  = last(INFVOL),
      inf_d3  = INFVOL[which.min(abs(time - 3))],
      tisch   = sum((ISCHo > 0.15) & (time > 3)) * 0.1,
      .groups = "drop"
    ) %>%
    left_join(pop[, c("ID", "AGE", "EBIB", "MFI", "HP22", "AREGB", "HGB0")],
              by = "ID") %>%
    mutate(
      dci_inf = infvol - inf_d3,
      dci     = as.numeric(dci_inf >= 3 | tisch >= 0.25),
      ang_any = as.numeric(spt_max >= 0.25),
      ang_ms  = as.numeric(spt_max >= 0.33),
      ang_sev = as.numeric(spt_max >= 0.50),
      harm    = 0.85 * pmax(0, fluid_max - 0.35),
      logit   = b$B0 + b$B_EBI * EBIB + b$B_INF * log1p(infvol) +
                b$B_HYD * pmin(pmax((icp_max - 15) / 15, 0), 2) +
                b$B_AGE * (AGE - 56) + b$B_HARM * harm +
                b$B_SD * log1p(sdcum) / 3 + b$B_NA * (na_min < 130),
      ppoor   = 1 / (1 + exp(-logit))
    )
}

## =============================================================================
## RESULT 1 — SCENARIO TABLE
## =============================================================================
sah_scenarios <- function(pop = sah_population()) {
  out <- lapply(names(sah_scenario_defs), function(k) {
    sc <- sah_scenario_defs[[k]]
    ep <- sah_endpoints(sah_simulate(pop, sc), pop)
    data.frame(
      arm      = k,
      label    = sc$label,
      ang_any  = 100 * mean(ep$ang_any),
      ang_ms   = 100 * mean(ep$ang_ms),
      ang_sev  = 100 * mean(ep$ang_sev),
      dci      = 100 * mean(ep$dci),
      inf_dci  = ifelse(sum(ep$dci) > 0, median(ep$infvol[ep$dci > 0.5]), 0),
      pbto2    = median(ep$pbt_min),
      poor     = 100 * mean(ep$ppoor)
    )
  })
  do.call(rbind, out)
}

## =============================================================================
## RESULT 2 — RISK RATIOS AGAINST STANDARD OF CARE, WITH TRIAL COMPARATORS
## =============================================================================
sah_risk_ratios <- function(pop = sah_population()) {
  eps <- lapply(sah_scenario_defs, function(sc) sah_endpoints(sah_simulate(pop, sc), pop))
  ref <- eps$S2
  rr <- do.call(rbind, lapply(names(eps), function(k) {
    e <- eps[[k]]
    data.frame(arm = k, label = sah_scenario_defs[[k]]$label,
               rr_ang_ms = mean(e$ang_ms) / mean(ref$ang_ms),
               rr_dci    = mean(e$dci)    / mean(ref$dci),
               rr_poor   = mean(e$ppoor)  / mean(ref$ppoor))
  }))
  targets <- data.frame(
    arm    = c("S1", "S4", "S6", "S7", "S8"),
    target = c("nimodipine RR poor 0.67 (inverse: Cochrane)",
               "clazosentan RR angMS 0.35 / RR poor 1.05 (CONSCIOUS-1/2)",
               "cilostazol RR DCI 0.47 (meta-analysis)",
               "simvastatin RR poor ~1.00 (STASH)",
               "lumbar drainage RR poor 0.76 (EARLYDRAIN)"))
  merge(rr, targets, by = "arm", all.x = TRUE)
}

## =============================================================================
## RESULT 3 — THE DISSOCIATION, AND TCD vs PbtO2 AS DIAGNOSTIC TESTS
## =============================================================================
sah_dissociation <- function(pop = sah_population()) {
  ep <- sah_endpoints(sah_simulate(pop, sah_scenario_defs$S2), pop)
  tests <- function(x, thr, above = TRUE) {
    pos <- if (above) x > thr else x < thr
    c(sens = 100 * sum(pos & ep$dci > 0.5) / sum(ep$dci > 0.5),
      spec = 100 * sum(!pos & ep$dci < 0.5) / sum(ep$dci < 0.5))
  }
  list(
    r2_caliber_infarct = cor(ep$spt_max, ep$infvol)^2,
    r2_tcd_infarct     = cor(ep$tcd_max, ep$infvol)^2,
    dci_without_modsev_spasm = 100 * sum(ep$dci > 0.5 & ep$ang_ms < 0.5) / sum(ep$dci > 0.5),
    modsev_spasm_without_dci = 100 * sum(ep$ang_ms > 0.5 & ep$dci < 0.5) / sum(ep$ang_ms > 0.5),
    tcd120  = tests(ep$tcd_max, 120),
    tcd160  = tests(ep$tcd_max, 160),
    tcd200  = tests(ep$tcd_max, 200),
    pbto2_20 = tests(ep$pbt_min, 20, above = FALSE),
    pbto2_15 = tests(ep$pbt_min, 15, above = FALSE)
  )
}

## =============================================================================
## RESULT 4 — WHICH CONSUMER IS EATING THE RESERVE, BY DAY
## =============================================================================
sah_reserve_decomp <- function(pop = sah_population(), days = c(2, 4, 6, 8, 10, 14)) {
  sim <- sah_simulate(pop, sah_scenario_defs$S2)
  do.call(rbind, lapply(days, function(d) {
    s <- sim[abs(sim$time - d) < 1e-6, ]
    comp <- c(large = mean(s$DLARGE), micro = mean(s$DMICRO),
              thromb = mean(s$DTHROMB), cpp = mean(s$DCPPd), sd = mean(s$DSDd))
    data.frame(day = d, t(round(100 * comp / sum(comp), 1)),
               demand_over_reserve = round(mean((s$DLARGE + s$DMICRO + s$DTHROMB +
                                                 s$DCPPd + s$DSDd) / s$RESERVE), 2),
               large_share_of_total_CVR =
                 round(100 * mean(s$RLo / (s$CPPo / s$CBFo)), 1))
  }))
}

## =============================================================================
## RESULT 5 — THE EMERGENT WINDOW (there is no day-4 switch in the model)
## =============================================================================
sah_window <- function(pop = sah_population()) {
  sim <- sah_simulate(pop, sah_scenario_defs$S2)
  ep  <- sah_endpoints(sim, pop)
  grp <- function(ids, lab) {
    s <- sim[sim$ID %in% ids, ]
    a <- s %>% group_by(time) %>%
      summarise(hb = mean(OXYHB), et = mean(ET1), sp = mean(SPTOT),
                isch = mean(ISCHo), .groups = "drop")
    data.frame(group = lab,
               oxyhb_peak_day = a$time[which.max(a$hb)],
               oxyhb_peak     = round(max(a$hb), 2),
               et1_peak_day   = a$time[which.max(a$et)],
               caliber_nadir_day = a$time[which.max(a$sp)],
               caliber_loss_pct  = round(100 * max(a$sp)),
               ischaemia_peak_day = a$time[which.max(a$isch)])
  }
  bind_rows(
    grp(pop$ID, "all"),
    grp(pop$ID[pop$MFI <= 2], "mFisher 1-2"),
    grp(pop$ID[pop$MFI >= 3], "mFisher 3-4"),
    grp(pop$ID[pop$HP22 == 0], "Hp1-1/1-2"),
    grp(pop$ID[pop$HP22 == 1], "Hp2-2")
  ) %>%
    mutate(dci_note = "see sah_window_incidence()") ->
    tab
  attr(tab, "dci_by_mfisher") <- sapply(1:4, function(g)
    100 * mean(ep$dci[ep$MFI == g]))
  attr(tab, "dci_by_hp") <- c(
    `Hp1-1/1-2` = 100 * mean(ep$dci[ep$HP22 == 0]),
    `Hp2-2`     = 100 * mean(ep$dci[ep$HP22 == 1]))
  attr(tab, "dci_by_areg_tertile") <- {
    q <- quantile(ep$AREGB, c(1/3, 2/3))
    c(low  = 100 * mean(ep$dci[ep$AREGB <= q[1]]),
      mid  = 100 * mean(ep$dci[ep$AREGB > q[1] & ep$AREGB < q[2]]),
      high = 100 * mean(ep$dci[ep$AREGB >= q[2]]))
  }
  tab
}

## =============================================================================
## RESULT 6 — PATH KNOCKOUTS AND SUPERADDITIVITY
## =============================================================================
sah_knockouts <- function(pop = sah_population()) {
  eps <- lapply(sah_knockout_defs, function(sc) sah_endpoints(sah_simulate(pop, sc), pop))
  k0 <- eps$K0
  tab <- do.call(rbind, lapply(names(eps), function(k) {
    e <- eps[[k]]
    data.frame(knockout = k, label = sah_knockout_defs[[k]]$label,
               dci = 100 * mean(e$dci),
               rr_dci = mean(e$dci) / mean(k0$dci),
               inf_dci = ifelse(sum(e$dci) > 0, median(e$infvol[e$dci > 0.5]), 0),
               rr_poor = mean(e$ppoor) / mean(k0$ppoor))
  }))
  a  <- 1 - mean(eps$K1$dci) / mean(k0$dci)
  b  <- 1 - mean(eps$K3$dci) / mean(k0$dci)
  ab <- 1 - mean(eps$K13$dci) / mean(k0$dci)
  attr(tab, "superadditivity_large_plus_thrombosis") <-
    c(single_large = a, single_thromb = b, sum = a + b, joint = ab, excess = ab - (a + b))
  tab
}

## =============================================================================
## RESULT 7 — PRESSURE-BASED RESCUE ONLY WORKS WHEN THE RESERVE IS PASSIVE
## =============================================================================
sah_pressure_test <- function(pop = sah_population()) {
  ref <- sah_endpoints(sah_simulate(pop, sah_scenario_defs$S2), pop)
  q   <- quantile(pop$AREGB, c(0.33, 0.67))
  lo  <- pop$ID[pop$AREGB <= q[1]]
  hi  <- pop$ID[pop$AREGB >= q[2]]
  do.call(rbind, lapply(c("S10", "S11"), function(k) {
    e <- sah_endpoints(sah_simulate(pop, sah_scenario_defs[[k]]), pop)
    data.frame(arm = k, label = sah_scenario_defs[[k]]$label,
      rr_dci_impaired_autoreg =
        mean(e$dci[e$ID %in% lo]) / mean(ref$dci[ref$ID %in% lo]),
      rr_dci_intact_autoreg =
        mean(e$dci[e$ID %in% hi]) / mean(ref$dci[ref$ID %in% hi]))
  }))
}

## =============================================================================
## RESULT 8 — CLAZOSENTAN: WHERE THE ANGIOGRAPHIC WIN IS LOST
## =============================================================================
sah_clazosentan <- function(pop = sah_population()) {
  ref <- sah_endpoints(sah_simulate(pop, sah_scenario_defs$S2), pop)
  clz <- sah_endpoints(sah_simulate(pop, sah_scenario_defs$S4), pop)
  ## counterfactual: identical ETA blockade, harm channel switched off
  sc_b <- sah_scenario_defs$S4
  sc_b$DMAP_CLZ <- 0; sc_b$DHGB_CLZ <- 0; sc_b$KFLIN <- 0
  ben <- sah_endpoints(sah_simulate(pop, sc_b), pop)
  list(
    ang_ms       = c(soc = 100 * mean(ref$ang_ms), clz = 100 * mean(clz$ang_ms),
                     rr = mean(clz$ang_ms) / mean(ref$ang_ms)),
    caliber_loss = c(soc = 100 * median(ref$spt_max), clz = 100 * median(clz$spt_max)),
    dci          = c(soc = 100 * mean(ref$dci), clz = 100 * mean(clz$dci),
                     rr = mean(clz$dci) / mean(ref$dci)),
    poor         = c(soc = 100 * mean(ref$ppoor), clz = 100 * mean(clz$ppoor),
                     rr = mean(clz$ppoor) / mean(ref$ppoor)),
    fluid        = c(soc = median(ref$fluid_max), clz = median(clz$fluid_max),
                     pulm_oedema_range_soc = 100 * mean(ref$fluid_max > 0.45),
                     pulm_oedema_range_clz = 100 * mean(clz$fluid_max > 0.45)),
    harm_off_counterfactual =
      c(poor = 100 * mean(ben$ppoor), rr = mean(ben$ppoor) / mean(ref$ppoor)),
    interpretation = paste("the angiographic effect is real;",
                           "the endpoint effect is eaten by redundancy first",
                           "and by the harm channel second")
  )
}

## =============================================================================
## RESULT 9 — MONITORING LEAD TIMES
## =============================================================================
sah_monitoring <- function(pop = sah_population()) {
  sim <- sah_simulate(pop, sah_scenario_defs$S2)
  ep  <- sah_endpoints(sim, pop)
  dci_ids <- ep$ID[ep$dci > 0.5]
  s <- sim[sim$ID %in% dci_ids, ]
  t_tcd <- s %>% group_by(ID) %>%
    summarise(t = { i <- which(TCD > 120 & time > 1); if (length(i)) time[i[1]] else NA },
              .groups = "drop") %>% pull(t)
  t_pbt <- s %>% group_by(ID) %>%
    summarise(t = { i <- which(PBTO2 < 20 & time > 1); if (length(i)) time[i[1]] else NA },
              .groups = "drop") %>% pull(t)
  t_inf <- s %>% group_by(ID) %>%
    summarise(t = { d <- c(0, diff(INFVOL)); i <- which(d > 0.05 & time > 1)
                    if (length(i)) time[i[1]] else NA }, .groups = "drop") %>% pull(t)
  summ <- function(x, lab) {
    x <- x[!is.na(x)]
    data.frame(signal = lab, median_day = median(x),
               q25 = quantile(x, .25), q75 = quantile(x, .75), n = length(x))
  }
  rbind(summ(t_tcd, "TCD > 120 cm/s"),
        summ(t_pbt, "PbtO2 < 20 mmHg"),
        summ(t_inf, "infarct starts"))
}

## =============================================================================
## RESULT 10 — SYSTEMIC CHANNELS
## =============================================================================
sah_systemic <- function(pop = sah_population()) {
  ep <- sah_endpoints(sah_simulate(pop, sah_scenario_defs$S2), pop)
  m <- ep$dci > 0.5
  list(
    hyponatraemia = c(`<135` = 100 * mean(ep$na_min < 135),
                      `<130` = 100 * mean(ep$na_min < 130)),
    icp_over_20   = 100 * mean(ep$icp_max > 20),
    sd_count      = c(all = median(ep$sdcum), dci = median(ep$sdcum[m]),
                      no_dci = median(ep$sdcum[!m])),
    autoreg_nadir = c(all = median(ep$areg_min), dci = median(ep$areg_min[m]),
                      no_dci = median(ep$areg_min[!m]))
  )
}

## =============================================================================
## RESULT 11 — ANAEMIA AND OXYGEN CONTENT ACT ON THE SAME NODE
## =============================================================================
sah_anaemia <- function(pop = sah_population()) {
  ref <- sah_endpoints(sah_simulate(pop, sah_scenario_defs$S2), pop)
  trf <- sah_endpoints(sah_simulate(pop, sah_scenario_defs$S13), pop)
  strata <- list(`Hb < 10` = ref$HGB0 < 10,
                 `Hb 10-13` = ref$HGB0 >= 10 & ref$HGB0 < 13,
                 `Hb >= 13` = ref$HGB0 >= 13)
  tab <- do.call(rbind, lapply(names(strata), function(k) {
    m <- strata[[k]]
    data.frame(stratum = k, n = sum(m), dci = 100 * mean(ref$dci[m]),
               median_infarct = median(ref$infvol[m]))
  }))
  m <- ref$HGB0 < 10
  attr(tab, "transfusing_only_the_anaemic") <-
    c(n = sum(m), dci_before = 100 * mean(ref$dci[m]),
      dci_after = 100 * mean(trf$dci[m]),
      rr = mean(trf$dci[m]) / max(mean(ref$dci[m]), 1e-9))
  tab
}

## =============================================================================
## RESULT 12 — ROUTE MATTERS BECAUSE THE HARM CHANNEL IS CPP
## =============================================================================
sah_route <- function(pop = sah_population()) {
  arms <- c("S1", "S2", "S3")
  do.call(rbind, lapply(arms, function(k) {
    sim <- sah_simulate(pop, sah_scenario_defs[[k]])
    ep  <- sah_endpoints(sim, pop)
    w   <- sim$time >= 4 & sim$time <= 11
    data.frame(arm = k, label = sah_scenario_defs[[k]]$label,
               ang_ms = 100 * mean(ep$ang_ms), dci = 100 * mean(ep$dci),
               poor = 100 * mean(ep$ppoor),
               mean_CPP_d4_11 = mean(sim$CPPo[w]))
  }))
}

## =============================================================================
## RUN EVERYTHING
## =============================================================================
sah_run_all <- function(n = 900) {
  pop <- sah_population(n)
  list(
    scenarios     = sah_scenarios(pop),
    risk_ratios   = sah_risk_ratios(pop),
    dissociation  = sah_dissociation(pop),
    reserve       = sah_reserve_decomp(pop),
    window        = sah_window(pop),
    knockouts     = sah_knockouts(pop),
    pressure_test = sah_pressure_test(pop),
    clazosentan   = sah_clazosentan(pop),
    monitoring    = sah_monitoring(pop),
    systemic      = sah_systemic(pop),
    anaemia       = sah_anaemia(pop),
    route         = sah_route(pop)
  )
}

## Example:
##   res <- sah_run_all(300)
##   res$scenarios
##   res$reserve
##   res$clazosentan
##
## The authoritative numbers, and the virtual-population sizes used to produce
## them, are in sah_reference_check_output.txt (see sah_reference_check.py).
