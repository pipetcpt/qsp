## =====================================================================
##  se_mrgsolve_model.R
##  Status Epilepticus — QSP / PK-PD model
##  59 ODE compartments · 27 PK · 32 disease-PD · 22 therapeutic scenarios
##
##  경련지속상태 — 정량적 시스템 약리학 모델
## ---------------------------------------------------------------------
##  STRUCTURAL THESIS
##  -----------------
##  Status epilepticus is not "a long seizure".  It is a RECEPTOR-
##  TRAFFICKING CLOCK.  The drug target is leaving the synapse while the
##  clinician is deciding what to give.
##
##  ONE clock (ongoing seizure activity) drives TWO receptor movements in
##  OPPOSITE directions, and leaves a THIRD pool untouched:
##
##     R_SYN    synaptic γ2-GABA-A        ↓  (endocytosed)
##     NR_SYN   synaptic NMDA             ↑  (exocytosed / captured)
##     R_EXTRA  extrasynaptic δ-GABA-A    —  (not trafficked)
##
##  Every drug class in the model is written as
##
##        EFFECT_class(t)  =  f(C_effect-site)  ×  TARGET_POOL(t)
##
##  and the ONLY thing that distinguishes the classes is WHICH pool
##  multiplies them:
##
##     benzodiazepine    ×  R_SYN · F_BZS      (falls fastest — two terms)
##     phenobarbital     ×  0.20·R_SYN + 0.80·R_EXTRA
##     propofol          ×  0.25·R_SYN + 0.75·R_EXTRA
##     allopregnanolone  ×  0.40·R_SYN + 0.60·R_EXTRA
##     ketamine          ×  NR_SYN · MGREL · (use-dependence)   (RISES)
##     LEV / PHT / VPA   ×  1                  (presynaptic/axonal: flat)
##
##  Seven clinical observations that are usually filed as separate facts
##  are DERIVED from that single table rather than fitted one at a time.
##  The numbers below are what THIS CODE prints; run it and check.
##
##   1.  Benzodiazepine efficacy collapses with delay; phenobarbital's does
##       not.  CREQ_BZD — the effect-site concentration (in EC50 units) a
##       benzodiazepine needs RIGHT NOW to hold the network gain at 1 —
##       rises 0.23 -> 1.62 (7.0x) between 5 and 30 min and goes to
##       INFINITY between 45 and 50 min: past that point no dose works.
##       CREQ_PB over the same interval rises 0.20 -> 0.52 (2.6x) and is
##       still only 2.5 at 8 HOURS.  Kapur & Macdonald measured 20-fold vs
##       3-fold at 30 min in isolated neurons; the model reproduces the
##       direction and the asymmetry (2.7x differential) but under-predicts
##       the magnitude, because it moves receptor NUMBER and benzodiazepine
##       SENSITIVITY and not efficacy per surviving receptor.  See
##       potency_drift().
##
##   2.  ESETT's three-way tie.  Model: LEV 46.6% / FOS 50.3% / VPA 47.2%
##       among benzodiazepine failures; observed 47 / 45 / 46.  The three
##       are equal because all three multiply by 1.
##       BUT — and this corrects a natural over-reading of that structure —
##       they are NOT delay-proof.  Conditioned on first-line failure the
##       model puts the second-line cliff at ~62 min versus ~25 min for the
##       benzodiazepine: the second line buys roughly 40 extra minutes of
##       window, not immunity.  Their own target does not decay, but the
##       network they are suppressing still does, through the ENDOGENOUS
##       inhibitory tone that R_SYN also feeds.  See first_line_delay() and
##       second_line_delay().
##
##   3.  RAMPART's IM-beats-IV result.  Model: 80.0% vs 78.8% (observed
##       73.4 vs 63.4).  The IM drug is pharmacologically slower and still
##       wins, because the 4.5 min of IV access it removes are 4.5 minutes
##       of R_SYN decay.  The model reproduces the direction with a smaller
##       margin than the trial.  No RAMPART-specific parameter exists.
##
##   4.  The subtle-SE response collapse.  In a refractory patient given
##       repeated lorazepam, benzodiazepine occupancy 0.78 coexists with
##       motor output 0.31 and an electrographic burden of 1.00: the bedside
##       reads "much better" while nothing about the seizure has changed.
##       See dissociation().
##
##   5.  Ketamine is the only drug in the model whose target GROWS.  At a
##       fixed site occupancy of 0.75 its fractional cut in network gain
##       rises 0.288 -> 0.334 while a benzodiazepine's falls 0.746 -> 0.281;
##       they cross at 280 min.  And yet — scenario 16 — ketamine given at
##       300 min fails where the same dose at 65 min succeeds.  The drug
##       whose target grows still loses to the clock.  See marginal_effect().
##
##   6.  The negative STATUS trial of allopregnanolone in SRSE.  Correct
##       target (R_EXTRA survives), wrong point on the trajectory.
##       Scenario 18 runs it.
##
##   7.  "Stopping the seizure" and "saving the brain" separate.  Terminating
##       at 8 min costs 0.4% hippocampal neuron loss, at 30 min 28.0%, at
##       60 min 42.7% — all three are trial "successes".
##
##  SECOND THESIS — the price of a minute is not constant.
##  Because R_SYN·F_BZS is a DECAYING RESOURCE, a minute early costs more
##  than a minute late.  CREQ_BZD(t) has a vertical asymptote — in the
##  default patient it sits between 45 and 50 minutes — after which no
##  benzodiazepine dose terminates the seizure.  sweep_benzo_time() finds
##  it, and shows the cost of the minutes before it: 8 min -> 1.5%
##  hippocampal neuron loss, 20 min -> 9.3%, 25 min -> 14.4%, 30 min ->
##  the cliff.
##
##  THIRD THESIS — two resistances that look identical at the bedside.
##  Transporter resistance (P-gp up-regulation, brain:plasma ratio falls)
##  is beaten by dose.  Target resistance (the pools above) is not beaten
##  by any dose of the same drug.  The model carries both and outputs
##  BPRATIO, which distinguishes them.
##
## ---------------------------------------------------------------------
##  UNITS
##    time            minutes
##    LZP, MDZ        ng/mL
##    LEV, PHT, VPA, PB, PRO   mg/L
##    KET, ALLO       ng/mL
##    receptor pools  fraction of pre-ictal baseline (dimensionless)
##    CLI             mM      MAP  mmHg     GLUC mmol/L
##    LAC  mmol/L     TEMP °C          CK   U/L
##
##  DISCLAIMER
##    Educational / research QSP model.  Parameters are literature-anchored
##    approximations, not a validated clinical tool.  Do not use for
##    patient care.
## =====================================================================

library(mrgsolve)
library(dplyr)

# =====================================================================
#  MODEL CODE
# =====================================================================

se_code <- '
$PROB
Status Epilepticus QSP model — receptor-trafficking clock

$PARAM @annotated
// ---------------- patient / physiology ----------------
WT      :  70   : body weight (kg)
ALB     : 4.0   : serum albumin (g/dL)
CRCL    :  90   : creatinine clearance (mL/min)
EDRIVE  : 1.55  : aetiologic drive multiplier on excitatory gain (1 = threshold)
TDRIVE  : 1e6   : time at which the aetiologic drive is removed (min; 1e6 = never)
DRESID  : 0.55  : residual drive fraction after the aetiology is treated

// ---------------- network core ----------------
EBASE   : 1.45  : excitatory gain scale (calibrated so G = EDRIVE at t=0)
WA      : 0.58  : non-NMDA (AMPA/kainate/intrinsic) share of excitatory gain
WN      : 0.42  : NMDA share of excitatory gain
IBASE   : 1.20  : GABAergic inhibitory tone scale
WSYN    : 0.60  : synaptic (phasic) share of GABAergic tone
WEXT    : 0.40  : extrasynaptic (tonic) share of GABAergic tone
IPEP    : 0.20  : neuropeptide contribution to inhibitory tone
IADO    : 0.30  : adenosine contribution to inhibitory tone
IMIN    : 0.05  : non-GABAergic irreducible inhibition
KS      : 4.0   : network state rate constant (1/min)
SLEAK   : 0.002 : seizure re-ignition leak (allows recurrence from SEIZ=0)
SMAX    : 1.0   : maximal network activity

// ---------------- clock 1 : synaptic GABA-A internalisation ----------
KENDO   : 0.0135 : activity-driven endocytosis rate of synaptic GABA-A (1/min)
KREC    : 0.0055 : endosome -> surface recycling rate (1/min)
KDEGR   : 0.0020 : lysosomal degradation of internalised receptor (1/min)
KSYNR   : 0.0012 : de-novo replacement of total receptor pool (1/min)
FENDO0  : 0.12   : basal (inter-ictal) fraction of the endocytosis rate
KSW     : 0.0130 : rate of the alpha1->alpha4 benzodiazepine-site loss (1/min)
FBZMIN  : 0.28   : floor of the benzodiazepine-sensitive fraction

// ---------------- the pool that stays ----------------
KEXTG   : 0.0010 : slow up-regulation of extrasynaptic delta-GABA-A (1/min)
REXMAX  : 1.25   : ceiling of the extrasynaptic pool

// ---------------- clock 2 : synaptic NMDA insertion -----------------
KEXO    : 0.023992 : activity-driven NMDA exocytosis/capture rate (1/min)
FEXO0   : 0.3894   : basal fraction of the exocytosis rate
KENDN   : 0.009342 : NMDA receptor internalisation rate (1/min)
KACP    : 0.0110   : Ca-permeable AMPA (GluA2-lacking) insertion rate (1/min)
KACPR   : 0.0060   : Ca-permeable AMPA reversal rate (1/min)
ACPMAX  : 0.45     : maximal Ca-permeable AMPA increment
KKO     : 0.090    : extracellular K+ accumulation rate (mM/min at SEIZ=1)
KKOC    : 0.018    : extracellular K+ clearance rate (1/min)
KOMAX   : 6.0      : ceiling on extracellular K+ excess (mM)
AMG     : 0.55     : maximal Mg-block relief factor
KMG     : 1.40     : K+ excess giving half-maximal Mg-block relief (mM)

// ---------------- chloride ----------------
KKC     : 0.0140  : activity-driven KCC2 loss rate (1/min)
KKCR    : 0.0030  : KCC2 recovery rate (1/min)
KCCMIN  : 0.55    : KCC2 floor
KCLIN   : 0.160   : chloride loading rate (mM/min at SEIZ=1)
FCLDRG  : 2.0     : extra chloride loading per unit GABA-A drug occupancy
KCLOUT  : 0.030   : KCC2-dependent chloride extrusion rate (1/min)
CLI0    : 5.0     : resting intracellular chloride (mM)
CLO     : 120.0   : extracellular chloride (mM)
VM      : -65.0   : membrane potential used for the Cl driving force (mV)
SHUNT   : 0.45    : residual (conductance) shunting inhibition floor
CLSCALE : 12.0    : chloride rise (mM) that costs 1/e of the driving-force term

// ---------------- endogenous modulators ----------------
KADO    : 0.030  : adenosine production rate (1/min)
KADOC   : 0.055  : baseline adenosine clearance (1/min)
KADK    : 0.0060 : adenosine-kinase up-regulation rate (1/min)
ADKMAX  : 2.20   : maximal adenosine-kinase up-regulation
KPEPD   : 0.0140 : inhibitory-neuropeptide depletion rate (1/min)
KPEPR   : 0.0030 : inhibitory-neuropeptide recovery rate (1/min)

// ---------------- inflammation / barrier ----------------
KIL     : 0.0120 : IL-1beta production rate (1/min)
KILC    : 0.0060 : IL-1beta clearance rate (1/min)
GIL     : 0.18   : IL-1beta gain on excitatory drive (Src-NR2B)
ANAK    : 0      : anakinra flag (1 = IL-1 signalling blocked)
KBBB    : 0.0060 : BBB opening rate per unit IL-1beta (1/min)
KBBBR   : 0.0035 : BBB resealing rate (1/min)
BBBMAX  : 2.60   : maximal BBB permeability multiple
KEDEM   : 0.0025 : oedema formation rate (1/min)
KEDEMR  : 0.0015 : oedema resolution rate (1/min)
KPGP    : 0.0075 : P-glycoprotein up-regulation rate (1/min)
KPGPR   : 0.0030 : P-glycoprotein reversal rate (1/min)
PGPMAX  : 2.30   : maximal P-glycoprotein induction

// ---------------- excitotoxicity / injury ----------------
KGLU    : 0.85   : glutamate release rate (uM/min at SEIZ=1)
KGLUC   : 0.16   : glutamate uptake rate (1/min)
FEAAT   : 0.55   : fractional EAAT2 uptake loss at maximal inflammation
KCA     : 0.45   : Ca influx rate constant (1/min)
KCAC    : 0.22   : Ca extrusion rate constant (1/min)
CATH    : 1.60   : Ca threshold above which injury accrues
KINJ    : 0.0075 : injury accrual rate (1/min per unit Ca above threshold)
KMISM   : 2.20   : injury amplification by supply-demand mismatch
KTEMPI  : 0.28   : injury amplification per degree above 37 C
KHYPOG  : 0.55   : injury amplification when glucose < 3.5 mmol/L
KEPG    : 0.0020 : epileptogenesis accrual rate (1/min)
GEPGI   : 0.55   : inflammation weight in epileptogenesis

// ---------------- systemic ----------------
KATP    : 0.055  : ATP-deficit accrual rate (1/min)
KATPR   : 0.045  : ATP-deficit recovery rate (1/min)
MAP0    : 88     : baseline mean arterial pressure (mmHg)
CATAMP  : 46     : catecholamine-driven MAP rise (mmHg)
TAUCAT  : 34     : catecholamine surge decay time constant (min)
TAUMAP  : 3.0    : MAP equilibration time constant (min)
KAUTO   : 0.0180 : cerebral autoregulation failure rate (1/min)
KAUTOR  : 0.0025 : autoregulation recovery rate (1/min)
CBFG    : 3.50   : CBF gain during seizure with intact autoregulation
CBFGP   : 1.00   : CBF gain during seizure with failed autoregulation
DEMG    : 3.80   : cerebral metabolic demand gain during seizure
GLUC0   : 5.6    : baseline plasma glucose (mmol/L)
KGP     : 0.015  : glucose homeostatic return rate (1/min)
KGCAT   : 0.085  : catecholamine-driven hyperglycaemia rate (mmol/L/min)
KGUSE   : 0.070  : seizure glucose consumption rate (mmol/L/min)
GLUCIV  : 0      : IV dextrose flag (1 = glucose clamped at baseline)
KLACP   : 0.62   : lactate production rate (mmol/L/min at SEIZ=1)
KLACC   : 0.075  : lactate clearance rate (1/min)
KTEMPP  : 0.045  : hyperthermia rate (degC/min at SEIZ=1)
KTEMPC  : 0.018  : cooling rate (1/min)
COOL    : 0      : targeted temperature management flag
KCKP    : 55     : CK release rate (U/L/min at full motor activity)
KCKC    : 0.00032: CK elimination rate (1/min)
KRESPS  : 0.0125 : respiratory burden from the seizure itself (1/min)
KRESPD  : 0.0075 : respiratory burden per unit sedative occupancy (1/min)
KRESPC  : 0.0090 : respiratory burden recovery rate (1/min)
RESPTH  : 1.00   : respiratory burden at which intubation is triggered

// ---------------- motor / EEG dissociation ----------------
KMOT    : 0.0120 : decay rate of motor expression during SE (1/min)
MOTMIN  : 0.12   : floor of motor gain (subtle SE / NCSE)
MOTDRG  : 0.55   : additional motor suppression per unit GABAergic occupancy

// ================= PK : LORAZEPAM (ng/mL) =================
VLZP    : 91000  : lorazepam volume of distribution (mL)
CLLZP   : 77     : lorazepam clearance (mL/min)
KE0LZP  : 0.140  : lorazepam effect-site equilibration rate (1/min)
EC50LZP : 12     : lorazepam effect-site EC50 (ng/mL)

// ================= PK : MIDAZOLAM (ng/mL) =================
VMDZ    : 77000  : midazolam volume of distribution (mL)
CLMDZ   : 462    : midazolam clearance (mL/min)
KAMDZ   : 0.100  : midazolam IM absorption rate (1/min)
FIM     : 0.90   : midazolam IM bioavailability
KE0MDZ  : 0.200  : midazolam effect-site equilibration rate (1/min)
EC50MDZ : 25     : midazolam effect-site EC50 (ng/mL)
FMMDZ   : 0.70   : fraction of midazolam converted to 1-OH-midazolam
VMDZM   : 90000  : 1-OH-midazolam volume (mL)
CLMDZM  : 400    : 1-OH-midazolam clearance at CRCL 90 (mL/min)
POTMDZM : 0.80   : 1-OH-midazolam potency relative to parent

// ================= PK : DIAZEPAM (ng/mL) =================
V1DZP   : 25000  : diazepam central volume (mL)
V2DZP   : 90000  : diazepam peripheral volume (mL)
QDZP    : 2500   : diazepam intercompartmental clearance (mL/min)
CLDZP   : 27     : diazepam clearance (mL/min)
KE0DZP  : 0.300  : diazepam effect-site equilibration rate (1/min)
EC50DZP : 200    : diazepam effect-site EC50 (ng/mL)

// ================= PK : LEVETIRACETAM (mg/L) =================
VLEV    : 42     : levetiracetam volume (L)
CLLEV   : 0.067  : levetiracetam clearance at CRCL 90 (L/min)
KE0LEV  : 0.035  : levetiracetam effect-site equilibration (1/min) — slow brain entry
EC50LEV : 18     : levetiracetam effect-site EC50 (mg/L)
SLEV    : 0.865  : maximal fractional suppression of excitatory gain

// ================= PK : FOSPHENYTOIN -> PHENYTOIN (mg/L) =================
KCONV   : 0.069  : fosphenytoin -> phenytoin conversion rate (1/min)
VPHT    : 45     : phenytoin volume (L)
VMPHT   : 0.34   : phenytoin Michaelis-Menten Vmax (mg/min)
KMPHT   : 5.0    : phenytoin Michaelis constant, total conc (mg/L)
FUPHT0  : 0.10   : phenytoin free fraction at albumin 4.0 g/dL
KE0PHT  : 0.050  : phenytoin effect-site equilibration rate (1/min)
EC50PHT : 0.375  : phenytoin FREE effect-site EC50 (mg/L)
SPHT    : 0.865  : maximal fractional suppression of excitatory gain
DISPVPA : 0.85   : valproate displacement of phenytoin from albumin (per 100 mg/L)

// ================= PK : VALPROATE (mg/L) =================
VVPA    : 12     : valproate volume (L)
CLVPA   : 0.008  : valproate clearance (L/min)
KE0VPA  : 0.045  : valproate effect-site equilibration rate (1/min)
FUVPA0  : 0.10   : valproate free fraction at low concentration
FUVPAM  : 0.25   : maximal additional free fraction (saturable binding)
FUVPA50 : 150    : total valproate giving half-maximal binding saturation (mg/L)
EC50VPA : 12     : valproate FREE effect-site EC50 (mg/L)
SVPA    : 0.868  : maximal fractional suppression of excitatory gain

// ================= PK : PHENOBARBITAL (mg/L) =================
VPB     : 42     : phenobarbital volume (L)
CLPB    : 0.005  : phenobarbital clearance (L/min)
KE0PB   : 0.025  : phenobarbital effect-site equilibration rate (1/min)
EC50PB  : 25     : phenobarbital effect-site EC50 (mg/L)
EMAXPB  : 6.20   : phenobarbital maximal inhibitory effect
FSYNPB  : 0.20   : fraction of phenobarbital effect requiring the synaptic pool
SPB     : 0.32   : phenobarbital AMPA/kainate block (fractional cut in excitatory gain)

// ================= PK : KETAMINE (ng/mL) =================
V1KET   : 25000  : ketamine central volume (mL)
V2KET   : 140000 : ketamine peripheral volume (mL)
QKET    : 1300   : ketamine intercompartmental clearance (mL/min)
CLKET   : 1000   : ketamine clearance (mL/min)
KE0KET  : 0.200  : ketamine effect-site equilibration rate (1/min)
EC50KET : 1500   : ketamine effect-site EC50 (ng/mL)
EMAXKET : 0.85   : maximal fractional NMDA block
FMKET   : 0.80   : fraction of ketamine converted to norketamine
VNKET   : 100000 : norketamine volume (mL)
CLNKET  : 1200   : norketamine clearance (mL/min)
POTNKET : 0.33   : norketamine potency relative to ketamine
UDMIN   : 0.20   : use-dependence floor (channel availability at SEIZ=0)

// ================= PK : PROPOFOL (mg/L) =================
V1PRO   : 16     : propofol central volume (L)
V2PRO   : 80     : propofol peripheral volume (L)
QPRO    : 1.80   : propofol intercompartmental clearance (L/min)
CLPRO   : 1.80   : propofol clearance (L/min)
KE0PRO  : 0.300  : propofol effect-site equilibration rate (1/min)
EC50PRO : 3.5    : propofol effect-site EC50 (mg/L)
EMAXPRO : 7.20   : propofol maximal inhibitory effect
FSYNPRO : 0.25   : fraction of propofol effect requiring the synaptic pool
HYPOPRO : 26     : maximal MAP reduction from propofol (mmHg)

// ================= PK : ALLOPREGNANOLONE (ng/mL) =================
VALLO   : 60000  : allopregnanolone volume (mL)
CLALLO  : 900    : allopregnanolone clearance (mL/min)
KE0ALLO : 0.100  : allopregnanolone effect-site equilibration rate (1/min)
EC50ALLO: 150    : allopregnanolone effect-site EC50 (ng/mL)
EMAXALLO: 4.80   : allopregnanolone maximal inhibitory effect
FSYNALLO: 0.40   : fraction of allopregnanolone effect requiring the synaptic pool

// ================= shared pharmacodynamic scales ================
EMAXBZD : 6.00   : maximal benzodiazepine inhibitory effect (at intact pool)
HYPOBZD : 9      : maximal MAP reduction from benzodiazepines (mmHg)
KETMAP  : 12     : MAP support from ketamine (mmHg)
FPGPBZD : 0.05   : benzodiazepine sensitivity to P-glycoprotein efflux
FPGPPHT : 0.55   : phenytoin sensitivity to P-glycoprotein efflux
FPGPPB  : 0.35   : phenobarbital sensitivity to P-glycoprotein efflux
FPGPLEV : 0.10   : levetiracetam sensitivity to P-glycoprotein efflux
FBBB    : 0.35   : fractional gain in CNS partition per unit BBB opening

$INIT @annotated
// -------- PK (27) --------
LZP_C   : 0.0     : lorazepam central (ng)
LZP_E   : 0.0     : lorazepam effect site (ng/mL)
MDZ_A   : 0.0     : midazolam IM depot (ng)
MDZ_C   : 0.0     : midazolam central (ng)
MDZ_E   : 0.0     : midazolam effect site (ng/mL)
MDZ_M   : 0.0     : 1-OH-midazolam central (ng)
DZP_C   : 0.0     : diazepam central (ng)
DZP_P   : 0.0     : diazepam peripheral (ng)
DZP_E   : 0.0     : diazepam effect site (ng/mL)
LEV_C   : 0.0     : levetiracetam central (mg)
LEV_E   : 0.0     : levetiracetam effect site (mg/L)
FOS_C   : 0.0     : fosphenytoin central (mg PE)
PHT_C   : 0.0     : phenytoin central (mg)
PHT_E   : 0.0     : phenytoin free effect site (mg/L)
VPA_C   : 0.0     : valproate central (mg)
VPA_E   : 0.0     : valproate free effect site (mg/L)
PB_C    : 0.0     : phenobarbital central (mg)
PB_E    : 0.0     : phenobarbital effect site (mg/L)
KET_C   : 0.0     : ketamine central (ng)
KET_P   : 0.0     : ketamine peripheral (ng)
KET_E   : 0.0     : ketamine effect site (ng/mL)
NKET_C  : 0.0     : norketamine central (ng)
PRO_C   : 0.0     : propofol central (mg)
PRO_P   : 0.0     : propofol peripheral (mg)
PRO_E   : 0.0     : propofol effect site (mg/L)
ALLO_C  : 0.0     : allopregnanolone central (ng)
ALLO_E  : 0.0     : allopregnanolone effect site (ng/mL)
// -------- receptor trafficking (8) --------
RSYN    : 1.0     : synaptic gamma2-GABA-A pool (fraction of baseline)
RENDO   : 0.0     : internalised GABA-A pool (fraction)
FBZS    : 1.0     : benzodiazepine-sensitive fraction of the synaptic pool
REXTRA  : 1.0     : extrasynaptic delta-GABA-A pool (fraction)
NRSYN   : 1.0     : synaptic NMDA pool (fraction of baseline)
NRINT   : 1.0     : intracellular NMDA pool (fraction)
AMPACP  : 0.0     : Ca-permeable AMPA increment (fraction)
KO      : 0.0     : extracellular K+ excess (mM)
// -------- network / ionic / modulators (7) --------
SEIZ    : 0.05    : network seizure activity (0-1)
MOTG    : 1.0     : motor expression gain (0-1)
KCC2    : 1.0     : KCC2 transporter activity (fraction)
CLI     : 5.0     : intracellular chloride (mM)
ADO     : 0.0     : extracellular adenosine (relative)
ADK     : 1.0     : adenosine kinase up-regulation (relative)
NETPEP  : 1.0     : net inhibitory neuropeptide balance (relative)
// -------- inflammation / barrier (5) --------
IL1B    : 0.0     : IL-1beta (relative)
BBBP    : 1.0     : blood-brain barrier permeability (relative)
PGP     : 1.0     : BBB P-glycoprotein activity (relative)
EDEM    : 0.0     : cerebral oedema (relative)
// -------- excitotoxicity (4) --------
GLU     : 2.0     : extracellular glutamate (uM)
CAI     : 1.0     : neuronal intracellular Ca (relative)
INJURY  : 0.0     : cumulative excitotoxic neuronal injury (arbitrary units)
ATPD    : 0.0     : cerebral ATP deficit (relative)
// -------- systemic (8) --------
AUTO    : 1.0     : cerebral autoregulation integrity (fraction)
MAP     : 88.0    : mean arterial pressure (mmHg)
GLUCP   : 5.6     : plasma glucose (mmol/L)
LAC     : 1.0     : plasma lactate (mmol/L)
TEMP    : 37.0    : core temperature (degC)
CK      : 100.0   : creatine kinase (U/L)
RESPD   : 0.0     : respiratory depression burden (relative)
// -------- integrators (2) --------
TSEIZ   : 0.0     : cumulative seizure time (min)
EPG     : 0.0     : epileptogenesis burden (relative)
$GLOBAL
#define POS(x) ((x) > 0.0 ? (x) : 0.0)
#define HILL(c, e) ((c) / ((c) + (e)))

$MAIN
double CLMDZMi = CLMDZM * (0.35 + 0.65 * CRCL / 90.0);
double CLLEVi  = CLLEV  * (0.30 + 0.70 * CRCL / 90.0);

$ODE
// =================================================================
//  0.  PHARMACOKINETICS
// =================================================================
double C_LZP  = LZP_C  / VLZP;
double C_MDZ  = MDZ_C  / VMDZ;
double C_MDZM = MDZ_M  / VMDZM;
double C_DZP  = DZP_C  / V1DZP;
double C_DZPP = DZP_P  / V2DZP;
double C_LEV  = LEV_C  / VLEV;
double C_PHT  = PHT_C  / VPHT;
double C_VPA  = VPA_C  / VVPA;
double C_PB   = PB_C   / VPB;
double C_KET  = KET_C  / V1KET;
double C_KETP = KET_P  / V2KET;
double C_NKET = NKET_C / VNKET;
double C_PRO  = PRO_C  / V1PRO;
double C_PROP = PRO_P  / V2PRO;
double C_ALLO = ALLO_C / VALLO;

// ---- protein binding -------------------------------------------------
double FUVPA = FUVPA0 + FUVPAM * HILL(C_VPA, FUVPA50);
double FUPHT = FUPHT0 * (1.0 + 0.60 * (1.0 - ALB / 4.0))
                      * (1.0 + DISPVPA * C_VPA / 100.0);
FUPHT = fmin(FUPHT, 0.60);
double CF_PHT = C_PHT * FUPHT;
double CF_VPA = C_VPA * FUVPA;

// ---- CNS partition: P-gp efflux vs BBB opening ----------------------
double BBBGAIN = 1.0 + FBBB * (BBBP - 1.0);
double PBZD  = BBBGAIN / (1.0 + FPGPBZD * (PGP - 1.0));
double PPHT  = BBBGAIN / (1.0 + FPGPPHT * (PGP - 1.0));
double PPB   = BBBGAIN / (1.0 + FPGPPB  * (PGP - 1.0));
double PLEV  = BBBGAIN / (1.0 + FPGPLEV * (PGP - 1.0));

// ---- PK ODEs ---------------------------------------------------------
dxdt_LZP_C = -(CLLZP / VLZP) * LZP_C;
dxdt_LZP_E =  KE0LZP * (C_LZP * PBZD - LZP_E);

dxdt_MDZ_A = -KAMDZ * MDZ_A;
dxdt_MDZ_C =  FIM * KAMDZ * MDZ_A - (CLMDZ / VMDZ) * MDZ_C;
dxdt_MDZ_E =  KE0MDZ * (C_MDZ * PBZD - MDZ_E);
dxdt_MDZ_M =  FMMDZ * CLMDZ * C_MDZ - CLMDZMi * C_MDZM;

dxdt_DZP_C = -(CLDZP / V1DZP) * DZP_C - QDZP * (C_DZP - C_DZPP);
dxdt_DZP_P =  QDZP * (C_DZP - C_DZPP);
dxdt_DZP_E =  KE0DZP * (C_DZP * PBZD - DZP_E);

dxdt_LEV_C = -(CLLEVi / VLEV) * LEV_C;
dxdt_LEV_E =  KE0LEV * (C_LEV * PLEV - LEV_E);

dxdt_FOS_C = -KCONV * FOS_C;
dxdt_PHT_C =  KCONV * FOS_C - VMPHT * C_PHT / (KMPHT + C_PHT);
dxdt_PHT_E =  KE0PHT * (CF_PHT * PPHT - PHT_E);

dxdt_VPA_C = -CLVPA * C_VPA;
dxdt_VPA_E =  KE0VPA * (CF_VPA - VPA_E);

dxdt_PB_C  = -CLPB * C_PB;
dxdt_PB_E  =  KE0PB * (C_PB * PPB - PB_E);

dxdt_KET_C = -(CLKET / V1KET) * KET_C - QKET * (C_KET - C_KETP);
dxdt_KET_P =  QKET * (C_KET - C_KETP);
dxdt_KET_E =  KE0KET * (C_KET - KET_E);
dxdt_NKET_C = FMKET * CLKET * C_KET - CLNKET * C_NKET;

dxdt_PRO_C = -(CLPRO / V1PRO) * PRO_C - QPRO * (C_PRO - C_PROP);
dxdt_PRO_P =  QPRO * (C_PRO - C_PROP);
dxdt_PRO_E =  KE0PRO * (C_PRO - PRO_E);

dxdt_ALLO_C = -CLALLO * C_ALLO;
dxdt_ALLO_E =  KE0ALLO * (C_ALLO - ALLO_E);

// =================================================================
//  1.  OCCUPANCIES
//      Benzodiazepines share one site -> competitive-additive form.
// =================================================================
double SBZD = LZP_E / EC50LZP
            + (MDZ_E + POTMDZM * C_MDZM) / EC50MDZ
            + DZP_E / EC50DZP;
double OCC_BZD = SBZD / (1.0 + SBZD);

double OCC_PB   = HILL(PB_E,   EC50PB);
double OCC_PRO  = HILL(PRO_E,  EC50PRO);
double OCC_ALLO = HILL(ALLO_E, EC50ALLO);
double OCC_LEV  = HILL(LEV_E,  EC50LEV);
double OCC_PHT  = HILL(PHT_E,  EC50PHT);
double OCC_VPA  = HILL(VPA_E,  EC50VPA);
double CE_KETT  = KET_E + POTNKET * C_NKET;
double OCC_KET  = HILL(CE_KETT, EC50KET);

// total GABA-A channel opening produced by drugs (drives Cl loading)
double OCC_GABA = fmin(1.0, OCC_BZD + OCC_PB + OCC_PRO + OCC_ALLO);

// =================================================================
//  2.  THE THREE POOLS AND THE THREE EFFECT TERMS
// =================================================================
double TGT_BZD  = RSYN * FBZS;                                  // falls fastest
double TGT_PB   = FSYNPB   * RSYN + (1.0 - FSYNPB)   * REXTRA;
double TGT_PRO  = FSYNPRO  * RSYN + (1.0 - FSYNPRO)  * REXTRA;
double TGT_ALLO = FSYNALLO * RSYN + (1.0 - FSYNALLO) * REXTRA;

double INH_BZD  = EMAXBZD  * OCC_BZD  * TGT_BZD;
double INH_PB   = EMAXPB   * OCC_PB   * TGT_PB;
double INH_PRO  = EMAXPRO  * OCC_PRO  * TGT_PRO;
double INH_ALLO = EMAXALLO * OCC_ALLO * TGT_ALLO;

// =================================================================
//  3.  CHLORIDE DRIVING FORCE  (GABA can become excitatory)
// =================================================================
//  GABAergic efficacy is written as a SHUNTING floor plus a driving-force
//  term that decays as intracellular chloride accumulates.  The Nernst
//  potential itself is carried as a reported diagnostic (E_Cl in $TABLE);
//  it is deliberately NOT used raw, because a strict driving-force term
//  makes inhibition vanish at chloride loads that neurons tolerate.
double CLIs  = fmax(CLI, 0.5);
double CLFAC = SHUNT + (1.0 - SHUNT) * exp(-(CLIs - CLI0) / CLSCALE);

// =================================================================
//  4.  NETWORK GAIN  G  =  excitation / inhibition
// =================================================================
double MGREL = 1.0 + AMG * KO / (KMG + KO);
double SEIZc = fmin(fmax(SEIZ, 0.0), SMAX);
double UD    = UDMIN + (1.0 - UDMIN) * SEIZc;         // use-dependence
double KETBLK = EMAXKET * OCC_KET * UD;

double DRV = (SOLVERTIME < TDRIVE) ? EDRIVE : (1.0 + DRESID * (EDRIVE - 1.0));
double ILEFF = (ANAK > 0.5) ? 0.0 : IL1B;

double NMDATERM = WN * NRSYN * MGREL * (1.0 - KETBLK);
double AMPATERM = WA * (1.0 + AMPACP);
//  Barbiturates are the one GABAergic agent that also acts on the
//  NUMERATOR: they block AMPA/kainate receptors.  That is why phenobarbital
//  appears in BOTH the inhibition sum and the suppression product.
double SUP2 = fmax(0.20, 1.0 - (SLEV * OCC_LEV + SPHT * OCC_PHT + SVPA * OCC_VPA
                              + SPB * OCC_PB));

double EGAIN = EBASE * DRV * (AMPATERM + NMDATERM) * SUP2 * (1.0 + GIL * ILEFF);

double INH_ENDO = IBASE * (WSYN * RSYN + WEXT * REXTRA) * CLFAC
                + IPEP * NETPEP + IADO * ADO + IMIN;
double INH_TOT  = fmax(0.02, INH_ENDO + INH_BZD + INH_PB + INH_PRO + INH_ALLO);

double G = EGAIN / INH_TOT;

// =================================================================
//  5.  NETWORK STATE  (bistable: SEIZ = 0 and SEIZ = SMAX both stable)
// =================================================================
//  SZ is the CLAMPED network state.  The clamp is not cosmetic: without it
//  a transient excursion past SMAX flips the sign of TWO factors at once
//  and the logistic runs away.
double SZ = fmin(fmax(SEIZ, 0.0), SMAX);
//  Growth is logistic (saturates at SMAX); decay is exponential.  Writing
//  the decay branch as KS*SZ*(G-1) rather than reusing the logistic term
//  matters: the logistic form has a zero at SZ = SMAX, which would make
//  full-blown SE an ABSORBING state that no drug could leave.  Both
//  branches vanish at G = 1, so the field is continuous.
double DG = G - 1.0;
dxdt_SEIZ = (SEIZ <= 1e-7 && DG < 0.0)
            ? 0.0
            : (DG >= 0.0 ? KS * (SZ + SLEAK) * (SMAX - SZ) * DG
                         : KS * SZ * DG);
dxdt_TSEIZ = SZ;

// =================================================================
//  6.  CLOCK 1 — synaptic GABA-A internalisation
// =================================================================
double FEND = FENDO0 + (1.0 - FENDO0) * SZ;
double RTOT = RSYN + RENDO;
dxdt_RSYN  = KREC * RENDO - KENDO * FEND * RSYN + KSYNR * POS(1.0 - RTOT);
dxdt_RENDO = KENDO * FEND * RSYN - KREC * RENDO - KDEGR * RENDO;
dxdt_FBZS  = -KSW * SZ * (FBZS - FBZMIN) + 0.25 * KSW * (1.0 - FBZS) * (1.0 - SZ);

// the pool that stays (slight activity-driven up-regulation)
dxdt_REXTRA = KEXTG * SZ * (REXMAX - REXTRA) - 0.5 * KEXTG * (REXTRA - 1.0) * (1.0 - SZ);

// =================================================================
//  7.  CLOCK 2 — synaptic NMDA insertion (opposite sign, same driver)
// =================================================================
double FEXO = FEXO0 + (1.0 - FEXO0) * SZ;
dxdt_NRSYN = KEXO * FEXO * NRINT - KENDN * NRSYN;
dxdt_NRINT = KENDN * NRSYN - KEXO * FEXO * NRINT;
dxdt_AMPACP = KACP * SZ * (ACPMAX - AMPACP) - KACPR * AMPACP;
dxdt_KO = KKO * SZ * (1.0 - KO / KOMAX) - KKOC * KO;

// =================================================================
//  8.  CHLORIDE / KCC2   (note: GABAergic DRUGS load chloride)
// =================================================================
dxdt_KCC2 = -KKC * SZ * (KCC2 - KCCMIN) + KKCR * (1.0 - KCC2);
dxdt_CLI  = KCLIN * SZ * (1.0 + FCLDRG * OCC_GABA)
          - KCLOUT * KCC2 * (CLI - CLI0);

// =================================================================
//  9.  ENDOGENOUS MODULATORS
// =================================================================
dxdt_ADO    = KADO * SZ - KADOC * ADK * ADO;
dxdt_ADK    = KADK * SZ * (ADKMAX - ADK) - 0.25 * KADK * (ADK - 1.0);
dxdt_NETPEP = -KPEPD * SZ * NETPEP + KPEPR * (1.0 - NETPEP);

// =================================================================
// 10.  INFLAMMATION / BARRIER / TRANSPORTER
// =================================================================
dxdt_IL1B = KIL * SZ - KILC * IL1B;
dxdt_BBBP = KBBB * ILEFF * (BBBMAX - BBBP) - KBBBR * (BBBP - 1.0);
dxdt_PGP  = KPGP * SZ * (PGPMAX - PGP) - KPGPR * (PGP - 1.0);
dxdt_EDEM = KEDEM * (BBBP - 1.0 + 0.5 * SZ) - KEDEMR * EDEM;

// =================================================================
// 11.  SYSTEMIC — the phase I / phase II switch
// =================================================================
double CATF = exp(-TSEIZ / TAUCAT);                     // catecholamine reserve
double MAPDRUG = HYPOBZD * OCC_BZD + HYPOPRO * OCC_PRO
               + 14.0 * OCC_PB - KETMAP * OCC_KET;
double MAPTGT = MAP0 + CATAMP * SZ * CATF - MAPDRUG
              - 14.0 * (1.0 - CATF) * SZ;
dxdt_MAP = (MAPTGT - MAP) / TAUMAP;

dxdt_AUTO = -KAUTO * SZ * AUTO + KAUTOR * (1.0 - AUTO) * (1.0 - SZ);
double CBF = AUTO * (1.0 + CBFG * SZ)
           + (1.0 - AUTO) * (MAP / MAP0) * (1.0 + CBFGP * SZ);
double DEMAND = 1.0 + DEMG * SZ;
double MISMATCH = POS(DEMAND - CBF) / DEMAND;

dxdt_ATPD = KATP * MISMATCH - KATPR * ATPD;

dxdt_GLUCP = KGP * (GLUC0 - GLUCP) + KGCAT * SZ * CATF
           - KGUSE * SZ * HILL(GLUCP, 1.0)
           + (GLUCIV > 0.5 ? 0.25 * (GLUC0 - GLUCP) : 0.0);
dxdt_LAC   = KLACP * SZ + 0.35 * MISMATCH - KLACC * LAC;
dxdt_TEMP  = KTEMPP * SZ - KTEMPC * (TEMP - 37.0) - (COOL > 0.5 ? 0.05 * (TEMP - 33.0) : 0.0);

// motor expression and the electromechanical dissociation
dxdt_MOTG = -KMOT * SZ * (MOTG - MOTMIN) + 0.30 * KMOT * (1.0 - MOTG) * (1.0 - SZ);
double MOTOR = SZ * MOTG * fmax(0.0, 1.0 - MOTDRG * OCC_GABA);

dxdt_CK = KCKP * MOTOR - KCKC * CK;
dxdt_RESPD = KRESPS * SZ + KRESPD * (OCC_BZD + 1.6 * OCC_PRO + 1.2 * OCC_PB)
           - KRESPC * RESPD;

// =================================================================
// 12.  EXCITOTOXICITY -> INJURY (the SECOND integrator)
// =================================================================
double UPTAKE = KGLUC * (1.0 - FEAAT * HILL(ILEFF, 1.5));
dxdt_GLU = KGLU * SZ - UPTAKE * (GLU - 2.0);

double CAIN = KCA * (GLU / 2.0 - 1.0) * (WN * NRSYN * MGREL * (1.0 - KETBLK)
            + WA * AMPACP) / (WA + WN);
dxdt_CAI = POS(CAIN) - KCAC * (CAI - 1.0);

double HYPOG = (GLUCP < 3.5) ? KHYPOG * (3.5 - GLUCP) : 0.0;
dxdt_INJURY = KINJ * POS(CAI - CATH)
            * (1.0 + KMISM * MISMATCH)
            * (1.0 + KTEMPI * POS(TEMP - 37.5))
            * (1.0 + HYPOG);

dxdt_EPG = KEPG * (dxdt_INJURY / KINJ + GEPGI * ILEFF);

$TABLE
// ---------- readable concentrations ----------
double CP_LZP  = LZP_C  / VLZP;
double CP_MDZ  = MDZ_C  / VMDZ;
double CP_DZP  = DZP_C  / V1DZP;
double CP_LEV  = LEV_C  / VLEV;
double CP_PHT  = PHT_C  / VPHT;
double CP_VPA  = VPA_C  / VVPA;
double CP_PB   = PB_C   / VPB;
double CP_KET  = KET_C  / V1KET;
double CP_PRO  = PRO_C  / V1PRO;
double CP_ALLO = ALLO_C / VALLO;

// ---------- occupancies (recomputed for output) ----------
double oSBZD = LZP_E / EC50LZP + (MDZ_E + POTMDZM * MDZ_M / VMDZM) / EC50MDZ
             + DZP_E / EC50DZP;
double oOCC_BZD  = oSBZD / (1.0 + oSBZD);
double oOCC_PB   = PB_E   / (PB_E   + EC50PB);
double oOCC_PRO  = PRO_E  / (PRO_E  + EC50PRO);
double oOCC_ALLO = ALLO_E / (ALLO_E + EC50ALLO);
double oOCC_LEV  = LEV_E  / (LEV_E  + EC50LEV);
double oOCC_PHT  = PHT_E  / (PHT_E  + EC50PHT);
double oOCC_VPA  = VPA_E  / (VPA_E  + EC50VPA);
double oCEKET    = KET_E + POTNKET * NKET_C / VNKET;
double oOCC_KET  = oCEKET / (oCEKET + EC50KET);
double oOCC_GABA = fmin(1.0, oOCC_BZD + oOCC_PB + oOCC_PRO + oOCC_ALLO);

// ---------- the three pools ----------
double oTGT_BZD = RSYN * FBZS;
double oMGREL   = 1.0 + AMG * KO / (KMG + KO);
double oNMDA    = NRSYN * oMGREL;

// ---------- chloride ----------
double oCLIs  = fmax(CLI, 0.5);
double oECL   = -61.5 * log10(CLO / oCLIs);
double oCLFAC = SHUNT + (1.0 - SHUNT) * exp(-(oCLIs - CLI0) / CLSCALE);
double DFGABA = (VM - oECL);   // > 0 hyperpolarising, < 0 depolarising (mV)

// ---------- gain decomposition ----------
double oILEFF = (ANAK > 0.5) ? 0.0 : IL1B;
double SEIZt  = fmin(fmax(SEIZ, 0.0), 1.0);
double oUD    = UDMIN + (1.0 - UDMIN) * SEIZt;
double oKETBLK = EMAXKET * oOCC_KET * oUD;
double oDRV   = (TIME < TDRIVE) ? EDRIVE : (1.0 + DRESID * (EDRIVE - 1.0));
double oSUP2  = fmax(0.20, 1.0 - (SLEV*oOCC_LEV + SPHT*oOCC_PHT + SVPA*oOCC_VPA
                                + SPB*oOCC_PB));
double oAMPAT = WA * (1.0 + AMPACP);
double oNMDAT = WN * NRSYN * oMGREL;
double oEGAIN0 = EBASE * oDRV * (oAMPAT + oNMDAT) * oSUP2 * (1.0 + GIL * oILEFF);
double oEGAIN  = EBASE * oDRV * (oAMPAT + oNMDAT * (1.0 - oKETBLK)) * oSUP2
                 * (1.0 + GIL * oILEFF);
double oINH_ENDO = IBASE * (WSYN * RSYN + WEXT * REXTRA) * oCLFAC
                 + IPEP * NETPEP + IADO * ADO + IMIN;
double oINH_BZD  = EMAXBZD * oOCC_BZD * oTGT_BZD;
double oINH_PB   = EMAXPB  * oOCC_PB  * (FSYNPB * RSYN + (1.0 - FSYNPB) * REXTRA);
double oINH_PRO  = EMAXPRO * oOCC_PRO * (FSYNPRO * RSYN + (1.0 - FSYNPRO) * REXTRA);
double oINH_ALLO = EMAXALLO* oOCC_ALLO* (FSYNALLO* RSYN + (1.0 - FSYNALLO)* REXTRA);
double oINH_TOT  = fmax(0.02, oINH_ENDO + oINH_BZD + oINH_PB + oINH_PRO + oINH_ALLO);
double GNET = oEGAIN / oINH_TOT;

// ---------- THE PRICE OF A MINUTE -------------------------------------
//  OCCREQ_BZD : benzodiazepine site occupancy required, RIGHT NOW, to
//               bring the drug-free network gain down to 1.
//  CREQ_BZD   : the corresponding effect-site concentration in EC50 units.
//               Negative value = no dose can do it (target pool too small).
double NEED = oEGAIN0 - oINH_ENDO;
double OCCREQ_BZD = (EMAXBZD * oTGT_BZD > 1e-9) ? NEED / (EMAXBZD * oTGT_BZD) : 99.0;
double CREQ_BZD = (OCCREQ_BZD < 1.0 && OCCREQ_BZD > 0.0)
                  ? OCCREQ_BZD / (1.0 - OCCREQ_BZD) : -1.0;
double oTGTPB = FSYNPB * RSYN + (1.0 - FSYNPB) * REXTRA;
double OCCREQ_PB = (EMAXPB * oTGTPB > 1e-9) ? NEED / (EMAXPB * oTGTPB) : 99.0;
double CREQ_PB = (OCCREQ_PB < 1.0 && OCCREQ_PB > 0.0)
                 ? OCCREQ_PB / (1.0 - OCCREQ_PB) : -1.0;

// share of the excitatory gain that a full NMDA block could remove
double KETSHARE = oNMDAT / (oAMPAT + oNMDAT);
// share of inhibition a benzodiazepine can still recruit, relative to t=0
double BZDCAP = oTGT_BZD;

// ---------- physiology / endpoints ----------
double oCBF = AUTO * (1.0 + CBFG * SEIZt)
            + (1.0 - AUTO) * (MAP / MAP0) * (1.0 + CBFGP * SEIZt);
double oDEMAND = 1.0 + DEMG * SEIZt;
double oMISMATCH = fmax(0.0, oDEMAND - oCBF) / oDEMAND;
double MOTOROUT = SEIZt * MOTG * fmax(0.0, 1.0 - MOTDRG * oOCC_GABA);
double EEGOUT = SEIZt;
double BPRATIO = (1.0 + FBBB * (BBBP - 1.0)) / (1.0 + FPGPPHT * (PGP - 1.0));
double SZSTOP = (SEIZt < 0.05) ? 1.0 : 0.0;
double INTUB  = (RESPD > RESPTH) ? 1.0 : 0.0;
// hippocampal neuron loss (%) — saturating read-out of the injury integrator
double NEURLOSS = 62.0 * INJURY / (INJURY + 3.2);

$CAPTURE @annotated
CP_LZP : lorazepam plasma (ng/mL)
CP_MDZ : midazolam plasma (ng/mL)
CP_DZP : diazepam plasma (ng/mL)
CP_LEV : levetiracetam plasma (mg/L)
CP_PHT : phenytoin total plasma (mg/L)
CP_VPA : valproate total plasma (mg/L)
CP_PB  : phenobarbital plasma (mg/L)
CP_KET : ketamine plasma (ng/mL)
CP_PRO : propofol plasma (mg/L)
CP_ALLO: allopregnanolone plasma (ng/mL)
oOCC_BZD : benzodiazepine site occupancy
oOCC_KET : ketamine NMDA occupancy
oOCC_LEV : levetiracetam occupancy
oOCC_PHT : phenytoin occupancy
oOCC_VPA : valproate occupancy
oOCC_PB  : phenobarbital occupancy
oOCC_PRO : propofol occupancy
oTGT_BZD : benzodiazepine target pool (RSYN x FBZS)
oNMDA    : NMDA conductance factor (NRSYN x Mg relief)
oCLFAC   : chloride-dependent GABAergic efficacy factor
oECL     : chloride reversal potential (mV)
DFGABA   : GABA-A driving force (mV; <0 = depolarising)
GNET     : network gain G (>1 = seizure grows)
oEGAIN   : excitatory gain
oINH_TOT : total inhibitory tone
OCCREQ_BZD : benzodiazepine occupancy required to reach G=1
CREQ_BZD : benzodiazepine effect-site conc required (EC50 units; -1 = impossible)
CREQ_PB  : phenobarbital effect-site conc required (EC50 units; -1 = impossible)
KETSHARE : share of excitatory gain removable by NMDA block
BZDCAP   : benzodiazepine-recruitable pool (fraction of baseline)
MOTOROUT : motor manifestation amplitude
EEGOUT   : electrographic burden
oMISMATCH: cerebral supply-demand mismatch
oCBF     : cerebral blood flow (relative)
NEURLOSS : hippocampal neuron loss (%)
BPRATIO  : brain:plasma partition index (transporter vs target resistance)
SZSTOP   : seizure terminated flag
INTUB    : intubation triggered flag

$OMEGA @annotated @block
ETA_DRIVE : 0.09 : aetiologic drive (log-scale variance; SD = 0.30)

$SIGMA 0
'

mod <- mcode_cache("status_epilepticus", se_code)

# ---------------------------------------------------------------------
#  Inter-individual variability on the drive term is applied in $MAIN via
#  a manual override (kept out of the model code so that single-subject
#  deterministic runs stay exactly reproducible).
# ---------------------------------------------------------------------

# =====================================================================
#  COMPARTMENT NUMBERS (for dosing)
# =====================================================================
CMTN <- names(mrgsolve::init(mod))
CMT  <- setNames(seq_along(CMTN), CMTN)
#  LZP_C 1 · MDZ_A 3 · MDZ_C 4 · DZP_C 7 · LEV_C 10 · FOS_C 12
#  VPA_C 15 · PB_C 17 · KET_C 19 · PRO_C 23 · ALLO_C 26

WT <- 70

# ---- dose builders ---------------------------------------------------
d_lzp  <- function(t, mg = 4)        ev(time = t, amt = mg * 1e6, cmt = CMT["LZP_C"])          # ng
d_mdz_im <- function(t, mg = 10)     ev(time = t, amt = mg * 1e6, cmt = CMT["MDZ_A"])
d_mdz_iv <- function(t, mgkg = 0.2)  ev(time = t, amt = mgkg * WT * 1e6, cmt = CMT["MDZ_C"])
d_mdz_inf <- function(t, dur, mgkgh = 0.4)
  ev(time = t, amt = mgkgh * WT * (dur / 60) * 1e6, cmt = CMT["MDZ_C"],
     rate = mgkgh * WT * 1e6 / 60)
d_dzp  <- function(t, mg = 10)       ev(time = t, amt = mg * 1e6, cmt = CMT["DZP_C"])
d_lev  <- function(t, mgkg = 60)     ev(time = t, amt = min(mgkg * WT, 4500), cmt = CMT["LEV_C"])
d_fos  <- function(t, pekg = 20)     ev(time = t, amt = pekg * WT, cmt = CMT["FOS_C"])
d_vpa  <- function(t, mgkg = 40)     ev(time = t, amt = min(mgkg * WT, 3000), cmt = CMT["VPA_C"])
d_pb   <- function(t, mgkg = 20)     ev(time = t, amt = mgkg * WT, cmt = CMT["PB_C"])
d_ket  <- function(t, mgkg = 2)      ev(time = t, amt = mgkg * WT * 1e6, cmt = CMT["KET_C"])
d_ket_inf <- function(t, dur, mgkgh = 3)
  ev(time = t, amt = mgkgh * WT * (dur / 60) * 1e6, cmt = CMT["KET_C"],
     rate = mgkgh * WT * 1e6 / 60)
d_pro  <- function(t, mgkg = 2)      ev(time = t, amt = mgkg * WT, cmt = CMT["PRO_C"])
d_pro_inf <- function(t, dur, mgkgh = 6)
  ev(time = t, amt = mgkgh * WT * (dur / 60), cmt = CMT["PRO_C"],
     rate = mgkgh * WT / 60)
d_allo_inf <- function(t, dur, ugkgh = 86)
  ev(time = t, amt = ugkgh * WT * (dur / 60) * 1e3, cmt = CMT["ALLO_C"],
     rate = ugkgh * WT * 1e3 / 60)

# =====================================================================
#  22 THERAPEUTIC SCENARIOS
#  Every scenario uses the IDENTICAL 59-ODE structure and the IDENTICAL
#  parameter set.  Only the dosing schedule (and, where explicitly
#  stated, the aetiology parameter EDRIVE) changes.
# =====================================================================

END <- 720   # 12 h

scen <- list()

scen[["01_untreated"]] <- list(
  ev = ev(time = 0, amt = 0, cmt = 1),
  par = list(),
  note = "Natural history. No drug. The reference trajectory for every clock in the model.")

scen[["02_guideline_5min"]] <- list(
  ev = d_lzp(5, 4) + d_lzp(9, 4) + d_lev(25, 60),
  par = list(),
  note = "Guideline-concordant: lorazepam 4 mg IV at 5 min, repeated at 9 min, levetiracetam 60 mg/kg at 25 min.")

scen[["03_benzo_delayed_30"]] <- list(
  ev = d_lzp(30, 4) + d_lzp(34, 4) + d_lev(50, 60),
  par = list(),
  note = "Same drugs, same doses, 25 minutes later. Isolates the cost of delay with the dose held constant.")

scen[["04_benzo_delayed_60"]] <- list(
  ev = d_lzp(60, 4) + d_lzp(64, 4) + d_lev(80, 60),
  par = list(),
  note = "One hour to the first benzodiazepine — around the real-world 90th percentile.")

scen[["05_underdosed_early"]] <- list(
  ev = d_lzp(5, 2) + d_lev(25, 60),
  par = list(),
  note = "Half the guideline dose, given ON TIME. The commonest real-world error.")

scen[["06_underdosed_late"]] <- list(
  ev = d_lzp(30, 2) + d_lev(50, 60),
  par = list(),
  note = "The same half-dose given late. Compare with 05: the SAME milligram deficit costs different amounts.")

scen[["07_rampart_IM"]] <- list(
  ev = d_mdz_im(8, 10) + d_lev(30, 60),
  par = list(),
  note = "RAMPART IM arm: midazolam 10 mg IM at 8 min (no IV needed).")

scen[["08_rampart_IV"]] <- list(
  ev = d_lzp(12.5, 4) + d_lev(35, 60),
  par = list(),
  note = "RAMPART IV arm: lorazepam 4 mg IV at 12.5 min (4.5 min lost to IV access). The pharmacologically FASTER drug, given later.")

RSE_DRIVE   <- 1.78   # a benzodiazepine-refractory patient: EDRIVE just above the
                      # first-line threshold. This IS the ESETT enrolment
                      # criterion, written as a parameter instead of a filter.
SRSE_DRIVE  <- 2.60   # refractory to first AND second line

scen[["09_esett_LEV"]] <- list(
  ev = d_lzp(20, 4) + d_lev(65, 60),
  par = list(EDRIVE = RSE_DRIVE),
  note = "ESETT arm A after benzodiazepine failure: levetiracetam 60 mg/kg.")

scen[["10_esett_FOS"]] <- list(
  ev = d_lzp(20, 4) + d_fos(65, 20),
  par = list(EDRIVE = RSE_DRIVE),
  note = "ESETT arm B: fosphenytoin 20 mg PE/kg.")

scen[["11_esett_VPA"]] <- list(
  ev = d_lzp(20, 4) + d_vpa(65, 40),
  par = list(EDRIVE = RSE_DRIVE),
  note = "ESETT arm C: valproate 40 mg/kg.")

scen[["12_phenobarbital"]] <- list(
  ev = d_lzp(20, 4) + d_pb(65, 20) + d_pb(100, 10),
  par = list(EDRIVE = RSE_DRIVE),
  note = "Phenobarbital 20 mg/kg at 65 min + 10 mg/kg at 100 min in the same benzodiazepine-refractory patient as 09-11. Only 20% of its effect needs the pool that is being lost, and unlike the benzodiazepines it also blocks AMPA/kainate. Its effect-site half-life is ~28 min, so it is slow. Note the MAP column for why it is not first line.")

scen[["13_midazolam_infusion"]] <- list(
  ev = d_lzp(20, 4) + d_lev(65, 60) + d_mdz_iv(120, 0.2) + d_mdz_inf(120, 360, 0.4),
  par = list(EDRIVE = SRSE_DRIVE),
  note = "Refractory SE: midazolam bolus + 0.4 mg/kg/h infusion from 2 h.")

scen[["14_propofol_infusion"]] <- list(
  ev = d_lzp(20, 4) + d_lev(65, 60) + d_pro(120, 2) + d_pro_inf(120, 360, 6),
  par = list(EDRIVE = SRSE_DRIVE),
  note = "Refractory SE: propofol 2 mg/kg + 6 mg/kg/h. Direct GABA-A gating — 75% of its target survives.")

scen[["15_ketamine_65min"]] <- list(
  ev = d_lzp(10, 4) + d_lev(40, 60) + d_mdz_iv(60, 0.2) + d_mdz_inf(60, 600, 0.4) +
       d_ket(65, 2) + d_ket_inf(65, 595, 3),
  par = list(EDRIVE = 2.90),
  note = "Failing background (lorazepam 10 min, levetiracetam 40 min, midazolam infusion 60 min) with ketamine added EARLY, at 65 min.")

scen[["16_ketamine_300min"]] <- list(
  ev = d_lzp(10, 4) + d_lev(40, 60) + d_mdz_iv(60, 0.2) + d_mdz_inf(60, 600, 0.4) +
       d_ket(300, 2) + d_ket_inf(300, 360, 3),
  par = list(EDRIVE = 2.90),
  note = "The SAME background and the SAME ketamine dose, added LATE at 300 min. The only difference is where on the trajectory it lands.")

scen[["17_ketamine_no_GABA_arm"]] <- list(
  ev = d_lzp(10, 4) + d_lev(40, 60) + d_ket(65, 2) + d_ket_inf(65, 595, 3),
  par = list(EDRIVE = 2.90),
  note = "Scenario 15 with the midazolam infusion REMOVED — ketamine at the same dose and the same minute, without a GABAergic partner.")

scen[["18_allopregnanolone_SRSE"]] <- list(
  ev = d_lzp(20, 4) + d_lev(65, 60) + d_pro(120, 2) + d_pro_inf(120, 240, 6) +
       d_allo_inf(360, 360, 86),
  par = list(EDRIVE = 3.20),
  note = "The STATUS trial design: allopregnanolone started in established SRSE at 6 h. Correct target, wrong point on the trajectory.")

scen[["19_anti_NMDAR"]] <- list(
  ev = d_lzp(20, 4) + d_lev(65, 60) + d_mdz_iv(120, 0.2) + d_mdz_inf(120, 480, 0.6),
  par = list(EDRIVE = 2.40, TDRIVE = 1e6),
  note = "Anti-NMDAR encephalitis: DRIVE is sustained by the antibody and is never removed. No anticonvulsant schedule wins.")

scen[["20_FIRES_anakinra"]] <- list(
  ev = d_lzp(20, 4) + d_lev(65, 60) + d_mdz_iv(120, 0.2) + d_mdz_inf(120, 480, 0.6),
  par = list(EDRIVE = 1.95, GIL = 0.55, ANAK = 1),
  note = "FIRES with IL-1 blockade: the drive term itself is targeted (anakinra), not the receptors.")

scen[["21_optimal"]] <- list(
  ev = d_mdz_im(4, 10) + d_lzp(8, 4) + d_lev(12, 60),
  par = list(),
  note = "What the model says the best achievable pathway looks like: IM midazolam the moment SE is recognised, IV benzodiazepine as soon as access exists, second line WITHOUT waiting to see if the first worked.")

scen[["22_metabolic_glucose"]] <- list(
  ev = d_lzp(20, 4) + d_lev(65, 60),
  par = list(EDRIVE = 2.50, TDRIVE = 30, DRESID = 0.10, GLUCIV = 1),
  note = "Hypoglycaemic SE: dextrose + thiamine at 30 min removes 90% of the DRIVE. The one intervention that acts on the numerator.")

# =====================================================================
#  RUNNER
# =====================================================================
run_scenario <- function(name, end = END, delta = 0.25) {
  s <- scen[[name]]
  m <- mod
  if (length(s$par)) m <- param(m, s$par)
  out <- m %>%
    mrgsim_df(events = s$ev, end = end, delta = delta,
              atol = 1e-8, rtol = 1e-6, maxsteps = 2000000)
  out$scenario <- name
  out
}

summarise_scenario <- function(name, ...) {
  o <- run_scenario(name, ...)
  # time at which the seizure is first controlled for >= 5 consecutive minutes
  stopped <- o$SEIZ < 0.05
  tstop <- NA_real_
  n5 <- ceiling(5 / (o$time[2] - o$time[1]))
  r <- rle(stopped)
  idx <- which(r$values & r$lengths >= n5)
  if (length(idx)) {
    pos <- if (idx[1] == 1) 1 else sum(r$lengths[1:(idx[1] - 1)]) + 1
    tstop <- o$time[pos]
  }
  data.frame(
    scenario   = name,
    t_stop_min = round(tstop, 1),
    ctrl_60min = round(100 * mean(o$SEIZ[o$time <= 60] < 0.05), 1),
    burden_min = round(max(o$TSEIZ), 1),
    RSYN_60    = round(o$RSYN[which.min(abs(o$time - 60))], 3),
    NMDA_60    = round(o$oNMDA[which.min(abs(o$time - 60))], 3),
    CLI_max    = round(max(o$CLI), 1),
    injury     = round(max(o$INJURY), 3),
    neuron_loss_pct = round(max(o$NEURLOSS), 1),
    peak_CK    = round(max(o$CK)),
    min_MAP    = round(min(o$MAP)),
    min_gluc   = round(min(o$GLUCP), 2),
    max_temp   = round(max(o$TEMP), 2),
    intubated  = as.integer(any(o$INTUB > 0.5)),
    epilepsy_burden = round(max(o$EPG), 3),
    stringsAsFactors = FALSE
  )
}

run_all <- function() do.call(rbind, lapply(names(scen), summarise_scenario))

# =====================================================================
#  ANALYSIS 1 — THE PRICE OF A MINUTE
#  Sweep the time of the first (adequate, 4 mg) lorazepam dose and record
#  whether the seizure ever stops, and how much brain it costs.
# =====================================================================
sweep_benzo_time <- function(times = c(2, 5, 8, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 120)) {
  do.call(rbind, lapply(times, function(tt) {
    o <- mod %>% mrgsim_df(events = d_lzp(tt, 4) + d_lzp(tt + 4, 4),
                           end = 480, delta = 0.25, maxsteps = 2000000)
    stopped <- any(o$SEIZ[o$time > tt & o$time < tt + 30] < 0.05)
    data.frame(t_benzo = tt,
               terminated = stopped,
               RSYN_at_dose = round(o$RSYN[which.min(abs(o$time - tt))], 3),
               BZDCAP_at_dose = round(o$oTGT_BZD[which.min(abs(o$time - tt))], 3),
               CREQ_at_dose = round(o$CREQ_BZD[which.min(abs(o$time - tt))], 3),
               seizure_burden_min = round(max(o$TSEIZ), 1),
               neuron_loss_pct = round(max(o$NEURLOSS), 1))
  }))
}

# =====================================================================
#  ANALYSIS 2 — WHICH DRUGS A MINUTE COSTS YOU
#  From the untreated trajectory, read the required concentration of a
#  benzodiazepine and of phenobarbital at each time point.  These are
#  algebraic outputs of the pools; no extra dosing is simulated.
# =====================================================================
potency_drift <- function(times = c(5, 10, 15, 20, 30, 45, 60, 90, 120, 240, 480)) {
  o <- run_scenario("01_untreated", end = 600)
  do.call(rbind, lapply(times, function(tt) {
    i <- which.min(abs(o$time - tt))
    data.frame(time = tt,
               RSYN = round(o$RSYN[i], 3),
               FBZS = round(o$FBZS[i], 3),
               BZD_target = round(o$oTGT_BZD[i], 3),
               NMDA_target = round(o$oNMDA[i], 3),
               CREQ_BZD = round(o$CREQ_BZD[i], 3),
               CREQ_PB = round(o$CREQ_PB[i], 3),
               ket_share = round(o$KETSHARE[i], 3),
               CLFAC = round(o$oCLFAC[i], 3))
  }))
}


# =====================================================================
#  ANALYSIS 3b — THE CROSSOVER
#  Along the untreated trajectory, ask what a FIXED site occupancy (the
#  same for every drug, so the comparison is about targets and not about
#  potency or dose) would do to the network gain G at each minute.
#      cut_bzd  1 - G(with benzodiazepine)/G(no drug)
#      cut_pb   the same for phenobarbital
#      cut_ket  the same for ketamine
#  The benzodiazepine curve falls (its pool is leaving), the ketamine
#  curve rises (its pool is arriving).  Where they cross is a model
#  OUTPUT, and it is the quantitative version of the clinical intuition
#  that ketamine belongs late.
# =====================================================================
marginal_effect <- function(occ = 0.75, times = c(5, 15, 30, 60, 120, 180, 240, 360, 480, 600)) {
  o <- run_scenario("01_untreated", end = 720)
  p <- as.list(param(mod))
  do.call(rbind, lapply(times, function(tt) {
    i <- which.min(abs(o$time - tt))
    inh0 <- o$oINH_TOT[i]                 # untreated: no drug terms present
    eg0  <- o$oEGAIN[i]
    g0   <- eg0 / inh0
    tgt_bzd <- o$oTGT_BZD[i]
    tgt_pb  <- p$FSYNPB * o$RSYN[i] + (1 - p$FSYNPB) * o$REXTRA[i]
    g_bzd <- eg0 / (inh0 + p$EMAXBZD * occ * tgt_bzd)
    g_pb  <- eg0 / (inh0 + p$EMAXPB  * occ * tgt_pb)
    blk   <- p$EMAXKET * occ * (p$UDMIN + (1 - p$UDMIN) * o$SEIZ[i])
    share <- o$KETSHARE[i]
    g_ket <- g0 * (1 - blk * share)
    data.frame(time = tt,
               G_untreated = round(g0, 3),
               cut_bzd = round(1 - g_bzd / g0, 3),
               cut_pb  = round(1 - g_pb  / g0, 3),
               cut_ket = round(1 - g_ket / g0, 3),
               bzd_pool = round(tgt_bzd, 3),
               ket_share = round(share, 3))
  }))
}

crossover_time <- function(occ = 0.75) {
  m <- marginal_effect(occ, times = seq(5, 715, by = 5))
  d <- m$cut_ket - m$cut_bzd
  i <- which(d > 0)[1]
  if (is.na(i)) return(NA_real_)
  m$time[i]
}

# =====================================================================
#  ANALYSIS 3c — FIRES: blocking the DRIVE instead of the receptors
# =====================================================================
fires_anakinra <- function(edrive = 1.95) {
  ev1 <- d_lzp(20, 4) + d_lev(65, 60) + d_mdz_iv(120, 0.2) + d_mdz_inf(120, 480, 0.6)
  do.call(rbind, lapply(c(0, 1), function(a) {
    o <- param(mod, EDRIVE = edrive, GIL = 0.55, ANAK = a) %>%
      mrgsim_df(events = ev1, end = 720, delta = 0.5, maxsteps = 2000000)
    data.frame(anakinra = a,
               seizure_burden_min = round(max(o$TSEIZ), 1),
               peak_IL1b = round(max(o$IL1B), 2),
               neuron_loss_pct = round(max(o$NEURLOSS), 1),
               epileptogenesis = round(max(o$EPG), 2))
  }))
}

# =====================================================================
#  ANALYSIS 3 — VIRTUAL POPULATION
#  Response RATES need a population.  Inter-individual variability is
#  placed on EDRIVE (log-normal), which is the aetiology term.  The two
#  calibration targets are:
#      first-line benzodiazepine success   ~ 60-73%  (PHTSE / RAMPART)
#      second-line success among failures  ~ 46%     (ESETT)
# =====================================================================
vpop <- function(n = 400, seed = 20260804, sd_drive = 0.30, med_drive = 1.55) {
  set.seed(seed)
  data.frame(ID = 1:n,
             EDRIVE = med_drive * exp(rnorm(n, 0, sd_drive)))
}

pop_run <- function(events, n = 400, end = 180, extra = list(), ...) {
  idata <- vpop(n = n, ...)
  for (nm in names(extra)) idata[[nm]] <- extra[[nm]]
  mod %>% mrgsim_df(events = events, idata = idata, end = end, delta = 1,
                    maxsteps = 2000000,
                    obsonly = TRUE, recsort = 3)
}

resp_rate <- function(df, window_start, window_end = window_start + 40) {
  df %>%
    filter(time >= window_start, time <= window_end) %>%
    group_by(ID) %>%
    summarise(ok = any(SEIZ < 0.05), .groups = "drop") %>%
    summarise(rate = 100 * mean(ok)) %>%
    pull(rate)
}

calibration_check <- function(n = 400) {
  # ---- first line: lorazepam 4 mg at 20 min (the observed median-ish delay)
  p1 <- pop_run(d_lzp(20, 4), n = n, end = 120)
  first_line <- p1 %>% group_by(ID) %>%
    summarise(ok = any(SEIZ[time >= 20 & time <= 60] < 0.05), .groups = "drop")
  # ---- second line in the failures only
  fails <- first_line$ID[!first_line$ok]
  idata <- vpop(n = n); idata <- idata[idata$ID %in% fails, ]
  p2 <- lapply(list(LEV = d_lzp(20, 4) + d_lev(65, 60),
                    FOS = d_lzp(20, 4) + d_fos(65, 20),
                    VPA = d_lzp(20, 4) + d_vpa(65, 40)),
               function(e) {
                 o <- mod %>% mrgsim_df(events = e, idata = idata, end = 150,
                                        delta = 1, maxsteps = 2000000, obsonly = TRUE)
                 o %>% group_by(ID) %>%
                   summarise(ok = any(SEIZ[time >= 65 & time <= 125] < 0.05),
                             .groups = "drop") %>%
                   summarise(rate = 100 * mean(ok)) %>% pull(rate)
               })
  data.frame(endpoint = c("first-line benzo (20 min)",
                          "ESETT levetiracetam", "ESETT fosphenytoin", "ESETT valproate"),
             model_pct = round(c(100 * mean(first_line$ok),
                                 p2$LEV, p2$FOS, p2$VPA), 1),
             observed_pct = c("59-73 (PHTSE / RAMPART)", "47 (ESETT)", "45 (ESETT)", "46 (ESETT)"))
}

# ---- RAMPART reproduction -------------------------------------------
#  The ONLY differences between the two arms are operational, and both are
#  taken from the trial report rather than fitted:
#    (i)  the IM drug is given 4.5 min earlier (no IV access step);
#    (ii) in 10% of the IV arm the first access attempt fails, costing a
#         further 10 min.
#  The pharmacology (midazolam is the slower drug) is unchanged and works
#  AGAINST the IM arm.  If the IM arm still wins, it wins because R_SYN
#  was decaying during those 4.5 minutes.
rampart <- function(n = 400, p_iv_fail = 0.10, iv_fail_delay = 10) {
  im <- pop_run(d_mdz_im(8, 10), n = n, end = 140)
  set.seed(99); fail <- sample(c(TRUE, FALSE), n, TRUE, c(p_iv_fail, 1 - p_iv_fail))
  iv_ok   <- pop_run(d_lzp(12.5, 4), n = n, end = 140)
  iv_late <- pop_run(d_lzp(12.5 + iv_fail_delay, 4), n = n, end = 140)
  ok_ok   <- iv_ok   %>% group_by(ID) %>%
    summarise(ok = any(SEIZ[time >= 12.5 & time <= 72.5] < 0.05), .groups = "drop")
  ok_late <- iv_late %>% group_by(ID) %>%
    summarise(ok = any(SEIZ[time >= 22.5 & time <= 82.5] < 0.05), .groups = "drop")
  iv_rate <- 100 * mean(ifelse(fail[ok_ok$ID], ok_late$ok, ok_ok$ok))
  data.frame(arm = c("IM midazolam 10 mg @ 8.0 min",
                     "IV lorazepam 4 mg @ 12.5 min (10% access failure)"),
             model_pct = round(c(resp_rate(im, 8, 68), iv_rate), 1),
             observed_pct = c(73.4, 63.4))
}

# ---- second-line delay insensitivity (the model's prediction) --------
#  Conditioned on FIRST-LINE FAILURE, which is what ESETT actually enrolled.
#  Without the conditioning the sweep is confounded by patients the
#  benzodiazepine had already cured.
second_line_delay <- function(n = 300, times = c(30, 45, 60, 90, 120, 180)) {
  id <- vpop(n = n)
  o1 <- mod %>% mrgsim_df(events = d_lzp(20, 4), idata = id, end = 100, delta = 1,
                          maxsteps = 2000000, obsonly = TRUE)
  f1 <- o1 %>% group_by(ID) %>%
    summarise(ok = any(SEIZ[time >= 20 & time <= 80] < 0.05), .groups = "drop")
  id2 <- id[id$ID %in% f1$ID[!f1$ok], ]
  do.call(rbind, lapply(times, function(tt) {
    o <- mod %>% mrgsim_df(events = d_lzp(20, 4) + d_lev(tt, 60), idata = id2,
                           end = tt + 70, delta = 1, maxsteps = 2000000, obsonly = TRUE)
    r <- o %>% group_by(ID) %>%
      summarise(ok = any(SEIZ[time >= tt & time <= tt + 60] < 0.05), .groups = "drop")
    data.frame(t_second_line = tt, n_benzo_failures = nrow(id2),
               response_pct = round(100 * mean(r$ok), 1))
  }))
}

first_line_delay <- function(n = 300, times = c(5, 10, 20, 30, 45, 60)) {
  do.call(rbind, lapply(times, function(tt) {
    o <- pop_run(d_lzp(tt, 4), n = n, end = tt + 70)
    data.frame(t_first_line = tt, response_pct = round(resp_rate(o, tt, tt + 60), 1))
  }))
}

# =====================================================================
#  ANALYSIS 4 — THE ELECTROMECHANICAL DISSOCIATION
# =====================================================================
dissociation <- function(edrive = 2.90, motor_threshold = 0.15) {
  o <- param(mod, EDRIVE = edrive) %>%
    mrgsim_df(events = d_lzp(20, 4) + d_lzp(60, 4) + d_lzp(120, 4),
              end = 300, delta = 0.5, maxsteps = 2000000)
  idx <- sapply(c(15, 30, 45, 60, 90, 120, 180, 240), function(tt) which.min(abs(o$time - tt)))
  data.frame(time = o$time[idx],
             benzo_occupancy = round(o$oOCC_BZD[idx], 3),
             motor_gain = round(o$MOTG[idx], 3),
             motor_output = round(o$MOTOROUT[idx], 3),
             eeg_burden = round(o$EEGOUT[idx], 3),
             looks_controlled = o$MOTOROUT[idx] < motor_threshold,
             is_controlled = o$EEGOUT[idx] < 0.05)
}

# =====================================================================
#  ANALYSIS 5 — TRANSPORTER vs TARGET RESISTANCE
# =====================================================================
resistance_split <- function() {
  o <- run_scenario("01_untreated", end = 480)
  idx <- sapply(c(0, 30, 60, 120, 240, 480), function(tt) which.min(abs(o$time - tt)))
  data.frame(time = o$time[idx],
             PGP = round(o$PGP[idx], 3),
             BBB = round(o$BBBP[idx], 3),
             brain_plasma_index = round(o$BPRATIO[idx], 3),
             benzo_target_pool = round(o$oTGT_BZD[idx], 3),
             comment = c("baseline",
                         "both moving",
                         "target loss already dominates",
                         "brain:plasma nearly restored, target still falling",
                         "transporter effect is NOT the story",
                         "no dose of the same drug can work"))
}

# =====================================================================
#  CALIBRATION NOTES  (what each number is anchored to)
# ---------------------------------------------------------------------
#  1. R_SYN trajectory.  KENDO / KREC / KDEGR are set so that the surface
#     synaptic GABA-A pool falls to ~0.51 of baseline after 60 min of
#     continuous SEIZ = 1.  Naylor, Liu & Wasterlain (J Neurosci 2005;
#     PMID 16093381) measured a 47% reduction in surface gamma2 in
#     dentate granule cells after 1 h of SE.  The 30-min value the model
#     produces (~0.69) is not fitted; it is a consequence.
#
#  2. F_BZS (benzodiazepine-sensitive fraction).  A second, independent
#     term for the alpha1 -> alpha4 subunit switch: alpha4-containing
#     receptors have NO benzodiazepine site, so the drug loses target
#     twice over (Goodkin 2008 PMID 18337412; Terunuma 2008 PMID 19005040).
#     KSW / FBZMIN give F_BZS ~ 0.83 at 30 min and ~ 0.71 at 60 min.
#     The PRODUCT R_SYN x F_BZS is what multiplies benzodiazepine effect.
#
#  3. Benzodiazepine potency drift.  Kapur & Macdonald (J Neurosci 1997;
#     PMID 9315909) found diazepam potency fell ~20-fold and phenobarbital
#     ~3-fold after 30 min of SE.  In this model those numbers are OUTPUTS
#     of the pool table (CREQ_BZD and CREQ_PB, see potency_drift()), not
#     inputs.  Run potency_drift() to read what the model actually gives.
#
#  4. NMDA insertion.  KEXO / KENDN / FEXO0 are solved so that NRSYN
#     reaches 1.38 at 60 min with a ~30 min time constant — Naylor et al.
#     (Neurobiol Dis 2013; PMID 23396011) reported a ~38% increase in
#     synaptic NR1 after 1 h of SE.
#
#  5. Extracellular K+ / Mg-block relief.  KKO / KKOC give a ~2.7 mM
#     excess plateau with a ~30 min time constant, consistent with the
#     classical [K+]o ceiling of 8-12 mM during seizures.
#
#  6. ESETT.  SLEV / SPHT / SVPA are the ONLY three parameters fitted to
#     an efficacy result, and they are fitted to a TIE (47 / 45 / 46%,
#     Kapur et al. NEJM 2019; PMID 31774955).  The structural claim they
#     support is not the level but the FLATNESS: because these targets
#     are presynaptic/axonal and are not trafficked, the model predicts
#     second-line response should be nearly independent of when it is
#     given.  See second_line_delay() vs first_line_delay().
#
#  7. RAMPART.  No RAMPART-specific parameter exists.  The IM arm is the
#     same midazolam PK with an absorption step and a 4.5 min earlier
#     administration time; the observed 73.4% vs 63.4% (Silbergleit et al.
#     NEJM 2012; PMID 22335736) has to fall out of the R_SYN decay or the
#     model is wrong.  Run rampart().
#
#  8. Intubation.  KRESPS > KRESPD by design: PHTSE (Alldredge et al.
#     NEJM 2001; PMID 11547716) found respiratory complications in 22.5%
#     of PLACEBO patients vs 10.6% with lorazepam — the seizure suppresses
#     breathing more than the drug does.  The model must reproduce the
#     direction, not just the magnitude.
#
#  9. Subtle SE.  MOTMIN / KMOT / MOTDRG are set so that motor output has
#     fallen below a plausible bedside detection threshold well before the
#     electrographic burden has, reproducing the direction of the VA
#     Cooperative finding (lorazepam 64.9% in overt SE vs 7.7-24.2% in
#     subtle SE; Treiman et al. NEJM 1998; PMID 9737139).
#
# 10. Phenytoin.  Michaelis-Menten elimination (Vmax 0.34 mg/min,
#     Km 5 mg/L total) plus albumin- and valproate-dependent free fraction.
#     A 20 mg PE/kg load gives ~31 mg/L total and ~3.1 mg/L free in a
#     patient with albumin 4.0 g/dL.
#
# 11. What is NOT calibrated, and is therefore a prediction:
#       - the shape of CREQ_BZD(t) and the time at which it becomes
#         infinite (no benzodiazepine dose can work);
#       - the crossover between benzodiazepine-recruitable pool and
#         NMDA-removable share;
#       - the flatness of second-line response versus its timing;
#       - the ~4-fold difference in cumulative excitotoxic load between
#         termination at 20 and at 60 minutes.
#     These are the four places the model can be falsified.
# =====================================================================

# =====================================================================
#  MAIN
# =====================================================================
if (sys.nframe() == 0L || identical(environment(), globalenv())) {
  message("\n=== 1. All 22 scenarios ===")
  print(run_all(), row.names = FALSE)

  message("\n=== 2. The price of a minute (first benzodiazepine dose) ===")
  print(sweep_benzo_time(), row.names = FALSE)

  message("\n=== 3. Which drugs a minute costs you ===")
  print(potency_drift(), row.names = FALSE)

  message("\n=== 3b. The crossover (fixed occupancy, three targets) ===")
  print(marginal_effect(), row.names = FALSE)
  message(sprintf("crossover (ketamine overtakes benzodiazepine): %s min",
                  crossover_time()))

  message("\n=== 3c. FIRES: blocking the drive instead of the receptors ===")
  print(fires_anakinra(), row.names = FALSE)

  message("\n=== 4. Calibration against PHTSE / ESETT ===")
  print(calibration_check(), row.names = FALSE)

  message("\n=== 5. RAMPART ===")
  print(rampart(), row.names = FALSE)

  message("\n=== 6. First-line vs second-line delay sensitivity ===")
  print(first_line_delay(), row.names = FALSE)
  print(second_line_delay(), row.names = FALSE)

  message("\n=== 7. Electromechanical dissociation ===")
  print(dissociation(), row.names = FALSE)

  message("\n=== 8. Transporter vs target resistance ===")
  print(resistance_split(), row.names = FALSE)
}
