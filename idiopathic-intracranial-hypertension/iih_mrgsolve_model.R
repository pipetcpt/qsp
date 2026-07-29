# =============================================================================
# Idiopathic Intracranial Hypertension (IIH / pseudotumour cerebri)
# 44-state Quantitative Systems Pharmacology model for mrgsolve
# =============================================================================
#
# WHAT THIS MODEL IS BUILT AROUND
# -------------------------------
# Intracranial pressure is not written here as a tank that fills.  It is the
# FIXED POINT of the Davson relation
#
#       ICP  =  P_sss  +  R_out * F_form                                   (1)
#
# closed through a mechanically collapsible transverse sinus
#
#       P_sss =  P_cv  +  Q_v * ( RTS0 + RTSF + RTSC(ICP) )                (2)
#
# The transverse sinus is a thin-walled tube sitting INSIDE the cranium, so its
# transmural pressure is (sinus pressure - ICP).  When ICP exceeds the local
# sinus pressure by more than the wall can support (PCRIT) the lumen narrows,
# RTSC rises, P_sss rises, and ICP rises further.  (1) + (2) is a loop with a
# gain
#
#       g  =  dP_sss / dICP                                                (3)
#
# and essentially every clinically puzzling feature of this disease is a
# property of (3) rather than of any drug:
#
#   * an intervention that acts on R_out or F_form is amplified by 1/(1-g);
#   * an intervention that acts on P_cv is NOT amplified (net coefficient ~1),
#     because central venous pressure raises the sinus pressure and the
#     intracranial pressure equally and therefore cancels out of the collapse
#     drive;
#   * as g -> 1 the fixed point is annihilated in a saddle-node: there is no
#     stable low-pressure state left, which is what "fulminant IIH" is;
#   * for a band of parameters BELOW that fold the system is bistable, so a
#     drug can move the STATE rather than a parameter and remission can outlast
#     withdrawal;
#   * cerebral blood flow MULTIPLIES the stenosis resistance in (2), so a drug
#     that vasodilates (acetazolamide, through the metabolic acidosis it
#     causes) pays a penalty proportional to the very stenosis that amplifies
#     its benefit.
#
# Three further structural commitments are worth stating before the code:
#
#   * R_out is NOT an independent lesion.  The arachnoid granulations and the
#     meningeal/cervical lymphatics discharge INTO the venous compartment, so
#     R_out carries a venous-pressure-dependent term (KROUTP).  This is the
#     second, amplified route by which body habitus acts, and it is what lets
#     weight loss out-perform a pure hydrostatic account of itself.
#   * the optic nerve reads PONS - IOP, not ICP.  Acetazolamide inhibits
#     carbonic anhydrase in the ciliary epithelium as well as the choroid
#     plexus, so it lowers IOP while lowering ICP and erodes part of its own
#     benefit AT THE DISC while none of it at the headache endpoint.
#   * the OCT RNFL number is a SUM of two states with opposite meanings
#     (surviving axons + prelaminar swelling), and swelling requires living
#     axons to swell.  A recovering disc and a dying disc therefore pass
#     through the same measured thickness.
#
# UNITS
# -----
#   time      day
#   pressure  mmHg          (clinical opening pressures are in cmH2O:
#                            1 mmHg = 1.35951 cmH2O; both are tabulated)
#   flow      mL/min        (converted inside $ODE by *1440)
#   RNFL      micrometre
#   weight    kg
#
# VERIFICATION
# ------------
# Every number quoted in README.md was produced by iih_reference_check.py, an
# independent dependency-free Python/RK4 transcription of the equations below
# (same parameter names, same right-hand sides).  Its full output is committed
# as iih_reference_output.txt, and its ANALYSIS 0 prints a state/derivative
# fingerprint at a deliberately awkward point (all seven drugs on board, stent
# 0.55, shunt 0.70, ONSF 0.80, lumbar puncture running) so that this file can
# be checked term by term rather than trusted.
#
# =============================================================================

library(mrgsolve)
library(dplyr)
library(ggplot2)

iih_code <- '
$PROB
# Idiopathic intracranial hypertension: 44-state QSP model
# ICP as the fixed point of the Davson relation closed through a
# collapsible transverse sinus.

$PARAM @annotated
// ---- cerebral venous outflow / transverse sinus -------------------------
QV0     : 600.0   : cerebral venous outflow through the transverse sinuses (mL/min)
RTS0    : 0.00333 : baseline sinus resistance (mmHg.min/mL; 2 mmHg gradient)
RTSCMAX : 0.0200  : maximum reversible compressive resistance (mmHg.min/mL)
PCRIT   : 10.95   : transmural pressure the sinus wall supports (mmHg) [SUSCEPTIBILITY]
PSTIFF  : 4.14    : width of the collapse transition (mmHg) [SUSCEPTIBILITY]
FSEG    : 0.50    : fraction of the gradient upstream of the collapsible segment
TAURTS  : 0.0097  : collapse time constant (day; 14 min)
KREM    : 0.0016  : compression -> fixed remodelled stenosis (1/day)
KREV    : 0.0025  : reversal of remodelling (1/day)
RTSFMAX : 0.0100  : ceiling of the fixed stenosis component (mmHg.min/mL)
KRESTEN : 0.00035 : loss of stent patency, restenosis (1/day)

// ---- CSF hydrodynamics ---------------------------------------------------
FFORM0  : 0.350   : choroid plexus CSF secretion (mL/min; 504 mL/day)
ROUTB   : 0.00    : pressure-independent CSF outflow resistance (mmHg.min/mL)
KROUTP  : 1.69    : outflow resistance per mmHg of central venous pressure
TAUROUT : 1.00    : outflow-resistance equilibration (day)
E1      : 0.200   : elastance coefficient, C = 1/(E1*ICP) (1/mL)
RSHUNT  : 2.50    : shunt resistance (mmHg.min/mL)
PSHUNT  : 7.40    : shunt opening pressure (mmHg; 10 cmH2O)
KSHFAIL : 0.0011  : shunt obstruction hazard (1/day)

// ---- body habitus --------------------------------------------------------
HT      : 1.65    : height (m)
TAUWT   : 60.0    : weight trajectory time constant, diet (day)
TAUWTB  : 25.0    : weight trajectory time constant, bariatric surgery (day)
KIAP    : 0.350   : intra-abdominal pressure per BMI unit above BMIREF (mmHg)
BMIREF  : 22.0    : BMI at which IAP is taken as zero
KPCV    : 1.00    : transmission of IAP to central venous pressure
PCV0    : 4.50    : central venous pressure at BMIREF (mmHg)

// ---- choroid plexus pharmacology ----------------------------------------
TAUSEC  : 0.150   : secretory machinery adaptation (day; 3.6 h)
EMAXCA  : 0.420   : maximum fractional fall in F_form from CA inhibition
EC50ACZ : 12.0    : acetazolamide plasma EC50 (mg/L)
EC50TPM : 9.00    : topiramate plasma EC50 (mg/L)
POTTPM  : 0.45    : relative carbonic-anhydrase potency of topiramate
EMAXGLP : 0.300   : maximum fractional fall in F_form from GLP-1R agonism
EC50EX  : 0.055   : exenatide EC50 (ng/mL)
EC50SEM : 12.0    : semaglutide EC50 (ng/mL)
EMAXFUR : 0.200   : maximum fractional fall in F_form from NKCC1 inhibition
EC50FUR : 0.25    : furosemide EC50 (mg/L)
EMAXGC  : 0.100   : glucocorticoid drive on secretion (11b-HSD1 route)
EC50GC  : 1.00    : half-maximal local cortisol tone
KANDR   : 0.0120  : adiposity/androgen drive on secretion per BMI unit above 25

// ---- acid-base and cerebral blood flow ----------------------------------
HCO30   : 25.0    : baseline plasma bicarbonate (mEq/L)
FHCO3   : 0.280   : maximum fractional fall in bicarbonate
TAUHCO3 : 3.00    : bicarbonate time constant (day)
KACID   : 0.550   : CBF gain per fractional fall in bicarbonate
TAUQV   : 0.050   : CBF time constant (day)

// ---- 11beta-HSD1 --------------------------------------------------------
TAUCORT : 0.500   : local cortisol tone time constant (day)
KCORTBMI: 0.035   : cortisol tone per BMI unit above 25
EMAXAZD : 0.850   : AZD4017 target engagement
EC50AZD : 55.0    : AZD4017 EC50 (ng/mL)

// ---- optic nerve head ---------------------------------------------------
TAUONS  : 0.020   : sheath equilibration with ICP (day)
KFEN    : 22.0    : fenestration conductance after ONSF (1/day)
PORB    : 6.00    : orbital tissue pressure (mmHg)
IOP0    : 15.5    : intraocular pressure at BMI 25 (mmHg)
KIOPBMI : 0.100   : IOP per BMI unit above 25 (mmHg)
EMAXIOP : 2.50    : maximum IOP fall from carbonic-anhydrase inhibition (mmHg)
TAUIOP  : 0.100   : IOP time constant (day)
TLPDTHR : 3.00    : translaminar gradient at which axoplasmic stasis begins (mmHg)
TAUSTAS : 2.00    : stasis time constant (day)
SSTAS   : 2.20    : logistic width of the stasis switch (mmHg)
KSW     : 1.20    : swelling rate per mmHg supra-threshold gradient (um/day)
RNFLSMAX: 320.0   : ceiling of prelaminar swelling (um)
KRES    : 0.100   : resorption of swelling (1/day)
AXON0   : 98.0    : baseline RNFL axonal thickness (um)
AXONFL  : 38.0    : residual RNFL floor (um)
KAXL    : 0.00120 : stasis-driven axon loss (1/day)
KISCH   : 0.00600 : ischaemic axon loss above TLPDISCH (1/day)
TLPDISCH: 22.0    : translaminar gradient at which ONH perfusion fails (mmHg)
SMD     : 0.300   : mean deviation per um of lost axon (dB/um)
SREV    : 2.00    : maximum reversible field loss (dB)
TAUMD   : 10.0    : mean-deviation time constant (day)

// ---- headache -----------------------------------------------------------
ICPHTHR : 13.0    : ICP above which pressure drives headache (mmHg)
W1HA    : 0.550   : weight of the ICP term in the headache drive
W2HA    : 0.300   : weight of the venous-gradient term
W3HA    : 0.250   : weight of the CGRP term
KSENS   : 0.120   : central sensitisation rate (1/day)
KDESENS : 0.045   : desensitisation rate (1/day)
FMOHDES : 0.600   : fractional block of desensitisation by medication overuse
KMOH    : 0.030   : medication-overuse accrual (1/day)
KMOHOFF : 0.020   : medication-overuse resolution (1/day)
TAUMHD  : 3.00    : monthly-headache-day time constant (day)
TAUCGRP : 1.00    : CGRP tone time constant (day)
KCG     : 0.700   : CGRP gain from sensitisation
KCGICP  : 0.030   : CGRP gain from pressure

// ---- other symptoms -----------------------------------------------------
TINNREF : 12.0    : venous gradient giving mid-scale pulsatile tinnitus (mmHg)
TAUTINN : 0.200   : tinnitus time constant (day)
KDIPL   : 0.020   : abducens palsy accrual (1/day/mmHg)
KDIPLR  : 0.060   : abducens palsy recovery (1/day)
ICPDIPL : 26.0    : ICP threshold for abducens stretch (mmHg)

// ---- tolerability / adherence -------------------------------------------
EC50PAR : 6.00    : acetazolamide concentration for half-maximal paraesthesia (mg/L)
EC50COG : 5.50    : topiramate concentration for half-maximal cognitive slowing (mg/L)
TAUAE   : 0.500   : adverse-effect time constant (day)
TAUADH  : 7.00    : adherence time constant (day)
ADHMAX  : 0.950   : adherence ceiling
KADHPAR : 0.450   : adherence loss per unit paraesthesia
KADHCOG : 0.500   : adherence loss per unit cognitive slowing
TAUSTONE: 10.0    : stone-risk time constant (day)

// ---- PK ------------------------------------------------------------------
KAACZ   : 28.8    : acetazolamide absorption (1/day)
VCACZ   : 14.0    : acetazolamide central volume (L)
CLACZ   : 46.0    : acetazolamide clearance (L/day)
KONACZ  : 2.60    : red-cell carbonic-anhydrase association (L/mg/day)
KOFFACZ : 0.90    : red-cell dissociation (1/day)
RBCMAX  : 130.0   : red-cell binding capacity (mg)
KATPM   : 24.0    : topiramate absorption (1/day)
VCTPM   : 55.0    : topiramate volume (L)
CLTPM   : 43.0    : topiramate clearance (L/day)
FTPM    : 0.80    : topiramate bioavailability
KAFUR   : 36.0    : furosemide absorption (1/day)
VCFUR   : 12.0    : furosemide volume (L)
CLFUR   : 130.0   : furosemide clearance (L/day)
FFUR    : 0.55    : furosemide bioavailability
KAEX    : 8.30    : exenatide absorption (1/day)
VCEX    : 28.0    : exenatide volume (L)
CLEX    : 190.0   : exenatide clearance (L/day)
KASEM   : 0.55    : semaglutide absorption (1/day)
VCSEM   : 12.5    : semaglutide volume (L)
CLSEM   : 0.87    : semaglutide clearance (L/day)
KAPRED  : 24.0    : prednisolone absorption (1/day)
VCPRED  : 35.0    : prednisolone volume (L)
CLPRED  : 190.0   : prednisolone clearance (L/day)
KAAZD   : 20.0    : AZD4017 absorption (1/day)
VCAZD   : 60.0    : AZD4017 volume (L)
CLAZD   : 165.0   : AZD4017 clearance (L/day)
FAZD    : 0.65    : AZD4017 bioavailability
FWTTPM  : 0.100   : maximum fractional weight loss on topiramate
ED50WTPM: 90.0    : topiramate dose for half-maximal weight loss (mg/day)

$PARAM @annotated
// ---- prescribed regimen (continuous daily input scaled by adherence) -----
// Dosing is written as a continuous daily rate rather than as discrete events
// so that the tolerability -> adherence -> exposure loop closes INSIDE the ODE
// system.  For acetazolamide the difference between a continuous input and
// 500 mg TID is < 4 % in the daily-average concentration that drives the
// pharmacology, because the secretory state FSEC has a 3.6 h time constant.
DOSACZ  : 0.0  : acetazolamide (mg/day)
DOSTPM  : 0.0  : topiramate (mg/day)
DOSFUR  : 0.0  : furosemide (mg/day)
DOSEX   : 0.0  : exenatide (ug/day; 10 ug BID = 20)
DOSSEM  : 0.0  : semaglutide (mg/week)
DOSPRED : 0.0  : prednisolone equivalent (mg/day)
DOSAZD  : 0.0  : AZD4017 (mg/day)
TSTART  : 0.0  : day therapy starts
TSTOP   : 1e9  : day therapy stops
FDIET   : 0.0  : fractional weight-loss target, diet
TDIET   : 0.0  : day the weight programme starts
FBAR    : 0.0  : fractional weight-loss target, bariatric surgery
TBAR    : 1e9  : day of bariatric surgery
TSTENT  : 1e9  : day of venous sinus stenting
TSHUNT  : 1e9  : day of shunt insertion
TONSF   : 1e9  : day of optic nerve sheath fenestration
TLP     : 1e9  : day of a 25 mL lumbar puncture
ANALG   : 0.25 : analgesic exposure (0-1)
CITRATE : 0.0  : potassium citrate co-prescription (0/1)
WT0     : 103.4: weight at presentation (kg)

$CMT @annotated
// ---- PK (13) ------------------------------------------------------------
AACZG  : acetazolamide, gut (mg)
AACZP  : acetazolamide, plasma (mg)
AACZR  : acetazolamide, red-cell carbonic anhydrase (mg)
ATPMG  : topiramate, gut (mg)
ATPMC  : topiramate, plasma (mg)
AFURC  : furosemide, plasma (mg)
AEXSC  : exenatide, subcutaneous depot (ug)
AEXC   : exenatide, plasma (ug)
ASEMS  : semaglutide, subcutaneous depot (ug)
ASEMC  : semaglutide, plasma (ug)
APRED  : prednisolone, plasma (mg)
AAZDG  : AZD4017, gut (mg)
AAZDC  : AZD4017, plasma (mg)
// ---- CSF secretion and systemic consequences (5) ------------------------
FSEC   : relative choroid-plexus secretory capacity
HCO3   : plasma bicarbonate (mEq/L)
QVREL  : relative cerebral blood flow
ROUT   : CSF outflow resistance (mmHg.min/mL)
SHUNTP : shunt patency (0-1)
// ---- the venous loop (4) ------------------------------------------------
ICP    : intracranial pressure (mmHg)
RTSC   : reversible compressive sinus resistance (mmHg.min/mL)
RTSF   : remodelled fixed sinus stenosis (mmHg.min/mL)
STENT  : stent patency (0-1)
// ---- body habitus (4) ---------------------------------------------------
WT     : body weight (kg)
IAP    : intra-abdominal pressure (mmHg)
PCV    : central venous pressure (mmHg)
CORTL  : local cortisol tone at the choroid plexus
// ---- optic nerve and vision (7) -----------------------------------------
PONS   : optic nerve sheath CSF pressure (mmHg)
IOP    : intraocular pressure (mmHg)
STASIS : axoplasmic transport stasis (0-1)
RNFLS  : prelaminar swelling (um)
AXON   : surviving RNFL axonal thickness (um)
MD     : perimetric mean deviation (dB)
CUMEXP : cumulative supra-threshold translaminar exposure (mmHg.day)
// ---- headache and other symptoms (6) -----------------------------------
STG    : central trigeminal sensitisation (0-1)
CGRP   : CGRP tone (relative)
MOH    : medication-overuse state (0-1)
MHD    : monthly headache days (0-30)
TINN   : pulsatile tinnitus (0-10)
DIPL   : abducens palsy / diplopia burden (0-1)
// ---- tolerability (5) ---------------------------------------------------
PARES  : paraesthesia burden (0-1)
COGSL  : cognitive slowing (0-1)
ADH    : adherence (0-1)
KSTONE : nephrolithiasis risk index (0-1)
ONSFS  : optic nerve sheath fenestration state (0-1)

$GLOBAL
#define POSF(x) ((x) > 0.0 ? (x) : 0.0)
#define SIGF(x) (1.0 / (1.0 + exp(-((x) > 40.0 ? 40.0 : ((x) < -40.0 ? -40.0 : (x))))))
#define CMH2O 1.35951

$MAIN
// the untreated secretory drive is a function of habitus, so the initial
// condition is written from WT0 rather than hard-coded
double _bmi0  = WT0 / (HT * HT);
double _cort0 = 1.0 + KCORTBMI * POSF(_bmi0 - 25.0);
double _egc0  = EMAXGC * _cort0 / (EC50GC + _cort0);
double _fs0   = (1.0 + _egc0) * (1.0 + KANDR * POSF(_bmi0 - 25.0));
double _iap0  = KIAP * POSF(_bmi0 - BMIREF);
double _pcv0  = PCV0 + KPCV * _iap0;

FSEC_0   = _fs0;
HCO3_0   = HCO30;
QVREL_0  = 1.0;
ROUT_0   = ROUTB + KROUTP * _pcv0;
SHUNTP_0 = 0.0;
ICP_0    = 25.8;          // overwritten by the steady-state pre-run (see $CAPTURE note)
RTSC_0   = 0.5 * RTSCMAX;
RTSF_0   = 0.0;
STENT_0  = 0.0;
WT_0     = WT0;
IAP_0    = _iap0;
PCV_0    = _pcv0;
CORTL_0  = _cort0;
PONS_0   = 25.8;
IOP_0    = IOP0 + KIOPBMI * POSF(_bmi0 - 25.0);
STASIS_0 = 0.0;
RNFLS_0  = 0.0;
AXON_0   = AXON0;
MD_0     = 0.0;
CGRP_0   = 1.0;
MHD_0    = 1.5;
ADH_0    = ADHMAX;

$ODE
double ON = ((SOLVERTIME >= TSTART) && (SOLVERTIME < TSTOP)) ? 1.0 : 0.0;
double BMI = WT / (HT * HT);

// ---------------------------------------------------------------- PK ------
double CACZ = AACZP / VCACZ;
double RBCIN  = KONACZ * CACZ * (RBCMAX - AACZR);
double RBCOUT = KOFFACZ * AACZR;
dxdt_AACZG = ON * DOSACZ * ADH - KAACZ * AACZG;
dxdt_AACZP = KAACZ * AACZG - CLACZ / VCACZ * AACZP - RBCIN + RBCOUT;
dxdt_AACZR = RBCIN - RBCOUT;

double CTPM = ATPMC / VCTPM;
dxdt_ATPMG = ON * DOSTPM * FTPM * ADH - KATPM * ATPMG;
dxdt_ATPMC = KATPM * ATPMG - CLTPM / VCTPM * ATPMC;

double CFUR = AFURC / VCFUR;
dxdt_AFURC = ON * DOSFUR * FFUR * ADH - CLFUR / VCFUR * AFURC;

double CEX = AEXC / VCEX;                      // ng/mL
dxdt_AEXSC = ON * DOSEX - KAEX * AEXSC;
dxdt_AEXC  = KAEX * AEXSC - CLEX / VCEX * AEXC;

double CSEM = ASEMC / VCSEM;                   // ng/mL
dxdt_ASEMS = ON * DOSSEM / 7.0 * 1000.0 - KASEM * ASEMS;
dxdt_ASEMC = KASEM * ASEMS - CLSEM / VCSEM * ASEMC;

double CPRED = APRED / VCPRED;
dxdt_APRED = ON * DOSPRED - CLPRED / VCPRED * APRED;

double CAZD = AAZDC / VCAZD * 1000.0;          // ng/mL
dxdt_AAZDG = ON * DOSAZD * FAZD - KAAZD * AAZDG;
dxdt_AAZDC = KAAZD * AAZDG - CLAZD / VCAZD * AAZDC;

// ------------------------------------------- carbonic anhydrase engagement -
double XCA   = CACZ / EC50ACZ + POTTPM * CTPM / EC50TPM;
double INHCA = EMAXCA * XCA / (1.0 + XCA);     // fractional fall in F_form
double OCCCA = XCA / (1.0 + XCA);              // 0-1 target occupancy

double EGLP = EMAXGLP * (CEX / (EC50EX + CEX) + CSEM / (EC50SEM + CSEM));
if (EGLP > EMAXGLP) EGLP = EMAXGLP;
double EFUR = EMAXFUR * CFUR / (EC50FUR + CFUR);
double EGC  = EMAXGC * CORTL / (EC50GC + CORTL);
double EANDR = KANDR * POSF(BMI - 25.0);

double FSECT = (1.0 - INHCA) * (1.0 - EGLP) * (1.0 - EFUR)
             * (1.0 + EGC) * (1.0 + EANDR);
dxdt_FSEC = (FSECT - FSEC) / TAUSEC;
double FFORM = FFORM0 * FSEC;

// --------------------------------------- acidosis -> cerebral blood flow ---
dxdt_HCO3 = (HCO30 * (1.0 - FHCO3 * OCCCA) - HCO3) / TAUHCO3;
double QVRELT = 1.0 + KACID * (HCO30 - HCO3) / HCO30;
dxdt_QVREL = (QVRELT - QVREL) / TAUQV;
double QV = QV0 * QVREL;

// ------------------------------------------------------- 11beta-HSD1 -------
double EAZD  = EMAXAZD * CAZD / (EC50AZD + CAZD);
double CORTT = (1.0 + KCORTBMI * POSF(BMI - 25.0)) * (1.0 - EAZD)
             + 0.18 * CPRED;
dxdt_CORTL = (CORTT - CORTL) / TAUCORT;

// ================= the venous loop: equations (2) and (3) =================
double RTSCMX = RTSCMAX * (1.0 - STENT);
double RFIX   = RTS0 + RTSF * (1.0 - STENT);
double RTOT   = RFIX + RTSC;
double GRAD   = QV * RTOT;                        // SSS -> jugular gradient
double PSSS   = PCV + GRAD;
double PSEG   = PCV + (1.0 - FSEG) * GRAD;        // pressure at the collapse
double DRIVE  = (ICP - PSEG - PCRIT) / PSTIFF;
double RTSCT  = RTSCMX * SIGF(DRIVE);
dxdt_RTSC = (RTSCT - RTSC) / TAURTS;
dxdt_RTSF = KREM * (RTSC / RTSCMAX) * (RTSFMAX - RTSF) * (1.0 - STENT)
          - KREV * RTSF;
dxdt_STENT = -KRESTEN * STENT;

// ============================= equation (1) ==============================
double FABS = POSF(ICP - PSSS) / ROUT;
double FSH  = SHUNTP * POSF(ICP - PSHUNT) / RSHUNT;
double FLP  = ((SOLVERTIME >= TLP) && (SOLVERTIME < TLP + 25.0/2.5/1440.0))
              ? 2.5 : 0.0;
dxdt_ICP  = (FFORM - FABS - FSH - FLP) * 1440.0 * E1 * ICP;
dxdt_ROUT = (ROUTB + KROUTP * PCV - ROUT) / TAUROUT;
dxdt_SHUNTP = -KSHFAIL * SHUNTP;

// --------------------------------------------------------- body habitus ---
double FDT  = (SOLVERTIME >= TDIET) ? FDIET : 0.0;
double FBR  = (SOLVERTIME >= TBAR)  ? FBAR  : 0.0;
double FSTER = 0.0016 * DOSPRED * ON;
double FTPMW = (FWTTPM * DOSTPM / (ED50WTPM + DOSTPM)) * ON;
double FGLPW = (0.055 * DOSEX / (14.0 + DOSEX)
              + 0.150 * DOSSEM / (1.4 + DOSSEM)) * ON;
double WTT, TAUW;
if (FBR > 0.0) {
  WTT  = WT0 * (1.0 - FBR - FTPMW - FGLPW + FSTER);
  TAUW = TAUWTB;
} else {
  WTT  = WT0 * (1.0 - FDT - FTPMW - FGLPW + FSTER);
  TAUW = TAUWT;
}
dxdt_WT  = (WTT - WT) / TAUW;
dxdt_IAP = (KIAP * POSF(BMI - BMIREF) - IAP) / 1.0;
dxdt_PCV = (PCV0 + KPCV * IAP - PCV) / 1.0;

// ------------------------------------------------------ optic nerve head ---
double ONSFT = (SOLVERTIME >= TONSF) ? 1.0 : 0.0;
dxdt_ONSFS = (ONSFT - ONSFS) / 0.5;
dxdt_PONS  = (ICP - PONS) / TAUONS - ONSFS * KFEN * (PONS - PORB);
double IOPT = IOP0 + KIOPBMI * POSF(BMI - 25.0) - EMAXIOP * OCCCA;
dxdt_IOP   = (IOPT - IOP) / TAUIOP;

double TLPD  = PONS - IOP;
double STAST = SIGF((TLPD - TLPDTHR) / SSTAS);
dxdt_STASIS = (STAST - STASIS) / TAUSTAS;
dxdt_RNFLS  = KSW * POSF(TLPD - TLPDTHR) * (1.0 - RNFLS / RNFLSMAX)
            * (AXON / AXON0) - KRES * RNFLS;
dxdt_AXON   = -(KAXL * STASIS * STASIS * STASIS
              + KISCH * POSF(TLPD - TLPDISCH)) * (AXON - AXONFL);
double MDT = -(SMD * (AXON0 - AXON)
             + SREV * (POSF(TLPD - TLPDTHR) / 15.0 > 1.0
                       ? 1.0 : POSF(TLPD - TLPDTHR) / 15.0));
dxdt_MD     = (MDT - MD) / TAUMD;
dxdt_CUMEXP = POSF(TLPD - TLPDTHR);

// ------------------------------------------------------------- headache ---
double DHA = W1HA * POSF(ICP - ICPHTHR) / 10.0
           + W2HA * POSF(GRAD - 4.0) / 10.0
           + W3HA * POSF(CGRP - 1.0);
dxdt_STG  = KSENS * DHA * (1.0 - STG)
          - KDESENS * STG * (1.0 - FMOHDES * MOH);
dxdt_CGRP = (1.0 + KCG * STG + KCGICP * POSF(ICP - 15.0) - CGRP) / TAUCGRP;
dxdt_MOH  = KMOH * ANALG * (1.0 - MOH) - KMOHOFF * MOH * (1.0 - ANALG);
double MHDF = 0.04 + 0.40 * STG + 0.20 * MOH + 0.26 * DHA;
if (MHDF > 1.0) MHDF = 1.0;
dxdt_MHD = (30.0 * MHDF - MHD) / TAUMHD;

// ----------------------------------------------- tinnitus and diplopia ----
double TF = pow(GRAD / TINNREF, 1.5);
if (TF > 1.0) TF = 1.0;
dxdt_TINN = (10.0 * TF - TINN) / TAUTINN;
dxdt_DIPL = KDIPL * POSF(ICP - ICPDIPL) * (1.0 - DIPL) - KDIPLR * DIPL;

// ---------------------------------------------- tolerability / adherence --
dxdt_PARES = ((CACZ / (EC50PAR + CACZ) + 0.45 * CTPM / (EC50COG + CTPM))
              - PARES) / TAUAE;
dxdt_COGSL = ((CTPM / (EC50COG + CTPM)) - COGSL) / TAUAE;
double PAR1 = (PARES > 1.0) ? 1.0 : PARES;
double COG1 = (COGSL > 1.0) ? 1.0 : COGSL;
double ADHT = ADHMAX * (1.0 - KADHPAR * PAR1 - KADHCOG * COG1);
if (ADHT < 0.05) ADHT = 0.05;
dxdt_ADH = (ADHT - ADH) / TAUADH;
dxdt_KSTONE = ((OCCCA * (1.0 - 0.60 * CITRATE)) - KSTONE) / TAUSTONE;

$TABLE
double _bmi = WT / (HT * HT);
double _qv  = QV0 * QVREL;
double _rtot = RTS0 + RTSF * (1.0 - STENT) + RTSC;
double _grad = _qv * _rtot;
double _psss = PCV + _grad;
double _tlpd = PONS - IOP;
double _rnfl = AXON + RNFLS;

capture ICPCM   = ICP * CMH2O;          // opening pressure in cmH2O
capture GRADmm  = _grad;                // transverse sinus gradient (mmHg)
capture PSSSmm  = _psss;
capture FFORMml = FFORM0 * FSEC;
capture TLPDmm  = _tlpd;
capture RNFLum  = _rnfl;                // what OCT reports
capture FRISEN  = (5.0 * RNFLS / (RNFLS + 70.0) > 5.0
                   ? 5.0 : 5.0 * RNFLS / (RNFLS + 70.0));
capture BMIout  = _bmi;
capture CACZmgL = AACZP / VCACZ;
capture CTPMmgL = ATPMC / VCTPM;
capture CEXngmL = AEXC / VCEX;
// intracranial compliance and the LP relaxation time
capture COMPL   = 1.0 / (E1 * ICP);
capture TAURELm = (1.0 / (E1 * ICP)) * ROUT;
'

iih <- mcode("iih", iih_code)

# =============================================================================
# 1.  STEADY STATE: the patient IS a fixed point, so find it before dosing
# =============================================================================
# The single most important practical point about this model: you must not
# start a simulation from an arbitrary ICP.  ICP is the fixed point of
# equations (1)+(2), so the initial condition has to be SOLVED for.  Because
# the loop can be bistable, the solution is not always unique, and which root
# you land on is itself clinically meaningful.
#
# iih_steady() finds every stable root by scanning, exactly as
# steady_icp() does in iih_reference_check.py.

iih_sinus <- function(p, icp, rtsf = 0, stent = 0, pcv, qv) {
  rtscmax <- p$RTSCMAX * (1 - stent)
  rfix    <- p$RTS0 + rtsf * (1 - stent)
  # target(r) - r is strictly decreasing in r, so bisect
  target <- function(r) {
    grad <- qv * (rfix + r)
    pseg <- pcv + (1 - p$FSEG) * grad
    rtscmax / (1 + exp(-(icp - pseg - p$PCRIT) / p$PSTIFF))
  }
  lo <- 0; hi <- rtscmax
  if (rtscmax <= 0) {
    r <- 0
  } else {
    for (i in 1:44) {
      m <- (lo + hi) / 2
      if (target(m) - m > 0) lo <- m else hi <- m
    }
    r <- (lo + hi) / 2
  }
  grad <- qv * (rfix + r)
  list(rtsc = r, grad = grad, psss = pcv + grad,
       pseg = pcv + (1 - p$FSEG) * grad)
}

iih_steady <- function(bmi = 38, pars = list(), n = 2000) {
  p <- as.list(param(iih)); p[names(pars)] <- pars
  iap  <- p$KIAP * max(0, bmi - p$BMIREF)
  pcv  <- p$PCV0 + p$KPCV * iap
  rout <- p$ROUTB + p$KROUTP * pcv
  cort <- 1 + p$KCORTBMI * max(0, bmi - 25)
  egc  <- p$EMAXGC * cort / (p$EC50GC + cort)
  fsec <- (1 + egc) * (1 + p$KANDR * max(0, bmi - 25))
  ff   <- p$FFORM0 * fsec
  G <- function(icp) {
    s <- iih_sinus(p, icp, 0, 0, pcv, p$QV0)
    ff - max(0, icp - s$psss) / rout
  }
  x <- seq(0.5, 140, length.out = n); g <- vapply(x, G, 0)
  roots <- c()
  for (i in seq_len(n - 1)) {
    if (g[i] * g[i + 1] < 0)
      roots <- c(roots, uniroot(G, c(x[i], x[i + 1]), tol = 1e-9)$root)
  }
  stab <- vapply(roots, function(r)
    (G(r + 1e-4) - G(r - 1e-4)) / 2e-4 < 0, TRUE)
  s <- iih_sinus(p, max(roots[stab]), 0, 0, pcv, p$QV0)
  gain <- {
    h <- 1e-3
    (iih_sinus(p, max(roots[stab]) + h, 0, 0, pcv, p$QV0)$psss -
     iih_sinus(p, max(roots[stab]) - h, 0, 0, pcv, p$QV0)$psss) / (2 * h)
  }
  list(roots = roots, stable = roots[stab], nstable = sum(stab),
       icp = max(roots[stab]), icp_cmH2O = max(roots[stab]) * 1.35951,
       grad = s$grad, pcv = pcv, rout = rout, fform = ff,
       g = gain, amplification = 1 / (1 - gain), fsec0 = fsec,
       wt0 = bmi * p$HT^2)
}

# initialise the model at a patient's own untreated fixed point
iih_patient <- function(bmi = 38, pars = list()) {
  s  <- iih_steady(bmi, pars)
  p  <- as.list(param(iih)); p[names(pars)] <- pars
  sn <- iih_sinus(p, s$icp, 0, 0, s$pcv, p$QV0)
  m <- iih %>% param(c(pars, list(WT0 = s$wt0)))
  m <- m %>% init(ICP = s$icp, PONS = s$icp, RTSC = sn$rtsc,
                  ROUT = s$rout, PCV = s$pcv, IAP = s$pcv - p$PCV0,
                  WT = s$wt0, FSEC = s$fsec0, HCO3 = p$HCO30,
                  QVREL = 1, CORTL = 1 + p$KCORTBMI * max(0, bmi - 25),
                  IOP = p$IOP0 + p$KIOPBMI * max(0, bmi - 25),
                  AXON = p$AXON0, CGRP = 1, MHD = 1.5, ADH = p$ADHMAX)
  attr(m, "steady") <- s
  m
}

# =============================================================================
# 2.  SIXTEEN TREATMENT SCENARIOS
# =============================================================================
# Every scenario starts from a 120-day untreated run-in ("presentation"), so
# that papilloedema, central sensitisation and sinus-wall remodelling have the
# values they would actually have on the day the patient is first seen.  ICP
# alone equilibrates in ~30 min; the states that matter clinically do not.

iih_present <- function(mod, days = 120) {
  out <- mod %>% mrgsim(end = days, delta = 1, atol = 1e-8, rtol = 1e-8)
  as.list(tail(as.data.frame(out), 1))
}

iih_run <- function(mod, days = 365, delta = 1, runin = 120, ...) {
  pres <- iih_present(mod, runin)
  cmts <- names(init(mod))
  ic <- pres[cmts]; names(ic) <- cmts
  mod %>% init(unlist(ic)) %>% param(...) %>%
    mrgsim(end = days, delta = delta, atol = 1e-8, rtol = 1e-8)
}

SCENARIOS <- list(

  # --- 1. natural history --------------------------------------------------
  s01_untreated = list(
    label = "untreated natural history",
    pars  = list()),

  # --- 2. the two most-used drugs, at the doses actually used ---------------
  s02_acz_low = list(
    label = "acetazolamide 500 mg/day",
    pars  = list(DOSACZ = 500)),

  s03_acz_std = list(
    label = "acetazolamide 2 g/day (mean tolerated dose in IIHTT)",
    pars  = list(DOSACZ = 2000)),

  s04_acz_max = list(
    label = "acetazolamide 4 g/day (IIHTT target dose)",
    pars  = list(DOSACZ = 4000)),

  # --- 3. IIHTT reconstruction --------------------------------------------
  s05_iihtt_placebo = list(
    label = "IIHTT placebo arm: diet -3.3 %",
    pars  = list(FDIET = 0.033)),

  s06_iihtt_active = list(
    label = "IIHTT active arm: acetazolamide + diet -7 %",
    pars  = list(DOSACZ = 2000, FDIET = 0.070)),

  # --- 4. weight as therapy ------------------------------------------------
  s07_diet = list(
    label = "supervised very-low-energy diet -15.7 % (Sinclair 2010)",
    pars  = list(FDIET = 0.157)),

  s08_bariatric = list(
    label = "bariatric surgery -21 % (IIH:WT)",
    pars  = list(FBAR = 0.21, TBAR = 0)),

  # --- 5. incretin pharmacology -------------------------------------------
  s09_exenatide = list(
    label = "exenatide 10 ug BID (Mitchell 2023)",
    pars  = list(DOSEX = 20)),

  s10_semaglutide = list(
    label = "semaglutide 2.4 mg weekly",
    pars  = list(DOSSEM = 2.4)),

  # --- 6. a mechanistically reasonable target that failed ------------------
  s11_azd4017 = list(
    label = "AZD4017 400 mg BID (11b-HSD1; Markey 2020, negative)",
    pars  = list(DOSAZD = 800)),

  # --- 7. topiramate: the right answer by the wrong route ------------------
  s12_topiramate = list(
    label = "topiramate 150 mg/day",
    pars  = list(DOSTPM = 150)),

  # --- 8. procedures -------------------------------------------------------
  s13_stent = list(
    label = "venous sinus stenting (removes the amplifier)",
    pars  = list(TSTENT = 0)),

  s14_shunt = list(
    label = "CSF shunt, 10 cmH2O valve",
    pars  = list(TSHUNT = 0)),

  s15_onsf = list(
    label = "optic nerve sheath fenestration (decouples the disc)",
    pars  = list(TONSF = 0)),

  # --- 9. the self-defeating bridge ---------------------------------------
  s16_steroid = list(
    label = "prednisolone 40 mg/day bridge (weight gain arm active)",
    pars  = list(DOSPRED = 40, TSTOP = 60)),

  # --- 10. treat-then-withdraw: does remission outlast the drug? ----------
  s17_acz_withdraw = list(
    label = "acetazolamide 2 g/day for 180 d, then stop",
    pars  = list(DOSACZ = 2000, TSTOP = 180)),

  # --- 11. combination -----------------------------------------------------
  s18_combo = list(
    label = "acetazolamide 2 g/day + weight programme -12 %",
    pars  = list(DOSACZ = 2000, FDIET = 0.12))
)

run_scenarios <- function(bmi = 38, days = 365, pars = list()) {
  base <- iih_patient(bmi, pars)
  bind_rows(lapply(names(SCENARIOS), function(nm) {
    sc <- SCENARIOS[[nm]]
    out <- do.call(iih_run, c(list(base, days = days), sc$pars))
    as.data.frame(out) %>% mutate(scenario = nm, label = sc$label)
  }))
}

# =============================================================================
# 3.  ANALYSIS FUNCTIONS
# =============================================================================
# Each of these answers a question that the model can answer and a
# dose-response curve cannot.

## 3.1 -- the loop, decomposed.  Where does ICP come from, term by term?
iih_decompose <- function(bmis = seq(22, 48, by = 2), pars = list()) {
  bind_rows(lapply(bmis, function(b) {
    s <- iih_steady(b, pars)
    data.frame(BMI = b, ICP_cmH2O = s$icp_cmH2O, PCV = s$pcv,
               venous_gradient = s$grad, Rout_Fform = s$rout * s$fform,
               Rout = s$rout, g = s$g, amplification = s$amplification,
               n_stable = s$nstable)
  }))
}

## 3.2 -- identical target engagement, different patients.
##        Only PCRIT (sinus wall support) differs between rows.
iih_responder_spectrum <- function(pcrits = seq(8.5, 16, by = 0.5),
                                   secretion_cut = 0.30, bmi = 38) {
  bind_rows(lapply(pcrits, function(pc) {
    b <- try(iih_steady(bmi, list(PCRIT = pc)), silent = TRUE)
    if (inherits(b, "try-error")) return(NULL)
    # the treated fixed point: same equations, F_form scaled
    p <- as.list(param(iih)); p$PCRIT <- pc
    tgt <- iih_steady(bmi, list(PCRIT = pc, FFORM0 = p$FFORM0 *
                                  (1 - secretion_cut)))
    data.frame(PCRIT = pc, ICP_cmH2O = b$icp_cmH2O, gradient = b$grad,
               g = b$g, amplification = b$amplification,
               dICP_cmH2O = tgt$icp_cmH2O - b$icp_cmH2O,
               n_stable = b$nstable)
  }))
}

## 3.3 -- continuation: find the folds.  Sweep the CSF-side parameter and
##        report how many stable pressures the patient has.
iih_continuation <- function(fmuls = seq(1.30, 0.44, by = -0.02),
                             bmi = 38, pcrit = 9.6) {
  p0 <- as.list(param(iih))$FFORM0
  bind_rows(lapply(fmuls, function(fm) {
    s <- iih_steady(bmi, list(PCRIT = pcrit, FFORM0 = p0 * fm))
    data.frame(f_multiplier = fm, K_mmHg = s$rout * s$fform,
               n_stable = s$nstable,
               lowest_cmH2O = min(s$stable) * 1.35951,
               highest_cmH2O = max(s$stable) * 1.35951)
  }))
}

## 3.4 -- four targets, one equation.  Matched -20 % on each term.
iih_four_targets <- function(bmi = 38) {
  b <- iih_steady(bmi)
  p <- as.list(param(iih))
  rows <- list(
    c("F_form  -20 %  (carbonic anhydrase inhibition)",
      iih_steady(bmi, list(FFORM0 = p$FFORM0 * 0.8))$icp_cmH2O),
    c("R_out   -20 %  (shunt, lymphatic egress)",
      iih_steady(bmi, list(KROUTP = p$KROUTP * 0.8,
                           ROUTB = p$ROUTB * 0.8))$icp_cmH2O),
    c("P_cv    -20 %  (weight, posture, jugular)",
      iih_steady(bmi, list(PCV0 = p$PCV0 - 0.2 * b$pcv))$icp_cmH2O),
    c("R_ts -> 0      (venous sinus stent)",
      iih_steady(bmi, list(RTSCMAX = 1e-9))$icp_cmH2O),
    c("Q_v +15 %      (acidosis vasodilation)",
      iih_steady(bmi, list(QV0 = p$QV0 * 1.15))$icp_cmH2O))
  data.frame(intervention = vapply(rows, `[`, "", 1),
             ICP_cmH2O = as.numeric(vapply(rows, `[`, "", 2)),
             dICP_cmH2O = as.numeric(vapply(rows, `[`, "", 2)) - b$icp_cmH2O)
}

## 3.5 -- acetazolamide's three arms, isolated by counterfactual
iih_acz_arms <- function(bmi = 38, dose = 2000, days = 180) {
  variants <- list(
    "full model"                        = list(),
    "no acidosis-vasodilation (KACID=0)" = list(KACID = 0),
    "no ocular arm (EMAXIOP=0)"          = list(EMAXIOP = 0),
    "no adherence loss (KADHPAR=0)"      = list(KADHPAR = 0),
    "neither penalty arm"                = list(KACID = 0, EMAXIOP = 0))
  bind_rows(lapply(names(variants), function(nm) {
    m <- iih_patient(bmi, variants[[nm]])
    b <- iih_present(m, 120)
    o <- as.data.frame(iih_run(m, days = days, DOSACZ = dose))
    e <- tail(o, 1)
    data.frame(variant = nm,
               dICP_cmH2O = e$ICPCM - b$ICPCM,
               dTLPD = e$TLPDmm - b$TLPDmm,
               dFrisen = e$FRISEN - b$FRISEN,
               dMD = e$MD - b$MD)
  }))
}

## 3.6 -- the acetazolamide dose-response, with adherence in the loop
iih_dose_response <- function(bmi = 38, doses = c(0, 250, 500, 750, 1000,
                                                  1500, 2000, 3000, 4000, 6000),
                              days = 180) {
  m <- iih_patient(bmi); b <- iih_present(m, 120)
  bind_rows(lapply(doses, function(d) {
    e <- tail(as.data.frame(iih_run(m, days = days, DOSACZ = d)), 1)
    data.frame(dose_mg_day = d, Cacz_mgL = e$CACZmgL, adherence = e$ADH,
               HCO3 = e$HCO3, QVREL = e$QVREL, ICP_cmH2O = e$ICPCM,
               IOP = e$IOP, TLPD = e$TLPDmm, Frisen = e$FRISEN, MD = e$MD,
               dICP = e$ICPCM - b$ICPCM)
  }))
}

## 3.7 -- the OCT trap: two eyes converging on one measured thickness
iih_rnfl_trap <- function(bmi = 38, days = 360) {
  m <- iih_patient(bmi)
  treated <- as.data.frame(iih_run(m, days = days, DOSACZ = 2000, FDIET = 0.10))
  untreat <- as.data.frame(iih_run(m, days = days))
  bind_rows(
    treated %>% transmute(time, arm = "treated at presentation",
                          RNFL = RNFLum, axons = AXON, swelling = RNFLS,
                          MD, ICP_cmH2O = ICPCM),
    untreat %>% transmute(time, arm = "untreated",
                          RNFL = RNFLum, axons = AXON, swelling = RNFLS,
                          MD, ICP_cmH2O = ICPCM))
}

## 3.8 -- delay versus dose: the injury endpoint is an integral
iih_delay_vs_dose <- function(bmi = 38, horizon = 730) {
  m <- iih_patient(bmi)
  policies <- list(
    "full therapy, week 0"  = list(DOSACZ = 2000, FDIET = 0.10, TSTART = 0,   TDIET = 0),
    "full therapy, week 2"  = list(DOSACZ = 2000, FDIET = 0.10, TSTART = 14,  TDIET = 14),
    "full therapy, week 4"  = list(DOSACZ = 2000, FDIET = 0.10, TSTART = 28,  TDIET = 28),
    "full therapy, week 8"  = list(DOSACZ = 2000, FDIET = 0.10, TSTART = 56,  TDIET = 56),
    "full therapy, week 12" = list(DOSACZ = 2000, FDIET = 0.10, TSTART = 84,  TDIET = 84),
    "full therapy, week 26" = list(DOSACZ = 2000, FDIET = 0.10, TSTART = 182, TDIET = 182),
    "HALF therapy, week 0"  = list(DOSACZ = 500,  FDIET = 0.05, TSTART = 0,   TDIET = 0),
    "stent, week 0"         = list(TSTENT = 0),
    "stent, week 12"        = list(TSTENT = 84),
    "no therapy"            = list())
  bind_rows(lapply(names(policies), function(nm) {
    e <- tail(as.data.frame(do.call(iih_run,
              c(list(m, days = horizon), policies[[nm]]))), 1)
    data.frame(policy = nm, MD = e$MD, axons = e$AXON, RNFL = e$RNFLum,
               exposure = e$CUMEXP)
  }))
}

## 3.9 -- the lumbar puncture: 30 minutes of pressure, days of symptom
iih_lp <- function(bmi = 38) {
  m <- iih_patient(bmi)
  pres <- iih_present(m, 120)
  cmts <- names(init(m)); ic <- unlist(pres[cmts]); names(ic) <- cmts
  fast <- m %>% init(ic) %>% param(TLP = 0) %>%
    mrgsim(end = 0.3, delta = 0.002, atol = 1e-10, rtol = 1e-10)
  slow <- m %>% init(ic) %>% param(TLP = 0) %>%
    mrgsim(end = 30, delta = 0.5)
  list(minutes = as.data.frame(fast) %>%
         transmute(min = time * 1440, ICP_cmH2O = ICPCM,
                   gradient = GRADmm, TLPD = TLPDmm),
       days = as.data.frame(slow) %>%
         transmute(day = time, ICP_cmH2O = ICPCM, headache_days = MHD,
                   sensitisation = STG))
}

## 3.10 -- headache does not read ICP
iih_headache <- function(bmi = 38, days = 270) {
  bind_rows(lapply(c(0.15, 0.85), function(a) {
    m <- iih_patient(bmi) %>% param(ANALG = a)
    bind_rows(lapply(list(
        "no ICP therapy" = list(),
        "shunt"          = list(TSHUNT = 0),
        "stent"          = list(TSTENT = 0),
        "acetazolamide"  = list(DOSACZ = 2000)), function(pp) {
      e <- tail(as.data.frame(do.call(iih_run,
                c(list(m, days = days), pp))), 1)
      data.frame(analgesic_use = a, ICP_cmH2O = e$ICPCM,
                 headache_days = e$MHD, sensitisation = e$STG,
                 overuse = e$MOH, tinnitus = e$TINN)
    }), .id = "therapy")
  }))
}

## 3.11 -- a virtual cohort: the problem is allocation, not efficacy
iih_cohort <- function(n = 200, seed = 20260728, horizon = 425) {
  set.seed(seed)
  pcrit <- pmax(9.15, pmin(13.5, rnorm(n, 11.0, 0.9)))
  bmi   <- pmax(27, pmin(52, rnorm(n, 38, 5)))
  analg <- runif(n, 0.1, 0.9)
  delay <- pmin(300, rexp(n, 1 / 70))
  policies <- list(
    "no therapy"        = function(d) list(TSTART = 0),
    "acetazolamide 2 g" = function(d) list(DOSACZ = 2000),
    "weight -12 %"      = function(d) list(FDIET = 0.12, TDIET = 0),
    "both"              = function(d) list(DOSACZ = 2000, FDIET = 0.12, TDIET = 0),
    "stent if gradient > 8 mmHg" = function(d) list(TSTENT = 0))
  bind_rows(lapply(seq_len(n), function(i) {
    m <- iih_patient(bmi[i], list(PCRIT = pcrit[i])) %>% param(ANALG = analg[i])
    b <- iih_present(m, max(1, delay[i]))
    cmts <- names(init(m)); ic <- unlist(b[cmts]); names(ic) <- cmts
    m2 <- m %>% init(ic)
    bind_rows(lapply(names(policies), function(nm) {
      pp <- policies[[nm]](delay[i])
      if (nm == "stent if gradient > 8 mmHg" && b$GRADmm <= 8) pp <- list()
      e <- tail(as.data.frame(m2 %>% param(pp) %>%
                mrgsim(end = horizon - max(1, delay[i]), delta = 5)), 1)
      data.frame(id = i, policy = nm, PCRIT = pcrit[i], BMI0 = bmi[i],
                 delay = delay[i], g0 = attr(m, "steady")$g,
                 gradient0 = b$GRADmm,
                 dMD = e$MD - b$MD, dICP = e$ICPCM - b$ICPCM,
                 dHA = e$MHD - b$MHD, axons = e$AXON)
    }))
  }))
}

# =============================================================================
# 4.  PARAMETER PROVENANCE AND CALIBRATION NOTES
# =============================================================================
#
# CSF hydrodynamics
#   FFORM0 0.35 mL/min (500 mL/day) is the classical human CSF formation rate.
#   E1 = 0.2 /mL gives a compliance of 0.5 mL/mmHg at ICP 10 mmHg and 0.19 at
#   ICP 26 mmHg, matching the Marmarou exponential pressure-volume relation and
#   reproducing the clinical observation that ICP pulse amplitude grows as mean
#   pressure rises.
#   R_out is NOT fitted as a free lesion.  ROUTB = 0 and KROUTP = 1.69
#   reproduce a measured R_out of ~8.8 mmHg.min/mL in a lean subject and
#   ~17.1 at BMI 38, which is the range reported by infusion studies in IIH.
#
# Venous loop
#   RTS0 gives the normal 2 mmHg superior-sagittal-sinus-to-jugular gradient.
#   PCRIT = 10.95 mmHg and PSTIFF = 4.14 mmHg were solved (not searched) from
#   three requirements simultaneously: an ICP of 25.9 mmHg (35.1 cmH2O) at
#   BMI 38, a loop gain of 0.53 there, and the resulting gradient of 8.4 mmHg.
#   The model then predicts, without further tuning, a normal opening pressure
#   of 15.7-17.6 cmH2O at BMI 22-24 (Whiteley 2006 reports a mean adult CSF
#   opening pressure near 18 cmH2O) and a gradient rising from 4.5 to 10.5 mmHg
#   across the BMI range.
#
# Trial anchors reproduced (see iih_reference_output.txt for the full table)
#   Sinclair 2010 (BMJ), very-low-energy diet, -15.7 % weight, 3 months:
#       published ICP -8.4 cmH2O ; simulated -7.1 cmH2O
#   Mitchell 2023 (Brain), exenatide 10 ug BID, 12 weeks, no weight change:
#       published ICP -5.6 cmH2O ; simulated -4.9 cmH2O
#   Markey 2020 (Brain Commun), AZD4017 400 mg BID, 12 weeks:
#       published: no significant ICP change ; simulated -0.3 cmH2O
#   IIHTT 2014 (JAMA / Wall), acetazolamide vs placebo, both with diet, 6 mo:
#       published papilloedema grade difference -0.70 ; simulated -0.51
#       published mean deviation difference +0.71 dB ; simulated +0.89 dB
#       published ICP difference -4.4 cmH2O ; simulated -6.1 cmH2O  <- the
#       model over-predicts the between-arm ICP difference by ~1.4x.  This is
#       recorded rather than tuned away: the IIHTT placebo arm fell by nearly
#       10 cmH2O on a 3 % weight loss, which no mechanism in this model
#       reproduces, and the 6-month opening pressure was a single measurement
#       in a subset of participants.
#
# Deliberate near-nulls
#   AZD4017 comes out flat because the adiposity drive on secretion is split
#   into an 11b-HSD1-dependent part (EMAXGC 0.10, blockable) and an
#   androgen-dependent part (KANDR, not blockable by 11b-HSD1 inhibition).
#   85 % target engagement of the smaller arm cannot move the pressure.  The
#   negative trial is therefore a prediction of the structure, not an input.
#   Furosemide is weak for the same kind of reason: NKCC1 inhibition at
#   achievable plasma concentrations removes ~8 % of secretion.
#
# Known biases, kept rather than fitted
#   1. The IIHTT ICP over-prediction above.
#   2. The simulated transverse sinus gradient at presentation (8-10 mmHg) sits
#      at the lower end of the 10-25 mmHg reported in stented cohorts, because
#      those cohorts are selected for severe stenosis.
#   3. The unamplified status of body weight follows from FSEG = 0.5, i.e. from
#      taking the collapse to be driven by the mid-segment sinus pressure.  If
#      the collapse were driven at the jugular end (FSEG = 0) weight would be
#      amplified as well.  ANALYSIS 1 of the Python reference prints the whole
#      FSEG sensitivity instead of hiding it.
#
# =============================================================================
# 5.  QUICK START
# =============================================================================
if (interactive()) {
  m <- iih_patient(bmi = 38)
  print(attr(m, "steady")[c("icp_cmH2O", "grad", "g", "amplification",
                            "nstable")])
  print(iih_decompose())
  print(iih_four_targets())
  print(iih_dose_response())
  print(iih_responder_spectrum())
  out <- run_scenarios(days = 365)
  ggplot(out, aes(time, ICPCM, colour = label)) +
    geom_line() + geom_hline(yintercept = 25, linetype = 2) +
    labs(x = "day", y = "opening pressure (cmH2O)",
         title = "IIH: sixteen scenarios from one fixed-point model") +
    theme_bw() + theme(legend.position = "bottom",
                       legend.title = element_blank())
}
