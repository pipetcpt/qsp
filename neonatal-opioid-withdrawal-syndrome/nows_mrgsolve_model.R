## =============================================================================
##  Neonatal opioid withdrawal syndrome (NOWS / neonatal abstinence syndrome)
##  Quantitative systems pharmacology model — mrgsolve
##  ---------------------------------------------------------------------------
##  THE ORGANISING IDEA
##
##  NOWS is a NORADRENERGIC disease that is treated by titrating a MU-OPIOID
##  and scored with an instrument that measures neither.  Everything below is
##  built on one difference:
##
##       GAP = A - ITONE
##
##       A      neuroadaptive set-point of the locus coeruleus, written by
##              months of continuous in-utero mu-receptor occupancy
##       ITONE  the inhibitory tone sitting on the LC GIRK conductance right now
##
##  GAP > 0  -> LC disinhibition -> noradrenaline efflux -> the signs
##  GAP < 0  -> over-sedation, which the Finnegan sheet scores as ZERO while it
##              suppresses feeding.  The feed volume, not the score, is the
##              over-treatment alarm.
##
##  Five structural commitments follow, and each one produces a testable number:
##
##  (1) ITONE is a UNION, not a sum:  ITONE = EMU + EA2 - EMU*EA2, because
##      mu-opioid and alpha-2 receptors open THE SAME GIRK channel on the LC
##      neuron.  Clonidine therefore substitutes for morphine, sub-additively,
##      with most of its marginal value at the END of a wean.
##      Phenobarbital is deliberately NOT in this union: it acts on GABA-A,
##      downstream, so it lowers the SCORE without lowering the GAP.
##
##  (2) A = AD + AT, two pools with different owners:
##        AD  durable, set in utero, t1/2 ~ 12 d, NOT re-inducible after birth
##        AT  the classical adenylyl-cyclase superactivation pool, tau ~ 60 h,
##            re-driven by ANY mu agonist on board -- including the treatment.
##      AD is the floor under treatment duration.  AT is why suppressing the
##      score harder is not free.
##
##  (3) A0 saturates in occupancy.  A0 = AGAIN*FGA*EAD/(EA50A+EAD) with
##      EAD ~ 0.97 at any clinical maternal methadone dose, so 40 mg/d and
##      160 mg/d differ in A0 by about 1%.  The maternal dose changes WHEN
##      withdrawal starts (log-linear washout, ~1.3 d over a four-fold dose
##      range), not HOW BAD it gets.  This is the model's reading of a
##      literature that has never settled the dose-severity question.
##
##  (4) Buprenorphine enters TWICE and with different numbers: partial agonism
##      (EMAXB = 0.80 for withdrawal suppression) in the acute effect, and poor
##      beta-arrestin recruitment (BADRV = 0.40) in the adaptation drive.  The
##      milder NOWS of buprenorphine-exposed infants is therefore derived from
##      receptor pharmacology rather than fitted as a severity constant.
##
##  (5) A "fixed" mg/kg/day dose is never fixed in effect.  Two hidden trends
##      run in opposite directions: clearance rises with post-menstrual age
##      (UGT2B7 sigmoid, PMA50 54.2 wk) so exposure falls, while the
##      tolerance-shifted EC50 = EC50*(1+KSHIFT*A) falls as A decays so the
##      same exposure buys more effect.
##
##  UNITS: time h, amounts mg, volumes L, concentrations ng/mL, weight kg,
##         doses mg/kg/day (morphine equivalents) unless stated.
##
##  CALIBRATION AND PROVENANCE are documented at the bottom of this file and in
##  README.md.  Every clinical anchor is named with its trial.  Where the model
##  reproduces a number it says so; where it cannot, it says that too.
##
##  NOT FOR CLINICAL USE.
## =============================================================================

library(mrgsolve)

code <- '
$PROB
# Neonatal opioid withdrawal syndrome QSP model (43 states)

$PLUGIN base

$PARAM @annotated
// ---------------- patient / exposure covariates ----------------
WT0     : 3.30  : birth weight (kg)
GA      : 39.0  : gestational age at birth (wk)
MDRUG   : 1     : maternal opioid (0 none, 1 methadone, 2 buprenorphine, 3 short-acting)
MDOSE   : 90.0  : maternal daily dose (mg/day)
FBZD    : 0     : maternal benzodiazepine co-exposure (0/1)
FNIC    : 0     : maternal cigarettes/day
CARE    : 0.55  : non-pharmacologic care bundle index (0-1)
BF      : 0     : breastfeeding (0/1)
ADFORCE : -1    : if >=0, override the durable adaptation pool AD (experiment only)

// ---------------- maturation ----------------
PMA50U  : 54.2  : PMA at half-maximal UGT2B7 capacity (wk)
HILLU   : 3.92  : Hill coefficient for UGT2B7 maturation
PMA50R  : 47.7  : PMA at half-maximal GFR (wk)
HILLR   : 3.40  : Hill coefficient for GFR maturation
PMA50C  : 52.0  : PMA at half-maximal CYP3A4/2B6 (wk)
HILLC   : 3.20  : Hill coefficient for CYP maturation
PNAIND  : 0.35  : fractional postnatal CYP induction
TIND    : 240   : time constant of postnatal CYP induction (h)

// ---------------- morphine (treatment) ----------------
FORAL   : 0.30  : oral bioavailability of morphine solution
KAM     : 1.10  : morphine absorption rate constant (1/h)
CLMREF  : 0.95  : morphine clearance at WT 3.4 kg, PMA 40 wk (L/h)
V1MREF  : 3.40  : morphine central volume at 3.4 kg (L)
V2MREF  : 5.44  : morphine peripheral volume at 3.4 kg (L)
QMREF   : 0.60  : morphine intercompartmental clearance (L/h)
FM6G    : 0.30  : fraction of morphine clearance forming M6G
CLGREF  : 0.10  : M6G renal clearance at reference (L/h)
VGREF   : 1.70  : M6G volume at reference (L)
KE0M    : 0.45  : morphine effect-site equilibration (1/h)
KE0G    : 0.25  : M6G effect-site equilibration (1/h)
RPM     : 1.00  : morphine relative mu potency (reference)
RPG     : 0.035 : M6G effective relative potency (BBB-limited)

// ---------------- methadone ----------------
FDOR    : 0.75  : methadone oral bioavailability
KAD     : 0.35  : methadone absorption rate constant (1/h)
CLDREF  : 0.737 : methadone clearance at reference (L/h)
VDREF   : 17.0  : methadone volume at reference (L)
KE0D    : 0.30  : methadone effect-site equilibration (1/h)
RPD     : 3.00  : methadone relative mu potency
CORD_D  : 0.55  : cord/maternal methadone plasma ratio
KMPD    : 3.60  : maternal methadone conc per mg/day (ng/mL)
RIDD    : 0.025 : relative infant dose via breast milk (methadone)

// ---------------- buprenorphine / norbuprenorphine ----------------
FBSL    : 0.55  : sublingual buprenorphine bioavailability
KAB     : 0.60  : buprenorphine absorption rate constant (1/h)
CLBREF  : 1.60  : buprenorphine clearance at reference (L/h)
VBREF   : 40.0  : buprenorphine volume at reference (L)
KE0B    : 0.25  : buprenorphine effect-site equilibration (1/h)
RPB     : 1.00  : buprenorphine relative potency scaler
EMAXB   : 0.80  : buprenorphine intrinsic activity at mu (partial agonist ceiling)
BADRV   : 0.40  : buprenorphine relative efficacy for DRIVING adaptation
FNORB   : 0.35  : fraction of buprenorphine clearance forming norbuprenorphine
CLNREF  : 0.28  : norbuprenorphine clearance at reference (L/h)
VNREF   : 15.0  : norbuprenorphine volume at reference (L)
KE0N    : 0.20  : norbuprenorphine effect-site equilibration (1/h)
RPN     : 0.30  : norbuprenorphine relative mu potency
CORD_B  : 0.65  : cord/maternal buprenorphine ratio
KMPB    : 0.28  : maternal buprenorphine conc per mg/day (ng/mL)
CNB0    : 1.50  : cord norbuprenorphine:buprenorphine ratio
RIDB    : 0.005 : relative infant dose via breast milk (buprenorphine)

// ---------------- clonidine ----------------
FCLO    : 0.85  : clonidine oral bioavailability
KACL    : 1.50  : clonidine absorption rate constant (1/h)
CLCREF  : 0.80  : clonidine clearance at reference (L/h)
VCREF   : 10.2  : clonidine volume at reference (L)
KE0C    : 0.60  : clonidine effect-site equilibration (1/h)
EC50A2  : 1.60  : clonidine alpha-2 EC50 (ng/mL)
EMAXA2  : 0.34  : maximal alpha-2 inhibitory tone
KTACH   : 0.010 : alpha-2 tachyphylaxis build rate (1/h)
KTOFF   : 0.004 : alpha-2 tachyphylaxis offset rate (1/h)

// ---------------- phenobarbital ----------------
FPHB    : 0.90  : phenobarbital bioavailability
KAP     : 0.80  : phenobarbital absorption rate constant (1/h)
CLPREF  : 0.0212: phenobarbital clearance at reference (L/h)
VPREF   : 3.06  : phenobarbital volume at reference (L)
EC50P   : 22.0  : phenobarbital sedation EC50 (mg/L)
EMAXP   : 0.55  : maximal phenobarbital sedation

// ---------------- mu receptor ----------------
EC50MU  : 12.0  : full-agonist mu EC50 in the opioid-naive neonate (ng/mL morphine-eq)
EC50BU  : 0.22  : buprenorphine EC50 (ng/mL)
EMAXMU  : 1.00  : full-agonist intrinsic activity
KSHIFT  : 3.00  : tolerance shift of EC50 per unit of A
KRDOWN  : 0.004 : mu receptor down-regulation rate (1/h)
KRREC   : 0.006 : mu receptor recovery rate (1/h)

// ---------------- neuroadaptation ----------------
AGAIN   : 1.123 : gain of the in-utero adaptation set-point map
EA50A   : 0.160 : adaptation-drive EC50 for the set-point map
GA50A   : 33.5  : gestational age at half-maximal adaptation capacity (wk)
HGAA    : 8.0   : Hill coefficient for gestational maturity of adaptation
ATMAX   : 0.180 : ceiling of the re-inducible adaptation pool
EAT50   : 0.350 : adaptation-drive EC50 for the re-inducible pool
KATON   : 0.0167: rate constant of the re-inducible pool (1/h, tau 60 h)
KADOFF  : 0.002406 : decay of the durable pool (1/h, t1/2 12 d)
AMAX    : 1.00  : ceiling on the total set-point

// ---------------- locus coeruleus ----------------
LCB0    : 0.70  : baseline LC drive
KENV    : 0.55  : environmental arousal gain on baseline drive
KLC     : 11.9  : gap-to-LC gain
KCG     : 0.30  : fraction of the gap gain abolished by full care
KINNE   : 1.10  : NE production gain
KOUTNE  : 1.00  : NE elimination (1/h)

// ---------------- sign pools ----------------
KINC    : 0.35  : CNS pool input (1/h)
KOUTC   : 0.35  : CNS pool output (1/h)
KINA    : 0.25  : autonomic pool input (1/h)
KOUTA   : 0.25  : autonomic pool output (1/h)
KING    : 0.12  : GI pool input (1/h)
KOUTG   : 0.12  : GI pool output (1/h)
WBZ     : 1.60  : weight of benzodiazepine withdrawal on the CNS pool
WNIC    : 0.45  : weight of nicotine withdrawal on the autonomic pool
KSEDO   : 1.30  : over-sedation produced per unit of negative gap
SCNS    : 1.20  : Finnegan weight of the CNS pool
SANS    : 0.80  : Finnegan weight of the autonomic pool
SGI     : 0.60  : Finnegan weight of the GI pool
KFENV   : 1.55  : Finnegan points contributed by an arousing environment
KFS     : 0.35  : score-smoothing rate, q3-4h scoring as a filter (1/h)

// ---------------- ESC functional items ----------------
THE     : 3.30  : NE at which feeding is half-impaired
SE      : 0.55  : steepness of the feeding response
THS     : 3.00  : NE at which sleep consolidation is half-lost
SS      : 0.55  : steepness of the sleep response
THC     : 3.60  : NE at which consolability is half-lost
SC      : 0.55  : steepness of the consolability response
KCARE   : 1.30  : rightward shift of consolability from caregiver presence
TAUE    : 6.0   : feeding response time constant (h)
TAUS    : 6.0   : sleep response time constant (h)
TAUC    : 4.0   : consolability response time constant (h)

// ---------------- growth ----------------
KCALMAX : 115   : maximal enteral intake (kcal/kg/day)
KCALBASE: 58    : basal expenditure (kcal/kg/day)
KACT    : 6.0   : extra expenditure per unit of NE (kcal/kg/day)
KCALGAIN: 5200  : kcal per kg of weight gain

// ---------------- co-exposures ----------------
KELBZ   : 0.0069: neonatal benzodiazepine elimination (1/h, t1/2 100 h)
KBZON   : 0.020 : benzodiazepine withdrawal onset rate (1/h)
KBZOFF  : 0.005 : benzodiazepine withdrawal offset rate (1/h)
KELNIC  : 0.0058: nicotine withdrawal decay (1/h)

// ---------------- protocol / controller ----------------
FTARGET : 6.0   : score the protocol titrates toward
THRSTART: 8.0   : score at which pharmacotherapy is started
ESCFAIL : 0.95  : ESC functional-failure units that trigger treatment
KUP     : 0.0030: escalation gain (mg/kg/day per point per h)
RUP     : 0.020 : maximum escalation rate (mg/kg/day per h)
DMAX    : 1.00  : maximum dose (mg/kg/day morphine equivalent)
DMIN    : 0.04  : dose below which treatment is stopped
TAPMODE : 1     : 0 = -10%/d of current, 1 = of stabilisation dose, 2 = of standard dose, 3 = A-tracking
TAPFRAC : 0.10  : fractional daily wean step
KTRK    : 0.05  : A-tracking controller gain
GTOL    : 0.10  : gap the A-tracking controller tolerates
STABREQ : 48    : hours of stability required before weaning
MINOBS  : 120   : minimum observation before discharge readiness (h)
KLATCH  : 1.0   : treatment-initiation latch rate (1/h)
KLM     : 0.50  : stabilisation-dose latch rate (1/h)
DINIT   : 0.32  : starting dose (mg/kg/day morphine equivalent)
TRTDRUG : 1     : 1 morphine, 2 methadone, 3 sublingual buprenorphine
CDOSE   : 0     : clonidine dose (ug/kg/day)
PDOSE   : 0     : phenobarbital maintenance (mg/kg/day)
PLOAD   : 0     : phenobarbital loading dose (mg/kg)
ESCMODE : 0     : 0 = Finnegan-driven protocol, 1 = Eat-Sleep-Console

$CMT @annotated
AGUT_M : morphine gut (mg)
AC_M   : morphine central (mg)
AP_M   : morphine peripheral (mg)
CE_M   : morphine effect site (ng/mL)
A_G    : M6G central (mg)
CE_G   : M6G effect site (ng/mL)
AGUT_D : methadone gut (mg)
AC_D   : methadone central (mg)
CE_D   : methadone effect site (ng/mL)
AGUT_B : buprenorphine depot (mg)
AC_B   : buprenorphine central (mg)
CE_B   : buprenorphine effect site (ng/mL)
AC_N   : norbuprenorphine central (mg)
CE_N   : norbuprenorphine effect site (ng/mL)
AGUT_C : clonidine gut (mg)
AC_C   : clonidine central (mg)
CE_C   : clonidine effect site (ng/mL)
TACH   : alpha-2 tachyphylaxis (0-1)
AGUT_P : phenobarbital gut (mg)
AC_P   : phenobarbital central (mg)
AD     : durable neuroadaptation pool
AT     : re-inducible neuroadaptation pool
RMU    : mu receptor availability (0-1)
NE     : locus coeruleus noradrenaline output
CNS    : CNS sign pool
ANS    : autonomic sign pool
GI     : gastrointestinal sign pool
SLP    : sleep consolidation (0-1)
EAT    : feeding effectiveness (0-1)
CONS   : consolability (0-1)
WT     : body weight (kg)
BZD    : neonatal benzodiazepine concentration (relative)
BZW    : GABAergic withdrawal drive
NICW   : nicotine withdrawal drive
CUMM   : cumulative morphine equivalents (mg/kg)
CUMD   : cumulative treatment days
SEIZH  : cumulative seizure hazard
AUCGAP : integral of positive gap
AUCSED : integral of negative gap
FNAS_S : smoothed Finnegan score
DOSE   : current dose (mg/kg/day morphine equivalent)
TRTON  : treatment-initiation latch (0-1)
STAB   : hours of stability accumulated
DSTAB  : latched stabilisation dose (mg/kg/day)

$GLOBAL
#define MATF(x, p50, h) (pow((x),(h)) / (pow((p50),(h)) + pow((x),(h))))

// competitive occupancy of the mu receptor by full agonists and buprenorphine,
// with a tolerance-shifted EC50.  Returns the two fractional occupancies.
struct occ_t { double thf; double thb; double ec50; };
occ_t mu_occ(double UF, double UB, double A,
             double ec50mu, double ec50bu, double kshift) {
  occ_t o;
  o.ec50 = ec50mu * (1.0 + kshift * A);
  double ec50b = ec50bu * (1.0 + kshift * A);
  double xf = UF / o.ec50;
  double xb = UB / ec50b;
  double den = 1.0 + xf + xb;
  o.thf = xf / den;
  o.thb = xb / den;
  return o;
}

$MAIN
// ------------------------------------------------------------------
// Birth: solve the in-utero fixed point.  A, the receptor pool and the
// tolerance shift are mutually dependent, so the set-point at delivery is
// found by iteration rather than written down.
// ------------------------------------------------------------------
double cd0 = 0.0, cb0 = 0.0, cn0 = 0.0;
if (MDRUG == 1) {
  cd0 = CORD_D * KMPD * MDOSE;
} else if (MDRUG == 2) {
  cb0 = CORD_B * KMPB * MDOSE;
  cn0 = CNB0 * cb0;
} else if (MDRUG == 3) {
  cd0 = 40.0;                       // short-acting opioid, morphine-equivalent
}
double FGA = MATF(GA, GA50A, HGAA);
double Ai = 0.5, ri = 0.8, emu_i = 0.0, ead_i = 0.0;
for (int k = 0; k < 500; ++k) {
  double UF0 = RPD * cd0 + RPN * cn0;
  double UB0 = RPB * cb0;
  occ_t o0 = mu_occ(UF0, UB0, Ai, EC50MU, EC50BU, KSHIFT);
  emu_i = ri * (EMAXMU * o0.thf + EMAXB * o0.thb);
  ead_i = ri * (EMAXMU * o0.thf + EMAXB * BADRV * o0.thb);
  ri = KRREC / (KRREC + KRDOWN * emu_i);
  Ai = AGAIN * FGA * ead_i / (EA50A + ead_i);
}
double AT_init = ATMAX * ead_i / (EAT50 + ead_i);
double AD_init = Ai - AT_init;
if (AD_init < 0.0) AD_init = 0.0;
if (ADFORCE >= 0.0) AD_init = ADFORCE;

AD_0  = AD_init;
AT_0  = AT_init;
RMU_0 = ri;
WT_0  = WT0;

double WS0 = WT0 / 3.4;
if (cd0 > 0.0) {
  AC_D_0 = cd0 / 1000.0 * VDREF * WS0 * (MDRUG == 3 ? 0.25 : 1.0);
  CE_D_0 = cd0;
}
if (cb0 > 0.0) {
  AC_B_0 = cb0 / 1000.0 * VBREF * WS0;  CE_B_0 = cb0;
  AC_N_0 = cn0 / 1000.0 * VNREF * WS0;  CE_N_0 = cn0;
}
AC_P_0 = PLOAD * WT0;

// neuro / sign states start at the gap-free equilibrium for this environment
double lc0 = LCB0 * (1.0 + KENV * (1.0 - CARE));
double ne0 = KINNE / KOUTNE * lc0;
NE_0  = ne0;  CNS_0 = ne0;  ANS_0 = ne0;  GI_0 = ne0;
EAT_0  = 1.0 / (1.0 + exp((ne0 - THE) / SE));
SLP_0  = 1.0 / (1.0 + exp((ne0 - THS) / SS));
CONS_0 = 1.0 / (1.0 + exp((ne0 - THC - KCARE * CARE) / SC));
FNAS_S_0 = (SCNS + SANS + SGI) * ne0 + KFENV * (1.0 - CARE);
if (FBZD > 0.0) BZD_0 = 1.0;
if (FNIC > 0.0) NICW_0 = (FNIC / 20.0 < 1.0) ? FNIC / 20.0 : 1.0;

$ODE
double WTx  = (WT > 0.5) ? WT : 0.5;
double PMA  = GA + SOLVERTIME / 168.0;
double SZ   = pow(WTx / 3.4, 0.75);
double WS   = WTx / 3.4;

double RUGT = MATF(PMA, PMA50U, HILLU) / MATF(40.0, PMA50U, HILLU);
double RGFR = MATF(PMA, PMA50R, HILLR) / MATF(40.0, PMA50R, HILLR);
double RCYP = MATF(PMA, PMA50C, HILLC) / MATF(40.0, PMA50C, HILLC)
              * (1.0 + PNAIND * (1.0 - exp(-SOLVERTIME / TIND)));

double CLM = CLMREF * SZ * RUGT;  double V1M = V1MREF * WS;  double V2M = V2MREF * WS;
double QM  = QMREF * SZ;
double CLG = CLGREF * SZ * RGFR;  double VG  = VGREF * WS;
double CLD = CLDREF * SZ * RCYP;  double VD  = VDREF * WS;
double CLB = CLBREF * SZ * RCYP;  double VB  = VBREF * WS;
double CLN = CLNREF * SZ * RGFR;  double VN  = VNREF * WS;
double CLC = CLCREF * SZ * (0.4 + 0.6 * RGFR);  double VC = VCREF * WS;
double CLP = CLPREF * SZ * (0.5 + 0.5 * RCYP);  double VP = VPREF * WS;

double CM  = 1000.0 * AC_M / V1M;   double CPM = 1000.0 * AP_M / V2M;
double CG  = 1000.0 * A_G  / VG;    double CD  = 1000.0 * AC_D / VD;
double CB  = 1000.0 * AC_B / VB;    double CN  = 1000.0 * AC_N / VN;
double CC  = 1000.0 * AC_C / VC;    double CP  = AC_P / VP;   // mg/L

// ---------------- dose delivery (continuous approximation of q3-4h) --------
double DOSEc = (DOSE > 0.0) ? DOSE : 0.0;
double rate_m = 0.0, rate_d = 0.0, rate_b = 0.0;
if (TRTDRUG == 1)       rate_m = DOSEc * WTx / 24.0;
else if (TRTDRUG == 2)  rate_d = DOSEc * WTx / 24.0 / 6.0;    // 1 mg methadone ~ 6 mg morphine
else                    rate_b = DOSEc * WTx / 24.0 / 20.0;   // 1 mg SL bup ~ 20 mg morphine
double rate_c = CDOSE * WTx / 1000.0 / 24.0;
double rate_p = PDOSE * WTx / 24.0;
double bm_d = (MDRUG == 1) ? RIDD * (MDOSE / 70.0) * WTx / 24.0 * BF : 0.0;
double bm_b = (MDRUG == 2) ? RIDB * (MDOSE / 70.0) * WTx / 24.0 * BF : 0.0;

// ---------------- pharmacokinetics ----------------
dxdt_AGUT_M = -KAM * AGUT_M + rate_m;
dxdt_AC_M   = FORAL * KAM * AGUT_M - CLM * CM / 1000.0
              - QM * CM / 1000.0 + QM * CPM / 1000.0;
dxdt_AP_M   = QM * CM / 1000.0 - QM * CPM / 1000.0;
dxdt_CE_M   = KE0M * (CM - CE_M);
dxdt_A_G    = FM6G * 1.62 * CLM * CM / 1000.0 - CLG * CG / 1000.0;
dxdt_CE_G   = KE0G * (CG - CE_G);

dxdt_AGUT_D = -KAD * AGUT_D + rate_d + bm_d;
dxdt_AC_D   = FDOR * KAD * AGUT_D - CLD * CD / 1000.0;
dxdt_CE_D   = KE0D * (CD - CE_D);

dxdt_AGUT_B = -KAB * AGUT_B + rate_b + bm_b;
dxdt_AC_B   = FBSL * KAB * AGUT_B - CLB * CB / 1000.0;
dxdt_CE_B   = KE0B * (CB - CE_B);
dxdt_AC_N   = FNORB * CLB * CB / 1000.0 - CLN * CN / 1000.0;
dxdt_CE_N   = KE0N * (CN - CE_N);

dxdt_AGUT_C = -KACL * AGUT_C + rate_c;
dxdt_AC_C   = FCLO * KACL * AGUT_C - CLC * CC / 1000.0;
dxdt_CE_C   = KE0C * (CC - CE_C);
dxdt_TACH   = KTACH * (CC / (EC50A2 + CC)) * (1.0 - TACH) - KTOFF * TACH;

dxdt_AGUT_P = -KAP * AGUT_P + rate_p;
dxdt_AC_P   = FPHB * KAP * AGUT_P - CLP * CP;

// ---------------- receptor occupancy and the two effects ----------------
double A = AD + AT;  if (A > AMAX) A = AMAX;
double UF = RPM * CE_M + RPG * CE_G + RPD * CE_D + RPN * CE_N;
double UB = RPB * CE_B;
occ_t oc = mu_occ(UF, UB, A, EC50MU, EC50BU, KSHIFT);
double EMU = RMU * (EMAXMU * oc.thf + EMAXB * oc.thb);
double EAD = RMU * (EMAXMU * oc.thf + EMAXB * BADRV * oc.thb);
dxdt_RMU = KRREC * (1.0 - RMU) - KRDOWN * EMU * RMU;

double EA2 = EMAXA2 * CE_C / (EC50A2 * (1.0 + 1.2 * TACH) + CE_C);

// the union: mu and alpha-2 open the same GIRK conductance
double ITONE = EMU + EA2 - EMU * EA2;

// ---------------- neuroadaptation ----------------
dxdt_AD = -KADOFF * AD;
dxdt_AT = KATON * (ATMAX * EAD / (EAT50 + EAD) - AT);

// ---------------- the gap ----------------
double GAP = A - ITONE;
double GP  = (GAP > 0.0) ? GAP : 0.0;
double GN  = (GAP < 0.0) ? -GAP : 0.0;

double SED = EMAXP * CP / (EC50P + CP) + KSEDO * GN;
if (SED > 0.95) SED = 0.95;

double LC = LCB0 * (1.0 + KENV * (1.0 - CARE)) + KLC * GP * (1.0 - KCG * CARE);
dxdt_NE = KINNE * LC - KOUTNE * NE;

// ---------------- non-opioid co-exposures ----------------
dxdt_BZD  = -KELBZ * BZD;
dxdt_BZW  = KBZON * FBZD * (1.0 - BZD) * (1.0 - BZW) - KBZOFF * BZW;
dxdt_NICW = -KELNIC * NICW;

// ---------------- sign domains ----------------
dxdt_CNS = KINC * (NE + WBZ * BZW) * (1.0 - 0.80 * SED) - KOUTC * CNS;
dxdt_ANS = KINA * (NE + WNIC * NICW) * (1.0 - 0.50 * SED) - KOUTA * ANS;
dxdt_GI  = KING * NE * (1.0 - 0.30 * SED) - KOUTG * GI;

double FRAW = SCNS * CNS + SANS * ANS + SGI * GI + KFENV * (1.0 - CARE);
dxdt_FNAS_S = KFS * (FRAW - FNAS_S);

// ---------------- ESC functional items ----------------
double Etar = (1.0 / (1.0 + exp((NE - THE) / SE))) * (1.0 - 0.60 * SED);
double Star = (1.0 / (1.0 + exp((NE - THS) / SS))) * (1.0 - 0.20 * SED) + 0.20 * SED;
double Ctar = 1.0 / (1.0 + exp((NE - THC - KCARE * CARE) / SC));
dxdt_EAT  = (Etar - EAT) / TAUE;
dxdt_SLP  = (Star - SLP) / TAUS;
dxdt_CONS = (Ctar - CONS) / TAUC;

// ---------------- growth ----------------
dxdt_WT = WTx * (KCALMAX * EAT - KCALBASE - KACT * NE) / KCALGAIN / 24.0;

// ---------------- the protocol, as a feedback controller ----------------
double sig, thr, tgt;
if (ESCMODE < 0.5) {
  sig = FNAS_S;  thr = THRSTART;  tgt = FTARGET;
} else {
  double fail = 3.0 - (EAT + SLP + CONS);
  double x = fail - ESCFAIL;  if (x < 0.0) x = 0.0;
  sig = 4.0 + 12.0 * x;  thr = 8.0;  tgt = 6.0;
}
double dTRTON = (sig > thr) ? KLATCH * (1.0 - TRTON) : 0.0;
dxdt_TRTON = dTRTON;

double ERR = sig - tgt;
double dDOSE = 0.0, dSTAB = 0.0;
if (TRTON > 0.5) {
  if (ERR > 0.0) {
    dDOSE = KUP * ERR;  if (dDOSE > RUP) dDOSE = RUP;
    dSTAB = -0.5 * STAB;
  } else {
    dSTAB = 1.0;
    if (STAB > STABREQ) {
      if (TAPMODE < 0.5)        dDOSE = -TAPFRAC * DOSEc / 24.0;
      else if (TAPMODE < 1.5)   dDOSE = -TAPFRAC * DSTAB / 24.0;
      else if (TAPMODE < 2.5)   dDOSE = -TAPFRAC * DINIT / 24.0;
      else {
        // A-tracking: uses the unobservable set-point.  Present as an upper
        // bound on what any protocol could achieve, not as a proposal.
        double atar = A - GTOL;  if (atar < 0.0) atar = 0.0;
        double head = EMAXMU * RMU - atar;  if (head < 1e-3) head = 1e-3;
        double ctar = oc.ec50 * atar / head;
        dDOSE = KTRK * (ctar / 44.7 - DOSEc);
      }
    }
  }
}
dDOSE += DINIT * dTRTON;                        // the starting dose is the latch step
if (DOSEc <= DMIN && dDOSE < 0.0) dDOSE = 0.0;
if (DOSEc >= DMAX && dDOSE > 0.0) dDOSE = 0.0;
dxdt_DOSE  = dDOSE;
dxdt_STAB  = dSTAB;
dxdt_DSTAB = (DOSEc > DSTAB) ? KLM * (DOSEc - DSTAB) : 0.0;

// ---------------- counters ----------------
dxdt_CUMM   = (rate_m + 6.0 * rate_d + 20.0 * rate_b) / WTx;
dxdt_CUMD   = (DOSEc > DMIN) ? 1.0 / 24.0 : 0.0;
double gz   = GP - 0.55;  if (gz < 0.0) gz = 0.0;
dxdt_SEIZH  = 0.004 * gz * gz;
dxdt_AUCGAP = GP;
dxdt_AUCSED = GN;

$TABLE
double A_out    = AD + AT;  if (A_out > AMAX) A_out = AMAX;
double UF_o     = RPM * CE_M + RPG * CE_G + RPD * CE_D + RPN * CE_N;
double UB_o     = RPB * CE_B;
occ_t oco       = mu_occ(UF_o, UB_o, A_out, EC50MU, EC50BU, KSHIFT);
double EMU_o    = RMU * (EMAXMU * oco.thf + EMAXB * oco.thb);
double EA2_o    = EMAXA2 * CE_C / (EC50A2 * (1.0 + 1.2 * TACH) + CE_C);
double ITONE_o  = EMU_o + EA2_o - EMU_o * EA2_o;
double GAP_o    = A_out - ITONE_o;
double CP_o     = AC_P / (VPREF * (WT / 3.4));
double CM_o     = 1000.0 * AC_M / (V1MREF * (WT / 3.4));
double CD_o     = 1000.0 * AC_D / (VDREF  * (WT / 3.4));
double CB_o     = 1000.0 * AC_B / (VBREF  * (WT / 3.4));
double CC_o     = 1000.0 * AC_C / (VCREF  * (WT / 3.4));
double CG_o     = 1000.0 * A_G  / (VGREF  * (WT / 3.4));
double ESCFAILo = 3.0 - (EAT + SLP + CONS);
double ONTRT    = (DOSE > DMIN) ? 1.0 : 0.0;
double ESCPASS  = (EAT > 0.55 && SLP > 0.55 && CONS > 0.55) ? 1.0 : 0.0;

$CAPTURE @annotated
A_out    : total neuroadaptive set-point
EMU_o    : mu-mediated inhibitory effect
EA2_o    : alpha-2 mediated inhibitory effect
ITONE_o  : total inhibitory tone on the LC GIRK conductance
GAP_o    : the gap (A - ITONE); positive = withdrawal, negative = over-sedation
CM_o     : morphine plasma concentration (ng/mL)
CD_o     : methadone plasma concentration (ng/mL)
CB_o     : buprenorphine plasma concentration (ng/mL)
CC_o     : clonidine plasma concentration (ng/mL)
CG_o     : M6G plasma concentration (ng/mL)
CP_o     : phenobarbital plasma concentration (mg/L)
ESCFAILo : Eat-Sleep-Console functional failure units (0-3)
ONTRT    : on pharmacotherapy (0/1)
ESCPASS  : all three ESC functions passing (0/1)
'

mod <- mcode("nows", code)

## =============================================================================
##  SIMULATION HELPERS
## =============================================================================

sim_nows <- function(model = mod, tmax = 45 * 24, delta = 1, ...) {
  p <- list(...)
  m <- model
  if (length(p)) m <- param(m, p)
  mrgsim(m, end = tmax, delta = delta, atol = 1e-9, rtol = 1e-6,
         maxsteps = 500000) |> as.data.frame()
}

## endpoints exactly as they are defined in the trials
nows_endpoints <- function(out, dmin = 0.04, minobs = 120) {
  on   <- out$DOSE > dmin
  trt  <- any(out$TRTON > 0.5)
  days <- sum(on) * (out$time[2] - out$time[1]) / 24
  stop_t <- if (any(on)) max(out$time[on]) else 0
  ok   <- out$ESCPASS > 0.5 & !on & out$time >= minobs
  rl   <- rle(ok)
  ready <- NA_real_
  idx <- cumsum(rl$lengths)
  hit <- which(rl$values & rl$lengths >= 48)
  if (length(hit)) ready <- out$time[idx[hit[1]]] / 24
  data.frame(
    A0            = out$A_out[1],
    peak_finnegan = max(out$FNAS_S),
    day_of_peak   = out$time[which.max(out$FNAS_S)] / 24,
    onset_day     = if (any(out$FNAS_S > 8)) out$time[which(out$FNAS_S > 8)[1]] / 24 else NA,
    treated       = trt,
    max_dose      = max(out$DOSE),
    treat_days    = days,
    stop_day      = stop_t / 24,
    ready_day     = ready,
    cum_morphine  = tail(out$CUMM, 1),
    auc_gap       = tail(out$AUCGAP, 1),
    auc_sedation  = tail(out$AUCSED, 1),
    weight_gain_g = (tail(out$WT, 1) - out$WT[1]) * 1000
  )
}

## =============================================================================
##  SCENARIOS
##  Each one is a question the model was built to answer, not a demonstration.
## =============================================================================

scenarios <- list(

  ## --- 1. What the disease looks like with nothing done ---------------------
  s01_untreated_methadone = list(
    label = "Methadone-exposed term infant, no pharmacotherapy",
    par = list(MDRUG = 1, MDOSE = 90, THRSTART = 999)),

  ## --- 2. The reference patient under a Finnegan-driven protocol ------------
  s02_methadone_morphine = list(
    label = "Methadone-exposed, oral morphine, Finnegan protocol, usual care",
    par = list(MDRUG = 1, MDOSE = 90, CARE = 0.55)),

  ## --- 3-5. The maternal-dose question --------------------------------------
  ## The literature has never agreed whether maternal methadone dose predicts
  ## NOWS severity.  The model says: the dose sets the ONSET TIME, not the
  ## severity, because the adaptation map is saturated at every clinical dose.
  s03_maternal_dose_40  = list(label = "Maternal methadone 40 mg/day",
                               par = list(MDRUG = 1, MDOSE = 40)),
  s04_maternal_dose_90  = list(label = "Maternal methadone 90 mg/day",
                               par = list(MDRUG = 1, MDOSE = 90)),
  s05_maternal_dose_160 = list(label = "Maternal methadone 160 mg/day",
                               par = list(MDRUG = 1, MDOSE = 160)),

  ## --- 6-7. Buprenorphine exposure (MOTHER trial) ---------------------------
  s06_bup_exposed = list(
    label = "Buprenorphine-exposed 16 mg/day, oral morphine",
    par = list(MDRUG = 2, MDOSE = 16)),
  s07_bup_exposed_untreated = list(
    label = "Buprenorphine-exposed, no pharmacotherapy",
    par = list(MDRUG = 2, MDOSE = 16, THRSTART = 999)),

  ## --- 8. Short-acting opioid: the gap opens before anyone is watching ------
  s08_short_acting = list(
    label = "Heroin / short-acting opioid exposure",
    par = list(MDRUG = 3)),

  ## --- 9-10. The preterm experiment ----------------------------------------
  ## Preterm NOWS is milder for TWO reasons that are always confounded in
  ## cohort data: a less mature adaptation capacity, and an immature clearance
  ## that makes weight-based dosing relatively generous.  Scenario 10 forces a
  ## 34-week infant to carry a TERM-level durable adaptation pool, isolating
  ## the pharmacokinetic half of the explanation.
  s09_preterm_34 = list(
    label = "Preterm 34 wk, own (lower) adaptation",
    par = list(GA = 34, WT0 = 2.10)),
  s10_preterm_34_termA = list(
    label = "Preterm 34 wk carrying a TERM adaptation pool",
    par = list(GA = 34, WT0 = 2.10, ADFORCE = 0.578)),

  ## --- 11-13. The four weaning rules ---------------------------------------
  ## Not better and worse: four points on one trade-off curve between
  ## treatment days and cumulative withdrawal burden.
  s11_wean_current = list(label = "Wean -10%/day of the CURRENT dose",
                          par = list(TAPMODE = 0)),
  s12_wean_standard = list(label = "Wean -10%/day of a STANDARD dose",
                           par = list(TAPMODE = 2)),
  s13_wean_Atrack = list(label = "Wean by tracking A (unobservable oracle)",
                         par = list(TAPMODE = 3, GTOL = 0.10)),

  ## --- 14. How hard to suppress: the U-shaped harm curve --------------------
  s14_suppress_hard = list(label = "Titrate to a Finnegan score of 3",
                           par = list(FTARGET = 3.0)),
  s15_suppress_light = list(label = "Titrate to a Finnegan score of 7.8",
                            par = list(FTARGET = 7.8)),

  ## --- 16-17. Clonidine: the same channel, a different receptor -------------
  s16_clonidine_adjunct = list(label = "Morphine + clonidine 6 ug/kg/day",
                               par = list(CDOSE = 6)),
  s17_clonidine_high = list(label = "Morphine + clonidine 12 ug/kg/day",
                            par = list(CDOSE = 12)),

  ## --- 18-20. Phenobarbital: score versus gap -------------------------------
  s18_pheno_pure_opioid = list(
    label = "Pure opioid exposure + phenobarbital (score falls, gap does not)",
    par = list(PDOSE = 5, PLOAD = 20)),
  s19_polysubstance = list(
    label = "Benzodiazepine co-exposure, morphine only",
    par = list(FBZD = 1)),
  s20_polysubstance_pheno = list(
    label = "Benzodiazepine co-exposure + phenobarbital",
    par = list(FBZD = 1, PDOSE = 5, PLOAD = 20)),

  ## --- 21-23. Care, breastfeeding, Eat-Sleep-Console ------------------------
  s21_rooming_in = list(label = "Rooming-in / full care bundle, Finnegan protocol",
                        par = list(CARE = 0.90)),
  s22_esc = list(label = "Eat-Sleep-Console with full care bundle",
                 par = list(CARE = 0.90, ESCMODE = 1)),
  s23_breastfeeding = list(label = "Breastfeeding on methadone",
                           par = list(BF = 1)),

  ## --- 24-25. Alternative treatment agents ---------------------------------
  s24_treat_bup = list(label = "Sublingual buprenorphine as treatment (Kraft protocol)",
                       par = list(TRTDRUG = 3)),
  s25_treat_methadone = list(label = "Oral methadone as treatment",
                             par = list(TRTDRUG = 2))
)

run_all <- function() {
  res <- lapply(names(scenarios), function(nm) {
    s <- scenarios[[nm]]
    o <- do.call(sim_nows, s$par)
    cbind(scenario = nm, label = s$label, nows_endpoints(o))
  })
  do.call(rbind, res)
}

## =============================================================================
##  CALIBRATION AND PROVENANCE
##  ---------------------------------------------------------------------------
##  PHARMACOKINETICS
##   - Morphine in the term neonate: CL 0.28 L/h/kg, Vss 2.6 L/kg, terminal
##     t1/2 6-9 h; oral bioavailability ~0.30.  UGT2B7 maturation uses the
##     PMA-driven sigmoid (PMA50 54.2 wk, Hill 3.92) that the neonatal morphine
##     population analyses converged on; at PMA 40 wk this is 23% of the adult
##     per-kg^0.75 capacity, and it rises ~31% over four weeks at term and ~44%
##     over four weeks from 34 wk.  That difference is HIDDEN TAPER #1.
##   - Methadone in the neonate: t1/2 16-25 h, i.e. per-kg clearance roughly a
##     sixth of the adult value despite a similar half-life.  Cord/maternal
##     ratio ~0.5; maternal plasma ~3.6 ng/mL per mg/day.
##   - Buprenorphine and norbuprenorphine: sublingual bioavailability ~0.55 in
##     the neonate, cord/maternal ~0.65, cord norbuprenorphine exceeding
##     buprenorphine.  Buprenorphine EC50 is set 20-fold below morphine on a
##     molar-equivalent scale to reflect its affinity, with Emax capped at 0.60.
##   - M6G: formed by the same UGT2B7 step, cleared renally with a maturation
##     sigmoid of its own, neonatal t1/2 ~12 h.  Effective mu potency is set to
##     3.5% of morphine because the glucuronide crosses the blood-brain barrier
##     poorly -- but it is the LAST mu signal to leave, which is why the
##     terminal hours of a wean are not as bare as the parent drug suggests.
##   - Clonidine: neonatal t1/2 ~8-12 h; EC50 1.6 ng/mL with Emax 0.34, chosen
##     so that 6-12 ug/kg/day is an ADJUNCT and not a substitute for morphine.
##   - Phenobarbital: neonatal t1/2 ~100 h, sedation EC50 22 mg/L.
##
##  CLINICAL ANCHORS THE MODEL IS ASKED TO REPRODUCE
##   - Onset: short-acting <24 h, buprenorphine 36-60 h, methadone 48-72 h,
##     peak in the first week.
##   - Untreated peak Finnegan score in a methadone-exposed term infant:
##     mid-teens.
##   - Treatment concentrations of morphine during NOWS therapy: 15-45 ng/mL.
##   - MOTHER (Jones 2010): buprenorphine-exposed neonates needed less morphine
##     and a shorter course than methadone-exposed neonates, with a similar
##     proportion treated (47% vs 57%).
##   - Kraft 2017: sublingual buprenorphine shortened treatment versus oral
##     morphine (median 15 vs 28 days).
##   - Agthe 2009: clonidine as an adjunct shortened treatment (11 vs 15 days).
##   - ESC-NOW (Young 2023): Eat-Sleep-Console cut pharmacotherapy from 52.0%
##     to 19.5% and time to medically-ready-for-discharge from 14.9 to 8.2 days.
##
##  WHAT THE MODEL DOES NOT CLAIM
##   - It does not predict neurodevelopmental outcome.  The AUCGAP and AUCSED
##     integrals are offered as exposure metrics, not as outcome predictors;
##     the observational literature linking NOWS to later development is too
##     confounded by the social environment to calibrate against.
##   - The absolute treatment duration is a property of the WEANING RULE, not
##     of the infant.  The model reproduces the 3-fold spread in published
##     durations by changing TAPFRAC and TAPMODE alone, at fixed biology; it
##     should not be read as predicting any single centre s length of stay.
##   - The A-tracking controller (TAPMODE 3) is an oracle: it reads a state no
##     bedside measurement provides.  It is included to bound what is
##     achievable, and the notable finding is how LITTLE it buys.
## =============================================================================
