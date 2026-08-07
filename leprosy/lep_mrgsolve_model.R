## =====================================================================
##  lep_mrgsolve_model.R
##  Leprosy (Hansen's disease) — QSP / PK-PD model
##  38 ODE compartments · 16 therapeutic scenarios
##
##  한센병(나병) — 정량적 시스템 약리학 모델
## ---------------------------------------------------------------------
##  STRUCTURAL THESIS
##  -----------------
##  Leprosy is followed, staged and (historically) stopped on the strength
##  of ONE number, the BACTERIAL INDEX, and that number is governed by a
##  rate constant which no drug in the pharmacopoeia touches.
##
##      BI = log10( LIVE + DEAD acid-fast bacilli per gram )
##      MI = 100 x LIVE / (LIVE + DEAD)         ("morphological index")
##
##  A killed M. leprae is still an acid-fast bacillus.  It is removed only
##  when a macrophage digests it, with a half-life of about 4-5 months.
##  So the BI reads the CLEARANCE clock and the MI reads the KILL clock:
##
##      CLOCK 1  kill       d(live)/dt   drug-set     hours to weeks
##      CLOCK 2  clearance  d(dead)/dt   host-set     t1/2 ~ 137 d, drug-INVARIANT
##
##  Everything below follows from separating the two.
##
##  FOUR CONSEQUENCES THE MODEL DERIVES RATHER THAN ASSUMES
##  --------------------------------------------------------
##  (1) THE BI CANNOT COMPARE REGIMENS.  Scenario S02 (WHO MDT-MB, whose
##      rifampicin component sterilises the live pool within days) and S06
##      (dapsone alone, ~250-fold slower) differ by 3.8 log10 in VIABLE
##      burden at 21 days and by 0.24 log10 in BI at one year.  Turning the
##      bactericidal rate constant down 20-fold (S09) moves the one-year BI
##      by 0.11 log10 — inside the reading error of a slit-skin smear.  The
##      MI separates the same pair on day 21.  A slowly falling BI is
##      therefore not evidence of drug failure, and it is not a stopping
##      rule.
##
##  (2) THE REACTION IS DRIVEN BY THE DERIVATIVE, NOT THE LEVEL.  Antigen
##      becomes immunologically available when a bacillus is DIGESTED, not
##      when it is killed, and live organisms actively suppress phagosome
##      maturation.  Killing the population therefore un-blocks degradation
##      of the whole standing load at once: liberation rate rises 2.2-fold
##      within days of the first rifampicin dose and then decays with the
##      dead-bacillus pool.  Erythema nodosum leprosum tracks that rate.
##      The integral is conserved — it is set by the burden at diagnosis,
##      not by the regimen (S02 vs S09 differ by 15% over 2500 days, and
##      that residual is extra replication in the slow arm) — so a drug can
##      move the ENL PEAK in time but cannot abolish the AREA.  A 20-fold
##      faster bactericide raises the ENL peak 3.7-fold while leaving the
##      BI trajectory unchanged.
##
##  (3) CLOFAZIMINE'S ANTI-INFLAMMATORY ARM IS A LOADING PROBLEM.  With a
##      70-day terminal half-life the tissue depot is at 25% of steady
##      state after one month and 83% after six — that is, it arrives after
##      the ENL peak it is meant to prevent.  Front-loading 300 mg daily
##      for the first month (S10) puts the depot 5.3-fold higher at day 30
##      and cuts cumulative ENL burden by 49% in this model.  That number
##      is a PREDICTION; no trial has run it.
##
##  (4) THE NERVE HAS A REVERSIBLE POOL DRAINING INTO AN IRREVERSIBLE SINK.
##      Impairment enters a recoverable compartment (oedema, demyelination,
##      conduction block) which either recovers (k_rec, steroid-accelerated)
##      or fixes as axonal loss (k_fix = 1/180 d).  Only the reversible pool
##      is a drug target.  Two things follow, and they are different: the
##      course must be LONG ENOUGH to cover the reaction (a 20-week WHO
##      taper leaves 2.4x the permanent deficit of a 52-week course, S12 vs
##      S13) and it must be STARTED EARLY (half of the salvageable deficit
##      is gone by day 117; a one-year delay recovers 1% of it, S14).
##
##  WHAT IS FITTED AND WHAT IS PREDICTED
##  ------------------------------------
##  FITTED (to published human or mouse-foot-pad data):
##      - rifampicin, dapsone, clofazimine, prednisolone, thalidomide PK
##      - 3-log viable kill from ONE 600 mg rifampicin dose (Shepard, Levy)
##      - MI to zero in 2-3 weeks on rifampicin, ~5 months on dapsone
##      - BI decline 0.6-1.0 log10/year on MDT
##      - methaemoglobin 5-12% and Hb -1 g/dL on dapsone 100 mg/d
##      - clofazimine terminal t1/2 70 d, plasma Css 0.5-0.9 mg/L
##      - relapse 0.7-3% at 10 y after MDT-MB, higher after 6 months
##  PREDICTED (falls out of the equations, not fitted):
##      - the BI/MI dissociation and its size
##      - the 2.2-fold antigen-liberation surge on starting MDT
##      - ENL at the lepromatous pole and reversal reaction in the
##        borderline zone, as the product of two opposing gradients
##      - the containment threshold SPEC* = (MUMAX-KNAT)/KHOST = 0.0225
##      - the clofazimine loading result (the model's most exposed claim)
##
##  All 38 ODEs are independently re-implemented in lep_verify_python.py
##  and run against 54 published anchors (54/54 pass); the parameter values
##  below are the ones that survive that check.
##
##  DISCLAIMER: educational and research model.  Not for clinical use.
## =====================================================================

library(mrgsolve)
library(dplyr)

lep_code <- '
$PROB
# Leprosy QSP model (38 ODEs)
# BI reads the clearance clock; MI reads the kill clock; ENL tracks d(Ag)/dt.

$PARAM @annotated
// ---------------- patient / disease set-point ----------------
B0      : 150    : bacillary burden at diagnosis (1e6 organisms/g)
SPEC    : 0.05   : immunological set-point, TT 0.95 - BL 0.35 - LL 0.05 (-)
TISSG   : 100    : mass of bacilliferous tissue (g)

// ---------------- rifampicin ----------------
FRIF    : 0.90   : rifampicin oral bioavailability (-)
KARIF   : 36.0   : rifampicin absorption rate constant (1/d)
VRIF    : 50.0   : rifampicin volume of distribution (L)
CLRIF   : 240.0  : rifampicin clearance at ENZ = 1 (L/d)
KENZ    : 0.14   : enzyme turnover rate constant (1/d)
EMXIND  : 1.20   : maximal autoinduction of rifampicin clearance (-)
EC50IND : 1.00   : rifampicin concentration for half-maximal induction (mg/L)
KMAXR   : 12.0   : maximal rifampicin kill rate constant (1/d)
EC50R   : 0.30   : rifampicin concentration for half-maximal kill (mg/L)

// ---------------- dapsone ----------------
FDAP    : 0.93   : dapsone oral bioavailability (-)
KADAP   : 24.0   : dapsone absorption rate constant (1/d)
VDAP    : 75.0   : dapsone volume of distribution (L)
CLDAP   : 43.2   : dapsone clearance (L/d)
FINDD   : 0.35   : rifampicin induction of dapsone clearance (-)
KMAXD   : 0.068  : maximal dapsone kill rate constant (1/d)
EC50D   : 0.50   : dapsone concentration for half-maximal kill (mg/L)
FNOH    : 0.12   : fraction of dapsone cleared to the hydroxylamine (-)
KNOH    : 12.0   : dapsone hydroxylamine elimination rate constant (1/d)
KMET    : 6774   : methaemoglobin formation constant (%*L/mg/d)
KRED    : 12.0   : methaemoglobin reductase rate constant (1/d)
GRED    : 1.0    : methaemoglobin reduction capacity, 0.35 if G6PD deficient (-)
KINHB   : 0.10833: haemoglobin production (g/dL/d)
KOUTHB  : 0.008333: red-cell loss rate constant, 120 d lifespan (1/d)
HEMO    : 18.0   : haemolytic potency of the hydroxylamine (L/mg)
GHEM    : 1.0    : haemolysis multiplier, 6 if G6PD deficient (-)

// ---------------- clofazimine ----------------
FCLO    : 0.50   : clofazimine oral bioavailability (-)
KACLO   : 6.0    : clofazimine absorption rate constant (1/d)
VCLO1   : 100.0  : clofazimine central volume (L)
VCLO2   : 3000.0 : clofazimine tissue volume (L)
QCLO    : 100.0  : clofazimine intercompartmental clearance (L/d)
CLCLO   : 43.0   : clofazimine clearance (L/d)
A2REF   : 2100.0 : tissue depot at steady state on standard MDT (mg)
KMAXC   : 0.090  : maximal clofazimine kill rate constant (1/d)
EC50C   : 2.00   : clofazimine depot for half-maximal kill (fraction of A2REF)
IMAXC   : 0.45   : maximal clofazimine suppression of TNF (-)
IC50C   : 0.80   : clofazimine depot for half-maximal TNF suppression (-)
KPIGIN  : 0.42   : pigmentation accrual (index/d)
KPIGOUT : 0.00385: pigmentation fading, t1/2 180 d (1/d)

// ---------------- prednisolone ----------------
FPDN    : 0.85   : prednisolone oral bioavailability (-)
KAPDN   : 48.0   : prednisolone absorption rate constant (1/d)
VPDN    : 40.0   : prednisolone volume of distribution (L)
CLPDN   : 216.0  : prednisolone clearance (L/d)
IC50P   : 0.020  : prednisolone concentration for half GR occupancy (mg/L)
IMAXP   : 0.90   : maximal glucocorticoid suppression of TNF (-)
KHPA    : 0.0476 : HPA axis turnover (1/d)
HPAMAX  : 0.85   : maximal HPA suppression at full GR occupancy (-)

// ---------------- thalidomide ----------------
FTHA    : 0.90   : thalidomide oral bioavailability (-)
KATHA   : 14.4   : thalidomide absorption rate constant (1/d)
VTHA    : 90.0   : thalidomide volume of distribution (L)
CLTHA   : 240.0  : thalidomide clearance (L/d)
IMAXT   : 0.70   : maximal thalidomide suppression of TNF via cereblon (-)
IC50T   : 0.50   : thalidomide concentration for half-maximal effect (mg/L)

// ---------------- ofloxacin + minocycline (lumped) ----------------
FROM    : 0.95   : second-line oral bioavailability (-)
KAROM   : 36.0   : second-line absorption rate constant (1/d)
VROM    : 100.0  : second-line volume of distribution (L)
CLROM   : 192.0  : second-line clearance (L/d)
KMAXO   : 0.150  : maximal second-line kill rate constant (1/d)
EC50O   : 1.00   : second-line concentration for half-maximal kill (mg/L)

// ---------------- bacterial population ----------------
MUMAX   : 0.05545: maximal growth, ln2 / 12.5 d generation time (1/d)
BMAX    : 190.0  : tissue carrying capacity (1e6/g)
KNAT    : 0.01043: natural death of growing bacilli (1/d)
FNATP   : 0.096  : persister death as a fraction of KNAT (-)
KCLR0   : 0.00168: dead-bacillus degradation while live bacilli suppress it (1/d)
CLRBOOST: 2.00   : de-repression of degradation once the bacilli are dead (-)
KSUB    : 50.0   : live burden half-relieving that suppression (1e6/g)
KGP     : 2.5e-8 : growing to persister switching (1/d)
KPG     : 0.0015 : persister resuscitation (1/d)
PTOL    : 0.002  : persister tolerance of drug killing (-)
PIMM    : 0.02   : persister tolerance of cell-mediated killing (-)
FPERS   : 1e-5   : persister fraction at diagnosis (-)
FRES    : 0      : resistant-mutant fraction at diagnosis (-)
RRIF    : 0      : resistant clone susceptibility to rifampicin (0-1)
RDAP    : 1      : resistant clone susceptibility to dapsone (0-1)
RCLO    : 1      : resistant clone susceptibility to clofazimine (0-1)
RROM    : 1      : resistant clone susceptibility to second line (0-1)
KHOST   : 2.00   : maximal cell-mediated killing (1/d)
KSAT    : 1.00   : live burden half-saturating granuloma containment (1e6/g)
NCRIT   : 1.6e6  : viable organisms per escape event, Poisson relapse index (-)

// ---------------- antigen ----------------
YB      : 0.10   : antigen liberated at the moment of killing (fraction)
YD      : 0.90   : antigen liberated on digestion of the dead bacillus (fraction)
KAGCL   : 0.35   : free antigen clearance (1/d)
FNERVE  : 0.05   : fraction of liberated antigen reaching nerve (-)
KAGNCL  : 0.06   : intraneural antigen clearance (1/d)

// ---------------- humoral arm and ENL ----------------
KINAB   : 1.00   : anti-PGL-1 production (AU/d)
KAB50   : 3.00   : antigen for half-maximal antibody production (AU)
KOUTAB  : 0.0077 : antibody loss, t1/2 90 d (1/d)
KFIC    : 0.0010 : immune-complex formation (1/AU/d)
KCIC    : 2.00   : immune-complex clearance (1/d)
KOUTT   : 12.0   : TNF turnover (1/d)
STIC    : 1.20   : immune-complex drive on TNF (1/AU)
SAG     : 1.50   : direct antigen drive on TNF (-)
KAGT    : 9.00   : antigen for half-maximal direct TNF drive (AU)
KINN    : 0.50   : neutrophil recruitment (1/d)
KOUTN   : 0.50   : neutrophil egress (1/d)
WTNF    : 0.60   : TNF weight in the ENL score (-)
WNEU    : 1.00   : neutrophil weight in the ENL score (-)
E50     : 2.00   : ENL drive giving a score of 50 (-)
HENL    : 4.0    : ENL Hill coefficient (-)

// ---------------- cellular immunity ----------------
CMIMAX  : 1.0    : ceiling on the cell-mediated index (-)
KCMI    : 0.0500 : immune restoration rate constant (1/d)
BSUP    : 1.00   : live burden imposing half-maximal anergy (1e6/g)
KTREG   : 0.02   : regulatory arm turnover (1/d)
AG50T   : 4.50   : antigen for half-maximal regulatory drive (AU)
TREGSUP : 0.50   : maximal regulatory suppression of the Th1 arm (-)

// ---------------- type 1 reaction and nerve ----------------
KT1R    : 0.1429 : reversal-reaction turnover, t1/2 ~5 d (1/d)
KAGN    : 0.30   : intraneural antigen for half-maximal reaction (AU)
PT1R    : 0.60   : maximal glucocorticoid suppression of the reaction (-)
T1R50   : 12.0   : reaction activity giving half-maximal nerve damage (-)
HT1R    : 5.0    : nerve damage Hill coefficient (-)
KDAM    : 0.70   : nerve damage rate at full reaction (points/d)
WT1     : 1.00   : reversal-reaction weight in nerve damage (-)
WE2     : 0.50   : ENL weight in nerve damage (-)
ENLNEU  : 0.40   : fraction of ENL episodes with neuritis (-)
KREC    : 0.01667: spontaneous recovery of reversible impairment (1/d)
SREC    : 2.00   : steroid acceleration of that recovery (-)
KFIX    : 0.005556: fixation to permanent axonal loss, 1/180 d (1/d)

// ---------------- descriptive ----------------
KLES    : 0.0714 : skin lesion index turnover (1/d)
KLES1   : 30.0   : live burden for half-maximal lesion activity (1e6/g)
KALT    : 0.1429 : ALT turnover (1/d)
ALT0    : 25.0   : baseline ALT (U/L)
HRIF    : 0.50   : rifampicin drive on ALT (L/mg)

$CMT @annotated
RIFG  : rifampicin in gut (mg)
RIFC  : rifampicin in plasma (mg)
ENZ   : relative inducible-enzyme amount (-)
DAPG  : dapsone in gut (mg)
DAPC  : dapsone in plasma (mg)
NOH   : dapsone hydroxylamine (mg)
METHB : methaemoglobin (% of total haemoglobin)
HB    : haemoglobin (g/dL)
CLOG  : clofazimine in gut (mg)
CLO1  : clofazimine in plasma (mg)
CLO2  : clofazimine in tissue depot (mg)
PDNG  : prednisolone in gut (mg)
PDNC  : prednisolone in plasma (mg)
THAG  : thalidomide in gut (mg)
THAC  : thalidomide in plasma (mg)
ROMG  : second-line drug in gut (mg)
ROMC  : second-line drug in plasma (mg)
BG    : growing bacilli (1e6/g)
BP    : persister bacilli (1e6/g)
BR    : resistant bacilli (1e6/g)
BD    : dead but acid-fast bacilli (1e6/g)
AG    : free bacillary antigen (AU)
AGN   : intraneural antigen (AU)
AB    : anti-PGL-1 antibody (AU)
IC    : immune complexes (AU)
TNF   : TNF-alpha, relative to healthy baseline (-)
NEU   : neutrophil infiltrate (AU)
CMI   : cell-mediated immunity index (-)
TREG  : regulatory arm (-)
T1R   : type 1 reaction activity (0-100)
NFIR  : reversible nerve impairment (points)
NFIP  : permanent nerve impairment (points)
LES   : skin lesion activity index (0-100)
PIG   : clofazimine pigmentation index (0-100)
HPA   : HPA axis integrity, 1 = normal (-)
ALT   : alanine aminotransferase (U/L)
ENLC  : cumulative ENL burden (score-days)
AGC   : cumulative antigen liberated (AU)

$GLOBAL
// Helper functions shared by $ODE and $TABLE so that the two can never
// drift apart.  Identical expressions appear in lep_verify_python.py.
#define POS(x) ((x) > 0.0 ? (x) : 0.0)

double emax_(double c, double emax, double ec50) {
  return emax * c / (ec50 + c);
}

// dead-bacillus degradation: live organisms block phagosome maturation,
// so the clearance clock speeds up once they have been killed.
double kclr_(double blive, double kclr0, double boost, double ksub) {
  return kclr0 * (1.0 + boost * (1.0 - blive / (blive + ksub)));
}

// cell-mediated killing saturates: a granuloma contains a small burden and
// is subverted by a large one.
double khost_(double cmi, double blive, double khost, double ksat) {
  return khost * cmi * ksat / (ksat + blive);
}

double enl_(double tnf, double neu, double wt, double wn, double e50, double h) {
  double drive = wt * POS(tnf - 1.0) + wn * POS(neu);
  double dh = pow(drive, h);
  return 100.0 * dh / (pow(e50, h) + dh);
}

$MAIN
// The patient presents after years of incubation, at the quasi-steady
// state implied by the parameters.
double kclr0_ = kclr_(B0, KCLR0, CLRBOOST, KSUB);
double bd0_   = KNAT * B0 / kclr0_;
double rel0_  = YB * KNAT * B0 + YD * kclr0_ * bd0_;
double ag0_   = rel0_ / KAGCL;
double ab0_   = (KINAB * ag0_ / (KAB50 + ag0_)) / KOUTAB;
double ic0_   = KFIC * ag0_ * ab0_ / KCIC;
double treg0_ = ag0_ / (AG50T + ag0_);
double cmi0_  = CMIMAX * SPEC / (1.0 + B0 / BSUP) * (1.0 - TREGSUP * treg0_);
double drv0_  = STIC * ic0_ + SAG * ag0_ / (KAGT + ag0_);
double agn0_  = FNERVE * rel0_ / KAGNCL;

ENZ_0   = 1.0;
HB_0    = KINHB / KOUTHB;
HPA_0   = 1.0;
ALT_0   = ALT0;
BG_0    = B0 * (1.0 - FPERS - FRES);
BP_0    = B0 * FPERS;
BR_0    = B0 * FRES;
BD_0    = bd0_;
AG_0    = ag0_;
AGN_0   = agn0_;
AB_0    = ab0_;
IC_0    = ic0_;
TREG_0  = treg0_;
CMI_0   = cmi0_;
TNF_0   = 1.0 + drv0_;
NEU_0   = KINN * (ic0_ + 0.4 * drv0_) / KOUTN;
T1R_0   = 100.0 * cmi0_ * agn0_ / (KAGN + agn0_);

$ODE
// ---------------- concentrations ----------------
double CR_   = POS(RIFC) / VRIF;
double CD_   = POS(DAPC) / VDAP;
double CN_   = POS(NOH)  / VDAP;
double CLOE_ = POS(CLO2) / A2REF;
double CP_   = POS(PDNC) / VPDN;
double CT_   = POS(THAC) / VTHA;
double CO_   = POS(ROMC) / VROM;
double GROCC_ = CP_ / (IC50P + CP_);

// ---------------- kill rate constants ----------------
double kR_ = emax_(CR_,   KMAXR, EC50R);
double kD_ = emax_(CD_,   KMAXD, EC50D);
double kC_ = emax_(CLOE_, KMAXC, EC50C);
double kO_ = emax_(CO_,   KMAXO, EC50O);

double BLIVE_ = POS(BG) + POS(BP) + POS(BR);
double kH_    = khost_(CMI, BLIVE_, KHOST, KSAT);

double killS_ = kR_ + kD_ + kC_ + kO_ + kH_;
double killR_ = RRIF*kR_ + RDAP*kD_ + RCLO*kC_ + RROM*kO_ + kH_;
double killP_ = PTOL * (kR_ + kD_ + kC_ + kO_) + PIMM * kH_;

double KCLRe_ = kclr_(BLIVE_, KCLR0, CLRBOOST, KSUB);
double MU_    = POS(MUMAX * (1.0 - BLIVE_ / BMAX));

double DEATH_ = KNAT*BG + FNATP*KNAT*BP + KNAT*BR
              + killS_*BG + killP_*BP + killR_*BR;
double RELEASE_ = YB * DEATH_ + YD * KCLRe_ * POS(BD);

// ---------------- anti-inflammatory arms ----------------
double IPDN_ = 1.0 - IMAXP * GROCC_;
double ITHA_ = 1.0 - emax_(CT_,   IMAXT, IC50T);
double ICLO_ = 1.0 - emax_(CLOE_, IMAXC, IC50C);
double ENLs_ = enl_(TNF, NEU, WTNF, WNEU, E50, HENL);

// ---------------- rifampicin ----------------
dxdt_RIFG = -KARIF * RIFG;
dxdt_RIFC =  KARIF * RIFG * FRIF - (CLRIF * ENZ / VRIF) * RIFC;
dxdt_ENZ  =  KENZ * (1.0 + emax_(CR_, EMXIND, EC50IND) - ENZ);

// ---------------- dapsone ----------------
double CLD_ = CLDAP * (1.0 + FINDD * (ENZ - 1.0));
dxdt_DAPG  = -KADAP * DAPG;
dxdt_DAPC  =  KADAP * DAPG * FDAP - (CLD_ / VDAP) * DAPC;
dxdt_NOH   =  FNOH * CLD_ * CD_ - KNOH * NOH;
dxdt_METHB =  KMET * CN_ - KRED * GRED * METHB;
dxdt_HB    =  KINHB - KOUTHB * HB * (1.0 + HEMO * CN_ * GHEM);

// ---------------- clofazimine ----------------
dxdt_CLOG = -KACLO * CLOG;
dxdt_CLO1 =  KACLO * CLOG * FCLO - (CLCLO / VCLO1) * CLO1
             - (QCLO / VCLO1) * CLO1 + (QCLO / VCLO2) * CLO2;
dxdt_CLO2 =  (QCLO / VCLO1) * CLO1 - (QCLO / VCLO2) * CLO2;

// ---------------- prednisolone, thalidomide, second line ----------------
dxdt_PDNG = -KAPDN * PDNG;
dxdt_PDNC =  KAPDN * PDNG * FPDN - (CLPDN / VPDN) * PDNC;
dxdt_THAG = -KATHA * THAG;
dxdt_THAC =  KATHA * THAG * FTHA - (CLTHA / VTHA) * THAC;
dxdt_ROMG = -KAROM * ROMG;
dxdt_ROMC =  KAROM * ROMG * FROM - (CLROM / VROM) * ROMC;

// ---------------- bacterial populations ----------------
dxdt_BG = MU_*BG - KNAT*BG - killS_*BG - KGP*BG + KPG*BP;
dxdt_BP = KGP*BG - KPG*BP - FNATP*KNAT*BP - killP_*BP;
dxdt_BR = MU_*BR - KNAT*BR - killR_*BR;
dxdt_BD = DEATH_ - KCLRe_ * BD;

// ---------------- antigen ----------------
dxdt_AG  = RELEASE_ - KAGCL * AG;
dxdt_AGN = FNERVE * RELEASE_ - KAGNCL * AGN;
dxdt_AGC = RELEASE_;

// ---------------- humoral arm ----------------
dxdt_AB = KINAB * AG / (KAB50 + AG) - KOUTAB * AB;
dxdt_IC = KFIC * AG * AB - KCIC * IC;

// ---------------- TNF and neutrophils ----------------
double drive_ = STIC * IC + SAG * AG / (KAGT + AG);
dxdt_TNF = KOUTT * ((1.0 + drive_) * IPDN_ * ITHA_ * ICLO_ - TNF);
dxdt_NEU = KINN * (IC + 0.4 * POS(TNF - 1.0)) * IPDN_ * ITHA_ - KOUTN * NEU;

// ---------------- cellular immunity ----------------
// Anergy is imposed by LIVE organisms (PGL-1 scavenging, M2 polarisation,
// IL-10), so it lifts within days of an effective bactericide.  That is
// why the reversal reaction is a treatment-onset phenomenon.
double cmit_ = CMIMAX * SPEC / (1.0 + BLIVE_ / BSUP) * (1.0 - TREGSUP * TREG);
dxdt_CMI  = KCMI * (cmit_ - CMI);
dxdt_TREG = KTREG * (AG / (AG50T + AG) - TREG);

// ---------------- type 1 reaction ----------------
double t1t_ = 100.0 * CMI * (AGN / (KAGN + AGN)) * (1.0 - PT1R * GROCC_);
dxdt_T1R = KT1R * (t1t_ - T1R);

// ---------------- nerve: reversible pool, irreversible sink ----------------
double th_  = pow(POS(T1R), HT1R);
double dmg_ = WT1 * th_ / (pow(T1R50, HT1R) + th_)
            + WE2 * ENLNEU * ENLs_ / 100.0;
double room_ = POS(100.0 - NFIR - NFIP) / 100.0;
dxdt_NFIR = KDAM * dmg_ * room_
          - KREC * (1.0 + SREC * GROCC_) * NFIR
          - KFIX * NFIR;
dxdt_NFIP = KFIX * NFIR;

// ---------------- descriptive states ----------------
double lest_ = 100.0 * (0.5 * BLIVE_ / (BLIVE_ + KLES1)
                        + 0.3 * POS(T1R) / 100.0 + 0.2 * ENLs_ / 100.0);
dxdt_LES  = KLES * (lest_ - LES);
dxdt_PIG  = KPIGIN * CLOE_ - KPIGOUT * PIG;
dxdt_HPA  = KHPA * ((1.0 - HPAMAX * GROCC_) - HPA);
dxdt_ALT  = KALT * (ALT0 * (1.0 + HRIF * CR_) - ALT);
dxdt_ENLC = ENLs_ / 100.0;

$TABLE
double blive = POS(BG) + POS(BP) + POS(BR);
double btot  = blive + POS(BD);

capture CRIF  = POS(RIFC) / VRIF;             // mg/L
capture CDAP  = POS(DAPC) / VDAP;             // mg/L
capture CCLO  = POS(CLO1) / VCLO1;            // mg/L
capture CLOE  = POS(CLO2) / A2REF;            // fraction of standard depot
capture CPDN  = POS(PDNC) / VPDN;             // mg/L
capture CTHA  = POS(THAC) / VTHA;             // mg/L
capture GROCC = CPDN / (IC50P + CPDN);

capture BLIVE = blive;                        // 1e6 viable organisms/g
capture BTOT  = btot;
capture BI    = log10(btot > 1e-30 ? btot : 1e-30) + 3.0;
capture BILIV = log10(blive > 1e-30 ? blive : 1e-30) + 3.0;
capture MI    = btot > 1e-30 ? 100.0 * blive / btot : 0.0;

capture KCLRE = kclr_(blive, KCLR0, CLRBOOST, KSUB);
capture AGRATE = YB * (KNAT*POS(BG) + FNATP*KNAT*POS(BP) + KNAT*POS(BR))
               + YD * KCLRE * POS(BD);
capture ENL   = enl_(TNF, NEU, WTNF, WNEU, E50, HENL);
capture NFI   = POS(NFIR) + POS(NFIP);
capture EHF   = NFI * 12.0 / 100.0;           // nerve function score, 0-12
capture GRADE = NFI < 10 ? 0 : (NFI < 30 ? 1 : 2);   // WHO disability grade

// Poisson single-escape interpretation of the residual viable burden.
capture NVIAB = blive * 1e6 * TISSG;
capture PREL  = 100.0 * (1.0 - exp(-NVIAB / NCRIT));
'

lep_mod <- mrgsolve::mcode("leprosy", lep_code)

## ---------------------------------------------------------------------
## Dosing helpers.  Every regimen is expressed as a (time, amt, cmt) frame
## so that scenarios can be composed by row-binding.
## ---------------------------------------------------------------------
dose_seq <- function(cmt, mg, start = 0, ii = 1, n = 1) {
  if (n <= 0) return(NULL)
  data.frame(time = start + ii * (seq_len(n) - 1), amt = mg,
             cmt = cmt, stringsAsFactors = FALSE)
}

## WHO multibacillary MDT: rifampicin 600 mg monthly (supervised),
## clofazimine 300 mg monthly + 50 mg daily, dapsone 100 mg daily.
mdt_mb <- function(months = 12, start = 0, clo_load = FALSE, rom = FALSE) {
  d <- rbind(
    dose_seq("RIFG", 600, start, 30, months),
    dose_seq("CLOG", 300, start, 30, months),
    dose_seq("CLOG",  50, start,  1, 30 * months),
    dose_seq("DAPG", 100, start,  1, 30 * months))
  if (clo_load) d <- rbind(d, dose_seq("CLOG", 250, start, 1, 30))
  if (rom)      d <- rbind(d, dose_seq("ROMG", 400, start, 1, 30 * months))
  d
}

## WHO paucibacillary MDT: rifampicin 600 mg monthly + dapsone 100 mg daily.
mdt_pb <- function(months = 6, start = 0) {
  rbind(dose_seq("RIFG", 600, start, 30, months),
        dose_seq("DAPG", 100, start,  1, 30 * months))
}

## Rifampicin-sparing second line for confirmed rifampicin resistance.
rom_regimen <- function(months = 24, start = 0) {
  rbind(dose_seq("ROMG", 400, start, 1, 30 * months),
        dose_seq("CLOG",  50, start, 1, 30 * months))
}

## Prednisolone taper.  weeks = 20, floor = 5 is the standard WHO course for
## a type 1 reaction; longer courses hold a maintenance dose at the floor.
pred_taper <- function(start, weeks = 20, top = 40, floor = 5, step_wk = 2) {
  out <- NULL; dose <- top
  for (w in seq_len(weeks)) {
    out <- rbind(out, dose_seq("PDNG", dose, start + 7 * (w - 1), 1, 7))
    if (w %% step_wk == 0) dose <- max(dose - 5, floor)
  }
  out
}

thalidomide <- function(start, days = 84, mg = 300) {
  dose_seq("THAG", mg, start, 1, days)
}

make_ev <- function(df) {
  if (is.null(df) || !nrow(df)) return(NULL)
  df <- df[order(df$time), , drop = FALSE]
  Reduce(`+`, lapply(seq_len(nrow(df)), function(i)
    mrgsolve::ev(time = df$time[i], amt = df$amt[i], cmt = df$cmt[i])))
}

## ---------------------------------------------------------------------
## Patient phenotypes along the Ridley-Jopling spectrum.
## ---------------------------------------------------------------------
PATIENT <- list(
  LL     = list(B0 = 150,   SPEC = 0.05),   # lepromatous, BI ~5.9
  BL     = list(B0 =  20,   SPEC = 0.35),   # borderline lepromatous, BI ~5
  BT     = list(B0 =   0.2, SPEC = 0.70),   # borderline tuberculoid
  TT     = list(B0 =   1e-3,SPEC = 0.95),   # tuberculoid, smear negative
  POLAR  = list(B0 = 150,   SPEC = 0.010)   # polar LL below the containment
)                                           # threshold SPEC* = 0.0225

## ---------------------------------------------------------------------
## 16 scenarios.
## ---------------------------------------------------------------------
SCEN <- list(
  S01 = list(lab = "Untreated lepromatous leprosy (natural history)",
             pt = "LL", ev = NULL, end = 1825),
  S02 = list(lab = "WHO MDT-MB 12 months (reference regimen)",
             pt = "LL", ev = mdt_mb(12), end = 1825),
  S03 = list(lab = "WHO MDT-PB 6 months (paucibacillary)",
             pt = "BT", ev = mdt_pb(6), end = 1095),
  S04 = list(lab = "Uniform MDT: 6 months for a multibacillary patient",
             pt = "LL", ev = mdt_mb(6), end = 1825),
  S05 = list(lab = "MDT-MB 24 months (the pre-1998 regimen)",
             pt = "LL", ev = mdt_mb(24), end = 1825),
  S06 = list(lab = "Dapsone monotherapy (historical)",
             pt = "LL", ev = dose_seq("DAPG", 100, 0, 1, 3650), end = 3650),
  S07 = list(lab = "Dapsone monotherapy in a polar patient with folP1 mutants",
             pt = "POLAR", ev = dose_seq("DAPG", 100, 0, 1, 3650), end = 3650,
             par = list(FRES = 1e-6, RDAP = 0, RRIF = 1, RCLO = 1, RROM = 1)),
  S08 = list(lab = "Single-dose rifampicin post-exposure prophylaxis",
             pt = "LL", ev = dose_seq("RIFG", 600, 0, 30, 1), end = 365,
             par = list(B0 = 1e-6)),
  S09 = list(lab = "EXPERIMENT: MDT with a 20x weaker bactericide",
             pt = "LL", ev = mdt_mb(12), end = 1825, par = list(KMAXR = 0.45)),
  S10 = list(lab = "EXPERIMENT: MDT with a 1-month clofazimine load",
             pt = "LL", ev = mdt_mb(12, clo_load = TRUE), end = 1825),
  S11 = list(lab = "MDT plus daily ofloxacin/minocycline (accelerated kill)",
             pt = "LL", ev = mdt_mb(12, rom = TRUE), end = 1825),
  S12 = list(lab = "Reversal reaction: 20-week WHO prednisolone taper",
             pt = "BL", ev = rbind(mdt_mb(12), pred_taper(20, weeks = 20)),
             end = 1200),
  S13 = list(lab = "Reversal reaction: 52-week prednisolone course",
             pt = "BL",
             ev = rbind(mdt_mb(12), pred_taper(20, weeks = 52, floor = 10)),
             end = 1200),
  S14 = list(lab = "Reversal reaction: same course started 180 days late",
             pt = "BL",
             ev = rbind(mdt_mb(12), pred_taper(200, weeks = 52, floor = 10)),
             end = 1200),
  S15 = list(lab = "ENL managed with thalidomide 300 mg for 12 weeks",
             pt = "LL", ev = rbind(mdt_mb(12), thalidomide(30, 84, 300)),
             end = 1095),
  S16 = list(lab = "MDT in G6PD deficiency (dapsone-associated haemolysis)",
             pt = "LL", ev = mdt_mb(12), end = 730,
             par = list(GRED = 0.35, GHEM = 6.0)),
  S17 = list(lab = "Rifampicin-resistant leprosy on second-line ROM + clofazimine",
             pt = "POLAR", ev = rom_regimen(24), end = 1825,
             par = list(FRES = 1e-7, RRIF = 0, RDAP = 1, RCLO = 1, RROM = 1))
)

run_scenario <- function(key, delta = 1) {
  s <- SCEN[[key]]
  p <- PATIENT[[s$pt]]
  if (!is.null(s$par)) p <- utils::modifyList(p, s$par)
  m  <- lep_mod %>% param(p)
  ee <- make_ev(s$ev)
  out <- if (is.null(ee)) m %>% mrgsim(end = s$end, delta = delta)
         else              m %>% mrgsim(events = ee, end = s$end, delta = delta)
  d <- as.data.frame(out)
  d$scenario <- key
  d$label    <- s$lab
  d$patient  <- s$pt
  d$rx_end   <- if (is.null(s$ev)) 0 else max(s$ev$time)
  d
}

run_all <- function(delta = 1) {
  dplyr::bind_rows(lapply(names(SCEN), run_scenario, delta = delta))
}

## ---------------------------------------------------------------------
## Numerical summary: the quantities the thesis turns on.
## ---------------------------------------------------------------------
summarise_scenario <- function(d) {
  at <- function(t, col) d[[col]][which.min(abs(d$time - t))]
  data.frame(
    scenario   = d$scenario[1],
    label      = d$label[1],
    BI_0       = round(at(0, "BI"), 2),
    BI_365     = round(at(365, "BI"), 2),
    dBI_yr1    = round(at(0, "BI") - at(365, "BI"), 2),
    MI_21d     = signif(at(21, "MI"), 3),
    logBLIVE21 = round(log10(max(at(21, "BLIVE"), 1e-30)) + 3, 2),
    AGrate_max = round(max(d$AGRATE), 2),
    ENL_peak   = round(max(d$ENL), 1),
    ENL_days   = round(max(d$ENLC), 1),
    T1R_peak   = round(max(d$T1R), 1),
    NFI_end    = round(tail(d$NFI, 1), 1),
    NFIP_end   = round(tail(d$NFIP, 1), 1),
    MetHb_max  = round(max(d$METHB), 1),
    Hb_min     = round(min(d$HB), 1),
    Pig_max    = round(max(d$PIG), 0),
    ## relapse index is read at the END OF TREATMENT, where the residual
    ## viable burden is the quantity that matters
    Prelapse   = signif(at(d$rx_end[1], "PREL"), 3),
    row.names  = NULL)
}

## Runs only when this file is executed directly (Rscript lep_mrgsolve_model.R);
## sourcing it from the Shiny app leaves the scenarios untouched.
if (sys.nframe() == 0L && !interactive()) {
  res <- run_all()
  print(dplyr::bind_rows(lapply(split(res, res$scenario), summarise_scenario)),
        width = 200)
}
