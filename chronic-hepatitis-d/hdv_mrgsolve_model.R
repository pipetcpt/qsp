## =============================================================================
##  CHRONIC HEPATITIS D (HDV) -- QSP / PK-PD MODEL  (mrgsolve)
##  Bulevirtide - peg-IFN alfa/lambda - lonafarnib/ritonavir - HBsAg-siRNA
## =============================================================================
##
##  FILE          hdv_mrgsolve_model.R
##  COMPANIONS    hdv_qsp_model.dot/.svg/.png   mechanistic map (18 clusters)
##                hdv_reference_model.py        pure-Python reference + report
##                hdv_model_report.txt          the numbers this file must match
##                hdv_shiny_app.R               interactive dashboard
##                hdv_references.md             source for every parameter
##
##  RUN
##      library(mrgsolve); library(dplyr); library(ggplot2)
##      mod <- mread("hdv_mrgsolve_model.R")
##      out <- scenario_all(mod)        # all 13 scenarios, bottom of this file
##
## -----------------------------------------------------------------------------
##  THE ONE IDEA
## -----------------------------------------------------------------------------
##  HDV is not "a virus you kill".  It is a TWO-SUBSTRATE ASSEMBLY LINE:
##
##      substrate 1  the HDV genome    made by the HOST (RNA Pol II, double
##                                     rolling circle) -- entirely independent
##                                     of HBV replication, which is why a
##                                     nucleos(t)ide analogue does nothing
##      substrate 2  the HBsAg envelope made by HBV cccDNA, in the same cell or
##                                     supplied by the surrounding HBsAg+ pool
##
##  So every drug is classified by WHICH FLUX it cuts, not by potency:
##
##      (E) entry / re-infection       NTCP-dependent      bulevirtide
##      (A) envelopment / egress       FTase, HBsAg        lonafarnib, siRNA
##      (C) intracellular + cell loss  Pol II, immunity    peg-IFN alfa/lambda
##
##  and the infected-hepatocyte pool Id is ALSO held up by two fluxes that no
##  entry inhibitor can touch:
##
##      division-mediated spread   an Id cell that divides gives two Id cells
##                                 (HDV replicates through mitosis; Giersch 2019)
##      cell-to-cell spread        only partly NTCP-dependent (fraction PHINT)
##
##  That pair is the FLOOR.  The central quantitative claim of this model is
##  that the FLOOR, not target occupancy, is what limits bulevirtide.
##
## -----------------------------------------------------------------------------
##  WHAT IS FITTED (4 numbers) AND WHAT IS PREDICTED (everything else)
## -----------------------------------------------------------------------------
##  FITTED
##    ROATP, KDNTCP     <- the two published total-bile-acid fold-rises on
##                         bulevirtide 2 mg (x3.2) and 10 mg (x13.0).
##                         NTCP OCCUPANCY IS THEN DERIVED, NEVER ASSUMED:
##                            2 mg  -> Css 1.72 nM -> occupancy 0.710
##                           10 mg  -> Css 14.5 nM -> occupancy 0.954
##                         and the fitted KDNTCP of 0.70 nM lands inside the
##                         reported in-vitro NTCP affinity window for
##                         myrcludex B / bulevirtide WITHOUT having been shown
##                         a single binding experiment.
##    DIMM              <- MYR301 virologic response rate at week 48 on 2 mg (71%)
##                         fitted value 0.02636 /day.  NOTE the sign: a LARGER
##                         baseline killing rate gives a SHALLOWER on-treatment
##                         decline, because the untreated steady state must
##                         still balance, so the re-infection flux holding it up
##                         is larger and the residual flux left after blocking a
##                         fixed fraction of entry is larger too.
##    ALTBASE           <- MYR301 ALT normalisation rate at week 48 on 2 mg (51%)
##
##  BACK-SOLVED (not fitted): every parameter that can be pinned by requiring
##    d/dt = 0 at the observed untreated chronic steady state is pinned that
##    way, so the untreated arm is a steady state BY CONSTRUCTION.  See
##    hdv_reference_model.py::back_solve().
##
##  PREDICTED / HELD OUT: the 10 mg arm, the untreated control arm, week 96,
##    peg-IFN on- and off-treatment response, bulevirtide + peg-IFN synergy,
##    lonafarnib/ritonavir, peg-IFN lambda, the entry-inhibition ceiling, the
##    bile-acid therapeutic-index curve, and the fibrosis / HCC projections.
##
## -----------------------------------------------------------------------------
##  KEY RESULTS THIS FILE REPRODUCES (all computed in hdv_model_report.txt)
## -----------------------------------------------------------------------------
##  1. BULEVIRTIDE HAS NO FIRST PHASE, AND THAT IS STRUCTURAL.  An entry
##     inhibitor cannot touch a cell that is already infected, so it cannot
##     lower secretion on day 1; interferon (replication + production block)
##     and lonafarnib (assembly block) both can.  The model was never told
##     this -- it follows from which flux each drug cuts.
##  2. 2 mg DOES NOT SATURATE NTCP (occupancy 0.710, i.e. 29% of receptor left
##     free) and yet 10 mg buys almost nothing on the combined endpoint,
##     because entry supplies only ~53% of the inflow maintaining Id and the
##     rest is division-mediated and cell-to-cell spread.
##  3. THE CEILING OF THE WHOLE ENTRY-INHIBITOR CLASS is a few points above
##     the approved 2 mg dose.  Dose escalation is not the lever.
##  4. ALT AND HDV RNA READ DIFFERENT THINGS.  ALT tracks the INFLOW of newly
##     infected hepatocytes (a term that falls within days of blocking entry);
##     HDV RNA tracks the SIZE of the infected pool (which decays over months).
##     That is the model's explanation for ALT normalisation in virologic
##     non-responders -- and it is a HYPOTHESIS with a stated falsifier.
##  5. THE LONAFARNIB PARADOX.  Blocking farnesyltransferase stops envelopment,
##     not replication, so INTRACELLULAR HDV RNA RISES (~+50%) while serum HDV
##     RNA falls.  Serum RNA measures the assembly flux, not the reservoir, so
##     withdrawal releases the trapped pool and rebound is fast.
##  6. DURABILITY COMES FROM EXHAUSTION REVERSAL, NOT FROM SUPPRESSION.
##     Interferon is the only current mechanism that raises the DEATH RATE of
##     infected cells and reverses exhaustion, so it is the only one whose
##     response survives withdrawal.  Bulevirtide is a suppressive therapy by
##     construction, however deep its on-treatment response.
##  7. AND ONE HONEST DISAGREEMENT: constrained by the bile-acid anchors, the
##     model separates 2 mg from 10 mg on the VIROLOGIC endpoint much more than
##     MYR301 did.  Three candidate resolutions, one of them testable, are set
##     out in hdv_model_report.txt section A13.  It is reported, not absorbed.
##
## -----------------------------------------------------------------------------
##  UNITS
##      time      days
##      pools     T, Ib, Id are FRACTIONS of the hepatocyte population (sum = 1)
##      Rg        intracellular HDV genomic RNA, copies per infected cell
##      Vd, Sser  serum HDV RNA and HBsAg, IU/mL
##      Cblv      nM     Cifn  ng/mL     Clnf, Crtv  uM
##      TBA       umol/L        ALT  U/L        Fib  Ishak 0-6
## =============================================================================

$PROB
# Chronic hepatitis D QSP model: entry (bulevirtide) x envelopment (lonafarnib,
# siRNA) x intracellular + infected-cell loss (peg-IFN), 33 ODE compartments.

$PARAM @annotated
// ---------------------------------------------------------------- bulevirtide
MWBLV   : 5398.9  : bulevirtide molecular weight (Da)
KABLV   : 20.0    : SC absorption rate constant (1/day)
FBLV    : 0.85    : SC bioavailability (-)
VCBLV   : 8.0     : central volume (L)
CLBLV   : 90.0    : linear clearance (L/day)
VMBLV   : 300.0   : saturable target-mediated clearance Vmax (nmol/day)
KMBLV   : 1.5     : saturable clearance Km (nM)
KOUTLIV : 2.0     : degradation of internalised drug (1/day)
KDNTCP  : 0.700652 : FITTED NTCP dissociation constant (nM)
OCCFIX  : -1.0    : if >=0, override occupancy (use 1.0 for a perfect blocker)

// ------------------------------------------------ bile acids = occupancy read-out
TBA0    : 5.0     : baseline total serum bile acids (umol/L)
CLNTCPB : 1.0     : NTCP-mediated bile acid uptake clearance (1/day, normalised)
ROATP   : 0.0333097 : FITTED OATP1B residual uptake capacity, CL_OATP/CL_NTCP (-)
INBA    : 5.16655 : bile acid synthesis + ileal return (umol/L/day)

// ------------------------------------------------------- peg-interferon alfa/lambda
KAIFN   : 1.44    : SC absorption rate constant (1/day)
FIFN    : 0.80    : SC bioavailability (-)
VCIFN   : 10.0    : central volume (L)
CLIFN   : 1.40    : clearance (L/day)
EC50IFN : 6.0     : ISG induction EC50 (ng/mL)
HIFN    : 1.5     : ISG induction Hill coefficient (-)
KINISG  : 1.0     : ISG production rate constant (1/day)
KOUTISG : 1.0     : ISG loss rate constant (1/day)
KASOCS  : 0.5     : SOCS1 production rate constant (1/day)
KDSOCS  : 0.5     : SOCS1 loss rate constant (1/day)
KSSOCS  : 1.0     : SOCS1 level that halves ISG production (-)
STIMEND : 0.25    : maximal endogenous MDA5-driven ISG drive (-)
KENDO   : 0.05    : infected fraction giving half the endogenous drive (-)
ISG0    : 0.145497 : endogenous ISG tone at baseline; ZERO of all drug effect
EMAXREP : 0.75    : max block of intracellular HDV replication (-)
EC50REP : 0.25    : replication block EC50 on the ISG increment (-)
EMAXPRD : 0.90    : max block of virion PRODUCTION/export; makes the first phase
EC50PRD : 0.25    : production block EC50 on the ISG increment (-)
DMAXIFN : 0.014   : max ADDED infected-cell loss rate (1/day)
EC50DIF : 0.35    : added-loss EC50 on the ISG increment (-)
EMAXHBS : 0.50    : max reduction of HBsAg secretion (-)
EC50HBS : 0.40    : HBsAg reduction EC50 (-)
EADAR   : 1.20    : max fold-increase of amber/W editing (ADAR1 is IFN-inducible)
EC50ADR : 0.40    : ADAR1 induction EC50 (-)
ECURE   : 0.50    : ISG boost of non-cytolytic intracellular clearance (-)
BIFNE   : 1.50    : ISG boost of effector activity (-)
CIFNX   : 3.00    : ISG-driven reversal of exhaustion (-)
LAMBDA  : 0.0     : 1 = peg-IFN lambda (no myelosuppression, hepatic AE instead)
ALTTOXL : 5.0     : peg-IFN lambda hepatic ALT toxicity gain (U/L per ISG unit)

// ---------------------------------------------------- lonafarnib + ritonavir
KALNF   : 6.0     : oral absorption rate constant (1/day)
FLNF    : 0.30    : oral bioavailability (-)
VCLNF   : 100.0   : central volume (L)
VPLNF   : 200.0   : peripheral volume (L)
QLNF    : 30.0    : intercompartmental clearance (L/day)
CLLNF   : 480.0   : unboosted CYP3A4 clearance (L/day)
KIRTV   : 0.40    : ritonavir CYP3A4 inhibition constant (uM)
IC50FT  : 0.10    : farnesyltransferase inhibition IC50 (uM)
VCRTV   : 100.0   : ritonavir central volume (L)
CLRTV   : 130.0   : ritonavir clearance (L/day)

// ------------------------------------------------------- HBsAg-directed siRNA
CLSIRN  : 40.0    : siRNA plasma clearance (L/day)
VCSIRN  : 5.0     : siRNA central volume (L)
KINSIR  : 0.20    : transfer into the hepatic effect site (1/day)
KOUTSIR : 0.030   : hepatic effect-site loss (1/day), t1/2 ~ 23 d
EMAXSIR : 0.85    : max HBsAg transcript knockdown (-)
EC50SIR : 0.50    : knockdown EC50 (effect-site a.u.)

// ------------------------------- hepatocyte pools and the maintenance fluxes
DHEP    : 0.0040  : background hepatocyte death (1/day)
DHBV    : 0.0020  : extra death of HBsAg+ HDV-negative cells (1/day)
DINN    : 0.0040  : innate (non-CTL) loss of infected cells (1/day)
DIMM    : 0.02636 : FITTED baseline immune killing of infected cells (1/day)
KCURE   : 0.0020  : non-cytolytic intracellular HDV clearance (1/day)
PHINT   : 0.50    : NTCP-DEPENDENT fraction of cell-to-cell spread (-)
LOCREN  : 0.999   : regeneration locality within HBsAg+ / HBsAg- classes (-)
BETAD   : 1.73924e-08 : HDV de-novo infection rate constant (mL/IU/day)
BETAB   : 2.0e-09 : HBV de-novo infection rate constant (mL/IU/day)
KCC     : 0.0137499 : cell-to-cell spread rate constant (1/day)
GAMMAX  : 0.15    : killing gain per unit of REVERSED exhaustion (-)
NHEP    : 2.0e11  : hepatocytes in an adult liver (-)

// --------------------------------------- intracellular replication / assembly
ALPHAR  : 2.816   : rolling-circle replication rate constant (1/day)
RMAX    : 6000.0  : intracellular HDV RNA carrying capacity (copies/cell)
MUR     : 0.30    : intracellular HDV RNA turnover (1/day)
KSS     : 100.0   : S-HDAg support of replication (a.u.)
KIL     : 3000.0  : L-HDAg trans-dominant repression constant (a.u.)
KTL     : 1.00    : HDAg translation scale (1/day)
FEDIT0  : 0.25    : baseline amber/W editing fraction (-)
MUS     : 1.50    : S-HDAg turnover (1/day)
KPREN   : 2.00    : farnesylation of L-HDAg (1/day)
MUL     : 0.50    : unprenylated L-HDAg turnover (1/day)
MULP    : 1.50    : prenylated L-HDAg turnover (1/day)
KEX     : 2.57142 : per-cell export rate constant (1/day)
KPLP    : 400.0   : prenyl-L-HDAg half-saturation for envelopment (a.u.)
KSENV   : 3034.0  : serum HBsAg half-saturation for envelopment (IU/mL)
CONVSER : 644.168 : copies/cell/day -> IU/mL/day conversion (-)
CD      : 0.55    : serum HDV RNA clearance (1/day), t1/2 ~ 30 h

// ------------------------------------------------ HBV side = envelope factory
CCCMAX  : 20.0    : cccDNA carrying capacity (copies/cell)
RCCC    : 0.050   : cccDNA replenishment rate constant (1/day)
MUCCC   : 0.025   : cccDNA loss rate constant (1/day)
KHBS    : 550.624 : HBsAg production rate constant (IU/mL/day per cccDNA)
CS      : 0.35    : serum HBsAg clearance (1/day), t1/2 ~ 2 d
KHBV    : 26.6667 : HBV DNA production rate constant (IU/mL/day per cccDNA)
CB      : 0.60    : serum HBV DNA clearance (1/day)
KSUP    : 0.10    : HDV suppression of HBV replication (-)
EPSNUC  : 0.999   : NUC block of HBV DNA production (-)
NUC     : 1.0     : 1 = on a nucleos(t)ide analogue, 0 = off

// ---------------------------------------------------------------- immunity
AE      : 0.060   : effector production rate constant (a.u./day)
KEHALF  : 0.020   : infected fraction giving half effector drive (-)
DE      : 0.050   : effector loss rate constant (1/day)
AX      : 0.030   : exhaustion accumulation rate constant (1/day)
KXHALF  : 0.0005  : infected fraction giving half exhaustion drive (-)
DX      : 0.0199005 : exhaustion reversal rate constant (1/day)

// -------------------------------------------- injury, fibrosis, organ, outcome
KAPPAF  : 1.05144 : injury per unit ENTRY flux -- the ALT/RNA decoupling term
ALTBASE : 42.0    : FITTED median NON-HDV component of ALT (U/L)
KELALT  : 0.35    : ALT elimination rate constant (1/day), t1/2 ~ 47 h
KALT    : 4703.56 : ALT gain per unit killing flux (U/L per fraction/day)
ULNALT  : 40.0    : ALT upper limit of normal used by the trials (U/L)
AH      : 0.00595041 : stellate cell activation rate constant (1/day)
DH      : 0.020   : stellate cell deactivation rate constant (1/day)
ALTREF  : 40.0    : ALT reference for stellate activation (U/L)
KF      : 0.001795 : fibrogenesis rate constant (Ishak/day)
KREV    : 0.0002  : fibrolysis rate constant (Ishak/day)
KHREV   : 0.15    : activated HSC level that halves fibrolysis (-)
PLT0    : 250.0   : platelet count with no fibrosis (10^9/L)
GPLT    : 0.62    : portal-hypertension platelet effect at Ishak 6 (-)
MYEPLT  : 0.42    : peg-IFN alfa myelosuppression on platelets (-)
KPLT    : 0.15    : platelet turnover rate constant (1/day)
NEU0    : 4.0     : neutrophil count off treatment (10^9/L)
MYENEU  : 0.55    : peg-IFN alfa myelosuppression on neutrophils (-)
KNEU    : 0.20    : neutrophil turnover rate constant (1/day)
H0HCC   : 1.36986e-05 : HCC hazard at Ishak 2 (1/day)
BFHCC   : 0.4302  : HCC log-hazard per Ishak stage (-)
BVHCC   : 0.15    : HCC log-hazard per log10 HDV RNA (-)

$CMT @annotated
// ---- drug PK / PD (13) ----
ASC     : bulevirtide SC depot (nmol)
ACEN    : bulevirtide central (nmol)
ALIV    : bulevirtide internalised hepatic pool (nmol)
ISC     : peg-IFN SC depot (ug)
ICEN    : peg-IFN central (ug)
ISG     : interferon-stimulated-gene activity (-)
SOCS    : SOCS1/USP18 desensitisation signal (-)
LGUT    : lonafarnib gut (umol)
LCEN    : lonafarnib central (umol)
LPER    : lonafarnib peripheral (umol)
RCEN    : ritonavir central (umol)
QCEN    : siRNA central (nmol)
QEFF    : siRNA hepatic effect site (a.u.)
// ---- hepatocyte pools (3) ----
T       : HBsAg-negative hepatocytes (fraction)
IB      : HBsAg+ HDV-negative hepatocytes, the ENVELOPE DONOR pool (fraction)
ID      : HBsAg+ HDV+ hepatocytes, productively co-infected (fraction)
// ---- intracellular HDV (4) ----
RG      : intracellular HDV genomic RNA (copies/cell)
SAG     : S-HDAg (a.u./cell)
LAG     : L-HDAg, unprenylated (a.u./cell)
LAGP    : L-HDAg, prenylated and assembly-competent (a.u./cell)
// ---- serum virology (4) ----
VD      : serum HDV RNA (IU/mL)
CCC     : HBV cccDNA per HBsAg+ cell (copies/cell)
SSER    : serum HBsAg (IU/mL)
VB      : serum HBV DNA (IU/mL)
// ---- immunity (2) ----
EC      : HDV-specific CD8 effector activity (a.u.)
EXH     : exhaustion level (-)
// ---- injury / organ / clinical (7) ----
ALT     : serum ALT (U/L)
HSCA    : activated hepatic stellate cells (-)
FIB     : fibrosis stage, Ishak (-)
PLT     : platelets (10^9/L)
NEU     : neutrophils (10^9/L)
TBA     : total serum bile acids (umol/L)
CH      : cumulative HCC hazard (-)

$MAIN
// ---- baseline hepatocyte pools and the observed chronic steady state -------
T_0    = 0.55;
IB_0   = 0.35;
ID_0   = 0.10;
RG_0   = 3000.0;
SAG_0  = 1500.0;
LAG_0  = 300.0;
LAGP_0 = 400.0;
VD_0   = 316228.0;      // 10^5.5 IU/mL
CCC_0  = 10.0;
SSER_0 = 7079.46;       // 10^3.85 IU/mL
VB_0   = 100.0;
EC_0   = 1.0;
EXH_0  = 0.60;
ALT_0  = 110.0;
HSCA_0 = 0.45;
FIB_0  = 2.0;
PLT_0  = PLT0 * (1.0 - GPLT * 2.0 / 6.0);
NEU_0  = NEU0;
TBA_0  = TBA0;
ISG_0  = ISG0;
SOCS_0 = ISG0;

$ODE
// ===========================================================================
//  1.  BULEVIRTIDE PK AND NTCP OCCUPANCY
//      Clearance is split into a linear route and a SATURABLE route (NTCP-
//      mediated internalisation = target-mediated clearance).  The saturable
//      route is what makes exposure MORE than dose proportional -- Css goes
//      1.72 -> 14.47 nM for a 5-fold dose step -- and that supra-proportional
//      exposure is what the bile-acid inversion needs in order to place the
//      two doses at different points on the same occupancy curve.
// ===========================================================================
double Cblv = ACEN / VCBLV;                                  // nM
double upt  = VMBLV * Cblv / (KMBLV + Cblv);
dxdt_ASC  = -KABLV * ASC;
dxdt_ACEN = FBLV * KABLV * ASC - CLBLV * Cblv - upt;
dxdt_ALIV = upt - KOUTLIV * ALIV;

double OCC = Cblv / (Cblv + KDNTCP);
if (OCCFIX >= 0.0) OCC = OCCFIX;      // hypothetical perfect entry blocker

// ===========================================================================
//  2.  PEG-INTERFERON PK -> ISG, WITH SOCS1 TACHYPHYLAXIS
//      Chronic hepatitis D already runs a strong endogenous ISG signature and
//      the virus persists anyway, so ALL interferon effects below are driven by
//      the INCREMENT above that endogenous tone.  If they were written on
//      absolute ISG, the untreated patient would already be partly treated and
//      could not be a steady state.
// ===========================================================================
double Cifn = ICEN / VCIFN;                                  // ng/mL
dxdt_ISC  = -KAIFN * ISC;
dxdt_ICEN = FIFN * KAIFN * ISC - CLIFN * Cifn;

double stimD = (Cifn > 0.0)
  ? pow(Cifn, HIFN) / (pow(Cifn, HIFN) + pow(EC50IFN, HIFN)) : 0.0;
double stimE = STIMEND * ID / (ID + KENDO);
double stim  = stimD + stimE; if (stim > 1.0) stim = 1.0;
dxdt_ISG  = KINISG * stim / (1.0 + SOCS / KSSOCS) - KOUTISG * ISG;
dxdt_SOCS = KASOCS * ISG - KDSOCS * SOCS;

double dis     = ISG - ISG0; if (dis < 0.0) dis = 0.0;   // the drug increment
double epsRep  = EMAXREP * dis / (EC50REP + dis);
double epsProd = EMAXPRD * dis / (EC50PRD + dis);        // -> the first phase
double dltIfn  = DMAXIFN * dis / (EC50DIF + dis);        // -> durability
double epsHbs  = EMAXHBS * dis / (EC50HBS + dis);
double fEdit   = FEDIT0 * (1.0 + EADAR * dis / (EC50ADR + dis));
if (fEdit > 0.95) fEdit = 0.95;
double kCure   = KCURE * (1.0 + ECURE * dis);
double hem     = (LAMBDA > 0.5) ? 0.0 : dis;            // myelosuppression

// ===========================================================================
//  3.  LONAFARNIB + RITONAVIR -> FARNESYLTRANSFERASE INHIBITION
// ===========================================================================
double Clnf = LCEN / VCLNF;                                  // uM
double Crtv = RCEN / VCRTV;                                  // uM
double CLeff = CLLNF / (1.0 + Crtv / KIRTV);                 // CYP3A4 boost
dxdt_LGUT = -KALNF * LGUT;
dxdt_LCEN = FLNF * KALNF * LGUT - CLeff * Clnf
            - QLNF * (Clnf - LPER / VPLNF);
dxdt_LPER = QLNF * (Clnf - LPER / VPLNF);
dxdt_RCEN = -CLRTV * Crtv;
double IFT = Clnf / (Clnf + IC50FT);

// ===========================================================================
//  4.  HBsAg-DIRECTED siRNA
// ===========================================================================
double Cq = QCEN / VCSIRN;
dxdt_QCEN = -CLSIRN * Cq;
dxdt_QEFF = KINSIR * Cq - KOUTSIR * QEFF;
double epsS = EMAXSIR * QEFF / (EC50SIR + QEFF);

// ===========================================================================
//  5.  INTRACELLULAR HDV -- replication, editing, prenylation, export
//      Export is the DOMINANT drain on the intracellular pool.  That is what
//      makes an assembly block RAISE intracellular HDV RNA instead of merely
//      lowering serum RNA: the lonafarnib paradox in result (5) above.
// ===========================================================================
double gsupp = SAG / (SAG + KSS);
double grepr = 1.0 / (1.0 + LAG / KIL);
double gLP   = LAGP / (LAGP + KPLP);
double gS    = SSER / (SSER + KSENV);
double grow  = 1.0 - RG / RMAX; if (grow < 0.0) grow = 0.0;
double expo  = KEX * RG * gLP * gS * (1.0 - IFT) * (1.0 - epsProd);

dxdt_RG   = ALPHAR * RG * grow * gsupp * grepr * (1.0 - epsRep)
            - MUR * RG - expo;
dxdt_SAG  = KTL * RG * (1.0 - fEdit) - MUS * SAG;
dxdt_LAG  = KTL * RG * fEdit - KPREN * (1.0 - IFT) * LAG - MUL * LAG;
dxdt_LAGP = KPREN * (1.0 - IFT) * LAG - MULP * LAGP;

// ===========================================================================
//  6.  HEPATOCYTE POOLS -- the three maintenance fluxes and the FLOOR
//
//      (E) F_entry   NTCP-dependent, fully blockable
//      (C) F_cc      only the PHINT fraction is NTCP-dependent
//      (D) renId     an Id cell that divides gives Id daughters: UNBLOCKABLE
//
//      Regeneration is LOCAL.  A slot vacated by a dying hepatocyte is
//      refilled by a neighbour, and the neighbours of an HBsAg+ cell are
//      mostly HBsAg+ too.  That is what makes (a) the HBsAg+ pool persist for
//      decades and (b) HDV+ cells get DILUTED by dividing HDV-negative
//      HBsAg+ neighbours -- the mechanism the HBV field calls cccDNA dilution
//      by hepatocyte turnover.  Without it the whole HBsAg+ compartment drains
//      within a year, which patients plainly do not do.
// ===========================================================================
double N  = T + IB + ID;   if (N  < 1e-9)  N  = 1e-9;
double Hb = IB + ID;       if (Hb < 1e-12) Hb = 1e-12;

double Fentry = BETAD * VD * IB * (1.0 - OCC);
double Fcc    = KCC * ID * IB * (1.0 - PHINT * OCC);
double Fhbv   = BETAB * VB * T * (1.0 - OCC);

// exhaustion REVERSED below baseline is what keeps killing capacity up after
// the drug stops -- see result (6).
double revx = (0.60 - EXH) / 0.60; if (revx < 0.0) revx = 0.0;
double fx   = 1.0 + GAMMAX * revx;

double dthT  = DHEP;
double dthIb = DHEP + DHBV;
double dthId = DHEP + DINN + DIMM * fx + dltIfn;

double Hd    = dthIb * IB + dthId * ID;      // deaths in the HBsAg+ class
double Dt    = dthT * T;                     // deaths in the HBsAg- class
double poolR = (1.0 - LOCREN) * (Dt + Hd);
double renT  = LOCREN * Dt + poolR * T  / N;
double renIb = LOCREN * Hd * IB / Hb + poolR * IB / N;
double renId = LOCREN * Hd * ID / Hb + poolR * ID / N;   // (D) UNBLOCKABLE

dxdt_T  = renT  - dthT  * T  - Fhbv;
dxdt_IB = renIb - dthIb * IB + Fhbv - Fentry - Fcc + kCure * ID;
dxdt_ID = renId - dthId * ID         + Fentry + Fcc - kCure * ID;

// ===========================================================================
//  7.  SERUM VIROLOGY
// ===========================================================================
dxdt_VD   = CONVSER * expo * ID - CD * VD;
double cgrow = 1.0 - CCC / CCCMAX; if (cgrow < 0.0) cgrow = 0.0;
dxdt_CCC  = RCCC * CCC * cgrow - MUCCC * CCC * (1.0 + 0.6 * dis);
dxdt_SSER = KHBS * CCC * Hb * (1.0 - epsS) * (1.0 - epsHbs) - CS * SSER;
double nucf = (NUC > 0.5) ? (1.0 - EPSNUC) : 1.0;
dxdt_VB   = KHBV * CCC * Hb * nucf / (1.0 + ID / KSUP) - CB * VB;

// ===========================================================================
//  8.  IMMUNITY.  EC and EXH are reported biomarkers; only the REVERSAL of
//      exhaustion feeds back on killing (see fx above).  Exhaustion is made
//      almost insensitive to the size of the infected pool (KXHALF is tiny)
//      and sensitive to ISG instead: that is the discriminating assumption --
//      interferon reverses exhaustion, an entry inhibitor does not, even
//      though both lower antigen.
// ===========================================================================
dxdt_EC  = AE * (ID / (ID + KEHALF)) * (1.0 + BIFNE * dis) - DE * EC;
dxdt_EXH = AX * (ID / (ID + KXHALF)) * (1.0 - EXH)
           - DX * EXH * (1.0 + CIFNX * dis);

// ===========================================================================
//  9.  INJURY -> FIBROSIS -> OUTCOME
//      NOTE the second injury term.  ALT is driven not only by the SIZE of the
//      infected pool but by the INFLOW of newly infected hepatocytes, a term
//      that collapses within days of blocking entry while the pool itself
//      decays over months.  That is the model's account of ALT normalisation
//      in virologic non-responders, and it is a HYPOTHESIS: the falsifier is
//      in hdv_model_report.txt A13.
// ===========================================================================
double kill = (DINN + DIMM * fx + dltIfn) * ID + KAPPAF * Fentry;
double altTox = (LAMBDA > 0.5) ? ALTTOXL * dis : 0.0;
dxdt_ALT = KALT * kill + KELALT * ALTBASE + altTox - KELALT * ALT;

dxdt_HSCA = AH * (ALT / ALTREF) * (1.0 - HSCA) - DH * HSCA;
double fgrow = 1.0 - FIB / 6.0; if (fgrow < 0.0) fgrow = 0.0;
dxdt_FIB  = KF * HSCA * fgrow - KREV * FIB / (1.0 + HSCA / KHREV);

double pltT = PLT0 * (1.0 - GPLT * FIB / 6.0) * (1.0 - MYEPLT * hem);
dxdt_PLT = KPLT * (pltT - PLT);
double neuT = NEU0 * (1.0 - MYENEU * hem);
dxdt_NEU = KNEU * (neuT - NEU);

// bile acids: the same NTCP, read as an occupancy biomarker
double CLupt = CLNTCPB * (1.0 - OCC) + CLNTCPB * ROATP;
dxdt_TBA = INBA - CLupt * TBA;

double lgVD = log10(VD > 1e-12 ? VD : 1e-12);
dxdt_CH = H0HCC * exp(BFHCC * (FIB - 2.0) + BVHCC * (lgVD - 5.5));

$TABLE
double Cblv_o = ACEN / VCBLV;
double OCC_o  = (OCCFIX >= 0.0) ? OCCFIX : Cblv_o / (Cblv_o + KDNTCP);
capture CBLV   = Cblv_o;
capture OCCUP  = OCC_o;                       // DERIVED NTCP occupancy
capture FREENT = 1.0 - OCC_o;                 // residual free NTCP
capture CIFN   = ICEN / VCIFN;
capture CLNF   = LCEN / VCLNF;
capture CRTV   = RCEN / VCRTV;
capture IFTC   = (LCEN / VCLNF) / ((LCEN / VCLNF) + IC50FT);
capture LGVD   = log10(VD > 1e-12 ? VD : 1e-12);
capture LGSAG  = log10(SSER > 1e-12 ? SSER : 1e-12);
capture LGVB   = log10(VB > 1e-12 ? VB : 1e-12);
capture DLGVD  = log10(VD > 1e-12 ? VD : 1e-12) - 5.5;
// NOTE: DLGVD is referenced to the DEFAULT baseline of 5.5 log10 IU/mL.  If you
// override VD_0 (e.g. from the Shiny sidebar), compute LGVD minus your own
// baseline instead -- the app does exactly that.
capture TBAFLD = TBA / TBA0;                  // fold rise = occupancy read-out
capture ALTN   = (ALT <= ULNALT) ? 1.0 : 0.0; // ALT normalised
capture VR     = ((VD < 6.0) || (LGVD <= 3.5)) ? 1.0 : 0.0;  // >=2 log or <LOD
capture COMB   = (((VD < 6.0) || (LGVD <= 3.5)) && (ALT <= ULNALT)) ? 1.0 : 0.0;
capture HDAGPC = 100.0 * ID;                  // % of hepatocytes HDAg-positive
capture HCCINC = 100.0 * (1.0 - exp(-CH));    // cumulative HCC incidence (%)
capture RGFOLD = RG / 3000.0;                 // intracellular RNA vs baseline

## =============================================================================
##  R-SIDE HELPERS -- 13 THERAPEUTIC SCENARIOS
##  Everything below is ordinary R and is ignored by mread().
## =============================================================================
$ENV

## --- dosing builders --------------------------------------------------------
## bulevirtide: mg SC once daily -> nmol into ASC
blv_dose <- function(mg, start = 0, dur = 336) {
  if (mg <= 0) return(NULL)
  mrgsolve::ev(time = start, amt = mg * 1e-3 / 5398.9 * 1e9,
               cmt = "ASC", ii = 1, addl = dur - 1)
}
## peg-IFN alfa-2a or lambda-1a: ug SC once weekly
ifn_dose <- function(ug, start = 0, dur = 336) {
  if (ug <= 0) return(NULL)
  mrgsolve::ev(time = start, amt = ug, cmt = "ISC",
               ii = 7, addl = floor(dur / 7) - 1)
}
## lonafarnib mg PO BID (+ ritonavir mg PO BID as a booster)
lnf_dose <- function(mg, rtv = 100, start = 0, dur = 336) {
  if (mg <= 0) return(NULL)
  e <- mrgsolve::ev(time = start, amt = mg * 1e-3 / 638.6 * 1e6,
                    cmt = "LGUT", ii = 0.5, addl = dur * 2 - 1)
  if (rtv > 0) {
    e <- e + mrgsolve::ev(time = start, amt = 0.7 * rtv * 1e-3 / 720.9 * 1e6,
                          cmt = "RCEN", ii = 0.5, addl = dur * 2 - 1)
  }
  e
}
## HBsAg-directed siRNA mg SC every 4 weeks
sirna_dose <- function(mg, start = 0, dur = 336) {
  if (mg <= 0) return(NULL)
  mrgsolve::ev(time = start, amt = mg * 1e-3 / 16000 * 1e9, cmt = "QCEN",
               ii = 28, addl = floor(dur / 28) - 1)
}
combine_ev <- function(...) {
  es <- Filter(Negate(is.null), list(...))
  if (!length(es)) return(mrgsolve::ev(time = 0, amt = 0, cmt = "ASC"))
  Reduce(`+`, es)
}

## --- the 13 scenarios -------------------------------------------------------
## Each entry: label, event object, total simulated horizon (days), parameter
## overrides, and what the scenario is FOR.
hdv_scenarios <- function() {
  list(
    list(id = "S01", lab = "Natural history (NUC only, no anti-HDV therapy)",
         ev = NULL, tend = 5 * 365, par = list(),
         why = "the control arm; must stay flat on HDV RNA and ALT while fibrosis progresses"),

    list(id = "S02", lab = "Bulevirtide 2 mg SC qd, 48 weeks",
         ev = blv_dose(2, 0, 336), tend = 336, par = list(),
         why = "the approved regimen; MYR301 week-48 anchors live here"),

    list(id = "S03", lab = "Bulevirtide 2 mg SC qd, 96 weeks",
         ev = blv_dose(2, 0, 672), tend = 672, par = list(),
         why = "does the response keep deepening?  MYR301 says 45% -> 55%"),

    list(id = "S04", lab = "Bulevirtide 10 mg SC qd, 48 weeks",
         ev = blv_dose(10, 0, 336), tend = 336, par = list(),
         why = "PREDICTION: 6x less residual entry than 2 mg -- does it matter?"),

    list(id = "S05", lab = "Bulevirtide 2 mg, 48 weeks, then STOP",
         ev = blv_dose(2, 0, 336), tend = 504, par = list(),
         why = "an entry inhibitor is suppressive by construction: expect relapse"),

    list(id = "S06", lab = "Peg-IFN alfa-2a 180 ug qw, 48 weeks, then STOP",
         ev = ifn_dose(180, 0, 336), tend = 504, par = list(),
         why = "HIDIT-1 comparator; the only monotherapy with off-treatment durability"),

    list(id = "S07", lab = "Peg-IFN alfa-2a 180 ug qw, 96 weeks, then STOP",
         ev = ifn_dose(180, 0, 672), tend = 840, par = list(),
         why = "HIDIT-2: does doubling the duration buy durability?"),

    list(id = "S08", lab = "Bulevirtide 2 mg + peg-IFN alfa 180 ug, 48 wk, then STOP",
         ev = combine_ev(blv_dose(2, 0, 336), ifn_dose(180, 0, 336)),
         tend = 504, par = list(),
         why = "MYR204: inflow-cut x outflow-raise, the structural synergy"),

    list(id = "S09", lab = "Bulevirtide 10 mg + peg-IFN alfa 180 ug, 48 wk, then STOP",
         ev = combine_ev(blv_dose(10, 0, 336), ifn_dose(180, 0, 336)),
         tend = 504, par = list(),
         why = "the best combination tested clinically"),

    list(id = "S10", lab = "Lonafarnib 50 mg BID + ritonavir 100 mg BID, 48 wk, then STOP",
         ev = lnf_dose(50, 100, 0, 336), tend = 504, par = list(),
         why = "the assembly block: watch RGFOLD rise while LGVD falls, then rebound"),

    list(id = "S11", lab = "Lonafarnib/ritonavir + peg-IFN alfa, 48 weeks, then STOP",
         ev = combine_ev(lnf_dose(50, 100, 0, 336), ifn_dose(180, 0, 336)),
         tend = 504, par = list(),
         why = "D-LIVR triple therapy"),

    list(id = "S12", lab = "Peg-IFN lambda-1a 180 ug qw, 48 weeks, then STOP",
         ev = ifn_dose(180, 0, 336), tend = 504, par = list(LAMBDA = 1),
         why = "LIMT-1: hepatocyte-restricted receptor, no myelosuppression, hepatic AE"),

    list(id = "S13", lab = "Bulevirtide 2 mg + HBsAg siRNA 200 mg q4w, 48 weeks",
         ev = combine_ev(blv_dose(2, 0, 336), sirna_dose(200, 0, 336)),
         tend = 336, par = list(),
         why = "starve the envelope: the only non-IFN mechanism that moves the FLOOR"),

    list(id = "S14", lab = "HYPOTHETICAL perfect entry inhibitor (100% occupancy)",
         ev = NULL, tend = 336, par = list(OCCFIX = 1.0),
         why = "the CEILING of the entire entry-inhibitor class")
  )
}

## --- run one scenario -------------------------------------------------------
run_scenario <- function(mod, sc, delta = 1) {
  m <- mod
  if (length(sc$par)) m <- mrgsolve::param(m, sc$par)
  m <- mrgsolve::update(m, end = sc$tend, delta = delta)
  out <- if (is.null(sc$ev)) mrgsolve::mrgsim(m) else mrgsolve::mrgsim(m, events = sc$ev)
  d <- as.data.frame(out)
  d$scenario <- sc$id
  d$label <- sc$lab
  d
}

## --- run everything and summarise -------------------------------------------
scenario_all <- function(mod, delta = 1) {
  scs <- hdv_scenarios()
  do.call(rbind, lapply(scs, function(s) run_scenario(mod, s, delta)))
}

## week-48 / end-of-horizon summary table
scenario_summary <- function(mod) {
  scs <- hdv_scenarios()
  rows <- lapply(scs, function(s) {
    d <- run_scenario(mod, s)
    pick <- function(tt) d[which.min(abs(d$time - tt)), ]
    w48 <- pick(min(336, s$tend))
    eos <- d[nrow(d), ]
    data.frame(
      id = s$id, label = s$lab,
      wk48_dlog10_HDV = round(w48$DLGVD, 2),
      wk48_ALT = round(w48$ALT, 1),
      wk48_ALT_normal = w48$ALTN == 1,
      wk48_HDAg_pos_pct = round(w48$HDAGPC, 3),
      wk48_intracell_RNA_fold = round(w48$RGFOLD, 2),
      wk48_NTCP_occupancy = round(w48$OCCUP, 3),
      wk48_bile_acid_fold = round(w48$TBAFLD, 1),
      end_day = round(eos$time),
      end_dlog10_HDV = round(eos$DLGVD, 2),
      end_ALT = round(eos$ALT, 1),
      end_Ishak = round(eos$FIB, 2),
      end_platelets = round(eos$PLT),
      end_HCC_pct = round(eos$HCCINC, 2),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

## --- the two analyses that carry the argument -------------------------------
## A. bile acids invert to NTCP occupancy, and benefit saturates while cost
##    does not.  This is the dose-finding instrument of section A12.
blv_dose_response <- function(mod, doses = c(0.5, 1, 2, 5, 10, 20)) {
  do.call(rbind, lapply(doses, function(mg) {
    d <- as.data.frame(mrgsolve::mrgsim(
      mrgsolve::update(mod, end = 336, delta = 7),
      events = blv_dose(mg, 0, 336)))
    e <- d[nrow(d), ]
    data.frame(dose_mg = mg,
               Css_nM = round(mean(tail(d$CBLV, 48)), 3),
               NTCP_occupancy = round(e$OCCUP, 4),
               free_NTCP_pct = round(100 * e$FREENT, 1),
               bile_acid_fold = round(e$TBAFLD, 1),
               dlog10_HDV_wk48 = round(e$DLGVD, 2),
               ALT_wk48 = round(e$ALT, 1))
  }))
}

## B. what sets the FLOOR?  Perturb each parameter by +/-30% at FIXED 2 mg
##    dosing and see what moves the week-48 response.  If the model is right,
##    host cell biology outranks target affinity.
floor_sensitivity <- function(mod, pars = c("DHEP", "KCC", "PHINT", "KCURE",
                                            "DIMM", "KDNTCP"), f = 0.3) {
  base <- as.data.frame(mrgsolve::mrgsim(
    mrgsolve::update(mod, end = 336, delta = 14), events = blv_dose(2, 0, 336)))
  b <- base[nrow(base), "DLGVD"]
  do.call(rbind, lapply(pars, function(pn) {
    v <- unlist(mrgsolve::param(mod))[[pn]]
    got <- sapply(c(1 - f, 1 + f), function(k) {
      m <- mrgsolve::param(mod, setNames(list(v * k), pn))
      d <- as.data.frame(mrgsolve::mrgsim(
        mrgsolve::update(m, end = 336, delta = 14), events = blv_dose(2, 0, 336)))
      d[nrow(d), "DLGVD"]
    })
    data.frame(parameter = pn, base_dlog10 = round(b, 2),
               minus30 = round(got[1], 2), plus30 = round(got[2], 2),
               abs_effect = round(mean(abs(got - b)), 3))
  }))
}

## --- quick look -------------------------------------------------------------
## mod <- mread("hdv_mrgsolve_model.R")
## print(scenario_summary(mod))
## print(blv_dose_response(mod))
## print(floor_sensitivity(mod))
##
## library(ggplot2)
## out <- scenario_all(mod)
## ggplot(out, aes(time / 7, DLGVD, colour = scenario)) + geom_line(linewidth = .7) +
##   labs(x = "week", y = expression(Delta*log[10]*" HDV RNA (IU/mL)"),
##        title = "Chronic hepatitis D: entry vs envelopment vs infected-cell loss") +
##   theme_minimal()
##
## Expected shape, and the whole point of the figure: bulevirtide falls with NO
## first phase, peg-IFN and lonafarnib both drop steeply in week 1 and then
## flatten onto the same cell-loss-limited slope, and only the peg-IFN-containing
## arms hold their response after the drug stops.
