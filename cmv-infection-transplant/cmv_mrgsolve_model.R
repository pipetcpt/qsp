## ===========================================================================
##  cmv_mrgsolve_model.R
##  Cytomegalovirus infection and disease in transplant recipients
##  이식 후 거대세포바이러스(CMV) 감염 — 48-ODE mrgsolve QSP model
##
##  This file MIRRORS cmv_python_reference.py equation for equation.  Every
##  number in the calibration notes below was produced by RUNNING that Python
##  reference (output in cmv_reference_output.txt), because no R/mrgsolve
##  toolchain was available in the environment where this model was written.
##  That is stated here rather than implied: the ODE system has been integrated
##  and verified, but not by this file.
##
##  ORGANISING THESIS -- ONE NUMBER, THEN A SECOND
##  ---------------------------------------------------------------------------
##  Write the systemic virus as a target-cell-limited process.  At quasi-steady
##  state in the virion pool the net exponent of the infected-cell compartment is
##
##      r = KPROD*(1-e_pol)*(1-e_pack) - DELI - KE8*E8eff - KENK*NKA
##
##  Two measured clinical numbers fix KPROD:
##      untreated DNAemia doubling time  ~1.2 d  ->  r0   = 0.5776 /d
##      on-treatment decline half-life   ~2.4 d  ->  DELI = 0.2888 /d
##      KPROD = r0 + DELI                       =  0.8664 /d
##
##  THRESHOLD 1 (the drug threshold)
##      e* = r0/KPROD = 0.6667
##      Potency, resistance, renal dose banding and drug-drug interaction are
##      all the same question asked of this one number.  Verified values:
##          valGCV 900 od      e = 0.890      margin +0.223   CONTROL
##          valGCV 900 BID     e = 0.970      margin +0.303   CONTROL
##          valGCV 450 q2d     e = 0.770      margin +0.104   CONTROL  (CrCl 30)
##          valGCV 450 2x/wk   e = 0.657      margin -0.010   BREAKTHROUGH
##          letermovir 480 od  e = 0.969      margin +0.303   CONTROL
##          letermovir 240 od  e = 0.940      margin +0.274   CONTROL  (+CsA)
##          maribavir 400 BID  e = 0.792      margin +0.125   CONTROL
##          foscarnet 90 q12h  e = 0.719      margin +0.053   CONTROL
##      and on the resistant strains:
##          valGCV 900 od  / UL97 mut   e = 0.112   (12.6% of the effect left)
##          valGCV 900 BID / UL97 mut   e = 0.335   still BELOW e*
##          letermovir     / UL56 C325Y e = 0.010   ( 1.1% of the effect left)
##      -> the GCV-TP needed to clear e* on the UL97 mutant is 265 uM-eq, which
##         implies a plasma Cavg of 20.4 uM = 2.0x the licensed treatment dose.
##         Dose escalation cannot rescue it; the drug must be changed.
##
##  THRESHOLD 2 (the immune threshold)
##      E8* = r0/KE8 = 4.81 CMV-specific CD8 per uL -- control with NO drug.
##      This is the endpoint of the illness and no drug moves it.
##
##  THE COUPLING (why the two fight each other)
##      d(E8)/dt is proportional to ANTIGEN = the infected-cell pool the drug
##      just removed.  Clearing threshold 1 therefore delays reaching threshold
##      2.  Verified: at the moment 200-day valganciclovir prophylaxis stops,
##      E8eff = 0.094 /uL against E8* = 4.81 -- a 51-fold shortfall -- and
##      DNAemia recrosses 1000 IU/mL 19 days later.  Late-onset CMV is not bolted
##      on; it is what this coupling does when the drug stops.
##
##  THE THIRD AXIS (monitoring interval is a dose)
##      Pre-emptive therapy does not act on e.  It sets the HEIGHT at which the
##      exponential is caught: fold-rise between draws = 2^(dt/1.2 d).
##          dt = 3.5 d ->   7.6x    a nominal 1000 IU/mL trigger is really 7551
##          dt = 7   d ->  57.0x    ... 57,018 IU/mL
##          dt = 14  d -> 3251x     ... 3,250,997 IU/mL
##      Verified consequence: peak DNAemia 3.44 / 3.94 / 4.62 log10 for q3.5 /
##      q7 / q14, and P(CMV disease) 0.12 / 0.17 / 0.34.
##
##  FOUR STRUCTURAL DECISIONS
##  ---------------------------------------------------------------------------
##  [1] MEASURED DNAemia IS TWO SPECIES.  VVx = encapsidated virion DNA (the
##      only infectious species, the thing BETA multiplies).  VLY = free DNA
##      shed by dying infected cells, ~67% of what a plasma PCR reports.
##      Polymerase inhibitors cut both.  Terminase/kinase inhibitors cut only
##      VVx and, because unit-length genomes are no longer excised from the
##      concatemer, RAISE the DNA per dying cell by AMPPK.  Verified: at a fixed
##      infected-cell pool the immediate change in measured DNAemia is
##          ganciclovir 900 BID  1.52 log10 fall
##          foscarnet   90 q12h  0.55 log10 fall
##          letermovir  480 od   0.00 log10 fall   (0.17 to -0.12 over AMPPK 0-1)
##          maribavir   400 BID  0.00 log10 fall
##      -> on-letermovir DNAemia is not an early efficacy read-out, and
##         letermovir is a prophylaxis drug rather than a treatment drug.
##
##  [2] GANCICLOVIR'S ACTIVATION STEP IS A VIRAL GENE PRODUCT.  pUL97 makes
##      GCV-MP; maribavir inhibits pUL97; so maribavir enters the ganciclovir
##      arm through the ACTIVATION term (F_KIN below) and the pair is
##      antagonistic by construction.  Verified at treatment doses:
##          GCV alone 0.970 · MBV alone 0.792 · both 0.913
##      i.e. the combination is 0.057 WORSE in e than ganciclovir alone.
##      UL97 GCV-resistance is modelled as loss of KINASE efficiency (FK_A =
##      0.125), not as a polymerase EC50 shift, which is why it spares maribavir
##      (MBV resistance sits at different residues: T409M/H411Y/C480F).
##
##  [3] THE NEUTROPENIA LOOP CLOSES THROUGH THE KIDNEY.  CrCl sets ganciclovir
##      exposure (90% renal), exposure drives a Friberg 5-compartment marrow
##      chain, the ANC rule cuts the dose, and the cut can drop e below e*.
##      Verified, correctly dose-banded: e 0.880 -> 0.657 as CrCl falls 100 ->
##      20, with ANC steady state staying above 1.0 throughout.  Verified with
##      the renal adjustment OMITTED: e rises to 0.989 but ANC steady state
##      falls to 0.72, and the simulated arm spends 53 days with ANC < 1.0.
##      The two errors fail in opposite directions.
##
##  [4] RESISTANCE IS A SHAPE, NOT A POTENCY, AND ITS RISK IS A DURATION.
##      A deterministic mutation flux is NOT used: it lets 1e-9 of an infected
##      cell exist and grow at the mutant's positive exponent, which makes
##      resistance certain in every arm by ~day 70 (flatly wrong).  Instead the
##      model integrates DAYS SPENT IN THE SELECTION WINDOW -- the state in which
##      the wild type is held (e_W > e*) while the mutant's own R_eff exceeds 1.
##      Verified: 46 d (pre-emptive q3.5) · 61 d (q7) · 83 d (q14) · 200 d
##      (valGCV prophylaxis 200 d) · 202 d in the UL56 window (letermovir 200 d).
##      Prophylaxis multiplies time-under-selection by 3-4x versus pre-emption.
##
##  CALIBRATION LEDGER (what was fitted, and to what)
##  ---------------------------------------------------------------------------
##    KPROD, DELI      <- doubling time 1.2 d and decline half-life 2.4 d
##    LAMT   = 2.00    <- ONE number fitted to the observed peak DNAemia of
##                        untreated primary D+/R- infection (model 10^4.83)
##    KDIS   = 0.0300  <- ONE number fitted so untreated D+/R- gives
##                        P(CMV disease) = 0.60 at one year
##    GAM/EMAX_MYE/E_MPA  <- chronic-exposure ANC, not the usual Friberg
##                        chemotherapy fits (see note at the Friberg block)
##    everything else  <- taken from published PK/PD or in-vitro EC50 values
##
##  RUN
##    library(mrgsolve); library(dplyr)
##    mod <- mread("cmv_mrgsolve_model.R")
##    out <- cmv_scenario(mod, "S7")      # helper at the bottom of this file
## ===========================================================================

$PROB
# CMV in transplant recipients — 48-ODE QSP model (mirror of cmv_python_reference.py)

$SET end = 365, delta = 0.25, atol = 1e-9, rtol = 1e-7, maxsteps = 100000

## ---------------------------------------------------------------------------
$PARAM @annotated
// ---- 1. systemic viral dynamics -----------------------------------------
DBL0     :  1.20  : untreated plasma DNAemia doubling time (d)
THALFTX  :  2.40  : on-treatment DNAemia decline half-life (d)
CVIR     :  1.20  : encapsidated virion DNA clearance (1/d)
CLYS     :  6.00  : free (lysis-derived) DNA clearance (1/d)
PVIR     :  1e4   : virion DNA produced per unit I per day (IU/mL/d)
RLYS     :  2.00  : Vl:Vv ratio at untreated quasi-steady state
AMPPK    :  0.50  : extra DNA per dying cell when packaging is blocked
LAMT     :  2.00  : permissive target-cell renewal (1/d) [fitted to peak VL]
TC0      :  1.00  : permissive pool (normalised)

// ---- 2. latency and reactivation ----------------------------------------
KREACT   :  3e-5  : baseline reactivation flux per unit reservoir (1/d)
AREACT   :  6.00  : immunosuppression multiplier on reactivation
AINFL    :  3.00  : inflammation multiplier on reactivation
KLAT     :  1e-3  : reservoir re-seeding from active infection

// ---- 3. strain fitness and resistance ----------------------------------
FIT_A    :  0.85  : UL97 activation-mutant replicative fitness
FIT_B    :  0.95  : UL56 C325Y replicative fitness
FK_A     :  0.125 : UL97-mutant kinase efficiency (= 8-fold GCV resistance)
RESBLTV  :  3000  : UL56 C325Y letermovir fold-shift
MU_A     :  3e-5  : per-infection probability of generating UL97 resistance
MU_B     :  3e-5  : per-infection probability of generating UL56 resistance
IEXT     :  1e-8  : extinction floor (one infected cell body-wide)

// ---- 4. adaptive and innate immunity ------------------------------------
KE8      :  0.12  : CMV-specific CD8 killing (per (cell/uL) per d)
KENK     :  0.010 : adaptive NK killing (per (cell/uL) per d)
RHO8     :  0.65  : max CD8 expansion rate (1/d)
KAG      :  0.05  : antigen half-saturation for CD8 expansion (I units)
KHELP    :  1.20  : CD4 help half-saturation (cells/uL)
KISE     :  1.20  : immunosuppression index halving CD8 expansion
DE8      :  0.05  : CD8 effector contraction (1/d)
KMEM     :  0.35  : fraction of contracting effectors entering memory
DEM8     :  0.002 : memory CD8 loss (1/d)
FMEM     :  0.50  : cytotoxic weight of a memory cell vs an effector
E8MAX    : 60.0   : CD8 carrying capacity (cells/uL)
SPRE8    :  0.002 : CD8 precursor influx (cells/uL/d)
RHO4     :  0.55  : max CD4 expansion rate (1/d)
DE4      :  0.012 : CD4 loss (1/d)
E4MAX    : 25.0   : CD4 carrying capacity (cells/uL)
SPRE4    :  0.004 : CD4 precursor influx (cells/uL/d)
KISE4    :  0.80  : immunosuppression index halving CD4 expansion
RHONK    :  0.06  : adaptive NK expansion (1/d)
DNK      :  0.02  : adaptive NK loss (1/d)
NKMAX    : 40.0   : adaptive NK ceiling (cells/uL)
KAGNK    :  0.20  : antigen half-saturation for NK expansion
KISNK    :  4.00  : immunosuppression index halving NK expansion
KAB      :  0.30  : neutralising antibody production
DAB      :  0.02  : neutralising antibody loss (1/d)
KNAB     :  8.00  : antibody titre halving effective infectivity

// ---- 5. immunosuppression ----------------------------------------------
LYM0     :  1.80  : baseline lymphocyte count (10^3/uL)
KDEPL    :  3.00  : ATG-driven lymphodepletion (per (ug/mL) per d)
KREPOP   :  0.004 : lymphocyte repopulation (1/d)
KELATG   :  0.115 : ATG elimination (1/d)
VTAC     : 1200   : tacrolimus apparent volume (L)
CLTAC    : 30.0   : tacrolimus apparent clearance (L/h)
DDITLTV  :  2.40  : letermovir effect on tacrolimus AUC
DDITMBV  :  1.50  : maribavir effect on tacrolimus AUC
KDDILTV  :  0.012 : free letermovir for half-maximal CYP3A/OATP inhibition (uM)
KDDIMBV  :  0.050 : free maribavir for half-maximal CYP3A inhibition (uM)
KSTER    :  0.035 : steroid taper rate (1/d)
WTAC     :  0.55  : ISI weight, tacrolimus
WLYM     :  0.80  : ISI weight, lymphodepletion
WSTER    :  0.25  : ISI weight, corticosteroid
WMPA     :  0.40  : ISI weight, mycophenolate
WMTOR    :  0.35  : ISI credit, mTOR inhibitor

// ---- 6. ganciclovir / valganciclovir -----------------------------------
MWVGCV   : 354.4  : valganciclovir MW
FVGCV    :  0.60  : valganciclovir -> ganciclovir systemic availability
KAGCV    :  8.00  : valganciclovir absorption (1/d)
VGCV     : 60.0   : ganciclovir volume (L)
CLGCV    : 13.0   : ganciclovir clearance at CrCl 100 (L/h)
FRENGCV  :  0.90  : renal fraction of ganciclovir clearance
KPHOS    : 12.0   : pUL97 phosphorylation of GCV (1/d)
KDEGTP   :  0.924 : GCV-TP loss (1/d), t1/2 18 h
EC50GTP  : 23.4   : GCV-TP for 50% polymerase block (uM-eq)
HGTP     :  2.00  : Hill coefficient of chain termination

// ---- 7. letermovir -----------------------------------------------------
MWLTV    : 572.6  : letermovir MW
FLTV     :  0.94  : letermovir bioavailability
KALTV    :  6.00  : letermovir absorption (1/d)
VLTV     : 45.0   : letermovir volume (L)
CLLTV    :  2.60  : letermovir clearance (L/h)
FULTV    :  0.010 : letermovir free fraction
EC50LTV  :  0.0040: free letermovir EC50 vs pUL56 terminase (uM)
PENLTVS  :  0.010 : sanctuary penetration (free drug)

// ---- 8. maribavir ------------------------------------------------------
MWMBV    : 376.2  : maribavir MW
FMBV     :  0.90  : maribavir bioavailability
KAMBV    : 10.0   : maribavir absorption (1/d)
VMBV     : 27.0   : maribavir volume (L)
CLMBV    :  3.50  : maribavir clearance (L/h)
FUMBV    :  0.020 : maribavir free fraction
EC50MBV  :  0.120 : free maribavir EC50 vs pUL97 (uM)
KIMBVK   :  0.120 : Ki for blocking GCV phosphorylation (uM)
PENMBVS  :  0.400 : sanctuary penetration

// ---- 9. foscarnet and cidofovir ---------------------------------------
VFOS     : 40.0   : foscarnet volume (L)
CLFOS    :  6.50  : foscarnet clearance (L/h)
FRENFOS  :  0.95  : renal fraction of foscarnet clearance
EC50FOS  : 250.0  : foscarnet EC50 (uM)
PENFOSS  :  0.200 : sanctuary penetration
VCDV     : 25.0   : cidofovir volume (L)
CLCDV    : 10.0   : cidofovir clearance (L/h)
FRENCDV  :  0.90  : renal fraction of cidofovir clearance
KPPCDV   :  6.00  : cidofovir -> CDV-PP (1/d)
KDEGPP   :  0.30  : CDV-PP loss (1/d)
EC50CPP  :  1.20  : CDV-PP EC50 (uM-eq)

// ---- 10. myelosuppression (Friberg) -----------------------------------
MTT      : 125.0  : mean transit time (h)
GAM      :  0.35  : feedback exponent (chronic-exposure calibration)
CIRC0    :  4.50  : baseline ANC on MMF-based maintenance (10^9/L)
EMAXMYE  :  0.60  : maximal ganciclovir marrow effect
EC50MYE  : 12.0   : plasma GCV for half-maximal marrow effect (uM)
EMPA     :  0.075 : additive mycophenolate marrow effect
KGCSF    :  0.60  : G-CSF acceleration of the transit chain

// ---- 11. kidney -------------------------------------------------------
GFR0     : 86.36  : allograft eGFR set-point (mL/min/1.73)
KRECG    :  0.10  : eGFR return to set-point (1/d)
KINJG    :  0.344 : tubular injury -> eGFR loss
KREPT    :  0.05  : tubular repair (1/d)
AFOS     :  0.50  : foscarnet tubular injury
ACDV     :  0.08  : cidofovir tubular injury
ATAC     :  0.03  : tacrolimus tubular injury above trough 12
KMGIN    :  1.00  : magnesium input
KMGOUT   :  0.50  : magnesium output (1/d)
FMG      :  1.20  : foscarnet magnesium wasting
CFIB     :  0.030 : damage -> eGFR set-point coupling

// ---- 12. tissue and sanctuary ----------------------------------------
LAMTT    :  0.04  : tissue target-cell renewal (1/d)
TT0      :  1.00  : tissue permissive pool
BETATREL :  0.35  : tissue infectivity relative to systemic
LOSST    :  0.22  : tissue infected-cell loss (1/d)
KE8T     :  0.045 : CD8 killing in tissue
PVIRT    :  3e3   : tissue viral DNA production (IU/g/d)
CVIRT    :  1.00  : tissue viral DNA clearance (1/d)
RS       :  0.12  : sanctuary replication rate (1/d)
KSCLR    :  0.05  : sanctuary clearance (1/d)
KE8S     :  0.030 : CD8 clearance in sanctuary

// ---- 13. clinical hazards -------------------------------------------
KDIS     :  0.0300: CMV-disease hazard scale (1/d) [fitted: untreated 0.60]
KD50     : 250.0  : tissue load for half-maximal disease hazard (IU/g)
HDIS     :  1.50  : Hill coefficient of the disease hazard
EDIS     :  6.00  : CD8 level halving the disease hazard
H0REJ    :  3.2e-4: baseline acute-rejection hazard (1/d)
CREJV    :  0.55  : CMV amplification of the rejection hazard
CREJIS   :  0.55  : immunosuppressive protection from rejection

// ---- 14. regimen inputs (set by the scenario helper) -----------------
RVGCV    :  0.0   : valganciclovir input (umol/d)
RGCVIV   :  0.0   : ganciclovir IV input (umol/d)
RLTV     :  0.0   : letermovir input (umol/d)
RMBV     :  0.0   : maribavir input (umol/d)
RFOS     :  0.0   : foscarnet input (umol/d)
RCDV     :  0.0   : cidofovir input (umol/d)
RATG     :  0.0   : ATG input (ug/mL/d)
RTAC     : 6480   : tacrolimus input (ng/d equivalent)
STERSS   :  5.0   : steroid target (mg prednisone equivalent)
MPA      :  1.0   : mycophenolate on/off
MTOR     :  0.0   : mTOR inhibitor on/off
INFL     :  0.0   : inflammation forcing (sepsis / rejection)
GCSF     :  0.0   : G-CSF on/off
ACT      :  0.0   : adoptive CMV-specific T-cell infusion rate
CMVIG    :  0.0   : CMV hyperimmune globulin input
COSTRATE :  0.0   : running cost (currency/d)

// ---- 15. costs (illustrative list-price order of magnitude) ----------
CSTVGCV  : 42.0   : valganciclovir per day
CSTLTV   : 225.0  : letermovir per day
CSTMBV   : 690.0  : maribavir per day
CSTFOS   : 1450.0 : foscarnet per day
CSTCDV   : 980.0  : cidofovir per day
CSTPCR   : 155.0  : one plasma CMV PCR
CSTGCSF  : 340.0  : G-CSF per day

## ---------------------------------------------------------------------------
$CMT @annotated
LAT   : latent reservoir (normalised genome copies)
TC    : permissive target cells (normalised)
IW    : infected cells, wild type
IA    : infected cells, UL97 activation-mutant
IB    : infected cells, UL56 C325Y mutant
VVW   : plasma ENCAPSIDATED virion DNA, wild type (IU/mL)
VVA   : plasma encapsidated virion DNA, UL97 mutant (IU/mL)
VVB   : plasma encapsidated virion DNA, UL56 mutant (IU/mL)
VLY   : plasma LYSIS-DERIVED free DNA, all strains (IU/mL)
TT    : tissue permissive cells (normalised)
IT    : tissue infected cells
VT    : tissue viral load (IU/g)
VS    : sanctuary (retina/CNS) viral load (IU/mL)
E8    : CMV-specific CD8 effectors (cells/uL)
EM8   : CMV-specific CD8 memory (cells/uL)
E4    : CMV-specific CD4 (cells/uL)
NKA   : adaptive NKG2C+ NK cells (cells/uL)
AB    : neutralising antibody (relative units)
LYM   : total lymphocytes (10^3/uL)
TAC   : tacrolimus whole-blood trough (ng/mL)
ATGC  : ATG concentration (ug/mL)
STER  : prednisone-equivalent dose (mg)
AGCV  : valganciclovir gut depot (umol)
GCV   : plasma ganciclovir (uM)
GTP   : intracellular GCV-triphosphate (uM-eq)
ALTV  : letermovir gut depot (umol)
LTV   : plasma letermovir, total (uM)
AMBV  : maribavir gut depot (umol)
MBV   : plasma maribavir, total (uM)
FOS   : plasma foscarnet (uM)
CDV   : plasma cidofovir (uM)
CTP   : intracellular CDV-diphosphate (uM-eq)
PROL  : Friberg proliferating pool
TR1   : Friberg transit 1
TR2   : Friberg transit 2
TR3   : Friberg transit 3
ANC   : circulating neutrophils (10^9/L)
GFR   : allograft eGFR (mL/min/1.73)
TUBI  : proximal tubular injury burden
MG    : serum magnesium (mg/dL)
AUCV  : integral of measured plasma DNAemia (IU*d/mL)
AUCT  : integral of tissue viral load
HZD   : cumulative CMV-disease hazard
HZR   : cumulative acute-rejection hazard
TNEU  : cumulative days with ANC < 1.0
COST  : cumulative cost
MUTA  : days inside the UL97 selection window
MUTB  : days inside the UL56 selection window

## ---------------------------------------------------------------------------
$GLOBAL
#define LN2 0.6931471805599453

// smooth logistic gate; a BOOLEAN "drug present" test is a bug here, because
// 1e-12 of solver round-off latches it on for ever (see the DDI note below)
double sig(double x, double sc) { return 1.0 / (1.0 + exp(-x / sc)); }

// fraction of viral DNA SYNTHESIS removed by polymerase-directed agents
double epol(double gtp, double fk, double fosf, double ctp,
            double ec50g, double hg, double ec50f, double ec50c) {
  double x  = (gtp * fk) / ec50g;
  double xh = pow(x, hg);
  double eg = (x > 0.0) ? xh / (1.0 + xh) : 0.0;
  double ef = (fosf > 0.0) ? fosf / (fosf + ec50f) : 0.0;
  double ec = (ctp  > 0.0) ? ctp  / (ctp  + ec50c) : 0.0;
  return 1.0 - (1.0 - eg) * (1.0 - ef) * (1.0 - ec);
}

// fraction of ENCAPSIDATION / nuclear egress removed (terminase, kinase)
double epack(double ltvf, double mbvf, double resl, double resm,
             double ec50l, double ec50m) {
  double el = (ltvf > 0.0) ? ltvf / (ltvf + ec50l * resl) : 0.0;
  double em = (mbvf > 0.0) ? mbvf / (mbvf + ec50m * resm) : 0.0;
  return 1.0 - (1.0 - el) * (1.0 - em);
}

## ---------------------------------------------------------------------------
$MAIN
// ---- derived constants: the two thresholds the whole model turns on -----
double r0     = LN2 / DBL0;                 // 0.5776 /d
double DELI   = LN2 / THALFTX;              // 0.2888 /d
double KPROD  = r0 + DELI;                  // 0.8664 /d
double BETA0  = KPROD * CVIR / (PVIR * TC0);
double QLYS   = RLYS * PVIR * CLYS / (CVIR * DELI);
double EPSTAR = r0 / KPROD;                 // 0.6667  THRESHOLD 1
double E8STAR = r0 / KE8;                   // 4.81     THRESHOLD 2
double KTR    = 4.0 / (MTT / 24.0);

// ---- default initial conditions: D+/R- kidney recipient after ATG -------
LAT_0  = 0.05;    // graft-borne inoculum; a seronegative recipient has no
TC_0   = TC0;     // memory pool at all, which is the whole of the D+/R- risk
IW_0   = 2.0e-4;
TT_0   = TT0;
E8_0   = 0.004;   // 0.04 naive precursors x 0.10 ATG depletion
EM8_0  = 0.0;
E4_0   = 0.008;
LYM_0  = 0.12;
TAC_0  = 9.0;
ATGC_0 = 60.0;
STER_0 = 20.0;
PROL_0 = CIRC0; TR1_0 = CIRC0; TR2_0 = CIRC0; TR3_0 = CIRC0; ANC_0 = CIRC0;
GFR_0  = GFR0;
MG_0   = 2.0;

## ---------------------------------------------------------------------------
$ODE
// =========================== immunosuppression ==========================
double ISI = WTAC * TAC / 8.0
           + WLYM * (1.0 - LYM / LYM0 > 0.0 ? 1.0 - LYM / LYM0 : 0.0)
           + WSTER * STER / 20.0
           + WMPA * MPA
           - WMTOR * MTOR;
if (ISI < 0.0) ISI = 0.0;

// =========================== free drug =================================
double fosf = FOS;
double ltvf = LTV * FULTV;
double mbvf = MBV * FUMBV;

// per-strain inhibition.  UL97 GCV-resistance enters as loss of KINASE
// efficiency (FK_A), not as a polymerase EC50 shift -- which is exactly why it
// leaves maribavir and foscarnet untouched.
double epW = epol(GTP, 1.0,  fosf, CTP, EC50GTP, HGTP, EC50FOS, EC50CPP);
double epA = epol(GTP, FK_A, fosf, CTP, EC50GTP, HGTP, EC50FOS, EC50CPP);
double epB = epW;
double ekW = epack(ltvf, mbvf, 1.0,     1.0, EC50LTV, EC50MBV);
double ekA = ekW;
double ekB = epack(ltvf, mbvf, RESBLTV, 1.0, EC50LTV, EC50MBV);

// =========================== immune killing ============================
double E8EFF = E8 + FMEM * EM8;
double kill  = KE8 * E8EFF + KENK * NKA;
double lossW = DELI + kill;
double lossA = lossW;
double lossB = lossW;

// =========================== infection =================================
double beta = BETA0 / (1.0 + AB / KNAB);   // antibody neutralisation
double infW = beta * TC * VVW;
double infA = beta * TC * VVA;
double infB = beta * TC * VVB;
double react = KREACT * LAT * (1.0 + AREACT * ISI) * (1.0 + AINFL * INFL);

dxdt_LAT = -1.0e-4 * react + KLAT * (IW + IA + IB);
dxdt_TC  = LAMT * (TC0 - TC) - (infW + infA + infB);

// NOTE ON MUTATION.  A continuous mutation flux into IA/IB is deliberately NOT
// written here.  With MU ~ 3e-5 and a normalisation in which 1 unit of I is
// ~1e8 cells, a deterministic flux lets 1e-9 of an infected cell exist and then
// grow at the mutant's positive exponent, so resistance becomes CERTAIN in
// every arm within 60-90 d -- against ~5% observed on valGCV prophylaxis.  The
// strain compartments are seeded only by an explicit scenario event, and the
// resistance RISK is carried by MUTA/MUTB below.
dxdt_IW = infW + react - lossW * IW;
dxdt_IA = infA - lossA * IA;
dxdt_IB = infB - lossB * IB;

// =============== the two measured DNA species [decision 1] =============
dxdt_VVW = PVIR * IW * (1 - epW) * (1 - ekW) - CVIR * VVW;
dxdt_VVA = PVIR * FIT_A * IA * (1 - epA) * (1 - ekA) - CVIR * VVA;
dxdt_VVB = PVIR * FIT_B * IB * (1 - epB) * (1 - ekB) - CVIR * VVB;

double lysflux = lossW * IW * (1 - epW) * (1 + AMPPK * ekW)
               + lossA * IA * (1 - epA) * (1 + AMPPK * ekA)
               + lossB * IB * (1 - epB) * (1 + AMPPK * ekB);
dxdt_VLY = QLYS * lysflux - CLYS * VLY;

double VVTOT = VVW + VVA + VVB;
double VMEAS = VVTOT + VLY;

// =========================== tissue compartment ========================
double betat = beta * BETATREL;
double inft  = betat * TT * VVTOT;
double losst = LOSST + KE8T * E8EFF;
dxdt_TT = LAMTT * (TT0 - TT) - inft;
dxdt_IT = inft - losst * IT;
dxdt_VT = PVIRT * IT * (1 - epW) * (1 - ekW) - CVIRT * VT;

// ======================== sanctuary (retina / CNS) =====================
double eps_s = epol(GTP * 0.30, 1.0, fosf * PENFOSS, CTP * 0.15,
                    EC50GTP, HGTP, EC50FOS, EC50CPP);
double ekk_s = epack(ltvf * PENLTVS, mbvf * PENMBVS, 1.0, 1.0,
                     EC50LTV, EC50MBV);
double gs = RS * (1 - eps_s) * (1 - ekk_s);
dxdt_VS = gs * VS * (1.0 - VS / 1.0e6) - (KSCLR + KE8S * E8EFF) * VS
          + 1.0e-3 * VVTOT;

// ==================== adaptive and innate immunity =====================
// THE COUPLING: the expansion term is proportional to ANTIGEN, i.e. to the
// infected-cell pool the antiviral has just removed.
double AG   = IW + IA + IB + 0.5 * IT;
double sat  = AG / (AG + KAG);
double hlp  = E4 / (E4 + KHELP);
double fis  = 1.0 / (1.0 + ISI / KISE);
double fis4 = 1.0 / (1.0 + ISI / KISE4);
double fnk  = 1.0 / (1.0 + ISI / KISNK);

dxdt_E8  = RHO8 * E8EFF * sat * hlp * fis * (1 - E8EFF / E8MAX)
           - DE8 * E8 + SPRE8 + 0.25 * ACT;
dxdt_EM8 = KMEM * DE8 * E8 - DEM8 * EM8;
dxdt_E4  = RHO4 * E4 * sat * fis4 * (1 - E4 / E4MAX) - DE4 * E4 + SPRE4;
dxdt_NKA = RHONK * (AG / (AG + KAGNK)) * fnk * (NKMAX - NKA) - DNK * NKA;
dxdt_AB  = KAB * sat * hlp * fis4 - DAB * AB + CMVIG;

// ==================== immunosuppression pharmacology ===================
dxdt_LYM  = -KDEPL * ATGC * LYM + KREPOP * (LYM0 - LYM);
dxdt_ATGC = RATG - KELATG * ATGC;

// CYP3A / OATP1B1 inhibition MUST be concentration-dependent.  Written as a
// boolean "if (LTV > 0) clearance /= 2.4" it is latched permanently by solver
// round-off and silently trebles the tacrolimus trough.
double ddi = 1.0
           + (DDITLTV - 1.0) * ltvf / (ltvf + KDDILTV)
           + (DDITMBV - 1.0) * mbvf / (mbvf + KDDIMBV);
double cltac = (CLTAC * 24.0) / ddi;
dxdt_TAC  = RTAC / VTAC - cltac / VTAC * TAC;
dxdt_STER = -KSTER * (STER - STERSS);

// ============================ antiviral PK =============================
double crcl = GFR * 1.10;  if (crcl < 8.0) crcl = 8.0;

double clgcv = (CLGCV * (FRENGCV * crcl / 100.0 + (1 - FRENGCV))) * 24.0;
dxdt_AGCV = RVGCV - KAGCV * AGCV;
dxdt_GCV  = KAGCV * AGCV * FVGCV / VGCV + RGCVIV / VGCV - clgcv / VGCV * GCV;

// [decision 2] pUL97 makes GCV-MP, and maribavir inhibits pUL97, so maribavir
// enters the ganciclovir arm through the ACTIVATION term
double fkin = 1.0 / (1.0 + mbvf / KIMBVK);
dxdt_GTP  = KPHOS * GCV * fkin - KDEGTP * GTP;

dxdt_ALTV = RLTV - KALTV * ALTV;
dxdt_LTV  = KALTV * ALTV * FLTV / VLTV - CLLTV * 24.0 / VLTV * LTV;
dxdt_AMBV = RMBV - KAMBV * AMBV;
dxdt_MBV  = KAMBV * AMBV * FMBV / VMBV - CLMBV * 24.0 / VMBV * MBV;

double clfos = (CLFOS * (FRENFOS * crcl / 100.0 + (1 - FRENFOS))) * 24.0;
dxdt_FOS = RFOS / VFOS - clfos / VFOS * FOS;
double clcdv = (CLCDV * (FRENCDV * crcl / 100.0 + (1 - FRENCDV))) * 24.0;
dxdt_CDV = RCDV / VCDV - clcdv / VCDV * CDV;
dxdt_CTP = KPPCDV * CDV - KDEGPP * CTP;

// ============== myelosuppression: Friberg 5-compartment ================
// GAM = 0.35 rather than the usual 0.16-0.24.  With CHRONIC dosing the steady
// state is ANC = CIRC0*(1-Edrug)^(1/GAM), so GAM = 0.17 puts a patient on
// mycophenolate alone at ANC 0.83 -- the classic Friberg fits were made on
// transient chemotherapy exposures and do not transfer.  GAM, EMAXMYE and EMPA
// were re-fitted against chronic-exposure ANC instead: valGCV 900 od + MMF ->
// ANC 1.99, letermovir + MMF -> ANC 3.61, valGCV 900 od at CrCl 30 without
// renal adjustment -> ANC 0.97.
double edrug = EMAXMYE * GCV / (EC50MYE + GCV) + EMPA * MPA;
if (edrug > 0.95) edrug = 0.95;
double ktr = KTR * (1.0 + KGCSF * GCSF);
double ancf = (ANC > 0.05 ? ANC : 0.05);
dxdt_PROL = ktr * PROL * ((1.0 - edrug) * pow(CIRC0 / ancf, GAM) - 1.0);
dxdt_TR1  = ktr * (PROL - TR1);
dxdt_TR2  = ktr * (TR1 - TR2);
dxdt_TR3  = ktr * (TR2 - TR3);
dxdt_ANC  = ktr * TR3 - ktr * ANC;

// ================================ kidney ===============================
// [decision 3, second half] crcl above is what closes the loop back onto the
// ganciclovir exposure that drives edrug
double tacex = (TAC > 12.0 ? TAC - 12.0 : 0.0);
dxdt_TUBI = AFOS * FOS / 1000.0 + ACDV * CTP / 10.0 + ATAC * tacex / 4.0
            - KREPT * TUBI;
double damage = HZR + 0.15 * log10(1.0 + VMEAS) / 10.0;
double gfrset = GFR0 * exp(-CFIB * damage * 10.0);
dxdt_GFR = KRECG * (gfrset - GFR) - KINJG * TUBI;
dxdt_MG  = KMGIN - KMGOUT * MG * (1.0 + FMG * FOS / 400.0);

// =========================== endpoint accumulators =====================
dxdt_AUCV = VMEAS;
dxdt_AUCT = VT;
double vth = pow(VT, HDIS);
dxdt_HZD = KDIS * (vth / (vth + pow(KD50, HDIS))) / (1.0 + E8EFF / EDIS);
dxdt_HZR = H0REJ * (1.0 + CREJV * log10(1.0 + VMEAS)) * exp(-CREJIS * ISI);
dxdt_TNEU = 1.0 / (1.0 + exp((ANC - 1.0) / 0.05));
dxdt_COST = COSTRATE;

// ============ [decision 4] resistance risk = TIME UNDER SELECTION ======
// The only state that converts a minority variant into the dominant one is:
// the wild type is held (e_W > e*) while the mutant's own R_eff still exceeds 1.
double ReffA = (KPROD * FIT_A * (1 - epA) * (1 - ekA)) / (lossA > 1e-9 ? lossA : 1e-9);
double ReffB = (KPROD * FIT_B * (1 - epB) * (1 - ekB)) / (lossB > 1e-9 ? lossB : 1e-9);
double eWtot = 1.0 - (1 - epW) * (1 - ekW);
double gateW = sig(eWtot - EPSTAR, 0.02);
dxdt_MUTA = gateW * sig(ReffA - 1.0, 0.05);
dxdt_MUTB = gateW * sig(ReffB - 1.0, 0.05);

## ---------------------------------------------------------------------------
$TABLE
double r0o     = LN2 / DBL0;
double DELIo   = LN2 / THALFTX;
double EPSTARo = r0o / (r0o + DELIo);
double E8STARo = r0o / KE8;

double fosfo = FOS;
double ltvfo = LTV * FULTV;
double mbvfo = MBV * FUMBV;
double epWo = epol(GTP, 1.0, fosfo, CTP, EC50GTP, HGTP, EC50FOS, EC50CPP);
double ekWo = epack(ltvfo, mbvfo, 1.0, 1.0, EC50LTV, EC50MBV);

double VMEASo = VVW + VVA + VVB + VLY;
double E8EFFo = E8 + FMEM * EM8;
double ISIo = WTAC * TAC / 8.0
            + WLYM * (1.0 - LYM / LYM0 > 0.0 ? 1.0 - LYM / LYM0 : 0.0)
            + WSTER * STER / 20.0 + WMPA * MPA - WMTOR * MTOR;
if (ISIo < 0.0) ISIo = 0.0;

capture DNAEMIA  = VMEASo;                        // IU/mL, what a PCR reports
capture LOG10VL  = log10(VMEASo > 1.0 ? VMEASo : 1.0);
capture VIRIONFR = (VMEASo > 0 ? (VVW+VVA+VVB)/VMEASo : 0.0);
capture EPS_TOT  = 1.0 - (1 - epWo) * (1 - ekWo); // the number vs EPSTAR
capture EPSTARC  = EPSTARo;
capture MARGIN   = (1.0 - (1 - epWo) * (1 - ekWo)) - EPSTARo;
capture E8EFFC   = E8EFFo;
capture E8STARC  = E8STARo;
capture IMMMARG  = E8EFFo - E8STARo;
capture ISIC     = ISIo;
capture CRCLC    = (GFR * 1.10 < 8.0 ? 8.0 : GFR * 1.10);
capture PDIS     = 1.0 - exp(-HZD);               // P(CMV end-organ disease)
capture PREJ     = 1.0 - exp(-HZR);               // P(acute rejection)
capture RESFRACA = ((VVW+VVA+VVB) > 1.0 ? VVA/(VVW+VVA+VVB) : 0.0);
capture RESFRACB = ((VVW+VVA+VVB) > 1.0 ? VVB/(VVW+VVA+VVB) : 0.0);
capture SELWINA  = MUTA;                          // days under UL97 selection
capture SELWINB  = MUTB;                          // days under UL56 selection

## ---------------------------------------------------------------------------
$CAPTURE DNAEMIA LOG10VL EPS_TOT MARGIN E8EFFC IMMMARG ISIC ANC GFR MG PDIS PREJ

## ===========================================================================
##  R DRIVER
##  ---------------------------------------------------------------------------
##  The discrete parts of the clinical strategy -- when a PCR is drawn, what the
##  trigger does, the ANC dose-modification rule, the guideline definition of
##  refractory CMV -- are decisions, not differential equations.  They live here,
##  in a loop that advances the ODE system one 6-hour window at a time and
##  updates the parameter vector between windows.  This is exactly what
##  cmv_python_reference.py does, so the two implementations stay comparable.
## ===========================================================================
/*
library(mrgsolve); library(dplyr)
mod <- mread("cmv_mrgsolve_model.R")

MW <- list(VGCV = 354.4, LTV = 572.6, MBV = 376.2, FOS = 126.0, CDV = 279.2)

## label-consistent valganciclovir dosing, mg/day averaged over the interval
vgcv_mg_per_day <- function(crcl, level = 1, treat = FALSE) {
  if (treat) {
    b <- if (crcl >= 60) c(900, 0.5) else if (crcl >= 40) c(450, 0.5) else
         if (crcl >= 25) c(450, 1.0) else c(450, 2.0)
  } else {
    b <- if (crcl >= 60) c(900, 1.0) else if (crcl >= 40) c(450, 1.0) else
         if (crcl >= 25) c(450, 2.0) else c(450, 3.5)
  }
  b[1] * level / b[2]
}

## initial condition sets by donor/recipient serostatus and induction
init_cmv <- function(mod, sero = "D+/R-", induction = "ATG", crcl = 95) {
  i <- as.list(init(mod))
  i$GFR <- crcl / 1.10
  if (sero == "D+/R-") {
    i$LAT <- 0.05; i$IW <- 2e-4; i$E8 <- 0.04; i$EM8 <- 0;   i$E4 <- 0.08; i$AB <- 0
  } else if (sero == "D+/R+") {
    i$LAT <- 1.00; i$IW <- 5e-5; i$E8 <- 1.00; i$EM8 <- 20;  i$E4 <- 7.00; i$AB <- 12
  } else if (sero == "D-/R+") {
    i$LAT <- 1.00; i$IW <- 2e-5; i$E8 <- 1.00; i$EM8 <- 22;  i$E4 <- 8.00; i$AB <- 12
  } else {                                     # D-/R-
    i$LAT <- 0.00; i$IW <- 0.00; i$E8 <- 0.04; i$EM8 <- 0;   i$E4 <- 0.08; i$AB <- 0
  }
  f <- switch(induction, ATG = 0.10, alemtuzumab = 0.04, 1.00)
  if (induction %in% c("ATG", "alemtuzumab")) {
    i$ATGC <- if (induction == "ATG") 60 else 90
    i$LYM  <- if (induction == "ATG") 0.12 else 0.05
    i$E8 <- i$E8 * f; i$EM8 <- i$EM8 * f; i$E4 <- i$E4 * f
  } else i$ATGC <- 0
  i
}

cmv_sim <- function(mod,
                    strategy      = "prophylaxis",  # none | preemptive | prophylaxis
                    proph_drug    = "VGCV",         # VGCV | LTV
                    proph_days    = 200,
                    monitor_int   = 14,             # PCR interval, days
                    trigger       = 1000,           # IU/mL
                    lloq          = 137,
                    tx_drug       = "VGCV",
                    tx_drug2      = NULL,           # salvage on refractory CMV
                    tx_min_days   = 0,              # secondary prophylaxis
                    stop_rule     = 2,              # consecutive negatives
                    post_monitor  = TRUE,
                    renal_adjust  = TRUE,
                    tac_cut_ltv   = TRUE,
                    adhere        = 1.0,
                    resist_seed   = 0,              # minority UL97 variant
                    gcsf_policy   = TRUE,
                    mpa = 1, mtor = 0,
                    sero = "D+/R-", induction = "ATG", crcl = 95,
                    tend = 365, dt = 0.25) {

  y <- init_cmv(mod, sero, induction, crcl)
  p <- as.list(param(mod)); p$GFR0 <- crcl / 1.10
  st <- list(on_tx = FALSE, agent = NA, neg = 0, ok = 0, level = 1, gcsf = 0,
             npcr = 0, txd = 0, seeded = FALSE, switched = NA,
             tx_t0 = 0, tx_v0 = 0, nextpcr = 0, log = character(0))
  keep <- list(); t <- 0

  while (t < tend - 1e-9) {
    ## ---------------- monitoring visit ---------------------------------
    monitoring <- (strategy == "preemptive") ||
      (strategy == "prophylaxis" && post_monitor && t >= proph_days) || st$on_tx
    if (monitoring && t >= st$nextpcr - 1e-9) {
      st$npcr <- st$npcr + 1
      y$COST  <- y$COST + p$CSTPCR
      vobs <- y$VVW + y$VVA + y$VVB + y$VLY
      if (vobs < lloq) vobs <- 0
      if (!st$on_tx) {
        if (vobs >= trigger) {
          st$on_tx <- TRUE; st$agent <- tx_drug; st$neg <- 0
          st$tx_t0 <- t; st$tx_v0 <- vobs
          if (resist_seed > 0 && !st$seeded) {
            ## a minority variant selected during an EARLIER course.  Seeding at
            ## transplant does nothing: at the drug-free target-limited
            ## equilibrium the fitness-costed mutant is competitively excluded.
            y$IA <- y$IA + resist_seed * y$IW
            y$VVA <- y$VVA + resist_seed * y$VVW
            st$seeded <- TRUE
          }
          st$log <- c(st$log, sprintf("d%.0f start %s at %.0f IU/mL",
                                      t, st$agent, vobs))
        }
      } else {
        if (vobs == 0) {
          st$neg <- st$neg + 1
          if (st$neg >= stop_rule && (t - st$tx_t0) >= tx_min_days) {
            st$on_tx <- FALSE; st$agent <- NA; st$level <- 1
            st$log <- c(st$log, sprintf("d%.0f stop antiviral", t))
          }
        } else st$neg <- 0
        ## refractory CMV, guideline definition: >= 14 d of therapy with
        ## < 0.5 log10 fall from the load at which therapy started
        if (is.na(st$switched) && !is.null(tx_drug2) && vobs >= trigger &&
            (t - st$tx_t0) >= 14 && vobs > 0.316 * st$tx_v0) {
          st$switched <- t
          if (tx_drug2 != st$agent) {
            st$agent <- tx_drug2
            st$log <- c(st$log, sprintf("d%.0f REFRACTORY -> switch to %s",
                                        t, st$agent))
          } else st$log <- c(st$log, sprintf("d%.0f REFRACTORY -> stay on %s",
                                             t, st$agent))
        }
      }
      st$nextpcr <- t + monitor_int
    }

    ## ---------------- weekly CBC and the ANC dose rule ------------------
    if (abs(t / 7 - round(t / 7)) < 1e-9) {
      if (y$ANC < 0.5) {
        st$level <- 0; st$ok <- 0; st$gcsf <- as.numeric(gcsf_policy)
        st$log <- c(st$log, sprintf("d%.0f ANC %.2f -> hold + G-CSF", t, y$ANC))
      } else if (y$ANC < 1.0) {
        st$level <- 0.5; st$ok <- 0; st$gcsf <- as.numeric(gcsf_policy)
      } else if (y$ANC > 1.5) {
        st$ok <- st$ok + 1
        if (st$ok >= 2) { st$level <- 1; st$gcsf <- 0 }
      } else st$ok <- 0
    }

    ## ---------------- build the parameter vector for this window --------
    crcl_now <- max(8, y$GFR * 1.10)
    tgt <- if (t < 90) 9 else 6
    ddi_exp <- if (tac_cut_ltv && strategy == "prophylaxis" &&
                   proph_drug == "LTV" && t < proph_days) p$DDITLTV else 1
    pp <- p
    pp$RTAC   <- tgt * (p$CLTAC * 24) / ddi_exp
    pp$STERSS <- if (t < 7) 20 else if (t < 30) 10 else 5
    pp$MPA <- mpa; pp$MTOR <- mtor; pp$GCSF <- st$gcsf
    pp[c("RVGCV","RLTV","RMBV","RFOS","RCDV","COSTRATE")] <- 0

    drug <- NA; treat <- FALSE
    if (strategy == "prophylaxis" && t < proph_days) drug <- proph_drug
    else if (st$on_tx) { drug <- st$agent; treat <- TRUE }
    adh <- if (treat) 1 else adhere

    if (!is.na(drug)) {
      if (drug == "VGCV" || drug == "GCV+MBV") {
        mgd <- adh * vgcv_mg_per_day(if (renal_adjust) crcl_now else 100,
                                     st$level, treat)
        pp$RVGCV <- mgd / MW$VGCV * 1000
        pp$COSTRATE <- pp$COSTRATE +
          p$CSTVGCV * (if (treat && crcl_now >= 60) 2 else 1)
        if (drug == "GCV+MBV") {
          pp$RMBV <- 800 / MW$MBV * 1000
          pp$COSTRATE <- pp$COSTRATE + p$CSTMBV
        }
      } else if (drug == "LTV") {
        pp$RLTV <- adh * 480 / MW$LTV * 1000
        pp$COSTRATE <- pp$COSTRATE + p$CSTLTV
      } else if (drug == "MBV") {
        pp$RMBV <- 800 / MW$MBV * 1000
        pp$COSTRATE <- pp$COSTRATE + p$CSTMBV
      } else if (drug == "FOS") {
        pp$RFOS <- 90 * 70 * 2 / MW$FOS * 1000 * min(1, crcl_now / 100)
        pp$COSTRATE <- pp$COSTRATE + p$CSTFOS
      } else if (drug == "CDV") {
        pp$RCDV <- (5 * 70 / 7) / MW$CDV * 1000
        pp$COSTRATE <- pp$COSTRATE + p$CSTCDV
      }
      st$txd <- st$txd + dt
    }
    if (st$gcsf > 0) pp$COSTRATE <- pp$COSTRATE + p$CSTGCSF

    ## ---------------- advance the ODEs one window ----------------------
    o <- mod %>% param(pp) %>% init(unlist(y)) %>%
      mrgsim(start = 0, end = dt, delta = dt, obsonly = FALSE) %>% as.data.frame()
    last <- o[nrow(o), ]
    for (nm in names(y)) if (nm %in% names(last)) y[[nm]] <- max(0, last[[nm]])

    ## EXTINCTION FLOOR.  A strain compartment at exactly zero is an UNSTABLE
    ## equilibrium whenever that strain's exponent is positive, so ODE solver
    ## round-off (~1e-16) grows at 0.4-0.5/d and is macroscopic in ~70 days --
    ## the model then "discovers" resistance out of floating point.  A strain
    ## below one infected cell body-wide is extinct, not rare.
    if (y$IA < p$IEXT && y$VVA < 1) { y$IA <- 0; y$VVA <- 0 }
    if (y$IB < p$IEXT && y$VVB < 1) { y$IB <- 0; y$VVB <- 0 }
    if (y$IW < p$IEXT && y$VVW < 1) { y$IW <- 0; y$VVW <- 0 }

    last$time <- t + dt
    keep[[length(keep) + 1]] <- last
    t <- t + dt
  }
  list(sim = bind_rows(keep), state = y, meta = st)
}

## ---- the sixteen reported arms ------------------------------------------
CMV_SCENARIOS <- list(
  S1  = list(strategy="none"),
  S2  = list(strategy="preemptive", monitor_int=7),
  S3  = list(strategy="preemptive", monitor_int=3.5),
  S4  = list(strategy="preemptive", monitor_int=14),
  S5  = list(strategy="prophylaxis", proph_drug="VGCV", proph_days=100,
             post_monitor=FALSE),
  S6  = list(strategy="prophylaxis", proph_drug="VGCV", proph_days=200,
             post_monitor=FALSE),
  S7  = list(strategy="prophylaxis", proph_drug="VGCV", proph_days=200,
             monitor_int=14),
  S8  = list(strategy="prophylaxis", proph_drug="LTV",  proph_days=200,
             monitor_int=14),
  S9  = list(strategy="prophylaxis", proph_drug="LTV",  proph_days=200,
             monitor_int=14, tac_cut_ltv=FALSE),
  S10 = list(strategy="prophylaxis", proph_drug="VGCV", proph_days=200,
             monitor_int=14, crcl=30),
  S11 = list(strategy="prophylaxis", proph_drug="VGCV", proph_days=200,
             monitor_int=14, crcl=30, renal_adjust=FALSE),
  S12 = list(strategy="preemptive", monitor_int=7, sero="D+/R+"),
  S13 = list(strategy="preemptive", monitor_int=7, sero="D-/R+",
             induction="basiliximab"),
  S14 = list(strategy="preemptive", monitor_int=7, sero="D+/R+",
             induction="alemtuzumab"),
  S15 = list(strategy="prophylaxis", proph_drug="VGCV", proph_days=200,
             monitor_int=14, mpa=0, mtor=1),
  S16 = list(strategy="preemptive", monitor_int=7, resist_seed=0.20,
             tx_min_days=42, tx_drug2="MBV"),
  S17 = list(strategy="preemptive", monitor_int=7, resist_seed=0.20,
             tx_min_days=42, tx_drug2="FOS"),
  S18 = list(strategy="preemptive", monitor_int=7, resist_seed=0.20,
             tx_min_days=42, tx_drug2="GCV+MBV"),
  S19 = list(strategy="preemptive", monitor_int=7, resist_seed=0.20,
             tx_min_days=42, tx_drug2="VGCV")
)

cmv_scenario <- function(mod, id, ...) {
  do.call(cmv_sim, c(list(mod = mod), CMV_SCENARIOS[[id]], list(...)))
}

cmv_summary <- function(r) {
  s <- r$sim; y <- r$state
  data.frame(
    peak_log10VL = max(s$LOG10VL),
    day_of_peak  = s$time[which.max(s$LOG10VL)],
    P_disease    = 1 - exp(-y$HZD),
    P_rejection  = 1 - exp(-y$HZR),
    E8eff_d200   = s$E8EFFC[which.min(abs(s$time - 200))],
    E8star       = unname(param(mod)$KE8) * 0 + 0.6931472/1.2/0.12,
    ANC_nadir    = min(s$ANC),
    days_ANC_lt1 = y$TNEU,
    selwin_UL97  = y$MUTA,
    selwin_UL56  = y$MUTB,
    eGFR_nadir   = min(s$GFR),
    Mg_nadir     = min(s$MG),
    n_PCR        = r$meta$npcr,
    drug_days    = r$meta$txd,
    cost         = y$COST
  )
}

## ---- the three algebraic read-outs, no integration needed ---------------
cmv_thresholds <- function(mod) {
  p <- as.list(param(mod))
  r0 <- log(2) / p$DBL0; deli <- log(2) / p$THALFTX
  epsstar <- r0 / (r0 + deli)
  cat(sprintf("THRESHOLD 1 (drug)   e*  = %.4f\n", epsstar))
  cat(sprintf("THRESHOLD 2 (immune) E8* = %.2f CMV-specific CD8/uL\n", r0 / p$KE8))
  cat("\nmonitoring interval expressed as a dose:\n")
  for (d in c(2, 3.5, 7, 10, 14, 21))
    cat(sprintf("  q%-5.1fd  %10.1f-fold rise between draws  ->  a nominal"
                , d, 2^(d / p$DBL0)),
        sprintf(" 1000 IU/mL trigger is really %10.0f\n", 1000 * 2^(d / p$DBL0)))
  invisible(epsstar)
}
*/
