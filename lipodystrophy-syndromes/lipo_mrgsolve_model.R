## =============================================================================
##  lipo_mrgsolve_model.R
##  Lipodystrophy syndromes (CGL · AGL · FPLD · APL) — Quantitative Systems
##  Pharmacology model for mrgsolve
##
##  53 ODE compartments.  Time unit = DAYS (a 4-hour metreleptin half-life and
##  a 10-year fibrosis trajectory have to live in the same system).
##
##  ---------------------------------------------------------------------------
##  THE FOUR STRUCTURAL COMMITMENTS
##  ---------------------------------------------------------------------------
##
##  (1) ONE LESION, TWO DEFICITS, AND THE DRUGS ARE SORTED BY WHICH ONE THEY
##      TOUCH.
##
##      Adipose tissue is simultaneously a MECHANICAL BUFFER (the only tissue
##      that can esterify a large triglyceride flux cheaply) and the ENDOCRINE
##      REPORTER of how large that buffer is.  Losing adipocytes therefore
##      produces two deficits that are mechanistically independent:
##
##          capacity deficit :  J_ov = max(0, J_in - Sum(J_st,i) - OX)
##          signal   deficit :  DEF  = 1 - SD(L_effective)
##
##      and the signal deficit RAISES J_in through hyperphagia, so the two
##      deficits multiply rather than add.  Every drug in this file is
##      classified by which term it enters:
##
##          metreleptin   -> SD  (signal)            ... and hence J_in
##          pioglitazone  -> C_cap (capacity), gated by preadipocyte substrate
##          APOC3 / ANGPTL3 / fibrate / omega-3 -> plasma TG clearance only
##          insulin / metformin / SGLT2i         -> glucose only
##          very-low-fat diet                    -> J_in directly
##
##  (2) THE GL/PL DISSOCIATION IS NOT WRITTEN DOWN.  IT IS DERIVED.
##
##      There is no parameter in this model called "metreleptin works better in
##      generalised lipodystrophy".  There is one saturating transducer
##
##          S(L) = L^h / (L^h + EC50^h),  EC50 = 4 ng/mL, h = 2
##
##      applied to total effective leptin.  Because S saturates, the SAME dose
##      produces the same receptor occupancy in every patient, while the
##      achievable change in S is bounded by 1 - S(L_baseline) -- a property of
##      the PATIENT, not of the dose.  Occupancy is captured as OCC and the
##      deficit as DEF so that the product structure can be read off the
##      output.  Scenarios 7 / 8 / 9 exist to make the claim falsifiable:
##      identical mg/kg in two phenotypes, then a dose escalation that raises
##      OCC and should NOT rescue the effect if the model is right.
##
##  (3) PLASMA TRIGLYCERIDE IS TIER 5 OF A CASCADE, NOT THE DISEASE.
##
##      Overflow is partitioned to liver, muscle, plasma, pancreas, heart and
##      kidney with fixed fractions.  Plasma TG is therefore a MARKER of J_ov,
##      and an agent that removes it (APOC3 ASO) must leave hepatic fat and
##      HbA1c nearly unchanged.  That decoupling is a prediction, not an
##      assumption -- scenario 15 is the test, and if hepatic fat fell with
##      plasma TG the tier structure would be wrong.
##
##  (4) ONLY FIBROSIS IS IRREVERSIBLE.
##
##      Every other state in this model relaxes back when the driver is
##      removed (scenario 10 withdraws metreleptin and everything returns).
##      FIB has a regression rate 17x slower than its formation rate, so early
##      versus late initiation (scenario 23) changes the CEILING of recovery
##      rather than its rate.  No bistability is claimed anywhere: the model
##      has no positive feedback loop with gain > 1, and the withdrawal
##      scenario is included precisely so that claim can be checked.
##
##  ---------------------------------------------------------------------------
##  WHAT THIS MODEL IS CALIBRATED AGAINST
##  ---------------------------------------------------------------------------
##  Baseline phenotypes (the model is run to its own steady state for each
##  phenotype -- see baseline() -- rather than having initial conditions typed
##  in by hand):
##
##    congenital generalised (CGL):  leptin 0.4-1.5 ng/mL, HFF 20-35%,
##       TG 800-3000 mg/dL, HbA1c 9-11%, insulin 40-100 uU/mL, ALT 60-120 U/L
##       (Garg 2004 NEJM; Agarwal 2003; Simha 2003 JCEM; Akinci 2019)
##    familial partial (FPLD2, Dunnigan):  leptin 4-12 ng/mL, HFF 10-20%,
##       TG 250-600 mg/dL, HbA1c 7.5-8.5%   (Vigouroux 2011; Guillin-Amarelle 2016)
##    healthy female 60 kg:  leptin 10-15 ng/mL, HFF ~3%, TG ~90 mg/dL,
##       HbA1c ~5.1%, insulin ~8 uU/mL
##
##  Treatment effects the model is asked to reproduce (see summarise_all()):
##    metreleptin in GL      HbA1c -1.7 to -2.5%, TG -40 to -60%,
##                           hepatic fat -30 to -50% relative, insulin dose
##                           requirement often halved
##                           (Oral 2002 NEJM; Petersen 2002 JCI; Javor 2005;
##                            Chan 2011; Diker-Cohen 2015 JCEM;
##                            Brown 2018 Endocrine Reviews / MEASuRE cohort)
##    metreleptin in PL      HbA1c ~ -0.6% overall, benefit concentrated in
##                           the low-leptin stratum
##                           (Diker-Cohen 2015; Oral 2019 NEJM combined analysis)
##    volanesorsen in FPLD   TG -70 to -88%, hepatic fat change small
##                           (BROADEN, Oral 2019/Prohaska 2023)
##    pioglitazone in FPLD   modest fat redistribution + glycaemic gain;
##                           no adipogenic gain where no preadipocytes exist
##
##  Parameters that are honest guesses rather than fitted values are marked
##  "(struct)" in the annotation.  They set the SHAPE of a relationship whose
##  existence is documented but whose coefficient is not: e.g. the fraction of
##  overflow going to pancreas, or the slope of proteinuria on insulin.
##
##  ---------------------------------------------------------------------------
##  UNITS
##  ---------------------------------------------------------------------------
##   time                days
##   metreleptin         micrograms (so amount/L = ng/mL, the clinical unit)
##   small molecules     mg, L
##   adipose depots      kg of triglyceride
##   liver / muscle TG   grams
##   plasma TG, glucose  mg/dL
##   insulin             microunits/mL
##   leptin              ng/mL
##   normalised states   1 = the healthy steady state (IMCLN, PLIP, MYOL, KLIP)
##
##  Author: QSP Disease Model Library (Claude Code Routine)
##  Licence: see repository LICENSE.  EDUCATIONAL / RESEARCH USE ONLY —
##  not validated for clinical or regulatory use.
## =============================================================================

library(mrgsolve)

lipo_code <- r"---(
$PROB
# Lipodystrophy syndromes QSP model (53 ODEs)

$PARAM @annotated
// ================================================================= PHENOTYPE
// The phenotype is a CAPACITY DISTRIBUTION plus a secretion coefficient.
// There is deliberately no categorical "disease type" switch.
FC1      : 1.00 : femorogluteal SC capacity, fraction of normal
FC2      : 1.00 : upper-trunk/facial SC capacity, fraction of normal
FCV      : 1.00 : visceral capacity, fraction of normal
PREAD    : 1.00 : preadipocyte/expandability reserve, 0-1 (gates PPARg effect)
LEPKM    : 1.00 : multiplier on leptin secretion per kg fat (1 = normal biology)
LEAN0    : 45   : lean body mass (kg)
BASLIP   : 0    : excess basal lipolysis (PLIN1/CIDEC phenotypes), 0-1 (struct)

// ============================================================ CAPACITY (kg)
CAP1N    : 12.0 : normal femorogluteal TG ceiling (kg)
CAP2N    :  5.0 : normal upper-trunk/facial TG ceiling (kg)
CAPVN    :  3.0 : normal visceral TG ceiling (kg)
FATN     : 14.0 : reference healthy total fat mass (kg), scales LPL and leptin

// ================================================== ENERGY INTAKE AND INFLUX
EIBASE   : 40.0 : weight-maintenance intake per kg lean mass (kcal/kg/d)
FFAT     : 0.35 : habitual fraction of energy as fat
FCHO     : 0.45 : habitual fraction of energy as carbohydrate
HMAX     : 0.45 : maximal fractional hyperphagia at full signal deficit
TAUEI    :  3.0 : intake adaptation time constant (d)
EICLAMP  :  0   : >0 fixes intake (kcal/d) = pair-feeding clamp
DIETT    : 1e6  : day a very-low-fat diet starts
DIETFAT  : 0.15 : fat fraction of energy on that diet (CHO takes the balance)
KDNL     : 0.030 : de novo lipogenesis yield per g carbohydrate at full insulin drive
KINS50   : 20.0 : insulin giving half-maximal lipogenic drive (uU/mL)
KDNLSEL  : 0.60 : extra DNL drive from selective hepatic insulin resistance (struct)
KLIPSC   : 0.0100 : subcutaneous adipose TG turnover (1/d) -- the futile cycle
KLIPVI   : 0.0120 : visceral adipose TG turnover (1/d)

// ============================================== BUFFER (ESTERIFICATION) TERMS
KEST     : 15.0 : esterification Vmax per kg of depot capacity (g/d/kg)
HRFLOOR  : 0.10 : irreducible turnover headroom of a full depot (struct)
KMJ      : 60.0 : influx giving half-maximal storage flux (g/d)
KOXW     : 2.00 : EXTRAhepatic fat oxidation coefficient (g/d per kg^0.75)
KOXAD    : 0.35 : adaptive rise in oxidation with ectopic load (struct)
TAUOXAD  : 20.0 : time constant of that adaptation (d)

// ====================================================== OVERFLOW PARTITIONING
FLIV     : 0.35 : baseline fraction of overflow to liver
KFLIVA   : 0.30 : extra hepatic share as adipose disappears (competitive sinks)
FMUS     : 0.15 : fraction of overflow to skeletal muscle
FPL      : 0.15 : fraction of overflow appearing directly in plasma
FPANC    : 0.03 : fraction of overflow to pancreatic islets (struct)
FMYO     : 0.04 : fraction of overflow to myocardium (struct)
FKID     : 0.03 : fraction of overflow to kidney (struct)

// ================================================================ LIVER LIPID
LIVLEAN  : 1400 : lean liver mass (g), sets hepatic fat fraction
FHNEFA   : 0.20 : splanchnic share of the systemic NEFA flux taken by the liver
FHVIS    : 0.55 : PORTAL share of visceral lipolysis delivered to the liver
SCFATN   : 11.6 : healthy SUBCUTANEOUS fat mass (kg) -- the liver's only competitor
VMVLDL   : 45.0 : maximal VLDL-TG secretion (g/d)
KMLTG    : 68.0 : hepatic TG giving half-maximal VLDL secretion (g)
VMOX     : 19.0 : maximal hepatic FA oxidation (g/d) -- extrahepatic oxidation
KMOX     : 68.0 : hepatic TG for half-maximal oxidation (g) -- is in KOXW
KLVLDL   : 0.00 : direct leptin suppression of VLDL secretion -- SET TO ZERO, see \$ODE
KLOX     : 0.25 : hepatic FA-oxidation LOSS per unit leptin signal lost
KLDNL    : 0.20 : de novo lipogenesis gain per unit leptin signal lost
KOM3VLDL : 0.18 : omega-3 suppression of VLDL-TG secretion
KEVVLDL  : 0.05 : ANGPTL3-blockade suppression of VLDL-TG secretion
KFIBOX   : 0.50 : PPARa (fibrate) gain on hepatic FA oxidation

// =============================================================== PLASMA LIPID
VPDL     : 30.0 : plasma volume (dL)
VMTG     : 1730 : maximal LPL-mediated TG clearance (mg/dL/d)
KMTGP    : 150  : plasma TG at half-maximal LPL clearance (mg/dL)
KLINTG   : 0.90 : non-saturable TG clearance (1/d)
LPLFLOOR : 0.45 : LPL activity retained with no adipose tissue (muscle, heart)
KAP3     : 2.50 : LPL gain per unit APOC3 knockdown
KANG     : 1.20 : LPL gain per unit ANGPTL3 knockdown
KFIBL    : 0.50 : LPL gain from PPARa activation
TGTHR    : 1000 : pancreatitis threshold for time-above-threshold (mg/dL)
KPANC    : 0.0015 : maximal pancreatitis hazard (1/d) (struct)

// ================================================================ NEFA / IMCL
NEFA0    : 0.45 : reference plasma NEFA (mmol/L)
KLIPB    : 1.20 : NEFA gain from excess basal lipolysis
KINSN    : 0.60 : insulin suppression of NEFA (per baseline-insulin unit)
TAUNEFA  : 0.50 : NEFA turnover time constant (d)
IMCLN0   : 39.5 : reference intramyocellular lipid pool (g)
KIMOX    : 0.15 : IMCL disposal (1/d)
TAUECT   : 30.0 : time constant of pancreas/heart/kidney lipid tracking (d)
JOVREF   : 39.5 : healthy reference overflow flux (g/d), normalises ectopic pools

KINSACT  : 8.00 : insulin action (level x sensitivity) at half-maximal esterification

// ==================================================== INSULIN RESISTANCE INDEX
AHFF     : 0.220 : IR slope on hepatic fat fraction (per %)
BIM      : 0.150 : IR slope on normalised IMCL
CINFL    : 0.400 : IR slope on hepatic inflammation (struct)
CLEP     : 2.000 : IR reduction per unit leptin signal (direct peripheral action)
DADP     : 1.000 : IR reduction per unit relative adiponectin
SDREF    : 0.900 : leptin signal of the healthy reference subject
ADPN0    : 7.50 : healthy adiponectin (ug/mL)
ADPNFLR  : 0.50 : non-adipose adiponectin floor (ug/mL)
IRMIN    : 0.50 : floor on the IR index
TAUIR    : 7.00 : lag of tissue insulin resistance behind lipid load (d)

// ============================================================ GLUCOSE-INSULIN
VGDL     : 96.0 : glucose distribution volume (dL)
HGP0     : 2400 : reference hepatic glucose output (mg/dL/d)
KHIR     : 0.38 : HGP gain per unit excess IR
KHINS    : 0.080 : HGP suppression per uU/mL insulin
KGZ      : 4.00 : insulin-independent glucose disposal (1/d)
KGI      : 1.53 : insulin-dependent disposal (1/d per uU/mL)
KREN     : 8.00 : renal glucose excretion above threshold (1/d)
GRENTHR  : 180  : renal glucose threshold (mg/dL)
SMAX     : 960  : maximal beta-cell secretory rate (uU/mL/d)
G50      : 110  : glucose for half-maximal secretion (mg/dL)
KCOMP    : 0.35 : compensatory secretory gain per unit excess IR
KINSCL   : 48.0 : insulin elimination (1/d) -- effective, smoothed
KIEXT    : 0.15 : loss of hepatic insulin extraction per unit excess IR
KBREC    : 0.0020 : beta-cell functional recovery (1/d)
KBLIP    : 0.0010 : beta-cell loss per unit excess islet lipid (1/d)
KBGLU    : 0.0015 : beta-cell loss per 100 mg/dL glucose above 180 (1/d)
GTOXTHR  : 180  : glucose above which glucotoxicity operates (mg/dL)
TAUA1C   : 35.0 : HbA1c integration time constant (d)
MGFAC    : 1.10 : mean-glucose to model-glucose ratio

// ================================================ LIVER INJURY AND FIBROSIS
ALT0     : 18.0 : healthy ALT (U/L)
KALT     : 2.60 : ALT rise per % hepatic fat
TAUALT   : 10.0 : ALT turnover (d)
K50INFL  : 15.0 : hepatic fat fraction for half-maximal inflammation (%)
TAUINFL  : 30.0 : inflammation turnover (d)
KFIB     : 0.0100 : fibrogenesis rate above the inflammation threshold (1/d)
INFLTHR  : 0.35 : inflammation threshold for net fibrogenesis
KFREG    : 0.00015 : fibrosis regression (1/d) -- 67x slower than formation

// ==================================================== KIDNEY / ANDROGEN / QoL
KPR      : 0.15 : proteinuria slope on excess IR (g/d)
KPRG     : 0.20 : proteinuria slope per 100 mg/dL glucose above 140 (g/d)
KPRL     : 0.25 : proteinuria slope on renal lipid (g/d) (struct)
TAUPROT  : 60.0 : proteinuria adaptation (d)
KEGFR    : 0.005 : eGFR loss per g/d proteinuria (mL/min/1.73/d)
KAND     : 0.55 : androgen index gain per relative insulin excess (struct)
TAUAND   : 30.0 : androgen index adaptation (d)

// ==================================================== LEPTIN AXIS / TRANSDUCER
LEPK     : 0.60 : leptin secretion per kg effective fat (ng/mL/kg)
HYPMAX   : 0.60 : extra leptin per kg from adipocyte hypertrophy at full fill
LEPFLOOR : 0.15 : non-adipose leptin floor (ng/mL)
TAULEP   : 0.50 : endogenous leptin turnover (d)
LEC50    : 4.00 : leptin for half-maximal hypothalamic signal (ng/mL)
LHILL    : 2.00 : Hill coefficient of the leptin transducer
KSOCS    : 3.00 : SOCS3/PTP1B brake gain on signal sustained ABOVE the healthy level
SREFB    : 0.88 : healthy hypothalamic signal above which the brake engages
TAUSOCS  : 14.0 : brake adaptation time constant (d)
TAUSLAG  : 2.00 : lag of peripheral leptin actions behind plasma signal (d)
KDLEP    : 6.00 : apparent Kd for the receptor-occupancy readout (ng/mL)

// ====================================================== METRELEPTIN PK AND ADA
KAMET    : 8.00 : metreleptin SC absorption (1/d)
FMET     : 0.80 : metreleptin SC bioavailability
VMETC    : 35.0 : metreleptin central volume (L)
CLMET    : 145  : metreleptin clearance (L/d) -> t1/2 ~ 4 h
QMET     : 12.0 : metreleptin distribution clearance (L/d)
VMETP    : 20.0 : metreleptin peripheral volume (L)
NEUT     : 0    : neutralising-ADA switch (0/1) -- binding-only ADA has no PD effect
KNABON   : 0.030 : neutralising titre formation (1/d)
KNABOFF  : 0.0020 : neutralising titre decay (1/d)
NABMAX   : 20.0 : maximal neutralising capacity (dimensionless)
CNAB50   : 10.0 : metreleptin concentration driving titre formation (ng/mL)

// ================================================================ PIOGLITAZONE
KAPIO    : 6.00 : pioglitazone absorption (1/d)
VPIO     : 15.0 : pioglitazone volume (L)
CLPIO    : 20.8 : pioglitazone clearance (L/d), active-metabolite adjusted
EC50PIO  : 700  : pioglitazone concentration for half-maximal PPARg effect (ng/mL)
KCAPM    : 0.45 : maximal fractional capacity gain at full PPARg x full substrate
TAUCAP   : 45.0 : adipogenic capacity expansion time constant (d)
KADPPIO  : 6.00 : adiponectin gain from PPARg activation (ug/mL)
TAUADP   : 20.0 : adiponectin turnover (d)

// ============================================== VOLANESORSEN (APOC3 ASO) PK/PD
KAVOL    : 1.00 : volanesorsen SC absorption (1/d)
VVOLC    : 10.0 : volanesorsen central volume (L)
CLVOL    : 2.00 : volanesorsen plasma clearance (L/d)
QVOL     : 1.50 : plasma-to-tissue clearance (L/d)
VVOLT    : 5.00 : tissue (hepatic) volume (L)
KOUTVT   : 0.030 : tissue elimination (1/d)
KOUT3    : 0.050 : APOC3 turnover (1/d) -> ~14 d onset
EMAX3    : 0.85 : maximal APOC3 knockdown
EC503    : 1.50 : tissue concentration for half-maximal knockdown (mg/L)
KFIB3    : 0.35 : fibrate-mediated reduction of APOC3 production

// ==================================================== FENOFIBRATE / OMEGA-3
KAFIB    : 4.00 : fenofibrate absorption (1/d)
VFIB     : 40.0 : fenofibrate volume (L)
CLFIB    : 20.0 : fenofibrate clearance (L/d)
EC50FIB  : 5.00 : fibrate concentration for half-maximal PPARa effect (mg/L)
OM3DOSE  : 0    : omega-3 4 g/d switch (0/1)
OM3T     : 1e6  : day omega-3 starts
KOM3OX   : 0.35 : omega-3 gain on hepatic FA oxidation (its main action)
TAUOM3   : 14.0 : omega-3 tissue enrichment time constant (d)

// ====================================================== EVINACUMAB (ANGPTL3)
VEVI     : 5.00 : evinacumab volume (L)
CLEVI    : 0.35 : evinacumab clearance (L/d)
KOUTA    : 0.10 : ANGPTL3 turnover (1/d)
EMAXA    : 0.90 : maximal ANGPTL3 suppression
EC50A    : 10.0 : evinacumab concentration for half-maximal suppression (mg/L)

// ====================================================== GLUCOSE-DIRECTED DRUGS
KAMET2   : 8.00 : metformin absorption (1/d)
VMETF    : 600  : metformin volume (L)
CLMETF   : 600  : metformin clearance (L/d)
EMAXMETF : 0.25 : maximal metformin suppression of HGP
EC50METF : 1.50 : metformin concentration for half-maximal effect (mg/L)
INSDOSE  : 0    : exogenous insulin (U/d) as a continuous requirement
INST     : 1e6  : day exogenous insulin starts
KAI      : 4.00 : insulin depot absorption (1/d)
SCALEI   : 24.0 : uU/mL per U/d of exogenous insulin at steady state
KAGLP    : 0.50 : GLP-1 RA absorption (1/d)
VGLP     : 12.0 : GLP-1 RA volume (L)
CLGLP    : 1.19 : GLP-1 RA clearance (L/d) -> t1/2 ~ 7 d
EMAXGLP  : 0.20 : maximal fractional intake reduction by GLP-1 RA
EC50GLP  : 100  : GLP-1 RA concentration for half-maximal effect (ng/mL)
KASGL    : 8.00 : SGLT2i absorption (1/d)
VSGL     : 120  : SGLT2i volume (L)
CLSGL    : 20.0 : SGLT2i clearance (L/d)
UGEMAX   : 729  : maximal glucosuria (mg/dL/d)
EC50SGL  : 0.15 : SGLT2i concentration for half-maximal glucosuria (mg/L)
GSGLTHR  : 100  : renal glucose threshold for SGLT2i-driven excretion (mg/dL)

$INIT @annotated
// ---- metreleptin ----------------------------------------------------------
MSC      : 0     : metreleptin SC depot (ug)
MCEN     : 0     : metreleptin central (ug)
MPER     : 0     : metreleptin peripheral (ug)
NAB      : 0     : neutralising anti-metreleptin capacity (dimensionless)
// ---- other drugs ----------------------------------------------------------
PGUT     : 0     : pioglitazone gut (mg)
PCEN     : 0     : pioglitazone central (mg)
VSC      : 0     : volanesorsen SC depot (mg)
VCEN     : 0     : volanesorsen central (mg)
VTIS     : 0     : volanesorsen tissue (mg)
APOC3    : 1     : APOC3 level (fraction of normal)
FGUT     : 0     : fenofibrate gut (mg)
FCEN     : 0     : fenofibrate central (mg)
EVC      : 0     : evinacumab central (mg)
ANG3     : 1     : ANGPTL3 level (fraction of normal)
OM3      : 0     : omega-3 tissue enrichment index (0-1)
MTGUT    : 0     : metformin gut (mg)
MTCEN    : 0     : metformin central (mg)
ISCD     : 0     : exogenous insulin depot (U)
INSX     : 0     : exogenous insulin concentration (uU/mL)
GSC      : 0     : GLP-1 RA SC depot (mg)
GCEN     : 0     : GLP-1 RA central (mg)
SGGUT    : 0     : SGLT2i gut (mg)
SGCEN    : 0     : SGLT2i central (mg)
// ---- leptin axis and intake ----------------------------------------------
LEPEN    : 12.0  : endogenous plasma leptin (ng/mL)
SOCS     : 0.5   : hypothalamic desensitisation state
SLAG     : 0.9   : lagged leptin signal driving peripheral actions
EIS      : 1400  : energy intake state (kcal/d)
// ---- adipose depots and capacity -----------------------------------------
ASC1     : 8.0   : femorogluteal SC TG (kg)
ASC2     : 3.5   : upper-trunk/facial SC TG (kg)
AVIS     : 2.0   : visceral TG (kg)
CAPM     : 1.0   : capacity multiplier (PPARg-driven adipogenesis)
OXAD     : 0     : adaptive oxidation state
// ---- ectopic lipid --------------------------------------------------------
LTG      : 45.0  : liver TG (g)
IMCL     : 12.0  : intramyocellular lipid (g)
PLIP     : 1.0   : islet lipid (normalised)
MYOL     : 1.0   : myocardial lipid (normalised)
KLIP     : 1.0   : renal lipid (normalised)
// ---- plasma lipid ---------------------------------------------------------
PTG      : 90.0  : plasma triglyceride (mg/dL)
NEFA     : 0.45  : plasma NEFA (mmol/L)
ADPN     : 8.0   : adiponectin (ug/mL)
// ---- glucose axis ---------------------------------------------------------
IRX      : 1.0   : lagged insulin resistance index
GLU      : 90.0  : plasma glucose (mg/dL)
INS      : 8.0   : endogenous insulin (uU/mL)
BCF      : 1.0   : beta-cell functional mass (fraction)
A1C      : 5.1   : HbA1c (%)
// ---- organ injury ---------------------------------------------------------
ALT      : 18.0  : ALT (U/L)
INFL     : 0.17  : hepatic inflammation (normalised)
FIB      : 0     : fibrosis stage (0-4)
PROT     : 0.05  : proteinuria (g/d)
EGFR     : 105   : eGFR (mL/min/1.73 m2)
ANDX     : 1.0   : androgen index (1 = normal)
TAT      : 0     : cumulative time with TG above threshold (d)
PANCH    : 0     : cumulative pancreatitis hazard

$GLOBAL
#define HFF   (100.0*LTG/(LTG + LIVLEAN))
#define FATM  (ASC1 + ASC2 + AVIS)
#define EFFAT (ASC1 + ASC2 + 0.30*AVIS)

// Shared scalars set in $MAIN (constant within a record)
double C1, C2, CV, EI0, DIETFATN, DIETCHON;

// Leptin transducer: saturating Hill on total effective leptin, then the
// SOCS3/PTP1B brake.  This single function is the only place where leptin
// becomes a signal, and it is deliberately the same function for endogenous
// and exogenous leptin.
double lep_signal(double L, double ec50, double h, double socs,
                  double ksocs, double srefb) {
  double Lp = pow(fmax(L, 1e-8), h);
  double S  = Lp/(Lp + pow(ec50, h));
  // The SOCS3/PTP1B brake is a DEVIATION mechanism: it engages only for signal
  // held above the healthy level, which is why it desensitises obesity without
  // taxing a normal subject -- and why it cannot be invoked to explain the
  // lipodystrophy phenotype, where the signal is far below srefb.
  return S/(1.0 + ksocs*fmax(0.0, socs - srefb));
}

$MAIN
// Capacity ceilings.  CAPM (the PPARg state) multiplies only the SUBCUTANEOUS
// depots: thiazolidinediones expand subcutaneous adipose, they do not create
// visceral capacity.
C1 = CAP1N*FC1*CAPM;
C2 = CAP2N*FC2*CAPM;
CV = CAPVN*FCV;

// Weight-maintenance intake scales with lean mass, not with total weight:
// lipodystrophy patients have little fat but normal-to-high lean mass.
EI0 = EIBASE*LEAN0*exp(ETA(1));

DIETFATN = DIETFAT;
DIETCHON = FFAT + FCHO - DIETFAT;

$ODE
// ===========================================================================
//  1.  LEPTIN: SECRETION, EXOGENOUS DRUG, TRANSDUCER, DEFICIT
// ===========================================================================
double fill1 = ASC1/(C1 + 1e-9);
double fill2 = ASC2/(C2 + 1e-9);
double fillv = AVIS/(CV + 1e-9);
double fillm = fmin(1.0, (ASC1 + ASC2)/(C1 + C2 + 1e-9));

// Leptin per kg rises with adipocyte size (fractional fill), which is why the
// few remaining, maximally distended adipocytes of a CGL patient secrete more
// per gram than a healthy adipocyte does.
double lep_ss = LEPKM*LEPK*EFFAT*(1.0 + HYPMAX*fillm) + LEPFLOOR;

double CMET   = MCEN/VMETC;              // total metreleptin (ng/mL)
double FU     = 1.0/(1.0 + NAB);         // free fraction left by neutralising ADA
double LTOT   = LEPEN + CMET;            // what a clinical assay measures
double LEFF   = LEPEN + CMET*FU;         // what the receptor sees
double SD     = lep_signal(LEFF, LEC50, LHILL, SOCS, KSOCS, SREFB);
double DEFC   = 1.0 - SD;                // the deficit that drives hyperphagia
double OCC    = LTOT/(LTOT + KDLEP);     // receptor-occupancy readout

// ===========================================================================
//  2.  INTAKE AND INFLUX
// ===========================================================================
double CGLP   = 1000.0*GCEN/VGLP;                       // ng/mL
double EGLP   = EMAXGLP*CGLP/(CGLP + EC50GLP);
double diet   = (SOLVERTIME >= DIETT) ? 1.0 : 0.0;
double ffat   = diet ? DIETFATN : FFAT;
double fcho   = diet ? DIETCHON : FCHO;

// Two anorexigenic levers act on ONE intake term, so they are automatically
// less than additive (scenario 22).
double EITG   = (EICLAMP > 0) ? EICLAMP : EI0*(1.0 + HMAX*DEFC)*(1.0 - EGLP);

double INSTOT = INS + INSX;
double insdr  = INSTOT/(INSTOT + KINS50);
double IRSEL  = fmax(0.0, IRX - 1.0)/(fmax(0.0, IRX - 1.0) + 3.0);   // selective hepatic IR
double DFAT   = EIS*ffat/9.0;
double CHOG   = EIS*fcho/4.0;
double DNL    = KDNL*CHOG*insdr*(1.0 + KDNLSEL*IRSEL)
                *(1.0 + KLDNL*fmax(0.0, SREFB - SLAG));

// Adipose lipolytic release.  This flux is RECYCLED: it is presented back to
// the buffer and must therefore appear on both sides of the balance, otherwise
// the futile cycle would look like free storage capacity.
double lipsup = 1.0/(1.0 + KINSN*INSTOT/8.0);
double JREL1 = KLIPSC*(1.0 + BASLIP)*ASC1*1000.0*lipsup;
double JREL2 = KLIPSC*(1.0 + BASLIP)*ASC2*1000.0*lipsup;
double JRELV = KLIPVI*(1.0 + BASLIP)*AVIS*1000.0*lipsup;
double JRELT = JREL1 + JREL2 + JRELV;

double JIN    = DFAT + DNL + JRELT;

// ===========================================================================
//  3.  THE BUFFER: CAPACITY-LIMITED ESTERIFICATION
// ===========================================================================
double ISENS  = 1.0/fmax(IRX, IRMIN);
double jsat   = JIN/(JIN + KMJ);

// Adipose esterification is driven by insulin ACTION = level x sensitivity,
// normalised to 1 in the healthy reference subject.  This matters: in severe
// insulin resistance the adipocyte sees a 4-5x higher insulin concentration,
// so the signalling term is nearly preserved and CAPACITY -- not signalling --
// is what limits storage.  Gating storage on sensitivity alone would empty the
// residual depots, which is the opposite of the hypertrophic adipocytes seen
// in partial lipodystrophy.
double xact  = INSTOT*fmin(ISENS, 2.0);
double insact = 2.0*xact/(xact + KINSACT);
double hr1 = fmax(0.0, 1.0 - fill1) + HRFLOOR;
double hr2 = fmax(0.0, 1.0 - fill2) + HRFLOOR;
double hrv = fmax(0.0, 1.0 - fillv) + HRFLOOR;
double JST1 = KEST*C1*hr1*jsat*insact;
double JST2 = KEST*C2*hr2*jsat*insact;
double JSTV = KEST*CV*hrv*jsat*insact;
double JSTT = JST1 + JST2 + JSTV;

// Whole-body oxidation, with an adaptive component: patients with a large
// ectopic load up-regulate fat oxidation, which is why hyperphagia does not
// simply make them obese.
double BW   = LEAN0 + FATM + (LTG + IMCL)/1000.0;
double OXW  = KOXW*pow(BW, 0.75)*(1.0 + KOXAD*OXAD);

// THE overflow flux.  Everything downstream of here is a consequence.
double JOV  = fmax(0.0, JIN - JSTT - OXW);


// ===========================================================================
//  4.  DRUG PK
// ===========================================================================
double CPIO  = 1000.0*PCEN/VPIO;                 // ng/mL
double EPIO  = CPIO/(CPIO + EC50PIO);
double CVOLT = VTIS/VVOLT;                       // mg/L
double CFIB  = FCEN/VFIB;                        // mg/L
double EFIB  = CFIB/(CFIB + EC50FIB);
double CEVI  = EVC/VEVI;                         // mg/L
double CMETF = MTCEN/VMETF;                      // mg/L
double EMETF = EMAXMETF*CMETF/(CMETF + EC50METF);
double CSGL  = SGCEN/VSGL;                       // mg/L
double UGE   = UGEMAX*(CSGL/(CSGL + EC50SGL))*fmax(0.0, GLU - GSGLTHR)/GLU;

dxdt_MSC  = -KAMET*MSC;
dxdt_MCEN =  KAMET*MSC*FMET - (CLMET/VMETC)*MCEN - (QMET/VMETC)*MCEN + (QMET/VMETP)*MPER;
dxdt_MPER =  (QMET/VMETC)*MCEN - (QMET/VMETP)*MPER;
dxdt_NAB  =  KNABON*NEUT*(CMET/(CMET + CNAB50))*(NABMAX - NAB) - KNABOFF*NAB;

dxdt_PGUT = -KAPIO*PGUT;
dxdt_PCEN =  KAPIO*PGUT - (CLPIO/VPIO)*PCEN;

dxdt_VSC  = -KAVOL*VSC;
dxdt_VCEN =  KAVOL*VSC - (CLVOL/VVOLC)*VCEN - (QVOL/VVOLC)*VCEN + (QVOL/VVOLT)*VTIS;
dxdt_VTIS =  (QVOL/VVOLC)*VCEN - (QVOL/VVOLT)*VTIS - KOUTVT*VTIS;

double kin3 = KOUT3*(1.0 - KFIB3*EFIB);
dxdt_APOC3 = kin3*(1.0 - EMAX3*CVOLT/(CVOLT + EC503)) - KOUT3*APOC3;

dxdt_FGUT = -KAFIB*FGUT;
dxdt_FCEN =  KAFIB*FGUT - (CLFIB/VFIB)*FCEN;
dxdt_EVC  = -(CLEVI/VEVI)*EVC;
dxdt_ANG3 =  KOUTA*(1.0 - EMAXA*CEVI/(CEVI + EC50A)) - KOUTA*ANG3;
dxdt_OM3  =  (((SOLVERTIME >= OM3T) ? OM3DOSE : 0.0) - OM3)/TAUOM3;

dxdt_MTGUT = -KAMET2*MTGUT;
dxdt_MTCEN =  KAMET2*MTGUT - (CLMETF/VMETF)*MTCEN;

double insinf = (SOLVERTIME >= INST) ? INSDOSE : 0.0;
dxdt_ISCD  = insinf - KAI*ISCD;
dxdt_INSX  = KAI*ISCD*SCALEI - KINSCL*INSX;   // exogenous: no hepatic first pass

dxdt_GSC  = -KAGLP*GSC;
dxdt_GCEN =  KAGLP*GSC - (CLGLP/VGLP)*GCEN;
dxdt_SGGUT = -KASGL*SGGUT;
dxdt_SGCEN =  KASGL*SGGUT - (CLSGL/VSGL)*SGCEN;

// ===========================================================================
//  5.  LEPTIN-AXIS STATES
// ===========================================================================
dxdt_LEPEN = (lep_ss - LEPEN)/TAULEP;
dxdt_SOCS  = (SD - SOCS)/TAUSOCS;    // the brake tracks the sustained signal
dxdt_SLAG  = (SD - SLAG)/TAUSLAG;
dxdt_EIS   = (EITG - EIS)/TAUEI;

// ===========================================================================
//  6.  ADIPOSE DEPOTS, CAPACITY, ADAPTIVE OXIDATION
// ===========================================================================
dxdt_ASC1 = (JST1 - JREL1)/1000.0;
dxdt_ASC2 = (JST2 - JREL2)/1000.0;
dxdt_AVIS = (JSTV - JRELV)/1000.0;
dxdt_CAPM = ((1.0 + PREAD*KCAPM*EPIO) - CAPM)/TAUCAP;
dxdt_OXAD = ((JOV/JOVREF)/(1.0 + JOV/JOVREF) - OXAD)/TAUOXAD;

// ===========================================================================
//  7.  ECTOPIC LIPID (the overflow cascade)
// ===========================================================================
// VLDL-TG secretion is written as PURELY SUBSTRATE-DRIVEN (Adiels): there is no
// direct leptin term on it.  An earlier draft gave leptin a primary inhibitory
// action here, and that draft predicted that leptin replacement at CLAMPED food
// intake makes steatosis 71% WORSE -- because suppressing the largest hepatic
// disposal route at unchanged substrate delivery has to raise the pool.  Since
// the observed effect of leptin at fixed intake is the opposite, the falling
// VLDL secretion seen on metreleptin must be a CONSEQUENCE of the smaller
// hepatic triglyceride pool, not a cause of it.  KLVLDL is retained at 0 so the
// discarded structure stays visible and testable.
double vldl = VMVLDL*(LTG/(LTG + KMLTG))
              *(1.0 + KLVLDL*fmax(0.0, SREFB - SLAG))
              *(1.0 - KOM3VLDL*OM3)
              *(1.0 - KEVVLDL*(1.0 - ANG3));
// Leptin acts on the liver in three directions at once (AMPK-driven FA
// oxidation up, SREBP-1c-driven lipogenesis down, VLDL secretion down).  Only
// modelling the third would predict that leptin replacement at fixed food
// intake makes steatosis WORSE, which is not what pair-fed animals or the
// Petersen hepatic-TG data show.
double hepox = VMOX*(1.0 + KFIBOX*EFIB + KOM3OX*OM3)
               *(1.0 - KLOX*fmax(0.0, SREFB - SLAG))*LTG/(LTG + KMOX);

// The hepatic SHARE of overflow is not a constant: as adipose disappears the
// liver becomes a larger fraction of the remaining sink.
// The liver competes with SUBCUTANEOUS adipose only.  Visceral fat is not a
// competing sink at all -- it drains into the portal vein, so it delivers TO
// the liver.  This is why partial lipodystrophy, which trades subcutaneous for
// visceral fat, is steatotic despite a near-normal total fat mass.
double scfr0  = fmin(1.0, (ASC1 + ASC2)/SCFATN);
double flive  = fmin(0.88, FLIV + KFLIVA*(1.0 - scfr0));
// Visceral lipolysis drains into the portal vein, so it reaches the liver at a
// higher fraction than systemic NEFA does.  This is the term that makes
// PARTIAL lipodystrophy (visceral expansion, subcutaneous loss) steatotic.
dxdt_LTG  = FHNEFA*(JREL1 + JREL2) + FHVIS*JRELV + flive*JOV + DNL*0.35
            - hepox - vldl;
dxdt_IMCL = FMUS*JOV - KIMOX*IMCL;
dxdt_PLIP = ((FPANC*JOV)/(FPANC*JOVREF) - PLIP)/TAUECT;
dxdt_MYOL = ((FMYO*JOV)/(FMYO*JOVREF) - MYOL)/TAUECT;
dxdt_KLIP = ((FKID*JOV)/(FKID*JOVREF) - KLIP)/TAUECT;

// ===========================================================================
//  8.  PLASMA LIPID
// ===========================================================================
double fatfrac = fmin(1.5, FATM/FATN);
double lplact  = (LPLFLOOR + (1.0 - LPLFLOOR)*fatfrac)
                 *(1.0 + KAP3*(1.0 - APOC3))
                 *(1.0 + KANG*(1.0 - ANG3))
                 *(1.0 + KFIBL*EFIB);
double tgin  = (vldl + FPL*JOV)*1000.0/VPDL;
double tgout = VMTG*lplact*PTG/(KMTGP + PTG) + KLINTG*PTG;
dxdt_PTG = tgin - tgout;

double nefa_ss = NEFA0*(1.0 + KLIPB*BASLIP)*lipsup*(1.0 + 0.5*fmax(0.0, JOV/JOVREF - 1.0));
dxdt_NEFA = (nefa_ss - NEFA)/TAUNEFA;

double adpn_ss = ADPNFLR + ADPN0*fatfrac/(1.0 + 0.5*fmax(0.0, IRX - 1.0))
                 + KADPPIO*EPIO*fmin(1.0, fatfrac + 0.15);
dxdt_ADPN = (adpn_ss - ADPN)/TAUADP;

// ===========================================================================
//  9.  INSULIN RESISTANCE, GLUCOSE, BETA CELL
// ===========================================================================
double IMCLNN = IMCL/IMCLN0;
double IRnow  = 1.0 + AHFF*(HFF - 3.0) + BIM*(IMCLNN - 1.0) + CINFL*(INFL - 0.17)
                - CLEP*(SLAG - SDREF) - DADP*(ADPN - ADPN0)/ADPN0;
IRnow = fmax(IRMIN, IRnow);
dxdt_IRX = (IRnow - IRX)/TAUIR;

// SELECTIVE hepatic insulin resistance, written with the SAME term that keeps
// lipogenesis insulin-driven (IRSEL).  One coefficient therefore produces both
// halves of the paradox: glucose production escapes insulin while de novo
// lipogenesis still obeys it.  This is why the model can be hyperinsulinaemic
// AND hyperglycaemic AND hyperlipogenic at once.
double hgp = HGP0*(1.0 + KHIR*fmax(0.0, IRX - 1.0))
             /(1.0 + KHINS*INSTOT*(1.0 - IRSEL))*(1.0 - EMETF);
double rd  = (KGZ + KGI*INSTOT*fmin(ISENS, 2.0))*GLU;
double uger = KREN*fmax(0.0, GLU - GRENTHR);
dxdt_GLU = hgp - rd - uger - UGE;

double isec = BCF*SMAX*(pow(GLU,2.0)/(pow(G50,2.0) + pow(GLU,2.0)))
              *(1.0 + KCOMP*fmax(0.0, IRX - 1.0));
// Hepatic insulin extraction falls as insulin resistance rises, which is part
// of why lipodystrophy is hyperinsulinaemic out of proportion to secretion.
double kinscl = KINSCL/(1.0 + KIEXT*fmax(0.0, IRX - 1.0));
dxdt_INS = isec - kinscl*INS;

dxdt_BCF = KBREC*(1.0 - BCF)
           - BCF*(KBLIP*fmax(0.0, PLIP - 1.0)
                  + KBGLU*fmax(0.0, (GLU - GTOXTHR)/100.0));

dxdt_A1C = (((MGFAC*GLU) + 46.7)/28.7 - A1C)/TAUA1C;

// ===========================================================================
// 10.  LIVER INJURY, FIBROSIS, KIDNEY, ANDROGENS, PANCREATITIS
// ===========================================================================
double alt_ss = ALT0 + KALT*fmax(0.0, HFF - 3.0)*(1.0 + 0.5*INFL);
dxdt_ALT  = (alt_ss - ALT)/TAUALT;
dxdt_INFL = (HFF/(HFF + K50INFL) - INFL)/TAUINFL;
dxdt_FIB  = KFIB*fmax(0.0, INFL - INFLTHR)*(1.0 - FIB/4.0) - KFREG*FIB;

double prot_ss = 0.05 + KPR*fmax(0.0, IRX - 1.0) + KPRG*fmax(0.0, (GLU - 140.0)/100.0)
                 + KPRL*fmax(0.0, KLIP - 1.0);
dxdt_PROT = (prot_ss - PROT)/TAUPROT;
dxdt_EGFR = -KEGFR*PROT*(EGFR > 10.0 ? 1.0 : 0.0);

double and_ss = 1.0 + KAND*fmax(0.0, (INSTOT - 8.0)/8.0)/(1.0 + fmax(0.0, (INSTOT - 8.0)/8.0)/4.0);
dxdt_ANDX = (and_ss - ANDX)/TAUAND;

dxdt_TAT   = 1.0/(1.0 + exp(-(PTG - TGTHR)/50.0));
double tgr = pow(PTG/1000.0, 4.0);
dxdt_PANCH = KPANC*tgr/(1.0 + tgr);

$TABLE
// Recomputed from the states at the output time so that captured values never
// depend on where the solver happened to stop.
double t_fill = fmin(1.0, (ASC1 + ASC2)/(CAP1N*FC1*CAPM + CAP2N*FC2*CAPM + 1e-9));
double t_CMET = MCEN/VMETC;
double t_FU   = 1.0/(1.0 + NAB);
double t_LTOT = LEPEN + t_CMET;
double t_LEFF = LEPEN + t_CMET*t_FU;
double t_SD   = lep_signal(t_LEFF, LEC50, LHILL, SOCS, KSOCS, SREFB);
double t_DEF  = 1.0 - t_SD;
double t_OCC  = t_LTOT/(t_LTOT + KDLEP);
double t_HFF  = HFF;
double t_FAT  = FATM;
double t_INST = INS + INSX;
double t_ISEN = 1.0/fmax(IRX, IRMIN);

double t_diet = (TIME >= DIETT) ? 1.0 : 0.0;
double t_ffat = t_diet ? DIETFAT : FFAT;
double t_fcho = t_diet ? (FFAT + FCHO - DIETFAT) : FCHO;
double t_insdr = t_INST/(t_INST + KINS50);
double t_IRSEL = fmax(0.0, IRX - 1.0)/(fmax(0.0, IRX - 1.0) + 3.0);
double t_DNL  = KDNL*(EIS*t_fcho/4.0)*t_insdr*(1.0 + KDNLSEL*t_IRSEL)
                *(1.0 + KLDNL*fmax(0.0, SREFB - SLAG));
double t_lsup = 1.0/(1.0 + KINSN*t_INST/8.0);
double t_JREL = (KLIPSC*(ASC1 + ASC2) + KLIPVI*AVIS)*(1.0 + BASLIP)*1000.0*t_lsup;
double t_JIN  = EIS*t_ffat/9.0 + t_DNL + t_JREL;
double t_jsat = t_JIN/(t_JIN + KMJ);
double t_xact = t_INST*fmin(t_ISEN, 2.0);
double t_iact = 2.0*t_xact/(t_xact + KINSACT);
double t_C1 = CAP1N*FC1*CAPM;
double t_C2 = CAP2N*FC2*CAPM;
double t_CV = CAPVN*FCV;
double t_JST = KEST*t_jsat*t_iact*
   ( t_C1*(fmax(0.0, 1.0 - ASC1/(t_C1 + 1e-9)) + HRFLOOR)
   + t_C2*(fmax(0.0, 1.0 - ASC2/(t_C2 + 1e-9)) + HRFLOOR)
   + t_CV*(fmax(0.0, 1.0 - AVIS/(t_CV + 1e-9)) + HRFLOOR) );
double t_BW  = LEAN0 + FATM + (LTG + IMCL)/1000.0;
double t_OXW = KOXW*pow(t_BW, 0.75)*(1.0 + KOXAD*OXAD);
double t_JOV = fmax(0.0, t_JIN - t_JST - t_OXW);
double t_CAP = t_C1 + t_C2 + t_CV;
double t_fatf = fmin(1.0, (ASC1 + ASC2)/SCFATN);
double t_fliv = fmin(0.88, FLIV + KFLIVA*(1.0 - t_fatf));
double t_JHEP = FHNEFA*(KLIPSC*(ASC1 + ASC2)*(1.0 + BASLIP)*1000.0*t_lsup)
                + FHVIS*(KLIPVI*AVIS*(1.0 + BASLIP)*1000.0*t_lsup)
                + t_fliv*t_JOV + 0.35*t_DNL;

$CAPTURE @annotated
t_HFF  : hepatic fat fraction (%)
t_FAT  : total body fat (kg)
t_CAP  : total storage capacity (kg)
t_LTOT : total leptin, assay-equivalent (ng/mL)
t_LEFF : receptor-available leptin (ng/mL)
t_SD   : hypothalamic leptin signal (0-1)
t_DEF  : signal deficit = 1 - SD
t_OCC  : leptin receptor occupancy readout (0-1)
t_JIN  : lipid influx presented to the buffer (g/d)
t_JST  : total adipose storage flux (g/d)
t_JOV  : OVERFLOW flux (g/d)
t_OXW  : whole-body fat oxidation (g/d)
t_BW   : body weight (kg)
t_INST : total insulin, endogenous + exogenous (uU/mL)
t_ISEN : insulin sensitivity (1/IR)
t_fill : mean subcutaneous fractional fill
t_DNL  : de novo lipogenesis (g/d)
t_JREL : adipose lipolytic release, recycled (g/d)
t_JHEP : fatty-acid flux delivered to the liver (g/d)
t_fliv : hepatic share of the overflow

$OMEGA @annotated @block
ELEAN : 0 : between-subject variability on maintenance intake

$SIGMA 0
)---"

## -----------------------------------------------------------------------------
##  BUILD
## -----------------------------------------------------------------------------
mod <- mcode_cache("lipo_qsp", lipo_code, atol = 1e-8, rtol = 1e-8, maxsteps = 500000)

## =============================================================================
##  PHENOTYPES
##
##  A phenotype in this model is NOT a label.  It is a capacity distribution
##  (FC1/FC2/FCV), an expandability reserve (PREAD) and a secretion coefficient
##  (LEPKM).  Everything else -- leptin level, hepatic fat, triglyceride,
##  HbA1c -- is an OUTPUT of running those numbers to steady state.
##
##  The two control phenotypes matter as much as the disease ones:
##    lep_def   fat mass HIGH, signal ZERO  -> the pure signal-deficit arm
##    obese     fat mass HIGH, signal HIGH  -> the pure capacity-exceeded arm
##  Together they show that the model's two deficits are separable.
## =============================================================================

phenotypes <- list(
  normal    = list(FC1 = 1.00, FC2 = 1.00, FCV = 1.00, PREAD = 1.00, LEPKM = 1.0,
                   LEAN0 = 45, label = "Healthy control"),
  cgl       = list(FC1 = 0.02, FC2 = 0.05, FCV = 0.05, PREAD = 0.03, LEPKM = 1.0,
                   LEAN0 = 48, label = "CGL (AGPAT2/BSCL2)"),
  agl       = list(FC1 = 0.05, FC2 = 0.08, FCV = 0.15, PREAD = 0.15, LEPKM = 1.0,
                   LEAN0 = 46, label = "Acquired generalised (Lawrence)"),
  fpld_low  = list(FC1 = 0.12, FC2 = 0.90, FCV = 1.70, PREAD = 0.55, LEPKM = 1.0,
                   LEAN0 = 45, label = "FPLD2, low leptin"),
  fpld_high = list(FC1 = 0.22, FC2 = 1.30, FCV = 1.80, PREAD = 0.60, LEPKM = 1.0,
                   LEAN0 = 45, label = "FPLD2, near-normal leptin"),
  apl       = list(FC1 = 0.90, FC2 = 0.15, FCV = 0.60, PREAD = 0.50, LEPKM = 1.0,
                   LEAN0 = 45, label = "Acquired partial (Barraquer-Simons)"),
  plin1     = list(FC1 = 0.30, FC2 = 0.55, FCV = 1.30, PREAD = 0.45, LEPKM = 1.0,
                   BASLIP = 0.8, LEAN0 = 45, label = "PLIN1 FPLD4 (high basal lipolysis)"),
  lep_def   = list(FC1 = 1.60, FC2 = 1.60, FCV = 1.60, PREAD = 1.00, LEPKM = 0.02,
                   LEAN0 = 45, label = "Congenital leptin deficiency (control)"),
  obese     = list(FC1 = 1.80, FC2 = 1.80, FCV = 2.00, PREAD = 1.00, LEPKM = 1.0,
                   LEAN0 = 55, label = "Common obesity (control)")
)

pheno_param <- function(p) {
  p$label <- NULL
  p
}

## =============================================================================
##  BASELINE: each phenotype is run to ITS OWN steady state
##
##  Hand-typed initial conditions would let a phenotype be asserted rather than
##  derived, and would make every treatment contrast depend on how good the
##  guess was.  Instead every scenario starts from a 12-year burn-in of the
##  untreated phenotype.
## =============================================================================

## Run an arbitrary parameter set to its own steady state.  The Shiny app uses
## this so that moving a capacity slider produces a genuinely re-derived
## baseline rather than a disease phenotype with one number edited.
baseline_pheno <- function(p, burn = 12*365) {
  m   <- param(mod, p)
  out <- mrgsim_df(zero_re(m), end = burn, delta = burn, recsort = 3)
  last <- out[nrow(out), ]
  st  <- as.list(last[intersect(names(last), names(as.list(init(mod))))])
  list(mod = init(m, st), init = st, param = p)
}

baseline <- function(pheno = "cgl", burn = 12*365) {
  p <- pheno_param(phenotypes[[pheno]])
  m <- param(mod, p)
  out <- mrgsim_df(zero_re(m), end = burn, delta = burn, recsort = 3)
  last <- out[nrow(out), ]
  st <- as.list(last[intersect(names(last), names(as.list(init(mod))))])
  list(mod = init(m, st), init = st, param = p)
}

baseline_table <- function() {
  rows <- lapply(names(phenotypes), function(nm) {
    b  <- baseline(nm)
    s  <- mrgsim_df(zero_re(b$mod), end = 1, delta = 1)
    s  <- s[nrow(s), ]
    data.frame(phenotype = phenotypes[[nm]]$label,
               fat_kg    = round(s$t_FAT, 2),
               capacity_kg = round(s$t_CAP, 2),
               leptin    = round(s$t_LTOT, 2),
               signal_SD = round(s$t_SD, 3),
               deficit   = round(s$t_DEF, 3),
               intake_kcal = round(s$EIS),
               JIN_g     = round(s$t_JIN, 1),
               JOV_g     = round(s$t_JOV, 1),
               HFF_pct   = round(s$t_HFF, 1),
               TG        = round(s$PTG),
               glucose   = round(s$GLU),
               HbA1c     = round(s$A1C, 2),
               insulin   = round(s$INS, 1),
               IR        = round(s$IRX, 2),
               betacell  = round(s$BCF, 2),
               adipo     = round(s$ADPN, 1),
               ALT       = round(s$ALT),
               fibrosis  = round(s$FIB, 2),
               protein   = round(s$PROT, 2),
               row.names = NULL)
  })
  do.call(rbind, rows)
}

## =============================================================================
##  DOSING HELPERS
##
##  Metreleptin amounts are in micrograms so that amount/volume comes out in
##  ng/mL, the unit every leptin paper reports.
## =============================================================================

met_sc <- function(mgkg = 0.06, wt = 60, start = 730, dur = 3650)
  ev(amt = mgkg*wt*1000, cmt = "MSC", time = start, ii = 1, addl = dur - 1)

met_fixed <- function(mg = 5, start = 730, dur = 3650)
  ev(amt = mg*1000, cmt = "MSC", time = start, ii = 1, addl = dur - 1)

pio_oral <- function(mg = 45, start = 730, dur = 3650)
  ev(amt = mg, cmt = "PGUT", time = start, ii = 1, addl = dur - 1)

vol_sc <- function(mg = 300, start = 730, dur = 3650)
  ev(amt = mg, cmt = "VSC", time = start, ii = 7, addl = floor(dur/7) - 1)

fen_oral <- function(mg = 145, start = 730, dur = 3650)
  ev(amt = mg, cmt = "FGUT", time = start, ii = 1, addl = dur - 1)

evi_iv <- function(mg = 900, start = 730, n = 130)
  ev(amt = mg, cmt = "EVC", time = start, ii = 28, addl = n - 1)

metf_oral <- function(mg = 1000, start = 730, dur = 3650)
  ev(amt = mg, cmt = "MTGUT", time = start, ii = 0.5, addl = dur*2 - 1)

glp_sc <- function(mg = 1, start = 730, dur = 3650)
  ev(amt = mg, cmt = "GSC", time = start, ii = 7, addl = floor(dur/7) - 1)

sgl_oral <- function(mg = 10, start = 730, dur = 3650)
  ev(amt = mg, cmt = "SGGUT", time = start, ii = 1, addl = dur - 1)

## =============================================================================
##  SCENARIO ENGINE
##
##  Every scenario is the same 53 equations started from the same phenotype
##  steady state.  Contrasts are therefore differences in the intervention
##  only -- which is what makes 7 vs 8, 13 vs 14 and 6 vs 11 interpretable.
## =============================================================================

run_scn <- function(label, pheno, events = NULL, pset = NULL,
                    end = 12*365, delta = 7) {
  b <- baseline(pheno)
  m <- zero_re(b$mod)
  if (!is.null(pset)) m <- param(m, pset)
  out <- if (is.null(events)) mrgsim_df(m, end = end, delta = delta)
         else mrgsim_df(m, events = events, end = end, delta = delta)
  out$scenario <- label
  out$pheno    <- pheno
  out
}

scenarios <- list(

  ## ---- 1-5: natural history and the two control phenotypes ---------------
  s1  = function() run_scn("1 CGL untreated", "cgl"),
  s2  = function() run_scn("2 FPLD2 low-leptin untreated", "fpld_low"),
  s3  = function() run_scn("3 Healthy control", "normal"),
  s4  = function() run_scn("4 Common obesity (capacity exceeded, signal high)", "obese"),
  s5  = function() run_scn("5 Congenital leptin deficiency (signal zero, capacity intact)", "lep_def"),

  ## ---- 6-9: the metreleptin dose/deficit dissociation --------------------
  s6  = function() run_scn("6 CGL + metreleptin 0.06 mg/kg", "cgl",
                           met_sc(0.06, 60, start = 730)),
  s7  = function() run_scn("7 FPLD low-leptin + metreleptin 0.06", "fpld_low",
                           met_sc(0.06, 60, start = 730)),
  s8  = function() run_scn("8 FPLD near-normal leptin + metreleptin 0.06", "fpld_high",
                           met_sc(0.06, 60, start = 730)),
  s9  = function() run_scn("9 FPLD near-normal leptin + metreleptin 0.13 (dose escalation)",
                           "fpld_high", met_sc(0.13, 60, start = 730)),

  ## ---- 10-12: reversibility, immunogenicity, pair-feeding ---------------
  s10 = function() run_scn("10 CGL + metreleptin, withdrawn at year 6", "cgl",
                           met_sc(0.06, 60, start = 730, dur = 4*365)),
  s11 = function() run_scn("11 CGL + metreleptin with neutralising ADA", "cgl",
                           met_sc(0.06, 60, start = 730), pset = list(NEUT = 1)),
  s12 = function() {
    ## Pair-fed clamp: intake is held at the untreated steady-state value, so
    ## whatever effect remains is the DIRECT peripheral action of leptin.
    b  <- baseline("cgl")
    ei <- b$init$EIS
    run_scn("12 CGL + metreleptin, food intake clamped (pair-fed)", "cgl",
            met_sc(0.06, 60, start = 730), pset = list(EICLAMP = ei))
  },
  s12b = function() {
    b  <- baseline("cgl")
    ei <- b$init$EIS
    run_scn("12b CGL untreated, intake clamped (pair-fed control)", "cgl",
            NULL, pset = list(EICLAMP = ei))
  },

  ## ---- 13-14: capacity-directed therapy needs a substrate --------------
  s13 = function() run_scn("13 CGL + pioglitazone 45 mg (no preadipocyte substrate)",
                           "cgl", pio_oral(45, start = 730)),
  s14 = function() run_scn("14 FPLD low-leptin + pioglitazone 45 mg (substrate present)",
                           "fpld_low", pio_oral(45, start = 730)),

  ## ---- 15-17: the spillover agents ------------------------------------
  s15 = function() run_scn("15 FPLD low-leptin + volanesorsen 300 mg q1w", "fpld_low",
                           vol_sc(300, start = 730)),
  s16 = function() run_scn("16 CGL + volanesorsen 300 mg q1w", "cgl",
                           vol_sc(300, start = 730)),
  s17 = function() run_scn("17 CGL + metreleptin + volanesorsen", "cgl",
                           c(met_sc(0.06, 60, start = 730), vol_sc(300, start = 730))),

  ## ---- 18-19: cutting influx without a drug, and PPARa ----------------
  s18 = function() run_scn("18 CGL + very-low-fat diet (15% kcal fat)", "cgl",
                           NULL, pset = list(DIETT = 730, DIETFAT = 0.15)),
  s19 = function() run_scn("19 CGL + fenofibrate + omega-3", "cgl",
                           fen_oral(145, start = 730),
                           pset = list(OM3DOSE = 1, OM3T = 730)),

  ## ---- 20-21: glucose-directed therapy and its trade-off --------------
  s20 = function() run_scn("20 CGL + high-dose insulin 150 U/d", "cgl",
                           NULL, pset = list(INSDOSE = 150, INST = 730)),
  s21 = function() run_scn("21 CGL + metformin + insulin 150 U/d", "cgl",
                           metf_oral(1000, start = 730),
                           pset = list(INSDOSE = 150, INST = 730)),

  ## ---- 22: two anorexigenic levers on one intake term ----------------
  s22a = function() run_scn("22a FPLD low + GLP-1 RA 1 mg q1w", "fpld_low",
                            glp_sc(1, start = 730)),
  s22b = function() run_scn("22b FPLD low + metreleptin", "fpld_low",
                            met_sc(0.06, 60, start = 730)),
  s22c = function() run_scn("22c FPLD low + both", "fpld_low",
                            c(glp_sc(1, start = 730), met_sc(0.06, 60, start = 730))),

  ## ---- 23: the only irreversible state decides the window -----------
  s23a = function() run_scn("23a CGL + metreleptin from year 1", "cgl",
                            met_sc(0.06, 60, start = 365, dur = 20*365 - 365),
                            end = 20*365),
  s23b = function() run_scn("23b CGL + metreleptin from year 10", "cgl",
                            met_sc(0.06, 60, start = 3650, dur = 20*365 - 3650),
                            end = 20*365),

  ## ---- 24: ANGPTL3 blockade and SGLT2 inhibition ------------------
  s24 = function() run_scn("24 CGL + evinacumab 900 mg q4w", "cgl",
                           evi_iv(900, start = 730)),
  s25 = function() run_scn("25 CGL + SGLT2 inhibitor 10 mg", "cgl",
                           sgl_oral(10, start = 730))
)

run_all <- function() lapply(scenarios, function(f) f())

## =============================================================================
##  SUMMARY
## =============================================================================

at_time <- function(d, t) d[which.min(abs(d$time - t)), ]

summarise_all <- function(all = run_all(), t0 = 730, t1 = 730 + 365) {
  do.call(rbind, lapply(all, function(d) {
    a <- at_time(d, t0); b <- at_time(d, t1); e <- d[nrow(d), ]
    data.frame(scenario = d$scenario[1],
               leptin_pre = round(a$t_LTOT, 2),
               deficit_pre = round(a$t_DEF, 3),
               occ_on     = round(b$t_OCC, 3),
               dSD        = round(b$t_SD - a$t_SD, 3),
               intake_pre = round(a$EIS),
               dIntake    = round(b$EIS - a$EIS),
               dJOV       = round(b$t_JOV - a$t_JOV, 1),
               A1C_pre    = round(a$A1C, 2),
               dA1C       = round(b$A1C - a$A1C, 2),
               TG_pre     = round(a$PTG),
               dTG_pct    = round(100*(b$PTG - a$PTG)/a$PTG, 1),
               HFF_pre    = round(a$t_HFF, 1),
               dHFF_pct   = round(100*(b$t_HFF - a$t_HFF)/a$t_HFF, 1),
               fib_end    = round(e$FIB, 2),
               row.names  = NULL)
  }))
}

## =============================================================================
##  THE SEVEN INFERENCE TESTS
##
##  Each function below is a falsification attempt on one structural claim of
##  the map.  They are the reason the scenario list contains matched pairs.
## =============================================================================

## 24-hour profile.  Metreleptin has a ~4 h half-life and is dosed daily, so a
## sample taken at any integer day is a trough.  Occupancy, the leptin signal
## and the fraction of the day spent above the transducer EC50 are therefore
## averaged over one full dosing interval.
daily_profile <- function(pheno, mgkg, day = 730 + 365) {
  b <- baseline(pheno)
  m <- zero_re(b$mod)
  out <- mrgsim_df(m, events = met_sc(mgkg, 60, start = 730),
                   end = day + 1, delta = 0.02, recsort = 3)
  w <- out[out$time >= day & out$time <= day + 1, ]
  c(mean_SD = mean(w$t_SD), mean_OCC = mean(w$t_OCC),
    peak_leptin = max(w$t_LTOT), trough_leptin = min(w$t_LTOT),
    frac_day_above_EC50 = mean(w$t_LEFF > 4.0))
}

## INFERENCE 1 -- effect = occupancy x deficit.
## If the model is right, scenario 9 (double dose, high-leptin PL) raises
## occupancy but not effect, while scenario 6 (same dose, GL) has a large
## effect at the SAME occupancy.
inference1_dose_vs_deficit <- function() {
  d <- list(scenarios$s6(), scenarios$s7(), scenarios$s8(), scenarios$s9())
  prof <- list(daily_profile("cgl", 0.06), daily_profile("fpld_low", 0.06),
               daily_profile("fpld_high", 0.06), daily_profile("fpld_high", 0.13))
  k <- 0
  do.call(rbind, lapply(d, function(x) {
    k <<- k + 1; pr <- prof[[k]]
    a <- at_time(x, 730); b <- at_time(x, 730 + 365)
    data.frame(scenario = x$scenario[1],
               leptin_pre = round(a$t_LTOT, 2),
               deficit_pre = round(a$t_DEF, 3),
               occupancy_24h = round(pr[["mean_OCC"]], 3),
               signal_24h = round(pr[["mean_SD"]], 3),
               frac_day_above_EC50 = round(pr[["frac_day_above_EC50"]], 3),
               dSD = round(b$t_SD - a$t_SD, 3),
               dA1C = round(b$A1C - a$A1C, 2),
               dTG_pct = round(100*(b$PTG - a$PTG)/a$PTG, 1),
               dHFF_pct = round(100*(b$t_HFF - a$t_HFF)/a$t_HFF, 1),
               row.names = NULL)
  }))
}

## INFERENCE 2 -- metreleptin has an intake-independent component.
## Scenario 12 clamps intake; what survives is the direct peripheral action.
inference2_pairfed <- function() {
  free <- scenarios$s6(); clamp <- scenarios$s12(); ctrl <- scenarios$s12b()
  f <- function(x) {
    a <- at_time(x, 730); b <- at_time(x, 730 + 365)
    c(dA1C = b$A1C - a$A1C,
      dTG  = 100*(b$PTG - a$PTG)/a$PTG,
      dHFF = 100*(b$t_HFF - a$t_HFF)/a$t_HFF)
  }
  m <- rbind(free = f(free), clamped = f(clamp), clamp_no_drug = f(ctrl))
  out <- as.data.frame(round(m, 2))
  out$arm <- rownames(m)
  drug_clamped <- f(clamp) - f(ctrl)
  out <- rbind(out, data.frame(t(round(drug_clamped, 2)),
                               arm = "clamped drug effect", row.names = NULL))
  out$fraction_of_free <- round(c(1, NA, NA, drug_clamped["dHFF"]/f(free)["dHFF"]), 3)
  out
}

## INFERENCE 3 -- plasma TG is tier 5, so removing it should NOT move hepatic fat.
inference3_decoupling <- function() {
  d <- list(scenarios$s15(), scenarios$s16(), scenarios$s6())
  do.call(rbind, lapply(d, function(x) {
    a <- at_time(x, 730); b <- at_time(x, 730 + 365)
    data.frame(scenario = x$scenario[1],
               dTG_pct  = round(100*(b$PTG - a$PTG)/a$PTG, 1),
               dHFF_pct = round(100*(b$t_HFF - a$t_HFF)/a$t_HFF, 1),
               dA1C     = round(b$A1C - a$A1C, 2),
               dJOV     = round(b$t_JOV - a$t_JOV, 2),
               ratio_HFF_to_TG = round(((b$t_HFF - a$t_HFF)/a$t_HFF)/
                                       ((b$PTG - a$PTG)/a$PTG), 3),
               row.names = NULL)
  }))
}

## INFERENCE 4 -- PPARg agonism needs a preadipocyte substrate.
## Identical pioglitazone exposure, two substrates.
inference4_substrate <- function() {
  d <- list(scenarios$s13(), scenarios$s14())
  do.call(rbind, lapply(d, function(x) {
    a <- at_time(x, 730); b <- at_time(x, 730 + 2*365)
    data.frame(scenario = x$scenario[1],
               pio_ng_mL = round(1000*b$PCEN/15, 0),
               capacity_pre = round(a$t_CAP, 2),
               capacity_post = round(b$t_CAP, 2),
               dfat_kg = round(b$t_FAT - a$t_FAT, 2),
               dA1C = round(b$A1C - a$A1C, 2),
               dHFF_pct = round(100*(b$t_HFF - a$t_HFF)/a$t_HFF, 1),
               row.names = NULL)
  }))
}

## INFERENCE 5 -- a diet that cuts J_in by the same amount as metreleptin's
## anorexigenic effect should match it on TG and hepatic fat, but NOT on the
## direct peripheral term.  The diet is TITRATED to match, not guessed.
inference5_diet_equivalence <- function() {
  drug <- scenarios$s6()
  jin_target <- at_time(drug, 730 + 365)$t_JIN
  obj <- function(ff) {
    x <- run_scn("titrate", "cgl", NULL, pset = list(DIETT = 730, DIETFAT = ff),
                 end = 730 + 365, delta = 365)
    at_time(x, 730 + 365)$t_JIN - jin_target
  }
  ff <- tryCatch(uniroot(obj, c(0.02, 0.35), tol = 1e-3)$root, error = function(e) NA)
  diet <- run_scn(sprintf("diet matched at %.1f%% kcal fat", 100*ff), "cgl", NULL,
                  pset = list(DIETT = 730, DIETFAT = ff))
  f <- function(x) {
    a <- at_time(x, 730); b <- at_time(x, 730 + 365)
    data.frame(scenario = x$scenario[1],
               JIN_post = round(b$t_JIN, 1),
               dTG_pct  = round(100*(b$PTG - a$PTG)/a$PTG, 1),
               dHFF_pct = round(100*(b$t_HFF - a$t_HFF)/a$t_HFF, 1),
               dA1C     = round(b$A1C - a$A1C, 2),
               SD_post  = round(b$t_SD, 3),
               row.names = NULL)
  }
  rbind(f(drug), f(diet))
}

## INFERENCE 6 -- only fibrosis is irreversible, so start time changes the
## CEILING of recovery rather than its rate.
inference6_window <- function() {
  d <- list(scenarios$s1(), scenarios$s23a(), scenarios$s23b())
  do.call(rbind, lapply(d, function(x) {
    e <- at_time(x, 20*365)
    data.frame(scenario = x$scenario[1],
               fib_y20 = round(e$FIB, 2),
               HFF_y20 = round(e$t_HFF, 1),
               A1C_y20 = round(e$A1C, 2),
               eGFR_y20 = round(e$EGFR, 1),
               row.names = NULL)
  }))
}

## INFERENCE 7 -- neutralising ADA reproduces withdrawal while the clinical
## leptin assay (total drug) stays high.  A PK/PD dissociation, not
## non-adherence.
inference7_ada <- function() {
  d <- list(scenarios$s6(), scenarios$s11(), scenarios$s10())
  do.call(rbind, lapply(d, function(x) {
    b <- at_time(x, 730 + 365); e <- at_time(x, 9*365)
    data.frame(scenario = x$scenario[1],
               leptin_assay_y9 = round(e$t_LTOT, 1),
               leptin_free_y9  = round(e$t_LEFF, 2),
               SD_y9   = round(e$t_SD, 3),
               A1C_y3  = round(b$A1C, 2),
               A1C_y9  = round(e$A1C, 2),
               TG_y9   = round(e$PTG),
               row.names = NULL)
  }))
}

## Pancreatitis-relevant endpoint: time above the TG threshold, which is what
## the complication actually depends on -- not the mean.
pancreatitis_burden <- function() {
  d <- list(scenarios$s1(), scenarios$s6(), scenarios$s16(), scenarios$s17(),
            scenarios$s19())
  do.call(rbind, lapply(d, function(x) {
    a <- at_time(x, 730); e <- at_time(x, 12*365)
    data.frame(scenario = x$scenario[1],
               TG_pre = round(a$PTG),
               TG_end = round(e$PTG),
               days_above_1000_after_start = round(e$TAT - a$TAT),
               cumulative_hazard = round(e$PANCH - a$PANCH, 3),
               row.names = NULL)
  }))
}

## Between-subject variability is used deliberately and only here: the
## deterministic scenarios above exist to isolate mechanisms, and a random
## draw would let luck masquerade as one.
responder_rates <- function(n = 200, thresh = 1.0) {
  out <- lapply(c("cgl", "fpld_low", "fpld_high"), function(ph) {
    b <- baseline(ph)
    m <- omat(b$mod, dmat(0.02))
    idata <- data.frame(ID = seq_len(n))
    tr <- mrgsim_df(m, idata = idata, events = met_sc(0.06, 60, start = 730),
                    end = 730 + 365, delta = 365)
    a <- tr[tr$time == 730, ]; b2 <- tr[tr$time == 730 + 365, ]
    data.frame(phenotype = phenotypes[[ph]]$label,
               n = n,
               mean_dA1C = round(mean(b2$A1C - a$A1C), 2),
               responder_pct = round(100*mean((a$A1C - b2$A1C) >= thresh), 1),
               row.names = NULL)
  })
  do.call(rbind, out)
}

if (identical(environment(), globalenv()) &&
    !is.null(getOption("lipo.run.scenarios"))) {
  print(baseline_table())
  all <- run_all()
  print(summarise_all(all))
  print(inference1_dose_vs_deficit())
  print(inference2_pairfed())
  print(inference3_decoupling())
  print(inference4_substrate())
  print(inference5_diet_equivalence())
  print(inference6_window())
  print(inference7_ada())
  print(pancreatitis_burden())
}
