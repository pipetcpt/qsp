## =====================================================================
##  Cancer Anorexia-Cachexia Syndrome (CACS) — mrgsolve QSP model
##  cacs_mrgsolve_model.R
## ---------------------------------------------------------------------
##  THE PREMISE
##
##  Cachexia is not starvation. In starvation the body defends itself:
##  REE falls with intake, ketogenesis spares protein, and refeeding
##  restores what was lost. In cancer cachexia all three compensations
##  break, and it is the BREAKAGE — not the calorie deficit — that makes
##  the syndrome refractory to food.
##
##  The model is built as TWO ARMS converging on muscle MASS, plus a
##  THIRD AXIS carrying muscle QUALITY that neither nutrition nor most
##  drugs ever touch:
##
##    ARM A  intake      GDF-15 -> GFRAL (area postrema) -> CALCR/PBN
##                       aversive drive -> melanocortin tone -> INTAKE.
##                       The leptin-fall rescue signal is present but is
##                       BLOCKED by central IL-1b/PGE2, which is why
##                       cachectic patients do not get hungry the way
##                       starving people do.
##    ARM B  catabolism  IL-6/LIF -> STAT3 and activin A -> SMAD2/3
##                       -> FoxO -> MAFbx/MuRF1 -> UPS + autophagy, with
##                       IGF-1/Akt/mTORC1 suppressed (ANABOLIC RESISTANCE)
##                       and PTHrP/IL-6 browning WAT so REE RISES while
##                       intake falls.
##    AXIS C quality     inflammation -> ROS -> myosin nitration, RyR1
##                       leak, PGC-1a collapse -> FORCE PER KG falls; and
##                       mass added FASTER than mitochondria are built is
##                       DILUTED, low-quality mass.
##
##      GRIP = MASS x QUALITY
##
##  That one line is why this model exists. It is the mechanistic
##  explanation for the most reproducible result in cachexia drug
##  development:
##
##      anamorelin (ROMANA 1/2)  LBM endpoint MET, handgrip MISSED
##      enobosarm  (POWER 1/2)   LBM endpoint MET, stair-climb MISSED
##
##  Both add mass quickly through anabolic signalling while the
##  inflammatory drive — and therefore quality — is untouched. The
##  dilution term (FDIL) reproduces this. FDIL is the single parameter
##  that carries the dissociation and it is called out as such in the
##  sensitivity diagnostic; everything else is anchored independently.
##
##  Two further structural features:
##
##  1. THE DISEASE IS GENERATED, NOT ASSUMED. Set TVOL0 = 0 and the
##     patient is weight-stable forever. Every state is initialised at the
##     analytic fixed point of its own ODE at zero tumour burden, and the
##     baseline diagnostic verifies ~0% drift over 365 days. Anorexia,
##     hypermetabolism, wasting and refractoriness all EMERGE from tumour
##     burden plus one susceptibility parameter (CXSENS).
##
##  2. THE POINT OF NO RETURN. Sustained STAT3 depletes the Pax7+
##     satellite-cell / myonuclear pool (SATC) with a recovery rate ~25x
##     slower than its loss rate. Once depleted, removing the tumour no
##     longer restores mass. That hysteresis IS Fearon's refractory stage,
##     and it is why in this model timing beats potency.
##
## ---------------------------------------------------------------------
##  STRUCTURE — 78 ODE compartments
##
##    Tumour and mediators (10)  TUMOR, IL6, TNFA, GDF15, ACTA, MSTN,
##                               LIF, CRP, ALB, LIVM
##    Arm A — central (7)        BSSIG, NAUS, AGRP, POMC, ANOR, LEPT, GHRL
##    Intake and energy (3)      INTK, REE, PACT
##    Arm B — muscle (9)         STAT3, SMAD, NFKB, FOXO, ATRO, UPS,
##                               AUTOP, ROS, AAPL
##    Anabolic (8)               AKT, MTOR, ARES, GH, IGF1, TEST, CORT, INSR
##    Axis C — quality (4)       PGC1, MITO, MYOQ, SATC
##    Body composition (11)      MUSC (what DXA sees), MUSCC (what pulls),
##                               LNW (GH-driven lean water), FATM, OLBM,
##                               EDEM, BWL1-3, UCP1, ATGL
##    Function and outcome (2)   ECOG, CHZ
##    Drug PK (24)               anamorelin 2, megestrol 2, dexamethasone 2,
##                               olanzapine 2, enobosarm 2, espindolol 2,
##                               EPA 2, celecoxib 2, ponsegromab 3 + complex,
##                               tocilizumab 2, bimagrumab 2
##
##  307 annotated parameters. Baseline drift with TVOL0 = 0 is 0.000000%
##  over 365 days on every reported state (diagnostic 1 below).
##
## ---------------------------------------------------------------------
##  CALIBRATION ANCHORS (full citations in cacs_references.md)
##
##   Baseline    70 kg male, 1.73 m, advanced NSCLC, pre-morbid
##               weight-stable: 25 kg skeletal muscle, 17.5 kg fat,
##               REE 1550, TEE 2250 kcal/d, grip 40 kg, CRP 2 mg/L,
##               albumin 42 g/L, GDF-15 0.5 ng/mL, IGF-1 150 ng/mL,
##               FAACT A/CS 44, muscle FSR 1.55 %/d.
##   Natural hx  Fearon 2011 Lancet Oncol staging; Martin 2015 JCO
##               %WL x BMI survival grid; untreated advanced NSCLC loses
##               roughly 5% of body weight per 3 months, accelerating.
##   ROMANA 1/2  anamorelin 100 mg qd x 12 wk: LBM +0.99 / +0.65 kg vs
##               -0.47 / -0.98 kg placebo (difference ~ +1.1 to +1.6 kg);
##               body weight ~ +2.2 kg vs placebo; HANDGRIP NO DIFFERENCE
##               (Temel 2016 Lancet Oncol).
##   Ponsegromab Groarke 2024 NEJM phase 2, 12 wk: 400 mg q4w +5.6% body
##               weight over placebo; appetite, cachexia symptoms and
##               actigraphy-measured physical activity all improved; free
##               GDF-15 driven below the limit of quantification.
##   POWER 1/2   enobosarm 3 mg: LBM co-primary met, stair-climb power
##               co-primary missed (Dobs 2013 Lancet Oncol; Crawford 2016).
##   ACT-ONE     espindolol 10 mg bid x 16 wk: weight +0.9 kg vs -1.2 kg
##               placebo AND handgrip improved (Stewart Coats 2016 JCSM).
##   Megestrol   weight gain that is fat and fluid, no lean-mass or QoL
##               benefit, oedema/VTE and adrenal suppression
##               (Ruiz-Garcia 2013 Cochrane).
##   Olanzapine  2.5 mg qd with chemotherapy: >5% weight gain in ~60% vs
##               ~9% on placebo, appetite improved (Sandhya 2023 JCO).
##   Steroids    appetite benefit for ~2-4 weeks then tachyphylaxis, with
##               proximal myopathy emerging thereafter (Yennurajalingam
##               2013 JCO; Paulsen 2014 JCO).
##   Tocilizumab CRP normalises within 1-2 weeks and weight stabilises,
##               while PLASMA IL-6 RISES (receptor blockade, not ligand
##               removal) — reproduced explicitly below.
##
## ---------------------------------------------------------------------
##  Requires: mrgsolve (>= 1.0). Run:  Rscript cacs_mrgsolve_model.R
## =====================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
})

code <- '
$PARAM @annotated
// ---------------- tumour ------------------------------------------------
TVOL0   :  25   : Tumour burden at t=0 (g; 0 = healthy control)
TMAX    : 1200  : Gompertz carrying capacity (g)
KGROW   : 0.006 : Gompertz growth rate constant (1/d)
FCACHEX :  1.0  : Histology cachexigenicity (pancreas 1.6, GI 1.3, NSCLC 1.0, breast 0.5)
KAAFEED : 0.0   : Amino-acid feed-forward on tumour growth (1/(mM*d)); >0 closes the vicious cycle

// ---------------- host susceptibility -----------------------------------
CXSENS  :  1.0  : Host cachexia sensitivity multiplier (mediator production)
SEXM    :  1.0  : 1 = male, 0 = female (body composition and grip scaling)

// ---------------- circulating mediators ---------------------------------
IL6B    :  2.0  : Baseline plasma IL-6 (pg/mL)
KOUTIL6 :  24   : IL-6 elimination rate (1/d; t1/2 ~ 40 min)
FRMED6  : 0.60  : Receptor-mediated share of IL-6 clearance (blocked by tocilizumab)
SIL6    : 1.10  : Tumour -> IL-6 secretion scaler (pg/mL/d per g)
TNFB    :  1.6  : Baseline TNF-alpha (pg/mL)
KOUTTNF :  30   : TNF elimination rate (1/d)
STNF    : 0.20  : Tumour -> TNF secretion scaler
GDFB    :  0.5  : Baseline GDF-15 (ng/mL)
KOUTGDF : 0.55  : GDF-15 elimination rate (1/d; t1/2 ~ 30 h)
SGDF    : 0.016 : Tumour -> GDF-15 secretion scaler (ng/mL/d per g)
ACTB    :  300  : Baseline activin A (pg/mL)
KOUTACT :  2.8  : Activin A elimination rate (1/d)
SACT    :  2.6  : Tumour -> activin A secretion scaler
LIFB    :  8    : Baseline LIF (pg/mL)
KOUTLIF :  20   : LIF elimination rate (1/d)
SLIF    : 0.30  : Tumour -> LIF secretion scaler
MSTNB   :  3.0  : Baseline plasma myostatin (ng/mL)
KOUTMST :  1.4  : Myostatin elimination rate (1/d)

// ---------------- hepatic acute-phase response --------------------------
CRPB    :  2.0  : Baseline CRP (mg/L)
KOUTCRP : 1.25  : CRP elimination rate (1/d; t1/2 ~ 19 h)
EMXCRP  :  30   : Emax of IL-6-driven CRP synthesis (fold)
EC50CRP :  14   : IL-6 excess for half-maximal CRP synthesis (pg/mL)
ALBB    :  42   : Baseline albumin (g/L)
KOUTALB : 0.038 : Albumin fractional turnover (1/d; t1/2 ~ 18 d)
IMXALB  : 0.42  : Maximal fractional suppression of albumin synthesis
IC50ALB :  16   : IL-6 excess for half-maximal albumin suppression (pg/mL)
LIVMB   :  1.6  : Baseline liver mass (kg)
KLIVM   : 0.03  : Liver mass turnover (1/d)
SLIVM   : 0.22  : Maximal APR-driven hepatomegaly (fraction)
EC50APR :  30   : CRP for half-maximal hepatomegaly (mg/L)

// ---------------- ARM A: brainstem GFRAL and appetite -------------------
KBS     :  1.6  : Brainstem aversive-drive turnover (1/d)
EMXBS   :  1.0  : Emax of GDF-15 on brainstem aversive drive
EC50BS  :  5.0  : Free GDF-15 for half-maximal brainstem drive (ng/mL)
HBS     :  1.8  : Hill coefficient, GDF-15 -> brainstem drive
KNAUS   :  1.2  : Nausea turnover (1/d)
SNAUS   :  3.5  : Brainstem drive -> nausea gain (0-10 scale)
KAGRP   :  2.0  : AgRP/NPY tone turnover (1/d)
KPOMC   :  2.0  : POMC tone turnover (1/d)
SLEPAG  : 1.00  : Leptin-fall rescue drive on AgRP tone
SIL1AG  : 0.62  : Central IL-1b/PGE2 BLOCKADE of the AgRP rescue
EC50IL1 :  20   : IL-6 excess for half-maximal central IL-1b effect (pg/mL)
SGHSAG  : 0.55  : GHSR-1a drive on AgRP tone
SLEPPO  : 0.40  : Leptin drive on POMC tone
SIL1PO  : 0.85  : Central inflammation drive on POMC tone
KANOR   :  1.5  : Anorexia integrator turnover (1/d)
WBSANO  : 0.70  : Weight of brainstem drive in the anorexia integrator
WMCANO  : 0.30  : Weight of melanocortin imbalance in the anorexia integrator
KLEPT   :  1.0  : Leptin turnover (1/d)
LEPKG   :  1.1  : Leptin per kg fat mass (ng/mL/kg)
KGHRL   :  6.0  : Ghrelin turnover (1/d)
GHRLB   :  550  : Baseline acyl-ghrelin (pg/mL)
SGHRWL  :  1.5  : Maximal fold rise in ghrelin with weight loss
WLGHR   : 0.20  : Fractional weight loss giving the maximal ghrelin rise
SLEAP2  : 0.45  : Inflammation-driven LEAP2 antagonism of GHSR signalling

// ---------------- intake ------------------------------------------------
INTKB   : 2250  : Weight-stable energy intake (kcal/d)
KINTK   :  1.0  : Intake adaptation rate (1/d)
IMXANO  : 0.40  : Maximal fractional intake suppression by the anorexia drive
IC50ANO : 0.55  : Anorexia drive for half-maximal intake suppression
SNAUSI  : 0.025 : Additional fractional intake loss per nausea unit
IMXBOOS : 0.10  : Maximal fractional intake GAIN from orexigenic drugs
EC50BOO : 0.60  : Orexigenic drug drive for half-maximal intake gain
GIBAR   :  0.0  : Mechanical/mucosal GI barrier (0-1 fractional intake loss)
PROTB   : 1.00  : Baseline protein intake (g/kg/d)
PROTTG  : 1.00  : Prescribed protein intake (g/kg/d)
ONSKCAL :  0.0  : ONS/EN/PN prescribed energy (kcal/d)
ONSADH  : 0.65  : Adherence to prescribed ONS (fraction)

// ---------------- energy expenditure ------------------------------------
REEB    : 1550  : Baseline resting energy expenditure (kcal/d)
KREE    : 0.35  : REE adaptation rate (1/d)
SREEIL6 : 0.16  : Maximal fractional REE rise from inflammation
EC50RE6 :  14   : IL-6 excess for half-maximal REE rise (pg/mL)
SREEUCP : 0.08  : Maximal fractional REE rise from WAT browning
SCORI   : 0.18  : Maximal fractional REE rise from Cori-cycle futile flux
SREEGDF : 0.08  : Maximal fractional REE rise from GFRAL-driven sympathetic output
EC50COR :  420  : Tumour burden for half-maximal Cori flux (g)
SREEMUS : 0.55  : Elasticity of REE to lean mass
AEEB    :  500  : Nominal activity energy expenditure (kcal/d) - the value actually used is solved in $MAIN so the energy budget closes
KPACT   : 0.20  : Physical-activity adaptation rate (1/d)
DITF    : 0.089 : Diet-induced thermogenesis (fraction of intake)
KCALFAT : 9440  : Energy density of adipose tissue (kcal/kg)
KCALMUS : 1050  : Energy density of muscle wet mass (kcal/kg)
FPARTF  :  1.0  : Fraction of residual energy balance partitioned to fat

// ---------------- ARM B: muscle catabolic signalling --------------------
KSTAT3  :  4.0  : Muscle STAT3 turnover (1/d)
SST6    : 0.85  : IL-6 contribution to STAT3 activation
EC50ST6 :  11   : IL-6 for half-maximal STAT3 (pg/mL)
SSTLIF  : 0.35  : LIF contribution to STAT3
EC50STL :  26   : LIF for half-maximal STAT3 (pg/mL)
KSMAD   :  4.0  : SMAD2/3 turnover (1/d)
SSMACT  : 0.90  : Activin A contribution to SMAD
EC50SMA :  760  : Activin A for half-maximal SMAD (pg/mL)
SSMMST  : 0.30  : Myostatin contribution to SMAD
EC50SMM :  3.4  : Myostatin for half-maximal SMAD (ng/mL)
KNFKB   :  4.0  : Muscle NF-kB turnover (1/d)
SNFTNF  : 0.80  : TNF contribution to muscle NF-kB
EC50NFT :  5.5  : TNF for half-maximal NF-kB (pg/mL)
KFOXO   :  2.0  : FoxO turnover (1/d)
WFXSTAT : 0.40  : Weight of STAT3 deviation on FoxO
WFXSMAD : 0.34  : Weight of SMAD deviation on FoxO
WFXNFKB : 0.16  : Weight of NF-kB deviation on FoxO
WFXGR   : 0.42  : Weight of glucocorticoid tone (KLF15/REDD1) on FoxO
WFXAKT  : 0.75  : Weight of Akt LOSS on FoxO (nuclear exclusion relieved)
KATRO   :  1.2  : Atrogene (MAFbx/MuRF1) transcript turnover (1/d)
EMXATR  :  1.2  : Emax of FoxO on atrogene transcription (fold above baseline)
EC50ATR : 0.50  : FoxO excess for half-maximal atrogene induction
KUPS    :  1.5  : Proteasome-flux turnover (1/d)
SUPS    : 0.20  : Atrogene -> proteasome flux gain
KAUTOP  :  1.5  : Autophagy-flux turnover (1/d)
SAUTOP  : 0.20  : FoxO -> autophagy gain
KROS    :  2.0  : Muscle ROS turnover (1/d)
SROS    : 1.25  : STAT3/NF-kB -> ROS gain
KAAPL   :  12   : Plasma amino-acid pool turnover (1/d)
AAPLB   :  3.5  : Baseline plasma BCAA-equivalent pool (mM)
SAAMUS  : 0.65  : Weight of muscle proteolysis in the plasma AA pool
SAAPRO  : 0.35  : Weight of dietary protein in the plasma AA pool

// ---------------- ARM B: anabolic signalling ----------------------------
KAKT    :  4.0  : Akt turnover (1/d)
SAKTIGF : 0.85  : IGF-1 drive on Akt
EC50AKI :  105  : IGF-1 for half-maximal Akt (ng/mL)
SAKTINS : 0.25  : Insulin-sensitivity contribution to Akt
KMTOR   :  4.0  : mTORC1 turnover (1/d)
WMTAKT  : 0.55  : Weight of Akt on mTORC1
WMTLEU  : 0.45  : Weight of the leucine/amino-acid signal on mTORC1
SREDD1  : 0.45  : Glucocorticoid REDD1 brake on mTORC1
KARES   : 0.20  : Anabolic-resistance index turnover (1/d)
SARSTAT : 0.50  : STAT3 contribution to anabolic resistance
SARSMAD : 0.25  : SMAD contribution to anabolic resistance
SARINAC : 0.30  : Inactivity (disuse) contribution to anabolic resistance
ARESMAX : 0.62  : Maximal fractional waste of an anabolic stimulus
AREX    : 0.70  : Fraction of anabolic resistance reversed by resistance training
KGH     :  6.0  : GH secretion-state turnover (1/d)
GHB     :  1.0  : Baseline GH secretory tone (relative)
SGHWL   : 0.85  : Rise in GH tone driven by hepatic GH resistance
KIGF1   : 0.35  : IGF-1 turnover (1/d; ternary-complex t1/2 ~ 2 d)
IGF1B   :  150  : Baseline IGF-1 (ng/mL)
SGHRES  : 0.52  : Maximal hepatic GH resistance (SOCS3)
EC50GHR :  13   : IL-6 excess for half-maximal GH resistance (pg/mL)
SIGFPRO : 0.22  : Nutritional contribution to IGF-1
KTEST   : 0.30  : Testosterone turnover (1/d)
TESTB   :  480  : Baseline total testosterone (ng/dL, male)
STHYPO  : 0.55  : Maximal inflammation-driven hypogonadism
EC50HYP :  12   : IL-6 excess for half-maximal hypogonadism (pg/mL)
KCORT   :  6.0  : Cortisol turnover (1/d)
CORTB   :  12   : Baseline cortisol (ug/dL)
SCORTIN : 0.55  : Maximal inflammation-driven cortisol rise
EC50COT :  14   : IL-6 excess for half-maximal cortisol rise (pg/mL)
KINSR   : 0.10  : Insulin-sensitivity index turnover (1/d)
SINSR   : 0.45  : Maximal inflammation-driven insulin resistance
EC50INS :  15   : IL-6 excess for half-maximal insulin resistance (pg/mL)

// ---------------- muscle mass flux gains --------------------------------
MUSCB   :  25   : Baseline skeletal muscle mass (kg, 70 kg male)
KDMUSB  : 0.0155: Baseline muscle protein turnover rate (1/d; FSR 1.55 %/d)
GKS     : 0.16  : Gain from the integrated anabolic index to synthesis
GENO    : 0.12  : Direct androgen-receptor (enobosarm) gain on synthesis
GBIM    : 0.18  : Direct ActRIIB-blockade (bimagrumab) gain on synthesis
KSMIN   : 0.35  : Lower bound on relative synthesis
KSMAX   : 1.80  : Upper bound on relative synthesis
HYPMAX  : 0.35  : Hypertrophy ceiling as a fraction above baseline muscle mass
GHOMEO  : 0.45  : Strength of the homeostatic drive back towards the muscle set point
GNUT    : 0.08  : Gain from nutritional adequacy to synthesis (gated by anabolic resistance)
SETFLR  : 0.55  : Fraction of the muscle set point that survives total myonuclear loss
KDMIN   : 0.60  : Lower bound on relative degradation
KDMAX   : 1.80  : Upper bound on relative degradation

// ---------------- AXIS C: muscle quality --------------------------------
KPGC1   : 0.09  : PGC-1a turnover (1/d)
SPGCINF : 0.25  : Maximal inflammatory suppression of PGC-1a
EC50PGC :  10   : IL-6 excess for half-maximal PGC-1a suppression (pg/mL)
SPGCEX  : 0.35  : Maximal exercise induction of PGC-1a
EXKPGC  :  1.0  : Exercise dose giving half-maximal PGC-1a induction
KMITO   : 0.045 : Mitochondrial content turnover (1/d; ~15 d half-life)
KMYOQ   : 0.030 : Muscle-quality turnover (1/d; ~23 d half-life)
WQMITO  : 0.55  : Weight of mitochondrial content on quality
WQOXID  : 0.45  : Weight of oxidative myosin/RyR1 integrity on quality
SQROS   : 0.28  : Maximal ROS-driven loss of contractile quality
ECROS   :  1.4  : ROS excess for half-maximal contractile damage
QNEW    : 0.30  : Contractile fraction of UNLOADED pharmacological hypertrophy
KLNW    : 0.15  : Turnover of GH-driven lean body water (1/d)
SLNW    :  2.8  : Maximal GH-driven lean body water (kg) - DXA scores this as lean mass
EC50LNW :  1.0  : GH tone above baseline giving half-maximal water retention
SQEX    : 0.06  : Direct resistance-training gain in contractile quality

// ---------------- other body composition --------------------------------
FATB    : 17.5  : Baseline fat mass (kg)
OLBMB   : 26.0  : Baseline non-muscle, non-liver lean mass (kg)
HTM     : 1.73  : Height (m)
KEDEM   : 0.10  : Oedema turnover (1/d)
SEDEM   :  9.0  : Oedema (kg) at complete loss of oncotic pressure
KWLAG   : 0.0167: Transit rate of the 3-stage 6-month body-weight memory (1/d)
FATMIN  :  1.6  : Non-mobilisable (essential) fat mass (kg)
KOLBM   : 0.004 : Non-muscle lean mass turnover (1/d)
SOLBM   : 0.18  : Maximal fractional loss of non-muscle lean mass
KATGL   : 0.30  : Adipose lipolytic tone turnover (1/d)
SATGL   :  1.9  : Maximal inflammatory/ZAG induction of lipolytic tone
KLIPEX  : 0.0022: Extra lipolysis per unit of lipolytic tone above baseline (1/d)
KUCP1   : 0.045 : WAT browning turnover (1/d)
SUCP1   :  1.0  : Maximal browning signal (PTHrP + IL-6 + sympathetic)
EC50UCP :  330  : Tumour burden for half-maximal PTHrP browning (g)

// ---------------- satellite cells / point of no return ------------------
KSATLS  : 0.004 : Satellite-cell depletion rate at full STAT3 drive (1/d)
KSATRC  : 0.00016: Satellite-cell recovery rate (1/d) — 25x slower than loss
SATFLR  : 0.20  : Floor of the regenerative contribution to synthesis
SATEX   :  3.5  : Fold acceleration of satellite recovery by resistance exercise
STATSCL : 0.60  : STAT3 excess corresponding to full catabolic drive

// ---------------- function and outcome ----------------------------------
GRIPB   :  40   : Baseline handgrip strength (kg, male)
SMIREF  :  50   : Baseline L3 skeletal muscle index (cm2/m2, male; x0.80 female)
KECOG   : 0.035 : ECOG performance-status adaptation rate (1/d)
SECOG   :  4.0  : Maximal ECOG contribution from loss of muscle function
SECOGN  : 0.09  : ECOG contribution per nausea unit
SECOGT  :  1.4  : Maximal ECOG contribution from tumour burden itself
ECT50   :  400  : Tumour burden for half-maximal ECOG contribution (g)
HAZB    : 0.00075: Baseline daily hazard in advanced cancer
HZWL    : 0.055 : Log-hazard per 1% weight loss
HZMGPS  : 0.33  : Log-hazard per mGPS point
HZECOG  : 0.30  : Log-hazard per ECOG point
HZTUM   : 0.0011: Log-hazard per gram of tumour burden

// ---------------- non-pharmacologic interventions -----------------------
EXRES   :  0.0  : Resistance-exercise dose (0 = none, 1 = full 3x/week programme)
EXAER   :  0.0  : Aerobic-exercise dose (0-1)
LEUBOL  :  0.0  : Leucine/HMB-enriched protein bolus flag (0/1)
ACTEX   : 0.30  : Fractional rise in activity energy expenditure at full training

// ---------------- anticancer therapy ------------------------------------
ONCEFF  :  0.0  : Anticancer therapy efficacy (fractional tumour kill rate, 1/d)
ONCSTART:  1e6  : Day anticancer therapy starts
MYOTOX  :  0.0  : Direct chemotherapy myotoxicity (0-1)

// ---------------- anamorelin PK/PD --------------------------------------
KAANA   : 24.0  : Anamorelin absorption rate (1/d; Tmax ~1 h)
CLANA   :  770  : Anamorelin clearance (L/d)
VANA    : 7800  : Anamorelin volume of distribution (L)
FANA    : 0.35  : Anamorelin oral bioavailability (fasted)
EMXANAG : 0.90  : Anamorelin Emax on GHSR-1a signalling
EC50ANA :  6.5  : Anamorelin for half-maximal GHSR effect (ng/mL)
EMXANAH : 0.55  : Anamorelin Emax on GH secretory tone
EC50ANH :  9.0  : Anamorelin for half-maximal GH effect (ng/mL)
EMXANAP : 0.35  : Anamorelin prokinetic relief of the GI barrier

// ---------------- megestrol acetate PK/PD -------------------------------
KAMEG   :  6.0  : Megestrol absorption rate (1/d)
CLMEG   :  240  : Megestrol clearance (L/d)
VMEG    : 5400  : Megestrol volume (L)
FMEG    : 0.30  : Megestrol bioavailability
EMXMEGA : 0.55  : Megestrol Emax on AgRP/NPY orexigenic tone
EC50MEG :  180  : Megestrol for half-maximal effect (ng/mL)
EMXMEGC : 0.30  : Megestrol Emax on cytokine production
EMXMEGG : 0.22  : Megestrol Emax as an OFF-TARGET glucocorticoid-receptor agonist
EMXMEGT : 0.70  : Megestrol Emax on gonadal testosterone suppression
EMXMEGE :  2.2  : Megestrol-driven fluid retention (kg at Emax)

// ---------------- dexamethasone PK/PD -----------------------------------
KADEX   : 12.0  : Dexamethasone absorption rate (1/d)
CLDEX   :  230  : Dexamethasone clearance (L/d)
VDEX    :  60   : Dexamethasone volume (L)
FDEX    : 0.80  : Dexamethasone bioavailability
EMXDEXA : 0.55  : Dexamethasone Emax on appetite (transient)
EC50DEX :  4.0  : Dexamethasone for half-maximal effect (ng/mL)
TAUDEXA :  16   : Time constant of appetite tachyphylaxis (d)
EMXDEXG : 0.45  : Dexamethasone Emax as a glucocorticoid-receptor agonist
EMXDEXI : 0.55  : Dexamethasone Emax on cytokine production

// ---------------- olanzapine PK/PD --------------------------------------
KAOLZ   :  4.0  : Olanzapine absorption rate (1/d)
CLOLZ   :  620  : Olanzapine clearance (L/d)
VOLZ    : 16000 : Olanzapine volume (L)
FOLZ    : 0.60  : Olanzapine bioavailability
EMXOLZA : 0.75  : Olanzapine Emax on appetite (H1 / 5-HT2C)
EC50OLZ :  3.0  : Olanzapine for half-maximal effect (ng/mL)
EMXOLZN : 0.75  : Olanzapine Emax on nausea (D2 / 5-HT3 / muscarinic)
EMXOLZI : 0.10  : Olanzapine-driven insulin resistance

// ---------------- enobosarm PK/PD ---------------------------------------
KAENO   :  6.0  : Enobosarm absorption rate (1/d)
CLENO   :  100  : Enobosarm clearance (L/d)
VENO    : 2900  : Enobosarm volume (L)
FENO    : 0.60  : Enobosarm bioavailability
EMXENO  : 0.85  : Enobosarm Emax on androgen-receptor anabolic drive
EC50ENO :  35   : Enobosarm for half-maximal effect (ng/mL)

// ---------------- espindolol PK/PD --------------------------------------
KAESP   : 12.0  : Espindolol absorption rate (1/d)
CLESP   :  580  : Espindolol clearance (L/d)
VESP    : 3600  : Espindolol volume (L)
FESP    : 0.85  : Espindolol bioavailability
EMXESPR : 0.70  : Espindolol Emax on the sympathetic/browning component of REE
EC50ESP :  3.5  : Espindolol for half-maximal effect (ng/mL)
EMXESPU : 0.40  : Espindolol Emax on proteasome flux (anti-catabolic)

// ---------------- EPA and celecoxib PK/PD -------------------------------
KAEPA   :  1.0  : EPA absorption rate (1/d)
KOUTEPA : 0.055 : EPA membrane washout rate (1/d)
VEPA    :  40   : EPA distribution volume (L)
EMXEPA  : 0.32  : EPA Emax on NF-kB-driven cytokine production
EC50EPA :  55   : EPA for half-maximal effect (ug/mL)
KACEL   : 12.0  : Celecoxib absorption rate (1/d)
CLCEL   :  660  : Celecoxib clearance (L/d)
VCEL    : 11000 : Celecoxib volume (L)
FCEL    : 0.40  : Celecoxib bioavailability
EMXCEL  : 0.30  : Celecoxib Emax on PGE2-driven central anorexia and IL-6
EC50CEL :  180  : Celecoxib for half-maximal effect (ng/mL)

// ---------------- ponsegromab PK/PD -------------------------------------
KAPON   : 0.28  : Ponsegromab SC absorption rate (1/d)
FPON    : 0.65  : Ponsegromab SC bioavailability
CLPON   : 0.22  : Ponsegromab clearance (L/d)
VPON    :  3.2  : Ponsegromab central volume (L)
QPON    : 0.40  : Ponsegromab intercompartmental clearance (L/d)
V2PON   :  2.6  : Ponsegromab peripheral volume (L)
KONPON  :  190  : Ponsegromab-GDF15 association (1/(nM*d))
KOFFPON : 0.06  : Ponsegromab-GDF15 dissociation (1/d)
KDEGCPX : 0.10  : Ponsegromab-GDF15 complex elimination rate (1/d)
KCPXEQ  :  24   : Equilibration rate of the reported complex state (1/d)
MWGDF   : 25000 : GDF-15 dimer molecular weight (g/mol)
FCRNALB : 0.35  : Fractional rise in mAb clearance at severe hypoalbuminaemia

// ---------------- tocilizumab PK/PD -------------------------------------
KATCZ   : 0.25  : Tocilizumab SC absorption rate (1/d)
FTCZ    : 0.80  : Tocilizumab SC bioavailability
CLTCZ   : 0.30  : Tocilizumab linear clearance (L/d)
VTCZ    :  3.5  : Tocilizumab volume (L)
VMTCZ   :  2.6  : Tocilizumab target-mediated Vmax (mg/d)
KMTCZ   :  1.2  : Tocilizumab target-mediated Km (mg/L)
IMXTCZ  : 0.92  : Tocilizumab maximal blockade of IL-6 signalling
IC50TCZ :  0.9  : Tocilizumab for half-maximal blockade (mg/L)

// ---------------- bimagrumab PK/PD --------------------------------------
CLBIM   : 0.26  : Bimagrumab clearance (L/d)
VBIM    :  3.0  : Bimagrumab central volume (L)
QBIM    : 0.45  : Bimagrumab intercompartmental clearance (L/d)
V2BIM   :  2.8  : Bimagrumab peripheral volume (L)
IMXBIM  : 0.90  : Bimagrumab maximal ActRIIB blockade
IC50BIM :  1.8  : Bimagrumab for half-maximal blockade (mg/L)

$CMT @annotated
TUMOR  : Tumour burden (g)
IL6    : Plasma IL-6 (pg/mL)
TNFA   : Plasma TNF-alpha (pg/mL)
GDF15  : FREE plasma GDF-15 (ng/mL)
ACTA   : Plasma activin A (pg/mL)
MSTN   : Plasma myostatin (ng/mL)
LIF    : Plasma LIF (pg/mL)
CRP    : C-reactive protein (mg/L)
ALB    : Serum albumin (g/L)
LIVM   : Liver mass (kg)
BSSIG  : Brainstem aversive drive (0-1)
NAUS   : Nausea score (0-10)
AGRP   : AgRP/NPY orexigenic tone (relative)
POMC   : POMC anorexigenic tone (relative)
ANOR   : Integrated anorexia drive
LEPT   : Plasma leptin (ng/mL)
GHRL   : Plasma acyl-ghrelin (pg/mL)
INTK   : Energy intake (kcal/d)
REE    : Resting energy expenditure (kcal/d)
PACT   : Physical activity energy expenditure (kcal/d)
STAT3  : Muscle STAT3 activation
SMAD   : Muscle SMAD2/3 activation
NFKB   : Muscle NF-kB activation
FOXO   : FoxO nuclear activity (1 = baseline)
ATRO   : Atrogene (MAFbx/MuRF1) transcript (1 = baseline)
UPS    : Proteasome flux (1 = baseline)
AUTOP  : Autophagy flux (1 = baseline)
ROS    : Muscle oxidative stress (1 = baseline)
AAPL   : Plasma amino-acid pool (mM)
AKT    : Akt activation
MTOR   : mTORC1 activity
ARES   : Anabolic-resistance index (0-1)
GH     : GH secretory tone (relative)
IGF1   : Plasma IGF-1 (ng/mL)
TEST   : Testosterone (ng/dL)
CORT   : Cortisol (ug/dL)
INSR   : Insulin sensitivity (1 = normal)
PGC1   : PGC-1alpha (1 = baseline)
MITO   : Mitochondrial content (1 = baseline)
MYOQ   : Muscle quality, force per kg (1 = baseline)
SATC   : Satellite-cell / myonuclear pool (1 = intact)
MUSC   : Skeletal muscle mass, DXA/CT-measurable (kg)
MUSCC  : CONTRACTILE muscle mass (kg) - the part that actually generates force
LNW    : GH-driven lean body water (kg) - counted as lean mass by DXA
FATM   : Fat mass (kg)
OLBM   : Non-muscle, non-liver lean mass (kg)
EDEM   : Oedema / excess extracellular fluid (kg)
BWL1   : Body-weight memory transit 1 (kg)
BWL2   : Body-weight memory transit 2 (kg)
BWL3   : Body-weight memory transit 3 = weight ~6 months ago (kg)
UCP1   : WAT browning / UCP1 (relative)
ATGL   : Adipose lipolytic tone (1 = baseline)
ECOG   : ECOG performance status (0-4)
CHZ    : Cumulative mortality hazard
ANAD   : Anamorelin gut depot (mg)
ANAC   : Anamorelin central (mg)
MEGD   : Megestrol gut depot (mg)
MEGC   : Megestrol central (mg)
DEXD   : Dexamethasone gut depot (mg)
DEXC   : Dexamethasone central (mg)
OLZD   : Olanzapine gut depot (mg)
OLZC   : Olanzapine central (mg)
ENOD   : Enobosarm gut depot (mg)
ENOC   : Enobosarm central (mg)
ESPD   : Espindolol gut depot (mg)
ESPC   : Espindolol central (mg)
EPAD   : EPA gut depot (mg)
EPAC   : EPA membrane pool (mg)
CELD   : Celecoxib gut depot (mg)
CELC   : Celecoxib central (mg)
POND   : Ponsegromab SC depot (nmol)
PONC   : Ponsegromab central (nmol)
PONP   : Ponsegromab peripheral (nmol)
CPX    : Ponsegromab-GDF15 complex (nmol)
TCZD   : Tocilizumab SC depot (mg)
TCZC   : Tocilizumab central (mg)
BIMC   : Bimagrumab central (mg)
BIMP   : Bimagrumab peripheral (mg)

$GLOBAL
#define CP_ANA   (ANAC/VANA*1e6)        // ng/mL
#define CP_MEG   (MEGC/VMEG*1e6)        // ng/mL
#define CP_DEX   (DEXC/VDEX*1e6)        // ng/mL
#define CP_OLZ   (OLZC/VOLZ*1e6)        // ng/mL
#define CP_ENO   (ENOC/VENO*1e6)        // ng/mL
#define CP_ESP   (ESPC/VESP*1e6)        // ng/mL
#define CP_EPA   (EPAC/VEPA)            // ug/mL
#define CP_CEL   (CELC/VCEL*1e6)        // ng/mL
#define CP_TCZ   (TCZC/VTCZ)            // mg/L
#define CP_BIM   (BIMC/VBIM)            // mg/L
#define CP_PON   (PONC/VPON)            // nM
#define POS(x)   ((x) > 0.0 ? (x) : 0.0)

double emaxf(double c, double emax, double ec50) {
  double cc = c > 0.0 ? c : 0.0;
  return emax * cc / (ec50 + cc);
}
double emaxh(double c, double emax, double ec50, double h) {
  double cc = c > 0.0 ? c : 0.0;
  double num = pow(cc, h);
  return emax * num / (pow(ec50, h) + num);
}
double clampd(double x, double lo, double hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

$MAIN
// =======================================================================
//  Initial conditions.
//  Every state is set to the ANALYTIC fixed point of its own ODE at
//  TUMOR = 0, so a TVOL0 = 0 run has no numerical drift. This is what
//  makes "the disease is generated, not assumed" a checkable claim
//  rather than a slogan.
// =======================================================================
TUMOR_0 = TVOL0;
IL6_0   = IL6B;
TNFA_0  = TNFB;
GDF15_0 = GDFB;
ACTA_0  = ACTB;
MSTN_0  = MSTNB;
LIF_0   = LIFB;
CRP_0   = CRPB;
ALB_0   = ALBB;
LIVM_0  = LIVMB;

BSSIG_0 = emaxh(GDFB, EMXBS, EC50BS, HBS);
NAUS_0  = SNAUS * BSSIG_0;
AGRP_0  = 1.0;
POMC_0  = 1.0;
ANOR_0  = WBSANO * BSSIG_0;          // melanocortin imbalance is exactly 0 at baseline
LEPT_0  = LEPKG * FATB;
GHRL_0  = GHRLB;

// The baseline anorexia drive is small but NOT zero (GDF-15 is 0.5 ng/mL
// even in health), so INTK_0 is the fixed point of the intake ODE, and the
// baseline activity expenditure is then SOLVED so that TEE == INTK exactly.
// Without this the healthy control silently loses ~1.6% of its fat per year
// and every treatment effect is measured against a moving baseline.
double SUP0 = IMXANO * ANOR_0 / (IC50ANO + fabs(ANOR_0)) + SNAUSI * NAUS_0;
INTK_0  = INTKB * (1.0 - SUP0);
REE_0   = REEB;
double INTKREF = INTK_0;
double PROTREF = PROTB * (INTK_0 / INTKB);
PACT_0  = INTK_0 * (1.0 - DITF) - REEB;   // solved: TEE == INTK exactly at t = 0
double PACTREF = PACT_0;

STAT3_0 = emaxf(IL6B, SST6, EC50ST6) + emaxf(LIFB, SSTLIF, EC50STL);
SMAD_0  = emaxf(ACTB, SSMACT, EC50SMA) + emaxf(MSTNB, SSMMST, EC50SMM);
NFKB_0  = emaxf(TNFB, SNFTNF, EC50NFT);
AKT_0   = emaxf(IGF1B, SAKTIGF, EC50AKI) + SAKTINS;
MTOR_0  = WMTAKT * AKT_0 + WMTLEU * 1.0;
FOXO_0  = 1.0;
ATRO_0  = 1.0;
UPS_0   = 1.0;
AUTOP_0 = 1.0;
ROS_0   = 1.0;
AAPL_0  = AAPLB;
ARES_0  = 0.0;
GH_0    = GHB;
IGF1_0  = IGF1B;
TEST_0  = TESTB * (SEXM > 0.5 ? 1.0 : 0.06);
CORT_0  = CORTB;
INSR_0  = 1.0;

PGC1_0  = 1.0;
MITO_0  = 1.0;
MYOQ_0  = 1.0;
SATC_0  = 1.0;

MUSC_0  = MUSCB * (SEXM > 0.5 ? 1.0 : 0.72);
MUSCC_0 = MUSC_0;
LNW_0   = 0.0;
FATM_0  = FATB  * (SEXM > 0.5 ? 1.0 : 1.55);
OLBM_0  = OLBMB * (SEXM > 0.5 ? 1.0 : 0.86);
EDEM_0  = 0.0;

double BW0 = MUSC_0 + FATM_0 + OLBM_0 + LIVMB;
BWL1_0  = BW0;
BWL2_0  = BW0;
BWL3_0  = BW0;

UCP1_0  = 0.0;
ATGL_0  = 1.0;
ECOG_0  = SECOGN * NAUS_0 + SECOGT * TVOL0 / (ECT50 + TVOL0);
double ECOGREF = ECOG_0;
CHZ_0   = 0.0;

// Reference (baseline) values used as the zero point for every deviation
double MUSCREF = MUSC_0;
double OLBMREF = OLBM_0;
double GRIPREF = GRIPB * (SEXM > 0.5 ? 1.0 : 0.60);
double BWREF   = BW0;
double AKTREF  = AKT_0;
double MTORREF = MTOR_0;
double STATREF = STAT3_0;
double SMADREF = SMAD_0;
double NFKBREF = NFKB_0;
double LEPREF  = LEPT_0;
double BSSREF  = BSSIG_0;
double TESTREF = TEST_0;

$ODE
// =======================================================================
//  0. Derived quantities
// =======================================================================
double BW    = MUSC + FATM + OLBM + LIVM + EDEM + LNW;
double PCTWL = 100.0 * (BWL3 - BW) / (BWL3 > 1.0 ? BWL3 : 1.0);
double BMI   = BW / (HTM * HTM);
double MGPS  = (CRP > 10.0 ? 1.0 : 0.0) + ((CRP > 10.0 && ALB < 35.0) ? 1.0 : 0.0);
double ONCON = (SOLVERTIME >= ONCSTART) ? 1.0 : 0.0;

// -------- drug effects --------------------------------------------------
double E_ANA_G = emaxf(CP_ANA, EMXANAG, EC50ANA);
double E_ANA_H = emaxf(CP_ANA, EMXANAH, EC50ANH);
double E_ANA_P = emaxf(CP_ANA, EMXANAP, EC50ANA);
double E_MEG_A = emaxf(CP_MEG, EMXMEGA, EC50MEG);
double E_MEG_C = emaxf(CP_MEG, EMXMEGC, EC50MEG);
double E_MEG_G = emaxf(CP_MEG, EMXMEGG, EC50MEG);
double E_MEG_T = emaxf(CP_MEG, EMXMEGT, EC50MEG);
double E_MEG_E = emaxf(CP_MEG, EMXMEGE, EC50MEG);
double TACHY   = exp(-POS(SOLVERTIME - 3.0) / TAUDEXA);
double E_DEX_A = emaxf(CP_DEX, EMXDEXA, EC50DEX) * TACHY;
double E_DEX_G = emaxf(CP_DEX, EMXDEXG, EC50DEX);
double E_DEX_I = emaxf(CP_DEX, EMXDEXI, EC50DEX);
double E_OLZ_A = emaxf(CP_OLZ, EMXOLZA, EC50OLZ);
double E_OLZ_N = emaxf(CP_OLZ, EMXOLZN, EC50OLZ);
double E_OLZ_I = emaxf(CP_OLZ, EMXOLZI, EC50OLZ);
double E_ENO   = emaxf(CP_ENO, EMXENO,  EC50ENO);
double E_ESP_R = emaxf(CP_ESP, EMXESPR, EC50ESP);
double E_ESP_U = emaxf(CP_ESP, EMXESPU, EC50ESP);
double E_EPA   = emaxf(CP_EPA, EMXEPA,  EC50EPA);
double E_CEL   = emaxf(CP_CEL, EMXCEL,  EC50CEL);
double E_TCZ   = emaxf(CP_TCZ, IMXTCZ,  IC50TCZ);
double E_BIM   = emaxf(CP_BIM, IMXBIM,  IC50BIM);

// Anti-inflammatory pressure on cytokine PRODUCTION (not on the receptor)
double ANTIINF = clampd(1.0 - E_EPA - E_CEL - E_MEG_C - E_DEX_I, 0.10, 1.0);

// Total glucocorticoid-receptor tone: endogenous + megestrol + dexamethasone
double GRTONE  = (CORT / CORTB) + E_DEX_G + E_MEG_G;

// =======================================================================
//  1. Tumour — the single upstream driver
// =======================================================================
double AAFEED = KAAFEED * POS(AAPL - AAPLB);
double GOMP   = (TUMOR > 1e-6 && TMAX > TUMOR)
                ? (KGROW + AAFEED) * TUMOR * log(TMAX / TUMOR) : 0.0;
dxdt_TUMOR = GOMP - ONCON * ONCEFF * TUMOR;

double TUMD = TUMOR * FCACHEX * CXSENS;

// =======================================================================
//  2. Circulating mediators
// =======================================================================
// NOTE ON TOCILIZUMAB: it blocks the RECEPTOR, not the ligand, so plasma
// IL-6 RISES on treatment while every IL-6-dependent readout falls. Only
// IL6EFF (the signal actually delivered) is blocked. This is real, and it
// is routinely misread in the clinic as treatment failure.
// IL-6 is cleared substantially through its own receptor, so blocking IL-6R
// SLOWS IL-6 elimination: plasma IL-6 RISES several-fold on tocilizumab while
// every IL-6-dependent readout falls. A rising IL-6 on treatment is target
// engagement, not treatment failure, and the model must say so.
dxdt_IL6  = KOUTIL6 * IL6B + SIL6 * TUMD * ANTIINF
            - KOUTIL6 * (1.0 - FRMED6 * E_TCZ) * IL6;
dxdt_TNFA = KOUTTNF * TNFB + STNF * TUMD * ANTIINF - KOUTTNF * TNFA;
dxdt_LIF  = KOUTLIF * LIFB + SLIF * TUMD * ANTIINF - KOUTLIF * LIF;
dxdt_ACTA = KOUTACT * ACTB + SACT * TUMD           - KOUTACT * ACTA;

// GDF15 is carried as TOTAL circulating GDF-15 (ng/mL). Ponsegromab binds
// with sub-picomolar affinity (Kd = KOFFPON/KONPON ~ 0.3 pM), which is orders
// of magnitude faster than any other process here, so the bound fraction is
// solved at quasi-equilibrium rather than integrated — the explicit on/off
// form is numerically stiff and buys nothing.
//   1 ng/mL of a 25 kDa dimer = 1000/MWGDF nM.
double GDFTOTNM = GDF15 * 1000.0 / MWGDF;             // ng/mL -> nM
double KDPON    = KOFFPON / KONPON;                   // nM
double FREEFR   = KDPON / (KDPON + CP_PON);           // free fraction
double GDFREE   = GDF15 * FREEFR;                     // ng/mL, the ligand GFRAL sees
dxdt_GDF15 = KOUTGDF * GDFB + SGDF * TUMD
             - KOUTGDF * GDF15 * FREEFR               // free ligand cleared normally
             - KDEGCPX * GDF15 * (1.0 - FREEFR);      // bound ligand cleared with the mAb

// The human myostatin paradox: muscle is the source, so it FALLS
dxdt_MSTN = KOUTMST * MSTNB * (MUSC / MUSCREF) - KOUTMST * MSTN;

double IL6EFF = IL6 * (1.0 - E_TCZ);
double IL6EXC = POS(IL6EFF - IL6B);

// =======================================================================
//  3. Hepatic acute-phase response
// =======================================================================
dxdt_CRP  = KOUTCRP * CRPB * (1.0 + emaxf(IL6EXC, EMXCRP, EC50CRP))
            - KOUTCRP * CRP;
dxdt_ALB  = KOUTALB * ALBB * (1.0 - emaxf(IL6EXC, IMXALB, IC50ALB))
            - KOUTALB * ALB;
dxdt_LIVM = KLIVM * LIVMB * (1.0 + SLIVM * emaxf(POS(CRP - CRPB), 1.0, EC50APR))
            - KLIVM * LIVM;

// =======================================================================
//  4. ARM A(i) — GDF-15 -> GFRAL -> brainstem aversive drive
// =======================================================================
dxdt_BSSIG = KBS * emaxh(GDFREE, EMXBS, EC50BS, HBS) - KBS * BSSIG;
dxdt_NAUS  = KNAUS * SNAUS * BSSIG * (1.0 - E_OLZ_N) - KNAUS * NAUS;

// =======================================================================
//  5. ARM A(ii) — melanocortin tone, and the rescue that fails
// =======================================================================
double CENTINF = emaxf(IL6EXC, 1.0, EC50IL1) * (1.0 - E_CEL);
double LEPREL  = LEPT / (LEPREF > 1e-9 ? LEPREF : 1e-9);
double GHSRSIG = (GHRL / GHRLB) * (1.0 - SLEAP2 * CENTINF) + E_ANA_G + E_MEG_A;

// AgRP SHOULD rise as leptin falls. Central IL-1b/PGE2 blocks exactly that.
double AGRPTGT = 1.0 + SLEPAG * (1.0 - LEPREL)
                     - SIL1AG * CENTINF
                     + SGHSAG * (GHSRSIG - 1.0)
                     + E_OLZ_A + E_DEX_A;
dxdt_AGRP = KAGRP * POS(AGRPTGT) - KAGRP * AGRP;

double POMCTGT = 1.0 + SIL1PO * CENTINF - SLEPPO * (1.0 - LEPREL);
dxdt_POMC = KPOMC * POS(POMCTGT) - KPOMC * POMC;

// Signed melanocortin imbalance: negative = net orexigenic (drug-driven)
double MCIMB   = (POMC - AGRP) / (POMC + AGRP + 1e-9);
double ANORTGT = clampd(WBSANO * BSSIG + WMCANO * MCIMB, -0.30, 1.0);
dxdt_ANOR = KANOR * ANORTGT - KANOR * ANOR;

dxdt_LEPT = KLEPT * LEPKG * FATM - KLEPT * LEPT;
double WLFRAC = clampd(POS(1.0 - BW / BWREF) / WLGHR, 0.0, 1.0);
dxdt_GHRL = KGHRL * GHRLB * (1.0 + SGHRWL * WLFRAC) - KGHRL * GHRL;

// =======================================================================
//  6. Energy intake
// =======================================================================
double ANOSUP  = IMXANO * ANOR / (IC50ANO + fabs(ANOR));
double GIEFF   = clampd(GIBAR * (1.0 - E_ANA_P), 0.0, 0.9);
double SUPTOT  = clampd(ANOSUP + SNAUSI * NAUS + GIEFF, -0.20, 0.90);
double OREXDRV = E_ANA_G + E_OLZ_A + E_DEX_A + E_MEG_A;
double BOOST   = emaxf(OREXDRV, IMXBOOS, EC50BOO);
double INTKTGT = INTKB * (1.0 - SUPTOT + BOOST) + ONSKCAL * ONSADH;
dxdt_INTK = KINTK * INTKTGT - KINTK * INTK;

double PROTACT = PROTTG * (INTK / INTKB);
double LEUSIG  = pow(PROTACT / PROTREF, 0.55) * (1.0 + 0.25 * LEUBOL);

// =======================================================================
//  7. Energy expenditure — REE rises while intake falls
// =======================================================================
// GFRAL signalling raises sympathetic output and thermogenesis as well as
// suppressing intake, so blocking GDF-15 moves BOTH sides of the energy
// budget. That is why an anti-GDF-15 antibody can produce weight gain out of
// proportion to the appetite change it produces.
double REEUP  = (emaxf(IL6EXC, SREEIL6, EC50RE6)
                 + SREEUCP * UCP1
                 + SREEGDF * POS(BSSIG - BSSREF)
                 + SCORI * TUMOR / (EC50COR + TUMOR)) * (1.0 - E_ESP_R);
double LEANF  = pow(MUSC / MUSCREF, SREEMUS);
dxdt_REE = KREE * REEB * LEANF * (1.0 + REEUP) - KREE * REE;

double ACTTGT = POS(PACTREF * (1.0 - 0.24 * ECOG) / (1.0 - 0.24 * ECOGREF)
                    * (1.0 + ACTEX * (EXRES + EXAER)));
if (ACTTGT < 60.0) ACTTGT = 60.0;
dxdt_PACT = KPACT * ACTTGT - KPACT * PACT;

double TEE  = REE + PACT + DITF * INTK;
double EBAL = INTK - TEE;

// =======================================================================
//  8. ARM B(i) — muscle catabolic signalling
//     Everything is expressed as a DEVIATION from its own baseline, so
//     the catabolic drive is exactly zero in a tumour-free host.
// =======================================================================
dxdt_STAT3 = KSTAT3 * (emaxf(IL6EFF, SST6, EC50ST6)
                       + emaxf(LIF, SSTLIF, EC50STL)) - KSTAT3 * STAT3;

double ACTAEFF = ACTA * (1.0 - E_BIM);
double MSTNEFF = MSTN * (1.0 - E_BIM);
dxdt_SMAD = KSMAD * (emaxf(ACTAEFF, SSMACT, EC50SMA)
                     + emaxf(MSTNEFF, SSMMST, EC50SMM)) - KSMAD * SMAD;

dxdt_NFKB = KNFKB * (emaxf(TNFA, SNFTNF, EC50NFT) * (1.0 - E_EPA)
                     + 0.35 * MYOTOX) - KNFKB * NFKB;

// Signed: ActRIIB blockade pushes SMAD BELOW baseline and that is the
// entire point of bimagrumab, so the term must be allowed to go negative.
double CATDRV = WFXSTAT * (STAT3 - STATREF)
              + WFXSMAD * (SMAD  - SMADREF)
              + WFXNFKB * (NFKB  - NFKBREF)
              + WFXGR   * (GRTONE - 1.0)
              + WFXAKT  * (AKTREF - AKT);
dxdt_FOXO = KFOXO * clampd(1.0 + CATDRV, 0.40, 3.0) - KFOXO * FOXO;

double FDEV = FOXO - 1.0;
dxdt_ATRO  = KATRO * (1.0 + EMXATR * FDEV / (EC50ATR + fabs(FDEV))) - KATRO * ATRO;
dxdt_UPS   = KUPS * (1.0 + SUPS * (ATRO - 1.0) * (1.0 - E_ESP_U)) - KUPS * UPS;
dxdt_AUTOP = KAUTOP * (1.0 + SAUTOP * FDEV) - KAUTOP * AUTOP;
dxdt_ROS   = KROS * (1.0 + SROS * (POS(STAT3 - STATREF)
                                   + 0.6 * POS(NFKB - NFKBREF))
                     + 0.5 * MYOTOX) - KROS * ROS;

// =======================================================================
//  9. ARM B(ii) — anabolic signalling and ANABOLIC RESISTANCE
// =======================================================================
double GHRESIST = emaxf(IL6EXC, SGHRES, EC50GHR);
dxdt_GH = KGH * GHB * (1.0 + SGHWL * GHRESIST + E_ANA_H) - KGH * GH;

double NUTADQ = (INTK / INTKREF) * pow(PROTACT / PROTREF, 0.5);
dxdt_IGF1 = KIGF1 * IGF1B * (GH / GHB) * (1.0 - GHRESIST)
            * (1.0 - SIGFPRO * POS(1.0 - NUTADQ)) - KIGF1 * IGF1;

dxdt_TEST = KTEST * TESTREF * (1.0 - emaxf(IL6EXC, STHYPO, EC50HYP))
            * (1.0 - E_MEG_T) - KTEST * TEST;
dxdt_CORT = KCORT * CORTB * (1.0 + emaxf(IL6EXC, SCORTIN, EC50COT)) - KCORT * CORT;
dxdt_INSR = KINSR * clampd(1.0 - emaxf(IL6EXC, SINSR, EC50INS)
                           - E_OLZ_I - 0.15 * E_DEX_G, 0.15, 1.0) - KINSR * INSR;

dxdt_AKT  = KAKT * (emaxf(IGF1, SAKTIGF, EC50AKI) + SAKTINS * INSR) - KAKT * AKT;
dxdt_MTOR = KMTOR * POS((WMTAKT * AKT + WMTLEU * LEUSIG)
                        * (1.0 - SREDD1 * POS(GRTONE - 1.0))) - KMTOR * MTOR;

double INACT   = POS(1.0 - PACT / PACTREF);
double ARESTGT = ARESMAX * clampd(SARSTAT * POS(STAT3 - STATREF) / STATSCL
                                + SARSMAD * POS(SMAD - SMADREF) / 0.50
                                + SARINAC * INACT, 0.0, 1.0)
                 * (1.0 - AREX * EXRES);
dxdt_ARES = KARES * ARESTGT - KARES * ARES;

// =======================================================================
// 10. Satellite-cell pool — THE POINT OF NO RETURN
// =======================================================================
dxdt_SATC = KSATRC * (1.0 + SATEX * EXRES) * (1.0 - SATC)
            - KSATLS * POS(STAT3 - STATREF) / STATSCL * SATC;
double SATEFF = SATFLR + (1.0 - SATFLR) * clampd(SATC, 0.0, 1.0);   // reported only

// =======================================================================
// 11. Muscle mass — a SMALL imbalance in a LARGE flux
//     Baseline FSR is 1.55 %/d, so a 2 kg loss over 12 weeks is only a
//     ~6% imbalance between synthesis and degradation. Getting that scale
//     right is the difference between a model and a cartoon.
// =======================================================================
// SATEFF is deliberately NOT applied here: the satellite-cell pool acts
// through the SET POINT (MSET) below. Applying it in both places would
// double-count one lesion and block recovery for the wrong reason.
double XANAB = (MTOR / MTORREF) * (1.0 - ARES);
// DRUGANA is the part of the synthesis drive supplied directly by a drug
// acting on the myonucleus (androgen receptor, ActRIIB blockade) rather than
// by relief of the disease. It is the part that arrives without a matching
// mechanical stimulus, and it is the part that carries the QNEW penalty.
double DRUGANA = GENO * E_ENO + GBIM * E_BIM;

// THE MUSCLE SET POINT. Without one, a proportional turnover model has no
// memory: any mass is an equilibrium and nothing ever recovers after the
// tumour is removed. The set point is carried by the myonuclear/satellite
// pool, so as SATC is depleted the set point itself falls and cannot be
// argued back up — that is the point of no return, expressed as a number.
double MSET  = MUSCREF * (SETFLR + (1.0 - SETFLR) * clampd(SATC, 0.0, 1.0));
double HOMEO = GHOMEO * (MSET / (MUSC > 1e-6 ? MUSC : 1e-6) - 1.0);
// Substrate availability contributes to synthesis, but the contribution is
// GATED BY ANABOLIC RESISTANCE — which is the whole reason feeding a
// cachectic patient more protein does not simply rebuild the muscle.
double NUTGAIN = GNUT * (NUTADQ - 1.0) * (1.0 - ARES / ARESMAX);
double KSREL = clampd(1.0 + GKS * (XANAB - 1.0) + DRUGANA + HOMEO + NUTGAIN,
                      KSMIN, KSMAX);
double KDREL = clampd(0.72 * UPS + 0.28 * AUTOP, KDMIN, KDMAX);
// Synthesis is proportional to the muscle mass that is actually there, so
// the whole equation is a small fractional imbalance in a large flux:
//   dM/dt = M * FSR_baseline * (relative synthesis - relative degradation)
// A 2 kg loss over 12 weeks is a ~6% imbalance, not a collapse. HYPCAP is a
// hypertrophy ceiling so that anabolic agents cannot grow muscle without end.
double HYPCAP = (MUSC <= MUSCREF) ? 1.0
                : clampd(1.0 - (MUSC / MUSCREF - 1.0) / HYPMAX, 0.0, 1.0);
double KSEFF  = KSREL * HYPCAP;
double KDRATE = KDMUSB * KDREL;                // 1/d
double KSFLX  = KDMUSB * MUSC * KSEFF;         // kg/d
dxdt_MUSC = KSFLX - KDRATE * MUSC;

// MASS THAT DXA SEES vs MASS THAT PULLS.
// Hypertrophy driven pharmacologically but WITHOUT a mechanical loading
// stimulus is only QNEW (~30%) as contractile per kilogram as normal
// muscle: myofibrillar protein is added faster than myonuclei, capillaries
// and neural drive can be matched to it. Resistance training supplies the
// missing stimulus, so EXRES restores the quality of the new tissue.
// Degradation removes contractile and non-contractile mass alike.
double QNEWEFF = clampd(QNEW + (1.0 - QNEW) * EXRES, 0.0, 1.0);
dxdt_MUSCC = KDMUSB * MUSC * (KSEFF - (1.0 - QNEWEFF) * DRUGANA * HYPCAP)
             - KDRATE * MUSCC;

// GH-DRIVEN LEAN BODY WATER.
// GHSR-1a agonists raise GH, and GH retains sodium and water. DXA cannot
// tell that water from muscle and reports both as "lean body mass", which
// is exactly why a GH-axis agent can win an LBM endpoint that a hand
// dynamometer will not confirm. Oedema was a real adverse event in ROMANA.
dxdt_LNW = KLNW * (SLNW * emaxf(POS(GH / GHB - 1.0), 1.0, EC50LNW) - LNW);

// Plasma amino-acid pool: 35% dietary, 65% from proteolysis. This is the
// currency of the vicious cycle — the muscle literally feeds the tumour
// and the hepatic acute-phase response.
dxdt_AAPL = KAAPL * AAPLB * (SAAPRO * (PROTACT / PROTREF)
                             + SAAMUS * (KDRATE * MUSC) / (KDMUSB * MUSCREF))
            - KAAPL * AAPL;

// =======================================================================
// 12. AXIS C — muscle QUALITY (where mass and function part company)
// =======================================================================
double EXTOT  = EXAER + 0.6 * EXRES;
dxdt_PGC1 = KPGC1 * POS(1.0 - emaxf(IL6EXC, SPGCINF, EC50PGC)
                        + SPGCEX * EXTOT / (EXKPGC + EXTOT)) - KPGC1 * PGC1;
dxdt_MITO = KMITO * PGC1 - KMITO * MITO;

double OXDAM = emaxf(POS(ROS - 1.0), SQROS, ECROS);
double QTGT  = POS(WQMITO * MITO + WQOXID * (1.0 - OXDAM) + SQEX * EXRES);

// MYOQ is the force per kilogram of CONTRACTILE tissue: mitochondrial
// capacity and the integrity of myosin and RyR1. It is driven down by
// inflammation and up by loading, and no purely anabolic agent touches it.
dxdt_MYOQ = KMYOQ * QTGT - KMYOQ * MYOQ;

// =======================================================================
// 13. Adipose tissue and the energy partition
// =======================================================================
dxdt_ATGL = KATGL * (1.0 + SATGL * (emaxf(IL6EXC, 0.6, 12.0)
                                    + emaxf(TUMD, 0.4, 350.0))) - KATGL * ATGL;
dxdt_UCP1 = KUCP1 * SUCP1 * (TUMOR / (EC50UCP + TUMOR)) * (1.0 - E_ESP_R)
            - KUCP1 * UCP1;

// Muscle turnover has already claimed its own energy; the residual balance
// is partitioned to fat, plus a directly cytokine/ZAG-driven lipolysis.
double DEMUS   = (KSFLX - KDRATE * MUSC) * KCALMUS;
double FATFLUX = (EBAL - DEMUS) * FPARTF / KCALFAT;
double LIPEX   = KLIPEX * POS(ATGL - 1.0) * POS(FATM - FATMIN);
dxdt_FATM = FATFLUX - LIPEX;
if (FATM <= FATMIN && dxdt_FATM < 0.0) dxdt_FATM = 0.0;

dxdt_OLBM = KOLBM * (OLBMREF * (1.0 - SOLBM * clampd(FDEV, 0.0, 1.0)) - OLBM);

// Oedema from hypoalbuminaemia (and megestrol) MASKS true tissue loss
dxdt_EDEM = KEDEM * (SEDEM * POS(1.0 - ALB / ALBB) + E_MEG_E) - KEDEM * EDEM;

// Six-month body-weight memory (3-stage transit; mean 3/KWLAG = 180 d)
dxdt_BWL1 = KWLAG * (BW   - BWL1);
dxdt_BWL2 = KWLAG * (BWL1 - BWL2);
dxdt_BWL3 = KWLAG * (BWL2 - BWL3);

// =======================================================================
// 14. Function, performance status, hazard
// =======================================================================
double GRIPNOW  = GRIPREF * (MUSCC / MUSCREF) * MYOQ;
double FUNCLOSS = clampd(1.0 - GRIPNOW / GRIPREF, 0.0, 1.0);
double ECOGTGT  = clampd(SECOG * pow(FUNCLOSS, 1.25) + SECOGN * NAUS
                         + SECOGT * TUMOR / (ECT50 + TUMOR), 0.0, 4.0);
dxdt_ECOG = KECOG * ECOGTGT - KECOG * ECOG;

dxdt_CHZ = HAZB * exp(HZWL * POS(PCTWL) + HZMGPS * MGPS
                      + HZECOG * ECOG + HZTUM * TUMOR);

// =======================================================================
// 15. Pharmacokinetics
// =======================================================================
dxdt_ANAD = -KAANA * ANAD;
dxdt_ANAC =  KAANA * ANAD - (CLANA / VANA) * ANAC;
dxdt_MEGD = -KAMEG * MEGD;
dxdt_MEGC =  KAMEG * MEGD - (CLMEG / VMEG) * MEGC;
dxdt_DEXD = -KADEX * DEXD;
dxdt_DEXC =  KADEX * DEXD - (CLDEX / VDEX) * DEXC;
dxdt_OLZD = -KAOLZ * OLZD;
dxdt_OLZC =  KAOLZ * OLZD - (CLOLZ / VOLZ) * OLZC;
dxdt_ENOD = -KAENO * ENOD;
dxdt_ENOC =  KAENO * ENOD - (CLENO / VENO) * ENOC;
dxdt_ESPD = -KAESP * ESPD;
dxdt_ESPC =  KAESP * ESPD - (CLESP / VESP) * ESPC;
dxdt_EPAD = -KAEPA * EPAD;
dxdt_EPAC =  KAEPA * EPAD - KOUTEPA * EPAC;
dxdt_CELD = -KACEL * CELD;
dxdt_CELC =  KACEL * CELD - (CLCEL / VCEL) * CELC;

// FcRn recycling is albumin-dependent, so hypoalbuminaemia SPEEDS UP mAb
// clearance: the sickest patients get the least drug exposure.
double FCRNF = 1.0 + FCRNALB * clampd(POS(1.0 - ALB / ALBB) / 0.40, 0.0, 1.0);

dxdt_POND = -KAPON * POND;
dxdt_PONC =  KAPON * POND - (CLPON * FCRNF / VPON) * PONC
             - QPON * (PONC / VPON - PONP / V2PON)
             - KDEGCPX * CPX;                 // antibody lost with the complex
dxdt_PONP =  QPON * (PONC / VPON - PONP / V2PON);
// CPX simply reports the bound pool (nmol) implied by the equilibrium
dxdt_CPX  =  KCPXEQ * (GDFTOTNM * (1.0 - FREEFR) * VPON - CPX);

dxdt_TCZD = -KATCZ * TCZD;
dxdt_TCZC =  KATCZ * TCZD - CLTCZ * FCRNF * CP_TCZ
             - VMTCZ * CP_TCZ / (KMTCZ + CP_TCZ);

dxdt_BIMC = -(CLBIM / VBIM) * BIMC - QBIM * (BIMC / VBIM - BIMP / V2BIM);
dxdt_BIMP =  QBIM * (BIMC / VBIM - BIMP / V2BIM);

$TABLE
// mrgsolve hoists $MAIN/$ODE/$TABLE locals into one scope, so names already
// declared in $ODE are re-assigned here rather than re-declared.
BW             = MUSC + FATM + OLBM + LIVM + EDEM + LNW;
double LBM     = MUSC + OLBM + LIVM + LNW;   // as DXA would report it
BMI            = BW / (HTM * HTM);
PCTWL          = 100.0 * (BWL3 - BW) / (BWL3 > 1.0 ? BWL3 : 1.0);
double GRIP    = GRIPB * (SEXM > 0.5 ? 1.0 : 0.60) * (MUSCC / MUSC_0) * MYOQ;
double SMI     = SMIREF * (SEXM > 0.5 ? 1.0 : 0.80) * (MUSC / MUSC_0);  // L3 SMI (cm2/m2)
MGPS           = (CRP > 10.0 ? 1.0 : 0.0) + ((CRP > 10.0 && ALB < 35.0) ? 1.0 : 0.0);
TEE            = REE + PACT + DITF * INTK;
double PREE    = REEB * pow(MUSC / MUSC_0, SREEMUS);
double RQREE   = REE / PREE;
double SURV    = exp(-CHZ);
double EBALO   = INTK - TEE;
double KDOUT   = KDMUSB * clampd(0.72 * UPS + 0.28 * AUTOP, KDMIN, KDMAX);
double GDFTOT  = GDF15;                       // the state IS total GDF-15
double GDFFREE = GDF15 * (KOFFPON/KONPON) / ((KOFFPON/KONPON) + CP_PON);
double NONCON  = (MUSC - MUSCC) + LNW;

// FAACT anorexia/cachexia subscale (0-48); <= 37 defines anorexia
double FAACT = clampd(44.0 - 22.0 * ANOR - 0.85 * NAUS, 0.0, 48.0);

// Fearon 2011 stage: 0 none, 1 precachexia, 2 cachexia, 3 refractory
double SARCO = (SMI < (SEXM > 0.5 ? 43.0 : 41.0)) ? 1.0 : 0.0;
double STAGE = 0.0;
if (PCTWL > 2.0) STAGE = 1.0;
if (PCTWL > 5.0 || (PCTWL > 2.0 && BMI < 20.0) || (PCTWL > 2.0 && SARCO > 0.5)) STAGE = 2.0;
if (ECOG >= 3.0 && SATC < 0.45) STAGE = 3.0;

// Martin 2015 %WL x BMI grade (0-4)
double MGRADE = 0.0;
if (PCTWL > 2.5  && BMI < 28.0) MGRADE = 1.0;
if (PCTWL > 6.0  && BMI < 25.0) MGRADE = 2.0;
if (PCTWL > 11.0 && BMI < 22.0) MGRADE = 3.0;
if (PCTWL > 15.0 && BMI < 20.0) MGRADE = 4.0;

double CANA = CP_ANA;
double CPON = CP_PON;
double CTCZ = CP_TCZ;
double CBIM = CP_BIM;
double CMEG = CP_MEG;
double COLZ = CP_OLZ;
double CENO = CP_ENO;
double CESP = CP_ESP;
double CDEX = CP_DEX;

$CAPTURE @annotated
BW     : Body weight (kg)
LBM    : Lean body mass (kg)
BMI    : Body mass index (kg/m2)
PCTWL  : Weight loss over the preceding ~6 months (%)
GRIP   : Handgrip strength (kg)
SMI    : L3 skeletal muscle index (cm2/m2)
MGPS   : modified Glasgow Prognostic Score (0-2)
TEE    : Total energy expenditure (kcal/d)
PREE   : Predicted resting energy expenditure (kcal/d)
RQREE  : Measured/predicted REE ratio (>1.10 = hypermetabolic)
EBALO  : Net energy balance (kcal/d)
KDOUT  : Muscle protein degradation rate constant (1/d)
NONCON : Non-contractile share of measured lean gain (kg)
GDFTOT : TOTAL GDF-15, free + antibody-bound (ng/mL)
GDFFREE: FREE (GFRAL-accessible) GDF-15 (ng/mL)
FAACT  : FAACT anorexia-cachexia subscale (0-48)
STAGE  : Fearon stage (0 none, 1 pre, 2 cachexia, 3 refractory)
MGRADE : Martin %WL x BMI grade (0-4)
SARCO  : CT-defined sarcopenia flag
SURV   : Survival probability
CANA   : Anamorelin concentration (ng/mL)
CPON   : Ponsegromab concentration (nM)
CTCZ   : Tocilizumab concentration (mg/L)
CBIM   : Bimagrumab concentration (mg/L)
CMEG   : Megestrol concentration (ng/mL)
COLZ   : Olanzapine concentration (ng/mL)
CENO   : Enobosarm concentration (ng/mL)
CESP   : Espindolol concentration (ng/mL)
CDEX   : Dexamethasone concentration (ng/mL)
'

mod <- mcode_cache("cacs", code, soloc = tempdir())

CMT_NAMES <- as.list(mod)$cmt
CMT <- function(name) {
  i <- which(CMT_NAMES == name)
  if (length(i) != 1) stop("unknown compartment: ", name)
  i
}

## =====================================================================
##  SCENARIO BUILDERS
## =====================================================================
DAYS  <- 168
TGRID <- c(seq(0, 14, 0.25), seq(15, DAYS, 1))

pat_nsclc   <- list(TVOL0 = 25,  KGROW = 0.0060, FCACHEX = 1.00, CXSENS = 1.00)
pat_panc    <- list(TVOL0 = 40,  KGROW = 0.0075, FCACHEX = 1.60, CXSENS = 1.15)
pat_healthy <- list(TVOL0 = 0,   KGROW = 0.0060, FCACHEX = 1.00, CXSENS = 1.00)
pat_early   <- list(TVOL0 = 8,   KGROW = 0.0060, FCACHEX = 1.00, CXSENS = 1.00)
pat_refr    <- list(TVOL0 = 320, KGROW = 0.0070, FCACHEX = 1.35, CXSENS = 1.25)

## Registration-trial populations are NOT untreated natural history: everyone
## in them receives oncological care, symptom control and dietetic input, so
## their placebo arms are close to weight-stable over 12 weeks. pat_trial is
## calibrated to the ROMANA/POWER placebo arms (LBM about -0.5 kg and grip
## about -1.5 kg at week 12) so that model-vs-published comparisons are made
## against the right control, not against a sicker one.
pat_trial   <- list(TVOL0 = 12,  KGROW = 0.0050, FCACHEX = 0.62, CXSENS = 0.85)

## The ponsegromab phase 2 ENRICHED for elevated GDF-15 (>= 1.5 ng/mL), so its
## control arm is not the ROMANA control arm. Testing an anti-GDF-15 antibody
## against a GDF-15-low population would understate it for a reason that has
## nothing to do with the drug.
pat_gdfhi   <- list(TVOL0 = 45,  KGROW = 0.0055, FCACHEX = 1.10, CXSENS = 1.00)

## Oral doses are given into the gut depot AFTER bioavailability, so the
## F parameters are applied here rather than with F_CMT.
ev_oral <- function(cmtname, dose, F1, ii = 1, dur = DAYS, start = 0)
  ev(amt = dose * F1, cmt = CMT(cmtname), ii = ii,
     addl = ceiling(dur / ii) - 1, time = start)

ev_anamorelin    <- function(dur = DAYS, dose = 100) ev_oral("ANAD", dose, 0.35, 1,   dur)
ev_megestrol     <- function(dur = DAYS, dose = 800) ev_oral("MEGD", dose, 0.30, 1,   dur)
ev_dexamethasone <- function(dur = DAYS, dose = 4)   ev_oral("DEXD", dose, 0.80, 1,   dur)
ev_olanzapine    <- function(dur = DAYS, dose = 2.5) ev_oral("OLZD", dose, 0.60, 1,   dur)
ev_enobosarm     <- function(dur = DAYS, dose = 3)   ev_oral("ENOD", dose, 0.60, 1,   dur)
ev_espindolol    <- function(dur = DAYS, dose = 10)  ev_oral("ESPD", dose, 0.85, 0.5, dur)
ev_epa           <- function(dur = DAYS, dose = 2000)ev_oral("EPAD", dose, 0.90, 1,   dur)
ev_celecoxib     <- function(dur = DAYS, dose = 200) ev_oral("CELD", dose, 0.40, 0.5, dur)

ev_ponsegromab <- function(dur = DAYS, dose_mg = 400)
  ev(amt = dose_mg / 147000 * 1e6 * 0.65, cmt = CMT("POND"),
     ii = 28, addl = floor(dur / 28), time = 0)

ev_tocilizumab <- function(dur = DAYS, dose = 162)
  ev(amt = dose * 0.80, cmt = CMT("TCZD"), ii = 7, addl = floor(dur / 7), time = 0)

ev_bimagrumab <- function(dur = DAYS, dose_mgkg = 10, wt = 70)
  ev(amt = dose_mgkg * wt, cmt = CMT("BIMC"), ii = 28, addl = floor(dur / 28), time = 0)

combine_ev <- function(...) {
  parts <- Filter(Negate(is.null), list(...))
  if (length(parts) == 0) return(NULL)
  Reduce(function(a, b) a + b, parts)
}

run <- function(patient = pat_nsclc, params = list(), events = NULL,
                days = DAYS, tgrid = TGRID) {
  m <- param(mod, modifyList(patient, params))
  if (is.null(events)) {
    mrgsim_df(m, end = days, add = tgrid, delta = 1)
  } else {
    mrgsim_df(m, events = events, end = days, add = tgrid, delta = 1)
  }
}

## =====================================================================
##  THE 20 SCENARIOS
## =====================================================================
scenarios <- list(
  list(id = "S00", label = "Healthy control (no tumour)",
       patient = pat_healthy, params = list(), ev = NULL),
  list(id = "S01", label = "Natural history - advanced NSCLC, untreated",
       patient = pat_nsclc, params = list(), ev = NULL),
  list(id = "S02", label = "Natural history - pancreatic cancer (high cachexigenicity)",
       patient = pat_panc, params = list(), ev = NULL),
  list(id = "S03", label = "Nutrition only - ONS +600 kcal/d, protein 1.5 g/kg/d",
       patient = pat_nsclc, params = list(ONSKCAL = 600, PROTTG = 1.5), ev = NULL),
  list(id = "S04", label = "Nutrition + resistance and aerobic exercise",
       patient = pat_nsclc,
       params = list(ONSKCAL = 600, PROTTG = 1.5, EXRES = 1, EXAER = 0.6, LEUBOL = 1),
       ev = NULL),
  list(id = "S05", label = "Anamorelin 100 mg PO qd (ROMANA)",
       patient = pat_nsclc, params = list(), ev = ev_anamorelin()),
  list(id = "S06", label = "Ponsegromab 400 mg SC q4w",
       patient = pat_nsclc, params = list(), ev = ev_ponsegromab(dose_mg = 400)),
  list(id = "S07", label = "Ponsegromab 100 mg SC q4w",
       patient = pat_nsclc, params = list(), ev = ev_ponsegromab(dose_mg = 100)),
  list(id = "S08", label = "Megestrol acetate 800 mg/d",
       patient = pat_nsclc, params = list(), ev = ev_megestrol()),
  list(id = "S09", label = "Dexamethasone 4 mg/d continuous",
       patient = pat_nsclc, params = list(), ev = ev_dexamethasone()),
  list(id = "S10", label = "Olanzapine 2.5 mg PO qhs",
       patient = pat_nsclc, params = list(), ev = ev_olanzapine()),
  list(id = "S11", label = "Enobosarm 3 mg PO qd (POWER)",
       patient = pat_nsclc, params = list(), ev = ev_enobosarm()),
  list(id = "S12", label = "Bimagrumab 10 mg/kg IV q4w",
       patient = pat_nsclc, params = list(), ev = ev_bimagrumab()),
  list(id = "S13", label = "Tocilizumab 162 mg SC weekly",
       patient = pat_nsclc, params = list(), ev = ev_tocilizumab()),
  list(id = "S14", label = "Espindolol 10 mg PO bid (ACT-ONE)",
       patient = pat_nsclc, params = list(), ev = ev_espindolol()),
  list(id = "S15", label = "EPA 2 g/d + celecoxib 200 mg bid",
       patient = pat_nsclc, params = list(), ev = combine_ev(ev_epa(), ev_celecoxib())),
  list(id = "S16", label = "MULTIMODAL - ponsegromab + nutrition + exercise + EPA",
       patient = pat_nsclc,
       params = list(ONSKCAL = 600, PROTTG = 1.5, EXRES = 1, EXAER = 0.8, LEUBOL = 1),
       ev = combine_ev(ev_ponsegromab(dose_mg = 400), ev_epa())),
  list(id = "S17", label = "Effective anticancer therapy alone (RECIST responder)",
       patient = pat_nsclc, params = list(ONCEFF = 0.030, ONCSTART = 14, MYOTOX = 0.35),
       ev = NULL),
  list(id = "S18", label = "EARLY multimodal, started in precachexia",
       patient = pat_early,
       params = list(ONSKCAL = 600, PROTTG = 1.5, EXRES = 1, EXAER = 0.8, LEUBOL = 1),
       ev = combine_ev(ev_ponsegromab(dose_mg = 400), ev_epa())),
  list(id = "S19", label = "LATE multimodal, same package in refractory cachexia",
       patient = pat_refr,
       params = list(ONSKCAL = 600, PROTTG = 1.5, EXRES = 1, EXAER = 0.8, LEUBOL = 1),
       ev = combine_ev(ev_ponsegromab(dose_mg = 400), ev_epa()))
)

run_scenario <- function(s, days = DAYS) {
  out <- run(s$patient, s$params, s$ev, days = days)
  out$id <- s$id
  out
}

## =====================================================================
##  DIAGNOSTICS
## =====================================================================

## 1. Baseline stability. If TVOL0 = 0 does not hold perfectly still, then
##    every treatment effect below is contaminated by numerical drift.
check_baseline <- function(days = 365) {
  o <- run(pat_healthy, days = days, tgrid = seq(0, days, 1))
  f <- function(v) 100 * (tail(v, 1) - head(v, 1)) / head(v, 1)
  data.frame(
    state = c("BW", "MUSC", "FATM", "GRIP", "MYOQ", "CRP", "ALB",
              "INTK", "REE", "IGF1", "ECOG"),
    drift_pct = round(c(f(o$BW), f(o$MUSC), f(o$FATM), f(o$GRIP), f(o$MYOQ),
                        f(o$CRP), f(o$ALB), f(o$INTK), f(o$REE), f(o$IGF1),
                        tail(o$ECOG, 1)), 6))
}

pick <- function(d, t) d[which.min(abs(d$time - t)), ]

## 2. Week-12 endpoint grid — the readout that trials report
week12 <- function(res_list, t = 84) {
  do.call(rbind, lapply(res_list, function(d) {
    a <- pick(d, 0); b <- pick(d, t)
    data.frame(id = d$id[1],
               dBW    = round(b$BW   - a$BW, 2),
               dLBM   = round(b$LBM  - a$LBM, 2),
               dMUSC  = round(b$MUSC - a$MUSC, 2),
               dFAT   = round(b$FATM - a$FATM, 2),
               dGRIP  = round(b$GRIP - a$GRIP, 2),
               MYOQ   = round(b$MYOQ, 3),
               dFAACT = round(b$FAACT - a$FAACT, 1),
               INTK   = round(b$INTK),
               RQREE  = round(b$RQREE, 3),
               CRP    = round(b$CRP, 1),
               GDF15  = round(b$GDF15, 2),
               ECOG   = round(b$ECOG, 2),
               STAGE  = b$STAGE)
  }))
}

## 3. THE SIGNATURE DISSOCIATION.
##    Mass-only drugs must raise LBM without raising grip; only
##    quality-touching interventions may raise both. This is the model's
##    central claim, so it is tested rather than asserted.
check_dissociation <- function(res_list, t = 84) {
  g <- week12(res_list, t)
  pl <- g[g$id == "S01", ]
  ids <- c("S05", "S11", "S12", "S06", "S08", "S09", "S14", "S04", "S16")
  labs <- c("anamorelin", "enobosarm", "bimagrumab", "ponsegromab",
            "megestrol", "dexamethasone", "espindolol",
            "nutrition+exercise", "multimodal")
  do.call(rbind, Map(function(id, lab) {
    r <- g[g$id == id, ]
    dl <- r$dLBM - pl$dLBM; dg <- r$dGRIP - pl$dGRIP
    data.frame(id = id, intervention = lab,
               dLBM_vs_untreated  = round(dl, 2),
               dGRIP_vs_untreated = round(dg, 2),
               grip_per_kg_LBM    = ifelse(abs(dl) < 0.05, NA, round(dg / dl, 2)),
               pct_of_pure_mass   = ifelse(abs(dl) < 0.05, NA,
                                     round(100 * (dg / dl) / (40 / 25))),
               MYOQ = r$MYOQ)
  }, ids, labs))
}

## 4. Sensitivity of the dissociation to the two parameters that carry it:
##    QNEW (contractile fraction of unloaded hypertrophy) and SLNW (GH-driven
##    lean water). If the dissociation only exists at implausible values,
##    that should be visible rather than buried.
check_dissoc_sens <- function() {
  base <- function(pars, evx) {
    a <- run(pat_nsclc, pars, evx); b <- run(pat_nsclc, pars, NULL)
    c(dLBM  = (pick(a, 84)$LBM  - pick(a, 0)$LBM)  - (pick(b, 84)$LBM  - pick(b, 0)$LBM),
      dGRIP = (pick(a, 84)$GRIP - pick(a, 0)$GRIP) - (pick(b, 84)$GRIP - pick(b, 0)$GRIP))
  }
  g1 <- data.frame(drug = "anamorelin", knob = "SLNW", value = c(0, 1.4, 2.8, 4.2))
  g1 <- cbind(g1, t(sapply(g1$value, function(v)
    round(base(list(SLNW = v), ev_anamorelin()), 2))))
  g2 <- data.frame(drug = "enobosarm", knob = "QNEW", value = c(0.0, 0.3, 0.6, 1.0))
  g2 <- cbind(g2, t(sapply(g2$value, function(v)
    round(base(list(QNEW = v), ev_enobosarm()), 2))))
  g3 <- data.frame(drug = "bimagrumab", knob = "QNEW", value = c(0.0, 0.3, 0.6, 1.0))
  g3 <- cbind(g3, t(sapply(g3$value, function(v)
    round(base(list(QNEW = v), ev_bimagrumab()), 2))))
  rbind(g1, g2, g3)
}

## 4b. Head-to-head against published trial results, in the trial population.
validate_trials <- function() {
  d <- function(pars, evx, t = 84) {
    o <- run(pat_trial, pars, evx, days = max(t, 84))
    a <- pick(o, 0); b <- pick(o, t)
    c(BW = b$BW - a$BW, LBM = b$LBM - a$LBM, GRIP = b$GRIP - a$GRIP)
  }
  dg <- function(pars, evx, t = 84) {
    o <- run(pat_gdfhi, pars, evx, days = t)
    a <- pick(o, 0); b <- pick(o, t)
    c(BW = b$BW - a$BW, LBM = b$LBM - a$LBM, GRIP = b$GRIP - a$GRIP)
  }
  pl84  <- d(list(), NULL, 84)
  pl112 <- d(list(), NULL, 112)
  plg84 <- dg(list(), NULL, 84)
  pg    <- dg(list(), ev_ponsegromab(dose_mg = 400))
  row <- function(lab, v, ref, wk = 12) {
    base <- if (wk == 16) pl112 else pl84
    data.frame(trial = lab, week = wk,
               dBW_model   = round(v[["BW"]]   - base[["BW"]], 2),
               dLBM_model  = round(v[["LBM"]]  - base[["LBM"]], 2),
               dGRIP_model = round(v[["GRIP"]] - base[["GRIP"]], 2),
               published   = ref)
  }
  rbind(
    data.frame(trial = "PLACEBO ARM (absolute change)", week = 12,
               dBW_model = round(pl84[["BW"]], 2),
               dLBM_model = round(pl84[["LBM"]], 2),
               dGRIP_model = round(pl84[["GRIP"]], 2),
               published = "ROMANA placebo: BW +0.14, LBM -0.47, grip -1.58"),
    row("ROMANA 1/2 anamorelin 100 mg qd", d(list(), ev_anamorelin()),
        "LBM +1.1 to +1.6 kg, BW ~ +2.2 kg, GRIP no difference"),
    row("POWER 1/2 enobosarm 3 mg qd", d(list(), ev_enobosarm()),
        "LBM endpoint met (~+1.3 kg), stair-climb power MISSED"),
    data.frame(trial = "Ponsegromab 400 mg q4w (GDF-15-high population)", week = 12,
               dBW_model   = round(pg[["BW"]]   - plg84[["BW"]], 2),
               dLBM_model  = round(pg[["LBM"]]  - plg84[["LBM"]], 2),
               dGRIP_model = round(pg[["GRIP"]] - plg84[["GRIP"]], 2),
               published   = "BW +5.6% over placebo (~+2.8 kg), appetite and activity up"),
    row("Megestrol acetate 800 mg/d", d(list(), ev_megestrol()),
        "weight gain (fat and fluid), NO lean-mass or QoL benefit"),
    row("Olanzapine 2.5 mg qhs", d(list(), ev_olanzapine()),
        ">5% weight gain in ~60% vs ~9% placebo"),
    row("ACT-ONE espindolol 10 mg bid", d(list(), ev_espindolol(), 112),
        "BW +0.9 vs -1.2 kg placebo AND handgrip improved", 16),
    row("Dexamethasone 4 mg/d", d(list(), ev_dexamethasone()),
        "appetite up for 2-4 wk, then tachyphylaxis and myopathy")
  )
}

## 5. Hysteresis: identical curative therapy, started later and later.
check_hysteresis <- function(delays = c(0, 60, 120, 180, 240, 300)) {
  do.call(rbind, lapply(delays, function(d) {
    o <- run(pat_nsclc, list(ONCEFF = 0.045, ONCSTART = d), NULL,
             days = 540, tgrid = seq(0, 540, 2))
    a <- o[1, ]; nad <- o[which.min(o$MUSC), ]; b <- o[nrow(o), ]
    data.frame(start_day = d,
               nadir_MUSC = round(nad$MUSC, 2), nadir_day = nad$time,
               MUSC_d540 = round(b$MUSC, 2),
               pct_of_premorbid = round(100 * b$MUSC / a$MUSC, 1),
               regained_kg = round(b$MUSC - nad$MUSC, 2),
               permanent_deficit_kg = round(a$MUSC - b$MUSC, 2),
               SATC_d540 = round(b$SATC, 3),
               GRIP_d540 = round(b$GRIP, 1))
  }))
}

## 6. Arm decomposition — which arm does each intervention actually cut?
check_arms <- function(res_list, t = 84) {
  do.call(rbind, lapply(res_list, function(d) {
    b <- pick(d, t)
    data.frame(id = d$id[1],
               A_intake = round(b$INTK), A_anorexia = round(b$ANOR, 3),
               B_STAT3 = round(b$STAT3, 3), B_kd = round(b$KDOUT, 5),
               B_ARES = round(b$ARES, 3), B_REEratio = round(b$RQREE, 3),
               C_MYOQ = round(b$MYOQ, 3), C_MITO = round(b$MITO, 3),
               SATC = round(b$SATC, 3))
  }))
}

## 7. Tocilizumab: plasma IL-6 must RISE while CRP falls.
check_tcz_paradox <- function() {
  o <- run(pat_nsclc, list(), ev_tocilizumab())
  p <- run(pat_nsclc, list(), NULL)
  rbind(
    data.frame(arm = "untreated",  day = c(0, 28, 84),
               IL6 = round(sapply(c(0, 28, 84), function(t) pick(p, t)$IL6), 2),
               CRP = round(sapply(c(0, 28, 84), function(t) pick(p, t)$CRP), 1)),
    data.frame(arm = "tocilizumab", day = c(0, 28, 84),
               IL6 = round(sapply(c(0, 28, 84), function(t) pick(o, t)$IL6), 2),
               CRP = round(sapply(c(0, 28, 84), function(t) pick(o, t)$CRP), 1)))
}

## 8. Ponsegromab target engagement: FREE GDF-15 collapses while TOTAL rises.
check_pons_te <- function() {
  o <- run(pat_nsclc, list(), ev_ponsegromab(dose_mg = 400))
  p <- run(pat_nsclc, list(), NULL)
  data.frame(day = c(0, 7, 28, 56, 84),
             free_untreated = round(sapply(c(0, 7, 28, 56, 84), function(t) pick(p, t)$GDFFREE), 2),
             free_ponse     = round(sapply(c(0, 7, 28, 56, 84), function(t) pick(o, t)$GDFFREE), 4),
             total_ponse    = round(sapply(c(0, 7, 28, 56, 84), function(t) pick(o, t)$GDFTOT), 2))
}

## =====================================================================
##  MAIN
## =====================================================================
if (identical(environment(), globalenv())) {

  message("\n=== 1. Baseline stability (TVOL0 = 0, 365 days, % drift) ========")
  print(check_baseline(), row.names = FALSE)

  message("\n=== 2. Running the 20 scenarios =================================")
  res <- lapply(scenarios, run_scenario)
  names(res) <- sapply(scenarios, `[[`, "id")
  labs <- data.frame(id = sapply(scenarios, `[[`, "id"),
                     scenario = sapply(scenarios, `[[`, "label"))
  print(labs, row.names = FALSE)

  message("\n=== 3. Week-12 endpoint grid ====================================")
  print(week12(res), row.names = FALSE)

  message("\n=== 4. Week-24 endpoint grid ====================================")
  print(week12(res, t = 168), row.names = FALSE)

  message("\n=== 5. THE SIGNATURE DISSOCIATION (vs untreated, week 12) =======")
  print(check_dissociation(res), row.names = FALSE)

  message("\n=== 6a. Validation against published trials (trial population) ==")
  print(validate_trials(), row.names = FALSE)

  message("\n=== 6b. Sensitivity of the dissociation (anamorelin, week 12) ===")
  print(check_dissoc_sens(), row.names = FALSE)

  message("\n=== 7. Arm decomposition (week 12) ==============================")
  print(check_arms(res), row.names = FALSE)

  message("\n=== 8. Tocilizumab paradox: IL-6 UP, CRP DOWN ===================")
  print(check_tcz_paradox(), row.names = FALSE)

  message("\n=== 9. Ponsegromab target engagement ============================")
  print(check_pons_te(), row.names = FALSE)

  message("\n=== 10. Hysteresis / point of no return =========================")
  print(check_hysteresis(), row.names = FALSE)

  message("\n=== 11. Week-24 staging and survival ============================")
  print(do.call(rbind, lapply(res, function(d) {
    b <- d[nrow(d), ]
    data.frame(id = b$id, PCTWL = round(b$PCTWL, 1), BMI = round(b$BMI, 1),
               SMI = round(b$SMI, 1), MGRADE = b$MGRADE, STAGE = b$STAGE,
               MGPS = b$MGPS, ECOG = round(b$ECOG, 2),
               FAACT = round(b$FAACT, 1), SURV = round(b$SURV, 3))
  })), row.names = FALSE)

  invisible(res)
}
