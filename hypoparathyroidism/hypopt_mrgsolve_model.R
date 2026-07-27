## ============================================================================
##  Chronic Hypoparathyroidism (HypoPT) — QSP / PK-PD Model
##  ============================================================================
##  Disease  : Absent or inappropriately low parathyroid hormone. The index
##             phenotype modelled here is POST-SURGICAL HYPOPARATHYROIDISM
##             (75-80% of cases), but the same structure covers autoimmune
##             disease, genetic forms, and — through the CaSR set-point
##             parameter SETPT — autosomal dominant hypocalcaemia type 1
##             (ADH1, activating CASR mutation), which behaves very
##             differently under therapy and is therefore given its own
##             scenarios.
##
##             Physiology in one paragraph. Extracellular ionised calcium is
##             held inside a ~5% band by a single fast negative-feedback loop:
##             the parathyroid chief cell reads iCa through the calcium-sensing
##             receptor (CaSR) and secretes PTH on an inverse sigmoid. PTH then
##             does three things, and only three: (i) it turns on distal
##             tubular calcium reabsorption (TRPV5/calbindin-D28k/NCX1), (ii) it
##             turns on renal CYP27B1 so that 1,25(OH)2D rises and the gut
##             absorbs more calcium, and (iii) it raises bone remodelling so the
##             skeletal reservoir can be tapped. It also dumps phosphate by
##             internalising NPT2a. Remove PTH and all four of those go at once.
##
##             The clinically decisive point is NOT that serum calcium falls.
##             It is that the ONLY remaining lever is the gut, and the gut lever
##             is driven by exogenous calcitriol and oral calcium — neither of
##             which can be sensed and neither of which restores the renal arm.
##             So conventional therapy raises serum calcium by pushing MORE
##             calcium into a kidney that has LOST the ability to reclaim it.
##             The consequence is the central non-obvious behaviour of this
##             disease, and the one this model exists to reproduce:
##
##                 *** Filtered load rises, fractional reabsorption does not,
##                     therefore urinary calcium rises steeply with serum
##                     calcium. Patients are made eucalcaemic at the price of
##                     hypercalciuria, nephrocalcinosis and a 2-17x excess
##                     risk of CKD. PTH replacement is the only intervention
##                     that lowers urine calcium AT THE SAME serum calcium,
##                     because it restores the renal arm rather than
##                     overwhelming it. ***
##
##             That mechanism is why the trial endpoint in this field is a
##             TRIPLE COMPOSITE (normal albumin-corrected calcium AND oral
##             calcium <=600 mg/day AND active vitamin D = 0) rather than
##             serum calcium alone — a serum calcium target is trivially
##             hit by conventional therapy while the kidney is being damaged.
##             RESP3 in $TABLE is that composite.
##
##  Pharmacology modelled
##             * Oral calcium salts: gut-lumen compartment with SATURABLE
##               1,25D-dependent active transport (TRPV6/calbindin-D9k/PMCA1b)
##               plus a non-saturable paracellular path. Saturation is what
##               makes fractional absorption fall as the dose rises, and it is
##               why splitting the dose matters clinically.
##             * Calcitriol (1,25(OH)2D3) PO — enters the same 1,25D pool the
##               kidney would have filled, which is exactly its pharmacology.
##             * Alfacalcidol PO — 1-alpha-OH-D3, needs hepatic 25-hydroxylation,
##               modelled as an extra transit step (slower onset/offset).
##             * Cholecalciferol — fills the 25(OH)D substrate pool only.
##             * Hydrochlorothiazide — raises distal fractional calcium
##               reabsorption; the only conventional agent that touches urine Ca.
##             * rhPTH(1-84) (Natpara; withdrawn 2024) — BIPHASIC subcutaneous
##               absorption (fast + slow depots), which is why its profile is
##               peaky and why calcium can swing within a dosing interval.
##             * Teriparatide PTH(1-34) BID — short, tall peaks.
##             * Palopegteriparatide / TransCon PTH — an inert PEG-linked
##               prodrug depot that releases free PTH(1-34) by pH-dependent
##               linker hydrolysis with a ~60 h release half-life, producing
##               near-flat 24 h exposure. Modelled as depot -> circulating
##               prodrug -> free PTH, because the flatness IS the drug.
##             * Eneboparatide (AZP-3601) — RG-conformation-selective agonist;
##               represented by shifting the bone-resorption EC50 up relative to
##               the renal EC50 (less resorptive drive per unit renal effect).
##             * Encaleret — CaSR negative allosteric modulator. Acts on the
##               SET-POINT, not on PTH levels directly, so it only works when
##               chief cells survive (ADH1). It simultaneously de-represses
##               parathyroid secretion AND relieves CaSR-mediated inhibition of
##               tubular reabsorption, which is why it fixes hypercalciuria.
##             * IV calcium gluconate — acute post-thyroidectomy crisis.
##             * Magnesium repletion — hypomagnesaemia blocks both PTH SECRETION
##               (Galpha activation) and PTH ACTION (adenylyl cyclase). Both are
##               modelled, which is why scenario 12 is refractory until Mg is
##               given and then corrects abruptly.
##
##  Structure : 36 ODE compartments
##              calcium (gut lumen, ECF, labile bone, structural bone),
##              phosphate (ECF, bone), magnesium,
##              endogenous PTH, exogenous PTH (2 SC depots + central),
##              TransCon prodrug (depot + circulating),
##              calcitriol depot, alfacalcidol (depot + hepatic),
##              1,25(OH)2D, cholecalciferol depot, 25(OH)D,
##              thiazide (depot + central), encaleret (depot + central),
##              PTH1R signalling state (endosomal cAMP persistence),
##              osteoclast state, osteoblast state, anabolic pulse memory,
##              FGF23, nephrocalcinosis burden, GFR, BMD Z-score,
##              calcium effect site, quality of life,
##              cumulative urine Ca, cumulative urine Pi, basal-ganglia calcium
##
##  Calibration anchors (see hypopt_references.md for the citations)
##              healthy      : Ca 9.5 mg/dL, iCa 1.19 mmol/L, PTH 35 pg/mL,
##                             1,25D 40 pg/mL, Pi 3.5 mg/dL, TmP/GFR 3.1 mg/dL,
##                             urine Ca 150 mg/d, gut Ca absorption ~30%
##              untreated    : Ca 6.3-6.8 mg/dL, PTH <5 pg/mL, 1,25D ~17 pg/mL,
##                 HypoPT      Pi 4.6-4.8 mg/dL, bone turnover markers ~40-55%
##                             of normal, BMD Z-score +1.3 to +1.6
##              conventional : Ca ~8.6 mg/dL on 1000 mg elemental Ca + 0.5 ug
##                 therapy     calcitriol/day, 24 h urine Ca ~240 mg
##              PTH therapy  : same serum calcium, 24 h urine Ca ~150 mg,
##                             phosphate normalised, oral calcium withdrawn
##
##  DISCLAIMER: educational / research QSP model. Not validated for clinical
##  decision-making, dose selection or regulatory use.
## ============================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
})

hypopt_code <- '
$PROB
# Chronic Hypoparathyroidism — QSP model
- 36 compartments; Ca-Pi-Mg-PTH-vitamin D homeostasis with explicit renal arm
- Conventional therapy, PTH replacement (3 molecules), calcilytic, IV rescue
- Long-horizon renal outcome module (hypercalciuria -> nephrocalcinosis -> CKD)

$PARAM @annotated
// ------------------------- patient / anatomy -----------------------------
WT      :  70.0   : Body weight (kg)
ALB     :   4.0   : Serum albumin (g/dL)
PHART   :   7.40  : Arterial pH (alkalosis lowers ionised calcium)
VECF    :  14.0   : Extracellular fluid volume (L)
GFRC0   : 120.0   : Reference GFR (mL/min/1.73m2)
UVOL    :   2.0   : Daily urine volume (L/day)

// ------------------------- disease drivers -------------------------------
PTMASS  :   1.00  : Functional chief-cell mass (1 = normal, 0.05 = post-surgical)
SETPT   :   1.20  : CaSR set-point for PTH secretion (mmol/L ionised Ca)
HILLP   :   6.0   : Hill coefficient of the Ca-PTH sigmoid
SBAS    :   0.10  : Non-suppressible fraction of PTH secretion
SMAXP   : 934.0   : Maximal PTH secretion rate (pg/mL/h at full chief-cell mass)
KEPTH   :  10.4   : PTH elimination rate constant (1/h, t1/2 ~4 min)
ID50D   :  90.0   : 1,25D concentration suppressing PTH transcription 50% (pg/mL)
KMGP    :   1.20  : Mg for half-maximal PTH secretion (mg/dL, Hill 3)

// ------------------------- PTH1R potency ---------------------------------
EC50CA  :  35.0   : PTH EC50, distal tubular Ca reabsorption (pg/mL eq)
EC50PI  :  35.0   : PTH EC50, NPT2a internalisation (pg/mL eq)
EC50D   :  35.0   : PTH EC50, renal CYP27B1 induction (pg/mL eq)
EC50OC  :  45.0   : PTH EC50, osteoclast drive (pg/mL eq)
EC50OB  :  60.0   : PTH EC50, osteoblast drive (pg/mL eq)
POTX    :   1.00  : Potency of exogenous PTH per pg vs PTH(1-84) (2.29 for 1-34)
KSIG    :   0.30  : PTH1R signal decay (1/h); low = R0-biased, sustained endosomal cAMP

// ------------------------- renal calcium handling ------------------------
FUFCA   :   0.60  : Ultrafilterable fraction of total serum calcium
FR0     :   0.96727 : PTH-independent fractional Ca reabsorption
FRPTH   :   0.035 : Maximal PTH-dependent increment in fractional reabsorption
CASRSL  :   0.050 : TAL CaSR slope (fall in FR per mmol/L rise in ionised Ca)
CASRGN  :   1.00  : Renal CaSR gain (>1 for activating CASR mutations)
CASROFF :   0.00  : Tonic CaSR-driven loss of fractional reabsorption (ADH1)
CAION0  :   1.185 : Reference ionised calcium (mmol/L)

// ------------------------- renal phosphate handling ----------------------
FUFPI   :   0.95  : Ultrafilterable fraction of serum phosphate
TMPMAX  :   4.65  : TmP/GFR with no PTH and no FGF23 effect (mg/dL)
ETMPP   :   0.55  : Maximal fractional fall in TmP/GFR from PTH
ETMPF   :   0.15  : Maximal fractional fall in TmP/GFR from FGF23

// ------------------------- gut -------------------------------------------
DIETCA  : 1000.0  : Dietary calcium (mg/day, enters gut lumen continuously)
DIETPI  : 1200.0  : Dietary phosphate (mg/day)
DIETMG  :  300.0  : Dietary magnesium (mg/day)
KTRG    :   0.25  : Gut lumen transit rate constant (1/h)
VMCA    :  42.5   : Vmax of active (TRPV6) calcium absorption (mg/h)
KMCA    : 150.0   : Km of active calcium absorption (mg in lumen)
KDGUT   :  45.0   : 1,25D producing half-maximal active transport (pg/mL)
FPASS   :   0.10  : Paracellular fraction of transiting luminal calcium
ENDFEC  :   6.25  : Endogenous faecal calcium secretion (mg/h)
FPIP    :   0.50  : Passive fraction of dietary phosphate absorbed
FPIA    :   0.32  : 1,25D-dependent increment in phosphate absorption
ENDFPI  :   6.25  : Endogenous faecal phosphate secretion (mg/h)
BINDER  :   0.0   : Phosphate binder on (1) / off (0)
EBIND   :   0.35  : Fractional reduction of Pi absorption by binder
FMG     :   0.35  : Fractional magnesium absorption
MGSUP   :   0.0   : Elemental magnesium supplement (mg/day)
KMGCL   :   0.016889 : Renal magnesium clearance constant (1/h)
EMGP    :   0.15  : PTH-dependent reduction in magnesium excretion
FUFMG   :   0.70  : Ultrafilterable fraction of serum magnesium
MGTHR   :   2.20  : Serum Mg above which renal escape begins (mg/dL)
KMGEX   :   1.00  : Fractional renal escape of magnesium above MGTHR

// ------------------------- bone ------------------------------------------
KRES    :  16.67  : Bone resorption flux at unit osteoclast activity (mg Ca/h)
KFORM   :  16.67  : Bone formation flux at unit osteoblast activity (mg Ca/h)
KEX     :   0.050 : Labile bone -> ECF exchange rate constant (1/h)
KEXB    :   0.04511: ECF -> labile bone exchange rate constant (1/h)
CAPRAT  :   2.15  : Ca:P mass ratio of bone mineral (hydroxyapatite)
BASOC   :   0.25  : PTH-independent osteoclast activity
EMOC    :   1.70  : Maximal PTH-driven osteoclast activity
KOUTOC  :   0.01389 : Osteoclast state turnover (1/h, ~72 h)
BASOB   :   0.45  : PTH-independent osteoblast activity
EMOB    :   1.50  : Maximal PTH-driven osteoblast activity
KOUTOB  :   0.004167 : Osteoblast state turnover (1/h, ~240 h)
EMPK    :   0.55  : Anabolic gain at saturated pulse memory (intermittent PTH)
PKTHR   :   0.55  : Fractional PTH occupancy above which a peak is anabolic
KPKIN   :   2.0   : Pulse-memory accumulation rate (1/h)
KPKOUT  :   0.020 : Pulse-memory decay rate (1/h)
KZ      :   2.95e-4 : BMD Z-score gain per unit turnover suppression (1/h)
KZA     :   1.0e-4  : BMD Z-score gain from anabolic pulses (1/h)
KZO     :   1.14e-4 : BMD Z-score reversion rate (1/h)
P1NP0   :  45.0   : Reference P1NP (ug/L)
CTX0    :   0.35  : Reference CTX (ng/mL)
BSAP0   :  12.0   : Reference bone-specific ALP (ug/L)

// ------------------------- vitamin D -------------------------------------
RIN25   :   0.04125 : Baseline 25(OH)D input (ng/mL/h)
K25     :   0.001375: 25(OH)D elimination (1/h, t1/2 ~3 weeks)
KAD3    :   0.050 : Cholecalciferol absorption rate constant (1/h)
FD3     :   3.3e-4: 25(OH)D rise per absorbed IU of cholecalciferol (ng/mL per IU)
VMAX1A  :  43.6   : Maximal renal 1,25D production (pg/mL/h)
BAS1A   :   0.15  : PTH-independent fraction of CYP27B1 activity
K25M    :  15.0   : 25(OH)D Km for 1-alpha-hydroxylation (ng/mL)
KPI1A   :  12.0   : Phosphate inhibition constant for CYP27B1 (mg/dL)
KF1A    :   1.0   : FGF23 inhibition constant for CYP27B1 (relative units)
KO125   :   0.1155: Baseline 1,25D elimination (1/h, t1/2 ~6 h)
IND24   :   1.00  : Maximal CYP24A1 autoinduction factor
K24     :  60.0   : 1,25D for half-maximal CYP24A1 autoinduction (pg/mL)
VD125   :   5.0   : Apparent volume of the 1,25D pool (L)
KACTR   :   1.50  : Calcitriol absorption rate constant (1/h)
FCTR    :   0.70  : Calcitriol oral bioavailability
KAALF   :   1.00  : Alfacalcidol absorption rate constant (1/h)
FALF    :   0.75  : Alfacalcidol oral bioavailability
K25ALF  :   0.10  : Hepatic 25-hydroxylation of alfacalcidol (1/h)

// ------------------------- FGF23 -----------------------------------------
KINF    :   0.70  : FGF23 production (relative units/h)
KOUTF   :   0.70  : FGF23 elimination (1/h, t1/2 ~1 h)
AP      :   1.50  : FGF23 sensitivity to fractional phosphate deviation
AD      :   0.80  : FGF23 sensitivity to fractional 1,25D deviation
PI0     :   3.50  : Reference serum phosphate (mg/dL)
D1250   :  40.0   : Reference 1,25D (pg/mL)
CAT0    :   9.50  : Reference total serum calcium (mg/dL)
FION0   :   0.50  : Ionised fraction of total calcium at pH 7.40

// ------------------------- exogenous PTH PK ------------------------------
FPTHSC  :   0.53  : SC bioavailability of injected PTH
FFAST   :   0.40  : Fraction of the SC dose in the fast depot
KA1F    :   4.00  : Fast SC absorption rate constant (1/h)
KA1S    :   0.20  : Slow SC absorption rate constant (1/h)
VPTHX   : 130.0   : Apparent volume of distribution of exogenous PTH (L)
CLPTHX  :  78.0   : Apparent clearance of exogenous PTH (L/h)
KATC    :   0.060 : TransCon PTH depot absorption rate constant (1/h)
FTC     :   0.70  : TransCon PTH SC bioavailability (prodrug)
KRELTC  :   0.0116: TransCon linker hydrolysis / PTH release (1/h, t1/2 ~60 h)
KCLTC   :   0.0015: Clearance of intact TransCon prodrug (1/h)

// ------------------------- thiazide --------------------------------------
KATHZ   :   1.20  : Hydrochlorothiazide absorption rate constant (1/h)
FTHZ    :   0.70  : Hydrochlorothiazide bioavailability
VTHZ    : 200.0   : Hydrochlorothiazide V/F (L)
CLTHZ   :  20.0   : Hydrochlorothiazide CL/F (L/h)
EMXTHZ  :   0.012 : Maximal thiazide increment in fractional Ca reabsorption
EC50THZ :  25.0   : Thiazide EC50 (ng/mL)

// ------------------------- encaleret (CaSR NAM) --------------------------
KAENC   :   1.00  : Encaleret absorption rate constant (1/h)
FENC    :   0.60  : Encaleret bioavailability
VENC    : 300.0   : Encaleret V/F (L)
CLENC   :  30.0   : Encaleret CL/F (L/h)
EMXENC  :   0.40  : Maximal rightward shift of the CaSR set-point
EMXENCK :   0.022 : Maximal renal (TAL) increment in fractional Ca reabsorption
EC50ENC : 100.0   : Encaleret EC50 (ng/mL)

// ------------------------- renal outcome module --------------------------
UCTHR   : 100.0   : Urine calcium concentration threshold (mg/L)
KNC     :   5.0e-5: Nephrocalcinosis accretion rate (units/h per unit risk)
KNCR    :   1.0e-5: Nephrocalcinosis resolution rate (1/h)
KGFR0   :   9.1e-5: Age-related GFR decline (mL/min/1.73m2 per h)
KGFRNC  :   3.0e-4: Additional GFR decline per unit nephrocalcinosis
KBG     :   2.0e-6: Basal-ganglia calcification accretion (units per mg2/dL2/h)

// ------------------------- symptoms / QoL --------------------------------
KEO     :   0.15  : Equilibration rate to the calcium effect site (1/h)
SMAX    : 100.0   : Maximal symptom score
SC50    :   0.98  : Effect-site ionised Ca at half-maximal symptoms (mmol/L)
SGAM    :   0.08  : Steepness of the symptom curve (mmol/L)
KQ      :   0.002 : Quality-of-life adaptation rate (1/h)
QTC0    : 400.0   : QTc at reference calcium (ms)
KQT     :  14.0   : QTc prolongation per mg/dL fall in total calcium (ms)

// ------------------------- regimen bookkeeping ---------------------------
DCASUP  :   0.0   : Prescribed elemental oral calcium (mg/day, for composite)
DCTR    :   0.0   : Prescribed active vitamin D (ug/day, for composite)

$CMT @annotated
GUTCA  : Calcium in gut lumen (mg)
CAE    : Extracellular calcium (mg)
CABONE : Rapidly exchangeable bone calcium (mg)
CAMIN  : Structural bone mineral calcium (mg)
PIE    : Extracellular phosphate (mg)
PIBONE : Bone mineral phosphate (mg)
MGE    : Extracellular magnesium (mg)
PTHE   : Endogenous plasma PTH (pg/mL)
PTHD1  : Exogenous PTH, fast SC depot (ug)
PTHD2  : Exogenous PTH, slow SC depot (ug)
PTHX   : Exogenous PTH, central (ug)
TCDEP  : TransCon PTH prodrug, SC depot (ug PTH equivalents)
TCPRO  : TransCon PTH prodrug, circulating (ug PTH equivalents)
CTRDEP : Calcitriol in gut (ug)
ALFDEP : Alfacalcidol in gut (ug)
ALFC   : Alfacalcidol awaiting hepatic 25-hydroxylation (ug)
D125   : 1,25-dihydroxyvitamin D (pg/mL)
D3DEP  : Cholecalciferol in gut (IU)
D25    : 25-hydroxyvitamin D (ng/mL)
THZDEP : Hydrochlorothiazide in gut (mg)
THZC   : Hydrochlorothiazide, central (mg)
ENCDEP : Encaleret in gut (mg)
ENCC   : Encaleret, central (mg)
OC     : Osteoclast activity (relative to normal)
OB     : Osteoblast activity (relative to normal)
PSIG   : PTH1R signalling state (pg/mL equivalents)
PKS    : Anabolic pulse memory (0-1, saturating)
FGF23  : FGF23 (relative units, 1 = normal)
NC     : Nephrocalcinosis burden (arbitrary units, 1 = imaging threshold)
GFRC   : Glomerular filtration rate (mL/min/1.73m2)
BMDZ   : Lumbar spine BMD Z-score
CAEFF  : Ionised calcium at the neuromuscular effect site (mmol/L)
QOL    : Quality-of-life score (0-100, higher better)
UCACUM : Cumulative urinary calcium (mg)
UPICUM : Cumulative urinary phosphate (mg)
BGCALC : Basal-ganglia calcification burden (arbitrary units)

$GLOBAL
// --------------------------------------------------------------------------
// Algebraic layer. Written as macros so that $ODE and $TABLE see exactly the
// same expressions — a QSP model that reports one thing and integrates
// another is worse than no model at all.
// --------------------------------------------------------------------------
#define GFRLH    (GFRC*0.06)
#define CATOT    (CAE/(VECF*10.0))
#define PICONC   (PIE/(VECF*10.0))
#define MGCONC   (MGE/(VECF*10.0))
#define FIONX    (FION0*(1.0 - 0.40*(PHART - 7.40)))
#define CAIONX   (CATOT*FIONX/4.008)
#define CACORRX  (CATOT + 0.8*(4.0 - ALB))

// ---- drug concentrations
#define CPTHX    (PTHX*1000.0/VPTHX)
#define CTHZ     (THZC/VTHZ*1000.0)
#define CENC     (ENCC/VENC*1000.0)

// ---- magnesium modulation of PTH ACTION (adenylyl cyclase coupling)
#define MGRES    (fmin(1.10, (MGCONC*MGCONC/(0.64 + MGCONC*MGCONC))/0.86207))
#define MGSEC    (pow(MGCONC,3.0)/(pow(KMGP,3.0) + pow(MGCONC,3.0)))

// ---- total bioactive PTH signal
#define PTHTOT   (PTHE + POTX*CPTHX)
#define PTHEFF   (PTHTOT*MGRES)

// Downstream organ effects read the SIGNALLING state, not the plasma
// concentration. PTH1R internalises and keeps producing cAMP from the
// endosome, so the pharmacodynamic duration outlasts the pharmacokinetic
// one — markedly so for R0-biased ligands such as PTH(1-84) and LA-PTH,
// barely at all for RG-selective ones such as eneboparatide. KSIG carries
// that distinction. Pulse memory alone reads the instantaneous ligand,
// because it is the PEAK that opens the anabolic window.
#define EPCA     (PSIG/(EC50CA + PSIG))
#define EPPI     (PSIG/(EC50PI + PSIG))
#define EPD      (PSIG/(EC50D  + PSIG))
#define EPOC     (PSIG/(EC50OC + PSIG))
#define EPOB     (PSIG/(EC50OB + PSIG))
#define EPOBP    (PTHEFF/(EC50OB + PTHEFF))

// ---- CaSR set-point, shifted rightwards by a calcilytic
#define SETPTE   (SETPT*(1.0 + EMXENC*CENC/(EC50ENC + CENC)))
#define SECFRAC  (SBAS + (1.0-SBAS)/(1.0 + pow(CAIONX/SETPTE, HILLP)))
#define VDFB     (1.0/(1.0 + pow(D125/ID50D, 2.0)))
#define PTHSEC   (PTMASS*SMAXP*SECFRAC*MGSEC*VDFB)

// ---- renal calcium
#define FRTHZM   (EMXTHZ*CTHZ/(EC50THZ + CTHZ))
#define FRENCM   (EMXENCK*CENC/(EC50ENC + CENC))
#define FRRAW    (FR0 + FRPTH*EPCA + FRTHZM + FRENCM - CASROFF - CASRSL*CASRGN*(CAIONX - CAION0))
#define FRCA     (FRRAW > 0.997 ? 0.997 : (FRRAW < 0.900 ? 0.900 : FRRAW))
#define CAFILT   (GFRLH*CATOT*10.0*FUFCA)
#define UCARATE  (CAFILT*(1.0 - FRCA))

// ---- renal phosphate
#define EFGFTM   (FGF23/(1.0 + FGF23))
#define TMPGFR   (TMPMAX*(1.0 - ETMPP*EPPI)*(1.0 - ETMPF*EFGFTM))
#define PIFILT   (GFRLH*PICONC*10.0*FUFPI)
#define PIREABS  (fmin(PIFILT, TMPGFR*GFRLH*10.0*FUFPI))
#define UPIRATE  (PIFILT - PIREABS)

// ---- gut
#define DFAC     (D125/(KDGUT + D125))
#define ACTCA    (VMCA*DFAC*GUTCA/(KMCA + GUTCA))
#define PASCA    (FPASS*KTRG*GUTCA)
#define ABSCA    (ACTCA + PASCA)
#define ABSPI    ((FPIP + FPIA*DFAC)*(1.0 - EBIND*BINDER)*DIETPI/24.0)
#define ABSMG    (FMG*(DIETMG + MGSUP)/24.0)
// Magnesium leaves by two routes: an avid, PTH-modulated conservation path
// that sets the level during depletion, and a threshold escape that dumps
// any surplus. Without the escape term, oral repletion drives serum Mg to
// frankly toxic levels; without the conservation path, dietary depletion
// cannot produce hypomagnesaemia at all.
#define UMGRATE  (KMGCL*MGE*(1.0 - EMGP*EPCA) + KMGEX*fmax(0.0, MGCONC - MGTHR)*GFRLH*10.0*FUFMG)

// ---- bone
#define RESCA    (KRES*OC)
#define FORMCA   (KFORM*OB*fmin(1.30, sqrt(fmax(0.05, CATOT/CAT0))))
#define FB2E     (KEX*CABONE)
#define FE2B     (KEXB*CAE)

// ---- vitamin D
#define D25FAC   (D25/(K25M + D25))
#define PIINH1A  (1.0/(1.0 + PICONC/KPI1A))
#define FGFINH1A (1.0/(1.0 + FGF23/KF1A))
#define PROD125  (VMAX1A*(BAS1A + (1.0-BAS1A)*EPD)*D25FAC*(GFRC/GFRC0)*PIINH1A*FGFINH1A)
#define KEL125   (KO125*(1.0 + IND24*D125/(K24 + D125)))

// ---- renal risk
#define CAXP     (CATOT*PICONC)
#define UCADAY   (UCARATE*24.0)
#define UCACONC  (UCADAY/UVOL)
#define NCRISK   ((UCACONC/UCTHR)*(1.0 + 0.5*fmax(0.0, (CAXP - 40.0)/15.0)))

// ---- symptoms
#define SYMV     (SMAX/(1.0 + exp((CAEFF - SC50)/SGAM)))

$MAIN
// Healthy steady state as the default initial condition. Disease and therapy
// scenarios are started from an equilibrated state produced in R by
// hypopt_baseline(), which runs this model forward with the disease
// parameters until the derivatives vanish.
CAE_0    = CAT0*VECF*10.0;      // 1330 mg  -> 9.5 mg/dL
GUTCA_0  = 129.6;               // steady luminal calcium on a 1000 mg/day diet
CABONE_0 = 1200.0;
CAMIN_0  = 1.0e6;
PIE_0    = PI0*VECF*10.0;       // 490 mg   -> 3.5 mg/dL
PIBONE_0 = 5.0e5;
MGE_0    = 2.0*VECF*10.0;       // 280 mg   -> 2.0 mg/dL
PTHE_0   = 35.0;
PSIG_0   = 35.0;
D125_0   = 40.0;
D25_0    = 30.0;
OC_0     = 1.0;
OB_0     = 1.0;
FGF23_0  = 1.0;
GFRC_0   = GFRC0;
CAEFF_0  = CAT0*FION0/4.008;
QOL_0    = 92.85;

// SC bioavailability is applied to the two parallel PTH depots. The dosing
// helper puts the SAME nominal amount in both; F splits it.
F_PTHD1  = FPTHSC*FFAST;
F_PTHD2  = FPTHSC*(1.0 - FFAST);
F_TCDEP  = FTC;
F_CTRDEP = FCTR;
F_ALFDEP = FALF;
F_THZDEP = FTHZ;
F_ENCDEP = FENC;
F_MGE    = FMG;   // oral magnesium salts are dosed straight into the Mg pool

$ODE
// ---------------------------------------------------------------- gut lumen
dxdt_GUTCA = DIETCA/24.0 - KTRG*GUTCA - ACTCA;

// ------------------------------------------------------ extracellular calcium
dxdt_CAE = ABSCA - ENDFEC - UCARATE
           + RESCA - FORMCA
           + FB2E - FE2B;

// --------------------------------------------------------------------- bone
dxdt_CABONE = FE2B - FB2E;
dxdt_CAMIN  = FORMCA - RESCA;

// ---------------------------------------------------------------- phosphate
dxdt_PIE    = ABSPI - ENDFPI - UPIRATE + RESCA/CAPRAT - FORMCA/CAPRAT;
dxdt_PIBONE = FORMCA/CAPRAT - RESCA/CAPRAT;

// ---------------------------------------------------------------- magnesium
dxdt_MGE = ABSMG - UMGRATE;

// ------------------------------------------------------------ endogenous PTH
dxdt_PTHE = PTHSEC - KEPTH*PTHE;

// ------------------------------------------------------- exogenous PTH (PK)
dxdt_PTHD1 = -KA1F*PTHD1;
dxdt_PTHD2 = -KA1S*PTHD2;
dxdt_TCDEP = -KATC*TCDEP;
dxdt_TCPRO =  KATC*TCDEP - KRELTC*TCPRO - KCLTC*TCPRO;
dxdt_PTHX  =  KA1F*PTHD1 + KA1S*PTHD2 + KRELTC*TCPRO - (CLPTHX/VPTHX)*PTHX;

// -------------------------------------------------------------- vitamin D
dxdt_CTRDEP = -KACTR*CTRDEP;
dxdt_ALFDEP = -KAALF*ALFDEP;
dxdt_ALFC   =  KAALF*ALFDEP - K25ALF*ALFC;
dxdt_D3DEP  = -KAD3*D3DEP;
dxdt_D25    =  RIN25 + FD3*KAD3*D3DEP - K25*D25;
dxdt_D125   =  PROD125 - KEL125*D125
               + (KACTR*CTRDEP + K25ALF*ALFC)*1000.0/VD125;

// ------------------------------------------------------------------ thiazide
dxdt_THZDEP = -KATHZ*THZDEP;
dxdt_THZC   =  KATHZ*THZDEP - (CLTHZ/VTHZ)*THZC;

// ----------------------------------------------------------------- encaleret
dxdt_ENCDEP = -KAENC*ENCDEP;
dxdt_ENCC   =  KAENC*ENCDEP - (CLENC/VENC)*ENCC;

// ------------------------------------------------------------- bone cells
// Osteoclasts follow the time-averaged PTH signal (slow state), osteoblasts
// follow PTH plus a pulse-memory term. That asymmetry is the anabolic window:
// short tall peaks build bone, flat exposure does not.
dxdt_PSIG = KSIG*(PTHEFF - PSIG);
dxdt_PKS  = KPKIN*fmax(0.0, EPOBP - PKTHR)*(1.0 - PKS) - KPKOUT*PKS;
dxdt_OC  = KOUTOC*(BASOC + EMOC*EPOC - OC);
dxdt_OB  = KOUTOB*(BASOB + EMOB*EPOB + EMPK*PKS - OB);

// ---------------------------------------------------------------- FGF23
dxdt_FGF23 = KINF*fmax(0.20, 1.0 + AP*(PICONC/PI0 - 1.0) + AD*(D125/D1250 - 1.0))
             - KOUTF*FGF23;

// ------------------------------------------------------- renal outcome arm
dxdt_NC   = KNC*fmax(0.0, NCRISK - 1.0) - KNCR*NC;
dxdt_GFRC = (GFRC > 10.0 ? -(KGFR0 + KGFRNC*NC) : 0.0);

// -------------------------------------------------------------------- bone
dxdt_BMDZ = KZ*(1.0 - OC) + KZA*PKS - KZO*BMDZ;

// ------------------------------------------------------- symptoms and QoL
dxdt_CAEFF = KEO*(CAIONX - CAEFF);
dxdt_QOL   = KQ*((100.0 - SYMV) - QOL);

// ----------------------------------------------------------- accumulators
dxdt_UCACUM = UCARATE;
dxdt_UPICUM = UPIRATE;
dxdt_BGCALC = KBG*fmax(0.0, CAXP - 45.0);

$TABLE
double CATOTAL = CATOT;
double CAIONMM = CAIONX;
double CACORRD = CACORRX;
double PISER   = PICONC;
double MGSER   = MGCONC;
double PTHPG   = PTHTOT;
double PTHENDO = PTHE;
double PTHEXOG = CPTHX;
double PTHSIG  = PSIG;
double D125PG  = D125;
double D25NG   = D25;
double FGF23RU = FGF23;
double TMPPERG = TMPGFR;
double FRCAPCT = FRCA*100.0;
double FEXCA   = (1.0 - FRCA)*100.0;
double UCA24   = UCADAY;
double UPI24   = UPIRATE*24.0;
double UCAMGKG = UCADAY/WT;
double CAABS   = ABSCA*24.0;
double FABSPCT = (DIETCA > 0.0 ? 100.0*ABSCA*24.0/DIETCA : 0.0);
double CAPROD  = CAXP;
double EGFR    = GFRC;
double NCBURD  = NC;
double BMDZS   = BMDZ;
double P1NP    = P1NP0*OB;
double CTX     = CTX0*OC;
double BSAP    = BSAP0*OB;
double TURNOV  = 0.5*(OB + OC);
double SYMSCOR = SYMV;
double QTCMS   = QTC0 + KQT*(CAT0 - CATOT);
double QOLSC   = QOL;
double BGCALCI = BGCALC;

// Tetany-risk flag: effect-site ionised calcium below 1.00 mmol/L
double TETANY  = (CAEFF < 1.00 ? 1.0 : 0.0);

// PaTHway/PARALLAX-style triple composite response
double INDEP   = ((DCASUP <= 600.0) && (DCTR <= 0.0)) ? 1.0 : 0.0;
double RESP3   = ((CACORRD >= 8.0) && (CACORRD <= 10.6) && (INDEP > 0.5)) ? 1.0 : 0.0;

$CAPTURE @annotated
CATOTAL : Total serum calcium (mg/dL)
CACORRD : Albumin-corrected serum calcium (mg/dL)
CAIONMM : Ionised calcium (mmol/L)
PISER   : Serum phosphate (mg/dL)
MGSER   : Serum magnesium (mg/dL)
PTHPG   : Total bioactive PTH signal (pg/mL PTH(1-84) equivalents)
PTHENDO : Endogenous PTH (pg/mL)
PTHEXOG : Exogenous PTH concentration (pg/mL)
PTHSIG  : PTH1R signalling state (pg/mL equivalents)
D125PG  : 1,25-dihydroxyvitamin D (pg/mL)
D25NG   : 25-hydroxyvitamin D (ng/mL)
FGF23RU : FGF23 (relative units)
TMPPERG : TmP/GFR (mg/dL)
FRCAPCT : Fractional tubular calcium reabsorption (%)
FEXCA   : Fractional excretion of calcium (%)
UCA24   : 24-hour urinary calcium (mg/day)
UPI24   : 24-hour urinary phosphate (mg/day)
UCAMGKG : Urinary calcium (mg/kg/day)
CAABS   : Absorbed calcium (mg/day)
FABSPCT : Fractional calcium absorption (% of dietary intake)
CAPROD  : Calcium-phosphate product (mg2/dL2)
EGFR    : eGFR (mL/min/1.73m2)
NCBURD  : Nephrocalcinosis burden (units)
BMDZS   : Lumbar spine BMD Z-score
P1NP    : P1NP (ug/L)
CTX     : CTX (ng/mL)
BSAP    : Bone-specific alkaline phosphatase (ug/L)
TURNOV  : Composite bone turnover index (1 = normal)
SYMSCOR : Neuromuscular symptom score (0-100)
QTCMS   : QTc interval (ms)
QOLSC   : Quality-of-life score (0-100)
BGCALCI : Basal-ganglia calcification burden (units)
TETANY  : Tetany-risk flag (1 = effect-site iCa < 1.00 mmol/L)
RESP3   : Triple composite response (1 = met)
'

## ---------------------------------------------------------------------------
##  Build
## ---------------------------------------------------------------------------
HYPOPT_build <- function() {
  mrgsolve::mcode_cache("hypopt", hypopt_code)
}

## ---------------------------------------------------------------------------
##  Equilibration.
##
##  A disease state is not an initial condition you can write down — it is
##  where the system settles once the parameter is changed. This runs the model
##  forward with the disease parameters and no therapy until the derivatives
##  vanish, then hands back the state vector to be used as `init` for the
##  therapy simulation. Cumulative compartments are reset so that on-treatment
##  urinary load starts from zero.
## ---------------------------------------------------------------------------
hypopt_baseline <- function(mod = HYPOPT_build(), pars = list(), days = 400) {
  m <- mod
  if (length(pars)) m <- mrgsolve::update(m, param = pars)
  out  <- mrgsolve::mrgsim_df(m, end = days*24, delta = 24, hmax = 2)
  last <- out[nrow(out), , drop = FALSE]
  cmts <- names(mrgsolve::init(mod))
  iv   <- as.numeric(last[, cmts]); names(iv) <- cmts
  iv[c("UCACUM", "UPICUM")] <- 0
  iv
}

## ---------------------------------------------------------------------------
##  Number of additional doses to carry an every-`ii`-hours regimen for `days`
## ---------------------------------------------------------------------------
hypopt_addl <- function(days, ii) ceiling(days*24/ii) - 1

## ---------------------------------------------------------------------------
##  Dosing helpers
## ---------------------------------------------------------------------------

#' Oral elemental calcium, mg per dose
hypopt_calcium <- function(mg_per_dose, ii = 8, days = 180, start = 0) {
  mrgsolve::ev(time = start, amt = mg_per_dose, cmt = "GUTCA",
               ii = ii, addl = hypopt_addl(days, ii))
}

#' Oral calcitriol, ug per dose
hypopt_calcitriol <- function(ug_per_dose, ii = 12, days = 180, start = 0) {
  mrgsolve::ev(time = start, amt = ug_per_dose, cmt = "CTRDEP",
               ii = ii, addl = hypopt_addl(days, ii))
}

#' Oral alfacalcidol, ug per dose
hypopt_alfacalcidol <- function(ug_per_dose, ii = 24, days = 180, start = 0) {
  mrgsolve::ev(time = start, amt = ug_per_dose, cmt = "ALFDEP",
               ii = ii, addl = hypopt_addl(days, ii))
}

#' Cholecalciferol, IU per dose
hypopt_vitd3 <- function(iu_per_dose, ii = 24, days = 180, start = 0) {
  mrgsolve::ev(time = start, amt = iu_per_dose, cmt = "D3DEP",
               ii = ii, addl = hypopt_addl(days, ii))
}

#' Hydrochlorothiazide, mg per dose
hypopt_thiazide <- function(mg_per_dose = 25, ii = 12, days = 180, start = 0) {
  mrgsolve::ev(time = start, amt = mg_per_dose, cmt = "THZDEP",
               ii = ii, addl = hypopt_addl(days, ii))
}

#' Oral elemental magnesium, mg per dose (absorbed fraction FMG applied)
hypopt_magnesium <- function(mg_per_dose = 400, ii = 8, days = 60, start = 0) {
  mrgsolve::ev(time = start, amt = mg_per_dose, cmt = "MGE",
               ii = ii, addl = hypopt_addl(days, ii))
}

#' Encaleret, mg per dose
hypopt_encaleret <- function(mg_per_dose = 60, ii = 12, days = 180, start = 0) {
  mrgsolve::ev(time = start, amt = mg_per_dose, cmt = "ENCDEP",
               ii = ii, addl = hypopt_addl(days, ii))
}

#' Injected PTH (rhPTH(1-84) or teriparatide). The same nominal amount goes to
#' both SC depots; F_PTHD1 / F_PTHD2 split it into the fast and slow fractions,
#' which is what reproduces the biphasic subcutaneous profile.
hypopt_pth_sc <- function(ug_per_dose, ii = 24, days = 180, start = 0) {
  a <- mrgsolve::ev(time = start, amt = ug_per_dose, cmt = "PTHD1",
                    ii = ii, addl = hypopt_addl(days, ii))
  b <- mrgsolve::ev(time = start, amt = ug_per_dose, cmt = "PTHD2",
                    ii = ii, addl = hypopt_addl(days, ii))
  mrgsolve::as.ev(rbind(as.data.frame(a), as.data.frame(b)))
}

#' Palopegteriparatide (TransCon PTH), ug PTH(1-34) equivalents per dose
hypopt_transcon <- function(ug_per_dose = 18, ii = 24, days = 180, start = 0) {
  mrgsolve::ev(time = start, amt = ug_per_dose, cmt = "TCDEP",
               ii = ii, addl = hypopt_addl(days, ii))
}

#' IV calcium gluconate. 1 g calcium gluconate = 93 mg elemental calcium.
#' `g_gluconate` is infused over `hours`.
hypopt_iv_calcium <- function(g_gluconate = 2, hours = 1, start = 0,
                              ii = 0, addl = 0) {
  amt <- g_gluconate*93
  mrgsolve::ev(time = start, amt = amt, cmt = "CAE",
               rate = amt/hours, ii = ii, addl = addl)
}

## ---------------------------------------------------------------------------
##  Parameter sets for the exogenous PTH molecules.
##
##  rhPTH(1-84)      : biphasic SC absorption, F 0.53, potency 1.0 per pg
##  teriparatide     : fast single-phase SC, F 0.95, potency 2.29 per pg
##                     (PTH(1-34) MW 4118 vs PTH(1-84) MW 9425)
##  palopegteriparatide: prodrug depot; the free PTH(1-34) it liberates uses
##                     the teriparatide disposition parameters
## ---------------------------------------------------------------------------
HYPOPT_pth_params <- list(
  rhpth84 = list(FPTHSC = 0.53, FFAST = 0.40, KA1F = 4.0, KA1S = 0.20,
                 VPTHX = 130, CLPTHX = 78, POTX = 1.00, KSIG = 0.15),
  terip   = list(FPTHSC = 0.95, FFAST = 1.00, KA1F = 6.0, KA1S = 0.20,
                 VPTHX = 100, CLPTHX = 62, POTX = 2.29, KSIG = 0.30),
  transcon= list(VPTHX = 100, CLPTHX = 62, POTX = 2.29, KSIG = 0.30),
  eneb    = list(FPTHSC = 0.80, FFAST = 1.00, KA1F = 3.0, KA1S = 0.20,
                 VPTHX = 100, CLPTHX = 62, POTX = 2.29,
                 KSIG = 0.45, EC50OC = 160)  # RG-selective: renal effect kept, resorption blunted
)

## ---------------------------------------------------------------------------
##  Fifteen prebuilt scenarios
## ---------------------------------------------------------------------------
#' @param mod   model object from HYPOPT_build()
#' @param which integer vector of scenario numbers (default all 15)
#' @return data.frame with one row per output time and a `scenario` column
HYPOPT_simulate_scenarios <- function(mod = HYPOPT_build(), which = 1:15) {

  D  <- 180            # standard therapy horizon (days)
  H  <- D*24
  LT <- 10*365*24      # long-horizon renal outcome (10 years)

  ## disease parameter blocks
  P_HEALTHY <- list(PTMASS = 1.00, SETPT = 1.20)
  P_SURG    <- list(PTMASS = 0.05, SETPT = 1.20)   # post-surgical HypoPT
  ## ADH1: the mutation activates CaSR in BOTH tissues — the parathyroid
  ## set-point shifts left AND the thick ascending limb tonically wastes
  ## calcium, which is why urine calcium is high for the serum calcium.
  P_ADH1    <- list(PTMASS = 1.00, SETPT = 0.78, CASRGN = 1.50, CASROFF = 0.012)
  P_MGLOW   <- list(PTMASS = 1.00, SETPT = 1.20, DIETMG = 120) # severe Mg depletion

  run <- function(label, pars = list(), base_pars = NULL, events = NULL,
                  end = H, delta = 2, init = NULL) {
    m <- mrgsolve::update(mod, param = utils::modifyList(pars, list()))
    if (length(pars)) m <- mrgsolve::update(mod, param = pars)
    if (is.null(init)) {
      bp   <- if (is.null(base_pars)) pars else base_pars
      init <- hypopt_baseline(mod, bp)
    }
    m <- mrgsolve::update(m, init = as.list(init))
    d <- if (is.null(events)) {
      mrgsolve::mrgsim_df(m, end = end, delta = delta, hmax = 2)
    } else {
      mrgsolve::mrgsim_df(m, events = events, end = end, delta = delta, hmax = 2)
    }
    d$scenario <- label
    d
  }

  S <- list()

  ## 1 — healthy control. The reference against which everything else is read.
  S[[1]] <- function() run("01 정상 대조군 (Healthy control)", P_HEALTHY)

  ## 2 — post-surgical HypoPT, untreated. Serum Ca settles ~6.3-6.8 mg/dL,
  ##     phosphate ~4.7, 1,25D ~17, turnover markers ~half normal, BMD Z +1.5.
  S[[2]] <- function() run("02 수술후 HypoPT — 무치료 (Untreated)", P_SURG)

  ## 3 — conventional therapy. Calcium carbonate 500 mg elemental TID with
  ##     meals + calcitriol 0.25 ug BID. Serum Ca reaches target; urine Ca does
  ##     not. This is the scenario the whole model exists to contrast against.
  S[[3]] <- function() run(
    "03 기존 치료 (Ca 1000 mg/일 + 칼시트리올 0.5 ug/일)",
    utils::modifyList(P_SURG, list(DCASUP = 1000, DCTR = 0.5)),
    base_pars = P_SURG,
    events = c(hypopt_calcium(500, ii = 12, days = D),
               hypopt_calcitriol(0.25, ii = 12, days = D),
               hypopt_vitd3(1000, ii = 24, days = D)))

  ## 4 — conventional therapy + hydrochlorothiazide 25 mg BID on a low-sodium
  ##     diet. The only conventional lever on urine calcium; expect roughly a
  ##     25-30% fall in 24 h urine Ca at a slightly higher serum calcium.
  S[[4]] <- function() run(
    "04 기존 치료 + 티아지드 (칼슘 절감: Ca 500 mg/일 + HCTZ 25 mg BID)",
    utils::modifyList(P_SURG, list(DCASUP = 500, DCTR = 0.5)),
    base_pars = P_SURG,
    events = c(hypopt_calcium(500, ii = 24, days = D),
               hypopt_calcitriol(0.25, ii = 12, days = D),
               hypopt_vitd3(1000, ii = 24, days = D),
               hypopt_thiazide(25, ii = 12, days = D)))

  ## 5 — over-treatment. Ca 1000 mg TID + calcitriol 0.5 ug BID. Serum calcium
  ##     drifts high-normal, Ca x P crosses 45, urine Ca and nephrocalcinosis
  ##     risk climb sharply. The classic iatrogenic trajectory.
  S[[5]] <- function() run(
    "05 과치료 (Ca 3000 mg/일 + 칼시트리올 1.0 ug/일)",
    utils::modifyList(P_SURG, list(DCASUP = 3000, DCTR = 1.0)),
    base_pars = P_SURG,
    events = c(hypopt_calcium(1000, ii = 8, days = D),
               hypopt_calcitriol(0.5, ii = 12, days = D)))

  ## 6 — rhPTH(1-84) 50 ug SC daily added to halved conventional therapy, then
  ##     conventional therapy withdrawn at week 4. Biphasic absorption gives a
  ##     visible within-day calcium swing — the known weakness of once-daily
  ##     PTH(1-84) and the reason for the pump and prodrug programmes.
  S[[6]] <- function() run(
    "06 rhPTH(1-84) 50→100 ug SC QD + 기존치료 감량",
    utils::modifyList(c(P_SURG, HYPOPT_pth_params$rhpth84),
                      list(DCASUP = 0, DCTR = 0)),
    base_pars = P_SURG,
    events = c(hypopt_pth_sc(50,  ii = 24, days = 28),
               hypopt_pth_sc(100, ii = 24, days = D - 28, start = 28*24),
               hypopt_calcium(500, ii = 24, days = 28),
               hypopt_calcitriol(0.25, ii = 24, days = 28)))

  ## 7 — teriparatide 20 ug SC BID. Tall short peaks: good calcium control but
  ##     the pulse-memory term drives an anabolic bone response and BMD falls
  ##     from the artificially high HypoPT baseline toward normal.
  S[[7]] <- function() run(
    "07 테리파라타이드 30 ug SC BID",
    utils::modifyList(c(P_SURG, HYPOPT_pth_params$terip),
                      list(DCASUP = 0, DCTR = 0)),
    base_pars = P_SURG,
    events = hypopt_pth_sc(30, ii = 12, days = D))

  ## 8 — palopegteriparatide 18 ug daily, conventional therapy stopped on day 1.
  ##     Near-flat exposure, phosphate normalises, urine calcium falls BELOW the
  ##     conventional-therapy value at the same serum calcium. Triple composite
  ##     response (RESP3) becomes attainable here and in scenario 9 only.
  S[[8]] <- function() run(
    "08 팔로페그테리파라타이드 18 ug QD (기존치료 중단)",
    utils::modifyList(c(P_SURG, HYPOPT_pth_params$transcon),
                      list(DCASUP = 0, DCTR = 0)),
    base_pars = P_SURG,
    events = hypopt_transcon(18, ii = 24, days = D))

  ## 9 — palopegteriparatide titrated to 30 ug daily.
  S[[9]] <- function() run(
    "09 팔로페그테리파라타이드 30 ug QD (증량)",
    utils::modifyList(c(P_SURG, HYPOPT_pth_params$transcon),
                      list(DCASUP = 0, DCTR = 0)),
    base_pars = P_SURG,
    events = hypopt_transcon(30, ii = 24, days = D))

  ## 10 — ADH1, untreated. Chief cells are intact but the set-point is shifted
  ##      left, so PTH is low-normal — not absent — at a low serum calcium, and
  ##      the TAL CaSR is simultaneously over-active. Urine calcium is therefore
  ##      HIGH for the serum calcium, which is why conventional therapy is
  ##      especially dangerous in ADH1.
  S[[10]] <- function() run("10 ADH1 (CASR 기능획득) — 무치료", P_ADH1)

  ## 11 — ADH1 + encaleret 60 mg BID. The set-point is pushed back rightwards,
  ##      endogenous PTH rises, and tubular reabsorption is de-repressed:
  ##      serum calcium up AND urine calcium down at the same time. No other
  ##      agent in this model does both.
  S[[11]] <- function() run(
    "11 ADH1 + 엔칼레렛 60 mg BID",
    utils::modifyList(P_ADH1, list(DCASUP = 0, DCTR = 0)),
    base_pars = P_ADH1,
    events = hypopt_encaleret(60, ii = 12, days = D))

  ## 12 — functional hypoparathyroidism from profound magnesium depletion.
  ##      Chief-cell mass is normal but both secretion and PTH action are
  ##      blocked; calcium is refractory to calcium and calcitriol and corrects
  ##      abruptly once magnesium is replaced on day 30.
  S[[12]] <- function() run(
    "12 저마그네슘혈증성 기능성 HypoPT → Mg 보충 (30일)",
    P_MGLOW, base_pars = P_MGLOW,
    ## Magnesium repletion is DOSED, not switched on with a parameter. A single
    ## ev() record carrying a covariate column is applied by mrgsolve from t=0,
    ## which would silently make the patient magnesium-replete for the entire
    ## run — the opposite of the scenario.
    events = hypopt_magnesium(400, ii = 8, days = 60, start = 30*24),
    end = 90*24, delta = 2)

  ## 13 — acute post-thyroidectomy hypocalcaemic crisis. Chief-cell mass drops
  ##      to 2% at hour 0 in a previously normal patient; IV calcium gluconate
  ##      2 g over 1 h at hour 18 and 24, then an infusion, then oral therapy.
  S[[13]] <- function() run(
    "13 급성 저칼슘혈증 위기 — IV 칼슘 글루콘산",
    list(PTMASS = 0.02, SETPT = 1.20, DCASUP = 1500, DCTR = 1.0),
    base_pars = P_HEALTHY,
    events = c(hypopt_iv_calcium(2, hours = 1, start = 72),
               hypopt_iv_calcium(2, hours = 1, start = 78),
               hypopt_iv_calcium(10, hours = 12, start = 84),
               hypopt_calcium(500, ii = 8, days = 9, start = 72),
               hypopt_calcitriol(0.5, ii = 12, days = 9, start = 72)),
    end = 288, delta = 0.25)

  ## 14 — ten years of conventional therapy. This is where the renal cost is
  ##      actually paid: sustained urine calcium above the crystallisation
  ##      threshold accretes nephrocalcinosis, which accelerates GFR loss.
  S[[14]] <- function() run(
    "14 장기 10년 — 기존 치료 (신장 결과)",
    utils::modifyList(P_SURG, list(DCASUP = 1000, DCTR = 0.5)),
    base_pars = P_SURG,
    events = c(hypopt_calcium(500, ii = 12, days = 3650),
               hypopt_calcitriol(0.25, ii = 12, days = 3650)),
    end = LT, delta = 24*14)

  ## 15 — ten years of palopegteriparatide. Same disease, same patient, renal
  ##      arm restored. The GFR trajectories are the headline comparison.
  S[[15]] <- function() run(
    "15 장기 10년 — 팔로페그테리파라타이드 (신장 결과)",
    utils::modifyList(c(P_SURG, HYPOPT_pth_params$transcon),
                      list(DCASUP = 0, DCTR = 0)),
    base_pars = P_SURG,
    events = hypopt_transcon(24, ii = 24, days = 3650),
    end = LT, delta = 24*14)

  do.call(rbind, lapply(which, function(i) S[[i]]()))
}

## ---------------------------------------------------------------------------
##  Compact numeric summary of a scenario set
## ---------------------------------------------------------------------------
HYPOPT_summary <- function(sim) {
  sp <- split(sim, sim$scenario)
  last_frac <- function(d, col, frac = 0.2) {
    n <- nrow(d); k <- max(1, floor(n*(1 - frac)))
    mean(d[[col]][k:n], na.rm = TRUE)
  }
  out <- do.call(rbind, lapply(names(sp), function(nm) {
    d <- sp[[nm]]
    data.frame(
      scenario   = nm,
      Ca_mgdL    = round(last_frac(d, "CACORRD"), 2),
      iCa_mM     = round(last_frac(d, "CAIONMM"), 3),
      Pi_mgdL    = round(last_frac(d, "PISER"), 2),
      PTH_pgmL   = round(last_frac(d, "PTHPG"), 1),
      D125_pgmL  = round(last_frac(d, "D125PG"), 1),
      uCa_mgday  = round(last_frac(d, "UCA24"), 0),
      FECa_pct   = round(last_frac(d, "FEXCA"), 2),
      TmPGFR     = round(last_frac(d, "TMPPERG"), 2),
      turnover   = round(last_frac(d, "TURNOV"), 2),
      BMD_Z      = round(last_frac(d, "BMDZS"), 2),
      eGFR_end   = round(d$EGFR[nrow(d)], 1),
      NC_end     = round(d$NCBURD[nrow(d)], 2),
      symptom    = round(last_frac(d, "SYMSCOR"), 1),
      QTc_ms     = round(last_frac(d, "QTCMS"), 0),
      RESP3_pct  = round(100*mean(d$RESP3, na.rm = TRUE), 0),
      stringsAsFactors = FALSE)
  }))
  out[order(out$scenario), ]
}

## ---------------------------------------------------------------------------
##  Serum-calcium excursion within a dosing interval — the metric that
##  separates once-daily PTH(1-84) from a prodrug with flat exposure.
## ---------------------------------------------------------------------------
HYPOPT_swing <- function(sim, from_day = 150) {
  sp <- split(sim, sim$scenario)
  do.call(rbind, lapply(names(sp), function(nm) {
    d <- sp[[nm]]
    d <- d[d$time >= from_day*24, ]
    if (!nrow(d)) return(NULL)
    data.frame(scenario = nm,
               Ca_min   = round(min(d$CATOTAL), 2),
               Ca_max   = round(max(d$CATOTAL), 2),
               Ca_swing = round(max(d$CATOTAL) - min(d$CATOTAL), 2),
               PTH_min  = round(min(d$PTHPG), 1),
               PTH_max  = round(max(d$PTHPG), 1),
               PT_ratio = round(max(d$PTHPG)/pmax(1e-6, min(d$PTHPG)), 1),
               stringsAsFactors = FALSE)
  }))
}

## ---------------------------------------------------------------------------
##  Dose-titration helper: find the oral calcium / calcitriol combination that
##  lands albumin-corrected calcium in the 8.0-9.0 mg/dL guideline window, and
##  report what urinary calcium that costs.
## ---------------------------------------------------------------------------
HYPOPT_titrate_conventional <- function(mod = HYPOPT_build(),
                                        ca_doses  = c(250, 500, 750, 1000),
                                        ctr_doses = c(0.125, 0.25, 0.5, 0.75),
                                        days = 120) {
  init <- hypopt_baseline(mod, list(PTMASS = 0.05))
  grid <- expand.grid(ca = ca_doses, ctr = ctr_doses)
  res  <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    ev <- c(hypopt_calcium(grid$ca[i], ii = 8, days = days),
            hypopt_calcitriol(grid$ctr[i], ii = 12, days = days))
    m  <- mrgsolve::update(mod, param = list(PTMASS = 0.05),
                           init = as.list(init))
    d  <- mrgsolve::mrgsim_df(m, events = ev, end = days*24, delta = 12, hmax = 2)
    d  <- d[d$time >= (days - 14)*24, ]
    data.frame(ca_per_dose = grid$ca[i], ca_per_day = grid$ca[i]*3,
               calcitriol_per_day = grid$ctr[i]*2,
               Ca_mgdL  = round(mean(d$CACORRD), 2),
               uCa_mgd  = round(mean(d$UCA24), 0),
               in_target = mean(d$CACORRD) >= 8.0 && mean(d$CACORRD) <= 9.0,
               stringsAsFactors = FALSE)
  }))
  res[order(res$uCa_mgd), ]
}

## ---------------------------------------------------------------------------
##  Demo
## ---------------------------------------------------------------------------
if (identical(environment(), globalenv()) && !interactive()) {
  ## nothing runs on source(); call the functions explicitly
}

if (FALSE) {
  mod <- HYPOPT_build()
  sim <- HYPOPT_simulate_scenarios(mod, which = 1:13)
  print(HYPOPT_summary(sim))
  print(HYPOPT_swing(sim, from_day = 150))

  long <- HYPOPT_simulate_scenarios(mod, which = 14:15)
  print(HYPOPT_summary(long))

  print(HYPOPT_titrate_conventional(mod))
}
